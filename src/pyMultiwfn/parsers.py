"""Output parsers for Multiwfn results, especially those with perhaps more complciated or less well formatted outputs."""

import re
from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional, Tuple, Type

from .config import MultiwfnError


# Shared float pattern for all numeric parsing
# Matches: 1, -1, 1., .5, -0.5, 1e10, .5E-1, -1.2e+03, 1.23E-02
FLOAT_PATTERN = r"[-+]?(?:\d+\.\d*|\d*\.\d+|\d+)(?:[Ee][-+]?\d+)?"


class OutputParser(ABC):
    """Abstract base class for output parsers."""
    
    @abstractmethod
    def parse(self, stdout: str, **kwargs) -> Any:
        """Parse Multiwfn output and return structured data."""
        pass


class ChargeParser(OutputParser):
    """Parser for atomic charges."""
    
    def parse(self, stdout: str, method: str = "Hirshfeld") -> Dict[int, float]:
        """
        Extract atomic charges from Multiwfn output.
        
        Parameters
        ----------
        stdout : str
            Multiwfn standard output
        method : str
            Charge method name (e.g., "Hirshfeld", "ADCH", "RESP")
            
        Returns
        -------
        dict
            Dictionary mapping atom indices to charges
        """
        charges = {}
        # Match atom index, element in parentheses, then a floating-point charge
        # Supports optional sign, optional leading zero, and scientific notation
        pattern = rf"^\s*(\d+)\s*\([A-Za-z]+\s*\)\s+({FLOAT_PATTERN})"

        in_charge_section = False
        for line in stdout.split('\n'):
            if method.lower() in line.lower() or "atomic charges" in line.lower():
                in_charge_section = True
                continue

            if in_charge_section:
                if match := re.match(pattern, line):
                    atom_idx = int(match[1])
                    charge = float(match[2])
                    charges[atom_idx] = charge
                elif line.strip() == "" or "sum of atomic charges" in line.lower():
                    if charges:
                        break

        return charges


class BondOrderParser(OutputParser):
    """Parser for bond orders."""
    
    def parse(self, stdout: str, **kwargs) -> Dict[Tuple[int, int], float]:
        """
        Extract bond orders from Multiwfn output.
        
        Parameters
        ----------
        stdout : str
            Multiwfn standard output
            
        Returns
        -------
        dict
            Dictionary mapping atom pairs to bond orders
        """
        bond_orders = {}
        # Pattern 1: "   1  -    2    1.4523"
        # Pattern 2: "   1(C ) -    2(C ):  1.4523"
        patterns = [
            rf"^\s*(\d+)\s*-\s*(\d+)\s+({FLOAT_PATTERN})",
            rf"^\s*(\d+)\([^)]+\)\s*-\s*(\d+)\([^)]+\)\s*:\s*({FLOAT_PATTERN})",
        ]

        for line in stdout.split('\n'):
            for pattern in patterns:
                if match := re.match(pattern, line):
                    atom1 = int(match[1])
                    atom2 = int(match[2])
                    bo = float(match[3])
                    if atom1 > atom2:
                        atom1, atom2 = atom2, atom1
                    bond_orders[(atom1, atom2)] = bo
                    break

        return bond_orders


class CriticalPointParser(OutputParser):
    """Parser for critical point information."""
    
    CP_TYPE_NAMES = {
        '(3,-3)': 'nuclear',
        '(3,-1)': 'bond',
        '(3,+1)': 'ring',
        '(3,+3)': 'cage',
    }
    
    def parse(self, stdout: str, **kwargs) -> List[Dict[str, Any]]:
        """
        Extract critical point information from topology analysis.
        
        Parameters
        ----------
        stdout : str
            Multiwfn standard output
            
        Returns
        -------
        list
            List of dictionaries containing CP information.
            Position is None if not found in output.
        """
        cps = []

        # Pattern for format: "CP   1 (3,-3) Nuclear critical point"
        # followed by "Position (Bohr):    0.000000    0.000000    0.000000"
        pattern = r"CP\s+(\d+)\s+\((\d+),([+-]?\d+)\)"
        pos_pattern = rf"Position.*?:\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})"

        lines = stdout.split('\n')
        for i, line in enumerate(lines):
            if match := re.search(pattern, line):
                cp_index = int(match[1])
                cp_type = f"({match[2]},{match[3]})"

                position = None
                for j in range(i, min(i + 5, len(lines))):
                    if pos_match := re.search(pos_pattern, lines[j]):
                        position = float(pos_match[1]), float(pos_match[2]), float(pos_match[3])
                        break

                cp = {
                    'index': cp_index,
                    'type': cp_type,
                    'cp_type': self.CP_TYPE_NAMES.get(cp_type, 'unknown'),
                    'position': position
                }
                cps.append(cp)

        return cps


class SpectrumParser(OutputParser):
    """Parser for spectrum data."""
    
    def parse(self, stdout: str, **kwargs) -> Dict[str, List[float]]:
        """
        Extract spectrum data (frequencies, intensities).
        
        Parameters
        ----------
        stdout : str
            Multiwfn standard output
            
        Returns
        -------
        dict
            Dictionary with 'frequencies' and 'intensities' lists
        """
        spectrum = {'frequencies': [], 'intensities': []}

        pattern1 = rf"({FLOAT_PATTERN})\s+cm\^?-1.*?Intensity:\s+({FLOAT_PATTERN})"
        pattern2 = rf"^\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s*$"

        for match in re.finditer(pattern1, stdout):
            spectrum['frequencies'].append(float(match.group(1)))
            spectrum['intensities'].append(float(match.group(2)))

        if not spectrum['frequencies']:
            for line in stdout.split('\n'):
                if match := re.match(pattern2, line):
                    spectrum['frequencies'].append(float(match[1]))
                    spectrum['intensities'].append(float(match[2]))

        return spectrum


class ParserRegistry:
    """
    Registry for output parsers.
    
    Allows registration of custom parsers and provides a unified interface
    for parsing Multiwfn output.
    """
    
    _parsers: Dict[str, OutputParser] = {}
    
    @classmethod
    def register(cls, name: str, parser: Optional[OutputParser] = None):
        """
        Register a parser with the registry.
        
        Can be used as a decorator or called directly.
        
        Parameters
        ----------
        name : str
            Name to register the parser under
        parser : OutputParser, optional
            Parser instance (if not using as decorator)
        """
        def decorator(parser_cls: Type[OutputParser]):
            cls._parsers[name] = parser_cls()
            return parser_cls
        
        if parser is not None:
            cls._parsers[name] = parser
            return None
        return decorator
    
    @classmethod
    def get(cls, name: str) -> OutputParser:
        """Get a parser by name."""
        if name not in cls._parsers:
            raise MultiwfnError(f"Unknown parser: {name}")
        return cls._parsers[name]
    
    @classmethod
    def parse(cls, name: str, stdout: str, **kwargs) -> Any:
        """
        Parse output using a registered parser.
        
        Parameters
        ----------
        name : str
            Name of the parser to use
        stdout : str
            Multiwfn output to parse
        **kwargs
            Additional arguments to pass to the parser
            
        Returns
        -------
        Any
            Parsed data
        """
        return cls.get(name).parse(stdout, **kwargs)
    
    @classmethod
    def list_parsers(cls) -> List[str]:
        """List all registered parser names."""
        return list(cls._parsers.keys())


ParserRegistry.register('charges', ChargeParser())
ParserRegistry.register('bond_orders', BondOrderParser())
ParserRegistry.register('critical_points', CriticalPointParser())
ParserRegistry.register('spectrum', SpectrumParser())