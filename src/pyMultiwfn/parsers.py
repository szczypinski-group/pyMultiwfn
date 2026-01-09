"""Output parsers for Multiwfn results."""

import re
from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional, Tuple, Type

from .config import MultiwfnError


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
        pattern = r"^\s*(\d+)\s*\([A-Za-z]+\s*\)\s+([-+]?\d+\.\d+)"
        
        in_charge_section = False
        for line in stdout.split('\n'):
            if method.lower() in line.lower() or "atomic charges" in line.lower():
                in_charge_section = True
                continue
            
            if in_charge_section:
                match = re.match(pattern, line)
                if match:
                    atom_idx = int(match.group(1))
                    charge = float(match.group(2))
                    charges[atom_idx] = charge
                elif line.strip() == "" or "sum" in line.lower():
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
            r"^\s*(\d+)\s*-\s*(\d+)\s+([-+]?\d+\.\d+)",
            r"^\s*(\d+)\([^)]+\)\s*-\s*(\d+)\([^)]+\)\s*:\s*([-+]?\d+\.\d+)",
        ]
        
        for line in stdout.split('\n'):
            for pattern in patterns:
                match = re.match(pattern, line)
                if match:
                    atom1 = int(match.group(1))
                    atom2 = int(match.group(2))
                    bo = float(match.group(3))
                    # Ensure consistent ordering (smaller index first)
                    if atom1 > atom2:
                        atom1, atom2 = atom2, atom1
                    bond_orders[(atom1, atom2)] = bo
                    break
        
        return bond_orders


class CriticalPointParser(OutputParser):
    """Parser for critical point information."""
    
    # Map CP types to names
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
            List of dictionaries containing CP information
        """
        cps = []
        
        # Pattern for format: "CP   1 (3,-3) Nuclear critical point"
        # followed by "Position (Bohr):    0.000000    0.000000    0.000000"
        pattern = r"CP\s+(\d+)\s+\((\d+),([+-]?\d+)\)"
        pos_pattern = r"Position.*?:\s+([-+]?\d+\.?\d*[Ee]?[+-]?\d*)\s+([-+]?\d+\.?\d*[Ee]?[+-]?\d*)\s+([-+]?\d+\.?\d*[Ee]?[+-]?\d*)"
        
        lines = stdout.split('\n')
        for i, line in enumerate(lines):
            match = re.search(pattern, line)
            if match:
                cp_index = int(match.group(1))
                cp_type = f"({match.group(2)},{match.group(3)})"
                
                # Look for position in next few lines
                position = (0.0, 0.0, 0.0)
                for j in range(i, min(i + 5, len(lines))):
                    pos_match = re.search(pos_pattern, lines[j])
                    if pos_match:
                        position = (
                            float(pos_match.group(1)),
                            float(pos_match.group(2)),
                            float(pos_match.group(3))
                        )
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
        
        # Pattern 1: "500.00 cm^-1 ... Intensity: 12.34"
        # Pattern 2: "    500.00           12.34" (two columns)
        pattern1 = r"(\d+\.?\d*)\s+cm\^?-1.*?Intensity:\s+([-+]?\d+\.?\d*)"
        pattern2 = r"^\s+(\d+\.?\d+)\s+([-+]?\d+\.?\d+)\s*$"
        
        for match in re.finditer(pattern1, stdout):
            spectrum['frequencies'].append(float(match.group(1)))
            spectrum['intensities'].append(float(match.group(2)))
        
        # If pattern1 didn't match, try pattern2
        if not spectrum['frequencies']:
            for line in stdout.split('\n'):
                match = re.match(pattern2, line)
                if match:
                    spectrum['frequencies'].append(float(match.group(1)))
                    spectrum['intensities'].append(float(match.group(2)))
        
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


# Register default parsers
ParserRegistry.register('charges', ChargeParser())
ParserRegistry.register('bond_orders', BondOrderParser())
ParserRegistry.register('critical_points', CriticalPointParser())
ParserRegistry.register('spectrum', SpectrumParser())