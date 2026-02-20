"""Output parsers for Multiwfn results."""

import re
from typing import Any

# Shared float pattern for all numeric parsing
FLOAT_PATTERN = r"[-+]?(?:\d+\.\d*|\d*\.\d+|\d+)(?:[Ee][-+]?\d+)?"


class OutputParser:
    """Base class for output parsers."""

    pass


class ChargeParser(OutputParser):
    """Parser for atomic charges."""

    @staticmethod
    def parse(stdout: str, method: str = "Hirshfeld") -> dict[int, float]:
        """Extract atomic charges from Multiwfn output.

        Parameters
        ----------
        stdout : str
            Multiwfn standard output
        method : str
            Charge method name (e.g., "Hirshfeld", "ADCH", "RESP")

        Returns
        -------
        dict[int, float]
            Dictionary mapping atom indices to charges
        """
        charges: dict[int, float] = {}

        # Pattern 1: "Hirshfeld charge of atom     1(C ) is  0.03208687"
        pattern1 = (
            rf"{method}\s+charge of atom\s+(\d+)\s*\([A-Za-z]+\s*\)\s+is\s+"
            rf"({FLOAT_PATTERN})"
        )
        # Pattern 2: "Atom    1(C ):     0.03209323" (Final atomic charges)
        pattern2 = rf"Atom\s+(\d+)\s*\([A-Za-z]+\s*\)\s*:\s*({FLOAT_PATTERN})"
        # Pattern 3: "    1(C )   -0.0523" (summary table format)
        pattern3 = rf"^\s*(\d+)\s*\([A-Za-z]+\s*\)\s+({FLOAT_PATTERN})"
        # Pattern 4: "   1  C      -0.052300" (column format)
        pattern4 = rf"^\s*(\d+)\s+[A-Za-z]+\s+({FLOAT_PATTERN})"

        in_charge_section = False
        in_final_section = False

        for line in stdout.split("\n"):
            line_lower = line.lower()

            # Check for section markers
            if "final atomic charges" in line_lower:
                in_final_section = True
                in_charge_section = True
                charges.clear()  # Prefer final charges
                continue
            elif method.lower() in line_lower and "charge" in line_lower:
                in_charge_section = True
                continue

            # Try pattern 1 first (most specific - includes method name)
            if match := re.search(pattern1, line, re.IGNORECASE):
                atom_idx = int(match[1])
                charge = float(match[2])
                charges[atom_idx] = charge
                continue

            # Try pattern 2 (Final atomic charges format)
            if match := re.search(pattern2, line):
                if in_charge_section or in_final_section:
                    atom_idx = int(match[1])
                    charge = float(match[2])
                    charges[atom_idx] = charge
                continue

            # Try pattern 3 (table format)
            if in_charge_section and (match := re.match(pattern3, line)):
                atom_idx = int(match[1])
                charge = float(match[2])
                charges[atom_idx] = charge
                continue

            # Try pattern 4 (column format)
            if in_charge_section and (match := re.match(pattern4, line)):
                atom_idx = int(match[1])
                charge = float(match[2])
                charges[atom_idx] = charge
                continue

            # End of section detection
            if in_final_section and charges and "calculation took" in line_lower:
                break

        return charges


class BondOrderParser(OutputParser):
    """Parser for bond orders."""

    @staticmethod
    def parse(
        stdout: str, **kwargs: OutputParser
    ) -> dict[tuple[int, int], float]:
        """Extract bond orders from Multiwfn output.

        Parameters
        ----------
        stdout : str
            Multiwfn standard output

        Returns
        -------
        dict[tuple[int, int], float]
            Dictionary mapping atom pairs to bond orders
        """
        bond_orders: dict[tuple[int, int], float] = {}

        # Pattern 1: "#    1:         1(C )    2(N )    1.95394965"
        pattern1 = (
            rf"#\s*\d+:\s+(\d+)\s*\([^)]+\)\s+(\d+)\s*\([^)]+\)\s+"
            rf"({FLOAT_PATTERN})"
        )

        # Pattern 2: "1(C ) -    2(C ):  1.4523"
        pattern2 = (
            rf"(\d+)\s*\([^)]+\)\s*-\s*(\d+)\s*\([^)]+\)\s*:\s*"
            rf"({FLOAT_PATTERN})"
        )

        # Pattern 3: Simple "1 - 2   1.4523"
        pattern3 = rf"^\s*(\d+)\s*-\s*(\d+)\s+({FLOAT_PATTERN})"

        # Pattern 4: "   1    2    1.4523" (column format)
        pattern4 = rf"^\s*(\d+)\s+(\d+)\s+({FLOAT_PATTERN})\s*$"

        in_bond_section = False

        for line in stdout.split("\n"):
            line_lower = line.lower()

            # Check for section markers
            if "bond order" in line_lower or "bond-order" in line_lower:
                in_bond_section = True
                continue

            match = None

            # Try all patterns
            if match := re.search(pattern1, line):
                pass
            if match:
                atom1 = int(match[1])
                atom2 = int(match[2])
                bo = float(match[3])
                # Ensure consistent ordering (smaller index first)
                if atom1 > atom2:
                    atom1, atom2 = atom2, atom1
                bond_orders[(atom1, atom2)] = bo

        return bond_orders


class CriticalPointParser(OutputParser):
    """Parser for critical point information."""

    CP_TYPE_NAMES: dict[str, str] = {
        "(3,-3)": "nuclear",
        "(3,-1)": "bond",
        "(3,+1)": "ring",
        "(3,+3)": "cage",
    }

    @staticmethod
    def parse(stdout: str, **kwargs: OutputParser) -> list[dict[str, Any]]:
        """Extract critical point information from topology analysis.

        Parameters
        ----------
        stdout : str
            Multiwfn standard output

        Returns
        -------
        list[dict[str, Any]]
            List of dictionaries containing CP information
        """
        cps: list[dict[str, Any]] = []

        pattern = r"CP\s+(\d+)\s+\((\d+),([+-]?\d+)\)"
        pos_pattern = (
            rf"Position.*?:\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})"
            rf"\s+({FLOAT_PATTERN})"
        )

        lines = stdout.split("\n")
        for i, line in enumerate(lines):
            if match := re.search(pattern, line):
                cp_index = int(match[1])
                cp_type = f"({match[2]},{match[3]})"

                position: tuple[float, float, float] | None = None
                for j in range(i, min(i + 5, len(lines))):
                    if pos_match := re.search(pos_pattern, lines[j]):
                        position = (
                            float(pos_match[1]),
                            float(pos_match[2]),
                            float(pos_match[3]),
                        )
                        break

                cp = {
                    "index": cp_index,
                    "type": cp_type,
                    "cp_type": CriticalPointParser.CP_TYPE_NAMES.get(
                        cp_type, "unknown"
                    ),
                    "position": position,
                }
                cps.append(cp)

        return cps


class SpectrumParser(OutputParser):
    """Parser for spectrum data."""

    @staticmethod
    def parse(stdout: str, **kwargs: OutputParser) -> dict[str, list[float]]:
        """Extract spectrum data (frequencies, intensities).

        Parameters
        ----------
        stdout : str
            Multiwfn standard output

        Returns
        -------
        dict[str, list[float]]
            Dictionary with 'frequencies' and 'intensities' lists
        """
        spectrum: dict[str, list[float]] = {
            "frequencies": [],
            "intensities": [],
        }

        pattern1 = (
            rf"({FLOAT_PATTERN})\s+cm\^?-1.*?Intensity:\s+({FLOAT_PATTERN})"
        )
        pattern2 = rf"^\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s*$"

        for match in re.finditer(pattern1, stdout):
            spectrum["frequencies"].append(float(match.group(1)))
            spectrum["intensities"].append(float(match.group(2)))

        if not spectrum["frequencies"]:
            for line in stdout.split("\n"):
                if match2 := re.match(pattern2, line):
                    spectrum["frequencies"].append(float(match2[1]))
                    spectrum["intensities"].append(float(match2[2]))

        return spectrum