"""Output parsers for Multiwfn result.

Especially those with perhaps more complicated or worse formatted outputs.

"""

import re
from typing import Any

# Shared float pattern for all numeric parsing
# Matches: 1, -1, 1., .5, -0.5, 1e10, .5E-1, -1.2e+03, 1.23E-02
FLOAT_PATTERN = r"[-+]?(?:\d+\.\d*|\d*\.\d+|\d+)(?:[Ee][-+]?\d+)?"


class OutputParser:
    """Parent class for output parsers."""

    pass


class ChargeParser(OutputParser):
    """Parser for atomic charges."""

    def parse(
        self, stdout: str, method: str = "Hirshfeld"
    ) -> dict[int, float]:
        """Extract atomic charges from Multiwfn output.

        Parameters
        ----------
        stdout : str
            Multiwfn standard output
        method : str
            Charge method name (e.g., "Hirshfeld", "ADCH", "RESP")

        Returns:
        -------
        dict
            Dictionary mapping atom indices to charges
        """
        charges = {}
        # Match atom id, element in parentheses, then a floating-point charge
        # Supports optional sign, leading zero, and scientific notation
        pattern = rf"^\s*(\d+)\s*\([A-Za-z]+\s*\)\s+({FLOAT_PATTERN})"

        in_charge_section = False
        for line in stdout.split("\n"):
            if (
                method.lower() in line.lower()
                or "atomic charges" in line.lower()
            ):
                in_charge_section = True
                continue

            if in_charge_section:
                if match := re.match(pattern, line):
                    atom_idx = int(match[1])
                    charge = float(match[2])
                    charges[atom_idx] = charge
                elif (
                    line.strip() == ""
                    or "sum of atomic charges" in line.lower()
                ):
                    if charges:
                        break

        return charges


class BondOrderParser(OutputParser):
    """Parser for bond orders."""

    def parse(
        self,
        stdout: str,
    ) -> dict[tuple[int, int], float]:
        """Extract bond orders from Multiwfn output.

        Parameters
        ----------
        stdout : str
            Multiwfn standard output

        Returns:
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

        for line in stdout.split("\n"):
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
        "(3,-3)": "nuclear",
        "(3,-1)": "bond",
        "(3,+1)": "ring",
        "(3,+3)": "cage",
    }

    def parse(
        self,
        stdout: str,
    ) -> list[dict[str, Any]]:
        """Extract critical point information from topology analysis.

        Parameters
        ----------
        stdout : str
            Multiwfn standard output

        Returns:
        -------
        list
            List of dictionaries containing CP information.
            Position is None if not found in output.
        """
        cps = []

        # Pattern for format: "CP   1 (3,-3) Nuclear critical point"
        # followed by "Position (Bohr):    0.000000    0.000000    0.000000"
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

                position = None
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
                    "cp_type": self.CP_TYPE_NAMES.get(cp_type, "unknown"),
                    "position": position,
                }
                cps.append(cp)

        return cps


class SpectrumParser(OutputParser):
    """Parser for spectrum data."""

    def parse(
        self,
        stdout: str,
    ) -> dict[str, list[float]]:
        """Extract spectrum data (frequencies, intensities).

        Parameters
        ----------
        stdout : str
            Multiwfn standard output

        Returns:
        -------
        dict
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
                match2 = re.match(pattern2, line)
                if match2 is not None:
                    spectrum["frequencies"].append(float(match2[1]))
                    spectrum["intensities"].append(float(match2[2]))

        return spectrum
