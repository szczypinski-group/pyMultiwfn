"""Result container for Multiwfn job execution."""

from dataclasses import dataclass
from typing import Any

from pymultiwfn.analysis.parsers import (
    BondOrderParser,
    ChargeParser,
    CriticalPointParser,
    SpectrumParser,
)


@dataclass
class MultiwfnResult:
    """Encapsulates parsed results from a Multiwfn execution."""

    type: str
    value: Any

    def parse_charges(
        self,
        stdout: str,
    ) -> dict[int, float]:
        """Parse atomic charges from output.

        Parameters
        ----------
        method : str
            Charge method name (e.g., "hirshfeld", "mulliken", "adch")

        Returns
        -------
        dict[int, float]
            Dictionary mapping atom indices to charges
        """
        return ChargeParser.parse(stdout)

    def parse_bond_orders(
        self,
        stdout: str,
    ) -> dict[tuple[int, int], float]:
        """Parse bond orders from output.

        Returns
        -------
        dict[tuple[int, int], float]
            Dictionary mapping atom pairs to bond orders
        """
        return BondOrderParser.parse(stdout)

    def parse_critical_points(
        self,
        stdout: str,
    ) -> list[dict[str, Any]]:
        """Parse critical points from output.

        Returns
        -------
        list[dict[str, Any]]
            List of critical point dictionaries
        """
        return CriticalPointParser.parse(stdout)

    def parse_spectrum(
        self,
        stdout: str,
    ) -> dict[str, list[float]]:
        """Parse spectrum data from output.

        Returns
        -------
        dict[str, list[float]]
            Dictionary with 'frequencies' and 'intensities' lists
        """
        return SpectrumParser.parse(stdout)
