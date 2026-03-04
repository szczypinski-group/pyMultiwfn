"""Result container for Multiwfn job execution."""

from dataclasses import dataclass
from typing import Any

from pymultiwfn.api.parsers import (
    BondOrderParser,
    ChargeParser,
    CriticalPointParser,
    SpectrumParser,
)


@dataclass
class MultiwfnResult:
    """Encapsulates results from a Multiwfn execution.

    Attributes
    ----------
    stdout : str
        Standard output from Multiwfn
    stderr : str
        Standard error from Multiwfn
    returncode : int
        Process return code
    execution_time : float
        Execution time in seconds
    commands : list of str
        Commands that were executed
    input_file : Path
        Input file that was processed
    """

    type: str
    value: Any

    def parse_charges(
        self,
        stdout: str,
        method: str = "hirshfeld",
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
        return ChargeParser.parse(stdout, method=method)

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
