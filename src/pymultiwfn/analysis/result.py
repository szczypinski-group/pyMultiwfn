"""Result container for Multiwfn job execution."""

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from pymultiwfn.interface.parsers import (
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

    stdout: str
    stderr: str
    returncode: int
    execution_time: float
    commands: list[str]
    input_file: Path

    @property
    def success(self) -> bool:
        """Check if execution was successful.

        Notes
        -----
        Multiwfn often returns non-zero codes in batch mode even on success.
        We check for actual error indicators instead of just return code.
        """
        if self.returncode is None or self.returncode < 0:
            return False
        error_indicators = [
            "error",
            "fatal",
            "cannot open",
            "not found",
            "failed to",
        ]
        has_error = any(ind in self.stderr.lower() for ind in error_indicators)
        return not has_error

    def parse_charges(self, method: str = "hirshfeld") -> dict[int, float]:
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
        return ChargeParser.parse(self.stdout, method=method)

    def parse_bond_orders(self) -> dict[tuple[int, int], float]:
        """Parse bond orders from output.

        Returns
        -------
        dict[tuple[int, int], float]
            Dictionary mapping atom pairs to bond orders
        """
        return BondOrderParser.parse(self.stdout)

    def parse_critical_points(self) -> list[dict[str, Any]]:
        """Parse critical points from output.

        Returns
        -------
        list[dict[str, Any]]
            List of critical point dictionaries
        """
        return CriticalPointParser.parse(self.stdout)

    def parse_spectrum(self) -> dict[str, list[float]]:
        """Parse spectrum data from output.

        Returns
        -------
        dict[str, list[float]]
            Dictionary with 'frequencies' and 'intensities' lists
        """
        return SpectrumParser.parse(self.stdout)

    def save_output(self, filename: str | Path) -> None:
        """Save stdout to file.

        Parameters
        ----------
        filename : str or Path
            Output file path
        """
        Path(filename).write_text(self.stdout, encoding="utf-8")

    def __str__(self) -> str:
        """Human-readable string representation."""
        status = (
            "SUCCESS" if self.success else f"FAILED (code {self.returncode})"
        )
        return (
            f"MultiwfnResult({status}, {self.execution_time:.2f}s, "
            f"{len(self.commands)} commands)"
        )

    def __repr__(self) -> str:
        """Detailed string representation for debugging."""
        return (
            f"MultiwfnResult("
            f"returncode={self.returncode}, "
            f"execution_time={self.execution_time:.3f}, "
            f"input_file={self.input_file!r}, "
            f"commands={self.commands!r}, "
            f"stdout_len={len(self.stdout)}, "
            f"stderr_len={len(self.stderr)})"
        )
