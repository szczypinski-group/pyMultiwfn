"""Result container for Multiwfn job execution."""

from dataclasses import dataclass
from pathlib import Path


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
        # Check stderr for error indicators
        error_indicators = [
            "error",
            "fatal",
            "cannot open",
            "not found",
            "failed to",
        ]
        has_error = any(ind in self.stderr.lower() for ind in error_indicators)
        return not has_error

    # def parse(
    #     self,
    #     parser_name: str,
    # ) -> Any:
    #     """Parse output using a registered parser.

    #     Parameters
    #     ----------
    #     parser_name : str
    #         Name of the parser to use
    #     **kwargs
    #         Additional arguments to pass to the parser

    #     Returns:
    #     -------
    #     Any
    #         Parsed data
    #     """
    #     return ParserRegistry.parse(parser_name, self.stdout, **kwargs)

    # def parse_charges(self, method: str = "hirshfeld") -> dict[int, float]:
    #     """Parse atomic charges from output."""
    #     return self.parse("charges", method=str(method))

    # def parse_bond_orders(self) -> dict[tuple[int, int], float]:
    #     """Parse bond orders from output."""
    #     return self.parse("bond_orders")

    # def parse_critical_points(self) -> list[dict[str, Any]]:
    #     """Parse critical points from output."""
    #     return self.parse("critical_points")

    # def parse_spectrum(self) -> dict[str, list[float]]:
    #     """Parse spectrum data from output."""
    #     return self.parse("spectrum")

    # def save_output(self, filename: str | Path) -> None:
    #     """Save output to file."""
    #     with open(filename, "w", encoding="utf-8") as f:
    #         f.write(self.stdout)

    # def __str__(self) -> str:
    #     """Human-readable string representation."""
    #     status = (
    #         "SUCCESS" if self.success else f"FAILED (code {self.returncode})"
    #     )
    #     return (
    #               f"MultiwfnResult({status}, {self.execution_time:.2f}s, F
    #               f"{len(self.commands)} commands)"
    #     )

    # def __repr__(self) -> str:
    #     """Detailed string representation for debugging."""
    #     return (
    #         f"MultiwfnResult("
    #         f"returncode={self.returncode}, "
    #         f"execution_time={self.execution_time:.3f}, "
    #         f"input_file={self.input_file!r}, "
    #         f"commands={self.commands!r}, "
    #         f"stdout_len={len(self.stdout)}, "
    #         f"stderr_len={len(self.stderr)})"
    #     )
