"""Raw results container for MultiwfnJob execution."""

from dataclasses import dataclass
from pathlib import Path


@dataclass
class MultiwfnJobOutcome:
    """Encapsulates results from a Multiwfn execution.

    Attributes
    ----------
    stdout
        Standard output from Multiwfn
    stderr
        Standard error from Multiwfn
    return_code
        Process return code
    execution_time
        Execution time in seconds

    """

    stdout: str
    stderr: str
    return_code: int | None
    execution_time: float

    @property
    def success(self) -> bool:
        """Check if execution was successful.

        Notes
        -----
        Multiwfn often returns non-zero codes in batch mode even on success.
        We check for actual error indicators instead of just return code.
        """
        if self.return_code is None or self.return_code < 0:
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

    def save_output(self, filename: str | Path) -> None:
        """Save stdout to file.

        Parameters
        ----------
        filename : str or Path
            Output file path
        """
        Path(filename).write_text(self.stdout, encoding="utf-8")
