"""Configuration interaction with Multiwfn.

This class relates to the actual Multiwfn executable (e.g., executable path),
rather than the details of a specific job/analysis (e.g., timeout and parsing).

"""

import platform
from pathlib import Path

from pymultiwfn.api.exceptions import MultiwfnError


class Multiwfn:
    """Configuration for Multiwfn execution.

    Attributes
    ----------
    exe_path : Path, optional
        Explicit path to Multiwfn executable. If not provided,
        looks for executable in bin/ directory relative to this package.

    Examples
    --------
    >>> # Simple configuration
    >>> multiwfn = Multiwfn()

    """

    def __init__(
        self,
        exe_path: Path | None = None,
    ) -> None:
        self.exe_path = self._parse_exe(exe_path)

    def _parse_exe(self, exe_path: Path | None) -> Path:
        """Parse the executable path."""
        if exe_path:
            if exe_path.resolve().exists():
                return exe_path.resolve()
            raise MultiwfnError(
                f"Specified executable not found: {self.exe_path.resolve()}"
            )

        if platform.system() == "Windows":
            bin_path = Path(__package__).resolve().parent / "bin/Multifwn.exe"

        elif platform.system() == "Linux":
            bin_path = Path(__package__).resolve().parent / "bin/Multifwn"

        else:
            raise MultiwfnError(
                f"No Multiwfn executable found for {platform.system()} OS."
            )

        if bin_path.exists():
            return bin_path

        else:
            raise MultiwfnError(
                f"Multiwfn executable not found at {bin_path}. "
                "Please specify the executable path."
            )

    def __str__(self) -> str:
        return f"{self.exe_path}"

    def __repr__(self) -> str:
        return f"pymultiwfn.Multiwfn(exe_path={self.exe_path})"
