"""Configuration management for pyMultiwfn."""

import platform
from dataclasses import dataclass, field
from pathlib import Path

from pymultiwfn.api.exceptions import MultiwfnError


@dataclass
class MultiwfnConfig:
    """Configuration for Multiwfn execution.

    Attributes
    ----------
    exe_path : Path, optional
        Explicit path to Multiwfn executable. If not provided,
        looks for executable in bin/ directory relative to this package.
    working_dir : Path
        Working directory for output files
    timeout : int or TimeoutConfig, optional
        Timeout configuration. Can be:
        - None: Use default TimeoutConfig
        - int: Use this value for all analyses
        - TimeoutConfig: Use specific timeouts for different analyses
    verbose : bool
        Whether to print output during execution

    Examples
    --------
    >>> # Simple configuration
    >>> config = MultiwfnConfig(timeout=300, verbose=True)

    >>> # Advanced timeout configuration
    >>> timeout = TimeoutConfig(default=120, topology=600)
    >>> config = MultiwfnConfig(timeout=timeout)
    """

    exe_path: Path | None = None
    working_dir: Path = field(default_factory=Path.cwd)
    timeout: int | None = 120
    verbose: bool = False

    _resolved_exe: Path | None = field(default=None, init=False, repr=False)

    def __post_init__(self) -> None:
        """Validate and normalize configuration after initialization."""
        # Normalize working_dir to Path
        if isinstance(self.working_dir, str):
            object.__setattr__(self, "working_dir", Path(self.working_dir))

        # Normalize exe_path to Path
        if isinstance(self.exe_path, str):
            object.__setattr__(self, "exe_path", Path(self.exe_path))

        # Set up timeout configuration

    @property
    def executable(self) -> Path:
        """Get the resolved Multiwfn executable path."""
        if self._resolved_exe is None:
            object.__setattr__(self, "_resolved_exe", self._find_executable())
        return self._resolved_exe  # type: ignore[return-value]

    def _find_executable(self) -> Path:
        """Find Multiwfn executable in bin/ directory."""
        if self.exe_path:
            path = Path(self.exe_path)
            if path.exists():
                return path
            raise MultiwfnError(
                f"Specified executable not found: {self.exe_path}"
            )

        is_windows = platform.system() == "Windows"
        exe_name = "Multiwfn.exe" if is_windows else "Multiwfn"

        # Check bin/ directory relative to this package
        bin_exe = Path(__package__).resolve().parent / "bin" / exe_name
        if bin_exe.exists():
            return bin_exe

        raise MultiwfnError(
            f"Multiwfn executable not found at {bin_exe}. "
            "Please place the executable in the bin/ directory "
            "or specify exe_path."
        )
