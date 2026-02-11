"""Configuration management for pyMultiwfn."""

import platform
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


class MultiwfnError(Exception):
    """Base exception for pyMultiwfn errors."""
    pass


@dataclass(slots=True)
class MultiwfnConfig:
    """
    Configuration for Multiwfn execution.
    
    Attributes
    ----------
    exe_path : Path, optional
        Explicit path to Multiwfn executable. If not provided,
        looks for executable in bin/ directory relative to this package.
    working_dir : Path
        Working directory for output files
    timeout : int, optional
        Default timeout in seconds for job execution
    verbose : bool
        Whether to print output during execution
    """
    exe_path: Optional[Path] = None
    working_dir: Path = field(default_factory=Path.cwd)
    timeout: Optional[int] = None
    verbose: bool = False
    
    _resolved_exe: Optional[Path] = field(default=None, init=False, repr=False)
    
    @property
    def executable(self) -> Path:
        """Get the resolved Multiwfn executable path."""
        if self._resolved_exe is None:
            object.__setattr__(self, '_resolved_exe', self._find_executable())
        return self._resolved_exe
    
    def _find_executable(self) -> Path:
        """Find Multiwfn executable in bin/ directory."""
        if self.exe_path:
            path = Path(self.exe_path)
            if path.exists():
                return path
            raise MultiwfnError(f"Specified executable not found: {self.exe_path}")
        
        is_windows = platform.system() == "Windows"
        exe_name = "Multiwfn.exe" if is_windows else "Multiwfn"
        
        # Check bin/ directory relative to this package
        bin_exe = Path(__file__).resolve().parent / "bin" / exe_name
        if bin_exe.exists():
            return bin_exe
        
        raise MultiwfnError(
            f"Multiwfn executable not found at {bin_exe}. "
            "Please place the executable in the bin/ directory or specify exe_path."
        )