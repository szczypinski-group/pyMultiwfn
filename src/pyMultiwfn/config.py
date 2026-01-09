"""Configuration management for pyMultiwfn."""

import platform
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


class MultiwfnError(Exception):
    """Base exception for pyMultiwfn errors."""
    pass


@dataclass(slots=True)
class MultiwfnConfig:
    """
    Configuration and executable discovery for Multiwfn.
    
    Attributes
    ----------
    exe_path : Path, optional
        Explicit path to Multiwfn executable
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
        """Find Multiwfn executable in common locations."""
        # Check explicit path
        if self.exe_path:
            if Path(self.exe_path).exists():
                return Path(self.exe_path)
            raise MultiwfnError(f"Specified executable not found: {self.exe_path}")
        
        # Check PATH
        if (found := shutil.which("Multiwfn")):
            return Path(found)
        
        is_windows = platform.system() == "Windows"
        exe_name = "Multiwfn.exe" if is_windows else "Multiwfn"
        
        # Check project bin/ directory (walk up to 5 levels)
        current = Path(__file__).resolve().parent
        for _ in range(5):
            bin_exe = current / "bin" / exe_name
            if bin_exe.exists():
                return bin_exe
            current = current.parent
        
        # Check common system locations
        candidates = [
            Path.home() / "Multiwfn" / exe_name,
            Path(f"C:/Multiwfn/{exe_name}") if is_windows else Path("/usr/local/bin/Multiwfn"),
            Path(f"C:/Program Files/Multiwfn/{exe_name}") if is_windows else Path("/opt/Multiwfn/Multiwfn"),
        ]
        
        for path in candidates:
            if path.exists() and path.is_file():
                return path
        
        raise MultiwfnError(
            "Multiwfn executable not found. Please specify path with exe_path parameter."
        )