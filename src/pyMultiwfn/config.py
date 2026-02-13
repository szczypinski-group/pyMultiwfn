"""Configuration management for pyMultiwfn."""

import platform
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Union


class MultiwfnError(Exception):
    """Base exception for pyMultiwfn errors."""
    pass


class TimeoutConfig:
    """
    Timeout configuration with sensible defaults for different analysis types.
    
    Attributes
    ----------
    default : int
        Default timeout for simple analyses (seconds)
    topology : int
        Timeout for topology analyses (longer due to CP search)
    batch : int
        Timeout for batch/multiple analyses
    cube : int
        Timeout for cube file generation
    
    Examples
    --------
    >>> timeout = TimeoutConfig(default=120, topology=600)
    >>> timeout.get_for_analysis('topology_search_cps')
    600
    """
    
    # Analysis types that require longer timeouts
    TOPOLOGY_ANALYSES = frozenset({
        'topology_search_cps', 'topology_generate_paths', 
        'topology_interbasin_surfaces', 'topology_analysis_complete',
        'topology_esp_analysis', 'topology_lol_analysis',
        'complete_qtaim_analysis', 'basin_analysis_aim', 'basin_analysis_elf'
    })
    
    CUBE_ANALYSES = frozenset({
        'cube_density', 'cube_spin_density', 'cube_elf', 'cube_lol',
        'cube_esp', 'cube_laplacian', 'cube_fukui_minus', 'cube_fukui_plus',
        'cube_dual_descriptor', 'cube_orbital', 'batch_cube_generation'
    })
    
    BATCH_ANALYSES = frozenset({
        'run_all', 'run_all_charges', 'run_all_bond_orders',
        'run_all_weak_interactions', 'run_all_topology', 'run_all_cubes',
        'run_all_spectra', 'run_all_surfaces', 'run_all_aromaticity',
        'run_all_cdft', 'weak_interaction_suite'
    })
    
    def __init__(
        self,
        default: int = 120,
        topology: int = 300,
        batch: int = 600,
        cube: int = 180
    ):
        """
        Initialize timeout configuration.
        
        Parameters
        ----------
        default : int
            Default timeout in seconds (default: 120)
        topology : int
            Timeout for topology analyses in seconds (default: 300)
        batch : int
            Timeout for batch analyses in seconds (default: 600)
        cube : int
            Timeout for cube generation in seconds (default: 180)
        """
        self._validate_timeout(default, "default")
        self._validate_timeout(topology, "topology")
        self._validate_timeout(batch, "batch")
        self._validate_timeout(cube, "cube")
        
        self.default = default
        self.topology = topology
        self.batch = batch
        self.cube = cube
    
    @staticmethod
    def _validate_timeout(value: int, name: str) -> None:
        """Validate a timeout value."""
        if not isinstance(value, int):
            raise TypeError(f"{name} timeout must be an integer, got {type(value).__name__}")
        if value <= 0:
            raise ValueError(f"{name} timeout must be positive, got {value}")
    
    def get_for_analysis(self, analysis_name: str) -> int:
        """
        Get the appropriate timeout for a given analysis type.
        
        Parameters
        ----------
        analysis_name : str
            Name of the analysis function
            
        Returns
        -------
        int
            Timeout in seconds
        """
        if analysis_name in self.BATCH_ANALYSES:
            return self.batch
        elif analysis_name in self.TOPOLOGY_ANALYSES:
            return self.topology
        elif analysis_name in self.CUBE_ANALYSES:
            return self.cube
        else:
            return self.default
    
    def __repr__(self) -> str:
        return (
            f"TimeoutConfig(default={self.default}, topology={self.topology}, "
            f"batch={self.batch}, cube={self.cube})"
        )


@dataclass
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
    exe_path: Optional[Path] = None
    working_dir: Path = field(default_factory=Path.cwd)
    timeout: Optional[Union[int, TimeoutConfig]] = None
    verbose: bool = False
    
    _resolved_exe: Optional[Path] = field(default=None, init=False, repr=False)
    _timeout_config: Optional[TimeoutConfig] = field(default=None, init=False, repr=False)
    
    def __post_init__(self):
        """Validate and normalize configuration after initialization."""
        # Normalize working_dir to Path
        if isinstance(self.working_dir, str):
            object.__setattr__(self, 'working_dir', Path(self.working_dir))
        
        # Normalize exe_path to Path
        if isinstance(self.exe_path, str):
            object.__setattr__(self, 'exe_path', Path(self.exe_path))
        
        # Set up timeout configuration
        if self.timeout is None:
            self._timeout_config = TimeoutConfig()
        elif isinstance(self.timeout, int):
            if self.timeout <= 0:
                raise ValueError("Timeout must be positive")
            self._timeout_config = TimeoutConfig(
                default=self.timeout,
                topology=self.timeout,
                batch=self.timeout * 3,
                cube=self.timeout
            )
        elif isinstance(self.timeout, TimeoutConfig):
            self._timeout_config = self.timeout
        else:
            raise TypeError(
                f"timeout must be int, TimeoutConfig, or None, got {type(self.timeout).__name__}"
            )
    
    @property
    def timeout_config(self) -> TimeoutConfig:
        """Get the timeout configuration."""
        if self._timeout_config is None:
            self._timeout_config = TimeoutConfig()
        return self._timeout_config
    
    def get_timeout(self, analysis_name: Optional[str] = None) -> int:
        """
        Get timeout for a specific analysis type.
        
        Parameters
        ----------
        analysis_name : str, optional
            Name of the analysis. If None, returns default timeout.
            
        Returns
        -------
        int
            Timeout in seconds
        """
        if analysis_name is None:
            return self.timeout_config.default
        return self.timeout_config.get_for_analysis(analysis_name)
    
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