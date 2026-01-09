"""Job management for Multiwfn execution."""

import os
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Callable, List, Optional, Union

from .config import MultiwfnConfig, MultiwfnError
from .result import MultiwfnResult


class MultiwfnJobBuilder:
    """
    Builder for constructing MultiwfnJob instances.
    
    Provides a fluent interface for configuring and building jobs.
    
    Examples
    --------
    >>> job = (MultiwfnJobBuilder("molecule.wfn")
    ...        .with_working_dir(Path("/output"))
    ...        .with_timeout(300)
    ...        .with_menu_sequence(menu.hirshfeld_charge)
    ...        .with_menu_sequence(menu.mayer_bond_order)
    ...        .build())
    """
    
    __slots__ = ('_input_file', '_commands', '_working_dir', '_timeout', '_exe_path', '_verbose')
    
    def __init__(self, input_file: Union[str, Path]):
        """
        Initialize the builder.
        
        Parameters
        ----------
        input_file : str or Path
            Path to wavefunction file
        """
        self._input_file = Path(input_file)
        if not self._input_file.exists():
            raise FileNotFoundError(f"Input file not found: {input_file}")
        
        self._commands: List[str] = []
        self._working_dir: Optional[Path] = None
        self._timeout: Optional[int] = None
        self._exe_path: Optional[Path] = None
        self._verbose: bool = False
    
    @property
    def input_file(self) -> Path:
        """Get the input file path."""
        return self._input_file
    
    @property
    def timeout(self) -> Optional[int]:
        """Get the timeout value."""
        return self._timeout
    
    @timeout.setter
    def timeout(self, value: Optional[int]) -> None:
        """Set timeout with validation."""
        if value is not None and (not isinstance(value, int) or value <= 0):
            raise ValueError("Timeout must be a positive integer")
        self._timeout = value
    
    def with_working_dir(self, path: Union[str, Path]) -> 'MultiwfnJobBuilder':
        """Set the working directory for output files."""
        self._working_dir = Path(path)
        return self
    
    def with_timeout(self, seconds: int) -> 'MultiwfnJobBuilder':
        """Set execution timeout in seconds."""
        self.timeout = seconds
        return self
    
    def with_executable(self, path: Union[str, Path]) -> 'MultiwfnJobBuilder':
        """Set the path to Multiwfn executable."""
        self._exe_path = Path(path)
        return self
    
    def with_verbose(self, verbose: bool = True) -> 'MultiwfnJobBuilder':
        """Enable or disable verbose output."""
        self._verbose = verbose
        return self
    
    def with_config(self, config: MultiwfnConfig) -> 'MultiwfnJobBuilder':
        """
        Set configuration from a MultiwfnConfig object.
        
        Parameters
        ----------
        config : MultiwfnConfig
            Configuration object
        """
        self._exe_path = config.exe_path
        self._working_dir = config.working_dir
        self._timeout = config.timeout
        self._verbose = config.verbose
        return self
    
    def with_menu_sequence(self, menu_func: Callable, **kwargs) -> 'MultiwfnJobBuilder':
        """
        Add a menu sequence from a menu function.
        
        Parameters
        ----------
        menu_func : callable
            Menu function from menu.py
        **kwargs
            Arguments to pass to the menu function
        """
        sequence = menu_func(**kwargs) if kwargs else menu_func()
        if sequence and sequence[-1] == 'q':
            sequence = sequence[:-1]
        self._commands.extend(sequence)
        return self
    
    def with_custom_commands(self, commands: List[str]) -> 'MultiwfnJobBuilder':
        """
        Add custom command sequence.
        
        Parameters
        ----------
        commands : list of str
            List of commands to add
        """
        self._commands.extend(commands)
        return self
    
    def build(self) -> 'MultiwfnJob':
        """
        Build the MultiwfnJob instance.
        
        Returns
        -------
        MultiwfnJob
            Configured job ready for execution
        """
        config = MultiwfnConfig(
            exe_path=self._exe_path,
            working_dir=self._working_dir or Path.cwd(),
            timeout=self._timeout,
            verbose=self._verbose
        )
        
        job = MultiwfnJob(self._input_file, config=config)
        job._commands = self._commands.copy()
        return job
    
    def __repr__(self) -> str:
        return (
            f"MultiwfnJobBuilder("
            f"input_file={self._input_file!r}, "
            f"commands={len(self._commands)}, "
            f"working_dir={self._working_dir!r}, "
            f"timeout={self._timeout}, "
            f"verbose={self._verbose})"
        )


class MultiwfnJob:
    """
    Encapsulates a Multiwfn job with input file, menu sequence, and execution.
    
    Parameters
    ----------
    input_file : str or Path
        Path to wavefunction file
    config : MultiwfnConfig, optional
        Configuration object
        
    Examples
    --------
    >>> job = MultiwfnJob("molecule.wfn")
    >>> job.add_menu_sequence(menu.hirshfeld_charge)
    >>> result = job.run()
    >>> charges = result.parse_charges()
    """
    
    __slots__ = ('_input_file', '_config', '_commands', '_result')
    
    def __init__(
        self,
        input_file: Union[str, Path],
        config: Optional[MultiwfnConfig] = None
    ):
        self._input_file = Path(input_file).resolve()  # Make absolute
        if not self._input_file.exists():
            raise FileNotFoundError(f"Input file not found: {input_file}")
        
        self._config = config or MultiwfnConfig()
        self._commands: List[str] = []
        self._result: Optional[MultiwfnResult] = None
    
    @classmethod
    def builder(cls, input_file: Union[str, Path]) -> MultiwfnJobBuilder:
        """
        Create a builder for this job.
        
        Parameters
        ----------
        input_file : str or Path
            Path to wavefunction file
            
        Returns
        -------
        MultiwfnJobBuilder
            Builder instance
        """
        return MultiwfnJobBuilder(input_file)
    
    @property
    def input_file(self) -> Path:
        """Get the input file path (read-only)."""
        return self._input_file
    
    @property
    def config(self) -> MultiwfnConfig:
        """Get the job configuration (read-only)."""
        return self._config
    
    @property
    def commands(self) -> List[str]:
        """Get a copy of the current command sequence."""
        return self._commands.copy()
    
    @property
    def result(self) -> Optional[MultiwfnResult]:
        """Get the execution result, if available (read-only)."""
        return self._result
    
    @property
    def has_result(self) -> bool:
        """Check if job has been executed."""
        return self._result is not None
    
    def add_menu_sequence(self, menu_func: Callable, **kwargs) -> 'MultiwfnJob':
        """
        Add a menu sequence from a menu function.
        
        Parameters
        ----------
        menu_func : callable
            Menu function from menu.py
        **kwargs
            Arguments to pass to the menu function
            
        Returns
        -------
        self
            For method chaining
        """
        sequence = menu_func(**kwargs) if kwargs else menu_func()
        if sequence and sequence[-1] == 'q':
            sequence = sequence[:-1]
        self._commands.extend(sequence)
        return self
    
    def add_custom_commands(self, commands: List[str]) -> 'MultiwfnJob':
        """
        Add custom command sequence.
        
        Parameters
        ----------
        commands : list of str
            List of commands to add
            
        Returns
        -------
        self
            For method chaining
        """
        self._commands.extend(commands)
        return self
    
    def clear_commands(self) -> 'MultiwfnJob':
        """
        Clear all commands.
        
        Returns
        -------
        self
            For method chaining
        """
        self._commands.clear()
        return self
    
    def run(self, verbose: Optional[bool] = None, timeout: Optional[int] = None) -> MultiwfnResult:
        """
        Execute the Multiwfn job.
        
        Parameters
        ----------
        verbose : bool, optional
            Print stdout during execution (overrides config)
        timeout : int, optional
            Timeout in seconds (overrides config)
            
        Returns
        -------
        MultiwfnResult
            Execution results
        """
        verbose = verbose if verbose is not None else self._config.verbose
        timeout = timeout if timeout is not None else self._config.timeout
        
        exe = self._config.executable
        
        commands = self._commands.copy()
        if not commands or commands[-1] != "q":
            commands.append("q")
        
        with tempfile.NamedTemporaryFile(
            mode='w',
            delete=False,
            suffix=".inp",
            dir=self._config.working_dir,
            encoding='utf-8'
        ) as batch_file:
            batch_file.write("\n".join(commands) + "\n")
            batch_path = batch_file.name
        
        start_time = time.time()
        original_dir = Path.cwd()
        
        try:
            os.chdir(self._config.working_dir)
            
            with open(batch_path, 'r', encoding='utf-8') as batch:
                proc = subprocess.Popen(
                    [str(exe), str(self._input_file)],
                    stdin=batch,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True
                )
                
                try:
                    stdout, stderr = proc.communicate(timeout=timeout)
                    returncode = proc.returncode
                except subprocess.TimeoutExpired:
                    proc.kill()
                    stdout, stderr = proc.communicate()
                    raise MultiwfnError(f"Multiwfn execution timed out after {timeout}s")
            
            execution_time = time.time() - start_time
            
            # Decode bytes to string
            if isinstance(stdout, bytes):
                stdout = stdout.decode('utf-8', errors='replace')
            if isinstance(stderr, bytes):
                stderr = stderr.decode('utf-8', errors='replace')
            
            if verbose:
                print(stdout)
            
            self._result = MultiwfnResult(
                stdout=stdout,
                stderr=stderr,
                returncode=returncode,
                execution_time=execution_time,
                commands=commands,
                input_file=self._input_file
            )
            
            # Multiwfn often returns non-zero even on success in batch mode
            # Only fail on obvious errors (negative codes or very high codes)
            # or if stderr contains clear error messages
            error_indicators = ['error', 'fatal', 'cannot open', 'not found', 'failed to']
            has_error = any(ind in stderr.lower() for ind in error_indicators)
            
            if (returncode is not None and returncode < 0) or has_error:
                raise MultiwfnError(
                    f"Multiwfn failed with return code {returncode}\n"
                    f"STDERR: {stderr}"
                )
            
            return self._result
            
        finally:
            os.chdir(original_dir)
            try:
                os.unlink(batch_path)
            except OSError:
                pass
    
    def __str__(self) -> str:
        status = "executed" if self.has_result else "pending"
        return f"MultiwfnJob({self._input_file.name}, {len(self._commands)} commands, {status})"
    
    def __repr__(self) -> str:
        return (
            f"MultiwfnJob("
            f"input_file={self._input_file!r}, "
            f"commands={self._commands!r}, "
            f"config={self._config!r}, "
            f"has_result={self.has_result})"
        )