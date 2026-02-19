"""Job management for Multiwfn execution."""

import contextlib
import os
import subprocess
import tempfile
import time
from pathlib import Path

from pymultiwfn.config import (
    MultiwfnConfig,
    MultiwfnError,
    TimeoutConfig,
)
from pymultiwfn.menu import Menu
from pymultiwfn.result import MultiwfnResult


class MultiwfnJobBuilder:
    """Builder for constructing MultiwfnJob instances.

    Provides a fluent interface for configuring and building jobs.

    Examples
    --------
    >>> from pyMultiwfn.menu import Menu
    >>> job = (MultiwfnJobBuilder("molecule.wfn")
    ...        .with_working_dir(Path("/output"))
    ...        .with_timeout(300)
    ...        .with_menu(Menu.HIRSHFELD_CHARGE)
    ...        .with_menu(Menu.MAYER_BOND_ORDER)
    ...        .build())
    """

    def __init__(self, input_file: str | Path) -> None:
        """Initialize the builder.

        Parameters
        ----------
        input_file : str or Path
            Path to wavefunction file
        """
        self._input_file = Path(input_file)
        if not self._input_file.exists():
            raise FileNotFoundError(f"Input file not found: {input_file}")

        self._commands: list[str] = []
        self._working_dir: Path | None = None
        self._timeout: int | TimeoutConfig | None = None
        self._exe_path: Path | None = None
        self._verbose: bool = False
        self._analysis_names: list[
            str
        ] = []  # Track analysis names for timeout

    @property
    def input_file(self) -> Path:
        """Get the input file path."""
        return self._input_file

    @property
    def timeout(self) -> int | TimeoutConfig | None:
        """Get the timeout value."""
        return self._timeout

    @timeout.setter
    def timeout(self, value: int | TimeoutConfig | None) -> None:
        """Set timeout with validation."""
        if value is not None:
            if isinstance(value, int):
                if value <= 0:
                    raise ValueError("Timeout must be a positive integer")
            elif not isinstance(value, TimeoutConfig):
                raise TypeError("Timeout must be int, TimeoutConfig, or None")
        self._timeout = value

    def with_working_dir(self, path: str | Path) -> "MultiwfnJobBuilder":
        """Set the working directory for output files."""
        self._working_dir = Path(path)
        return self

    def with_timeout(
        self, seconds: int | TimeoutConfig
    ) -> "MultiwfnJobBuilder":
        """Set execution timeout.

        Parameters
        ----------
        seconds : int or TimeoutConfig
            Either a simple timeout in seconds, or a TimeoutConfig
            for analysis-specific timeouts
        """
        self.timeout = seconds
        return self

    def with_executable(self, path: str | Path) -> "MultiwfnJobBuilder":
        """Set the path to Multiwfn executable."""
        self._exe_path = Path(path)
        return self

    def with_verbose(self, verbose: bool = True) -> "MultiwfnJobBuilder":
        """Enable or disable verbose output."""
        self._verbose = verbose
        return self

    def with_config(self, config: MultiwfnConfig) -> "MultiwfnJobBuilder":
        """Set configuration from a MultiwfnConfig object.

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

    def with_menu(self, menu_item: Menu) -> "MultiwfnJobBuilder":
        """Add a menu sequence from a Menu enum member.

        Parameters
        ----------
        menu_item : Menu
            Menu enum member
        *args
            Additional arguments for the menu sequence
        """
        sequence = menu_item.get_sequence()
        if sequence and sequence[-1] == "q":
            sequence = sequence[:-1]
        self._commands.extend(sequence)
        self._analysis_names.append(menu_item.name.lower())
        return self

    def with_custom_commands(
        self, commands: list[str]
    ) -> "MultiwfnJobBuilder":
        """Add custom command sequence.

        Parameters
        ----------
        commands : list of str
            List of commands to add
        """
        self._commands.extend(commands)
        return self

    def build(self) -> "MultiwfnJob":
        """Build the MultiwfnJob instance.

        Returns
        -------
        MultiwfnJob
            Configured job ready for execution
        """
        config = MultiwfnConfig(
            exe_path=self._exe_path,
            working_dir=self._working_dir or Path.cwd(),
            timeout=self._timeout,
            verbose=self._verbose,
        )

        job = MultiwfnJob(self._input_file, config=config)
        job._commands = self._commands.copy()
        job._analysis_names = self._analysis_names.copy()
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
    """Encapsulates a Multiwfn job.

    Parameters
    ----------
    input_file : str or Path
        Path to wavefunction file
    config : MultiwfnConfig, optional
        Configuration object

    Examples
    --------
    >>> from pyMultiwfn import MultiwfnJob
    >>> from pyMultiwfn.menu import Menu
    >>> job = MultiwfnJob("molecule.wfn")
    >>> job.add_menu(Menu.HIRSHFELD_CHARGE)
    >>> job.add_menu(Menu.MAYER_BOND_ORDER)
    >>> result = job.run()
    >>> charges = result.parse_charges()
    """

    def __init__(
        self,
        input_file: str | Path,
        config: MultiwfnConfig | None = None,
    ) -> None:
        """Initialise the mutliwfn job."""
        self._input_file = Path(input_file).resolve()  # Make absolute
        if not self._input_file.exists():
            raise FileNotFoundError(f"Input file not found: {input_file}")

        self._config = config or MultiwfnConfig()
        self._commands: list[str] = []
        self._result: MultiwfnResult | None = None
        self._analysis_names: list[
            str
        ] = []  # Track analysis names for timeout

    @classmethod
    def builder(cls, input_file: str | Path) -> MultiwfnJobBuilder:
        """Create a builder for this job.

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
    def commands(self) -> list[str]:
        """Get a copy of the current command sequence."""
        return self._commands.copy()

    @property
    def result(self) -> MultiwfnResult | None:
        """Get the execution result, if available (read-only)."""
        return self._result

    @property
    def has_result(self) -> bool:
        """Check if job has been executed."""
        return self._result is not None

    def add_menu(self, menu_item: Menu) -> "MultiwfnJob":
        """Add a menu sequence from a Menu enum member.

        Parameters
        ----------
        menu_item : Menu
            Menu enum member
        *args
            Additional arguments for the menu sequence

        Returns
        -------
        self
            For method chaining
        """
        sequence = menu_item.get_sequence()
        if sequence and sequence[-1] == "q":
            sequence = sequence[:-1]
        self._commands.extend(sequence)
        self._analysis_names.append(menu_item.name.lower())
        return self

    def add_custom_commands(self, commands: list[str]) -> "MultiwfnJob":
        """Add custom command sequence.

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

    def clear_commands(self) -> "MultiwfnJob":
        """Clear all commands.

        Returns
        -------
        self
            For method chaining
        """
        self._commands.clear()
        self._analysis_names.clear()
        return self

    def _calculate_timeout(self, override: int | None = None) -> int:
        """Calculate the appropriate timeout for this job.

        Parameters
        ----------
        override : int, optional
            If provided, use this timeout value directly

        Returns
        -------
        int
            Timeout in seconds
        """
        if override is not None:
            return override

        # If no analysis names tracked, use default
        if not self._analysis_names:
            return self._config.get_timeout()

        # Get the maximum timeout across all analyses
        max_timeout = 0
        for name in self._analysis_names:
            timeout = self._config.get_timeout(name)
            max_timeout = max(max_timeout, timeout)

        # If multiple analyses, add buffer time
        if len(self._analysis_names) > 1:
            # Add 30 seconds per additional analysis as buffer
            buffer = (len(self._analysis_names) - 1) * 30
            max_timeout += buffer

        return max_timeout

    def run(
        self, verbose: bool | None = None, timeout: int | None = None
    ) -> MultiwfnResult:
        """Execute the Multiwfn job.

        Parameters
        ----------
        verbose : bool, optional
            Print stdout during execution (overrides config)
        timeout : int, optional
            Timeout in seconds (overrides config and automatic calculation)

        Returns
        -------
        MultiwfnResult
            Execution results

        Raises
        ------
        MultiwfnError
            If execution times out or fails with errors
        """
        verbose = verbose if verbose is not None else self._config.verbose
        effective_timeout = self._calculate_timeout(timeout)

        exe = self._config.executable

        commands = self._commands.copy()
        if not commands or commands[-1] != "q":
            commands.append("q")

        with tempfile.NamedTemporaryFile(
            mode="w",
            delete=False,
            suffix=".inp",
            dir=self._config.working_dir,
            encoding="utf-8",
        ) as batch_file:
            batch_file.write("\n".join(commands) + "\n")
            batch_path = Path(batch_file.name)

        start_time = time.time()
        original_dir = Path.cwd()

        try:
            os.chdir(self._config.working_dir)

            with batch_path.open(encoding="utf-8") as batch:
                proc = subprocess.Popen(
                    [str(exe), str(self._input_file)],
                    stdin=batch,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )

                try:
                    stdout, stderr = proc.communicate(
                        timeout=effective_timeout
                    )
                    returncode = proc.returncode
                except subprocess.TimeoutExpired as err:
                    proc.kill()
                    stdout, stderr = proc.communicate()
                    raise MultiwfnError(
                        "Multiwfn execution timed out after "
                        f"{effective_timeout}s. "
                        f"Analyses: {', '.join(self._analysis_names)}. "
                        "Consider increasing timeout or using "
                        "TimeoutConfig for complex analyses."
                    ) from err

            execution_time = time.time() - start_time

            # Decode bytes to string
            if isinstance(stdout, bytes):
                stdout = stdout.decode("utf-8", errors="replace")
            if isinstance(stderr, bytes):
                stderr = stderr.decode("utf-8", errors="replace")

            if verbose:
                print(stdout)

            self._result = MultiwfnResult(
                stdout=stdout,
                stderr=stderr,
                returncode=returncode,
                execution_time=execution_time,
                commands=commands,
                input_file=self._input_file,
            )

            # Multiwfn often returns non-zero even on success in batch mode
            # Only fail on obvious errors (negative codes or very high codes)
            # or if stderr contains clear error messages
            error_indicators = [
                "error",
                "fatal",
                "cannot open",
                "not found",
                "failed to",
            ]
            has_error = any(ind in stderr.lower() for ind in error_indicators)

            if (returncode is not None and returncode < 0) or has_error:
                raise MultiwfnError(
                    f"Multiwfn failed with return code {returncode}\n"
                    f"STDERR: {stderr}"
                )

            return self._result

        finally:
            os.chdir(original_dir)
            with contextlib.suppress(OSError):
                batch_path.unlink()

    def __str__(self) -> str:
        status = "executed" if self.has_result else "pending"
        return (
            f"MultiwfnJob({self._input_file.name}, {len(self._commands)} "
            f"commands, {status})"
        )

    def __repr__(self) -> str:
        return (
            f"MultiwfnJob("
            f"input_file={self._input_file!r}, "
            f"commands={self._commands!r}, "
            f"config={self._config!r}, "
            f"has_result={self.has_result})"
        )
