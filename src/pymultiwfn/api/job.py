"""Job management for Multiwfn execution."""

import subprocess
import tempfile
import time
from datetime import datetime
from pathlib import Path
from uuid import uuid4

from pymultiwfn.analysis.analysis import MultiwfnAnalysis
from pymultiwfn.analysis.result import MultiwfnResult
from pymultiwfn.api.exceptions import MultiwfnError
from pymultiwfn.api.multiwfn import Multiwfn
from pymultiwfn.enums.menu import Menu


class MultiwfnJob:
    """Encapsulates a Multiwfn job.

    This is the core class for actually interacting with Multiwfn. It connects
    the requested analyses with menu items, manages the execution of Multiwfn,
    and captures the results.

    However, this is not the intended entry point for most users. Instead,
    users should use the MultiwfnAnalysis class, which provides a friendly
    interface for job creations and analysis.

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
        analysis: MultiwfnAnalysis | None,
        multiwfn: Multiwfn | None = None,
        timeout: int | None = None,
        work_dir: Path | None = None,
        verbose: bool = False,
        cached: bool = True,
    ) -> None:
        """Initialise the Mutliwfn job.

        Parameters
        ----------
        input_file
            Path to wavefunction file (e.g., .wfn, .wfx, .fchk).

        multiwfn
            Multiwfn instance with executable configuration. If None, a default
            one will be created.

        analysis
            MultiwfnAnalysis to perform. Each analysis entry will be
            converted into a sequence of menu commands. If None, an empty job
            will be created for manual command addition.

        timeout
            Optional timeout in seconds for the Multiwfn execution. If None,
            there will be noe timeout, which might lead to hanging for complex
            analysed (e.g., elaborate cube generation).

        work_dir
            Optional working directory for execution. If None, a temporary
            location will be used in the current directory.

        verbose
            If True, print Multiwfn stdout during execution. Defaults to False.

        Notes
        -----
        While not recommended, it is possible to create an empty MultiwfnJob
        without going through setting up any MultiwfnAnalysis. This allows for
        manual interaction with the Multiwfn executable, but can lead to errors
        if wrong menu commands are added.

        """
        if not Path(input_file).exists():
            raise FileNotFoundError(f"Input file not found: {input_file}")
        self._input_file = Path(input_file).resolve()

        self._multiwfn = multiwfn if multiwfn is not None else Multiwfn()
        self._commands: list[str] = []
        self._result: MultiwfnResult | None = None
        self._timeout = self._validate_timeout(timeout)

        if work_dir is None:
            today_str = datetime.today().strftime("%Y-%m-%d")
            today_str += f"_{uuid4()}"
            work_dir = Path.cwd() / f"{today_str}"

        self._work_dir = work_dir.resolve()
        self._verbose = verbose
        self._cached = cached

        if analysis is not None:
            self._analysis = analysis
            for menu_item in analysis.analyses:
                self.add_menu(menu_item)

        else:
            self._analysis = MultiwfnAnalysis(
                input_file=input_file,
                analyses=[],
                cached=cached,
            )

    @classmethod
    def from_analysis(
        cls,
        analysis: MultiwfnAnalysis,
        multiwfn: Multiwfn | None = None,
        timeout: int | None = None,
        work_dir: Path | None = None,
        verbose: bool = False,
        cached: bool | None = True,
    ) -> "MultiwfnJob":
        """Create a MultiwfnJob directly from a MultiwfnAnalysis.

        Parameters
        ----------
        analysis
            MultiwfnAnalysis to perform.

        multiwfn
            Multiwfn instance with executable configuration. If None, a default
            one will be created.

        timeout
            Optional timeout in seconds for the Multiwfn execution. If None,
            there will be noe timeout, which might lead to hanging for complex
            analysed (e.g., elaborate cube generation).

        work_dir
            Optional working directory for execution. If None, a temporary
            location will be used in the current directory.

        verbose
            If True, print Multiwfn stdout during execution. Defaults to False.

        Return
        ------
        A MultiwfnJob instance ready to be executed, with menu commands
        generated from the analysis configuration.

        Notes
        -----
        This is the intended entry point for most users. The input file and
        menu sequences will be deduced from the provided MultiwfnAnalysis,
        which will subsequently be updated with the results.

        """
        return cls(
            input_file=analysis.input_file,
            analysis=analysis,
            multiwfn=multiwfn,
            timeout=timeout,
            work_dir=work_dir,
            verbose=verbose,
            cached=cached if cached is not None else analysis.cached,
        )

    @classmethod
    def from_file(
        cls,
        input_file: str | Path,
        analyses: list[Menu] | None = None,
        multiwfn: Multiwfn | None = None,
        timeout: int | None = None,
        work_dir: Path | None = None,
        verbose: bool = False,
        cached: bool = True,
    ) -> "MultiwfnJob":
        """Create a MultiwfnJob directly from an input file.

        Parameters
        ----------
        input_file
            Path to wavefunction file (e.g., .wfn, .wfx, .fchk).

        analyses
            A list of Menu enum members representing the analyses to perform.
            If None, no analyses will be added and an empty job will be created
            for manual command addition.

        multiwfn
            Multiwfn instance with executable configuration. If None, a default
            one will be created.

        timeout
            Optional timeout in seconds for the Multiwfn execution. If None,
            there will be noe timeout, which might lead to hanging for complex
            analysed (e.g., elaborate cube generation).

        work_dir
            Optional working directory for execution. If None, a temporary
            location will be used in the current directory.

        verbose
            If True, print Multiwfn stdout during execution. Defaults to False.

        Return
        ------
        A MultiwfnJob instance ready to be executed, with menu commands
        generated from the analysis configuration.

        Notes
        -----
        This is *not* the intended entry point for most users. The input file
        and menu sequences have to be provided manually; However, this method
        can be useful for more direct integration with other software.

        """
        analysis = MultiwfnAnalysis(
            input_file,
            analyses,
            cached=cached,
        )
        return cls(
            input_file=input_file,
            analysis=analysis,
            multiwfn=multiwfn,
            timeout=timeout,
            work_dir=work_dir,
            verbose=verbose,
        )

    def _validate_timeout(self, value: int | None) -> int | None:
        """Validate and set the timeout value."""
        if value is not None and value <= 0:
            raise ValueError("Timeout must be a positive integer or None")
        return value

    @property
    def timeout(self) -> int | None:
        """Get the timeout value."""
        return self._timeout

    @timeout.setter
    def timeout(self, value: int | None) -> None:
        """Set the timeout value."""
        self._timeout = self._validate_timeout(value)

    @property
    def input_file(self) -> Path:
        """Get the input file path (read-only)."""
        return self._input_file

    @property
    def multiwfn(self) -> Multiwfn:
        """Get the Multiwfn configuration."""
        return self._multiwfn

    @multiwfn.setter
    def multiwfn(self, value: Multiwfn) -> None:
        """Set the Multiwfn configuration."""
        if not isinstance(value, Multiwfn):
            raise ValueError("multiwfn must be an instance of Multiwfn")
        self._multiwfn = value

    @property
    def verbose(self) -> bool:
        """Get the verbosity setting."""
        return self._verbose

    @verbose.setter
    def verbose(self, value: bool) -> None:
        """Set the verbosity setting."""
        self._verbose = value

    @property
    def cached(self) -> bool:
        """Get the cache setting."""
        return self._cached

    @cached.setter
    def cached(self, value: bool) -> None:
        """Set the verbosity setting."""
        self._verbose = value
        self._analysis.cached = value

    @property
    def work_dir(self) -> Path:
        """Get the working directory."""
        return self._work_dir

    @work_dir.setter
    def work_dir(self, value: Path) -> None:
        """Set the working directory."""
        if not isinstance(value, Path):
            raise ValueError("work_dir must be a Path object")
        self._work_dir = value.resolve()

    @property
    def commands(self) -> list[str]:
        """Get a copy of the current command sequence."""
        return self._commands.copy()

    def add_menu(self, menu_item: Menu) -> None:
        """Add a menu sequence from a Menu enum member.

        Parameters
        ----------
        menu_item
            Menu enum member

        """
        sequence = menu_item.get_sequence()
        if sequence and sequence[-1] == "q":
            sequence = sequence[:-1]
        self._commands.extend(sequence)

    def add_custom_commands(self, commands: list[str]) -> None:
        """Add custom command sequence.

        Parameters
        ----------
        commands
            List of commands to add

        """
        self._commands.extend(commands)

    def run(
        self,
    ) -> MultiwfnResult:
        """Execute the Multiwfn job.

        Returns
        -------
        MultiwfnResult
            Execution results

        Raises
        ------
        MultiwfnError
            If execution times out or fails with errors

        """
        # Ensure we start with a newline for Multiwfn input
        commands = ["0", ""] + self._commands
        if commands[-1] != "q":
            commands.append("q")

        if not self.work_dir.exists():
            self.work_dir.mkdir(parents=True, exist_ok=True)

        with tempfile.NamedTemporaryFile(
            mode="w",
            newline="\n",
            encoding="utf-8",
            delete=False,
            suffix=".inp",
            dir=self.work_dir,
        ) as batch_file:
            batch_file.write("\n".join(commands) + "\n")

            start_time = time.time()

            proc = subprocess.Popen(
                [str(self.multiwfn.exe_path), str(self.input_file)],
                stdin=batch_file,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                cwd=self.work_dir,
            )

            try:
                stdout, stderr = proc.communicate(timeout=self.timeout)
                returncode = proc.returncode

            except subprocess.TimeoutExpired as err:
                proc.kill()
                stdout, stderr = proc.communicate()
                raise MultiwfnError(
                    "Multiwfn execution timed out after "
                    f"{self._timeout}s. "
                    "Consider increasing timeout or using "
                    "TimeoutConfig for complex analyses."
                ) from err

            execution_time = time.time() - start_time

            # Decode bytes to string
            if isinstance(stdout, bytes):
                stdout = stdout.decode("utf-8", errors="replace")
            if isinstance(stderr, bytes):
                stderr = stderr.decode("utf-8", errors="replace")

            if self.verbose:
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

    def __str__(self) -> str:
        return f"MultiwfnJob on {self._input_file.name}."

    def __repr__(self) -> str:
        return (
            f"MultiwfnJob("
            f"input_file={self._input_file!r}, "
            f"analysis={self._analysis!r}, "
            f"multiwfn={self._multiwfn!r}, "
            f"timeout={self._timeout!r}, "
            f"work_dir={self._work_dir!r}, "
            f"verbose={self._verbose!r})"
        )
