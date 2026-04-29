"""Job management for Multiwfn execution."""

import subprocess
import time
from pathlib import Path

from pymultiwfn.api.exceptions import InvalidInputFileError, MultiwfnError
from pymultiwfn.api.multiwfn import Multiwfn
from pymultiwfn.api.outcome import MultiwfnJobOutcome
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

    _SUPPORTED_INPUT_EXTENSIONS: tuple[str, ...] = (
        ".mwfn",
        ".wfn",
        ".wfx",
        ".fch",
        ".fchk",
        ".molden",
        ".31",
        ".32",
        ".33",
        ".34",
        ".35",
        ".36",
        ".37",
        ".38",
        ".39",
        ".40",
        ".pdb",
        ".xyz",
        ".chg",
        ".cub",
        ".cube",
        ".grd",
        ".mol",
        ".mol2",
        ".sdf",
        ".gro",
        ".cif",
        ".log",
        ".out",
        ".gjf",
        ".com",
        ".inp",
        ".mop",
    )

    def __init__(
        self,
        input_file: str | Path,
        analysis: Menu | None,
        multiwfn: Multiwfn | None = None,
        timeout: int | None = None,
        work_dir: Path | None = None,
        verbose: bool = False,
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
            A Menu enum members representing the analyses to perform.
            If None, no analyses will be added and an empty job will be created
            for manual command addition.

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
        input_path = Path(input_file)
        if not input_path.exists():
            raise FileNotFoundError(f"Input file not found: {input_file}")
        self._validate_input_file(input_path)
        self._input_file = input_path.resolve()

        self._multiwfn = multiwfn if multiwfn is not None else Multiwfn()
        self._commands: list[str] = []
        self._timeout = self._validate_timeout(timeout)

        if work_dir is None:
            work_dir = Path.cwd()

        self._work_dir = work_dir.resolve()
        self._verbose = verbose

        if analysis is not None:
            self._analysis = analysis
            self._parse_menu(analysis)

        # Attributes below are populated on execution.

        self._executed = False
        self._result: MultiwfnJobOutcome | None = None

    @classmethod
    def _validate_input_file(cls, input_path: Path) -> None:
        """Validate input file extension against Multiwfn-supported formats."""
        extension = input_path.suffix.lower()
        if extension not in cls._SUPPORTED_INPUT_EXTENSIONS:
            raise InvalidInputFileError(
                input_path=str(input_path),
                extension=extension,
                supported_extensions=cls._SUPPORTED_INPUT_EXTENSIONS,
            )

    @classmethod
    def from_file(
        cls,
        input_file: str | Path,
        analysis: Menu | None = None,
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

        analysis
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

    @property
    def executed(self) -> bool:
        """Get a flag whether the job has been executed (read-only)."""
        return self._executed

    @property
    def stderr(self) -> str | None:
        """Get the standard error from execution (read-only)."""
        return self._result.stderr if self._result is not None else None

    @property
    def stdout(self) -> str:
        """Get the standard output from execution (read-only)."""
        return self._result.stdout if self._result is not None else ""

    @property
    def return_code(self) -> int | None:
        """Get the return code from execution (read-only)."""
        return self._result.return_code if self._result is not None else None

    @property
    def execution_time(self) -> float | None:
        """Get the execution time (read-only)."""
        return (
            self._result.execution_time if self._result is not None else None
        )

    @property
    def success(self) -> bool | None:
        if self._result is not None:
            return self._result.is_successful()
        else:
            return None

    def _parse_menu(self, menu_item: Menu) -> None:
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

    def add_commands(self, commands: list[str]) -> None:
        """Add custom command sequence.

        Parameters
        ----------
        commands
            List of commands to add

        """
        self._commands.extend(commands)

    # add argument for input file name
    def run(
        self,
    ) -> "MultiwfnJob":
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
        commands = [
            "0",
            "",
        ] + self._commands
        if commands[-1] != "q":
            commands.append("q")

        self.work_dir.mkdir(parents=True, exist_ok=True)

        start_time = time.time()

        try:
            proc = subprocess.run(
                [str(self.multiwfn.exe_path), str(self.input_file)],
                input="\n".join(commands),
                capture_output=True,
                text=True,
                cwd=self.work_dir,
                timeout=self.timeout,
            )

        except subprocess.TimeoutExpired as e:
            raise MultiwfnError(
                "Multiwfn execution timed out after "
                f"{self._timeout}s. "
                "Consider increasing timeout or using "
                "setting to None."
            ) from e

        execution_time = time.time() - start_time

        if self.verbose:
            print(proc.stdout)

        result = MultiwfnJobOutcome(
            stderr=proc.stderr,
            stdout=proc.stdout,
            return_code=proc.returncode,
            execution_time=execution_time,
        )

        self._result = result
        self._executed = True

        return self

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
