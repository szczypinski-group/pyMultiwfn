"""Module to plan and execute Multiwfn analyses."""

from pathlib import Path

from pymultiwfn.analysis.result import MultiwfnResult
from pymultiwfn.api.exceptions import MultiwfnError
from pymultiwfn.api.job import MultiwfnJob
from pymultiwfn.api.logging import BatchLogger
from pymultiwfn.api.multiwfn import Multiwfn
from pymultiwfn.api.result_store import ResultStore
from pymultiwfn.enums.analyses import AnalysisClasses
from pymultiwfn.enums.menu import Menu


class MultiwfnAnalysis:
    """Base class for Multiwfn analyses.

    This class serves as a template for specific analyses. Each analysis type
    should inherit from this base class and implement the get_menu_sequence()
    method to provide the appropriate menu commands for that analysis.

    Attributes
    ----------
    input_file : str or Path
        Path to the wavefunction file to be analyzed.

    """

    def __init__(
        self,
        input_file: str | Path | list[str | Path],
        analyses: Menu | list[Menu] | AnalysisClasses | None = None,
        cached: bool = True,
    ) -> None:
        # Accept a single file or a list for batch runs.
        if isinstance(input_file, (str, Path)):
            self.input_files: list[Path] = [Path(input_file)]
        else:
            self.input_files = [Path(f) for f in input_file]

        # Keep a single-file alias for backwards compatibility.
        self.input_file: Path = self.input_files[0]

        self.analyses: list[Menu] = []
        if analyses is not None:
            self.add_menu(analyses)
        self.results: dict[Menu, MultiwfnResult] = {}
        self.jobs: list[MultiwfnJob] = []
        self.cached = cached

        # Populated after run() completes.
        self._logger: BatchLogger | None = None
        self._stores: dict[Path, ResultStore] = {}

    @property
    def log_path(self) -> Path | None:
        """Path to the batch log file, available after :meth:`run`."""
        if self._logger is not None:
            return self._logger.log_path
        return None

    def get_store(self, input_file: Path | str) -> ResultStore | None:
        """Return the :class:`ResultStore` for *input_file*, if available."""
        return self._stores.get(Path(input_file))

    def run(
        self,
        multiwfn: Multiwfn | None = None,
        timeout: int | None = None,
        work_dir: Path | None = None,
        verbose: bool = False,
    ) -> None:
        """Run all queued Multiwfn analyses across all input files.

        A batch log file is automatically written into *work_dir* (or the
        current directory if *work_dir* is None). Parsed results are
        accumulated into a per-molecule ``.json`` file in the same
        directory. After completion the log path is accessible via
        :attr:`log_path`.

        Parameters
        ----------
        multiwfn
            Multiwfn instance with executable configuration. If None, a
            default one will be created.

        timeout
            Optional timeout in seconds for each Multiwfn execution.

        work_dir
            Optional working directory for execution. If None, a temporary
            location will be used in the current directory.

        verbose
            If True, print Multiwfn stdout during execution.

        """
        log_dir = work_dir if work_dir is not None else Path.cwd()

        logger = BatchLogger(log_dir=log_dir)
        self._logger = logger

        logger.start_batch(
            files=self.input_files,
            analyses=self.analyses,
        )

        # Create a ResultStore per input file.
        for input_file in self.input_files:
            store = ResultStore(
                input_file=input_file,
                work_dir=log_dir,
            )
            self._stores[input_file] = store

        try:
            for input_file in self.input_files:
                store = self._stores[input_file]
                for menu in self.analyses:
                    # Check the JSON store for a cached result.
                    if self.cached and store.has_result(menu):
                        logger.log_job_skipped(
                            input_file=str(input_file),
                            analysis_name=menu.name,
                            reason="cached result exists in JSON store",
                        )
                    else:
                        self._create_and_run(
                            input_file=input_file,
                            analysis=menu,
                            multiwfn=multiwfn,
                            timeout=timeout,
                            work_dir=work_dir,
                            verbose=verbose,
                            logger=logger,
                            store=store,
                        )
        finally:
            logger.end_batch()

    def _create_and_run(
        self,
        input_file: Path,
        analysis: Menu,
        multiwfn: Multiwfn | None = None,
        timeout: int | None = None,
        work_dir: Path | None = None,
        verbose: bool = False,
        logger: BatchLogger | None = None,
        store: ResultStore | None = None,
    ) -> None:
        """Create and run a single MultiwfnJob with automatic logging."""
        job = MultiwfnJob(
            input_file=input_file,
            analysis=analysis,
            multiwfn=multiwfn,
            timeout=timeout,
            work_dir=work_dir,
            verbose=verbose,
        )

        log_entry = None
        if logger is not None:
            log_entry = logger.log_job_start(
                input_file=str(input_file),
                analysis_name=analysis.name,
                commands=job.commands,
            )

        error: Exception | None = None
        try:
            job = job.run()
        except (MultiwfnError, Exception) as exc:
            error = exc

        self.jobs.append(job)

        # Log the outcome.
        if logger is not None and log_entry is not None:
            logger.log_job_end(
                entry=log_entry,
                outcome=job._result,
                error=error,
            )

        # Parse and persist the result into the per-molecule JSON.
        if error is None and store is not None and job.stdout is not None:
            store.store_result(analysis=analysis, stdout=job.stdout)

        if error is not None:
            raise error

    def add_menu(
        self,
        menu: Menu | list[Menu] | AnalysisClasses,
    ) -> None:
        """Add a Menu enum member to the analysis.

        Parameters
        ----------
        menu
            Menu enum member or list of Menu members.

        """
        if isinstance(menu, list):
            self.analyses.extend(menu)
        elif isinstance(menu, Menu):
            self.analyses.append(menu)
        elif isinstance(menu, AnalysisClasses):
            self.add_menu(menu.value)
        else:
            raise TypeError(
                "Menu has to be a valid Menu enum or a list of Menu enums."
            )

    #######################################################################
    ### Convenience methods for bulk adding pre-defined AnalysisClasses ###
    #######################################################################

    def _add_charges_menus(self) -> None:
        """Run all charge analyses."""
        self.add_menu(AnalysisClasses.CHARGES)

    def _add_bond_orders_menus(self) -> None:
        """Run all bond order analyses."""
        self.add_menu(AnalysisClasses.BOND_ORDERS)

    def _add_topology_menus(self) -> None:
        """Run all topology analyses."""
        self.add_menu(AnalysisClasses.TOPOLOGY)

    def _add_weak_interactions_menus(self) -> None:
        """Run all weak interaction analyses."""
        self.add_menu(AnalysisClasses.WEAK_INTERACTIONS)

    def _add_spectra_menus(self) -> None:
        """Run all spectrum analyses."""
        self.add_menu(AnalysisClasses.SPECTRA)

    def _add_surfaces_menus(self) -> None:
        """Run all surface analyses."""
        self.add_menu(AnalysisClasses.SURFACES)

    def _add_aromaticity_menus(self) -> None:
        """Run all aromaticity analyses."""
        self.add_menu(AnalysisClasses.AROMATICITY)

    def _add_cdft_menus(self) -> None:
        """Run all CDFT analyses."""
        self.add_menu(AnalysisClasses.CDFT)

    def _add_dos_menus(self) -> None:
        """Run all density of states analyses."""
        self.add_menu(AnalysisClasses.DOS)

    def _add_basin_menus(self) -> None:
        """Run all basin analyses."""
        self.add_menu(AnalysisClasses.BASIN)

    def _add_excitation_menus(self) -> None:
        """Run all electron excitation analyses."""
        self.add_menu(AnalysisClasses.EXCITATION)

    def _add_cubes_menus(self) -> None:
        """Run all cube generation analyses."""
        self.add_menu(AnalysisClasses.CUBES)

    def _add_orbital_composition_menus(self) -> None:
        """Run all orbital composition analyses."""
        self.add_menu(AnalysisClasses.ORBITAL_COMPOSITION)

    def _add_orbital_localization_menus(self) -> None:
        """Run all orbital localization analyses."""
        self.add_menu(AnalysisClasses.ORBITAL_LOCALIZATION)

    def _add_fuzzy_space_menus(self) -> None:
        """Run all fuzzy atomic space analyses."""
        self.add_menu(AnalysisClasses.FUZZY_SPACE)

    def _add_eda_menus(self) -> None:
        """Run all energy decomposition analyses."""
        self.add_menu(AnalysisClasses.EDA)

    def _add_polarizability_menus(self) -> None:
        """Run all polarizability analyses."""
        self.add_menu(AnalysisClasses.POLARIZABILITY)

    def _add_wavefunction_menus(self) -> None:
        """Run all wavefunction check/modify analyses."""
        self.add_menu(AnalysisClasses.WAVEFUNCTION)

    def _add_line_plots_menus(self) -> None:
        """Run all line property plot analyses."""
        self.add_menu(AnalysisClasses.LINE_PLOTS)

    def _add_plane_maps_menus(self) -> None:
        """Run all plane property map analyses."""
        self.add_menu(AnalysisClasses.PLANE_MAPS)