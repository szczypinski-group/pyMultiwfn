"""Module to plan and execute Multiwfn analyses."""

import logging
from pathlib import Path

from pymultiwfn.analysis.result import MultiwfnResult, ResultStore
from pymultiwfn.api.exceptions import MultiwfnError
from pymultiwfn.api.job import MultiwfnJob
from pymultiwfn.api.multiwfn import Multiwfn
from pymultiwfn.enums.analyses import AnalysisClasses
from pymultiwfn.enums.menu import Menu

logger = logging.getLogger(__name__)


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
        input_file: str | Path,
        analyses: Menu | list[Menu] | AnalysisClasses | None = None,
        cached: bool = True,
        json_path: Path | None = None,
    ) -> None:
        self.input_file: Path = Path(input_file)

        self.analyses: list[Menu] = []
        if analyses is not None:
            self.add_menu(analyses)
        self.results: list[MultiwfnResult] = []
        self.jobs: list[MultiwfnJob] = []
        self.cached = cached
        self._json_path: Path | None = (
            Path(json_path) if json_path is not None else None
        )
        self._store: ResultStore | None = None

    @property
    def json_path(self) -> Path | None:
        """Path to the JSON output file, or ``None`` to disable.

        When set to a :class:`~pathlib.Path`, every completed analysis
        is persisted to that file.  When ``None`` (the default),
        results are still cached in memory and retrievable via the
        internal store, but no file is created on disk.

        Setting this from ``None`` to a ``Path`` after results have
        already been collected will immediately flush all cached data
        to that file.
        """
        return self._json_path

    @json_path.setter
    def json_path(self, value: Path | None) -> None:
        self._json_path = Path(value) if value is not None else None
        # Propagate to an already-initialised store so that subsequent
        # (or deferred) writes respect the new setting immediately.
        if self._store is not None:
            self._store.json_path = self._json_path

    def _get_store(self, work_dir: Path | None = None) -> ResultStore:
        """Lazily initialise or return the per-molecule result store."""
        if self._store is None:
            wd = work_dir if work_dir is not None else Path.cwd()
            self._store = ResultStore(
                input_file=self.input_file,
                work_dir=wd,
                json_path=self._json_path,
            )
        return self._store

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
        store = self._get_store(work_dir)

        for menu in self.analyses:
            # Check the JSON store for a cached result.
            if self.cached and store.has_result(menu):
                logger.info(
                    f"Cached result found for {menu.name} analysis. "
                    "Loading stored result instead of re-running Multiwfn."
                )
                # Reconstruct a MultiwfnResult from stored data
                cached_data = store.get_result(menu)
                if cached_data is not None:
                    result = MultiwfnResult(analysis=menu)
                    # The cached entry is already parsed; attach it
                    # as a lightweight marker so downstream code sees
                    # a non-empty result list.
                    self.results.append(result)
            else:
                self._create_and_run(
                    analysis=menu,
                    multiwfn=multiwfn,
                    timeout=timeout,
                    work_dir=work_dir,
                    verbose=verbose,
                )

    def _create_and_run(
        self,
        analysis: Menu,
        multiwfn: Multiwfn | None = None,
        timeout: int | None = None,
        work_dir: Path | None = None,
        verbose: bool = False,
    ) -> None:
        """Create and run a single MultiwfnJob with automatic logging."""
        job = MultiwfnJob(
            input_file=self.input_file,
            analysis=analysis,
            multiwfn=multiwfn,
            timeout=timeout,
            work_dir=work_dir,
            verbose=verbose,
        )

        error: Exception | None = None
        try:
            job = job.run()
        except (MultiwfnError, Exception) as exc:
            error = exc

        # Parse stdout into a MultiwfnResult.
        result = MultiwfnResult(analysis=analysis)
        result.parse(job.stdout)
        self.results.append(result)
        self.jobs.append(job)

        # Persist parsed results to the per-molecule JSON file.
        store = self._get_store(work_dir)
        store.store(result)

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

    def _has_result(self, analysis: Menu) -> bool:
        """Check whether a parsed result already exists for *analysis*."""
        return any(result.analysis == analysis for result in self.results)

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
