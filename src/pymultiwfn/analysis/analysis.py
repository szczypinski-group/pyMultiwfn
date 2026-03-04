"""Module to plan and execute Multiwfn analyses."""

from pathlib import Path

from pymultiwfn.analysis.result import MultiwfnResult
from pymultiwfn.api.job import MultiwfnJob
from pymultiwfn.api.multiwfn import Multiwfn
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
        input_file: str | Path,
        analyses: Menu | list[Menu] | AnalysisClasses | None = None,
        cached: bool = True,
    ) -> None:
        self.input_file = input_file
        self.analyses: list[Menu] = []
        if analyses is not None:
            self.add_menu(analyses)
        self.results: dict[Menu, MultiwfnResult] = {}
        self.jobs: list[MultiwfnJob] = []
        self.cached = cached

    def run(
        self,
        multiwfn: Multiwfn | None = None,
        timeout: int | None = None,
        work_dir: Path | None = None,
        verbose: bool = False,
    ) -> MultiwfnResult:
        """Run Multiwfn job specified by the analysis.

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

        cached,
            If True, will used cached results for already calculated analyses;
            otherwise, previous analysis will be overwritten.

        Return
        ------
            Unparsed result of the Multiwfn calculation.

        """
        self._job = MultiwfnJob.from_analysis(
            analysis=self,
            multiwfn=multiwfn,
            timeout=timeout,
            work_dir=work_dir,
            verbose=verbose,
            cached=self.cached,
        )

        self._result = self._job.run()

        return self._result

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
