"""Module to plan and execute Multiwfn analyses."""

from pathlib import Path

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

    Examples
    --------
    >>> # Example of a specific analysis inheriting from MultiwfnAnalysis
    >>> class ElectronDensityAnalysis(MultiwfnAnalysis):
    ...     def get_menu_sequence(self) -> list[str]:
    ...         return ["1", "2", "3"]

    """

    def __init__(
        self,
        input_file: str | Path,
        analyses: list[Menu] | None = None,
    ) -> None:
        self.input_file = input_file
        self.analyses = analyses if analyses is not None else []

    def _add_charges_menus(self) -> None:
        """Run all charge analyses."""
        self.analyses.extend(AnalysisClasses.CHARGES.value)

    def _add_bond_orders_menus(self) -> None:
        """Run all bond order analyses."""
        self.analyses.extend(AnalysisClasses.BOND_ORDERS.value)

    def _add_topology_menus(self) -> None:
        """Run all topology analyses."""
        self.analyses.extend(AnalysisClasses.TOPOLOGY.value)

    def _add_weak_interactions_menus(self) -> None:
        """Run all weak interaction analyses."""
        self.analyses.extend(AnalysisClasses.WEAK_INTERACTIONS.value)

    def _add_spectra_menus(self) -> None:
        """Run all spectrum analyses."""
        self.analyses.extend(AnalysisClasses.SPECTRA.value)

    def _add_surfaces_menus(self) -> None:
        """Run all surface analyses."""
        self.analyses.extend(AnalysisClasses.SURFACES.value)

    def _add_aromaticity_menus(self) -> None:
        """Run all aromaticity analyses."""
        self.analyses.extend(AnalysisClasses.AROMATICITY.value)

    def _add_cdft_menus(self) -> None:
        """Run all CDFT analyses."""
        self.analyses.extend(AnalysisClasses.CDFT.value)

    def _add_dos_menus(self) -> None:
        """Run all density of states analyses."""
        self.analyses.extend(AnalysisClasses.DOS.value)

    def _add_basin_menus(self) -> None:
        """Run all basin analyses."""
        self.analyses.extend(AnalysisClasses.BASIN.value)

    def _add_excitation_menus(self) -> None:
        """Run all electron excitation analyses."""
        self.analyses.extend(AnalysisClasses.EXCITATION.value)

    def _add_cubes_menus(self) -> None:
        """Run all cube generation analyses."""
        self.analyses.extend(AnalysisClasses.CUBES.value)

    def _add_orbital_composition_menus(self) -> None:
        """Run all orbital composition analyses."""
        self.analyses.extend(AnalysisClasses.ORBITAL_COMPOSITION.value)

    def _add_orbital_localization_menus(self) -> None:
        """Run all orbital localization analyses."""
        self.analyses.extend(AnalysisClasses.ORBITAL_LOCALIZATION.value)

    def _add_fuzzy_space_menus(self) -> None:
        """Run all fuzzy atomic space analyses."""
        self.analyses.extend(AnalysisClasses.FUZZY_SPACE.value)

    def _add_eda_menus(self) -> None:
        """Run all energy decomposition analyses."""
        self.analyses.extend(AnalysisClasses.EDA.value)

    def _add_polarizability_menus(self) -> None:
        """Run all polarizability analyses."""
        self.analyses.extend(AnalysisClasses.POLARIZABILITY.value)

    def _add_wavefunction_menus(self) -> None:
        """Run all wavefunction check/modify analyses."""
        self.analyses.extend(AnalysisClasses.WAVEFUNCTION.value)

    def _add_line_plots_menus(self) -> None:
        """Run all line property plot analyses."""
        self.analyses.extend(AnalysisClasses.LINE_PLOTS.value)

    def _add_plane_maps_menus(self) -> None:
        """Run all plane property map analyses."""
        self.analyses.extend(AnalysisClasses.PLANE_MAPS.value)
