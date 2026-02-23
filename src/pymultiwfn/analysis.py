"""Convenience API for running Multiwfn calculations.

Provides the Analysis class for running calculations without manually
building MultiwfnJob objects.

Usage
-----
>>> from pymultiwfn import Analysis, Menu
>>> analysis = Analysis("molecule.wfn")
>>> result = analysis.run(Menu.HIRSHFELD_CHARGE)
>>> charges = result.parse_charges()

>>> # Run category
>>> results = analysis.run_charges()

>>> # Run everything
>>> all_results = analysis.run_all()
"""

from pathlib import Path

from pymultiwfn.config import MultiwfnConfig
from pymultiwfn.job import MultiwfnJob
from pymultiwfn.menu import Menu
from pymultiwfn.result import MultiwfnResult


class MultiwfnAnalysis:
    """Convenience class for running Multiwfn analyses.

    Parameters
    ----------
    filepath : str or Path
        Path to wavefunction file
    config : MultiwfnConfig, optional
        Configuration for Multiwfn execution

    Examples
    --------
    >>> from pymultiwfn import Analysis, Menu
    >>> analysis = Analysis("molecule.wfn")

    >>> # Single analysis
    >>> result = analysis.run(Menu.HIRSHFELD_CHARGE)
    >>> charges = result.parse_charges()

    >>> # Multiple analyses in one session
    >>> menus = [Menu.HIRSHFELD_CHARGE, Menu.MAYER_BOND_ORDER]
    >>> result = analysis.run_multiple(menus)

    >>> # All charges
    >>> results = analysis.run_charges()
    >>> hirshfeld = results["hirshfeld_charge"].parse_charges()

    >>> # Everything
    >>> all_results = analysis.run_all()
    """

    # =========================================================================
    # Category definitions
    # =========================================================================

    CHARGES = [
        Menu.HIRSHFELD_CHARGE,
        Menu.VDD_POPULATION,
        Menu.MULLIKEN_POPULATION,
        Menu.LOWDIN_POPULATION,
        Menu.SCPA_POPULATION,
        Menu.STOUT_POLITZER_POPULATION,
        Menu.BICKELHAUPT_POPULATION,
        Menu.BECKE_CHARGE,
        Menu.ADCH_CHARGE,
        Menu.CHELPG_CHARGE,
        Menu.MK_CHARGE,
        Menu.AIM_CHARGE,
        Menu.HIRSHFELD_I_CHARGE,
        Menu.CM5_CHARGE,
        Menu.EEM_CHARGE,
        Menu.RESP_CHARGE,
        Menu.GASTEIGER_CHARGE,
        Menu.MBIS_CHARGE,
        Menu.DDEC_CHARGE,
    ]

    BOND_ORDERS = [
        Menu.MAYER_BOND_ORDER,
        Menu.MULTICENTER_BOND_ORDER,
        Menu.MULTICENTER_BOND_ORDER_NAO,
        Menu.WIBERG_BOND_ORDER,
        Menu.MULLIKEN_BOND_ORDER,
        Menu.MULLIKEN_BOND_ORDER_DECOMPOSE,
        Menu.ORBITAL_PERTURBED_MAYER,
        Menu.FUZZY_BOND_ORDER,
        Menu.LAPLACIAN_BOND_ORDER,
        Menu.WIBERG_DECOMPOSITION,
        Menu.IBSI_ANALYSIS,
        Menu.AV1245_INDEX,
    ]

    TOPOLOGY = [
        Menu.TOPOLOGY_SEARCH_CPS,
        Menu.TOPOLOGY_GENERATE_PATHS,
        Menu.TOPOLOGY_INTERBASIN_SURFACES,
        Menu.TOPOLOGY_ANALYSIS_COMPLETE,
        Menu.TOPOLOGY_ESP_ANALYSIS,
        Menu.TOPOLOGY_LOL_ANALYSIS,
        Menu.TOPOLOGY_ELF_ANALYSIS,
        Menu.TOPOLOGY_LAPLACIAN_ANALYSIS,
        Menu.TOPOLOGY_SEARCH_BCP,
        Menu.TOPOLOGY_SEARCH_RCP,
        Menu.TOPOLOGY_SEARCH_CCP,
    ]

    WEAK_INTERACTIONS = [
        Menu.NCI_ANALYSIS,
        Menu.NCI_PROMOLECULAR,
        Menu.ANCI_ANALYSIS,
        Menu.IRI_ANALYSIS,
        Menu.DORI_ANALYSIS,
        Menu.VDW_POTENTIAL,
        Menu.IGM_ANALYSIS,
        Menu.IGMH_ANALYSIS,
        Menu.AIGM_ANALYSIS,
        Menu.MIGM_ANALYSIS,
        Menu.AMIGM_ANALYSIS,
    ]

    SPECTRA = [
        Menu.PLOT_IR_SPECTRUM,
        Menu.PLOT_RAMAN_SPECTRUM,
        Menu.PLOT_UV_VIS_SPECTRUM,
        Menu.PLOT_ECD_SPECTRUM,
        Menu.PLOT_VCD_SPECTRUM,
        Menu.PLOT_ROA_SPECTRUM,
        Menu.PLOT_NMR_SPECTRUM,
        Menu.PLOT_FLUORESCENCE_SPECTRUM,
        Menu.PLOT_PVS,
        Menu.PREDICT_COLOR,
    ]

    SURFACES = [
        Menu.SURFACE_ANALYSIS_ESP,
        Menu.SURFACE_ANALYSIS_ALIE,
        Menu.SURFACE_AREA_VOLUME,
        Menu.BECKE_SURFACE,
        Menu.HIRSHFELD_SURFACE,
        Menu.SURFACE_EXTREMA,
        Menu.HIRSHFELD_SURFACE_FINGERPRINT,
    ]

    AROMATICITY = [
        Menu.AICD_ANALYSIS,
        Menu.NICS_POINT,
        Menu.ICSS_ANALYSIS,
        Menu.NICS_SCAN,
        Menu.BIRD_INDEX,
        Menu.HOMA_INDEX,
        Menu.HOMAC_HOMER,
        Menu.STANGER_INDEX,
        Menu.NICS_1D_SCAN,
        Menu.NICS_2D_MAP,
    ]

    CDFT = [
        Menu.CDFT_ANALYSIS,
        Menu.FUKUI_FUNCTION,
        Menu.DUAL_DESCRIPTOR,
        Menu.CONDENSED_FUKUI,
        Menu.LOCAL_HARDNESS,
        Menu.LOCAL_IONIZATION_ENERGY,
    ]

    DOS = [
        Menu.PLOT_DOS,
        Menu.PLOT_PDOS,
        Menu.PLOT_OPDOS,
        Menu.PLOT_LDOS,
        Menu.PLOT_PHOTOELECTRON_SPECTRUM,
        Menu.PLOT_COHP,
    ]

    BASIN = [
        Menu.BASIN_ANALYSIS_AIM,
        Menu.BASIN_ANALYSIS_ELF,
        Menu.BASIN_INTEGRATE_PROPERTY,
        Menu.BASIN_ANALYSIS_ESP,
        Menu.BASIN_ANALYSIS_LOL,
        Menu.BASIN_ANALYSIS_LOL_ALPHA,
        Menu.BASIN_ANALYSIS_CUSTOM,
    ]

    EXCITATION = [
        Menu.HOLE_ELECTRON_ANALYSIS,
        Menu.TRANSITION_DENSITY_MATRIX,
        Menu.CHARGE_TRANSFER_ANALYSIS,
        Menu.DELTA_R_INDEX,
        Menu.TRANSITION_DIPOLE_MOMENTS,
        Menu.GENERATE_NTO,
        Menu.IFCT_ANALYSIS,
        Menu.LAMBDA_INDEX,
        Menu.CTS_ANALYSIS,
        Menu.CONDITIONAL_DENSITY,
    ]

    CUBES = [
        Menu.CUBE_DENSITY,
        Menu.CUBE_SPIN_DENSITY,
        Menu.CUBE_ELF,
        Menu.CUBE_LOL,
        Menu.CUBE_ESP,
        Menu.CUBE_LAPLACIAN,
        Menu.CUBE_GRADIENT_NORM,
        Menu.CUBE_KINETIC_G,
        Menu.CUBE_KINETIC_K,
        Menu.CUBE_ALIE,
        Menu.CUBE_RDG,
        Menu.CUBE_SIGN_LAMBDA2_RHO,
        Menu.CUBE_SOURCE_FUNCTION,
        Menu.CUBE_ORBITAL_WAVEFUNCTION,
        Menu.CUBE_FUKUI_MINUS,
        Menu.CUBE_FUKUI_PLUS,
        Menu.CUBE_DUAL_DESCRIPTOR,
        Menu.CUBE_PROMOLECULAR_DENSITY,
        Menu.CUBE_DEFORMATION_DENSITY,
        Menu.CUBE_DENSITY_HIGH,
        Menu.CUBE_ESP_HIGH,
        Menu.CUBE_ELF_HIGH,
    ]

    ORBITAL_COMPOSITION = [
        Menu.ORBITAL_COMPOSITION_MULLIKEN,
        Menu.ORBITAL_COMPOSITION_SCPA,
        Menu.ORBITAL_COMPOSITION_STOUT_POLITZER,
        Menu.ORBITAL_COMPOSITION_FRAGMENT_MULLIKEN,
        Menu.ORBITAL_COMPOSITION_FRAGMENT_STOUT,
        Menu.ORBITAL_COMPOSITION_FRAGMENT_SCPA,
        Menu.ORBITAL_COMPOSITION_NAO,
        Menu.ORBITAL_COMPOSITION_HIRSHFELD,
        Menu.ORBITAL_COMPOSITION_BECKE,
        Menu.LOBA_OXIDATION_STATE,
    ]

    ORBITAL_LOCALIZATION = [
        Menu.BOYS_LOCALIZATION,
        Menu.PIPEK_MEZEY_LOCALIZATION,
    ]

    FUZZY_SPACE = [
        Menu.FUZZY_INTEGRATE_PROPERTY,
        Menu.ATOMIC_DIPOLE_MOMENTS,
        Menu.ATOMIC_OVERLAP_MATRIX,
        Menu.LOCALIZATION_DELOCALIZATION_INDEX,
        Menu.PDI_AROMATICITY,
        Menu.FLU_AROMATICITY,
        Menu.FLU_PI_AROMATICITY,
        Menu.FUZZY_INTEGRATE_OVERLAP,
        Menu.CONDENSED_LINEAR_RESPONSE,
        Menu.PARA_LINEAR_RESPONSE,
        Menu.MULTICENTER_DI,
        Menu.ITA_AROMATICITY,
        Menu.ATOMIC_VOLUME_POLARIZABILITY,
        Menu.IFDI_ANALYSIS,
    ]

    EDA = [
        Menu.EDA_FF,
        Menu.EDA_SBL,
        Menu.SOBEDA_ANALYSIS,
        Menu.DISPERSION_ATOMIC_CONTRIBUTION,
    ]

    POLARIZABILITY = [
        Menu.PARSE_POLARIZABILITY,
        Menu.SOS_POLARIZABILITY,
        Menu.POLARIZABILITY_DENSITY,
        Menu.UNIT_SPHERE_POLARIZABILITY,
    ]

    WAVEFUNCTION = [
        Menu.CHECK_WAVEFUNCTION,
        Menu.PRINT_ORBITAL_INFO,
        Menu.PRINT_GTF_INFO,
        Menu.PRINT_BASIS_INFO,
        Menu.PRINT_DENSITY_MATRIX,
        Menu.PRINT_OVERLAP_MATRIX,
        Menu.MODIFY_OCCUPATION,
        Menu.SAVE_WAVEFUNCTION_WFN,
        Menu.DELETE_INNER_ORBITALS,
    ]

    LINE_PLOTS = [
        Menu.LINE_ESP,
        Menu.LINE_ELECTRON_DENSITY,
        Menu.LINE_LAPLACIAN,
        Menu.LINE_ELF,
        Menu.LINE_LOL,
        Menu.LINE_RDG,
        Menu.LINE_SPIN_DENSITY,
        Menu.LINE_GRADIENT_NORM,
        Menu.LINE_KINETIC_G,
        Menu.LINE_KINETIC_K,
        Menu.LINE_ALIE,
        Menu.LINE_SOURCE_FUNCTION,
    ]

    PLANE_MAPS = [
        Menu.PLANE_MAP_DENSITY,
        Menu.PLANE_MAP_ESP,
        Menu.PLANE_MAP_ELF,
        Menu.PLANE_MAP_LOL,
        Menu.PLANE_MAP_GRADIENT,
        Menu.PLANE_MAP_LAPLACIAN,
        Menu.PLANE_MAP_SPIN_DENSITY,
        Menu.PLANE_MAP_RDG,
        Menu.PLANE_MAP_SIGN_LAMBDA2_RHO,
        Menu.PLANE_MAP_ALIE,
        Menu.PLANE_MAP_KINETIC_G,
        Menu.PLANE_MAP_KINETIC_K,
        Menu.PLANE_MAP_SOURCE_FUNCTION,
        Menu.PLANE_MAP_ORBITAL_WAVEFUNCTION,
        Menu.PLANE_MAP_FUKUI_MINUS,
        Menu.PLANE_MAP_FUKUI_PLUS,
        Menu.PLANE_MAP_DUAL_DESCRIPTOR,
        Menu.PLANE_MAP_DEFORMATION_DENSITY,
        Menu.PLANE_MAP_PROMOLECULAR_DENSITY,
    ]

    CATEGORIES: dict[str, list[Menu]] = {
        "charges": CHARGES,
        "bond_orders": BOND_ORDERS,
        "topology": TOPOLOGY,
        "weak_interactions": WEAK_INTERACTIONS,
        "spectra": SPECTRA,
        "surfaces": SURFACES,
        "aromaticity": AROMATICITY,
        "cdft": CDFT,
        "dos": DOS,
        "basin": BASIN,
        "excitation": EXCITATION,
        "cubes": CUBES,
        "orbital_composition": ORBITAL_COMPOSITION,
        "orbital_localization": ORBITAL_LOCALIZATION,
        "fuzzy_space": FUZZY_SPACE,
        "eda": EDA,
        "polarizability": POLARIZABILITY,
        "wavefunction": WAVEFUNCTION,
        "line_plots": LINE_PLOTS,
        "plane_maps": PLANE_MAPS,
    }

    def __init__(
        self,
        filepath: str | Path,
        config: MultiwfnConfig | None = None,
    ) -> None:
        """Initialize Analysis with a wavefunction file."""
        self._filepath = Path(filepath)
        if not self._filepath.exists():
            raise FileNotFoundError(f"File not found: {filepath}")
        self._config = config or MultiwfnConfig()

    @property
    def filepath(self) -> Path:
        """Get the wavefunction file path."""
        return self._filepath

    @property
    def config(self) -> MultiwfnConfig:
        """Get the configuration."""
        return self._config

    @config.setter
    def config(self, value: MultiwfnConfig) -> None:
        """Set the configuration."""
        self._config = value

    # =========================================================================
    # Core run methods
    # =========================================================================

    def run(
        self,
        menu: Menu,
        verbose: bool = False,
    ) -> MultiwfnResult:
        """Run a single analysis.

        Parameters
        ----------
        menu : Menu
            Menu item to run
        verbose : bool
            Print output during execution

        Returns
        -------
        MultiwfnResult
            Execution result
        """
        job = MultiwfnJob(self._filepath, config=self._config)
        job.add_menu(menu)
        return job.run(verbose=verbose)

    def run_multiple(
        self,
        menus: list[Menu],
        verbose: bool = False,
    ) -> MultiwfnResult:
        """Run multiple analyses in a single Multiwfn session.

        Parameters
        ----------
        menus : list of Menu
            Menu items to run
        verbose : bool
            Print output

        Returns
        -------
        MultiwfnResult
            Combined result from all analyses
        """
        job = MultiwfnJob(self._filepath, config=self._config)
        for menu in menus:
            job.add_menu(menu)
        return job.run(verbose=verbose)

    def run_batch(
        self,
        menus: list[Menu],
        verbose: bool = False,
    ) -> dict[str, MultiwfnResult]:
        """Run analyses separately, returning individual results.

        Parameters
        ----------
        menus : list of Menu
            Menu items to run
        verbose : bool
            Print output

        Returns
        -------
        dict[str, MultiwfnResult]
            Map of menu name to result
        """
        return {menu.name.lower(): self.run(menu, verbose) for menu in menus}

    # =========================================================================
    # Category runners
    # =========================================================================

    def run_charges(self, verbose: bool = False) -> dict[str, MultiwfnResult]:
        """Run all charge analyses."""
        return self.run_batch(self.CHARGES, verbose)

    def run_bond_orders(
        self, verbose: bool = False
    ) -> dict[str, MultiwfnResult]:
        """Run all bond order analyses."""
        return self.run_batch(self.BOND_ORDERS, verbose)

    def run_topology(self, verbose: bool = False) -> dict[str, MultiwfnResult]:
        """Run all topology analyses."""
        return self.run_batch(self.TOPOLOGY, verbose)

    def run_weak_interactions(
        self, verbose: bool = False
    ) -> dict[str, MultiwfnResult]:
        """Run all weak interaction analyses."""
        return self.run_batch(self.WEAK_INTERACTIONS, verbose)

    def run_spectra(self, verbose: bool = False) -> dict[str, MultiwfnResult]:
        """Run all spectrum analyses."""
        return self.run_batch(self.SPECTRA, verbose)

    def run_surfaces(self, verbose: bool = False) -> dict[str, MultiwfnResult]:
        """Run all surface analyses."""
        return self.run_batch(self.SURFACES, verbose)

    def run_aromaticity(
        self, verbose: bool = False
    ) -> dict[str, MultiwfnResult]:
        """Run all aromaticity analyses."""
        return self.run_batch(self.AROMATICITY, verbose)

    def run_cdft(self, verbose: bool = False) -> dict[str, MultiwfnResult]:
        """Run all CDFT analyses."""
        return self.run_batch(self.CDFT, verbose)

    def run_dos(self, verbose: bool = False) -> dict[str, MultiwfnResult]:
        """Run all density of states analyses."""
        return self.run_batch(self.DOS, verbose)

    def run_basin(self, verbose: bool = False) -> dict[str, MultiwfnResult]:
        """Run all basin analyses."""
        return self.run_batch(self.BASIN, verbose)

    def run_excitation(
        self, verbose: bool = False
    ) -> dict[str, MultiwfnResult]:
        """Run all electron excitation analyses."""
        return self.run_batch(self.EXCITATION, verbose)

    def run_cubes(self, verbose: bool = False) -> dict[str, MultiwfnResult]:
        """Run all cube generation analyses."""
        return self.run_batch(self.CUBES, verbose)

    def run_orbital_composition(
        self, verbose: bool = False
    ) -> dict[str, MultiwfnResult]:
        """Run all orbital composition analyses."""
        return self.run_batch(self.ORBITAL_COMPOSITION, verbose)

    def run_orbital_localization(
        self, verbose: bool = False
    ) -> dict[str, MultiwfnResult]:
        """Run all orbital localization analyses."""
        return self.run_batch(self.ORBITAL_LOCALIZATION, verbose)

    def run_fuzzy_space(
        self, verbose: bool = False
    ) -> dict[str, MultiwfnResult]:
        """Run all fuzzy atomic space analyses."""
        return self.run_batch(self.FUZZY_SPACE, verbose)

    def run_eda(self, verbose: bool = False) -> dict[str, MultiwfnResult]:
        """Run all energy decomposition analyses."""
        return self.run_batch(self.EDA, verbose)

    def run_polarizability(
        self, verbose: bool = False
    ) -> dict[str, MultiwfnResult]:
        """Run all polarizability analyses."""
        return self.run_batch(self.POLARIZABILITY, verbose)

    def run_wavefunction(
        self, verbose: bool = False
    ) -> dict[str, MultiwfnResult]:
        """Run all wavefunction check/modify analyses."""
        return self.run_batch(self.WAVEFUNCTION, verbose)

    def run_line_plots(
        self, verbose: bool = False
    ) -> dict[str, MultiwfnResult]:
        """Run all line property plot analyses."""
        return self.run_batch(self.LINE_PLOTS, verbose)

    def run_plane_maps(
        self, verbose: bool = False
    ) -> dict[str, MultiwfnResult]:
        """Run all plane property map analyses."""
        return self.run_batch(self.PLANE_MAPS, verbose)

    def run_category(
        self,
        category: str,
        verbose: bool = False,
    ) -> dict[str, MultiwfnResult]:
        """Run all analyses in a named category.

        Parameters
        ----------
        category : str
            Category name (charges, bond_orders, topology, etc.)
        verbose : bool
            Print output

        Returns
        -------
        dict[str, MultiwfnResult]
            Results keyed by analysis name
        """
        if category not in self.CATEGORIES:
            available = ", ".join(self.CATEGORIES.keys())
            raise ValueError(
                f"Unknown category: {category}. Available: {available}"
            )
        return self.run_batch(self.CATEGORIES[category], verbose)

    def run_all(
        self,
        verbose: bool = False,
        categories: list[str] | None = None,
    ) -> dict[str, dict[str, MultiwfnResult]]:
        """Run all available analyses.

        Parameters
        ----------
        verbose : bool
            Print output
        categories : list of str, optional
            Specific categories to run (default: all)

        Returns
        -------
        dict[str, dict[str, MultiwfnResult]]
            Nested dict: category -> analysis_name -> result
        """
        cats = categories or list(self.CATEGORIES.keys())
        results: dict[str, dict[str, MultiwfnResult]] = {
            cat: self.run_batch(self.CATEGORIES[cat], verbose)
            for cat in cats
            if cat in self.CATEGORIES
        }
        return results

    # =========================================================================
    # Helpers
    # =========================================================================

    @classmethod
    def list_categories(cls) -> list[str]:
        """List available analysis categories."""
        return list(cls.CATEGORIES.keys())

    @classmethod
    def list_analyses(cls, category: str | None = None) -> list[str]:
        """List available analyses, optionally filtered by category."""
        if category:
            if category not in cls.CATEGORIES:
                raise ValueError(f"Unknown category: {category}")
            return [m.name.lower() for m in cls.CATEGORIES[category]]
        return [m.name.lower() for m in Menu]

    def __repr__(self) -> str:
        return f"MultiwfnAnalysis({self._filepath!r})"