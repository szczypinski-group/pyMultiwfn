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

    # Category definitions
    CHARGES = [
        Menu.HIRSHFELD_CHARGE,
        Menu.VDD_POPULATION,
        Menu.MULLIKEN_POPULATION,
        Menu.LOWDIN_POPULATION,
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
    ]

    BOND_ORDERS = [
        Menu.MAYER_BOND_ORDER,
        Menu.MULTICENTER_BOND_ORDER,
        Menu.WIBERG_BOND_ORDER,
        Menu.MULLIKEN_BOND_ORDER,
        Menu.FUZZY_BOND_ORDER,
        Menu.LAPLACIAN_BOND_ORDER,
    ]

    TOPOLOGY = [
        Menu.TOPOLOGY_SEARCH_CPS,
        Menu.TOPOLOGY_GENERATE_PATHS,
        Menu.TOPOLOGY_ANALYSIS_COMPLETE,
    ]

    WEAK_INTERACTIONS = [
        Menu.NCI_ANALYSIS,
        Menu.IRI_ANALYSIS,
        Menu.DORI_ANALYSIS,
        Menu.IGM_ANALYSIS,
        Menu.IGMH_ANALYSIS,
    ]

    SPECTRA = [
        Menu.PLOT_IR_SPECTRUM,
        Menu.PLOT_RAMAN_SPECTRUM,
        Menu.PLOT_UV_VIS_SPECTRUM,
        Menu.PLOT_ECD_SPECTRUM,
        Menu.PLOT_VCD_SPECTRUM,
        Menu.PLOT_ROA_SPECTRUM,
        Menu.PLOT_NMR_SPECTRUM,
    ]

    SURFACES = [
        Menu.SURFACE_ANALYSIS_ESP,
        Menu.SURFACE_ANALYSIS_ALIE,
        Menu.SURFACE_AREA_VOLUME,
    ]

    AROMATICITY = [
        Menu.ICSS_ANALYSIS,
        Menu.NICS_SCAN,
        Menu.HOMA_INDEX,
    ]

    CDFT = [
        Menu.FUKUI_FUNCTION,
        Menu.DUAL_DESCRIPTOR,
        Menu.CONDENSED_FUKUI,
    ]

    DOS = [
        Menu.PLOT_DOS,
        Menu.PLOT_PDOS,
    ]

    BASIN = [
        Menu.BASIN_ANALYSIS_AIM,
        Menu.BASIN_ANALYSIS_ELF,
    ]

    EXCITATION = [
        Menu.HOLE_ELECTRON_ANALYSIS,
        Menu.CHARGE_TRANSFER_ANALYSIS,
    ]

    CUBES = [
        Menu.CUBE_DENSITY,
        Menu.CUBE_SPIN_DENSITY,
        Menu.CUBE_ELF,
        Menu.CUBE_LOL,
        Menu.CUBE_ESP,
        Menu.CUBE_LAPLACIAN,
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

    def run_cubes(self, verbose: bool = False) -> dict[str, MultiwfnResult]:
        """Run all cube generation analyses."""
        return self.run_batch(self.CUBES, verbose)

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
