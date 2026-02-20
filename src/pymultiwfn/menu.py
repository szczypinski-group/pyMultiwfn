"""Complete Multiwfn menu mapping for batch automation.

Notes
-----
Based on Multiwfn 3.8 (dev) manual.

"""

from enum import Enum


class Menu(Enum):
    """Enumeration of all Multiwfn menu functions.

    Examples
    --------
    >>> Menu.HIRSHFELD_CHARGE.get_sequence()
    ['7', '1', 'q']

    >>> job.add_menu(Menu.HIRSHFELD_CHARGE)
    >>> job.add_menu(Menu.MAYER_BOND_ORDER)
    """

    # Main Menu 0: Visualization
    VIEW_STRUCTURE = ("0",)

    # Main Menu 2: Topology analysis
    TOPOLOGY_SEARCH_CPS = ("2", "2", "-1")
    TOPOLOGY_GENERATE_PATHS = ("2", "3")
    TOPOLOGY_INTERBASIN_SURFACES = ("2", "4")
    TOPOLOGY_ANALYSIS_COMPLETE = ("2", "2", "-1", "3", "4")
    TOPOLOGY_ESP_ANALYSIS = ("2", "-2", "2", "-1")
    TOPOLOGY_LOL_ANALYSIS = ("2", "-10", "2", "-1")

    # Main Menu 3: Line analysis
    LINE_ESP = ("3", "12")

    # Main Menu 4: Plane maps
    PLANE_MAP_DENSITY = ("4", "1", "1")
    PLANE_MAP_ESP = ("4", "12", "1")
    PLANE_MAP_ELF = ("4", "4", "1")
    PLANE_MAP_LOL = ("4", "5", "1")
    PLANE_MAP_GRADIENT = ("4", "2", "1")
    PLANE_MAP_LAPLACIAN = ("4", "3", "1")

    # Main Menu 5: Cube generation
    CUBE_DENSITY = ("5", "1", "2", "0")
    CUBE_SPIN_DENSITY = ("5", "2", "2", "0")
    CUBE_ELF = ("5", "4", "2", "0")
    CUBE_LOL = ("5", "5", "2", "0")
    CUBE_ESP = ("5", "12", "2", "0")
    CUBE_LAPLACIAN = ("5", "3", "2", "0")
    CUBE_FUKUI_MINUS = ("5", "0", "-1", "2", "0")
    CUBE_FUKUI_PLUS = ("5", "0", "1", "2", "0")
    CUBE_DUAL_DESCRIPTOR = ("5", "0", "3", "2", "0")

    # Main Menu 6: Wavefunction
    # Menu 6 shows wavefunction info immediately, then shows submenu
    # "0" returns to main menu
    CHECK_WAVEFUNCTION = ("6", "0")
    PRINT_ORBITAL_INFO = ("6", "1", "0", "0")  # 1=show orbitals, first 0=back, second 0=main
    PRINT_BASIS_INFO = ("6", "2", "0", "0")
    MODIFY_OCCUPATION = ("6", "3", "0", "0")

    # Main Menu 7: Population analysis and atomic charges
    # Note: Many charge methods require additional inputs:
    # - "1" to select built-in atomic densities
    # - "n" to skip outputting .chg file
    # - "0" to return to main menu
    HIRSHFELD_CHARGE = ("7", "1", "1", "n", "0")
    VDD_POPULATION = ("7", "2", "1", "n", "0")
    MULLIKEN_POPULATION = ("7", "5", "0")
    LOWDIN_POPULATION = ("7", "6", "0")
    BECKE_CHARGE = ("7", "10", "1", "n", "0")
    ADCH_CHARGE = ("7", "11", "1", "n", "0")
    CHELPG_CHARGE = ("7", "12", "n", "0")
    MK_CHARGE = ("7", "13", "n", "0")
    AIM_CHARGE = ("7", "14", "0")
    HIRSHFELD_I_CHARGE = ("7", "15", "1", "n", "0")
    CM5_CHARGE = ("7", "16", "1", "n", "0")
    EEM_CHARGE = ("7", "17", "n", "0")
    RESP_CHARGE = ("7", "18", "n", "0")
    GASTEIGER_CHARGE = ("7", "19", "n", "0")
    MBIS_CHARGE = ("7", "20", "1", "n", "0")

    # Main Menu 8: Orbital composition analysis
    LOBA_OXIDATION_STATE = ("8", "100", "0")

    # Main Menu 9: Bond order analysis
    # Bond order analyses typically need "0" to return to main menu
    MAYER_BOND_ORDER = ("9", "1", "0")
    MULTICENTER_BOND_ORDER = ("9", "2", "0")
    WIBERG_BOND_ORDER = ("9", "3", "0")
    MULLIKEN_BOND_ORDER = ("9", "4", "0")
    FUZZY_BOND_ORDER = ("9", "7", "0")
    LAPLACIAN_BOND_ORDER = ("9", "8", "0")
    WIBERG_DECOMPOSITION = ("9", "9", "0")
    IBSI_ANALYSIS = ("9", "10", "0")
    AV1245_INDEX = ("9", "11", "0")

    # Main Menu 10: DOS
    PLOT_DOS = ("10", "1")
    PLOT_PDOS = ("10", "2")
    PLOT_OPDOS = ("10", "3")
    PLOT_LDOS = ("10", "4")
    PLOT_PHOTOELECTRON_SPECTRUM = ("10", "5")
    PLOT_COHP = ("10", "6")

    # Main Menu 11: Spectra
    PLOT_IR_SPECTRUM = ("11", "1")
    PLOT_RAMAN_SPECTRUM = ("11", "2")
    PLOT_UV_VIS_SPECTRUM = ("11", "3")
    PLOT_ECD_SPECTRUM = ("11", "4")
    PLOT_VCD_SPECTRUM = ("11", "5")
    PLOT_ROA_SPECTRUM = ("11", "6")
    PLOT_NMR_SPECTRUM = ("11", "7")
    PLOT_FLUORESCENCE_SPECTRUM = ("11", "8")
    PLOT_PVS = ("11", "9")
    PREDICT_COLOR = ("11", "10")

    # Main Menu 12: Surface analysis
    SURFACE_ANALYSIS_ESP = ("12", "1")
    SURFACE_ANALYSIS_ALIE = ("12", "2")
    SURFACE_AREA_VOLUME = ("12", "3")
    BECKE_SURFACE = ("12", "4")
    HIRSHFELD_SURFACE = ("12", "5")
    SURFACE_EXTREMA = ("12", "6")

    # Main Menu 13: Grid processing
    EXPORT_CUBE = ("13", "0")
    GRID_MATH_OPERATIONS = ("13", "11")
    GRID_EXTRACT_PLANE = ("13", "2")
    GRID_PLOT_INTEGRAL_CURVE = ("13", "18")

    # Main Menu 14: AdNDP
    ADNDP_ANALYSIS = ("14",)

    # Main Menu 15: Fuzzy atomic space
    FUZZY_INTEGRATE_PROPERTY = ("15", "1")
    ATOMIC_DIPOLE_MOMENTS = ("15", "2")
    ATOMIC_OVERLAP_MATRIX = ("15", "3")
    LOCALIZATION_DELOCALIZATION_INDEX = ("15", "4")
    PDI_AROMATICITY = ("15", "5")
    FLU_AROMATICITY = ("15", "6")
    FLU_PI_AROMATICITY = ("15", "7")
    MULTICENTER_DI = ("15", "11")
    ITA_AROMATICITY = ("15", "12")
    ATOMIC_VOLUME_POLARIZABILITY = ("15", "13")

    # Main Menu 16: CDA
    CDA_ANALYSIS = ("16",)

    # Main Menu 17: Basin analysis
    BASIN_ANALYSIS_AIM = ("17", "1")
    BASIN_ANALYSIS_ELF = ("17", "2")
    BASIN_INTEGRATE_PROPERTY = ("17", "3")

    # Main Menu 18: Excitation analysis
    HOLE_ELECTRON_ANALYSIS = ("18", "1")
    TRANSITION_DENSITY_MATRIX = ("18", "2")
    CHARGE_TRANSFER_ANALYSIS = ("18", "3")
    DELTA_R_INDEX = ("18", "4")
    TRANSITION_DIPOLE_MOMENTS = ("18", "5")
    GENERATE_NTO = ("18", "6")
    IFCT_ANALYSIS = ("18", "8")
    LAMBDA_INDEX = ("18", "14")
    CTS_ANALYSIS = ("18", "16")

    # Main Menu 19: Localization
    BOYS_LOCALIZATION = ("19", "1")
    PIPEK_MEZEY_LOCALIZATION = ("19", "2")

    # Main Menu 20: Weak interactions
    NCI_ANALYSIS = ("20", "1")
    NCI_PROMOLECULAR = ("20", "2")
    ANCI_ANALYSIS = ("20", "3")
    IRI_ANALYSIS = ("20", "4")
    DORI_ANALYSIS = ("20", "5")
    VDW_POTENTIAL = ("20", "6")
    IGM_ANALYSIS = ("20", "10")
    IGMH_ANALYSIS = ("20", "11")
    AIGM_ANALYSIS = ("20", "12")
    MIGM_ANALYSIS = ("20", "-10")

    # Main Menu 21: EDA
    EDA_FF = ("21", "1")
    EDA_SBL = ("21", "2")
    SOBEDA_ANALYSIS = ("21", "3")
    DISPERSION_ATOMIC_CONTRIBUTION = ("21", "4")

    # Main Menu 22: CDFT
    CDFT_ANALYSIS = ("22",)
    FUKUI_FUNCTION = ("22", "1")
    DUAL_DESCRIPTOR = ("22", "2")
    CONDENSED_FUKUI = ("22", "3")

    # Main Menu 23: ETS-NOCV
    ETS_NOCV_ANALYSIS = ("23",)

    # Main Menu 24: Polarizability
    PARSE_POLARIZABILITY = ("24", "1")
    SOS_POLARIZABILITY = ("24", "2")
    POLARIZABILITY_DENSITY = ("24", "3")
    UNIT_SPHERE_POLARIZABILITY = ("24", "5")

    # Main Menu 25: Aromaticity
    ICSS_ANALYSIS = ("25", "3")
    NICS_SCAN = ("25", "4")
    HOMA_INDEX = ("25", "6")
    HOMAC_HOMER = ("25", "7")
    NICS_1D_SCAN = ("25", "13")
    NICS_2D_MAP = ("25", "14")

    # Main Menu 100: Utilities Part 1
    SCATTER_GRAPH_TWO_FUNCTIONS = ("100", "1")
    EXPORT_VARIOUS_FILES = ("100", "2")
    VDW_VOLUME = ("100", "3")
    INTEGRATE_WHOLE_SPACE = ("100", "4")
    ORBITAL_OVERLAP_INTEGRAL = ("100", "5")
    MONITOR_SCF_CONVERGENCE = ("100", "6")
    FRAGMENT_GUESS_INPUT = ("100", "8")
    ATOMIC_COORDINATION = ("100", "9")
    ORBITAL_OVERLAP_CENTROID = ("100", "11")
    BIORTHOGONALIZATION = ("100", "12")
    LOLIPOP_INDEX = ("100", "14")
    INTERMOLECULAR_OVERLAP = ("100", "15")
    GENERATE_FOCK_MATRIX = ("100", "17")
    ELECTRON_TRANSPORT_ROUTE = ("100", "18")
    COMBINE_FRAGMENTS = ("100", "19")
    HELLMANN_FEYNMAN_FORCES = ("100", "20")
    GEOMETRY_PROPERTIES = ("100", "21")
    DETECT_PI_ORBITALS = ("100", "22")
    FIT_FUNCTION_TO_ATOMS = ("100", "23")

    # Main Menu 200: Utilities Part 2
    CVB_INDEX = ("200", "1")
    ATOMIC_BOND_DIPOLES = ("200", "2")
    MULTIPLE_ORBITAL_CUBES = ("200", "3")
    RADIAL_DISTRIBUTION = ("200", "5")
    ORBITAL_CORRESPONDENCE = ("200", "6")
    AVERAGE_BOND_LENGTH = ("200", "9")
    ORBITAL_INTEGRALS = ("200", "10")
    FUNCTION_MOMENTS = ("200", "11")
    ENERGY_INDEX = ("200", "12")
    ORBITAL_CONTRIBUTIONS_TO_GRID = ("200", "13")
    DOMAIN_ANALYSIS = ("200", "14")
    CORRELATION_INDEX = ("200", "15")
    NATURAL_ORBITALS = ("200", "16")
    COULOMB_EXCHANGE_INTEGRALS = ("200", "17")
    BLA_BOA_ANALYSIS = ("200", "18")
    SPATIAL_DELOCALIZATION_INDEX = ("200", "19")
    BOD_NADO_ANALYSIS = ("200", "20")
    LOWDIN_ORTHOGONALIZATION = ("200", "21")

    # Main Menu 300: Utilities Part 3
    FREE_VOLUME_IN_CELL = ("300", "1")
    FIT_ATOMIC_RADIAL_DENSITY = ("300", "2")
    STM_IMAGE = ("300", "4")
    ELECTRIC_MULTIPOLE_MOMENTS = ("300", "5")
    ORBITAL_ENERGIES_FROM_FOCK = ("300", "6")
    GEOMETRY_OPERATIONS = ("300", "7")
    SURFACE_DISTANCE_PROJECTION = ("300", "8")
    DETERMINE_FERMI_LEVEL = ("300", "9")

    def get_sequence(self) -> list[str]:
        """Get the menu sequence with 'q' appended.

        Parameters
        ----------
        *args
            Additional arguments to append (e.g., orbital index)
        """
        result = list(self.value)
        # if args:
        #     result.extend(str(arg) for arg in args)
        result.append("q")
        return result

    @classmethod
    def search(cls, keyword: str) -> list["Menu"]:
        """Search menu items by keyword in name."""
        keyword = keyword.lower()
        return [m for m in cls if keyword in m.name.lower()]

    @classmethod
    def list_all(cls) -> list[str]:
        """List all menu item names."""
        return [m.name for m in cls]


def custom_sequence(seq: list[str]) -> list[str]:
    """Create a custom menu sequence.

    Parameters
    ----------
    seq : list of str
        Menu choices

    Returns
    -------
    list of str
        Sequence with 'q' appended if not present
    """
    if not isinstance(seq, list):
        raise TypeError("seq must be a list of menu choices")
    seq = list(seq)
    if not seq or seq[-1] != "q":
        seq.append("q")
    return seq