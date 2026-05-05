"""Convenience groups of Multifwn analyses.

This is a strongly opinionated set of analyses that one could run on their
wavefunction files in order to obtain related results.

"""

from enum import Enum

from pymultiwfn.enums.menu import Menu


class AnalysisClasses(Enum):
    """Enum containing various Multiwfn analysis groups."""

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
        Menu.CM5_CHARGE,
        Menu.EEM_CHARGE,
        Menu.RESP_CHARGE,
        Menu.GASTEIGER_CHARGE,
        Menu.MBIS_CHARGE,
        Menu.MULLIKEN_DECOMPOSE_ATOMIC_POPULATION,
        Menu.MULLIKEN_DECOMPOSE_BASIS_FUNCTION,
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
        Menu.TOPOLOGY_VISUALISE_CPS,
        Menu.TOPOLOGY_CP_NUCLEAR_POSITION,
        Menu.TOPOLOGY_CP_MIDPOINTS,
        Menu.TOPOLOGY_CP_TRIANGLE_CENTRES,
        Menu.TOPOLOGY_CP_PYRAMID_CENTRES,
        Menu.TOPOLOGY_CP_SPHERE_POINTS,
        Menu.TOPOLOGY_CP_REAL_SPACE_POINTS,
        Menu.TOPOLOGY_CP_PATHS_3MINUS3_3MINUS1,
        Menu.TOPOLOGY_CP_PATHS_3PLUS1_3PLUS3,
        Menu.TOPOLOGY_SEARCH_CPS,
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
        Menu.QMSA_ESP,
        Menu.QMSA_ALIE,
        Menu.QMSA_LEA,
        Menu.QMSA_LEAE,
        Menu.QMSA_EDR,
        Menu.QMSA_MAXEDR,
        Menu.QMSA_EDENSITY,
        Menu.QMSA_LAMBDA2_RHO,
        Menu.SURFACE_DISTANCE_PROJECTION,
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
        Menu.HOMA_BIRD_AROMATICITY,
        Menu.NICS_ZZ_NONPLANAR,
        Menu.LOLIPOP_INDEX,
    ]

    CDFT = [
        Menu.CDFT_ANALYSIS,
        Menu.CONDENSED_FUKUI,
        Menu.ORBITAL_WEIGHTED_FUKUI,
        Menu.GRID_FUKUI,
        Menu.LOCAL_HARDNESS,
        Menu.ORBITAL_WEIGHTS,
        Menu.SUPERDELOCALIZABILITIES_NUC_E,
    ]

    DOS = [
        Menu.PLOT_TDOS,
        Menu.PLOT_TDOS_OPDOS,
    ]

    BASIN = [
        Menu.BASIN_ANALYSIS_AIM,
        Menu.BASIN_ANALYSIS_ELF,
        Menu.BASIN_INTEGRATE_PROPERTY,
        Menu.BASIN_ANALYSIS_ESP,
        Menu.BASIN_ANALYSIS_LOL,
        Menu.BASIN_ANALYSIS_LOL_ALPHA,
        Menu.BASIN_ANALYSIS_CUSTOM,
        Menu.BASIN_ANALYSIS_RHO,
        Menu.BASIN_EDENSITY,
        Menu.BASIN_NORM_RHO,
        Menu.BASIN_LAPLACIAN,
        Menu.BASIN_ORB_WFN_HOMO,
        Menu.BASIN_ORB_WFN_LUMO,
        Menu.BASIN_ESPIN_DENSITY,
        Menu.BASIN_KR,
        Menu.BASIN_GR,
        Menu.BASIN_ESP_CHARGES,
        Menu.BASIN_ELF,
        Menu.BASIN_LOL,
        Menu.BASIN_LOCAL_ENTROPY,
        Menu.BASIN_ESP,
        Menu.BASIN_RDG,
        Menu.BASIN_RDG_PROMOLECULAR,
        Menu.BASIN_LAMBDA2RHO,
        Menu.BASIN_LAMBDA2RGO_PROMOLECULAR,
        Menu.BASIN_ALIE,
        Menu.BASIN_EDR,
        Menu.BASIN_ORB_OVERLAP_DR,
        Menu.BASIN_DELTAG_PROMOLECULAR,
        Menu.BASIN_DELTAG_HIRSHFELD,
        Menu.BASIN_IRI,
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
        Menu.CUBE_LOCAL_INFORMATION_ENTROPY,
        Menu.CUBE_ELECTROSTATIC_POTENTIAL_FROM_CHARGE,
        Menu.CUBE_RDG_QUICK,
        Menu.CUBE_CORRELATION_HOLE_ALPHA,
        Menu.CUBE_EDR,
        Menu.CUBE_DELTAG_PROMOLECULAR_APPROX,
        Menu.CUBE_DELTAG_HIRSHFELD_PAER_APPROX,
        Menu.CUBE_IRI,
        Menu.CUBE_VDW_POTENTIAL,
    ]

    ORBITAL_COMPOSITION = [
        Menu.ORBITAL_COMPOSITION_MULLIKEN_ALL,
        Menu.ORBITAL_COMPOSITION_SCPA_ALL,
        Menu.ORBITAL_COMPOSITION_STOUT_POLITZER_ALL,
        Menu.ORBITAL_COMPOSITION_FRAGMENT_MULLIKEN,
        Menu.ORBITAL_COMPOSITION_FRAGMENT_STOUT,
        Menu.ORBITAL_COMPOSITION_FRAGMENT_SCPA,
        Menu.ORBITAL_COMPOSITION_NAO,
        Menu.ORBITAL_COMPOSITION_HIRSHFELD,
        Menu.ORBITAL_COMPOSITION_BECKE,
        Menu.LOBA_OXIDATION_STATE,
        Menu.ORBITAL_COMPOSITION_MULLIKEN_HOMO,
        Menu.ORBITAL_COMPOSITION_MULLIKEN_LUMO,
        Menu.ORBITAL_COMPOSITION_STOUT_POLITZER_HOMO,
        Menu.ORBITAL_COMPOSITION_STOUT_POLITZER_LUMO,
        Menu.ORBITAL_COMPOSITION_SCPA_HOMO,
        Menu.ORBITAL_COMPOSITION_SCPA_LUMO,
        Menu.FRAGMENT_CONTRIBUTION_HIRSHFELD_HOMO,
        Menu.FRAGMENT_CONTRIBUTION_HIRSHFELD_LUMO,
        Menu.FRAGMENT_CONTRIBUTION_HIRSHFELD_ALL,
        Menu.ATOM_CONTRIBUTION_HIRSHFELD,
        Menu.FRAGMENT_CONTRIBUTION_BECKE_HOMO,
        Menu.FRAGMENT_CONTRIBUTION_BECKE_LUMO,
        Menu.FRAGMENT_CONTRIBUTION_BECKE_ALL,
        Menu.ATOM_CONTRIBUTION_BECKE,
    ]

    ORBITAL_LOCALIZATION = [
        Menu.BOYS_LOCALIZATION_OCCUPIED,
        Menu.BOYS_LOCALIZATION_ALL,
        Menu.PIPEK_MEZEY_LOCALIZATION_HIRSHFELD_OCCUPIED,
        Menu.PIPEK_MEZEY_LOCALIZATION_HIRSHFELD_ALL,
        Menu.PIPEK_MEZEY_LOCALIZATION_LOWDIN_OCUPIED,
        Menu.PIPEK_MEZEY_LOCALIZATION_LOWDIN_ALL,
        Menu.PIPEK_MEZEY_LOCALIZATION_BECKE_OCCUPIED,
        Menu.PIPEK_MEZEY_LOCALIZATION_BECKE_ALL,
    ]

    FUZZY_SPACE = [
        Menu.FUZZY_INTEGRATE_EDENSITY,
        Menu.FUZZY_INTEGRATE_NORM_RHO,
        Menu.FUZZY_INTEGRATE_LAPLACIAN,
        Menu.FUZZY_INTEGRATE_ORB_WFN_HOMO,
        Menu.FUZZY_INTEGRATE_ORB_WFN_LUMO,
        Menu.FUZZY_INTEGRATE_ESPIN_DENSITY,
        Menu.FUZZY_INTEGRATE_KR,
        Menu.FUZZY_INTEGRATE_GR,
        Menu.FUZZY_INTEGRATE_ESP_CHARGES,
        Menu.FUZZY_INTEGRATE_ELF,
        Menu.FUZZY_INTEGRATE_LOL,
        Menu.FUZZY_INTEGRATE_LOCAL_ENTROPY,
        Menu.FUZZY_INTEGRATE_ESP,
        Menu.FUZZY_INTEGRATE_RDG,
        Menu.FUZZY_INTEGRATE_RDG_PROMOLECULAR,
        Menu.FUZZY_INTEGRATE_LAMBDA2RHO,
        Menu.FUZZY_INTEGRATE_LAMBDA2RGO_PROMOLECULAR,
        Menu.FUZZY_INTEGRATE_ALIE,
        Menu.FUZZY_INTEGRATE_EDR,
        Menu.FUZZY_INTEGRATE_ORB_OVERLAP_DR,
        Menu.FUZZY_INTEGRATE_DELTAG_PROMOLECULAR,
        Menu.FUZZY_INTEGRATE_DELTAG_HIRSHFELD,
        Menu.FUZZY_INTEGRATE_IRI,
        Menu.FUZZY_MULTIPOLE,
        Menu.FUZZY_OVERLAP_EDENSITY,
        Menu.FUZZY_OVERLAP_LAPLACIAN,
        Menu.FUZZY_OVERLAP_ORB_WFN_HOMO,
        Menu.FUZZY_OVERLAP_ORB_WFN_LUMO,
        Menu.FUZZY_OVERLAP_ESPIN_DENSITY,
        Menu.FUZZY_OVERLAP_KR,
        Menu.FUZZY_OVERLAP_GR,
        Menu.FUZZY_OVERLAP_ESP_CHARGES,
        Menu.FUZZY_OVERLAP_ELF,
        Menu.FUZZY_OVERLAP_LOL,
        Menu.FUZZY_OVERLAP_LOCAL_ENTROPY,
        Menu.FUZZY_OVERLAP_ESP,
        Menu.FUZZY_OVERLAP_RDG,
        Menu.FUZZY_OVERLAP_RDG_PROMOLECULAR,
        Menu.FUZZY_OVERLAP_LAMBDA2RHO,
        Menu.FUZZY_OVERLAP_LAMBDA2RGO_PROMOLECULAR,
        Menu.FUZZY_OVERLAP_ALIE,
        Menu.FUZZY_OVERLAP_EDR,
        Menu.FUZZY_OVERLAP_ORB_OVERLAP_DR,
        Menu.FUZZY_OVERLAP_DELTAG_PROMOLECULAR,
        Menu.FUZZY_OVERLAP_DELTAG_HIRSHFELD,
        Menu.FUZZY_OVERLAP_IRI,
        Menu.CLRK_MATRIX,
        Menu.ATOMIC_OVERLAP_MATRIX,
        Menu.LOCALIZATION_DELOCALIZATION_INDEX,
        Menu.PDI_AROMATICITY,
        Menu.FLU_AROMATICITY,
        Menu.FLU_PI_AROMATICITY,
        Menu.CONDENSED_LINEAR_RESPONSE,
        Menu.PARA_LINEAR_RESPONSE,
        Menu.SPATIAL_DELOCALIZATION_INDEX,
    ]

    EDA = [
        Menu.EDA_FF,
        Menu.EDA_SBL,
        Menu.SOBEDA_ANALYSIS,
        Menu.DISPERSION_ATOMIC_CONTRIBUTION,
        Menu.ETS_NOCV_ANALYSIS,
        Menu.CDA_ANALYSIS,
    ]

    POLARIZABILITY = [
        Menu.PARSE_POLARIZABILITY,
        Menu.SOS_POLARIZABILITY,
        Menu.POLARIZABILITY_DENSITY,
        Menu.UNIT_SPHERE_POLARIZABILITY,
        Menu.PARSE_POLARIZABILITY_GAUSSIAN,
        Menu.SOS_HYPERPOLARIZABILITY,
    ]

    WAVEFUNCTION = [
        Menu.PRINT_DENSITY_MATRIX,
        Menu.PRINT_INTEGRAL_MATRIX_OVERLAP,
        Menu.PRINT_INTEGRAL_MATRIX_ELECTRIC_DIPOLE,
        Menu.PRINT_INTEGRAL_MATRIX_QUADRUPOLE,
        Menu.PRINT_COEFFICIENT_MATRIX,
        Menu.PRINT_INTEGRAL_MATRIX_FOCK,
        Menu.PRINT_INTEGRAL_MATRIX_MAGNETIC_DIPOLE,
        Menu.PRINT_INTEGRAL_MATRIX_VELOCITY,
        Menu.PRINT_INTEGRAL_MATRIX_EKINETIC,
        Menu.PRINT_INTEGRAL_MATRIX_OCTOPOLE,
        Menu.PRINT_INTEGRAL_MATRIX_HEXADECAPOLE,
        Menu.NATURAL_ORBITALS,
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

    GRID_PROCESSING = [
        Menu.EXPORT_CUBE,
        Menu.EXPORT_GRID_ALL_POINTS,
        Menu.GRID_EXTRACT_PLANE_XY,
        Menu.GRID_EXTRACT_PLANE_XZ,
        Menu.GRID_EXTRACT_PLANE_YZ,
        Menu.GRID_AVERAGE_XY,
        Menu.GRID_AVERAGE_XZ,
        Menu.GRID_AVERAGE_YZ,
        Menu.GRID_EXTRACT_PLANE_3ATOMS,
        Menu.GRID_EXTRACT_PLANE_3POINTS,
        Menu.GRID_EXTRACT_VALUE_RANGE,
        Menu.GRID_MATH_OPERATIONS,
        Menu.GRID_MAP_TO_ISOSURFACE,
        Menu.GRID_SET_VALUE_DISTANCE,
        Menu.GRID_SET_VALUE_FRAGMENT,
        Menu.GRID_SET_VALUE_RANGE,
        Menu.GRID_SCALE_RANGE,
        Menu.GRID_STATISTIC_DATA,
        Menu.GRID_PLOT_INTEGRAL_CURVE,
        Menu.GRID_VISUALIZE_ISOSURFACE,
    ]

    UTILITIES = [
        Menu.VIEW_STRUCTURE,
        Menu.PROPERTIES_AT_POINT,
        Menu.GEOMETRY_PROPERTIES,
        Menu.GEOMETRY_OPERATIONS,
        Menu.ELECTRIC_MULTIPOLE_MOMENTS,
        Menu.BLA_BOA_ANALYSIS,
        Menu.ATOMIC_COORDINATION,
        Menu.AVERAGE_BOND_LENGTH,
        Menu.HELLMANN_FEYNMAN_FORCES,
        Menu.VDW_VOLUME,
        Menu.INTEGRATE_WHOLE_SPACE,
        Menu.RING_AREA_PERIMETER,
        Menu.DETECT_PI_ORBITALS,
        Menu.FIT_FUNCTION_TO_ATOMS,
        Menu.FUNCTION_MOMENTS,
        Menu.ENERGY_INDEX,
        Menu.CORRELATION_INDEX,
        Menu.DOMAIN_ANALYSIS,
        Menu.BOD_NADO_ANALYSIS,
        Menu.CVB_INDEX,
        Menu.ATOMIC_BOND_DIPOLES,
        Menu.RADIAL_DISTRIBUTION,
        Menu.DETERMINE_FERMI_LEVEL,
        Menu.FREE_VOLUME_IN_CELL,
        Menu.FIT_ATOMIC_RADIAL_DENSITY,
        Menu.STM_IMAGE,
    ]

    ORBITAL_ANALYSIS = [
        Menu.ORBITAL_OVERLAP_INTEGRAL,
        Menu.ORBITAL_OVERLAP_CENTROID,
        Menu.ORBITAL_CORRESPONDENCE,
        Menu.ORBITAL_INTEGRALS,
        Menu.ORBITAL_INTEGRAL_ELECTRIC_DIPOLE,
        Menu.ORBITAL_INTEGRAL_MAGNETIC_DIPOLE,
        Menu.ORBITAL_INTEGRAL_VELOCITY,
        Menu.ORBITAL_INTEGRAL_KINETIC_ENERGY,
        Menu.ORBITAL_INTEGRAL_OVERLAP,
        Menu.ORBITAL_CONTRIBUTIONS_TO_GRID,
        Menu.ORBITAL_ENERGIES_FROM_FOCK,
        Menu.MULTIPLE_ORBITAL_CUBES,
        Menu.ICSS_CUBES,
        Menu.BIORTHOGONALIZATION,
        Menu.LOWDIN_ORTHOGONALIZATION,
        Menu.COULOMB_EXCHANGE_INTEGRALS,
    ]

    SPATIAL_DELOCALIZATION = [
        Menu.SPACIAL_DELOCALISATION_EDENSITY,
        Menu.SPACIAL_DELOCALISATIOn_NORM_RHO,
        Menu.SPATIAL_DELOCALISATION_LAPLACIAN,
        Menu.SPATIAL_DELOCALISATION_ORB_WFN,
        Menu.SPATIAL_DELOCALISATION_ESPIN_DENSITY,
        Menu.SPATIAL_DELOCALISATION_KR,
        Menu.SPATIAL_DELOCALISATION_GR,
        Menu.SPATIAL_DELOCALISATION_ESP_CHARGES,
        Menu.SPATIAL_DELOCALISATION_ELF,
        Menu.SPATIAL_DELOCALISATION_LOL,
        Menu.SPATIAL_DELOCALISATION_LOCAL_ENTROPY,
        Menu.SPATIAL_DELOCALISATION_ESP_TOTAL,
        Menu.SPATIAL_DELOCALISATION_RDG,
        Menu.SPATIAL_DELOCALISATION_RDG_PROMOLECULAR,
        Menu.SPATIAL_DELOCALISATION_LAMBDA2RHO,
        Menu.SPATIAL_DELOCALISATION_LAMBDA2RHO_PROMOLECULAR,
        Menu.SPATIAL_DELOCALISATION_ALIE,
        Menu.SPATIAL_DELOCALISATION_EDR,
        Menu.SPATIAL_DELOCALISATION_ORB_OVERLAP_DR,
        Menu.SPATIAL_DELOCALISATION_DELTAG_PROMOLECULAR,
        Menu.SPATIAL_DELOCALISATION_DELTAG_HIRSHFELD,
        Menu.SPATIAL_DELOCALISATION_IRI,
    ]

    FILE_EXPORT = [
        Menu.EXPORT_VARIOUS_FILES,
        Menu.SCATTER_GRAPH_TWO_FUNCTIONS,
        Menu.MONITOR_SCF_CONVERGENCE,
        Menu.GAUSSIAN_INITIAL_GUESS,
        Menu.FRAGMENT_GUESS_INPUT,
        Menu.COMBINE_FRAGMENTS,
        Menu.GENERATE_FOCK_MATRIX,
        Menu.ELECTRON_TRANSPORT_ROUTE,
        Menu.INTERMOLECULAR_OVERLAP,
        Menu.ADNDP_ANALYSIS,
    ]

    # =========================================================================
    # Helpers for listing
    # =========================================================================

    @classmethod
    def list_categories(cls) -> list[str]:
        """List available analysis categories."""
        return [category.name for category in cls]

    @classmethod
    def list_analyses(cls, category: str | None = None) -> list[Menu]:
        """List available analyses, optionally filtered by category."""
        if category is not None:
            if category not in [cat.name for cat in cls]:
                raise ValueError(f"Unknown category: {category}")
            else:
                analyses = []
                for analysis in cls:
                    if analysis.name.upper() in category.upper():
                        analyses.extend(analysis.value)
                return analyses

        else:
            analyses = []
            for analysis in cls:
                analyses.extend(analysis.value)
            return analyses

    @classmethod
    def find_category(
        cls,
        menu: Menu,
    ) -> str | None:
        """Find which category a menu belongs to."""
        for analysis in cls:
            if menu in analysis.value:
                return analysis.name
        return None
