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
        Menu.CONDENSED_FUKUI,
        Menu.ORBITAL_WEIGHTED_FUKUI,
        Menu.GRID_FUKUI,
        Menu.LOCAL_HARDNESS,
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
    ]

    ORBITAL_LOCALIZATION = [
        Menu.BOYS_LOCALIZATION_OCCUPIED,
        Menu.PIPEK_MEZEY_LOCALIZATION_HIRSHFELD_OCCUPIED,
    ]

    FUZZY_SPACE = [
        Menu.FUZZY_INTEGRATE_EDENSITY,
        Menu.FUZZY_MULTIPOLE,
        Menu.ATOMIC_OVERLAP_MATRIX,
        Menu.LOCALIZATION_DELOCALIZATION_INDEX,
        Menu.PDI_AROMATICITY,
        Menu.FLU_AROMATICITY,
        Menu.FLU_PI_AROMATICITY,
        Menu.FUZZY_OVERLAP_EDENSITY,
        Menu.CONDENSED_LINEAR_RESPONSE,
        Menu.PARA_LINEAR_RESPONSE,
        Menu.SPATIAL_DELOCALIZATION_INDEX,
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
        Menu.PRINT_DENSITY_MATRIX,
        Menu.PRINT_INTEGRAL_MATRIX_OVERLAP,
        Menu.PRINT_INTEGRAL_MATRIX_ELECTRIC_DIPOLE,
        Menu.PRINT_INTEGRAL_MATRIX_QUADRUPOLE,
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
        menu: str,
    ) -> str | None:
        """Find which category a menu belongs to."""
        for analysis in cls:
            if menu in analysis.value:
                return analysis.name
        return None
