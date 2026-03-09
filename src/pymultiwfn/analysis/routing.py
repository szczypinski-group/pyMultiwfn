"""Routing for parsers based on Multiwfn menu options."""

from typing import Any

from pymultiwfn.analysis.parsers import (
    AromaticityParser,
    BasinParser,
    BondOrderParser,
    CDFTParser,
    ChargeParser,
    CriticalPointParser,
    CubeParser,
    DOSParser,
    EDAParser,
    ExcitationParser,
    FuzzySpaceParser,
    OrbitalCompositionParser,
    PolarizabilityParser,
    SpectrumParser,
    SurfaceParser,
    UtilityParser,
    WavefunctionParser,
    WeakInteractionParser,
)
from pymultiwfn.enums.menu import Menu


class ParserRoute:
    """Routing for Multiwfn output parsers.

    This class manages a routing table that maps Multiwfn menu options to the
    appropriate parser functions. It also provides helper methods that call
    a combination of parsers from the `analysis.parsers.py` module.

    """

    @staticmethod
    def _parse_charges(
        stdout: str,
    ) -> dict[str, dict[str, float] | dict[int, float]]:
        """Parse charges and optionally the dipole moment."""
        charges = ChargeParser.parse(stdout)
        result: dict[str, dict[str, float] | dict[int, float]] = {
            "charges": charges,
        }
        dipole = ChargeParser.parse_dipole(stdout)
        if dipole is not None:
            result["dipole"] = dipole
        return result

    @staticmethod
    def _parse_bond_orders(
        stdout: str,
    ) -> dict[str, dict[tuple[int, int], float] | dict[int, dict[str, float]]]:
        bos = BondOrderParser.parse(stdout)
        valences = BondOrderParser.parse_valence(stdout)
        result: dict[
            str, dict[tuple[int, int], float] | dict[int, dict[str, float]]
        ] = {
            "bond_orders": bos,
        }
        if valences:
            result["valences"] = valences
        return result

    @staticmethod
    def _parse_multicenter_bo(
        stdout: str,
    ) -> dict[str, list[dict[str, list[int] | float]]]:
        mc = BondOrderParser.parse_multicenter(stdout)
        return {"multicenter_bond_orders": mc}

    @staticmethod
    def _parse_bo_decomposition(
        stdout: str,
    ) -> dict[str, list[dict[str, int | float]]]:
        decomp = BondOrderParser.parse_decomposition(stdout)
        return {"bond_order_decomposition": decomp}

    @staticmethod
    def _parse_critical_points(
        stdout: str,
    ) -> dict[str, int | float | str | tuple[float, float, float] | None]:
        cps = CriticalPointParser.parse(stdout)
        paths = CriticalPointParser.parse_bond_paths(stdout)
        summary = CriticalPointParser.summary(cps)
        result: dict[str, Any] = {
            "critical_points": cps,
            "summary": summary,
        }
        if paths:
            result["bond_paths"] = paths
        return result

    @staticmethod
    def _parse_spectrum(stdout: str) -> dict[str, Any]:
        peaks = SpectrumParser.parse(stdout)
        transitions = SpectrumParser.parse_transitions(stdout)
        result: dict[str, Any] = {"peaks": peaks}
        if transitions:
            result["transitions"] = transitions
        color = SpectrumParser.parse_color(stdout)
        if color is not None:
            result["predicted_color"] = color
        return result

    @staticmethod
    def _parse_dos(stdout: str) -> dict[str, Any]:
        data = DOSParser.parse(stdout)
        orbitals = DOSParser.parse_orbital_energies(stdout)
        result: dict[str, Any] = {"dos_curves": data}
        if orbitals:
            result["orbital_energies"] = orbitals
        return result

    @staticmethod
    def _parse_surface(stdout: str) -> dict[str, Any]:
        stats = SurfaceParser.parse(stdout)
        extrema = SurfaceParser.parse_extrema(stdout)
        result: dict[str, Any] = {"surface_statistics": stats}
        if extrema:
            result["surface_extrema"] = extrema
        return result

    @staticmethod
    def _parse_orbital_composition(stdout: str) -> dict[str, Any]:
        orbitals = OrbitalCompositionParser.parse(stdout)
        return {"orbital_compositions": orbitals}

    @staticmethod
    def _parse_oxidation_states(stdout: str) -> dict[str, Any]:
        states = OrbitalCompositionParser.parse_oxidation_states(stdout)
        return {"oxidation_states": {str(k): v for k, v in states.items()}}

    @staticmethod
    def _parse_fuzzy_properties(stdout: str) -> dict[str, Any]:
        props = FuzzySpaceParser.parse_atomic_properties(stdout)
        return {
            "fuzzy_atomic_properties": {str(k): v for k, v in props.items()},
        }

    @staticmethod
    def _parse_delocalization(stdout: str) -> dict[str, Any]:
        data = FuzzySpaceParser.parse_delocalization_indices(stdout)
        result: dict[str, Any] = {}
        if data.get("localization"):
            result["localization_indices"] = {
                str(k): v for k, v in data["localization"].items()
            }
        if data.get("delocalization"):
            result["delocalization_indices"] = {
                f"{a}-{b}": v for (a, b), v in data["delocalization"].items()
            }
        return result

    @staticmethod
    def _parse_fuzzy_aromaticity(stdout: str) -> dict[str, Any]:
        return {
            "aromaticity_indices": FuzzySpaceParser.parse_aromaticity_index(
                stdout
            )
        }

    @staticmethod
    def _parse_basin(stdout: str) -> dict[str, Any]:
        basins = BasinParser.parse(stdout)
        charges = BasinParser.parse_charges(stdout)
        result: dict[str, Any] = {"basins": basins}
        if charges:
            result["basin_charges"] = {str(k): v for k, v in charges.items()}
        return result

    @staticmethod
    def _parse_hole_electron(stdout: str) -> dict[str, Any]:
        return {"hole_electron": ExcitationParser.parse_hole_electron(stdout)}

    @staticmethod
    def _parse_charge_transfer(stdout: str) -> dict[str, Any]:
        return {
            "charge_transfer": ExcitationParser.parse_charge_transfer(stdout)
        }

    @staticmethod
    def _parse_delta_r(stdout: str) -> dict[str, Any]:
        return {"delta_r": ExcitationParser.parse_delta_r(stdout)}

    @staticmethod
    def _parse_lambda_index(stdout: str) -> dict[str, Any]:
        return {"lambda_index": ExcitationParser.parse_lambda_index(stdout)}

    @staticmethod
    def _parse_weak_interactions(stdout: str) -> dict[str, Any]:
        return {
            "weak_interaction_summary": WeakInteractionParser.parse(stdout)
        }

    @staticmethod
    def _parse_eda(stdout: str) -> dict[str, Any]:
        components = EDAParser.parse(stdout)
        disp = EDAParser.parse_dispersion_contributions(stdout)
        result: dict[str, Any] = {"eda_components": components}
        if disp:
            result["dispersion_contributions"] = {
                str(k): v for k, v in disp.items()
            }
        return result

    @staticmethod
    def _parse_cdft_global(stdout: str) -> dict[str, Any]:
        return {"cdft_global_indices": CDFTParser.parse_global_indices(stdout)}

    @staticmethod
    def _parse_condensed_fukui(stdout: str) -> dict[str, Any]:
        fukui = CDFTParser.parse_condensed_fukui(stdout)
        return {"condensed_fukui": {str(k): v for k, v in fukui.items()}}

    @staticmethod
    def _parse_dual_descriptor(stdout: str) -> dict[str, Any]:
        dd = CDFTParser.parse_dual_descriptor(stdout)
        return {"dual_descriptor": {str(k): v for k, v in dd.items()}}

    @staticmethod
    def _parse_polarizability(stdout: str) -> dict[str, Any]:
        return {"polarizability": PolarizabilityParser.parse(stdout)}

    @staticmethod
    def _parse_aromaticity(stdout: str) -> dict[str, Any]:
        indices = AromaticityParser.parse(stdout)
        scan = AromaticityParser.parse_nics_scan(stdout)
        result: dict[str, Any] = {"aromaticity": indices}
        if scan.get("distances"):
            result["nics_scan"] = scan
        return result

    @staticmethod
    def _parse_wavefunction_info(stdout: str) -> dict[str, Any]:
        return {"orbital_info": WavefunctionParser.parse_orbital_info(stdout)}

    @staticmethod
    def _parse_cube_info(stdout: str) -> dict[str, Any]:
        return {"cube_info": CubeParser.parse(stdout)}

    @staticmethod
    def _parse_geometry(stdout: str) -> dict[str, Any]:
        return {"geometry": UtilityParser.parse_geometry(stdout)}

    @staticmethod
    def _parse_multipole_moments(stdout: str) -> dict[str, Any]:
        return {
            "multipole_moments": UtilityParser.parse_multipole_moments(stdout)
        }

    @staticmethod
    def _parse_bla_boa(stdout: str) -> dict[str, Any]:
        return {"bla_boa": UtilityParser.parse_bla_boa(stdout)}

    ROUTE_TABLE = {
        # Menu 7 — charges
        Menu.HIRSHFELD_CHARGE: _parse_charges,
        Menu.VDD_POPULATION: _parse_charges,
        Menu.MULLIKEN_POPULATION: _parse_charges,
        Menu.LOWDIN_POPULATION: _parse_charges,
        Menu.SCPA_POPULATION: _parse_charges,
        Menu.STOUT_POLITZER_POPULATION: _parse_charges,
        Menu.BICKELHAUPT_POPULATION: _parse_charges,
        Menu.BECKE_CHARGE: _parse_charges,
        Menu.ADCH_CHARGE: _parse_charges,
        Menu.CHELPG_CHARGE: _parse_charges,
        Menu.MK_CHARGE: _parse_charges,
        Menu.AIM_CHARGE: _parse_charges,
        Menu.HIRSHFELD_I_CHARGE: _parse_charges,
        Menu.CM5_CHARGE: _parse_charges,
        Menu.EEM_CHARGE: _parse_charges,
        Menu.RESP_CHARGE: _parse_charges,
        Menu.GASTEIGER_CHARGE: _parse_charges,
        Menu.MBIS_CHARGE: _parse_charges,
        Menu.DDEC_CHARGE: _parse_charges,
        # Menu 8 — orbital composition
        Menu.ORBITAL_COMPOSITION_MULLIKEN: _parse_orbital_composition,
        Menu.ORBITAL_COMPOSITION_SCPA: _parse_orbital_composition,
        Menu.ORBITAL_COMPOSITION_STOUT_POLITZER: _parse_orbital_composition,
        Menu.ORBITAL_COMPOSITION_FRAGMENT_MULLIKEN: _parse_orbital_composition,
        Menu.ORBITAL_COMPOSITION_FRAGMENT_STOUT: _parse_orbital_composition,
        Menu.ORBITAL_COMPOSITION_FRAGMENT_SCPA: _parse_orbital_composition,
        Menu.ORBITAL_COMPOSITION_NAO: _parse_orbital_composition,
        Menu.ORBITAL_COMPOSITION_HIRSHFELD: _parse_orbital_composition,
        Menu.ORBITAL_COMPOSITION_BECKE: _parse_orbital_composition,
        Menu.LOBA_OXIDATION_STATE: _parse_oxidation_states,
        # Menu 9 — bond orders
        Menu.MAYER_BOND_ORDER: _parse_bond_orders,
        Menu.WIBERG_BOND_ORDER: _parse_bond_orders,
        Menu.MULLIKEN_BOND_ORDER: _parse_bond_orders,
        Menu.FUZZY_BOND_ORDER: _parse_bond_orders,
        Menu.LAPLACIAN_BOND_ORDER: _parse_bond_orders,
        Menu.IBSI_ANALYSIS: _parse_bond_orders,
        Menu.MULTICENTER_BOND_ORDER: _parse_multicenter_bo,
        Menu.MULTICENTER_BOND_ORDER_NAO: _parse_multicenter_bo,
        Menu.MULLIKEN_BOND_ORDER_DECOMPOSE: _parse_bo_decomposition,
        Menu.ORBITAL_PERTURBED_MAYER: _parse_bo_decomposition,
        Menu.WIBERG_DECOMPOSITION: _parse_bo_decomposition,
        Menu.AV1245_INDEX: _parse_fuzzy_aromaticity,
        # Menu 2 — topology
        Menu.TOPOLOGY_SEARCH_CPS: _parse_critical_points,
        Menu.TOPOLOGY_GENERATE_PATHS: _parse_critical_points,
        Menu.TOPOLOGY_ANALYSIS_COMPLETE: _parse_critical_points,
        Menu.TOPOLOGY_ESP_ANALYSIS: _parse_critical_points,
        Menu.TOPOLOGY_LOL_ANALYSIS: _parse_critical_points,
        Menu.TOPOLOGY_ELF_ANALYSIS: _parse_critical_points,
        Menu.TOPOLOGY_LAPLACIAN_ANALYSIS: _parse_critical_points,
        Menu.TOPOLOGY_SEARCH_BCP: _parse_critical_points,
        Menu.TOPOLOGY_SEARCH_RCP: _parse_critical_points,
        Menu.TOPOLOGY_SEARCH_CCP: _parse_critical_points,
        # Menu 10 — density of states
        Menu.PLOT_DOS: _parse_dos,
        Menu.PLOT_PDOS: _parse_dos,
        Menu.PLOT_OPDOS: _parse_dos,
        Menu.PLOT_LDOS: _parse_dos,
        Menu.PLOT_PHOTOELECTRON_SPECTRUM: _parse_dos,
        Menu.PLOT_COHP: _parse_dos,
        # Menu 11 — spectra
        Menu.PLOT_IR_SPECTRUM: _parse_spectrum,
        Menu.PLOT_RAMAN_SPECTRUM: _parse_spectrum,
        Menu.PLOT_UV_VIS_SPECTRUM: _parse_spectrum,
        Menu.PLOT_ECD_SPECTRUM: _parse_spectrum,
        Menu.PLOT_VCD_SPECTRUM: _parse_spectrum,
        Menu.PLOT_ROA_SPECTRUM: _parse_spectrum,
        Menu.PLOT_NMR_SPECTRUM: _parse_spectrum,
        Menu.PLOT_FLUORESCENCE_SPECTRUM: _parse_spectrum,
        Menu.PLOT_PVS: _parse_spectrum,
        Menu.PREDICT_COLOR: _parse_spectrum,
        # Menu 12 — surface analysis
        Menu.SURFACE_ANALYSIS_ESP: _parse_surface,
        Menu.SURFACE_ANALYSIS_ALIE: _parse_surface,
        Menu.SURFACE_AREA_VOLUME: _parse_surface,
        Menu.BECKE_SURFACE: _parse_surface,
        Menu.HIRSHFELD_SURFACE: _parse_surface,
        Menu.SURFACE_EXTREMA: _parse_surface,
        Menu.HIRSHFELD_SURFACE_FINGERPRINT: _parse_surface,
        # Menu 15 — fuzzy atomic space
        Menu.FUZZY_INTEGRATE_PROPERTY: _parse_fuzzy_properties,
        Menu.ATOMIC_DIPOLE_MOMENTS: _parse_fuzzy_properties,
        Menu.ATOMIC_OVERLAP_MATRIX: _parse_fuzzy_properties,
        Menu.ATOMIC_VOLUME_POLARIZABILITY: _parse_fuzzy_properties,
        Menu.LOCALIZATION_DELOCALIZATION_INDEX: _parse_delocalization,
        Menu.PDI_AROMATICITY: _parse_fuzzy_aromaticity,
        Menu.FLU_AROMATICITY: _parse_fuzzy_aromaticity,
        Menu.FLU_PI_AROMATICITY: _parse_fuzzy_aromaticity,
        Menu.MULTICENTER_DI: _parse_fuzzy_aromaticity,
        Menu.ITA_AROMATICITY: _parse_fuzzy_aromaticity,
        Menu.CONDENSED_LINEAR_RESPONSE: _parse_fuzzy_properties,
        Menu.PARA_LINEAR_RESPONSE: _parse_fuzzy_aromaticity,
        Menu.IFDI_ANALYSIS: _parse_delocalization,
        Menu.FUZZY_INTEGRATE_OVERLAP: _parse_fuzzy_properties,
        # Menu 17 — basin analysis
        Menu.BASIN_ANALYSIS_AIM: _parse_basin,
        Menu.BASIN_ANALYSIS_ELF: _parse_basin,
        Menu.BASIN_INTEGRATE_PROPERTY: _parse_basin,
        Menu.BASIN_ANALYSIS_ESP: _parse_basin,
        Menu.BASIN_ANALYSIS_LOL: _parse_basin,
        Menu.BASIN_ANALYSIS_LOL_ALPHA: _parse_basin,
        Menu.BASIN_ANALYSIS_CUSTOM: _parse_basin,
        # Menu 18 — excitation analysis
        Menu.HOLE_ELECTRON_ANALYSIS: _parse_hole_electron,
        Menu.CHARGE_TRANSFER_ANALYSIS: _parse_charge_transfer,
        Menu.DELTA_R_INDEX: _parse_delta_r,
        Menu.LAMBDA_INDEX: _parse_lambda_index,
        Menu.IFCT_ANALYSIS: _parse_charge_transfer,
        Menu.CTS_ANALYSIS: _parse_charge_transfer,
        # Menu 20 — weak interactions
        Menu.NCI_ANALYSIS: _parse_weak_interactions,
        Menu.NCI_PROMOLECULAR: _parse_weak_interactions,
        Menu.ANCI_ANALYSIS: _parse_weak_interactions,
        Menu.IRI_ANALYSIS: _parse_weak_interactions,
        Menu.DORI_ANALYSIS: _parse_weak_interactions,
        Menu.VDW_POTENTIAL: _parse_weak_interactions,
        Menu.IGM_ANALYSIS: _parse_weak_interactions,
        Menu.IGMH_ANALYSIS: _parse_weak_interactions,
        Menu.AIGM_ANALYSIS: _parse_weak_interactions,
        Menu.MIGM_ANALYSIS: _parse_weak_interactions,
        Menu.AMIGM_ANALYSIS: _parse_weak_interactions,
        # Menu 21 — EDA
        Menu.EDA_FF: _parse_eda,
        Menu.EDA_SBL: _parse_eda,
        Menu.SOBEDA_ANALYSIS: _parse_eda,
        Menu.DISPERSION_ATOMIC_CONTRIBUTION: _parse_eda,
        # Menu 22 — CDFT
        Menu.CDFT_ANALYSIS: _parse_cdft_global,
        Menu.FUKUI_FUNCTION: _parse_condensed_fukui,
        Menu.DUAL_DESCRIPTOR: _parse_dual_descriptor,
        Menu.CONDENSED_FUKUI: _parse_condensed_fukui,
        Menu.LOCAL_HARDNESS: _parse_cdft_global,
        Menu.LOCAL_IONIZATION_ENERGY: _parse_cdft_global,
        # Menu 24 — polarizability
        Menu.PARSE_POLARIZABILITY: _parse_polarizability,
        Menu.SOS_POLARIZABILITY: _parse_polarizability,
        Menu.POLARIZABILITY_DENSITY: _parse_polarizability,
        Menu.UNIT_SPHERE_POLARIZABILITY: _parse_polarizability,
        # Menu 25 — aromaticity
        Menu.AICD_ANALYSIS: _parse_aromaticity,
        Menu.NICS_POINT: _parse_aromaticity,
        Menu.ICSS_ANALYSIS: _parse_aromaticity,
        Menu.NICS_SCAN: _parse_aromaticity,
        Menu.BIRD_INDEX: _parse_aromaticity,
        Menu.HOMA_INDEX: _parse_aromaticity,
        Menu.HOMAC_HOMER: _parse_aromaticity,
        Menu.STANGER_INDEX: _parse_aromaticity,
        Menu.NICS_1D_SCAN: _parse_aromaticity,
        Menu.NICS_2D_MAP: _parse_aromaticity,
        # Menu 6 — wavefunction info
        Menu.PRINT_ORBITAL_INFO: _parse_wavefunction_info,
        Menu.PRINT_GTF_INFO: _parse_wavefunction_info,
        Menu.PRINT_BASIS_INFO: _parse_wavefunction_info,
        # Menu 5 — cube generation (store metadata only, not the cube)
        Menu.CUBE_DENSITY: _parse_cube_info,
        Menu.CUBE_SPIN_DENSITY: _parse_cube_info,
        Menu.CUBE_ELF: _parse_cube_info,
        Menu.CUBE_LOL: _parse_cube_info,
        Menu.CUBE_ESP: _parse_cube_info,
        Menu.CUBE_LAPLACIAN: _parse_cube_info,
        Menu.CUBE_GRADIENT_NORM: _parse_cube_info,
        Menu.CUBE_KINETIC_G: _parse_cube_info,
        Menu.CUBE_KINETIC_K: _parse_cube_info,
        Menu.CUBE_ALIE: _parse_cube_info,
        Menu.CUBE_RDG: _parse_cube_info,
        Menu.CUBE_SIGN_LAMBDA2_RHO: _parse_cube_info,
        Menu.CUBE_SOURCE_FUNCTION: _parse_cube_info,
        Menu.CUBE_ORBITAL_WAVEFUNCTION: _parse_cube_info,
        Menu.CUBE_FUKUI_MINUS: _parse_cube_info,
        Menu.CUBE_FUKUI_PLUS: _parse_cube_info,
        Menu.CUBE_DUAL_DESCRIPTOR: _parse_cube_info,
        Menu.CUBE_PROMOLECULAR_DENSITY: _parse_cube_info,
        Menu.CUBE_DEFORMATION_DENSITY: _parse_cube_info,
        Menu.CUBE_DENSITY_HIGH: _parse_cube_info,
        Menu.CUBE_ESP_HIGH: _parse_cube_info,
        Menu.CUBE_ELF_HIGH: _parse_cube_info,
        # Menu 100 — utilities
        Menu.GEOMETRY_PROPERTIES: _parse_geometry,
        Menu.ELECTRIC_MULTIPOLE_MOMENTS: _parse_multipole_moments,
        Menu.BLA_BOA_ANALYSIS: _parse_bla_boa,
    }
