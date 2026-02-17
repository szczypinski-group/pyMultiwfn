"""
Functional API for Multiwfn analysis.

This module provides simple functions that execute Multiwfn analyses
and return parsed results.

Usage
-----
>>> import pymultiwfn
>>> charges = pymultiwfn.hirshfeld_charges("molecule.wfn")
>>> bonds = pymultiwfn.mayer_bond_orders("molecule.wfn")
"""

# from pathlib import Path
# from typing import Any

# from pymultiwfn.config import MultiwfnConfig
# from pymultiwfn.job import MultiwfnJob
# from pymultiwfn.menu import Menu
# from pymultiwfn.parsers import (
#     BondOrderParser,
#     ChargeParser,
#     CriticalPointParser,
#     SpectrumParser,
# )

# # Module-level default configuration
# _default_config: MultiwfnConfig | None = None

# # Reusable parser instances
# _charge_parser = ChargeParser()
# _bond_order_parser = BondOrderParser()
# _cp_parser = CriticalPointParser()
# _spectrum_parser = SpectrumParser()


# # ===========================================================================
# # Configuration
# # ===========================================================================


# def set_default_config(config: MultiwfnConfig | None) -> None:
#     """Set the default configuration for all analysis functions."""
#     global _default_config
#     _default_config = config


# def get_default_config() -> MultiwfnConfig | None:
#     """Get the current default configuration."""
#     return _default_config


# def _cfg(config: MultiwfnConfig | None) -> MultiwfnConfig:
#     """Get effective config."""
#     return config or _default_config or MultiwfnConfig()


# def _run(
#     filepath: str | Path,
#     menu: Menu,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> str:
#     """Execute analysis and return stdout."""
#     job = MultiwfnJob(filepath, config=_cfg(config))
#     job.add_menu(menu)
#     return job.run(verbose=verbose).stdout


# # ===========================================================================
# # Charge Analysis (Menu 7)
# # ===========================================================================


# def hirshfeld_charges(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[int, float]:
#     """Calculate Hirshfeld charges."""
#     out = _run(filepath, Menu.HIRSHFELD_CHARGE, config, verbose)
#     return _charge_parser.parse(out, "Hirshfeld")


# def vdd_charges(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[int, float]:
#     """Calculate VDD (Voronoi Deformation Density) charges."""
#     out = _run(filepath, Menu.VDD_POPULATION, config, verbose)
#     return _charge_parser.parse(out, "VDD")


# def mulliken_charges(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[int, float]:
#     """Calculate Mulliken charges."""
#     out = _run(filepath, Menu.MULLIKEN_POPULATION, config, verbose)
#     return _charge_parser.parse(out, "Mulliken")


# def lowdin_charges(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[int, float]:
#     """Calculate Lowdin charges."""
#     out = _run(filepath, Menu.LOWDIN_POPULATION, config, verbose)
#     return _charge_parser.parse(out, "Lowdin")


# def becke_charges(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[int, float]:
#     """Calculate Becke charges."""
#     out = _run(filepath, Menu.BECKE_CHARGE, config, verbose)
#     return _charge_parser.parse(out, "Becke")


# def adch_charges(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[int, float]:
#     """Calculate ADCH (Atomic Dipole Corrected Hirshfeld) charges."""
#     out = _run(filepath, Menu.ADCH_CHARGE, config, verbose)
#     return _charge_parser.parse(out, "ADCH")


# def chelpg_charges(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[int, float]:
#     """Calculate CHELPG charges."""
#     out = _run(filepath, Menu.CHELPG_CHARGE, config, verbose)
#     return _charge_parser.parse(out, "CHELPG")


# def mk_charges(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[int, float]:
#     """Calculate Merz-Kollman (MK) charges."""
#     out = _run(filepath, Menu.MK_CHARGE, config, verbose)
#     return _charge_parser.parse(out, "MK")


# def aim_charges(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[int, float]:
#     """Calculate AIM (QTAIM) charges."""
#     out = _run(filepath, Menu.AIM_CHARGE, config, verbose)
#     return _charge_parser.parse(out, "AIM")


# def hirshfeld_i_charges(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[int, float]:
#     """Calculate Hirshfeld-I charges."""
#     out = _run(filepath, Menu.HIRSHFELD_I_CHARGE, config, verbose)
#     return _charge_parser.parse(out, "Hirshfeld-I")


# def cm5_charges(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[int, float]:
#     """Calculate CM5 charges."""
#     out = _run(filepath, Menu.CM5_CHARGE, config, verbose)
#     return _charge_parser.parse(out, "CM5")


# def eem_charges(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[int, float]:
#     """Calculate EEM charges."""
#     out = _run(filepath, Menu.EEM_CHARGE, config, verbose)
#     return _charge_parser.parse(out, "EEM")


# def resp_charges(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[int, float]:
#     """Calculate RESP charges."""
#     out = _run(filepath, Menu.RESP_CHARGE, config, verbose)
#     return _charge_parser.parse(out, "RESP")


# def gasteiger_charges(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[int, float]:
#     """Calculate Gasteiger charges."""
#     out = _run(filepath, Menu.GASTEIGER_CHARGE, config, verbose)
#     return _charge_parser.parse(out, "Gasteiger")


# def mbis_charges(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[int, float]:
#     """Calculate MBIS charges."""
#     out = _run(filepath, Menu.MBIS_CHARGE, config, verbose)
#     return _charge_parser.parse(out, "MBIS")


# # ===========================================================================
# # Bond Order Analysis (Menu 9)
# # ===========================================================================


# def mayer_bond_orders(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[tuple[int, int], float]:
#     """Calculate Mayer bond orders."""
#     out = _run(filepath, Menu.MAYER_BOND_ORDER, config, verbose)
#     return _bond_order_parser.parse(out)


# def multicenter_bond_orders(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[tuple[int, int], float]:
#     """Calculate multicenter bond orders."""
#     out = _run(filepath, Menu.MULTICENTER_BOND_ORDER, config, verbose)
#     return _bond_order_parser.parse(out)


# def wiberg_bond_orders(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[tuple[int, int], float]:
#     """Calculate Wiberg bond orders."""
#     out = _run(filepath, Menu.WIBERG_BOND_ORDER, config, verbose)
#     return _bond_order_parser.parse(out)


# def mulliken_bond_orders(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[tuple[int, int], float]:
#     """Calculate Mulliken bond orders."""
#     out = _run(filepath, Menu.MULLIKEN_BOND_ORDER, config, verbose)
#     return _bond_order_parser.parse(out)


# def fuzzy_bond_orders(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[tuple[int, int], float]:
#     """Calculate fuzzy bond orders."""
#     out = _run(filepath, Menu.FUZZY_BOND_ORDER, config, verbose)
#     return _bond_order_parser.parse(out)


# def laplacian_bond_orders(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[tuple[int, int], float]:
#     """Calculate Laplacian bond orders."""
#     out = _run(filepath, Menu.LAPLACIAN_BOND_ORDER, config, verbose)
#     return _bond_order_parser.parse(out)


# # ===========================================================================
# # Topology Analysis (Menu 2)
# # ===========================================================================


# def critical_points(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> list[dict[str, Any]]:
#     """Find QTAIM critical points."""
#     out = _run(filepath, Menu.TOPOLOGY_SEARCH_CPS, config, verbose)
#     return _cp_parser.parse(out)


# def topology_analysis(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> list[dict[str, Any]]:
#     """Run complete topology analysis."""
#     out = _run(filepath, Menu.TOPOLOGY_ANALYSIS_COMPLETE, config, verbose)
#     return _cp_parser.parse(out)


# # ===========================================================================
# # Weak Interaction Analysis (Menu 20)
# # ===========================================================================


# def nci_analysis(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> str:
#     """Run NCI (Non-Covalent Interaction) analysis."""
#     return _run(filepath, Menu.NCI_ANALYSIS, config, verbose)


# def igm_analysis(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> str:
#     """Run IGM (Independent Gradient Model) analysis."""
#     return _run(filepath, Menu.IGM_ANALYSIS, config, verbose)


# def igmh_analysis(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> str:
#     """Run IGMH analysis."""
#     return _run(filepath, Menu.IGMH_ANALYSIS, config, verbose)


# def iri_analysis(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> str:
#     """Run IRI (Interaction Region Indicator) analysis."""
#     return _run(filepath, Menu.IRI_ANALYSIS, config, verbose)


# def dori_analysis(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> str:
#     """Run DORI analysis."""
#     return _run(filepath, Menu.DORI_ANALYSIS, config, verbose)


# # ===========================================================================
# # Spectrum Analysis (Menu 11)
# # ===========================================================================


# def ir_spectrum(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[str, list[float]]:
#     """Calculate IR spectrum."""
#     out = _run(filepath, Menu.PLOT_IR_SPECTRUM, config, verbose)
#     return _spectrum_parser.parse(out)


# def raman_spectrum(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[str, list[float]]:
#     """Calculate Raman spectrum."""
#     out = _run(filepath, Menu.PLOT_RAMAN_SPECTRUM, config, verbose)
#     return _spectrum_parser.parse(out)


# def uv_vis_spectrum(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[str, list[float]]:
#     """Calculate UV-Vis spectrum."""
#     out = _run(filepath, Menu.PLOT_UV_VIS_SPECTRUM, config, verbose)
#     return _spectrum_parser.parse(out)


# def ecd_spectrum(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[str, list[float]]:
#     """Calculate ECD spectrum."""
#     out = _run(filepath, Menu.PLOT_ECD_SPECTRUM, config, verbose)
#     return _spectrum_parser.parse(out)


# def vcd_spectrum(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[str, list[float]]:
#     """Calculate VCD spectrum."""
#     out = _run(filepath, Menu.PLOT_VCD_SPECTRUM, config, verbose)
#     return _spectrum_parser.parse(out)


# def roa_spectrum(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[str, list[float]]:
#     """Calculate ROA spectrum."""
#     out = _run(filepath, Menu.PLOT_ROA_SPECTRUM, config, verbose)
#     return _spectrum_parser.parse(out)


# def nmr_spectrum(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[str, list[float]]:
#     """Calculate NMR spectrum."""
#     out = _run(filepath, Menu.PLOT_NMR_SPECTRUM, config, verbose)
#     return _spectrum_parser.parse(out)


# # ===========================================================================
# # Surface Analysis (Menu 12)
# # ===========================================================================


# def surface_esp(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> str:
#     """Analyze ESP on molecular surface."""
#     return _run(filepath, Menu.SURFACE_ANALYSIS_ESP, config, verbose)


# def surface_alie(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> str:
#     """Analyze ALIE on molecular surface."""
#     return _run(filepath, Menu.SURFACE_ANALYSIS_ALIE, config, verbose)


# def surface_area_volume(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> str:
#     """Calculate molecular surface area and volume."""
#     return _run(filepath, Menu.SURFACE_AREA_VOLUME, config, verbose)


# # ===========================================================================
# # Aromaticity Analysis (Menu 25)
# # ===========================================================================


# def icss_analysis(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> str:
#     """Run ICSS aromaticity analysis."""
#     return _run(filepath, Menu.ICSS_ANALYSIS, config, verbose)


# def nics_scan(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> str:
#     """Run NICS scan."""
#     return _run(filepath, Menu.NICS_SCAN, config, verbose)


# def homa_index(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> str:
#     """Calculate HOMA aromaticity index."""
#     return _run(filepath, Menu.HOMA_INDEX, config, verbose)


# # ===========================================================================
# # CDFT Analysis (Menu 22)
# # ===========================================================================


# def fukui_function(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> str:
#     """Calculate Fukui function."""
#     return _run(filepath, Menu.FUKUI_FUNCTION, config, verbose)


# def dual_descriptor(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> str:
#     """Calculate dual descriptor."""
#     return _run(filepath, Menu.DUAL_DESCRIPTOR, config, verbose)


# def condensed_fukui(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> str:
#     """Calculate condensed Fukui functions."""
#     return _run(filepath, Menu.CONDENSED_FUKUI, config, verbose)


# # ===========================================================================
# # Excitation Analysis (Menu 18)
# # ===========================================================================


# def hole_electron_analysis(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> str:
#     """Run hole-electron analysis."""
#     return _run(filepath, Menu.HOLE_ELECTRON_ANALYSIS, config, verbose)


# def charge_transfer_analysis(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> str:
#     """Run charge transfer analysis."""
#     return _run(filepath, Menu.CHARGE_TRANSFER_ANALYSIS, config, verbose)


# # ===========================================================================
# # DOS Analysis (Menu 10)
# # ===========================================================================


# def dos_analysis(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> str:
#     """Plot density of states."""
#     return _run(filepath, Menu.PLOT_DOS, config, verbose)


# def pdos_analysis(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> str:
#     """Plot partial density of states."""
#     return _run(filepath, Menu.PLOT_PDOS, config, verbose)


# # ===========================================================================
# # Basin Analysis (Menu 17)
# # ===========================================================================


# def basin_aim(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> str:
#     """Run AIM basin analysis."""
#     return _run(filepath, Menu.BASIN_ANALYSIS_AIM, config, verbose)


# def basin_elf(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> str:
#     """Run ELF basin analysis."""
#     return _run(filepath, Menu.BASIN_ANALYSIS_ELF, config, verbose)


# # ===========================================================================
# # Fuzzy Atomic Space (Menu 15)
# # ===========================================================================


# def localization_index(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> str:
#     """Calculate localization/delocalization index."""
#     return _run(
#         filepath, Menu.LOCALIZATION_DELOCALIZATION_INDEX, config, verbose
#     )


# def atomic_dipoles(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> str:
#     """Calculate atomic dipole moments."""
#     return _run(filepath, Menu.ATOMIC_DIPOLE_MOMENTS, config, verbose)


# # ===========================================================================
# # EDA Analysis (Menu 21)
# # ===========================================================================


# def eda_analysis(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> str:
#     """Run EDA-FF analysis."""
#     return _run(filepath, Menu.EDA_FF, config, verbose)


# # ===========================================================================
# # Wavefunction Info (Menu 6)
# # ===========================================================================


# def wavefunction_info(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> str:
#     """Check wavefunction information."""
#     return _run(filepath, Menu.CHECK_WAVEFUNCTION, config, verbose)


# def orbital_info(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> str:
#     """Print orbital information."""
#     return _run(filepath, Menu.PRINT_ORBITAL_INFO, config, verbose)


# # ===========================================================================
# # Convinience Functions(runnning multiple categories of calculations at once)
# # ===========================================================================


# def all_topology(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> list[dict[str, Any]]:
#     """Run complete topology analysis."""
#     return topology_analysis(filepath, config, verbose)


# def all_planeMap(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[str, list[float]]:
#     """Calculate all plane maps."""
#     return {
#         "density": surface_esp(filepath, config, verbose),
#         "esp": surface_esp(filepath, config, verbose),
#         "elf": surface_alie(filepath, config, verbose),
#         "lol": surface_alie(filepath, config, verbose),
#         "gradient": surface_alie(filepath, config, verbose),
#         "laplacian": surface_alie(filepath, config, verbose),
#     }


# def all_cube(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[str, list[float]]:
#     """Calculate all cubes."""
#     return {
#         "density": surface_esp(filepath, config, verbose),
#         "spin_density": surface_esp(filepath, config, verbose),
#         "elf": surface_alie(filepath, config, verbose),
#         "lol": surface_alie(filepath, config, verbose),
#         "esp": surface_alie(filepath, config, verbose),
#         "laplacian": surface_alie(filepath, config, verbose),
#         "fukui_minus": surface_alie(filepath, config, verbose),
#         "fukui_plus": surface_alie(filepath, config, verbose),
#         "dual_descriptor": surface_alie(filepath, config, verbose),
#     }


# def all_populations(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[str, dict[int, float]]:
#     """Calculate all population analyses."""
#     return {
#         "hirshfeld": hirshfeld_charges(filepath, config, verbose),
#         "vdd": vdd_charges(filepath, config, verbose),
#         "mulliken": mulliken_charges(filepath, config, verbose),
#         "lowdin": lowdin_charges(filepath, config, verbose),
#         "becke": becke_charges(filepath, config, verbose),
#         "adch": adch_charges(filepath, config, verbose),
#         "chelpg": chelpg_charges(filepath, config, verbose),
#         "mk": mk_charges(filepath, config, verbose),
#         "aim": aim_charges(filepath, config, verbose),
#         "hirshfeld_i": hirshfeld_i_charges(filepath, config, verbose),
#         "cm5": cm5_charges(filepath, config, verbose),
#         "eem": eem_charges(filepath, config, verbose),
#         "resp": resp_charges(filepath, config, verbose),
#         "gasteiger": gasteiger_charges(filepath, config, verbose),
#         "mbis": mbis_charges(filepath, config, verbose),
#     }


# def all_bondOrders(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[str, dict[tuple[int, int], float]]:
#     """Calculate all bond orders."""
#     return {
#         "mayer": mayer_bond_orders(filepath, config, verbose),
#         "multicenter": multicenter_bond_orders(filepath, config, verbose),
#         "wiberg": wiberg_bond_orders(filepath, config, verbose),
#         "mulliken": mulliken_bond_orders(filepath, config, verbose),
#         "fuzzy": fuzzy_bond_orders(filepath, config, verbose),
#         "laplacian": laplacian_bond_orders(filepath, config, verbose),
#     }


# def all_dos(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[str, str]:
#     """Calculate all DOS analyses."""
#     return {
#         "dos": dos_analysis(filepath, config, verbose),
#         "pdos": pdos_analysis(filepath, config, verbose),
#     }


# def all_spectra(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[str, dict[str, list[float]]]:
#     """Calculate all spectra."""
#     return {
#         "ir": ir_spectrum(filepath, config, verbose),
#         "raman": raman_spectrum(filepath, config, verbose),
#         "uv_vis": uv_vis_spectrum(filepath, config, verbose),
#         "ecd": ecd_spectrum(filepath, config, verbose),
#         "vcd": vcd_spectrum(filepath, config, verbose),
#         "roa": roa_spectrum(filepath, config, verbose),
#         "nmr": nmr_spectrum(filepath, config, verbose),
#     }


# def all_surfaceAnalysis(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[str, str]:
#     """Calculate all surface analyses."""
#     return {
#         "esp": surface_esp(filepath, config, verbose),
#         "alie": surface_alie(filepath, config, verbose),
#         "area_volume": surface_area_volume(filepath, config, verbose),
#     }


# def all_gridProcessing(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[str, str]:
#     """Calculate all grid processing analyses."""
#     return {
#         "icss": icss_analysis(filepath, config, verbose),
#         "nics_scan": nics_scan(filepath, config, verbose),
#         "homa_index": homa_index(filepath, config, verbose),
#     }


# def all_fuzzyAtomicSpace(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[str, str]:
#     """Calculate all fuzzy atomic space analyses."""
#     return {
#         "localization_index": localization_index(filepath, config, verbose),
#         "atomic_dipoles": atomic_dipoles(filepath, config, verbose),
#     }


# def all_basinAnalysis(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[str, str]:
#     """Calculate all basin analyses."""
#     return {
#         "aim": basin_aim(filepath, config, verbose),
#         "elf": basin_elf(filepath, config, verbose),
#     }


# def all_excitationAnalysis(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[str, str]:
#     """Calculate all excitation analyses."""
#     return {
#         "hole_electron": hole_electron_analysis(filepath, config, verbose),
#         "charge_transfer": charge_transfer_analysis(
#           filepath, config, verbose),
#     }


# def all_localisationAnalysis(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[str, str]:
#     """Calculate all localization analyses."""
#     return {
#         "fukui_function": fukui_function(filepath, config, verbose),
#         "dual_descriptor": dual_descriptor(filepath, config, verbose),
#         "condensed_fukui": condensed_fukui(filepath, config, verbose),
#     }


# def all_weakInteractions(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[str, str]:
#     """Calculate all weak interaction analyses."""
#     return {
#         "nci": nci_analysis(filepath, config, verbose),
#         "igm": igm_analysis(filepath, config, verbose),
#         "igmh": igmh_analysis(filepath, config, verbose),
#         "iri": iri_analysis(filepath, config, verbose),
#         "dori": dori_analysis(filepath, config, verbose),
#     }


# def all_eda(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[str, str]:
#     """Calculate all EDA analyses."""
#     return {
#         "eda": eda_analysis(filepath, config, verbose),
#     }


# def all_cfdt(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[str, str]:
#     """Calculate all CDFT analyses."""
#     return {
#         "fukui_function": fukui_function(filepath, config, verbose),
#         "dual_descriptor": dual_descriptor(filepath, config, verbose),
#         "condensed_fukui": condensed_fukui(filepath, config, verbose),
#     }


# def all_polarisability(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[str, str]:
#     """Calculate all polarizability analyses."""
#     return {
#         "ir": ir_spectrum(filepath, config, verbose),
#         "raman": raman_spectrum(filepath, config, verbose),
#     }


# def all_aromaticity(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[str, str]:
#     """Calculate all aromaticity analyses."""
#     return {
#         "icss": icss_analysis(filepath, config, verbose),
#         "nics_scan": nics_scan(filepath, config, verbose),
#         "homa_index": homa_index(filepath, config, verbose),
#     }

# def all_utilities1(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[str, str]:
#     """Calculate all utility analyses."""
#     return {
#         "wavefunction_info": wavefunction_info(filepath, config, verbose),
#         "orbital_info": orbital_info(filepath, config, verbose),
#     }


# def all_utilities2(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[str, str]:
#     """Calculate all utility analyses."""
#     return {
#         "wavefunction_info": wavefunction_info(filepath, config, verbose),
#         "orbital_info": orbital_info(filepath, config, verbose),
#     }


# # TO DO: make sure functions have relevant calculations
# def all_utilities3(
#     filepath: str | Path,
#     config: MultiwfnConfig | None = None,
#     verbose: bool = False,
# ) -> dict[str, str]:
#     """Calculate all utility analyses."""
#     return {
#         "wavefunction_info": wavefunction_info(filepath, config, verbose),
#         "orbital_info": orbital_info(filepath, config, verbose),
#     }


# # =========================================================================
# # Exports
# # =========================================================================

# __all__ = [
#     # Config
#     "set_default_config",
#     "get_default_config",
#     # Charges
#     "hirshfeld_charges",
#     "vdd_charges",
#     "mulliken_charges",
#     "lowdin_charges",
#     "becke_charges",
#     "adch_charges",
#     "chelpg_charges",
#     "mk_charges",
#     "aim_charges",
#     "hirshfeld_i_charges",
#     "cm5_charges",
#     "eem_charges",
#     "resp_charges",
#     "gasteiger_charges",
#     "mbis_charges",
#     # Bond orders
#     "mayer_bond_orders",
#     "multicenter_bond_orders",
#     "wiberg_bond_orders",
#     "mulliken_bond_orders",
#     "fuzzy_bond_orders",
#     "laplacian_bond_orders",
#     # Topology
#     "critical_points",
#     "topology_analysis",
#     # Weak interactions
#     "nci_analysis",
#     "igm_analysis",
#     "igmh_analysis",
#     "iri_analysis",
#     "dori_analysis",
#     # Spectra
#     "ir_spectrum",
#     "raman_spectrum",
#     "uv_vis_spectrum",
#     "ecd_spectrum",
#     "vcd_spectrum",
#     "roa_spectrum",
#     "nmr_spectrum",
#     # Surface
#     "surface_esp",
#     "surface_alie",
#     "surface_area_volume",
#     # Aromaticity
#     "icss_analysis",
#     "nics_scan",
#     "homa_index",
#     # CDFT
#     "fukui_function",
#     "dual_descriptor",
#     "condensed_fukui",
#     # Excitation
#     "hole_electron_analysis",
#     "charge_transfer_analysis",
#     # DOS
#     "dos_analysis",
#     "pdos_analysis",
#     # Basin
#     "basin_aim",
#     "basin_elf",
#     # Fuzzy
#     "localization_index",
#     "atomic_dipoles",
#     # EDA
#     "eda_analysis",
#     # Wavefunction
#     "wavefunction_info",
#     "orbital_info",
# ]
