"""Common analysis sets for pyMultiwfn."""

TOPOLOGY_ANALYSES = frozenset(
    {
        "topology_search_cps",
        "topology_generate_paths",
        "topology_interbasin_surfaces",
        "topology_analysis_complete",
        "topology_esp_analysis",
        "topology_lol_analysis",
        "complete_qtaim_analysis",
        "basin_analysis_aim",
        "basin_analysis_elf",
    }
)

CUBE_ANALYSES = frozenset(
    {
        "cube_density",
        "cube_spin_density",
        "cube_elf",
        "cube_lol",
        "cube_esp",
        "cube_laplacian",
        "cube_fukui_minus",
        "cube_fukui_plus",
        "cube_dual_descriptor",
        "cube_orbital",
        "batch_cube_generation",
    }
)

BATCH_ANALYSES = frozenset(
    {
        "run_all",
        "run_all_charges",
        "run_all_bond_orders",
        "run_all_weak_interactions",
        "run_all_topology",
        "run_all_cubes",
        "run_all_spectra",
        "run_all_surfaces",
        "run_all_aromaticity",
        "run_all_cdft",
        "weak_interaction_suite",
    }
)
