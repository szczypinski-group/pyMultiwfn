"""
Complete Multiwfn menu mapping for batch automation.
Based on Multiwfn 3.8 (dev) manual.
Each function returns a menu sequence (list of strings) compatible with run_multiwfn.
"""

# =========================================================================
# Main Menu 0: Showing molecular structure and viewing orbitals/isosurfaces
# =========================================================================
def view_structure():
    """View molecular structure in GUI"""
    return ["0", "q"]

def view_orbital(orbital_index=None):
    """View specific orbital or all orbitals
    
    Args:
        orbital_index: Orbital number to view, or None to view all
    """
    if orbital_index is None:
        return ["0", "q"]
    return ["0", str(orbital_index), "q"]

# =========================================================================
# Main Menu 1: Output all properties at a point
# =========================================================================
def properties_at_point(x, y, z):
    """Output all properties at specified coordinates
    
    Args:
        x, y, z: Coordinates in Angstrom
    """
    return ["1", f"{x},{y},{z}", "q"]

# =========================================================================
# Main Menu 2: Topology analysis
# =========================================================================
def topology_search_cps():
    """Search for critical points"""
    return ["2", "2", "-1", "q"]

def topology_generate_paths():
    """Generate topology paths"""
    return ["2", "3", "q"]

def topology_interbasin_surfaces():
    """Generate interbasin surfaces"""
    return ["2", "4", "q"]

def topology_analysis_complete():
    """Complete topology analysis (search CPs + paths + surfaces)"""
    return ["2", "2", "-1", "3", "4", "q"]

def topology_esp_analysis():
    """Topology analysis of electrostatic potential"""
    return ["2", "-2", "2", "-1", "q"]

def topology_lol_analysis():
    """Topology analysis of LOL (Localized Orbital Locator)"""
    return ["2", "-10", "2", "-1", "q"]

# =========================================================================
# Main Menu 3: Output and plot specific property in a line
# =========================================================================
def line_scan(func_type="1"):
    """Plot property along a line
    
    Args:
        func_type: Function type (1=density, 2=gradient, etc.)
    """
    return ["3", func_type, "q"]

def line_esp():
    """ESP along a line"""
    return ["3", "12", "q"]

# =========================================================================
# Main Menu 4: Output and plot specific property in a plane
# =========================================================================
def plane_map_density():
    """Electron density plane map"""
    return ["4", "1", "1", "q"]

def plane_map_esp():
    """ESP plane map"""
    return ["4", "12", "1", "q"]

def plane_map_elf():
    """ELF plane map"""
    return ["4", "4", "1", "q"]

def plane_map_lol():
    """LOL plane map"""
    return ["4", "5", "1", "q"]

def plane_map_gradient():
    """Gradient norm plane map"""
    return ["4", "2", "1", "q"]

def plane_map_laplacian():
    """Laplacian plane map"""
    return ["4", "3", "1", "q"]

def plane_map_custom(func_type, plane_type="1"):
    """Custom plane map
    
    Args:
        func_type: Function type code
        plane_type: Plane type (1=filled, 2=contour, etc.)
    """
    return ["4", str(func_type), str(plane_type), "q"]

# =========================================================================
# Main Menu 5: Output and plot property within spatial region (cube files)
# =========================================================================
def cube_density():
    """Generate electron density cube file"""
    return ["5", "1", "2", "0", "q"]

def cube_spin_density():
    """Generate spin density cube file"""
    return ["5", "2", "2", "0", "q"]

def cube_elf():
    """Generate ELF cube file"""
    return ["5", "4", "2", "0", "q"]

def cube_lol():
    """Generate LOL cube file"""
    return ["5", "5", "2", "0", "q"]

def cube_esp():
    """Generate ESP cube file"""
    return ["5", "12", "2", "0", "q"]

def cube_laplacian():
    """Generate Laplacian cube file"""
    return ["5", "3", "2", "0", "q"]

def cube_fukui_minus():
    """Generate Fukui f- function cube"""
    return ["5", "0", "-1", "2", "0", "q"]

def cube_fukui_plus():
    """Generate Fukui f+ function cube"""
    return ["5", "0", "1", "2", "0", "q"]

def cube_dual_descriptor():
    """Generate dual descriptor cube"""
    return ["5", "0", "3", "2", "0", "q"]

def cube_orbital(orbital_index):
    """Generate specific orbital cube file"""
    return ["5", "0", "5", str(orbital_index), "2", "0", "q"]

# =========================================================================
# Main Menu 6: Check and modify wavefunction
# =========================================================================
def check_wavefunction():
    """Check wavefunction information"""
    return ["6", "q"]

def print_orbital_info():
    """Print orbital information"""
    return ["6", "1", "q"]

def print_basis_info():
    """Print basis function information"""
    return ["6", "2", "q"]

def modify_occupation():
    """Modify orbital occupation numbers"""
    return ["6", "3", "q"]

# =========================================================================
# Main Menu 7: Population analysis and atomic charges
# =========================================================================
def hirshfeld_charge():
    """Hirshfeld atomic charges"""
    return ["7", "1", "q"]

def vdd_population():
    """VDD (Voronoi Deformation Density) population"""
    return ["7", "2", "q"]

def mulliken_population():
    """Mulliken population analysis"""
    return ["7", "5", "q"]

def lowdin_population():
    """Löwdin population analysis"""
    return ["7", "6", "q"]

def becke_charge():
    """Becke atomic charges with dipole correction"""
    return ["7", "10", "q"]

def adch_charge():
    """ADCH (Atomic Dipole Corrected Hirshfeld) charges"""
    return ["7", "11", "1", "q"]

def chelpg_charge():
    """CHELPG ESP fitting charges"""
    return ["7", "12", "q"]

def mk_charge():
    """Merz-Kollman ESP fitting charges"""
    return ["7", "13", "q"]

def aim_charge():
    """AIM (QTAIM) atomic charges"""
    return ["7", "14", "q"]

def hirshfeld_i_charge():
    """Hirshfeld-I iterative charges"""
    return ["7", "15", "q"]

def cm5_charge():
    """CM5 charges"""
    return ["7", "16", "q"]

def eem_charge():
    """EEM (Electronegativity Equalization Method) charges"""
    return ["7", "17", "q"]

def resp_charge():
    """RESP (Restrained ESP) charges"""
    return ["7", "18", "q"]

def gasteiger_charge():
    """PEOE/Gasteiger charges"""
    return ["7", "19", "q"]

def mbis_charge():
    """MBIS (Minimal Basis Iterative Stockholder) charges"""
    return ["7", "20", "q"]

# =========================================================================
# Main Menu 8: Orbital composition analysis
# =========================================================================
def orbital_composition_mulliken(orbital_index):
    """Orbital composition by Mulliken method"""
    return ["8", "1", str(orbital_index), "q"]

def orbital_composition_nao(orbital_index):
    """Orbital composition by NAO method"""
    return ["8", "7", str(orbital_index), "q"]

def orbital_composition_hirshfeld(orbital_index):
    """Orbital composition by Hirshfeld method"""
    return ["8", "8", str(orbital_index), "q"]

def orbital_composition_becke(orbital_index):
    """Orbital composition by Becke method"""
    return ["8", "9", str(orbital_index), "q"]

def loba_oxidation_state():
    """LOBA/mLOBA oxidation state analysis"""
    return ["8", "100", "q"]

# =========================================================================
# Main Menu 9: Bond order analysis
# =========================================================================
def mayer_bond_order():
    """Mayer bond order analysis"""
    return ["9", "1", "q"]

def multicenter_bond_order():
    """Multi-center bond order (MCBO) analysis"""
    return ["9", "2", "q"]

def wiberg_bond_order():
    """Wiberg bond order in Löwdin basis"""
    return ["9", "3", "q"]

def mulliken_bond_order():
    """Mulliken bond order"""
    return ["9", "4", "q"]

def fuzzy_bond_order():
    """Fuzzy bond order"""
    return ["9", "7", "q"]

def laplacian_bond_order():
    """Laplacian bond order (LBO)"""
    return ["9", "8", "q"]

def wiberg_decomposition():
    """Wiberg bond order decomposition in NAO"""
    return ["9", "9", "q"]

def ibsi_analysis():
    """IBSI (Intrinsic Bond Strength Index)"""
    return ["9", "10", "q"]

def av1245_index():
    """AV1245 aromaticity index"""
    return ["9", "11", "q"]

# =========================================================================
# Main Menu 10: Plot DOS (Density of States)
# =========================================================================
def plot_dos():
    """Plot total DOS"""
    return ["10", "1", "q"]

def plot_pdos():
    """Plot partial DOS"""
    return ["10", "2", "q"]

def plot_opdos():
    """Plot overlap population DOS"""
    return ["10", "3", "q"]

def plot_ldos():
    """Plot local DOS"""
    return ["10", "4", "q"]

def plot_photoelectron_spectrum():
    """Plot photoelectron spectrum (PES)"""
    return ["10", "5", "q"]

def plot_cohp():
    """Plot COHP (Crystal Orbital Hamilton Population)"""
    return ["10", "6", "q"]

# =========================================================================
# Main Menu 11: Plot spectra (IR, Raman, UV-Vis, ECD, VCD, ROA, NMR)
# =========================================================================
def plot_ir_spectrum():
    """Plot IR spectrum"""
    return ["11", "1", "q"]

def plot_raman_spectrum():
    """Plot Raman spectrum"""
    return ["11", "2", "q"]

def plot_uv_vis_spectrum():
    """Plot UV-Vis absorption spectrum"""
    return ["11", "3", "q"]

def plot_ecd_spectrum():
    """Plot ECD (Electronic Circular Dichroism) spectrum"""
    return ["11", "4", "q"]

def plot_vcd_spectrum():
    """Plot VCD (Vibrational Circular Dichroism) spectrum"""
    return ["11", "5", "q"]

def plot_roa_spectrum():
    """Plot ROA (Raman Optical Activity) spectrum"""
    return ["11", "6", "q"]

def plot_nmr_spectrum():
    """Plot NMR spectrum"""
    return ["11", "7", "q"]

def plot_fluorescence_spectrum():
    """Plot fluorescence emission spectrum"""
    return ["11", "8", "q"]

def plot_pvs():
    """Plot partial vibrational spectrum (PVS)"""
    return ["11", "9", "q"]

def predict_color():
    """Predict color from UV-Vis spectrum"""
    return ["11", "10", "q"]

# =========================================================================
# Main Menu 12: Quantitative molecular surface analysis
# =========================================================================
def surface_analysis_esp():
    """ESP on molecular surface"""
    return ["12", "1", "q"]

def surface_analysis_alie():
    """Average Local Ionization Energy on surface"""
    return ["12", "2", "q"]

def surface_area_volume():
    """Calculate surface area and enclosed volume"""
    return ["12", "3", "q"]

def becke_surface():
    """Becke surface analysis"""
    return ["12", "4", "q"]

def hirshfeld_surface():
    """Hirshfeld surface analysis"""
    return ["12", "5", "q"]

def surface_extrema():
    """Locate surface minima and maxima"""
    return ["12", "6", "q"]

# =========================================================================
# Main Menu 13: Process grid data (cube files)
# =========================================================================
def export_cube():
    """Export current grid data to cube file"""
    return ["13", "0", "q"]

def grid_math_operations():
    """Mathematical operations on grid data"""
    return ["13", "11", "q"]

def grid_extract_plane():
    """Extract data in a plane"""
    return ["13", "2", "q"]

def grid_plot_integral_curve():
    """Plot integral curve along axis"""
    return ["13", "18", "q"]

# =========================================================================
# Main Menu 14: AdNDP (Adaptive Natural Density Partitioning)
# =========================================================================
def adndp_analysis():
    """Perform AdNDP analysis"""
    return ["14", "q"]

# =========================================================================
# Main Menu 15: Fuzzy atomic space analysis
# =========================================================================
def fuzzy_integrate_property():
    """Integrate property in fuzzy atomic spaces"""
    return ["15", "1", "q"]

def atomic_dipole_moments():
    """Calculate atomic dipole moments"""
    return ["15", "2", "q"]

def atomic_overlap_matrix():
    """Calculate atomic overlap matrix (AOM)"""
    return ["15", "3", "q"]

def localization_delocalization_index():
    """Calculate LI and DI indices"""
    return ["15", "4", "q"]

def pdi_aromaticity():
    """Para-delocalization index (PDI)"""
    return ["15", "5", "q"]

def flu_aromaticity():
    """Aromatic fluctuation index (FLU)"""
    return ["15", "6", "q"]

def flu_pi_aromaticity():
    """FLU-π aromaticity index"""
    return ["15", "7", "q"]

def multicenter_di():
    """Multi-center delocalization index"""
    return ["15", "11", "q"]

def ita_aromaticity():
    """Information-theoretic aromaticity index"""
    return ["15", "12", "q"]

def atomic_volume_polarizability():
    """Atomic effective volume and polarizability"""
    return ["15", "13", "q"]

# =========================================================================
# Main Menu 16: Charge Decomposition Analysis (CDA)
# =========================================================================
def cda_analysis():
    """Perform CDA analysis"""
    return ["16", "q"]

# =========================================================================
# Main Menu 17: Basin analysis
# =========================================================================
def basin_analysis_aim():
    """AIM basin analysis"""
    return ["17", "1", "q"]

def basin_analysis_elf():
    """ELF basin analysis"""
    return ["17", "2", "q"]

def basin_integrate_property():
    """Integrate properties in basins"""
    return ["17", "3", "q"]

# =========================================================================
# Main Menu 18: Electron excitation analysis
# =========================================================================
def hole_electron_analysis(state_num="1"):
    """Hole-electron distribution analysis"""
    return ["18", "1", str(state_num), "q"]

def transition_density_matrix():
    """Plot transition density matrix as heatmap"""
    return ["18", "2", "q"]

def charge_transfer_analysis():
    """Analyze charge transfer in excitation"""
    return ["18", "3", "q"]

def delta_r_index():
    """Calculate Δr charge-transfer index"""
    return ["18", "4", "q"]

def transition_dipole_moments():
    """Calculate transition dipole moments"""
    return ["18", "5", "q"]

def generate_nto():
    """Generate Natural Transition Orbitals"""
    return ["18", "6", "q"]

def ifct_analysis():
    """Interfragment charge transfer (IFCT) analysis"""
    return ["18", "8", "q"]

def lambda_index():
    """Calculate Λ index for excitation character"""
    return ["18", "14", "q"]

def cts_analysis():
    """Charge-transfer spectrum analysis"""
    return ["18", "16", "q"]

# =========================================================================
# Main Menu 19: Orbital localization
# =========================================================================
def boys_localization():
    """Boys orbital localization"""
    return ["19", "1", "q"]

def pipek_mezey_localization():
    """Pipek-Mezey orbital localization"""
    return ["19", "2", "q"]

# =========================================================================
# Main Menu 20: Weak interaction analysis
# =========================================================================
def nci_analysis():
    """NCI (Non-Covalent Interaction) analysis"""
    return ["20", "1", "q"]

def nci_promolecular():
    """NCI analysis with promolecular density"""
    return ["20", "2", "q"]

def anci_analysis():
    """Averaged NCI (aNCI) analysis"""
    return ["20", "3", "q"]

def iri_analysis():
    """IRI (Interaction Region Indicator) analysis"""
    return ["20", "4", "q"]

def dori_analysis():
    """DORI analysis"""
    return ["20", "5", "q"]

def vdw_potential():
    """van der Waals potential visualization"""
    return ["20", "6", "q"]

def igm_analysis():
    """IGM (Independent Gradient Model) analysis"""
    return ["20", "10", "q"]

def igmh_analysis():
    """IGMH analysis (IGM with Hirshfeld partition)"""
    return ["20", "11", "q"]

def aigm_analysis():
    """Averaged IGM (aIGM) analysis"""
    return ["20", "12", "q"]

def migm_analysis():
    """Modified IGM (mIGM) analysis"""
    return ["20", "-10", "q"]

# =========================================================================
# Main Menu 21: Energy Decomposition Analysis
# =========================================================================
def eda_ff():
    """EDA based on force field (EDA-FF)"""
    return ["21", "1", "q"]

def eda_sbl():
    """Shubin Liu's energy decomposition"""
    return ["21", "2", "q"]

def sobeda_analysis():
    """sobEDA analysis"""
    return ["21", "3", "q"]

def dispersion_atomic_contribution():
    """Atomic contribution to dispersion energy"""
    return ["21", "4", "q"]

# =========================================================================
# Main Menu 22: Conceptual DFT (CDFT)
# =========================================================================
def cdft_analysis():
    """Conceptual DFT analysis (Fukui, dual descriptor, etc.)"""
    return ["22", "q"]

def fukui_function():
    """Calculate Fukui functions"""
    return ["22", "1", "q"]

def dual_descriptor():
    """Calculate dual descriptor"""
    return ["22", "2", "q"]

def condensed_fukui():
    """Condensed Fukui function"""
    return ["22", "3", "q"]

# =========================================================================
# Main Menu 23: ETS-NOCV Analysis
# =========================================================================
def ets_nocv_analysis():
    """ETS-NOCV analysis"""
    return ["23", "q"]

# =========================================================================
# Main Menu 24: (Hyper)polarizability analysis
# =========================================================================
def parse_polarizability():
    """Parse Gaussian polar output"""
    return ["24", "1", "q"]

def sos_polarizability():
    """Sum-over-states polarizability calculation"""
    return ["24", "2", "q"]

def polarizability_density():
    """(Hyper)polarizability density"""
    return ["24", "3", "q"]

def unit_sphere_polarizability():
    """Unit sphere representation of polarizability"""
    return ["24", "5", "q"]

# =========================================================================
# Main Menu 25: Aromaticity analysis
# =========================================================================
def icss_analysis():
    """Iso-Chemical Shielding Surface (ICSS)"""
    return ["25", "3", "q"]

def nics_scan():
    """NICS scan (non-planar systems)"""
    return ["25", "4", "q"]

def homa_index():
    """HOMA aromaticity index"""
    return ["25", "6", "q"]

def homac_homer():
    """HOMAc and HOMER indices"""
    return ["25", "7", "q"]

def nics_1d_scan():
    """NICS-1D scan curve"""
    return ["25", "13", "q"]

def nics_2d_map():
    """NICS-2D plane map"""
    return ["25", "14", "q"]

# =========================================================================
# Main Menu 100: Other functions (Part 1)
# =========================================================================
def scatter_graph_two_functions():
    """Scatter graph between two functions"""
    return ["100", "1", "q"]

def export_various_files():
    """Export various file formats"""
    return ["100", "2", "q"]

def vdw_volume():
    """Calculate van der Waals volume"""
    return ["100", "3", "q"]

def integrate_whole_space():
    """Integrate function over whole space"""
    return ["100", "4", "q"]

def orbital_overlap_integral():
    """Overlap integral between alpha/beta orbitals"""
    return ["100", "5", "q"]

def monitor_scf_convergence():
    """Monitor Gaussian SCF convergence"""
    return ["100", "6", "q"]

def fragment_guess_input():
    """Generate fragment guess Gaussian input"""
    return ["100", "8", "q"]

def atomic_coordination():
    """Evaluate atomic coordination number"""
    return ["100", "9", "q"]

def orbital_overlap_centroid():
    """Orbital overlap and centroid distance"""
    return ["100", "11", "q"]

def biorthogonalization():
    """Biorthogonalization of alpha/beta orbitals"""
    return ["100", "12", "q"]

def lolipop_index():
    """LOLIPOP (LOL Integrated Pi Over Plane)"""
    return ["100", "14", "q"]

def intermolecular_overlap():
    """Intermolecular orbital overlap"""
    return ["100", "15", "q"]

def generate_fock_matrix():
    """Generate Fock/KS matrix"""
    return ["100", "17", "q"]

def electron_transport_route():
    """Yoshizawa's electron transport analysis"""
    return ["100", "18", "q"]

def combine_fragments():
    """Generate wavefunction from fragments"""
    return ["100", "19", "q"]

def hellmann_feynman_forces():
    """Calculate Hellmann-Feynman forces"""
    return ["100", "20", "q"]

def geometry_properties():
    """Geometry-based property calculations"""
    return ["100", "21", "q"]

def detect_pi_orbitals():
    """Detect π orbitals and composition"""
    return ["100", "22", "q"]

def fit_function_to_atoms():
    """Fit function distribution to atomic values"""
    return ["100", "23", "q"]

# =========================================================================
# Main Menu 200: Other functions (Part 2)
# =========================================================================
def cvb_index():
    """Core-Valence Bifurcation index"""
    return ["200", "1", "q"]

def atomic_bond_dipoles():
    """Atomic and bond dipole moments"""
    return ["200", "2", "q"]

def multiple_orbital_cubes():
    """Generate cube for multiple orbitals"""
    return ["200", "3", "q"]

def radial_distribution():
    """Radial distribution function"""
    return ["200", "5", "q"]

def orbital_correspondence():
    """Analyze orbital correspondence between wavefunctions"""
    return ["200", "6", "q"]

def average_bond_length():
    """Average bond length and coordination"""
    return ["200", "9", "q"]

def orbital_integrals():
    """Various orbital integrals"""
    return ["200", "10", "q"]

def function_moments():
    """Center, moments, radius of gyration"""
    return ["200", "11", "q"]

def energy_index():
    """Energy index (EI) and bond polarity index"""
    return ["200", "12", "q"]

def orbital_contributions_to_grid():
    """Orbital contributions to density difference"""
    return ["200", "13", "q"]

def domain_analysis():
    """Domain analysis within isosurfaces"""
    return ["200", "14", "q"]

def correlation_index():
    """Electron correlation index"""
    return ["200", "15", "q"]

def natural_orbitals():
    """Generate natural orbitals, NSO, SNO"""
    return ["200", "16", "q"]

def coulomb_exchange_integrals():
    """Coulomb and exchange integrals"""
    return ["200", "17", "q"]

def bla_boa_analysis():
    """Bond Length/Order Alternation"""
    return ["200", "18", "q"]

def spatial_delocalization_index():
    """Spatial delocalization index (SDI)"""
    return ["200", "19", "q"]

def bod_nado_analysis():
    """Bond Order Density and NAdO analysis"""
    return ["200", "20", "q"]

def lowdin_orthogonalization():
    """Löwdin orthogonalization of orbitals"""
    return ["200", "21", "q"]

# =========================================================================
# Main Menu 300: Other functions (Part 3)
# =========================================================================
def free_volume_in_cell():
    """View free regions and calculate free volume"""
    return ["300", "1", "q"]

def fit_atomic_radial_density():
    """Fit atomic radial density as STOs/GTFs"""
    return ["300", "2", "q"]

def stm_image():
    """Simulate STM image"""
    return ["300", "4", "q"]

def electric_multipole_moments():
    """Electric dipole/multipole moments"""
    return ["300", "5", "q"]

def orbital_energies_from_fock():
    """Calculate orbital energies from Fock matrix"""
    return ["300", "6", "q"]

def geometry_operations():
    """Geometry manipulation operations"""
    return ["300", "7", "q"]

def surface_distance_projection():
    """Surface distance projection map"""
    return ["300", "8", "q"]

def determine_fermi_level():
    """Determine Fermi level"""
    return ["300", "9", "q"]

# =========================================================================
# Composite/Advanced Functions
# =========================================================================
def complete_qtaim_analysis():
    """Complete QTAIM analysis workflow"""
    return ["2", "2", "-1", "3", "4", "q"]

def fukui_dual_descriptor_cubes():
    """Generate both Fukui and dual descriptor cubes"""
    seq = []
    # f- function
    seq.extend(["5", "0", "-1", "2", "fukui_minus.cub"])
    # f+ function  
    seq.extend(["5", "0", "1", "2", "fukui_plus.cub"])
    # dual descriptor
    seq.extend(["5", "0", "3", "2", "dual_descriptor.cub"])
    seq.append("q")
    return seq

def complete_aromaticity_analysis():
    """Comprehensive aromaticity analysis"""
    return ["25", "6", "13", "14", "q"]

def molecular_properties_summary():
    """Calculate multiple molecular properties"""
    seq = []
    seq.extend(["100", "3"])  # vdW volume
    seq.extend(["100", "21"])  # geometry properties
    seq.append("q")
    return seq

# =========================================================================
# Custom menu builder
# =========================================================================
def custom_menu(seq):
    """Provide arbitrary menu sequence for batch mode
    
    Args:
        seq: List of menu choices
        
    Returns:
        Menu sequence with proper exit
    """
    if isinstance(seq, str):
        seq = [seq]
    if not seq or seq[-1] != "q":
        seq.append("q")
    return seq

def batch_cube_generation(func_list):
    """Generate multiple cube files in one run
    
    Args:
        func_list: List of function codes to generate
        
    Example:
        batch_cube_generation([1, 4, 5, 12])  # density, ELF, LOL, ESP
    """
    seq = []
    for func in func_list:
        seq.extend(["5", str(func), "2", "0"])
    seq.append("q")
    return seq

# =========================================================================
# Helper function for common analyses
# =========================================================================
def weak_interaction_suite():
    """Run multiple weak interaction analyses"""
    analyses = [
        nci_analysis(),
        iri_analysis(),
        igm_analysis()
    ]
    # Combine all but remove redundant 'q's
    combined = []
    for analysis in analyses:
        combined.extend(analysis[:-1])
    combined.append("q")
    return combined

def charge_analysis_suite():
    """Run multiple charge analysis methods"""
    return {
        'hirshfeld': hirshfeld_charge(),
        'adch': adch_charge(),
        'chelpg': chelpg_charge(),
        'resp': resp_charge(),
        'cm5': cm5_charge()
    }

def bond_analysis_suite():
    """Run multiple bond order analyses"""
    return {
        'mayer': mayer_bond_order(),
        'wiberg': wiberg_bond_order(),
        'fuzzy': fuzzy_bond_order(),
        'laplacian': laplacian_bond_order()
    }

# =========================================================================
# Documentation helper
# =========================================================================
def list_all_functions():
    """Return dictionary of all available functions with descriptions"""
    import inspect
    
    functions = {}
    for name, obj in globals().items():
        if callable(obj) and not name.startswith('_'):
            doc = inspect.getdoc(obj)
            if doc:
                functions[name] = doc.split('\n')[0]
    
    return functions

def get_function_by_category():
    """Return functions organized by category"""
    categories = {
        'Visualization': ['view_structure', 'view_orbital'],
        'Topology': ['topology_search_cps', 'topology_generate_paths', 'topology_interbasin_surfaces'],
        'Population': ['hirshfeld_charge', 'mulliken_population', 'adch_charge', 'resp_charge'],
        'Orbital': ['orbital_composition_mulliken', 'orbital_composition_nao', 'boys_localization'],
        'Bond Analysis': ['mayer_bond_order', 'wiberg_bond_order', 'fuzzy_bond_order', 'ibsi_analysis'],
        'Spectroscopy': ['plot_ir_spectrum', 'plot_uv_vis_spectrum', 'plot_nmr_spectrum'],
        'Weak Interactions': ['nci_analysis', 'igm_analysis', 'iri_analysis', 'dori_analysis'],
        'Aromaticity': ['homa_index', 'nics_scan', 'flu_aromaticity', 'pdi_aromaticity'],
        'Excitation': ['hole_electron_analysis', 'generate_nto', 'ifct_analysis'],
        'Energy Decomposition': ['eda_ff', 'ets_nocv_analysis', 'cda_analysis'],
        'Surface Analysis': ['surface_analysis_esp', 'becke_surface', 'hirshfeld_surface'],
        'Grid/Cube': ['cube_density', 'cube_elf', 'cube_esp', 'grid_math_operations']
    }
    return categories