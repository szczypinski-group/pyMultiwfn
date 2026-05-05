"""Complete Multiwfn menu mapping for batch automation.

Notes
-----
Based on Multiwfn 3.8 (dev) manual.

"""

from enum import Enum


class Menu(Enum):
    """Enumeration of all Multiwfn menu functions.

    Sequences represent the keystrokes to navigate Multiwfn interactively.
    Each entry can be run standalone or composed sequentially via
    job.add_menu().

    Key conventions:
      - "0"  → return to main menu (from most submenus)
      - "-1" → trigger search-all or global option within a submenu
      - "q"  → quit (used in some contexts, e.g. batch mode)
      - "n"  → decline file-output prompt (skip writing .chg file)
      - "1"  → select built-in atomic densities (Hirshfeld-family methods)

    Examples
    --------
    >>> Menu.HIRSHFELD_CHARGE.get_sequence()
    ('7', '1', '1', 'n', '0')

    >>> job.add_menu(Menu.HIRSHFELD_CHARGE)
    >>> job.add_menu(Menu.MAYER_BOND_ORDER)
    """

    @classmethod
    def search(cls, query: str) -> list["Menu"]:
        """Search menu items by name (case-insensitive).

        Parameters
        ----------
        query : str
            Search string to match against member names

        Returns
        -------
        list[Menu]
            Matching menu items
        """
        query_upper = query.upper()
        return [item for item in cls if query_upper in item.name]

    @classmethod
    def list_all(cls) -> list[str]:
        """Return names of all menu items.

        Returns
        -------
        list[str]
            All member names
        """
        return [item.name for item in cls]

    def get_sequence(self) -> tuple[str, ...]:
        return self.value

    # TODO(fs): make sure all the menu sequences work as expected,
    # I will expand on testing each sequence
    # to make sure they work as expected.
    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 0: Show molecular structure / view orbitals
    # ─────────────────────────────────────────────────────────────────────────
    # Display 3D molecular structure and orbital isosurfaces in the GUI; prints
    # atom coordinates to screen
    VIEW_STRUCTURE = ("0",)

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 1: Output all properties at a point
    # Interactive only – prompts for a coordinate or atom index then prints
    # all supported real-space functions at that point (rho, ESP, ELF, ...).
    # INTERACTIVE ONLY - Ignore for now, output will not be parsed corectly or
    # output expected results
    # ─────────────────────────────────────────────────────────────────────────
    # Print every real-space function value (rho, ESP, ELF, G, K, ...) at a
    # user-supplied coordinate or nucleus
    # Requires xyz input and user interaction, ignore for now
    PROPERTIES_AT_POINT = ("1",)

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 2: Topology analysis
    # Uses Newton iteration to locate critical points (CPs) of a chosen
    # real-space function, then traces gradient paths between them.
    # CP types: (3,-3) nuclear/max, (3,-1) bond, (3,+1) ring, (3,+3) cage.
    # ─────────────────────────────────────────────────────────────────────────
    # NEW SEQUENCES TO BE TESTED
    TOPOLOGY_VISUALISE_CPS = ("2", "0")
    TOPOLOGY_CP_NUCLEAR_POSITION = ("2", "2")
    TOPOLOGY_CP_MIDPOINTS = ("2", "3")
    TOPOLOGY_CP_TRIANGLE_CENTRES = ("2", "4")
    TOPOLOGY_CP_PYRAMID_CENTRES = ("2", "5")
    TOPOLOGY_CP_SPHERE_POINTS = ("2", "6", "0")
    # Outputs CPout.txt
    TOPOLOGY_CP_REAL_SPACE_POINTS = ("2", "7", "0")
    TOPOLOGY_CP_PATHS_3MINUS3_3MINUS1 = ("2", "8")
    TOPOLOGY_CP_PATHS_3PLUS1_3PLUS3 = ("2", "9")

    # ALL FOLLOWING SEQUENCES ARE DEPRECATED
    # Search ALL critical points of rho starting from every nuclear position;
    # finds NCPs, BCPs, RCPs, and CCPs in one pass
    TOPOLOGY_SEARCH_CPS = ("2", "2", "-1")
    # Full AIM workflow in one sequence: search all CPs, generate bond paths,
    # then generate interbasin surfaces
    TOPOLOGY_ANALYSIS_COMPLETE = ("2", "2", "-1", "3", "4")
    # Topology analysis of the electrostatic potential (ESP): locate ESP
    # critical points starting from all nuclei
    TOPOLOGY_ESP_ANALYSIS = ("2", "-2", "2", "-1")
    # Topology analysis of the Localized Orbital Locator (LOL): find LOL maxima
    # corresponding to bonding and lone-pair basins
    TOPOLOGY_LOL_ANALYSIS = ("2", "-10", "2", "-1")
    # Topology analysis of the Electron Localization Function (ELF): find ELF
    # attractors and characterise bonding basins
    TOPOLOGY_ELF_ANALYSIS = ("2", "9", "2", "-1")
    # Topology analysis of nabla^2 rho: locate charge-concentration (3,-3) and
    # charge-depletion (3,+3) critical points
    TOPOLOGY_LAPLACIAN_ANALYSIS = ("2", "3", "2", "-1")
    # Search for bond critical points (BCPs) using atom-pair midpoints as
    # Newton starting guesses
    TOPOLOGY_SEARCH_BCP = ("2", "3", "-1")
    # Search for ring critical points (RCPs) using triangle centres of atom
    # triplets as starting guesses
    TOPOLOGY_SEARCH_RCP = ("2", "4", "-1")
    # Search for cage critical points (CCPs) using pyramid centres of atom
    # quartets as starting guesses
    TOPOLOGY_SEARCH_CCP = ("2", "5", "-1")

    # ────────────────────────────────────────────────────────────────────────
    # Main Menu 3: Output / plot property along a line
    # Evaluates a real-space function at 3000 points between two atoms or
    # coordinates and produces a curve map (line.txt exported on request).
    # FULLY INTERACTIVE - Ignore for now, output will not be parsed corectly
    # or output expected results
    # ────────────────────────────────────────────────────────────────────────
    # Plot total electrostatic potential (ESP) along a line; reveals
    # electrophilic and nucleophilic interaction sites between atoms
    LINE_ESP = ("3", "12")
    # Plot electron density rho(r) along a line; shows charge accumulation and
    # depletion between bonded atoms
    LINE_ELECTRON_DENSITY = ("3", "1")
    # Plot Laplacian of electron density nabla^2 rho along a line; negative
    # regions indicate charge concentration (covalent bonds, lone pairs)
    LINE_LAPLACIAN = ("3", "3")
    # Plot Electron Localization Function (ELF, range 0-1) along a line; peaks
    # identify bonding pairs and lone pairs
    LINE_ELF = ("3", "9")
    # Plot Localized Orbital Locator (LOL, range 0-1) along a line; similar to
    # ELF but with sharper basin boundary definition
    LINE_LOL = ("3", "10")
    # Plot Reduced Density Gradient (RDG) along a line; low-RDG regions reveal
    # non-covalent interaction zones
    LINE_RDG = ("3", "13")
    # Plot spin density rho_alpha - rho_beta along a line; used for radical and
    # open-shell systems to locate unpaired electrons
    LINE_SPIN_DENSITY = ("3", "5")
    # Plot magnitude of the electron density gradient |nabla rho| along a line;
    # maxima mark interatomic surface boundaries
    LINE_GRADIENT_NORM = ("3", "2")
    # Plot Lagrangian (positive-definite) kinetic energy density G(r) along a
    # line
    LINE_KINETIC_G = ("3", "7")
    # Plot Hamiltonian kinetic energy density K(r) along a line; related to G:
    # K = G - (1/4) nabla^2 rho
    LINE_KINETIC_K = ("3", "6")
    # Plot Average Local Ionization Energy (ALIE) along a line; low values
    # indicate weakly bound, reactive electrons susceptible to electrophilic
    # attack
    LINE_ALIE = ("3", "18")
    # Plot Source Function SF(r', r) along a line relative to a fixed reference
    # point r set in settings.ini; shows electron-density contributions from
    # each spatial region
    LINE_SOURCE_FUNCTION = ("3", "19")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 4: Output / plot property in a plane
    # The third tuple element encodes the plane-definition mode:
    #   "1" = XY plane (input Z value)
    #   "2" = XZ plane (input Y value)
    #   "3" = YZ plane (input X value)
    #   "4" = plane defined by three atom indices
    #   "5" = plane defined by three Cartesian points
    # Produces colour-filled, contour, relief, or gradient/vector maps.
    # FULLY INTERACTIVE - Ignore for now, output will not be parsed corectly
    # or output expected results
    # ─────────────────────────────────────────────────────────────────────────
    # Colour-filled map of electron density rho(r) in the XY plane; visualises
    # how charge is distributed across the molecule
    PLANE_MAP_DENSITY = ("4", "1", "1")
    # Colour-filled map of total ESP in the XY plane; red regions are
    # nucleophilic (-), blue regions are electrophilic (+)
    PLANE_MAP_ESP = ("4", "12", "1")
    # Colour-filled map of ELF in the XY plane; red/orange islands mark bonding
    # pairs and lone pairs, blue marks depleted regions
    PLANE_MAP_ELF = ("4", "9", "1")
    # Colour-filled map of LOL in the XY plane; chemically equivalent to ELF
    # with cleaner inter-basin boundaries
    PLANE_MAP_LOL = ("4", "10", "1")
    # Colour-filled map of |nabla rho| in the XY plane; highlights interatomic
    # surface and shell-structure regions
    PLANE_MAP_GRADIENT = ("4", "2", "1")
    # Colour-filled or contour map of nabla^2 rho in the XY plane; negative =
    # charge concentration (bonds/lone pairs), positive = depletion
    PLANE_MAP_LAPLACIAN = ("4", "3", "1")
    # Colour-filled map of spin density rho_alpha - rho_beta in the XY plane;
    # shows location of unpaired electrons in open-shell species
    PLANE_MAP_SPIN_DENSITY = ("4", "5", "1")
    # Colour-filled map of RDG in the XY plane; low-value regions reveal non-
    # covalent interactions when combined with sign(lambda2)rho colouring
    PLANE_MAP_RDG = ("4", "13", "1")
    # Colour-filled map of sign(lambda2)*rho in the XY plane; negative =
    # attractive NCI (H-bond), positive = repulsive steric clash
    PLANE_MAP_SIGN_LAMBDA2_RHO = ("4", "15", "1")
    # Colour-filled map of ALIE in the XY plane;
    # low-ALIE pockets on the surface
    # predict electrophilic and radical attack sites
    PLANE_MAP_ALIE = ("4", "18", "1")
    # Colour-filled map of Lagrangian kinetic energy density G(r) in the XY
    # plane
    PLANE_MAP_KINETIC_G = ("4", "7", "1")
    # Colour-filled map of Hamiltonian kinetic energy density K(r) in the XY
    # plane
    PLANE_MAP_KINETIC_K = ("4", "6", "1")
    # Colour-filled map of Source Function in the XY plane relative to the
    # reference point defined in settings.ini
    PLANE_MAP_SOURCE_FUNCTION = ("4", "19", "1")
    # Colour-filled or contour map of a single MO wavefunction psi_i in the XY
    # plane; user is prompted to supply the orbital index
    PLANE_MAP_ORBITAL_WAVEFUNCTION = ("4", "4", "1")
    # Colour-filled map of Fukui f-(r) in the XY plane via custom-operation
    # rho_N - rho_{N-1}; predicts electrophilic attack sites
    PLANE_MAP_FUKUI_MINUS = ("4", "0", "-1")
    # Colour-filled map of Fukui f+(r) in the XY plane via rho_{N+1} - rho_N;
    # predicts nucleophilic attack sites
    PLANE_MAP_FUKUI_PLUS = ("4", "0", "1")
    # Colour-filled map of dual descriptor Delta_f = f+ - f- in the XY plane;
    # positive = nucleophilic centre, negative = electrophilic centre
    PLANE_MAP_DUAL_DESCRIPTOR = ("4", "0", "3")
    # Colour-filled map of deformation density rho_mol - rho_promol in the XY
    # plane; shows electron redistribution upon chemical bond formation
    PLANE_MAP_DEFORMATION_DENSITY = ("4", "-2")
    # Colour-filled map of promolecular density (superposition of free-atom
    # densities) in the XY plane; reference state before bonding
    PLANE_MAP_PROMOLECULAR_DENSITY = ("4", "-1")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 5: Cube / grid generation
    # Evaluates a real-space function on a 3D grid and exports a Gaussian
    # .cube file compatible with VMD, GaussView, ChemCraft, Molekel, etc.
    # Grid quality codes:
    #   "1" ~ 50^3 points (low / preview)
    #   "2" ~ 80^3 points (medium, default for most workflows)
    #   "3" ~ 120^3 points (high, for publication-quality figures)
    # "0" at the post-process prompt triggers cube export and returns.
    # ─────────────────────────────────────────────────────────────────────────
    # NEW SEQUENCES - TO BE TESTED
    CUBE_LOCAL_INFORMATION_ENTROPY = ("5", "11", "2", "0")
    CUBE_ELECTROSTATIC_POTENTIAL_FROM_CHARGE = ("5", "8", "2", "0")
    CUBE_RDG_QUICK = ("5", "14", "2", "0")
    CUBE_CORRELATION_HOLE_ALPHA = ("5", "17", "2", "1", "2", "0")
    CUBE_EDR = ("5", "20", "2", "2", "2", "0")
    CUBE_DELTAG_PROMOLECULAR_APPROX = ("5", "22", "2", "2", "0")
    CUBE_DELTAG_HIRSHFELD_PAER_APPROX = ("5", "23", "2", "2", "0")
    CUBE_IRI = ("5", "24", "2", "2", "0")
    CUBE_VDW_POTENTIAL = ("5", "25", "2", "2", "0")

    # EDITED SEQUENCES - TO BE TESTED
    CUBE_ALIE = ("5", "18", "2", "2", "0")
    CUBE_SOURCE_FUNCTION = ("5", "19", "2", "2", "0")

    # CORRECT SEQUENCES - TO BE TESTED
    # Generate medium-quality .cube of electron density rho(r); standard input
    # for VMD/GaussView isosurface rendering
    CUBE_DENSITY = ("5", "1", "2", "2", "0")
    # Generate medium-quality .cube of spin density rho_alpha - rho_beta;
    # visualise unpaired electrons in radicals and open-shell molecules
    CUBE_SPIN_DENSITY = ("5", "5", "2", "2", "0")
    # Generate medium-quality .cube of ELF; bonding pairs, lone pairs, and core
    # shells visible as isosurfaces
    CUBE_ELF = ("5", "9", "2", "2", "0")
    # Generate medium-quality .cube of LOL; electron-pair localisation function
    # with sharper basin boundaries than ELF
    CUBE_LOL = ("5", "10", "2", "2", "0")
    # Generate medium-quality .cube of total ESP; map onto rho=0.001 isosurface
    # to produce a surface electrostatic potential map
    CUBE_ESP = ("5", "12", "2", "2", "0")
    # Generate medium-quality .cube of nabla^2 rho; negative isosurfaces encode
    # charge concentration (covalent bonds, lone pairs, shells)
    CUBE_LAPLACIAN = ("5", "3", "2", "2", "0")
    # Generate medium-quality .cube of |nabla rho|; high-value isosurfaces mark
    # interatomic and atomic-shell boundaries
    CUBE_GRADIENT_NORM = ("5", "2", "2", "2", "0")
    # Generate medium-quality .cube of G(r), the positive-definite Lagrangian
    # kinetic energy density
    CUBE_KINETIC_G = ("5", "7", "2", "2", "0")
    # Generate medium-quality .cube of K(r), the Hamiltonian kinetic energy
    # density; K = G - (1/4) nabla^2 rho
    CUBE_KINETIC_K = ("5", "6", "2", "2", "0")
    # Generate medium-quality .cube of RDG; render at low isovalue (~0.5)
    # together with sign(lambda2)rho cube for NCI isosurface colouring
    CUBE_RDG = ("5", "13", "2", "2", "0")
    # Generate medium-quality .cube of sign(lambda2)*rho; colour-maps onto RDG
    # isosurface: blue=H-bond, green=vdW, red=steric repulsion
    CUBE_SIGN_LAMBDA2_RHO = ("5", "15", "2", "2", "0")

    # INTERACTIVE OR IMPLICITLY INTERACTIVE SEQUENCES - TO BE TESTED
    # Generate medium-quality .cube of a single MO wavefunction psi_i; user is
    # prompted for the orbital index; standard for orbital visualisation
    # INTERACTIVE
    CUBE_ORBITAL_WAVEFUNCTION = ("5", "4", "2", "2", "0")
    # Generate medium-quality .cube of Fukui f-(r) via two-wavefunction
    # subtraction rho_N - rho_{N-1}; maps electrophilic attack susceptibility
    # NOT POSSIBLE - REQUIRES MULTIPLE WAVEFUNCTIONS
    CUBE_FUKUI_MINUS = ("5", "0", "-1", "2", "0")
    # Generate medium-quality .cube of Fukui f+(r) via rho_{N+1} - rho_N; maps
    # nucleophilic attack susceptibility
    CUBE_FUKUI_PLUS = ("5", "0", "1", "2", "0")
    # Generate medium-quality .cube of dual descriptor f+ - f-; positive
    # isosurfaces are nucleophilic sites, negative are electrophilic
    CUBE_DUAL_DESCRIPTOR = ("5", "0", "3", "2", "0")
    # Generate medium-quality .cube of promolecular density (superposition of
    # free-atom densities); used as a reference before bond formation
    CUBE_PROMOLECULAR_DENSITY = ("5", "-1", "1", "2", "0")
    # Generate medium-quality .cube of deformation density Delta_rho = rho_mol
    #     # rho_promol; shows how bonding redistributes electron density
    CUBE_DEFORMATION_DENSITY = ("5", "-2", "1", "2", "0")
    # Generate HIGH-quality .cube of electron density (~120^3 grid); finer
    # isosurfaces suitable for publication figures
    CUBE_DENSITY_HIGH = ("5", "1", "3", "0")
    # Generate HIGH-quality .cube of total ESP (~120^3 grid); needed for
    # accurate surface-ESP colour maps on large or complex molecules
    CUBE_ESP_HIGH = ("5", "12", "3", "0")
    # Generate HIGH-quality .cube of ELF (~120^3 grid); resolves fine bonding
    # features in systems requiring detailed electron-pair analysis
    CUBE_ELF_HIGH = ("5", "9", "3", "0")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 6: Check & modify wavefunction
    # Provides subfunctions to inspect, edit, and save the loaded wavefunction.
    # "0" from the submenu returns to main menu.
    # ─────────────────────────────────────────────────────────────────────────
    # NEW SEQUENCES - TO BE TESTED
    PRINT_COEFFICIENT_MATRIX = ("6", "5", "1")

    PRINT_DENSITY_MATRIX = ("6", "6", "1")

    PRINT_INTEGRAL_MATRIX_FOCK = ("6", "7", "0", "1", "1")
    PRINT_INTEGRAL_MATRIX_OVERLAP = ("6", "7", "1", "1")
    PRINT_INTEGRAL_MATRIX_ELECTRIC_DIPOLE = ("6", "7", "2", "1")
    PRINT_INTEGRAL_MATRIX_MAGNETIC_DIPOLE = ("6", "7", "3", "1")
    PRINT_INTEGRAL_MATRIX_VELOCITY = ("6", "7", "4", "1")
    PRINT_INTEGRAL_MATRIX_EKINETIC = ("6", "7", "5", "1")
    PRINT_INTEGRAL_MATRIX_QUADRUPOLE = ("6", "7", "6", "1")
    PRINT_INTEGRAL_MATRIX_OCTOPOLE = ("6", "7", "7", "1")
    PRINT_INTEGRAL_MATRIX_HEXADECAPOLE = ("6", "7", "8", "1")

    # #INTERACTIVE SEQUENCES - to be tested
    # #CORRECT SEQUENCES
    # # Save the (possibly modified) wavefunction to new.wfn; also converts
    # # fch/molden to .wfn format; zero-occupation orbitals are dropped
    # # automatically
    # SAVE_WFN = ("6", "0")
    # Print centre atom, angular-momentum type (s/p/d/f/...), and exponent for
    # # every Gaussian-type function (GTF) in the basis
    # PRINT_ALL_GTF = ("6", "1")
    # Print shell assignments, contracted function types, and GTF index ranges
    # # for every basis function
    # PRINT_ALL_BASIS_FUNCTIONS = ("6", "2")
    # # Print index, energy, occupation number, and spin type for every o
    # # rbital in the loaded wavefunction
    # PRINT_ORBITAL_INFO = ("6", "3")

    # # Print the one-particle density matrix P expressed in the basis-function
    # # representation
    # PRINT_DENSITY_MATRIX = ("6", "6", "0")
    # # Manually set the occupation number of selected orbitals; set to 0 to
    # remove their contribution from subsequent real-space function evaluations
    # MODIFY_OCCUPATION = ("6", "26", "0")
    # # Remove all core (inner-shell) orbitals from the wavefunction, retaining
    # # only valence-shell orbitals for subsequent analyses
    # DELETE_INNER_ORBITALS = ("6", "34", "0")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 7: Population analysis & atomic charges
    # Hirshfeld-family methods require atomic reference densities:
    #   "1" → use Multiwfn's built-in sphericalised free-atom densities
    #   "n" → skip writing the .chg output file to disk
    #   "0" → return to main menu after printing charges
    # ─────────────────────────────────────────────────────────────────────────
    # NEW SEQUENCES - TO BE TESTED
    MULLIKEN_DECOMPOSE_ATOMIC_POPULATION = ("7", "5", "2", "n")
    MULLIKEN_DECOMPOSE_BASIS_FUNCTION = ("7", "5", "3", "n")

    # INTERACTIVE SEQUENCES - TO BE TESTED
    # Lowdin atomic charges via Lowdin orthogonalisation; slightly more basis-
    # set stable than Mulliken but still dependent on basis choice
    # INTERACTIVE
    LOWDIN_POPULATION = ("7", "6", "ENTER", "0")

    # INGNORE FOR NOW
    # Electronegativity Equalization Method (EEM) charges; fast geometry-based
    # empirical model requiring no wavefunction evaluation
    # NEEDS MOL2
    EEM_CHARGE = ("7", "17", "n", "0")

    # CORRECT SEQUENCES - TO BE TESTED
    # Hirshfeld atomic charges from deformation-density partitioning;
    # qualitatively correct but systematically underestimates charge transfer
    HIRSHFELD_CHARGE = ("7", "1", "1", "n", "0")
    # Voronoi Deformation Density charges; Hirshfeld variant using Voronoi
    # cell weights instead of Hirshfeld weights; results similar to Hirshfeld
    VDD_POPULATION = ("7", "2", "1", "n", "0")
    # Mulliken atomic charges and basis-function populations; oldest method,
    # highly basis-set dependent, avoid with diffuse basis functions
    MULLIKEN_POPULATION = ("7", "5", "1", "n", "0")

    # Ros-Schuit C-squared Population Analysis (SCPA); modified Mulliken scheme
    # that prevents negative population numbers
    SCPA_POPULATION = ("7", "7", "n", "0")
    # Stout-Politzer modified Mulliken charges; cross-terms partitioned by the
    # ratio of squared orbital coefficients
    STOUT_POLITZER_POPULATION = ("7", "8", "n", "0")
    # Bickelhaupt modified Mulliken charges; cross-terms weighted by the total
    # local populations summed across all orbitals
    BICKELHAUPT_POPULATION = ("7", "9", "n", "0")
    # Becke atomic charges with atomic-dipole-moment correction applied;
    # reasonable for typical organic systems using default CSD radii
    BECKE_CHARGE = ("7", "10", "0", "n", "0")
    # Atomic Dipole moment Corrected Hirshfeld (ADCH) charges; exactly
    # reproduces the molecular dipole moment and gives reliable ESP; highly
    # recommended
    ADCH_CHARGE = ("7", "11", "1", "n", "0")
    # CHELPG ESP-fitting charges on a cubic grid of points; best rotational
    # invariance among ESP-fit methods; widely used for force-field
    # parametrisation
    CHELPG_CHARGE = ("7", "12", "1", "n", "0")
    # Merz-Kollmann (MK) ESP-fitting charges on concentric shells at 1.4-2.0x
    # vdW radii; standard for AMBER/GAFF parametrisation
    MK_CHARGE = ("7", "13", "1", "n", "0")

    # CM5 charges (charge model 5); Hirshfeld-based with empirical correction
    # terms for improved dipole-moment reproduction across diverse molecules
    CM5_CHARGE = ("7", "16", "1", "n", "0")

    # RESP (Restrained ESP) charges; ESP fit with hyperbolic
    # restraint toward zero; the standard charge model for AMBER/GAFF force
    # fields
    # FORCED DEFAULTS
    RESP_CHARGE = ("7", "18", "1n", "0")
    # Gasteiger-Marsili empirical charges; purely connectivity-based, extremely
    # fast, no wavefunction needed; suitable for large databases
    GASTEIGER_CHARGE = ("7", "19", "n", "0")
    # Minimal Basis Iterative Stockholder (MBIS) charges; information-theoretic
    # partitioning giving excellent dipole and higher-multipole reproduction
    MBIS_CHARGE = ("7", "20", "1", "n", "0")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 8: Orbital composition analysis
    # Decomposes each MO into percentage contributions from basis functions,
    # shells, atoms, or user-defined fragments.
    # ─────────────────────────────────────────────────────────────────────────
    # NEW SEQUENCES - TO BE TESTED
    ORBITAL_COMPOSITION_MULLIKEN_HOMO = ("8", "1", "h")
    ORBITAL_COMPOSITION_MULLIKEN_LUMO = ("8", "1", "l")
    ORBITAL_COMPOSITION_MULLIKEN_ALL = ("8", "1", "a")

    ORBITAL_COMPOSITION_STOUT_POLITZER_HOMO = ("8", "2", "h")
    ORBITAL_COMPOSITION_STOUT_POLITZER_LUMO = ("8", "2", "l")
    ORBITAL_COMPOSITION_STOUT_POLITZER_ALL = ("8", "2", "a")

    ORBITAL_COMPOSITION_SCPA_HOMO = ("8", "3", "h")
    ORBITAL_COMPOSITION_SCPA_LUMO = ("8", "3", "l")
    ORBITAL_COMPOSITION_SCPA_ALL = ("8", "3", "a")

    FRAGMENT_CONTRIBUTION_HIRSHFELD_HOMO = ("8", "8", "1", "h")
    FRAGMENT_CONTRIBUTION_HIRSHFELD_LUMO = ("8", "8", "1", "l")
    FRAGMENT_CONTRIBUTION_HIRSHFELD_ALL = ("8", "8", "1", "-1")
    ATOM_CONTRIBUTION_HIRSHFELD = ("8", "8", "1", "-4")

    FRAGMENT_CONTRIBUTION_BECKE_HOMO = ("8", "9", "1", "h")
    FRAGMENT_CONTRIBUTION_BECKE_LUMO = ("8", "9", "1", "l")
    FRAGMENT_CONTRIBUTION_BECKE_ALL = ("8", "9", "1", "-1")
    ATOM_CONTRIBUTION_BECKE = ("8", "9", "1", "-4")

    # INTERACTIVE SEQUENCES - IGNORE
    # Mulliken fragment composition: prints the total percentage of a pre-
    # defined fragment in every occupied MO as a table
    ORBITAL_COMPOSITION_FRAGMENT_MULLIKEN = ("8", "4", "0")
    # Stout-Politzer fragment composition including inter-fragment cross-term
    # breakdown
    ORBITAL_COMPOSITION_FRAGMENT_STOUT = ("8", "5", "0")
    # SCPA fragment composition: sums C^2 contributions of all fragment basis
    # functions per orbital; no negative values
    ORBITAL_COMPOSITION_FRAGMENT_SCPA = ("8", "6", "0")
    # Natural Atomic Orbital (NAO) composition using the MO-in-NAO coefficient
    # matrix from NBO output; excellent basis-set stability for occupied MOs
    # REQUIRES SPECIFIC FILE TYPE
    ORBITAL_COMPOSITION_NAO = ("8", "7", "0")

    # Hirshfeld orbital composition: integral of |psi_i|^2 * w_A(r) for each
    # atom A; highly stable regardless of basis-set choice
    ORBITAL_COMPOSITION_HIRSHFELD = ("8", "8", "0")
    # Becke orbital composition: Becke-partition integral of |psi_i|^2 per atom
    # no atomic reference density files required
    ORBITAL_COMPOSITION_BECKE = ("8", "9", "0")
    # Modified LOBA (Localized Orbital Bonding Analysis): assigns formal
    # oxidation states by analysing orbital-population partitioning
    # among bonded atoms
    LOBA_OXIDATION_STATE = ("8", "100", "0")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 9: Bond order analysis
    # Prints bond orders between all atom pairs above a threshold, plus
    # total and free valences for each atom.  "0" returns to main menu.
    # ─────────────────────────────────────────────────────────────────────────
    # NEW SEQUENCES - TO BE TESTED

    # INTERACTIVE - REQUIRES USER INPUT
    MULTICENTER_BOND_ORDER = ("9", "2", "0")
    # Multi-centre bond order in the NAO basis (S becomes identity); much more
    # basis-set stable, recommended when diffuse functions are present
    MULTICENTER_BOND_ORDER_NAO = ("9", "-2", "0")
    # Wiberg bond order in the Lowdin-orthogonalised basis (WL); more stable
    # than Mayer for large basis sets but can overestimate polar bonds
    # Decompose Mulliken bond order between a chosen atom pair into per-orbital
    # contributions to identify bonding vs. antibonding MOs
    MULLIKEN_BOND_ORDER_DECOMPOSE = ("9", "5", "0")
    # Orbital-occupancy-perturbed Mayer bond order: prints the contribution of
    # each occupied MO to the Mayer bond order for a selected A-B pair
    ORBITAL_PERTURBED_MAYER = ("9", "6", "0")
    # Decompose Wiberg bond order between two atoms into per-orbital
    # contributions; identifies which MOs are responsible for the bond
    # NEEDS NAO BASIS
    WIBERG_DECOMPOSITION = ("9", "9", "0")
    # AV1245 multicentric aromaticity index computed from Mayer bond orders
    # along a ring; robust against ring size and basis-set choice
    AV1245_INDEX = ("9", "11", "0")

    # CORRECT SEQUENCES - TO BE TESTED
    # Mayer bond orders from PS-matrix products; values approximately 1/2/3 for
    # single/double/triple bonds; the best all-round bond-order method
    MAYER_BOND_ORDER = ("9", "1", "n0")
    # Multi-centre bond order (up to 12 centres) in the original basis-function
    # representation; sensitive to diffuse basis functions
    WIBERG_BOND_ORDER = ("9", "3", "n", "0")
    # Mulliken bond orders (2*PS off-diagonal elements); positive = bonding
    # character, negative = antibonding; qualitative indicator only
    MULLIKEN_BOND_ORDER = ("9", "4", "n", "0")
    # Fuzzy bond order (Becke-space integration of PS products); more basis-set
    # stable than Mayer; essentially equivalent to the AIM delocalization index
    FUZZY_BOND_ORDER = ("9", "7", "n", "0")
    # Laplacian Bond Order (LBO): integral of -nabla^2 rho in the fuzzy overlap
    # space; correlates with BDE and vibrational frequency; independent of
    # wavefunction type
    LAPLACIAN_BOND_ORDER = ("9", "8", "n", "0")
    # Intrinsic Bond Strength Index (IBSI): geometry-free bond-order measure
    # derived from the electron density at the bond critical point
    IBSI_ANALYSIS = ("9", "10", "1", "1", "0")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 10: Density of states
    # Reads orbital energies from .fch/.molden/Gaussian-output/plain-text
    # and produces broadened TDOS, PDOS (per fragment), and OPDOS curves.
    # ─────────────────────────────────────────────────────────────────────────
    # Plot Total Density of States (TDOS): broadened orbital energy spectrum
    # showing how densely electronic states are distributed in energy
    PLOT_TDOS = ("10", "0", "2")

    PLOT_TDOS_OPDOS = ("10", "00", "2")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 11: Spectra simulation
    # Reads transition data from Gaussian/ORCA output or plain-text files
    # and broadens discrete transitions into simulated spectra using
    # Gaussian, Lorentzian, or pseudo-Voigt broadening functions.
    # REQUIRES .OUT FILES - IGNORE FOR NOW
    # ─────────────────────────────────────────────────────────────────────────
    # Simulate IR absorption spectrum by broadening harmonic (or anharmonic)
    # vibrational frequencies weighted by IR intensities in km/mol
    PLOT_IR_SPECTRUM = ("11", "1", "0", "1")
    # Simulate Raman spectrum from Raman activities; optionally converts
    # activities to intensities given a laser wavelength and temperature
    PLOT_RAMAN_SPECTRUM = ("11", "2", "0", "1")
    # Simulate UV-Vis absorption spectrum by broadening TD-DFT/CIS excitation
    # energies weighted by oscillator strengths; area calibrated to molar
    # absorptivity
    PLOT_UV_VIS_SPECTRUM = ("11", "3", "0", "1")
    # Simulate Electronic Circular Dichroism (ECD) spectrum from rotatory
    # strengths in either length or velocity gauge representation
    PLOT_ECD_SPECTRUM = ("11", "4", "0", "1")
    # Simulate Vibrational Circular Dichroism (VCD) spectrum by broadening
    # vibrational rotatory strengths
    PLOT_VCD_SPECTRUM = ("11", "5", "0", "1")
    # Simulate Raman Optical Activity (ROA) spectrum from computed ROA
    # intensities
    PLOT_ROA_SPECTRUM = ("11", "6", "0", "1")
    # Simulate NMR spectrum from calculated chemical shielding tensors; convert
    # to chemical shifts by subtracting a reference shielding value
    PLOT_NMR_SPECTRUM = ("11", "7", "0", "1")
    # Simulate fluorescence/phosphorescence spectrum; applies Kasha rule so
    # only the first excited state contributes when appropriate
    PLOT_FLUORESCENCE_SPECTRUM = ("11", "8", "0", "1")
    # Plot photoelectron valence spectrum  by broadening orbital ionisation
    # energies; analogous to TDOS but focused on valence region
    PLOT_PVS = ("11", "9", "0", "1")
    # Predict the perceived colour of a compound from its computed UV-Vis
    # absorption spectrum using CIE colour matching functions
    PREDICT_COLOR = ("11", "0")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 12: Quantitative molecular surface analysis
    # Generates the rho = 0.001 a.u. vdW isosurface via Marching Tetrahedra,
    # maps a chosen function onto it, and computes surface statistical
    # descriptors (V_S, sigma^2, Pi, etc.) and locates surface extrema.
    # ─────────────────────────────────────────────────────────────────────────
    # NEW SEQUENCES - TO BE TESTED
    QMSA_ESP = ("12", "0", "-2")
    QMSA_ALIE = ("12", "2", "2", "0", "-2")
    QMSA_LEA = ("12", "2", "4", "0", "-2")
    QMSA_LEAE = ("12", "2", "-4", "0", "-2")
    QMSA_EDR = ("12", "2", "5", "0", "-2")
    QMSA_MAXEDR = ("12", "2", "6", "0", "-2")
    QMSA_EDENSITY = ("12", "2", "11", "0", "-2")
    QMSA_LAMBDA2_RHO = ("12", "2", "12", "0", "-2")

    # DEPRECATED SEQUENCES - TO BE TESTED
    # Map ESP onto the vdW surface (rho=0.001); compute V_S+, V_S-, sigma^2, Pi
    # and locate surface ESP minima/maxima; outputs GIPF descriptors for
    # property prediction
    SURFACE_ANALYSIS_ESP = ("12", "1")
    # Map ALIE onto the vdW surface; locate surface minima which predict the
    # most reactive electrophilic and radical attack sites
    SURFACE_ANALYSIS_ALIE = ("12", "2")
    # Compute the molecular vdW surface area and enclosed volume from the
    # rho=0.001 isosurface
    SURFACE_AREA_VOLUME = ("12", "3")
    # Becke-partition surface analysis: map a function onto each atomic Becke
    # surface and compute per-atom surface descriptors
    BECKE_SURFACE = ("12", "4")
    # Hirshfeld surface analysis for crystal packing studies: generate the
    # Hirshfeld surface and compute shape-index, curvedness, and contact-
    # distance properties
    HIRSHFELD_SURFACE = ("12", "5")
    # Locate and list all local minima and maxima of the mapped function on the
    # molecular surface, with their coordinates and values
    SURFACE_EXTREMA = ("12", "6")
    # Generate a Hirshfeld surface fingerprint plot (2D histogram of d_i vs.
    # d_e distances); reveals intermolecular contact patterns in crystal
    # structures
    HIRSHFELD_SURFACE_FINGERPRINT = ("12", "5", "2")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 13: Process grid data
    # Works on grid data already held in memory (from Menu 5 cube generation)
    # or loaded from an external .cube/.grd file at startup.
    # ─────────────────────────────────────────────────────────────────────────
    # NEW SEQUENCES - TO BE TESTED
    # Export the current in-memory grid data to a Gaussian .cube file in the
    # working directory
    EXPORT_CUBE = ("13", "0")
    # Export all grid-point coordinates together with their function values to
    # plain-text file (output.txt)
    EXPORT_GRID_ALL_POINTS = ("13", "1")
    # Extract a 2D slice of the 3D grid in the XY plane at a user-specified Z
    # index and save to a text file
    GRID_EXTRACT_PLANE_XY = ("13", "2")
    # Extract a 2D slice of the 3D grid in the XZ plane at a user-specified Y
    # index and save to a text file
    GRID_EXTRACT_PLANE_XZ = ("13", "3")
    # Extract a 2D slice of the 3D grid in the YZ plane at a user-specified X
    # index and save to a text file
    GRID_EXTRACT_PLANE_YZ = ("13", "4")
    # Compute the planar average of grid values over all XY planes in a
    # specified Z range; outputs a 1-D average profile along Z
    GRID_AVERAGE_XY = ("13", "5")
    # Compute the planar average of grid values over all XZ planes in a
    # specified X range; outputs a 1-D average profile along X
    GRID_AVERAGE_XZ = ("13", "6")
    # Compute the planar average of grid values over all YZ planes in a
    # specified Y range; outputs a 1-D average profile along Y
    GRID_AVERAGE_YZ = ("13", "7")
    # Extract a 2D grid slice in the plane defined by the nuclear positions of
    # three user-specified atoms
    GRID_EXTRACT_PLANE_3ATOMS = ("13", "8")
    # Extract a 2D grid slice in the plane defined by three user-specified
    # Cartesian coordinates
    GRID_EXTRACT_PLANE_3POINTS = ("13", "9")
    # Export all grid points whose function value falls within a user-specified
    # [min, max] interval
    GRID_EXTRACT_VALUE_RANGE = ("13", "10")
    # Perform element-wise arithmetic (+, -, *, /) between the current grid and
    # a second cube file; used to construct difference densities, Fukui
    # functions, etc.
    GRID_MATH_OPERATIONS = ("13", "11")
    # Map function values from a second cube file onto the isosurface of the
    # current grid (e.g., colour an ELF isosurface by ESP)
    GRID_MAP_TO_ISOSURFACE = ("13", "12")
    # Set all grid points farther than (or closer than) a distance cutoff from
    # chosen atoms to a fixed value; used to mask unwanted isosurface regions
    GRID_SET_VALUE_DISTANCE = ("13", "13")
    # Set grid points outside the fuzzy overlap region of two user-defined
    # fragments to a fixed value; isolates the inter-fragment interaction zone
    GRID_SET_VALUE_FRAGMENT = ("13", "14")
    # Replace all grid values within a specified numerical range with a single
    # fixed value; applies thresholding or clamping to the data
    GRID_SET_VALUE_RANGE = ("13", "15")
    # Linearly rescale all grid values from their current [min, max] to a user-
    # specified [new_min, new_max]
    GRID_SCALE_RANGE = ("13", "16")
    # Print statistics (min, max, mean, std dev, and spatial integral) for grid
    # points in a user-defined spatial and/or value range
    GRID_STATISTIC_DATA = ("13", "17")
    # Compute and plot the running integral of the grid data along X, Y, or Z;
    # useful for generating charge-displacement curves in EDA
    GRID_PLOT_INTEGRAL_CURVE = ("13", "18")
    # Open the interactive GUI to visualise the isosurface of the currently
    # loaded grid data at an adjustable isovalue slider
    GRID_VISUALIZE_ISOSURFACE = ("13", "-2")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 14: Adaptive Natural Density Partitioning (AdNDP)
    # Decomposes the electron density into n-centre two-electron (nc-2e)
    # bonding elements interactively; visualises results in Multiwfn GUI.
    # NEED SPECIFIC FILE
    # ─────────────────────────────────────────────────────────────────────────
    # Launch the AdNDP interactive interface to search for 1c-2e lone pairs,
    # 2c-2e bonds, and multi-centre bonding elements; widely used for cluster
    # and aromatic-system bonding analysis
    ADNDP_ANALYSIS = ("14",)

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 15: Fuzzy atomic space analysis
    # Numerical integration of real-space functions in Becke or Hirshfeld
    # fuzzy atomic spaces; computes delocalization indices and aromaticity.
    # ─────────────────────────────────────────────────────────────────────────
    # NEW SEQUENCES - TO BE TESTED
    FUZZY_INTEGRATE_EDENSITY = ("15", "1", "1")
    FUZZY_INTEGRATE_NORM_RHO = ("15", "1", "2")
    FUZZY_INTEGRATE_LAPLACIAN = ("15", "1", "3")
    FUZZY_INTEGRATE_ORB_WFN_HOMO = ("15", "1", "4", "h")
    FUZZY_INTEGRATE_ORB_WFN_LUMO = ("15", "1", "4", "l")
    FUZZY_INTEGRATE_ESPIN_DENSITY = ("15", "1", "5")
    FUZZY_INTEGRATE_KR = ("15", "1", "6")
    FUZZY_INTEGRATE_GR = ("15", "1", "7")
    FUZZY_INTEGRATE_ESP_CHARGES = ("15", "1", "8")
    FUZZY_INTEGRATE_ELF = ("15", "1", "9")
    FUZZY_INTEGRATE_LOL = ("15", "1", "10")
    FUZZY_INTEGRATE_LOCAL_ENTROPY = ("15", "1", "11")
    FUZZY_INTEGRATE_ESP = ("15", "1", "12")
    FUZZY_INTEGRATE_RDG = ("15", "1", "13")
    FUZZY_INTEGRATE_RDG_PROMOLECULAR = ("15", "1", "14")
    FUZZY_INTEGRATE_LAMBDA2RHO = ("15", "1", "15")
    FUZZY_INTEGRATE_LAMBDA2RGO_PROMOLECULAR = ("15", "1", "16")
    FUZZY_INTEGRATE_ALIE = ("15", "1", "18")
    FUZZY_INTEGRATE_EDR = ("15", "1", "20")
    FUZZY_INTEGRATE_ORB_OVERLAP_DR = ("15", "1", "21")
    FUZZY_INTEGRATE_DELTAG_PROMOLECULAR = ("15", "1", "22")
    FUZZY_INTEGRATE_DELTAG_HIRSHFELD = ("15", "1", "23")
    FUZZY_INTEGRATE_IRI = ("15", "1", "24")

    FUZZY_MULTIPOLE = ("15", "2", "1")

    FUZZY_OVERLAP_EDENSITY = ("15", "8", "1", "n")
    FUZZY_OVERLAP_NORM_RHO = ("15", "8", "1", "n")
    FUZZY_OVERLAP_LAPLACIAN = ("15", "8", "3", "n")
    FUZZY_OVERLAP_ORB_WFN_HOMO = ("15", "8", "4", "h", "n")
    FUZZY_OVERLAP_ORB_WFN_LUMO = ("15", "8", "4", "l", "n")
    FUZZY_OVERLAP_ESPIN_DENSITY = ("15", "8", "5", "n")
    FUZZY_OVERLAP_KR = ("15", "8", "6", "n")
    FUZZY_OVERLAP_GR = ("15", "8", "7", "n")
    FUZZY_OVERLAP_ESP_CHARGES = ("15", "8", "8", "n")
    FUZZY_OVERLAP_ELF = ("15", "8", "9", "n")
    FUZZY_OVERLAP_LOL = ("15", "8", "10", "n")
    FUZZY_OVERLAP_LOCAL_ENTROPY = ("15", "8", "11", "n")
    FUZZY_OVERLAP_ESP = ("15", "8", "12", "n")
    FUZZY_OVERLAP_RDG = ("15", "8", "13", "n")
    FUZZY_OVERLAP_RDG_PROMOLECULAR = ("15", "8", "14", "n")
    FUZZY_OVERLAP_LAMBDA2RHO = ("15", "8", "15", "n")
    FUZZY_OVERLAP_LAMBDA2RGO_PROMOLECULAR = ("15", "8", "16", "n")
    FUZZY_OVERLAP_ALIE = ("15", "8", "18", "n")
    FUZZY_OVERLAP_EDR = ("15", "8", "20", "n")
    FUZZY_OVERLAP_ORB_OVERLAP_DR = ("15", "8", "21", "n")
    FUZZY_OVERLAP_DELTAG_PROMOLECULAR = ("15", "8", "22", "n")
    FUZZY_OVERLAP_DELTAG_HIRSHFELD = ("15", "8", "23", "n")
    FUZZY_OVERLAP_IRI = ("15", "8", "24", "n")

    CLRK_MATRIX = ("15", "9", "n")

    # Compute the Atomic Overlap Matrix (AOM) S_AB = integral(psi_i * psi_j *
    # w_A) dr; required input for NOCV, ETS-NOCV, and delocalization-index
    # calculations
    ATOMIC_OVERLAP_MATRIX = ("15", "3")
    # Compute localization index and pairwise delocalization index (DI) for
    # all atom pairs in fuzzy spaces; DI is a bond-order analogue grounded in
    # density-matrix theory
    LOCALIZATION_DELOCALIZATION_INDEX = ("15", "4", "n")
    # Para Delocalization Index (PDI): average DI between para-related carbon
    # atoms in a 6-membered ring; larger PDI indicates stronger aromaticity
    PDI_AROMATICITY = ("15", "5", "q")
    # Aromatic Fluctuation Index (FLU): measures deviation of DI from reference
    # bond values; value near 0 indicates a fully aromatic ring
    FLU_AROMATICITY = ("15", "6", "q")
    # FLU-pi: FLU computed using only pi-orbital contributions to the
    # delocalization index; more selective for pi-electron aromaticity
    # INTERACTIVE - REQUIRES USER INPUT
    FLU_PI_AROMATICITY = ("15", "7")
    # Compute the condensed linear response kernel chi_AB; measures how much
    # electron density at atom B responds to an external perturbation at atom A
    CONDENSED_LINEAR_RESPONSE = ("15", "9")
    # Compute the para linear response (PLR) aromaticity index from the
    # condensed linear response kernel; large PLR indicates strong aromatic
    # delocalization
    PARA_LINEAR_RESPONSE = ("15", "10", "q")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 16: Charge Decomposition Analysis (CDA)
    # Generalised CDA (GCDA): decomposes orbital interactions between two or
    # more fragments into donation, back-donation, repulsion, and residual.
    # INTERACTIVE
    # ─────────────────────────────────────────────────────────────────────────
    # Launch the CDA/GCDA: decompose charge transfer between fragments
    # into donation (d), back-donation (b), repulsion (r), and residual (Delta)
    # terms; plot the orbital interaction diagram
    # INTERACTIVE - REQUIRES USER INPUT
    CDA_ANALYSIS = ("16",)

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 17: Basin analysis
    # Locates attractors of a real-space function, builds gradient-following
    # basins, and integrates properties (rho, multipoles, DI) within each.
    # ─────────────────────────────────────────────────────────────────────────
    # NEW SEQUENCES - TO BE TESTED
    BASIN_ANALYSIS_RHO = (
        "17",
        "1",
        "1",
        "2",
    )

    BASIN_EDENSITY = (
        "17",
        "1",
        "1",
        "1",
    )
    BASIN_NORM_RHO = ("17", "1", "2", "n")
    BASIN_LAPLACIAN = ("17", "1", "3", "n")
    BASIN_ORB_WFN_HOMO = ("17", "1", "4", "h", "n")
    BASIN_ORB_WFN_LUMO = ("17", "1", "4", "l", "n")
    BASIN_ESPIN_DENSITY = ("17", "1", "5", "n")
    BASIN_KR = ("17", "1", "6", "n")
    BASIN_GR = ("17", "1", "7", "n")
    BASIN_ESP_CHARGES = ("17", "1", "8", "n")
    BASIN_ELF = ("17", "1", "9", "n")
    BASIN_LOL = ("17", "1", "10", "n")
    BASIN_LOCAL_ENTROPY = ("17", "1", "11", "n")
    BASIN_ESP = ("17", "1", "12", "n")
    BASIN_RDG = ("17", "1", "13", "n")
    BASIN_RDG_PROMOLECULAR = ("17", "1", "14", "n")
    BASIN_LAMBDA2RHO = ("17", "1", "15", "n")
    BASIN_LAMBDA2RGO_PROMOLECULAR = ("17", "1", "16", "n")
    BASIN_ALIE = ("17", "1", "18", "n")
    BASIN_EDR = ("17", "1", "20", "n")
    BASIN_ORB_OVERLAP_DR = ("17", "1", "21", "n")
    BASIN_DELTAG_PROMOLECULAR = ("17", "1", "22", "n")
    BASIN_DELTAG_HIRSHFELD = ("17", "1", "23", "n")
    BASIN_IRI = ("17", "1", "24", "n")

    # AIM basin analysis: attractors are nuclei; integrates rho in each QTAIM
    # atomic basin to yield Bader/AIM charges, atomic multipoles, LI, and DI
    BASIN_ANALYSIS_AIM = ("17", "1")
    # ELF basin analysis: finds ELF attractors (bonding pairs, lone pairs, core
    # shells) and integrates rho per basin to obtain basin electron populations
    BASIN_ANALYSIS_ELF = ("17", "2")
    # Integrate any user-chosen real-space function over basins that have
    # already been generated in the current Multiwfn session
    BASIN_INTEGRATE_PROPERTY = ("17", "3")
    # ESP basin analysis: ESP minima serve as attractors; integrates rho in
    # corresponding basins to characterise nucleophilic pockets and lone-pair
    # regions
    BASIN_ANALYSIS_ESP = ("17", "4")
    # LOL basin analysis: LOL maxima serve as attractors (bonding and lone-pair
    # regions); integrates rho for sharp-boundary basin populations
    BASIN_ANALYSIS_LOL = ("17", "5")
    # LOL-alpha basin analysis restricted to alpha-spin electrons; isolates pi-
    # electron basins for aromaticity and pi-conjugation studies
    BASIN_ANALYSIS_LOL_ALPHA = ("17", "6")
    # Basin analysis using a user-defined real-space function specified via the
    # iuserfunc parameter in settings.ini
    BASIN_ANALYSIS_CUSTOM = ("17", "0")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 18: Electron excitation analysis
    # Analyses TD-DFT/CIS excited states loaded from Gaussian or ORCA output.
    # REQUIRES BOTH FCH AND OUT FILES - IGNORE FOR NOW
    # ─────────────────────────────────────────────────────────────────────────
    # Analyse hole and electron density distributions for a TD-DFT transition:
    # compute centroid distance, t index, and Sr overlap integral; visualise
    # hole/electron isosurfaces
    HOLE_ELECTRON_ANALYSIS = ("18", "1")
    # Plot the transition density matrix as a colour-filled 2D map in the atom-
    # atom/MO-MO representation; reveals coherence length and charge-transfer
    # character
    TRANSITION_DENSITY_MATRIX = ("18", "2")
    # Analyse charge transfer from an electron-density-difference grid: compute
    # CT distance, amount of transferred charge per fragment using the Plasser-
    # Lischka method
    CHARGE_TRANSFER_ANALYSIS = ("18", "3")
    # Compute the Delta_r index (charge-transfer length) for each excited state
    # to quantify local (small Delta_r) vs. charge-transfer (large Delta_r)
    # excitation character
    DELTA_R_INDEX = ("18", "4")
    # Calculate and print transition electric (and optionally magnetic) dipole
    # moments between all pairs of excited states in the loaded TD-DFT output
    TRANSITION_DIPOLE_MOMENTS = ("18", "5")
    # Generate Natural Transition Orbitals (NTOs) for a chosen excitation; the
    # dominant particle-hole NTO pair visually captures the essential character
    # of the transition
    GENERATE_NTO = ("18", "6")
    # Inter-Fragment Charge Transfer (IFCT): quantify hole and electron
    # populations transferred between user-defined fragments upon
    # photoexcitation
    IFCT_ANALYSIS = ("18", "8")
    # Compute the Lambda diagnostic for each TD-DFT transition; Lambda < 0.3
    # flags charge-transfer states that may be poorly described by standard
    # functionals
    LAMBDA_INDEX = ("18", "14")
    # Charge Transfer Spectrum (CTS): run batch IFCT over all excited states,
    # then plot the result as an inter-fragment CT contribution spectrum
    CTS_ANALYSIS = ("18", "16")
    # Compute the conditional electron density: given one electron is fixed
    # at reference point r0, maps the conditional probability of finding a
    # second electron at r
    CONDITIONAL_DENSITY = ("18", "17")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 19: Orbital localization
    # Transforms canonical delocalized MOs into spatially localised MOs
    # for chemical interpretation as bonds, lone pairs, and core orbitals.
    # ─────────────────────────────────────────────────────────────────────────
    # NEW SEQUENCES - TO BE TESTED - FORCED DEFAULT TO HIRSHFELD COMPOSITION
    PIPEK_MEZEY_LOCALIZATION_HIRSHFELD_OCCUPIED = ("19", "1")
    PIPEK_MEZEY_LOCALIZATION_HIRSHFELD_ALL = ("19", "2")
    PIPEK_MEZEY_LOCALIZATION_LOWDIN_OCUPIED = ("19", "-6", "2", "1")
    PIPEK_MEZEY_LOCALIZATION_LOWDIN_ALL = ("19", "-6", "2", "2")
    PIPEK_MEZEY_LOCALIZATION_BECKE_OCCUPIED = ("19", "-6", "3", "1")
    PIPEK_MEZEY_LOCALIZATION_BECKE_ALL = ("19", "-6", "3", "2")
    BOYS_LOCALIZATION_OCCUPIED = ("19", "-6", "10", "1")
    BOYS_LOCALIZATION_ALL = ("19", "-6", "10", "2")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 20: Weak interaction analysis
    # Methods based on the Reduced Density Gradient (RDG) and related
    # functions to visualise and quantify non-covalent interactions (NCI).
    # ─────────────────────────────────────────────────────────────────────────
    # NCI analysis from wavefunction: generate RDG and sign(lambda2)*rho cubes;
    # render RDG isosurface coloured by sign(lambda2)*rho in VMD to visualise
    # H-bonds, vdW contacts, and steric clashes
    NCI_ANALYSIS = ("20", "1", "2", "1", "2")
    # Promolecular NCI: same as NCI but built from superposition of free-atom
    # densities; very fast and suitable for macromolecules and protein-ligand
    # binding without a full wavefunction
    NCI_PROMOLECULAR = ("20", "2", "2", "1", "2")
    # Interaction Region Indicator (IRI): improved NCI variant that decays
    # smoothly at atomic cores and better resolves weak intermolecular contacts
    # recommended over standard RDG
    IRI_ANALYSIS = ("20", "4", "2", "1", "2")
    # Density Overlap Regions Indicator (DORI): highlights regions where two
    # electron densities overlap, clearly revealing both covalent and non-
    # covalent interaction zones
    DORI_ANALYSIS = ("20", "5", "2", "1", "2")
    # Compute and visualise the van der Waals interaction potential landscape
    # around the molecule
    VDW_POTENTIAL = ("20", "6", "2", "1", "2")

    # Averaged NCI (ANCI): NCI analysis averaged over an ensemble of MD
    # trajectory snapshots; reveals which non-covalent interactions persist
    # statistically over time
    # INTERACTIVE - REQUIRES USER INPUT
    ANCI_ANALYSIS = ("20", "3", "2", "1", "2")

    # INTERACTIVE - REQUIRES USER INPUT
    # Independent Gradient Model (IGM): decomposes the density gradient into
    # intra- and inter-fragment parts; isolates and visualises inter-fragment
    # non-covalent interactions
    # INTEACTIVE - REQUIRES USER INPUT
    IGM_ANALYSIS = ("20", "10")
    # IGM-H (Hirshfeld-based IGM): uses Hirshfeld atomic densities as the
    # reference for gradient decomposition; more accurate inter-fragment
    # gradient isolation than standard IGM
    # INTERACTIVE - REQUIRES USER INPUT
    IGMH_ANALYSIS = ("20", "11")
    # Averaged IGM (aIGM): IGM analysis averaged over an MD trajectory ensemble
    # identifies which inter-fragment interactions are persistent in dynamic
    # systems
    # INTERACTIVE - REQUIRES USER INPUT
    AIGM_ANALYSIS = ("20", "12")
    # Modified IGM (mIGM): coordinate-only IGM variant that avoids the need for
    # a wavefunction; nearly identical results to IGMH for weak interactions at
    # orders-of-magnitude lower cost
    # INTERACTIVE - REQUIRES USER INPUT
    MIGM_ANALYSIS = ("20", "-10")
    # Averaged mIGM (amIGM): mIGM averaged over MD snapshots; extends mIGM to
    # fluctuation environments; more robust than aNCI and recommended for MD-
    # based NCI studies
    # INTERACTIVE - REQUIRES USER INPUT
    AMIGM_ANALYSIS = ("20", "-11")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 21: Energy Decomposition Analysis (EDA)
    # INTERACTIVE OR REQUIRES SPECIFIC FILE TYPES
    # ─────────────────────────────────────────────────────────────────────────
    # Simple EDA using combined fragment wavefunctions: decomposes interaction
    # energy into electrostatic, Pauli exchange-repulsion, polarisation, and
    # dispersion components
    EDA_FF = ("21", "1")
    # Symmetry-based localisation EDA (SBL-EDA): energy partitioning scheme
    # using symmetry constraints on localised fragment orbitals
    EDA_SBL = ("21", "2")
    # Second-Order Bond Energy Analysis (SOBEDA): decomposes pairwise bond
    # energies into individual occupied-MO contributions via second-order
    # perturbation theory
    SOBEDA_ANALYSIS = ("21", "3")
    # Compute per-atom contributions to the DFT-D3 or D4 empirical dispersion
    # correction energy; identifies which atoms drive van der Waals attraction
    # most strongly
    DISPERSION_ATOMIC_CONTRIBUTION = ("21", "4")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 22: Conceptual DFT (CDFT)
    # Chemical reactivity descriptors derived from DFT response theory.
    # NEED MULTIPLE FILES
    # ─────────────────────────────────────────────────────────────────────────
    # Launch the CDFT module: compute global reactivity indices including
    # chemical potential mu, chemical hardness eta, softness S, and
    # electrophilicity index omega
    CDFT_ANALYSIS = ("22",)
    # Condense Fukui functions to atomic values using Mulliken, Hirshfeld,
    # Becke, or AIM partitioning; tabulates f+, f-, and f0 for each atom
    CONDENSED_FUKUI = ("22", "3")
    # Compute local hardness eta(r) from the local chemical potential;
    # complements the dual descriptor for predicting hard-soft reactivity site
    # preferences
    LOCAL_HARDNESS = ("22", "4")
    # Compute orbital weights as a reactivity map; low-weight regions indicate
    # preferred electrophilic attack sites
    ORBITAL_WEIGHTS = ("22", "5")

    ORBITAL_WEIGHTED_FUKUI = ("22", "6")
    GRID_FUKUI = ("22", "7", "2", "5", "6", "7", "8")
    SUPERDELOCALIZABILITIES_NUC_E = ("22", "8")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 23: ETS-NOCV
    # Extended Transition State combined with Natural Orbitals for Chemical
    # Valence: energy-decomposed orbital interaction analysis.
    # INTERACTIVE OR REQUIRES SPECIFIC FILE TYPES
    # ─────────────────────────────────────────────────────────────────────────
    # ETS-NOCV: diagonalise the deformation density matrix to obtain NOCV
    # pairs; decompose the orbital interaction energy Delta_E_orb into
    # individual pairwise NOCV channel contributions; useful for characterising
    # sigma, pi, and delta bond formation
    ETS_NOCV_ANALYSIS = ("23",)

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 24: Polarizability
    # INTERACTIVE
    # ─────────────────────────────────────────────────────────────────────────
    # Parse and print polarizability (alpha) and hyperpolarizability (beta,
    # gamma) tensors from a Gaussian frequency or finite-field task output file
    PARSE_POLARIZABILITY = ("24", "1")
    # Compute (hyper)polarizability via the Sum-Over-States (SOS) method from
    # TD-DFT excited-state data; useful when a direct field-perturbation
    # calculation is impractical
    SOS_POLARIZABILITY = ("24", "2")
    # Compute and visualise the polarizability density p(r): the spatial
    # distribution of where molecular polarisability originates in real space
    POLARIZABILITY_DENSITY = ("24", "3")
    # Project the polarizability tensor onto a unit sphere to produce a 3D
    # directional surface showing polarisability anisotropy
    UNIT_SPHERE_POLARIZABILITY = ("24", "5")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 25: Aromaticity
    # Magnetic, geometric, and electronic indices for aromaticity assessment.
    # REQUIRES GAUSSIAN
    # ─────────────────────────────────────────────────────────────────────────
    # Anisotropy of the Induced Current Density: compute the magnetically
    # induced ring-current density and export as a 3D isosurface; diatropic
    # current = aromatic
    AICD_ANALYSIS = ("25", "1")
    # Compute Nucleus-Independent Chemical Shift (NICS) at a single user-
    # specified point; negative NICS = diatropic / aromatic, positive =
    # paratropic / antiaromatic
    NICS_POINT = ("25", "2")
    # Iso-Chemical Shielding Surface (ICSS): compute NMR shielding on a 3D grid
    # and generate isosurfaces to visualise the spatial extent and shape of
    # ring-current shielding cones
    ICSS_ANALYSIS = ("25", "3")
    # NICS scan along a user-defined path perpendicular to a ring plane; plots
    # the NICS(z) profile to show how aromaticity decays with distance from the
    # ring
    NICS_SCAN = ("25", "4")
    # Bird aromaticity index (I5 for 5-membered, I6 for 6-membered rings):
    # geometry-based index from bond-length uniformity; 100 = fully aromatic
    # reference
    BIRD_INDEX = ("25", "5")
    # Harmonic Oscillator Model of Aromaticity: geometry-based index; 1 =
    # fully aromatic, 0 = nonaromatic, negative values indicate antiaromaticity
    HOMA_INDEX = ("25", "6")
    # HOMAC (corrected HOMA) and HOMER aromaticity indices; improved geometric
    # models that account for updated bond-length reference values from crystal
    # data
    HOMAC_HOMER = ("25", "7")
    # Stanger EN geometric aromaticity index: separates bond-length alternation
    # (EN_GEO) from uniform bond-length deviation (EN_BLA) for a cleaner
    # aromaticity measure
    STANGER_INDEX = ("25", "8")
    # Automated 1-D NICS scan perpendicular to a ring: places ghost atoms at
    # incremental heights, runs Gaussian, reads back shielding tensors, and
    # plots the NICS(z) curve
    NICS_1D_SCAN = ("25", "13")
    # Automated 2D NICS shielding map in or perpendicular to the ring plane;
    # produces a colour-filled map to spatially visualise the magnetic ring-
    # current delocalization
    NICS_2D_MAP = ("25", "14")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 100: Utilities Part 1
    # INTERACTIVE
    # ─────────────────────────────────────────────────────────────────────────
    # Plot a 2D scatter graph of two real-space functions evaluated on the same
    # grid and export both as .cube files; classic use is RDG vs.
    # sign(lambda2)*rho for the NCI scatter plot
    SCATTER_GRAPH_TWO_FUNCTIONS = ("100", "1")
    # Export the loaded wavefunction or geometry to .pdb, .xyz, .wfn, .molden,
    # .fch, Gaussian input, GAMESS input, or CP2K input format
    EXPORT_VARIOUS_FILES = ("100", "2")
    # Calculate the molecular van der Waals volume as the space enclosed within
    # the rho = 0.001 a.u. electron-density isosurface
    VDW_VOLUME = ("100", "3")
    # Numerically integrate a chosen real-space function over all space using
    # Becke multi-centre quadrature; used to verify total electron count
    # (integral of rho = N)
    INTEGRATE_WHOLE_SPACE = ("100", "4")
    # Compute the spatial overlap integral <psi_i_alpha | psi_j_beta> between
    # each alpha-beta orbital pair; measures how well alpha and beta orbitals
    # correspond spatially in UHF/UKS
    ORBITAL_OVERLAP_INTEGRAL = ("100", "5")
    # Parse a Gaussian output file and plot the SCF energy and DIIS error
    # convergence as a function of iteration number
    MONITOR_SCF_CONVERGENCE = ("100", "6")
    # Generate a Gaussian input file that uses the already-converged MO
    # coefficients as initial guess; saves SCF iterations when restarting or
    # changing basis
    GAUSSIAN_INITIAL_GUESS = ("100", "7")
    # Generate a Gaussian input file whose initial MO guess is assembled from
    # separately converged fragment wavefunctions; used to prepare CDA and EDA
    # calculations
    FRAGMENT_GUESS_INPUT = ("100", "8")
    # Evaluate and print the coordination number for every atom using an
    # empirical distance criterion scaled from vdW or covalent radii
    ATOMIC_COORDINATION = ("100", "9")
    # Calculate integral of |psi_i(r)| * |psi_j(r)| over all space; measures
    # how much two orbital densities spatially overlap without phase
    # cancellation
    ORBITAL_OVERLAP_CENTROID = ("100", "11")
    # Biorthogonalise a set of orbitals (e.g., localised MOs) and compute their
    # one-electron energies from the Fock matrix without re-running the SCF
    BIORTHOGONALIZATION = ("100", "12")
    # Compute HOMA and Bird aromaticity indices for a user-specified ring
    # directly from the molecular geometry coordinates
    HOMA_BIRD_AROMATICITY = ("100", "13")
    # Compute the LOLIPOP index (LOL Integrated Pi Over Plane): integral of LOL
    # in the pi-symmetry plane above and below the ring; quantifies pi-electron
    # delocalization
    LOLIPOP_INDEX = ("100", "14")
    # Calculate intermolecular orbital overlap integrals between two molecular
    # species; predicts charge-transfer coupling strengths relevant to organic
    # semiconductor design
    INTERMOLECULAR_OVERLAP = ("100", "15")
    # Construct and export the Fock matrix F = C^{-1T} * epsilon * C^{-1} from
    # orbital energies and coefficients; needed for biorthogonalisation and
    # ETS-NOCV
    GENERATE_FOCK_MATRIX = ("100", "17")
    # Yoshizawa electron-transport route analysis: identify the dominant
    # through-bond tunnelling pathway between two terminal atoms using squared
    # orbital-amplitude products
    ELECTRON_TRANSPORT_ROUTE = ("100", "18")
    # Combine separately computed fragment wavefunctions at their molecular
    # geometry to form a promolecular .wfn file;
    # reference state for EDA and CDA
    COMBINE_FRAGMENTS = ("100", "19")
    # Compute Hellmann-Feynman electrostatic forces on each nucleus from
    # the electron density distribution and all other nuclear charges
    HELLMANN_FEYNMAN_FORCES = ("100", "20")
    # Print selected geometric properties: bond lengths, bond angles, dihedral
    # angles, and distances for user-specified atoms
    GEOMETRY_PROPERTIES = ("100", "21")
    # Automatically identify pi-symmetry orbitals in planar systems based on
    # nodal-plane criteria and optionally zero their occupation numbers to
    # isolate the sigma frame
    DETECT_PI_ORBITALS = ("100", "22")
    # Fit the distribution of a real-space function onto per-atom parameter
    # values using a least-squares procedure; useful for atom-centred property
    # assignments
    FIT_FUNCTION_TO_ATOMS = ("100", "23")
    # Compute the out-of-plane NICS_ZZ shielding component for non-planar rings
    # by projecting the full shielding tensor onto the local ring-normal
    # direction
    NICS_ZZ_NONPLANAR = ("100", "24")
    # Calculate the area enclosed by and the perimeter of a ring defined by
    # user-specified atom indices from the molecular geometry
    RING_AREA_PERIMETER = ("100", "25")
    # Generate a CP2K periodic-DFT input file from the loaded geometry and cell
    # information; enter 'cp2k' at the export-type prompt to activate
    GENERATE_CP2K_INPUT = ("100", "2")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 200: Utilities Part 2
    # ─────────────────────────────────────────────────────────────────────────
    # NEW SEQUENCES - TO BE TESTED - FORCED DEFAULT TO HIRSHFELD COMPOSITION
    ORBITAL_INTEGRAL_ELECTRIC_DIPOLE = ("200", "10", "1", "3")
    ORBITAL_INTEGRAL_MAGNETIC_DIPOLE = ("200", "10", "2", "3")
    ORBITAL_INTEGRAL_VELOCITY = ("200", "10", "3", "3")
    ORBITAL_INTEGRAL_KINETIC_ENERGY = ("200", "10", "4", "3")
    ORBITAL_INTEGRAL_OVERLAP = ("200", "10", "5", "3")

    SPACIAL_DELOCALISATION_EDENSITY = ("200", "19", "1", "1")
    SPACIAL_DELOCALISATIOn_NORM_RHO = ("200", "19", "1", "2")
    SPATIAL_DELOCALISATION_LAPLACIAN = ("200", "19", "1", "3")
    SPATIAL_DELOCALISATION_ORB_WFN = ("200", "19", "1", "4")
    SPATIAL_DELOCALISATION_ESPIN_DENSITY = ("200", "19", "1", "5")
    SPATIAL_DELOCALISATION_KR = ("200", "19", "1", "6")
    SPATIAL_DELOCALISATION_GR = ("200", "19", "1", "7")
    SPATIAL_DELOCALISATION_ESP_CHARGES = ("200", "19", "1", "8")
    SPATIAL_DELOCALISATION_ELF = ("200", "19", "1", "9")
    SPATIAL_DELOCALISATION_LOL = ("200", "19", "1", "10")
    SPATIAL_DELOCALISATION_LOCAL_ENTROPY = ("200", "19", "1", "11")
    SPATIAL_DELOCALISATION_ESP_TOTAL = ("200", "19", "1", "12")
    SPATIAL_DELOCALISATION_RDG = ("200", "19", "1", "13")
    SPATIAL_DELOCALISATION_RDG_PROMOLECULAR = ("200", "19", "1", "14")
    SPATIAL_DELOCALISATION_LAMBDA2RHO = ("200", "19", "1", "15")
    SPATIAL_DELOCALISATION_LAMBDA2RHO_PROMOLECULAR = ("200", "19", "1", "16")
    SPATIAL_DELOCALISATION_ALIE = ("200", "19", "1", "18")
    SPATIAL_DELOCALISATION_EDR = ("200", "19", "1", "20")
    SPATIAL_DELOCALISATION_ORB_OVERLAP_DR = ("200", "19", "1", "21")
    SPATIAL_DELOCALISATION_DELTAG_PROMOLECULAR = ("200", "19", "1", "22")
    SPATIAL_DELOCALISATION_DELTAG_HIRSHFELD = ("    200", "19", "1", "23")
    SPATIAL_DELOCALISATION_IRI = ("200", "19", "1", "24")
    # Fluctuation NCI (aRDG): RDG-based weak-interaction analysis time-averaged
    # over an MD trajectory; reveals which non-covalent contacts persist
    # throughout the simulation
    CVB_INDEX = ("200", "1")
    # Compute atomic and bond dipole moments in Hilbert space using Mulliken
    # partitioning; decomposes the total molecular dipole into atomic
    # and inter-atomic bond contributions
    ATOMIC_BOND_DIPOLES = ("200", "2")
    # Generate Gaussian .cube files for multiple orbitals in a single automated
    # run; each selected orbital wavefunction is exported as a separate file
    MULTIPLE_ORBITAL_CUBES = ("200", "3")
    # Generate 3D NMR shielding (ICSS) data as .cube files; visualise
    # magnetically induced ring-current effects as shielding isosurfaces in VMD
    ICSS_CUBES = ("200", "4")
    # Compute and plot the spherical radial distribution function RDF(r) =
    # 4*pi*r^2 * f(r) of a chosen real-space function; useful for spherically
    # symmetric systems such as fullerenes or atoms
    RADIAL_DISTRIBUTION = ("200", "5")
    # Analyse the overlap correspondence between orbitals in two different
    # wavefunctions; tracks how orbitals evolve along a reaction path or change
    # between two levels of theory
    ORBITAL_CORRESPONDENCE = ("200", "6")
    # Parse and tabulate polarizability (alpha) and first/second/third
    # hyperpolarizability (beta, gamma, delta) tensors from a Gaussian polar-
    # keyword output file
    PARSE_POLARIZABILITY_GAUSSIAN = ("200", "7")
    # Compute polarizability and 1st/2nd/3rd hyperpolarizabilities by the Sum-
    # Over-States (SOS) method using TD-DFT excited-state energies and
    # transition moments
    SOS_HYPERPOLARIZABILITY = ("200", "8")
    # Calculate and print average bond lengths and average coordination numbers
    # for all element-type combinations in the molecule
    AVERAGE_BOND_LENGTH = ("200", "9")
    # Compute one-electron integrals between pairs of orbitals: kinetic energy
    # <i|T|j>, nuclear attraction <i|V|j>, Coulomb, and exchange integrals
    ORBITAL_INTEGRALS = ("200", "10")
    # Compute the centroid (first moment), second moments, and radius of
    # gyration of a chosen real-space function; characterises its spatial
    # spread and shape
    FUNCTION_MOMENTS = ("200", "11")
    # Compute the Energy Index (EI) and Bond Polarity Index (BPI) from kinetic
    # and potential energy densities evaluated at bond critical points
    ENERGY_INDEX = ("200", "12")
    # Decompose a grid-based property (e.g., electron density) into
    # contributions from individual MOs and export per-orbital .cube files for
    # detailed analysis
    ORBITAL_CONTRIBUTIONS_TO_GRID = ("200", "13")
    # Identify and characterise topologically connected electron-density
    # domains: contiguous regions above a chosen isovalue threshold in the 3D
    # grid
    DOMAIN_ANALYSIS = ("200", "14")
    # Compute a spatial correlation index between two real-space functions
    # evaluated on the same grid; quantifies how similarly their distributions
    # are arranged in space
    CORRELATION_INDEX = ("200", "15")
    # Diagonalise the one-particle density matrix to obtain natural orbitals
    # (NOs) and their occupation numbers; useful for characterising electron
    # correlation and multi-reference character
    NATURAL_ORBITALS = ("200", "16")
    # Evaluate two-electron Coulomb (J_ij) and exchange (K_ij) integrals
    # between all pairs of molecular orbitals; provides input for EDA,
    # perturbation theory, and excited-state analysis
    COULOMB_EXCHANGE_INTEGRALS = ("200", "17")
    # Compute Bond Length Alternation (BLA) and Bond Order Alternation (BOA)
    # along a user-defined conjugated chain; quantifies the degree of
    # alternation in polyene or polyacetylene-type pi systems
    BLA_BOA_ANALYSIS = ("200", "18")
    # Compute the spatial delocalization index (SDI): measures how broadly the
    # electron density is spread across space; indicator of overall
    # delocalisation
    SPATIAL_DELOCALIZATION_INDEX = ("200", "19")
    # Bond Order Decomposition (BOD) and Natural Atomic Dipole Orbital (NADO)
    # analysis: decomposes Mayer bond orders and atomic dipole moments into
    # natural orbital pair contributions
    BOD_NADO_ANALYSIS = ("200", "20")
    # Perform Lowdin symmetric orthogonalization of the current orbital set and
    # optionally evaluate the resulting orbital energies directly from the Fock
    # matrix
    LOWDIN_ORTHOGONALIZATION = ("200", "21")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 300: Utilities Part 3
    # ─────────────────────────────────────────────────────────────────────────
    # Calculate the free void volume in a periodic unit cell not occupied by
    # atomic vdW spheres; relevant for porosity characterisation in MOFs,
    # zeolites, and porous organic cages
    FREE_VOLUME_IN_CELL = ("300", "1")
    # Fit the spherically averaged radial electron density of each atom in the
    # molecule to a set of Slater-type or Gaussian-type radial functions;
    # generates transferable atomic density parameters
    FIT_ATOMIC_RADIAL_DENSITY = ("300", "2")
    # Simulate a Scanning Tunnelling Microscopy (STM) image using the Tersoff-
    # Hamann approximation from the local density of states near the Fermi
    # energy
    STM_IMAGE = ("300", "4")
    # Compute and print the molecular electric multipole moment tensor up to
    # hexadecapole order from the electron density and nuclear charges
    ELECTRIC_MULTIPOLE_MOMENTS = ("300", "5")
    # Evaluate orbital energies epsilon_i = <psi_i|F|psi_i> directly from the
    # Fock matrix without running a new SCF; particularly useful after orbital
    # localisation or biorthogonalisation
    ORBITAL_ENERGIES_FROM_FOCK = ("300", "6")
    # Perform geometric transformations on the molecular structure: translate,
    # rotate, reflect, invert, or generate symmetry-equivalent atoms
    GEOMETRY_OPERATIONS = ("300", "7")
    # Generate a molecular surface distance projection map plotting d_i
    # (distance to nearest internal atom) vs. d_e (distance to nearest external
    # atom) for every surface point; reveals steric complementarity and packing
    # motifs
    SURFACE_DISTANCE_PROJECTION = ("300", "8", "0")
    # Determine the Fermi energy level from the orbital energy distribution and
    # occupation numbers; particularly useful for periodic or large cluster
    # systems
    DETERMINE_FERMI_LEVEL = ("300", "9")


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
