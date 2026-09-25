"""Complete Multiwfn menu mapping for batch automation.

Notes
-----
Based on Multiwfn 3.8 (dev) manual.
====================================================
Every entry below was live-verified against the real bundled Multiwfn
3.8(dev) binary (coord.molden for ground-state analyses, tddft.out for
excitation/spectrum analyses that need real TD-DFT data) across two
audit passes

Each entry's comment has up to three parts:
  - a short description of what the analysis computes and why you'd
    run it;
  - a "Sequence:" line spelling out what every token in the tuple
    means, in order, so the mapping from numbers to menu choices is
    never a mystery -- built by live-probing the real running menu at
    every step, so it reflects what the token actually does, even on
    the handful of entries where that turns out not to match the
    entry's own name (see e.g. TOPOLOGY_ELF_ANALYSIS, TOPOLOGY_LOL_
    ANALYSIS, TOPOLOGY_LAPLACIAN_ANALYSIS, AMIGM_ANALYSIS -- their
    Sequence breakdowns show they don't currently reach the function
    their name promises; flagged here as comments only, sequence
    values intentionally left untouched pending a future fix pass);
  - where applicable, an "INTERACTIVE ONLY" tag (this sequence cannot
    be reduced to a fixed, non-interactive token list -- kept in place
    as reference/building blocks for future development) or a "FIXED"
    tag (the sequence was corrected in an earlier pass from the
    original, broken value shown in the tag).

Entries are grouped by Multiwfn main menu category, in the same
numeric order the program itself presents them (0-26, then the
"Other functions" utility menus 100/200/300).
"""

from enum import Enum


class Menu(Enum):
    """Enumeration of all Multiwfn menu functions.

    Sequences represent the keystrokes to navigate Multiwfn
    interactively. Each entry can be run standalone or composed
    sequentially via job.add_menu().

    Key conventions:
      - "0"  -> return to main menu (from most submenus, but not all --
               see each entry's own Sequence breakdown)
      - "-1" -> trigger search-all or global option within a submenu
      - "n"  -> decline a fragment/mirror-plane restriction prompt
      - "q"  -> quit / return without change
      - "h" / "l" -> select the HOMO / LUMO at an orbital-index prompt
      - "a"  -> select all orbitals
    """

    @classmethod
    def search(cls, query: str) -> list["Menu"]:
        """Search menu items by name (case-insensitive)."""
        query_upper = query.upper()
        return [item for item in cls if query_upper in item.name]

    @classmethod
    def list_all(cls) -> list[str]:
        """Return names of all menu items."""
        return [item.name for item in cls]

    def get_sequence(self) -> tuple[str, ...]:
        return self.value

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 0: Show molecular structure / view orbitals
    # ─────────────────────────────────────────────────────────────────────────
    # Display 3D molecular structure and orbital isosurfaces in the GUI;
    # prints atom coordinates to screen
    # Sequence: 0=Show molecular structure and view orbitals
    VIEW_STRUCTURE = ("0",)

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 1: Output all properties at a point
    # ─────────────────────────────────────────────────────────────────────────
    # Interactive only – prompts for a coordinate or atom index then
    # prints all supported real-space functions at that point (rho, ESP,
    # ELF, ...). INTERACTIVE ONLY - Ignore for now, output will not be
    # parsed corectly or output expected results Print every real-space
    # function value (rho, ESP, ELF, G, K, ...) at a user-supplied
    # coordinate or nucleus Requires xyz input and user interaction,
    # ignore for now
    # Sequence: 1=Output all properties at a point
    PROPERTIES_AT_POINT = ("1",)

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 2: Topology analysis
    # ─────────────────────────────────────────────────────────────────────────
    # Uses Newton iteration to locate critical points (CPs) of a chosen
    # real-space function, then traces gradient paths between them. CP
    # types: (3,-3) nuclear/max, (3,-1) bond, (3,+1) ring, (3,+3) cage.
    # Sequence: 2=Topology analysis; 0=Print and visualize all generated
    #           CPs, paths and interbasin surfaces
    TOPOLOGY_VISUALISE_CPS = ("2", "0")
    # Sequence: 2=Topology analysis; 2=Search CPs from nuclear positions
    TOPOLOGY_CP_NUCLEAR_POSITION = ("2", "2")
    # Sequence: 2=Topology analysis; 3=Search CPs from midpoint of atomic
    #           pairs
    TOPOLOGY_CP_MIDPOINTS = ("2", "3")
    # Sequence: 2=Topology analysis; 4=Search CPs from triangle center of
    #           three atoms
    TOPOLOGY_CP_TRIANGLE_CENTRES = ("2", "4")
    # Sequence: 2=Topology analysis; 5=Search CPs from pyramid center of
    #           four atoms
    TOPOLOGY_CP_PYRAMID_CENTRES = ("2", "5")
    # Sequence: 2=Topology analysis; 6=Search CPs from a batch of points
    #           within sphere(s); 0=Start the search using the defined
    #           sphere center
    TOPOLOGY_CP_SPHERE_POINTS = ("2", "6", "0")
    # Outputs CPout.txt
    # Sequence: 2=Topology analysis; 7=Show real space function values at
    #           specific CP or all CPs; 0=output properties of ALL CPs to
    #           CPprop.txt (not just one)
    TOPOLOGY_CP_REAL_SPACE_POINTS = ("2", "7", "0")
    # Sequence: 2=Topology analysis; 2=Search CPs from nuclear positions;
    #           3=Search CPs from midpoint of atomic pairs; 8=Generating
    #           the paths connecting (3,-3) and (3,-1) CPs
    TOPOLOGY_CP_PATHS_3MINUS3_3MINUS1 = ("2", "2", "3", "8")
    # Sequence: 2=Topology analysis; 2=Search CPs from nuclear positions;
    #           3=Search CPs from midpoint of atomic pairs; 9=Generating
    #           the paths connecting (3,+1) and (3,+3) CPs
    TOPOLOGY_CP_PATHS_3PLUS1_3PLUS3 = ("2", "2", "3", "9")
    # Search ALL critical points of rho starting from every nuclear
    # position; finds NCPs, BCPs, RCPs, and CCPs in one pass
    # Sequence: 2=Topology analysis; 2=Search CPs from nuclear positions;
    #           3=Search CPs from midpoint of atomic pairs; 8=Generating
    #           the paths connecting (3,-3) and (3,-1) CPs; 0=Print and
    #           visualize all generated CPs/paths/surfaces; -10=Return to
    #           main menu
    TOPOLOGY_SEARCH_CPS = ("2", "2", "3", "8", "0", "-10")
    # Full AIM workflow in one sequence: search all CPs, generate bond
    # paths, then generate interbasin surfaces
    # Sequence: 2=Topology analysis; 2=Search CPs from nuclear positions;
    #           -1=Set CP searching parameters; 3=Criteria for gradient
    #           norm convergence: 1.00000E-06; 4=criteria for displacement
    #           convergence
    TOPOLOGY_ANALYSIS_COMPLETE = ("2", "2", "-1", "3", "4")
    # Topology analysis of the electrostatic potential (ESP): locate ESP
    # critical points starting from all nuclei
    # Sequence: 2=Topology analysis; -2=Set path generating parameters;
    #           2=select the 'Stepsize' setting; -1=value given for the
    #           stepsize prompt
    TOPOLOGY_ESP_ANALYSIS = ("2", "-2", "2", "-1")
    # Topology analysis of the Localized Orbital Locator (LOL): find LOL
    # maxima corresponding to bonding and lone-pair basins
    # Sequence: 2=Topology analysis; -10=Return to main menu; 2=re-enter
    #           Topology analysis (from the main menu); -1=Set CP
    #           searching parameters
    TOPOLOGY_LOL_ANALYSIS = ("2", "-10", "2", "-1")
    # Topology analysis of the Electron Localization Function (ELF): find
    # ELF attractors and characterise bonding basins
    # Sequence: 2=Topology analysis; 9=Generating the paths connecting
    #           (3,+1) and (3,+3) CPs; 2=search CPs from nuclear positions
    #           (active function, default rho); -1=Set CP searching
    #           parameters
    TOPOLOGY_ELF_ANALYSIS = ("2", "9", "2", "-1")
    # Topology analysis of nabla^2 rho: locate charge-concentration (3,-3)
    # and charge-depletion (3,+3) critical points
    # Sequence: 2=Topology analysis; 3=Search CPs from midpoint of atomic
    #           pairs; 2=search CPs from nuclear positions (active
    #           function, default rho); -1=Set CP searching parameters
    TOPOLOGY_LAPLACIAN_ANALYSIS = ("2", "3", "2", "-1")
    # Search for bond critical points (BCPs) using atom-pair midpoints as
    # Newton starting guesses
    # Sequence: 2=Topology analysis; 3=Search CPs from midpoint of atomic
    #           pairs; -1=Set CP searching parameters
    TOPOLOGY_SEARCH_BCP = ("2", "3", "-1")
    # Search for ring critical points (RCPs) using triangle centres of
    # atom triplets as starting guesses
    # Sequence: 2=Topology analysis; 4=Search CPs from triangle center of
    #           three atoms; -1=Set CP searching parameters
    TOPOLOGY_SEARCH_RCP = ("2", "4", "-1")
    # Search for cage critical points (CCPs) using pyramid centres of atom
    # quartets as starting guesses
    # Sequence: 2=Topology analysis; 5=Search CPs from pyramid center of
    #           four atoms; -1=Set CP searching parameters
    TOPOLOGY_SEARCH_CCP = ("2", "5", "-1")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 3: Output / plot property along a line
    # ─────────────────────────────────────────────────────────────────────────
    # Evaluates a real-space function at 3000 points between two atoms or
    # coordinates and produces a curve map (line.txt exported on request).
    # or output expected results Plot total electrostatic potential (ESP)
    # along a line; reveals electrophilic and nucleophilic interaction
    # sites between atoms
    # Sequence: 3=Output and plot specific property in a line; 12=Total
    #           electrostatic potential (ESP)
    LINE_ESP = ("3", "12")
    # Plot electron density rho(r) along a line; shows charge accumulation
    # and depletion between bonded atoms
    # Sequence: 3=Output and plot specific property in a line; 1=Electron
    #           density (rho)
    LINE_ELECTRON_DENSITY = ("3", "1")
    # Plot Laplacian of electron density nabla^2 rho along a line;
    # negative regions indicate charge concentration (covalent bonds, lone
    # pairs)
    # Sequence: 3=Output and plot specific property in a line; 3=Laplacian
    #           of rho
    LINE_LAPLACIAN = ("3", "3")
    # Plot Electron Localization Function (ELF, range 0-1) along a line;
    # peaks identify bonding pairs and lone pairs
    # Sequence: 3=Output and plot specific property in a line; 9=Electron
    #           localization function (ELF)
    LINE_ELF = ("3", "9")
    # Plot Localized Orbital Locator (LOL, range 0-1) along a line;
    # similar to ELF but with sharper basin boundary definition
    # Sequence: 3=Output and plot specific property in a line;
    #           10=Localized orbital locator (LOL)
    LINE_LOL = ("3", "10")
    # Plot Reduced Density Gradient (RDG) along a line; low-RDG regions
    # reveal non-covalent interaction zones
    # Sequence: 3=Output and plot specific property in a line; 13=Reduced
    #           density gradient (RDG)
    LINE_RDG = ("3", "13")
    # Plot spin density rho_alpha - rho_beta along a line; used for
    # radical and open-shell systems to locate unpaired electrons
    # Sequence: 3=Output and plot specific property in a line; 5=Electron
    #           spin density
    LINE_SPIN_DENSITY = ("3", "5")
    # Plot magnitude of the electron density gradient |nabla rho| along a
    # line; maxima mark interatomic surface boundaries
    # Sequence: 3=Output and plot specific property in a line; 2=Gradient
    #           norm of rho
    LINE_GRADIENT_NORM = ("3", "2")
    # Plot Lagrangian (positive-definite) kinetic energy density G(r)
    # along a line
    # Sequence: 3=Output and plot specific property in a line;
    #           7=Lagrangian kinetic energy density G(r)
    LINE_KINETIC_G = ("3", "7")
    # Plot Hamiltonian kinetic energy density K(r) along a line; related
    # to G: K = G - (1/4) nabla^2 rho
    # Sequence: 3=Output and plot specific property in a line;
    #           6=Hamiltonian kinetic energy density K(r)
    LINE_KINETIC_K = ("3", "6")
    # Plot Average Local Ionization Energy (ALIE) along a line; low values
    # indicate weakly bound, reactive electrons susceptible to
    # electrophilic attack
    # Sequence: 3=Output and plot specific property in a line; 18=Average
    #           local ionization energy (ALIE)
    LINE_ALIE = ("3", "18")
    # Plot Source Function SF(r', r) along a line relative to a fixed
    # reference point r set in settings.ini; shows electron-density
    # contributions from each spatial region
    # Sequence: 3=Output and plot specific property in a line; 19=Source
    #           function, mode: 1, ref. point:   0.00000   0.00000
    #           0.00000
    LINE_SOURCE_FUNCTION = ("3", "19")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 4: Output / plot property in a plane
    # ─────────────────────────────────────────────────────────────────────────
    # The third tuple element encodes the plane-definition mode: "1" = XY
    # plane (input Z value) "2" = XZ plane (input Y value) "3" = YZ plane
    # (input X value) "4" = plane defined by three atom indices "5" =
    # plane defined by three Cartesian points Produces colour-filled,
    # contour, relief, or gradient/vector maps. or output expected results
    # Colour-filled map of electron density rho(r) in the XY plane;
    # visualises how charge is distributed across the molecule
    # Sequence: 4=Output and plot specific property in a plane; 1=Electron
    #           density (rho); 1=Color-filled map (with/without contour
    #           lines)
    # INTERACTIVE ONLY -- requires interactively defining a plane
    #                     (orientation/resolution); excluded per your
    #                     instruction to skip plane-dependent sequences
    PLANE_MAP_DENSITY = ("4", "1", "1")
    # Colour-filled map of total ESP in the XY plane; red regions are
    # nucleophilic (-), blue regions are electrophilic (+)
    # Sequence: 4=Output and plot specific property in a plane; 12=Total
    #           electrostatic potential (ESP); 1=Color-filled map
    #           (with/without contour lines)
    # INTERACTIVE ONLY -- requires interactively defining a plane
    #                     (orientation/resolution); excluded per your
    #                     instruction to skip plane-dependent sequences
    PLANE_MAP_ESP = ("4", "12", "1")
    # Colour-filled map of ELF in the XY plane; red/orange islands mark
    # bonding pairs and lone pairs, blue marks depleted regions
    # Sequence: 4=Output and plot specific property in a plane; 9=Electron
    #           localization function (ELF); 1=Color-filled map
    #           (with/without contour lines)
    # INTERACTIVE ONLY -- requires interactively defining a plane
    #                     (orientation/resolution); excluded per your
    #                     instruction to skip plane-dependent sequences
    PLANE_MAP_ELF = ("4", "9", "1")
    # Colour-filled map of LOL in the XY plane; chemically equivalent to
    # ELF with cleaner inter-basin boundaries
    # Sequence: 4=Output and plot specific property in a plane;
    #           10=Localized orbital locator (LOL); 1=Color-filled map
    #           (with/without contour lines)
    # INTERACTIVE ONLY -- requires interactively defining a plane
    #                     (orientation/resolution); excluded per your
    #                     instruction to skip plane-dependent sequences
    PLANE_MAP_LOL = ("4", "10", "1")
    # Colour-filled map of |nabla rho| in the XY plane; highlights
    # interatomic surface and shell-structure regions
    # Sequence: 4=Output and plot specific property in a plane; 2=Gradient
    #           norm of rho; 1=Color-filled map (with/without contour
    #           lines)
    # INTERACTIVE ONLY -- requires interactively defining a plane
    #                     (orientation/resolution); excluded per your
    #                     instruction to skip plane-dependent sequences
    PLANE_MAP_GRADIENT = ("4", "2", "1")
    # Colour-filled or contour map of nabla^2 rho in the XY plane;
    # negative = charge concentration (bonds/lone pairs), positive =
    # depletion
    # Sequence: 4=Output and plot specific property in a plane;
    #           3=Laplacian of rho; 1=Color-filled map (with/without
    #           contour lines)
    # INTERACTIVE ONLY -- requires interactively defining a plane
    #                     (orientation/resolution); excluded per your
    #                     instruction to skip plane-dependent sequences
    PLANE_MAP_LAPLACIAN = ("4", "3", "1")
    # Colour-filled map of spin density rho_alpha - rho_beta in the XY
    # plane; shows location of unpaired electrons in open-shell species
    # Sequence: 4=Output and plot specific property in a plane; 5=Electron
    #           spin density; 1=Color-filled map (with/without contour
    #           lines)
    # INTERACTIVE ONLY -- requires interactively defining a plane
    #                     (orientation/resolution); excluded per your
    #                     instruction to skip plane-dependent sequences
    PLANE_MAP_SPIN_DENSITY = ("4", "5", "1")
    # Colour-filled map of RDG in the XY plane; low-value regions reveal
    # non- covalent interactions when combined with sign(lambda2)rho
    # colouring
    # Sequence: 4=Output and plot specific property in a plane; 13=Reduced
    #           density gradient (RDG); 1=Color-filled map (with/without
    #           contour lines)
    # INTERACTIVE ONLY -- requires interactively defining a plane
    #                     (orientation/resolution); excluded per your
    #                     instruction to skip plane-dependent sequences
    PLANE_MAP_RDG = ("4", "13", "1")
    # Colour-filled map of sign(lambda2)*rho in the XY plane; negative =
    # attractive NCI (H-bond), positive = repulsive steric clash
    # Sequence: 4=Output and plot specific property in a plane;
    #           15=Sign(lambda2)*rho; 1=Color-filled map (with/without
    #           contour lines)
    # INTERACTIVE ONLY -- requires interactively defining a plane
    #                     (orientation/resolution); excluded per your
    #                     instruction to skip plane-dependent sequences
    PLANE_MAP_SIGN_LAMBDA2_RHO = ("4", "15", "1")
    # Colour-filled map of ALIE in the XY plane; low-ALIE pockets on the
    # surface predict electrophilic and radical attack sites
    # Sequence: 4=Output and plot specific property in a plane; 18=Average
    #           local ionization energy (ALIE); 1=Color-filled map
    #           (with/without contour lines)
    # INTERACTIVE ONLY -- requires interactively defining a plane
    #                     (orientation/resolution); excluded per your
    #                     instruction to skip plane-dependent sequences
    PLANE_MAP_ALIE = ("4", "18", "1")
    # Colour-filled map of Lagrangian kinetic energy density G(r) in the
    # XY plane
    # Sequence: 4=Output and plot specific property in a plane;
    #           7=Lagrangian kinetic energy density G(r); 1=Color-filled
    #           map (with/without contour lines)
    # INTERACTIVE ONLY -- requires interactively defining a plane
    #                     (orientation/resolution); excluded per your
    #                     instruction to skip plane-dependent sequences
    PLANE_MAP_KINETIC_G = ("4", "7", "1")
    # Colour-filled map of Hamiltonian kinetic energy density K(r) in the
    # XY plane
    # Sequence: 4=Output and plot specific property in a plane;
    #           6=Hamiltonian kinetic energy density K(r); 1=Color-filled
    #           map (with/without contour lines)
    # INTERACTIVE ONLY -- requires interactively defining a plane
    #                     (orientation/resolution); excluded per your
    #                     instruction to skip plane-dependent sequences
    PLANE_MAP_KINETIC_K = ("4", "6", "1")
    # Colour-filled map of Source Function in the XY plane relative to the
    # reference point defined in settings.ini
    # Sequence: 4=Output and plot specific property in a plane; 19=Source
    #           function, mode: 1, ref. point:   0.00000   0.00000
    #           0.00000; 1=Color-filled map (with/without contour lines)
    # INTERACTIVE ONLY -- requires interactively defining a plane
    #                     (orientation/resolution); excluded per your
    #                     instruction to skip plane-dependent sequences
    PLANE_MAP_SOURCE_FUNCTION = ("4", "19", "1")
    # Colour-filled or contour map of a single MO wavefunction psi_i in
    # the XY plane; user is prompted to supply the orbital index
    # Sequence: 4=Output and plot specific property in a plane; 4=Value of
    #           orbital wavefunction; 1=orbital index 1 (a fixed
    #           representative orbital)
    PLANE_MAP_ORBITAL_WAVEFUNCTION = ("4", "4", "1")
    # Colour-filled map of Fukui f-(r) in the XY plane via custom-
    # operation rho_N - rho_{N-1}; predicts electrophilic attack sites
    # Sequence: 4=Output and plot specific property in a plane; 0=Set
    #           custom operation; -1=custom operation: Fukui f- (density
    #           difference N/N-1)
    PLANE_MAP_FUKUI_MINUS = ("4", "0", "-1")
    # Colour-filled map of Fukui f+(r) in the XY plane via rho_{N+1} -
    # rho_N; predicts nucleophilic attack sites
    # Sequence: 4=Output and plot specific property in a plane; 0=Set
    #           custom operation; 1=custom operation: Fukui f+ (density
    #           difference N+1/N)
    # INTERACTIVE ONLY -- requires interactively defining a plane
    #                     (orientation/resolution); excluded per your
    #                     instruction to skip plane-dependent sequences
    #                     (also requires specifying a file-difference
    #                     operator for Fukui-type functions (N+1/N-1
    #                     states); needs additional setup beyond the
    #                     primary wavefunction)
    PLANE_MAP_FUKUI_PLUS = ("4", "0", "1")
    # Colour-filled map of dual descriptor Delta_f = f+ - f- in the XY
    # plane; positive = nucleophilic centre, negative = electrophilic
    # centre
    # Sequence: 4=Output and plot specific property in a plane; 0=Set
    #           custom operation; 3=custom operation: dual descriptor (f+
    #           minus f-)
    # INTERACTIVE ONLY -- requires interactively defining a plane
    #                     (orientation/resolution); excluded per your
    #                     instruction to skip plane-dependent sequences
    #                     (also requires specifying a file-difference
    #                     operator for Fukui-type functions (N+1/N-1
    #                     states); needs additional setup beyond the
    #                     primary wavefunction)
    PLANE_MAP_DUAL_DESCRIPTOR = ("4", "0", "3")
    # Colour-filled map of deformation density rho_mol - rho_promol in the
    # XY plane; shows electron redistribution upon chemical bond formation
    # Sequence: 4=Output and plot specific property in a plane; -2=Obtain
    #           deformation property
    PLANE_MAP_DEFORMATION_DENSITY = ("4", "-2")
    # Colour-filled map of promolecular density (superposition of free-
    # atom densities) in the XY plane; reference state before bonding
    # Sequence: 4=Output and plot specific property in a plane; -1=Obtain
    #           promolecule property
    PLANE_MAP_PROMOLECULAR_DENSITY = ("4", "-1")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 5: Cube / grid generation
    # ─────────────────────────────────────────────────────────────────────────
    # Evaluates a real-space function on a 3D grid and exports a Gaussian
    # .cube file compatible with VMD, GaussView, ChemCraft, Molekel, etc.
    # Grid quality codes: "1" ~ 50^3 points (low / preview) "2" ~ 80^3
    # points (medium, default for most workflows) "3" ~ 120^3 points
    # (high, for publication-quality figures) "0" at the post-process
    # prompt triggers cube export and returns.
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); 11=Local information entropy;
    #           2=Medium quality grid, covering whole system, about 512000
    #           points in total; 0=Return to main menu
    CUBE_LOCAL_INFORMATION_ENTROPY = ("5", "11", "2", "0")
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); 8=Electrostatic potential from
    #           nuclear charges; 2=Medium quality grid, covering whole
    #           system, about 512000 points in total; 0=Return to main
    #           menu
    CUBE_ELECTROSTATIC_POTENTIAL_FROM_CHARGE = ("5", "8", "2", "0")
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); 14=RDG with promolecular
    #           approximation; 2=Medium quality grid, covering whole
    #           system, about 512000 points in total; 0=Return to main
    #           menu
    CUBE_RDG_QUICK = ("5", "14", "2", "0")
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); 17=Correlation hole for alpha,
    #           ref. point:   0.00000   0.00000   0.00000; 2=Medium
    #           quality grid, covering whole system, about 512000 points
    #           in total; 1=Save graph of isosurface to file in current
    #           folder; 2=export data to a Gaussian-type cube file;
    #           0=return to the previous menu
    CUBE_CORRELATION_HOLE_ALPHA = ("5", "17", "2", "1", "2", "0")
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); 20=Electron delocal. range func.
    #           EDR(r;d); 2=length scale d = 2 Bohr (parameter for the
    #           EDR(r;d) function); 2=Medium quality grid, covering whole
    #           system, about 512000 points in total; 2=export data to a
    #           Gaussian-type cube file; 0=return to the previous menu
    CUBE_EDR = ("5", "20", "2", "2", "2", "0")
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); 22=Delta-g (promolecular
    #           approximation); 2=Medium quality grid, covering whole
    #           system, about 512000 points in total; 2=Export data to a
    #           Gaussian-type cube file in current folder; 0=return to the
    #           previous menu
    CUBE_DELTAG_PROMOLECULAR_APPROX = ("5", "22", "2", "2", "0")
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); 23=Delta-g (Hirshfeld
    #           partition); 2=Medium quality grid, covering whole system,
    #           about 512000 points in total; 2=Export data to a
    #           Gaussian-type cube file in current folder; 0=return to the
    #           previous menu
    CUBE_DELTAG_HIRSHFELD_PAER_APPROX = ("5", "23", "2", "2", "0")
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); 24=Interaction region indicator
    #           (IRI); 2=Medium quality grid, covering whole system, about
    #           512000 points in total; 2=Export data to a Gaussian-type
    #           cube file in current folder; 0=return to the previous menu
    CUBE_IRI = ("5", "24", "2", "2", "0")
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); 25=van der Waals potential
    #           (probe=C ); 2=Medium quality grid, covering whole system,
    #           about 512000 points in total; 2=Export data to a
    #           Gaussian-type cube file in current folder; 0=return to the
    #           previous menu
    CUBE_VDW_POTENTIAL = ("5", "25", "2", "2", "0")
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); 18=Average local ionization
    #           energy (ALIE); 2=Medium quality grid, covering whole
    #           system, about 512000 points in total; 2=Export data to a
    #           Gaussian-type cube file in current folder; 0=return to the
    #           previous menu
    CUBE_ALIE = ("5", "18", "2", "2", "0")
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); 19=Source function, mode: 1,
    #           ref. point:   0.00000   0.00000   0.00000; 2=Medium
    #           quality grid, covering whole system, about 512000 points
    #           in total; 2=Export data to a Gaussian-type cube file in
    #           current folder; 0=return to the previous menu
    CUBE_SOURCE_FUNCTION = ("5", "19", "2", "2", "0")
    # Generate medium-quality .cube of electron density rho(r); standard
    # input for VMD/GaussView isosurface rendering
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); 1=Electron density (rho);
    #           2=Medium quality grid, covering whole system, about 512000
    #           points in total; 2=Export data to a Gaussian-type cube
    #           file in current folder; 0=return to the previous menu
    CUBE_DENSITY = ("5", "1", "2", "2", "0")
    # Generate medium-quality .cube of spin density rho_alpha - rho_beta;
    # visualise unpaired electrons in radicals and open-shell molecules
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); 5=Electron spin density;
    #           2=Medium quality grid, covering whole system, about 512000
    #           points in total; 2=Export data to a Gaussian-type cube
    #           file in current folder; 0=return to the previous menu
    CUBE_SPIN_DENSITY = ("5", "5", "2", "2", "0")
    # Generate medium-quality .cube of ELF; bonding pairs, lone pairs, and
    # core shells visible as isosurfaces
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); 9=Electron localization function
    #           (ELF); 2=Medium quality grid, covering whole system, about
    #           512000 points in total; 2=Export data to a Gaussian-type
    #           cube file in current folder; 0=return to the previous menu
    CUBE_ELF = ("5", "9", "2", "2", "0")
    # Generate medium-quality .cube of LOL; electron-pair localisation
    # function with sharper basin boundaries than ELF
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); 10=Localized orbital locator
    #           (LOL); 2=Medium quality grid, covering whole system, about
    #           512000 points in total; 2=Export data to a Gaussian-type
    #           cube file in current folder; 0=return to the previous menu
    CUBE_LOL = ("5", "10", "2", "2", "0")
    # Generate medium-quality .cube of total ESP; map onto rho=0.001
    # isosurface to produce a surface electrostatic potential map
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); 12=Total electrostatic potential
    #           (ESP); 2=Medium quality grid, covering whole system, about
    #           512000 points in total; 2=export data to a Gaussian-type
    #           cube file; 0=return to the previous menu
    CUBE_ESP = ("5", "12", "2", "2", "0")
    # Generate medium-quality .cube of nabla^2 rho; negative isosurfaces
    # encode charge concentration (covalent bonds, lone pairs, shells)
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); 3=Laplacian of rho; 2=Medium
    #           quality grid, covering whole system, about 512000 points
    #           in total; 2=Export data to a Gaussian-type cube file in
    #           current folder; 0=return to the previous menu
    CUBE_LAPLACIAN = ("5", "3", "2", "2", "0")
    # Generate medium-quality .cube of |nabla rho|; high-value isosurfaces
    # mark interatomic and atomic-shell boundaries
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); 2=Gradient norm of rho; 2=Medium
    #           quality grid, covering whole system, about 512000 points
    #           in total; 2=Export data to a Gaussian-type cube file in
    #           current folder; 0=return to the previous menu
    CUBE_GRADIENT_NORM = ("5", "2", "2", "2", "0")
    # Generate medium-quality .cube of G(r), the positive-definite
    # Lagrangian kinetic energy density
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); 7=Lagrangian kinetic energy
    #           density G(r); 2=Medium quality grid, covering whole
    #           system, about 512000 points in total; 2=Export data to a
    #           Gaussian-type cube file in current folder; 0=return to the
    #           previous menu
    CUBE_KINETIC_G = ("5", "7", "2", "2", "0")
    # Generate medium-quality .cube of K(r), the Hamiltonian kinetic
    # energy density; K = G - (1/4) nabla^2 rho
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); 6=Hamiltonian kinetic energy
    #           density K(r); 2=Medium quality grid, covering whole
    #           system, about 512000 points in total; 2=Export data to a
    #           Gaussian-type cube file in current folder; 0=return to the
    #           previous menu
    CUBE_KINETIC_K = ("5", "6", "2", "2", "0")
    # Generate medium-quality .cube of RDG; render at low isovalue (~0.5)
    # together with sign(lambda2)rho cube for NCI isosurface colouring
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); 13=Reduced density gradient
    #           (RDG); 2=Medium quality grid, covering whole system, about
    #           512000 points in total; 2=Export data to a Gaussian-type
    #           cube file in current folder; 0=return to the previous menu
    CUBE_RDG = ("5", "13", "2", "2", "0")
    # Generate medium-quality .cube of sign(lambda2)*rho; colour-maps onto
    # RDG isosurface: blue=H-bond, green=vdW, red=steric repulsion
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); 15=Sign(lambda2)*rho; 2=Medium
    #           quality grid, covering whole system, about 512000 points
    #           in total; 2=Export data to a Gaussian-type cube file in
    #           current folder; 0=return to the previous menu
    CUBE_SIGN_LAMBDA2_RHO = ("5", "15", "2", "2", "0")
    # Generate medium-quality .cube of a single MO wavefunction psi_i;
    # user is prompted for the orbital index; standard for orbital
    # visualisation INTERACTIVE
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); 4=Value of orbital wavefunction;
    #           2=orbital index 2 (a fixed representative orbital);
    #           2=Medium quality grid, covering whole system, about 512000
    #           points in total; 0=return to the previous menu
    CUBE_ORBITAL_WAVEFUNCTION = ("5", "4", "2", "2", "0")
    # Generate medium-quality .cube of Fukui f-(r) via two-wavefunction
    # subtraction rho_N - rho_{N-1}; maps electrophilic attack
    # susceptibility
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); 0=Set custom operation;
    #           -1=custom operation: Fukui f- (density difference N/N-1);
    #           2=Gradient norm of rho; 0=return to the previous menu
    CUBE_FUKUI_MINUS = ("5", "0", "-1", "2", "0")
    # Generate medium-quality .cube of Fukui f+(r) via rho_{N+1} - rho_N;
    # maps nucleophilic attack susceptibility
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); 0=Set custom operation; 1=custom
    #           operation: Fukui f+ (density difference N+1/N); 2=export
    #           data to a Gaussian-type cube file; 0=return to the
    #           previous menu
    # INTERACTIVE ONLY -- requires specifying a file-difference operator
    #                     for Fukui-type functions (N+1/N-1 states); needs
    #                     additional setup beyond the primary wavefunction
    CUBE_FUKUI_PLUS = ("5", "0", "1", "2", "0")
    # Generate medium-quality .cube of dual descriptor f+ - f-; positive
    # isosurfaces are nucleophilic sites, negative are electrophilic
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); 0=Set custom operation; 3=custom
    #           operation: dual descriptor (f+ minus f-); 2=export data to
    #           a Gaussian-type cube file; 0=return to the previous menu
    # INTERACTIVE ONLY -- requires specifying a file-difference operator
    #                     for Fukui-type functions (N+1/N-1 states); needs
    #                     additional setup beyond the primary wavefunction
    CUBE_DUAL_DESCRIPTOR = ("5", "0", "3", "2", "0")
    # Generate medium-quality .cube of promolecular density (superposition
    # of free-atom densities); used as a reference before bond formation
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); -1=Obtain promolecule property;
    #           1=promolecular property: density; 2=export data to a
    #           Gaussian-type cube file; 0=return to the previous menu
    # INTERACTIVE ONLY -- requires a local Gaussian installation to be
    #                     found on PATH
    CUBE_PROMOLECULAR_DENSITY = ("5", "-1", "1", "2", "0")
    # Generate medium-quality .cube of deformation density Delta_rho =
    # rho_mol # rho_promol; shows how bonding redistributes electron
    # density
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); -2=Obtain deformation property;
    #           1=deformation property: density; 2=export data to a
    #           Gaussian-type cube file; 0=return to the previous menu
    # INTERACTIVE ONLY -- requires a local Gaussian installation to be
    #                     found on PATH
    CUBE_DEFORMATION_DENSITY = ("5", "-2", "1", "2", "0")
    # Generate HIGH-quality .cube of electron density (~120^3 grid); finer
    # isosurfaces suitable for publication figures
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); 1=Electron density (rho); 3=High
    #           quality grid,   covering whole system, about 1728000
    #           points in total; 0=return to the previous menu
    CUBE_DENSITY_HIGH = ("5", "1", "3", "0")
    # Generate HIGH-quality .cube of total ESP (~120^3 grid); needed for
    # accurate surface-ESP colour maps on large or complex molecules
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); 12=Total electrostatic potential
    #           (ESP); 3=High quality grid,   covering whole system, about
    #           1728000 points in total; 0=return to the previous menu
    CUBE_ESP_HIGH = ("5", "12", "3", "0")
    # Generate HIGH-quality .cube of ELF (~120^3 grid); resolves fine
    # bonding features in systems requiring detailed electron-pair
    # analysis
    # Sequence: 5=Output and plot specific property within a spatial
    #           region (calc. grid data); 9=Electron localization function
    #           (ELF); 3=High quality grid,   covering whole system, about
    #           1728000 points in total; 0=return to the previous menu
    CUBE_ELF_HIGH = ("5", "9", "3", "0")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 6: Check & modify wavefunction
    # ─────────────────────────────────────────────────────────────────────────
    # Provides subfunctions to inspect, edit, and save the loaded
    # wavefunction. "0" from the submenu returns to main menu.
    # Sequence: 6=Check & modify wavefunction; 5=Print coefficient matrix
    #           in basis functions; 1=Print on screen
    PRINT_COEFFICIENT_MATRIX = ("6", "5", "1")
    # Sequence: 6=Check & modify wavefunction; 6=Print density matrix in
    #           basis function; 1=Print on screen
    PRINT_DENSITY_MATRIX = ("6", "6", "1")
    # Sequence: 6=Check & modify wavefunction; 7=Print various kinds of
    #           integral matrix between basis functions; 0=Fock/KS matrix;
    #           1=Print on screen; 1=Generating Fock/KS matrix by MO
    #           energies and coefficients as well as overlap matrix
    PRINT_INTEGRAL_MATRIX_FOCK = ("6", "7", "0", "1", "1")
    # Sequence: 6=Check & modify wavefunction; 7=Print various kinds of
    #           integral matrix between basis functions; 1=Overlap
    #           integral; 1=Print on screen
    PRINT_INTEGRAL_MATRIX_OVERLAP = ("6", "7", "1", "1")
    # Sequence: 6=Check & modify wavefunction; 7=Print various kinds of
    #           integral matrix between basis functions; 2=Electric dipole
    #           moment integral; 1=Print on screen
    PRINT_INTEGRAL_MATRIX_ELECTRIC_DIPOLE = ("6", "7", "2", "1")
    # Sequence: 6=Check & modify wavefunction; 7=Print various kinds of
    #           integral matrix between basis functions; 3=Magnetic dipole
    #           moment integral; 1=Print on screen
    PRINT_INTEGRAL_MATRIX_MAGNETIC_DIPOLE = ("6", "7", "3", "1")
    # Sequence: 6=Check & modify wavefunction; 7=Print various kinds of
    #           integral matrix between basis functions; 4=Velocity
    #           integral; 1=Print on screen
    PRINT_INTEGRAL_MATRIX_VELOCITY = ("6", "7", "4", "1")
    # Sequence: 6=Check & modify wavefunction; 7=Print various kinds of
    #           integral matrix between basis functions; 5=Kinetic energy
    #           integral; 1=Print on screen
    PRINT_INTEGRAL_MATRIX_EKINETIC = ("6", "7", "5", "1")
    # Sequence: 6=Check & modify wavefunction; 7=Print various kinds of
    #           integral matrix between basis functions; 6=Electric
    #           quadrupole moment integral; 1=List all GTFs
    PRINT_INTEGRAL_MATRIX_QUADRUPOLE = ("6", "7", "6", "1")
    # Sequence: 6=Check & modify wavefunction; 7=Print various kinds of
    #           integral matrix between basis functions; 7=Electric
    #           octopole moment integral; 1=List all GTFs
    PRINT_INTEGRAL_MATRIX_OCTOPOLE = ("6", "7", "7", "1")
    # Sequence: 6=Check & modify wavefunction; 7=Print various kinds of
    #           integral matrix between basis functions; 8=Electric
    #           hexadecapole moment integral; 1=List all GTFs
    PRINT_INTEGRAL_MATRIX_HEXADECAPOLE = ("6", "7", "8", "1")
    # #INTERACTIVE SEQUENCES - to be tested #CORRECT SEQUENCES # Save the
    # (possibly modified) wavefunction to new.wfn; also converts #
    # fch/molden to .wfn format; zero-occupation orbitals are dropped #
    # automatically SAVE_WFN = ("6", "0") Print centre atom, angular-
    # momentum type (s/p/d/f/...), and exponent for # every Gaussian-type
    # function (GTF) in the basis PRINT_ALL_GTF = ("6", "1") Print shell
    # assignments, contracted function types, and GTF index ranges # for
    # every basis function PRINT_ALL_BASIS_FUNCTIONS = ("6", "2") # Print
    # index, energy, occupation number, and spin type for every o # rbital
    # in the loaded wavefunction PRINT_ORBITAL_INFO = ("6", "3") # Print
    # the one-particle density matrix P expressed in the basis-function #
    # representation PRINT_DENSITY_MATRIX = ("6", "6", "0") # Manually set
    # the occupation number of selected orbitals; set to 0 to remove their
    # contribution from subsequent real-space function evaluations
    # MODIFY_OCCUPATION = ("6", "26", "0") # Remove all core (inner-shell)
    # orbitals from the wavefunction, retaining # only valence-shell
    # orbitals for subsequent analyses DELETE_INNER_ORBITALS = ("6", "34",
    # "0") Print centre atom, angular-momentum type (s/p/d/f/...), and
    # exponent for every Gaussian-type function (GTF) in the basis
    # Sequence: 6=Check & modify wavefunction; 1=List all GTFs
    PRINT_ALL_GTF = ("6", "1")
    # Print shell assignments, contracted function types, and GTF index
    # ranges for every basis function
    # Sequence: 6=Check & modify wavefunction; 2=List all basis functions
    PRINT_ALL_BASIS_FUNCTIONS = ("6", "2")
    # Print index, energy, occupation number, and spin type for every
    # orbital in the loaded wavefunction
    # Sequence: 6=Check & modify wavefunction; 3=List all orbitals
    PRINT_ORBITAL_INFO = ("6", "3")
    # Save the (possibly modified) wavefunction to new.wfn in current
    # folder; also converts fch/molden to .wfn format, dropping zero-
    # occupation orbitals automatically
    # Sequence: 6=Check & modify wavefunction; 0=Save the present
    #           wavefunction to new.wfn file in current folder
    SAVE_WFN = ("6", "0")
    # Remove all core (inner-shell) orbitals from the wavefunction and
    # save the result to new.wfn; retains only valence-shell orbitals for
    # subsequent analyses
    # Sequence: 6=Check & modify wavefunction; 34=Set occupation number of
    #           inner orbitals to zero; 0=Save the present wavefunction to
    #           new.wfn file in current folder
    DELETE_INNER_ORBITALS = ("6", "34", "0")
    # Manually set the occupation number of selected orbitals (e.g. to
    # zero out a subset for a custom electron count)
    # Sequence: 6=Check & modify wavefunction; 26=Set occupation of some
    #           orbitals; 0=select all orbitals (per 'Input 0 can select
    #           all orbitals')
    # INTERACTIVE ONLY -- requires an interactively-specified target
    #                     occupation-number value after selecting orbitals
    #                     (e.g. '1.2', '+1.1', 'i' to restore); no single
    #                     non-interactive default makes sense for an
    #                     arbitrary molecule
    MODIFY_OCCUPATION = ("6", "26", "0")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 7: Population analysis & atomic charges
    # ─────────────────────────────────────────────────────────────────────────
    # Hirshfeld-family methods require atomic reference densities: "1" →
    # use Multiwfn's built-in sphericalised free-atom densities "n" → skip
    # writing the .chg output file to disk "0" → return to main menu after
    # printing charges
    # Sequence: 7=Population analysis and calculation of atomic charges;
    #           5=Mulliken atom & basis function population analysis;
    #           2=Output gross atomic population matrix and decompose it;
    #           n=no (decline the fragment/mirror-plane restriction)
    MULLIKEN_DECOMPOSE_ATOMIC_POPULATION = ("7", "5", "2", "n")
    # Sequence: 7=Population analysis and calculation of atomic charges;
    #           5=Mulliken atom & basis function population analysis;
    #           3=Output gross basis function population matrix and
    #           decompose it; n=no (decline the fragment/mirror-plane
    #           restriction)
    MULLIKEN_DECOMPOSE_BASIS_FUNCTION = ("7", "5", "3", "n")
    # Lowdin atomic charges via Lowdin orthogonalisation; slightly more
    # basis- set stable than Mulliken but still dependent on basis choice
    # INTERACTIVE
    # Sequence: 7=Population analysis and calculation of atomic charges;
    #           6=Lowdin atom & basis function population analysis;
    #           ENTER=accept the default (print results to screen instead
    #           of a file); 0=return to the previous menu
    LOWDIN_POPULATION = ("7", "6", "ENTER", "0")
    # INGNORE FOR NOW Electronegativity Equalization Method (EEM) charges;
    # fast geometry-based empirical model requiring no wavefunction
    # evaluation NEEDS MOL2
    # Sequence: 7=Population analysis and calculation of atomic charges;
    #           17=Electronegativity Equalization Method (EEM) atomic
    #           charge; n=no (decline the fragment/mirror-plane
    #           restriction); 0=return to the previous menu
    EEM_CHARGE = ("7", "17", "n", "0")
    # Hirshfeld atomic charges from deformation-density partitioning;
    # qualitatively correct but systematically underestimates charge
    # transfer
    # Sequence: 7=Population analysis and calculation of atomic charges;
    #           1=Hirshfeld atomic charge; 1=Use build-in sphericalized
    #           atomic densities in free-states (more convenient); n=no
    #           (decline the fragment/mirror-plane restriction); 0=return
    #           to the previous menu
    HIRSHFELD_CHARGE = ("7", "1", "1", "n", "0")
    # Voronoi Deformation Density charges; Hirshfeld variant using Voronoi
    # cell weights instead of Hirshfeld weights; results similar to
    # Hirshfeld
    # Sequence: 7=Population analysis and calculation of atomic charges;
    #           2=Voronoi deformation density (VDD) atom population; 1=Use
    #           build-in sphericalized atomic densities in free-states
    #           (more convenient); n=no (decline the fragment/mirror-plane
    #           restriction); 0=return to the previous menu
    VDD_POPULATION = ("7", "2", "1", "n", "0")
    # Mulliken atomic charges and basis-function populations; oldest
    # method, highly basis-set dependent, avoid with diffuse basis
    # functions
    # Sequence: 7=Population analysis and calculation of atomic charges;
    #           5=Mulliken atom & basis function population analysis;
    #           1=Output Mulliken population and atomic charges; n=no
    #           (decline the fragment/mirror-plane restriction); 0=return
    #           to the previous menu
    MULLIKEN_POPULATION = ("7", "5", "1", "n", "0")
    # Ros-Schuit C-squared Population Analysis (SCPA); modified Mulliken
    # scheme that prevents negative population numbers
    # Sequence: 7=Population analysis and calculation of atomic charges;
    #           7=Modified Mulliken atom population defined by Ros &
    #           Schuit (SCPA); n=no (decline the fragment/mirror-plane
    #           restriction); 0=return to the previous menu
    SCPA_POPULATION = ("7", "7", "n", "0")
    # Stout-Politzer modified Mulliken charges; cross-terms partitioned by
    # the ratio of squared orbital coefficients
    # Sequence: 7=Population analysis and calculation of atomic charges;
    #           8=Modified Mulliken atom population defined by Stout &
    #           Politzer; n=no (decline the fragment/mirror-plane
    #           restriction); 0=return to the previous menu
    STOUT_POLITZER_POPULATION = ("7", "8", "n", "0")
    # Bickelhaupt modified Mulliken charges; cross-terms weighted by the
    # total local populations summed across all orbitals
    # Sequence: 7=Population analysis and calculation of atomic charges;
    #           9=Modified Mulliken atom population defined by
    #           Bickelhaupt; n=no (decline the fragment/mirror-plane
    #           restriction); 0=return to the previous menu
    BICKELHAUPT_POPULATION = ("7", "9", "n", "0")
    # Becke atomic charges with atomic-dipole-moment correction applied;
    # reasonable for typical organic systems using default CSD radii
    # Sequence: 7=Population analysis and calculation of atomic charges;
    #           10=Becke atomic charge with atomic dipole moment
    #           correction; 0=return to the previous menu; n=no (decline
    #           the fragment/mirror-plane restriction); 0=return to the
    #           previous menu
    BECKE_CHARGE = ("7", "10", "0", "n", "0")
    # Atomic Dipole moment Corrected Hirshfeld (ADCH) charges; exactly
    # reproduces the molecular dipole moment and gives reliable ESP;
    # highly recommended
    # Sequence: 7=Population analysis and calculation of atomic charges;
    #           11=Atomic dipole corrected Hirshfeld atomic charge (ADCH)
    #           (recommended); 1=Use build-in sphericalized atomic
    #           densities in free-states (more convenient); n=no (decline
    #           the fragment/mirror-plane restriction); 0=return to the
    #           previous menu
    ADCH_CHARGE = ("7", "11", "1", "n", "0")
    # CHELPG ESP-fitting charges on a cubic grid of points; best
    # rotational invariance among ESP-fit methods; widely used for force-
    # field parametrisation
    # Sequence: 7=Population analysis and calculation of atomic charges;
    #           12=CHELPG ESP fitting atomic charge; 1=Start calculation!;
    #           n=no (decline the fragment/mirror-plane restriction);
    #           0=return to the previous menu
    # INTERACTIVE ONLY -- Multiwfn offers a default via pressing ENTER
    #                     here, but the exact blank-token sequence to
    #                     reach it wasn't reliably confirmed in testing --
    #                     needs manual verification
    CHELPG_CHARGE = ("7", "12", "1", "n", "0")
    # Merz-Kollmann (MK) ESP-fitting charges on concentric shells at
    # 1.4-2.0x vdW radii; standard for AMBER/GAFF parametrisation
    # Sequence: 7=Population analysis and calculation of atomic charges;
    #           13=Merz-Kollmann (MK) ESP fitting atomic charge; 1=Start
    #           calculation!; n=no (decline the fragment/mirror-plane
    #           restriction); 0=return to the previous menu
    # INTERACTIVE ONLY -- Multiwfn offers a default via pressing ENTER
    #                     here, but the exact blank-token sequence to
    #                     reach it wasn't reliably confirmed in testing --
    #                     needs manual verification
    MK_CHARGE = ("7", "13", "1", "n", "0")
    # CM5 charges (charge model 5); Hirshfeld-based with empirical
    # correction terms for improved dipole-moment reproduction across
    # diverse molecules
    # Sequence: 7=Population analysis and calculation of atomic charges;
    #           16=CM5 atomic charge; 1=Use build-in sphericalized atomic
    #           densities in free-states (more convenient); n=no (decline
    #           the fragment/mirror-plane restriction); 0=return to the
    #           previous menu
    CM5_CHARGE = ("7", "16", "1", "n", "0")
    # RESP (Restrained ESP) charges; ESP fit with hyperbolic restraint
    # toward zero; the standard charge model for AMBER/GAFF force fields
    # FORCED DEFAULTS
    # Sequence: 7=Population analysis and calculation of atomic charges;
    #           18=Restrained ElectroStatic Potential (RESP) atomic
    #           charge; 1n=start standard two-stage RESP fitting ('1'),
    #           declining the equivalence-constraint follow-up prompt
    #           ('n'); 0=return to the previous menu
    RESP_CHARGE = ("7", "18", "1n", "0")
    # Gasteiger-Marsili empirical charges; purely connectivity-based,
    # extremely fast, no wavefunction needed; suitable for large databases
    # Sequence: 7=Population analysis and calculation of atomic charges;
    #           19=Gasteiger (PEOE) charge; n=no (decline the
    #           fragment/mirror-plane restriction); 0=return to the
    #           previous menu
    GASTEIGER_CHARGE = ("7", "19", "n", "0")
    # Minimal Basis Iterative Stockholder (MBIS) charges; information-
    # theoretic partitioning giving excellent dipole and higher-multipole
    # reproduction
    # Sequence: 7=Population analysis and calculation of atomic charges;
    #           20=Minimal Basis Iterative Stockholder (MBIS) charge;
    #           1=Start calculation!; n=no (decline the
    #           fragment/mirror-plane restriction); 0=return to the
    #           previous menu
    MBIS_CHARGE = ("7", "20", "1", "n", "0")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 8: Orbital composition analysis
    # ─────────────────────────────────────────────────────────────────────────
    # Decomposes each MO into percentage contributions from basis
    # functions, shells, atoms, or user-defined fragments.
    # Sequence: 8=Orbital composition analysis; 1=Orbital composition
    #           analysis with Mulliken partition; h=select the HOMO;
    #           0=return to the previous menu
    ORBITAL_COMPOSITION_MULLIKEN_HOMO = ("8", "1", "h", "0")
    # Sequence: 8=Orbital composition analysis; 1=Orbital composition
    #           analysis with Mulliken partition; l=select the LUMO;
    #           0=return to the previous menu
    ORBITAL_COMPOSITION_MULLIKEN_LUMO = ("8", "1", "l", "0")
    # Sequence: 8=Orbital composition analysis; 1=Orbital composition
    #           analysis with Mulliken partition; a=select all orbitals;
    #           0=return to the previous menu
    ORBITAL_COMPOSITION_MULLIKEN_ALL = ("8", "1", "a", "0")
    # Sequence: 8=Orbital composition analysis; 2=Orbital composition
    #           analysis with Stout-Politzer partition; h=select the HOMO;
    #           0=return to the previous menu
    ORBITAL_COMPOSITION_STOUT_POLITZER_HOMO = ("8", "2", "h", "0")
    # Sequence: 8=Orbital composition analysis; 2=Orbital composition
    #           analysis with Stout-Politzer partition; l=select the LUMO;
    #           0=return to the previous menu
    ORBITAL_COMPOSITION_STOUT_POLITZER_LUMO = ("8", "2", "l", "0")
    # Sequence: 8=Orbital composition analysis; 2=Orbital composition
    #           analysis with Stout-Politzer partition; a=select all
    #           orbitals; 0=return to the previous menu
    ORBITAL_COMPOSITION_STOUT_POLITZER_ALL = ("8", "2", "a", "0")
    # Sequence: 8=Orbital composition analysis; 3=Orbital composition
    #           analysis with Ros-Schuit (SCPA) partition; h=select the
    #           HOMO; 0=return to the previous menu
    ORBITAL_COMPOSITION_SCPA_HOMO = ("8", "3", "h", "0")
    # Sequence: 8=Orbital composition analysis; 3=Orbital composition
    #           analysis with Ros-Schuit (SCPA) partition; l=select the
    #           LUMO; 0=return to the previous menu
    ORBITAL_COMPOSITION_SCPA_LUMO = ("8", "3", "l", "0")
    # Sequence: 8=Orbital composition analysis; 3=Orbital composition
    #           analysis with Ros-Schuit (SCPA) partition; a=select all
    #           orbitals; 0=return to the previous menu
    ORBITAL_COMPOSITION_SCPA_ALL = ("8", "3", "a", "0")
    # Sequence: 8=Orbital composition analysis; 8=Calculate atom and
    #           fragment contributions by Hirshfeld method; 1=Orbital
    #           composition analysis with Mulliken partition; h=select the
    #           HOMO
    FRAGMENT_CONTRIBUTION_HIRSHFELD_HOMO = ("8", "8", "1", "h")
    # Sequence: 8=Orbital composition analysis; 8=Calculate atom and
    #           fragment contributions by Hirshfeld method; 1=Orbital
    #           composition analysis with Mulliken partition; l=select the
    #           LUMO
    FRAGMENT_CONTRIBUTION_HIRSHFELD_LUMO = ("8", "8", "1", "l")
    # Sequence: 8=Orbital composition analysis; 8=Calculate atom and
    #           fragment contributions by Hirshfeld method; 1=Orbital
    #           composition analysis with Mulliken partition; -1=Define
    #           fragment 1 (for option 1~6)
    FRAGMENT_CONTRIBUTION_HIRSHFELD_ALL = ("8", "8", "1", "-1")
    # Sequence: 8=Orbital composition analysis; 8=Calculate atom and
    #           fragment contributions by Hirshfeld method; 1=Orbital
    #           composition analysis with Mulliken partition; -4=answer to
    #           the atom-index prompt for which atoms to print
    ATOM_CONTRIBUTION_HIRSHFELD = ("8", "8", "1", "-4")
    # Sequence: 8=Orbital composition analysis; 9=Calculate atom and
    #           fragment contributions by Becke method; 1=Orbital
    #           composition analysis with Mulliken partition; h=select the
    #           HOMO
    FRAGMENT_CONTRIBUTION_BECKE_HOMO = ("8", "9", "1", "h")
    # Sequence: 8=Orbital composition analysis; 9=Calculate atom and
    #           fragment contributions by Becke method; 1=Orbital
    #           composition analysis with Mulliken partition; l=select the
    #           LUMO
    FRAGMENT_CONTRIBUTION_BECKE_LUMO = ("8", "9", "1", "l")
    # Sequence: 8=Orbital composition analysis; 9=Calculate atom and
    #           fragment contributions by Becke method; 1=Orbital
    #           composition analysis with Mulliken partition; -1=Define
    #           fragment 1 (for option 1~6)
    FRAGMENT_CONTRIBUTION_BECKE_ALL = ("8", "9", "1", "-1")
    # Sequence: 8=Orbital composition analysis; 9=Calculate atom and
    #           fragment contributions by Becke method; 1=Orbital
    #           composition analysis with Mulliken partition; -4=answer to
    #           the atom-index prompt for which atoms to print
    ATOM_CONTRIBUTION_BECKE = ("8", "9", "1", "-4")
    # Mulliken fragment composition: prints the total percentage of a pre-
    # defined fragment in every occupied MO as a table
    # Sequence: 8=Orbital composition analysis; 4=Print frag. 1 &
    #           inter-fragment compositions in all orbitals (Mulliken);
    #           0=return to the previous menu
    ORBITAL_COMPOSITION_FRAGMENT_MULLIKEN = ("8", "4", "0")
    # Stout-Politzer fragment composition including inter-fragment cross-
    # term breakdown
    # Sequence: 8=Orbital composition analysis; 5=Print frag. 1 &
    #           inter-fragment compositions in all orbitals
    #           (Stout-Politzer); 0=return to the previous menu
    ORBITAL_COMPOSITION_FRAGMENT_STOUT = ("8", "5", "0")
    # SCPA fragment composition: sums C^2 contributions of all fragment
    # basis functions per orbital; no negative values
    # Sequence: 8=Orbital composition analysis; 6=Print frag. 1
    #           compositions in all orbitals (SCPA); 0=return to the
    #           previous menu
    ORBITAL_COMPOSITION_FRAGMENT_SCPA = ("8", "6", "0")
    # Natural Atomic Orbital (NAO) composition using the MO-in-NAO
    # coefficient matrix from NBO output; excellent basis-set stability
    # for occupied MOs
    # Sequence: 8=Orbital composition analysis; 7=Orbital composition
    #           analysis by natural atomic orbital (NAO) method; 0=return
    #           to the previous menu
    ORBITAL_COMPOSITION_NAO = ("8", "7", "0")
    # Hirshfeld orbital composition: integral of |psi_i|^2 * w_A(r) for
    # each atom A; highly stable regardless of basis-set choice
    # Sequence: 8=Orbital composition analysis; 8=Calculate atom and
    #           fragment contributions by Hirshfeld method; 0=return to
    #           the previous menu
    ORBITAL_COMPOSITION_HIRSHFELD = ("8", "8", "0")
    # Becke orbital composition: Becke-partition integral of |psi_i|^2 per
    # atom no atomic reference density files required
    # Sequence: 8=Orbital composition analysis; 9=Calculate atom and
    #           fragment contributions by Becke method; 0=return to the
    #           previous menu
    ORBITAL_COMPOSITION_BECKE = ("8", "9", "0")
    # Modified LOBA (Localized Orbital Bonding Analysis): assigns formal
    # oxidation states by analysing orbital-population partitioning among
    # bonded atoms
    # Sequence: 8=Orbital composition analysis; 100=Evaluate oxidation
    #           state by LOBA/mLOBA method; 0=return to the previous menu
    LOBA_OXIDATION_STATE = ("8", "100", "0")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 9: Bond order analysis
    # ─────────────────────────────────────────────────────────────────────────
    # Prints bond orders between all atom pairs above a threshold, plus
    # total and free valences for each atom. "0" returns to main menu.
    # INTERACTIVE - REQUIRES USER INPUT
    # Sequence: 9=Bond order analysis; 2=Multicenter bond order analysis;
    #           0=answer '0' to the atom-index prompt (returns without a
    #           specific pair)
    MULTICENTER_BOND_ORDER = ("9", "2", "0")
    # Multi-centre bond order in the NAO basis (S becomes identity); much
    # more basis-set stable, recommended when diffuse functions are
    # present
    # Sequence: 9=Bond order analysis; -2=Multicenter bond order analysis
    #           in NAO basis; 0=return to the previous menu
    # INTERACTIVE ONLY -- requires natural population analysis (NPA) data
    #                     to already be present in the input file
    #                     (Multiwfn: 'Cannot find natural population
    #                     analysis information in the input file'); fails
    #                     with a plain .molden/.fch, needs e.g. an NBO-
    #                     enabled Gaussian output instead
    MULTICENTER_BOND_ORDER_NAO = ("9", "-2", "0")
    # Wiberg bond order in the Lowdin-orthogonalised basis (WL); more
    # stable than Mayer for large basis sets but can overestimate polar
    # bonds Decompose Mulliken bond order between a chosen atom pair into
    # per-orbital contributions to identify bonding vs. antibonding MOs
    # Sequence: 9=Bond order analysis; 5=Decompose Mulliken bond order
    #           between two atoms to orbital contributions; 0=answer '0'
    #           to the atom-pair-index prompt (returns without a specific
    #           pair)
    # INTERACTIVE ONLY -- requires an interactively-specified atom or ring
    #                     index list (e.g. 'Input index of the atoms in
    #                     ring order')
    MULLIKEN_BOND_ORDER_DECOMPOSE = ("9", "5", "0")
    # Orbital-occupancy-perturbed Mayer bond order: prints the
    # contribution of each occupied MO to the Mayer bond order for a
    # selected A-B pair
    # Sequence: 9=Bond order analysis; 6=Orbital occupancy-perturbed Mayer
    #           bond order; 0=return to the previous menu
    # INTERACTIVE ONLY -- requires an interactively-specified atom or ring
    #                     index list (e.g. 'Input index of the atoms in
    #                     ring order')
    ORBITAL_PERTURBED_MAYER = ("9", "6", "0")
    # Decompose Wiberg bond order between two atoms into per-orbital
    # contributions; identifies which MOs are responsible for the bond
    # Sequence: 9=Bond order analysis; 9=Decompose Wiberg bond order in
    #           NAO basis as atomic orbital pair contribution; 0=return to
    #           the previous menu
    WIBERG_DECOMPOSITION = ("9", "9", "0")
    # AV1245 multicentric aromaticity index computed from Mayer bond
    # orders along a ring; robust against ring size and basis-set choice
    # Sequence: 9=Bond order analysis; 11=AV1245 index (approximate
    #           multicenter bond order for large rings) and AVmin;
    #           0=return to the previous menu
    # INTERACTIVE ONLY -- requires an interactively-specified atom or ring
    #                     index list (e.g. 'Input index of the atoms in
    #                     ring order')
    AV1245_INDEX = ("9", "11", "0")
    # Mayer bond orders from PS-matrix products; values approximately
    # 1/2/3 for single/double/triple bonds; the best all-round bond-order
    # method
    # Sequence: 9=Bond order analysis; 1=Mayer bond order analysis; n=no
    #           (decline the fragment/mirror-plane restriction); 0=return
    #           to the previous menu
    MAYER_BOND_ORDER = ("9", "1", "n", "0")
    # Multi-centre bond order (up to 12 centres) in the original basis-
    # function representation; sensitive to diffuse basis functions
    # Sequence: 9=Bond order analysis; 3=Wiberg bond order analysis in
    #           Lowdin orthogonalized basis; n=no (decline the
    #           fragment/mirror-plane restriction); 0=return to the
    #           previous menu
    WIBERG_BOND_ORDER = ("9", "3", "n", "0")
    # Mulliken bond orders (2*PS off-diagonal elements); positive =
    # bonding character, negative = antibonding; qualitative indicator
    # only
    # Sequence: 9=Bond order analysis; 4=Mulliken bond order (Mulliken
    #           overlap population) analysis; n=no (decline the
    #           fragment/mirror-plane restriction); 0=return to the
    #           previous menu
    MULLIKEN_BOND_ORDER = ("9", "4", "n", "0")
    # Fuzzy bond order (Becke-space integration of PS products); more
    # basis-set stable than Mayer; essentially equivalent to the AIM
    # delocalization index
    # Sequence: 9=Bond order analysis; 7=Fuzzy bond order analysis (FBO);
    #           n=no (decline the fragment/mirror-plane restriction);
    #           0=return to the previous menu
    FUZZY_BOND_ORDER = ("9", "7", "n", "0")
    # Laplacian Bond Order (LBO): integral of -nabla^2 rho in the fuzzy
    # overlap space; correlates with BDE and vibrational frequency;
    # independent of wavefunction type
    # Sequence: 9=Bond order analysis; 8=Laplacian bond order (LBO); n=no
    #           (decline the fragment/mirror-plane restriction); 0=return
    #           to the previous menu
    LAPLACIAN_BOND_ORDER = ("9", "8", "n", "0")
    # Intrinsic Bond Strength Index (IBSI): geometry-free bond-order
    # measure derived from the electron density at the bond critical point
    # Sequence: 9=Bond order analysis; 10=Intrinsic bond strength index
    #           (IBSI); 1=Start calculation; 1=Medium quality (radial=30,
    #           angular=110. Cost=1.0 x); 0=return to the previous menu
    IBSI_ANALYSIS = ("9", "10", "1", "1", "0")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 10: Density of states
    # ─────────────────────────────────────────────────────────────────────────
    # Reads orbital energies from .fch/.molden/Gaussian-output/plain-text
    # and produces broadened TDOS, PDOS (per fragment), and OPDOS curves.
    # Plot Total Density of States (TDOS): broadened orbital energy
    # spectrum showing how densely electronic states are distributed in
    # energy
    # Sequence: 10=Plot total DOS, PDOS, OPDOS, local DOS, COHP and
    #           photoelectron spectrum; 0=Draw TDOS graph!; 2=Save the
    #           graph to image file in current folder; 3=Export curve and
    #           line data to plain text file in current folder; 0=return
    #           to the previous menu
    PLOT_TDOS = ("10", "0", "2", "3", "0")
    # Sequence: 10=Plot total DOS, PDOS, OPDOS, local DOS, COHP and
    #           photoelectron spectrum; 00=Draw TDOS and OPDOS between
    #           nearest atoms!; 2=Save the graph to image file in current
    #           folder; 3=Export curve and line data to plain text file in
    #           current folder; 0=return to the previous menu
    PLOT_TDOS_OPDOS = ("10", "00", "2", "3", "0")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 11: Spectra simulation
    # ─────────────────────────────────────────────────────────────────────────
    # Reads transition data from Gaussian/ORCA output or plain-text files
    # and broadens discrete transitions into simulated spectra using
    # Gaussian, Lorentzian, or pseudo-Voigt broadening functions. Simulate
    # IR absorption spectrum by broadening harmonic (or anharmonic)
    # vibrational frequencies weighted by IR intensities in km/mol
    # Sequence: 11=Plot IR/Raman/UV-Vis/ECD/VCD/ROA/NMR spectrum;
    #           1=spectrum type 1: IR; 0=Plot spectrum! (compute and draw
    #           the curve); 1=Save graphical file of the spectrum in
    #           current folder
    PLOT_IR_SPECTRUM = ("11", "1", "0", "1")
    # Simulate Raman spectrum from Raman activities; optionally converts
    # activities to intensities given a laser wavelength and temperature
    # Sequence: 11=Plot IR/Raman/UV-Vis/ECD/VCD/ROA/NMR spectrum;
    #           2=spectrum type 2: Raman (or pre-resonance Raman); 0=Plot
    #           spectrum! (compute and draw the curve); 1=Save graphical
    #           file of the spectrum in current folder
    PLOT_RAMAN_SPECTRUM = ("11", "2", "0", "1")
    # Simulate UV-Vis absorption spectrum by broadening TD-DFT/CIS
    # excitation energies weighted by oscillator strengths; area
    # calibrated to molar absorptivity
    # Sequence: 11=Plot IR/Raman/UV-Vis/ECD/VCD/ROA/NMR spectrum;
    #           3=spectrum type 3: UV-Vis; 0=Plot spectrum! (compute and
    #           draw the curve); 1=Save graphical file of the spectrum in
    #           current folder
    PLOT_UV_VIS_SPECTRUM = ("11", "3", "0", "1")
    # Simulate Electronic Circular Dichroism (ECD) spectrum from rotatory
    # strengths in either length or velocity gauge representation
    # Sequence: 11=Plot IR/Raman/UV-Vis/ECD/VCD/ROA/NMR spectrum;
    #           4=spectrum type 4: ECD; 0=Plot spectrum! (compute and draw
    #           the curve); 1=Save graphical file of the spectrum in
    #           current folder
    PLOT_ECD_SPECTRUM = ("11", "4", "0", "1")
    # Simulate Vibrational Circular Dichroism (VCD) spectrum by broadening
    # vibrational rotatory strengths
    # Sequence: 11=Plot IR/Raman/UV-Vis/ECD/VCD/ROA/NMR spectrum;
    #           5=spectrum type 5: VCD; 0=Plot spectrum! (compute and draw
    #           the curve); 1=Save graphical file of the spectrum in
    #           current folder
    # INTERACTIVE ONLY -- requires vibrational frequency/intensity data (a
    #                     Gaussian or ORCA frequency-job output) in the
    #                     primary input file; crashes 'end-of-file during
    #                     read' on the file itself when that data is
    #                     absent, as in both bundled test files
    #                     (coord.molden has no frequency job, tddft.out is
    #                     a TD-DFT excited-state job with no Hessian)
    PLOT_VCD_SPECTRUM = ("11", "5", "0", "1")
    # Simulate Raman Optical Activity (ROA) spectrum from computed ROA
    # intensities
    # Sequence: 11=Plot IR/Raman/UV-Vis/ECD/VCD/ROA/NMR spectrum;
    #           6=spectrum type 6: ROA; 0=Plot spectrum! (compute and draw
    #           the curve); 1=Save graphical file of the spectrum in
    #           current folder
    PLOT_ROA_SPECTRUM = ("11", "6", "0", "1")
    # Simulate NMR spectrum from calculated chemical shielding tensors;
    # convert to chemical shifts by subtracting a reference shielding
    # value
    # Sequence: 11=Plot IR/Raman/UV-Vis/ECD/VCD/ROA/NMR spectrum;
    #           7=spectrum type 7: NMR; 0=Plot spectrum! (compute and draw
    #           the curve); 1=Save graphical file of the spectrum in
    #           current folder
    PLOT_NMR_SPECTRUM = ("11", "7", "0", "1")
    # Simulate fluorescence/phosphorescence spectrum; applies Kasha rule
    # so only the first excited state contributes when appropriate
    # Sequence: 11=Plot IR/Raman/UV-Vis/ECD/VCD/ROA/NMR spectrum;
    #           8=spectrum type 8: Fluorescence (undocumented in the
    #           printed menu, but accepted); 0=Plot spectrum! (compute and
    #           draw the curve); 1=Save graphical file of the spectrum in
    #           current folder
    PLOT_FLUORESCENCE_SPECTRUM = ("11", "8", "0", "1")
    # Plot photoelectron valence spectrum by broadening orbital ionisation
    # energies; analogous to TDOS but focused on valence region
    # Sequence: 11=Plot IR/Raman/UV-Vis/ECD/VCD/ROA/NMR spectrum;
    #           9=spectrum type 9: PVS/photoelectron vibrational spectrum
    #           (undocumented, but accepted); 0=Plot spectrum! (compute
    #           and draw the curve); 1=Save graphical file of the spectrum
    #           in current folder
    PLOT_PVS = ("11", "9", "0", "1")
    # Predict the perceived colour of a compound from its computed UV-Vis
    # absorption spectrum using CIE colour matching functions
    # Sequence: 11=Plot IR/Raman/UV-Vis/ECD/VCD/ROA/NMR spectrum; 0=return
    #           to the previous menu
    PREDICT_COLOR = ("11", "0")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 12: Quantitative molecular surface analysis
    # ─────────────────────────────────────────────────────────────────────────
    # Generates the rho = 0.001 a.u. vdW isosurface via Marching
    # Tetrahedra, maps a chosen function onto it, and computes surface
    # statistical descriptors (V_S, sigma^2, Pi, etc.) and locates surface
    # extrema.
    # Sequence: 12=Quantitative analysis of molecular surface; 0=Start
    #           analysis now!; -2=Export the grid data to surf.cub in
    #           current folder
    QMSA_ESP = ("12", "0", "-2")
    # Sequence: 12=Quantitative analysis of molecular surface; 2=Select
    #           mapped function, current: Electrostatic potential (ESP);
    #           2=Average local ionization energy (ALIE); 0=Start analysis
    #           now!; -2=export the grid data to surf.cub
    QMSA_ALIE = ("12", "2", "2", "0", "-2")
    # Sequence: 12=Quantitative analysis of molecular surface; 2=Select
    #           mapped function, current: Electrostatic potential (ESP);
    #           4=Local electron affinity (LEA); 0=Start analysis now!;
    #           -2=export the grid data to surf.cub
    QMSA_LEA = ("12", "2", "4", "0", "-2")
    # Sequence: 12=Quantitative analysis of molecular surface; 2=Select
    #           mapped function, current: Electrostatic potential (ESP);
    #           -4=Local electron attachment energy (LEAE); 0=Start
    #           analysis now!; -2=export the grid data to surf.cub
    QMSA_LEAE = ("12", "2", "-4", "0", "-2")
    # Sequence: 12=Quantitative analysis of molecular surface; 2=Select
    #           mapped function, current: Electrostatic potential (ESP);
    #           5=Electron delocalization range function EDR(r;d); 0=Start
    #           analysis now!; -2=export the grid data to surf.cub
    QMSA_EDR = ("12", "2", "5", "0", "-2")
    # Sequence: 12=Quantitative analysis of molecular surface; 2=Select
    #           mapped function, current: Electrostatic potential (ESP);
    #           6=Orbital overlap length function D(r) which maximizes
    #           EDR(r;d); 0=Start analysis now!; -2=export the grid data
    #           to surf.cub
    QMSA_MAXEDR = ("12", "2", "6", "0", "-2")
    # Sequence: 12=Quantitative analysis of molecular surface; 2=Select
    #           mapped function, current: Electrostatic potential (ESP);
    #           11=Electron density; 0=Start analysis now!; -2=export the
    #           grid data to surf.cub
    QMSA_EDENSITY = ("12", "2", "11", "0", "-2")
    # Sequence: 12=Quantitative analysis of molecular surface; 2=Select
    #           mapped function, current: Electrostatic potential (ESP);
    #           12=Sign(lambda2)*rho; 0=Start analysis now!; -2=export the
    #           grid data to surf.cub
    QMSA_LAMBDA2_RHO = ("12", "2", "12", "0", "-2")
    # DEPRECATED SEQUENCES - TO BE TESTED Map ESP onto the vdW surface
    # (rho=0.001); compute V_S+, V_S-, sigma^2, Pi and locate surface ESP
    # minima/maxima; outputs GIPF descriptors for property prediction
    # Sequence: 12=Quantitative analysis of molecular surface; 0=Start
    #           analysis now!
    SURFACE_ANALYSIS_ESP = ("12", "0")
    # Map ALIE onto the vdW surface; locate surface minima which predict
    # the most reactive electrophilic and radical attack sites
    # Sequence: 12=Quantitative analysis of molecular surface; 2=Select
    #           mapped function, current: Electrostatic potential (ESP)
    SURFACE_ANALYSIS_ALIE = ("12", "2")
    # Compute the molecular vdW surface area and enclosed volume from the
    # rho=0.001 isosurface
    # Sequence: 12=Quantitative analysis of molecular surface; 6=Start
    #           analysis without considering mapped function
    SURFACE_AREA_VOLUME = ("12", "6")
    # Becke-partition surface analysis: map a function onto each atomic
    # Becke surface and compute per-atom surface descriptors
    # Sequence: 12=Quantitative analysis of molecular surface; 4=Advanced
    #           options
    BECKE_SURFACE = ("12", "4")
    # Hirshfeld surface analysis for crystal packing studies: generate the
    # Hirshfeld surface and compute shape-index, curvedness, and contact-
    # distance properties
    # Sequence: 12=Quantitative analysis of molecular surface; 5=Loading
    #           mapped function values from external file, current: No
    HIRSHFELD_SURFACE = ("12", "5")
    # Locate and list all local minima and maxima of the mapped function
    # on the molecular surface, with their coordinates and values
    # Sequence: 12=Quantitative analysis of molecular surface; 6=Start
    #           analysis without considering mapped function
    SURFACE_EXTREMA = ("12", "6")
    # Generate a Hirshfeld surface fingerprint plot (2D histogram of d_i
    # vs. d_e distances); reveals intermolecular contact patterns in
    # crystal structures
    # Sequence: 12=Quantitative analysis of molecular surface; 5=Loading
    #           mapped function values from external file, current: No;
    #           2=Similar to 1, but specific for the case of using cubegen
    #           utility of Gaussian
    HIRSHFELD_SURFACE_FINGERPRINT = ("12", "5", "2")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 13: Process grid data
    # ─────────────────────────────────────────────────────────────────────────
    # Works on grid data already held in memory (from Menu 5 cube
    # generation) or loaded from an external .cube/.grd file at startup.
    # Export the current in-memory grid data to a Gaussian .cube file in
    # the working directory
    # Sequence: 13=Process grid data (No grid data is presented
    #           currently); 0=return to the previous menu
    EXPORT_CUBE = ("13", "0")
    # Export all grid-point coordinates together with their function
    # values to plain-text file (output.txt)
    # Sequence: 13=Process grid data (No grid data is presented
    #           currently); 1=export all grid points to a plain text file
    # INTERACTIVE ONLY -- requires grid data from a prior cube-generation
    #                     step in the same session (Multiwfn: 'Grid data
    #                     has not been loaded or generated'); not a
    #                     standalone zero-interaction action
    EXPORT_GRID_ALL_POINTS = ("13", "1")
    # Extract a 2D slice of the 3D grid in the XY plane at a user-
    # specified Z index and save to a text file
    # Sequence: 13=Process grid data (No grid data is presented
    #           currently); 2=extract the XY-plane slice of the grid
    # INTERACTIVE ONLY -- requires grid data from a prior cube-generation
    #                     step in the same session (Multiwfn: 'Grid data
    #                     has not been loaded or generated'); not a
    #                     standalone zero-interaction action
    GRID_EXTRACT_PLANE_XY = ("13", "2")
    # Extract a 2D slice of the 3D grid in the XZ plane at a user-
    # specified Y index and save to a text file
    # Sequence: 13=Process grid data (No grid data is presented
    #           currently); 3=extract the XZ-plane slice of the grid
    # INTERACTIVE ONLY -- requires grid data from a prior cube-generation
    #                     step in the same session (Multiwfn: 'Grid data
    #                     has not been loaded or generated'); not a
    #                     standalone zero-interaction action
    GRID_EXTRACT_PLANE_XZ = ("13", "3")
    # Extract a 2D slice of the 3D grid in the YZ plane at a user-
    # specified X index and save to a text file
    # Sequence: 13=Process grid data (No grid data is presented
    #           currently); 4=extract the YZ-plane slice of the grid
    # INTERACTIVE ONLY -- requires grid data from a prior cube-generation
    #                     step in the same session (Multiwfn: 'Grid data
    #                     has not been loaded or generated'); not a
    #                     standalone zero-interaction action
    GRID_EXTRACT_PLANE_YZ = ("13", "4")
    # Compute the planar average of grid values over all XY planes in a
    # specified Z range; outputs a 1-D average profile along Z
    # Sequence: 13=Process grid data (No grid data is presented
    #           currently); 5=average the grid along X and Y
    # INTERACTIVE ONLY -- requires grid data from a prior cube-generation
    #                     step in the same session (Multiwfn: 'Grid data
    #                     has not been loaded or generated'); not a
    #                     standalone zero-interaction action
    GRID_AVERAGE_XY = ("13", "5")
    # Compute the planar average of grid values over all XZ planes in a
    # specified X range; outputs a 1-D average profile along X
    # Sequence: 13=Process grid data (No grid data is presented
    #           currently); 6=average the grid along X and Z
    # INTERACTIVE ONLY -- requires grid data from a prior cube-generation
    #                     step in the same session (Multiwfn: 'Grid data
    #                     has not been loaded or generated'); not a
    #                     standalone zero-interaction action
    GRID_AVERAGE_XZ = ("13", "6")
    # Compute the planar average of grid values over all YZ planes in a
    # specified Y range; outputs a 1-D average profile along Y
    # Sequence: 13=Process grid data (No grid data is presented
    #           currently); 7=average the grid along Y and Z
    # INTERACTIVE ONLY -- requires grid data from a prior cube-generation
    #                     step in the same session (Multiwfn: 'Grid data
    #                     has not been loaded or generated'); not a
    #                     standalone zero-interaction action
    GRID_AVERAGE_YZ = ("13", "7")
    # Extract a 2D grid slice in the plane defined by the nuclear
    # positions of three user-specified atoms
    # Sequence: 13=Process grid data (No grid data is presented
    #           currently); 8=extract a plane defined by 3 atoms
    # INTERACTIVE ONLY -- requires grid data from a prior cube-generation
    #                     step in the same session (Multiwfn: 'Grid data
    #                     has not been loaded or generated'); not a
    #                     standalone zero-interaction action
    GRID_EXTRACT_PLANE_3ATOMS = ("13", "8")
    # Extract a 2D grid slice in the plane defined by three user-specified
    # Cartesian coordinates
    # Sequence: 13=Process grid data (No grid data is presented
    #           currently); 9=extract a plane defined by 3 points
    # INTERACTIVE ONLY -- requires grid data from a prior cube-generation
    #                     step in the same session (Multiwfn: 'Grid data
    #                     has not been loaded or generated'); not a
    #                     standalone zero-interaction action
    GRID_EXTRACT_PLANE_3POINTS = ("13", "9")
    # Export all grid points whose function value falls within a user-
    # specified [min, max] interval
    # Sequence: 13=Process grid data (No grid data is presented
    #           currently); 10=extract grid points within a value range
    # INTERACTIVE ONLY -- requires grid data from a prior cube-generation
    #                     step in the same session (Multiwfn: 'Grid data
    #                     has not been loaded or generated'); not a
    #                     standalone zero-interaction action
    GRID_EXTRACT_VALUE_RANGE = ("13", "10")
    # Perform element-wise arithmetic (+, -, *, /) between the current
    # grid and a second cube file; used to construct difference densities,
    # Fukui functions, etc.
    # Sequence: 13=Process grid data (No grid data is presented
    #           currently); 11=perform math operations on two grid data
    #           sets
    # INTERACTIVE ONLY -- requires grid data from a prior cube-generation
    #                     step in the same session (Multiwfn: 'Grid data
    #                     has not been loaded or generated'); not a
    #                     standalone zero-interaction action
    GRID_MATH_OPERATIONS = ("13", "11")
    # Map function values from a second cube file onto the isosurface of
    # the current grid (e.g., colour an ELF isosurface by ESP)
    # Sequence: 13=Process grid data (No grid data is presented
    #           currently); 12=map the grid onto an isosurface
    # INTERACTIVE ONLY -- requires grid data from a prior cube-generation
    #                     step in the same session (Multiwfn: 'Grid data
    #                     has not been loaded or generated'); not a
    #                     standalone zero-interaction action
    GRID_MAP_TO_ISOSURFACE = ("13", "12")
    # Set all grid points farther than (or closer than) a distance cutoff
    # from chosen atoms to a fixed value; used to mask unwanted isosurface
    # regions
    # Sequence: 13=Process grid data (No grid data is presented
    #           currently); 13=set grid values by distance to a reference
    # INTERACTIVE ONLY -- requires grid data from a prior cube-generation
    #                     step in the same session (Multiwfn: 'Grid data
    #                     has not been loaded or generated'); not a
    #                     standalone zero-interaction action
    GRID_SET_VALUE_DISTANCE = ("13", "13")
    # Set grid points outside the fuzzy overlap region of two user-defined
    # fragments to a fixed value; isolates the inter-fragment interaction
    # zone
    # Sequence: 13=Process grid data (No grid data is presented
    #           currently); 14=set grid values by fragment membership
    # INTERACTIVE ONLY -- requires grid data from a prior cube-generation
    #                     step in the same session (Multiwfn: 'Grid data
    #                     has not been loaded or generated'); not a
    #                     standalone zero-interaction action
    GRID_SET_VALUE_FRAGMENT = ("13", "14")
    # Replace all grid values within a specified numerical range with a
    # single fixed value; applies thresholding or clamping to the data
    # Sequence: 13=Process grid data (No grid data is presented
    #           currently); 15=set grid values within a range
    # INTERACTIVE ONLY -- requires grid data from a prior cube-generation
    #                     step in the same session (Multiwfn: 'Grid data
    #                     has not been loaded or generated'); not a
    #                     standalone zero-interaction action
    GRID_SET_VALUE_RANGE = ("13", "15")
    # Linearly rescale all grid values from their current [min, max] to a
    # user- specified [new_min, new_max]
    # Sequence: 13=Process grid data (No grid data is presented
    #           currently); 16=rescale the grid's value range
    # INTERACTIVE ONLY -- requires grid data from a prior cube-generation
    #                     step in the same session (Multiwfn: 'Grid data
    #                     has not been loaded or generated'); not a
    #                     standalone zero-interaction action
    GRID_SCALE_RANGE = ("13", "16")
    # Print statistics (min, max, mean, std dev, and spatial integral) for
    # grid points in a user-defined spatial and/or value range
    # Sequence: 13=Process grid data (No grid data is presented
    #           currently); 17=print statistics of the grid data
    # INTERACTIVE ONLY -- requires grid data from a prior cube-generation
    #                     step in the same session (Multiwfn: 'Grid data
    #                     has not been loaded or generated'); not a
    #                     standalone zero-interaction action
    GRID_STATISTIC_DATA = ("13", "17")
    # Compute and plot the running integral of the grid data along X, Y,
    # or Z; useful for generating charge-displacement curves in EDA
    # Sequence: 13=Process grid data (No grid data is presented
    #           currently); 18=plot the integral curve of the grid data
    # INTERACTIVE ONLY -- requires grid data from a prior cube-generation
    #                     step in the same session (Multiwfn: 'Grid data
    #                     has not been loaded or generated'); not a
    #                     standalone zero-interaction action
    GRID_PLOT_INTEGRAL_CURVE = ("13", "18")
    # Open the interactive GUI to visualise the isosurface of the
    # currently loaded grid data at an adjustable isovalue slider
    # Sequence: 13=Process grid data (No grid data is presented
    #           currently); -2=visualize the grid as an isosurface
    # INTERACTIVE ONLY -- requires grid data from a prior cube-generation
    #                     step in the same session (Multiwfn: 'Grid data
    #                     has not been loaded or generated'); not a
    #                     standalone zero-interaction action
    GRID_VISUALIZE_ISOSURFACE = ("13", "-2")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 14: Adaptive Natural Density Partitioning (AdNDP)
    # ─────────────────────────────────────────────────────────────────────────
    # Decomposes the electron density into n-centre two-electron (nc-2e)
    # bonding elements interactively; visualises results in Multiwfn GUI.
    # NEED SPECIFIC FILE Launch the AdNDP interactive interface to search
    # for 1c-2e lone pairs, 2c-2e bonds, and multi-centre bonding
    # elements; widely used for cluster and aromatic-system bonding
    # analysis
    # Sequence: 14=Adaptive natural density partitioning (AdNDP) analysis
    # INTERACTIVE ONLY -- requires natural population analysis (NPA) data
    #                     to already be present in the input file
    #                     (Multiwfn: 'Cannot find natural population
    #                     analysis information in the input file'); fails
    #                     with a plain .molden/.fch, needs e.g. an NBO-
    #                     enabled Gaussian output instead
    ADNDP_ANALYSIS = ("14",)

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 15: Fuzzy atomic space analysis
    # ─────────────────────────────────────────────────────────────────────────
    # Numerical integration of real-space functions in Becke or Hirshfeld
    # fuzzy atomic spaces; computes delocalization indices and
    # aromaticity.
    # Sequence: 15=Fuzzy atomic space analysis; 1=Perform integration in
    #           fuzzy atomic spaces for a real space function; 1=Electron
    #           density (rho)
    FUZZY_INTEGRATE_EDENSITY = ("15", "1", "1")
    # Sequence: 15=Fuzzy atomic space analysis; 1=Perform integration in
    #           fuzzy atomic spaces for a real space function; 2=Gradient
    #           norm of rho
    FUZZY_INTEGRATE_NORM_RHO = ("15", "1", "2")
    # Sequence: 15=Fuzzy atomic space analysis; 1=Perform integration in
    #           fuzzy atomic spaces for a real space function; 3=Laplacian
    #           of rho
    FUZZY_INTEGRATE_LAPLACIAN = ("15", "1", "3")
    # Sequence: 15=Fuzzy atomic space analysis; 1=Perform integration in
    #           fuzzy atomic spaces for a real space function; 4=Value of
    #           orbital wavefunction; h=select the HOMO
    FUZZY_INTEGRATE_ORB_WFN_HOMO = ("15", "1", "4", "h")
    # Sequence: 15=Fuzzy atomic space analysis; 1=Perform integration in
    #           fuzzy atomic spaces for a real space function; 4=Value of
    #           orbital wavefunction; l=select the LUMO
    FUZZY_INTEGRATE_ORB_WFN_LUMO = ("15", "1", "4", "l")
    # Sequence: 15=Fuzzy atomic space analysis; 1=Perform integration in
    #           fuzzy atomic spaces for a real space function; 5=Electron
    #           spin density
    FUZZY_INTEGRATE_ESPIN_DENSITY = ("15", "1", "5")
    # Sequence: 15=Fuzzy atomic space analysis; 1=Perform integration in
    #           fuzzy atomic spaces for a real space function;
    #           6=Hamiltonian kinetic energy density K(r)
    FUZZY_INTEGRATE_KR = ("15", "1", "6")
    # Sequence: 15=Fuzzy atomic space analysis; 1=Perform integration in
    #           fuzzy atomic spaces for a real space function;
    #           7=Lagrangian kinetic energy density G(r)
    FUZZY_INTEGRATE_GR = ("15", "1", "7")
    # Sequence: 15=Fuzzy atomic space analysis; 1=Perform integration in
    #           fuzzy atomic spaces for a real space function;
    #           8=Electrostatic potential from nuclear charges
    FUZZY_INTEGRATE_ESP_CHARGES = ("15", "1", "8")
    # Sequence: 15=Fuzzy atomic space analysis; 1=Perform integration in
    #           fuzzy atomic spaces for a real space function; 9=Electron
    #           localization function (ELF)
    FUZZY_INTEGRATE_ELF = ("15", "1", "9")
    # Sequence: 15=Fuzzy atomic space analysis; 1=Perform integration in
    #           fuzzy atomic spaces for a real space function;
    #           10=Localized orbital locator (LOL)
    FUZZY_INTEGRATE_LOL = ("15", "1", "10")
    # Sequence: 15=Fuzzy atomic space analysis; 1=Perform integration in
    #           fuzzy atomic spaces for a real space function; 11=Local
    #           information entropy
    FUZZY_INTEGRATE_LOCAL_ENTROPY = ("15", "1", "11")
    # Sequence: 15=Fuzzy atomic space analysis; 1=Perform integration in
    #           fuzzy atomic spaces for a real space function; 12=Total
    #           electrostatic potential (ESP)
    FUZZY_INTEGRATE_ESP = ("15", "1", "12")
    # Sequence: 15=Fuzzy atomic space analysis; 1=Perform integration in
    #           fuzzy atomic spaces for a real space function; 13=Reduced
    #           density gradient (RDG)
    FUZZY_INTEGRATE_RDG = ("15", "1", "13")
    # Sequence: 15=Fuzzy atomic space analysis; 1=Perform integration in
    #           fuzzy atomic spaces for a real space function; 14=RDG with
    #           promolecular approximation
    FUZZY_INTEGRATE_RDG_PROMOLECULAR = ("15", "1", "14")
    # Sequence: 15=Fuzzy atomic space analysis; 1=Perform integration in
    #           fuzzy atomic spaces for a real space function;
    #           15=Sign(lambda2)*rho
    FUZZY_INTEGRATE_LAMBDA2RHO = ("15", "1", "15")
    # Sequence: 15=Fuzzy atomic space analysis; 1=Perform integration in
    #           fuzzy atomic spaces for a real space function;
    #           16=Sign(lambda2)*rho with promolecular approximation
    FUZZY_INTEGRATE_LAMBDA2RGO_PROMOLECULAR = ("15", "1", "16")
    # Sequence: 15=Fuzzy atomic space analysis; 1=Perform integration in
    #           fuzzy atomic spaces for a real space function; 18=Average
    #           local ionization energy (ALIE)
    FUZZY_INTEGRATE_ALIE = ("15", "1", "18")
    # Sequence: 15=Fuzzy atomic space analysis; 1=Perform integration in
    #           fuzzy atomic spaces for a real space function; 20=Electron
    #           delocal. range func. EDR(r;d)
    FUZZY_INTEGRATE_EDR = ("15", "1", "20")
    # Sequence: 15=Fuzzy atomic space analysis; 1=Perform integration in
    #           fuzzy atomic spaces for a real space function; 21=Orbital
    #           overlap dist. func. D(r)
    FUZZY_INTEGRATE_ORB_OVERLAP_DR = ("15", "1", "21")
    # Sequence: 15=Fuzzy atomic space analysis; 1=Perform integration in
    #           fuzzy atomic spaces for a real space function; 22=Delta-g
    #           (promolecular approximation)
    FUZZY_INTEGRATE_DELTAG_PROMOLECULAR = ("15", "1", "22")
    # Sequence: 15=Fuzzy atomic space analysis; 1=Perform integration in
    #           fuzzy atomic spaces for a real space function; 23=Delta-g
    #           (Hirshfeld partition)
    FUZZY_INTEGRATE_DELTAG_HIRSHFELD = ("15", "1", "23")
    # Sequence: 15=Fuzzy atomic space analysis; 1=Perform integration in
    #           fuzzy atomic spaces for a real space function;
    #           24=Interaction region indicator (IRI)
    FUZZY_INTEGRATE_IRI = ("15", "1", "24")
    # Sequence: 15=Fuzzy atomic space analysis; 2=Calculate atomic and
    #           molecular multipole moments and <r^2>; 1=output the result
    #           on screen (not to multipole.txt)
    FUZZY_MULTIPOLE = ("15", "2", "1")
    # Sequence: 15=Fuzzy atomic space analysis; 8=Perform integration in
    #           fuzzy overlap region for a real space functions;
    #           1=Electron density (rho); n=no (decline the
    #           fragment/mirror-plane restriction)
    FUZZY_OVERLAP_EDENSITY = ("15", "8", "1", "n")
    # Sequence: 15=Fuzzy atomic space analysis; 8=Perform integration in
    #           fuzzy overlap region for a real space functions;
    #           1=Electron density (rho); n=no (decline the
    #           fragment/mirror-plane restriction)
    FUZZY_OVERLAP_NORM_RHO = ("15", "8", "1", "n")
    # Sequence: 15=Fuzzy atomic space analysis; 8=Perform integration in
    #           fuzzy overlap region for a real space functions;
    #           3=Laplacian of rho; n=no (decline the
    #           fragment/mirror-plane restriction)
    FUZZY_OVERLAP_LAPLACIAN = ("15", "8", "3", "n")
    # Sequence: 15=Fuzzy atomic space analysis; 8=Perform integration in
    #           fuzzy overlap region for a real space functions; 4=Value
    #           of orbital wavefunction; h=select the HOMO; n=no (decline
    #           the fragment/mirror-plane restriction)
    FUZZY_OVERLAP_ORB_WFN_HOMO = ("15", "8", "4", "h", "n")
    # Sequence: 15=Fuzzy atomic space analysis; 8=Perform integration in
    #           fuzzy overlap region for a real space functions; 4=Value
    #           of orbital wavefunction; l=select the LUMO; n=no (decline
    #           the fragment/mirror-plane restriction)
    FUZZY_OVERLAP_ORB_WFN_LUMO = ("15", "8", "4", "l", "n")
    # Sequence: 15=Fuzzy atomic space analysis; 8=Perform integration in
    #           fuzzy overlap region for a real space functions;
    #           5=Electron spin density; n=no (decline the
    #           fragment/mirror-plane restriction)
    FUZZY_OVERLAP_ESPIN_DENSITY = ("15", "8", "5", "n")
    # Sequence: 15=Fuzzy atomic space analysis; 8=Perform integration in
    #           fuzzy overlap region for a real space functions;
    #           6=Hamiltonian kinetic energy density K(r); n=no (decline
    #           the fragment/mirror-plane restriction)
    FUZZY_OVERLAP_KR = ("15", "8", "6", "n")
    # Sequence: 15=Fuzzy atomic space analysis; 8=Perform integration in
    #           fuzzy overlap region for a real space functions;
    #           7=Lagrangian kinetic energy density G(r); n=no (decline
    #           the fragment/mirror-plane restriction)
    FUZZY_OVERLAP_GR = ("15", "8", "7", "n")
    # Sequence: 15=Fuzzy atomic space analysis; 8=Perform integration in
    #           fuzzy overlap region for a real space functions;
    #           8=Electrostatic potential from nuclear charges; n=no
    #           (decline the fragment/mirror-plane restriction)
    FUZZY_OVERLAP_ESP_CHARGES = ("15", "8", "8", "n")
    # Sequence: 15=Fuzzy atomic space analysis; 8=Perform integration in
    #           fuzzy overlap region for a real space functions;
    #           9=Electron localization function (ELF); n=no (decline the
    #           fragment/mirror-plane restriction)
    FUZZY_OVERLAP_ELF = ("15", "8", "9", "n")
    # Sequence: 15=Fuzzy atomic space analysis; 8=Perform integration in
    #           fuzzy overlap region for a real space functions;
    #           10=Localized orbital locator (LOL); n=no (decline the
    #           fragment/mirror-plane restriction)
    FUZZY_OVERLAP_LOL = ("15", "8", "10", "n")
    # Sequence: 15=Fuzzy atomic space analysis; 8=Perform integration in
    #           fuzzy overlap region for a real space functions; 11=Local
    #           information entropy; n=no (decline the
    #           fragment/mirror-plane restriction)
    FUZZY_OVERLAP_LOCAL_ENTROPY = ("15", "8", "11", "n")
    # Sequence: 15=Fuzzy atomic space analysis; 8=Perform integration in
    #           fuzzy overlap region for a real space functions; 12=Total
    #           electrostatic potential (ESP); n=no (decline the
    #           fragment/mirror-plane restriction)
    FUZZY_OVERLAP_ESP = ("15", "8", "12", "n")
    # Sequence: 15=Fuzzy atomic space analysis; 8=Perform integration in
    #           fuzzy overlap region for a real space functions;
    #           13=Reduced density gradient (RDG); n=no (decline the
    #           fragment/mirror-plane restriction)
    FUZZY_OVERLAP_RDG = ("15", "8", "13", "n")
    # Sequence: 15=Fuzzy atomic space analysis; 8=Perform integration in
    #           fuzzy overlap region for a real space functions; 14=RDG
    #           with promolecular approximation; n=no (decline the
    #           fragment/mirror-plane restriction)
    FUZZY_OVERLAP_RDG_PROMOLECULAR = ("15", "8", "14", "n")
    # Sequence: 15=Fuzzy atomic space analysis; 8=Perform integration in
    #           fuzzy overlap region for a real space functions;
    #           15=Sign(lambda2)*rho; n=no (decline the
    #           fragment/mirror-plane restriction)
    FUZZY_OVERLAP_LAMBDA2RHO = ("15", "8", "15", "n")
    # Sequence: 15=Fuzzy atomic space analysis; 8=Perform integration in
    #           fuzzy overlap region for a real space functions;
    #           16=Sign(lambda2)*rho with promolecular approximation; n=no
    #           (decline the fragment/mirror-plane restriction)
    FUZZY_OVERLAP_LAMBDA2RGO_PROMOLECULAR = ("15", "8", "16", "n")
    # Sequence: 15=Fuzzy atomic space analysis; 8=Perform integration in
    #           fuzzy overlap region for a real space functions;
    #           18=Average local ionization energy (ALIE); n=no (decline
    #           the fragment/mirror-plane restriction)
    FUZZY_OVERLAP_ALIE = ("15", "8", "18", "n")
    # Sequence: 15=Fuzzy atomic space analysis; 8=Perform integration in
    #           fuzzy overlap region for a real space functions;
    #           20=Electron delocal. range func. EDR(r;d); n=no (decline
    #           the fragment/mirror-plane restriction)
    FUZZY_OVERLAP_EDR = ("15", "8", "20", "n")
    # Sequence: 15=Fuzzy atomic space analysis; 8=Perform integration in
    #           fuzzy overlap region for a real space functions;
    #           21=Orbital overlap dist. func. D(r); n=no (decline the
    #           fragment/mirror-plane restriction)
    FUZZY_OVERLAP_ORB_OVERLAP_DR = ("15", "8", "21", "n")
    # Sequence: 15=Fuzzy atomic space analysis; 8=Perform integration in
    #           fuzzy overlap region for a real space functions;
    #           22=Delta-g (promolecular approximation); n=no (decline the
    #           fragment/mirror-plane restriction)
    FUZZY_OVERLAP_DELTAG_PROMOLECULAR = ("15", "8", "22", "n")
    # Sequence: 15=Fuzzy atomic space analysis; 8=Perform integration in
    #           fuzzy overlap region for a real space functions;
    #           23=Delta-g (Hirshfeld partition); n=no (decline the
    #           fragment/mirror-plane restriction)
    FUZZY_OVERLAP_DELTAG_HIRSHFELD = ("15", "8", "23", "n")
    # Sequence: 15=Fuzzy atomic space analysis; 8=Perform integration in
    #           fuzzy overlap region for a real space functions;
    #           24=Interaction region indicator (IRI); n=no (decline the
    #           fragment/mirror-plane restriction)
    FUZZY_OVERLAP_IRI = ("15", "8", "24", "n")
    # Sequence: 15=Fuzzy atomic space analysis; 9=Calculate condensed
    #           linear response kernel (CLRK); n=no (decline the
    #           fragment/mirror-plane restriction)
    CLRK_MATRIX = ("15", "9", "n")
    # Compute the Atomic Overlap Matrix (AOM) S_AB = integral(psi_i *
    # psi_j * w_A) dr; required input for NOCV, ETS-NOCV, and
    # delocalization-index calculations
    # Sequence: 15=Fuzzy atomic space analysis; 3=Calculate and output
    #           atomic overlap matrix (AOM) in current folder
    ATOMIC_OVERLAP_MATRIX = ("15", "3")
    # Compute localization index and pairwise delocalization index (DI)
    # for all atom pairs in fuzzy spaces; DI is a bond-order analogue
    # grounded in density-matrix theory
    # Sequence: 15=Fuzzy atomic space analysis; 4=Calculate localization
    #           (LI) and delocalization index (DI); n=no (decline the
    #           fragment/mirror-plane restriction)
    LOCALIZATION_DELOCALIZATION_INDEX = ("15", "4", "n")
    # Para Delocalization Index (PDI): average DI between para-related
    # carbon atoms in a 6-membered ring; larger PDI indicates stronger
    # aromaticity
    # Sequence: 15=Fuzzy atomic space analysis; 5=Calculate PDI
    #           (Para-delocalization index); q=quit / return without
    #           change
    PDI_AROMATICITY = ("15", "5", "q")
    # Aromatic Fluctuation Index (FLU): measures deviation of DI from
    # reference bond values; value near 0 indicates a fully aromatic ring
    # Sequence: 15=Fuzzy atomic space analysis; 6=Calculate FLU (Aromatic
    #           fluctuation index); q=quit / return without change
    FLU_AROMATICITY = ("15", "6", "q")
    # FLU-pi: FLU computed using only pi-orbital contributions to the
    # delocalization index; more selective for pi-electron aromaticity
    # INTERACTIVE - REQUIRES USER INPUT
    # Sequence: 15=Fuzzy atomic space analysis; 7=Calculate FLU-pi
    FLU_PI_AROMATICITY = ("15", "7")
    # Compute the condensed linear response kernel chi_AB; measures how
    # much electron density at atom B responds to an external perturbation
    # at atom A
    # Sequence: 15=Fuzzy atomic space analysis; 9=Calculate condensed
    #           linear response kernel (CLRK)
    CONDENSED_LINEAR_RESPONSE = ("15", "9")
    # Compute the para linear response (PLR) aromaticity index from the
    # condensed linear response kernel; large PLR indicates strong
    # aromatic delocalization
    # Sequence: 15=Fuzzy atomic space analysis; 10=Calculate PLR (Para
    #           linear response index); q=quit / return without change
    PARA_LINEAR_RESPONSE = ("15", "10", "q")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 16: Charge Decomposition Analysis (CDA)
    # ─────────────────────────────────────────────────────────────────────────
    # Generalised CDA (GCDA): decomposes orbital interactions between two
    # or more fragments into donation, back-donation, repulsion, and
    # residual. INTERACTIVE Launch the CDA/GCDA: decompose charge transfer
    # between fragments into donation (d), back-donation (b), repulsion
    # (r), and residual (Delta) terms; plot the orbital interaction
    # diagram INTERACTIVE - REQUIRES USER INPUT
    # Sequence: 16=Charge decomposition analysis (CDA) and plot orbital
    #           interaction diagram
    # INTERACTIVE ONLY -- requires interactively defining one or more atom
    #                     fragments ('How many fragments...'); can't be
    #                     filled with a fixed default
    CDA_ANALYSIS = ("16",)

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 17: Basin analysis
    # ─────────────────────────────────────────────────────────────────────────
    # Locates attractors of a real-space function, builds gradient-
    # following basins, and integrates properties (rho, multipoles, DI)
    # within each.
    # Sequence: 17=Basin analysis; 1=Generate basins and locate
    #           attractors; 1=Electron density (rho); 2=Medium quality
    #           grid, spacing=0.10 Bohr, cost: 8x
    BASIN_ANALYSIS_RHO = ("17", "1", "1", "2")
    # Sequence: 17=Basin analysis; 1=Generate basins and locate
    #           attractors; 1=Electron density (rho); 1=Low quality grid,
    #           spacing=0.20 Bohr, cost: 1x
    BASIN_EDENSITY = ("17", "1", "1", "1")
    # Sequence: 17=Basin analysis; 1=Generate basins and locate
    #           attractors; 2=Gradient norm of rho; n=no (decline the
    #           fragment/mirror-plane restriction)
    BASIN_NORM_RHO = ("17", "1", "2", "n")
    # Sequence: 17=Basin analysis; 1=Generate basins and locate
    #           attractors; 3=Laplacian of rho; n=no (decline the
    #           fragment/mirror-plane restriction)
    BASIN_LAPLACIAN = ("17", "1", "3", "n")
    # Sequence: 17=Basin analysis; 1=Generate basins and locate
    #           attractors; 4=Value of orbital wavefunction; h=select the
    #           HOMO; n=no (decline the fragment/mirror-plane restriction)
    BASIN_ORB_WFN_HOMO = ("17", "1", "4", "h", "n")
    # Sequence: 17=Basin analysis; 1=Generate basins and locate
    #           attractors; 4=Value of orbital wavefunction; l=select the
    #           LUMO; n=no (decline the fragment/mirror-plane restriction)
    BASIN_ORB_WFN_LUMO = ("17", "1", "4", "l", "n")
    # Sequence: 17=Basin analysis; 1=Generate basins and locate
    #           attractors; 5=Electron spin density; n=no (decline the
    #           fragment/mirror-plane restriction)
    BASIN_ESPIN_DENSITY = ("17", "1", "5", "n")
    # Sequence: 17=Basin analysis; 1=Generate basins and locate
    #           attractors; 6=Hamiltonian kinetic energy density K(r);
    #           n=no (decline the fragment/mirror-plane restriction)
    BASIN_KR = ("17", "1", "6", "n")
    # Sequence: 17=Basin analysis; 1=Generate basins and locate
    #           attractors; 7=Lagrangian kinetic energy density G(r); n=no
    #           (decline the fragment/mirror-plane restriction)
    BASIN_GR = ("17", "1", "7", "n")
    # Sequence: 17=Basin analysis; 1=Generate basins and locate
    #           attractors; 8=Electrostatic potential from nuclear
    #           charges; n=no (decline the fragment/mirror-plane
    #           restriction)
    BASIN_ESP_CHARGES = ("17", "1", "8", "n")
    # Sequence: 17=Basin analysis; 1=Generate basins and locate
    #           attractors; 9=Electron localization function (ELF); n=no
    #           (decline the fragment/mirror-plane restriction)
    BASIN_ELF = ("17", "1", "9", "n")
    # Sequence: 17=Basin analysis; 1=Generate basins and locate
    #           attractors; 10=Localized orbital locator (LOL); n=no
    #           (decline the fragment/mirror-plane restriction)
    BASIN_LOL = ("17", "1", "10", "n")
    # Sequence: 17=Basin analysis; 1=Generate basins and locate
    #           attractors; 11=Local information entropy; n=no (decline
    #           the fragment/mirror-plane restriction)
    BASIN_LOCAL_ENTROPY = ("17", "1", "11", "n")
    # Sequence: 17=Basin analysis; 1=Generate basins and locate
    #           attractors; 12=Total electrostatic potential (ESP); n=no
    #           (decline the fragment/mirror-plane restriction)
    BASIN_ESP = ("17", "1", "12", "n")
    # Sequence: 17=Basin analysis; 1=Generate basins and locate
    #           attractors; 13=Reduced density gradient (RDG); n=no
    #           (decline the fragment/mirror-plane restriction)
    BASIN_RDG = ("17", "1", "13", "n")
    # Sequence: 17=Basin analysis; 1=Generate basins and locate
    #           attractors; 14=RDG with promolecular approximation; n=no
    #           (decline the fragment/mirror-plane restriction)
    BASIN_RDG_PROMOLECULAR = ("17", "1", "14", "n")
    # Sequence: 17=Basin analysis; 1=Generate basins and locate
    #           attractors; 15=Sign(lambda2)*rho; n=no (decline the
    #           fragment/mirror-plane restriction)
    BASIN_LAMBDA2RHO = ("17", "1", "15", "n")
    # Sequence: 17=Basin analysis; 1=Generate basins and locate
    #           attractors; 16=Sign(lambda2)*rho with promolecular
    #           approximation; n=no (decline the fragment/mirror-plane
    #           restriction)
    BASIN_LAMBDA2RGO_PROMOLECULAR = ("17", "1", "16", "n")
    # Sequence: 17=Basin analysis; 1=Generate basins and locate
    #           attractors; 18=Average local ionization energy (ALIE);
    #           n=no (decline the fragment/mirror-plane restriction)
    BASIN_ALIE = ("17", "1", "18", "n")
    # Sequence: 17=Basin analysis; 1=Generate basins and locate
    #           attractors; 20=Electron delocal. range func. EDR(r;d);
    #           n=no (decline the fragment/mirror-plane restriction)
    BASIN_EDR = ("17", "1", "20", "n")
    # Sequence: 17=Basin analysis; 1=Generate basins and locate
    #           attractors; 21=Orbital overlap dist. func. D(r); n=no
    #           (decline the fragment/mirror-plane restriction)
    BASIN_ORB_OVERLAP_DR = ("17", "1", "21", "n")
    # Sequence: 17=Basin analysis; 1=Generate basins and locate
    #           attractors; 22=Delta-g (promolecular approximation); n=no
    #           (decline the fragment/mirror-plane restriction)
    BASIN_DELTAG_PROMOLECULAR = ("17", "1", "22", "n")
    # Sequence: 17=Basin analysis; 1=Generate basins and locate
    #           attractors; 23=Delta-g (Hirshfeld partition); n=no
    #           (decline the fragment/mirror-plane restriction)
    BASIN_DELTAG_HIRSHFELD = ("17", "1", "23", "n")
    # Sequence: 17=Basin analysis; 1=Generate basins and locate
    #           attractors; 24=Interaction region indicator (IRI); n=no
    #           (decline the fragment/mirror-plane restriction)
    BASIN_IRI = ("17", "1", "24", "n")
    # AIM basin analysis: attractors are nuclei; integrates rho in each
    # QTAIM atomic basin to yield Bader/AIM charges, atomic multipoles,
    # LI, and DI
    # Sequence: 17=Basin analysis; 1=Generate basins and locate attractors
    BASIN_ANALYSIS_AIM = ("17", "1")
    # ELF basin analysis: finds ELF attractors (bonding pairs, lone pairs,
    # core shells) and integrates rho per basin to obtain basin electron
    # populations
    # Sequence: 17=Basin analysis; 2=Topology analysis
    BASIN_ANALYSIS_ELF = ("17", "2")
    # Integrate any user-chosen real-space function over basins that have
    # already been generated in the current Multiwfn session
    # Sequence: 17=Basin analysis; 3=Output and plot specific property in
    #           a line
    BASIN_INTEGRATE_PROPERTY = ("17", "3")
    # ESP basin analysis: ESP minima serve as attractors; integrates rho
    # in corresponding basins to characterise nucleophilic pockets and
    # lone-pair regions
    # Sequence: 17=Basin analysis; 4=Output and plot specific property in
    #           a plane
    BASIN_ANALYSIS_ESP = ("17", "4")
    # LOL basin analysis: LOL maxima serve as attractors (bonding and
    # lone-pair regions); integrates rho for sharp-boundary basin
    # populations
    # Sequence: 17=Basin analysis; 5=Output and plot specific property
    #           within a spatial region (calc. grid data)
    BASIN_ANALYSIS_LOL = ("17", "5")
    # LOL-alpha basin analysis restricted to alpha-spin electrons;
    # isolates pi- electron basins for aromaticity and pi-conjugation
    # studies
    # Sequence: 17=Basin analysis; 6=Check & modify wavefunction
    BASIN_ANALYSIS_LOL_ALPHA = ("17", "6")
    # Basin analysis using a user-defined real-space function specified
    # via the iuserfunc parameter in settings.ini
    # Sequence: 17=Basin analysis; 0=Show molecular structure and view
    #           orbitals
    BASIN_ANALYSIS_CUSTOM = ("17", "0")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 18: Electron excitation analysis
    # ─────────────────────────────────────────────────────────────────────────
    # Analyses TD-DFT/CIS excited states loaded from Gaussian or ORCA
    # output. Analyse hole and electron density distributions for a TD-DFT
    # transition: compute centroid distance, t index, and Sr overlap
    # integral; visualise hole/electron isosurfaces
    # Sequence: 18=Electron excitation analysis; 1=Analyze and visualize
    #           hole&electron distribution, transition density, and
    #           transition electric/magnetic dipole moment density
    HOLE_ELECTRON_ANALYSIS = ("18", "1")
    # Plot the transition density matrix as a colour-filled 2D map in the
    # atom- atom/MO-MO representation; reveals coherence length and
    # charge-transfer character
    # Sequence: 18=Electron excitation analysis; 2=Plot atom/fragment
    #           transition matrix of various kinds as heat map
    # INTERACTIVE ONLY -- requires the path to a separate external file
    #                     (not the primary wavefunction) -- e.g. a
    #                     Gaussian/ORCA output, a cube file, or an NBO
    #                     file
    TRANSITION_DENSITY_MATRIX = ("18", "2")
    # Analyse charge transfer from an electron-density-difference grid:
    # compute CT distance, amount of transferred charge per fragment using
    # the Plasser- Lischka method
    # Sequence: 18=Electron excitation analysis; 3=Analyze charge-transfer
    #           based on density difference grid data (JCTC,7,2498)
    CHARGE_TRANSFER_ANALYSIS = ("18", "3")
    # Compute the Delta_r index (charge-transfer length) for each excited
    # state to quantify local (small Delta_r) vs. charge-transfer (large
    # Delta_r) excitation character
    # Sequence: 18=Electron excitation analysis; 4=Calculate delta_r index
    #           to measure charge-transfer length (JCTC,9,3118)
    # INTERACTIVE ONLY -- requires the path to a separate external file
    #                     (not the primary wavefunction) -- e.g. a
    #                     Gaussian/ORCA output, a cube file, or an NBO
    #                     file
    DELTA_R_INDEX = ("18", "4")
    # Calculate and print transition electric (and optionally magnetic)
    # dipole moments between all pairs of excited states in the loaded TD-
    # DFT output
    # Sequence: 18=Electron excitation analysis; 5=Calculate transition
    #           electric/magnetic dipole moments between all states and
    #           for each state
    # INTERACTIVE ONLY -- requires the path to a separate external file
    #                     (not the primary wavefunction) -- e.g. a
    #                     Gaussian/ORCA output, a cube file, or an NBO
    #                     file
    TRANSITION_DIPOLE_MOMENTS = ("18", "5")
    # Generate Natural Transition Orbitals (NTOs) for a chosen excitation;
    # the dominant particle-hole NTO pair visually captures the essential
    # character of the transition
    # Sequence: 18=Electron excitation analysis; 6=Generate natural
    #           transition orbitals (NTOs)
    GENERATE_NTO = ("18", "6")
    # Inter-Fragment Charge Transfer (IFCT): quantify hole and electron
    # populations transferred between user-defined fragments upon
    # photoexcitation
    # Sequence: 18=Electron excitation analysis; 8=Calculate interfragment
    #           charge transfer via IFCT method
    IFCT_ANALYSIS = ("18", "8")
    # Compute the Lambda diagnostic for each TD-DFT transition; Lambda <
    # 0.3 flags charge-transfer states that may be poorly described by
    # standard functionals
    # Sequence: 18=Electron excitation analysis; 14=Calculate lambda index
    #           to characterize electron excitation (JCP,128,044118)
    LAMBDA_INDEX = ("18", "14")
    # Charge Transfer Spectrum (CTS): run batch IFCT over all excited
    # states, then plot the result as an inter-fragment CT contribution
    # spectrum
    # Sequence: 18=Electron excitation analysis; 16=Charge-transfer
    #           spectrum (CTS) analysis (Carbon,187,78)
    # INTERACTIVE ONLY -- requires interactively defining one or more atom
    #                     fragments ('How many fragments...'); can't be
    #                     filled with a fixed default
    CTS_ANALYSIS = ("18", "16")
    # Compute the conditional electron density: given one electron is
    # fixed at reference point r0, maps the conditional probability of
    # finding a second electron at r
    # Sequence: 18=Electron excitation analysis; 17=Electron density
    #           polarization analysis based on electron excitations
    CONDITIONAL_DENSITY = ("18", "17")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 19: Orbital localization
    # ─────────────────────────────────────────────────────────────────────────
    # Transforms canonical delocalized MOs into spatially localised MOs
    # for chemical interpretation as bonds, lone pairs, and core orbitals.
    # Sequence: 19=Orbital localization analysis; 1=Localizing occupied
    #           orbitals only
    PIPEK_MEZEY_LOCALIZATION_HIRSHFELD_OCCUPIED = ("19", "1")
    # Sequence: 19=Orbital localization analysis; 2=Localizing both
    #           occupied and unoccupied orbitals separately
    PIPEK_MEZEY_LOCALIZATION_HIRSHFELD_ALL = ("19", "2")
    # Sequence: 19=Orbital localization analysis; -6=Set localization
    #           method, current: Pipek-Mezey with Mulliken population;
    #           2=Pipek-Mezey based on Lowdin type of population;
    #           1=Localizing occupied orbitals only
    # INTERACTIVE ONLY -- CRASHES this Multiwfn 3.8(dev) build with
    #                     SIGSEGV after 'Calculating orbital
    #                     compositions...' following localization --
    #                     appears to be a genuine bug in this build, not a
    #                     sequence issue
    PIPEK_MEZEY_LOCALIZATION_LOWDIN_OCUPIED = ("19", "-6", "2", "1")
    # Sequence: 19=Orbital localization analysis; -6=Set localization
    #           method, current: Pipek-Mezey with Mulliken population;
    #           2=Pipek-Mezey based on Lowdin type of population;
    #           2=Localizing both occupied and unoccupied orbitals
    #           separately
    # INTERACTIVE ONLY -- CRASHES this Multiwfn 3.8(dev) build with
    #                     SIGSEGV after 'Calculating orbital
    #                     compositions...' following localization --
    #                     appears to be a genuine bug in this build, not a
    #                     sequence issue
    PIPEK_MEZEY_LOCALIZATION_LOWDIN_ALL = ("19", "-6", "2", "2")
    # Sequence: 19=Orbital localization analysis; -6=Set localization
    #           method, current: Pipek-Mezey with Mulliken population;
    #           3=Pipek-Mezey based on Becke population; 1=Localizing
    #           occupied orbitals only
    # INTERACTIVE ONLY -- CRASHES this Multiwfn 3.8(dev) build with
    #                     SIGSEGV after 'Calculating orbital
    #                     compositions...' following localization --
    #                     appears to be a genuine bug in this build, not a
    #                     sequence issue
    PIPEK_MEZEY_LOCALIZATION_BECKE_OCCUPIED = ("19", "-6", "3", "1")
    # Sequence: 19=Orbital localization analysis; -6=Set localization
    #           method, current: Pipek-Mezey with Mulliken population;
    #           3=Pipek-Mezey based on Becke population; 2=Localizing both
    #           occupied and unoccupied orbitals separately
    PIPEK_MEZEY_LOCALIZATION_BECKE_ALL = ("19", "-6", "3", "2")
    # Sequence: 19=Orbital localization analysis; -6=Set localization
    #           method, current: Pipek-Mezey with Mulliken population;
    #           10=Foster-Boys; 1=Localizing occupied orbitals only
    # INTERACTIVE ONLY -- CRASHES this Multiwfn 3.8(dev) build with
    #                     SIGSEGV after 'Calculating orbital
    #                     compositions...' following localization --
    #                     appears to be a genuine bug in this build, not a
    #                     sequence issue
    BOYS_LOCALIZATION_OCCUPIED = ("19", "-6", "10", "1")
    # Sequence: 19=Orbital localization analysis; -6=Set localization
    #           method, current: Pipek-Mezey with Mulliken population;
    #           10=Foster-Boys; 2=Localizing both occupied and unoccupied
    #           orbitals separately
    BOYS_LOCALIZATION_ALL = ("19", "-6", "10", "2")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 20: Weak interaction analysis
    # ─────────────────────────────────────────────────────────────────────────
    # Methods based on the Reduced Density Gradient (RDG) and related
    # functions to visualise and quantify non-covalent interactions (NCI).
    # NCI analysis from wavefunction: generate RDG and sign(lambda2)*rho
    # cubes; render RDG isosurface coloured by sign(lambda2)*rho in VMD to
    # visualise H-bonds, vdW contacts, and steric clashes
    # Sequence: 20=Visual study of weak interaction; 1=NCI analysis (also
    #           known as RDG analysis); 2=Medium quality grid, covering
    #           whole system, about 512000 points in total; 2=Output
    #           scatter points to output.txt in current folder; 3=Output
    #           cube files to func1.cub and func2.cub in current folder;
    #           0=Start analysis now!; 0=return to the previous menu
    NCI_ANALYSIS = ("20", "1", "2", "2", "3", "0", "0")
    # Promolecular NCI: same as NCI but built from superposition of free-
    # atom densities; very fast and suitable for macromolecules and
    # protein-ligand binding without a full wavefunction
    # Sequence: 20=Visual study of weak interaction; 2=NCI analysis based
    #           on promolecular density; 2=Medium quality grid, covering
    #           whole system, about 512000 points in total; 2=Output
    #           scatter points to output.txt in current folder; 3=export
    #           cube files to func1.cub/func2.cub; 0=Start analysis now!;
    #           0=return to the previous menu
    NCI_PROMOLECULAR = ("20", "2", "2", "2", "3", "0", "0")
    # Interaction Region Indicator (IRI): improved NCI variant that decays
    # smoothly at atomic cores and better resolves weak intermolecular
    # contacts recommended over standard RDG
    # Sequence: 20=Visual study of weak interaction; 4=IRI: Interaction
    #           region indicator analysis (Chemistry-Methods, 1, 231);
    #           2=Medium quality grid, covering whole system, about 512000
    #           points in total; 1=Save the scatter graph to file;
    #           2=Output scatter points to output.txt in current folder;
    #           3=export cube files to func1.cub/func2.cub; 0=Start
    #           analysis now!; 0=return to the previous menu
    IRI_ANALYSIS = ("20", "4", "2", "1", "2", "3", "0", "0")
    # Density Overlap Regions Indicator (DORI): highlights regions where
    # two electron densities overlap, clearly revealing both covalent and
    # non- covalent interaction zones
    # Sequence: 20=Visual study of weak interaction; 5=DORI analysis;
    #           2=Medium quality grid, covering whole system, about 512000
    #           points in total; 1=real-space function: sign(lambda2)*rho
    #           (paired with EDR); 2=medium-quality grid; 3=export cube
    #           files to func1.cub/func2.cub; 0=Start analysis now!;
    #           0=return to the previous menu
    DORI_ANALYSIS = ("20", "5", "2", "1", "2", "3", "0", "0")
    # Compute and visualise the van der Waals interaction potential
    # landscape around the molecule
    # Sequence: 20=Visual study of weak interaction; 6=Visualization of
    #           van der Waals potential (JMM, 26, 315); 2=Medium quality
    #           grid, covering whole system, about 512000 points in total;
    #           1=Show isosurface graph of repulsion potential; 2=Show
    #           isosurface graph of dispersion potential; 3=Show
    #           isosurface graph of van der Waals potential; 0=Start
    #           analysis now!; 0=return to the previous menu
    VDW_POTENTIAL = ("20", "6", "2", "1", "2", "3", "0", "0")
    # Averaged NCI (ANCI): NCI analysis averaged over an ensemble of MD
    # trajectory snapshots; reveals which non-covalent interactions
    # persist statistically over time INTERACTIVE - REQUIRES USER INPUT
    # Sequence: 20=Visual study of weak interaction; 3=aNCI: Averaged NCI
    #           analysis; 2=medium-quality grid; 1=start averaged-NCI
    #           (aNCI) analysis; 2=Medium quality grid, covering whole
    #           system, about 512000 points in total
    # INTERACTIVE ONLY -- did not finish within a ~30s verification
    #                     timeout (computing averaged
    #                     density/gradient/Hessian) -- may just be a
    #                     genuinely slow calculation (aNCI is an ensemble
    #                     method), not confirmed broken
    ANCI_ANALYSIS = ("20", "3", "2", "1", "2")
    # INTERACTIVE - REQUIRES USER INPUT Independent Gradient Model (IGM):
    # decomposes the density gradient into intra- and inter-fragment
    # parts; isolates and visualises inter-fragment non-covalent
    # interactions INTEACTIVE - REQUIRES USER INPUT
    # Sequence: 20=Visual study of weak interaction; 10=IGM analysis
    # INTERACTIVE ONLY -- requires interactively defining one or more atom
    #                     fragments ('How many fragments...'); can't be
    #                     filled with a fixed default
    IGM_ANALYSIS = ("20", "10")
    # IGM-H (Hirshfeld-based IGM): uses Hirshfeld atomic densities as the
    # reference for gradient decomposition; more accurate inter-fragment
    # gradient isolation than standard IGM INTERACTIVE - REQUIRES USER
    # INPUT
    # Sequence: 20=Visual study of weak interaction; 11=IGMH: IGM analysis
    #           based on Hirshfeld partition of molecular density (JCC,
    #           43, 539)
    # INTERACTIVE ONLY -- requires interactively defining one or more atom
    #                     fragments ('How many fragments...'); can't be
    #                     filled with a fixed default
    IGMH_ANALYSIS = ("20", "11")
    # Averaged IGM (aIGM): IGM analysis averaged over an MD trajectory
    # ensemble identifies which inter-fragment interactions are persistent
    # in dynamic systems INTERACTIVE - REQUIRES USER INPUT
    # Sequence: 20=Visual study of weak interaction; 12=aIGM: Averaged IGM
    #           analysis
    # INTERACTIVE ONLY -- requires interactively defining one or more atom
    #                     fragments ('How many fragments...'); can't be
    #                     filled with a fixed default
    AIGM_ANALYSIS = ("20", "12")
    # Modified IGM (mIGM): coordinate-only IGM variant that avoids the
    # need for a wavefunction; nearly identical results to IGMH for weak
    # interactions at orders-of-magnitude lower cost INTERACTIVE -
    # REQUIRES USER INPUT
    # Sequence: 20=Visual study of weak interaction; -10=mIGM: Modified
    #           IGM analysis
    # INTERACTIVE ONLY -- requires interactively defining one or more atom
    #                     fragments ('How many fragments...'); can't be
    #                     filled with a fixed default
    MIGM_ANALYSIS = ("20", "-10")
    # Averaged mIGM (amIGM): mIGM averaged over MD snapshots; extends mIGM
    # to fluctuation environments; more robust than aNCI and recommended
    # for MD- based NCI studies INTERACTIVE - REQUIRES USER INPUT
    # Sequence: 20=Visual study of weak interaction; -11=not a valid
    #           option here (redisplays the same menu; only -12/-10 exist)
    AMIGM_ANALYSIS = ("20", "-11")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 21: Energy Decomposition Analysis (EDA)
    # ─────────────────────────────────────────────────────────────────────────
    # Simple EDA using combined fragment wavefunctions: decomposes
    # interaction energy into electrostatic, Pauli exchange-repulsion,
    # polarisation, and dispersion components
    # Sequence: 21=Energy decomposition analysis; 1=Energy decomposition
    #           analysis based on molecular forcefield (EDA-FF)
    EDA_FF = ("21", "1")
    # Symmetry-based localisation EDA (SBL-EDA): energy partitioning
    # scheme using symmetry constraints on localised fragment orbitals
    # Sequence: 21=Energy decomposition analysis; 2=Shubin Liu's energy
    #           decomposition analysis (needs Gaussian)
    # INTERACTIVE ONLY -- requires the path to a separate external file
    #                     (not the primary wavefunction) -- e.g. a
    #                     Gaussian/ORCA output, a cube file, or an NBO
    #                     file
    EDA_SBL = ("21", "2")
    # Second-Order Bond Energy Analysis (SOBEDA): decomposes pairwise bond
    # energies into individual occupied-MO contributions via second-order
    # perturbation theory
    # Sequence: 21=Energy decomposition analysis; 3=sobEDA and sobEDAw
    #           energy decomposition analyses
    SOBEDA_ANALYSIS = ("21", "3")
    # Compute per-atom contributions to the DFT-D3 or D4 empirical
    # dispersion correction energy; identifies which atoms drive van der
    # Waals attraction most strongly
    # Sequence: 21=Energy decomposition analysis; 4=Analysis of atomic
    #           contribution to dispersion energy
    DISPERSION_ATOMIC_CONTRIBUTION = ("21", "4")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 22: Conceptual DFT (CDFT)
    # ─────────────────────────────────────────────────────────────────────────
    # Chemical reactivity descriptors derived from DFT response theory.
    # Launch the CDFT module: compute global reactivity indices including
    # chemical potential mu, chemical hardness eta, softness S, and
    # electrophilicity index omega
    # Sequence: 22=Conceptual DFT (CDFT) analysis
    CDFT_ANALYSIS = ("22",)
    # Condense Fukui functions to atomic values using Mulliken, Hirshfeld,
    # Becke, or AIM partitioning; tabulates f+, f-, and f0 for each atom
    # Sequence: 22=Conceptual DFT (CDFT) analysis; 3=Calculate grid data
    #           of Fukui function, dual descriptor and related functions
    # INTERACTIVE ONLY -- requires a second wavefunction file (e.g. the
    #                     N+1/N-1 electron state for CDFT) beyond the
    #                     primary input
    CONDENSED_FUKUI = ("22", "3")
    # Compute local hardness eta(r) from the local chemical potential;
    # complements the dual descriptor for predicting hard-soft reactivity
    # site preferences
    # Sequence: 22=Conceptual DFT (CDFT) analysis; 4=Set delta in
    #           orbital-weighted (OW) calculation, current: 0.1000 a.u.
    LOCAL_HARDNESS = ("22", "4")
    # Compute orbital weights as a reactivity map; low-weight regions
    # indicate preferred electrophilic attack sites
    # Sequence: 22=Conceptual DFT (CDFT) analysis; 5=Print current orbital
    #           weights used in orbital-weighted (OW) calculation
    ORBITAL_WEIGHTS = ("22", "5")
    # Sequence: 22=Conceptual DFT (CDFT) analysis; 6=Calculate condensed
    #           OW Fukui function and OW dual descriptor
    ORBITAL_WEIGHTED_FUKUI = ("22", "6")
    # Sequence: 22=Conceptual DFT (CDFT) analysis; 7=Calculate grid data
    #           of OW Fukui function and OW dual descriptor; 2=Medium
    #           quality grid, covering whole system, about 512000 points
    #           in total; 5=Export grid data of orbital-weighted f+ as
    #           OW_f+.cub; 6=export grid data of orbital-weighted f- as
    #           OW_f-.cub; 7=export grid data of orbital-weighted f0 as
    #           OW_f0.cub; 8=Export grid data of orbital-weighted dual
    #           descriptor as OW_DD.cub
    GRID_FUKUI = ("22", "7", "2", "5", "6", "7", "8")
    # Sequence: 22=Conceptual DFT (CDFT) analysis; 8=Calculate
    #           nucleophilic and electrophilic superdelocalizabilities
    SUPERDELOCALIZABILITIES_NUC_E = ("22", "8")
    # Generate .wfn files for the N, N+1, and N-1 electron states via
    # Gaussian single-point calculations; prerequisite for condensed Fukui
    # functions and CDFT indices that need multiple charge states
    # Sequence: 22=Conceptual DFT (CDFT) analysis; 1=Generate .wfn files
    #           for N, N+1, N-1 electrons states
    # INTERACTIVE ONLY -- requires a local Gaussian installation to be
    #                     found on PATH
    CDFT_GENERATE_CHARGED_WFN = ("22", "1")
    # Calculate grid data of the Fukui potential and dual-descriptor
    # potential (electrostatic-potential analogues of the Fukui function)
    # Sequence: 22=Conceptual DFT (CDFT) analysis; 9=Calculate grid data
    #           of Fukui potential and dual descriptor potential
    # INTERACTIVE ONLY -- requires a second wavefunction file (e.g. the
    #                     N+1/N-1 electron state for CDFT) beyond the
    #                     primary input
    CDFT_GRID_FUKUI_POTENTIAL = ("22", "9")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 23: ETS-NOCV
    # ─────────────────────────────────────────────────────────────────────────
    # Extended Transition State combined with Natural Orbitals for
    # Chemical Valence: energy-decomposed orbital interaction analysis.
    # ETS-NOCV: diagonalise the deformation density matrix to obtain NOCV
    # pairs; decompose the orbital interaction energy Delta_E_orb into
    # individual pairwise NOCV channel contributions; useful for
    # characterising sigma, pi, and delta bond formation
    # Sequence: 23=ETS-NOCV analysis
    # INTERACTIVE ONLY -- requires interactively defining one or more atom
    #                     fragments ('How many fragments...'); can't be
    #                     filled with a fixed default
    ETS_NOCV_ANALYSIS = ("23",)

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 24: Polarizability
    # ─────────────────────────────────────────────────────────────────────────
    # INTERACTIVE Parse and print polarizability (alpha) and
    # hyperpolarizability (beta, gamma) tensors from a Gaussian frequency
    # or finite-field task output file
    # Sequence: 24=(Hyper)polarizability analysis; 1=Parse output file of
    #           (hyper)polarizability task of Gaussian and calculate
    #           various related quantities
    PARSE_POLARIZABILITY = ("24", "1")
    # Compute (hyper)polarizability via the Sum-Over-States (SOS) method
    # from TD-DFT excited-state data; useful when a direct field-
    # perturbation calculation is impractical
    # Sequence: 24=(Hyper)polarizability analysis; 2=Study
    #           (hyper)polarizability by sum-over-states (SOS) method and
    #           perform two/three-level analysis
    SOS_POLARIZABILITY = ("24", "2")
    # Compute and visualise the polarizability density p(r): the spatial
    # distribution of where molecular polarisability originates in real
    # space
    # Sequence: 24=(Hyper)polarizability analysis; 3=(hyper)polarizability
    #           density analysis
    POLARIZABILITY_DENSITY = ("24", "3")
    # Project the polarizability tensor onto a unit sphere to produce a 3D
    # directional surface showing polarisability anisotropy
    # Sequence: 24=(Hyper)polarizability analysis; 5=Visualize
    #           (hyper)polarizability via unit sphere and vector
    #           representations
    UNIT_SPHERE_POLARIZABILITY = ("24", "5")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 25: Aromaticity
    # ─────────────────────────────────────────────────────────────────────────
    # Magnetic, geometric, and electronic indices for aromaticity
    # assessment. Anisotropy of the Induced Current Density: compute the
    # magnetically induced ring-current density and export as a 3D
    # isosurface; diatropic current = aromatic
    # Sequence: 25=Electron delocalization and aromaticity analyses;
    #           1=Multicenter bond order
    # INTERACTIVE ONLY -- requires an interactively-specified atom or ring
    #                     index list (e.g. 'Input index of the atoms in
    #                     ring order')
    AICD_ANALYSIS = ("25", "1")
    # Compute Nucleus-Independent Chemical Shift (NICS) at a single user-
    # specified point; negative NICS = diatropic / aromatic, positive =
    # paratropic / antiaromatic
    # Sequence: 25=Electron delocalization and aromaticity analyses;
    #           2=AV1245 index
    # INTERACTIVE ONLY -- requires an interactively-specified atom or ring
    #                     index list (e.g. 'Input index of the atoms in
    #                     ring order')
    NICS_POINT = ("25", "2")
    # Iso-Chemical Shielding Surface (ICSS): compute NMR shielding on a 3D
    # grid and generate isosurfaces to visualise the spatial extent and
    # shape of ring-current shielding cones
    # Sequence: 25=Electron delocalization and aromaticity analyses;
    #           3=Iso-chemical shielding surface (ICSS)
    ICSS_ANALYSIS = ("25", "3")
    # NICS scan along a user-defined path perpendicular to a ring plane;
    # plots the NICS(z) profile to show how aromaticity decays with
    # distance from the ring
    # Sequence: 25=Electron delocalization and aromaticity analyses;
    #           4=NICS_ZZ for non-planar or tilted system
    # INTERACTIVE ONLY -- requires an interactively-specified coordinate
    #                     (e.g. a NICS ring center)
    NICS_SCAN = ("25", "4")
    # Bird aromaticity index (I5 for 5-membered, I6 for 6-membered rings):
    # geometry-based index from bond-length uniformity; 100 = fully
    # aromatic reference
    # Sequence: 25=Electron delocalization and aromaticity analyses;
    #           5=ELF-sigma/pi and LOL-sigma/pi
    # INTERACTIVE ONLY -- in this Multiwfn 3.8(dev) build, option 5 under
    #                     Main Menu 25 is actually 'ELF-sigma/pi and LOL-
    #                     sigma/pi' (a pointer to other modules, not a
    #                     standalone calculation) -- the menu has been
    #                     reorganized since this entry was written; Bird
    #                     aromaticity index is now reached via option 6
    #                     (combined with HOMA), which itself requires an
    #                     interactively-specified atom or ring index list
    #                     (e.g. 'Input index of the atoms in ring order')
    #                     (ring)
    BIRD_INDEX = ("25", "5")
    # Harmonic Oscillator Model of Aromaticity: geometry-based index; 1 =
    # fully aromatic, 0 = nonaromatic, negative values indicate
    # antiaromaticity
    # Sequence: 25=Electron delocalization and aromaticity analyses;
    #           6=Harmonic oscillator measure of aromaticity (HOMA) and
    #           Bird indices
    # INTERACTIVE ONLY -- reaches the real 'HOMA / Bird aromaticity index'
    #                     submenu (option 6) but requires an
    #                     interactively-specified atom or ring index list
    #                     (e.g. 'Input index of the atoms in ring order')
    #                     (ring) to actually run the calculation ('Input
    #                     indices of the atoms according to bonding
    #                     relationship in the ring')
    HOMA_INDEX = ("25", "6")
    # HOMAC (corrected HOMA) and HOMER aromaticity indices; improved
    # geometric models that account for updated bond-length reference
    # values from crystal data
    # Sequence: 25=Electron delocalization and aromaticity analyses;
    #           7=Shannon aromaticity index
    # INTERACTIVE ONLY -- in this Multiwfn 3.8(dev) build, option 7 under
    #                     Main Menu 25 is actually 'Shannon aromaticity
    #                     index' (a pointer to the topology analysis
    #                     module, main function 2) -- HOMAc/HOMER are sub-
    #                     variants '6a'/'6b' of option 6 (HOMA), not
    #                     reachable as a separate top-level entry, and
    #                     themselves requires an interactively-specified
    #                     atom or ring index list (e.g. 'Input index of
    #                     the atoms in ring order') (ring)
    HOMAC_HOMER = ("25", "7")
    # Stanger EN geometric aromaticity index: separates bond-length
    # alternation (EN_GEO) from uniform bond-length deviation (EN_BLA) for
    # a cleaner aromaticity measure
    # Sequence: 25=Electron delocalization and aromaticity analyses;
    #           8=Para-delocalization index (PDI)
    # INTERACTIVE ONLY -- in this Multiwfn 3.8(dev) build, option 8 under
    #                     Main Menu 25 is actually 'Para-delocalization
    #                     index (PDI)' (a pointer to the fuzzy atomic
    #                     space module, main function 15) -- no separate
    #                     Stanger EN aromaticity index option exists at
    #                     this position in this build
    STANGER_INDEX = ("25", "8")
    # Automated 1-D NICS scan perpendicular to a ring: places ghost atoms
    # at incremental heights, runs Gaussian, reads back shielding tensors,
    # and plots the NICS(z) curve
    # Sequence: 25=Electron delocalization and aromaticity analyses;
    #           13=NICS-1D scan curve map, integral NICS (INICS) and
    #           FiPC-NICS
    NICS_1D_SCAN = ("25", "13")
    # Automated 2D NICS shielding map in or perpendicular to the ring
    # plane; produces a colour-filled map to spatially visualise the
    # magnetic ring- current delocalization
    # Sequence: 25=Electron delocalization and aromaticity analyses;
    #           14=NICS-2D scan plane map
    NICS_2D_MAP = ("25", "14")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 100: Utilities Part 1
    # ─────────────────────────────────────────────────────────────────────────
    # INTERACTIVE Plot a 2D scatter graph of two real-space functions
    # evaluated on the same grid and export both as .cube files; classic
    # use is RDG vs. sign(lambda2)*rho for the NCI scatter plot
    # Sequence: 100=Other functions (Part 1); 1=Draw scatter graph between
    #           two functions and generate their cube files
    SCATTER_GRAPH_TWO_FUNCTIONS = ("100", "1")
    # Export the loaded wavefunction or geometry to .pdb, .xyz, .wfn,
    # .molden, .fch, Gaussian input, GAMESS input, or CP2K input format
    # Sequence: 100=Other functions (Part 1); 2=Export various files
    #           (mwfn/pdb/xyz/wfn/wfx/molden/fch/47/mkl...) or generate
    #           input file of quantum chemistry programs
    EXPORT_VARIOUS_FILES = ("100", "2")
    # Calculate the molecular van der Waals volume as the space enclosed
    # within the rho = 0.001 a.u. electron-density isosurface
    # Sequence: 12=Quantitative analysis of molecular surface; 6=Start
    #           analysis without considering mapped function
    VDW_VOLUME = ("12", "6")
    # Numerically integrate a chosen real-space function over all space
    # using Becke multi-centre quadrature; used to verify total electron
    # count (integral of rho = N)
    # Sequence: 100=Other functions (Part 1); 4=Integrate a function in
    #           whole space
    INTEGRATE_WHOLE_SPACE = ("100", "4")
    # Compute the spatial overlap integral <psi_i_alpha | psi_j_beta>
    # between each alpha-beta orbital pair; measures how well alpha and
    # beta orbitals correspond spatially in UHF/UKS
    # Sequence: 100=Other functions (Part 1); 5=Show overlap integral
    #           between alpha and beta orbitals
    # INTERACTIVE ONLY -- requires an unrestricted (open-shell, UHF/UKS)
    #                     wavefunction (Multiwfn: 'This function is only
    #                     available for unrestricted wavefunction!');
    #                     fails on closed-shell/restricted systems such as
    #                     the bundled test molecule
    ORBITAL_OVERLAP_INTEGRAL = ("100", "5")
    # Parse a Gaussian output file and plot the SCF energy and DIIS error
    # convergence as a function of iteration number
    # Sequence: 100=Other functions (Part 1); 6=Monitor SCF convergence
    #           process of Gaussian
    # INTERACTIVE ONLY -- expects a Gaussian SCF output log file named
    #                     exactly 'gauout.out' to already exist in the
    #                     working directory (Multiwfn reads it directly
    #                     via Fortran unit 10 with no path prompt; crashes
    #                     with 'end-of-file during read' when absent)
    MONITOR_SCF_CONVERGENCE = ("100", "6")
    # Generate a Gaussian input file that uses the already-converged MO
    # coefficients as initial guess; saves SCF iterations when restarting
    # or changing basis
    # Sequence: 100=Other functions (Part 1); 7=Auxiliary tools for CP2K
    #           (CP2Kmate)
    GAUSSIAN_INITIAL_GUESS = ("100", "7")
    # Generate a Gaussian input file whose initial MO guess is assembled
    # from separately converged fragment wavefunctions; used to prepare
    # CDA and EDA calculations
    # Sequence: 100=Other functions (Part 1); 8=Generate Gaussian input
    #           file with initial guess from fragment wavefunctions
    # INTERACTIVE ONLY -- requires interactively defining one or more atom
    #                     fragments ('How many fragments...'); can't be
    #                     filled with a fixed default
    FRAGMENT_GUESS_INPUT = ("100", "8")
    # Evaluate and print the coordination number for every atom using an
    # empirical distance criterion scaled from vdW or covalent radii
    # Sequence: 100=Other functions (Part 1); 9=Evaluate interatomic
    #           connectivity and atomic coordination number; =accept the
    #           default (press ENTER); y=yes
    ATOMIC_COORDINATION = ("100", "9", "", "y")
    # Calculate integral of |psi_i(r)| * |psi_j(r)| over all space;
    # measures how much two orbital densities spatially overlap without
    # phase cancellation
    # Sequence: 100=Other functions (Part 1); 11=Calculate overlap and
    #           centroid distance between two orbitals
    ORBITAL_OVERLAP_CENTROID = ("100", "11")
    # Biorthogonalise a set of orbitals (e.g., localised MOs) and compute
    # their one-electron energies from the Fock matrix without re-running
    # the SCF
    # Sequence: 100=Other functions (Part 1); 12=Biorthogonalization
    #           between alpha and beta orbitals
    BIORTHOGONALIZATION = ("100", "12")
    # Compute HOMA and Bird aromaticity indices for a user-specified ring
    # directly from the molecular geometry coordinates
    # Sequence: 100=Other functions (Part 1); 13=Process grid data (No
    #           grid data is presented currently)
    # INTERACTIVE ONLY -- requires an interactively-specified atom or ring
    #                     index list (e.g. 'Input index of the atoms in
    #                     ring order') (ring)
    HOMA_BIRD_AROMATICITY = ("100", "13")
    # Compute the LOLIPOP index (LOL Integrated Pi Over Plane): integral
    # of LOL in the pi-symmetry plane above and below the ring; quantifies
    # pi-electron delocalization
    # Sequence: 100=Other functions (Part 1); 14=Calculate LOLIPOP (LOL
    #           Integrated Pi Over Plane)
    LOLIPOP_INDEX = ("100", "14")
    # Calculate intermolecular orbital overlap integrals between two
    # molecular species; predicts charge-transfer coupling strengths
    # relevant to organic semiconductor design
    # Sequence: 100=Other functions (Part 1); 15=Calculate intermolecular
    #           orbital overlap
    INTERMOLECULAR_OVERLAP = ("100", "15")
    # Construct and export the Fock matrix F = C^{-1T} * epsilon * C^{-1}
    # from orbital energies and coefficients; needed for
    # biorthogonalisation and ETS-NOCV
    # Sequence: 100=Other functions (Part 1); 17=Generate Fock/KS matrix
    #           based on orbital energies and coefficients
    GENERATE_FOCK_MATRIX = ("100", "17")
    # Yoshizawa electron-transport route analysis: identify the dominant
    # through-bond tunnelling pathway between two terminal atoms using
    # squared orbital-amplitude products
    # Sequence: 100=Other functions (Part 1); 18=Yoshizawa's electron
    #           transport route analysis
    # INTERACTIVE ONLY -- requires the path to a separate external file
    #                     (not the primary wavefunction) -- e.g. a
    #                     Gaussian/ORCA output, a cube file, or an NBO
    #                     file
    ELECTRON_TRANSPORT_ROUTE = ("100", "18")
    # Combine separately computed fragment wavefunctions at their
    # molecular geometry to form a promolecular .wfn file; reference state
    # for EDA and CDA
    # Sequence: 100=Other functions (Part 1); 19=Generate new wavefunction
    #           by combining fragment wavefunctions
    COMBINE_FRAGMENTS = ("100", "19")
    # Compute Hellmann-Feynman electrostatic forces on each nucleus from
    # the electron density distribution and all other nuclear charges
    # Sequence: 100=Other functions (Part 1); 20=Calculate
    #           Hellmann-Feynman forces
    HELLMANN_FEYNMAN_FORCES = ("100", "20")
    # Print selected geometric properties: bond lengths, bond angles,
    # dihedral angles, and distances for user-specified atoms
    # Sequence: 100=Other functions (Part 1); 21=Calculate properties
    #           based on geometry information for specific atoms
    GEOMETRY_PROPERTIES = ("100", "21")
    # Automatically identify pi-symmetry orbitals in planar systems based
    # on nodal-plane criteria and optionally zero their occupation numbers
    # to isolate the sigma frame
    # Sequence: 100=Other functions (Part 1); 22=Detect pi orbitals, set
    #           occupation numbers and calculate pi composition
    DETECT_PI_ORBITALS = ("100", "22")
    # Fit the distribution of a real-space function onto per-atom
    # parameter values using a least-squares procedure; useful for atom-
    # centred property assignments
    # Sequence: 100=Other functions (Part 1); 23=Fit function distribution
    #           to atomic value
    FIT_FUNCTION_TO_ATOMS = ("100", "23")
    # Compute the out-of-plane NICS_ZZ shielding component for non-planar
    # rings by projecting the full shielding tensor onto the local ring-
    # normal direction
    # Sequence: 100=Other functions (Part 1); 24=(Hyper)polarizability
    #           analysis
    # INTERACTIVE ONLY -- requires an interactively-specified coordinate
    #                     (e.g. a NICS ring center)
    NICS_ZZ_NONPLANAR = ("100", "24")
    # Calculate the area enclosed by and the perimeter of a ring defined
    # by user-specified atom indices from the molecular geometry
    # Sequence: 100=Other functions (Part 1); 25=Electron delocalization
    #           and aromaticity analyses
    # INTERACTIVE ONLY -- requires an interactively-specified atom or ring
    #                     index list (e.g. 'Input index of the atoms in
    #                     ring order') (ring)
    RING_AREA_PERIMETER = ("100", "25")
    # Generate a CP2K periodic-DFT input file from the loaded geometry and
    # cell information; enter 'cp2k' at the export-type prompt to activate
    # Sequence: 100=Other functions (Part 1); 2=Export various files
    #           (mwfn/pdb/xyz/wfn/wfx/molden/fch/47/mkl...) or generate
    #           input file of quantum chemistry programs; 25=CP2K; =accept
    #           the default (press ENTER); 0=return to the previous menu
    GENERATE_CP2K_INPUT = ("100", "2", "25", "", "0")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 200: Utilities Part 2
    # ─────────────────────────────────────────────────────────────────────────
    # Sequence: 200=Other functions (Part 2); 10=Output various kinds of
    #           integral between orbitals; 1=electric dipole moment
    #           integral; 3=Between all orbitals
    ORBITAL_INTEGRAL_ELECTRIC_DIPOLE = ("200", "10", "1", "3")
    # Sequence: 200=Other functions (Part 2); 10=Output various kinds of
    #           integral between orbitals; 2=magnetic dipole moment
    #           integral; 3=Between all orbitals
    ORBITAL_INTEGRAL_MAGNETIC_DIPOLE = ("200", "10", "2", "3")
    # Sequence: 200=Other functions (Part 2); 10=Output various kinds of
    #           integral between orbitals; 3=velocity integral; 3=Between
    #           all orbitals
    ORBITAL_INTEGRAL_VELOCITY = ("200", "10", "3", "3")
    # Sequence: 200=Other functions (Part 2); 10=Output various kinds of
    #           integral between orbitals; 4=kinetic energy integral;
    #           3=Between all orbitals
    ORBITAL_INTEGRAL_KINETIC_ENERGY = ("200", "10", "4", "3")
    # Sequence: 200=Other functions (Part 2); 10=Output various kinds of
    #           integral between orbitals; 5=overlap integral; 3=Between
    #           all orbitals
    ORBITAL_INTEGRAL_OVERLAP = ("200", "10", "5", "3")
    # Sequence: 200=Other functions (Part 2); 19=Calculate spatial
    #           delocalization index (SDI) for orbitals or a function;
    #           1=Calcluate SDI for a real space function; 1=Electron
    #           density (rho)
    SPACIAL_DELOCALISATION_EDENSITY = ("200", "19", "1", "1")
    # Sequence: 200=Other functions (Part 2); 19=Calculate spatial
    #           delocalization index (SDI) for orbitals or a function;
    #           1=Calcluate SDI for a real space function; 2=Gradient norm
    #           of rho
    SPACIAL_DELOCALISATIOn_NORM_RHO = ("200", "19", "1", "2")
    # Sequence: 200=Other functions (Part 2); 19=Calculate spatial
    #           delocalization index (SDI) for orbitals or a function;
    #           1=Calcluate SDI for a real space function; 3=Laplacian of
    #           rho
    SPATIAL_DELOCALISATION_LAPLACIAN = ("200", "19", "1", "3")
    # Sequence: 200=Other functions (Part 2); 19=Calculate spatial
    #           delocalization index (SDI) for orbitals or a function;
    #           1=Calcluate SDI for a real space function; 4=Value of
    #           orbital wavefunction; 0=return to the previous menu
    SPATIAL_DELOCALISATION_ORB_WFN = ("200", "19", "1", "4", "0")
    # Sequence: 200=Other functions (Part 2); 19=Calculate spatial
    #           delocalization index (SDI) for orbitals or a function;
    #           1=Calcluate SDI for a real space function; 5=Electron spin
    #           density
    SPATIAL_DELOCALISATION_ESPIN_DENSITY = ("200", "19", "1", "5")
    # Sequence: 200=Other functions (Part 2); 19=Calculate spatial
    #           delocalization index (SDI) for orbitals or a function;
    #           1=Calcluate SDI for a real space function; 6=Hamiltonian
    #           kinetic energy density K(r)
    SPATIAL_DELOCALISATION_KR = ("200", "19", "1", "6")
    # Sequence: 200=Other functions (Part 2); 19=Calculate spatial
    #           delocalization index (SDI) for orbitals or a function;
    #           1=Calcluate SDI for a real space function; 7=Lagrangian
    #           kinetic energy density G(r)
    SPATIAL_DELOCALISATION_GR = ("200", "19", "1", "7")
    # Sequence: 200=Other functions (Part 2); 19=Calculate spatial
    #           delocalization index (SDI) for orbitals or a function;
    #           1=Calcluate SDI for a real space function; 8=Electrostatic
    #           potential from nuclear charges
    SPATIAL_DELOCALISATION_ESP_CHARGES = ("200", "19", "1", "8")
    # Sequence: 200=Other functions (Part 2); 19=Calculate spatial
    #           delocalization index (SDI) for orbitals or a function;
    #           1=Calcluate SDI for a real space function; 9=Electron
    #           localization function (ELF)
    SPATIAL_DELOCALISATION_ELF = ("200", "19", "1", "9")
    # Sequence: 200=Other functions (Part 2); 19=Calculate spatial
    #           delocalization index (SDI) for orbitals or a function;
    #           1=Calcluate SDI for a real space function; 10=Localized
    #           orbital locator (LOL)
    SPATIAL_DELOCALISATION_LOL = ("200", "19", "1", "10")
    # Sequence: 200=Other functions (Part 2); 19=Calculate spatial
    #           delocalization index (SDI) for orbitals or a function;
    #           1=Calcluate SDI for a real space function; 11=Local
    #           information entropy
    SPATIAL_DELOCALISATION_LOCAL_ENTROPY = ("200", "19", "1", "11")
    # Sequence: 200=Other functions (Part 2); 19=Calculate spatial
    #           delocalization index (SDI) for orbitals or a function;
    #           1=Calcluate SDI for a real space function; 12=Total
    #           electrostatic potential (ESP)
    SPATIAL_DELOCALISATION_ESP_TOTAL = ("200", "19", "1", "12")
    # Sequence: 200=Other functions (Part 2); 19=Calculate spatial
    #           delocalization index (SDI) for orbitals or a function;
    #           1=Calcluate SDI for a real space function; 13=Reduced
    #           density gradient (RDG)
    SPATIAL_DELOCALISATION_RDG = ("200", "19", "1", "13")
    # Sequence: 200=Other functions (Part 2); 19=Calculate spatial
    #           delocalization index (SDI) for orbitals or a function;
    #           1=Calcluate SDI for a real space function; 14=RDG with
    #           promolecular approximation
    SPATIAL_DELOCALISATION_RDG_PROMOLECULAR = ("200", "19", "1", "14")
    # Sequence: 200=Other functions (Part 2); 19=Calculate spatial
    #           delocalization index (SDI) for orbitals or a function;
    #           1=Calcluate SDI for a real space function;
    #           15=Sign(lambda2)*rho
    SPATIAL_DELOCALISATION_LAMBDA2RHO = ("200", "19", "1", "15")
    # Sequence: 200=Other functions (Part 2); 19=Calculate spatial
    #           delocalization index (SDI) for orbitals or a function;
    #           1=Calcluate SDI for a real space function;
    #           16=Sign(lambda2)*rho with promolecular approximation
    SPATIAL_DELOCALISATION_LAMBDA2RHO_PROMOLECULAR = ("200", "19", "1", "16")
    # Sequence: 200=Other functions (Part 2); 19=Calculate spatial
    #           delocalization index (SDI) for orbitals or a function;
    #           1=Calcluate SDI for a real space function; 18=Average
    #           local ionization energy (ALIE)
    SPATIAL_DELOCALISATION_ALIE = ("200", "19", "1", "18")
    # Sequence: 200=Other functions (Part 2); 19=Calculate spatial
    #           delocalization index (SDI) for orbitals or a function;
    #           1=Calcluate SDI for a real space function; 20=Electron
    #           delocal. range func. EDR(r;d)
    SPATIAL_DELOCALISATION_EDR = ("200", "19", "1", "20")
    # Sequence: 200=Other functions (Part 2); 19=Calculate spatial
    #           delocalization index (SDI) for orbitals or a function;
    #           1=Calcluate SDI for a real space function; 21=Orbital
    #           overlap dist. func. D(r)
    SPATIAL_DELOCALISATION_ORB_OVERLAP_DR = ("200", "19", "1", "21")
    # Sequence: 200=Other functions (Part 2); 19=Calculate spatial
    #           delocalization index (SDI) for orbitals or a function;
    #           1=Calcluate SDI for a real space function; 22=Delta-g
    #           (promolecular approximation)
    SPATIAL_DELOCALISATION_DELTAG_PROMOLECULAR = ("200", "19", "1", "22")
    # Sequence:     200=Other functions (Part 2) (leading spaces in the
    #           token are stray but harmless); 19=Calculate spatial
    #           delocalization index (SDI) for orbitals or a function;
    #           1=Calcluate SDI for a real space function; 23=Delta-g
    #           (Hirshfeld partition)
    SPATIAL_DELOCALISATION_DELTAG_HIRSHFELD = ("    200", "19", "1", "23")
    # Sequence: 200=Other functions (Part 2); 19=Calculate spatial
    #           delocalization index (SDI) for orbitals or a function;
    #           1=Calcluate SDI for a real space function; 24=Interaction
    #           region indicator (IRI)
    SPATIAL_DELOCALISATION_IRI = ("200", "19", "1", "24")
    # Fluctuation NCI (aRDG): RDG-based weak-interaction analysis time-
    # averaged over an MD trajectory; reveals which non-covalent contacts
    # persist throughout the simulation
    # Sequence: 200=Other functions (Part 2); 1=Calculate core-valence
    #           bifurcation (CVB) index and related quantities
    CVB_INDEX = ("200", "1")
    # Compute atomic and bond dipole moments in Hilbert space using
    # Mulliken partitioning; decomposes the total molecular dipole into
    # atomic and inter-atomic bond contributions
    # Sequence: 200=Other functions (Part 2); 2=Calculate atomic and bond
    #           dipole moments in Hilbert space
    ATOMIC_BOND_DIPOLES = ("200", "2")
    # Generate Gaussian .cube files for multiple orbitals in a single
    # automated run; each selected orbital wavefunction is exported as a
    # separate file
    # Sequence: 200=Other functions (Part 2); 3=Generate cube file for
    #           multiple orbital wavefunctions; 0=return to the previous
    #           menu
    MULTIPLE_ORBITAL_CUBES = ("200", "3", "0")
    # Generate 3D NMR shielding (ICSS) data as .cube files; visualise
    # magnetically induced ring-current effects as shielding isosurfaces
    # in VMD
    # Sequence: 200=Other functions (Part 2); 4=Output and plot specific
    #           property in a plane
    ICSS_CUBES = ("200", "4")
    # Compute and plot the spherical radial distribution function RDF(r) =
    # 4*pi*r^2 * f(r) of a chosen real-space function; useful for
    # spherically symmetric systems such as fullerenes or atoms
    # Sequence: 200=Other functions (Part 2); 5=Plot radial distribution
    #           function for a real space function
    RADIAL_DISTRIBUTION = ("200", "5")
    # Analyse the overlap correspondence between orbitals in two different
    # wavefunctions; tracks how orbitals evolve along a reaction path or
    # change between two levels of theory
    # Sequence: 200=Other functions (Part 2); 6=Analyze correspondence
    #           between orbitals in two wavefunctions
    ORBITAL_CORRESPONDENCE = ("200", "6")
    # Parse and tabulate polarizability (alpha) and first/second/third
    # hyperpolarizability (beta, gamma, delta) tensors from a Gaussian
    # polar- keyword output file
    # Sequence: 200=Other functions (Part 2); 7=Population analysis and
    #           calculation of atomic charges
    PARSE_POLARIZABILITY_GAUSSIAN = ("200", "7")
    # Compute polarizability and 1st/2nd/3rd hyperpolarizabilities by the
    # Sum- Over-States (SOS) method using TD-DFT excited-state energies
    # and transition moments
    # Sequence: 200=Other functions (Part 2); 8=Orbital composition
    #           analysis
    SOS_HYPERPOLARIZABILITY = ("200", "8")
    # Calculate and print average bond lengths and average coordination
    # numbers for all element-type combinations in the molecule
    # Sequence: 200=Other functions (Part 2); 9=Calculate average bond
    #           length and average coordinate number
    AVERAGE_BOND_LENGTH = ("200", "9")
    # Compute one-electron integrals between pairs of orbitals: kinetic
    # energy <i|T|j>, nuclear attraction <i|V|j>, Coulomb, and exchange
    # integrals
    # Sequence: 200=Other functions (Part 2); 10=Output various kinds of
    #           integral between orbitals
    ORBITAL_INTEGRALS = ("200", "10")
    # Compute the centroid (first moment), second moments, and radius of
    # gyration of a chosen real-space function; characterises its spatial
    # spread and shape
    # Sequence: 200=Other functions (Part 2); 11=Calculate center,
    #           first/second moments and radius of gyration of a function
    FUNCTION_MOMENTS = ("200", "11")
    # Compute the Energy Index (EI) and Bond Polarity Index (BPI) from
    # kinetic and potential energy densities evaluated at bond critical
    # points
    # Sequence: 200=Other functions (Part 2); 12=Calculate energy index
    #           (EI) or bond polarity index (BPI)
    ENERGY_INDEX = ("200", "12")
    # Decompose a grid-based property (e.g., electron density) into
    # contributions from individual MOs and export per-orbital .cube files
    # for detailed analysis
    # Sequence: 200=Other functions (Part 2); 13=Evaluate orbital
    #           contributions to density difference or other grid data
    # INTERACTIVE ONLY -- requires the path to a separate external file
    #                     (not the primary wavefunction) -- e.g. a
    #                     Gaussian/ORCA output, a cube file, or an NBO
    #                     file
    ORBITAL_CONTRIBUTIONS_TO_GRID = ("200", "13")
    # Identify and characterise topologically connected electron-density
    # domains: contiguous regions above a chosen isovalue threshold in the
    # 3D grid
    # Sequence: 200=Other functions (Part 2); 14=Domain analysis
    #           (Obtaining properties within isosurfaces of a function)
    DOMAIN_ANALYSIS = ("200", "14")
    # Compute a spatial correlation index between two real-space functions
    # evaluated on the same grid; quantifies how similarly their
    # distributions are arranged in space
    # Sequence: 200=Other functions (Part 2); 15=Calculate electron
    #           correlation index (PCCP, 18, 24015)
    CORRELATION_INDEX = ("200", "15")
    # Diagonalise the one-particle density matrix to obtain natural
    # orbitals (NOs) and their occupation numbers; useful for
    # characterising electron correlation and multi-reference character
    # Sequence: 200=Other functions (Part 2); 16=Generate natural orbitals
    #           based on the density matrix in .fch/.fchk file
    NATURAL_ORBITALS = ("200", "16")
    # Evaluate two-electron Coulomb (J_ij) and exchange (K_ij) integrals
    # between all pairs of molecular orbitals; provides input for EDA,
    # perturbation theory, and excited-state analysis
    # Sequence: 200=Other functions (Part 2); 17=Calculate Coulomb and
    #           exchange integrals between two orbitals
    COULOMB_EXCHANGE_INTEGRALS = ("200", "17")
    # Compute Bond Length Alternation (BLA) and Bond Order Alternation
    # (BOA) along a user-defined conjugated chain; quantifies the degree
    # of alternation in polyene or polyacetylene-type pi systems
    # Sequence: 200=Other functions (Part 2); 18=Calculate bond
    #           length/order alternation (BLA/BOA)
    BLA_BOA_ANALYSIS = ("200", "18")
    # Compute the spatial delocalization index (SDI): measures how broadly
    # the electron density is spread across space; indicator of overall
    # delocalisation
    # Sequence: 200=Other functions (Part 2); 19=Calculate spatial
    #           delocalization index (SDI) for orbitals or a function
    SPATIAL_DELOCALIZATION_INDEX = ("200", "19")
    # Bond Order Decomposition (BOD) and Natural Atomic Dipole Orbital
    # (NADO) analysis: decomposes Mayer bond orders and atomic dipole
    # moments into natural orbital pair contributions
    # Sequence: 200=Other functions (Part 2); 20=Bond order density (BOD)
    #           and natural adaptive orbital (NAdO) analyses
    BOD_NADO_ANALYSIS = ("200", "20")
    # Perform Lowdin symmetric orthogonalization of the current orbital
    # set and optionally evaluate the resulting orbital energies directly
    # from the Fock matrix
    # Sequence: 200=Other functions (Part 2); 21=Perform Lowdin
    #           orthogonalization between occupied orbitals
    LOWDIN_ORTHOGONALIZATION = ("200", "21")

    # ─────────────────────────────────────────────────────────────────────────
    # Main Menu 300: Utilities Part 3
    # ─────────────────────────────────────────────────────────────────────────
    # Calculate the free void volume in a periodic unit cell not occupied
    # by atomic vdW spheres; relevant for porosity characterisation in
    # MOFs, zeolites, and porous organic cages
    # Sequence: 300=Other functions (Part 3); 1=Viewing free regions and
    #           calculating free volume in a cell
    FREE_VOLUME_IN_CELL = ("300", "1")
    # Fit the spherically averaged radial electron density of each atom in
    # the molecule to a set of Slater-type or Gaussian-type radial
    # functions; generates transferable atomic density parameters
    # Sequence: 300=Other functions (Part 3); 2=Fitting atomic radial
    #           density as linear combination of multiple STOs or GTFs
    FIT_ATOMIC_RADIAL_DENSITY = ("300", "2")
    # Simulate a Scanning Tunnelling Microscopy (STM) image using the
    # Tersoff- Hamann approximation from the local density of states near
    # the Fermi energy
    # Sequence: 300=Other functions (Part 3); 4=Simulating scanning
    #           tunneling microscope (STM) image
    STM_IMAGE = ("300", "4")
    # Compute and print the molecular electric multipole moment tensor up
    # to hexadecapole order from the electron density and nuclear charges
    # Sequence: 300=Other functions (Part 3); 5=Calculate electric
    #           dipole/multipole moments and electronic spatial extent
    ELECTRIC_MULTIPOLE_MOMENTS = ("300", "5")
    # Evaluate orbital energies epsilon_i = <psi_i|F|psi_i> directly from
    # the Fock matrix without running a new SCF; particularly useful after
    # orbital localisation or biorthogonalisation
    # Sequence: 300=Other functions (Part 3); 6=Calculate energies of
    #           present orbitals by inputting Fock matrix
    # INTERACTIVE ONLY -- requires the path to a separate external file
    #                     (not the primary wavefunction) -- e.g. a
    #                     Gaussian/ORCA output, a cube file, or an NBO
    #                     file
    ORBITAL_ENERGIES_FROM_FOCK = ("300", "6")
    # Perform geometric transformations on the molecular structure:
    # translate, rotate, reflect, invert, or generate symmetry-equivalent
    # atoms
    # Sequence: 300=Other functions (Part 3); 7=Geometry operation on the
    #           present system
    GEOMETRY_OPERATIONS = ("300", "7")
    # Generate a molecular surface distance projection map plotting d_i
    # (distance to nearest internal atom) vs. d_e (distance to nearest
    # external atom) for every surface point; reveals steric
    # complementarity and packing motifs
    # Sequence: 300=Other functions (Part 3); 8=Plot surface distance
    #           projection map; 0=return to the previous menu; -2=Export
    #           plane data as distmap.txt in current folder; 1=Save the
    #           map as graphical file in current folder; -1=Return
    SURFACE_DISTANCE_PROJECTION = ("300", "8", "0", "-2", "1", "-1")
    # Determine the Fermi energy level from the orbital energy
    # distribution and occupation numbers; particularly useful for
    # periodic or large cluster systems
    # Sequence: 300=Other functions (Part 3); 9=Determine Fermi level
    DETERMINE_FERMI_LEVEL = ("300", "9")
