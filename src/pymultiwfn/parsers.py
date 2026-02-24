"""Output parsers for Multiwfn results."""

import re
from typing import Any

# Shared float pattern for all numeric parsing
FLOAT_PATTERN = r"[-+]?(?:\d+\.\d*|\d*\.\d+|\d+)(?:[Ee][-+]?\d+)?"


class OutputParser:
    """Base class for output parsers."""

    pass

#TODO(fs): maybe worth looking into making this code less long? I think i went
# a bit overkill with all the different parsign methods and regex. Maybe some 
# of the methods can be combined?
# =============================================================================
# Menu 7: Population analysis & atomic charges
# =============================================================================


class ChargeParser(OutputParser):
    """Parser for atomic charges.

    Handles output from all Menu 7 charge methods including Hirshfeld,
    VDD, Mulliken, Lowdin, SCPA, Stout-Politzer, Bickelhaupt, Becke,
    ADCH, CHELPG, MK, AIM, Hirshfeld-I, CM5, EEM, RESP, Gasteiger,
    MBIS, and DDEC.
    """

    @staticmethod
    def parse(stdout: str, method: str = "Hirshfeld") -> dict[int, float]:
        """Extract atomic charges from Multiwfn output.

        Parameters
        ----------
        stdout : str
            Multiwfn standard output
        method : str
            Charge method name (e.g., "Hirshfeld", "ADCH", "RESP",
            "Mulliken", "Lowdin", "CHELPG", "MK", "CM5", "MBIS",
            "DDEC", "Gasteiger", "EEM", "Becke", "VDD", "SCPA",
            "Stout-Politzer", "Bickelhaupt", "AIM", "Hirshfeld-I")

        Returns
        -------
        dict[int, float]
            Dictionary mapping atom indices to charges
        """
        charges: dict[int, float] = {}

        # Pattern 1: "Hirshfeld charge of atom     1(C ) is  0.03208687"
        pattern1 = (
            rf"{re.escape(method)}\s+charge of atom\s+(\d+)\s*\([A-Za-z]+"
            rf"\s*\)\s+is\s+({FLOAT_PATTERN})"
        )
        # Pattern 2: "Atom    1(C ):     0.03209323" (Final atomic charges)
        pattern2 = rf"Atom\s+(\d+)\s*\([A-Za-z]+\s*\)\s*:\s*({FLOAT_PATTERN})"
        # Pattern 3: "    1(C )   -0.0523" (summary table format)
        pattern3 = rf"^\s*(\d+)\s*\([A-Za-z]+\s*\)\s+({FLOAT_PATTERN})"
        # Pattern 4: "   1  C      -0.052300" (column format)
        pattern4 = rf"^\s*(\d+)\s+[A-Za-z]+\s+({FLOAT_PATTERN})"
        # Pattern 5: "Population of atom  1(C ) :   5.9679" (Mulliken/Lowdin)
        pattern5 = (
            rf"Population of atom\s+(\d+)\s*\([A-Za-z]+\s*\)\s*:\s+"
            rf"({FLOAT_PATTERN})"
        )
        # Pattern 6: "Charge of atom  1(C ) :  0.0321" (generic)
        pattern6 = (
            rf"Charge of atom\s+(\d+)\s*\([A-Za-z]+\s*\)\s*:\s+"
            rf"({FLOAT_PATTERN})"
        )
        # Pattern 7: "Net charge:  0.0321" after atom header
        # (some methods print atom header then net charge on next line)

        in_charge_section = False
        in_final_section = False

        for line in stdout.split("\n"):
            line_lower = line.lower()

            # Check for section markers
            #TODO(fs): lots of if statements here, worth refactoring i think?
            if "final atomic charges" in line_lower:
                in_final_section = True
                in_charge_section = True
                charges.clear()  # Prefer final charges
                continue
            elif method.lower() in line_lower and (
                "charge" in line_lower or "population" in line_lower
            ):
                in_charge_section = True
                continue
            elif "summary of" in line_lower and "charge" in line_lower:
                in_charge_section = True
                charges.clear()
                continue

            # Try pattern 1 first (most specific - includes method name)
            if match := re.search(pattern1, line, re.IGNORECASE):
                charges[int(match[1])] = float(match[2])
                continue

            # Try pattern 5 (population format for Mulliken/Lowdin)
            if match := re.search(pattern5, line):
                if in_charge_section:
                    charges[int(match[1])] = float(match[2])
                continue

            # Try pattern 6 (generic charge format)
            if match := re.search(pattern6, line):
                if in_charge_section:
                    charges[int(match[1])] = float(match[2])
                continue

            # Try pattern 2 (Final atomic charges format)
            if match := re.search(pattern2, line):
                if in_charge_section or in_final_section:
                    charges[int(match[1])] = float(match[2])
                continue

            # Try pattern 3 (table format)
            if in_charge_section and (match := re.match(pattern3, line)):
                charges[int(match[1])] = float(match[2])
                continue

            # Try pattern 4 (column format)
            if in_charge_section and (match := re.match(pattern4, line)):
                charges[int(match[1])] = float(match[2])
                continue

            # End of section detection
            if (
                in_final_section
                and charges
                and "calculation took" in line_lower
            ):
                break

        return charges

    @staticmethod
    def parse_dipole(stdout: str) -> dict[str, float] | None:
        """Extract molecular dipole moment from charge output.

        Returns
        -------
        dict[str, float] | None
            Dictionary with 'x', 'y', 'z', 'total' dipole components
            in Debye, or None if not found.
        """
        pattern = (
            rf"Dipole moment.*?X=\s*({FLOAT_PATTERN})\s+"
            rf"Y=\s*({FLOAT_PATTERN})\s+Z=\s*({FLOAT_PATTERN})\s+"
            rf"Tot=\s*({FLOAT_PATTERN})"
        )
        if match := re.search(pattern, stdout, re.IGNORECASE):
            return {
                "x": float(match[1]),
                "y": float(match[2]),
                "z": float(match[3]),
                "total": float(match[4]),
            }
        return None


# =============================================================================
# Menu 8: Orbital composition analysis
# =============================================================================


class OrbitalCompositionParser(OutputParser):
    """Parser for orbital composition analysis (Menu 8).

    Handles Mulliken, SCPA, Stout-Politzer, NAO, Hirshfeld, and Becke
    orbital composition methods, as well as fragment composition and
    LOBA oxidation state analysis.
    """

    @staticmethod
    def parse(stdout: str) -> list[dict[str, Any]]:
        """Extract orbital composition data.

        Returns
        -------
        list[dict[str, Any]]
            List of dicts with keys: 'orbital', 'energy', 'occupation',
            'contributions' (dict mapping atom/fragment labels to percentages)
        """
        orbitals: list[dict[str, Any]] = []

        # "Orbital    5  Occ= 2.000000  E= -0.72340 a.u."
        orb_pattern = (
            rf"Orbital\s+(\d+)\s+Occ=\s*({FLOAT_PATTERN})\s+"
            rf"E=\s*({FLOAT_PATTERN})"
        )
        # "  C  1    s :   23.45%"  or  "  1(C )    45.23%"
        contrib_pattern = (
            rf"([A-Za-z]+)\s+(\d+)\s+.*?:\s+({FLOAT_PATTERN})\s*%"
        )
        contrib_pattern2 = rf"(\d+)\s*\([A-Za-z]+\s*\)\s+({FLOAT_PATTERN})\s*%"

        current_orb: dict[str, Any] | None = None
        for line in stdout.split("\n"):
            if match := re.search(orb_pattern, line):
                if current_orb is not None:
                    orbitals.append(current_orb)
                current_orb = {
                    "orbital": int(match[1]),
                    "occupation": float(match[2]),
                    "energy": float(match[3]),
                    "contributions": {},
                }
            elif current_orb is not None:
                if match := re.search(contrib_pattern, line):
                    label = f"{match[1]}{match[2]}"
                    current_orb["contributions"][label] = float(match[3])
                elif match := re.search(contrib_pattern2, line):
                    label = f"atom_{match[1]}"
                    current_orb["contributions"][label] = float(match[2])

        if current_orb is not None:
            orbitals.append(current_orb)
        return orbitals

    @staticmethod
    def parse_oxidation_states(stdout: str) -> dict[int, int]:
        """Extract LOBA oxidation states.

        Returns
        -------
        dict[int, int]
            Mapping of atom index to formal oxidation state
        """
        pattern = (
            rf"Atom\s+(\d+)\s*\([A-Za-z]+\s*\).*?"
            rf"oxidation state.*?({FLOAT_PATTERN})"
        )
        states: dict[int, int] = {
            int(match[1]): int(round(float(match[2])))
            for match in re.finditer(pattern, stdout, re.IGNORECASE)
        }
        return states


# =============================================================================
# Menu 9: Bond order analysis
# =============================================================================


class BondOrderParser(OutputParser):
    """Parser for bond orders.

    Handles Mayer, Wiberg, Mulliken, fuzzy, Laplacian, IBSI, and
    multicenter bond orders, as well as AV1245 index and bond order
    decompositions.
    """

    @staticmethod
    def parse(stdout: str, **kwargs: Any) -> dict[tuple[int, int], float]:
        """Extract bond orders from Multiwfn output.

        Parameters
        ----------
        stdout : str
            Multiwfn standard output

        Returns
        -------
        dict[tuple[int, int], float]
            Dictionary mapping atom pairs to bond orders
        """
        bond_orders: dict[tuple[int, int], float] = {}

        # Pattern 1: "#    1:         1(C )    2(N )    1.95394965"
        pattern1 = (
            rf"#\s*\d+:\s+(\d+)\s*\([^)]+\)\s+(\d+)\s*\([^)]+\)\s+"
            rf"({FLOAT_PATTERN})"
        )

        # Pattern 2: "1(C ) -    2(C ):  1.4523"
        pattern2 = (
            rf"(\d+)\s*\([^)]+\)\s*-\s*(\d+)\s*\([^)]+\)\s*:\s*"
            rf"({FLOAT_PATTERN})"
        )

        # Pattern 3: Simple "1 - 2   1.4523"
        pattern3 = rf"^\s*(\d+)\s*-\s*(\d+)\s+({FLOAT_PATTERN})"

        # Pattern 4: "   1    2    1.4523" (column format)
        pattern4 = rf"^\s*(\d+)\s+(\d+)\s+({FLOAT_PATTERN})\s*$"

        for line in stdout.split("\n"):
            match = None

            if (
                (match := re.search(pattern1, line))
                or (match := re.search(pattern2, line))
                or (match := re.match(pattern3, line))
                or (match := re.match(pattern4, line))
            ):
                pass

            if match:
                atom1 = int(match[1])
                atom2 = int(match[2])
                bo = float(match[3])
                if atom1 > atom2:
                    atom1, atom2 = atom2, atom1
                bond_orders[(atom1, atom2)] = bo

        return bond_orders

    @staticmethod
    def parse_valence(stdout: str) -> dict[int, dict[str, float]]:
        """Extract total valence and free valence for each atom.

        Returns
        -------
        dict[int, dict[str, float]]
            Mapping of atom index to {'total_valence', 'free_valence'}
        """
        valences: dict[int, dict[str, float]] = {}

        # "Total valence of atom    1(C ):   3.9412"
        total_pattern = (
            rf"Total valence of atom\s+(\d+)\s*\([^)]+\)\s*:\s+"
            rf"({FLOAT_PATTERN})"
        )
        # "Free valence of atom     1(C ):   0.0588"
        free_pattern = (
            rf"Free valence of atom\s+(\d+)\s*\([^)]+\)\s*:\s+"
            rf"({FLOAT_PATTERN})"
        )

        for match in re.finditer(total_pattern, stdout):
            idx = int(match[1])
            valences.setdefault(idx, {})["total_valence"] = float(match[2])
        for match in re.finditer(free_pattern, stdout):
            idx = int(match[1])
            valences.setdefault(idx, {})["free_valence"] = float(match[2])

        return valences

    @staticmethod
    def parse_multicenter(stdout: str) -> list[dict[str, Any]]:
        """Extract multicenter bond order results.

        Returns
        -------
        list[dict[str, Any]]
            List of dicts with 'atoms' (list[int]) and 'bond_order' (float)
        """
        results: list[dict[str, Any]] = []
        # "Multi-center bond order of atoms  1  2  3 :  0.12345"
        pattern = (
            rf"Multi-center bond order of atoms\s+([\d\s]+):\s+"
            rf"({FLOAT_PATTERN})"
        )
        for match in re.finditer(pattern, stdout, re.IGNORECASE):
            atoms = [int(x) for x in match[1].split()]
            results.append({"atoms": atoms, "bond_order": float(match[2])})
        return results

    @staticmethod
    def parse_decomposition(
        stdout: str,
    ) -> list[dict[str, Any]]:
        """Extract per-orbital bond order decomposition.

        Returns
        -------
        list[dict[str, Any]]
            List of dicts with 'orbital', 'contribution' keys
        """
        decomp: list[dict[str, Any]] = []
        # "Orbital   5:   0.23456"
        pattern = rf"Orbital\s+(\d+)\s*:\s+({FLOAT_PATTERN})"
        decomp.extend(
            {
                "orbital": int(match[1]),
                "contribution": float(match[2]),
            }
            for match in re.finditer(pattern, stdout)
        )
        return decomp


# =============================================================================
# Menu 2: Topology analysis
# =============================================================================


class CriticalPointParser(OutputParser):
    """Parser for critical point information from topology analysis.

    Handles output from TOPOLOGY_SEARCH_CPS, TOPOLOGY_GENERATE_PATHS,
    TOPOLOGY_ANALYSIS_COMPLETE, and specialised topology analyses
    (ESP, LOL, ELF, Laplacian, BCP/RCP/CCP searches).
    """

    CP_TYPE_NAMES: dict[str, str] = {
        "(3,-3)": "nuclear",
        "(3,-1)": "bond",
        "(3,+1)": "ring",
        "(3,+3)": "cage",
    }

    @staticmethod
    def parse(stdout: str, **kwargs: Any) -> list[dict[str, Any]]:
        """Extract critical point information from topology analysis.

        Parameters
        ----------
        stdout : str
            Multiwfn standard output

        Returns
        -------
        list[dict[str, Any]]
            List of dicts with 'index', 'type', 'cp_type', 'position',
            and optionally 'rho', 'laplacian', 'ellipticity'
        """
        cps: list[dict[str, Any]] = []

        pattern = r"CP\s+(\d+)\s+\((\d+),([+-]?\d+)\)"
        pos_pattern = (
            rf"Position.*?:\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})"
            rf"\s+({FLOAT_PATTERN})"
        )
        rho_pattern = rf"Density of.*?:\s+({FLOAT_PATTERN})"
        lap_pattern = rf"Laplacian.*?:\s+({FLOAT_PATTERN})"
        ell_pattern = rf"Ellipticity.*?:\s+({FLOAT_PATTERN})"

        lines = stdout.split("\n")
        for i, line in enumerate(lines):
            if match := re.search(pattern, line):
                cp_index = int(match[1])
                cp_type = f"({match[2]},{match[3]})"

                cp: dict[str, Any] = {
                    "index": cp_index,
                    "type": cp_type,
                    "cp_type": CriticalPointParser.CP_TYPE_NAMES.get(
                        cp_type, "unknown"
                    ),
                    "position": None,
                    "rho": None,
                    "laplacian": None,
                    "ellipticity": None,
                }

                # Scan following lines for properties
                for j in range(i, min(i + 10, len(lines))):
                    sub = lines[j]
                    if pos_match := re.search(pos_pattern, sub):
                        cp["position"] = (
                            float(pos_match[1]),
                            float(pos_match[2]),
                            float(pos_match[3]),
                        )
                    if rho_match := re.search(rho_pattern, sub):
                        cp["rho"] = float(rho_match[1])
                    if lap_match := re.search(lap_pattern, sub):
                        cp["laplacian"] = float(lap_match[1])
                    if ell_match := re.search(ell_pattern, sub):
                        cp["ellipticity"] = float(ell_match[1])

                cps.append(cp)

        return cps

    @staticmethod
    def parse_bond_paths(
        stdout: str,
    ) -> list[dict[str, Any]]:
        """Extract bond path information.

        Returns
        -------
        list[dict[str, Any]]
            List of dicts with 'atom1', 'atom2', 'bcp_index', 'path_length'
        """
        paths: list[dict[str, Any]] = []
        # "Bond path between atom  1(C ) and atom  2(N ), BCP  3, length  2.456"
        pattern = (
            rf"Bond path between atom\s+(\d+).*?and atom\s+(\d+).*?"
            rf"BCP\s+(\d+).*?length\s+({FLOAT_PATTERN})"
        )
        paths.extend(
            {
                "atom1": int(match[1]),
                "atom2": int(match[2]),
                "bcp_index": int(match[3]),
                "path_length": float(match[4]),
            }
            for match in re.finditer(pattern, stdout, re.IGNORECASE)
        )
        return paths

    @staticmethod
    def summary(cps: list[dict[str, Any]]) -> dict[str, int]:
        """Count critical points by type.

        Returns
        -------
        dict[str, int]
            Counts per CP type name
        """
        counts: dict[str, int] = {}
        for cp in cps:
            name = cp.get("cp_type", "unknown")
            counts[name] = counts.get(name, 0) + 1
        return counts


# =============================================================================
# Menu 10: Density of states
# =============================================================================


class DOSParser(OutputParser):
    """Parser for density of states output.

    Handles TDOS, PDOS, OPDOS, LDOS, COHP, and photoelectron spectrum.
    """

    @staticmethod
    def parse(stdout: str) -> dict[str, list[float]]:
        """Extract DOS curve data.

        Returns
        -------
        dict[str, list[float]]
            Dictionary with 'energies' and 'dos' lists, plus optional
            fragment-resolved 'pdos_<n>' lists.
        """
        data: dict[str, list[float]] = {"energies": [], "dos": []}

        # Two-column: "  -10.2345    0.5678"
        pattern2 = rf"^\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s*$"
        # Multi-column for PDOS: " -10.2345  0.567  0.234  0.100"
        pattern_multi = rf"^\s+({FLOAT_PATTERN})((?:\s+{FLOAT_PATTERN})+)\s*$"

        in_data = False
        for line in stdout.split("\n"):
            if "TDOS" in line or "PDOS" in line or "OPDOS" in line:
                in_data = True
                continue
            if not in_data:
                continue

            if match := re.match(pattern_multi, line):
                energy = float(match[1])
                vals = [float(v) for v in match[2].split()]
                data["energies"].append(energy)
                if vals:
                    data["dos"].append(vals[0])
                for k, v in enumerate(vals[1:], start=1):
                    key = f"pdos_{k}"
                    data.setdefault(key, []).append(v)
            elif match := re.match(pattern2, line):
                data["energies"].append(float(match[1]))
                data["dos"].append(float(match[2]))

        return data

    @staticmethod
    def parse_orbital_energies(stdout: str) -> list[dict[str, Any]]:
        """Extract orbital energies used in DOS.

        Returns
        -------
        list[dict[str, Any]]
            List of dicts with 'index', 'energy_eV', 'occupation'
        """
        pattern = rf"(\d+)\s+({FLOAT_PATTERN})\s+eV\s+Occ=\s*({FLOAT_PATTERN})"
        orbitals: list[dict[str, Any]] = [
            {
                "index": int(match[1]),
                "energy_eV": float(match[2]),
                "occupation": float(match[3]),
            }
            for match in re.finditer(pattern, stdout)
        ]
        return orbitals


# =============================================================================
# Menu 11: Spectra simulation
# =============================================================================


class SpectrumParser(OutputParser):
    """Parser for spectrum data.

    Handles IR, Raman, UV-Vis, ECD, VCD, ROA, NMR, fluorescence, PVS,
    and colour prediction output.
    """

    @staticmethod
    def parse(stdout: str, **kwargs: Any) -> dict[str, list[float]]:
        """Extract spectrum data (frequencies/wavelengths, intensities).

        Parameters
        ----------
        stdout : str
            Multiwfn standard output

        Returns
        -------
        dict[str, list[float]]
            Dictionary with 'frequencies' and 'intensities' lists.
            For UV-Vis/ECD the keys are 'wavelengths' and
            'intensities'.
        """
        spectrum: dict[str, list[float]] = {
            "frequencies": [],
            "intensities": [],
        }

        # IR/Raman: "  1234.56 cm^-1  Intensity:  45.678"
        pattern1 = (
            rf"({FLOAT_PATTERN})\s+cm\^?-1.*?Intensity:\s+({FLOAT_PATTERN})"
        )
        # UV-Vis/ECD: "  345.67 nm  f= 0.1234"
        pattern_uv = (
            rf"({FLOAT_PATTERN})\s+nm.*?(?:f=|Str[.=])\s*({FLOAT_PATTERN})"
        )
        # NMR: "  Atom  1(C )  shift:  123.45 ppm"
        pattern_nmr = rf"Atom\s+(\d+)\s*\([^)]+\)\s+shift:\s+({FLOAT_PATTERN})"
        # Generic two-column
        pattern2 = rf"^\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s*$"

        # Try IR/Raman first
        for match in re.finditer(pattern1, stdout):
            spectrum["frequencies"].append(float(match[1]))
            spectrum["intensities"].append(float(match[2]))

        if spectrum["frequencies"]:
            return spectrum

        # Try UV-Vis/ECD
        uv_data: dict[str, list[float]] = {
            "wavelengths": [],
            "intensities": [],
        }
        for match in re.finditer(pattern_uv, stdout):
            uv_data["wavelengths"].append(float(match[1]))
            uv_data["intensities"].append(float(match[2]))
        if uv_data["wavelengths"]:
            return uv_data

        # Try NMR
        nmr_data: dict[str, list[float]] = {
            "atom_indices": [],
            "chemical_shifts": [],
        }
        for match in re.finditer(pattern_nmr, stdout):
            nmr_data["atom_indices"].append(float(match[1]))
            nmr_data["chemical_shifts"].append(float(match[2]))
        if nmr_data["atom_indices"]:
            return nmr_data

        # Fallback: generic two-column
        for line in stdout.split("\n"):
            if match2 := re.match(pattern2, line):
                spectrum["frequencies"].append(float(match2[1]))
                spectrum["intensities"].append(float(match2[2]))

        return spectrum

    @staticmethod
    def parse_transitions(stdout: str) -> list[dict[str, Any]]:
        """Extract discrete transition data (excitation energies, strengths).

        Returns
        -------
        list[dict[str, Any]]
            List of dicts with 'state', 'energy_eV', 'wavelength_nm',
            'osc_strength', and optionally 'rot_strength'
        """
        transitions: list[dict[str, Any]] = []
        # "Excited state   1:  E= 3.4567 eV  lam= 358.7 nm  f= 0.0123"
        pattern = (
            rf"Excited state\s+(\d+).*?E=\s*({FLOAT_PATTERN})\s*eV.*?"
            rf"lam=\s*({FLOAT_PATTERN})\s*nm.*?"
            rf"f=\s*({FLOAT_PATTERN})"
        )
        rot_pattern = rf"R.*?=\s*({FLOAT_PATTERN})"

        lines = stdout.split("\n")
        for i, line in enumerate(lines):
            if match := re.search(pattern, line):
                t: dict[str, Any] = {
                    "state": int(match[1]),
                    "energy_eV": float(match[2]),
                    "wavelength_nm": float(match[3]),
                    "osc_strength": float(match[4]),
                }
                # Check for rotatory strength on same or next line
                combined = line
                if i + 1 < len(lines):
                    combined += lines[i + 1]
                if rot_match := re.search(rot_pattern, combined):
                    t["rot_strength"] = float(rot_match[1])
                transitions.append(t)
        return transitions

    @staticmethod
    def parse_color(stdout: str) -> dict[str, Any] | None:
        """Extract predicted colour from PREDICT_COLOR output.

        Returns
        -------
        dict[str, Any] | None
            Dictionary with CIE coordinates and RGB values
        """
        result: dict[str, Any] = {}
        for pat_name, pat in [
            ("X", rf"X=\s*({FLOAT_PATTERN})"),
            ("Y", rf"Y=\s*({FLOAT_PATTERN})"),
            ("Z", rf"Z=\s*({FLOAT_PATTERN})"),
            ("R", r"R=\s*(\d+)"),
            ("G", r"G=\s*(\d+)"),
            ("B", r"B=\s*(\d+)"),
        ]:
            if m := re.search(pat, stdout):
                result[pat_name] = (
                    int(m[1]) if pat_name in ("R", "G", "B") else float(m[1])
                )
        return result or None


# =============================================================================
# Menu 12: Quantitative molecular surface analysis
# =============================================================================


class SurfaceParser(OutputParser):
    """Parser for molecular surface analysis output.

    Handles ESP surface analysis, ALIE surface analysis, surface area/
    volume, Becke surface, Hirshfeld surface, surface extrema, and
    Hirshfeld surface fingerprint.
    """

    @staticmethod
    def parse(stdout: str) -> dict[str, Any]:
        """Extract surface analysis statistics.

        Returns
        -------
        dict[str, Any]
            Dictionary with surface descriptors: 'area', 'volume',
            'V_S_plus', 'V_S_minus', 'sigma2_total', 'nu', 'Pi',
            'n_extrema_max', 'n_extrema_min', and lists of extrema.
        """
        result: dict[str, Any] = {}

        patterns: dict[str, str] = {
            "area": rf"(?:Overall|Molecular)\s+surface\s+area.*?:\s+({FLOAT_PATTERN})",
            "volume": rf"(?:Enclosed|Molecular)\s+volume.*?:\s+({FLOAT_PATTERN})",
            "V_S_plus": rf"V_S\+.*?:\s+({FLOAT_PATTERN})",
            "V_S_minus": rf"V_S-.*?:\s+({FLOAT_PATTERN})",
            "sigma2_total": rf"sigma\^?2_?tot.*?:\s+({FLOAT_PATTERN})",
            "nu": rf"[Nn]u.*?:\s+({FLOAT_PATTERN})",
            "Pi": rf"Pi\s*:\s+({FLOAT_PATTERN})",
            "V_S_max": rf"Global.*?maximum.*?:\s+({FLOAT_PATTERN})",
            "V_S_min": rf"Global.*?minimum.*?:\s+({FLOAT_PATTERN})",
            "balance": rf"Balance.*?:\s+({FLOAT_PATTERN})",
        }

        for key, pat in patterns.items():
            if match := re.search(pat, stdout, re.IGNORECASE):
                result[key] = float(match[1])

        return result

    @staticmethod
    def parse_extrema(stdout: str) -> list[dict[str, Any]]:
        """Extract surface extrema (minima and maxima).

        Returns
        -------
        list[dict[str, Any]]
            List of dicts with 'type' ('min'/'max'), 'value',
            'position' (x, y, z)
        """
        extrema: list[dict[str, Any]] = []
        # "Local  minimum   1:   -45.678  at   1.234   2.345   3.456"
        pattern = (
            rf"Local\s+(min\w*|max\w*)\s+(\d+)\s*:\s+({FLOAT_PATTERN})"
            rf"\s+at\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+"
            rf"({FLOAT_PATTERN})"
        )
        for match in re.finditer(pattern, stdout, re.IGNORECASE):
            ext_type = "min" if "min" in match[1].lower() else "max"
            extrema.append(
                {
                    "type": ext_type,
                    "index": int(match[2]),
                    "value": float(match[3]),
                    "position": (
                        float(match[4]),
                        float(match[5]),
                        float(match[6]),
                    ),
                }
            )
        return extrema


# =============================================================================
# Menu 15: Fuzzy atomic space analysis
# =============================================================================


class FuzzySpaceParser(OutputParser):
    """Parser for fuzzy atomic space analysis output.

    Handles integration results, atomic dipole moments, AOM,
    localization/delocalization indices, PDI, FLU, FLU-pi, MCI,
    ITA, atomic volumes, and IFDI.
    """

    @staticmethod
    def parse_atomic_properties(
        stdout: str,
    ) -> dict[int, dict[str, float]]:
        """Extract per-atom integrated properties.

        Returns
        -------
        dict[int, dict[str, float]]
            Mapping of atom index to property dict (e.g., 'population',
            'dipole_x', 'dipole_y', 'dipole_z', 'quadrupole')
        """
        atoms: dict[int, dict[str, float]] = {}

        # "Atom   1(C ):  population=  5.9678  dipole=  0.1234"
        pop_pattern = (
            rf"Atom\s+(\d+)\s*\([^)]+\)\s*:?\s+(?:population|integral)"
            rf"[=:\s]+({FLOAT_PATTERN})"
        )
        for match in re.finditer(pop_pattern, stdout, re.IGNORECASE):
            idx = int(match[1])
            atoms.setdefault(idx, {})["population"] = float(match[2])

        # Atomic dipole: "Atomic dipole of atom  1(C ):  X= 0.12  Y= 0.34  Z= 0.56"
        dip_pattern = (
            rf"Atomic dipole of atom\s+(\d+).*?"
            rf"X=\s*({FLOAT_PATTERN})\s+"
            rf"Y=\s*({FLOAT_PATTERN})\s+"
            rf"Z=\s*({FLOAT_PATTERN})"
        )
        for match in re.finditer(dip_pattern, stdout, re.IGNORECASE):
            idx = int(match[1])
            atoms.setdefault(idx, {}).update(
                {
                    "dipole_x": float(match[2]),
                    "dipole_y": float(match[3]),
                    "dipole_z": float(match[4]),
                }
            )

        # Atomic volume: "Atom  1  volume:  23.456"
        vol_pattern = rf"Atom\s+(\d+).*?volume[=:\s]+({FLOAT_PATTERN})"
        for match in re.finditer(vol_pattern, stdout, re.IGNORECASE):
            idx = int(match[1])
            atoms.setdefault(idx, {})["volume"] = float(match[2])

        return atoms

    @staticmethod
    def parse_delocalization_indices(
        stdout: str,
    ) -> dict[str, Any]:
        """Extract localization and delocalization indices.

        Returns
        -------
        dict[str, Any]
            'localization': dict[int, float] (LI per atom)
            'delocalization': dict[tuple[int,int], float] (DI per pair)
        """
        result: dict[str, Any] = {"localization": {}, "delocalization": {}}

        li_pattern = (
            rf"Localization index of atom\s+(\d+).*?:\s+({FLOAT_PATTERN})"
        )
        di_pattern = (
            rf"Delocalization index.*?atom\s+(\d+).*?atom\s+(\d+)"
            rf".*?:\s+({FLOAT_PATTERN})"
        )

        for match in re.finditer(li_pattern, stdout, re.IGNORECASE):
            result["localization"][int(match[1])] = float(match[2])

        for match in re.finditer(di_pattern, stdout, re.IGNORECASE):
            a1, a2 = int(match[1]), int(match[2])
            if a1 > a2:
                a1, a2 = a2, a1
            result["delocalization"][(a1, a2)] = float(match[3])

        return result

    @staticmethod
    def parse_aromaticity_index(stdout: str) -> dict[str, Any]:
        """Extract aromaticity indices (PDI, FLU, FLU-pi, MCI, ITA).

        Returns
        -------
        dict[str, Any]
            Dictionary with named index values found in the output
        """
        result: dict[str, Any] = {}

        index_patterns: dict[str, str] = {
            "PDI": rf"PDI[=:\s]+({FLOAT_PATTERN})",
            "FLU": rf"FLU[=:\s]+({FLOAT_PATTERN})",
            "FLU_pi": rf"FLU.*?pi[=:\s]+({FLOAT_PATTERN})",
            "MCI": rf"MCI[=:\s]+({FLOAT_PATTERN})",
            "Iring": rf"Iring[=:\s]+({FLOAT_PATTERN})",
            "ITA": rf"ITA[=:\s]+({FLOAT_PATTERN})",
            "PLR": rf"PLR[=:\s]+({FLOAT_PATTERN})",
        }

        for name, pat in index_patterns.items():
            if match := re.search(pat, stdout, re.IGNORECASE):
                result[name] = float(match[1])

        return result


# =============================================================================
# Menu 17: Basin analysis
# =============================================================================


class BasinParser(OutputParser):
    """Parser for basin analysis output.

    Handles AIM, ELF, ESP, LOL, and custom basin analyses.
    """

    @staticmethod
    def parse(stdout: str) -> list[dict[str, Any]]:
        """Extract basin integration results.

        Returns
        -------
        list[dict[str, Any]]
            List of dicts with 'basin', 'attractor', 'population',
            and optionally 'volume', 'charge'
        """
        basins: list[dict[str, Any]] = []

        # "Basin   1  attractor at atom  3(O )  population:  9.2345"
        pattern = (
            rf"Basin\s+(\d+).*?(?:attractor.*?atom\s+(\d+)\s*\(([^)]+)\))?"
            rf".*?population[=:\s]+({FLOAT_PATTERN})"
        )
        for match in re.finditer(pattern, stdout, re.IGNORECASE):
            basin: dict[str, Any] = {
                "basin": int(match[1]),
                "population": float(match[4]),
            }
            if match[2]:
                basin["attractor_atom"] = int(match[2])
                basin["attractor_element"] = match[3].strip()
            basins.append(basin)

        if not basins:
            # Simpler format: "  1   C    6.0123   23.456"
            simple = rf"^\s*(\d+)\s+([A-Za-z]+)\s+({FLOAT_PATTERN})"
            in_basin = False
            for line in stdout.split("\n"):
                if "basin" in line.lower() and (
                    "population" in line.lower() or "integral" in line.lower()
                ):
                    in_basin = True
                    continue
                if in_basin and (m := re.match(simple, line)):
                    basins.append(
                        {
                            "basin": int(m[1]),
                            "attractor_element": m[2],
                            "population": float(m[3]),
                        }
                    )

        return basins

    @staticmethod
    def parse_charges(stdout: str) -> dict[int, float]:
        """Extract AIM/Bader charges from basin analysis.

        Returns
        -------
        dict[int, float]
            Atom index to AIM charge
        """
        pattern = (
            rf"(?:AIM|Bader)\s+charge.*?atom\s+(\d+).*?:\s+({FLOAT_PATTERN})"
        )
        charges: dict[int, float] = {
            int(match[1]): float(match[2])
            for match in re.finditer(pattern, stdout, re.IGNORECASE)
        }
        return charges


# =============================================================================
# Menu 18: Electron excitation analysis
# =============================================================================


class ExcitationParser(OutputParser):
    """Parser for electron excitation analysis output.

    Handles hole-electron analysis, transition density matrix,
    charge transfer analysis, Delta_r index, NTOs, IFCT, Lambda
    index, CTS, and conditional density.
    """

    @staticmethod
    def parse_hole_electron(stdout: str) -> dict[str, Any]:
        """Extract hole-electron analysis descriptors.

        Returns
        -------
        dict[str, Any]
            Dictionary with 'D_index' (centroid distance), 'Sr'
            (overlap integral), 't_index', 'H_index', 'E_index',
            'hole_centroid', 'electron_centroid'
        """
        result: dict[str, Any] = {}

        patterns: dict[str, str] = {
            "D_index": rf"D\s+index[=:\s]+({FLOAT_PATTERN})",
            "Sr": rf"Sr[=:\s]+({FLOAT_PATTERN})",
            "t_index": rf"t\s+index[=:\s]+({FLOAT_PATTERN})",
            "H_index": rf"H\s+index[=:\s]+({FLOAT_PATTERN})",
            "E_index": rf"E\s+index[=:\s]+({FLOAT_PATTERN})",
            "HDI": rf"HDI[=:\s]+({FLOAT_PATTERN})",
            "EDI": rf"EDI[=:\s]+({FLOAT_PATTERN})",
        }

        for key, pat in patterns.items():
            if match := re.search(pat, stdout, re.IGNORECASE):
                result[key] = float(match[1])

        # Centroids
        centroid_pat = (
            rf"(hole|electron)\s+centroid.*?"
            rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})"
        )
        for match in re.finditer(centroid_pat, stdout, re.IGNORECASE):
            key = f"{match[1].lower()}_centroid"
            result[key] = (
                float(match[2]),
                float(match[3]),
                float(match[4]),
            )

        return result

    @staticmethod
    def parse_charge_transfer(stdout: str) -> dict[str, Any]:
        """Extract charge transfer analysis results.

        Returns
        -------
        dict[str, Any]
            Dictionary with 'CT_distance', 'CT_amount', and per-fragment
            transferred charge data
        """
        result: dict[str, Any] = {}

        ct_dist = rf"CT\s+distance[=:\s]+({FLOAT_PATTERN})"
        ct_amt = (
            rf"(?:transferred|CT)\s+(?:charge|amount)[=:\s]+({FLOAT_PATTERN})"
        )

        if match := re.search(ct_dist, stdout, re.IGNORECASE):
            result["CT_distance"] = float(match[1])
        if match := re.search(ct_amt, stdout, re.IGNORECASE):
            result["CT_amount"] = float(match[1])

        # Per-fragment: "Fragment 1:  hole= 0.85  electron= 0.15"
        frag_pattern = (
            rf"Fragment\s+(\d+).*?hole[=:\s]+({FLOAT_PATTERN})"
            rf".*?electron[=:\s]+({FLOAT_PATTERN})"
        )
        fragments: list[dict[str, Any]] = []
        fragments.extend(
            {
                "fragment": int(match[1]),
                "hole": float(match[2]),
                "electron": float(match[3]),
            }
            for match in re.finditer(frag_pattern, stdout, re.IGNORECASE)
        )
        if fragments:
            result["fragments"] = fragments

        return result

    @staticmethod
    def parse_delta_r(stdout: str) -> list[dict[str, float]]:
        """Extract Delta_r index for each excited state.

        Returns
        -------
        list[dict[str, float]]
            List of dicts with 'state' and 'delta_r'
        """
        pattern = (
            rf"(?:State|Excited state)\s+(\d+).*?"
            rf"Delta_?r[=:\s]+({FLOAT_PATTERN})"
        )
        results: list[dict[str, float]] = [
            {
                "state": float(match[1]),
                "delta_r": float(match[2]),
            }
            for match in re.finditer(pattern, stdout, re.IGNORECASE)
        ]
        return results

    @staticmethod
    def parse_lambda_index(stdout: str) -> list[dict[str, float]]:
        """Extract Lambda diagnostic for each excited state.

        Returns
        -------
        list[dict[str, float]]
            List of dicts with 'state' and 'lambda_index'
        """
        pattern = (
            rf"(?:State|Excited state)\s+(\d+).*?"
            rf"Lambda[=:\s]+({FLOAT_PATTERN})"
        )
        results: list[dict[str, float]] = [
            {
                "state": float(match[1]),
                "lambda_index": float(match[2]),
            }
            for match in re.finditer(pattern, stdout, re.IGNORECASE)
        ]
        return results


# =============================================================================
# Menu 20: Weak interaction analysis
# =============================================================================


class WeakInteractionParser(OutputParser):
    """Parser for weak interaction analysis output.

    Handles NCI, IRI, DORI, IGM, IGMH, and related analyses. These
    methods primarily produce cube files; this parser extracts any
    summary statistics printed to stdout.
    """

    @staticmethod
    def parse(stdout: str) -> dict[str, Any]:
        """Extract summary statistics from weak interaction analysis.

        Returns
        -------
        dict[str, Any]
            Dictionary with descriptors like 'delta_g_inter',
            'delta_g_intra', 'intrinsic_bond_strength', and any
            reported isosurface integral values.
        """
        result: dict[str, Any] = {}

        patterns: dict[str, str] = {
            "delta_g_inter": (rf"delta_?g_?inter.*?[=:\s]+({FLOAT_PATTERN})"),
            "delta_g_intra": (rf"delta_?g_?intra.*?[=:\s]+({FLOAT_PATTERN})"),
            "isosurface_integral": (
                rf"(?:Integral|integral).*?isosurface.*?[=:\s]+({FLOAT_PATTERN})"
            ),
        }

        for key, pat in patterns.items():
            if match := re.search(pat, stdout, re.IGNORECASE):
                result[key] = float(match[1])

        # Cube files produced
        cubes: list[str] = []
        cubes.extend(
            match[1]
            for match in re.finditer(
                r"(\S+\.cube)\s+has been generated", stdout
            )
        )
        if cubes:
            result["cube_files"] = cubes

        return result


# =============================================================================
# Menu 21: Energy Decomposition Analysis (EDA)
# =============================================================================


class EDAParser(OutputParser):
    """Parser for energy decomposition analysis output.

    Handles EDA_FF, EDA_SBL, SOBEDA, and dispersion contributions.
    """

    @staticmethod
    def parse(stdout: str) -> dict[str, float]:
        """Extract EDA energy components.

        Returns
        -------
        dict[str, float]
            Dictionary with energy components in kcal/mol:
            'electrostatic', 'exchange', 'repulsion', 'polarization',
            'dispersion', 'total_interaction'
        """
        result: dict[str, float] = {}

        patterns: dict[str, str] = {
            "electrostatic": (rf"[Ee]lectrostatic.*?[=:\s]+({FLOAT_PATTERN})"),
            "exchange": rf"[Ee]xchange.*?[=:\s]+({FLOAT_PATTERN})",
            "repulsion": (
                rf"(?:[Rr]epulsion|Pauli).*?[=:\s]+({FLOAT_PATTERN})"
            ),
            "polarization": (
                rf"[Pp]olari[sz]ation.*?[=:\s]+({FLOAT_PATTERN})"
            ),
            "dispersion": rf"[Dd]ispersion.*?[=:\s]+({FLOAT_PATTERN})",
            "orbital_interaction": (
                rf"[Oo]rbital\s+interaction.*?[=:\s]+({FLOAT_PATTERN})"
            ),
            "total_interaction": (
                rf"[Tt]otal\s+interaction.*?[=:\s]+({FLOAT_PATTERN})"
            ),
        }

        for key, pat in patterns.items():
            if match := re.search(pat, stdout):
                result[key] = float(match[1])

        return result

    @staticmethod
    def parse_dispersion_contributions(
        stdout: str,
    ) -> dict[int, float]:
        """Extract per-atom dispersion energy contributions.

        Returns
        -------
        dict[int, float]
            Atom index to dispersion contribution
        """
        pattern = (
            rf"Atom\s+(\d+).*?(?:dispersion|D[34]).*?[=:\s]+({FLOAT_PATTERN})"
        )
        contributions: dict[int, float] = {
            int(match[1]): float(match[2])
            for match in re.finditer(pattern, stdout, re.IGNORECASE)
        }
        return contributions


# =============================================================================
# Menu 22: Conceptual DFT (CDFT)
# =============================================================================


class CDFTParser(OutputParser):
    """Parser for conceptual DFT output.

    Handles global reactivity indices, Fukui functions, dual descriptor,
    condensed Fukui, local hardness, and local ionisation energy.
    """

    @staticmethod
    def parse_global_indices(stdout: str) -> dict[str, float]:
        """Extract global CDFT reactivity descriptors.

        Returns
        -------
        dict[str, float]
            Dictionary with 'chemical_potential', 'hardness',
            'softness', 'electrophilicity', 'nucleophilicity',
            'IP', 'EA'
        """
        result: dict[str, float] = {}

        patterns: dict[str, str] = {
            "chemical_potential": (
                rf"[Cc]hemical potential.*?[=:\s]+({FLOAT_PATTERN})"
            ),
            "hardness": rf"[Hh]ardness.*?[=:\s]+({FLOAT_PATTERN})",
            "softness": rf"[Ss]oftness.*?[=:\s]+({FLOAT_PATTERN})",
            "electrophilicity": (
                rf"[Ee]lectrophilicity.*?[=:\s]+({FLOAT_PATTERN})"
            ),
            "nucleophilicity": (
                rf"[Nn]ucleophilicity.*?[=:\s]+({FLOAT_PATTERN})"
            ),
            "IP": rf"(?:IP|Ionization potential).*?[=:\s]+({FLOAT_PATTERN})",
            "EA": rf"(?:EA|Electron affinity).*?[=:\s]+({FLOAT_PATTERN})",
        }

        for key, pat in patterns.items():
            if match := re.search(pat, stdout, re.IGNORECASE):
                result[key] = float(match[1])

        return result

    @staticmethod
    def parse_condensed_fukui(
        stdout: str,
    ) -> dict[int, dict[str, float]]:
        """Extract condensed Fukui function values per atom.

        Returns
        -------
        dict[int, dict[str, float]]
            Mapping of atom index to dict with 'f_plus', 'f_minus',
            'f_zero' values
        """
        # "Atom   1(C ):  f+= 0.1234  f-= 0.0567  f0= 0.0900"
        pattern = (
            rf"Atom\s+(\d+)\s*\([^)]+\)\s*:?\s+"
            rf"f\+[=:\s]+({FLOAT_PATTERN})\s+"
            rf"f-[=:\s]+({FLOAT_PATTERN})\s+"
            rf"f0[=:\s]+({FLOAT_PATTERN})"
        )
        fukui: dict[int, dict[str, float]] = {
            int(match[1]): {
                "f_plus": float(match[2]),
                "f_minus": float(match[3]),
                "f_zero": float(match[4]),
            }
            for match in re.finditer(pattern, stdout)
        }
        if not fukui:
            # Try separate patterns
            for label, key in [
                (r"f\+", "f_plus"),
                (r"f-", "f_minus"),
                (r"f0", "f_zero"),
            ]:
                pat = (
                    rf"Atom\s+(\d+)\s*\([^)]+\)\s*:?\s+"
                    rf"{label}[=:\s]+({FLOAT_PATTERN})"
                )
                for m in re.finditer(pat, stdout):
                    fukui.setdefault(int(m[1]), {})[key] = float(m[2])

        return fukui

    @staticmethod
    def parse_dual_descriptor(
        stdout: str,
    ) -> dict[int, float]:
        """Extract condensed dual descriptor per atom.

        Returns
        -------
        dict[int, float]
            Atom index to dual descriptor value (positive = nucleophilic)
        """
        pattern = (
            rf"Atom\s+(\d+)\s*\([^)]+\)\s*:?\s+"
            rf"(?:dual|Delta_?f|f\+\s*-\s*f-)[=:\s]+({FLOAT_PATTERN})"
        )
        dd: dict[int, float] = {
            int(match[1]): float(match[2])
            for match in re.finditer(pattern, stdout, re.IGNORECASE)
        }
        return dd


# =============================================================================
# Menu 24: Polarizability
# =============================================================================


class PolarizabilityParser(OutputParser):
    """Parser for polarizability output.

    Handles parsed polarizability/hyperpolarizability tensors and SOS
    results.
    """

    @staticmethod
    def parse(stdout: str) -> dict[str, Any]:
        """Extract polarizability tensor and related quantities.

        Returns
        -------
        dict[str, Any]
            Dictionary with 'alpha_iso', 'alpha_aniso', tensor
            components ('alpha_xx', etc.), and optionally
            'beta_total', 'gamma_total'
        """
        result: dict[str, Any] = {}

        # Isotropic polarizability
        iso_pat = rf"[Ii]sotropic.*?(?:alpha|polarizability).*?[=:\s]+({FLOAT_PATTERN})"
        if match := re.search(iso_pat, stdout):
            result["alpha_iso"] = float(match[1])

        # Anisotropy
        aniso_pat = rf"[Aa]nisotropy.*?[=:\s]+({FLOAT_PATTERN})"
        if match := re.search(aniso_pat, stdout):
            result["alpha_aniso"] = float(match[1])

        # Tensor components
        for comp in ["xx", "xy", "xz", "yy", "yz", "zz"]:
            pat = rf"alpha[_\s]*{comp}[=:\s]+({FLOAT_PATTERN})"
            if match := re.search(pat, stdout, re.IGNORECASE):
                result[f"alpha_{comp}"] = float(match[1])

        # Hyperpolarizability
        beta_pat = rf"[Bb]eta.*?total[=:\s]+({FLOAT_PATTERN})"
        if match := re.search(beta_pat, stdout):
            result["beta_total"] = float(match[1])

        gamma_pat = rf"[Gg]amma.*?(?:total|average)[=:\s]+({FLOAT_PATTERN})"
        if match := re.search(gamma_pat, stdout):
            result["gamma_total"] = float(match[1])

        return result


# =============================================================================
# Menu 25: Aromaticity
# =============================================================================


class AromaticityParser(OutputParser):
    """Parser for aromaticity analysis output.

    Handles NICS, ICSS, HOMA, HOMAC, HOMER, Bird, Stanger, and AICD.
    """

    @staticmethod
    def parse(stdout: str) -> dict[str, Any]:
        """Extract aromaticity indices.

        Returns
        -------
        dict[str, Any]
            Dictionary with named aromaticity index values
        """
        result: dict[str, Any] = {}

        patterns: dict[str, str] = {
            "NICS": rf"NICS\s*(?:\(0\))?[=:\s]+({FLOAT_PATTERN})",
            "NICS_1": rf"NICS\s*\(1\)[=:\s]+({FLOAT_PATTERN})",
            "NICS_ZZ": rf"NICS_?ZZ[=:\s]+({FLOAT_PATTERN})",
            "HOMA": rf"HOMA[=:\s]+({FLOAT_PATTERN})",
            "HOMAC": rf"HOMAC[=:\s]+({FLOAT_PATTERN})",
            "HOMER": rf"HOMER[=:\s]+({FLOAT_PATTERN})",
            "Bird": rf"Bird.*?index[=:\s]+({FLOAT_PATTERN})",
            "EN_GEO": rf"EN[_\s]*GEO[=:\s]+({FLOAT_PATTERN})",
            "EN_BLA": rf"EN[_\s]*BLA[=:\s]+({FLOAT_PATTERN})",
        }

        for key, pat in patterns.items():
            if match := re.search(pat, stdout, re.IGNORECASE):
                result[key] = float(match[1])

        return result

    @staticmethod
    def parse_nics_scan(
        stdout: str,
    ) -> dict[str, list[float]]:
        """Extract NICS scan profile data.

        Returns
        -------
        dict[str, list[float]]
            Dictionary with 'distances' and 'nics_values' lists
        """
        data: dict[str, list[float]] = {
            "distances": [],
            "nics_values": [],
        }
        pattern = rf"^\s*({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s*$"
        in_scan = False
        for line in stdout.split("\n"):
            if "nics" in line.lower() and "scan" in line.lower():
                in_scan = True
                continue
            if in_scan and (match := re.match(pattern, line)):
                data["distances"].append(float(match[1]))
                data["nics_values"].append(float(match[2]))
        return data


# =============================================================================
# Menu 6: Wavefunction info
# =============================================================================


class WavefunctionParser(OutputParser):
    """Parser for wavefunction check/modify output.

    Handles orbital info, GTF info, basis info, density matrix,
    and overlap matrix output.
    """

    @staticmethod
    def parse_orbital_info(stdout: str) -> list[dict[str, Any]]:
        """Extract orbital information.

        Returns
        -------
        list[dict[str, Any]]
            List of dicts with 'index', 'energy', 'occupation', 'spin'
        """
        # "   5   Alpha   Occ= 2.000000   E=  -0.72340 a.u.  -19.684 eV"
        pattern = (
            rf"(\d+)\s+(Alpha|Beta)\s+Occ=\s*({FLOAT_PATTERN})\s+"
            rf"E=\s*({FLOAT_PATTERN})\s*a\.u\.\s+({FLOAT_PATTERN})\s*eV"
        )
        orbitals: list[dict[str, Any]] = [
            {
                "index": int(match[1]),
                "spin": match[2].lower(),
                "occupation": float(match[3]),
                "energy_au": float(match[4]),
                "energy_eV": float(match[5]),
            }
            for match in re.finditer(pattern, stdout)
        ]
        return orbitals


# =============================================================================
# Menu 5 / 13: Cube / grid operations
# =============================================================================


class CubeParser(OutputParser):
    """Parser for cube generation and grid processing output.

    Extracts file paths and grid metadata from cube generation, and
    statistics from grid processing operations.
    """

    @staticmethod
    def parse(stdout: str) -> dict[str, Any]:
        """Extract cube file generation info and grid statistics.

        Returns
        -------
        dict[str, Any]
            Dictionary with 'cube_files' (list of generated files),
            'grid_points' (tuple of nx, ny, nz), 'min', 'max', 'mean',
            'integral'
        """
        result: dict[str, Any] = {}

        # Cube files generated
        cubes: list[str] = []
        cubes.extend(
            match[1]
            for match in re.finditer(
                r"(\S+\.cube)\s+has been generated", stdout
            )
        )
        if cubes:
            result["cube_files"] = cubes

        # Grid dimensions
        grid_pat = r"Grid dimensions:\s*(\d+)\s*x\s*(\d+)\s*x\s*(\d+)"
        if match := re.search(grid_pat, stdout):
            result["grid_points"] = (
                int(match[1]),
                int(match[2]),
                int(match[3]),
            )

        # Statistics
        stat_patterns: dict[str, str] = {
            "min": rf"[Mm]inimum.*?[=:\s]+({FLOAT_PATTERN})",
            "max": rf"[Mm]aximum.*?[=:\s]+({FLOAT_PATTERN})",
            "mean": rf"[Mm]ean.*?[=:\s]+({FLOAT_PATTERN})",
            "integral": rf"[Ii]ntegral.*?[=:\s]+({FLOAT_PATTERN})",
            "std_dev": rf"[Ss]td.*?dev.*?[=:\s]+({FLOAT_PATTERN})",
        }

        for key, pat in stat_patterns.items():
            if match := re.search(pat, stdout):
                result[key] = float(match[1])

        return result


# =============================================================================
# Menu 100/200/300: Utilities
# =============================================================================
class UtilityParser(OutputParser):
    """Parser for utility function outputs.

    Handles geometry properties, electric multipole moments, Hellmann-
    Feynman forces, coordination numbers, BLA/BOA, and various other
    utility results.
    """

    @staticmethod
    def parse_geometry(stdout: str) -> dict[str, Any]:
        """Extract geometry properties (bond lengths, angles, dihedrals).

        Returns
        -------
        dict[str, Any]
            Dictionary with 'bond_lengths', 'angles', 'dihedrals' lists
        """
        result: dict[str, Any] = {
            "bond_lengths": [],
            "angles": [],
            "dihedrals": [],
        }

        # "Bond length between atom  1(C ) and  2(N ):  1.3456 Angstrom"
        bl_pattern = (
            rf"Bond length.*?atom\s+(\d+).*?(?:and|atom)\s+(\d+).*?:\s+"
            rf"({FLOAT_PATTERN})"
        )
        for match in re.finditer(bl_pattern, stdout, re.IGNORECASE):
            result["bond_lengths"].append(
                {
                    "atom1": int(match[1]),
                    "atom2": int(match[2]),
                    "length": float(match[3]),
                }
            )

        # "Angle  1-2-3 :  120.345 degree"
        ang_pattern = (
            rf"[Aa]ngle\s+(\d+)\s*-\s*(\d+)\s*-\s*(\d+)\s*:\s+"
            rf"({FLOAT_PATTERN})"
        )
        for match in re.finditer(ang_pattern, stdout):
            result["angles"].append(
                {
                    "atoms": (int(match[1]), int(match[2]), int(match[3])),
                    "angle": float(match[4]),
                }
            )

        # "Dihedral  1-2-3-4 :  -45.678 degree"
        dih_pattern = (
            rf"[Dd]ihedral\s+(\d+)\s*-\s*(\d+)\s*-\s*(\d+)\s*-\s*(\d+)"
            rf"\s*:\s+({FLOAT_PATTERN})"
        )
        for match in re.finditer(dih_pattern, stdout):
            result["dihedrals"].append(
                {
                    "atoms": (
                        int(match[1]),
                        int(match[2]),
                        int(match[3]),
                        int(match[4]),
                    ),
                    "dihedral": float(match[5]),
                }
            )

        return result

    @staticmethod
    def parse_multipole_moments(stdout: str) -> dict[str, Any]:
        """Extract electric multipole moments.

        Returns
        -------
        dict[str, Any]
            Dictionary with 'dipole', 'quadrupole', 'octapole',
            'hexadecapole' components
        """
        result: dict[str, Any] = {}

        dip_pat = (
            rf"[Dd]ipole.*?X[=:\s]+({FLOAT_PATTERN})\s+"
            rf"Y[=:\s]+({FLOAT_PATTERN})\s+Z[=:\s]+({FLOAT_PATTERN})"
        )
        if match := re.search(dip_pat, stdout):
            result["dipole"] = {
                "x": float(match[1]),
                "y": float(match[2]),
                "z": float(match[3]),
            }

        quad_components = ["XX", "XY", "XZ", "YY", "YZ", "ZZ"]
        quad: dict[str, float] = {}
        for comp in quad_components:
            pat = rf"[Qq]uadrupole.*?{comp}[=:\s]+({FLOAT_PATTERN})"
            if m := re.search(pat, stdout):
                quad[comp] = float(m[1])
        if quad:
            result["quadrupole"] = quad

        return result

    @staticmethod
    def parse_coordination_numbers(
        stdout: str,
    ) -> dict[int, float]:
        """Extract atomic coordination numbers.

        Returns
        -------
        dict[int, float]
            Atom index to coordination number
        """
        pattern = rf"Atom\s+(\d+).*?coordination.*?[=:\s]+({FLOAT_PATTERN})"
        coords: dict[int, float] = {
            int(match[1]): float(match[2])
            for match in re.finditer(pattern, stdout, re.IGNORECASE)
        }
        return coords

    @staticmethod
    def parse_bla_boa(stdout: str) -> dict[str, float]:
        """Extract BLA and BOA values.

        Returns
        -------
        dict[str, float]
            Dictionary with 'BLA' and 'BOA' values
        """
        result: dict[str, float] = {}
        for key in ["BLA", "BOA"]:
            pat = rf"{key}[=:\s]+({FLOAT_PATTERN})"
            if match := re.search(pat, stdout, re.IGNORECASE):
                result[key] = float(match[1])
        return result

    @staticmethod
    def parse_generated_files(stdout: str) -> list[str]:
        """Extract paths of files generated during the analysis.

        Returns
        -------
        list[str]
            List of generated file paths
        """
        files: list[str] = []
        for pattern in [
            r"(\S+\.cube)\s+has been generated",
            r"(\S+\.wfn)\s+has been generated",
            r"(\S+\.molden)\s+has been generated",
            r"(\S+\.fch)\s+has been generated",
            r"(\S+\.xyz)\s+has been generated",
            r"(\S+\.pdb)\s+has been generated",
            r"(\S+\.txt)\s+has been generated",
            r"(\S+\.gjf)\s+has been generated",
        ]:
            files.extend(match[1] for match in re.finditer(pattern, stdout))
        return files