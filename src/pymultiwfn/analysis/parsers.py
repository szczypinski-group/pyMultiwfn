"""Output parsers for Multiwfn results.

ClaudeCode-extracted regex patterns from Multiwfn 3.8 manual (6 Jan 2026).

"""

import re
from typing import Any, Literal

from pymultiwfn.analysis.result import (
    BLA_BOA,
    Aromaticity,
    AromaticityIndex,
    Basin,
    BondAngle,
    BondLength,
    BondOrder,
    BondOrderDecomposition,
    BondPath,
    Charge,
    ChargeTransfer,
    ChargeTransferFragment,
    Color,
    CondensedFukui,
    CoordinationNumber,
    CriticalPoint,
    Cube,
    DelocalizationIndex,
    DeltaR,
    DensityOfStates,
    DihedralAngle,
    Dipole,
    DipoleMoment,
    DispersionContribution,
    DualDescriptor,
    EnergyDecompositionAnalysis,
    FuzzyAtomicProperty,
    HoleElectron,
    LambdaIndex,
    MultiCenterBondOrder,
    MultiwfnResult,
    NICSScan,
    Orbital,
    OrbitalComponent,
    OrbitalContribution,
    OrbitalEnergy,
    OxidationState,
    Polarizability,
    PolarizabilityTensor,
    QuadrupoleMoment,
    Reactivity,
    Spectrum,
    SurfaceAnalysis,
    SurfaceExtremum,
    Transition,
    Valence,
    WeakInteraction,
)
from pymultiwfn.enums.menu import Menu

# Shared float pattern for all numeric parsing
FLOAT_PATTERN = r"[-+]?(?:\d+\.\d*|\d*\.\d+|\d+)(?:[Ee][-+]?\d+)?"


class OutputParser:
    """Base class for output parsers."""

    pass


# =============================================================================
# Menu 7: Population analysis & atomic charges
# =============================================================================


class ChargeParser(
    OutputParser,
    MultiwfnResult,
):
    """Parser for atomic charges.

    Handles output from all Menu 7 charge methods including Hirshfeld,
    VDD, Mulliken, Lowdin, SCPA, Stout-Politzer, Bickelhaupt, Becke,
    ADCH, CHELPG, MK, AIM, Hirshfeld-I, CM5, EEM, RESP, Gasteiger,
    MBIS, and DDEC.
    """

    def _parse(self, stdout: str) -> None:
        self.result.extend(ChargeParser.parse_charges(stdout))
        if (dipole := ChargeParser.parse_dipole(stdout)) is not None:
            self.result.extend([dipole])

    @staticmethod
    def parse_charges(stdout: str) -> list[Charge]:
        """Extract atomic charges from Multiwfn output.

        Parameters
        ----------
        stdout
            Multiwfn standard output

        Returns
        -------
        ChargeResult
            Result object containing list of charge objects

        """
        charges: list[Charge] = []

        # Pattern 1: "Hirshfeld charge of atom     1(C ) is  0.03208687"
        pattern1 = (
            rf"charge of atom\s+(\d+)\s*\([A-Za-z]+"
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
            # TODO(fs): lots of if statements here, worth refactoring i think?
            if "final atomic charges" in line_lower:
                in_final_section = True
                in_charge_section = True
                charges.clear()  # Prefer final charges
                continue
            elif "charge" in line_lower or "population" in line_lower:
                in_charge_section = True
                continue
            elif "summary of" in line_lower and "charge" in line_lower:
                in_charge_section = True
                charges.clear()
                continue

            # Try pattern 1 first (most specific - includes method name)
            if match := re.search(pattern1, line, re.IGNORECASE):
                charges.append(
                    Charge(atom_id=int(match[1]), charge=float(match[2]))
                )
                continue

            # Try pattern 5 (population format for Mulliken/Lowdin)
            if match := re.search(pattern5, line):
                if in_charge_section:
                    charges.append(
                        Charge(atom_id=int(match[1]), charge=float(match[2]))
                    )
                continue

            # Try pattern 6 (generic charge format)
            if match := re.search(pattern6, line):
                if in_charge_section:
                    charges.append(
                        Charge(atom_id=int(match[1]), charge=float(match[2]))
                    )
                continue

            # Try pattern 2 (Final atomic charges format)
            if match := re.search(pattern2, line):
                if in_charge_section or in_final_section:
                    charges.append(
                        Charge(atom_id=int(match[1]), charge=float(match[2]))
                    )
                continue

            # Try pattern 3 (table format)
            if in_charge_section and (match := re.match(pattern3, line)):
                charges.append(
                    Charge(atom_id=int(match[1]), charge=float(match[2]))
                )
                continue

            # Try pattern 4 (column format)
            if in_charge_section and (match := re.match(pattern4, line)):
                charges.append(
                    Charge(atom_id=int(match[1]), charge=float(match[2]))
                )
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
    def parse_dipole(stdout: str) -> Dipole | None:
        """Extract molecular dipole moment from charge output.

        Returns
        -------
        Dictionary with 'x', 'y', 'z', 'total' dipole components
        in Debye, or None if not found.

        """
        pattern = (
            rf"Dipole moment.*?X=\s*({FLOAT_PATTERN})\s+"
            rf"Y=\s*({FLOAT_PATTERN})\s+Z=\s*({FLOAT_PATTERN})\s+"
            rf"Tot=\s*({FLOAT_PATTERN})"
        )
        if match := re.search(pattern, stdout, re.IGNORECASE):
            return Dipole(
                x=float(match[1]),
                y=float(match[2]),
                z=float(match[3]),
                total=float(match[4]),
            )
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
    def parse(stdout: str) -> list[OrbitalComponent]:
        """Extract orbital composition data.

        Returns
        -------
        List of dicts with keys: 'orbital', 'energy', 'occupation',
        'contributions' (dict mapping atom/fragment labels to percentages)

        """
        orbitals: list[OrbitalComponent] = []

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

        current_orb: OrbitalComponent | None = None
        for line in stdout.split("\n"):
            if match := re.search(orb_pattern, line):
                if current_orb is not None:
                    orbitals.append(current_orb)
                current_orb = OrbitalComponent(
                    orbital_id=int(match[1]),
                    occupation=float(match[2]),
                    energy=float(match[3]),
                    contributions=[],
                )
            elif current_orb is not None:
                if match := re.search(contrib_pattern, line):
                    label = f"{match[1]}{match[2]}"
                    current_orb.contributions.append(
                        OrbitalContribution(
                            label=label, percentage=float(match[3])
                        )
                    )
                elif match := re.search(contrib_pattern2, line):
                    current_orb.contributions.append(
                        OrbitalContribution(
                            label=f"atom_{match[1]}",
                            percentage=float(match[2]),
                        )
                    )

        if current_orb is not None:
            orbitals.append(current_orb)
        return orbitals

    @staticmethod
    def parse_oxidation_states(stdout: str) -> list[OxidationState]:
        """Extract LOBA oxidation states.

        Returns
        -------
        List of oxidation states for each atom.

        """
        pattern = (
            rf"Atom\s+(\d+)\s*\([A-Za-z]+\s*\).*?"
            rf"oxidation state.*?({FLOAT_PATTERN})"
        )
        states: list[OxidationState] = [
            OxidationState(
                atom_id=int(match[1]),
                oxidation_state=int(round(float(match[2]))),
            )
            for match in re.finditer(pattern, stdout, re.IGNORECASE)
        ]
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
    def parse(
        stdout: str,
    ) -> list[BondOrder]:
        """Extract bond orders from Multiwfn output.

        Parameters
        ----------
        stdout
            Multiwfn standard output

        Returns
        -------
        List of bond orders

        """
        bond_orders: list[BondOrder] = []

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
                bond_orders.append(
                    BondOrder(atom1_id=atom1, atom2_id=atom2, bond_order=bo)
                )

        return bond_orders

    @staticmethod
    def parse_valence(stdout: str) -> list[Valence]:
        """Extract total valence and free valence for each atom.

        Returns
        -------
        List of valence information for each atom

        """
        valences: list[Valence] = []

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
            valences.append(
                Valence(
                    atom_id=idx, type="total_valence", valence=float(match[2])
                )
            )
        for match in re.finditer(free_pattern, stdout):
            idx = int(match[1])
            valences.append(
                Valence(
                    atom_id=idx, type="free_valence", valence=float(match[2])
                )
            )

        return valences

    @staticmethod
    def parse_multicenter(stdout: str) -> list[MultiCenterBondOrder]:
        """Extract multicenter bond order results.

        Returns
        -------
        List of multicenter bond order results

        """
        results: list[MultiCenterBondOrder] = []
        # "Multi-center bond order of atoms  1  2  3 :  0.12345"
        pattern = (
            rf"Multi-center bond order of atoms\s+([\d\s]+):\s+"
            rf"({FLOAT_PATTERN})"
        )
        for match in re.finditer(pattern, stdout, re.IGNORECASE):
            atoms = [int(x) for x in match[1].split()]
            results.append(
                MultiCenterBondOrder(
                    atom_ids=atoms, bond_order=float(match[2])
                )
            )
        return results

    @staticmethod
    def parse_decomposition(
        stdout: str,
    ) -> list[BondOrderDecomposition]:
        """Extract per-orbital bond order decomposition.

        Returns
        -------
        List of bond order decomposition results

        """
        decomp: list[BondOrderDecomposition] = []
        # "Orbital   5:   0.23456"
        pattern = rf"Orbital\s+(\d+)\s*:\s+({FLOAT_PATTERN})"
        decomp.extend(
            BondOrderDecomposition(
                orbital_id=int(match[1]),
                contribution=float(match[2]),
            )
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

    CP_TYPE_NAMES: dict[
        str, Literal["nuclear", "bond", "ring", "cage", "unknown"]
    ] = {
        "(3,-3)": "nuclear",
        "(3,-1)": "bond",
        "(3,+1)": "ring",
        "(3,+3)": "cage",
    }

    @staticmethod
    def parse(
        stdout: str,
    ) -> list[CriticalPoint]:
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
        cps: list[CriticalPoint] = []

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

                # Scan following lines for properties
                for j in range(i, min(i + 10, len(lines))):
                    sub = lines[j]
                    if pos_match := re.search(pos_pattern, sub):
                        position = (
                            float(pos_match[1]),
                            float(pos_match[2]),
                            float(pos_match[3]),
                        )
                    if rho_match := re.search(rho_pattern, sub):
                        rho = float(rho_match[1])
                    if lap_match := re.search(lap_pattern, sub):
                        lap = float(lap_match[1])
                    if ell_match := re.search(ell_pattern, sub):
                        ell = float(ell_match[1])

                    cp = CriticalPoint(
                        index=cp_index,
                        x=position[0],
                        y=position[1],
                        z=position[2],
                        rho=rho if rho is not None else None,
                        laplacian=lap if lap is not None else None,
                        ellipticity=ell if ell is not None else None,
                        type=CriticalPointParser.CP_TYPE_NAMES.get(
                            cp_type, "unknown"
                        ),
                    )

                cps.append(cp)

        return cps

    @staticmethod
    def parse_bond_paths(
        stdout: str,
    ) -> list[BondPath]:
        """Extract bond path information.

        Returns
        -------
        list[BondPath]
            List of BondPath instances
        """
        paths: list[BondPath] = []
        # Example:
        # Bond path between atom  1(C ) and atom  2(N ), BCP  3, length  2.456
        pattern = (
            rf"Bond path between atom\s+(\d+).*?and atom\s+(\d+).*?"
            rf"BCP\s+(\d+).*?length\s+({FLOAT_PATTERN})"
        )
        for match in re.finditer(pattern, stdout, re.IGNORECASE):
            if match is not None:
                paths.append(
                    BondPath(
                        atom1_id=int(match[1]),
                        atom2_id=int(match[2]),
                        bcp_id=int(match[3]),
                        path_length=float(match[4]),
                    )
                )

        return paths


# =============================================================================
# Menu 10: Density of states
# =============================================================================


class DOSParser(OutputParser):
    """Parser for density of states output.

    Handles TDOS, PDOS, OPDOS, LDOS, COHP, and photoelectron spectrum.
    """

    @staticmethod
    def parse(stdout: str) -> DensityOfStates:
        """Extract DOS curve data.

        Returns
        -------
        Density of states results object.

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
                data["energies_eV"].append(energy)
                if vals:
                    data["dos"].append(vals[0])
                for k, v in enumerate(vals[1:], start=1):
                    key = f"pdos_{k}"
                    data.setdefault(key, []).append(v)
            elif match := re.match(pattern2, line):
                data["energies_eV"].append(float(match[1]))
                data["dos"].append(float(match[2]))

        pdos = {k: v for k, v in data.items() if k.startswith("pdos_")}
        return DensityOfStates(
            energies_eV=data["energies_eV"],
            dos=data["dos"],
            projected_dos=pdos if pdos else None,
        )

    @staticmethod
    def parse_orbital_energies(stdout: str) -> list[OrbitalEnergy]:
        """Extract orbital energies used in DOS.

        Returns
        -------
        list[OrbitalEnergy]
            List of orbital energy objects
        """
        pattern = rf"(\d+)\s+({FLOAT_PATTERN})\s+eV\s+Occ=\s*({FLOAT_PATTERN})"
        orbitals: list[OrbitalEnergy] = [
            OrbitalEnergy(
                index=int(match[1]),
                energy_eV=float(match[2]),
                occupation=float(match[3]),
            )
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
    def parse(
        stdout: str,
    ) -> Spectrum:
        """Extract spectrum data (frequencies/wavelengths, intensities).

        Parameters
        ----------
        stdout : str
            Multiwfn standard output

        Returns
        -------
        Spectrum
            Spectrum result object
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
            return Spectrum(
                frequencies=spectrum["frequencies"],
                intensities=spectrum["intensities"],
            )

        # Try UV-Vis/ECD
        uv_data: dict[str, list[float]] = {
            "wavelengths": [],
            "intensities": [],
        }
        for match in re.finditer(pattern_uv, stdout):
            uv_data["wavelengths"].append(float(match[1]))
            uv_data["intensities"].append(float(match[2]))
        if uv_data["wavelengths"]:
            return Spectrum(
                wavelengths=uv_data["wavelengths"],
                intensities=uv_data["intensities"],
            )

        # Try NMR
        nmr_data: dict[str, list[float]] = {
            "atom_indices": [],
            "chemical_shifts": [],
        }
        for match in re.finditer(pattern_nmr, stdout):
            nmr_data["atom_indices"].append(float(match[1]))
            nmr_data["chemical_shifts"].append(float(match[2]))
        if nmr_data["atom_indices"]:
            return Spectrum(
                atom_indices=[int(x) for x in nmr_data["atom_indices"]],
                chemical_shifts=nmr_data["chemical_shifts"],
            )

        # Fallback: generic two-column
        for line in stdout.split("\n"):
            if match2 := re.match(pattern2, line):
                spectrum["frequencies"].append(float(match2[1]))
                spectrum["intensities"].append(float(match2[2]))

        return Spectrum(
            frequencies=spectrum["frequencies"],
            intensities=spectrum["intensities"],
        )

    @staticmethod
    def parse_transitions(stdout: str) -> list[Transition]:
        """Extract discrete transition data (excitation energies, strengths).

        Returns
        -------
        list[Transition]
            List of transition objects
        """
        transitions: list[Transition] = []
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
                rot_strength = None
                # Check for rotatory strength on same or next line
                combined = line
                if i + 1 < len(lines):
                    combined += lines[i + 1]
                if rot_match := re.search(rot_pattern, combined):
                    rot_strength = float(rot_match[1])
                transitions.append(
                    Transition(
                        state=int(match[1]),
                        energy_eV=float(match[2]),
                        wavelength_nm=float(match[3]),
                        osc_strength=float(match[4]),
                        rot_strength=rot_strength,
                    )
                )
        return transitions

    @staticmethod
    def parse_color(stdout: str) -> Color | None:
        """Extract predicted colour from PREDICT_COLOR output.

        Returns
        -------
        Color | None
            Color result object or None
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
        if result:
            return Color(**result)
        return None


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
    def parse(stdout: str) -> SurfaceAnalysis:
        """Extract surface analysis statistics.

        Returns
        -------
        SurfaceAnalysis
            Surface analysis result object
        """
        result: dict[str, Any] = {}

        patterns: dict[str, str] = {
            "area": (
                rf"(?:Overall|Molecular)\s+"
                rf"surface\s+area.*?:\s+({FLOAT_PATTERN})"
            ),
            "volume": (
                rf"(?:Enclosed|Molecular)\s+"
                rf"volume.*?:\s+({FLOAT_PATTERN})"
            ),
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

        return SurfaceAnalysis(**result)

    @staticmethod
    def parse_extrema(stdout: str) -> list[SurfaceExtremum]:
        """Extract surface extrema (minima and maxima).

        Returns
        -------
        list[SurfaceExtremum]
            List of surface extremum objects
        """
        extrema: list[SurfaceExtremum] = []
        # "Local  minimum   1:   -45.678  at   1.234   2.345   3.456"
        pattern = (
            rf"Local\s+(min\w*|max\w*)\s+(\d+)\s*:\s+({FLOAT_PATTERN})"
            rf"\s+at\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+"
            rf"({FLOAT_PATTERN})"
        )
        for match in re.finditer(pattern, stdout, re.IGNORECASE):
            ext_type: Literal["min", "max"] = (
                "min" if "min" in match[1].lower() else "max"
            )
            extrema.append(
                SurfaceExtremum(
                    type=ext_type,
                    index=int(match[2]),
                    value=float(match[3]),
                    x=float(match[4]),
                    y=float(match[5]),
                    z=float(match[6]),
                )
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
    ) -> list[FuzzyAtomicProperty]:
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

        # "Atomic dipole of atom  1(C ):  X= 0.12  Y= 0.34  Z= 0.56"
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

        properties = []
        for atom_id, props in atoms.items():
            properties.append(FuzzyAtomicProperty(atom_id=atom_id, **props))
        return properties

    @staticmethod
    def parse_delocalization_indices(
        stdout: str,
    ) -> list[DelocalizationIndex]:
        """Extract delocalization indices.

        Returns
        -------
        list[DelocalizationIndex]
            List of delocalization index objects
        """
        indices: list[DelocalizationIndex] = []

        di_pattern = (
            rf"Delocalization index.*?atom\s+(\d+).*?atom\s+(\d+)"
            rf".*?:\s+({FLOAT_PATTERN})"
        )

        for match in re.finditer(di_pattern, stdout, re.IGNORECASE):
            a1, a2 = int(match[1]), int(match[2])
            if a1 > a2:
                a1, a2 = a2, a1
            indices.append(
                DelocalizationIndex(
                    atom1_id=a1, atom2_id=a2, index=float(match[3])
                )
            )

        return indices

    @staticmethod
    def parse_aromaticity_index(stdout: str) -> list[AromaticityIndex]:
        """Extract aromaticity indices (PDI, FLU, FLU-pi, MCI, ITA).

        Returns
        -------
        list[AromaticityIndex]
            List of aromaticity index objects
        """
        result: list[AromaticityIndex] = []

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
                result.append(
                    AromaticityIndex(index_name=name, value=float(match[1]))
                )

        return result


# =============================================================================
# Menu 17: Basin analysis
# =============================================================================


class BasinParser(OutputParser):
    """Parser for basin analysis output.

    Handles AIM, ELF, ESP, LOL, and custom basin analyses.
    """

    @staticmethod
    def parse(stdout: str) -> list[Basin]:
        """Extract basin integration results.

        Returns
        -------
        list[Basin]
            Basin result objects with 'basin', 'population', and optionally
            'attractor_atom' and 'attractor_element'.
        """
        basins: list[Basin] = []

        # "Basin   1  attractor at atom  3(O )  population:  9.2345"
        pattern = (
            rf"Basin\s+(\d+).*?(?:attractor.*?atom\s+(\d+)\s*\(([^)]+)\))?"
            rf".*?population[=:\s]+({FLOAT_PATTERN})"
        )
        for match in re.finditer(pattern, stdout, re.IGNORECASE):
            basin: Basin = Basin(
                basin_id=int(match[1]),
                population=float(match[4]),
            )
            if match[2]:
                basin.attractor_atom = int(match[2])
                basin.attractor_element = match[3].strip()
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
                        Basin(
                            basin_id=int(m[1]),
                            population=float(m[3]),
                            attractor_element=m[2],
                        )
                    )

        return basins

    @staticmethod
    def parse_charges(stdout: str) -> list[Charge]:
        """Extract AIM/Bader charges from basin analysis.

        Returns
        -------
        list[Charge]
            List of charge objects
        """
        pattern = (
            rf"(?:AIM|Bader)\s+charge.*?atom\s+(\d+).*?:\s+({FLOAT_PATTERN})"
        )
        charges: list[Charge] = [
            Charge(
                atom_id=int(match[1]),
                charge=float(match[2]),
            )
            for match in re.finditer(pattern, stdout, re.IGNORECASE)
        ]
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
    def parse_hole_electron(stdout: str) -> HoleElectron | None:
        """Extract hole-electron analysis descriptors.

        Returns
        -------
        HoleElectron | None
            Hole-electron analysis result or None
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
        centroids: dict[str, tuple[float, float, float]] = {}
        for match in re.finditer(centroid_pat, stdout, re.IGNORECASE):
            key = f"{match[1].lower()}_centroid"
            centroids[key] = (
                float(match[2]),
                float(match[3]),
                float(match[4]),
            )

        if "hole_centroid" in centroids and "electron_centroid" in centroids:
            return HoleElectron(
                hole_id=result["H_index"],
                electron_id=result["E_index"],
                transition_index=result["t_index"],
                electron_delocalisation_index=result["EDI"],
                hole_delocalisation_index=result["HDI"],
                Sr=result["Sr"],
                d_index=result["D_index"],
                hole_centroid=centroids["hole_centroid"],
                electron_centroid=centroids["electron_centroid"],
            )
        return None

    @staticmethod
    def parse_charge_transfer(stdout: str) -> ChargeTransfer:
        """Extract charge transfer analysis results.

        Returns
        -------
        Charge transfer result data.

        """
        ct_dist = rf"CT\s+distance[=:\s]+({FLOAT_PATTERN})"
        ct_amt = (
            rf"(?:transferred|CT)\s+(?:charge|amount)[=:\s]+({FLOAT_PATTERN})"
        )

        if match := re.search(ct_dist, stdout, re.IGNORECASE):
            ct_distance = float(match[1])
        if match := re.search(ct_amt, stdout, re.IGNORECASE):
            ct_amount = float(match[1])

        result = ChargeTransfer(
            distance=ct_distance if ct_distance is not None else None,
            transfer_amount=ct_amount if ct_amount is not None else None,
        )

        # Per-fragment: "Fragment 1:  hole= 0.85  electron= 0.15"
        frag_pattern = (
            rf"Fragment\s+(\d+).*?hole[=:\s]+({FLOAT_PATTERN})"
            rf".*?electron[=:\s]+({FLOAT_PATTERN})"
        )
        fragments: list[ChargeTransferFragment] = []
        for match in re.finditer(frag_pattern, stdout, re.IGNORECASE):
            fragments.append(
                ChargeTransferFragment(
                    fragment_id=int(match[1]),
                    hole_contribution=float(match[2]),
                    electron_contribution=float(match[3]),
                )
            )

        if fragments:
            result.fragments = fragments

        return result

    @staticmethod
    def parse_delta_r(stdout: str) -> list[DeltaR]:
        """Extract Delta_r index for each excited state.

        Returns
        -------
        list[DeltaR]
            List of DeltaR result objects
        """
        pattern = (
            rf"(?:State|Excited state)\s+(\d+).*?"
            rf"Delta_?r[=:\s]+({FLOAT_PATTERN})"
        )
        results: list[DeltaR] = [
            DeltaR(state_id=int(match[1]), delta_r=float(match[2]))
            for match in re.finditer(pattern, stdout, re.IGNORECASE)
        ]
        return results

    @staticmethod
    def parse_lambda_index(stdout: str) -> list[LambdaIndex]:
        """Extract Lambda diagnostic for each excited state.

        Returns
        -------
        list[LambdaIndex]
            List of LambdaIndex result objects
        """
        pattern = (
            rf"(?:State|Excited state)\s+(\d+).*?"
            rf"Lambda[=:\s]+({FLOAT_PATTERN})"
        )
        results: list[LambdaIndex] = [
            LambdaIndex(state_id=int(match[1]), lambda_index=float(match[2]))
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
    def parse(stdout: str) -> WeakInteraction:
        """Extract summary statistics from weak interaction analysis.

        Returns
        -------
        WeakInteraction
            Weak interaction analysis result
        """
        result: dict[str, float] = {}

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

        if result["delta_g_inter"] and result["delta_g_intra"] is not None:
            interaction = WeakInteraction(
                delta_g_inter=result["delta_g_inter"],
                delta_g_intra=result["delta_g_intra"],
            )

        if result["isosurface_integral"] is not None:
            interaction.isosurface_integral = result["isosurface_integral"]

        # Cube files produced
        cubes: list[str] = []
        cubes.extend(
            match[1]
            for match in re.finditer(
                r"(\S+\.cube)\s+has been generated", stdout
            )
        )
        if cubes:
            interaction.cube_names = cubes

        return interaction


# =============================================================================
# Menu 21: Energy Decomposition Analysis (EDA)
# =============================================================================


class EDAParser(OutputParser):
    """Parser for energy decomposition analysis output.

    Handles EDA_FF, EDA_SBL, SOBEDA, and dispersion contributions.
    """

    @staticmethod
    def parse(stdout: str) -> EnergyDecompositionAnalysis:
        """Extract EDA energy components.

        Returns
        -------
        EnergyDecompositionAnalysis
            EDA analysis result
        """
        components: dict[str, float] = {}

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
                components[key] = float(match[1])

        return EnergyDecompositionAnalysis(**components)

    @staticmethod
    def parse_dispersion_contributions(
        stdout: str,
    ) -> list[DispersionContribution]:
        """Extract per-atom dispersion energy contributions.

        Returns
        -------
        list[DispersionContribution]
            List of dispersion contribution objects
        """
        pattern = (
            rf"Atom\s+(\d+).*?(?:dispersion|D[34]).*?[=:\s]+({FLOAT_PATTERN})"
        )
        contributions = [
            DispersionContribution(
                atom_id=int(match[1]), contribution=float(match[2])
            )
            for match in re.finditer(pattern, stdout, re.IGNORECASE)
        ]
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
    def parse_global_indices(stdout: str) -> Reactivity:
        """Extract global CDFT indices.

        Returns
        -------
        ReactivityIndices
            Reactivity indices result object
        """
        indices: dict[str, float] = {}

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
                indices[key] = float(match[1])

        result = Reactivity(
            chemical_potential=indices.get("chemical_potential"),
            hardness=indices.get("hardness"),
            softness=indices.get("softness"),
            electrophilicity=indices.get("electrophilicity"),
            nucleophilicity=indices.get("nucleophilicity"),
            ionization_potential=indices.get("IP"),
            electron_affinity=indices.get("EA"),
        )

        return result

    @staticmethod
    def parse_condensed_fukui(
        stdout: str,
    ) -> list[CondensedFukui]:
        """Extract condensed Fukui function values per atom.

        Returns
        -------
        list[CondensedFukui]
            List of condensed Fukui result objects
        """
        # "Atom   1(C ):  f+= 0.1234  f-= 0.0567  f0= 0.0900"
        pattern = (
            rf"Atom\s+(\d+)\s*\([^)]+\)\s*:?\s+"
            rf"f\+[=:\s]+({FLOAT_PATTERN})\s+"
            rf"f-[=:\s]+({FLOAT_PATTERN})\s+"
            rf"f0[=:\s]+({FLOAT_PATTERN})"
        )
        fukui_list: list[CondensedFukui] = [
            CondensedFukui(
                atom_id=int(match[1]),
                fukui_plus=float(match[2]),
                fukui_minus=float(match[3]),
                fukui_zero=float(match[4]),
            )
            for match in re.finditer(pattern, stdout)
        ]
        if not fukui_list:
            # Try separate patterns
            atoms_dict: dict[int, dict[str, float]] = {}
            for label, key in [
                (r"f\+", "fukui_plus"),
                (r"f-", "fukui_minus"),
                (r"f0", "fukui_zero"),
            ]:
                pat = (
                    rf"Atom\s+(\d+)\s*\([^)]+\)\s*:?\s+"
                    rf"{label}[=:\s]+({FLOAT_PATTERN})"
                )
                for m in re.finditer(pat, stdout):
                    atoms_dict.setdefault(int(m[1]), {})[key] = float(m[2])
            fukui_list = [
                CondensedFukui(atom_id=aid, **props)
                for aid, props in atoms_dict.items()
            ]

        return fukui_list

    @staticmethod
    def parse_dual_descriptor(
        stdout: str,
    ) -> list[DualDescriptor]:
        """Extract condensed dual descriptor per atom.

        Returns
        -------
        list[DualDescriptor]
            List of dual descriptor result objects
        """
        pattern = (
            rf"Atom\s+(\d+)\s*\([^)]+\)\s*:?\s+"
            rf"(?:dual|Delta_?f|f\+\s*-\s*f-)[=:\s]+({FLOAT_PATTERN})"
        )
        dd_list: list[DualDescriptor] = [
            DualDescriptor(atom_id=int(match[1]), value=float(match[2]))
            for match in re.finditer(pattern, stdout, re.IGNORECASE)
        ]
        return dd_list


# =============================================================================
# Menu 24: Polarizability
# =============================================================================


class PolarizabilityParser(OutputParser):
    """Parser for polarizability output.

    Handles parsed polarizability/hyperpolarizability tensors and SOS
    results.
    """

    @staticmethod
    def parse(stdout: str) -> Polarizability:
        """Extract polarizability tensor data.

        Returns
        -------
        dict[str, Any]
            Dictionary with 'alpha_iso', 'alpha_aniso', tensor
            components ('alpha_xx', etc.), and optionally
            'beta_total', 'gamma_total'
        """
        result = Polarizability()

        # Isotropic polarizability
        iso_pat = (
            rf"[Ii]sotropic.*?(?:alpha|polarizability)"
            rf".*?[=:\s]+({FLOAT_PATTERN})"
        )
        if match := re.search(iso_pat, stdout):
            result.isotripic = float(match[1])

        # Anisotropy
        aniso_pat = rf"[Aa]nisotropy.*?[=:\s]+({FLOAT_PATTERN})"
        if match := re.search(aniso_pat, stdout):
            result.anisotropic = float(match[1])

        # Hyperpolarizability
        beta_pat = rf"[Bb]eta.*?total[=:\s]+({FLOAT_PATTERN})"
        if match := re.search(beta_pat, stdout):
            result.beta_total = float(match[1])

        gamma_pat = rf"[Gg]amma.*?(?:total|average)[=:\s]+({FLOAT_PATTERN})"
        if match := re.search(gamma_pat, stdout):
            result.gamma_total = float(match[1])

        # Tensor components
        components: dict[str, float] = {}
        for comp in ["xx", "xy", "xz", "yy", "yz", "zz"]:
            pat = rf"alpha[_\s]*{comp}[=:\s]+({FLOAT_PATTERN})"
            if match := re.search(pat, stdout, re.IGNORECASE):
                components[f"alpha_{comp}"] = float(match[1])

        # Build tensor if we have components
        if components:
            tensor = PolarizabilityTensor(
                alpha_xx=components.get("alpha_xx", 0.0),
                alpha_xy=components.get("alpha_xy", 0.0),
                alpha_xz=components.get("alpha_xz", 0.0),
                alpha_yy=components.get("alpha_yy", 0.0),
                alpha_yz=components.get("alpha_yz", 0.0),
                alpha_zz=components.get("alpha_zz", 0.0),
            )

            result.tensor = tensor

        return result


# =============================================================================
# Menu 25: Aromaticity
# =============================================================================


class AromaticityParser(OutputParser):
    """Parser for aromaticity analysis output.

    Handles NICS, ICSS, HOMA, HOMAC, HOMER, Bird, Stanger, and AICD.
    """

    @staticmethod
    def parse(stdout: str) -> Aromaticity:
        """Extract aromaticity indices.

        Returns
        -------
        Aromaticity
            Aromaticity result object
        """
        result = Aromaticity()

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
                setattr(result, key, float(match[1]))

        return result

    @staticmethod
    def parse_nics_scan(
        stdout: str,
    ) -> NICSScan:
        """Extract NICS scan profile data.

        Returns
        -------
        NICS
            NICS scan result object
        """
        distances: list[float] = []
        values: list[float] = []
        pattern = rf"^\s*({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s*$"
        in_scan = False
        for line in stdout.split("\n"):
            if "nics" in line.lower() and "scan" in line.lower():
                in_scan = True
                continue
            if in_scan and (match := re.match(pattern, line)):
                distances.append(float(match[1]))
                values.append(float(match[2]))
        return NICSScan(distances=distances, values=values)


# =============================================================================
# Menu 6: Wavefunction info
# =============================================================================


class WavefunctionParser(OutputParser):
    """Parser for wavefunction check/modify output.

    Handles orbital info, GTF info, basis info, density matrix,
    and overlap matrix output.
    """

    @staticmethod
    def parse_orbital_info(stdout: str) -> list[Orbital]:
        """Extract orbital information.

        Returns
        -------
        list[OrbitalInfo]
            List of orbital info result objects.
        """
        # "   5   Alpha   Occ= 2.000000   E=  -0.72340 a.u.  -19.684 eV"
        pattern = (
            rf"(\d+)\s+(Alpha|Beta)\s+Occ=\s*({FLOAT_PATTERN})\s+"
            rf"E=\s*({FLOAT_PATTERN})\s*a\.u\.\s+({FLOAT_PATTERN})\s*eV"
        )
        orbitals: list[Orbital] = []

        for match in re.finditer(pattern, stdout):
            if match:
                if match[2].lower() == "alpha":
                    spin: Literal["alpha", "beta"] | None = "alpha"
                elif match[2].lower() == "beta":
                    spin = "beta"
                else:
                    spin = None

                orbitals.append(
                    Orbital(
                        orbital_id=int(match[1]),
                        spin=spin,
                        occupation=float(match[3]),
                        energy_au=float(match[4]),
                        energy_eV=float(match[5]),
                    )
                )

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
    def parse(stdout: str) -> Cube:
        """Extract cube file generation info and grid statistics.

        Returns
        -------
        Cube
            Cube operations result object

        """
        # Statistics patterns
        stat_patterns: dict[str, str] = {
            "min": rf"[Mm]inimum.*?[=:\s]+({FLOAT_PATTERN})",
            "max": rf"[Mm]aximum.*?[=:\s]+({FLOAT_PATTERN})",
            "mean": rf"[Mm]ean.*?[=:\s]+({FLOAT_PATTERN})",
            "integral": rf"[Ii]ntegral.*?[=:\s]+({FLOAT_PATTERN})",
            "std_dev": rf"[Ss]td.*?dev.*?[=:\s]+({FLOAT_PATTERN})",
        }

        # Cube files generated
        for fname in re.finditer(r"(\S+\.cube)\s+has been generated", stdout):
            grid_pat = r"Grid dimensions:\s*(\d+)\s*x\s*(\d+)\s*x\s*(\d+)"
            if grid := re.search(grid_pat, stdout):
                cube = Cube(
                    file_name=fname[1],
                    x_dim=int(grid[1]),
                    y_dim=int(grid[2]),
                    z_dim=int(grid[3]),
                )

        for key, pat in stat_patterns.items():
            if match := re.search(pat, stdout):
                setattr(cube, key, float(match[1]))

        return cube


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
    def parse_bond_lengths(stdout: str) -> list[BondLength]:
        """Extract bond length information.

        Returns
        -------
        list[BondLength]
            Geometry result object with bond lengths, angles, and dihedrals
        """
        bond_lengths: list[BondLength] = []
        # "Bond length between atom  1(C ) and  2(N ):  1.3456 Angstrom"
        bl_pattern = (
            rf"Bond length.*?atom\s+(\d+).*?(?:and|atom)\s+(\d+).*?:\s+"
            rf"({FLOAT_PATTERN})"
        )
        for match in re.finditer(bl_pattern, stdout, re.IGNORECASE):
            bond_lengths.append(
                BondLength(
                    atom1_id=int(match[1]),
                    atom2_id=int(match[2]),
                    length=float(match[3]),
                )
            )
        return bond_lengths

    @staticmethod
    def parse_bond_angles(stdout: str) -> list[BondAngle]:
        """Extract bond angle information.

        Returns
        -------
        list[BondAngle]
            List of bond angle result objects.

        """
        bond_angles: list[BondAngle] = []

        # "Angle  1-2-3 :  120.345 degree"
        ang_pattern = (
            rf"[Aa]ngle\s+(\d+)\s*-\s*(\d+)\s*-\s*(\d+)\s*:\s+"
            rf"({FLOAT_PATTERN})"
        )
        for match in re.finditer(ang_pattern, stdout):
            bond_angles.append(
                BondAngle(
                    atom1_id=int(match[1]),
                    atom2_id=int(match[2]),
                    atom3_id=int(match[3]),
                    angle=float(match[4]),
                )
            )

        return bond_angles

    @staticmethod
    def parse_dihedral_angles(stdout: str) -> list[DihedralAngle]:
        """Extract dihedral angle information.

        Returns
        -------
        list[DihedralAngle]
            List of dihedral angle result objects.

        """
        dihedrals: list[DihedralAngle] = []
        # "Dihedral  1-2-3-4 :  -45.678 degree"
        dih_pattern = (
            rf"[Dd]ihedral\s+(\d+)\s*-\s*(\d+)\s*-\s*(\d+)\s*-\s*(\d+)"
            rf"\s*:\s+({FLOAT_PATTERN})"
        )
        for match in re.finditer(dih_pattern, stdout):
            dihedrals.append(
                DihedralAngle(
                    atom1_id=int(match[1]),
                    atom2_id=int(match[2]),
                    atom3_id=int(match[3]),
                    atom4_id=int(match[4]),
                    angle=float(match[5]),
                )
            )

        return dihedrals

    @staticmethod
    def parse_dipole_moments(stdout: str) -> DipoleMoment | None:
        """Extract electric dipole moment.

        Returns
        -------
        DipoleMoment
            Dipole moment result object with x, y, z components and total
        """
        dip_pat = (
            rf"[Dd]ipole.*?X[=:\s]+({FLOAT_PATTERN})\s+"
            rf"Y[=:\s]+({FLOAT_PATTERN})\s+Z[=:\s]+({FLOAT_PATTERN})"
        )
        if match := re.search(dip_pat, stdout):
            return DipoleMoment(
                x=float(match[1]),
                y=float(match[2]),
                z=float(match[3]),
            )

        return None

    @staticmethod
    def parse_quadrupole_moments(stdout: str) -> QuadrupoleMoment | None:
        """Extract electric quadrupole moment.

        Returns
        -------
        QuadrupoleMoment
            Quadrupole moment result object.
        """
        moments: dict[str, float] = {}
        quad_components = ["XX", "XY", "XZ", "YY", "YZ", "ZZ"]
        for comp in quad_components:
            pat = rf"[Qq]uadrupole.*?{comp}[=:\s]+({FLOAT_PATTERN})"
            if m := re.search(pat, stdout):
                moments[f"quadrupole_{comp}"] = float(m[1])

        return QuadrupoleMoment(
            xx=moments.get("quadrupole_XX"),
            xy=moments.get("quadrupole_XY"),
            xz=moments.get("quadrupole_XZ"),
            yy=moments.get("quadrupole_YY"),
            yz=moments.get("quadrupole_YZ"),
            zz=moments.get("quadrupole_ZZ"),
        )

    @staticmethod
    def parse_coordination_numbers(
        stdout: str,
    ) -> list[CoordinationNumber]:
        """Extract atomic coordination numbers.

        Returns
        -------
        list[CoordinationNumber]
            List of coordination number result objects
        """
        pattern = rf"Atom\s+(\d+).*?coordination.*?[=:\s]+({FLOAT_PATTERN})"
        return [
            CoordinationNumber(
                atom_id=int(match[1]),
                coordination_number=float(match[2]),
            )
            for match in re.finditer(pattern, stdout, re.IGNORECASE)
        ]

    @staticmethod
    def parse_bla_boa(stdout: str) -> BLA_BOA:
        """Extract Bond length alternation and bond order alternation.

        Returns
        -------
        BLA_BOA
            Bond length alternation and bond order alternation values

        """
        bla = None
        boa = None

        bla_pat = rf"BLA[=:\s]+({FLOAT_PATTERN})"
        if match := re.search(bla_pat, stdout, re.IGNORECASE):
            bla = float(match[1])

        boa_pat = rf"BOA[=:\s]+({FLOAT_PATTERN})"
        if match := re.search(boa_pat, stdout, re.IGNORECASE):
            boa = float(match[1])

        return BLA_BOA(bla=bla, boa=boa)


class ParserRoute:
    """Routing from Menu enums to OutputParser classes."""

    ROUTE_TABLE = {
        # Menu 7 — charges
        Menu.HIRSHFELD_CHARGE: ChargeParser,
        Menu.VDD_POPULATION: ChargeParser,
        Menu.MULLIKEN_POPULATION: ChargeParser,
        Menu.LOWDIN_POPULATION: ChargeParser,
        Menu.SCPA_POPULATION: ChargeParser,
        Menu.STOUT_POLITZER_POPULATION: ChargeParser,
        Menu.BICKELHAUPT_POPULATION: ChargeParser,
        Menu.BECKE_CHARGE: ChargeParser,
        Menu.ADCH_CHARGE: ChargeParser,
        Menu.CHELPG_CHARGE: ChargeParser,
        Menu.MK_CHARGE: ChargeParser,
        Menu.AIM_CHARGE: ChargeParser,
        Menu.HIRSHFELD_I_CHARGE: ChargeParser,
        Menu.CM5_CHARGE: ChargeParser,
        Menu.EEM_CHARGE: ChargeParser,
        Menu.RESP_CHARGE: ChargeParser,
        Menu.GASTEIGER_CHARGE: ChargeParser,
        Menu.MBIS_CHARGE: ChargeParser,
        Menu.DDEC_CHARGE: ChargeParser,
        # Menu 8 — orbital composition
        Menu.ORBITAL_COMPOSITION_MULLIKEN: OrbitalCompositionParser,
        Menu.ORBITAL_COMPOSITION_SCPA: OrbitalCompositionParser,
        Menu.ORBITAL_COMPOSITION_STOUT_POLITZER: OrbitalCompositionParser,
        Menu.ORBITAL_COMPOSITION_FRAGMENT_MULLIKEN: OrbitalCompositionParser,
        Menu.ORBITAL_COMPOSITION_FRAGMENT_STOUT: OrbitalCompositionParser,
        Menu.ORBITAL_COMPOSITION_FRAGMENT_SCPA: OrbitalCompositionParser,
        Menu.ORBITAL_COMPOSITION_NAO: OrbitalCompositionParser,
        Menu.ORBITAL_COMPOSITION_HIRSHFELD: OrbitalCompositionParser,
        Menu.ORBITAL_COMPOSITION_BECKE: OrbitalCompositionParser,
        Menu.LOBA_OXIDATION_STATE: OrbitalCompositionParser,
        # Menu 9 — bond orders
        Menu.MAYER_BOND_ORDER: BondOrderParser,
        Menu.WIBERG_BOND_ORDER: BondOrderParser,
        Menu.MULLIKEN_BOND_ORDER: BondOrderParser,
        Menu.FUZZY_BOND_ORDER: BondOrderParser,
        Menu.LAPLACIAN_BOND_ORDER: BondOrderParser,
        Menu.IBSI_ANALYSIS: BondOrderParser,
        Menu.MULTICENTER_BOND_ORDER: BondOrderParser,
        Menu.MULTICENTER_BOND_ORDER_NAO: BondOrderParser,
        Menu.MULLIKEN_BOND_ORDER_DECOMPOSE: BondOrderParser,
        Menu.ORBITAL_PERTURBED_MAYER: BondOrderParser,
        Menu.WIBERG_DECOMPOSITION: BondOrderParser,
        Menu.AV1245_INDEX: FuzzySpaceParser,
        # Menu 2 — topology
        Menu.TOPOLOGY_SEARCH_CPS: CriticalPointParser,
        Menu.TOPOLOGY_GENERATE_PATHS: CriticalPointParser,
        Menu.TOPOLOGY_ANALYSIS_COMPLETE: CriticalPointParser,
        Menu.TOPOLOGY_ESP_ANALYSIS: CriticalPointParser,
        Menu.TOPOLOGY_LOL_ANALYSIS: CriticalPointParser,
        Menu.TOPOLOGY_ELF_ANALYSIS: CriticalPointParser,
        Menu.TOPOLOGY_LAPLACIAN_ANALYSIS: CriticalPointParser,
        Menu.TOPOLOGY_SEARCH_BCP: CriticalPointParser,
        Menu.TOPOLOGY_SEARCH_RCP: CriticalPointParser,
        Menu.TOPOLOGY_SEARCH_CCP: CriticalPointParser,
        # Menu 10 — density of states
        Menu.PLOT_DOS: DOSParser,
        Menu.PLOT_PDOS: DOSParser,
        Menu.PLOT_OPDOS: DOSParser,
        Menu.PLOT_LDOS: DOSParser,
        Menu.PLOT_PHOTOELECTRON_SPECTRUM: DOSParser,
        Menu.PLOT_COHP: DOSParser,
        # Menu 11 — spectra
        Menu.PLOT_IR_SPECTRUM: SpectrumParser,
        Menu.PLOT_RAMAN_SPECTRUM: SpectrumParser,
        Menu.PLOT_UV_VIS_SPECTRUM: SpectrumParser,
        Menu.PLOT_ECD_SPECTRUM: SpectrumParser,
        Menu.PLOT_VCD_SPECTRUM: SpectrumParser,
        Menu.PLOT_ROA_SPECTRUM: SpectrumParser,
        Menu.PLOT_NMR_SPECTRUM: SpectrumParser,
        Menu.PLOT_FLUORESCENCE_SPECTRUM: SpectrumParser,
        Menu.PLOT_PVS: SpectrumParser,
        Menu.PREDICT_COLOR: SpectrumParser,
        # Menu 12 — surface analysis
        Menu.SURFACE_ANALYSIS_ESP: SurfaceParser,
        Menu.SURFACE_ANALYSIS_ALIE: SurfaceParser,
        Menu.SURFACE_AREA_VOLUME: SurfaceParser,
        Menu.BECKE_SURFACE: SurfaceParser,
        Menu.HIRSHFELD_SURFACE: SurfaceParser,
        Menu.SURFACE_EXTREMA: SurfaceParser,
        Menu.HIRSHFELD_SURFACE_FINGERPRINT: SurfaceParser,
        # Menu 15 — fuzzy atomic space
        Menu.FUZZY_INTEGRATE_PROPERTY: FuzzySpaceParser,
        Menu.ATOMIC_DIPOLE_MOMENTS: FuzzySpaceParser,
        Menu.ATOMIC_OVERLAP_MATRIX: FuzzySpaceParser,
        Menu.ATOMIC_VOLUME_POLARIZABILITY: FuzzySpaceParser,
        Menu.LOCALIZATION_DELOCALIZATION_INDEX: FuzzySpaceParser,
        Menu.PDI_AROMATICITY: FuzzySpaceParser,
        Menu.FLU_AROMATICITY: FuzzySpaceParser,
        Menu.FLU_PI_AROMATICITY: FuzzySpaceParser,
        Menu.MULTICENTER_DI: FuzzySpaceParser,
        Menu.ITA_AROMATICITY: FuzzySpaceParser,
        Menu.CONDENSED_LINEAR_RESPONSE: FuzzySpaceParser,
        Menu.PARA_LINEAR_RESPONSE: FuzzySpaceParser,
        Menu.IFDI_ANALYSIS: FuzzySpaceParser,
        Menu.FUZZY_INTEGRATE_OVERLAP: FuzzySpaceParser,
        # Menu 17 — basin analysis
        Menu.BASIN_ANALYSIS_AIM: BasinParser,
        Menu.BASIN_ANALYSIS_ELF: BasinParser,
        Menu.BASIN_INTEGRATE_PROPERTY: BasinParser,
        Menu.BASIN_ANALYSIS_ESP: BasinParser,
        Menu.BASIN_ANALYSIS_LOL: BasinParser,
        Menu.BASIN_ANALYSIS_LOL_ALPHA: BasinParser,
        Menu.BASIN_ANALYSIS_CUSTOM: BasinParser,
        # Menu 18 — excitation analysis
        Menu.HOLE_ELECTRON_ANALYSIS: ExcitationParser,
        Menu.CHARGE_TRANSFER_ANALYSIS: ExcitationParser,
        Menu.DELTA_R_INDEX: ExcitationParser,
        Menu.LAMBDA_INDEX: ExcitationParser,
        Menu.IFCT_ANALYSIS: ExcitationParser,
        Menu.CTS_ANALYSIS: ExcitationParser,
        # Menu 20 — weak interactions
        Menu.NCI_ANALYSIS: WeakInteractionParser,
        Menu.NCI_PROMOLECULAR: WeakInteractionParser,
        Menu.ANCI_ANALYSIS: WeakInteractionParser,
        Menu.IRI_ANALYSIS: WeakInteractionParser,
        Menu.DORI_ANALYSIS: WeakInteractionParser,
        Menu.VDW_POTENTIAL: WeakInteractionParser,
        Menu.IGM_ANALYSIS: WeakInteractionParser,
        Menu.IGMH_ANALYSIS: WeakInteractionParser,
        Menu.AIGM_ANALYSIS: WeakInteractionParser,
        Menu.MIGM_ANALYSIS: WeakInteractionParser,
        Menu.AMIGM_ANALYSIS: WeakInteractionParser,
        # Menu 21 — EDA
        Menu.EDA_FF: EDAParser,
        Menu.EDA_SBL: EDAParser,
        Menu.SOBEDA_ANALYSIS: EDAParser,
        Menu.DISPERSION_ATOMIC_CONTRIBUTION: EDAParser,
        # Menu 22 — CDFT
        Menu.CDFT_ANALYSIS: CDFTParser,
        Menu.FUKUI_FUNCTION: CDFTParser,
        Menu.DUAL_DESCRIPTOR: CDFTParser,
        Menu.CONDENSED_FUKUI: CDFTParser,
        Menu.LOCAL_HARDNESS: CDFTParser,
        Menu.LOCAL_IONIZATION_ENERGY: CDFTParser,
        # Menu 24 — polarizability
        Menu.PARSE_POLARIZABILITY: PolarizabilityParser,
        Menu.SOS_POLARIZABILITY: PolarizabilityParser,
        Menu.POLARIZABILITY_DENSITY: PolarizabilityParser,
        Menu.UNIT_SPHERE_POLARIZABILITY: PolarizabilityParser,
        # Menu 25 — aromaticity
        Menu.AICD_ANALYSIS: AromaticityParser,
        Menu.NICS_POINT: AromaticityParser,
        Menu.ICSS_ANALYSIS: AromaticityParser,
        Menu.NICS_SCAN: AromaticityParser,
        Menu.BIRD_INDEX: AromaticityParser,
        Menu.HOMA_INDEX: AromaticityParser,
        Menu.HOMAC_HOMER: AromaticityParser,
        Menu.STANGER_INDEX: AromaticityParser,
        Menu.NICS_1D_SCAN: AromaticityParser,
        Menu.NICS_2D_MAP: AromaticityParser,
        # Menu 6 — wavefunction info
        Menu.PRINT_ORBITAL_INFO: WavefunctionParser,
        Menu.PRINT_GTF_INFO: WavefunctionParser,
        Menu.PRINT_BASIS_INFO: WavefunctionParser,
        # Menu 5 — cube generation (store metadata only, not the cube)
        Menu.CUBE_DENSITY: CubeParser,
        Menu.CUBE_SPIN_DENSITY: CubeParser,
        Menu.CUBE_ELF: CubeParser,
        Menu.CUBE_LOL: CubeParser,
        Menu.CUBE_ESP: CubeParser,
        Menu.CUBE_LAPLACIAN: CubeParser,
        Menu.CUBE_GRADIENT_NORM: CubeParser,
        Menu.CUBE_KINETIC_G: CubeParser,
        Menu.CUBE_KINETIC_K: CubeParser,
        Menu.CUBE_ALIE: CubeParser,
        Menu.CUBE_RDG: CubeParser,
        Menu.CUBE_SIGN_LAMBDA2_RHO: CubeParser,
        Menu.CUBE_SOURCE_FUNCTION: CubeParser,
        Menu.CUBE_ORBITAL_WAVEFUNCTION: CubeParser,
        Menu.CUBE_FUKUI_MINUS: CubeParser,
        Menu.CUBE_FUKUI_PLUS: CubeParser,
        Menu.CUBE_DUAL_DESCRIPTOR: CubeParser,
        Menu.CUBE_PROMOLECULAR_DENSITY: CubeParser,
        Menu.CUBE_DEFORMATION_DENSITY: CubeParser,
        Menu.CUBE_DENSITY_HIGH: CubeParser,
        Menu.CUBE_ESP_HIGH: CubeParser,
        Menu.CUBE_ELF_HIGH: CubeParser,
        # Menu 100 — utilities
        Menu.GEOMETRY_PROPERTIES: UtilityParser,
        Menu.ELECTRIC_MULTIPOLE_MOMENTS: UtilityParser,
        Menu.BLA_BOA_ANALYSIS: UtilityParser,
    }
