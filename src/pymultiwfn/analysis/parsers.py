"""Output parsers for Multiwfn results.

ClaudeCode-extracted regex patterns from Multiwfn 3.8 manual (6 Jan 2026).

Each parser class inherits from :class:`OutputParser` and exposes a
:meth:`parse_for_result` classmethod that the :class:`MultiwfnResult`
container calls.  ``parse_for_result`` dispatches to the appropriate
static parsing helpers and returns a flat
``list[ParsedMultiwfnResult]``.
"""

from __future__ import annotations

import re
from typing import Any, Literal

from pymultiwfn.analysis.result import (
    BLA_BOA,
    AOMDiagnostics,
    Aromaticity,
    AromaticityIndex,
    AtomicMultipole,
    AtomInfo,
    Basin,
    BasisFunction,
    BondAngle,
    BondLength,
    BondOrder,
    BondOrderDecomposition,
    BondOrderSet,
    Charge,
    ChargeSet,
    ChargeTransfer,
    ChargeTransferFragment,
    CLRKMatrix,
    CoefficientMatrix,
    Color,
    CondensedFukui,
    CoordinationNumber,
    CriticalPoint,
    Cube,
    DelocalizationIndex,
    DelocalizationIndexMatrix,
    DeltaR,
    DensityMatrix,
    DensityOfStates,
    DihedralAngle,
    Dipole,
    DipoleMoment,
    DispersionContribution,
    DOSMetadata,
    DualDescriptor,
    ElectricMultipoleMomentReport,
    EnergyDecompositionAnalysis,
    ExportedMatrix,
    FLUReferenceParameter,
    FuzzyIntegrationEntry,
    FuzzyIntegrationResult,
    GaussianTypeFunction,
    GridExtremum,
    HoleElectron,
    HOMOLUMOGap,
    IBSIAnalysis,
    IBSIEntry,
    LambdaIndex,
    LocalizationIndex,
    MolecularMultipole,
    MultiCenterBondOrder,
    NICSScan,
    Orbital,
    OrbitalAtomComposition,
    OrbitalAtomEntry,
    OrbitalBasisComposition,
    OrbitalBasisEntry,
    OrbitalEnergy,
    OrbitalShellEntry,
    OrbitalShellTypeComposition,
    OrbitalWeightDecomposition,
    OrbitalWeightedFukuiEntry,
    OrbitalWeightedFukuiResult,
    OrbitalWeightEntry,
    OverlapIntegrationMatrix,
    OxidationState,
    ParsedMultiwfnResult,
    PoincareHopfCounts,
    Polarizability,
    PolarizabilityTensor,
    QuadrupoleMoment,
    Reactivity,
    Spectrum,
    SpectrumCurveExtrema,
    SpectrumExtremum,
    SuperdelocalizabilityEntry,
    SuperdelocalizabilityResult,
    SurfaceAnalysisResult,
    SurfaceExtremum,
    SurfaceGeometry,
    SurfaceStatistics,
    TopologyPath,
    Transition,
    Valence,
    WeakInteraction,
)
from pymultiwfn.enums.menu import Menu

# Shared float pattern for all numeric parsing
FLOAT_PATTERN = r"[-+]?(?:\d+\.\d*|\d*\.\d+|\d+)(?:[Ee][-+]?\d+)?"


class OutputParser:
    """Base class for output parsers.

    Subclasses must implement :meth:`parse_for_result` which receives
    the :class:`Menu` enum and raw stdout, and returns a list of
    :class:`ParsedMultiwfnResult` instances.
    """

    # Type alias for callables used in case maps.  Each callable
    # accepts ``stdout`` and returns either a single result, a list of
    # results, or ``None``.
    _Parser = Any  # Callable[[str], ParsedMultiwfnResult | list | None]

    @classmethod
    def parse_for_result(
        cls,
        analysis: Menu,
        stdout: str,
    ) -> list[ParsedMultiwfnResult]:
        """Dispatch to the appropriate static helpers and return results.

        Subclasses override this to call their specific ``parse_*``
        methods and return a flat list.
        """
        raise NotImplementedError

    @staticmethod
    def _collect(
        stdout: str,
        parsers: list[Any],
    ) -> list[ParsedMultiwfnResult]:
        """Run *parsers* against *stdout* and flatten into a single list.

        Each entry in *parsers* is a callable ``(str) -> T`` where *T*
        may be:

        * a single :class:`ParsedMultiwfnResult` — appended directly
        * a ``list`` of results — extended into the output
        * ``None`` — skipped

        For dataclasses whose "empty" state is meaningful (e.g.
        ``Spectrum`` with no frequencies), the caller should wrap the
        parser in a lambda that returns ``None`` when empty.
        """
        results: list[ParsedMultiwfnResult] = []
        for fn in parsers:
            value = fn(stdout)
            if value is None:
                continue
            if isinstance(value, list):
                results.extend(value)
            else:
                results.append(value)
        return results


# =============================================================================
# Menu 0: View structure
# =============================================================================


class AtomListParser(OutputParser):
    """Parser for atom list and HOMO-LUMO gap output."""

    @classmethod
    def parse_for_result(
        cls,
        analysis: Menu,
        stdout: str,
    ) -> list[ParsedMultiwfnResult]:
        return cls._collect(stdout, [cls.parse_atoms, cls.parse_homo_lumo_gap])

    @staticmethod
    def parse_atoms(stdout: str) -> list[AtomInfo]:
        """Extract atom list entries."""
        pattern = (
            rf"^\s*(\d+)\s*\(([A-Za-z]+)\s*\)\s*-->\s*Charge:\s*({FLOAT_PATTERN})"
            rf"\s+x,y,z\(Bohr\):\s*({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})"
        )
        return [
            AtomInfo(
                atom_id=int(match[1]),
                element=match[2],
                nuclear_charge=float(match[3]),
                x_bohr=float(match[4]),
                y_bohr=float(match[5]),
                z_bohr=float(match[6]),
            )
            for match in re.finditer(pattern, stdout, re.MULTILINE)
        ]

    @staticmethod
    def parse_homo_lumo_gap(stdout: str) -> HOMOLUMOGap | None:
        """Extract HOMO, LUMO, and gap information."""
        homo_pattern = (
            rf"Orbital\s+(\d+)\s+is\s+HOMO.*?energy:\s*({FLOAT_PATTERN})\s*a\.u\.\s+"
            rf"({FLOAT_PATTERN})\s*eV"
        )
        lumo_pattern = (
            rf"Orbital\s+(\d+)\s+is\s+LUMO.*?energy:\s*({FLOAT_PATTERN})\s*a\.u\.\s+"
            rf"({FLOAT_PATTERN})\s*eV"
        )
        gap_pattern = (
            rf"HOMO-LUMO\s+gap:\s*({FLOAT_PATTERN})\s*a\.u\.\s+"
            rf"({FLOAT_PATTERN})\s*eV\s+({FLOAT_PATTERN})\s*kJ/mol"
        )

        homo = re.search(homo_pattern, stdout)
        lumo = re.search(lumo_pattern, stdout)
        gap = re.search(gap_pattern, stdout)

        if homo and lumo and gap:
            return HOMOLUMOGap(
                homo_index=int(homo[1]),
                homo_energy_au=float(homo[2]),
                homo_energy_eV=float(homo[3]),
                lumo_index=int(lumo[1]),
                lumo_energy_au=float(lumo[2]),
                lumo_energy_eV=float(lumo[3]),
                gap_au=float(gap[1]),
                gap_eV=float(gap[2]),
                gap_kJ_mol=float(gap[3]),
            )
        return None


# =============================================================================
# Menu 2: Topology analysis
# =============================================================================


class CriticalPointParser(OutputParser):
    """Parser for critical point information from topology analysis."""

    CP_TYPE_NAMES: dict[
        str, Literal["nuclear", "bond", "ring", "cage", "unknown"]
    ] = {
        "(3,-3)": "nuclear",
        "(3,-1)": "bond",
        "(3,+1)": "ring",
        "(3,+3)": "cage",
    }

    @classmethod
    def parse_for_result(
        cls,
        analysis: Menu,
        stdout: str,
    ) -> list[ParsedMultiwfnResult]:
        results: list[ParsedMultiwfnResult] = []
        results.extend(cls.parse(stdout))
        results.extend(cls.parse_topology_paths(stdout))
        results.extend(cls.parse_bond_paths(stdout))
        ph = cls.parse_poincare_hopf(stdout)
        if ph is not None:
            results.append(ph)
        return results

    @classmethod
    def parse(cls, stdout: str) -> list[CriticalPoint]:
        """Extract critical points from both short and long summaries.

        The long summary (from path generation) includes nucleus and
        bond labels.  The short summary (from initial search) does
        not.  If both are present the long summary takes precedence.
        """
        cps: list[CriticalPoint] = []

        # ── Long summary ──
        # "  1    0.000   -4.680   -0.000   (3,-3)   Nucleus:   10(H )"
        # " 13    0.000   -3.966   -0.000   (3,-1)   10(H ) --    4(C )"
        # " 19    0.000    0.000   -0.000   (3,+1)"
        long_pattern = (
            rf"^\s*(\d+)\s+"
            rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+"
            rf"\((\d+),([+-]?\d+)\)"
            rf"(?:\s+Nucleus:\s*(\d+)\s*\(([A-Za-z]+)\s*\))?"
            rf"(?:\s+(\d+)\s*\([A-Za-z]+\s*\)\s*--\s*(\d+)\s*\([A-Za-z]+\s*\))?"
        )

        in_long_summary = False
        for line in stdout.split("\n"):
            if "Summary of found CPs" in line:
                in_long_summary = True
                cps.clear()
                continue
            if in_long_summary and (
                "number of critical points" in line.lower()
            ):
                in_long_summary = False
                continue

            if in_long_summary:
                if match := re.match(long_pattern, line):
                    cp_type_str = f"({match[5]},{match[6]})"
                    cp = CriticalPoint(
                        index=int(match[1]),
                        x=float(match[2]),
                        y=float(match[3]),
                        z=float(match[4]),
                        type=cls.CP_TYPE_NAMES.get(cp_type_str, "unknown"),
                    )
                    if match[7]:
                        cp.nucleus_atom_id = int(match[7])
                        cp.nucleus_element = match[8].strip()
                    if match[9]:
                        cp.bonded_atom1_id = int(match[9])
                        cp.bonded_atom2_id = int(match[10])
                    cps.append(cp)

        if cps:
            return cps

        # ── Short summary (fallback) ──
        short_pattern = (
            rf"^\s*(\d+)\s+"
            rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+"
            rf"\((\d+),([+-]?\d+)\)"
        )
        in_short_summary = False
        for line in stdout.split("\n"):
            if re.search(r"----\s*Summary\s*----", line):
                in_short_summary = True
                continue
            if in_short_summary and "Totally find" in line:
                in_short_summary = False
                continue

            if in_short_summary:
                if match := re.match(short_pattern, line):
                    cp_type_str = f"({match[5]},{match[6]})"
                    cps.append(
                        CriticalPoint(
                            index=int(match[1]),
                            x=float(match[2]),
                            y=float(match[3]),
                            z=float(match[4]),
                            type=cls.CP_TYPE_NAMES.get(cp_type_str, "unknown"),
                        )
                    )

        return cps

    @classmethod
    def parse_topology_paths(cls, stdout: str) -> list[TopologyPath]:
        """Extract topology gradient paths from path summary."""
        pattern = (
            rf"Path\s+(\d+),\s+CP:\s+(\d+)\s+\((\d+),([+-]?\d+)\)\s+"
            rf"-->\s+CP:\s+(\d+)\s+\((\d+),([+-]?\d+)\)\s+"
            rf"Length:\s+({FLOAT_PATTERN})"
        )
        return [
            TopologyPath(
                path_id=int(match[1]),
                bcp_index=int(match[2]),
                bcp_type=cls.CP_TYPE_NAMES.get(
                    f"({match[3]},{match[4]})", "unknown"
                ),
                target_cp_index=int(match[5]),
                target_cp_type=cls.CP_TYPE_NAMES.get(
                    f"({match[6]},{match[7]})", "unknown"
                ),
                length_bohr=float(match[8]),
            )
            for match in re.finditer(pattern, stdout)
        ]

    @staticmethod
    def parse_poincare_hopf(stdout: str) -> PoincareHopfCounts | None:
        """Extract Poincare-Hopf critical point counts."""
        count_pattern = (
            r"\(3,-3\):\s*(\d+).*?"
            r"\(3,-1\):\s*(\d+).*?"
            r"\(3,\+1\):\s*(\d+).*?"
            r"\(3,\+3\):\s*(\d+)"
        )
        if match := re.search(count_pattern, stdout):
            satisfied = bool(
                re.search(r"Poincare-Hopf relationship is satisfied", stdout)
            )
            return PoincareHopfCounts(
                nuclear=int(match[1]),
                bond=int(match[2]),
                ring=int(match[3]),
                cage=int(match[4]),
                satisfied=satisfied,
            )
        return None

    @staticmethod
    def parse_bond_paths(stdout: str) -> list[TopologyPath]:
        """Extract bond path information."""
        pattern = (
            rf"Bond path between atom\s+(\d+).*?and atom\s+(\d+).*?"
            rf"BCP\s+(\d+).*?length\s+({FLOAT_PATTERN})"
        )
        return [
            TopologyPath(
                atom1_id=int(match[1]),
                atom2_id=int(match[2]),
                bcp_id=int(match[3]),
                path_length=float(match[4]),
            )
            for match in re.finditer(pattern, stdout, re.IGNORECASE)
        ]


# =============================================================================
# Menu 5 / 13: Cube / grid operations
# =============================================================================
# ── parsers.py ──


class CubeParser(OutputParser):
    """Parser for cube generation and grid processing output.

    Each cube export sequence in the stdout produces a separate
    ``Cube`` instance.  The parser splits on cube file markers so
    that per-sequence grid metadata, statistics, and the output
    filename are kept together.
    """

    @classmethod
    def parse_for_result(
        cls,
        analysis: Menu,
        stdout: str,
    ) -> list[ParsedMultiwfnResult]:
        return cls.parse(stdout)

    @classmethod
    def parse(cls, stdout: str) -> list[Cube]:
        """Extract one ``Cube`` per exported cube file sequence."""
        # Split stdout into per-sequence chunks by cube file markers.
        file_pattern = r"(\S+\.cub[e]?)\s+in current folder"
        file_matches = list(re.finditer(file_pattern, stdout))

        if not file_matches:
            return []

        # Build chunk boundaries: each chunk ends at a file marker
        # and starts either at the beginning or after the previous
        # file marker.
        chunks: list[tuple[str, str]] = []
        for i, fm in enumerate(file_matches):
            start = file_matches[i - 1].end() if i > 0 else 0
            end = fm.end()
            chunks.append((stdout[start:end], fm[1]))

        results: list[Cube] = []
        for chunk, file_name in chunks:
            results.append(cls._parse_chunk(chunk, file_name))

        return results

    @classmethod
    def _parse_chunk(cls, chunk: str, file_name: str) -> Cube:
        """Parse a single sequence chunk into a ``Cube``."""
        cube = Cube(file_name=file_name)

        # ── Origin ──
        origin_pat = (
            rf"Coordinate of origin.*?X,Y,Z\s+is\s+"
            rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})"
        )
        if match := re.search(origin_pat, chunk):
            cube.origin_x = float(match[1])
            cube.origin_y = float(match[2])
            cube.origin_z = float(match[3])

        # ── End point ──
        end_pat = (
            rf"Coordinate of end point.*?X,Y,Z\s+is\s+"
            rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})"
        )
        if match := re.search(end_pat, chunk):
            cube.end_x = float(match[1])
            cube.end_y = float(match[2])
            cube.end_z = float(match[3])

        # ── Spacing ──
        spacing_pat = (
            rf"Grid spacing.*?X,Y,Z\s+is\s+"
            rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})"
        )
        if match := re.search(spacing_pat, chunk):
            cube.spacing_x = float(match[1])
            cube.spacing_y = float(match[2])
            cube.spacing_z = float(match[3])

        # ── Grid dimensions ──
        dim_pat = (
            r"Number of points.*?X,Y,Z\s+is\s+"
            r"(\d+)\s+(\d+)\s+(\d+)\s+Total:\s+(\d+)"
        )
        if match := re.search(dim_pat, chunk):
            cube.x_dim = int(match[1])
            cube.y_dim = int(match[2])
            cube.z_dim = int(match[3])
            cube.total_points = int(match[4])

        # ── Minimum ──
        min_pat = (
            rf"The minimum is\s+({FLOAT_PATTERN})\s+at\s+"
            rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})"
        )
        if match := re.search(min_pat, chunk):
            cube.minimum = GridExtremum(
                value=float(match[1]),
                x_bohr=float(match[2]),
                y_bohr=float(match[3]),
                z_bohr=float(match[4]),
            )

        # ── Maximum ──
        max_pat = (
            rf"The maximum is\s+({FLOAT_PATTERN})\s+at\s+"
            rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})"
        )
        if match := re.search(max_pat, chunk):
            cube.maximum = GridExtremum(
                value=float(match[1]),
                x_bohr=float(match[2]),
                y_bohr=float(match[3]),
                z_bohr=float(match[4]),
            )

        # ── Integrals ──
        total_pat = (
            rf"Summing up all value.*?multiply differential element:\s+"
            rf"({FLOAT_PATTERN})"
        )
        if match := re.search(total_pat, chunk, re.DOTALL):
            cube.integral_total = float(match[1])

        pos_pat = (
            rf"Summing up positive value.*?multiply differential element:\s+"
            rf"({FLOAT_PATTERN})"
        )
        if match := re.search(pos_pat, chunk, re.DOTALL):
            cube.integral_positive = float(match[1])

        neg_pat = (
            rf"Summing up negative value.*?multiply differential element:\s+"
            rf"({FLOAT_PATTERN})"
        )
        if match := re.search(neg_pat, chunk, re.DOTALL):
            cube.integral_negative = float(match[1])

        return cube


# =============================================================================
# Menu 6: Wavefunction info
# =============================================================================


class WavefunctionParser(OutputParser):
    """Parser for wavefunction check/modify output (Menu 6)."""

    @classmethod
    def parse_for_result(
        cls,
        analysis: Menu,
        stdout: str,
    ) -> list[ParsedMultiwfnResult]:
        case_map: dict[Menu, list] = {}
        default = [
            cls.parse_orbital_info,
            cls.parse_gtf_info,
            cls.parse_basis_info,
            cls.parse_coefficient_matrices,
            cls.parse_density_matrices,
            cls.parse_exported_matrices,
        ]
        parsers = case_map.get(analysis, default)
        return cls._collect(stdout, parsers)

    @staticmethod
    def parse_orbital_info(stdout: str) -> list[Orbital]:
        """Extract orbital information.

        Handles both formats:
        - "N  Alpha/Beta  Occ=  E=  a.u.  eV"
        - "Orb: N Ene(au/eV): E_au E_eV Occ: occ Type: type"
        """
        orbitals: list[Orbital] = []

        # Format 1: detailed listing
        pattern1 = (
            rf"(\d+)\s+(Alpha|Beta)\s+Occ=\s*({FLOAT_PATTERN})\s+"
            rf"E=\s*({FLOAT_PATTERN})\s*a\.u\.\s+({FLOAT_PATTERN})\s*eV"
        )
        for match in re.finditer(pattern1, stdout):
            spin: Literal["alpha", "beta"] | None = None
            if match[2].lower() == "alpha":
                spin = "alpha"
            elif match[2].lower() == "beta":
                spin = "beta"
            orbitals.append(
                Orbital(
                    orbital_id=int(match[1]),
                    spin=spin,
                    occupation=float(match[3]),
                    energy_au=float(match[4]),
                    energy_eV=float(match[5]),
                )
            )

        if orbitals:
            return orbitals

        # Format 2: "Orb: N Ene(au/eV): E_au E_eV Occ: occ Type: type"
        pattern2 = (
            rf"Orb:\s+(\d+)\s+Ene\(au/eV\):\s+({FLOAT_PATTERN})\s+"
            rf"({FLOAT_PATTERN})\s+Occ:\s+({FLOAT_PATTERN})\s+"
            rf"Type:\s*(\S+)"
        )
        for match in re.finditer(pattern2, stdout):
            orbitals.append(
                Orbital(
                    orbital_id=int(match[1]),
                    energy_au=float(match[2]),
                    energy_eV=float(match[3]),
                    occupation=float(match[4]),
                )
            )

        return orbitals

    @staticmethod
    def parse_gtf_info(stdout: str) -> list[GaussianTypeFunction]:
        """Extract Gaussian-type function listing."""
        pattern = (
            rf"(\d+)\s+Center:\s+(\d+)\s*\(([A-Za-z]+)\s*\)\s+"
            rf"Type:\s+(\S+)\s+Exponent:\s+({FLOAT_PATTERN})"
        )
        return [
            GaussianTypeFunction(
                gtf_index=int(match[1]),
                center_atom_id=int(match[2]),
                center_element=match[3].strip(),
                function_type=match[4].strip(),
                exponent=float(match[5]),
            )
            for match in re.finditer(pattern, stdout)
        ]

    @staticmethod
    def parse_basis_info(stdout: str) -> list[BasisFunction]:
        """Extract basis function to shell/GTF mapping."""
        pattern = (
            r"Basis:\s+(\d+)\s+Shell:\s+(\d+)\s+Center:\s+(\d+)"
            r"\s*\(([A-Za-z]+)\s*\)\s+Type:\s*(\S+)\s+"
            r"GTF:\s+(\d+)\s+to\s+(\d+)"
        )
        return [
            BasisFunction(
                basis_index=int(match[1]),
                shell_index=int(match[2]),
                center_atom_id=int(match[3]),
                center_element=match[4].strip(),
                function_type=match[5].strip(),
                gtf_start=int(match[6]),
                gtf_end=int(match[7]),
            )
            for match in re.finditer(pattern, stdout)
        ]

    @staticmethod
    def parse_coefficient_matrices(
        stdout: str,
    ) -> list[CoefficientMatrix]:
        """Extract MO coefficient matrices printed in column blocks.

        Multiple matrices may appear in a single stdout (e.g. alpha
        and beta).  Each is delimited by a header or section break.
        The block format prints 5 orbital columns at a time with row
        indices on the left.
        """
        results: list[CoefficientMatrix] = []

        # Split into sections by density matrix headers or end markers
        # to avoid mixing coefficient and density data.
        # Coefficient data has no "density matrix" header — it's just
        # rows of basis-index followed by floats, preceded by a
        # column-header row of orbital indices.
        col_header_pat = re.compile(r"^\s+(\d+(?:\s+\d+)+)\s*$")
        row_pat = re.compile(rf"^\s+(\d+)((?:\s+{FLOAT_PATTERN})+)\s*$")

        # Don't parse inside density matrix sections
        in_density = False
        current_cols: list[int] = []
        data: dict[tuple[int, int], float] = {}
        max_row = 0
        max_col = 0

        for line in stdout.split("\n"):
            if "density matrix" in line.lower():
                in_density = True
                continue
            if in_density:
                if "Trace of density matrix" in line:
                    in_density = False
                continue

            if "Basic information of all orbitals" in line:
                continue
            if "Orb:" in line or "Largest exponent" in line:
                continue
            if "Center:" in line and "Type:" in line and "Exponent:" in line:
                continue
            if "Basis:" in line and "Shell:" in line and "GTF:" in line:
                continue

            if match := col_header_pat.match(line):
                vals = match[1].split()
                # Column headers are orbital indices (usually 1-based)
                # Only treat as column header if values are sequential-ish
                # integers and there are <= 10 of them
                if len(vals) <= 10:
                    try:
                        cols = [int(v) for v in vals]
                        if all(c > 0 for c in cols):
                            current_cols = cols
                            continue
                    except ValueError:
                        pass

            if current_cols and (match := row_pat.match(line)):
                row_idx = int(match[1])
                vals = [float(v) for v in match[2].split()]
                for k, val in enumerate(vals):
                    if k < len(current_cols):
                        col_idx = current_cols[k]
                        data[(row_idx, col_idx)] = val
                        if row_idx > max_row:
                            max_row = row_idx
                        if col_idx > max_col:
                            max_col = col_idx

        if data:
            results.append(
                CoefficientMatrix(
                    label="coefficient_matrix",
                    n_basis=max_row,
                    n_orbitals=max_col,
                    data=data,
                )
            )

        return results

    @staticmethod
    def parse_density_matrices(stdout: str) -> list[DensityMatrix]:
        """Extract density matrices (symmetric, lower triangle).

        Each density matrix section starts with a header like
        ``"*** Total density matrix ***"`` and ends at the trace
        lines.  Multiple sections produce multiple results.
        """
        results: list[DensityMatrix] = []

        dm_header_pat = re.compile(
            r"\*+\s*(.*?density matrix.*?)\s*\*+", re.IGNORECASE
        )
        col_header_pat = re.compile(r"^\s+(\d+(?:\s+\d+)*)\s*$")
        row_pat = re.compile(rf"^\s+(\d+)((?:\s+{FLOAT_PATTERN})+)\s*$")
        trace_pat = re.compile(
            rf"Trace of density matrix:\s+({FLOAT_PATTERN})"
        )
        trace_overlap_pat = re.compile(
            rf"Trace of density matrix multiplied by overlap matrix:\s+"
            rf"({FLOAT_PATTERN})"
        )

        in_dm = False
        label = ""
        current_cols: list[int] = []
        data: dict[tuple[int, int], float] = {}
        max_idx = 0
        trace: float | None = None
        trace_overlap: float | None = None

        def _flush() -> None:
            nonlocal data, max_idx, trace, trace_overlap, label
            if data:
                results.append(
                    DensityMatrix(
                        label=label,
                        n_basis=max_idx,
                        data=dict(data),
                        trace=trace,
                        trace_overlap=trace_overlap,
                    )
                )
            data = {}
            max_idx = 0
            trace = None
            trace_overlap = None
            label = ""

        for line in stdout.split("\n"):
            if match := dm_header_pat.search(line):
                if in_dm:
                    _flush()
                in_dm = True
                label = match[1].strip()
                current_cols = []
                continue

            if not in_dm:
                continue

            if match := trace_pat.search(line):
                trace = float(match[1])
                continue
            if match := trace_overlap_pat.search(line):
                trace_overlap = float(match[1])
                _flush()
                in_dm = False
                continue

            if match := col_header_pat.match(line):
                current_cols = [int(v) for v in match[1].split()]
                continue

            if current_cols and (match := row_pat.match(line)):
                row_idx = int(match[1])
                vals = [float(v) for v in match[2].split()]
                for k, val in enumerate(vals):
                    if k < len(current_cols):
                        col_idx = current_cols[k]
                        data[(row_idx, col_idx)] = val
                        if row_idx > max_idx:
                            max_idx = row_idx
                        if col_idx > max_idx:
                            max_idx = col_idx

        if in_dm and data:
            _flush()

        return results

    @staticmethod
    def parse_exported_matrices(stdout: str) -> list[ExportedMatrix]:
        """Extract matrix export file references."""
        pattern = (
            r"(?:matrix|The matrix)\s+has been\s+(?:outputted|exported)\s+"
            r"to\s+(\S+)\s+in current folder"
        )
        return [
            ExportedMatrix(file_name=match[1])
            for match in re.finditer(pattern, stdout, re.IGNORECASE)
        ]


# =============================================================================
# Menu 7: Population analysis & atomic charges
# =============================================================================


class ChargeParser(OutputParser):
    """Parser for atomic charges (Menu 7).

    Each distinct charge block is stored as a ``ChargeSet`` with a
    standardized ``method`` name and optional ``stage``, so users
    can filter by method unambiguously.
    """

    MENU_METHOD: dict[Menu, str] = {
        Menu.HIRSHFELD_CHARGE: "hirshfeld",
        Menu.VDD_POPULATION: "vdd",
        Menu.MULLIKEN_POPULATION: "mulliken",
        Menu.LOWDIN_POPULATION: "lowdin",
        Menu.SCPA_POPULATION: "scpa",
        Menu.STOUT_POLITZER_POPULATION: "stout_politzer",
        Menu.BICKELHAUPT_POPULATION: "bickelhaupt",
        Menu.BECKE_CHARGE: "becke",
        Menu.ADCH_CHARGE: "adch",
        Menu.CHELPG_CHARGE: "chelpg",
        Menu.MK_CHARGE: "mk",
        Menu.CM5_CHARGE: "cm5",
        Menu.EEM_CHARGE: "eem",
        Menu.RESP_CHARGE: "resp",
        Menu.GASTEIGER_CHARGE: "gasteiger",
        Menu.MBIS_CHARGE: "mbis",
    }
    _HEADER_METHOD: list[tuple[re.Pattern[str], str]] = []

    @classmethod
    def parse_for_result(
        cls,
        analysis: Menu,
        stdout: str,
    ) -> list[ParsedMultiwfnResult]:
        results: list[ParsedMultiwfnResult] = []
        base_method = cls.MENU_METHOD.get(analysis, "unknown")
        results.extend(cls.parse_charge_sets(stdout, base_method))
        if (dipole := cls.parse_dipole(stdout)) is not None:
            results.append(dipole)
        return results

    @classmethod
    def _detect_method_from_context(
        cls, context_lines: list[str]
    ) -> str | None:
        """Scan recent context lines to detect the charge method."""
        for line in reversed(context_lines):
            for pat, method in cls._HEADER_METHOD:
                if pat.search(line):
                    return method
        return None

    @classmethod
    def parse_charge_sets(
        cls,
        stdout: str,
        base_method: str = "unknown",
    ) -> list[ChargeSet]:
        """Extract all distinct charge blocks with method labels.

        Each block produces a separate ``ChargeSet`` with a
        standardized ``method`` and ``stage``.
        """
        sets: list[ChargeSet] = []

        # Track recent context for method detection
        context: list[str] = []
        current_method = base_method
        current_stage: str | None = None
        current_charges: list[Charge] = []
        current_total: float | None = None

        # ── Header patterns ──
        final_header = re.compile(r"Final atomic charges:", re.I)
        adc_header = re.compile(
            r"-+\s*(.*?(?:ADC|ADCH)\s+.*?charges)\s*-+", re.I
        )
        cm5_header = re.compile(r"-+\s*(.*?CM5\s+charges)\s*-+", re.I)
        resp_stage_pat = re.compile(
            r"\*+\s*(Stage\s+(\d+).*?RESP.*?)\s*$", re.I
        )
        center_header = re.compile(r"^\s*Center\s+Charge\s*$", re.I)
        atom_header = re.compile(r"^\s*Atom\s+Charge\s*$", re.I)
        pop_atoms_header = re.compile(r"Population of atoms:", re.I)

        # ── Charge line patterns ──
        # "Atom    1(C ):    -0.04351137"
        atom_colon_pat = re.compile(
            rf"Atom\s+(\d+)\s*\([A-Za-z]+\s*\)\s*:\s+({FLOAT_PATTERN})"
        )
        # "Atom:    1C   Corrected charge:   -0.041271"
        corrected_pat = re.compile(
            rf"Atom:\s+(\d+)[A-Za-z]+\s+Corrected\s+charge:\s+"
            rf"({FLOAT_PATTERN})",
            re.I,
        )
        # "Atom:    1C   CM5 charge:   -0.097209"
        cm5_pat = re.compile(
            rf"Atom:\s+(\d+)[A-Za-z]+\s+CM5\s+charge:\s+"
            rf"({FLOAT_PATTERN})",
            re.I,
        )
        # "     1(C )  -0.0912959843"
        center_pat = re.compile(
            rf"^\s+(\d+)\s*\([A-Za-z]+\s*\)\s+({FLOAT_PATTERN})\s*$"
        )
        # "Atom     1(C )  Population:  6.128  Net charge: -0.128"
        pop_charge_pat = re.compile(
            rf"Atom\s+(\d+)\s*\([A-Za-z]+\s*\)\s+Population:\s+"
            rf"{FLOAT_PATTERN}\s+(?:Atomic charge|Net charge):\s+"
            rf"({FLOAT_PATTERN})"
        )

        # ── Total patterns ──
        total_pat = re.compile(
            rf"(?:Total\s+(?:net\s+)?charge|Sum\s+of\s+charges|"
            rf"Total\s+charge)\s*:\s+({FLOAT_PATTERN})",
            re.I,
        )
        sum_pat = re.compile(
            rf"Summing up all.*?charges:\s+({FLOAT_PATTERN})", re.I
        )

        # State tracking
        in_center_block = False
        in_atom_block = False
        in_pop_atoms = False

        def _flush_current() -> None:
            nonlocal current_charges, current_total
            nonlocal current_method, current_stage
            if current_charges:
                sets.append(
                    ChargeSet(
                        method=current_method,
                        stage=current_stage,
                        charges=list(current_charges),
                        total_charge=current_total,
                    )
                )
            current_charges = []
            current_total = None
            current_stage = None

        for line in stdout.split("\n"):
            context.append(line)
            if len(context) > 20:
                context.pop(0)

            # ── Detect new section headers ──

            # ADC/ADCH header
            if match := adc_header.search(line):
                _flush_current()
                header_text = match[1].lower()
                if "adch" in header_text:
                    current_method = "adch"
                else:
                    detected = cls._detect_method_from_context(context[:-1])
                    if detected == "becke":
                        current_method = "becke"
                    else:
                        current_method = detected or base_method
                current_stage = "corrected"
                in_center_block = False
                in_atom_block = False
                in_pop_atoms = False
                continue

            # CM5 header
            if cm5_header.search(line):
                _flush_current()
                current_method = "cm5"
                current_stage = "raw"
                in_center_block = False
                in_atom_block = False
                in_pop_atoms = False
                continue

            # RESP stage header
            if match := resp_stage_pat.search(line):
                _flush_current()
                current_method = "resp"
                current_stage = f"stage_{match[2]}"
                in_center_block = False
                in_atom_block = False
                continue

            # "Final atomic charges:" header
            if final_header.search(line):
                _flush_current()
                detected = cls._detect_method_from_context(context[:-1])
                current_method = detected or base_method
                current_stage = "final"
                in_center_block = False
                in_atom_block = False
                in_pop_atoms = False
                continue

            # "Center  Charge" header (CHELPG/MK/RESP)
            if center_header.match(line):
                _flush_current()
                detected = cls._detect_method_from_context(context[:-1])
                current_method = detected or base_method
                current_stage = "esp_fit"
                in_center_block = True
                in_atom_block = False
                in_pop_atoms = False
                continue

            # "Atom  Charge" header (Gasteiger)
            if atom_header.match(line):
                _flush_current()
                detected = cls._detect_method_from_context(context[:-1])
                current_method = detected or base_method
                current_stage = "raw"
                in_atom_block = True
                in_center_block = False
                in_pop_atoms = False
                continue

            # "Population of atoms:" header
            if pop_atoms_header.search(line):
                _flush_current()
                current_method = base_method
                current_stage = "population"
                in_pop_atoms = True
                in_center_block = False
                in_atom_block = False
                continue

            # ── Parse charge lines ──

            # Center block format
            if in_center_block:
                if match := center_pat.match(line):
                    current_charges.append(
                        Charge(
                            atom_id=int(match[1]),
                            charge=float(match[2]),
                        )
                    )
                    continue
                if match := total_pat.search(line):
                    current_total = float(match[1])
                    in_center_block = False
                    continue
                if match := sum_pat.search(line):
                    current_total = float(match[1])
                    in_center_block = False
                    continue

            # Atom block format (Gasteiger)
            if in_atom_block:
                if match := center_pat.match(line):
                    current_charges.append(
                        Charge(
                            atom_id=int(match[1]),
                            charge=float(match[2]),
                        )
                    )
                    continue
                if match := total_pat.search(line):
                    current_total = float(match[1])
                    in_atom_block = False
                    continue

            # Population of atoms section
            if in_pop_atoms:
                if match := pop_charge_pat.search(line):
                    current_charges.append(
                        Charge(
                            atom_id=int(match[1]),
                            charge=float(match[2]),
                        )
                    )
                    continue
                if match := total_pat.search(line):
                    current_total = float(match[1])
                    in_pop_atoms = False
                    continue

            # ADC corrected charges
            if match := corrected_pat.search(line):
                current_charges.append(
                    Charge(
                        atom_id=int(match[1]),
                        charge=float(match[2]),
                    )
                )
                continue

            # CM5 charges
            if match := cm5_pat.search(line):
                current_charges.append(
                    Charge(
                        atom_id=int(match[1]),
                        charge=float(match[2]),
                    )
                )
                continue

            # "Atom N(X): charge" format (final charges)
            if match := atom_colon_pat.search(line):
                if not in_center_block and not in_atom_block:
                    current_charges.append(
                        Charge(
                            atom_id=int(match[1]),
                            charge=float(match[2]),
                        )
                    )
                continue

            # Totals outside blocks
            if match := total_pat.search(line):
                current_total = float(match[1])
                continue
            if match := sum_pat.search(line):
                current_total = float(match[1])
                continue

        _flush_current()
        return sets

    @staticmethod
    def parse_dipole(stdout: str) -> Dipole | None:
        """Extract molecular dipole moment from charge output."""
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

    @classmethod
    def parse_charges(cls, stdout: str) -> list[Charge]:
        """Backward-compatible flat charge parser."""
        sets = cls.parse_charge_sets(stdout, base_method="unknown")
        if not sets:
            return []
        return sets[-1].charges


# =============================================================================
# Menu 8: Orbital composition analysis
# =============================================================================


class OrbitalCompositionParser(OutputParser):
    """Parser for orbital composition analysis (Menu 8).

    Produces ``OrbitalBasisComposition`` for methods that give
    per-basis detail (Mulliken/SCPA/Stout-Politzer), and
    ``OrbitalAtomComposition`` for methods that give only per-atom
    contributions (Hirshfeld/Becke/fragment).

    Each orbital is stored separately with its orbital_id, energy,
    occupation, and type so that multiple orbitals in one stdout
    do not overwrite each other.
    """

    @classmethod
    def parse_for_result(
        cls,
        analysis: Menu,
        stdout: str,
    ) -> list[ParsedMultiwfnResult]:
        if analysis == Menu.LOBA_OXIDATION_STATE:
            return cls._collect(stdout, [cls.parse_oxidation_states])

        # Try full basis composition first
        basis_results = cls.parse_basis_compositions(stdout, method)
        if basis_results:
            return basis_results

        # Fall back to atom-only composition
        atom_results = cls.parse_atom_compositions(stdout, method)
        if atom_results:
            return atom_results

        return []

    # ── Full basis composition (Mulliken/SCPA/Stout-Politzer) ────────

    @staticmethod
    def parse_basis_compositions(
        stdout: str,
        method: str,
    ) -> list[OrbitalBasisComposition]:
        """Extract per-orbital full basis compositions.

        Each orbital block starts with an orbital header, followed by
        basis contributions, shell compositions, shell-type
        composition, atom compositions, and a delocalization index.
        """
        results: list[OrbitalBasisComposition] = []

        # Orbital header:
        # "Orbital:    21  Energy(a.u.):     -0.246939  Occ:  2.000000  Type: Alpha&Beta"  # noqa: E501
        orb_header_pat = re.compile(
            rf"Orbital:\s+(\d+)\s+Energy\(a\.u\.\):\s+({FLOAT_PATTERN})"
            rf"\s+Occ:\s+({FLOAT_PATTERN})\s+Type:\s+(\S+)"
        )

        # Basis contribution:
        # "    20   Z        2(C )    9      8.67678 %      5.61233 %     14.28911 %"  # noqa: E501
        basis_pat = re.compile(
            rf"^\s+(\d+)\s+([A-Z]+)\s+(\d+)\s*\(([A-Za-z]+)\s*\)\s+"
            rf"(\d+)\s+({FLOAT_PATTERN})\s*%\s+({FLOAT_PATTERN})\s*%"
            rf"\s+({FLOAT_PATTERN})\s*%"
        )

        # Shell composition:
        # "Shell     9 Type: P    in atom    2(C ) :    14.28911 %"
        shell_pat = re.compile(
            rf"Shell\s+(\d+)\s+Type:\s+(\S+)\s+in\s+atom\s+(\d+)"
            rf"\s*\(([A-Za-z]+)\s*\)\s*:\s+({FLOAT_PATTERN})\s*%"
        )

        # Shell-type composition:
        # "s:   0.000  p:  99.010  d:   0.990  f:   0.000  g:   0.000  h:   0.000"  # noqa: E501
        shell_type_pat = re.compile(
            rf"^\s*s:\s+({FLOAT_PATTERN})\s+p:\s+({FLOAT_PATTERN})\s+"
            rf"d:\s+({FLOAT_PATTERN})\s+f:\s+({FLOAT_PATTERN})\s+"
            rf"g:\s+({FLOAT_PATTERN})\s+h:\s+({FLOAT_PATTERN})"
        )

        # Atom composition:
        # "Atom     1(C ) :     0.30819 %"
        atom_pat = re.compile(
            rf"Atom\s+(\d+)\s*\(([A-Za-z]+)\s*\)\s*:\s+"
            rf"({FLOAT_PATTERN})\s*%"
        )

        # Delocalization index:
        # "Orbital delocalization index:   24.69"
        deloc_pat = re.compile(
            rf"Orbital delocalization index:\s+({FLOAT_PATTERN})"
        )

        # "Composition of each shell" / "Composition of different types"
        # / "Composition of each atom:" are section markers
        section_shell = re.compile(r"Composition of each shell")
        section_type = re.compile(r"Composition of different types of shells")
        section_atom = re.compile(r"Composition of each atom:")

        current: OrbitalBasisComposition | None = None
        in_basis = False
        in_shell = False
        in_type = False
        in_atom = False

        for line in stdout.split("\n"):
            # New orbital header
            if match := orb_header_pat.search(line):
                # Check if this is followed by basis data (not just
                # a listing).  We can't know yet, so start a new one
                # and discard later if it stays empty.
                if current is not None and (
                    current.basis_contributions or current.atom_contributions
                ):
                    results.append(current)
                current = OrbitalBasisComposition(
                    method=method,
                    orbital_id=int(match[1]),
                    energy_au=float(match[2]),
                    occupation=float(match[3]),
                    orbital_type=match[4],
                )
                in_basis = True
                in_shell = False
                in_type = False
                in_atom = False
                continue

            if current is None:
                continue

            # Section markers
            if section_shell.search(line):
                in_basis = False
                in_shell = True
                in_type = False
                in_atom = False
                continue
            if section_type.search(line):
                in_basis = False
                in_shell = False
                in_type = True
                in_atom = False
                continue
            if section_atom.search(line):
                in_basis = False
                in_shell = False
                in_type = False
                in_atom = True
                continue

            # Basis contributions
            if in_basis and (match := basis_pat.match(line)):
                current.basis_contributions.append(
                    OrbitalBasisEntry(
                        basis_index=int(match[1]),
                        function_type=match[2].strip(),
                        atom_id=int(match[3]),
                        atom_element=match[4].strip(),
                        shell_index=int(match[5]),
                        local_pct=float(match[6]),
                        cross_pct=float(match[7]),
                        total_pct=float(match[8]),
                    )
                )
                continue

            # Shell compositions
            if in_shell and (match := shell_pat.search(line)):
                current.shell_contributions.append(
                    OrbitalShellEntry(
                        shell_index=int(match[1]),
                        shell_type=match[2].strip(),
                        atom_id=int(match[3]),
                        atom_element=match[4].strip(),
                        composition_pct=float(match[5]),
                    )
                )
                continue

            # Shell-type composition
            if in_type and (match := shell_type_pat.search(line)):
                current.shell_type_composition = OrbitalShellTypeComposition(
                    s=float(match[1]),
                    p=float(match[2]),
                    d=float(match[3]),
                    f=float(match[4]),
                    g=float(match[5]),
                    h=float(match[6]),
                )
                in_type = False
                continue

            # Atom compositions
            if in_atom and (match := atom_pat.search(line)):
                current.atom_contributions.append(
                    OrbitalAtomEntry(
                        atom_id=int(match[1]),
                        atom_element=match[2].strip(),
                        composition_pct=float(match[3]),
                    )
                )
                continue

            # Delocalization index (ends the block)
            if match := deloc_pat.search(line):
                if current is not None:
                    current.delocalization_index = float(match[1])
                    if (
                        current.basis_contributions
                        or current.atom_contributions
                    ):
                        results.append(current)
                    current = None
                    in_basis = False
                    in_shell = False
                    in_type = False
                    in_atom = False
                continue

        # Flush last orbital if not already flushed
        if current is not None and (
            current.basis_contributions or current.atom_contributions
        ):
            results.append(current)

        return results

    # ── Atom-only composition (Hirshfeld/Becke/fragment) ─────────────

    @staticmethod
    def parse_atom_compositions(
        stdout: str,
        method: str,
    ) -> list[OrbitalAtomComposition]:
        """Extract per-orbital atom-only compositions.

        These outputs have an orbital header, optionally a sum before
        normalization, then "Contributions after normalization:" with
        per-atom percentages, and a delocalization index.
        """
        results: list[OrbitalAtomComposition] = []

        # Orbital header (same format but may have different spacing)
        orb_header_pat = re.compile(
            rf"Orbital:\s+(\d+)\s+Energy\(a\.u\.\):\s+({FLOAT_PATTERN})"
            rf"\s+Occ:\s+({FLOAT_PATTERN})\s+Type:\s+(\S+)"
        )

        # "The sum of contributions before normalization   99.999199 %"
        sum_before_pat = re.compile(
            rf"sum of contributions before normalization\s+"
            rf"({FLOAT_PATTERN})\s*%",
            re.I,
        )

        # "Contributions after normalization:"
        contrib_header = re.compile(
            r"Contributions after normalization:", re.I
        )

        # "Atom     1(C ) :      4.109 %"
        atom_pat = re.compile(
            rf"Atom\s+(\d+)\s*\(([A-Za-z]+)\s*\)\s*:\s+"
            rf"({FLOAT_PATTERN})\s*%"
        )

        # Delocalization index
        deloc_pat = re.compile(
            rf"Orbital delocalization index:\s+({FLOAT_PATTERN})"
        )

        current: OrbitalAtomComposition | None = None
        in_contrib = False

        for line in stdout.split("\n"):
            # New orbital header
            if match := orb_header_pat.search(line):
                if current is not None and current.atom_contributions:
                    results.append(current)
                current = OrbitalAtomComposition(
                    method=method,
                    orbital_id=int(match[1]),
                    energy_au=float(match[2]),
                    occupation=float(match[3]),
                    orbital_type=match[4],
                )
                in_contrib = False
                continue

            if current is None:
                continue

            # Sum before normalization
            if match := sum_before_pat.search(line):
                current.sum_before_normalization = float(match[1])
                continue

            # Contributions header
            if contrib_header.search(line):
                in_contrib = True
                continue

            # Atom contributions
            if in_contrib and (match := atom_pat.search(line)):
                current.atom_contributions.append(
                    OrbitalAtomEntry(
                        atom_id=int(match[1]),
                        atom_element=match[2].strip(),
                        composition_pct=float(match[3]),
                    )
                )
                continue

            # Delocalization index (ends the block)
            if match := deloc_pat.search(line):
                if current is not None:
                    current.delocalization_index = float(match[1])
                    if current.atom_contributions:
                        results.append(current)
                    current = None
                    in_contrib = False
                continue

        # Flush last
        if current is not None and current.atom_contributions:
            results.append(current)

        return results

    # ── LOBA oxidation states ────────────────────────────────────────

    @staticmethod
    def parse_oxidation_states(stdout: str) -> list[OxidationState]:
        """Extract LOBA oxidation states."""
        pattern = (
            rf"Atom\s+(\d+)\s*\([A-Za-z]+\s*\).*?"
            rf"oxidation state.*?({FLOAT_PATTERN})"
        )
        return [
            OxidationState(
                atom_id=int(match[1]),
                oxidation_state=int(round(float(match[2]))),
            )
            for match in re.finditer(pattern, stdout, re.IGNORECASE)
        ]


# =============================================================================
# Menu 9: Bond order analysis
# =============================================================================


_BOND_ORDER_METHOD: dict[Menu, str] = {
    Menu.MAYER_BOND_ORDER: "mayer",
    Menu.WIBERG_BOND_ORDER: "wiberg",
    Menu.MULLIKEN_BOND_ORDER: "mulliken",
    Menu.FUZZY_BOND_ORDER: "fuzzy",
    Menu.LAPLACIAN_BOND_ORDER: "laplacian",
    Menu.MULTICENTER_BOND_ORDER: "multicenter",
    Menu.MULTICENTER_BOND_ORDER_NAO: "multicenter_nao",
    Menu.MULLIKEN_BOND_ORDER_DECOMPOSE: "mulliken_decomposition",
    Menu.WIBERG_DECOMPOSITION: "wiberg_decomposition",
    Menu.ORBITAL_PERTURBED_MAYER: "orbital_perturbed_mayer",
    Menu.IBSI_ANALYSIS: "ibsi",
}


class BondOrderParser(OutputParser):
    """Parser for bond orders (Menu 9).

    Each bond order block is stored as a ``BondOrderSet`` with a
    standardized ``method`` name so multiple methods coexist
    without overwriting.
    """

    @classmethod
    def parse_for_result(
        cls,
        analysis: Menu,
        stdout: str,
    ) -> list[ParsedMultiwfnResult]:
        method = _BOND_ORDER_METHOD.get(analysis, "unknown")

        case_map: dict[Menu, list] = {
            Menu.MULTICENTER_BOND_ORDER: [cls.parse_multicenter],
            Menu.MULTICENTER_BOND_ORDER_NAO: [cls.parse_multicenter],
            Menu.MULLIKEN_BOND_ORDER_DECOMPOSE: [cls.parse_decomposition],
            Menu.WIBERG_DECOMPOSITION: [cls.parse_decomposition],
            Menu.IBSI_ANALYSIS: [
                lambda s: cls.parse_ibsi(s),
                lambda s: cls.parse_bond_order_sets(s, method),
            ],
        }
        default = [
            lambda s: cls.parse_bond_order_sets(s, method),
            cls.parse_valence,
        ]
        parsers = case_map.get(analysis, default)
        return cls._collect(stdout, parsers)

    # ── Bond order sets ──────────────────────────────────────────────

    @staticmethod
    def parse_bond_order_sets(
        stdout: str,
        method: str,
    ) -> list[BondOrderSet]:
        """Extract all bond order blocks, each as a BondOrderSet.

        Handles both formats:
        - "#  N:  atom1(X)  atom2(X)  value"
        - "atom1(X)  atom2(X)  ..."  (IBSI-style, but the bond order
          part is handled here for the "total bond order" blocks)

        Multiple blocks (e.g. different thresholds, or "The total
        bond order" after individual ones) each get their own set.
        """
        sets: list[BondOrderSet] = []

        # Threshold header:
        # "Bond orders with absolute value >=  0.050000"
        # "The total bond order >=  0.050000"
        threshold_pat = re.compile(
            rf"(?:Bond orders|bond order|total bond order)"
            rf".*?>=?\s+({FLOAT_PATTERN})",
            re.I,
        )

        # Bond order line:
        # "#    1:         1(C )    2(C )    1.45268106"
        # or without leading #:
        # "#    1:         1(C )    2(C )    1.45268106"
        numbered_pat = re.compile(
            rf"#\s*(\d+):\s+(\d+)\s*\([^)]+\)\s+(\d+)\s*\([^)]+\)\s+"
            rf"({FLOAT_PATTERN})"
        )

        current_threshold: float | None = None
        current_bonds: list[BondOrder] = []
        current_label = method

        def _flush() -> None:
            nonlocal current_bonds, current_threshold, current_label
            if current_bonds:
                sets.append(
                    BondOrderSet(
                        method=current_label,
                        threshold=current_threshold,
                        bond_orders=list(current_bonds),
                    )
                )
            current_bonds = []
            current_threshold = None

        for line in stdout.split("\n"):
            # New threshold header starts a new set
            if match := threshold_pat.search(line):
                _flush()
                current_threshold = float(match[1])
                # Detect if this is "total bond order" vs regular
                if "total bond order" in line.lower():
                    current_label = f"{method}_total"
                else:
                    current_label = method
                continue

            # Numbered bond order line
            if match := numbered_pat.search(line):
                atom1 = int(match[2])
                atom2 = int(match[3])
                bo = float(match[4])
                if atom1 > atom2:
                    atom1, atom2 = atom2, atom1
                current_bonds.append(
                    BondOrder(
                        atom1_id=atom1,
                        atom2_id=atom2,
                        bond_order=bo,
                    )
                )
                continue

        _flush()
        return sets

    # ── Valences ─────────────────────────────────────────────────────

    @staticmethod
    def parse_valence(stdout: str) -> list[Valence]:
        """Extract total valence and free valence for each atom.

        Handles the two-column format:
        "Atom     1(C ) :    3.94042232    0.00000000"
        where first value is total valence and second is free valence.

        Also handles the single-value format:
        "Total valence of atom  1(C ): 3.940"
        "Free valence of atom  1(C ): 0.000"
        """
        valences: list[Valence] = []

        # Two-column format (Mayer style)
        two_col_pat = re.compile(
            rf"Atom\s+(\d+)\s*\([^)]+\)\s*:\s+"
            rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})"
        )
        in_valence_section = False
        for line in stdout.split("\n"):
            if "valences" in line.lower() and (
                "total" in line.lower() or "free" in line.lower()
            ):
                in_valence_section = True
                continue
            if in_valence_section:
                if match := two_col_pat.search(line):
                    atom_id = int(match[1])
                    valences.append(
                        Valence(
                            atom_id=atom_id,
                            type="total_valence",
                            valence=float(match[2]),
                        )
                    )
                    valences.append(
                        Valence(
                            atom_id=atom_id,
                            type="free_valence",
                            valence=float(match[3]),
                        )
                    )
                    continue
                # End of section if we hit a non-matching non-blank line
                if line.strip() and not two_col_pat.search(line):
                    in_valence_section = False

        if valences:
            return valences

        # Single-value format (fallback)
        total_pattern = (
            rf"Total valence of atom\s+(\d+)\s*\([^)]+\)\s*:\s+"
            rf"({FLOAT_PATTERN})"
        )
        free_pattern = (
            rf"Free valence of atom\s+(\d+)\s*\([^)]+\)\s*:\s+"
            rf"({FLOAT_PATTERN})"
        )

        valences.extend(
            Valence(
                atom_id=int(match[1]),
                type="total_valence",
                valence=float(match[2]),
            )
            for match in re.finditer(total_pattern, stdout)
        )
        valences.extend(
            Valence(
                atom_id=int(match[1]),
                type="free_valence",
                valence=float(match[2]),
            )
            for match in re.finditer(free_pattern, stdout)
        )
        return valences

    # ── IBSI analysis ────────────────────────────────────────────────

    @staticmethod
    def parse_ibsi(stdout: str) -> IBSIAnalysis | None:
        """Extract IBSI analysis entries.

        Format:
        "    1(C )    2(C )  Dist:  1.3950   Int(dg_pair): 0.87601   IBSI: 0.79442"
        """  # noqa: E501
        pattern = re.compile(
            rf"(\d+)\s*\([^)]+\)\s+(\d+)\s*\([^)]+\)\s+"
            rf"Dist:\s+({FLOAT_PATTERN})\s+"
            rf"Int\(dg_pair\):\s+({FLOAT_PATTERN})\s+"
            rf"IBSI:\s+({FLOAT_PATTERN})"
        )

        entries: list[IBSIEntry] = []
        for match in re.finditer(pattern, stdout):
            atom1 = int(match[1])
            atom2 = int(match[2])
            if atom1 > atom2:
                atom1, atom2 = atom2, atom1
            entries.append(
                IBSIEntry(
                    atom1_id=atom1,
                    atom2_id=atom2,
                    distance=float(match[3]),
                    int_dg_pair=float(match[4]),
                    ibsi=float(match[5]),
                )
            )

        if entries:
            return IBSIAnalysis(entries=entries)
        return None

    # ── Multicenter bond orders ──────────────────────────────────────

    @staticmethod
    def parse_multicenter(stdout: str) -> list[MultiCenterBondOrder]:
        """Extract multicenter bond order results."""
        results: list[MultiCenterBondOrder] = []
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

    # ── Per-orbital decomposition ────────────────────────────────────

    @staticmethod
    def parse_decomposition(stdout: str) -> list[BondOrderDecomposition]:
        """Extract per-orbital bond order decomposition."""
        pattern = rf"Orbital\s+(\d+)\s*:\s+({FLOAT_PATTERN})"
        return [
            BondOrderDecomposition(
                orbital_id=int(match[1]),
                contribution=float(match[2]),
            )
            for match in re.finditer(pattern, stdout)
        ]


# =============================================================================
# Menu 10: Density of states
# =============================================================================


class DOSParser(OutputParser):
    """Parser for density of states output."""

    @classmethod
    def parse_for_result(
        cls,
        analysis: Menu,
        stdout: str,
    ) -> list[ParsedMultiwfnResult]:
        results: list[ParsedMultiwfnResult] = []
        dos = cls.parse(stdout)
        if dos.energies_eV:
            results.append(dos)
        orb_energies = cls.parse_orbital_energies(stdout)
        results.extend(orb_energies)
        metadata = cls.parse_metadata(stdout)
        if metadata is not None:
            results.append(metadata)
        return results

    @staticmethod
    def parse(stdout: str) -> DensityOfStates:
        """Extract DOS curve data."""
        energies: list[float] = []
        dos_vals: list[float] = []
        pdos: dict[str, list[float]] = {}

        pattern2 = rf"^\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s*$"
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
                energies.append(energy)
                if vals:
                    dos_vals.append(vals[0])
                for k, v in enumerate(vals[1:], start=1):
                    key = f"pdos_{k}"
                    pdos.setdefault(key, []).append(v)
            elif match := re.match(pattern2, line):
                energies.append(float(match[1]))
                dos_vals.append(float(match[2]))

        return DensityOfStates(
            energies_eV=energies,
            dos=dos_vals,
            projected_dos=pdos if pdos else None,
        )

    @staticmethod
    def parse_orbital_energies(stdout: str) -> list[OrbitalEnergy]:
        """Extract orbital energies used in DOS."""
        pattern = rf"(\d+)\s+({FLOAT_PATTERN})\s+eV\s+Occ=\s*({FLOAT_PATTERN})"
        return [
            OrbitalEnergy(
                index=int(match[1]),
                energy_eV=float(match[2]),
                occupation=float(match[3]),
            )
            for match in re.finditer(pattern, stdout)
        ]

    @staticmethod
    def parse_metadata(stdout: str) -> DOSMetadata | None:
        """Extract TDOS center and HOMO level metadata."""
        tdos_center: float | None = None
        homo_level: float | None = None

        tdos_pat = re.compile(rf"Center of TDOS:\s+({FLOAT_PATTERN})\s*a\.u\.")
        homo_pat = re.compile(
            rf"vertical dash line corresponds to HOMO level at\s+"
            rf"({FLOAT_PATTERN})\s*a\.u\."
        )

        if match := tdos_pat.search(stdout):
            tdos_center = float(match[1])
        if match := homo_pat.search(stdout):
            homo_level = float(match[1])

        if tdos_center is not None or homo_level is not None:
            return DOSMetadata(
                tdos_center_au=tdos_center,
                homo_level_au=homo_level,
            )
        return None


# =============================================================================
# Menu 11: Spectra simulation
# =============================================================================


class SpectrumParser(OutputParser):
    """Parser for spectrum data."""

    @classmethod
    def _parse_spectrum_or_none(cls, stdout: str) -> Spectrum | None:
        """Return parsed spectrum only if it contains data."""
        s = cls.parse(stdout)
        return (
            s if (s.frequencies or s.wavelengths or s.atom_indices) else None
        )

    @classmethod
    def parse_for_result(
        cls,
        analysis: Menu,
        stdout: str,
    ) -> list[ParsedMultiwfnResult]:
        base: list = [
            cls._parse_spectrum_or_none,
            cls.parse_spectrum_curve_extrema,
            cls.parse_transitions,
        ]
        case_map: dict[Menu, list] = {
            Menu.PREDICT_COLOR: [*base, cls.parse_color],
        }
        parsers = case_map.get(analysis, base)
        return cls._collect(stdout, parsers)

    @staticmethod
    def parse(stdout: str) -> Spectrum:
        """Extract spectrum data (frequencies/wavelengths, intensities)."""
        spectrum: dict[str, list[float]] = {
            "frequencies": [],
            "intensities": [],
        }

        pattern1 = (
            rf"({FLOAT_PATTERN})\s+cm\^?-1.*?Intensity:\s+({FLOAT_PATTERN})"
        )
        pattern_uv = (
            rf"({FLOAT_PATTERN})\s+nm.*?(?:f=|Str[.=])\s*({FLOAT_PATTERN})"
        )
        pattern_nmr = rf"Atom\s+(\d+)\s*\([^)]+\)\s+shift:\s+({FLOAT_PATTERN})"
        pattern2 = rf"^\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s*$"

        for match in re.finditer(pattern1, stdout):
            spectrum["frequencies"].append(float(match[1]))
            spectrum["intensities"].append(float(match[2]))

        if spectrum["frequencies"]:
            return Spectrum(
                frequencies=spectrum["frequencies"],
                intensities=spectrum["intensities"],
            )

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
        """Extract discrete transition data."""
        transitions: list[Transition] = []
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
    def parse_spectrum_curve_extrema(
        stdout: str,
    ) -> SpectrumCurveExtrema | None:
        """Extract extrema reported under spectrum-curve extrema blocks."""
        if "Extrema on the spectrum curve" not in stdout:
            return None

        pattern = (
            rf"(Maximum|Minimum)\s+(\d+)\s+X:\s*({FLOAT_PATTERN})\s+"
            rf"Value:\s*({FLOAT_PATTERN})\.?"
        )

        extrema: list[SpectrumExtremum] = []
        for match in re.finditer(pattern, stdout, flags=re.IGNORECASE):
            kind = match[1].lower()
            extrema.append(
                SpectrumExtremum(
                    kind="maximum" if kind == "maximum" else "minimum",
                    index=int(match[2]),
                    x=float(match[3]),
                    value=float(match[4]),
                )
            )

        if not extrema:
            return None
        return SpectrumCurveExtrema(extrema=extrema)

    @staticmethod
    def parse_color(stdout: str) -> Color | None:
        """Extract predicted colour from PREDICT_COLOR output."""
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
        return Color(**result) if result else None


# =============================================================================
# Menu 12: Quantitative molecular surface analysis
# =============================================================================


class SurfaceParser(OutputParser):
    """Parser for molecular surface analysis output (Menu 12)."""

    @classmethod
    def parse_for_result(
        cls,
        analysis: Menu,
        stdout: str,
    ) -> list[ParsedMultiwfnResult]:
        mapped_property = _SURFACE_PROPERTY.get(analysis, "unknown")
        result = cls.parse(stdout, mapped_property)
        return [result] if result is not None else []

    @classmethod
    def parse(
        cls,
        stdout: str,
        mapped_property: str,
    ) -> SurfaceAnalysisResult | None:
        """Parse full surface analysis output into one result."""
        geometry = cls.parse_geometry(stdout)
        minima = cls.parse_extrema(stdout, "min")
        maxima = cls.parse_extrema(stdout, "max")
        statistics = cls.parse_statistics(stdout, mapped_property)

        # Only return if we found something meaningful
        if (
            geometry is None
            and not minima
            and not maxima
            and statistics is None
        ):
            return None

        return SurfaceAnalysisResult(
            mapped_property=mapped_property,
            geometry=geometry,
            minima=minima,
            maxima=maxima,
            statistics=statistics,
        )

    # ── Geometry ─────────────────────────────────────────────────────

    @staticmethod
    def parse_geometry(stdout: str) -> SurfaceGeometry | None:
        """Extract isosurface geometry data."""
        geo = SurfaceGeometry()
        found = False

        # Volume
        vol_pat = re.compile(
            rf"Volume.*?:\s+({FLOAT_PATTERN})\s+Bohr\^3\s+"
            rf"\(\s*({FLOAT_PATTERN})\s+Angstrom\^3\)",
            re.I,
        )
        if match := vol_pat.search(stdout):
            geo.volume_bohr3 = float(match[1])
            geo.volume_angstrom3 = float(match[2])
            found = True

        # Area
        area_pat = re.compile(
            rf"(?:Isosurface|Overall surface)\s+area:\s+"
            rf"({FLOAT_PATTERN})\s+Bohr\^2\s+"
            rf"\(\s*({FLOAT_PATTERN})\s+Angstrom\^2\)",
            re.I,
        )
        if match := area_pat.search(stdout):
            geo.area_bohr2 = float(match[1])
            geo.area_angstrom2 = float(match[2])
            found = True

        # Sphericity
        sph_pat = re.compile(rf"Sphericity:\s+({FLOAT_PATTERN})")
        if match := sph_pat.search(stdout):
            geo.sphericity = float(match[1])
            found = True

        # Density
        dens_pat = re.compile(
            rf"Estimated density.*?:\s+({FLOAT_PATTERN})\s+g/cm\^3",
            re.I,
        )
        if match := dens_pat.search(stdout):
            geo.density_g_cm3 = float(match[1])
            found = True

        # Vertices/edges/facets (after elimination)
        vef_pat = re.compile(
            r"After elimination.*?V=\s*(\d+).*?E=\s*(\d+).*?"
            r"F=\s*(\d+)"
        )
        if match := vef_pat.search(stdout):
            geo.n_vertices = int(match[1])
            geo.n_edges = int(match[2])
            geo.n_facets = int(match[3])
            found = True

        return geo if found else None

    # ── Extrema ──────────────────────────────────────────────────────

    @staticmethod
    def parse_extrema(
        stdout: str,
        ext_type: Literal["min", "max"],
    ) -> list[SurfaceExtremum]:
        """Extract surface minima or maxima table."""
        extrema: list[SurfaceExtremum] = []

        # Section header
        if ext_type == "min":
            header_pat = re.compile(r"The number of surface minima:\s+(\d+)")
        else:
            header_pat = re.compile(r"The number of surface maxima:\s+(\d+)")

        # Entry line (with optional * for global):
        # "*    1 -0.02750965   -0.748576  -17.262580   -0.225  0.366  -1.831"
        # "     1 -0.02750962   -0.748575  -17.262563   -0.431  0.028  -1.831"
        entry_pat = re.compile(
            rf"^\s*(\*?)\s*(\d+)\s+({FLOAT_PATTERN})\s+"
            rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+"
            rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+"
            rf"({FLOAT_PATTERN})"
        )

        in_section = False
        for line in stdout.split("\n"):
            if header_pat.search(line):
                in_section = True
                continue
            if in_section:
                # End of section: blank line or next section header
                if (
                    line.strip() == ""
                    or "The number of surface" in line
                    or "Summary of surface" in line
                ):
                    if extrema:
                        break
                    continue
                # Skip the column header line
                if "#" in line and "a.u." in line:
                    continue
                if match := entry_pat.match(line):
                    extrema.append(
                        SurfaceExtremum(
                            type=ext_type,
                            index=int(match[2]),
                            value_au=float(match[3]),
                            value_eV=float(match[4]),
                            value_kcal_mol=float(match[5]),
                            x_angstrom=float(match[6]),
                            y_angstrom=float(match[7]),
                            z_angstrom=float(match[8]),
                            is_global=match[1] == "*",
                        )
                    )

        return extrema

    # ── Statistics ───────────────────────────────────────────────────

    @staticmethod
    def parse_statistics(
        stdout: str,
        mapped_property: str,
    ) -> SurfaceStatistics | None:
        """Extract the Summary of surface analysis block."""
        # Find the summary section
        summary_start = stdout.find("Summary of surface analysis")
        if summary_start == -1:
            return None

        block = stdout[summary_start:]
        stats = SurfaceStatistics(mapped_property=mapped_property)
        found = False

        # Helper to search within block
        def _get(pattern: str, flags: int = re.I) -> re.Match[str] | None:
            return re.search(pattern, block, flags)

        # ── Global extrema from summary ──
        # "Minimal value:    -17.26258 kcal/mol   Maximal value:     11.82103 kcal/mol"  # noqa: E501
        minmax_pat = (
            rf"Minimal value:\s+({FLOAT_PATTERN})\s+kcal/mol\s+"
            rf"Maximal value:\s+({FLOAT_PATTERN})\s+kcal/mol"
        )
        if match := _get(minmax_pat):
            stats.global_min_kcal_mol = float(match[1])
            stats.global_max_kcal_mol = float(match[2])
            found = True

        # Global min/max in a.u. from the pre-summary lines
        gmin_pat = rf"Global surface minimum:\s+({FLOAT_PATTERN})\s+a\.u\."
        gmax_pat = rf"Global surface maximum:\s+({FLOAT_PATTERN})\s+a\.u\."
        if match := re.search(gmin_pat, stdout):
            stats.global_min_au = float(match[1])
            found = True
        if match := re.search(gmax_pat, stdout):
            stats.global_max_au = float(match[1])
            found = True

        # ── Areas ──
        # "Overall surface area:  448.293 Bohr^2  ( 125.535 Angstrom^2)"
        area_both = (
            rf"({FLOAT_PATTERN})\s+Bohr\^2\s+"
            rf"\(\s*({FLOAT_PATTERN})\s+Angstrom\^2\)"
        )
        if match := _get(r"Overall surface area:\s+" + area_both):
            stats.overall_area_bohr2 = float(match[1])
            stats.overall_area_angstrom2 = float(match[2])
            found = True
        if match := _get(r"Positive surface area:\s+" + area_both):
            stats.positive_area_bohr2 = float(match[1])
            stats.positive_area_angstrom2 = float(match[2])
            found = True
        if match := _get(r"Negative surface area:\s+" + area_both):
            stats.negative_area_bohr2 = float(match[1])
            stats.negative_area_angstrom2 = float(match[2])
            found = True

        # ── Averages ──
        # "Overall average value:  0.00002044 a.u. (  0.01282 kcal/mol)"
        avg_both = (
            rf"({FLOAT_PATTERN})\s+a\.u\.\s+"
            rf"\(\s*({FLOAT_PATTERN})\s+kcal/mol\)"
        )
        if match := _get(r"Overall average value:\s+" + avg_both):
            stats.overall_average_au = float(match[1])
            stats.overall_average_kcal_mol = float(match[2])
            found = True
        if match := _get(r"Positive average value:\s+" + avg_both):
            stats.positive_average_au = float(match[1])
            stats.positive_average_kcal_mol = float(match[2])
            found = True
        if match := _get(r"Negative average value:\s+" + avg_both):
            stats.negative_average_au = float(match[1])
            stats.negative_average_kcal_mol = float(match[2])
            found = True

        # ── Variances ──
        # "Overall variance (sigma^2_tot):  0.00009608 a.u.^2 (  37.83 (kcal/mol)^2)"  # noqa: E501
        var_both = (
            rf"({FLOAT_PATTERN})\s+a\.u\.\^2\s+"
            rf"\(\s*({FLOAT_PATTERN})\s+\(kcal/mol\)\^2\)"
        )
        if match := _get(r"Overall variance.*?sigma.*?:\s+" + var_both):
            stats.sigma2_total_au2 = float(match[1])
            stats.sigma2_total_kcal_mol2 = float(match[2])
            found = True
        if match := _get(r"Positive variance:\s+" + var_both):
            stats.positive_variance_au2 = float(match[1])
            stats.positive_variance_kcal_mol2 = float(match[2])
            found = True
        if match := _get(r"Negative variance:\s+" + var_both):
            stats.negative_variance_au2 = float(match[1])
            stats.negative_variance_kcal_mol2 = float(match[2])
            found = True

        # ── Nu ──
        # "Balance of charges (nu):   0.19601755"
        if match := _get(rf"Balance of charges.*?:\s+({FLOAT_PATTERN})"):
            stats.nu = float(match[1])
            found = True

        # ── sigma^2_tot * nu ──
        # "Product of sigma^2_tot and nu:  0.00001883 a.u.^2 (  7.416 (kcal/mol)^2)"  # noqa: E501
        if match := _get(r"Product of sigma.*?nu:\s+" + var_both):
            stats.sigma2_tot_times_nu_au2 = float(match[1])
            stats.sigma2_tot_times_nu_kcal_mol2 = float(match[2])
            found = True

        # ── Pi ──
        # "Internal charge separation (Pi):  0.01236 a.u. (  7.757 kcal/mol)"
        if match := _get(r"Internal charge separation.*?:\s+" + avg_both):
            stats.pi_au = float(match[1])
            stats.pi_kcal_mol = float(match[2])
            found = True

        # ── MPI ──
        # "Molecular polarity index (MPI):  0.33644 eV (  7.758 kcal/mol)"
        mpi_pat = (
            rf"Molecular polarity index.*?:\s+"
            rf"({FLOAT_PATTERN})\s+eV\s+"
            rf"\(\s*({FLOAT_PATTERN})\s+kcal/mol\)"
        )
        if match := _get(mpi_pat):
            stats.mpi_eV = float(match[1])
            stats.mpi_kcal_mol = float(match[2])
            found = True

        # ── Polar/nonpolar area ──
        # "Nonpolar surface area (|ESP| <= 10 kcal/mol):  88.18 Angstrom^2  ( 70.25 %)"  # noqa: E501
        nonpolar_pat = (
            rf"Nonpolar surface area.*?:\s+"
            rf"({FLOAT_PATTERN})\s+Angstrom\^2\s+"
            rf"\(\s*({FLOAT_PATTERN})\s*%\s*\)"
        )
        if match := _get(nonpolar_pat):
            stats.nonpolar_area_angstrom2 = float(match[1])
            stats.nonpolar_area_pct = float(match[2])
            found = True

        polar_pat = (
            rf"Polar surface area.*?:\s+"
            rf"({FLOAT_PATTERN})\s+Angstrom\^2\s+"
            rf"\(\s*({FLOAT_PATTERN})\s*%\s*\)"
        )
        if match := _get(polar_pat):
            stats.polar_area_angstrom2 = float(match[1])
            stats.polar_area_pct = float(match[2])
            found = True

        # ── Skewness ──
        if match := _get(rf"Overall skewness:\s+({FLOAT_PATTERN})"):
            stats.overall_skewness = float(match[1])
            found = True
        if match := _get(rf"Positive skewness:\s+({FLOAT_PATTERN})"):
            stats.positive_skewness = float(match[1])
            found = True
        if match := _get(rf"Negative skewness:\s+({FLOAT_PATTERN})"):
            stats.negative_skewness = float(match[1])
            found = True

        return stats if found else None


# =============================================================================
# Menu 15: Fuzzy atomic space analysis
# =============================================================================

# Legacy placeholder kept for backwards compatibility with older mappings.
_FUZZY_PROPERTY: dict[Menu, str] = {}

# For the specific sub-menus, we'll detect from Menu name
_FUZZY_PROPERTY_FROM_NAME: dict[str, str] = {
    "FUZZY_INTEGRATE_EDENSITY": "edensity",
    "FUZZY_INTEGRATE_NORM_RHO": "norm_rho",
    "FUZZY_INTEGRATE_LAPLACIAN": "laplacian",
    "FUZZY_INTEGRATE_ORB_WFN_HOMO": "orb_wfn_homo",
    "FUZZY_INTEGRATE_ORB_WFN_LUMO": "orb_wfn_lumo",
    "FUZZY_INTEGRATE_ESPIN_DENSITY": "espin_density",
    "FUZZY_INTEGRATE_KR": "kr",
    "FUZZY_INTEGRATE_GR": "gr",
    "FUZZY_INTEGRATE_ESP_CHARGES": "esp_charges",
    "FUZZY_INTEGRATE_ELF": "elf",
    "FUZZY_INTEGRATE_LOL": "lol",
    "FUZZY_INTEGRATE_LOCAL_ENTROPY": "local_entropy",
    "FUZZY_INTEGRATE_ESP": "esp",
    "FUZZY_INTEGRATE_RDG": "rdg",
    "FUZZY_INTEGRATE_RDG_PROMOLECULAR": "rdg_promolecular",
    "FUZZY_INTEGRATE_LAMBDA2RHO": "lambda2rho",
    "FUZZY_INTEGRATE_LAMBDA2RGO_PROMOLECULAR": "lambda2rho_promolecular",
    "FUZZY_INTEGRATE_ALIE": "alie",
    "FUZZY_INTEGRATE_EDR": "edr",
    "FUZZY_INTEGRATE_ORB_OVERLAP_DR": "orb_overlap_dr",
    "FUZZY_INTEGRATE_DELTAG_PROMOLECULAR": "deltag_promolecular",
    "FUZZY_INTEGRATE_DELTAG_HIRSHFELD": "deltag_hirshfeld",
    "FUZZY_INTEGRATE_IRI": "iri",
    "FUZZY_OVERLAP_EDENSITY": "edensity",
    "FUZZY_OVERLAP_NORM_RHO": "norm_rho",
    "FUZZY_OVERLAP_LAPLACIAN": "laplacian",
    "FUZZY_OVERLAP_ORB_WFN_HOMO": "orb_wfn_homo",
    "FUZZY_OVERLAP_ORB_WFN_LUMO": "orb_wfn_lumo",
    "FUZZY_OVERLAP_ESPIN_DENSITY": "espin_density",
    "FUZZY_OVERLAP_KR": "kr",
    "FUZZY_OVERLAP_GR": "gr",
    "FUZZY_OVERLAP_ESP_CHARGES": "esp_charges",
    "FUZZY_OVERLAP_ELF": "elf",
    "FUZZY_OVERLAP_LOL": "lol",
    "FUZZY_OVERLAP_LOCAL_ENTROPY": "local_entropy",
    "FUZZY_OVERLAP_ESP": "esp",
    "FUZZY_OVERLAP_RDG": "rdg",
    "FUZZY_OVERLAP_RDG_PROMOLECULAR": "rdg_promolecular",
    "FUZZY_OVERLAP_LAMBDA2RHO": "lambda2rho",
    "FUZZY_OVERLAP_LAMBDA2RGO_PROMOLECULAR": "lambda2rho_promolecular",
    "FUZZY_OVERLAP_ALIE": "alie",
    "FUZZY_OVERLAP_EDR": "edr",
    "FUZZY_OVERLAP_ORB_OVERLAP_DR": "orb_overlap_dr",
    "FUZZY_OVERLAP_DELTAG_PROMOLECULAR": "deltag_promolecular",
    "FUZZY_OVERLAP_DELTAG_HIRSHFELD": "deltag_hirshfeld",
    "FUZZY_OVERLAP_IRI": "iri",
}


class FuzzySpaceParser(OutputParser):
    """Parser for fuzzy atomic space analysis output (Menu 15)."""

    @classmethod
    def parse_for_result(
        cls,
        analysis: Menu,
        stdout: str,
    ) -> list[ParsedMultiwfnResult]:
        results: list[ParsedMultiwfnResult] = []

        # Determine integrated property from Menu name
        prop = _FUZZY_PROPERTY_FROM_NAME.get(analysis.name, "unknown")

        # Fuzzy integration table
        integ = cls.parse_fuzzy_integration(stdout, prop)
        if integ is not None:
            results.append(integ)

        # Atomic multipoles
        results.extend(cls.parse_atomic_multipoles(stdout))

        # Molecular multipoles
        mol = cls.parse_molecular_multipole(stdout)
        if mol is not None:
            results.append(mol)

        # AOM diagnostics
        aom = cls.parse_aom_diagnostics(stdout)
        if aom is not None:
            results.append(aom)

        # Delocalization index matrices
        results.extend(cls.parse_di_matrices(stdout))

        # Overlap integration matrices
        results.extend(cls.parse_overlap_matrices(stdout, prop))

        # CLRK matrix
        clrk = cls.parse_clrk_matrix(stdout)
        if clrk is not None:
            results.append(clrk)

        # FLU reference parameters
        results.extend(cls.parse_flu_references(stdout))

        # Aromaticity indices
        results.extend(cls.parse_aromaticity_index(stdout))

        # Pairwise delocalization indices
        results.extend(cls.parse_delocalization_indices(stdout))

        return results

    # ── Fuzzy integration table ──────────────────────────────────────

    @staticmethod
    def parse_fuzzy_integration(
        stdout: str,
        integrated_property: str,
    ) -> FuzzyIntegrationResult | None:
        """Extract per-atom fuzzy integration values."""
        # "     1(C )            6.21199090            14.790461            14.790461"  # noqa: E501
        entry_pat = re.compile(
            rf"^\s+(\d+)\s*\(([A-Za-z]+)\s*\)\s+"
            rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+"
            rf"({FLOAT_PATTERN})"
        )
        sum_pat = re.compile(rf"Summing up above values:\s+({FLOAT_PATTERN})")
        sum_abs_pat = re.compile(
            rf"Summing up absolute value.*?:\s+({FLOAT_PATTERN})"
        )

        entries: list[FuzzyIntegrationEntry] = []
        total_sum: float | None = None
        total_sum_abs: float | None = None

        in_table = False
        for line in stdout.split("\n"):
            if "Atomic space" in line and "Value" in line:
                in_table = True
                continue
            if in_table:
                if match := entry_pat.match(line):
                    entries.append(
                        FuzzyIntegrationEntry(
                            atom_id=int(match[1]),
                            atom_element=match[2].strip(),
                            value=float(match[3]),
                            pct_of_sum=float(match[4]),
                            pct_of_sum_abs=float(match[5]),
                        )
                    )
                    continue
                if match := sum_pat.search(line):
                    total_sum = float(match[1])
                    continue
                if match := sum_abs_pat.search(line):
                    total_sum_abs = float(match[1])
                    in_table = False
                    continue

        if entries:
            return FuzzyIntegrationResult(
                integrated_property=integrated_property,
                entries=entries,
                total_sum=total_sum,
                total_sum_abs=total_sum_abs,
            )
        return None

    # ── Atomic multipoles ────────────────────────────────────────────

    @staticmethod
    def parse_atomic_multipoles(
        stdout: str,
    ) -> list[AtomicMultipole]:
        """Extract per-atom multipole moments."""
        results: list[AtomicMultipole] = []

        # Section header: "*****  Atom     1(C )  *****"
        atom_header = re.compile(
            r"\*+\s*Atom\s+(\d+)\s*\(([A-Za-z]+)\s*\)\s*\*+"
        )

        charge_pat = re.compile(rf"Atomic charge:\s+({FLOAT_PATTERN})")
        monopole_pat = re.compile(
            rf"Atomic monopole moment.*?:\s+({FLOAT_PATTERN})"
        )
        dipole_pat = re.compile(
            rf"X=\s+({FLOAT_PATTERN})\s+Y=\s+({FLOAT_PATTERN})\s+"
            rf"Z=\s+({FLOAT_PATTERN})\s+Norm=\s+({FLOAT_PATTERN})"
        )
        quad_cart_pat = re.compile(rf"([XY][XYZ])=\s+({FLOAT_PATTERN})")
        quad_mag_pat = re.compile(
            rf"Magnitude of the traceless quadrupole.*?:\s+"
            rf"({FLOAT_PATTERN})"
        )
        extent_pat = re.compile(
            rf"Atomic electronic spatial extent.*?:\s+"
            rf"({FLOAT_PATTERN})"
        )
        extent_comp_pat = re.compile(
            rf"Components of <r\^2>:\s+X=\s+({FLOAT_PATTERN})\s+"
            rf"Y=\s+({FLOAT_PATTERN})\s+Z=\s+({FLOAT_PATTERN})"
        )
        octo_mag_pat = re.compile(
            rf"Magnitude:\s+\|Q_3\|=\s+({FLOAT_PATTERN})"
        )

        current: AtomicMultipole | None = None
        in_dipole = False
        in_contrib = False
        in_quad_traceless = False

        for line in stdout.split("\n"):
            if match := atom_header.search(line):
                if current is not None:
                    results.append(current)
                current = AtomicMultipole(
                    atom_id=int(match[1]),
                    atom_element=match[2].strip(),
                )
                in_dipole = False
                in_contrib = False
                in_quad_traceless = False
                continue

            if current is None:
                continue

            if "Molecular dipole" in line or "Molecular quadrupole" in line:
                if current is not None:
                    results.append(current)
                    current = None
                continue

            if match := charge_pat.search(line):
                current.atomic_charge = float(match[1])
                continue

            if match := monopole_pat.search(line):
                current.monopole_moment = float(match[1])
                continue

            if "Atomic dipole moments:" in line:
                in_dipole = True
                in_contrib = False
                continue
            if "Contribution to molecular dipole" in line:
                in_dipole = False
                in_contrib = True
                continue

            if (in_dipole or in_contrib) and (
                match := dipole_pat.search(line)
            ):
                if in_dipole:
                    current.dipole_x = float(match[1])
                    current.dipole_y = float(match[2])
                    current.dipole_z = float(match[3])
                    current.dipole_norm = float(match[4])
                    in_dipole = False
                elif in_contrib:
                    current.mol_dipole_contrib_x = float(match[1])
                    current.mol_dipole_contrib_y = float(match[2])
                    current.mol_dipole_contrib_z = float(match[3])
                    current.mol_dipole_contrib_norm = float(match[4])
                    in_contrib = False
                continue

            if "Traceless Cartesian form" in line:
                in_quad_traceless = True
                continue
            if in_quad_traceless:
                for qm in quad_cart_pat.finditer(line):
                    comp = qm[1]
                    val = float(qm[2])
                    if comp == "XX":
                        current.quadrupole_xx = val
                    elif comp == "XY":
                        current.quadrupole_xy = val
                    elif comp == "XZ":
                        current.quadrupole_xz = val
                    elif comp == "YY":
                        current.quadrupole_yy = val
                    elif comp == "YZ":
                        current.quadrupole_yz = val
                    elif comp == "ZZ":
                        current.quadrupole_zz = val
                        in_quad_traceless = False

            if match := quad_mag_pat.search(line):
                current.quadrupole_magnitude = float(match[1])
                continue

            if match := extent_pat.search(line):
                current.spatial_extent_r2 = float(match[1])
                continue
            if match := extent_comp_pat.search(line):
                current.spatial_extent_x = float(match[1])
                current.spatial_extent_y = float(match[2])
                current.spatial_extent_z = float(match[3])
                continue

            if match := octo_mag_pat.search(line):
                current.octopole_magnitude = float(match[1])
                continue

        if current is not None:
            results.append(current)

        return results

    # ── Molecular multipoles ─────────────────────────────────────────

    @staticmethod
    def parse_molecular_multipole(
        stdout: str,
    ) -> MolecularMultipole | None:
        """Extract molecular multipole moments."""
        if "Molecular dipole" not in stdout:
            return None

        mol = MolecularMultipole()
        found = False

        electrons_pat = re.compile(
            rf"Total number of electrons:\s+({FLOAT_PATTERN})\s+"
            rf"Net charge:\s+({FLOAT_PATTERN})"
        )
        if match := electrons_pat.search(stdout):
            mol.total_electrons = float(match[1])
            mol.net_charge = float(match[2])
            found = True

        dip_au_pat = re.compile(
            rf"Molecular dipole moment \(a\.u\.\):\s+"
            rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+"
            rf"({FLOAT_PATTERN})"
        )
        if match := dip_au_pat.search(stdout):
            mol.dipole_x_au = float(match[1])
            mol.dipole_y_au = float(match[2])
            mol.dipole_z_au = float(match[3])
            found = True

        dip_debye_pat = re.compile(
            rf"Molecular dipole moment \(Debye\):\s+"
            rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+"
            rf"({FLOAT_PATTERN})"
        )
        if match := dip_debye_pat.search(stdout):
            mol.dipole_x_debye = float(match[1])
            mol.dipole_y_debye = float(match[2])
            mol.dipole_z_debye = float(match[3])
            found = True

        dip_mag_pat = re.compile(
            rf"Magnitude of molecular dipole.*?a\.u\.&Debye\):\s+"
            rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})"
        )
        if match := dip_mag_pat.search(stdout):
            mol.dipole_magnitude_au = float(match[1])
            mol.dipole_magnitude_debye = float(match[2])
            found = True

        # Traceless quadrupole
        traceless_start = stdout.find(
            "Molecular quadrupole moments (Traceless"
        )
        if traceless_start != -1:
            block = stdout[traceless_start : traceless_start + 500]
            comp_pat = re.compile(rf"([XY][XYZ])=\s+({FLOAT_PATTERN})")
            for qm in comp_pat.finditer(block):
                comp = qm[1]
                val = float(qm[2])
                if comp == "XX":
                    mol.quadrupole_xx = val
                elif comp == "XY":
                    mol.quadrupole_xy = val
                elif comp == "XZ":
                    mol.quadrupole_xz = val
                elif comp == "YY":
                    mol.quadrupole_yy = val
                elif comp == "YZ":
                    mol.quadrupole_yz = val
                elif comp == "ZZ":
                    mol.quadrupole_zz = val
            found = True

        quad_mag_pat = re.compile(
            rf"Magnitude of the traceless quadrupole moment tensor:\s+"
            rf"({FLOAT_PATTERN})"
        )
        # Find the one after "Molecular"
        mol_section = stdout[stdout.find("Molecular dipole") :]
        if match := quad_mag_pat.search(mol_section):
            mol.quadrupole_magnitude = float(match[1])

        mol_extent_pat = re.compile(
            rf"Molecular electronic spatial extent.*?:\s+"
            rf"({FLOAT_PATTERN})"
        )
        if match := mol_extent_pat.search(stdout):
            mol.spatial_extent_r2 = float(match[1])
            found = True

        mol_extent_comp = re.compile(
            rf"Components of <r\^2>:\s+X=\s+({FLOAT_PATTERN})\s+"
            rf"Y=\s+({FLOAT_PATTERN})\s+Z=\s+({FLOAT_PATTERN})"
        )
        # Take the last one (molecular, not atomic)
        matches = list(mol_extent_comp.finditer(stdout))
        if matches:
            m = matches[-1]
            mol.spatial_extent_x = float(m[1])
            mol.spatial_extent_y = float(m[2])
            mol.spatial_extent_z = float(m[3])

        octo_mag = re.compile(rf"Magnitude:\s+\|Q_3\|=\s+({FLOAT_PATTERN})")
        mol_octo_matches = list(octo_mag.finditer(mol_section))
        if mol_octo_matches:
            mol.octopole_magnitude = float(mol_octo_matches[-1][1])

        return mol if found else None

    # ── AOM diagnostics ──────────────────────────────────────────────

    @staticmethod
    def parse_aom_diagnostics(stdout: str) -> AOMDiagnostics | None:
        """Extract AOM quality diagnostics."""
        error_pat = re.compile(rf"Error of AOM is\s+({FLOAT_PATTERN})")
        diag_pat = re.compile(
            rf"Maximum diagonal deviation to 1:\s+({FLOAT_PATTERN})"
            rf"\s+at orbital\s+(\d+)"
        )
        nondiag_pat = re.compile(
            rf"Maximum nondiagonal deviation to 0:\s+"
            rf"({FLOAT_PATTERN})\s+between orbitals\s+(\d+)\s+(\d+)"
        )
        export_pat = re.compile(r"exported to\s+(\S+)\s+in current folder")

        aom = AOMDiagnostics()
        found = False

        if match := error_pat.search(stdout):
            aom.error = float(match[1])
            found = True
        if match := diag_pat.search(stdout):
            aom.max_diagonal_deviation = float(match[1])
            aom.max_diagonal_orbital = int(match[2])
            found = True
        if match := nondiag_pat.search(stdout):
            aom.max_nondiagonal_deviation = float(match[1])
            aom.max_nondiagonal_orbitals = (
                int(match[2]),
                int(match[3]),
            )
            found = True
        # Find AOM-specific export
        for m in export_pat.finditer(stdout):
            if "AOM" in m[1]:
                aom.exported_file = m[1]
                found = True

        return aom if found else None

    # ── Delocalization index matrices ────────────────────────────────

    @classmethod
    def parse_di_matrices(
        cls,
        stdout: str,
    ) -> list[DelocalizationIndexMatrix]:
        """Extract delocalization/localization index matrices."""
        results: list[DelocalizationIndexMatrix] = []

        header_pat = re.compile(
            r"\*+\s*(.*?(?:delocalization|localization)\s+index"
            r"\s+matrix)\s*\*+",
            re.I,
        )
        col_header_pat = re.compile(r"^\s+(\d+(?:\s+\d+)*)\s*$")
        row_pat = re.compile(rf"^\s+(\d+)((?:\s+{FLOAT_PATTERN})+)\s*$")
        li_pat = re.compile(
            rf"(\d+)\s*\(([A-Za-z]+)\s*\)\s*:\s+({FLOAT_PATTERN})"
        )

        in_matrix = False
        label = ""
        current_cols: list[int] = []
        data: dict[tuple[int, int], float] = {}
        max_idx = 0
        loc_indices: list[LocalizationIndex] = []

        def _flush() -> None:
            nonlocal data, max_idx, label, in_matrix, loc_indices
            if data:
                results.append(
                    DelocalizationIndexMatrix(
                        label=label,
                        n_atoms=max_idx,
                        data=dict(data),
                        localization_indices=list(loc_indices),
                    )
                )
            data = {}
            max_idx = 0
            label = ""
            in_matrix = False
            loc_indices = []

        for line in stdout.split("\n"):
            if match := header_pat.search(line):
                if in_matrix:
                    _flush()
                in_matrix = True
                label = match[1].strip()
                current_cols = []
                continue

            if not in_matrix:
                # Localization index lines outside matrix
                if "Localization index:" in line:
                    for m in li_pat.finditer(line):
                        loc_indices.append(
                            LocalizationIndex(
                                atom_id=int(m[1]),
                                atom_element=m[2].strip(),
                                index=float(m[3]),
                            )
                        )
                    # Also check subsequent lines
                    continue
                if loc_indices and (
                    match := re.search(
                        rf"(\d+)\s*\([A-Za-z]+\s*\)\s*:\s+"
                        rf"({FLOAT_PATTERN})",
                        line,
                    )
                ):
                    for m in li_pat.finditer(line):
                        loc_indices.append(
                            LocalizationIndex(
                                atom_id=int(m[1]),
                                atom_element=m[2].strip(),
                                index=float(m[3]),
                            )
                        )
                    continue
                continue

            # End markers
            if (
                "Localization index:" in line
                or "Note:" in line
                or "Current FLU" in line
                or "Integration of" in line
                or "Condensed linear" in line
            ):
                _flush()
                # Parse LI from this line if present
                if "Localization index:" in line:
                    for m in li_pat.finditer(line):
                        loc_indices.append(
                            LocalizationIndex(
                                atom_id=int(m[1]),
                                atom_element=m[2].strip(),
                                index=float(m[3]),
                            )
                        )
                continue

            if match := col_header_pat.match(line):
                current_cols = [int(v) for v in match[1].split()]
                continue

            if current_cols and (match := row_pat.match(line)):
                row_idx = int(match[1])
                vals = [float(v) for v in match[2].split()]
                for k, val in enumerate(vals):
                    if k < len(current_cols):
                        col_idx = current_cols[k]
                        data[(row_idx, col_idx)] = val
                        if row_idx > max_idx:
                            max_idx = row_idx
                        if col_idx > max_idx:
                            max_idx = col_idx

        if in_matrix and data:
            _flush()

        # Attach any trailing localization indices to last matrix
        if loc_indices and results:
            results[-1].localization_indices = loc_indices

        return results

    # ── Overlap integration matrices ─────────────────────────────────

    @classmethod
    def parse_overlap_matrices(
        cls,
        stdout: str,
        integrated_property: str,
    ) -> list[OverlapIntegrationMatrix]:
        """Extract positive/negative/all overlap integration matrices."""
        results: list[OverlapIntegrationMatrix] = []

        header_pat = re.compile(
            r"\*+\s*Integration of\s+(positive|negative|all)\s+"
            r"values in overlap region\s*\*+",
            re.I,
        )
        col_header_pat = re.compile(r"^\s+(\d+(?:\s+\d+)*)\s*$")
        row_pat = re.compile(rf"^\s+(\d+)((?:\s+{FLOAT_PATTERN})+)\s*$")
        sum_diag_pat = re.compile(
            rf"Summing up diagonal.*?:\s+({FLOAT_PATTERN})"
        )
        sum_nondiag_pat = re.compile(
            rf"Summing up non-diagonal.*?:\s+({FLOAT_PATTERN})"
        )
        sum_all_pat = re.compile(rf"Summing up all.*?:\s+({FLOAT_PATTERN})")

        in_matrix = False
        category = ""
        current_cols: list[int] = []
        data: dict[tuple[int, int], float] = {}
        max_idx = 0
        sum_d: float | None = None
        sum_nd: float | None = None
        sum_a: float | None = None

        def _flush() -> None:
            nonlocal data, max_idx, category, in_matrix
            nonlocal sum_d, sum_nd, sum_a
            if data or sum_d is not None:
                results.append(
                    OverlapIntegrationMatrix(
                        category=category,
                        integrated_property=integrated_property,
                        n_atoms=max_idx,
                        data=dict(data),
                        sum_diagonal=sum_d,
                        sum_nondiagonal=sum_nd,
                        sum_all=sum_a,
                    )
                )
            data = {}
            max_idx = 0
            sum_d = None
            sum_nd = None
            sum_a = None
            in_matrix = False

        for line in stdout.split("\n"):
            if match := header_pat.search(line):
                if in_matrix:
                    _flush()
                in_matrix = True
                category = match[1].lower()
                current_cols = []
                continue

            if not in_matrix:
                continue

            if match := sum_diag_pat.search(line):
                sum_d = float(match[1])
                continue
            if match := sum_nondiag_pat.search(line):
                sum_nd = float(match[1])
                continue
            if match := sum_all_pat.search(line):
                sum_a = float(match[1])
                _flush()
                continue

            if match := col_header_pat.match(line):
                current_cols = [int(v) for v in match[1].split()]
                continue

            if current_cols and (match := row_pat.match(line)):
                row_idx = int(match[1])
                vals = [float(v) for v in match[2].split()]
                for k, val in enumerate(vals):
                    if k < len(current_cols):
                        col_idx = current_cols[k]
                        data[(row_idx, col_idx)] = val
                        if row_idx > max_idx:
                            max_idx = row_idx
                        if col_idx > max_idx:
                            max_idx = col_idx

        if in_matrix:
            _flush()

        return results

    # ── CLRK matrix ──────────────────────────────────────────────────

    @staticmethod
    def parse_clrk_matrix(stdout: str) -> CLRKMatrix | None:
        """Extract condensed linear response kernel matrix."""
        header_pat = re.compile(
            r"\*+\s*.*?(?:CLRK|[Cc]ondensed linear response "
            r"kernel).*?matrix\s*\*+",
            re.I,
        )
        col_header_pat = re.compile(r"^\s+(\d+(?:\s+\d+)*)\s*$")
        row_pat = re.compile(rf"^\s+(\d+)((?:\s+{FLOAT_PATTERN})+)\s*$")

        in_matrix = False
        current_cols: list[int] = []
        data: dict[tuple[int, int], float] = {}
        max_idx = 0

        for line in stdout.split("\n"):
            if header_pat.search(line):
                in_matrix = True
                current_cols = []
                data = {}
                max_idx = 0
                continue

            if not in_matrix:
                continue

            # End on next section or blank
            if line.strip() and not re.match(r"^\s+[\d-]", line):
                if match := col_header_pat.match(line):
                    current_cols = [int(v) for v in match[1].split()]
                    continue
                if data:
                    break

            if current_cols and (match := row_pat.match(line)):
                row_idx = int(match[1])
                vals = [float(v) for v in match[2].split()]
                for k, val in enumerate(vals):
                    if k < len(current_cols):
                        col_idx = current_cols[k]
                        data[(row_idx, col_idx)] = val
                        if row_idx > max_idx:
                            max_idx = row_idx
                        if col_idx > max_idx:
                            max_idx = col_idx

        if data:
            return CLRKMatrix(n_atoms=max_idx, data=data)
        return None

    # ── FLU references ───────────────────────────────────────────────

    @staticmethod
    def parse_flu_references(
        stdout: str,
    ) -> list[FLUReferenceParameter]:
        """Extract FLU reference bond order parameters."""
        pattern = re.compile(
            rf"([A-Za-z]+)\s*-\s*([A-Za-z]+)\s*:\s+({FLOAT_PATTERN})"
        )
        results: list[FLUReferenceParameter] = []
        in_section = False
        for line in stdout.split("\n"):
            if "FLU reference parameters" in line:
                in_section = True
                continue
            if in_section:
                if match := pattern.search(line):
                    results.append(
                        FLUReferenceParameter(
                            element1=match[1].strip(),
                            element2=match[2].strip(),
                            reference_value=float(match[3]),
                        )
                    )
                elif line.strip() and not pattern.search(line):
                    in_section = False
        return results

    # ── Aromaticity indices ──────────────────────────────────────────

    @staticmethod
    def parse_aromaticity_index(
        stdout: str,
    ) -> list[AromaticityIndex]:
        """Extract aromaticity indices (PDI, FLU, etc.)."""
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

    # ── Pairwise delocalization indices ──────────────────────────────

    @staticmethod
    def parse_delocalization_indices(
        stdout: str,
    ) -> list[DelocalizationIndex]:
        """Extract pairwise delocalization indices."""
        di_pattern = (
            rf"Delocalization index.*?atom\s+(\d+).*?atom\s+(\d+)"
            rf".*?:\s+({FLOAT_PATTERN})"
        )
        indices: list[DelocalizationIndex] = []
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


# =============================================================================
# Menu 17: Basin analysis
# =============================================================================


class BasinParser(OutputParser):
    """Parser for basin analysis output."""

    @classmethod
    def parse_for_result(
        cls,
        analysis: Menu,
        stdout: str,
    ) -> list[ParsedMultiwfnResult]:
        results: list[ParsedMultiwfnResult] = []
        results.extend(cls.parse(stdout))
        results.extend(cls.parse_charges(stdout))
        return results

    @staticmethod
    def parse(stdout: str) -> list[Basin]:
        """Extract basin integration results."""
        basins: list[Basin] = []

        pattern = (
            rf"Basin\s+(\d+).*?(?:attractor.*?atom\s+(\d+)\s*\(([^)]+)\))?"
            rf".*?population[=:\s]+({FLOAT_PATTERN})"
        )
        for match in re.finditer(pattern, stdout, re.IGNORECASE):
            basin = Basin(
                basin_id=int(match[1]),
                population=float(match[4]),
            )
            if match[2]:
                basin.attractor_atom = int(match[2])
                basin.attractor_element = match[3].strip()
            basins.append(basin)

        if not basins:
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
        """Extract AIM/Bader charges from basin analysis."""
        pattern = (
            rf"(?:AIM|Bader)\s+charge.*?atom\s+(\d+).*?:\s+({FLOAT_PATTERN})"
        )
        return [
            Charge(atom_id=int(match[1]), charge=float(match[2]))
            for match in re.finditer(pattern, stdout, re.IGNORECASE)
        ]


# =============================================================================
# Menu 18: Electron excitation analysis
# =============================================================================


class ExcitationParser(OutputParser):
    """Parser for electron excitation analysis output."""

    @classmethod
    def parse_for_result(
        cls,
        analysis: Menu,
        stdout: str,
    ) -> list[ParsedMultiwfnResult]:
        default = [
            cls.parse_hole_electron,
            cls.parse_delta_r,
            cls.parse_lambda_index,
        ]
        case_map: dict[Menu, list] = {
            Menu.HOLE_ELECTRON_ANALYSIS: [cls.parse_hole_electron],
            Menu.CHARGE_TRANSFER_ANALYSIS: [cls.parse_charge_transfer],
            Menu.IFCT_ANALYSIS: [cls.parse_charge_transfer],
            Menu.CTS_ANALYSIS: [cls.parse_charge_transfer],
            Menu.DELTA_R_INDEX: [cls.parse_delta_r],
            Menu.LAMBDA_INDEX: [cls.parse_lambda_index],
        }
        parsers = case_map.get(analysis, default)
        return cls._collect(stdout, parsers)

    @staticmethod
    def parse_hole_electron(stdout: str) -> HoleElectron | None:
        """Extract hole-electron analysis descriptors."""
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

        required = {
            "H_index",
            "E_index",
            "t_index",
            "EDI",
            "HDI",
            "Sr",
            "D_index",
        }
        if (
            required.issubset(result)
            and "hole_centroid" in centroids
            and "electron_centroid" in centroids
        ):
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
        """Extract charge transfer analysis results."""
        ct_distance: float | None = None
        ct_amount: float | None = None

        ct_dist = rf"CT\s+distance[=:\s]+({FLOAT_PATTERN})"
        ct_amt = (
            rf"(?:transferred|CT)\s+(?:charge|amount)[=:\s]+({FLOAT_PATTERN})"
        )

        if match := re.search(ct_dist, stdout, re.IGNORECASE):
            ct_distance = float(match[1])
        if match := re.search(ct_amt, stdout, re.IGNORECASE):
            ct_amount = float(match[1])

        result = ChargeTransfer(
            distance=ct_distance,
            transfer_amount=ct_amount,
        )

        frag_pattern = (
            rf"Fragment\s+(\d+).*?hole[=:\s]+({FLOAT_PATTERN})"
            rf".*?electron[=:\s]+({FLOAT_PATTERN})"
        )
        fragments: list[ChargeTransferFragment] = []
        fragments.extend(
            ChargeTransferFragment(
                fragment_id=int(match[1]),
                hole_contribution=float(match[2]),
                electron_contribution=float(match[3]),
            )
            for match in re.finditer(frag_pattern, stdout, re.IGNORECASE)
        )
        if fragments:
            result.fragments = fragments

        return result

    @staticmethod
    def parse_delta_r(stdout: str) -> list[DeltaR]:
        """Extract Delta_r index for each excited state."""
        pattern = (
            rf"(?:State|Excited state)\s+(\d+).*?"
            rf"Delta_?r[=:\s]+({FLOAT_PATTERN})"
        )
        return [
            DeltaR(state_id=int(match[1]), delta_r=float(match[2]))
            for match in re.finditer(pattern, stdout, re.IGNORECASE)
        ]

    @staticmethod
    def parse_lambda_index(stdout: str) -> list[LambdaIndex]:
        """Extract Lambda diagnostic for each excited state."""
        pattern = (
            rf"(?:State|Excited state)\s+(\d+).*?"
            rf"Lambda[=:\s]+({FLOAT_PATTERN})"
        )
        return [
            LambdaIndex(state_id=int(match[1]), lambda_index=float(match[2]))
            for match in re.finditer(pattern, stdout, re.IGNORECASE)
        ]


# =============================================================================
# Menu 20: Weak interaction analysis
# =============================================================================


class WeakInteractionParser(OutputParser):
    """Parser for weak interaction analysis output."""

    @classmethod
    def parse_for_result(
        cls,
        analysis: Menu,
        stdout: str,
    ) -> list[ParsedMultiwfnResult]:
        results: list[ParsedMultiwfnResult] = []
        interaction = cls.parse(stdout)
        if interaction is not None:
            results.append(interaction)
        return results

    @staticmethod
    def parse(stdout: str) -> WeakInteraction | None:
        """Extract summary statistics from weak interaction analysis."""
        vals: dict[str, float] = {}

        patterns: dict[str, str] = {
            "delta_g_inter": rf"delta_?g_?inter.*?[=:\s]+({FLOAT_PATTERN})",
            "delta_g_intra": rf"delta_?g_?intra.*?[=:\s]+({FLOAT_PATTERN})",
            "isosurface_integral": (
                rf"(?:Integral|integral).*?isosurface.*?[=:\s]+({FLOAT_PATTERN})"
            ),
        }

        for key, pat in patterns.items():
            if match := re.search(pat, stdout, re.IGNORECASE):
                vals[key] = float(match[1])

        if "delta_g_inter" not in vals and "delta_g_intra" not in vals:
            return None

        interaction = WeakInteraction(
            delta_g_inter=vals.get("delta_g_inter", 0.0),
            delta_g_intra=vals.get("delta_g_intra", 0.0),
        )

        if "isosurface_integral" in vals:
            interaction.isosurface_integral = vals["isosurface_integral"]

        if cubes := [
            match[1]
            for match in re.finditer(
                r"(\S+\.cube)\s+has been generated", stdout
            )
        ]:
            interaction.cube_names = cubes

        return interaction


# =============================================================================
# Menu 21: Energy Decomposition Analysis (EDA)
# =============================================================================


class EDAParser(OutputParser):
    """Parser for energy decomposition analysis output."""

    @classmethod
    def parse_for_result(
        cls,
        analysis: Menu,
        stdout: str,
    ) -> list[ParsedMultiwfnResult]:
        case_map: dict[Menu, list] = {
            Menu.DISPERSION_ATOMIC_CONTRIBUTION: [
                cls.parse_dispersion_contributions
            ],
        }
        parsers = case_map.get(analysis, [cls.parse])
        return cls._collect(stdout, parsers)

    @staticmethod
    def parse(stdout: str) -> EnergyDecompositionAnalysis:
        """Extract EDA energy components."""
        components: dict[str, float] = {}

        patterns: dict[str, str] = {
            "electrostatic": rf"[Ee]lectrostatic.*?[=:\s]+({FLOAT_PATTERN})",
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
        """Extract per-atom dispersion energy contributions."""
        pattern = (
            rf"Atom\s+(\d+).*?(?:dispersion|D[34]).*?[=:\s]+({FLOAT_PATTERN})"
        )
        return [
            DispersionContribution(
                atom_id=int(match[1]), contribution=float(match[2])
            )
            for match in re.finditer(pattern, stdout, re.IGNORECASE)
        ]


# =============================================================================
# Menu 22: Conceptual DFT (CDFT)
# =============================================================================


class CDFTParser(OutputParser):
    """Parser for conceptual DFT output (Menu 22)."""

    @classmethod
    def parse_for_result(
        cls,
        analysis: Menu,
        stdout: str,
    ) -> list[ParsedMultiwfnResult]:
        results: list[ParsedMultiwfnResult] = []

        reactivity = cls.parse_reactivity(stdout)
        if reactivity is not None:
            results.append(reactivity)

        results.extend(cls.parse_condensed_fukui(stdout))
        results.extend(cls.parse_dual_descriptor(stdout))

        superdeloc = cls.parse_superdelocalizability(stdout)
        if superdeloc is not None:
            results.append(superdeloc)

        ow_fukui = cls.parse_orbital_weighted_fukui(stdout)
        if ow_fukui is not None:
            results.append(ow_fukui)

        results.extend(cls.parse_orbital_weight_decomposition(stdout))

        return results

    # ── Global reactivity indices ────────────────────────────────────

    @staticmethod
    def parse_reactivity(stdout: str) -> Reactivity | None:
        """Extract global reactivity indices."""
        r = Reactivity()
        found = False

        homo_pat = re.compile(
            rf"HOMO energy:\s+({FLOAT_PATTERN})\s+a\.u\.\s+"
            rf"({FLOAT_PATTERN})\s+eV"
        )
        lumo_pat = re.compile(
            rf"LUMO energy:\s+({FLOAT_PATTERN})\s+a\.u\.\s+"
            rf"({FLOAT_PATTERN})\s+eV"
        )
        mu_pat = re.compile(
            rf"Chemical potential:\s+({FLOAT_PATTERN})\s+a\.u\.\s+"
            rf"({FLOAT_PATTERN})\s+eV"
        )
        delta_pat = re.compile(
            rf"Delta parameter:\s+({FLOAT_PATTERN})\s+a\.u\.\s+"
            rf"({FLOAT_PATTERN})\s+eV"
        )
        hard_pat = re.compile(
            rf"(?:Chemical |Global )?[Hh]ardness:\s+"
            rf"({FLOAT_PATTERN})"
        )
        soft_pat = re.compile(
            rf"(?:Chemical |Global )?[Ss]oftness:\s+"
            rf"({FLOAT_PATTERN})"
        )
        electro_pat = re.compile(
            rf"[Ee]lectrophilicity.*?:\s+({FLOAT_PATTERN})"
        )
        nucleo_pat = re.compile(rf"[Nn]ucleophilicity.*?:\s+({FLOAT_PATTERN})")
        ip_pat = re.compile(
            rf"[Ii]onization potential.*?:\s+({FLOAT_PATTERN})"
        )
        ea_pat = re.compile(rf"[Ee]lectron affinity.*?:\s+({FLOAT_PATTERN})")

        if match := homo_pat.search(stdout):
            r.homo_energy_au = float(match[1])
            r.homo_energy_eV = float(match[2])
            found = True
        if match := lumo_pat.search(stdout):
            r.lumo_energy_au = float(match[1])
            r.lumo_energy_eV = float(match[2])
            found = True
        if match := mu_pat.search(stdout):
            r.chemical_potential = float(match[1])
            r.chemical_potential_eV = float(match[2])
            found = True
        if match := delta_pat.search(stdout):
            r.delta_parameter_au = float(match[1])
            r.delta_parameter_eV = float(match[2])
            found = True
        if match := hard_pat.search(stdout):
            r.hardness = float(match[1])
            found = True
        if match := soft_pat.search(stdout):
            r.softness = float(match[1])
            found = True
        if match := electro_pat.search(stdout):
            r.electrophilicity = float(match[1])
            found = True
        if match := nucleo_pat.search(stdout):
            r.nucleophilicity = float(match[1])
            found = True
        if match := ip_pat.search(stdout):
            r.ionization_potential = float(match[1])
            found = True
        if match := ea_pat.search(stdout):
            r.electron_affinity = float(match[1])
            found = True

        return r if found else None

    # ── Condensed Fukui ──────────────────────────────────────────────

    @staticmethod
    def parse_condensed_fukui(stdout: str) -> list[CondensedFukui]:
        """Extract condensed Fukui functions."""
        results: list[CondensedFukui] = []

        # Pattern for Fukui table rows
        # "     1(C )        0.14827        0.15436        0.15132"
        fukui_pat = re.compile(
            rf"^\s+(\d+)\s*\([A-Za-z]+\s*\)\s+"
            rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+"
            rf"({FLOAT_PATTERN})\s*$"
        )

        in_fukui = False
        for line in stdout.split("\n"):
            if (
                re.search(r"Atom.*?f\+.*?f-.*?f0", line, re.I)
                and "OW" not in line
            ):
                in_fukui = True
                continue
            if in_fukui:
                if match := fukui_pat.match(line):
                    results.append(
                        CondensedFukui(
                            atom_id=int(match[1]),
                            fukui_plus=float(match[2]),
                            fukui_minus=float(match[3]),
                            fukui_zero=float(match[4]),
                        )
                    )
                    continue
                if line.strip() and not fukui_pat.match(line):
                    in_fukui = False

        return results

    # ── Dual descriptor ──────────────────────────────────────────────

    @staticmethod
    def parse_dual_descriptor(stdout: str) -> list[DualDescriptor]:
        """Extract dual descriptor values."""
        results: list[DualDescriptor] = []

        dd_pat = re.compile(
            rf"^\s+(\d+)\s*\([A-Za-z]+\s*\)\s+.*?"
            rf"({FLOAT_PATTERN})\s*$"
        )

        in_dd = False
        for line in stdout.split("\n"):
            if re.search(r"Atom.*?[Dd]ual\s+[Dd]escriptor", line):
                in_dd = True
                continue
            if in_dd:
                if match := dd_pat.match(line):
                    results.append(
                        DualDescriptor(
                            atom_id=int(match[1]),
                            value=float(match[2]),
                        )
                    )
                    continue
                if line.strip() and "Sum" not in line:
                    in_dd = False

        return results

    # ── Superdelocalizability ────────────────────────────────────────

    @staticmethod
    def parse_superdelocalizability(
        stdout: str,
    ) -> SuperdelocalizabilityResult | None:
        """Extract superdelocalizability analysis."""
        alpha_pat = re.compile(
            rf"alpha parameter:\s+({FLOAT_PATTERN})\s+Hartree"
        )
        entry_pat = re.compile(
            rf"^\s+(\d+)\s*\(([A-Za-z]+)\s*\)\s+"
            rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+"
            rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})"
        )
        sum_dn_pat = re.compile(rf"Sum of D_N:\s+({FLOAT_PATTERN})")
        sum_de_pat = re.compile(rf"Sum of D_E:\s+({FLOAT_PATTERN})")
        sum_dn0_pat = re.compile(rf"Sum of D_N_0:\s+({FLOAT_PATTERN})")
        sum_de0_pat = re.compile(rf"Sum of D_E_0:\s+({FLOAT_PATTERN})")

        # Check if superdelocalizability output exists
        if "superdelocalizability" not in stdout.lower():
            return None

        result = SuperdelocalizabilityResult()

        if match := alpha_pat.search(stdout):
            result.alpha_parameter = float(match[1])

        in_table = False
        for line in stdout.split("\n"):
            if "Atom" in line and "D_N" in line and "D_E" in line:
                in_table = True
                continue
            if in_table:
                if match := entry_pat.match(line):
                    result.entries.append(
                        SuperdelocalizabilityEntry(
                            atom_id=int(match[1]),
                            atom_element=match[2].strip(),
                            d_n=float(match[3]),
                            d_e=float(match[4]),
                            d_n_0=float(match[5]),
                            d_e_0=float(match[6]),
                        )
                    )
                    continue
                if match := sum_dn_pat.search(line):
                    result.sum_d_n = float(match[1])
                    continue
                if match := sum_de_pat.search(line):
                    result.sum_d_e = float(match[1])
                    continue
                if match := sum_dn0_pat.search(line):
                    result.sum_d_n_0 = float(match[1])
                    continue
                if match := sum_de0_pat.search(line):
                    result.sum_d_e_0 = float(match[1])
                    in_table = False
                    continue

        if result.entries:
            return result
        return None

    # ── Orbital-weighted Fukui ───────────────────────────────────────

    @staticmethod
    def parse_orbital_weighted_fukui(
        stdout: str,
    ) -> OrbitalWeightedFukuiResult | None:
        """Extract orbital-weighted Fukui functions."""
        entry_pat = re.compile(
            rf"^\s+(\d+)\s*\(([A-Za-z]+)\s*\)\s+"
            rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+"
            rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})"
        )
        sum_fplus_pat = re.compile(
            rf"Sum of orbital weighted f\+\s+({FLOAT_PATTERN})"
        )
        sum_fminus_pat = re.compile(
            rf"Sum of orbital weighted f-\s+({FLOAT_PATTERN})"
        )

        entries: list[OrbitalWeightedFukuiEntry] = []
        sum_fplus: float | None = None
        sum_fminus: float | None = None

        in_table = False
        for line in stdout.split("\n"):
            if "Atom" in line and "OW f+" in line:
                in_table = True
                continue
            if in_table:
                if match := entry_pat.match(line):
                    entries.append(
                        OrbitalWeightedFukuiEntry(
                            atom_id=int(match[1]),
                            atom_element=match[2].strip(),
                            ow_f_plus=float(match[3]),
                            ow_f_minus=float(match[4]),
                            ow_f_zero=float(match[5]),
                            ow_dd=float(match[6]),
                        )
                    )
                    continue
                if match := sum_fplus_pat.search(line):
                    sum_fplus = float(match[1])
                    continue
                if match := sum_fminus_pat.search(line):
                    sum_fminus = float(match[1])
                    in_table = False
                    continue

        if entries:
            return OrbitalWeightedFukuiResult(
                entries=entries,
                sum_ow_f_plus=sum_fplus,
                sum_ow_f_minus=sum_fminus,
            )
        return None

    # ── Orbital weight decomposition ─────────────────────────────────

    @staticmethod
    def parse_orbital_weight_decomposition(
        stdout: str,
    ) -> list[OrbitalWeightDecomposition]:
        """Extract orbital weight decompositions for f+ and f-."""
        results: list[OrbitalWeightDecomposition] = []

        # "10 Highest weights in orbital-weighted f+"
        header_pat = re.compile(
            r"Highest weights in orbital-weighted (f[+-])", re.I
        )
        # "Orbital    22 (LUMO  )   Weight:  48.16 %   E_diff:     3.410 eV"
        entry_pat = re.compile(
            rf"Orbital\s+(\d+)\s+\(([^)]+)\)\s+Weight:\s+"
            rf"({FLOAT_PATTERN})\s+%\s+E_diff:\s+({FLOAT_PATTERN})\s+eV"
        )
        # "Total weight of above listed orbitals: 100.00 %"
        total_pat = re.compile(
            rf"Total weight of above.*?:\s+({FLOAT_PATTERN})\s+%"
        )

        current_type: str | None = None
        current_entries: list[OrbitalWeightEntry] = []
        current_total: float | None = None

        def _flush() -> None:
            nonlocal current_entries, current_total, current_type
            if current_entries and current_type is not None:
                results.append(
                    OrbitalWeightDecomposition(
                        fukui_type=current_type,
                        entries=list(current_entries),
                        total_weight_pct=current_total,
                    )
                )
            current_entries = []
            current_total = None
            current_type = None

        for line in stdout.split("\n"):
            if match := header_pat.search(line):
                _flush()
                current_type = match[1]
                continue

            if current_type is not None:
                if match := entry_pat.search(line):
                    current_entries.append(
                        OrbitalWeightEntry(
                            orbital_id=int(match[1]),
                            orbital_label=match[2].strip(),
                            weight_pct=float(match[3]),
                            e_diff_eV=float(match[4]),
                        )
                    )
                    continue
                if match := total_pat.search(line):
                    current_total = float(match[1])
                    _flush()
                    continue

        _flush()
        return results


# =============================================================================
# Menu 24: Polarizability
# =============================================================================


class PolarizabilityParser(OutputParser):
    """Parser for polarizability output."""

    @classmethod
    def parse_for_result(
        cls,
        analysis: Menu,
        stdout: str,
    ) -> list[ParsedMultiwfnResult]:
        return [cls.parse(stdout)]

    @staticmethod
    def parse(stdout: str) -> Polarizability:
        """Extract polarizability tensor data."""
        result = Polarizability()

        iso_pat = (
            rf"[Ii]sotropic.*?(?:alpha|polarizability)"
            rf".*?[=:\s]+({FLOAT_PATTERN})"
        )
        if match := re.search(iso_pat, stdout):
            result.isotropic = float(match[1])

        aniso_pat = rf"[Aa]nisotropy.*?[=:\s]+({FLOAT_PATTERN})"
        if match := re.search(aniso_pat, stdout):
            result.anisotropic = float(match[1])

        beta_pat = rf"[Bb]eta.*?total[=:\s]+({FLOAT_PATTERN})"
        if match := re.search(beta_pat, stdout):
            result.beta_total = float(match[1])

        gamma_pat = rf"[Gg]amma.*?(?:total|average)[=:\s]+({FLOAT_PATTERN})"
        if match := re.search(gamma_pat, stdout):
            result.gamma_total = float(match[1])

        components: dict[str, float] = {}
        for comp in ["xx", "xy", "xz", "yy", "yz", "zz"]:
            pat = rf"alpha[_\s]*{comp}[=:\s]+({FLOAT_PATTERN})"
            if match := re.search(pat, stdout, re.IGNORECASE):
                components[f"alpha_{comp}"] = float(match[1])

        if components:
            result.tensor = PolarizabilityTensor(
                alpha_xx=components.get("alpha_xx", 0.0),
                alpha_xy=components.get("alpha_xy", 0.0),
                alpha_xz=components.get("alpha_xz", 0.0),
                alpha_yy=components.get("alpha_yy", 0.0),
                alpha_yz=components.get("alpha_yz", 0.0),
                alpha_zz=components.get("alpha_zz", 0.0),
            )

        return result


# =============================================================================
# Menu 25: Aromaticity
# =============================================================================


class AromaticityParser(OutputParser):
    """Parser for aromaticity analysis output."""

    @classmethod
    def _parse_nics_scan_or_none(cls, stdout: str) -> NICSScan | None:
        """Return parsed NICS scan only if it contains data."""
        scan = cls.parse_nics_scan(stdout)
        return scan if scan.distances else None

    @classmethod
    def parse_for_result(
        cls,
        analysis: Menu,
        stdout: str,
    ) -> list[ParsedMultiwfnResult]:
        case_map: dict[Menu, list] = {
            Menu.NICS_SCAN: [cls.parse, cls._parse_nics_scan_or_none],
            Menu.NICS_1D_SCAN: [cls.parse, cls._parse_nics_scan_or_none],
        }
        parsers = case_map.get(analysis, [cls.parse])
        return cls._collect(stdout, parsers)

    @staticmethod
    def parse(stdout: str) -> Aromaticity:
        """Extract aromaticity indices."""
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
    def parse_nics_scan(stdout: str) -> NICSScan:
        """Extract NICS scan profile data."""
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
# Menu 100/200/300: Utilities
# =============================================================================


class UtilityParser(OutputParser):
    """Parser for utility function outputs."""

    @classmethod
    def parse_for_result(
        cls,
        analysis: Menu,
        stdout: str,
    ) -> list[ParsedMultiwfnResult]:
        case_map: dict[Menu, list] = {
            Menu.GEOMETRY_PROPERTIES: [
                cls.parse_bond_lengths,
                cls.parse_bond_angles,
                cls.parse_dihedral_angles,
            ],
            Menu.ELECTRIC_MULTIPOLE_MOMENTS: [
                cls.parse_electric_multipole_moment_report,
                cls.parse_dipole_moments,
                cls.parse_quadrupole_moments,
            ],
            Menu.BLA_BOA_ANALYSIS: [cls.parse_bla_boa],
        }
        default = [
            cls.parse_bond_lengths,
            cls.parse_bond_angles,
            cls.parse_dihedral_angles,
            cls.parse_dipole_moments,
            cls.parse_coordination_numbers,
        ]
        parsers = case_map.get(analysis, default)
        return cls._collect(stdout, parsers)

    @staticmethod
    def parse_bond_lengths(stdout: str) -> list[BondLength]:
        """Extract bond length information."""
        bl_pattern = (
            rf"Bond length.*?atom\s+(\d+).*?(?:and|atom)\s+(\d+).*?:\s+"
            rf"({FLOAT_PATTERN})"
        )
        return [
            BondLength(
                atom1_id=int(match[1]),
                atom2_id=int(match[2]),
                length=float(match[3]),
            )
            for match in re.finditer(bl_pattern, stdout, re.IGNORECASE)
        ]

    @staticmethod
    def parse_bond_angles(stdout: str) -> list[BondAngle]:
        """Extract bond angle information."""
        ang_pattern = (
            rf"[Aa]ngle\s+(\d+)\s*-\s*(\d+)\s*-\s*(\d+)\s*:\s+"
            rf"({FLOAT_PATTERN})"
        )
        return [
            BondAngle(
                atom1_id=int(match[1]),
                atom2_id=int(match[2]),
                atom3_id=int(match[3]),
                angle=float(match[4]),
            )
            for match in re.finditer(ang_pattern, stdout)
        ]

    @staticmethod
    def parse_dihedral_angles(stdout: str) -> list[DihedralAngle]:
        """Extract dihedral angle information."""
        dih_pattern = (
            rf"[Dd]ihedral\s+(\d+)\s*-\s*(\d+)\s*-\s*(\d+)\s*-\s*(\d+)"
            rf"\s*:\s+({FLOAT_PATTERN})"
        )
        return [
            DihedralAngle(
                atom1_id=int(match[1]),
                atom2_id=int(match[2]),
                atom3_id=int(match[3]),
                atom4_id=int(match[4]),
                angle=float(match[5]),
            )
            for match in re.finditer(dih_pattern, stdout)
        ]

    @staticmethod
    def parse_dipole_moments(stdout: str) -> DipoleMoment | None:
        """Extract electric dipole moment."""
        xyz_pat = (
            rf"[Dd]ipole.*?X[=:\s]+({FLOAT_PATTERN})\s+"
            rf"Y[=:\s]+({FLOAT_PATTERN})\s+Z[=:\s]+({FLOAT_PATTERN})"
        )
        plain_triplet_pat = (
            rf"[Dd]ipole moment\s*\((?:a\.u\.|Debye)\)\s*:\s*"
            rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})"
        )
        if match := re.search(xyz_pat, stdout):
            return DipoleMoment(
                x=float(match[1]),
                y=float(match[2]),
                z=float(match[3]),
            )
        if match := re.search(plain_triplet_pat, stdout, re.IGNORECASE):
            return DipoleMoment(
                x=float(match[1]),
                y=float(match[2]),
                z=float(match[3]),
            )
        return None

    @staticmethod
    def parse_quadrupole_moments(stdout: str) -> QuadrupoleMoment | None:
        """Extract electric quadrupole moment."""
        moments: dict[str, float] = {}
        for comp in ["XX", "XY", "XZ", "YY", "YZ", "ZZ"]:
            pat = rf"[Qq]uadrupole.*?{comp}[=:\s]+({FLOAT_PATTERN})"
            if m := re.search(pat, stdout):
                moments[comp] = float(m[1])

        if moments:
            return QuadrupoleMoment(
                xx=moments.get("XX"),
                xy=moments.get("XY"),
                xz=moments.get("XZ"),
                yy=moments.get("YY"),
                yz=moments.get("YZ"),
                zz=moments.get("ZZ"),
            )
        return None

    @staticmethod
    def parse_coordination_numbers(
        stdout: str,
    ) -> list[CoordinationNumber]:
        """Extract atomic coordination numbers."""
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
        """Extract Bond length alternation and bond order alternation."""
        bla = None
        boa = None

        bla_pat = rf"BLA[=:\s]+({FLOAT_PATTERN})"
        if match := re.search(bla_pat, stdout, re.IGNORECASE):
            bla = float(match[1])

        boa_pat = rf"BOA[=:\s]+({FLOAT_PATTERN})"
        if match := re.search(boa_pat, stdout, re.IGNORECASE):
            boa = float(match[1])

        return BLA_BOA(bla=bla, boa=boa)

    @staticmethod
    def parse_electric_multipole_moment_report(
        stdout: str,
    ) -> ElectricMultipoleMomentReport | None:
        """Extract Menu 300 electric multipole analysis report."""
        if (
            "Quadrupole moments" not in stdout
            and "Dipole moment" not in stdout
        ):
            return None

        def _vec(pattern: str) -> tuple[float, float, float] | None:
            if m := re.search(pattern, stdout, re.IGNORECASE):
                return (float(m[1]), float(m[2]), float(m[3]))
            return None

        def _section(start: str, end_markers: list[str]) -> str:
            start_idx = stdout.find(start)
            if start_idx < 0:
                return ""
            end_idx = len(stdout)
            for marker in end_markers:
                idx = stdout.find(marker, start_idx + len(start))
                if idx >= 0:
                    end_idx = min(end_idx, idx)
            return stdout[start_idx:end_idx]

        def _pairs(section_text: str, key_pattern: str) -> dict[str, float]:
            return {
                m[1].replace(" ", ""): float(m[2])
                for m in re.finditer(
                    rf"({key_pattern})\s*=\s*({FLOAT_PATTERN})",
                    section_text,
                    re.IGNORECASE,
                )
            }

        quadrupole_standard_text = _section(
            "Quadrupole moments (Standard Cartesian form):",
            ["Quadrupole moments (Traceless Cartesian form):"],
        )
        quadrupole_traceless_text = _section(
            "Quadrupole moments (Traceless Cartesian form):",
            [
                "Magnitude of the traceless quadrupole moment tensor:",
                "Quadrupole moments (Spherical harmonic form):",
            ],
        )
        quadrupole_spherical_text = _section(
            "Quadrupole moments (Spherical harmonic form):",
            ["Magnitude: |Q_2|="],
        )
        octopole_cartesian_text = _section(
            "Octopole moments (Cartesian form):",
            ["Octopole moments (Spherical harmonic form):"],
        )
        octopole_spherical_text = _section(
            "Octopole moments (Spherical harmonic form):",
            ["Magnitude: |Q_3|="],
        )
        hexadecapole_text = _section(
            "Hexadecapole moments:",
            ["Electronic spatial extent <r^2>:"],
        )

        quadrupole_traceless_magnitude = None
        if m := re.search(
            rf"Magnitude of the traceless quadrupole moment tensor:\s*"
            rf"({FLOAT_PATTERN})",
            stdout,
            re.IGNORECASE,
        ):
            quadrupole_traceless_magnitude = float(m[1])

        quadrupole_spherical_magnitude = None
        if m := re.search(
            rf"Magnitude:\s*\|Q_2\|=\s*({FLOAT_PATTERN})",
            stdout,
            re.IGNORECASE,
        ):
            quadrupole_spherical_magnitude = float(m[1])

        octopole_spherical_magnitude = None
        if m := re.search(
            rf"Magnitude:\s*\|Q_3\|=\s*({FLOAT_PATTERN})",
            stdout,
            re.IGNORECASE,
        ):
            octopole_spherical_magnitude = float(m[1])

        dipole_magnitude_au = None
        dipole_magnitude_debye = None
        if m := re.search(
            rf"Magnitude of dipole moment:\s*({FLOAT_PATTERN})\s*a\.u\.\s*"
            rf"({FLOAT_PATTERN})\s*Debye",
            stdout,
            re.IGNORECASE,
        ):
            dipole_magnitude_au = float(m[1])
            dipole_magnitude_debye = float(m[2])

        spatial_extent_r2 = None
        if m := re.search(
            rf"Electronic spatial extent\s*<r\^2>:\s*({FLOAT_PATTERN})",
            stdout,
            re.IGNORECASE,
        ):
            spatial_extent_r2 = float(m[1])

        spatial_extent_components = None
        if m := re.search(
            rf"Components of <r\^2>:\s*X=\s*({FLOAT_PATTERN})\s*"
            rf"Y=\s*({FLOAT_PATTERN})\s*Z=\s*({FLOAT_PATTERN})",
            stdout,
            re.IGNORECASE,
        ):
            spatial_extent_components = (
                float(m[1]),
                float(m[2]),
                float(m[3]),
            )

        report = ElectricMultipoleMomentReport(
            positive_charge_center_angstrom=_vec(
                rf"X,\s*Y,\s*Z of center of positive charges.*?\n\s*"
                rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})"
            ),
            negative_charge_center_angstrom=_vec(
                rf"X,\s*Y,\s*Z of center of negative charges.*?\n\s*"
                rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})"
            ),
            nuclear_dipole_au=_vec(
                rf"Dipole moment from nuclear charges\s*\(a\.u\.\):\s*"
                rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})"
            ),
            electronic_dipole_au=_vec(
                rf"Dipole moment from electrons\s*\(a\.u\.\):\s*"
                rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})"
            ),
            total_dipole_au=_vec(
                rf"Dipole moment\s*\(a\.u\.\):\s*"
                rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})"
            ),
            total_dipole_debye=_vec(
                rf"Dipole moment\s*\(Debye\):\s*"
                rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})"
            ),
            dipole_magnitude_au=dipole_magnitude_au,
            dipole_magnitude_debye=dipole_magnitude_debye,
            quadrupole_standard_cartesian=_pairs(
                quadrupole_standard_text,
                r"XX|XY|XZ|YX|YY|YZ|ZX|ZY|ZZ",
            ),
            quadrupole_traceless_cartesian=_pairs(
                quadrupole_traceless_text,
                r"XX|XY|XZ|YX|YY|YZ|ZX|ZY|ZZ",
            ),
            quadrupole_traceless_magnitude=quadrupole_traceless_magnitude,
            quadrupole_spherical_harmonic=_pairs(
                quadrupole_spherical_text,
                r"Q_2,\s*-?\d",
            ),
            quadrupole_spherical_magnitude=quadrupole_spherical_magnitude,
            octopole_cartesian=_pairs(
                octopole_cartesian_text,
                r"XXX|YYY|ZZZ|XYY|XXY|XXZ|XZZ|YZZ|YYZ|XYZ",
            ),
            octopole_spherical_harmonic=_pairs(
                octopole_spherical_text,
                r"Q_3,\s*-?\d",
            ),
            octopole_spherical_magnitude=octopole_spherical_magnitude,
            hexadecapole=_pairs(
                hexadecapole_text,
                r"XXXX|YYYY|ZZZZ|XXXY|XXXZ|YYYX|YYYZ|ZZZX|ZZZY|XXYY|XXZZ|YYZZ|XXYZ|YYXZ|ZZXY",
            ),
            electronic_spatial_extent_r2=spatial_extent_r2,
            electronic_spatial_extent_components=spatial_extent_components,
        )

        has_data = any(
            [
                report.positive_charge_center_angstrom is not None,
                report.negative_charge_center_angstrom is not None,
                report.nuclear_dipole_au is not None,
                report.electronic_dipole_au is not None,
                report.total_dipole_au is not None,
                report.total_dipole_debye is not None,
                report.dipole_magnitude_au is not None,
                bool(report.quadrupole_standard_cartesian),
                bool(report.quadrupole_traceless_cartesian),
                report.quadrupole_traceless_magnitude is not None,
                bool(report.quadrupole_spherical_harmonic),
                report.quadrupole_spherical_magnitude is not None,
                bool(report.octopole_cartesian),
                bool(report.octopole_spherical_harmonic),
                report.octopole_spherical_magnitude is not None,
                bool(report.hexadecapole),
                report.electronic_spatial_extent_r2 is not None,
                report.electronic_spatial_extent_components is not None,
            ]
        )
        return report if has_data else None


# =============================================================================
# Parser routing table
# =============================================================================


class ParserRoute:
    """Routing from Menu enums to OutputParser classes."""

    ROUTE_TABLE: dict[Menu, type[OutputParser]] = {
        # Menu 0 - View Sttructure
        Menu.VIEW_STRUCTURE: AtomListParser,
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
        Menu.CM5_CHARGE: ChargeParser,
        Menu.EEM_CHARGE: ChargeParser,
        Menu.RESP_CHARGE: ChargeParser,
        Menu.GASTEIGER_CHARGE: ChargeParser,
        Menu.MBIS_CHARGE: ChargeParser,
        # Menu 8 — orbital composition
        Menu.ORBITAL_COMPOSITION_MULLIKEN_ALL: OrbitalCompositionParser,
        Menu.ORBITAL_COMPOSITION_SCPA_ALL: OrbitalCompositionParser,
        Menu.ORBITAL_COMPOSITION_STOUT_POLITZER_ALL: OrbitalCompositionParser,
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
        Menu.TOPOLOGY_VISUALISE_CPS: CriticalPointParser,
        Menu.TOPOLOGY_CP_NUCLEAR_POSITION: CriticalPointParser,
        Menu.TOPOLOGY_CP_MIDPOINTS: CriticalPointParser,
        Menu.TOPOLOGY_CP_TRIANGLE_CENTRES: CriticalPointParser,
        Menu.TOPOLOGY_CP_PYRAMID_CENTRES: CriticalPointParser,
        Menu.TOPOLOGY_CP_SPHERE_POINTS: CriticalPointParser,
        Menu.TOPOLOGY_CP_REAL_SPACE_POINTS: CriticalPointParser,
        Menu.TOPOLOGY_CP_PATHS_3MINUS3_3MINUS1: CriticalPointParser,
        Menu.TOPOLOGY_CP_PATHS_3PLUS1_3PLUS3: CriticalPointParser,
        # Menu 10 — density of states
        Menu.PLOT_TDOS: DOSParser,
        Menu.PLOT_TDOS_OPDOS: DOSParser,
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
        Menu.FUZZY_INTEGRATE_EDENSITY: FuzzySpaceParser,
        Menu.FUZZY_MULTIPOLE: FuzzySpaceParser,
        Menu.ATOMIC_OVERLAP_MATRIX: FuzzySpaceParser,
        Menu.LOCALIZATION_DELOCALIZATION_INDEX: FuzzySpaceParser,
        Menu.PDI_AROMATICITY: FuzzySpaceParser,
        Menu.FLU_AROMATICITY: FuzzySpaceParser,
        Menu.FLU_PI_AROMATICITY: FuzzySpaceParser,
        Menu.CONDENSED_LINEAR_RESPONSE: FuzzySpaceParser,
        Menu.PARA_LINEAR_RESPONSE: FuzzySpaceParser,
        Menu.FUZZY_OVERLAP_EDENSITY: FuzzySpaceParser,
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
        Menu.CONDENSED_FUKUI: CDFTParser,
        Menu.ORBITAL_WEIGHTED_FUKUI: CDFTParser,
        Menu.GRID_FUKUI: CDFTParser,
        Menu.LOCAL_HARDNESS: CDFTParser,
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
        # Menu 5 — cube generation
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
