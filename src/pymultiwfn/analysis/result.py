"""Result container for Multiwfn job execution.

Provides the parsed-result dataclasses, the :class:`MultiwfnResult`
container that accumulates parser output, and per-molecule JSON
persistence (replacing the former ``storage.py``).
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Literal

from pymultiwfn.enums.menu import Menu

# ═════════════════════════════════════════════════════════════════════════════
# Parsed result dataclasses
# ═════════════════════════════════════════════════════════════════════════════


@dataclass
class ParsedMultiwfnResult:
    """Base class for all Multiwfn parsed result types."""

    def to_dict(self) -> dict[str, Any]:
        """Serialise to a plain dict (JSON-safe)."""
        return asdict(self)


# ── Menu 0: View structure ──────────────────────────────────────────────


@dataclass
class AtomInfo(ParsedMultiwfnResult):
    """Atom entry from the atom list output."""

    atom_id: int
    element: str
    nuclear_charge: float
    x_bohr: float
    y_bohr: float
    z_bohr: float


@dataclass
class HOMOLUMOGap(ParsedMultiwfnResult):
    """HOMO-LUMO gap information."""

    homo_index: int
    homo_energy_au: float
    homo_energy_eV: float  # noqa: N815
    lumo_index: int
    lumo_energy_au: float
    lumo_energy_eV: float  # noqa: N815
    gap_au: float
    gap_eV: float  # noqa: N815
    gap_kJ_mol: float  # noqa: N815


# ── Menu 2: Topology ────────────────────────────────────────────────────


@dataclass
class CriticalPoint(ParsedMultiwfnResult):
    """Critical point from topology summary output."""

    index: int
    x: float
    y: float
    z: float
    type: Literal["nuclear", "bond", "ring", "cage", "unknown"] = "unknown"
    nucleus_atom_id: int | None = None
    nucleus_element: str | None = None
    bonded_atom1_id: int | None = None
    bonded_atom2_id: int | None = None


@dataclass
class TopologyPath(ParsedMultiwfnResult):
    """A gradient path connecting two critical points."""

    path_id: int
    bcp_index: int
    bcp_type: Literal["nuclear", "bond", "ring", "cage", "unknown"]
    target_cp_index: int
    target_cp_type: Literal["nuclear", "bond", "ring", "cage", "unknown"]
    length_bohr: float

    @property
    def atom1_id(self) -> int:
        return self.bcp_index

    @property
    def atom2_id(self) -> int:
        return self.target_cp_index

    @property
    def path_length(self) -> float:
        return self.length_bohr


# Backward-compatible alias for older parser/test code.
BondPath = TopologyPath


@dataclass
class PoincareHopfCounts(ParsedMultiwfnResult):
    """Critical point counts and Poincare-Hopf verification."""

    nuclear: int = 0
    bond: int = 0
    ring: int = 0
    cage: int = 0
    satisfied: bool = False


# ── Menu 5: Cube ────────────────────────────────────────────────────────


@dataclass
class GridExtremum:
    """Position and value of a grid minimum or maximum."""

    value: float
    x_bohr: float
    y_bohr: float
    z_bohr: float


@dataclass
class Cube(ParsedMultiwfnResult):
    """Result for a single cube/grid calculation sequence."""

    file_name: str = ""
    origin_x: float | None = None
    origin_y: float | None = None
    origin_z: float | None = None
    end_x: float | None = None
    end_y: float | None = None
    end_z: float | None = None
    spacing_x: float | None = None
    spacing_y: float | None = None
    spacing_z: float | None = None
    x_dim: int = 0
    y_dim: int = 0
    z_dim: int = 0
    total_points: int = 0
    minimum: GridExtremum | None = None
    maximum: GridExtremum | None = None
    integral_total: float | None = None
    integral_positive: float | None = None
    integral_negative: float | None = None


# ── Menu 6: Wavefunction ────────────────────────────────────────────────
@dataclass
class GaussianTypeFunction(ParsedMultiwfnResult):
    """Single Gaussian-type function (GTF) from basis set info."""

    gtf_index: int
    center_atom_id: int
    center_element: str
    function_type: str
    exponent: float


@dataclass
class BasisFunction(ParsedMultiwfnResult):
    """Single basis function mapping to shell and GTF range."""

    basis_index: int
    shell_index: int
    center_atom_id: int
    center_element: str
    function_type: str
    gtf_start: int
    gtf_end: int


@dataclass
class CoefficientMatrix(ParsedMultiwfnResult):
    """MO coefficient matrix (basis x orbitals) for one sequence.

    Rows are basis functions, columns are orbital indices.
    Stored as a flat dict mapping ``(row, col)`` to the coefficient
    value, to avoid importing numpy.
    """

    label: str = ""
    n_basis: int = 0
    n_orbitals: int = 0
    data: dict[tuple[int, int], float] = field(default_factory=dict)


@dataclass
class DensityMatrix(ParsedMultiwfnResult):
    """Density matrix (symmetric, basis x basis) for one sequence.

    Lower-triangle storage: keys ``(row, col)`` with ``row >= col``.
    """

    label: str = ""
    n_basis: int = 0
    data: dict[tuple[int, int], float] = field(default_factory=dict)
    trace: float | None = None
    trace_overlap: float | None = None


@dataclass
class ExportedMatrix(ParsedMultiwfnResult):
    """Record that a matrix was exported to a file."""

    file_name: str = ""
    label: str = ""


# ── Menu 7: Charges ──────────────────────────────────────────────────────


@dataclass
class Charge(ParsedMultiwfnResult):
    """Single atomic charge entry."""

    atom_id: int
    charge: float


@dataclass
class ChargeSet(ParsedMultiwfnResult):
    """A complete set of atomic charges from one method.

    Parameters
    ----------
    method
        Standardized method identifier.  One of: ``"mulliken"``,
        ``"lowdin"``, ``"hirshfeld"``, ``"vdd"``, ``"scpa"``,
        ``"stout_politzer"``, ``"bickelhaupt"``, ``"becke"``,
        ``"adch"``, ``"chelpg"``, ``"mk"``, ``"aim"``,
        ``"hirshfeld_i"``, ``"cm5"``, ``"eem"``, ``"resp"``,
        ``"gasteiger"``, ``"mbis"``, ``"ddec"``,
        ``"esp"`` (for raw ESP fit), ``"population"`` (for
        Mulliken/Lowdin gross populations).
    stage
        For methods with multiple stages (e.g. RESP stage 1/2),
        the stage identifier.  ``"final"`` for the final normalized
        charges.  ``"corrected"`` for ADC-corrected intermediate.
        ``"raw"`` for the raw/unnormalized charges.  ``None`` if
        there is only one stage.
    charges
        Per-atom charges.
    total_charge
        Sum of all charges, if reported in the output.
    """

    method: str
    stage: str | None = None
    charges: list[Charge] = field(default_factory=list)
    total_charge: float | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "method": self.method,
            "stage": self.stage,
            "charges": [asdict(c) for c in self.charges],
            "total_charge": self.total_charge,
        }


@dataclass
class Dipole(ParsedMultiwfnResult):
    """Result for dipole moment analysis."""

    x: float
    y: float
    z: float
    total: float


# ── Menu 8: Orbital composition ─────────────────────────────────────────
@dataclass
class OrbitalBasisComposition(ParsedMultiwfnResult):
    """Per-basis-function contribution to an orbital.

    Each orbital that has detailed composition output produces one
    of these containers.  The ``method`` field distinguishes which
    decomposition scheme was used (Mulliken/SCPA/Stout-Politzer).
    """

    method: str
    orbital_id: int
    energy_au: float
    occupation: float
    orbital_type: str = ""

    # Per-basis contributions
    basis_contributions: list[OrbitalBasisEntry] = field(default_factory=list)

    # Per-shell contributions
    shell_contributions: list[OrbitalShellEntry] = field(default_factory=list)

    # Per-shell-type (angular momentum) contributions
    shell_type_composition: OrbitalShellTypeComposition | None = None

    # Per-atom contributions
    atom_contributions: list[OrbitalAtomEntry] = field(default_factory=list)

    # Summary
    delocalization_index: float | None = None


@dataclass
class OrbitalBasisEntry:
    """Single basis function contribution to an orbital."""

    basis_index: int
    function_type: str
    atom_id: int
    atom_element: str
    shell_index: int
    local_pct: float
    cross_pct: float
    total_pct: float


@dataclass
class OrbitalShellEntry:
    """Single shell contribution to an orbital."""

    shell_index: int
    shell_type: str
    atom_id: int
    atom_element: str
    composition_pct: float


@dataclass
class OrbitalShellTypeComposition:
    """Angular momentum type composition of an orbital."""

    s: float = 0.0
    p: float = 0.0
    d: float = 0.0
    f: float = 0.0
    g: float = 0.0
    h: float = 0.0


@dataclass
class OrbitalAtomEntry:
    """Single atom contribution to an orbital."""

    atom_id: int
    atom_element: str
    composition_pct: float


@dataclass
class OrbitalAtomComposition(ParsedMultiwfnResult):
    """Per-atom contribution to an orbital (Hirshfeld/Becke/fragment methods).

    Simpler than ``OrbitalBasisComposition`` — only atom-level
    contributions are available, no basis/shell breakdown.
    """

    method: str
    orbital_id: int
    energy_au: float
    occupation: float
    orbital_type: str = ""

    atom_contributions: list[OrbitalAtomEntry] = field(default_factory=list)
    delocalization_index: float | None = None
    sum_before_normalization: float | None = None


# Backward-compatible alias for older integration tests.
OrbitalComponent = OrbitalBasisComposition

# ── Menu 9: Bond orders ─────────────────────────────────────────────────


@dataclass
class BondOrder(ParsedMultiwfnResult):
    """Single bond order entry between two atoms."""

    atom1_id: int
    atom2_id: int
    bond_order: float


@dataclass
class BondOrderSet(ParsedMultiwfnResult):
    """A complete set of bond orders from one method.

    Parameters
    ----------
    method
        Standardized method identifier.  One of: ``"mayer"``,
        ``"wiberg"``, ``"mulliken"``, ``"fuzzy"``,
        ``"laplacian"``.
    threshold
        The printing threshold used, if reported (e.g. 0.05).
    bond_orders
        Per-pair bond orders.
    """

    method: str
    threshold: float | None = None
    bond_orders: list[BondOrder] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "method": self.method,
            "threshold": self.threshold,
            "bond_orders": [asdict(b) for b in self.bond_orders],
        }


@dataclass
class BondOrderDecomposition(ParsedMultiwfnResult):
    """Result for bond order decomposition (per orbital) analysis."""

    orbital_id: int
    contribution: float


@dataclass
class Valence(ParsedMultiwfnResult):
    """Result for valence analysis."""

    atom_id: int
    type: Literal["total_valence", "free_valence"]
    valence: float


@dataclass
class MultiCenterBondOrder(ParsedMultiwfnResult):
    """Result for multi-center bond order analysis."""

    atom_ids: list[int] = field(default_factory=list)
    bond_order: float = 0.0


@dataclass
class IBSIEntry(ParsedMultiwfnResult):
    """Single IBSI analysis entry between two atoms."""

    atom1_id: int
    atom2_id: int
    distance: float
    int_dg_pair: float
    ibsi: float


@dataclass
class IBSIAnalysis(ParsedMultiwfnResult):
    """Complete IBSI analysis result set."""

    entries: list[IBSIEntry] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "entries": [asdict(e) for e in self.entries],
        }


# ── Menu 10: Density of states ──────────────────────────────────────────


@dataclass
class DOSMetadata(ParsedMultiwfnResult):
    """Metadata from density of states calculation.

    Parameters
    ----------
    tdos_center_au
        Center of TDOS in Hartree, if reported.
    homo_level_au
        HOMO energy level used for the vertical dash line, in Hartree.
    """

    tdos_center_au: float | None = None
    homo_level_au: float | None = None


@dataclass
class DensityOfStates(ParsedMultiwfnResult):
    """Discrete DOS/PDOS arrays sampled on an energy grid."""

    energies_eV: list[float] = field(default_factory=list)  # noqa: N815
    dos: list[float] = field(default_factory=list)
    projected_dos: dict[str, list[float]] | None = None


@dataclass
class OrbitalEnergy(ParsedMultiwfnResult):
    """Orbital energy and occupation entry."""

    index: int
    energy_eV: float  # noqa: N815
    occupation: float


# ── Menu 11: Spectra ────────────────────────────────────────────────────


@dataclass
class Spectrum(ParsedMultiwfnResult):
    """Result for spectrum analysis."""

    maximum: list[int] | None = None
    X: list[float] | None = None
    value: list[float] | None = None
    frequencies: list[float] | None = None
    intensities: list[float] | None = None
    wavelengths: list[float] | None = None
    atom_indices: list[int] | None = None
    chemical_shifts: list[float] | None = None


@dataclass
class SpectrumExtremum(ParsedMultiwfnResult):
    """Single extremum point on a simulated spectrum curve."""

    kind: Literal["maximum", "minimum"]
    index: int
    x: float
    value: float


@dataclass
class SpectrumCurveExtrema(ParsedMultiwfnResult):
    """Collection of extrema reported for a spectrum curve."""

    extrema: list[SpectrumExtremum] = field(default_factory=list)


@dataclass
class Transition(ParsedMultiwfnResult):
    """Result for transition analysis."""

    state: int
    energy_eV: float  # noqa: N815
    wavelength_nm: float
    osc_strength: float
    rot_strength: float | None = None


@dataclass
class Orbital(ParsedMultiwfnResult):
    """Orbital descriptor used by wavefunction and DOS parsers."""

    orbital_id: int
    energy_au: float
    energy_eV: float  # noqa: N815
    occupation: float | None = None
    spin: Literal["alpha", "beta"] | None = None


@dataclass
class OxidationState(ParsedMultiwfnResult):
    """Integer oxidation state inferred for one atom."""

    atom_id: int
    oxidation_state: int


@dataclass
class Color(ParsedMultiwfnResult):
    """Result for color prediction in CIE XYZ color space and RGB values."""

    X: float
    Y: float  # noqa: E741
    Z: float
    R: int
    G: int
    B: int


# ── Menu 12: Surface analysis ───────────────────────────────────────────


@dataclass
class SurfaceGeometry(ParsedMultiwfnResult):
    """Isosurface geometry statistics.

    Extracted once per surface calculation regardless of mapped
    property.
    """

    volume_bohr3: float | None = None
    volume_angstrom3: float | None = None
    area_bohr2: float | None = None
    area_angstrom2: float | None = None
    sphericity: float | None = None
    density_g_cm3: float | None = None
    n_vertices: int | None = None
    n_edges: int | None = None
    n_facets: int | None = None


@dataclass
class SurfaceExtremum(ParsedMultiwfnResult):
    """A single surface extremum (minimum or maximum)."""

    type: Literal["min", "max"]
    index: int
    value_au: float
    value_eV: float  # noqa: N815
    value_kcal_mol: float
    x_angstrom: float
    y_angstrom: float
    z_angstrom: float
    is_global: bool = False


@dataclass
class SurfaceStatistics(ParsedMultiwfnResult):
    """Summary statistics for a mapped surface property.

    Parameters
    ----------
    mapped_property
        Standardized identifier for what was mapped onto the
        surface.  One of: ``"esp"``, ``"alie"``, ``"lea"``,
        ``"leae"``, ``"edr"``, ``"maxedr"``, ``"edensity"``,
        ``"lambda2_rho"``.
    """

    mapped_property: str

    # Global extrema
    global_min_au: float | None = None
    global_max_au: float | None = None
    global_min_kcal_mol: float | None = None
    global_max_kcal_mol: float | None = None

    # Area breakdown
    overall_area_bohr2: float | None = None
    overall_area_angstrom2: float | None = None
    positive_area_bohr2: float | None = None
    positive_area_angstrom2: float | None = None
    negative_area_bohr2: float | None = None
    negative_area_angstrom2: float | None = None

    # Averages
    overall_average_au: float | None = None
    overall_average_kcal_mol: float | None = None
    positive_average_au: float | None = None
    positive_average_kcal_mol: float | None = None
    negative_average_au: float | None = None
    negative_average_kcal_mol: float | None = None

    # Variances
    sigma2_total_au2: float | None = None
    sigma2_total_kcal_mol2: float | None = None
    positive_variance_au2: float | None = None
    positive_variance_kcal_mol2: float | None = None
    negative_variance_au2: float | None = None
    negative_variance_kcal_mol2: float | None = None

    # Descriptors
    nu: float | None = None
    sigma2_tot_times_nu_au2: float | None = None
    sigma2_tot_times_nu_kcal_mol2: float | None = None
    pi_au: float | None = None
    pi_kcal_mol: float | None = None
    mpi_eV: float | None = None  # noqa: N815
    mpi_kcal_mol: float | None = None

    # Polar surface area
    nonpolar_area_angstrom2: float | None = None
    nonpolar_area_pct: float | None = None
    polar_area_angstrom2: float | None = None
    polar_area_pct: float | None = None

    # Skewness
    overall_skewness: float | None = None
    positive_skewness: float | None = None
    negative_skewness: float | None = None


@dataclass
class SurfaceAnalysisResult(ParsedMultiwfnResult):
    """Complete surface analysis for one mapped property.

    Bundles geometry, extrema, and statistics so that each
    sequence (ESP, ALIE, etc.) is stored as a single coherent
    result that won't overwrite other sequences.

    Parameters
    ----------
    mapped_property
        Standardized identifier: ``"esp"``, ``"alie"``, ``"lea"``,
        ``"leae"``, ``"edr"``, ``"maxedr"``, ``"edensity"``,
        ``"lambda2_rho"``.
    """

    mapped_property: str
    geometry: SurfaceGeometry | None = None
    minima: list[SurfaceExtremum] = field(default_factory=list)
    maxima: list[SurfaceExtremum] = field(default_factory=list)
    statistics: SurfaceStatistics | None = None

    @property
    def area(self) -> float | None:
        if self.geometry is None:
            return None
        return self.geometry.area_angstrom2

    @property
    def volume(self) -> float | None:
        if self.geometry is None:
            return None
        return self.geometry.volume_angstrom3

    @property
    def V_S_plus(self) -> float | None:  # noqa: N802
        if self.statistics is None:
            return None
        return self.statistics.positive_surface_volume

    @property
    def V_S_minus(self) -> float | None:  # noqa: N802
        if self.statistics is None:
            return None
        return self.statistics.negative_surface_volume


# Backward-compatible alias for older integration tests.
SurfaceAnalysis = SurfaceAnalysisResult

# ── Menu 15: Fuzzy atomic space ─────────────────────────────────────────


@dataclass
class FuzzyIntegrationEntry:
    """Single atom entry in a fuzzy integration table."""

    atom_id: int
    atom_element: str
    value: float
    pct_of_sum: float
    pct_of_sum_abs: float


@dataclass
class FuzzyIntegrationResult(ParsedMultiwfnResult):
    """Result of fuzzy space integration for one property.

    Parameters
    ----------
    integrated_property
        Standardized identifier: ``"edensity"``, ``"norm_rho"``,
        ``"laplacian"``, ``"orb_wfn_homo"``, ``"orb_wfn_lumo"``,
        ``"espin_density"``, ``"kr"``, ``"gr"``, ``"esp_charges"``,
        ``"elf"``, ``"lol"``, ``"local_entropy"``, ``"esp"``,
        ``"rdg"``, ``"rdg_promolecular"``, ``"lambda2rho"``,
        ``"lambda2rho_promolecular"``, ``"alie"``, ``"edr"``,
        ``"orb_overlap_dr"``, ``"deltag_promolecular"``,
        ``"deltag_hirshfeld"``, ``"iri"``.
    """

    integrated_property: str
    entries: list[FuzzyIntegrationEntry] = field(default_factory=list)
    total_sum: float | None = None
    total_sum_abs: float | None = None


# ── Atomic multipoles ───────────────────────────────────────────────


@dataclass
class AtomicMultipole(ParsedMultiwfnResult):
    """Full atomic multipole moments for one atom."""

    atom_id: int
    atom_element: str

    # Charge / monopole
    atomic_charge: float | None = None
    monopole_moment: float | None = None

    # Dipole
    dipole_x: float | None = None
    dipole_y: float | None = None
    dipole_z: float | None = None
    dipole_norm: float | None = None

    # Contribution to molecular dipole
    mol_dipole_contrib_x: float | None = None
    mol_dipole_contrib_y: float | None = None
    mol_dipole_contrib_z: float | None = None
    mol_dipole_contrib_norm: float | None = None

    # Quadrupole (traceless Cartesian)
    quadrupole_xx: float | None = None
    quadrupole_xy: float | None = None
    quadrupole_xz: float | None = None
    quadrupole_yy: float | None = None
    quadrupole_yz: float | None = None
    quadrupole_zz: float | None = None
    quadrupole_magnitude: float | None = None

    # Spatial extent
    spatial_extent_r2: float | None = None
    spatial_extent_x: float | None = None
    spatial_extent_y: float | None = None
    spatial_extent_z: float | None = None

    # Octopole (spherical harmonic magnitudes)
    octopole_magnitude: float | None = None


@dataclass
class MolecularMultipole(ParsedMultiwfnResult):
    """Molecular dipole and multipole moments."""

    total_electrons: float | None = None
    net_charge: float | None = None

    # Dipole
    dipole_x_au: float | None = None
    dipole_y_au: float | None = None
    dipole_z_au: float | None = None
    dipole_magnitude_au: float | None = None
    dipole_x_debye: float | None = None
    dipole_y_debye: float | None = None
    dipole_z_debye: float | None = None
    dipole_magnitude_debye: float | None = None

    # Quadrupole (traceless Cartesian)
    quadrupole_xx: float | None = None
    quadrupole_xy: float | None = None
    quadrupole_xz: float | None = None
    quadrupole_yy: float | None = None
    quadrupole_yz: float | None = None
    quadrupole_zz: float | None = None
    quadrupole_magnitude: float | None = None

    # Spatial extent
    spatial_extent_r2: float | None = None
    spatial_extent_x: float | None = None
    spatial_extent_y: float | None = None
    spatial_extent_z: float | None = None

    # Octopole magnitude
    octopole_magnitude: float | None = None


# ── AOM diagnostics ──────────────────────────────────────────────────


@dataclass
class AOMDiagnostics(ParsedMultiwfnResult):
    """Atomic overlap matrix quality diagnostics."""

    error: float | None = None
    max_diagonal_deviation: float | None = None
    max_diagonal_orbital: int | None = None
    max_nondiagonal_deviation: float | None = None
    max_nondiagonal_orbitals: tuple[int, int] | None = None
    exported_file: str | None = None


# ── Delocalization / Localization index ──────────────────────────────


@dataclass
class LocalizationIndex(ParsedMultiwfnResult):
    """Localization index for one atom."""

    atom_id: int
    index: float
    atom_element: str | None = None


@dataclass
class DelocalizationIndexMatrix(ParsedMultiwfnResult):
    """Full delocalization index matrix (atom × atom).

    Diagonal elements are localization indices (or atomic valence
    for closed-shell).  The ``label`` distinguishes "Total" from
    regular DI matrices.
    """

    label: str = ""
    n_atoms: int = 0
    data: dict[tuple[int, int], float] = field(default_factory=dict)
    localization_indices: list[LocalizationIndex] = field(default_factory=list)


# ── Overlap integration matrices ─────────────────────────────────────


@dataclass
class OverlapIntegrationMatrix(ParsedMultiwfnResult):
    """Overlap region integration matrix for one sign category.

    Parameters
    ----------
    category
        One of ``"positive"``, ``"negative"``, ``"all"``.
    """

    category: str
    integrated_property: str = ""
    n_atoms: int = 0
    data: dict[tuple[int, int], float] = field(default_factory=dict)
    sum_diagonal: float | None = None
    sum_nondiagonal: float | None = None
    sum_all: float | None = None


# ── CLRK matrix ──────────────────────────────────────────────────────


@dataclass
class CLRKMatrix(ParsedMultiwfnResult):
    """Condensed linear response kernel matrix (atom × atom)."""

    n_atoms: int = 0
    data: dict[tuple[int, int], float] = field(default_factory=dict)


# ── FLU reference parameters ─────────────────────────────────────────


@dataclass
class FLUReferenceParameter(ParsedMultiwfnResult):
    """FLU reference bond order for one element pair."""

    element1: str
    element2: str
    reference_value: float


# ── Keep existing AromaticityIndex / DelocalizationIndex ─────────────


@dataclass
class AromaticityIndex(ParsedMultiwfnResult):
    """Result for aromaticity indices (Menu 15)."""

    index_name: str
    value: float


@dataclass
class DelocalizationIndex(ParsedMultiwfnResult):
    """Pairwise delocalization index between two atoms."""

    atom1_id: int
    atom2_id: int
    index: float


# ── Menu 17: Basin analysis ─────────────────────────────────────────────


@dataclass
class Basin(ParsedMultiwfnResult):
    """Result for basin analysis."""

    basin_id: int
    population: float
    attractor_atom: int | None = None
    attractor_element: str | None = None
    volume: float | None = None
    charge: float | None = None


# ── Menu 18: Excitation analysis ────────────────────────────────────────


@dataclass
class HoleElectron(ParsedMultiwfnResult):
    """Result for hole-electron analysis."""

    hole_id: float
    electron_id: float
    transition_index: float
    electron_delocalisation_index: float
    hole_delocalisation_index: float
    Sr: float
    d_index: float
    hole_centroid: tuple[float, float, float] = (0.0, 0.0, 0.0)
    electron_centroid: tuple[float, float, float] = (0.0, 0.0, 0.0)


@dataclass
class ChargeTransferFragment(ParsedMultiwfnResult):
    """Individual fragment contribution to charge transfer analysis."""

    fragment_id: int
    hole_contribution: float
    electron_contribution: float


@dataclass
class ChargeTransfer(ParsedMultiwfnResult):
    """Result for charge transfer analysis."""

    distance: float | None = None
    transfer_amount: float | None = None
    fragments: list[ChargeTransferFragment] | None = None


@dataclass
class DeltaR(ParsedMultiwfnResult):
    """Result for Delta_r index."""

    state_id: int
    delta_r: float


@dataclass
class LambdaIndex(ParsedMultiwfnResult):
    """Result for Lambda index."""

    state_id: int
    lambda_index: float


# -- Menu 19: Orbital localisation ────────────────────────────────────────


@dataclass
class LMOAtomContribution:
    """Single atom contribution to a localized molecular orbital."""

    atom_id: int
    atom_element: str
    contribution_pct: float


@dataclass
class LocalizedOrbital:
    """Character of a single localized molecular orbital.

    Parameters
    ----------
    orbital_id
        1-based LMO index.
    category
        One of ``"single_center"``, ``"two_center"``,
        ``"delocalized"``.
    contributions
        Atom contributions listed for this LMO, ordered by
        decreasing percentage.
    """

    orbital_id: int
    category: Literal["single_center", "two_center", "delocalized"]
    contributions: list[LMOAtomContribution] = field(default_factory=list)


@dataclass
class LocalizationConvergence:
    """Convergence history for one localization run.

    Parameters
    ----------
    orbital_set
        Which orbitals were localized: ``"occupied"`` or
        ``"unoccupied"`` or ``"all"``.
    n_cycles
        Number of cycles to convergence.
    final_P
        Final value of the localization functional P.
    converged
        Whether the localization converged successfully.
    """

    orbital_set: str
    n_cycles: int = 0
    final_P: float | None = None
    converged: bool = False


@dataclass
class OrbitalLocalizationResult(ParsedMultiwfnResult):
    """Complete result for one orbital localization run.

    Each ``Menu`` sequence (e.g. Pipek-Mezey/Hirshfeld/occupied)
    produces one of these.  The ``method``, ``population_scheme``,
    and ``orbital_set`` together uniquely identify the calculation.

    Parameters
    ----------
    method
        Localization method: ``"pipek_mezey"`` or ``"boys"``.
    population_scheme
        Population analysis used for Pipek-Mezey: ``"hirshfeld"``,
        ``"lowdin"``, ``"becke"``.  ``None`` for Boys.
    orbital_set
        Which orbitals: ``"occupied"``, ``"all"``.
    convergence
        Convergence history for each sub-run (occupied and/or
        unoccupied).
    occupied_lmos
        Character of occupied LMOs.
    unoccupied_lmos
        Character of unoccupied LMOs (only present for "all"
        runs).
    exported_file
        Name of the exported file, if any.
    """

    method: str
    population_scheme: str | None = None
    orbital_set: str = ""

    convergence: list[LocalizationConvergence] = field(default_factory=list)
    occupied_lmos: list[LocalizedOrbital] = field(default_factory=list)
    unoccupied_lmos: list[LocalizedOrbital] = field(default_factory=list)
    exported_file: str | None = None


# ── Menu 20: Weak interactions ──────────────────────────────────────────


@dataclass
class WeakInteraction(ParsedMultiwfnResult):
    """Result for weak interaction analysis."""

    delta_g_inter: float = 0.0
    delta_g_intra: float = 0.0
    isosurface_integral: float | None = None
    cube_names: list[str] | None = None


# ── Menu 21: EDA ────────────────────────────────────────────────────────


@dataclass
class EnergyDecompositionAnalysis(ParsedMultiwfnResult):
    """Result for energy decomposition analysis."""

    electrostatic: float | None = None
    exchange: float | None = None
    repulsion: float | None = None
    polarization: float | None = None
    dispersion: float | None = None
    orbital_interaction: float | None = None
    total_interaction: float | None = None


@dataclass
class DispersionContribution(ParsedMultiwfnResult):
    """Result for dispersion contributions."""

    atom_id: int
    contribution: float


# ── Menu 22: CDFT ───────────────────────────────────────────────────────


@dataclass
class Reactivity(ParsedMultiwfnResult):
    """Result for global reactivity indices."""

    chemical_potential: float | None = None
    chemical_potential_eV: float | None = None  # noqa: N815
    hardness: float | None = None
    softness: float | None = None
    electrophilicity: float | None = None
    nucleophilicity: float | None = None
    ionization_potential: float | None = None
    electron_affinity: float | None = None
    homo_energy_au: float | None = None
    homo_energy_eV: float | None = None  # noqa: N815
    lumo_energy_au: float | None = None
    lumo_energy_eV: float | None = None  # noqa: N815
    delta_parameter_au: float | None = None
    delta_parameter_eV: float | None = None  # noqa: N815


# ── Condensed Fukui ──────────────────────────────────────────────────


@dataclass
class CondensedFukui(ParsedMultiwfnResult):
    """Result for condensed Fukui functions."""

    atom_id: int
    fukui_plus: float | None = None
    fukui_minus: float | None = None
    fukui_zero: float | None = None


@dataclass
class DualDescriptor(ParsedMultiwfnResult):
    """Result for dual descriptor."""

    atom_id: int
    value: float


# ── Superdelocalizability ────────────────────────────────────────────


@dataclass
class SuperdelocalizabilityEntry:
    """Per-atom superdelocalizability values."""

    atom_id: int
    atom_element: str
    d_n: float
    d_e: float
    d_n_0: float
    d_e_0: float


@dataclass
class SuperdelocalizabilityResult(ParsedMultiwfnResult):
    """Complete superdelocalizability analysis.

    Contains the alpha parameter, per-atom values, and sums.
    """

    alpha_parameter: float | None = None
    entries: list[SuperdelocalizabilityEntry] = field(default_factory=list)
    sum_d_n: float | None = None
    sum_d_e: float | None = None
    sum_d_n_0: float | None = None
    sum_d_e_0: float | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "alpha_parameter": self.alpha_parameter,
            "entries": [asdict(e) for e in self.entries],
            "sum_d_n": self.sum_d_n,
            "sum_d_e": self.sum_d_e,
            "sum_d_n_0": self.sum_d_n_0,
            "sum_d_e_0": self.sum_d_e_0,
        }


# ── Orbital-weighted Fukui ───────────────────────────────────────────


@dataclass
class OrbitalWeightedFukuiEntry:
    """Per-atom orbital-weighted Fukui values."""

    atom_id: int
    atom_element: str
    ow_f_plus: float
    ow_f_minus: float
    ow_f_zero: float
    ow_dd: float


@dataclass
class OrbitalWeightedFukuiResult(ParsedMultiwfnResult):
    """Complete orbital-weighted Fukui analysis.

    Contains per-atom values and sums of f+ and f-.
    """

    entries: list[OrbitalWeightedFukuiEntry] = field(default_factory=list)
    sum_ow_f_plus: float | None = None
    sum_ow_f_minus: float | None = None


# ── Orbital weight decomposition ────────────────────────────────────


@dataclass
class OrbitalWeightEntry:
    """Single orbital's weight in f+ or f- decomposition."""

    orbital_id: int
    orbital_label: str
    weight_pct: float
    e_diff_eV: float  # noqa: N815


@dataclass
class OrbitalWeightDecomposition(ParsedMultiwfnResult):
    """Decomposition of orbital-weighted Fukui into orbital contributions.

    Parameters
    ----------
    fukui_type
        Which Fukui function: ``"f+"``, ``"f-"``.
    """

    fukui_type: str
    entries: list[OrbitalWeightEntry] = field(default_factory=list)
    total_weight_pct: float | None = None


# ── Menu 24: Polarizability ─────────────────────────────────────────────


@dataclass
class PolarizabilityTensor(ParsedMultiwfnResult):
    """Result for polarizability tensor."""

    alpha_xx: float = 0.0
    alpha_xy: float = 0.0
    alpha_xz: float = 0.0
    alpha_yy: float = 0.0
    alpha_yz: float = 0.0
    alpha_zz: float = 0.0


@dataclass
class Polarizability(ParsedMultiwfnResult):
    """Result for polarizability."""

    isotropic: float | None = None
    anisotropic: float | None = None
    beta_total: float | None = None
    gamma_total: float | None = None
    tensor: PolarizabilityTensor | None = None


# ── Menu 25: Aromaticity ────────────────────────────────────────────────


@dataclass
class Aromaticity(ParsedMultiwfnResult):
    """Result for aromaticity analysis (Menu 25)."""

    NICS: float | None = None
    NICS_ZZ: float | None = None
    NICS_1: float | None = None
    HOMA: float | None = None
    HOMAC: float | None = None
    HOMER: float | None = None
    Bird: float | None = None
    EN_GEO: float | None = None
    EN_BLA: float | None = None


@dataclass
class NICSScan(ParsedMultiwfnResult):
    """Result for Nucleus Independent Chemical Shift scan."""

    distances: list[float] = field(default_factory=list)
    values: list[float] = field(default_factory=list)


# ── Menu 100/200/300: Utilities ──────────────────────────────────────────


@dataclass
class BondLength(ParsedMultiwfnResult):
    """Result for bond length analysis."""

    atom1_id: int
    atom2_id: int
    length: float


@dataclass
class BondAngle(ParsedMultiwfnResult):
    """Result for bond angle analysis."""

    atom1_id: int
    atom2_id: int
    atom3_id: int
    angle: float


@dataclass
class DihedralAngle(ParsedMultiwfnResult):
    """Result for dihedral angle analysis."""

    atom1_id: int
    atom2_id: int
    atom3_id: int
    atom4_id: int
    angle: float


@dataclass
class DipoleMoment(Dipole):
    """Utility-menu dipole vector (same data model as :class:`Dipole`)."""


@dataclass
class QuadrupoleMoment(ParsedMultiwfnResult):
    """Result for quadrupole moment analysis."""

    xx: float | None = None
    xy: float | None = None
    xz: float | None = None
    yy: float | None = None
    yz: float | None = None
    zz: float | None = None


@dataclass
class MultipoleMoments(ParsedMultiwfnResult):
    """Result for multipole moments."""

    moments: dict[str, float] = field(default_factory=dict)


@dataclass
class ElectricMultipoleMomentReport(ParsedMultiwfnResult):
    """Comprehensive electric multipole report for Menu 300."""

    positive_charge_center_angstrom: tuple[float, float, float] | None = None
    negative_charge_center_angstrom: tuple[float, float, float] | None = None
    nuclear_dipole_au: tuple[float, float, float] | None = None
    electronic_dipole_au: tuple[float, float, float] | None = None
    total_dipole_au: tuple[float, float, float] | None = None
    total_dipole_debye: tuple[float, float, float] | None = None
    dipole_magnitude_au: float | None = None
    dipole_magnitude_debye: float | None = None
    quadrupole_standard_cartesian: dict[str, float] = field(
        default_factory=dict
    )
    quadrupole_traceless_cartesian: dict[str, float] = field(
        default_factory=dict
    )
    quadrupole_traceless_magnitude: float | None = None
    quadrupole_spherical_harmonic: dict[str, float] = field(
        default_factory=dict
    )
    quadrupole_spherical_magnitude: float | None = None
    octopole_cartesian: dict[str, float] = field(default_factory=dict)
    octopole_spherical_harmonic: dict[str, float] = field(default_factory=dict)
    octopole_spherical_magnitude: float | None = None
    hexadecapole: dict[str, float] = field(default_factory=dict)
    electronic_spatial_extent_r2: float | None = None
    electronic_spatial_extent_components: tuple[float, float, float] | None = (
        None
    )


@dataclass
class CoordinationNumber(ParsedMultiwfnResult):
    """Result for coordination numbers."""

    atom_id: int
    coordination_number: float


@dataclass
class CorrelationIndex(ParsedMultiwfnResult):
    """Result for correlation index analysis."""

    nondynamic_correlation_index: float
    dynamic_correlation_index: float
    total_correlation_index: float


@dataclass
class BLA_BOA(ParsedMultiwfnResult):  # noqa: N801
    """Result for BLA/BOA."""

    bla: float | None = None
    boa: float | None = None


# ═════════════════════════════════════════════════════════════════════════════
# MultiwfnResult — container that accumulates parsed results per analysis
# ═════════════════════════════════════════════════════════════════════════════


class MultiwfnResult:
    """Container for parsed Multiwfn analysis results.

    Each instance is bound to a single :class:`Menu` analysis type.
    After :meth:`parse` is called with raw stdout, the parsed
    :class:`ParsedMultiwfnResult` objects are available in
    :attr:`result`.

    Parameters
    ----------
    analysis
        The Menu enum member identifying this analysis.

    """

    def __init__(self, analysis: Menu) -> None:
        self.analysis: Menu = analysis
        self.result: list[ParsedMultiwfnResult] = []

    def parse(self, stdout: str) -> None:
        """Parse *stdout* using the router and store results.

        Imports :class:`ParserRoute` lazily to avoid circular imports
        (parsers.py imports dataclasses from this module).
        """
        from pymultiwfn.analysis.parsers import ParserRoute

        parser_cls = ParserRoute.ROUTE_TABLE.get(self.analysis)
        if parser_cls is None:
            return

        parsed = parser_cls.parse_for_result(self.analysis, stdout)
        if parsed:
            self.result.extend(parsed)

    # ── serialisation ────────────────────────────────────────────────────

    def to_dict(self) -> dict[str, Any]:
        """Serialise the result list to a JSON-safe dict."""
        return {
            "analysis": self.analysis.name,
            "results": [r.to_dict() for r in self.result],
        }

    def __repr__(self) -> str:
        return (
            f"MultiwfnResult(analysis={self.analysis.name!r}, "
            f"n_results={len(self.result)})"
        )


# ═════════════════════════════════════════════════════════════════════════════
# ResultStore — per-molecule JSON persistence (replaces storage.py)
# ═════════════════════════════════════════════════════════════════════════════


class ResultStore:
    """Per-molecule result store with optional JSON persistence.

    Parsed results are always cached in memory so that
    :meth:`has_result` and :meth:`get_result` work regardless of
    whether a JSON file exists on disk.

    Parameters
    ----------
    input_file
        Path to the wavefunction input file.
    work_dir
        Directory used for the working data.  Created if it does not
        exist.
    json_path
        Controls JSON persistence:

        * ``None`` (the default) — results are cached in memory only;
          no file is written to disk.
        * A :class:`~pathlib.Path` — results are written to that exact
          file every time :meth:`store` is called.  If the file already
          exists it is loaded on construction so that previous results
          are available immediately.

        The value can be changed at any time via the :attr:`json_path`
        property.  Setting it from ``None`` to a ``Path`` will
        immediately flush the current in-memory data to disk.

    """

    def __init__(
        self,
        input_file: Path,
        work_dir: Path,
        json_path: Path | None = None,
    ) -> None:
        self._input_file = Path(input_file)
        self._work_dir = Path(work_dir)
        self._work_dir.mkdir(parents=True, exist_ok=True)

        self._json_path: Path | None = (
            Path(json_path) if json_path is not None else None
        )
        self._data: dict[str, Any] = self._load()

    @property
    def json_path(self) -> Path | None:
        """Path to the JSON file, or ``None`` if persistence is disabled."""
        return self._json_path

    @json_path.setter
    def json_path(self, value: Path | None) -> None:
        self._json_path = Path(value) if value is not None else None
        if self._json_path is not None:
            # Flush current in-memory data to disk immediately.
            self._save()

    @property
    def data(self) -> dict[str, Any]:
        """Return a *copy* of the stored data."""
        return dict(self._data)

    # ── persistence ──────────────────────────────────────────────────────

    def _load(self) -> dict[str, Any]:
        if self._json_path is not None and self._json_path.exists():
            with Path.open(self._json_path, encoding="utf-8") as f:
                return json.load(f)
        return {
            "input_file": str(self._input_file),
            "analyses": {},
        }

    def _save(self) -> None:
        if self._json_path is None:
            return
        self._json_path.parent.mkdir(parents=True, exist_ok=True)
        with Path.open(self._json_path, "w", encoding="utf-8") as f:
            json.dump(self._data, f, indent=2, default=str)

    # ── read / write ─────────────────────────────────────────────────────

    def has_result(self, analysis: Menu) -> bool:
        """Check whether a parsed result already exists for *analysis*."""
        return analysis.name in self._data.get("analyses", {})

    def get_result(self, analysis: Menu) -> dict[str, Any] | None:
        """Retrieve a previously stored parsed result, or ``None``."""
        return self._data.get("analyses", {}).get(analysis.name)

    def store(self, mwfn_result: MultiwfnResult) -> None:
        """Cache a :class:`MultiwfnResult` in memory.

        If :attr:`json_path` is not ``None``, the result is also
        written to the JSON file on disk.

        Parameters
        ----------
        mwfn_result
            A fully parsed ``MultiwfnResult`` whose ``.result`` list
            is non-empty.

        """
        if not mwfn_result.result:
            return

        entry: dict[str, Any] = {
            "parsed": mwfn_result.to_dict(),
            "timestamp": datetime.now().isoformat(),
        }
        self._data.setdefault("analyses", {})[mwfn_result.analysis.name] = (
            entry
        )
        self._save()

    def store_from_stdout(
        self,
        analysis: Menu,
        stdout: str,
    ) -> MultiwfnResult | None:
        """Parse *stdout*, persist, and return the :class:`MultiwfnResult`.

        Convenience method combining :meth:`MultiwfnResult.parse` and
        :meth:`store` in a single call.

        Returns ``None`` if no parser is available or parsing yields no
        results.
        """
        result = MultiwfnResult(analysis=analysis)
        result.parse(stdout)
        if not result.result:
            return None
        self.store(result)
        return result
