"""Result container for Multiwfn job execution."""

from dataclasses import dataclass
from typing import Literal

from pymultiwfn.analysis.parsers import ParserRoute
from pymultiwfn.enums.menu import Menu


@dataclass
class ParsedMultiwfnResult:
    """Base class for Multiwfn result types."""

    pass


@dataclass
class MultiwfnResult:
    """Base class for Multiwfn result types."""

    analysis: Menu
    result: list[ParsedMultiwfnResult] = []

    def parse(
        self,
        stdout: str,
    ) -> None:
        ParserRoute.ROUTE_TABLE[self.analysis](stdout)


@dataclass
class Charge(ParsedMultiwfnResult):
    """Result for atomic charge analysis."""

    atom_id: int
    charge: float


@dataclass
class Dipole(ParsedMultiwfnResult):
    """Result for dipole moment analysis."""

    x: float
    y: float
    z: float
    total: float


@dataclass
class OrbitalContribution:
    """Individual contribution to an orbital."""

    label: str
    percentage: float


@dataclass
class OrbitalComponent(ParsedMultiwfnResult):
    """Result for orbital composition analysis."""

    orbital_id: int
    occupation: float
    energy: float
    contributions: list[OrbitalContribution]


@dataclass
class BondOrder(ParsedMultiwfnResult):
    """Result for bond order analysis."""

    atom1_id: int
    atom2_id: int
    bond_order: float


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

    atom_ids: list[int]
    bond_order: float


@dataclass
class OxidationState(ParsedMultiwfnResult):
    """Result for oxidation state analysis."""

    atom_id: int
    oxidation_state: int


@dataclass
class CriticalPoint:
    """Individual critical point."""

    index: int
    x: float
    y: float
    z: float
    rho: float | None = None
    laplacian: float | None = None
    ellipticity: float | None = None
    type: Literal["nuclear", "bond", "ring", "cage", "unknown"] = "unknown"


@dataclass
class BondPath(ParsedMultiwfnResult):
    """Result for bond path analysis."""

    atom1_id: int
    atom2_id: int
    bcp_id: int
    path_length: float


@dataclass
class ProjectedDensityOfStates(ParsedMultiwfnResult):
    """Result for projected density of states analysis."""

    orbital_id: int
    dos: float


@dataclass
class DensityOfStates(ParsedMultiwfnResult):
    """Result for density of states analysis."""

    energies_eV: list[float]  # noqa: N815
    dos: list[float]
    projected_dos: dict[str, list[float]] | None = None


@dataclass
class OrbitalEnergy(ParsedMultiwfnResult):
    """Result for orbital energy analysis."""

    index: int
    energy_eV: float  # noqa: N815
    occupation: float


@dataclass
class Spectrum(ParsedMultiwfnResult):
    """Result for spectrum analysis."""

    frequencies: list[float] | None = None
    intensities: list[float] | None = None
    wavelengths: list[float] | None = None
    atom_indices: list[int] | None = None
    chemical_shifts: list[float] | None = None


@dataclass
class Transition(ParsedMultiwfnResult):
    """Result for transition analysis."""

    state: int
    energy_eV: float  # noqa: N815
    wavelength_nm: float
    osc_strength: float
    rot_strength: float | None = None


@dataclass
class Color(ParsedMultiwfnResult):
    """Result for color prediction in CIE XYZ color space and RGB values."""

    X: float
    Y: float  # noqa: E741
    Z: float
    R: int
    G: int
    B: int


@dataclass
class SurfaceAnalysis(ParsedMultiwfnResult):
    """Result for surface analysis."""

    area: float | None = None
    volume: float | None = None
    V_S_plus: float | None = None
    V_S_minus: float | None = None
    sigma2_total: float | None = None
    nu: float | None = None
    pi: float | None = None
    V_S_max: float | None = None
    V_S_min: float | None = None
    balance: float | None = None


@dataclass
class SurfaceExtremum(ParsedMultiwfnResult):
    """Result for surface extrema."""

    type: Literal["min", "max"]
    index: int
    value: float
    x: float
    y: float
    z: float


@dataclass
class FuzzyAtomicProperty(ParsedMultiwfnResult):
    """Result for atomic properties in fuzzy space."""

    atom_id: int
    population: float | None = None
    dipole_x: float | None = None
    dipole_y: float | None = None
    dipole_z: float | None = None
    quadrupole: float | None = None


@dataclass
class DelocalizationIndex(ParsedMultiwfnResult):
    """Result for delocalization index."""

    atom1_id: int
    atom2_id: int
    index: float


@dataclass
class Basin(ParsedMultiwfnResult):
    """Result for basin analysis."""

    basin_id: int
    population: float
    attractor_atom: int | None = None
    attractor_element: str | None = None
    volume: float | None = None
    charge: float | None = None


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
    hole_centroid: tuple[float, float, float]
    electron_centroid: tuple[float, float, float]


@dataclass
class ChargeTransferFragment(ParsedMultiwfnResult):
    """Individual fragment contribution to charge transfer analysis."""

    fragment_id: int
    hole_contribution: float
    electron_contribution: float


@dataclass
class ChargeTransfer(ParsedMultiwfnResult):
    """Result for charge transfer analysis."""

    distance: float | None
    transfer_amount: float | None
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


@dataclass
class WeakInteraction(ParsedMultiwfnResult):
    """Result for weak interaction analysis."""

    delta_g_inter: float
    delta_g_intra: float
    isosurface_integral: float | None = None
    cube_names: list[str] | None = None


@dataclass
class EnergyDecompositionAnalysis(ParsedMultiwfnResult):
    """Result for energy decomposition analysis."""

    electrostatic: float | None = None
    exchange: float | None = None
    repulsion: float | None = None
    polarization: float | None = None
    dispersion: float | None = None
    total_interaction: float | None = None


@dataclass
class DispersionContribution(ParsedMultiwfnResult):
    """Result for dispersion contributions."""

    atom_id: int
    contribution: float


@dataclass
class Reactivity(ParsedMultiwfnResult):
    """Result for global reactivity indices."""

    chemical_potential: float | None = None
    hardness: float | None = None
    softness: float | None = None
    electrophilicity: float | None = None
    nucleophilicity: float | None = None
    ionization_potential: float | None = None
    electron_affinity: float | None = None


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


@dataclass
class PolarizabilityTensor(ParsedMultiwfnResult):
    """Result for polarizability tensor."""

    alpha_xx: float
    alpha_xy: float
    alpha_xz: float
    alpha_yy: float
    alpha_yz: float
    alpha_zz: float


@dataclass
class Polarizability(ParsedMultiwfnResult):
    """Result for polarizability."""

    isotripic: float | None = None
    anisotropic: float | None = None
    beta_total: float | None = None
    gamma_total: float | None = None
    tensor: PolarizabilityTensor | None = None


@dataclass
class AromaticityIndex(ParsedMultiwfnResult):
    """Result for aromaticity indices (Menu 15)."""

    index_name: str
    value: float


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

    distances: list[float]
    values: list[float]


@dataclass
class Orbital(ParsedMultiwfnResult):
    """Result for orbital info."""

    orbital_id: int
    energy_eV: float  # noqa: N815
    energy_au: float
    occupation: float
    spin: Literal["alpha", "beta"] | None


@dataclass
class Cube(ParsedMultiwfnResult):
    """Result for cube operations."""

    file_name: str
    x_dim: int
    y_dim: int
    z_dim: int
    min: float | None = None
    max: float | None = None
    mean: float | None = None
    integral: float | None = None
    std_dev: float | None = None


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
class DipoleMoment(ParsedMultiwfnResult):
    """Result for dipole moment analysis."""

    x: float
    y: float
    z: float


@dataclass
class QuadrupoleMoment(ParsedMultiwfnResult):
    """Result for quadrupole moment analysis."""

    xx: float | None
    xy: float | None
    xz: float | None
    yy: float | None
    yz: float | None
    zz: float | None


@dataclass
class MultipoleMoments(ParsedMultiwfnResult):
    """Result for multipole moments."""

    moments: dict[str, float]


@dataclass
class CoordinationNumber(ParsedMultiwfnResult):
    """Result for coordination numbers."""

    atom_id: int
    coordination_number: float


@dataclass
class BLA_BOA(ParsedMultiwfnResult):  # noqa: N801
    """Result for BLA/BOA."""

    bla: float | None = None
    boa: float | None = None


@dataclass
class LocalizationIndex(ParsedMultiwfnResult):
    """Result for localization index."""

    atom_id: int
    index: float
