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


# ── Menu 7: Charges ──────────────────────────────────────────────────────


# @dataclass
# class Orbital(ParsedMultiwfnResult):
#     """Individual orbital contribution to charge analysis."""

#     orbital_id: int
#     occupation: float
#     energy: float
#     symmetry: str


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


# ── Menu 8: Orbital composition ─────────────────────────────────────────


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
    contributions: list[OrbitalContribution] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["contributions"] = [asdict(c) for c in self.contributions]
        return d


@dataclass
class OxidationState(ParsedMultiwfnResult):
    """Result for oxidation state analysis."""

    atom_id: int
    oxidation_state: int


# ── Menu 9: Bond orders ─────────────────────────────────────────────────


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

    atom_ids: list[int] = field(default_factory=list)
    bond_order: float = 0.0


# ── Menu 2: Topology ────────────────────────────────────────────────────


@dataclass
class CriticalPoint(ParsedMultiwfnResult):
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


# ── Menu 10: Density of states ──────────────────────────────────────────


@dataclass
class DensityOfStates(ParsedMultiwfnResult):
    """Result for density of states analysis."""

    energies_eV: list[float] = field(default_factory=list)  # noqa: N815
    dos: list[float] = field(default_factory=list)
    projected_dos: dict[str, list[float]] | None = None


@dataclass
class OrbitalEnergy(ParsedMultiwfnResult):
    """Result for orbital energy analysis."""

    index: int
    energy_eV: float  # noqa: N815
    occupation: float


# ── Menu 11: Spectra ────────────────────────────────────────────────────


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


# ── Menu 12: Surface analysis ───────────────────────────────────────────


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


# ── Menu 15: Fuzzy atomic space ─────────────────────────────────────────


@dataclass
class FuzzyAtomicProperty(ParsedMultiwfnResult):
    """Result for atomic properties in fuzzy space."""

    atom_id: int
    population: float | None = None
    dipole_x: float | None = None
    dipole_y: float | None = None
    dipole_z: float | None = None
    quadrupole: float | None = None
    volume: float | None = None


@dataclass
class DelocalizationIndex(ParsedMultiwfnResult):
    """Result for delocalization index."""

    atom1_id: int
    atom2_id: int
    index: float


@dataclass
class AromaticityIndex(ParsedMultiwfnResult):
    """Result for aromaticity indices (Menu 15)."""

    index_name: str
    value: float


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

    def to_dict(self) -> dict[str, Any]:
        d: dict[str, Any] = {
            "distance": self.distance,
            "transfer_amount": self.transfer_amount,
        }
        if self.fragments is not None:
            d["fragments"] = [asdict(f) for f in self.fragments]
        else:
            d["fragments"] = None
        return d


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

    def to_dict(self) -> dict[str, Any]:
        d: dict[str, Any] = {
            "isotropic": self.isotropic,
            "anisotropic": self.anisotropic,
            "beta_total": self.beta_total,
            "gamma_total": self.gamma_total,
        }
        if self.tensor is not None:
            d["tensor"] = asdict(self.tensor)
        else:
            d["tensor"] = None
        return d


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


# ── Menu 6: Wavefunction ────────────────────────────────────────────────


@dataclass
class Orbital(ParsedMultiwfnResult):
    """Result for orbital info."""

    orbital_id: int
    energy_eV: float  # noqa: N815
    energy_au: float
    occupation: float
    spin: Literal["alpha", "beta"] | None = None


# ── Menu 5: Cube ────────────────────────────────────────────────────────


@dataclass
class Cube(ParsedMultiwfnResult):
    """Result for cube operations."""

    file_name: str = ""
    x_dim: int = 0
    y_dim: int = 0
    z_dim: int = 0
    min: float | None = None
    max: float | None = None
    mean: float | None = None
    integral: float | None = None
    std_dev: float | None = None


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
class DipoleMoment(ParsedMultiwfnResult):
    """Result for dipole moment analysis."""

    x: float
    y: float
    z: float


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
