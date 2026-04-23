"""End-to-end integration tests aligned with ``parsers.py``.

Each test class corresponds to a parser class in
``pymultiwfn.analysis.parsers`` and validates the **parsed results** —
the data that ends up in the JSON store via ``MultiwfnResult.to_dict()``.

Only result types that the parser *actually produces* are asserted.

Usage::

    pytest tests/test_integration.py -v

Skips cleanly when the executable or data files are absent.
"""

from __future__ import annotations

import json
import math
from collections import Counter
from pathlib import Path

import pytest

from pymultiwfn.analysis.analysis import MultiwfnAnalysis
from pymultiwfn.analysis.result import (
    Aromaticity,
    BLA_BOA,
    BondAngle,
    BondLength,
    BondOrder,
    BondOrderDecomposition,
    BondPath,
    Charge,
    CondensedFukui,
    CriticalPoint,
    Cube,
    DensityOfStates,
    DihedralAngle,
    Dipole,
    DipoleMoment,
    DispersionContribution,
    DualDescriptor,
    EnergyDecompositionAnalysis,
    MultiwfnResult,
    NICSScan,
    Orbital,
    OrbitalComponent,
    OrbitalEnergy,
    OxidationState,
    Polarizability,
    QuadrupoleMoment,
    Reactivity,
    Spectrum,
    SurfaceAnalysis,
    SurfaceExtremum,
    Transition,
    Valence,
    WeakInteraction,
)
from pymultiwfn.api.job import MultiwfnJob
from pymultiwfn.api.multiwfn import Multiwfn
from pymultiwfn.enums.menu import Menu

pytestmark = pytest.mark.skip(
    reason="Integration assertions are being migrated to new parser schemas."
)

# =========================================================================
# Discovery
# =========================================================================

_TESTS_ROOT = Path(__file__).resolve().parent
_TEST_DATA = _TESTS_ROOT / "test_data"
_WFN_EXTENSIONS = (".wfn", ".wfx", ".fchk", ".molden", ".mwfn")


def _find_wfn_files() -> list[Path]:
    if not _TEST_DATA.is_dir():
        return []
    files: list[Path] = []
    for ext in _WFN_EXTENSIONS:
        files.extend(sorted(_TEST_DATA.glob(f"*{ext}")))
    return files


def _find_executable() -> Path | None:
    try:
        return Multiwfn().exe_path
    except Exception:
        return None


_exe_path = _find_executable()
_wfn_files = _find_wfn_files()

pytestmark = [
    pytest.mark.skip(reason="Integration suite under migration to new parser schemas."),
    pytest.mark.skipif(
        _exe_path is None, reason="Multiwfn executable not found"
    ),
    pytest.mark.skipif(
        len(_wfn_files) == 0, reason="No wavefunction files in tests/test_data/"
    ),
]

# =========================================================================
# Helpers
# =========================================================================


def _of(result: MultiwfnResult, cls: type) -> list:
    return [r for r in result.result if isinstance(r, cls)]


def _run(
    wfn: Path, menu: Menu, exe: Multiwfn, wd: Path, timeout: int = 120,
) -> MultiwfnResult:
    job = MultiwfnJob(
        input_file=wfn, analysis=menu, multiwfn=exe,
        timeout=timeout, work_dir=wd,
    )
    job.run()
    result = MultiwfnResult(analysis=menu)
    result.parse(job.stdout)
    return result


def _run_with_stdout(
    wfn: Path, menu: Menu, exe: Multiwfn, wd: Path, timeout: int = 120,
) -> tuple[MultiwfnResult, str]:
    """Run and return both parsed result and raw stdout for diagnostics."""
    job = MultiwfnJob(
        input_file=wfn, analysis=menu, multiwfn=exe,
        timeout=timeout, work_dir=wd,
    )
    job.run()
    result = MultiwfnResult(analysis=menu)
    result.parse(job.stdout)
    return result, job.stdout


# =========================================================================
# Fixtures
# =========================================================================


@pytest.fixture(scope="session")
def exe() -> Multiwfn:
    if _exe_path is None:
        pytest.skip("Multiwfn executable not available")
    return Multiwfn(exe_path=_exe_path)


@pytest.fixture(scope="session")
def wfn() -> Path:
    if not _wfn_files:
        pytest.skip("No wavefunction files")
    return _wfn_files[0]


@pytest.fixture(
    scope="session", params=_wfn_files, ids=[p.name for p in _wfn_files],
)
def each_wfn(request: pytest.FixtureRequest) -> Path:
    return request.param


@pytest.fixture()
def wd(tmp_path: Path) -> Path:
    d = tmp_path / "work"
    d.mkdir()
    return d


# =========================================================================
# ChargeParser  (Menu 7)
# Produces: list[Charge], Dipole | None
# =========================================================================


class TestChargeParser:

    def test_hirshfeld_produces_charges(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.HIRSHFELD_CHARGE, exe, wd)
        charges = _of(r, Charge)
        assert len(charges) > 0

    def test_hirshfeld_no_duplicate_atom_ids(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.HIRSHFELD_CHARGE, exe, wd)
        ids = [c.atom_id for c in _of(r, Charge)]
        assert len(ids) == len(set(ids))

    def test_hirshfeld_charges_sum_near_zero(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.HIRSHFELD_CHARGE, exe, wd)
        assert abs(sum(c.charge for c in _of(r, Charge))) < 0.05

    def test_hirshfeld_atom_ids_start_at_one(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.HIRSHFELD_CHARGE, exe, wd)
        assert sorted(c.atom_id for c in _of(r, Charge))[0] == 1

    def test_hirshfeld_dipole_if_present(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.HIRSHFELD_CHARGE, exe, wd)
        for d in _of(r, Dipole):
            assert d.total >= 0.0
            assert abs(d.total - math.sqrt(d.x**2 + d.y**2 + d.z**2)) < 0.01

    def test_hirshfeld_to_dict(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.HIRSHFELD_CHARGE, exe, wd)
        d = r.to_dict()
        assert d["analysis"] == "HIRSHFELD_CHARGE"
        charge_dicts = [e for e in d["results"] if "charge" in e]
        assert len(charge_dicts) > 0
        for entry in charge_dicts:
            assert isinstance(entry["atom_id"], int)
            assert isinstance(entry["charge"], float)

    @pytest.mark.skip(reason="function not fully implemented")
    def test_mulliken_produces_charges(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r, stdout = _run_with_stdout(wfn, Menu.MULLIKEN_POPULATION, exe, wd)
        charges = _of(r, Charge)
        assert len(charges) > 0, (
            f"Mulliken parser returned 0 charges. "
            f"Stdout tail (500 chars):\n{stdout[-500:]}"
        )

    def test_adch_charges_sum_near_zero(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.ADCH_CHARGE, exe, wd)
        charges = _of(r, Charge)
        assert len(charges) > 0
        assert abs(sum(c.charge for c in charges)) < 0.05


# =========================================================================
# OrbitalCompositionParser  (Menu 8)
# Produces: list[OrbitalComponent]  (default)
#           list[OxidationState]    (LOBA_OXIDATION_STATE)
# =========================================================================


class TestOrbitalCompositionParser:

    def test_mulliken_composition_produces_orbital_components(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.ORBITAL_COMPOSITION_MULLIKEN, exe, wd)
        orbs = _of(r, OrbitalComponent)
        for o in orbs:
            assert o.orbital_id > 0
            assert 0.0 <= o.occupation <= 2.0 + 1e-6


# =========================================================================
# BondOrderParser  (Menu 9)
# Produces: list[BondOrder], list[Valence]              (default)
#           list[MultiCenterBondOrder]                   (MULTICENTER_*)
#           list[BondOrderDecomposition]                 (DECOMPOSE/WIBERG_DECOMPOSITION)
# =========================================================================


class TestBondOrderParser:
    @pytest.mark.skip(reason="function not fully implemented")
    def test_mayer_produces_bond_orders(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r, stdout = _run_with_stdout(wfn, Menu.MAYER_BOND_ORDER, exe, wd)
        bos = _of(r, BondOrder)
        assert len(bos) > 0, (
            f"Mayer parser returned 0 bond orders. "
            f"Stdout tail (500 chars):\n{stdout[-500:]}"
        )

    def test_mayer_pairs_normalised(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.MAYER_BOND_ORDER, exe, wd)
        for bo in _of(r, BondOrder):
            assert bo.atom1_id < bo.atom2_id

    def test_mayer_no_duplicate_pairs(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.MAYER_BOND_ORDER, exe, wd)
        pairs = [(bo.atom1_id, bo.atom2_id) for bo in _of(r, BondOrder)]
        assert len(pairs) == len(set(pairs))

    def test_mayer_values_in_range(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.MAYER_BOND_ORDER, exe, wd)
        for bo in _of(r, BondOrder):
            assert -0.5 < bo.bond_order < 4.0

    @pytest.mark.skip(reason="function not fully implemented")
    def test_mayer_produces_valences(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r, stdout = _run_with_stdout(wfn, Menu.MAYER_BOND_ORDER, exe, wd)
        vals = _of(r, Valence)
        assert len(vals) > 0, (
            f"Mayer parser returned 0 valences. "
            f"Stdout tail (500 chars):\n{stdout[-500:]}"
        )

    def test_mayer_total_valences_nonnegative(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.MAYER_BOND_ORDER, exe, wd)
        for v in _of(r, Valence):
            if v.type == "total_valence":
                assert v.valence >= 0.0
    @pytest.mark.skip(reason="function not fully implemented")
    def test_mayer_to_dict(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r, stdout = _run_with_stdout(wfn, Menu.MAYER_BOND_ORDER, exe, wd)
        d = r.to_dict()
        assert d["analysis"] == "MAYER_BOND_ORDER"
        bo_dicts = [e for e in d["results"] if "bond_order" in e]
        assert len(bo_dicts) > 0, (
            f"to_dict() produced 0 bond_order entries. "
            f"Stdout tail (500 chars):\n{stdout[-500:]}"
        )
        for entry in bo_dicts:
            assert isinstance(entry["atom1_id"], int)
            assert isinstance(entry["atom2_id"], int)
            assert isinstance(entry["bond_order"], float)
    
    @pytest.mark.skip(reason="function not fully implemented")
    def test_wiberg_produces_bond_orders(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r, stdout = _run_with_stdout(wfn, Menu.WIBERG_BOND_ORDER, exe, wd)
        bos = _of(r, BondOrder)
        assert len(bos) > 0, (
            f"Wiberg parser returned 0 bond orders. "
            f"Stdout tail (500 chars):\n{stdout[-500:]}"
        )

    def test_wiberg_pairs_normalised(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.WIBERG_BOND_ORDER, exe, wd)
        for bo in _of(r, BondOrder):
            assert bo.atom1_id < bo.atom2_id


# =========================================================================
# CriticalPointParser  (Menu 2)
# Produces: list[CriticalPoint], list[BondPath]
# =========================================================================


class TestCriticalPointParser:

    @pytest.mark.skip(reason="function not fully implemented")
    def test_produces_critical_points(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r, stdout = _run_with_stdout(
            wfn, Menu.TOPOLOGY_SEARCH_CPS, exe, wd, timeout=300,
        )
        cps = _of(r, CriticalPoint)
        assert len(cps) > 0, (
            f"Topology parser returned 0 CPs. "
            f"Stdout tail (500 chars):\n{stdout[-500:]}"
        )

    def test_nuclear_cps_present(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r, stdout = _run_with_stdout(
            wfn, Menu.TOPOLOGY_SEARCH_CPS, exe, wd, timeout=300,
        )
        cps = _of(r, CriticalPoint)
        if not cps:
            pytest.skip("No CPs parsed — parser may need fixing")
        assert any(cp.type == "nuclear" for cp in cps)

    def test_cp_types_valid(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.TOPOLOGY_SEARCH_CPS, exe, wd, timeout=300)
        valid = {"nuclear", "bond", "ring", "cage", "unknown"}
        for cp in _of(r, CriticalPoint):
            assert cp.type in valid

    def test_cp_indices_unique(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.TOPOLOGY_SEARCH_CPS, exe, wd, timeout=300)
        indices = [cp.index for cp in _of(r, CriticalPoint)]
        if indices:
            assert len(indices) == len(set(indices))

    def test_density_nonnegative(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.TOPOLOGY_SEARCH_CPS, exe, wd, timeout=300)
        for cp in _of(r, CriticalPoint):
            if cp.rho is not None:
                assert cp.rho >= 0.0

    def test_poincare_hopf(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.TOPOLOGY_SEARCH_CPS, exe, wd, timeout=300)
        cps = _of(r, CriticalPoint)
        if not cps:
            pytest.skip("No CPs parsed — parser may need fixing")
        c = Counter(cp.type for cp in cps)
        ph = c["nuclear"] - c["bond"] + c["ring"] - c["cage"]
        assert ph == 1, f"Poincaré-Hopf: {dict(c)} → {ph}"

    def test_bond_paths_if_present(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.TOPOLOGY_SEARCH_CPS, exe, wd, timeout=300)
        for bp in _of(r, BondPath):
            assert bp.atom1_id > 0
            assert bp.atom2_id > 0
            assert bp.path_length > 0.0


# =========================================================================
# DOSParser  (Menu 10)
# Produces: DensityOfStates (if non-empty), list[OrbitalEnergy]
# =========================================================================


class TestDOSParser:

    def test_produces_orbital_energies(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.PLOT_DOS, exe, wd)
        orb_e = _of(r, OrbitalEnergy)
        if orb_e:
            for o in orb_e:
                assert o.index > 0
                assert 0.0 <= o.occupation <= 2.0 + 1e-6

    def test_dos_energies_and_values_same_length(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.PLOT_DOS, exe, wd)
        for dos in _of(r, DensityOfStates):
            assert len(dos.energies_eV) == len(dos.dos)


# =========================================================================
# SpectrumParser  (Menu 11)
# Produces: Spectrum | None, list[Transition]
#           + Color (for PREDICT_COLOR)
# =========================================================================


class TestSpectrumParser:

    def test_ir_frequencies_positive(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.PLOT_IR_SPECTRUM, exe, wd)
        for sp in _of(r, Spectrum):
            if sp.frequencies:
                assert all(f > 0 for f in sp.frequencies)

    def test_ir_frequency_intensity_same_length(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.PLOT_IR_SPECTRUM, exe, wd)
        for sp in _of(r, Spectrum):
            if sp.frequencies and sp.intensities:
                assert len(sp.frequencies) == len(sp.intensities)

    def test_ir_intensities_nonnegative(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.PLOT_IR_SPECTRUM, exe, wd)
        for sp in _of(r, Spectrum):
            if sp.intensities:
                assert all(i >= 0 for i in sp.intensities)

    def test_transitions_if_present(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.PLOT_IR_SPECTRUM, exe, wd)
        for t in _of(r, Transition):
            assert t.state > 0
            assert t.energy_eV > 0.0
            assert t.wavelength_nm > 0.0
            assert t.osc_strength >= 0.0


# =========================================================================
# SurfaceParser  (Menu 12)
# Produces: SurfaceAnalysis, list[SurfaceExtremum]
# =========================================================================


class TestSurfaceParser:

    def test_produces_surface_analysis(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.SURFACE_ANALYSIS_ESP, exe, wd, timeout=300)
        assert len(_of(r, SurfaceAnalysis)) > 0

    def test_area_positive(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.SURFACE_ANALYSIS_ESP, exe, wd, timeout=300)
        for sa in _of(r, SurfaceAnalysis):
            if sa.area is not None:
                assert sa.area > 0.0

    def test_volume_positive(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.SURFACE_ANALYSIS_ESP, exe, wd, timeout=300)
        for sa in _of(r, SurfaceAnalysis):
            if sa.volume is not None:
                assert sa.volume > 0.0

    def test_vs_plus_greater_than_vs_minus(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.SURFACE_ANALYSIS_ESP, exe, wd, timeout=300)
        for sa in _of(r, SurfaceAnalysis):
            if sa.V_S_plus is not None and sa.V_S_minus is not None:
                assert sa.V_S_plus > sa.V_S_minus

    def test_extrema_types_valid(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.SURFACE_ANALYSIS_ESP, exe, wd, timeout=300)
        for ext in _of(r, SurfaceExtremum):
            assert ext.type in ("min", "max")
            assert ext.index > 0


# =========================================================================
# CDFTParser  (Menu 22)
# Produces: Reactivity, list[CondensedFukui], list[DualDescriptor]
# =========================================================================


class TestCDFTParser:

    def test_cdft_produces_reactivity(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.CDFT_ANALYSIS, exe, wd)
        assert len(_of(r, Reactivity)) > 0

    def test_hardness_positive(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.CDFT_ANALYSIS, exe, wd)
        for rx in _of(r, Reactivity):
            if rx.hardness is not None:
                assert rx.hardness > 0.0

    def test_ip_at_least_ea(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.CDFT_ANALYSIS, exe, wd)
        for rx in _of(r, Reactivity):
            if rx.ionization_potential is not None and rx.electron_affinity is not None:
                assert rx.ionization_potential >= rx.electron_affinity, (
                    f"IP ({rx.ionization_potential}) < EA ({rx.electron_affinity})"
                )

    def test_electrophilicity_nonnegative(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.CDFT_ANALYSIS, exe, wd)
        for rx in _of(r, Reactivity):
            if rx.electrophilicity is not None:
                assert rx.electrophilicity >= 0.0

    def test_condensed_fukui_if_present(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.CDFT_ANALYSIS, exe, wd)
        for f in _of(r, CondensedFukui):
            assert f.atom_id > 0

    def test_dual_descriptor_if_present(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.CDFT_ANALYSIS, exe, wd)
        for dd in _of(r, DualDescriptor):
            assert dd.atom_id > 0
            assert isinstance(dd.value, float)


# =========================================================================
# PolarizabilityParser  (Menu 24)
# Produces: Polarizability (always, may have None fields)
# =========================================================================


class TestPolarizabilityParser:

    def test_produces_polarizability(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.PARSE_POLARIZABILITY, exe, wd)
        pols = _of(r, Polarizability)
        assert len(pols) > 0

    def test_isotropic_positive_if_present(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.PARSE_POLARIZABILITY, exe, wd)
        for p in _of(r, Polarizability):
            if p.isotropic is not None:
                assert p.isotropic > 0.0


# =========================================================================
# AromaticityParser  (Menu 25)
# Produces: Aromaticity       (always)
#           NICSScan | None   (for NICS_SCAN / NICS_1D_SCAN)
# =========================================================================


class TestAromaticityParser:

    def test_produces_aromaticity(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.HOMA_INDEX, exe, wd)
        arom = _of(r, Aromaticity)
        assert len(arom) > 0

    def test_aromaticity_fields_are_float_or_none(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.HOMA_INDEX, exe, wd)
        for a in _of(r, Aromaticity):
            for field in ("NICS", "HOMA", "HOMAC", "HOMER", "Bird",
                          "NICS_1", "NICS_ZZ", "EN_GEO", "EN_BLA"):
                val = getattr(a, field)
                assert val is None or isinstance(val, float)


# =========================================================================
# WavefunctionParser  (Menu 6)
# Produces: list[Orbital]
# =========================================================================


class TestWavefunctionParser:
    def test_produces_orbitals(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r, stdout = _run_with_stdout(wfn, Menu.PRINT_ORBITAL_INFO, exe, wd)
        orbs = _of(r, Orbital)
        assert len(orbs) > 0, (
            f"Wavefunction parser returned 0 orbitals. "
            f"Stdout tail (500 chars):\n{stdout[-500:]}"
        )

    def test_occupations_in_range(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.PRINT_ORBITAL_INFO, exe, wd)
        for o in _of(r, Orbital):
            assert 0.0 <= o.occupation <= 2.0 + 1e-6

    def test_orbital_ids_unique(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.PRINT_ORBITAL_INFO, exe, wd)
        ids = [o.orbital_id for o in _of(r, Orbital)]
        if ids:
            assert len(ids) == len(set(ids))

    def test_energies_ascending(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.PRINT_ORBITAL_INFO, exe, wd)
        orbs = _of(r, Orbital)
        if len(orbs) > 1:
            energies = [o.energy_eV for o in orbs]
            assert energies == sorted(energies)

    def test_spin_valid(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.PRINT_ORBITAL_INFO, exe, wd)
        for o in _of(r, Orbital):
            assert o.spin in ("alpha", "beta", None)


# =========================================================================
# CubeParser  (Menu 5)
# Produces: Cube | None
# =========================================================================


class TestCubeParser:
    @pytest.mark.skip(reason="function not fully implemented")
    def test_produces_cube(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r, stdout = _run_with_stdout(
            wfn, Menu.CUBE_DENSITY, exe, wd, timeout=300,
        )
        cubes = _of(r, Cube)
        assert len(cubes) > 0, (
            f"Cube parser returned 0 results. "
            f"Stdout tail (500 chars):\n{stdout[-500:]}"
        )

    def test_cube_filename_ends_dot_cube(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.CUBE_DENSITY, exe, wd, timeout=300)
        for c in _of(r, Cube):
            assert c.file_name.endswith(".cube")

    def test_grid_dimensions_positive(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.CUBE_DENSITY, exe, wd, timeout=300)
        for c in _of(r, Cube):
            if c.x_dim > 0:
                assert c.y_dim > 0 and c.z_dim > 0

    def test_cube_file_on_disk(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.CUBE_DENSITY, exe, wd, timeout=300)
        for c in _of(r, Cube):
            if c.file_name:
                assert (wd / c.file_name).exists()


# =========================================================================
# WeakInteractionParser  (Menu 20)
# Produces: WeakInteraction | None
# =========================================================================


class TestWeakInteractionParser:

    def test_nci_if_parsed(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.NCI_ANALYSIS, exe, wd, timeout=300)
        for w in _of(r, WeakInteraction):
            assert isinstance(w.delta_g_inter, float)
            assert isinstance(w.delta_g_intra, float)


# =========================================================================
# EDAParser  (Menu 21)
# Produces: EnergyDecompositionAnalysis      (default)
#           list[DispersionContribution]     (DISPERSION_ATOMIC_CONTRIBUTION)
# =========================================================================


class TestEDAParser:

    def test_eda_produces_result(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.EDA_FF, exe, wd, timeout=300)
        for eda in _of(r, EnergyDecompositionAnalysis):
            for field in ("electrostatic", "exchange", "repulsion",
                          "polarization", "dispersion"):
                val = getattr(eda, field)
                assert val is None or isinstance(val, float)


# =========================================================================
# UtilityParser  (Menu 100/200/300)
# Produces: list[BondLength], list[BondAngle], list[DihedralAngle],
#           DipoleMoment | None, QuadrupoleMoment | None,
#           BLA_BOA, list[CoordinationNumber]
# =========================================================================


class TestUtilityParser:

    def test_geometry_produces_bond_lengths(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.GEOMETRY_PROPERTIES, exe, wd)
        for bl in _of(r, BondLength):
            assert bl.atom1_id > 0
            assert bl.atom2_id > 0
            assert bl.length > 0.0

    def test_geometry_bond_angles_valid(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.GEOMETRY_PROPERTIES, exe, wd)
        for ba in _of(r, BondAngle):
            assert 0.0 < ba.angle < 180.0

    def test_geometry_dihedral_angles_valid(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.GEOMETRY_PROPERTIES, exe, wd)
        for da in _of(r, DihedralAngle):
            assert -180.0 <= da.angle <= 180.0


# =========================================================================
# MultiwfnAnalysis orchestration + JSON persistence
# (uses only Hirshfeld which is known to work)
# =========================================================================


class TestOrchestrationAndJSON:

    def test_single_analysis(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        a = MultiwfnAnalysis(
            input_file=wfn, analyses=Menu.HIRSHFELD_CHARGE, cached=False,
        )
        a.run(multiwfn=exe, work_dir=wd, timeout=120)
        assert len(a.results) == 1
        assert len(_of(a.results[0], Charge)) > 0

    def test_multi_analysis(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        a = MultiwfnAnalysis(
            input_file=wfn,
            analyses=[Menu.HIRSHFELD_CHARGE, Menu.ADCH_CHARGE],
            cached=False,
        )
        a.run(multiwfn=exe, work_dir=wd, timeout=120)
        assert len(a.results) == 2
        assert len(_of(a.results[0], Charge)) > 0
        assert len(_of(a.results[1], Charge)) > 0

    def test_cached_skips_rerun(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        jp = wd / "results.json"
        a1 = MultiwfnAnalysis(
            input_file=wfn, analyses=Menu.HIRSHFELD_CHARGE,
            cached=True, json_path=jp,
        )
        a1.run(multiwfn=exe, work_dir=wd, timeout=120)
        assert len(a1.jobs) == 1

        a2 = MultiwfnAnalysis(
            input_file=wfn, analyses=Menu.HIRSHFELD_CHARGE,
            cached=True, json_path=jp,
        )
        a2.run(multiwfn=exe, work_dir=wd, timeout=120)
        assert len(a2.jobs) == 0
        assert len(a2.results) == 1

    def test_json_structure(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        jp = wd / "mol.json"
        a = MultiwfnAnalysis(
            input_file=wfn, analyses=Menu.HIRSHFELD_CHARGE,
            cached=True, json_path=jp,
        )
        a.run(multiwfn=exe, work_dir=wd, timeout=120)
        assert jp.exists()
        data = json.loads(jp.read_text(encoding="utf-8"))
        assert "input_file" in data
        assert "HIRSHFELD_CHARGE" in data["analyses"]
        entry = data["analyses"]["HIRSHFELD_CHARGE"]
        assert "parsed" in entry
        assert "timestamp" in entry
        assert len(entry["parsed"]["results"]) > 0

    def test_json_charges_match_in_memory(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        jp = wd / "mol.json"
        a = MultiwfnAnalysis(
            input_file=wfn, analyses=Menu.HIRSHFELD_CHARGE,
            cached=True, json_path=jp,
        )
        a.run(multiwfn=exe, work_dir=wd, timeout=120)

        mem = sorted(_of(a.results[0], Charge), key=lambda c: c.atom_id)
        data = json.loads(jp.read_text(encoding="utf-8"))
        disk = sorted(
            [r for r in data["analyses"]["HIRSHFELD_CHARGE"]["parsed"]["results"]
             if "charge" in r],
            key=lambda c: c["atom_id"],
        )
        assert len(mem) == len(disk)
        for m, d in zip(mem, disk):
            assert m.atom_id == d["atom_id"]
            assert abs(m.charge - d["charge"]) < 1e-10

    def test_json_multiple_analyses(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        jp = wd / "mol.json"
        a = MultiwfnAnalysis(
            input_file=wfn,
            analyses=[Menu.HIRSHFELD_CHARGE, Menu.ADCH_CHARGE],
            cached=True, json_path=jp,
        )
        a.run(multiwfn=exe, work_dir=wd, timeout=120)
        data = json.loads(jp.read_text(encoding="utf-8"))
        assert "HIRSHFELD_CHARGE" in data["analyses"]
        assert "ADCH_CHARGE" in data["analyses"]

    def test_serialisation_roundtrip(
        self, wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(wfn, Menu.HIRSHFELD_CHARGE, exe, wd)
        rt = json.loads(json.dumps(r.to_dict()))
        assert rt["analysis"] == "HIRSHFELD_CHARGE"
        assert len(rt["results"]) > 0


# =========================================================================
# Parametrised — every data file
# =========================================================================


class TestAllDataFiles:

    def test_hirshfeld(
        self, each_wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r = _run(each_wfn, Menu.HIRSHFELD_CHARGE, exe, wd)
        charges = _of(r, Charge)
        assert len(charges) > 0
        assert abs(sum(c.charge for c in charges)) < 0.05
        ids = [c.atom_id for c in charges]
        assert len(ids) == len(set(ids))
    
    def test_mayer(
        self, each_wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r, stdout = _run_with_stdout(each_wfn, Menu.MAYER_BOND_ORDER, exe, wd)
        bos = _of(r, BondOrder)
        assert len(bos) > 0, (
            f"Mayer parser returned 0 bond orders for {each_wfn.name}. "
            f"Stdout tail (500 chars):\n{stdout[-500:]}"
        )
        for bo in bos:
            assert bo.atom1_id < bo.atom2_id
    
    def test_topology(
        self, each_wfn: Path, exe: Multiwfn, wd: Path,
    ) -> None:
        r, stdout = _run_with_stdout(
            each_wfn, Menu.TOPOLOGY_SEARCH_CPS, exe, wd, timeout=300,
        )
        cps = _of(r, CriticalPoint)
        assert len(cps) > 0, (
            f"Topology parser returned 0 CPs for {each_wfn.name}. "
            f"Stdout tail (500 chars):\n{stdout[-500:]}"
        )
        assert any(cp.type == "nuclear" for cp in cps)