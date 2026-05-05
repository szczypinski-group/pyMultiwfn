"""Tests for pymultiwfn.analysis.result — MultiwfnResult and dataclasses."""

import json
from pathlib import Path

import pytest

from pymultiwfn.analysis.result import (
    BondOrder,
    BondOrderSet,
    Charge,
    ChargeSet,
    CriticalPoint,
    Cube,
    GridExtremum,
    MultiwfnResult,
    Reactivity,
    ResultStore,
    Spectrum,
    SurfaceAnalysisResult,
    SurfaceGeometry,
    TopologyPath,
)
from pymultiwfn.enums.menu import Menu


class TestParsedMultiwfnResultBase:
    """Tests for the base ParsedMultiwfnResult class."""

    def test_to_dict(self) -> None:
        charge = Charge(atom_id=1, charge=-0.05)
        d = charge.to_dict()
        assert d["atom_id"] == 1
        assert d["charge"] == pytest.approx(-0.05)


class TestChargeSet:
    """Tests for ChargeSet dataclass."""

    def test_create(self) -> None:
        cs = ChargeSet(
            method="hirshfeld",
            stage="final",
            charges=[Charge(atom_id=1, charge=-0.05)],
            total_charge=-0.05,
        )
        assert cs.method == "hirshfeld"
        assert len(cs.charges) == 1

    def test_to_dict(self) -> None:
        cs = ChargeSet(
            method="mulliken",
            charges=[Charge(atom_id=1, charge=0.1)],
        )
        d = cs.to_dict()
        assert d["method"] == "mulliken"
        assert len(d["charges"]) == 1


class TestBondOrderSet:
    """Tests for BondOrderSet dataclass."""

    def test_create(self) -> None:
        bos = BondOrderSet(
            method="mayer",
            threshold=0.05,
            bond_orders=[BondOrder(atom1_id=1, atom2_id=2, bond_order=1.5)],
        )
        assert bos.method == "mayer"
        assert len(bos.bond_orders) == 1

    def test_to_dict(self) -> None:
        bos = BondOrderSet(
            method="wiberg",
            bond_orders=[BondOrder(atom1_id=1, atom2_id=2, bond_order=1.0)],
        )
        d = bos.to_dict()
        assert d["method"] == "wiberg"
        assert len(d["bond_orders"]) == 1


class TestCriticalPoint:
    """Tests for CriticalPoint dataclass."""

    def test_create_with_defaults(self) -> None:
        cp = CriticalPoint(index=1, x=0.0, y=0.0, z=0.0)
        assert cp.type == "unknown"
        assert cp.nucleus_atom_id is None

    def test_create_nuclear(self) -> None:
        cp = CriticalPoint(
            index=1,
            x=0.0,
            y=0.0,
            z=0.0,
            type="nuclear",
            nucleus_atom_id=1,
            nucleus_element="C",
        )
        assert cp.type == "nuclear"
        assert cp.nucleus_atom_id == 1


class TestTopologyPath:
    """Tests for TopologyPath dataclass."""

    def test_backward_compat_properties(self) -> None:
        tp = TopologyPath(
            path_id=1,
            bcp_index=5,
            bcp_type="bond",
            target_cp_index=1,
            target_cp_type="nuclear",
            length_bohr=2.456,
        )
        assert tp.atom1_id == 5
        assert tp.atom2_id == 1
        assert tp.path_length == pytest.approx(2.456)


class TestCube:
    """Tests for Cube dataclass."""

    def test_create_minimal(self) -> None:
        cube = Cube(file_name="density.cube")
        assert cube.file_name == "density.cube"
        assert cube.x_dim == 0

    def test_create_with_extrema(self) -> None:
        cube = Cube(
            file_name="test.cube",
            minimum=GridExtremum(
                value=0.001, x_bohr=1.0, y_bohr=2.0, z_bohr=3.0
            ),
            maximum=GridExtremum(
                value=0.999, x_bohr=4.0, y_bohr=5.0, z_bohr=6.0
            ),
        )
        assert cube.minimum.value == pytest.approx(0.001)
        assert cube.maximum.value == pytest.approx(0.999)


class TestSurfaceAnalysisResult:
    """Tests for SurfaceAnalysisResult dataclass."""

    def test_create(self) -> None:
        sar = SurfaceAnalysisResult(mapped_property="esp")
        assert sar.mapped_property == "esp"
        assert sar.geometry is None
        assert sar.minima == []

    def test_area_property(self) -> None:
        geo = SurfaceGeometry(area_angstrom2=125.5)
        sar = SurfaceAnalysisResult(mapped_property="esp", geometry=geo)
        assert sar.area == pytest.approx(125.5)

    def test_area_property_no_geometry(self) -> None:
        sar = SurfaceAnalysisResult(mapped_property="esp")
        assert sar.area is None

    def test_volume_property(self) -> None:
        geo = SurfaceGeometry(volume_angstrom3=34.7)
        sar = SurfaceAnalysisResult(mapped_property="esp", geometry=geo)
        assert sar.volume == pytest.approx(34.7)


class TestReactivity:
    """Tests for Reactivity dataclass."""

    def test_create_with_defaults(self) -> None:
        r = Reactivity()
        assert r.chemical_potential is None
        assert r.hardness is None

    def test_create_with_values(self) -> None:
        r = Reactivity(
            chemical_potential=-0.15,
            hardness=0.21,
            electrophilicity=0.05,
        )
        assert r.chemical_potential == pytest.approx(-0.15)
        assert r.hardness == pytest.approx(0.21)


class TestMultiwfnResult:
    """Tests for the MultiwfnResult container."""

    def test_create(self) -> None:
        r = MultiwfnResult(analysis=Menu.HIRSHFELD_CHARGE)
        assert r.analysis == Menu.HIRSHFELD_CHARGE
        assert r.result == []

    def test_parse_charges(self) -> None:
        r = MultiwfnResult(analysis=Menu.HIRSHFELD_CHARGE)
        stdout = "Final atomic charges:\nAtom    1(C ):    -0.05230000\n"
        r.parse(stdout)
        assert len(r.result) > 0
        # Should contain ChargeSet objects
        charge_sets = [x for x in r.result if isinstance(x, ChargeSet)]
        assert len(charge_sets) >= 1
        assert len(charge_sets[0].charges) > 0

    def test_parse_bond_orders(self) -> None:
        r = MultiwfnResult(analysis=Menu.MAYER_BOND_ORDER)
        stdout = (
            "Bond orders with absolute value >=  0.050000\n"
            "#    1:         1(C )    2(N )    1.45230000\n"
        )
        r.parse(stdout)
        bo_sets = [x for x in r.result if isinstance(x, BondOrderSet)]
        assert len(bo_sets) >= 1

    def test_parse_spectrum(self) -> None:
        r = MultiwfnResult(analysis=Menu.PLOT_IR_SPECTRUM)
        stdout = "  500.0 cm^-1  Intensity:  12.34\n"
        r.parse(stdout)
        spectrum_results = [x for x in r.result if isinstance(x, Spectrum)]
        if not None:
            assert len(spectrum_results) == 1
            assert len(spectrum_results[0].frequencies) == 1

    def test_parse_charges_empty(self) -> None:
        r = MultiwfnResult(analysis=Menu.HIRSHFELD_CHARGE)
        r.parse("")
        # Dipole is also checked, but empty stdout means no results
        assert r.result == []

    def test_to_dict(self) -> None:
        r = MultiwfnResult(analysis=Menu.HIRSHFELD_CHARGE)
        stdout = "Final atomic charges:\nAtom    1(C ):    -0.05230000\n"
        r.parse(stdout)
        d = r.to_dict()
        assert d["analysis"] == "HIRSHFELD_CHARGE"
        assert len(d["results"]) > 0

    def test_repr(self) -> None:
        r = MultiwfnResult(analysis=Menu.HIRSHFELD_CHARGE)
        assert "MultiwfnResult" in repr(r)
        assert "HIRSHFELD_CHARGE" in repr(r)


class TestResultStore:
    """Tests for the ResultStore class."""

    def test_create_without_json(self, temp_dir: Path) -> None:
        store = ResultStore(
            input_file=Path("test.wfn"),
            work_dir=temp_dir,
        )
        assert store.json_path is None

    def test_create_with_json(self, temp_dir: Path) -> None:
        jp = temp_dir / "results.json"
        store = ResultStore(
            input_file=Path("test.wfn"),
            work_dir=temp_dir,
            json_path=jp,
        )
        assert store.json_path == jp

    def test_has_result_false(self, temp_dir: Path) -> None:
        store = ResultStore(
            input_file=Path("test.wfn"),
            work_dir=temp_dir,
        )
        assert store.has_result(Menu.HIRSHFELD_CHARGE) is False

    def test_store_and_retrieve(self, temp_dir: Path) -> None:
        store = ResultStore(
            input_file=Path("test.wfn"),
            work_dir=temp_dir,
        )
        result = MultiwfnResult(analysis=Menu.HIRSHFELD_CHARGE)
        result.result.append(Charge(atom_id=1, charge=-0.05))
        store.store(result)
        assert store.has_result(Menu.HIRSHFELD_CHARGE) is True

    def test_store_persists_to_json(self, temp_dir: Path) -> None:
        jp = temp_dir / "results.json"
        store = ResultStore(
            input_file=Path("test.wfn"),
            work_dir=temp_dir,
            json_path=jp,
        )
        result = MultiwfnResult(analysis=Menu.HIRSHFELD_CHARGE)
        result.result.append(Charge(atom_id=1, charge=-0.05))
        store.store(result)
        assert jp.exists()
        data = json.loads(jp.read_text(encoding="utf-8"))
        assert "HIRSHFELD_CHARGE" in data["analyses"]

    def test_store_empty_result_skipped(self, temp_dir: Path) -> None:
        store = ResultStore(
            input_file=Path("test.wfn"),
            work_dir=temp_dir,
        )
        result = MultiwfnResult(analysis=Menu.HIRSHFELD_CHARGE)
        store.store(result)
        assert store.has_result(Menu.HIRSHFELD_CHARGE) is False

    def test_store_from_stdout(self, temp_dir: Path) -> None:
        store = ResultStore(
            input_file=Path("test.wfn"),
            work_dir=temp_dir,
        )
        stdout = "Final atomic charges:\nAtom    1(C ):    -0.05230000\n"
        result = store.store_from_stdout(Menu.HIRSHFELD_CHARGE, stdout)
        assert result is not None
        assert store.has_result(Menu.HIRSHFELD_CHARGE) is True

    def test_store_from_stdout_empty(self, temp_dir: Path) -> None:
        store = ResultStore(
            input_file=Path("test.wfn"),
            work_dir=temp_dir,
        )
        result = store.store_from_stdout(Menu.HIRSHFELD_CHARGE, "")
        assert result is None

    def test_json_path_setter_flushes(self, temp_dir: Path) -> None:
        store = ResultStore(
            input_file=Path("test.wfn"),
            work_dir=temp_dir,
        )
        result = MultiwfnResult(analysis=Menu.HIRSHFELD_CHARGE)
        result.result.append(Charge(atom_id=1, charge=-0.05))
        store.store(result)

        jp = temp_dir / "new_results.json"
        store.json_path = jp
        assert jp.exists()

    def test_get_result(self, temp_dir: Path) -> None:
        store = ResultStore(
            input_file=Path("test.wfn"),
            work_dir=temp_dir,
        )
        result = MultiwfnResult(analysis=Menu.HIRSHFELD_CHARGE)
        result.result.append(Charge(atom_id=1, charge=-0.05))
        store.store(result)

        retrieved = store.get_result(Menu.HIRSHFELD_CHARGE)
        assert retrieved is not None
        assert "parsed" in retrieved

    def test_get_result_missing(self, temp_dir: Path) -> None:
        store = ResultStore(
            input_file=Path("test.wfn"),
            work_dir=temp_dir,
        )
        assert store.get_result(Menu.HIRSHFELD_CHARGE) is None

    def test_data_property(self, temp_dir: Path) -> None:
        store = ResultStore(
            input_file=Path("test.wfn"),
            work_dir=temp_dir,
        )
        data = store.data
        assert "input_file" in data
        assert "analyses" in data

    def test_load_existing_json(self, temp_dir: Path) -> None:
        jp = temp_dir / "existing.json"
        jp.write_text(
            json.dumps(
                {
                    "input_file": "test.wfn",
                    "analyses": {
                        "HIRSHFELD_CHARGE": {
                            "parsed": {
                                "analysis": "HIRSHFELD_CHARGE",
                                "results": [],
                            },
                            "timestamp": "2025-01-01T00:00:00",
                        }
                    },
                }
            ),
            encoding="utf-8",
        )

        store = ResultStore(
            input_file=Path("test.wfn"),
            work_dir=temp_dir,
            json_path=jp,
        )
        assert store.has_result(Menu.HIRSHFELD_CHARGE) is True
