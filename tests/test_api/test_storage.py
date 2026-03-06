"""Tests for pymultiwfn.api.storage — ResultStore and parser routing."""

import json
from pathlib import Path

import pytest

from pymultiwfn.api.storage import ResultStore, _PARSER_ROUTE
from pymultiwfn.enums.menu import Menu


@pytest.fixture
def store(temp_dir: Path, mock_wfn_file: Path) -> ResultStore:
    return ResultStore(input_file=mock_wfn_file, work_dir=temp_dir)


class TestResultStoreInit:
    def test_json_path_naming(self, store: ResultStore) -> None:
        assert store.json_path.name == "mock.wfn.json"

    def test_empty_data_on_fresh(self, store: ResultStore) -> None:
        data = store.data
        assert "input_file" in data
        assert data["analyses"] == {}

    def test_creates_work_dir(self, temp_dir: Path, mock_wfn_file: Path) -> None:
        sub = temp_dir / "deep" / "nested"
        store = ResultStore(input_file=mock_wfn_file, work_dir=sub)
        assert sub.exists()


class TestStoreAndRetrieve:
    def test_store_charges(self, store: ResultStore) -> None:
        stdout = (
            "Final atomic charges:\n"
            "Atom    1(C ):    -0.05230000\n"
            "Atom    2(H ):     0.05230000\n"
        )
        parsed = store.store_result(Menu.HIRSHFELD_CHARGE, stdout)

        assert parsed is not None
        assert "charges" in parsed
        assert parsed["charges"]["1"] == pytest.approx(-0.0523)

    def test_has_result_after_store(self, store: ResultStore) -> None:
        stdout = "Final atomic charges:\nAtom    1(C ):    -0.05\n"
        store.store_result(Menu.HIRSHFELD_CHARGE, stdout)
        assert store.has_result(Menu.HIRSHFELD_CHARGE) is True

    def test_has_result_false_before_store(self, store: ResultStore) -> None:
        assert store.has_result(Menu.HIRSHFELD_CHARGE) is False

    def test_get_result(self, store: ResultStore) -> None:
        stdout = "Final atomic charges:\nAtom    1(C ):    -0.05\n"
        store.store_result(Menu.HIRSHFELD_CHARGE, stdout)
        result = store.get_result(Menu.HIRSHFELD_CHARGE)
        assert result is not None
        assert "parsed" in result
        assert "timestamp" in result

    def test_get_result_none_for_missing(self, store: ResultStore) -> None:
        assert store.get_result(Menu.HIRSHFELD_CHARGE) is None


class TestPersistence:
    def test_json_written_to_disk(self, store: ResultStore) -> None:
        stdout = "Final atomic charges:\nAtom    1(C ):    -0.05\n"
        store.store_result(Menu.HIRSHFELD_CHARGE, stdout)
        assert store.json_path.exists()

        with open(store.json_path) as f:
            data = json.load(f)
        assert "HIRSHFELD_CHARGE" in data["analyses"]

    def test_reload_existing_json(
        self, temp_dir: Path, mock_wfn_file: Path
    ) -> None:
        store1 = ResultStore(input_file=mock_wfn_file, work_dir=temp_dir)
        stdout = "Final atomic charges:\nAtom    1(C ):    -0.05\n"
        store1.store_result(Menu.HIRSHFELD_CHARGE, stdout)

        store2 = ResultStore(input_file=mock_wfn_file, work_dir=temp_dir)
        assert store2.has_result(Menu.HIRSHFELD_CHARGE) is True

    def test_multiple_analyses_accumulate(self, store: ResultStore) -> None:
        store.store_result(
            Menu.HIRSHFELD_CHARGE,
            "Final atomic charges:\nAtom    1(C ):    -0.05\n",
        )
        store.store_result(
            Menu.MAYER_BOND_ORDER,
            "#    1:         1(C )    2(H )    0.95000000\n",
        )
        data = store.data
        assert "HIRSHFELD_CHARGE" in data["analyses"]
        assert "MAYER_BOND_ORDER" in data["analyses"]


class TestParserRouting:
    def test_no_parser_returns_none(self, store: ResultStore) -> None:
        """Menu items without a parser route return None."""
        result = store.store_result(Menu.VIEW_STRUCTURE, "some output")
        assert result is None

    def test_bond_order_routing(self, store: ResultStore) -> None:
        stdout = "#    1:         1(C )    2(N )    1.45230000\n"
        parsed = store.store_result(Menu.MAYER_BOND_ORDER, stdout)
        assert parsed is not None
        assert "bond_orders" in parsed

    def test_topology_routing(self, store: ResultStore) -> None:
        stdout = (
            "CP  1 (3,-3)\n"
            "Position (Bohr):   0.000000   0.000000   0.000000\n"
            "Density of all electrons:   0.29876500\n"
        )
        parsed = store.store_result(Menu.TOPOLOGY_SEARCH_CPS, stdout)
        assert parsed is not None
        assert "critical_points" in parsed

    def test_spectrum_routing(self, store: ResultStore) -> None:
        stdout = "  500.0 cm^-1  Intensity:  12.34\n"
        parsed = store.store_result(Menu.PLOT_IR_SPECTRUM, stdout)
        assert parsed is not None
        assert "peaks" in parsed

    def test_surface_routing(self, store: ResultStore) -> None:
        stdout = "Overall surface area:  234.5678\nEnclosed volume:  345.6789\n"
        parsed = store.store_result(Menu.SURFACE_ANALYSIS_ESP, stdout)
        assert parsed is not None
        assert "surface_statistics" in parsed

    def test_empty_stdout_returns_none(self, store: ResultStore) -> None:
        """Empty stdout produces empty parse → None."""
        parsed = store.store_result(Menu.HIRSHFELD_CHARGE, "")
        assert parsed is None

    def test_all_charge_menus_routed(self) -> None:
        charge_menus = [
            Menu.HIRSHFELD_CHARGE, Menu.VDD_POPULATION,
            Menu.MULLIKEN_POPULATION, Menu.LOWDIN_POPULATION,
            Menu.ADCH_CHARGE, Menu.CHELPG_CHARGE, Menu.MK_CHARGE,
            Menu.RESP_CHARGE, Menu.CM5_CHARGE, Menu.MBIS_CHARGE,
            Menu.DDEC_CHARGE,
        ]
        for menu in charge_menus:
            assert menu in _PARSER_ROUTE, f"{menu.name} missing from route"

    def test_all_bond_order_menus_routed(self) -> None:
        bo_menus = [
            Menu.MAYER_BOND_ORDER, Menu.WIBERG_BOND_ORDER,
            Menu.MULLIKEN_BOND_ORDER, Menu.FUZZY_BOND_ORDER,
            Menu.LAPLACIAN_BOND_ORDER,
        ]
        for menu in bo_menus:
            assert menu in _PARSER_ROUTE, f"{menu.name} missing from route"