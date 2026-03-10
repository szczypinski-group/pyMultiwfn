"""Tests for pymultiwfn.analysis.analysis — MultiwfnAnalysis."""

from pathlib import Path

import pytest

from pymultiwfn.analysis.analysis import MultiwfnAnalysis
from pymultiwfn.enums.analyses import AnalysisClasses
from pymultiwfn.enums.menu import Menu


class TestMultiwfnAnalysisInit:
    """Tests for initializing MultiwfnAnalysis."""

    def test_single_file(self, mock_wfn_file: Path) -> None:
        mol = MultiwfnAnalysis(mock_wfn_file)
        assert mol.input_file == mock_wfn_file

    def test_single_file_string(self, mock_wfn_file: Path) -> None:
        mol = MultiwfnAnalysis(str(mock_wfn_file))
        assert mol.input_file == Path(str(mock_wfn_file))

    def test_with_single_menu(self, mock_wfn_file: Path) -> None:
        mol = MultiwfnAnalysis(mock_wfn_file, analyses=Menu.HIRSHFELD_CHARGE)
        assert mol.analyses == [Menu.HIRSHFELD_CHARGE]

    def test_with_menu_list(self, mock_wfn_file: Path) -> None:
        menus = [Menu.HIRSHFELD_CHARGE, Menu.MAYER_BOND_ORDER]
        mol = MultiwfnAnalysis(mock_wfn_file, analyses=menus)
        assert mol.analyses == menus

    def test_with_analysis_class(self, mock_wfn_file: Path) -> None:
        mol = MultiwfnAnalysis(mock_wfn_file, analyses=AnalysisClasses.CHARGES)
        assert len(mol.analyses) == len(AnalysisClasses.CHARGES.value)

    def test_no_analyses(self, mock_wfn_file: Path) -> None:
        mol = MultiwfnAnalysis(mock_wfn_file)
        assert mol.analyses == []

    def test_cached_default_true(self, mock_wfn_file: Path) -> None:
        mol = MultiwfnAnalysis(mock_wfn_file)
        assert mol.cached is True

    def test_cached_false(self, mock_wfn_file: Path) -> None:
        mol = MultiwfnAnalysis(mock_wfn_file, cached=False)
        assert mol.cached is False


class TestAddMenu:
    """Tests for the add_menu() method of MultiwfnAnalysis."""

    def test_add_single(self, mock_wfn_file: Path) -> None:
        mol = MultiwfnAnalysis(mock_wfn_file)
        mol.add_menu(Menu.HIRSHFELD_CHARGE)
        assert Menu.HIRSHFELD_CHARGE in mol.analyses

    def test_add_list(self, mock_wfn_file: Path) -> None:
        mol = MultiwfnAnalysis(mock_wfn_file)
        mol.add_menu([Menu.HIRSHFELD_CHARGE, Menu.MAYER_BOND_ORDER])
        assert len(mol.analyses) == 2

    def test_add_analysis_class(self, mock_wfn_file: Path) -> None:
        mol = MultiwfnAnalysis(mock_wfn_file)
        mol.add_menu(AnalysisClasses.BOND_ORDERS)
        assert len(mol.analyses) == len(AnalysisClasses.BOND_ORDERS.value)

    def test_add_invalid_type(self, mock_wfn_file: Path) -> None:
        mol = MultiwfnAnalysis(mock_wfn_file)
        with pytest.raises(TypeError, match="valid Menu enum"):
            mol.add_menu("not_a_menu")  # type: ignore[arg-type]

    def test_add_multiple_times(self, mock_wfn_file: Path) -> None:
        mol = MultiwfnAnalysis(mock_wfn_file)
        mol.add_menu(Menu.HIRSHFELD_CHARGE)
        mol.add_menu(Menu.MAYER_BOND_ORDER)
        assert len(mol.analyses) == 2


class TestConvenienceMethods:
    """Test the _add_*_menus convenience methods."""

    def test_add_charges(self, mock_wfn_file: Path) -> None:
        mol = MultiwfnAnalysis(mock_wfn_file)
        mol._add_charges_menus()
        assert len(mol.analyses) == len(AnalysisClasses.CHARGES.value)

    def test_add_bond_orders(self, mock_wfn_file: Path) -> None:
        mol = MultiwfnAnalysis(mock_wfn_file)
        mol._add_bond_orders_menus()
        assert len(mol.analyses) == len(AnalysisClasses.BOND_ORDERS.value)

    def test_add_topology(self, mock_wfn_file: Path) -> None:
        mol = MultiwfnAnalysis(mock_wfn_file)
        mol._add_topology_menus()
        assert len(mol.analyses) == len(AnalysisClasses.TOPOLOGY.value)

    def test_add_spectra(self, mock_wfn_file: Path) -> None:
        mol = MultiwfnAnalysis(mock_wfn_file)
        mol._add_spectra_menus()
        assert len(mol.analyses) == len(AnalysisClasses.SPECTRA.value)

    def test_add_surfaces(self, mock_wfn_file: Path) -> None:
        mol = MultiwfnAnalysis(mock_wfn_file)
        mol._add_surfaces_menus()
        assert len(mol.analyses) == len(AnalysisClasses.SURFACES.value)

    def test_add_aromaticity(self, mock_wfn_file: Path) -> None:
        mol = MultiwfnAnalysis(mock_wfn_file)
        mol._add_aromaticity_menus()
        assert len(mol.analyses) == len(AnalysisClasses.AROMATICITY.value)

    def test_add_cdft(self, mock_wfn_file: Path) -> None:
        mol = MultiwfnAnalysis(mock_wfn_file)
        mol._add_cdft_menus()
        assert len(mol.analyses) == len(AnalysisClasses.CDFT.value)

    def test_add_cubes(self, mock_wfn_file: Path) -> None:
        mol = MultiwfnAnalysis(mock_wfn_file)
        mol._add_cubes_menus()
        assert len(mol.analyses) == len(AnalysisClasses.CUBES.value)


class TestLogPath:
    """Tests for the _store attribute of MultiwfnAnalysis."""

    def test_store_none_before_run(self, mock_wfn_file: Path) -> None:
        mol = MultiwfnAnalysis(mock_wfn_file)
        assert mol._store is None

    def test_get_store_returns_store(self, mock_wfn_file: Path) -> None:
        mol = MultiwfnAnalysis(mock_wfn_file)
        store = mol._get_store()
        assert store is not None
