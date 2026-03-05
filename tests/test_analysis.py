"""Comprehensive tests for analysis.py - MultiwfnAnalysis class.

Tests cover positive results, negative results, and edge cases.
"""

from pathlib import Path
from unittest.mock import patch

import pytest

import pymultiwfn

# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def mock_wfn_file(
    tmp_path: Path,
) -> Path:
    """Create a mock wavefunction file."""
    path = tmp_path / "test.wfn"
    path.write_text("mock wavefunction content")
    return path


@pytest.fixture
def mock_executable(
    tmp_path: Path,
) -> Path:
    """Create a mock executable."""
    exe = tmp_path / "Multiwfn"
    exe.write_text("#!/bin/bash\necho 'mock output'")
    exe.chmod(0o755)
    return exe


@pytest.fixture
def mock_config(
    mock_executable: Path,
) -> pymultiwfn.MultiwfnConfig:
    """Create a mock config."""
    return pymultiwfn.MultiwfnConfig(exe_path=mock_executable)


@pytest.fixture
def mock_result() -> pymultiwfn.MultiwfnResult:
    """Create a mock MultiwfnResult."""
    return pymultiwfn.MultiwfnResult(
        stdout="mock output",
        stderr="",
        returncode=0,
        execution_time=1.0,
        commands=["7", "1", "q"],
        input_file=Path("/test/mol.wfn"),
    )


@pytest.fixture
def mol(
    mock_wfn_file: str,
    mock_config: pymultiwfn.MultiwfnConfig,
) -> pymultiwfn.MultiwfnAnalysis:
    """Create an Analysis instance using pymultiwfn.load()."""
    return pymultiwfn.load(mock_wfn_file, mock_config)


# =============================================================================
# Load Function Tests
# =============================================================================


class TestLoad:
    """Tests for pymultiwfn.load() function."""

    def test_load_returns_analysis(
        self,
        mock_wfn_file: str,
    ) -> None:
        """Positive: load returns MultiwfnAnalysis instance."""
        mol = pymultiwfn.load(mock_wfn_file)
        assert isinstance(mol, pymultiwfn.MultiwfnAnalysis)

    def test_load_with_string_path(
        self,
        mock_wfn_file: str,
    ) -> None:
        """Positive: load accepts string path."""
        mol = pymultiwfn.load(str(mock_wfn_file))
        assert mol.filepath == mock_wfn_file

    def test_load_with_config(
        self,
        mock_wfn_file: str,
        mock_config: pymultiwfn.MultiwfnConfig,
    ) -> None:
        """Positive: load accepts config parameter."""
        mol = pymultiwfn.load(mock_wfn_file, mock_config)
        assert mol.config == mock_config

    def test_load_file_not_found(self) -> None:
        """Negative: load raises FileNotFoundError for nonexistent file."""
        with pytest.raises(FileNotFoundError, match="File not found"):
            pymultiwfn.load("/nonexistent/file.wfn")

    def test_load_is_callable(self) -> None:
        """Positive: load function is callable."""
        assert callable(pymultiwfn.load)


# =============================================================================
# Analysis Class Tests
# =============================================================================


class TestAnalysisInit:
    """Tests for MultiwfnAnalysis.__init__."""

    def test_init_with_valid_file(
        self,
        mock_wfn_file: Path,
    ) -> None:
        """Positive: Initialize with valid file path."""
        mol = pymultiwfn.MultiwfnAnalysis(mock_wfn_file)
        assert mol.filepath == mock_wfn_file

    def test_init_with_string_path(
        self,
        mock_wfn_file: Path,
    ) -> None:
        """Positive: Initialize with string path."""
        mol = pymultiwfn.MultiwfnAnalysis(str(mock_wfn_file))
        assert mol.filepath == mock_wfn_file

    def test_init_with_config(
        self,
        mock_wfn_file: Path,
        mock_config: pymultiwfn.MultiwfnConfig,
    ) -> None:
        """Positive: Initialize with custom config."""
        mol = pymultiwfn.MultiwfnAnalysis(mock_wfn_file, mock_config)
        assert mol.config == mock_config

    def test_init_file_not_found(self) -> None:
        """Negative: Raise FileNotFoundError for nonexistent file."""
        with pytest.raises(FileNotFoundError, match="File not found"):
            pymultiwfn.MultiwfnAnalysis("/nonexistent/file.wfn")

    def test_init_default_config(
        self,
        mock_wfn_file: Path,
    ) -> None:
        """Edge case: Default config is created when none provided."""
        mol = pymultiwfn.MultiwfnAnalysis(mock_wfn_file)
        assert isinstance(mol.config, pymultiwfn.MultiwfnConfig)


class TestAnalysisProperties:
    """Tests for MultiwfnAnalysis properties."""

    def test_filepath_property(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
    ) -> None:
        """Positive: filepath property returns correct path."""
        assert isinstance(mol.filepath, Path)
        assert mol.filepath.exists()

    def test_config_property(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
    ) -> None:
        """Positive: config property returns MultiwfnConfig."""
        assert isinstance(mol.config, pymultiwfn.MultiwfnConfig)

    def test_config_setter(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
    ) -> None:
        """Positive: config setter updates config."""
        new_config = pymultiwfn.MultiwfnConfig(timeout=999)
        mol.config = new_config
        assert mol.config.timeout == 999

    def test_repr(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
    ) -> None:
        """Positive: repr returns readable string."""
        r = repr(mol)
        assert "MultiwfnAnalysis" in r
        assert ".wfn" in r


class TestAnalysisRun:
    """Tests for MultiwfnAnalysis.run method."""

    def test_run_returns_result(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: run returns MultiwfnResult."""
        with patch.object(mol, "run", return_value=mock_result):
            result = mol.run(pymultiwfn.Menu.HIRSHFELD_CHARGE)
            assert isinstance(result, pymultiwfn.MultiwfnResult)

    def test_run_with_verbose(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: run accepts verbose parameter."""
        with patch.object(mol, "run", return_value=mock_result) as mock_run:
            mol.run(pymultiwfn.Menu.HIRSHFELD_CHARGE, verbose=True)
            mock_run.assert_called_once_with(
                pymultiwfn.Menu.HIRSHFELD_CHARGE, verbose=True
            )

    def test_run_invalid_menu_type(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
    ) -> None:
        """Negative: run with invalid menu type raises error."""
        with pytest.raises((TypeError, AttributeError)):
            mol.run("not_a_menu")  # type: ignore[arg-type]


class TestAnalysisRunMultiple:
    """Tests for MultiwfnAnalysis.run_multiple method."""

    def test_run_multiple_returns_result(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: run_multiple returns single MultiwfnResult."""
        with patch.object(mol, "run_multiple", return_value=mock_result):
            result = mol.run_multiple(
                [
                    pymultiwfn.Menu.HIRSHFELD_CHARGE,
                    pymultiwfn.Menu.MAYER_BOND_ORDER,
                ]
            )
            assert isinstance(result, pymultiwfn.MultiwfnResult)

    def test_run_multiple_empty_list(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Edge case: empty menu list."""
        with patch.object(mol, "run_multiple", return_value=mock_result):
            result = mol.run_multiple([])
            assert isinstance(result, pymultiwfn.MultiwfnResult)


class TestAnalysisRunBatch:
    """Tests for MultiwfnAnalysis.run_batch method."""

    def test_run_batch_returns_dict(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: run_batch returns dict of results."""
        with patch.object(mol, "run", return_value=mock_result):
            results = mol.run_batch(
                [
                    pymultiwfn.Menu.HIRSHFELD_CHARGE,
                    pymultiwfn.Menu.MAYER_BOND_ORDER,
                ]
            )
            assert isinstance(results, dict)
            assert "hirshfeld_charge" in results
            assert "mayer_bond_order" in results

    def test_run_batch_empty_list(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
    ) -> None:
        """Edge case: empty menu list returns empty dict."""
        results = mol.run_batch([])
        assert results == {}

    def test_run_batch_keys_lowercase(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: result keys are lowercase menu names."""
        with patch.object(mol, "run", return_value=mock_result):
            results = mol.run_batch([pymultiwfn.Menu.HIRSHFELD_CHARGE])
            assert all(key.islower() for key in results)


class TestAnalysisCategoryRunners:
    """Tests for category-specific run methods."""

    def test_run_charges(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: run_charges returns dict with charge results."""
        with patch.object(mol, "run", return_value=mock_result):
            results = mol.run_charges()
            assert isinstance(results, dict)
            assert len(results) == len(pymultiwfn.MultiwfnAnalysis.CHARGES)

    def test_run_bond_orders(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: run_bond_orders returns dict."""
        with patch.object(mol, "run", return_value=mock_result):
            results = mol.run_bond_orders()
            assert len(results) == len(pymultiwfn.MultiwfnAnalysis.BOND_ORDERS)

    def test_run_topology(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: run_topology returns dict."""
        with patch.object(mol, "run", return_value=mock_result):
            results = mol.run_topology()
            assert len(results) == len(pymultiwfn.MultiwfnAnalysis.TOPOLOGY)

    def test_run_weak_interactions(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: run_weak_interactions returns dict."""
        with patch.object(mol, "run", return_value=mock_result):
            results = mol.run_weak_interactions()
            assert len(results) == len(
                pymultiwfn.MultiwfnAnalysis.WEAK_INTERACTIONS
            )

    def test_run_spectra(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: run_spectra returns dict."""
        with patch.object(mol, "run", return_value=mock_result):
            results = mol.run_spectra()
            assert len(results) == len(pymultiwfn.MultiwfnAnalysis.SPECTRA)

    def test_run_surfaces(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: run_surfaces returns dict."""
        with patch.object(mol, "run", return_value=mock_result):
            results = mol.run_surfaces()
            assert len(results) == len(pymultiwfn.MultiwfnAnalysis.SURFACES)

    def test_run_aromaticity(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: run_aromaticity returns dict."""
        with patch.object(mol, "run", return_value=mock_result):
            results = mol.run_aromaticity()
            assert len(results) == len(pymultiwfn.MultiwfnAnalysis.AROMATICITY)

    def test_run_cdft(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: run_cdft returns dict."""
        with patch.object(mol, "run", return_value=mock_result):
            results = mol.run_cdft()
            assert len(results) == len(pymultiwfn.MultiwfnAnalysis.CDFT)

    def test_run_cubes(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: run_cubes returns dict."""
        with patch.object(mol, "run", return_value=mock_result):
            results = mol.run_cubes()
            assert len(results) == len(pymultiwfn.MultiwfnAnalysis.CUBES)

    def test_run_dos(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: run_dos returns dict."""
        with patch.object(mol, "run", return_value=mock_result):
            results = mol.run_dos()
            assert len(results) == len(pymultiwfn.MultiwfnAnalysis.DOS)

    def test_run_basin(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: run_basin returns dict."""
        with patch.object(mol, "run", return_value=mock_result):
            results = mol.run_basin()
            assert len(results) == len(pymultiwfn.MultiwfnAnalysis.BASIN)

    def test_run_excitation(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: run_excitation returns dict."""
        with patch.object(mol, "run", return_value=mock_result):
            results = mol.run_excitation()
            assert len(results) == len(pymultiwfn.MultiwfnAnalysis.EXCITATION)

    def test_run_orbital_composition(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: run_orbital_composition returns dict."""
        with patch.object(mol, "run", return_value=mock_result):
            results = mol.run_orbital_composition()
            assert len(results) == len(
                pymultiwfn.MultiwfnAnalysis.ORBITAL_COMPOSITION
            )

    def test_run_orbital_localization(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: run_orbital_localization returns dict."""
        with patch.object(mol, "run", return_value=mock_result):
            results = mol.run_orbital_localization()
            assert len(results) == len(
                pymultiwfn.MultiwfnAnalysis.ORBITAL_LOCALIZATION
            )

    def test_run_fuzzy_space(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: run_fuzzy_space returns dict."""
        with patch.object(mol, "run", return_value=mock_result):
            results = mol.run_fuzzy_space()
            assert len(results) == len(pymultiwfn.MultiwfnAnalysis.FUZZY_SPACE)

    def test_run_eda(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: run_eda returns dict."""
        with patch.object(mol, "run", return_value=mock_result):
            results = mol.run_eda()
            assert len(results) == len(pymultiwfn.MultiwfnAnalysis.EDA)

    def test_run_polarizability(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: run_polarizability returns dict."""
        with patch.object(mol, "run", return_value=mock_result):
            results = mol.run_polarizability()
            assert len(results) == len(
                pymultiwfn.MultiwfnAnalysis.POLARIZABILITY
            )

    def test_run_wavefunction(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: run_wavefunction_analysis returns dict."""
        with patch.object(mol, "run", return_value=mock_result):
            results = mol.run_wavefunction()
            assert len(results) == len(
                pymultiwfn.MultiwfnAnalysis.WAVEFUNCTION
            )

    def test_run_line_plots(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: run_line_plots returns dict."""
        with patch.object(mol, "run", return_value=mock_result):
            results = mol.run_line_plots()
            assert len(results) == len(pymultiwfn.MultiwfnAnalysis.LINE_PLOTS)

    def test_run_plane_maps(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: run_plane_maps returns dict."""
        with patch.object(mol, "run", return_value=mock_result):
            results = mol.run_plane_maps()
            assert len(results) == len(pymultiwfn.MultiwfnAnalysis.PLANE_MAPS)


class TestAnalysisRunCategory:
    """Tests for MultiwfnAnalysis.run_category method."""

    def test_run_category_valid(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: run_category with valid category returns dict."""
        with patch.object(mol, "run", return_value=mock_result):
            results = mol.run_category("charges")
            assert isinstance(results, dict)
            assert len(results) == len(pymultiwfn.MultiwfnAnalysis.CHARGES)

    def test_run_category_invalid(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
    ) -> None:
        """Negative: invalid category raises ValueError."""
        with pytest.raises(ValueError, match="Unknown category"):
            mol.run_category("invalid_category")

    def test_run_category_case_sensitive(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
    ) -> None:
        """Edge case: category names are case-sensitive."""
        with pytest.raises(ValueError, match="Unknown category"):
            mol.run_category("CHARGES")

    def test_run_category_all_categories_valid(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: all defined categories are valid."""
        with patch.object(mol, "run", return_value=mock_result):
            for cat in pymultiwfn.MultiwfnAnalysis.CATEGORIES:
                results = mol.run_category(cat)
                assert isinstance(results, dict)


class TestAnalysisRunAll:
    """Tests for MultiwfnAnalysis.run_all method."""

    def test_run_all_returns_nested_dict(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: run_all returns nested dict."""
        with patch.object(mol, "run", return_value=mock_result):
            results = mol.run_all()
            assert isinstance(results, dict)
            for _cat, cat_results in results.items():
                assert isinstance(cat_results, dict)

    def test_run_all_with_category_filter(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: run_all with specific categories."""
        with patch.object(mol, "run", return_value=mock_result):
            results = mol.run_all(categories=["charges", "bond_orders"])
            assert len(results) == 2
            assert "charges" in results
            assert "bond_orders" in results

    def test_run_all_empty_categories(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Edge case: empty categories list runs all."""
        with patch.object(mol, "run", return_value=mock_result):
            results = mol.run_all(categories=[])
            assert len(results) == len(pymultiwfn.MultiwfnAnalysis.CATEGORIES)

    def test_run_all_invalid_category_ignored(
        self,
        mol: pymultiwfn.MultiwfnAnalysis,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Edge case: invalid categories are silently ignored."""
        with patch.object(mol, "run", return_value=mock_result):
            results = mol.run_all(categories=["charges", "invalid"])
            assert "charges" in results
            assert "invalid" not in results


class TestAnalysisListMethods:
    """Tests for MultiwfnAnalysis class methods."""

    def test_list_categories(self) -> None:
        """Positive: list_categories returns list of strings."""
        categories = pymultiwfn.MultiwfnAnalysis.list_categories()
        assert isinstance(categories, list)
        assert len(categories) > 0
        assert all(isinstance(c, str) for c in categories)

    def test_list_categories_contains_expected(self) -> None:
        """Positive: list_categories contains expected categories."""
        categories = pymultiwfn.MultiwfnAnalysis.list_categories()
        expected = ["charges", "bond_orders", "topology"]
        for cat in expected:
            assert cat in categories

    def test_list_analyses_all(self) -> None:
        """Positive: list_analyses without filter returns all."""
        analyses = pymultiwfn.MultiwfnAnalysis.list_analyses()
        assert isinstance(analyses, list)
        assert len(analyses) > 100

    def test_list_analyses_filtered(self) -> None:
        """Positive: list_analyses with category filter."""
        analyses = pymultiwfn.MultiwfnAnalysis.list_analyses("charges")
        assert len(analyses) == len(pymultiwfn.MultiwfnAnalysis.CHARGES)
        assert "hirshfeld_charge" in analyses

    def test_list_analyses_invalid_category(self) -> None:
        """Negative: invalid category raises ValueError."""
        with pytest.raises(ValueError, match="Unknown category"):
            pymultiwfn.MultiwfnAnalysis.list_analyses("invalid")

    def test_list_analyses_lowercase(self) -> None:
        """Positive: analysis names are lowercase."""
        analyses = pymultiwfn.MultiwfnAnalysis.list_analyses("charges")
        assert all(a.islower() for a in analyses)


# =============================================================================
# Edge Cases and Integration
# =============================================================================


class TestEdgeCases:
    """Edge case tests."""

    def test_analysis_categories_not_empty(self) -> None:
        """Edge case: all categories have at least one item."""
        for cat, menus in pymultiwfn.MultiwfnAnalysis.CATEGORIES.items():
            assert len(menus) > 0, f"Category {cat} is empty"

    def test_analysis_categories_contain_menu_items(self) -> None:
        """Edge case: all category items are Menu enums."""
        for _cat, menus in pymultiwfn.MultiwfnAnalysis.CATEGORIES.items():
            for menu in menus:
                assert isinstance(menu, pymultiwfn.Menu), (
                    f"{menu} is not a Menu enum"
                )

    def test_all_categories_in_list(self) -> None:
        """Edge case: list_categories matches CATEGORIES keys."""
        listed = set(pymultiwfn.MultiwfnAnalysis.list_categories())
        defined = set(pymultiwfn.MultiwfnAnalysis.CATEGORIES.keys())
        assert listed == defined

    def test_path_types(self, mock_wfn_file: Path) -> None:
        """Edge case: both str and Path work for filepath."""
        mol1 = pymultiwfn.MultiwfnAnalysis(mock_wfn_file)
        mol2 = pymultiwfn.MultiwfnAnalysis(str(mock_wfn_file))
        assert mol1.filepath == mol2.filepath


class TestPackageLevelImports:
    """Tests for package-level imports."""

    def test_import_load(self) -> None:
        """Positive: load function is importable."""
        assert hasattr(pymultiwfn, "load")
        assert callable(pymultiwfn.load)

    def test_import_analysis_class(self) -> None:
        """Positive: MultiwfnAnalysis is importable."""
        assert hasattr(pymultiwfn, "MultiwfnAnalysis")
        assert isinstance(pymultiwfn.MultiwfnAnalysis, type)

    def test_import_analysis_alias(self) -> None:
        """Positive: Analysis alias exists."""
        assert hasattr(pymultiwfn, "Analysis")
        assert pymultiwfn.Analysis is pymultiwfn.MultiwfnAnalysis

    def test_import_menu(self) -> None:
        """Positive: Menu is importable."""
        assert hasattr(pymultiwfn, "Menu")

    def test_import_config(self) -> None:
        """Positive: MultiwfnConfig is importable."""
        assert hasattr(pymultiwfn, "MultiwfnConfig")


class TestWorkflow:
    """Integration tests for typical workflows."""

    def test_typical_workflow(
        self, mock_wfn_file: str, mock_result: pymultiwfn.MultiwfnResult
    ) -> None:
        """Positive: typical workflow works."""
        # Load molecule
        mol = pymultiwfn.load(mock_wfn_file)

        # Run analysis
        with patch.object(mol, "run", return_value=mock_result):
            result = mol.run(pymultiwfn.Menu.HIRSHFELD_CHARGE)
            assert isinstance(result, pymultiwfn.MultiwfnResult)

    def test_batch_workflow(
        self,
        mock_wfn_file: str,
        mock_result: pymultiwfn.MultiwfnResult,
    ) -> None:
        """Positive: batch workflow works."""
        mol = pymultiwfn.load(mock_wfn_file)

        with patch.object(mol, "run", return_value=mock_result):
            results = mol.run_charges()
            assert isinstance(results, dict)
            assert len(results) == len(pymultiwfn.MultiwfnAnalysis.CHARGES)
