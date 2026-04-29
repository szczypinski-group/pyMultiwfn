"""Tests for pymultiwfn.enums.analyses — AnalysisClasses."""

import pytest

from pymultiwfn.enums.analyses import AnalysisClasses
from pymultiwfn.enums.menu import Menu


class TestAnalysisClassesStructure:
    """Tests for the structure and contents of AnalysisClasses."""

    def test_all_categories_non_empty(self) -> None:
        for category in AnalysisClasses:
            assert len(category.value) > 0, f"{category.name} is empty"

    def test_all_values_are_menu_items(self) -> None:
        for category in AnalysisClasses:
            for item in category.value:
                assert isinstance(item, Menu), (
                    f"{item} in {category.name} is not a Menu"
                )

    def test_charges_contains_hirshfeld(self) -> None:
        assert Menu.HIRSHFELD_CHARGE in AnalysisClasses.CHARGES.value

    def test_bond_orders_contains_mayer(self) -> None:
        assert Menu.MAYER_BOND_ORDER in AnalysisClasses.BOND_ORDERS.value

    def test_every_menu_member_is_in_at_least_one_category(self) -> None:
        """Every Menu enum member should appear in at least one category."""
        all_categorised: set[Menu] = set()
        for category in AnalysisClasses:
            for item in category.value:
                all_categorised.add(item)
        missing = [m for m in Menu if m not in all_categorised]
        assert missing == [], (
            f"Menu members not in any AnalysisClasses category: "
            f"{[m.name for m in missing]}"
        )

    def test_no_duplicate_within_category(self) -> None:
        """No Menu member should appear twice in the same category."""
        for category in AnalysisClasses:
            seen: set[Menu] = set()
            for item in category.value:
                assert item not in seen, (
                    f"{item.name} appears twice in {category.name}"
                )
                seen.add(item)


class TestListCategories:
    """Tests for AnalysisClasses.list_categories()."""

    def test_returns_list(self) -> None:
        cats = AnalysisClasses.list_categories()
        assert isinstance(cats, list)
        assert len(cats) > 0

    def test_contains_expected(self) -> None:
        cats = AnalysisClasses.list_categories()
        assert "CHARGES" in cats
        assert "BOND_ORDERS" in cats
        assert "TOPOLOGY" in cats

    def test_contains_new_categories(self) -> None:
        cats = AnalysisClasses.list_categories()
        assert "GRID_PROCESSING" in cats
        assert "UTILITIES" in cats
        assert "ORBITAL_ANALYSIS" in cats
        assert "SPATIAL_DELOCALIZATION" in cats
        assert "FILE_EXPORT" in cats


class TestListAnalyses:
    """Tests for AnalysisClasses.list_analyses()."""

    def test_all_analyses(self) -> None:
        analyses = AnalysisClasses.list_analyses()
        assert len(analyses) > 50

    def test_filtered_by_category(self) -> None:
        analyses = AnalysisClasses.list_analyses("CHARGES")
        assert len(analyses) == len(AnalysisClasses.CHARGES.value)

    def test_invalid_category(self) -> None:
        with pytest.raises(ValueError, match="Unknown category"):
            AnalysisClasses.list_analyses("NONEXISTENT")


class TestFindCategory:
    """Tests for AnalysisClasses.find_category()."""

    def test_find_existing(self) -> None:
        result = AnalysisClasses.find_category(Menu.HIRSHFELD_CHARGE)
        assert result == "CHARGES"

    def test_find_bond_order(self) -> None:
        result = AnalysisClasses.find_category(Menu.MAYER_BOND_ORDER)
        assert result == "BOND_ORDERS"

    def test_find_new_categories(self) -> None:
        assert (
            AnalysisClasses.find_category(Menu.GRID_EXTRACT_PLANE_XY)
            == "GRID_PROCESSING"
        )
        assert (
            AnalysisClasses.find_category(Menu.GEOMETRY_PROPERTIES)
            == "UTILITIES"
        )
        assert (
            AnalysisClasses.find_category(Menu.ORBITAL_OVERLAP_INTEGRAL)
            == "ORBITAL_ANALYSIS"
        )
        assert (
            AnalysisClasses.find_category(Menu.SPACIAL_DELOCALISATION_EDENSITY)
            == "SPATIAL_DELOCALIZATION"
        )
        assert (
            AnalysisClasses.find_category(Menu.EXPORT_VARIOUS_FILES)
            == "FILE_EXPORT"
        )
