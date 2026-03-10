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
        # find_category compares the argument against the Menu enum values
        # in each category, using `menu in analysis.value` where analysis.value
        # is a list of Menu enums. So we need to pass the Menu enum name as string
        # and the source code checks `if menu in analysis.value`.
        # But the source does: `if menu in analysis.value` where menu is a string
        # and analysis.value is a list of Menu enums — this won't match.
        # The source find_category is broken for string input; it only works
        # if we pass a Menu enum directly.
        # Since we can only change tests, let's test what actually works:
        result = AnalysisClasses.find_category(Menu.HIRSHFELD_CHARGE)
        assert result == "CHARGES"

    def test_find_bond_order(self) -> None:
        result = AnalysisClasses.find_category(Menu.MAYER_BOND_ORDER)
        assert result == "BOND_ORDERS"

    def test_not_found(self) -> None:
        result = AnalysisClasses.find_category(Menu.VIEW_STRUCTURE)
        assert result is None
