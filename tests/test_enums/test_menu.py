"""Tests for pymultiwfn.enums.menu — Menu enum and custom_sequence."""

import pytest

from pymultiwfn.enums.menu import Menu, custom_sequence


class TestMenuStructure:
    """Tests for the structure of the Menu enum."""

    def test_is_enum(self) -> None:
        assert hasattr(Menu, "__members__")
        assert len(Menu) > 100

    def test_all_values_are_tuples(self) -> None:
        for item in Menu:
            assert isinstance(item.value, tuple)
            assert all(isinstance(x, str) for x in item.value)


class TestGetSequence:
    """Tests for Menu.get_sequence()."""

    def test_returns_tuple(self) -> None:
        seq = Menu.HIRSHFELD_CHARGE.get_sequence()
        assert isinstance(seq, tuple)

    def test_hirshfeld_sequence(self) -> None:
        seq = Menu.HIRSHFELD_CHARGE.get_sequence()
        assert seq == ("7", "1", "1", "n", "0")

    def test_mayer_starts_with_9(self) -> None:
        assert Menu.MAYER_BOND_ORDER.get_sequence()[0] == "9"


class TestMenuCategories:
    """Tests for Menu category properties."""

    @pytest.mark.parametrize(
        "item",
        [
            Menu.HIRSHFELD_CHARGE,
            Menu.MULLIKEN_POPULATION,
            Menu.ADCH_CHARGE,
            Menu.RESP_CHARGE,
            Menu.CM5_CHARGE,
        ],
    )
    def test_charges_start_with_7(self, item: Menu) -> None:
        assert item.get_sequence()[0] == "7"

    @pytest.mark.parametrize(
        "item",
        [Menu.MAYER_BOND_ORDER, Menu.WIBERG_BOND_ORDER, Menu.FUZZY_BOND_ORDER],
    )
    def test_bond_orders_start_with_9(self, item: Menu) -> None:
        assert item.get_sequence()[0] == "9"

    @pytest.mark.parametrize(
        "item",
        [Menu.NCI_ANALYSIS, Menu.IRI_ANALYSIS, Menu.IGM_ANALYSIS],
    )
    def test_weak_start_with_20(self, item: Menu) -> None:
        assert item.get_sequence()[0] == "20"

    @pytest.mark.parametrize(
        "item",
        [Menu.CUBE_DENSITY, Menu.CUBE_ESP, Menu.CUBE_ELF],
    )
    def test_cubes_start_with_5(self, item: Menu) -> None:
        assert item.get_sequence()[0] == "5"

    @pytest.mark.parametrize(
        "item",
        [Menu.PLOT_IR_SPECTRUM, Menu.PLOT_UV_VIS_SPECTRUM],
    )
    def test_spectra_start_with_11(self, item: Menu) -> None:
        assert item.get_sequence()[0] == "11"


class TestSearch:
    """Tests for Menu.search()."""

    def test_search_returns_matches(self) -> None:
        results = Menu.search("hirshfeld")
        assert Menu.HIRSHFELD_CHARGE in results

    def test_search_case_insensitive(self) -> None:
        results = Menu.search("BoNd")
        assert any("BOND" in item.name for item in results)

    def test_search_no_match(self) -> None:
        assert Menu.search("zzzznonexistent") == []


class TestListAll:
    """Tests for Menu.list_all()."""

    def test_returns_names(self) -> None:
        names = Menu.list_all()
        assert "HIRSHFELD_CHARGE" in names
        assert "MAYER_BOND_ORDER" in names

    def test_count_matches_enum(self) -> None:
        assert len(Menu.list_all()) == len(Menu)


class TestCustomSequence:
    """Tests for the custom_sequence() function."""

    def test_appends_q(self) -> None:
        assert custom_sequence(["7", "1"]) == ["7", "1", "q"]

    def test_keeps_existing_q(self) -> None:
        seq = custom_sequence(["7", "1", "q"])
        assert seq.count("q") == 1

    def test_requires_list(self) -> None:
        with pytest.raises(TypeError, match="list"):
            custom_sequence("7")  # type: ignore[arg-type]

    def test_empty_list(self) -> None:
        assert custom_sequence([]) == ["q"]