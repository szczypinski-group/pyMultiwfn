"""Tests for Enum-based Menu implementation (pv-dev branch)."""

import pytest

from pymultiwfn.menu import Menu, custom_sequence


class TestEnumStructure:
    """Basic structural tests for Menu Enum."""

    def test_menu_is_enum(self) -> None:
        """Menu should behave as an Enum."""
        assert hasattr(Menu, "__members__")
        assert len(Menu) > 100  # Should contain many items

    def test_all_members_have_tuple_values(self) -> None:
        """Each Enum member should store a tuple of strings."""
        for item in Menu:
            assert isinstance(item.value, tuple)
            assert all(isinstance(x, str) for x in item.value)


class TestGetSequence:
    """Tests for get_sequence()."""

    def test_all_sequences_return_list(self) -> None:
        """All menu items should return list[str] ending with 'q'."""
        for item in Menu:
            seq = item.get_sequence()
            assert isinstance(seq, list)
            assert all(isinstance(x, str) for x in seq)
            assert seq[-1] == "q"

    def test_sequence_preserves_original_values(self) -> None:
        """get_sequence should preserve tuple values before 'q'."""
        item = Menu.HIRSHFELD_CHARGE
        seq = item.get_sequence()
        assert seq[:-1] == list(item.value)


class TestMenuCategories:
    """Tests that specific categories start with correct main menu numbers."""

    @pytest.mark.parametrize(
        "item",
        [
            Menu.HIRSHFELD_CHARGE,
            Menu.MULLIKEN_POPULATION,
            Menu.LOWDIN_POPULATION,
            Menu.ADCH_CHARGE,
            Menu.RESP_CHARGE,
            Menu.CM5_CHARGE,
        ],
    )
    def test_population_menu_starts_with_7(self, item: Menu) -> None:
        """
        Test to ensure population analysis menu items start with '7'.

        :param self: Description
        :param item: Description
        :type item: Menu
        """
        seq = item.get_sequence()
        assert seq[0] == "7"

    @pytest.mark.parametrize(
        "item",
        [
            Menu.MAYER_BOND_ORDER,
            Menu.WIBERG_BOND_ORDER,
            Menu.MULTICENTER_BOND_ORDER,
            Menu.FUZZY_BOND_ORDER,
        ],
    )
    def test_bond_order_menu_starts_with_9(
        self,
        item: Menu,
    ) -> None:
        """
        Test to ensure bond order analysis menu items start with '9'.

        :param self: Description
        :param item: Description
        :type item: Menu
        """
        seq = item.get_sequence()
        assert seq[0] == "9"

    @pytest.mark.parametrize(
        "item",
        [
            Menu.NCI_ANALYSIS,
            Menu.IRI_ANALYSIS,
            Menu.DORI_ANALYSIS,
            Menu.IGM_ANALYSIS,
        ],
    )
    def test_weak_interaction_menu_starts_with_20(
        self,
        item: Menu,
    ) -> None:
        """
        Test to ensure weak interaction analysis menu items start with '20'.

        :param self: Description
        :param item: Description
        :type item: Menu
        """
        seq = item.get_sequence()
        assert seq[0] == "20"

    @pytest.mark.parametrize(
        "item",
        [
            Menu.CUBE_DENSITY,
            Menu.CUBE_ESP,
            Menu.CUBE_ELF,
            Menu.CUBE_LAPLACIAN,
        ],
    )
    def test_cube_menu_starts_with_5(self, item: Menu) -> None:
        """
        Test to ensure cube generation menu items start with '5'.

        :param self: Description
        :param item: Description
        :type item: Menu
        """
        seq = item.get_sequence()
        assert seq[0] == "5"

    @pytest.mark.parametrize(
        "item",
        [
            Menu.PLOT_IR_SPECTRUM,
            Menu.PLOT_RAMAN_SPECTRUM,
            Menu.PLOT_UV_VIS_SPECTRUM,
        ],
    )
    def test_spectrum_menu_starts_with_11(self, item: Menu) -> None:
        """
        Test to ensure spectrum menu items start with '11'.

        :param self: Description
        :param item: Description
        :type item: Menu
        """
        seq = item.get_sequence()
        assert seq[0] == "11"


class TestSearchAndListing:
    """Tests for search() and list_all()."""

    def test_search_returns_matching_items(self) -> None:
        """
        Test to ensure search returns menu items matching the query.

        :param self: Description
        """
        results = Menu.search("hirshfeld")
        assert isinstance(results, list)
        assert Menu.HIRSHFELD_CHARGE in results

    def test_search_case_insensitive(self) -> None:
        """
        Test to ensure search is case-insensitive.

        :param self: Description
        """
        results = Menu.search("BoNd")
        assert any("BOND" in item.name for item in results)

    def test_list_all_returns_names(self) -> None:
        """
        Test to ensure list_all returns a list of all menu item names.

        :param self: Description
        """
        names = Menu.list_all()
        assert isinstance(names, list)
        assert "HIRSHFELD_CHARGE" in names
        assert "MAYER_BOND_ORDER" in names


class TestCustomSequence:
    """Tests for custom_sequence helper."""

    def test_custom_sequence_appends_q(self) -> None:
        """
        Test to ensure custom_sequence appends 'q' to the end of the sequence.

        :param self: Description
        """
        seq = custom_sequence(["7", "1"])
        assert seq == ["7", "1", "q"]

    def test_custom_sequence_keeps_existing_q(self) -> None:
        """
        Test to ensure custom_sequence does not append 'q' if already present.

        :param self: Description
        """
        seq = custom_sequence(["7", "1", "q"])
        assert seq == ["7", "1", "q"]
        assert seq.count("q") == 1

    def test_custom_sequence_requires_list(self) -> None:
        """
        Docstring for test_custom_sequence_requires_list.

        :param self: Description
        """
        with pytest.raises(TypeError, match="seq must be a list"):
            custom_sequence("7")  # type: ignore[arg-type]
