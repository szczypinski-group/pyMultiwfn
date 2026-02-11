"""Tests for menu.py - Menu sequences and run_all functions."""

import pytest

from pyMultiwfn import menu


class TestMenuFunctionStructure:
    """Test that all menu functions return valid sequences."""
    
    def test_all_functions_return_lists(self):
        """All no-arg menu functions should return lists."""
        # Get all callable functions that don't start with _
        funcs = [
            getattr(menu, name) 
            for name in dir(menu) 
            if callable(getattr(menu, name)) 
            and not name.startswith('_')
            and name not in ('list_all_functions', 'get_function_by_category')
        ]
        
        for func in funcs:
            try:
                # Try calling with no args
                result = func()
                if isinstance(result, list):
                    assert all(isinstance(item, str) for item in result), f"{func.__name__} returned non-string items"
            except TypeError:
                # Function requires arguments, skip
                pass
    
    def test_sequences_end_with_q(self):
        """Most sequences should end with 'q' to exit."""
        simple_funcs = [
            menu.hirshfeld_charge,
            menu.mayer_bond_order,
            menu.nci_analysis,
            menu.topology_search_cps,
            menu.cube_density,
        ]
        
        for func in simple_funcs:
            seq = func()
            assert seq[-1] == 'q', f"{func.__name__} doesn't end with 'q'"


class TestChargeMenuFunctions:
    """Tests for charge/population menu functions."""
    
    @pytest.mark.parametrize("func,expected_menu", [
        (menu.hirshfeld_charge, "7"),
        (menu.mulliken_population, "7"),
        (menu.lowdin_population, "7"),
        (menu.adch_charge, "7"),
        (menu.resp_charge, "7"),
        (menu.cm5_charge, "7"),
        (menu.chelpg_charge, "7"),
        (menu.mk_charge, "7"),
        (menu.becke_charge, "7"),
        (menu.vdd_population, "7"),
        (menu.aim_charge, "7"),
        (menu.eem_charge, "7"),
        (menu.gasteiger_charge, "7"),
        (menu.mbis_charge, "7"),
    ])
    def test_charge_functions_use_menu_7(self, func, expected_menu):
        """Charge functions should start with menu 7."""
        seq = func()
        assert seq[0] == expected_menu


class TestBondOrderMenuFunctions:
    """Tests for bond order menu functions."""
    
    @pytest.mark.parametrize("func,expected_menu", [
        (menu.mayer_bond_order, "9"),
        (menu.wiberg_bond_order, "9"),
        (menu.fuzzy_bond_order, "9"),
        (menu.laplacian_bond_order, "9"),
        (menu.mulliken_bond_order, "9"),
        (menu.multicenter_bond_order, "9"),
    ])
    def test_bond_functions_use_menu_9(self, func, expected_menu):
        """Bond order functions should start with menu 9."""
        seq = func()
        assert seq[0] == expected_menu


class TestTopologyMenuFunctions:
    """Tests for topology menu functions."""
    
    @pytest.mark.parametrize("func", [
        menu.topology_search_cps,
        menu.topology_generate_paths,
        menu.topology_interbasin_surfaces,
        menu.topology_analysis_complete,
    ])
    def test_topology_functions_use_menu_2(self, func):
        """Topology functions should start with menu 2."""
        seq = func()
        assert seq[0] == "2"


class TestWeakInteractionMenuFunctions:
    """Tests for weak interaction menu functions."""
    
    @pytest.mark.parametrize("func", [
        menu.nci_analysis,
        menu.igm_analysis,
        menu.iri_analysis,
        menu.dori_analysis,
    ])
    def test_weak_interaction_functions_use_menu_20(self, func):
        """Weak interaction functions should start with menu 20."""
        seq = func()
        assert seq[0] == "20"


class TestCubeMenuFunctions:
    """Tests for cube generation menu functions."""
    
    @pytest.mark.parametrize("func", [
        menu.cube_density,
        menu.cube_esp,
        menu.cube_elf,
        menu.cube_lol,
        menu.cube_laplacian,
    ])
    def test_cube_functions_use_menu_5(self, func):
        """Cube functions should start with menu 5."""
        seq = func()
        assert seq[0] == "5"


class TestSpectrumMenuFunctions:
    """Tests for spectrum menu functions."""
    
    @pytest.mark.parametrize("func", [
        menu.plot_ir_spectrum,
        menu.plot_raman_spectrum,
        menu.plot_uv_vis_spectrum,
    ])
    def test_spectrum_functions_use_menu_11(self, func):
        """Spectrum functions should start with menu 11."""
        seq = func()
        assert seq[0] == "11"


class TestFunctionsWithArguments:
    """Tests for menu functions that require arguments."""
    
    def test_cube_orbital(self):
        """cube_orbital requires orbital index."""
        seq = menu.cube_orbital(5)
        assert "5" in seq  # menu 5
        assert seq[-1] == "q"
    
    def test_properties_at_point(self):
        """properties_at_point requires coordinates."""
        seq = menu.properties_at_point(0.0, 1.0, 2.0)
        assert seq[0] == "1"
        assert "0.0,1.0,2.0" in seq
    
    def test_view_orbital(self):
        """view_orbital with orbital index."""
        seq = menu.view_orbital(10)
        assert seq[0] == "0"
        assert "10" in seq
    
    def test_hole_electron_analysis(self):
        """hole_electron_analysis with state number."""
        seq = menu.hole_electron_analysis("2")
        assert seq[0] == "18"


class TestCustomMenu:
    """Tests for custom_menu function."""
    
    def test_custom_menu_with_list(self):
        """custom_menu accepts list."""
        seq = menu.custom_menu(["7", "1", "1"])
        assert seq == ["7", "1", "1", "q"]
    
    def test_custom_menu_with_string_raises(self):
        """custom_menu rejects string input."""
        with pytest.raises(TypeError, match="must be a list"):
            menu.custom_menu("7 1 1")
    
    def test_custom_menu_already_has_q(self):
        """custom_menu doesn't add extra q if present."""
        seq = menu.custom_menu(["7", "1", "q"])
        assert seq == ["7", "1", "q"]
        assert seq.count("q") == 1


class TestRunAllFunctions:
    """Tests for run_all* functions."""
    
    def test_run_all_returns_list(self):
        """run_all returns a list."""
        seq = menu.run_all()
        assert isinstance(seq, list)
        assert len(seq) > 0
    
    def test_run_all_ends_with_q(self):
        """run_all ends with single 'q'."""
        seq = menu.run_all()
        assert seq[-1] == 'q'
        assert seq.count('q') == 1
    
    def test_run_all_charges(self):
        """run_all_charges returns combined charge sequences."""
        seq = menu.run_all_charges()
        assert isinstance(seq, list)
        assert seq[-1] == 'q'
        # Should contain menu 7 entries
        assert '7' in seq
    
    def test_run_all_bond_orders(self):
        """run_all_bond_orders returns combined bond order sequences."""
        seq = menu.run_all_bond_orders()
        assert isinstance(seq, list)
        assert seq[-1] == 'q'
        # Should contain menu 9 entries
        assert '9' in seq
    
    def test_run_all_weak_interactions(self):
        """run_all_weak_interactions returns combined sequences."""
        seq = menu.run_all_weak_interactions()
        assert isinstance(seq, list)
        assert seq[-1] == 'q'
        # Should contain menu 20 entries
        assert '20' in seq
    
    def test_run_all_topology(self):
        """run_all_topology returns combined topology sequences."""
        seq = menu.run_all_topology()
        assert isinstance(seq, list)
        assert seq[-1] == 'q'
        assert '2' in seq
    
    def test_run_all_cubes(self):
        """run_all_cubes returns combined cube sequences."""
        seq = menu.run_all_cubes()
        assert isinstance(seq, list)
        assert seq[-1] == 'q'
        assert '5' in seq
    
    def test_run_all_spectra(self):
        """run_all_spectra returns combined spectrum sequences."""
        seq = menu.run_all_spectra()
        assert isinstance(seq, list)
        assert seq[-1] == 'q'
    
    def test_run_all_surfaces(self):
        """run_all_surfaces returns combined surface sequences."""
        seq = menu.run_all_surfaces()
        assert isinstance(seq, list)
        assert seq[-1] == 'q'
    
    def test_run_all_aromaticity(self):
        """run_all_aromaticity returns combined aromaticity sequences."""
        seq = menu.run_all_aromaticity()
        assert isinstance(seq, list)
        assert seq[-1] == 'q'
    
    def test_run_all_cdft(self):
        """run_all_cdft returns combined CDFT sequences."""
        seq = menu.run_all_cdft()
        assert isinstance(seq, list)
        assert seq[-1] == 'q'


class TestSuiteFunctions:
    """Tests for analysis suite functions."""
    
    def test_charge_analysis_suite(self):
        """charge_analysis_suite returns dict of sequences."""
        suite = menu.charge_analysis_suite()
        assert isinstance(suite, dict)
        assert 'hirshfeld' in suite
        assert 'adch' in suite
    
    def test_bond_analysis_suite(self):
        """bond_analysis_suite returns dict of sequences."""
        suite = menu.bond_analysis_suite()
        assert isinstance(suite, dict)
        assert 'mayer' in suite
        assert 'wiberg' in suite
    
    def test_weak_interaction_suite(self):
        """weak_interaction_suite returns combined list of sequences."""
        suite = menu.weak_interaction_suite()
        assert isinstance(suite, list)
        assert suite[-1] == 'q'
        assert '20' in suite  # weak interaction menu


class TestDocumentationHelpers:
    """Tests for documentation helper functions."""
    
    def test_list_all_functions(self):
        """list_all_functions returns dict."""
        funcs = menu.list_all_functions()
        assert isinstance(funcs, dict)
        assert len(funcs) > 100  # Should have many functions
        assert 'hirshfeld_charge' in funcs
    
    def test_get_function_by_category(self):
        """get_function_by_category returns categorized dict."""
        categories = menu.get_function_by_category()
        assert isinstance(categories, dict)
        assert 'Population' in categories or 'Visualization' in categories