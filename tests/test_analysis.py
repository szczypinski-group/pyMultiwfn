"""Tests for analysis.py - Analysis helper functions."""

import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock

from pyMultiwfn import (
    list_functions,
    get_function,
    run_analysis,
    run_menu_sequence,
    run_custom_sequence,
    quick_analysis,
    charge_analysis,
    bond_order_analysis,
    topology_analysis,
    nci_analysis,
    menu,
)


class TestListFunctions:
    """Tests for list_functions."""
    
    def test_returns_list(self):
        """list_functions returns a list."""
        funcs = list_functions()
        assert isinstance(funcs, list)
    
    def test_contains_expected_functions(self):
        """List contains expected function names."""
        funcs = list_functions()
        
        assert 'hirshfeld_charge' in funcs
        assert 'mayer_bond_order' in funcs
        assert 'nci_analysis' in funcs
        assert 'topology_search_cps' in funcs
    
    def test_sorted(self):
        """List is sorted alphabetically."""
        funcs = list_functions()
        assert funcs == sorted(funcs)
    
    def test_no_private_functions(self):
        """List excludes private functions."""
        funcs = list_functions()
        assert not any(f.startswith('_') for f in funcs)


class TestGetFunction:
    """Tests for get_function."""
    
    def test_get_existing_function(self):
        """get_function returns callable for valid name."""
        func = get_function('hirshfeld_charge')
        assert callable(func)
        assert func is menu.hirshfeld_charge
    
    def test_get_nonexistent_function(self):
        """get_function raises ValueError for invalid name."""
        with pytest.raises(ValueError, match="Unknown function"):
            get_function('nonexistent_function')


class TestRunAnalysis:
    """Tests for run_analysis function."""
    
    def test_run_analysis_by_name(self, mock_wfn_file, mock_executable):
        """run_analysis executes named function."""
        with patch('subprocess.Popen') as mock_popen:
            mock_process = MagicMock()
            mock_process.communicate.return_value = (b"output", b"")
            mock_process.returncode = 0
            mock_popen.return_value = mock_process
            
            result = run_analysis(
                mock_wfn_file, 
                'hirshfeld_charge',
                exe_path=mock_executable
            )
        
        assert result is not None
        assert result.returncode == 0
    
    def test_run_analysis_with_func_kwargs(self, mock_wfn_file, mock_executable):
        """run_analysis passes func_kwargs to menu function."""
        with patch('subprocess.Popen') as mock_popen:
            mock_process = MagicMock()
            mock_process.communicate.return_value = (b"output", b"")
            mock_process.returncode = 0
            mock_popen.return_value = mock_process
            
            result = run_analysis(
                mock_wfn_file,
                'cube_orbital',
                func_kwargs={'orbital_index': 5},
                exe_path=mock_executable
            )
        
        assert result is not None
    
    def test_run_analysis_invalid_name(self, mock_wfn_file, mock_executable):
        """run_analysis raises error for invalid function name."""
        with pytest.raises(ValueError, match="Unknown function"):
            run_analysis(mock_wfn_file, 'nonexistent', exe_path=mock_executable)


class TestRunMenuSequence:
    """Tests for run_menu_sequence function."""
    
    def test_run_menu_sequence(self, mock_wfn_file, mock_executable):
        """run_menu_sequence executes menu function."""
        with patch('subprocess.Popen') as mock_popen:
            mock_process = MagicMock()
            mock_process.communicate.return_value = (b"output", b"")
            mock_process.returncode = 0
            mock_popen.return_value = mock_process
            
            result = run_menu_sequence(
                mock_wfn_file,
                menu.hirshfeld_charge,
                exe_path=mock_executable
            )
        
        assert result is not None
        assert result.returncode == 0


class TestRunCustomSequence:
    """Tests for run_custom_sequence function."""
    
    def test_run_custom_sequence(self, mock_wfn_file, mock_executable):
        """run_custom_sequence executes custom commands."""
        with patch('subprocess.Popen') as mock_popen:
            mock_process = MagicMock()
            mock_process.communicate.return_value = (b"output", b"")
            mock_process.returncode = 0
            mock_popen.return_value = mock_process
            
            result = run_custom_sequence(
                mock_wfn_file,
                ["7", "1", "1", "q"],
                exe_path=mock_executable
            )
        
        assert result is not None


class TestQuickAnalysis:
    """Tests for quick_analysis function."""
    
    def test_quick_analysis_charges(self, mock_wfn_file, mock_executable, sample_hirshfeld_output):
        """quick_analysis parses charge results."""
        with patch('subprocess.Popen') as mock_popen:
            mock_process = MagicMock()
            mock_process.communicate.return_value = (sample_hirshfeld_output.encode(), b"")
            mock_process.returncode = 0
            mock_popen.return_value = mock_process
            
            results = quick_analysis(
                mock_wfn_file,
                'charges',
                method='hirshfeld',
                exe_path=mock_executable
            )
        
        assert 'charges' in results
        assert 'stdout' in results
        assert 'success' in results
        assert results['success'] is True
    
    def test_quick_analysis_bonds(self, mock_wfn_file, mock_executable, sample_mayer_output):
        """quick_analysis parses bond order results."""
        with patch('subprocess.Popen') as mock_popen:
            mock_process = MagicMock()
            mock_process.communicate.return_value = (sample_mayer_output.encode(), b"")
            mock_process.returncode = 0
            mock_popen.return_value = mock_process
            
            results = quick_analysis(
                mock_wfn_file,
                'bonds',
                method='mayer',
                exe_path=mock_executable
            )
        
        assert 'bond_orders' in results
    
    def test_quick_analysis_topology(self, mock_wfn_file, mock_executable, sample_topology_output):
        """quick_analysis parses topology results."""
        with patch('subprocess.Popen') as mock_popen:
            mock_process = MagicMock()
            mock_process.communicate.return_value = (sample_topology_output.encode(), b"")
            mock_process.returncode = 0
            mock_popen.return_value = mock_process
            
            results = quick_analysis(
                mock_wfn_file,
                'topology',
                exe_path=mock_executable
            )
        
        assert 'critical_points' in results
    
    def test_quick_analysis_nci(self, mock_wfn_file, mock_executable):
        """quick_analysis returns NCI output."""
        with patch('subprocess.Popen') as mock_popen:
            mock_process = MagicMock()
            mock_process.communicate.return_value = (b"NCI output", b"")
            mock_process.returncode = 0
            mock_popen.return_value = mock_process
            
            results = quick_analysis(
                mock_wfn_file,
                'nci',
                exe_path=mock_executable
            )
        
        assert 'output' in results
    
    def test_quick_analysis_invalid_type(self, mock_wfn_file, mock_executable):
        """quick_analysis raises error for invalid type."""
        with pytest.raises(ValueError, match="Unknown analysis type"):
            quick_analysis(mock_wfn_file, 'invalid_type', exe_path=mock_executable)


class TestConvenienceFunctions:
    """Tests for convenience analysis functions."""
    
    def test_charge_analysis(self, mock_wfn_file, mock_executable, sample_hirshfeld_output):
        """charge_analysis returns charges dict."""
        with patch('subprocess.Popen') as mock_popen:
            mock_process = MagicMock()
            mock_process.communicate.return_value = (sample_hirshfeld_output.encode(), b"")
            mock_process.returncode = 0
            mock_popen.return_value = mock_process
            
            charges = charge_analysis(
                mock_wfn_file,
                method='hirshfeld',
                exe_path=mock_executable
            )
        
        assert isinstance(charges, dict)
        assert len(charges) == 6
    
    def test_bond_order_analysis(self, mock_wfn_file, mock_executable, sample_mayer_output):
        """bond_order_analysis returns bond orders dict."""
        with patch('subprocess.Popen') as mock_popen:
            mock_process = MagicMock()
            mock_process.communicate.return_value = (sample_mayer_output.encode(), b"")
            mock_process.returncode = 0
            mock_popen.return_value = mock_process
            
            bonds = bond_order_analysis(
                mock_wfn_file,
                method='mayer',
                exe_path=mock_executable
            )
        
        assert isinstance(bonds, dict)
    
    def test_topology_analysis(self, mock_wfn_file, mock_executable, sample_topology_output):
        """topology_analysis returns critical points list."""
        with patch('subprocess.Popen') as mock_popen:
            mock_process = MagicMock()
            mock_process.communicate.return_value = (sample_topology_output.encode(), b"")
            mock_process.returncode = 0
            mock_popen.return_value = mock_process
            
            cps = topology_analysis(mock_wfn_file, exe_path=mock_executable)
        
        assert isinstance(cps, list)
    
    def test_nci_analysis(self, mock_wfn_file, mock_executable):
        """nci_analysis returns output string."""
        with patch('subprocess.Popen') as mock_popen:
            mock_process = MagicMock()
            mock_process.communicate.return_value = (b"NCI analysis output", b"")
            mock_process.returncode = 0
            mock_popen.return_value = mock_process
            
            output = nci_analysis(mock_wfn_file, exe_path=mock_executable)
        
        assert output is not None
        assert isinstance(output, str)