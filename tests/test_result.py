"""Tests for result.py - MultiwfnResult."""

import pytest
from pathlib import Path

from pyMultiwfn import MultiwfnResult


class TestMultiwfnResult:
    """Tests for MultiwfnResult."""
    
    @pytest.fixture
    def sample_result(self, sample_hirshfeld_output):
        """Create a sample result."""
        return MultiwfnResult(
            stdout=sample_hirshfeld_output,
            stderr="",
            returncode=0,
            execution_time=1.5,
            commands=["7", "1", "1", "q"],
            input_file=Path("/test/mol.wfx")
        )
    
    @pytest.fixture
    def failed_result(self):
        """Create a failed result."""
        return MultiwfnResult(
            stdout="",
            stderr="Error occurred",
            returncode=1,
            execution_time=0.5,
            commands=["7", "1", "q"],
            input_file=Path("/test/mol.wfx")
        )
    
    def test_success_property_true(self, sample_result):
        """Success is True when returncode is 0."""
        assert sample_result.success is True
    
    def test_success_property_false(self, failed_result):
        """Success is False when returncode is non-zero."""
        assert failed_result.success is False
    
    def test_parse_charges(self, sample_result):
        """Result can parse charges."""
        charges = sample_result.parse_charges("hirshfeld")
        
        assert len(charges) == 6
        assert charges[1] == pytest.approx(-0.0523)
    
    def test_parse_bond_orders(self, sample_mayer_output):
        """Result can parse bond orders."""
        result = MultiwfnResult(
            stdout=sample_mayer_output,
            stderr="",
            returncode=0,
            execution_time=1.0,
            commands=["9", "1", "q"],
            input_file=Path("/test/mol.wfx")
        )
        
        bonds = result.parse_bond_orders()
        assert len(bonds) == 5
    
    def test_parse_critical_points(self, sample_topology_output):
        """Result can parse critical points."""
        result = MultiwfnResult(
            stdout=sample_topology_output,
            stderr="",
            returncode=0,
            execution_time=5.0,
            commands=["2", "2", "-1", "q"],
            input_file=Path("/test/mol.wfx")
        )
        
        cps = result.parse_critical_points()
        assert len(cps) == 4
    
    def test_parse_spectrum(self, sample_spectrum_output):
        """Result can parse spectrum."""
        result = MultiwfnResult(
            stdout=sample_spectrum_output,
            stderr="",
            returncode=0,
            execution_time=2.0,
            commands=["11", "1", "q"],
            input_file=Path("/test/mol.wfx")
        )
        
        spectrum = result.parse_spectrum()
        assert 'frequencies' in spectrum
        assert len(spectrum['frequencies']) == 5
    
    def test_save_output(self, sample_result, temp_dir):
        """Result can save output to file."""
        output_file = temp_dir / "output.log"
        sample_result.save_output(output_file)
        
        assert output_file.exists()
        content = output_file.read_text()
        assert "Hirshfeld" in content
    
    def test_str_representation(self, sample_result):
        """Result has string representation."""
        s = str(sample_result)
        
        assert "SUCCESS" in s
        assert "1.50s" in s
        assert "4 commands" in s
    
    def test_str_failed(self, failed_result):
        """Failed result shows FAILED."""
        s = str(failed_result)
        assert "FAILED" in s
    
    def test_repr(self, sample_result):
        """Result has detailed repr."""
        r = repr(sample_result)
        
        assert "MultiwfnResult" in r
        assert "returncode=0" in r
        assert "execution_time=" in r


class TestMultiwfnResultGenericParse:
    """Tests for generic parse method."""
    
    def test_parse_with_parser_name(self, sample_hirshfeld_output):
        """Parse using parser name."""
        result = MultiwfnResult(
            stdout=sample_hirshfeld_output,
            stderr="",
            returncode=0,
            execution_time=1.0,
            commands=[],
            input_file=Path("/test/mol.wfx")
        )
        
        charges = result.parse('charges', method='hirshfeld')
        assert len(charges) == 6
    
    def test_parse_unknown_parser(self):
        """Unknown parser raises error."""
        from pyMultiwfn import MultiwfnError
        
        result = MultiwfnResult(
            stdout="",
            stderr="",
            returncode=0,
            execution_time=0,
            commands=[],
            input_file=Path("/test/mol.wfx")
        )
        
        with pytest.raises(MultiwfnError, match="Unknown parser"):
            result.parse('nonexistent_parser')