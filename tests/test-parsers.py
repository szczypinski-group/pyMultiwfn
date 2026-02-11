"""Tests for parsers.py - Output parsers and ParserRegistry."""

import pytest

from pyMultiwfn import (
    ParserRegistry,
    ChargeParser,
    BondOrderParser,
    CriticalPointParser,
    SpectrumParser,
    MultiwfnError,
)


class TestChargeParser:
    """Tests for ChargeParser."""
    
    def test_parse_hirshfeld(self, sample_hirshfeld_output):
        """Parse Hirshfeld charges."""
        charges = ChargeParser.parse(sample_hirshfeld_output, method="Hirshfeld")
        
        assert len(charges) == 6
        assert charges[1] == pytest.approx(-0.0523)
        assert charges[5] == pytest.approx(-0.3245)
        assert charges[6] == pytest.approx(0.1823)
    
    def test_parse_mulliken(self, sample_mulliken_output):
        """Parse Mulliken charges."""
        charges = ChargeParser.parse(sample_mulliken_output, method="Mulliken")
        
        assert len(charges) == 6
        assert charges[1] == pytest.approx(-0.1234)
        assert charges[5] == pytest.approx(-0.4567)
    
    def test_parse_empty_output(self):
        """Empty output returns empty dict."""
        charges = ChargeParser.parse("", method="Hirshfeld")
        assert charges == {}
    
    def test_parse_no_charges_section(self):
        """Output without charges section returns empty dict."""
        output = "Some random output\nNo charges here"
        charges = ChargeParser.parse(output, method="Hirshfeld")
        assert charges == {}
    
    def test_parse_case_insensitive(self, sample_hirshfeld_output):
        """Method name is case-insensitive."""
        charges1 = ChargeParser.parse(sample_hirshfeld_output, method="hirshfeld")
        charges2 = ChargeParser.parse(sample_hirshfeld_output, method="HIRSHFELD")
        
        assert charges1 == charges2


class TestBondOrderParser:
    """Tests for BondOrderParser."""
    
    def test_parse_mayer(self, sample_mayer_output):
        """Parse Mayer bond orders."""
        bonds = BondOrderParser.parse(sample_mayer_output)
        
        assert len(bonds) == 5
        assert bonds[(1, 2)] == pytest.approx(1.4523)
        assert bonds[(2, 5)] == pytest.approx(1.8934)
        assert bonds[(5, 6)] == pytest.approx(0.8234)
    
    def test_parse_wiberg(self, sample_wiberg_output):
        """Parse Wiberg bond orders."""
        bonds = BondOrderParser.parse(sample_wiberg_output)
        
        assert len(bonds) == 3
        assert bonds[(1, 2)] == pytest.approx(1.5234)
    
    def test_parse_empty_output(self):
        """Empty output returns empty dict."""
        bonds = BondOrderParser.parse("")
        assert bonds == {}
    
    def test_bond_tuple_ordering(self, sample_mayer_output):
        """Bond tuples should have smaller index first."""
        bonds = BondOrderParser.parse(sample_mayer_output)
        
        for (a, b) in bonds.keys():
            assert a < b


class TestCriticalPointParser:
    """Tests for CriticalPointParser."""
    
    def test_parse_all_cp_types(self, sample_topology_output):
        """Parse all critical point types."""
        cps = CriticalPointParser.parse(sample_topology_output)
        
        assert len(cps) == 4
    
    def test_nuclear_critical_point(self, sample_topology_output):
        """Parse nuclear (3,-3) critical point."""
        cps = CriticalPointParser.parse(sample_topology_output)
        nuclear = [cp for cp in cps if cp['type'] == '(3,-3)'][0]
        
        assert nuclear['cp_type'] == 'nuclear'
        assert nuclear['position'] == (0.0, 0.0, 0.0)
    
    def test_bond_critical_point(self, sample_topology_output):
        """Parse bond (3,-1) critical point."""
        cps = CriticalPointParser.parse(sample_topology_output)
        bond = [cp for cp in cps if cp['type'] == '(3,-1)'][0]
        
        assert bond['cp_type'] == 'bond'
        assert bond['position'][0] == pytest.approx(1.234567)
    
    def test_ring_critical_point(self, sample_topology_output):
        """Parse ring (3,+1) critical point."""
        cps = CriticalPointParser.parse(sample_topology_output)
        ring = [cp for cp in cps if cp['type'] == '(3,+1)'][0]
        
        assert ring['cp_type'] == 'ring'
    
    def test_cage_critical_point(self, sample_topology_output):
        """Parse cage (3,+3) critical point."""
        cps = CriticalPointParser.parse(sample_topology_output)
        cage = [cp for cp in cps if cp['type'] == '(3,+3)'][0]
        
        assert cage['cp_type'] == 'cage'
    
    def test_parse_empty_output(self):
        """Empty output returns empty list."""
        cps = CriticalPointParser.parse("")
        assert cps == []


class TestSpectrumParser:
    """Tests for SpectrumParser."""
    
    def test_parse_spectrum(self, sample_spectrum_output):
        """Parse IR spectrum data."""
        spectrum = SpectrumParser.parse(sample_spectrum_output)
        
        assert 'frequencies' in spectrum
        assert 'intensities' in spectrum
        assert len(spectrum['frequencies']) == 5
        assert len(spectrum['intensities']) == 5
    
    def test_spectrum_values(self, sample_spectrum_output):
        """Verify spectrum values."""
        spectrum = SpectrumParser.parse(sample_spectrum_output)
        
        assert spectrum['frequencies'][0] == pytest.approx(500.0)
        assert spectrum['intensities'][0] == pytest.approx(12.34)
        assert spectrum['frequencies'][-1] == pytest.approx(3000.0)
        assert spectrum['intensities'][-1] == pytest.approx(156.78)
    
    def test_parse_empty_output(self):
        """Empty output returns empty spectrum."""
        spectrum = SpectrumParser.parse("")
        assert spectrum == {'frequencies': [], 'intensities': []}


class TestParserRegistry:
    """Tests for ParserRegistry."""
    
    def test_list_parsers(self):
        """List available parsers."""
        parsers = ParserRegistry.list_parsers()
        
        assert 'charges' in parsers
        assert 'bond_orders' in parsers
        assert 'critical_points' in parsers
        assert 'spectrum' in parsers
    
    def test_get_parser(self):
        """Get parser by name."""
        parser = ParserRegistry.get_parser('charges')
        assert parser == ChargeParser
    
    def test_get_unknown_parser(self):
        """Unknown parser raises error."""
        with pytest.raises(MultiwfnError, match="Unknown parser"):
            ParserRegistry.get_parser('nonexistent')
    
    def test_parse_through_registry(self, sample_hirshfeld_output):
        """Parse using registry."""
        charges = ParserRegistry.parse('charges', sample_hirshfeld_output, method="Hirshfeld")
        
        assert len(charges) == 6
        assert charges[1] == pytest.approx(-0.0523)
    
    def test_register_custom_parser(self):
        """Register and use custom parser."""
        class CustomParser:
            @staticmethod
            def parse(output, **kwargs):
                return {"custom": True}
        
        ParserRegistry.register('custom_test', CustomParser)
        
        assert 'custom_test' in ParserRegistry.list_parsers()
        result = ParserRegistry.parse('custom_test', "any output")
        assert result == {"custom": True}
        
        # Cleanup
        del ParserRegistry._parsers['custom_test']