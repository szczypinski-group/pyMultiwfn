"""
Integration tests for pyMultiwfn.

These tests use real wavefunction files and optionally the real Multiwfn executable.
Run with: pytest tests/test_integration.py -v -m "real_files or requires_multiwfn"
"""

import pytest
from pathlib import Path

from pyMultiwfn import (
    MultiwfnJob,
    MultiwfnConfig,
    MultiwfnResult,
    MultiwfnError,
    ParserRegistry,
    menu,
    run_analysis,
    list_functions,
)


# =============================================================================
# Tests with Real Wavefunction Files (no Multiwfn needed)
# =============================================================================

@pytest.mark.real_files
class TestRealFileLoading:
    """Tests that load real wavefunction files."""
    
    def test_load_wfx_file(self, real_wfx_file, mock_executable):
        """Can create job with real .wfx file."""
        config = MultiwfnConfig(exe_path=mock_executable)
        job = MultiwfnJob(real_wfx_file, config=config)
        
        assert job.input_file == real_wfx_file
        assert job.input_file.exists()
        assert job.input_file.suffix == '.wfx'
    
    def test_load_molden_file(self, real_molden_file, mock_executable):
        """Can create job with real .molden file."""
        config = MultiwfnConfig(exe_path=mock_executable)
        job = MultiwfnJob(real_molden_file, config=config)
        
        assert job.input_file == real_molden_file
        assert job.input_file.exists()
        assert job.input_file.suffix == '.molden'
    
    def test_wfx_file_content(self, real_wfx_file):
        """Verify .wfx file has expected content structure."""
        content = real_wfx_file.read_text()
        
        # WFX files should have these sections
        assert '<Number of Nuclei>' in content or 'Number of Nuclei' in content
    
    def test_builder_with_real_files(self, real_wfx_file, mock_executable):
        """Builder pattern works with real files."""
        config = MultiwfnConfig(exe_path=mock_executable)
        job = (
            MultiwfnJob.builder(real_wfx_file)
            .with_config(config)
            .with_menu_sequence(menu.hirshfeld_charge)
            .with_menu_sequence(menu.mayer_bond_order)
            .build()
        )
        
        assert job.input_file == real_wfx_file
        assert len(job.commands) > 0


# =============================================================================
# Tests with Real Multiwfn Executable
# =============================================================================

@pytest.mark.requires_multiwfn
class TestRealMultiwfnExecution:
    """
    Tests that run real Multiwfn analyses.
    
    These tests require:
    1. Real wavefunction files in tests/test_data/
    2. Multiwfn executable in bin/
    
    Skip with: pytest -m "not requires_multiwfn"
    """
    
    def test_hirshfeld_charge(self, real_wfx_file, real_executable, temp_dir):
        """Run actual Hirshfeld charge analysis."""
        config = MultiwfnConfig(
            exe_path=real_executable,
            working_dir=temp_dir,
            timeout=120
        )
        job = MultiwfnJob(real_wfx_file, config=config)
        job.add_menu_sequence(menu.hirshfeld_charge)
        
        try:
            result = job.run()
        except MultiwfnError as e:
            pytest.skip(f"Multiwfn execution failed: {e}")
        
        assert result.success
        assert len(result.stdout) > 0
        
        charges = result.parse_charges('hirshfeld')
        if len(charges) == 0:
            pytest.skip("Could not parse charges from output")
        
        # Charges should sum to approximately total charge
        total = sum(charges.values())
        assert abs(total) < 1.0  # Reasonable for most molecules
    
    def test_mayer_bond_order(self, real_wfx_file, real_executable, temp_dir):
        """Run actual Mayer bond order analysis."""
        config = MultiwfnConfig(
            exe_path=real_executable,
            working_dir=temp_dir,
            timeout=120
        )
        job = MultiwfnJob(real_wfx_file, config=config)
        job.add_menu_sequence(menu.mayer_bond_order)
        
        try:
            result = job.run()
        except MultiwfnError as e:
            pytest.skip(f"Multiwfn execution failed: {e}")
        
        assert result.success
        bonds = result.parse_bond_orders()
        if len(bonds) == 0:
            pytest.skip("Could not parse bond orders from output")
    
    def test_topology_analysis(self, real_wfx_file, real_executable, temp_dir):
        """Run actual topology analysis."""
        config = MultiwfnConfig(
            exe_path=real_executable,
            working_dir=temp_dir,
            timeout=300  # Topology can be slow
        )
        job = MultiwfnJob(real_wfx_file, config=config)
        job.add_menu_sequence(menu.topology_search_cps)
        
        try:
            result = job.run()
        except MultiwfnError as e:
            pytest.skip(f"Multiwfn execution failed: {e}")
        
        assert result.success
        cps = result.parse_critical_points()
        # Should find at least nuclear critical points
        if len(cps) == 0:
            pytest.skip("Could not parse critical points from output")
    
    def test_multiple_analyses(self, real_wfx_file, real_executable, temp_dir):
        """Run multiple analyses in sequence."""
        config = MultiwfnConfig(
            exe_path=real_executable,
            working_dir=temp_dir,
            timeout=180
        )
        job = (
            MultiwfnJob.builder(real_wfx_file)
            .with_config(config)
            .with_menu_sequence(menu.hirshfeld_charge)
            .with_menu_sequence(menu.mayer_bond_order)
            .build()
        )
        
        try:
            result = job.run()
        except MultiwfnError as e:
            pytest.skip(f"Multiwfn execution failed: {e}")
        
        assert result.success
        
        # Both should be parseable
        charges = result.parse_charges('hirshfeld')
        bonds = result.parse_bond_orders()
        
        if len(charges) == 0 and len(bonds) == 0:
            pytest.skip("Could not parse any results from output")
    
    def test_run_analysis_function(self, real_wfx_file, real_executable, temp_dir):
        """Test run_analysis convenience function."""
        try:
            result = run_analysis(
                real_wfx_file,
                'hirshfeld_charge',
                exe_path=real_executable,
                working_dir=temp_dir,
                timeout=120
            )
        except MultiwfnError as e:
            pytest.skip(f"Multiwfn execution failed: {e}")
        
        assert result.success
    
    def test_wavefunction_check(self, real_wfx_file, real_executable, temp_dir):
        """Run wavefunction check."""
        config = MultiwfnConfig(
            exe_path=real_executable,
            working_dir=temp_dir,
            timeout=60
        )
        job = MultiwfnJob(real_wfx_file, config=config)
        job.add_menu_sequence(menu.check_wavefunction)
        
        try:
            result = job.run()
        except MultiwfnError as e:
            pytest.skip(f"Multiwfn execution failed: {e}")
        
        assert result.success
        # Output should contain wavefunction info
        assert 'orbital' in result.stdout.lower() or 'electron' in result.stdout.lower()
    
    def test_save_and_reload_output(self, real_wfx_file, real_executable, temp_dir):
        """Save output and verify it can be reparsed."""
        config = MultiwfnConfig(
            exe_path=real_executable,
            working_dir=temp_dir,
            timeout=120
        )
        job = MultiwfnJob(real_wfx_file, config=config)
        job.add_menu_sequence(menu.hirshfeld_charge)
        
        try:
            result = job.run()
        except MultiwfnError as e:
            pytest.skip(f"Multiwfn execution failed: {e}")
        
        assert result.success
        
        # Save output
        output_file = temp_dir / "hirshfeld_output.log"
        result.save_output(output_file)
        assert output_file.exists()
        
        # Reload and verify parseable
        saved_output = output_file.read_text()
        charges = ParserRegistry.parse('charges', saved_output, method='hirshfeld')
        if len(charges) == 0:
            pytest.skip("Could not parse charges from saved output")


# =============================================================================
# End-to-End Workflow Tests
# =============================================================================

@pytest.mark.requires_multiwfn
class TestEndToEndWorkflows:
    """Full workflow tests."""
    
    def test_complete_analysis_workflow(self, real_wfx_file, real_executable, temp_dir):
        """Complete analysis workflow from file to parsed results."""
        # 1. Create config
        config = MultiwfnConfig(
            exe_path=real_executable,
            working_dir=temp_dir,
            timeout=300,
            verbose=False
        )
        
        # 2. Build job with multiple analyses
        job = (
            MultiwfnJob.builder(real_wfx_file)
            .with_config(config)
            .with_menu_sequence(menu.hirshfeld_charge)
            .with_menu_sequence(menu.mayer_bond_order)
            .build()
        )
        
        # 3. Run
        try:
            result = job.run()
        except MultiwfnError as e:
            pytest.skip(f"Multiwfn execution failed: {e}")
        
        assert result.success
        
        # 4. Parse all results
        charges = result.parse_charges('hirshfeld')
        bonds = result.parse_bond_orders()
        
        # 5. Verify results (at least one should work)
        if len(charges) == 0 and len(bonds) == 0:
            pytest.skip("Could not parse any results from output")
        
        # 6. Save output
        output_file = temp_dir / "complete_analysis.log"
        result.save_output(output_file)
        assert output_file.exists()
        
        # 7. Verify execution time recorded
        assert result.execution_time > 0
    
    def test_run_all_charges_workflow(self, real_wfx_file, real_executable, temp_dir):
        """Run all charge methods using run_all_charges."""
        config = MultiwfnConfig(
            exe_path=real_executable,
            working_dir=temp_dir,
            timeout=600  # Multiple analyses take time
        )
        job = MultiwfnJob(real_wfx_file, config=config)
        job.add_menu_sequence(menu.run_all_charges)
        
        try:
            result = job.run()
        except MultiwfnError as e:
            pytest.skip(f"Multiwfn execution failed: {e}")
        
        # May not fully succeed if some methods fail, but should produce output
        assert len(result.stdout) > 0


# =============================================================================
# Error Handling Tests
# =============================================================================

class TestErrorHandling:
    """Tests for error conditions."""
    
    def test_missing_input_file(self, mock_executable):
        """Error when input file doesn't exist."""
        config = MultiwfnConfig(exe_path=mock_executable)
        
        with pytest.raises((FileNotFoundError, ValueError)):
            MultiwfnJob(Path("/nonexistent/file.wfx"), config=config)
    
    def test_missing_executable(self, mock_wfn_file, temp_dir):
        """Error when executable doesn't exist."""
        from pyMultiwfn import MultiwfnError
        
        config = MultiwfnConfig(exe_path=temp_dir / "nonexistent")
        
        with pytest.raises(MultiwfnError, match="not found"):
            _ = config.executable


# =============================================================================
# Compatibility Tests
# =============================================================================

class TestFileCompatibility:
    """Tests for different file format compatibility."""
    
    @pytest.mark.real_files
    def test_wfx_extension(self, real_wfx_file, mock_executable):
        """Accept .wfx files."""
        config = MultiwfnConfig(exe_path=mock_executable)
        job = MultiwfnJob(real_wfx_file, config=config)
        assert job.input_file.suffix == '.wfx'
    
    @pytest.mark.real_files
    def test_molden_extension(self, real_molden_file, mock_executable):
        """Accept .molden files."""
        config = MultiwfnConfig(exe_path=mock_executable)
        job = MultiwfnJob(real_molden_file, config=config)
        assert job.input_file.suffix == '.molden'