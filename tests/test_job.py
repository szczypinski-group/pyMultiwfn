"""Tests for job.py - MultiwfnJob and MultiwfnJobBuilder."""

import platform
import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock

from pyMultiwfn import MultiwfnJob, MultiwfnJobBuilder, MultiwfnConfig, MultiwfnError, menu


class TestMultiwfnJobBuilder:
    """Tests for MultiwfnJobBuilder."""
    
    def test_builder_creation(self, mock_wfn_file):
        """Builder can be created with input file."""
        builder = MultiwfnJobBuilder(mock_wfn_file)
        assert builder._input_file == mock_wfn_file
    
    def test_builder_from_job(self, mock_wfn_file):
        """Builder can be created via MultiwfnJob.builder()."""
        builder = MultiwfnJob.builder(mock_wfn_file)
        assert isinstance(builder, MultiwfnJobBuilder)
    
    def test_with_config(self, mock_wfn_file, mock_executable, temp_dir):
        """Builder accepts config."""
        config = MultiwfnConfig(exe_path=mock_executable, working_dir=temp_dir)
        builder = MultiwfnJobBuilder(mock_wfn_file).with_config(config)
        
        assert builder._exe_path == mock_executable
        assert builder._working_dir == temp_dir
    
    def test_with_timeout(self, mock_wfn_file):
        """Builder accepts timeout."""
        builder = MultiwfnJobBuilder(mock_wfn_file).with_timeout(300)
        assert builder._timeout == 300
    
    def test_with_invalid_timeout(self, mock_wfn_file):
        """Builder rejects invalid timeout."""
        with pytest.raises(ValueError, match="positive"):
            MultiwfnJobBuilder(mock_wfn_file).with_timeout(-1)
        
        with pytest.raises(ValueError, match="positive"):
            MultiwfnJobBuilder(mock_wfn_file).with_timeout(0)
    
    def test_with_menu_sequence(self, mock_wfn_file):
        """Builder accepts menu sequence."""
        builder = MultiwfnJobBuilder(mock_wfn_file)
        builder.with_menu_sequence(menu.hirshfeld_charge)
        
        assert len(builder._commands) >= 1
    
    def test_chaining(self, mock_wfn_file):
        """Builder methods can be chained."""
        builder = (
            MultiwfnJobBuilder(mock_wfn_file)
            .with_timeout(300)
            .with_menu_sequence(menu.hirshfeld_charge)
            .with_menu_sequence(menu.mayer_bond_order)
        )
        
        assert builder._timeout == 300
        assert len(builder._commands) >= 2  # Has commands from both sequences
    
    def test_build(self, mock_wfn_file, mock_executable):
        """Builder creates job."""
        config = MultiwfnConfig(exe_path=mock_executable)
        job = (
            MultiwfnJobBuilder(mock_wfn_file)
            .with_config(config)
            .with_menu_sequence(menu.hirshfeld_charge)
            .build()
        )
        
        assert isinstance(job, MultiwfnJob)
        assert job.input_file == mock_wfn_file
    
    def test_build_requires_input_file(self):
        """Builder requires valid input file."""
        with pytest.raises((ValueError, FileNotFoundError)):
            MultiwfnJobBuilder(Path("/nonexistent/file.wfn")).build()


class TestMultiwfnJob:
    """Tests for MultiwfnJob."""
    
    def test_job_creation(self, mock_wfn_file, mock_executable):
        """Job can be created."""
        config = MultiwfnConfig(exe_path=mock_executable)
        job = MultiwfnJob(mock_wfn_file, config=config)
        
        assert job.input_file == mock_wfn_file
    
    def test_add_menu_sequence(self, mock_wfn_file, mock_executable):
        """Job accepts menu sequences."""
        config = MultiwfnConfig(exe_path=mock_executable)
        job = MultiwfnJob(mock_wfn_file, config=config)
        
        job.add_menu_sequence(menu.hirshfeld_charge)
        
        assert len(job.commands) > 0
        # Commands should contain menu entries (7 for charges)
    
    def test_add_multiple_sequences(self, mock_wfn_file, mock_executable):
        """Job combines multiple sequences."""
        config = MultiwfnConfig(exe_path=mock_executable)
        job = MultiwfnJob(mock_wfn_file, config=config)
        
        job.add_menu_sequence(menu.hirshfeld_charge)
        job.add_menu_sequence(menu.mayer_bond_order)
        
        # Should have commands from both
        assert len(job.commands) > 0
        assert '7' in job.commands  # charges menu
        assert '9' in job.commands  # bond order menu
    
    def test_add_custom_commands(self, mock_wfn_file, mock_executable):
        """Job accepts custom command list."""
        config = MultiwfnConfig(exe_path=mock_executable)
        job = MultiwfnJob(mock_wfn_file, config=config)
        
        job.add_custom_commands(["7", "1", "1", "q"])
        
        assert "7" in job.commands
    
    def test_clear_commands(self, mock_wfn_file, mock_executable):
        """Job can clear commands."""
        config = MultiwfnConfig(exe_path=mock_executable)
        job = MultiwfnJob(mock_wfn_file, config=config)
        
        job.add_menu_sequence(menu.hirshfeld_charge)
        assert len(job.commands) > 0
        
        job.clear_commands()
        assert len(job.commands) == 0
    
    def test_method_chaining(self, mock_wfn_file, mock_executable):
        """Job methods can be chained."""
        config = MultiwfnConfig(exe_path=mock_executable)
        job = MultiwfnJob(mock_wfn_file, config=config)
        
        result = job.add_menu_sequence(menu.hirshfeld_charge).add_menu_sequence(menu.mayer_bond_order)
        
        assert result is job


class TestMultiwfnJobExecution:
    """Tests for job execution (mocked)."""
    
    def test_run_returns_result(self, mock_wfn_file, mock_executable):
        """Run returns MultiwfnResult."""
        config = MultiwfnConfig(exe_path=mock_executable)
        job = MultiwfnJob(mock_wfn_file, config=config)
        job.add_menu_sequence(menu.hirshfeld_charge)
        
        with patch('pyMultiwfn.job.subprocess.Popen') as mock_popen:
            mock_process = MagicMock()
            mock_process.communicate.return_value = (b"output", b"")
            mock_process.returncode = 0
            mock_popen.return_value = mock_process
            
            result = job.run()
        
        assert result is not None
        assert result.returncode == 0
        assert "output" in result.stdout  # stdout is decoded to string
    
    def test_run_captures_stderr(self, mock_wfn_file, mock_executable):
        """Run captures stderr on success."""
        config = MultiwfnConfig(exe_path=mock_executable)
        job = MultiwfnJob(mock_wfn_file, config=config)
        job.add_menu_sequence(menu.hirshfeld_charge)
        
        with patch('pyMultiwfn.job.subprocess.Popen') as mock_popen:
            mock_process = MagicMock()
            mock_process.communicate.return_value = (b"output", b"warning message")
            mock_process.returncode = 0  # Success with warnings
            mock_popen.return_value = mock_process
            
            result = job.run()
        
        assert "warning" in result.stderr
        assert result.returncode == 0
    
    def test_run_timeout(self, mock_wfn_file, temp_dir):
        """Run raises error on timeout."""
        import subprocess as sp
        
        # Create a fake executable that sleeps
        if platform.system() == "Windows":
            fake_exe = temp_dir / "fake_multiwfn.bat"
            fake_exe.write_text("@ping 127.0.0.1 -n 10 > nul")
        else:
            fake_exe = temp_dir / "fake_multiwfn"
            fake_exe.write_text("#!/bin/bash\nsleep 10")
            fake_exe.chmod(0o755)
        
        config = MultiwfnConfig(exe_path=fake_exe, timeout=1)
        job = MultiwfnJob(mock_wfn_file, config=config)
        job.add_menu_sequence(menu.hirshfeld_charge)
        
        with pytest.raises(MultiwfnError, match="timed out"):
            job.run()
    
    def test_run_stores_result(self, mock_wfn_file, mock_executable):
        """Run stores result on job."""
        config = MultiwfnConfig(exe_path=mock_executable)
        job = MultiwfnJob(mock_wfn_file, config=config)
        job.add_menu_sequence(menu.hirshfeld_charge)
        
        with patch('pyMultiwfn.job.subprocess.Popen') as mock_popen:
            mock_process = MagicMock()
            mock_process.communicate.return_value = (b"output", b"")
            mock_process.returncode = 0
            mock_popen.return_value = mock_process
            
            result = job.run()
        
        assert job.result is result


class TestMultiwfnJobWithRealFiles:
    """Tests using real wavefunction files."""
    
    @pytest.mark.real_files
    def test_job_with_real_wfx(self, real_wfx_file, mock_executable):
        """Job accepts real .wfx file."""
        config = MultiwfnConfig(exe_path=mock_executable)
        job = MultiwfnJob(real_wfx_file, config=config)
        
        assert job.input_file == real_wfx_file
        assert job.input_file.exists()
    
    @pytest.mark.real_files
    def test_job_with_real_molden(self, real_molden_file, mock_executable):
        """Job accepts real .molden file."""
        config = MultiwfnConfig(exe_path=mock_executable)
        job = MultiwfnJob(real_molden_file, config=config)
        
        assert job.input_file == real_molden_file
        assert job.input_file.exists()
    
    @pytest.mark.real_files
    def test_builder_with_real_file(self, real_wfx_file, mock_executable):
        """Builder works with real file."""
        config = MultiwfnConfig(exe_path=mock_executable)
        job = (
            MultiwfnJob.builder(real_wfx_file)
            .with_config(config)
            .with_menu_sequence(menu.hirshfeld_charge)
            .build()
        )
        
        assert job.input_file == real_wfx_file