"""Pytest configuration and shared fixtures for pymultiwfn tests."""

import platform
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))


# =============================================================================
# Path Fixtures
# =============================================================================


@pytest.fixture(scope="session")
def project_root() -> Path:
    return Path(__file__).parent.parent


@pytest.fixture(scope="session")
def test_data_dir(project_root: Path) -> Path:
    return project_root / "tests" / "test_data"


@pytest.fixture
def temp_dir(tmp_path: Path) -> Path:
    return tmp_path


@pytest.fixture
def mock_wfn_file(temp_dir: Path) -> Path:
    path = temp_dir / "mock.wfn"
    path.write_text("mock wavefunction content")
    return path


@pytest.fixture
def mock_executable(temp_dir: Path) -> Path:
    if platform.system() == "Windows":
        exe_path = temp_dir / "Multiwfn.exe"
        exe_path.write_text("@echo mock")
    else:
        exe_path = temp_dir / "Multiwfn"
        exe_path.write_text("#!/bin/bash\necho 'mock'")
        exe_path.chmod(0o755)
    return exe_path


# =============================================================================
# Sample Output Fixtures
# =============================================================================


@pytest.fixture
def sample_hirshfeld_output() -> str:
    return (
        "Final atomic charges:\n"
        "Atom    1(C ):    -0.05230000\n"
        "Atom    2(C ):    -0.08920000\n"
        "Atom    3(H ):     0.05340000\n"
        "Atom    4(H ):     0.05210000\n"
        "Atom    5(O ):    -0.32450000\n"
        "Atom    6(H ):     0.18230000\n"
    )


@pytest.fixture
def sample_mayer_output() -> str:
    return (
        "Mayer bond order analysis\n"
        "Bond order list:\n"
        "#    1:         1(C )    2(C )    1.45230000\n"
        "#    2:         1(C )    3(H )    0.92340000\n"
        "#    3:         2(C )    4(H )    0.95670000\n"
        "#    4:         2(C )    5(O )    1.89340000\n"
        "#    5:         5(O )    6(N )    0.82340000\n"
    )


@pytest.fixture
def sample_topology_output() -> str:
    return (
        "Topology analysis\n"
        "CP  1 (3,-3)\n"
        "Position (Bohr):   0.000000   0.000000   0.000000\n"
        "Density of all electrons:   0.29876500\n"
        "Laplacian of electron density:  -1.12345600\n"
        "Ellipticity:   0.00000000\n"
        "\n\n\n\n\n"
        "CP  2 (3,-1)\n"
        "Position (Bohr):   1.234567   0.567890   0.000000\n"
        "Density of all electrons:   0.25432100\n"
        "Laplacian of electron density:  -0.87654300\n"
        "Ellipticity:   0.04523000\n"
        "\n\n\n\n\n"
        "CP  3 (3,+1)\n"
        "Position (Bohr):   0.500000   0.500000   0.500000\n"
        "Density of all electrons:   0.01234500\n"
        "Laplacian of electron density:   0.05678900\n"
        "\n\n\n\n\n"
        "CP  4 (3,+3)\n"
        "Position (Bohr):   1.000000   1.000000   1.000000\n"
        "Density of all electrons:   0.00123400\n"
        "Laplacian of electron density:   0.00987600\n"
    )


@pytest.fixture
def sample_spectrum_output() -> str:
    return (
        "IR spectrum data\n"
        "  500.0 cm^-1  Intensity:  12.34\n"
        "  800.0 cm^-1  Intensity:  45.67\n"
        " 1200.0 cm^-1  Intensity:  89.01\n"
        " 1500.0 cm^-1  Intensity:  23.45\n"
        " 3000.0 cm^-1  Intensity: 156.78\n"
    )


# =============================================================================
# Markers
# =============================================================================


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line("markers", "slow: marks tests as slow")
    config.addinivalue_line("markers", "integration: marks integration tests")
    config.addinivalue_line(
        "markers", "requires_multiwfn: marks tests requiring Multiwfn"
    )
