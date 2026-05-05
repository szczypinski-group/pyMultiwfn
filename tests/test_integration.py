"""End-to-end integration tests aligned with current parser schemas.

These tests require the Multiwfn executable and test data files.
They are automatically skipped when either is unavailable.

Usage::

    pytest tests/test_integration.py -v

"""

from __future__ import annotations

from pathlib import Path

import pytest

from pymultiwfn.api.multiwfn import Multiwfn

# =========================================================================
# Discovery
# =========================================================================

_TESTS_ROOT = Path(__file__).resolve().parent
_TEST_DATA = _TESTS_ROOT / "test_data"
_WFN_EXTENSIONS = (".wfn", ".wfx", ".fchk", ".molden", ".mwfn")


def _find_wfn_files() -> list[Path]:
    if not _TEST_DATA.is_dir():
        return []
    files: list[Path] = []
    for ext in _WFN_EXTENSIONS:
        files.extend(sorted(_TEST_DATA.glob(f"*{ext}")))
    return files


def _find_executable() -> Path | None:
    try:
        return Multiwfn().exe_path
    except Exception:
        return None


_exe_path = _find_executable()
_wfn_files = _find_wfn_files()

# Skip the entire module if prerequisites are missing
pytestmark = [
    pytest.mark.skipif(
        _exe_path is None, reason="Multiwfn executable not found"
    ),
    pytest.mark.skipif(
        len(_wfn_files) == 0,
        reason="No wavefunction files in tests/test_data/",
    ),
    pytest.mark.integration,
]


# =========================================================================
# Fixtures
# =========================================================================


@pytest.fixture(scope="session")
def exe() -> Multiwfn:
    if _exe_path is None:
        pytest.skip("Multiwfn executable not available")
    return Multiwfn(exe_path=_exe_path)


@pytest.fixture(scope="session")
def wfn() -> Path:
    if not _wfn_files:
        pytest.skip("No wavefunction files")
    return _wfn_files[0]


@pytest.fixture()
def wd(tmp_path: Path) -> Path:
    d = tmp_path / "work"
    d.mkdir()
    return d


# =========================================================================
# Placeholder test — verifies the integration setup works
# =========================================================================


class TestIntegrationSetup:
    """Verify integration test infrastructure is correctly set up."""

    def test_executable_exists(self, exe: Multiwfn) -> None:
        assert exe.exe_path.exists()

    def test_wavefunction_file_exists(self, wfn: Path) -> None:
        assert wfn.exists()
