"""Tests for pymultiwfn.analysis.result — MultiwfnResult."""

import pytest

from pymultiwfn.analysis.result import MultiwfnResult
from pymultiwfn.enums.menu import Menu


class TestMultiwfnResult:
    """Tests for the MultiwfnResult class."""

    def test_create(self) -> None:
        r = MultiwfnResult(analysis=Menu.HIRSHFELD_CHARGE)
        assert r.analysis == Menu.HIRSHFELD_CHARGE
        assert r.result == []

    def test_parse_charges(self) -> None:
        r = MultiwfnResult(analysis=Menu.HIRSHFELD_CHARGE)
        stdout = "Final atomic charges:\nAtom    1(C ):    -0.05230000\n"
        r.parse(stdout)
        # Result should contain Charge dataclass(es)
        assert len(r.result) > 0
        charge_results = [x for x in r.result if hasattr(x, "charge")]
        assert charge_results[0].charge == pytest.approx(-0.0523)

    def test_parse_bond_orders(self) -> None:
        r = MultiwfnResult(analysis=Menu.MAYER_BOND_ORDER)
        stdout = "#    1:         1(C )    2(N )    1.45230000\n"
        r.parse(stdout)
        bo_results = [x for x in r.result if hasattr(x, "bond_order")]
        assert bo_results[0].bond_order == pytest.approx(1.4523)

    def test_parse_critical_points(self) -> None:
        r = MultiwfnResult(analysis=Menu.TOPOLOGY_SEARCH_CPS)
        stdout = (
            "CP  1 (3,-3)\n"
            "Position (Bohr):   0.0   0.0   0.0\n"
            "Density of all electrons:   0.298\n"
        )
        r.parse(stdout)
        cp_results = [x for x in r.result if hasattr(x, "rho")]
        assert len(cp_results) == 1

    def test_parse_spectrum(self) -> None:
        r = MultiwfnResult(analysis=Menu.PLOT_IR_SPECTRUM)
        stdout = "  500.0 cm^-1  Intensity:  12.34\n"
        r.parse(stdout)
        spectrum_results = [x for x in r.result if hasattr(x, "frequencies")]
        assert len(spectrum_results) == 1
        assert len(spectrum_results[0].frequencies) == 1

    def test_parse_charges_empty(self) -> None:
        r = MultiwfnResult(analysis=Menu.HIRSHFELD_CHARGE)
        r.parse("")
        assert r.result == []

    def test_parse_bond_orders_empty(self) -> None:
        r = MultiwfnResult(analysis=Menu.MAYER_BOND_ORDER)
        r.parse("")
        assert r.result == []
