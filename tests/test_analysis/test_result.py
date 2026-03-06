"""Tests for pymultiwfn.analysis.result — MultiwfnResult."""

import pytest

from pymultiwfn.analysis.result import MultiwfnResult


class TestMultiwfnResult:
    """Tests for the MultiwfnResult class."""

    def test_create(self) -> None:
        r = MultiwfnResult(type="charges", value={"1": 0.05})
        assert r.type == "charges"
        assert r.value == {"1": 0.05}

    def test_parse_charges(self) -> None:
        r = MultiwfnResult(type="charges", value=None)
        stdout = "Final atomic charges:\nAtom    1(C ):    -0.05230000\n"
        charges = r.parse_charges(stdout)
        assert charges[1] == pytest.approx(-0.0523)

    def test_parse_bond_orders(self) -> None:
        r = MultiwfnResult(type="bond_orders", value=None)
        stdout = "#    1:         1(C )    2(N )    1.45230000\n"
        bos = r.parse_bond_orders(stdout)
        assert bos[(1, 2)] == pytest.approx(1.4523)

    def test_parse_critical_points(self) -> None:
        r = MultiwfnResult(type="topology", value=None)
        stdout = (
            "CP  1 (3,-3)\n"
            "Position (Bohr):   0.0   0.0   0.0\n"
            "Density of all electrons:   0.298\n"
        )
        cps = r.parse_critical_points(stdout)
        assert len(cps) == 1

    def test_parse_spectrum(self) -> None:
        r = MultiwfnResult(type="spectrum", value=None)
        stdout = "  500.0 cm^-1  Intensity:  12.34\n"
        sp = r.parse_spectrum(stdout)
        assert len(sp["frequencies"]) == 1

    def test_parse_charges_empty(self) -> None:
        r = MultiwfnResult(type="charges", value=None)
        assert r.parse_charges("") == {}

    def test_parse_bond_orders_empty(self) -> None:
        r = MultiwfnResult(type="bond_orders", value=None)
        assert r.parse_bond_orders("") == {}
