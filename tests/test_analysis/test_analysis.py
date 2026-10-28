"""Tests for pymultiwfn.analysis.analysis — MultiwfnAnalysis."""

import json
from pathlib import Path

import pytest

from pymultiwfn.analysis.analysis import MultiwfnAnalysis
from pymultiwfn.api.exceptions import MultiwfnError
from pymultiwfn.api.job import MultiwfnJob
from pymultiwfn.api.multiwfn import Multiwfn
from pymultiwfn.enums.analyses import AnalysisClasses
from pymultiwfn.enums.menu import Menu


@pytest.fixture
def multiwfn(mock_executable: Path) -> Multiwfn:
    return Multiwfn(exe_path=mock_executable)


class TestOutputDir:
    """Tests for MultiwfnAnalysis._output_dir()."""

    def test_named_after_full_input_file_name(self, temp_dir: Path) -> None:
        analysis = MultiwfnAnalysis(input_file=temp_dir / "coord.molden")
        output_dir = analysis._output_dir(temp_dir)
        assert output_dir == temp_dir / "coord.molden.output"

    def test_distinguishes_same_stem_different_extension(
        self, temp_dir: Path
    ) -> None:
        wfn_analysis = MultiwfnAnalysis(input_file=temp_dir / "mol.wfn")
        fch_analysis = MultiwfnAnalysis(input_file=temp_dir / "mol.fch")
        assert wfn_analysis._output_dir(temp_dir) != fch_analysis._output_dir(
            temp_dir
        )
        assert wfn_analysis._output_dir(temp_dir).name == "mol.wfn.output"
        assert fch_analysis._output_dir(temp_dir).name == "mol.fch.output"

    def test_defaults_to_cwd_when_work_dir_none(self, temp_dir: Path) -> None:
        analysis = MultiwfnAnalysis(input_file=temp_dir / "mock.wfn")
        assert analysis._output_dir(None) == Path.cwd() / "mock.wfn.output"


class TestAddMenu:
    """Tests for MultiwfnAnalysis.add_menu()."""

    def test_add_single_menu(self, temp_dir: Path) -> None:
        analysis = MultiwfnAnalysis(input_file=temp_dir / "mock.wfn")
        analysis.add_menu(Menu.HIRSHFELD_CHARGE)
        assert analysis.analyses == [Menu.HIRSHFELD_CHARGE]

    def test_add_list_of_menus(self, temp_dir: Path) -> None:
        analysis = MultiwfnAnalysis(input_file=temp_dir / "mock.wfn")
        analysis.add_menu([Menu.HIRSHFELD_CHARGE, Menu.MAYER_BOND_ORDER])
        assert analysis.analyses == [
            Menu.HIRSHFELD_CHARGE,
            Menu.MAYER_BOND_ORDER,
        ]

    def test_add_analysis_class_expands_to_its_menus(self, temp_dir: Path) -> None:
        analysis = MultiwfnAnalysis(input_file=temp_dir / "mock.wfn")
        analysis.add_menu(AnalysisClasses.CHARGES)
        assert analysis.analyses == AnalysisClasses.CHARGES.value

    def test_add_invalid_type_raises(self, temp_dir: Path) -> None:
        analysis = MultiwfnAnalysis(input_file=temp_dir / "mock.wfn")
        with pytest.raises(TypeError):
            analysis.add_menu("not a menu")  # type: ignore[arg-type]

    def test_constructor_accepts_initial_analyses(self, temp_dir: Path) -> None:
        analysis = MultiwfnAnalysis(
            input_file=temp_dir / "mock.wfn",
            analyses=Menu.HIRSHFELD_CHARGE,
        )
        assert analysis.analyses == [Menu.HIRSHFELD_CHARGE]


class TestJsonPath:
    """Tests for MultiwfnAnalysis.json_path."""

    def test_defaults_to_none(self, temp_dir: Path) -> None:
        analysis = MultiwfnAnalysis(input_file=temp_dir / "mock.wfn")
        assert analysis.json_path is None

    def test_explicit_json_path_is_respected(self, temp_dir: Path) -> None:
        json_path = temp_dir / "custom.json"
        analysis = MultiwfnAnalysis(
            input_file=temp_dir / "mock.wfn", json_path=json_path
        )
        assert analysis.json_path == json_path


class TestRun:
    """Integration tests for MultiwfnAnalysis.run() against a mock exe."""

    def test_run_creates_output_dir_named_after_full_input_file_name(
        self,
        mock_wfn_file: Path,
        multiwfn: Multiwfn,
        temp_dir: Path,
    ) -> None:
        analysis = MultiwfnAnalysis(
            input_file=mock_wfn_file,
            analyses=Menu.HIRSHFELD_CHARGE,
        )
        analysis.run(multiwfn=multiwfn, work_dir=temp_dir)

        expected_output_dir = temp_dir / f"{mock_wfn_file.name}.output"
        assert expected_output_dir.is_dir()

    def test_run_writes_stdout_named_after_menu_member(
        self,
        mock_wfn_file: Path,
        multiwfn: Multiwfn,
        temp_dir: Path,
    ) -> None:
        analysis = MultiwfnAnalysis(
            input_file=mock_wfn_file,
            analyses=Menu.HIRSHFELD_CHARGE,
        )
        analysis.run(multiwfn=multiwfn, work_dir=temp_dir)

        output_dir = temp_dir / f"{mock_wfn_file.name}.output"
        stdout_file = output_dir / "HIRSHFELD_CHARGE.txt"
        assert stdout_file.is_file()
        assert "mock" in stdout_file.read_text(encoding="utf-8")

    def test_run_writes_one_stdout_file_per_queued_analysis(
        self,
        mock_wfn_file: Path,
        multiwfn: Multiwfn,
        temp_dir: Path,
    ) -> None:
        analysis = MultiwfnAnalysis(
            input_file=mock_wfn_file,
            analyses=[Menu.HIRSHFELD_CHARGE, Menu.MAYER_BOND_ORDER],
        )
        analysis.run(multiwfn=multiwfn, work_dir=temp_dir)

        output_dir = temp_dir / f"{mock_wfn_file.name}.output"
        assert (output_dir / "HIRSHFELD_CHARGE.txt").is_file()
        assert (output_dir / "MAYER_BOND_ORDER.txt").is_file()

    def test_run_auto_json_path_lives_inside_output_dir(
        self,
        mock_wfn_file: Path,
        multiwfn: Multiwfn,
        temp_dir: Path,
    ) -> None:
        analysis = MultiwfnAnalysis(
            input_file=mock_wfn_file,
            analyses=Menu.HIRSHFELD_CHARGE,
        )
        analysis.run(multiwfn=multiwfn, work_dir=temp_dir)

        output_dir = temp_dir / f"{mock_wfn_file.name}.output"
        json_path = output_dir / f"{mock_wfn_file.stem}.json"
        assert json_path.is_file()

        data = json.loads(json_path.read_text(encoding="utf-8"))
        assert "generated_files" in data
        assert "HIRSHFELD_CHARGE.txt" in data["generated_files"]

    def test_run_respects_explicit_json_path(
        self,
        mock_wfn_file: Path,
        multiwfn: Multiwfn,
        temp_dir: Path,
    ) -> None:
        custom_json = temp_dir / "elsewhere" / "results.json"
        analysis = MultiwfnAnalysis(
            input_file=mock_wfn_file,
            analyses=Menu.HIRSHFELD_CHARGE,
            json_path=custom_json,
        )
        analysis.run(multiwfn=multiwfn, work_dir=temp_dir)

        assert custom_json.is_file()
        output_dir = temp_dir / f"{mock_wfn_file.name}.output"
        assert not (output_dir / f"{mock_wfn_file.stem}.json").exists()

    def test_run_populates_in_memory_results_and_jobs(
        self,
        mock_wfn_file: Path,
        multiwfn: Multiwfn,
        temp_dir: Path,
    ) -> None:
        analysis = MultiwfnAnalysis(
            input_file=mock_wfn_file,
            analyses=Menu.HIRSHFELD_CHARGE,
        )
        analysis.run(multiwfn=multiwfn, work_dir=temp_dir)

        assert len(analysis.results) == 1
        assert len(analysis.jobs) == 1
        assert analysis.jobs[0].executed

    def test_run_uses_cached_result_and_skips_rerun(
        self,
        mock_wfn_file: Path,
        multiwfn: Multiwfn,
        temp_dir: Path,
    ) -> None:
        output_dir = temp_dir / f"{mock_wfn_file.name}.output"
        output_dir.mkdir(parents=True)
        json_path = output_dir / f"{mock_wfn_file.stem}.json"
        json_path.write_text(
            json.dumps(
                {
                    "input_file": str(mock_wfn_file),
                    "analyses": {
                        "HIRSHFELD_CHARGE": {
                            "parsed": {},
                            "timestamp": "2024-01-01T00:00:00",
                        }
                    },
                }
            ),
            encoding="utf-8",
        )

        analysis = MultiwfnAnalysis(
            input_file=mock_wfn_file,
            analyses=Menu.HIRSHFELD_CHARGE,
            cached=True,
        )
        analysis.run(multiwfn=multiwfn, work_dir=temp_dir)

        # The cached branch appends a lightweight result but never
        # creates/runs a MultiwfnJob for it.
        assert len(analysis.results) == 1
        assert len(analysis.jobs) == 0

    def test_create_and_run_propagates_job_error(
        self,
        mock_wfn_file: Path,
        multiwfn: Multiwfn,
        temp_dir: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        def _boom(self: MultiwfnJob) -> MultiwfnJob:
            raise MultiwfnError("simulated Multiwfn failure")

        monkeypatch.setattr(MultiwfnJob, "run", _boom)

        analysis = MultiwfnAnalysis(
            input_file=mock_wfn_file,
            analyses=Menu.HIRSHFELD_CHARGE,
        )
        with pytest.raises(MultiwfnError):
            analysis.run(multiwfn=multiwfn, work_dir=temp_dir)

        # Even though the job failed, a result/job pair is still
        # recorded before the error is re-raised.
        assert len(analysis.results) == 1
        assert len(analysis.jobs) == 1


class TestHasResult:
    """Tests for MultiwfnAnalysis._has_result()."""

    def test_false_when_no_results(self, temp_dir: Path) -> None:
        analysis = MultiwfnAnalysis(input_file=temp_dir / "mock.wfn")
        assert analysis._has_result(Menu.HIRSHFELD_CHARGE) is False

    def test_true_after_matching_result_recorded(
        self,
        mock_wfn_file: Path,
        multiwfn: Multiwfn,
        temp_dir: Path,
    ) -> None:
        analysis = MultiwfnAnalysis(
            input_file=mock_wfn_file,
            analyses=Menu.HIRSHFELD_CHARGE,
        )
        analysis.run(multiwfn=multiwfn, work_dir=temp_dir)
        assert analysis._has_result(Menu.HIRSHFELD_CHARGE) is True
        assert analysis._has_result(Menu.MAYER_BOND_ORDER) is False


class TestJsonPathPropagation:
    """Tests for json_path setter propagation to an initialised store."""

    def test_setter_propagates_to_existing_store(self, temp_dir: Path) -> None:
        analysis = MultiwfnAnalysis(input_file=temp_dir / "mock.wfn")
        analysis._get_store(temp_dir)
        assert analysis._store is not None

        new_json = temp_dir / "new_results.json"
        analysis.json_path = new_json

        assert analysis._store.json_path == new_json

    def test_setter_before_store_init_has_no_store_to_update(
        self, temp_dir: Path
    ) -> None:
        analysis = MultiwfnAnalysis(input_file=temp_dir / "mock.wfn")
        assert analysis._store is None
        analysis.json_path = temp_dir / "new_results.json"
        assert analysis._store is None


class TestConvenienceMenuMethods:
    """Each _add_*_menus() helper adds exactly its AnalysisClasses."""

    @pytest.mark.parametrize(
        ("method_name", "category"),
        [
            ("_add_charges_menus", AnalysisClasses.CHARGES),
            ("_add_bond_orders_menus", AnalysisClasses.BOND_ORDERS),
            ("_add_topology_menus", AnalysisClasses.TOPOLOGY),
            (
                "_add_weak_interactions_menus",
                AnalysisClasses.WEAK_INTERACTIONS,
            ),
            ("_add_spectra_menus", AnalysisClasses.SPECTRA),
            ("_add_surfaces_menus", AnalysisClasses.SURFACES),
            ("_add_aromaticity_menus", AnalysisClasses.AROMATICITY),
            ("_add_cdft_menus", AnalysisClasses.CDFT),
            ("_add_dos_menus", AnalysisClasses.DOS),
            ("_add_basin_menus", AnalysisClasses.BASIN),
            ("_add_excitation_menus", AnalysisClasses.EXCITATION),
            ("_add_cubes_menus", AnalysisClasses.CUBES),
            (
                "_add_orbital_composition_menus",
                AnalysisClasses.ORBITAL_COMPOSITION,
            ),
            (
                "_add_orbital_localization_menus",
                AnalysisClasses.ORBITAL_LOCALIZATION,
            ),
            ("_add_fuzzy_space_menus", AnalysisClasses.FUZZY_SPACE),
            ("_add_eda_menus", AnalysisClasses.EDA),
            ("_add_polarizability_menus", AnalysisClasses.POLARIZABILITY),
            ("_add_wavefunction_menus", AnalysisClasses.WAVEFUNCTION),
            ("_add_line_plots_menus", AnalysisClasses.LINE_PLOTS),
            ("_add_plane_maps_menus", AnalysisClasses.PLANE_MAPS),
        ],
    )
    def test_convenience_method_adds_expected_category(
        self,
        temp_dir: Path,
        method_name: str,
        category: AnalysisClasses,
    ) -> None:
        analysis = MultiwfnAnalysis(input_file=temp_dir / "mock.wfn")
        getattr(analysis, method_name)()
        assert analysis.analyses == category.value
