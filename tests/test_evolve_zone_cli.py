"""evolve-zone on large networks and in explosive burning.

Given a complete JINA database, ``nucnetpy evolve-zone`` must let the user cut
the network and set the tolerances and the step control, report how far a long
integration has got, and say what to try when it advances too slowly to finish.
When the run ends it prints a short summary and saves the history.
"""
from pathlib import Path

import numpy as np
import pytest

import nucnetpy.cli as cli
from nucnetpy import evolve_zone, read_xml, time_grid
from nucnetpy.cli import main

GOLDEN = str(Path(__file__).parent / "golden" / "golden_network.xml")
SHORT = ["evolve-zone", GOLDEN, "--t1", "1e-3", "--steps", "3"]


@pytest.fixture(autouse=True)
def _run_in_a_scratch_folder(tmp_path, monkeypatch):
    """evolve-zone saves every run in a new folder of the current directory."""
    monkeypatch.chdir(tmp_path)


@pytest.mark.parametrize("method", ["bdf", "adaptive_implicit_euler"])
def test_progress_hook_follows_the_integration_and_changes_nothing(method):
    net = read_xml(GOLDEN)
    times = time_grid(0.0, 1.0, 5)
    seen = []
    watched = evolve_zone(net, net.zone(0), times, method=method, progress=seen.append)
    plain = evolve_zone(net, net.zone(0), times, method=method)
    assert seen and max(seen) == pytest.approx(1.0)
    assert np.array_equal(watched.y, plain.y)


def test_cut_options_restrict_the_network(capsys):
    assert main(SHORT + ["--zmax", "6", "--progress", "0", "--output", "run"]) == 0
    assert capsys.readouterr().out.startswith("evolve-zone: finished in")
    names = set(read_xml("run/output.xml").zones[-1].abundances)
    assert names and names <= {"n", "h1", "he4", "be8", "c12"}


def test_progress_lines_go_to_standard_error(capsys):
    main(SHORT + ["--progress", "1e-9"])
    captured = capsys.readouterr()
    assert "evolve-zone: t = " in captured.err and " dt = " in captured.err
    assert " dt = " not in captured.out


def test_progress_can_be_turned_off(capsys):
    main(SHORT + ["--progress", "0"])
    assert " dt = " not in capsys.readouterr().err


def test_a_large_network_is_announced(monkeypatch, capsys):
    monkeypatch.setattr(cli, "LARGE_NETWORK", 3)
    main(SHORT + ["--progress", "0"])
    err = capsys.readouterr().err
    assert "note: the network has" in err
    assert cli.EXPLOSIVE_STEP_CONTROL in err and "--zmax" in err


def test_a_small_network_is_not_announced(capsys):
    main(SHORT + ["--progress", "0"])
    assert "note" not in capsys.readouterr().err


def test_tolerances_reach_the_integrator(monkeypatch):
    seen = {}
    real = cli.evolve_zone

    def spy(*args, **kwargs):
        seen.update(kwargs)
        return real(*args, **kwargs)

    monkeypatch.setattr(cli, "evolve_zone", spy)
    main(SHORT + ["--rtol", "1e-8", "--atol", "1e-12", "--progress", "0"])
    assert seen["rtol"] == 1e-8 and seen["atol"] == 1e-12


def test_step_control_reaches_the_integrator(monkeypatch):
    seen = {}
    real = cli.evolve_zone

    def spy(*args, **kwargs):
        seen.update(kwargs)
        return real(*args, **kwargs)

    monkeypatch.setattr(cli, "evolve_zone", spy)
    assert main(SHORT + ["--method", "adaptive_implicit_euler", "--max-change", "0.05",
                         "--min-mass-fraction", "1e-10", "--progress", "0"]) == 0
    assert seen["max_change"] == 0.05 and seen["min_mass_fraction"] == 1e-10


def test_step_control_defaults_are_those_of_evolve_zone():
    import inspect
    defaults = inspect.signature(evolve_zone).parameters
    args = cli.build_parser().parse_args(SHORT)
    assert args.max_change == defaults["max_change"].default
    assert args.min_mass_fraction == defaults["min_mass_fraction"].default


@pytest.mark.parametrize("method, advice", [("bdf", cli.EXPLOSIVE_STEP_CONTROL),
                                            ("adaptive_implicit_euler", "--max-change")])
def test_an_integration_that_will_not_finish_gets_advice(capsys, method, advice):
    args = cli.build_parser().parse_args(SHORT + ["--method", method])
    progress = cli._Progress(time_grid(0.0, 1.0, 5), 0.0, args)
    progress(1e-12)
    progress(2e-12)
    err = capsys.readouterr().err
    assert err.count("hint:") == 1 and advice in err


def test_an_integration_that_advances_gets_no_advice(capsys):
    args = cli.build_parser().parse_args(SHORT)
    progress = cli._Progress(time_grid(0.0, 1.0, 5), 0.0, args)
    progress(0.3)
    progress(0.4)
    assert "hint" not in capsys.readouterr().err


def test_a_transient_whose_steps_grow_gets_no_advice(capsys):
    # Far from the end, but the reached time grows geometrically: the
    # integration is getting through a fast transient, not crawling.
    args = cli.build_parser().parse_args(SHORT + ["--method", "adaptive_implicit_euler"])
    progress = cli._Progress(time_grid(0.0, 1.0, 5), 0.0, args)
    for t in (1e-13, 1e-11, 1e-9):
        progress(t)
    assert "hint" not in capsys.readouterr().err


def test_advice_does_not_repeat_the_settings_in_use(capsys):
    args = cli.build_parser().parse_args(SHORT + ["--method", "adaptive_implicit_euler",
                                                  "--max-change", "0.05",
                                                  "--min-mass-fraction", "1e-10"])
    progress = cli._Progress(time_grid(0.0, 1.0, 5), 0.0, args)
    progress(1e-12)
    progress(2e-12)
    err = capsys.readouterr().err
    assert "hint:" in err and "--zmax" in err
    assert "--max-change" not in err and "--min-mass-fraction" not in err


# --- run conditions from the zone, as in NucNet Tools ------------------------

def _zone_file(tmp_path, **props):
    """The golden network with its zone's properties replaced by ``props``."""
    from nucnetpy import write_xml
    net = read_xml(GOLDEN)
    zone = net.zone(0)
    zone.properties.clear()
    zone.optional_properties.clear()
    zone.optional_properties.update({k: str(v) for k, v in props.items()})
    path = tmp_path / "zone_network.xml"
    write_xml(net, path)
    return str(path)


def _spy(monkeypatch):
    seen = {}
    real = cli.evolve_zone

    def spy(net, zone, times, **kwargs):
        seen.update(kwargs, times=times)
        return real(net, zone, times, **kwargs)

    monkeypatch.setattr(cli, "evolve_zone", spy)
    return seen


def test_conditions_come_from_the_zone(tmp_path, monkeypatch, capsys):
    path = _zone_file(tmp_path, t9_0=0.4, rho_0=2e4, tend=1e-3, dt=1e-7)
    seen = _spy(monkeypatch)
    assert main(["evolve-zone", path, "--steps", "3", "--progress", "0"]) == 0
    assert seen["times"][-1] == pytest.approx(1e-3)
    assert seen["first_step"] == pytest.approx(1e-7)
    assert seen["thermo"](0.0, {}) == (0.4, 2e4)
    err = capsys.readouterr().err
    for name in ("zone property t9_0", "zone property rho_0", "zone property tend",
                 "zone property dt"):
        assert name in err


def test_options_override_the_zone(tmp_path, monkeypatch, capsys):
    path = _zone_file(tmp_path, t9_0=0.4, rho_0=2e4, tend=1e-3, dt=1e-7)
    seen = _spy(monkeypatch)
    assert main(["evolve-zone", path, "--steps", "3", "--progress", "0", "--t9", "0.3",
                 "--rho", "1e3", "--t1", "2e-3", "--dt", "1e-9"]) == 0
    assert seen["times"][-1] == pytest.approx(2e-3)
    assert seen["first_step"] == pytest.approx(1e-9)
    assert seen["thermo"](0.0, {}) == (0.3, 1e3)
    err = capsys.readouterr().err
    assert "(--t9)" in err and "(--t1)" in err and "(--dt)" in err


def test_without_end_time_or_step_the_defaults_are_said(tmp_path, monkeypatch, capsys):
    path = _zone_file(tmp_path, t9=0.4, rho=2e4)
    seen = _spy(monkeypatch)
    assert main(["evolve-zone", path, "--steps", "3", "--progress", "0"]) == 0
    assert seen["times"][-1] == 1.0 and seen["first_step"] is None
    err = capsys.readouterr().err
    assert "(default)" in err and "chosen by the integrator" in err


def test_a_run_without_temperature_or_density_is_refused(tmp_path, capsys):
    path = _zone_file(tmp_path, tend=1e-3)
    assert main(["evolve-zone", path, "--progress", "0"]) == 1
    err = capsys.readouterr().err
    assert "no temperature" in err and "--t9" in err
    assert "no density" in err and "--rho" in err


@pytest.mark.parametrize("method", ["bdf", "adaptive_implicit_euler"])
def test_first_step_reaches_the_integrator(monkeypatch, method):
    import scipy.integrate
    net = read_xml(GOLDEN)
    times = time_grid(0.0, 1e-2, 4)
    seen = {}
    real = scipy.integrate.solve_ivp

    def spy(*args, **kwargs):
        seen.update(kwargs)
        return real(*args, **kwargs)

    monkeypatch.setattr(scipy.integrate, "solve_ivp", spy)
    reached = []
    result = evolve_zone(net, net.zone(0), times, method=method, first_step=1e-9,
                         progress=reached.append)
    assert result.success
    if method == "bdf":
        assert seen["first_step"] == 1e-9
    else:
        assert reached[0] == pytest.approx(1e-9)


def test_a_first_step_longer_than_the_run_is_cut():
    net = read_xml(GOLDEN)
    result = evolve_zone(net, net.zone(0), time_grid(0.0, 1e-6, 3), first_step=1.0)
    assert result.success


def test_zone_reads_end_time_and_time_step():
    from nucnetpy import Zone
    zone = Zone(optional_properties={"tend": "3.15e7", "dt": "1e-5", "t9_0": "0.2"},
                properties={"dt": "1e-4"})
    assert zone.end_time() == 3.15e7
    assert zone.time_step() == 1e-4          # a property on <zone> wins
    assert zone.find_property(("t9", "t9_0")) == ("t9_0", "0.2")
    assert Zone().end_time() is None and Zone().end_time(5.0) == 5.0


# --- what a run reports while it goes, and what it leaves behind -------------

@pytest.mark.parametrize("method", ["bdf", "radau", "lsoda", "adaptive_implicit_euler", "rk4"])
def test_progress_is_called_once_per_accepted_step(method):
    net = read_xml(GOLDEN)
    end = 1e-6 if method == "rk4" else 1.0
    seen = []
    result = evolve_zone(net, net.zone(0), time_grid(0.0, end, 5), method=method,
                         progress=seen.append)
    assert result.steps == len(seen) > 0
    assert all(a < b for a, b in zip(seen, seen[1:]))
    assert seen[-1] == pytest.approx(end)


def test_progress_reports_the_step_count_and_size(capsys):
    args = cli.build_parser().parse_args(SHORT)
    progress = cli._Progress(time_grid(0.0, 1.0, 5), 0.0, args)
    progress(1e-3)
    progress(3e-3)
    assert "t = 0.003 s of 1 s | step 2, dt = 0.002 s | 0 of 4 output times" in \
        capsys.readouterr().err


def test_a_run_is_saved_in_a_new_folder(tmp_path, capsys):
    assert main(SHORT + ["--progress", "0"]) == 0
    folders = [p for p in tmp_path.iterdir() if p.name.startswith("evolve-zone_")]
    assert len(folders) == 1
    saved = read_xml(folders[0] / "output.xml")
    assert [z.label[0] for z in saved.zones] == ["0", "1", "2"]
    times = [float(z.optional_properties["time"]) for z in saved.zones]
    assert times == pytest.approx([0.0, 5e-4, 1e-3])
    assert {"t9", "rho"} <= set(saved.zones[-1].optional_properties)
    out = capsys.readouterr().out
    assert f"saved: {folders[0].name}/output.xml (3 zones" in out
    assert "--zone-index 2" in out


def test_a_second_run_gets_its_own_folder(tmp_path, capsys):
    main(SHORT + ["--progress", "0"])
    main(SHORT + ["--progress", "0"])
    assert len([p for p in tmp_path.iterdir() if p.name.startswith("evolve-zone_")]) == 2


def test_no_output_saves_nothing(tmp_path, capsys):
    assert main(SHORT + ["--progress", "0", "--no-output"]) == 0
    assert list(tmp_path.iterdir()) == []
    out = capsys.readouterr().out
    assert out.startswith("evolve-zone: finished in") and "saved:" not in out


def test_output_names_the_folder(tmp_path, capsys):
    assert main(SHORT + ["--progress", "0", "--output", "runs/first"]) == 0
    assert (tmp_path / "runs" / "first" / "output.xml").is_file()


def test_saved_results_are_never_overwritten(tmp_path, monkeypatch, capsys):
    called = []
    monkeypatch.setattr(cli, "evolve_zone", lambda *a, **k: called.append(1))
    (tmp_path / "run").mkdir()
    (tmp_path / "run" / "output.xml").write_text("an earlier run")
    assert main(SHORT + ["--progress", "0", "--output", "run"]) == 1
    assert "already exists" in capsys.readouterr().err
    assert (tmp_path / "run" / "output.xml").read_text() == "an earlier run"
    assert not called, "the run started although its results could not be saved"


def test_the_summary_replaces_the_list_of_abundances(capsys):
    assert main(SHORT + ["--progress", "0", "--no-output"]) == 0
    lines = capsys.readouterr().out.splitlines()
    assert lines[0].startswith("evolve-zone: finished in") and lines[0].endswith("success")
    assert lines[1].strip().endswith("Jacobians")
    assert lines[2].strip().startswith("sum X = ") and "Ye = " in lines[2]
    assert lines[3].strip().startswith("largest mass fractions at t = 0.001 s")
    assert len(lines) < 10


def test_an_interrupted_run_leaves_no_empty_folder(tmp_path, monkeypatch, capsys):
    def interrupted(*args, **kwargs):
        raise KeyboardInterrupt

    monkeypatch.setattr(cli, "evolve_zone", interrupted)
    with pytest.raises(KeyboardInterrupt):
        main(SHORT + ["--progress", "0"])
    assert list(tmp_path.iterdir()) == []
