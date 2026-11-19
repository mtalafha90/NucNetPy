"""evolve-zone on large networks and in explosive burning.

Given a complete JINA database, ``nucnetpy evolve-zone`` must let the user cut
the network and set the tolerances and the step control, report how far a long
integration has got, and say what to try when it advances too slowly to finish.
"""
from pathlib import Path

import numpy as np
import pytest

import nucnetpy.cli as cli
from nucnetpy import evolve_zone, read_xml, time_grid
from nucnetpy.cli import main

GOLDEN = str(Path(__file__).parent / "golden" / "golden_network.xml")
SHORT = ["evolve-zone", GOLDEN, "--t1", "1e-3", "--steps", "3"]


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
    assert main(SHORT + ["--zmax", "6", "--progress", "0"]) == 0
    lines = capsys.readouterr().out.splitlines()
    assert lines[0].startswith("success True")
    names = {line.split()[0] for line in lines[1:]}
    assert names and names <= {"n", "h1", "he4", "be8", "c12"}


def test_progress_lines_go_to_standard_error(capsys):
    main(SHORT + ["--progress", "1e-9"])
    captured = capsys.readouterr()
    assert "evolve-zone: reached t =" in captured.err
    assert "reached" not in captured.out


def test_progress_can_be_turned_off(capsys):
    main(SHORT + ["--progress", "0"])
    assert "reached" not in capsys.readouterr().err


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
