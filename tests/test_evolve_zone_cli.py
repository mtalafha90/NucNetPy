"""evolve-zone on networks too large to integrate quickly.

Given a complete JINA database, ``nucnetpy evolve-zone`` used to work for hours
without a word: it integrated every species in the file, offered no way to cut
the network or to set the tolerances, and printed nothing until it finished.
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
    assert "warning: the network has" in err and "--zmax" in err


def test_a_small_network_is_not_announced(capsys):
    main(SHORT + ["--progress", "0"])
    assert "warning" not in capsys.readouterr().err


def test_tolerances_reach_the_integrator(monkeypatch):
    seen = {}
    real = cli.evolve_zone

    def spy(*args, **kwargs):
        seen.update(kwargs)
        return real(*args, **kwargs)

    monkeypatch.setattr(cli, "evolve_zone", spy)
    main(SHORT + ["--rtol", "1e-8", "--atol", "1e-12", "--progress", "0"])
    assert seen["rtol"] == 1e-8 and seen["atol"] == 1e-12
