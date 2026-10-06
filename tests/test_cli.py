"""Command-line behaviour that a user relies on but the library tests do not see."""
from pathlib import Path

import pytest

from nucnetpy import read_xml, evolve_zone, time_grid
from nucnetpy.cli import main

GOLDEN = str(Path(__file__).parent / "golden" / "golden_network.xml")


def test_zone_properties_include_optional_properties(capsys):
    # libnucnet zone files keep T9 and rho among the optional properties.
    assert main(["zone-properties", GOLDEN]) in (None, 0)
    out = capsys.readouterr().out.split()
    assert out == ["rho", "1e5", "t9", "0.5"]


def test_energy_generation_is_reported_in_erg_per_gram_per_second(capsys):
    main(["energy-generation", GOLDEN, "--t9", "2", "--rho", "1e5"])
    lines = capsys.readouterr().out.splitlines()
    assert len(lines) == 2 and all(line.endswith("erg/g/s") for line in lines)
    masses, q = (float(line.split()[-2]) for line in lines)
    assert masses == pytest.approx(3.7657e22, rel=1e-4)
    assert q == pytest.approx(masses, rel=1e-6)


@pytest.mark.parametrize("argv, message", [
    (["summary", "no_such_file.xml"], "No such file"),
    (["zone-abundances", GOLDEN, "--zone-index", "5"], "zone index 5 is out of range"),
])
def test_input_errors_give_one_line_and_exit_status_one(argv, message, capsys):
    assert main(argv) == 1
    err = capsys.readouterr().err
    assert message in err and "Traceback" not in err


def test_unknown_method_is_rejected_by_the_command_line(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["evolve-zone", GOLDEN, "--method", "bfd"])
    assert exc.value.code == 2
    assert "invalid choice" in capsys.readouterr().err


def test_unknown_method_is_rejected_by_the_library():
    net = read_xml(GOLDEN)
    with pytest.raises(ValueError, match="unknown method"):
        evolve_zone(net, net.zone(0), time_grid(0.0, 1e-6, 3), method="nonsense")
    with pytest.raises(ValueError, match="unknown jac_mode"):
        evolve_zone(net, net.zone(0), time_grid(0.0, 1e-6, 3), jac_mode="numeric")


def test_method_names_are_case_insensitive():
    net = read_xml(GOLDEN)
    a = evolve_zone(net, net.zone(0), time_grid(0.0, 1e-6, 3), method="RK4")
    b = evolve_zone(net, net.zone(0), time_grid(0.0, 1e-6, 3), method="rk4")
    assert (a.y == b.y).all()
