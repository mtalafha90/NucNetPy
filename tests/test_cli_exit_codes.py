"""A command that failed must say so and exit non-zero.

The package's own advice is that a solver result is meaningless until its
success flag has been checked. The command line has to follow that advice
itself: a failed solve still carries a full abundance vector, so printing the
numbers without the verdict makes a diverged run look like an ordinary one.
"""
import subprocess
import sys

import pytest

from nucnetpy import Network, Reaction, Species, Zone
from nucnetpy.io.xml import write_xml


def run(*argv, cwd=None):
    return subprocess.run([sys.executable, "-m", "nucnetpy.cli", *argv],
                          capture_output=True, text=True, timeout=300, cwd=cwd)


@pytest.fixture
def failing_network(tmp_path):
    """A rate large enough that the stiff solver cannot take its first step."""
    net = Network()
    net.add_species(Species("he4", 2, 4, mass_excess=2.425))
    net.add_species(Species("c12", 6, 12, mass_excess=0.0))
    net.reactions.add(Reaction.from_names(["he4", "he4", "he4"], ["c12"],
                                          constant_rate=1.0e300))
    net.add_zone(Zone(abundances={"he4": 0.25}))
    path = tmp_path / "fails.xml"
    write_xml(net, path)
    return path


@pytest.fixture
def working_network(tmp_path):
    net = Network()
    net.add_species(Species("he4", 2, 4, mass_excess=2.425))
    net.add_species(Species("c12", 6, 12, mass_excess=0.0))
    net.reactions.add(Reaction.from_names(["he4", "he4", "he4"], ["c12"],
                                          constant_rate=1.0e-6))
    net.add_zone(Zone(abundances={"he4": 0.25}))
    path = tmp_path / "works.xml"
    write_xml(net, path)
    return path


def test_failed_evolution_exits_non_zero(failing_network):
    proc = run("evolve-zone", str(failing_network), "--t9", "1.0",
               "--rho", "1e8", "--t0", "0", "--t1", "1e6", "--steps", "5")
    assert proc.returncode != 0, (
        "a failed evolution exited 0; a diverged run is indistinguishable "
        "from a successful one at the command line")


def test_failed_evolution_reports_the_failure(failing_network):
    proc = run("evolve-zone", str(failing_network), "--t9", "1.0",
               "--rho", "1e8", "--t0", "0", "--t1", "1e6", "--steps", "5")
    assert "success False" in proc.stdout
    assert proc.stderr.strip(), "no diagnostic was written to stderr"


def test_successful_evolution_exits_zero_and_prints_abundances(working_network):
    proc = run("evolve-zone", str(working_network), "--t9", "1.0",
               "--rho", "1e6", "--t0", "0", "--t1", "1.0", "--steps", "20")
    assert proc.returncode == 0, proc.stderr
    assert "success True" in proc.stdout
    assert "he4" in proc.stdout


def test_nse_reports_status_and_exit_code(working_network):
    proc = run("nse", str(working_network), "--t9", "5", "--rho", "1e8",
               "--ye", "0.5")
    assert proc.stdout.startswith("success "), proc.stdout[:80]
    expected = 0 if proc.stdout.startswith("success True") else 1
    assert proc.returncode == expected
