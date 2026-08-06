"""The archived silicon-burning reference case must keep reproducing.

validation/reference/ holds the extracted network, its initial zone and the
expected results for the principal demonstration of the accompanying article.
These tests are what stop that archive drifting away from the code: if a change
alters the reference result, it fails here rather than in a reader's hands.
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from nucnetpy.io.xml import read_xml
from nucnetpy.nse import solve_nse
from nucnetpy.solver import constant_thermo, evolve_zone, time_grid

REFERENCE = Path(__file__).resolve().parents[1] / "validation" / "reference"
pytestmark = pytest.mark.skipif(not REFERENCE.is_dir(),
                                reason="reference archive not present")


@pytest.fixture(scope="module")
def spec():
    return json.loads((REFERENCE / "si_burning_expected.json").read_text())


@pytest.fixture(scope="module")
def network(spec):
    return read_xml(str(REFERENCE / spec["network"]["file"]))


def test_archived_network_matches_its_checksum(spec):
    import hashlib
    path = REFERENCE / spec["network"]["file"]
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    assert digest == spec["network"]["sha256"], (
        "the archived network has changed; the expectations beside it were "
        "produced from a different file")


def test_archived_network_has_the_documented_size(network, spec):
    assert len(network.species) == spec["network"]["species"]
    assert len(network.reactions.reactions) == spec["network"]["reactions"]
    assert network.zones, "the archive must carry its initial zone"


def test_reference_case_reproduces(network, spec):
    c, exp, tol = spec["conditions"], spec["expected"], spec["tolerances"]
    zone = network.zone(0)
    ye = zone.ye(network.species)
    assert np.isclose(ye, c["initial_ye"], rtol=tol["rtol"], atol=tol["atol"])

    result = evolve_zone(network, zone, time_grid(0.0, c["t_end"], c["steps"]),
                         thermo=constant_thermo(c["t9"], c["rho"]),
                         method=c["method"], rtol=c["rtol"], atol=c["atol"])
    assert result.success, result.message

    final = result.final_abundances
    for name, want in exp["final_abundances"].items():
        got = float(final.get(name, 0.0))
        assert np.isclose(got, want, rtol=tol["rtol"], atol=tol["atol"]), name

    xsum = sum(network.species[k].a * v for k, v in final.items()
               if k in network.species)
    assert np.isclose(xsum, exp["final_xsum"], rtol=0.0, atol=tol["xsum_atol"])


def test_reference_nse_reproduces(network, spec):
    c, exp, tol = spec["conditions"], spec["expected"], spec["tolerances"]
    sp = [n for n in network.species if n not in {"n", "h1"}]
    nse = solve_nse(network, t9=c["t9"], rho=c["rho"], ye=c["initial_ye"],
                    species=sp)
    assert nse.success
    assert np.isclose(nse.mu_p, exp["nse_mu_p"], rtol=tol["rtol"], atol=tol["atol"])
    assert np.isclose(nse.mu_n, exp["nse_mu_n"], rtol=tol["rtol"], atol=tol["atol"])
    for name, want in exp["nse_abundances"].items():
        got = float(nse.abundances.get(name, 0.0))
        assert np.isclose(got, want, rtol=tol["rtol"], atol=tol["atol"]), name


def test_reproduction_script_runs_and_passes():
    """The script a reader is told to run must actually exit 0."""
    script = REFERENCE.parent / "reproduce_si_burning.py"
    proc = subprocess.run([sys.executable, str(script), "--quiet"],
                          capture_output=True, text=True, timeout=600)
    assert proc.returncode == 0, proc.stdout + proc.stderr


def test_provenance_records_the_source_database():
    """The database is not redistributed, so it must at least be identified."""
    prov = json.loads((REFERENCE / "si_burning_provenance.json").read_text())
    assert prov["files"], "no source files recorded"
    for entry in prov["files"]:
        assert len(entry["sha256"]) == 64
        assert entry["bytes"] > 0
        assert entry["source_labels"], "no data-release labels recorded"
    parsed = prov["parsed"]
    # The distinction that makes the species counts intelligible.
    assert (parsed["nuclide_file_entries"] + parsed["placeholders_without_nuclear_data"]
            == parsed["species_in_memory"])


def test_detailed_balance_result_reproduces(network, spec):
    """The article's headline figure must keep coming out of the archive."""
    from nucnetpy.detailed_balance import consistent_reverse_network
    exp = spec["expected"].get("detailed_balance")
    if not exp:
        pytest.skip("archive predates the detailed-balance expectations")
    c, tol = spec["conditions"], spec["tolerances"]
    rtol = tol.get("detailed_balance_rtol", 1.0e-3)
    for label, want in exp.items():
        net = consistent_reverse_network(
            read_xml(str(REFERENCE / spec["network"]["file"])),
            tabulate=(label == "tabulated"))
        res = evolve_zone(net, net.zone(0),
                          time_grid(0.0, c["t_end"], c["steps"]),
                          thermo=constant_thermo(c["t9"], c["rho"]),
                          method=c["method"], rtol=c["rtol"], atol=c["atol"])
        assert res.success, res.message
        sp = [n for n in net.species if n not in {"n", "h1"}]
        nse = solve_nse(net, t9=c["t9"], rho=c["rho"], ye=c["initial_ye"],
                        species=sp)
        xn = {k: net.species[k].a * v for k, v in res.final_abundances.items()
              if k in net.species}
        xq = {k: net.species[k].a * v for k, v in nse.abundances.items()
              if k in net.species}
        d = [abs(xn.get(k, 0.0) - v) / v for k, v in xq.items()
             if k not in {"n", "h1"} and v >= 1.0e-6]
        assert np.isclose(float(np.median(d)), want["median_rel_diff_vs_nse"],
                          rtol=rtol), label
        assert np.isclose(float(max(d)), want["max_rel_diff_vs_nse"],
                          rtol=rtol), label
