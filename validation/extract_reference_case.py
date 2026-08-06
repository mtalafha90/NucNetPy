"""Extract the silicon-burning reference case into a self-contained archive.

The demonstration of the accompanying article is driven by a third-party
JINA/libnucnet database that cannot be redistributed. Everything the principal
result depends on, however, is a 15-nuclide, 32-reaction subnetwork of it. This
script writes that subnetwork out together with its initial zone, the exact
conditions, the expected results and a provenance record of the database it
came from, so that the result can be reproduced without the database.

Run it once, against the same database used for the article:

    python validation/extract_reference_case.py \\
        --nuclides nuclides.xml --reactions reaction_data.xml

It writes into validation/reference/:

    si_burning_network.xml    the extracted network and its initial zone
    si_burning_expected.json  expected outputs, with the tolerances to use
    si_burning_provenance.json  checksums and provenance of the source files

Reproducing the case afterwards needs only the first of those; see
validation/reproduce_si_burning.py.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import platform
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Dict

import numpy as np

from nucnetpy.core import Zone
from nucnetpy.io.jina import read_jina_xml
from nucnetpy.io.xml import write_xml
from nucnetpy.network_limiter import limit_network
from nucnetpy.nse import solve_nse
from nucnetpy.solver import constant_thermo, evolve_zone, time_grid

# Must match validation/demonstration.py exactly, or the archived case would
# not be the case the article reports.
ALPHA_CHAIN = ["gamma", "n", "h1", "he4", "c12", "o16", "ne20", "mg24", "si28",
               "s32", "ar36", "ca40", "ti44", "cr48", "fe52", "ni56"]

T9 = 5.0
RHO = 1.0e8
T_END = 10.0
STEPS = 120
RTOL = 1.0e-8
ATOL = 1.0e-14

REFERENCE_DIR = Path(__file__).resolve().parent / "reference"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def source_labels(path: Path, limit: int = 12) -> Dict[str, int]:
    """Histogram of <source> labels, which identifies the data release."""
    import re
    text = Path(path).read_text(errors="replace")
    counts = Counter(re.findall(r"<source>([^<]*)</source>", text))
    return dict(counts.most_common(limit))


def environment_record() -> Dict[str, str]:
    import scipy
    import nucnetpy
    return {
        "nucnetpy": nucnetpy.__version__,
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
        "processor": platform.processor() or platform.machine(),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nuclides", required=True)
    ap.add_argument("--reactions", required=True)
    ap.add_argument("--outdir", default=str(REFERENCE_DIR))
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    nuc, rea = Path(args.nuclides), Path(args.reactions)
    print(f"reading {nuc} and {rea} ...", flush=True)
    t0 = time.perf_counter()
    base = read_jina_xml(str(nuc), str(rea))
    load_seconds = time.perf_counter() - t0

    full_species = len(base.species)
    full_reactions = len(base.reactions.reactions)
    placeholders = base.species_without_nuclear_data()

    # ---- the extracted case -------------------------------------------------
    net = copy.deepcopy(base)
    limit_network(net, ALPHA_CHAIN)
    zone = Zone(abundances={"si28": 1.0 / 28.0})
    net.zones = [zone]
    ye = zone.ye(net.species)
    print(f"extracted {len(net.species)} species, "
          f"{len(net.reactions.reactions)} reactions", flush=True)

    network_path = outdir / "si_burning_network.xml"
    write_xml(net, network_path)

    # Read it back and run from the archived file, not from memory: the archive
    # is only useful if it is what reproduces the numbers.
    from nucnetpy.io.xml import read_xml
    check = read_xml(str(network_path))
    assert len(check.species) == len(net.species), "round trip lost species"
    assert len(check.reactions.reactions) == len(net.reactions.reactions), \
        "round trip lost reactions"
    assert check.zones, "round trip lost the initial zone"

    zone_r = check.zone(0)
    times = time_grid(0.0, T_END, STEPS)
    result = evolve_zone(check, zone_r, times, thermo=constant_thermo(T9, RHO),
                         method="bdf", rtol=RTOL, atol=ATOL)
    if not result.success:
        print(f"evolution failed: {result.message}", file=sys.stderr)
        return 1
    nse = solve_nse(check, t9=T9, rho=RHO, ye=ye)

    final = result.final_abundances
    xsum = sum(check.species[k].a * v for k, v in final.items() if k in check.species)

    # The detailed-balance variant is the article's headline result, so the
    # archive carries it too.  Reverse rates are rebuilt from the equilibrium
    # constant, which shares its prefactor with solve_nse; the agreement below
    # therefore measures mutual consistency, not the correctness of that
    # shared formulation.
    from nucnetpy.detailed_balance import consistent_reverse_network
    db = {}
    for label, tabulate in (("function", False), ("tabulated", True)):
        dnet = consistent_reverse_network(read_xml(str(network_path)),
                                          tabulate=tabulate)
        dzone = dnet.zone(0)
        dres = evolve_zone(dnet, dzone, times, thermo=constant_thermo(T9, RHO),
                           method="bdf", rtol=RTOL, atol=ATOL)
        dnse = solve_nse(dnet, t9=T9, rho=RHO, ye=ye)
        xn = {k: dnet.species[k].a * v for k, v in dres.final_abundances.items()
              if k in dnet.species}
        xq = {k: dnet.species[k].a * v for k, v in dnse.abundances.items()
              if k in dnet.species}
        diffs = [abs(xn.get(k, 0.0) - v) / v for k, v in xq.items()
                 if k not in {"n", "h1"} and v >= 1.0e-6]
        db[label] = {"evolve_success": bool(dres.success),
                     "median_rel_diff_vs_nse": float(np.median(diffs)),
                     "max_rel_diff_vs_nse": float(max(diffs))}

    expected = {
        "case": "silicon burning",
        "description": "Pure Si-28 burned to its stationary composition and "
                       "compared with a separately solved NSE.",
        "conditions": {"t9": T9, "rho": RHO, "t_end": T_END, "steps": STEPS,
                       "method": "bdf", "rtol": RTOL, "atol": ATOL,
                       "initial_abundances": {"si28": 1.0 / 28.0},
                       "initial_ye": ye},
        "network": {"file": network_path.name,
                    "sha256": sha256(network_path),
                    "species": len(check.species),
                    "reactions": len(check.reactions.reactions)},
        "expected": {
            "evolve_success": True,
            "final_abundances": {k: float(v) for k, v in sorted(final.items())},
            "final_xsum": float(xsum),
            "nse_success": bool(nse.success),
            "nse_abundances": {k: float(v) for k, v in sorted(nse.abundances.items())},
            "nse_mu_p": float(nse.mu_p),
            "nse_mu_n": float(nse.mu_n),
            "nse_xsum": float(nse.xsum),
            "nse_ye": float(nse.computed_ye),
            "detailed_balance": db,
        },
        "tolerances": {
            "comment": "Abundances are compared with numpy.isclose using these "
                       "values. They are loose enough to absorb BLAS and "
                       "library differences between platforms and tight enough "
                       "that a change in the physics fails the check.",
            "rtol": 1.0e-6,
            "atol": 1.0e-18,
            "xsum_atol": 1.0e-9,
            "detailed_balance_rtol": 1.0e-3,
        },
        "environment": environment_record(),
    }
    (outdir / "si_burning_expected.json").write_text(json.dumps(expected, indent=2) + "\n")

    provenance = {
        "comment": "Provenance of the third-party database the reference case "
                   "was extracted from. The database itself is not "
                   "redistributed; these records identify it.",
        "schema": "libnucnet 2019-01-15",
        "files": [
            {"role": "nuclides", "name": nuc.name, "bytes": nuc.stat().st_size,
             "sha256": sha256(nuc), "source_labels": source_labels(nuc)},
            {"role": "reactions", "name": rea.name, "bytes": rea.stat().st_size,
             "sha256": sha256(rea), "source_labels": source_labels(rea)},
        ],
        "parsed": {
            "nuclide_file_entries": full_species - len(placeholders),
            "species_in_memory": full_species,
            "placeholders_without_nuclear_data": len(placeholders),
            "reactions": full_reactions,
            "load_seconds": load_seconds,
        },
        "extraction": {"species_kept": ALPHA_CHAIN,
                       "note": "gamma is retained so photodisintegration "
                               "records survive the cut"},
        "environment": environment_record(),
    }
    (outdir / "si_burning_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")

    print(f"wrote {network_path}")
    print(f"wrote {outdir/'si_burning_expected.json'}")
    print(f"wrote {outdir/'si_burning_provenance.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
