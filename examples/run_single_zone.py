"""Evolve one zone of a JINA/libnucnet database read from separate XML files.

Unlike ``use_nucnetpy.py``, which builds a toy network in code, this script
works from real data files, which is the usual starting point.

Run with:

    python examples/run_single_zone.py \\
        --nuclides nuclides.xml \\
        --reactions reaction_data.xml \\
        --zone zone.xml

The thermodynamic conditions and the end time are taken from the zone's own
properties where it supplies them, so the same command works unchanged on a
different database.
"""
import argparse
import sys

from nucnetpy import read_jina_xml
from nucnetpy.solver import constant_thermo, evolve_zone, time_grid


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nuclides", required=True, help="JINA nuclide XML file")
    ap.add_argument("--reactions", required=True, help="JINA reaction XML file")
    ap.add_argument("--zone", help="zone XML file supplying the initial composition")
    ap.add_argument("--steps", type=int, default=200)
    ap.add_argument("--rtol", type=float, default=1.0e-8)
    ap.add_argument("--atol", type=float, default=1.0e-20)
    ap.add_argument("--min-x", type=float, default=1.0e-20,
                    help="do not print species below this mass fraction")
    args = ap.parse_args()

    net = read_jina_xml(nuclides_xml=args.nuclides,
                        reactions_xml=args.reactions,
                        zones_xml=args.zone)

    print(f"species:   {len(net.species)}")
    print(f"reactions: {len(net.reactions.reactions)}")
    print(f"zones:     {len(net.zones)}")

    if not net.zones:
        print("No zone in the input, so there is no composition to evolve.",
              file=sys.stderr)
        return 1

    zone = net.zone(0)

    # A libnucnet zone file carries its metadata inside <optional_properties>,
    # which the reader keeps in zone.optional_properties.  zone.properties holds
    # only properties written directly on the <zone> element, so a script that
    # consults one and not the other silently gets nothing and falls back to
    # its defaults.
    props = {**zone.optional_properties, **zone.properties}
    print("\nInitial zone properties:")
    for k, v in sorted(props.items()):
        print(f"  {k} = {v}")
    if not props:
        print("  (none recorded in the file)")

    def prop(*names, default):
        for n in names:
            if n in props:
                try:
                    return float(props[n])
                except ValueError:
                    pass
        return default

    t9 = prop("t9_0", "t9", default=0.20)
    rho = prop("rho_0", "rho", default=1.5e4)
    t_end = prop("tend", default=100.0)
    print(f"\nEvolving at T9 = {t9}, rho = {rho} g/cm^3 to t = {t_end} s")

    result = evolve_zone(net, zone, time_grid(0.0, t_end, args.steps),
                         thermo=constant_thermo(t9=t9, rho=rho),
                         method="bdf", rtol=args.rtol, atol=args.atol)

    # A failed solve, or a positivity projection that had to invent baryon
    # number, both surface here.  A composition from a failed run looks
    # perfectly plausible, so this check is not optional.
    if not result.success:
        print(f"\nEvolution failed: {result.message}", file=sys.stderr)
        return 1

    final = result.final_abundances
    xsum = sum(net.species[k].a * v for k, v in final.items() if k in net.species)
    print(f"\nSolver: {result.message}")
    print(f"Right-hand-side evaluations: {result.nfev}, Jacobians: {result.njev}")
    print(f"Sum of mass fractions: {xsum:.12f}  (1 - sum = {1.0 - xsum:.3e})")

    print(f"\nFinal composition above X = {args.min_x:g}:")
    rows = [(k, v, net.species[k].a * v) for k, v in final.items()
            if k in net.species and net.species[k].a * v > args.min_x]
    for name, y, x in sorted(rows, key=lambda r: -r[2]):
        print(f"  {name:8s} Y = {y:.8e}   X = {x:.8e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
