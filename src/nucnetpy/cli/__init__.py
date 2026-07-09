from __future__ import annotations
import argparse, sys
from pathlib import Path

from ..io.xml import read_xml, write_xml
from ..io.jina import read_jina_xml, combine_jina_xml, jina_database_summary
from ..analysis import largest_mass_fractions, species_history, flows, ydot, element_abundances, abundance_moment, energy_generation_rate, compare_rates, charge_changing_flows, system_timescales, entropy_generation_rate
from ..detailed_balance import net_flows as detailed_net_flows
from ..solver import evolve_zone, constant_thermo, time_grid
from ..graph import reaction_network_dot
from ..nse import solve_nse
from ..validation import validate_network, validate_zone


def _load(path):
    return read_xml(path)


def cmd_summary(args):
    n = _load(args.xml)
    print(f"species: {len(n.species)}")
    print(f"reactions: {len(n.reactions.reactions)}")
    print(f"zones: {len(n.zones)}")
    v = n.validate()
    print(f"missing_species: {len(v['missing_species'])}")
    print(f"invalid_reactions: {len(v['invalid_reactions'])}")


def cmd_print_output(args):
    n = _load(args.xml)
    cmd_summary(args)
    for i, z in enumerate(n.zones[:args.max_zones]):
        print(f"zone {i} label={z.label} props={z.properties}")
        for name, y in sorted(z.abundances.items()):
            if y >= args.min_abundance:
                print(f"  {name:8s} {y:.8e}")


def cmd_largest_x(args):
    n = _load(args.xml); z = n.zone(args.zone_index)
    for name, x in largest_mass_fractions(z, n.species, args.n, args.min_x):
        print(f"{name:8s} {x:.12e}")


def cmd_zone_abundances(args):
    n = _load(args.xml); z = n.zone(args.zone_index)
    for name, y in sorted(z.abundances.items()):
        if y >= args.min_abundance: print(f"{name} {y:.12e}")


def cmd_zone_properties(args):
    n = _load(args.xml); z = n.zone(args.zone_index)
    for k, v in sorted(z.properties.items()): print(k, v)


def cmd_element(args):
    n = _load(args.xml); z = n.zone(args.zone_index)
    for el, y in element_abundances(z, n.species).items():
        if el == args.element.lower() or args.element == '*': print(el, y)


def cmd_rates(args):
    n = _load(args.xml)
    for r, val in n.reactions.rates(args.t9, args.rho).items():
        if val >= args.min_rate: print(f"{r:60s} {val:.12e}")


def cmd_flows(args):
    n = _load(args.xml)
    for r, val in flows(n, args.zone_index, args.t9, args.rho).items():
        if val >= args.min_flow: print(f"{r:60s} {val:.12e}")


def cmd_ydot(args):
    n = _load(args.xml)
    for s, val in sorted(ydot(n, args.zone_index, args.t9, args.rho).items()):
        if abs(val) >= args.min_abs: print(f"{s:8s} {val:.12e}")


def cmd_net_flows(args):
    n = _load(args.xml); z = n.zone(args.zone_index)
    t9 = args.t9 or z.temperature9(); rho = args.rho or z.density()
    for r, (f, rev, net) in detailed_net_flows(n, z.abundances, t9=t9, rho=rho).items():
        if max(abs(f), abs(rev)) >= args.min_flow:
            print(f"{r:60s} fwd={f:.6e} rev={rev:.6e} net={net:.6e}")


def cmd_charge_flows(args):
    n = _load(args.xml)
    contributions = charge_changing_flows(n, args.zone_index, args.t9, args.rho)
    for r, val in sorted(contributions.items(), key=lambda kv: -abs(kv[1])):
        if abs(val) >= args.min_flow: print(f"{r:60s} {val:.12e}")
    print(f"dYe/dt {sum(contributions.values()):.12e}")


def cmd_timescales(args):
    n = _load(args.xml)
    ts = system_timescales(n, args.zone_index, args.t9, args.rho)
    for s, val in sorted(ts.items(), key=lambda kv: kv[1])[:args.n]:
        print(f"{s:8s} {val:.12e}")


def cmd_entropy(args):
    n = _load(args.xml)
    print(entropy_generation_rate(n, args.zone_index, args.t9, args.rho))


def cmd_conservation(args):
    n = _load(args.xml)
    bad = n.reactions.invalid_reactions(n.species)
    print(f"invalid_reactions {len(bad)}")
    for r in bad[:args.max]:
        ok, da, dz = r.conserves_a_z(n.species)
        print(f"{r.string} dA={da} dZ={dz}")


def cmd_remove_duplicates(args):
    n = _load(args.xml); before = len(n.reactions.reactions)
    n.reactions.remove_duplicates(); write_xml(n, args.output)
    print(f"removed {before - len(n.reactions.reactions)} duplicate reactions")


def cmd_remove_invalid(args):
    n = _load(args.xml); before = len(n.reactions.reactions)
    n.reactions.filter_valid(n.species); write_xml(n, args.output)
    print(f"removed {before - len(n.reactions.reactions)} invalid reactions")


def cmd_export_zone_xml(args):
    n = _load(args.xml); out = type(n)()
    out.species = n.species; out.reactions = n.reactions
    out.zones = [n.zone(args.zone_index)]
    write_xml(out, args.output)


def cmd_reactions_latex(args):
    n = _load(args.xml)
    with open(args.output, 'w') as f:
        f.write('\\begin{tabular}{ll}\nReaction & $Q$ [MeV]\\\\\n')
        for r in n.reactions.reactions:
            f.write(f"${r.latex}$ & {r.q_value:.6g} \\\\\n")
        f.write('\\end{tabular}\n')


def cmd_species_history(args):
    n = _load(args.xml)
    for i, label, y in species_history(n, args.species): print(i, '|'.join(label), f"{y:.12e}")


def cmd_evolve(args):
    n = _load(args.xml); z = n.zone(args.zone_index)
    times = time_grid(args.t0, args.t1, args.steps, args.log_time)
    res = evolve_zone(n, z, times, thermo=constant_thermo(args.t9 or z.temperature9(), args.rho or z.density()), method=args.method)
    for s, y in sorted(res.final_abundances.items()):
        if y >= args.min_abundance: print(s, f"{y:.12e}")


def cmd_energy(args):
    n = _load(args.xml)
    print(energy_generation_rate(n, args.zone_index, args.t9, args.rho))


def cmd_dot(args):
    n = _load(args.xml)
    dot = reaction_network_dot(n, min_rate=args.min_rate, t9=args.t9, rho=args.rho)
    if args.output: Path(args.output).write_text(dot)
    else: print(dot)



def cmd_validate(args):
    n = _load(args.xml)
    issues = validate_network(n, strict=args.strict)
    for i, z in enumerate(n.zones[:args.max_zones]):
        issues.extend(validate_zone(z, n))
    print(f"issues {len(issues)}")
    for issue in issues[:args.max]:
        print(f"{issue.level} {issue.code}: {issue.message}")


def cmd_nse(args):
    n = _load(args.xml)
    corr = None
    if args.coulomb:
        from ..coulomb import nse_correction as corr
    res = solve_nse(n, t9=args.t9, rho=args.rho, ye=args.ye, nse_correction=corr)
    print(f"success {res.success} mu_p {res.mu_p:.8e} mu_n {res.mu_n:.8e} xsum {res.xsum:.8e} ye {res.computed_ye:.8e}")
    for name, y in sorted(res.abundances.items()):
        x = n.species.get(name).a * y if name in n.species else y
        if y >= args.min_abundance or x >= args.min_x:
            print(f"{name:8s} Y={y:.12e} X={x:.12e}")


def cmd_qse(args):
    from ..qse import solve_qse, QSECluster
    n = _load(args.xml)
    clusters = []
    for spec in args.cluster or []:
        names, _, constraint = spec.rpartition(':')
        clusters.append(QSECluster([s for s in names.split(',') if s], float(constraint)))
    corr = None
    if args.coulomb:
        from ..coulomb import nse_correction as corr
    res = solve_qse(n, t9=args.t9, rho=args.rho, ye=args.ye, clusters=clusters, nse_correction=corr)
    lam = ' '.join(f"{v:.6e}" for v in res.lambdas)
    print(f"success {res.success} mu_p {res.mu_p:.8e} mu_n {res.mu_n:.8e} lambdas [{lam}] xsum {res.xsum:.8e} ye {res.computed_ye:.8e}")
    for name, y in sorted(res.abundances.items()):
        x = n.species.get(name).a * y if name in n.species else y
        if y >= args.min_abundance or x >= args.min_x:
            print(f"{name:8s} Y={y:.12e} X={x:.12e}")


def cmd_jina_summary(args):
    n = read_jina_xml(args.nuclides_xml, args.reactions_xml, zones_xml=args.zones_xml)
    summary = jina_database_summary(n)
    for key, value in summary.items():
        print(f"{key}: {value}")
    if args.show_invalid:
        for r in n.reactions.invalid_reactions(n.species)[:args.max_invalid]:
            ok, da, dz = r.conserves_a_z(n.species)
            print(f"INVALID {r.string} dA={da} dZ={dz}")


def cmd_jina_combine(args):
    n = combine_jina_xml(args.nuclides_xml, args.reactions_xml, args.output_xml, zones_xml=args.zones_xml)
    print(f"wrote {args.output_xml}")
    print(f"species: {len(n.species)}")
    print(f"reactions: {len(n.reactions.reactions)}")
    print(f"zones: {len(n.zones)}")

def build_parser():
    p = argparse.ArgumentParser(prog='nucnetpy', description='Pure-Python NucNet Tools replacement commands')
    sub = p.add_subparsers(dest='cmd', required=True)
    def add_xml(sp): sp.add_argument('xml')
    s = sub.add_parser('summary'); add_xml(s); s.set_defaults(func=cmd_summary)
    s = sub.add_parser('print-output'); add_xml(s); s.add_argument('--max-zones', type=int, default=3); s.add_argument('--min-abundance', type=float, default=0.0); s.set_defaults(func=cmd_print_output)
    s = sub.add_parser('largest-x'); add_xml(s); s.add_argument('-n', type=int, default=10); s.add_argument('--zone-index', type=int, default=0); s.add_argument('--min-x', type=float, default=0.0); s.set_defaults(func=cmd_largest_x)
    s = sub.add_parser('zone-abundances'); add_xml(s); s.add_argument('--zone-index', type=int, default=0); s.add_argument('--min-abundance', type=float, default=0.0); s.set_defaults(func=cmd_zone_abundances)
    s = sub.add_parser('zone-properties'); add_xml(s); s.add_argument('--zone-index', type=int, default=0); s.set_defaults(func=cmd_zone_properties)
    s = sub.add_parser('element-abundances'); add_xml(s); s.add_argument('element'); s.add_argument('--zone-index', type=int, default=0); s.set_defaults(func=cmd_element)
    s = sub.add_parser('rates'); add_xml(s); s.add_argument('--t9', type=float, default=1.0); s.add_argument('--rho', type=float, default=1.0); s.add_argument('--min-rate', type=float, default=0.0); s.set_defaults(func=cmd_rates)
    s = sub.add_parser('flows'); add_xml(s); s.add_argument('--zone-index', type=int, default=0); s.add_argument('--t9', type=float); s.add_argument('--rho', type=float); s.add_argument('--min-flow', type=float, default=0.0); s.set_defaults(func=cmd_flows)
    s = sub.add_parser('ydot'); add_xml(s); s.add_argument('--zone-index', type=int, default=0); s.add_argument('--t9', type=float); s.add_argument('--rho', type=float); s.add_argument('--min-abs', type=float, default=0.0); s.set_defaults(func=cmd_ydot)
    s = sub.add_parser('net-flows', help='forward, detailed-balance reverse, and net fluxes'); add_xml(s); s.add_argument('--zone-index', type=int, default=0); s.add_argument('--t9', type=float); s.add_argument('--rho', type=float); s.add_argument('--min-flow', type=float, default=0.0); s.set_defaults(func=cmd_net_flows)
    s = sub.add_parser('charge-flows', help='per-reaction dYe/dt contributions'); add_xml(s); s.add_argument('--zone-index', type=int, default=0); s.add_argument('--t9', type=float); s.add_argument('--rho', type=float); s.add_argument('--min-flow', type=float, default=0.0); s.set_defaults(func=cmd_charge_flows)
    s = sub.add_parser('timescales', help='shortest species timescales Y/|dY/dt|'); add_xml(s); s.add_argument('--zone-index', type=int, default=0); s.add_argument('--t9', type=float); s.add_argument('--rho', type=float); s.add_argument('-n', type=int, default=20); s.set_defaults(func=cmd_timescales)
    s = sub.add_parser('entropy-generation', help='dS/dt per nucleon in k_B/s'); add_xml(s); s.add_argument('--zone-index', type=int, default=0); s.add_argument('--t9', type=float); s.add_argument('--rho', type=float); s.set_defaults(func=cmd_entropy)
    s = sub.add_parser('conservation'); add_xml(s); s.add_argument('--max', type=int, default=20); s.set_defaults(func=cmd_conservation)
    s = sub.add_parser('remove-duplicates'); add_xml(s); s.add_argument('output'); s.set_defaults(func=cmd_remove_duplicates)
    s = sub.add_parser('remove-invalid'); add_xml(s); s.add_argument('output'); s.set_defaults(func=cmd_remove_invalid)
    s = sub.add_parser('export-zone-xml'); add_xml(s); s.add_argument('output'); s.add_argument('--zone-index', type=int, default=0); s.set_defaults(func=cmd_export_zone_xml)
    s = sub.add_parser('reactions-latex'); add_xml(s); s.add_argument('output'); s.set_defaults(func=cmd_reactions_latex)
    s = sub.add_parser('species-history'); add_xml(s); s.add_argument('species'); s.set_defaults(func=cmd_species_history)
    s = sub.add_parser('evolve-zone'); add_xml(s); s.add_argument('--zone-index', type=int, default=0); s.add_argument('--t0', type=float, default=0.0); s.add_argument('--t1', type=float, default=1.0); s.add_argument('--steps', type=int, default=50); s.add_argument('--t9', type=float); s.add_argument('--rho', type=float); s.add_argument('--method', default='bdf'); s.add_argument('--log-time', action='store_true'); s.add_argument('--min-abundance', type=float, default=0.0); s.set_defaults(func=cmd_evolve)
    s = sub.add_parser('energy-generation'); add_xml(s); s.add_argument('--zone-index', type=int, default=0); s.add_argument('--t9', type=float); s.add_argument('--rho', type=float); s.set_defaults(func=cmd_energy)
    s = sub.add_parser('net-dot'); add_xml(s); s.add_argument('-o','--output'); s.add_argument('--t9', type=float, default=1.0); s.add_argument('--rho', type=float, default=1.0); s.add_argument('--min-rate', type=float, default=0.0); s.set_defaults(func=cmd_dot)
    s = sub.add_parser('validate'); add_xml(s); s.add_argument('--strict', action='store_true'); s.add_argument('--max', type=int, default=50); s.add_argument('--max-zones', type=int, default=10); s.set_defaults(func=cmd_validate)
    s = sub.add_parser('nse'); add_xml(s); s.add_argument('--t9', type=float, required=True); s.add_argument('--rho', type=float, required=True); s.add_argument('--ye', type=float, required=True); s.add_argument('--coulomb', action='store_true', help='apply Bravo & Garcia-Senz Coulomb corrections'); s.add_argument('--min-abundance', type=float, default=0.0); s.add_argument('--min-x', type=float, default=0.0); s.set_defaults(func=cmd_nse)
    s = sub.add_parser('qse', help='constrained cluster equilibrium (libnuceq)'); add_xml(s); s.add_argument('--t9', type=float, required=True); s.add_argument('--rho', type=float, required=True); s.add_argument('--ye', type=float, required=True); s.add_argument('--cluster', action='append', metavar='SP1,SP2,...:Y', help='cluster species and constrained total abundance; repeatable'); s.add_argument('--coulomb', action='store_true'); s.add_argument('--min-abundance', type=float, default=0.0); s.add_argument('--min-x', type=float, default=0.0); s.set_defaults(func=cmd_qse)

    s = sub.add_parser('jina-summary', help='summarize separate JINA nuclide and reaction XML files')
    s.add_argument('nuclides_xml')
    s.add_argument('reactions_xml')
    s.add_argument('--zones-xml')
    s.add_argument('--show-invalid', action='store_true')
    s.add_argument('--max-invalid', type=int, default=20)
    s.set_defaults(func=cmd_jina_summary)

    s = sub.add_parser('jina-combine', help='combine JINA nuclide and reaction XML files into one nucnetpy XML file')
    s.add_argument('nuclides_xml')
    s.add_argument('reactions_xml')
    s.add_argument('output_xml')
    s.add_argument('--zones-xml')
    s.set_defaults(func=cmd_jina_combine)
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    return args.func(args)
