from __future__ import annotations
import argparse
import sys
from pathlib import Path

from ..io.xml import read_xml, write_xml, write_evolution_xml
from ..io.jina import read_jina_xml, combine_jina_xml, jina_database_summary
from ..analysis import (
    charge_changing_flows,
    element_abundances,
    energy_generation_rate,
    entropy_generation_rate,
    flows,
    largest_mass_fractions,
    species_history,
    system_timescales,
    ydot,
)
from ..detailed_balance import net_flows as detailed_net_flows
from ..solver import evolve_zone, constant_thermo, time_grid, METHODS
from ..graph import reaction_network_dot
from ..nse import solve_nse
from ..validation import validate_network, validate_zone


def _load(path):
    return read_xml(path)


def _properties(zone):
    """All of a zone's properties.

    libnucnet zone files keep T9 and rho among the optional properties, so
    showing only ``zone.properties`` prints nothing for them.  Where a name
    appears in both, ``zone.properties`` wins, as in ``Zone.temperature9``.
    """
    return {**zone.optional_properties, **zone.properties}


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
        print(f"zone {i} label={z.label} props={_properties(z)}")
        for name, y in sorted(z.abundances.items()):
            if y >= args.min_abundance:
                print(f"  {name:8s} {y:.8e}")


def cmd_largest_x(args):
    n = _load(args.xml)
    z = n.zone(args.zone_index)
    for name, x in largest_mass_fractions(z, n.species, args.n, args.min_x):
        print(f"{name:8s} {x:.12e}")


def cmd_zone_abundances(args):
    n = _load(args.xml)
    z = n.zone(args.zone_index)
    for name, y in sorted(z.abundances.items()):
        if y >= args.min_abundance:
            print(f"{name} {y:.12e}")


def cmd_zone_properties(args):
    n = _load(args.xml)
    z = n.zone(args.zone_index)
    for k, v in sorted(_properties(z).items()):
        print(k, v)


def cmd_element(args):
    n = _load(args.xml)
    z = n.zone(args.zone_index)
    for el, y in element_abundances(z, n.species).items():
        if el == args.element.lower() or args.element == '*':
            print(el, y)


def cmd_rates(args):
    n = _load(args.xml)
    for r, val in n.reactions.rates(args.t9, args.rho).items():
        if val >= args.min_rate:
            print(f"{r:60s} {val:.12e}")


def cmd_flows(args):
    n = _load(args.xml)
    for r, val in flows(n, args.zone_index, args.t9, args.rho).items():
        if val >= args.min_flow:
            print(f"{r:60s} {val:.12e}")


def cmd_ydot(args):
    n = _load(args.xml)
    for s, val in sorted(ydot(n, args.zone_index, args.t9, args.rho).items()):
        if abs(val) >= args.min_abs:
            print(f"{s:8s} {val:.12e}")


def cmd_net_flows(args):
    n = _load(args.xml)
    z = n.zone(args.zone_index)
    t9 = args.t9 or z.temperature9()
    rho = args.rho or z.density()
    for r, (f, rev, net) in detailed_net_flows(n, z.abundances, t9=t9, rho=rho).items():
        if max(abs(f), abs(rev)) >= args.min_flow:
            print(f"{r:60s} fwd={f:.6e} rev={rev:.6e} net={net:.6e}")


def cmd_charge_flows(args):
    n = _load(args.xml)
    contributions = charge_changing_flows(n, args.zone_index, args.t9, args.rho)
    for r, val in sorted(contributions.items(), key=lambda kv: -abs(kv[1])):
        if abs(val) >= args.min_flow:
            print(f"{r:60s} {val:.12e}")
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
    n = _load(args.xml)
    before = len(n.reactions.reactions)
    n.reactions.remove_duplicates()
    write_xml(n, args.output)
    print(f"removed {before - len(n.reactions.reactions)} duplicate reactions")


def cmd_remove_invalid(args):
    n = _load(args.xml)
    before = len(n.reactions.reactions)
    n.reactions.filter_valid(n.species)
    write_xml(n, args.output)
    print(f"removed {before - len(n.reactions.reactions)} invalid reactions")


def cmd_export_zone_xml(args):
    n = _load(args.xml)
    out = type(n)()
    out.species = n.species
    out.reactions = n.reactions
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
    for i, label, y in species_history(n, args.species):
        print(i, '|'.join(label), f"{y:.12e}")


#: Above this many species ``evolve-zone`` says up front how long one
#: evaluation of the right-hand side takes, and what to do if the integration
#: advances slowly.
LARGE_NETWORK = 1000

#: Settings of ``adaptive_implicit_euler`` that get through explosive burning
#: in a large network; ``evolve-zone`` suggests them when an integration
#: advances slowly.
EXPLOSIVE_STEP_CONTROL = "--method adaptive_implicit_euler --max-change 0.05 --min-mass-fraction 1e-10"

#: ``evolve-zone`` gives its advice when, at a progress report after the
#: first, the integration has covered less than ``SLOW_FRACTION`` of the
#: requested time and the time it has reached grew by less than
#: ``SLOW_GROWTH`` since the previous report.  An integration getting through
#: a fast transient takes steps that grow geometrically, and passes this test;
#: one held to a fixed tiny step does not.
SLOW_FRACTION = 1.0e-6
SLOW_GROWTH = 3.0


def _duration(seconds: float) -> str:
    seconds = int(round(seconds))
    if seconds < 60:
        return f"{seconds} s"
    minutes, seconds = divmod(seconds, 60)
    if minutes < 60:
        return f"{minutes} min {seconds} s"
    hours, minutes = divmod(minutes, 60)
    return f"{hours} h {minutes} min"


class _Progress:
    """Report how far an integration has got, on standard error.

    ``evolve_zone`` calls it after every accepted step with the time reached,
    so it counts the steps and takes the step size from successive times.  It
    prints at most once per ``interval`` seconds of computing time, so a short
    run prints nothing.  If two successive reports find the integration
    crawling (see ``SLOW_FRACTION``), it also says what to try instead, once.
    """

    def __init__(self, times, interval: float, args=None):
        import time
        self.times = times
        self.interval = interval
        self.args = args
        self.start = self.last = time.perf_counter()
        self.reached = float(times[0])
        self.steps = 0
        self.step = 0.0
        self.previous = None
        self.advised = False

    def __call__(self, t: float) -> None:
        import time
        import numpy as np
        self.steps += 1
        self.step = float(t) - self.reached
        self.reached = max(self.reached, float(t))
        now = time.perf_counter()
        if now - self.last < self.interval:
            return
        self.last = now
        done = int(np.searchsorted(self.times, self.reached, side="right")) - 1
        print(f"evolve-zone: t = {self.reached:.4g} s of {self.times[-1]:.4g} s | "
              f"step {self.steps}, dt = {self.step:.2g} s | "
              f"{done} of {len(self.times) - 1} output times | "
              f"after {_duration(now - self.start)}", file=sys.stderr, flush=True)
        covered = self.reached - float(self.times[0])
        span = float(self.times[-1]) - float(self.times[0])
        if (not self.advised and self.previous is not None and span > 0
                and covered < SLOW_FRACTION * span and covered < SLOW_GROWTH * self.previous):
            self.advised = True
            print(f"evolve-zone: hint: {_advice(self.args)}", file=sys.stderr, flush=True)
        self.previous = covered


def _advice(args) -> str:
    """What to try when an integration advances too slowly to finish."""
    cut = "cut the network to the nuclides the burning reaches with --zmax and --amax"
    if getattr(args, "method", "bdf") != "adaptive_implicit_euler":
        return ("the integration is advancing too slowly to finish. Explosive burning "
                f"makes the network very stiff; try {EXPLOSIVE_STEP_CONTROL}, or {cut}.")
    looser = []
    if getattr(args, "max_change", 0.05) < 0.05:
        looser.append("--max-change 0.05")
    if getattr(args, "min_mass_fraction", 1e-10) < 1e-10:
        looser.append("--min-mass-fraction 1e-10")
    if looser:
        return ("the integration is advancing too slowly to finish. Larger steps need "
                f"{' and '.join(looser)}, or {cut}.")
    return f"the integration is advancing too slowly to finish; {cut}."


def _warn_if_large(n, z, thermo) -> None:
    """Say up front how costly a large network is, and what to do if it is slow."""
    import time
    import numpy as np
    from ..solver import rhs
    from ..species import is_massless
    species = [s for s in n.species_names() if not is_massless(s)]
    if len(species) <= LARGE_NETWORK:
        return
    f = rhs(n, species, thermo)
    y0 = np.array([z.get_abundance(s) for s in species])
    f(0.0, y0)  # the first call also evaluates the rate parts that depend on T9 alone
    samples = []
    for _ in range(3):
        t0 = time.perf_counter()
        f(0.0, y0)
        samples.append(time.perf_counter() - t0)
    one = sorted(samples)[1]
    print(f"evolve-zone: note: the network has {len(species)} species and "
          f"{len(n.reactions.reactions)} reactions; one evaluation of the "
          f"right-hand side takes {one:.2g} s here. A stiff integration needs "
          f"thousands of evaluations, and explosive burning far more. If the "
          f"integration advances slowly, use {EXPLOSIVE_STEP_CONTROL}, or cut "
          f"the network with --zmax and --amax.",
          file=sys.stderr, flush=True)


def _run_conditions(args, z):
    """T9, rho, end time and first step of an evolve-zone run, with their sources.

    An option wins over the zone's property, as NucNet Tools reads its zone
    files.  Without either, the end time is 1 s and the integrator chooses the
    first step.  The temperature and the density have no default: a run at an
    assumed temperature gives a plausible but meaningless composition.
    """
    from ..core import (DENSITY_PROPERTIES, END_TIME_PROPERTIES,
                        TEMPERATURE_PROPERTIES, TIME_STEP_PROPERTIES)

    def pick(value, option, names, read):
        if value is not None:
            return value, option
        key, _ = z.find_property(names)
        return (None, None) if key is None else (read(), f"zone property {key}")

    t9, t9_from = pick(args.t9, "--t9", TEMPERATURE_PROPERTIES, z.temperature9)
    rho, rho_from = pick(args.rho, "--rho", DENSITY_PROPERTIES, z.density)
    t1, t1_from = pick(args.t1, "--t1", END_TIME_PROPERTIES, z.end_time)
    dt, dt_from = pick(args.dt, "--dt", TIME_STEP_PROPERTIES, z.time_step)
    missing = []
    if t9 is None:
        missing.append("no temperature: give --t9, or a t9 or t9_0 property in the zone")
    if rho is None:
        missing.append("no density: give --rho, or a rho or rho_0 property in the zone")
    if missing:
        raise ValueError("; ".join(missing))
    if t1 is None:
        t1, t1_from = 1.0, "default"
    return t9, t9_from, rho, rho_from, t1, t1_from, dt, dt_from


def cmd_evolve(args):
    from ..network_limiter import limit_network, select_species
    n = _load(args.xml)
    z = n.zone(args.zone_index)
    t9, t9_from, rho, rho_from, t1, t1_from, dt, dt_from = _run_conditions(args, z)
    first = (f"first step {dt:g} s ({dt_from})" if dt is not None
             else "first step chosen by the integrator")
    print(f"evolve-zone: T9 = {t9:g} ({t9_from}), rho = {rho:g} g/cm^3 ({rho_from}), "
          f"t = {args.t0:g} to {t1:g} s ({t1_from}), {first}", file=sys.stderr, flush=True)
    if args.zmax is not None or args.amax is not None:
        # Keep the nuclides inside the cut and the reactions among them;
        # photons and leptons stay, so photodisintegrations and decays survive.
        limit_network(n, select_species(n, zmax=args.zmax, amax=args.amax))
    times = time_grid(args.t0, t1, args.steps, args.log_time)
    thermo = constant_thermo(t9, rho)
    _warn_if_large(n, z, thermo)
    # The folder is made before the integration, so that a long run never ends
    # unable to save its results.
    folder, created = (None, False) if args.no_output else _output_folder(args.output)
    target = None if folder is None else folder / "output.xml"
    if target is not None:
        print(f"evolve-zone: the results will be saved in {target}", file=sys.stderr, flush=True)
    progress = _Progress(times, args.progress, args) if args.progress > 0 else None
    import time
    start = time.perf_counter()
    try:
        res = evolve_zone(n, z, times, thermo=thermo, method=args.method,
                          rtol=args.rtol, atol=args.atol, max_change=args.max_change,
                          min_mass_fraction=args.min_mass_fraction, progress=progress,
                          first_step=dt)
        elapsed = time.perf_counter() - start
        saved = None
        if target is not None and len(res.time):
            saved = (target, write_evolution_xml(n, res, target, thermo=thermo,
                                                 min_abundance=args.min_abundance))
    finally:
        if created and not any(folder.iterdir()):
            folder.rmdir()
    # A failed solve still carries a full abundance vector, and a positivity
    # projection that had to invent baryon number also reports failure here.
    # The summary therefore opens with the verdict, which also sets the exit code.
    for line in _summary(n, res, args.method, elapsed, saved):
        print(line)
    if not res.success:
        print(f"error: {res.message}", file=sys.stderr)
    return 0 if res.success else 1


def _output_folder(requested):
    """The folder a run saves its results in, and whether it was made now.

    Without ``requested``, a new folder named after the date and time, in the
    current directory.  An ``output.xml`` already in the folder is refused,
    rather than overwritten after the run.
    """
    from datetime import datetime
    if requested:
        folder = Path(requested)
    else:
        stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        folder, k = Path(f"evolve-zone_{stamp}"), 2
        while folder.exists():
            folder, k = Path(f"evolve-zone_{stamp}_{k}"), k + 1
    if (folder / "output.xml").exists():
        raise FileExistsError(f"{folder / 'output.xml'} already exists; choose another "
                              "--output folder, or move the file away")
    created = not folder.exists()
    folder.mkdir(parents=True, exist_ok=True)
    return folder, created


def _summary(network, res, method, elapsed, saved):
    """The lines that end a run: the verdict, the cost, conservation and the result."""
    import numpy as np
    from ..species import Species
    lines = [f"evolve-zone: finished in {_duration(elapsed)}: "
             + ("success" if res.success else "FAILED")]
    if not len(res.time):
        lines.append("  the integration produced no output, so nothing was saved")
        return lines
    counts = [f"{res.steps} steps" + (f" ({res.rejected_steps} rejected)"
                                       if method == "adaptive_implicit_euler" else "")]
    if res.nfev:
        counts.append(f"{res.nfev} right-hand sides")
    if res.njev:
        counts.append(f"{res.njev} Jacobians")
    lines.append("  " + ", ".join(counts))
    nuclides = [network.species.get(s) or Species.parse(s) for s in res.species]
    a = np.array([sp.a for sp in nuclides], dtype=float)
    charge = np.array([sp.z for sp in nuclides], dtype=float)
    first, last = np.asarray(res.y[0], dtype=float), np.asarray(res.y[-1], dtype=float)
    xsum = float(a @ last)
    lines.append(f"  sum X = {xsum:.12g} (changed by {xsum - float(a @ first):.1e} during "
                 f"the run), Ye = {float(charge @ last):.6g}")
    x = a * last
    order = [i for i in np.argsort(-x, kind="stable")[:10] if x[i] > 0]
    lines.append(f"  largest mass fractions at t = {float(res.time[-1]):.4g} s:")
    for row in range(0, len(order), 3):
        lines.append("    " + "   ".join(f"{res.species[i]:<6} {x[i]:.3e}"
                                         for i in order[row:row + 3]))
    if saved is not None:
        path, zones = saved
        lines.append(f"  saved: {path} ({zones} zones, one per output time"
                     + (")" if res.success else "; the run failed, so do not trust them)"))
        lines.append(f"  final abundances: nucnetpy zone-abundances {path} --zone-index {zones - 1}")
    return lines


def cmd_energy(args):
    from ..analysis import nuclear_energy_generation_rate
    from ..constants import AVOGADRO, MEV_TO_ERG
    n = _load(args.xml)
    # The mass-excess form is exact and consistent with the equilibrium
    # solvers.  The Q-value form uses the rates' own Q-values, which matters
    # for a network whose nuclides carry no mass excesses.  Both are printed in
    # erg/g/s; energy_generation_rate itself returns MeV per nucleon per second.
    from_masses = nuclear_energy_generation_rate(n, args.zone_index, args.t9, args.rho)
    from_q = energy_generation_rate(n, args.zone_index, args.t9, args.rho) * AVOGADRO * MEV_TO_ERG
    print(f"from mass excesses  {from_masses:.12e} erg/g/s")
    print(f"from reaction Q     {from_q:.12e} erg/g/s")


def cmd_dot(args):
    n = _load(args.xml)
    dot = reaction_network_dot(n, min_rate=args.min_rate, t9=args.t9, rho=args.rho)
    if args.output:
        Path(args.output).write_text(dot)
    else:
        print(dot)


def cmd_validate(args):
    n = _load(args.xml)
    issues = validate_network(n, strict=args.strict)
    for z in n.zones[:args.max_zones]:
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
    return 0 if res.success else 1


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
    return 0 if res.success else 1


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

    def command(name, func, help):
        """Add a command that reads one network XML file, given first."""
        s = sub.add_parser(name, help=help, description=help)
        s.add_argument('xml', help='network XML file')
        s.set_defaults(func=func)
        return s

    def zone_and_conditions(s):
        """Options shared by the commands that evaluate one zone of the network."""
        s.add_argument('--zone-index', type=int, default=0)
        s.add_argument('--t9', type=float)
        s.add_argument('--rho', type=float)

    # Inspecting a network and its zones
    command('summary', cmd_summary, 'count the species, reactions and zones, and check the network')

    s = command('print-output', cmd_print_output, 'print the properties and abundances of each zone')
    s.add_argument('--max-zones', type=int, default=3)
    s.add_argument('--min-abundance', type=float, default=0.0)

    s = command('largest-x', cmd_largest_x, 'list the largest mass fractions in a zone')
    s.add_argument('-n', type=int, default=10)
    s.add_argument('--zone-index', type=int, default=0)
    s.add_argument('--min-x', type=float, default=0.0)

    s = command('zone-abundances', cmd_zone_abundances, 'list the abundances in a zone')
    s.add_argument('--zone-index', type=int, default=0)
    s.add_argument('--min-abundance', type=float, default=0.0)

    s = command('zone-properties', cmd_zone_properties, 'list the properties of a zone')
    s.add_argument('--zone-index', type=int, default=0)

    s = command('element-abundances', cmd_element, 'list the abundances in a zone summed by element')
    s.add_argument('element', help="element symbol, such as fe; '*' lists every element")
    s.add_argument('--zone-index', type=int, default=0)

    # Rates, flows and derived quantities
    s = command('rates', cmd_rates, 'list the rate of every reaction at a given T9 and rho')
    s.add_argument('--t9', type=float, default=1.0)
    s.add_argument('--rho', type=float, default=1.0)
    s.add_argument('--min-rate', type=float, default=0.0)

    s = command('flows', cmd_flows, 'list the reaction flows in a zone')
    zone_and_conditions(s)
    s.add_argument('--min-flow', type=float, default=0.0)

    s = command('ydot', cmd_ydot, 'list dY/dt for every species in a zone')
    zone_and_conditions(s)
    s.add_argument('--min-abs', type=float, default=0.0)

    s = command('net-flows', cmd_net_flows, 'list the forward, detailed-balance reverse and net flows in a zone')
    zone_and_conditions(s)
    s.add_argument('--min-flow', type=float, default=0.0)

    s = command('charge-flows', cmd_charge_flows, "list each reaction's contribution to dYe/dt in a zone")
    zone_and_conditions(s)
    s.add_argument('--min-flow', type=float, default=0.0)

    s = command('timescales', cmd_timescales, 'list the shortest species timescales, Y/|dY/dt|, in a zone')
    zone_and_conditions(s)
    s.add_argument('-n', type=int, default=20)

    s = command('entropy-generation', cmd_entropy, 'print the entropy generation rate of a zone, dS/dt per nucleon in k_B/s')
    zone_and_conditions(s)

    # Checking and editing a reaction set
    s = command('conservation', cmd_conservation, 'list the reactions that do not conserve A and Z')
    s.add_argument('--max', type=int, default=20)

    s = command('remove-duplicates', cmd_remove_duplicates, 'write a copy of the network without duplicate reactions')
    s.add_argument('output', help='XML file to write')

    s = command('remove-invalid', cmd_remove_invalid, 'write a copy of the network without the reactions that do not conserve A and Z')
    s.add_argument('output', help='XML file to write')

    s = command('export-zone-xml', cmd_export_zone_xml, 'write a copy of the network holding only one zone')
    s.add_argument('output', help='XML file to write')
    s.add_argument('--zone-index', type=int, default=0)

    s = command('reactions-latex', cmd_reactions_latex, 'write the reactions and their Q values as a LaTeX table')
    s.add_argument('output', help='LaTeX file to write')

    s = command('species-history', cmd_species_history, 'list the abundance of one species in every zone')
    s.add_argument('species', help='species name, such as fe56')

    # Evolution
    s = command('evolve-zone', cmd_evolve, 'integrate a zone in time at constant T9 and rho, and save its history')
    s.add_argument('--zone-index', type=int, default=0)
    s.add_argument('--t0', type=float, default=0.0, help='start time in seconds (default 0)')
    s.add_argument('--t1', type=float, help="end time in seconds (default: the zone's tend property, else 1)")
    s.add_argument('--steps', type=int, default=50, help='number of output times (default 50)')
    s.add_argument('--t9', type=float, help="temperature in units of 10^9 K (default: the zone's t9 or t9_0 property)")
    s.add_argument('--rho', type=float, help="density in g/cm^3 (default: the zone's rho or rho_0 property)")
    s.add_argument('--dt', type=float, help="first time step in seconds (default: the zone's dt property; without either, the integrator chooses)")
    s.add_argument('--method', default='bdf', type=str.lower, choices=METHODS)
    s.add_argument('--log-time', action='store_true')
    s.add_argument('--min-abundance', type=float, default=0.0, metavar='Y', help='leave out of the saved file the abundances smaller than Y (default 0: every positive abundance is saved)')
    s.add_argument('--zmax', type=int, help='keep only nuclides with Z <= ZMAX, and the reactions among them')
    s.add_argument('--amax', type=int, help='keep only nuclides with A <= AMAX, and the reactions among them')
    s.add_argument('--rtol', type=float, default=1e-6, help='relative tolerance of the SciPy integrators (default 1e-6)')
    s.add_argument('--atol', type=float, default=1e-30, help='absolute tolerance on the abundances (default 1e-30, which makes the relative tolerance govern every abundance; 1e-12 is often enough and much faster)')
    s.add_argument('--max-change', type=float, default=0.005, metavar='FRACTION', help='adaptive_implicit_euler only: the largest fractional change of an abundance in one step (default 0.005; 0.05 suits explosive burning)')
    s.add_argument('--min-mass-fraction', type=float, default=1e-15, metavar='X', help='adaptive_implicit_euler only: abundances with a smaller mass fraction do not limit the step (default 1e-15)')
    s.add_argument('--progress', type=float, default=30.0, metavar='SECONDS', help='report progress on standard error every SECONDS of computing time (default 30; 0 turns it off)')
    saving = s.add_mutually_exclusive_group()
    saving.add_argument('--output', metavar='FOLDER', help='folder to save the results in, as output.xml: one zone per output time, in the libnucnet format of NucNet Tools (default: a new folder evolve-zone_DATE_TIME in the current directory)')
    saving.add_argument('--no-output', action='store_true', help='save nothing; only print the summary')

    s = command('energy-generation', cmd_energy, 'print the nuclear energy generation rate of a zone, in erg/g/s')
    zone_and_conditions(s)

    s = command('net-dot', cmd_dot, 'write the network as a Graphviz DOT graph')
    s.add_argument('-o', '--output')
    s.add_argument('--t9', type=float, default=1.0)
    s.add_argument('--rho', type=float, default=1.0)
    s.add_argument('--min-rate', type=float, default=0.0)

    s = command('validate', cmd_validate, 'check the network and its zones for problems')
    s.add_argument('--strict', action='store_true')
    s.add_argument('--max', type=int, default=50)
    s.add_argument('--max-zones', type=int, default=10)

    # Equilibrium
    s = command('nse', cmd_nse, 'solve for nuclear statistical equilibrium')
    s.add_argument('--t9', type=float, required=True)
    s.add_argument('--rho', type=float, required=True)
    s.add_argument('--ye', type=float, required=True)
    s.add_argument('--coulomb', action='store_true', help='apply Bravo & Garcia-Senz Coulomb corrections')
    s.add_argument('--min-abundance', type=float, default=0.0)
    s.add_argument('--min-x', type=float, default=0.0)

    s = command('qse', cmd_qse, 'solve for constrained cluster equilibrium (libnuceq)')
    s.add_argument('--t9', type=float, required=True)
    s.add_argument('--rho', type=float, required=True)
    s.add_argument('--ye', type=float, required=True)
    s.add_argument('--cluster', action='append', metavar='SP1,SP2,...:Y', help='cluster species and constrained total abundance; repeatable')
    s.add_argument('--coulomb', action='store_true')
    s.add_argument('--min-abundance', type=float, default=0.0)
    s.add_argument('--min-x', type=float, default=0.0)

    # Separate JINA nuclide, reaction and zone files
    help = 'summarise separate JINA nuclide and reaction XML files'
    s = sub.add_parser('jina-summary', help=help, description=help)
    s.add_argument('nuclides_xml', help='nuclide XML file')
    s.add_argument('reactions_xml', help='reaction XML file')
    s.add_argument('--zones-xml')
    s.add_argument('--show-invalid', action='store_true')
    s.add_argument('--max-invalid', type=int, default=20)
    s.set_defaults(func=cmd_jina_summary)

    help = 'combine JINA nuclide and reaction XML files into one nucnetpy XML file'
    s = sub.add_parser('jina-combine', help=help, description=help)
    s.add_argument('nuclides_xml', help='nuclide XML file')
    s.add_argument('reactions_xml', help='reaction XML file')
    s.add_argument('output_xml', help='combined XML file to write')
    s.add_argument('--zones-xml')
    s.set_defaults(func=cmd_jina_combine)
    return p


def main(argv=None):
    import xml.etree.ElementTree as ET
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except (OSError, ET.ParseError, IndexError, ValueError) as exc:
        # Input mistakes -- a missing or malformed file, a zone index that
        # does not exist, an unknown option value -- get one line, not a
        # traceback.  Anything else is a bug and keeps its traceback.
        print(f"nucnetpy {args.cmd}: error: {exc}", file=sys.stderr)
        return 1
