"""Build the network and the initial composition for the reproduction.

The paper ran NucNet Tools on JINA ReacLib rates for every nuclide from
hydrogen to nickel, starting from the solar composition of Lodders (2003).
This script assembles the same two inputs for NucNetPy:

``data/network_reaclib_default2_20220329_z28.xml.gz``
    Every ReacLib rate whose nuclides all have Z <= 28, written in the
    libnucnet XML format that ``nucnetpy.read_xml`` reads.  The rates come
    from the ReacLib snapshot of 29 March 2022 that pynucastro distributes
    (the paper was submitted on 15 March 2022).  Masses and spins are from
    AME2020 and NUBASE2020, and partition functions from Rauscher and
    Thielemann (2000), all as distributed by pynucastro.

``data/initial_composition_lodders2003.json``
    The Lodders (2003) isotopic composition, rescaled to the present-day
    photospheric X = 0.7491, Y = 0.2377 and Z = 0.0133 that the paper quotes.

``data/provenance.json``
    Versions, checksums and counts that identify both inputs.

It needs pynucastro (for ReacLib and the nuclear data) and the Lodders (2003)
table ``sollo03.dat`` from the StarFit package::

    pip install pynucastro==3.2.0
    pip download starfit==0.19.1 --no-deps --no-binary :all:
    tar xzf starfit-0.19.1.tar.gz starfit-0.19.1/src/starfit/data/ref/sollo03.dat
    python build_inputs.py --lodders starfit-0.19.1/src/starfit/data/ref/sollo03.dat

The outputs are committed, so the calculation itself needs only NucNetPy.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import math
import sys
from collections import Counter, deque
from pathlib import Path

from common import COMPOSITION_FILE, DATA, NETWORK_FILE, PROVENANCE_FILE

LIBRARY = "reaclib_default2_20220329"
MAX_Z = 28
PARTITION_T9_MAX = 10.0
#: Lodders (2003), Table 2 and section 4: present-day photospheric values,
#: quoted by the paper for its initial composition.
X_PHOTO, Y_PHOTO, Z_PHOTO = 0.7491, 0.2377, 0.0133

ELEMENTS = ["n", "h", "he", "li", "be", "b", "c", "n", "o", "f", "ne", "na",
            "mg", "al", "si", "p", "s", "cl", "ar", "k", "ca", "sc", "ti", "v",
            "cr", "mn", "fe", "co", "ni"]


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def species_name(z: int, a: int) -> str:
    """libnucnet name of a nuclide: n, h1, he4, al26, ..."""
    return "n" if (z, a) == (0, 1) else f"{ELEMENTS[z]}{a}"


def participants(rate):
    """Reactants and products of a pynucastro rate, with photons and leptons.

    libnucnet files list photons and leptons explicitly, so that charge is
    conserved and a capture can be paired with its photodisintegration:

    * beta-minus decay gains an electron and an antineutrino as products;
    * beta-plus decay (or a laboratory electron-capture decay) gains a
      positron and a neutrino;
    * an electron capture (ReacLib label ``ec``) gains an electron reactant
      and a neutrino product;
    * a radiative capture into one nucleus gains a photon product;
    * the photodisintegration that ReacLib derives from such a capture (its
      ``v`` flag) gains a photon reactant.
    """
    reactants = [species_name(n.Z, n.A) for n in rate.reactants]
    products = [species_name(n.Z, n.A) for n in rate.products]
    if rate.weak_type == "beta_neg":
        products += ["electron", "anti-neutrino_e"]
    elif rate.weak_type == "beta_pos":
        products += ["positron", "neutrino_e"]
    elif rate.weak_type == "electron_capture":
        reactants += ["electron"]
        products += ["neutrino_e"]
    elif len(products) == 1 and len(reactants) > 1:
        products += ["gamma"]
    elif len(reactants) == 1 and len(products) > 1 and rate.derived_from_inverse:
        reactants += ["gamma"]
    return reactants, products


def derive_missing_masses(nuclides, rates, known):
    """Mass excesses for nuclides AME2020 does not list, from ReacLib Q values.

    About a third of the nuclides ReacLib holds up to nickel lie near the drip
    lines, beyond the mass evaluation.  Their masses do not enter the
    calculation (the rates are fitted, not derived), but a libnucnet file
    needs one, and a placeholder of zero would make them look as bound as
    carbon-12.  Each is obtained from the Q value of a ReacLib rate that
    links it to nuclides whose masses are known, working outwards from the
    evaluated ones.  ReacLib quotes every Q value, weak ones included, as the
    difference of atomic mass excesses, so
    ``Q = sum(reactant excesses) - sum(product excesses)`` holds for all of
    them; strong and electromagnetic rates are used first, and weak rates
    only for the few nuclides that nothing else reaches.
    """
    masses = dict(known)
    source = {}
    for group in ([r for r in rates if not r.weak_type], list(rates)):
        pending = deque(group)
        progress = True
        while progress:
            progress = False
            for _ in range(len(pending)):
                r = pending.popleft()
                names = [(species_name(n.Z, n.A), +1) for n in r.reactants] + \
                        [(species_name(n.Z, n.A), -1) for n in r.products]
                unknown = {n for n, _ in names if n not in masses}
                if len(unknown) != 1:
                    if unknown:
                        pending.append(r)
                    continue
                (target,) = unknown
                coeff = sum(sign for n, sign in names if n == target)
                if coeff == 0:
                    continue
                rest = sum(sign * masses[n] for n, sign in names if n != target)
                masses[target] = (r.Q - rest) / coeff
                source[target] = f"ReacLib Q of {r.fname}"
                progress = True
    missing = sorted(set(nuclides) - set(masses))
    return masses, source, missing


def build_network(lib_path: Path):
    import pynucastro as pyna
    from nucnetpy.core import Network
    from nucnetpy.reactions import RateFit, Reaction, ReactionNetwork
    from nucnetpy.species import Species

    library = pyna.ReacLibLibrary(libfile=LIBRARY)
    rates = [r for r in library.get_rates()
             if all(n.Z <= MAX_Z for n in r.reactants + r.products)]
    nuclei = {}
    for r in rates:
        for n in r.reactants + r.products:
            nuclei[species_name(n.Z, n.A)] = n

    known = {name: float(n.dm) for name, n in nuclei.items() if n.dm is not None}
    masses, mass_source, missing = derive_missing_masses(nuclei, rates, known)
    if missing:
        raise RuntimeError(f"no mass could be found for {missing}")

    species = {}
    for name, n in sorted(nuclei.items(), key=lambda kv: (kv[1].Z, kv[1].A)):
        spin = None if n.spin_states is None else (float(n.spin_states) - 1.0) / 2.0
        partition = {}
        pf = n.partition_function
        if pf is not None:
            # pynucastro stores ln G, G being normalised to the ground state.
            # Only the points up to T9 = 10 are kept: the grid of the JINA
            # nuclide file, and well beyond the 3.6 GK of silicon burning.
            partition = {float(t): float(math.exp(lg))
                         for t, lg in zip(pf.T9_points, pf.log_pf_data)
                         if t <= PARTITION_T9_MAX}
        src = "AME2020" if name in known else mass_source[name]
        species[name] = Species(name, int(n.Z), int(n.A), mass_excess=masses[name],
                                spin=spin, source=src, partition=partition)

    reactions = []
    for r in rates:
        reactants, products = participants(r)
        fits = [RateFit(list(s.a), label=s.labelprops.strip()) for s in r.sets]
        reactions.append(Reaction.from_names(reactants, products, rate_fits=fits,
                                             q_value=float(r.Q),
                                             source=r.sets[0].labelprops[:4].strip()))
    net = Network(species=species, reactions=ReactionNetwork(reactions),
                  metadata={"library": LIBRARY})
    keys = Counter(r.key for r in reactions)
    duplicates = [k for k, c in keys.items() if c > 1]
    if duplicates:
        raise RuntimeError(f"duplicate reactions: {duplicates}")
    stats = {
        "library": LIBRARY,
        "library_sha256": sha256(lib_path.read_bytes()),
        "pynucastro": pyna.__version__,
        "max_z": MAX_Z,
        "nuclides": len(species),
        "reactions": len(reactions),
        "reverse_rates_flagged_v": sum(1 for r in rates if r.derived_from_inverse),
        "weak_rates": dict(Counter(r.weak_type for r in rates if r.weak_type)),
        "electron_captures": sorted(r.fname for r in rates
                                    if r.weak_type == "electron_capture"),
        "masses_from_ame2020": len(known),
        "masses_from_reaclib_q": len(mass_source),
        "spins_unknown": sum(1 for s in species.values() if s.spin is None),
        "partition_tables": sum(1 for s in species.values() if s.partition),
        "partition_source": "Rauscher & Thielemann (2000) FRDM, via pynucastro",
        "partition_t9_max": PARTITION_T9_MAX,
    }
    return net, stats


def build_composition(lodders_path: Path, network_species):
    """Lodders (2003) isotopes, rescaled to the paper's X, Y and Z.

    ``sollo03.dat`` holds Lodders' protosolar mass fractions (X = 0.711,
    Y = 0.274).  Lodders derives the present-day photosphere from them by
    gravitational settling, which changes the hydrogen, helium and metal
    totals but not the proportions within each group.  Each group is
    therefore scaled to the photospheric total.  Nuclides heavier than nickel
    (4e-6 of the mass) are outside the network; the metals are renormalised
    so that Z stays 0.0133.  Lodders' rounded totals sum to 1.0001, so the
    whole composition is finally divided by that sum to make the mass
    fractions add up to one.
    """
    raw = {}
    for line in lodders_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith(";"):
            continue
        name, value = line.split()[:2]
        raw[name.lower()] = float(value)
    hydrogen = {k: v for k, v in raw.items() if k in ("h1", "h2")}
    helium = {k: v for k, v in raw.items() if k in ("he3", "he4")}
    metals = {k: v for k, v in raw.items() if k not in hydrogen and k not in helium}
    kept = {k: v for k, v in metals.items() if k in network_species}
    dropped = {k: v for k, v in metals.items() if k not in network_species}

    def scaled(group, total):
        s = sum(group.values())
        return {k: v * total / s for k, v in group.items()}

    x = {**scaled(hydrogen, X_PHOTO), **scaled(helium, Y_PHOTO), **scaled(kept, Z_PHOTO)}
    total = sum(x.values())
    x = {k: v / total for k, v in x.items()}
    record = {
        "comment": "Initial mass fractions: Lodders (2003, ApJ 591, 1220) isotopic "
                   "composition from StarFit's sollo03.dat, with hydrogen, helium "
                   "and metals each scaled to the present-day photospheric "
                   "X = 0.7491, Y = 0.2377, Z = 0.0133 that the paper quotes, then "
                   "divided by their sum, 1.0001, so that the mass fractions add "
                   "up to one.",
        "source_file": lodders_path.name,
        "source_sha256": sha256(lodders_path.read_bytes()),
        "source_package": "starfit 0.19.1 (src/starfit/data/ref/sollo03.dat, GPL-3.0)",
        "protosolar_totals": {"X": sum(hydrogen.values()), "Y": sum(helium.values()),
                              "Z": sum(metals.values())},
        "photospheric_totals": {"X": X_PHOTO, "Y": Y_PHOTO, "Z": Z_PHOTO},
        "dropped_beyond_nickel_protosolar_mass": sum(dropped.values()),
        "mass_fractions": dict(sorted(x.items(), key=lambda kv: -kv[1])),
    }
    return record


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lodders", required=True, type=Path,
                    help="path to StarFit's sollo03.dat (Lodders 2003)")
    args = ap.parse_args()

    import pynucastro as pyna
    lib_path = Path(pyna.__file__).parent / "data" / LIBRARY
    net, stats = build_network(lib_path)
    print(f"{stats['nuclides']} nuclides, {stats['reactions']} reactions", flush=True)

    from nucnetpy.io.xml import write_xml
    DATA.mkdir(parents=True, exist_ok=True)
    tmp = io.BytesIO()
    write_xml(net, tmp)
    # A fixed mtime keeps the compressed file identical between runs.
    with open(NETWORK_FILE, "wb") as raw, \
            gzip.GzipFile(filename=NETWORK_FILE.stem, mode="wb", fileobj=raw,
                          mtime=0) as gz:
        gz.write(tmp.getvalue())
    stats["xml_sha256"] = sha256(tmp.getvalue())

    composition = build_composition(args.lodders, set(net.species))
    COMPOSITION_FILE.write_text(json.dumps(composition, indent=1) + "\n")

    provenance = {
        "comment": "Inputs of the NucNetPy reproduction of Abdallah & Ershaidat "
                   "(2024, Jordan J. Phys. 17, 273). Regenerate with build_inputs.py.",
        "network_file": NETWORK_FILE.name,
        "network": stats,
        "composition_file": COMPOSITION_FILE.name,
        "composition_source_sha256": composition["source_sha256"],
        "python": sys.version.split()[0],
    }
    PROVENANCE_FILE.write_text(json.dumps(provenance, indent=1) + "\n")
    print(f"wrote {NETWORK_FILE.name}, {COMPOSITION_FILE.name}, {PROVENANCE_FILE.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
