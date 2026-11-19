"""Fast, permissive NucNet/JINA/libnucnet XML reader and writer.

This module is tuned for real JINA XML databases, where files may contain
thousands of nuclides and tens of thousands of reactions.  It supports:

* separate ``<nuclear_data>`` and ``<reaction_data>`` files,
* combined ``<libnucnet_input>`` output files,
* zone files with mass fractions or abundances,
* ReacLib/JINA ``non_smoker_fit`` blocks with ``a1``--``a7`` coefficients,
* single-rate beta-decay style records,
* partition-function tables.
"""
from __future__ import annotations

import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from typing import List, Optional

from ..core import Network, Zone
from ..species import Species, normalize_species_name, species_from_za, abundance_from_mass_fraction
from ..reactions import Reaction, ReactionNetwork, ReactionParticipant, RateFit, TabularRate


def _tag(el) -> str:
    return el.tag.split('}', 1)[-1].lower() if isinstance(el.tag, str) else ''


def _text(el, default: str = '') -> str:
    return (el.text or default).strip() if el is not None else default


def _attr(el, *names, default=None):
    if el is None:
        return default
    for n in names:
        if n in el.attrib:
            return el.attrib[n]
    lower = {k.lower(): v for k, v in el.attrib.items()}
    for n in names:
        if n.lower() in lower:
            return lower[n.lower()]
    return default


def _child(el, *names):
    names = {n.lower() for n in names}
    for c in list(el):
        if _tag(c) in names:
            return c
    return None


def _children(el, *names):
    names = {n.lower() for n in names}
    return [c for c in list(el) if _tag(c) in names]


def _descendants(el, *names):
    names = {n.lower() for n in names}
    for c in el.iter():
        if c is not el and _tag(c) in names:
            yield c


def _float_or_none(x):
    try:
        if x is None or str(x).strip() == '' or str(x).lower() in {'nan', 'none'}:
            return None
        return float(x)
    except Exception:
        return None


def _int_or_none(x):
    v = _float_or_none(x)
    return None if v is None else int(v)


def read_xml(path: str | Path, reverse_partition_functions: bool = True) -> Network:
    """Read a libnucnet-style XML file holding nuclides, reactions and zones.

    With ``reverse_partition_functions`` (the default) the library reverse
    rates are multiplied by the partition functions they are fitted without;
    see :func:`nucnetpy.apply_reverse_partition_functions`.
    """
    return read_xml_root(ET.parse(path).getroot(), source=str(path),
                         reverse_partition_functions=reverse_partition_functions)


read_network_xml = read_xml


def read_xml_root(root, source: str = '<xml-root>',
                  reverse_partition_functions: bool = True) -> Network:
    net = Network(metadata={"source_file": str(source)})
    described = _parse_species(root, net)
    _parse_reactions(root, net)
    _parse_zones(root, net)
    _record_species_without_nuclear_data(net, described)
    if reverse_partition_functions:
        from ..detailed_balance import apply_reverse_partition_functions
        apply_reverse_partition_functions(net)
    return net


def _record_species_without_nuclear_data(net: Network, described) -> None:
    """Create and record the species that no nuclide record describes.

    A reaction or a zone can name a nuclide that the nuclide records leave out.
    :meth:`Network.ensure_species` then gives it a placeholder mass excess of
    zero.  That is not a physical value: in an equilibrium calculation it makes
    an unbound nuclide look as tightly bound as the most stable one.  The names
    are recorded in ``metadata["species_without_nuclear_data"]`` so that
    validation can report them, the equilibrium solvers can leave them out,
    and :func:`write_xml` does not write them out as if they had nuclear data.

    ``described`` holds the species that a nuclide record with data defined.
    Photons and leptons are participants, not nuclides, and are never recorded.
    """
    from ..species import is_massless
    names = set(net.species)
    names.update(net.reactions.species_names())
    for zone in net.zones:
        names.update(zone.abundances)
    placeholders = set(net.species_without_nuclear_data())
    for name in names:
        key = normalize_species_name(name)
        if key in described or is_massless(key):
            continue
        try:
            net.ensure_species(key)
        except Exception:
            continue
        placeholders.add(key)
    if placeholders:
        net.metadata["species_without_nuclear_data"] = ",".join(sorted(placeholders))


def _parse_species(root, net: Network) -> set:
    """Parse true nuclide/species records, avoiding zone mass-fraction records.

    Returns the names of the species defined by a record that carries nuclear
    data.  A bare reference such as ``<nuc name="h1" count="2"/>`` names a
    species without describing it, so a species known only from one is left
    out and later recorded as having no nuclear data.
    """
    described = set()
    for el in root.iter():
        tg = _tag(el)
        if tg not in {"nuclide", "nuc", "species"}:
            continue
        parent_like_zone_record = (_child(el, 'x', 'y', 'abundance') is not None) and (_child(el, 'mass_excess', 'spin', 'partf_table') is None)
        # Zone abundance records often use <nuclide name=...><z>..<a>..<x>..</x></nuclide>.
        # Do not use these to overwrite mass-table entries.
        if parent_like_zone_record and _attr(el, 'name', default=None) is not None:
            continue

        name = _attr(el, "name", "species", "id")
        z = _attr(el, "z", "Z")
        a = _attr(el, "a", "A")
        mass_excess = _attr(el, "mass_excess", "mass-excess", "massExcess", "mass excess", default=None)
        spin = _attr(el, "spin", "spin_factor", default=None)
        source = _attr(el, 'source', default=None)

        if name is None:
            ce = _child(el, 'name')
            name = _text(ce) if ce is not None else None
        if z is None:
            z = _text(_child(el, 'z'), None)
        if a is None:
            a = _text(_child(el, 'a'), None)
        if source is None:
            source = _text(_child(el, 'source'), None)
        if mass_excess is None:
            mass_excess = _text(_child(el, 'mass_excess', 'mass-excess', 'massExcess'), None)
        if spin is None:
            spin = _text(_child(el, 'spin'), None)

        # Avoid overwriting rich mass-table records with lightweight reaction-side
        # references such as <nuc name="h1" count="2"/>.
        lightweight = bool(name and z is None and a is None and mass_excess is None
                           and spin is None
                           and _child(el, 'partf_table', 'partition_function',
                                      'partition-function') is None)
        if lightweight and normalize_species_name(name) in net.species:
            continue

        try:
            kwargs = {"mass_excess": float(mass_excess or 0.0)}
            if spin is not None:
                sv = _float_or_none(spin)
                if sv is not None:
                    kwargs['spin'] = sv
            if source is not None:
                kwargs['source'] = source
            zi = _int_or_none(z)
            ai = _int_or_none(a)
            if name:
                sp = Species.parse(name, **kwargs)
                if zi is not None and ai is not None:
                    sp = Species(sp.name, zi, ai, **kwargs)
            elif zi is not None and ai is not None:
                sp = species_from_za(zi, ai, **kwargs)
            else:
                continue

            part = {}
            part_table = _child(el, 'partf_table', 'partition_function_table', 'partition-functions')
            search_root = part_table if part_table is not None else el
            for pf in search_root.iter():
                if pf is search_root and _tag(pf) not in {'point', 'row', 'partition_function', 'partition-function', 'partf'}:
                    continue
                if _tag(pf) not in {"partition_function", "partition-function", "partf", "partition", "point", "row"}:
                    continue
                t9 = _attr(pf, "t9", "T9", "temperature", default=None)
                val = _attr(pf, "value", "g", "partf", "partition_function", default=None)
                log_val = None
                if t9 is None:
                    t9 = _text(_child(pf, 't9', 'T9', 'temperature'), None)
                if val is None:
                    val = _text(_child(pf, 'partf', 'partition_function', 'value', 'g'), None)
                log_child = _child(pf, 'log10_partf')
                if val is None and log_child is not None:
                    log_val = _text(log_child, None)
                    val = log_val
                try:
                    if t9 is not None and val is not None:
                        fv = float(val)
                        if log_val is not None:
                            fv = 10.0 ** fv
                        part[float(t9)] = fv
                except Exception:
                    pass
            if part:
                sp = Species(sp.name, sp.z, sp.a, mass_excess=sp.mass_excess, spin=sp.spin, source=sp.source, partition=part)
            net.add_species(sp)
            if not lightweight:
                described.add(normalize_species_name(sp.name))
        except Exception:
            continue
    return described


def _participant_name(el) -> Optional[str]:
    v = _attr(el, 'name', 'species', 'nuc', 'id')
    if v:
        return v
    txt = _text(el)
    return txt or None


def _participants_from_names(names: List[str]) -> List[ReactionParticipant]:
    c = Counter(normalize_species_name(n) for n in names if str(n).strip())
    return [ReactionParticipant(k, v) for k, v in sorted(c.items())]


def _participants_from_side_text(text: str) -> List[ReactionParticipant]:
    return _participants_from_names([t for t in str(text).replace('+', ' ').split() if t])


def _parse_reactions(root, net: Network):
    reactions: List[Reaction] = []
    for el in root.iter():
        if _tag(el) != 'reaction':
            continue
        reactant_names: List[str] = []
        product_names: List[str] = []
        source = ''
        for c in list(el):
            tg = _tag(c)
            if tg == 'source':
                source = _text(c)
            elif tg == 'reactant':
                name = _participant_name(c)
                if name:
                    cnt = int(float(_attr(c, 'count', 'number', default='1') or 1))
                    reactant_names.extend([name] * cnt)
            elif tg == 'product':
                name = _participant_name(c)
                if name:
                    cnt = int(float(_attr(c, 'count', 'number', default='1') or 1))
                    product_names.extend([name] * cnt)
            elif tg in {'reactants', 'products'}:
                target = reactant_names if tg == 'reactants' else product_names
                if list(c):
                    for q in list(c):
                        if _tag(q) in {'nuc', 'nuclide', 'species', 'reactant', 'product'}:
                            name = _participant_name(q)
                            if name:
                                cnt = int(float(_attr(q, 'count', 'number', default='1') or 1))
                                target.extend([name] * cnt)
                else:
                    target.extend([t for t in _text(c).replace('+', ' ').split() if t])

        if not reactant_names or not product_names:
            s = _attr(el, 'string', 'reaction', default=None) or _text(_child(el, 'string'), '')
            if '->' in s:
                left, right = s.split('->', 1)
                reactants = _participants_from_side_text(left)
                products = _participants_from_side_text(right)
            else:
                continue
        else:
            reactants = _participants_from_names(reactant_names)
            products = _participants_from_names(product_names)

        fits: List[RateFit] = []
        constant = _attr(el, 'rate', 'constant_rate', default=None)
        for single in _children(el, 'single_rate', 'single-rate'):
            v = _float_or_none(_text(single))
            if v is not None:
                constant = v

        # non_smoker_fit may contain a1..a7 directly or one/more nested <fit> records.
        # Some libnucnet test files place a coefficient-style <single_rate>
        # inside <rate_data>, so scan descendants, not only direct children.
        for container in _descendants(el, 'non_smoker_fit', 'non-smoker-fit', 'rate_fit', 'rate-fit', 'reaclib', 'single_rate', 'single-rate'):
            if _tag(container) in {'single_rate', 'single-rate'} and not any(_child(container, f'a{i}') is not None for i in range(1,8)):
                continue
            fit_records = [f for f in list(container) if _tag(f) == 'fit'] or [container]
            for fit in fit_records:
                coeffs = []
                # JINA/libnucnet ReacLib XML uses a1..a7.  Compact files may use a0..a6.
                if _child(fit, 'a1') is not None or _attr(fit, 'a1', default=None) is not None:
                    keys = ['a1','a2','a3','a4','a5','a6','a7']
                else:
                    keys = ['a0','a1','a2','a3','a4','a5','a6']
                for key in keys:
                    val = _attr(fit, key, default=None)
                    if val is None:
                        val = _text(_child(fit, key), None)
                    fv = _float_or_none(val)
                    if fv is not None:
                        coeffs.append(fv)
                if len(coeffs) != 7:
                    # fallback for whitespace coefficient lists
                    toks = _text(fit).replace(',', ' ').split()
                    coeffs = []
                    for tok in toks:
                        fv = _float_or_none(tok)
                        if fv is not None:
                            coeffs.append(fv)
                    coeffs = coeffs[:7]
                if len(coeffs) == 7:
                    chapter = _int_or_none(_attr(fit, 'chapter', default=None))
                    label = _attr(fit, 'label', 'note', 'source', default='') or source
                    fits.append(RateFit(coeffs, label=label, chapter=chapter))

        tabular = None
        for tab in _children(el, 'tabular_rate', 'tabular-rate', 'rate_table', 'rate-table'):
            xs, ys = [], []
            for row in list(tab):
                if _tag(row) not in {'point', 'row', 'rate'}:
                    continue
                t9v = _attr(row, 't9', 'T9', default=None)
                rv = _attr(row, 'rate', 'value', default=None)
                toks = _text(row).replace(',', ' ').split()
                if t9v is None and toks:
                    t9v = toks[0]
                if rv is None and len(toks) >= 2:
                    rv = toks[1]
                ft, fr = _float_or_none(t9v), _float_or_none(rv)
                if ft is not None and fr is not None:
                    xs.append(ft)
                    ys.append(fr)
            if xs and ys:
                tabular = TabularRate(xs, ys)

        q = _attr(el, 'q', 'q_value', 'q-value', default=None)
        if q is None:
            q = _text(_child(el, 'q', 'q_value', 'q-value'), '0')
        qv = _float_or_none(q) or 0.0
        label = _attr(el, 'label', 'id', default='') or ''
        rxn = Reaction(reactants, products, rate_fits=fits, tabular_rate=tabular, q_value=qv, source=source, label=label)
        if constant is not None:
            cv = _float_or_none(constant)
            if cv is not None:
                rxn.constant_rate = cv
        reactions.append(rxn)
    net.reactions = ReactionNetwork(reactions)


def _parse_zones(root, net: Network):
    for el in root.iter():
        if _tag(el) != 'zone':
            continue
        label = (_attr(el, 'label1', default=None), _attr(el, 'label2', default=None), _attr(el, 'label3', default=None))
        if all(v is None for v in label):
            label = (_attr(el, 'id', default=str(len(net.zones))), '0', '0')
        label = tuple(str(v if v is not None else '0') for v in label)
        zone = Zone(label=label)

        # optional_properties may wrap many property records.
        for container in [el] + _children(el, 'optional_properties', 'optional-properties', 'properties'):
            for p in list(container):
                tg = _tag(p)
                if tg not in {'property', 'optional_property', 'optional-property'}:
                    continue
                name = _attr(p, 'name', 'tag1', default=None)
                if name:
                    value = _text(p) or _attr(p, 'value', default='')
                    if container is not el or tg in {'optional_property', 'optional-property'}:
                        zone.optional_properties[str(name)] = value
                    else:
                        zone.properties[str(name)] = value

        # Direct abundance records.
        for c in list(el):
            tg = _tag(c)
            if tg in {'abundance', 'species', 'nuc'}:
                name = _attr(c, 'name', 'species', 'nuc', 'id', default=None) or _text(_child(c, 'name'), None)
                val = _attr(c, 'y', 'abundance', 'value', default=None) or _text(_child(c, 'y', 'abundance'), None) or _text(c)
                if name:
                    fv = _float_or_none(val)
                    if fv is not None:
                        zone.set_abundance(name, fv)
            elif tg in {'mass_fractions', 'mass-fractions', 'abundances'}:
                for ael in list(c):
                    name = _attr(ael, 'name', 'species', 'nuc', 'id', default=None) or _text(_child(ael, 'name'), None)
                    y_val = _attr(ael, 'y', 'abundance', 'value', default=None) or _text(_child(ael, 'y', 'abundance'), None)
                    x_val = _attr(ael, 'x', 'mass_fraction', 'mass-fraction', default=None) or _text(_child(ael, 'x', 'mass_fraction', 'mass-fraction'), None)
                    z_val = _attr(ael, 'z', default=None) or _text(_child(ael, 'z'), None)
                    a_val = _attr(ael, 'a', default=None) or _text(_child(ael, 'a'), None)
                    if name is None and z_val is not None and a_val is not None:
                        try:
                            name = species_from_za(int(float(z_val)), int(float(a_val))).name
                        except Exception:
                            name = None
                    if not name:
                        continue
                    if y_val is not None:
                        fv = _float_or_none(y_val)
                    elif x_val is not None:
                        av = _int_or_none(a_val)
                        if av is None:
                            try:
                                av = net.species.get(normalize_species_name(name), Species.parse(name)).a
                            except Exception:
                                av = 1
                        xv = _float_or_none(x_val)
                        fv = None if xv is None else abundance_from_mass_fraction(xv, av)
                    else:
                        fv = _float_or_none(_text(ael))
                    if fv is not None:
                        zone.set_abundance(name, fv)
        if zone.abundances or zone.properties or zone.optional_properties:
            net.add_zone(zone)


def write_xml(network: Network, path: str | Path) -> None:
    from ..species import is_massless
    root = ET.Element('nucnetpy')
    nucs = ET.SubElement(root, 'nuclear_data')
    # A placeholder has no nuclear data, so writing it as a nuclide record would
    # turn its invented mass excess of zero into data.  It is left out, and the
    # reader recreates and records it from the reactions that name it.  Photons
    # and leptons are reaction participants, not nuclides.
    skip = set(network.species_without_nuclear_data())
    for sp in sorted(network.species.values(), key=lambda s: (s.z, s.a, s.name)):
        key = normalize_species_name(sp.name)
        if key in skip or is_massless(key):
            continue
        ne = ET.SubElement(nucs, 'nuclide')
        ET.SubElement(ne, 'z').text = str(sp.z)
        ET.SubElement(ne, 'a').text = str(sp.a)
        if sp.source:
            ET.SubElement(ne, 'source').text = str(sp.source)
        ET.SubElement(ne, 'mass_excess').text = f'{sp.mass_excess:.17g}'
        if sp.spin is not None:
            ET.SubElement(ne, 'spin').text = f'{sp.spin:.17g}'
        if sp.partition:
            pt = ET.SubElement(ne, 'partf_table')
            for t9, val in sorted(sp.partition.items()):
                pe = ET.SubElement(pt, 'point')
                ET.SubElement(pe, 't9').text = f'{float(t9):.17g}'
                ET.SubElement(pe, 'partf').text = f'{float(val):.17g}'
    rs = ET.SubElement(root, 'reaction_data')
    for r in network.reactions.reactions:
        rel = ET.SubElement(rs, 'reaction')
        if r.source:
            ET.SubElement(rel, 'source').text = str(r.source)
        for p in r.reactants:
            for _ in range(p.count):
                ET.SubElement(rel, 'reactant').text = p.species
        for p in r.products:
            for _ in range(p.count):
                ET.SubElement(rel, 'product').text = p.species
        if r.constant_rate is not None:
            ET.SubElement(rel, 'single_rate').text = f'{r.constant_rate:.17g}'
        for fit in r.rate_fits:
            fe = ET.SubElement(rel, 'non_smoker_fit')
            if fit.label:
                fe.set('source', fit.label)
            for i, c in enumerate(fit.coefficients, start=1):
                ET.SubElement(fe, f'a{i}').text = f'{float(c):.17g}'
    zd = ET.SubElement(root, 'zone_data')
    for zone in network.zones:
        ze = ET.SubElement(zd, 'zone', label1=zone.label[0], label2=zone.label[1], label3=zone.label[2])
        if zone.properties or zone.optional_properties:
            props = ET.SubElement(ze, 'optional_properties')
            for k, v in {**zone.properties, **zone.optional_properties}.items():
                pe = ET.SubElement(props, 'property', name=str(k))
                pe.text = str(v)
        abs_el = ET.SubElement(ze, 'abundances')
        for name, y in sorted(zone.abundances.items()):
            ET.SubElement(abs_el, 'abundance', name=name, y=f'{y:.17g}')
    ET.indent(root)
    ET.ElementTree(root).write(path, encoding='utf-8', xml_declaration=True)


write_zone_xml = write_xml


def write_evolution_xml(network: Network, result, path: str | Path, thermo=None,
                        min_abundance: float = 0.0) -> int:
    """Write an evolution as a libnucnet output file, one zone per output time.

    The file has the layout of a NucNet Tools output file: the nuclear data of
    the evolved nuclides, then a zone for each output time of ``result`` (an
    :class:`~nucnetpy.solver.EvolutionResult`).  The zones are labelled by the
    index of the output time, from 0 for the start, and carry the time and,
    when ``thermo`` is given, ``t9`` and ``rho`` as properties.  Each lists the
    mass fractions of the nuclides whose abundance is positive and at least
    ``min_abundance``.  wnutils reads the file, and so do :func:`read_xml` and
    the ``print-output``, ``species-history`` and ``zone-abundances`` commands.

    libnucnet requires a spin for every nuclide in the nuclear data, so a
    nuclide without a known spin, or without nuclear data at all, is left out
    of it; its mass fractions are still written, with its Z and A.  The
    partition functions are written as libnucnet stores them, as
    ``log10_partf``, the inverse of what :func:`read_xml` does.

    The file is written as it is built, so a long history of a large network
    needs little memory.  Returns the number of zones written.
    """
    import math
    from xml.sax.saxutils import escape, quoteattr
    from ..species import is_massless

    def num(value) -> str:
        # The shortest decimal that reads back as the same double.
        return repr(float(value))

    skip = set(network.species_without_nuclear_data())
    nuclides = {}
    for name in result.species:
        key = normalize_species_name(name)
        if is_massless(key):
            continue
        sp = network.species.get(key)
        if sp is None:
            try:
                sp = Species.parse(key)
            except ValueError:  # not a nuclide name: nothing to write
                continue
        nuclides[key] = sp
    with open(path, "w", encoding="utf-8") as out:
        out.write("<?xml version='1.0' encoding='utf-8'?>\n<libnucnet_input>\n")
        out.write("  <nuclear_data>\n")
        for key, sp in sorted(nuclides.items(), key=lambda kv: (kv[1].z, kv[1].a, kv[0])):
            if key in skip or key not in network.species or sp.spin is None:
                continue
            out.write(f"    <nuclide>\n      <z>{sp.z}</z>\n      <a>{sp.a}</a>\n")
            if sp.source:
                out.write(f"      <source>{escape(str(sp.source))}</source>\n")
            out.write(f"      <mass_excess>{num(sp.mass_excess)}</mass_excess>\n"
                      f"      <spin>{num(sp.spin)}</spin>\n")
            points = [(t9, g) for t9, g in sorted(sp.partition.items()) if g > 0]
            if points:
                out.write("      <partf_table>\n")
                for t9, g in points:
                    out.write(f"        <point><t9>{num(t9)}</t9>"
                              f"<log10_partf>{num(math.log10(g))}</log10_partf></point>\n")
                out.write("      </partf_table>\n")
            out.write("    </nuclide>\n")
        out.write("  </nuclear_data>\n  <zone_data>\n")
        names = [normalize_species_name(s) for s in result.species]
        written = 0
        for k, (t, y) in enumerate(zip(result.time, result.y)):
            out.write(f'    <zone label1="{k}" label2="0" label3="0">\n'
                      "      <optional_properties>\n"
                      f'        <property name="time">{num(t)}</property>\n')
            if thermo is not None:
                t9, rho = thermo(t, y)
                out.write(f'        <property name="t9">{num(t9)}</property>\n'
                          f'        <property name="rho">{num(rho)}</property>\n')
            out.write("      </optional_properties>\n      <mass_fractions>\n")
            for name, value in zip(names, y):
                sp = nuclides.get(name)
                if sp is None or not value > 0.0 or value < min_abundance:
                    continue
                out.write(f"        <nuclide name={quoteattr(name)}><z>{sp.z}</z><a>{sp.a}</a>"
                          f"<x>{num(sp.a * value)}</x></nuclide>\n")
            out.write("      </mass_fractions>\n    </zone>\n")
            written += 1
        out.write("  </zone_data>\n</libnucnet_input>\n")
    return written


def read_xml_string(text: str, reverse_partition_functions: bool = True) -> Network:
    root = ET.fromstring(text)
    return read_xml_root(root, source='string',
                         reverse_partition_functions=reverse_partition_functions)
