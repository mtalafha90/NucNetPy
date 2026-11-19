"""The history of an evolution, saved as a libnucnet output file.

``write_evolution_xml`` writes what NucNet Tools writes: the nuclear data, then
one zone per output time with its time, T9 and rho and the mass fractions.
NucNetPy, wnutils and NucNet Tools must all read it, so the layout is checked
element by element, and the numbers must read back exactly.
"""
import math
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pytest

from nucnetpy import (
    Network,
    Reaction,
    Species,
    Zone,
    constant_thermo,
    evolve_zone,
    read_xml,
    time_grid,
    write_evolution_xml,
)

GOLDEN = str(Path(__file__).parent / "golden" / "golden_network.xml")


@pytest.fixture
def run():
    net = read_xml(GOLDEN)
    thermo = constant_thermo(1.0, 1.0e5)
    result = evolve_zone(net, net.zone(0), time_grid(0.0, 1.0, 5), thermo=thermo)
    assert result.success
    return net, thermo, result


def test_every_output_time_is_a_zone_that_reads_back_exactly(run, tmp_path):
    net, thermo, result = run
    path = tmp_path / "output.xml"
    assert write_evolution_xml(net, result, path, thermo=thermo) == len(result.time)
    zones = read_xml(path).zones
    assert [z.label for z in zones] == [(str(k), "0", "0") for k in range(len(result.time))]
    for zone, t, y in zip(zones, result.time, result.y):
        props = zone.optional_properties
        assert float(props["time"]) == t
        assert float(props["t9"]) == 1.0 and float(props["rho"]) == 1.0e5
        for name, value in zip(result.species, y):
            if value > 0:
                assert zone.get_abundance(name) == pytest.approx(value, rel=1e-15, abs=0.0)
            else:
                assert name not in zone.abundances


def test_the_layout_is_that_of_libnucnet(run, tmp_path):
    net, thermo, result = run
    path = tmp_path / "output.xml"
    write_evolution_xml(net, result, path, thermo=thermo)
    root = ET.parse(path).getroot()
    assert root.tag == "libnucnet_input"
    assert [child.tag for child in root] == ["nuclear_data", "zone_data"]
    nuclide = root.find("zone_data/zone/mass_fractions/nuclide")
    assert nuclide.get("name") and [c.tag for c in nuclide] == ["z", "a", "x"]
    # libnucnet requires a spin for every nuclide of the nuclear data.
    for record in root.findall("nuclear_data/nuclide"):
        assert record.find("spin") is not None and record.find("mass_excess") is not None


def test_mass_fractions_are_written_not_abundances(run, tmp_path):
    net, _, result = run
    path = tmp_path / "output.xml"
    write_evolution_xml(net, result, path)
    last = ET.parse(path).getroot().findall("zone_data/zone")[-1]
    i = result.species.index("c12")
    record = last.find("mass_fractions/nuclide[@name='c12']")
    assert float(record.find("x").text) == pytest.approx(12 * result.y[-1][i], rel=1e-15)


def test_without_thermo_only_the_time_is_written(run, tmp_path):
    net, _, result = run
    path = tmp_path / "output.xml"
    write_evolution_xml(net, result, path)
    assert set(read_xml(path).zones[-1].optional_properties) == {"time"}


def test_small_abundances_can_be_left_out(run, tmp_path):
    net, _, result = run
    path = tmp_path / "output.xml"
    write_evolution_xml(net, result, path, min_abundance=1e-6)
    for zone in read_xml(path).zones:
        assert zone.abundances and min(zone.abundances.values()) >= 1e-6


def test_partition_functions_and_spins_read_back(tmp_path):
    net = Network()
    net.add_species(Species("he4", 2, 4, mass_excess=2.4249, spin=0.0,
                            partition={1.0: 1.0, 10.0: 1.02}))
    net.add_species(Species("c12", 6, 12, mass_excess=0.0, spin=0.0,
                            partition={1.0: 1.0, 10.0: 1.3}))
    net.reactions.add(Reaction.from_names(["he4", "he4", "he4"], ["c12"],
                                          constant_rate=1.0e-3))
    net.add_zone(Zone(abundances={"he4": 0.25}))
    result = evolve_zone(net, net.zone(0), time_grid(0.0, 1.0, 3),
                         thermo=constant_thermo(1.0, 1.0e5))
    path = tmp_path / "output.xml"
    write_evolution_xml(net, result, path)
    back = read_xml(path, reverse_partition_functions=False).species["c12"]
    assert back.spin == 0.0 and back.mass_excess == 0.0
    assert back.partition == pytest.approx({1.0: 1.0, 10.0: 1.3}, rel=1e-15)
    point = ET.parse(path).getroot().find("nuclear_data/nuclide[a='12']/partf_table/point")
    assert float(point.find("log10_partf").text) == pytest.approx(math.log10(1.0), abs=0)


def test_a_nuclide_without_a_spin_keeps_its_mass_fractions(tmp_path):
    net = Network()
    net.add_species(Species("he4", 2, 4, mass_excess=2.4249))
    net.add_species(Species("c12", 6, 12, mass_excess=0.0, spin=0.0))
    net.reactions.add(Reaction.from_names(["he4", "he4", "he4"], ["c12"],
                                          constant_rate=1.0e-3))
    net.add_zone(Zone(abundances={"he4": 0.25}))
    result = evolve_zone(net, net.zone(0), time_grid(0.0, 1.0, 3),
                         thermo=constant_thermo(1.0, 1.0e5))
    path = tmp_path / "output.xml"
    write_evolution_xml(net, result, path)
    root = ET.parse(path).getroot()
    assert [n.find("a").text for n in root.findall("nuclear_data/nuclide")] == ["12"]
    assert root.find("zone_data/zone/mass_fractions/nuclide[@name='he4']") is not None
    assert np.isclose(read_xml(path).zones[-1].get_abundance("he4"), result.y[-1][0],
                      rtol=1e-15, atol=0.0)
