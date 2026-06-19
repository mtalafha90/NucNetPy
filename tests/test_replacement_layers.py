import numpy as np
from nucnetpy import Network, Species, Zone, Reaction, ReactionNetwork, RateFit, time_grid, evolve_zone
from nucnetpy import read_xml_string, solve_nse, weak_screening_factor, WeakRateTable, validate_network
from nucnetpy.solver import jacobian_sparsity


def small_network():
    net = Network()
    for name in ["h1", "he4"]:
        net.add_species(Species.parse(name))
    net.reactions = ReactionNetwork([
        Reaction.from_names(["h1","h1"], ["he4"], constant_rate=1e-6)
    ])
    return net


def test_xml_edge_cases_partition_optional_and_tabular():
    xml = '''<nucnet><nuclides>
    <nuclide name="h1" z="1" a="1"><partition_function t9="1" value="2"/></nuclide>
    <nuclide name="he4" z="2" a="4"/>
    </nuclides><reactions><reaction><reactants><nuc name="h1" count="2"/></reactants><products><nuc name="he4"/></products><tabular_rate><point t9="1" rate="3"/></tabular_rate></reaction></reactions>
    <zones><zone label1="0"><optional_property name="note">x</optional_property><abundances><abundance name="h1" y="1"/></abundances></zone></zones></nucnet>'''
    net = read_xml_string(xml)
    assert net.species["h1"].partition[1.0] == 2.0
    assert net.reactions.reactions[0].rate(1.0) == 3.0
    assert net.zones[0].optional_properties["note"] == "x"


def test_sparse_jacobian_and_solver():
    net = small_network()
    z = Zone(abundances={"h1":1.0, "he4":0.0}, properties={"t9":"1", "rho":"1"})
    res = evolve_zone(net, z, time_grid(0, 10, 5), method="bdf")
    assert res.success
    assert res.final_abundances["he4"] >= 0
    sp = jacobian_sparsity(net, ["h1", "he4"])
    assert sp.shape == (2,2)


def test_screening_weakrate_validation_nse():
    assert weak_screening_factor(1,1,1.0,1e5,0.5) >= 1.0
    wr = WeakRateTable("h1", "he4", [1.0,2.0], [1e6,1e7], [[1e-3,2e-3],[3e-3,4e-3]])
    assert wr.rate(1.5, 3e6) > 0
    net = small_network()
    issues = validate_network(net)
    assert any(i.code == 'nonconserving_reaction' for i in issues)
    nse = solve_nse(net, 5.0, 1e7, 0.5)
    assert nse.abundances


def _iron_peak_network():
    data = [("n", 8.071), ("h1", 7.289), ("he4", 2.425), ("c12", 0.0),
            ("o16", -4.737), ("si28", -21.49), ("fe56", -60.6), ("ni56", -53.9)]
    net = Network()
    for name, me in data:
        net.add_species(Species.parse(name, mass_excess=me))
    return net


def test_nse_satisfies_constraints():
    net = _iron_peak_network()
    for t9, rho, ye in [(5.0, 1e8, 0.5), (3.0, 1e7, 0.5), (10.0, 1e9, 0.45)]:
        res = solve_nse(net, t9, rho, ye)
        assert res.success
        assert abs(res.xsum - 1.0) < 1e-6
        assert abs(res.computed_ye - ye) < 1e-6


def test_nse_symmetric_matter_favors_n_equals_z():
    # At Ye=0.5 the most-bound N=Z iron-peak nucleus (ni56) should dominate.
    net = _iron_peak_network()
    res = solve_nse(net, 5.0, 1e8, 0.5)
    x = {name: Species.parse(name).a * y for name, y in res.abundances.items()}
    dominant = max(x, key=x.get)
    assert dominant == "ni56"
    assert x["ni56"] > 0.5
