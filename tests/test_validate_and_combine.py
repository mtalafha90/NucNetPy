"""Photons, leptons and placeholder species through validation and file round trips.

Two behaviours are pinned here.  ``nucnetpy validate`` must not report the
photons and leptons a reaction names as missing species: they are reaction
participants, not nuclides.  And a combined file written by
``nucnetpy jina-combine`` must keep the record of which species have no nuclear
data.  Written out as nuclide records with a mass excess of zero, the
placeholders would let an unbound nuclide compete with real ones in an
equilibrium solve on the combined file.
"""
import pytest

from nucnetpy import (Zone, combine_jina_xml, read_jina_xml, read_xml, solve_nse,
                      validate_network)
from nucnetpy.cli import main
from nucnetpy.species import is_massless

NUCLIDES = """<nuclear_data>
  <nuclide><z>0</z><a>1</a><source>ame11</source><mass_excess>8.071317</mass_excess><spin>0.5</spin></nuclide>
  <nuclide><z>1</z><a>1</a><source>ame11</source><mass_excess>7.288971</mass_excess><spin>0.5</spin></nuclide>
  <nuclide><z>2</z><a>4</a><source>ame11</source><mass_excess>2.424916</mass_excess><spin>0</spin></nuclide>
</nuclear_data>"""

# A beta decay with its two leptons, and a capture that names a photon and a
# nuclide, li5, that the nuclide file leaves out.
REACTIONS = """<reaction_data>
  <reaction><source>wc12</source><reactant>n</reactant>
    <product>h1</product><product>electron</product><product>anti-neutrino_e</product>
    <single_rate>1.0e-3</single_rate></reaction>
  <reaction><source>test</source><reactant>h1</reactant><reactant>he4</reactant>
    <product>li5</product><product>gamma</product>
    <single_rate>1.0e-10</single_rate></reaction>
</reaction_data>"""


@pytest.fixture
def database(tmp_path):
    nuclides = tmp_path / "nuclides.xml"
    reactions = tmp_path / "reaction_data.xml"
    nuclides.write_text(NUCLIDES)
    reactions.write_text(REACTIONS)
    return nuclides, reactions


def test_validate_does_not_report_photons_and_leptons_as_missing(database):
    net = read_jina_xml(*database)
    assert not [i for i in validate_network(net) if i.code == "missing_species"]
    # The command-line validator and Network.validate agree.
    assert net.validate()["missing_species"] == []


def test_validate_still_reports_a_nuclide_missing_from_the_species_map(database):
    net = read_jina_xml(*database)
    del net.species["li5"]
    net.metadata.pop("species_without_nuclear_data")
    missing = [i.message for i in validate_network(net) if i.code == "missing_species"]
    assert missing == ["li5: appears in reactions but not in species list"]


def test_validate_reports_species_without_nuclear_data(database):
    net = read_jina_xml(*database)
    issues = [i for i in validate_network(net) if i.code == "no_nuclear_data"]
    assert len(issues) == 1
    assert "1 species has no nuclide record (li5)" in issues[0].message


def test_combined_file_keeps_the_species_without_nuclear_data(database, tmp_path):
    combined = tmp_path / "combined.xml"
    original = combine_jina_xml(*database, combined)
    assert original.species_without_nuclear_data() == ["li5"]

    reread = read_xml(combined)
    assert reread.species_without_nuclear_data() == ["li5"]
    # The placeholder is not written out as a nuclide record with data.
    assert "<z>3</z>" not in combined.read_text()

    # An equilibrium solve therefore leaves it out, as it does for the original
    # files, instead of treating its mass excess of zero as data.
    a = solve_nse(original, t9=3.0, rho=1.0e7, ye=0.5)
    b = solve_nse(reread, t9=3.0, rho=1.0e7, ye=0.5)
    assert a.success and b.success
    assert "li5" not in b.abundances
    for name, y in a.abundances.items():
        assert b.abundances[name] == pytest.approx(y, rel=1e-9)


def test_command_line_validate_on_a_combined_file(database, tmp_path, capsys):
    combined = tmp_path / "combined.xml"
    assert main(["jina-combine", str(database[0]), str(database[1]),
                 str(combined)]) in (None, 0)
    capsys.readouterr()
    main(["validate", str(combined)])
    out = capsys.readouterr().out.splitlines()
    assert out[0] == "issues 1"
    assert out[1].startswith("warning no_nuclear_data: 1 species has no nuclide record (li5)")


def test_zone_species_without_a_nuclide_record_are_recorded(tmp_path, database):
    zones = tmp_path / "zone.xml"
    zones.write_text("""<zone_data><zone label1="0"><mass_fractions>
      <nuclide name="h1"><z>1</z><a>1</a><x>0.7</x></nuclide>
      <nuclide name="li7"><z>3</z><a>7</a><x>0.3</x></nuclide>
    </mass_fractions></zone></zone_data>""")
    net = read_jina_xml(*database, zones_xml=zones)
    assert net.species_without_nuclear_data() == ["li5", "li7"]
    assert net.validate()["missing_species"] == []


def test_abundance_matrix_leaves_out_photons_and_leptons(database):
    net = read_jina_xml(*database)
    net.add_zone(Zone(abundances={"h1": 0.5, "he4": 0.125}))
    names, matrix = net.abundance_matrix()
    assert names and not [n for n in names if is_massless(n)]
    assert matrix.shape == (1, len(names))


def test_hdf5_round_trip_leaves_out_photons_and_leptons(database, tmp_path):
    pytest.importorskip("h5py")
    from nucnetpy.io.hdf5 import read_hdf5, write_hdf5
    net = read_jina_xml(*database)
    net.add_zone(Zone(abundances={"h1": 0.5, "he4": 0.125}))
    path = tmp_path / "net.h5"
    write_hdf5(net, path)
    back = read_hdf5(path)
    assert not [n for n in back.species if is_massless(n)]
    assert back.zones[0].get_abundance("he4") == pytest.approx(0.125)
