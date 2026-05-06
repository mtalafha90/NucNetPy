from pathlib import Path

from nucnetpy import read_jina_xml


def test_read_separate_jina_xml_files(tmp_path: Path):
    nuclides = tmp_path / "nuclides.xml"
    reactions = tmp_path / "reactions.xml"
    nuclides.write_text(
        """
        <libnucnet__nuc>
          <nuclide><name>h1</name><z>1</z><a>1</a><mass_excess>7.288970</mass_excess><spin>2</spin></nuclide>
          <nuclide><name>he4</name><z>2</z><a>4</a><mass_excess>2.424915</mass_excess>
            <partf_table><point><t9>1.0</t9><partf>1.0</partf></point></partf_table>
          </nuclide>
        </libnucnet__nuc>
        """
    )
    reactions.write_text(
        """
        <libnucnet__reac>
          <reaction>
            <reactant>h1</reactant><reactant>h1</reactant>
            <product>he4</product>
            <q>1.0</q>
            <rate_data>
              <single_rate>
                <a1>0</a1><a2>0</a2><a3>0</a3><a4>0</a4><a5>0</a5><a6>0</a6><a7>0</a7>
              </single_rate>
            </rate_data>
          </reaction>
        </libnucnet__reac>
        """
    )
    net = read_jina_xml(nuclides, reactions)
    assert "h1" in net.species
    assert "he4" in net.species
    assert net.species["he4"].partition[1.0] == 1.0
    assert len(net.reactions.reactions) == 1
    assert abs(net.reactions.reactions[0].rate(1.0) - 1.0) < 1e-12
