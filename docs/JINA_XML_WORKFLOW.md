# JINA XML workflow

NucNet Tools calculations usually begin from JINA/libnucnet XML files, and so
do NucNetPy's. The files are read into NucNetPy's Python objects: a `Network`
holding the species, the reactions and the zones.

## Files

- `nuclides.xml`: the nuclides, with Z, A, mass excess, spin and partition
  function.
- `reaction_data.xml`: the reactions and their ReacLib rate coefficients.
- `zone.xml` (optional): initial abundances, or the zones of an earlier
  calculation.

## Python usage

```python
from nucnetpy import read_jina_xml

net = read_jina_xml("nuclides.xml", "reaction_data.xml")
print(net.validate())
```

With zones:

```python
net = read_jina_xml("nuclides.xml", "reaction_data.xml", zones_xml="zone.xml")
```

## Command-line usage

```bash
nucnetpy jina-summary nuclides.xml reaction_data.xml
nucnetpy jina-combine nuclides.xml reaction_data.xml combined_network.xml
nucnetpy validate combined_network.xml
nucnetpy rates combined_network.xml --t9 1.0 --rho 1e5
nucnetpy conservation combined_network.xml
```

The combined file keeps the nuclide records exactly as the nuclide file gives
them. A species that the reactions name but the nuclide file leaves out has no
record in it either; the reader recreates it as a placeholder and records it as
having no nuclear data, as it does for the separate files. `nucnetpy validate`
reports those placeholders in one `no_nuclear_data` warning. It does not report
photons and leptons: they take part in reactions but are not nuclides.

## Comparing with NucNet Tools

To compare with the C++ tools, use exactly the same XML files, do not reorder
or filter the reactions unless that operation is what you are testing, and
compare rates and flows at fixed `T9`, `rho` and abundances.
`docs/EXACT_NUMERICAL_IDENTITY.md` describes the comparison in full.
