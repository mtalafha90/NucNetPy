# JINA XML workflow in nucnetpy

The original C++ NucNet Tools examples usually begin from JINA/libnucnet XML
files.  The Python port follows the same idea, but uses a Python object model.

## Files

Use your stored XML files as follows:

- `jina_nuclides.xml`: nuclides/species, including Z, A, mass excess, spin, and partition data.
- `jina_reactions.xml`: reactions and ReacLib rate coefficients.
- optional `zones.xml`: initial abundances or output zones.

## Python usage

```python
from nucnetpy import read_jina_xml

net = read_jina_xml("jina_nuclides.xml", "jina_reactions.xml")
print(net.validate())
```

With zones:

```python
net = read_jina_xml(
    "jina_nuclides.xml",
    "jina_reactions.xml",
    zones_xml="initial_zones.xml",
)
```

## Command-line usage

```bash
nucnetpy jina-summary jina_nuclides.xml jina_reactions.xml
nucnetpy jina-combine jina_nuclides.xml jina_reactions.xml combined_network.xml
nucnetpy validate combined_network.xml
nucnetpy rates combined_network.xml --t9 1.0 --rho 1e5
nucnetpy conservation combined_network.xml
```

The combined file keeps the nuclide records exactly as the nuclide file gives
them. A species that the reactions name but the nuclide file leaves out has no
record in it either; the reader recreates it as a placeholder and records it as
having no nuclear data, as it does for the separate files. `nucnetpy validate`
reports those placeholders in one `no_nuclear_data` warning. It does not report
photons and leptons: they take part in reactions but are not nuclides. Combined
files written before version 1.1.1 stored the placeholders as nuclide records
with a mass excess of zero, so regenerate them.

## Important numerical point

To approach numerical identity with the original C++ tools, use the exact same
JINA XML files, do not reorder or filter reactions unless intentionally testing
that operation, and compare rates/flows at fixed `T9`, `rho`, and abundance
values.
