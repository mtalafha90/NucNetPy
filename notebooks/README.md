# nucnetpy Jupyter tutorial notebooks

These notebooks demonstrate the pure-Python `nucnetpy` replacement package.

Recommended order:

1. `00_installation_and_first_network.ipynb`
2. `01_species_zones_and_abundances.ipynb`
3. `02_xml_read_write_and_cli.ipynb`
4. `03_reaction_rates_flows_and_conservation.ipynb`
5. `04_one_zone_evolution.ipynb`
6. `05_nse_screening_and_weak_rates.ipynb`
7. `06_validation_and_regression_workflow.ipynb`

From the package root, install with:

```bash
python -m pip install -e .[dev,hdf5]
python -m pip install notebook matplotlib pandas
jupyter notebook notebooks
```

The examples use a small artificial alpha-chain network so they run quickly. Replace the toy XML/network with your real NucNet XML files when you are ready.

## Real JINA/XML validation notebook

`08_validate_real_jina_files.ipynb` is designed for real libnucnet/JINA files. Copy your large XML files into `notebooks/data/` using these names:

```text
nuclides.xml
reaction_data.xml
zone.xml
output_Nova_exp.xml
```

Then run the notebook to check parser counts, mass-fraction normalization, reaction conservation, and sample rates.
