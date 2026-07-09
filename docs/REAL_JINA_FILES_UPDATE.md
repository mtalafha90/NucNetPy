# Real JINA XML support update

This version was tested against the uploaded real files:

- `nuclides.xml`
- `reaction_data.xml`
- `zone(1).xml`
- `output_Nova_exp(1).xml`

## Main fixes

1. Faster parsing of large JINA reaction XML files.
2. Correct parsing of JINA `non_smoker_fit` blocks, including nested `<fit>` records.
3. Correct parsing of `single_rate` beta-decay style records.
4. Correct conversion of zone mass fractions `X` to abundances `Y = X/A`.
5. Support for special particles in reaction sides: `gamma`, `electron`, `positron`, `neutrino_e`, and `anti-neutrino_e`.
6. Optional zone properties such as `t9`, `rho`, `t9_0`, and `rho_0` are now used by `Zone.temperature9()` and `Zone.density()`.
7. Added `validation/validate_real_jina_xml.py`.
8. Added notebook `notebooks/08_validate_real_jina_files.ipynb`.

## Validation result on the uploaded files

```text
Input JINA database:
  species:           8136
  reactions:         81758
  zones:             1
  missing species:   0
  invalid reactions: 0

Initial zone:
  abundances: 47
  sum X:      1.0
  T9:         0.2
  rho:        15000.0

Nova output XML:
  species:   480
  reactions: 3089
  zones:     344
  zone 0 sum X: 1.0
```

## Run the validation script

```bash
python validation/validate_real_jina_xml.py \
  --nuclides /path/to/nuclides.xml \
  --reactions /path/to/reaction_data.xml \
  --zone /path/to/zone.xml \
  --output /path/to/output_Nova_exp.xml
```
