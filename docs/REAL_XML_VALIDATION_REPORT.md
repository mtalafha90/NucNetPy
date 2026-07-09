# NucNetPy real JINA XML validation report
This report was generated from the uploaded real XML files.
## Input database
- **species**: `8136`
- **reactions**: `81758`
- **zones**: `1`
- **missing_species**: `0`
- **invalid_reactions**: `0`
- **load_seconds**: `15.131`

## Initial zone
- **abundance_count**: `47`
- **xsum**: `1.0`
- **ysum**: `0.7814645329200833`
- **Ye**: `0.8559863592831524`
- **T9**: `0.2`
- **rho**: `15000.0`
- **optional_properties**: `{'rho_0': '1.5e4', 'dt': '1.e-15', 'tau': '0.2', 'steps': '2', 't9_0': '0.20', 'munuekT': '-inf', 'tend': '100'}`

## Nova output XML
- **species**: `480`
- **reactions**: `3089`
- **zones**: `344`
- **zone0_abundance_count**: `83`
- **zone0_xsum**: `1.0`
- **zone0_T9**: `0.19999999999999968`
- **zone0_rho**: `14999.999999999925`
- **load_seconds**: `2.637`

## First 10 reaction-rate checks at T9=0.2, rho=1.5e4

| reaction | rate | fits | constant |
|---|---:|---:|---:|
| `gamma + zr85 -> he4 + sr81` | `1.738191e-149` | 1 | `None` |
| `n + p45 -> gamma + p46` | `1.217623e+01` | 1 | `None` |
| `he4 + tb151 -> ho154 + n` | `2.696862e-312` | 1 | `None` |
| `h1 + hg181 -> au178 + he4` | `2.632329e-59` | 1 | `None` |
| `h1 + pb269 -> bi269 + n` | `2.335002e-62` | 1 | `None` |
| `rf310 -> antineutrinoe + db310 + electron` | `2.032420e+00` | 0 | `2.03242` |
| `as106 + gamma -> as105 + n` | `9.858780e+06` | 1 | `None` |
| `gamma + sb123 -> he4 + in119` | `1.935840e-225` | 1 | `None` |
| `tb203 -> antineutrinoe + dy201 + electron + n + n` | `2.293870e+00` | 0 | `2.29387` |
| `hs338 -> antineutrinoe + electron + mt336 + n + n` | `5.221110e+00` | 0 | `5.22111` |
