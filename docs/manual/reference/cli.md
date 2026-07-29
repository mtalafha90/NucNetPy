# Command line

Installing the package puts a `nucnetpy` command on the path. This page is
generated from the argument parser, so it matches the code.

```
usage: nucnetpy [-h]
                {summary,print-output,largest-x,zone-abundances,zone-properties,element-abundances,rates,flows,ydot,net-flows,charge-flows,timescales,entropy-generation,conservation,remove-duplicates,remove-invalid,export-zone-xml,reactions-latex,species-history,evolve-zone,energy-generation,net-dot,validate,nse,qse,jina-summary,jina-combine}
                ...

Pure-Python NucNet Tools replacement commands

positional arguments:
  {summary,print-output,largest-x,zone-abundances,zone-properties,element-abundances,rates,flows,ydot,net-flows,charge-flows,timescales,entropy-generation,conservation,remove-duplicates,remove-invalid,export-zone-xml,reactions-latex,species-history,evolve-zone,energy-generation,net-dot,validate,nse,qse,jina-summary,jina-combine}
    net-flows           forward, detailed-balance reverse, and net fluxes
    charge-flows        per-reaction dYe/dt contributions
    timescales          shortest species timescales Y/|dY/dt|
    entropy-generation  dS/dt per nucleon in k_B/s
    qse                 constrained cluster equilibrium (libnuceq)
    jina-summary        summarize separate JINA nuclide and reaction XML files
    jina-combine        combine JINA nuclide and reaction XML files into one
                        nucnetpy XML file

options:
  -h, --help            show this help message and exit
```

## `summary`

```
usage: nucnetpy summary [-h] xml

positional arguments:
  xml

options:
  -h, --help  show this help message and exit
```

## `print-output`

```
usage: nucnetpy print-output [-h] [--max-zones MAX_ZONES]
                             [--min-abundance MIN_ABUNDANCE]
                             xml

positional arguments:
  xml

options:
  -h, --help            show this help message and exit
  --max-zones MAX_ZONES
  --min-abundance MIN_ABUNDANCE
```

## `largest-x`

```
usage: nucnetpy largest-x [-h] [-n N] [--zone-index ZONE_INDEX]
                          [--min-x MIN_X]
                          xml

positional arguments:
  xml

options:
  -h, --help            show this help message and exit
  -n N
  --zone-index ZONE_INDEX
  --min-x MIN_X
```

## `zone-abundances`

```
usage: nucnetpy zone-abundances [-h] [--zone-index ZONE_INDEX]
                                [--min-abundance MIN_ABUNDANCE]
                                xml

positional arguments:
  xml

options:
  -h, --help            show this help message and exit
  --zone-index ZONE_INDEX
  --min-abundance MIN_ABUNDANCE
```

## `zone-properties`

```
usage: nucnetpy zone-properties [-h] [--zone-index ZONE_INDEX] xml

positional arguments:
  xml

options:
  -h, --help            show this help message and exit
  --zone-index ZONE_INDEX
```

## `element-abundances`

```
usage: nucnetpy element-abundances [-h] [--zone-index ZONE_INDEX] xml element

positional arguments:
  xml
  element

options:
  -h, --help            show this help message and exit
  --zone-index ZONE_INDEX
```

## `rates`

```
usage: nucnetpy rates [-h] [--t9 T9] [--rho RHO] [--min-rate MIN_RATE] xml

positional arguments:
  xml

options:
  -h, --help           show this help message and exit
  --t9 T9
  --rho RHO
  --min-rate MIN_RATE
```

## `flows`

```
usage: nucnetpy flows [-h] [--zone-index ZONE_INDEX] [--t9 T9] [--rho RHO]
                      [--min-flow MIN_FLOW]
                      xml

positional arguments:
  xml

options:
  -h, --help            show this help message and exit
  --zone-index ZONE_INDEX
  --t9 T9
  --rho RHO
  --min-flow MIN_FLOW
```

## `ydot`

```
usage: nucnetpy ydot [-h] [--zone-index ZONE_INDEX] [--t9 T9] [--rho RHO]
                     [--min-abs MIN_ABS]
                     xml

positional arguments:
  xml

options:
  -h, --help            show this help message and exit
  --zone-index ZONE_INDEX
  --t9 T9
  --rho RHO
  --min-abs MIN_ABS
```

## `net-flows`

```
usage: nucnetpy net-flows [-h] [--zone-index ZONE_INDEX] [--t9 T9] [--rho RHO]
                          [--min-flow MIN_FLOW]
                          xml

positional arguments:
  xml

options:
  -h, --help            show this help message and exit
  --zone-index ZONE_INDEX
  --t9 T9
  --rho RHO
  --min-flow MIN_FLOW
```

## `charge-flows`

```
usage: nucnetpy charge-flows [-h] [--zone-index ZONE_INDEX] [--t9 T9]
                             [--rho RHO] [--min-flow MIN_FLOW]
                             xml

positional arguments:
  xml

options:
  -h, --help            show this help message and exit
  --zone-index ZONE_INDEX
  --t9 T9
  --rho RHO
  --min-flow MIN_FLOW
```

## `timescales`

```
usage: nucnetpy timescales [-h] [--zone-index ZONE_INDEX] [--t9 T9]
                           [--rho RHO] [-n N]
                           xml

positional arguments:
  xml

options:
  -h, --help            show this help message and exit
  --zone-index ZONE_INDEX
  --t9 T9
  --rho RHO
  -n N
```

## `entropy-generation`

```
usage: nucnetpy entropy-generation [-h] [--zone-index ZONE_INDEX] [--t9 T9]
                                   [--rho RHO]
                                   xml

positional arguments:
  xml

options:
  -h, --help            show this help message and exit
  --zone-index ZONE_INDEX
  --t9 T9
  --rho RHO
```

## `conservation`

```
usage: nucnetpy conservation [-h] [--max MAX] xml

positional arguments:
  xml

options:
  -h, --help  show this help message and exit
  --max MAX
```

## `remove-duplicates`

```
usage: nucnetpy remove-duplicates [-h] xml output

positional arguments:
  xml
  output

options:
  -h, --help  show this help message and exit
```

## `remove-invalid`

```
usage: nucnetpy remove-invalid [-h] xml output

positional arguments:
  xml
  output

options:
  -h, --help  show this help message and exit
```

## `export-zone-xml`

```
usage: nucnetpy export-zone-xml [-h] [--zone-index ZONE_INDEX] xml output

positional arguments:
  xml
  output

options:
  -h, --help            show this help message and exit
  --zone-index ZONE_INDEX
```

## `reactions-latex`

```
usage: nucnetpy reactions-latex [-h] xml output

positional arguments:
  xml
  output

options:
  -h, --help  show this help message and exit
```

## `species-history`

```
usage: nucnetpy species-history [-h] xml species

positional arguments:
  xml
  species

options:
  -h, --help  show this help message and exit
```

## `evolve-zone`

```
usage: nucnetpy evolve-zone [-h] [--zone-index ZONE_INDEX] [--t0 T0] [--t1 T1]
                            [--steps STEPS] [--t9 T9] [--rho RHO]
                            [--method METHOD] [--log-time]
                            [--min-abundance MIN_ABUNDANCE]
                            xml

positional arguments:
  xml

options:
  -h, --help            show this help message and exit
  --zone-index ZONE_INDEX
  --t0 T0
  --t1 T1
  --steps STEPS
  --t9 T9
  --rho RHO
  --method METHOD
  --log-time
  --min-abundance MIN_ABUNDANCE
```

## `energy-generation`

```
usage: nucnetpy energy-generation [-h] [--zone-index ZONE_INDEX] [--t9 T9]
                                  [--rho RHO]
                                  xml

positional arguments:
  xml

options:
  -h, --help            show this help message and exit
  --zone-index ZONE_INDEX
  --t9 T9
  --rho RHO
```

## `net-dot`

```
usage: nucnetpy net-dot [-h] [-o OUTPUT] [--t9 T9] [--rho RHO]
                        [--min-rate MIN_RATE]
                        xml

positional arguments:
  xml

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
  --t9 T9
  --rho RHO
  --min-rate MIN_RATE
```

## `validate`

```
usage: nucnetpy validate [-h] [--strict] [--max MAX] [--max-zones MAX_ZONES]
                         xml

positional arguments:
  xml

options:
  -h, --help            show this help message and exit
  --strict
  --max MAX
  --max-zones MAX_ZONES
```

## `nse`

```
usage: nucnetpy nse [-h] --t9 T9 --rho RHO --ye YE [--coulomb]
                    [--min-abundance MIN_ABUNDANCE] [--min-x MIN_X]
                    xml

positional arguments:
  xml

options:
  -h, --help            show this help message and exit
  --t9 T9
  --rho RHO
  --ye YE
  --coulomb             apply Bravo & Garcia-Senz Coulomb corrections
  --min-abundance MIN_ABUNDANCE
  --min-x MIN_X
```

## `qse`

```
usage: nucnetpy qse [-h] --t9 T9 --rho RHO --ye YE [--cluster SP1,SP2,...:Y]
                    [--coulomb] [--min-abundance MIN_ABUNDANCE]
                    [--min-x MIN_X]
                    xml

positional arguments:
  xml

options:
  -h, --help            show this help message and exit
  --t9 T9
  --rho RHO
  --ye YE
  --cluster SP1,SP2,...:Y
                        cluster species and constrained total abundance;
                        repeatable
  --coulomb
  --min-abundance MIN_ABUNDANCE
  --min-x MIN_X
```

## `jina-summary`

```
usage: nucnetpy jina-summary [-h] [--zones-xml ZONES_XML] [--show-invalid]
                             [--max-invalid MAX_INVALID]
                             nuclides_xml reactions_xml

positional arguments:
  nuclides_xml
  reactions_xml

options:
  -h, --help            show this help message and exit
  --zones-xml ZONES_XML
  --show-invalid
  --max-invalid MAX_INVALID
```

## `jina-combine`

```
usage: nucnetpy jina-combine [-h] [--zones-xml ZONES_XML]
                             nuclides_xml reactions_xml output_xml

positional arguments:
  nuclides_xml
  reactions_xml
  output_xml

options:
  -h, --help            show this help message and exit
  --zones-xml ZONES_XML
```

