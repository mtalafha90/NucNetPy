# Command line

Installing the package puts a `nucnetpy` command on the path. This page is
generated from the argument parser by `docs/manual/generate_cli_reference.py`,
so it matches the code.

```
usage: nucnetpy [-h]
                {summary,print-output,largest-x,zone-abundances,zone-properties,element-abundances,rates,flows,ydot,net-flows,charge-flows,timescales,entropy-generation,conservation,remove-duplicates,remove-invalid,export-zone-xml,reactions-latex,species-history,evolve-zone,energy-generation,net-dot,validate,nse,qse,jina-summary,jina-combine}
                ...

Pure-Python NucNet Tools replacement commands

positional arguments:
  {summary,print-output,largest-x,zone-abundances,zone-properties,element-abundances,rates,flows,ydot,net-flows,charge-flows,timescales,entropy-generation,conservation,remove-duplicates,remove-invalid,export-zone-xml,reactions-latex,species-history,evolve-zone,energy-generation,net-dot,validate,nse,qse,jina-summary,jina-combine}
    summary             count the species, reactions and zones, and check the
                        network
    print-output        print the properties and abundances of each zone
    largest-x           list the largest mass fractions in a zone
    zone-abundances     list the abundances in a zone
    zone-properties     list the properties of a zone
    element-abundances  list the abundances in a zone summed by element
    rates               list the rate of every reaction at a given T9 and rho
    flows               list the reaction flows in a zone
    ydot                list dY/dt for every species in a zone
    net-flows           list the forward, detailed-balance reverse and net
                        flows in a zone
    charge-flows        list each reaction's contribution to dYe/dt in a zone
    timescales          list the shortest species timescales, Y/|dY/dt|, in a
                        zone
    entropy-generation  print the entropy generation rate of a zone, dS/dt per
                        nucleon in k_B/s
    conservation        list the reactions that do not conserve A and Z
    remove-duplicates   write a copy of the network without duplicate
                        reactions
    remove-invalid      write a copy of the network without the reactions that
                        do not conserve A and Z
    export-zone-xml     write a copy of the network holding only one zone
    reactions-latex     write the reactions and their Q values as a LaTeX
                        table
    species-history     list the abundance of one species in every zone
    evolve-zone         integrate a zone in time at constant T9 and rho, and
                        save its history
    energy-generation   print the nuclear energy generation rate of a zone, in
                        erg/g/s
    net-dot             write the network as a Graphviz DOT graph
    validate            check the network and its zones for problems
    nse                 solve for nuclear statistical equilibrium
    qse                 solve for constrained cluster equilibrium (libnuceq)
    jina-summary        summarise separate JINA nuclide and reaction XML files
    jina-combine        combine JINA nuclide and reaction XML files into one
                        nucnetpy XML file

options:
  -h, --help            show this help message and exit
```

## `summary`

```
usage: nucnetpy summary [-h] xml

count the species, reactions and zones, and check the network

positional arguments:
  xml         network XML file

options:
  -h, --help  show this help message and exit
```

## `print-output`

```
usage: nucnetpy print-output [-h] [--max-zones MAX_ZONES]
                             [--min-abundance MIN_ABUNDANCE]
                             xml

print the properties and abundances of each zone

positional arguments:
  xml                   network XML file

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

list the largest mass fractions in a zone

positional arguments:
  xml                   network XML file

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

list the abundances in a zone

positional arguments:
  xml                   network XML file

options:
  -h, --help            show this help message and exit
  --zone-index ZONE_INDEX
  --min-abundance MIN_ABUNDANCE
```

## `zone-properties`

```
usage: nucnetpy zone-properties [-h] [--zone-index ZONE_INDEX] xml

list the properties of a zone

positional arguments:
  xml                   network XML file

options:
  -h, --help            show this help message and exit
  --zone-index ZONE_INDEX
```

## `element-abundances`

```
usage: nucnetpy element-abundances [-h] [--zone-index ZONE_INDEX] xml element

list the abundances in a zone summed by element

positional arguments:
  xml                   network XML file
  element               element symbol, such as fe; '*' lists every element

options:
  -h, --help            show this help message and exit
  --zone-index ZONE_INDEX
```

## `rates`

```
usage: nucnetpy rates [-h] [--t9 T9] [--rho RHO] [--min-rate MIN_RATE] xml

list the rate of every reaction at a given T9 and rho

positional arguments:
  xml                  network XML file

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

list the reaction flows in a zone

positional arguments:
  xml                   network XML file

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

list dY/dt for every species in a zone

positional arguments:
  xml                   network XML file

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

list the forward, detailed-balance reverse and net flows in a zone

positional arguments:
  xml                   network XML file

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

list each reaction's contribution to dYe/dt in a zone

positional arguments:
  xml                   network XML file

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

list the shortest species timescales, Y/|dY/dt|, in a zone

positional arguments:
  xml                   network XML file

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

print the entropy generation rate of a zone, dS/dt per nucleon in k_B/s

positional arguments:
  xml                   network XML file

options:
  -h, --help            show this help message and exit
  --zone-index ZONE_INDEX
  --t9 T9
  --rho RHO
```

## `conservation`

```
usage: nucnetpy conservation [-h] [--max MAX] xml

list the reactions that do not conserve A and Z

positional arguments:
  xml         network XML file

options:
  -h, --help  show this help message and exit
  --max MAX
```

## `remove-duplicates`

```
usage: nucnetpy remove-duplicates [-h] xml output

write a copy of the network without duplicate reactions

positional arguments:
  xml         network XML file
  output      XML file to write

options:
  -h, --help  show this help message and exit
```

## `remove-invalid`

```
usage: nucnetpy remove-invalid [-h] xml output

write a copy of the network without the reactions that do not conserve A and Z

positional arguments:
  xml         network XML file
  output      XML file to write

options:
  -h, --help  show this help message and exit
```

## `export-zone-xml`

```
usage: nucnetpy export-zone-xml [-h] [--zone-index ZONE_INDEX] xml output

write a copy of the network holding only one zone

positional arguments:
  xml                   network XML file
  output                XML file to write

options:
  -h, --help            show this help message and exit
  --zone-index ZONE_INDEX
```

## `reactions-latex`

```
usage: nucnetpy reactions-latex [-h] xml output

write the reactions and their Q values as a LaTeX table

positional arguments:
  xml         network XML file
  output      LaTeX file to write

options:
  -h, --help  show this help message and exit
```

## `species-history`

```
usage: nucnetpy species-history [-h] xml species

list the abundance of one species in every zone

positional arguments:
  xml         network XML file
  species     species name, such as fe56

options:
  -h, --help  show this help message and exit
```

## `evolve-zone`

```
usage: nucnetpy evolve-zone [-h] [--zone-index ZONE_INDEX] [--t0 T0] [--t1 T1]
                            [--steps STEPS] [--t9 T9] [--rho RHO] [--dt DT]
                            [--method {bdf,radau,lsoda,rk45,dop853,adaptive_implicit_euler,rk4,euler,implicit_euler,backward_euler}]
                            [--log-time] [--min-abundance Y] [--zmax ZMAX]
                            [--amax AMAX] [--rtol RTOL] [--atol ATOL]
                            [--max-change FRACTION] [--min-mass-fraction X]
                            [--progress SECONDS]
                            [--output FOLDER | --no-output]
                            xml

integrate a zone in time at constant T9 and rho, and save its history

positional arguments:
  xml                   network XML file

options:
  -h, --help            show this help message and exit
  --zone-index ZONE_INDEX
  --t0 T0               start time in seconds (default 0)
  --t1 T1               end time in seconds (default: the zone's tend
                        property, else 1)
  --steps STEPS         number of output times (default 50)
  --t9 T9               temperature in units of 10^9 K (default: the zone's t9
                        or t9_0 property)
  --rho RHO             density in g/cm^3 (default: the zone's rho or rho_0
                        property)
  --dt DT               first time step in seconds (default: the zone's dt
                        property; without either, the integrator chooses)
  --method {bdf,radau,lsoda,rk45,dop853,adaptive_implicit_euler,rk4,euler,implicit_euler,backward_euler}
  --log-time
  --min-abundance Y     leave out of the saved file the abundances smaller
                        than Y (default 0: every positive abundance is saved)
  --zmax ZMAX           keep only nuclides with Z <= ZMAX, and the reactions
                        among them
  --amax AMAX           keep only nuclides with A <= AMAX, and the reactions
                        among them
  --rtol RTOL           relative tolerance of the SciPy integrators (default
                        1e-6)
  --atol ATOL           absolute tolerance on the abundances (default 1e-30,
                        which makes the relative tolerance govern every
                        abundance; 1e-12 is often enough and much faster)
  --max-change FRACTION
                        adaptive_implicit_euler only: the largest fractional
                        change of an abundance in one step (default 0.005;
                        0.05 suits explosive burning)
  --min-mass-fraction X
                        adaptive_implicit_euler only: abundances with a
                        smaller mass fraction do not limit the step (default
                        1e-15)
  --progress SECONDS    report progress on standard error every SECONDS of
                        computing time (default 30; 0 turns it off)
  --output FOLDER       folder to save the results in, as output.xml: one zone
                        per output time, in the libnucnet format of NucNet
                        Tools (default: a new folder evolve-zone_DATE_TIME in
                        the current directory)
  --no-output           save nothing; only print the summary
```

## `energy-generation`

```
usage: nucnetpy energy-generation [-h] [--zone-index ZONE_INDEX] [--t9 T9]
                                  [--rho RHO]
                                  xml

print the nuclear energy generation rate of a zone, in erg/g/s

positional arguments:
  xml                   network XML file

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

write the network as a Graphviz DOT graph

positional arguments:
  xml                   network XML file

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

check the network and its zones for problems

positional arguments:
  xml                   network XML file

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

solve for nuclear statistical equilibrium

positional arguments:
  xml                   network XML file

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

solve for constrained cluster equilibrium (libnuceq)

positional arguments:
  xml                   network XML file

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

summarise separate JINA nuclide and reaction XML files

positional arguments:
  nuclides_xml          nuclide XML file
  reactions_xml         reaction XML file

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

combine JINA nuclide and reaction XML files into one nucnetpy XML file

positional arguments:
  nuclides_xml          nuclide XML file
  reactions_xml         reaction XML file
  output_xml            combined XML file to write

options:
  -h, --help            show this help message and exit
  --zones-xml ZONES_XML
```
