"""Shared definitions for the reproduction of Abdallah and Ershaidat (2024).

The paper is "A Study of the Abundance of Aluminium in Massive Stars",
Jordan Journal of Physics 17 (3) 273-282, doi:10.47011/17.3.2.  It follows the
core of a 25 solar-mass star through six burning stages.  Each stage is held
at a constant temperature and density for a fixed time (its Table 1), and each
one starts from the composition the previous one left behind.

This module holds what the other scripts share: the six stages, the species
each figure of the paper shows, the numbers the paper quotes, the file paths
and the functions that read the inputs.
"""
from __future__ import annotations

import gzip
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
RESULTS = HERE / "results"
FIGURES = HERE / "figures"

NETWORK_FILE = DATA / "network_reaclib_default2_20220329_z28.xml.gz"
COMPOSITION_FILE = DATA / "initial_composition_lodders2003.json"
PROVENANCE_FILE = DATA / "provenance.json"

#: One Julian year in seconds.  The paper does not say which year it uses.
YEAR = 3.15576e7
UNITS = {"Myr": 1.0e6 * YEAR, "kyr": 1.0e3 * YEAR, "yr": YEAR, "d": 86400.0}


@dataclass(frozen=True)
class Stage:
    """One core burning stage of the paper's Table 1."""

    key: str
    title: str
    t9: float
    rho: float
    duration: float          #: in ``unit``
    unit: str
    first_output: float      #: first time shown on the paper's axis, in ``unit``
    figure: int
    species: Tuple[str, ...]  #: the species the paper's figure shows

    @property
    def seconds(self) -> float:
        return self.duration * UNITS[self.unit]


STAGES: Tuple[Stage, ...] = (
    Stage("hydrogen", "Hydrogen burning", 0.0381, 3.81, 6.7, "Myr", 1.0e-18, 1,
          ("h1", "he4", "c12", "c13", "n14", "o16", "o17", "o18")),
    Stage("helium", "Helium burning", 0.196, 7.6e2, 0.839, "Myr", 1.0e-20, 2,
          ("he4", "c12", "o16", "ne20", "mg24", "mg25", "mg26")),
    Stage("carbon", "Carbon burning", 0.9, 1.0e5, 3.2, "kyr", 1.0e-16, 3,
          ("c12", "o16", "ne20", "ne21", "ne22", "na23", "mg24", "mg25", "mg26",
           "al27")),
    Stage("neon", "Neon burning", 1.5, 5.0e6, 0.891, "yr", 1.0e-12, 4,
          ("o16", "ne20", "mg24", "mg25", "mg26", "al27", "si28", "si29", "si30")),
    Stage("oxygen", "Oxygen burning", 2.2, 3.0e6, 0.412, "yr", 1.0e-12, 5,
          ("o16", "al27", "si28", "si29", "si30", "p31", "s32", "s33", "s34",
           "cl35", "ar36")),
    # Table 1 of the paper gives 0.733 yr, but its Figs. 6, 7(d) and 8(d) run
    # to 0.733 d, and 0.733 d is what reproduces the iron-group composition
    # its text quotes.  The figures are what the calculation used.
    Stage("silicon", "Silicon burning", 3.6, 3.0e7, 0.733, "d", 1.0e-10, 6,
          ("mg24", "al27", "si28", "si29", "si30", "s32", "s34", "ar36", "cr52",
           "mn53", "fe54", "fe55", "fe56", "co57", "ni58")),
)
STAGE_BY_KEY = {s.key: s for s in STAGES}
ADVANCED = ("carbon", "neon", "oxygen", "silicon")

#: Final mass fractions quoted in the paper's text and captions.
PAPER_FINAL: Dict[str, Dict[str, float]] = {
    "hydrogen": {"h1": 2.73e-5, "he4": 0.985},
    # The text gives "20Ne (0.0141)", but the paper's Fig. 2 shows 20Ne near
    # 1e-3 throughout, and its Fig. 3 starts 22Ne at about 0.014: the quoted
    # value belongs to 22Ne.  See README.md.
    "helium": {"he4": 0.006, "c12": 0.3075, "o16": 0.6671, "ne22": 0.0141},
    "carbon": {"c12": 0.004, "o16": 0.497, "ne20": 0.409, "na23": 0.023,
               "mg24": 0.040, "mg25": 0.011},
    "neon": {"o16": 0.708, "mg24": 0.115, "al27": 0.015, "si28": 0.112,
             "si29": 0.015, "si30": 0.013},
    # The text gives 0.299 for 32S and the caption of Fig. 5 gives 0.290.
    "oxygen": {"o16": 0.001, "si28": 0.565, "s32": 0.299, "s34": 0.031,
               "ar36": 0.028},
    "silicon": {"cr52": 0.008, "mn53": 0.018, "fe54": 0.683, "fe55": 0.047,
                "fe56": 0.081, "co57": 0.020, "ni58": 0.103},
}

#: Values the paper quotes from Iliadis, Nuclear Physics of Stars (2015).
ILIADIS_FINAL: Dict[str, Dict[str, float]] = {
    "carbon": {"o16": 0.600, "ne20": 0.350},
    "neon": {"o16": 0.770, "mg24": 0.110},
    "oxygen": {"si28": 0.54, "s32": 0.280, "s34": 0.044, "ar36": 0.027},
}

#: Average 26Al/27Al mass-fraction ratio per stage, and over the star's life
#: (the paper's abstract).  The paper does not say how it averaged.
PAPER_RAL_MEAN = {"carbon": 2.85e-4, "neon": 2.23e-4, "oxygen": 1.29e-4,
                  "silicon": 3.46e-3}
PAPER_RAL_OVERALL = 1.68e-4


def load_network(path: Path = NETWORK_FILE, reverse_partition_functions: bool = True):
    """Read the gzipped libnucnet-style network with NucNetPy's XML reader.

    The reader puts the partition-function factor on the rates that ReacLib
    derives by detailed balance; it finds them by ReacLib's ``v`` flag, which
    the network file keeps on each fit label.
    """
    from nucnetpy.io.xml import read_xml
    with gzip.open(path, "rb") as fh:
        return read_xml(fh, reverse_partition_functions=reverse_partition_functions)


def load_composition(path: Path = COMPOSITION_FILE) -> Dict[str, float]:
    """Initial mass fractions, by NucNetPy species name."""
    return json.loads(Path(path).read_text())["mass_fractions"]
