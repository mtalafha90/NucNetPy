"""Validate NucNetPy against real JINA/libnucnet XML files.

Usage
-----
python validation/validate_real_jina_xml.py \
    --nuclides /path/to/nuclides.xml \
    --reactions /path/to/reaction_data.xml \
    --zone /path/to/zone.xml \
    --output /path/to/output_Nova_exp.xml
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from nucnetpy.io.jina import read_jina_xml
from nucnetpy.io.xml import read_xml


def summarize(nuclides, reactions, zone=None, output=None, sample_t9=0.2, sample_rho=1.5e4):
    t0 = time.time()
    net = read_jina_xml(nuclides, reactions, zone)
    dt_net = time.time() - t0
    val = net.validate()
    report = {
        "input_database": {
            "species": len(net.species),
            "reactions": len(net.reactions.reactions),
            "zones": len(net.zones),
            "missing_species": len(val["missing_species"]),
            "invalid_reactions": len(val["invalid_reactions"]),
            "load_seconds": dt_net,
        },
        "sample_rates": [],
    }
    if net.zones:
        z = net.zones[0]
        report["initial_zone"] = {
            "abundance_count": len(z.abundances),
            "xsum": z.xsum(net.species),
            "ysum": z.ysum(),
            "Ye": z.ye(net.species),
            "T9": z.temperature9(),
            "rho": z.density(),
            "properties": z.properties,
            "optional_properties": z.optional_properties,
        }
    for r in net.reactions.reactions[:10]:
        report["sample_rates"].append({
            "reaction": r.string,
            "rate": r.rate(sample_t9, sample_rho),
            "number_of_fits": len(r.rate_fits),
            "constant_rate": r.constant_rate,
        })
    if output:
        t0 = time.time()
        out = read_xml(output)
        dt_out = time.time() - t0
        report["output_xml"] = {
            "species": len(out.species),
            "reactions": len(out.reactions.reactions),
            "zones": len(out.zones),
            "load_seconds": dt_out,
        }
        if out.zones:
            z0 = out.zones[0]
            zf = out.zones[-1]
            report["output_xml"].update({
                "zone0_abundance_count": len(z0.abundances),
                "zone0_xsum": z0.xsum(out.species),
                "zone0_T9": z0.temperature9(),
                "zone0_rho": z0.density(),
                "final_zone_abundance_count": len(zf.abundances),
                "final_zone_xsum": zf.xsum(out.species),
                "final_zone_T9": zf.temperature9(),
                "final_zone_rho": zf.density(),
            })
    return report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nuclides", required=True)
    ap.add_argument("--reactions", required=True)
    ap.add_argument("--zone")
    ap.add_argument("--output")
    ap.add_argument("--json", default="real_xml_validation_summary.json")
    args = ap.parse_args()
    report = summarize(args.nuclides, args.reactions, args.zone, args.output)
    Path(args.json).write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
