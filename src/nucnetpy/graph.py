"""Graph export helpers replacing the Graphviz-oriented examples."""
from __future__ import annotations
from typing import Mapping, Optional
from .core import Network


def reaction_network_dot(network: Network, min_rate: float = 0.0, t9: float = 1.0, rho: float = 1.0) -> str:
    lines = ["digraph nucnet {", "  rankdir=LR;"]
    for r in network.reactions.reactions:
        rate = r.rate(t9, rho=rho)
        if rate < min_rate:
            continue
        rnode = f"r_{abs(hash(r.string))}"
        lines.append(f'  "{rnode}" [shape=point,label=""];')
        for p in r.reactants:
            lines.append(f'  "{p.species}" -> "{rnode}" [label="{p.count}"];')
        for p in r.products:
            lines.append(f'  "{rnode}" -> "{p.species}" [label="{p.count}"];')
    lines.append("}")
    return "\n".join(lines)
