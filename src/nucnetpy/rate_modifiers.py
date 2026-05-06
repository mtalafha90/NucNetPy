"""Pure-Python rate modification hooks replacing user/rate_modifiers.*."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Callable, Dict, Mapping, Optional

RateModifier = Callable[[object, float, float, Optional[float], Mapping[str, float]], float]

@dataclass
class RateModifierRegistry:
    modifiers: Dict[str, RateModifier] = field(default_factory=dict)

    def register(self, key: str, func: RateModifier) -> None:
        self.modifiers[str(key)] = func

    def factor(self, reaction, t9: float, rho: float, ye: Optional[float], abundances: Mapping[str, float]) -> float:
        fac = 1.0
        for key in getattr(reaction, 'metadata', {}) or {}:
            if key in self.modifiers:
                fac *= float(self.modifiers[key](reaction, t9, rho, ye, abundances))
        if getattr(reaction, 'label', '') in self.modifiers:
            fac *= float(self.modifiers[reaction.label](reaction, t9, rho, ye, abundances))
        return float(fac)


def constant_factor(value: float) -> RateModifier:
    return lambda reaction, t9, rho, ye, abundances: float(value)


def exp_temperature_factor(alpha: float, power: float = 1.0) -> RateModifier:
    import math
    return lambda reaction, t9, rho, ye, abundances: math.exp(float(alpha) * (float(t9) ** float(power)))
