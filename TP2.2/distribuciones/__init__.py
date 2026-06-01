"""Paquete de generadores de distribuciones de probabilidad (TP 2.2).

Construye nueve distribuciones a partir de uniformes U(0,1) generados por
el GCL del TP 2.1, siguiendo el capitulo 4 de Naylor.

Uso tipico::

    from distribuciones import LCG, Exponential
    rng = LCG(seed=12345)
    expo = Exponential(rng, mean=2.0)
    valores = expo.sample(10000)
"""

from __future__ import annotations

from .rng import LCG
from .distributions import (
    Uniform,
    Exponential,
    Gamma,
    Normal,
    Pascal,
    Binomial,
    Hypergeometric,
    Poisson,
    EmpiricalDiscrete,
)

__all__ = [
    "LCG",
    "Uniform",
    "Exponential",
    "Gamma",
    "Normal",
    "Pascal",
    "Binomial",
    "Hypergeometric",
    "Poisson",
    "EmpiricalDiscrete",
]
