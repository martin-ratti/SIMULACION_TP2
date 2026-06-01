"""Testeo de los generadores de distribuciones (TP 2.2).

El enunciado (Tabla 1) marca con "si" en la columna Testeo a:
Uniforme, Exponencial, Normal, Binomial, Poisson y Empirica discreta.
Para esas se aplica una prueba de bondad de ajuste:

    - Continuas (Uniforme, Exponencial, Normal): Kolmogorov-Smirnov contra
      la CDF teorica.
    - Discretas (Binomial, Poisson, Empirica): chi-cuadrado de frecuencias
      observadas vs esperadas.

Para todas las distribuciones (incluso las que el enunciado no exige
testear) se calcula ademas la comparacion de media y varianza empirica
contra la teorica, como verificacion minima.

Las distribuciones de chi-cuadrado y normal se aproximan con ``math`` para
no depender de ``scipy`` (igual criterio que el TP 2.1).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Sequence


# ---------------------------------------------------------------------------
# Utilidades de distribucion (sin scipy)
# ---------------------------------------------------------------------------


def _norm_cdf(x: float, mu: float = 0.0, sigma: float = 1.0) -> float:
    """CDF de la normal N(mu, sigma) via funcion error."""
    return 0.5 * (1.0 + math.erf((x - mu) / (sigma * math.sqrt(2.0))))


def _lower_incomplete_gamma(s: float, x: float) -> float:
    """Gamma incompleta inferior regularizada por serie."""
    if x <= 0:
        return 0.0
    term = 1.0 / s
    total = term
    k = 1
    while True:
        term *= x / (s + k)
        total += term
        if term < total * 1e-12 or k > 2000:
            break
        k += 1
    return total * math.exp(-x + s * math.log(x) - math.lgamma(s))


def _chi2_sf(x: float, dof: int) -> float:
    """p-valor de una chi-cuadrado con dof grados de libertad."""
    if x <= 0:
        return 1.0
    return max(0.0, 1.0 - _lower_incomplete_gamma(dof / 2.0, x / 2.0))


def _ks_p_value(d: float, n: int) -> float:
    """Aproximacion asintotica del p-valor de Kolmogorov-Smirnov."""
    t = (math.sqrt(n) + 0.12 + 0.11 / math.sqrt(n)) * d
    if t <= 0:
        return 1.0
    total = 0.0
    for j in range(1, 101):
        total += (-1) ** (j - 1) * math.exp(-2.0 * (j ** 2) * (t ** 2))
    return min(1.0, max(0.0, 2.0 * total))


# ---------------------------------------------------------------------------
# Resultado
# ---------------------------------------------------------------------------


@dataclass
class GoodnessOfFit:
    """Resultado del testeo de una distribucion."""

    dist_name: str
    test_name: str
    statistic: float
    p_value: float
    passed: bool
    mean_emp: float
    mean_theo: float
    var_emp: float
    var_theo: float
    detail: str = ""
    extra: Dict[str, float] = field(default_factory=dict)

    @property
    def verdict(self) -> str:
        return "PASA" if self.passed else "FALLA"


# ---------------------------------------------------------------------------
# Estadisticos basicos
# ---------------------------------------------------------------------------


def _mean(sample: Sequence[float]) -> float:
    return sum(sample) / len(sample)


def _variance(sample: Sequence[float]) -> float:
    m = _mean(sample)
    return sum((x - m) ** 2 for x in sample) / len(sample)


# ---------------------------------------------------------------------------
# Bondad de ajuste para continuas: Kolmogorov-Smirnov
# ---------------------------------------------------------------------------


def ks_test(
    sample: Sequence[float],
    cdf: Callable[[float], float],
    alpha: float = 0.05,
) -> tuple:
    """KS contra una CDF teorica dada. Devuelve (D, p_value, passed)."""
    n = len(sample)
    ordered = sorted(sample)
    d = 0.0
    for i, x in enumerate(ordered):
        fx = cdf(x)
        d = max(d, (i + 1) / n - fx, fx - i / n)
    p = _ks_p_value(d, n)
    return d, p, (p >= alpha)


# ---------------------------------------------------------------------------
# Bondad de ajuste para discretas: chi-cuadrado
# ---------------------------------------------------------------------------


def chi2_discrete_test(
    sample: Sequence[int],
    pmf: Callable[[int], float],
    support: Sequence[int],
    alpha: float = 0.05,
) -> tuple:
    """Chi-cuadrado de frecuencias para una variable discreta.

    Agrupa las categorias con frecuencia esperada < 5 en las colas para
    que el estadistico sea valido. Devuelve (chi2, p_value, passed, dof).
    """
    n = len(sample)
    # Frecuencias observadas por categoria del soporte.
    observed = {k: 0 for k in support}
    extra_low = 0  # valores por debajo del soporte
    extra_high = 0  # valores por encima del soporte
    smin, smax = min(support), max(support)
    for x in sample:
        if x < smin:
            extra_low += 1
        elif x > smax:
            extra_high += 1
        else:
            observed[x] = observed.get(x, 0) + 1

    # Construir celdas (categoria, observado, esperado) agrupando colas.
    cells = []
    for k in support:
        cells.append([observed[k], n * pmf(k)])
    # Sumar masa fuera del soporte a las celdas extremas.
    if extra_low:
        cells[0][0] += extra_low
    if extra_high:
        cells[-1][0] += extra_high

    # Agrupar celdas con esperado < 5 para validez del chi2.
    grouped = []
    acc_o, acc_e = 0.0, 0.0
    for o, e in cells:
        acc_o += o
        acc_e += e
        if acc_e >= 5.0:
            grouped.append((acc_o, acc_e))
            acc_o, acc_e = 0.0, 0.0
    if acc_e > 0:  # ultima celda residual
        if grouped:
            o_last, e_last = grouped[-1]
            grouped[-1] = (o_last + acc_o, e_last + acc_e)
        else:
            grouped.append((acc_o, acc_e))

    chi2 = sum((o - e) ** 2 / e for o, e in grouped if e > 0)
    dof = max(1, len(grouped) - 1)
    p = _chi2_sf(chi2, dof)
    return chi2, p, (p >= alpha), dof


# ---------------------------------------------------------------------------
# CDF / PMF teoricas (para los tests)
# ---------------------------------------------------------------------------


def uniform_cdf(a: float, b: float) -> Callable[[float], float]:
    def cdf(x: float) -> float:
        if x <= a:
            return 0.0
        if x >= b:
            return 1.0
        return (x - a) / (b - a)
    return cdf


def exponential_cdf(mean: float) -> Callable[[float], float]:
    alpha = 1.0 / mean
    def cdf(x: float) -> float:
        return 0.0 if x < 0 else 1.0 - math.exp(-alpha * x)
    return cdf


def normal_cdf(mu: float, sigma: float) -> Callable[[float], float]:
    def cdf(x: float) -> float:
        return _norm_cdf(x, mu, sigma)
    return cdf


def binomial_pmf(n: int, p: float) -> Callable[[int], float]:
    def pmf(k: int) -> float:
        if k < 0 or k > n:
            return 0.0
        return math.comb(n, k) * (p ** k) * ((1 - p) ** (n - k))
    return pmf


def poisson_pmf(lam: float) -> Callable[[int], float]:
    def pmf(k: int) -> float:
        if k < 0:
            return 0.0
        return math.exp(-lam) * (lam ** k) / math.factorial(k)
    return pmf


def empirical_pmf(values: Sequence, probs: Sequence[float]) -> Callable:
    table = dict(zip(values, probs))
    def pmf(k) -> float:
        return table.get(k, 0.0)
    return pmf
