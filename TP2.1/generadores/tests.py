"""Tests estadisticos para evaluar generadores pseudoaleatorios (TP 2.1).

Se implementan cinco pruebas, todas operando sobre una muestra de
flotantes ``u_i in [0, 1)``:

    1. Chi-cuadrado de uniformidad : reparte la muestra en k intervalos
       y compara las frecuencias observadas con las esperadas.
    2. Test de rachas (runs)       : cuenta las rachas ascendentes y
       descendentes y las contrasta con su distribucion asintotica normal.
    3. Kolmogorov-Smirnov          : maxima distancia entre la CDF empirica
       y la CDF uniforme teorica.
    4. Autocorrelacion (lag k)     : correlacion serial de la muestra
       consigo misma desplazada k posiciones.
    5. Test de poker               : agrupa digitos y clasifica las "manos"
       (todos distintos, un par, etc.) contra sus probabilidades teoricas.

Cada test devuelve un ``TestResult`` con el estadistico, el valor critico
o p-valor, y un veredicto PASA/FALLA a un nivel de significacion ``alpha``.

Las funciones de distribucion (chi2 y normal) se calculan con
aproximaciones de ``math`` para no depender de ``scipy``; son suficientes
para los tamanos de muestra y niveles de significacion habituales.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Sequence

# ---------------------------------------------------------------------------
# Utilidades de distribucion (sin scipy)
# ---------------------------------------------------------------------------


def _norm_cdf(x: float) -> float:
    """CDF de la normal estandar via funcion error."""
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def _norm_sf(x: float) -> float:
    """Funcion de supervivencia (cola superior) de la normal estandar."""
    return 1.0 - _norm_cdf(x)


def _lower_incomplete_gamma(s: float, x: float) -> float:
    """Gamma incompleta inferior gamma(s, x) por serie (x < s + 1)."""
    if x < 0:
        raise ValueError("x debe ser no negativo")
    if x == 0:
        return 0.0
    term = 1.0 / s
    total = term
    n = 1
    while True:
        term *= x / (s + n)
        total += term
        if term < total * 1e-12 or n > 1000:
            break
        n += 1
    return total * math.exp(-x + s * math.log(x) - math.lgamma(s))


def _chi2_cdf(x: float, k: int) -> float:
    """CDF de una chi-cuadrado con k grados de libertad."""
    if x <= 0:
        return 0.0
    return _lower_incomplete_gamma(k / 2.0, x / 2.0)


def _chi2_sf(x: float, k: int) -> float:
    """p-valor: probabilidad de superar x en una chi2 con k g.l."""
    return max(0.0, 1.0 - _chi2_cdf(x, k))


# ---------------------------------------------------------------------------
# Resultado de un test
# ---------------------------------------------------------------------------


@dataclass
class TestResult:
    """Resultado de aplicar un test a una muestra."""

    name: str
    statistic: float
    p_value: float
    passed: bool
    detail: str = ""
    extra: Dict[str, float] = field(default_factory=dict)

    @property
    def verdict(self) -> str:
        return "PASA" if self.passed else "FALLA"


# ---------------------------------------------------------------------------
# 1. Chi-cuadrado de uniformidad
# ---------------------------------------------------------------------------


def chi_square_test(
    sample: Sequence[float], bins: int = 10, alpha: float = 0.05
) -> TestResult:
    """Test chi-cuadrado de bondad de ajuste a la uniforme [0, 1).

    Divide [0, 1) en ``bins`` intervalos de igual ancho y compara las
    frecuencias observadas O_i con la esperada E = n / bins:

        chi2 = sum_i (O_i - E)^2 / E

    Se rechaza la hipotesis de uniformidad si el p-valor < alpha.
    """
    n = len(sample)
    if n == 0:
        raise ValueError("la muestra no puede estar vacia")

    observed = [0] * bins
    for u in sample:
        idx = int(u * bins)
        if idx == bins:  # u == 1.0 por redondeo
            idx = bins - 1
        observed[idx] += 1

    expected = n / bins
    chi2 = sum((o - expected) ** 2 / expected for o in observed)
    dof = bins - 1
    p_value = _chi2_sf(chi2, dof)
    passed = p_value >= alpha

    return TestResult(
        name="Chi-cuadrado (uniformidad)",
        statistic=chi2,
        p_value=p_value,
        passed=passed,
        detail=f"bins={bins}, g.l.={dof}, E={expected:.1f}",
        extra={"dof": dof, "bins": bins},
    )


# ---------------------------------------------------------------------------
# 2. Test de rachas (runs up/down)
# ---------------------------------------------------------------------------


def runs_test(sample: Sequence[float], alpha: float = 0.05) -> TestResult:
    """Test de rachas ascendentes y descendentes.

    Se construye la secuencia de signos de las diferencias consecutivas y
    se cuenta el numero de rachas ``R``. Bajo independencia, para una
    muestra de tamano n:

        E[R]   = (2n - 1) / 3
        Var[R] = (16n - 29) / 90

    El estadistico Z = (R - E[R]) / sqrt(Var[R]) es aproximadamente
    normal estandar. Se rechaza si |Z| supera el valor critico bilateral.
    """
    n = len(sample)
    if n < 2:
        raise ValueError("se necesitan al menos 2 valores")

    runs = 1
    for i in range(1, n - 1):
        prev_up = sample[i] >= sample[i - 1]
        next_up = sample[i + 1] >= sample[i]
        if prev_up != next_up:
            runs += 1

    exp_runs = (2 * n - 1) / 3.0
    var_runs = (16 * n - 29) / 90.0
    z = (runs - exp_runs) / math.sqrt(var_runs)
    p_value = 2.0 * _norm_sf(abs(z))
    passed = p_value >= alpha

    return TestResult(
        name="Rachas (runs up/down)",
        statistic=z,
        p_value=p_value,
        passed=passed,
        detail=f"R={runs}, E[R]={exp_runs:.1f}, sd={math.sqrt(var_runs):.2f}",
        extra={"runs": runs, "expected_runs": exp_runs},
    )


# ---------------------------------------------------------------------------
# 3. Kolmogorov-Smirnov
# ---------------------------------------------------------------------------


def _ks_p_value(d: float, n: int) -> float:
    """Aproximacion del p-valor de KS (formula de Kolmogorov asintotica)."""
    t = (math.sqrt(n) + 0.12 + 0.11 / math.sqrt(n)) * d
    if t <= 0:
        return 1.0
    total = 0.0
    for j in range(1, 101):
        total += (-1) ** (j - 1) * math.exp(-2.0 * (j ** 2) * (t ** 2))
    p = 2.0 * total
    return min(1.0, max(0.0, p))


def kolmogorov_smirnov_test(
    sample: Sequence[float], alpha: float = 0.05
) -> TestResult:
    """Test de Kolmogorov-Smirnov contra la uniforme [0, 1).

    Calcula la maxima discrepancia entre la CDF empirica y la teorica
    F(x) = x:

        D = max_i  max( i/n - u_(i),  u_(i) - (i-1)/n )

    donde u_(i) son los valores ordenados. Se rechaza si el p-valor < alpha.
    """
    n = len(sample)
    if n == 0:
        raise ValueError("la muestra no puede estar vacia")

    ordered = sorted(sample)
    d_plus = max((i + 1) / n - ordered[i] for i in range(n))
    d_minus = max(ordered[i] - i / n for i in range(n))
    d = max(d_plus, d_minus)
    p_value = _ks_p_value(d, n)
    passed = p_value >= alpha

    return TestResult(
        name="Kolmogorov-Smirnov",
        statistic=d,
        p_value=p_value,
        passed=passed,
        detail=f"D+={d_plus:.4f}, D-={d_minus:.4f}",
        extra={"d_plus": d_plus, "d_minus": d_minus},
    )


# ---------------------------------------------------------------------------
# 4. Autocorrelacion serial
# ---------------------------------------------------------------------------


def autocorrelation_test(
    sample: Sequence[float], lag: int = 1, alpha: float = 0.05
) -> TestResult:
    """Test de autocorrelacion serial al desplazamiento ``lag``.

    Estima la correlacion entre la muestra y su version desplazada:

        rho_k = ( (1/(n-k)) sum u_i u_{i+k} - 0.25 ) / (1/12)

    Bajo independencia, rho_k ~ N(0, 1/(n-k)) aproximadamente, por lo que
    Z = rho_k * sqrt(n - k) es normal estandar. Se rechaza si |Z| es
    demasiado grande.
    """
    n = len(sample)
    if n <= lag:
        raise ValueError("la muestra es mas corta que el lag pedido")

    m = n - lag
    prod = sum(sample[i] * sample[i + lag] for i in range(m)) / m
    # Bajo H0 (uniformes independientes): E[u_i u_{i+k}] = 1/4, Var(u) = 1/12.
    rho = (prod - 0.25) / (1.0 / 12.0)
    z = rho * math.sqrt(m)
    p_value = 2.0 * _norm_sf(abs(z))
    passed = p_value >= alpha

    return TestResult(
        name=f"Autocorrelacion (lag={lag})",
        statistic=z,
        p_value=p_value,
        passed=passed,
        detail=f"rho={rho:.4f}, pares={m}",
        extra={"rho": rho, "lag": lag},
    )


# ---------------------------------------------------------------------------
# 5. Test de poker
# ---------------------------------------------------------------------------


def poker_test(
    sample: Sequence[float], digits: int = 5, alpha: float = 0.05
) -> TestResult:
    """Test de poker sobre grupos de ``digits`` digitos decimales.

    De cada flotante se toman ``digits`` digitos (por defecto 5) y se
    clasifica la "mano" segun el patron de repeticiones. Se comparan las
    frecuencias observadas de cada categoria con las teoricas y se aplica
    una chi-cuadrado.

    Para el caso clasico de 5 digitos las categorias y probabilidades son:

        Todos distintos   0.30240
        Un par            0.50400
        Dos pares         0.10800
        Terna             0.07200
        Full              0.00900
        Poker             0.00450
        Quintilla         0.00010
    """
    if digits != 5:
        raise ValueError("este test de poker esta calibrado para 5 digitos")

    categories = [
        "todos_distintos",
        "un_par",
        "dos_pares",
        "terna",
        "full",
        "poker",
        "quintilla",
    ]
    probs = {
        "todos_distintos": 0.30240,
        "un_par": 0.50400,
        "dos_pares": 0.10800,
        "terna": 0.07200,
        "full": 0.00900,
        "poker": 0.00450,
        "quintilla": 0.00010,
    }

    def classify(u: float) -> str:
        ds = [int(d) for d in f"{u:.5f}"[2:2 + digits]]
        counts = sorted(_digit_counts(ds).values(), reverse=True)
        if counts == [1, 1, 1, 1, 1]:
            return "todos_distintos"
        if counts == [2, 1, 1, 1]:
            return "un_par"
        if counts == [2, 2, 1]:
            return "dos_pares"
        if counts == [3, 1, 1]:
            return "terna"
        if counts == [3, 2]:
            return "full"
        if counts == [4, 1]:
            return "poker"
        return "quintilla"

    observed = {c: 0 for c in categories}
    n = len(sample)
    for u in sample:
        observed[classify(u)] += 1

    chi2 = 0.0
    for c in categories:
        exp = n * probs[c]
        if exp > 0:
            chi2 += (observed[c] - exp) ** 2 / exp
    dof = len(categories) - 1
    p_value = _chi2_sf(chi2, dof)
    passed = p_value >= alpha

    detail = ", ".join(f"{c}={observed[c]}" for c in categories if observed[c])
    return TestResult(
        name="Poker (5 digitos)",
        statistic=chi2,
        p_value=p_value,
        passed=passed,
        detail=detail,
        extra={"dof": dof},
    )


def _digit_counts(digits: List[int]) -> Dict[int, int]:
    counts: Dict[int, int] = {}
    for d in digits:
        counts[d] = counts.get(d, 0) + 1
    return counts


# ---------------------------------------------------------------------------
# Orquestador de todos los tests
# ---------------------------------------------------------------------------


def run_all_tests(
    sample: Sequence[float], alpha: float = 0.05
) -> List[TestResult]:
    """Aplica las cinco pruebas a una muestra y devuelve los resultados."""
    return [
        chi_square_test(sample, bins=10, alpha=alpha),
        runs_test(sample, alpha=alpha),
        kolmogorov_smirnov_test(sample, alpha=alpha),
        autocorrelation_test(sample, lag=1, alpha=alpha),
        poker_test(sample, digits=5, alpha=alpha),
    ]
