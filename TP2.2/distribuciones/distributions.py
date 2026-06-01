"""Generadores de distintas distribuciones de probabilidad (TP 2.2).

Cada distribucion se construye a partir de numeros uniformes U(0, 1)
producidos por el GCL del TP 2.1 (modulo ``rng``). Los metodos siguen el
capitulo 4 de Naylor, "Tecnicas de simulacion en computadoras" (1982).

Metodo empleado por distribucion (segun la Tabla 1 del enunciado):

    Continuas
        Uniforme       : transformacion inversa  x = a + (b-a)*r       (4-23)
        Exponencial    : transformacion inversa  x = -EX*ln(r)         (4-30)
        Gamma (Erlang) : convolucion de k exponenciales                (4-54)
        Normal         : metodo directo Box-Muller                     (4-81)

    Discretas
        Pascal         : suma de k geometricas                         (4-133)
        Binomial       : n ensayos de Bernoulli                        (4-136)
        Hipergeometrica: muestreo sin reemplazo (p se actualiza)       (4-144)
        Poisson        : producto de uniformes < e^{-lambda}           (4-153)
        Empirica disc. : transformacion inversa sobre la acumulada     (4-154)

Cada clase expone ``generate()`` (un valor) y ``sample(n)`` (lista de n
valores), y atributos ``mean_theoretical`` / ``var_theoretical`` con la
media y varianza teoricas para validacion.
"""

from __future__ import annotations

import math
from typing import List, Sequence

from .rng import LCG


# ===========================================================================
# DISTRIBUCIONES CONTINUAS
# ===========================================================================


class Uniform:
    """Distribucion uniforme continua en [a, b].

    Densidad: f(x) = 1/(b-a) para a <= x <= b.
    Generacion por transformacion inversa (Naylor 4-23):

        x = a + (b - a) * r,    r ~ U(0, 1).
    """

    name = "Uniforme"

    def __init__(self, rng: LCG, a: float = 0.0, b: float = 1.0) -> None:
        if b <= a:
            raise ValueError("se requiere b > a")
        self.rng = rng
        self.a = a
        self.b = b
        self.mean_theoretical = (a + b) / 2.0
        self.var_theoretical = (b - a) ** 2 / 12.0

    def generate(self) -> float:
        r = self.rng.uniform()
        return self.a + (self.b - self.a) * r

    def sample(self, n: int) -> List[float]:
        return [self.generate() for _ in range(n)]


class Exponential:
    """Distribucion exponencial con media EX = 1/alpha.

    Densidad: f(x) = alpha * e^{-alpha x}, x >= 0.
    Generacion por transformacion inversa (Naylor 4-30):

        x = -EX * ln(r),    r ~ U(0, 1).

    Se usa ``uniform_open`` para evitar ln(0).
    """

    name = "Exponencial"

    def __init__(self, rng: LCG, mean: float = 1.0) -> None:
        if mean <= 0:
            raise ValueError("la media debe ser positiva")
        self.rng = rng
        self.mean = mean  # EX = 1/alpha
        self.mean_theoretical = mean
        self.var_theoretical = mean ** 2

    def generate(self) -> float:
        r = self.rng.uniform_open()
        return -self.mean * math.log(r)

    def sample(self, n: int) -> List[float]:
        return [self.generate() for _ in range(n)]


class Gamma:
    """Distribucion gamma (Erlang) de parametro entero k y media EX.

    Para k entero, la gamma es la suma de k exponenciales de igual media.
    Naylor usa la forma computacional (4-54):

        x = -(1/alpha) * ln( producto_{i=1..k} r_i ),

    donde 1/alpha = EX/k es la media de cada exponencial componente. Asi la
    media total es EX y la varianza EX^2 / k.
    """

    name = "Gamma (Erlang)"

    def __init__(self, rng: LCG, k: int = 2, mean: float = 1.0) -> None:
        if k < 1:
            raise ValueError("k debe ser un entero >= 1")
        if mean <= 0:
            raise ValueError("la media debe ser positiva")
        self.rng = rng
        self.k = k
        self.mean = mean
        self.alpha = k / mean  # tasa de cada exponencial componente
        self.mean_theoretical = mean
        self.var_theoretical = mean ** 2 / k

    def generate(self) -> float:
        prod = 1.0
        for _ in range(self.k):
            prod *= self.rng.uniform_open()
        return -(1.0 / self.alpha) * math.log(prod)

    def sample(self, n: int) -> List[float]:
        return [self.generate() for _ in range(n)]


class Normal:
    """Distribucion normal de media mu y desvio sigma.

    Metodo directo de Box-Muller (Naylor 4-81/4-82): a partir de dos
    uniformes independientes r1, r2 se obtienen dos normales estandar

        z1 = sqrt(-2 ln r1) * cos(2 pi r2)
        z2 = sqrt(-2 ln r1) * sin(2 pi r2)

    y luego x = mu + sigma * z. Se cachea z2 para no desperdiciar la
    segunda normal que produce cada par.
    """

    name = "Normal"

    def __init__(self, rng: LCG, mu: float = 0.0, sigma: float = 1.0) -> None:
        if sigma <= 0:
            raise ValueError("sigma debe ser positivo")
        self.rng = rng
        self.mu = mu
        self.sigma = sigma
        self.mean_theoretical = mu
        self.var_theoretical = sigma ** 2
        self._cached = None  # segunda normal estandar pendiente

    def generate(self) -> float:
        if self._cached is not None:
            z = self._cached
            self._cached = None
            return self.mu + self.sigma * z
        r1 = self.rng.uniform_open()
        r2 = self.rng.uniform()
        radius = math.sqrt(-2.0 * math.log(r1))
        z1 = radius * math.cos(2.0 * math.pi * r2)
        z2 = radius * math.sin(2.0 * math.pi * r2)
        self._cached = z2
        return self.mu + self.sigma * z1

    def sample(self, n: int) -> List[float]:
        return [self.generate() for _ in range(n)]


# ===========================================================================
# DISTRIBUCIONES DISCRETAS
# ===========================================================================


class Pascal:
    """Distribucion de Pascal (binomial negativa con k entero).

    Cuenta el numero de fallas antes de obtener k exitos en ensayos de
    Bernoulli con probabilidad p. Es la suma de k variables geometricas.

    Naylor (4-125) genera una geometrica por transformacion inversa como
    x_g = floor( ln(r) / ln(q) ), q = 1 - p, redondeando "al entero menor".
    La Pascal es la suma de k de estas geometricas:

        x = sum_{i=1..k} floor( ln(r_i) / ln(q) ).

    Importante: el redondeo (floor) se aplica a CADA geometrica por
    separado; aplicarlo una sola vez a la suma (como sugiere la forma
    compacta 4-133) sesga la media, porque floor(a)+floor(b) != floor(a+b).

    Media EX = k*q/p, varianza VX = k*q/p^2.
    """

    name = "Pascal"

    def __init__(self, rng: LCG, k: int = 3, p: float = 0.5) -> None:
        if k < 1:
            raise ValueError("k debe ser un entero >= 1")
        if not 0 < p < 1:
            raise ValueError("p debe estar en (0, 1)")
        self.rng = rng
        self.k = k
        self.p = p
        self.q = 1.0 - p
        self._log_q = math.log(self.q)
        self.mean_theoretical = k * self.q / p
        self.var_theoretical = k * self.q / p ** 2

    def generate(self) -> int:
        x = 0
        for _ in range(self.k):
            r = self.rng.uniform_open()
            x += int(math.floor(math.log(r) / self._log_q))
        return x

    def sample(self, n: int) -> List[int]:
        return [self.generate() for _ in range(n)]


class Binomial:
    """Distribucion binomial: numero de exitos en n ensayos de Bernoulli.

    Naylor (figuras 4-19/4-20) la genera reproduciendo n ensayos: por cada
    ensayo se genera r ~ U(0,1) y se cuenta un exito si r <= p.

    Media EX = n*p, varianza VX = n*p*q.
    """

    name = "Binomial"

    def __init__(self, rng: LCG, n: int = 10, p: float = 0.5) -> None:
        if n < 1:
            raise ValueError("n debe ser un entero >= 1")
        if not 0 <= p <= 1:
            raise ValueError("p debe estar en [0, 1]")
        self.rng = rng
        self.n = n
        self.p = p
        self.mean_theoretical = n * p
        self.var_theoretical = n * p * (1.0 - p)

    def generate(self) -> int:
        exitos = 0
        for _ in range(self.n):
            if self.rng.uniform() <= self.p:
                exitos += 1
        return exitos

    def sample(self, n: int) -> List[int]:
        return [self.generate() for _ in range(n)]


class Hypergeometric:
    """Distribucion hipergeometrica: muestreo sin reemplazo.

    Poblacion de N elementos, de los cuales N*p son de la clase I. Se
    extrae una muestra de ``draws`` elementos sin reemplazo y se cuenta
    cuantos son de la clase I. Naylor (figura 4-22) actualiza p tras cada
    extraccion:

        si el elemento extraido es de clase I (r <= p): S = 1, x += 1
        si no:                                          S = 0
        p_nuevo = (N*p - S) / (N - 1);   N_nuevo = N - 1

    Media EX = n*p, varianza VX = n*p*q*(N-n)/(N-1).
    """

    name = "Hipergeometrica"

    def __init__(
        self, rng: LCG, N: int = 100, p: float = 0.5, draws: int = 10
    ) -> None:
        if N < 1:
            raise ValueError("N debe ser un entero >= 1")
        if not 0 <= p <= 1:
            raise ValueError("p debe estar en [0, 1]")
        if not 1 <= draws <= N:
            raise ValueError("draws debe estar entre 1 y N")
        self.rng = rng
        self.N = N
        self.p = p
        self.draws = draws
        q = 1.0 - p
        self.mean_theoretical = draws * p
        self.var_theoretical = draws * p * q * (N - draws) / (N - 1) \
            if N > 1 else 0.0

    def generate(self) -> int:
        N = float(self.N)
        p = self.p
        x = 0
        for _ in range(self.draws):
            r = self.rng.uniform()
            if r <= p:
                s = 1.0
                x += 1
            else:
                s = 0.0
            # Actualizar p y N para la siguiente extraccion (sin reemplazo).
            if N > 1:
                p = (N * p - s) / (N - 1.0)
                p = min(1.0, max(0.0, p))  # evitar derivas numericas
            N -= 1.0
        return x

    def sample(self, n: int) -> List[int]:
        return [self.generate() for _ in range(n)]


class Poisson:
    """Distribucion de Poisson de parametro lambda.

    Cuenta el numero de eventos en un intervalo unitario. Naylor (4-153)
    usa la relacion con la exponencial: se multiplican uniformes hasta que
    el producto cae por debajo de e^{-lambda}:

        se incrementa x mientras  producto_{i} r_i >= e^{-lambda}.

    Media EX = lambda, varianza VX = lambda.
    """

    name = "Poisson"

    def __init__(self, rng: LCG, lam: float = 4.0) -> None:
        if lam <= 0:
            raise ValueError("lambda debe ser positivo")
        self.rng = rng
        self.lam = lam
        self._threshold = math.exp(-lam)  # B = e^{-lambda}
        self.mean_theoretical = lam
        self.var_theoretical = lam

    def generate(self) -> int:
        x = 0
        prod = self.rng.uniform()
        while prod >= self._threshold:
            prod *= self.rng.uniform()
            x += 1
        return x

    def sample(self, n: int) -> List[int]:
        return [self.generate() for _ in range(n)]


class EmpiricalDiscrete:
    """Distribucion empirica discreta arbitraria.

    Dada una tabla de valores ``values`` con probabilidades ``probs`` (que
    suman 1), genera por transformacion inversa sobre la acumulada
    (Naylor 4-154): se genera r ~ U(0,1) y se devuelve el valor b_i tal que

        P_1 + ... + P_{i-1} < r <= P_1 + ... + P_i.
    """

    name = "Empirica discreta"

    def __init__(
        self, rng: LCG, values: Sequence, probs: Sequence[float]
    ) -> None:
        if len(values) != len(probs):
            raise ValueError("values y probs deben tener igual longitud")
        total = sum(probs)
        if abs(total - 1.0) > 1e-6:
            raise ValueError(f"las probabilidades deben sumar 1 (suman {total})")
        self.rng = rng
        self.values = list(values)
        self.probs = list(probs)
        # Acumulada para la busqueda inversa.
        self.cumulative: List[float] = []
        acc = 0.0
        for pr in self.probs:
            acc += pr
            self.cumulative.append(acc)
        self.mean_theoretical = sum(v * pr for v, pr in zip(values, probs))
        ex2 = sum((v ** 2) * pr for v, pr in zip(values, probs))
        self.var_theoretical = ex2 - self.mean_theoretical ** 2

    def generate(self):
        r = self.rng.uniform()
        for value, cum in zip(self.values, self.cumulative):
            if r <= cum:
                return value
        return self.values[-1]  # salvaguarda por redondeo

    def sample(self, n: int) -> List:
        return [self.generate() for _ in range(n)]
