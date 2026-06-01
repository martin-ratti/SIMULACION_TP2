"""Generadores de numeros pseudoaleatorios (TP 2.1).

Se implementan tres generadores deterministas mas un envoltorio sobre el
generador nativo de Python (Mersenne Twister) usado como referencia:

    - ``LCG``          : Generador Congruencial Lineal (obligatorio).
    - ``MiddleSquare`` : Metodo de los Cuadrados Medios (Von Neumann).
    - ``XorShift``     : Generador XorShift de 32 bits.
    - ``PythonRNG``    : ``random.Random`` de la libreria estandar.

Todos producen flotantes uniformes en ``[0, 1)`` via ``next_float`` y una
muestra de tamano ``n`` via ``sample(n)``. La normalizacion a ``[0, 1)``
divide el entero crudo por el modulo (o rango) correspondiente, de modo que
los tests estadisticos pueden tratarlos de forma homogenea.
"""

from __future__ import annotations

import random
from abc import ABC, abstractmethod
from typing import List, Optional


class Generator(ABC):
    """Interfaz comun a todos los generadores pseudoaleatorios."""

    name: str = "generic"
    kind: str = "pseudoaleatorio"

    @abstractmethod
    def next_int(self) -> int:
        """Devuelve el siguiente entero crudo del estado interno."""

    @abstractmethod
    def next_float(self) -> float:
        """Devuelve el siguiente flotante uniforme en [0, 1)."""

    def sample(self, n: int) -> List[float]:
        """Genera ``n`` flotantes uniformes en [0, 1)."""
        if n <= 0:
            raise ValueError("n debe ser un entero positivo")
        return [self.next_float() for _ in range(n)]


class LCG(Generator):
    """Generador Congruencial Lineal (GCL / LCG).

    Recurrencia:

        X_{n+1} = (a * X_n + c) mod m

    El flotante en [0, 1) se obtiene como ``X_{n+1} / m``.

    Los parametros por defecto son los del generador ``MINSTD`` revisado
    (Park & Miller, 1993): ``a = 48271``, ``c = 0``, ``m = 2**31 - 1``.
    Es un generador congruencial multiplicativo de periodo completo
    ``m - 1`` y de calidad estadistica conocida, ampliamente documentado.
    """

    name = "GCL (LCG)"
    kind = "pseudoaleatorio"

    def __init__(
        self,
        seed: int = 12345,
        a: int = 48271,
        c: int = 0,
        m: int = 2 ** 31 - 1,
    ) -> None:
        if m <= 0:
            raise ValueError("el modulo m debe ser positivo")
        self.a = a
        self.c = c
        self.m = m
        self._state = seed % m
        if self._state == 0 and c == 0:
            # Un estado 0 en un GCL multiplicativo queda atrapado en 0.
            self._state = 1

    def next_int(self) -> int:
        self._state = (self.a * self._state + self.c) % self.m
        return self._state

    def next_float(self) -> float:
        return self.next_int() / self.m


class MiddleSquare(Generator):
    """Metodo de los Cuadrados Medios de Von Neumann.

    Se eleva al cuadrado la semilla, se rellena con ceros a la izquierda
    hasta ``2 * digits`` cifras y se extraen los ``digits`` digitos
    centrales, que constituyen el nuevo estado. El flotante en [0, 1) se
    obtiene dividiendo por ``10 ** digits``.

    Es un generador historico y deliberadamente debil: tiende a degenerar
    (caer en ciclos cortos o en cero) cuando los digitos centrales se
    vuelven pequenos. Se incluye precisamente para evidenciar como un mal
    generador falla los tests estadisticos.
    """

    name = "Cuadrados Medios"
    kind = "pseudoaleatorio"

    def __init__(self, seed: int = 675248, digits: int = 6) -> None:
        if digits % 2 != 0:
            raise ValueError("digits debe ser par para tomar el centro")
        self.digits = digits
        self._state = seed % (10 ** digits)
        if self._state == 0:
            self._state = 10 ** (digits // 2) + 1

    def next_int(self) -> int:
        squared = self._state ** 2
        # Rellenar a 2*digits cifras y tomar el bloque central.
        text = str(squared).zfill(2 * self.digits)
        start = (len(text) - self.digits) // 2
        center = text[start:start + self.digits]
        self._state = int(center)
        return self._state

    def next_float(self) -> float:
        return self.next_int() / (10 ** self.digits)


class XorShift(Generator):
    """Generador XorShift de 32 bits (Marsaglia, 2003).

    Aplica una secuencia de desplazamientos y XOR sobre un estado de 32
    bits. Es rapido y de mejor calidad que un GCL simple, con periodo
    ``2**32 - 1``. Se usa como tercer generador propio para la comparacion.
    """

    name = "XorShift32"
    kind = "pseudoaleatorio"
    MASK = 0xFFFFFFFF
    DIVISOR = float(0x100000000)  # 2**32

    def __init__(self, seed: int = 2463534242) -> None:
        self._state = seed & self.MASK
        if self._state == 0:
            self._state = 2463534242  # estado 0 es punto fijo: evitarlo

    def next_int(self) -> int:
        x = self._state
        x ^= (x << 13) & self.MASK
        x ^= (x >> 17)
        x ^= (x << 5) & self.MASK
        self._state = x & self.MASK
        return self._state

    def next_float(self) -> float:
        return self.next_int() / self.DIVISOR


class PythonRNG(Generator):
    """Envoltorio sobre ``random.Random`` (Mersenne Twister) de Python.

    Sirve como generador de referencia de alta calidad para contrastar con
    los generadores implementados a mano.
    """

    name = "Python (MT19937)"
    kind = "referencia"

    def __init__(self, seed: Optional[int] = 12345) -> None:
        self._rng = random.Random(seed)

    def next_int(self) -> int:
        return self._rng.getrandbits(32)

    def next_float(self) -> float:
        return self._rng.random()


def build_generators(seed: int = 12345) -> List[Generator]:
    """Construye la lista estandar de generadores a comparar.

    Se usa una semilla comun (donde aplica) para que las corridas sean
    reproducibles. Los cuadrados medios reciben una semilla de 6 cifras
    distinta de cero para arrancar con suficientes digitos significativos.
    """
    return [
        LCG(seed=seed),
        MiddleSquare(seed=675248),
        XorShift(seed=(seed * 2654435761) & 0xFFFFFFFF or 2463534242),
        PythonRNG(seed=seed),
    ]
