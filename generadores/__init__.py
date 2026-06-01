"""Paquete de generadores pseudoaleatorios y tests estadisticos (TP 2.1).

Expone los generadores implementados y el motor de tests para poder
importarlos comodamente desde ``programa.py``::

    from generadores import LCG, MiddleSquare, XorShift, PythonRNG
    from generadores import run_all_tests

Todos los generadores implementan la misma interfaz minima:

    - ``next_int()``   -> entero crudo del estado interno
    - ``next_float()`` -> flotante uniforme en [0, 1)
    - ``sample(n)``    -> lista de ``n`` flotantes en [0, 1)
"""

from __future__ import annotations

from .generators import (
    Generator,
    LCG,
    MiddleSquare,
    XorShift,
    PythonRNG,
    build_generators,
)
from .tests import (
    TestResult,
    chi_square_test,
    runs_test,
    kolmogorov_smirnov_test,
    autocorrelation_test,
    poker_test,
    run_all_tests,
)

__all__ = [
    "Generator",
    "LCG",
    "MiddleSquare",
    "XorShift",
    "PythonRNG",
    "build_generators",
    "TestResult",
    "chi_square_test",
    "runs_test",
    "kolmogorov_smirnov_test",
    "autocorrelation_test",
    "poker_test",
    "run_all_tests",
]
