"""Graficos de las distribuciones generadas (TP 2.2).

Por cada distribucion se superpone el histograma de la muestra generada
con la funcion teorica:
    - continuas: histograma normalizado + curva de densidad f(x).
    - discretas: barras de frecuencia relativa + puntos de la PMF teorica.
"""

from __future__ import annotations

import os
from typing import Callable, List, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def _ensure_dir(path: str) -> None:
    if path:
        os.makedirs(path, exist_ok=True)


def _save(fig, output_dir: str, filename: str) -> str:
    _ensure_dir(output_dir)
    target = os.path.join(output_dir, filename)
    fig.savefig(target, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return target


def plot_continuous(
    sample: Sequence[float],
    pdf: Callable[[float], float],
    title: str,
    output_dir: str,
    filename: str,
    bins: int = 40,
) -> str:
    """Histograma normalizado vs densidad teorica f(x)."""
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.hist(sample, bins=bins, density=True, color="steelblue",
            edgecolor="white", alpha=0.8, label="muestra generada")

    lo, hi = min(sample), max(sample)
    xs = [lo + (hi - lo) * i / 400 for i in range(401)]
    ys = [pdf(x) for x in xs]
    ax.plot(xs, ys, color="crimson", linewidth=2.0, label="densidad teorica f(x)")

    ax.set_xlabel("x")
    ax.set_ylabel("densidad")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=9)
    return _save(fig, output_dir, filename)


def plot_discrete(
    sample: Sequence[int],
    pmf: Callable[[int], float],
    title: str,
    output_dir: str,
    filename: str,
) -> str:
    """Frecuencia relativa observada vs PMF teorica."""
    smin, smax = min(sample), max(sample)
    support = list(range(smin, smax + 1))
    n = len(sample)

    counts = {k: 0 for k in support}
    for x in sample:
        counts[x] += 1
    rel_freq = [counts[k] / n for k in support]
    theo = [pmf(k) for k in support]

    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.bar(support, rel_freq, color="steelblue", edgecolor="white",
           alpha=0.8, label="frecuencia observada")
    ax.plot(support, theo, "o-", color="crimson", linewidth=1.5,
            markersize=5, label="PMF teorica")
    ax.set_xlabel("x")
    ax.set_ylabel("probabilidad")
    ax.set_title(title)
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(loc="best", fontsize=9)
    return _save(fig, output_dir, filename)
