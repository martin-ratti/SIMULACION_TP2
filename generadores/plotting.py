"""Graficos para el analisis de generadores pseudoaleatorios (TP 2.1).

Genera, por cada generador:
    - Histograma de la muestra en [0, 1) frente a la densidad uniforme.
    - Scatter de pares (u_i, u_{i+1}) para evidenciar estructura/planos.

Y un grafico comparativo de p-valores del test chi-cuadrado entre todos
los generadores.
"""

from __future__ import annotations

import os
from typing import Dict, List, Sequence

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


def plot_histogram(
    sample: Sequence[float],
    title: str,
    output_dir: str,
    filename: str,
    bins: int = 20,
) -> str:
    """Histograma normalizado de la muestra vs. densidad uniforme = 1."""
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.hist(sample, bins=bins, range=(0, 1), density=True,
            color="steelblue", edgecolor="white", alpha=0.85,
            label="frecuencia observada")
    ax.axhline(1.0, color="crimson", linestyle="--", linewidth=1.4,
               label="densidad uniforme = 1")
    ax.set_xlabel("u (valor en [0, 1))")
    ax.set_ylabel("densidad")
    ax.set_title(title)
    ax.set_ylim(0, max(2.0, ax.get_ylim()[1]))
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8)
    return _save(fig, output_dir, filename)


def plot_scatter_pairs(
    sample: Sequence[float],
    title: str,
    output_dir: str,
    filename: str,
    max_points: int = 3000,
) -> str:
    """Scatter de pares consecutivos (u_i, u_{i+1}).

    Si los puntos se alinean en pocas rectas/planos, el generador tiene
    estructura latticial (planos de Marsaglia), sintoma de baja calidad.
    """
    xs = list(sample[:-1])[:max_points]
    ys = list(sample[1:])[:max_points]
    fig, ax = plt.subplots(figsize=(5.2, 5.2))
    ax.scatter(xs, ys, s=4, color="darkslategray", alpha=0.5)
    ax.set_xlabel("u_i")
    ax.set_ylabel("u_{i+1}")
    ax.set_title(title)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.25)
    return _save(fig, output_dir, filename)


def plot_pvalue_comparison(
    pvalues: Dict[str, float],
    title: str,
    output_dir: str,
    filename: str,
    alpha: float = 0.05,
) -> str:
    """Barras de p-valores (un test) por generador, con linea de alpha."""
    names = list(pvalues.keys())
    values = [pvalues[n] for n in names]
    colors = ["seagreen" if v >= alpha else "indianred" for v in values]

    fig, ax = plt.subplots(figsize=(8, 4.5))
    bars = ax.bar(names, values, color=colors, edgecolor="black", alpha=0.85)
    ax.axhline(alpha, color="black", linestyle="--", linewidth=1.3,
               label=f"alpha = {alpha}")
    ax.set_ylabel("p-valor")
    ax.set_title(title)
    ax.set_ylim(0, max(1.0, max(values) * 1.1) if values else 1.0)
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(loc="best", fontsize=8)
    for bar, v in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                f"{v:.3f}", ha="center", va="bottom", fontsize=8)
    plt.setp(ax.get_xticklabels(), rotation=15, ha="right")
    return _save(fig, output_dir, filename)
