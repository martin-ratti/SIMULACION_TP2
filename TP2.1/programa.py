"""TP 2.1 - Generadores pseudoaleatorios.

Genera muestras con varios generadores (GCL, Cuadrados Medios, XorShift y
el RNG nativo de Python), les aplica cinco tests estadisticos y produce
una tabla comparativa por consola mas los graficos del informe.

CLI:

    python programa.py -n MUESTRA [-s SEMILLA] [-a ALPHA] [--output-dir DIR]

Donde:
    -n MUESTRA      tamano de la muestra por generador (default 10000)
    -s SEMILLA      semilla comun para los generadores (default 12345)
    -a ALPHA        nivel de significacion de los tests (default 0.05)
    --output-dir    carpeta de salida para graficos (default ./output)
    --no-plots      no generar graficos (solo tabla por consola)
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, List, Optional

from generadores import build_generators, run_all_tests
from generadores.generators import Generator
from generadores.tests import TestResult


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="programa.py",
        description="Comparacion de generadores pseudoaleatorios (TP 2.1).",
    )
    parser.add_argument("-n", "--muestra", type=int, default=10000,
                        help="tamano de la muestra por generador")
    parser.add_argument("-s", "--semilla", type=int, default=12345,
                        help="semilla comun para los generadores")
    parser.add_argument("-a", "--alpha", type=float, default=0.05,
                        help="nivel de significacion de los tests")
    parser.add_argument("--output-dir", default="output", dest="output_dir",
                        help="carpeta de salida para graficos")
    parser.add_argument("--no-plots", action="store_true", dest="no_plots",
                        help="no generar graficos")
    return parser.parse_args(argv)


def _print_generator_block(gen: Generator, results: List[TestResult]) -> None:
    print(f"\n=== {gen.name}  [{gen.kind}] ===")
    header = f"  {'Test':<28} {'Estadistico':>12} {'p-valor':>10}  Veredicto"
    print(header)
    print("  " + "-" * (len(header) - 2))
    for r in results:
        print(f"  {r.name:<28} {r.statistic:>12.4f} {r.p_value:>10.4f}"
              f"  {r.verdict}")


def _print_summary(table: Dict[str, List[TestResult]]) -> None:
    """Tabla resumen: cuantos tests paso cada generador."""
    print("\n=== Resumen (tests superados / total) ===")
    for name, results in table.items():
        passed = sum(1 for r in results if r.passed)
        print(f"  {name:<22} {passed}/{len(results)}")


def _make_plots(
    samples: Dict[str, List[float]],
    table: Dict[str, List[TestResult]],
    output_dir: str,
    alpha: float,
) -> None:
    # Importar aca para no exigir matplotlib si se usa --no-plots.
    from generadores import plotting

    os.makedirs(output_dir, exist_ok=True)
    slug = {
        "GCL (LCG)": "gcl",
        "Cuadrados Medios": "midsq",
        "XorShift32": "xorshift",
        "Python (MT19937)": "python",
    }
    for name, sample in samples.items():
        key = slug.get(name, name.lower().replace(" ", "_"))
        plotting.plot_histogram(
            sample, f"Histograma - {name}", output_dir, f"hist_{key}.png")
        plotting.plot_scatter_pairs(
            sample, f"Pares (u_i, u_i+1) - {name}", output_dir,
            f"scatter_{key}.png")

    # Comparacion de p-valores del chi-cuadrado (primer test de la lista).
    chi_pvalues = {name: results[0].p_value for name, results in table.items()}
    plotting.plot_pvalue_comparison(
        chi_pvalues, "p-valor Chi-cuadrado por generador",
        output_dir, "comparacion_chi2.png", alpha=alpha)

    # Comparacion de p-valores de Kolmogorov-Smirnov (tercer test).
    ks_pvalues = {name: results[2].p_value for name, results in table.items()}
    plotting.plot_pvalue_comparison(
        ks_pvalues, "p-valor Kolmogorov-Smirnov por generador",
        output_dir, "comparacion_ks.png", alpha=alpha)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    if args.muestra <= 0:
        raise SystemExit("ERROR: -n debe ser un entero positivo")
    if not 0 < args.alpha < 1:
        raise SystemExit("ERROR: -a (alpha) debe estar entre 0 y 1")

    print("== TP 2.1 Generadores pseudoaleatorios ==")
    print(f"Muestra por generador : {args.muestra}")
    print(f"Semilla comun         : {args.semilla}")
    print(f"Nivel alpha           : {args.alpha}")

    generators = build_generators(seed=args.semilla)
    samples: Dict[str, List[float]] = {}
    table: Dict[str, List[TestResult]] = {}

    for gen in generators:
        sample = gen.sample(args.muestra)
        results = run_all_tests(sample, alpha=args.alpha)
        samples[gen.name] = sample
        table[gen.name] = results
        _print_generator_block(gen, results)

    _print_summary(table)

    if not args.no_plots:
        _make_plots(samples, table, args.output_dir, args.alpha)
        print(f"\nGraficos en: {os.path.abspath(args.output_dir)}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
