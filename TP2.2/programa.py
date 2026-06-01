"""TP 2.2 - Generadores de distribuciones de probabilidad.

Genera valores de nueve distribuciones (uniforme, exponencial, gamma,
normal, pascal, binomial, hipergeometrica, poisson y empirica discreta) a
partir del GCL del TP 2.1, las testea segun lo que pide el enunciado y
produce los graficos comparando la muestra con la distribucion teorica.

CLI:

    python programa.py -n MUESTRA [-s SEMILLA] [-a ALPHA] [--output-dir DIR]
                       [--no-plots]
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from typing import List, Optional

from distribuciones import (
    LCG,
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
from distribuciones import tests as T
from distribuciones import plotting as P


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="programa.py",
        description="Generadores de distribuciones de probabilidad (TP 2.2).",
    )
    parser.add_argument("-n", "--muestra", type=int, default=10000,
                        help="tamano de muestra por distribucion")
    parser.add_argument("-s", "--semilla", type=int, default=12345,
                        help="semilla del GCL fuente de uniformes")
    parser.add_argument("-a", "--alpha", type=float, default=0.05,
                        help="nivel de significacion de los tests")
    parser.add_argument("--output-dir", default="output", dest="output_dir",
                        help="carpeta de salida para graficos")
    parser.add_argument("--no-plots", action="store_true", dest="no_plots",
                        help="no generar graficos")
    return parser.parse_args(argv)


def _fmt(value: float) -> str:
    return f"{value:.4f}"


def _print_row(name: str, gof) -> None:
    if gof is None:
        return
    print(f"\n--- {name} ---")
    print(f"  Media   : empirica={_fmt(gof.mean_emp):>10} | "
          f"teorica={_fmt(gof.mean_theo):>10}")
    print(f"  Varianza: empirica={_fmt(gof.var_emp):>10} | "
          f"teorica={_fmt(gof.var_theo):>10}")
    if gof.test_name != "(sin test exigido)":
        print(f"  {gof.test_name}: estadistico={_fmt(gof.statistic)} "
              f"| p-valor={_fmt(gof.p_value)} -> {gof.verdict}")
    else:
        print("  Testeo formal no exigido por el enunciado (Tabla 1).")
    if gof.detail:
        print(f"  {gof.detail}")


def _mean(sample) -> float:
    return sum(sample) / len(sample)


def _variance(sample) -> float:
    m = _mean(sample)
    return sum((x - m) ** 2 for x in sample) / len(sample)


def _build_gof(dist, sample, test_name, statistic, p_value, passed,
               detail=""):
    return T.GoodnessOfFit(
        dist_name=dist.name,
        test_name=test_name,
        statistic=statistic,
        p_value=p_value,
        passed=passed,
        mean_emp=_mean(sample),
        mean_theo=dist.mean_theoretical,
        var_emp=_variance(sample),
        var_theo=dist.var_theoretical,
        detail=detail,
    )


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    if args.muestra <= 0:
        raise SystemExit("ERROR: -n debe ser un entero positivo")
    if not 0 < args.alpha < 1:
        raise SystemExit("ERROR: -a (alpha) debe estar entre 0 y 1")

    n = args.muestra
    alpha = args.alpha
    out = args.output_dir

    print("== TP 2.2 Distribuciones de probabilidad ==")
    print(f"Muestra por distribucion : {n}")
    print(f"Semilla GCL (fuente U)   : {args.semilla}")
    print(f"Nivel alpha              : {alpha}")
    if not args.no_plots:
        os.makedirs(out, exist_ok=True)

    # Un GCL independiente por distribucion para que cada una arranque desde
    # el mismo estado y los resultados sean reproducibles y comparables.
    def rng():
        return LCG(seed=args.semilla)

    results = []

    # ---- CONTINUAS ----

    # Uniforme [2, 8] - testeo: SI (KS)
    d = Uniform(rng(), a=2.0, b=8.0)
    s = d.sample(n)
    dstat, pval, ok = T.ks_test(s, T.uniform_cdf(2.0, 8.0), alpha)
    gof = _build_gof(d, s, "Kolmogorov-Smirnov", dstat, pval, ok,
                     "Metodo: transformacion inversa x=a+(b-a)r")
    results.append(("uniforme", gof))
    if not args.no_plots:
        P.plot_continuous(s, lambda x: 1.0 / (8.0 - 2.0),
                          "Uniforme [2, 8]", out, "uniforme.png")

    # Exponencial media=2 - testeo: SI (KS)
    d = Exponential(rng(), mean=2.0)
    s = d.sample(n)
    dstat, pval, ok = T.ks_test(s, T.exponential_cdf(2.0), alpha)
    gof = _build_gof(d, s, "Kolmogorov-Smirnov", dstat, pval, ok,
                     "Metodo: transformacion inversa x=-EX*ln(r)")
    results.append(("exponencial", gof))
    if not args.no_plots:
        alpha_e = 1.0 / 2.0
        P.plot_continuous(s, lambda x: alpha_e * math.exp(-alpha_e * x),
                          "Exponencial (media=2)", out, "exponencial.png")

    # Gamma k=3, media=6 - testeo: NO exigido
    d = Gamma(rng(), k=3, mean=6.0)
    s = d.sample(n)
    gof = _build_gof(d, s, "(sin test exigido)", 0.0, 1.0, True,
                     "Metodo: convolucion de k=3 exponenciales")
    results.append(("gamma", gof))
    if not args.no_plots:
        k, al = 3, 3 / 6.0
        def gamma_pdf(x, k=k, al=al):
            if x <= 0:
                return 0.0
            return (al ** k) * (x ** (k - 1)) * math.exp(-al * x) \
                / math.factorial(k - 1)
        P.plot_continuous(s, gamma_pdf, "Gamma (k=3, media=6)", out,
                          "gamma.png")

    # Normal mu=10, sigma=2 - testeo: SI (KS)
    d = Normal(rng(), mu=10.0, sigma=2.0)
    s = d.sample(n)
    dstat, pval, ok = T.ks_test(s, T.normal_cdf(10.0, 2.0), alpha)
    gof = _build_gof(d, s, "Kolmogorov-Smirnov", dstat, pval, ok,
                     "Metodo: Box-Muller")
    results.append(("normal", gof))
    if not args.no_plots:
        mu, sg = 10.0, 2.0
        def normal_pdf(x, mu=mu, sg=sg):
            return math.exp(-0.5 * ((x - mu) / sg) ** 2) \
                / (sg * math.sqrt(2 * math.pi))
        P.plot_continuous(s, normal_pdf, "Normal (mu=10, sigma=2)", out,
                          "normal.png")

    # ---- DISCRETAS ----

    # Pascal k=3, p=0.4 - testeo: NO exigido
    d = Pascal(rng(), k=3, p=0.4)
    s = d.sample(n)
    gof = _build_gof(d, s, "(sin test exigido)", 0.0, 1.0, True,
                     "Metodo: suma de k=3 geometricas")
    results.append(("pascal", gof))
    if not args.no_plots:
        k, p = 3, 0.4
        def pascal_pmf(x, k=k, p=p):
            if x < 0:
                return 0.0
            return math.comb(k + x - 1, x) * (p ** k) * ((1 - p) ** x)
        P.plot_discrete(s, pascal_pmf, "Pascal (k=3, p=0.4)", out,
                        "pascal.png")

    # Binomial n=20, p=0.3 - testeo: SI (chi2)
    d = Binomial(rng(), n=20, p=0.3)
    s = d.sample(n)
    pmf = T.binomial_pmf(20, 0.3)
    chi2, pval, ok, dof = T.chi2_discrete_test(s, pmf, list(range(0, 21)),
                                               alpha)
    gof = _build_gof(d, s, "Chi-cuadrado", chi2, pval, ok,
                     f"Metodo: n ensayos Bernoulli | g.l.={dof}")
    results.append(("binomial", gof))
    if not args.no_plots:
        P.plot_discrete(s, pmf, "Binomial (n=20, p=0.3)", out, "binomial.png")

    # Hipergeometrica N=100, p=0.4, draws=20 - testeo: NO exigido
    d = Hypergeometric(rng(), N=100, p=0.4, draws=20)
    s = d.sample(n)
    gof = _build_gof(d, s, "(sin test exigido)", 0.0, 1.0, True,
                     "Metodo: muestreo sin reemplazo")
    results.append(("hipergeometrica", gof))
    if not args.no_plots:
        N, pp, dr = 100, 0.4, 20
        Np, Nq = int(N * pp), int(N * (1 - pp))
        def hyper_pmf(x, Np=Np, Nq=Nq, N=N, dr=dr):
            if x < 0 or x > dr or x > Np or (dr - x) > Nq:
                return 0.0
            return (math.comb(Np, x) * math.comb(Nq, dr - x)) / math.comb(N, dr)
        P.plot_discrete(s, hyper_pmf, "Hipergeometrica (N=100, p=0.4, n=20)",
                        out, "hipergeometrica.png")

    # Poisson lambda=5 - testeo: SI (chi2)
    d = Poisson(rng(), lam=5.0)
    s = d.sample(n)
    pmf = T.poisson_pmf(5.0)
    smax = max(max(s), 20)
    chi2, pval, ok, dof = T.chi2_discrete_test(s, pmf, list(range(0, smax + 1)),
                                               alpha)
    gof = _build_gof(d, s, "Chi-cuadrado", chi2, pval, ok,
                     f"Metodo: producto de uniformes | g.l.={dof}")
    results.append(("poisson", gof))
    if not args.no_plots:
        P.plot_discrete(s, pmf, "Poisson (lambda=5)", out, "poisson.png")

    # Empirica discreta - testeo: SI (chi2)
    values = [1, 2, 3, 4, 5]
    probs = [0.10, 0.25, 0.30, 0.20, 0.15]
    d = EmpiricalDiscrete(rng(), values=values, probs=probs)
    s = d.sample(n)
    pmf = T.empirical_pmf(values, probs)
    chi2, pval, ok, dof = T.chi2_discrete_test(s, pmf, values, alpha)
    gof = _build_gof(d, s, "Chi-cuadrado", chi2, pval, ok,
                     f"Metodo: inversa sobre acumulada | g.l.={dof}")
    results.append(("empirica", gof))
    if not args.no_plots:
        P.plot_discrete(s, pmf, "Empirica discreta", out, "empirica.png")

    # ---- REPORTE ----
    for _, gof in results:
        _print_row(gof.dist_name, gof)

    print("\n== Resumen de tests exigidos por el enunciado ==")
    for _, gof in results:
        if gof.test_name != "(sin test exigido)":
            print(f"  {gof.dist_name:<18} {gof.test_name:<20} {gof.verdict}")

    if not args.no_plots:
        print(f"\nGraficos en: {os.path.abspath(out)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
