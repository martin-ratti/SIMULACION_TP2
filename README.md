# SIMULACION 2026 - TP 2

Trabajos prácticos de la cátedra de Simulación, UTN FRRO, 2026.

## Estructura del repositorio

```
SIMULACION_TP2/
  TP2.1/        # Generadores pseudoaleatorios + tests estadisticos
    programa.py
    generadores/
    output/     # PNGs generados (gitignored)
    informe/    # informe LaTeX
    requirements.txt
    README.md

  TP2.2/        # Generadores de distribuciones de probabilidad
    programa.py
    distribuciones/
    output/     # PNGs generados (gitignored)
    informe/    # informe LaTeX (incluye codigo embebido)
    requirements.txt
    README.md
```

Cada TP es independiente: tiene su propio `programa.py`, paquete,
`requirements.txt`, `output/` e `informe/`.

## TP 2.1 - Generadores pseudoaleatorios

Implementa cuatro generadores (GCL, Cuadrados Medios, XorShift y el
Mersenne Twister de Python) y los evalúa con cinco tests estadísticos
(chi-cuadrado, rachas, Kolmogorov-Smirnov, autocorrelación, póker).

```bash
cd TP2.1
pip install -r requirements.txt
python programa.py -n 10000 -s 12345 -e 17
```

Detalles en `TP2.1/README.md`.

## TP 2.2 - Distribuciones de probabilidad

Construye nueve distribuciones (uniforme, exponencial, gamma, normal,
Pascal, binomial, hipergeométrica, Poisson y empírica discreta) a partir
del GCL del TP 2.1, siguiendo el capítulo 4 de Naylor. Testea con
Kolmogorov-Smirnov (continuas) y chi-cuadrado (discretas).

```bash
cd TP2.2
pip install -r requirements.txt
python programa.py -n 10000 -s 12345 -a 0.05
```

Detalles en `TP2.2/README.md`.

## Referencias

- Enunciados:
  - `SIMULACIÓN_2026_TP2_1.pdf`
  - `SIMULACIÓN_2026_TP2_2.pdf`
- Naylor, T. H. (1982). *Técnicas de simulación en computadoras*, cap. 4.
- Cátedra Simulación, UTN FRRO, 2026.
```
