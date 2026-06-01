# SIMULACION 2026 - TP 2.1: Generadores Pseudoaleatorios

Trabajo práctico de la cátedra de Simulación, UTN FRRO, 2026.

Implementación en Python 3 de varios generadores de números
pseudoaleatorios y una batería de cinco tests estadísticos para evaluar la
calidad de su generación, con comparación contra el RNG nativo de Python
(Mersenne Twister).

## Estructura

```
SIMULACION_TP2/
  programa.py          # CLI principal: corre generadores + tests + graficos
  generadores/
    __init__.py
    generators.py      # GCL, Cuadrados Medios, XorShift, PythonRNG
    tests.py           # Chi2, Rachas, KS, Autocorrelacion, Poker
    plotting.py        # Histogramas, scatter de pares, comparativas
  requirements.txt
  output/              # PNGs generados (gitignored)
  informe/             # Informe LaTeX
  README.md
```

## Generadores implementados

| Generador        | Tipo            | Descripción                                   |
|------------------|-----------------|-----------------------------------------------|
| GCL (LCG)        | Pseudoaleatorio | Congruencial lineal `X = (a·X + c) mod m`     |
| Cuadrados Medios | Pseudoaleatorio | Método de Von Neumann (histórico, débil)      |
| XorShift32       | Pseudoaleatorio | Desplazamientos y XOR de 32 bits (Marsaglia)  |
| Python (MT19937) | Referencia      | `random.Random`, Mersenne Twister             |

## Tests estadísticos

1. **Chi-cuadrado de uniformidad** — frecuencias por intervalos.
2. **Rachas (runs up/down)** — independencia de la secuencia.
3. **Kolmogorov-Smirnov** — ajuste a la CDF uniforme.
4. **Autocorrelación (lag 1)** — correlación serial.
5. **Póker (5 dígitos)** — patrones de repetición de dígitos.

Todos los tests se calculan sin dependencias externas (las distribuciones
chi-cuadrado y normal se aproximan con `math`).

## Uso

```bash
pip install -r requirements.txt
python programa.py -n 10000 -s 12345 -a 0.05
```

Opciones:

```
-n / --muestra      tamaño de muestra por generador (default 10000)
-s / --semilla      semilla común (default 12345)
-a / --alpha        nivel de significación (default 0.05)
--output-dir DIR    carpeta de gráficos (default ./output)
--no-plots          omitir gráficos
```

## Informe

El informe en LaTeX está en `informe/informe.tex`. Para compilarlo
(requiere los PNGs en `output/`, generados al correr `programa.py`):

```bash
cd informe
pdflatex informe.tex
```

## Referencias

- Enunciado: `SIMULACIÓN_2026_TP2_1.pdf` — Cátedra Simulación, UTN FRRO.
```
