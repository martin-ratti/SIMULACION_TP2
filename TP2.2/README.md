# TP 2.2 - Generadores de distribuciones de probabilidad

Genera valores de nueve distribuciones de probabilidad a partir del GCL
validado en el TP 2.1, siguiendo los métodos del capítulo 4 de Naylor
(*Técnicas de simulación en computadoras*, 1982).

## Estructura

```
TP2.2/
  programa.py               # CLI: genera, testea y grafica las 9 distribuciones
  distribuciones/
    __init__.py
    rng.py                  # GCL (fuente de uniformes, reusado del TP2.1)
    distributions.py        # las 9 distribuciones
    tests.py                # KS, chi-cuadrado, media/varianza teórica
    plotting.py             # histogramas + curva/PMF teórica
  requirements.txt
  output/                   # PNGs generados (gitignored)
  informe/                  # informe LaTeX
```

## Distribuciones y métodos

| Distribución      | Tipo     | Método (Naylor)                    | Testeo exigido |
|-------------------|----------|------------------------------------|----------------|
| Uniforme          | continua | Inversa `x=a+(b-a)r`               | KS             |
| Exponencial       | continua | Inversa `x=-EX·ln(r)`             | KS             |
| Gamma (Erlang)    | continua | Convolución de k exponenciales     | no             |
| Normal            | continua | Box-Muller                         | KS             |
| Pascal            | discreta | Suma de k geométricas              | no             |
| Binomial          | discreta | n ensayos de Bernoulli             | chi²           |
| Hipergeométrica   | discreta | Muestreo sin reemplazo             | no             |
| Poisson           | discreta | Producto de uniformes              | chi²           |
| Empírica discreta | discreta | Inversa sobre la acumulada         | chi²           |

El testeo formal sigue la columna "Testeo" de la Tabla 1 del enunciado.
Para todas las distribuciones se compara además media y varianza empírica
contra la teórica.

## Uso

```bash
pip install -r requirements.txt
python programa.py -n 10000 -s 12345 -a 0.05
```

Opciones:

```
-n / --muestra      tamaño de muestra por distribución (default 10000)
-s / --semilla      semilla del GCL fuente de uniformes (default 12345)
-a / --alpha        nivel de significación (default 0.05)
--output-dir DIR    carpeta de gráficos (default ./output)
--no-plots          omitir gráficos
```

## Informe

`informe/informe.tex`. Requiere los PNG de `output/` (generados al correr
`programa.py`). El informe incluye el código fuente embebido (paquete
`listings`).

## Referencias

- Naylor, T. H. (1982). *Técnicas de simulación en computadoras*, cap. 4.
- Enunciado: `SIMULACIÓN_2026_TP2_2.pdf` — Cátedra Simulación, UTN FRRO.
```
