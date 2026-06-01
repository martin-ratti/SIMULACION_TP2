# Informe LaTeX - TP 2.2

Informe del TP 2.2 (Generadores de distribuciones de probabilidad).

## Archivos necesarios para compilar

El `.tex` referencia:

1. **Código fuente** vía `\lstinputlisting{codigo/...}` → la carpeta
   `codigo/` con `rng.py` y `distributions.py` (copiados del paquete
   `distribuciones/`). Ya están incluidos.
2. **Imágenes** → los 9 PNG de `../output/` (generados al correr
   `programa.py`), referenciados **sin prefijo de carpeta**.

## Compilación local

```bash
cd ..
python programa.py -n 10000 -s 12345 -a 0.05   # genera los PNG en output/
cp output/*.png informe/                        # copiar imagenes junto al .tex
cd informe
pdflatex -shell-escape informe.tex
```

> El paquete `listings` no requiere `-shell-escape`; se incluye por las
> dudas según el motor.

## Para Overleaf

Subir a la raíz del proyecto Overleaf:

- `informe.tex`
- la carpeta `codigo/` (con `rng.py` y `distributions.py`)
- los 9 PNG de `output/`: `uniforme.png`, `exponencial.png`, `gamma.png`,
  `normal.png`, `pascal.png`, `binomial.png`, `hipergeometrica.png`,
  `poisson.png`, `empirica.png`

Lo más cómodo es importar todo como ZIP (*New Project → Upload Project*).

## Contenido

- Introducción: del uniforme a distribuciones arbitrarias.
- Métodos generales: transformación inversa, rechazo, convolución.
- 9 distribuciones: teoría + derivación de F(x)/inversa + código + testeo.
- Resultados: tablas de media/varianza y de p-valores + gráficos.
- Discusión (incluye el bug de la Pascal y su corrección) y conclusiones.
