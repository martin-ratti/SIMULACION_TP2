# Informe LaTeX - TP 2.1

Informe del TP 2.1 (Generadores Pseudoaleatorios).

## Compilación

El informe incluye las imágenes de `../output/`, por lo que primero hay
que generarlas corriendo el programa desde la raíz del proyecto:

```bash
cd ..
python programa.py -n 10000 -s 12345 -a 0.05
```

Luego compilar el informe:

```bash
cd informe
pdflatex informe.tex
```

Las rutas de las figuras (`output/...`) son relativas a la raíz del
proyecto, así que conviene compilar desde la raíz o copiar la carpeta
`output/` junto al `.tex`. Alternativamente, en Overleaf subir el `.tex`
y la carpeta `output/` con los PNG.

## Contenido

- Introducción: pseudoaleatoriedad vs. aleatoriedad real.
- Materiales y métodos: GCL, Cuadrados Medios, XorShift, Mersenne Twister
  y las cinco pruebas estadísticas con sus fórmulas.
- Resultados: tablas de p-valores e imágenes (histogramas, scatter,
  comparativas).
- Discusión y conclusiones.
- Referencias.
