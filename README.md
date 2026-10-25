# Antena Coil Generator

Generador del soporte imprimible (STL) de una **antena solenoide**: un cilindro hueco con un
surco helicoidal en U donde se aloja el cable. Las vueltas se reparten de forma no uniforme
(más juntas en los extremos que en el centro) para compensar el efecto de borde y obtener un
campo axial homogéneo en la región que se elija. Pensado para antenas de RF de resonancia
magnética de campo bajo, pero sirve para cualquier solenoide que necesite campo uniforme.

## Uso rápido (Windows)

1. Clone el repositorio o descárguelo como ZIP y descomprímalo.
2. Haga doble clic en **`iniciar_gui.bat`**.
3. La primera vez:
   - si no hay **Python 3.10 o superior**, una ventana pregunta si se instala Python 3.12
     (con `winget`, el instalador de Windows);
   - otra ventana pregunta si se instalan las dependencias (numpy, scipy, matplotlib y
     manifold3d, unos 100 MB). Se instalan en la carpeta `.venv` del proyecto, sin tocar el
     resto del sistema.
4. Se abre la interfaz en el navegador. La ventana de consola debe quedar abierta mientras
   se usa; al cerrarla se detiene el generador.

La interfaz sólo escucha en `127.0.0.1`: no es accesible desde otros equipos.

## La interfaz

Cada generación exige un **nombre de proyecto** y deja todos los archivos en
`proyectos/<nombre>/` (si la carpeta ya existe, se reemplazan sus archivos). Al pasar el
cursor sobre cada parámetro aparece una nota con su explicación y su valor por defecto.

| Parámetro | Por defecto | Qué controla |
|---|---|---|
| Diámetro interno / externo | 100 / 104 mm | Cilindro impreso; deben quedar ≥ 0.5 mm de pared bajo el surco |
| Largo | 170 mm | Largo total del soporte |
| Margen en los extremos | 5 mm | Distancia mínima del borde a la primera y última vuelta |
| Diámetro del surco | 1.0 mm | Surco en U al ras de la superficie |
| Diámetro del cable | 1.0 mm | Conductor real; apoya en el fondo del surco |
| Largo máximo del cable | 7.5 m | Fija el número de vueltas (se usa todo: la SNR crece con él) |
| Separación mínima | 2.0 mm | Entre centros de vueltas vecinas; ≥ surco + 0.4 mm |
| Región homogénea | 70 × 100 mm | Cilindro centrado (diámetro × largo) |
| Criterio de homogeneidad | RMS | RMS o pico a pico (ver abajo) |
| Largo del bobinado libre | activado | El bobinado más corto que cumpla la inhomogeneidad máxima |
| Inhomogeneidad máxima | 1 % | Límite según el criterio elegido |
| Segmentos por vuelta / del arco | 360 / 24 | Resolución de la malla (tamaño del STL) |

### Criterio de homogeneidad: RMS o pico a pico

*RMS* es la raíz del error cuadrático medio de Bz en la región (minimizar el RMS o el error
cuadrático medio da el mismo diseño). *Pico a pico* es la diferencia entre el punto más
fuerte y el más débil. Comparados sobre cuatro geometrías distintas, optimizar el pico a pico
baja el peor caso sólo un 4–8 %, pero empeora el RMS un 33–55 % y no mejora el campo mínimo:
el pico a pico lo deciden unos pocos puntos en las esquinas de la región, junto al cable,
donde manda la asimetría propia de la hélice. Por eso el valor por defecto es RMS.

### Largo del bobinado libre

Con el cable fijo, acortar el bobinado junta las vueltas y sube el campo por amperio, pero
empeora la homogeneidad. Con el largo libre se busca por bisección el bobinado más corto que
cumple la inhomogeneidad máxima; el resultado se verifica con Biot-Savart sobre la hélice
real y, si la asimetría de la hélice lo deja por encima del límite, se ajusta hasta que
cumpla. Si el límite no se alcanza ni con el largo máximo, se usa la máxima homogeneidad y
la interfaz lo avisa. Con el largo fijo se usa todo el largo disponible.

## Archivos de cada proyecto

| Archivo | Contenido |
|---|---|
| `bobina_solenoide.stl` | Soporte en mm, sólido cerrado listo para el laminador |
| `Resumen.txt` | Geometría, distribución, homogeneidad, inductancia y resistencia |
| `vueltas.csv` | Altura de cada vuelta en θ = 0 y paso local |
| `trayectoria_cable.csv` | Centro del cable en cada segmento (θ, x, y, z) |
| `parametros.json` | Parámetros y métricas en formato máquina |
| `distribucion_vueltas.png` | Densidad de vueltas y paso a lo largo del eje |
| `campo_perfiles.png` | Bz sobre el eje y en el borde de la región homogénea |
| `mapa_homogeneidad.png` | Desviación de Bz en el plano (ρ, z): uniforme vs elegida |
| `corte_stl.png` | Corte del STL para revisar el surco |

## Cómo se distribuyen las vueltas

La densidad de vueltas por unidad de largo es un polinomio par en la coordenada normalizada
u (−1 en la primera vuelta, 0 en el centro, +1 en la última):

    n(u) = 1 + c1·u² + c2·u⁴ + c3·u⁶

El número de vueltas queda fijado por el largo del cable. Los coeficientes se optimizan
(Nelder-Mead) con tres restricciones: densidad creciente desde el centro hacia los extremos,
densidad positiva y separación mínima entre vueltas. Para optimizar, el bobinado se modela
como una lámina de anillos coaxiales (campo exacto con integrales elípticas); después se
verifica con Biot-Savart sobre la hélice real. La hélice se construye invirtiendo la
fracción acumulada de vueltas F(u): θ(u) = 2π·N·F(u), z(u) = z_centro + u·L/2.

## Geometría y malla

- Perfil en U: fondo semicircular del diámetro del surco y paredes rectas hasta la
  superficie. El perfil se barre en planos meridianos y se escala por 1/cos(α) para que el
  ancho perpendicular al cable sea exactamente el diámetro del surco.
- Terminales: canal recto paralelo al eje desde cada extremo del bobinado hasta el borde.
- Las operaciones booleanas usan [manifold3d](https://github.com/elalish/manifold). Para que
  el STL no tenga triángulos degenerados al exportarse en float32, la hélice se muestrea en
  ángulos exactos, el facetado del tubo se gira lejos de las aristas del surco y, si aun así
  queda alguna astilla, la operación se repite con un giro mínimo.
- Si la impresora cierra los huecos, aumente el diámetro del surco 0.1–0.2 mm.

## Línea de comandos

Con las dependencias instaladas (`pip install -r requirements.txt`), los valores por defecto
de `Parametros` en `generar_bobina.py` se usan así:

```bash
python generar_bobina.py nombre_del_proyecto
```

`python servidor_gui.py --puerto 8765 --no-abrir` levanta la interfaz sin abrir el navegador.

## Estructura

| Archivo | Rol |
|---|---|
| `iniciar_gui.bat`, `iniciar.ps1` | Lanzador: Python, dependencias en `.venv` e interfaz |
| `servidor_gui.py`, `gui/index.html` | Interfaz web local |
| `generar_bobina.py` | Parámetros, diseño, gráficos, resumen y archivos |
| `distribucion.py` | Distribución de vueltas, hélice y optimización |
| `campo.py` | Campo de espiras, Biot-Savart e inductancia |
| `geometria.py` | Malla del soporte y exportación STL |

## Notas

- El campo y la inductancia no incluyen los terminales ni el cableado externo.
- La inductancia es una estimación (una espira coaxial por vuelta, aproximación de alta
  frecuencia para la autoinductancia de cada espira).
- La generación tarda entre 25 y 50 s según la geometría.
- Pendiente: campo B0 o frecuencia de trabajo para calcular resistencia en RF, Q y SNR;
  tipo de salida de terminales y forma de la región homogénea desde la interfaz.
