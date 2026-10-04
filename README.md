# Relativistic raytracer

Simulador en tiempo real de un agujero negro de Kerr con disco de acreción, y un prototipo en CPU.
Unidades G = c = M = 1: las distancias están en masas del agujero.

## Física del simulador

| Pieza | Modelo |
| --- | --- |
| Luz | Geodésicas nulas de Kerr, hamiltoniano en coordenadas Kerr-Schild salientes, RK4 trazado hacia el pasado |
| Cámara | Observador ZAMO (estático si el spin es 0) con base ortonormal local |
| Disco | Fino y opaco hasta 25M (transparente hasta 50M), desde la ISCO de Bardeen-Press-Teukolsky |
| Emisión | Flujo de Page-Thorne, cuerpo negro local T = Tmax · (F/Fmax)^1/4 |
| Color y brillo | Cuerpo negro observado a g·T (factor de corrimiento exacto), integrado con las funciones CIE 1931 |
| Cielo | Estrellas de cuerpo negro desplazadas al azul en la posición de la cámara |
| Animación | Turbulencia arrastrada a la velocidad orbital, evaluada en el instante de emisión (retraso de la luz) |

La física está escrita dos veces: en GLSL (`simulator/blackhole.fs`) y en C++ doble precisión
(`common/KerrPhysics.h`). Las pruebas comprueban la versión C++ contra resultados analíticos, y
`tools/validate.py` comprueba que el shader da lo mismo píxel a píxel.

## Compilar

```
cmake -S . -B build
cmake --build build
```

En Windows con MinGW usa el raylib incluido en `simulator/lib`. En Linux o macOS busca raylib
instalado y, si no está, lo descarga. El shader se copia junto al ejecutable `build/simulator`.

## Controles

| Tecla | Acción |
| --- | --- |
| WASD + ratón | Mover la cámara |
| R | Spin activado / desactivado (vuelve al último valor, 0,9 por defecto) |
| Izquierda / Derecha | Spin −0,05 / +0,05 (entre −0,95 y 0,95; negativo = disco retrógrado) |
| 1 / 2 | Temperatura máxima del disco ÷1,25 / ×1,25 |
| Arriba / Abajo | Exposición |
| V | Modo de vista: imagen, mapa de g, mapa de temperatura |
| T | Turbulencia del disco |
| RePág / AvPág | Velocidad del tiempo (M por segundo) |

## Comprobar

```
./build/physics_tests                                   # física contra valores analíticos
xvfb-run -a python3 tools/validate.py --build build     # shader contra valores exactos y la referencia en CPU
```

`tools/validate.py` y `tools/render_shader.py` necesitan `pip install moderngl numpy pillow`;
`xvfb-run` solo hace falta en una máquina sin pantalla.
