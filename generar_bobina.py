"""Antena solenoide con surco helicoidal y vueltas no uniformes.

Genera el STL del soporte (cilindro hueco con surco en U para el cable), elige la
distribución de vueltas que compensa el efecto de borde (más vueltas por centímetro en
los extremos que en el centro) y verifica el campo con Biot-Savart sobre la hélice real.

Uso:   python generar_bobina.py [nombre_proyecto]
       (o la interfaz: iniciar_gui.bat)
Salida: proyectos/<nombre_proyecto>/ (STL, trayectoria del cable, gráficos y Resumen.txt)
"""
import json
import re
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

import geometria as geo
from campo import campo_polilinea, inductancia_espiras
from distribucion import (Bobinado, Distribucion, bz_lamina, homogeneidad, optimizar,
                          optimizar_largo, puntos_roi)

CARPETA_PROYECTOS = Path(__file__).resolve().parent / "proyectos"

# Colores de los gráficos
COLOR_OPT = "#2a78d6"
COLOR_UNI = "#eb6834"
TEXTO = "#0b0b0b"
TEXTO_2 = "#52514e"
REJILLA = "#e4e3df"
FONDO = "#fcfcfb"

ARCHIVOS = [
    ("bobina_solenoide.stl", "Soporte para imprimir (mm)"),
    ("Resumen.txt", "Resumen de geometría, distribución y campo"),
    ("vueltas.csv", "Altura de cada vuelta en theta = 0 y paso local"),
    ("trayectoria_cable.csv", "Centro del cable (theta, x, y, z en mm)"),
    ("parametros.json", "Parámetros y métricas"),
    ("distribucion_vueltas.png", "Densidad de vueltas y paso"),
    ("campo_perfiles.png", "Bz sobre el eje y en el borde de la ROI"),
    ("mapa_homogeneidad.png", "Mapa de homogeneidad de Bz"),
    ("corte_stl.png", "Corte del STL"),
]


# ----------------------------------------------------------------------------------------
# PARÁMETROS
# ----------------------------------------------------------------------------------------
@dataclass
class Parametros:
    # Editables desde la interfaz (por defecto: la antena de 100/104 x 170 mm)
    diametro_interno_mm: float = 100.0
    diametro_externo_mm: float = 104.0
    largo_cilindro_mm: float = 170.0
    margen_extremos_mm: float = 5.0       # del borde del cilindro al bobinado más largo posible
    diametro_surco_mm: float = 1.0        # súmele holgura (1.1-1.2) si la impresora cierra huecos
    diametro_cable_mm: float = 1.0        # cable real (cobre); descansa en el fondo del surco
    largo_cable_m: float = 7.5            # cable alojado en el surco helicoidal (sin terminales)
    paso_minimo_mm: float = 2.0           # separación mínima entre vueltas (centro a centro)
    roi_diametro_mm: float = 70.0         # región homogénea: cilindro centrado
    roi_largo_mm: float = 100.0
    objetivo: str = "rms"                 # criterio de homogeneidad: "rms" o "pico" (pico a pico)
    largo_libre: bool = True              # bobinado más corto que cumpla la inhomogeneidad máxima
    inhomogeneidad_max_pct: float = 1.0   # según el criterio (sólo con largo libre)
    segmentos_por_vuelta: int = 360       # resolución angular de la hélice
    segmentos_arco_surco: int = 24        # resolución del fondo semicircular del surco

    # Avanzados (por ahora sólo desde código)
    profundidad_centro_cable_mm: float | None = None   # None: surco al ras de la superficie
    distribucion: str = "optimizada"      # "optimizada" o "uniforme"
    orden_polinomio: int = 3              # coeficientes c_k de n(u) = 1 + c1 u^2 + c2 u^4 + ...
    densidad_monotona: bool = True        # densidad creciente desde el centro a los extremos
    coeficientes_manuales: list | None = None
    roi_forma: str = "cilindro"           # "cilindro" o "esfera" (usa roi_diametro_mm)
    salida_terminales: str = "axial"      # "axial", "radial" o "ninguna"
    diametro_agujero_mm: float = 1.4      # sólo para "radial"
    segmentos_cilindro: int | None = None  # None: el doble de segmentos_por_vuelta
    resistividad_ohm_m: float = 1.68e-8   # cobre a 20 °C

    @property
    def r_int(self):
        return self.diametro_interno_mm / 2.0

    @property
    def r_ext(self):
        return self.diametro_externo_mm / 2.0

    @property
    def cable_mm(self):
        return self.diametro_cable_mm

    @property
    def profundidad_centro_mm(self):
        if self.profundidad_centro_cable_mm is None:
            return self.diametro_surco_mm / 2.0
        return self.profundidad_centro_cable_mm

    @property
    def r_surco_mm(self):
        """Radio del centro del semicírculo del surco."""
        return self.r_ext - self.profundidad_centro_mm

    @property
    def r_cable_mm(self):
        """Radio del centro del cable, que descansa en el fondo del surco."""
        return self.r_surco_mm - (self.diametro_surco_mm - self.diametro_cable_mm) / 2.0

    @property
    def paso_min_mm(self):
        return self.paso_minimo_mm

    @property
    def nombre_criterio(self):
        return "RMS" if self.objetivo == "rms" else "pico a pico"

    @property
    def seg_cilindro(self):
        return self.segmentos_cilindro or 2 * self.segmentos_por_vuelta

    def triangulos_estimados(self):
        """Tamaño aproximado de la malla (ajustado con proyectos generados)."""
        perfil = self.segmentos_arco_surco + 3
        return int(self.vueltas_aprox() * self.segmentos_por_vuelta * (2 * perfil + 12))

    @property
    def largo_bobinado_max_mm(self):
        return self.largo_cilindro_mm - 2.0 * self.margen_extremos_mm

    @property
    def roi_radio_mm(self):
        return self.roi_diametro_mm / 2.0

    @property
    def roi_semialto_mm(self):
        return self.roi_radio_mm if self.roi_forma == "esfera" else self.roi_largo_mm / 2.0

    def vueltas_aprox(self):
        return self.largo_cable_m * 1e3 / (2.0 * np.pi * self.r_cable_mm)

    def validar(self):
        """Lista de errores en lenguaje llano (vacía si los parámetros son válidos)."""
        errores = []
        numeros = {"diametro_interno_mm": "El diámetro interno", "diametro_externo_mm": "El diámetro externo",
                   "largo_cilindro_mm": "El largo del cilindro", "margen_extremos_mm": "El margen en los extremos",
                   "diametro_surco_mm": "El diámetro del surco", "diametro_cable_mm": "El diámetro del cable",
                   "largo_cable_m": "El largo del cable", "paso_minimo_mm": "La separación mínima",
                   "roi_diametro_mm": "El diámetro de la región homogénea",
                   "roi_largo_mm": "El largo de la región homogénea"}
        for nombre, texto in numeros.items():
            valor = getattr(self, nombre)
            if not isinstance(valor, (int, float)) or not np.isfinite(valor) or valor <= 0:
                errores.append(f"{texto} debe ser un número mayor que 0.")
        otros = []
        for nombre, texto, minimo, maximo in (("segmentos_por_vuelta", "Los segmentos por vuelta", 36, 1440),
                                              ("segmentos_arco_surco", "Los segmentos del arco del surco", 4, 64)):
            valor = getattr(self, nombre)
            if not isinstance(valor, int) or not minimo <= valor <= maximo:
                otros.append(f"{texto} deben ser un entero entre {minimo} y {maximo}.")
        malla_valida = not otros
        if self.objetivo not in ("rms", "pico"):
            otros.append("El criterio de homogeneidad debe ser 'rms' o 'pico'.")
        if errores:
            return errores + otros  # con números no válidos las comprobaciones siguientes no tienen sentido
        errores = otros
        if self.diametro_cable_mm > self.diametro_surco_mm:
            errores.append("El cable debe ser igual o más delgado que el surco.")
        if self.paso_minimo_mm < self.diametro_surco_mm + 0.4:
            errores.append(f"La separación mínima debe ser al menos el surco + 0.4 mm "
                           f"({self.diametro_surco_mm + 0.4:.2f} mm) para que quede pared entre vueltas.")
        if self.margen_extremos_mm < self.diametro_surco_mm:
            errores.append("El margen en los extremos debe ser al menos el diámetro del surco.")
        pared = self.r_ext - self.r_int
        fondo = self.profundidad_centro_mm + self.diametro_surco_mm / 2.0
        if pared <= 0:
            errores.append("El diámetro externo debe ser mayor que el interno.")
        elif pared - fondo < 0.5:
            errores.append(f"La pared ({pared:.2f} mm) es muy delgada para un surco de "
                           f"{fondo:.2f} mm de profundidad: deben quedar al menos 0.5 mm bajo el surco.")
        if self.roi_diametro_mm >= self.diametro_interno_mm:
            errores.append("El diámetro de la región homogénea debe ser menor que el diámetro interno.")
        if self.roi_forma == "cilindro" and self.roi_largo_mm > self.largo_cilindro_mm:
            errores.append("El largo de la región homogénea no puede superar el largo del cilindro.")
        if self.largo_bobinado_max_mm <= 10.0:
            errores.append("El cilindro es demasiado corto para bobinar.")
        elif self.largo_bobinado_max_mm / self.vueltas_aprox() < self.paso_min_mm:
            cabe = self.largo_bobinado_max_mm / self.paso_min_mm * 2 * np.pi * self.r_cable_mm / 1e3
            errores.append(f"El cable no cabe: con vueltas separadas {self.paso_min_mm:.1f} mm caben "
                           f"como máximo {cabe:.2f} m de cable.")
        if self.vueltas_aprox() > 400:
            errores.append("Demasiadas vueltas (más de 400): la malla sería enorme.")
        elif malla_valida and self.triangulos_estimados() > 5e6:
            errores.append(f"La malla sería demasiado grande (≈{self.triangulos_estimados() / 1e6:.1f} M "
                           f"triángulos): baje los segmentos por vuelta o del arco del surco.")
        if self.largo_libre and not (0 < self.inhomogeneidad_max_pct < 100):
            errores.append("La inhomogeneidad máxima debe estar entre 0 y 100 %.")
        return errores


def carpeta_proyecto(nombre):
    """Carpeta proyectos/<nombre>, validando que el nombre sirva como carpeta de Windows."""
    nombre = (nombre or "").strip()
    if not nombre:
        raise ValueError("Falta el nombre del proyecto.")
    if (len(nombre) > 60 or not re.fullmatch(r"[\w\- .]+", nombre) or nombre.strip(".") == ""
            or nombre.endswith((".", " "))
            or nombre.upper().split(".")[0] in {"CON", "PRN", "AUX", "NUL", *(f"COM{i}" for i in range(1, 10)),
                                                *(f"LPT{i}" for i in range(1, 10))}):
        raise ValueError("Nombre de proyecto no válido: use letras, números, espacios, '-', '_' o '.' "
                         "(máximo 60 caracteres).")
    return CARPETA_PROYECTOS / nombre


# ----------------------------------------------------------------------------------------
# Diseño
# ----------------------------------------------------------------------------------------
def disenar(p, log=print):
    """Bobinado uniforme de referencia (largo máximo) y bobinado elegido."""
    radio = p.r_cable_mm * 1e-3
    largo_max = p.largo_bobinado_max_mm * 1e-3
    z_centro = p.largo_cilindro_mm / 2.0 * 1e-3
    roi = puntos_roi(p.roi_forma, p.roi_radio_mm * 1e-3, p.roi_semialto_mm * 1e-3, z_centro)
    paso_min = p.paso_min_mm * 1e-3
    info = {"largo_libre": False, "cumple_objetivo": None}

    uniforme = Bobinado(Distribucion([]), radio, largo_max, z_centro, p.largo_cable_m)
    largo = largo_max
    if p.distribucion == "uniforme":
        coeficientes = []
    elif p.coeficientes_manuales is not None:
        coeficientes = list(p.coeficientes_manuales)
    elif p.distribucion == "optimizada" and p.largo_libre:
        log(f"Buscando el bobinado más corto con inhomogeneidad {p.nombre_criterio} "
            f"<= {p.inhomogeneidad_max_pct} %...")
        objetivo = p.inhomogeneidad_max_pct / 100.0
        meta = objetivo
        for _ in range(4):
            # La optimización usa el modelo axisimétrico; la hélice real suma su propia
            # asimetría, así que se verifica sobre ella y se ajusta la meta si se pasa.
            largo, coeficientes, _, cumple = optimizar_largo(
                radio, largo_max, z_centro, p.largo_cable_m, roi, p.orden_polinomio, paso_min,
                meta, monotona=p.densidad_monotona, objetivo=p.objetivo)
            real = _inhomogeneidad_helice(p, Bobinado(Distribucion(coeficientes), radio, largo,
                                                      z_centro, p.largo_cable_m))
            log(f"  Bobinado de {largo * 1e3:.1f} mm: {real * 100:.3f} % sobre la hélice real")
            if not cumple or real <= objetivo * 1.005:
                break
            meta *= 0.98 * objetivo / real
        cumple = real <= objetivo * 1.005
        info.update(largo_libre=True, cumple_objetivo=bool(cumple))
        if not cumple:
            log("  El objetivo no se alcanza ni con el largo máximo: se usa la máxima homogeneidad.")
    elif p.distribucion == "optimizada":
        log("Optimizando la distribución con el largo máximo...")
        coeficientes, _ = optimizar(radio, largo_max, z_centro, p.largo_cable_m, roi,
                                    p.orden_polinomio, paso_min, monotona=p.densidad_monotona,
                                    objetivo=p.objetivo)
    else:
        raise ValueError(f"distribución desconocida: {p.distribucion!r}")
    elegido = Bobinado(Distribucion(coeficientes), radio, largo, z_centro, p.largo_cable_m)

    paso = elegido.paso(np.linspace(-1, 1, 2001)).min() * 1e3
    if paso < max(p.paso_min_mm, p.diametro_surco_mm) - 1e-6:
        raise ValueError(f"vueltas demasiado juntas: paso mínimo {paso:.2f} mm")
    info["largo_bobinado_mm"] = largo * 1e3
    log(f"  Bobinado de {largo * 1e3:.1f} mm, {elegido.vueltas:.2f} vueltas, "
        f"n(u) = {np.round(elegido.distribucion.coeficientes, 3).tolist()}")
    return uniforme, elegido, roi, info


def _inhomogeneidad_helice(p, bobinado):
    """Inhomogeneidad (según p.objetivo) de Bz en la ROI calculada sobre la hélice real."""
    puntos, pesos = puntos_roi_3d(p, bobinado)
    bz = campo_polilinea(bobinado.puntos_cable(360), puntos)[:, 2]
    return homogeneidad(bz, pesos)[2 if p.objetivo == "pico" else 1]


def puntos_roi_3d(p, bobinado, n_rho=11, n_phi=16, n_z=21):
    """Puntos cartesianos de la ROI para la verificación con la hélice real."""
    radio = p.roi_radio_mm * 1e-3
    rho = np.linspace(0.0, radio, n_rho)
    phi = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)
    z = np.linspace(-1.0, 1.0, n_z) * p.roi_semialto_mm * 1e-3
    r, f, zz = np.meshgrid(rho, phi, z, indexing="ij")
    if p.roi_forma == "esfera":
        dentro = r**2 + zz**2 <= radio**2 * (1 + 1e-9)
    else:
        dentro = np.ones_like(r, dtype=bool)
    r, f, zz = r[dentro], f[dentro], zz[dentro]
    pesos = np.where(r > 0, r, (rho[1] - rho[0]) / 8.0)
    puntos = np.column_stack([r * np.cos(f), r * np.sin(f), zz + bobinado.z_centro])
    return puntos, pesos


def metricas(p, bobinado, roi, puntos_3d, pesos_3d):
    rho, z, pesos = roi
    media, rms, pp = homogeneidad(bz_lamina(bobinado, rho, z), pesos)
    cable = bobinado.puntos_cable(360)
    b = campo_polilinea(cable, puntos_3d)
    media_h, rms_h, pp_h = homogeneidad(b[:, 2], pesos_3d)
    centro = campo_polilinea(cable, np.array([[0.0, 0.0, bobinado.z_centro]]))[0]
    z_v, frac = bobinado.vueltas_individuales()
    inductancia = inductancia_espiras(bobinado.radio, z_v, frac, p.cable_mm * 0.5e-3)
    return {
        "vueltas": bobinado.vueltas,
        "largo_bobinado_mm": bobinado.largo * 1e3,
        "largo_cable_m": bobinado.largo_cable(bobinado.vueltas),
        "paso_extremo_mm": float(bobinado.paso(1.0) * 1e3),
        "paso_centro_mm": float(bobinado.paso(0.0) * 1e3),
        "bz_centro_uT_por_A": float(centro[2] * 1e6),
        "lamina": {"bz_medio_uT_por_A": media * 1e6, "rms_pct": rms * 100, "pico_pico_pct": pp * 100},
        "helice": {"bz_medio_uT_por_A": media_h * 1e6, "bz_min_uT_por_A": float(b[:, 2].min() * 1e6),
                   "rms_pct": rms_h * 100, "pico_pico_pct": pp_h * 100,
                   "b_transversal_max_pct": float(np.hypot(b[:, 0], b[:, 1]).max() / media_h * 100)},
        "inductancia_uH": inductancia * 1e6,
    }


def construir_solido(p, bobinado):
    theta, z = bobinado.trayectoria(p.segmentos_por_vuelta)
    z_mm = z * 1e3
    r_c = p.r_surco_mm  # el surco se talla en su radio; el cable apoya en su fondo
    d = p.diametro_surco_mm
    cortadores = [geo.cortador_helicoidal(theta, z_mm, r_c, d, p.r_ext, p.segmentos_arco_surco)]
    # Ángulos donde los cortadores tienen aristas que atraviesan las paredes del tubo
    angulos = [np.arange(p.segmentos_por_vuelta) * 2.0 * np.pi / p.segmentos_por_vuelta, [theta[-1]]]
    if p.salida_terminales == "axial":
        # El canal sobrepasa la pared del surco en 0.2 d: codo despejado y sin caras coplanares
        a = 0.7 * d
        cortadores.append(geo.canal_axial(theta[0], -1.0, z_mm[0] + a, r_c, d, p.r_ext,
                                          p.segmentos_arco_surco))
        cortadores.append(geo.canal_axial(theta[-1], z_mm[-1] - a, p.largo_cilindro_mm + 1.0, r_c, d,
                                          p.r_ext, p.segmentos_arco_surco))
        w = np.arcsin(d / 2.0 / p.r_ext)
        angulos.append([theta[0] - w, theta[0] + w, theta[-1] - w, theta[-1] + w])
    elif p.salida_terminales == "radial":
        # El agujero queda entero dentro del surco, con su borde un cuarto de surco antes del
        # extremo: si el borde tocara justo la tapa del surco, la booleana dejaría astillas
        delta = (p.diametro_agujero_mm / 2.0 + 0.25 * d) / r_c
        w = np.arcsin(p.diametro_agujero_mm / 2.0 / p.r_int)
        # El agujero termina en la parte recta de las paredes del surco: si terminara justo
        # en r_c coincidiría con vértices del perfil y la booleana dejaría astillas
        r_fin = r_c + 0.15 * d
        for t in (theta[0] + delta, theta[-1] - delta):
            z_t = bobinado.z_de_u(bobinado.u_de_theta(t)) * 1e3
            cortadores.append(geo.agujero_radial(t, z_t, p.r_int - 1.0, r_fin, p.diametro_agujero_mm))
            angulos.append([t - w, t, t + w])
    elif p.salida_terminales != "ninguna":
        raise ValueError(f"salida de terminales desconocida: {p.salida_terminales!r}")
    giro = geo.giro_facetado(np.concatenate(angulos), p.seg_cilindro)
    paso = 2.0 * np.pi / p.segmentos_por_vuelta
    for intento in range(6):
        # Si la booleana deja astillas (una intersección pasó casi por un vértice), se giran
        # levemente las salidas y el facetado del tubo (centésimas de mm) y se repite.
        tubo = geo.tubo(p.r_int, p.r_ext, p.largo_cilindro_mm, p.seg_cilindro, giro + intento * 0.031 * paso)
        salidas = [c.rotate([0.0, 0.0, np.degrees(intento * 0.037 * paso)]) for c in cortadores[1:]]
        solido = geo.mf.Manifold.batch_boolean([tubo, cortadores[0]] + salidas, geo.mf.OpType.Subtract)
        if solido.status() != geo.mf.Error.NoError:
            raise RuntimeError(f"fallo en la operación booleana: {solido.status()}")
        if geo.astillas(solido) == 0:
            break
    return solido


def guardar_tablas(p, bobinado, carpeta):
    theta, z = bobinado.trayectoria(p.segmentos_por_vuelta)
    xyz = np.column_stack([bobinado.radio * np.cos(theta), bobinado.radio * np.sin(theta), z]) * 1e3
    np.savetxt(carpeta / "trayectoria_cable.csv", np.column_stack([np.degrees(theta), xyz]),
               delimiter=",", fmt="%.5f", header="theta_grados,x_mm,y_mm,z_mm", comments="")

    k, z_k = bobinado.cruces()
    u_k = (z_k - bobinado.z_centro) / (bobinado.largo / 2.0)
    paso = bobinado.paso(u_k) * 1e3
    filas = np.column_stack([k, z_k * 1e3, (z_k - bobinado.z_centro) * 1e3, paso,
                             np.append(np.diff(z_k) * 1e3, np.nan)])
    np.savetxt(carpeta / "vueltas.csv", filas, delimiter=",", fmt="%.4f", comments="",
               header="vuelta,z_mm,z_desde_centro_mm,paso_local_mm,distancia_a_siguiente_mm")


# ----------------------------------------------------------------------------------------
# Gráficos (API orientada a objetos: segura fuera del hilo principal)
# ----------------------------------------------------------------------------------------
def _figura(**kwargs):
    import matplotlib
    from matplotlib.figure import Figure

    matplotlib.rcParams.update({
        "figure.facecolor": FONDO, "axes.facecolor": FONDO, "savefig.facecolor": FONDO,
        "axes.edgecolor": TEXTO_2, "axes.labelcolor": TEXTO, "text.color": TEXTO,
        "xtick.color": TEXTO_2, "ytick.color": TEXTO_2, "axes.grid": True,
        "grid.color": REJILLA, "grid.linewidth": 0.8, "axes.spines.top": False,
        "axes.spines.right": False, "font.size": 10, "axes.titlesize": 11,
        "axes.titleweight": "bold", "legend.frameon": False, "lines.linewidth": 2,
    })
    return Figure(layout="constrained", **kwargs)


def _sombrear_roi(ax, p):
    zc = p.largo_cilindro_mm / 2.0
    ax.axvspan(zc - p.roi_semialto_mm, zc + p.roi_semialto_mm, color=REJILLA, alpha=0.6, lw=0)


def grafico_distribucion(p, uniforme, elegido, carpeta):
    fig = _figura(figsize=(8, 6.5))
    ax1, ax2 = fig.subplots(2, 1, sharex=True)
    u = np.linspace(-1, 1, 801)
    for b, color, nombre in ((uniforme, COLOR_UNI, "Uniforme (largo máximo)"),
                             (elegido, COLOR_OPT, "Elegida")):
        z = b.z_de_u(u) * 1e3
        densidad = b.vueltas * b.distribucion.densidad_normalizada(u) / (b.largo * 50.0)  # vueltas/cm
        ax1.plot(z, densidad, color=color, label=nombre)
        ax2.plot(z, b.paso(u) * 1e3, color=color, label=nombre)
    _, z_k = elegido.cruces()
    for zk in z_k * 1e3:
        ax1.axvline(zk, color=COLOR_OPT, alpha=0.25, linewidth=1)
    for ax in (ax1, ax2):
        _sombrear_roi(ax, p)
    ax1.set_ylabel("Densidad de vueltas [vueltas/cm]")
    fig.suptitle(f"Distribución de vueltas ({elegido.vueltas:.2f} vueltas, {p.largo_cable_m:g} m de cable, "
                 f"bobinado de {elegido.largo * 1e3:.1f} mm)", fontweight="bold")
    ax1.set_title("líneas verticales: paso del cable por θ = 0 · sombra: región homogénea",
                  fontsize=9, fontweight="normal", color=TEXTO_2)
    ax1.legend(loc="upper center", ncol=2)
    paso_max = max(uniforme.paso(0.0), elegido.paso(np.linspace(-1, 1, 201)).max()) * 1e3
    ax2.set_ylim(0, paso_max * 1.3)
    ax2.set_ylabel("Paso local [mm]")
    ax2.set_xlabel("Altura z sobre la base del cilindro [mm]")
    ax2.axhline(p.paso_min_mm, color=TEXTO_2, linestyle="--", linewidth=1)
    ax2.text(p.largo_cilindro_mm / 2, p.paso_min_mm + paso_max * 0.03,
             f"paso mínimo ({p.paso_min_mm:g} mm)", ha="center", color=TEXTO_2, fontsize=9)
    ax2.set_xlim(0, p.largo_cilindro_mm)
    fig.savefig(carpeta / "distribucion_vueltas.png", dpi=150)


def grafico_campo(p, uniforme, elegido, carpeta):
    z = np.linspace(0.0, p.largo_cilindro_mm, 341) * 1e-3
    radios = [0.0, p.roi_radio_mm * 1e-3]
    fig = _figura(figsize=(11, 4.4))
    ejes = fig.subplots(1, 2, sharey=True)
    for ax, r in zip(ejes, radios):
        for b, color, nombre in ((uniforme, COLOR_UNI, "Uniforme (largo máximo)"),
                                 (elegido, COLOR_OPT, "Elegida")):
            cable = b.puntos_cable(360)
            puntos = np.column_stack([np.full_like(z, r), np.zeros_like(z), z])
            bz = campo_polilinea(cable, puntos)[:, 2]
            b0 = campo_polilinea(cable, np.array([[0.0, 0.0, b.z_centro]]))[0, 2]
            ax.plot(z * 1e3, (bz / b0 - 1.0) * 100, color=color, label=nombre)
        _sombrear_roi(ax, p)
        ax.set_title("Sobre el eje (ρ = 0)" if r == 0 else f"A ρ = {r * 1e3:.0f} mm del eje (borde de la ROI)")
        ax.set_xlabel("Altura z [mm]")
        ax.set_xlim(0, p.largo_cilindro_mm)
    ejes[0].set_ylabel("Bz respecto del centro [%]")
    ejes[0].legend(loc="lower center")
    ejes[1].text(p.largo_cilindro_mm / 2, 0.96, "ROI", ha="center", va="top", color=TEXTO_2, fontsize=9,
                 transform=ejes[1].get_xaxis_transform())
    fig.suptitle("Campo axial Bz (Biot-Savart sobre la hélice real)", fontweight="bold")
    fig.savefig(carpeta / "campo_perfiles.png", dpi=150)


def grafico_mapa(p, uniforme, elegido, carpeta):
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

    rho_max = 0.9 * p.r_int
    rho = np.linspace(0.0, rho_max, 91) * 1e-3
    z = np.linspace(0.0, p.largo_cilindro_mm, 171) * 1e-3
    r, zz = np.meshgrid(rho, z)
    mapa = LinearSegmentedColormap.from_list("div", [COLOR_OPT, "#f0efec", COLOR_UNI])
    niveles_lineas = [-5, -2, -1, -0.5, 0.5, 1, 2, 5]
    fig = _figura(figsize=(10, 5.6))
    ejes = fig.subplots(1, 2, sharey=True)
    zc = p.largo_cilindro_mm / 2
    alto = p.roi_semialto_mm
    for ax, b, nombre in ((ejes[0], uniforme, "Uniforme (largo máximo)"), (ejes[1], elegido, "Elegida")):
        bz = bz_lamina(b, r.ravel(), zz.ravel()).reshape(r.shape)
        b0 = bz_lamina(b, np.array([0.0]), np.array([b.z_centro]))[0]
        desv = np.clip((bz / b0 - 1.0) * 100, -10, 10)
        malla = ax.contourf(r * 1e3, zz * 1e3, desv, levels=np.linspace(-10, 10, 41), cmap=mapa,
                            norm=TwoSlopeNorm(0, -10, 10))
        cs = ax.contour(r * 1e3, zz * 1e3, desv, levels=niveles_lineas, colors=TEXTO_2, linewidths=0.7)
        ax.clabel(cs, fmt="%g %%", fontsize=7)
        if p.roi_forma == "esfera":
            t = np.linspace(-np.pi / 2, np.pi / 2, 100)
            ax.plot(p.roi_radio_mm * np.cos(t), zc + p.roi_radio_mm * np.sin(t), color=TEXTO, lw=1.2)
        else:
            ax.plot([0, p.roi_radio_mm, p.roi_radio_mm, 0], [zc - alto, zc - alto, zc + alto, zc + alto],
                    color=TEXTO, lw=1.2)
        _, z_k = b.cruces()
        ax.plot(np.full_like(z_k, 0.97 * rho_max), z_k * 1e3, ">", color=TEXTO, ms=5)
        ax.set_title(nombre)
        ax.set_xlabel("ρ [mm]")
        ax.set_xlim(0, rho_max)
        ax.grid(False)
    ejes[0].set_ylabel("Altura z [mm]")
    barra = fig.colorbar(malla, ax=ejes, shrink=0.85, ticks=[-10, -5, 0, 5, 10])
    barra.set_label("Bz respecto del centro [%] (saturado en ±10 %)")
    fig.suptitle("Homogeneidad de Bz (recuadro: ROI; ▸ altura de cada vuelta en θ = 0)",
                 fontweight="bold", fontsize=10)
    fig.savefig(carpeta / "mapa_homogeneidad.png", dpi=150)


def grafico_corte(p, solido, carpeta):
    poligonos = geo.seccion_meridiana(solido)
    fig = _figura(figsize=(10, 6))
    ax1, ax2 = fig.subplots(1, 2, gridspec_kw={"width_ratios": [1.2, 1]})
    for ax in (ax1, ax2):
        for pol in poligonos:
            ax.fill(pol[:, 0], pol[:, 1], color="#c3c2b7", edgecolor=TEXTO, linewidth=0.6)
        ax.set_aspect("equal")
        ax.grid(False)
    ax1.set_title("Corte por el plano y = 0")
    ax1.set_xlabel("x [mm]")
    ax1.set_ylabel("z [mm]")
    ax2.set_xlim(p.r_int - 1.5, p.r_ext + 1.5)
    ax2.set_ylim(-1, min(30.0, p.largo_cilindro_mm))
    ax2.set_title("Detalle de la pared (x > 0, abajo):\nsurco en U y salida")
    ax2.set_xlabel("x [mm]")
    fig.savefig(carpeta / "corte_stl.png", dpi=150)


# ----------------------------------------------------------------------------------------
def escribir_resumen(p, carpeta, info, elegido, m_uni, m_opt, info_stl):
    c = elegido.distribucion.coeficientes
    terminos = " ".join(f"{ck:+.4f} u^{2 * k + 2}" for k, ck in enumerate(c)) or "(uniforme)"
    resistencia = p.resistividad_ohm_m * p.largo_cable_m / (np.pi * (p.cable_mm * 0.5e-3) ** 2)
    roi = (f"esfera de {p.roi_diametro_mm:g} mm de diámetro" if p.roi_forma == "esfera" else
           f"cilindro de {p.roi_diametro_mm:g} mm de diámetro x {p.roi_largo_mm:g} mm de largo")
    z_ini = (p.largo_cilindro_mm - elegido.largo * 1e3) / 2.0
    terminales = {"axial": f"canal axial de {z_ini:.1f} mm hasta cada borde",
                  "radial": f"agujero radial de {p.diametro_agujero_mm} mm al interior en cada extremo",
                  "ninguna": "sin salida"}[p.salida_terminales]
    if info["largo_libre"]:
        estado = "cumple" if info["cumple_objetivo"] else "NO se alcanza; se usa el largo máximo"
        largo_txt = (f"libre: el más corto con inhomogeneidad {p.nombre_criterio} <= "
                     f"{p.inhomogeneidad_max_pct:g} % ({estado})")
    else:
        largo_txt = "largo máximo disponible"

    def fila(nombre, clave, sub=None, fmt="{:.3f}"):
        a = m_uni[clave] if sub is None else m_uni[sub][clave]
        b = m_opt[clave] if sub is None else m_opt[sub][clave]
        return f"  {nombre:<40}{fmt.format(a):>12}{fmt.format(b):>12}"

    lineas = [
        "ANTENA SOLENOIDE CON VUELTAS NO UNIFORMES",
        f"Proyecto: {carpeta.name}",
        f"Generado: {time.strftime('%Y-%m-%d %H:%M')}",
        "",
        "SOPORTE (STL en mm)",
        f"  Diámetro interno / externo: {p.diametro_interno_mm:g} / {p.diametro_externo_mm:g} mm "
        f"(pared {p.r_ext - p.r_int:.2f} mm)",
        f"  Largo del cilindro: {p.largo_cilindro_mm:g} mm",
        f"  Surco en U: {p.diametro_surco_mm:g} mm de diámetro, fondo a "
        f"{p.profundidad_centro_mm + p.diametro_surco_mm / 2:.2f} mm de la superficie",
        f"  Cable de {p.diametro_cable_mm:g} mm apoyado en el fondo; radio de su centro: {elegido.radio * 1e3:.2f} mm",
        f"  Separación mínima entre vueltas: {p.paso_minimo_mm:g} mm; margen en los extremos: {p.margen_extremos_mm:g} mm",
        f"  Bobinado: {elegido.largo * 1e3:.1f} mm, entre z = {z_ini:.1f} y "
        f"{p.largo_cilindro_mm - z_ini:.1f} mm ({largo_txt})",
        f"  Terminales: {terminales}",
        f"  Triángulos: {info_stl['triangulos']}, volumen: {info_stl['volumen_cm3']:.2f} cm3, "
        f"género: {info_stl['genero']}, sólido válido: {info_stl['valido']}",
        f"  Malla: {p.segmentos_por_vuelta} segmentos por vuelta, {p.segmentos_arco_surco} en el arco del surco, "
        f"{p.seg_cilindro} en el cilindro",
        "",
        "DISTRIBUCIÓN DE VUELTAS",
        f"  Tipo: {p.distribucion}" + (" (coeficientes manuales)" if p.coeficientes_manuales else ""),
        f"  n(u) = 1 {terminos}",
        "  u = -1 en la primera vuelta (abajo), 0 en el centro, +1 en la última (arriba)",
        f"  Densidad extremo / centro: {elegido.distribucion.densidad(1.0):.2f}",
        f"  Región homogénea (ROI): {roi}",
        f"  Criterio de homogeneidad: {p.nombre_criterio}",
        "",
        f"  {'':<40}{'Uniforme*':>12}{'Elegida':>12}",
        fila("Largo del bobinado [mm]", "largo_bobinado_mm", fmt="{:.1f}"),
        fila("Vueltas", "vueltas"),
        fila("Largo de cable helicoidal [m]", "largo_cable_m", fmt="{:.4f}"),
        fila("Paso en el centro [mm]", "paso_centro_mm", fmt="{:.2f}"),
        fila("Paso en los extremos [mm]", "paso_extremo_mm", fmt="{:.2f}"),
        fila("Bz en el centro [uT/A]", "bz_centro_uT_por_A", fmt="{:.2f}"),
        fila("Bz medio en la ROI [uT/A]", "bz_medio_uT_por_A", "helice", fmt="{:.2f}"),
        fila("Bz mínimo en la ROI [uT/A]", "bz_min_uT_por_A", "helice", fmt="{:.2f}"),
        fila("Inhomogeneidad RMS en ROI [%]", "rms_pct", "helice"),
        fila("Inhomogeneidad pico a pico en ROI [%]", "pico_pico_pct", "helice"),
        fila("|B perpendicular al eje| máx / Bz [%]", "b_transversal_max_pct", "helice"),
        fila("RMS modelo axisimétrico (lámina) [%]", "rms_pct", "lamina"),
        fila("Inductancia estimada [uH]", "inductancia_uH", fmt="{:.2f}"),
        "  * Uniforme: el mismo cable bobinado parejo en todo el largo máximo, como referencia.",
        "",
        f"  Resistencia DC del cable ({p.largo_cable_m:g} m, cobre, {p.cable_mm:g} mm): "
        f"{resistencia * 1e3:.1f} mOhm",
        "  La inductancia es una estimación (una espira coaxial por vuelta, sin terminales).",
        "  Campo: Biot-Savart sobre la hélice real (sin terminales). La diferencia con el modelo",
        "  axisimétrico viene de la vuelta incompleta y de la componente axial de la hélice.",
        "",
        "ARCHIVOS",
        *(f"  {nombre:<26}{descripcion}" for nombre, descripcion in ARCHIVOS),
    ]
    texto = "\n".join(lineas) + "\n"
    (carpeta / "Resumen.txt").write_text(texto, encoding="utf-8")
    return texto


def generar(p, carpeta, log=print):
    """Diseña la bobina y escribe todos los archivos en `carpeta`. Devuelve un resumen."""
    errores = p.validar()
    if errores:
        raise ValueError("\n".join(errores))
    t0 = time.time()
    carpeta = Path(carpeta)
    carpeta.mkdir(parents=True, exist_ok=True)

    log("Calculando la distribución de vueltas...")
    uniforme, elegido, roi, info = disenar(p, log)

    log("Verificando el campo con Biot-Savart sobre la hélice...")
    puntos_3d, pesos_3d = puntos_roi_3d(p, elegido)
    m_uni = metricas(p, uniforme, roi, puntos_3d, pesos_3d)
    m_opt = metricas(p, elegido, roi, puntos_3d, pesos_3d)

    log("Construyendo la malla y las operaciones booleanas...")
    solido = construir_solido(p, elegido)
    triangulos = geo.guardar_stl(solido, carpeta / "bobina_solenoide.stl",
                                 "antena solenoide - surco helicoidal")
    info_stl = {"triangulos": triangulos, "volumen_cm3": solido.volume() / 1000.0,
                "genero": solido.genus(), "valido": solido.status() == geo.mf.Error.NoError}

    log("Guardando tablas y gráficos...")
    guardar_tablas(p, elegido, carpeta)
    grafico_distribucion(p, uniforme, elegido, carpeta)
    grafico_campo(p, uniforme, elegido, carpeta)
    grafico_mapa(p, uniforme, elegido, carpeta)
    grafico_corte(p, solido, carpeta)

    resistencia = p.resistividad_ohm_m * p.largo_cable_m / (np.pi * (p.cable_mm * 0.5e-3) ** 2)
    info["resistencia_dc_mohm"] = resistencia * 1e3
    resultado = {"proyecto": carpeta.name, "parametros": asdict(p), "diseno": info,
                 "coeficientes": elegido.distribucion.coeficientes.tolist(),
                 "uniforme": m_uni, "elegida": m_opt, "stl": info_stl,
                 "segundos": round(time.time() - t0, 1)}
    with open(carpeta / "parametros.json", "w", encoding="utf-8") as f:
        json.dump(resultado, f, indent=2, ensure_ascii=False)
    resultado["resumen"] = escribir_resumen(p, carpeta, info, elegido, m_uni, m_opt, info_stl)
    log(f"Listo en {resultado['segundos']} s -> {carpeta}")
    return resultado


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    nombre = sys.argv[1] if len(sys.argv) > 1 else "por_defecto"
    resultado = generar(Parametros(), carpeta_proyecto(nombre))
    print()
    print(resultado["resumen"])


if __name__ == "__main__":
    main()
