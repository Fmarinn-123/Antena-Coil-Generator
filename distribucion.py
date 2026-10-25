"""Distribución no uniforme de las vueltas a lo largo del eje y su optimización.

La densidad de vueltas (vueltas por unidad de largo) se describe con un polinomio par

    n(u) = 1 + c1 u^2 + c2 u^4 + ...,   u = 2 (z - z_centro) / L_bobinado  en [-1, 1]

de modo que n(0) = 1 en el centro y los coeficientes positivos concentran vueltas en los
extremos. Se buscan los coeficientes que hacen más uniforme el campo axial B_z dentro de
la región de interés (ROI), modelando el bobinado como una lámina de corriente formada por
muchos anillos coaxiales (cálculo exacto con integrales elípticas).
"""
import numpy as np
from scipy.optimize import brentq, minimize

from campo import campo_espira


class Distribucion:
    """Densidad relativa de vueltas n(u) = 1 + sum_k c_k u^(2k)."""

    def __init__(self, coeficientes):
        self.coeficientes = np.atleast_1d(np.asarray(coeficientes, dtype=float))

    def densidad(self, u):
        """Densidad relativa al centro (n(0) = 1)."""
        u = np.asarray(u, dtype=float)
        return np.ones_like(u) + sum(c * u ** (2 * k + 2) for k, c in enumerate(self.coeficientes))

    def derivada(self, u):
        """dn/du."""
        u = np.asarray(u, dtype=float)
        return np.zeros_like(u) + sum(c * (2 * k + 2) * u ** (2 * k + 1) for k, c in enumerate(self.coeficientes))

    def _primitiva(self, u):
        u = np.asarray(u, dtype=float)
        return u + sum(c * u ** (2 * k + 3) / (2 * k + 3) for k, c in enumerate(self.coeficientes))

    def densidad_normalizada(self, u):
        """dF/du: densidad que integra 1 en [-1, 1]."""
        return self.densidad(u) / (2.0 * self._primitiva(1.0))

    def acumulada(self, u):
        """F(u): fracción de las vueltas totales comprendida entre u = -1 y u."""
        g1 = self._primitiva(1.0)
        return (self._primitiva(u) + g1) / (2.0 * g1)


class Bobinado:
    """Hélice de paso variable sobre un cilindro, con largo de cable fijo.

    radio:      radio del centro del cable [m]
    largo:      largo axial ocupado por el bobinado [m]
    z_centro:   altura del centro del bobinado [m]
    largo_cable: largo del cable helicoidal [m]; fija el número de vueltas
    """

    def __init__(self, distribucion, radio, largo, z_centro, largo_cable):
        self.distribucion = distribucion
        self.radio = radio
        self.largo = largo
        self.z_centro = z_centro
        self.largo_cable_objetivo = largo_cable
        self.vueltas = self._resolver_vueltas(largo_cable)
        self._u_tabla = np.linspace(-1.0, 1.0, 20001)
        self._f_tabla = distribucion.acumulada(self._u_tabla)

    def largo_cable(self, vueltas):
        """Largo exacto de la hélice: integral de sqrt((r dtheta)^2 + dz^2)."""
        u = np.linspace(-1.0, 1.0, 8001)
        r_dtheta = self.radio * 2.0 * np.pi * vueltas * self.distribucion.densidad_normalizada(u)
        return float(np.trapezoid(np.hypot(r_dtheta, self.largo / 2.0), u))

    def _resolver_vueltas(self, largo_cable):
        n0 = largo_cable / (2.0 * np.pi * self.radio)
        return brentq(lambda n: self.largo_cable(n) - largo_cable, 0.5 * n0, n0)

    def paso(self, u):
        """Paso local de la hélice [m] (avance axial por vuelta)."""
        return (self.largo / 2.0) / (self.vueltas * self.distribucion.densidad_normalizada(u))

    def z_de_u(self, u):
        return self.z_centro + np.asarray(u) * self.largo / 2.0

    def u_de_theta(self, theta):
        return np.interp(np.asarray(theta) / (2.0 * np.pi * self.vueltas), self._f_tabla, self._u_tabla)

    def trayectoria(self, puntos_por_vuelta):
        """Ángulo theta [rad] y altura z [m].

        theta toma múltiplos exactos de 2*pi/puntos_por_vuelta (más el punto final), así los
        vértices de la malla caen en los mismos ángulos en todas las vueltas.
        """
        paso = 2.0 * np.pi / puntos_por_vuelta
        theta_fin = 2.0 * np.pi * self.vueltas
        theta = np.append(np.arange(0.0, theta_fin - 0.25 * paso, paso), theta_fin)
        return theta, self.z_de_u(self.u_de_theta(theta))

    def puntos_cable(self, puntos_por_vuelta):
        theta, z = self.trayectoria(puntos_por_vuelta)
        return np.column_stack([self.radio * np.cos(theta), self.radio * np.sin(theta), z])

    def anillos(self, n=600):
        """Modelo de lámina: alturas de n anillos y sus vueltas (fraccionarias)."""
        bordes = np.linspace(-1.0, 1.0, n + 1)
        pesos = self.vueltas * np.diff(self.distribucion.acumulada(bordes))
        return self.z_de_u(0.5 * (bordes[1:] + bordes[:-1])), pesos

    def vueltas_individuales(self, muestras=720):
        """Altura media y fracción de cada vuelta (la última puede ser incompleta)."""
        z_media, fraccion = [], []
        k = 0
        while k < self.vueltas:
            fin = min(k + 1.0, self.vueltas)
            theta = 2.0 * np.pi * np.linspace(k, fin, muestras)
            z_media.append(float(np.mean(self.z_de_u(self.u_de_theta(theta)))))
            fraccion.append(fin - k)
            k += 1
        return np.array(z_media), np.array(fraccion)

    def cruces(self):
        """Altura donde el cable pasa por theta = 0 (cada vuelta completa) y por el final."""
        k = np.arange(0, int(np.floor(self.vueltas)) + 1, dtype=float)
        k = np.append(k, self.vueltas) if self.vueltas - k[-1] > 1e-9 else k
        return k, self.z_de_u(self.u_de_theta(2.0 * np.pi * k))


def puntos_roi(forma, radio, semilargo, z_centro, n_rho=21, n_z=41):
    """Puntos (rho, z) de la región de interés y pesos de volumen (simetría axial)."""
    rho = np.linspace(0.0, radio, n_rho)
    alto = radio if forma == "esfera" else semilargo
    z = np.linspace(-alto, alto, n_z)
    r, zz = np.meshgrid(rho, z)
    if forma == "esfera":
        dentro = r**2 + zz**2 <= radio**2 * (1 + 1e-9)
    elif forma == "cilindro":
        dentro = np.ones_like(r, dtype=bool)
    else:
        raise ValueError(f"forma de ROI desconocida: {forma!r}")
    pesos = np.where(r > 0, r, (rho[1] - rho[0]) / 8.0)
    return r[dentro], zz[dentro] + z_centro, pesos[dentro]


def homogeneidad(bz, pesos):
    """Campo medio, desviación RMS relativa y pico a pico relativo."""
    media = np.sum(pesos * bz) / np.sum(pesos)
    rms = np.sqrt(np.sum(pesos * (bz / media - 1.0) ** 2) / np.sum(pesos))
    return media, rms, (bz.max() - bz.min()) / media


def bz_lamina(bobinado, rho, z, n_anillos=600):
    z_anillos, pesos = bobinado.anillos(n_anillos)
    _, bz = campo_espira(bobinado.radio, z_anillos[:, None], rho[None, :], z[None, :])
    return pesos @ bz


def optimizar(radio, largo, z_centro, largo_cable, roi, orden, paso_minimo,
              monotona=True, objetivo="rms", densidad_minima=0.05, n_anillos=600):
    """Coeficientes de la distribución que minimizan la inhomogeneidad de B_z en la ROI.

    objetivo: "rms" (desviación cuadrática media) o "pico" (máximo menos mínimo).

    roi: tupla (rho, z, pesos) de `puntos_roi`.
    Restricciones (penalizadas): n(u) >= densidad_minima, paso local >= paso_minimo y,
    si `monotona`, densidad creciente desde el centro hacia cada extremo.
    """
    rho, z, pesos = roi
    bordes = np.linspace(-1.0, 1.0, n_anillos + 1)
    z_anillos = z_centro + 0.5 * (bordes[1:] + bordes[:-1]) * largo / 2.0
    _, matriz = campo_espira(radio, z_anillos[:, None], rho[None, :], z[None, :])
    vueltas_aprox = largo_cable / (2.0 * np.pi * radio)
    u_control = np.linspace(-1.0, 1.0, 401)

    def costo(c):
        d = Distribucion(c)
        penal = max(0.0, densidad_minima - d.densidad(u_control).min())
        paso_min = (largo / 2.0) / (vueltas_aprox * d.densidad_normalizada(u_control).max())
        penal += max(0.0, (paso_minimo - paso_min) / paso_minimo)
        if monotona:
            penal += max(0.0, -d.derivada(u_control[u_control > 0]).min())
        if penal > 0:
            return 1.0 + penal
        bz = np.diff(d.acumulada(bordes)) @ matriz
        return homogeneidad(bz, pesos)[2 if objetivo == "pico" else 1]

    x0 = np.zeros(orden)
    mejor = None
    for _ in range(3):  # reinicios de Nelder-Mead desde el mejor punto
        simplex = np.vstack([x0, x0 + 2.0 * np.eye(orden)])
        res = minimize(costo, x0, method="Nelder-Mead",
                       options={"initial_simplex": simplex, "xatol": 1e-6, "fatol": 1e-12,
                                "maxiter": 20000, "maxfev": 20000})
        if mejor is None or res.fun < mejor.fun:
            mejor = res
        x0 = mejor.x
    return mejor.x, mejor.fun


def optimizar_largo(radio, largo_max, z_centro, largo_cable, roi, orden, paso_minimo,
                    inhomogeneidad_max, monotona=True, objetivo="rms", tolerancia=0.5e-3):
    """Largo del bobinado libre: el bobinado más corto que cumple `inhomogeneidad_max`.

    Con el largo de cable fijo, acortar el bobinado junta las vueltas y sube el campo por
    amperio, pero empeora la homogeneidad. La inhomogeneidad mínima alcanzable crece al
    acortar, así que se busca por bisección el largo más corto que todavía la cumple.

    Devuelve (largo, coeficientes, inhomogeneidad, cumple). Si ni con `largo_max` se cumple
    el objetivo, devuelve el largo máximo con la máxima homogeneidad y cumple = False.
    """
    def resolver(largo):
        return optimizar(radio, largo, z_centro, largo_cable, roi, orden, paso_minimo,
                         monotona=monotona, objetivo=objetivo)

    c_max, f_max = resolver(largo_max)
    if f_max > inhomogeneidad_max:
        return largo_max, c_max, f_max, False
    vueltas = largo_cable / (2.0 * np.pi * radio)
    corto = max(1.02 * vueltas * paso_minimo, 0.2 * largo_max)
    if corto >= largo_max:
        return largo_max, c_max, f_max, True
    c, f = resolver(corto)
    if f <= inhomogeneidad_max:
        return corto, c, f, True
    largo, mejor_c, mejor_f = largo_max, c_max, f_max
    while largo - corto > tolerancia:
        medio = 0.5 * (largo + corto)
        c, f = resolver(medio)
        if f <= inhomogeneidad_max:
            largo, mejor_c, mejor_f = medio, c, f
        else:
            corto = medio
    return largo, mejor_c, mejor_f, True
