"""Campo magnético e inductancia de la bobina (unidades SI: metros, amperios, teslas).

- campo_espira: campo exacto de una espira circular (integrales elípticas).
- campo_polilinea: Biot-Savart de un conductor filamentario descrito como polilínea.
- inductancia_espiras: autoinductancia aproximando cada vuelta por una espira coaxial.
"""
import numpy as np
from scipy.special import ellipe, ellipk

MU0 = 4e-7 * np.pi


def campo_espira(a, z0, rho, z, corriente=1.0):
    """Campo (B_rho, B_z) de una espira de radio `a` centrada en el eje a la altura `z0`.

    `rho` y `z` son arreglos con las coordenadas cilíndricas de los puntos de observación.
    """
    rho = np.asarray(rho, dtype=float)
    dz = np.asarray(z, dtype=float) - z0
    alfa2 = (a - rho) ** 2 + dz**2
    beta2 = (a + rho) ** 2 + dz**2
    beta = np.sqrt(beta2)
    m = 4.0 * a * rho / beta2
    K = ellipk(m)
    E = ellipe(m)
    c = MU0 * corriente / (2.0 * np.pi * beta)
    bz = c * (K + (a**2 - rho**2 - dz**2) / alfa2 * E)
    rho_seguro = np.where(rho > 0, rho, 1.0)
    br = np.where(rho > 0, c * dz / rho_seguro * (-K + (a**2 + rho**2 + dz**2) / alfa2 * E), 0.0)
    return br, bz


def campo_polilinea(cable, puntos, corriente=1.0, bloque=128):
    """Campo B (N x 3) en `puntos` (N x 3) producido por la polilínea `cable` (M x 3).

    Usa la expresión cerrada de Biot-Savart para cada segmento recto.
    """
    a = cable[:-1]
    b = cable[1:]
    salida = np.zeros((len(puntos), 3))
    for i in range(0, len(puntos), bloque):
        p = puntos[i : i + bloque, None, :]
        r1 = a[None] - p
        r2 = b[None] - p
        n1 = np.linalg.norm(r1, axis=-1)
        n2 = np.linalg.norm(r2, axis=-1)
        producto = np.cross(r1, r2)
        denominador = n1 * n2 * (n1 * n2 + np.sum(r1 * r2, axis=-1))
        factor = (n1 + n2) / denominador
        salida[i : i + bloque] = np.sum(producto * factor[..., None], axis=1)
    return salida * MU0 * corriente / (4.0 * np.pi)


def inductancia_mutua_coaxial(a, b, d):
    """Inductancia mutua entre dos espiras coaxiales de radios a, b separadas d."""
    k2 = 4.0 * a * b / ((a + b) ** 2 + d**2)
    k = np.sqrt(k2)
    return MU0 * np.sqrt(a * b) * ((2.0 / k - k) * ellipk(k2) - 2.0 / k * ellipe(k2))


def inductancia_espiras(radio, z_vueltas, pesos, radio_hilo):
    """Autoinductancia [H] de un conjunto de vueltas coaxiales de igual radio.

    `pesos` permite contar una vuelta incompleta (fracción de vuelta). La autoinductancia
    de cada espira usa la aproximación de alta frecuencia (corriente superficial):
    L = mu0 R (ln(8R/a) - 2).
    """
    z = np.asarray(z_vueltas, dtype=float)
    w = np.asarray(pesos, dtype=float)
    propia = MU0 * radio * (np.log(8.0 * radio / radio_hilo) - 2.0)
    d = np.abs(z[:, None] - z[None, :])
    np.fill_diagonal(d, 1.0)
    m = inductancia_mutua_coaxial(radio, radio, d)
    np.fill_diagonal(m, propia)
    return float(w @ m @ w)
