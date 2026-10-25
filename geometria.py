"""Malla del soporte: cilindro hueco con surco helicoidal en U y salidas de terminales.

Todas las dimensiones en milímetros (unidad habitual de los STL).
"""
import numpy as np
import manifold3d as mf

EXCESO = 0.6  # cuánto sobresale el cortador del surco por fuera del cilindro [mm]
# El canal de salida queda apenas más hondo que el surco: si sus fondos coinciden, la
# operación booleana deja astillas que se aplastan al exportar el STL en float32 [mm]
HONDURA_EXTRA_CANAL = 0.03


def _solido(vertices, triangulos):
    """Crea un Manifold desde arreglos numpy y verifica que sea un sólido válido."""
    malla = mf.Mesh64(vert_properties=np.ascontiguousarray(vertices, dtype=np.float64),
                      tri_verts=np.ascontiguousarray(triangulos, dtype=np.uint32))
    solido = mf.Manifold(malla)
    if solido.status() != mf.Error.NoError:
        raise RuntimeError(f"malla inválida: {solido.status()}")
    return solido


def perfil_surco(r_centro, diametro, r_exterior, n_arco):
    """Perfil en U en el plano meridiano (rho, zeta), recorrido antihorario.

    Fondo semicircular de radio diametro/2 centrado en r_centro y paredes rectas que
    salen por fuera de la superficie exterior.
    """
    a = diametro / 2.0
    phi = np.linspace(np.pi / 2.0, 1.5 * np.pi, n_arco + 1)
    rho = np.append(r_centro + a * np.cos(phi), [r_exterior + EXCESO] * 2)
    zeta = np.append(a * np.sin(phi), [-a, a])
    return rho, zeta


def cortador_helicoidal(theta, z, r_centro, diametro, r_exterior, n_arco):
    """Sólido que se resta al cilindro para tallar el surco a lo largo de (theta, z).

    El perfil se coloca en cada plano meridiano; su altura se escala por 1/cos(alfa)
    (alfa = ángulo de la hélice) para que el ancho perpendicular al cable sea `diametro`.
    """
    rho, zeta = perfil_surco(r_centro, diametro, r_exterior, n_arco)
    p = len(rho)
    m = len(theta)
    dz = np.gradient(z, theta)
    escala = np.hypot(r_centro, dz) / r_centro
    x = rho[None, :] * np.cos(theta)[:, None]
    y = rho[None, :] * np.sin(theta)[:, None]
    zz = z[:, None] + zeta[None, :] * escala[:, None]
    vertices = np.column_stack([x.ravel(), y.ravel(), zz.ravel()])

    j = np.arange(m - 1)[:, None]
    k = np.arange(p)[None, :]
    k1 = (k + 1) % p
    a, b = j * p + k, j * p + k1
    c, d = (j + 1) * p + k, (j + 1) * p + k1
    laterales = np.concatenate([np.stack([a, d, b], -1).reshape(-1, 3),
                                np.stack([a, c, d], -1).reshape(-1, 3)])
    abanico = np.arange(1, p - 1)
    inicio = np.column_stack([np.zeros_like(abanico), abanico, abanico + 1])
    fin = (m - 1) * p + np.column_stack([np.zeros_like(abanico), abanico + 1, abanico])
    return _solido(vertices, np.concatenate([laterales, inicio, fin]))


def canal_axial(theta, z_desde, z_hasta, r_centro, diametro, r_exterior, n_arco):
    """Canal recto paralelo al eje (mismo perfil en U) en el ángulo theta [rad]."""
    a = diametro / 2.0
    phi = np.linspace(np.pi, 2.0 * np.pi, n_arco + 1)
    t = np.append(a * np.cos(phi), [a, -a])
    r = np.append(r_centro - HONDURA_EXTRA_CANAL + a * np.sin(phi), [r_exterior + EXCESO] * 2)
    seccion = mf.CrossSection([np.column_stack([t, r])])
    canal = mf.Manifold.extrude(seccion, z_hasta - z_desde).translate([0.0, 0.0, z_desde])
    return canal.rotate([0.0, 0.0, np.degrees(theta) - 90.0])


def agujero_radial(theta, z, r_desde, r_hasta, diametro, segmentos=32):
    """Agujero cilíndrico radial entre r_desde y r_hasta en (theta [rad], z)."""
    agujero = mf.Manifold.cylinder(r_hasta - r_desde, diametro / 2.0, diametro / 2.0, segmentos)
    agujero = agujero.rotate([0.0, 90.0, 0.0]).translate([r_desde, 0.0, 0.0])
    return agujero.rotate([0.0, 0.0, np.degrees(theta)]).translate([0.0, 0.0, z])


def astillas(solido, umbral=1e-4):
    """Triángulos con alguna arista menor que `umbral` [mm].

    Las operaciones booleanas pueden dejar aristas de ~1e-6 mm cuando una intersección pasa
    casi por un vértice; al exportar el STL en float32 esos triángulos se aplastan.
    """
    malla = solido.to_mesh64()
    p = np.asarray(malla.vert_properties)[:, :3][np.asarray(malla.tri_verts)]
    lados = np.linalg.norm(p - np.roll(p, 1, axis=1), axis=2)
    return int(np.sum(lados.min(axis=1) < umbral))


def tubo(r_interior, r_exterior, largo, segmentos, giro=0.0):
    """Cilindro hueco facetado; `giro` [rad] rota los vértices del facetado."""
    exterior = mf.Manifold.cylinder(largo, r_exterior, r_exterior, segmentos)
    interior = mf.Manifold.cylinder(largo + 2.0, r_interior, r_interior, segmentos)
    return (exterior - interior.translate([0.0, 0.0, -1.0])).rotate([0.0, 0.0, np.degrees(giro)])


def giro_facetado(angulos, segmentos):
    """Giro del facetado del tubo que lo aleja lo más posible de `angulos` [rad].

    Los vértices del tubo que caen casi sobre aristas del cortador generan astillas de
    ~1e-6 mm en la operación booleana, que se aplastan al exportar el STL en float32.
    """
    paso = 2.0 * np.pi / segmentos
    candidatos = np.linspace(0.0, paso, 512, endpoint=False)
    resto = np.mod(np.asarray(angulos)[None, :] - candidatos[:, None], paso)
    return candidatos[np.argmax(np.minimum(resto, paso - resto).min(axis=1))]


def guardar_stl(solido, ruta, nombre="bobina"):
    """Escribe un STL binario."""
    malla = solido.to_mesh64()
    v = np.asarray(malla.vert_properties)[:, :3]
    t = np.asarray(malla.tri_verts)
    p = v[t]
    normales = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    normales /= np.maximum(np.linalg.norm(normales, axis=1, keepdims=True), 1e-30)
    registro = np.dtype([("normal", "<f4", 3), ("vertices", "<f4", (3, 3)), ("atributo", "<u2")])
    datos = np.zeros(len(t), dtype=registro)
    datos["normal"] = normales
    datos["vertices"] = p
    with open(ruta, "wb") as f:
        f.write(nombre.encode("ascii", "replace")[:80].ljust(80, b" "))
        f.write(np.uint32(len(t)).tobytes())
        f.write(datos.tobytes())
    return len(t)


def seccion_meridiana(solido):
    """Polígonos (x, z) del corte del sólido por el plano y = 0."""
    corte = solido.rotate([90.0, 0.0, 0.0]).slice(0.0)
    return [np.column_stack([pol[:, 0], -pol[:, 1]]) for pol in corte.to_polygons()]
