"""Interfaz web local para generar bobinas.

Uso:  python servidor_gui.py [--puerto 8765] [--no-abrir]
      (o iniciar_gui.bat). Abre http://127.0.0.1:<puerto> en el navegador.

Sólo escucha en 127.0.0.1: la interfaz no es accesible desde otros equipos.
"""
import argparse
import json
import mimetypes
import os
import sys
import threading
import time
import traceback
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, unquote, urlparse

import matplotlib

matplotlib.use("Agg")

import generar_bobina as gb  # noqa: E402

PAGINA = Path(__file__).resolve().parent / "gui" / "index.html"

# Parámetros que la interfaz puede modificar (el resto queda con su valor por defecto)
CAMPOS_GUI = {
    "diametro_interno_mm": float, "diametro_externo_mm": float, "largo_cilindro_mm": float,
    "margen_extremos_mm": float, "diametro_surco_mm": float, "diametro_cable_mm": float,
    "largo_cable_m": float, "paso_minimo_mm": float, "roi_diametro_mm": float,
    "roi_largo_mm": float, "objetivo": str, "largo_libre": bool, "inhomogeneidad_max_pct": float,
    "segmentos_por_vuelta": int, "segmentos_arco_surco": int,
}


class Trabajo:
    """Estado de la generación en curso (una a la vez)."""

    def __init__(self):
        self.lock = threading.Lock()
        self.estado = "inactivo"
        self.proyecto = None
        self.log = []
        self.error = None
        self.inicio = None
        self.fin = None

    def agregar(self, linea):
        with self.lock:
            self.log.append(f"[{time.time() - self.inicio:5.1f} s] {linea}")

    def foto(self):
        with self.lock:
            fin = self.fin or time.time()
            return {"estado": self.estado, "proyecto": self.proyecto, "log": list(self.log),
                    "error": self.error,
                    "segundos": round(fin - self.inicio, 1) if self.inicio else None}

    def lanzar(self, proyecto, carpeta, parametros):
        with self.lock:
            if self.estado == "trabajando":
                return False
            self.estado, self.proyecto, self.log, self.error = "trabajando", proyecto, [], None
            self.inicio, self.fin = time.time(), None
        threading.Thread(target=self._correr, args=(carpeta, parametros), daemon=True).start()
        return True

    def _correr(self, carpeta, parametros):
        try:
            gb.generar(parametros, carpeta, log=self.agregar)
            estado, error = "listo", None
        except Exception as exc:  # se informa en la interfaz
            estado, error = "error", str(exc)
            self.agregar("ERROR: " + str(exc))
            traceback.print_exc()
        with self.lock:
            self.estado, self.error, self.fin = estado, error, time.time()


TRABAJO = Trabajo()


def leer_parametros(datos):
    """Parametros a partir del JSON de la interfaz; lanza ValueError con mensajes legibles."""
    valores = {}
    for campo, tipo in CAMPOS_GUI.items():
        if campo not in datos:
            continue
        valor = datos[campo]
        if tipo is bool:
            valores[campo] = bool(valor)
        elif tipo is str:
            valores[campo] = str(valor)
        else:
            try:
                numero = float(valor)
            except (TypeError, ValueError):
                raise ValueError(f"{campo}: '{valor}' no es un número.") from None
            if tipo is int:
                if not numero.is_integer():
                    raise ValueError(f"{campo}: '{valor}' debe ser un número entero.")
                numero = int(numero)
            valores[campo] = numero
    p = gb.Parametros(**valores)
    errores = p.validar()
    if errores:
        raise ValueError("\n".join(errores))
    return p


def listar_proyectos():
    if not gb.CARPETA_PROYECTOS.exists():
        return []
    proyectos = []
    for carpeta in gb.CARPETA_PROYECTOS.iterdir():
        resumen = carpeta / "parametros.json"
        if carpeta.is_dir():
            # la fecha de la carpeta no cambia al sobrescribir archivos: se usa la del resultado
            marca = (resumen if resumen.exists() else carpeta).stat().st_mtime
            proyectos.append({"nombre": carpeta.name, "completo": resumen.exists(),
                              "fecha": time.strftime("%Y-%m-%d %H:%M", time.localtime(marca))})
    return sorted(proyectos, key=lambda x: x["fecha"], reverse=True)


def info_proyecto(nombre):
    carpeta = gb.carpeta_proyecto(nombre)
    datos_json = carpeta / "parametros.json"
    if not datos_json.exists():
        raise FileNotFoundError(f"El proyecto '{nombre}' no tiene resultados.")
    datos = json.loads(datos_json.read_text(encoding="utf-8"))
    resumen = carpeta / "Resumen.txt"
    datos["resumen"] = resumen.read_text(encoding="utf-8") if resumen.exists() else ""
    datos["archivos"] = [
        {"nombre": nombre_arch, "descripcion": descripcion,
         "bytes": (carpeta / nombre_arch).stat().st_size,
         "version": int((carpeta / nombre_arch).stat().st_mtime)}
        for nombre_arch, descripcion in gb.ARCHIVOS if (carpeta / nombre_arch).exists()
    ]
    datos["carpeta"] = str(carpeta)
    return datos


class Manejador(BaseHTTPRequestHandler):
    server_version = "BobinaGUI/1.0"

    def log_message(self, formato, *args):  # silencia el registro de cada petición
        pass

    def _responder(self, codigo, cuerpo, tipo="application/json; charset=utf-8", extra=None):
        if isinstance(cuerpo, (dict, list)):
            cuerpo = json.dumps(cuerpo, ensure_ascii=False).encode("utf-8")
        elif isinstance(cuerpo, str):
            cuerpo = cuerpo.encode("utf-8")
        self.send_response(codigo)
        self.send_header("Content-Type", tipo)
        self.send_header("Content-Length", str(len(cuerpo)))
        self.send_header("Cache-Control", "no-store")
        for clave, valor in (extra or {}).items():
            self.send_header(clave, valor)
        self.end_headers()
        self.wfile.write(cuerpo)

    def _error(self, codigo, mensaje):
        self._responder(codigo, {"error": mensaje})

    def do_GET(self):
        url = urlparse(self.path)
        ruta = unquote(url.path)
        try:
            if ruta in ("/", "/index.html"):
                self._responder(200, PAGINA.read_bytes(), "text/html; charset=utf-8")
            elif ruta == "/api/defaults":
                p = gb.Parametros()
                self._responder(200, {"valores": {c: getattr(p, c) for c in CAMPOS_GUI},
                                      "margen_extremos_mm": p.margen_extremos_mm})
            elif ruta == "/api/estado":
                self._responder(200, TRABAJO.foto())
            elif ruta == "/api/proyectos":
                self._responder(200, listar_proyectos())
            elif ruta.startswith("/api/proyecto/"):
                self._responder(200, info_proyecto(ruta[len("/api/proyecto/"):]))
            elif ruta.startswith("/proyectos/"):
                self._archivo(ruta, parse_qs(url.query))
            else:
                self._error(404, "No encontrado")
        except (ValueError, FileNotFoundError) as exc:
            self._error(404, str(exc))

    def _archivo(self, ruta, consulta):
        partes = ruta.split("/")
        if len(partes) != 4:
            return self._error(404, "No encontrado")
        carpeta = gb.carpeta_proyecto(partes[2])
        nombres = {n for n, _ in gb.ARCHIVOS}
        archivo = carpeta / partes[3]
        if partes[3] not in nombres or not archivo.is_file():
            return self._error(404, "No encontrado")
        tipo = {".stl": "model/stl", ".csv": "text/csv; charset=utf-8",
                ".txt": "text/plain; charset=utf-8", ".json": "application/json; charset=utf-8"}.get(
            archivo.suffix, mimetypes.guess_type(archivo.name)[0] or "application/octet-stream")
        extra = {}
        if "descargar" in consulta:
            extra["Content-Disposition"] = f'attachment; filename="{carpeta.name}_{archivo.name}"'
        self._responder(200, archivo.read_bytes(), tipo, extra)

    def do_POST(self):
        ruta = urlparse(self.path).path
        try:
            largo = int(self.headers.get("Content-Length", 0))
            datos = json.loads(self.rfile.read(largo) or b"{}")
        except (ValueError, json.JSONDecodeError):
            return self._error(400, "Petición no válida")
        try:
            if ruta == "/api/generar":
                carpeta = gb.carpeta_proyecto(datos.get("proyecto"))
                parametros = leer_parametros(datos)
                existia = carpeta.exists()
                if not TRABAJO.lanzar(carpeta.name, carpeta, parametros):
                    return self._error(409, "Ya hay una generación en curso; espere a que termine.")
                self._responder(202, {"proyecto": carpeta.name, "existia": existia})
            elif ruta == "/api/abrir":
                carpeta = gb.carpeta_proyecto(datos.get("proyecto"))
                if not carpeta.exists():
                    return self._error(404, "La carpeta no existe.")
                if hasattr(os, "startfile"):
                    os.startfile(carpeta)
                self._responder(200, {"carpeta": str(carpeta)})
            else:
                self._error(404, "No encontrado")
        except ValueError as exc:
            self._error(400, str(exc))


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    analizador = argparse.ArgumentParser(description=__doc__)
    analizador.add_argument("--puerto", type=int, default=8765)
    analizador.add_argument("--no-abrir", action="store_true", help="no abrir el navegador")
    args = analizador.parse_args()

    servidor = None
    for puerto in range(args.puerto, args.puerto + 20):
        try:
            servidor = ThreadingHTTPServer(("127.0.0.1", puerto), Manejador)
            break
        except OSError:
            continue
    if servidor is None:
        sys.exit(f"No hay puertos libres entre {args.puerto} y {args.puerto + 19}.")
    url = f"http://127.0.0.1:{servidor.server_address[1]}/"
    print(f"Interfaz de bobinas en {url}")
    print("Cierre esta ventana (o Ctrl+C) para detenerla.")
    if not args.no_abrir and not os.environ.get("BOBINA_NO_ABRIR"):
        threading.Timer(1.0, webbrowser.open, [url]).start()
    try:
        servidor.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        servidor.server_close()


if __name__ == "__main__":
    main()
