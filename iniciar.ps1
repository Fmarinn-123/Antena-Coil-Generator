# Inicia el generador de antenas: prepara Python y las dependencias (preguntando antes en
# una ventana) y abre la interfaz en el navegador. Se ejecuta desde iniciar_gui.bat.
#
# Variables de entorno opcionales:
#   BOBINA_AUTO_INSTALAR=1  responde "Sí" a las preguntas (instalación desatendida)
#   BOBINA_NO_ABRIR=1       no abre el navegador
param([Parameter(ValueFromRemainingArguments = $true)] $Resto)

$ErrorActionPreference = 'Continue'
Set-Location -LiteralPath $PSScriptRoot
$Titulo = 'Generador de antena solenoide'
$VenvPython = Join-Path $PSScriptRoot '.venv\Scripts\python.exe'
$Modulos = 'import numpy, scipy, matplotlib, manifold3d'

function Preguntar([string]$Texto) {
    if ($env:BOBINA_AUTO_INSTALAR) { return $true }
    Add-Type -AssemblyName System.Windows.Forms
    $respuesta = [System.Windows.Forms.MessageBox]::Show($Texto, $Titulo, 'YesNo', 'Question')
    return $respuesta -eq 'Yes'
}

function Avisar([string]$Texto) {
    Write-Host $Texto -ForegroundColor Yellow
    if ($env:BOBINA_AUTO_INSTALAR) { return }
    Add-Type -AssemblyName System.Windows.Forms
    [void][System.Windows.Forms.MessageBox]::Show($Texto, $Titulo, 'OK', 'Warning')
}

function Es-PythonValido([string]$Exe, [string[]]$Previos) {
    # Rechaza también el acceso directo de Microsoft Store, que no es un Python real
    & $Exe @Previos -c 'import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)' 2>$null | Out-Null
    return $LASTEXITCODE -eq 0
}

function Buscar-Python {
    $candidatos = @()
    if (Get-Command py -ErrorAction SilentlyContinue) { $candidatos += , @('py', @('-3')) }
    if (Get-Command python -ErrorAction SilentlyContinue) { $candidatos += , @('python', @()) }
    $carpetas = @("$env:LOCALAPPDATA\Programs\Python\Python3*\python.exe",
                  "$env:LOCALAPPDATA\Python\pythoncore-3*\python.exe",
                  "$env:ProgramFiles\Python3*\python.exe")
    foreach ($patron in $carpetas) {
        Get-ChildItem -Path $patron -ErrorAction SilentlyContinue | Sort-Object FullName -Descending |
            ForEach-Object { $candidatos += , @($_.FullName, @()) }
    }
    foreach ($c in $candidatos) {
        if (Es-PythonValido $c[0] $c[1]) { return $c }
    }
    return $null
}

function Dependencias-Listas {
    if (-not (Test-Path -LiteralPath $VenvPython)) { return $false }
    & $VenvPython -c $Modulos 2>$null | Out-Null
    return $LASTEXITCODE -eq 0
}

if (-not (Dependencias-Listas)) {
    $python = Buscar-Python
    if (-not $python) {
        $texto = "No se encontró Python 3.10 o superior en este equipo.`n`n" +
                 "¿Desea instalar Python 3.12 ahora? Se usa winget (el instalador de Windows) y tarda unos minutos."
        if (-not (Preguntar $texto)) {
            Avisar "Sin Python no se puede abrir el generador. Puede instalarlo desde https://www.python.org/downloads/ y volver a ejecutar iniciar_gui.bat."
            exit 1
        }
        if (-not (Get-Command winget -ErrorAction SilentlyContinue)) {
            Avisar "winget no está disponible en este equipo. Instale Python desde https://www.python.org/downloads/ (marque 'Add python.exe to PATH') y vuelva a ejecutar iniciar_gui.bat."
            Start-Process 'https://www.python.org/downloads/'
            exit 1
        }
        Write-Host 'Instalando Python 3.12 con winget...'
        winget install --exact --id Python.Python.3.12 --scope user --accept-package-agreements --accept-source-agreements
        $python = Buscar-Python
        if (-not $python) {
            Avisar "Python se instaló pero no se encuentra todavía. Cierre esta ventana y vuelva a ejecutar iniciar_gui.bat."
            exit 1
        }
    }

    $texto = "Falta instalar las dependencias del generador: numpy, scipy, matplotlib y manifold3d.`n`n" +
             "Se instalarán en la carpeta .venv de este proyecto, sin tocar el resto del sistema. " +
             "Se descargan unos 100 MB y tarda de 1 a 3 minutos.`n`n¿Instalarlas ahora?"
    if (-not (Preguntar $texto)) {
        Avisar 'Instalación cancelada: el generador necesita esas dependencias para funcionar.'
        exit 1
    }
    Write-Host 'Creando el entorno local .venv ...'
    & $python[0] @($python[1]) -m venv --clear .venv
    if (-not (Test-Path -LiteralPath $VenvPython)) {
        Avisar 'No se pudo crear el entorno .venv. Revise que Python esté bien instalado.'
        exit 1
    }
    Write-Host 'Instalando dependencias (puede tardar unos minutos)...'
    & $VenvPython -m pip install --disable-pip-version-check -r requirements.txt
    if (-not (Dependencias-Listas)) {
        Avisar 'No se pudieron instalar las dependencias. Revise la conexión a internet y vuelva a intentarlo.'
        exit 1
    }
    Write-Host 'Dependencias instaladas.' -ForegroundColor Green
}

Write-Host 'Iniciando la interfaz. Cierre esta ventana para detenerla.'
& $VenvPython servidor_gui.py @Resto
exit $LASTEXITCODE
