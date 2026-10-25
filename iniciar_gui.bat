@echo off
rem Doble clic para abrir el generador de antenas.
rem Si faltan Python o las dependencias, pregunta en una ventana y las instala (ver iniciar.ps1).
cd /d "%~dp0"
title Generador de antena solenoide
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0iniciar.ps1" %*
if errorlevel 1 pause
