@echo off
rem Run the BRep modeling ground-truth generator. The exe is /MD, so it is
rem copied next to the OCCT DLLs (win64\vc14\bin) before it starts.
rem THIRDPARTY_DIR must be ABSOLUTE (env.bat's default is relative).
rem Usage: run_probe.bat <output-dir>   (default: data\modeling)
setlocal
set "OCCT=D:\source\occt-8.0.0"
set "THIRDPARTY_DIR=%OCCT%\3rdparty-vc14-64"
copy /y "%~dp0modeling_probe.exe" "%OCCT%\win64\vc14\bin\zz_modeling_probe.exe" >nul
if errorlevel 1 (echo copy failed & exit /b 1)
cd /d "%OCCT%"
call env.bat vc14 64 >nul
if "%~1"=="" (set "OUT=D:/source/repos/dogs/data/modeling") else (set "OUT=%~1")
"%OCCT%\win64\vc14\bin\zz_modeling_probe.exe" "%OUT%"
