@echo off
rem Run the polygonal HLR ground-truth generator. The exe is /MD, so it is
rem copied next to the OCCT DLLs (win64\vc14\bin) before it starts.
rem THIRDPARTY_DIR must be ABSOLUTE (env.bat's default is relative).
rem Usage: run_probe.bat <step-file> [deflection] [angle-rad] [xoy|xoz|yoz]
setlocal
set "OCCT=D:\source\occt-8.0.0"
set "THIRDPARTY_DIR=%OCCT%\3rdparty-vc14-64"
copy /y "%~dp0hlr_probe.exe" "%OCCT%\win64\vc14\bin\zz_hlr_probe.exe" >nul
if errorlevel 1 (echo copy failed & exit /b 1)
cd /d "%OCCT%"
call env.bat vc14 64 >nul
"%OCCT%\win64\vc14\bin\zz_hlr_probe.exe" %*
