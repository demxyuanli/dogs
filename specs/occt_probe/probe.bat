@echo off
rem Run the OCCT ground-truth probe. Requires the OCCT 8.0.0 distribution and its
rem third-party DLLs (tbb12.dll / jemalloc.dll come from THIRDPARTY_DIR).
rem THIRDPARTY_DIR must be ABSOLUTE: env.bat's default is relative to its cwd.
setlocal
set "OCCT=D:\source\occt-8.0.0"
set "THIRDPARTY_DIR=%OCCT%\3rdparty-vc14-64"
cd /d "%OCCT%"
call env.bat vc14 64 >nul
"%~dp0occt_probe.exe" %*