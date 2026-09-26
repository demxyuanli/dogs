@echo off
rem Build (cl.exe + OCCT libs) then run the probe:  build_run.bat <step-file> [mode]
call "%~dp0build.bat"
if errorlevel 1 exit /b 1
call "%~dp0run_probe.bat" %*
if errorlevel 1 exit /b 1
call "%~dp0probe.bat" %*