@echo off
rem The probe'\\''s CRT is /MD: copy it next to the OCCT DLLs so it can start.
set "OCCT=D:\source\occt-8.0.0"
copy /y "%~dp0occt_probe.exe" "%OCCT%\win64\vc14\bin\zz_occt_probe.exe" >nul