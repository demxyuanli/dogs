@echo off
rem Run bop_split_probe.exe: copy it next to the OCCT DLLs and set the env.
set "OCCT=D:\source\occt-8.0.0"
copy /y "%~dp0bop_split_probe.exe" "%OCCT%\win64\vc14\bin\zz_bop_split_probe.exe" >nul
pushd "%OCCT%"
call env.bat vc14 64 >nul
popd
"%OCCT%\win64\vc14\bin\zz_bop_split_probe.exe" %*
