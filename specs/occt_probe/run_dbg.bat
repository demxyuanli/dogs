@echo off
set "OCCT=D:\source\occt-8.0.0"
set "THIRDPARTY_DIR=%OCCT%\3rdparty-vc14-64"
pushd "%OCCT%"
call env.bat vc14 64 >nul
popd
"%~dp0zz_wires_probe.exe" %*
