@echo off
rem Run zz_iges_obj.exe with the OCCT + third-party DLLs on PATH.
setlocal
set "OCCT=D:\source\occt-8.0.0"
set "THIRDPARTY_DIR=%OCCT%\3rdparty-vc14-64"
pushd "%OCCT%"
call env.bat vc14 64 >nul
popd
"%~dp0zz_iges_obj.exe" %*
exit /b %ERRORLEVEL%
