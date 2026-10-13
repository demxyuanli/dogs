@echo off
rem Run DRAWEXE (the OCCT DRAW interpreter) with the OCCT DLLs on PATH.
setlocal
set "OCCT=D:\source\occt-8.0.0"
set "THIRDPARTY_DIR=%OCCT%\3rdparty-vc14-64"
pushd "%OCCT%"
call env.bat vc14 64 >nul
popd
"%OCCT%\win64\vc14\bin\DRAWEXE.exe" %*
exit /b %ERRORLEVEL%
