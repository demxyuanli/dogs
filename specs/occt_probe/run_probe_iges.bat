@echo off
rem Run the OCCT ground-truth probe with the OCCT + third-party DLLs on PATH
rem (same setup as run_here.bat), keeping the caller's working directory.
setlocal
set "OCCT=D:\source\occt-8.0.0"
set "THIRDPARTY_DIR=%OCCT%\3rdparty-vc14-64"
pushd "%OCCT%"
call env.bat vc14 64 >nul
popd
"%~dp0occt_probe.exe" %*
exit /b %ERRORLEVEL%
