@echo off
rem Run the OCCT ground-truth probe while KEEPING the caller's working directory.
rem `probe.bat` does `cd /d %OCCT%`, which forces absolute input paths; this one
rem pushes into the OCCT dir only long enough to source env.bat (PATH for the
rem OCCT + third-party DLLs), then returns and runs the probe in the caller's cwd.
setlocal
set "OCCT=D:\source\occt-8.0.0"
set "THIRDPARTY_DIR=%OCCT%\3rdparty-vc14-64"
pushd "%OCCT%"
call env.bat vc14 64 >nul
popd
"%~dp0occt_probe.exe" %*
exit /b %ERRORLEVEL%
