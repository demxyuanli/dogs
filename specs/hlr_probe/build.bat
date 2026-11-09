@echo off
rem Build the polygonal HLR ground-truth generator against the installed
rem OCCT 8.0.0 (vc14 x64). Produces hlr_probe.exe next to this script.
call "C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvars64.bat" >nul
if errorlevel 1 (echo vcvars failed & exit /b 1)
cd /d "%~dp0"
cl /nologo /std:c++17 /EHsc /MD /O2 /I "D:\source\occt-8.0.0\inc" hlr_probe.cpp ^
  /link /LIBPATH:"D:\source\occt-8.0.0\win64\vc14\lib" ^
  TKernel.lib TKMath.lib TKG2d.lib TKG3d.lib TKGeomBase.lib TKBRep.lib ^
  TKTopAlgo.lib TKMesh.lib TKGeomAlgo.lib TKHLR.lib TKXSBase.lib ^
  TKDESTEP.lib TKDE.lib ^
  /OUT:hlr_probe.exe
if errorlevel 1 (echo link failed & exit /b 1)
echo BUILD OK
