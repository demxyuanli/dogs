@echo off
rem Build the BRep modeling ground-truth generator against the installed
rem OCCT 8.0.0 (vc14 x64). Produces modeling_probe.exe next to this script.
call "C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvars64.bat" >nul
if errorlevel 1 (echo vcvars failed & exit /b 1)
cd /d "%~dp0"
cl /nologo /std:c++17 /EHsc /MD /O2 /I "D:\source\occt-8.0.0\inc" modeling_probe.cpp ^
  /link /LIBPATH:"D:\source\occt-8.0.0\win64\vc14\lib" ^
  TKernel.lib TKMath.lib TKG2d.lib TKG3d.lib TKGeomBase.lib TKBRep.lib ^
  TKTopAlgo.lib TKPrim.lib TKMesh.lib TKGeomAlgo.lib TKBO.lib TKBool.lib ^
  TKFillet.lib TKOffset.lib TKShHealing.lib ^
  /OUT:modeling_probe.exe
if errorlevel 1 (echo link failed & exit /b 1)
echo BUILD OK
