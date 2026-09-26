@echo off
rem Build the OCCT ground-truth probe against the installed OCCT 8.0.0 (vc14 x64).
call "C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvars64.bat" >nul
if errorlevel 1 (echo vcvars failed & exit /b 1)
cd /d "%~dp0"
cl /nologo /std:c++17 /EHsc /MD /O2 /I "D:\source\occt-8.0.0\inc" occt_probe.cpp ^
  /link /LIBPATH:"D:\source\occt-8.0.0\win64\vc14\lib" ^
  TKernel.lib TKMath.lib TKG2d.lib TKG3d.lib TKGeomBase.lib TKBRep.lib ^
  TKBO.lib TKBool.lib TKPrim.lib TKTopAlgo.lib TKGeomAlgo.lib TKShHealing.lib TKXSBase.lib TKDESTEP.lib ^
  /OUT:occt_probe.exe
if errorlevel 1 (echo link failed & exit /b 1)
echo BUILD OK
