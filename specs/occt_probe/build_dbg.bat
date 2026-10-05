@echo off
call "C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvars64.bat" >nul
if errorlevel 1 (echo vcvars failed & exit /b 1)
cd /d "%~dp0"
cl /nologo /std:c++17 /EHsc /MD /O2 /I "D:\source\occt-8.0.0\inc" wires_probe.cpp _dbg\ZZ_ComposeShell.cxx _dbg\ZZ_ShapeFix_Face.cxx _dbg\ZZ_ActorRead.cxx _dbg\ZZ_Classifier.cxx ^
  /link /LIBPATH:"D:\source\occt-8.0.0\win64\vc14\lib" ^
  TKernel.lib TKMath.lib TKG2d.lib TKG3d.lib TKGeomBase.lib TKBRep.lib ^
  TKBO.lib TKBool.lib TKPrim.lib TKTopAlgo.lib TKGeomAlgo.lib TKShHealing.lib TKXSBase.lib TKDESTEP.lib TKMesh.lib ^
  /MAP:zz_wires_probe.map /OUT:zz_wires_probe.exe
if errorlevel 1 (echo link failed & exit /b 1)
echo BUILD OK
