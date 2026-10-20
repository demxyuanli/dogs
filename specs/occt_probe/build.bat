@echo off
rem Build the OCCT ground-truth probe against the installed OCCT 8.0.0 (vc14 x64).
rem The probe is one executable built from several translation units:
rem   occt_probe.cpp      entry point and the `run*` dispatch loop
rem   modes_*.cpp         the `argv[2]` command groups
rem   *_probe.cpp         structural dumps / per-stage replays they call
rem   occt_probe_common.h shared includes, `ProbeCtx` and cross-file declarations
call "C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvars64.bat" >nul
if errorlevel 1 (echo vcvars failed & exit /b 1)
cd /d "%~dp0"
set "SRCS=occt_probe.cpp modes_mesh.cpp modes_pcurve.cpp modes_uv.cpp modes_faceids.cpp modes_edges.cpp modes_fix.cpp modes_report.cpp modes_boolean.cpp modes_classify.cpp modes_default.cpp delaun_probe.cpp boundary_probe.cpp structdump_probe.cpp srfparams_probe.cpp meshdbg_probe.cpp wiredbg_probe.cpp iges_probe.cpp extrema_probe.cpp"
cl /nologo /std:c++17 /EHsc /MD /O2 /MP /I "D:\source\occt-8.0.0\inc" %SRCS% _dbg\ZZ_TorusRangeSplitter.cxx ^
  /link /LIBPATH:"D:\source\occt-8.0.0\win64\vc14\lib" ^
  TKernel.lib TKMath.lib TKG2d.lib TKG3d.lib TKGeomBase.lib TKBRep.lib ^
  TKBO.lib TKBool.lib TKPrim.lib TKTopAlgo.lib TKGeomAlgo.lib TKShHealing.lib TKXSBase.lib TKDESTEP.lib TKDEIGES.lib TKMesh.lib TKOffset.lib ^
  /OUT:occt_probe.exe
if errorlevel 1 (echo link failed & exit /b 1)
echo BUILD OK
