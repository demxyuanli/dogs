# Regenerates an OCCT ground-truth OBJ reference for a data/*.step model,
# matching the header style of the existing data/occ-*.obj files
# ("# Exported by Open CASCADE Technology [dev.opencascade.org]").
#
# ---------------------------------------------------------------------------
# Provenance (established 2026-09-14)
# ---------------------------------------------------------------------------
# The legacy data/occ-*.obj references were produced by DRAWEXE (OCCT 8.0.0)
# with this exact pipeline:
#
#   ReadStep D <in.step>            ;# XDE/XCAF read: keeps product name + colors
#   XGetOneShape s D                ;# single compound of the document free shapes
#   bounding s -nodraw -print -noTriangulation -finitePart
#                                   ;# AABB, same algorithm as Prs3d::GetDeflection
#   incmesh s <lin> -angular 20     ;# BRepMesh_IncrementalMesh, absolute deflection
#   WriteObj D <out.obj>            ;# RWObj_CafWriter: group name = product name
#
# <lin> is the Prs3d relative deflection used by the visualization pipeline
# (StdPrs_ToolTriangulatedShape::Tessellate):
#   lin = max(maxComp(bbox) * DeviationCoefficient * 4, Precision::Confusion())
# with Prs3d_Drawer defaults DeviationCoefficient = 0.001 and DeviationAngle
# = 20 deg.  maxComp is the largest of the AABB side lengths (dx, dy, dz).
# References (OCCT source in this checkout):
#   Prs3d.hxx:70-72                     GetDeflection = maxComp*coeff*4
#   StdPrs_ToolTriangulatedShape.cxx:118 GetDeflection (BRepBndLib::Add(...,false))
#   StdPrs_ToolTriangulatedShape.cxx:160 Tessellate (relative deflection + angle)
#   Prs3d_Drawer.cxx                     DeviationCoefficient/DeviationAngle defaults
#   MeshTest.cxx:1717                    DRAWEXE "incmesh" command
#   XSDRAWOBJ.cxx:233                    DRAWEXE "WriteObj" command
#   BRepTest_BasicCommands.cxx:1830      DRAWEXE "bounding" command
#
# NOTE: do not use "incmesh s 0.001 -prs": that DRAW option also forces
# IMeshTools_Parameters::DeflectionInterior to 0.001, which over-refines the
# interior and does not match the references (e.g. Shape 20815/40700 instead of
# 6150/11372).  Pass the absolute linear deflection explicitly, as above.
#
# ---------------------------------------------------------------------------
# Invocation (Windows / PowerShell)
# ---------------------------------------------------------------------------
#   $root = "D:\source\occt-8.0.0"
#   $tp   = "$root\3rdparty-vc14-64"
#   $env:PATH = "$root\win64\vc14\bin;$tp\msvc-vc14-64;$tp\tcltk-8.6.15-x64\bin;" +
#               "$tp\freetype-2.13.3-x64\bin;$tp\freeimage-3.18.0-x64\bin;" +
#               "$tp\angle-gles2-2.1.0-vc14-64\bin;$tp\tbb-2021.13.0-x64\bin;" +
#               "$tp\vtk-9.4.1-x64\bin;$tp\ffmpeg-3.3.4-64\bin;" +
#               "$tp\jemalloc-vc14-64\bin;$tp\openvr-1.14.15-64\bin\win64;$env:PATH"
#   $env:CASROOT = $root
#   $env:CSF_OCCTResourcePath = "$root\src"
#   $env:DRAWHOME = "$root\src\DrawResources"
#   $env:CSF_DrawPluginDefaults = "$root\src\DrawResources"
#   $env:DRAWDEFAULT = "$root\src\DrawResources\DrawDefault"
#   $env:CSF_LANGUAGE = "us"; $env:MMGT_CLEAR = "1"
#   $env:OCC_IN  = "D:/source/repos/dogs/data/rev.step"
#   $env:OCC_OUT = "D:/source/repos/dogs/data/occ-rev.obj"
#   & "$root\win64\vc14\bin\DRAWEXE.exe" -b -f "D:/source/repos/dogs/data/_occ_ref_export.tcl"
#
# Optional overrides: OCC_COEFF (default 0.001), OCC_ANG (default 20).
#
# Reproduction check (2026-09-14, all exact against the stored data/occ-*.obj):
#   Cube 24/12, Cone 195/301, Cylinder 146/140, Sphere 642/1244, Torus 1369/2592,
#   rev 104/92, screw 600/790, linkrods 3494/5078, HoledPlate 180/128,
#   OffsetPlaneHoleEdge 48/32, Shape 6150/11372, Shape-1 3343/4336, Shape-2 3098/4778.
#   ATU01038 17720/22135 (regenerated with this recipe; the earlier file was
#   17726/22147 from a lost pipeline and has been overwritten).
#
# The group name and materials in the output come from the STEP product itself
# (XDE reader): files without a product name yield the default group
# "Open CASCADE STEP translator 7.6" and no mtllib, while named/colored
# products (e.g. ATU01038) yield the product name group plus an .mtl sidecar.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Provenance addendum (2026-09-24)
# ---------------------------------------------------------------------------
# The references for the four remaining `data/occ/*.step` models were added on
# 2026-09-24 20:59-21:00 (occ-bottom / occ-T0M / occ-TDB / occ-top, plus the
# IGES pair data/iges/occ-bearing.obj and data/iges/occ-hammer.obj at 20:58).
# They are real RWObj output -- header, product-name group, per-vertex `vn`,
# and a colour sidecar where the STEP carries colours -- in the same vertex
# format and the same frame as the references that were already stored
# (`data/occ/occ-*.obj`).  The batch that wrote them left no trace in this
# repository, so their exact producer build is unrecorded.  Header counts as
# committed: top 16905/22800, bottom 17587/23557, T0M 60548/67366,
# TDB 73544/79977.
#
# Frame caveat -- do NOT overwrite a stored reference with a blind re-run of
# this script: run with the OCCT 8.0.0 install above it reproduces the *same*
# mesh (top: 16905/22800, every stored vertex equal to the run's vertex rotated
# by R_x(90); bottom: identical counts; T0M/TDB: +-1 vertex; the .mtl sidecars
# are byte-identical) but in a frame rotated 90 degrees about X relative to
# every stored reference.  That is the same signature as board card T-10
# (`data/occ-ATU01038.obj` was once overwritten in place by such a run).  All
# stored references, and the Rust port's own exports, share the un-rotated
# frame, so a fresh run has to be validated against a stored reference before
# it is committed.  The "Reproduction check (2026-09-14)" list above only covers
# the non-assembly `data/*.step` models, which carry no root placement, so it
# did not exercise this difference.
# ---------------------------------------------------------------------------

set in    $::env(OCC_IN)
set out   $::env(OCC_OUT)
set coeff 0.001
set ang   20.0
if {[info exists ::env(OCC_COEFF)]} { set coeff $::env(OCC_COEFF) }
if {[info exists ::env(OCC_ANG)]}   { set ang   $::env(OCC_ANG) }

pload ALL
ReadStep D $in
XGetOneShape s D

set bbox [bounding s -nodraw -print -noTriangulation -finitePart]
set dx 0.0
set dy 0.0
set dz 0.0
foreach l [split $bbox "\n"] {
  if {[regexp {X-range:\s+(\S+)\s+(\S+)} $l m a b]} { set dx [expr {$b - $a}] }
  if {[regexp {Y-range:\s+(\S+)\s+(\S+)} $l m a b]} { set dy [expr {$b - $a}] }
  if {[regexp {Z-range:\s+(\S+)\s+(\S+)} $l m a b]} { set dz [expr {$b - $a}] }
}
set mc [expr {$dx > $dy ? ($dx > $dz ? $dx : $dz) : ($dy > $dz ? $dy : $dz)}]
set lin [expr {$mc * $coeff * 4.0}]
if {$lin < 1.0e-7} { set lin 1.0e-7 }
puts "OCC_LIN=$lin maxComp=$mc dx=$dx dy=$dy dz=$dz"

incmesh s $lin -angular $ang
WriteObj D $out
puts "OCC_DONE=$out"
exit
