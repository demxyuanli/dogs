# Regenerates an OCCT ground-truth glTF reference for a data/*.step model.
#
# Same pipeline as data/_occ_ref_export.tcl up to the mesh, then
# RWGltf_CafWriter instead of RWObj_CafWriter:
#
#   ReadStep D <in.step>            ;# XDE/XCAF read
#   XGetOneShape s D                ;# document free shapes
#   bounding s -nodraw -print -noTriangulation -finitePart
#   incmesh s <lin> -angular 20     ;# BRepMesh_IncrementalMesh
#   WriteGltf D <out.gltf>          ;# RWGltf_CafWriter, defaults
#
# WriteGltf defaults (XSDRAWGLTF.cxx:300-503): non-binary .gltf (a .glb
# extension switches to binary), trsfFormat compact, systemCoordSys Zup (so no
# coordinate conversion — same frame as the OBJ references), nodeNameFormat
# instOrProd, meshNameFormat product, no forced UV export, no face merging.
#
# Invocation (Windows / PowerShell): same environment as _occ_ref_export.tcl,
# with OCC_IN / OCC_OUT for the .gltf path. The .bin buffer is written next to
# it as <stem>.bin.

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
if {[info exists ::env(OCC_FORCE_UV)] && $::env(OCC_FORCE_UV) ne "" && $::env(OCC_FORCE_UV) ne "0"} {
  WriteGltf D $out -forceUVExport
} else {
  WriteGltf D $out
}
puts "OCC_DONE=$out"
exit
