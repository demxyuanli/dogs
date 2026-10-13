puts "=== mesh-ish commands ==="
puts [lsort [info commands *mesh*]]
puts [lsort [info commands *tri*]]
puts "=== read ==="
set x [ReadIges _igesdump/oracle/Torus.iges]
puts "returned=<$x>"
nbshapes $x
puts "=== incmesh ==="
incmesh $x 0.096 0.3490658503988659
puts "=== bbox after incmesh ==="
bounding $x
puts "=== done ==="
