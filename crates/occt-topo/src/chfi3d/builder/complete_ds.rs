//! Port of the `CompleteDS` static of `ChFi3d_Builder.cxx:80-136`.
//!
//! It is not a member of `ChFi3d_Builder` in OCCT either: the file-local static
//! takes the DS and the input shape, adds the ancestor wire of every DS edge
//! and the ancestor shell of every DS face, then sets the range of each DS
//! curve from the curve interferences it holds.

use crate::chfi_ds::ChFiDSMap;
use crate::shape::TopoShape;
use crate::topo_tools::explorer_shapes;
use crate::abs::ShapeType;
use crate::topopebrep::ds::HDataStructure;

/// `CompleteDS(DStr, S)` (`ChFi3d_Builder.cxx:80-136`).
pub fn complete_ds(dstr: &mut HDataStructure, s: &TopoShape) {
    let mut map_ew = ChFiDSMap::new();
    map_ew.fill(s, ShapeType::Edge, ShapeType::Wire);
    let mut map_fs = ChFiDSMap::new();
    map_fs.fill(s, ShapeType::Face, ShapeType::Shell);

    // `ChFi3d_Builder.cxx:87-103`: the ancestor wires of a DS edge.
    for e in explorer_shapes(s, ShapeType::Edge) {
        if !dstr.has_geometry(&e) {
            continue;
        }
        for wire_anc in map_ew.get_from_key(&e) {
            dstr.change_ds().add_shape(wire_anc);
        }
    }

    // `ChFi3d_Builder.cxx:106-122`: the ancestor shells of a DS face.
    for f in explorer_shapes(s, ShapeType::Face) {
        if !dstr.has_geometry(&f) {
            continue;
        }
        for shell_anc in map_fs.get_from_key(&f) {
            dstr.change_ds().add_shape(shell_anc);
        }
    }

    // `ChFi3d_Builder.cxx:124-135`: the range of every DS curve. An empty
    // interference list keeps `RealLast()` / `RealFirst()` as the range,
    // exactly like the OCCT loop that never enters its body.
    for ic in 1..=dstr.nb_curves() {
        let interferences = dstr.ds().curve_interferences(ic).to_vec();
        let range = dstr.min_max_on_parameter(&interferences).ok().flatten();
        let (parmin, parmax) = range.unwrap_or((f64::MAX, -f64::MAX));
        dstr.change_curve(ic).set_range(parmin, parmax);
    }
}
