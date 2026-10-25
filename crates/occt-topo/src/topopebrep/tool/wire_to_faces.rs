//! `TopOpeBRepTool_TOOL::WireToFace` — rebuild faces from a wire map.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepTool/TopOpeBRepTool_TOOL.cxx`
//! (OCCT 8.0.0, 1943-1974). The map binds a boundary wire to the wires that
//! describe the holes inside it (`TopOpeBRepTool_RegularizeW.cxx:855-877`).

use crate::abs::Orientation;
use crate::builder::TopoBuilder;
use crate::shape::TopoShape;
use crate::tgeometry::GeometryRegistry;

/// `TopOpeBRepTool_TOOL::WireToFace(Fref, mapWlow, lFs)` (`cxx:1943-1974`).
///
/// `map_wlow` binds each boundary wire to its holes; `l_fs` is cleared and
/// filled with one face per boundary. The faces are empty copies of `f_ref`
/// (so they carry its surface, tolerance and natural-restriction flag) with
/// the wires added, reversed when `f_ref` is reversed.
pub fn wire_to_face(
    f_ref: &TopoShape,
    map_wlow: &[(TopoShape, Vec<TopoShape>)],
    l_fs: &mut Vec<TopoShape>,
) {
    l_fs.clear();

    // `cxx:1950-1953`
    let mut f = f_ref.clone();
    f.set_orientation(Orientation::Forward);
    let to_reverse = f_ref.orientation().is_reversed();

    let reg = GeometryRegistry::global();
    let builder = TopoBuilder::new();
    // `cxx:1956-1973`
    for (wi, holes) in map_wlow {
        let mut ff = f.empty_copied();
        // `F.EmptyCopied()` on a face copies `BRep_TFace`'s surface and
        // tolerance (`BRep_TFace.cxx:36-43`).
        if let Some(geom) = reg.face_geom(&f) {
            reg.set_face(&ff, geom);
        }
        builder.add(&mut ff, wi);
        for wwi in holes {
            builder.add(&mut ff, wwi);
        }
        if to_reverse {
            ff.set_orientation(Orientation::Reversed);
        }
        l_fs.push(ff);
    }
}
