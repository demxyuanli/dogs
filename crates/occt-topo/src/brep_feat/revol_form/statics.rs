//! The two file-local helpers of `BRepFeat_MakeRevolutionForm`.
//!
//! `MajMap` (`cxx:1937-1981`) and `SetGluedFaces` (`cxx:1983-2031`) are the
//! revolution form's own versions: they key the map on the profile wires of
//! `FirstShape()` / `LastShape()` and on the basis edges, unlike the linear
//! form's helpers (`BRepFeat_MakeLinearForm.cxx:1292-1366`), which key it on
//! the basis faces.

use crate::abs::ShapeType;
use crate::locope::shape_maps::{ShapeDataMap, ShapeValueMap};
use crate::locope::RevolutionForm;
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, is_same};

/// `MajMap(theB, theP, theMap, theFShape, theLShape)` (`cxx:1937-1981`).
///
/// `theB` is `myPbase`, `theP` the finished `LocOpe_RevolutionForm`, `theMap`
/// `myMap`, and the two `TopoDS_Shape&` are `myFShape` / `myLShape`.
pub fn maj_map(
    the_b: &TopoShape,
    the_p: &RevolutionForm,
    the_map: &mut ShapeValueMap<Vec<TopoShape>>,
    the_f_shape: &mut Option<TopoShape>,
    the_l_shape: &mut Option<TopoShape>,
) {
    // `cxx:1947-1957`: the first wire of `FirstShape()` heads the map, followed
    // by every face of `FirstShape()`.
    if let Some(first) = the_p.first_shape() {
        if let Some(w) = explorer_shapes(first, ShapeType::Wire).first() {
            *the_f_shape = Some(w.clone());
            the_map.bind(w, explorer_shapes(first, ShapeType::Face));
        }
    }

    // `cxx:1959-1968`: the same for `LastShape()`.
    if let Some(last) = the_p.last_shape() {
        if let Some(w) = explorer_shapes(last, ShapeType::Wire).first() {
            *the_l_shape = Some(w.clone());
            the_map.bind(w, explorer_shapes(last, ShapeType::Face));
        }
    }

    // `cxx:1971-1979`: every basis edge not yet in the map answers the shapes
    // `LocOpe_RevolutionForm` generated from it. `theMap(K) = theP.Shapes(K)`
    // replaces the whole list, so a `bind` is the same operation.
    for e in explorer_shapes(the_b, ShapeType::Edge) {
        if !the_map.is_bound(&e) {
            let images = the_p.shapes(&e).map(|v| v.to_vec()).unwrap_or_default();
            the_map.bind(&e, images);
        }
    }
}

/// `SetGluedFaces(theSlmap, thePrism, SlidingMap, theMap)` (`cxx:1983-2031`).
///
/// `theSlmap` is `mySlface`, `thePrism` the finished `LocOpe_RevolutionForm`,
/// `SlidingMap` the face-to-faces map `Perform` builds and `theMap` `myGluedF`.
pub fn set_glued_faces(
    the_slmap: &ShapeValueMap<Vec<TopoShape>>,
    the_prism: &RevolutionForm,
    sliding_map: &ShapeValueMap<Vec<TopoShape>>,
    the_map: &mut ShapeDataMap,
) {
    // `cxx:1995`: an empty `theSlmap` leaves `myGluedF` untouched.
    if the_slmap.is_empty() {
        return;
    }
    for (fac, ledg) in the_slmap.iter() {
        for e in ledg {
            // `cxx:2004-2009`: the reference traces "Pb SetGluedFace" when the
            // prism made other than one shape from this edge and then reads
            // `First()` regardless.
            let gfac = the_prism.shapes(e).unwrap_or(&[]);
            for (ff, lfaces) in sliding_map.iter() {
                if lfaces.is_empty() {
                    continue;
                }
                let fff = &lfaces[0];
                if let Some(g) = gfac.first() {
                    if is_same(g, ff) {
                        the_map.bind(fff, fac);
                    }
                }
            }
        }
    }
}
