//! The two file-local helpers of `BRepFeat_MakeLinearForm.cxx`: `MajMap`
//! (`cxx:1292-1334`) and `SetGluedFaces` (`cxx:1336-1366`).

use crate::abs::ShapeType;
use crate::locope::shape_maps::{ShapeDataMap, ShapeValueMap};
use crate::locope::LinearForm;
use crate::shape::TopoShape;
use crate::topo_tools::explorer_shapes;

/// `MajMap(theB, theP, theMap, theFShape, theLShape)` (`cxx:1292-1334`).
///
/// The first / last wire of the prism become the first / last shape, and their
/// face children are their descendants. Every edge of the basis that the map
/// does not know yet is bound to the prism's images of that edge.
pub fn maj_map(
    base: &TopoShape,
    form: &LinearForm,
    map: &mut ShapeValueMap<Vec<TopoShape>>,
    f_shape: &mut Option<TopoShape>,
    l_shape: &mut Option<TopoShape>,
) {
    // `cxx:1300-1310`: the first wire and its faces.
    if let Some(first) = form.first_shape() {
        if let Some(w) = explorer_shapes(first, ShapeType::Wire).first() {
            *f_shape = Some(w.clone());
            map.bind(w, explorer_shapes(first, ShapeType::Face));
        }
    }

    // `cxx:1312-1323`: the last wire and its faces.
    if let Some(last) = form.last_shape() {
        if let Some(w) = explorer_shapes(last, ShapeType::Wire).first() {
            *l_shape = Some(w.clone());
            map.bind(w, explorer_shapes(last, ShapeType::Face));
        }
    }

    // `cxx:1325-1333`.
    for e in explorer_shapes(base, ShapeType::Edge) {
        if !map.is_bound(&e) {
            let images = form.shapes(&e).map(|v| v.to_vec()).unwrap_or_default();
            map.bind(&e, images);
        }
    }
}

/// `SetGluedFaces(theSlmap, thePrism, theMap)` (`cxx:1336-1366`): every
/// section edge of a sliding face is glued to the single prism face it
/// generated.
pub fn set_glued_faces(
    slmap: &ShapeValueMap<Vec<TopoShape>>,
    form: &LinearForm,
    map: &mut ShapeDataMap,
) {
    for (fac, ledg) in slmap.iter() {
        for e in ledg {
            // `cxx:1356-1362`: the reference prints on a size other than one
            // and then dereferences `First()` regardless.
            let gfac = form.shapes(e).unwrap_or(&[]);
            if let Some(first) = gfac.first() {
                map.bind(first, fac);
            }
        }
    }
}
