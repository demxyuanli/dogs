//! The descendant bookkeeping of `BRepFeat_Form`.
//!
//! Source: `BRepFeat_Form.cxx:1310-1336` (`UpdateDescendants(Gluer)`),
//! `cxx:1489-1508` (the file-local `Descendants`), `cxx:1512-1580`
//! (`UpdateDescendants(BooleanOperation, S, SkipFace)`).

use crate::abs::ShapeType;
use crate::brep_feat::algo_api::BooleanOperation;
use crate::brep_feat::builder::BRepFeatBuilder;
use crate::locope::{Gluer, ShapeSet};
use crate::shape::{Face, TopoShape};
use crate::topo_tools::{faces_of, is_same};

use super::BRepFeatForm;

/// `Descendants(S, theFB, mapF)` (`cxx:1489-1508`): every face of `S` is
/// looked up in the history of `theFB` and its images collected.
pub fn descendants(s: &TopoShape, fb: &BRepFeatBuilder, map_f: &mut ShapeSet) {
    map_f.clear();
    for fdsc in faces_of(s) {
        for v in fb.modified(&fdsc.0) {
            map_f.add(&v);
        }
    }
}

/// `UpdateDescendants(const LocOpe_Gluer&)` (`cxx:1310-1336`): every value of
/// `myMap` is replaced by the gluer descendants of its faces.
///
/// `Err` is the `Standard_TypeMismatch` of `TopoDS::Face` on a value that is
/// not a face.
pub fn update_descendants_gluer(form: &mut BRepFeatForm, g: &Gluer) -> Result<(), String> {
    let keys: Vec<TopoShape> = form.my_map.iter().map(|(k, _)| k.clone()).collect();
    for orig in keys {
        let list = form.my_map.find(&orig).cloned().unwrap_or_default();
        let mut newdsc = ShapeSet::new();
        for v in &list {
            let fdsc = Face::wrap(v.clone())
                .ok_or_else(|| "Standard_TypeMismatch: TopoDS::Face".to_string())?;
            for d in g.descendant_faces(&fdsc).map_err(|e| e.to_string())? {
                newdsc.add(d);
            }
        }
        let slot = form
            .my_map
            .change(&orig)
            .ok_or_else(|| "Standard_NoSuchObject: myMap".to_string())?;
        slot.clear();
        for k in newdsc.iter() {
            slot.push(k.clone());
        }
    }
    Ok(())
}

/// `UpdateDescendants(const BRepAlgoAPI_BooleanOperation&, S, SkipFace)`
/// (`cxx:1512-1580`).
///
/// A value that disappeared from `S` is replaced by the `Modified` images of
/// the Boolean; a value still present in `S` is kept as-is; only the results
/// that are faces of `S` stay in the map.
pub fn update_descendants_bop(
    form: &mut BRepFeatForm,
    bop: &BooleanOperation,
    s: &TopoShape,
    skip_face: bool,
) -> Result<(), String> {
    let s_faces = faces_of(s);
    let keys: Vec<TopoShape> = form.my_map.iter().map(|(k, _)| k.clone()).collect();
    for orig in keys {
        // `cxx:1524-1528`.
        if skip_face && orig.shape_type() == ShapeType::Face {
            continue;
        }
        // `cxx:1532-1535`: an empty list is seeded with the key itself.
        if let Some(slot) = form.my_map.change(&orig) {
            if slot.is_empty() {
                slot.push(orig.clone());
            }
        }
        let mut newdsc = ShapeSet::new();
        // `cxx:1537-1563`.
        let list = form.my_map.find(&orig).cloned().unwrap_or_default();
        for sh in &list {
            if sh.shape_type() != ShapeType::Face {
                continue;
            }
            if s_faces.iter().any(|f| is_same(&f.0, sh)) {
                // preserved
                newdsc.add(sh);
                continue;
            }
            for a_s in bop.modified(sh) {
                newdsc.add(&a_s);
            }
        }
        // `cxx:1564-1578`: only the images that belong to `S` stay.
        let slot = form
            .my_map
            .change(&orig)
            .ok_or_else(|| "Standard_NoSuchObject: myMap".to_string())?;
        slot.clear();
        for k in newdsc.iter() {
            if s_faces.iter().any(|f| is_same(&f.0, k)) {
                slot.push(k.clone());
            }
        }
    }
    Ok(())
}
