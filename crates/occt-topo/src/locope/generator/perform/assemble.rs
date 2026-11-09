//! `LocOpe_Generator.cxx:1092-1218`: the final assembly.
//!
//! The surviving base faces and the generated faces go into one loose face list
//! that `LocOpe_BuildShape` turns into shells and solids. The reference then
//! restores, on a single-solid answer, the orientation the first base face had:
//! when `BuildShape` came back with that face reversed, the solid is rebuilt
//! around a reversed shell (`cxx:1194-1213`).

use crate::abs::{Orientation, ShapeType};
use crate::builder::TopoBuilder;
use crate::locope::build_shape::BuildShape;
use crate::locope::shape_maps::ShapeValueMap;
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, is_same};

use super::Stage;

/// `cxx:1133-1218`.
pub(super) fn build(stage: &Stage, mod_shapes: &mut ShapeValueMap<Vec<TopoShape>>) -> Option<TopoShape> {
    let mut lfres: Vec<TopoShape> = Vec::new();
    let mut face_ref_ori: Option<TopoShape> = None;

    // `cxx:1137-1164`.
    for fac in explorer_shapes(&stage.shape, ShapeType::Face) {
        // `cxx:1141`.
        if stage.left.contains(&fac) {
            continue;
        }
        // `cxx:1143-1157`.
        if !stage.to_remove.contains(&fac) {
            lfres.push(fac.clone());
            if !mod_shapes.is_bound(&fac) {
                mod_shapes.bind(&fac, Vec::new());
            }
            if let Some(list) = mod_shapes.change(&fac) {
                list.push(fac.clone());
            }
            if face_ref_ori.is_none() {
                face_ref_ori = Some(fac.clone());
            }
        } else if let Some(first) = mod_shapes.find(&fac).and_then(|l| l.first().cloned()) {
            // `cxx:1158-1161`. A bound-but-empty list raises in the reference's
            // `First()`; the guard keeps such a face out of the result.
            lfres.push(first.oriented(fac.orientation()));
        }
    }

    // `cxx:1166-1187`.
    let orsolid = stage.shape.orientation();
    for fac in stage.g.oriented_faces() {
        if stage.to_remove.contains(fac) {
            continue;
        }
        if orsolid == Orientation::Forward {
            lfres.push(fac.clone());
        } else {
            lfres.push(fac.oriented(fac.orientation().reversed()));
        }
        if !mod_shapes.is_bound(fac) {
            mod_shapes.bind(fac, Vec::new());
        }
        if let Some(list) = mod_shapes.change(fac) {
            list.push(fac.clone());
        }
    }

    // `cxx:1189-1190`.
    let mut bs = BuildShape::new();
    let _ = bs.perform(&lfres);
    let my_res = bs.shape().cloned()?;

    // `cxx:1192-1214`: give a single solid back the orientation of the first
    // base face.
    if my_res.shape_type() == ShapeType::Solid {
        let found = explorer_shapes(&my_res, ShapeType::Face)
            .into_iter()
            .find(|f| face_ref_ori.as_ref().is_some_and(|r| is_same(f, r)));
        if let (Some(found), Some(reference)) = (found, face_ref_ori.clone()) {
            if found.orientation() != reference.orientation() {
                let shell = explorer_shapes(&my_res, ShapeType::Shell)
                    .into_iter()
                    .next()
                    .map(|s| s.oriented(s.orientation().reversed()));
                if let Some(shell) = shell {
                    // `B.MakeSolid(NewSol); B.Add(NewSol, shell)`: the new solid
                    // is FORWARD, so `myRes.Orientation` resets to FORWARD.
                    let new_solid = TopoBuilder::new().make_solid(&[crate::shape::Shell::wrap(shell)
                        .expect("a shell child of a solid")]);
                    return Some(new_solid.0);
                }
            }
        }
    }
    Some(my_res)
}
