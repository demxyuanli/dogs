//! `BRepFeat_Form::TransformShapeFU` (`cxx:1378-1480`): limits the `From` /
//! `Until` shape when it is a single infinite face.

use crate::abs::ShapeType;
use crate::brep_feat::statics::face_until;
use crate::brep_tool::BRepTool;
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, faces_of};

use super::BRepFeatForm;

/// `TransformShapeFU(flag)` (`cxx:1378-1480`).
///
/// `flag` 0 is `mySFrom` and 1 is `mySUntil`; anything else answers `false`
/// (`cxx:1395-1398`). The returned flag is the reference's `Trf`: true when the
/// single face was infinite and [`face_until`] trimmed it.
pub fn transform_shape_fu(form: &mut BRepFeatForm, flag: i32) -> bool {
    let mut trf = false;

    // `cxx:1385-1398`.
    let shapefu = match flag {
        0 => form.my_sfrom.clone(),
        1 => form.my_suntil.clone(),
        _ => return trf,
    };
    // `cxx:1400-1407`: no face means an invalid shape. The port has no null
    // `TopoDS_Shape`, so an unset member is the empty group of the reference.
    let Some(shapefu) = shapefu else {
        return trf;
    };

    let faces = faces_of(&shapefu);
    if faces.is_empty() {
        return trf;
    }

    if faces.len() == 1 {
        // `cxx:1411-1457`: the only face. Is it infinite?
        let mut fac = faces[0].clone();

        let surface = BRepTool::face_surface(&fac);
        let surface = match surface {
            Some(s) => s,
            // `cxx:1414-1416` dereferences the null handle.
            None => return trf,
        };
        let surface = surface.rectangular_trimmed_basis().unwrap_or(surface);

        // `cxx:1422-1438`: only the three analytic families can be trimmed.
        if surface.gp_pln().is_some()
            || surface.gp_cylinder().is_some()
            || surface.gp_cone().is_some()
        {
            if explorer_shapes(&fac.0, ShapeType::Wire).is_empty() {
                trf = true;
            } else {
                trf = BRepTool::natural_restriction(&fac);
            }
        }
        // `cxx:1437-1439`: `FaceUntil` rewrites `fac` in place. It nullifies
        // the face only for the non-analytic families, which `trf` excludes.
        if trf {
            let base = form
                .my_sbase
                .clone()
                .unwrap_or_else(|| TopoShape::new(ShapeType::Compound));
            if let Some(trimmed) = face_until(&base, &fac) {
                fac = trimmed;
            }
        }

        // `cxx:1441-1456`.
        if flag == 0 {
            if let Some(key) = form.my_sfrom.clone() {
                form.my_map.bind(&key, vec![fac.0.clone()]);
            }
            form.my_sfrom = Some(fac.0.clone());
        } else {
            if let Some(key) = form.my_suntil.clone() {
                form.my_map.bind(&key, vec![fac.0.clone()]);
            }
            form.my_suntil = Some(fac.0.clone());
        }
    } else {
        // `cxx:1459-1468`.
        for fac in faces {
            form.my_map.bind(&fac.0, vec![fac.0.clone()]);
        }
    }

    trf
}
