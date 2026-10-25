//! Port of `TopOpeBRepTool_tol` (`TopOpeBRepTool_tol.hxx:23-72`,
//! `TopOpeBRepTool_tol.cxx:41-311`): the face/face intersection tolerance
//! computation `FTOL_FaceTolerances` plus its three wrappers.
//!
//! `FBOX_GetHBoxTool()` (the global `TopOpeBRepTool_HBoxTool` in
//! `TopOpeBRepTool_box.cxx`) is passed in as `&HBoxTool`, as the rest of this
//! port does for `BoxSort`.

use occt_core::bnd::BndBox;
use occt_core::precision::CONFUSION;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face};
use crate::topo_tools::explorer_shapes;

use super::adaptor::SurfaceAdaptor;
use super::h_box_tool::HBoxTool;

/// `FTOL_FaceTolerances` (`cxx:41-219`).
///
/// `tol1` / `tol2` / `deflection` / `max_uv` are the by-reference outputs.
#[allow(clippy::too_many_arguments)]
pub fn ftol_face_tolerances(
    b1: &BndBox,
    b2: &BndBox,
    face1: &Face,
    face2: &Face,
    surface1: &SurfaceAdaptor,
    surface2: &SurfaceAdaptor,
    tol1: &mut f64,
    tol2: &mut f64,
    deflection: &mut f64,
    max_uv: &mut f64,
) {
    let a_tol_f1 = BRepTool::face_tolerance(face1);
    let a_tol_f2 = BRepTool::face_tolerance(face2);
    *tol1 = a_tol_f1 + a_tol_f2;
    *tol2 = *tol1;

    let box1_ok = !b1.is_open_xmin()
        && !b1.is_open_xmax()
        && !b1.is_open_ymin()
        && !b1.is_open_ymax()
        && !b1.is_open_zmin()
        && !b1.is_open_zmax()
        && !b1.is_void();
    let box2_ok = !b2.is_open_xmin()
        && !b2.is_open_xmax()
        && !b2.is_open_ymin()
        && !b2.is_open_ymax()
        && !b2.is_open_zmin()
        && !b2.is_open_zmax()
        && !b2.is_void();

    let mut dx;
    let mut dy;
    let mut dz;
    if box1_ok {
        let (x0, y0, z0, x1, y1, z1) = b1.get().unwrap_or((0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
        dx = x1 - x0;
        dy = y1 - y0;
        dz = z1 - z0;

        if box2_ok {
            let (x0, y0, z0, x1, y1, z1) = b2.get().unwrap_or((0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
            let (dx2, dy2, dz2) = (x1 - x0, y1 - y0, z1 - z0);
            if dx2 > dx {
                dx = dx2;
            }
            if dy2 > dy {
                dy = dy2;
            }
            if dz2 > dz {
                dz = dz2;
            }
        }
    } else if box2_ok {
        let (x0, y0, z0, x1, y1, z1) = b2.get().unwrap_or((0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
        dx = x1 - x0;
        dy = y1 - y0;
        dz = z1 - z0;
    } else {
        dx = 1.0;
        dy = 1.0;
        dz = 1.0;
    }

    if dx < dy {
        dx = dy;
    }
    if dx < dz {
        dx = dz;
    }
    if dx > 1000000.0 {
        dx = 1000000.0;
    }

    let mut tolef1 = CONFUSION;
    for shape in explorer_shapes(&face1.0, ShapeType::Edge) {
        let Some(edge) = Edge::wrap(shape) else {
            continue;
        };
        let tole = BRepTool::edge_tolerance(&edge);
        if tole > tolef1 {
            tolef1 = tole;
        }
    }
    let mut tolef2 = CONFUSION;
    for shape in explorer_shapes(&face2.0, ShapeType::Edge) {
        let Some(edge) = Edge::wrap(shape) else {
            continue;
        };
        let tole = BRepTool::edge_tolerance(&edge);
        if tole > tolef2 {
            tolef2 = tole;
        }
    }
    let mut _tolef = tolef1;
    if tolef2 > _tolef {
        _tolef = tolef2;
    }

    *deflection = 0.01;
    *max_uv = 0.01;
    *deflection *= dx;

    let mut mdeflection = *deflection;
    let mut mmaxuv = 0.01;

    let mu0 = surface1.first_u_parameter();
    let mu1 = surface1.last_u_parameter();
    let mv0 = surface1.first_v_parameter();
    let mv1 = surface1.last_v_parameter();
    let dp1 = mu1 - mu0;
    let dp2 = mv1 - mv0;
    let mut dp = if dp2 < dp1 { dp2 } else { dp1 };

    let mu0 = surface2.first_u_parameter();
    let mu1 = surface2.last_u_parameter();
    let mv0 = surface2.first_v_parameter();
    let mv1 = surface2.last_v_parameter();
    let mut dp1 = mu1 - mu0;
    let dp2 = mv1 - mv0;
    if dp2 > dp1 {
        dp1 = dp2;
    }
    if dp1 < dp {
        dp = dp1;
    }

    let _ = dp * 0.01;

    if mmaxuv < 1e-3 {
        mmaxuv = 1e-3;
    }
    if mdeflection < 1e-3 {
        mdeflection = 1e-3;
    }

    if mmaxuv > 0.01 {
        mmaxuv = 0.01;
    }
    if mdeflection > 0.1 {
        mdeflection = 0.1;
    }

    *deflection = mdeflection;
    *max_uv = mmaxuv;
}

/// `FTOL_FaceTolerances3d(F1, F2, Tol)` (`cxx:221-266`): the boxing tool
/// supplies both bounding boxes, a missing box falling back to
/// `Bnd_Box::Update(0,0,0,1,1,1)`.
pub fn ftol_face_tolerances3d_faces(
    hbt: &HBoxTool,
    face1: &Face,
    face2: &Face,
    tol: &mut f64,
) {
    let mut b1 = BndBox::new();
    let mut b2 = BndBox::new();
    if hbt.has_box(&face1.0) {
        b1 = hbt.box_of(&face1.0).cloned().unwrap_or_default();
    } else {
        b1.update(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    }
    if hbt.has_box(&face2.0) {
        b2 = hbt.box_of(&face2.0).cloned().unwrap_or_default();
    } else {
        b2.update(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    }

    let surface1 = SurfaceAdaptor::of_face(face1);
    let surface2 = SurfaceAdaptor::of_face(face2);
    let mut deflection = 0.01;
    let mut max_uv = 0.01;
    let mut tol1 = 0.0;
    let mut tol2 = 0.0;
    ftol_face_tolerances(
        &b1,
        &b2,
        face1,
        face2,
        &surface1,
        &surface2,
        &mut tol1,
        &mut tol2,
        &mut deflection,
        &mut max_uv,
    );
    tol1 = if tol1 > 1.0e-4 { 1.0e-4 } else { tol1 };
    tol2 = if tol2 > 1.0e-4 { 1.0e-4 } else { tol2 };
    *tol = tol1.max(tol2);
}

/// `FTOL_FaceTolerances3d(B1, B2, F1, F2, S1, S2, Tol1, Tol2, Deflection, MaxUV)`
/// (`cxx:268-286`): the same body as [`ftol_face_tolerances`], kept as a named
/// entry point because the intersection code calls it by that name.
#[allow(clippy::too_many_arguments)]
pub fn ftol_face_tolerances3d_boxes(
    b1: &BndBox,
    b2: &BndBox,
    face1: &Face,
    face2: &Face,
    surface1: &SurfaceAdaptor,
    surface2: &SurfaceAdaptor,
    tol1: &mut f64,
    tol2: &mut f64,
    deflection: &mut f64,
    max_uv: &mut f64,
) {
    ftol_face_tolerances(
        b1, b2, face1, face2, surface1, surface2, tol1, tol2, deflection, max_uv,
    );
}

/// `FTOL_FaceTolerances2d(B1, B2, F1, F2, S1, S2, Tol1, Tol2)` (`cxx:288-311`).
#[allow(clippy::too_many_arguments)]
pub fn ftol_face_tolerances2d(
    b1: &BndBox,
    b2: &BndBox,
    face1: &Face,
    face2: &Face,
    surface1: &SurfaceAdaptor,
    surface2: &SurfaceAdaptor,
    tol1: &mut f64,
    tol2: &mut f64,
) {
    let mut bid_deflection = 0.01;
    let mut bid_max_uv = 0.01;
    ftol_face_tolerances(
        b1,
        b2,
        face1,
        face2,
        surface1,
        surface2,
        tol1,
        tol2,
        &mut bid_deflection,
        &mut bid_max_uv,
    );
}
