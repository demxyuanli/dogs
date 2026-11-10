use std::f64::consts::PI;

use super::*;

use occt_core::gp::{
    GpAx2d, GpDir2d, GpPnt2d, GpTrsf2d, GpVec2d,
};
use occt_core::precision::CONFUSION;

use occt_topo::brep_tool::BRepTool;
use occt_topo::shape::TopoShape;

use super::super::iges_data_iges_model::EntityData;
use super::super::iges_to_brep_curve_and_surface::Session;

/// Parametric map applied to a 142 or 141 UV curve.
/// `ParamSurface` (`IGESToBRep_TopoSurface.cxx:1846-2070`).
///
/// Entity 120 stored as `Geom_SurfaceOfRevolution` does not enter the
/// cylinder/cone/torus/sphere shift (`cxx:1859-1862`). `TheULength` is reset
/// to 1 at `TransferTopoSurface` (`cxx:191`) and replaced by the line length
/// only for a line generatrix (`cxx:732`).
pub(in crate::iges_read) fn uv_map(
    session: &Session<'_>,
    basis: usize,
    face: &TopoShape,
) -> Result<(GpTrsf2d, f64), String> {
    let mut trans = GpTrsf2d::identity();
    let mut cscale = session.u_length;
    let mut uscale = 1.0;
    let surf = BRepTool::face_surface_world(&occt_topo::shape::Face(face.clone())).ok_or_else(|| {
        "BRep_Tool::Surface: IGESToBRep_TopoSurface.cxx:1851".to_string()
    })?;
    // One rectangular unwrap (`cxx:1853-1857`). A nested trim stays a trim.
    let surf = surf.rectangular_trimmed_basis().unwrap_or(surf);
    let (mut paramu, mut paramv) = elementary_param_shift(session, basis, face, surf.as_ref());
    if paramu.abs() <= CONFUSION {
        paramu = 0.0;
    }
    if paramv.abs() <= CONFUSION {
        paramv = 0.0;
    }
    // Offset unwrap is one level (`cxx:1988-1991`).
    let mut isrf = basis;
    if session.model.entity(isrf)?.type_number == 140 {
        if let EntityData::OffsetSurface { surface, .. } = &session.model.entity(isrf)?.data {
            isrf = *surface;
        }
    }
    let isrf_type = session.model.entity(isrf)?.type_number;
    if isrf_type == 120 {
        if let EntityData::Revolution { end_angle, .. } = &session.model.entity(isrf)?.data {
            let mut tmp = GpTrsf2d::identity();
            tmp.set_translation_vec(&GpVec2d::new(0.0, -2.0 * PI));
            trans.pre_multiply(&tmp);
            tmp = GpTrsf2d::identity();
            tmp.set_mirror_ax2d(&GpAx2d::new(
                GpPnt2d::new(0.0, 0.0),
                GpDir2d::new(1.0, 0.0).expect("ox"),
            ));
            trans.pre_multiply(&tmp);
            tmp = GpTrsf2d::identity();
            tmp.set_mirror_ax2d(&GpAx2d::new(
                GpPnt2d::new(0.0, 0.0),
                GpDir2d::new(1.0, 1.0).expect("diag"),
            ));
            trans.pre_multiply(&tmp);
            uscale = 1.0 / cscale;
            paramu = -(2.0 * PI - *end_angle);
        }
    } else {
        paramu = 0.0; // cxx:2010-2012
    }
    if isrf_type == 118 {
        uscale = 1.0 / cscale; // cxx:2015-2018
    }
    if isrf_type == 122 {
        // `Surf->Bounds` after the rectangular unwrap (`cxx:2021-2032`).
        let (umin, umax) = surf.u_range();
        let uln = (umax - umin).abs();
        uscale = uln / cscale;
        paramu = umin / uln;
    }
    if isrf_type == 192 || isrf_type == 194 {
        uscale = PI / 180.0; // cxx:2035-2038
    }
    if isrf_type == 196 {
        cscale = PI / 180.0; // cxx:2041-2044
        uscale = 1.0;
    }
    if isrf_type == 198 {
        let mut tmp = GpTrsf2d::identity();
        tmp.set_translation_vec(&GpVec2d::new(0.0, -360.0));
        trans.pre_multiply(&tmp);
        tmp = GpTrsf2d::identity();
        tmp.set_mirror_ax2d(&GpAx2d::new(
            GpPnt2d::new(0.0, 0.0),
            GpDir2d::new(1.0, 0.0).expect("ox"),
        ));
        trans.pre_multiply(&tmp);
        tmp = GpTrsf2d::identity();
        tmp.set_mirror_ax2d(&GpAx2d::new(
            GpPnt2d::new(0.0, 0.0),
            GpDir2d::new(1.0, 1.0).expect("diag"),
        ));
        trans.pre_multiply(&tmp);
        if paramv > 0.0 {
            paramv = paramv * 180.0 / PI; // cxx:2056-2058
        }
        cscale = PI / 180.0;
        uscale = 1.0;
    }
    let mut tmp = GpTrsf2d::identity();
    tmp.set_translation_pnts(&GpPnt2d::new(0.0, 0.0), &GpPnt2d::new(paramu, paramv));
    trans.pre_multiply(&tmp);
    tmp = GpTrsf2d::identity();
    tmp.set_scale(&GpPnt2d::new(0.0, 0.0), cscale)
        .map_err(|e| e.to_string())?;
    trans.pre_multiply(&tmp);
    Ok((trans, uscale))
}
