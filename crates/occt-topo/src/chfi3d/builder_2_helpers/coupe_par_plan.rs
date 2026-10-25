//! `ChFi3d_CoupeParPlan` — cut the fillet surface by the plane through the two
//! common points and their edge tangents.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_2.cxx:99-188`
//! (OCCT 8.0.0).
//!
//! The plane carries `P1`, `D12` and `T1`; when `T2` is within `0.01` of that
//! plane the fillet surface is intersected with it and the resulting
//! intersection curve plus its p-curve on the fillet become the stripe
//! extremity. Otherwise `plane` answers `false` and the caller falls back to
//! another extremity construction.
//!
//! Translation notes:
//! - `Pdeb` / `Pfin` are the `NCollection_Array1<double>(1, 4)` pair OCCT
//!   fills with the `UV` of the fillet at the extremity and the `UV` of the
//!   projected point on the plane; the port passes `[f64; 4]`,
//! - `gp_Dir` construction throws on a null / parallel input. OCCT would abort
//!   the whole build there; the port answers `plane = false`, the same
//!   outcome as the non-planar branch, so the caller keeps its fallback,
//! - `Extrema_ExtPS(..., Extrema_ExtFlag_MIN)` is
//!   [`ExtPs`] with `set_flag(Min)` / `set_algo(Grad)` set before
//!   `initialize`, the order of the OCCT constructor (`Extrema_ExtPS.cxx:192-193`).

use std::sync::Arc;

use occt_core::gp::{GpAx3, GpDir, GpPln, GpPnt2d, GpVec};
use occt_geom::curve::Curve;
use occt_geom::extrema_surf::{ExtPs, ExtremaExtAlgo, ExtremaExtFlag};
use occt_geom::{GeomPlane, Surface};
use occt_geom2d::curve::Curve2d;

use crate::chfi3d::builder_0_helpers::{ch_fi3d_compute_curves, ComputeCurves};
use crate::chfi_ds::CommonPoint;
use crate::shape::Edge;
use crate::topopebrep::tool::adaptor::{CurveAdaptor, SurfaceAdaptor};

/// The out-parameters of [`ch_fi3d_coupe_par_plan`].
#[derive(Default, Clone)]
pub struct CoupeParPlan {
    /// `C3d`.
    pub c3d: Option<Arc<dyn Curve>>,
    /// `pcurve` — the same intersection seen from the fillet surface.
    pub pcurve: Option<Arc<dyn Curve2d>>,
    /// `tolreached`.
    pub tolreached: f64,
    /// `Pardeb`.
    pub pardeb: f64,
    /// `Parfin`.
    pub parfin: f64,
    /// `plane`: `false` when the analytic cut does not apply.
    pub plane: bool,
}

/// `ChFi3d_CoupeParPlan(compoint1, compoint2, HConge, UV1, UV2, tol3d, tol2d,
/// C3d, pcurve, tolreached, Pardeb, Parfin, plane)` (`cxx:99-188`).
#[allow(clippy::too_many_arguments)]
pub fn ch_fi3d_coupe_par_plan(
    compoint1: &CommonPoint,
    compoint2: &CommonPoint,
    hconge: &SurfaceAdaptor,
    uv1: &GpPnt2d,
    uv2: &GpPnt2d,
    tol3d: f64,
    tol2d: f64,
    out: &mut CoupeParPlan,
) {
    // `plane = true;` (`cxx:115`).
    out.plane = true;
    // `if (compoint1.IsOnArc() && compoint2.IsOnArc())` (`cxx:116`).
    if !(compoint1.is_on_arc() && compoint2.is_on_arc()) {
        out.plane = false;
        return;
    }
    let (Some(arc1), Some(arc2)) = (compoint1.arc(), compoint2.arc()) else {
        out.plane = false;
        return;
    };
    let (Some(edge1), Some(edge2)) = (Edge::wrap(arc1.clone()), Edge::wrap(arc2.clone())) else {
        out.plane = false;
        return;
    };

    // `BRepAdaptor_Curve BCurv1(compoint1.Arc()); ... parE1 =
    // compoint1.ParameterOnArc(); BCurv1.D1(parE1, P1, t1);` (`cxx:118-123`).
    let bcurv1 = CurveAdaptor::of_edge(&edge1);
    let bcurv2 = CurveAdaptor::of_edge(&edge2);
    let (Some(par_e1), Some(par_e2)) = (
        compoint1.parameter_on_arc(),
        compoint2.parameter_on_arc(),
    ) else {
        out.plane = false;
        return;
    };
    let (Some((p1, t1)), Some((p2, t2))) = (bcurv1.d1(par_e1), bcurv2.d1(par_e2)) else {
        out.plane = false;
        return;
    };

    // `gp_Dir tgt1(t1); gp_Dir tgt2(t2); gp_Vec v12(P2 - P1); gp_Dir d12(v12);
    // gp_Dir nor = tgt1.Crossed(d12);` (`cxx:126-129`).
    let Ok(tgt1) = GpDir::from_vec(&t1) else {
        out.plane = false;
        return;
    };
    let Ok(tgt2) = GpDir::from_vec(&t2) else {
        out.plane = false;
        return;
    };
    let v12 = GpVec::new(p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z());
    let Ok(d12) = GpDir::from_vec(&v12) else {
        out.plane = false;
        return;
    };
    let Ok(nor) = tgt1.crossed(&d12) else {
        out.plane = false;
        return;
    };

    // `occ::handle<Geom_Plane> Plan = new Geom_Plane(P1, nor);` (`cxx:130`).
    let Ok(ax3) = GpAx3::from_point_dir(p1, nor) else {
        out.plane = false;
        return;
    };
    let plan = GeomPlane::new(GpPln::new(ax3));

    // `scal = std::abs(nor.Dot(tgt2)); if (scal < 0.01)` (`cxx:131-134`).
    let scal = nor.dot(&tgt2).abs();
    if scal >= 0.01 {
        // `else { plane = false; }` (`cxx:178-181`).
        out.plane = false;
        return;
    }

    // `occ::handle<GeomAdaptor_Surface> HPlan = new GeomAdaptor_Surface(Plan);`
    // (`cxx:136`).
    let surface: Arc<dyn Surface> = Arc::new(plan);
    let hplan = SurfaceAdaptor::of_surface(surface.clone());

    // `GeomAdaptor_Surface AS(Plan); Extrema_ExtPS anExtPS(P1, AS, ..., 1.e-3,
    // 1.e-3, Extrema_ExtFlag_MIN); anExtPS.Point(1).Parameter(u1, v1);`
    // (`cxx:138-149`).
    let mut ext_ps = ExtPs::new();
    ext_ps.set_flag(ExtremaExtFlag::Min);
    ext_ps.set_algo(ExtremaExtAlgo::Grad);
    ext_ps.initialize(
        surface.as_ref(),
        hplan.first_u_parameter(),
        hplan.last_u_parameter(),
        hplan.first_v_parameter(),
        hplan.last_v_parameter(),
        1.0e-3,
        1.0e-3,
    );
    ext_ps.perform(&p1);
    // OCCT dereferences `Point(1)` unconditionally; a failed search raises
    // `StdFail_NotDone` there. The port answers the non-planar branch.
    if !ext_ps.is_done() || ext_ps.nb_ext() < 1 {
        out.plane = false;
        return;
    }
    let (u1, v1, _) = ext_ps.point(1);

    // `Pdeb(1) = UV1.X(); Pdeb(2) = UV1.Y(); Pdeb(3) = u1; Pdeb(4) = v1;`
    // (`cxx:150-154`).
    let pdeb = [uv1.x(), uv1.y(), u1, v1];

    // `anExtPS.Perform(P2); anExtPS.Point(1).Parameter(u1, v1);` (`cxx:155-157`).
    ext_ps.perform(&p2);
    if !ext_ps.is_done() || ext_ps.nb_ext() < 1 {
        out.plane = false;
        return;
    }
    let (u1, v1, _) = ext_ps.point(1);
    // `Pfin(1) = UV2.X(); Pfin(2) = UV2.Y(); Pfin(3) = u1; Pfin(4) = v1;`
    // (`cxx:157-159`).
    let pfin = [uv2.x(), uv2.y(), u1, v1];

    // `if (ChFi3d_ComputeCurves(HConge, HPlan, Pdeb, Pfin, C3d, pcurve,
    // C2dint2, tol3d, tol2d, tolreached))` (`cxx:160-169`).
    let mut curves = ComputeCurves::default();
    if ch_fi3d_compute_curves(hconge, &hplan, pdeb, pfin, &mut curves, tol3d, tol2d, false) {
        // `Pardeb = C3d->FirstParameter(); Parfin = C3d->LastParameter();`
        // (`cxx:170-171`).
        out.pardeb = curves.c3d.as_ref().map(|c| c.first_parameter()).unwrap_or(0.0);
        out.parfin = curves.c3d.as_ref().map(|c| c.last_parameter()).unwrap_or(0.0);
        out.c3d = curves.c3d.clone();
        out.pcurve = curves.pc1.clone();
        out.tolreached = curves.tolreached;
    } else {
        // `else { plane = false; }` (`cxx:173-176`).
        out.plane = false;
    }
}
