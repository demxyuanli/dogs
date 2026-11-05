//! `IntPlanEdge` / `ComputeEdgeParameter` — the two file statics of
//! `FilletSurf_InternalBuilder.cxx`.
//!
//! Source: `ModelingAlgorithms/TKFillet/FilletSurf/FilletSurf_InternalBuilder.cxx:59-153`
//! (OCCT 8.0.0).
//!
//! `IntPlanEdge` intersects the elementary spine of one edge with the plane
//! through the guideline point, normal to the guideline tangent, and answers
//! the parameter of the intersection closest to the plane origin:
//!
//! ```text
//! Intersection.Perform(Ed, Plan);
//! for (iip = 1; iip <= nbp; iip++) {
//!   d = Point(iip).Pnt().Distance(Or);
//!   if (d < dist) { done = true; w = Point(iip).W(); dist = d; }
//! }
//! // then the two extremities of Ed are tested the same way, but only when
//! // they lie on the plane within tol3d
//! ```
//!
//! Divergences that are not heuristics:
//! - `IntCurveSurface_HInter::Perform(Ed, Plan)` is
//!   [`intcurvesurface::perform_curve_surface`], the port of the same class.
//!   That port needs a finite `(u0, v0, u1, v1)` window while OCCT's
//!   `GeomAdaptor_Surface(Geom_Plane)` is unbounded; [`unbounded_plane_window`]
//!   is the representable stand-in, exactly as in
//!   [`crate::chfi3d::builder_0_helpers::inter_plane_edge`] (the same adaptor
//!   pair, from `ChFi3d_InterPlaneEdge`),
//! - `ElSLib::Parameters(P, Pnt, u, v)` / `ElSLib::Value(u, v, P)` are
//!   [`slib::plane_parameters`] / [`slib::plane_value`],
//! - `Spine->ElSpine(ind)` is [`Spine::el_spine_of_edge_index`], and the
//!   `Standard_NoSuchObject` of an out-of-range `ind` is the `false` of
//!   [`compute_edge_parameter`].

use occt_core::elib::slib;
use occt_core::gp::{GpAx3, GpDir, GpPln};
use occt_geom::{Curve, GeomPlane};

use crate::chfi_ds::SpineHandle;
use crate::intcurvesurface::perform_curve_surface;

/// Stand-in for the infinite `U` / `V` bounds of `GeomAdaptor_Surface` on a
/// `Geom_Plane`. `1e100` keeps every `ComputeAppendPoint` window test inside
/// the window while leaving the Euclidean arithmetic untouched.
pub const UNBOUNDED_PLANE_WINDOW: f64 = 1.0e100;

/// The plane window used for an unbounded plane surface.
pub fn unbounded_plane_window() -> (f64, f64, f64, f64) {
    let r = UNBOUNDED_PLANE_WINDOW;
    (-r, r, -r, r)
}

/// `IntPlanEdge(Ed, P, w, tol3d)` (`cxx:72-138`).
///
/// `w` is the parameter on `ed` of the intersection closest to `P.Location()`;
/// `false` is OCCT's `done == false`, and `w` is then left untouched.
pub fn int_plan_edge(ed: &dyn Curve, p: &GpPln, w: &mut f64, tol3d: f64) -> bool {
    let f = ed.first_parameter();
    let l = ed.last_parameter();
    let op = p.location();
    let pln = GeomPlane::new(p.clone());

    let mut done = false;
    // `double dist = RealLast();` (`cxx:87`).
    let mut dist = f64::MAX;

    // `Intersection.Perform(Ed, Plan);` (`cxx:89`).
    if let Ok(res) = perform_curve_surface(ed, &pln, (f, l), unbounded_plane_window()) {
        // `nbp = Intersection.NbPoints(); for (iip = 1; iip <= nbp; iip++)`
        // (`cxx:93-107`).
        for ip in res.points() {
            let pint = ip.pnt();
            let d = pint.distance(&op);
            if d < dist {
                done = true;
                *w = ip.param();
                dist = d;
            }
        }
    }

    // `gp_Pnt pdeb = Ed->Value(f); gp_Pnt pfin = Ed->Value(l);` (`cxx:107-108`).
    let pdeb = ed.d0(f);
    let pfin = ed.d0(l);

    // Check if the extremities are not solution (`cxx:109-138`).
    let pos: GpAx3 = p.position();
    let (u, v) = slib::plane_parameters(&pos, &pdeb);
    let projdeb = slib::plane_value(p, u, v);
    if pdeb.distance(&projdeb) < tol3d {
        let d = pdeb.distance(&op);
        if d < dist {
            done = true;
            *w = f;
            dist = d;
        }
    }
    let (u, v) = slib::plane_parameters(&pos, &pfin);
    let projfin = slib::plane_value(p, u, v);
    if pfin.distance(&projfin) < tol3d {
        let d = pfin.distance(&op);
        if d < dist {
            done = true;
            *w = l;
            // `dist = d;` (`cxx:135`) is the last write of the C++ body and
            // nothing reads it afterwards; the branch is otherwise identical.
        }
    }
    done
}

/// `ComputeEdgeParameter(Spine, ind, pelsp, ped, tol3d)` (`cxx:140-153`).
///
/// The plane is the one through the guideline point `Guide->D1(pelsp)`,
/// normal to the guideline tangent there.
pub fn compute_edge_parameter(
    spine: &SpineHandle,
    ind: i32,
    pelsp: f64,
    ped: &mut f64,
    tol3d: f64,
) -> bool {
    // `occ::handle<ChFiDS_ElSpine> Guide = Spine->ElSpine(ind);` (`cxx:146`).
    let guide = {
        let Ok(guard) = spine.read() else {
            return false;
        };
        guard.el_spine_of_edge_index(ind)
    };
    let Some(guide) = guide else {
        return false;
    };
    // `gp_Pnt P; gp_Vec V; Guide->D1(pelsp, P, V);` (`cxx:147-149`).
    let (p, v) = guide.d1(pelsp);
    // `gp_Pln pln(P, V);` (`cxx:150`).
    let Ok(dir) = GpDir::new(v.x(), v.y(), v.z()) else {
        return false;
    };
    let Ok(ax3) = GpAx3::from_point_dir(p, dir) else {
        return false;
    };
    let pln = GpPln::new(ax3);

    // `occ::handle<BRepAdaptor_Curve> ed =
    // new BRepAdaptor_Curve(Spine->CurrentElementarySpine(ind));` (`cxx:151`).
    let ed = {
        let Ok(mut guard) = spine.write() else {
            return false;
        };
        guard.current_elementary_spine(ind)
    };
    let Some(ed) = ed else {
        return false;
    };
    int_plan_edge(&ed, &pln, ped, tol3d)
}
