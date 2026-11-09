//! File-local `PutPCurve(Edg, Fac)` (`LocOpe_WiresOnShape.cxx:722-910`).
//!
//! Projects the 3d curve of an edge onto a face and stores the p-curve, after
//! shifting it inside the face's UV box when the surface is periodic. Only the
//! 3d projection is reproduced here (the reference never calls
//! `ShapeConstruct_ProjectCurveOnSurface::PerformByCurve` in this helper).

use std::sync::Arc;

use occt_core::gp::{GpTrsf2d, GpVec2d};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom2d::curve::Curve2d;

use crate::brep_lib_same_parameter::same_parameter_edge_inplace;
use crate::brep_offset::update_edge_pcurve;
use crate::brep_tool::BRepTool;
use crate::brep_tools_uv::uv_bounds_face;
use crate::pcurve_full::{project_curve_on_surface_perform, ProjectorCache};
use crate::shape::{Edge, Face};
use crate::tgeometry::GeometryRegistry;

/// `gp_Geom2d_Curve::Translated(gp_Vec2d)` (`cxx:869`, `cxx:902`): a translated
/// copy, never an in-place move.
fn translated(curve: &Arc<dyn Curve2d>, dx: f64, dy: f64) -> Arc<dyn Curve2d> {
    let mut t = GpTrsf2d::identity();
    t.set_translation_vec(&GpVec2d::new(dx, dy));
    Arc::from(curve.transformed(&t))
}

/// `PutPCurve(Edg, Fac)` (`LocOpe_WiresOnShape.cxx:722-910`).
///
/// `the_umin` / `the_umax` / `the_vmin` / `the_vmax` are only loop counters in
/// the two period blocks - the reference never reads them afterwards - hence the
/// `unused_assignments` allowance.
#[allow(unused_assignments)]
pub(super) fn put_pcurve(edg: &Edge, fac: &Face) {
    // `cxx:726-737`: a rectangular trim is replaced by its basis surface, and a
    // plane needs no p-curve at all.
    let Some(mut s) = BRepTool::face_surface_world(fac) else {
        return;
    };
    if let Some(basis) = s.rectangular_trimmed_basis() {
        s = basis;
    }
    if s.gp_pln().is_some() {
        return;
    }

    // `cxx:741-742`.
    let (umin, umax, vmin, vmax) = uv_bounds_face(fac);

    // `cxx:745-764`.
    let (f, l) = BRepTool::edge_parameters(edg);
    if let Some((a_c2d, cf, cl)) = crate::boptools_2d::curve_on_surface_range(edg, fac) {
        let p2d = a_c2d.d0((cf + cl) * 0.5);
        let mut is_in = true;
        if p2d.x() < umin - PCONFUSION || p2d.x() > umax + PCONFUSION {
            is_in = false;
        }
        if p2d.y() < vmin - PCONFUSION || p2d.y() > vmax + PCONFUSION {
            is_in = false;
        }
        if is_in {
            return;
        }
    }

    // `cxx:766-775`: `BRep_Tool::Curve(Edg, Loc, f, l)` with the location folded
    // into the curve; `edge_curve_world` is exactly that, and `f`/`l` are the
    // edge range read just above.
    let Some(c) = BRepTool::edge_curve_world(edg) else {
        return;
    };
    // `cxx:777-780`: `Geom_TrimmedCurve(C, f, l)`. The port's trimmed curve
    // rescales its parameter to [0, 1] (`Geom_TrimmedCurve::D0`), unlike OCCT's,
    // so the wrapper is skipped and the projector below gets the explicit
    // [f, l] window: the same basis parameters the OCCT trimmed curve answers.
    //
    // `cxx:781`: `S = BRep_Tool::Surface(Fac)` - the world surface, not the
    // rectangular-trim basis read at `cxx:729`.
    let Some(s_world) = BRepTool::face_surface_world(fac) else {
        return;
    };

    // `cxx:783-792`: `TolFirst` / `TolLast` from the two edge vertices.
    let (v_first, v_last) = crate::topo_tools::edge_vertices(edg);
    let mut tol_first = -1.0;
    let mut tol_last = -1.0;
    if let Some(v1) = v_first.as_ref() {
        tol_first = BRepTool::vertex_tolerance(v1);
    }
    if let Some(v2) = v_last.as_ref() {
        tol_last = BRepTool::vertex_tolerance(v2);
    }

    // `cxx:794-799`: `ShapeConstruct_ProjectCurveOnSurface::Init(S, tol2d)` with
    // `tol2d = Precision::Confusion()`, then `Perform(C, f, l, C2d, TolFirst,
    // TolLast)`.
    const TOL2D: f64 = CONFUSION;
    let mut cache = ProjectorCache::default();
    let Some(c2d) = project_curve_on_surface_perform(
        c.as_ref(),
        s_world.as_ref(),
        f,
        l,
        TOL2D,
        tol_first,
        tol_last,
        &mut cache,
    ) else {
        return;
    };

    // `cxx:804-809`. Both points feed only the unported "internal vertices"
    // arm below (see the note there), so they stay named for the reference's
    // control flow without a reader.
    let pf = c2d.d0(f);
    let pl = c2d.d0(l);
    let _p_f = s_world.d0(pf.x(), pf.y());
    let _p_l = s_world.d0(pl.x(), pl.y());

    // `cxx:811-827`: `V1` / `V2` are re-read with `TopExp::FirstVertex` /
    // `LastVertex` (CumOri = true) and flipped on a REVERSED edge. They feed only
    // the arm below.
    //
    // `cxx:829-841`: "Handling of internal vertices". The reference enters this
    // arm only when `V2` is null and then still calls `BRep_Tool::Tolerance(V2)`
    // and `BRep_Tool::Pnt(V2)`, i.e. it dereferences the null vertex. That is
    // unreachable for an edge with two endpoints and undefined otherwise, so the
    // arm is UNPORTED rather than given an invented non-null substitute.

    let mut c2d = c2d;
    // `cxx:843-871`: the U period shift against `[Umin, Umax]`.
    if s_world.is_u_periodic() {
        let up = s_world.u_period();
        let tolu = PCONFUSION;
        let mut nbtra = 0_i32;
        let mut the_umin = pf.x().min(pl.x());
        let mut the_umax = pf.x().max(pl.x());
        if the_umin < umin - tolu {
            while the_umin < umin - tolu {
                the_umin += up;
                the_umax += up;
                nbtra += 1;
            }
        } else if the_umax > umax + tolu {
            while the_umax > umax + tolu {
                the_umax -= up;
                the_umin -= up;
                nbtra -= 1;
            }
        }
        if nbtra != 0 {
            c2d = translated(&c2d, nbtra as f64 * up, 0.0);
        }
    }

    // `cxx:873-904`: the V period shift against `[Vmin, Vmax]`.
    if s_world.is_v_periodic() {
        let vp = s_world.v_period();
        let tolv = PCONFUSION;
        let mut nbtra = 0_i32;
        let mut the_vmin = pf.y().min(pl.y());
        let mut the_vmax = pf.y().max(pl.y());
        if the_vmin < vmin - tolv {
            while the_vmin < vmin - tolv {
                the_vmin += vp;
                the_vmax += vp;
                nbtra += 1;
            }
        } else if the_vmax > vmax + tolv {
            while the_vmax > vmax + tolv {
                the_vmax -= vp;
                the_vmin -= vp;
                nbtra -= 1;
            }
        }
        if nbtra != 0 {
            c2d = translated(&c2d, 0.0, nbtra as f64 * vp);
        }
    }
    // `cxx:905-908`.
    update_edge_pcurve(edg, &c2d, fac, BRepTool::edge_tolerance(edg));
    GeometryRegistry::global().set_same_parameter(&edg.0, false);
    same_parameter_edge_inplace(edg, TOL2D);
}
