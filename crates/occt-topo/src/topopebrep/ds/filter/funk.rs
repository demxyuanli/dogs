//! `TopOpeBRepDS_funk.cxx`: the 3d state of a point relative to a face, used to
//! resolve `(UNKNOWN, UNKNOWN)` face transitions.
//!
//! Sources (OCCT 8.0.0):
//! - `FUN_staPinF3d` (`TopOpeBRepDS_funk.cxx:44-72`)
//! - `FUN_UNKFstasta` (`TopOpeBRepDS_funk.cxx:74-279`)
//!
//! Two points at `0.41237118973` and its symmetric partners of the edge's
//! parameter are classified on `FS`, either through the caller's
//! `TopOpeBRepTool_ShapeClassifier` (`pClassif`) or, when there is none, through
//! `FUN_staPinF3d` (which classifies against the face restriction via
//! `BRepExtrema_ExtPF`).

use occt_core::gp::{GpDir, GpPnt, GpVec};

use crate::abs::{Orientation, ShapeType, State};
use crate::brep_tool::BRepTool;
use crate::brep_tools_uv;
use crate::boptools_2d;
use crate::shape::{Edge, Face};
use crate::topexp::Explorer;
use crate::topopebrep::tool::shape_classifier::ShapeClassifier;
use crate::topopebrep::tool::topology::{fun_tool_ng_geom_f, fun_tool_value_uv};

/// `FUN_staPinF3d(P, F)` (`cxx:44-73`).
///
/// `UNKNOWN` is both the initial value (`cxx:49`) and the `!ok` arm of
/// `FUN_tool_projPonboundedF` (`cxx:51-54`).
pub fn fun_sta_pin_f_3d(p: &GpPnt, f: &Face) -> State {
    let mut st = State::Unknown;
    // `cxx:47-56`
    let Some((pu, pv, d)) = crate::topopebrep::tool::project::fun_tool_proj_pon_bounded_f(p, f) else {
        return st;
    };
    let uv = occt_core::gp::GpPnt2d::new(pu, pv);
    // `cxx:57-60`
    let tol_f = BRepTool::face_tolerance(f);
    if d < tol_f {
        // `TopAbs_IN` (the C++ comments out `TopAbs_ON`).
        return State::In;
    }
    // `cxx:61-62`
    let Some(p_f) = fun_tool_value_uv(&uv, f) else {
        return st;
    };
    // `cxx:63-67`
    let Some(mut nt_f) = fun_tool_ng_geom_f(&uv, f) else {
        return st;
    };
    if f.0.orientation() == Orientation::Reversed {
        nt_f.reverse();
    }
    // `cxx:68-70`: `gp_Dir PpF(gp_Vec(P, pF))`, which raises on a zero vector.
    let Ok(pp_f) = GpDir::from_vec(&GpVec::from_pnts(p, &p_f)) else {
        return st;
    };
    let dot = nt_f.dot(&pp_f);
    st = if dot > 0.0 { State::In } else { State::Out };
    st
}

/// Whether `e` is one of `f`'s edges (`TopExp_Explorer` over `TopAbs_EDGE`).
fn edge_of_face(e: &Edge, f: &Face) -> bool {
    let mut ex = Explorer::new(&f.0, ShapeType::Edge);
    while ex.more() {
        let cur = ex.current().clone();
        ex.next();
        if crate::topo_tools::is_same(&cur, &e.0) {
            return true;
        }
    }
    false
}

/// `FUN_UNKFstasta(FF, FS, EE, EEofFF, stateb, statea, pClassif)`
/// (`cxx:75-279`).
///
/// `ee_of_ff` is OCCT's `EEofFF`: whether `EE` is an edge of `FF`.
/// `p_classif` is the classifier `FUN_resolveFUNKNOWN` passes when the face
/// pair did not produce a transition; it is mutated (`StateP3DReference`) and
/// reused across interferences, matching OCCT.
pub fn fun_unk_fstasta(
    ff: &Face,
    fs: &Face,
    ee: &Edge,
    ee_of_ff: bool,
    p_classif: Option<&mut ShapeClassifier>,
) -> (State, State) {
    let mut stateb = State::Unknown;
    let mut statea = State::Unknown;

    let (fe, le) = BRepTool::edge_parameters(ee);
    // NYI (`cxx:88-91`): get the points from the 2d curve.
    let Some(cee) = BRepTool::edge_curve_world(ee) else {
        return (stateb, statea);
    };
    let Some(sff) = BRepTool::face_surface_world(ff) else {
        return (stateb, statea);
    };

    let tt_e = 0.41237118973;
    let par_e = (1.0 - tt_e) * fe + tt_e * le;
    let (pe, ve) = cee.d1(par_e);

    // `GeomAPI_ProjectPointOnSurf PonS(PE, SFF)` (`cxx:100-110`).
    let Some(proj) =
        occt_geom::geom_api::project_point_on_surface(sff.as_ref(), &pe, occt_core::precision::CONFUSION)
    else {
        return (stateb, statea);
    };
    let (u, v) = (proj.u, proj.v);
    let (_puv, d1u, d1v) = sff.d1(u, v);
    let n = d1u.crossed(&d1v);

    // `BRepTools::UVBounds(FF, ...)` (`cxx:118-119`).
    let (fu_min, fu_max, fv_min, fv_max) = brep_tools_uv::uv_bounds_face(ff);

    // `BRepTools::UVBounds(FF, EE, ...)` (`cxx:121-181`).
    let (eu_min, eu_max, ev_min, ev_max) = if ee_of_ff {
        brep_tools_uv::uv_bounds_face_edge(ff, ee)
    } else {
        // `EE` is not an edge of `FF`, i.e. it is an edge of `FS`.
        //
        // `cxx:128-131` (the `CEE.IsNull()` arm building a
        // `ProjLib_ProjectedCurve`) is unreachable: `CEE` was rejected as null
        // at `cxx:88-91`.
        let pc = if edge_of_face(ee, ff) {
            // `cxx:146-157`: `BRep_Tool::CurveOnSurface(EE, FF, fE, lE)`.
            match boptools_2d::curve_on_surface(ee, ff) {
                Some(pc) => pc,
                None => return (stateb, statea),
            }
        } else {
            // `cxx:158-161`: `GeomProjLib::Curve2d(CEE, fE, lE, SFF)`, whose port
            // stand-in is `BOPTools_AlgoTools2D::Make2D` (the port's existing
            // replacement for `GeomProjLib`, see `boptools_2d::make_2d`).
            match boptools_2d::make_2d(ee, ff) {
                Ok(pc) => pc,
                Err(_) => return (stateb, statea),
            }
        };
        // `BRepTools::UVBounds(FFx, EEx, ...)` (`cxx:168-181`): the temporary
        // face carries `SFF` itself, so its bounds are those of `FF` and only
        // the p-curve box differs.
        let mut box2d = occt_core::bnd::BndBox2d::new();
        brep_tools_uv::add_uv_bounds_pcurve(ff, pc.as_ref(), fe, le, &mut box2d);
        match box2d.get() {
            Some((umin, vmin, umax, vmax)) => (umin, umax, vmin, vmax),
            None => (0.0, 0.0, 0.0, 0.0),
        }
    };

    // `cxx:183-186`
    let e_iso_u = (ev_max - ev_min).abs() < occt_core::precision::CONFUSION;
    let e_iso_v = (eu_max - eu_min).abs() < occt_core::precision::CONFUSION;
    // `cxx:187-193`
    let ttu = 1.0e-2;
    let paru = (ttu * (fu_max - fu_min)).abs();
    let ttv = 1.0e-2;
    let parv = (ttv * (fv_max - fv_min)).abs();

    // `cxx:195-210`
    let mut up = u;
    let mut vp = v;
    if e_iso_v {
        up += paru;
    } else if e_iso_u {
        vp += parv;
    } else {
        up += paru;
        vp += parv;
    }
    let pb = sff.d0(up, vp);

    // `cxx:212-228`
    let mut um = u;
    let mut vm = v;
    if e_iso_v {
        um -= paru;
    } else if e_iso_u {
        vm -= parv;
    } else {
        um -= paru;
        vm -= parv;
    }
    let pa = sff.d0(um, vm);

    // `cxx:230-263`: on the `EisoU` arm and on the general arm OCCT probes
    // `d1v`; only the `EisoV` arm probes `d1u`.
    let mut permute = false;
    let ve_cro_n = ve.crossed(&n);
    let dot = if e_iso_v { ve_cro_n.dot(&d1u) } else { ve_cro_n.dot(&d1v) };
    if dot < 0.0 {
        permute = true;
    }
    let (pa, pb) = if permute { (pb, pa) } else { (pa, pb) };

    match p_classif {
        Some(psc) => {
            // `cxx:265-272`
            if psc.state_p3d_reference(&pb).is_err() {
                return (stateb, statea);
            }
            stateb = psc.state();
            if psc.state_p3d_reference(&pa).is_err() {
                return (stateb, statea);
            }
            statea = psc.state();
        }
        None => {
            // `cxx:274-278`
            stateb = fun_sta_pin_f_3d(&pb, fs);
            statea = fun_sta_pin_f_3d(&pa, fs);
        }
    }
    (stateb, statea)
}
