//! `TopOpeBRep_ProcessGR.cxx` — the geometry/restriction helpers that
//! `ProcessSectionEdges`, `KeepRLine` and the VPoint processing read.
//!
//! Sources (OCCT 8.0.0, `ModelingAlgorithms/TKBool/TopOpeBRep`):
//! - `StBipVPonF` (`TopOpeBRep_ProcessGR.cxx:59-186`)
//! - `StateVPonFace` (`ProcessGR.cxx:192-224`)
//! - `LSameDomainERL` (`ProcessGR.cxx:262-312`)
//! - `VPParamOnER` (`ProcessGR.cxx:368-408`)
//! - `FUN_EqualPonR` (`ProcessGR.cxx:410-423`)
//! - `FUN_EqualponR` (`ProcessGR.cxx:425-432`)
//! - `EqualpPonR` (`ProcessGR.cxx:436-443`)
//!
//! `Lminmax` / `IsVPtransLok` / `TransvpOK` already live in
//! [`super::FacesFiller`].
//!
//! Translation notes:
//! - `myPShapeClassifier` is a pointer in OCCT, so `StBipVPonF` and
//!   `StateVPonFace` are `const` while still writing the classifier's state.
//!   The port passes `&mut ShapeClassifier` explicitly.
//! - The two `TopoDS_Shape F` members are passed as the two faces, because a
//!   method on `FacesFiller` cannot borrow `self.f1`/`self.f2` and
//!   `self.pshape_classifier` at the same time.
//! - `GLOBAL_tolFF` (`ProcessGR.cxx:55`, defined at
//!   `TopOpeBRep_FacesIntersector.cxx:37`, set at `:212`) is
//!   [`FacesIntersector::global_tol_ff`]; `LSameDomainERL` takes it as an
//!   argument.

use occt_core::gp::GpPnt2d;
use occt_core::precision::PCONFUSION;

use crate::abs::{Orientation, ShapeType, State};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, TopoShape};
use crate::topopebrep::ff_transition_tool::FFTransitionTool;
use crate::topopebrep::geom_tool::GeomTool;
use crate::topopebrep::interf::VPointInter;
use crate::topopebrep::line_inter::LineInter;
use crate::topopebrep::tool::adaptor::CurveAdaptor;
use crate::topopebrep::tool::project_curve::{
    fun_tool_bounds, fun_tool_proj_pon_e_tol,
};
use crate::topopebrep::tool::sc::fsc_state_pon_face;
use crate::topopebrep::tool::shape_classifier::ShapeClassifier;
use crate::topopebrep::tool::topology::fun_tool_pin_c_range;
use crate::topopebrep::type_line_curve::TypeLineCurve;

/// `StBipVPonF(vpf, vpl, Lrest, isonedge1)` (`ProcessGR.cxx:59-186`): the
/// state of the middle of the restriction arc `Lrest` on the *other* face.
///
/// `f1` / `f2` are `myF1` / `myF2`; `f2` is classified when the arc is on edge
/// 1 and `f1` otherwise (`cxx:144-151`).
pub fn st_bip_vp_on_f(
    psc: &mut ShapeClassifier,
    f1: &TopoShape,
    f2: &TopoShape,
    vpf: &VPointInter,
    vpl: &VPointInter,
    l: &LineInter,
    is_on_edge1: bool,
) -> Result<State, String> {
    let sind = if is_on_edge1 { 2 } else { 1 };
    let stf = vpf.state(sind)?;
    let stl = vpl.state(sind)?;
    if stf == State::Out || stl == State::Out {
        return Ok(State::Out);
    }
    if stf == State::In || stl == State::In {
        return Ok(State::In);
    }

    let arc = l
        .arc()
        .cloned()
        .ok_or_else(|| "TopOpeBRep_FacesFiller::StBipVPonF : no arc".to_string())?;
    let bac = CurveAdaptor::of_edge(&Edge(arc.clone()));
    let ct = bac.get_type();
    let isperiodic = ct == crate::geom_abs_curve_type::GeomAbsCurveType::Circle
        || ct == crate::geom_abs_curve_type::GeomAbsCurveType::Ellipse;

    let mut vpff = vpf.clone();
    let mut vpll = vpl.clone();

    // cxx:95-141: on a periodic restriction the bounds may be described in the
    // reverse order; swap them when the IN orientations are REVERSED/FORWARD.
    if isperiodic {
        let iarc = if l.arc_is_edge(1) {
            1
        } else if l.arc_is_edge(2) {
            2
        } else {
            0
        };
        if iarc == 0 {
            // cxx:116-122: `throw Standard_Failure("StBipVPonF")` in debug,
            // `return TopAbs_UNKNOWN` otherwise.
            return Ok(State::Unknown);
        }
        let isi = if iarc == 1 { 2 } else { 1 };
        let sif = vpf.shape_index();
        let sil = vpl.shape_index();
        let act = (sif == 3 || sif == isi) && (sil == 3 || sil == isi);
        if act {
            let ori_f = vpf
                .edge(isi)
                .map(|e| e.orientation())
                .unwrap_or(Orientation::Forward);
            let ori_l = vpl
                .edge(isi)
                .map(|e| e.orientation())
                .unwrap_or(Orientation::Forward);
            let tf = FFTransitionTool::process_line_transition_with_index(vpf, isi, ori_f);
            let tl = FFTransitionTool::process_line_transition_with_index(vpl, isi, ori_l);
            let toreverse = tf.orientation(State::In, ShapeType::Face) == Orientation::Reversed;
            let toreverse = toreverse
                && tl.orientation(State::In, ShapeType::Face) == Orientation::Forward;
            if toreverse {
                vpff = vpl.clone();
                vpll = vpf.clone();
            }
        }
    }

    let f = if is_on_edge1 { f2 } else { f1 };
    let mut uf = vp_param_on_er(&vpff, l)?;
    let mut ul = vp_param_on_er(&vpll, l)?;

    let badparametrized = uf > ul;
    let fpar = bac.first_parameter();
    let lpar = bac.last_parameter();
    if badparametrized && isperiodic {
        if uf == lpar {
            uf = fpar;
        }
        if ul == fpar {
            ul = lpar;
        }
    }

    let x = 0.789;
    let parmil = (1.0 - x) * uf + x * ul;
    let pmil = bac
        .d0(parmil)
        .ok_or_else(|| "TopOpeBRep_FacesFiller::StBipVPonF : no curve".to_string())?;

    Ok(fsc_state_pon_face(&pmil, f, psc))
}

/// `StateVPonFace(VP)` (`ProcessGR.cxx:192-224`): the state of `VP` on the
/// face it is *not* on.
pub fn state_vp_on_face(
    psc: &mut ShapeClassifier,
    f1: &TopoShape,
    f2: &TopoShape,
    vp: &VPointInter,
) -> Result<State, String> {
    let ivp = vp.shape_index();
    if ivp == 3 {
        return Ok(State::On);
    }
    let iother = if ivp == 1 { 2 } else { 1 };
    let f = if iother == 1 { f1 } else { f2 };
    let (u, v) = if iother == 1 {
        vp.parameters_on_s1()?
    } else {
        vp.parameters_on_s2()?
    };
    psc.set_reference(f);
    psc.state_p2d_reference(&GpPnt2d::new(u, v))?;
    Ok(psc.state())
}

/// `LSameDomainERL(L, ERL)` (`ProcessGR.cxx:262-312`): `true` when the line
/// shares a geometric domain with at least one edge of `ERL`.
///
/// `global_tol_ff` is `GLOBAL_tolFF` (`cxx:301`). `l` is `&mut` because
/// `TopOpeBRep_GeomTool::MakeCurve` reads the non-const `LineInter::Curve()`.
pub fn lsame_domain_erl(
    l: &mut LineInter,
    erl: &[TopoShape],
    global_tol_ff: f64,
) -> Result<bool, String> {
    if l.type_line_curve() == TypeLineCurve::Walking {
        return Ok(false);
    }

    let (f0, l0) = super::FacesFiller::lminmax(l);
    let d = (f0 - l0).abs();

    if l.inl() && d == 0.0 {
        // null length line, made of VPoints only
        return Ok(false);
    }
    if d <= PCONFUSION {
        return Ok(false);
    }

    let cl = GeomTool::make_curve(f0, l0, l)?
        .ok_or_else(|| "TopOpeBRep_FacesFiller::LSameDomainERL : no curve".to_string())?;
    let t = 0.417789;
    let p = (1.0 - t) * f0 + t * l0;
    let pm = cl.d0(p);

    for e in erl {
        let tole = BRepTool::edge_tolerance(&Edge(e.clone()));
        let maxtol = if tole > global_tol_ff {
            tole
        } else {
            global_tol_ff
        };
        let bac = CurveAdaptor::of_edge(&Edge(e.clone()));
        let f = bac.first_parameter();
        let lp = bac.last_parameter();
        if fun_tool_pin_c_range(&pm, &bac, f, lp, maxtol) {
            return Ok(true);
        }
    }
    Ok(false)
}

/// `VPParamOnER(vp, Lrest)` (`ProcessGR.cxx:368-408`): the parameter of `vp`
/// on the restriction edge of `Lrest`.
///
/// The three "edge boundary / on domain" arms return directly; otherwise the
/// point is projected with `tole = BRep_Tool::Tolerance(E) * 1.e2`
/// (`cxx:400-402`). A failed projection raises `Standard_ProgramError` in OCCT
/// (`cxx:405-407`) and is the `Err` here.
pub fn vp_param_on_er(vp: &VPointInter, l: &LineInter) -> Result<f64, String> {
    let e = l
        .arc()
        .ok_or_else(|| "TopOpeBRep_FacesFiller::VPParamOnER : no arc".to_string())?;
    let isedge1 = l.arc_is_edge(1);
    let isedge2 = l.arc_is_edge(2);

    if isedge1 && vp.is_vertex_on_s1()? {
        let v1 = vp.vertex_on_s1()?;
        return Ok(BRepTool::parameter_of_vertex(e, v1));
    }
    if isedge2 && vp.is_vertex_on_s2()? {
        let v2 = vp.vertex_on_s2()?;
        return Ok(BRepTool::parameter_of_vertex(e, v2));
    }
    if isedge1 && vp.is_on_dom_s1()? {
        return vp.parameter_on_arc1();
    }
    if isedge2 && vp.is_on_dom_s2()? {
        return vp.parameter_on_arc2();
    }

    let mut tolee = BRepTool::edge_tolerance(&Edge(e.clone()));
    tolee *= 1.0e2; // xpu290998 : PRO15369
    let p = vp.value()?;
    match fun_tool_proj_pon_e_tol(p, tolee, e) {
        Some((param, _dist)) => Ok(param),
        None => Err("TopOpeBRep_FacesFiller::VPParamOnER".to_string()),
    }
}

/// `FUN_EqualPonR(Lrest, VP1, VP2)` (`ProcessGR.cxx:410-423`): equal 3D point
/// at the larger of the two VPoint tolerances.
pub fn equal_p_on_r(vp1: &VPointInter, vp2: &VPointInter) -> Result<bool, String> {
    let p1 = vp1.value()?;
    let p2 = vp2.value()?;
    let ptol1 = vp1.tolerance()?;
    let ptol2 = vp2.tolerance()?;
    let ptol = if ptol1 > ptol2 { ptol1 } else { ptol2 };
    Ok(p1.is_equal_with_tolerance(p2, ptol))
}

/// `FUN_EqualponR(Lrest, VP1, VP2)` (`ProcessGR.cxx:425-432`): equal parameter
/// on the restriction edge, within `Precision::PConfusion()`.
pub fn equal_pon_r(
    vp1: &VPointInter,
    vp2: &VPointInter,
    l: &LineInter,
) -> Result<bool, String> {
    let p1 = vp_param_on_er(vp1, l)?;
    let p2 = vp_param_on_er(vp2, l)?;
    Ok((p1 - p2).abs() < PCONFUSION)
}

/// `EqualpPonR(Lrest, VP1, VP2)` (`ProcessGR.cxx:436-443`): both the 3D point
/// and the parameter on the restriction edge are equal.
pub fn equal_p_pon_r(
    vp1: &VPointInter,
    vp2: &VPointInter,
    l: &LineInter,
) -> Result<bool, String> {
    let pequal = equal_p_on_r(vp1, vp2)?;
    let pequal_par = equal_pon_r(vp1, vp2, l)?;
    Ok(pequal && pequal_par)
}
