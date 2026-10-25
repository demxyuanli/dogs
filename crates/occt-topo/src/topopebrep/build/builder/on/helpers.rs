//! `TopOpeBRepBuild_BuilderON` — the free helpers of `TopOpeBRepBuild_BuilderON.cxx`.
//!
//! Sources (`ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_BuilderON.cxx`,
//! OCCT 8.0.0):
//! - `FUN_keepEON` (`:303-410`)
//! - `FUN_build_TB` (`:412-436`)
//! - `FUN_Kpart0` / `FUN_Kpart` (`:438-476`)
//! - `FUN_findeSD` (`:479-516`)
//! - `ComputeFaceCrvtInSec` (`:519-552`)
//!
//! The `#ifdef OCCT_DEBUG` arms of `FUN_keepEON` (`:331-348`, `:395-407`) only
//! print the three `keep` values; the port keeps the `keep3` computation the
//! function returns.

use occt_core::gp::{GpDir, GpPnt2d};
use occt_geom::geom_lprop::sl_props::SlProps;

use crate::abs::{Orientation, State};
use crate::brep_tool::BRepTool;
use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::shape::{Edge, Face, TopoShape};
use crate::topopebrep::ds::data_structure::DataStructure;
use crate::topopebrep::ds::export::{fds_has_same_domain2d, fds_has_same_domain3d};
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::tool::{
    fun_tool_find_p_in_e, fun_tool_orient_ein_f, fun_tool_sta_pin_e,
};

use super::super::Builder;

/// `FUN_build_TB(PB, rank)` (`cxx:412-436`): the state of `FOR` implied by the
/// operation, in the `Opefus` / `Opecom` / `Opec12` / `Opec21` order of the
/// C++ `if`s (each later hit overwrites `sta`).
pub(super) fn fun_build_tb(b: &Builder, rank: i32) -> State {
    let mut sta = State::Unknown;
    if b.opefus() {
        sta = State::Out;
    }
    if b.opecom() {
        sta = State::In;
    }
    if b.opec12() {
        sta = if rank == 1 { State::Out } else { State::In };
    }
    if b.opec21() {
        sta = if rank == 2 { State::Out } else { State::In };
    }
    sta
}

/// `FUN_Kpart0(ST1, ST2)` (`cxx:438-462`).
pub(super) fn fun_kpart0(st1: GeomAbsSurfaceType, st2: GeomAbsSurfaceType) -> bool {
    let plane1 = st1 == GeomAbsSurfaceType::Plane;
    let cyli1 = st1 == GeomAbsSurfaceType::Cylinder;
    let cone1 = st1 == GeomAbsSurfaceType::Cone;

    let cyli2 = st2 == GeomAbsSurfaceType::Cylinder;
    let cone2 = st2 == GeomAbsSurfaceType::Cone;
    let sphe2 = st2 == GeomAbsSurfaceType::Sphere;
    let quad2 = cyli2 || cone2 || sphe2 || (st2 == GeomAbsSurfaceType::Torus);

    if plane1 && quad2 {
        return true;
    }
    if cyli1 && sphe2 {
        return true;
    }
    if cone1 && sphe2 {
        return true;
    }
    false
}

/// `FUN_Kpart(ST1, ST2)` (`cxx:464-476`): `1` when `(ST1, ST2)` matches,
/// `2` when the swapped pair matches, `0` otherwise.
pub(super) fn fun_kpart(st1: GeomAbsSurfaceType, st2: GeomAbsSurfaceType) -> i32 {
    if fun_kpart0(st1, st2) {
        return 1;
    }
    if fun_kpart0(st2, st1) {
        return 2;
    }
    0
}

/// `FUN_findeSD(BDS, EspON, EG, FOR, oeSD, D)` (`cxx:479-516`): the
/// same-domain partner of `EG` on `FOR` (3-D for `d == 3`, 2-D for `d == 2`)
/// that contains `EspON`.
///
/// The probe point is `FUN_tool_findPinE(EspON)` (`cxx:489`), a partner is
/// kept only when it is a bound of `FOR` (`cxx:503-506`) and the probe lies
/// `IN` it (`cxx:508`). `oe_sd` receives that bound's orientation in `FOR`,
/// and the answer is the DS index (`0` when nothing matched).
pub(super) fn fun_find_e_sd(
    bds: &DataStructure,
    esp_on: &TopoShape,
    eg: &TopoShape,
    for_: &TopoShape,
    oe_sd: &mut Orientation,
    d: i32,
) -> i32 {
    let Some((pt_on, _par_on)) = fun_tool_find_p_in_e(esp_on) else {
        return 0;
    };
    let lesd = if d == 3 {
        fds_has_same_domain3d(bds, eg, true).1
    } else if d == 2 {
        fds_has_same_domain2d(bds, eg, true).1
    } else {
        Vec::new()
    };

    for e_sd in lesd {
        // `cxx:503`: `FUN_tool_orientEinFFORWARD(eSD, FOR, oesd)`.
        let Some(oesd) = crate::topopebrep::tool::fun_tool_orient_ein_f_forward(&e_sd, for_)
        else {
            continue;
        };
        if fun_tool_sta_pin_e(&pt_on, &e_sd) == State::In {
            *oe_sd = oesd;
            return bds.shape_index(&e_sd, true);
        }
    }
    0
}

/// `FUN_keepEON(B, sEG, sFOR, sFS, EGBoundFOR, TFE, TB1, TB2)` (`cxx:303-410`).
///
/// The `Builder` argument and `TB2` are read only inside `#ifdef OCCT_DEBUG`,
/// so they are not parameters here (`sFS` is unused by the C++ body as well;
/// it is kept for the call sites).
pub(super) fn fun_keep_eon(
    s_eg: &TopoShape,
    s_for: &TopoShape,
    tfe: &Transition,
    tb1: State,
) -> bool {
    // `FF = FOR` with a FORWARD orientation (`cxx:328`); `FUN_tool_orientEinF`
    // reads the stored occurrence orientation, which `Oriented(FORWARD)`
    // leaves untouched.
    let ff = s_for;
    let is_closed_ff = BRepTool::is_closed_edge_face(&Edge(s_eg.clone()), &Face(ff.clone()));
    if is_closed_ff {
        return true;
    }

    let Some(o_eg_ff) = fun_tool_orient_ein_f(s_eg, ff) else {
        // `FUN_tool_orientEinF` leaves `oriEinF = FORWARD` when the edge is not
        // a bound of `FF` (`TOPOLOGY.cxx:404-415`) and the caller ignores `ok`
        // (`BuilderON.cxx:365`).
        return tfe.after_state() == tb1;
    };
    let tfe_b = tfe.before_state();
    let tfe_a = tfe.after_state();
    match o_eg_ff {
        Orientation::Forward => tfe_a == tb1,
        Orientation::Reversed => tfe_b == tb1,
        Orientation::External | Orientation::Internal => tfe_a == tb1 || tfe_b == tb1,
    }
}

/// `ComputeFaceCrvtInSec(aFace, aP2d, aSecDir, aCrvt)` (`cxx:519-552`): the
/// normal curvature of `aFace` at `aP2d` along `aSecDir`.
///
/// `false` is `cxx:526-529` (continuity below `C2`) and `cxx:531-534`
/// (curvature undefined).
pub(super) fn compute_face_crvt_in_sec(
    face: &Face,
    p2d: &GpPnt2d,
    sec_dir: &GpDir,
) -> Option<f64> {
    let surf = BRepTool::face_surface_world(face)?;
    // `cxx:524-528`: `BRepLProp_SurfaceTool::Continuity` is
    // `min(UContinuity(), VContinuity())`, both `Geom_Surface::Continuity()`.
    if (surf.continuity() as i32) < 2 {
        return None;
    }
    let mut prop = SlProps::new(
        surf.as_ref(),
        p2d.x(),
        p2d.y(),
        2,
        occt_core::precision::CONFUSION,
    )
    .ok()?;
    // `cxx:530-534`
    if !prop.is_curvature_defined(surf.as_ref()) {
        return None;
    }
    if prop.is_umbilic(surf.as_ref()).ok()? {
        return prop.max_curvature(surf.as_ref()).ok();
    }
    let max_crv = prop.max_curvature(surf.as_ref()).ok()?;
    let min_crv = prop.min_curvature(surf.as_ref()).ok()?;
    let (max_dir, min_dir) = prop.curvature_directions(surf.as_ref()).ok()?;
    let cos_max = sec_dir.dot(&max_dir);
    let cos_min = sec_dir.dot(&min_dir);
    Some(max_crv * cos_max + min_crv * cos_min)
}

/// `M_FORWARD` / `M_REVERSED` / `M_INTERNAL` / `M_EXTERNAL` (`cxx:134-137`).
pub(super) const M_FORWARD: Orientation = Orientation::Forward;
/// `M_REVERSED` (`cxx:135`).
pub(super) const M_REVERSED: Orientation = Orientation::Reversed;
/// `M_INTERNAL` (`cxx:136`).
pub(super) const M_INTERNAL: Orientation = Orientation::Internal;
/// `M_EXTERNAL` (`cxx:137`).
pub(super) const M_EXTERNAL: Orientation = Orientation::External;
