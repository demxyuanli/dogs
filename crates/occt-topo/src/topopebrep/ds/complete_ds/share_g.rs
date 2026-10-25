//! `FUN_ds_shareG`, `FUN_ds_mkTonFsdm` and `FUN_ds_oriEinF`.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_EXPORT.cxx:1423-1791`.

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};
use crate::topopebrep::ds::data_structure::DataStructure;
use crate::topopebrep::ds::export::fun_ds_sdm;
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::tool::curvature::get_xx;
use crate::topopebrep::tool::tangent::tg_geom_e;
use crate::topopebrep::tool::project_curve::{
    fun_tool_bounds, fun_tool_proj_pon_e, fun_tool_value,
};
use crate::topopebrep::tool::topology::{fun_tool_ng_geom_f_par, fun_tool_orient_ein_f_forward};

/// `#define UNKNOWN (0)` (`cxx:1669`).
pub const ORI_E_IN_F_UNKNOWN: i32 = 0;
/// `#define ONSAMESHA (1)` (`cxx:1670`).
pub const ORI_E_IN_F_ON_SAME_SHA: i32 = 1;
/// `#define CLOSESAME (11)` (`cxx:1671`).
pub const ORI_E_IN_F_CLOSE_SAME: i32 = 11;
/// `#define ONOPPOSHA (2)` (`cxx:1672`).
pub const ORI_E_IN_F_ON_OPPO_SHA: i32 = 2;
/// `#define CLOSEOPPO (22)` (`cxx:1673`).
pub const ORI_E_IN_F_CLOSE_OPPO: i32 = 22;
/// `#define FORREVOPPO (222)` (`cxx:1674`).
pub const ORI_E_IN_F_FOR_REV_OPPO: i32 = 222;

/// `FUN_ds_shareG(HDS, iF1, iF2, iE2, Esp, shareG)` (`cxx:1423-1522`).
///
/// `F1` same-domain `F2`, sharing the split `Esp` of `E2`. `None` is the OCCT
/// `false` return; `Some(share)` is `true` with the `shareG` output.
pub fn fun_ds_share_g(
    hds: &HDataStructure,
    if1: i32,
    if2: i32,
    ie2: i32,
    esp: &TopoShape,
) -> Option<bool> {
    let bds = hds.ds();
    // `cxx:1433-1437`
    let e2_shape = bds.shape(ie2, true)?;
    if !hds.has_same_domain(&e2_shape, true) {
        return None;
    }

    // `cxx:1439-1441`
    let f1 = bds.shape(if1, true)?;
    let f2 = bds.shape(if2, true)?;
    let e2 = bds.shape(ie2, true)?;

    // `cxx:1443-1466`
    let tol = occt_core::precision::CONFUSION * 1.0e3;
    let (f, l) = fun_tool_bounds(esp);
    let x = 0.45678;
    let par = (1.0 - x) * f + x * l;
    let p = fun_tool_value(par, esp)?;
    let (par2, d2) = fun_tool_proj_pon_e(&p, &e2)?;
    if d2 > tol {
        return None;
    }

    // `cxx:1468-1476`
    let ng_geom_f2 = fun_tool_ng_geom_f_par(par2, &e2, &f2)?;
    let nxx2 = get_xx(&e2, &f2, par2, &ng_geom_f2)?;

    // `cxx:1478-1519`
    let mut mf1: Vec<TopoShape> = Vec::new();
    let mut ex = crate::topexp::Explorer::new(&f1, ShapeType::Edge);
    while ex.more() {
        mf1.push(ex.current().clone());
        ex.next();
    }

    for e1 in bds.shape_same_domain(&e2_shape) {
        if !mf1.iter().any(|e| crate::topo_tools::is_same(e, e1)) {
            continue;
        }
        let Some((par1, d1)) = fun_tool_proj_pon_e(&p, e1) else {
            continue;
        };
        if d1 > tol {
            continue;
        }
        // E1 on F1, E2 on F2, E1 and E2 share the split Esp (nyi: check it)
        let Some(ng_geom_f1) = fun_tool_ng_geom_f_par(par1, e1, &f1) else {
            return None;
        };
        let Some(nxx1) = get_xx(e1, &f1, par1, &ng_geom_f1) else {
            return None;
        };
        let prod = nxx1.dot(&nxx2);
        return Some(prod > 0.0);
    }
    None
}

/// `FUN_ds_mkTonFsdm(HDS, iF1, iF2, iE2, iEG, paronEG, Esp, pardef, T)`
/// (`cxx:1523-1668`): the 2d transition of `F1` when crossing `EG` relative to
/// the 2d matter of `F2` limited by `E2`.
///
/// When `pardef` is false OCCT projects the still default-constructed `gp_Pnt`
/// (`cxx:1565-1577`, `P` is only assigned at `cxx:1585`), i.e. the origin; the
/// port keeps that call verbatim.
pub fn fun_ds_mk_ton_fsdm(
    hds: &HDataStructure,
    if1: i32,
    if2: i32,
    ie2: i32,
    ieg: i32,
    paron_eg: f64,
    esp: &TopoShape,
    pardef: bool,
) -> Option<crate::topopebrep::ds::transition::Transition> {
    let bds = hds.ds();
    // `cxx:1547-1552`
    let e2_shape = bds.shape(ie2, true)?;
    if !hds.has_same_domain(&e2_shape, true) {
        return None;
    }

    // `cxx:1554-1558`
    let f1 = bds.shape(if1, true)?;
    let f2 = bds.shape(if2, true)?;
    let e2 = bds.shape(ie2, true)?;
    let eg = bds.shape(ieg, true)?;

    let eg_is_e2 = ieg == ie2;
    let tol = occt_core::precision::CONFUSION * 1.0e3;

    // `cxx:1562-1590`
    let mut p = occt_core::gp::GpPnt::new(0.0, 0.0, 0.0);
    let par_eg;
    if pardef {
        par_eg = paron_eg;
    } else {
        let (_f, _l) = fun_tool_bounds(esp);
        let (par, d_eg) = fun_tool_proj_pon_e(&p, &eg)?;
        if d_eg > tol {
            return None;
        }
        par_eg = par;
    }
    p = fun_tool_value(par_eg, &eg)?;

    // `cxx:1591-1605`
    let tgt_eg = tg_geom_e(par_eg, &eg)?;
    let ng_f1 = fun_tool_ng_geom_f_par(par_eg, &eg, &f1)?;
    let beafter =
        occt_core::gp::GpVec::new(ng_f1.x(), ng_f1.y(), ng_f1.z()).crossed(&tgt_eg);

    // `cxx:1607-1636`
    let par2 = if eg_is_e2 {
        par_eg
    } else {
        let (par, d2) = fun_tool_proj_pon_e(&p, &e2)?;
        if d2 > tol {
            return None;
        }
        par
    };
    let ng_f2 = fun_tool_ng_geom_f_par(par2, &e2, &f2)?;
    let nxx2 = get_xx(&e2, &f2, par2, &ng_f2)?;

    // `cxx:1638-1648`
    let mut sdm_eg_e2 = eg_is_e2;
    if !sdm_eg_e2 {
        sdm_eg_e2 = fun_ds_sdm(bds, &eg, &e2);
    }
    if !sdm_eg_e2 {
        return None;
    }

    // `cxx:1650-1667`
    let n2 = occt_core::gp::GpVec::new(nxx2.x(), nxx2.y(), nxx2.z());
    let prod = beafter.dot(&n2);
    let tola = occt_core::precision::ANGULAR * 1.0e3;
    if (1.0 - prod.abs()).abs() >= tola {
        return None;
    }
    let mut t = crate::topopebrep::ds::transition::Transition::default();
    if prod > 0.0 {
        t.set_state_before(crate::abs::State::Out);
        t.set_state_after(crate::abs::State::In);
    } else {
        t.set_state_before(crate::abs::State::In);
        t.set_state_after(crate::abs::State::Out);
    }
    Some(t)
}

/// `FUN_ds_oriEinF(BDS, E, F, O)` (`cxx:1676-1791`): the orientation of `E`
/// in `F`, following the same-domain face of the other ancestor rank.
///
/// `O` is the `TopAbs_Orientation` output; the `int` result is one of the
/// `ORI_E_IN_F_*` constants.
pub fn fun_ds_ori_e_in_f(
    bds: &DataStructure,
    e: &TopoShape,
    f: &TopoShape,
) -> (i32, Orientation) {
    // `cxx:1687-1692`
    let mut o = Orientation::External;
    let rk_f = bds.ancestor_rank_of(f);
    let rk_e = bds.ancestor_rank_of(e);

    let ee = e.clone();
    let ff = f.clone();
    let i_f = bds.shape_index(f, true);
    let o_f = bds.shape(i_f, true).map(|s| s.orientation()).unwrap_or(Orientation::Forward);

    if rk_f == rk_e {
        // `cxx:1694-1713`
        let Some(o_ee) = fun_tool_orient_ein_f_forward(&ee, &ff) else {
            return (ORI_E_IN_F_UNKNOWN, o);
        };
        o = o_ee;
        let isc_e = BRepTool::is_closed_edge_face(&Edge(ee.clone()), &Face(ff.clone()));
        if isc_e {
            return (ORI_E_IN_F_CLOSE_SAME, o);
        }
        return (ORI_E_IN_F_ON_SAME_SHA, o);
    }

    // `cxx:1714-1789`
    let sdm_fs = bds.shape_same_domain(&ff).to_vec();
    if sdm_fs.is_empty() {
        return (ORI_E_IN_F_UNKNOWN, o);
    }

    let mut has_for = false;
    let mut has_rev = false;
    let c = bds.same_domain_ori_of(&ff);

    for fsdm in &sdm_fs {
        let i_fsdm = bds.shape_index(fsdm, true);
        let rksdm = bds.ancestor_rank_of(fsdm);
        if rksdm == rk_f {
            continue;
        }

        let Some(o_fsdm) = fun_tool_orient_ein_f_forward(&ee, fsdm) else {
            continue;
        };
        o = o_fsdm;

        let isc_e = BRepTool::is_closed_edge_face(&Edge(ee.clone()), &Face(fsdm.clone()));
        if isc_e {
            return (ORI_E_IN_F_CLOSE_OPPO, o);
        }

        let c_sdm = bds.same_domain_ori_of(fsdm);
        let to_reverse1 = c_sdm != c && o != Orientation::Internal && o != Orientation::External;
        if to_reverse1 {
            o = o.complement();
        }

        let o_fsdm_shape = bds
            .shape(i_fsdm, true)
            .map(|s| s.orientation())
            .unwrap_or(Orientation::Forward);
        let to_reverse2 = o_f != o_fsdm_shape
            && o_fsdm_shape != Orientation::Internal
            && o_fsdm_shape != Orientation::External;
        if to_reverse2 {
            o = o.complement();
        }
        if !has_for {
            has_for = o == Orientation::Forward;
        }
        if !has_rev {
            has_rev = o == Orientation::Reversed;
        }
    }

    // `cxx:1774-1785`
    if has_for && has_rev {
        return (ORI_E_IN_F_FOR_REV_OPPO, o);
    }
    if has_for || has_rev {
        o = if has_for {
            Orientation::Forward
        } else {
            Orientation::Reversed
        };
        return (ORI_E_IN_F_ON_OPPO_SHA, o);
    }
    (ORI_E_IN_F_UNKNOWN, o)
}
