//! `TopOpeBRepDS_EXPORT` — the free functions of the classic boolean engine
//! (`FDS_*` / `FUN_ds_*`).
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_EXPORT.hxx:24-91`,
//! `TopOpeBRepDS_EXPORT.cxx:1-3952`.
//!
//! Ported here: the pure list / interference accessors and the same-domain
//! queries that need no geometry evaluation (`cxx:61-360`, `396-462`,
//! `464-495`, `647-712`, `1874-1995`). The `FUN_ds_completefor*` /
//! `FUN_ds_PURGEforE9` / `FUN_ds_PointToVertex` / `FUN_ds_redusamsha`
//! family lives in [`super::complete_ds`], keyed by its own `.cxx` lines.
//!
//! Still UNPORTED in this file's source range: `FDS_HasSameDomain3d`'s
//! geometry branch is in [`fds_has_same_domain3d`]; `FDS_SIisGIofIofSBAofTofI`
//! is [`fds_si_is_gi_of_i_of_sba_of_t_of_i`]. `FUN_ds_hasFEI` (`:3279`),
//! `FUN_ds_ONesd` (`:3308`), `FUN_ds_unkeepEVIonGb1` (`:3364`),
//! `FDS_LOIinfsup`'s 2d half (`:3427`), `FDS_parbefaft` (`:3571`),
//! `FDS_stateEwithF2d` (`:3598`) and `FUN_ds_FEIGb1TO0` (`:3647`) are the
//! remaining entry points of the file; `FUN_ds_unkeepEVIonGb1` is only
//! called from commented-out code in `TopOpeBRep_DSFiller.cxx:769`.
//!
//! OCCT's `bool(a, b)` output pairs become `Option`/tuple returns, and the
//! `bool FDS_data(it, I, ...)` overload becomes [`fds_data_is_unknown`].

use crate::shape::TopoShape;
use occt_core::gp::GpDir;

use super::data_structure::DataStructure;
use super::h_data_structure::HDataStructure;
use super::interference::Interference;
use super::interference_iterator::InterferenceList;
use super::kind::{Config, Kind};
use super::tki::TKI;
use super::transition::Transition;

/// `FDS_SetT(T, T0)` (`cxx:61-68`).
pub fn fds_set_t(t: &mut Transition, t0: &Transition) {
    let stb = t0.before_state();
    let sta = t0.after_state();
    t.set_state_before(stb);
    t.set_state_after(sta);
}

/// `FDS_hasUNK(T)` (`cxx:70-77`): either side is `UNKNOWN`.
pub fn fds_has_unk(t: &Transition) -> bool {
    t.before_state() == crate::abs::State::Unknown || t.after_state() == crate::abs::State::Unknown
}

/// `FDS_copy(LI, LII)` (`cxx:79-88`): append every element of `LI` to `LII`.
pub fn fds_copy(li: &InterferenceList, lii: &mut InterferenceList) {
    lii.extend(li.iter().cloned());
}

/// `FDS_copy(LI, LII)` for `List<TopoDS_Shape>` (`cxx:91-99`).
pub fn fds_copy_shapes(li: &[TopoShape], lii: &mut Vec<TopoShape>) {
    lii.extend(li.iter().cloned());
}

/// `FDS_assign(LI, LII)` (`cxx:101-108`): clear then copy.
pub fn fds_assign(li: &InterferenceList, lii: &mut InterferenceList) {
    lii.clear();
    fds_copy(li, lii);
}

/// `FDS_assign(LI, LII)` for `List<TopoDS_Shape>` (`cxx:110-117`).
pub fn fds_assign_shapes(li: &[TopoShape], lii: &mut Vec<TopoShape>) {
    lii.clear();
    fds_copy_shapes(li, lii);
}

/// `FUN_ds_samRk(BDS, Rk, LI, LIsrk)` (`cxx:119-142`): move the elements of
/// `LI` whose `AncestorRank` is `Rk` into `LIsrk` (which is cleared first).
pub fn fun_ds_sam_rk(bds: &DataStructure, rk: i32, li: &mut Vec<TopoShape>, lisrk: &mut Vec<TopoShape>) {
    lisrk.clear();
    let mut i = 0;
    while i < li.len() {
        let crk = bds.ancestor_rank_of(&li[i]);
        if crk == rk {
            let moved = li.remove(i);
            lisrk.push(moved);
        } else {
            i += 1;
        }
    }
}

/// `FDS_Tdata(I, SB, IB, SA, IA)` (`cxx:144-162`): the four transition sides.
pub fn fds_tdata(t: &Transition) -> (crate::abs::ShapeType, i32, crate::abs::ShapeType, i32) {
    (
        t.shape_before(),
        t.index_before(),
        t.shape_after(),
        t.index_after(),
    )
}

/// `FDS_data(I, GT1, G1, ST1, S1)` (`cxx:164-181`).
pub fn fds_data(i: &Interference) -> (Kind, i32, Kind, i32) {
    (i.geometry_type(), i.geometry(), i.support_type(), i.support())
}

/// `FDS_data(it1, I1, GT1, G1, ST1, S1)` (`cxx:184-203`): fills the four
/// outputs and reports whether the transition is `(UNKNOWN, UNKNOWN)`.
pub fn fds_data_is_unknown(i: &Interference) -> bool {
    i.transition().is_unknown()
}

/// `FDS_Idata(I, SB, IB, SA, IA, GT1, G1, ST1, S1)` (`cxx:205-223`): the
/// transition and data quadruples in one call.
pub fn fds_idata(
    i: &Interference,
) -> (
    crate::abs::ShapeType,
    i32,
    crate::abs::ShapeType,
    i32,
    Kind,
    i32,
    Kind,
    i32,
) {
    let (sb, ib, sa, ia) = fds_tdata(&i.transition());
    let (gt, g, st, s) = fds_data(i);
    (sb, ib, sa, ia, gt, g, st, s)
}

/// `FUN_ds_hasFEI(pDS2d, F, GI, ITRA)` (`EXPORT.cxx:3279-3306`): `pDS2d`
/// already carries an interference of `F` whose geometry is `GI` and whose
/// transition index is `ITRA`.
pub fn fun_ds_has_fei(p_ds2d: &DataStructure, f: &TopoShape, gi: i32, itra: i32) -> bool {
    if !p_ds2d.has_shape(f, true) {
        return false;
    }
    for i in p_ds2d.shape_interferences(f, true) {
        let Ok(guard) = i.read() else {
            continue;
        };
        let (_gt, g, _st, _s) = fds_data(&guard);
        if g == gi && guard.transition().index().unwrap_or(-1) == itra {
            return true;
        }
    }
    false
}

/// `FUN_ds_getVsdm(BDS, iV, iVsdm)` (`cxx:226-252`): the index of the first
/// shape same-domain with `BDS.Shape(iV)` that is not the shape itself.
pub fn fun_ds_get_vsdm(bds: &DataStructure, iv: i32) -> Option<i32> {
    let imax = bds.nb_shapes();
    if iv < 1 || iv > imax {
        return None;
    }
    let v = bds.shape(iv, true)?;
    for vv in bds.shape_same_domain(&v) {
        if crate::topo_tools::is_same(&v, vv) {
            continue;
        }
        return Some(bds.shape_index(vv, true));
    }
    None
}

/// `FUN_ds_sdm(BDS, s1, s2)` (`cxx:255-274`): `s2` is in `s1`'s same-domain
/// list. `HasShape` defaults to `FindKeep = true`.
pub fn fun_ds_sdm(bds: &DataStructure, s1: &TopoShape, s2: &TopoShape) -> bool {
    if !bds.has_shape(s1, true) || !bds.has_shape(s2, true) {
        return false;
    }
    bds.shape_same_domain(s1)
        .iter()
        .any(|s| crate::topo_tools::is_same(s, s2))
}

/// `FUN_ds_getoov(v, BDS, oov)` (`cxx:647-668`): the "other" vertex of a
/// same-domain vertex pair (the DS binds at most two).
pub fn fun_ds_getoov(bds: &DataStructure, v: &TopoShape) -> Option<TopoShape> {
    for vcur in bds.shape_same_domain(v) {
        if crate::topo_tools::is_same(vcur, v) {
            continue;
        }
        return Some(vcur.clone());
    }
    None
}

/// `FUN_ds_getoov(v, HDS, oov)` (`cxx:671-694`): the same on the HDS, which
/// first checks `HasSameDomain(v)`.
pub fn fun_ds_getoov_in_hds(hds: &HDataStructure, v: &TopoShape) -> Option<TopoShape> {
    let hsd = hds.has_same_domain(v, true);
    if !hsd {
        return None;
    }
    fun_ds_getoov(hds.ds(), v)
}

/// `FDS_Parameter(I, par)` (`cxx:464-475`): true only for an
/// `EdgeVertexInterference` or a `CurvePointInterference`.
pub fn fds_parameter(i: &Interference) -> Option<f64> {
    match i.data() {
        super::interference::InterferenceData::CurvePoint { .. }
        | super::interference::InterferenceData::EdgeVertex(_) => i.parameter(),
        _ => None,
    }
}

/// `FDS_Parameter(I)` (`cxx:478-495`): the `Standard_Failure` for any other
/// subclass is `None`.
pub fn fds_parameter_or_fail(i: &Interference) -> Result<f64, String> {
    fds_parameter(i).ok_or_else(|| "FDS_Parameter".to_string())
}

/// `FUN_selectTRAINTinterference(li, liINTERNAL)` (`cxx:696-712`): collect the
/// interferences whose `Orientation(TopAbs_IN)` is `INTERNAL`.
pub fn fun_select_train_tinterference(
    li: &InterferenceList,
    li_internal: &mut InterferenceList,
) -> bool {
    li_internal.clear();
    for i in li {
        let ori = i
            .read()
            .expect("interference lock")
            .transition()
            .orientation(crate::abs::State::In, crate::abs::ShapeType::Face);
        if ori == crate::abs::Orientation::Internal {
            li_internal.push(i.clone());
        }
    }
    !li_internal.is_empty()
}

/// `FDS_aresamdom(BDS, ES, F1, F2)` (`cxx:282-321`): true when the faces
/// `F1`, `F2`, issued from the same origin shape and connected by the edge
/// `ES`, are tangent.
///
/// 1/ the SameDomain coding of the DS; 2/ when the DS carries none — the
/// SameDomain information is not looked for between the faces of one origin
/// shape — the geometric test on the two plane normals at the UV points of
/// `ES` (`cxx:295-318`).
pub fn fds_aresamdom(
    bds: &DataStructure,
    es: &TopoShape,
    f1: &TopoShape,
    f2: &TopoShape,
) -> bool {
    for ss in bds.shape_same_domain(f1) {
        if crate::topo_tools::is_same(ss, f2) {
            return true;
        }
    }

    let face1 = crate::shape::Face(f1.clone());
    let face2 = crate::shape::Face(f2.clone());
    let (Some(su1), Some(su2)) = (
        crate::brep_tool::BRepTool::face_surface(&face1),
        crate::brep_tool::BRepTool::face_surface(&face2),
    ) else {
        return false;
    };
    let su1 = crate::topopebrep::tool::basis_surface(&su1);
    let su2 = crate::topopebrep::tool::basis_surface(&su2);
    let st1 = crate::geom_bnd_lib_surface3d::surface_kind(&*su1);
    let st2 = crate::geom_bnd_lib_surface3d::surface_kind(&*su2);
    let plpl = st1 == crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType::Plane
        && st2 == crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType::Plane;
    if !plpl {
        return false;
    }

    let edge = crate::shape::Edge(es.clone());
    let (p2d1, p2d2) = crate::brep_tool::BRepTool::uv_points(&edge, &face1);
    let d1 = crate::topopebrep::tool::fun_tool_ng_s(&p2d1, &*su1);
    let d2 = crate::topopebrep::tool::fun_tool_ng_s(&p2d2, &*su2);
    let tola = occt_core::precision::ANGULAR;
    let dot = d1.dot(&d2);
    (1.0 - dot.abs()).abs() < tola
}

/// `FDS_aresamdom(BDS, SI, isb1, isb2)` (`cxx:330-357`): the shape-index
/// overload, which filters the three shapes by type first.
pub fn fds_aresamdom_from_indices(bds: &DataStructure, si: i32, isb1: i32, isb2: i32) -> bool {
    let Some(es) = bds.shape(si, true) else {
        return false;
    };
    if es.shape_type() != crate::abs::ShapeType::Edge {
        return false;
    }
    let Some(f1) = bds.shape(isb1, true) else {
        return false;
    };
    if f1.shape_type() != crate::abs::ShapeType::Face {
        return false;
    }
    let Some(f2) = bds.shape(isb2, true) else {
        return false;
    };
    if f2.shape_type() != crate::abs::ShapeType::Face {
        return false;
    }
    fds_aresamdom(bds, &es, &f1, &f2)
}

/// `FDS_EdgeIsConnexToSameDomainFaces(E, HDS)` (`cxx:357-391`): at least two
/// faces connected through `E`, one pair of which `FDS_aresamdom` reports as
/// same-domain.
///
/// The inner loop starts at the outer cursor (`cxx:376`), so the first pair
/// tested is `(f1, f1)`; `samdom` is overwritten on every inner iteration and
/// the loops stop on the first `true` (`cxx:379-388`). `HDS` is unused by
/// OCCT (`cxx:359`) beyond reaching the DS.
pub fn fds_edge_is_connex_to_same_domain_faces(
    connex: &super::connex::Connexity,
    e: &TopoShape,
    hds: &HDataStructure,
) -> bool {
    let bds = hds.ds();
    let lf = connex.edge_connexity_same_shape(e, hds);
    if lf.len() < 2 {
        return false;
    }
    for i1 in 0..lf.len() {
        let f1 = &lf[i1];
        for j in i1..lf.len() {
            let f2 = &lf[j];
            if fds_aresamdom(bds, e, f1, f2) {
                return true;
            }
        }
    }
    false
}

/// `FDS_HasSameDomain3d(BDS, E, PLSD)` (`cxx:498-537`): `E` has a
/// same-domain partner whose `SameDomainOri` is `UNSHGEOMETRY` (3-D geometry
/// not shared).
///
/// The C++ `PLSD` out-parameter collects every such partner; `None` is a null
/// `PLSD`, in which case OCCT `break`s on the first hit (`cxx:529-532`). The
/// port always collects, and the `break`/`continue` split is reproduced by the
/// `collect_all` flag so a caller emitting only the boolean reads the same
/// early exit.
pub fn fds_has_same_domain3d(
    bds: &DataStructure,
    e: &TopoShape,
    collect_all: bool,
) -> (bool, Vec<TopoShape>) {
    let lssd = bds.shape_same_domain(e);
    let mut out = Vec::new();
    if lssd.is_empty() {
        return (false, out);
    }
    let mut hsd3d = false;
    for esd in lssd.iter() {
        let c = bds.same_domain_ori_of(esd);
        if c == Config::UnshGeometry {
            hsd3d = true;
            out.push(esd.clone());
            if !collect_all {
                break;
            }
        }
    }
    (hsd3d, out)
}

/// `FDS_Config3d(E1, E2, c)` (`cxx:539-575`): `SAMEORIENTED` /
/// `DIFFORIENTED` from the sign of the dot product of the two 3-D tangents at
/// the midpoint of `E1` and its projection on `E2`. `None` is `cxx:563-566`
/// (`!ok1 || !ok2`).
pub fn fds_config3d(e1: &TopoShape, e2: &TopoShape) -> Option<Config> {
    let (pe1, par_e1) = crate::topopebrep::tool::fun_tool_find_p_in_e(e1)?;
    let ve1 = crate::topopebrep::tool::tg_geom_e(par_e1, e1)?;
    let (par_e2, _d_e2) = crate::topopebrep::tool::fun_tool_proj_pon_e(&pe1, e2)?;
    let ve2 = crate::topopebrep::tool::tg_geom_e(par_e2, e2)?;
    let de1 = GpDir::from_vec(&ve1).ok()?;
    let de2 = GpDir::from_vec(&ve2).ok()?;
    let dot = de1.dot(&de2);
    Some(if dot > 0.0 {
        Config::SameOriented
    } else {
        Config::DiffOriented
    })
}

/// `FDS_HasSameDomain2d(BDS, E, PLSD)` (`cxx:577-615`): the 2-D counterpart,
/// keeping partners whose `SameDomainOri` is `SAMEORIENTED` or `DIFFORIENTED`.
pub fn fds_has_same_domain2d(
    bds: &DataStructure,
    e: &TopoShape,
    collect_all: bool,
) -> (bool, Vec<TopoShape>) {
    let lssd = bds.shape_same_domain(e);
    let mut out = Vec::new();
    if lssd.is_empty() {
        return (false, out);
    }
    let mut hsd2d = false;
    for esd in lssd.iter() {
        let c = bds.same_domain_ori_of(esd);
        if c == Config::SameOriented || c == Config::DiffOriented {
            hsd2d = true;
            out.push(esd.clone());
            if !collect_all {
                break;
            }
        }
    }
    (hsd2d, out)
}

/// `FUN_ds_ONesd(BDS, IE, EspON, IEsd)` (`cxx:3308-3345`): the same-domain
/// partner of `E = BDS.Shape(IE)` that carries the ON part `EspON`.
///
/// The probe point is taken at `x = 0.456789` of `EspON`'s range
/// (`cxx:3316-3319`); a partner is kept when the probe projects on it within
/// `1.e3 * tolerance` (`cxx:3332-3339`). The answer is the DS index of that
/// partner, written through `iesd`.
pub fn fun_ds_on_esd(
    bds: &DataStructure,
    ie: i32,
    esp_on: &TopoShape,
    iesd: &mut i32,
) -> bool {
    // `cxx:3313`: `BDS.Shape(IE)` — the default `findWithoutKeep = true`.
    let Some(e) = bds.shape(ie, true) else {
        return false;
    };
    let (f, l) = crate::topopebrep::tool::fun_tool_bounds(esp_on);
    let x = 0.456789;
    let par = (1.0 - x) * f + x * l;
    let Some(p3d) = crate::topopebrep::tool::fun_tool_value(par, esp_on) else {
        return false;
    };

    let lesd = bds.shape_same_domain(&e).to_vec();
    for esd in lesd {
        let Some((_parp, d)) = crate::topopebrep::tool::fun_tool_proj_pon_e(&p3d, &esd) else {
            continue;
        };
        let tolesd = crate::brep_tool::BRepTool::edge_tolerance(&crate::shape::Edge(esd.clone()));
        if !(d < tolesd * 1.0e3) {
            continue;
        }
        // `cxx:3341`: `IEsd = BDS.Shape(esd)` — again `findWithoutKeep = true`.
        *iesd = bds.shape_index(&esd, true);
        return true;
    }
    false
}

/// `FDS_LOIinfsup(BDS, E, pE, KDS, GDS, LOI, pbef, paft, isonboundper)`
/// (`cxx:3427-3566`): the parameters `pbef < pE < paft` bracketing `pE` on
/// `E`, taken from the interferences of `LOI` on a vertex/point other than
/// `(KDS, GDS)` whose transition is `(FACE, FACE)`.
///
/// The C++ tests `n == 0` *twice* (`cxx:3444`, `:3503`): the second one is
/// after the grouping loop and skips to the same answer. Both return `true`
/// (`cxx:3447`, `:3505`), so the port's `Err` only covers the
/// `FDS_Parameter` throw (`cxx:3508`, `:3555`).
pub fn fds_loi_inf_sup(
    e: &TopoShape,
    p_e: f64,
    k_ds: Kind,
    g_ds: i32,
    loi: &InterferenceList,
) -> Result<(f64, f64, bool), String> {
    let (f, l) = crate::topopebrep::tool::fun_tool_bounds(e);
    let mut pbef = f;
    let mut paft = l;
    if loi.is_empty() {
        return Ok((pbef, paft, false));
    }

    let mut loi_sans_gds: InterferenceList = Vec::new();
    let mut tki = TKI::new();
    tki.fill_on_geometry(loi);
    tki.init();
    while tki.more() {
        let group = tki.value().map(|(k, g, l)| (k, g, l.clone()));
        tki.next();
        let Some((k, g, list)) = group else {
            continue;
        };
        let pv = k == Kind::Point || k == Kind::Vertex;
        if !pv {
            continue;
        }
        let mk = k == k_ds;
        let mg = g == g_ds;
        if mk && mg {
            continue;
        }
        for h in list.iter() {
            let keep = {
                let i = h.read().expect("interference lock");
                let (tsb1, _isb1, tsa1, _isa1, _gt1, _g1, _st1, _s1) = fds_idata(&i);
                if tsb1 != crate::abs::ShapeType::Face {
                    false
                } else if tsa1 != crate::abs::ShapeType::Face {
                    false
                } else {
                    matches!(
                        i.data(),
                        super::interference::InterferenceData::CurvePoint { .. }
                            | super::interference::InterferenceData::EdgeVertex(_)
                    )
                }
            };
            if keep {
                loi_sans_gds.push(h.clone());
                break;
            }
        }
    }

    if loi_sans_gds.is_empty() {
        return Ok((pbef, paft, false));
    }

    let v = crate::topopebrep::tool::closed_e(e);
    let e_closed = v.is_some();
    let tole = crate::brep_tool::BRepTool::edge_tolerance(&crate::shape::Edge(e.clone()));
    let mut tol = occt_core::precision::Precision::parametric(tole, 1.0e2);
    let mut isonboundper = false;
    if e_closed {
        let v = v.expect("closed_e returned a vertex");
        let tolv = crate::brep_tool::BRepTool::vertex_tolerance(&crate::shape::Vertex(v));
        let tolv = occt_core::precision::Precision::parametric(tolv, 1.0e2);
        if tolv > tol {
            tol = tolv;
        }
        let p_e_is_ef = (p_e - f).abs() <= tol;
        let p_e_is_el = (p_e - l).abs() <= tol;
        isonboundper = p_e_is_ef || p_e_is_el;
    }

    for h in loi_sans_gds.iter() {
        let p = fds_parameter_or_fail(&h.read().expect("interference lock"))?;
        if isonboundper {
            // pbef > paft
            if p > pbef {
                pbef = p;
            }
            if p < paft {
                paft = p;
            }
        } else {
            if p > pbef && p < p_e {
                pbef = p;
            }
            if p < paft && p > p_e {
                paft = p;
            }
        }
    }
    Ok((pbef, paft, isonboundper))
}

/// `FDS_parbefaft(BDS, E, pE, pbef, paft, isonboundper, p1, p2)`
/// (`cxx:3572-3596`): the two interior parameters `p1` in `[pbef, pE[` and
/// `p2` in `]pE, paft]` with the fixed ratio `t = 0.3178947713`.
pub fn fds_parbefaft(e: &TopoShape, p_e: f64, pbef: f64, paft: f64, isonboundper: bool) -> (f64, f64) {
    let t = 0.3178947713;
    let (f, l) = crate::topopebrep::tool::fun_tool_bounds(e);
    if isonboundper {
        ((1.0 - t) * pbef + t * l, (1.0 - t) * f + t * paft)
    } else {
        ((1.0 - t) * pbef + t * p_e, (1.0 - t) * p_e + t * paft)
    }
}

/// `FDS_stateEwithF2d(BDS, E, pE, KDS, GDS, F1, TrmemeS)`
/// (`cxx:3598-3643`): the 2d transition of `E` on `F1` at `pE`, i.e. the
/// states of the two bracketing points.
///
/// `ON` is mapped to `IN` (`cxx:3629-3636`, "xpu190898 : cto014I2").
pub fn fds_state_e_with_f2d(
    bds: &DataStructure,
    e: &TopoShape,
    p_e: f64,
    k_ds: Kind,
    g_ds: i32,
    f1: &TopoShape,
    tr_meme_s: &mut Transition,
) -> Result<bool, String> {
    let loi = bds.shape_interferences(e, true).to_vec();
    let (pbef, paft, isonper) = fds_loi_inf_sup(e, p_e, k_ds, g_ds, &loi)?;
    let (t1, t2) = fds_parbefaft(e, p_e, pbef, paft, isonper);
    let (Some(p1), Some(p2)) = (
        crate::topopebrep::tool::fun_tool_value(t1, e),
        crate::topopebrep::tool::fun_tool_value(t2, e),
    ) else {
        return Ok(false);
    };

    let mut psc = crate::topopebrep::tool::fsc_get_psc(f1);
    let mut sta1 = crate::topopebrep::tool::fsc_state_pon_face(&p1, f1, &mut psc);
    let mut sta2 = crate::topopebrep::tool::fsc_state_pon_face(&p2, f1, &mut psc);

    // xpu190898 : cto014I2 (e5,fTRASHA14,vG7,eS10)
    if sta1 == crate::abs::State::On {
        sta1 = crate::abs::State::In;
    }
    if sta2 == crate::abs::State::On {
        sta2 = crate::abs::State::In;
    }

    tr_meme_s.before(sta1, crate::abs::ShapeType::Face);
    tr_meme_s.after(sta2, crate::abs::ShapeType::Face);
    Ok(true)
}

/// `FDS_getupperlower(HDS, edgeIndex, paredge, p1, p2)`
/// (`TopOpeBRepDS_EXPORT.cxx:617-643`): `p1`/`p2` start at the edge bounds and
/// are narrowed to the point-interference parameters bracketing `paredge`.
///
/// The C++ `EdgePoints(E)` walk yields only point/vertex interferences, whose
/// `Parameter()` never raises; the port still propagates the `Err` of the four
/// other interference types in case a list holds one.
pub fn fds_getupperlower(
    hds: &HDataStructure,
    e: &TopoShape,
    paredge: f64,
) -> Result<(f64, f64), String> {
    let (mut p1, mut p2) = crate::topopebrep::tool::fun_tool_bounds(e);
    let mut it = hds.edge_points(e);
    while it.more() {
        let par = it.parameter()?;
        it.next();
        let parsup1 = par > p1;
        let parinfe = par < paredge;
        let parinf2 = par < p2;
        let parsupe = par > paredge;
        if parsup1 && parinfe {
            p1 = par;
        }
        if parinf2 && parsupe {
            p2 = par;
        }
    }
    Ok((p1, p2))
}

// UNPORTED in `TopOpeBRepDS_EXPORT.cxx` (all call sites outside this crate's
// reached packages, or dead in OCCT itself):
// - `FUN_ds_hasFEI` (`:3279-3306`) and `FUN_ds_ONesd` (`:3308-3345`): the only
//   callers are `TopOpeBRepBuild_GridFF.cxx:227/246/281` and
//   `TopOpeBRepBuild_BuilderON.cxx:672`, i.e. the `TopOpeBRepBuild` package.
// - `FUN_ds_unkeepEVIonGb1` (`:3364-3425`) and `FUN_ds_FEIGb1TO0`
//   (`:3647-3774`): their only call sites are commented out in OCCT
//   (`TopOpeBRep_DSFiller.cxx:768`, `TopOpeBRepDS_Reducer.cxx:44`), and the
//   first one's body is itself inside `/* */`.
//
// The `FUN_ds_completefor*` / `FUN_ds_PURGEforE9` / `FUN_ds_redusamsha` /
// `FUN_ds_PointToVertex` / `FUN_ds_complete1dForSESDM` family (`:715-3952`)
// lives in [`super::complete_ds`]; the 1d/2d reduction helpers of
// `TopOpeBRepDS_redu.cxx` live in [`super::redu`].

/// `FUN_EisSE2` return codes (`cxx:1792-1795`).
pub const EISSE2_NONE: i32 = 0;
/// `EREFandESUP`
pub const EISSE2_EREF_AND_ESUP: i32 = 2;

/// `FUN_EisSE2(BDS, Fsdm, E, I)` (`cxx:1797-1871`).
///
/// `Fsdm` is accepted but unused by OCCT (`cxx:1800`).
pub fn fun_eis_se2(bds: &DataStructure, e: &TopoShape, i: &Interference) -> i32 {
    let (sb, _ib, _sa, _ia, gt, g, _st, _s) = fds_idata(i);
    if gt != Kind::Vertex {
        return EISSE2_NONE;
    }
    // `cxx:1811-1818`: the orientation of the transition's state before, and
    // no `EXTERNAL` interference.
    let o = i
        .transition()
        .orientation(crate::abs::State::In, crate::abs::ShapeType::Shape);
    if o == crate::abs::Orientation::External {
        return EISSE2_NONE;
    }
    // `cxx:1820-1825`: E same-domain edge
    if sb != crate::abs::ShapeType::Edge {
        return EISSE2_NONE;
    }

    let e_index = bds.shape_index(e, true);
    let rk_e = bds.ancestor_rank(e_index);
    let rk_g = bds.ancestor_rank(g);
    let Some(vg) = bds.shape(g, true) else {
        return EISSE2_NONE;
    };
    let gsd = super::export::fun_ds_get_vsdm(bds, g);
    let o_g_in_e = if let Some(gsd) = gsd {
        let vgsd = bds.shape(gsd, true).unwrap_or_else(|| vg.clone());
        if rk_e == rk_g {
            crate::topopebrep::tool::fun_tool_orient_vin_e(&vg, e)
        } else {
            crate::topopebrep::tool::fun_tool_orient_vin_e(&vgsd, e)
        }
    } else {
        crate::topopebrep::tool::fun_tool_orient_vin_e(&vg, e)
    };

    // `cxx:1850-1853`
    if o_g_in_e == 0 {
        return EISSE2_EREF_AND_ESUP; // G in IN E
    }
    // `cxx:1855-1861`
    let noshare = (o_g_in_e == crate::topopebrep::tool::LAST
        && o == crate::abs::Orientation::Forward)
        || (o_g_in_e == crate::topopebrep::tool::FIRST
            && o == crate::abs::Orientation::Reversed);
    if noshare {
        return EISSE2_NONE;
    }
    // `cxx:1863-1869`: E same-domain edge(ST); both arms answer
    // `EREFandESUP`.
    EISSE2_EREF_AND_ESUP
}

/// `FUN_ds_FillSDMFaces(HDS)` (`cxx:1874-1947`): every edge of a face that has
/// a same-domain face of the other ancestor rank becomes a section edge when
/// [`fun_eis_se2`] accepts one of its interferences.
pub fn fun_ds_fill_sdm_faces(hds: &mut HDataStructure) -> Result<(), String> {
    let ns = hds.change_ds().nb_shapes();
    for i in 1..=ns {
        let Some(s) = hds.ds().shape(i, true) else {
            continue;
        };
        if s.shape_type() != crate::abs::ShapeType::Face {
            continue;
        }
        let rk_s = hds.ds().ancestor_rank(i);
        // `cxx:1888-1899`: the same-domain faces of another rank.
        let mut fsdm: Vec<TopoShape> = Vec::new();
        for f in hds.ds().shape_same_domain(&s) {
            let rkf = hds.ds().ancestor_rank_of(f);
            if rkf != rk_s {
                fsdm.push(f.clone());
            }
        }
        if fsdm.is_empty() {
            continue;
        }

        let mut edges: Vec<TopoShape> = Vec::new();
        let mut ex = crate::topexp::Explorer::new(&s, crate::abs::ShapeType::Edge);
        while ex.more() {
            edges.push(ex.current().clone());
            ex.next();
        }
        for e in edges {
            // `cxx:1908-1916`
            if !hds.ds().has_shape(&e, true) {
                continue;
            }
            if hds.ds().is_section_edge(&e, true) {
                continue;
            }
            let li = hds.ds().shape_interferences(&e, true).to_vec();
            if li.len() < 1 {
                continue;
            }
            // `cxx:1926-1943`
            for handle in &li {
                let itf = handle.read().expect("interference lock").clone();
                let is_se2 = fun_eis_se2(hds.ds(), &e, &itf);
                if is_se2 == EISSE2_NONE {
                    continue;
                }
                hds.change_ds().add_section_edge(&e);
                if is_se2 == EISSE2_EREF_AND_ESUP {
                    let st = itf.support();
                    if let Some(es) = hds.ds().shape(st, true) {
                        hds.change_ds().add_section_edge(&es);
                    }
                }
                break;
            }
        }
    }
    Ok(())
}

/// `FUN_ds_addSEsdm1d(HDS)` (`cxx:1949-1995`): every non-degenerated, non
/// section edge of rank 1 that shares geometry with a same-domain edge becomes
/// a section edge along with its `EShareG` partners.
pub fn fun_ds_add_sesdm1d(hds: &mut HDataStructure) -> Result<(), String> {
    let ns = hds.change_ds().nb_shapes();
    for i in 1..=ns {
        let Some(s) = hds.ds().shape(i, true) else {
            continue;
        };
        if s.shape_type() != crate::abs::ShapeType::Edge {
            continue;
        }
        // `cxx:1961-1965`
        if crate::brep_tool::BRepTool::is_degenerated(&crate::shape::Edge(s.clone())) {
            continue;
        }
        // `cxx:1967-1971`
        if hds.ds().is_section_edge(&s, true) {
            continue;
        }
        // `cxx:1972-1976`
        if hds.ds().ancestor_rank(i) != 1 {
            continue;
        }
        let mut lsd: Vec<TopoShape> = Vec::new();
        super::tool::e_share_g(hds, &s, &mut lsd);
        let share_g = !lsd.is_empty();
        for e in &lsd {
            hds.change_ds().add_section_edge(e);
        }
        if share_g {
            hds.change_ds().add_section_edge(&s);
        }
    }
    Ok(())
}

/// `FDS_SIisGIofIofSBAofTofI(BDS, SI, I)` (`cxx:396-462`): the face that `I`
/// names (before side first, then after side) carries an interference whose
/// geometry is the edge `SI`.
pub fn fds_si_is_gi_of_i_of_sba_of_t_of_i(
    bds: &DataStructure,
    si: i32,
    i: &Interference,
) -> bool {
    if si == 0 {
        return false;
    }
    let (sb1, ib1, sa1, ia1, _gt1, _g1, _st1, _s1) = fds_idata(i);

    let face_index = if sb1 == crate::abs::ShapeType::Face {
        ib1
    } else if sa1 == crate::abs::ShapeType::Face {
        ia1
    } else {
        return false;
    };

    for ib in bds.shape_interferences_at(face_index, true) {
        let (gtb, gb) = {
            let other = ib.read().expect("interference lock");
            let (_sbb, _ibb, _sab, _iab, gtb, gb, _stb, _sb) = fds_idata(&other);
            (gtb, gb)
        };
        if gtb == Kind::Edge && gb == si {
            // the face carries an interference whose geometry is the edge SI
            return true;
        }
    }
    false
}

