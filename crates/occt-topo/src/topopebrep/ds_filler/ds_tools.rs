//! The package-level helpers `TopOpeBRep_DSFiller` runs at the end of its
//! insert passes.
//!
//! Source: `TKBool/TopOpeBRep/TopOpeBRep_DSFiller.cxx:174-424`,
//! `TopOpeBRep_sort.cxx:26-104`, `TopOpeBRep_mergePDS.cxx:29-196`
//! (OCCT 8.0.0).
//!
//! * `BREP_correctgbound` (`DSFiller.cxx:174-278`) re-derives the `GBound`
//!   flag of every `VERTEX`-geometry shape-shape interference of an edge from
//!   the edge's own vertex map, then rewrites the geometry index to the
//!   same-domain reference vertex.
//! * `BREP_UnfillSameDomain` (`:282-300`) drops a same-domain couple the
//!   classifier reports as fully `OUT` or fully `UNKNOWN`.
//! * `FUN_shareNOG` (`:305-424`) decides whether two same-domain faces share
//!   no geometric domain.
//! * `BREP_sortonparameter` (`sort.cxx:90-104`) groups an edge's interferences
//!   by geometry and orders the groups by the parameter of their first
//!   `CurvePointInterference`.
//! * `BREP_mergePDS` (`mergePDS.cxx:61-196`) turns curve-point interferences
//!   whose point coincides with a vertex of one of the curve's faces into
//!   vertex interferences.
//!
//! The `qsort` of `sort.cxx:81` is not stable; the port uses a stable sort, so
//! groups whose comparator answers `0` (`sort.cxx:27-60`: a first element that
//! is not a `CurvePointInterference`) keep their relative order.

use occt_core::gp::{GpDir, GpVec};
use occt_core::precision::ANGULAR;

use crate::abs::{ShapeType, State};
use crate::shape::{Face, TopoShape};
use crate::topexp::Explorer;
use crate::topo_tools::is_same;
use crate::topopebrep::ds::curve::Curve;
use crate::topopebrep::ds::curve_explorer::CurveExplorer;
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::interference_iterator::InterferenceList;
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::ds::point::Point;
use crate::topopebrep::ds::tki::TKI;
use crate::topopebrep::tool::{
    fun_tool_ng_geom_f, fun_tool_par_e, fun_tool_par_on_e_f, fun_tool_tg_geom_e, get_geom_xx,
    ShapeClassifier,
};

/// `BREP_correctgbound(HDS)` (`DSFiller.cxx:174-278`).
pub fn brep_correct_g_bound(hds: &mut HDataStructure) {
    let n = hds.ds().nb_shapes();
    for i in 1..=n {
        // `cxx:180-186`
        let Some(s) = hds.ds().shape(i, true) else {
            continue;
        };
        if s.shape_type() != ShapeType::Edge {
            continue;
        }
        let interferences = hds.ds().shape_interferences(&s, true).to_vec();
        // `cxx:188-191`
        if interferences.is_empty() {
            continue;
        }

        // `cxx:192-212`: does the edge carry a shape-shape interference whose
        // geometry is a VERTEX?
        let mut ehassiv = false;
        for handle in &interferences {
            let Ok(itf) = handle.read() else {
                continue;
            };
            if itf.shape_shape_data().is_none() {
                continue;
            }
            if itf.geometry_type() != Kind::Vertex {
                continue;
            }
            ehassiv = true;
            break;
        }
        if !ehassiv {
            continue;
        }

        // `cxx:218-219`: `TopExp::MapShapes(s, TopAbs_VERTEX, imev)`.
        let imev = vertices_of(&s);
        let contains = |v: &TopoShape| imev.iter().any(|x| is_same(x, v));

        // `cxx:221-241`: GBound of every VERTEX geometry.
        for handle in &interferences {
            let Ok(mut itf) = handle.write() else {
                continue;
            };
            if itf.shape_shape_data().is_none() {
                continue;
            }
            if itf.geometry_type() != Kind::Vertex {
                continue;
            }
            let Some(v) = hds.ds().shape(itf.geometry(), true) else {
                continue;
            };
            let vofe = contains(&v);
            itf.set_g_bound(vofe);
        }

        // `cxx:244-276`: replace the geometry by its same-domain reference.
        for handle in &interferences {
            let Ok(mut itf) = handle.write() else {
                continue;
            };
            if itf.shape_shape_data().is_none() {
                continue;
            }
            if itf.geometry_type() != Kind::Vertex {
                continue;
            }
            let gi = itf.geometry();
            let Some(v) = hds.ds().shape(gi, true) else {
                continue;
            };
            // `cxx:261-265`
            if !hds.has_same_domain(&v, true) {
                continue;
            }
            let ivref = hds.ds().same_domain_ref(gi);
            if ivref == gi {
                continue;
            }
            let Some(vref) = hds.ds().shape(ivref, true) else {
                continue;
            };
            let vrefofe = contains(&vref);
            itf.set_geometry(ivref);
            itf.set_g_bound(vrefofe);
        }
    }
}

/// `BREP_UnfillSameDomain(F1, F2, HDS, SC)` (`DSFiller.cxx:282-300`).
pub fn brep_unfill_same_domain(
    f1: &TopoShape,
    f2: &TopoShape,
    hds: &mut HDataStructure,
    sc: &mut ShapeClassifier,
) -> Result<bool, String> {
    let samdom = 1;
    let st1 = sc.state_shape_shape(f1, f2, samdom)?;
    let st2 = sc.state_shape_shape(f2, f1, samdom)?;
    let unfill = (st1 == State::Out && st2 == State::Out)
        || (st1 == State::Unknown && st2 == State::Unknown);
    if unfill {
        hds.change_ds().unfill_shapes_same_domain(f1, f2);
    }
    Ok(unfill)
}

/// `FUN_shareNOG(HDS, lFF1, lFF2)` (`DSFiller.cxx:305-424`).
pub fn fun_share_nog(hds: &HDataStructure, lff1: &TopoShape, lff2: &TopoShape) -> bool {
    let bds = hds.ds();
    let f1 = Face(lff1.clone());
    let f2 = Face(lff2.clone());
    let map1 = edges_of(lff1);
    let map2 = edges_of(lff2);
    let tola = ANGULAR;

    // `cxx:317-333`
    for handle in bds.shape_interferences(lff1, true) {
        let Ok(itf) = handle.read() else {
            continue;
        };
        if itf.geometry_type() != Kind::Edge {
            continue;
        }
        let Some(eg) = bds.shape(itf.geometry(), true) else {
            continue;
        };
        if map2.iter().any(|e| is_same(e, &eg)) {
            return false;
        }
    }
    // `cxx:336-350`
    for handle in bds.shape_interferences(lff2, true) {
        let Ok(itf) = handle.read() else {
            continue;
        };
        if itf.geometry_type() != Kind::Edge {
            continue;
        }
        let Some(eg) = bds.shape(itf.geometry(), true) else {
            continue;
        };
        if map1.iter().any(|e| is_same(e, &eg)) {
            return false;
        }
    }

    // `cxx:352-421`
    let mut ex1 = Explorer::new(lff1, ShapeType::Edge);
    while ex1.more() {
        let e1 = ex1.current().clone();
        // `cxx:356-367`
        if hds.has_same_domain(&e1, true) {
            for e2 in hds.same_domain(&e1) {
                if map2.iter().any(|e| is_same(e, e2)) {
                    return false;
                }
            }
        }
        // `cxx:370-420`
        for handle in bds.shape_interferences(&e1, true) {
            let Ok(itf) = handle.read() else {
                continue;
            };
            // `cxx:377-380`
            if itf.support_type() != Kind::Edge {
                continue;
            }
            let Some(e2) = bds.shape(itf.support(), true) else {
                continue;
            };
            if !map2.iter().any(|e| is_same(e, &e2)) {
                continue;
            }
            // `cxx:388-390`
            if itf.geometry_type() != Kind::Point || !itf.is_curve_point() {
                continue;
            }
            let Some(par1) = itf.parameter() else {
                continue;
            };
            // `cxx:395-403`
            let Some(par2) = fun_tool_par_e(&e1, par1, &e2) else {
                return false;
            };
            // `cxx:405-410`
            let tge1 = fun_tool_tg_geom_e(par1, &e1);
            let tge2 = fun_tool_tg_geom_e(par2, &e2);
            let (Ok(d1), Ok(d2)) = (GpDir::from_vec(&tge1), GpDir::from_vec(&tge2)) else {
                return false;
            };
            let dot: f64 = d1.dot(&d2);
            let x = (1.0f64 - dot.abs()).abs();
            if x > tola {
                return false; // e1,e2 not tangent
            }
            // `cxx:412-417`
            let xx_f1 = get_geom_xx_on_face(lff1, &e1, par1);
            let xx_f2 = get_geom_xx_on_face(lff2, &e2, par2);
            let (Ok(dx1), Ok(dx2)) = (GpDir::from_vec(&xx_f1), GpDir::from_vec(&xx_f2)) else {
                return false;
            };
            if dx1.dot(&dx2) > 0.0 {
                return false; // F1,F2 share geometric domain near G1
            }
        }
        ex1.next();
    }
    let _ = (&f1, &f2);
    true
}

/// `FUN_tool_getgeomxx(Fi, Ei, parOnEi)` (`TopOpeBRepTool_matter.cxx:89-104`):
/// the 3-argument overload. A failed `FUN_tool_paronEF` keeps the initial
/// `gp_Vec(1,0,0)` (`cxx:91-98`); the 4-argument overload answers a null
/// vector when the tangent lookup fails (`cxx:76-86`).
fn get_geom_xx_on_face(f: &TopoShape, e: &TopoShape, par: f64) -> GpVec {
    let mut xx = GpVec::new(1.0, 0.0, 0.0);
    let Some(uv) = fun_tool_par_on_e_f(e, par, f) else {
        return xx;
    };
    let Some(ng) = fun_tool_ng_geom_f(&uv, &Face(f.clone())) else {
        return GpVec::new(0.0, 0.0, 0.0);
    };
    if let Some(v) = get_geom_xx(e, par, &ng) {
        xx = v;
    } else {
        xx = GpVec::new(0.0, 0.0, 0.0);
    }
    xx
}

/// `BREP_sortonparameter(HDS)` (`sort.cxx:90-104`).
pub fn brep_sort_on_parameter(hds: &mut HDataStructure) {
    let n = hds.ds().nb_shapes();
    for i in 1..=n {
        let Some(s) = hds.ds().shape(i, true) else {
            continue;
        };
        if s.shape_type() != ShapeType::Edge {
            continue;
        }
        let mut list = hds.ds().shape_interferences(&s, true).to_vec();
        brep_sort_on_parameter_2(&mut list);
        *hds.change_ds().change_shape_interferences(&s) = list;
    }
}

/// `compll(l1, l2)` (`sort.cxx:27-60`): the parameter of the group's first
/// element when that element is a `CurvePointInterference`; otherwise the
/// groups are incomparable.
fn group_parameter(list: &InterferenceList) -> Option<f64> {
    let first = list.first()?;
    let itf = first.read().ok()?;
    if !itf.is_curve_point() {
        return None;
    }
    itf.parameter()
}

/// `BREP_sortonparameter2(LOI)` (`sort.cxx:62-88`).
fn brep_sort_on_parameter_2(loi: &mut InterferenceList) {
    let mut tki = TKI::new();
    tki.fill_on_geometry(loi);
    let mut groups: Vec<(Option<f64>, InterferenceList)> = Vec::new();
    tki.init();
    while tki.more() {
        let list = tki.value().map(|(_, _, l)| l.clone());
        if let Some(list) = list {
            groups.push((group_parameter(&list), list));
        }
        tki.next();
    }
    groups.sort_by(|a, b| match (a.0, b.0) {
        (Some(x), Some(y)) => x.partial_cmp(&y).unwrap_or(std::cmp::Ordering::Equal),
        _ => std::cmp::Ordering::Equal,
    });
    loi.clear();
    for (_, list) in groups {
        loi.extend(list);
    }
}

/// `BREP_mergePDS(HDS)` (`mergePDS.cxx:61-196`).
pub fn brep_merge_pds(hds: &mut HDataStructure) -> Result<(), String> {
    // `cxx:64-68`: an empty curve explorer returns immediately.
    let curves: Vec<(i32, Option<Curve>)> = {
        let bds = hds.ds();
        let mut cex = CurveExplorer::with_ds(bds, false);
        let mut out = Vec::new();
        while cex.more() {
            let ic = cex.index();
            // `cxx:78-82`: `!ChangeCurveInterferences(ic).More()` -> continue
            if bds.curve_interferences(ic).is_empty() {
                cex.next();
                continue;
            }
            out.push((ic, cex.curve().cloned()));
            cex.next();
        }
        out
    };

    for (ic, curve) in curves {
        let curve = curve.ok_or_else(|| "BREP_mergePDS : no curve".to_string())?;
        // `cxx:85-97`
        let f1 = curve
            .shape1()
            .cloned()
            .ok_or_else(|| "BREP_mergePDS : no shape1".to_string())?;
        let f2 = curve
            .shape2()
            .cloned()
            .ok_or_else(|| "BREP_mergePDS : no shape2".to_string())?;
        let mvp1 = make_idmovp(&f1);
        let mvp2 = make_idmovp(&f2);

        let list = hds.ds().curve_interferences(ic).to_vec();
        for handle in list {
            let Ok(mut itf) = handle.write() else {
                continue;
            };
            // `cxx:104-107`: `occ::down_cast<TopOpeBRepDS_CurvePointInterference>`
            if !itf.is_curve_point() {
                continue;
            }
            // `cxx:110-113`
            if itf.geometry_type() != Kind::Point {
                continue;
            }
            // `cxx:115-121`
            let gi = itf.geometry();
            if gi > hds.ds().nb_points() {
                continue;
            }
            let Some(pds) = hds.ds().point(gi).cloned() else {
                continue;
            };

            // `cxx:125-137`
            let ivp1 = find_pds_among_idmovp(&pds, &mvp1);
            let (k1, iv1) = match ivp1 {
                Some(idx) => {
                    let v = mvp1[idx].0.clone();
                    let i = hds.change_ds().add_shape(&v);
                    (Kind::Vertex, i)
                }
                None => (Kind::Unknown, 0),
            };
            // `cxx:138-148`
            let ivp2 = find_pds_among_idmovp(&pds, &mvp2);
            let (k2, iv2) = match ivp2 {
                Some(idx) => {
                    let v = mvp2[idx].0.clone();
                    let i = hds.change_ds().add_shape(&v);
                    (Kind::Vertex, i)
                }
                None => (Kind::Unknown, 0),
            };

            // `cxx:150-153`
            if let (Some(i1), Some(i2)) = (ivp1, ivp2) {
                hds.change_ds()
                    .fill_shapes_same_domain(&mvp1[i1].0, &mvp2[i2].0, true)?;
            }

            // `cxx:155-168`
            if ivp1.is_some() {
                itf.set_geometry_type(k1);
                itf.set_geometry(iv1);
            } else if ivp2.is_some() {
                itf.set_geometry_type(k2);
                itf.set_geometry(iv2);
            }
        }
    }
    Ok(())
}

/// `BREP_makeIDMOVP(S, IDMOVP)` (`mergePDS.cxx:47-59`): the vertices of `S`
/// with their `TopOpeBRepDS_Point`, in `TopExp_Explorer` order.
fn make_idmovp(s: &TopoShape) -> Vec<(TopoShape, Point)> {
    let mut out = Vec::new();
    let mut ex = Explorer::new(s, ShapeType::Vertex);
    while ex.more() {
        let v = ex.current().clone();
        // `TopOpeBRepDS_Point PDS(v)` (`mergePDS.cxx:56`): built on the vertex
        // alone, independent of the DS.
        match Point::of_shape(&v) {
            Ok(p) => out.push((v, p)),
            Err(_) => break,
        }
        ex.next();
    }
    out
}

/// `BREP_findPDSamongIDMOVP(PDS, IDMOVP)` (`mergePDS.cxx:29-45`): the map
/// position (1-based in OCCT, 0-based here) of the first point equal to
/// `PDS`. `None` is the `0` answer.
fn find_pds_among_idmovp(pds: &Point, idmovp: &[(TopoShape, Point)]) -> Option<usize> {
    idmovp.iter().position(|(_, pm)| pds.is_equal(pm))
}

/// `TopExp::MapShapes(S, TopAbs_VERTEX)` as a list.
fn vertices_of(s: &TopoShape) -> Vec<TopoShape> {
    let mut out = Vec::new();
    let mut ex = Explorer::new(s, ShapeType::Vertex);
    while ex.more() {
        out.push(ex.current().clone());
        ex.next();
    }
    out
}

/// `TopExp::MapShapes(S, TopAbs_EDGE)` as a list.
fn edges_of(s: &TopoShape) -> Vec<TopoShape> {
    let mut out = Vec::new();
    let mut ex = Explorer::new(s, ShapeType::Edge);
    while ex.more() {
        out.push(ex.current().clone());
        ex.next();
    }
    out
}
