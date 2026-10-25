//! `TopOpeBRep_ProcessSectionEdges.cxx` — the section-edge pass of
//! `TopOpeBRep_FacesFiller::Insert`.
//!
//! Sources (OCCT 8.0.0, `ModelingAlgorithms/TKBool/TopOpeBRep`):
//! - `GetESL` (`TopOpeBRep_ProcessSectionEdges.cxx:54-103`)
//! - `KeepRLine` (`ProcessSectionEdges.cxx:107-283`)
//! - `FUN_brep_sdmRE` (`ProcessSectionEdges.cxx:285-330`)
//! - `ProcessSectionEdges` (`ProcessSectionEdges.cxx:334-468`)
//!
//! Translation notes:
//! - `KeepRLine` is `const` in OCCT but mutates the shape classifier through
//!   `myPShapeClassifier`; the port takes `&mut self` and drives
//!   `StBipVPonF` through [`super::process_gr::st_bip_vp_on_f`].
//! - OCCT's `VPBounds(f, l, n)` is `const` and caches the bounds inside the
//!   line. The port's is `&mut`, so `KeepRLine` clones the line out of the
//!   intersector: every read in this pass is read-only, and OCCT's cache only
//!   serves the reads made inside this same function.
//! - `myDS->Shape(E)` / `myDS->AncestorRank(E)` (`cxx:354-355`) are pure
//!   lookups (`DataStructure.cxx:1159-1168`, `:...`); the port calls them and
//!   discards the result, as OCCT does.
//! - The `toappend` local at `cxx:385` is a `double` initialised to `true`;
//!   OCCT never changes it. The port keeps the unconditional append.
//! - `LES` / `LOI` are paired positionally at `cxx:399-403`: an edge that got
//!   neither rank shifts the pairing. That is OCCT's behaviour and is kept.
//! - `mapELE` is a `NCollection_DataMap` iterated in insertion order here (the
//!   OCCT `DataMap` order is unspecified but stable within a run); the port
//!   keeps a vector of `(key, values)`.

use crate::abs::State;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, TopoShape};
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::Config;
use crate::topopebrep::faces_intersector::FacesIntersector;
use crate::topopebrep::interf::VPointInter;
use crate::topopebrep::line_inter::LineInter;
use crate::topopebrep::tool::adaptor::CurveAdaptor;
use crate::topopebrep::tool::closed::{closed_e, vertices};
use crate::topopebrep::tool::project_curve::fun_tool_bounds;
use crate::topopebrep::tool::shape_tool::{pnt_of, tolerance_of};
use crate::topopebrep::tool::topology::fun_tool_pin_c;
use crate::topopebrep::type_line_curve::TypeLineCurve;

/// First `i >= from` whose `VPoint` passes the `checkkeep` filter
/// (`TopOpeBRep_VPointInterIterator::Init`, `cxx:49-70`).
fn first_kept(l: &LineInter, from: i32, check_keep: bool) -> Option<i32> {
    let nb = l.nb_v_point();
    let mut i = from;
    while i <= nb {
        let ok = match l.v_point(i) {
            Some(vp) => !check_keep || vp.keep(),
            None => false,
        };
        if ok {
            return Some(i);
        }
        i += 1;
    }
    None
}

/// `FUN_brep_sdmRE(E1, E2)` (`ProcessSectionEdges.cxx:285-330`): `true` when
/// both end vertices of `E1` or both end vertices of `E2` lie on the other
/// edge within `max(Tolerance(E), max(Tolerance(v3), Tolerance(v4)))`.
///
/// The four `BRep_Tool::Tolerance` calls are kept in OCCT's order and the same
/// tolerance is reused for both points of an edge (`cxx:294-309`, `:314-328`).
pub fn fun_brep_sdm_re(e1: &TopoShape, e2: &TopoShape) -> Result<bool, String> {
    let (v1, v2) = vertices(e1);
    let (v3, v4) = vertices(e2);
    let (Some(v1), Some(v2)) = (v1, v2) else {
        return Err("FUN_brep_sdmRE : E1 has no end vertex".to_string());
    };
    let (Some(v3), Some(v4)) = (v3, v4) else {
        return Err("FUN_brep_sdmRE : E2 has no end vertex".to_string());
    };

    let mut ok = false;

    let bac1 = CurveAdaptor::of_edge(&Edge(e1.clone()));
    let tol1 = BRepTool::edge_tolerance(&Edge(e1.clone()));
    let tol2 = tolerance_of(&v3)?;
    let tol3 = tolerance_of(&v4)?;
    let tol4 = tol1.max(tol2).max(tol3);
    if !ok {
        let p3 = pnt_of(&v3)?;
        ok = fun_tool_pin_c(&p3, &bac1, tol4);
    }
    if !ok {
        let p4 = pnt_of(&v4)?;
        ok = fun_tool_pin_c(&p4, &bac1, tol4);
    }

    if !ok {
        let bac2 = CurveAdaptor::of_edge(&Edge(e2.clone()));
        let tol1 = BRepTool::edge_tolerance(&Edge(e2.clone()));
        let tol2 = tolerance_of(&v1)?;
        let tol3 = tolerance_of(&v2)?;
        let tol4 = tol1.max(tol2).max(tol3);
        if !ok {
            let p1 = pnt_of(&v1)?;
            ok = fun_tool_pin_c(&p1, &bac2, tol4);
        }
        if !ok {
            let p2 = pnt_of(&v2)?;
            ok = fun_tool_pin_c(&p2, &bac2, tol4);
        }
    }
    Ok(ok)
}

impl super::FacesFiller {
    /// `GetESL(LES)` (`ProcessSectionEdges.cxx:54-103`): the restriction edges
    /// of the current intersector that `KeepRLine` accepts, without repeats.
    ///
    /// OCCT's `mapES` is a `TopTools_ShapeMapHasher` set (`IsSame`); the port
    /// compares with [`crate::topo_tools::is_same`].
    pub fn get_esl(&mut self, facint: &FacesIntersector) -> Result<Vec<TopoShape>, String> {
        let mut les: Vec<TopoShape> = Vec::new();
        let mut map_es: Vec<TopoShape> = Vec::new();

        for i in 1..=facint.nb_lines() {
            let l = facint
                .lines()
                .and_then(|h| h.value(i))
                .ok_or_else(|| "TopOpeBRep_FacesFiller::GetESL : no line".to_string())?;
            if l.type_line_curve() != TypeLineCurve::Restriction {
                continue;
            }
            let e = l
                .arc()
                .cloned()
                .ok_or_else(|| "TopOpeBRep_FacesFiller::GetESL : no arc".to_string())?;

            let mut add = !map_es.iter().any(|x| crate::topo_tools::is_same(x, &e));
            if add {
                add = self.keep_r_line(facint, i, false)?;
            }
            if add {
                map_es.push(e.clone());
                les.push(e);
            }
        }
        Ok(les)
    }

    /// `KeepRLine(L, checkkeep)` (`ProcessSectionEdges.cxx:107-283`): `true`
    /// when the restriction line `L` carries a kept edge portion with a
    /// `IN`/`ON` middle point on the other face.
    pub fn keep_r_line(
        &mut self,
        facint: &FacesIntersector,
        index: i32,
        check_keep: bool,
    ) -> Result<bool, String> {
        let mut l = facint
            .lines()
            .and_then(|h| h.value(index))
            .cloned()
            .ok_or_else(|| "TopOpeBRep_FacesFiller::KeepRLine : no line".to_string())?;

        if l.type_line_curve() != TypeLineCurve::Restriction {
            return Ok(false);
        }
        let el = l
            .arc()
            .cloned()
            .ok_or_else(|| "TopOpeBRep_FacesFiller::KeepRLine : no arc".to_string())?;
        if BRepTool::is_degenerated(&Edge(el.clone())) {
            return Ok(false);
        }

        let isedge1 = l.arc_is_edge(1);

        // VPI.Init(L, checkkeep)
        let Some(first) = first_kept(&l, 1, check_keep) else {
            return Ok(false);
        };
        let vpf = l
            .v_point(first)
            .cloned()
            .ok_or_else(|| "TopOpeBRep_FacesFiller::KeepRLine : no VPoint".to_string())?;

        let mut samevp = true;
        let mut vpl = VPointInter::new();

        let mut middle = false; // xpu011098 : cto012U1
        if closed_e(&el).is_some() {
            let (parf, parl) = fun_tool_bounds(&el);
            let mut i = first_kept(&l, 1, check_keep);
            while let Some(k) = i {
                let vp = l
                    .v_point(k)
                    .cloned()
                    .ok_or_else(|| "TopOpeBRep_FacesFiller::KeepRLine : no VPoint".to_string())?;
                let pf = super::process_gr::vp_param_on_er(&vp, &l)?;
                if parf < pf && pf < parl {
                    middle = true;
                    samevp = false;
                    vpl = vp;
                    break;
                }
                i = first_kept(&l, k + 1, check_keep);
            }
            if middle {
                // VPI.Init(L, checkkeep)
                let mut i = first_kept(&l, 1, check_keep);
                while let Some(k) = i {
                    let vp = l.v_point(k).cloned().ok_or_else(|| {
                        "TopOpeBRep_FacesFiller::KeepRLine : no VPoint".to_string()
                    })?;
                    vpl = vp;
                    let samept = super::process_gr::equal_p_on_r(&vpf, &vpl)?;
                    if !samept {
                        break;
                    }
                    i = first_kept(&l, k + 1, check_keep);
                }
            }
        }

        if !middle {
            // VPI.Init(L, checkkeep); if (More) Next();
            let mut i = match first_kept(&l, 1, check_keep) {
                Some(k) => first_kept(&l, k + 1, check_keep),
                None => None,
            };
            while let Some(k) = i {
                let vp = l
                    .v_point(k)
                    .cloned()
                    .ok_or_else(|| "TopOpeBRep_FacesFiller::KeepRLine : no VPoint".to_string())?;
                vpl = vp;
                samevp = super::process_gr::equal_pon_r(&vpf, &vpl, &l)?;
                if !samevp {
                    break;
                }
                i = first_kept(&l, k + 1, check_keep);
            }
        }

        if !samevp {
            // xpu151098 : cto 904 C8
            let samept = super::process_gr::equal_p_on_r(&vpf, &vpl)?;
            if samept {
                match closed_e(&el) {
                    Some(vclo) => {
                        let tolvclo = tolerance_of(&vclo)?;
                        let ptclo = pnt_of(&vclo)?;
                        let tolf = vpf.tolerance()?;
                        let ptf = vpf.value()?;
                        let d = ptf.distance(&ptclo);
                        let sameclo = d < tolvclo.max(tolf);
                        if !sameclo {
                            return Ok(false);
                        }
                    }
                    None => return Ok(false),
                }
            }
        }

        let mut out = false;
        if samevp {
            let isper = BRepTool::edge_curve(&Edge(el.clone()))
                .map(|c| c.is_periodic())
                .unwrap_or(false);
            match l.vp_bounds() {
                Some((f, ll, n)) => {
                    if isper && n == 2 {
                        let ioo = if isedge1 { 2 } else { 1 };
                        let vpf1 = l.v_point(f).cloned().ok_or_else(|| {
                            "TopOpeBRep_FacesFiller::KeepRLine : no VPoint".to_string()
                        })?;
                        let vpl1 = l.v_point(ll).cloned().ok_or_else(|| {
                            "TopOpeBRep_FacesFiller::KeepRLine : no VPoint".to_string()
                        })?;
                        let sf = vpf1.state(ioo)?;
                        let sl = vpl1.state(ioo)?;
                        let bf = sf == State::In || sf == State::On;
                        let bl = sl == State::In || sl == State::On;
                        // xpu120898 : when projection fails we get unknown status
                        let bfl = if sf == State::Unknown || sl == State::Unknown {
                            bf || bl
                        } else {
                            bf && bl
                        };
                        out = !bfl;
                    } else {
                        out = true;
                    }
                }
                // cxx:236-258: the `else` arm when the bounds are not `isper &&
                // n == 2`.
                None => out = true,
            }
        }
        if out {
            return Ok(false);
        }

        let f1 = self
            .f1
            .as_ref()
            .ok_or_else(|| "TopOpeBRep_FacesFiller::KeepRLine : no F1".to_string())?
            .0
            .clone();
        let f2 = self
            .f2
            .as_ref()
            .ok_or_else(|| "TopOpeBRep_FacesFiller::KeepRLine : no F2".to_string())?
            .0
            .clone();
        let psc = self
            .pshape_classifier
            .clone()
            .ok_or_else(|| "TopOpeBRep_FacesFiller::KeepRLine : no classifier".to_string())?;
        let mut psc = psc.borrow_mut();
        let st = super::process_gr::st_bip_vp_on_f(&mut psc, &f1, &f2, &vpf, &vpl, &l, isedge1)?;
        Ok(st == State::In || st == State::On) // REST1
    }

    /// `ProcessSectionEdges()` (`ProcessSectionEdges.cxx:334-468`).
    pub fn process_section_edges(
        &mut self,
        facint: &mut FacesIntersector,
        hds: &mut HDataStructure,
    ) -> Result<(), String> {
        // recuperation des aretes d'intersection LES
        let les = self.get_esl(facint)?;

        // add LES edges as section edges in the DS.
        {
            let ds = hds.change_ds();
            for e in &les {
                if BRepTool::is_degenerated(&Edge(e.clone())) {
                    continue; // xpu290698
                }
                ds.add_section_edge(e);
                // cxx:354-355: `myDS->Shape(E); myDS->AncestorRank(E);` — pure
                // lookups whose result is discarded.
                let _ = ds.shape_index(e, true);
                let _ = ds.ancestor_rank_of(e);
            }
        }

        // LOI = liste des rank (1 ou 2) des aretes de section (liste LES)
        let mut loi: Vec<i32> = Vec::new();
        for eles in &les {
            let mut is1 = false;
            let mut is2 = false;
            for i in 1..=facint.nb_lines() {
                let l = facint
                    .lines()
                    .and_then(|h| h.value(i))
                    .ok_or_else(|| {
                        "TopOpeBRep_FacesFiller::ProcessSectionEdges : no line".to_string()
                    })?;
                if l.type_line_curve() != TypeLineCurve::Restriction {
                    continue;
                }
                if let Some(arc) = l.arc() {
                    if crate::topo_tools::is_equal(arc, eles) {
                        is1 = l.arc_is_edge(1);
                        is2 = l.arc_is_edge(2);
                        break;
                    }
                }
            }
            if is1 {
                loi.push(1);
            } else if is2 {
                loi.push(2);
            }
        }

        // ajout des aretes de section dans la DS de shape, connaissant leur rank
        {
            let ds = hds.change_ds();
            for (e, r) in les.iter().zip(loi.iter()) {
                ds.add_shape_with_ancestor(e, *r);
            }
        }

        // mapELE(arete(1)) -> {arete(2)} ; mapELE(arete(2)) -> {arete(1)}
        let mut map_ele: Vec<(TopoShape, Vec<TopoShape>)> = Vec::new();
        {
            let ds = hds.change_ds();
            for e1 in &les {
                let ie1 = ds.shape_index(e1, true);
                let re1 = ds.ancestor_rank(ie1);
                if re1 != 1 {
                    continue;
                }
                map_ele.push((e1.clone(), Vec::new()));

                for e2 in &les {
                    let ie2 = ds.shape_index(e2, true);
                    let re2 = ds.ancestor_rank(ie2);
                    if re2 == 0 || ie1 == ie2 || re2 == re1 {
                        continue;
                    }
                    if fun_brep_sdm_re(e1, e2)? {
                        map_ele
                            .last_mut()
                            .expect("just pushed")
                            .1
                            .push(e2.clone());
                    }
                }
            }
        }

        for (e1, list) in &map_ele {
            let ds = hds.change_ds();
            let ie1 = ds.shape_index(e1, true);
            let re1 = ds.ancestor_rank(ie1);
            let face1 = facint.face(re1)?.clone();
            let is_closing1 =
                BRepTool::is_closed_edge_face(&Edge(e1.clone()), &face1);
            for e2 in list {
                let ie2 = ds.shape_index(e2, true);
                let re2 = ds.ancestor_rank(ie2);
                let face2 = facint.face(re2)?;
                let is_closing2 = BRepTool::is_closed_edge_face(&Edge(e2.clone()), face2);
                let ref_first = is_closing1 || !is_closing2;
                ds.fill_shapes_same_domain_with_config(
                    e1,
                    e2,
                    Config::UnshGeometry,
                    Config::UnshGeometry,
                    ref_first,
                )?;
            }
        }
        Ok(())
    }
}
