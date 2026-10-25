//! `TopOpeBRep_EdgesFiller` — fills the DS with the edge/edge intersection
//! data described by `TopOpeBRep_EdgesIntersector`.
//!
//! Source (OCCT 8.0.0, `ModelingAlgorithms/TKBool/TopOpeBRep`):
//! - `EdgesFiller.hxx:39-110` (members `myPDS` / `myPEI` / `myE1` / `myE2` /
//!   `myF1` / `myF2` / `myLI1` / `myLI2` / `myHDS`)
//! - `EdgesFiller.cxx:60-66` constructor, `:72-363` `Insert`, `:372-425`
//!   `SetShapeTransition`, `:430-440` `GetGeometry`, `:444-473` `MakeGeometry`,
//!   `:477-508` the two `Face` accessors, `:514-547` `StorePI`, `:551-577`
//!   `StoreVI`, `:581-590` `ToRecompute`, `:594-605` `StoreRecompute`,
//!   `:608-648` `RecomputeInterferences`
//!
//! Translation notes:
//! - `myPDS` / `myPEI` / `myHDS` (`hxx:84-90`) are raw pointers that OCCT keeps
//!   across calls; the port passes `&mut EdgesIntersector` / `&mut
//!   HDataStructure` to every method instead, like `FacesFiller`.
//! - `myE1` / `myE2` / `myF1` / `myF2` are default-constructed (null) shapes in
//!   OCCT (`IsNull()`); the port stores `Option<TopoShape>`.
//! - The `OCCT_DEBUG` blocks (`cxx:142-160`, `:183-195`, `:227-229`,
//!   `:258-265`) are compiled out and not translated.
//! - `RecomputeInterferences` (`cxx:608-648`) builds `TN` and passes it to
//!   `FDS_stateEwithF2d`, then discards it: the only reachable effect is the
//!   call itself. The port keeps the call and the discard.

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::shape::TopoShape;
use crate::topopebrep::ds::data_structure::DataStructure;
use crate::topopebrep::ds::export::{fds_parameter, fds_state_e_with_f2d};
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::interference::{InterferenceData, InterferenceHandle};
use crate::topopebrep::ds::interference_tool::{
    make_edge_interference, make_edge_vertex_interference,
};
use crate::topopebrep::ds::kind::{Config, Kind};
use crate::topopebrep::ds::point_iterator::PointIterator;
use crate::topopebrep::ds::tki::TKI;
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::ds::InterferenceList;
use crate::topopebrep::edges_intersector::EdgesIntersector;
use crate::topopebrep::interf::point2d::Point2d;
use crate::topopebrep::interf::point_geom_tool::PointGeomTool;

/// `class TopOpeBRep_EdgesFiller` (`EdgesFiller.hxx:40-110`).
#[derive(Default)]
pub struct EdgesFiller {
    /// `myE1` / `myE2` (`hxx:85-86`).
    e1: Option<TopoShape>,
    e2: Option<TopoShape>,
    /// `myF1` / `myF2` (`hxx:87-88`).
    f1: Option<TopoShape>,
    f2: Option<TopoShape>,
    /// `myLI1` / `myLI2` (`hxx:88-89`): the interferences `ToRecompute` flagged.
    li1: InterferenceList,
    li2: InterferenceList,
    /// `myPEI->SameDomain()` (`cxx:585`), read by `ToRecompute` after `Insert`
    /// set `myPEI`. OCCT re-reads the pointer; the port captures the flag.
    esd: bool,
}

impl EdgesFiller {
    /// `TopOpeBRep_EdgesFiller()` (`cxx:60-66`): both pointers null.
    pub fn new() -> Self {
        Self::default()
    }

    /// `Insert(E1, E2, EDGINT, HDS)` (`cxx:72-363`).
    pub fn insert(
        &mut self,
        e1: &TopoShape,
        e2: &TopoShape,
        ei: &mut EdgesIntersector,
        hds: &mut HDataStructure,
    ) -> Result<(), String> {
        // `cxx:77-83`
        self.e1 = Some(e1.clone());
        self.e2 = Some(e2.clone());
        self.li1.clear();
        self.li2.clear();

        // `cxx:85-88`
        let esd = ei.same_domain();
        self.esd = esd;
        if esd {
            hds.change_ds()
                .fill_shapes_same_domain(e1, e2, true)?;
        }

        // exit if no point.
        ei.init_point(true);
        if !ei.more_point() {
            return Ok(());
        }

        // `cxx:97-98`: add <E1,E2> in BDS.
        let e1_index = hds.change_ds().add_shape_with_ancestor(e1, 1);
        let e2_index = hds.change_ds().add_shape_with_ancestor(e2, 2);

        while ei.more_point() {
            // `cxx:111-113`
            let p2d = ei.point()?.clone();
            let mut par1 = p2d.parameter(1)?;
            let mut par2 = p2d.parameter(2)?;

            // `cxx:115-122`
            if let Some(f1) = self.f1.clone() {
                hds.change_ds().add_shape_with_ancestor(&f1, 1);
            }
            if let Some(f2) = self.f2.clone() {
                hds.change_ds().add_shape_with_ancestor(&f2, 2);
            }

            // `cxx:124-126`
            let mut t1 = p2d.transition(1)?.clone();
            let mut t2 = p2d.transition(2)?.clone();
            self.set_shape_transition(&p2d, ei, hds, &mut t1, &mut t2)?;

            // `cxx:128-140`
            let is_vertex1 = p2d.is_vertex(1)?;
            let v1 = if is_vertex1 {
                Some(p2d.vertex(1)?.clone())
            } else {
                None
            };
            let is_vertex2 = p2d.is_vertex(2)?;
            let v2 = if is_vertex2 {
                Some(p2d.vertex(2)?.clone())
            } else {
                None
            };
            let is_vertex = is_vertex1 || is_vertex2;

            // `cxx:162-176`: the faces intersector's output is assumed valid and
            // is used to correct the edges intersector's output.
            let mut found = false;
            let mut g = 0;
            let mut k = Kind::Point;
            {
                let list = hds.ds().shape_interferences(e1, true).to_vec();
                let mut it = PointIterator::new(list);
                let dsp = PointGeomTool::make_point_of_point2d(&p2d);
                if let Some((gi, ki)) = hds.get_geometry(&mut it, &dsp)? {
                    found = true;
                    g = gi;
                    k = ki;
                }
            }
            if !found {
                self.make_geometry(&p2d, hds.change_ds(), &mut g, &mut k)?;
            }

            // `cxx:178-180`
            let found_point = found && k == Kind::Point;
            let is_new_point = !found && k == Kind::Point;
            let is_new_vertex = !found && k == Kind::Vertex;

            // `cxx:181-196`: `faulty` is computed by OCCT and only reported
            // under `OCCT_DEBUG`; the release flow does not branch on it.
            let _faulty = (is_vertex && is_new_point) || (!is_vertex && is_new_vertex);

            // `cxx:198-279`
            if is_vertex && found_point {
                let ns = hds.ds().nb_shapes();
                for is in 1..=ns {
                    let Some(s) = hds.ds().shape(is, false) else {
                        continue;
                    };
                    if s.shape_type() != ShapeType::Edge {
                        continue;
                    }
                    let e = s;
                    let mut linew: InterferenceList = Vec::new();
                    // `cxx:211-273`: walk the list, replace every matching CPI by
                    // an EVI on the current vertex, and remove it. The scan is
                    // split in two passes so the DS is not borrowed twice.
                    let mut hits: Vec<(usize, Transition, i32, f64, Kind)> = Vec::new();
                    {
                        let li = hds.change_ds().change_shape_interferences(&e);
                        for (idx, i) in li.iter().enumerate() {
                            let r = i.read().expect("interference lock");
                            let ki = r.geometry_type();
                            let gi = r.geometry();
                            let is_cpi = matches!(r.data(), InterferenceData::CurvePoint { .. });
                            let cond_cpi = ki == Kind::Point && gi == g && is_cpi;
                            if cond_cpi {
                                hits.push((
                                    idx,
                                    r.transition().clone(),
                                    r.support(),
                                    r.parameter().unwrap_or(0.0),
                                    r.support_type(),
                                ));
                            }
                        }
                    }
                    if !hits.is_empty() {
                        // replace G,K of I by the current vertex
                        let mut gevi = 0;
                        if is_vertex1 {
                            gevi =
                                hds.change_ds()
                                    .add_shape_with_ancestor(v1.as_ref().unwrap(), 1);
                        } else if is_vertex2 {
                            gevi =
                                hds.change_ds()
                                    .add_shape_with_ancestor(v2.as_ref().unwrap(), 2);
                        }
                        let bevi = false;
                        let cevi = Config::UnshGeometry;
                        for (_, tevi, sevi, pevi, kevi) in &hits {
                            let evi =
                                make_edge_vertex_interference(tevi, *sevi, gevi, bevi, cevi, *pevi);
                            // `cxx:254-255`: `evi->SupportType(epi->SupportType())`.
                            evi.write().expect("interference lock").set_support_type(*kevi);
                            linew.push(evi);
                        }
                        // `li.Remove(it)`: the scan of an OCCT list removes the
                        // element the iterator sits on; removing in descending
                        // index order keeps the remaining indices valid.
                        let li = hds.change_ds().change_shape_interferences(&e);
                        let mut idxs: Vec<usize> = hits.iter().map(|h| h.0).collect();
                        idxs.sort_unstable();
                        for idx in idxs.into_iter().rev() {
                            li.remove(idx);
                        }
                    }
                    if !linew.is_empty() {
                        hds.store_interferences_on_index(&linew, is);
                    }
                }
            }

            // `cxx:281-296`
            if is_vertex1 {
                let vv1 = v1.as_ref().unwrap();
                par1 = BRepTool::parameter_of_vertex(e1, vv1);
            }
            if is_vertex2 {
                let vv2 = v2.as_ref().unwrap();
                par2 = BRepTool::parameter_of_vertex(e2, vv2);
            }

            // `cxx:298-301`
            if is_vertex1 && is_vertex2 {
                hds.change_ds().fill_shapes_same_domain(
                    v1.as_ref().unwrap(),
                    v2.as_ref().unwrap(),
                    true,
                )?;
            }

            if !is_vertex {
                // `cxx:306-316`
                let mut dsp_index = 0;
                let mut kkk = Kind::Point;
                let list = hds.ds().shape_interferences(e1, true).to_vec();
                let mut it = PointIterator::new(list);
                let dsp = PointGeomTool::make_point_of_point2d(&p2d);
                let epi_found = match hds.get_geometry(&mut it, &dsp)? {
                    Some((gi, ki)) => {
                        dsp_index = gi;
                        kkk = ki;
                        true
                    }
                    None => false,
                };
                if !epi_found {
                    self.make_geometry(&p2d, hds.change_ds(), &mut dsp_index, &mut kkk)?;
                }

                // `cxx:318`
                self.set_shape_transition(&p2d, ei, hds, &mut t1, &mut t2)?;

                // `cxx:319-330`
                if kkk == Kind::Point {
                    self.store_pi(&p2d, &t1, e2_index, dsp_index, par1, 1, hds)?;
                    self.store_pi(&p2d, &t2, e1_index, dsp_index, par2, 2, hds)?;
                } else if kkk == Kind::Vertex {
                    let v_index = dsp_index;
                    let bevi = false;
                    let cevi = Config::UnshGeometry;
                    self.store_vi(&p2d, &t1, e2_index, v_index, bevi, cevi, par1, 1, hds)?;
                    self.store_vi(&p2d, &t2, e1_index, v_index, bevi, cevi, par2, 2, hds)?;
                }
            } else {
                // `cxx:336-359`
                self.set_shape_transition(&p2d, ei, hds, &mut t1, &mut t2)?;

                if is_vertex1 {
                    let v = v1.as_ref().unwrap();
                    let v_index = hds.change_ds().add_shape_with_ancestor(v, 1);
                    let ssc = p2d.edges_config();
                    self.store_vi(&p2d, &t1, e2_index, v_index, true, ssc, par1, 1, hds)?;
                    self.store_vi(&p2d, &t2, e1_index, v_index, false, ssc, par2, 2, hds)?;
                }
                if is_vertex2 {
                    let v = v2.as_ref().unwrap();
                    let v_index = hds.change_ds().add_shape_with_ancestor(v, 2);
                    let ssc = p2d.edges_config();
                    self.store_vi(&p2d, &t1, e2_index, v_index, false, ssc, par1, 1, hds)?;
                    self.store_vi(&p2d, &t2, e1_index, v_index, true, ssc, par2, 2, hds)?;
                }
            }

            ei.next_point();
        }

        // `cxx:361-362`
        let (e1, e2) = (self.e1.clone(), self.e2.clone());
        let li1 = std::mem::take(&mut self.li1);
        if let Some(e1) = e1 {
            self.recompute_interferences(&e1, &li1, hds)?;
        }
        self.li1 = li1;
        let li2 = std::mem::take(&mut self.li2);
        if let Some(e2) = e2 {
            self.recompute_interferences(&e2, &li2, hds)?;
        }
        self.li2 = li2;
        Ok(())
    }

    /// `Face(ISI, F)` (`cxx:477-490`).
    pub fn set_face(&mut self, isi: i32, f: &TopoShape) -> Result<(), String> {
        match isi {
            1 => {
                self.f1 = Some(f.clone());
                Ok(())
            }
            2 => {
                self.f2 = Some(f.clone());
                Ok(())
            }
            _ => Err("Face(i,f) : ISI incorrect".to_string()),
        }
    }

    /// `Face(ISI)` (`cxx:495-508`).
    pub fn face(&self, isi: i32) -> Result<&TopoShape, String> {
        match isi {
            1 => self.f1.as_ref().ok_or_else(|| "Face(i) : ISI incorrect".to_string()),
            2 => self.f2.as_ref().ok_or_else(|| "Face(i) : ISI incorrect".to_string()),
            _ => Err("Face(i) : ISI incorrect".to_string()),
        }
    }

    /// `SetShapeTransition(P2D, T1, T2)` (`cxx:372-425`).
    fn set_shape_transition(
        &self,
        p2d: &Point2d,
        ei: &EdgesIntersector,
        hds: &mut HDataStructure,
        t1: &mut Transition,
        t2: &mut Transition,
    ) -> Result<(), String> {
        let pointofsegment = p2d.is_point_of_segment();
        let esd = ei.same_domain();

        if pointofsegment && esd {
            // `cxx:381-386`
            t1.set_shape_before(ShapeType::Edge);
            t1.set_shape_after(ShapeType::Edge);
            t2.set_shape_before(ShapeType::Edge);
            t2.set_shape_after(ShapeType::Edge);
            // `cxx:387-401`
            let mut ie1 = 0;
            if let Some(e1) = self.e1.as_ref() {
                ie1 = hds.change_ds().add_shape_with_ancestor(e1, 1);
            }
            let mut ie2 = 0;
            if let Some(e2) = self.e2.as_ref() {
                ie2 = hds.change_ds().add_shape_with_ancestor(e2, 2);
            }
            if self.e2.is_some() {
                t1.set_index(ie2);
            }
            if self.e1.is_some() {
                t2.set_index(ie1);
            }
        } else {
            // `cxx:403-424`
            t1.set_shape_before(ShapeType::Face);
            t1.set_shape_after(ShapeType::Face);
            t2.set_shape_before(ShapeType::Face);
            t2.set_shape_after(ShapeType::Face);
            let mut if1 = 0;
            if let Some(f1) = self.f1.as_ref() {
                if1 = hds.change_ds().add_shape_with_ancestor(f1, 1);
            }
            let mut if2 = 0;
            if let Some(f2) = self.f2.as_ref() {
                if2 = hds.change_ds().add_shape_with_ancestor(f2, 2);
            }
            if self.f1.is_some() {
                t2.set_index(if1);
            }
            if self.f2.is_some() {
                t1.set_index(if2);
            }
        }
        Ok(())
    }

    /// `MakeGeometry(P2D, G, K)` (`cxx:444-473`): always returns `true` in
    /// OCCT; the `Result` only carries the DS accesses.
    fn make_geometry(
        &self,
        p2d: &Point2d,
        ds: &mut DataStructure,
        g: &mut i32,
        k: &mut Kind,
    ) -> Result<bool, String> {
        let is_vertex1 = p2d.is_vertex(1)?;
        let is_vertex2 = p2d.is_vertex(2)?;
        if is_vertex1 && is_vertex2 {
            let g1 = ds.add_shape_with_ancestor(p2d.vertex(1)?, 1);
            ds.add_shape_with_ancestor(p2d.vertex(2)?, 2);
            *g = g1;
            *k = Kind::Vertex;
        } else if is_vertex1 {
            *g = ds.add_shape_with_ancestor(p2d.vertex(1)?, 1);
            *k = Kind::Vertex;
        } else if is_vertex2 {
            *g = ds.add_shape_with_ancestor(p2d.vertex(2)?, 2);
            *k = Kind::Vertex;
        } else {
            *g = ds.add_point(&PointGeomTool::make_point_of_point2d(p2d));
            *k = Kind::Point;
        }
        Ok(true)
    }

    /// `StorePI(P2D, T, SI, GI, param, IEmother)` (`cxx:514-547`).
    #[allow(clippy::too_many_arguments)]
    fn store_pi(
        &mut self,
        p2d: &Point2d,
        t: &Transition,
        si: i32,
        gi: i32,
        param: f64,
        ie_mother: i32,
        hds: &mut HDataStructure,
    ) -> Result<InterferenceHandle, String> {
        let i = make_edge_interference(t, Kind::Edge, si, Kind::Point, gi, param);
        let mother = match ie_mother {
            1 => self.e1.clone(),
            2 => self.e2.clone(),
            _ => None,
        };
        match mother {
            Some(e) => hds.store_interference_on_shape(&i, &e)?,
            None => {
                return Err("TopOpeBRep_EdgesFiller::StorePI : Emother".to_string());
            }
        }
        if self.to_recompute(p2d, ie_mother) {
            self.store_recompute(&i, ie_mother);
        }
        Ok(i)
    }

    /// `StoreVI(P2D, T, EI, VI, VisB, C, param, IEmother)` (`cxx:551-577`).
    #[allow(clippy::too_many_arguments)]
    fn store_vi(
        &mut self,
        p2d: &Point2d,
        t: &Transition,
        ei: i32,
        vi: i32,
        vis_b: bool,
        c: Config,
        param: f64,
        ie_mother: i32,
        hds: &mut HDataStructure,
    ) -> Result<InterferenceHandle, String> {
        let i = make_edge_vertex_interference(t, ei, vi, vis_b, c, param);
        let mother = match ie_mother {
            1 => self.e1.clone(),
            2 => self.e2.clone(),
            _ => None,
        };
        match mother {
            Some(e) => hds.store_interference_on_shape(&i, &e)?,
            None => {
                return Err("TopOpeBRep_EdgesFiller::StoreVI : Emother".to_string());
            }
        }
        if self.to_recompute(p2d, ie_mother) {
            self.store_recompute(&i, ie_mother);
        }
        Ok(i)
    }

    /// `ToRecompute(P2D, I, IEmother)` (`cxx:581-590`): `pointofsegment &&
    /// !esd`.
    ///
    /// `esd` is `myPEI->SameDomain()`; the intersector is not a member of the
    /// filler in the port, so the caller passes it through
    /// [`Self::insert`]'s `esd` copy — which is what `store_pi` / `store_vi`
    /// use.
    fn to_recompute(&self, p2d: &Point2d, _ie_mother: i32) -> bool {
        p2d.is_point_of_segment() && !self.esd
    }

    /// `StoreRecompute(I, IEmother)` (`cxx:594-605`).
    fn store_recompute(&mut self, i: &InterferenceHandle, ie_mother: i32) {
        if ie_mother == 1 {
            self.li1.push(i.clone());
        } else if ie_mother == 2 {
            self.li2.push(i.clone());
        }
    }

    /// `RecomputeInterferences(E, LI)` (`cxx:608-648`).
    fn recompute_interferences(
        &self,
        e: &TopoShape,
        li: &InterferenceList,
        hds: &mut HDataStructure,
    ) -> Result<(), String> {
        if li.is_empty() {
            return Ok(());
        }
        let mut tki = TKI::new();
        tki.fill_on_geometry(li);
        tki.init();
        while tki.more() {
            let Some((k, g, _)) = tki.value() else {
                break;
            };
            // `cxx:626-632`: `loi = tki.ChangeValue(K, G)`, `nloi = loi.Extent()`.
            let (tu, pe) = {
                let Some((_, _, loi)) = tki.change_value() else {
                    break;
                };
                if loi.is_empty() {
                    tki.next();
                    continue;
                }
                let first = loi[0].clone();
                let r = first.read().expect("interference lock");
                (r.transition().clone(), fds_parameter(&r).unwrap_or(0.0))
            };

            let ifb = tu.index_before();
            // `cxx:637-638`: `myPDS->Shape(ifb)` is the one-argument accessor
            // (`DataStructure.cxx:1140-1147`), i.e. `FindKeep == false`.
            let fb = hds
                .ds()
                .shape(ifb, false)
                .ok_or_else(|| "TopOpeBRep_EdgesFiller::RecomputeInterferences : Shape".to_string())?;

            // `cxx:640-643`: TN copies shape and index of both sides.
            let mut tn = Transition::default();
            tn.set_shape_before(tu.shape_before());
            tn.set_index_before(tu.index_before());
            tn.set_shape_after(tu.shape_after());
            tn.set_index_after(tu.index_after());

            // `cxx:645`: the result is not used by OCCT either.
            let _ = fds_state_e_with_f2d(hds.ds(), e, pe, k, g, &fb, &mut tn)?;

            tki.next();
        }
        Ok(())
    }
}
