//! `TopOpeBRep_FaceEdgeFiller` — fills the DS with the face/edge intersection
//! data described by `TopOpeBRep_FaceEdgeIntersector`.
//!
//! Source (OCCT 8.0.0, `ModelingAlgorithms/TKBool/TopOpeBRep`):
//! - `FaceEdgeFiller.hxx:30-70` (members `myPDS` / `myHDS` are only reached
//!   through the `HDS` argument, so the port keeps no state at all)
//! - `FaceEdgeFiller.cxx:36` constructur, `:40-169` `Insert`, `:187-203`
//!   `ScanInterfList`, `:206-220` `GetGeometry(IT,FEINT,G,BDS)`, `:223-241`
//!   `GetGeometry(L,DSP,G,BDS)`, `:248-286` `StoreInterference`, `:289-316`
//!   `MakeGeometry`
//!
//! Translation notes:
//! - The `OCCT_DEBUG` blocks (`cxx:30-33` declaration, `:75-82` inside the
//!   point loop) and `TopOpeBRep_FaceEdgeFiller_DEB.cxx` are compiled out and
//!   not translated.
//! - `FEINT.UVPoint(pUV)` (`cxx:72-73`) only feeds the two blocks of `NYI`
//!   `FacePointInterference` / `FaceVertexInterference` code that OCCT keeps
//!   commented out (`cxx:107-116`, `:151-158`), so the port does not call it.
//! - `StoreInterference(I, LI, BDS)` (`cxx:248-286`) takes `LI` by reference;
//!   both call sites (`cxx:121`, `:144`, `:166`) pass `EIL2`, which is the live
//!   `ChangeShapeInterferences(EE)` list taken at `cxx:66`. The port therefore
//!   stores through [`HDataStructure::store_interference_on_shape`], which
//!   appends to that same list and then performs the same geometry append.

use crate::shape::TopoShape;
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::interference::InterferenceHandle;
use crate::topopebrep::ds::interference_tool::make_curve_interference;
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::ds::point_iterator::PointIterator;
use crate::topopebrep::ds::Point;
use crate::topopebrep::face_edge_intersector::FaceEdgeIntersector;
use crate::topopebrep::interf::point_geom_tool::PointGeomTool;

/// `class TopOpeBRep_FaceEdgeFiller` (`FaceEdgeFiller.hxx:30-70`).
///
/// OCCT stores `myPDS` / `myHDS` as pointers that `Insert` refreshes; every
/// method of this class takes `HDS` explicitly, so the port is stateless.
#[derive(Debug, Default)]
pub struct FaceEdgeFiller;

impl FaceEdgeFiller {
    /// `TopOpeBRep_FaceEdgeFiller()` (`cxx:36`): both pointers null.
    pub fn new() -> Self {
        Self
    }

    /// `Insert(F, E, FEINT, HDS)` (`cxx:40-169`).
    pub fn insert(
        &mut self,
        f: &TopoShape,
        e: &TopoShape,
        feint: &mut FaceEdgeIntersector,
        hds: &mut HDataStructure,
    ) -> Result<(), String> {
        // exit if no point.
        feint.init_point();
        if !feint.more_point() {
            return Ok(());
        }

        let ff_ori = f.orientation();

        // `cxx:60-62`: add <FF,EE> in BDS.
        let ff_index = hds.change_ds().add_shape_with_ancestor(f, 1);
        hds.change_ds().add_shape_with_ancestor(e, 2);

        while feint.more_point() {
            // `cxx:74`: `parE = FEINT.Parameter()`.
            let par_e = feint.parameter();

            // `cxx:82-83`: `TEE = FEINT.Transition(2, FFori); TEE.Index(FFindex)`.
            let mut tee = feint.transition(2, ff_ori)?;
            tee.set_index(ff_index);

            // `cxx:85-90`
            let v1 = feint.is_vertex(1);
            let v2 = feint.is_vertex(2);
            let is_vertex_f = v1.is_some();
            let is_vertex_e = v2.is_some();
            let is_vertex = is_vertex_f || is_vertex_e;

            if !is_vertex {
                // `cxx:99-101`: reuse the geometry of an equal point among the
                // interferences of <FF>, else create one.
                let list = hds.ds().shape_interferences(f, true).to_vec();
                let mut it = PointIterator::new(list);
                let dsp_index = match self.get_geometry(&mut it, feint, hds)? {
                    Some(g) => g,
                    None => self.make_geometry(feint, hds),
                };

                // `cxx:118-121`: EPI on EE, `MakeCurveInterference` (a
                // `CurvePointInterference` with a FACE/POINT geometry pair).
                let epi = make_curve_interference(
                    &tee,
                    Kind::Face,
                    ff_index,
                    Kind::Point,
                    dsp_index,
                    par_e,
                );
                self.store_interference(&epi, e, hds)?;
            } else if is_vertex_f {
                // `cxx:126-128`: a face vertex wins over an edge vertex.
                let v_index = hds
                    .change_ds()
                    .add_shape_with_ancestor(v1.as_ref().unwrap().shape(), 1);

                // `cxx:139-144`
                let epi = make_curve_interference(
                    &tee,
                    Kind::Face,
                    ff_index,
                    Kind::Vertex,
                    v_index,
                    par_e,
                );
                self.store_interference(&epi, e, hds)?;
            } else if is_vertex_e {
                // `cxx:147-149`
                let v_index = hds
                    .change_ds()
                    .add_shape_with_ancestor(v2.as_ref().unwrap().shape(), 2);

                // `cxx:160-166`
                let epi = make_curve_interference(
                    &tee,
                    Kind::Face,
                    ff_index,
                    Kind::Vertex,
                    v_index,
                    par_e,
                );
                self.store_interference(&epi, e, hds)?;
            }

            feint.next_point();
        }
        Ok(())
    }

    /// `ScanInterfList(IT, DSP, BDS)` (`cxx:187-203`): walk the list until an
    /// interference whose `Geometry()` names a DS point equal to `DSP`.
    ///
    /// OCCT reads `BDS.Point(DSPG)` unconditionally; a `Geometry()` that is not
    /// a DS point index raises inside `NCollection_DataMap::FindKey`, so the
    /// port returns an error there.
    fn scan_interf_list(
        &self,
        it: &mut PointIterator,
        dsp: &Point,
        hds: &HDataStructure,
    ) -> Result<bool, String> {
        while it.more() {
            let g = {
                let i = it.current();
                i.geometry()
            };
            // `const TopOpeBRepDS_Point& otherDSP = BDS.Point(DSPG);`
            let other = hds
                .ds()
                .point(g)
                .ok_or_else(|| "TopOpeBRepDS_DataStructure::Point".to_string())?;
            if PointGeomTool::is_equal(dsp, other) {
                return Ok(true);
            }
            it.next();
        }
        Ok(false)
    }

    /// `GetGeometry(IT, FEINT, G, BDS)` (`cxx:206-220`): `MakePoint(FEINT)` plus
    /// [`Self::scan_interf_list`], then `G = IT.Value()->Geometry()`.
    ///
    /// OCCT also has `GetGeometry(L, DSP, G, BDS)` (`cxx:223-241`), whose
    /// `else` arm adds `DSP` to the DS. It has no caller in the package; the
    /// port keeps the reachable half only and records the rest as `UNPORTED`.
    fn get_geometry(
        &self,
        it: &mut PointIterator,
        feint: &FaceEdgeIntersector,
        hds: &HDataStructure,
    ) -> Result<Option<i32>, String> {
        let dsp = PointGeomTool::make_point_of_face_edge_int(feint);
        if self.scan_interf_list(it, &dsp, hds)? {
            return Ok(Some(it.geometry()));
        }
        Ok(None)
    }

    /// `StoreInterference(I, LI, BDS)` (`cxx:248-286`).
    ///
    /// `LI` is `EIL2 == BDS.ChangeShapeInterferences(EE)` at every call site, so
    /// the port stores on the edge shape. The `switch (I->GeometryType())`
    /// body of OCCT is exactly [`HDataStructure::store_interference_on_shape`]'s
    /// geometry append (`SOLID/FACE/EDGE/VERTEX/POINT` are no-ops, `SURFACE`
    /// and `CURVE` append to the matching DS list).
    fn store_interference(
        &self,
        i: &InterferenceHandle,
        e: &TopoShape,
        hds: &mut HDataStructure,
    ) -> Result<(), String> {
        hds.store_interference_on_shape(i, e)
    }

    /// `MakeGeometry(FEINT, BDS)` (`cxx:289-316`).
    fn make_geometry(&self, feint: &FaceEdgeIntersector, hds: &mut HDataStructure) -> i32 {
        let v1 = feint.is_vertex(1);
        let v2 = feint.is_vertex(2);
        let is_vertex = v1.is_some() || v2.is_some();

        let mut g = 0;
        if is_vertex {
            if let Some(v) = v1 {
                g = hds.change_ds().add_shape_with_ancestor(v.shape(), 1);
            } else if let Some(v) = v2 {
                g = hds.change_ds().add_shape_with_ancestor(v.shape(), 2);
            }
        } else {
            g = hds
                .change_ds()
                .add_point(&PointGeomTool::make_point_of_face_edge_int(feint));
        }
        g
    }
}
