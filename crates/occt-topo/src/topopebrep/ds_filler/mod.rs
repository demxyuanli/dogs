//! `TopOpeBRep_DSFiller` — fills a `TopOpeBRepDS_HDataStructure` with the
//! results of a shape/shape intersection.
//!
//! Source: `TKBool/TopOpeBRep/TopOpeBRep_DSFiller.hxx:44-160`,
//! `TopOpeBRep_DSFiller.cxx:135-1215` (OCCT 8.0.0).
//!
//! The class owns one intersector of each family (`3d`, `2d`) and one filler
//! per intersecting pair kind (face/face, edge/edge, face/edge). The insert
//! pass ([`insert`]) walks the intersector and dispatches every pair to the
//! matching filler; the complete pass ([`complete`]) deduces the parent
//! shells and wires from the recorded interferences.
//!
//! The two package-level helpers the insert pass ends with —
//! `BREP_sortonparameter` (`TopOpeBRep_sort.cxx:90-104`) and
//! `BREP_mergePDS` (`TopOpeBRep_mergePDS.cxx:61-196`) — plus
//! `BREP_correctgbound` / `BREP_UnfillSameDomain` / `FUN_shareNOG`
//! (`TopOpeBRep_DSFiller.cxx:174-427`) live in [`ds_tools`].
//!
//! UNPORTED, called out where the call sites are (source file:line):
//! - `FBOX_Prepare()` (`TopOpeBRepTool_box.cxx:22`) — the global
//!   `TopOpeBRepTool_HBoxTool` is a process-wide singleton in OCCT; the port's
//!   box tool is per-intersector;
//! - `FC2D_Prepare(S1, S2)` (`TopOpeBRepTool_2d.cxx:96`) — sets the global
//!   `TopOpeBRepTool_CurveTool` reference surface;
//! - `FDSSDM_prepare(HDS)` (`TopOpeBRepDS_samdom.cxx:50`) — the same-domain
//!   edge map of `TopOpeBRepDS_samdom`.

mod complete;
mod ds_tools;
mod insert;

pub use ds_tools::*;

use crate::abs::ShapeType;
use crate::shape::TopoShape;
use crate::topexp::Explorer;
use crate::topo_tools::is_equal;
use crate::topopebrep::ds::connex::Connexity;
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::edges_filler::EdgesFiller;
use crate::topopebrep::face_edge_filler::FaceEdgeFiller;
use crate::topopebrep::faces_filler::FacesFiller;
use crate::topopebrep::shape_intersector::ShapeIntersector;
use crate::topopebrep::shape_intersector2d::ShapeIntersector2d;
use crate::topopebrep::tool::ShapeClassifier;

/// `TopOpeBRep_DSFiller` (`hxx:44-160`).
pub struct DsFiller {
    /// `myShapeIntersector` (`hxx:154`).
    shape_intersector: ShapeIntersector,
    /// `myShapeIntersector2d` (`hxx:155`).
    shape_intersector2d: ShapeIntersector2d,
    /// `myFacesFiller` (`hxx:156`).
    faces_filler: FacesFiller,
    /// `myEdgesFiller` (`hxx:157`).
    edges_filler: EdgesFiller,
    /// `myFaceEdgeFiller` (`hxx:158`).
    face_edge_filler: FaceEdgeFiller,
    /// `myPShapeClassifier` (`hxx:159`). OCCT holds a raw pointer shared with
    /// `TopOpeBRep_FacesFiller`; the port owns one `Rc<RefCell<..>>` and hands
    /// a clone to the faces filler on every insert (`cxx:451-454`,
    /// `cxx:911-914`).
    p_shape_classifier: Option<std::rc::Rc<std::cell::RefCell<ShapeClassifier>>>,
    /// `FDSCNX_Prepare(S1, S2, HDS)`'s global tables (`TopOpeBRepDS_connex.cxx:22-36`).
    /// OCCT keeps them in `thread_local` globals; the port keeps the prepared
    /// tables on the filler so `GapFiller` can query them (`gxx:127-220`).
    connex: Connexity,
}

impl Default for DsFiller {
    fn default() -> Self {
        Self::new()
    }
}

impl DsFiller {
    /// `TopOpeBRep_DSFiller()` (`cxx:137-141`): `myPShapeClassifier` null.
    pub fn new() -> Self {
        Self {
            shape_intersector: ShapeIntersector::new(),
            shape_intersector2d: ShapeIntersector2d::new(),
            faces_filler: FacesFiller::new(),
            edges_filler: EdgesFiller::new(),
            face_edge_filler: FaceEdgeFiller::new(),
            p_shape_classifier: None,
            connex: Connexity::new(),
        }
    }

    /// `PShapeClassifier()` (`cxx:155-158`).
    pub fn p_shape_classifier(&self) -> Option<std::rc::Rc<std::cell::RefCell<ShapeClassifier>>> {
        self.p_shape_classifier.clone()
    }

    /// The `cxx:451-454` / `cxx:911-914` lazy creation followed by
    /// `SetPShapeClassifier(myPShapeClassifier)`.
    fn ensure_p_shape_classifier(&mut self) -> std::rc::Rc<std::cell::RefCell<ShapeClassifier>> {
        if self.p_shape_classifier.is_none() {
            self.p_shape_classifier = Some(std::rc::Rc::new(std::cell::RefCell::new(
                ShapeClassifier::new(),
            )));
        }
        self.p_shape_classifier
            .clone()
            .expect("p_shape_classifier was just created")
    }

    /// `ChangeShapeIntersector()` (`cxx:1182-1185`).
    pub fn change_shape_intersector(&mut self) -> &mut ShapeIntersector {
        &mut self.shape_intersector
    }

    /// `ChangeShapeIntersector2d()`.
    pub fn change_shape_intersector2d(&mut self) -> &mut ShapeIntersector2d {
        &mut self.shape_intersector2d
    }

    /// `ChangeFacesFiller()`.
    pub fn change_faces_filler(&mut self) -> &mut FacesFiller {
        &mut self.faces_filler
    }

    /// `ChangeEdgesFiller()`.
    pub fn change_edges_filler(&mut self) -> &mut EdgesFiller {
        &mut self.edges_filler
    }

    /// `ChangeFaceEdgeFiller()`.
    pub fn change_face_edge_filler(&mut self) -> &mut FaceEdgeFiller {
        &mut self.face_edge_filler
    }

    /// `CheckInsert(S1, S2)` (`cxx:1136-1139`): `!S1.IsEqual(S2)`.
    fn check_insert(&self, s1: &TopoShape, s2: &TopoShape) -> bool {
        !is_equal(s1, s2)
    }

    /// `IsMadeOf1d(S)` (`cxx:1026-1066`): true for a WIRE or an EDGE, and for a
    /// non-empty COMPOUND all of whose children are `IsMadeOf1d`.
    pub fn is_made_of_1d(&self, s: &TopoShape) -> bool {
        match s.shape_type() {
            ShapeType::Wire | ShapeType::Edge => true,
            ShapeType::Compound => {
                // `TopoDS_Iterator it(aS)` (`cxx:1040`): direct children only.
                let children = crate::bop_occt_util::iter_children(s);
                if children.is_empty() {
                    false
                } else {
                    let mut res = true;
                    for child in &children {
                        if !self.is_made_of_1d(child) {
                            res = false;
                            break;
                        }
                    }
                    res
                }
            }
            // `cxx:1061-1064`: != COMPOUND, WIRE, EDGE
            _ => false,
        }
    }

    /// `IsContext1d(S)` (`cxx:1072-1079`).
    pub fn is_context_1d(&self, s: &TopoShape) -> bool {
        self.is_made_of_1d(s)
    }

    /// `ClearShapeSameDomain(S1, S2, HDS)` (`cxx:1145-1179`).
    fn clear_shape_same_domain(
        &self,
        s1: &TopoShape,
        s2: &TopoShape,
        hds: &mut HDataStructure,
    ) -> bool {
        if !self.check_insert(s1, s2) {
            return false;
        }
        let mut ex1 = Explorer::new(s1, ShapeType::Face);
        while ex1.more() {
            let shape1 = ex1.current().clone();
            // `cxx:1156-1164`: the C++ test is `!HasShape(Shape1)`, inverted
            // with respect to the second loop (`cxx:1167-1175`).
            if !hds.has_shape(&shape1, true) && hds.has_same_domain(&shape1, true) {
                if let Some(list) = hds.change_ds().change_shape_same_domain(&shape1) {
                    list.clear();
                }
            }
            ex1.next();
        }
        let mut ex2 = Explorer::new(s2, ShapeType::Face);
        while ex2.more() {
            let shape2 = ex2.current().clone();
            if hds.has_shape(&shape2, true) && hds.has_same_domain(&shape2, true) {
                if let Some(list) = hds.change_ds().change_shape_same_domain(&shape2) {
                    list.clear();
                }
            }
            ex2.next();
        }
        true
    }

    /// `FDSCNX_Prepare(aS1, aS2, HDS)` (`TopOpeBRepDS_connex.cxx:127-220`):
    /// rebuilds the edge/connexity tables. `S1`/`S2` are unused in OCCT.
    fn fdscnx_prepare(&mut self, _s1: &TopoShape, _s2: &TopoShape, hds: &mut HDataStructure) {
        self.connex.prepare(hds);
    }

    /// `FDSSDM_prepare(HDS)` (`TopOpeBRepDS_samdom.cxx:50`): UNPORTED — no
    /// same-domain edge map in the port yet.
    fn fdssdm_prepare(&self, _hds: &mut HDataStructure) {}
}
