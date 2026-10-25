//! `TopOpeBRep_ShapeIntersector2d` — intersects two shapes, keeping only the
//! face/face SameDomain couples and the edge/edge couples inside them.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRep/TopOpeBRep_ShapeIntersector2d.hxx:47-120`,
//! `.cxx:28-411`.
//!
//! Differences against [`ShapeIntersector`](crate::topopebrep::ShapeIntersector)
//! that the `.cxx` keeps:
//! - there is no FE / EF / plain-EE family: only `FF` and `EEFF`;
//! - `FindFFIntersection` (`cxx:242-271`) does **not** run
//!   `FacesIntersector::Perform`; it only boxes the couple and sets
//!   `myFFDone = true`, i.e. every face couple is reported as an intersection
//!   and the faces are assumed SameDomain;
//! - `InitEEFFIntersection` (`cxx:306-331`) calls the two-argument
//!   `EdgesIntersector::SetFaces(F1, F2)` (`EdgesIntersector.cxx:98-106`), not
//!   the `SetFaces(F1, F2, B1, B2)` overload.
//!
//! Translation notes:
//! - `myHBoxTool = FBOX_GetHBoxTool()` and the two `SetHBoxTool` calls
//!   (`cxx:30-31`) are UNPORTED for the same reason as in `ShapeIntersector`:
//!   each `BoxSort` owns its `HBoxTool` and computes a missing box on demand.
//! - `DumpCurrent` (`cxx:130-151`) and `Index` (`cxx:156-181`) exist only under
//!   `#ifdef OCCT_DEBUG`; the release `Index` returns 0 and is the one noted.
//! - The `Standard_ProgramError` raises of `CurrentGeomShape` (`cxx:97`) and
//!   `Shape` (`cxx:404`) become `Option::None`.

use crate::abs::ShapeType;
use crate::shape::TopoShape;
use crate::topopebrep::edges_intersector::EdgesIntersector;
use crate::topopebrep::shape_scanner::ShapeScanner;
use crate::topopebrep::tool::shape_explorer::ShapeExplorer;

/// `class TopOpeBRep_ShapeIntersector2d` (`hxx:47`).
pub struct ShapeIntersector2d {
    /// `myShape1` (`hxx:100`).
    shape1: Option<TopoShape>,
    /// `myShape2` (`hxx:101`).
    shape2: Option<TopoShape>,
    // `myHBoxTool` (`hxx:102`): UNPORTED, see the module note.
    /// `myFaceExplorer` (`hxx:103`).
    face_explorer: ShapeExplorer,
    /// `myFaceScanner` (`hxx:104`).
    face_scanner: ShapeScanner,
    /// `myEdgeExplorer` (`hxx:105`).
    edge_explorer: ShapeExplorer,
    /// `myEdgeScanner` (`hxx:106`).
    edge_scanner: ShapeScanner,
    /// `myEEIntersector` (`hxx:107`).
    ee_intersector: EdgesIntersector,
    /// `myIntersectionDone` (`hxx:108`).
    intersection_done: bool,
    /// `myFFDone` (`hxx:109`).
    ff_done: bool,
    /// `myEEFFDone` (`hxx:110`).
    eeff_done: bool,
    /// `myFFInit` (`hxx:111`).
    ff_init: bool,
    /// `myEEFFInit` (`hxx:112`).
    eeff_init: bool,
}

impl ShapeIntersector2d {
    /// `TopOpeBRep_ShapeIntersector2d()` (`cxx:28-32`).
    pub fn new() -> Self {
        let mut si = Self {
            shape1: None,
            shape2: None,
            face_explorer: ShapeExplorer::new(),
            face_scanner: ShapeScanner::new(),
            edge_explorer: ShapeExplorer::new(),
            edge_scanner: ShapeScanner::new(),
            ee_intersector: EdgesIntersector::new(),
            intersection_done: false,
            ff_done: false,
            eeff_done: false,
            ff_init: false,
            eeff_init: false,
        };
        si.reset();
        si
    }

    /// `Reset()` (`cxx:38-48`).
    pub fn reset(&mut self) {
        self.intersection_done = false;
        self.ff_done = false;
        self.eeff_done = false;
        self.ff_init = false;
        self.eeff_init = false;
    }

    /// `Init(S1, S2)` (`cxx:52-57`). `myHBoxTool->Clear()` is UNPORTED.
    pub fn init(&mut self, s1: &TopoShape, s2: &TopoShape) {
        self.reset();
        self.shape1 = Some(s1.clone());
        self.shape2 = Some(s2.clone());
    }

    /// `SetIntersectionDone()` (`cxx:62-65`).
    fn set_intersection_done(&mut self) {
        self.intersection_done = self.ff_done || self.eeff_done;
    }

    /// `CurrentGeomShape(Index)` (`cxx:69-97`).
    ///
    /// `cxx:97` throws `Standard_ProgramError("CurrentGeomShape : no
    /// intersection 2d")`; the port answers `None`.
    pub fn current_geom_shape(&self, index: i32) -> Option<TopoShape> {
        if self.intersection_done {
            if self.ff_done {
                if index == 1 {
                    return self.face_scanner_current().ok();
                } else if index == 2 {
                    return self.face_explorer_current().ok();
                }
            } else if self.eeff_done {
                if index == 1 {
                    return self.edge_scanner_current().ok();
                } else if index == 2 {
                    return self.edge_explorer_current().ok();
                }
            }
        }
        None
    }

    /// `InitIntersection(S1, S2)` (`cxx:102-106`).
    pub fn init_intersection(&mut self, s1: &TopoShape, s2: &TopoShape) -> Result<(), String> {
        self.init(s1, s2);
        self.init_ff_intersection()
    }

    /// `MoreIntersection()` (`cxx:109-125`).
    pub fn more_intersection(&self) -> bool {
        self.intersection_done
    }

    /// `NextIntersection()` (`cxx:187-217`).
    pub fn next_intersection(&mut self) -> Result<(), String> {
        self.intersection_done = false;

        if self.ff_done {
            // Previous state of More(): two faces.
            self.ff_done = false;
            self.init_eeff_intersection()?;
            self.find_eeff_intersection()?;
            if !self.intersection_done {
                self.next_ff_couple()?;
                self.find_ff_intersection()?;
            }
        } else if self.eeff_done {
            self.next_eeff_couple()?;
            self.find_eeff_intersection()?;
            if !self.intersection_done {
                self.next_ff_couple()?;
                self.find_ff_intersection()?;
            }
        }

        if !self.intersection_done {
            self.init_ff_intersection()?;
        }
        Ok(())
    }

    /// `Shape(Index)` (`cxx:392-405`); the raise of `cxx:404` is `None`.
    pub fn shape(&self, index: i32) -> Option<&TopoShape> {
        if index == 1 {
            self.shape1.as_ref()
        } else if index == 2 {
            self.shape2.as_ref()
        } else {
            None
        }
    }

    /// `ChangeEdgesIntersector()` (`cxx:409-411`).
    pub fn change_edges_intersector(&mut self) -> &mut EdgesIntersector {
        &mut self.ee_intersector
    }

    // ========
    // FFFFFFFF
    // ========

    /// `InitFFIntersection()` (`cxx:226-239`).
    pub fn init_ff_intersection(&mut self) -> Result<(), String> {
        if !self.ff_init {
            let shape1 = self
                .shape1
                .clone()
                .ok_or_else(|| "InitFFIntersection : no shape 1".to_string())?;
            let shape2 = self
                .shape2
                .clone()
                .ok_or_else(|| "InitFFIntersection : no shape 2".to_string())?;
            self.face_scanner.clear();
            self.face_scanner
                .add_boxes_make_cob(&shape1, ShapeType::Face, ShapeType::Shape)?;
            self.face_explorer
                .init(&shape2, ShapeType::Face, ShapeType::Shape);
            self.face_scanner
                .init_from_explorer(&mut self.face_explorer)?;
            self.find_ff_intersection()?;
        }
        self.ff_init = true;
        Ok(())
    }

    /// `FindFFIntersection()` (`cxx:242-271`).
    ///
    /// The two `BS.Box(GS1)` / `BS.Box(GS2)` calls (`cxx:266-268`) only warm the
    /// box cache: their results are not read, and this overload of
    /// `FindFFIntersection` never runs `FacesIntersector::Perform`.
    pub fn find_ff_intersection(&mut self) -> Result<(), String> {
        self.ff_done = false;

        if self.more_ff_couple() {
            let gs1 = self.face_scanner_current()?;
            let gs2 = self.face_explorer_current()?;

            let bs = self.face_scanner.change_box_sort();
            let _ = bs.box_of(&gs1);
            let _ = bs.box_of(&gs2);

            self.ff_done = true;
        }

        self.set_intersection_done();
        Ok(())
    }

    /// `MoreFFCouple()` (`cxx:276-280`).
    pub fn more_ff_couple(&self) -> bool {
        self.face_scanner.more() && self.face_explorer.more()
    }

    /// `NextFFCouple()` (`cxx:285-298`).
    pub fn next_ff_couple(&mut self) -> Result<(), String> {
        self.face_scanner.next();
        let mut b1 = !self.face_scanner.more();
        let mut b2 = self.face_explorer.more();
        while b1 && b2 {
            self.face_explorer.next();
            self.face_scanner
                .init_from_explorer(&mut self.face_explorer)?;
            b1 = !self.face_scanner.more();
            b2 = self.face_explorer.more();
        }
        Ok(())
    }

    // ========
    // EEFFEEFF
    // ========

    /// `InitEEFFIntersection()` (`cxx:306-331`).
    pub fn init_eeff_intersection(&mut self) -> Result<(), String> {
        // Prepare exploration of the edges of the two current SameDomain faces.
        let face1 = self.face_scanner_current()?;
        let face2 = self.face_explorer_current()?;

        self.ee_intersector.set_faces(&face1, &face2);

        self.edge_scanner.clear();
        self.edge_scanner
            .add_boxes_make_cob(&face1, ShapeType::Edge, ShapeType::Shape)?;
        self.edge_explorer
            .init(&face2, ShapeType::Edge, ShapeType::Shape);
        self.edge_scanner
            .init_from_explorer(&mut self.edge_explorer)?;

        self.eeff_init = true;
        Ok(())
    }

    /// `FindEEFFIntersection()` (`cxx:336-367`).
    pub fn find_eeff_intersection(&mut self) -> Result<(), String> {
        self.eeff_done = false;
        while self.more_eeff_couple() {
            let gs1 = self.edge_scanner_current()?;
            let gs2 = self.edge_explorer_current()?;
            self.ee_intersector.perform(&gs1, &gs2);

            self.eeff_done = !self.ee_intersector.is_empty();
            if self.eeff_done {
                break;
            }
            self.next_eeff_couple()?;
        }
        self.set_intersection_done();
        Ok(())
    }

    /// `MoreEEFFCouple()` (`cxx:372-377`).
    pub fn more_eeff_couple(&self) -> bool {
        self.edge_scanner.more() && self.edge_explorer.more()
    }

    /// `NextEEFFCouple()` (`cxx:381-388`).
    pub fn next_eeff_couple(&mut self) -> Result<(), String> {
        self.edge_scanner.next();
        while !self.edge_scanner.more() && self.edge_explorer.more() {
            self.edge_explorer.next();
            self.edge_scanner
                .init_from_explorer(&mut self.edge_explorer)?;
        }
        Ok(())
    }

    // ========
    // helpers
    // ========

    /// `myFaceScanner.Current()` (`cxx:252`).
    fn face_scanner_current(&self) -> Result<TopoShape, String> {
        match self.face_scanner.current() {
            Some(Ok(s)) => Ok(s.clone()),
            Some(Err(e)) => Err(e),
            None => Err("ShapeIntersector2d : no current face".to_string()),
        }
    }

    /// `myFaceExplorer.Current()` (`cxx:253`).
    fn face_explorer_current(&self) -> Result<TopoShape, String> {
        self.face_explorer
            .current()
            .cloned()
            .ok_or_else(|| "ShapeIntersector2d : no current explored face".to_string())
    }

    /// `myEdgeScanner.Current()` (`cxx:340`).
    fn edge_scanner_current(&self) -> Result<TopoShape, String> {
        match self.edge_scanner.current() {
            Some(Ok(s)) => Ok(s.clone()),
            Some(Err(e)) => Err(e),
            None => Err("ShapeIntersector2d : no current edge".to_string()),
        }
    }

    /// `myEdgeExplorer.Current()` (`cxx:341`).
    fn edge_explorer_current(&self) -> Result<TopoShape, String> {
        self.edge_explorer
            .current()
            .cloned()
            .ok_or_else(|| "ShapeIntersector2d : no current explored edge".to_string())
    }
}

impl Default for ShapeIntersector2d {
    fn default() -> Self {
        Self::new()
    }
}
