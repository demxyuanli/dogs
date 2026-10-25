//! The six couple families of `TopOpeBRep_ShapeIntersector`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRep/TopOpeBRep_ShapeIntersector.cxx:624-1032`.
//!
//! Every family has the same shape: `Init<X>Intersection` prepares the scanner
//! (over `myShape1`) and the explorer (over `myShape2`), then
//! `Find<X>Intersection` walks the couples until one produces a non-empty
//! intersection, `More<X>Couple`/`Next<X>Couple` advancing the pair.
//!
//! Translation notes:
//! - `AddBoxesMakeCOB(S, TS)` (`hxx:45-47`) defaults `TA` to `TopAbs_SHAPE`;
//!   the two-argument call sites below pass `ShapeType::Shape`.
//! - `Init(S, ToFind)` of `TopOpeBRepTool_ShapeExplorer` defaults `ToAvoid` to
//!   `TopAbs_SHAPE` for the same reason.
//! - OCCT reads `Current()` as a reference; the port clones it, because the
//!   scanners borrow `self` mutably right after.

use crate::abs::{Orientation, ShapeType};
use crate::shape::{Face, TopoShape};

use super::ShapeIntersector;

impl ShapeIntersector {
    /// `myFaceScanner.Current()` (`cxx:654`).
    ///
    /// OCCT assumes the iterator is positioned (`More()` is checked first); the
    /// port turns a missing / failed position into an `Err` instead of reading
    /// a null `TopoDS_Shape`.
    fn face_scanner_current(&self) -> Result<TopoShape, String> {
        match self.face_scanner.current() {
            Some(Ok(s)) => Ok(s.clone()),
            Some(Err(e)) => Err(e),
            None => Err("ShapeIntersector : no current face".to_string()),
        }
    }

    /// `myEdgeScanner.Current()` (`cxx:783`).
    fn edge_scanner_current(&self) -> Result<TopoShape, String> {
        match self.edge_scanner.current() {
            Some(Ok(s)) => Ok(s.clone()),
            Some(Err(e)) => Err(e),
            None => Err("ShapeIntersector : no current edge".to_string()),
        }
    }

    /// `myFaceExplorer.Current()` (`cxx:655`).
    fn face_explorer_current(&self) -> Result<TopoShape, String> {
        match self.face_explorer.current() {
            Some(s) => Ok(s.clone()),
            None => Err("ShapeIntersector : no current explored face".to_string()),
        }
    }

    /// `myEdgeExplorer.Current()` (`cxx:784`).
    fn edge_explorer_current(&self) -> Result<TopoShape, String> {
        match self.edge_explorer.current() {
            Some(s) => Ok(s.clone()),
            None => Err("ShapeIntersector : no current explored edge".to_string()),
        }
    }

    // ========
    // FFFFFFFF
    // ========

    /// `InitFFIntersection()` (`cxx:628-641`).
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
            self.face_scanner.init_from_explorer(&mut self.face_explorer)?;
            self.find_ff_intersection()?;
        }
        self.ff_init = true;
        Ok(())
    }

    /// `FindFFIntersection()` (`cxx:646-711`).
    pub fn find_ff_intersection(&mut self) -> Result<(), String> {
        self.ff_done = false;
        self.ff_same_domain = false;

        while self.more_ff_couple() {
            // The two candidate intersecting GeomShapes GS1, GS2.
            let gs1 = self.face_scanner_current()?;
            let gs2 = self.face_explorer_current()?;

            let (b1, b2) = {
                let bs = self.face_scanner.change_box_sort();
                (bs.box_of(&gs1), bs.box_of(&gs2))
            };
            let (Ok(b1), Ok(b2)) = (b1, b2) else {
                self.next_ff_couple()?;
                continue;
            };
            let (Some(f1), Some(f2)) = (Face::wrap(gs1), Face::wrap(gs2)) else {
                self.next_ff_couple()?;
                continue;
            };
            self.ff_intersector.perform(&f1, &f2, &b1, &b2);

            // `bool ok = myFFIntersector.IsDone();` (`cxx:675`).
            if !self.ff_intersector.is_done() {
                self.next_ff_couple()?;
                continue;
            }

            self.ff_same_domain = self.ff_intersector.same_domain()?;

            if self.ff_same_domain {
                self.ff_done = true;
                break;
            } else {
                self.ff_done = !self.ff_intersector.is_empty();

                // Update face/face intersection tolerances.
                if self.ff_done {
                    let (tol1, tol2) = self.ff_intersector.get_tolerances();
                    self.tol1 = self.tol1.max(tol1);
                    self.tol2 = self.tol2.max(tol2);
                }

                if self.ff_done {
                    break;
                }
            }

            self.next_ff_couple()?;
        }

        self.set_intersection_done();
        Ok(())
    }

    /// `MoreFFCouple()` (`cxx:716-719`).
    pub fn more_ff_couple(&self) -> bool {
        self.face_scanner.more() && self.face_explorer.more()
    }

    /// `NextFFCouple()` (`cxx:724-736`).
    pub fn next_ff_couple(&mut self) -> Result<(), String> {
        self.face_scanner.next();
        let mut b1 = !self.face_scanner.more();
        let mut b2 = self.face_explorer.more();
        while b1 && b2 {
            self.face_explorer.next();
            self.face_scanner.init_from_explorer(&mut self.face_explorer)?;
            b1 = !self.face_scanner.more();
            b2 = self.face_explorer.more();
        }
        Ok(())
    }

    // ========
    // EEFFEEFF
    // ========

    /// `InitEEFFIntersection()` (`cxx:747-773`).
    pub fn init_eeff_intersection(&mut self) -> Result<(), String> {
        // Prepare exploration of the edges of the two current SameDomain faces.
        let face1 = self.face_scanner_current()?;
        let face2 = self.face_explorer_current()?;

        let (b1, b2) = {
            let bs = self.face_scanner.change_box_sort();
            (bs.box_of(&face1), bs.box_of(&face2))
        };
        let (Ok(b1), Ok(b2)) = (b1, b2) else {
            return Err("InitEEFFIntersection : box".to_string());
        };
        self.ee_intersector
            .set_faces_with_boxes(&face1, &face2, &b1, &b2);

        self.edge_scanner.clear();
        self.edge_scanner
            .add_boxes_make_cob(&face1, ShapeType::Edge, ShapeType::Shape)?;
        self.edge_explorer
            .init(&face2, ShapeType::Edge, ShapeType::Shape);
        self.edge_scanner.init_from_explorer(&mut self.edge_explorer)?;

        self.eeff_init = true;
        Ok(())
    }

    /// `FindEEFFIntersection()` (`cxx:778-808`).
    pub fn find_eeff_intersection(&mut self) -> Result<(), String> {
        self.eeff_done = false;
        while self.more_eeff_couple() {
            let gs1 = self.edge_scanner_current()?;
            let gs2 = self.edge_explorer_current()?;
            self.ee_intersector.perform(&gs1, &gs2);

            self.eeff_done = !self.ee_intersector.is_empty();
            if self.eeff_done {
                break;
            } else {
                self.next_eeff_couple()?;
            }
        }
        self.set_intersection_done();
        Ok(())
    }

    /// `MoreEEFFCouple()` (`cxx:813-816`).
    pub fn more_eeff_couple(&self) -> bool {
        self.edge_scanner.more() && self.edge_explorer.more()
    }

    /// `NextEEFFCouple()` (`cxx:821-830`).
    pub fn next_eeff_couple(&mut self) -> Result<(), String> {
        self.edge_scanner.next();
        while !self.edge_scanner.more() && self.edge_explorer.more() {
            self.edge_explorer.next();
            self.edge_scanner.init_from_explorer(&mut self.edge_explorer)?;
        }
        Ok(())
    }

    // ========
    // FEFEFEFE
    // ========

    /// `InitFEIntersection()` (`cxx:837-852`).
    pub fn init_fe_intersection(&mut self) -> Result<(), String> {
        if !self.fe_init {
            let shape1 = self
                .shape1
                .clone()
                .ok_or_else(|| "InitFEIntersection : no shape 1".to_string())?;
            let shape2 = self
                .shape2
                .clone()
                .ok_or_else(|| "InitFEIntersection : no shape 2".to_string())?;
            self.face_scanner.clear();
            self.face_scanner
                .add_boxes_make_cob(&shape1, ShapeType::Face, ShapeType::Shape)?;
            // `NYI defaut de spec` (`cxx:846`): the edge explorer avoids faces.
            self.edge_explorer
                .init(&shape2, ShapeType::Edge, ShapeType::Face);
            self.face_scanner.init_from_explorer(&mut self.edge_explorer)?;
            self.find_fe_intersection()?;
        }
        self.fe_init = true;
        Ok(())
    }

    /// `FindFEIntersection()` (`cxx:857-874`).
    pub fn find_fe_intersection(&mut self) -> Result<(), String> {
        self.fe_done = false;
        while self.more_fe_couple() {
            let gs1 = self.face_scanner_current()?;
            let gs2 = self.edge_explorer_current()?;
            self.fe_intersector.perform(&gs1, &gs2);
            self.fe_done = !self.fe_intersector.is_empty();
            if self.fe_done {
                break;
            } else {
                self.next_fe_couple()?;
            }
        }
        self.set_intersection_done();
        Ok(())
    }

    /// `MoreFECouple()` (`cxx:879-882`).
    pub fn more_fe_couple(&self) -> bool {
        self.face_scanner.more() && self.edge_explorer.more()
    }

    /// `NextFECouple()` (`cxx:887-896`).
    pub fn next_fe_couple(&mut self) -> Result<(), String> {
        self.face_scanner.next();
        while !self.face_scanner.more() && self.edge_explorer.more() {
            self.edge_explorer.next();
            self.face_scanner.init_from_explorer(&mut self.edge_explorer)?;
        }
        Ok(())
    }

    // ========
    // EFEFEFEF
    // ========

    /// `InitEFIntersection()` (`cxx:904-917`).
    pub fn init_ef_intersection(&mut self) -> Result<(), String> {
        if !self.ef_init {
            let shape1 = self
                .shape1
                .clone()
                .ok_or_else(|| "InitEFIntersection : no shape 1".to_string())?;
            let shape2 = self
                .shape2
                .clone()
                .ok_or_else(|| "InitEFIntersection : no shape 2".to_string())?;
            self.edge_scanner.clear();
            // `NYI defaut de spec` (`cxx:910`): the edge scanner avoids faces.
            self.edge_scanner
                .add_boxes_make_cob(&shape1, ShapeType::Edge, ShapeType::Face)?;
            self.face_explorer
                .init(&shape2, ShapeType::Face, ShapeType::Shape);
            self.edge_scanner.init_from_explorer(&mut self.face_explorer)?;
            self.find_ef_intersection()?;
        }
        self.ef_init = true;
        Ok(())
    }

    /// `FindEFIntersection()` (`cxx:922-939`).
    pub fn find_ef_intersection(&mut self) -> Result<(), String> {
        self.ef_done = false;
        while self.more_ef_couple() {
            let gs1 = self.edge_scanner_current()?;
            let gs2 = self.face_explorer_current()?;
            // `myFEIntersector.Perform(GS2, GS1)` (`cxx:928`).
            self.fe_intersector.perform(&gs2, &gs1);
            self.ef_done = !self.fe_intersector.is_empty();
            if self.ef_done {
                break;
            } else {
                self.next_ef_couple()?;
            }
        }
        self.set_intersection_done();
        Ok(())
    }

    /// `MoreEFCouple()` (`cxx:944-947`).
    pub fn more_ef_couple(&self) -> bool {
        self.edge_scanner.more() && self.face_explorer.more()
    }

    /// `NextEFCouple()` (`cxx:952-960`).
    pub fn next_ef_couple(&mut self) -> Result<(), String> {
        self.edge_scanner.next();
        while !self.edge_scanner.more() && self.face_explorer.more() {
            self.face_explorer.next();
            self.edge_scanner.init_from_explorer(&mut self.face_explorer)?;
        }
        Ok(())
    }

    // ========
    // EEEEEEEE
    // ========

    /// `InitEEIntersection()` (`cxx:967-988`).
    pub fn init_ee_intersection(&mut self) -> Result<(), String> {
        if !self.ee_init {
            let face1 = self
                .ee_face1
                .as_ref()
                .ok_or_else(|| "InitEEIntersection : no face 1".to_string())?
                .0
                .clone()
                .oriented(Orientation::Forward);
            let face2 = self
                .ee_face2
                .as_ref()
                .ok_or_else(|| "InitEEIntersection : no face 2".to_string())?
                .0
                .clone()
                .oriented(Orientation::Forward);
            let (b1, b2) = {
                let bs = self.face_scanner.change_box_sort();
                (bs.box_of(&face1), bs.box_of(&face2))
            };
            let (Ok(b1), Ok(b2)) = (b1, b2) else {
                return Err("InitEEIntersection : box".to_string());
            };
            self.ee_intersector
                .set_faces_with_boxes(&face1, &face2, &b1, &b2);

            let shape1 = self
                .shape1
                .clone()
                .ok_or_else(|| "InitEEIntersection : no shape 1".to_string())?;
            let shape2 = self
                .shape2
                .clone()
                .ok_or_else(|| "InitEEIntersection : no shape 2".to_string())?;
            self.edge_scanner.clear();
            self.edge_scanner
                .add_boxes_make_cob(&shape1, ShapeType::Edge, ShapeType::Shape)?;
            self.edge_explorer
                .init(&shape2, ShapeType::Edge, ShapeType::Shape);
            self.edge_scanner.init_from_explorer(&mut self.edge_explorer)?;
            self.find_ee_intersection()?;
        }
        self.ee_init = true;
        Ok(())
    }

    /// `FindEEIntersection()` (`cxx:993-1012`).
    pub fn find_ee_intersection(&mut self) -> Result<(), String> {
        self.ee_done = false;
        while self.more_ee_couple() {
            let gs1 = self.edge_scanner_current()?;
            let gs2 = self.edge_explorer_current()?;
            self.ee_intersector.perform(&gs1, &gs2);
            self.ee_done = !self.ee_intersector.is_empty();
            if self.ee_done {
                break;
            } else {
                self.next_ee_couple()?;
            }
        }
        self.set_intersection_done();
        Ok(())
    }

    /// `MoreEECouple()` (`cxx:1017-1020`).
    pub fn more_ee_couple(&self) -> bool {
        self.edge_scanner.more() && self.edge_explorer.more()
    }

    /// `NextEECouple()` (`cxx:1025-1032`).
    pub fn next_ee_couple(&mut self) -> Result<(), String> {
        self.edge_scanner.next();
        while !self.edge_scanner.more() && self.edge_explorer.more() {
            self.edge_explorer.next();
            self.edge_scanner.init_from_explorer(&mut self.edge_explorer)?;
        }
        Ok(())
    }
}
