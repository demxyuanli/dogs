//! `TopOpeBRep_ShapeIntersector` — intersects two shapes.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRep/TopOpeBRep_ShapeIntersector.hxx:57-180`,
//! `.cxx:72-1200`.
//!
//! Layout of the port:
//! - this file — state, `Init`, `InitIntersection`, `NextIntersection`, the
//!   accessors (`Shape`, `CurrentGeomShape`, `Change*Intersector`,
//!   `GetTolerances`);
//! - [`couples`] — the six couple families (`FF`, `EEFF`, `FE`, `EF`, `EE`)
//!   that drive the three sub-intersectors;
//! - [`rejected_faces`] — `RejectedFaces` together with the two file-local
//!   statics it needs (`OneShapeIsHalfSpace`, `GetNewSolid`).
//!
//! Translation notes:
//! - `myHBoxTool = FBOX_GetHBoxTool()` (`cxx:76-79`) shares one `HBoxTool`
//!   between the scanners, the `FacesIntersector` and `BRepTopAdaptor`; the
//!   port gives each `BoxSort` its own `HBoxTool` and computes a missing box
//!   on demand (`tool/box_sort.rs:269-283`), so the `SetHBoxTool` sharing is
//!   UNPORTED (it changes how many boxes are cached, not which are compared).
//! - `Index` (`cxx:500-546`) and `DumpCurrent` (`cxx:456-495`) exist only under
//!   `#ifdef OCCT_DEBUG`; the release definitions (`cxx:543-547`,
//!   `cxx:495`) are the ones ported here.
//! - The `Standard_Failure` raises of `CurrentGeomShape` (`cxx:180`) and
//!   `Shape` (`cxx:1047`) are `None` here.
//! - `Result<_, String>` carries the failures of the `BoxSort` /
//!   `ShapeScanner` helpers, whose OCCT counterparts cannot fail.

use crate::shape::{Face, TopoShape};
use crate::topopebrep::edges_intersector::EdgesIntersector;
use crate::topopebrep::face_edge_intersector::FaceEdgeIntersector;
use crate::topopebrep::faces_intersector::FacesIntersector;
use crate::topopebrep::shape_scanner::ShapeScanner;
use crate::topopebrep::tool::shape_explorer::ShapeExplorer;

mod couples;
mod rejected_faces;

/// `class TopOpeBRep_ShapeIntersector` (`hxx:57`).
pub struct ShapeIntersector {
    /// `myShape1` (`hxx:154`).
    shape1: Option<TopoShape>,
    /// `myShape2` (`hxx:155`).
    shape2: Option<TopoShape>,
    // `myHBoxTool` (`hxx:156`): UNPORTED, see the module note.
    /// `myFaceExplorer` (`hxx:157`).
    face_explorer: ShapeExplorer,
    /// `myFaceScanner` (`hxx:158`).
    face_scanner: ShapeScanner,
    /// `myFFIntersector` (`hxx:159`).
    ff_intersector: FacesIntersector,
    /// `myFFSameDomain` (`hxx:160`).
    ff_same_domain: bool,
    /// `myEdgeExplorer` (`hxx:161`).
    edge_explorer: ShapeExplorer,
    /// `myEdgeScanner` (`hxx:162`).
    edge_scanner: ShapeScanner,
    /// `myEEIntersector` (`hxx:163`).
    ee_intersector: EdgesIntersector,
    /// `myFEIntersector` (`hxx:165`).
    fe_intersector: FaceEdgeIntersector,
    /// `myEEFace1` (`hxx:166`).
    ee_face1: Option<Face>,
    /// `myEEFace2` (`hxx:167`).
    ee_face2: Option<Face>,
    /// `myIntersectionDone` (`hxx:168`).
    intersection_done: bool,
    /// `myTol1` (`hxx:169`).
    tol1: f64,
    /// `myTol2` (`hxx:170`).
    tol2: f64,
    /// `myFFDone` (`hxx:171`).
    ff_done: bool,
    /// `myEEFFDone` (`hxx:172`).
    eeff_done: bool,
    /// `myEFDone` (`hxx:173`).
    ef_done: bool,
    /// `myFEDone` (`hxx:174`).
    fe_done: bool,
    /// `myEEDone` (`hxx:175`).
    ee_done: bool,
    /// `myFFInit` (`hxx:176`).
    ff_init: bool,
    /// `myEEFFInit` (`hxx:177`).
    eeff_init: bool,
    /// `myEFInit` (`hxx:178`).
    ef_init: bool,
    /// `myFEInit` (`hxx:179`).
    fe_init: bool,
    /// `myEEInit` (`hxx:180`).
    ee_init: bool,
}

impl ShapeIntersector {
    /// `TopOpeBRep_ShapeIntersector()` (`cxx:72-80`).
    pub fn new() -> Self {
        let mut si = Self {
            shape1: None,
            shape2: None,
            face_explorer: ShapeExplorer::new(),
            face_scanner: ShapeScanner::new(),
            ff_intersector: FacesIntersector::new(),
            ff_same_domain: false,
            edge_explorer: ShapeExplorer::new(),
            edge_scanner: ShapeScanner::new(),
            ee_intersector: EdgesIntersector::new(),
            fe_intersector: FaceEdgeIntersector::new(),
            ee_face1: None,
            ee_face2: None,
            intersection_done: false,
            tol1: 0.0,
            tol2: 0.0,
            ff_done: false,
            eeff_done: false,
            ef_done: false,
            fe_done: false,
            ee_done: false,
            ff_init: false,
            eeff_init: false,
            ef_init: false,
            fe_init: false,
            ee_init: false,
        };
        si.reset();
        let (tol1, tol2) = si.ff_intersector.get_tolerances();
        si.tol1 = tol1;
        si.tol2 = tol2;
        si
    }

    /// `Reset()` (`cxx:86-99`).
    fn reset(&mut self) {
        self.intersection_done = false;

        self.ff_done = false;
        self.ff_same_domain = false;
        self.eeff_done = false;
        self.ef_done = false;
        self.fe_done = false;
        self.ee_done = false;

        self.ff_init = false;
        self.eeff_init = false;
        self.ef_init = false;
        self.fe_init = false;
        self.ee_init = false;
    }

    /// `Init(S1, S2)` (`cxx:105-109`).
    pub fn init(&mut self, s1: &TopoShape, s2: &TopoShape) {
        self.reset();
        self.shape1 = Some(s1.clone());
        self.shape2 = Some(s2.clone());
    }

    /// `SetIntersectionDone()` (`cxx:114-117`).
    fn set_intersection_done(&mut self) {
        self.intersection_done =
            self.ff_done || self.eeff_done || self.fe_done || self.ef_done || self.ee_done;
    }

    /// `InitIntersection(S1, S2)` (`cxx:373-394`).
    pub fn init_intersection(&mut self, s1: &TopoShape, s2: &TopoShape) -> Result<(), String> {
        self.init(s1, s2);

        self.init_ff_intersection()?;
        if self.more_ff_couple() {
            return Ok(());
        }

        self.init_fe_intersection()?;
        if self.more_fe_couple() {
            return Ok(());
        }

        self.init_ef_intersection()?;
        if self.more_ef_couple() {
            return Ok(());
        }

        Ok(())
    }

    /// `InitIntersection(S1, S2, F1, F2)` (`cxx:397-408`).
    pub fn init_intersection_with_faces(
        &mut self,
        s1: &TopoShape,
        s2: &TopoShape,
        f1: &Face,
        f2: &Face,
    ) -> Result<(), String> {
        self.init(s1, s2);

        self.ee_face1 = Some(f1.clone());
        self.ee_face2 = Some(f2.clone());

        self.init_ee_intersection()
    }

    /// `MoreIntersection()` (`cxx:412-450`).
    pub fn more_intersection(&self) -> bool {
        self.intersection_done
    }

    /// `NextIntersection()` (`cxx:550-622`).
    pub fn next_intersection(&mut self) -> Result<(), String> {
        self.intersection_done = false;

        if self.ff_same_domain {
            // Previous state of More(): two SameDomain faces.
            self.ff_done = false;
            self.ff_same_domain = false;
            self.init_eeff_intersection()?;
            self.find_eeff_intersection()?;
            if !self.intersection_done {
                self.next_ff_couple()?;
                self.find_ff_intersection()?;
            }
        } else if self.ff_done {
            self.next_ff_couple()?;
            self.find_ff_intersection()?;
        } else if self.eeff_done {
            self.next_eeff_couple()?;
            self.find_eeff_intersection()?;
            if !self.intersection_done {
                self.next_ff_couple()?;
                self.find_ff_intersection()?;
            }
        } else if self.fe_done {
            self.next_fe_couple()?;
            self.find_fe_intersection()?;
        } else if self.ef_done {
            self.next_ef_couple()?;
            self.find_ef_intersection()?;
        } else if self.ee_done {
            self.next_ee_couple()?;
            self.find_ee_intersection()?;
        }

        if !self.intersection_done {
            self.init_ff_intersection()?;
        }

        if !self.intersection_done {
            self.init_fe_intersection()?;
        }

        if !self.intersection_done {
            self.init_ef_intersection()?;
        }

        if !self.intersection_done && self.ee_face1.is_some() && self.ee_face2.is_some() {
            self.init_ee_intersection()?;
        }

        Ok(())
    }

    /// `CurrentGeomShape(Index)` (`cxx:119-181`).
    ///
    /// `None` is the port of the `Standard_Failure("CurrentGeomShape : no
    /// intersection")` raise at `cxx:180`.
    pub fn current_geom_shape(&self, index: i32) -> Option<TopoShape> {
        if self.intersection_done {
            if self.ff_done {
                if index == 1 {
                    return self.face_scanner.current().and_then(|r| r.ok()).cloned();
                } else if index == 2 {
                    return self.face_explorer.current().cloned();
                }
            } else if self.eeff_done {
                if index == 1 {
                    return self.edge_scanner.current().and_then(|r| r.ok()).cloned();
                } else if index == 2 {
                    return self.edge_explorer.current().cloned();
                }
            } else if self.fe_done {
                if index == 1 {
                    return self.face_scanner.current().and_then(|r| r.ok()).cloned();
                } else if index == 2 {
                    return self.edge_explorer.current().cloned();
                }
            } else if self.ef_done {
                if index == 1 {
                    return self.edge_scanner.current().and_then(|r| r.ok()).cloned();
                } else if index == 2 {
                    return self.face_explorer.current().cloned();
                }
            } else if self.ee_done {
                if index == 1 {
                    return self.edge_scanner.current().and_then(|r| r.ok()).cloned();
                } else if index == 2 {
                    return self.edge_explorer.current().cloned();
                }
            }
        }

        None
    }

    /// `Shape(Index)` (`cxx:1036-1048`).
    ///
    /// `None` is the port of the `Standard_Failure("ShapeIntersector : no
    /// shape")` raise at `cxx:1047`.
    pub fn shape(&self, index: i32) -> Option<&TopoShape> {
        if index == 1 {
            self.shape1.as_ref()
        } else if index == 2 {
            self.shape2.as_ref()
        } else {
            None
        }
    }

    /// `ChangeFacesIntersector()` (`cxx:1054-1057`).
    pub fn change_faces_intersector(&mut self) -> &mut FacesIntersector {
        &mut self.ff_intersector
    }

    /// `ChangeEdgesIntersector()` (`cxx:1061-1064`).
    pub fn change_edges_intersector(&mut self) -> &mut EdgesIntersector {
        &mut self.ee_intersector
    }

    /// `ChangeFaceEdgeIntersector()` (`cxx:1068-1071`).
    pub fn change_face_edge_intersector(&mut self) -> &mut FaceEdgeIntersector {
        &mut self.fe_intersector
    }

    /// `GetTolerances(double&, double&)` (`cxx:1075-1079`).
    pub fn get_tolerances(&self) -> (f64, f64) {
        (self.tol1, self.tol2)
    }

    /// `Index(K)` (`cxx:543-547`): the release definition, which is
    /// `return 0`. The `#ifdef OCCT_DEBUG` body (`cxx:500-541`) is UNPORTED.
    pub fn index(&self, _k: i32) -> i32 {
        0
    }

    /// `DumpCurrent(K)` (`cxx:495`): the release definition, an empty body.
    /// The `#ifdef OCCT_DEBUG` body (`cxx:456-492`) is UNPORTED.
    pub fn dump_current(&self, _k: i32) {}
}
