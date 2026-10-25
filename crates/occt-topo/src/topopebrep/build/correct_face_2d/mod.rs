//! `TopOpeBRepBuild_CorrectFace2d` — glue the 2d contours of a periodic face
//! back together.
//!
//! Source: `TKBool/TopOpeBRepBuild/TopOpeBRepBuild_CorrectFace2d.hxx:29-101`,
//! `.cxx:49-1128` (OCCT 8.0.0). Split by subject:
//!
//! - [`mod.rs`] — the state machine: `Perform` (`cxx:587-642`), `CheckFace`
//!   (`:226-336`), `MakeRightWire` (`:646-683`), the accessors and the
//!   `ErrorStatus` table (`:1120-1128`).
//! - [`copy_data`] — `BuildCopyData` (`:116-222`), `ConnectWire` (`:340-585`),
//!   `GetP2dFL` (`:687-702`), `MakeHeadList` (`:706-744`), `CheckList`
//!   (`:748-847`).
//! - [`wires_2d`] — `TranslateCurve2d` (`:851-862`), `OuterWire` (`:866-890`),
//!   `BndBoxWire` (`:894-907`), `MoveWire2d` (`:911-977`), `MoveWires2d`
//!   (`:982-1030`) and the two `UpdateEdge` overloads (`:1033-1117`).
//!
//! `myMap` is a `void*` in OCCT (`hxx:87`): the caller hands in an
//! `IndexedDataMap<TopoDS_Shape, TopoDS_Shape>` and the class writes the new
//! edges created for the *original* face into it. The port keeps the same
//! aliasing with an [`Rc<RefCell<..>>`](std::rc::Rc); when no map was set the
//! `pMap` arms are skipped (OCCT dereferences a null `void*` there, i.e. UB).

use std::cell::RefCell;
use std::rc::Rc;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::shape::{Face, TopoShape};
use crate::topo_tools::explorer_shapes;
use crate::topopebrep::build::tools2d;
use crate::topopebrep::ds::data_structure::{IndexedShapeData, IndexedShapeSet};

pub mod copy_data;
pub mod wires_2d;

/// `ErrorStatus` table (`cxx:1120-1128`).
///
/// 0 — all is done; 1 — nothing is done; 2 — the wire does not contain any
/// licit edge and the wire is not closed (documented only, no `return` sets
/// it); 3 — a wire is not connectable; 4 — the base surface is not periodic
/// (also set by `MakeRightWire` when the path walk fails, `cxx:667`); 5 — no
/// outer wire found; 6 — wrong connectability of a wire in 2d.
pub const ERROR_NOTHING_DONE: i32 = 1;
pub const ERROR_NOT_CONNECTABLE: i32 = 3;
pub const ERROR_NOT_PERIODIC: i32 = 4;
pub const ERROR_NO_OUTER_WIRE: i32 = 5;
pub const ERROR_WRONG_CONNECTABILITY: i32 = 6;

/// `class TopOpeBRepBuild_CorrectFace2d` (`hxx:38-100`).
pub struct CorrectFace2d {
    /// `myFace` (`hxx:78`).
    my_face: Face,
    /// `myCorrectedFace` (`hxx:79`).
    my_corrected_face: Face,
    /// `myIsDone` (`hxx:80`).
    my_is_done: bool,
    /// `myErrorStatus` (`hxx:81`).
    my_error_status: i32,
    /// `myFaceTolerance` (`hxx:82`).
    my_face_tolerance: f64,
    /// `myCurrentWire` (`hxx:83`): a null `TopoDS_Wire` until `CheckFace` /
    /// `Perform` set it.
    my_current_wire: TopoShape,
    /// `myOrderedWireList` (`hxx:84`).
    my_ordered_wire_list: Vec<TopoShape>,
    /// `myAvoidMap` (`hxx:85`).
    my_avoid_map: IndexedShapeSet,
    /// `myMap` (`hxx:87`), the `void*` aliasing the caller's
    /// `IndexedDataMap<TopoDS_Shape, TopoDS_Shape>`.
    my_map: Option<Rc<RefCell<IndexedShapeData<TopoShape>>>>,
    /// `myCopyFace` (`hxx:88`).
    my_copy_face: Face,
    /// `myCopyAvoidMap` (`hxx:89`).
    my_copy_avoid_map: IndexedShapeSet,
    /// `myEdMapInversed` (`hxx:90`): copy edge -> original edge.
    my_ed_map_inversed: IndexedShapeData<TopoShape>,
}

impl Default for CorrectFace2d {
    fn default() -> Self {
        Self::new()
    }
}

impl CorrectFace2d {
    /// `TopOpeBRepBuild_CorrectFace2d()` (`cxx:49-54`).
    ///
    /// OCCT leaves `myFaceTolerance` and `myMap` uninitialized here; the port
    /// uses `0.` and `None` (see the module note on `myMap`).
    pub fn new() -> Self {
        Self {
            my_face: Face::new(),
            my_corrected_face: Face::new(),
            my_is_done: false,
            my_error_status: ERROR_NOTHING_DONE,
            my_face_tolerance: 0.0,
            my_current_wire: TopoShape::new(ShapeType::Wire),
            my_ordered_wire_list: Vec::new(),
            my_avoid_map: IndexedShapeSet::new(),
            my_map: None,
            my_copy_face: Face::new(),
            my_copy_avoid_map: IndexedShapeSet::new(),
            my_ed_map_inversed: IndexedShapeData::new(),
        }
    }

    /// `TopOpeBRepBuild_CorrectFace2d(aFace, anAvoidMap, aMap)` (`cxx:56-68`).
    pub fn with_face(
        a_face: &Face,
        an_avoid_map: &IndexedShapeSet,
        a_map: Rc<RefCell<IndexedShapeData<TopoShape>>>,
    ) -> Self {
        Self {
            my_face: a_face.clone(),
            my_avoid_map: an_avoid_map.clone(),
            my_map: Some(a_map),
            ..Self::new()
        }
    }

    /// `Face()` (`cxx:86-90`).
    pub fn face(&self) -> &Face {
        &self.my_face
    }

    /// `IsDone()` (`cxx:92-96`).
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }

    /// `ErrorStatus()` (`cxx:100-104`).
    pub fn error_status(&self) -> i32 {
        self.my_error_status
    }

    /// `CorrectedFace()` (`cxx:106-111`).
    pub fn corrected_face(&self) -> &Face {
        &self.my_corrected_face
    }

    /// `SetMapOfTrans2dInfo(aMap)` (`cxx:70-74`).
    pub fn set_map_of_trans2d_info(&mut self, a_map: Rc<RefCell<IndexedShapeData<TopoShape>>>) {
        self.my_map = Some(a_map);
    }

    /// `MapOfTrans2dInfo()` (`cxx:78-84`).
    pub fn map_of_trans2d_info(&self) -> Option<Rc<RefCell<IndexedShapeData<TopoShape>>>> {
        self.my_map.clone()
    }

    /// `Perform()` (`cxx:587-642`).
    pub fn perform(&mut self) {
        // 0. (`cxx:589-597`)
        self.my_copy_avoid_map.clear();
        let face = self.my_face.clone();
        let avoid = self.my_avoid_map.clone();
        let mut copy_face = Face::new();
        let mut copy_avoid = IndexedShapeSet::new();
        self.build_copy_data(&face, &avoid, &mut copy_face, &mut copy_avoid, true);
        self.my_copy_face = copy_face;
        self.my_copy_avoid_map = copy_avoid;

        self.my_corrected_face = self.my_copy_face.clone();
        self.my_corrected_face
            .0
            .set_orientation(self.my_copy_face.0.orientation());
        self.my_face_tolerance = BRepTool::face_tolerance(&self.my_copy_face);

        // 1. Check the input face first (`cxx:598-604`).
        self.check_face();
        if self.my_is_done {
            return;
        }

        // 2. Make all wires connected (`cxx:606-626`).
        let wires = explorer_shapes(&self.my_copy_face.0, ShapeType::Wire);
        for w in &wires {
            self.my_current_wire = w.clone();
            if !self.my_copy_avoid_map.contains(w) {
                let priz = self.connect_wire_copy(false);
                if priz != 0 {
                    // This face contains a wire that is not closed.
                    self.my_is_done = false;
                    self.my_error_status = ERROR_NOT_CONNECTABLE;
                    return;
                }
            }
        }

        // 3. Define the outer wire (`cxx:628-636`).
        let Some(an_outer_wire) = self.outer_wire() else {
            self.my_is_done = false;
            // Can't find the outer wire.
            self.my_error_status = ERROR_NO_OUTER_WIRE;
            return;
        };

        // 4. Moving the outer wire and the other wires in 2d space
        //    (`cxx:638-639`).
        self.move_wires_2d(&an_outer_wire);

        self.my_is_done = true;
        self.my_error_status = 0;
    }

    /// `CheckFace()` (`cxx:226-336`).
    fn check_face(&mut self) {
        // I. Is the surface periodic (`cxx:228-247`).
        let Some(surf) = BRepTool::face_surface_world(&self.my_copy_face) else {
            // `BRep_Tool::Surface` on a face without a surface dereferences a
            // null handle in OCCT; answer the "not periodic" status.
            self.my_is_done = true;
            self.my_error_status = ERROR_NOT_PERIODIC;
            return;
        };
        if !(surf.is_u_periodic() || surf.is_v_periodic()) {
            self.my_is_done = true;
            self.my_error_status = ERROR_NOT_PERIODIC;
            return;
        }
        // Don't treat torus surfaces (NIZHNY-MZV, `cxx:239-247`).
        if surf.is_u_periodic() && surf.is_v_periodic() {
            self.my_is_done = true;
            self.my_error_status = ERROR_NOT_PERIODIC;
            return;
        }

        // II. Has each wire at least one "licit" edge (`cxx:249-292`).
        let mut a_void_wires = IndexedShapeSet::new();
        let wires = explorer_shapes(&self.my_copy_face.0, ShapeType::Wire);
        for a_wire in &wires {
            if self.my_copy_avoid_map.contains(a_wire) {
                continue;
            }
            let mut an_edge_exists = false;
            for an_edge in explorer_shapes(a_wire, ShapeType::Edge) {
                if self.my_copy_avoid_map.contains(&an_edge) {
                    an_edge_exists = true;
                    break;
                }
            }
            if !an_edge_exists {
                // This is the wire without any "Right" information.
                a_void_wires.add(a_wire);
            }
        }

        // For void wires we assume the first edge is the non-movable one
        // (PKV 24-Feb-2000, `cxx:282-292`).
        let a_nb_void_wires = a_void_wires.extent();
        for i in 1..=a_nb_void_wires {
            let Some(w) = a_void_wires.find_key(i).cloned() else {
                continue;
            };
            let mut a_em = IndexedShapeSet::new();
            for e in explorer_shapes(&w, ShapeType::Edge) {
                a_em.add(&e);
            }
            if a_em.extent() > 0 {
                if let Some(first) = a_em.find_key(1).cloned() {
                    self.my_copy_avoid_map.add(&first);
                }
            }
        }

        // III. Check all wires to know whether they are closed (`cxx:294-310`).
        for a_wire in &wires {
            self.my_current_wire = a_wire.clone();
            let priz = self.make_right_wire();
            if priz != 0 {
                // This face contains a wire that is not closed.
                self.my_is_done = true;
                self.my_error_status = ERROR_NOT_CONNECTABLE;
                return;
            }
        }

        // VI. Check the connectability of the wires in 2d (`cxx:312-335`).
        let mut a_copy_face = Face::new();
        let mut a_copy_avoid_map = IndexedShapeSet::new();
        let copy_face = self.my_copy_face.clone();
        let copy_avoid = self.my_copy_avoid_map.clone();
        self.build_copy_data(
            &copy_face,
            &copy_avoid,
            &mut a_copy_face,
            &mut a_copy_avoid_map,
            false,
        );

        let copy_wires = explorer_shapes(&a_copy_face.0, ShapeType::Wire);
        for w in &copy_wires {
            self.my_current_wire = w.clone();
            if !a_copy_avoid_map.contains(w) {
                let priz = self.connect_wire(&mut a_copy_face, &a_copy_avoid_map, true);
                if priz != 0 {
                    self.my_is_done = true;
                    self.my_error_status = ERROR_WRONG_CONNECTABILITY;
                    return;
                }
            }
        }

        // Face seems to be OK.
    }

    /// `ConnectWire(myCopyFace, myCopyAvoidMap, false)` — the `Perform` caller
    /// (`cxx:617`). `myCopyFace` is a member there, so it is passed by value
    /// here to satisfy the borrow of `&mut self`.
    fn connect_wire_copy(&mut self, a_try_flag: bool) -> i32 {
        let mut a_copy_face = self.my_copy_face.clone();
        let avoid = self.my_copy_avoid_map.clone();
        let priz = self.connect_wire(&mut a_copy_face, &avoid, a_try_flag);
        self.my_copy_face = a_copy_face;
        priz
    }

    /// `MakeRightWire()` (`cxx:646-683`).
    pub(super) fn make_right_wire(&mut self) -> i32 {
        self.my_ordered_wire_list.clear();

        // 1. Real amount of the edges on aWire (`cxx:651-659`).
        let a_nb_edges_really = explorer_shapes(&self.my_current_wire, ShapeType::Edge).len() as i32;

        // 2. `TopOpeBRepBuild_Tools2d::Path` (`cxx:661-668`).
        let a_l = tools2d::path(&self.my_current_wire);
        if a_l.len() as i32 != a_nb_edges_really {
            self.my_error_status = ERROR_NOT_PERIODIC;
            return 1;
        }

        // The path list comes back in reverse order; prepend to restore it
        // (NIZNHY-PKV, `cxx:669-681`).
        let mut a_fl: Vec<TopoShape> = Vec::with_capacity(a_l.len());
        for v in &a_l {
            a_fl.insert(0, v.clone());
        }
        self.my_ordered_wire_list = a_fl;
        0
    }
}
