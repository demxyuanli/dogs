//! `TopOpeBRepTool_BoxSort` — the box sorter of the classic boolean engine.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepTool/TopOpeBRepTool_BoxSort.hxx:34-89`,
//! `.cxx:53-356`.
//!
//! Translation notes:
//! - `Bnd_BoundSortBox` is [`occt_core::bnd::bound_sort_box::BoundSortBox`];
//!   both index their box array 1-based with a void slot 0, so `myHAI` /
//!   `myHAB` are plain `Vec`s here and `Compare` answers 1-based positions.
//! - `NCollection_List<int>::Iterator myIterator` (`hxx:84`) holds a reference
//!   into `myBSB`'s last result list. The port copies that result into
//!   [`BoxSort::iter_list`] and keeps a cursor; the list is replaced only by
//!   the next [`BoxSort::compare`], exactly like OCCT.
//! - `Compare`'s plane arm needs the `GeomAdaptor_Surface::GetType() == Plane`
//!   test and `GeomAdaptor_Surface::Plane()`; both recurse into a
//!   `Geom_RectangularTrimmedSurface` basis, which
//!   [`crate::geom_bnd_lib_surface3d::surface_kind`] and
//!   [`adaptor_plane`] do.
//! - `BS::Box3` raises `Standard_ProgramError`; the port returns that message.
//! - `Box(S)` is non-const in OCCT (`cxx:333`) because it may lazily attach the
//!   `HBoxTool` (`cxx:335-338`); the port keeps `&mut self` for that reason.

use occt_core::bnd::bound_sort_box::BoundSortBox;
use occt_core::bnd::box3d::BndBox;
use occt_core::gp::GpPln;
use occt_geom::surface::Surface;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::geom_bnd_lib_surface3d::{surface_kind, GeomAbsSurfaceType};
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{edge_vertices, explorer_shapes_avoiding};
use crate::topopebrep::tool::h_box_tool::HBoxTool;

/// `GeomAdaptor_Surface::Plane()`: the plane of a plane, of a trimmed plane, or
/// of a plane reached through any chain of `Geom_RectangularTrimmedSurface`
/// bases (`GeomAdaptor_Surface::load` recurses the same way, `cxx:423-430`).
fn adaptor_plane(s: &dyn Surface) -> Option<GpPln> {
    if let Some(p) = s.gp_pln() {
        return Some(p);
    }
    adaptor_plane(s.rectangular_trimmed_basis()?.as_ref())
}

/// `class TopOpeBRepTool_BoxSort` (`hxx:34-89`).
#[derive(Default)]
pub struct BoxSort {
    /// `myCOB` (`hxx:81`).
    cob: BndBox,
    /// `myBSB` (`hxx:82`).
    bsb: BoundSortBox,
    /// The list `myIterator` walks (`hxx:84`), as 1-based box indexes.
    iter_list: Vec<usize>,
    /// The iterator cursor into `iter_list`.
    iter_pos: usize,
    /// `myLastCompareShape` (`hxx:85`).
    last_compare_shape: Option<TopoShape>,
    /// `myLastCompareShapeBox` (`hxx:86`).
    last_compare_shape_box: BndBox,
    /// `myHBT` (`hxx:87`).
    hbt: Option<HBoxTool>,
    /// `myHAB` (`hxx:88`), 1-based through `BoundSortBox`.
    hab: Vec<BndBox>,
    /// `myHAI` (`hxx:89`), the `HBoxTool` index of each `hab` entry.
    hai: Vec<i32>,
}

impl BoxSort {
    /// `TopOpeBRepTool_BoxSort()` (`cxx:53`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `TopOpeBRepTool_BoxSort(HBT)` (`cxx:59-62`).
    pub fn with_h_box_tool(hbt: HBoxTool) -> Self {
        let mut s = Self::new();
        s.set_h_box_tool(hbt);
        s
    }

    /// `SetHBoxTool(HBT)` (`cxx:78-81`).
    pub fn set_h_box_tool(&mut self, hbt: HBoxTool) {
        self.hbt = Some(hbt);
    }

    /// `HBoxTool()` (`cxx:85-88`).
    pub fn h_box_tool(&self) -> Option<&HBoxTool> {
        self.hbt.as_ref()
    }

    /// `Clear()` (`cxx:92-96`): `myCOB.SetVoid()` only; the box map is kept.
    pub fn clear(&mut self) {
        self.cob.set_void();
        self.iter_list.clear();
    }

    /// `myListIterator.Initialize(anEmptyList)` (`TopOpeBRep_ShapeScanner.cxx:60-62`):
    /// reset the iterator without touching `myCOB` / `myBSB`.
    pub fn reset_iterator(&mut self) {
        self.iter_list.clear();
        self.iter_pos = 0;
    }

    /// `AddBoxes(S, TS, TA)` (`cxx:102-110`), creating `myHBT` when absent.
    pub fn add_boxes(&mut self, s: &TopoShape, ts: ShapeType, ta: ShapeType) -> Result<(), String> {
        if self.hbt.is_none() {
            self.hbt = Some(HBoxTool::new());
        }
        let hbt = self.hbt.as_mut().expect("HBoxTool set above");
        hbt.add_boxes(s, ts, ta)
    }

    /// `MakeHAB(S, TS, TA)` (`cxx:115-160`): build `myHAB` / `myHAI` from the
    /// `TS` sub-shapes of `S`, adding a box to `myHBT` when one is missing.
    pub fn make_hab(&mut self, s: &TopoShape, ts: ShapeType, ta: ShapeType) -> Result<(), String> {
        let subs = explorer_shapes_avoiding(s, ts, ta);
        self.hab.clear();
        self.hai.clear();
        for ss in subs {
            if self.hbt.is_none() {
                self.hbt = Some(HBoxTool::new());
            }
            let hbt = self.hbt.as_mut().expect("HBoxTool set above");
            if !hbt.has_box(&ss) {
                hbt.add_box(&ss)?;
            }
            let im = hbt.index_of(&ss);
            let b = *hbt.box_of(&ss)?;
            self.hai.push(im);
            self.hab.push(b);
        }
        Ok(())
    }

    /// `HAB()` (`cxx:165-168`).
    pub fn hab(&self) -> &[BndBox] {
        &self.hab
    }

    /// `MakeHABCOB(HAB, COB)` (`cxx:172-183`): `COB` is the union of the boxes.
    pub fn make_hab_cob(hab: &[BndBox], cob: &mut BndBox) {
        cob.set_void();
        for b in hab {
            cob.add_box(b);
        }
    }

    /// `HABShape(I)` (`cxx:186-197`): the shape of box `i`, `BS::Box3` outside
    /// `[1, |myHAI|]`.
    pub fn hab_shape(&self, i: i32) -> Result<&TopoShape, String> {
        let iu = self.hai.len() as i32;
        if i < 1 || i > iu {
            return Err("BS::Box3".to_string());
        }
        let im = self.hai[(i - 1) as usize];
        match self.hbt.as_ref() {
            Some(hbt) => hbt.shape(im),
            None => Err("BS::Box3".to_string()),
        }
    }

    /// `MakeCOB(S, TS, TA)` (`cxx:200-215`).
    pub fn make_cob(&mut self, s: &TopoShape, ts: ShapeType, ta: ShapeType) -> Result<(), String> {
        self.make_hab(s, ts, ta)?;
        Self::make_hab_cob(&self.hab, &mut self.cob);
        self.bsb.initialize_with_enclosing(&self.cob, &self.hab);
        Ok(())
    }

    /// `AddBoxesMakeCOB(S, TS, TA)` (`cxx:219-224`).
    pub fn add_boxes_make_cob(
        &mut self,
        s: &TopoShape,
        ts: ShapeType,
        ta: ShapeType,
    ) -> Result<(), String> {
        self.add_boxes(s, ts, ta)?;
        self.make_cob(s, ts, ta)
    }

    /// `Compare(S)` (`cxx:228-318`): re-initialise the iterator over the boxes
    /// touched by `S`.
    pub fn compare(&mut self, s: &TopoShape) -> Result<(), String> {
        if self.hbt.is_none() {
            self.hbt = Some(HBoxTool::new());
        }

        let mut plane = GpPln::default();
        let mut is_plane = false;
        let t = s.shape_type();
        {
            let hbt = self.hbt.as_mut().expect("HBoxTool set above");
            if !hbt.has_box(s) {
                hbt.add_box(s)?;
            }
        }

        self.last_compare_shape = Some(s.clone());
        self.last_compare_shape_box.set_void();

        if t == ShapeType::Face {
            // TopoDS::Face(S), then BRep_Tool::NaturalRestriction(F).
            if let Some(f) = Face::wrap(s.clone()) {
                let natu = BRepTool::natural_restriction(&f);
                if natu {
                    if let Some(surf) = BRepTool::face_surface(&f) {
                        is_plane = surface_kind(surf.as_ref()) == GeomAbsSurfaceType::Plane;
                        if is_plane {
                            plane = adaptor_plane(surf.as_ref()).unwrap_or_default();
                        } else {
                            self.last_compare_shape_box = self.face_box(&f)?;
                        }
                    } else {
                        // OCCT builds a GeomAdaptor_Surface on a null handle
                        // here; the port falls back to the box arm.
                        self.last_compare_shape_box = self.face_box(&f)?;
                    }
                } else {
                    self.last_compare_shape_box = self.face_box(&f)?;
                }
            }
        } else if t == ShapeType::Edge {
            if let Some(e) = Edge::wrap(s.clone()) {
                let (v1, v2) = edge_vertices(&e);
                let perso = v1.is_none() || v2.is_none();
                if perso {
                    HBoxTool::compute_box_on_vertices(
                        &e.0,
                        &mut self.last_compare_shape_box,
                    );
                } else {
                    self.last_compare_shape_box = self.edge_box(&e)?;
                }
            }
        }

        if is_plane {
            self.iter_list = self.bsb.compare_plane(&plane).to_vec();
        } else {
            self.iter_list = self.bsb.compare(&self.last_compare_shape_box).to_vec();
        }
        self.iter_pos = 0;
        Ok(())
    }

    /// `myHBT->Box(F)` for the face arm of [`Self::compare`].
    fn face_box(&self, f: &Face) -> Result<BndBox, String> {
        self.hbt
            .as_ref()
            .ok_or_else(|| "HBT::Box1".to_string())?
            .box_of(&f.0)
            .map(|b| *b)
    }

    /// `myHBT->Box(E)` for the edge arm of [`Self::compare`].
    fn edge_box(&self, e: &Edge) -> Result<BndBox, String> {
        self.hbt
            .as_ref()
            .ok_or_else(|| "HBT::Box1".to_string())?
            .box_of(&e.0)
            .map(|b| *b)
    }

    /// `TouchedShape(LI)` (`cxx:322-327`).
    pub fn touched_shape(&self, i: i32) -> Result<&TopoShape, String> {
        self.hab_shape(i)
    }

    /// `Box(S)` (`cxx:332-356`).
    pub fn box_of(&mut self, s: &TopoShape) -> Result<BndBox, String> {
        if self.hbt.is_none() {
            self.hbt = Some(HBoxTool::new());
        }
        let hbt = self.hbt.as_ref().expect("HBoxTool set above");
        if hbt.has_box(s) {
            return hbt.box_of(s).map(|b| *b);
        }
        if let Some(last) = self.last_compare_shape.as_ref() {
            if last.same_tshape(s) && !self.last_compare_shape_box.is_void() {
                return Ok(self.last_compare_shape_box);
            }
        }
        hbt.box_of(s).map(|b| *b)
    }

    /// `myListIterator.More()` (`cxx:87-90`).
    pub fn more(&self) -> bool {
        self.iter_pos < self.iter_list.len()
    }

    /// `myListIterator.Next()` (`cxx:94-97`).
    pub fn next(&mut self) {
        self.iter_pos += 1;
    }

    /// `TouchedShape(myListIterator)` / `Current()` (`cxx:102-106`).
    pub fn current(&self) -> Option<Result<&TopoShape, String>> {
        let i = *self.iter_list.get(self.iter_pos)?;
        Some(self.touched_shape(i as i32))
    }

    /// `myListIterator.Value()` (`cxx:123-131`), 0 when the iterator is done.
    pub fn index(&self) -> i32 {
        match self.iter_list.get(self.iter_pos) {
            Some(&i) => i as i32,
            None => 0,
        }
    }
}
