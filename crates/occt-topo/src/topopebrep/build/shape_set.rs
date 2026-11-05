//! `TopOpeBRepBuild_ShapeSet` — an exploration of a set of shapes to build
//! faces or solids.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_ShapeSet.hxx:27-152`,
//! `.cxx:22-499`.
//!
//! To build faces the shapes are wires and the elements are edges; to build
//! solids the shapes are shells and the elements are faces. The set stores a
//! list of shapes, a list of elements to start reconstructions, and a map to
//! search neighbours: the map stores the connection between elements through
//! sub-shapes of type `mySubShapeType`.
//!
//! Translation notes:
//! - `NCollection_List<TopoDS_Shape>` is [`Vec<TopoShape>`] plus a cursor
//!   field standing in for the member `NCollection_List::Iterator`
//!   (`hxx:97-110`), which is what `NCollection_List::Iterator::Initialize`
//!   reset to does.
//! - `mySubShapeMap` is an
//!   `NCollection_IndexedDataMap<Shape, List<Shape>, TopTools_ShapeMapHasher>`
//!   and maps to [`IndexedShapeData<Vec<TopoShape>>`]; the three
//!   `NCollection_IndexedMap<Shape>` sets (`myOMSS`, `myOMES`, `myOMSH`) map
//!   to [`IndexedShapeSet`], whose key is `TopTools_ShapeMapHasher`
//!   (`TShape` + `Location`).
//! - `MakeNeighboursList` returns `const List&` into `mySubShapeMap`. The base
//!   implementation has no other reader of the map during the neighbour walk,
//!   so the port copies the list into `myCurrentShapeNeighbours` (the field
//!   `hxx:110` already declares) and iterates that.

use crate::abs::ShapeType;
use crate::brep_check::Analyzer;
use crate::shape::TopoShape;
use crate::topopebrep::ds::data_structure::{IndexedShapeData, IndexedShapeSet};
use crate::topopebrep::tool::shape_explorer::ShapeExplorer;

/// `class TopOpeBRepBuild_ShapeSet` (`hxx:51`).
pub struct ShapeSet {
    /// `myShapeType` (`hxx:115`).
    pub(super) my_shape_type: ShapeType,
    /// `mySubShapeType` (`hxx:116`).
    pub(super) my_sub_shape_type: ShapeType,
    /// `mySubShapeExplorer` (`hxx:117`).
    pub(super) my_sub_shape_explorer: ShapeExplorer,
    /// `myStartShapes` (`hxx:118`).
    pub(super) my_start_shapes: Vec<TopoShape>,
    /// `myStartShapesIter` (`hxx:119`).
    pub(super) my_start_shapes_iter: usize,
    /// `mySubShapeMap` (`hxx:120-121`).
    pub(super) my_sub_shape_map: IndexedShapeData<Vec<TopoShape>>,
    /// `myIncidentShapesIter` (`hxx:122`) — a cursor into
    /// `myCurrentShapeNeighbours`.
    pub(super) my_incident_shapes_iter: usize,
    /// `myShapes` (`hxx:123`).
    pub(super) my_shapes: Vec<TopoShape>,
    /// `myShapesIter` (`hxx:124`).
    pub(super) my_shapes_iter: usize,
    /// `myCurrentShape` (`hxx:125`).
    pub(super) my_current_shape: Option<TopoShape>,
    /// `myCurrentShapeNeighbours` (`hxx:126`).
    pub(super) my_current_shape_neighbours: Vec<TopoShape>,
    /// `myDEBNumber` (`hxx:127`).
    pub(super) my_deb_number: i32,
    /// `myDEBName` (`hxx:128`).
    pub(super) my_deb_name: String,
    /// `myOMSS` (`hxx:129`).
    pub(super) my_omss: IndexedShapeSet,
    /// `myOMES` (`hxx:130`).
    pub(super) my_omes: IndexedShapeSet,
    /// `myOMSH` (`hxx:131`).
    pub(super) my_omsh: IndexedShapeSet,
    /// `myCheckShape` (`hxx:132`).
    pub(super) my_check_shape: bool,
}

impl ShapeSet {
    /// `TopOpeBRepBuild_ShapeSet(SubShapeType, checkshape)` (`cxx:32-49`).
    ///
    /// `cxx:47` overwrites `myCheckShape` with `false` ("temporary NYI"), so
    /// the `checkshape` argument has no effect, exactly as in OCCT.
    pub fn new(sub_shape_type: ShapeType, checkshape: bool) -> Self {
        let my_shape_type = if sub_shape_type == ShapeType::Edge {
            ShapeType::Face
        } else if sub_shape_type == ShapeType::Vertex {
            ShapeType::Edge
        } else {
            // cxx:43-46: Standard_ProgramError("ShapeSet : bad ShapeType")
            panic!("ShapeSet : bad ShapeType");
        };
        let _ = checkshape;
        Self {
            my_shape_type,
            my_sub_shape_type: sub_shape_type,
            my_sub_shape_explorer: ShapeExplorer::new(),
            my_start_shapes: Vec::new(),
            my_start_shapes_iter: 0,
            my_sub_shape_map: IndexedShapeData::new(),
            my_incident_shapes_iter: 0,
            my_shapes: Vec::new(),
            my_shapes_iter: 0,
            my_current_shape: None,
            my_current_shape_neighbours: Vec::new(),
            my_deb_number: 0,
            my_deb_name: String::new(),
            my_omss: IndexedShapeSet::new(),
            my_omes: IndexedShapeSet::new(),
            my_omsh: IndexedShapeSet::new(),
            my_check_shape: false, // cxx:47 — temporary NYI
        }
    }

    /// `AddShape(S)` (`cxx:53-66`).
    pub fn add_shape(&mut self, s: &TopoShape) {
        if !self.check_shape_of(s, false) {
            return;
        }
        self.process_add_shape(s);
    }

    /// `AddStartElement(S)` (`cxx:70-83`).
    pub fn add_start_element(&mut self, s: &TopoShape) {
        if !self.check_shape_of(s, false) {
            return;
        }
        self.process_add_start_element(s);
    }

    /// `AddElement(S)` (`cxx:87-100`).
    pub fn add_element(&mut self, s: &TopoShape) {
        if !self.check_shape_of(s, false) {
            return;
        }
        self.process_add_element(s);
    }

    /// `ProcessAddShape(S)` (`cxx:104-111`).
    pub fn process_add_shape(&mut self, s: &TopoShape) {
        if !self.my_omsh.contains(s) {
            self.my_omsh.add(s);
            self.my_shapes.push(s.clone());
        }
    }

    /// `ProcessAddStartElement(S)` (`cxx:115-123`).
    pub fn process_add_start_element(&mut self, s: &TopoShape) {
        if !self.my_omss.contains(s) {
            self.my_omss.add(s);
            self.my_start_shapes.push(s.clone());
            self.process_add_element(s);
        }
    }

    /// `ProcessAddElement(S)` (`cxx:127-146`).
    pub fn process_add_element(&mut self, s: &TopoShape) {
        if self.my_omes.contains(s) {
            return;
        }
        self.my_omes.add(s);
        let mut ex = ShapeExplorer::with_shape(s, self.my_sub_shape_type, ShapeType::Shape);
        while ex.more() {
            let subshape = match ex.current() {
                Some(c) => c.clone(),
                None => break,
            };
            if !self.my_sub_shape_map.contains(&subshape) {
                self.my_sub_shape_map.add(&subshape, Vec::new());
            }
            if let Some(l) = self.my_sub_shape_map.change_from_key(&subshape) {
                l.push(s.clone());
            }
            ex.next();
        }
    }

    /// `StartElements()` (`cxx:150-153`).
    pub fn start_elements(&self) -> &[TopoShape] {
        &self.my_start_shapes
    }

    /// `InitShapes()` (`cxx:157-160`).
    pub fn init_shapes(&mut self) {
        self.my_shapes_iter = 0;
    }

    /// `MoreShapes()` (`cxx:164-168`).
    pub fn more_shapes(&self) -> bool {
        self.my_shapes_iter < self.my_shapes.len()
    }

    /// `NextShape()` (`cxx:172-175`).
    pub fn next_shape(&mut self) {
        self.my_shapes_iter += 1;
    }

    /// `Shape()` (`cxx:179-183`).
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_shapes.get(self.my_shapes_iter)
    }

    /// `InitStartElements()` (`cxx:187-190`).
    pub fn init_start_elements(&mut self) {
        self.my_start_shapes_iter = 0;
    }

    /// `MoreStartElements()` (`cxx:194-198`).
    pub fn more_start_elements(&self) -> bool {
        self.my_start_shapes_iter < self.my_start_shapes.len()
    }

    /// `NextStartElement()` (`cxx:202-205`).
    pub fn next_start_element(&mut self) {
        self.my_start_shapes_iter += 1;
    }

    /// `StartElement()` (`cxx:209-213`).
    pub fn start_element(&self) -> Option<&TopoShape> {
        self.my_start_shapes.get(self.my_start_shapes_iter)
    }

    /// `InitNeighbours(S)` (`cxx:217-222`).
    pub fn init_neighbours(&mut self, s: &TopoShape) {
        self.my_sub_shape_explorer
            .init(s, self.my_sub_shape_type, ShapeType::Shape);
        self.my_current_shape = Some(s.clone());
        self.find_neighbours();
    }

    /// `MoreNeighbours()` (`cxx:226-230`).
    pub fn more_neighbours(&self) -> bool {
        self.my_incident_shapes_iter < self.my_current_shape_neighbours.len()
    }

    /// `NextNeighbour()` (`cxx:234-247`).
    ///
    /// `NextNeighbour` is non-virtual but calls the virtual `FindNeighbours`
    /// at `cxx:243`; callers that hold a subclass of `ShapeSet` must go through
    /// [`ShapeSetSource::next_neighbour`], which dispatches that call.
    pub fn next_neighbour_base(&mut self) {
        self.my_incident_shapes_iter += 1;
        let no_is_more = !self.more_neighbours();
        if no_is_more {
            if self.my_sub_shape_explorer.more() {
                self.my_sub_shape_explorer.next();
                self.find_neighbours();
            }
        }
    }

    /// `Neighbour()` (`cxx:251-255`).
    pub fn neighbour(&self) -> Option<&TopoShape> {
        self.my_current_shape_neighbours
            .get(self.my_incident_shapes_iter)
    }

    /// `ChangeStartShapes()` (`cxx:259-262`).
    pub fn change_start_shapes(&mut self) -> &mut Vec<TopoShape> {
        &mut self.my_start_shapes
    }

    /// `FindNeighbours()` (`cxx:266-296`).
    pub fn find_neighbours(&mut self) {
        while self.my_sub_shape_explorer.more() {
            let v = match self.my_sub_shape_explorer.current() {
                Some(c) => c.clone(),
                None => break,
            };
            let e = self.my_current_shape.clone();
            let l = self.make_neighbours_list(e.as_ref(), &v);
            self.my_current_shape_neighbours = l;
            self.my_incident_shapes_iter = 0;
            if self.more_neighbours() {
                break;
            } else {
                self.my_sub_shape_explorer.next();
            }
        }
    }

    /// `MakeNeighboursList(Earg, Varg)` (`cxx:302-308`): the neighbours of
    /// `Earg` through `Varg`, i.e. the incident list stored under `Varg`.
    pub fn make_neighbours_list(
        &self,
        _e_arg: Option<&TopoShape>,
        v_arg: &TopoShape,
    ) -> Vec<TopoShape> {
        // cxx:306 `mySubShapeMap.FindFromKey(Varg)`; a shape that was never
        // added has no entry and answers the empty list.
        self.my_sub_shape_map
            .find_from_key(v_arg)
            .cloned()
            .unwrap_or_default()
    }

    /// `MaxNumberSubShape(Shape)` (`cxx:312-335`): the largest number of
    /// elements incident to one of the sub-shapes of `Shape`.
    pub fn max_number_sub_shape(&self, shape: &TopoShape) -> i32 {
        let mut m = 0;
        let mut se = ShapeExplorer::with_shape(shape, self.my_sub_shape_type, ShapeType::Shape);
        while se.more() {
            let sub_shape = match se.current() {
                Some(c) => c.clone(),
                None => break,
            };
            if !self.my_sub_shape_map.contains(&sub_shape) {
                se.next();
                continue;
            }
            let i = match self.my_sub_shape_map.find_from_key(&sub_shape) {
                Some(l) => l.len() as i32,
                None => 0,
            };
            if i > m {
                m = i;
            }
            se.next();
        }
        m
    }

    /// `CheckShape(const bool checkshape)` (`cxx:339-352`).
    pub fn set_check_shape(&mut self, checkshape: bool) {
        self.my_check_shape = checkshape;
    }

    /// `CheckShape() const` (`cxx:356-359`).
    pub fn check_shape(&self) -> bool {
        self.my_check_shape
    }

    /// `CheckShape(S, checkgeom)` (`cxx:363-373`).
    ///
    /// `BRepCheck_Analyzer ana(S, checkgeom)` (`cxx:370`); the port's
    /// [`Analyzer::new`] has no geometry-controls switch yet, so the
    /// `checkgeom` flag is currently not forwarded.
    pub fn check_shape_of(&mut self, s: &TopoShape, _check_geom: bool) -> bool {
        if !self.my_check_shape {
            return true;
        }
        Analyzer::new(s.clone()).is_valid()
    }

    /// `DumpName(OS, str)` (`cxx:377-381`) — debug output only.
    pub fn dump_name(&self, _str: &str) {}

    /// `DumpCheck(OS, str, S, chk)` (`cxx:386-430`) — debug output only.
    pub fn dump_check(&self, _str: &str, _s: &TopoShape, _chk: bool) {}

    /// `DumpSS()` (`cxx:434`).
    pub fn dump_ss(&self) {}

    /// `DumpBB()` (`cxx:438`).
    pub fn dump_bb(&self) {}

    /// `DEBName(N)` (`cxx:442-445`).
    pub fn set_deb_name(&mut self, n: &str) {
        self.my_deb_name = n.to_string();
    }

    /// `DEBName()` (`cxx:449-452`).
    pub fn deb_name(&self) -> &str {
        &self.my_deb_name
    }

    /// `DEBNumber(I)` (`cxx:456-459`).
    pub fn set_deb_number(&mut self, i: i32) {
        self.my_deb_number = i;
    }

    /// `DEBNumber()` (`cxx:463-466`).
    pub fn deb_number(&self) -> i32 {
        self.my_deb_number
    }

    /// `SName(S, sb, sa)` (`cxx:470-475`) — empty without `OCCT_DEBUG`.
    pub fn s_name(&self, _s: &TopoShape, _sb: &str, _sa: &str) -> String {
        String::new()
    }

    /// `SNameori(S, sb, sa)` (`cxx:479-485`) — empty without `OCCT_DEBUG`.
    pub fn s_name_ori(&self, _s: &TopoShape, _sb: &str, _sa: &str) -> String {
        String::new()
    }

    /// `SName(L, sb, sa)` (`cxx:489-495`) — empty without `OCCT_DEBUG`.
    pub fn s_name_list(&self, _l: &[TopoShape], _sb: &str, _sa: &str) -> String {
        String::new()
    }

    /// `SNameori(L, sb, sa)` (`cxx:499-506`) — empty without `OCCT_DEBUG`.
    pub fn s_name_ori_list(&self, _l: &[TopoShape], _sb: &str, _sa: &str) -> String {
        String::new()
    }
}

/// The virtual surface of `TopOpeBRepBuild_ShapeSet` (`hxx:104-190`).
///
/// `BlockBuilder::MakeBlock` (`BlockBuilder.cxx:65-127`) drives a `ShapeSet`
/// by reference and relies on the virtual calls `AddStartElement`,
/// `AddElement`, `AddShape`, `InitNeighbours` and `FindNeighbours` to reach the
/// subclass (`TopOpeBRepBuild_WireEdgeSet`, `TopOpeBRepBuild_ShellFaceSet`).
/// The non-virtual members `MoreNeighbours`, `NextNeighbour` and `Neighbour`
/// are the provided methods below; `NextNeighbour` re-dispatches to
/// `FindNeighbours` exactly as `cxx:243` does.
pub trait ShapeSetSource {
    fn base(&self) -> &ShapeSet;
    fn base_mut(&mut self) -> &mut ShapeSet;

    fn add_shape(&mut self, s: &TopoShape);
    fn add_start_element(&mut self, s: &TopoShape);
    fn add_element(&mut self, s: &TopoShape);
    fn init_neighbours(&mut self, s: &TopoShape);
    fn find_neighbours(&mut self);
    fn make_neighbours_list(&mut self, e: Option<&TopoShape>, v: &TopoShape) -> Vec<TopoShape>;

    /// `ChangeStartShapes()` (`cxx:259-262`).
    fn change_start_shapes(&mut self) -> &mut Vec<TopoShape> {
        ShapeSet::change_start_shapes(self.base_mut())
    }

    /// `StartElements()` (`cxx:65-70`).
    fn start_elements(&self) -> &[TopoShape] {
        ShapeSet::start_elements(self.base())
    }

    /// `InitStartElements()` (`cxx:108-111`).
    fn init_start_elements(&mut self) {
        ShapeSet::init_start_elements(self.base_mut());
    }

    /// `MoreStartElements()` (`cxx:113-116`).
    fn more_start_elements(&mut self) -> bool {
        ShapeSet::more_start_elements(self.base_mut())
    }

    /// `NextStartElement()` (`cxx:118-121`).
    fn next_start_element(&mut self) {
        ShapeSet::next_start_element(self.base_mut());
    }

    /// `StartElement()` (`cxx:123-129`).
    fn start_element(&self) -> Option<&TopoShape> {
        ShapeSet::start_element(self.base())
    }

    /// `MaxNumberSubShape(S)` (`cxx:216-219`).
    fn max_number_sub_shape(&mut self, s: &TopoShape) -> usize {
        ShapeSet::max_number_sub_shape(self.base_mut(), s) as usize
    }

    /// `InitNeighbours(S)` is virtual; the base body is `cxx:222-229`.

    /// `MoreNeighbours()` (`cxx:231-233`).
    fn more_neighbours(&self) -> bool {
        ShapeSet::more_neighbours(self.base())
    }

    /// `NextNeighbour()` (`cxx:234-247`) with the virtual `FindNeighbours` at
    /// `cxx:243`.
    fn next_neighbour(&mut self) {
        self.base_mut().my_incident_shapes_iter += 1;
        if !ShapeSet::more_neighbours(self.base()) && self.base().my_sub_shape_explorer.more() {
            self.base_mut().my_sub_shape_explorer.next();
            self.find_neighbours();
        }
    }

    /// `Neighbour()` (`cxx:251-255`).
    fn neighbour(&self) -> Option<&TopoShape> {
        ShapeSet::neighbour(self.base())
    }

    /// `InitShapes()` (`cxx:358-362`).
    fn init_shapes(&mut self) {
        ShapeSet::init_shapes(self.base_mut());
    }

    /// `MoreShapes()` (`cxx:364-367`).
    fn more_shapes(&mut self) -> bool {
        ShapeSet::more_shapes(self.base_mut())
    }

    /// `NextShape()` (`cxx:369-372`).
    fn next_shape(&mut self) {
        ShapeSet::next_shape(self.base_mut());
    }

    /// `Shape()` (`cxx:374-377`).
    fn shape(&self) -> Option<&TopoShape> {
        ShapeSet::shape(self.base())
    }
}

impl ShapeSetSource for ShapeSet {
    fn base(&self) -> &ShapeSet {
        self
    }

    fn base_mut(&mut self) -> &mut ShapeSet {
        self
    }

    fn add_shape(&mut self, s: &TopoShape) {
        ShapeSet::add_shape(self, s);
    }

    fn add_start_element(&mut self, s: &TopoShape) {
        ShapeSet::add_start_element(self, s);
    }

    fn add_element(&mut self, s: &TopoShape) {
        ShapeSet::add_element(self, s);
    }

    fn init_neighbours(&mut self, s: &TopoShape) {
        ShapeSet::init_neighbours(self, s);
    }

    fn find_neighbours(&mut self) {
        ShapeSet::find_neighbours(self);
    }

    fn make_neighbours_list(&mut self, e: Option<&TopoShape>, v: &TopoShape) -> Vec<TopoShape> {
        ShapeSet::make_neighbours_list(self, e, v)
    }
}
