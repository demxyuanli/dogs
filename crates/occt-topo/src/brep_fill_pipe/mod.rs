//! `BRepFill_Pipe` — sweep a profile shape along a wire spine
//! (`BRepFill_Pipe.cxx:1-1221`, `.hxx:45-150`).
//!
//! Port status, keeping the OCCT control flow:
//!
//! * Ported: both constructors (`cxx:148-190`), `Perform` (`cxx:193-311`),
//!   `Spine` / `Profile` / `Shape` / `ErrorOnSurface` / `FirstShape` /
//!   `LastShape` (`cxx:317-356`), `Generated` (`cxx:359-366`), `Face`
//!   (`cxx:372-410`), `Edge` (`cxx:416-453`), `Section` (`cxx:458-493`),
//!   `PipeLine` (`cxx:502-526`), the recursive `MakeShape` (`cxx:530-851`),
//!   `FindEdge` (`cxx:860-910`), `FindVertex` (`cxx:918-990`),
//!   `DefineRealSegmax` (`cxx:999-1049`), `RebuildTopOrBottomFace`
//!   (`cxx:1058-1090`), `UpdateTolFromTopOrBottomPCurve` (`cxx:96-146`) and
//!   `BuildHistory` (`cxx:1100-1220`).
//!
//! Nothing in `BRepFill_Pipe.cxx` is left unported. Three reference behaviours
//! have no directly translatable counterpart; each is documented at its call
//! site instead of being replaced by a different rule:
//!
//! * `TopLoc_Location::operator!=` (`cxx:298`) compares `TopLoc_Datum3D`
//!   handles, so the two datums built in `Perform` are always distinct and the
//!   `else myLast = myFirst` arm (`cxx:305-307`) is unreachable.
//! * `Section` (`cxx:487-490`) and `RebuildTopOrBottomFace` (`cxx:1060`,
//!   `:1075`) index `NCollection_Array2` with the column bound in the row
//!   position (and vice versa). Out-of-range reads are undefined in OCCT;
//!   [`ShapeTable`] answers `None` there and the loops skip the slot.
//! * `TopoDS_Iterator::Value()` on an exhausted iterator (`cxx:644-659`,
//!   `:1091`) is undefined in OCCT; the port passes `None`.
//!
//! Translation notes:
//! - `handle(BRepFill_LocationLaw)` is `Option<Arc<dyn LocationLaw>>`; a null
//!   handle is `None`. `mySections` / `myFaces` / `myEdges` are
//!   `Option<ShapeTable>`, i.e. `handle(NCollection_HArray2)` with a null
//!   handle.
//! - `mySpine` / `myProfile` are `Option<Wire>` / `Option<TopoShape>`: a
//!   default-constructed `BRepFill_Pipe` stores null shapes in OCCT.

pub mod history;
pub mod make_shape;
pub mod perform;
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use occt_core::gp::{GpPnt, GpTrsf};
use occt_core::kernel::geomabs::Shape;

use occt_geom::geom_fill::GeomFillTrihedron;

use crate::abs::ShapeType;
use crate::brep_fill_location_law::LocationLaw;
use crate::brep_fill_sweep::helpers::ShapeTable;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Compound, Edge, Face, TopoShape, Vertex, Wire};
use crate::topo_tools::{is_same, shape_map_key, ShapeMapKey};

/// `BRepFill_Pipe` (`BRepFill_Pipe.hxx:48-150`): the class data.
pub struct Pipe {
    /// `mySpine` (`hxx:121`), a null wire in OCCT before `Perform`.
    my_spine: Option<Wire>,
    /// `myProfile` (`hxx:122`).
    my_profile: Option<TopoShape>,
    /// `myShape` (`hxx:123`).
    my_shape: Option<TopoShape>,
    /// `myTrsf` (`hxx:124`): the section placement transform.
    my_trsf: GpTrsf,
    /// `myLoc` (`hxx:125`).
    my_loc: Option<Arc<dyn LocationLaw>>,
    /// `mySections` (`hxx:126`).
    my_sections: Option<ShapeTable>,
    /// `myFaces` (`hxx:127`).
    my_faces: Option<ShapeTable>,
    /// `myEdges` (`hxx:128`).
    my_edges: Option<ShapeTable>,
    /// `myReversedEdges` (`hxx:129`).
    my_reversed_edges: HashSet<ShapeMapKey>,
    /// `myTapes` (`hxx:130-132`).
    my_tapes: HashMap<ShapeMapKey, ShapeTable>,
    /// `myRails` (`hxx:133-135`).
    my_rails: HashMap<ShapeMapKey, ShapeTable>,
    /// `myCurIndexOfSectionEdge` (`hxx:136`).
    my_cur_index_of_section_edge: usize,
    /// `myFirst` (`hxx:137`).
    my_first: Option<TopoShape>,
    /// `myLast` (`hxx:138`).
    my_last: Option<TopoShape>,
    /// `myGenMap` (`hxx:139-140`).
    my_gen_map: HashMap<ShapeMapKey, Vec<TopoShape>>,
    /// `myDegmax` (`hxx:141`).
    my_degmax: i32,
    /// `mySegmax` (`hxx:142`).
    my_segmax: i32,
    /// `myContinuity` (`hxx:143`).
    my_continuity: Shape,
    /// `myMode` (`hxx:144`).
    my_mode: GeomFillTrihedron,
    /// `myForceApproxC1` (`hxx:145`).
    my_force_approx_c1: bool,
    /// `myErrorOnSurf` (`hxx:146`).
    my_error_on_surf: f64,
}

impl Default for Pipe {
    /// `BRepFill_Pipe()` (`cxx:148-158`).
    fn default() -> Self {
        Self {
            my_spine: None,
            my_profile: None,
            my_shape: None,
            my_trsf: GpTrsf::identity(),
            my_loc: None,
            my_sections: None,
            my_faces: None,
            my_edges: None,
            my_reversed_edges: HashSet::new(),
            my_tapes: HashMap::new(),
            my_rails: HashMap::new(),
            // `myCurIndexOfSectionEdge = 1` (`cxx:157`).
            my_cur_index_of_section_edge: 1,
            my_first: None,
            my_last: None,
            my_gen_map: HashMap::new(),
            my_degmax: 11,
            my_segmax: 100,
            my_continuity: Shape::C2,
            my_mode: GeomFillTrihedron::CorrectedFrenet,
            my_force_approx_c1: false,
            my_error_on_surf: 0.0,
        }
    }
}

impl Pipe {
    /// `BRepFill_Pipe()` (`cxx:148-158`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `BRepFill_Pipe(Spine, Profile, aMode, ForceApproxC1, GeneratePartCase)`
    /// (`cxx:161-189`).
    ///
    /// `aMode` only keeps the three trihedra the class can sweep with; every
    /// other value falls back to `CorrectedFrenet` (`cxx:172-189`). The last
    /// parameter is forwarded to `Perform`, which ignores it (`cxx:195`).
    pub fn with_spine(
        spine: &Wire,
        profile: &TopoShape,
        mode: GeomFillTrihedron,
        force_approx_c1: bool,
        generate_part_case: bool,
    ) -> Result<Self, String> {
        let mut this = Self::default();
        this.my_mode = match mode {
            GeomFillTrihedron::Frenet
            | GeomFillTrihedron::CorrectedFrenet
            | GeomFillTrihedron::DiscreteTrihedron => mode,
            _ => GeomFillTrihedron::CorrectedFrenet,
        };
        this.my_continuity = if this.my_mode == GeomFillTrihedron::DiscreteTrihedron {
            Shape::C0
        } else {
            Shape::C2
        };
        this.my_force_approx_c1 = force_approx_c1;
        this.perform(spine, profile, generate_part_case)?;
        Ok(this)
    }

    /// `Spine()` (`cxx:317-320`): a null wire before `Perform` in OCCT.
    pub fn spine(&self) -> Option<&Wire> {
        self.my_spine.as_ref()
    }

    /// `Profile()` (`cxx:323-326`).
    pub fn profile(&self) -> Option<&TopoShape> {
        self.my_profile.as_ref()
    }

    /// `Shape()` (`cxx:331-334`).
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_shape.as_ref()
    }

    /// `ErrorOnSurface()` (`cxx:338-341`).
    pub fn error_on_surface(&self) -> f64 {
        self.my_error_on_surf
    }

    /// `FirstShape()` (`cxx:344-347`).
    pub fn first_shape(&self) -> Option<&TopoShape> {
        self.my_first.as_ref()
    }

    /// `LastShape()` (`cxx:351-354`).
    pub fn last_shape(&self) -> Option<&TopoShape> {
        self.my_last.as_ref()
    }

    /// `Generated(theShape, theList)` (`cxx:359-366`): the shapes generated
    /// from `the_shape`, `IsSame`-keyed as `NCollection_DataMap` is.
    pub fn generated(&self, the_shape: &TopoShape) -> Vec<TopoShape> {
        self.my_gen_map.get(&shape_map_key(the_shape)).cloned().unwrap_or_default()
    }

    /// `Face(ESpine, EProfile)` (`cxx:372-410`): the face swept from an edge of
    /// the spine and an edge of the profile.
    pub fn face(&self, e_spine: &Edge, e_profile: &Edge) -> Result<Face, String> {
        // `cxx:376-379`: a degenerated profile edge answers a null face; the
        // port answers the empty face `TopoShape` builds for that type.
        if BRepTool::is_degenerated(e_profile) {
            return Ok(Face::new());
        }
        // `cxx:385-390`: "Search if EProfile is an edge of myProfile".
        let mut count = 0i32;
        let profile = self.profile().cloned().unwrap_or_else(|| TopoShape::new(ShapeType::Shape));
        let iprof = self.find_edge(&profile, e_profile, &mut count)?;
        if iprof == 0 {
            return Err("Standard_DomainError: BRepFill_Pipe::Face : Edge not in the Profile"
                .to_string());
        }

        // `cxx:397-407`.
        let loc = self.location()?;
        let mut ispin = 0usize;
        for ii in 1..=loc.base().nb_law() {
            if ispin != 0 {
                break;
            }
            if let Some(edge) = loc.base().edge(ii) {
                if is_same(&edge.0, &e_spine.0) {
                    ispin = ii;
                }
            }
        }
        if ispin == 0 {
            return Err("Standard_DomainError: BRepFill_Pipe::Edge  : Edge not in the Spine"
                .to_string());
        }

        // `cxx:409-410`.
        let faces = self
            .my_faces
            .as_ref()
            .ok_or_else(|| "BRepFill_Pipe::Face : the sweep produced no faces".to_string())?;
        let the_face = faces
            .value(iprof as usize, ispin)
            .cloned()
            .ok_or_else(|| "BRepFill_Pipe::Face : null face".to_string())?;
        Face::wrap(the_face)
            .ok_or_else(|| "Standard_TypeMismatch: BRepFill_Pipe::Face".to_string())
    }

    /// `Edge(ESpine, VProfile)` (`cxx:416-453`): the edge swept from an edge of
    /// the spine and a vertex of the profile.
    pub fn edge(&self, e_spine: &Edge, v_profile: &Vertex) -> Result<Edge, String> {
        // `cxx:422-427`.
        let mut count = 0i32;
        let profile = self.profile().cloned().unwrap_or_else(|| TopoShape::new(ShapeType::Shape));
        let iprof = self.find_vertex(&profile, v_profile, &mut count)?;
        if iprof == 0 {
            return Err("Standard_DomainError: BRepFill_Pipe::Edge : Vertex not in the Profile"
                .to_string());
        }

        // `cxx:433-445`.
        let loc = self.location()?;
        let mut ispin = 0usize;
        for ii in 1..=loc.base().nb_law() {
            if ispin != 0 {
                break;
            }
            if let Some(edge) = loc.base().edge(ii) {
                if is_same(&edge.0, &e_spine.0) {
                    ispin = ii;
                }
            }
        }
        if ispin == 0 {
            return Err("Standard_DomainError: BRepFill_Pipe::Edge  : Edge not in the Spine"
                .to_string());
        }

        // `cxx:450-451`.
        let edges = self
            .my_edges
            .as_ref()
            .ok_or_else(|| "BRepFill_Pipe::Edge : the sweep produced no rails".to_string())?;
        let the_edge = edges
            .value(iprof as usize, ispin)
            .cloned()
            .ok_or_else(|| "BRepFill_Pipe::Edge : null edge".to_string())?;
        Edge::wrap(the_edge)
            .ok_or_else(|| "Standard_TypeMismatch: BRepFill_Pipe::Edge".to_string())
    }

    /// `Section(VSpine)` (`cxx:458-493`): the section at a vertex of the spine,
    /// as a compound of the stored section edges.
    pub fn section(&self, v_spine: &Vertex) -> Result<TopoShape, String> {
        // `cxx:470-483`.
        let loc = self.location()?;
        let mut ispin = 0usize;
        for ii in 1..=loc.base().nb_law() + 1 {
            if ispin != 0 {
                break;
            }
            if let Some(vertex) = loc.base().vertex(ii) {
                if is_same(&vertex.0, &v_spine.0) {
                    ispin = ii;
                }
            }
        }
        if ispin == 0 {
            return Err("Standard_DomainError: BRepFill_Pipe::Section  : Vertex not in the Spine"
                .to_string());
        }

        // `cxx:485-491`. The loop bound is `ColLength()` while the *row* index
        // is the loop variable; see the module note on the reference's
        // out-of-range reads.
        let sections = self
            .my_sections
            .as_ref()
            .ok_or_else(|| "BRepFill_Pipe::Section : the sweep produced no sections".to_string())?;
        let builder = TopoBuilder::new();
        let mut comp = Compound::new();
        for ii in 1..=sections.nb_cols() {
            if let Some(section) = sections.value(ii, ispin) {
                builder.add_compound(&mut comp, section);
            }
        }
        Ok(comp.0)
    }

    /// `PipeLine(Point)` (`cxx:502-526`): the wire swept by one point.
    pub fn pipe_line(&mut self, point: &GpPnt) -> Result<Wire, String> {
        // `cxx:503-508`: "Positioning".
        let p = point.transformed(&self.my_trsf);
        let vertex_section = TopoBuilder::new().make_vertex(p, occt_core::precision::CONFUSION);
        let section: Arc<dyn crate::brep_fill_section_law::SectionLaw> =
            Arc::new(crate::brep_fill_shape_law::ShapeLaw::with_vertex(&vertex_section, true)?);

        // `cxx:510-521`: "Sweeping".
        let location = self.location()?;
        let mut mk_sw = crate::brep_fill_sweep::BrepFillSweep::new(section, location, true);
        mk_sw.set_force_approx_c1(self.my_force_approx_c1);
        let _ = make_shape::run_sweep(self, &mut mk_sw);

        // `cxx:522-525`.
        let a_local_shape = mk_sw.shape().cloned();
        self.my_error_on_surf = mk_sw.error_on_surface();
        self.build_history(&mk_sw, &vertex_section.0)?;
        let shape = a_local_shape
            .ok_or_else(|| "StdFail_NotDone: BRepFill_Pipe::PipeLine".to_string())?;
        Wire::wrap(shape)
            .ok_or_else(|| "Standard_TypeMismatch: BRepFill_Pipe::PipeLine".to_string())
    }

    /// The stored location law, or the error a missing one raises.
    fn location(&self) -> Result<Arc<dyn LocationLaw>, String> {
        self.my_loc
            .clone()
            .ok_or_else(|| "BRepFill_Pipe: no location law, Perform did not run".to_string())
    }
}

/// Accessors the split translation files share (`Perform`, `MakeShape`,
/// `BuildHistory`, `PipeLine`). OCCT reaches these fields directly; the port
/// keeps them private to the module tree.
impl Pipe {
    /// `myMode` (`hxx:144`).
    pub(crate) fn mode(&self) -> GeomFillTrihedron {
        self.my_mode
    }

    /// `myTrsf` (`hxx:124`).
    pub(crate) fn trsf(&self) -> &GpTrsf {
        &self.my_trsf
    }

    /// `myTrsf` (`hxx:124`).
    pub(crate) fn set_trsf(&mut self, trsf: GpTrsf) {
        self.my_trsf = trsf;
    }

    /// `myLoc` (`hxx:125`).
    pub(crate) fn loc(&self) -> Option<Arc<dyn LocationLaw>> {
        self.my_loc.clone()
    }

    /// `myLoc` (`hxx:125`).
    pub(crate) fn set_loc(&mut self, loc: Option<Arc<dyn LocationLaw>>) {
        self.my_loc = loc;
    }

    /// `mySpine` (`hxx:121`).
    pub(crate) fn set_spine(&mut self, spine: Option<Wire>) {
        self.my_spine = spine;
    }

    /// `mySpine.IsClosed()` (`cxx:298`, `:599`, `:820`). A missing spine is not
    /// closed.
    pub(crate) fn spine_closed(&self) -> bool {
        self.my_spine.as_ref().map(|wire| wire.is_closed()).unwrap_or(false)
    }

    /// `myProfile` (`hxx:122`).
    pub(crate) fn set_profile(&mut self, profile: Option<TopoShape>) {
        self.my_profile = profile;
    }

    /// `myShape` (`hxx:123`).
    pub(crate) fn set_shape(&mut self, shape: Option<TopoShape>) {
        self.my_shape = shape;
    }

    /// `mySections` (`hxx:126`).
    pub(crate) fn sections(&self) -> Option<&ShapeTable> {
        self.my_sections.as_ref()
    }

    /// `mySections` (`hxx:126`).
    pub(crate) fn set_sections(&mut self, sections: Option<ShapeTable>) {
        self.my_sections = sections;
    }

    /// `myFaces` (`hxx:127`).
    pub(crate) fn faces(&self) -> Option<&ShapeTable> {
        self.my_faces.as_ref()
    }

    /// `myFaces` (`hxx:127`).
    pub(crate) fn set_faces(&mut self, faces: Option<ShapeTable>) {
        self.my_faces = faces;
    }

    /// `myEdges` (`hxx:128`).
    pub(crate) fn edges(&self) -> Option<&ShapeTable> {
        self.my_edges.as_ref()
    }

    /// `myEdges` (`hxx:128`).
    pub(crate) fn set_edges(&mut self, edges: Option<ShapeTable>) {
        self.my_edges = edges;
    }

    /// `myCurIndexOfSectionEdge` (`hxx:136`).
    pub(crate) fn cur_index_of_section_edge(&self) -> usize {
        self.my_cur_index_of_section_edge
    }

    /// `myCurIndexOfSectionEdge` (`hxx:136`).
    pub(crate) fn set_cur_index_of_section_edge(&mut self, index: usize) {
        self.my_cur_index_of_section_edge = index;
    }

    /// `myFirst` (`hxx:137`).
    pub(crate) fn set_first(&mut self, first: Option<TopoShape>) {
        self.my_first = first;
    }

    /// `myLast` (`hxx:138`).
    pub(crate) fn set_last(&mut self, last: Option<TopoShape>) {
        self.my_last = last;
    }

    /// `myGenMap` (`hxx:139-140`).
    pub(crate) fn gen_map(&self) -> &HashMap<ShapeMapKey, Vec<TopoShape>> {
        &self.my_gen_map
    }

    /// `myGenMap` (`hxx:139-140`).
    pub(crate) fn gen_map_mut(&mut self) -> &mut HashMap<ShapeMapKey, Vec<TopoShape>> {
        &mut self.my_gen_map
    }

    /// `mySegmax` (`hxx:142`).
    pub(crate) fn segmax(&self) -> i32 {
        self.my_segmax
    }

    /// `mySegmax` (`hxx:142`).
    pub(crate) fn set_segmax(&mut self, segmax: i32) {
        self.my_segmax = segmax;
    }

    /// `myForceApproxC1` (`hxx:145`).
    pub(crate) fn force_approx_c1(&self) -> bool {
        self.my_force_approx_c1
    }

    /// `myErrorOnSurf` (`hxx:146`).
    pub(crate) fn set_error_on_surf(&mut self, error: f64) {
        self.my_error_on_surf = error;
    }
}
