//! `BRepFill_Sweep` — a swept shell from a generating wire, a section law and
//! a location law (`BRepFill_Sweep.cxx:1-4151`, `.hxx:37-190`).
//!
//! Port status, keeping the OCCT control flow:
//!
//! * Ported (batch A): the data layer, `SetBounds`, `SetTolerance`,
//!   `SetAngularControl`, `SetForceApproxC1`, `CorrectApproxParameters`, the
//!   five accessors, and the static helpers in [`helpers`]
//!   (`NumberOfPoles`, `HasPCurves`, `Translate`, `Couture`, `IsDegen`,
//!   `KeepEdge`, `ReverseEdgeInFirstOrLastWire`, `ReverseModifiedEdges`,
//!   `Oriente`, `BuildVertex`, `NullEdge`).
//! * Ported (batch B): the same-parameter family and the face/edge builders —
//!   [`same_parameter`] (`CheckSameParameter` `cxx:251-283`,
//!   `CheckSameParameterExact` `cxx:288-311`, `SameParameter` `cxx:319-391`,
//!   `CorrectSameParameter` `cxx:394-448`, `UpdateEdgeOnPlane` `cxx:494-517`),
//!   [`build_face`] (`cxx:526-762`), [`build_edge`] (`cxx:766-862` and
//!   `cxx:1476-1655`) and [`substitute`] (`cxx:1341-1367`).
//!   `BRepCheck_Edge::Tolerance` (`BRepCheck_Edge.cxx:598-707`) they need lives
//!   in [`crate::brep_check::edge_tolerance`].
//! * Ported (batch C): [`build_wire`] (`cxx:2078-2209`),
//!   [`update_edge_iso`] (`cxx:1659-1820`), [`build_shell`] (`cxx:2213-3208`),
//!   [`build`] (`cxx:3213-3530`) and the corner helpers of [`corner`] —
//!   `EvalExtrapol` (`cxx:3910-3978`) fully, `PerformCorner`
//!   (`cxx:3587-3906`) up to the `BRepFill_TrimShellCorner` call.
//! * Ported (batch D): [`filling`] — `Filling` (`cxx:867-1338`), the corner
//!   face of the `Round` transition.
//! * Not ported yet, each still absent as a callable method so no call site can
//!   silently take a different branch:
//!   - `BRepFill_TrimShellCorner` (`BRepFill_TrimShellCorner.cxx`, used by
//!     `PerformCorner` at `cxx:3715-3766`).
//!   - `PerformCorner`'s `Round` filling call site (`cxx:3781-3887`).
//!
//! Translation notes:
//! - `GeomAbs_Shape` is [`occt_core::kernel::geomabs::Shape`] and
//!   `GeomFill_ApproxStyle` is [`ApproxStyle`].
//! - `handle(BRepFill_SectionLaw)` / `handle(BRepFill_LocationLaw)` are
//!   `Arc<dyn SectionLaw>` / `Arc<dyn LocationLaw>`; the concrete laws embed
//!   their base data and answer through `base()`.
//! - `NCollection_HArray2<TopoDS_Shape>` is [`helpers::ShapeTable`].
//! - A null `TopoDS_Wire` is `None`; `BRepLib::CheckSameRange` is
//!   [`crate::brep_lib_same_range::check_same_range`].

pub mod build;
pub mod build_edge;
pub mod build_face;
pub mod build_shell;
pub mod build_wire;
pub mod corner;
pub mod filling;
pub mod helpers;
pub mod same_parameter;
pub mod shell_tables;
pub mod substitute;
pub mod trim_shell_corner;
pub mod update_edge_iso;

use std::sync::Arc;

use occt_core::kernel::geomabs::Shape;
use occt_geom::geom_fill::sweep::ApproxStyle;

use crate::brep_fill_location_law::LocationLaw;
use crate::brep_fill_section_law::SectionLaw;
use crate::brep_lib_same_range::check_same_range;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::shape::{TopoShape, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{shape_map_key, ShapeMapKey};

use helpers::ShapeTable;

/// `BRepFill_TransitionStyle` (`BRepFill_TransitionStyle.hxx`): the corner
/// transition between two sweepings.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TransitionStyle {
    /// `BRepFill_Modified`: the corner is the intersection of the two sweeps.
    Modified,
    /// `BRepFill_Right`: the corner keeps the section plane.
    Right,
    /// `BRepFill_Round`: a fillet-like corner.
    Round,
    /// `BRepFill_Transformed`: the transformed corner.
    Transformed,
}

/// `SetTolerance` defaults (`hxx:63-66`).
pub const DEFAULT_BOUND_TOL: f64 = 1.0;
/// `SetTolerance` defaults (`hxx:64`).
pub const DEFAULT_TOL_2D: f64 = 1.0e-5;
/// `SetTolerance` defaults (`hxx:66`).
pub const DEFAULT_TOL_ANGULAR: f64 = 1.0e-2;
/// `SetAngularControl` defaults (`hxx:71`).
pub const DEFAULT_ANG_MIN: f64 = 0.01;
/// `SetAngularControl` defaults (`hxx:71`).
pub const DEFAULT_ANG_MAX: f64 = 6.0;

/// `BRepFill_Sweep` (`BRepFill_Sweep.hxx:44-189`): the class data.
pub struct BrepFillSweep {
    /// `isDone` (`hxx:164`).
    is_done: bool,
    /// `KPart` (`hxx:165`).
    k_part: bool,
    /// `myTol3d` (`hxx:166`).
    my_tol_3d: f64,
    /// `myBoundTol` (`hxx:167`).
    my_bound_tol: f64,
    /// `myTol2d` (`hxx:168`).
    my_tol_2d: f64,
    /// `myTolAngular` (`hxx:169`).
    my_tol_angular: f64,
    /// `myAngMin` (`hxx:170`).
    my_ang_min: f64,
    /// `myAngMax` (`hxx:171`).
    my_ang_max: f64,
    /// `myApproxStyle` (`hxx:172`).
    my_approx_style: ApproxStyle,
    /// `myContinuity` (`hxx:173`).
    my_continuity: Shape,
    /// `myDegmax` (`hxx:174`).
    my_degmax: i32,
    /// `mySegmax` (`hxx:175`).
    my_segmax: i32,
    /// `myForceApproxC1` (`hxx:176`).
    my_force_approx_c1: bool,
    /// `myShape` (`hxx:177`), a null shape is `None`.
    my_shape: Option<TopoShape>,
    /// `myLoc` (`hxx:178`).
    my_loc: Arc<dyn LocationLaw>,
    /// `mySec` (`hxx:179`).
    my_sec: Arc<dyn SectionLaw>,
    /// `myUEdges` (`hxx:180`): the 2-column inter-face table.
    my_u_edges: ShapeTable,
    /// `myVEdges` (`hxx:181`): the section table.
    my_v_edges: ShapeTable,
    /// `myVEdgesModified` (`hxx:182`).
    my_v_edges_modified: std::collections::HashMap<ShapeMapKey, TopoShape>,
    /// `myFaces` (`hxx:183`).
    my_faces: ShapeTable,
    /// `myAuxShape` (`hxx:184`).
    my_aux_shape: Vec<TopoShape>,
    /// `myTapes` (`hxx:185`), 1-based in OCCT and indexed as `index - 1` here.
    my_tapes: Vec<Option<TopoShape>>,
    /// `Error` (`hxx:186`).
    error: f64,
    /// `FirstShape` (`hxx:187`).
    first_shape: Option<Wire>,
    /// `LastShape` (`hxx:188`).
    last_shape: Option<Wire>,
}

impl BrepFillSweep {
    /// `BRepFill_Sweep::BRepFill_Sweep(Section, Location, WithKPart)`
    /// (`cxx:1940-1958`).
    pub fn new(section: Arc<dyn SectionLaw>, location: Arc<dyn LocationLaw>, with_k_part: bool) -> Self {
        let mut this = Self {
            is_done: false,
            k_part: with_k_part,
            my_tol_3d: 0.0,
            my_bound_tol: 0.0,
            my_tol_2d: 0.0,
            my_tol_angular: 0.0,
            my_ang_min: 0.0,
            my_ang_max: 0.0,
            my_approx_style: ApproxStyle::Location,
            my_continuity: Shape::C2,
            my_degmax: 11,
            my_segmax: 30,
            my_force_approx_c1: false,
            my_shape: None,
            my_loc: location,
            my_sec: section,
            my_u_edges: ShapeTable::default(),
            my_v_edges: ShapeTable::default(),
            my_v_edges_modified: std::collections::HashMap::new(),
            my_faces: ShapeTable::default(),
            my_aux_shape: Vec::new(),
            my_tapes: Vec::new(),
            error: 0.0,
            first_shape: None,
            last_shape: None,
        };
        this.set_tolerance(1.0e-4, DEFAULT_BOUND_TOL, DEFAULT_TOL_2D, DEFAULT_TOL_ANGULAR);
        this.set_angular_control(DEFAULT_ANG_MIN, DEFAULT_ANG_MAX);
        this.my_aux_shape.clear();
        this.my_approx_style = ApproxStyle::Location;
        this.my_continuity = Shape::C2;
        this.my_degmax = 11;
        this.my_segmax = 30;
        this.my_force_approx_c1 = false;
        this
    }

    /// `BRepFill_Sweep::SetBounds` (`cxx:1962-2007`): the two limiting wires,
    /// with the `SameRange` flags cleared on the edges that fail
    /// `BRepLib::CheckSameRange` (PRO13551).
    pub fn set_bounds(&mut self, first: Option<Wire>, last: Option<Wire>) {
        self.first_shape = first;
        self.last_shape = last;
        let reg = GeometryRegistry::global();
        for w in [self.first_shape.as_ref(), self.last_shape.as_ref()].into_iter().flatten() {
            let mut wexp = WireExplorer::of_wire(w);
            while wexp.more() {
                if let Some(e) = wexp.current_edge() {
                    if !check_same_range(&e.0, 1.0e-12) {
                        reg.set_same_range(&e.0, false);
                        reg.set_same_parameter(&e.0, false);
                    }
                }
                wexp.next();
            }
        }
    }

    /// `BRepFill_Sweep::SetTolerance` (`cxx:2011-2021`).
    pub fn set_tolerance(&mut self, tol3d: f64, bound_tol: f64, tol2d: f64, tol_angular: f64) {
        self.my_tol_3d = tol3d;
        self.my_bound_tol = bound_tol;
        self.my_tol_2d = tol2d;
        self.my_tol_angular = tol_angular;
    }

    /// `BRepFill_Sweep::SetTolerance(Tol3d, BoundTol = 1., Tol2d = 1.e-5,
    /// TolAngular = 1.e-2)` (`hxx:63-66`) through its one-argument call site
    /// (`BRepFill_Draft.cxx:470`): only `Tol3d` is replaced, the other three
    /// keep the constructor's values.
    pub fn set_tolerance_3d(&mut self, tol3d: f64) {
        self.my_tol_3d = tol3d;
    }

    /// `BRepFill_Sweep::SetAngularControl` (`cxx:2024-2028`).
    pub fn set_angular_control(&mut self, min_angle: f64, max_angle: f64) {
        self.my_ang_min = min_angle.max(occt_core::precision::ANGULAR);
        self.my_ang_max = max_angle.min(6.28);
    }

    /// `BRepFill_Sweep::SetForceApproxC1` (`cxx:2036-2039`).
    pub fn set_force_approx_c1(&mut self, force_approx_c1: bool) {
        self.my_force_approx_c1 = force_approx_c1;
    }

    /// `BRepFill_Sweep::CorrectApproxParameters` (`cxx:2043-2075`): lower the
    /// requested continuity to the path's worst edge and raise `Segmax` to its
    /// widest edge, answering `true` when either value changed.
    pub fn correct_approx_parameters(&mut self) -> Result<bool, String> {
        let path = self.my_loc.base().wire().clone();
        let mut new_cont = self.my_continuity;
        let mut new_segmax = self.my_segmax;
        for e in crate::topo_tools::edges_of(path.shape()) {
            let Some(c) = crate::brep_tool::BRepTool::edge_curve(&e) else {
                return Err("BRepAdaptor_Curve : null curve".to_string());
            };
            let a_continuity = shape_from_continuity(c.continuity());
            let a_nb_interv = c.nb_intervals(Shape::CN as u8);
            if a_continuity < new_cont {
                new_cont = a_continuity;
            }
            if a_nb_interv > new_segmax {
                new_segmax = a_nb_interv;
            }
        }
        let corrected = new_cont != self.my_continuity || new_segmax != self.my_segmax;
        self.my_continuity = new_cont;
        self.my_segmax = new_segmax;
        Ok(corrected)
    }

    /// `BRepFill_Sweep::IsDone` (`cxx:3534-3537`).
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// `BRepFill_Sweep::Shape` (`cxx:3541-3544`): a null shape is `None`.
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_shape.as_ref()
    }

    /// `BRepFill_Sweep::ErrorOnSurface` (`cxx:3548-3551`).
    pub fn error_on_surface(&self) -> f64 {
        self.error
    }

    /// `BRepFill_Sweep::SubShape` (`cxx:3555-3558`): `myFaces`.
    pub fn sub_shape(&self) -> &ShapeTable {
        &self.my_faces
    }

    /// `BRepFill_Sweep::InterFaces` (`cxx:3562-3565`): `myUEdges`.
    pub fn inter_faces(&self) -> &ShapeTable {
        &self.my_u_edges
    }

    /// `BRepFill_Sweep::Sections` (`cxx:3571-3574`): `myVEdges`.
    pub fn sections(&self) -> &ShapeTable {
        &self.my_v_edges
    }

    /// `BRepFill_Sweep::Tape(Index)` (`cxx:3580-3583`): 1-based index.
    pub fn tape(&self, index: usize) -> Option<&TopoShape> {
        if index == 0 {
            return None;
        }
        self.my_tapes.get(index - 1).and_then(|t| t.as_ref())
    }

    /// The `ReversedEdges` map a future `Build` fills (`hxx:93`): keyed by the
    /// pre-modification edge, valued with its reversed copy.
    pub fn reversed_edges(&self) -> &std::collections::HashMap<ShapeMapKey, TopoShape> {
        &self.my_v_edges_modified
    }

    /// `KPart` (`hxx:165`).
    pub fn k_part(&self) -> bool {
        self.k_part
    }
}

/// `GeomAbs_Shape` from a `Continuity()` code (`GeomAbs_Shape.hxx`).
fn shape_from_continuity(c: u8) -> Shape {
    match c {
        0 => Shape::C0,
        1 => Shape::G1,
        2 => Shape::C1,
        3 => Shape::G2,
        4 => Shape::C2,
        5 => Shape::C3,
        _ => Shape::CN,
    }
}

/// `TopoDS_Shape` hash key used by the `myVEdgesModified` data map.
pub fn edge_map_key(s: &TopoShape) -> ShapeMapKey {
    shape_map_key(s)
}
