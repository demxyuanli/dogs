use super::super::prelude::*;
use super::super::*;

/// Reference resolver with memoization; unknown/unsupported entities are
/// recorded as warnings and skipped.
pub(in crate::step) struct Resolver<'a> {
    pub(in crate::step) records: &'a HashMap<usize, Record>,
    /// `STEPControl_ActorRead::myPrecision` for this file
    /// (`STEPControl_ActorRead.cxx:2370-2384`). Reaches the pcurve post-pass as
    /// `StepToTopoDS_TranslateEdgeLoop::Precision()` (`cxx:236`, `cxx:875`).
    pub(in crate::step) precision: f64,
    /// `StepData_Factors::PlaneAngleFactor()` (`StepControl_ActorRead.cxx:2347-2355`,
    /// `StepConstruct_UnitContext.cxx:439-444`): the file's plane angle unit
    /// expressed in radians. Reaches `StepToGeom::MakeConicalSurface`
    /// (`StepToGeom.cxx:1316`) and `GeomConvert_Units::DegreeToRadian`
    /// (`StepToTopoDS_TranslateEdge.cxx:579`).
    pub(in crate::step) plane_angle_factor: f64,
    /// `StepData_Factors::LengthFactor()` (`StepControl_ActorRead.cxx:2347-2355`,
    /// `StepConstruct_UnitContext.cxx:424-428`): the file's length unit expressed
    /// in the cascade unit (millimetres). `1.` for a millimetre file. Reaches
    /// every `StepToGeom::Make*` that turns a STEP coordinate, vector magnitude,
    /// radius or distance into geometry (`StepToGeom.cxx:1179`, `:1222`, `:1315`,
    /// `:1452`, `:1549`, `:1616`, `:1706`, `:1899`, `:1947`, `:2111`, `:2577`) and
    /// `GeomConvert_Units::DegreeToRadian` (`StepToTopoDS_TranslateEdge.cxx:573-580`).
    pub(in crate::step) length_factor: f64,
    pub(in crate::step) b: TopoBuilder,
    pub(in crate::step) shape_cache: RefCell<HashMap<usize, TopoShape>>,
    /// `StepToTopoDS_TranslateTool::Bind` for vertices: a `VERTEX_POINT` record
    /// that another record was bound onto resolves to the bound `TopoDS_Vertex`
    /// (`StepToTopoDS_TranslateEdgeLoop.cxx:384-396`, `:466-477`).
    pub(in crate::step) vertex_bind: RefCell<HashMap<usize, TopoShape>>,
    pub(in crate::step) point_cache: RefCell<HashMap<usize, GpPnt>>,
    pub(in crate::step) dir_cache: RefCell<HashMap<usize, GpDir>>,
    pub(in crate::step) axis_cache: RefCell<HashMap<usize, GpAx2>>,
    pub(in crate::step) curve_cache: RefCell<HashMap<usize, Arc<dyn Curve>>>,
    pub(in crate::step) surface_cache: RefCell<HashMap<usize, Arc<dyn Surface>>>,
    pub(in crate::step) curve2d_cache: RefCell<HashMap<usize, (Arc<dyn Curve2d>, (f64, f64))>>,
    /// SURFACE_CURVE / SEAM_CURVE id → its `associated_geometry` (pcurve_or_surface)
    /// reference list, carried for the face-level pcurve association.
    pub(in crate::step) surface_curve_pcurves: RefCell<HashMap<usize, Vec<usize>>>,
    /// Edge TShape pointer → the curve entity id of its geometry (an EDGE_CURVE's
    /// `curve` attribute, which may be a SURFACE_CURVE carrying pcurves).
    pub(in crate::step) edge_curve_ref: RefCell<HashMap<usize, usize>>,
    pub(in crate::step) resolving: RefCell<HashSet<usize>>,
    pub(in crate::step) warnings: RefCell<Vec<String>>,
    /// `ShapeFix_Shape::myContext` (`ShapeFix_Shape.cxx:75`): the
    /// `ShapeBuild_ReShape` handle shared by every fix tool of one `FixShape`
    /// pass. `FromSTEP.exec.op` is exactly `FixShape`
    /// (`STEPControl_Controller.cxx:201`) and
    /// `ShapeProcess_OperLibrary.cxx:830` runs it once per transferred root,
    /// so one context covers a whole root shape. `ShapeFix_Shape::Perform`
    /// materialises it with `myResult = Context()->Apply(S)`
    /// (`ShapeFix_Shape.cxx:257`, `ShapeFix_Shell.cxx:139`), which is what
    /// rewrites a face whose *neighbor* was split (T0M face 1693 vs 1695).
    /// The stack is pushed per root in `resolve_shell` / `resolve_solid` and
    /// popped together with that `Apply`.
    pub(in crate::step) heal_context: RefCell<Vec<crate::shape_fix_compose_shell::SharedReShape>>,
}
