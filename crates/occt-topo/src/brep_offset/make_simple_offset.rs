//! `BRepOffset_MakeSimpleOffset` (`BRepOffset_MakeSimpleOffset.cxx:1-703`,
//! `BRepOffset_MakeSimpleOffset.hxx:34-177`): the simple offset algorithm.
//!
//! The algorithm is the one the header describes (`hxx:41-51`): build the
//! source -> image maps through [`BRepOffsetSimpleOffset`] (a
//! [`Modification`]), let [`Modifier`] apply it, then restore the topological
//! integrity of the answer. [`BRepOffsetMakeSimpleOffset::perform`] is
//! `cxx:155-207`:
//!
//! * `Clear()` (`cxx:158`), the null-input test (`cxx:161-165`), the lazily
//!   cached [`tgtfaces`] angle (`cxx:167-170`).
//! * `myBuilder.Init` + `new BRepOffset_SimpleOffset` + `myBuilder.Perform`
//!   (`cxx:172-178`), `ErrorOffsetComputation` when the modifier fails
//!   (`cxx:178-182`).
//! * `ModifiedShape(myInputShape)` (`cxx:184`) and the degeneracy fix
//!   (`cxx:186-200`): an edge the input marks as degenerated stays degenerated
//!   on its image.
//! * `BuildMissingWalls` when the solid flag is set (`cxx:202-206`).
//!
//! Port notes:
//!
//! * `BRepOffsetSimple_Status` (`hxx:30-38`) is
//!   [`BRepOffsetSimpleStatus`]; `GetErrorMessage` (`cxx:86-117`) answers the
//!   same five strings, empty for `OK`.
//! * `TopExp::MapShapesAndAncestors` (`TopExp.cxx:80-121`) of
//!   `ComputeMaxAngleOnShape` (`cxx:331-332`) is
//!   [`crate::topo_tools::AncestorMap`]; indices are 0-based.
//! * `BRepAdaptor_Surface(F, false)->D1(u, v)` (`cxx:235-238`, `:286`) is
//!   [`crate::brep_tool::BRepTool::face_surface`] + `Surface::d1`.
//! * `BRepAdaptor_Curve2d(E, F)->D0(u)` (`cxx:253-258`, `:284`, `:300`) is the
//!   pcurve of [`crate::boptools_2d::curve_on_surface_range`] evaluated at `u`.
//!   `E.Orientation(FORWARD / REVERSED)` (`cxx:252`, `:256-258`) does not
//!   change that lookup: `BRep_Tool::CurveOnSurface` ignores the occurrence
//!   orientation, exactly as the port's helper does.
//! * `BRep_Tool::Range(E, f, l)` (`cxx:264`) is
//!   [`crate::brep_tool::BRepTool::edge_parameters`].
//! * `gp_Vec::Crossed` / `Magnitude` / `Reverse` / `Angle` (`cxx:285-318`) are
//!   the same-named `GpVec` methods.
//! * `BRep_Tool::Continuity(E, F1, F2) <= GeomAbs_C0` (`cxx:376`) is
//!   [`crate::brep_lib_continuity::stored_continuity`] answering `0`
//!   (`GeomAbs_C0`); an edge that was never encoded reads `C0`, as in OCCT.
//! * The `try` / `catch (Standard_Failure)` around `tgtfaces`
//!   (`cxx:377-386`) has no counterpart: `Surface::d1` and `Curve2d::value`
//!   answer values instead of throwing, and a degenerate normal is the
//!   `continue` of `cxx:290-293` / `:309-312`.
//! * `BRep_Tool::MaxTolerance` of `GetSafeOffset` (`cxx:144`) is
//!   [`crate::brep_tool::BRepTool::max_tolerance`].
//! * `myReShape` (`hxx:163`) is a [`SharedReShape`], so the `const`
//!   `Generated` / `Modified` (`cxx:671-703`) can still run `Apply` / `Status`.
//! * `BuildMissingWalls` (`cxx:404-509`) starts from
//!   `ShapeAnalysis_FreeBounds(myInputShape)`, which is the three-argument
//!   constructor with its defaults (`hxx:106-110`: `splitclosed = false`,
//!   `splitopen = true`, `checkinternaledges = false`); the port calls
//!   [`ShapeAnalysisFreeBounds::of_shells`] with exactly those arguments.
//!   `GetClosedWires()` is [`ShapeAnalysisFreeBounds::closed_wires`].
//! * `ShapeFix_Edge::SetContext(myReShape)` (`cxx:440`) has no effect on the
//!   following `FixSameParameter(edge)`: the one-argument overload
//!   (`ShapeFix_Edge.cxx:790-794`) forwards to the three-argument one, which
//!   never reads `myContext`. The call is [`fix_same_parameter_no_face`].
//! * `BRepTools_Quilt` (`cxx:479-481`) is [`crate::brep_sweep::Quilt`];
//!   `BRep_Tool::IsClosed(shell)` (`cxx:496`) is
//!   [`crate::brep_tool::BRepTool::is_closed`].
//! * `BuildWallFace` (`cxx:513-667`): `BRepLib_MakeEdge(V1, V2)`
//!   (`cxx:547`, `:568`) is [`make_edge_between`], which fails like
//!   `BRepLib_LineThroughIdenticPoints`; `BRepLib::BuildCurves3d(wire)`
//!   (`cxx:592`) is [`build_curves3d`]; `BRepLib_MakeFace(wire, true)`
//!   (`cxx:599`) is [`make_face_from_wire`] with `only_plane`, and its
//!   `try` / `catch` arm is the `None` both OCCT paths answer with;
//!   `GeomFill_Generator` (`cxx:620-624`) is
//!   [`generator_surface`]; `theSurf->UIso(u)` (`cxx:641`, `:655`, `:658`) is
//!   [`Surface::u_iso_curve`], absent on a surface that has no iso curve, in
//!   which case OCCT would store a null 3D curve and the port keeps the
//!   previous one; `BRepLib_MakeFace(theSurf, wire)` (`cxx:663`) is
//!   [`make_face_from_surface_and_wire`]; the `BRep_Builder` writes are
//!   [`GeometryRegistry`] (see `make_offset/missing_walls/perform.rs`, which maps the
//!   same family of calls).
//! * `TopoDS_Shape::Oriented(o)` (`cxx:519`, `:534-535`, `:547-548`,
//!   `:563`, `:568-569`) *replaces* the orientation
//!   (`TopoDS_Shape.hxx`), which is what [`TopoShape::oriented`] does, so no
//!   compose helper is used here. `Orientation(o)` (`cxx:519`, `:529`) is
//!   [`TopoShape::set_orientation`].
//! * UNPORTED: nothing on the `BuildMissingWalls` / `BuildWallFace` path. The
//!   reference dereferences a null shell at `cxx:496` (there would be no shell
//!   if every face failed to sew); the port answers
//!   `BRepOffsetSimple_ErrorNonClosedShell` instead of throwing
//!   `Standard_NullObject`.

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::gp::{GpDir2d, GpPnt2d};
use occt_core::precision::CONFUSION;
use occt_geom::curve::Curve;
use occt_geom::geom_fill::generator::generator_surface;
use occt_geom::surface::Surface;
use occt_geom::trimmed::GeomTrimmedCurve;
use occt_geom2d::line::Geom2dLine;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_lib_build_curves3d::build_curves3d;
use crate::brep_lib_continuity::stored_continuity;
use crate::brep_lib_make_face::{make_face_from_surface_and_wire, make_face_from_wire};
use crate::brep_offset::simple_offset::BRepOffsetSimpleOffset;
use crate::brep_sweep::Quilt;
use crate::brep_tool::BRepTool;
use crate::brep_tools_modifier::{Modification, Modifier};
use crate::builder::TopoBuilder;
use crate::shape::{Compound, Edge, Face, TopoShape, Vertex, Wire};
use crate::shape_fix_compose_shell::{ReShape, SharedReShape};
use crate::shhealing::{fix_same_parameter_no_face, ShapeAnalysisFreeBounds};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{
    edge_vertices, explorer_shapes, is_same, shape_map_key, AncestorMap, ShapeMapKey,
};

/// `BRepOffsetSimple_Status` (`BRepOffset_MakeSimpleOffset.hxx:30-38`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BRepOffsetSimpleStatus {
    /// `BRepOffsetSimple_OK`.
    Ok,
    /// `BRepOffsetSimple_NullInputShape`.
    NullInputShape,
    /// `BRepOffsetSimple_ErrorOffsetComputation`.
    ErrorOffsetComputation,
    /// `BRepOffsetSimple_ErrorWallFaceComputation`.
    ErrorWallFaceComputation,
    /// `BRepOffsetSimple_ErrorInvalidNbShells`.
    ErrorInvalidNbShells,
    /// `BRepOffsetSimple_ErrorNonClosedShell`.
    ErrorNonClosedShell,
}

/// `tgtfaces(Ed, F1, F2, couture, theResAngle)` (`cxx:215-322`): the angle at
/// the border between two squares that share a front edge.
///
/// The surface normals are sampled at `23 + 1` parameters of the edge range,
/// shrunk by `1/100` at both ends (`cxx:266-268`), and the maximum angle is
/// accumulated into `the_res_angle` (`cxx:318-320`). A null pcurve on either
/// face answers immediately (`cxx:225-233`); a degenerate normal (magnitude at
/// most `1e-12`) skips that sample (`cxx:290-293`, `:309-312`).
fn tgtfaces(edge: &Edge, f1: &Face, f2: &Face, couture: bool, the_res_angle: &mut f64) {
    // `cxx:225-233`: the pcurves have to exist on both faces.
    let Some((c2d1, _, _)) = curve_on_surface_range(edge, f1) else {
        return;
    };
    let Some((c2d2, _, _)) = curve_on_surface_range(edge, f2) else {
        return;
    };

    // `cxx:235-249`: the two surface adaptors. `couture` makes both read the
    // same face, hence the same surface.
    let Some(surf1) = BRepTool::face_surface(f1) else {
        return;
    };
    let surf2: Arc<dyn Surface> = if couture {
        surf1.clone()
    } else {
        match BRepTool::face_surface(f2) {
            Some(s) => s,
            None => return,
        }
    };

    // `cxx:260-261`.
    let rev1 = f1.0.orientation() == Orientation::Reversed;
    let rev2 = f2.0.orientation() == Orientation::Reversed;

    // `cxx:262-268`: `BRep_Tool::Range(E, f, l)` shrunk by `(l - f) / 100`.
    let (a_first, a_last) = BRepTool::edge_parameters(edge);
    let eps = (a_last - a_first) / 100.0;
    let first = a_first + eps;
    let last = a_last - eps;

    const NB_PNT: i32 = 23;
    for i in 0..=NB_PNT {
        // `cxx:280`: the parameter of this sample.
        let u = first + (last - first) * f64::from(i) / f64::from(NB_PNT);

        // `cxx:283-286`: the normal of the first face at the pcurve point.
        let p: GpPnt2d = c2d1.value(u);
        let (_, du1, dv1) = surf1.d1(p.x(), p.y());
        let mut d1 = du1.crossed(&dv1);
        let norm = d1.magnitude();
        if norm > 1.0e-12 {
            d1 = d1.divided(norm);
        } else {
            continue; // `cxx:290-293`: skip degenerated point.
        }
        // `cxx:294-297`.
        if rev1 {
            d1.reverse();
        }

        // `cxx:300-303`.
        let p: GpPnt2d = c2d2.value(u);
        let (_, du2, dv2) = surf2.d1(p.x(), p.y());
        let mut d2 = du2.crossed(&dv2);
        let norm = d2.magnitude();
        if norm > 1.0e-12 {
            d2 = d2.divided(norm);
        } else {
            continue; // `cxx:309-312`: skip degenerated point.
        }
        // `cxx:313-316`.
        if rev2 {
            d2.reverse();
        }

        // `cxx:318-320`.
        *the_res_angle = the_res_angle.max(d1.angle(&d2));
    }
}

/// `ComputeMaxAngleOnShape(S, theResAngle)` (`cxx:329-388`): the maximum angle
/// over every edge whose two ancestor faces meet with continuity at most
/// `GeomAbs_C0`.
fn compute_max_angle_on_shape(shape: &TopoShape, the_res_angle: &mut f64) {
    // `cxx:331-332`.
    let edge_face_map =
        AncestorMap::map_shapes_and_ancestors(shape, ShapeType::Edge, ShapeType::Face);

    for idx in 0..edge_face_map.length() {
        // `cxx:340`.
        let Some(edge) = Edge::wrap(edge_face_map.find_key(idx).clone()) else {
            continue;
        };
        let mut found = false;
        let mut couture = false;
        // `cxx:343`: `F1.Nullify()`.
        let mut f1: Option<Face> = None;
        let mut f2: Option<Face> = None;

        // `cxx:345-357`: the first two distinct ancestor faces.
        for anc in edge_face_map.ancestors(idx) {
            if found {
                break; // `cxx:345`: `It.More() && !found`.
            }
            let Some(face) = Face::wrap(anc.clone()) else {
                continue;
            };
            match f1.as_ref() {
                None => f1 = Some(face), // `cxx:347-349`.
                Some(cur) => {
                    // `cxx:352-356`.
                    if !is_same(&cur.0, &face.0) {
                        found = true;
                        f2 = Some(face);
                    }
                }
            }
        }

        // `cxx:359-372`: a sewing edge, i.e. the same face on both sides.
        if !found {
            if let Some(cur_f1) = f1.as_ref() {
                let or_e = edge.0.orientation();
                for cur_e in explorer_shapes(&cur_f1.0, ShapeType::Edge) {
                    if found {
                        break; // `cxx:360`: `Ex.More() && !found`.
                    }
                    if is_same(&edge.0, &cur_e) && or_e != cur_e.orientation() {
                        found = true;
                        couture = true;
                        f2 = Some(cur_f1.clone());
                    }
                }
            }
        }

        // `cxx:375-387`.
        if found {
            let (Some(cur_f1), Some(cur_f2)) = (f1.as_ref(), f2.as_ref()) else {
                continue;
            };
            if stored_continuity(&edge, cur_f1, cur_f2) == 0 {
                tgtfaces(&edge, cur_f1, cur_f2, couture, the_res_angle);
            }
        }
    }
}

/// `BRepLib_MakeEdge(V1, V2)` (`BRepLib_MakeEdge.cxx:185-198`): the segment
/// between two existing vertices. Fails with
/// `BRepLib_LineThroughIdenticPoints` when the two points are within
/// `gp::Resolution()`.
fn make_edge_between(vertex1: &Vertex, vertex2: &Vertex) -> Option<Edge> {
    let p1 = BRepTool::vertex_point(vertex1);
    let p2 = BRepTool::vertex_point(vertex2);
    if p1.distance(&p2) <= occt_core::precision::REAL_SMALL {
        return None;
    }
    Some(TopoBuilder::new().make_edge_segment_with_vertices(&p1, &p2, vertex1, vertex2))
}

/// `BRepOffset_MakeSimpleOffset` (`BRepOffset_MakeSimpleOffset.hxx:61-171`).
pub struct BRepOffsetMakeSimpleOffset {
    /// `myInputShape`.
    input_shape: Option<TopoShape>,
    /// `myOffsetValue`.
    offset_value: f64,
    /// `myTolerance` (for singularities).
    tolerance: f64,
    /// `myIsBuildSolid`.
    is_build_solid: bool,
    /// `myMaxAngle`.
    max_angle: f64,
    /// `myError`.
    error: BRepOffsetSimpleStatus,
    /// `myIsDone`.
    is_done: bool,
    /// `myMapVE`: vertex -> wall edge, shared by the adjacent wall faces
    /// (`NCollection_DataMap<TopoDS_Vertex, TopoDS_Edge>`, `hxx:158`).
    map_ve: HashMap<ShapeMapKey, Edge>,
    /// `myBuilder`.
    builder: Modifier,
    /// `myReShape`.
    reshape: SharedReShape,
    /// `myResShape`.
    res_shape: Option<TopoShape>,
}

impl Default for BRepOffsetMakeSimpleOffset {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepOffsetMakeSimpleOffset {
    /// `BRepOffset_MakeSimpleOffset()` (`cxx:47-57`).
    pub fn new() -> Self {
        Self {
            input_shape: None,
            offset_value: 0.0,
            tolerance: CONFUSION,
            is_build_solid: false,
            max_angle: 0.0,
            error: BRepOffsetSimpleStatus::Ok,
            is_done: false,
            map_ve: HashMap::new(),
            builder: Modifier::new(false),
            reshape: SharedReShape::new(),
            res_shape: None,
        }
    }

    /// `BRepOffset_MakeSimpleOffset(theInputShape, theOffsetValue)`
    /// (`cxx:61-72`).
    pub fn with_shape(input_shape: &TopoShape, offset_value: f64) -> Self {
        Self {
            input_shape: Some(input_shape.clone()),
            offset_value,
            ..Self::new()
        }
    }

    /// `Initialize(theInputShape, theOffsetValue)` (`cxx:76-82`).
    pub fn initialize(&mut self, input_shape: &TopoShape, offset_value: f64) {
        self.input_shape = Some(input_shape.clone()); // `cxx:79`.
        self.offset_value = offset_value; // `cxx:80`.
        self.clear(); // `cxx:81`.
    }

    /// `GetErrorMessage()` (`cxx:86-117`): the message of the current error,
    /// empty for `BRepOffsetSimple_OK`.
    pub fn error_message(&self) -> &'static str {
        match self.error {
            BRepOffsetSimpleStatus::NullInputShape => "Null input shape",
            BRepOffsetSimpleStatus::ErrorOffsetComputation => "Error during offset construction",
            BRepOffsetSimpleStatus::ErrorWallFaceComputation => "Error during building wall face",
            BRepOffsetSimpleStatus::ErrorInvalidNbShells => "Result contains two or more shells",
            BRepOffsetSimpleStatus::ErrorNonClosedShell => "Result shell is not closed",
            BRepOffsetSimpleStatus::Ok => "",
        }
    }

    /// `GetError()` (`hxx:81`).
    pub fn error(&self) -> BRepOffsetSimpleStatus {
        self.error
    }

    /// `GetBuildSolidFlag()` (`hxx:85`).
    pub fn build_solid_flag(&self) -> bool {
        self.is_build_solid
    }

    /// `SetBuildSolidFlag(theBuildFlag)` (`hxx:88`).
    pub fn set_build_solid_flag(&mut self, the_build_flag: bool) {
        self.is_build_solid = the_build_flag;
    }

    /// `GetOffsetValue()` (`hxx:91`).
    pub fn offset_value(&self) -> f64 {
        self.offset_value
    }

    /// `SetOffsetValue(theOffsetValue)` (`hxx:94`).
    pub fn set_offset_value(&mut self, the_offset_value: f64) {
        self.offset_value = the_offset_value;
    }

    /// `GetTolerance()` (`hxx:97`).
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// `SetTolerance(theValue)` (`hxx:101`).
    pub fn set_tolerance(&mut self, the_value: f64) {
        self.tolerance = the_value;
    }

    /// `IsDone()` (`hxx:104`).
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// `GetResultShape()` (`hxx:107`).
    pub fn result_shape(&self) -> Option<&TopoShape> {
        self.res_shape.as_ref()
    }

    /// `Clear()` (`cxx:121-128`): drops the previous result.
    pub fn clear(&mut self) {
        self.is_done = false; // `cxx:123`.
        self.error = BRepOffsetSimpleStatus::Ok; // `cxx:124`.
        self.max_angle = 0.0; // `cxx:125`.
        self.map_ve.clear(); // `cxx:126`.
        // `cxx:127`: `myReShape->Clear()` drops the stored modifications.
        self.reshape = SharedReShape::new();
    }

    /// `GetSafeOffset(theExpectedToler)` (`cxx:132-150`): the offset distance
    /// at which the tolerance growth of a face junction stays inside
    /// `theExpectedToler`.
    pub fn get_safe_offset(&mut self, the_expected_toler: f64) -> f64 {
        // `cxx:134-137`.
        let Some(input) = self.input_shape.clone() else {
            return 0.0; // Input shape is null.
        };

        // `cxx:139-143`: non-initialized.
        if self.max_angle == 0.0 {
            self.compute_max_angle();
        }

        // `cxx:145-146`.
        let a_max_tol = BRepTool::max_tolerance(&input, ShapeType::Vertex);

        // `cxx:148-149`: the minimal distance cannot be lower than `0.0`.
        ((the_expected_toler - a_max_tol) / (2.0 * self.max_angle)).max(0.0)
    }

    /// `ComputeMaxAngle()` (`cxx:393-397`).
    pub fn compute_max_angle(&mut self) {
        let Some(input) = self.input_shape.clone() else {
            return;
        };
        compute_max_angle_on_shape(&input, &mut self.max_angle);
    }

    /// `Perform()` (`cxx:155-207`).
    pub fn perform(&mut self) {
        // `cxx:158`: clear the result of previous computations.
        self.clear();

        // `cxx:161-165`.
        let Some(input) = self.input_shape.clone() else {
            self.error = BRepOffsetSimpleStatus::NullInputShape;
            return;
        };

        // `cxx:167-170`: non-initialized.
        if self.max_angle == 0.0 {
            self.compute_max_angle();
        }

        // `cxx:172-178`.
        self.builder.init(&input);
        let mapper = BRepOffsetSimpleOffset::new(&input, self.offset_value, self.tolerance);
        let _ = self.builder.perform(&mapper as &dyn Modification);

        // `cxx:178-182`.
        if !self.builder.is_done() {
            self.error = BRepOffsetSimpleStatus::ErrorOffsetComputation;
            return;
        }

        // `cxx:184`.
        self.res_shape = self.builder.modified_shape(&input).ok();

        // `cxx:186-200`: fix degeneracy, a degenerated edge has to map to a
        // degenerated one.
        let reg = GeometryRegistry::global();
        for edge_shape in explorer_shapes(&input, ShapeType::Edge) {
            // `cxx:192-195`.
            if !reg.is_degenerated_edge(&edge_shape) {
                continue;
            }
            // `cxx:197-198`.
            let Ok(image) = self.builder.modified_shape(&edge_shape) else {
                continue;
            };
            // `cxx:199`: `aBB.Degenerated(anEdge, true)`.
            reg.set_degenerated(&image, true);
        }

        // `cxx:202-206`: a failed wall build leaves `myIsDone` unset.
        if self.is_build_solid && !self.build_missing_walls() {
            return;
        }

        // `cxx:207`.
        self.is_done = true;
    }

    /// `BuildMissingWalls()` (`cxx:403-509`): closes the offset surface with
    /// a wall on every outer edge and stitches the lot into a solid.
    fn build_missing_walls(&mut self) -> bool {
        // `cxx:406-408`: the compound of wall faces.
        let builder = TopoBuilder::new();
        let mut new_faces = Compound::new();

        // `cxx:411-412`: `ShapeAnalysis_FreeBounds(myInputShape)` is the
        // three-argument constructor with its defaults (`hxx:106-110`:
        // `splitclosed = false`, `splitopen = true`, `checkinternaledges =
        // false`); its `GetClosedWires()` answers the outer contours.
        let Some(input) = self.input_shape.clone() else {
            return false;
        };
        let free_bounds = ShapeAnalysisFreeBounds::of_shells(&input, false, true, false);
        let free_wires = free_bounds.closed_wires().clone();

        // `cxx:415-436`: one wall face per outer edge, in explorer order.
        for wire_shape in explorer_shapes(&free_wires.0, ShapeType::Wire) {
            let Some(cur_wire) = Wire::wrap(wire_shape) else {
                continue;
            };
            for edge_shape in explorer_shapes(&cur_wire.0, ShapeType::Edge) {
                let Some(cur_edge) = Edge::wrap(edge_shape) else {
                    continue;
                };
                // `cxx:426-432`.
                let Some(new_face) = self.build_wall_face(&cur_edge) else {
                    self.error = BRepOffsetSimpleStatus::ErrorWallFaceComputation;
                    return false;
                };
                builder.add_compound(&mut new_faces, &new_face.0); // `cxx:434`.
            }
        }

        // `cxx:439-447`: `ShapeFix_Edge::FixSameParameter` on every wall edge.
        // `SetContext(myReShape)` (`cxx:440`) has no effect on this call: the
        // one-argument overload (`ShapeFix_Edge.cxx:790-794`) forwards straight
        // to the three-argument one, which never reads `myContext`.
        for edge_shape in explorer_shapes(&new_faces.0, ShapeType::Edge) {
            if let Some(edge) = Edge::wrap(edge_shape) {
                fix_same_parameter_no_face(&edge);
            }
        }

        // `cxx:450-473`.
        let mut res_compound = Compound::new();
        for face_shape in explorer_shapes(&input, ShapeType::Face) {
            builder.add_compound(&mut res_compound, &face_shape); // `cxx:454-458`.
        }
        let Some(res_shape) = self.res_shape.clone() else {
            return false;
        };
        for face_shape in explorer_shapes(&res_shape, ShapeType::Face) {
            builder.add_compound(&mut res_compound, &face_shape); // `cxx:461-465`.
        }
        for face_shape in explorer_shapes(&new_faces.0, ShapeType::Face) {
            builder.add_compound(&mut res_compound, &face_shape); // `cxx:468-473`.
        }

        // `cxx:475-476`: `myReShape->Apply(aResCompound)`. `Apply` mutates the
        // stored modification state, hence the clone - exactly what
        // `Generated` (`cxx:684`) does.
        let mut reshape = self.reshape.clone();
        let applied = reshape.apply(&res_compound.0);

        // `cxx:478-481`: `BRepTools_Quilt` sews the shared wall edges.
        let mut quilt = Quilt::new();
        quilt.add(&applied);
        let shells = quilt.shells();

        // `cxx:483-494`: exactly one shell, otherwise
        // `BRepOffsetSimple_ErrorInvalidNbShells`.
        let mut res_shell: Option<crate::shape::Shell> = None;
        for shell_shape in explorer_shapes(&shells, ShapeType::Shell) {
            if res_shell.is_some() {
                self.error = BRepOffsetSimpleStatus::ErrorInvalidNbShells;
                return false;
            }
            res_shell = crate::shape::Shell::wrap(shell_shape);
        }
        let Some(res_shell) = res_shell else {
            // `cxx:496` would run `BRep_Tool::IsClosed` on a null shell, which
            // throws `Standard_NullObject`; the port reports the same class of
            // failure instead of dereferencing it.
            self.error = BRepOffsetSimpleStatus::ErrorNonClosedShell;
            return false;
        };

        // `cxx:496-500`.
        if !BRepTool::is_closed(&res_shell.0) {
            self.error = BRepOffsetSimpleStatus::ErrorNonClosedShell;
            return false;
        }

        // `cxx:502-506`: `MakeSolid` + `Add(solid, shell)`.
        let res_solid = builder.make_solid(&[res_shell]);
        self.res_shape = Some(res_solid.0);

        true
    }

    /// `BuildWallFace(theOrigEdge)` (`cxx:513-667`): the ruled face between an
    /// outer edge and its offset image, closed by the two connecting edges that
    /// adjacent wall faces share through [`Self::map_ve`].
    fn build_wall_face(&mut self, the_orig_edge: &Edge) -> Option<Face> {
        // `cxx:518-519`: the offset image, reversed so the wire reads as
        // orig (forward) -> wall1 -> image (reversed) -> wall2.
        let mut new_edge = Edge::wrap(self.builder.modified_shape(&the_orig_edge.0).ok()?)?;
        new_edge.0.set_orientation(Orientation::Reversed);

        // `cxx:521-522`: `TopExp::Vertices(aNewEdge, aNewV1, aNewV2)`.
        let (new_v1, new_v2) = edge_vertices(&new_edge);

        // `cxx:528-531`: a FORWARD copy of the original, so that its vertices
        // come out in the stored order, and their `TopExp::Vertices` pair.
        let orig_copy = Edge(the_orig_edge.0.oriented(Orientation::Forward));
        let (v1, v2) = edge_vertices(&orig_copy);
        let (v1, v2) = (v1?, v2?);

        // `cxx:534-535`: the FORWARD-oriented vertices are the `myMapVE` keys
        // (`TopTools_ShapeMapHasher` ignores the orientation anyway).
        let forward_v1 = v1.0.oriented(Orientation::Forward);
        let forward_v2 = v2.0.oriented(Orientation::Forward);

        // `cxx:538-556`: the wall edge at the original edge's last vertex.
        let wall1 = match self.map_ve.get(&shape_map_key(&forward_v2)) {
            // `cxx:541-542`: the stored edge is already FORWARD.
            Some(wall) => wall.clone(),
            None => {
                // `cxx:547-548`: `BRepLib_MakeEdge(V2, aNewV2.Reversed())`.
                let start = Vertex(v2.0.oriented(Orientation::Forward));
                let end = Vertex(new_v2?.0.oriented(Orientation::Reversed));
                let wall = make_edge_between(&start, &end)?;
                self.map_ve
                    .insert(shape_map_key(&forward_v2), wall.clone()); // `cxx:555`.
                wall
            }
        };

        // `cxx:559-580`: the wall edge at the original edge's first vertex.
        let wall2 = match self.map_ve.get(&shape_map_key(&forward_v1)) {
            // `cxx:563`: retrieved already reversed.
            Some(wall) => Edge(wall.0.oriented(Orientation::Reversed)),
            None => {
                // `cxx:568-569`: `BRepLib_MakeEdge(V1, aNewV1.Reversed())`.
                let start = Vertex(v1.0.oriented(Orientation::Forward));
                let end = Vertex(new_v1?.0.oriented(Orientation::Reversed));
                let wall = make_edge_between(&start, &end)?;
                self.map_ve
                    .insert(shape_map_key(&forward_v1), wall.clone()); // `cxx:576`.
                Edge(wall.0.oriented(Orientation::Reversed)) // `cxx:579`.
            }
        };

        // `cxx:584-589`: `MakeWire` + the four `Add`s, each edge keeping the
        // orientation written above.
        let builder = TopoBuilder::new();
        let wire = builder.make_wire(&[orig_copy, wall1.clone(), new_edge.clone(), wall2.clone()]);

        // `cxx:592`: `BRepLib::BuildCurves3d(aWire)`.
        let _ = build_curves3d(&wire.0);

        // `cxx:595-607`: the planar attempt. `BRepLib_MakeFace(aWire, true)`
        // answers a null face both when it throws and when it fails, so the
        // port's `None` covers both arms.
        if let Some(face) = make_face_from_wire(&wire, true) {
            return Some(face);
        }

        // `cxx:611-663`: the thrusections fallback.
        // `cxx:613-618`: the two trimmed section curves. `ToReverse` is
        // `false` (`cxx:612`) and stays so.
        let edge_curve = BRepTool::edge_curve(the_orig_edge)?;
        let (fpar, lpar) = BRepTool::edge_parameters(the_orig_edge);
        let tr_edge_curve: Arc<dyn Curve> =
            Arc::new(GeomTrimmedCurve::new(edge_curve, fpar, lpar));
        let offset_curve = BRepTool::edge_curve(&new_edge)?;
        let (fpar_oe, lpar_oe) = BRepTool::edge_parameters(&new_edge);
        let tr_offset_curve: Arc<dyn Curve> =
            Arc::new(GeomTrimmedCurve::new(offset_curve, fpar_oe, lpar_oe));

        // `cxx:620-624`: `GeomFill_Generator` + `Perform(PConfusion)`.
        let the_surf = generator_surface(&[tr_edge_curve, tr_offset_curve]).ok()?;

        // `cxx:626-627`: `theSurf->Bounds(Uf, Ul, Vf, Vl)`.
        let (uf, ul) = the_surf.u_range();
        let (vf, vl) = the_surf.v_range();

        let reg = GeometryRegistry::global();
        let dir_x = GpDir2d::new(1.0, 0.0).ok()?;

        // `cxx:629-631`: the original edge's pcurve is the line `v = Vf`.
        let edge_line2d = Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, vf), dir_x);
        reg.update_edge_pcurve_on_surface(
            &the_orig_edge.0,
            &the_surf,
            Some(Arc::new(edge_line2d)),
            CONFUSION,
        );

        // `cxx:632-633`: the offset edge's pcurve is the line `v = Vl`.
        let oe_line2d = Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, vl), dir_x);
        reg.update_edge_pcurve_on_surface(
            &new_edge.0,
            &the_surf,
            Some(Arc::new(oe_line2d)),
            CONFUSION,
        );

        // `cxx:634-635`: `ToReverse` is false, so `UonV1 = Uf`, `UonV2 = Ul`.
        let (uon_v1, uon_v2) = (uf, ul);
        let dir_y = GpDir2d::new(0.0, 1.0).ok()?;
        // `cxx:636-637`.
        let a_line2d = Geom2dLine::from_pnt_dir(GpPnt2d::new(uon_v2, 0.0), dir_y);
        let a_line2d2 = Geom2dLine::from_pnt_dir(GpPnt2d::new(uon_v1, 0.0), dir_y);

        if is_same(&wall1.0, &wall2.0) {
            // `cxx:638-644`: one wall edge carrying both pcurves (a seam).
            reg.set_edge_pcurves_on_surface(
                &wall1.0,
                &the_surf,
                vec![Arc::new(a_line2d), Arc::new(a_line2d2)],
                CONFUSION,
            );
            // `cxx:641-642`: `UpdateEdge(aWall1, theSurf->UIso(Uf), Confusion)`.
            if let Some(bspl_c34) = the_surf.u_iso_curve(uf) {
                reg.set_edge_curve3d(&wall1.0, bspl_c34, CONFUSION);
            }
            // `cxx:643`: `Range(aWall1, Vf, Vl)` writes every representation.
            reg.set_ranges_all(&wall1.0, vf, vl);
        } else {
            // `cxx:647-650`.
            reg.set_same_parameter(&wall1.0, false);
            reg.set_same_range(&wall1.0, false);
            reg.set_same_parameter(&wall2.0, false);
            reg.set_same_range(&wall2.0, false);
            // `cxx:651-652`.
            reg.update_edge_pcurve_on_surface(
                &wall1.0,
                &the_surf,
                Some(Arc::new(a_line2d)),
                CONFUSION,
            );
            reg.set_pcurve_range_on_surface(&wall1.0, &the_surf, vf, vl);
            // `cxx:653-654`.
            reg.update_edge_pcurve_on_surface(
                &wall2.0,
                &the_surf,
                Some(Arc::new(a_line2d2)),
                CONFUSION,
            );
            reg.set_pcurve_range_on_surface(&wall2.0, &the_surf, vf, vl);
            // `cxx:655-657`.
            if let Some(bspl_c3) = the_surf.u_iso_curve(uon_v2) {
                reg.set_edge_curve3d(&wall1.0, bspl_c3, CONFUSION);
            }
            reg.set_edge_range(&wall1.0, vf, vl); // `Range(..., true)`: only 3D.
            // `cxx:658-660`.
            if let Some(bspl_c4) = the_surf.u_iso_curve(uon_v1) {
                reg.set_edge_curve3d(&wall2.0, bspl_c4, CONFUSION);
            }
            reg.set_edge_range(&wall2.0, vf, vl);
        }

        // `cxx:663`.
        Some(make_face_from_surface_and_wire(the_surf, &wire))
    }

    /// `Generated(theShape)` (`cxx:671-686`): the shape the modification made
    /// of `theShape`, after the shape-healing modifications.
    pub fn generated(&self, shape: &TopoShape) -> Option<TopoShape> {
        // `cxx:674-681`: a null image stays null.
        let res = self.builder.modified_shape(shape).ok()?;
        // `cxx:684`: `aRes = myReShape->Apply(aRes)`.
        let mut reshape = self.reshape.clone();
        Some(reshape.apply(&res))
    }

    /// `Modified(theShape)` (`cxx:690-703`): the shape when a modification is
    /// recorded for it or its sub-shapes, `None` otherwise.
    pub fn modified(&mut self, shape: &TopoShape) -> Option<TopoShape> {
        // `cxx:694-700`: a status of `0` means no modification.
        self.reshape.status_last(shape)
    }
}
