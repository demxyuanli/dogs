//! `BRepOffsetAPI_MakeOffset` — the 2D offset of a planar face or wire.
//!
//! Source: `ModelingAlgorithms/TKOffset/BRepOffsetAPI/BRepOffsetAPI_MakeOffset.cxx:93-505`
//! and `.hxx` (OCCT 8.0.0).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `NeedsConvertion(Wire)` | `cxx:47-61` | [`needs_convertion`] |
//! | `ConvertFace(F, tol)` | `cxx:63-89` | [`convert_face`] |
//! | `BRepOffsetAPI_MakeOffset()` | `cxx:93-99` | [`BRepOffsetApiMakeOffset::new`] |
//! | `(Face, Join, IsOpenResult)` | `cxx:104-108` | [`BRepOffsetApiMakeOffset::with_face`] |
//! | `Init(Face, Join, IsOpenResult)` | `cxx:112-127` | [`BRepOffsetApiMakeOffset::init`] |
//! | `(Wire, Join, IsOpenResult)` | `cxx:131-139` | [`BRepOffsetApiMakeOffset::with_wire`] |
//! | `Init(Join, IsOpenResult)` | `cxx:145-148` | [`BRepOffsetApiMakeOffset::init_join`] |
//! | `SetApprox(ToApprox)` | `cxx:157-160` | [`BRepOffsetApiMakeOffset::set_approx`] |
//! | `AddWire(Wire)` | `cxx:164-168` | [`BRepOffsetApiMakeOffset::add_wire`] |
//! | `BuildDomains(...)` | `cxx:172-338` | [`BRepOffsetApiMakeOffset::build_domains`] |
//! | `Perform(Offset, Alt)` | `cxx:342-474` | [`BRepOffsetApiMakeOffset::perform`] |
//! | `Build(Range)` | `cxx:478-481` | [`BRepOffsetApiMakeOffset::build`] |
//! | `Generated(S)` | `cxx:485-505` | [`BRepOffsetApiMakeOffset::generated`] |
//!
//! # Translation notes
//!
//! `NeedsConvertion`'s edge iterator is `TopoDS_Iterator` in OCCT
//! (`cxx:48-49`) and `TopExp_Explorer` here; a wire has only edges, so the two
//! agree. `ConvertFace` is declared `static` in OCCT (`hxx:95`) and touches no
//! instance state, so the port keeps it as the module-level [`convert_face`].
//! The `myIsToApprox` branch's `aFace` can be left null when no wire
//! yields a face (`cxx:352-364`); OCCT then dereferences the null
//! `Geom_Surface` inside `ConvertWire` (`BRepAlgo.cxx:103`), which is undefined
//! behaviour — the port skips those wires instead.
//!
//! * `BuildDomains` `throw StdFail_NotDone` at `cxx:190` and `cxx:241` is
//!   caught by `Perform`'s `catch (Standard_Failure)`
//!   (`cxx:463-474`) and only produces `NotDone()`. The port returns `false`
//!   from [`BRepOffsetApiMakeOffset::build_domains`] and takes the same
//!   `NotDone` path.
//! * `StdFail_NotDone_Raise_if(!myIsInitialized, ...)` (`cxx:343`) happens
//!   *before* the `try`, so OCCT really raises it. The port leaves
//!   `IsDone() == false` instead of raising, which is how the rest of the port
//!   reports `StdFail_NotDone`.
//! * `BRepAlgo_FaceRestrictor FR; FR.Init(myFace, true)` (`cxx:214`) passes
//!   `Proj = true` with the *default* `ControlOrientation = false`
//!   (`BRepAlgo_FaceRestrictor.hxx:47-49`), i.e. the `WireToFace` path of
//!   `Perform`, not `PerformWithCorrection`.
//! * `BRepBuilderAPI_MakeFace(Wire, OnlyPlane)` (`cxx:186`) is
//!   [`crate::brep_lib_make_face::make_face_from_wire`].
//! * `BRepAdaptor_Surface S(F, false)` + `Extrema_ExtPS(P, S, Tol, Tol,
//!   Extrema_ExtFlag_MIN)` (`cxx:287-305`) is
//!   [`occt_geom::extrema_surf::ExtPs`] over the face's surface with
//!   `ExtremaExtFlag::Min` and the surface's own `u_range` / `v_range`, which
//!   is what the unrestricted adaptor reports.
//! * `BRepTopAdaptor_FClass2d` (`cxx:290`) is
//!   [`crate::brep_top_adaptor_fclass2d::BrepTopAdaptorFClass2d`].
//! * The `TopoDS_Compound Res` of `cxx:391-459` is collected as a `Vec` and
//!   built with `TopoBuilder::make_compound_of`; OCCT turns `myShape` into the
//!   compound once two or more domain shapes are produced (`cxx:449-451`).

use occt_core::gp::GpPnt2d;
use occt_core::kernel::geomabs::JoinType;
use occt_core::precision::Precision;
use occt_geom::extrema_surf::{ExtPs, ExtremaExtFlag};

use crate::abs::{Orientation, ShapeType};
use crate::brep_algo_face_restrictor::BRepAlgoFaceRestrictor;
use crate::brep_fill_offset_wire::BRepFillOffsetWire;
use crate::brep_lib_make_face::make_face_from_wire;
use crate::brep_top_adaptor_fclass2d::BrepTopAdaptorFClass2d;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::fclass2d::FaceState;
use crate::geom_abs_curve_type::{curve_kind, GeomAbsCurveType};
use crate::shape::{Face, TopoShape, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{explorer_shapes, is_same, wire_vertices};

/// `NeedsConvertion(theWire)` (`BRepOffsetAPI_MakeOffset.cxx:47-61`): true as
/// soon as one edge is neither a line nor a circle.
///
/// Callers: `ConvertFace` (`cxx:76`) and the `myIsToApprox` branch of `Perform`
/// (`cxx:367`).
fn needs_convertion(w: &TopoShape) -> bool {
    // `cxx:48-58`
    for e in explorer_shapes(w, ShapeType::Edge) {
        let kind = match BRepTool::edge_curve(&crate::shape::Edge(e.clone())) {
            Some(c) => curve_kind(c.as_ref()),
            None => GeomAbsCurveType::OtherCurve,
        };
        // `cxx:55-58`: `BRepAdaptor_Curve(anEdge).GetType()`.
        if kind != GeomAbsCurveType::Line && kind != GeomAbsCurveType::Circle {
            return true;
        }
    }
    // `cxx:60`
    false
}

/// `BRepOffsetAPI_MakeOffset::ConvertFace(theFace, theAngleTolerance)`
/// (`cxx:63-89`).
///
/// The face is copied with `EmptyCopied` (surface and tolerance only), every
/// child wire that [`needs_convertion`] flags is run through
/// `BRepAlgo::ConvertWire` with the face FORWARD, `BRepLib::BuildCurves3d`
/// materialises the 3D curves, the wire's own orientation is restored and the
/// result is added back; a wire that needs no conversion is added unchanged.
/// `cxx:86` then restores the face's orientation.
pub fn convert_face(the_face: &Face, angle_tolerance: f64) -> Face {
    // `cxx:65-68`: `anOr` saved, `aFace` FORWARD.
    let an_or = the_face.0.orientation();
    let a_face = Face(the_face.0.oriented(Orientation::Forward));
    let builder = TopoBuilder::new();
    let reg = GeometryRegistry::global();

    // `cxx:70`: `TopoDS_Face aNewFace = TopoDS::Face(aFace.EmptyCopied())`.
    let mut a_new_face = a_face.0.empty_copied();
    if let Some(g) = reg.face_geom(&a_face.0) {
        reg.set_face(&a_new_face, g);
    }
    reg.set_face_tolerance(&a_new_face, reg.face_tolerance(&a_face.0));

    // `cxx:72-85`: `TopoDS_Iterator anIter(aFace)`.
    let mut builder_face = Face(a_new_face.clone());
    for sub in explorer_shapes(&a_face.0, ShapeType::Wire) {
        let mut wire = Wire(sub);
        let wire = if needs_convertion(&wire.0) {
            // `cxx:78-83`.
            let wire_or = wire.0.orientation();
            wire.0.set_orientation(Orientation::Forward);
            let Some(mut converted) = crate::brep_algo::convert_wire(&wire, angle_tolerance, &a_face)
            else {
                // `cxx:84`: `aBB.Add(aNewFace, aWire)` with a null wire is a
                // no-op.
                continue;
            };
            crate::brep_lib_build_curves3d::build_curves3d(&converted.0);
            converted.0.set_orientation(wire_or);
            converted
        } else {
            wire
        };
        // `cxx:84`.
        builder.add_wire(&mut builder_face, &wire);
    }
    a_new_face = builder_face.0;

    // `cxx:86`.
    a_new_face.set_orientation(an_or);
    Face(a_new_face)
}

/// `class BRepOffsetAPI_MakeOffset` (`hxx:36-110`).
pub struct BRepOffsetApiMakeOffset {
    /// `myFace` (`hxx:88`).
    my_face: Option<Face>,
    /// `myWires` (`hxx:89`).
    my_wires: Vec<TopoShape>,
    /// `myLeft` (`hxx:91`).
    my_left: Vec<BRepFillOffsetWire>,
    /// `myRight` (`hxx:93`).
    my_right: Vec<BRepFillOffsetWire>,
    /// `myJoin` (`hxx:95`).
    my_join: JoinType,
    /// `myIsOpenResult` (`hxx:96`).
    my_is_open_result: bool,
    /// `myIsInitialized` (`hxx:97`).
    my_is_initialized: bool,
    /// `myIsToApprox` (`hxx:98`).
    my_is_to_approx: bool,
    /// `myLastIsLeft` (`hxx:99`).
    my_last_is_left: bool,
    /// `myGenerated` (`hxx:101`).
    my_generated: Vec<TopoShape>,
    /// `BRepBuilderAPI_MakeShape::myShape` (`BRepBuilderAPI_MakeShape.hxx`).
    my_shape: Option<TopoShape>,
    /// `BRepBuilderAPI_MakeShape::myDone`.
    my_done: bool,
}

impl Default for BRepOffsetApiMakeOffset {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepOffsetApiMakeOffset {
    /// `BRepOffsetAPI_MakeOffset()` (`cxx:93-99`): `myJoin` starts at
    /// `GeomAbs_Arc`, the other flags are false.
    pub fn new() -> Self {
        // `cxx:95-97`
        Self {
            my_face: None,
            my_wires: Vec::new(),
            my_left: Vec::new(),
            my_right: Vec::new(),
            my_join: JoinType::Arc,
            my_is_open_result: false,
            my_is_initialized: false,
            my_is_to_approx: false,
            my_last_is_left: false,
            my_generated: Vec::new(),
            my_shape: None,
            my_done: false,
        }
    }

    /// `BRepOffsetAPI_MakeOffset(Spine, Join, IsOpenResult)` (`cxx:104-108`).
    pub fn with_face(spine: &Face, join: JoinType, is_open_result: bool) -> Self {
        let mut s = Self::new();
        s.init(spine, join, is_open_result);
        s
    }

    /// `BRepOffsetAPI_MakeOffset(Spine, Join, IsOpenResult)` for a wire
    /// (`cxx:131-139`).
    pub fn with_wire(spine: &TopoShape, join: JoinType, is_open_result: bool) -> Self {
        let mut s = Self::new();
        // `cxx:135-138`
        s.my_wires.push(spine.clone());
        s.my_is_initialized = true;
        s.my_join = join;
        s.my_is_open_result = is_open_result;
        s.my_is_to_approx = false;
        s
    }

    /// `Init(Spine, Join, IsOpenResult)` (`cxx:112-127`).
    pub fn init(&mut self, spine: &Face, join: JoinType, is_open_result: bool) {
        // `cxx:117-120`
        self.my_face = Some(spine.clone());
        self.my_is_initialized = true;
        self.my_join = join;
        self.my_is_open_result = is_open_result;
        self.my_is_to_approx = false;
        // `cxx:121-126`: every wire of the face, in `TopExp_Explorer` order.
        for w in explorer_shapes(&spine.0, ShapeType::Wire) {
            self.my_wires.push(w);
        }
    }

    /// `Init(Join, IsOpenResult)` (`cxx:145-148`).
    pub fn init_join(&mut self, join: JoinType, is_open_result: bool) {
        self.my_join = join;
        self.my_is_open_result = is_open_result;
    }

    /// `SetApprox(ToApprox)` (`cxx:157-160`).
    ///
    /// The flag drives the `myIsToApprox` branch of `Perform`
    /// (`cxx:347-387`), see [`BRepOffsetApiMakeOffset::perform_body`].
    pub fn set_approx(&mut self, to_approx: bool) {
        self.my_is_to_approx = to_approx;
    }

    /// `AddWire(Spine)` (`cxx:164-168`).
    pub fn add_wire(&mut self, spine: &TopoShape) {
        // `cxx:166-167`
        self.my_is_initialized = true;
        self.my_wires.push(spine.clone());
    }

    /// `IsDone()` (`BRepBuilderAPI_MakeShape`).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `Shape()` (`BRepBuilderAPI_MakeShape`).
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_shape.as_ref()
    }

    /// `JoinType()` (`hxx:74`).
    pub fn join_type(&self) -> JoinType {
        self.my_join
    }

    /// `IsOpenResult()` (`hxx:79`).
    pub fn is_open_result(&self) -> bool {
        self.my_is_open_result
    }

    /// `Build(Range)` (`cxx:478-481`): `Done()`.
    pub fn build(&mut self) {
        self.my_done = true;
    }

    /// `Generated(S)` (`cxx:485-505`).
    pub fn generated(&mut self, s: &TopoShape) -> &[TopoShape] {
        // `cxx:487`
        self.my_generated.clear();
        // `cxx:488-494`: the list `Perform` last walked.
        let algos: &mut Vec<BRepFillOffsetWire> = if self.my_last_is_left {
            &mut self.my_left
        } else {
            &mut self.my_right
        };
        // `cxx:495-503`
        let forward = s.oriented(Orientation::Forward);
        let reversed = s.oriented(Orientation::Reversed);
        for ow in algos.iter_mut() {
            let mut l = ow.generated_shapes(&forward);
            self.my_generated.append(&mut l);
            let mut l = ow.generated_shapes(&reversed);
            self.my_generated.append(&mut l);
        }
        &self.my_generated
    }

    /// `BuildDomains(myFace, WorkWires, myAlgos, myJoin, myIsOpenResult,
    /// isPositive, isWasReversed)` (`cxx:172-338`).
    ///
    /// Answers `false` where OCCT throws `StdFail_NotDone` (`cxx:190`,
    /// `cxx:241`); `Perform` turns that into `NotDone()`.
    fn build_domains(
        my_face: &mut Option<Face>,
        work_wires: &mut Vec<TopoShape>,
        my_algos: &mut Vec<BRepFillOffsetWire>,
        join: JoinType,
        is_open_result: bool,
        is_positive: bool,
        is_was_reversed: &mut bool,
    ) -> bool {
        let reg = GeometryRegistry::global();
        let builder = TopoBuilder::new();
        let mut fr = BRepAlgoFaceRestrictor::new();
        let mut low: Vec<TopoShape> = Vec::new();

        // `cxx:186-192`: build a face from the first work wire when none was
        // given; a wire that is not planar throws.
        if my_face.is_none() {
            let Some(first) = work_wires.first().cloned() else {
                return false;
            };
            let Some(f) = make_face_from_wire(&Wire(first), true) else {
                return false;
            };
            *my_face = Some(f);
        }
        let face = my_face.clone().expect("face set above");
        let face_geom = reg.face_geom(&face.0);

        // `cxx:194-212`: reverse every work wire when the first one is
        // oriented like the face's own first wire (or unlike it, for
        // `isPositive` false).
        let Some(a_wire1) = work_wires.first().cloned() else {
            return false;
        };
        if let Some(a_wire2) = explorer_shapes(&face.0, ShapeType::Wire).into_iter().next() {
            let same = a_wire1.orientation() == a_wire2.orientation();
            let complement = a_wire1.orientation() == a_wire2.orientation().complement();
            if (same && is_positive) || (complement && !is_positive) {
                *work_wires = work_wires
                    .iter()
                    .map(|w| w.oriented(w.orientation().reversed()))
                    .collect();
                *is_was_reversed = true;
            }
        }

        // `cxx:214`: `FR.Init(myFace, true)` — `Proj = true`,
        // `ControlOrientation` stays false.
        fr.init(Some(face.clone()), true, false);

        // `cxx:219-237`
        for w in work_wires.iter() {
            if w.is_closed() {
                // `cxx:223-227`
                fr.add(w);
                continue;
            }
            let (vf, vl) = wire_vertices(&Wire(w.clone()));
            let closed_by_vertices = match (vf, vl) {
                (Some(a), Some(b)) => is_same(&a.0, &b.0),
                _ => false,
            };
            if closed_by_vertices {
                // `cxx:229-232`
                fr.add(w);
            } else {
                // `cxx:234-236`
                low.push(w.clone());
            }
        }

        // `cxx:238-242`
        fr.perform();
        if !fr.is_done() {
            return false;
        }
        // `cxx:243-258`
        let mut faces: Vec<TopoShape> = Vec::new();
        while fr.more() {
            if let Some(f) = fr.current() {
                faces.push(f.0);
            }
            fr.next();
        }

        // `cxx:260-277`: no closed wire means a single domain.
        if faces.is_empty() {
            let mut f = face.0.empty_copied();
            if let Some(g) = &face_geom {
                reg.set_face(&f, g.clone());
            }
            // `cxx:269-273`
            for w in &low {
                builder.add(&mut f, w);
            }
            // `cxx:274`
            my_algos.push(BRepFillOffsetWire::with_face(&Face(f), join, is_open_result));
            return true;
        }

        // `cxx:279-327`: give each open wire to the face that contains it.
        let tol = BRepTool::face_tolerance(&face);
        for i in 0..faces.len() {
            let Some(surface) = reg.face_surface(&faces[i]) else {
                continue;
            };
            let f_key = Face(faces[i].clone());
            // `cxx:290`
            let cl = BrepTopAdaptorFClass2d::new(&f_key, Precision::CONFUSION);

            // `cxx:292-326`
            let mut k = 0;
            while k < low.len() {
                let w = low[k].clone();
                // `cxx:299-301`: `TopExp_Explorer exp(W, TopAbs_VERTEX)` then
                // `P3d = BRep_Tool::Pnt(V)`. A wire without vertices cannot be
                // classified, so the port leaves it in `LOW` where OCCT would
                // dereference a null vertex.
                let Some(v) = explorer_shapes(&w, ShapeType::Vertex).into_iter().next() else {
                    k += 1;
                    continue;
                };
                let p3d = BRepTool::vertex_point(&Vertex(v));

                // `cxx:302-305`: `Extrema_ExtPS(P3d, S, Tol, Tol, MIN)`.
                let mut ext = ExtPs::new();
                ext.set_flag(ExtremaExtFlag::Min);
                let (u0, u1) = surface.u_range();
                let (v0, v1) = surface.v_range();
                ext.initialize(surface.as_ref(), u0, u1, v0, v1, tol, tol);
                ext.perform(&p3d);

                // `cxx:306-318`
                let mut dist2_min = Precision::INFINITE;
                let mut found = false;
                let mut pv = GpPnt2d::new(0.0, 0.0);
                for ie in 1..=ext.nb_ext() {
                    if ext.square_distance(ie) < dist2_min {
                        dist2_min = ext.square_distance(ie);
                        found = true;
                        let (x, y, _) = ext.point(ie);
                        pv = GpPnt2d::new(x, y);
                    }
                }

                // `cxx:319-325`
                if found && cl.perform(pv, false) == FaceState::In {
                    builder.add(&mut faces[i], &w);
                    low.remove(k);
                } else {
                    k += 1;
                }
            }
        }

        // `cxx:329-337`
        for f in &faces {
            my_algos.push(BRepFillOffsetWire::with_face(
                &Face(f.clone()),
                join,
                is_open_result,
            ));
        }
        true
    }

    /// `Perform(Offset, Alt)` (`cxx:342-474`).
    ///
    /// `cxx:345-462` is a `try` / `catch (Standard_Failure const&)` that leaves
    /// `myShape` null and the generator not done (`cxx:463-472`). The OCCT
    /// exceptions of that body are panics in the port, so the body runs under
    /// `catch_unwind` and the catch state is restored here.
    pub fn perform(&mut self, offset: f64, alt: f64) {
        // `cxx:343`: `StdFail_NotDone_Raise_if(!myIsInitialized, ...)`, raised
        // before the `try`.
        if !self.my_is_initialized {
            self.my_done = false;
            return;
        }
        let prev_hook = std::panic::take_hook();
        std::panic::set_hook(Box::new(|_| {}));
        let land = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            self.perform_body(offset, alt)
        }));
        std::panic::set_hook(prev_hook);
        if land.is_err() {
            self.my_done = false;
            self.my_shape = None;
        }
    }

    /// The `try` block of `Perform` (`cxx:345-462`).
    fn perform_body(&mut self, offset: f64, alt: f64) {
        // `cxx:347-387`: the `myIsToApprox` branch.
        if self.my_is_to_approx {
            // `cxx:349`: `double aTol = 0.1`.
            let a_tol = 0.1;
            if self.my_face.is_none() {
                // `cxx:352-364`: the first wire that yields a face.
                let mut a_face: Option<Face> = None;
                for w in &self.my_wires {
                    // `cxx:357-358`: `BRepBuilderAPI_MakeFace(Wire, OnlyPlane)`,
                    // `Error() == BRepBuilderAPI_FaceDone`.
                    if let Some(f) =
                        crate::brep_lib_make_face::make_face_from_wire(&Wire(w.clone()), true)
                    {
                        a_face = Some(f);
                        break;
                    }
                }
                // `cxx:365-374`. Without a face OCCT dereferences a null
                // `Geom_Surface` in `ConvertWire` (`BRepAlgo.cxx:103`); the
                // port leaves the wires alone instead.
                if let Some(a_face) = a_face {
                    for shape in self.my_wires.iter_mut() {
                        let wire = Wire(shape.clone());
                        if !needs_convertion(&wire.0) {
                            continue;
                        }
                        // `cxx:369-372`.
                        let Some(mut new_wire) =
                            crate::brep_algo::convert_wire(&wire, a_tol, &a_face)
                        else {
                            continue;
                        };
                        crate::brep_lib_build_curves3d::build_curves3d(&new_wire.0);
                        new_wire.0.set_orientation(wire.0.orientation());
                        *shape = new_wire.0;
                    }
                }
            } else if let Some(face) = self.my_face.clone() {
                // `cxx:378-385`.
                let new_face = convert_face(&face, a_tol);
                crate::brep_lib_build_curves3d::build_curves3d(&new_face.0);
                self.my_face = Some(new_face);
                self.my_wires.clear();
                let face = self.my_face.clone().unwrap();
                for sub in explorer_shapes(&face.0, ShapeType::Wire) {
                    self.my_wires.push(sub);
                }
            }
        }

        // `cxx:389-394`
        let mut i = 1;
        let mut result: Vec<TopoShape> = Vec::new();
        self.my_last_is_left = offset <= 0.0;
        let mut is_was_reversed = false;

        // `cxx:396-420`: `Offset <= 0` walks `myLeft`.
        if offset <= 0.0 {
            if self.my_left.is_empty() {
                if !Self::build_domains(
                    &mut self.my_face,
                    &mut self.my_wires,
                    &mut self.my_left,
                    self.my_join,
                    self.my_is_open_result,
                    false,
                    &mut is_was_reversed,
                ) {
                    // `cxx:463-472`: `NotDone(); myShape.Nullify();`
                    self.my_done = false;
                    self.my_shape = None;
                    return;
                }
            }
            for algo in self.my_left.iter_mut() {
                // `cxx:407`
                algo.perform(offset.abs(), alt);
                // `cxx:408-418`
                if algo.is_done() {
                    if let Some(shape) = algo.shape() {
                        let s = if is_was_reversed {
                            shape.oriented(shape.orientation().reversed())
                        } else {
                            shape.clone()
                        };
                        result.push(s.clone());
                        if i == 1 {
                            self.my_shape = Some(s);
                        }
                        i += 1;
                    }
                }
            }
        } else {
            // `cxx:422-447`: `Offset > 0` walks `myRight`.
            if self.my_right.is_empty() {
                if !Self::build_domains(
                    &mut self.my_face,
                    &mut self.my_wires,
                    &mut self.my_right,
                    self.my_join,
                    self.my_is_open_result,
                    true,
                    &mut is_was_reversed,
                ) {
                    self.my_done = false;
                    self.my_shape = None;
                    return;
                }
            }
            for algo in self.my_right.iter_mut() {
                // `cxx:433`
                algo.perform(offset, alt);
                // `cxx:435-445`
                if algo.is_done() {
                    if let Some(shape) = algo.shape() {
                        let s = if is_was_reversed {
                            shape.oriented(shape.orientation().reversed())
                        } else {
                            shape.clone()
                        };
                        result.push(s.clone());
                        if i == 1 {
                            self.my_shape = Some(s);
                        }
                        i += 1;
                    }
                }
            }
        }

        // `cxx:449-451`: two or more domain shapes are returned as one
        // compound.
        if i > 2 {
            let compound = TopoBuilder::new().make_compound_of(&result);
            self.my_shape = Some(compound.0);
        }

        // `cxx:453-461`
        self.my_done = self.my_shape.is_some();
    }
}
