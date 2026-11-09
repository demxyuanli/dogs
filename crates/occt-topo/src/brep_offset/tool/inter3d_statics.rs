//! File-local helpers that `BRepOffset_Tool::Inter3D` (`BRepOffset_Tool.cxx:
//! 1441-1941`) and its `AssembleEdge` (`:1345-1433`) call, plus the two
//! `GeomProjLib` entry points they depend on.
//!
//! Ported here (control flow read from
//! `src/ModelingAlgorithms/TKOffset/BRepOffset/BRepOffset_Tool.cxx`):
//!
//! - `UpdateVertexTolerances` (`:4613-4658`)
//! - `IsAutonomVertex(V, PDS, F1, F2)` (`:811-832`)
//! - `IsAutonomVertex(V, PDS)` (`:839-911`)
//! - `AreConnex` (`:918-928`), `AreClosed` (`:930-939`)
//! - `BSplineEdges` (`:942-987`), `AngleWireEdge` (`:991-1032`)
//! - `ReconstructPCurves` (`:1035-1057`) and `GeomProjLib::Curve2d(C, First,
//!   Last, S, Tol)` (`GeomProjLib.cxx:54-132`) as
//!   [`BRepOffsetTool::geom_projlib_curve2d`]
//! - `ConcatPCurves` (`:1062-1160`), `Glue` (`:1164-1259`) and
//!   `AssembleEdge` (`:1345-1433`)
//!
//! Still UNPORTED in this group: none. `PerformPlanes` (`:4537-4602`) lives in
//! `perform_planes.rs` (`IntTools_FaceFace` is `crate::int_face_face`), and
//! `CheckIntersFF` (`:1264-1343`) in `check_inters_ff.rs`; `Inter3D` itself is
//! `inter3d.rs`.
//!
//! Port guards, each mirroring a read OCCT performs on an undefined value:
//! `BSplineEdges` (`cxx:976-985`) returns `false` without writing `angle`, and
//! `AngleWireEdge` (`cxx:1011-1028`) then returns that uninitialised `Angle`;
//! the port answers `None` there instead of inventing an angle. `AngleWireEdge`
//! likewise leaves the first-edge vertices unset when no edge of the wire
//! touches `CV` (`cxx:995-1010`), so the port answers `None` as well.

use std::sync::Arc;

use occt_core::convert::ParameterisationType;
use occt_core::elib::clib2d::{
    parameter_circ2d, parameter_elips2d, parameter_hypr2d, parameter_lin2d, parameter_parab2d,
};
use occt_core::kernel::geomabs::{CurveType, Shape};
use occt_core::precision::{CONFUSION, PCONFUSION, REAL_SMALL};
use occt_geom::convert_approx_curve::GeomConvertApproxCurve;
use occt_geom::convert_bspl::CompCurveToBSplineCurve;
use occt_geom::projlib::ProjectedCurve;
use occt_geom::{Curve, GeomTrimmedCurveBasis, Surface};
use occt_geom2d::comp_curve_to_bspline::CompCurveToBSplineCurve2d;
use occt_geom2d::convert_approx_curve::Geom2dConvertApproxCurve;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;
use occt_geom2d::{
    Curve2d, Geom2dCircle, Geom2dEllipse, Geom2dHyperbola, Geom2dLine, Geom2dParabola,
};

use crate::bopds::BopdsDS;
use crate::brep_lib_make_edge;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::geom_abs_curve_type::{curve_kind, GeomAbsCurveType};
use crate::geom_bnd_lib_curve2d::curve2d_abs_type;
use crate::inttools_data::CommonPartType;
use crate::shape::{Edge, Face, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{is_same, wire_vertices};

use super::BRepOffsetTool;

/// `GeomAbs_Shape` code of `C1` (`GeomAbs_Shape.hxx`); `Continuity() < C1`
/// is OCCT's `Continuity() < GeomAbs_C1` in `ConcatPCurves` / `Glue`.
const GEOM_ABS_C1: u8 = 2;

/// `gp::Resolution()`, the derivative-magnitude floor of `BSplineEdges`
/// (`cxx:976`).
const GP_RESOLUTION: f64 = REAL_SMALL;

/// True for `V1.IsSame(V2)` with OCCT's null-handle rule: two null shapes are
/// `IsSame` (`TopoDS_Shape::IsSame`).
fn same_vertex(a: &Option<Vertex>, b: &Option<Vertex>) -> bool {
    match (a.as_ref(), b.as_ref()) {
        (Some(x), Some(y)) => is_same(&x.0, &y.0),
        (None, None) => true,
        _ => false,
    }
}

/// `TopoDS_Vertex::IsSame(V)` for one optional and one concrete vertex.
fn is_same_vertex(v: &Option<Vertex>, w: &Vertex) -> bool {
    v.as_ref().map_or(false, |x| is_same(&x.0, &w.0))
}

/// `BB.UpdateVertex(V, Tol)` (`BRep_Builder.cxx:1113-1125`): the stored
/// tolerance is the larger of `tol` and the current one.
fn raise_vertex_tolerance(v: &Vertex, tol: f64) {
    if tol > BRepTool::vertex_tolerance(v) {
        v.set_tolerance(tol);
    }
}

/// The trimmed-unwrapped view of a 3D curve, as the `Geom_TrimmedCurve`
/// down-casts of `cxx:952-963` / `cxx:1190-1198` do. A curve that is not
/// trimmed is returned unchanged.
fn basis_3d(c: &Arc<dyn Curve>) -> Arc<dyn Curve> {
    match c.untrimmed_basis() {
        Some((basis, _, _)) => basis,
        None => c.clone(),
    }
}

/// The trimmed-unwrapped view of a 2D curve (`cxx:1078-1087`).
fn basis_2d(c: &Arc<dyn Curve2d>) -> Arc<dyn Curve2d> {
    match c.trimmed_basis() {
        Some(basis) => Arc::from(basis.clone_dyn()),
        None => c.clone(),
    }
}

/// `Geom_Line` / `Geom_Conic` membership, the `IsInstance(Geom_Line) ||
/// IsKind(Geom_Conic)` test of `cxx:1207-1208` and `cxx:1096-1098`.
fn is_line_or_conic(t: GeomAbsCurveType) -> bool {
    matches!(
        t,
        GeomAbsCurveType::Line
            | GeomAbsCurveType::Circle
            | GeomAbsCurveType::Ellipse
            | GeomAbsCurveType::Parabola
            | GeomAbsCurveType::Hyperbola
    )
}

/// The same test on the 2D side (`Geom2d_Line` / `Geom2d_Conic`).
fn is_line_or_conic_2d(c: &dyn Curve2d) -> bool {
    matches!(
        curve2d_abs_type(c),
        crate::geom_bnd_lib_curve2d::GeomAbsCurve2d::Line
            | crate::geom_bnd_lib_curve2d::GeomAbsCurve2d::Circle
            | crate::geom_bnd_lib_curve2d::GeomAbsCurve2d::Ellipse
            | crate::geom_bnd_lib_curve2d::GeomAbsCurve2d::Parabola
            | crate::geom_bnd_lib_curve2d::GeomAbsCurve2d::Hyperbola
    )
}

impl BRepOffsetTool {
    /// `UpdateVertexTolerances(theFace)` (`cxx:4613-4658`): raises the
    /// tolerance of every vertex of `face` to the largest distance between the
    /// vertex and the vertex parameters of its incident edges, in 3D and on
    /// the face's surface (through the pcurves).
    ///
    /// `TopExp::MapShapesAndAncestors(theFace, VERTEX, EDGE, VEmap)` is
    /// [`crate::topo_tools::map_shapes_and_ancestors`] in OCCT order; the port
    /// walks the face's own vertices and the face's edges, which yields the
    /// same `(vertex, incident edges)` pairs.
    ///
    /// A missing 3D curve or pcurve is the null handle `BRepAdaptor_Curve`
    /// would dereference (it raises in OCCT); the port skips that evaluation.
    pub fn update_vertex_tolerances(face: &Face) {
        let edges = crate::topo_tools::edges_of(&face.0);
        let surface = BRepTool::face_surface(face);
        for v in crate::topo_tools::vertices_of(&face.0) {
            let pnt_vtx = BRepTool::vertex_point_world(&v);
            for e in &edges {
                let (v1, v2) = Self::edge_vertices(e);
                if !is_same_vertex(&v1, &v) && !is_same_vertex(&v2, &v) {
                    continue;
                }
                // `cxx:4633-4635`.
                let (fpar, lpar) = BRepTool::edge_parameters(e);
                let a_param = if is_same_vertex(&v1, &v) { fpar } else { lpar };
                let closed = same_vertex(&v1, &v2);

                // `cxx:4636-4649`: the 3D curve, skipped when degenerated.
                if !BRepTool::is_degenerated(e) {
                    if let Some(c) = BRepTool::edge_curve_world(e) {
                        let p = c.d0(a_param);
                        raise_vertex_tolerance(&v, pnt_vtx.distance(&p));
                        if closed {
                            let p = c.d0(lpar);
                            raise_vertex_tolerance(&v, pnt_vtx.distance(&p));
                        }
                    }
                }

                // `cxx:4650-4658`: `BRepAdaptor_Curve(anEdge, theFace)`, i.e.
                // `S.Value(C2d.Value(U))` at the same parameters.
                if let (Some(s), Some((pc, _, _))) =
                    (surface.as_ref(), crate::boptools_2d::curve_on_surface_range(e, face))
                {
                    let uv = pc.d0(a_param);
                    let p = s.d0(uv.x(), uv.y());
                    raise_vertex_tolerance(&v, pnt_vtx.distance(&p));
                    if closed {
                        let uv = pc.d0(lpar);
                        let p = s.d0(uv.x(), uv.y());
                        raise_vertex_tolerance(&v, pnt_vtx.distance(&p));
                    }
                }
            }
        }
    }

    /// `IsAutonomVertex(theVertex, thePDS, theFace1, theFace2)` (`cxx:811-832`):
    /// the vertex is created by the pave filler and is not a `VerticesOn`
    /// vertex of either face.
    ///
    /// OCCT answers `false` when the face info contains the vertex and `true`
    /// when neither does. A vertex or face missing from the DS is an
    /// out-of-range `pDS->Index` / `FaceInfo` in OCCT; the port answers `false`
    /// (not autonom, so nothing is glued on unindexed input).
    pub fn is_autonom_vertex_on_faces(
        vertex: &Vertex,
        pds: &BopdsDS,
        f1: &Face,
        f2: &Face,
    ) -> bool {
        let Some(nv) = pds.index(&vertex.0) else {
            return false;
        };
        for face in [f1, f2] {
            let Some(nf) = pds.index(&face.0) else {
                return false;
            };
            let Some(info) = pds.face_info(nf) else {
                return false;
            };
            if info.verts_on().contains(&nv) {
                return false;
            }
        }
        true
    }

    /// `IsAutonomVertex(aVertex, pDS)` (`cxx:839-911`): `true` only for a new
    /// shape that is the `IndexNew` of no vertex/vertex, edge/edge (vertex
    /// common part) or edge/face (vertex common part) interference.
    ///
    /// `cxx:843-859` falls back to a linear scan over
    /// `[NbSourceShapes, NbShapes)` when `Index` answers `-1`; a vertex that is
    /// still unindexed makes OCCT read `Shape(-1)`, so the port answers
    /// `false` there.
    pub fn is_autonom_vertex(vertex: &Vertex, pds: &BopdsDS) -> bool {
        let index = match pds.index(&vertex.0) {
            Some(i) => i,
            None => {
                let mut found = None;
                for i in pds.nb_source_shapes()..pds.nb_shapes() {
                    if let Some(s) = pds.shape(i) {
                        if is_same(s, &vertex.0) {
                            found = Some(i);
                            break;
                        }
                    }
                }
                match found {
                    Some(i) => i,
                    None => return false,
                }
            }
        };
        // `cxx:862-865`.
        if !pds.is_new_shape(index) {
            return false;
        }
        let as_new = index as i64;
        // `cxx:867-879`: VV.
        for vv in pds.interf_vv() {
            if vv.has_index_new() && vv.index_new() == as_new {
                return false;
            }
        }
        // `cxx:881-893`: EE with a vertex common part.
        for ee in pds.interf_ee() {
            let is_vertex = ee
                .common_part
                .as_ref()
                .map_or(false, |cp| cp.part_type() == CommonPartType::Vertex);
            if is_vertex && ee.index_new() == as_new {
                return false;
            }
        }
        // `cxx:895-907`: EF with a vertex common part.
        for ef in pds.interf_ef() {
            let is_vertex = ef
                .common_part
                .as_ref()
                .map_or(false, |cp| cp.part_type() == CommonPartType::Vertex);
            if is_vertex && ef.index_new() == as_new {
                return false;
            }
        }
        true
    }

    /// `AreConnex(W1, W2)` (`cxx:918-928`): the two wires share an endpoint
    /// vertex.
    pub fn are_connex(w1: &Wire, w2: &Wire) -> bool {
        let (v11, v12) = wire_vertices(w1);
        let (v21, v22) = wire_vertices(w2);
        same_vertex(&v11, &v21)
            || same_vertex(&v11, &v22)
            || same_vertex(&v12, &v21)
            || same_vertex(&v12, &v22)
    }

    /// `AreClosed(E1, E2)` (`cxx:930-939`): the two edges join the same two
    /// vertices, in either order.
    pub fn are_closed(e1: &Edge, e2: &Edge) -> bool {
        let (v11, v12) = Self::edge_vertices(e1);
        let (v21, v22) = Self::edge_vertices(e2);
        (same_vertex(&v11, &v21) && same_vertex(&v12, &v22))
            || (same_vertex(&v11, &v22) && same_vertex(&v12, &v21))
    }

    /// `BSplineEdges(E1, E2, par1, par2, angle)` (`cxx:942-987`): the angle
    /// between the first derivatives of the two B-spline edges at the chosen
    /// ends (`0` = first, anything else = last), or `PI/2` when a derivative
    /// vanishes.
    ///
    /// `None` is OCCT's `false` (a curve that is not a `Geom_BSplineCurve`
    /// after the trimmed unwrap): `cxx:976-985` leaves `angle` untouched, so no
    /// value exists to return.
    pub fn bspline_edges(e1: &Edge, e2: &Edge, par1: i32, par2: i32) -> Option<f64> {
        // `cxx:951-963`: the stored curves, then their trimmed unwrap.
        let c1 = basis_3d(&BRepTool::edge_curve_world(e1)?);
        let c2 = basis_3d(&BRepTool::edge_curve_world(e2)?);
        let (first1, last1) = BRepTool::edge_parameters(e1);
        let (first2, last2) = BRepTool::edge_parameters(e2);
        // `cxx:965-968`.
        if curve_kind(c1.as_ref()) != GeomAbsCurveType::BSplineCurve
            || curve_kind(c2.as_ref()) != GeomAbsCurveType::BSplineCurve
        {
            return None;
        }
        // `cxx:969-970`.
        let param1 = if par1 == 0 { first1 } else { last1 };
        let param2 = if par2 == 0 { first2 } else { last2 };
        // `cxx:972-975`: `C->D1(Param)` on the *basis* curves.
        let (_, der1) = c1.d1(param1);
        let (_, der2) = c2.d1(param2);
        // `cxx:977-985`.
        if der1.magnitude() <= GP_RESOLUTION || der2.magnitude() <= GP_RESOLUTION {
            Some(std::f64::consts::FRAC_PI_2)
        } else {
            Some(der1.angle(&der2))
        }
    }

    /// `AngleWireEdge(aWire, anEdge)` (`cxx:991-1032`): the angle between the
    /// wire edge incident to the shared vertex `CV` and `anEdge`, measured at
    /// `CV`.
    ///
    /// `None` is OCCT's undefined read: the first-edge search (`:995-1010`)
    /// leaves the endpoint pair unset when no edge of `aWire` touches `CV`, and
    /// `BSplineEdges` (`:1011-1028`) can leave `Angle` unwritten.
    pub fn angle_wire_edge(wire: &Wire, edge: &Edge) -> Option<f64> {
        // `cxx:993-996`.
        let (w11, w12) = wire_vertices(wire);
        let (v21, v22) = Self::edge_vertices(edge);
        let cv = if same_vertex(&w11, &v21) || same_vertex(&w11, &v22) {
            w11.clone()?
        } else {
            w12.clone()?
        };

        // `cxx:997-1010`: the first edge of the wire that touches `CV`.
        let mut found: Option<(Edge, Option<Vertex>, Option<Vertex>)> = None;
        for fe in crate::topo_tools::edges_of_wire(wire) {
            let (v1, v2) = Self::edge_vertices(&fe);
            if is_same_vertex(&v1, &cv) || is_same_vertex(&v2, &cv) {
                found = Some((fe, v1, v2));
                break;
            }
        }
        let (first_edge, v11, v12) = found?;

        // `cxx:1011-1029`.
        if is_same_vertex(&v11, &cv) && is_same_vertex(&v21, &cv) {
            Some(std::f64::consts::PI - Self::bspline_edges(&first_edge, edge, 0, 0)?)
        } else if is_same_vertex(&v11, &cv) && is_same_vertex(&v22, &cv) {
            Self::bspline_edges(&first_edge, edge, 0, 1)
        } else if is_same_vertex(&v12, &cv) && is_same_vertex(&v21, &cv) {
            Self::bspline_edges(&first_edge, edge, 1, 0)
        } else {
            Some(std::f64::consts::PI - Self::bspline_edges(&first_edge, edge, 1, 1)?)
        }
    }

    /// `GeomProjLib::Curve2d(C, First, Last, S, Tolerance)`
    /// (`GeomProjLib.cxx:54-132`): the 2D image of `C` on `S`.
    ///
    /// The five analytic arms (`:80-107`) build the `Geom2d_*` object from the
    /// projector data, the Bezier/BSpline arms hand back the projector's own
    /// curve, and the `default` arm (`:108-109`) returns a null handle. The
    /// trailing `Geom_TrimmedCurve` arm (`:119-129`) clips the image to the
    /// trim range on non-periodic results.
    ///
    /// `None` groups the null handle (`:108-109`), the projector-type mismatch
    /// of the `Line()`..`Hyperbola()` accessors (`Standard_NoSuchObject`), and
    /// a trim range the `Geom2d_TrimmedCurve` constructor rejects.
    pub fn geom_projlib_curve2d(
        c: &Arc<dyn Curve>,
        first: f64,
        last: f64,
        s: &Arc<dyn Surface>,
        tol: f64,
    ) -> Option<Arc<dyn Curve2d>> {
        // `cxx:64`.
        let tol = tol.max(PCONFUSION);
        // `cxx:66-71`: `GeomAdaptor_Curve(C, First, Last)` is the bounded view
        // over the same parameterisation.
        let view: Arc<dyn Curve> = if c.first_parameter() == first && c.last_parameter() == last {
            c.clone()
        } else {
            Arc::new(GeomTrimmedCurveBasis::new(c.clone(), first, last))
        };
        // `cxx:72`: `ProjLib_ProjectedCurve(HS, HC, Tolerance)`.
        let proj = ProjectedCurve::with_surface_curve_tol(s.clone(), view, tol);
        // `cxx:77-110`.
        let c2d: Arc<dyn Curve2d> = match proj.get_type() {
            CurveType::Line => {
                let lin = proj.line().ok()?;
                Arc::new(Geom2dLine::from_pnt_dir(lin.location(), *lin.direction()))
            }
            CurveType::Circle => Arc::new(Geom2dCircle::new(proj.circle().ok()?)),
            CurveType::Ellipse => Arc::new(Geom2dEllipse::new(proj.ellipse().ok()?)),
            CurveType::Parabola => Arc::new(Geom2dParabola::new(proj.parabola().ok()?)),
            CurveType::Hyperbola => Arc::new(Geom2dHyperbola::new(proj.hyperbola().ok()?)),
            CurveType::BezierCurve => proj.bezier()?,
            CurveType::BSplineCurve => proj.bspline()?,
            _ => return None,
        };
        // `cxx:119-129`.
        if c.is_geom_trimmed() {
            let mut u1 = c.first_parameter();
            let mut u2 = c.last_parameter();
            if !c2d.is_periodic() {
                u1 = u1.max(c2d.first_parameter());
                u2 = u2.min(c2d.last_parameter());
            }
            if u1 > u2 {
                // Port guard: `Geom2d_TrimmedCurve` raises on an inverted range.
                return None;
            }
            return Some(Arc::new(Geom2dTrimmedCurve::new(c2d, u1, u2)));
        }
        Some(c2d)
    }

    /// `ReconstructPCurves(anEdge)` (`cxx:1035-1057`): reprojects the edge's 3D
    /// curve onto every surface it already carries a pcurve for, replacing the
    /// pcurve in place.
    ///
    /// `cxx:1047-1050` builds `theSurf = Surface->Transformed(Location)` from
    /// the edge's location times the representation's. The port keeps face
    /// surfaces in world coordinates and carries no per-representation
    /// location, so the transform applies only when the edge itself has a
    /// non-identity location; the replacement is keyed by the stored surface
    /// handle, which is OCCT's representation identity.
    pub fn reconstruct_pcurves(edge: &Edge) {
        let Some(c3d) = BRepTool::edge_curve_world(edge) else {
            return;
        };
        let (f, l) = BRepTool::edge_parameters(edge);
        let reg = GeometryRegistry::global();
        let loc = edge.0.location();
        let trsf = loc.transformation();
        for (surf, _pc, _a, _b) in reg.edge_pcurve_reps(&edge.0) {
            // `cxx:1048-1050`.
            let target: Arc<dyn Surface> = if loc.is_identity() {
                surf.clone()
            } else {
                Arc::from(surf.transformed(&trsf))
            };
            // `cxx:1051`.
            if let Some(pc) = Self::geom_projlib_curve2d(&c3d, f, l, &target, PCONFUSION) {
                // `cxx:1052-1055`: `CurveRep->PCurve(ProjPCurve)`.
                reg.update_edge_pcurve_on_surface(&edge.0, &surf, Some(pc), 0.0);
            }
        }
    }

    /// `ConcatPCurves(E1, E2, F, After, newFirst, newLast)`
    /// (`cxx:1062-1160`): the pcurve joining the images of `E1` and `E2` on
    /// `F`, with the resulting COS range.
    ///
    /// `None` is the null `newPCurve` OCCT's `Glue` would append (`cxx:1248`):
    /// a failed 2D concatenation (`Add` false in OCCT is dereferenced at
    /// `:1149`) or an inverted trim range.
    pub fn concat_pcurves(
        e1: &Edge,
        e2: &Edge,
        face: &Face,
        after: bool,
    ) -> Option<(Arc<dyn Curve2d>, f64, f64)> {
        // `cxx:1073-1087`.
        let (pc1, first1, last1) = crate::boptools_2d::curve_on_surface_range(e1, face)?;
        let (pc2, mut first2, mut last2) = crate::boptools_2d::curve_on_surface_range(e2, face)?;
        let pcurve1 = basis_2d(&pc1);
        let pcurve2 = basis_2d(&pc2);

        // `cxx:1089-1092`: the same representation. OCCT compares the handles
        // after the trimmed unwrap, so the comparison is on the curve object
        // identity (`addr_eq` of the `dyn` data pointer).
        if std::ptr::addr_eq(pcurve1.as_ref(), pcurve2.as_ref()) {
            return Some((
                pcurve1,
                first1.min(first2),
                last1.max(last2),
            ));
        }

        // `cxx:1093-1139`: the same concrete class, `Geom2d_Line` or a conic.
        if curve2d_abs_type(pcurve1.as_ref()) == curve2d_abs_type(pcurve2.as_ref())
            && is_line_or_conic_2d(pcurve1.as_ref())
        {
            // `cxx:1100-1139`.
            let p1 = pcurve2.d0(first2);
            let p2 = pcurve2.d0(last2);
            if let Some(lin) = pcurve1.gp_lin2d() {
                first2 = parameter_lin2d(&lin, &p1);
                last2 = parameter_lin2d(&lin, &p2);
            } else if let Some(circ) = pcurve1.gp_circ2d() {
                first2 = parameter_circ2d(&circ, &p1);
                last2 = parameter_circ2d(&circ, &p2);
            } else if let Some(elips) = pcurve1.gp_elips2d() {
                first2 = parameter_elips2d(&elips, &p1);
                last2 = parameter_elips2d(&elips, &p2);
            } else if let Some(parab) = pcurve1.gp_parab2d() {
                first2 = parameter_parab2d(&parab, &p1);
                last2 = parameter_parab2d(&parab, &p2);
            } else if let Some(hypr) = pcurve1.gp_hypr2d() {
                first2 = parameter_hypr2d(&hypr, &p1);
                last2 = parameter_hypr2d(&hypr, &p2);
            }
            return Some((pcurve1, first1.min(first2), last1.max(last2)));
        }

        // `cxx:1141-1158`: concatenate the trimmed pieces and approximate when
        // the result is not at least C1.
        let tc1: Arc<dyn Curve2d> = Arc::new(Geom2dTrimmedCurve::new(pcurve1, first1, last1));
        let tc2: Arc<dyn Curve2d> = Arc::new(Geom2dTrimmedCurve::new(pcurve2, first2, last2));
        let mut concat =
            CompCurveToBSplineCurve2d::from_curve(tc1.as_ref(), ParameterisationType::TgtThetaOver2)?;
        // `cxx:1146`: `Add(TC2, Precision::Confusion(), After)`; `false` leaves
        // the accumulated curve as it was, and OCCT still reads
        // `BSplineCurve()` at `:1148`.
        concat.add(tc2.as_ref(), CONFUSION, after).ok()?;
        let bspl = concat.into_curve()?;
        let mut new_pcurve: Arc<dyn Curve2d> = Arc::new(bspl);
        // `cxx:1149-1155`.
        if new_pcurve.continuity() < GEOM_ABS_C1 {
            let approx = Geom2dConvertApproxCurve::new(
                new_pcurve.as_ref(),
                1.0e-7,
                Shape::C1,
                16,
                14,
            );
            if approx.has_result() {
                if let Some(c) = approx.curve() {
                    new_pcurve = Arc::new(c.clone());
                }
            }
        }
        // `cxx:1157-1158`.
        let new_first = new_pcurve.first_parameter();
        let new_last = new_pcurve.last_parameter();
        Some((new_pcurve, new_first, new_last))
    }

    /// `Glue(E1, E2, Vfirst, Vlast, After, F1, addPCurve1, F2, addPCurve2,
    /// theGlueTol)` (`cxx:1164-1259`): one edge over the union of the two
    /// curves, carrying the concatenated pcurve on each requested face.
    ///
    /// `None` is the null edge OCCT returns at `cxx:1219` (a failed 3D
    /// concatenation) and at `cxx:1240` (`BRepLib_MakeEdge` not done).
    #[allow(clippy::too_many_arguments)]
    pub fn glue(
        e1: &Edge,
        e2: &Edge,
        vfirst: &Vertex,
        vlast: &Vertex,
        after: bool,
        f1: &Face,
        add_pcurve1: bool,
        f2: &Face,
        add_pcurve2: bool,
        glue_tol: f64,
    ) -> Option<Edge> {
        let builder = TopoBuilder::new();
        // `cxx:1188-1198`.
        let c1 = basis_3d(&BRepTool::edge_curve_world(e1)?);
        let c2 = basis_3d(&BRepTool::edge_curve_world(e2)?);
        let (first1, last1) = BRepTool::edge_parameters(e1);
        let (first2, last2) = BRepTool::edge_parameters(e2);

        let mut is_canonic = false;
        let mut fparam = 0.0;
        let mut lparam = 0.0;
        let new_curve: Arc<dyn Curve>;
        // `cxx:1200-1205`: `C1 == C2` on the unwrapped curve handles.
        if std::ptr::addr_eq(c1.as_ref(), c2.as_ref()) {
            // `cxx:1200-1205`.
            new_curve = c1;
            fparam = first1.min(first2);
            lparam = last1.max(last2);
        } else if curve_kind(c1.as_ref()) == curve_kind(c2.as_ref())
            && is_line_or_conic(curve_kind(c1.as_ref()))
        {
            // `cxx:1206-1211`.
            is_canonic = true;
            new_curve = c1;
        } else {
            // `cxx:1212-1231`.
            let tc1: Arc<dyn Curve> = Arc::new(GeomTrimmedCurveBasis::new(c1, first1, last1));
            let tc2: Arc<dyn Curve> = Arc::new(GeomTrimmedCurveBasis::new(c2, first2, last2));
            let mut concat =
                CompCurveToBSplineCurve::from_curve(tc1.as_ref(), ParameterisationType::TgtThetaOver2)
                    .ok()?;
            // `cxx:1217-1220`: `!Add(...)` returns the null edge.
            if !concat
                .add(tc2.as_ref(), glue_tol, after, true, 0)
                .ok()?
            {
                return None;
            }
            let mut bspl: Arc<dyn Curve> = Arc::new(concat.into_curve()?);
            // `cxx:1221-1229`.
            if bspl.continuity() < GEOM_ABS_C1 {
                let approx =
                    GeomConvertApproxCurve::new(bspl.as_ref(), 1.0e-7, Shape::C1, 16, 14);
                if approx.has_result() {
                    if let Some(c) = approx.curve() {
                        bspl = Arc::new(c.clone());
                    }
                }
            }
            fparam = bspl.first_parameter();
            lparam = bspl.last_parameter();
            new_curve = bspl;
        }

        // `cxx:1233-1244`: the two `BRepLib_MakeEdge` arms.
        let build = if is_canonic {
            brep_lib_make_edge::init_3d_projecting(
                &new_curve,
                Some(vfirst),
                Some(vlast),
            )?
        } else {
            brep_lib_make_edge::init_3d_vertices(
                &new_curve,
                Some(vfirst),
                Some(vlast),
                fparam,
                lparam,
            )?
        };
        let mut new_edge = builder.make_edge(new_curve, build.first, build.last);
        brep_lib_make_edge::add_endpoint_vertices(&builder, &mut new_edge, &build);
        if build.closed {
            new_edge.0.set_closed(true);
        }

        // `cxx:1247-1257`.
        if add_pcurve1 {
            let (pc, new_first, new_last) = Self::concat_pcurves(e1, e2, f1, after)?;
            let key = GeometryRegistry::shape_key(&f1.0);
            let reg = GeometryRegistry::global();
            reg.set_edge_pcurve(&new_edge.0, key, pc);
            reg.set_pcurve_range(&new_edge.0, key, new_first, new_last);
        }
        if add_pcurve2 {
            let (pc, new_first, new_last) = Self::concat_pcurves(e1, e2, f2, after)?;
            let key = GeometryRegistry::shape_key(&f2.0);
            let reg = GeometryRegistry::global();
            reg.set_edge_pcurve(&new_edge.0, key, pc);
            reg.set_pcurve_range(&new_edge.0, key, new_first, new_last);
        }
        Some(new_edge)
    }

    /// `AssembleEdge(pDS, F1, F2, addPCurve1, addPCurve2, EdgesForConcat)`
    /// (`cxx:1345-1433`): glues the sequence into one edge, or answers the null
    /// edge when a joint is not autonom (`cxx:1380`, `cxx:1420`) or the glue
    /// fails (`cxx:1424`).
    ///
    /// The two null-vertex/edge-unset reads of OCCT (an empty
    /// `EdgesForConcat`, a `CommonVertex` that is null with `CV.IsNull()`
    /// handled, but the four `IsSame` branches assuming non-null `V11`..`V22`)
    /// become `None` here.
    pub fn assemble_edge(
        pds: &BopdsDS,
        f1: &Face,
        f2: &Face,
        add_pcurve1: bool,
        add_pcurve2: bool,
        edges_for_concat: &[Edge],
    ) -> Option<Edge> {
        // `cxx:1355`: `TopoDS::Edge(EdgesForConcat(1))`.
        let mut cur_edge = edges_for_concat.first()?.clone();
        // `cxx:1356`.
        let mut glue_tol = CONFUSION;

        // `cxx:1358-1432`.
        for an_edge in &edges_for_concat[1..] {
            let mut after = false;
            let (vfirst, vlast);
            // `cxx:1361-1419`.
            if Self::are_closed(&cur_edge, an_edge) {
                // `cxx:1363-1381`.
                let (v1, v2) = Self::edge_vertices(&cur_edge);
                let is_autonom_v1 = v1
                    .as_ref()
                    .map_or(false, |v| Self::is_autonom_vertex_on_faces(v, pds, f1, f2));
                let is_autonom_v2 = v2
                    .as_ref()
                    .map_or(false, |v| Self::is_autonom_vertex_on_faces(v, pds, f1, f2));
                if is_autonom_v1 {
                    after = false;
                    let v = v2.clone()?;
                    vfirst = v.clone();
                    vlast = v;
                } else if is_autonom_v2 {
                    after = true;
                    let v = v1.clone()?;
                    vfirst = v.clone();
                    vlast = v;
                } else {
                    return None;
                }
            } else {
                // `cxx:1384-1420`.
                let cv = crate::topo_tools::common_vertex(&cur_edge, an_edge)?;
                let cv = Vertex(cv);
                let is_autonom_cv = Self::is_autonom_vertex_on_faces(&cv, pds, f1, f2);
                if !is_autonom_cv {
                    return None;
                }
                // `cxx:1396`.
                glue_tol = BRepTool::vertex_tolerance(&cv);
                let (v11, v12) = Self::edge_vertices(&cur_edge);
                let (v21, v22) = Self::edge_vertices(an_edge);
                if is_same_vertex(&v11, &cv) && is_same_vertex(&v21, &cv) {
                    vfirst = v22.clone()?;
                    vlast = v12.clone()?;
                } else if is_same_vertex(&v11, &cv) && is_same_vertex(&v22, &cv) {
                    vfirst = v21.clone()?;
                    vlast = v12.clone()?;
                } else if is_same_vertex(&v12, &cv) && is_same_vertex(&v21, &cv) {
                    vfirst = v11.clone()?;
                    vlast = v22.clone()?;
                } else {
                    vfirst = v11.clone()?;
                    vlast = v21.clone()?;
                }
            }

            // `cxx:1422-1428`.
            let new_edge = Self::glue(
                &cur_edge,
                an_edge,
                &vfirst,
                &vlast,
                after,
                f1,
                add_pcurve1,
                f2,
                add_pcurve2,
                glue_tol,
            )?;
            cur_edge = new_edge;
        }
        // `cxx:1432`.
        Some(cur_edge)
    }
}
