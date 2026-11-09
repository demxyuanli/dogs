//! `BRepOffset_SimpleOffset` (`BRepOffset_SimpleOffset.cxx:1-427`,
//! `BRepOffset_SimpleOffset.hxx:44-190`): the topology-preserving offset as a
//! [`Modification`], i.e. without intersections.
//!
//! The mapping the header describes (`hxx:37-43`) is built once, in the
//! constructor, by [`FillOffsetData`] (`cxx:168-202`):
//!
//! * `FillFaceData` (`cxx:206-233`) collapses the B-spline singularities of the
//!   face (`BRepOffset::CollapseSingularities`) and replaces the surface by
//!   `BRepOffset::Surface(A, +-myOffsetValue, status, true)`; the face
//!   orientation picks the sign (`cxx:221-225`).
//! * `FillEdgeData` (`cxx:237-314`) re-approximates the pcurve of the edge on
//!   the first adjacent face's offset surface: `BRepBuilderAPI_MakeEdge(C2d, S,
//!   f, l)` (`cxx:267-268`), then `BRepLib::BuildCurves3d(E, MaxTolerance(E,
//!   VERTEX))` (`cxx:272-273`), then the maximum of
//!   `BRepLib_ValidateEdge(3d curve, pcurve on each adjacent offset surface,
//!   sameParameter = true)` (`cxx:279-309`).
//! * `FillVertexData` (`cxx:318-427`) averages the ends of the offset curves of
//!   every adjacent edge (`cxx:346-396`), closes the closed-edge case with a
//!   second sample (`cxx:387-393`), and takes `max(edge tolerance, spread) *
//!   1.001` as the new tolerance (`cxx:407-423`).
//!
//! Port notes:
//!
//! * `NCollection_IndexedDataMap<Shape, List<Shape>>` is [`AncestorMap`]: key
//!   insertion order follows `TopExp::MapShapesAndAncestors`
//!   (`TopExp.cxx:80-110`) and indices here are 0-based (the reference's are
//!   1-based).
//! * `BRep_Tool::CurveOnSurface(E, F, f, l)` is
//!   [`crate::boptools_2d::curve_on_surface_range`].
//! * `Adaptor3d_CurveOnSurface(Geom2dAdaptor_Curve(C2d, f, l),
//!   GeomAdaptor_Surface(S))` (`cxx:292-296`) is
//!   [`crate::meshing::edge_discret::CurveOnSurface::new`], which carries the
//!   same `[f, l]` window; `GeomAdaptor_Curve(C3d, f, l)` (`cxx:299`) is the
//!   built 3D curve itself (its own range is that window).
//! * `BRepLib::BuildCurves3d(E, Tol)` (`cxx:273`) is
//!   [`crate::brep_lib_build_curves3d::build_curve3d`] with the header defaults
//!   `GeomAbs_C1`, degree 14, `MaxSegment` 0.
//! * The pcurve of `cxx:267` is stored on the new edge through
//!   `BRep_Builder::UpdateEdge(E, C2d, S, L, Tol)`
//!   ([`GeometryRegistry::update_edge_pcurve_on_surface`]) so the 3D build has
//!   the same input the reference's curve-on-surface representation holds.
//!   The 3D curve slot the edge was created with is cleared first
//!   ([`GeometryRegistry::remove_edge_curve3d`]), because `BuildCurve3d`
//!   returns immediately when a 3D curve is already present
//!   (`BRepLib.cxx:321-325`).
//! * `TopLoc_Location` (`NewFaceData::myL`, always the null location,
//!   `cxx:229`) has no counterpart: the port keeps geometry in the shape's own
//!   frame, so the reported location is the identity.
//! * `BRep_Tool::MaxTolerance` (`cxx:272`) is
//!   [`crate::brep_tool::BRepTool::max_tolerance`].
//! * UNPORTED: none.

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_core::kernel::geomabs::Shape as GeomAbsShape;
use occt_core::precision::CONFUSION;
use occt_core::toploc::TopLocLocation;
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_geom2d::curve::Curve2d;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_lib_build_curves3d::build_curve3d;
use crate::brep_lib_make_edge::{add_endpoint_vertices, init_2d};
use crate::brep_lib_validate_edge::ValidateEdge;
use crate::brep_offset::{brep_offset_surface, collapse_singularities};
use crate::brep_tool::BRepTool;
use crate::brep_tools_modifier::{
    geom_abs_shape_from_u8, Modification, NewCurveInfo, NewPointInfo, NewParameterInfo,
    NewSurfaceInfo,
};
use crate::builder::TopoBuilder;
use crate::meshing::edge_discret::CurveOnSurface;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{edge_vertices, explorer_shapes, shape_map_key, AncestorMap, ShapeMapKey};

/// `NewFaceData` (`BRepOffset_SimpleOffset.hxx:147-154`).
#[derive(Clone)]
struct NewFaceData {
    /// `myOffsetS`.
    offset_s: Arc<dyn Surface>,
    /// `myL`: always the null location (`cxx:229`).
    l: TopLocLocation,
    /// `myTol`.
    tol: f64,
    /// `myRevWires`.
    rev_wires: bool,
    /// `myRevFace`.
    rev_face: bool,
}

/// `NewEdgeData` (`BRepOffset_SimpleOffset.hxx:156-161`).
#[derive(Clone)]
struct NewEdgeData {
    /// `myOffsetC`: `None` is OCCT's null handle, which
    /// `cxx:281` uses to stop the deviation loop.
    offset_c: Option<Arc<dyn Curve>>,
    /// `myL`.
    l: TopLocLocation,
    /// `myTol`.
    tol: f64,
}

/// `NewVertexData` (`BRepOffset_SimpleOffset.hxx:163-167`).
#[derive(Clone, Copy)]
struct NewVertexData {
    /// `myP`.
    p: GpPnt,
    /// `myTol`.
    tol: f64,
}

/// `BRepOffset_SimpleOffset` (`BRepOffset_SimpleOffset.hxx:44-190`).
pub struct BRepOffsetSimpleOffset {
    /// `myFaceInfo`.
    face_info: HashMap<ShapeMapKey, NewFaceData>,
    /// `myEdgeInfo`.
    edge_info: HashMap<ShapeMapKey, NewEdgeData>,
    /// `myVertexInfo`.
    vertex_info: HashMap<ShapeMapKey, NewVertexData>,
    /// `myOffsetValue`.
    offset_value: f64,
    /// `myTolerance`.
    tolerance: f64,
}

impl BRepOffsetSimpleOffset {
    /// `BRepOffset_SimpleOffset(theInputShape, theOffsetValue, theTolerance)`
    /// (`cxx:40-47`).
    pub fn new(input_shape: &TopoShape, offset_value: f64, tolerance: f64) -> Self {
        let mut this = Self {
            face_info: HashMap::new(),
            edge_info: HashMap::new(),
            vertex_info: HashMap::new(),
            offset_value,
            tolerance,
        };
        this.fill_offset_data(input_shape);
        this
    }

    /// `FillOffsetData(theShape)` (`cxx:168-202`).
    fn fill_offset_data(&mut self, shape: &TopoShape) {
        // `cxx:170-173`: clears old data.
        self.face_info.clear();
        self.edge_info.clear();
        self.vertex_info.clear();

        // `cxx:175-181`: the faces.
        for face_shape in explorer_shapes(shape, ShapeType::Face) {
            if let Some(face) = Face::wrap(face_shape) {
                self.fill_face_data(&face);
            }
        }

        // `cxx:183-191`: the edges with their adjacent faces.
        let edge_face_map =
            AncestorMap::map_shapes_and_ancestors(shape, ShapeType::Edge, ShapeType::Face);
        for idx in 0..edge_face_map.length() {
            let Some(edge) = Edge::wrap(edge_face_map.find_key(idx).clone()) else {
                continue;
            };
            let faces = edge_face_map.ancestors(idx).to_vec();
            self.fill_edge_data(&edge, &faces);
        }

        // `cxx:193-201`: the vertices with their adjacent edges.
        let vertex_edge_map =
            AncestorMap::map_shapes_and_ancestors(shape, ShapeType::Vertex, ShapeType::Edge);
        for idx in 0..vertex_edge_map.length() {
            let Some(vertex) = Vertex::wrap(vertex_edge_map.find_key(idx).clone()) else {
                continue;
            };
            let edges = vertex_edge_map.ancestors(idx).to_vec();
            self.fill_vertex_data(&vertex, &edges);
        }
    }

    /// `FillFaceData(theFace)` (`cxx:206-233`).
    fn fill_face_data(&mut self, face: &Face) {
        // `cxx:208-211`.
        let tol = BRepTool::face_tolerance(face);

        // `cxx:217`: `BRep_Tool::Surface(theFace)`; a null surface has no
        // offset, and the reference would dereference it, so the face is
        // skipped here.
        let Some(a_s) = BRepTool::face_surface(face) else {
            return;
        };
        // `cxx:218`: `BRepOffset::CollapseSingularities(aS, theFace, myTolerance)`.
        let a_s = collapse_singularities(&a_s, face, self.tolerance);

        // `cxx:221-225`: take the face orientation into account.
        let a_mult = if face.0.orientation() == Orientation::Reversed {
            -1.0
        } else {
            1.0
        };

        // `cxx:228`: `BRepOffset::Surface(aS, aMult * myOffsetValue, aStatus, true)`.
        let Ok((offset, _status)) = brep_offset_surface(a_s.as_ref(), a_mult * self.offset_value, true)
        else {
            return;
        };
        let nfd = NewFaceData {
            offset_s: Arc::from(offset),
            // `cxx:229`: `aNFD.myL = TopLoc_Location();` (null transformation).
            l: TopLocLocation::identity(),
            tol,
            // `cxx:209-210`.
            rev_wires: false,
            rev_face: false,
        };
        // `cxx:232`: `myFaceInfo.Bind(theFace, aNFD)`.
        self.face_info.insert(shape_map_key(&face.0), nfd);
    }

    /// `FillEdgeData(theEdge, theEdgeFaceMap, theIdx)` (`cxx:237-314`).
    fn fill_edge_data(&mut self, edge: &Edge, faces: &[TopoShape]) {
        // `cxx:246-249`: free edges are skipped.
        if faces.is_empty() {
            return;
        }

        // `cxx:252`: `aCurrFace = TopoDS::Face(aFacesList.First())`.
        let Some(curr_face) = Face::wrap(faces[0].clone()) else {
            return;
        };

        // `cxx:254-257`.
        let Some(nfd) = self.face_info.get(&shape_map_key(&curr_face.0)).cloned() else {
            return;
        };

        // `cxx:259-261`: the offset surface of that face.
        let offset_surf = nfd.offset_s;

        // `cxx:264-265`: `BRep_Tool::CurveOnSurface(theEdge, aCurrFace, aF, aL)`.
        let Some((pc, a_f, a_l)) = curve_on_surface_range(edge, &curr_face) else {
            return;
        };

        // `cxx:267-268`: `BRepBuilderAPI_MakeEdge(aC2d, anOffsetSurf, aF, aL)`
        // and `anEdgeMaker.Edge()`; `!IsDone()` is the reference's null edge.
        let Some(build) = init_2d(&pc, &offset_surf, a_f, a_l) else {
            return;
        };
        let builder = TopoBuilder::new();
        // The edge is created with the curve-on-surface view as its 3D slot;
        // the slot is cleared right after so `BuildCurve3d` runs the build.
        let cos: Arc<dyn Curve> = Arc::new(CurveOnSurface::new(
            pc.clone(),
            offset_surf.clone(),
            build.first,
            build.last,
        ));
        let mut new_edge = builder.make_edge(cos, build.first, build.last);
        add_endpoint_vertices(&builder, &mut new_edge, &build);

        let reg = GeometryRegistry::global();
        // `BRep_Builder::UpdateEdge(E, C2d, S, L, Tol)`: the reference stores
        // the pcurve in the edge's curve-on-surface representation.
        reg.update_edge_pcurve_on_surface(&new_edge.0, &offset_surf, Some(pc.clone()), CONFUSION);
        reg.remove_edge_curve3d(&new_edge.0);

        // `cxx:272-273`: `aTol = BRep_Tool::MaxTolerance(theEdge, TopAbs_VERTEX)`
        // and `BRepLib::BuildCurves3d(aNewEdge, aTol)`.
        let a_tol = BRepTool::max_tolerance(&edge.0, ShapeType::Vertex);
        let built = build_curve3d(
            &new_edge,
            a_tol,
            GeomAbsShape::C1,
            14,
            0,
        );
        let _ = built;

        // `cxx:275-276`: `anOffsetC = BRep_Tool::Curve(aNewEdge, L, aF, aL)`.
        let offset_c = reg.edge_curve(&new_edge.0);

        // `cxx:278-309`: the maximum deviation over the adjacent faces.
        let mut an_edge_tol: f64 = 0.0;
        if let Some(offset_c) = offset_c.as_ref() {
            for face_shape in faces {
                let Some(cur_face) = Face::wrap(face_shape.clone()) else {
                    continue;
                };
                let Some(cur_nfd) = self.face_info.get(&shape_map_key(&cur_face.0)).cloned() else {
                    continue;
                };
                // `cxx:291`: `BRep_Tool::CurveOnSurface(theEdge, aCurFace, aF, aL)`.
                let Some((a_c2d_new, f2, l2)) = curve_on_surface_range(edge, &cur_face) else {
                    return;
                };
                // `cxx:292-296`: the curve-on-surface adaptor on the offset
                // surface of that face.
                let cos_new: Arc<dyn Curve> = Arc::new(CurveOnSurface::new(
                    a_c2d_new,
                    cur_nfd.offset_s.clone(),
                    f2,
                    l2,
                ));
                // `cxx:299`: `GeomAdaptor_Curve(anOffsetC, aF, aL)`.
                let curve3d: Arc<dyn Curve> = offset_c.clone();
                // `cxx:302-308`. The exact arm of `BRepLib_ValidateEdge` is
                // never taken here (`SetExactMethod` is not called), so the
                // curve-on-surface 2D handle is not needed.
                let mut validate = ValidateEdge::new(
                    curve3d.as_ref(),
                    cos_new.as_ref(),
                    Some(cur_nfd.offset_s.as_ref()),
                    None,
                    true,
                );
                validate.process();
                if validate.is_done() {
                    let a_max_tol1 = validate.get_max_distance();
                    an_edge_tol = an_edge_tol.max(a_max_tol1);
                }
            }
        }
        // `cxx:310`: `myTol = max(BRep_Tool::Tolerance(aNewEdge), anEdgeTol)`.
        let tol = reg.edge_tolerance(&new_edge.0).max(an_edge_tol);

        // `cxx:313`: `myEdgeInfo.Bind(theEdge, aNED)`.
        self.edge_info.insert(
            shape_map_key(&edge.0),
            NewEdgeData {
                offset_c,
                l: TopLocLocation::identity(),
                tol,
            },
        );
    }

    /// `FillVertexData(theVertex, theVertexEdgeMap, theIdx)` (`cxx:318-427`).
    fn fill_vertex_data(&mut self, vertex: &Vertex, edges: &[TopoShape]) {
        // `cxx:331`: `gp_Pnt aCurrPnt = BRep_Tool::Pnt(theVertex)`.
        let curr_pnt = BRepTool::vertex_point_world(vertex);

        // `cxx:333-338`: free vertices are skipped.
        if edges.is_empty() {
            return;
        }

        // `cxx:341`: the offset points.
        let mut offset_points: Vec<GpPnt> = Vec::new();
        // `cxx:343`.
        let mut max_edge_tol: f64 = 0.0;

        // `cxx:346-396`.
        for edge_shape in edges {
            let Some(curr_edge) = Edge::wrap(edge_shape.clone()) else {
                continue;
            };
            // `cxx:351-354`.
            let Some(ned) = self.edge_info.get(&shape_map_key(&curr_edge.0)).cloned() else {
                continue;
            };

            // `cxx:357-364`: protection from degenerated edges.
            let reg = GeometryRegistry::global();
            let (a_f, a_l) = BRepTool::edge_parameters(&curr_edge);
            let Some(c3d) = reg.edge_curve(&curr_edge.0) else {
                continue;
            };

            // `cxx:366-370`.
            let pnt_f = c3d.value(a_f);
            let pnt_l = c3d.value(a_l);
            let sq_dist_f = pnt_f.square_distance(&curr_pnt);
            let sq_dist_l = pnt_l.square_distance(&curr_pnt);

            // `cxx:372-378`.
            let (min_param, max_param) = if sq_dist_l < sq_dist_f {
                (a_l, a_f)
            } else {
                (a_f, a_l)
            };

            // `cxx:381-384`.
            let Some(offset_curve) = ned.offset_c.as_ref() else {
                continue;
            };
            let offset_point = offset_curve.value(min_param);
            offset_points.push(offset_point);

            // `cxx:386-393`: a closed edge contributes both ends.
            let (v1, v2) = edge_vertices(&curr_edge);
            if let (Some(v1), Some(v2)) = (v1, v2) {
                if v1.0.same_tshape(&v2.0) {
                    offset_points.push(offset_curve.value(max_param));
                }
            }

            // `cxx:395`.
            max_edge_tol = max_edge_tol.max(ned.tol);
        }

        // `cxx:400-405`: the barycentre of the offset points.
        let n = offset_points.len() as f64;
        let mut center = GpPnt::zero();
        for p in &offset_points {
            let c = center.coord.add(&p.coord);
            center = GpPnt::from_xyz(&c);
        }
        center = GpPnt::from_xyz(&center.coord.divided(n));

        // `cxx:407-416`: the maximum distance from the barycentre.
        let mut sq_max_dist: f64 = 0.0;
        for p in &offset_points {
            let sq_dist = center.square_distance(p);
            if sq_dist > sq_max_dist {
                sq_max_dist = sq_dist;
            }
        }

        // `cxx:418`.
        let res_tol = max_edge_tol.max(sq_max_dist.sqrt());

        // `cxx:420-426`.
        const MULT_COEFF: f64 = 1.001; // Avoid tolerance problems.
        self.vertex_info.insert(
            shape_map_key(&vertex.0),
            NewVertexData { p: center, tol: res_tol * MULT_COEFF },
        );
    }
}

impl Modification for BRepOffsetSimpleOffset {
    /// `NewSurface(F, S, L, Tol, RevWires, RevFace)` (`cxx:51-72`).
    fn new_surface(&self, f: &TopoShape) -> Option<NewSurfaceInfo> {
        // `cxx:58-61`.
        let nfd = self.face_info.get(&shape_map_key(f))?;
        Some(NewSurfaceInfo {
            // `S = aNFD.myOffsetS` (`cxx:65`).
            surface: nfd.offset_s.clone(),
            // `L = aNFD.myL` (`cxx:66`).
            location: nfd.l.clone(),
            // `Tol = aNFD.myTol` (`cxx:67`).
            tolerance: nfd.tol,
            // `RevWires = aNFD.myRevWires` (`cxx:68`).
            rev_wires: nfd.rev_wires,
            // `RevFace = aNFD.myRevFace` (`cxx:69`).
            rev_face: nfd.rev_face,
        })
    }

    /// `NewCurve(E, C, L, Tol)` (`cxx:76-93`).
    fn new_curve(&self, e: &TopoShape) -> Option<NewCurveInfo> {
        // `cxx:81-84`.
        let ned = self.edge_info.get(&shape_map_key(e))?;
        Some(NewCurveInfo {
            // `C = aNED.myOffsetC` (`cxx:88`).
            curve: ned.offset_c.clone(),
            // `L = aNED.myL` (`cxx:89`).
            location: ned.l.clone(),
            // `Tol = aNED.myTol` (`cxx:90`).
            tolerance: ned.tol,
        })
    }

    /// `NewPoint(V, P, Tol)` (`cxx:97-110`).
    fn new_point(&self, v: &TopoShape) -> Option<NewPointInfo> {
        // `cxx:99-102`.
        let nvd = self.vertex_info.get(&shape_map_key(v))?;
        Some(NewPointInfo {
            // `P = aNVD.myP` (`cxx:106`).
            point: nvd.p,
            // `Tol = aNVD.myTol` (`cxx:107`).
            tolerance: nvd.tol,
        })
    }

    /// `NewCurve2d(E, F, NewE, NewF, C, Tol)` (`cxx:114-132`): always `true`.
    fn new_curve2d(
        &self,
        e: &TopoShape,
        f: &TopoShape,
        new_e: &TopoShape,
        new_f: &TopoShape,
    ) -> Option<Arc<dyn Curve2d>> {
        let _ = (new_e, new_f);
        let (Some(edge), Some(face)) = (Edge::wrap(e.clone()), Face::wrap(f.clone())) else {
            return None;
        };
        // `C = BRep_Tool::CurveOnSurface(E, F, aF, aL)` (`cxx:123`): "Use
        // original pcurve." `Tol = BRep_Tool::Tolerance(E)` (`cxx:124`) is the
        // three-argument `Tolerance(Edge)` clamp.
        let (c, ..) = curve_on_surface_range(&edge, &face)?;
        let _tol = self
            .edge_info
            .get(&shape_map_key(&edge.0))
            .map_or_else(|| BRepTool::edge_tolerance(&edge), |ned| ned.tol);
        Some(c)
    }

    /// `NewParameter(V, E, P, Tol)` (`cxx:136-151`): always `true`.
    fn new_parameter(&self, v: &TopoShape, e: &TopoShape) -> Option<NewParameterInfo> {
        let Some(vertex) = Vertex::wrap(v.clone()) else {
            return None;
        };
        // `P = BRep_Tool::Parameter(V, E)` (`cxx:142`).
        let p = BRepTool::parameter_of_vertex(e, v);
        // `Tol = BRep_Tool::Tolerance(V)` (`cxx:143`), replaced by the new
        // vertex tolerance when one was computed (`cxx:145-148`).
        let tol = self
            .vertex_info
            .get(&shape_map_key(&vertex.0))
            .map_or_else(|| BRepTool::vertex_tolerance(&vertex), |nvd| nvd.tol);
        Some(NewParameterInfo { parameter: p, tolerance: tol })
    }

    /// `Continuity(E, F1, F2, ...)` (`cxx:155-164`): the source continuity.
    fn continuity(
        &self,
        e: &TopoShape,
        f1: &TopoShape,
        f2: &TopoShape,
        new_e: &TopoShape,
        new_f1: &TopoShape,
        new_f2: &TopoShape,
    ) -> GeomAbsShape {
        let _ = (new_e, new_f1, new_f2);
        let (Some(edge), Some(face1), Some(face2)) =
            (Edge::wrap(e.clone()), Face::wrap(f1.clone()), Face::wrap(f2.clone()))
        else {
            return GeomAbsShape::C0;
        };
        geom_abs_shape_from_u8(crate::brep_lib_continuity::stored_continuity(
            &edge, &face1, &face2,
        ))
    }
}
