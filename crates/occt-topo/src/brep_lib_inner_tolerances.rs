//! `BRepLib::UpdateInnerTolerances` (`BRepLib.cxx:1983-2083`).
//!
//! `TShape::Locked` (`BRep_Builder.cxx:1446-1448`, `:1002-1004`) is not stored,
//! so `UpdateVertex` / `UpdateEdge` always raise through `UpdateTolerance`.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpTrsf};
use occt_core::precision::epsilon;
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_geom2d::curve::Curve2d;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{edges_of, faces_of};

/// `BRepLib::UpdateInnerTolerances` (`BRepLib.cxx:1983-2083`).
pub fn update_inner_tolerances(shape: &TopoShape) {
    let reg = GeometryRegistry::global();
    for (edge, faces) in edge_ancestors(shape) {
        if !is_geometric(&edge) {
            continue;
        }
        let (v1, v2) = edge_ends(&edge);
        let (fpar, lpar) = reg.edge_parameters(&edge);
        let view = Edge(edge.clone());
        let same = BRepTool::same_parameter(&view);
        let degenerated = BRepTool::is_degenerated(&view);
        // A geometric edge with no loadable adaptor throws
        // `Standard_NullObject` (`BRepAdaptor_Curve.cxx:109`).
        let Some(center) = center_curve(&edge) else {
            continue;
        };
        let pnt1 = v1.as_ref().map(vertex_pnt);
        let pnt2 = v2.as_ref().map(vertex_pnt);

        if !degenerated && !faces.is_empty() {
            let reps: Vec<OnFace> = faces.iter().filter_map(|f| on_face(&edge, f)).collect();
            let nb = if same { 23 } else { 2 };
            let delta = (lpar - fpar) / (nb as f64 - 1.0);
            let mut max_dist = 0.0;
            for rep in &reps {
                for k in 0..=nb {
                    let param_on_center = if k == nb {
                        lpar
                    } else {
                        fpar + (k as f64) * delta
                    };
                    let center_pt = center.at(param_on_center);
                    let param_on_curve = if same {
                        param_on_center
                    } else if k == 0 {
                        rep.first
                    } else {
                        rep.last
                    };
                    let point = rep.at(param_on_curve);
                    let mut dist = center_pt.distance(&point);
                    dist += 2.0 * epsilon(dist);
                    if dist > max_dist {
                        max_dist = dist;
                    }
                    if k == 0 {
                        if let (Some(v), Some(p)) = (&v1, pnt1) {
                            let mut d1 = p.distance(&point);
                            d1 += 2.0 * epsilon(d1);
                            raise_vertex(v, d1);
                        }
                    }
                    if k == nb {
                        if let (Some(v), Some(p)) = (&v2, pnt2) {
                            let mut d2 = p.distance(&point);
                            d2 += 2.0 * epsilon(d2);
                            raise_vertex(v, d2);
                        }
                    }
                }
            }
            reg.raise_edge_tolerance(&edge, max_dist);
        }

        let tol_edge = reg.edge_tolerance(&edge);
        if let (Some(v), Some(p)) = (&v1, pnt1) {
            let end = center.at(fpar);
            let mut dist = p.distance(&end);
            dist += 2.0 * epsilon(dist);
            raise_vertex(v, dist.max(tol_edge));
        }
        if let (Some(v), Some(p)) = (&v2, pnt2) {
            let end = center.at(lpar);
            let mut dist = p.distance(&end);
            dist += 2.0 * epsilon(dist);
            raise_vertex(v, dist.max(tol_edge));
        }
    }
}

enum Center {
    Curve(Arc<dyn Curve>),
    OnSurf(OnFace),
}

impl Center {
    fn at(&self, t: f64) -> GpPnt {
        match self {
            Center::Curve(curve) => curve.value(t),
            Center::OnSurf(on) => on.at(t),
        }
    }
}

struct OnFace {
    curve: Arc<dyn Curve2d>,
    surface: Arc<dyn Surface>,
    trsf: GpTrsf,
    first: f64,
    last: f64,
}

impl OnFace {
    fn at(&self, t: f64) -> GpPnt {
        let uv = self.curve.value(t);
        self.surface
            .value(uv.x(), uv.y())
            .transformed(&self.trsf)
    }
}

/// `BRepAdaptor_Curve::Initialize(E)` (`BRepAdaptor_Curve.cxx:80-112`).
fn center_curve(edge: &TopoShape) -> Option<Center> {
    if let Some(curve) = BRepTool::edge_curve_world(&Edge(edge.clone())) {
        return Some(Center::Curve(curve));
    }
    // No 3D curve: the first CurveOnSurface, located by `E.Location() * GC->Location`
    // (`BRep_Tool.cxx:454-476`). The representation location is not stored.
    let (surface, curve, first, last) =
        GeometryRegistry::global().edge_pcurve_reps(edge).into_iter().next()?;
    Some(Center::OnSurf(OnFace {
        curve,
        surface,
        trsf: edge.location().transformation(),
        first,
        last,
    }))
}

/// `BRepAdaptor_Curve::Initialize(E, F)` (`BRepAdaptor_Curve.cxx:117-135`).
fn on_face(edge: &TopoShape, face: &TopoShape) -> Option<OnFace> {
    let (curve, first, last) =
        crate::boptools_2d::curve_on_surface_range(&Edge(edge.clone()), &Face(face.clone()))?;
    let surface = BRepTool::face_surface(&Face(face.clone()))?;
    Some(OnFace {
        curve,
        surface,
        trsf: face.location().transformation(),
        first,
        last,
    })
}

/// `BRep_Tool::IsGeometric(E)` (`BRep_Tool.cxx:238-261`).
fn is_geometric(edge: &TopoShape) -> bool {
    BRepTool::edge_curve(&Edge(edge.clone())).is_some()
        || !GeometryRegistry::global().edge_pcurve_reps(edge).is_empty()
}

/// `TopExp::MapShapesAndAncestors(S, EDGE, FACE)` (`BRepLib.cxx:1985-1986`).
fn edge_ancestors(shape: &TopoShape) -> Vec<(TopoShape, Vec<TopoShape>)> {
    let mut keys: Vec<TopoShape> = Vec::new();
    let mut faces: Vec<Vec<TopoShape>> = Vec::new();
    for face in faces_of(shape) {
        for edge in edges_of(&face.0) {
            if let Some(i) = keys.iter().position(|k| same_located(k, &edge.0)) {
                if !faces[i].iter().any(|f| same_located(f, &face.0)) {
                    faces[i].push(face.0.clone());
                }
            } else {
                keys.push(edge.0.clone());
                faces.push(vec![face.0.clone()]);
            }
        }
    }
    keys.into_iter().zip(faces).collect()
}

/// `TopTools_ShapeMapHasher` / `IsSame`: same `TShape` and the same location.
fn same_located(a: &TopoShape, b: &TopoShape) -> bool {
    crate::topo_tools_full::is_same(a, b)
}

/// `TopExp::Vertices(E, V1, V2)` with `CumOri = false` (`TopExp.hxx:126-129`).
/// Location is still cumulated (`TopoDS_Iterator` default `cumLoc`).
fn edge_ends(edge: &TopoShape) -> (Option<TopoShape>, Option<TopoShape>) {
    let stored = edge
        .tshape
        .read()
        .expect("poisoned TShape lock")
        .children
        .clone();
    let parent_loc = edge.location().clone();
    let move_loc = !parent_loc.is_identity();
    let mut first = None;
    let mut last = None;
    for mut child in stored {
        if child.shape_type() != ShapeType::Vertex {
            continue;
        }
        if move_loc {
            child.move_location(&parent_loc);
        }
        if child.orientation() == Orientation::Forward {
            first = Some(child);
        } else if child.orientation() == Orientation::Reversed {
            last = Some(child);
        }
    }
    (first, last)
}

/// `BRep_Tool::Pnt` (`BRep_Tool.cxx`): the vertex point moved by its location.
fn vertex_pnt(v: &TopoShape) -> GpPnt {
    let local = GeometryRegistry::global().vertex_point(v);
    local.transformed(&v.location().transformation())
}

/// `BRep_TVertex::UpdateTolerance` (`BRep_TVertex.lxx:33-37`): raise only.
fn raise_vertex(v: &TopoShape, tol: f64) {
    let mut ts = v.tshape.write().expect("poisoned TShape lock");
    let core = ts.vertex_core_mut();
    if tol > core.tolerance {
        core.tolerance = tol;
    }
}
