//! High-level shape construction — a port of the BRepBuilderAPI package.
//! Source: `BRepBuilderAPI_MakePolygon`, `MakeEdge`, `MakeFace`, `MakeSolid`,
//! `MakePrism` (TKTopAlgo).
//!
//! Thin convenience wrappers over `TopoBuilder` plus polygon/arc/primitive
//! helpers used by downstream CAD workflows: closed/open polygonal wires,
//! planar faces from polygons, circular arcs through three points, solid
//! assembly from faces, and full cylinder/cone primitives with their faces
//! exposed.

use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::gp::{GpAx3, GpCirc, GpDir, GpLin, GpPln, GpPnt, GpVec, GpXyz};
use occt_geom::{Curve, GeomCircle, GeomLine, GeomPlane, GeomTrimmedCurve};

use crate::builder::TopoBuilder;
use crate::shape::{Compound, Edge, Face, Solid, TopoShape, Vertex, Wire};

/// Result of a closed-polygon construction: the wire plus its boundary
/// vertices and edges (`BRepBuilderAPI_MakePolygon`).
#[derive(Debug, Clone)]
pub struct PolygonBuilder {
    pub wire: Wire,
    pub vertices: Vec<Vertex>,
    pub edges: Vec<Edge>,
}

/// Cylinder primitive with its three faces exposed (`BRepPrimAPI_MakeCylinder`).
#[derive(Debug, Clone)]
pub struct CylinderPrim {
    pub solid: Solid,
    pub axis: GpAx3,
    pub radius: f64,
    pub height: f64,
    pub bottom: Face,
    pub top: Face,
    pub lateral: Face,
}

/// Cone primitive with its two faces exposed (`BRepPrimAPI_MakeCone`).
#[derive(Debug, Clone)]
pub struct ConePrim {
    pub solid: Solid,
    pub radius: f64,
    pub height: f64,
    pub base: Face,
    pub lateral: Face,
}

/// Straight line curve through `p1` and `p2`, parameterized on `[0, |p2−p1|]`.
fn line_curve(p1: &GpPnt, p2: &GpPnt) -> Arc<dyn Curve> {
    let d = GpDir::from_vec(&GpVec::from_pnts(p1, p2))
        .unwrap_or_else(|_| GpDir::new(1.0, 0.0, 0.0).unwrap());
    Arc::new(GeomLine::new(GpLin::from_pnt_dir(*p1, d)))
}

/// Plane through three non-collinear points.
fn plane_through3(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> Result<GpPln, String> {
    let n = GpVec::from_pnts(a, b).xyz().crossed(GpVec::from_pnts(a, c).xyz());
    if n.modulus() < 1e-12 {
        return Err("plane_through3: collinear points".into());
    }
    let d = GpDir::from_xyz(&n).map_err(|_| "plane_through3: collinear points")?;
    let z_axis = GpDir::new(0.0, 0.0, 1.0).unwrap();
    let x_dir = if d.is_normal(&z_axis) { z_axis } else { GpDir::new(1.0, 0.0, 0.0).unwrap() };
    GpAx3::new(*a, d, &x_dir).map(GpPln::new).map_err(|e| e.to_string())
}

/// Indices of the first three non-collinear points in `points`.
fn first_three_noncollinear(points: &[GpPnt]) -> Result<(usize, usize, usize), String> {
    for i in 0..points.len() {
        for j in (i + 1)..points.len() {
            for k in (j + 1)..points.len() {
                let v1 = GpVec::from_pnts(&points[i], &points[j]);
                let v2 = GpVec::from_pnts(&points[i], &points[k]);
                if v1.xyz().crossed(v2.xyz()).modulus() > 1e-12 {
                    return Ok((i, j, k));
                }
            }
        }
    }
    Err("make_face_from_polygon: all points are collinear".into())
}

/// Shared edge builder for a polyline: registers vertex points and attaches the
/// shared corner vertices to each line edge. `closed` adds the closing edge.
fn build_wire_from_points(
    points: &[GpPnt],
    closed: bool,
) -> Result<(Wire, Vec<Vertex>, Vec<Edge>), String> {
    let min = if closed { 3 } else { 2 };
    if points.len() < min {
        return Err(format!("need at least {min} points, got {}", points.len()));
    }
    let b = TopoBuilder::new();
    let verts: Vec<Vertex> = points.iter().map(|p| b.make_vertex(*p, 0.0)).collect();
    let n = points.len();
    let segs: Vec<(usize, usize)> = if closed {
        (0..n).map(|i| (i, (i + 1) % n)).collect()
    } else {
        (0..n - 1).map(|i| (i, i + 1)).collect()
    };
    let mut edges = Vec::with_capacity(segs.len());
    for (i, j) in segs {
        let mut e = b.make_edge(line_curve(&points[i], &points[j]), 0.0, points[i].distance(&points[j]));
        b.add(&mut e.0, &verts[i].0);
        b.add(&mut e.0, &verts[j].0);
        edges.push(e);
    }
    let wire = b.make_wire(&edges);
    wire.set_closed(closed);
    Ok((wire, verts, edges))
}

/// Closed polygonal wire of straight segments through `points` (`BRepBuilderAPI_MakePolygon`).
/// Corner vertices are shared between adjacent edges; the closing edge returns
/// to the first point.
pub fn make_polygon(points: &[GpPnt]) -> Result<PolygonBuilder, String> {
    let (wire, vertices, edges) = build_wire_from_points(points, true)?;
    Ok(PolygonBuilder { wire, vertices, edges })
}

/// Open polygonal wire of straight segments through `points` (no closing edge).
pub fn make_wire_from_points(points: &[GpPnt]) -> Result<Wire, String> {
    let (wire, _, _) = build_wire_from_points(points, false)?;
    Ok(wire)
}

/// Planar face bounded by the closed polygon `points`. The plane passes through
/// the first three non-collinear points; the polygon must be non-degenerate.
pub fn make_face_from_polygon(points: &[GpPnt]) -> Result<Face, String> {
    let poly = make_polygon(points)?;
    if polygon_area(points) <= 1e-12 {
        return Err("make_face_from_polygon: degenerate polygon (zero area)".into());
    }
    let (i, j, k) = first_three_noncollinear(points)?;
    let pln = plane_through3(&points[i], &points[j], &points[k])?;
    // Coplanarity check: every vertex must lie on the computed plane.
    let axis = pln.axis();
    let n = *axis.direction().xyz();
    let d = pln.location().coord.dot(&n);
    if points.iter().any(|p| (p.coord.dot(&n) - d).abs() > 1e-6) {
        return Err("make_face_from_polygon: points are not coplanar".into());
    }
    let b = TopoBuilder::new();
    Ok(b.make_face(Arc::new(GeomPlane::new(pln)), &[poly.wire]))
}

/// 3×3 determinant (Cramer's-rule helper for the circumcenter solve).
fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Solve a 3×3 linear system via Cramer's rule.
fn solve3(a: &[[f64; 3]; 3], rhs: &[f64; 3]) -> Option<[f64; 3]> {
    let d = det3(a);
    if d.abs() < 1e-20 {
        return None;
    }
    let mut x = [0.0; 3];
    for k in 0..3 {
        let mut m = *a;
        for i in 0..3 {
            m[i][k] = rhs[i];
        }
        x[k] = det3(&m) / d;
    }
    Some(x)
}

/// Circumcenter of three non-collinear 3D points, if it exists.
///
/// Solves the perpendicular-bisector system: the center O is equidistant from
/// a, b, c (two bisector equations) and lies in their plane.
fn circumcenter(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> Option<GpPnt> {
    let d1 = GpVec::from_pnts(a, b);
    let d2 = GpVec::from_pnts(a, c);
    let n = d1.xyz().crossed(d2.xyz());
    if n.modulus() < 1e-20 {
        return None;
    }
    let n2 = |p: &GpPnt| p.coord.dot(&p.coord);
    let mat = [
        [d1.xyz().x, d1.xyz().y, d1.xyz().z],
        [d2.xyz().x, d2.xyz().y, d2.xyz().z],
        [n.x, n.y, n.z],
    ];
    let rhs = [0.5 * (n2(b) - n2(a)), 0.5 * (n2(c) - n2(a)), a.coord.dot(&n)];
    let o = solve3(&mat, &rhs)?;
    Some(GpPnt::new(o[0], o[1], o[2]))
}

/// Angle of `p` around `ax`'s Z axis, measured CCW from the X direction, in
/// `[0, 2π)`.
fn angle_in_frame(ax: &GpAx3, p: &GpPnt) -> f64 {
    let d = GpVec::from_pnts(&ax.location(), p);
    let x = d.dot(&GpVec::from_xyz(ax.x_direction().xyz()));
    let y = d.dot(&GpVec::from_xyz(ax.y_direction().xyz()));
    let a = y.atan2(x);
    if a < 0.0 { a + 2.0 * PI } else { a }
}

/// Circular arc through three points (`BRepBuilderAPI_MakeEdge(P1, P2, P3)`).
/// The curve is a `GeomCircle` trimmed to the angular range that sweeps the
/// minor arc from `p1` to `p3` through `p2`; endpoint vertex children are
/// attached. Errors on collinear points.
pub fn make_edge_arc(p1: &GpPnt, p2: &GpPnt, p3: &GpPnt) -> Result<Edge, String> {
    let center = circumcenter(p1, p2, p3).ok_or("make_edge_arc: points are collinear")?;
    let radius = p1.distance(&center);
    // Circle-plane normal, flipped if needed so the arc through p2 is CCW.
    let n = GpVec::from_pnts(p1, p2).xyz().crossed(GpVec::from_pnts(p1, p3).xyz());
    let mut normal = GpDir::from_xyz(&n).map_err(|_| "make_edge_arc: collinear points")?;
    // X axis points at p1, so p1 sits at angle 0 on the circle.
    let x_dir = GpDir::from_vec(&GpVec::from_pnts(&center, p1))
        .map_err(|_| "make_edge_arc: p1 coincides with the circumcenter")?;
    let mut ax = GpAx3::new(center, normal, &x_dir).map_err(|e| e.to_string())?;
    let a2 = angle_in_frame(&ax, p2);
    let mut a3 = angle_in_frame(&ax, p3);
    if a2 > a3 {
        // The arc through p2 goes clockwise in this frame; flip the normal so
        // it becomes the CCW arc from p1 (angle 0) to p3 (angle a3).
        normal = normal.reversed();
        ax = GpAx3::new(center, normal, &x_dir).map_err(|e| e.to_string())?;
        a3 = angle_in_frame(&ax, p3);
    }
    let circle: Arc<dyn Curve> = Arc::new(GeomCircle::new(GpCirc::new(ax.ax2(), radius)));
    let curve: Arc<dyn Curve> = Arc::new(GeomTrimmedCurve::new(circle, 0.0, a3));
    let b = TopoBuilder::new();
    let mut e = b.make_edge(curve, 0.0, 1.0); // trimmed curve normalizes to [0, 1]
    let v1 = b.make_vertex(*p1, 0.0);
    let v3 = b.make_vertex(*p3, 0.0);
    b.add(&mut e.0, &v1.0);
    b.add(&mut e.0, &v3.0);
    Ok(e)
}

/// Convenience circular-arc edge on the circle of `radius` centered at
/// `center` in the plane with the given `normal`, over the angular range
/// `[a0, a1]` (radians, measured CCW from an auto-selected in-plane X axis).
pub fn make_edge_circle_arc(center: &GpPnt, normal: &GpDir, radius: f64, a0: f64, a1: f64) -> Edge {
    let z_axis = GpDir::new(0.0, 0.0, 1.0).unwrap();
    let x_dir = if normal.is_normal(&z_axis) { z_axis } else { GpDir::new(1.0, 0.0, 0.0).unwrap() };
    let ax = GpAx3::new(*center, *normal, &x_dir).expect("make_edge_circle_arc: invalid frame");
    let circle: Arc<dyn Curve> = Arc::new(GeomCircle::new(GpCirc::new(ax.ax2(), radius)));
    let curve: Arc<dyn Curve> = Arc::new(GeomTrimmedCurve::new(circle, a0, a1));
    TopoBuilder::new().make_edge(curve, 0.0, 1.0)
}

/// Solid from a non-empty set of faces: one shell containing them, then a solid
/// (`BRepBuilderAPI_MakeSolid`).
pub fn make_solid_from_faces(faces: &[Face]) -> Result<Solid, String> {
    if faces.is_empty() {
        return Err("make_solid_from_faces: no faces".into());
    }
    let b = TopoBuilder::new();
    let shell = b.make_shell(faces);
    Ok(b.make_solid(&[shell]))
}

/// Full cylinder primitive (`BRepPrimAPI_MakeCylinder`), also exposing the
/// bottom, top and lateral faces.
pub fn make_cylinder_full(radius: f64, height: f64) -> Result<CylinderPrim, String> {
    if radius <= 0.0 || height <= 0.0 {
        return Err("make_cylinder_full: radius and height must be positive".into());
    }
    let prim = crate::primitives::BRepPrimCylinder::make_cylinder(radius, height);
    let faces = crate::topo_tools_full::faces_of(&prim.solid.0);
    if faces.len() != 3 {
        return Err("make_cylinder_full: expected 3 faces".into());
    }
    let axis = GpAx3::new(
        GpPnt::zero(),
        GpDir::new(0.0, 0.0, 1.0).unwrap(),
        &GpDir::new(1.0, 0.0, 0.0).unwrap(),
    )
    .map_err(|e| e.to_string())?;
    Ok(CylinderPrim {
        solid: prim.solid,
        axis,
        radius,
        height,
        bottom: faces[0].clone(),
        top: faces[1].clone(),
        lateral: faces[2].clone(),
    })
}

/// Full cone primitive (`BRepPrimAPI_MakeCone`), also exposing the base and
/// lateral faces.
pub fn make_cone_full(radius: f64, height: f64) -> Result<ConePrim, String> {
    if radius <= 0.0 || height <= 0.0 {
        return Err("make_cone_full: radius and height must be positive".into());
    }
    let prim = crate::primitives::BRepPrimCone::make_cone(radius, height);
    let faces = crate::topo_tools_full::faces_of(&prim.solid.0);
    if faces.len() != 2 {
        return Err("make_cone_full: expected 2 faces".into());
    }
    Ok(ConePrim {
        solid: prim.solid,
        radius,
        height,
        base: faces[0].clone(),
        lateral: faces[1].clone(),
    })
}

/// Prism extruding a planar face by `height` along `dir` (`BRepPrimAPI_MakePrism`).
/// `dir` is a direction; the full displacement is `dir` scaled to `height`
/// (`prism_from_face` expects the full displacement vector).
pub fn make_prism(face: &Face, dir: &GpVec, height: f64) -> Result<crate::sweep::Prism, String> {
    let mag = dir.xyz().modulus();
    if mag < 1e-30 {
        return Err("make_prism: zero direction".into());
    }
    if height <= 0.0 {
        return Err("make_prism: height must be positive".into());
    }
    let d = GpVec::new(dir.x() / mag * height, dir.y() / mag * height, dir.z() / mag * height);
    Ok(crate::sweep::prism_from_face(face, &d))
}

/// Compound containing `shapes` (`BRepBuilderAPI_MakeCompound`).
pub fn make_compound_shapes(shapes: &[TopoShape]) -> Compound {
    TopoBuilder::new().make_compound_of(shapes)
}

/// Signed planar polygon area from the cross-product shoelace formula
/// (`½ |Σ Pi × Pi+1|`), used by MakeFace degenerate-polygon validation.
pub fn polygon_area(points: &[GpPnt]) -> f64 {
    let n = points.len();
    if n < 3 {
        return 0.0;
    }
    let mut acc = GpXyz::zero();
    for i in 0..n {
        let j = (i + 1) % n;
        acc = acc.added(&points[i].coord.crossed(&points[j].coord));
    }
    0.5 * acc.modulus()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_tool::BRepTool;
    use crate::brep_surface::face_is_planar;
    use crate::topo_tools_full::{vertex_position, wire_is_closed};

    fn approx(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9 * b.abs().max(1.0)
    }

    #[test]
    fn polygon_square() {
        let pts = [
            GpPnt::new(0., 0., 0.),
            GpPnt::new(1., 0., 0.),
            GpPnt::new(1., 1., 0.),
            GpPnt::new(0., 1., 0.),
        ];
        let poly = make_polygon(&pts).expect("square polygon");
        assert_eq!(poly.vertices.len(), 4);
        assert_eq!(poly.edges.len(), 4);
        assert!(poly.wire.closed(), "polygon wire is flagged closed");
        assert!(wire_is_closed(&poly.wire), "polygon wire chains end-to-end");
        for (v, p) in poly.vertices.iter().zip(&pts) {
            assert!(vertex_position(v).is_equal(p), "vertex point matches input");
        }
        let face = make_face_from_polygon(&pts).expect("square face");
        assert!(face_is_planar(&face), "polygon face is planar");
    }

    #[test]
    fn edge_arc_through_three_points() {
        let p1 = GpPnt::new(1., 0., 0.);
        let p2 = GpPnt::new(0., 1., 0.);
        let p3 = GpPnt::new(-1., 0., 0.);
        let e = make_edge_arc(&p1, &p2, &p3).expect("arc edge");
        let curve = BRepTool::edge_curve(&e).expect("edge curve");
        // All samples sit on the unit circle centered at the origin.
        for t in [0.0, 0.25, 0.5, 0.75, 1.0] {
            let pt = curve.d0(t);
            assert!((pt.distance(&GpPnt::zero()) - 1.0).abs() < 1e-9, "radius at t={t}");
        }
        assert!(curve.d0(0.0).distance(&p1) < 1e-9, "start at p1");
        assert!(curve.d0(1.0).distance(&p3) < 1e-9, "end at p3");
        assert!(curve.d0(0.5).distance(&p2) < 1e-9, "midpoint at p2");
        // Endpoint vertex children match p1/p3.
        let (a, z) = crate::topo_tools_full::edge_vertices(&e);
        assert!(vertex_position(&a.unwrap()).distance(&p1) < 1e-9);
        assert!(vertex_position(&z.unwrap()).distance(&p3) < 1e-9);
    }

    #[test]
    fn cylinder_full_three_faces() {
        let c = make_cylinder_full(2.0, 3.0).expect("cylinder");
        let faces = crate::topo_tools_full::faces_of(&c.solid.0);
        assert_eq!(faces.len(), 3);
        assert!(c.bottom.is_face() && c.top.is_face() && c.lateral.is_face());
        assert!(approx(c.radius, 2.0) && approx(c.height, 3.0));
        assert!(approx(PI * c.radius * c.radius * c.height, PI * 4.0 * 3.0));
        // Caps are planar; the lateral surface is not.
        assert!(face_is_planar(&c.bottom) && face_is_planar(&c.top));
        assert!(!face_is_planar(&c.lateral));
    }

    #[test]
    fn prism_of_square_face() {
        let pts = [
            GpPnt::new(0., 0., 0.),
            GpPnt::new(1., 0., 0.),
            GpPnt::new(1., 1., 0.),
            GpPnt::new(0., 1., 0.),
        ];
        let face = make_face_from_polygon(&pts).expect("base face");
        let prism = make_prism(&face, &GpVec::new(0., 0., 3.), 3.0).expect("prism");
        assert_eq!(crate::sweep::prism_counts(&prism), (8, 12, 6));
        let vol = crate::sweep::prism_volume(&prism);
        assert!(approx(vol, 3.0), "volume = base area 1 × height 3, got {vol}");
        // A unit direction scales identically.
        let prism2 = make_prism(&face, &GpVec::new(0., 0., 1.), 3.0).expect("prism unit dir");
        assert!(approx(crate::sweep::prism_volume(&prism2), 3.0));
    }

    #[test]
    fn polygon_area_square_and_triangle() {
        let sq = [
            GpPnt::new(0., 0., 0.),
            GpPnt::new(1., 0., 0.),
            GpPnt::new(1., 1., 0.),
            GpPnt::new(0., 1., 0.),
        ];
        assert!(approx(polygon_area(&sq), 1.0));
        let tri = [GpPnt::new(0., 0., 0.), GpPnt::new(2., 0., 0.), GpPnt::new(0., 2., 0.)];
        assert!(approx(polygon_area(&tri), 2.0));
    }

    #[test]
    fn wire_and_solid_helpers() {
        // Open wire of two segments.
        let w = make_wire_from_points(&[
            GpPnt::new(0., 0., 0.),
            GpPnt::new(1., 0., 0.),
            GpPnt::new(1., 1., 0.),
        ])
        .expect("open wire");
        assert!(!w.closed());
        assert_eq!(crate::topo_tools_full::edges_of_wire(&w).len(), 2);

        // Solid from faces, compound from shapes.
        let faces = crate::topo_tools_full::faces_of(&make_cylinder_full(1.0, 1.0).unwrap().solid.0);
        let solid = make_solid_from_faces(&faces).expect("solid");
        assert!(solid.is_solid());
        assert!(make_solid_from_faces(&[]).is_err());

        let comp = make_compound_shapes(&[faces[0].0.clone(), faces[1].0.clone()]);
        assert!(comp.is_compound());

        // Convenience arc edge: endpoints on the circle.
        let e = make_edge_circle_arc(&GpPnt::zero(), &GpDir::new(0.0, 0.0, 1.0).unwrap(), 1.0, 0.0, PI);
        let curve = BRepTool::edge_curve(&e).expect("curve");
        assert!(curve.d0(0.0).distance(&GpPnt::new(1.0, 0.0, 0.0)) < 1e-9);
        assert!(curve.d0(1.0).distance(&GpPnt::new(-1.0, 0.0, 0.0)) < 1e-9);
    }
}
