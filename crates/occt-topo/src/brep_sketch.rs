//! Phase 4 module: brep_sketch — 2D sketch → 3D profile → face.
//!
//! Port of the 2D-sketching front-end of `BRepBuilderAPI_MakeFace` /
//! `Geom2dAPI_PointsToBSpline`: a planar sketch (2D points + line segments)
//! embedded into a 3D plane, then turned into a wire or a planar face.

use std::sync::Arc;

use occt_core::gp::{GpDir, GpLin, GpPln, GpPnt, GpPnt2d, GpVec};
use occt_geom::GeomLine;

use crate::brep_builder_api;
use crate::builder::TopoBuilder;
use crate::shape::{Face, Vertex, Wire};

/// A 2D planar sketch: points plus line-segment connectivity.
#[derive(Debug, Clone, Default)]
pub struct Sketch {
    pub points: Vec<GpPnt2d>,
    pub lines: Vec<(usize, usize)>,
}

impl Sketch {
    pub fn new() -> Self {
        Self::default()
    }

    /// Append a 2D point and return its index.
    pub fn add_point(&mut self, p: GpPnt2d) -> usize {
        self.points.push(p);
        self.points.len() - 1
    }

    /// Append a line segment between two point indices.
    pub fn add_line(&mut self, a: usize, b: usize) {
        self.lines.push((a, b));
    }

    /// Number of points in the sketch.
    pub fn len(&self) -> usize {
        self.points.len()
    }
}

/// Map a 2D sketch point into the 3D plane: `loc + u·x_dir + v·y_dir`.
fn to_3d(p: &GpPnt2d, plane: &GpPln) -> GpPnt {
    let ax = plane.position();
    let loc = ax.location();
    let xd = ax.x_direction();
    let yd = ax.y_direction();
    GpPnt::new(
        loc.x() + p.x() * xd.x() + p.y() * yd.x(),
        loc.y() + p.x() * xd.y() + p.y() * yd.y(),
        loc.z() + p.x() * xd.z() + p.y() * yd.z(),
    )
}

/// Trace the closed boundary of a sketch starting at point 0. Returns an
/// ordered list of point indices around the cycle, or an empty list when the
/// line graph does not form a single closed cycle.
fn ordered_cycle(sketch: &Sketch) -> Vec<usize> {
    if sketch.points.is_empty() || sketch.lines.len() < 3 {
        return Vec::new();
    }
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); sketch.points.len()];
    for &(a, b) in &sketch.lines {
        if a >= sketch.points.len() || b >= sketch.points.len() {
            return Vec::new();
        }
        adj[a].push(b);
        adj[b].push(a);
    }
    let start = 0;
    let mut order = vec![start];
    let mut cur = start;
    let mut prev = usize::MAX;
    for _ in 0..=sketch.lines.len() {
        let mut next = None;
        for &nb in &adj[cur] {
            if nb != prev {
                next = Some(nb);
                break;
            }
        }
        match next {
            Some(n) if n == start => return order,
            Some(n) => {
                order.push(n);
                prev = cur;
                cur = n;
            }
            None => return Vec::new(),
        }
    }
    Vec::new()
}

/// Whether the sketch's line segments form a closed cycle.
fn lines_form_cycle(sketch: &Sketch) -> bool {
    !ordered_cycle(sketch).is_empty()
}

/// Build a wire from the sketch's line segments, embedded in `plane`.
/// Vertices are shared between adjacent edges (one per sketch point).
pub fn sketch_to_wire(sketch: &Sketch, plane: &GpPln) -> Result<Wire, String> {
    if sketch.lines.is_empty() {
        return Err("sketch_to_wire: sketch has no lines".into());
    }
    let pts3d: Vec<GpPnt> = sketch.points.iter().map(|p| to_3d(p, plane)).collect();
    let b = TopoBuilder::new();
    let verts: Vec<Vertex> = pts3d.iter().map(|p| b.make_vertex(*p, 0.0)).collect();
    let mut edges = Vec::with_capacity(sketch.lines.len());
    for &(i, j) in &sketch.lines {
        if i >= pts3d.len() || j >= pts3d.len() {
            return Err(format!("sketch_to_wire: line ({i},{j}) out of range"));
        }
        let dir = GpDir::from_vec(&GpVec::from_pnts(&pts3d[i], &pts3d[j]))
            .map_err(|_| "sketch_to_wire: degenerate line".to_string())?;
        let mut e = b.make_edge(
            Arc::new(GeomLine::new(GpLin::from_pnt_dir(pts3d[i], dir))),
            0.0,
            pts3d[i].distance(&pts3d[j]),
        );
        b.add(&mut e.0, &verts[i].0);
        b.add(&mut e.0, &verts[j].0);
        edges.push(e);
    }
    let wire = b.make_wire(&edges);
    wire.set_closed(lines_form_cycle(sketch));
    Ok(wire)
}

/// Closed sketch → planar face. The sketch points are taken in boundary order
/// (the closing edge back to the first point is added by the polygon face
/// builder).
pub fn sketch_to_face(sketch: &Sketch, plane: &GpPln) -> Result<Face, String> {
    let pts3d: Vec<GpPnt> = sketch.points.iter().map(|p| to_3d(p, plane)).collect();
    if pts3d.len() < 3 {
        return Err("sketch_to_face: need >= 3 points".into());
    }
    brep_builder_api::make_face_from_polygon(&pts3d)
}

/// Build a wire from an ordered list of 2D points embedded in `plane`.
/// `close` appends the closing edge back to the first point.
pub fn profile_from_points(points: &[GpPnt2d], plane: &GpPln, close: bool) -> Result<Wire, String> {
    let mut sketch = Sketch::new();
    for p in points {
        sketch.add_point(*p);
    }
    let n = points.len();
    for i in 0..n.saturating_sub(1) {
        sketch.add_line(i, i + 1);
    }
    if close && n >= 3 {
        sketch.add_line(n - 1, 0);
    }
    sketch_to_wire(&sketch, plane)
}

/// Axis-aligned 2D bounds of the sketch's points (empty sketch → the origin).
pub fn sketch_bbox2d(sketch: &Sketch) -> (GpPnt2d, GpPnt2d) {
    if sketch.points.is_empty() {
        return (GpPnt2d::zero(), GpPnt2d::zero());
    }
    let mut min = GpPnt2d::new(f64::INFINITY, f64::INFINITY);
    let mut max = GpPnt2d::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
    for p in &sketch.points {
        min = GpPnt2d::new(min.x().min(p.x()), min.y().min(p.y()));
        max = GpPnt2d::new(max.x().max(p.x()), max.y().max(p.y()));
    }
    (min, max)
}

/// Area of the closed boundary traced by the sketch's lines (shoelace
/// formula). The boundary is assumed to be a single closed cycle.
pub fn sketch_area2d(sketch: &Sketch) -> f64 {
    let order = ordered_cycle(sketch);
    if order.len() < 3 {
        return 0.0;
    }
    let n = order.len();
    let mut sum = 0.0;
    for i in 0..n {
        let a = sketch.points[order[i]];
        let b = sketch.points[order[(i + 1) % n]];
        sum += a.x() * b.y() - b.x() * a.y();
    }
    0.5 * sum.abs()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_surface::face_is_planar;
    use crate::topo_tools_full::{edges_of_wire, vertex_position, vertices_of};
    use occt_core::gp::GpAx3;

    fn standard_plane() -> GpPln {
        GpPln::new(GpAx3::standard())
    }

    fn square_sketch() -> Sketch {
        let mut s = Sketch::new();
        let a = s.add_point(GpPnt2d::new(0.0, 0.0));
        let b = s.add_point(GpPnt2d::new(1.0, 0.0));
        let c = s.add_point(GpPnt2d::new(1.0, 1.0));
        let d = s.add_point(GpPnt2d::new(0.0, 1.0));
        s.add_line(a, b);
        s.add_line(b, c);
        s.add_line(c, d);
        s.add_line(d, a);
        s
    }

    #[test]
    fn square_sketch_to_wire_and_face() {
        let sketch = square_sketch();
        assert_eq!(sketch.len(), 4);
        let plane = standard_plane();
        let wire = sketch_to_wire(&sketch, &plane).expect("wire");
        assert_eq!(edges_of_wire(&wire).len(), 4);
        assert!(wire.closed(), "square wire is flagged closed");
        // Standard plane → every vertex at z = 0.
        for v in vertices_of(&wire.0) {
            assert!((vertex_position(&v).z() - 0.0).abs() < 1e-12);
        }
        // Shoelace area and bbox of the closed boundary.
        assert!((sketch_area2d(&sketch) - 1.0).abs() < 1e-12);
        let (mn, mx) = sketch_bbox2d(&sketch);
        assert!((mn.x() - 0.0).abs() < 1e-12 && (mn.y() - 0.0).abs() < 1e-12);
        assert!((mx.x() - 1.0).abs() < 1e-12 && (mx.y() - 1.0).abs() < 1e-12);

        let face = sketch_to_face(&sketch, &plane).expect("face");
        assert!(face_is_planar(&face));
    }

    #[test]
    fn open_sketch_produces_open_wire() {
        let mut s = Sketch::new();
        for p in [
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(1.0, 0.0),
            GpPnt2d::new(2.0, 0.0),
            GpPnt2d::new(2.0, 1.0),
            GpPnt2d::new(1.0, 1.0),
        ] {
            s.add_point(p);
        }
        for i in 0..4 {
            s.add_line(i, i + 1);
        }
        let wire = sketch_to_wire(&s, &standard_plane()).expect("wire");
        assert!(!wire.closed(), "open wire is not closed");
        assert_eq!(edges_of_wire(&wire).len(), 4);
        assert_eq!(vertices_of(&wire.0).len(), 5);
    }

    #[test]
    fn profile_from_points_helpers() {
        let plane = standard_plane();
        let pts = [
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(1.0, 0.0),
            GpPnt2d::new(1.0, 1.0),
            GpPnt2d::new(0.0, 1.0),
        ];
        let open = profile_from_points(&pts, &plane, false).expect("open profile");
        assert!(!open.closed());
        assert_eq!(edges_of_wire(&open).len(), 3);
        let closed = profile_from_points(&pts, &plane, true).expect("closed profile");
        assert!(closed.closed());
        assert_eq!(edges_of_wire(&closed).len(), 4);
    }
}
