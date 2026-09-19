//! Phase 4 module: sweep_revolve — surfaces of revolution (lathe).
//! Source: `BRepPrimAPI_MakeRevol`, `BRepSweep_Revol`.
//!
//! Revolves a 2D profile (radius, height) in the (r, z) plane around the Z
//! axis. Each profile segment sweeps one face:
//! - horizontal segment → planar disc/annulus (`GpPln`);
//! - vertical segment   → cylindrical face (`GpCylinder`);
//! - diagonal segment   → conical face (`GpCone`, a frustum or a full cone);
//! - segments on the axis (r = 0) are degenerate and produce no face.
//!
//! Every ring around the axis has `steps` circular-arc edges; consecutive
//! faces share ring vertices, so a closed profile
//! `[(0,0),(r,0),(r,h),(0,h)]` sweeps a solid cylinder (bottom disc + lateral
//! cylinder + top disc) and `[(0,0),(r,h),(0,h)]` sweeps a cone.

use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::gp::{GpAx2, GpAx3, GpCone, GpCylinder, GpDir, GpLin, GpPln, GpPnt, GpPnt2d, GpVec};
use occt_geom::{GeomCone, GeomCylinder, GeomLine, GeomPlane, Surface};

use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, Solid, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;

/// Result of revolving a 2D profile around the Z axis.
#[derive(Debug, Clone)]
pub struct RevolvedSolid {
    pub solid: Solid,
    pub faces: Vec<Face>,
    pub rings: Vec<Vec<Vertex>>,
    pub profiles: Vec<Wire>,
}

const EPS: f64 = 1e-12;

fn z_dir() -> GpDir {
    GpDir::new(0.0, 0.0, 1.0).expect("sweep_revolve: Z axis")
}
fn x_dir() -> GpDir {
    GpDir::new(1.0, 0.0, 0.0).expect("sweep_revolve: X axis")
}

fn validate(profile: &[GpPnt2d], steps: usize, angle: f64) -> Result<(), String> {
    if profile.len() < 2 {
        return Err("sweep_revolve: profile needs at least 2 points".into());
    }
    if steps < 4 {
        return Err("sweep_revolve: steps must be >= 4".into());
    }
    if profile.iter().any(|p| p.x() < 0.0) {
        return Err("sweep_revolve: profile radii must be >= 0".into());
    }
    if !(angle > 0.0 && angle <= 2.0 * PI) {
        return Err("sweep_revolve: angle must be in (0, 2π]".into());
    }
    Ok(())
}

/// One ring of `steps` circular-arc edges at radius `r`, height `z`, spanning
/// the angular range `[0, total]`. Returns the ring vertices and edges.
///
/// A full revolution (total == 2π) closes the ring: the last vertex is the
/// first handle (`steps` distinct vertices). A partial angle leaves an open
/// arc of `steps + 1` distinct vertices. In both cases `verts[steps]` is the
/// angular end of the ring (== `verts[0]` for a closed ring), so the lateral
/// face wire can attach its seam edges uniformly.
fn build_ring(b: &TopoBuilder, r: f64, z: f64, steps: usize, total: f64) -> (Vec<Vertex>, Vec<Edge>) {
    let full = (total - 2.0 * PI).abs() < 1e-9;
    let mut verts: Vec<Vertex> = Vec::with_capacity(steps + 1);
    for i in 0..=steps {
        let theta = total * i as f64 / steps as f64;
        if full && i == steps {
            verts.push(verts[0].clone());
        } else {
            verts.push(b.make_vertex(GpPnt::new(r * theta.cos(), r * theta.sin(), z), 0.0));
        }
    }
    let ax2 = GpAx2::new(GpPnt::new(0.0, 0.0, z), z_dir(), x_dir()).expect("sweep_revolve: ring axis");
    let mut edges = Vec::with_capacity(steps);
    for j in 0..steps {
        let a0 = total * j as f64 / steps as f64;
        let a1 = total * (j + 1) as f64 / steps as f64;
        let mut e = b.make_edge_circle(&ax2, r, a0, a1);
        b.add_edge_vertices(&mut e, &verts[j], &verts[j + 1]);
        edges.push(e);
    }
    (verts, edges)
}

/// Line edge between two existing ring vertices (the seam of a lateral face).
fn edge_through(b: &TopoBuilder, v1: &Vertex, v2: &Vertex) -> Edge {
    let p1 = GeometryRegistry::global().vertex_point(&v1.0);
    let p2 = GeometryRegistry::global().vertex_point(&v2.0);
    let dir = GpDir::from_vec(&GpVec::from_pnts(&p1, &p2)).unwrap_or_else(|_| x_dir());
    let mut e = b.make_edge(Arc::new(GeomLine::new(GpLin::from_pnt_dir(p1, dir))), 0.0, p1.distance(&p2));
    b.add_edge_vertices(&mut e, v1, v2);
    e
}

/// Boundary wire of a cylindrical/conical face between two rings: lower ring
/// arcs, a seam up, upper ring arcs (conceptually traversed back), a seam
/// down. Mirrors OCCT's seam-edge convention for a single-face boundary.
fn lateral_wire(b: &TopoBuilder, ring1: &[Vertex], edges1: &[Edge], ring2: &[Vertex], edges2: &[Edge]) -> Wire {
    let last1 = ring1.len() - 1;
    let last2 = ring2.len() - 1;
    let seam_up = edge_through(b, &ring1[last1], &ring2[last2]);
    let seam_down = edge_through(b, &ring2[0], &ring1[0]);
    let mut edges = Vec::with_capacity(edges1.len() + edges2.len() + 2);
    edges.extend(edges1.iter().cloned());
    edges.push(seam_up);
    edges.extend(edges2.iter().rev().cloned());
    edges.push(seam_down);
    b.make_wire(&edges)
}

/// The analytic surface swept by a profile segment `(r1,z1)→(r2,z2)`.
fn segment_surface(r1: f64, z1: f64, r2: f64, z2: f64) -> Result<Arc<dyn Surface>, String> {
    let dr = r2 - r1;
    let dz = z2 - z1;
    if dz.abs() < EPS {
        // Horizontal segment → planar disc/annulus at height z1.
        let pln = GpPln::new(GpAx3::new(GpPnt::new(0.0, 0.0, z1), z_dir(), &x_dir())?);
        Ok(Arc::new(GeomPlane::new(pln)))
    } else if dr.abs() < EPS {
        // Vertical segment → cylinder of radius r1 between z1 and z2.
        let ax3 = GpAx3::new(GpPnt::new(0.0, 0.0, z1), z_dir(), &x_dir())?;
        Ok(Arc::new(GeomCylinder::new(GpCylinder::new(ax3, r1)?)))
    } else {
        // Diagonal segment → cone frustum (full cone when one radius is 0).
        // The radius is linear in z; the apex is where it vanishes.
        let semi_angle = (dr.abs() / dz.abs()).atan();
        let z_apex = z1 - r1 * dz / dr;
        let dir = if z_apex <= z1.min(z2) + 1e-9 {
            z_dir()
        } else {
            GpDir::new(0.0, 0.0, -1.0).expect("sweep_revolve: Z axis")
        };
        let ax3 = GpAx3::new(GpPnt::new(0.0, 0.0, z_apex), dir, &x_dir())?;
        Ok(Arc::new(GeomCone::new(GpCone::new(ax3, 0.0, semi_angle)?)))
    }
}

/// Boundary wire for a profile segment: the outer circle ring for a horizontal
/// segment, or the lateral panel wire (rings + seams) otherwise.
fn segment_wire(
    b: &TopoBuilder,
    ring1: &[Vertex],
    edges1: &[Edge],
    ring2: &[Vertex],
    edges2: &[Edge],
    r1: f64,
    z1: f64,
    r2: f64,
    z2: f64,
) -> Wire {
    if (z2 - z1).abs() < EPS {
        if r1 >= r2 {
            b.make_wire(edges1)
        } else {
            b.make_wire(edges2)
        }
    } else {
        lateral_wire(b, ring1, edges1, ring2, edges2)
    }
}

/// Build the face swept by one profile segment. Returns `None` for degenerate
/// segments (on the axis, or a single point).
fn build_segment_face(
    b: &TopoBuilder,
    ring1: &[Vertex],
    edges1: &[Edge],
    ring2: &[Vertex],
    edges2: &[Edge],
    r1: f64,
    z1: f64,
    r2: f64,
    z2: f64,
) -> Result<Option<(Face, Wire)>, String> {
    let dr = r2 - r1;
    let dz = z2 - z1;
    if (r1.abs() < EPS && r2.abs() < EPS) || (dr.abs() < EPS && dz.abs() < EPS) {
        return Ok(None);
    }
    let surface = segment_surface(r1, z1, r2, z2)?;
    let profile_wire;
    let mut wires: Vec<Wire> = Vec::new();
    if dz.abs() < EPS {
        // Horizontal: outer circle wire; when both radii are positive the inner
        // ring is added as a hole wire (planar annulus).
        let (outer_edges, inner_edges) = if r1 >= r2 { (edges1, edges2) } else { (edges2, edges1) };
        profile_wire = b.make_wire(outer_edges);
        wires.push(profile_wire.clone());
        if r1.min(r2) > EPS {
            wires.push(b.make_wire(inner_edges));
        }
    } else {
        profile_wire = lateral_wire(b, ring1, edges1, ring2, edges2);
        wires.push(profile_wire.clone());
    }
    let face = b.make_face(surface, &wires);
    Ok(Some((face, profile_wire)))
}

/// Revolve a 2D profile (x = radius ≥ 0, y = height) around the Z axis by a
/// full 360° revolution.
pub fn revolve_polyline_around_z(profile: &[GpPnt2d], steps: usize) -> Result<RevolvedSolid, String> {
    revolve_angle(profile, steps, 2.0 * PI)
}

/// Revolve a 2D profile around the Z axis by `angle` radians (0 < angle ≤ 2π).
/// The rings span the angular range `[0, angle]`.
pub fn revolve_angle(profile: &[GpPnt2d], steps: usize, angle: f64) -> Result<RevolvedSolid, String> {
    validate(profile, steps, angle)?;
    let b = TopoBuilder::new();
    let mut rings = Vec::with_capacity(profile.len());
    let mut ring_edges = Vec::with_capacity(profile.len());
    for p in profile {
        let (v, e) = build_ring(&b, p.x(), p.y(), steps, angle);
        rings.push(v);
        ring_edges.push(e);
    }
    let mut faces = Vec::new();
    let mut profiles = Vec::new();
    for i in 0..profile.len() - 1 {
        let (r1, z1) = (profile[i].x(), profile[i].y());
        let (r2, z2) = (profile[i + 1].x(), profile[i + 1].y());
        if let Some((face, wire)) = build_segment_face(
            &b,
            &rings[i],
            &ring_edges[i],
            &rings[i + 1],
            &ring_edges[i + 1],
            r1,
            z1,
            r2,
            z2,
        )? {
            faces.push(face);
            profiles.push(wire);
        }
    }
    let shell = b.make_shell(&faces);
    let solid = b.make_solid(&[shell]);
    Ok(RevolvedSolid { solid, faces, rings, profiles })
}

/// Revolve the closed rectangle profile `[(0,0),(r,0),(r,h),(0,h)]` → a solid
/// cylinder: bottom disc + lateral cylinder + top disc.
pub fn revolve_rectangle(radius: f64, height: f64, steps: usize) -> Result<RevolvedSolid, String> {
    if !(radius > 0.0 && height > 0.0) {
        return Err("sweep_revolve: revolve_rectangle needs positive radius and height".into());
    }
    let profile = [
        GpPnt2d::new(0.0, 0.0),
        GpPnt2d::new(radius, 0.0),
        GpPnt2d::new(radius, height),
        GpPnt2d::new(0.0, height),
    ];
    revolve_polyline_around_z(&profile, steps)
}

/// Revolve the triangle profile `[(0,0),(r,h),(0,h)]` → a cone (apex at the
/// origin): lateral conical surface + base disc.
pub fn revolve_triangle(radius: f64, height: f64, steps: usize) -> Result<RevolvedSolid, String> {
    if !(radius > 0.0 && height > 0.0) {
        return Err("sweep_revolve: revolve_triangle needs positive radius and height".into());
    }
    let profile = [
        GpPnt2d::new(0.0, 0.0),
        GpPnt2d::new(radius, height),
        GpPnt2d::new(0.0, height),
    ];
    revolve_polyline_around_z(&profile, steps)
}

/// Boundary wire of one profile segment (two points), helper for tests and
/// partial-angle work.
pub fn revolve_open_wire(points: &[GpPnt2d], steps: usize) -> Result<Wire, String> {
    if points.len() != 2 {
        return Err("sweep_revolve: revolve_open_wire expects exactly 2 profile points".into());
    }
    validate(points, steps, 2.0 * PI)?;
    let b = TopoBuilder::new();
    let (ring1, edges1) = build_ring(&b, points[0].x(), points[0].y(), steps, 2.0 * PI);
    let (ring2, edges2) = build_ring(&b, points[1].x(), points[1].y(), steps, 2.0 * PI);
    Ok(segment_wire(
        &b,
        &ring1,
        &edges1,
        &ring2,
        &edges2,
        points[0].x(),
        points[0].y(),
        points[1].x(),
        points[1].y(),
    ))
}

/// Radius and height of a ring, from its first (θ = 0) vertex.
fn ring_r_z(verts: &[Vertex]) -> (f64, f64) {
    match verts.first() {
        Some(v) => {
            let p = GeometryRegistry::global().vertex_point(&v.0);
            (p.x().hypot(p.y()), p.z())
        }
        None => (0.0, 0.0),
    }
}

/// Volume of the revolved solid by Pappus's centroid theorem.
///
/// For each profile segment the region between the segment and the axis is a
/// trapezoid of area `A = (r1 + r2)·|Δz|/2` and centroid radius
/// `r_c = (r1² + r1·r2 + r2²)/(3·(r1 + r2))`, so
/// `V = 2π·r_c·A = π·|Δz|·(r1² + r1·r2 + r2²)/3`. Horizontal disc segments
/// (|Δz| = 0) contribute zero.
pub fn revolved_volume(rev: &RevolvedSolid) -> f64 {
    let mut vol = 0.0;
    for i in 0..rev.rings.len().saturating_sub(1) {
        let (r1, z1) = ring_r_z(&rev.rings[i]);
        let (r2, z2) = ring_r_z(&rev.rings[i + 1]);
        vol += PI * (z2 - z1).abs() * (r1 * r1 + r1 * r2 + r2 * r2) / 3.0;
    }
    vol
}

/// Swept surface area by Pappus's second theorem: `Σ 2π·r_c·L` per segment,
/// with centroid radius `r_c = (r1 + r2)/2` and segment length `L`.
pub fn revolved_lateral_area(rev: &RevolvedSolid) -> f64 {
    let mut area = 0.0;
    for i in 0..rev.rings.len().saturating_sub(1) {
        let (r1, z1) = ring_r_z(&rev.rings[i]);
        let (r2, z2) = ring_r_z(&rev.rings[i + 1]);
        let l = (r2 - r1).hypot(z2 - z1);
        if l < EPS {
            continue;
        }
        area += PI * (r1 + r2) * l;
    }
    area
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_tool::BRepTool;
    use crate::topo_tools_full;

    fn approx(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9 * b.abs().max(1.0)
    }

    #[test]
    fn rectangle_revolves_to_cylinder() {
        let rev = revolve_rectangle(1.0, 3.0, 16).unwrap();
        // Bottom disc + lateral cylinder + top disc.
        assert_eq!(rev.faces.len(), 3);
        assert!(approx(revolved_volume(&rev), PI * 1.0 * 1.0 * 3.0));
        assert!(approx(revolved_lateral_area(&rev), 2.0 * PI * 1.0 * 3.0 + 2.0 * PI * 1.0 * 1.0));
        assert!(topo_tools_full::structure_is_valid(&rev.solid.0));
    }

    #[test]
    fn closed_profile_is_full_cylinder() {
        let profile = [
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(1.0, 0.0),
            GpPnt2d::new(1.0, 3.0),
            GpPnt2d::new(0.0, 3.0),
        ];
        let rev = revolve_polyline_around_z(&profile, 16).unwrap();
        assert_eq!(rev.faces.len(), 3);
        assert!(approx(revolved_volume(&rev), 3.0 * PI));
    }

    #[test]
    fn open_lateral_profile_is_one_face() {
        let rev = revolve_polyline_around_z(&[GpPnt2d::new(1.0, 0.0), GpPnt2d::new(1.0, 3.0)], 16).unwrap();
        assert_eq!(rev.faces.len(), 1);
        // Pappus treats the segment as closing along the axis: full cylinder.
        assert!(approx(revolved_volume(&rev), 3.0 * PI));
        assert!(approx(revolved_lateral_area(&rev), 2.0 * PI * 1.0 * 3.0));
    }

    #[test]
    fn triangle_revolves_to_cone() {
        let rev = revolve_triangle(1.0, 3.0, 16).unwrap();
        // Lateral conical surface + base disc.
        assert_eq!(rev.faces.len(), 2);
        assert!(approx(revolved_volume(&rev), PI * 1.0 * 1.0 * 3.0 / 3.0));
    }

    #[test]
    fn partial_angle_spans_quarter() {
        let rev = revolve_angle(&[GpPnt2d::new(1.0, 0.0), GpPnt2d::new(1.0, 3.0)], 8, PI / 2.0).unwrap();
        assert_eq!(rev.faces.len(), 1);
        // The first ring runs from θ = 0 to θ = π/2.
        let r0 = &rev.rings[0];
        let p0 = BRepTool::vertex_point(&r0[0]);
        let pend = BRepTool::vertex_point(&r0[r0.len() - 1]);
        assert!(p0.distance(&GpPnt::new(1.0, 0.0, 0.0)) < 1e-9);
        assert!(pend.distance(&GpPnt::new(0.0, 1.0, 0.0)) < 1e-9);
        // The profile wire's first edge is a ring arc ending at π/16; the last
        // ring arc ends at π/2.
        let edges = topo_tools_full::edges_of_wire(&rev.profiles[0]);
        let (a0, b0) = BRepTool::edge_parameters(&edges[0]);
        let (a7, b7) = BRepTool::edge_parameters(&edges[7]);
        assert!(approx(a0, 0.0) && approx(b0, PI / 16.0));
        assert!(approx(a7, 7.0 * PI / 16.0) && approx(b7, PI / 2.0));
    }

    #[test]
    fn open_wire_helper_builds_panel() {
        let w = revolve_open_wire(&[GpPnt2d::new(1.0, 0.0), GpPnt2d::new(1.0, 3.0)], 8).unwrap();
        let edges = topo_tools_full::edges_of_wire(&w);
        // 8 lower ring arcs + seam + 8 upper ring arcs + seam.
        assert_eq!(edges.len(), 18);
    }

    #[test]
    fn validation_rejects_bad_input() {
        assert!(revolve_polyline_around_z(&[GpPnt2d::new(0.0, 0.0), GpPnt2d::new(1.0, 1.0)], 3).is_err());
        assert!(revolve_polyline_around_z(&[GpPnt2d::new(-1.0, 0.0), GpPnt2d::new(1.0, 1.0)], 8).is_err());
        assert!(revolve_polyline_around_z(&[GpPnt2d::new(0.0, 0.0)], 8).is_err());
        assert!(revolve_angle(&[GpPnt2d::new(1.0, 0.0), GpPnt2d::new(1.0, 1.0)], 8, 0.0).is_err());
        assert!(revolve_angle(&[GpPnt2d::new(1.0, 0.0), GpPnt2d::new(1.0, 1.0)], 8, 3.0 * PI).is_err());
        assert!(revolve_rectangle(0.0, 3.0, 8).is_err());
    }
}
