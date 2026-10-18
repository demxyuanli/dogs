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
//!
//! Ring and seam structure follows `BRepPrim_OneAxis`, which is the base of
//! both `BRepPrim_Revolution` and `BRepSweep_Rotation`:
//! * the lateral wire order and occurrence orientations are
//!   `BRepPrim_OneAxis::LateralWire` (`BRepPrim_OneAxis.cxx:666-679` plus
//!   `BRepPrim_Builder::AddWireEdge`, `BRepPrim_Builder.cxx:184-192`): bottom
//!   ring FORWARD, seam at u = 0 REVERSED, top ring REVERSED, seam at
//!   u = angle FORWARD;
//! * a full revolution has `!HasSides()` (`cxx:331-334`), so `StartEdge()`
//!   aliases `EndEdge()` (`cxx:976-979`, `:1028-1031`) and **one** seam TShape
//!   is used twice;
//! * the pcurves are `BRepPrim_OneAxis::LateralFace` (`cxx:388-439`) — ring
//!   arcs as `gp_Lin2d((0, V), +X)` over their angular range, the seam as the
//!   closed-form pair `gp_Lin2d((angle, 0), +Y)` / `gp_Lin2d((0, 0), +Y)` —
//!   and `TopFace`/`BottomFace` (`cxx:465-468`, `:504-509`) — cap rings as
//!   `gp_Circ2d(gp_Ax2d((0, 0), +X), R)`;
//! * ring vertices carry `Precision::Confusion()` like
//!   `BRepPrim_Builder::MakeVertex` (`BRepPrim_Builder.cxx:129-132`).
//!
//! Measured 2026-10-16 (package P039/T-169): handing this analytic solid to
//! `crate::bop_builder::boolean` still fails — `box ∪ analytic neck` drops the
//! tool and `faceted cylinder − analytic groove ring` adds material — so the
//! remaining gap is not the OneAxis wire/seam/pcurve conventions alone.
//! `crate::brepfeat` therefore keeps its lathe-mesh tool (UNPORTED workaround).

use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::gp::{
    GpAx2, GpAx22d, GpAx3, GpCirc2d, GpCone, GpCylinder, GpDir, GpDir2d, GpLin, GpPln, GpPnt,
    GpPnt2d, GpVec,
};
use occt_geom::{GeomCone, GeomCylinder, GeomLine, GeomPlane, Surface};
use occt_geom2d::{Geom2dCircle, Geom2dLine};

use crate::abs::Orientation;
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
            // `BRepPrim_Builder::MakeVertex` (`BRepPrim_Builder.cxx:129-132`) is
            // `BRep_Builder::MakeVertex(V, P, Precision::Confusion())`; the 1e-7
            // tolerance is what `BOPAlgo_WireSplitter::Tolerance2D`
            // (`_1.cxx:859-881`) turns into the UV closing tolerance, so a zero
            // tolerance keeps the seam loops from closing (aD2 < 0).
            verts.push(b.make_vertex(
                GpPnt::new(r * theta.cos(), r * theta.sin(), z),
                occt_core::precision::CONFUSION,
            ));
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

/// Angular start/end of each of the `steps` arcs of a ring spanning `[0, total]`.
fn ring_arcs(steps: usize, total: f64) -> Vec<(f64, f64)> {
    (0..steps)
        .map(|j| (total * j as f64 / steps as f64, total * (j + 1) as f64 / steps as f64))
        .collect()
}

/// The face swept by one profile segment, plus the pieces whose pcurves have to
/// be registered on it.
struct LateralPanel {
    wire: Wire,
    /// Boundary arcs lying on the ring of the segment's first profile point.
    ring1: Vec<Edge>,
    /// Boundary arcs lying on the ring of the segment's second profile point.
    ring2: Vec<Edge>,
    /// Seam edges as *stored* (first ring point -> second ring point). One entry
    /// for a full revolution — `BRepPrim_OneAxis` then aliases `StartEdge()` and
    /// `EndEdge()` (`BRepPrim_OneAxis.cxx:976-979`, `:1028-1031`), so the same
    /// edge TShape is used twice — two entries (u = 0, u = total) otherwise.
    seams: Vec<Edge>,
}

impl LateralPanel {
    fn full(&self) -> bool {
        self.seams.len() == 1
    }
}

/// Boundary wire of a cylindrical/conical face between two rings.
///
/// Wire order and occurrence orientations follow
/// `BRepPrim_OneAxis::LateralWire` (`BRepPrim_OneAxis.cxx:666-679`):
/// `AddWireEdge(TopEdge(), false)`, `AddWireEdge(EndEdge(), true)`,
/// `AddWireEdge(BottomEdge(), true)`, `AddWireEdge(StartEdge(), false)`.
/// `BRepPrim_Builder::AddWireEdge` (`BRepPrim_Builder.cxx:184-192`) applies
/// `EE.Reverse()` when `direct == false`, and for a full revolution
/// `StartEdge()` aliases `EndEdge()`. Rotating the loop so the bottom ring
/// comes first (a wire is a cycle) gives
/// **bottom ring FORWARD, seam at u = 0 REVERSED, top ring REVERSED,
/// seam at u = total FORWARD**.
fn lateral_wire(
    b: &TopoBuilder,
    ring1: &[Vertex],
    edges1: &[Edge],
    ring2: &[Vertex],
    edges2: &[Edge],
    total: f64,
) -> LateralPanel {
    let last1 = ring1.len() - 1;
    let last2 = ring2.len() - 1;
    let full = (total - 2.0 * PI).abs() < 1e-9;
    let seam_start = edge_through(b, &ring1[0], &ring2[0]);
    let seam_end = if full {
        seam_start.clone()
    } else {
        edge_through(b, &ring1[last1], &ring2[last2])
    };
    let mut seam_start_back = seam_start.clone();
    seam_start_back.0.reverse();
    let mut edges: Vec<Edge> = Vec::with_capacity(edges1.len() + edges2.len() + 2);
    edges.extend(edges1.iter().cloned());
    edges.push(seam_start_back);
    edges.extend(edges2.iter().rev().map(|e| Edge(e.0.oriented(Orientation::Reversed))));
    edges.push(seam_end.clone());
    let seams = if full { vec![seam_start] } else { vec![seam_start, seam_end] };
    LateralPanel { wire: b.make_wire(&edges), ring1: edges1.to_vec(), ring2: edges2.to_vec(), seams }
}

/// The analytic surface swept by a profile segment `(r1,z1)→(r2,z2)`, together
/// with the V parameter of the meridian at each end.
///
/// `BRepPrim_OneAxis` uses the *meridian parameter* V: the pcurves of the
/// lateral face are `gp_Lin2d((0, V), +X)` for the ring edges
/// (`BRepPrim_OneAxis::LateralFace`, `BRepPrim_OneAxis.cxx:388-439`). For a
/// cylinder of revolution V is the axial abscissa from the surface origin; for
/// a cone it is the distance from the apex along the axis
/// (`ElSLib::ConeVIso`, `ElSLib.cxx:1793-1810`: radius `R + V*sin(A)`, center
/// `V*cos(A)` along the axis).
struct SegmentSurface {
    surface: Arc<dyn Surface>,
    v1: f64,
    v2: f64,
}

fn segment_surface(r1: f64, z1: f64, r2: f64, z2: f64) -> Result<SegmentSurface, String> {
    let dr = r2 - r1;
    let dz = z2 - z1;
    if dz.abs() < EPS {
        // Horizontal segment → planar disc/annulus at height z1.
        let pln = GpPln::new(GpAx3::new(GpPnt::new(0.0, 0.0, z1), z_dir(), &x_dir())?);
        Ok(SegmentSurface { surface: Arc::new(GeomPlane::new(pln)), v1: 0.0, v2: 0.0 })
    } else if dr.abs() < EPS {
        // Vertical segment → cylinder of radius r1 between z1 and z2.
        let ax3 = GpAx3::new(GpPnt::new(0.0, 0.0, z1), z_dir(), &x_dir())?;
        Ok(SegmentSurface {
            surface: Arc::new(GeomCylinder::new(GpCylinder::new(ax3, r1)?)),
            v1: 0.0,
            v2: dz.abs(),
        })
    } else {
        // Diagonal segment → cone frustum (full cone when one radius is 0).
        // The radius is linear in z; the apex is where it vanishes.
        let semi_angle = (dr.abs() / dz.abs()).atan();
        let z_apex = z1 - r1 * dz / dr;
        let upward = z_apex <= z1.min(z2) + 1e-9;
        let dir = if upward { z_dir() } else { GpDir::new(0.0, 0.0, -1.0).expect("sweep_revolve: Z axis") };
        let ax3 = GpAx3::new(GpPnt::new(0.0, 0.0, z_apex), dir, &x_dir())?;
        let sign = if upward { 1.0 } else { -1.0 };
        Ok(SegmentSurface {
            surface: Arc::new(GeomCone::new(GpCone::new(ax3, 0.0, semi_angle)?)),
            v1: sign * (z1 - z_apex),
            v2: sign * (z2 - z_apex),
        })
    }
}

/// Pcurves of a lateral face (`BRepPrim_OneAxis::LateralFace`,
/// `BRepPrim_OneAxis.cxx:388-439`):
///   * `ETOP` / `EBOTTOM` — the ring edges — carry `gp_Lin2d((0, VMax|VMin),
///     +X)`; the UV abscissa is the angle around the axis, so each arc's
///     pcurve is that line restricted to the arc's angular range;
///   * `ESTART` (u = 0) / `EEND` (u = angle) — the seam — carry
///     `gp_Lin2d((0, -offset), +Y)` and `gp_Lin2d((angle, -offset), +Y)`.
///     For a full revolution `SetPCurve(E, F, c1, c2)` (`cxx:434-439`) stores
///     *both* on the same edge, which is what makes
///     `BRep_Tool::IsClosed(seam, lateral)` true (`BRep_Tool.cxx:820-840`).
fn attach_lateral_pcurves(face: &Face, panel: &LateralPanel, srf: &SegmentSurface, steps: usize, total: f64) {
    let reg = GeometryRegistry::global();
    let fk = GeometryRegistry::shape_key(&face.0);
    let dir_u = GpDir2d::new(1.0, 0.0).expect("u direction");
    let dir_v = GpDir2d::new(0.0, 1.0).expect("v direction");
    let arcs = ring_arcs(steps, total);
    for (ring, v) in [(&panel.ring1, srf.v1), (&panel.ring2, srf.v2)] {
        for (e, arc) in ring.iter().zip(arcs.iter()) {
            reg.set_edge_pcurve(
                &e.0,
                fk,
                Arc::new(Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, v), dir_u)),
            );
            reg.set_pcurve_range(&e.0, fk, arc.0, arc.1);
        }
    }
    // The seam's 3-D curve is a plain line parameterized by arc length, so its
    // pcurve range is the meridian span `[0, v2 - v1]` (identical to the arc
    // length for a cylinder).
    let span = srf.v2 - srf.v1;
    if panel.full() {
        reg.set_edge_pcurves(
            &panel.seams[0].0,
            fk,
            vec![
                Arc::new(Geom2dLine::from_pnt_dir(GpPnt2d::new(total, 0.0), dir_v)),
                Arc::new(Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, 0.0), dir_v)),
            ],
        );
        reg.set_pcurve_range(&panel.seams[0].0, fk, 0.0, span);
    } else {
        reg.set_edge_pcurve(
            &panel.seams[0].0,
            fk,
            Arc::new(Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, 0.0), dir_v)),
        );
        reg.set_pcurve_range(&panel.seams[0].0, fk, 0.0, span);
        reg.set_edge_pcurve(
            &panel.seams[1].0,
            fk,
            Arc::new(Geom2dLine::from_pnt_dir(GpPnt2d::new(total, 0.0), dir_v)),
        );
        reg.set_pcurve_range(&panel.seams[1].0, fk, 0.0, span);
    }
}

/// Pcurves of a planar cap (`BRepPrim_OneAxis::TopFace` / `BottomFace`,
/// `BRepPrim_OneAxis.cxx:465-468`, `:504-509`): the ring edge carries a
/// `gp_Circ2d(gp_Ax2d(gp_Pnt2d(0, 0), +X), R)` — the cap plane's own UV frame,
/// whose abscissa is the angle around the axis.
fn attach_cap_pcurves(face: &Face, rings: &[(&[Edge], f64)], steps: usize, total: f64) {
    let reg = GeometryRegistry::global();
    let fk = GeometryRegistry::shape_key(&face.0);
    let dir_u = GpDir2d::new(1.0, 0.0).expect("u direction");
    let dir_v = GpDir2d::new(0.0, 1.0).expect("v direction");
    let ax22d = GpAx22d::new(GpPnt2d::new(0.0, 0.0), dir_u, dir_v).expect("cap frame");
    let arcs = ring_arcs(steps, total);
    for (ring, radius) in rings {
        let circle = Arc::new(Geom2dCircle::new(GpCirc2d::new(ax22d, *radius)));
        for (e, arc) in ring.iter().zip(arcs.iter()) {
            reg.set_edge_pcurve(&e.0, fk, circle.clone());
            reg.set_pcurve_range(&e.0, fk, arc.0, arc.1);
        }
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
    total: f64,
) -> Wire {
    if (z2 - z1).abs() < EPS {
        if r1 >= r2 {
            b.make_wire(edges1)
        } else {
            b.make_wire(edges2)
        }
    } else {
        lateral_wire(b, ring1, edges1, ring2, edges2, total).wire
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
    steps: usize,
    total: f64,
) -> Result<Option<(Face, Wire)>, String> {
    let dr = r2 - r1;
    let dz = z2 - z1;
    if (r1.abs() < EPS && r2.abs() < EPS) || (dr.abs() < EPS && dz.abs() < EPS) {
        return Ok(None);
    }
    let srf = segment_surface(r1, z1, r2, z2)?;
    let profile_wire;
    let mut wires: Vec<Wire> = Vec::new();
    if dz.abs() < EPS {
        // Horizontal: outer circle wire; when both radii are positive the inner
        // ring is added as a hole wire (planar annulus).
        let (outer_edges, outer_r, inner_edges, inner_r) = if r1 >= r2 {
            (edges1, r1, edges2, r2)
        } else {
            (edges2, r2, edges1, r1)
        };
        profile_wire = b.make_wire(outer_edges);
        wires.push(profile_wire.clone());
        let mut caps: Vec<(&[Edge], f64)> = vec![(outer_edges, outer_r)];
        if inner_r > EPS {
            wires.push(b.make_wire(inner_edges));
            caps.push((inner_edges, inner_r));
        }
        let face = b.make_face(srf.surface.clone(), &wires);
        attach_cap_pcurves(&face, &caps, steps, total);
        return Ok(Some((face, profile_wire)));
    }
    let panel = lateral_wire(b, ring1, edges1, ring2, edges2, total);
    profile_wire = panel.wire.clone();
    wires.push(profile_wire.clone());
    let face = b.make_face(srf.surface.clone(), &wires);
    attach_lateral_pcurves(&face, &panel, &srf, steps, total);
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
            steps,
            angle,
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
        2.0 * PI,
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
    use crate::topo_tools;

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
        assert!(topo_tools::structure_is_valid(&rev.solid.0));
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
        let edges = topo_tools::edges_of_wire(&rev.profiles[0]);
        let (a0, b0) = BRepTool::edge_parameters(&edges[0]);
        let (a7, b7) = BRepTool::edge_parameters(&edges[7]);
        assert!(approx(a0, 0.0) && approx(b0, PI / 16.0));
        assert!(approx(a7, 7.0 * PI / 16.0) && approx(b7, PI / 2.0));
    }

    #[test]
    fn open_wire_helper_builds_panel() {
        let w = revolve_open_wire(&[GpPnt2d::new(1.0, 0.0), GpPnt2d::new(1.0, 3.0)], 8).unwrap();
        let edges = topo_tools::edges_of_wire(&w);
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
