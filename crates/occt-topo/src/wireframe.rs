//! Wireframe & face tessellation.
//!
//! Port of the *meshing* half of `BRepMesh_IncrementalMesh`: an edge is turned
//! into a deflection-bounded polyline and a face into a UV-grid triangle soup.
//! Geometry is read from the side-table registry (`tgeometry::GeometryRegistry`)
//! — the same data `BRep_Tool` reads.
//!
//! Edge polylines go through the faithful tessellator
//! (`meshing::edge_discret::CurveTessellator` = `GCPnts_TangentialDeflection`
//! initialize/PerformCurve), matching OCCT's mesh pipeline where
//! `BRepMesh_CurveTessellator` uses `GCPnts_TangentialDeflection`.
//! `GCPnts_UniformDeflection` is **UNPORTED** (see `occt-core/src/gcpnts.rs`).

use occt_core::gp::GpPnt;
use occt_core::gcpnts::{CurveSample, UniformPoints};
use occt_core::poly::triangulation::Triangle;
use occt_geom::Curve;

use crate::meshing::edge_discret::{tessellator_min_points, CurveTessellator};
use crate::shape::{Edge, Face};
use crate::tgeometry::GeometryRegistry;

fn reg() -> &'static GeometryRegistry {
    GeometryRegistry::global()
}

/// Adapter making a `&dyn Curve` sampleable by the uniform fallback sampler.
struct CurveAdapter<'a> {
    curve: &'a dyn Curve,
}
impl CurveSample for CurveAdapter<'_> {
    fn point(&self, u: f64) -> GpPnt {
        self.curve.d0(u)
    }
}

/// Discretize an edge into a polyline of 3D points.
///
/// The curve is tessellated over the edge's parameter range with the faithful
/// `GCPnts_TangentialDeflection` engine (`GCPnts_TangentialDeflection::initialize`
/// `cxx:415-453` + `PerformCurve` `cxx:522-916`), i.e. what OCCT's
/// `BRepMesh_CurveTessellator` calls; `min_points` comes from
/// `BRepMesh_CurveTessellator.cxx:100-114` (`tessellator_min_points`).
/// Both endpoints are always included. If the edge has no registered curve (or
/// an unbounded range) an empty polyline is returned; a fallback of 64 uniform
/// samples covers the degenerate case.
pub fn edge_to_polyline(e: &Edge, deflection: f64) -> Vec<GpPnt> {
    let Some(curve) = reg().edge_curve(e) else { return Vec::new() };
    let (f0, f1) = reg().edge_parameters(e);
    let (a, b) = if f0.is_finite() && f1.is_finite() && f1 >= f0 { (f0, f1) } else { return Vec::new() };
    if b - a < 1e-12 {
        return vec![curve.d0(a)];
    }
    let tol = deflection.max(1e-9);
    let min_points = tessellator_min_points(curve.as_ref());
    let tess = CurveTessellator::from_range(curve.clone(), a, b, tol, min_points);
    if tess.points_nb() >= 2 {
        tess.points().to_vec()
    } else {
        UniformPoints::from_curve(&CurveAdapter { curve: curve.as_ref() }, a, b, 64).points
    }
}

/// Maximum deviation of a sampled polyline from the true edge curve.
///
/// Each polyline segment is mapped to a parameter span by chord-length
/// parametrization and the curve is sampled at 64 interior points per span;
/// the furthest distance to the chord is returned.
pub fn edge_chord_error(e: &Edge, polyline: &[GpPnt]) -> f64 {
    let n = polyline.len();
    if n < 2 { return 0.0; }
    let Some(curve) = reg().edge_curve(e) else { return 0.0 };
    let (f0, f1) = reg().edge_parameters(e);
    if !(f0.is_finite() && f1.is_finite()) || f1 <= f0 { return 0.0; }

    // Cumulative chord lengths → curve parameter per polyline point.
    let mut s = Vec::with_capacity(n);
    s.push(0.0);
    for i in 1..n {
        s.push(s[i - 1] + polyline[i].distance(&polyline[i - 1]));
    }
    let total = *s.last().unwrap();
    if total <= 1e-12 { return 0.0; }

    let mut max_dev = 0.0f64;
    for i in 0..n - 1 {
        let ua = f0 + (f1 - f0) * s[i] / total;
        let ub = f0 + (f1 - f0) * s[i + 1] / total;
        let pa = &polyline[i];
        let pb = &polyline[i + 1];
        for k in 1..64 {
            let u = ua + (ub - ua) * k as f64 / 64.0;
            let d = point_segment_dist(&curve.d0(u), pa, pb);
            if d > max_dev { max_dev = d; }
        }
    }
    max_dev
}

/// Tessellate a face into a triangle soup (vertices + triangles).
///
/// Delegates to the ported `BRepMesh` pipeline (`IncrementalMesh`, the same
/// control flow as `BRepMesh_IncrementalMesh`: `BRepMesh_FaceDiscret` builds the
/// boundary from the face pcurves, then `BRepMesh_Delaun` triangulates the
/// interior — see `meshing::incremental_mesh`). OCCT has no alternative
/// tessellator, so a face whose Delaunay fails keeps `IMeshData_Failure`
/// (`BRepMesh_BaseMeshAlgo.cxx:40-62`) and yields an empty soup. Winding follows
/// the face occurrence orientation (`RWMesh_FaceIterator::TriangleOriented`,
/// applied by `build_shape_mesh`).
pub fn face_to_triangles(f: &Face, deflection: f64) -> (Vec<GpPnt>, Vec<Triangle>) {
    let inc = crate::meshing::incremental_mesh::IncrementalMesh::new();
    inc.discretize_face(f, deflection)
}

/// Sum of triangle areas of a face's tessellation.
pub fn wireframe_face_area(f: &Face) -> f64 {
    let (pts, tris) = face_to_triangles(f, 0.01);
    tris.iter().map(|t| triangle_area(&pts, t)).sum()
}

/// UV domain of a face -- `BRepTools::UVBounds(F, UMin, UMax, VMin, VMax)`.
///
/// Delegates to the ported `BRepTools::AddUVBounds` chain
/// (`crate::brep_uv_bounds`): the 2-D box of the face's pcurve boundaries, with
/// the stored surface's natural range as OCCT's own fallback for a face whose
/// edge box is void (`BRepTools.cxx:126-157`). OCCT builds the box from the
/// `BndLib_Add2dCurve` range of each pcurve and clamps it to the surface bounds
/// unless a B-spline periodicity probe succeeds; there is no UV sampling and no
/// `[0,1]^2` window, so neither is reproduced here.
pub(crate) use crate::brep_surface::face_uv_bounds;

/// Triangle area by cross-product magnitude.
fn triangle_area(pts: &[GpPnt], t: &Triangle) -> f64 {
    let ab = pts[t.n1].coord.subtracted(&pts[t.n0].coord);
    let ac = pts[t.n2].coord.subtracted(&pts[t.n0].coord);
    0.5 * ab.crossed(&ac).modulus()
}

/// Distance from `p` to the segment `a..b`.
fn point_segment_dist(p: &GpPnt, a: &GpPnt, b: &GpPnt) -> f64 {
    let ab = b.coord.subtracted(&a.coord);
    let len2 = ab.square_modulus();
    if len2 <= f64::EPSILON {
        return p.distance(a);
    }
    let ap = p.coord.subtracted(&a.coord);
    let t = (ap.dot(&ab) / len2).clamp(0.0, 1.0);
    let proj = a.coord.added(&ab.multiplied(t));
    p.coord.subtracted(&proj).modulus()
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use std::sync::Arc;
    use occt_core::gp::{GpDir, GpLin, GpPln};
    use occt_geom::{GeomLine, GeomPlane};
    use crate::builder::TopoBuilder;
    use crate::shape::Wire;
    use crate::tgeometry::{EdgeGeom, FaceGeom};

    // ---- shared fixtures (also used by sibling modules' tests) ----

    /// Edge along `(x0,y0)` → `(x1,y1)` in the z=0 plane, parameterized [0,1].
    pub(crate) fn segment_edge(x0: f64, y0: f64, x1: f64, y1: f64) -> Edge {
        let e = Edge::new();
        let dir = GpDir::new(x1 - x0, y1 - y0, 0.0).unwrap();
        let lin = GpLin::from_pnt_dir(GpPnt::new(x0, y0, 0.0), dir);
        reg().set_edge(&e.0, EdgeGeom::new(Arc::new(GeomLine::new(lin)), 0.0, 1.0));
        e
    }

    /// Unit-square face in the z=0 plane (surface registered, 4-edge wire).
    pub(crate) fn square_face() -> Face {
        let b = TopoBuilder::new();
        let mut face = Face::new();
        reg().set_face(&face.0, FaceGeom::new(Arc::new(GeomPlane::new(GpPln::default()))));
        let mut wire = Wire::new();
        for (x0, y0, x1, y1) in [
            (0.0, 0.0, 1.0, 0.0),
            (1.0, 0.0, 1.0, 1.0),
            (1.0, 1.0, 0.0, 1.0),
            (0.0, 1.0, 0.0, 0.0),
        ] {
            let e = segment_edge(x0, y0, x1, y1);
            b.add_edge(&mut wire, &e);
        }
        b.add_wire(&mut face, &wire);
        face
    }

    /// A face with an outer 2×2 square and a 24-gon hole (a box-top ring cut by
    /// a boss cylinder) triangulates to the ring area, not the full square.
    #[test]
    fn face_with_hole_triangulates_ring_area() {
        use occt_core::gp::{GpAx3, GpPnt};
        let b = TopoBuilder::new();
        let mut face = Face::new();
        let pln = GpPln::new(GpAx3::new(GpPnt::new(0.0, 0.0, 1.0), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap());
        reg().set_face(&face.0, FaceGeom::new(Arc::new(GeomPlane::new(pln.clone()))));
        // Outer 2×2 square. Consecutive edges must share a vertex: a wire with
        // no common vertex between neighbors is marked `IMeshData_OpenWire` by
        // `BRepMesh_ModelHealer::fixFaceBoundaries` (`BRepMesh_ModelHealer.cxx:
        // 346-356` + `getCommonVertex`), and an inner open wire is skipped in
        // `BRepMesh_NodeInsertionMeshAlgo::initDataStructure` (`hxx:70-74`), so
        // its hole would not be subtracted.
        let mut outer = Wire::new();
        let outer_pts = [
            GpPnt::new(0.0, 0.0, 1.0),
            GpPnt::new(2.0, 0.0, 1.0),
            GpPnt::new(2.0, 2.0, 1.0),
            GpPnt::new(0.0, 2.0, 1.0),
        ];
        let outer_verts: Vec<_> = outer_pts.iter().map(|p| b.make_vertex(*p, 0.0)).collect();
        for i in 0..outer_pts.len() {
            let j = (i + 1) % outer_pts.len();
            let e = b.make_edge_segment_with_vertices(
                &outer_pts[i],
                &outer_pts[j],
                &outer_verts[i],
                &outer_verts[j],
            );
            b.add_edge(&mut outer, &e);
        }
        b.add_wire(&mut face, &outer);
        // Hole: 24-gon at (1,1), r=0.25, also with shared vertices.
        //
        // An inner wire is stored wound the other way round from the outer one
        // (TopoDS face convention, what the readers and `ShapeFix_Face`
        // produce), so the 24-gon is traversed clockwise here. The mesher
        // itself does not fix wire winding: `BRepMesh_ShapeVisitor::Visit`
        // (`BRepMesh_ShapeVisitor.cxx:58-87`) only does
        // `ShapeAnalysis::OuterWire` + `TopExp_Explorer` + `addWire`
        // (`cxx:93-134`), and `BRepMesh_Classifier::RegisterWire`
        // (`BRepMesh_Classifier.cxx:65-127`) orients each wire by the turning
        // angle of its UV polygon, `Perform` (`cxx:35-59`) returns OUT as soon
        // as one wire's oriented test says out. A same-wound inner wire would
        // therefore be read as a second outer bound (the hole area would be the
        // only IN region), not as a hole. `ShapeFix_Face::FixOrientation`
        // (`ShapeFix_Face.cxx:1274-1602`) belongs to the reader path, not here.
        let r = 0.25;
        let n = 24usize;
        let mut hole = Wire::new();
        let hole_pts: Vec<GpPnt> = (0..n)
            .map(|i| {
                let t = -2.0 * std::f64::consts::PI * i as f64 / n as f64;
                GpPnt::new(1.0 + r * t.cos(), 1.0 + r * t.sin(), 1.0)
            })
            .collect();
        let hole_verts: Vec<_> = hole_pts.iter().map(|p| b.make_vertex(*p, 0.0)).collect();
        for i in 0..n {
            let j = (i + 1) % n;
            let e = b.make_edge_segment_with_vertices(
                &hole_pts[i],
                &hole_pts[j],
                &hole_verts[i],
                &hole_verts[j],
            );
            b.add_edge(&mut hole, &e);
        }
        b.add_wire(&mut face, &hole);
        let (pts, tris) = face_to_triangles(&face, 0.05);
        let mut area = 0.0;
        for t in &tris {
            let a = pts[t.n0].coord;
            let x = pts[t.n1].coord.subtracted(&a);
            let y = pts[t.n2].coord.subtracted(&a);
            area += 0.5 * x.crossed(&y).modulus();
        }
        let expect = 4.0 - std::f64::consts::PI * r * r;
        assert!(
            (area - expect).abs() < 0.05,
            "ring area {area} vs expected {expect} (hole not subtracted?)"
        );
    }

    #[test]
    fn edge_polyline_endpoints_match() {
        let e = segment_edge(0.0, 0.0, 1.0, 0.0);
        let poly = edge_to_polyline(&e, 0.01);
        assert!(poly.len() >= 2, "poly len {}", poly.len());
        assert!(poly.first().unwrap().distance(&GpPnt::new(0.0, 0.0, 0.0)) < 1e-9);
        assert!(poly.last().unwrap().distance(&GpPnt::new(1.0, 0.0, 0.0)) < 1e-9);
    }

    #[test]
    fn edge_chord_error_straight_line_is_zero() {
        let e = segment_edge(0.0, 0.0, 1.0, 0.0);
        let poly = edge_to_polyline(&e, 0.01);
        assert!(edge_chord_error(&e, &poly) < 1e-9, "err {}", edge_chord_error(&e, &poly));
    }

    #[test]
    fn face_triangulates_to_unit_area() {
        let f = square_face();
        let (pts, tris) = face_to_triangles(&f, 0.25);
        // Planar faces triangulate exactly from their boundary polygon: a
        // square is 2 triangles (not the dense UV grid).
        assert_eq!(tris.len(), 2, "tris {}", tris.len());
        let area: f64 = tris.iter().map(|t| triangle_area(&pts, t)).sum();
        assert!((area - 1.0).abs() < 1e-6, "area {area}");
    }

    #[test]
    fn face_area_of_unit_square() {
        let f = square_face();
        assert!((wireframe_face_area(&f) - 1.0).abs() < 1e-6, "area {}", wireframe_face_area(&f));
    }
}
