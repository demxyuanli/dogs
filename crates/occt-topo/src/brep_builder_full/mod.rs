//! Unified BRepBuilderAPI construction surface — MakeEdge, MakeWire, MakeFace,
//! MakeShell, MakeSolid, MakePolygon and MakeVertex as ergonomic builders.
//! Source: `BRepBuilderAPI_MakeEdge`, `BRepBuilderAPI_MakeWire`,
//! `BRepBuilderAPI_MakeFace`, `BRepBuilderAPI_MakeShell`,
//! `BRepBuilderAPI_MakeSolid`, `BRepBuilderAPI_MakePolygon`,
//! `BRepBuilderAPI_MakeVertex` (TKTopAlgo).
//!
//! The port is organised as a small set of namespace structs, one per OCCT
//! builder class, plus a handful of free convenience functions. Each builder
//! mirrors the constructor overloads of its OCCT counterpart:
//!
//! * [`BRepBuilderEdge`] — `BRepBuilderAPI_MakeEdge` (points, curves, circles,
//!   ellipses, arcs through three points, B-splines).
//! * [`BRepBuilderWire`] — `BRepBuilderAPI_MakeWire` (validated chains of
//!   edges, open polylines, auto-closed wires).
//! * [`BRepBuilderFace`] — `BRepBuilderAPI_MakeFace` (planar faces from wires,
//!   faces from a surface and a trimming wire, bare surfaces).
//! * [`BRepBuilderShell`] — `BRepBuilderAPI_MakeShell` (faces → shell, closed
//!   box shells, closedness queries).
//! * [`BRepBuilderSolid`] — `BRepBuilderAPI_MakeSolid` (shells → solid,
//!   faces → solid, box primitives).
//! * [`BRepBuilderPolygon`] — `BRepBuilderAPI_MakePolygon` (open/closed
//!   polyline wires).
//!
//! Where an equivalent already exists in [`crate::brep_builder_api`] the new
//! builder delegates to it rather than duplicating geometry logic; this module
//! only *adds* the missing constructor surface (validated wire assembly, the
//! ellipse/B-spline edge constructors, box shells, two-vertex edges and the
//! auto-closing wire helper).
//!
//! # Worked example
//!
//! Building a closed square box solid, then recovering its geometry:
//!
//! ```text
//! use occt_core::gp::{GpAx2, GpPnt};
//! use occt_topo::brep_builder_full::{
//!     BRepBuilderEdge, BRepBuilderFace, BRepBuilderPolygon, BRepBuilderSolid,
//!     BRepBuilderWire,
//! };
//!
//! // Edge chain of a unit square.
//! let e1 = BRepBuilderEdge::from_points(GpPnt::new(0., 0., 0.), GpPnt::new(1., 0., 0.));
//! let e2 = BRepBuilderEdge::from_points(GpPnt::new(1., 0., 0.), GpPnt::new(1., 1., 0.));
//! let e3 = BRepBuilderEdge::from_points(GpPnt::new(1., 1., 0.), GpPnt::new(0., 1., 0.));
//! let e4 = BRepBuilderEdge::from_points(GpPnt::new(0., 1., 0.), GpPnt::new(0., 0., 0.));
//! let square = BRepBuilderWire::from_edges(&[e1, e2, e3, e4]);   // closed wire
//! ```
//!
//! Every builder validates its inputs and returns `Result<_, String>` with a
//! descriptive error, so a malformed construction (coincident points, a
//! disjoint wire chain, a degenerate box) fails loudly instead of producing a
//! broken shape. See [`wire_from_edges_auto`] for a best-effort variant that
//! bridges gaps automatically.
//!
//! # Delegation
//!
//! The low-level shape building is done through [`TopoBuilder`], so every
//! shape returned here carries registered geometry in the
//! [`GeometryRegistry`](crate::tgeometry::GeometryRegistry) side-table and can
//! be queried back through [`BRepTool`]. Existing polygon / arc / primitive
//! helpers in [`crate::brep_builder_api`] are reused verbatim — for example
//! [`BRepBuilderEdge::from_arc_3pts`] and [`BRepBuilderPolygon::from_points`]
//! are thin delegations to `brep_builder_api::make_edge_arc` /
//! `make_polygon`.
mod prelude {

pub(crate) use std::sync::Arc;

pub(crate) use occt_core::gp::{GpAx2, GpDir, GpElips, GpLin, GpPln, GpPnt, GpVec};
pub(crate) use occt_geom::{Curve, GeomBSplineCurve, GeomEllipse, GeomLine, GeomPlane, Surface};

pub(crate) use crate::brep_builder_api;
pub(crate) use crate::brep_tool::BRepTool;
pub(crate) use crate::builder::TopoBuilder;
pub(crate) use crate::shape::{Edge, Face, Shell, Solid, TopoShape, Vertex, Wire};
pub(crate) use crate::topo_tools_full;

}

use prelude::*;

// ============================================================================
// Compatibility notes
// ============================================================================
//
// The builders here map OCCT's `BRepBuilderAPI_*` classes to Rust `Result`
// returning free functions. Differences from OCCT worth knowing:
//
// * **No exception model.** OCCT's `MakeEdge(P1, P2)` raises `StdFail_NotDone`
//   when the points coincide; here the same condition is an `Err(String)`.
//   Every builder validates its preconditions up front and returns a
//   descriptive message instead of panicking (except where noted, e.g. the
//   low-level `BRepPrimBox` dimensions are asserted in `primitives.rs`).
//
// * **Validation is constructive.** `BRepBuilderWire` rejects a chain whose
//   consecutive edges are disjoint, whereas OCCT's `MakeWire` accepts any
//   edges and only reports `WireDone() == false` afterwards. The port
//   deliberately fails fast so a broken boundary is caught at construction
//   time.
//
// * **Vertex sharing.** Edges built from points share their endpoint vertex
//   `TShape`s across the wire (the OCCT shared-TShape model), so a box built
//   from `BRepBuilderWire` has 4 corner vertices, not 8 duplicate ones. This
//   is what makes `shell_euler_characteristic` come out to 2 for a closed box.
//
// * **Planes are unbounded.** A face built with [`BRepBuilderFace::from_wire`]
//   carries a `GeomPlane` whose natural parameter range is infinite; the
//   bounding wire defines the actual trimmed region, exactly as in OCCT.
//
// * **Closedness is a flag, not a guarantee.** [`BRepBuilderWire`] sets the
//   wire's `closed` flag based on endpoint coincidence within [`WIRE_TOL`].
//   The flag is advisory; use `topo_tools_full::wire_is_closed` for a
//   structural re-check on an arbitrary wire.

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_measure::edge_length;
    use crate::brep_surface::{classify_surface, face_is_planar, SurfaceKind};
    use crate::fillet_edge::classify_surface_full;
    use crate::primitives::BRepPrimBox;
    use crate::topo_tools_full::edges_of_wire;
    use occt_core::gp::{GpAx3, GpCylinder};
    use occt_geom::GeomCylinder;
    use std::f64::consts::PI;

    #[test]
    fn edge_from_points_and_curve() {
        // Straight segment of length 1.
        let e = BRepBuilderEdge::from_points(GpPnt::zero(), GpPnt::new(1.0, 0.0, 0.0)).expect("segment");
        assert!((edge_length(&e, 16) - 1.0).abs() < 1e-9, "segment length");
        let (f, l) = BRepTool::edge_parameters(&e);
        assert!((f - 0.0).abs() < 1e-12 && (l - 1.0).abs() < 1e-12);

        // Circle edge: parameters span first..last, curve sits on radius 2.
        let c = BRepBuilderEdge::from_circle(GpAx2::standard(), 2.0, 0.0, PI);
        let (f, l) = BRepTool::edge_parameters(&c);
        assert!((f - 0.0).abs() < 1e-12 && (l - PI).abs() < 1e-12);
        let curve = BRepTool::edge_curve(&c).expect("circle curve");
        assert!(curve.d0(0.0).distance(&GpPnt::new(2.0, 0.0, 0.0)) < 1e-9);
        assert!(curve.d0(PI).distance(&GpPnt::new(-2.0, 0.0, 0.0)) < 1e-9);
    }

    #[test]
    fn edge_from_ellipse() {
        let e = BRepBuilderEdge::from_ellipse(GpAx2::standard(), 4.0, 2.0, 0.0, 2.0 * PI).expect("ellipse edge");
        let curve = BRepTool::edge_curve(&e).expect("ellipse curve");
        // Major vertex at angle 0, minor vertex at π/2. The OCCT ellipse
        // convention is x = a·cos u, y = −b·sin u, so the π/2 point lies at
        // (0, −minor).
        assert!(curve.d0(0.0).distance(&GpPnt::new(4.0, 0.0, 0.0)) < 1e-9, "major vertex");
        assert!(curve.d0(0.5 * PI).distance(&GpPnt::new(0.0, -2.0, 0.0)) < 1e-9, "minor vertex");
        // Invalid radii rejected.
        assert!(BRepBuilderEdge::from_ellipse(GpAx2::standard(), 2.0, 3.0, 0.0, 1.0).is_err());
        assert!(BRepBuilderEdge::from_ellipse(GpAx2::standard(), 0.0, 2.0, 0.0, 1.0).is_err());
    }

    #[test]
    fn edge_arc_3pts() {
        let p1 = GpPnt::new(1.0, 0.0, 0.0);
        let p2 = GpPnt::new(0.0, 1.0, 0.0);
        let p3 = GpPnt::new(-1.0, 0.0, 0.0);
        let e = BRepBuilderEdge::from_arc_3pts(p1, p2, p3).expect("arc");
        let curve = BRepTool::edge_curve(&e).expect("arc curve");
        // The arc curve passes through the middle point at the half parameter.
        assert!(curve.d0(0.5).distance(&p2) < 1e-9, "arc midpoint at p2");
        assert!(curve.d0(0.0).distance(&p1) < 1e-9, "arc start at p1");
        assert!(curve.d0(1.0).distance(&p3) < 1e-9, "arc end at p3");
        // Endpoint vertex children match p1/p3.
        let (a, z) = crate::topo_tools_full::edge_vertices(&e);
        assert!(crate::topo_tools_full::vertex_position(&a.unwrap()).distance(&p1) < 1e-9);
        assert!(crate::topo_tools_full::vertex_position(&z.unwrap()).distance(&p3) < 1e-9);
    }

    #[test]
    fn edge_bspline() {
        let poles = vec![GpPnt::new(0., 0., 0.), GpPnt::new(1., 1., 0.), GpPnt::new(2., 0., 0.)];
        let knots = vec![0., 0., 0., 1., 1., 1.];
        let e = BRepBuilderEdge::from_bspline(&poles, &knots, 2).expect("bspline edge");
        let curve = BRepTool::edge_curve(&e).expect("bspline curve");
        let (f, l) = BRepTool::edge_parameters(&e);
        assert!(curve.d0(f).distance(&poles[0]) < 1e-6, "start at first pole");
        assert!(curve.d0(l).distance(&poles[2]) < 1e-6, "end at last pole");
        // Bad knot count rejected.
        assert!(BRepBuilderEdge::from_bspline(&poles, &[0., 0., 0., 1., 1.], 2).is_err());
    }

    #[test]
    fn wire_from_edges_validates() {
        let e1 = BRepBuilderEdge::from_points(GpPnt::zero(), GpPnt::new(1., 0., 0.)).unwrap();
        let e2 = BRepBuilderEdge::from_points(GpPnt::new(1., 0., 0.), GpPnt::new(1., 1., 0.)).unwrap();
        assert!(BRepBuilderWire::from_edges(&[e1.clone(), e2]).is_ok(), "connected chain ok");

        let e3 = BRepBuilderEdge::from_points(GpPnt::new(5., 0., 0.), GpPnt::new(6., 0., 0.)).unwrap();
        assert!(BRepBuilderWire::from_edges(&[e1, e3]).is_err(), "disjoint chain rejected");
        assert!(BRepBuilderWire::from_edges(&[]).is_err(), "empty rejected");
    }

    #[test]
    fn wire_from_points_and_closed() {
        let pts = [
            GpPnt::new(0., 0., 0.),
            GpPnt::new(1., 0., 0.),
            GpPnt::new(1., 1., 0.),
            GpPnt::new(0., 1., 0.),
        ];
        let open = BRepBuilderWire::from_points(&pts).expect("open wire");
        assert_eq!(edges_of_wire(&open).len(), 3, "open: n-1 edges");
        assert!(!open.closed());

        let closed = BRepBuilderPolygon::from_points(&pts, true).expect("closed wire");
        assert_eq!(edges_of_wire(&closed).len(), 4, "closed: n edges");
        assert!(closed.closed());
        assert!(crate::topo_tools_full::wire_is_closed(&closed), "closed chains end-to-end");
    }

    #[test]
    fn wire_incremental_builder() {
        let a = GpPnt::new(0., 0., 0.);
        let b = GpPnt::new(1., 0., 0.);
        let c = GpPnt::new(1., 1., 0.);
        let mut w = BRepBuilderWire::new();
        w.add_edge(BRepBuilderEdge::from_points(a, b).unwrap()).expect("add e1");
        w.add_edge(BRepBuilderEdge::from_points(b, c).unwrap()).expect("add e2");
        // Disjoint third edge is rejected.
        let far = BRepBuilderEdge::from_points(GpPnt::new(9., 0., 0.), GpPnt::new(8., 0., 0.)).unwrap();
        assert!(w.add_edge(far).is_err(), "disjoint edge rejected");
        let wire = w.build().expect("build wire");
        assert_eq!(edges_of_wire(&wire).len(), 2);
        assert!(!wire.closed());
    }

    #[test]
    fn face_from_wire_planar() {
        let pts = [
            GpPnt::new(0., 0., 0.),
            GpPnt::new(1., 0., 0.),
            GpPnt::new(1., 1., 0.),
            GpPnt::new(0., 1., 0.),
        ];
        let wire = BRepBuilderPolygon::from_points(&pts, true).expect("square wire");
        let plane = GpPln::new(GpAx3::standard());
        let face = BRepBuilderFace::from_wire(&wire, &plane).expect("planar face");
        assert!(face_is_planar(&face));
        let surf = BRepTool::face_surface(&face).expect("surface");
        assert_eq!(classify_surface(surf.as_ref()), SurfaceKind::Plane);
    }

    #[test]
    fn face_from_surface_and_wire() {
        // A cylinder face bounded by two circle wires (bottom + top).
        let b = TopoBuilder::new();
        let ax = GpAx2::standard();
        let r = 2.0;
        let height = 3.0;
        let bottom = b.make_edge_circle(&ax, r, 0.0, 2.0 * PI);
        let mut top_ax = ax;
        top_ax.set_location(GpPnt::new(0.0, 0.0, height));
        let top = b.make_edge_circle(&top_ax, r, 0.0, 2.0 * PI);
        let wire_b = b.make_wire(&[bottom]);
        let wire_t = b.make_wire(&[top]);
        let cyl = GpCylinder::new(GpAx3::standard(), r).expect("cylinder");
        let face = b.make_face(Arc::new(GeomCylinder::new(cyl)), &[wire_b.clone(), wire_t]);
        let surf = BRepTool::face_surface(&face).expect("surface");
        assert_eq!(classify_surface_full(surf.as_ref()), SurfaceKind::Cylinder);

        // Builder path with an explicit wire.
        let f2 = BRepBuilderFace::from_surface_and_wire(
            Arc::new(GeomCylinder::new(GpCylinder::new(GpAx3::standard(), r).unwrap())),
            &wire_b,
        )
        .expect("surface+wire face");
        assert!(BRepTool::face_surface(&f2).is_some());
    }

    #[test]
    fn shell_closed_box() {
        let shell = BRepBuilderShell::closed_box(1.0, 2.0, 3.0).expect("box shell");
        assert!(crate::shell_check::shell_is_closed(&shell), "closed box shell");
        assert_eq!(crate::topo_tools_full::faces_of(&shell.0).len(), 6);
        assert!(BRepBuilderShell::new(shell).is_closed());
        assert!(BRepBuilderShell::closed_box(0.0, 1.0, 1.0).is_err());
    }

    #[test]
    fn solid_from_faces() {
        let bx = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let faces = crate::topo_tools_full::faces_of(&bx.solid.0);
        let solid = BRepBuilderSolid::from_faces(&faces).expect("solid");
        let shell = Shell(solid.0.tshape.read().unwrap().children[0].clone());
        assert!(crate::shell_check::shell_is_closed(&shell), "solid shell closed");
        assert_eq!(crate::topo_tools_full::vertices_of(&solid.0).len(), 8);
        assert!(BRepBuilderSolid::from_faces(&[]).is_err());
    }

    #[test]
    fn polygon_open_closed() {
        let pts = [
            GpPnt::new(0., 0., 0.),
            GpPnt::new(1., 0., 0.),
            GpPnt::new(1., 1., 0.),
            GpPnt::new(0., 1., 0.),
            GpPnt::new(0., 0., 1.),
        ];
        let open = BRepBuilderPolygon::from_points(&pts, false).expect("open polygon");
        assert_eq!(edges_of_wire(&open).len(), pts.len() - 1, "open: n-1 edges");
        assert!(!open.closed());
        let closed = BRepBuilderPolygon::from_points(&pts, true).expect("closed polygon");
        assert_eq!(edges_of_wire(&closed).len(), pts.len(), "closed: n edges");
        assert!(closed.closed());
    }

    #[test]
    fn edge_between_vertices() {
        let b = TopoBuilder::new();
        let v1 = b.make_vertex(GpPnt::new(0., 0., 0.), 0.0);
        let v2 = b.make_vertex(GpPnt::new(2., 0., 0.), 0.0);
        let e = edge_from_two_vertices(&v1, &v2).expect("segment");
        assert!((edge_length(&e, 16) - 2.0).abs() < 1e-9, "length 2");
        let (a, z) = BRepTool::edge_vertices(&e).expect("endpoints");
        assert!(a.distance(&GpPnt::zero()) < 1e-9);
        assert!(z.distance(&GpPnt::new(2., 0., 0.)) < 1e-9);
        // The two input vertices are reused as the edge's boundary children.
        let (va, vz) = crate::topo_tools_full::edge_vertices(&e);
        assert!(crate::topo_tools_full::is_same(&va.unwrap().0, &v1.0));
        assert!(crate::topo_tools_full::is_same(&vz.unwrap().0, &v2.0));
    }

    #[test]
    fn prism_solid_from_base() {
        let tri = [GpPnt::new(0., 0., 0.), GpPnt::new(1., 0., 0.), GpPnt::new(0., 1., 0.)];
        let solid = solid_from_prism(&tri, &GpVec::new(0., 0., 1.), 3.0).expect("prism");
        let shell = Shell(solid.0.tshape.read().unwrap().children[0].clone());
        assert!(crate::shell_check::shell_is_closed(&shell), "triangle prism closed");
        // Boundary counts: 6 vertices, 9 edges, 5 faces.
        assert_eq!(crate::topo_tools_full::vertices_of(&solid.0).len(), 6);
        assert_eq!(crate::topo_tools_full::edges_of(&solid.0).len(), 9);
        assert_eq!(crate::topo_tools_full::faces_of(&solid.0).len(), 5);
        // Analytic prism volume = base area × sweep height, recovered from the
        // solid's own geometry (the mesh-based volume is unreliable on a
        // triangle base: the face grid overhangs the boundary).
        let faces = crate::topo_tools_full::faces_of(&solid.0);
        let zc = |f: &Face| crate::brep_surface::face_centroid(f, 8, 8).map(|p| p.z()).unwrap_or(f64::NAN);
        let base = faces.iter().min_by(|a, b| zc(a).partial_cmp(&zc(b)).unwrap()).unwrap();
        let top = faces.iter().max_by(|a, b| zc(a).partial_cmp(&zc(b)).unwrap()).unwrap();
        let base_plane = crate::brep_surface::face_plane(base).expect("base plane");
        let top_plane = crate::brep_surface::face_plane(top).expect("top plane");
        let height = top_plane.location().z() - base_plane.location().z();
        assert!((height - 3.0).abs() < 1e-9, "sweep height {height}");
        let wire = crate::topo_tools_full::wires_of_face(base).into_iter().next().expect("base wire");
        let ring: Vec<GpPnt> = crate::topo_tools_full::edges_of_wire(&wire)
            .iter()
            .filter_map(|e| crate::topo_tools_full::edge_vertices(e).0)
            .map(|v| crate::topo_tools_full::vertex_position(&v))
            .collect();
        let area = crate::brep_builder_api::polygon_area(&ring);
        assert!((area - 0.5).abs() < 1e-9, "base area {area}");
        let vol = area * height;
        assert!(vol > 0.0, "volume positive");
        assert!((vol - 1.5).abs() < 1e-9, "volume {vol} ~ 1.5");
        // Invalid inputs rejected.
        assert!(solid_from_prism(&tri, &GpVec::new(0., 0., 1.), 0.0).is_err());
        assert!(solid_from_prism(&[GpPnt::new(0., 0., 0.), GpPnt::new(1., 0., 0.)], &GpVec::new(0., 0., 1.), 1.0).is_err());
        assert!(solid_from_prism(&tri, &GpVec::new(0., 0., 0.), 1.0).is_err());
    }

    #[test]
    fn collinear_arc_errors() {
        let e = BRepBuilderEdge::from_arc_3pts(
            GpPnt::new(0., 0., 0.),
            GpPnt::new(1., 0., 0.),
            GpPnt::new(2., 0., 0.),
        );
        assert!(e.is_err(), "collinear points cannot define an arc");
    }
}

mod p01;
pub use p01::*;
