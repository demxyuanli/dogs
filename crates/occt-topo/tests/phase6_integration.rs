//! Phase 6 end-to-end: exact surface intersection (IntPatch), offset
//! operations, BinXCAF binary container, glTF export, polygon boolean and
//! curve/surface extrema.
//!
//! These tests exercise the Phase 6 modules together on real geometry and
//! verify each against analytic expectations.

use occt_core::geom::polygon_boolean::{polygon_difference, polygon_intersect, polygon_union};
use occt_core::gp::{GpAx1, GpAx3, GpCirc, GpDir, GpLin, GpPnt, GpTrsf, GpVec};
use occt_geom::extrema::{curve_surface_extrema, point_curve_extrema, surface_surface_extrema};
use occt_geom::{GeomCircle, GeomLine, GeomPlane, GeomSphere};
use occt_topo::brep_offset::{offset_polygon, offset_shell};
use occt_topo::gltf::{build_gltf_bin, GltfOptions};
use occt_topo::primitives::BRepPrimBox;

fn make_line(p: GpPnt, d: GpDir) -> GeomLine {
    GeomLine::new(GpLin::new(GpAx1::new(p, d)))
}

fn plane_at(z: f64) -> GeomPlane {
    let mut ax = GpAx3::standard();
    ax.set_location(GpPnt::new(0.0, 0.0, z));
    GeomPlane::new(occt_core::gp::GpPln::new(ax))
}

#[test]
fn polygon_boolean_ops_match_areas() {
    use occt_core::gp::GpPnt2d;
    let rect = |x0: f64, y0: f64, x1: f64, y1: f64| vec![
        GpPnt2d::new(x0, y0),
        GpPnt2d::new(x1, y0),
        GpPnt2d::new(x1, y1),
        GpPnt2d::new(x0, y1),
    ];
    let a = rect(0.0, 0.0, 2.0, 2.0);
    let b = rect(1.0, 1.0, 3.0, 3.0);
    let area = |p: &[GpPnt2d]| occt_core::geom::polygon_boolean::signed_area2d(p).abs();
    let u: f64 = polygon_union(&a, &b).iter().map(|p| area(p)).sum();
    let i: f64 = polygon_intersect(&a, &b).iter().map(|p| area(p)).sum();
    let d: f64 = polygon_difference(&a, &b).iter().map(|p| area(p)).sum();
    // Union 7, intersection 1, difference 3 for the two unit-overlap squares.
    assert!((u - 7.0).abs() < 1e-8, "union {u}");
    assert!((i - 1.0).abs() < 1e-8, "intersect {i}");
    assert!((d - 3.0).abs() < 1e-8, "difference {d}");
}

#[test]
fn geom_extrema_distances() {
    // Point to line: (3,4,0) to x-axis → 4.
    let line = make_line(GpPnt::new(0.0, 0.0, 0.0), GpDir::new(1.0, 0.0, 0.0).unwrap());
    let e = point_curve_extrema(&line, &GpPnt::new(3.0, 4.0, 0.0));
    assert!((e.distance - 4.0).abs() < 1e-7, "point-line {}", e.distance);

    // Curve–surface: line y=3 vs unit sphere at origin → 2.
    let sphere = GeomSphere::new(occt_core::gp::GpSphere::new(GpAx3::standard(), 1.0).unwrap());
    let l2 = make_line(GpPnt::new(0.0, 3.0, 0.0), GpDir::new(1.0, 0.0, 0.0).unwrap());
    let cs = curve_surface_extrema(&l2, &sphere, 16);
    assert!((cs.distance - 2.0).abs() < 1e-4, "curve-surface {}", cs.distance);

    // Surface–surface: sphere vs plane z=5 → 4.
    let plane = plane_at(5.0);
    let ss = surface_surface_extrema(&sphere, &plane, 4);
    assert!((ss.distance - 4.0).abs() < 1e-4, "surface-surface {}", ss.distance);
}

#[test]
fn offset_polygon_and_shell() {
    use occt_core::gp::GpPnt2d;
    let plane = occt_core::gp::GpPln::new(GpAx3::standard());
    let sq: Vec<GpPnt> = vec![
        GpPnt2d::new(-1.0, -1.0),
        GpPnt2d::new(1.0, -1.0),
        GpPnt2d::new(1.0, 1.0),
        GpPnt2d::new(-1.0, 1.0),
    ]
    .iter()
    .map(|p| GpPnt::new(p.x(), p.y(), 0.0))
    .collect();
    let out = offset_polygon(&sq, &plane, 0.5).expect("offset polygon");
    let mx = out.iter().map(|p| p.x().abs()).fold(0.0, f64::max);
    assert!((mx - 1.5).abs() < 1e-9, "offset square max x {mx}");

    // Box spans [0,2]³; offset outward by +0.5 → max x = 2.5.
    let b = BRepPrimBox::make_box(2.0, 2.0, 2.0);
    let grown = offset_shell(&b.solid.0, 0.5).expect("offset shell");
    let v = occt_topo::topo_tools_full::vertices_of(&grown);
    let mx = v.iter().map(|vv| occt_topo::brep_tool::BRepTool::vertex_point(vv).x()).fold(0.0, f64::max);
    assert!((mx - 2.5).abs() < 1e-9, "grown box max x {mx}");
}

#[test]
fn bincaf_roundtrip_preserves_tree() {
    use occt_topo::bincaf::{deserialize_bincaf, serialize_bincaf, BinXcaf, BinXcafEntry, XcafAttribute};
    let mut doc = BinXcaf::default();
    doc.root.attributes.push(XcafAttribute { kind: "name".into(), value: "root-part".into() });
    let b = BRepPrimBox::make_box(1.0, 2.0, 3.0);
    doc.root.shape = Some(b.solid.0.clone());
    let mut child = BinXcafEntry::default();
    child.attributes.push(XcafAttribute { kind: "color".into(), value: "1,0,0".into() });
    let s = BRepPrimBox::make_box(0.5, 0.5, 0.5);
    child.shape = Some(s.solid.0.clone());
    doc.root.children.push(child);

    let bytes = serialize_bincaf(&doc).expect("serialize");
    let back = deserialize_bincaf(&bytes).expect("deserialize");
    assert_eq!(back.root.attributes[0].value, "root-part");
    assert_eq!(back.root.children.len(), 1);
    assert_eq!(back.root.children[0].attributes[0].kind, "color");
    // Solid topology preserved.
    let nf = occt_topo::topo_tools_full::faces_of(back.root.shape.as_ref().unwrap()).len();
    assert_eq!(nf, 6, "box face count after round-trip");
}

#[test]
fn gltf_bin_is_wellformed() {
    let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
    let (json, bin) = build_gltf_bin(&b.solid.0, 0.2, &GltfOptions::default()).expect("gltf bin");
    assert!(json.contains("\"version\": \"2.0\""));
    assert!(json.contains("\"meshes\""));
    assert!(json.contains("\"scenes\""));
    assert!(!bin.is_empty());
    let marker = "\"byteLength\": ";
    let pos = json.find(marker).map(|p| p + marker.len()).expect("byteLength present");
    let rest = &json[pos..];
    let end = rest.find([',', '}']).map(|p| pos + p).unwrap_or(json.len());
    let reported: usize = json[pos..end].trim().parse().expect("numeric byteLength");
    assert_eq!(reported, bin.len(), "byteLength matches bin");
}
