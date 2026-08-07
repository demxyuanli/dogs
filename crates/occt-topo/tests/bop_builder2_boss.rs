//! Wave 6 gate: boss fuse via `bop_builder2` merges the box top hole with the
//! cylinder base into one solid.
//!
//! A boss (cylinder standing on the box top face) fused with the box must yield
//! **one** solid whose box-top is a 2-wire ring face (the hole the cylinder
//! base closes). Two cylinder fixtures cover both face representations:
//! - single-disc base (one 24-gon planar face, OCCT-style);
//! - triangular-fan base (faceted `shape_mesh_to_brep`, the brepfeat fixture).
//!
//! Regressed against OCCT `BOPAlgo_Builder::FillImagesSolids` → `FillIn3DParts`
//! + `BOPAlgo_SplitSolid`: a covering face (the base over the box-top hole) is
//! included as an internal face, and cross-solid pieces sharing a face collapse
//! to one union region.

use occt_core::gp::{GpAx3, GpDir, GpPln, GpPnt, GpVec};
use occt_topo::bop_builder2::{builder_bop, BoolOp2};
use occt_topo::builder::TopoBuilder;
use occt_topo::primitives::BRepPrimBox;
use occt_topo::shape::{Edge, Face, TopoShape};
use occt_topo::topo_tools_full::{faces_of, wires_of};

/// A faceted cylinder whose base and top are SINGLE 24-gon faces (OCCT-style).
fn single_disc_cylinder(radius: f64, height: f64, slices: usize) -> TopoShape {
    let b = TopoBuilder::new();
    let base_ax = GpAx3::new(GpPnt::new(0.0, 0.0, 0.0), GpDir::new(0.0, 0.0, 1.0).unwrap(),
        &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
    let top_ax = GpAx3::new(GpPnt::new(0.0, 0.0, height), GpDir::new(0.0, 0.0, -1.0).unwrap(),
        &GpDir::new(1.0, 0.0, 0.0).unwrap()).unwrap();
    let mut rim_bot: Vec<Edge> = Vec::new();
    let mut rim_top: Vec<Edge> = Vec::new();
    for i in 0..slices {
        let th = 2.0 * std::f64::consts::PI * i as f64 / slices as f64;
        let th2 = 2.0 * std::f64::consts::PI * (i + 1) as f64 / slices as f64;
        let a = GpPnt::new(radius * th.cos(), radius * th.sin(), 0.0);
        let c = GpPnt::new(radius * th2.cos(), radius * th2.sin(), 0.0);
        rim_bot.push(b.make_edge_segment(&a, &c));
        let at = GpPnt::new(radius * th.cos(), radius * th.sin(), height);
        let ct = GpPnt::new(radius * th2.cos(), radius * th2.sin(), height);
        rim_top.push(b.make_edge_segment(&ct, &at));
    }
    let mut side_faces: Vec<Face> = Vec::new();
    for i in 0..slices {
        let j = (i + 1) % slices;
        let th = 2.0 * std::f64::consts::PI * i as f64 / slices as f64;
        let thj = 2.0 * std::f64::consts::PI * j as f64 / slices as f64;
        let a = GpPnt::new(radius * th.cos(), radius * th.sin(), 0.0);
        let c = GpPnt::new(radius * thj.cos(), radius * thj.sin(), 0.0);
        let vert_i = b.make_edge_segment(&a, &GpPnt::new(a.x(), a.y(), height));
        let vert_j = b.make_edge_segment(&c, &GpPnt::new(c.x(), c.y(), height));
        let mut e1 = rim_bot[i].clone();
        e1.0.set_orientation(occt_topo::abs::Orientation::Forward);
        let mut e2 = vert_j.clone();
        e2.0.set_orientation(occt_topo::abs::Orientation::Forward);
        let mut e3 = rim_top[i].clone();
        e3.0.set_orientation(occt_topo::abs::Orientation::Forward);
        let mut e4 = vert_i.clone();
        e4.0.set_orientation(occt_topo::abs::Orientation::Reversed);
        let wire = b.make_wire(&[e1, e2, e3, e4]);
        let n = GpVec::new(c.x() - a.x(), c.y() - a.y(), 0.0).crossed(&GpVec::new(0.0, 0.0, height)).normalized();
        let ax = GpAx3::new(a, GpDir::from_vec(&n).unwrap(), &GpDir::new(0.0, 0.0, 1.0).unwrap()).unwrap();
        let mut face = b.make_face_plane(&GpPln::new(ax));
        b.add_wire(&mut face, &wire);
        side_faces.push(face);
    }
    let base_wire = b.make_wire(&rim_bot);
    let mut base = b.make_face_plane(&GpPln::new(base_ax));
    b.add_wire(&mut base, &base_wire);
    let top_wire = b.make_wire(&rim_top);
    let mut top = b.make_face_plane(&GpPln::new(top_ax));
    b.add_wire(&mut top, &top_wire);
    let mut all_faces = side_faces;
    all_faces.push(base);
    all_faces.push(top);
    b.make_solid(&[b.make_shell(&all_faces)]).0
}

/// Faceted cylinder with TRIANGULAR-fan base/top (the brepfeat::mesh_cylinder
/// fixture).
fn tri_fan_cylinder(radius: f64, height: f64, slices: usize) -> TopoShape {
    use occt_core::poly::triangulation::Triangle;
    let mut vertices = Vec::new();
    let mut triangles = Vec::new();
    for i in 0..=slices {
        let theta = 2.0 * std::f64::consts::PI * i as f64 / slices as f64;
        let (x, y) = (radius * theta.cos(), radius * theta.sin());
        vertices.push(GpPnt::new(x, y, 0.0));
        vertices.push(GpPnt::new(x, y, height));
    }
    for i in 0..slices {
        let (a, b, c, d) = (2 * i, 2 * i + 1, 2 * i + 2, 2 * i + 3);
        triangles.push(Triangle::new(a, b, c));
        triangles.push(Triangle::new(b, d, c));
    }
    let (top_idx, bot_idx) = (vertices.len(), vertices.len() + 1);
    vertices.push(GpPnt::new(0.0, 0.0, height));
    vertices.push(GpPnt::new(0.0, 0.0, 0.0));
    for i in 0..slices {
        let (tb, tt) = (2 * i + 1, 2 * i + 3);
        let (bb, bt) = (2 * i, 2 * i + 2);
        triangles.push(Triangle::new(top_idx, tt, tb));
        triangles.push(Triangle::new(bot_idx, bb, bt));
    }
    let mesh = occt_topo::mesh::ShapeMesh { vertices, triangles, source_shape: occt_topo::abs::ShapeType::Solid };
    occt_topo::mesh_to_brep::shape_mesh_to_brep(&mesh).solid.expect("closed").0
}

/// Fuse box + boss and assert one solid with a 2-wire ring box-top face.
fn assert_boss_one_solid(cyl: &TopoShape) {
    let box_s = BRepPrimBox::make_box(2.0, 2.0, 1.0);
    let cyl = occt_topo::shape_ops::translated_copy(cyl, &GpVec::new(1.0, 1.0, 1.0)).unwrap();
    let result = builder_bop(&[box_s.solid.0.clone()], &[cyl], BoolOp2::Fuse).expect("fuse");
    let counts = occt_topo::topo_tools_full::shape_counts(&result);
    let n_solids = counts.get(&occt_topo::abs::ShapeType::Solid).copied().unwrap_or(0);
    assert_eq!(n_solids, 1, "boss fuse must merge box+cylinder into one solid");
    let ring = faces_of(&result).iter().filter(|f| wires_of(&f.0).len() == 2).count();
    assert!(ring >= 1, "result must contain the 2-wire ring box-top face");
}

#[test]
fn boss_single_disc_base_merges_one_solid() {
    let cyl = single_disc_cylinder(0.25, 0.8, 24);
    assert_boss_one_solid(&cyl);
}

#[test]
fn boss_tri_fan_base_merges_one_solid() {
    let cyl = tri_fan_cylinder(0.25, 0.8, 24);
    assert_boss_one_solid(&cyl);
}
