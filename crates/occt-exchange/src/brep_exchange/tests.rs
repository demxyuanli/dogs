use super::*;
use occt_core::gp::{GpAx1, GpAx3, GpDir, GpPln, GpPnt};
use occt_topo::brep_tool::BRepTool;
use occt_topo::builder::TopoBuilder;
use occt_topo::shape::{TopoShape, Vertex};

/// Build a unit cube [0,1]^3 as a Solid with real geometry.
fn unit_box() -> TopoShape {
    let b = TopoBuilder::new();
    let mk = |x: f64, y: f64, z: f64| b.make_vertex(GpPnt::new(x, y, z), 0.0);
    let p0 = mk(0.,0.,0.); let p1 = mk(1.,0.,0.); let p2 = mk(1.,1.,0.); let p3 = mk(0.,1.,0.);
    let p4 = mk(0.,0.,1.); let p5 = mk(1.,0.,1.); let p6 = mk(1.,1.,1.); let p7 = mk(0.,1.,1.);

    let face = |c: [&Vertex; 4], normal: GpDir| {
        let e = [
            b.make_edge_segment(&BRepTool::vertex_point(c[0]), &BRepTool::vertex_point(c[1])),
            b.make_edge_segment(&BRepTool::vertex_point(c[1]), &BRepTool::vertex_point(c[2])),
            b.make_edge_segment(&BRepTool::vertex_point(c[2]), &BRepTool::vertex_point(c[3])),
            b.make_edge_segment(&BRepTool::vertex_point(c[3]), &BRepTool::vertex_point(c[0])),
        ];
        let wire = b.make_wire(&e);
        let pln = GpPln::new(GpAx3::from_ax1(&GpAx1::new(BRepTool::vertex_point(c[0]), normal)));
        let mut f = b.make_face_plane(&pln);
        b.add_wire(&mut f, &wire);
        f
    };

    let nz = GpDir::new(0.,0.,-1.).unwrap(); let pz = GpDir::new(0.,0.,1.).unwrap();
    let nx = GpDir::new(-1.,0.,0.).unwrap(); let px = GpDir::new(1.,0.,0.).unwrap();
    let ny = GpDir::new(0.,-1.,0.).unwrap(); let py = GpDir::new(0.,1.,0.).unwrap();

    let faces = vec![
        face([&p0, &p1, &p2, &p3], nz),
        face([&p4, &p5, &p6, &p7], pz),
        face([&p0, &p3, &p7, &p4], nx),
        face([&p1, &p5, &p6, &p2], px),
        face([&p0, &p1, &p5, &p4], ny),
        face([&p3, &p2, &p6, &p7], py),
    ];
    let shell = b.make_shell(&faces);
    let solid = b.make_solid(&[shell]);
    solid.into()
}

#[test]
fn obj_export_has_vertices_and_faces() {
    let shape = unit_box();
    let obj = brep_to_obj(&shape, 0.1);
    assert!(obj.contains("v 0.000000000 0.000000000 0.000000000"), "origin corner present");
    assert!(obj.lines().filter(|l| l.starts_with("f ")).count() > 0);
}

#[test]
fn stl_ascii_export() {
    let shape = unit_box();
    let stl = brep_to_stl_ascii(&shape, 0.1);
    assert!(stl.contains("solid"));
    assert!(stl.contains("endsolid"));
}

#[test]
fn stl_binary_and_ply_export() {
    let shape = unit_box();
    let bytes = brep_to_stl_binary(&shape, 0.1);
    assert!(bytes.len() >= 84);
    assert_eq!((bytes.len() - 84) % 50, 0);
    let ply = brep_to_ply(&shape, 0.1);
    assert!(ply.starts_with("ply"));
}

#[test]
fn write_obj_file_roundtrip() {
    let shape = unit_box();
    let path = { let d = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("target").join("test-out"); let _ = std::fs::create_dir_all(&d); d }.join(format!("occt_brep_obj_{}.obj", std::process::id()));
    let p = path.to_str().unwrap().to_string();
    brep_write_obj(&p, &shape, 0.1).unwrap();
    let m = occt_core::io::obj::read_obj_file(&p).unwrap();
    assert!(m.vertices.len() > 0);
    assert!(m.faces.len() > 0);
    let _ = std::fs::remove_file(&p);
}
