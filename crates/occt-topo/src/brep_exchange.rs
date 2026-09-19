//! Module 2: brep_exchange — export BRep shapes to OBJ / STL / PLY.
//! Source: `RWObj`, `RWStl`, `RWPLY` (simplified OCCT writers).

use crate::mesh::ShapeMesh;
use crate::shape::TopoShape;
use occt_core::io::obj::{ObjFace, ObjMesh};
use occt_core::io::ply::PlyMesh;
use occt_core::io::stl::StlMesh;

/// Convert a `ShapeMesh` to an OBJ mesh (each triangle becomes an OBJ face).
fn shape_mesh_to_obj(mesh: &ShapeMesh) -> ObjMesh {
    ObjMesh {
        vertices: mesh.vertices.clone(),
        texcoords: Vec::new(),
        normals: Vec::new(),
        faces: mesh.triangles.iter().map(|t| ObjFace {
            v: vec![t.n0 as i32, t.n1 as i32, t.n2 as i32],
            vt: None,
            vn: None,
        }).collect(),
    }
}

/// Convert a `ShapeMesh` to an STL mesh (independent per-triangle vertices).
fn shape_mesh_to_stl(mesh: &ShapeMesh) -> StlMesh {
    StlMesh {
        triangles: mesh.triangles.iter().map(|t| [
            mesh.vertices[t.n0],
            mesh.vertices[t.n1],
            mesh.vertices[t.n2],
        ]).collect(),
        normals: Vec::new(),
    }
}

/// Convert a `ShapeMesh` to a PLY mesh.
fn shape_mesh_to_ply(mesh: &ShapeMesh) -> PlyMesh {
    PlyMesh {
        vertices: mesh.vertices.clone(),
        faces: mesh.triangles.iter().map(|t| vec![t.n0, t.n1, t.n2]).collect(),
    }
}

/// Mesh a shape for export using the OCCT `BRepMesh_IncrementalMesh` pipeline:
/// the Delaunay path (`meshing::incremental_mesh_to_shape_mesh`) shares each
/// edge's 3D polyline between its adjacent faces, so the output is watertight on
/// shared edges and periodic seams. Falls back to the quadtree subdivision path
/// (`brepmesh`) when the Delaunay pipeline errors.
///
/// ponytail: B-spline faces whose STEP `SURFACE_CURVE` pcurves are dropped
/// (`step.rs`) get a wrong boundary, so the Delaunay bbox drifts past `EXACT_TOL`
/// on `Shape.step`/`Shape-2.step` — a STEP-import gap, not a BRepMesh one.
/// `Prs3d::GetDeflection(shape, drawer)` (`Prs3d.hxx:82-103`):
/// `BRepBndLib::Add(shape, box, false)` then
/// `maxComp(diag) * DeviationCoefficient * 4`, coefficient default 0.001.
pub fn prs3d_get_deflection(shape: &TopoShape, maximal_chordial: f64) -> f64 {
    const DEVIATION_COEFFICIENT: f64 = 0.001;
    const CONFUSION: f64 = 1e-7;
    let b = crate::brep_bnd_lib::shape_bnd_box(shape);
    if b.is_void() {
        return maximal_chordial;
    }
    let b = if b.is_open() {
        if !b.has_finite_part() {
            return maximal_chordial;
        }
        b.finite_part()
    } else {
        b
    };
    let mn = b.corner_min();
    let mx = b.corner_max();
    let max_comp = (mx.x() - mn.x())
        .max(mx.y() - mn.y())
        .max(mx.z() - mn.z());
    (max_comp * DEVIATION_COEFFICIENT * 4.0).max(CONFUSION)
}

fn export_mesh(shape: &TopoShape, deflection: f64) -> ShapeMesh {
    // `RWObj_CafWriter` writes the vis triangulation
    // (`StdPrs_ToolTriangulatedShape::Tessellate` →
    // `BRepMesh_DiscretFactory::Discret(shape, GetDeflection, DeviationAngle)`).
    // Drawer defaults: 20 deg (`Prs3d_Drawer.cxx:96`) and relative
    // `GetDeflection` (`Prs3d.hxx:71`). The `deflection` argument is only
    // the void-box fallback (`MaximalChordialDeviation`).
    const PRS3D_DEV_ANGLE: f64 = 20.0 * std::f64::consts::PI / 180.0;
    let lin = prs3d_get_deflection(shape, deflection);
    crate::meshing::incremental_mesh::IncrementalMesh::from_deflection(
        shape,
        lin,
        false,
        PRS3D_DEV_ANGLE,
    )
    .mesh()
    .cloned()
    .or_else(|| crate::brepmesh::incremental_mesh(shape, deflection).ok().map(|im| im.mesh))
    .unwrap_or_else(|| crate::shape_mesh::mesh_shape(shape, deflection))
}

/// Export a shape to Wavefront OBJ text.
///
/// Vertices stay per face (`RWObj_CafWriter` / `RWMesh_FaceIterator`):
/// a box is 24 nodes / 12 triangles, not a globally welded 8-node mesh.
pub fn brep_to_obj(shape: &TopoShape, deflection: f64) -> String {
    let mesh = export_mesh(shape, deflection);
    occt_core::io::obj::write_obj(&shape_mesh_to_obj(&mesh))
}

/// Export a shape to ASCII STL text.
pub fn brep_to_stl_ascii(shape: &TopoShape, deflection: f64) -> String {
    let mesh = export_mesh(shape, deflection);
    occt_core::io::stl::write_ascii_stl(&shape_mesh_to_stl(&mesh))
}

/// Export a shape to binary STL bytes (vertices welded first).
pub fn brep_to_stl_binary(shape: &TopoShape, deflection: f64) -> Vec<u8> {
    let mut mesh = export_mesh(shape, deflection);
    crate::shape_mesh::weld_vertices(&mut mesh, 1e-9);
    occt_core::io::stl::write_binary_stl(&shape_mesh_to_stl(&mesh))
}

/// Export a shape to ASCII PLY text (vertices welded first).
pub fn brep_to_ply(shape: &TopoShape, deflection: f64) -> String {
    let mut mesh = export_mesh(shape, deflection);
    crate::shape_mesh::weld_vertices(&mut mesh, 1e-9);
    occt_core::io::ply::write_ply(&shape_mesh_to_ply(&mesh))
}

/// Write a shape to an OBJ file (per-face nodes, same as [`brep_to_obj`]).
pub fn brep_write_obj(path: &str, shape: &TopoShape, deflection: f64) -> std::io::Result<()> {
    let mesh = export_mesh(shape, deflection);
    occt_core::io::obj::write_obj_file(path, &shape_mesh_to_obj(&mesh))
}

/// Write a shape to a (binary) STL file.
pub fn brep_write_stl(path: &str, shape: &TopoShape, deflection: f64) -> std::io::Result<()> {
    std::fs::write(path, brep_to_stl_binary(shape, deflection))
}

/// Write a shape to an ASCII PLY file.
pub fn brep_write_ply(path: &str, shape: &TopoShape, deflection: f64) -> std::io::Result<()> {
    std::fs::write(path, brep_to_ply(shape, deflection))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_tool::BRepTool;
    use crate::builder::TopoBuilder;
    use crate::shape::Vertex;
    use occt_core::gp::{GpAx1, GpAx3, GpDir, GpPln, GpPnt};

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
        let path = std::env::temp_dir().join(format!("occt_brep_obj_{}.obj", std::process::id()));
        let p = path.to_str().unwrap().to_string();
        brep_write_obj(&p, &shape, 0.1).unwrap();
        let m = occt_core::io::obj::read_obj_file(&p).unwrap();
        assert!(m.vertices.len() > 0);
        assert!(m.faces.len() > 0);
        let _ = std::fs::remove_file(&p);
    }
}
