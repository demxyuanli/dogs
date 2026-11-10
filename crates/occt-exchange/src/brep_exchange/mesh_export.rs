use occt_core::io::obj::{ObjFace, ObjMesh};
use occt_core::io::ply::PlyMesh;
use occt_core::io::stl::StlMesh;
use occt_topo::mesh::ShapeMesh;
use occt_topo::shape::TopoShape;

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
/// shared edges and periodic seams. A failed tessellation yields no triangles:
/// OCCT has no second mesher (audit A17/A18/A19).
///
/// ponytail: B-spline faces whose STEP `SURFACE_CURVE` pcurves are dropped
/// (`step.rs`) get a wrong boundary, so the Delaunay bbox drifts past `EXACT_TOL`
/// on `Shape.step`/`Shape-2.step` — a STEP-import gap, not a BRepMesh one.
/// `Prs3d::GetDeflection(shape, drawer)` (`Prs3d.hxx:82-103`):
/// `BRepBndLib::Add(shape, box, false)` then
/// `maxComp(diag) * DeviationCoefficient * 4`, coefficient default 0.001.
pub fn prs3d_get_deflection(shape: &TopoShape, maximal_chordial_deviation: f64) -> f64 {
    const DEVIATION_COEFFICIENT: f64 = 0.001;
    const CONFUSION: f64 = 1e-7;
    let b = occt_topo::brep_bnd_lib::shape_bnd_box(shape);
    if b.is_void() {
        return maximal_chordial_deviation;
    }
    let b = if b.is_open() {
        if !b.has_finite_part() {
            return maximal_chordial_deviation;
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

fn export_mesh(shape: &TopoShape, maximal_chordial_deviation: f64) -> ShapeMesh {
    // `RWObj_CafWriter` writes the vis triangulation
    // (`StdPrs_ToolTriangulatedShape::Tessellate` →
    // `BRepMesh_DiscretFactory::Discret(shape, GetDeflection, DeviationAngle)`).
    // Drawer defaults: 20 deg (`Prs3d_Drawer.cxx:96`) and relative
    // `GetDeflection` (`Prs3d.hxx:71`). The argument is the drawer's
    // `MaximalChordialDeviation` — it is **not** a linear deflection and only
    // matters for a void/unbounded bounding box (`Prs3d.hxx:82-103`).
    //
    // A failed tessellation yields an **empty** mesh: OCCT's OBJ writer reads
    // the triangulation already stored on the shape
    // (`RWMesh_FaceIterator.cxx:87`) and **skips** a face without one (`:89`) —
    // it never switches to another mesher. The two silent fallbacks that used to
    // sit here (legacy quadtree `brepmesh::incremental_mesh`, UV-grid
    // `shape_mesh::mesh_shape`) were non-OCCT algorithms and are removed
    // (audit A17); `brepmesh` itself is deleted as of A18/A19.
    const PRS3D_DEV_ANGLE: f64 = 20.0 * std::f64::consts::PI / 180.0;
    let lin = prs3d_get_deflection(shape, maximal_chordial_deviation);
    occt_topo::meshing::incremental_mesh::IncrementalMesh::from_deflection(
        shape,
        lin,
        false,
        PRS3D_DEV_ANGLE,
    )
    .mesh()
    .cloned()
    .unwrap_or_else(|| ShapeMesh {
        vertices: Vec::new(),
        triangles: Vec::new(),
        source_shape: occt_topo::abs::ShapeType::Shape,
    })
}

/// Export a shape to Wavefront OBJ text.
///
/// Vertices stay per face (`RWObj_CafWriter` / `RWMesh_FaceIterator`):
/// a box is 24 nodes / 12 triangles, not a globally welded 8-node mesh.
pub fn brep_to_obj(shape: &TopoShape, maximal_chordial_deviation: f64) -> String {
    let mesh = export_mesh(shape, maximal_chordial_deviation);
    occt_core::io::obj::obj_mesh_to_string(&shape_mesh_to_obj(&mesh))
}

/// Export a shape to ASCII STL text.
pub fn brep_to_stl_ascii(shape: &TopoShape, maximal_chordial_deviation: f64) -> String {
    let mesh = export_mesh(shape, maximal_chordial_deviation);
    occt_core::io::stl::write_ascii_stl(&shape_mesh_to_stl(&mesh))
}

/// Export a shape to binary STL bytes (vertices welded first).
pub fn brep_to_stl_binary(shape: &TopoShape, maximal_chordial_deviation: f64) -> Vec<u8> {
    let mut mesh = export_mesh(shape, maximal_chordial_deviation);
    occt_topo::shape_mesh::weld_vertices(&mut mesh, 1e-9);
    occt_core::io::stl::write_binary_stl(&shape_mesh_to_stl(&mesh))
}

/// Export a shape to ASCII PLY text.
///
/// `RWPly_PlyWriterContext` writes each **face's own node block** with a running
/// `myVertOffset` (`RWPly_PlyWriterContext.cxx:276`, `:300`) — it never merges
/// nodes across faces, so this path must **not** weld: a box is 24 vertex
/// records (6 faces × 4), exactly like the OBJ path ([`brep_to_obj`]). The
/// previous body called `weld_vertices(1e-9)` first, producing a shared 8-node
/// mesh that OCCT never writes (audit A26).
pub fn brep_to_ply(shape: &TopoShape, maximal_chordial_deviation: f64) -> String {
    let mesh = export_mesh(shape, maximal_chordial_deviation);
    occt_core::io::ply::ply_mesh_to_string(&shape_mesh_to_ply(&mesh))
}

/// Write a shape to an OBJ file (per-face nodes, same as [`brep_to_obj`]).
pub fn brep_write_obj(path: &str, shape: &TopoShape, maximal_chordial_deviation: f64) -> std::io::Result<()> {
    let mesh = export_mesh(shape, maximal_chordial_deviation);
    occt_core::io::obj::write_obj_file(path, &shape_mesh_to_obj(&mesh))
}

// ---------------------------------------------------------------------------
// RWObj_CafWriter: OBJ + MTL for a document (`BRepModel`)
// ---------------------------------------------------------------------------

pub fn brep_write_stl(path: &str, shape: &TopoShape, maximal_chordial_deviation: f64) -> std::io::Result<()> {
    std::fs::write(path, brep_to_stl_binary(shape, maximal_chordial_deviation))
}

/// Write a shape to an ASCII PLY file.
pub fn brep_write_ply(path: &str, shape: &TopoShape, maximal_chordial_deviation: f64) -> std::io::Result<()> {
    std::fs::write(path, brep_to_ply(shape, maximal_chordial_deviation))
}

