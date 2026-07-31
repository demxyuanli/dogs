//! PLY (Polygon File Format) ASCII reader/writer. Source: `RWPLY` (simplified).

use std::fs;

use crate::bnd::BndBox;
use crate::gp::GpPnt;
use crate::io::obj::ObjMesh;
use crate::poly::triangulation::Triangle;
use crate::poly::Triangulation;

/// PLY mesh: vertices plus polygon faces (0-based vertex indices).
#[derive(Debug, Clone, Default)]
pub struct PlyMesh {
    pub vertices: Vec<GpPnt>,
    pub faces: Vec<Vec<usize>>,
}

/// Parse ASCII PLY: header with `element vertex N` / `element face M` / `end_header`,
/// then N vertex lines `x y z` and M face lines `k i1 i2 ... ik`.
pub fn parse_ply(content: &str) -> Result<PlyMesh, String> {
    let mut lines = content.lines();
    let mut vertex_count = None;
    let mut face_count = None;
    let mut seen_end = false;
    for line in &mut lines {
        let t = line.trim();
        if t == "end_header" {
            seen_end = true;
            break;
        }
        let mut parts = t.split_whitespace();
        if parts.next() == Some("element") {
            let kind = parts.next().unwrap_or("");
            let n: usize = parts
                .next()
                .ok_or("element missing count")?
                .parse()
                .map_err(|_| format!("bad element count: {t}"))?;
            match kind {
                "vertex" => vertex_count = Some(n),
                "face" => face_count = Some(n),
                _ => {}
            }
        }
    }
    if !seen_end {
        return Err("missing end_header".into());
    }
    let vc = vertex_count.ok_or("no vertex element in header")?;
    let fc = face_count.unwrap_or(0);

    let mut mesh = PlyMesh::default();
    for line in lines {
        let t = line.trim();
        if t.is_empty() {
            continue;
        }
        if mesh.vertices.len() < vc {
            let mut p = t.split_whitespace();
            let x: f64 = p.next().ok_or("vertex missing x")?.parse().map_err(|e| format!("bad vertex x: {e}"))?;
            let y: f64 = p.next().ok_or("vertex missing y")?.parse().map_err(|e| format!("bad vertex y: {e}"))?;
            let z: f64 = p.next().ok_or("vertex missing z")?.parse().map_err(|e| format!("bad vertex z: {e}"))?;
            mesh.vertices.push(GpPnt::new(x, y, z));
        } else if mesh.faces.len() < fc {
            let mut p = t.split_whitespace();
            let k: usize = p.next().ok_or("face missing vertex count")?.parse().map_err(|e| format!("bad face count: {e}"))?;
            let mut face = Vec::with_capacity(k);
            for _ in 0..k {
                let idx: usize = p.next().ok_or("face missing index")?.parse().map_err(|e| format!("bad face index: {e}"))?;
                face.push(idx);
            }
            mesh.faces.push(face);
        } else {
            break;
        }
    }
    if mesh.vertices.len() != vc {
        return Err(format!("expected {vc} vertices, got {}", mesh.vertices.len()));
    }
    if mesh.faces.len() != fc {
        return Err(format!("expected {fc} faces, got {}", mesh.faces.len()));
    }
    Ok(mesh)
}

/// Read an ASCII PLY file.
pub fn read_ply_file(path: &str) -> Result<PlyMesh, String> {
    let content = fs::read_to_string(path).map_err(|e| format!("cannot read {}: {e}", path))?;
    parse_ply(&content)
}

/// Serialize mesh to ASCII PLY.
pub fn write_ply(mesh: &PlyMesh) -> String {
    let mut s = String::with_capacity(64 + mesh.vertices.len() * 32 + mesh.faces.len() * 24);
    s.push_str("ply\nformat ascii 1.0\n");
    s.push_str(&format!("element vertex {}\n", mesh.vertices.len()));
    s.push_str("property float x\nproperty float y\nproperty float z\n");
    s.push_str(&format!("element face {}\n", mesh.faces.len()));
    s.push_str("property list uchar int vertex_indices\nend_header\n");
    for v in &mesh.vertices {
        s.push_str(&format!("{} {} {}\n", v.x(), v.y(), v.z()));
    }
    for f in &mesh.faces {
        s.push_str(&f.len().to_string());
        for &i in f {
            s.push_str(&format!(" {i}"));
        }
        s.push('\n');
    }
    s
}

/// Build a `Triangulation`, fan-triangulating polygon faces (triangles pass through).
pub fn to_triangulation(mesh: &PlyMesh) -> Triangulation {
    let nodes = mesh.vertices.clone();
    let mut tris = Vec::new();
    for f in &mesh.faces {
        if f.len() < 3 {
            continue;
        }
        for i in 1..f.len() - 1 {
            tris.push(Triangle::new(f[0], f[i], f[i + 1]));
        }
    }
    Triangulation::new(nodes, tris)
}

/// Convert an OBJ mesh (faces already carry 0-based vertex indices).
pub fn from_obj(mesh: &ObjMesh) -> PlyMesh {
    PlyMesh {
        vertices: mesh.vertices.clone(),
        faces: mesh.faces.iter().map(|f| f.v.iter().map(|&i| i as usize).collect()).collect(),
    }
}

pub fn vertex_count(mesh: &PlyMesh) -> usize {
    mesh.vertices.len()
}

pub fn face_count(mesh: &PlyMesh) -> usize {
    mesh.faces.len()
}

/// Axis-aligned bounding box of all mesh vertices.
pub fn ply_bbox(mesh: &PlyMesh) -> BndBox {
    let mut b = BndBox::new();
    for v in &mesh.vertices {
        b.add_point(v);
    }
    b
}

#[cfg(test)]
mod tests {
    use super::*;

    fn quad_ply() -> PlyMesh {
        PlyMesh {
            vertices: vec![
                GpPnt::new(0.0, 0.0, 0.0),
                GpPnt::new(1.0, 0.0, 0.0),
                GpPnt::new(1.0, 1.0, 0.0),
                GpPnt::new(0.0, 1.0, 0.0),
            ],
            faces: vec![vec![0, 1, 2], vec![0, 2, 3]],
        }
    }

    #[test]
    fn parse_ply_two_faces() {
        let content = "ply
format ascii 1.0
element vertex 4
property float x
property float y
property float z
element face 2
property list uchar int vertex_indices
end_header
0 0 0
1 0 0
1 1 0
0 1 0
3 0 1 2
3 0 2 3
";
        let m = parse_ply(content).unwrap();
        assert_eq!(vertex_count(&m), 4);
        assert_eq!(face_count(&m), 2);
        assert_eq!(m.faces[1], vec![0, 2, 3]);
        assert_eq!(m.vertices[2].x(), 1.0);
    }

    #[test]
    fn write_parse_roundtrip() {
        let m = quad_ply();
        let s = write_ply(&m);
        let back = parse_ply(&s).unwrap();
        assert_eq!(back.vertices.len(), 4);
        assert_eq!(back.faces.len(), 2);
        assert_eq!(back.vertices[3].y(), 1.0);
        assert_eq!(back.faces[0], vec![0, 1, 2]);
    }

    #[test]
    fn quad_triangulates_to_two() {
        let m = PlyMesh { faces: vec![vec![0, 1, 2, 3]], ..quad_ply() };
        let t = to_triangulation(&m);
        assert_eq!(t.nodes.len(), 4);
        assert_eq!(t.triangles.len(), 2);
        assert_eq!(t.triangles[1].n2, 3);
    }
}
