//! Wavefront OBJ file reader/writer. Source: `RWObj_Reader` / `RWObj_Writer` (simplified)
//! Minimal but complete OBJ parser supporting v/vt/vn/f with negative indices.
use std::io::{self, Write};
use crate::gp::GpPnt;
use crate::poly::Triangulation;

/// Parsed OBJ mesh: vertices, texture coords, normals, faces.
#[derive(Debug, Clone, Default)]
pub struct ObjMesh {
    pub vertices: Vec<GpPnt>,
    pub texcoords: Vec<(f64, f64)>,
    pub normals: Vec<GpPnt>,
    pub faces: Vec<ObjFace>,
}

/// A face: indices into vertices (0-based), optional texcoord/normal indices.
#[derive(Debug, Clone)]
pub struct ObjFace {
    pub v: Vec<i32>, // vertex indices (0-based after parsing)
    pub vt: Option<Vec<i32>>,
    pub vn: Option<Vec<i32>>,
}

impl ObjMesh {
    /// Parse OBJ text content. Returns Err on malformed face or unknown line.
    pub fn parse(content: &str) -> Result<Self, String> {
        let mut mesh = ObjMesh::default();
        for (lineno, raw) in content.lines().enumerate() {
            let line = raw.trim();
            if line.is_empty() || line.starts_with('#') { continue; }
            let mut parts = line.split_whitespace();
            let tag = parts.next().unwrap_or("");
            match tag {
                "v" => {
                    let coords: Vec<f64> = parts.filter_map(|p| p.parse().ok()).collect();
                    if coords.len() < 3 { return Err(format!("line {}: bad vertex", lineno+1)); }
                    mesh.vertices.push(GpPnt::new(coords[0], coords[1], coords[2]));
                }
                "vt" => {
                    let c: Vec<f64> = parts.filter_map(|p| p.parse().ok()).collect();
                    if c.len() < 2 { return Err(format!("line {}: bad texcoord", lineno+1)); }
                    mesh.texcoords.push((c[0], c[1]));
                }
                "vn" => {
                    let c: Vec<f64> = parts.filter_map(|p| p.parse().ok()).collect();
                    if c.len() < 3 { return Err(format!("line {}: bad normal", lineno+1)); }
                    mesh.normals.push(GpPnt::new(c[0], c[1], c[2]));
                }
                "f" => {
                    let mut face = ObjFace { v: Vec::new(), vt: None, vn: None };
                    let mut vts = Vec::new();
                    let mut vns = Vec::new();
                    for tok in parts {
                        // Formats: v, v/vt, v//vn, v/vt/vn
                        let segs: Vec<&str> = tok.split('/').collect();
                        let vi: i32 = segs[0].parse().map_err(|_| format!("line {}: bad face index", lineno+1))?;
                        face.v.push(resolve_index(vi, mesh.vertices.len())?);
                        if segs.len() > 1 && !segs[1].is_empty() {
                            vts.push(resolve_index(segs[1].parse::<i32>().unwrap(), mesh.texcoords.len())?);
                        }
                        if segs.len() > 2 && !segs[2].is_empty() {
                            vns.push(resolve_index(segs[2].parse::<i32>().unwrap(), mesh.normals.len())?);
                        }
                    }
                    if !vts.is_empty() { face.vt = Some(vts); }
                    if !vns.is_empty() { face.vn = Some(vns); }
                    mesh.faces.push(face);
                }
                _ => { /* skip unsupported: o, g, s, usemtl, mtllib */ }
            }
        }
        Ok(mesh)
    }

    /// Convert triangulated faces into a Triangulation (fan triangulation).
    pub fn to_triangulation(&self) -> Triangulation {
        let mut triangles = Vec::new();
        for face in &self.faces {
            for i in 1..face.v.len().saturating_sub(1) {
                triangles.push(crate::poly::triangulation::Triangle::new(
                    face.v[0] as usize, face.v[i] as usize, face.v[i+1] as usize,
                ));
            }
        }
        Triangulation::new(self.vertices.clone(), triangles)
    }

    /// Compute per-vertex normals by averaging adjacent face normals.
    pub fn compute_normals(&mut self) {
        let n = self.vertices.len();
        let mut accum = vec![GpPnt::zero(); n];
        for face in &self.faces {
            if face.v.len() < 3 { continue; }
            let a = &self.vertices[face.v[0] as usize];
            let b = &self.vertices[face.v[1] as usize];
            let c = &self.vertices[face.v[2] as usize];
            let nrm = b.coord.subtracted(&a.coord).crossed(&c.coord.subtracted(&a.coord));
            let nrm = if nrm.modulus() > 1e-30 { GpPnt::from_xyz(&nrm.divided(nrm.modulus())) } else { GpPnt::zero() };
            for &vi in &face.v { accum[vi as usize] = GpPnt::new(accum[vi as usize].x()+nrm.x(), accum[vi as usize].y()+nrm.y(), accum[vi as usize].z()+nrm.z()); }
        }
        self.normals = accum.into_iter().map(|p| {
            let m = (p.x()*p.x()+p.y()*p.y()+p.z()*p.z()).sqrt();
            if m > 1e-30 { GpPnt::new(p.x()/m, p.y()/m, p.z()/m) } else { GpPnt::zero() }
        }).collect();
    }
}

/// Resolve 1-based (or negative = from-end) OBJ index to 0-based.
fn resolve_index(i: i32, len: usize) -> Result<i32, String> {
    let resolved = if i > 0 { i - 1 } else { len as i32 + i };
    if resolved < 0 || resolved >= len as i32 {
        Err(format!("index {i} out of range (len {len})"))
    } else { Ok(resolved) }
}

/// Write mesh to OBJ text.
pub fn write_obj(mesh: &ObjMesh) -> String {
    let mut out = String::new();
    for v in &mesh.vertices { out.push_str(&format!("v {:.9} {:.9} {:.9}\n", v.x(), v.y(), v.z())); }
    for (u, v) in &mesh.texcoords { out.push_str(&format!("vt {:.9} {:.9}\n", u, v)); }
    for n in &mesh.normals { out.push_str(&format!("vn {:.9} {:.9} {:.9}\n", n.x(), n.y(), n.z())); }
    for f in &mesh.faces {
        match (&f.vt, &f.vn) {
            (None, None) => out.push_str(&format!("f {}\n", f.v.iter().map(|i| (i+1).to_string()).collect::<Vec<_>>().join(" "))),
            (Some(vt), None) => out.push_str(&format!("f {}\n", (0..f.v.len()).map(|k| format!("{}/{}", f.v[k]+1, vt[k]+1)).collect::<Vec<_>>().join(" "))),
            (None, Some(vn)) => out.push_str(&format!("f {}\n", (0..f.v.len()).map(|k| format!("{}//{}", f.v[k]+1, vn[k]+1)).collect::<Vec<_>>().join(" "))),
            (Some(vt), Some(vn)) => out.push_str(&format!("f {}\n", (0..f.v.len()).map(|k| format!("{}/{}/{}", f.v[k]+1, vt[k]+1, vn[k]+1)).collect::<Vec<_>>().join(" "))),
        }
    }
    out
}

/// Write OBJ to file.
pub fn write_obj_file(path: &str, mesh: &ObjMesh) -> io::Result<()> {
    let mut f = std::fs::File::create(path)?;
    f.write_all(write_obj(mesh).as_bytes())?;
    Ok(())
}

/// Read OBJ from file.
pub fn read_obj_file(path: &str) -> Result<ObjMesh, String> {
    let content = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    ObjMesh::parse(&content)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_triangle() {
        let obj = "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n";
        let mesh = ObjMesh::parse(obj).unwrap();
        assert_eq!(mesh.vertices.len(), 3);
        assert_eq!(mesh.faces.len(), 1);
        assert_eq!(mesh.faces[0].v, vec![0, 1, 2]);
    }

    #[test]
    fn parse_negative_index() {
        let obj = "v 0 0 0\nv 1 0 0\nv 0 1 0\nf -3 -2 -1\n";
        let mesh = ObjMesh::parse(obj).unwrap();
        assert_eq!(mesh.faces[0].v, vec![0, 1, 2]);
    }

    #[test]
    fn roundtrip() {
        let mesh = ObjMesh {
            vertices: vec![GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.), GpPnt::new(0.,1.,0.)],
            texcoords: vec![], normals: vec![],
            faces: vec![ObjFace { v: vec![0,1,2], vt: None, vn: None }],
        };
        let out = write_obj(&mesh);
        let back = ObjMesh::parse(&out).unwrap();
        assert_eq!(back.vertices.len(), 3);
        assert_eq!(back.faces[0].v, vec![0, 1, 2]);
    }

    #[test]
    fn triangulate_quad() {
        let obj = "v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nf 1 2 3 4\n";
        let mesh = ObjMesh::parse(obj).unwrap();
        let tri = mesh.to_triangulation();
        assert_eq!(tri.triangles.len(), 2);
    }
}
