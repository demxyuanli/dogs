//! STL (STereoLithography) mesh reader/writer. Source: `RWStl` (simplified).
//! Supports both ASCII and little-endian binary STL.

use std::collections::HashMap;
use std::fs;

use crate::bnd::BndBox;
use crate::gp::GpPnt;
use crate::poly::triangulation::Triangle;
use crate::poly::Triangulation;

/// Triangle mesh with one normal per facet.
#[derive(Debug, Clone, Default)]
pub struct StlMesh {
    pub triangles: Vec<[GpPnt; 3]>,
    pub normals: Vec<[f64; 3]>,
}

/// Parse ASCII STL text (`facet normal nx ny nz` / `outer loop` / `vertex x y z` x3 / `endloop` / `endfacet`).
pub fn parse_ascii_stl(content: &str) -> Result<StlMesh, String> {
    let mut mesh = StlMesh::default();
    let mut normal = [0.0f64; 3];
    let mut verts: Vec<GpPnt> = Vec::new();
    for line in content.lines() {
        let mut parts = line.trim().split_whitespace();
        match parts.next() {
            Some("facet") => {
                normal = [0.0; 3];
                verts.clear();
                if parts.next() == Some("normal") {
                    let nx: f64 = parts.next().ok_or("facet normal missing x")?.parse().map_err(|e| format!("bad normal x: {e}"))?;
                    let ny: f64 = parts.next().ok_or("facet normal missing y")?.parse().map_err(|e| format!("bad normal y: {e}"))?;
                    let nz: f64 = parts.next().ok_or("facet normal missing z")?.parse().map_err(|e| format!("bad normal z: {e}"))?;
                    normal = [nx, ny, nz];
                }
            }
            Some("vertex") => {
                let x: f64 = parts.next().ok_or("vertex missing x")?.parse().map_err(|e| format!("bad vertex x: {e}"))?;
                let y: f64 = parts.next().ok_or("vertex missing y")?.parse().map_err(|e| format!("bad vertex y: {e}"))?;
                let z: f64 = parts.next().ok_or("vertex missing z")?.parse().map_err(|e| format!("bad vertex z: {e}"))?;
                verts.push(GpPnt::new(x, y, z));
            }
            Some("endfacet") => {
                if verts.len() != 3 {
                    return Err(format!("facet has {} vertices, expected 3", verts.len()));
                }
                mesh.triangles.push([verts[0], verts[1], verts[2]]);
                mesh.normals.push(normal);
                verts.clear();
            }
            // solid / endsolid / outer loop / endloop / blank lines: ignored
            _ => {}
        }
    }
    Ok(mesh)
}

/// Parse binary STL: 80-byte header, u32 LE triangle count, then 50-byte records
/// (12 f32 LE: normal + 3 vertices, then u16 LE attribute word).
pub fn parse_binary_stl(bytes: &[u8]) -> Result<StlMesh, String> {
    if bytes.len() < 84 {
        return Err(format!("binary STL too short: {} bytes", bytes.len()));
    }
    let count = u32::from_le_bytes([bytes[80], bytes[81], bytes[82], bytes[83]]) as usize;
    let need = 84usize
        .checked_add(count.checked_mul(50).ok_or("triangle count overflow")?)
        .ok_or("size overflow")?;
    if bytes.len() < need {
        return Err(format!("binary STL truncated: expected {need} bytes, got {}", bytes.len()));
    }
    let mut mesh = StlMesh::default();
    for i in 0..count {
        let rec = &bytes[84 + i * 50..84 + i * 50 + 50];
        let mut f = [0.0f32; 12];
        for k in 0..12 {
            f[k] = f32::from_le_bytes([rec[4 * k], rec[4 * k + 1], rec[4 * k + 2], rec[4 * k + 3]]);
        }
        mesh.normals.push([f[0] as f64, f[1] as f64, f[2] as f64]);
        mesh.triangles.push([
            GpPnt::new(f[3] as f64, f[4] as f64, f[5] as f64),
            GpPnt::new(f[6] as f64, f[7] as f64, f[8] as f64),
            GpPnt::new(f[9] as f64, f[10] as f64, f[11] as f64),
        ]);
    }
    Ok(mesh)
}

/// Detect STL format: the first five bytes are `solid` ⇒ ASCII, else binary.
///
/// `RWStl.cxx:523-531`: OCCT reads 5 bytes and sets
/// `isAscii = (strncmp(aHeader, "solid", 5) == 0)` — it does **not** additionally
/// look for `facet` (the previous body required both, an invented extra rule,
/// audit A26).
pub fn detect_stl_format(bytes: &[u8]) -> &'static str {
    let head = &bytes[..bytes.len().min(5)];
    if head == b"solid" {
        "ascii"
    } else {
        "binary"
    }
}

/// Read an STL file, auto-detecting ASCII vs binary.
pub fn read_stl_file(path: &str) -> Result<StlMesh, String> {
    let bytes = fs::read(path).map_err(|e| format!("cannot read {}: {e}", path))?;
    match detect_stl_format(&bytes) {
        "ascii" => parse_ascii_stl(&String::from_utf8_lossy(&bytes)),
        _ => parse_binary_stl(&bytes),
    }
}

/// Serialize mesh to little-endian binary STL.
pub fn write_binary_stl(mesh: &StlMesh) -> Vec<u8> {
    let mut out = Vec::with_capacity(84 + mesh.triangles.len() * 50);
    let mut header = [0u8; 80];
    // `RWStl.cxx:374-375`: `char aHeader[80] = "STL Exported by Open CASCADE
    // Technology [dev.opencascade.org]"; theStream.write(aHeader, 80);` — the
    // 80-byte header is that literal, zero-padded (the previous body wrote a
    // port-specific "solid generated by occt-core", audit A26).
    const OCCT_HEADER: &[u8] = b"STL Exported by Open CASCADE Technology [dev.opencascade.org]";
    let n = OCCT_HEADER.len().min(80);
    header[..n].copy_from_slice(&OCCT_HEADER[..n]);
    out.extend_from_slice(&header);
    out.extend_from_slice(&(mesh.triangles.len() as u32).to_le_bytes());
    for (i, tri) in mesh.triangles.iter().enumerate() {
        let n = if i < mesh.normals.len() { mesh.normals[i] } else { compute_normal(tri) };
        let coords = [
            n[0] as f32, n[1] as f32, n[2] as f32,
            tri[0].x() as f32, tri[0].y() as f32, tri[0].z() as f32,
            tri[1].x() as f32, tri[1].y() as f32, tri[1].z() as f32,
            tri[2].x() as f32, tri[2].y() as f32, tri[2].z() as f32,
        ];
        for c in coords {
            out.extend_from_slice(&c.to_le_bytes());
        }
        out.extend_from_slice(&0u16.to_le_bytes());
    }
    out
}

/// Serialize mesh to ASCII STL.
pub fn write_ascii_stl(mesh: &StlMesh) -> String {
    let mut s = String::with_capacity(16 + mesh.triangles.len() * 96);
    // `RWStl.cxx:304-305`: `// note that space after 'solid' is necessary for
    // many systems` / `theStream << "solid \n";` (no solid name).
    s.push_str("solid \n");
    for (i, tri) in mesh.triangles.iter().enumerate() {
        let n = if i < mesh.normals.len() { mesh.normals[i] } else { compute_normal(tri) };
        s.push_str(&format!("facet normal {} {} {}\n", n[0], n[1], n[2]));
        s.push_str("  outer loop\n");
        for p in tri {
            s.push_str(&format!("    vertex {} {} {}\n", p.x(), p.y(), p.z()));
        }
        s.push_str("  endloop\n");
        s.push_str("endfacet\n");
    }
    s.push_str("endsolid\n");
    s
}

/// Build a `Triangulation`, merging facet nodes whose coordinates are **exactly**
/// equal - `RWStl_Reader` merges "on the fly" through `Poly_MergeNodesTool`
/// (`RWStl_Reader.hxx:36` "The nodes with equal coordinates are merged
/// automatically on the fly") whose merge tolerance defaults to `0.0`
/// (`Poly_MergeNodesTool.hxx:50-56`: "0.0 by default (only 3D points with exactly
/// matching coordinates are merged)"; the merge angle defaults to `M_PI/2`, i.e.
/// all nodes merge regardless of the facet angle, `RWStl_Reader.hxx:93-95`).
/// The previous body quantised every coordinate onto a `1e-9` grid, an invented
/// tolerance (audit A26) that also merged near-coincident nodes OCCT keeps apart.
pub fn to_triangulation(mesh: &StlMesh) -> Triangulation {
    let mut nodes: Vec<GpPnt> = Vec::new();
    let mut map: HashMap<(u64, u64, u64), usize> = HashMap::new();
    // `0.0 == -0.0` in OCCT's coordinate comparison, so normalise the sign bit
    // before hashing.
    let key_of = |v: f64| if v == 0.0 { 0.0f64.to_bits() } else { v.to_bits() };
    let mut tris = Vec::with_capacity(mesh.triangles.len());
    for tri in &mesh.triangles {
        let mut idx = [0usize; 3];
        for (j, p) in tri.iter().enumerate() {
            let key = (key_of(p.x()), key_of(p.y()), key_of(p.z()));
            idx[j] = *map.entry(key).or_insert_with(|| {
                nodes.push(*p);
                nodes.len() - 1
            });
        }
        tris.push(Triangle::new(idx[0], idx[1], idx[2]));
    }
    Triangulation::new(nodes, tris)
}

/// Axis-aligned bounding box of all mesh vertices.
pub fn mesh_bbox(mesh: &StlMesh) -> BndBox {
    let mut b = BndBox::new();
    for tri in &mesh.triangles {
        for p in tri {
            b.add_point(p);
        }
    }
    b
}

/// Area of a triangle via half the cross-product magnitude.
pub fn triangle_area(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> f64 {
    let ab = (b.x() - a.x(), b.y() - a.y(), b.z() - a.z());
    let ac = (c.x() - a.x(), c.y() - a.y(), c.z() - a.z());
    let cx = ab.1 * ac.2 - ab.2 * ac.1;
    let cy = ab.2 * ac.0 - ab.0 * ac.2;
    let cz = ab.0 * ac.1 - ab.1 * ac.0;
    0.5 * (cx * cx + cy * cy + cz * cz).sqrt()
}

/// Total surface area of all triangles.
pub fn total_area(mesh: &StlMesh) -> f64 {
    mesh.triangles.iter().map(|t| triangle_area(&t[0], &t[1], &t[2])).sum()
}

/// Recompute all facet normals (unit length, zero for degenerate facets).
pub fn compute_facet_normals(mesh: &mut StlMesh) {
    mesh.normals = mesh.triangles.iter().map(compute_normal).collect();
}

/// Unit normal of a triangle via normalized cross product; the zero vector when
/// the cross product's **squared** magnitude does not exceed `gp::Resolution()`
/// (`RWStl.cxx:325` / `:407`: `if (aVNorm.SquareMagnitude() > gp::Resolution())`,
/// and `gp::Resolution() == RealSmall() == DBL_MIN`). The previous threshold was
/// the invented `|cross| < 1e-12` (audit A26).
fn compute_normal(tri: &[GpPnt; 3]) -> [f64; 3] {
    let (a, b, c) = (&tri[0], &tri[1], &tri[2]);
    let ab = (b.x() - a.x(), b.y() - a.y(), b.z() - a.z());
    let ac = (c.x() - a.x(), c.y() - a.y(), c.z() - a.z());
    let cx = ab.1 * ac.2 - ab.2 * ac.1;
    let cy = ab.2 * ac.0 - ab.0 * ac.2;
    let cz = ab.0 * ac.1 - ab.1 * ac.0;
    let sq = cx * cx + cy * cy + cz * cz;
    if sq > crate::precision::REAL_SMALL {
        let len = sq.sqrt();
        [cx / len, cy / len, cz / len]
    } else {
        [0.0, 0.0, 0.0]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn single_triangle() -> StlMesh {
        let tri = [
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
        ];
        StlMesh { triangles: vec![tri], normals: vec![[0.0, 0.0, 1.0]] }
    }

    #[test]
    fn parse_ascii_single_triangle() {
        let content = "solid demo
facet normal 0 0 1
  outer loop
    vertex 0 0 0
    vertex 1 0 0
    vertex 0 1 0
  endloop
endfacet
endsolid demo
";
        let m = parse_ascii_stl(content).unwrap();
        assert_eq!(m.triangles.len(), 1);
        assert_eq!(m.normals[0], [0.0, 0.0, 1.0]);
        let t = &m.triangles[0];
        assert_eq!(t[0].x(), 0.0);
        assert_eq!(t[1].y(), 0.0);
        assert_eq!(t[2].z(), 0.0);
    }

    #[test]
    fn binary_write_read_roundtrip() {
        let m = single_triangle();
        let bytes = write_binary_stl(&m);
        assert_eq!(bytes.len(), 84 + 50);
        let back = parse_binary_stl(&bytes).unwrap();
        assert_eq!(back.triangles.len(), 1);
        assert_eq!(back.triangles[0][1].x(), 1.0);
        assert_eq!(back.triangles[0][2].y(), 1.0);
        assert_eq!(back.normals[0], [0.0, 0.0, 1.0]);
    }

    #[test]
    fn to_triangulation_one_triangle() {
        let t = to_triangulation(&single_triangle());
        assert_eq!(t.nodes.len(), 3);
        assert_eq!(t.triangles.len(), 1);
        assert_eq!(t.triangles[0].n2, 2);
    }

    #[test]
    fn triangle_area_half() {
        let a = GpPnt::new(0.0, 0.0, 0.0);
        let b = GpPnt::new(1.0, 0.0, 0.0);
        let c = GpPnt::new(0.0, 1.0, 0.0);
        assert!((triangle_area(&a, &b, &c) - 0.5).abs() < 1e-12);
        assert!((total_area(&single_triangle()) - 0.5).abs() < 1e-12);
    }
}
