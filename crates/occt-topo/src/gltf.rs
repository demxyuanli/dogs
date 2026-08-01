//! glTF 2.0 mesh export (`RWMesh` / `XCAFGlTF`-lite).
//!
//! Meshes a `TopoShape` with `mesh_shape` and writes a valid glTF 2.0 `.gltf`
//! JSON document plus a `.bin` buffer holding `f32` positions, optional `f32`
//! area-weighted vertex normals, and `u32` triangle indices. The buffer can be
//! referenced externally (`mesh.bin`) or embedded as a base64 data URI.
//!
//! The JSON is assembled by hand (no external JSON dependency): a single
//! `buffer`, three `bufferView`s (positions / normals / indices), three
//! `accessor`s, one mesh primitive (`mode = 4` triangles), one node and one
//! scene. Vertex normals are computed as the area-weighted average of incident
//! triangle face normals; degenerate averages fall back to `(0, 0, 1)`.
//! Non-finite vertices are dropped so NaN/Inf can never leak into the buffer.

use occt_core::gp::{GpPnt, GpXyz};
use occt_core::poly::triangulation::Triangle;

use crate::mesh::ShapeMesh;
use crate::shape::TopoShape;

/// Export options for the glTF writer.
pub struct GltfOptions {
    /// Emit a `NORMAL` accessor with area-weighted vertex normals.
    pub include_normals: bool,
    /// Embed the buffer as a base64 data URI instead of an external `mesh.bin`.
    pub embed_buffer: bool,
    /// Pretty-print the JSON (multi-line, 2-space indent) vs compact.
    pub pretty: bool,
}

impl Default for GltfOptions {
    fn default() -> Self {
        Self {
            include_normals: true,
            embed_buffer: false,
            pretty: true,
        }
    }
}

/// Build the `.gltf` JSON text for `shape` (buffer bytes are discarded).
pub fn write_gltf(shape: &TopoShape, deflection: f64, options: &GltfOptions) -> Result<String, String> {
    Ok(build_gltf_bin(shape, deflection, options)?.0)
}

/// Build the complete glTF export: the `.gltf` JSON text and the `.bin` bytes.
///
/// Bin layout is `[positions f32 xyz × N][normals f32 xyz × N][indices u32 × T]`
/// with matching `byteOffset`s in the JSON. The position accessor carries
/// `min`/`max` arrays computed from the exact `f32` values written to the bin.
pub fn build_gltf_bin(
    shape: &TopoShape,
    deflection: f64,
    options: &GltfOptions,
) -> Result<(String, Vec<u8>), String> {
    let mesh = crate::shape_mesh::mesh_shape(shape, deflection);
    let (verts, tris) = sanitize_mesh(&mesh);
    if verts.is_empty() || tris.is_empty() {
        return Err("gltf: shape produced an empty mesh".into());
    }
    if verts.len() > u32::MAX as usize {
        return Err(format!("gltf: too many vertices ({}) for u32 indices", verts.len()));
    }
    let n = verts.len();
    let t = tris.len();

    let mut bin = Vec::with_capacity(24 * n + 12 * t);
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for p in &verts {
        let (fx, fy, fz) = (p.x() as f32, p.y() as f32, p.z() as f32);
        bin.extend_from_slice(&fx.to_le_bytes());
        bin.extend_from_slice(&fy.to_le_bytes());
        bin.extend_from_slice(&fz.to_le_bytes());
        let (x, y, z) = (fx as f64, fy as f64, fz as f64);
        min[0] = min[0].min(x);
        min[1] = min[1].min(y);
        min[2] = min[2].min(z);
        max[0] = max[0].max(x);
        max[1] = max[1].max(y);
        max[2] = max[2].max(z);
    }

    let has_normals = options.include_normals;
    if has_normals {
        for nx in &compute_normals(&verts, &tris) {
            bin.extend_from_slice(&nx[0].to_le_bytes());
            bin.extend_from_slice(&nx[1].to_le_bytes());
            bin.extend_from_slice(&nx[2].to_le_bytes());
        }
    }
    for tr in &tris {
        bin.extend_from_slice(&(tr.n0 as u32).to_le_bytes());
        bin.extend_from_slice(&(tr.n1 as u32).to_le_bytes());
        bin.extend_from_slice(&(tr.n2 as u32).to_le_bytes());
    }

    let uri = if options.embed_buffer {
        format!("data:application/octet-stream;base64,{}", base64_encode(&bin))
    } else {
        "mesh.bin".to_string()
    };
    let json = build_gltf_json(n, t, bin.len(), &uri, &min, &max, has_normals, options);
    Ok((json, bin))
}

/// Write `<path>.gltf` and, unless the buffer is embedded, `<path>.bin`.
pub fn write_gltf_file(
    shape: &TopoShape,
    deflection: f64,
    path: &str,
    options: &GltfOptions,
) -> Result<(), String> {
    let (gltf, bin) = build_gltf_bin(shape, deflection, options)?;
    std::fs::write(path, &gltf).map_err(|e| format!("gltf: write {path}: {e}"))?;
    if !options.embed_buffer {
        let bin_path = path.strip_suffix(".gltf").unwrap_or(path).to_string() + ".bin";
        std::fs::write(&bin_path, &bin).map_err(|e| format!("gltf: write {bin_path}: {e}"))?;
    }
    Ok(())
}

/// Standard base64 encoding (RFC 4648) without external dependencies.
pub fn base64_encode(bytes: &[u8]) -> String {
    const ALPHABET: &[u8; 64] =
        b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::with_capacity((bytes.len() + 2) / 3 * 4);
    let mut i = 0;
    while i + 3 <= bytes.len() {
        let b = [bytes[i], bytes[i + 1], bytes[i + 2]];
        out.push(ALPHABET[(b[0] >> 2) as usize] as char);
        out.push(ALPHABET[(((b[0] & 0x03) << 4) | (b[1] >> 4)) as usize] as char);
        out.push(ALPHABET[(((b[1] & 0x0F) << 2) | (b[2] >> 6)) as usize] as char);
        out.push(ALPHABET[(b[2] & 0x3F) as usize] as char);
        i += 3;
    }
    let rem = bytes.len() - i;
    if rem == 1 {
        let b = bytes[i];
        out.push(ALPHABET[(b >> 2) as usize] as char);
        out.push(ALPHABET[((b & 0x03) << 4) as usize] as char);
        out.push('=');
        out.push('=');
    } else if rem == 2 {
        let b = [bytes[i], bytes[i + 1]];
        out.push(ALPHABET[(b[0] >> 2) as usize] as char);
        out.push(ALPHABET[(((b[0] & 0x03) << 4) | (b[1] >> 4)) as usize] as char);
        out.push(ALPHABET[((b[1] & 0x0F) << 2) as usize] as char);
        out.push('=');
    }
    out
}

// ---------------------------------------------------------------------------
// Mesh helpers
// ---------------------------------------------------------------------------

/// Drop non-finite vertices and any triangle that references a dropped vertex,
/// remapping indices so the exported buffer contains only finite floats.
fn sanitize_mesh(mesh: &ShapeMesh) -> (Vec<GpPnt>, Vec<Triangle>) {
    let mut remap = vec![usize::MAX; mesh.vertices.len()];
    let mut verts = Vec::new();
    for (i, p) in mesh.vertices.iter().enumerate() {
        if p.x().is_finite() && p.y().is_finite() && p.z().is_finite() {
            remap[i] = verts.len();
            verts.push(*p);
        }
    }
    let mut tris = Vec::new();
    for t in &mesh.triangles {
        let (a, b, c) = (remap[t.n0], remap[t.n1], remap[t.n2]);
        if a != usize::MAX && b != usize::MAX && c != usize::MAX {
            tris.push(Triangle::new(a, b, c));
        }
    }
    (verts, tris)
}

/// Area-weighted vertex normals. Zero-length accumulations fall back to
/// `(0, 0, 1)` so every emitted normal is unit-length.
fn compute_normals(verts: &[GpPnt], tris: &[Triangle]) -> Vec<[f32; 3]> {
    let n = verts.len();
    let mut accum = vec![GpXyz::zero(); n];
    for t in tris {
        let a = verts[t.n0].coord;
        let b = verts[t.n1].coord;
        let c = verts[t.n2].coord;
        let nrm = b.subtracted(&a).crossed(&c.subtracted(&a));
        accum[t.n0] = accum[t.n0].added(&nrm);
        accum[t.n1] = accum[t.n1].added(&nrm);
        accum[t.n2] = accum[t.n2].added(&nrm);
    }
    accum
        .into_iter()
        .map(|v| {
            let m = v.modulus();
            if m > 1e-30 {
                let u = v.divided(m);
                [u.x as f32, u.y as f32, u.z as f32]
            } else {
                [0.0, 0.0, 1.0]
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// JSON assembly
// ---------------------------------------------------------------------------

/// Render an f64 as a valid JSON number (always carries a `.0` when integral).
fn json_num(x: f64) -> String {
    if !x.is_finite() {
        return "0.0".into();
    }
    let s = format!("{:?}", x);
    if s.contains('.') || s.contains('e') || s.contains('E') {
        s
    } else {
        format!("{s}.0")
    }
}

fn json_arr3(v: &[f64; 3]) -> String {
    format!("[{},{},{}]", json_num(v[0]), json_num(v[1]), json_num(v[2]))
}

/// Escape a string for a JSON string literal (quotes, backslash, control
/// characters). ASCII digits/letters pass through untouched, so the generator
/// name, `mesh.bin` uri and base64 data uris are emitted byte-for-byte.
fn json_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out
}

/// Remove all whitespace outside JSON string literals (safe here: the emitted
/// strings — generator, uri, type names — contain no whitespace).
fn minify_json(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut in_str = false;
    let mut chars = s.chars();
    while let Some(c) = chars.next() {
        if in_str {
            out.push(c);
            if c == '\\' {
                if let Some(n) = chars.next() {
                    out.push(n);
                }
            } else if c == '"' {
                in_str = false;
            }
        } else {
            match c {
                '"' => {
                    in_str = true;
                    out.push(c);
                }
                c if c.is_whitespace() => {}
                c => out.push(c),
            }
        }
    }
    out
}

/// Build the glTF 2.0 JSON document (pretty form; compacted when `!pretty`).
#[allow(clippy::too_many_arguments)]
fn build_gltf_json(
    n: usize,
    t: usize,
    bin_len: usize,
    uri: &str,
    min: &[f64; 3],
    max: &[f64; 3],
    has_normals: bool,
    opts: &GltfOptions,
) -> String {
    let pos_bytes = 12 * n;
    let norm_bytes = if has_normals { 12 * n } else { 0 };
    let idx_off = pos_bytes + norm_bytes;
    let idx_bytes = 12 * t;
    let idx_accessor = if has_normals { 2 } else { 1 };
    let min_s = json_arr3(min);
    let max_s = json_arr3(max);

    let mut s = String::new();
    s.push_str("{\n");
    s.push_str("  \"asset\": {\n");
    s.push_str("    \"version\": \"2.0\",\n");
    s.push_str(&format!(
        "    \"generator\": \"{}\",\n",
        json_escape("occt-rust")
    ));
    s.push_str("  },\n");
    s.push_str("  \"scene\": 0,\n");
    s.push_str("  \"scenes\": [\n");
    s.push_str("    {\"nodes\": [0]}\n");
    s.push_str("  ],\n");
    s.push_str("  \"nodes\": [\n");
    s.push_str("    {\"mesh\": 0}\n");
    s.push_str("  ],\n");
    s.push_str("  \"meshes\": [\n");
    let attrs = if has_normals {
        "\"POSITION\": 0, \"NORMAL\": 1"
    } else {
        "\"POSITION\": 0"
    };
    s.push_str(&format!(
        "    {{\"primitives\": [{{\"attributes\": {{{attrs}}}, \"indices\": {idx_accessor}, \"mode\": 4}}]}}\n"
    ));
    s.push_str("  ],\n");
    s.push_str("  \"buffers\": [\n");
    s.push_str(&format!(
        "    {{\"uri\": \"{}\", \"byteLength\": {bin_len}}}\n",
        json_escape(uri)
    ));
    s.push_str("  ],\n");
    s.push_str("  \"bufferViews\": [\n");
    s.push_str(&format!(
        "    {{\"buffer\": 0, \"byteOffset\": 0, \"byteLength\": {pos_bytes}, \"target\": 34962}},\n"
    ));
    if has_normals {
        s.push_str(&format!(
            "    {{\"buffer\": 0, \"byteOffset\": {pos_bytes}, \"byteLength\": {norm_bytes}, \"target\": 34962}},\n"
        ));
    }
    s.push_str(&format!(
        "    {{\"buffer\": 0, \"byteOffset\": {idx_off}, \"byteLength\": {idx_bytes}, \"target\": 34963}}\n"
    ));
    s.push_str("  ],\n");
    s.push_str("  \"accessors\": [\n");
    s.push_str(&format!(
        "    {{\"bufferView\": 0, \"componentType\": 5126, \"count\": {n}, \"type\": \"VEC3\", \"min\": {min_s}, \"max\": {max_s}}},\n"
    ));
    if has_normals {
        s.push_str(&format!(
            "    {{\"bufferView\": 1, \"componentType\": 5126, \"count\": {n}, \"type\": \"VEC3\"}},\n"
        ));
    }
    s.push_str(&format!(
        "    {{\"bufferView\": {idx_accessor}, \"componentType\": 5125, \"count\": {t}, \"type\": \"SCALAR\"}}\n"
    ));
    s.push_str("  ]\n");
    s.push_str("}\n");

    if opts.pretty {
        s
    } else {
        minify_json(&s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::BRepPrimBox;

    /// Balanced-brace/string sanity check for generated JSON.
    fn balanced_json(s: &str) -> bool {
        let mut depth = 0i32;
        let mut in_str = false;
        let mut chars = s.chars();
        while let Some(c) = chars.next() {
            if in_str {
                if c == '\\' {
                    chars.next();
                } else if c == '"' {
                    in_str = false;
                }
            } else {
                match c {
                    '"' => in_str = true,
                    '{' => depth += 1,
                    '}' => {
                        depth -= 1;
                        if depth < 0 {
                            return false;
                        }
                    }
                    _ => {}
                }
            }
        }
        depth == 0 && !in_str
    }

    /// All `u64` values following `"key":` (whitespace-tolerant).
    fn json_nums(s: &str, key: &str) -> Vec<u64> {
        let needle = format!("\"{key}\":");
        let mut rest = s;
        let mut out = Vec::new();
        while let Some(pos) = rest.find(&needle) {
            let tail = &rest[pos + needle.len()..];
            let num: String = tail
                .chars()
                .skip_while(|c| c.is_whitespace())
                .take_while(|c| c.is_ascii_digit())
                .collect();
            if let Ok(v) = num.parse() {
                out.push(v);
            }
            rest = &tail[num.len()..];
        }
        out
    }

    /// The first `[a,b,c]` array following `"key":`.
    fn json_arr3_field(s: &str, key: &str) -> Option<[f64; 3]> {
        let needle = format!("\"{key}\":");
        let pos = s.find(&needle)?;
        let tail = &s[pos + needle.len()..];
        let tail = tail.trim_start_matches(|c: char| c.is_whitespace());
        let start = tail.find('[')?;
        let rest = &tail[start..];
        let end = rest.find(']')?;
        let parts: Vec<f64> = rest[1..end]
            .split(',')
            .filter_map(|x| x.trim().parse().ok())
            .collect();
        if parts.len() == 3 {
            Some([parts[0], parts[1], parts[2]])
        } else {
            None
        }
    }

    fn box_shape() -> TopoShape {
        BRepPrimBox::make_box(2.0, 2.0, 2.0).solid.0
    }

    fn compact(s: &str) -> String {
        s.chars().filter(|c| !c.is_whitespace()).collect()
    }

    #[test]
    fn gltf_box_valid_json() {
        let json = write_gltf(&box_shape(), 0.5, &GltfOptions::default()).expect("write");
        assert!(json.trim_start().starts_with('{'));
        for needle in ["\"asset\"", "\"version\": \"2.0\"", "\"meshes\"", "\"scenes\""] {
            assert!(json.contains(needle), "missing {needle} in:\n{json}");
        }
        assert!(balanced_json(&json), "unbalanced braces");
    }

    #[test]
    fn gltf_accessor_counts() {
        let shape = box_shape();
        let opts = GltfOptions::default();
        let (json, bin) = build_gltf_bin(&shape, 0.5, &opts).expect("build");
        let mesh = crate::shape_mesh::mesh_shape(&shape, 0.5);
        let n = mesh.vertices.len();
        let t = mesh.triangles.len();
        assert_eq!(bin.len(), 12 * n + 12 * n + 12 * t, "bin length");
        let c = compact(&json);
        assert!(c.contains(&format!("\"byteLength\":{bin_len}", bin_len = bin.len())));
        assert!(c.contains(&format!("\"count\":{n}")), "positions count missing");
    }

    #[test]
    fn gltf_min_max_present() {
        let json = write_gltf(&box_shape(), 0.5, &GltfOptions::default()).expect("write");
        assert!(json.contains("\"min\":"));
        assert!(json.contains("\"max\":"));
    }

    #[test]
    fn gltf_embed_buffer_uri() {
        let opts = GltfOptions {
            embed_buffer: true,
            ..Default::default()
        };
        let json = write_gltf(&box_shape(), 0.5, &opts).expect("write");
        assert!(
            json.contains("data:application/octet-stream;base64,"),
            "missing data uri:\n{json}"
        );
        assert!(!json.contains("mesh.bin"));
    }

    #[test]
    fn gltf_bin_decodable() {
        let shape = box_shape();
        let (json, bin) = build_gltf_bin(&shape, 0.5, &GltfOptions::default()).expect("build");
        let offsets = json_nums(&json, "byteOffset");
        let lengths = json_nums(&json, "byteLength");
        assert!(!offsets.is_empty() && !lengths.is_empty());
        let total = lengths[0];
        assert_eq!(total, bin.len() as u64, "buffer byteLength must match bin");
        // Each bufferView: own offset + own length <= total.
        for (&off, &len) in offsets.iter().zip(lengths.iter().skip(1)) {
            assert!(off + len <= total, "view off {off} + len {len} > total {total}");
        }
    }

    #[test]
    fn gltf_normals_unit_length() {
        let shape = box_shape();
        let opts = GltfOptions::default();
        let (_, bin) = build_gltf_bin(&shape, 0.5, &opts).expect("build");
        let n = crate::shape_mesh::mesh_shape(&shape, 0.5).vertices.len();
        let mut k = 12 * n;
        while k + 12 <= 24 * n {
            let x = f32::from_le_bytes([bin[k], bin[k + 1], bin[k + 2], bin[k + 3]]) as f64;
            let y = f32::from_le_bytes([bin[k + 4], bin[k + 5], bin[k + 6], bin[k + 7]]) as f64;
            let z = f32::from_le_bytes([bin[k + 8], bin[k + 9], bin[k + 10], bin[k + 11]]) as f64;
            let len2 = x * x + y * y + z * z;
            assert!((len2 - 1.0).abs() < 1e-3, "normal len {len2} at vertex {}", k / 12);
            k += 12;
        }
    }

    #[test]
    fn gltf_box_bbox() {
        let json = write_gltf(&box_shape(), 0.5, &GltfOptions::default()).expect("write");
        let min = json_arr3_field(&json, "min").expect("min array");
        let max = json_arr3_field(&json, "max").expect("max array");
        for i in 0..3 {
            assert!((min[i] - 0.0).abs() < 1e-6, "min[{i}] = {}", min[i]);
            assert!((max[i] - 2.0).abs() < 1e-6, "max[{i}] = {}", max[i]);
        }
    }

    #[test]
    fn gltf_no_nan() {
        let (_, bin) = build_gltf_bin(&box_shape(), 0.5, &GltfOptions::default()).expect("build");
        let mut k = 0;
        while k + 4 <= bin.len() {
            let v = f32::from_le_bytes([bin[k], bin[k + 1], bin[k + 2], bin[k + 3]]);
            assert!(v.is_finite(), "non-finite f32 at byte {k}: {v}");
            k += 4;
        }
    }

    #[test]
    fn gltf_compact_pretty_equivalent() {
        let pretty = GltfOptions::default();
        let compact_opts = GltfOptions {
            pretty: false,
            ..Default::default()
        };
        let p = write_gltf(&box_shape(), 0.5, &pretty).expect("pretty");
        let c = write_gltf(&box_shape(), 0.5, &compact_opts).expect("compact");
        assert_eq!(compact(&p), c, "minified pretty JSON must equal compact JSON");
    }

    #[test]
    fn base64_roundtrip_known_vector() {
        assert_eq!(base64_encode(b"f"), "Zg==");
        assert_eq!(base64_encode(b"fo"), "Zm8=");
        assert_eq!(base64_encode(b"foo"), "Zm9v");
        assert_eq!(base64_encode(b"foobar"), "Zm9vYmFy");
    }
}
