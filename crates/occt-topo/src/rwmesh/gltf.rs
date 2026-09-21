use super::prelude::*;
use super::*;


/// Column-major glTF 4x4 matrix → `GpTrsf` (linear part + translation).
pub(super) fn trsf_from_matrix(m: &[f64]) -> GpTrsf {
    let mut mat = occt_core::gp::GpMat::zero();
    mat.m[0][0] = m[0];
    mat.m[0][1] = m[4];
    mat.m[0][2] = m[8];
    mat.m[1][0] = m[1];
    mat.m[1][1] = m[5];
    mat.m[1][2] = m[9];
    mat.m[2][0] = m[2];
    mat.m[2][1] = m[6];
    mat.m[2][2] = m[10];
    GpTrsf {
        scale: 1.0,
        shape: TrsfForm::CompoundTrsf,
        matrix: mat,
        loc: GpXyz::new(m[12], m[13], m[14]),
    }
}

/// Depth-first walk of the node graph; emits one `MeshNode` per mesh node with
/// the accumulated world transform applied to its vertices and the primitive's
/// material (and base-color texture path, if any) attached.
pub(super) fn collect_gltf_nodes(
    nodes: &[GlNode],
    idx: usize,
    world: &GpTrsf,
    ctx: &GltfCtx,
    out: &mut Vec<MeshNode>,
) -> Result<(), String> {
    let n = &nodes[idx];
    let world = world.multiplied(&n.trsf);
    if let Some(m) = n.mesh {
        let mesh = ctx.meshes.get(m).ok_or_else(|| format!("gltf: node {idx}: mesh {m} out of range"))?;
        for prim in &mesh.prims {
            let (verts, tris) = read_primitive(prim, ctx.buffers, ctx.views)?;
            let verts: Vec<GpPnt> = verts.iter().map(|p| p.transformed(&world)).collect();
            let name = if n.name.is_empty() { format!("mesh_{idx}") } else { n.name.clone() };

            let (material, texture_path) = prim.mat_idx
                .and_then(|mi| ctx.materials.get(mi))
                .map(|m| {
                    let mat = Material {
                        name: m.name.clone(),
                        diffuse: m.diffuse,
                        specular: (1.0, 1.0, 1.0),
                        opacity: m.opacity,
                    };
                    let tex = m.base_color_texture
                        .and_then(|ti| ctx.textures.get(ti).and_then(|o| *o))
                        .and_then(|ii| ctx.images.get(ii))
                        .map(|uri| resolve_uri(uri, ctx.base_dir));
                    (mat, tex)
                })
                .map(|(mat, tex)| (Some(mat), tex))
                .unwrap_or((None, None));

            out.push(MeshNode {
                name,
                vertices: verts,
                triangles: tris,
                transform: Some(world.clone()),
                material,
                uv: None,
                texture_path,
            });
        }
    }
    for &c in &n.children {
        collect_gltf_nodes(nodes, c, &world, ctx, out)?;
    }
    Ok(())
}

/// Resolve a glTF image URI against the scene directory (data: URIs pass through).
pub(super) fn resolve_uri(uri: &str, base_dir: &Path) -> String {
    if uri.starts_with("data:") {
        uri.to_string()
    } else {
        base_dir.join(uri).to_string_lossy().to_string()
    }
}

/// Read a mesh primitive's positions and indices from the buffer views.
pub(super) fn read_primitive(prim: &Prim, buffers: &[Vec<u8>], views: &[BufView]) -> Result<(Vec<GpPnt>, Vec<(usize, usize, usize)>), String> {
    let vbytes = view_bytes(prim.pos_view, buffers, views)?;
    if prim.pos_comp != 5126 {
        return Err(format!("gltf: POSITION componentType {} not supported", prim.pos_comp));
    }
    let mut verts = Vec::with_capacity(prim.pos_count);
    for i in 0..prim.pos_count {
        let off = prim.pos_offset + i * 12;
        if off + 12 > vbytes.len() {
            return Err(format!("gltf: POSITION accessor out of range (vertex {i})"));
        }
        let x = f32::from_le_bytes([vbytes[off], vbytes[off + 1], vbytes[off + 2], vbytes[off + 3]]) as f64;
        let y = f32::from_le_bytes([vbytes[off + 4], vbytes[off + 5], vbytes[off + 6], vbytes[off + 7]]) as f64;
        let z = f32::from_le_bytes([vbytes[off + 8], vbytes[off + 9], vbytes[off + 10], vbytes[off + 11]]) as f64;
        verts.push(GpPnt::new(x, y, z));
    }

    let mut tris = Vec::new();
    if let Some(vw) = prim.idx_view {
        let ibytes = view_bytes(Some(vw), buffers, views)?;
        let sz = match prim.idx_comp {
            5121 => 1,
            5123 => 2,
            _ => 4,
        };
        for i in 0..prim.idx_count / 3 {
            let o = prim.idx_offset + i * 3 * sz;
            if o + 3 * sz > ibytes.len() {
                return Err(format!("gltf: indices accessor out of range (tri {i})"));
            }
            let a = read_uint(prim.idx_comp, ibytes, o);
            let b = read_uint(prim.idx_comp, ibytes, o + sz);
            let c = read_uint(prim.idx_comp, ibytes, o + 2 * sz);
            tris.push((a, b, c));
        }
    } else {
        // No index accessor: vertices are a triangle list.
        for i in 0..prim.pos_count / 3 {
            tris.push((3 * i, 3 * i + 1, 3 * i + 2));
        }
    }
    Ok((verts, tris))
}

pub(super) fn view_bytes<'a>(view: Option<usize>, buffers: &'a [Vec<u8>], views: &[BufView]) -> Result<&'a [u8], String> {
    let v = views.get(view.ok_or("gltf: missing bufferView")?)
        .ok_or("gltf: bufferView index out of range")?;
    let buf = buffers.get(v.buffer).ok_or("gltf: buffer index out of range")?;
    if v.offset + v.len > buf.len() {
        return Err("gltf: bufferView extends past buffer end".into());
    }
    Ok(&buf[v.offset..v.offset + v.len])
}

pub(super) fn read_uint(comp: u32, b: &[u8], off: usize) -> usize {
    match comp {
        5121 => b[off] as usize,
        5123 => u16::from_le_bytes([b[off], b[off + 1]]) as usize,
        _ => u32::from_le_bytes([b[off], b[off + 1], b[off + 2], b[off + 3]]) as usize,
    }
}

/// Decode a standard base64 string (RFC 4648) with or without padding.
pub(super) fn base64_decode(s: &str) -> Result<Vec<u8>, String> {
    let s: String = s.chars().filter(|c| !c.is_whitespace()).collect();
    let s = s.trim_end_matches('=');
    let mut out = Vec::with_capacity(s.len() * 3 / 4);
    let mut bits: u32 = 0;
    let mut n: u32 = 0;
    for c in s.chars() {
        let v = match c {
            'A'..='Z' => c as u32 - 'A' as u32,
            'a'..='z' => c as u32 - 'a' as u32 + 26,
            '0'..='9' => c as u32 - '0' as u32 + 52,
            '+' => 62,
            '/' => 63,
            _ => return Err(format!("gltf: invalid base64 character '{c}'")),
        };
        bits = (bits << 6) | v;
        n += 6;
        if n >= 8 {
            n -= 8;
            out.push((bits >> n) as u8);
        }
    }
    Ok(out)
}

// ---------------------------------------------------------------------------
// Minimal JSON parser
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub(super) enum JVal {
    Null,
    Bool(bool),
    Num(f64),
    Str(String),
    Arr(Vec<JVal>),
    Obj(Vec<(String, JVal)>),
}

impl JVal {
    pub(super) fn as_arr(&self) -> Option<&[JVal]> {
        match self {
            JVal::Arr(a) => Some(a),
            _ => None,
        }
    }
    pub(super) fn as_obj(&self) -> Option<&[(String, JVal)]> {
        match self {
            JVal::Obj(o) => Some(o),
            _ => None,
        }
    }
    pub(super) fn as_str(&self) -> Option<&str> {
        match self {
            JVal::Str(s) => Some(s),
            _ => None,
        }
    }
    pub(super) fn as_num(&self) -> Option<f64> {
        match self {
            JVal::Num(n) => Some(*n),
            _ => None,
        }
    }
    pub(super) fn get(&self, key: &str) -> Option<&JVal> {
        self.as_obj().and_then(|o| o.iter().find(|(k, _)| k == key).map(|(_, v)| v))
    }
}

pub(super) struct JP<'a> {
    pub(super) s: &'a [u8],
    pub(super) i: usize,
}

impl<'a> JP<'a> {
    pub(super) fn ws(&mut self) {
        while self.i < self.s.len() && matches!(self.s[self.i], b' ' | b'\t' | b'\n' | b'\r') {
            self.i += 1;
        }
    }
    pub(super) fn peek(&self) -> Option<u8> {
        self.s.get(self.i).copied()
    }
    pub(super) fn eat(&mut self, c: u8) -> Result<(), String> {
        self.ws();
        if self.peek() == Some(c) {
            self.i += 1;
            Ok(())
        } else {
            Err(format!("json: expected '{}' at {}", c as char, self.i))
        }
    }
    pub(super) fn value(&mut self) -> Result<JVal, String> {
        self.ws();
        match self.peek() {
            Some(b'{') => self.object(),
            Some(b'[') => self.array(),
            Some(b'"') => Ok(JVal::Str(self.string()?)),
            Some(b't') => {
                self.lit("true")?;
                Ok(JVal::Bool(true))
            }
            Some(b'f') => {
                self.lit("false")?;
                Ok(JVal::Bool(false))
            }
            Some(b'n') => {
                self.lit("null")?;
                Ok(JVal::Null)
            }
            Some(c) if c == b'-' || c.is_ascii_digit() => Ok(JVal::Num(self.number()?)),
            _ => Err(format!("json: unexpected char at {}", self.i)),
        }
    }
    pub(super) fn lit(&mut self, word: &str) -> Result<(), String> {
        if self.s.get(self.i..).map(|x| x.starts_with(word.as_bytes())).unwrap_or(false) {
            self.i += word.len();
            Ok(())
        } else {
            Err(format!("json: expected '{word}'"))
        }
    }
    pub(super) fn string(&mut self) -> Result<String, String> {
        if self.peek() != Some(b'"') {
            return Err("json: expected string".into());
        }
        self.i += 1;
        let mut out = String::new();
        loop {
            let c = self.peek().ok_or("json: unterminated string")?;
            self.i += 1;
            match c {
                b'"' => break,
                b'\\' => {
                    let e = self.peek().ok_or("json: bad escape")?;
                    self.i += 1;
                    match e {
                        b'"' => out.push('"'),
                        b'\\' => out.push('\\'),
                        b'/' => out.push('/'),
                        b'b' => out.push('\u{8}'),
                        b'f' => out.push('\u{c}'),
                        b'n' => out.push('\n'),
                        b'r' => out.push('\r'),
                        b't' => out.push('\t'),
                        b'u' => {
                            let hex = self.s.get(self.i..self.i + 4).ok_or("json: bad \\u escape")?;
                            let n = u32::from_str_radix(std::str::from_utf8(hex).map_err(|_| "json: bad utf8")?, 16)
                                .map_err(|_| "json: bad \\u escape")?;
                            self.i += 4;
                            out.push(char::from_u32(n).ok_or("json: bad \\u value")?);
                        }
                        _ => return Err("json: unknown escape".into()),
                    }
                }
                0x00..=0x1F => return Err("json: control char in string".into()),
                c if c < 0x80 => out.push(c as char),
                _ => {
                    // multi-byte UTF-8
                    let extra = if c >= 0xF0 {
                        3
                    } else if c >= 0xE0 {
                        2
                    } else if c >= 0xC0 {
                        1
                    } else {
                        return Err("json: bad utf8".into());
                    };
                    let end = self.i + extra;
                    let seq = self.s.get(self.i..end).ok_or("json: bad utf8")?;
                    out.push_str(std::str::from_utf8(seq).map_err(|_| "json: bad utf8")?);
                    self.i = end;
                }
            }
        }
        Ok(out)
    }
    pub(super) fn number(&mut self) -> Result<f64, String> {
        let start = self.i;
        while self.i < self.s.len()
            && (self.s[self.i].is_ascii_digit()
                || matches!(self.s[self.i], b'-' | b'+' | b'.' | b'e' | b'E'))
        {
            self.i += 1;
        }
        std::str::from_utf8(&self.s[start..self.i])
            .map_err(|_| "json: bad number".to_string())?
            .parse::<f64>()
            .map_err(|_| "json: bad number".to_string())
    }
    pub(super) fn array(&mut self) -> Result<JVal, String> {
        self.i += 1; // '['
        let mut out = Vec::new();
        self.ws();
        if self.peek() == Some(b']') {
            self.i += 1;
            return Ok(JVal::Arr(out));
        }
        loop {
            out.push(self.value()?);
            self.ws();
            match self.peek() {
                Some(b',') => {
                    self.i += 1;
                }
                Some(b']') => {
                    self.i += 1;
                    break;
                }
                _ => return Err("json: expected ',' or ']'".into()),
            }
        }
        Ok(JVal::Arr(out))
    }
    pub(super) fn object(&mut self) -> Result<JVal, String> {
        self.i += 1; // '{'
        let mut out = Vec::new();
        self.ws();
        if self.peek() == Some(b'}') {
            self.i += 1;
            return Ok(JVal::Obj(out));
        }
        loop {
            self.ws();
            let key = self.string()?;
            self.eat(b':')?;
            let val = self.value()?;
            out.push((key, val));
            self.ws();
            match self.peek() {
                Some(b',') => {
                    self.i += 1;
                }
                Some(b'}') => {
                    self.i += 1;
                    break;
                }
                _ => return Err("json: expected ',' or '}'".into()),
            }
        }
        Ok(JVal::Obj(out))
    }
}

pub(super) fn parse_json(text: &str) -> Result<JVal, String> {
    let mut p = JP { s: text.as_bytes(), i: 0 };
    let v = p.value()?;
    p.ws();
    if p.i != p.s.len() {
        return Err("json: trailing content".into());
    }
    Ok(v)
}

pub(super) fn arr3(a: &[JVal]) -> [f64; 3] {
    let mut out = [0.0; 3];
    for (k, v) in a.iter().take(3).enumerate() {
        if let Some(n) = v.as_num() {
            out[k] = n;
        }
    }
    out
}

pub(super) fn arr4(a: &[JVal]) -> [f64; 4] {
    let mut out = [0.0; 4];
    for (k, v) in a.iter().take(4).enumerate() {
        if let Some(n) = v.as_num() {
            out[k] = n;
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Small utilities
// ---------------------------------------------------------------------------

pub(super) fn file_stem(path: &str) -> &str {
    Path::new(path).file_stem().and_then(|s| s.to_str()).unwrap_or("mesh")
}

pub(super) fn ext(path: &str) -> String {
    Path::new(path)
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_lowercase()
}
