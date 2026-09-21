use super::prelude::*;
use super::*;


/// Resolve a 1-based (or negative = from-end) OBJ index to a 0-based index.
pub(super) fn obj_index(i: i32, len: usize) -> Result<usize, String> {
    let r = if i > 0 {
        (i - 1) as usize
    } else if i < 0 {
        let v = len as i64 + i as i64;
        if v < 0 {
            return Err(format!("rwmesh: obj index {i} out of range (len {len})"));
        }
        v as usize
    } else {
        return Err("rwmesh: obj index 0 is invalid".into());
    };
    if r >= len {
        return Err(format!("rwmesh: obj index {i} out of range (len {len})"));
    }
    Ok(r)
}

// ---------------------------------------------------------------------------
// VRML 2.0 scene parser
// ---------------------------------------------------------------------------

/// Read a VRML 2.0 `.wrl` file into a scene. `Transform` nesting is accumulated
/// onto child vertices; `Shape`/`IndexedFaceSet` polygons are fan-triangulated;
/// a `Material` `diffuseColor` becomes the node material. A bare
/// `IndexedFaceSet` at the top level (no `Shape`/`Transform`) is also accepted.
pub fn read_vrml_scene(path: &str) -> Result<MeshScene, String> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| format!("rwmesh: read {path}: {e}"))?;
    parse_vrml_scene(&content, file_stem(path))
}

/// Parse VRML 2.0 text into a scene (see `read_vrml_scene`).
pub fn parse_vrml_scene(content: &str, name: &str) -> Result<MeshScene, String> {
    let toks = tokenize_vrml(content)?;
    let mut p = Vp { toks, i: 0 };
    let mut scene = MeshScene { name: name.to_string(), nodes: Vec::new(), materials: Vec::new() };
    while p.i < p.toks.len() {
        p.parse_scene_node(&GpTrsf::identity(), &mut scene)?;
    }
    Ok(scene)
}

#[derive(Debug, Clone, PartialEq)]
pub(super) enum VTok {
    Id(String),
    Num(f64),
    LBrace,
    RBrace,
    LBracket,
    RBracket,
}

/// Tokenize VRML 2.0: identifiers/numbers, `{}[]`, quoted strings become `Id`,
/// `#` starts a comment, commas and whitespace separate tokens.
pub(super) fn tokenize_vrml(content: &str) -> Result<Vec<VTok>, String> {
    let b = content.as_bytes();
    let mut toks = Vec::new();
    let mut i = 0;
    while i < b.len() {
        match b[i] {
            b' ' | b'\t' | b'\r' | b'\n' | b',' => i += 1,
            b'#' => {
                while i < b.len() && b[i] != b'\n' {
                    i += 1;
                }
            }
            b'{' => {
                toks.push(VTok::LBrace);
                i += 1;
            }
            b'}' => {
                toks.push(VTok::RBrace);
                i += 1;
            }
            b'[' => {
                toks.push(VTok::LBracket);
                i += 1;
            }
            b']' => {
                toks.push(VTok::RBracket);
                i += 1;
            }
            b'"' => {
                i += 1;
                let start = i;
                while i < b.len() && b[i] != b'"' {
                    i += 1;
                }
                if i >= b.len() {
                    return Err("vrml: unterminated string".into());
                }
                let s = std::str::from_utf8(&b[start..i]).map_err(|_| "vrml: bad string utf8")?;
                toks.push(VTok::Id(s.to_string()));
                i += 1;
            }
            c if c.is_ascii_alphabetic() || c == b'_' || c == b'.' || c == b'-' || c.is_ascii_digit() => {
                let start = i;
                while i < b.len()
                    && !matches!(b[i], b' ' | b'\t' | b'\r' | b'\n' | b',' | b'{' | b'}' | b'[' | b']' | b'"' | b'#')
                {
                    i += 1;
                }
                let w = std::str::from_utf8(&b[start..i]).map_err(|_| "vrml: bad token utf8")?;
                toks.push(match w.parse::<f64>() {
                    Ok(n) => VTok::Num(n),
                    Err(_) => VTok::Id(w.to_string()),
                });
            }
            c => return Err(format!("vrml: unexpected char '{}' at {i}", c as char)),
        }
    }
    Ok(toks)
}

/// Recursive-descent VRML parser. Understands `Transform`, `Shape`,
/// `Appearance`, `Material`, `IndexedFaceSet`, `Coordinate`,
/// `TextureCoordinate`; everything else is skipped structurally.
pub(super) struct Vp {
    pub(super) toks: Vec<VTok>,
    pub(super) i: usize,
}

impl Vp {
    pub(super) fn peek(&self) -> Option<&VTok> {
        self.toks.get(self.i)
    }
    pub(super) fn peek_is(&self, t: &VTok) -> bool {
        self.peek() == Some(t)
    }
    pub(super) fn peek_lbrace(&self) -> bool {
        self.peek_is(&VTok::LBrace)
    }
    pub(super) fn peek_rbrace(&self) -> bool {
        self.peek_is(&VTok::RBrace)
    }
    pub(super) fn peek_rbracket(&self) -> bool {
        self.peek_is(&VTok::RBracket)
    }
    pub(super) fn peek_id_eq(&self, s: &str) -> bool {
        matches!(self.peek(), Some(VTok::Id(x)) if x == s)
    }
    pub(super) fn take(&mut self) -> Option<VTok> {
        let t = self.toks.get(self.i).cloned();
        if t.is_some() {
            self.i += 1;
        }
        t
    }
    pub(super) fn expect_id(&mut self) -> Result<String, String> {
        match self.take() {
            Some(VTok::Id(s)) => Ok(s),
            other => Err(format!("vrml: expected identifier, got {other:?}")),
        }
    }
    pub(super) fn expect_num(&mut self) -> Result<f64, String> {
        match self.take() {
            Some(VTok::Num(n)) => Ok(n),
            other => Err(format!("vrml: expected number, got {other:?}")),
        }
    }
    pub(super) fn expect_lbrace(&mut self) -> Result<(), String> {
        if self.take() == Some(VTok::LBrace) {
            Ok(())
        } else {
            Err("vrml: expected '{{'".into())
        }
    }
    pub(super) fn expect_rbrace(&mut self) -> Result<(), String> {
        if self.take() == Some(VTok::RBrace) {
            Ok(())
        } else {
            Err("vrml: expected '}}'".into())
        }
    }
    pub(super) fn expect_lbracket(&mut self) -> Result<(), String> {
        if self.take() == Some(VTok::LBracket) {
            Ok(())
        } else {
            Err("vrml: expected '['".into())
        }
    }
    pub(super) fn expect_rbracket(&mut self) -> Result<(), String> {
        if self.take() == Some(VTok::RBracket) {
            Ok(())
        } else {
            Err("vrml: expected ']'".into())
        }
    }
    pub(super) fn parse_vec2(&mut self) -> Result<[f64; 2], String> {
        Ok([self.expect_num()?, self.expect_num()?])
    }
    pub(super) fn parse_vec3(&mut self) -> Result<[f64; 3], String> {
        Ok([self.expect_num()?, self.expect_num()?, self.expect_num()?])
    }
    pub(super) fn parse_vec4(&mut self) -> Result<[f64; 4], String> {
        Ok([self.expect_num()?, self.expect_num()?, self.expect_num()?, self.expect_num()?])
    }

    /// Skip one field value: `[ ... ]`, a `Node { ... }`, or a scalar run.
    pub(super) fn skip_value(&mut self) -> Result<(), String> {
        match self.peek().cloned() {
            Some(VTok::LBracket) => {
                self.take();
                self.skip_to_matching(&VTok::RBracket, &VTok::LBracket)?;
            }
            Some(VTok::LBrace) => {
                self.take();
                self.skip_to_matching(&VTok::RBrace, &VTok::LBrace)?;
            }
            Some(VTok::Id(_)) => {
                if self.i + 1 < self.toks.len() && self.toks[self.i + 1] == VTok::LBrace {
                    self.take();
                    self.take();
                    self.skip_to_matching(&VTok::RBrace, &VTok::LBrace)?;
                } else {
                    self.take();
                }
            }
            Some(VTok::Num(_)) => {
                while matches!(self.peek(), Some(VTok::Num(_))) {
                    self.take();
                }
            }
            Some(_) => {
                self.take();
            }
            None => return Err("vrml: unexpected end of file".into()),
        }
        Ok(())
    }

    pub(super) fn skip_to_matching(&mut self, close: &VTok, open: &VTok) -> Result<(), String> {
        let mut depth = 1;
        loop {
            match self.take() {
                Some(t) if &t == open => depth += 1,
                Some(t) if &t == close => {
                    depth -= 1;
                    if depth == 0 {
                        return Ok(());
                    }
                }
                Some(_) => {}
                None => return Err("vrml: unbalanced delimiters".into()),
            }
        }
    }

    /// Skip a whole unknown node `Name { ... }`.
    pub(super) fn skip_node(&mut self) -> Result<(), String> {
        self.expect_id()?;
        self.expect_lbrace()?;
        self.skip_to_matching(&VTok::RBrace, &VTok::LBrace)
    }

    pub(super) fn parse_scene_node(&mut self, world: &GpTrsf, scene: &mut MeshScene) -> Result<(), String> {
        match self.peek().cloned() {
            Some(VTok::Id(name)) if name == "Transform" => self.parse_transform(world, scene),
            Some(VTok::Id(name)) if name == "Shape" => self.parse_shape(world, scene),
            Some(VTok::Id(name)) if name == "IndexedFaceSet" => {
                // A bare IndexedFaceSet without a Shape/Transform wrapper.
                let (coord, coord_index, tex_coord) = self.parse_indexed_face_set()?;
                scene.nodes.push(build_vrml_node(&coord, &coord_index, tex_coord, world, None, "mesh"));
                Ok(())
            }
            Some(VTok::Id(_)) => self.skip_node(),
            other => Err(format!("vrml: expected a node, got {other:?}")),
        }
    }

    pub(super) fn parse_transform(&mut self, world: &GpTrsf, scene: &mut MeshScene) -> Result<(), String> {
        self.expect_id()?; // Transform
        self.expect_lbrace()?;
        let mut t = [0.0; 3];
        let mut axis = [0.0; 3];
        let mut angle = 0.0;
        let mut s = [1.0; 3];
        let mut children_at: Option<usize> = None;
        while !self.peek_rbrace() {
            let fname = self.expect_id()?;
            match fname.as_str() {
                "translation" => t = self.parse_vec3()?,
                "rotation" => {
                    let q = self.parse_vec4()?;
                    axis = [q[0], q[1], q[2]];
                    angle = q[3];
                }
                "scale" => s = self.parse_vec3()?,
                "children" => {
                    children_at = Some(self.i);
                    self.skip_value()?;
                }
                _ => self.skip_value()?,
            }
        }
        self.expect_rbrace()?;
        // Revisit `children` now that the local transform is fully known.
        if let Some(pos) = children_at {
            let saved = self.i;
            self.i = pos;
            let local = trsf_from_vrml(&t, &axis, angle, &s);
            let world2 = world.multiplied(&local);
            self.parse_children(&world2, scene)?;
            self.i = saved;
        }
        Ok(())
    }

    pub(super) fn parse_children(&mut self, world: &GpTrsf, scene: &mut MeshScene) -> Result<(), String> {
        if self.peek_is(&VTok::LBracket) {
            self.take();
            while !self.peek_rbracket() {
                self.parse_scene_node(world, scene)?;
            }
            self.take(); // ]
        } else {
            self.parse_scene_node(world, scene)?;
        }
        Ok(())
    }

    pub(super) fn parse_shape(&mut self, world: &GpTrsf, scene: &mut MeshScene) -> Result<(), String> {
        self.expect_id()?; // Shape
        self.expect_lbrace()?;
        let mut material: Option<Material> = None;
        let mut coord: Vec<GpPnt> = Vec::new();
        let mut coord_index: Vec<i32> = Vec::new();
        let mut tex_coord: Vec<GpPnt2d> = Vec::new();
        while !self.peek_rbrace() {
            let fname = self.expect_id()?;
            match fname.as_str() {
                "appearance" => {
                    if self.peek_id_eq("Appearance") {
                        material = self.parse_appearance()?;
                    } else {
                        self.skip_value()?;
                    }
                }
                "geometry" => {
                    if self.peek_id_eq("IndexedFaceSet") {
                        let (c, ci, tc) = self.parse_indexed_face_set()?;
                        coord = c;
                        coord_index = ci;
                        tex_coord = tc;
                    } else {
                        self.skip_value()?;
                    }
                }
                _ => self.skip_value()?,
            }
        }
        self.expect_rbrace()?;
        if !coord.is_empty() && !coord_index.is_empty() {
            scene.nodes.push(build_vrml_node(&coord, &coord_index, tex_coord, world, material, "mesh"));
        }
        Ok(())
    }

    pub(super) fn parse_appearance(&mut self) -> Result<Option<Material>, String> {
        self.expect_id()?; // Appearance
        self.expect_lbrace()?;
        let mut mat = None;
        while !self.peek_rbrace() {
            let fname = self.expect_id()?;
            match fname.as_str() {
                "material" => {
                    if self.peek_id_eq("Material") {
                        mat = Some(self.parse_material()?);
                    } else {
                        self.skip_value()?;
                    }
                }
                _ => self.skip_value()?,
            }
        }
        self.expect_rbrace()?;
        Ok(mat)
    }

    pub(super) fn parse_material(&mut self) -> Result<Material, String> {
        self.expect_id()?; // Material
        self.expect_lbrace()?;
        let mut m = Material::default();
        while !self.peek_rbrace() {
            let fname = self.expect_id()?;
            match fname.as_str() {
                "diffuseColor" => {
                    let c = self.parse_vec3()?;
                    m.diffuse = (c[0], c[1], c[2]);
                }
                "specularColor" => {
                    let c = self.parse_vec3()?;
                    m.specular = (c[0], c[1], c[2]);
                }
                "transparency" => {
                    let v = self.expect_num()?;
                    m.opacity = 1.0 - v.clamp(0.0, 1.0);
                }
                _ => self.skip_value()?,
            }
        }
        self.expect_rbrace()?;
        Ok(m)
    }

    pub(super) fn parse_indexed_face_set(&mut self) -> Result<(Vec<GpPnt>, Vec<i32>, Vec<GpPnt2d>), String> {
        self.expect_id()?; // IndexedFaceSet
        self.expect_lbrace()?;
        let mut coord = Vec::new();
        let mut coord_index = Vec::new();
        let mut tex_coord = Vec::new();
        while !self.peek_rbrace() {
            let fname = self.expect_id()?;
            match fname.as_str() {
                "coord" => {
                    if self.peek_id_eq("Coordinate") {
                        coord = self.parse_coordinate()?;
                    } else {
                        self.skip_value()?;
                    }
                }
                "coordIndex" => coord_index = self.parse_i32_array()?,
                "texCoord" => {
                    if self.peek_id_eq("TextureCoordinate") {
                        tex_coord = self.parse_texture_coordinate()?;
                    } else {
                        self.skip_value()?;
                    }
                }
                _ => self.skip_value()?,
            }
        }
        self.expect_rbrace()?;
        Ok((coord, coord_index, tex_coord))
    }

    pub(super) fn parse_coordinate(&mut self) -> Result<Vec<GpPnt>, String> {
        self.expect_id()?; // Coordinate
        self.expect_lbrace()?;
        let mut pts = Vec::new();
        while !self.peek_rbrace() {
            let fname = self.expect_id()?;
            match fname.as_str() {
                "point" => {
                    self.expect_lbracket()?;
                    while !self.peek_rbracket() {
                        let p = self.parse_vec3()?;
                        pts.push(GpPnt::new(p[0], p[1], p[2]));
                    }
                    self.expect_rbracket()?;
                }
                _ => self.skip_value()?,
            }
        }
        self.expect_rbrace()?;
        Ok(pts)
    }

    pub(super) fn parse_texture_coordinate(&mut self) -> Result<Vec<GpPnt2d>, String> {
        self.expect_id()?; // TextureCoordinate
        self.expect_lbrace()?;
        let mut pts = Vec::new();
        while !self.peek_rbrace() {
            let fname = self.expect_id()?;
            match fname.as_str() {
                "point" => {
                    self.expect_lbracket()?;
                    while !self.peek_rbracket() {
                        let p = self.parse_vec2()?;
                        pts.push(GpPnt2d::new(p[0], p[1]));
                    }
                    self.expect_rbracket()?;
                }
                _ => self.skip_value()?,
            }
        }
        self.expect_rbrace()?;
        Ok(pts)
    }

    pub(super) fn parse_i32_array(&mut self) -> Result<Vec<i32>, String> {
        self.expect_lbracket()?;
        let mut out = Vec::new();
        while !self.peek_rbracket() {
            match self.take() {
                Some(VTok::Num(n)) => out.push(n as i32),
                other => return Err(format!("vrml: coordIndex expects numbers, got {other:?}")),
            }
        }
        self.expect_rbracket()?;
        Ok(out)
    }
}

/// VRML `Transform` field values → `GpTrsf` (T·R·S). Rotation is axis-angle.
pub(super) fn trsf_from_vrml(t: &[f64; 3], axis: &[f64; 3], angle: f64, s: &[f64; 3]) -> GpTrsf {
    let (x, y, z) = (axis[0], axis[1], axis[2]);
    let len = (x * x + y * y + z * z).sqrt();
    let (x, y, z) = if len > 1e-12 {
        (x / len, y / len, z / len)
    } else {
        (0.0, 0.0, 1.0)
    };
    let (c, sn) = (angle.cos(), angle.sin());
    let one_c = 1.0 - c;
    let rot = GpMat::new(
        one_c * x * x + c,
        one_c * x * y - sn * z,
        one_c * x * z + sn * y,
        one_c * x * y + sn * z,
        one_c * y * y + c,
        one_c * y * z - sn * x,
        one_c * x * z - sn * y,
        one_c * y * z + sn * x,
        one_c * z * z + c,
    );
    // Fold scale in column-wise (S applied first: R·S).
    let mut m = GpMat::zero();
    for i in 0..3 {
        for j in 0..3 {
            m.m[i][j] = rot.m[i][j] * s[j];
        }
    }
    GpTrsf {
        scale: 1.0,
        shape: TrsfForm::CompoundTrsf,
        matrix: m,
        loc: GpXyz::new(t[0], t[1], t[2]),
    }
}

/// Assemble a node from an `IndexedFaceSet`: transform the coordinates by
/// `world`, fan-triangulate the `-1`-separated polygons, attach the material.
/// Texture coordinates map 1:1 to vertices when their counts match
/// (ponytail: `texCoordIndex` per-corner UVs are not handled).
pub(super) fn build_vrml_node(
    coord: &[GpPnt],
    coord_index: &[i32],
    tex_coord: Vec<GpPnt2d>,
    world: &GpTrsf,
    material: Option<Material>,
    name: &str,
) -> MeshNode {
    let verts: Vec<GpPnt> = coord.iter().map(|p| p.transformed(world)).collect();
    let mut faces: Vec<Vec<usize>> = Vec::new();
    let mut face: Vec<usize> = Vec::new();
    for &idx in coord_index {
        if idx < 0 {
            if face.len() >= 3 {
                faces.push(std::mem::take(&mut face));
            } else {
                face.clear();
            }
        } else {
            face.push(idx as usize);
        }
    }
    if face.len() >= 3 {
        faces.push(face);
    }
    let mut tris = Vec::new();
    for f in &faces {
        for i in 1..f.len().saturating_sub(1) {
            tris.push((f[0], f[i], f[i + 1]));
        }
    }
    let uv = if !tex_coord.is_empty() && tex_coord.len() == verts.len() {
        Some(tex_coord)
    } else {
        None
    };
    MeshNode {
        name: name.to_string(),
        vertices: verts,
        triangles: tris,
        transform: Some(world.clone()),
        material,
        uv,
        texture_path: None,
    }
}

// ---------------------------------------------------------------------------
// glTF internals
// ---------------------------------------------------------------------------

#[derive(Default, Clone, Copy)]
pub(super) struct BufView {
    pub(super) buffer: usize,
    pub(super) offset: usize,
    pub(super) len: usize,
}

#[derive(Default, Clone, Copy)]
pub(super) struct Accessor {
    pub(super) view: Option<usize>,
    pub(super) offset: usize,
    pub(super) comp: u32,
    pub(super) count: usize,
}

#[derive(Default, Clone, Copy)]
pub(super) struct Prim {
    pub(super) pos_view: Option<usize>,
    pub(super) pos_offset: usize,
    pub(super) pos_count: usize,
    pub(super) pos_comp: u32,
    pub(super) idx_view: Option<usize>,
    pub(super) idx_offset: usize,
    pub(super) idx_count: usize,
    pub(super) idx_comp: u32,
    pub(super) mat_idx: Option<usize>,
}

pub(super) struct GlMesh {
    pub(super) prims: Vec<Prim>,
}

pub(super) struct GlNode {
    pub(super) name: String,
    pub(super) mesh: Option<usize>,
    pub(super) trsf: GpTrsf,
    pub(super) children: Vec<usize>,
}

/// glTF material: base-color diffuse + opacity plus an optional base-color
/// texture reference (texture index → `textures[].source` → image).
pub(super) struct GlMat {
    pub(super) name: String,
    pub(super) diffuse: (f64, f64, f64),
    pub(super) opacity: f64,
    pub(super) base_color_texture: Option<usize>,
}

/// Shared resources handed to the node-graph walk.
pub(super) struct GltfCtx<'a> {
    pub(super) meshes: &'a [GlMesh],
    pub(super) buffers: &'a [Vec<u8>],
    pub(super) views: &'a [BufView],
    pub(super) materials: &'a [GlMat],
    pub(super) textures: &'a [Option<usize>],
    pub(super) images: &'a [String],
    pub(super) base_dir: &'a Path,
}

/// Build the local transform for a glTF node from `matrix` or `translation`/
/// `rotation`/`scale` (matrix wins per the spec).
pub(super) fn node_transform(nd: &JVal) -> GpTrsf {
    if let Some(mat) = nd.get("matrix").and_then(|v| v.as_arr()) {
        let m: Vec<f64> = mat.iter().filter_map(|v| v.as_num()).take(16).collect();
        if m.len() == 16 {
            return trsf_from_matrix(&m);
        }
    }
    let t = nd.get("translation").and_then(|v| v.as_arr()).map(arr3).unwrap_or([0.0; 3]);
    let r = nd.get("rotation").and_then(|v| v.as_arr()).map(arr4).unwrap_or([0.0, 0.0, 0.0, 1.0]);
    let s = nd.get("scale").and_then(|v| v.as_arr()).map(arr3).unwrap_or([1.0; 3]);
    trsf_from_trs(&t, &r, &s)
}

/// TRS → `GpTrsf`. Non-uniform scale is folded into the 3x3 linear part
/// (`scale = 1` so `transforms_xyz` never applies a second uniform scale).
pub(super) fn trsf_from_trs(t: &[f64; 3], r: &[f64; 4], s: &[f64; 3]) -> GpTrsf {
    let rot = GpQuaternion::new(r[0], r[1], r[2], r[3]).get_matrix();
    let mut m = occt_core::gp::GpMat::zero();
    for i in 0..3 {
        for j in 0..3 {
            m.m[i][j] = rot.m[i][j] * s[j];
        }
    }
    GpTrsf {
        scale: 1.0,
        shape: TrsfForm::CompoundTrsf,
        matrix: m,
        loc: GpXyz::new(t[0], t[1], t[2]),
    }
}
