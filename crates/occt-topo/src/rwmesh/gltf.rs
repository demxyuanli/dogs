use super::prelude::*;
use super::gltf_types::{GlMesh, GlNode, GltfCtx};
use super::gltf_triangulation_reader::GltfTriangulationReader;
use super::*;

use crate::tdata_std_named_data::TDataStdNamedData;
use crate::xcaf_doc_vis_material::XcafDocVisMaterial;

/// `RWGltf_GltfPrimitiveMode` (`RWGltf_GltfPrimitiveMode.hxx:20-30`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i32)]
pub enum RwGltfPrimitiveMode {
    /// `RWGltf_GltfPrimitiveMode_UNKNOWN` — unknown or invalid type.
    Unknown = -1,
    /// `GL_POINTS`.
    Points = 0,
    /// `GL_LINES`.
    Lines = 1,
    /// `GL_LINE_LOOP`.
    LineLoop = 2,
    /// `GL_LINE_STRIP`.
    LineStrip = 3,
    /// `GL_TRIANGLES` (the default when `mode` is absent).
    #[default]
    Triangles = 4,
    /// `GL_TRIANGLE_STRIP`.
    TriangleStrip = 5,
    /// `GL_TRIANGLE_FAN`.
    TriangleFan = 6,
}

impl RwGltfPrimitiveMode {
    /// `RWGltf_GltfJsonParser.cxx:1810-1814`: the accepted integer range.
    pub fn from_int(v: i64) -> Self {
        match v {
            0 => Self::Points,
            1 => Self::Lines,
            2 => Self::LineLoop,
            3 => Self::LineStrip,
            4 => Self::Triangles,
            5 => Self::TriangleStrip,
            6 => Self::TriangleFan,
            _ => Self::Unknown,
        }
    }

    /// `RWGltf_GltfJsonParser.cxx:1816-1822`: only points, lines and triangles
    /// are converted into shapes; the strip/loop/fan modes are skipped with a
    /// warning.
    pub fn is_supported(self) -> bool {
        matches!(self, Self::Points | Self::Lines | Self::Triangles)
    }

    /// `RWGltf_GltfJsonParser.cxx:1951-1958`: only triangles carry a
    /// triangulation.
    pub fn is_triangles(self) -> bool {
        self == Self::Triangles
    }
}


/// Column-major glTF 4x4 matrix → `GpTrsf`, following
/// `RWGltf_GltfJsonParser::parseTransformationMatrix`
/// (`RWGltf_GltfJsonParser.cxx:478-497`): an identity matrix leaves the
/// transform untouched, anything else goes through `gp_Trsf::SetValues`, which
/// recomputes the scale from the determinant and orthogonalizes the linear part
/// (so a matrix that is not a similarity loses its anisotropy, exactly as
/// `gp_Trsf` does).
pub(super) fn trsf_from_matrix(m: &[f64]) -> Result<GpTrsf, &'static str> {
    const IDENTITY: [f64; 16] = [
        1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0,
    ];
    if m == IDENTITY {
        return Ok(GpTrsf::identity());
    }
    let mut trsf = GpTrsf::identity();
    trsf.set_values(
        m[0], m[4], m[8], m[12],
        m[1], m[5], m[9], m[13],
        m[2], m[6], m[10], m[14],
    )?;
    Ok(trsf)
}

/// One placed primitive array: what `bindNamedShape` binds to a node's shape —
/// the mesh shape (here the deferred primitive array) plus the node location —
/// together with the attributes `RWMesh_NodeAttributes` carries
/// (`RWGltf_GltfJsonParser.cxx:2410-2505`).
pub(super) struct GltfNodeInst {
    pub(super) late: std::sync::Arc<GltfLatePrimitiveArray>,
    /// The accumulated world transform of the node. OCCT stores it in the
    /// shape's `TopLoc_Location`; this reader bakes it into the nodes when the
    /// scene is materialized.
    pub(super) world: GpTrsf,
    /// The glTF node and mesh this placement came from (`RWGltf_GltfJsonParser`
    /// builds one shape per (node, mesh, primitive array)); the document fill
    /// needs them to rebuild the node/mesh compound tree.
    pub(super) node: usize,
    pub(super) mesh: usize,
    pub(super) name: String,
    pub(super) material: Option<Material>,
    pub(super) vis_material: Option<XcafDocVisMaterial>,
    pub(super) texture_path: Option<String>,
    pub(super) named_data: Option<TDataStdNamedData>,
}

/// Depth-first walk of the node graph; records one [`GltfNodeInst`] per
/// (node, mesh, primitive array) with the accumulated world transform, the
/// node's attributes and the shared deferred array.
///
/// When `ctx.converter` is not empty, `TransformTransformation` is applied to
/// the accumulated location (`RWGltf_GltfJsonParser.cxx:493`, `:637`). The
/// positions themselves are converted by the triangulation reader
/// (`RWGltf_TriangulationReader.cxx:325`), not here.
pub(super) fn collect_gltf_nodes(
    nodes: &[GlNode],
    idx: usize,
    world: &GpTrsf,
    ctx: &GltfCtx,
    visited: &mut Vec<bool>,
    out: &mut Vec<GltfNodeInst>,
) -> Result<(), String> {
    // `findNodeShape` (`RWGltf_GltfJsonParser.cxx:1573-1576`): a node that was
    // already turned into a shape (another scene root, or a second parent in
    // the graph) reuses that shape; it is not placed twice.
    if visited[idx] {
        return Ok(());
    }
    visited[idx] = true;
    let n = &nodes[idx];
    // `RWGltf_GltfJsonParser.cxx:493`, `:637`: each node's *local* location is
    // conjugated by the converter, so a chain of converted locals still composes
    // to `C * W1 * W2 * ... * C^-1`.
    let mut local = n.trsf.clone();
    if !ctx.converter.is_empty() {
        ctx.converter.transform_transformation(&mut local);
    }
    let world = world.multiplied(&local);
    // glTF 2.0 carries one `mesh`; glTF 1.0 carries a `meshes` array and OCCT
    // walks both (`RWGltf_GltfJsonParser.cxx:1644-1707`), adding one mesh shape
    // per entry to the node compound.
    let mut mesh_ids: Vec<usize> = Vec::with_capacity(n.meshes.len() + 1);
    if let Some(m) = n.mesh {
        mesh_ids.push(m);
    }
    mesh_ids.extend(n.meshes.iter().copied());
    for m in mesh_ids {
        let mesh = ctx.meshes.get(m).ok_or_else(|| format!("gltf: node {idx}: mesh {m} out of range"))?;
        // `gltfParseMesh` (`cxx:1726-1776`): the mesh's primitive arrays are
        // parsed once and shared through `myShapeMap[ShapeMapGroup_PrimArray]`.
        let arrays = ctx.cache.prim_arrays(m, mesh, ctx)?;
        for (prim, late) in mesh.prims.iter().zip(arrays.iter()) {
            // `bindNamedShape` (`RWGltf_GltfJsonParser.cxx:2410-2447`):
            //  * the node's `name` wins; glTF 1.0 falls back to the node id;
            //  * an unnamed node takes the mesh name when
            //    `myUseMeshNameAsFallback` is set.
            // `mesh_{idx}` is this port's container-level uniqueness fallback for
            // a shape that would otherwise carry no name at all.
            let name = if !n.name.is_empty() {
                n.name.clone()
            } else if ctx.is_gltf1 && !n.id.is_empty() {
                n.id.clone()
            } else if ctx.use_mesh_name_as_fallback && !mesh.name.is_empty() {
                mesh.name.clone()
            } else {
                format!("mesh_{idx}")
            };

            let (material, vis_material, texture_path) = prim.mat_idx
                .and_then(|mi| ctx.materials.get(mi))
                .map(|m| {
                    // `RWGltf_CafReader` binds the `XCAFDoc_VisMaterial` to the
                    // node style; the `RWMesh_MaterialMap`-side view reads the
                    // same attribute back (`material_view`).
                    let vis = m.vis.clone();
                    let mat = material_view(&vis);
                    let tex = mat.diffuse_texture.clone();
                    (Some(mat), Some(vis), tex)
                })
                .unwrap_or((None, None, None));
            // `RWGltf_GltfJsonParser.cxx:2019-2027`: a line shape carries a
            // color (`Style.SetColorCurv`), not an `XCAFDoc_VisMaterial`.
            let (material, vis_material, texture_path) = if prim.mode == RwGltfPrimitiveMode::Lines
            {
                let color = material.as_ref().map(|m| {
                    let mut m = m.clone();
                    m.diffuse_texture = None;
                    m
                });
                (color, None, None)
            } else {
                (material, vis_material, texture_path)
            };

            out.push(GltfNodeInst {
                late: late.clone(),
                world: world.clone(),
                node: idx,
                mesh: m,
                name,
                material,
                vis_material,
                texture_path,
                // `bindNamedShape` stores `theExtras` on the shape attributes
                // (`cxx:2422-2424`); the node's own `extras` wins over the
                // mesh's, which the C++ binds to the Product label instead
                // (`cxx:1772-1774`).
                named_data: n.extras.clone().or_else(|| mesh.extras.clone()),
            });
        }
    }
    for &c in &n.children {
        collect_gltf_nodes(nodes, c, &world, ctx, visited, out)?;
    }
    Ok(())
}

/// Materialize one placed primitive array into a [`MeshNode`].
///
/// `RWGltf_CafReader::addShapeIntoDoc` (`cxx:470-688`) puts the shape in the
/// document with its location; the flat [`MeshScene`] has no location, so the
/// node's world transform is applied to the resolved triangulation's nodes —
/// the same point set (`RWGltf_CafReader.cxx:489-519` also scales nodes there).
/// An unresolved array (`myToSkipLateDataLoading`) yields the proxy shape: the
/// face exists (`cxx:2000-2007`) but holds no nodes yet.
pub(super) fn inst_to_node(inst: &GltfNodeInst) -> MeshNode {
    inst_to_node_at(inst, &inst.world)
}

/// The primitive array's own geometry, without the placement — the shape
/// `gltfParsePrimArray` builds and `myShapeMap[ShapeMapGroup_PrimArray]` shares
/// (`RWGltf_GltfJsonParser.cxx:1969-2013`). OCCT keeps it at the identity
/// location and places it through the node's `TopLoc_Location`; this crate bakes
/// the placement into the geometry, so the unplaced shape is what the document
/// fill keeps on the Product label (`XCAFDoc_ShapeTool::AddComponent(theAssembly,
/// theShape, ...)` strips the location, `cxx:999-1003`).
pub(super) fn inst_to_local_node(inst: &GltfNodeInst) -> MeshNode {
    inst_to_node_at(inst, &GpTrsf::identity())
}

fn inst_to_node_at(inst: &GltfNodeInst, world: &GpTrsf) -> MeshNode {
    let resolved = inst.late.resolved();
    let (vertices, triangles, uv, normals, polyline) = match &resolved {
        Some(r) => {
            let vertices = r
                .triangulation
                .nodes
                .iter()
                .map(|p| p.transformed(world))
                .collect();
            let triangles = r
                .triangulation
                .triangles
                .iter()
                .map(|t| (t.n0, t.n1, t.n2))
                .collect();
            let uv = r
                .triangulation
                .uv_nodes
                .as_ref()
                .map(|u| u.iter().map(|&(x, y)| GpPnt2d::new(x, y)).collect());
            // The placement is baked into the vertices, so the stored normals
            // take the same rotation (`gp_Trsf::Transforms` on a direction).
            let normals = r.triangulation.normals.as_ref().map(|ns| {
                ns.iter()
                    .map(|n| {
                        let mut v = n.coord;
                        world.transforms_xyz_dir(&mut v);
                        GpPnt::new(v.x(), v.y(), v.z())
                    })
                    .collect()
            });
            (vertices, triangles, uv, normals, r.edges.clone())
        }
        None => (Vec::new(), Vec::new(), None, None, Vec::new()),
    };
    MeshNode {
        name: inst.name.clone(),
        vertices,
        triangles,
        mode: inst.late.prim_mode(),
        polyline,
        transform: Some(inst.world.clone()),
        material: inst.material.clone(),
        vis_material: inst.vis_material.clone(),
        uv,
        normals,
        texture_path: inst.texture_path.clone(),
        named_data: inst.named_data.clone(),
        late: Some(inst.late.clone()),
    }
}

/// Shape-level reuse caches of `RWGltf_GltfJsonParser`: `findMeshShape`
/// (`RWGltf_GltfJsonParser.cxx:1733`) and `myShapeMap[ShapeMapGroup_PrimArray]`
/// (`cxx:1847-1877`). OCCT builds a mesh/primitive-array shape once and lets
/// every node that references it share the same `TopoDS_Shape`; here the shape
/// is the deferred primitive array, shared by `Arc`.
///
/// The primitive-array key is `aPrimArrayId` (`cxx:1847-1856`): the mode, the
/// indices accessor and every attribute accessor. OCCT additionally binds
/// `aPrimArrayIdWithMat` and copies the shape when only the material differs
/// (`cxx:1858-1877`); the deferred array does not depend on the material, and
/// every emitted node already binds its own material, so the two keys collapse
/// into one here.
#[derive(Default)]
pub(super) struct PrimCache {
    meshes: std::cell::RefCell<HashMap<usize, Vec<std::sync::Arc<GltfLatePrimitiveArray>>>>,
    prims: std::cell::RefCell<HashMap<String, std::sync::Arc<GltfLatePrimitiveArray>>>,
    /// `RWGltf_CafReader::myFaceList` (`RWGltf_CafReader.cxx:2000-2006`): the
    /// triangle primitive arrays, in creation order, handed to
    /// [`super::gltf_late_data::read_late_data`].
    faces: std::cell::RefCell<Vec<std::sync::Arc<GltfLatePrimitiveArray>>>,
}

impl PrimCache {
    /// `findMeshShape` (`cxx:1733-1736`): a mesh's primitive arrays are built
    /// once per mesh id.
    pub(super) fn prim_arrays(
        &self,
        mesh_idx: usize,
        mesh: &GlMesh,
        ctx: &GltfCtx,
    ) -> Result<Vec<std::sync::Arc<GltfLatePrimitiveArray>>, String> {
        if let Some(data) = self.meshes.borrow().get(&mesh_idx) {
            return Ok(data.clone());
        }
        let mut out = Vec::with_capacity(mesh.prims.len());
        for (prim_idx, prim) in mesh.prims.iter().enumerate() {
            let id = prim.array_id();
            let cached = self.prims.borrow().get(&id).cloned();
            let late = match cached {
                Some(late) => late,
                None => {
                    // `gltfParsePrimArray` (`cxx:1879-1885`).
                    let late = std::sync::Arc::new(gltf_types::build_late_array(
                        ctx, mesh_idx, prim_idx, &mesh.name, prim,
                    )?);
                    // `cxx:1960-1967`: deferred loading only works for
                    // triangulations; points and lines are filled right away
                    // because their shape needs the actual nodes.
                    if prim.mode != RwGltfPrimitiveMode::Triangles
                        && !late.data().is_empty()
                    {
                        let reader =
                            GltfTriangulationReader::new(ctx.converter);
                        let resolved = reader.load_all_data(&late)?;
                        late.set_resolved(resolved);
                    } else if prim.mode == RwGltfPrimitiveMode::Triangles {
                        self.faces.borrow_mut().push(late.clone());
                    }
                    self.prims.borrow_mut().insert(id, late.clone());
                    late
                }
            };
            out.push(late);
        }
        self.meshes.borrow_mut().insert(mesh_idx, out.clone());
        Ok(out)
    }

    /// `aDoc.FaceList()` (`RWGltf_CafReader.cxx:347`): the proxy faces whose
    /// triangulation the late-data pass fills.
    pub(super) fn faces(&self) -> Vec<std::sync::Arc<GltfLatePrimitiveArray>> {
        self.faces.borrow().clone()
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
    /// JSON integer literal. `RapidJSON`'s `IsInt()` and `IsDouble()` decide
    /// between `TDataStd_NamedData::SetInteger` and `SetReal`
    /// (`RWGltf_GltfJsonParser.cxx:213-233`) and between an integer and a real
    /// array (`cxx:273-296`), so the literal form is kept.
    Int(i64),
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
            JVal::Int(n) => Some(*n as f64),
            JVal::Num(n) => Some(*n),
            _ => None,
        }
    }
    pub(super) fn as_bool(&self) -> Option<bool> {
        match self {
            JVal::Bool(b) => Some(*b),
            _ => None,
        }
    }
    /// `RapidJSON::Value::IsInt()`: the literal has no fraction or exponent.
    pub(super) fn as_int(&self) -> Option<i64> {
        match self {
            JVal::Int(n) => Some(*n),
            _ => None,
        }
    }
    /// `RapidJSON::Value::IsDouble()`.
    pub(super) fn as_double(&self) -> Option<f64> {
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
            Some(c) if c == b'-' || c.is_ascii_digit() => self.number(),
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
    /// JSON number: an integer literal stays [`JVal::Int`] so callers can
    /// distinguish `IsInt()` from `IsDouble()` like RapidJSON does.
    pub(super) fn number(&mut self) -> Result<JVal, String> {
        let start = self.i;
        while self.i < self.s.len()
            && (self.s[self.i].is_ascii_digit()
                || matches!(self.s[self.i], b'-' | b'+' | b'.' | b'e' | b'E'))
        {
            self.i += 1;
        }
        let text = std::str::from_utf8(&self.s[start..self.i])
            .map_err(|_| "json: bad number".to_string())?;
        if !text.contains(['.', 'e', 'E']) {
            if let Ok(n) = text.parse::<i64>() {
                return Ok(JVal::Int(n));
            }
        }
        text.parse::<f64>()
            .map(JVal::Num)
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
