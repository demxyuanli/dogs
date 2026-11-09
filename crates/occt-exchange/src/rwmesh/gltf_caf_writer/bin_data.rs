//! Binary data stage of `RWGltf_CafWriter`: `RWGltf_GltfAccessor`,
//! `RWGltf_GltfBufferView`, `RWGltf_GltfFace` and the four array passes of
//! `writeBinData` (`RWGltf_CafWriter.cxx:811-1207`).
//!
//! Source: `RWGltf_GltfAccessor.hxx`, `RWGltf_GltfBufferView.hxx`,
//! `RWGltf_CafWriter.cxx:328-620` (`saveNodes`, `saveNormals`,
//! `saveTextCoords`, `saveIndices`, `saveTriangleIndices`) and `cxx:755-809`
//! (`writeShapesToBin`).
//!
//! The C++ writes into a `std::ostream` (a `.bin` sidecar, or a `.bin.tmp` the
//! GLB assembly later appends); the port accumulates the same byte stream in a
//! `Vec<u8>`, so `tellp()` is `bin.len()`.
//!
//! Only the face (`TopAbs_FACE`) arm is reachable here: the node lists this
//! crate's document readers hand to the writer carry triangulated faces only,
//! and the exported document is a STEP/shape document whose leaves are solids
//! and shells. The edge/vertex iterators (`RWMesh_EdgeIterator` /
//! `RWMesh_VertexIterator` arms of `writeBinData`, `cxx:1004-1046`) have no
//! source in the port's node list and are not reproduced.

use occt_core::gp::GpXyz;
use occt_core::poly::triangulation::Triangle;

use crate::rwmesh::coordinate_system::RwMeshCoordinateSystemConverter;
use crate::rwmesh::ObjFaceNode;
use crate::xcaf_prs::XcafPrsStyle;

/// `RWGltf_GltfAccessor::INVALID_ID`.
pub(crate) const INVALID_ID: i32 = -1;

/// `RWGltf_GltfBufferViewTarget_ARRAY_BUFFER` (`RWGltf_GltfBufferView.hxx`).
pub(crate) const TARGET_ARRAY_BUFFER: i32 = 34962;
/// `RWGltf_GltfBufferViewTarget_ELEMENT_ARRAY_BUFFER`.
pub(crate) const TARGET_ELEMENT_ARRAY_BUFFER: i32 = 34963;

/// `RWGltf_GltfAccessorCompType` (`RWGltf_GltfAccessorCompType.hxx`): the GL
/// component enumerants, which glTF stores verbatim.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum GltfAccessorCompType {
    /// `RWGltf_GltfAccessorCompType_UInt16` (5123).
    UInt16,
    /// `RWGltf_GltfAccessorCompType_UInt32` (5125).
    UInt32,
    /// `RWGltf_GltfAccessorCompType_Float32` (5126).
    Float32,
}

impl GltfAccessorCompType {
    /// The enumerant glTF writes into `componentType`.
    pub(crate) fn value(self) -> i32 {
        match self {
            // `GL_UNSIGNED_SHORT`.
            Self::UInt16 => 5123,
            // `GL_UNSIGNED_INT`.
            Self::UInt32 => 5125,
            // `GL_FLOAT`.
            Self::Float32 => 5126,
        }
    }
}

/// `RWGltf_GltfAccessorLayout` (`RWGltf_GltfAccessorLayout.hxx`): only the
/// three layouts the exporter emits.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum GltfAccessorLayout {
    /// `Vec3` — `"VEC3"`.
    Vec3,
    /// `Vec2` — `"VEC2"`.
    Vec2,
    /// `Scalar` — `"SCALAR"`.
    Scalar,
}

impl GltfAccessorLayout {
    /// The `type` string glTF writes.
    pub(crate) fn as_str(self) -> &'static str {
        match self {
            Self::Vec3 => "VEC3",
            Self::Vec2 => "VEC2",
            Self::Scalar => "SCALAR",
        }
    }
}

/// `RWGltf_GltfPrimitiveMode` (`RWGltf_GltfPrimitiveMode.hxx`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)] // `Points`/`Lines` need `RWMesh_EdgeIterator`/`VertexIterator` (unported)
pub(crate) enum GltfPrimitiveMode {
    /// `Points` (0).
    Points,
    /// `Lines` (1).
    Lines,
    /// `Triangles` (4).
    Triangles,
}

impl GltfPrimitiveMode {
    /// The `mode` value glTF writes.
    pub(crate) fn value(self) -> i32 {
        match self {
            Self::Points => 0,
            Self::Lines => 1,
            Self::Triangles => 4,
        }
    }
}

/// `RWGltf_GltfAccessor` (`RWGltf_GltfAccessor.hxx`): one `accessors` entry.
#[derive(Debug, Clone)]
pub(crate) struct GltfAccessor {
    /// `Id`.
    pub id: i32,
    /// `ByteOffset`.
    pub byte_offset: i64,
    /// `Count`.
    pub count: i64,
    /// `Type`.
    pub layout: GltfAccessorLayout,
    /// `ComponentType`.
    pub comp_type: GltfAccessorCompType,
    /// `BndBox` — the `Bnd_Box` the position arm fills; `None` while void.
    pub bnd_box: Option<([f64; 3], [f64; 3])>,
}

impl GltfAccessor {
    /// `RWGltf_GltfAccessor()` (`hxx`): every field zero/`INVALID_ID` and an
    /// empty box.
    pub(crate) fn new() -> Self {
        Self {
            id: INVALID_ID,
            byte_offset: 0,
            count: 0,
            layout: GltfAccessorLayout::Scalar,
            comp_type: GltfAccessorCompType::Float32,
            bnd_box: None,
        }
    }
}

/// `Bnd_Box::Add(gp_Pnt)` followed by `Bnd_Box::IsValid()`: the first point
/// seeds the box, later points extend it.
fn bnd_box_add(bnd_box: &mut Option<([f64; 3], [f64; 3])>, point: [f64; 3]) {
    match bnd_box {
        Some((min, max)) => {
            for i in 0..3 {
                min[i] = min[i].min(point[i]);
                max[i] = max[i].max(point[i]);
            }
        }
        None => *bnd_box = Some((point, point)),
    }
}

/// `RWGltf_GltfBufferView` (`RWGltf_GltfBufferView.hxx`).
#[derive(Debug, Clone, Copy)]
pub(crate) struct GltfBufferView {
    /// `Id`.
    pub id: i32,
    /// `ByteOffset`.
    pub byte_offset: i64,
    /// `ByteLength`.
    pub byte_length: i64,
    /// `ByteStride` — 0 when the view does not declare one (indices).
    pub byte_stride: i32,
    /// `Target`.
    pub target: i32,
}

impl GltfBufferView {
    /// One of the four array buffer views `writeBinData` sets up
    /// (`cxx:824-846`): `INVALID_ID`, zero offsets and the given stride/target.
    pub(crate) fn new(byte_stride: i32, target: i32) -> Self {
        Self { id: INVALID_ID, byte_offset: 0, byte_length: 0, byte_stride, target }
    }
}

/// `RWGltf_GltfFace` (`RWGltf_CafWriter.hxx:472-560`): the per-shape record the
/// four binary passes fill in, plus the identity of the source face.
#[derive(Debug, Clone)]
pub(crate) struct GltfFaceData {
    /// `Shape` — the face's `TShape` key (`GeometryRegistry::shape_key`). The
    /// map key holds the same value; kept here because `RWGltf_GltfFace` owns
    /// its `Shape` and `writePrimArray` switches on its type (`cxx:2077-2089`).
    #[allow(dead_code)]
    pub shape_key: usize,
    /// Index of the owning node in the writer's node list.
    pub node: usize,
    /// Index of the face inside that node (`GltfCafNode::faces`).
    pub face: usize,
    /// `NodePos`.
    pub node_pos: GltfAccessor,
    /// `NodeNorm`.
    pub node_norm: GltfAccessor,
    /// `NodeUV`.
    pub node_uv: GltfAccessor,
    /// `Indices`.
    pub indices: GltfAccessor,
    /// `NbIndexedNodes` — reset at the start of the position pass, then the
    /// running node offset of the indices pass.
    pub nb_indexed_nodes: i32,
}

impl GltfFaceData {
    /// Fresh record for one document face.
    pub(crate) fn new(shape_key: usize, node: usize, face: usize) -> Self {
        Self {
            shape_key,
            node,
            face,
            node_pos: GltfAccessor::new(),
            node_norm: GltfAccessor::new(),
            node_uv: GltfAccessor::new(),
            indices: GltfAccessor::new(),
            nb_indexed_nodes: 0,
        }
    }
}

/// One `myBinDataMap` entry: a styled shape and the `RWGltf_GltfFace` list it
/// serializes into. Without `toMergeFaces` the list always holds one element
/// (`RWGltf_CafWriter.cxx:733-749`).
#[derive(Debug, Clone)]
pub(crate) struct GltfFaceGroup {
    /// `RWGltf_StyledShape::Shape`.
    pub shape_key: usize,
    /// `RWGltf_StyledShape::Style` (`RWGltf_CafWriter.hxx:491-495`): the two
    /// key components of the map entry.
    pub style: XcafPrsStyle,
    /// `RWGltf_GltfFace` list, in document order.
    pub faces: Vec<GltfFaceData>,
}

/// `RWGltf_CafWriter::getShapeType` for a face-shaped `GltfFace`
/// (`cxx:286-316`): a document face is always `TopAbs_FACE`, so
/// [`Self::has_triangulation`] is TRUE and [`GltfPrimitiveMode::Triangles`] is
/// the primitive mode.
pub(crate) const PRIMITIVE_MODE: GltfPrimitiveMode = GltfPrimitiveMode::Triangles;

/// Append one little-endian `f32`.
fn write_f32(bin: &mut Vec<u8>, value: f32) {
    bin.extend_from_slice(&value.to_le_bytes());
}

/// `writeVec3(std::ostream&, const gp_XYZ&)` (`cxx:57-63`): three floats.
fn write_vec3(bin: &mut Vec<u8>, value: &GpXyz) {
    write_f32(bin, value.x() as f32);
    write_f32(bin, value.y() as f32);
    write_f32(bin, value.z() as f32);
}

/// `writeVec2(std::ostream&, const gp_XY&)` (`cxx:72-77`): two floats.
fn write_vec2(bin: &mut Vec<u8>, x: f64, y: f64) {
    write_f32(bin, x as f32);
    write_f32(bin, y as f32);
}

/// `RWGltf_CafWriter::saveNodes` (`RWGltf_CafWriter.cxx:328-368`).
pub(crate) fn save_nodes(
    bin: &mut Vec<u8>,
    data: &mut GltfFaceData,
    face: &ObjFaceNode,
    accessor_nb: &mut i32,
    cs_trsf: &RwMeshCoordinateSystemConverter,
    buff_view_offset: i64,
) {
    if data.node_pos.id == INVALID_ID {
        data.node_pos.id = *accessor_nb;
        *accessor_nb += 1;
        data.node_pos.byte_offset = bin.len() as i64 - buff_view_offset;
        data.node_pos.layout = GltfAccessorLayout::Vec3;
        data.node_pos.comp_type = GltfAccessorCompType::Float32;
    }
    data.node_pos.count += face.vertices.len() as i64;

    for vertex in &face.vertices {
        let mut node = GpXyz::new(vertex.x(), vertex.y(), vertex.z());
        cs_trsf.transform_position(&mut node);
        bnd_box_add(&mut data.node_pos.bnd_box, [node.x(), node.y(), node.z()]);
        write_vec3(bin, &node);
    }
}

/// `RWGltf_CafWriter::saveNormals` (`RWGltf_CafWriter.cxx:372-415`).
pub(crate) fn save_normals(
    bin: &mut Vec<u8>,
    data: &mut GltfFaceData,
    face: &ObjFaceNode,
    accessor_nb: &mut i32,
    cs_trsf: &RwMeshCoordinateSystemConverter,
    buff_view_offset: i64,
) {
    if !face.has_normals {
        return;
    }

    if data.node_norm.id == INVALID_ID {
        data.node_norm.id = *accessor_nb;
        *accessor_nb += 1;
        data.node_norm.byte_offset = bin.len() as i64 - buff_view_offset;
        data.node_norm.layout = GltfAccessorLayout::Vec3;
        data.node_norm.comp_type = GltfAccessorCompType::Float32;
    }
    data.node_norm.count += face.vertices.len() as i64;

    for i in 0..face.vertices.len() {
        // `RWMesh_FaceIterator::NormalTransformed` always has a value when
        // `HasNormals()` is set; the default keeps the stream complete.
        let normal = face.normals.get(i).copied().unwrap_or_default();
        let mut vec = GpXyz::new(
            normal.x() as f32 as f64,
            normal.y() as f32 as f64,
            normal.z() as f32 as f64,
        );
        cs_trsf.transform_normal(&mut vec);
        write_vec3(bin, &vec);
    }
}

/// `RWGltf_GltfMaterialMap::baseColorTexture` (`RWGltf_GltfMaterialMap.cxx:29-47`)
/// applied to the UV export gate: does the material reference any texture?
fn material_has_texture(material: &crate::xcaf_doc_vis_material::XcafDocVisMaterial) -> bool {
    let pbr = material.pbr_material();
    // `RWGltf_GltfMaterialMap::baseColorTexture`: PBR base color, else the
    // common material's diffuse texture.
    let base_color = pbr
        .base_color_texture
        .as_ref()
        .or_else(|| material.common_material().diffuse_texture.as_ref());
    base_color.is_some()
        || pbr.metallic_roughness_texture.is_some()
        || pbr.emissive_texture.is_some()
        || pbr.occlusion_texture.is_some()
        || pbr.normal_texture.is_some()
}

/// `RWGltf_CafWriter::saveTextCoords` (`RWGltf_CafWriter.cxx:420-481`).
pub(crate) fn save_text_coords(
    bin: &mut Vec<u8>,
    data: &mut GltfFaceData,
    face: &ObjFaceNode,
    accessor_nb: &mut i32,
    is_forced_uv_export: bool,
    buff_view_offset: i64,
) {
    if !face.has_tex_coords {
        return;
    }
    if !is_forced_uv_export {
        // `if (theFaceIter.Style().Material().IsNull()) return;` and the
        // "no texture at all" arm (`cxx:431-444`).
        let Some(material) = face.style.material() else {
            return;
        };
        if !material_has_texture(material) {
            return;
        }
    }

    if data.node_uv.id == INVALID_ID {
        data.node_uv.id = *accessor_nb;
        *accessor_nb += 1;
        data.node_uv.byte_offset = bin.len() as i64 - buff_view_offset;
        data.node_uv.layout = GltfAccessorLayout::Vec2;
        data.node_uv.comp_type = GltfAccessorCompType::Float32;
    }
    data.node_uv.count += face.vertices.len() as i64;

    for i in 0..face.vertices.len() {
        let (x, y) = match face.tex_coords.get(i) {
            // `aTexCoord.SetY(1.0 - aTexCoord.Y())` (`cxx:469`).
            Some(uv) => (uv.x(), 1.0 - uv.y()),
            None => (0.0, 1.0),
        };
        write_vec2(bin, x, y);
    }
}

/// `RWGltf_CafWriter::saveIndices` (`RWGltf_CafWriter.cxx:575-621`) together
/// with `saveTriangleIndices` (`cxx:484-517`), the only arm a face document
/// reaches.
pub(crate) fn save_indices(
    bin: &mut Vec<u8>,
    data: &mut GltfFaceData,
    face: &ObjFaceNode,
    accessor_nb: &mut i32,
    buff_view_offset: i64,
) {
    if data.indices.id == INVALID_ID {
        data.indices.id = *accessor_nb;
        *accessor_nb += 1;
        data.indices.byte_offset = bin.len() as i64 - buff_view_offset;
        data.indices.layout = GltfAccessorLayout::Scalar;
        // `NodePos.Count > numeric_limits<uint16_t>::max()` (`cxx:585-589`).
        data.indices.comp_type = if data.node_pos.count > u16::MAX as i64 {
            GltfAccessorCompType::UInt32
        } else {
            GltfAccessorCompType::UInt16
        };
    }

    // `const int aNodeFirst = theGltfFace.NbIndexedNodes - theFaceIter.ElemLower();`
    // (`cxx:489`). `Poly_Triangulation` elements are 1-based, so with the
    // 0-based triangles this port stores the shift is the running node count
    // itself.
    let node_first = data.nb_indexed_nodes;
    data.nb_indexed_nodes += face.vertices.len() as i32;
    data.indices.count += (face.triangles.len() * 3) as i64;

    for Triangle { n0, n1, n2 } in &face.triangles {
        let a = *n0 as i64 + node_first as i64;
        let b = *n1 as i64 + node_first as i64;
        let c = *n2 as i64 + node_first as i64;
        match data.indices.comp_type {
            GltfAccessorCompType::UInt16 => {
                for index in [a, b, c] {
                    bin.extend_from_slice(&(index as u16).to_le_bytes());
                }
            }
            _ => {
                for index in [a, b, c] {
                    bin.extend_from_slice(&(index as u32).to_le_bytes());
                }
            }
        }
    }
}
