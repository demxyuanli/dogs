//! glTF data structures and the primitive-array descriptor builder.
//!
//! Source: `RWGltf_GltfJsonParser`'s `gltfParseAccessor` (`cxx:2039-2198`),
//! `gltfParseBufferView` (`cxx:2204-2270`) and `gltfParseBuffer`
//! (`cxx:2276-2367`), plus the tree the parser builds for `nodes`, `meshes`
//! and `materials`.
//!
//! Nothing here decodes vertex bytes. Each primitive becomes a list of
//! [`GltfPrimArrayData`] descriptors that record where an accessor's bytes
//! start and whether they are already in memory (base64 buffer) or still have
//! to be read from a file (`RWGltf_GltfJsonParser.cxx:2294-2365`). Decoding is
//! the job of [`super::gltf_triangulation_reader`].

use std::cell::RefCell;
use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;

use occt_core::gp::{GpQuaternion, GpTrsf};

use super::coordinate_system::RwMeshCoordinateSystemConverter;
use super::gltf::{trsf_from_matrix, JVal, RwGltfPrimitiveMode};
use super::gltf_prim_array::{
    parse_attrib_type, GltfAccessor, GltfAccessorCompType, GltfAccessorLayout, GltfArrayType,
    GltfLatePrimitiveArray, GltfPrimArrayData,
};
use super::gltf::PrimCache;
use crate::tdata_std_named_data::TDataStdNamedData;
use crate::xcaf_doc_vis_material::XcafDocVisMaterial;

/// `RWGltf_GltfBufferView` (`RWGltf_GltfBufferView.hxx`): the fields
/// `gltfParseBufferView` validates.
#[derive(Default, Clone, Copy)]
pub(super) struct BufView {
    pub(super) buffer: usize,
    pub(super) offset: usize,
    pub(super) len: usize,
    /// `gltfParseBufferView`'s `aBuffView.ByteStride` (`cxx:2227`).
    pub(super) stride: usize,
    /// `gltfParseBufferView`'s `aBuffView.Target` (`cxx:2229-2237`), `0` when
    /// the bufferView does not define one.
    pub(super) target: i64,
}

/// The accessor fields `gltfParseAccessor` reads before it walks the
/// bufferView (`RWGltf_GltfJsonParser.cxx:2046-2117`).
#[derive(Default, Clone, Copy)]
pub(super) struct Accessor {
    pub(super) view: Option<usize>,
    pub(super) offset: usize,
    /// `RWGltf_GltfAccessor::ByteStride` from the accessor itself; glTF 1.0
    /// kept `byteStride` here (`cxx:2100`). The reader uses the bufferView
    /// stride, which gltfParseBuffer writes over it (`cxx:2300`).
    pub(super) stride: usize,
    pub(super) comp: i64,
    pub(super) count: usize,
    /// `RWGltf_GltfAccessor::Type` (`RWGltf_GltfAccessorLayout.hxx`).
    pub(super) layout: GltfAccessorLayout,
}

/// One mesh primitive array (`gltfParsePrimArray`, `cxx:1782-2035`).
#[derive(Default, Clone)]
pub(super) struct Prim {
    /// `RWGltf_GltfLatePrimitiveArray::PrimitiveMode`
    /// (`RWGltf_GltfJsonParser.cxx:1797-1822`).
    pub(super) mode: RwGltfPrimitiveMode,
    /// `anAttribs` in JSON member order, as (array type, accessor key). OCCT
    /// ignores an attribute whose name maps to `RWGltf_GltfArrayType_UNKNOWN`
    /// (`cxx:1899-1907`) and keeps the rest.
    pub(super) attribs: Vec<(GltfArrayType, usize)>,
    /// `indices` accessor key (`cxx:1931-1950`).
    pub(super) indices: Option<usize>,
    /// `POSITION` accessor key — part of `aPrimArrayId` (`cxx:1847-1857`)
    /// and the `hasPositions` check (`cxx:1887`, `:1920-1928`).
    pub(super) pos_acc: Option<usize>,
    pub(super) mat_idx: Option<usize>,
}

impl Prim {
    /// `aPrimArrayId` (`RWGltf_GltfJsonParser.cxx:1847-1857`): the mode, the
    /// indices key and every attribute's accessor key.
    pub(super) fn array_id(&self) -> String {
        let mut id = format!("{}:", self.mode as i32);
        if let Some(i) = self.indices {
            id.push_str(&i.to_string());
        }
        for (ty, acc) in &self.attribs {
            id.push(':');
            id.push_str(array_type_key(*ty));
            id.push(':');
            id.push_str(&acc.to_string());
        }
        id
    }
}

/// The glTF attribute name an array type came from; `indices` has no name.
fn array_type_key(ty: GltfArrayType) -> &'static str {
    match ty {
        GltfArrayType::Position => "POSITION",
        GltfArrayType::Normal => "NORMAL",
        GltfArrayType::Color => "COLOR",
        GltfArrayType::TCoord0 => "TEXCOORD_0",
        GltfArrayType::TCoord1 => "TEXCOORD_1",
        GltfArrayType::Joint => "JOINT",
        GltfArrayType::Weight => "WEIGHT",
        GltfArrayType::Indices => "indices",
        GltfArrayType::Unknown => "?",
    }
}

pub(super) struct GlMesh {
    /// The mesh's key in the `meshes` collection: its index in glTF 2.0, its id
    /// in glTF 1.0. glTF 1.0 falls back to it as the mesh name
    /// (`RWGltf_GltfJsonParser.cxx:2410-2418`).
    pub(super) id: String,
    /// glTF `mesh.name`. `RWGltf_GltfJsonParser` passes it as the mesh's
    /// `RawName`, which the node name falls back to when the node itself has no
    /// name (`cxx:2011-2028`, `cxx:2444-2447`).
    pub(super) name: String,
    /// `extras` of the mesh (`gltfParseMesh`, `cxx:1772-1774`): the
    /// `TDataStd_NamedData` `bindMeshShape` stores on the Product label
    /// (`cxx:2350-2365`).
    pub(super) extras: Option<TDataStdNamedData>,
    pub(super) prims: Vec<Prim>,
}

pub(super) struct GlNode {
    /// The node's `name`, empty when absent.
    pub(super) name: String,
    /// The node's key in the `nodes` collection — its index in glTF 2.0, its id
    /// in glTF 1.0. glTF 1.0 falls back to it as the node name
    /// (`RWGltf_GltfJsonParser.cxx:2412-2418`).
    pub(super) id: String,
    pub(super) mesh: Option<usize>,
    /// glTF 1.0 `meshes` array (`RWGltf_GltfJsonParser.cxx:1644-1676`): a node
    /// may reference several meshes; glTF 2.0 replaced this with `mesh`.
    pub(super) meshes: Vec<usize>,
    pub(super) trsf: GpTrsf,
    pub(super) children: Vec<usize>,
    /// `extras` of the node (`gltfParseSceneNode`, `cxx:1621-1622`): attached to
    /// the node's shape as `RWMesh_NodeAttributes::NamedData`.
    pub(super) extras: Option<TDataStdNamedData>,
}

/// glTF material table entry: the bound `XCAFDoc_VisMaterial` plus the id and
/// name glTF assigned. Source: `gltfBindMaterial`
/// (`RWGltf_GltfJsonParser.cxx:867-947`) and `myMaterials`.
pub(super) struct GlMat {
    pub(super) id: String,
    pub(super) name: String,
    /// `handle(XCAFDoc_VisMaterial)` — carries the PBR and/or common half.
    pub(super) vis: XcafDocVisMaterial,
}

/// The decoded bodies of base64 buffers, keyed by buffer key —
/// `RWGltf_GltfJsonParser::myDecodedBuffers` (`cxx:2319-2331`).
#[derive(Default)]
pub(super) struct DecodedBuffers {
    buffers: RefCell<HashMap<String, Arc<Vec<u8>>>>,
}

impl DecodedBuffers {
    fn or(&self, key: &str, body: &[u8]) -> Arc<Vec<u8>> {
        let mut cache = self.buffers.borrow_mut();
        cache
            .entry(key.to_string())
            .or_insert_with(|| Arc::new(body.to_vec()))
            .clone()
    }
}

/// Shared resources handed to the node-graph walk.
pub(super) struct GltfCtx<'a> {
    pub(super) meshes: &'a [GlMesh],
    pub(super) accessors: &'a [Accessor],
    /// The raw buffer bodies the parser read; used for base64 (`data:`) buffers
    /// and for the texture views. Primitive arrays that live in a file are
    /// re-read by the deferred loader.
    pub(super) buffers: &'a [Vec<u8>],
    /// The `uri` each buffer declared, empty for the GLB BIN body
    /// (`gltfParseBuffer`, `cxx:2287-2292`).
    pub(super) buffer_uris: &'a [String],
    /// Absolute offset of the GLB BIN body in the file — `myBinBodyOffset`
    /// (`cxx:2296`).
    pub(super) bin_offset: Option<usize>,
    /// The file being read — `myFilePath`, the `StreamUri` of a binary buffer
    /// (`cxx:2303`).
    pub(super) file_path: &'a str,
    pub(super) views: &'a [BufView],
    pub(super) materials: &'a [GlMat],
    /// `myCSTrsf` (`RWGltf_GltfJsonParser.hxx:490`): the glTF-to-OCCT converter.
    pub(super) converter: &'a RwMeshCoordinateSystemConverter,
    /// `myIsGltf1` (`RWGltf_GltfJsonParser.cxx:742`): enables the glTF 1.0 id
    /// name fallback.
    pub(super) is_gltf1: bool,
    /// `myUseMeshNameAsFallback` (`RWGltf_CafReader.cxx:175`, default true).
    pub(super) use_mesh_name_as_fallback: bool,
    /// `myShapeMap[ShapeMapGroup_Meshes]` / `[ShapeMapGroup_PrimArray]`: the
    /// mesh and primitive-array reuse caches (`RWGltf_GltfJsonParser.cxx:1733`,
    /// `:1847-1877`), shared by every scene root of one read.
    pub(super) cache: PrimCache,
    /// `myDecodedBuffers`.
    pub(super) decoded: DecodedBuffers,
}

/// `gltfParseAccessor` + `gltfParseBufferView` + `gltfParseBuffer`
/// (`RWGltf_GltfJsonParser.cxx:2039-2367`): validate one accessor and record
/// where its bytes are.
pub(super) fn prim_array_data(
    ctx: &GltfCtx<'_>,
    accessor_idx: usize,
    array_type: GltfArrayType,
) -> Result<GltfPrimArrayData, String> {
    let acc = ctx
        .accessors
        .get(accessor_idx)
        .ok_or_else(|| format!("gltf: accessor {accessor_idx} is out of range"))?;
    if acc.layout == GltfAccessorLayout::Unknown {
        return Err(format!("gltf: accessor '{accessor_idx}' has invalid type"));
    }
    let comp_type = GltfAccessorCompType::from_int(acc.comp);
    if !comp_type.is_valid() {
        return Err(format!(
            "gltf: accessor '{accessor_idx}' defines invalid componentType value"
        ));
    }
    if acc.count < 1 {
        return Err(format!("gltf: accessor '{accessor_idx}' defines invalid count"));
    }
    if acc.stride > 255 {
        return Err(format!("gltf: accessor '{accessor_idx}' defines invalid byteStride"));
    }
    let view_idx = acc
        .view
        .ok_or_else(|| format!("gltf: accessor '{accessor_idx}' does not define bufferView"))?;
    let view = ctx
        .views
        .get(view_idx)
        .ok_or_else(|| format!("gltf: accessor '{accessor_idx}' refers to non-existing bufferView"))?;
    if view.len == 0 {
        return Err(format!("gltf: bufferView '{view_idx}' defines invalid byteLength"));
    }
    if view.stride > 255 {
        return Err(format!("gltf: bufferView '{view_idx}' defines invalid byteStride"));
    }
    // `cxx:2229-2237`: only the two array targets are accepted.
    if view.target != 0 && view.target != 34962 && view.target != 34963 {
        return Err(format!("gltf: bufferView '{view_idx}' defines invalid target"));
    }
    let buffer_idx = view.buffer;
    let uri = ctx
        .buffer_uris
        .get(buffer_idx)
        .ok_or_else(|| format!("gltf: bufferView '{view_idx}' refers to non-existing buffer"))?;

    // `cxx:2292`: inside a GLB, a buffer without `uri` is the BIN body.
    let is_binary = ctx.bin_offset.is_some() && uri.is_empty();
    let offset = view.offset + acc.offset;
    let accessor = GltfAccessor {
        byte_offset: 0,
        count: acc.count,
        // `cxx:2300`, `:2317`: the bufferView stride (glTF 2.0) replaces the
        // accessor's; glTF 1.0 kept `byteStride` on the accessor
        // (`gltfParseAccessor`, cxx:2100-2112), so an empty bufferView stride
        // falls back to it.
        byte_stride: if view.stride != 0 { view.stride } else { acc.stride },
        layout: acc.layout,
        comp_type,
    };
    if is_binary {
        // `cxx:2296-2305`: the BIN body is re-read from the file, so the offset
        // carries the body's position in the GLB.
        return Ok(GltfPrimArrayData {
            stream_data: None,
            stream_uri: ctx.file_path.to_string(),
            stream_offset: ctx.bin_offset.unwrap_or(0) + offset,
            stream_length: view.len,
            accessor,
            array_type,
        });
    }
    // `cxx:2313`: only this exact prefix is a base64 buffer to OCCT; any other
    // `data:` URI is treated as a file name and fails the existence check.
    const BASE64_PREFIX: &str = "data:application/octet-stream;base64,";
    if let Some(rest) = uri.strip_prefix(BASE64_PREFIX) {
        let body = ctx
            .buffers
            .get(buffer_idx)
            .ok_or_else(|| format!("gltf: bufferView '{view_idx}' refers to non-existing buffer"))?;
        let _ = rest;
        // `cxx:2315-2331`: `StreamLength` stays 0; `StreamOffset` indexes the
        // decoded body.
        return Ok(GltfPrimArrayData {
            stream_data: Some(ctx.decoded.or(&buffer_idx.to_string(), body)),
            stream_uri: String::new(),
            stream_offset: offset,
            stream_length: 0,
            accessor,
            array_type,
        });
    }
    if uri.is_empty() {
        return Err("gltf: buffer does not define uri".to_string());
    }
    // `cxx:2343-2351`: a file buffer must exist; `StreamUri` is the resolved
    // path, read again by the deferred loader.
    let path = Path::new(ctx.file_path)
        .parent()
        .unwrap_or_else(|| Path::new("."))
        .join(uri);
    if !path.exists() {
        return Err(format!(
            "gltf: buffer refers to non-existing file '{uri}'"
        ));
    }
    Ok(GltfPrimArrayData {
        stream_data: None,
        stream_uri: path.to_string_lossy().to_string(),
        stream_offset: offset,
        stream_length: view.len,
        accessor,
        array_type,
    })
}

/// `gltfParsePrimArray` (`RWGltf_GltfJsonParser.cxx:1782-2035`) up to the shape:
/// collect the array's elements and its deferred node/triangle counts.
pub(super) fn build_late_array(
    ctx: &GltfCtx<'_>,
    mesh_idx: usize,
    prim_idx: usize,
    mesh_name: &str,
    prim: &Prim,
) -> Result<GltfLatePrimitiveArray, String> {
    let mesh_id = mesh_idx.to_string();
    let mut data: Vec<GltfPrimArrayData> = Vec::new();
    let mut has_positions = false;
    for (ty, acc_idx) in &prim.attribs {
        // `cxx:1889-1927`: attributes are walked in JSON order;
        // `RWGltf_GltfLatePrimitiveArray::AddPrimArrayData` keeps POSITION
        // first and places the index array right after it.
        let element = prim_array_data(ctx, *acc_idx, *ty)?;
        if *ty == GltfArrayType::Position {
            has_positions = true;
        }
        add_prim_array_data(&mut data, element);
    }
    if !has_positions {
        return Err(format!(
            "gltf: primitive array within mesh '{mesh_id}' does not define vertex positions"
        ));
    }
    if let Some(ii) = prim.indices {
        let element = prim_array_data(ctx, ii, GltfArrayType::Indices)?;
        add_prim_array_data(&mut data, element);
    }
    let mut late = GltfLatePrimitiveArray::new(
        mesh_id.clone(),
        mesh_name.to_string(),
        prim.mode,
        data,
    );
    // `cxx:2120-2126`, `:2173-2184`: the accessor counts give the deferred
    // estimates.
    let mut nb_nodes = 0usize;
    let mut nb_tris = 0usize;
    for element in late.data() {
        match element.array_type {
            GltfArrayType::Position => nb_nodes = element.accessor.count,
            GltfArrayType::Indices => {
                if prim.mode == RwGltfPrimitiveMode::Triangles {
                    nb_tris = element.accessor.count / 3;
                } else {
                    nb_nodes = element.accessor.count;
                    nb_tris = 0;
                }
            }
            _ => {}
        }
    }
    late.set_nb_deferred_nodes(nb_nodes);
    late.set_nb_deferred_triangles(nb_tris);
    if prim.indices.is_none() {
        // `cxx:1951-1958`: without an index accessor a triangle list is
        // reconstructed from consecutive nodes.
        if prim.mode == RwGltfPrimitiveMode::Triangles {
            let nb = late.nb_deferred_nodes() / 3;
            late.set_nb_deferred_triangles(nb);
        } else {
            late.set_nb_deferred_triangles(0);
        }
    }
    let _ = prim_idx;
    Ok(late)
}

/// `RWGltf_GltfLatePrimitiveArray::AddPrimArrayData` (`cxx:56-86`): POSITION
/// goes first, INDICES right after it, every other attribute is appended.
fn add_prim_array_data(data: &mut Vec<GltfPrimArrayData>, element: GltfPrimArrayData) {
    match element.array_type {
        GltfArrayType::Position => data.insert(0, element),
        GltfArrayType::Indices => {
            let after_position = data
                .first()
                .map(|d| d.array_type == GltfArrayType::Position)
                .unwrap_or(false);
            if after_position {
                data.insert(1, element);
            } else {
                data.insert(0, element);
            }
        }
        _ => data.push(element),
    }
}

/// The layout a glTF accessor `type` string names
/// (`RWGltf_GltfParseAccessorType`, `RWGltf_GltfAccessorLayout.hxx:33-67`).
pub(super) fn accessor_layout(s: &str) -> GltfAccessorLayout {
    GltfAccessorLayout::from_type_str(s)
}

/// `RWGltf_GltfParseAttribType` re-exported for the parser loop.
pub(super) fn attrib_type(name: &str) -> GltfArrayType {
    parse_attrib_type(name)
}

/// Build the local transform for a glTF node.
///
/// Source: `RWGltf_GltfJsonParser::parseTransformationMatrix`
/// (`RWGltf_GltfJsonParser.cxx:451-501`) and `parseTransformationComponents`
/// (`cxx:508-643`). The two forms are mutually exclusive (`cxx:1594-1600`:
/// "defines ambiguous transformation"), every component array is length and
/// type checked, and a zero scale component is an error (`cxx:585-588`).
///
/// `to_apply_scale` is `myToApplyScale` (`RWGltf_CafReader.cxx:180`, default
/// true). It only matters for a *non-uniform* scale (`cxx:590-591`: any two
/// components differ by more than `Precision::Confusion()`):
/// * true — OCCT hands the scale to the document fill, which multiplies the
///   primitive's triangulation nodes (`RWGltf_CafReader.cxx:501-519`) and keeps
///   the node location scale-free;
/// * false — the scale is folded into the location matrix and rebuilt through
///   `gp_Trsf::SetValues` (`cxx:596-616`), which orthogonalizes the linear part
///   (`gp_Trsf.cxx:862-938`), i.e. the anisotropy is lost.
///
/// A uniform scale is a similarity and goes straight into the location
/// (`cxx:631-634`: `aTrsf.SetScaleFactor`).
///
/// This reader emits world-space vertices per node, so "scale the
/// triangulation, keep the location scale-free" and "keep the scale in the
/// composed transform" are the same point set: the `true` arm folds the scale
/// into the node's linear part exactly, while the `false` arm reproduces
/// `SetValues`' orthogonalization.
pub(super) fn node_transform(nd: &JVal, id: &str, to_apply_scale: bool) -> Result<GpTrsf, String> {
    const CONFUSION: f64 = occt_core::precision::CONFUSION;
    const REAL_SMALL: f64 = occt_core::precision::REAL_SMALL;

    let mat_val = nd.get("matrix");
    let rot_val = nd.get("rotation");
    let scale_val = nd.get("scale");
    let tr_val = nd.get("translation");
    let has_components = rot_val.is_some() || scale_val.is_some() || tr_val.is_some();
    if has_components && mat_val.is_some() {
        return Err(format!("gltf: scene node '{id}' defines ambiguous transformation"));
    }

    if let Some(mv) = mat_val {
        // `cxx:455-476`: exactly 16 numbers, column-major.
        let arr = mv
            .as_arr()
            .filter(|a| a.len() == 16)
            .ok_or_else(|| format!("gltf: scene node '{id}' defines invalid transformation matrix array"))?;
        let mut m = [0.0f64; 16];
        for (i, v) in arr.iter().enumerate() {
            m[i] = v
                .as_num()
                .ok_or_else(|| format!("gltf: scene node '{id}' defines invalid transformation matrix"))?;
        }
        return trsf_from_matrix(&m).map_err(|e| format!("gltf: scene node '{id}': {e}"));
    }
    if !has_components {
        return Ok(GpTrsf::identity());
    }

    let mut trsf = GpTrsf::identity();
    if let Some(rv) = rot_val {
        // `cxx:516-543`: four numbers, identity quaternion leaves the transform alone.
        let arr = rv
            .as_arr()
            .filter(|a| a.len() == 4)
            .ok_or_else(|| format!("gltf: scene node '{id}' defines invalid rotation quaternion"))?;
        let mut r = [0.0f64; 4];
        for (i, v) in arr.iter().enumerate() {
            r[i] = v
                .as_num()
                .ok_or_else(|| format!("gltf: scene node '{id}' defines invalid rotation"))?;
        }
        if r[0].abs() > REAL_SMALL
            || r[1].abs() > REAL_SMALL
            || r[2].abs() > REAL_SMALL
            || (r[3] - 1.0).abs() > REAL_SMALL
        {
            trsf.set_rotation_quat(&GpQuaternion::new(r[0], r[1], r[2], r[3]));
        }
    }
    if let Some(tv) = tr_val {
        // `cxx:545-565`.
        let arr = tv
            .as_arr()
            .filter(|a| a.len() == 3)
            .ok_or_else(|| format!("gltf: scene node '{id}' defines invalid translation vector"))?;
        let mut t = [0.0f64; 3];
        for (i, v) in arr.iter().enumerate() {
            t[i] = v
                .as_num()
                .ok_or_else(|| format!("gltf: scene node '{id}' defines invalid translation"))?;
        }
        trsf.set_translation_part(&occt_core::gp::GpVec::new(t[0], t[1], t[2]));
    }
    if let Some(sv) = scale_val {
        // `cxx:567-634`.
        let arr = sv
            .as_arr()
            .filter(|a| a.len() == 3)
            .ok_or_else(|| format!("gltf: scene node '{id}' defines invalid scale vector"))?;
        let mut s = [0.0f64; 3];
        for (i, v) in arr.iter().enumerate() {
            s[i] = v
                .as_num()
                .ok_or_else(|| format!("gltf: scene node '{id}' defines invalid scale"))?;
            if s[i].abs() <= REAL_SMALL {
                return Err(format!("gltf: scene node '{id}' defines invalid scale"));
            }
        }
        let non_uniform = (s[0] - s[1]).abs() > CONFUSION
            || (s[1] - s[2]).abs() > CONFUSION
            || (s[0] - s[2]).abs() > CONFUSION;
        if non_uniform {
            if to_apply_scale {
                // `cxx:627-628` (`theHasScale`/`theScale`): the scale is applied
                // to the triangulation nodes; here it is folded into the node's
                // own linear part, giving the same world-space vertices.
                let mut m = trsf.h_vectorial_part();
                for i in 0..3 {
                    for j in 0..3 {
                        m.m[i][j] *= s[j];
                    }
                }
                trsf.matrix = m;
                trsf.shape = occt_core::gp::TrsfForm::CompoundTrsf;
            } else {
                // `cxx:599-616`: `aMat4 = aMat4 * aScaleMat` then
                // `gp_Trsf::SetValues(...)`.
                let mut m = trsf.h_vectorial_part();
                for i in 0..3 {
                    for j in 0..3 {
                        m.m[i][j] *= s[j];
                    }
                }
                let l = trsf.translation_part();
                let mut folded = GpTrsf::identity();
                folded
                    .set_values(
                        m.m[0][0], m.m[0][1], m.m[0][2], l.x(),
                        m.m[1][0], m.m[1][1], m.m[1][2], l.y(),
                        m.m[2][0], m.m[2][1], m.m[2][2], l.z(),
                    )
                    .map_err(|e| format!("gltf: scene node '{id}': {e}"))?;
                trsf = folded;
            }
        } else if (s[0] - 1.0).abs() > CONFUSION {
            // `cxx:631-634`.
            trsf.set_scale_factor(s[0])
                .map_err(|e| format!("gltf: scene node '{id}': {e}"))?;
        }
    }
    Ok(trsf)
}
