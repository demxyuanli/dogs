//! glTF file reader: buffers, accessors, materials, mesh/primitives, scene graph.
//!
//! Source: RWGltf_CafReader::performMesh (RWGltf_CafReader.cxx:186-410) and
//! RWGltf_GltfJsonParser (gltfParseScene cxx:1456-1516,
//! gltfParseSceneNodes cxx:1520-1562, gltfParseSceneNode
//! cxx:1566-1719). The container split lives in [`container`], the glTF 1.0
//! keyed-collection rewrite in [`gltf1_map`].

use super::gltf_late_data::{read_late_data, LateDataOptions};
use super::gltf_types::{
    node_transform, Accessor, BufView, DecodedBuffers, GlMat, GlMesh, GlNode, GltfCtx, Prim,
};
use super::prelude::*;
use super::*;
use crate::model::BRepModel;

mod container;
mod gltf1_map;

use container::*;
use gltf1_map::*;

/// Conversion to apply while reading (the port's stand-in for
/// `RWGltf_CafReader::SetCoordinateSystemConverter`,
/// `RWGltf_CafReader.hxx`).
#[derive(Debug, Clone, Default)]
pub struct GltfReadOptions {
    /// `myCoordSysConverter`. The default mirrors
    /// `RWGltf_CafReader::RWGltf_CafReader` (`RWGltf_CafReader.cxx:182-183`):
    /// input unit 1.0 (glTF is metres) and input coordinate system glTF
    /// (Y-up). The output side stays undefined, so the converter is empty and
    /// the file's coordinates are passed through unchanged — exactly the OCCT
    /// default.
    pub converter: RwMeshCoordinateSystemConverter,
    /// `myToSkipEmptyNodes` (`RWGltf_CafReader.cxx:174`, default true): a root
    /// node whose shape holds no face, edge or vertex is not added to the
    /// document (`RWGltf_GltfJsonParser.cxx:1554-1559`).
    ///
    /// In this port a scene node is only emitted from a primitive array, i.e.
    /// from a shape: a node whose mesh yields no array has no shape to add, so
    /// it never appears whether the flag is set or not. The flag is kept for
    /// reader-option parity.
    pub to_skip_empty_nodes: bool,
    /// `myToLoadAllScenes` (`RWGltf_CafReader.cxx:174`, default false): load
    /// every scene of a multi-scene glTF 2.0 file instead of the default one
    /// (`RWGltf_GltfJsonParser.cxx:1460-1481`).
    pub to_load_all_scenes: bool,
    /// `myUseMeshNameAsFallback` (`RWGltf_CafReader.cxx:175`, default true): an
    /// unnamed node takes the mesh name (`RWGltf_GltfJsonParser.cxx:2444-2487`).
    pub use_mesh_name_as_fallback: bool,
    /// `myToApplyScale` (`RWGltf_CafReader.cxx:180`, default true): a
    /// non-uniform node scale is applied to the primitive's triangulation
    /// directly instead of being folded into the node location
    /// (`RWGltf_GltfJsonParser.cxx:596-628`).
    pub to_apply_scale: bool,
    /// `myToSkipLateDataLoading` (`RWGltf_CafReader.cxx:174`, default false):
    /// leave the file-backed primitive arrays as proxies whose triangulation is
    /// loaded later (`RWGltf_CafReader::readLateData`, cxx:367-404).
    pub to_skip_late_data_loading: bool,
    /// `myToKeepLateData` (`RWGltf_CafReader.cxx:178`, default true): fill the
    /// source object in place instead of replacing the face triangulation
    /// (`BRep_Builder::UpdateFace`).
    pub to_keep_late_data: bool,
    /// `myToParallel` (`RWGltf_CafReader.cxx:173`, default false): decode the
    /// deferred triangulations on the default thread pool.
    pub to_parallel: bool,
}

impl GltfReadOptions {
    /// A reader configured like `RWGltf_CafReader`'s constructor
    /// (`cxx:171-184`).
    pub fn occt_defaults() -> Self {
        let mut converter = RwMeshCoordinateSystemConverter::new();
        converter.set_input_length_unit(1.0);
        converter.set_input_coordinate_system(RwMeshCoordinateSystem::GLTF);
        Self {
            converter,
            to_skip_empty_nodes: true,
            to_load_all_scenes: false,
            use_mesh_name_as_fallback: true,
            to_apply_scale: true,
            to_skip_late_data_loading: false,
            to_keep_late_data: true,
            to_parallel: false,
        }
    }
}

/// Read a glTF 2.0 `.gltf` or `.glb` file into a scene with the OCCT reader's
/// default options (`RWGltf_CafReader::RWGltf_CafReader`,
/// `RWGltf_CafReader.cxx:171-184`).
pub fn read_gltf_scene(path: &str) -> Result<MeshScene, String> {
    read_gltf_scene_with_options(path, &GltfReadOptions::occt_defaults())
}

/// Read a glTF file with an explicit coordinate-system converter.
///
/// Returns one node per mesh node, with the node's world transform (TRS or
/// matrix, accumulated through the scene graph) applied to its vertices.
/// External `.bin` buffers resolve relative to the glTF file, embedded `data:`
/// URIs are base64-decoded, and a buffer without `uri` inside a `.glb` is the
/// BIN chunk (`RWGltf_GltfJsonParser.cxx:2286-2306`). When the converter is not
/// empty its `TransformPosition` is applied to every vertex and
/// `TransformTransformation` to the node location, like
/// `RWGltf_TriangulationReader.cxx:325` and `RWGltf_GltfJsonParser.cxx:493`.
pub fn read_gltf_scene_with_options(
    path: &str,
    opts: &GltfReadOptions,
) -> Result<MeshScene, String> {
    read_gltf_impl(path, opts, None)
}

/// `RWGltf_CafReader::performMesh` followed by `fillDocument`
/// (`RWGltf_CafReader.cxx:186-410`, `:435-467`): the same read, with the
/// `XCAFDoc_ShapeTool` label tree filled into `model` on the way out.
pub(super) fn read_gltf_document_with_options(
    path: &str,
    opts: &GltfReadOptions,
) -> Result<BRepModel, String> {
    let mut model = BRepModel::new();
    read_gltf_impl(path, opts, Some(&mut model))?;
    Ok(model)
}

fn read_gltf_impl(
    path: &str,
    opts: &GltfReadOptions,
    doc: Option<&mut BRepModel>,
) -> Result<MeshScene, String> {
    let bytes = std::fs::read(path)
        .map_err(|e| format!("rwmesh: read {path}: {e}"))?;
    let container = gltf_container(&bytes)?;
    let text = container.json;
    let mut bin_body = container.bin;
    let bin_offset = container.bin_offset;
    let has_bin_body = bin_body.is_some();
    let root = parse_json(&text)?;
    let is_gltf1 = root.get("asset")
        .and_then(|a| a.get("version"))
        .and_then(|v| v.as_str())
        .map(|v| v.starts_with("1."))
        .unwrap_or(false);
    let mut root = root;
    // Keyed (glTF 1.0) collections and the glTF 1.0 id fallback names
    // (`RWGltf_GltfJsonParser.cxx:2410-2418`).
    let (node_ids, mesh_ids) = normalize_gltf1(&mut root);

    // buffers → raw bytes
    //
    // `RWGltf_GltfJsonParser::gltfParseBuffer` (`cxx:2286-2306`): inside a GLB a
    // buffer without `uri` (glTF 2.0) or named `binary_glTF` (glTF 1.0) is the
    // BIN chunk; everything else is a base64 `data:` URI or an external file.
    let mut buffers: Vec<Vec<u8>> = Vec::new();
    // Parallel to `buffers`: the `uri` each buffer declared (empty for the GLB
    // BIN chunk), needed to classify texture bufferViews.
    let mut buffer_uris: Vec<String> = Vec::new();
    if let Some(arr) = root.get("buffers").and_then(|v| v.as_arr()) {
        for (bi, buf) in arr.iter().enumerate() {
            let uri = buf.get("uri").and_then(|v| v.as_str()).unwrap_or("").to_string();
            let buf_name = buf.get("name").and_then(|v| v.as_str()).unwrap_or("");
            let want = buf.get("byteLength").and_then(|v| v.as_num()).map(|n| n as usize).unwrap_or(0);
            let is_binary = bin_body.is_some()
                && (uri.is_empty() || buf_name == "binary_glTF");
            let mut bytes = if is_binary {
                bin_body.take().expect("checked above")
            } else if let Some(data) = uri.strip_prefix("data:") {
                let b64 = data.split("base64,").nth(1)
                    .ok_or_else(|| format!("gltf: buffer {bi}: malformed data uri"))?;
                base64_decode(b64)?
            } else if uri.is_empty() {
                return Err(format!("gltf: buffer {bi} does not define uri"));
            } else {
                let dir = Path::new(path).parent().unwrap_or_else(|| Path::new("."));
                std::fs::read(dir.join(&uri))
                    .map_err(|e| format!("gltf: buffer {bi}: read {}: {e}", dir.join(&uri).display()))?
            };
            if want > 0 && bytes.len() < want {
                return Err(format!("gltf: buffer {bi}: expected {want} bytes, got {}", bytes.len()));
            }
            if want > 0 && bytes.len() > want {
                bytes.truncate(want);
            }
            buffers.push(bytes);
            buffer_uris.push(uri);
        }
    }
    if buffers.is_empty() {
        return Err("gltf: no buffers".into());
    }

    // bufferViews
    let mut views: Vec<BufView> = Vec::new();
    if let Some(arr) = root.get("bufferViews").and_then(|v| v.as_arr()) {
        for bv in arr {
            views.push(BufView {
                buffer: bv.get("buffer").and_then(|v| v.as_num()).map(|n| n as usize).unwrap_or(0),
                offset: bv.get("byteOffset").and_then(|v| v.as_num()).map(|n| n as usize).unwrap_or(0),
                len: bv.get("byteLength").and_then(|v| v.as_num()).map(|n| n as usize).unwrap_or(0),
                // `gltfParseBufferView` (`cxx:2227-2237`): `byteStride` and
                // `target` live here in glTF 2.0; the stride is what the
                // triangulation reader uses (`cxx:2300`).
                stride: bv
                    .get("byteStride")
                    .and_then(|v| v.as_int())
                    .map(|n| n as usize)
                    .unwrap_or(0),
                target: bv.get("target").and_then(|v| v.as_int()).unwrap_or(0),
            });
        }
    }

    // accessors
    let mut accessors: Vec<Accessor> = Vec::new();
    if let Some(arr) = root.get("accessors").and_then(|v| v.as_arr()) {
        for ac in arr {
            accessors.push(Accessor {
                view: ac.get("bufferView").and_then(|v| v.as_num()).map(|n| n as usize),
                offset: ac.get("byteOffset").and_then(|v| v.as_num()).map(|n| n as usize).unwrap_or(0),
                // `cxx:2100`: `ByteStride` is an integer, 0 when absent and
                // rejected above 255; glTF 1.0 kept it on the accessor.
                stride: ac
                    .get("byteStride")
                    .and_then(|v| v.as_int())
                    .map(|n| n as usize)
                    .unwrap_or(0),
                comp: ac.get("componentType").and_then(|v| v.as_int()).unwrap_or(0),
                count: ac.get("count").and_then(|v| v.as_num()).map(|n| n as usize).unwrap_or(0),
                // `cxx:2057-2068` / `RWGltf_GltfParseAccessorType`: the layout
                // decides whether an attribute's component count is acceptable.
                layout: ac
                    .get("type")
                    .and_then(|v| v.as_str())
                    .map(gltf_types::accessor_layout)
                    .unwrap_or(GltfAccessorLayout::Unknown),
            });
        }
    }

    // images / textures: `gltfParseTexture` (`cxx:1164-1311`) reads an image
    // either from a URI (file or `data:`) or from a bufferView slice.
    let mut images: Vec<GltfImageDef> = Vec::new();
    if let Some(obj) = root.get("images").and_then(|v| v.as_obj()) {
        // glTF 1.0 stores images as an object keyed by id.
        for (id, im) in obj {
            images.push(GltfImageDef {
                id: id.clone(),
                uri: im.get("uri").and_then(|v| v.as_str()).map(str::to_string),
                buffer_view: im.get("bufferView").and_then(|v| v.as_num()).map(|n| n as usize),
                khr_binary_buffer_view: im
                    .get("extensions")
                    .and_then(|e| e.get("KHR_binary_glTF"))
                    .and_then(|b| b.get("bufferView"))
                    .and_then(|v| v.as_num())
                    .map(|n| n as usize),
            });
        }
    } else if let Some(arr) = root.get("images").and_then(|v| v.as_arr()) {
        for (i, im) in arr.iter().enumerate() {
            images.push(GltfImageDef {
                id: i.to_string(),
                uri: im.get("uri").and_then(|v| v.as_str()).map(str::to_string),
                buffer_view: im.get("bufferView").and_then(|v| v.as_num()).map(|n| n as usize),
                khr_binary_buffer_view: None,
            });
        }
    }
    let mut textures: Vec<GltfTextureRef> = Vec::new();
    if let Some(obj) = root.get("textures").and_then(|v| v.as_obj()) {
        for (id, tx) in obj {
            textures.push(GltfTextureRef {
                id: id.clone(),
                source: tx.get("source").and_then(GltfKey::of),
                target: tx.get("target").and_then(|v| v.as_num()),
            });
        }
    } else if let Some(arr) = root.get("textures").and_then(|v| v.as_arr()) {
        for (i, tx) in arr.iter().enumerate() {
            textures.push(GltfTextureRef {
                id: i.to_string(),
                source: tx.get("source").and_then(GltfKey::of),
                target: tx.get("target").and_then(|v| v.as_num()),
            });
        }
    }

    let base_dir = Path::new(path).parent().unwrap_or_else(|| Path::new("."));
    let tex_src = GltfTextureSource {
        textures: &textures,
        images: &images,
        views: &views,
        buffers: &buffers,
        buffer_uris: &buffer_uris,
        base_dir,
        file_path: path,
        bin_offset: (bin_offset > 0).then_some(bin_offset),
        is_gltf1: root.get("asset").and_then(|a| a.get("version"))
            .and_then(|v| v.as_str())
            .map(|s| s.starts_with('1'))
            .unwrap_or(false),
    };

    // materials: `gltfParseMaterials` (`cxx:793-863`) parses a PBR and a common
    // half per entry, then `gltfBindMaterial` (`cxx:867-947`) merges both into
    // one `XCAFDoc_VisMaterial`.
    let mut gl_materials: Vec<GlMat> = Vec::new();
    let bind = |id: String, ma: &JVal, gl_materials: &mut Vec<GlMat>| {
        let name = ma.get("name").and_then(|v| v.as_str()).unwrap_or("").to_string();
        let mut pbr = parse_pbr_material(ma, &tex_src);
        let mut common = parse_common_material(ma, &tex_src)
            .or_else(|| parse_std_material(ma, &tex_src));
        if let Some(p) = pbr.as_mut() {
            p.id = id.clone();
            p.name = name.clone();
        }
        if let Some(c) = common.as_mut() {
            c.id = id.clone();
            c.name = name.clone();
        }
        if let Some(vis) = bind_material(pbr, common) {
            gl_materials.push(GlMat { id, name, vis });
        }
    };
    if let Some(obj) = root.get("materials").and_then(|v| v.as_obj()) {
        for (id, ma) in obj {
            bind(id.clone(), ma, &mut gl_materials);
        }
    } else if let Some(arr) = root.get("materials").and_then(|v| v.as_arr()) {
        for (i, ma) in arr.iter().enumerate() {
            bind(format!("mat_{i}"), ma, &mut gl_materials);
        }
    }

    // meshes → primitives
    let mut meshes: Vec<GlMesh> = Vec::new();
    if let Some(arr) = root.get("meshes").and_then(|v| v.as_arr()) {
        for (mi, ms) in arr.iter().enumerate() {
            let mut prims = Vec::new();
            if let Some(p_arr) = ms.get("primitives").and_then(|v| v.as_arr()) {
                for pr in p_arr {
                    // `RWGltf_GltfJsonParser.cxx:1797-1822`: `mode` defaults to
                    // triangles, must be an integer in 0..=6, and only points,
                    // lines and triangles produce a shape — the strip/loop/fan
                    // modes are skipped.
                    let mode = match pr.get("mode") {
                        None => RwGltfPrimitiveMode::Triangles,
                        Some(v) => v
                            .as_num()
                            .map(|n| RwGltfPrimitiveMode::from_int(n as i64))
                            .unwrap_or(RwGltfPrimitiveMode::Unknown),
                    };
                    if mode == RwGltfPrimitiveMode::Unknown {
                        return Err(format!(
                            "gltf: primitive array within mesh {mi} has unknown mode"
                        ));
                    }
                    if !mode.is_supported() {
                        continue;
                    }
                    let mut p = Prim::default();
                    p.mode = mode;
                    // `gltfParsePrimArray` (`cxx:1795-1803`): the primitive must
                    // define an `attributes` object.
                    let attrs = pr
                        .get("attributes")
                        .and_then(|v| v.as_obj())
                        .ok_or_else(|| {
                            format!("gltf: primitive array within mesh {mi} defines no attributes")
                        })?;
                    for (k, val) in attrs {
                        // `cxx:1899-1907`: an attribute name with no glTF array
                        // type is ignored; the accessor key must resolve.
                        let ty = gltf_types::attrib_type(k);
                        let ai = val
                            .as_num()
                            .map(|n| n as usize)
                            .ok_or_else(|| format!("gltf: attribute '{k}' accessor key is not a number"))?;
                        // `gltfParseAccessor` (`cxx:1909-1919`): the accessor
                        // must exist, whatever its array type.
                        if accessors.get(ai).is_none() {
                            return Err(format!(
                                "gltf: mesh {mi}: attribute '{k}' accessor {ai} out of range"
                            ));
                        }
                        if ty == GltfArrayType::Unknown {
                            continue;
                        }
                        if ty == GltfArrayType::Position {
                            p.pos_acc = Some(ai);
                        }
                        p.attribs.push((ty, ai));
                    }
                    if let Some(ii) = pr.get("indices") {
                        let ai = ii
                            .as_num()
                            .map(|n| n as usize)
                            .ok_or_else(|| format!("gltf: mesh {mi}: indices key is not a number"))?;
                        // `cxx:1935-1948`: the indices accessor must exist.
                        if accessors.get(ai).is_none() {
                            return Err(format!(
                                "gltf: mesh {mi}: indices accessor {ai} out of range"
                            ));
                        }
                        p.indices = Some(ai);
                    }
                    // glTF 2.0 references a material by index; glTF 1.0 by id
                    // (`RWGltf_GltfJsonParser.cxx:1823-1844` resolves the key
                    // through its material maps).
                    p.mat_idx = match pr.get("material") {
                        Some(JVal::Str(s)) => {
                            gl_materials.iter().position(|m| &m.id == s || &m.name == s)
                        }
                        Some(v) => v.as_num().map(|n| n as usize),
                        None => None,
                    };
                    prims.push(p);
                }
            }
            meshes.push(GlMesh {
                id: mesh_ids.get(mi).cloned().unwrap_or_else(|| mi.to_string()),
                name: ms.get("name").and_then(|v| v.as_str()).unwrap_or("").to_string(),
                // `cxx:1726`, `:1772-1774`: the mesh's `extras`.
                extras: parse_extras(ms.get("extras")),
                prims,
            });
        }
    }

    // nodes
    let mut nodes: Vec<GlNode> = Vec::new();
    if let Some(arr) = root.get("nodes").and_then(|v| v.as_arr()) {
        for (ni, nd) in arr.iter().enumerate() {
            let id = node_ids.get(ni).cloned().unwrap_or_else(|| ni.to_string());
            nodes.push(GlNode {
                name: nd.get("name").and_then(|v| v.as_str()).unwrap_or("").to_string(),
                id: id.clone(),
                mesh: nd.get("mesh").and_then(|v| v.as_num()).map(|n| n as usize),
                meshes: nd.get("meshes").and_then(|v| v.as_arr())
                    .map(|a| a.iter().filter_map(|c| c.as_num().map(|n| n as usize)).collect())
                    .unwrap_or_default(),
                trsf: node_transform(nd, &id, opts.to_apply_scale)?,
                children: nd.get("children").and_then(|v| v.as_arr())
                    .map(|a| a.iter().filter_map(|c| c.as_num().map(|n| n as usize)).collect())
                    .unwrap_or_default(),
                // `cxx:1621-1622`: the node's `extras`.
                extras: parse_extras(nd.get("extras")),
            });
        }
    }
    if nodes.is_empty() {
        return Err("gltf: no nodes".into());
    }

    // scene roots (`RWGltf_GltfJsonParser::gltfParseScene`, `cxx:1456-1516`)
    let scene_roots = gltf_scene_roots(&root, is_gltf1, opts.to_load_all_scenes)?;

    let mut scene_nodes = Vec::new();
    let gltf_ctx = GltfCtx {
        meshes: &meshes,
        accessors: &accessors,
        buffers: &buffers,
        buffer_uris: &buffer_uris,
        bin_offset: if has_bin_body { Some(bin_offset) } else { None },
        file_path: path,
        views: &views,
        materials: &gl_materials,
        converter: &opts.converter,
        is_gltf1,
        use_mesh_name_as_fallback: opts.use_mesh_name_as_fallback,
        cache: PrimCache::default(),
        decoded: DecodedBuffers::default(),
    };
    // `myShapeMap[ShapeMapGroup_Nodes]`: `findNodeShape` caches a node shape by
    // its id across every scene of one read.
    let mut visited = vec![false; nodes.len()];
    let mut insts: Vec<GltfNodeInst> = Vec::new();
    for roots in &scene_roots {
        for &r in roots {
            // `gltfParseSceneNodes` (`cxx:1536-1542`): a reference to a node
            // that does not exist is reported and stops the remaining roots of
            // that scene.
            if r >= nodes.len() {
                break;
            }
            collect_gltf_nodes(&nodes, r, &GpTrsf::identity(), &gltf_ctx, &mut visited, &mut insts)?;
        }
    }

    // `RWGltf_CafReader::readLateData` (`cxx:367-404`): every proxy face gets
    // its triangulation before the document is filled, because the fill needs
    // the nodes to build the planar wire and the `Poly_Triangulation`.
    let faces = gltf_ctx.cache.faces();
    read_late_data(
        &faces,
        &opts.converter,
        LateDataOptions {
            to_skip_late_data_loading: opts.to_skip_late_data_loading,
            to_keep_late_data: opts.to_keep_late_data,
            to_parallel: opts.to_parallel,
        },
    )?;

    // `gltfParseSceneNodes` (`RWGltf_GltfJsonParser.cxx:1554-1559`):
    // `myToSkipEmptyNodes` drops a root whose `findNodeShape` returned a nil
    // shape. A node gets a shape as soon as `gltfParsePrimArray` built a
    // primitive array for it, and that happens at parse time regardless of the
    // deferred triangulation, so a nil node shape means "the node carries no
    // mesh". Such a node produces no instance above, which is the same end
    // state; a proxy whose data is not loaded is *not* empty here.
    scene_nodes.extend(insts.iter().map(inst_to_node));

    // `RWGltf_CafReader::fillDocument` (`cxx:435-467`): replay the parsed shape
    // tree and the per-shape attributes into the document labels.
    if let Some(model) = doc {
        super::gltf_document::fill_document(
            model,
            &nodes,
            &meshes,
            &scene_roots,
            &insts,
            &gltf_ctx,
            opts.to_skip_empty_nodes,
        );
    }

    let materials = gl_materials.iter().map(|m| material_view(&m.vis)).collect();
    Ok(MeshScene { name: file_stem(path).to_string(), nodes: scene_nodes, materials })
}


// ---------------------------------------------------------------------------
// Scene selection
// ---------------------------------------------------------------------------

/// `RWGltf_GltfJsonParser::gltfParseScene` (`cxx:1456-1516`): the root node
/// lists to walk — one list per loaded scene.
///
/// * `myToLoadAllScenes` with a glTF 2.0 file of more than one scene loads every
///   scene (`cxx:1460-1481`);
/// * otherwise the default scene (`scene`) is loaded; a missing `scene` falls
///   back to the first scene in 2.0 (`cxx:1484-1500`) but stays an error in 1.0,
///   whose `scene` is a key;
/// * a scene without a `nodes` array is an error (`cxx:1506-1514`).
fn gltf_scene_roots(
    root: &JVal,
    is_gltf1: bool,
    to_load_all_scenes: bool,
) -> Result<Vec<Vec<usize>>, String> {
    let scenes = root.get("scenes").and_then(|v| v.as_arr());
    let scene_count = scenes.map(|a| a.len()).unwrap_or(0);
    if to_load_all_scenes && !is_gltf1 && scene_count > 1 {
        let mut out = Vec::with_capacity(scene_count);
        for sc in scenes.unwrap_or(&[]) {
            out.push(scene_node_ids(sc)?);
        }
        return Ok(out);
    }

    let selected = match root.get("scene") {
        // `scenes` may still be a keyed object in 1.0; `normalize_gltf1` turns
        // it (and the `scene` key) into array form before this point.
        Some(JVal::Str(s)) => root
            .get("scenes")
            .and_then(|v| v.as_obj())
            .and_then(|o| o.iter().find(|(k, _)| k == s).map(|(_, v)| v)),
        Some(v) => v
            .as_num()
            .and_then(|n| scenes.and_then(|a| a.get(n as usize))),
        None => None,
    };
    let selected = match selected {
        Some(sc) => sc,
        None if root.get("scene").is_some() => {
            return Err("gltf: default scene is not found".into())
        }
        None if !is_gltf1 => {
            // `cxx:1493-1499`: no `scene` — the first one is loaded.
            scenes.and_then(|a| a.first())
                .ok_or_else(|| "gltf: default scene is not found".to_string())?
        }
        None => return Err("gltf: default scene is not found".into()),
    };
    Ok(vec![scene_node_ids(selected)?])
}

/// `findObjectMember(scene, "nodes")` + `gltfParseSceneNodes` input check
/// (`cxx:1506-1514`): the root node list of one scene.
fn scene_node_ids(scene: &JVal) -> Result<Vec<usize>, String> {
    scene
        .get("nodes")
        .and_then(|v| v.as_arr())
        .map(|a| a.iter().filter_map(|c| c.as_num().map(|n| n as usize)).collect())
        .ok_or_else(|| "gltf: empty scene".to_string())
}
