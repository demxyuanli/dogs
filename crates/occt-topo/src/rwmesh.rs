//! RWMesh — mesh format import (OBJ/PLY/STL/glTF/VRML) and scene assembly.
//! Source: `RWMesh`.
//!
//! Reads triangle meshes from Wavefront OBJ (with MTL materials and UVs), PLY,
//! STL (binary/ASCII auto-detected), glTF 2.0 (external or embedded base64
//! buffers, node TRS/matrix transforms, PBR base-color materials) and VRML 2.0
//! (Transform nesting, IndexedFaceSet, diffuse materials), assembles them into a
//! `MeshScene`, and converts scenes back to B-Rep shapes via `mesh_to_brep`. A
//! tiny format router (`convert_mesh_format`) round-trips between the formats by
//! extension.

use std::collections::HashMap;
use std::path::Path;

use occt_core::gp::{GpMat, GpPnt, GpPnt2d, GpQuaternion, GpTrsf, GpXyz, TrsfForm};
use occt_core::io::ply::PlyMesh;
use occt_core::io::stl::StlMesh;
use occt_core::poly::Triangulation;
use occt_core::poly::triangulation::Triangle;

use crate::builder::TopoBuilder;
use crate::gltf::GltfOptions;
use crate::mesh_to_brep::triangulation_to_brep;
use crate::shape::{Compound, TopoShape};

/// A surface material. `diffuse`/`specular` are RGB triples in `[0,1]`;
/// `opacity` is `1.0` for fully opaque.
#[derive(Debug, Clone)]
pub struct Material {
    pub name: String,
    pub diffuse: (f64, f64, f64),
    pub specular: (f64, f64, f64),
    pub opacity: f64,
}

impl Default for Material {
    fn default() -> Self {
        Material {
            name: "default".to_string(),
            diffuse: (0.7, 0.7, 0.7),
            specular: (1.0, 1.0, 1.0),
            opacity: 1.0,
        }
    }
}

/// A triangulated mesh scene: a named collection of mesh nodes plus a shared
/// named material library (populated by OBJ `mtllib` / glTF `materials`).
#[derive(Debug, Clone)]
pub struct MeshScene {
    pub name: String,
    pub nodes: Vec<MeshNode>,
    pub materials: Vec<Material>,
}

/// One triangulated mesh within a scene. `triangles` are indices into `vertices`;
/// `transform` (when present) is the node's world transform that produced the
/// vertex positions. `material`/`uv`/`texture_path` are optional surface
/// attributes: `uv` holds one texture coordinate per vertex (same length as
/// `vertices`), `texture_path` is an image file reference for the material.
#[derive(Debug, Clone)]
pub struct MeshNode {
    pub name: String,
    pub vertices: Vec<GpPnt>,
    pub triangles: Vec<(usize, usize, usize)>,
    pub transform: Option<GpTrsf>,
    pub material: Option<Material>,
    pub uv: Option<Vec<GpPnt2d>>,
    pub texture_path: Option<String>,
}

// ---------------------------------------------------------------------------
// Format readers
// ---------------------------------------------------------------------------

/// Read a Wavefront OBJ file into a scene. `o`/`g` groups become separate nodes;
/// a file without groups yields a single node. Wavefront `.mtl` libraries
/// referenced by `mtllib` are loaded (diffuse/specular/opacity), `usemtl`
/// assigns them to groups, and `vt` + `f v/vt` texture coordinates are carried
/// into `node.uv`.
pub fn read_obj_scene(path: &str) -> Result<MeshScene, String> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| format!("rwmesh: read {path}: {e}"))?;
    let dir = Path::new(path).parent().map(|p| p.to_path_buf());
    parse_obj_scene_inner(&content, file_stem(path), dir.as_deref())
}

/// Read a PLY file (ASCII) into a scene with a single node.
pub fn read_ply_scene(path: &str) -> Result<MeshScene, String> {
    let mesh = occt_core::io::ply::read_ply_file(path)
        .map_err(|e| format!("rwmesh: ply: {e}"))?;
    let tri = occt_core::io::ply::to_triangulation(&mesh);
    let triangles = tri.triangles.iter().map(|t| (t.n0, t.n1, t.n2)).collect();
    Ok(MeshScene {
        name: file_stem(path).to_string(),
        nodes: vec![MeshNode {
            name: "mesh".to_string(),
            vertices: tri.nodes,
            triangles,
            transform: None,
            material: None,
            uv: None,
            texture_path: None,
        }],
        materials: Vec::new(),
    })
}

/// Read an STL file (binary or ASCII, auto-detected) into a scene. Shared
/// vertices are deduplicated by coordinate (1e-9 hash grid).
pub fn read_stl_scene(path: &str) -> Result<MeshScene, String> {
    let mesh = occt_core::io::stl::read_stl_file(path)
        .map_err(|e| format!("rwmesh: stl: {e}"))?;
    let tri = occt_core::io::stl::to_triangulation(&mesh);
    let triangles = tri.triangles.iter().map(|t| (t.n0, t.n1, t.n2)).collect();
    Ok(MeshScene {
        name: file_stem(path).to_string(),
        nodes: vec![MeshNode {
            name: "mesh".to_string(),
            vertices: tri.nodes,
            triangles,
            transform: None,
            material: None,
            uv: None,
            texture_path: None,
        }],
        materials: Vec::new(),
    })
}

/// Read a glTF 2.0 `.gltf` file into a scene. Returns one node per mesh node,
/// with the node's world transform (TRS or matrix, accumulated through the scene
/// graph) applied to its vertices. External `.bin` buffers are resolved relative
/// to the `.gltf` file; embedded `data:` URIs are base64-decoded.
pub fn read_gltf_scene(path: &str) -> Result<MeshScene, String> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("rwmesh: read {path}: {e}"))?;
    let root = parse_json(&text)?;

    // buffers → raw bytes
    let mut buffers: Vec<Vec<u8>> = Vec::new();
    if let Some(arr) = root.get("buffers").and_then(|v| v.as_arr()) {
        for (bi, buf) in arr.iter().enumerate() {
            let uri = buf.get("uri").and_then(|v| v.as_str()).unwrap_or("").to_string();
            let want = buf.get("byteLength").and_then(|v| v.as_num()).map(|n| n as usize).unwrap_or(0);
            let mut bytes = if let Some(data) = uri.strip_prefix("data:") {
                let b64 = data.split("base64,").nth(1)
                    .ok_or_else(|| format!("gltf: buffer {bi}: malformed data uri"))?;
                base64_decode(b64)?
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
                comp: ac.get("componentType").and_then(|v| v.as_num()).map(|n| n as u32).unwrap_or(0),
                count: ac.get("count").and_then(|v| v.as_num()).map(|n| n as usize).unwrap_or(0),
            });
        }
    }

    // images → uri (file path or data: URI), textures → source image index,
    // materials → PBR base-color factor + base-color texture reference.
    let mut images: Vec<String> = Vec::new();
    if let Some(arr) = root.get("images").and_then(|v| v.as_arr()) {
        for im in arr {
            images.push(im.get("uri").and_then(|v| v.as_str()).unwrap_or("").to_string());
        }
    }
    let mut textures: Vec<Option<usize>> = Vec::new();
    if let Some(arr) = root.get("textures").and_then(|v| v.as_arr()) {
        for tx in arr {
            textures.push(tx.get("source").and_then(|v| v.as_num()).map(|n| n as usize));
        }
    }
    let mut gl_materials: Vec<GlMat> = Vec::new();
    if let Some(arr) = root.get("materials").and_then(|v| v.as_arr()) {
        for ma in arr {
            let mut m = GlMat {
                name: ma.get("name").and_then(|v| v.as_str()).unwrap_or("").to_string(),
                diffuse: (0.7, 0.7, 0.7),
                opacity: 1.0,
                base_color_texture: None,
            };
            if let Some(pbr) = ma.get("pbrMetallicRoughness") {
                if let Some(bcf) = pbr.get("baseColorFactor").and_then(|v| v.as_arr()) {
                    let c = arr4(bcf);
                    m.diffuse = (c[0], c[1], c[2]);
                    m.opacity = c[3];
                }
                if let Some(bct) = pbr.get("baseColorTexture") {
                    m.base_color_texture = bct.get("index").and_then(|v| v.as_num()).map(|n| n as usize);
                }
            }
            gl_materials.push(m);
        }
    }

    // meshes → primitives
    let mut meshes: Vec<GlMesh> = Vec::new();
    if let Some(arr) = root.get("meshes").and_then(|v| v.as_arr()) {
        for ms in arr {
            let mut prims = Vec::new();
            if let Some(p_arr) = ms.get("primitives").and_then(|v| v.as_arr()) {
                for pr in p_arr {
                    let mut p = Prim::default();
                    if let Some(attrs) = pr.get("attributes").and_then(|v| v.as_obj()) {
                        for (k, val) in attrs {
                            if k == "POSITION" {
                                let ai = val.as_num().map(|n| n as usize).unwrap_or(0);
                                let a = &accessors[ai];
                                p.pos_view = a.view;
                                p.pos_offset = a.offset;
                                p.pos_count = a.count;
                                p.pos_comp = a.comp;
                            }
                        }
                    }
                    if let Some(ii) = pr.get("indices").and_then(|v| v.as_num()) {
                        let a = &accessors[ii as usize];
                        p.idx_view = a.view;
                        p.idx_offset = a.offset;
                        p.idx_count = a.count;
                        p.idx_comp = a.comp;
                    }
                    p.mat_idx = pr.get("material").and_then(|v| v.as_num()).map(|n| n as usize);
                    prims.push(p);
                }
            }
            meshes.push(GlMesh { prims });
        }
    }

    // nodes
    let mut nodes: Vec<GlNode> = Vec::new();
    if let Some(arr) = root.get("nodes").and_then(|v| v.as_arr()) {
        for nd in arr {
            nodes.push(GlNode {
                name: nd.get("name").and_then(|v| v.as_str()).unwrap_or("").to_string(),
                mesh: nd.get("mesh").and_then(|v| v.as_num()).map(|n| n as usize),
                trsf: node_transform(nd),
                children: nd.get("children").and_then(|v| v.as_arr())
                    .map(|a| a.iter().filter_map(|c| c.as_num().map(|n| n as usize)).collect())
                    .unwrap_or_default(),
            });
        }
    }
    if nodes.is_empty() {
        return Err("gltf: no nodes".into());
    }

    // scene roots
    let roots: Vec<usize> = root.get("scenes").and_then(|v| v.as_arr())
        .and_then(|a| a.first())
        .and_then(|s| s.get("nodes"))
        .and_then(|v| v.as_arr())
        .map(|a| a.iter().filter_map(|c| c.as_num().map(|n| n as usize)).collect())
        .unwrap_or_else(|| vec![0]);

    let mut scene_nodes = Vec::new();
    let gltf_ctx = GltfCtx {
        meshes: &meshes,
        buffers: &buffers,
        views: &views,
        materials: &gl_materials,
        textures: &textures,
        images: &images,
        base_dir: Path::new(path).parent().unwrap_or_else(|| Path::new(".")),
    };
    for &r in &roots {
        collect_gltf_nodes(&nodes, r, &GpTrsf::identity(), &gltf_ctx, &mut scene_nodes)?;
    }

    let materials = gl_materials
        .iter()
        .map(|m| Material {
            name: m.name.clone(),
            diffuse: m.diffuse,
            specular: (1.0, 1.0, 1.0),
            opacity: m.opacity,
        })
        .collect();
    Ok(MeshScene { name: file_stem(path).to_string(), nodes: scene_nodes, materials })
}

// ---------------------------------------------------------------------------
// Scene → shape
// ---------------------------------------------------------------------------

/// Convert one node's mesh to a B-Rep shape. A closed 2-manifold mesh becomes a
/// solid; otherwise the shell is returned. `tol` is unused (vertex welding is
/// fixed at 1e-9 by the underlying conversion).
pub fn scene_to_shape(scene: &MeshScene, node_index: usize, _tol: f64) -> Result<TopoShape, String> {
    let node = scene.nodes.get(node_index)
        .ok_or_else(|| format!("rwmesh: node index {node_index} out of range ({} nodes)", scene.nodes.len()))?;
    let nodes = node.vertices.clone();
    let triangles = node.triangles.iter().map(|&(a, b, c)| Triangle::new(a, b, c)).collect();
    let b = triangulation_to_brep(&Triangulation::new(nodes, triangles));
    Ok(b.solid.map(|s| s.0).unwrap_or(b.shell.0))
}

/// Convert every node of the scene to a shape and assemble them into a compound.
pub fn scene_to_compound(scene: &MeshScene, tol: f64) -> Result<TopoShape, String> {
    let b = TopoBuilder::new();
    let mut comp = Compound::new();
    for i in 0..scene.nodes.len() {
        let s = scene_to_shape(scene, i, tol)?;
        b.add_compound(&mut comp, &s);
    }
    Ok(comp.0)
}

/// Reverse direction: mesh a `TopoShape` into a one-node scene
/// (wraps `crate::shape_mesh::mesh_shape`).
pub fn mesh_from_shape(shape: &TopoShape, deflection: f64) -> MeshScene {
    let m = crate::shape_mesh::mesh_shape(shape, deflection);
    MeshScene {
        name: "mesh".to_string(),
        nodes: vec![MeshNode {
            name: "mesh".to_string(),
            vertices: m.vertices,
            triangles: m.triangles.iter().map(|t| (t.n0, t.n1, t.n2)).collect(),
            transform: None,
            material: None,
            uv: None,
            texture_path: None,
        }],
        materials: Vec::new(),
    }
}

// ---------------------------------------------------------------------------
// Format router
// ---------------------------------------------------------------------------

/// Read any supported mesh format (by source extension) and write it to
/// `dst_path` (by destination extension: `.obj`/`.ply`/`.stl`/`.gltf`).
pub fn convert_mesh_format(src_path: &str, dst_path: &str) -> Result<(), String> {
    let scene = match ext(src_path).as_str() {
        "obj" => read_obj_scene(src_path)?,
        "ply" => read_ply_scene(src_path)?,
        "stl" => read_stl_scene(src_path)?,
        "gltf" => read_gltf_scene(src_path)?,
        "wrl" => read_vrml_scene(src_path)?,
        e => return Err(format!("rwmesh: unsupported source format .{e}")),
    };
    write_scene(&scene, dst_path)
}

fn write_scene(scene: &MeshScene, path: &str) -> Result<(), String> {
    match ext(path).as_str() {
        "obj" => {
            let mut out = String::new();
            for node in &scene.nodes {
                out.push_str(&format!("o {}\n", node.name));
                for v in &node.vertices {
                    out.push_str(&format!("v {} {} {}\n", v.x(), v.y(), v.z()));
                }
                for &(a, b, c) in &node.triangles {
                    out.push_str(&format!("f {} {} {}\n", a + 1, b + 1, c + 1));
                }
            }
            std::fs::write(path, out).map_err(|e| format!("rwmesh: write {path}: {e}"))
        }
        "ply" => {
            let mut ply = PlyMesh::default();
            for node in &scene.nodes {
                let base = ply.vertices.len();
                ply.vertices.extend(node.vertices.iter().copied());
                for &(a, b, c) in &node.triangles {
                    ply.faces.push(vec![base + a, base + b, base + c]);
                }
            }
            std::fs::write(path, occt_core::io::ply::write_ply(&ply))
                .map_err(|e| format!("rwmesh: write {path}: {e}"))
        }
        "stl" => {
            let mut stl = StlMesh::default();
            for node in &scene.nodes {
                for &(a, b, c) in &node.triangles {
                    stl.triangles.push([node.vertices[a], node.vertices[b], node.vertices[c]]);
                }
            }
            std::fs::write(path, occt_core::io::stl::write_ascii_stl(&stl))
                .map_err(|e| format!("rwmesh: write {path}: {e}"))
        }
        "gltf" => {
            let shape = scene_to_compound(scene, 1e-7)?;
            crate::gltf::write_gltf_file(&shape, 0.1, path, &GltfOptions::default())
        }
        e => Err(format!("rwmesh: unsupported destination format .{e}")),
    }
}

// ---------------------------------------------------------------------------
// Scene helpers
// ---------------------------------------------------------------------------

/// Total number of triangles across all nodes.
pub fn mesh_triangle_count(scene: &MeshScene) -> usize {
    scene.nodes.iter().map(|n| n.triangles.len()).sum()
}

/// Total number of vertices across all nodes.
pub fn mesh_vertex_count(scene: &MeshScene) -> usize {
    scene.nodes.iter().map(|n| n.vertices.len()).sum()
}

/// Axis-aligned bounds of all scene vertices, or `None` for an empty scene.
pub fn scene_bounds(scene: &MeshScene) -> Option<(GpPnt, GpPnt)> {
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    let mut any = false;
    for n in &scene.nodes {
        for v in &n.vertices {
            any = true;
            min[0] = min[0].min(v.x());
            min[1] = min[1].min(v.y());
            min[2] = min[2].min(v.z());
            max[0] = max[0].max(v.x());
            max[1] = max[1].max(v.y());
            max[2] = max[2].max(v.z());
        }
    }
    if any {
        Some((GpPnt::new(min[0], min[1], min[2]), GpPnt::new(max[0], max[1], max[2])))
    } else {
        None
    }
}

/// Look up a named material in the scene's material library.
pub fn material_from_name<'a>(scene: &'a MeshScene, name: &str) -> Option<&'a Material> {
    scene.materials.iter().find(|m| m.name == name)
}

/// Render a node's texture space to a crude PPM checkerboard. With UVs present,
/// each pixel's quadrant color is chosen from its (u,v) coordinate; without UVs
/// the node's material diffuse (or gray) fills the image. Returns P3 (ASCII) PPM.
pub fn scene_with_texture_to_ppm(scene: &MeshScene, node_index: usize, width: u32, height: u32) -> Result<String, String> {
    let node = scene.nodes.get(node_index)
        .ok_or_else(|| format!("rwmesh: node index {node_index} out of range ({} nodes)", scene.nodes.len()))?;
    let mut out = String::new();
    out.push_str(&format!("P3\n{} {}\n255\n", width, height));
    let (w, h) = (width.max(1), height.max(1));
    for y in 0..height {
        for x in 0..width {
            let (r, g, b) = if node.uv.is_some() {
                let u = x as f64 / w as f64;
                let v = y as f64 / h as f64;
                match (u < 0.5, v < 0.5) {
                    (true, true) => (1.0, 0.0, 0.0),
                    (true, false) => (0.0, 1.0, 0.0),
                    (false, true) => (0.0, 0.0, 1.0),
                    (false, false) => (1.0, 1.0, 0.0),
                }
            } else {
                node.material.as_ref().map(|m| m.diffuse).unwrap_or((0.7, 0.7, 0.7))
            };
            let px = |c: f64| (c.clamp(0.0, 1.0) * 255.0).round() as u32;
            out.push_str(&format!("{} {} {}\n", px(r), px(g), px(b)));
        }
    }
    Ok(out)
}

// ---------------------------------------------------------------------------
// OBJ scene parser (group-aware)
// ---------------------------------------------------------------------------

/// Parse OBJ content without a base directory (any `mtllib` reference is
/// resolved against the current directory).
fn parse_obj_scene(content: &str, name: &str) -> Result<MeshScene, String> {
    parse_obj_scene_inner(content, name, None)
}

/// One OBJ group: a named collection of polygon faces; each corner is a
/// `(vertex_index, texcoord_index)` pair. `material` is the `usemtl` name in
/// effect when the group started.
struct ObjGroup {
    name: String,
    faces: Vec<Vec<(usize, Option<usize>)>>,
    material: Option<String>,
}

fn parse_obj_scene_inner(
    content: &str,
    name: &str,
    dir: Option<&Path>,
) -> Result<MeshScene, String> {
    let mut vertices: Vec<GpPnt> = Vec::new();
    let mut texcoords: Vec<GpPnt2d> = Vec::new();
    let mut groups: Vec<ObjGroup> = Vec::new();
    let mut current = 0usize;
    let mut mtllibs: Vec<String> = Vec::new();
    let mut cur_material: Option<String> = None;

    for (lineno, raw) in content.lines().enumerate() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut parts = line.split_whitespace();
        let tag = parts.next().unwrap_or("");
        match tag {
            "v" => {
                let c: Vec<f64> = parts.filter_map(|p| p.parse().ok()).collect();
                if c.len() < 3 {
                    return Err(format!("rwmesh: obj line {}: bad vertex", lineno + 1));
                }
                vertices.push(GpPnt::new(c[0], c[1], c[2]));
            }
            "vt" => {
                let c: Vec<f64> = parts.filter_map(|p| p.parse().ok()).collect();
                if c.len() >= 2 {
                    texcoords.push(GpPnt2d::new(c[0], c[1]));
                }
            }
            "mtllib" => {
                for f in parts {
                    mtllibs.push(f.to_string());
                }
            }
            "usemtl" => {
                if let Some(m) = parts.next() {
                    cur_material = Some(m.to_string());
                }
            }
            "o" | "g" => {
                let gname = parts.next().unwrap_or("mesh").to_string();
                groups.push(ObjGroup {
                    name: gname,
                    faces: Vec::new(),
                    material: cur_material.clone(),
                });
                current = groups.len() - 1;
            }
            "f" => {
                if groups.is_empty() {
                    groups.push(ObjGroup {
                        name: "mesh".to_string(),
                        faces: Vec::new(),
                        material: cur_material.clone(),
                    });
                    current = 0;
                }
                let mut corners = Vec::new();
                for tok in parts {
                    let mut seg = tok.split('/');
                    let vi: i32 = seg.next().unwrap_or("")
                        .parse()
                        .map_err(|_| format!("rwmesh: obj line {}: bad face index", lineno + 1))?;
                    let v_idx = obj_index(vi, vertices.len())?;
                    let vt_idx = seg.next()
                        .filter(|s| !s.is_empty())
                        .and_then(|s| s.parse::<i32>().ok())
                        .map(|t| obj_index(t, texcoords.len()))
                        .transpose()?;
                    corners.push((v_idx, vt_idx));
                }
                groups[current].faces.push(corners);
            }
            _ => {}
        }
    }

    // Load the material libraries referenced by `mtllib`.
    let mut materials: Vec<Material> = Vec::new();
    for ml in &mtllibs {
        let mtl_path = match dir {
            Some(d) => d.join(ml),
            None => Path::new(ml).to_path_buf(),
        };
        let mtl_content = std::fs::read_to_string(&mtl_path)
            .map_err(|e| format!("rwmesh: read mtl {}: {e}", mtl_path.display()))?;
        materials.extend(parse_mtl(&mtl_content));
    }

    let nodes = if groups.is_empty() {
        vec![MeshNode {
            name: "mesh".to_string(),
            vertices,
            triangles: Vec::new(),
            transform: None,
            material: None,
            uv: None,
            texture_path: None,
        }]
    } else {
        let mut out_nodes = Vec::with_capacity(groups.len());
        for g in groups {
            // Remap global OBJ (vertex, texcoord) pairs to a local list so a
            // shared vertex with different UVs (a texture seam) becomes two.
            let mut map: HashMap<(usize, Option<usize>), usize> = HashMap::new();
            let mut local: Vec<GpPnt> = Vec::new();
            let mut uv: Vec<GpPnt2d> = Vec::new();
            let mut tris = Vec::new();
            let has_uv = g.faces.iter().any(|f| f.iter().any(|&(_, vt)| vt.is_some()));
            for face in &g.faces {
                for i in 1..face.len().saturating_sub(1) {
                    let corners = [face[0], face[i], face[i + 1]];
                    let mut tri = [0usize; 3];
                    for (k, &(vi, vti)) in corners.iter().enumerate() {
                        let id = *map.entry((vi, vti)).or_insert_with(|| {
                            local.push(vertices[vi]);
                            uv.push(match vti {
                                Some(t) => texcoords[t],
                                None => GpPnt2d::new(0.0, 0.0),
                            });
                            local.len() - 1
                        });
                        tri[k] = id;
                    }
                    tris.push((tri[0], tri[1], tri[2]));
                }
            }
            let material = g.material.as_ref().and_then(|mn| find_material(&materials, mn)).cloned();
            out_nodes.push(MeshNode {
                name: g.name,
                vertices: local,
                triangles: tris,
                transform: None,
                material,
                uv: if has_uv { Some(uv) } else { None },
                texture_path: None,
            });
        }
        out_nodes
    };

    Ok(MeshScene { name: name.to_string(), nodes, materials })
}

/// Parse a Wavefront `.mtl` material library: `newmtl <name>` starts a material,
/// `Kd`/`Ks` set the diffuse/specular RGB, `d` sets opacity.
fn parse_mtl(content: &str) -> Vec<Material> {
    let mut mats: Vec<Material> = Vec::new();
    let mut cur: Option<Material> = None;
    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut parts = line.split_whitespace();
        match parts.next().unwrap_or("") {
            "newmtl" => {
                if let Some(m) = cur.take() {
                    mats.push(m);
                }
                cur = Some(Material {
                    name: parts.next().unwrap_or("").to_string(),
                    ..Material::default()
                });
            }
            "Kd" => {
                if let Some(m) = &mut cur {
                    let c: Vec<f64> = parts.filter_map(|p| p.parse().ok()).collect();
                    if c.len() >= 3 {
                        m.diffuse = (c[0], c[1], c[2]);
                    }
                }
            }
            "Ks" => {
                if let Some(m) = &mut cur {
                    let c: Vec<f64> = parts.filter_map(|p| p.parse().ok()).collect();
                    if c.len() >= 3 {
                        m.specular = (c[0], c[1], c[2]);
                    }
                }
            }
            "d" => {
                if let Some(m) = &mut cur {
                    if let Some(v) = parts.next().and_then(|p| p.parse().ok()) {
                        m.opacity = v;
                    }
                }
            }
            _ => {}
        }
    }
    if let Some(m) = cur {
        mats.push(m);
    }
    mats
}

fn find_material<'a>(materials: &'a [Material], name: &str) -> Option<&'a Material> {
    materials.iter().find(|m| m.name == name)
}

/// Resolve a 1-based (or negative = from-end) OBJ index to a 0-based index.
fn obj_index(i: i32, len: usize) -> Result<usize, String> {
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
enum VTok {
    Id(String),
    Num(f64),
    LBrace,
    RBrace,
    LBracket,
    RBracket,
}

/// Tokenize VRML 2.0: identifiers/numbers, `{}[]`, quoted strings become `Id`,
/// `#` starts a comment, commas and whitespace separate tokens.
fn tokenize_vrml(content: &str) -> Result<Vec<VTok>, String> {
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
struct Vp {
    toks: Vec<VTok>,
    i: usize,
}

impl Vp {
    fn peek(&self) -> Option<&VTok> {
        self.toks.get(self.i)
    }
    fn peek_is(&self, t: &VTok) -> bool {
        self.peek() == Some(t)
    }
    fn peek_lbrace(&self) -> bool {
        self.peek_is(&VTok::LBrace)
    }
    fn peek_rbrace(&self) -> bool {
        self.peek_is(&VTok::RBrace)
    }
    fn peek_rbracket(&self) -> bool {
        self.peek_is(&VTok::RBracket)
    }
    fn peek_id_eq(&self, s: &str) -> bool {
        matches!(self.peek(), Some(VTok::Id(x)) if x == s)
    }
    fn take(&mut self) -> Option<VTok> {
        let t = self.toks.get(self.i).cloned();
        if t.is_some() {
            self.i += 1;
        }
        t
    }
    fn expect_id(&mut self) -> Result<String, String> {
        match self.take() {
            Some(VTok::Id(s)) => Ok(s),
            other => Err(format!("vrml: expected identifier, got {other:?}")),
        }
    }
    fn expect_num(&mut self) -> Result<f64, String> {
        match self.take() {
            Some(VTok::Num(n)) => Ok(n),
            other => Err(format!("vrml: expected number, got {other:?}")),
        }
    }
    fn expect_lbrace(&mut self) -> Result<(), String> {
        if self.take() == Some(VTok::LBrace) {
            Ok(())
        } else {
            Err("vrml: expected '{{'".into())
        }
    }
    fn expect_rbrace(&mut self) -> Result<(), String> {
        if self.take() == Some(VTok::RBrace) {
            Ok(())
        } else {
            Err("vrml: expected '}}'".into())
        }
    }
    fn expect_lbracket(&mut self) -> Result<(), String> {
        if self.take() == Some(VTok::LBracket) {
            Ok(())
        } else {
            Err("vrml: expected '['".into())
        }
    }
    fn expect_rbracket(&mut self) -> Result<(), String> {
        if self.take() == Some(VTok::RBracket) {
            Ok(())
        } else {
            Err("vrml: expected ']'".into())
        }
    }
    fn parse_vec2(&mut self) -> Result<[f64; 2], String> {
        Ok([self.expect_num()?, self.expect_num()?])
    }
    fn parse_vec3(&mut self) -> Result<[f64; 3], String> {
        Ok([self.expect_num()?, self.expect_num()?, self.expect_num()?])
    }
    fn parse_vec4(&mut self) -> Result<[f64; 4], String> {
        Ok([self.expect_num()?, self.expect_num()?, self.expect_num()?, self.expect_num()?])
    }

    /// Skip one field value: `[ ... ]`, a `Node { ... }`, or a scalar run.
    fn skip_value(&mut self) -> Result<(), String> {
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

    fn skip_to_matching(&mut self, close: &VTok, open: &VTok) -> Result<(), String> {
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
    fn skip_node(&mut self) -> Result<(), String> {
        self.expect_id()?;
        self.expect_lbrace()?;
        self.skip_to_matching(&VTok::RBrace, &VTok::LBrace)
    }

    fn parse_scene_node(&mut self, world: &GpTrsf, scene: &mut MeshScene) -> Result<(), String> {
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

    fn parse_transform(&mut self, world: &GpTrsf, scene: &mut MeshScene) -> Result<(), String> {
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

    fn parse_children(&mut self, world: &GpTrsf, scene: &mut MeshScene) -> Result<(), String> {
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

    fn parse_shape(&mut self, world: &GpTrsf, scene: &mut MeshScene) -> Result<(), String> {
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

    fn parse_appearance(&mut self) -> Result<Option<Material>, String> {
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

    fn parse_material(&mut self) -> Result<Material, String> {
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

    fn parse_indexed_face_set(&mut self) -> Result<(Vec<GpPnt>, Vec<i32>, Vec<GpPnt2d>), String> {
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

    fn parse_coordinate(&mut self) -> Result<Vec<GpPnt>, String> {
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

    fn parse_texture_coordinate(&mut self) -> Result<Vec<GpPnt2d>, String> {
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

    fn parse_i32_array(&mut self) -> Result<Vec<i32>, String> {
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
fn trsf_from_vrml(t: &[f64; 3], axis: &[f64; 3], angle: f64, s: &[f64; 3]) -> GpTrsf {
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
fn build_vrml_node(
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
struct BufView {
    buffer: usize,
    offset: usize,
    len: usize,
}

#[derive(Default, Clone, Copy)]
struct Accessor {
    view: Option<usize>,
    offset: usize,
    comp: u32,
    count: usize,
}

#[derive(Default, Clone, Copy)]
struct Prim {
    pos_view: Option<usize>,
    pos_offset: usize,
    pos_count: usize,
    pos_comp: u32,
    idx_view: Option<usize>,
    idx_offset: usize,
    idx_count: usize,
    idx_comp: u32,
    mat_idx: Option<usize>,
}

struct GlMesh {
    prims: Vec<Prim>,
}

struct GlNode {
    name: String,
    mesh: Option<usize>,
    trsf: GpTrsf,
    children: Vec<usize>,
}

/// glTF material: base-color diffuse + opacity plus an optional base-color
/// texture reference (texture index → `textures[].source` → image).
struct GlMat {
    name: String,
    diffuse: (f64, f64, f64),
    opacity: f64,
    base_color_texture: Option<usize>,
}

/// Shared resources handed to the node-graph walk.
struct GltfCtx<'a> {
    meshes: &'a [GlMesh],
    buffers: &'a [Vec<u8>],
    views: &'a [BufView],
    materials: &'a [GlMat],
    textures: &'a [Option<usize>],
    images: &'a [String],
    base_dir: &'a Path,
}

/// Build the local transform for a glTF node from `matrix` or `translation`/
/// `rotation`/`scale` (matrix wins per the spec).
fn node_transform(nd: &JVal) -> GpTrsf {
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
fn trsf_from_trs(t: &[f64; 3], r: &[f64; 4], s: &[f64; 3]) -> GpTrsf {
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

/// Column-major glTF 4x4 matrix → `GpTrsf` (linear part + translation).
fn trsf_from_matrix(m: &[f64]) -> GpTrsf {
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
fn collect_gltf_nodes(
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
fn resolve_uri(uri: &str, base_dir: &Path) -> String {
    if uri.starts_with("data:") {
        uri.to_string()
    } else {
        base_dir.join(uri).to_string_lossy().to_string()
    }
}

/// Read a mesh primitive's positions and indices from the buffer views.
fn read_primitive(prim: &Prim, buffers: &[Vec<u8>], views: &[BufView]) -> Result<(Vec<GpPnt>, Vec<(usize, usize, usize)>), String> {
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

fn view_bytes<'a>(view: Option<usize>, buffers: &'a [Vec<u8>], views: &[BufView]) -> Result<&'a [u8], String> {
    let v = views.get(view.ok_or("gltf: missing bufferView")?)
        .ok_or("gltf: bufferView index out of range")?;
    let buf = buffers.get(v.buffer).ok_or("gltf: buffer index out of range")?;
    if v.offset + v.len > buf.len() {
        return Err("gltf: bufferView extends past buffer end".into());
    }
    Ok(&buf[v.offset..v.offset + v.len])
}

fn read_uint(comp: u32, b: &[u8], off: usize) -> usize {
    match comp {
        5121 => b[off] as usize,
        5123 => u16::from_le_bytes([b[off], b[off + 1]]) as usize,
        _ => u32::from_le_bytes([b[off], b[off + 1], b[off + 2], b[off + 3]]) as usize,
    }
}

/// Decode a standard base64 string (RFC 4648) with or without padding.
fn base64_decode(s: &str) -> Result<Vec<u8>, String> {
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
enum JVal {
    Null,
    Bool(bool),
    Num(f64),
    Str(String),
    Arr(Vec<JVal>),
    Obj(Vec<(String, JVal)>),
}

impl JVal {
    fn as_arr(&self) -> Option<&[JVal]> {
        match self {
            JVal::Arr(a) => Some(a),
            _ => None,
        }
    }
    fn as_obj(&self) -> Option<&[(String, JVal)]> {
        match self {
            JVal::Obj(o) => Some(o),
            _ => None,
        }
    }
    fn as_str(&self) -> Option<&str> {
        match self {
            JVal::Str(s) => Some(s),
            _ => None,
        }
    }
    fn as_num(&self) -> Option<f64> {
        match self {
            JVal::Num(n) => Some(*n),
            _ => None,
        }
    }
    fn get(&self, key: &str) -> Option<&JVal> {
        self.as_obj().and_then(|o| o.iter().find(|(k, _)| k == key).map(|(_, v)| v))
    }
}

struct JP<'a> {
    s: &'a [u8],
    i: usize,
}

impl<'a> JP<'a> {
    fn ws(&mut self) {
        while self.i < self.s.len() && matches!(self.s[self.i], b' ' | b'\t' | b'\n' | b'\r') {
            self.i += 1;
        }
    }
    fn peek(&self) -> Option<u8> {
        self.s.get(self.i).copied()
    }
    fn eat(&mut self, c: u8) -> Result<(), String> {
        self.ws();
        if self.peek() == Some(c) {
            self.i += 1;
            Ok(())
        } else {
            Err(format!("json: expected '{}' at {}", c as char, self.i))
        }
    }
    fn value(&mut self) -> Result<JVal, String> {
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
    fn lit(&mut self, word: &str) -> Result<(), String> {
        if self.s.get(self.i..).map(|x| x.starts_with(word.as_bytes())).unwrap_or(false) {
            self.i += word.len();
            Ok(())
        } else {
            Err(format!("json: expected '{word}'"))
        }
    }
    fn string(&mut self) -> Result<String, String> {
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
    fn number(&mut self) -> Result<f64, String> {
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
    fn array(&mut self) -> Result<JVal, String> {
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
    fn object(&mut self) -> Result<JVal, String> {
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

fn parse_json(text: &str) -> Result<JVal, String> {
    let mut p = JP { s: text.as_bytes(), i: 0 };
    let v = p.value()?;
    p.ws();
    if p.i != p.s.len() {
        return Err("json: trailing content".into());
    }
    Ok(v)
}

fn arr3(a: &[JVal]) -> [f64; 3] {
    let mut out = [0.0; 3];
    for (k, v) in a.iter().take(3).enumerate() {
        if let Some(n) = v.as_num() {
            out[k] = n;
        }
    }
    out
}

fn arr4(a: &[JVal]) -> [f64; 4] {
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

fn file_stem(path: &str) -> &str {
    Path::new(path).file_stem().and_then(|s| s.to_str()).unwrap_or("mesh")
}

fn ext(path: &str) -> String {
    Path::new(path)
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_lowercase()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh::mesh_box;
    use crate::primitives::{BRepPrimBox, BRepPrimSphere};

    fn temp_subdir(name: &str) -> std::path::PathBuf {
        let mut d = std::env::temp_dir();
        d.push(format!("rwmesh_{}_{}", std::process::id(), name));
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    fn triangle_bin() -> Vec<u8> {
        let mut bin = Vec::new();
        for p in [GpPnt::new(0., 0., 0.), GpPnt::new(1., 0., 0.), GpPnt::new(0., 1., 0.)] {
            for c in [p.x() as f32, p.y() as f32, p.z() as f32] {
                bin.extend_from_slice(&c.to_le_bytes());
            }
        }
        for i in [0u32, 1, 2] {
            bin.extend_from_slice(&i.to_le_bytes());
        }
        bin
    }

    /// Single-triangle glTF JSON; `node_extra` is spliced into the node object
    /// (e.g. a `,"matrix":[...]` attribute).
    fn tri_gltf_json(uri: &str, node_extra: &str) -> String {
        format!(
            r#"{{"asset":{{"version":"2.0"}},"scene":0,"scenes":[{{"nodes":[0]}}],"nodes":[{{"mesh":0{node_extra}}}],"meshes":[{{"primitives":[{{"attributes":{{"POSITION":0}},"indices":1}}]}}],"buffers":[{{"uri":"{uri}","byteLength":48}}],"bufferViews":[{{"buffer":0,"byteOffset":0,"byteLength":36,"target":34962}},{{"buffer":0,"byteOffset":36,"byteLength":12,"target":34963}}],"accessors":[{{"bufferView":0,"componentType":5126,"count":3,"type":"VEC3"}},{{"bufferView":1,"componentType":5125,"count":3,"type":"SCALAR"}}]}}"#
        )
    }

    #[test]
    fn read_obj_roundtrip() {
        use occt_core::io::obj::{ObjFace, ObjMesh};
        let mesh = ObjMesh {
            vertices: vec![
                GpPnt::new(0., 0., 0.),
                GpPnt::new(1., 0., 0.),
                GpPnt::new(0., 1., 0.),
                GpPnt::new(1., 1., 0.),
            ],
            texcoords: vec![],
            normals: vec![],
            faces: vec![
                ObjFace { v: vec![0, 1, 2], vt: None, vn: None },
                ObjFace { v: vec![1, 3, 2], vt: None, vn: None },
            ],
        };
        let path = temp_subdir("obj").join("m.obj");
        std::fs::write(&path, occt_core::io::obj::write_obj(&mesh)).unwrap();
        let scene = read_obj_scene(&path.to_string_lossy()).unwrap();
        assert_eq!(mesh_vertex_count(&scene), 4);
        assert_eq!(mesh_triangle_count(&scene), 2);
    }

    #[test]
    fn read_ply_roundtrip() {
        let ply = PlyMesh {
            vertices: vec![
                GpPnt::new(0., 0., 0.),
                GpPnt::new(1., 0., 0.),
                GpPnt::new(1., 1., 0.),
                GpPnt::new(0., 1., 0.),
            ],
            faces: vec![vec![0, 1, 2], vec![0, 2, 3]],
        };
        let path = temp_subdir("ply").join("m.ply");
        std::fs::write(&path, occt_core::io::ply::write_ply(&ply)).unwrap();
        let scene = read_ply_scene(&path.to_string_lossy()).unwrap();
        assert_eq!(mesh_vertex_count(&scene), 4);
        assert_eq!(mesh_triangle_count(&scene), 2);
    }

    #[test]
    fn read_stl_ascii() {
        let text = "solid demo\n\
            facet normal 0 0 1\n  outer loop\n    vertex 0 0 0\n    vertex 1 0 0\n    vertex 0 1 0\n  endloop\nendfacet\n\
            endsolid demo\n";
        let path = temp_subdir("stl_a").join("m.stl");
        std::fs::write(&path, text).unwrap();
        let scene = read_stl_scene(&path.to_string_lossy()).unwrap();
        assert_eq!(mesh_triangle_count(&scene), 1);
        assert_eq!(mesh_vertex_count(&scene), 3);
    }

    #[test]
    fn read_stl_binary() {
        let stl = StlMesh {
            triangles: vec![
                [GpPnt::new(0., 0., 0.), GpPnt::new(1., 0., 0.), GpPnt::new(0., 1., 0.)],
                [GpPnt::new(1., 0., 0.), GpPnt::new(1., 1., 0.), GpPnt::new(0., 1., 0.)],
            ],
            normals: vec![],
        };
        let path = temp_subdir("stl_b").join("m.stl");
        std::fs::write(&path, occt_core::io::stl::write_binary_stl(&stl)).unwrap();
        let scene = read_stl_scene(&path.to_string_lossy()).unwrap();
        assert_eq!(mesh_triangle_count(&scene), 2);
        // Shared vertices deduplicated.
        assert_eq!(mesh_vertex_count(&scene), 4);
    }

    #[test]
    fn read_gltf_embedded_base64() {
        let bin = triangle_bin();
        let uri = format!("data:application/octet-stream;base64,{}", crate::gltf::base64_encode(&bin));
        let path = temp_subdir("gltf_b64").join("m.gltf");
        std::fs::write(&path, tri_gltf_json(&uri, "")).unwrap();
        let scene = read_gltf_scene(&path.to_string_lossy()).unwrap();
        assert_eq!(mesh_vertex_count(&scene), 3);
        assert_eq!(mesh_triangle_count(&scene), 1);
    }

    #[test]
    fn read_gltf_external_bin() {
        let dir = temp_subdir("gltf_ext");
        std::fs::write(dir.join("mesh.bin"), triangle_bin()).unwrap();
        let path = dir.join("m.gltf");
        std::fs::write(&path, tri_gltf_json("mesh.bin", "")).unwrap();
        let scene = read_gltf_scene(&path.to_string_lossy()).unwrap();
        assert_eq!(mesh_vertex_count(&scene), 3);
        assert_eq!(mesh_triangle_count(&scene), 1);
    }

    #[test]
    fn scene_to_compound_box() {
        let shape = BRepPrimBox::make_box(2.0, 2.0, 2.0).solid.0;
        let scene = mesh_from_shape(&shape, 0.5);
        let comp = scene_to_compound(&scene, 1e-7).expect("compound");
        assert!(crate::shape_mesh::count_faces(&comp) > 0, "no faces in compound");
        assert!(mesh_vertex_count(&scene) > 0);

        // Closedness: a clean watertight box mesh converts to a solid.
        let m = mesh_box((GpPnt::zero(), GpPnt::new(2.0, 2.0, 2.0)));
        let scene2 = MeshScene {
            name: "box".into(),
            nodes: vec![MeshNode {
                name: "box".into(),
                vertices: m.vertices,
                triangles: m.triangles.iter().map(|t| (t.n0, t.n1, t.n2)).collect(),
                transform: None,
                material: None,
                uv: None,
                texture_path: None,
            }],
            materials: Vec::new(),
        };
        let s = scene_to_shape(&scene2, 0, 1e-7).expect("box shape");
        assert!(s.is_solid(), "watertight box should produce a solid");
    }

    #[test]
    fn scene_bounds_ok() {
        let scene = MeshScene {
            name: "s".into(),
            nodes: vec![MeshNode {
                name: "n".into(),
                vertices: vec![
                    GpPnt::new(1., 2., 3.),
                    GpPnt::new(-1., -2., -3.),
                    GpPnt::new(4., 0., 5.),
                ],
                triangles: vec![(0, 1, 2)],
                transform: None,
                material: None,
                uv: None,
                texture_path: None,
            }],
            materials: Vec::new(),
        };
        let (lo, hi) = scene_bounds(&scene).expect("bounds");
        for v in &scene.nodes[0].vertices {
            assert!(lo.x() <= v.x() && v.x() <= hi.x(), "x in bounds");
            assert!(lo.y() <= v.y() && v.y() <= hi.y(), "y in bounds");
            assert!(lo.z() <= v.z() && v.z() <= hi.z(), "z in bounds");
        }
        assert!(lo.is_equal(&GpPnt::new(-1., -2., -3.)));
        assert!(hi.is_equal(&GpPnt::new(4., 2., 5.)));
        assert!(scene_bounds(&MeshScene { name: "e".into(), nodes: vec![], materials: Vec::new() }).is_none());
    }

    #[test]
    fn convert_obj_to_ply() {
        let src = temp_subdir("conv").join("m.obj");
        let dst = temp_subdir("conv").join("m.ply");
        std::fs::write(
            &src,
            "v 0 0 0\nv 1 0 0\nv 0 1 0\nv 1 1 0\nf 1 2 3\nf 2 4 3\n",
        )
        .unwrap();
        convert_mesh_format(&src.to_string_lossy(), &dst.to_string_lossy()).unwrap();
        let scene = read_ply_scene(&dst.to_string_lossy()).unwrap();
        assert_eq!(mesh_triangle_count(&scene), 2);
        assert_eq!(mesh_vertex_count(&scene), 4);
    }

    #[test]
    fn roundtrip_shape_mesh() {
        let shape = BRepPrimSphere::make_sphere(1.0).solid.0;
        let scene = mesh_from_shape(&shape, 0.2);
        let n_before = mesh_vertex_count(&scene);
        assert!(n_before > 100, "sphere mesh too coarse: {n_before}");

        let s = scene_to_shape(&scene, 0, 1e-7).expect("sphere shape");
        // Each input triangle becomes one face.
        assert_eq!(crate::shape_mesh::count_faces(&s), mesh_triangle_count(&scene));

        // Vertex count preserved by the conversion (mesh is a single closed grid).
        let tri = Triangulation::new(
            scene.nodes[0].vertices.clone(),
            scene.nodes[0].triangles.iter().map(|&(a, b, c)| Triangle::new(a, b, c)).collect(),
        );
        let b = triangulation_to_brep(&tri);
        // The sphere grid collapses its pole rows (nu points → 1 each), so the
        // deduped conversion keeps slightly fewer vertices — but it must stay
        // close to the source count.
        let loss = n_before.saturating_sub(b.vertices.len());
        assert!(
            loss <= n_before / 4,
            "vertex count drifted: before {n_before} after {}",
            b.vertices.len()
        );
    }

    #[test]
    fn node_transform_applied_gltf() {
        let bin = triangle_bin();
        let uri = format!("data:application/octet-stream;base64,{}", crate::gltf::base64_encode(&bin));
        let dir = temp_subdir("gltf_trs");
        let path = dir.join("m.gltf");
        let json = tri_gltf_json(&uri, r#","matrix":[1,0,0,0,0,1,0,0,0,0,1,0,5,0,0,1]"#);
        std::fs::write(&path, json).unwrap();
        let scene = read_gltf_scene(&path.to_string_lossy()).unwrap();
        let node = &scene.nodes[0];
        assert!((node.vertices[0].x() - 5.0).abs() < 1e-6, "v0.x = {}", node.vertices[0].x());
        assert!((node.vertices[1].x() - 6.0).abs() < 1e-6, "v1.x = {}", node.vertices[1].x());
        assert!((node.vertices[2].y() - 1.0).abs() < 1e-6, "v2.y = {}", node.vertices[2].y());
        assert!(node.transform.is_some(), "node transform recorded");
    }

    #[test]
    fn obj_groups_become_nodes() {
        let obj = "o partA\nv 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\no partB\nv 5 5 5\nv 6 5 5\nv 5 6 5\nf 4 5 6\n";
        let scene = parse_obj_scene(obj, "groups").unwrap();
        assert_eq!(scene.nodes.len(), 2, "two object groups -> two nodes");
        assert_eq!(mesh_triangle_count(&scene), 2);
        assert_eq!(scene.nodes[0].name, "partA");
        assert_eq!(scene.nodes[1].name, "partB");
        assert_eq!(mesh_vertex_count(&scene), 6);
    }

    #[test]
    fn obj_mtl_material_parsed() {
        let dir = temp_subdir("obj_mtl");
        std::fs::write(
            dir.join("colors.mtl"),
            "newmtl red\nKd 1 0 0\nKs 0.2 0.2 0.2\nd 0.5\n",
        )
        .unwrap();
        let obj = "mtllib colors.mtl\nusemtl red\no part\nv 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n";
        let path = dir.join("m.obj");
        std::fs::write(&path, obj).unwrap();
        let scene = read_obj_scene(&path.to_string_lossy()).unwrap();
        assert_eq!(scene.materials.len(), 1, "mtl loaded into library");
        assert_eq!(scene.materials[0].name, "red");
        let (r, g, b) = scene.materials[0].diffuse;
        assert!((r - 1.0).abs() < 1e-6 && g.abs() < 1e-6 && b.abs() < 1e-6, "Kd 1 0 0");
        assert!((scene.materials[0].opacity - 0.5).abs() < 1e-6, "d 0.5");
        let mat = scene.nodes[0].material.as_ref().expect("usemtl applied to node");
        let (r, g, b) = mat.diffuse;
        assert!((r - 1.0).abs() < 1e-6 && g.abs() < 1e-6 && b.abs() < 1e-6, "node diffuse red");
    }

    #[test]
    fn obj_uv_coordinates() {
        let obj = "v 0 0 0\nv 1 0 0\nv 0 1 0\nvt 0 0\nvt 1 0\nvt 0 1\nf 1/1 2/2 3/3\n";
        let scene = parse_obj_scene(obj, "uv").unwrap();
        let node = &scene.nodes[0];
        let uv = node.uv.as_ref().expect("uv present");
        assert_eq!(uv.len(), 3, "one uv per vertex");
        assert!((uv[0].x() - 0.0).abs() < 1e-9 && (uv[0].y() - 0.0).abs() < 1e-9);
        assert!((uv[1].x() - 1.0).abs() < 1e-9 && (uv[1].y() - 0.0).abs() < 1e-9);
        assert!((uv[2].x() - 0.0).abs() < 1e-9 && (uv[2].y() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn gltf_material_color() {
        let bin = triangle_bin();
        let uri = format!("data:application/octet-stream;base64,{}", crate::gltf::base64_encode(&bin));
        let json = format!(
            r#"{{"asset":{{"version":"2.0"}},"scene":0,"scenes":[{{"nodes":[0]}}],"nodes":[{{"mesh":0}}],"meshes":[{{"primitives":[{{"attributes":{{"POSITION":0}},"indices":1,"material":0}}]}}],"materials":[{{"name":"redmat","pbrMetallicRoughness":{{"baseColorFactor":[1,0,0,1]}}}}],"buffers":[{{"uri":"{uri}","byteLength":48}}],"bufferViews":[{{"buffer":0,"byteOffset":0,"byteLength":36,"target":34962}},{{"buffer":0,"byteOffset":36,"byteLength":12,"target":34963}}],"accessors":[{{"bufferView":0,"componentType":5126,"count":3,"type":"VEC3"}},{{"bufferView":1,"componentType":5125,"count":3,"type":"SCALAR"}}]}}"#
        );
        let path = temp_subdir("gltf_mat").join("m.gltf");
        std::fs::write(&path, json).unwrap();
        let scene = read_gltf_scene(&path.to_string_lossy()).unwrap();
        let mat = scene.nodes[0].material.as_ref().expect("gltf material attached");
        let (r, g, b) = mat.diffuse;
        assert!((r - 1.0).abs() < 1e-6 && g.abs() < 1e-6 && b.abs() < 1e-6, "baseColorFactor [1,0,0,1]");
        assert_eq!(scene.materials.len(), 1, "gltf materials copied to library");
        assert_eq!(scene.materials[0].name, "redmat");
    }

    fn box_vrml() -> &'static str {
        "#VRML V2.0 utf8\n\
         Shape {\n\
           geometry IndexedFaceSet {\n\
             coord Coordinate {\n\
               point [\n\
                 0 0 0, 1 0 0, 1 1 0, 0 1 0,\n\
                 0 0 1, 1 0 1, 1 1 1, 0 1 1\n\
               ]\n\
             }\n\
             coordIndex [\n\
               0,1,2,3,-1,  4,7,6,5,-1,  0,4,5,1,-1,\n\
               1,5,6,2,-1,  2,6,7,3,-1,  3,7,4,0,-1\n\
             ]\n\
           }\n\
         }\n"
    }

    #[test]
    fn vrml_simple_box() {
        let scene = parse_vrml_scene(box_vrml(), "cube").unwrap();
        assert_eq!(scene.nodes.len(), 1, "one shape -> one node");
        assert_eq!(mesh_triangle_count(&scene), 12, "six quads fan to twelve tris");
        assert_eq!(mesh_vertex_count(&scene), 8, "cube has eight points");
    }

    #[test]
    fn vrml_transform_nesting() {
        let vrml = "#VRML V2.0 utf8\n\
            Transform {\n\
              translation 2 0 0\n\
              children [\n\
                Shape {\n\
                  geometry IndexedFaceSet {\n\
                    coord Coordinate { point [ 0 0 0, 1 0 0, 0 1 0 ] }\n\
                    coordIndex [ 0,1,2,-1 ]\n\
                  }\n\
                }\n\
              ]\n\
            }\n";
        let scene = parse_vrml_scene(vrml, "moved").unwrap();
        let node = &scene.nodes[0];
        assert_eq!(node.vertices.len(), 3);
        // Vertices (0,0,0), (1,0,0), (0,1,0) shifted by translation [2,0,0].
        assert!((node.vertices[0].x() - 2.0).abs() < 1e-9, "v0.x = {}", node.vertices[0].x());
        assert!((node.vertices[1].x() - 3.0).abs() < 1e-9, "v1.x = {}", node.vertices[1].x());
        assert!((node.vertices[2].x() - 2.0).abs() < 1e-9, "v2.x = {}", node.vertices[2].x());
        assert!((node.vertices[2].y() - 1.0).abs() < 1e-9, "v2.y = {}", node.vertices[2].y());
        assert_eq!(mesh_triangle_count(&scene), 1);
    }

    #[test]
    fn vrml_material_diffuse() {
        let vrml = "#VRML V2.0 utf8\n\
            Shape {\n\
              appearance Appearance {\n\
                material Material { diffuseColor 0 0 1 }\n\
              }\n\
              geometry IndexedFaceSet {\n\
                coord Coordinate { point [ 0 0 0, 1 0 0, 0 1 0 ] }\n\
                coordIndex [ 0,1,2,-1 ]\n\
              }\n\
            }\n";
        let scene = parse_vrml_scene(vrml, "colored").unwrap();
        let mat = scene.nodes[0].material.as_ref().expect("diffuse material attached");
        let (r, g, b) = mat.diffuse;
        assert!(r.abs() < 1e-6 && g.abs() < 1e-6 && (b - 1.0).abs() < 1e-6, "diffuseColor 0 0 1");
    }

    #[test]
    fn convert_wrl_to_obj() {
        let mesh = crate::mesh::mesh_box((GpPnt::zero(), GpPnt::new(1.0, 1.0, 1.0)));
        let wrl_text = crate::vrml::write_vrml(&mesh, "box");
        let dir = temp_subdir("wrl2obj");
        let src = dir.join("m.wrl");
        let dst = dir.join("m.obj");
        std::fs::write(&src, wrl_text).unwrap();
        convert_mesh_format(&src.to_string_lossy(), &dst.to_string_lossy()).unwrap();
        let scene = read_obj_scene(&dst.to_string_lossy()).unwrap();
        assert_eq!(mesh_triangle_count(&scene), 12, "box has 12 tris after obj roundtrip");
        assert_eq!(mesh_vertex_count(&scene), 8);
    }

    #[test]
    fn material_lookup() {
        let scene = MeshScene {
            name: "s".into(),
            nodes: Vec::new(),
            materials: vec![Material {
                name: "brass".into(),
                diffuse: (0.8, 0.6, 0.2),
                ..Material::default()
            }],
        };
        assert!(material_from_name(&scene, "brass").is_some(), "known material found");
        let m = material_from_name(&scene, "brass").unwrap();
        assert!((m.diffuse.0 - 0.8).abs() < 1e-9);
        assert!(material_from_name(&scene, "nope").is_none(), "unknown material -> None");
    }

    #[test]
    fn mtl_file_roundtrip_uv() {
        let dir = temp_subdir("mtl_uv");
        std::fs::write(dir.join("m.mtl"), "newmtl blue\nKd 0 0 1\n").unwrap();
        let obj = "mtllib m.mtl\nusemtl blue\no part\n\
                   v 0 0 0\nv 1 0 0\nv 0 1 0\n\
                   vt 0.25 0.25\nvt 0.75 0.25\nvt 0.25 0.75\n\
                   f 1/1 2/2 3/3\n";
        let path = dir.join("m.obj");
        std::fs::write(&path, obj).unwrap();
        let scene = read_obj_scene(&path.to_string_lossy()).unwrap();
        let node = &scene.nodes[0];
        let mat = node.material.as_ref().expect("material preserved");
        assert!((mat.diffuse.2 - 1.0).abs() < 1e-6, "blue diffuse");
        let uv = node.uv.as_ref().expect("uv preserved alongside material");
        assert_eq!(uv.len(), 3);
        assert!((uv[1].x() - 0.75).abs() < 1e-9);
        assert_eq!(scene.materials.len(), 1);
    }

    #[test]
    fn no_material_is_none() {
        let obj = "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n";
        let scene = parse_obj_scene(obj, "plain").unwrap();
        assert!(scene.nodes[0].material.is_none(), "no usemtl -> no material");
        assert!(scene.nodes[0].uv.is_none(), "no vt -> no uv");
    }

    #[test]
    fn ppm_checkerboard_uv() {
        let scene = MeshScene {
            name: "s".into(),
            nodes: vec![MeshNode {
                name: "n".into(),
                vertices: vec![GpPnt::new(0., 0., 0.)],
                triangles: Vec::new(),
                transform: None,
                material: None,
                uv: Some(vec![GpPnt2d::new(0.1, 0.1)]),
                texture_path: None,
            }],
            materials: Vec::new(),
        };
        let ppm = scene_with_texture_to_ppm(&scene, 0, 2, 2).unwrap();
        assert!(ppm.starts_with("P3\n2 2\n255\n"), "P3 header with size");
        // Top-left quadrant (u<0.5, v<0.5) is red.
        assert!(ppm.contains("255 0 0\n"), "uv quadrant red present");
        // Without UVs the image falls back to material diffuse.
        let plain = MeshScene {
            name: "p".into(),
            nodes: vec![MeshNode {
                name: "n".into(),
                vertices: vec![GpPnt::new(0., 0., 0.)],
                triangles: Vec::new(),
                transform: None,
                material: None,
                uv: None,
                texture_path: None,
            }],
            materials: Vec::new(),
        };
        let ppm2 = scene_with_texture_to_ppm(&plain, 0, 1, 1).unwrap();
        assert!(ppm2.contains("179 179 179\n"), "default gray 0.7 -> 179");
        assert!(scene_with_texture_to_ppm(&plain, 5, 1, 1).is_err(), "bad index errors");
    }
}
