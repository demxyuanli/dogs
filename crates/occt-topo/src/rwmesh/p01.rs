use super::prelude::*;
use super::*;

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

pub(super) fn write_scene(scene: &MeshScene, path: &str) -> Result<(), String> {
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
pub(super) fn parse_obj_scene(content: &str, name: &str) -> Result<MeshScene, String> {
    parse_obj_scene_inner(content, name, None)
}

/// One OBJ group: a named collection of polygon faces; each corner is a
/// `(vertex_index, texcoord_index)` pair. `material` is the `usemtl` name in
/// effect when the group started.
pub(super) struct ObjGroup {
    pub(super) name: String,
    pub(super) faces: Vec<Vec<(usize, Option<usize>)>>,
    pub(super) material: Option<String>,
}

pub(super) fn parse_obj_scene_inner(
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
pub(super) fn parse_mtl(content: &str) -> Vec<Material> {
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

pub(super) fn find_material<'a>(materials: &'a [Material], name: &str) -> Option<&'a Material> {
    materials.iter().find(|m| m.name == name)
}
