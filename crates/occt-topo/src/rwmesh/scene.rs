use super::prelude::*;
use super::*;

/// A surface material. `ambient`/`diffuse`/`specular` are RGB triples in
/// `[0,1]` (sRGB, as the `.mtl` text carries them); `shininess` is OCCT's
/// normalised specular exponent, `opacity` is `1.0` for fully opaque.
///
/// Source: `RWObj_Material` (`RWObj_Material.hxx:22-47`). That struct stores
/// `Transparency` where this one stores `opacity`, so `opacity = 1 -
/// transparency` throughout (`RWObj_MtlReader.cxx:210-240`). The three texture
/// fields are resolved to a path by `processTexturePath`
/// (`RWObj_MtlReader.cxx:270-306`), `None` when the material has no such map.
#[derive(Debug, Clone)]
pub struct Material {
    pub name: String,
    pub ambient: (f64, f64, f64),
    pub diffuse: (f64, f64, f64),
    pub specular: (f64, f64, f64),
    /// `Shininess` (`RWObj_Material.hxx:43`), `[0, 1]`; the MTL text carries
    /// `Ns = shininess * 1000` (`RWObj_ObjMaterialMap.cxx:84`).
    pub shininess: f64,
    pub opacity: f64,
    /// `DiffuseTexture` (`RWObj_Material.hxx:34`), from `map_Kd`.
    pub diffuse_texture: Option<String>,
    /// `SpecularTexture` (`RWObj_Material.hxx:35`), from `map_Ks`.
    pub specular_texture: Option<String>,
    /// `BumpTexture` (`RWObj_Material.hxx:36`), from `map_Bump`.
    pub bump_texture: Option<String>,
}

impl Default for Material {
    /// `RWObj_Material::RWObj_Material` (`RWObj_Material.hxx:45-46`):
    /// `AmbientColor(0.1)`, `DiffuseColor(0.8)`, `SpecularColor(0.2)`,
    /// `Shininess(1.0)`, `Transparency(0.0)`; the texture paths default empty
    /// (`TCollection_AsciiString` default-constructs empty).
    fn default() -> Self {
        Material {
            name: "default".to_string(),
            ambient: (0.1, 0.1, 0.1),
            diffuse: (0.8, 0.8, 0.8),
            specular: (0.2, 0.2, 0.2),
            shininess: 1.0,
            opacity: 1.0,
            diffuse_texture: None,
            specular_texture: None,
            bump_texture: None,
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
    /// `RWGltf_GltfLatePrimitiveArray::PrimitiveMode`
    /// (`RWGltf_GltfPrimitiveMode.hxx:20-30`). Formats that only produce
    /// triangulations (OBJ, PLY, STL, VRML) always report `Triangles`.
    pub mode: RwGltfPrimitiveMode,
    /// `RWGltf_GltfLatePrimitiveArray`'s edge array (`NbEdges`/`Edge`) for the
    /// `Lines` mode: the index buffer read in order, one node index per entry.
    /// The whole sequence becomes a single `Poly_Polygon3D`-backed edge
    /// (`RWGltf_GltfJsonParser.cxx:1985-1998`), not one edge per pair. Empty for
    /// the other modes.
    pub polyline: Vec<usize>,
    pub transform: Option<GpTrsf>,
    pub material: Option<Material>,
    /// The document-level material (`XCAFPrs_Style::SetMaterial`) the glTF
    /// reader bound to this node's shape (`RWGltf_GltfJsonParser.cxx:867-947`).
    /// `None` for formats and nodes that carry no such material.
    pub vis_material: Option<crate::xcaf_doc_vis_material::XcafDocVisMaterial>,
    pub uv: Option<Vec<GpPnt2d>>,
    /// `Poly_Triangulation::Normal` of the primitive's `NORMAL` accessor, after
    /// `RWGltf_TriangulationReader`'s coordinate conversion
    /// (`RWGltf_TriangulationReader.cxx:330-364`, `setNodeNormal`). The proxy
    /// shape's triangulation carries these back for export
    /// (`rwmesh::gltf_shape_reader`), and they are transformed by the node's
    /// rotation because this crate bakes the placement into the geometry.
    pub normals: Option<Vec<GpPnt>>,
    pub texture_path: Option<String>,
    /// `RWMesh_NodeAttributes::NamedData` for this node: the `extras` object
    /// parsed into `TDataStd_NamedData` (`RWGltf_GltfJsonParser.cxx:1621-1622`,
    /// `:2435`). `None` when neither the node nor its mesh defines `extras`.
    pub named_data: Option<crate::tdata_std_named_data::TDataStdNamedData>,
    /// The glTF primitive array this node's proxy shape carries in place of a
    /// triangulation — `RWGltf_GltfLatePrimitiveArray` put on the face by
    /// `BRep_Builder::MakeFace(face, theTriangulation)`
    /// (`RWGltf_GltfJsonParser.cxx:2000-2007`). `None` for formats that decode
    /// their data directly. `late.resolved()` is the triangulation currently in
    /// the slot; it stays `None` while `myToSkipLateDataLoading` deferred the
    /// array (`RWGltf_CafReader.cxx:376-388`).
    pub late: Option<std::sync::Arc<crate::rwmesh::GltfLatePrimitiveArray>>,
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
            mode: RwGltfPrimitiveMode::Triangles,
            polyline: Vec::new(),
            transform: None,
            material: None,
            vis_material: None,
            uv: None,
            normals: None,
            texture_path: None,
            named_data: None,
            late: None,
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
            mode: RwGltfPrimitiveMode::Triangles,
            polyline: Vec::new(),
            transform: None,
            material: None,
            vis_material: None,
            uv: None,
            normals: None,
            texture_path: None,
            named_data: None,
            late: None,
        }],
        materials: Vec::new(),
    })
}

// The glTF reader lives in gltf_file.rs.

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
            mode: RwGltfPrimitiveMode::Triangles,
            polyline: Vec::new(),
            transform: None,
            material: None,
            vis_material: None,
            uv: None,
            normals: None,
            texture_path: None,
            named_data: None,
            late: None,
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
        // `.glb` is the same reader: the container is sniffed from the `glTF`
        // magic (`RWGltf_CafReader.cxx:206-286`), not from the extension.
        "glb" => read_gltf_scene(src_path)?,
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
    // Fallback for a node with no material: `RWObj_Material`'s default diffuse
    // (`RWObj_Material.hxx:45-46`, carried by [`Material::default`]) so the
    // debug renderer and the reader agree on what "unstyled" looks like. This
    // whole function is port-only (OCCT has no PPM dump of a `MeshScene`).
    let unstyled = Material::default().diffuse;
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
                node.material.as_ref().map(|m| m.diffuse).unwrap_or(unstyled)
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
                // `RWObj_Tools::ReadName` (`RWObj_Tools.hxx:51-73`): the name is
                // the rest of the line, not the first token.
                let rest = line[tag.len()..].trim();
                cur_material = if rest.is_empty() { None } else { Some(rest.to_string()) };
            }
            "o" | "g" => {
                // `pushGroup` / `pushObject` use `ReadName` too, so an unquoted
                // group name may contain spaces (`g Open CASCADE STEP
                // translator 7.6 1` is one group).
                let rest = line[tag.len()..].trim();
                let gname = if rest.is_empty() { "mesh".to_string() } else { rest.to_string() };
                // `RWObj_SubMesh::Group` (`RWObj_Reader.cxx:628-646`): a new
                // group keeps the material of the sub-mesh being flushed; the
                // `f` branch below then replaces it with the material actually
                // in effect for the faces (`myLastFaceMaterial`,
                // `RWObj_TriangulationReader.cxx:39-55,61`).
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
                // `myLastFaceMaterial = theMesh.Material`
                // (`RWObj_TriangulationReader.cxx:61`): the group's bound
                // material is the one in effect for its **last** face, which is
                // why a `usemtl` after the `g` line still reaches the group.
                groups[current].material = cur_material.clone();
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
        // `RWObj_MtlReader::Read(theFolder, theFile)` (`cxx:117`): texture paths
        // in the library resolve against the folder the MTL lives in, not the
        // OBJ's folder.
        let mtl_dir = mtl_path.parent().unwrap_or_else(|| Path::new("."));
        materials.extend(parse_mtl(&mtl_content, mtl_dir));
    }

    let nodes = if groups.is_empty() {
        vec![MeshNode {
            name: "mesh".to_string(),
            vertices,
            triangles: Vec::new(),
            mode: RwGltfPrimitiveMode::Triangles,
            polyline: Vec::new(),
            transform: None,
            material: None,
            vis_material: None,
            uv: None,
            normals: None,
            texture_path: None,
            named_data: None,
            late: None,
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
            // `RWObj_CafReader::BindNamedShape` feeds `DiffuseTexture` into the
            // material's `XCAFDoc_VisMaterialCommon::DiffuseTexture`
            // (`RWObj_CafReader.cxx:63-65`); this crate keeps the same reference
            // on the node (`MeshNode::texture_path`), the slot the glTF reader
            // uses for its base-colour image.
            let texture_path = material.as_ref().and_then(|m| m.diffuse_texture.clone());
            out_nodes.push(MeshNode {
                name: g.name,
                vertices: local,
                triangles: tris,
                mode: RwGltfPrimitiveMode::Triangles,
                polyline: Vec::new(),
                transform: None,
                material,
                vis_material: None,
                uv: if has_uv { Some(uv) } else { None },
                normals: None,
                texture_path,
                named_data: None,
                late: None,
            });
        }
        out_nodes
    };

    Ok(MeshScene { name: name.to_string(), nodes, materials })
}

/// Parse a Wavefront `.mtl` material library: `newmtl <name>` starts a material,
/// `Kd`/`Ks`/`Ka` set the diffuse/specular/ambient RGB, `Ns` the specular
/// exponent, `d`/`Tr` the opacity, `map_Kd`/`map_Ks`/`map_Bump` the texture
/// paths.
///
/// `RWObj_MtlReader::Read` (`RWObj_MtlReader.cxx:135-262`) is the source: `Ka`
/// (`cxx:161-173`), `Kd` (`cxx:174-186`) and `Ks` (`cxx:187-199`) are the three
/// colours `Material` carries, `Ns` (`cxx:200-211`) is `Shininess =
/// min(Ns/1000, 1)`, `Tr` (`cxx:212-224`) stores transparency directly and `d`
/// (`cxx:225-236`) stores `1 - alpha`. The latter two collapse onto one field
/// here, so `Tr t` maps to `opacity = 1 - t` — the same value the writer derives
/// it from (`RWObj_ObjMaterialMap.cxx:118-121` writes `Tr = 1 - alpha`).
///
/// `folder` is `RWObj_MtlReader::Read`'s `theFolder`: texture names are resolved
/// against it (`processTexturePath`, `cxx:270-306`).
pub(super) fn parse_mtl(content: &str, folder: &Path) -> Vec<Material> {
    let mut mats: Vec<Material> = Vec::new();
    let mut cur: Option<Material> = None;
    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut parts = line.split_whitespace();
        let tag = parts.next().unwrap_or("");
        match tag {
            "newmtl" => {
                if let Some(m) = cur.take() {
                    mats.push(m);
                }
                cur = Some(Material {
                    name: parts.next().unwrap_or("").to_string(),
                    ..Material::default()
                });
            }
            "Ka" => {
                // `RWObj_MtlReader.cxx:161-173`: `AmbientColor`.
                if let Some(m) = &mut cur {
                    let c: Vec<f64> = parts.filter_map(|p| p.parse().ok()).collect();
                    if c.len() >= 3 && valid_color((c[0], c[1], c[2])) {
                        m.ambient = (c[0], c[1], c[2]);
                    }
                }
            }
            "Kd" => {
                if let Some(m) = &mut cur {
                    let c: Vec<f64> = parts.filter_map(|p| p.parse().ok()).collect();
                    if c.len() >= 3 && valid_color((c[0], c[1], c[2])) {
                        m.diffuse = (c[0], c[1], c[2]);
                    }
                }
            }
            "Ks" => {
                if let Some(m) = &mut cur {
                    let c: Vec<f64> = parts.filter_map(|p| p.parse().ok()).collect();
                    if c.len() >= 3 && valid_color((c[0], c[1], c[2])) {
                        m.specular = (c[0], c[1], c[2]);
                    }
                }
            }
            "Ns" => {
                // `RWObj_MtlReader.cxx:200-211`: `Shininess =
                // min(Ns / 1000.0, 1.0)`, rejected when negative.
                if let Some(m) = &mut cur {
                    if let Some(v) = parts.next().and_then(|p| p.parse::<f64>().ok()) {
                        if v >= 0.0 {
                            m.shininess = (v / 1000.0).min(1.0);
                        }
                    }
                }
            }
            "d" => {
                if let Some(m) = &mut cur {
                    if let Some(v) = parts.next().and_then(|p| p.parse::<f64>().ok()) {
                        // `cxx:228-240`: the value is a dissolve, `Transparency =
                        // 1 - alpha`, and `alpha < 0.01` is rejected.
                        if (0.0..=1.0).contains(&v) && v >= 0.01 {
                            m.opacity = v;
                        }
                    }
                }
            }
            "Tr" => {
                if let Some(m) = &mut cur {
                    if let Some(v) = parts.next().and_then(|p| p.parse::<f64>().ok()) {
                        // `cxx:212-224`: `Transparency = value`, rejected above
                        // 0.99 (fully transparent is not a colour).
                        if (0.0..=0.99).contains(&v) {
                            m.opacity = 1.0 - v;
                        }
                    }
                }
            }
            "map_Kd" | "map_Ks" | "map_Bump" => {
                // `cxx:242-267`: `RWObj_Tools::ReadName(aPos, aMat.*Texture)`
                // then `processTexturePath(aMat.*Texture, theFolder)`. A name
                // longer than `aLine[256]` cannot occur here because the reader
                // never truncates; an empty name is ignored like `ReadName`'s
                // `false` return.
                if let Some(m) = &mut cur {
                    let rest = line[tag.len()..].trim();
                    if !rest.is_empty() {
                        let resolved = process_texture_path(folder, rest);
                        match tag {
                            "map_Kd" => m.diffuse_texture = Some(resolved),
                            "map_Ks" => m.specular_texture = Some(resolved),
                            _ => m.bump_texture = Some(resolved),
                        }
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

/// `RWObj_MtlReader::processTexturePath` (`RWObj_MtlReader.cxx:270-306`): a
/// relative name is prefixed with the material library's folder; an absolute one
/// is kept as written, and when it does not exist the reader looks for the same
/// tail of the path next to the library (`findRelativePath`, `cxx:33-77`).
fn process_texture_path(folder: &Path, texture: &str) -> String {
    let p = Path::new(texture);
    if p.is_absolute() {
        if p.exists() {
            return texture.to_string();
        }
        // `findRelativePath(theAbsolutePath, theFolder, aRelativePath)` then
        // `theTexturePath = theFolder + aRelativePath` (`cxx:284-289`).
        if let Some(rel) = find_relative_path(p, folder) {
            return folder.join(rel).to_string_lossy().to_string();
        }
        return texture.to_string();
    }
    // `theTexturePath = theFolder + theTexturePath` (`cxx:302`).
    folder.join(p).to_string_lossy().to_string()
}

/// `findRelativePath` (`RWObj_MtlReader.cxx:33-77`): strip leading components
/// off the absolute path until the remaining tail exists under `folder`.
fn find_relative_path(absolute: &Path, folder: &Path) -> Option<std::path::PathBuf> {
    let mut tail: Vec<std::ffi::OsString> = Vec::new();
    for comp in absolute.components().rev() {
        if let std::path::Component::Normal(c) = comp {
            tail.insert(0, c.to_os_string());
            let candidate: std::path::PathBuf = tail.iter().collect();
            if folder.join(&candidate).exists() {
                return Some(candidate);
            }
        } else {
            break;
        }
    }
    None
}

/// `RWObj_MtlReader::validateColor` (`RWObj_MtlReader.cxx:336-350`): every
/// component within `[0, 1]`, otherwise the line is ignored with a warning.
/// The port has no message stream, so an invalid line simply leaves the field
/// at its default.
fn valid_color(c: (f64, f64, f64)) -> bool {
    (0.0..=1.0).contains(&c.0) && (0.0..=1.0).contains(&c.1) && (0.0..=1.0).contains(&c.2)
}

pub(super) fn find_material<'a>(materials: &'a [Material], name: &str) -> Option<&'a Material> {
    materials.iter().find(|m| m.name == name)
}
