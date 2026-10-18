// TEMPORARY probe: glTF 1.0 keyed maps, multi-scene, skip-empty, accessor stride.
use occt_topo::rwmesh::{read_gltf_scene, read_gltf_scene_with_options, GltfReadOptions};

fn b64(data: &[u8]) -> String {
    const T: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::new();
    for c in data.chunks(3) {
        let b = [c[0], *c.get(1).unwrap_or(&0), *c.get(2).unwrap_or(&0)];
        let n = (b[0] as u32) << 16 | (b[1] as u32) << 8 | b[2] as u32;
        for i in 0..4 {
            if i <= c.len() {
                out.push(T[((n >> (18 - 6 * i)) & 63) as usize] as char);
            } else {
                out.push('=');
            }
        }
    }
    out
}

fn tri_bytes() -> Vec<u8> {
    let mut bin = Vec::new();
    for p in [0.0f32, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0] {
        bin.extend_from_slice(&p.to_le_bytes());
    }
    for i in [0u16, 1, 2] {
        bin.extend_from_slice(&i.to_le_bytes());
    }
    bin
}

fn main() {
    std::fs::create_dir_all("_tmp/v1_probe").unwrap();
    let bin = tri_bytes();
    let uri = format!("data:application/octet-stream;base64,{}", b64(&bin));

    // glTF 1.0: every collection is an object keyed by id.
    let v1 = format!(
        r#"{{"asset":{{"version":"1.0"}},
"buffers":{{"buf":{{"byteLength":{},"uri":"{}"}}}},
"bufferViews":{{"bvpos":{{"buffer":"buf","byteOffset":0,"byteLength":36}},"bvidx":{{"buffer":"buf","byteOffset":36,"byteLength":6}}}},
"accessors":{{"accpos":{{"bufferView":"bvpos","componentType":5126,"count":3,"type":"VEC3"}},"accidx":{{"bufferView":"bvidx","componentType":5123,"count":3,"type":"SCALAR"}}}},
"meshes":{{"mesh0":{{"primitives":[{{"attributes":{{"POSITION":"accpos"}},"indices":"accidx"}}]}}}},
"nodes":{{"node0":{{"name":"v1node","mesh":"mesh0"}}}},
"scenes":{{"defaultScene":{{"nodes":["node0"]}}}},
"scene":"defaultScene"}}"#,
        bin.len(),
        uri
    );
    let p1 = "_tmp/v1_probe/v1.gltf";
    std::fs::write(p1, v1).unwrap();
    match read_gltf_scene(p1) {
        Ok(s) => println!(
            "v1: nodes={} name={:?} verts={} tris={}",
            s.nodes.len(),
            s.nodes.first().map(|n| n.name.clone()),
            s.nodes.first().map(|n| n.vertices.len()).unwrap_or(0),
            s.nodes.first().map(|n| n.triangles.len()).unwrap_or(0)
        ),
        Err(e) => println!("v1 ERR {e}"),
    }

    // glTF 2.0 with two scenes: scene 1 selected by default, both when asked.
    let bin_uri = format!("data:application/octet-stream;base64,{}", b64(&bin));
    let v2 = format!(
        r#"{{"asset":{{"version":"2.0"}},"scene":1,
"scenes":[{{"nodes":[0]}},{{"nodes":[1,2]}}],
"nodes":[{{"name":"a","mesh":0}},{{"name":"b","mesh":0}},{{"name":"empty"}}],
"meshes":[{{"primitives":[{{"attributes":{{"POSITION":0}},"indices":1}}]}}],
"buffers":[{{"byteLength":{},"uri":"{}"}}],
"bufferViews":[{{"buffer":0,"byteOffset":0,"byteLength":36}},{{"buffer":0,"byteOffset":36,"byteLength":6}}],
"accessors":[{{"bufferView":0,"componentType":5126,"count":3,"type":"VEC3"}},{{"bufferView":1,"componentType":5123,"count":3,"type":"SCALAR"}}]}}"#,
        bin.len(),
        bin_uri
    );
    let p2 = "_tmp/v1_probe/multi.gltf";
    std::fs::write(p2, v2).unwrap();
    let def = read_gltf_scene(p2).unwrap();
    println!("multi default: nodes={:?}", def.nodes.iter().map(|n| n.name.clone()).collect::<Vec<_>>());
    let mut all_opts = GltfReadOptions::occt_defaults();
    all_opts.to_load_all_scenes = true;
    let all = read_gltf_scene_with_options(p2, &all_opts).unwrap();
    println!("multi all: nodes={:?}", all.nodes.iter().map(|n| n.name.clone()).collect::<Vec<_>>());
    let mut no_skip = GltfReadOptions::occt_defaults();
    no_skip.to_skip_empty_nodes = false;
    all_opts.to_load_all_scenes = false;
    let _ = read_gltf_scene_with_options(p2, &no_skip).unwrap();
    println!("skip-empty keeps empty node: {}", def.nodes.len());

    // Interleaved POSITION with byteStride 16 (3 floats + 1 pad).
    let mut stride_bin = Vec::new();
    for p in [0.0f32, 0.0, 0.0, 9.0, 2.0, 0.0, 0.0, 9.0, 0.0, 2.0, 0.0, 9.0] {
        stride_bin.extend_from_slice(&p.to_le_bytes());
    }
    let suri = format!("data:application/octet-stream;base64,{}", b64(&stride_bin));
    let v2s = format!(
        r#"{{"asset":{{"version":"2.0"}},"scene":0,
"scenes":[{{"nodes":[0]}}],"nodes":[{{"name":"s","mesh":0}}],
"meshes":[{{"primitives":[{{"attributes":{{"POSITION":0}}}}]}}],
"buffers":[{{"byteLength":{},"uri":"{}"}}],
"bufferViews":[{{"buffer":0,"byteOffset":0,"byteLength":48}}],
"accessors":[{{"bufferView":0,"componentType":5126,"count":3,"type":"VEC3","byteStride":16}}]}}"#,
        stride_bin.len(),
        suri
    );
    let p3 = "_tmp/v1_probe/stride.gltf";
    std::fs::write(p3, v2s).unwrap();
    let s = read_gltf_scene(p3).unwrap();
    println!(
        "stride: v0=({:.1},{:.1},{:.1}) v2=({:.1},{:.1},{:.1}) tris={}",
        s.nodes[0].vertices[0].x(),
        s.nodes[0].vertices[0].y(),
        s.nodes[0].vertices[0].z(),
        s.nodes[0].vertices[2].x(),
        s.nodes[0].vertices[2].y(),
        s.nodes[0].vertices[2].z(),
        s.nodes[0].triangles.len()
    );

    // Missing default scene must be an error.
    let bad = r#"{"asset":{"version":"2.0"},"scene":7,"scenes":[{"nodes":[0]}],"nodes":[{"name":"a"}],"meshes":[],"buffers":[{"byteLength":0,"uri":"data:application/octet-stream;base64,"}]}"#;
    let p4 = "_tmp/v1_probe/bad.gltf";
    std::fs::write(p4, bad).unwrap();
    println!("bad scene: {:?}", read_gltf_scene(p4).err());

    // GLB: the buffer carries no uri and is the BIN chunk; its bytes are
    // addressed by the file offset (`RWGltf_GltfJsonParser.cxx:2286-2306`).
    let glb_json = format!(
        r#"{{"asset":{{"version":"2.0"}},"scene":0,"scenes":[{{"nodes":[0]}}],"nodes":[{{"name":"glb","mesh":0}}],"meshes":[{{"primitives":[{{"attributes":{{"POSITION":0}},"indices":1}}]}}],"buffers":[{{"byteLength":{}}}],"bufferViews":[{{"buffer":0,"byteOffset":0,"byteLength":36}},{{"buffer":0,"byteOffset":36,"byteLength":6}}],"accessors":[{{"bufferView":0,"componentType":5126,"count":3,"type":"VEC3"}},{{"bufferView":1,"componentType":5123,"count":3,"type":"SCALAR"}}]}}"#,
        bin.len()
    );
    let mut glb = Vec::new();
    let mut json = glb_json.into_bytes();
    while json.len() % 4 != 0 {
        json.push(b' ');
    }
    let mut body = bin.clone();
    while body.len() % 4 != 0 {
        body.push(0);
    }
    let total = 12 + 8 + json.len() + 8 + body.len();
    glb.extend_from_slice(b"glTF");
    glb.extend_from_slice(&2u32.to_le_bytes());
    glb.extend_from_slice(&(total as u32).to_le_bytes());
    glb.extend_from_slice(&(json.len() as u32).to_le_bytes());
    glb.extend_from_slice(&0x4E4F534Au32.to_le_bytes());
    glb.extend_from_slice(&json);
    glb.extend_from_slice(&(body.len() as u32).to_le_bytes());
    glb.extend_from_slice(&0x004E4942u32.to_le_bytes());
    glb.extend_from_slice(&body);
    let p5 = "_tmp/v1_probe/tri.glb";
    std::fs::write(p5, &glb).unwrap();
    match read_gltf_scene(p5) {
        Ok(s) => println!(
            "glb: nodes={} name={:?} verts={} tris={}",
            s.nodes.len(),
            s.nodes.first().map(|n| n.name.clone()),
            s.nodes.first().map(|n| n.vertices.len()).unwrap_or(0),
            s.nodes.first().map(|n| n.triangles.len()).unwrap_or(0)
        ),
        Err(e) => println!("glb ERR {e}"),
    }

    // Late-data branches on the same GLB: parallel decode and proxy-only read.
    let mut par_opts = GltfReadOptions::occt_defaults();
    par_opts.to_parallel = true;
    match read_gltf_scene_with_options(p5, &par_opts) {
        Ok(s) => println!("glb parallel: tris={}", s.nodes[0].triangles.len()),
        Err(e) => println!("glb parallel ERR {e}"),
    }
    let mut proxy_opts = GltfReadOptions::occt_defaults();
    proxy_opts.to_skip_late_data_loading = true;
    match read_gltf_scene_with_options(p5, &proxy_opts) {
        Ok(s) => println!(
            "glb proxy: nodes={} verts={} late={}",
            s.nodes.len(),
            s.nodes[0].vertices.len(),
            s.nodes[0].late.is_some()
        ),
        Err(e) => println!("glb proxy ERR {e}"),
    }

    // LINES / POINTS primitives (`RWGltf_CafReader.cxx:1972-1998`): an edge over
    // a polygon and a compound of vertices, not a face.
    for (mode, tag) in [(1u32, "lines"), (2, "lineloop"), (3, "linestrip"), (0, "points")] {
        let body = format!(
            r#"{{"asset":{{"version":"2.0"}},"scene":0,"scenes":[{{"nodes":[0]}}],"nodes":[{{"name":"{tag}","mesh":0}}],"meshes":[{{"primitives":[{{"attributes":{{"POSITION":0}},"indices":1,"mode":{mode}}}]}}],"buffers":[{{"byteLength":{},"uri":"{}"}}],"bufferViews":[{{"buffer":0,"byteOffset":0,"byteLength":36}},{{"buffer":0,"byteOffset":36,"byteLength":6}}],"accessors":[{{"bufferView":0,"componentType":5126,"count":3,"type":"VEC3"}},{{"bufferView":1,"componentType":5123,"count":3,"type":"SCALAR"}}]}}"#,
            bin.len(),
            bin_uri
        );
        let path = format!("_tmp/v1_probe/{tag}.gltf");
        std::fs::write(&path, body).unwrap();
        match read_gltf_scene(&path) {
            Ok(s) => println!(
                "{tag}: nodes={} verts={} tris={} polyline={}",
                s.nodes.len(),
                s.nodes.first().map(|n| n.vertices.len()).unwrap_or(0),
                s.nodes.first().map(|n| n.triangles.len()).unwrap_or(0),
                s.nodes.first().map(|n| n.polyline.len()).unwrap_or(0)
            ),
            Err(e) => println!("{tag} ERR {e}"),
        }
    }

    // Non-uniform node scale (`cxx:596-634`): the triangulation nodes are
    // scaled, so the world extent doubles on X only.
    if let Ok(s) = read_gltf_scene("_tmp/zz_probe_scale/s.gltf") {
        let xs: Vec<f64> = s.nodes[0].vertices.iter().map(|p| p.x()).collect();
        let ys: Vec<f64> = s.nodes[0].vertices.iter().map(|p| p.y()).collect();
        println!(
            "scale: x=[{:.1},{:.1}] y=[{:.1},{:.1}]",
            xs.iter().cloned().fold(f64::MAX, f64::min),
            xs.iter().cloned().fold(f64::MIN, f64::max),
            ys.iter().cloned().fold(f64::MAX, f64::min),
            ys.iter().cloned().fold(f64::MIN, f64::max)
        );
    }

    // `extras` -> `TDataStd_NamedData` (`setShapeNamedData`, cxx:658).
    let extras_json = format!(
        r#"{{"asset":{{"version":"2.0"}},"scene":0,"scenes":[{{"nodes":[0]}}],"nodes":[{{"name":"x","mesh":0,"extras":{{"a":2,"s":"v","arr":[1,2]}}}}],"meshes":[{{"primitives":[{{"attributes":{{"POSITION":0}}}}],"extras":{{"m":1.5}}}}],"buffers":[{{"byteLength":36,"uri":"data:application/octet-stream;base64,AAAAAAAAAAAAAAAAAACAPwAAAAAAAAAAAAAAAAAAgD8AAAAA"}}],"bufferViews":[{{"buffer":0,"byteLength":36}}],"accessors":[{{"bufferView":0,"componentType":5126,"count":3,"type":"VEC3"}}]}}"#
    );
    let p6 = "_tmp/v1_probe/extras.gltf";
    std::fs::write(p6, extras_json).unwrap();
    let s = read_gltf_scene(p6).unwrap();
    match &s.nodes[0].named_data {
        Some(nd) => println!(
            "extras: a={:?} s={:?} arr={:?} m={:?}",
            nd.integer("a"),
            nd.string("s"),
            nd.array_of_integers("arr"),
            nd.real("m")
        ),
        None => println!("extras: none"),
    }
    let model = occt_topo::rwmesh::read_gltf_model(p6).unwrap();
    println!(
        "extras model: {}",
        model
            .shapes
            .first()
            .and_then(|s| s.named_data.as_ref())
            .map(|nd| format!("a={:?}", nd.integer("a")))
            .unwrap_or_else(|| "none".to_string())
    );

    // Document labels (`fillDocument` / `addShapeIntoDoc`): the Product/Instance
    // split, assemblies and sub-shape labels.
    let dump = |tag: &str, m: &occt_topo::model::BRepModel| {
        println!(
            "{tag}: shapes={} labels={} roots={:?} mats={}",
            m.shapes.len(),
            m.labels.len(),
            m.roots,
            m.materials.len()
        );
        for (i, l) in m.labels.iter().enumerate() {
            println!(
                "  [{i}] {:?} name={:?} ref={:?} parent={:?} children={:?} shape={:?} color={:?} curv={:?} mat={:?} nd={:?}",
                l.role,
                l.name,
                l.referred,
                l.parent,
                l.children,
                l.shape.as_ref().map(|s| s.shape_type()),
                l.color.is_some(),
                l.color_curv.is_some(),
                l.material,
                l.named_data.as_ref().map(|n| n.integer("a"))
            );
        }
    };
    match occt_topo::rwmesh::read_gltf_model_with_options(p2, &all_opts) {
        Ok(m) => dump("doc all scenes", &m),
        Err(e) => println!("doc all scenes ERR {e}"),
    }
    match occt_topo::rwmesh::read_gltf_model(p6) {
        Ok(m) => dump("doc extras", &m),
        Err(e) => println!("doc extras ERR {e}"),
    }

    // Hierarchy: an assembly node instancing a named mesh.
    let asm = format!(
        r#"{{"asset":{{"version":"2.0"}},"scene":0,"scenes":[{{"nodes":[0]}}],"nodes":[{{"name":"Root","children":[1]}},{{"name":"Inst","mesh":0,"translation":[5,0,0]}}],"meshes":[{{"name":"Box","primitives":[{{"attributes":{{"POSITION":0}},"indices":1,"material":0}}]}}],"materials":[{{"name":"red","pbrMetallicRoughness":{{"baseColorFactor":[1,0,0,1]}}}}],"buffers":[{{"byteLength":{},"uri":"{}"}}],"bufferViews":[{{"buffer":0,"byteOffset":0,"byteLength":36}},{{"buffer":0,"byteOffset":36,"byteLength":6}}],"accessors":[{{"bufferView":0,"componentType":5126,"count":3,"type":"VEC3"}},{{"bufferView":1,"componentType":5123,"count":3,"type":"SCALAR"}}]}}"#,
        bin.len(),
        bin_uri
    );
    let p7 = "_tmp/v1_probe/asm.gltf";
    std::fs::write(p7, asm).unwrap();
    match occt_topo::rwmesh::read_gltf_model(p7) {
        Ok(m) => dump("doc asm", &m),
        Err(e) => println!("doc asm ERR {e}"),
    }

    // Input coordinate system: glTF is Y-up, OCCT is Z-up
    // (`RWMesh_CoordinateSystemConverter`, `RWGltf_CafReader.cxx:182-183`).
    let mut cs_opts = GltfReadOptions::occt_defaults();
    cs_opts
        .converter
        .set_output_coordinate_system(occt_topo::rwmesh::RwMeshCoordinateSystem::ZUP);
    match read_gltf_scene_with_options("_tmp/zz_probe_scale/s.gltf", &cs_opts) {
        Ok(s) => {
            let p = &s.nodes[0].vertices[2];
            println!("cs: v2=({:.1},{:.1},{:.1})", p.x(), p.y(), p.z());
        }
        Err(e) => println!("cs ERR {e}"),
    }

    // Document-style OBJ export (`RWObj_CafWriter::Perform` over
    // `XCAFPrs_DocumentExplorer`): group names from the label hierarchy,
    // per-face styles from `XCAFPrs::CollectStyleSettings`, triangles/normals
    // read back from the stored `Poly_Triangulation`.
    if let Ok(m) = occt_topo::rwmesh::read_gltf_model(p7) {
        if let Some(whole) = occt_topo::brep_exchange::document_shape(&m) {
            let (obj, mtl) = occt_topo::brep_exchange::model_to_obj_caf(&m, &whole, "asm.obj", 0.1);
            println!("doc obj:\n{obj}");
            println!("doc mtl:\n{mtl}");
        }
    }

    // PLY document export (`RWPly_CafWriter::Perform`). Both CAF writers walk
    // the same node list, so the vertex/face counts must equal the OCCT OBJ
    // baselines (`data/occ-*.obj`, `v` / `f` records).
    for stem in ["Cube", "Cone"] {
        let model = match occt_topo::step::read_step_file(&format!("data/{stem}.step")) {
            Ok(m) => m,
            Err(e) => {
                println!("ply {stem}: STEP ERR {e}");
                continue;
            }
        };
        let Some(whole) = occt_topo::brep_exchange::document_shape(&model) else { continue };
        let opts = occt_topo::rwmesh::PlyCafWriterOptions::default();
        let ply = match occt_topo::brep_exchange::model_to_ply_caf(
            &model,
            &whole,
            "probe.ply",
            0.1,
            &opts,
        ) {
            Ok(text) => text,
            Err(e) => {
                println!("ply {stem}: ERR {e}");
                continue;
            }
        };
        let header = |prefix: &str| -> usize {
            ply.lines()
                .find_map(|l| l.strip_prefix(prefix))
                .and_then(|s| s.trim().parse().ok())
                .unwrap_or(0)
        };
        let header_v = header("element vertex ");
        let header_f = header("element face ");
        let body = ply.split_once("end_header\n").map(|(_, b)| b).unwrap_or("");
        let body_v = body.lines().filter(|l| !l.starts_with("3 ")).count();
        let body_f = body.lines().filter(|l| l.starts_with("3 ")).count();
        let base = std::fs::read_to_string(format!("data/occ-{}.obj", stem.to_lowercase()))
            .unwrap_or_default();
        let base_v = base.lines().filter(|l| l.starts_with("v ")).count();
        let base_f = base.lines().filter(|l| l.starts_with("f ")).count();
        println!(
            "ply {stem}: header v={header_v} f={header_f} body v={body_v} f={body_f} obj v={base_v} f={base_f} match={}",
            header_v == base_v && header_f == base_f && body_v == base_v && body_f == base_f
        );
    }

    // glTF document export (`RWGltf_CafWriter::Perform`): structure and
    // numbers against the OCCT baselines `data/occ-*.gltf` + `.bin`.
    for stem in ["Cube", "Cone", "linkrods"] {
        let lower = stem.to_lowercase();
        let model = match occt_topo::step::read_step_file(&format!("data/{stem}.step")) {
            Ok(m) => m,
            Err(e) => {
                println!("gltf {stem}: STEP ERR {e}");
                continue;
            }
        };
        let Some(whole) = occt_topo::brep_exchange::document_shape(&model) else { continue };
        let path = format!("data/occ-{lower}.gltf");
        let opts = occt_topo::rwmesh::GltfCafWriterOptions::default();
        let out = match occt_topo::brep_exchange::model_to_gltf_caf(&model, &whole, &path, 0.1, &opts)
        {
            Ok(out) => out,
            Err(e) => {
                println!("gltf {stem}: ERR {e}");
                continue;
            }
        };
        let base_json = std::fs::read(&path).unwrap_or_default();
        let base_bin = std::fs::read(format!("data/occ-{lower}.bin")).unwrap_or_default();
        std::fs::create_dir_all("_tmp/gltf").unwrap();
        std::fs::write(format!("_tmp/gltf/{lower}.gltf"), &out.document).unwrap();
        std::fs::write(format!("_tmp/gltf/{lower}.bin"), &out.bin).unwrap();
        let doc: &[u8] = &out.document;
        let first_diff = doc
            .iter()
            .zip(&base_json)
            .position(|(a, b)| a != b)
            .unwrap_or_else(|| doc.len().min(base_json.len()));
        let window = |data: &[u8]| -> String {
            let start = first_diff.saturating_sub(40);
            let end = (first_diff + 80).min(data.len());
            String::from_utf8_lossy(&data[start..end]).replace('\n', "\\n")
        };
        println!(
            "gltf {stem}: json mine={} base={} equal={} first_diff={}",
            doc.len(),
            base_json.len(),
            doc == base_json.as_slice(),
            first_diff
        );
        if doc != base_json.as_slice() {
            println!("  mine: {}", window(doc));
            println!("  base: {}", window(&base_json));
        }
        let bin_diff = out
            .bin
            .iter()
            .zip(&base_bin)
            .position(|(a, b)| a != b)
            .unwrap_or_else(|| out.bin.len().min(base_bin.len()));
        println!(
            "gltf {stem}: bin mine={} base={} equal={} first_diff={}",
            out.bin.len(),
            base_bin.len(),
            out.bin == base_bin,
            bin_diff
        );
        if out.bin != base_bin {
            let diffs: Vec<usize> = out
                .bin
                .chunks_exact(4)
                .zip(base_bin.chunks_exact(4))
                .enumerate()
                .filter(|(_, (a, b))| a != b)
                .map(|(i, _)| i)
                .collect();
            println!("  diff words: {}", diffs.len());
            for &word in diffs.iter().take(8) {
                let a = &out.bin[word * 4..word * 4 + 4];
                let b = &base_bin[word * 4..word * 4 + 4];
                println!(
                    "  w{word} mine={:?}({:?}) base={:?}({:?})",
                    f32::from_le_bytes([a[0], a[1], a[2], a[3]]),
                    a,
                    f32::from_le_bytes([b[0], b[1], b[2], b[3]]),
                    b
                );
            }
        }
        // Forced UV export pass (PROBE_FORCE_UV=1): the same writer run with
        // `-forceUVExport`, dumped next to the reference glTF so the per-node
        // surface parameters can be compared against OCCT's.
        if std::env::var("PROBE_FORCE_UV").map(|v| v == "1").unwrap_or(false) {
            let mut uv_opts = occt_topo::rwmesh::GltfCafWriterOptions::default();
            uv_opts.is_forced_uv_export = true;
            let uv_path = format!("_tmp/gltf/{lower}_uv.gltf");
            match occt_topo::brep_exchange::model_to_gltf_caf(&model, &whole, &uv_path, 0.1, &uv_opts)
            {
                Ok(out) => {
                    std::fs::create_dir_all("_tmp/gltf").unwrap();
                    std::fs::write(&uv_path, &out.document).unwrap();
                    std::fs::write(format!("_tmp/gltf/{lower}_uv.bin"), &out.bin).unwrap();
                    println!("gltf {stem}: forced uv json={} bin={}", out.document.len(), out.bin.len());
                }
                Err(e) => println!("gltf {stem}: forced uv ERR {e}"),
            }
        }
    }
}
