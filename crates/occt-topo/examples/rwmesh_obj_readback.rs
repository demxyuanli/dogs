//! Read-back check for the OBJ + MTL pair the CAF exporter writes.
//!
//! Pipeline: `examples/export_data_obj.rs` (or `step_obj_gates`) writes
//! `data/output/<stem>.obj` + `<stem>.mtl`; this example reads the same files
//! back through `rwmesh::read_obj_scene` and reports, per node, whether the
//! `usemtl` material actually landed on the node. That closes the loop the
//! writer-side gates leave open: `step_to_obj` only round-trips vertices/faces
//! through `occt_core::io::obj`, so `mtllib` -> library -> node material was
//! never read back from the artifacts the exporter produces.
//!
//! Args: OBJ files, or directories scanned for `*.obj`. No argument defaults to
//! `data/output/gates` (the `step_obj_gates` dump). The `.mtl` referenced by
//! `mtllib` must sit next to the OBJ, which is what both producers write.
//!
//! Checks (exit code 1 on any failure):
//!   * `usemtl <name>` resolves to a material in the loaded library,
//!   * each node's material equals the **last** `usemtl` of its OBJ group —
//!     `RWObj_TriangulationReader::addMesh` binds the group's material from
//!     `myLastFaceMaterial` (`RWObj_TriangulationReader.cxx:46-50,61`), and the
//!     exporter writes `g <node>` *before* `usemtl` (`RWObj_CafWriter.cxx:294-308`),
//!     so a reader that binds the material at `g` time attaches nothing,
//!   * the attached material's `Kd` equals the `Kd` line of that material in the
//!     `.mtl` (the file is the oracle for what was written),
//!   * every triangle index is inside the node's vertex list.
//!
//! Usage:
//!   cargo run --example rwmesh_obj_readback
//!   cargo run --example rwmesh_obj_readback -- data/output/gates/Cone.obj
//!   cargo run --example rwmesh_obj_readback -- data/output

use std::path::{Path, PathBuf};
use std::process::ExitCode;

use occt_topo::rwmesh::{material_from_name, read_obj_scene, Material, MeshScene};

/// `data/output/gates` under the repository root.
fn default_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..").join("..").join("data").join("output").join("gates")
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let targets: Vec<PathBuf> = if args.is_empty() {
        vec![default_dir()]
    } else {
        args.iter().map(PathBuf::from).collect()
    };

    let mut files: Vec<PathBuf> = Vec::new();
    for t in &targets {
        if t.is_dir() {
            let mut found: Vec<PathBuf> = std::fs::read_dir(t)
                .unwrap_or_else(|e| panic!("read {}: {e}", t.display()))
                .filter_map(|e| e.ok())
                .map(|e| e.path())
                .filter(|p| {
                    p.extension().and_then(|e| e.to_str()).map(|e| e.eq_ignore_ascii_case("obj")).unwrap_or(false)
                })
                .collect();
            found.sort();
            files.extend(found);
        } else if t.is_file() {
            files.push(t.clone());
        } else {
            panic!("not a file or directory: {}", t.display());
        }
    }
    files.sort();
    files.dedup();

    if files.is_empty() {
        eprintln!("no .obj found under {:?}", targets);
        return ExitCode::FAILURE;
    }

    let mut checks = Checks::default();
    for obj in &files {
        checks.objects += 1;
        if let Err(e) = check_one(obj, &mut checks) {
            checks.failures += 1;
            eprintln!("err {}: {e}", obj.display());
        }
    }

    println!(
        "summary: objs={} nodes={} library={} attached={} unresolved={} kd_mismatch={} aspect_mismatch={} bad_index={} failures={}",
        checks.objects,
        checks.nodes,
        checks.library,
        checks.attached,
        checks.unresolved,
        checks.kd_mismatch,
        checks.aspect_mismatch,
        checks.bad_index,
        checks.failures
    );
    if checks.failures == 0 {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}

#[derive(Default)]
struct Checks {
    objects: usize,
    nodes: usize,
    library: usize,
    attached: usize,
    unresolved: usize,
    kd_mismatch: usize,
    /// Materials whose `Ka`/`Ks`/`Ns`/`Tr` did not survive the round trip.
    aspect_mismatch: usize,
    bad_index: usize,
    failures: usize,
}

/// One `newmtl` block as written, field by field. `None` means the line was not
/// emitted (`RWObj_ObjMaterialMap::DefineMaterial` writes the whole block or
/// nothing at all, `RWObj_ObjMaterialMap.cxx:78-121`).
#[derive(Default, Clone)]
struct WrittenMaterial {
    ka: Option<(f64, f64, f64)>,
    kd: Option<(f64, f64, f64)>,
    ks: Option<(f64, f64, f64)>,
    ns: Option<f64>,
    tr: Option<f64>,
    /// `map_Kd` / `map_Ks` / `map_Bump` value as written, unresolved.
    map_kd: Option<String>,
    map_ks: Option<String>,
    map_bump: Option<String>,
}

/// `mtllib` / `usemtl` / `g` records of an OBJ, plus every `newmtl` block of the
/// `.mtl` its `mtllib` names. Read here from the files themselves: the point of
/// the check is to compare what the reader produced against what was written, so
/// reusing `read_obj_scene`/`parse_mtl` for the expected side would compare the
/// code with itself.
#[derive(Default)]
struct Written {
    /// `usemtl` names in file order.
    used: Vec<String>,
    /// Group name -> the `usemtl` in effect for its last `f` line (`None` when
    /// the group has no faces or no material). This is `myLastFaceMaterial`.
    group_material: Vec<(String, Option<String>)>,
    /// Group name -> number of `f` lines.
    group_faces: Vec<(String, usize)>,
    /// `newmtl` name -> the block written for it.
    mats: Vec<(String, WrittenMaterial)>,
}

fn written(path: &Path) -> Result<Written, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("read {}: {e}", path.display()))?;
    let mut w = Written::default();
    let mut libs: Vec<String> = Vec::new();
    let mut group = String::from("mesh");
    let mut cur_mat: Option<String> = None;
    for line in text.lines() {
        let line = line.trim();
        if let Some(rest) = line.strip_prefix("usemtl ") {
            let name = rest.trim().to_string();
            w.used.push(name.clone());
            cur_mat = Some(name);
        } else if let Some(rest) = line.strip_prefix("mtllib ") {
            for f in rest.split_whitespace() {
                libs.push(f.to_string());
            }
        } else if line.starts_with("g ") || line.starts_with("o ") {
            group = line[2..].trim().to_string();
            set_group_material(&mut w, &group, cur_mat.clone());
        } else if line.starts_with("f ") {
            set_group_material(&mut w, &group, cur_mat.clone());
            match w.group_faces.iter_mut().find(|(n, _)| *n == group) {
                Some((_, n)) => *n += 1,
                None => w.group_faces.push((group.clone(), 1)),
            }
        }
    }
    let dir = path.parent().unwrap_or_else(|| Path::new("."));
    for lib in &libs {
        let mtl_path = dir.join(lib);
        let mtl = std::fs::read_to_string(&mtl_path)
            .map_err(|e| format!("read {}: {e}", mtl_path.display()))?;
        let mut cur: Option<String> = None;
        let mut block = WrittenMaterial::default();
        let flush = |name: Option<String>, block: &mut WrittenMaterial, w: &mut Written| {
            if let Some(name) = name {
                w.mats.push((name, std::mem::take(block)));
            }
        };
        for line in mtl.lines() {
            let line = line.trim();
            if let Some(rest) = line.strip_prefix("newmtl ") {
                flush(cur.take(), &mut block, &mut w);
                cur = Some(rest.trim().to_string());
                block = WrittenMaterial::default();
            } else if let Some(rest) = line.strip_prefix("Ka ") {
                block.ka = vec3(rest);
            } else if let Some(rest) = line.strip_prefix("Kd ") {
                block.kd = vec3(rest);
            } else if let Some(rest) = line.strip_prefix("Ks ") {
                block.ks = vec3(rest);
            } else if let Some(rest) = line.strip_prefix("Ns ") {
                block.ns = rest.trim().parse().ok();
            } else if let Some(rest) = line.strip_prefix("Tr ") {
                block.tr = rest.trim().parse().ok();
            } else if let Some(rest) = line.strip_prefix("map_Kd ") {
                block.map_kd = resolve_texture(mtl_path.parent(), rest);
            } else if let Some(rest) = line.strip_prefix("map_Ks ") {
                block.map_ks = resolve_texture(mtl_path.parent(), rest);
            } else if let Some(rest) = line.strip_prefix("map_Bump ") {
                block.map_bump = resolve_texture(mtl_path.parent(), rest);
            }
        }
        flush(cur.take(), &mut block, &mut w);
    }
    Ok(w)
}

fn vec3(rest: &str) -> Option<(f64, f64, f64)> {
    let c: Vec<f64> = rest.split_whitespace().filter_map(|p| p.parse().ok()).collect();
    (c.len() >= 3).then(|| (c[0], c[1], c[2]))
}

/// Expected value of a `map_*` line after `processTexturePath`
/// (`RWObj_MtlReader.cxx:270-306`): relative names are relative to the `.mtl`
/// folder, absolute names stay as written. Re-derived here rather than through
/// `parse_mtl` so the check does not compare the reader with itself.
fn resolve_texture(folder: Option<&Path>, rest: &str) -> Option<String> {
    let name = rest.trim();
    if name.is_empty() {
        return None;
    }
    let folder = folder.unwrap_or_else(|| Path::new("."));
    let p = Path::new(name);
    if p.is_absolute() {
        if p.exists() {
            return Some(name.to_string());
        }
        // `findRelativePath` (`RWObj_MtlReader.cxx:33-77`): retry the tail of the
        // absolute path under the library folder before giving up.
        let mut tail: Vec<std::ffi::OsString> = Vec::new();
        for comp in p.components().rev() {
            let std::path::Component::Normal(c) = comp else { break };
            tail.insert(0, c.to_os_string());
            let candidate: std::path::PathBuf = tail.iter().collect();
            if folder.join(&candidate).exists() {
                return Some(folder.join(&candidate).to_string_lossy().to_string());
            }
        }
        return Some(name.to_string());
    }
    Some(folder.join(p).to_string_lossy().to_string())
}

fn set_group_material(w: &mut Written, group: &str, material: Option<String>) {
    match w.group_material.iter_mut().find(|(n, _)| n == group) {
        Some(slot) => slot.1 = material,
        None => w.group_material.push((group.to_string(), material)),
    }
}

fn check_one(obj: &Path, c: &mut Checks) -> Result<(), String> {
    let scene = read_obj_scene(&obj.to_string_lossy())?;
    let w = written(obj)?;

    println!(
        "== {}\nscene {}  nodes={}  library={}  mtllib={}",
        obj.file_name().unwrap_or_default().to_string_lossy(),
        scene.name,
        scene.nodes.len(),
        scene.materials.len(),
        w.mats.len()
    );
    c.library += scene.materials.len();

    for name in &w.used {
        if material_from_name(&scene, name).is_none() {
            c.unresolved += 1;
            c.failures += 1;
            eprintln!("   usemtl {name}: not in library ({:?})", names(&scene));
        }
    }

    for node in &scene.nodes {
        c.nodes += 1;
        if let Some(t) = node.triangles.iter().find(|t| t.0 >= node.vertices.len() || t.1 >= node.vertices.len() || t.2 >= node.vertices.len()) {
            c.bad_index += 1;
            c.failures += 1;
            eprintln!("   node {}: triangle {t:?} out of range (v={})", node.name, node.vertices.len());
        }
        let expected = w.group_material.iter().find(|(n, _)| *n == node.name).and_then(|(_, m)| m.clone());
        if expected.is_none() && w.group_material.iter().all(|(n, _)| *n != node.name) {
            c.failures += 1;
            eprintln!("   node {}: no group of that name in the OBJ", node.name);
        }
        print_node(
            &node.name,
            node.vertices.len(),
            node.triangles.len(),
            node.material.as_ref(),
            node.uv.is_some(),
            node.texture_path.as_deref(),
        );

        match (&node.material, expected) {
            (Some(m), Some(want)) if m.name == want => {
                c.attached += 1;
                check_material(m, &w, node, c);
            }
            (Some(m), Some(want)) => {
                c.failures += 1;
                eprintln!("   node {}: material {} != last usemtl {want} of its group", node.name, m.name);
            }
            (None, Some(want)) => {
                c.failures += 1;
                eprintln!("   node {}: no material attached, group ends with usemtl {want}", node.name);
            }
            (Some(m), None) => {
                c.failures += 1;
                eprintln!("   node {}: material {} attached, but its group never sets usemtl", node.name, m.name);
            }
            (None, None) => {}
        }
    }

    // Every group that has faces must come back as a node with those faces, and
    // its face count must survive the fan triangulation (>= the OBJ polygon
    // count). This catches a reader that drops a group between `g` records.
    for (name, nf) in &w.group_faces {
        match scene.nodes.iter().find(|n| n.name == *name) {
            Some(node) if node.triangles.len() >= *nf && !node.vertices.is_empty() => {}
            Some(node) => {
                c.failures += 1;
                eprintln!("   group {name}: {nf} faces written, {} triangles read", node.triangles.len());
            }
            None => {
                c.failures += 1;
                eprintln!("   group {name}: {nf} faces written, node missing after read");
            }
        }
    }
    Ok(())
}

fn print_node(name: &str, nv: usize, nf: usize, mat: Option<&Material>, uv: bool, tex: Option<&str>) {
    let tex = tex.unwrap_or("-");
    match mat {
        Some(m) => println!(
            "   node {name:<24} v={nv:<7} f={nf:<7} uv={:<3} mat={} ka=({:.6},{:.6},{:.6}) kd=({:.6},{:.6},{:.6}) ks=({:.6},{:.6},{:.6}) ns={:.3} opacity={:.3} tex={tex}",
            if uv { "on" } else { "off" },
            m.name,
            m.ambient.0,
            m.ambient.1,
            m.ambient.2,
            m.diffuse.0,
            m.diffuse.1,
            m.diffuse.2,
            m.specular.0,
            m.specular.1,
            m.specular.2,
            m.shininess,
            m.opacity
        ),
        None => println!("   node {name:<24} v={nv:<7} f={nf:<7} uv={:<3} mat=- tex={tex}", if uv { "on" } else { "off" }),
    }
}

/// Compare one attached material against its written `newmtl` block.
///
/// `RWObj_ObjMaterialMap::DefineMaterial` writes `Ka`/`Kd`/`Ks`/`Ns` together
/// (`RWObj_ObjMaterialMap.cxx:85-121`) or none of them (a style with no surface
/// color), so a field the file does not carry must read back as the
/// `RWObj_Material` default (`RWObj_Material.hxx:45-46`), and `Tr` missing means
/// opaque. `Ns` is normalised by the reader (`Shininess = min(Ns/1000, 1)`,
/// `RWObj_MtlReader.cxx:200-211`), so 1000.000000 comes back as 1.0.
fn check_material(m: &Material, w: &Written, node: &occt_topo::rwmesh::MeshNode, c: &mut Checks) {
    let Some((_, block)) = w.mats.iter().find(|(n, _)| *n == m.name) else {
        c.failures += 1;
        eprintln!("   node {}: material {} has no newmtl in the written library", node.name, m.name);
        return;
    };
    let mut bad = Vec::new();
    match block.ka {
        Some(ka) if !same_rgb(m.ambient, ka) => bad.push(format!("Ka {:?} != {:?}", m.ambient, ka)),
        None if !same_rgb(m.ambient, (0.1, 0.1, 0.1)) => bad.push(format!("Ka {:?} != default", m.ambient)),
        _ => {}
    }
    match block.kd {
        Some(kd) if !same_rgb(m.diffuse, kd) => bad.push(format!("Kd {:?} != {:?}", m.diffuse, kd)),
        None if !same_rgb(m.diffuse, (0.8, 0.8, 0.8)) => bad.push(format!("Kd {:?} != default", m.diffuse)),
        _ => {}
    }
    match block.ks {
        Some(ks) if !same_rgb(m.specular, ks) => bad.push(format!("Ks {:?} != {:?}", m.specular, ks)),
        None if !same_rgb(m.specular, (0.2, 0.2, 0.2)) => bad.push(format!("Ks {:?} != default", m.specular)),
        _ => {}
    }
    match block.ns {
        Some(ns) => {
            let want = (ns / 1000.0).min(1.0);
            if (m.shininess - want).abs() > 1e-9 {
                bad.push(format!("Ns {ns} -> shininess {:?} != {want}", m.shininess));
            }
        }
        None if (m.shininess - 1.0).abs() > 1e-9 => bad.push(format!("shininess {:?} != default", m.shininess)),
        _ => {}
    }
    match block.tr {
        Some(tr) if (m.opacity - (1.0 - tr)).abs() > 1e-9 => {
            bad.push(format!("Tr {tr} -> opacity {:?} != {}", m.opacity, 1.0 - tr));
        }
        None if (m.opacity - 1.0).abs() > 1e-9 => bad.push(format!("opacity {:?} != opaque", m.opacity)),
        _ => {}
    }
    // `RWObj_MtlReader.cxx:242-267` reads all three maps through
    // `processTexturePath`; a missing line leaves the path empty (`None`).
    for (tag, written, got) in [
        ("map_Kd", &block.map_kd, &m.diffuse_texture),
        ("map_Ks", &block.map_ks, &m.specular_texture),
        ("map_Bump", &block.map_bump, &m.bump_texture),
    ] {
        match (written, got) {
            (Some(want), Some(have)) if want != have => bad.push(format!("{tag} {have} != {want}")),
            (Some(want), None) => bad.push(format!("{tag} {want} not read")),
            (None, Some(have)) => bad.push(format!("{tag} {have} read but not written")),
            _ => {}
        }
    }
    // The node's `texture_path` is the material's diffuse map, the same reference
    // `RWObj_CafReader` puts into `XCAFDoc_VisMaterialCommon::DiffuseTexture`.
    if node.texture_path != m.diffuse_texture {
        bad.push(format!("node texture_path {:?} != diffuse map {:?}", node.texture_path, m.diffuse_texture));
    }
    if bad.is_empty() {
        return;
    }
    if bad.iter().any(|b| b.starts_with("Kd")) {
        c.kd_mismatch += 1;
    }
    c.aspect_mismatch += 1;
    c.failures += 1;
    eprintln!("   node {}: material {}: {}", node.name, m.name, bad.join(", "));
}

fn names(scene: &MeshScene) -> Vec<&str> {
    scene.materials.iter().map(|m| m.name.as_str()).collect()
}

fn same_rgb(a: (f64, f64, f64), b: (f64, f64, f64)) -> bool {
    (a.0 - b.0).abs() < 1e-6 && (a.1 - b.1).abs() < 1e-6 && (a.2 - b.2).abs() < 1e-6
}
