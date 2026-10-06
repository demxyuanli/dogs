//! Overall-test dump: every `data/*.step` (and `.stp`) to `data/output/<stem>.obj`.
//!
//! Pipeline: `read_step_file` → compound if multiple roots → `brep_to_obj(..., 0.1)`.
//! The second argument is the drawer's `MaximalChordialDeviation` (the fallback
//! used only for a void/unbounded bounding box); the linear deflection actually
//! driving the mesh is `Prs3d::GetDeflection(shape, drawer)` =
//! `maxComp(bbox) * 0.001 * 4` (`Prs3d.hxx:82-103`), so changing `0.1` does not
//! change the export density. Re-run before visual / mesh overall checks. Do not
//! add a second exporter.
//!
//! Writes into the single `data/output/` dump directory shared with the
//! `step_obj_gates` test, so there is one artifact location to inspect.
//!
//! Optional positional arguments name extra source directories (relative to the
//! repository root) to scan instead of `data/`, so the OCCT test models under
//! `data/occ/` can be dumped with the same pipeline:
//! `cargo run --example export_data_obj -- data/occ`.
use std::fs;
use std::path::PathBuf;

use occt_topo::brep_exchange::brep_to_obj;
use occt_topo::builder::TopoBuilder;
use occt_topo::step::read_step_file;

fn main() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..").join("..");
    let out = root.join("data").join("output");
    fs::create_dir_all(&out).expect("create data/output/");

    let args: Vec<String> = std::env::args().skip(1).collect();
    let dirs: Vec<PathBuf> = if args.is_empty() {
        vec![root.join("data")]
    } else {
        args.iter().map(|a| root.join(a)).collect()
    };

    let mut steps: Vec<PathBuf> = Vec::new();
    for dir in &dirs {
        steps.extend(
            fs::read_dir(dir)
                .unwrap_or_else(|e| panic!("read {}: {e}", dir.display()))
                .filter_map(|e| e.ok())
                .map(|e| e.path())
                .filter(|p| {
                    p.extension()
                        .and_then(|e| e.to_str())
                        .map(|e| e.eq_ignore_ascii_case("step") || e.eq_ignore_ascii_case("stp"))
                        .unwrap_or(false)
                }),
        );
    }
    steps.sort();

    let mut failed = 0usize;
    for step in &steps {
        let stem = step.file_stem().unwrap().to_string_lossy();
        let dest = out.join(format!("{stem}.obj"));
        match export_one(step) {
            Ok((obj, exact)) => {
                let nv = obj.lines().filter(|l| l.starts_with("v ")).count();
                let nf = obj.lines().filter(|l| l.starts_with("f ")).count();
                fs::write(&dest, obj).unwrap_or_else(|e| panic!("write {}: {e}", dest.display()));
                println!("ok  {stem}.step -> data/output/{stem}.obj  v={nv} f={nf}{exact}");
            }
            Err(e) => {
                failed += 1;
                eprintln!("err {stem}.step: {e}");
            }
        }
    }
    if failed > 0 {
        std::process::exit(1);
    }
}

fn export_one(step: &PathBuf) -> Result<(String, String), String> {
    let model = read_step_file(&step.to_string_lossy())?;
    if model.shapes.is_empty() {
        return Err("parsed to no shapes".into());
    }
    let shape = if model.shapes.len() == 1 {
        model.shapes[0].shape.clone()
    } else {
        let parts: Vec<_> = model.shapes.iter().map(|s| s.shape.clone()).collect();
        TopoBuilder::new().make_compound_of(&parts).0
    };
    // Exact box (`BRepBndLib::Add` port: analytic extrema / control net, not the
    // tessellation), printed next to the mesh counts so a bbox drift can be
    // attributed to the parsed geometry or to the mesher at a glance.
    let exact = match occt_topo::brep_bnd_lib::shape_bnd_box(&shape).get() {
        Some((x0, x1, y0, y1, z0, z1)) => format!(
            "  exact=[{x0:.3},{y0:.3},{z0:.3}]~[{x1:.3},{y1:.3},{z1:.3}]"
        ),
        None => String::new(),
    };
    Ok((brep_to_obj(&shape, 0.1), exact))
}
