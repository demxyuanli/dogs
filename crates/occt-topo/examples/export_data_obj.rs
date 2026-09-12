//! Overall-test dump: every `data/*.step` (and `.stp`) to `output/<stem>.obj`.
//!
//! Pipeline: `read_step_file` → compound if multiple roots → `brep_to_obj(..., 0.1)`.
//! Re-run before visual / mesh overall checks. Do not add a second exporter.
use std::fs;
use std::path::PathBuf;

use occt_topo::brep_exchange::brep_to_obj;
use occt_topo::builder::TopoBuilder;
use occt_topo::step::read_step_file;

fn main() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..").join("..");
    let data = root.join("data");
    let out = root.join("output");
    fs::create_dir_all(&out).expect("create output/");

    let mut steps: Vec<PathBuf> = fs::read_dir(&data)
        .expect("read data/")
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| {
            p.extension()
                .and_then(|e| e.to_str())
                .map(|e| e.eq_ignore_ascii_case("step") || e.eq_ignore_ascii_case("stp"))
                .unwrap_or(false)
        })
        .collect();
    steps.sort();

    let mut failed = 0usize;
    for step in &steps {
        let stem = step.file_stem().unwrap().to_string_lossy();
        let dest = out.join(format!("{stem}.obj"));
        match export_one(step) {
            Ok(obj) => {
                let nv = obj.lines().filter(|l| l.starts_with("v ")).count();
                let nf = obj.lines().filter(|l| l.starts_with("f ")).count();
                fs::write(&dest, obj).unwrap_or_else(|e| panic!("write {}: {e}", dest.display()));
                println!("ok  {stem}.step -> output/{stem}.obj  v={nv} f={nf}");
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

fn export_one(step: &PathBuf) -> Result<String, String> {
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
    Ok(brep_to_obj(&shape, 0.1))
}
