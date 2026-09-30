//! TEMPORARY probe (T-94): does `fix_missing_seam` change the ORIGINAL face's
//! pcurves? The `swap_seam` step writes pcurves for the temporary face it builds
//! from `surf.clone()`, and `repr_key` resolves by surface pointer, so the write
//! can land on the original face's slot. This measures it.
//!
//! Usage: `cargo run --example zz_share_probe -- <file.step> [max-faces]`
use occt_topo::shape::Face;
use occt_topo::shhealing::ShapeFixFace;
use occt_topo::step::read_step_file;
use occt_topo::tgeometry::GeometryRegistry;
use occt_topo::topo_tools_full::{edges_of_wire, faces_of, wires_of_face};

/// One slot in the face key's table: how many pcurves, the stored range, and a
/// fingerprint of each pcurve (domain + evaluated end points), so a reorder is
/// visible rather than just "changed".
fn slot_desc(pcs: &[std::sync::Arc<dyn occt_geom2d::Curve2d>], rng: Option<(f64, f64)>) -> String {
    let mut s = format!("n={} rng={rng:?}", pcs.len());
    for (k, pc) in pcs.iter().enumerate() {
        let (a, b) = (pc.first_parameter(), pc.last_parameter());
        let (sa, sb) = if a.is_finite() && b.is_finite() {
            (format!("{a:.9}"), format!("{b:.9}"))
        } else {
            ("-inf".into(), "+inf".into())
        };
        // Evaluate inside the curve's own domain; fall back to 0/1 for an
        // unbounded curve (e.g. a Geom2d_Line, whose domain is (-inf, inf)).
        let (t0, t1) = if a.is_finite() && b.is_finite() { (a, b) } else { (0.0, 1.0) };
        let p0 = pc.d0(t0);
        let p1 = pc.d0(t1);
        let mid = pc.d0(0.5 * (t0 + t1));
        s += &format!(
            " | [{k}] dom=[{sa},{sb}] line={} p0=({:.9},{:.9}) p1=({:.9},{:.9}) mid=({:.9},{:.9})",
            pc.is_line(),
            p0.x(), p0.y(), p1.x(), p1.y(), mid.x(), mid.y()
        );
    }
    s
}

fn snapshot(face: &Face) -> Vec<(usize, usize, String)> {
    let reg = GeometryRegistry::global();
    let fk = GeometryRegistry::shape_key(&face.0);
    let mut out = Vec::new();
    for (wi, w) in wires_of_face(face).iter().enumerate() {
        for (ei, e) in edges_of_wire(w).iter().enumerate() {
            let pcs = reg.edge_pcurves(&e.0, fk);
            let rng = reg.pcurve_range(&e.0, fk);
            out.push((wi, ei, slot_desc(&pcs, rng)));
        }
    }
    out
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = args.get(1).expect("usage: zz_share_probe <file.step> [max-faces]");
    let max: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(usize::MAX);
    let model = read_step_file(path).expect("read step");
    let shape = model.shapes[0].shape.clone();

    let all = faces_of(&shape);
    let mut changed_face = 0usize;
    let mut changed_edges = 0usize;
    let mut checked = 0usize;
    for (i, f) in all.iter().take(max).enumerate() {
        let before = snapshot(f);
        if before.is_empty() {
            continue;
        }
        let mut sff = ShapeFixFace::with_face(f);
        let fixed = sff.fix_missing_seam();
        if !fixed {
            continue;
        }
        checked += 1;
        let after = snapshot(f);
        let mut new_edges = 0usize;
        let mut diffs = Vec::new();
        for (k, (wa, ea, va)) in before.iter().enumerate() {
            match after.get(k) {
                Some((_, _, vb)) if vb == va => {}
                Some((_, _, vb)) => {
                    diffs.push((k, *wa, *ea, va.clone(), vb.clone()));
                }
                None => {}
            }
        }
        // slots that only exist after (the seam step added edges)
        if after.len() > before.len() {
            new_edges = after.len() - before.len();
        }
        if !diffs.is_empty() {
            changed_face += 1;
            changed_edges += diffs.len();
            println!(
                "SHARE face={i} fix=true changed_slots={} (of {}) new_edges={new_edges}",
                diffs.len(),
                before.len()
            );
            for (k, wi, ei, va, vb) in &diffs {
                println!("SHARE   slot{k} wire{wi}.edge{ei}");
                println!("SHARE     before: {va}");
                println!("SHARE     after : {vb}");
            }
        }
    }
    println!(
        "SHARE TOTAL faces_with_fix={checked} faces_changed={changed_face} slots_changed={changed_edges}"
    );
}
