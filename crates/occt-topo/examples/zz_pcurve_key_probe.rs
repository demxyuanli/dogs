//! TEMPORARY probe (T-99): invariant check.
//!
//! For every face `F` of the read shape, verify that each of `F`'s edges has a
//! pcurve resolvable under `repr_key(F)` - i.e. call
//! `GeometryRegistry::edge_pcurves(edge, shape_key(F))`, which goes through
//! `repr_key` (`tgeometry.rs:524-531`) and therefore answers the question the
//! mesher asks: "does this face have a pcurve for this edge *on my surface*?"
//!
//! Reports the (face, edge) pairs that have none. Run against the baseline and
//! against the `unregister_face_surface(tmp_f)` experiment shape; if the
//! experiment violates the invariant on faces the baseline does not, the cause
//! of the area under-cover is "pcurve stored under a mismatched key".
//!
//! Usage: `cargo run --example zz_pcurve_key_probe -- <file.step>`
use occt_topo::shape::Face;
use occt_topo::step::read_step_file;
use occt_topo::tgeometry::GeometryRegistry;
use occt_topo::topo_tools_full::{edges_of_wire, faces_of, wires_of_face};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = args.get(1).expect("usage: zz_pcurve_key_probe <file.step>");
    let model = read_step_file(path).expect("read step");
    let shape = model.shapes[0].shape.clone();
    let reg = GeometryRegistry::global();

    let faces = faces_of(&shape);
    let mut total_edges = 0usize;
    let mut no_pcurve = 0usize;
    let mut faces_with_missing = 0usize;
    let mut worst: Vec<(usize, usize, usize, usize)> = Vec::new(); // face, edges, missing, wires

    for (fi, f) in faces.iter().enumerate() {
        let face = Face(f.0.clone());
        let fk = GeometryRegistry::shape_key(&face.0);
        let mut n_edges = 0usize;
        let mut n_missing = 0usize;
        for w in wires_of_face(&face) {
            for e in edges_of_wire(&w) {
                n_edges += 1;
                if reg.edge_pcurves(&e.0, fk).is_empty() {
                    n_missing += 1;
                }
            }
        }
        total_edges += n_edges;
        no_pcurve += n_missing;
        if n_missing > 0 {
            faces_with_missing += 1;
            worst.push((fi, n_edges, n_missing, wires_of_face(&face).len()));
        }
    }

    worst.sort_by_key(|t| std::cmp::Reverse(t.2));
    println!(
        "KEYPROBE faces={} edges={total_edges} edges_without_pcurve={no_pcurve} faces_with_missing={faces_with_missing}",
        faces.len()
    );
    for (fi, ne, nm, nw) in worst.iter().take(25) {
        println!("KEYPROBE   face={fi} edges={ne} missing={nm} wires={nw}");
    }

    // Per-face edge sets, so two runs (baseline vs experiment) can be diffed:
    // an edge is identified by its shape key, which is stable for a given input.
    if std::env::var_os("OCCT_TOPO_KEY_DUMP").is_some() {
        for (fi, f) in faces.iter().enumerate() {
            let face = Face(f.0.clone());
            let fk = GeometryRegistry::shape_key(&face.0);
            let mut keys: Vec<usize> = Vec::new();
            for w in wires_of_face(&face) {
                for e in edges_of_wire(&w) {
                    keys.push(GeometryRegistry::shape_key(&e.0));
                }
            }
            keys.sort_unstable();
            println!("KEYEDGE face={fi} fk={fk} n={} keys={keys:?}", keys.len());
        }
    }
}
