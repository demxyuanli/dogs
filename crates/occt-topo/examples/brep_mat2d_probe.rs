//! TEMPORARY probe: run `BRepMAT2d_Explorer` + `BRepMAT2d_BisectingLocus` over
//! every face of a STEP file and dump the resulting MAT graph sizes.
//!
//! Usage: cargo run --example brep_mat2d_probe -- data/HoledPlate.step
//!
//! Delete this file once `BRepFill_OffsetAncestors` supplies a real baseline.

use std::path::PathBuf;

use occt_core::kernel::geomabs::JoinType;
use occt_topo::brep_mat2d_bisecting_locus::BRepMat2dBisectingLocus;
use occt_topo::brep_mat2d_explorer::BRepMat2dExplorer;
use occt_topo::brep_tools_wire_explorer::WireExplorer;
use occt_topo::mat::MatSide;
use occt_topo::shape::{Face, TopoShape};
use occt_topo::step::read_step_file;
use occt_topo::topo_tools::wires_of_face;

fn collect_faces(shape: &TopoShape, out: &mut Vec<TopoShape>) {
    use occt_topo::abs::ShapeType;
    let kids: Vec<TopoShape> = shape.tshape.read().unwrap().children.clone();
    for k in kids {
        match k.shape_type() {
            ShapeType::Face => out.push(k),
            ShapeType::Shell
            | ShapeType::Solid
            | ShapeType::CompSolid
            | ShapeType::Compound => collect_faces(&k, out),
            _ => {}
        }
    }
}

fn main() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..").join("..");
    let arg = std::env::args().nth(1).unwrap_or_else(|| "data/HoledPlate.step".into());
    let path = if PathBuf::from(&arg).is_absolute() {
        PathBuf::from(&arg)
    } else {
        root.join(&arg)
    };
    let stack_kb: usize = std::env::var("PROBE_STACK_KB")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(2048);
    let handle = std::thread::Builder::new()
        .stack_size(stack_kb * 1024)
        .spawn(move || run(path))
        .expect("probe thread");
    handle.join().expect("probe thread join");
}

fn run(path: PathBuf) {
    println!("PROBE {}", path.display());

    let shape = read_step_file(&path.to_string_lossy()).expect("read_step_file");
    let Some(shape) = occt_topo::brep_exchange::document_shape(&shape) else {
        println!("no document shape");
        return;
    };
    let mut faces = Vec::new();
    collect_faces(&shape, &mut faces);
    println!("faces={}", faces.len());

    for (fi, f) in faces.iter().enumerate() {
        let face = Face(f.clone());
        for (wi, w) in wires_of_face(&face).iter().enumerate() {
            let mut e = WireExplorer::of_wire_on_face(w, &face);
            let mut steps = 0usize;
            let mut names = Vec::new();
            while e.more() {
                steps += 1;
                names.push(format!("{:?}", e.current().map(|x| x.0.orientation())));
                if steps > 64 {
                    println!("[{fi}/{wi}] EXPLORER CYCLES: {names:?}");
                    break;
                }
                e.next();
            }
            print!("[{fi}/{wi}] steps={steps} closed={} ", w.0.is_closed());
        }
        println!("| surface={}", occt_topo::tgeometry::GeometryRegistry::global().face_surface(&face.0).is_some());

        let mut explo = BRepMat2dExplorer::new();
        explo.perform(&face);
        let n_contours = explo.number_of_contours();
        let mut n_curves = 0usize;
        for c in 1..=n_contours {
            n_curves += explo.number_of_curves(c);
        }
        println!("[{fi}] perform done contours={n_contours} curves={n_curves}");

        for side in [MatSide::Left, MatSide::Right] {
            let mut locus = BRepMat2dBisectingLocus::new();
            let prev = std::panic::take_hook();
            std::panic::set_hook(Box::new(|_| {}));
            let res = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                locus.compute(&mut explo, 1, side, JoinType::Arc, false);
            }));
            std::panic::set_hook(prev);
            if res.is_err() {
                println!("[{fi}] {side:?} FAILED");
                continue;
            }
            let g = locus.graph();
            {
                let circ = locus.the_tool().circuit();
                eprintln!("RTCIRC items={}", circ.number_of_items());
                for k in 1..=circ.number_of_items() {
                    if circ.connexion_on(k) {
                        let c = circ.connexion(k).borrow();
                        eprintln!(
                            "RTCIRC conn key={k} first={},{} second={},{} l1={} i1={} l2={} i2={}",
                            c.point_on_first().x(),
                            c.point_on_first().y(),
                            c.point_on_second().x(),
                            c.point_on_second().y(),
                            c.index_first_line(),
                            c.index_item_on_first(),
                            c.index_second_line(),
                            c.index_item_on_second()
                        );
                    }
                }
            }
            println!(
                "[{fi}] {side:?} semi={} nbis={}",
                locus.the_mat().semi_infinite(),
                locus.the_mat().number_of_bisectors()
            );
            println!(
                "[{fi}] {side:?} done={} elts={} bis={} arcs={} nodes={}",
                locus.is_done(),
                g.number_of_basic_elts(),
                g.number_of_arcs() + g.number_of_infinite_nodes(),
                g.number_of_arcs(),
                g.number_of_nodes()
            );
        }
    }
    println!("PROBE END");
}
