//! TEMPORARY gate probe: does `BRepFill_OffsetWire` care about the spine wire
//! orientation? Builds the same rectangle face twice (wire FORWARD / REVERSED)
//! and reports the MAT2d locus and the offset result for both.

use std::path::PathBuf;

use occt_core::kernel::geomabs::JoinType;
use occt_topo::abs::{Orientation, ShapeType};
use occt_topo::brep_fill_offset_wire::BRepFillOffsetWire;
use occt_topo::brep_lib_find_surface::FindSurface;
use occt_topo::brep_mat2d_bisecting_locus::BRepMat2dBisectingLocus;
use occt_topo::brep_mat2d_explorer::BRepMat2dExplorer;
use occt_topo::builder::TopoBuilder;
use occt_topo::mat::MatSide;
use occt_topo::shape::{Face, TopoShape, Wire};
use occt_topo::step::read_step_file;
use occt_topo::topo_tools::{explorer_shapes, wires_of, wires_of_face};

fn collect_faces(shape: &TopoShape, out: &mut Vec<TopoShape>) {
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

fn face_from(wire: &Wire, like: &Face) -> Option<Face> {
    let surface = occt_topo::brep_tool::BRepTool::face_surface(like)?;
    let tol = occt_topo::brep_tool::BRepTool::face_tolerance(like);
    Some(TopoBuilder::new().make_face_with_tolerance(surface, &[wire.clone()], tol))
}

fn probe_open(tag: &str, face: &Face, off: f64) {
    let mut ow = BRepFillOffsetWire::with_face(face, JoinType::Arc, true);
    ow.perform(off, 0.0);
    let nw = ow.shape().map(|s| wires_of(s).len()).unwrap_or(0);
    println!(
        "OPEN {} off={:.1} done={} nwires={}",
        tag,
        off,
        if ow.is_done() { 1 } else { 0 },
        nw
    );
}

fn probe(tag: &str, face: &Face, off: f64) {
    let mut exp = BRepMat2dExplorer::new();
    exp.perform(face);
    let mut bilo = BRepMat2dBisectingLocus::new();
    bilo.compute(&mut exp, 1, MatSide::Left, JoinType::Arc, false);
    println!(
        "{} locus_done={} contours={}",
        tag,
        if bilo.is_done() { 1 } else { 0 },
        bilo.number_of_contours()
    );

    let mut ow = BRepFillOffsetWire::with_face(face, JoinType::Arc, false);
    ow.perform(off, 0.0);
    let nw = ow.shape().map(|s| wires_of(s).len()).unwrap_or(0);
    println!(
        "{} off={:.1} done={} nwires={}",
        tag,
        off,
        if ow.is_done() { 1 } else { 0 },
        nw
    );
}

fn main() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..").join("..");
    let arg = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "data/HoledPlate.step".into());
    let fi: usize = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(0);
    let path = if PathBuf::from(&arg).is_absolute() {
        PathBuf::from(&arg)
    } else {
        root.join(&arg)
    };
    let shape = read_step_file(&path.to_string_lossy()).expect("read_step_file");
    let Some(shape) = occt_topo::brep_exchange::document_shape(&shape) else {
        println!("no document shape");
        return;
    };
    let mut faces = Vec::new();
    collect_faces(&shape, &mut faces);
    let face = Face(faces[fi].clone());
    let outer = wires_of_face(&face)[0].clone();
    let outer_rev = Wire(outer.0.clone().oriented(Orientation::Reversed));

    let Some(fwd) = face_from(&outer, &face) else {
        println!("fwd face: none");
        return;
    };
    let Some(rev) = face_from(&outer_rev, &face) else {
        println!("rev face: none");
        return;
    };
    println!(
        "faces fwd_wires={} rev_wires={} fwd_ori={} rev_ori={}",
        wires_of_face(&fwd).len(),
        wires_of_face(&rev).len(),
        wires_of_face(&fwd)[0].0.orientation() as i32,
        wires_of_face(&rev)[0].0.orientation() as i32
    );

    probe_open("REV", &rev, 1.0);
    probe_open("FWD", &fwd, 1.0);
    probe("REV", &rev, 0.0);
    probe("FWD", &fwd, 0.0);
    probe("FWD", &fwd, 1.0);
    probe("REV", &rev, 1.0);
    probe("FWD", &fwd, -1.0);
    probe("REV", &rev, -1.0);
    println!("edges={}", explorer_shapes(&outer.0, ShapeType::Edge).len());
}
