//! TEMPORARY gate probe: dump `BRepAlgo_FaceRestrictor` results for one face,
//! line-for-line comparable with
//! `.target-gate/make_offset_gt/domains_probe.cpp`.
//!
//! Usage: cargo run --example face_restrict_probe -- data/HoledPlate.step 0

use std::path::PathBuf;

use occt_topo::abs::ShapeType;
use occt_topo::brep_algo_face_restrictor::BRepAlgoFaceRestrictor;
use occt_topo::brep_tool::BRepTool;
use occt_topo::brep_tools_wire_explorer::WireExplorer;
use occt_topo::shape::{Face, TopoShape};
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


/// Experiment: does flagging the face wires `Closed` change
/// `BRepOffsetAPI_MakeOffset` (OCCT reports `closed=1` for the STEP wires)?
fn offset_with_closed_flag(face: &Face, flag: bool, join: occt_core::kernel::geomabs::JoinType, off: f64) {
    use occt_topo::brep_offset_api_make_offset::BRepOffsetApiMakeOffset;
    if flag {
        for w in wires_of_face(face) {
            w.set_closed(true);
        }
    } else {
        let _ = face;
    }
    let mut mo = BRepOffsetApiMakeOffset::with_face(face, join, false);
    mo.perform(off, 0.0);
    let nw = mo.shape().map(|s| wires_of(s).len()).unwrap_or(0);
    println!(
        "OFF closed_flag={} done={} nw={}",
        if flag { 1 } else { 0 },
        if mo.is_done() { 1 } else { 0 },
        nw
    );
}
fn main() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..").join("..");
    let arg = std::env::args().nth(1).unwrap_or_else(|| "data/HoledPlate.step".into());
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

    let mut fr = BRepAlgoFaceRestrictor::new();
    fr.init(Some(face.clone()), true, false);

    for (wi, w) in wires_of_face(&face).into_iter().enumerate() {
        let ne = explorer_shapes(&w, ShapeType::Edge).len();
        println!(
            "IN {} ne={} closed={} ori={}",
            wi,
            ne,
            if w.is_closed() { 1 } else { 0 },
            w.orientation() as i32
        );
        fr.add(&w);
    }

    fr.perform();
    let mut out: Vec<TopoShape> = Vec::new();
    while fr.more() {
        if let Some(f) = fr.current() {
            out.push(f.0);
        }
        fr.next();
    }
    println!("DONE {} nfaces={}", if fr.is_done() { 1 } else { 0 }, out.len());

    for (k, f) in out.iter().enumerate() {
        let fw = explorer_shapes(f, ShapeType::Wire);
        println!("FACE {} nwires={}", k, fw.len());
        for (wj, w) in fw.iter().enumerate() {
            let mut e = WireExplorer::of_wire(&occt_topo::shape::Wire((*w).clone()));
            let mut ne = 0usize;
            let mut pts = Vec::new();
            while e.more() {
                ne += 1;
                if ne > 512 {
                    break;
                }
                if let Some(v) = e.current_vertex() {
                    pts.push(BRepTool::vertex_point(v));
                }
                e.next();
            }
            println!(
                "W {} {} ne={} closed={} ori={}",
                k,
                wj,
                ne,
                if w.is_closed() { 1 } else { 0 },
                w.orientation() as i32
            );
            for (p, pt) in pts.iter().enumerate() {
                println!(
                    "V {} {} {} {:.6} {:.6} {:.6}",
                    k,
                    wj,
                    p,
                    pt.x(),
                    pt.y(),
                    pt.z()
                );
            }
        }
    }
    {
        let w0 = wires_of_face(&face)[0].clone();
        let wr = w0.clone().oriented(occt_topo::abs::Orientation::Reversed);
        println!("ORI w0={} wr={}", w0.orientation() as i32, wr.orientation() as i32);
        for (tag, w) in [
            ("fwd", w0.clone()),
            ("rev", occt_topo::shape::Wire(wr.clone())),
        ] {
            let cum = occt_topo::iterator::cumulated_children(&w.0);
            println!(
                "CUM {} ori0={} ori1={} n={}",
                tag,
                cum.first().map(|c| c.orientation() as i32).unwrap_or(-1),
                cum.get(1).map(|c| c.orientation() as i32).unwrap_or(-1),
                cum.len()
            );
            let mut e = WireExplorer::of_wire(&w);
            let mut k = 0;
            while e.more() && k < 2 {
                if let Some(v) = e.current_vertex() {
                    let p = BRepTool::vertex_point(v);
                    println!("ORI {} {} {:.3} {:.3}", tag, k, p.x(), p.y());
                }
                k += 1;
                e.next();
            }
        }
    }
    println!("PROBE END");
    let f2 = Face(faces[fi].clone());
    offset_with_closed_flag(&f2, false, occt_core::kernel::geomabs::JoinType::Arc, 1.0);
    let f3 = Face(faces[fi].clone());
    offset_with_closed_flag(&f3, true, occt_core::kernel::geomabs::JoinType::Arc, 1.0);
    // Single-wire faces: is the outer offset broken on its own?
    for (wi, w) in wires_of_face(&f2).into_iter().enumerate() {
        let Some(single) = occt_topo::brep_lib_make_face::make_face_from_wire(
            &w.clone(),
            true,
        ) else {
            println!("SINGLE {} face=none", wi);
            continue;
        };
        let nw_in = wires_of_face(&single).len();
        let mut mo = occt_topo::brep_offset_api_make_offset::BRepOffsetApiMakeOffset::with_face(
            &single,
            occt_core::kernel::geomabs::JoinType::Arc,
            false,
        );
        mo.perform(1.0, 0.0);
        let nw = mo.shape().map(|s| wires_of(s).len()).unwrap_or(0);
        println!(
            "SINGLE {} inwires={} done={} nw={}",
            wi,
            nw_in,
            if mo.is_done() { 1 } else { 0 },
            nw
        );
    }
}
