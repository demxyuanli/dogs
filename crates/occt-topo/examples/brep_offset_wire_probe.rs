//! TEMPORARY gate probe: run `BRepFill_OffsetWire` over every planar face of a
//! STEP file and dump the resulting wires.
//!
//! Usage: cargo run --example brep_offset_wire_probe -- data/HoledPlate.step
//!
//! The output lines are line-for-line comparable with the OCCT side probe
//! (`.target-gate/offset_wire_gt/probe.cpp`).

use std::path::PathBuf;

use occt_core::gp::GpPnt;
use occt_core::kernel::geomabs::JoinType;
use occt_topo::abs::ShapeType;
use occt_topo::brep_fill_offset_wire::BRepFillOffsetWire;
use occt_topo::brep_surface::{face_is_plane_type, face_plane};
use occt_topo::brep_tool::BRepTool;
use occt_topo::brep_tools_wire_explorer::WireExplorer;
use occt_topo::shape::{Face, TopoShape};
use occt_exchange::step::read_step_file;
use occt_topo::topo_tools::{wires_of, wires_of_face};

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
        .unwrap_or(4096);
    let handle = std::thread::Builder::new()
        .stack_size(stack_kb * 1024)
        .spawn(move || run(path))
        .expect("probe thread");
    handle.join().expect("probe thread join");
}

fn join_name(j: JoinType) -> &'static str {
    match j {
        JoinType::Arc => "Arc",
        JoinType::Tangent => "Tangent",
        JoinType::Intersection => "Intersection",
    }
}

fn dump_vertex(prefix: &str, p: &GpPnt) {
    println!(
        "{} {:.6} {:.6} {:.6}",
        prefix,
        p.x(),
        p.y(),
        p.z()
    );
}

fn run(path: PathBuf) {
    println!("PROBE {}", path.display());

    let shape = read_step_file(&path.to_string_lossy()).expect("read_step_file");
    let Some(shape) = occt_exchange::brep_exchange::document_shape(&shape) else {
        println!("no document shape");
        return;
    };
    let mut faces = Vec::new();
    collect_faces(&shape, &mut faces);
    println!("faces={}", faces.len());
    let only_face: Option<usize> = std::env::var("PROBE_FACE").ok().and_then(|s| s.parse().ok());

    for (fi, f) in faces.iter().enumerate() {
        if let Some(only) = only_face {
            if fi != only {
                continue;
            }
        }
        let face = Face(f.clone());
        let is_plane = face_is_plane_type(&face);
        let nw = wires_of_face(&face).len();
        println!("FACE {} plane={} wires={}", fi, if is_plane { 1 } else { 0 }, nw);
        if !is_plane || nw == 0 {
            continue;
        }
        if let Some(pln) = face_plane(&face) {
            let n = pln.axis().direction().xyz().clone();
            println!("PLANE {} nx={:.6} ny={:.6} nz={:.6}", fi, n.x, n.y, n.z);
        } else {
            println!("PLANE {} none", fi);
        }

        for join in [JoinType::Arc, JoinType::Tangent, JoinType::Intersection] {
            if let Ok(want) = std::env::var("PROBE_JOIN") {
                if join_name(join) != want {
                    continue;
                }
            }
            for off in [1.0f64, -1.0, -5.0] {
                let prev = std::panic::take_hook();
                let direct = std::env::var("PROBE_NO_CATCH").is_ok();
                if !direct {
                    std::panic::set_hook(Box::new(|_| {}));
                }
                let body = || {
                    let mut ow = BRepFillOffsetWire::with_face(&face, join, false);
                    ow.perform(off, 0.0);
                    (ow.is_done(), ow.shape().cloned())
                };
                let res = if direct {
                    Ok(body())
                } else {
                    std::panic::catch_unwind(std::panic::AssertUnwindSafe(body))
                };
                std::panic::set_hook(prev);
                match res {
                    Err(e) => {
                        let msg = if let Some(s) = e.downcast_ref::<&str>() {
                            (*s).to_string()
                        } else if let Some(s) = e.downcast_ref::<String>() {
                            s.clone()
                        } else {
                            "?".to_string()
                        };
                        println!("OW {} {} {:.6} PANIC {}", fi, join_name(join), off, msg);
                        continue;
                    }
                    Ok((false, _)) => {
                        println!("OW {} {} {:.6} done=0", fi, join_name(join), off);
                        continue;
                    }
                    Ok((true, shape)) => {
                        let Some(shape) = shape else {
                            println!("OW {} {} {:.6} done=1 nw=0", fi, join_name(join), off);
                            continue;
                        };
                        let wires = wires_of(&shape);
                        println!("OW {} {} {:.6} done=1 nw={}", fi, join_name(join), off, wires.len());
                        for (wi, w) in wires.iter().enumerate() {
                            let mut e = WireExplorer::of_wire(w);
                            let mut ne = 0usize;
                            let mut pts: Vec<GpPnt> = Vec::new();
                            while e.more() {
                                ne += 1;
                                if let Some(v) = e.current_vertex() {
                                    pts.push(BRepTool::vertex_point(v));
                                }
                                e.next();
                                if ne > 512 {
                                    break;
                                }
                            }
                            println!(
                                "W {} {} {} ne={} closed={}",
                                fi,
                                wi,
                                join_name(join),
                                ne,
                                if w.0.is_closed() { 1 } else { 0 }
                            );
                            for (k, p) in pts.iter().enumerate() {
                                dump_vertex(&format!("V {} {} {} {}", fi, wi, k, join_name(join)), p);
                            }
                        }
                    }
                }
            }
        }
    }
    println!("PROBE END");
}
