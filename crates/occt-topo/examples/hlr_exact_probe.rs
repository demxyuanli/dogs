//! TEMPORARY probe: run the exact HLR chain
//! (`HLRBRep_Algo::Add` -> `HLRBRep_InternalAlgo::Update` ->
//! `HLRBRep_HLRToShape`) over a real STEP model and dump the resulting
//! compound sizes.
//!
//! Usage:
//!   cargo run --example hlr_exact_probe -- data/Cube.step [nbIso] [view]
//!
//! `view` is one of `xoy` (project along +Z, the default), `xoz`, `yoz`.

use std::path::PathBuf;
use std::sync::Arc;

use occt_core::gp::{GpAx2, GpDir, GpPnt};
use occt_topo::abs::ShapeType;
use occt_hlr::hlr_algo::HlrAlgoProjector;
use occt_hlr::hlrbrep_algo::HlrBRepAlgo;
use occt_hlr::hlrbrep_hlr_to_shape::HlrBRepHlrToShape;
use occt_topo::shape::TopoShape;
use occt_exchange::step::read_step_file;
use occt_topo::topo_tools::explorer_shapes;

fn projector_for(view: &str) -> HlrAlgoProjector {
    let (z, x) = match view {
        "xoz" => (
            GpDir::new(0.0, -1.0, 0.0).expect("dir"),
            GpDir::new(1.0, 0.0, 0.0).expect("dir"),
        ),
        "yoz" => (
            GpDir::new(-1.0, 0.0, 0.0).expect("dir"),
            GpDir::new(0.0, 1.0, 0.0).expect("dir"),
        ),
        _ => (
            GpDir::new(0.0, 0.0, 1.0).expect("dir"),
            GpDir::new(1.0, 0.0, 0.0).expect("dir"),
        ),
    };
    let ax2 = GpAx2::new(GpPnt::new(0.0, 0.0, 0.0), z, x).expect("ax2");
    HlrAlgoProjector::from_ax2(&ax2).expect("projector")
}

fn count(shape: &Option<TopoShape>, kind: ShapeType) -> usize {
    match shape {
        Some(s) => explorer_shapes(s, kind).len(),
        None => 0,
    }
}

fn main() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("..");
    let arg = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "data/Cube.step".into());
    let path = if PathBuf::from(&arg).is_absolute() {
        PathBuf::from(&arg)
    } else {
        root.join(&arg)
    };
    let nb_iso: i32 = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(0);
    let view = std::env::args().nth(3).unwrap_or_else(|| "xoy".into());

    let stack_kb: usize = std::env::var("PROBE_STACK_KB")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(65536);
    let handle = std::thread::Builder::new()
        .stack_size(stack_kb * 1024)
        .spawn(move || run(path, nb_iso, view))
        .expect("probe thread");
    handle.join().expect("probe thread join");
}

fn run(path: PathBuf, nb_iso: i32, view: String) {
    let file_name = path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    println!("EXACT {file_name} nbIso={nb_iso} view={view}");

    let model = read_step_file(&path.to_string_lossy()).expect("read_step_file");
    let shape: TopoShape =
        occt_exchange::brep_exchange::document_shape(&model).expect("document shape");

    let mut algo = HlrBRepAlgo::new();
    algo.base_mut()
        .set_projector(Arc::new(projector_for(&view)));
    algo.add(&shape, nb_iso);
    let update_result = algo.base_mut().update();
    if let Some(ds) = algo.base().data_structure() {
        println!(
            "ds vert={} edge={} face={} emap={} fmap={}",
            ds.nb_vertices(),
            ds.nb_edges(),
            ds.nb_faces(),
            ds.edge_map().len(),
            ds.face_map().len()
        );
    } else {
        println!("ds null");
    }
    match update_result {
        Ok(()) => println!("update ok"),
        Err(e) => {
            println!("update err: {e}");
            return;
        }
    }

    // `HLRBRep_Algo::Hide()` -> `HLRBRep_InternalAlgo::Hide()` covers the
    // `HLRBRep_Hider` pass; pass `nohide` as a fourth argument to keep the
    // pre-hider baseline (`H edge=0`).
    let do_hide = std::env::args().nth(4).as_deref() != Some("nohide");
    if do_hide {
        algo.base_mut().hide();
        println!("hide done");
        if let Some(ds) = algo.base().data_structure() {
            let mut all_hidden = 0;
            let mut partial = 0;
            for e in 1..=ds.nb_edges() {
                let st = ds.e_data(e).status();
                if st.all_hidden() {
                    all_hidden += 1;
                } else if !st.all_visible() {
                    partial += 1;
                }
            }
            println!(
                "status allhidden={all_hidden} partial={partial} of {}",
                ds.nb_edges()
            );
        }
    }

    let mut to_shape = HlrBRepHlrToShape::new(HlrBRepAlgo::from(&algo));
    let v = to_shape.v_compound();
    let h = to_shape.h_compound();
    let vo = to_shape.out_line_v_compound();
    let ho = to_shape.out_line_h_compound();
    let vr = to_shape.rgn_line_v_compound();
    let hr = to_shape.rgn_line_h_compound();
    println!(
        "V edge={} H edge={} Vout edge={} Hout edge={} Vrgn edge={} Hrgn edge={}",
        count(&v, ShapeType::Edge),
        count(&h, ShapeType::Edge),
        count(&vo, ShapeType::Edge),
        count(&ho, ShapeType::Edge),
        count(&vr, ShapeType::Edge),
        count(&hr, ShapeType::Edge)
    );
    println!(
        "V vert={} H vert={}",
        count(&v, ShapeType::Vertex),
        count(&h, ShapeType::Vertex)
    );
    println!("EXACT END");
}
