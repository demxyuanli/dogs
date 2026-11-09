//! TEMPORARY probe: run the polygonal HLR chain
//! (`HLRBRep_PolyAlgo::Update` -> `HLRBRep_PolyHLRToShape::Update` ->
//! `InternalCompound`) over a real STEP model and dump the resulting segment
//! sets.
//!
//! Usage:
//!   cargo run --example hlr_poly_probe -- data/Cube.step [deflection] [view]
//!
//! `view` is one of `xoy` (project along +Z, the default), `xoz`, `yoz`.
//!
//! Delete this file once the polygonal HLR chain has a stored reference.

use std::path::PathBuf;

use occt_core::gp::{GpAx2, GpDir, GpPnt};
use occt_topo::hlr_algo::HlrAlgoProjector;
use occt_topo::hlrbrep_poly_algo::HlrBRepPolyAlgo;
use occt_topo::hlrbrep_poly_hlr_to_shape::HlrBRepPolyHlrToShape;
use occt_topo::meshing::incremental_mesh::IncrementalMesh;
use occt_topo::shape::TopoShape;
use occt_topo::step::read_step_file;

fn main() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..").join("..");
    let arg = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "data/Cube.step".into());
    let path = if PathBuf::from(&arg).is_absolute() {
        PathBuf::from(&arg)
    } else {
        root.join(&arg)
    };
    let deflection: f64 = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(0.1);
    let view = std::env::args().nth(3).unwrap_or_else(|| "xoy".into());

    let stack_kb: usize = std::env::var("PROBE_STACK_KB")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(65536);
    let handle = std::thread::Builder::new()
        .stack_size(stack_kb * 1024)
        .spawn(move || run(path, deflection, view))
        .expect("probe thread");
    handle.join().expect("probe thread join");
}

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
    // Same frame as the C++ probe's `gp_Ax2(gp_Pnt(0,0,0), dir, vx)`.
    let ax2 = GpAx2::new(GpPnt::new(0.0, 0.0, 0.0), z, x).expect("ax2");
    HlrAlgoProjector::from_ax2(&ax2).expect("projector")
}

fn run(path: PathBuf, deflection: f64, view: String) {
    let file_name = path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    println!("PROBE {file_name} lin={deflection} view={view}");

    let model = read_step_file(&path.to_string_lossy()).expect("read_step_file");
    let shape: TopoShape = occt_topo::brep_exchange::document_shape(&model).expect("document shape");

    let inc = IncrementalMesh::from_deflection(&shape, deflection, false, 0.3490658503988659);
    let tris = inc.face_triangulations();
    let n_tri: usize = tris.iter().map(|t| t.triangles.len()).sum();
    println!("mesh faces={} triangles={n_tri}", tris.len());

    let mut committed = 0usize;
    let mut uncommitted = 0usize;
    for f in occt_topo::topo_tools::explorer_shapes(&shape, occt_topo::abs::ShapeType::Face) {
        if occt_topo::brep_tool::BRepTool::face_triangulation(&occt_topo::shape::Face(f)).is_some() {
            committed += 1;
        } else {
            uncommitted += 1;
        }
    }
    println!("committed_faces={committed} uncommitted_faces={uncommitted}");

    let mut algo = HlrBRepPolyAlgo::from_shape(shape.clone());
    let proj = projector_for(&view);
    algo.set_projector(proj);
    algo.update();
    println!("poly_algo more_hide={}", algo.more_hide());

    // Intermediate dump of every shell / face / triangle / bi-point, so the
    // OCCT transcript can be diffed stage by stage.
    {
        let shells = algo.algo().poly_shell();
        println!("shells={}", shells.len());
        for (si, shell) in shells.iter().enumerate() {
            println!(
                "SHELL {} hiding={} npoly={} nhidpoly={}",
                si + 1,
                u8::from(shell.hiding()),
                shell.poly_data().len(),
                shell.nb_hiding_poly_data()
            );
            for (fi, pd) in shell.poly_data().iter().enumerate() {
                println!(
                    " FACE {} nodes={} tris={} hid={}",
                    pd.face_index(),
                    pd.nodes().len(),
                    pd.t_data().len(),
                    u8::from(pd.hiding())
                );
                for (ni, n) in pd.nodes().iter().enumerate() {
                    println!("  NODE {} ({:.9},{:.9},{:.9})", ni + 1, n.x(), n.y(), n.z());
                }
                for (ti, t) in pd.t_data().iter().enumerate() {
                    println!(
                        "  TRI {} n=({},{},{}) flags=0x{:x}",
                        ti + 1,
                        t.node1,
                        t.node2,
                        t.node3,
                        t.flags as u32
                    );
                }
                if let Some(ph) = pd.ph_dat() {
                    for (hi, h) in ph.iter().enumerate() {
                        let i = h.indices();
                        let p = h.plane();
                        println!(
                            "  HIDE {} tri={} min={} max={} plane=({:.9},{:.9},{:.9},{:.9})",
                            hi + 1,
                            i.index,
                            i.min,
                            i.max,
                            p.normal.x(),
                            p.normal.y(),
                            p.normal.z(),
                            p.d
                        );
                    }
                }
                let _ = fi;
            }
            for (bi, bp) in shell.edges().iter().enumerate() {
                let p = bp.points();
                let i = bp.indices();
                println!(
                    " BP {} p1=({:.9},{:.9},{:.9}) p2=({:.9},{:.9},{:.9}) t1=({:.9},{:.9},{:.9}) t2=({:.9},{:.9},{:.9}) idx={} c1={} f1=({},{}) c2={} f2=({},{}) minseg={} maxseg={} hid={} rg1={} rgn={} outl={} intl={}",
                    bi,
                    p.pnt1.x(), p.pnt1.y(), p.pnt1.z(),
                    p.pnt2.x(), p.pnt2.y(), p.pnt2.z(),
                    p.pntp1.x(), p.pntp1.y(), p.pntp1.z(),
                    p.pntp2.x(), p.pntp2.y(), p.pntp2.z(),
                    i.shape_index,
                    i.face_conex1, i.face1_pt1, i.face1_pt2,
                    i.face_conex2, i.face2_pt1, i.face2_pt2,
                    i.min_seg,
                    i.max_seg,
                    u8::from(bp.hidden()),
                    u8::from(bp.rg1_line()),
                    u8::from(bp.rgn_line()),
                    u8::from(bp.out_line()),
                    u8::from(bp.int_line())
                );
            }
        }
    }

    // Edge polygon-on-triangulation dump, read back exactly the way
    // `InitBiPointsWithConnexity` reads it (HLRBRep_PolyAlgo.cxx:900-910).
    {
        let e_map = occt_topo::topo_tools::explorer_shapes(&shape, occt_topo::abs::ShapeType::Edge);
        let f_map = occt_topo::topo_tools::explorer_shapes(&shape, occt_topo::abs::ShapeType::Face);
        // Deduplicate the way `TopExp::MapShapes` does, so the printed `edge=`
        // ids match the C++ probe's `TopTools_IndexedMapOfShape` numbering.
        let emap_keys: Vec<_> = {
            let mut keys = Vec::new();
            for e in &e_map {
                let k = occt_topo::topo_tools::shape_map_key(e);
                if !keys.contains(&k) {
                    keys.push(k);
                }
            }
            keys
        };
        let edge_id = |e: &occt_topo::shape::TopoShape| -> usize {
            let k = occt_topo::topo_tools::shape_map_key(e);
            emap_keys.iter().position(|x| *x == k).map(|p| p + 1).unwrap_or(0)
        };
        println!("emap={} fmap={}", emap_keys.len(), f_map.len());
        for (fi, f) in f_map.iter().enumerate() {
            let face = occt_topo::shape::Face(f.clone());
            let tri = match occt_topo::brep_tool::BRepTool::face_triangulation(&face) {
                Some(t) => t,
                None => {
                    println!("POLY face={} no triangulation", fi + 1);
                    continue;
                }
            };
            for e in occt_topo::topo_tools::explorer_shapes(f, occt_topo::abs::ShapeType::Edge) {
                let ei = edge_id(&e);
                // `BRep_Tool::PolygonOnTriangulation` (`BRep_Tool.cxx:683-713`)
                // picks `PolygonOnTriangulation2` when the edge occurrence is
                // REVERSED, exactly as the occurrence that the explorer yields.
                let reversed = e.orientation().is_reversed();
                match occt_topo::brep_tool::BRepTool::polygon_on_triangulation(
                    &occt_topo::shape::Edge(e.clone()),
                    &tri,
                    reversed,
                ) {
                    Some(p) => {
                        let nodes: Vec<String> =
                            p.nodes.iter().map(|n| n.to_string()).collect();
                        let params: Vec<String> = p
                            .parameters
                            .as_ref()
                            .map(|ps| ps.iter().map(|v| format!("{v:.9}")).collect())
                            .unwrap_or_default();
                        println!(
                            "POLY face={} edge={} n={} nodes={} params={}",
                            fi + 1,
                            ei,
                            p.nodes.len(),
                            nodes.join(","),
                            params.join(",")
                        );
                    }
                    None => println!("POLY face={} edge={} null", fi + 1, ei),
                }
            }
        }
    }

    let mut out = HlrBRepPolyHlrToShape::new();
    out.update(&algo);
    let vis = out.bi_pnt_vis();
    let hid = out.bi_pnt_hid();
    println!("bi_pnt_vis={} bi_pnt_hid={}", vis.len(), hid.len());

    for (name, compound) in [
        ("VCompound", out.v_compound()),
        ("Rg1LineVCompound", out.rg1_line_v_compound()),
        ("RgNLineVCompound", out.rgn_line_v_compound()),
        ("OutLineVCompound", out.out_line_v_compound()),
        ("HCompound", out.h_compound()),
        ("Rg1LineHCompound", out.rg1_line_h_compound()),
        ("RgNLineHCompound", out.rgn_line_h_compound()),
        ("OutLineHCompound", out.out_line_h_compound()),
    ] {
        match compound {
            Some(c) => println!(
                "{name}: edges={}",
                occt_topo::topo_tools::explorer_shapes(&c, occt_topo::abs::ShapeType::Edge).len()
            ),
            None => println!("{name}: null"),
        }
    }

    for (kind, list) in [("VIS", vis), ("HID", hid)] {
        for seg in list {
            println!(
                "{kind} ({:.9},{:.9})-({:.9},{:.9}) rg1={} rgn={} outl={} intl={}",
                seg.p1().x(),
                seg.p1().y(),
                seg.p2().x(),
                seg.p2().y(),
                seg.rg1_line(),
                seg.rgn_line(),
                seg.out_line(),
                seg.int_line()
            );
        }
    }
    println!("PROBE END");
}
