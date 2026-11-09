//! TEMPORARY probe: dump `HLRBRep_Curve` state and evaluation over the edges
//! of a real STEP model, for the line-by-line diff against the OCCT 8.0 C++
//! generator (`specs/hlr_probe/hlr_probe.cpp`, `curve` mode).
//!
//! Usage:
//!   cargo run --example hlrbrep_curve_probe -- data/Cube.step [xoy|xoz|yoz]
//!
//! Delete this file once the exact-HLR curve layer has a stored reference.

use std::path::PathBuf;
use std::sync::Arc;

use occt_core::gp::{GpAx2, GpDir, GpPnt};
use occt_topo::geom_abs_curve_type::GeomAbsCurveType;
use occt_topo::hlr_algo::HlrAlgoProjector;
use occt_topo::hlrbrep_curve::HlrBRepCurve;
use occt_topo::shape::{Edge, TopoShape};
use occt_topo::step::read_step_file;

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

fn type_index(t: GeomAbsCurveType) -> i32 {
    match t {
        GeomAbsCurveType::Line => 0,
        GeomAbsCurveType::Circle => 1,
        GeomAbsCurveType::Ellipse => 2,
        GeomAbsCurveType::Hyperbola => 3,
        GeomAbsCurveType::Parabola => 4,
        GeomAbsCurveType::BezierCurve => 5,
        GeomAbsCurveType::BSplineCurve => 6,
        GeomAbsCurveType::OffsetCurve => 7,
        GeomAbsCurveType::OtherCurve => 8,
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
    let view = std::env::args().nth(2).unwrap_or_else(|| "xoy".into());

    let stack_kb: usize = std::env::var("PROBE_STACK_KB")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(65536);
    let handle = std::thread::Builder::new()
        .stack_size(stack_kb * 1024)
        .spawn(move || run(path, view))
        .expect("probe thread");
    handle.join().expect("probe thread join");
}

fn run(path: PathBuf, view: String) {
    let file_name = path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    println!("PROBE {file_name} view={view}");

    let model = read_step_file(&path.to_string_lossy()).expect("read_step_file");
    let shape: TopoShape =
        occt_topo::brep_exchange::document_shape(&model).expect("document shape");

    // Deduplicate the way `TopExp::MapShapes` does, so the `EDGE` ids match
    // the C++ probe's `TopTools_IndexedMapOfShape` numbering.
    let mut keys = Vec::new();
    let mut edges: Vec<TopoShape> = Vec::new();
    for e in occt_topo::topo_tools::explorer_shapes(&shape, occt_topo::abs::ShapeType::Edge) {
        let k = occt_topo::topo_tools::shape_map_key(&e);
        if !keys.contains(&k) {
            keys.push(k);
            edges.push(e);
        }
    }
    println!("emap={}", edges.len());

    let proj = Arc::new(projector_for(&view));
    for (i, e) in edges.iter().enumerate() {
        let edge = Edge(e.clone());
        let mut c = HlrBRepCurve::new();
        c.set_projector(Some(proj.clone()));
        c.set_curve(&edge);

        let mut mn = [f64::MAX; 16];
        let mut mx = [-f64::MAX; 16];
        let tol = c.update(&mut mn, &mut mx);
        println!(
            "EDGE {} type={} tol={:.9} first={:.9} last={:.9}",
            i + 1,
            type_index(c.get_type()),
            tol,
            c.first_parameter(),
            c.last_parameter()
        );
        print!(" MIN");
        for v in mn {
            print!(" {:.9}", v);
        }
        println!();
        print!(" MAX");
        for v in mx {
            print!(" {:.9}", v);
        }
        println!();

        let (u0, u1) = occt_topo::brep_tool::BRepTool::edge_parameters(&edge);
        for s in 0..=4 {
            let u = u0 + (u1 - u0) * (s as f64) / 4.0;
            let p3 = c.value3d(u);
            let p2 = c.d0(u);
            let (p1, v1) = c.d1(u);
            let (p2d, v21, v22) = c.d2(u);
            let z = c.z(u);
            let par2 = c.parameter2d(u);
            let par3 = c.parameter3d(par2);
            println!(
                " S{} u={:.9} p3=({:.9},{:.9},{:.9}) d0=({:.9},{:.9}) d1p=({:.9},{:.9}) d1v=({:.9},{:.9}) d2p=({:.9},{:.9}) d2v1=({:.9},{:.9}) d2v2=({:.9},{:.9}) z={:.9} par2={:.9} par3={:.9}",
                s, u,
                p3.x(), p3.y(), p3.z(),
                p2.x(), p2.y(),
                p1.x(), p1.y(),
                v1.x(), v1.y(),
                p2d.x(), p2d.y(),
                v21.x(), v21.y(),
                v22.x(), v22.y(),
                z, par2, par3
            );
        }
        match c.tangent(true) {
            Some((p, d)) => println!(" TAN1 p=({:.9},{:.9}) d=({:.9},{:.9})", p.x(), p.y(), d.x, d.y),
            None => println!(" TAN1 undef"),
        }
        match c.tangent(false) {
            Some((p, d)) => println!(" TAN2 p=({:.9},{:.9}) d=({:.9},{:.9})", p.x(), p.y(), d.x, d.y),
            None => println!(" TAN2 undef"),
        }
        if c.get_type() == GeomAbsCurveType::Line {
            match c.line() {
                Some(l) => println!(
                    " LINE p=({:.9},{:.9}) d=({:.9},{:.9})",
                    l.location().x(),
                    l.location().y(),
                    l.direction().x,
                    l.direction().y
                ),
                None => println!(" LINE null"),
            }
        }
    }
    println!("PROBE END");
}
