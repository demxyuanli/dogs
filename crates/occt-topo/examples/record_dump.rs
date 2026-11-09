//! RECORD v1 dump of one shape, in the exact layout
//! `specs/occt_probe/modes_record.cpp` writes, so a source-level parity run is
//! a field-by-field diff of two files rather than a comparison against a model
//! produced by an external application.
//!
//!   record_dump <input.step|.iges> <out.txt> [--defl <d>] [--angle <a>] [--statuses]
//!   record_dump <input.step|.iges> <out.txt> --rt iges <mid.iges> [--defl ..]
//!
//! `--rt iges` performs the write-side pipeline: read `input`, write it with
//! the port's `write_shape_iges`, read the written file back with
//! `read_iges_shape`, and dump the record of the re-read shape. The OCCT half
//! is `occt_probe <input> --record-rt iges <mid.iges> <out.txt>`, which runs the
//! same pipeline through `IGESControl_Writer` / `IGESControl_Reader`.
//!
//! Fields (see modes_record.cpp for the table):
//!   RECORD / STAGE / INPUT / TYPE / COUNTS / VALID / VOLUME / CENTROID /
//!   AREA / BBOX / FACE... / MESH / END

use std::fs;
use std::path::PathBuf;
use std::process::ExitCode;

use occt_core::gp::GpVec;
use occt_topo::brep_bnd_lib::shape_bnd_box;
use occt_topo::brep_check::Analyzer;
use occt_exchange::brep_exchange::document_shape;
use occt_topo::brep_gprop_full::{surface_properties, volume_properties};
use occt_topo::brep_tool::BRepTool;
use occt_exchange::iges::write_shape_iges;
use occt_exchange::iges_read::read_iges_shape;
use occt_topo::meshing::incremental_mesh::IncrementalMesh;
use occt_topo::meshing::range_splitter::{mesh_surface_type, MeshSurfaceType};
use occt_topo::shape::{Face, TopoShape};
use occt_exchange::step::{read_step_file, write_step};
use occt_topo::topo_tools::{edges_of_wire, faces_of, shape_counts, wires_of_face};
use occt_topo::tgeometry::GeometryRegistry;
use occt_topo::abs::ShapeType;

fn tname(t: MeshSurfaceType) -> &'static str {
    match t {
        MeshSurfaceType::Plane => "Plane",
        MeshSurfaceType::Cylinder => "Cylinder",
        MeshSurfaceType::Cone => "Cone",
        MeshSurfaceType::Sphere => "Sphere",
        MeshSurfaceType::Torus => "Torus",
        MeshSurfaceType::SurfaceOfRevolution => "Rev",
        MeshSurfaceType::SurfaceOfExtrusion => "Extr",
        MeshSurfaceType::BezierSurface => "Bezier",
        MeshSurfaceType::BSplineSurface => "BSpline",
        MeshSurfaceType::OffsetSurface => "Offset",
        MeshSurfaceType::OtherSurface => "Other",
    }
}

fn type_name(t: ShapeType) -> &'static str {
    match t {
        ShapeType::Vertex => "Vertex",
        ShapeType::Edge => "Edge",
        ShapeType::Wire => "Wire",
        ShapeType::Face => "Face",
        ShapeType::Shell => "Shell",
        ShapeType::Solid => "Solid",
        ShapeType::CompSolid => "CompSolid",
        ShapeType::Compound => "Compound",
        ShapeType::Shape => "Shape",
    }
}

fn g17(v: f64) -> String {
    if v == 0.0 {
        return "0".to_string();
    }
    format!("{v:.17e}")
}

fn bbox_text(shape: &TopoShape) -> [f64; 6] {
    match shape_bnd_box(shape).get() {
        Some((x0, x1, y0, y1, z0, z1)) => [x0, y0, z0, x1, y1, z1],
        None => [0.0; 6],
    }
}

fn render(
    shape: &TopoShape,
    stage: &str,
    input: &str,
    out: &PathBuf,
    defl: f64,
    ang: f64,
    statuses: bool,
    facemesh: bool,
) -> std::io::Result<()> {
    let mut s = String::new();
    s.push_str("RECORD v1\n");
    s.push_str(&format!("STAGE {stage}\n"));
    s.push_str(&format!("INPUT {input}\n"));
    s.push_str(&format!("TYPE {}\n", type_name(shape.shape_type())));

    // `TopExp::MapShapes` counts the root when its type matches, and the port's
    // own census (`all_subshapes`) already includes it.
    let counts = shape_counts(shape);
    let c = |t: ShapeType| counts.get(&t).copied().unwrap_or(0);
    s.push_str(&format!(
        "COUNTS V={} E={} W={} F={} SH={} SO={} CS={} C={}\n",
        c(ShapeType::Vertex),
        c(ShapeType::Edge),
        c(ShapeType::Wire),
        c(ShapeType::Face),
        c(ShapeType::Shell),
        c(ShapeType::Solid),
        c(ShapeType::CompSolid),
        c(ShapeType::Compound)
    ));

    let analyzer = Analyzer::new(shape.clone());
    s.push_str(&format!("VALID {}\n", if analyzer.is_valid() { 1 } else { 0 }));
    if statuses {
        // Diagnostic only (stdout, not the record): which sub-shape OCCT would
        // have accepted and the port rejects, and with which `BRepCheck_Status`.
        let reg = GeometryRegistry::global();
        for r in analyzer.results() {
            if r.is_valid() {
                continue;
            }
            let mut extra = String::new();
            if r.shape.shape_type() == ShapeType::Edge {
                let e = occt_topo::shape::Edge(r.shape.clone());
                let geom = reg.edge_geom(&e);
                let c3d = reg
                    .edge_curve(&e.0)
                    .map(|c| {
                        let kind = if c.is_line() {
                            "line"
                        } else if c.gp_circ().is_some() {
                            "circle"
                        } else if c.nurbs_degree().is_some() {
                            "bspline"
                        } else {
                            "other"
                        };
                        let (f, l) = reg.edge_parameters(&e.0);
                        let d = c.d0(f).distance(&c.d0(l));
                        let (basis, s, e2) = match c.untrimmed_basis() {
                            Some((b, s, e2)) => (b, s, e2),
                            None => (c.clone(), f, l),
                        };
                        format!(
                            "{kind}(range=[{f:.6},{l:.6}] span={d:.3e} win=[{:.6},{:.6}] per={} p={:.6} bwin=[{s:.6},{e2:.6}])",
                            c.first_parameter(),
                            c.last_parameter(),
                            basis.is_periodic(),
                            basis.period()
                        )
                    })
                    .unwrap_or_else(|| "None".to_string());
                extra = format!(
                    " degenerated={} same_parameter={} same_range={} geom={} curve3d_null={} c3d={} pcurves={}",
                    reg.is_degenerated_edge(&e),
                    reg.same_parameter(&e),
                    reg.same_range(&e),
                    geom.is_some(),
                    geom.as_ref().map(|g| g.curve3d_null).unwrap_or(true),
                    c3d,
                    geom.as_ref().map(|g| g.pcurves.len()).unwrap_or(0)
                );
            }
            println!(
                "INVALID type={} statuses={:?}{extra}",
                type_name(r.shape.shape_type()),
                r.statuses
            );
        }
    }

    match volume_properties(shape) {
        Ok(props) => {
            let c = props.center();
            s.push_str(&format!("VOLUME {}\n", g17(props.mass())));
            s.push_str(&format!(
                "CENTROID {} {} {}\n",
                g17(c.x()),
                g17(c.y()),
                g17(c.z())
            ));
        }
        Err(err) => {
            s.push_str("VOLUME nan\n");
            s.push_str(&format!("CENTROID nan nan nan # {err}\n"));
        }
    }
    match surface_properties(shape) {
        Ok((props, _)) => s.push_str(&format!("AREA {}\n", g17(props.mass()))),
        Err(err) => s.push_str(&format!("AREA nan # {err}\n")),
    }

    let b = bbox_text(shape);
    s.push_str(&format!(
        "BBOX {} {} {} {} {} {}\n",
        g17(b[0]),
        g17(b[1]),
        g17(b[2]),
        g17(b[3]),
        g17(b[4]),
        g17(b[5])
    ));

    // Face rows are compared as a multiset, so render then sort: traversal
    // order is not part of the record.
    let mut rows: Vec<String> = Vec::new();
    for face in faces_of(shape) {
        let ty = BRepTool::face_surface(&face)
            .as_ref()
            .map(|s| tname(mesh_surface_type(s.as_ref())))
            .unwrap_or("Other");
        let area = surface_properties(&face.0).map(|(p, _)| p.mass()).unwrap_or(f64::NAN);
        let wires = wires_of_face(&face);
        let edges: usize = wires.iter().map(|w| edges_of_wire(w).len()).sum();
        let fb = bbox_text(&face.0);
        rows.push(format!(
            "FACE type={ty} area={} wires={} edges={} bbox={} {} {} {} {} {}",
            g17(area),
            wires.len(),
            edges,
            g17(fb[0]),
            g17(fb[1]),
            g17(fb[2]),
            g17(fb[3]),
            g17(fb[4]),
            g17(fb[5])
        ));
    }
    rows.sort();
    for r in &rows {
        s.push_str(r);
        s.push('\n');
    }

    let inc = IncrementalMesh::from_deflection(shape, defl, false, ang);
    if facemesh {
        // Diagnostic only (stdout, not the record): per-face node/triangle
        // counts, the port half of the probe's `--mesh` list.
        for st in inc.face_stats() {
            println!(
                "FACE idx={} type={} nodes={} triangles={}",
                st.index,
                tname(st.surface),
                st.vertices,
                st.triangles
            );
        }
        println!("TOTAL faces={} nodes={} triangles={}", inc.face_stats().len(),
            inc.face_stats().iter().map(|s| s.vertices).sum::<usize>(),
            inc.face_stats().iter().map(|s| s.triangles).sum::<usize>());
    }
    let mut nv = 0usize;
    let mut nf = 0usize;
    let mut area = 0.0f64;
    for ft in inc.face_triangulations() {
        nv += ft.vertices.len();
        nf += ft.triangles.len();
        for t in &ft.triangles {
            let (Some(p1), Some(p2), Some(p3)) = (
                ft.vertices.get(t.n0),
                ft.vertices.get(t.n1),
                ft.vertices.get(t.n2),
            ) else {
                continue;
            };
            let v1 = GpVec::from_pnts(p1, p2);
            let v2 = GpVec::from_pnts(p1, p3);
            area += 0.5 * v1.cross_magnitude(&v2);
        }
    }
    s.push_str(&format!(
        "MESH nv={nv} nf={nf} area={} defl={} angle={}\n",
        g17(area),
        g17(defl),
        g17(ang)
    ));
    s.push_str("END\n");
    fs::write(out, s)?;
    println!("RECORD {stage} faces={} out={}", rows.len(), out.display());
    Ok(())
}

/// Diagnostic only (stdout, not the record): the `BRep_TEdge` representation
/// census `BRepCheck_Edge::Minimum` (`BRepCheck_Edge.cxx:88-148`) reads, in the
/// same shape as the probe's `--edgerep`, so the two sides can be diffed
/// edge by edge.
fn edge_rep_census(shape: &TopoShape) {
    let reg = GeometryRegistry::global();
    let mut n = 0usize;
    for s in occt_topo::topo_tools::all_subshapes(shape) {
        if s.shape_type() != ShapeType::Edge {
            continue;
        }
        n += 1;
        let e = occt_topo::shape::Edge(s);
        let geom = reg.edge_geom(&e);
        let (f, l) = reg.edge_parameters(&e.0);
        println!(
            "EDGEREP n={} deg={} spar={} srange={} c3drep={} c3dnull={} cos={} range=[{f:.6},{l:.6}]",
            n,
            reg.is_degenerated_edge(&e) as u8,
            reg.same_parameter(&e) as u8,
            reg.same_range(&e) as u8,
            geom.is_some() as u8,
            reg.edge_curve(&e.0).is_none() as u8,
            geom.as_ref().map(|g| g.pcurves.len()).unwrap_or(0)
        );
    }
}

/// Diagnostic only (stdout, not the record): every vertex with its point, its
/// tolerance and the (wire, edge) slots that reference it, so a junction that
/// one side shares and the other leaves split shows up as an identity column.
fn vertex_rep_census(shape: &TopoShape) {
    let mut verts: Vec<TopoShape> = Vec::new();
    let mut wires: Vec<TopoShape> = Vec::new();
    let mut eds: Vec<TopoShape> = Vec::new();
    for s in occt_topo::topo_tools::all_subshapes(shape) {
        match s.shape_type() {
            ShapeType::Vertex => verts.push(s),
            ShapeType::Wire => wires.push(s),
            ShapeType::Edge => eds.push(s),
            _ => {}
        }
    }
    let slot = |v: Option<occt_topo::shape::Vertex>, order: &[TopoShape]| match v {
        Some(v) => order
            .iter()
            .position(|o| occt_topo::topo_tools::is_same(o, &v.0))
            .map(|i| (i + 1) as i64)
            .unwrap_or(-1),
        None => -1,
    };
    for (i, v) in verts.iter().enumerate() {
        let vert = occt_topo::shape::Vertex(v.clone());
        let p = BRepTool::vertex_point(&vert);
        println!(
            "VERTEX n={} p=[{:.17e},{:.17e},{:.17e}] tol={:.17e}",
            i + 1,
            p.x(),
            p.y(),
            p.z(),
            BRepTool::vertex_tolerance(&vert)
        );
    }
    for (i, w) in wires.iter().enumerate() {
        let wire = occt_topo::shape::Wire(w.clone());
        let ids: Vec<String> = edges_of_wire(&wire)
            .iter()
            .map(|e| {
                eds.iter()
                    .position(|o| occt_topo::topo_tools::is_same(o, &e.0))
                    .map(|k| (k + 1).to_string())
                    .unwrap_or_else(|| "?".to_string())
            })
            .collect();
        println!("WIRE n={} edges=[{}]", i + 1, ids.join(","));
    }
    for (i, e) in eds.iter().enumerate() {
        let edge = occt_topo::shape::Edge(e.clone());
        println!(
            "EDGE n={} first={} last={}",
            i + 1,
            slot(occt_topo::shhealing::first_vertex(&edge), &verts),
            slot(occt_topo::shhealing::last_vertex(&edge), &verts)
        );
    }
}

/// Diagnostic only (stdout, not the record): the composed orientation of every
/// face and of the shell/solid above it, in the probe's default-mode shape
/// (`modes_default.cpp`), so the two sides can be diffed shape by shape.
/// `--orient raw` reads the stored child views instead of the composed ones.
fn orient_census(shape: &TopoShape, raw: bool) {
    let mut shell = 0usize;
    for s in occt_topo::topo_tools::all_subshapes(shape) {
        match s.shape_type() {
            ShapeType::Solid => println!("SOLID orient={}", s.orientation() as i32),
            ShapeType::Shell => {
                println!("SHELL orient={}", s.orientation() as i32);
                shell += 1;
            }
            _ => {}
        }
    }
    let _ = shell;
    if raw {
        for (i, f) in stored_faces(shape).iter().enumerate() {
            println!("RAWFACE n={} ori={}", i + 1, f.orientation() as i32);
        }
        return;
    }
    for (i, f) in faces_of(shape).iter().enumerate() {
        let ty = BRepTool::face_surface(f)
            .as_ref()
            .map(|s| tname(mesh_surface_type(s.as_ref())))
            .unwrap_or("Other");
        println!("FACE n={} ori={} type={ty}", i + 1, f.orientation() as i32);
    }
}

/// Faces read straight out of the stored `TShape::children` lists, i.e. the
/// view `TopoDS_Builder::Add` wrote (`TopoDS_Builder.cxx:74-92`), without the
/// composition `TopoDS_Iterator.cxx:74-80` applies on read.
fn stored_faces(shape: &TopoShape) -> Vec<Face> {
    let mut out = Vec::new();
    let mut stack = vec![shape.clone()];
    while let Some(s) = stack.pop() {
        if s.shape_type() == ShapeType::Face {
            out.push(Face(s));
            continue;
        }
        let children = s
            .tshape
            .read()
            .expect("poisoned TShape lock")
            .children
            .clone();
        for c in children.into_iter().rev() {
            stack.push(c);
        }
    }
    out
}

/// Diagnostic only (stdout, not the record): `BRepClass3d_SolidClassifier::
/// PerformInfinitePoint` on a solid built from each shell as stored, and on a
/// solid built from the same shell reversed. `ShapeFix_Solid::SolidFromShell`
/// (`ShapeFix_Solid.cxx:655-702`) reverses exactly when the stored shell comes
/// back `IN`, so the `stored` column must read OUT on an outward shell.
fn inf_solid_census(shape: &TopoShape) {
    let builder = occt_topo::builder::TopoBuilder::new();
    let mut i = 0usize;
    for sh in occt_topo::topo_tools::all_subshapes(shape) {
        if sh.shape_type() != ShapeType::Shell {
            continue;
        }
        for reversed in [false, true] {
            let mut s = sh.clone();
            if reversed {
                s.reverse();
            }
            let solid = builder.make_solid(&[occt_topo::shape::Shell(s.clone())]).0;
            let mut sc = occt_topo::brep_class3d::SolidClassifier::new();
            sc.load(solid);
            sc.perform_infinite_point(occt_core::precision::CONFUSION);
            println!(
                "INFSOLID shell={i} ori={} rev={} state={}",
                s.orientation() as i32,
                reversed as u8,
                sc.state() as i32
            );
        }
        i += 1;
    }
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let get = |flag: &str| -> Option<String> {
        args.iter()
            .position(|a| a == flag)
            .and_then(|i| args.get(i + 1))
            .cloned()
    };
    let rt = get("--rt");
    let mid = get("--mid");
    let defl: f64 = get("--defl").and_then(|v| v.parse().ok()).unwrap_or(0.1);
    let ang: f64 = get("--angle").and_then(|v| v.parse().ok()).unwrap_or(0.5);
    let statuses = args.iter().any(|a| a == "--statuses");
    let facemesh = args.iter().any(|a| a == "--facemesh");
    let edge_dump = args.iter().any(|a| a == "--edgerep");
    let vertex_dump = args.iter().any(|a| a == "--vertexrep");
    let orient = args.iter().any(|a| a == "--orient");
    let raw_orient = args.iter().any(|a| a == "--orientraw");
    let orient = orient || raw_orient;
    let infsolid = args.iter().any(|a| a == "--infsolid");
    let mut skip = vec![
        "--defl".to_string(),
        "--angle".to_string(),
        "--rt".to_string(),
        "--mid".to_string(),
        "--statuses".to_string(),
        "--facemesh".to_string(),
        "--edgerep".to_string(),
        "--vertexrep".to_string(),
        "--orient".to_string(),
        "--orientraw".to_string(),
        "--infsolid".to_string(),
    ];
    for f in ["--defl", "--angle", "--rt", "--mid"] {
        if let Some(v) = get(f) {
            skip.push(v);
        }
    }
    let pos: Vec<&String> = args.iter().filter(|a| !skip.iter().any(|s| s == *a)).collect();
    if pos.len() != 2 && !((edge_dump || orient || infsolid || vertex_dump) && pos.len() == 1) {
        eprintln!(
            "usage: record_dump <input.step|.iges> <out.txt> [--rt iges] [--mid <path>] [--defl d] [--angle a]"
        );
        return ExitCode::from(2);
    }
    let input = pos[0].clone();
    let out = pos.get(1).map(|s| PathBuf::from(s.as_str()));
    let lower = input.to_ascii_lowercase();
    let is_step = lower.ends_with(".step") || lower.ends_with(".stp");

    let (shape, model) = if is_step {
        match read_step_file(&input) {
            Ok(m) => match document_shape(&m) {
                Some(s) => (s, Some(m)),
                None => {
                    eprintln!("record_dump: no root shapes in {input}");
                    return ExitCode::from(1);
                }
            },
            Err(err) => {
                eprintln!("record_dump: {err}");
                return ExitCode::from(1);
            }
        }
    } else {
        match read_iges_shape(&input) {
            Ok(s) => (s, None),
            Err(err) => {
                eprintln!("record_dump: {err}");
                return ExitCode::from(1);
            }
        }
    };

    let (stage, read_back) = match rt.as_deref() {
        None => ("read", shape.clone()),
        Some("iges") => {
            let mid = match mid {
                Some(m) => PathBuf::from(m),
                None => {
                    eprintln!("record_dump: --rt iges needs --mid <path>");
                    return ExitCode::from(2);
                }
            };
            if let Err(err) = fs::write(&mid, write_shape_iges(&shape)) {
                eprintln!("record_dump: write {}: {err}", mid.display());
                return ExitCode::from(1);
            }
            match read_iges_shape(&mid.to_string_lossy()) {
                Ok(s) => ("rt-iges", s),
                Err(err) => {
                    eprintln!("record_dump: re-read {}: {err}", mid.display());
                    return ExitCode::from(1);
                }
            }
        }
        Some("step") => {
            let model = match model {
                Some(m) => m,
                None => {
                    eprintln!("record_dump: --rt step needs a STEP input");
                    return ExitCode::from(2);
                }
            };
            let mid = match mid {
                Some(m) => PathBuf::from(m),
                None => {
                    eprintln!("record_dump: --rt step needs --mid <path>");
                    return ExitCode::from(2);
                }
            };
            if let Err(err) = fs::write(&mid, write_step(&model)) {
                eprintln!("record_dump: write {}: {err}", mid.display());
                return ExitCode::from(1);
            }
            match read_step_file(&mid.to_string_lossy()) {
                Ok(m) => match document_shape(&m) {
                    Some(s) => ("rt-step", s),
                    None => {
                        eprintln!("record_dump: re-read {} produced no shapes", mid.display());
                        return ExitCode::from(1);
                    }
                },
                Err(err) => {
                    eprintln!("record_dump: re-read {}: {err}", mid.display());
                    return ExitCode::from(1);
                }
            }
        }
        Some(other) => {
            eprintln!("record_dump: unknown --rt {other}");
            return ExitCode::from(2);
        }
    };

    if edge_dump {
        edge_rep_census(&read_back);
        if out.is_none() && !orient {
            return ExitCode::SUCCESS;
        }
    }
    if vertex_dump {
        vertex_rep_census(&read_back);
        if out.is_none() && !orient && !edge_dump {
            return ExitCode::SUCCESS;
        }
    }
    if orient {
        orient_census(&read_back, raw_orient);
        if out.is_none() && !infsolid {
            return ExitCode::SUCCESS;
        }
    }
    if infsolid {
        inf_solid_census(&read_back);
        if out.is_none() {
            return ExitCode::SUCCESS;
        }
    }
    let out = match out {
        Some(o) => o,
        None => unreachable!("guarded by the argument count check above"),
    };
    match render(&read_back, stage, &input, &out, defl, ang, statuses, facemesh) {
        Ok(()) => ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("record_dump: write {}: {err}", out.display());
            ExitCode::from(1)
        }
    }
}
