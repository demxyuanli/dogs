//! TEMPORARY probe (T-111).
//!
//! Mirrors `specs/occt_probe/occt_probe.cpp --igespc`: for every face, in
//! `children_of_type(Face)` order (the order `write_iges` transfers them, which
//! matches OCCT's `TopExp_Explorer(aShape, TopAbs_FACE)`), print the basis
//! surface type / U-V bounds and, per edge, the 3D curve type plus the *stored*
//! p-curve type and range - the same facts `BRepToIGES_BRShell::TransferFace`
//! sees through `BRep_Tool::CurveOnSurface`.
//!
//! Usage:
//!   cargo run --example iges_pcurve_dump -- data/Shape-1.step [out.txt]

use std::io::Write;
use std::path::PathBuf;
use std::sync::Arc;

use occt_geom::{Curve, Surface};
use occt_geom2d::curve::Curve2d;
use occt_topo::abs::ShapeType;
use occt_topo::brep_tool::BRepTool;
use occt_topo::shape::{Edge, Face, TopoShape};
use occt_topo::iges_read::read_iges_shape;
use occt_topo::step::read_step_file;
use occt_topo::topo_tools::{edges_of_wire, wires_of_face};

fn c3_kind(c: &dyn Curve) -> String {
    if c.is_geom_trimmed() {
        "Geom_TrimmedCurve".into()
    } else if c.trimmed_basis_range().is_some() {
        "Geom_TrimmedCurve".into()
    } else if c.gp_circ().is_some() {
        "Geom_Circle".into()
    } else if c.gp_line().is_some() {
        "Geom_Line".into()
    } else if c.gp_ellipse().is_some() {
        "Geom_Ellipse".into()
    } else if c.bezier_poles().is_some() {
        "Geom_BezierCurve".into()
    } else if c.offset_curve().is_some() {
        "Geom_OffsetCurve".into()
    } else if c.bspline_poles().is_some() {
        "Geom_BSplineCurve".into()
    } else {
        "?".into()
    }
}

/// Faces of `shape` in the order `write_iges` transfers them: the same
/// depth-first descent as `TopExp_Explorer(aShape, TopAbs_FACE)`
/// (`solid -> shell -> face`), matching the oracle's `--igespc` face indices.
fn collect_faces(shape: &TopoShape, out: &mut Vec<TopoShape>) {
    let kids: Vec<TopoShape> = shape.tshape.read().unwrap().children.clone();
    for k in kids {
        match k.shape_type() {
            ShapeType::Face => out.push(k),
            ShapeType::Shell | ShapeType::Solid | ShapeType::CompSolid | ShapeType::Compound => {
                collect_faces(&k, out)
            }
            _ => {}
        }
    }
}

fn c2_kind(c: &dyn Curve2d) -> String {
    if c.trimmed_basis().is_some() {
        "Geom2d_TrimmedCurve".into()
    } else if c.gp_lin2d().is_some() {
        "Geom2d_Line".into()
    } else if c.gp_circ2d().is_some() {
        "Geom2d_Circle".into()
    } else if c.gp_elips2d().is_some() {
        "Geom2d_Ellipse".into()
    } else if c.bspline_poles2d().is_some() {
        "Geom2d_BSplineCurve".into()
    } else {
        "?".into()
    }
}

fn surf_kind(s: &dyn Surface) -> String {
    if s.gp_pln().is_some() {
        "Geom_Plane".into()
    } else if s.gp_cylinder().is_some() {
        "Geom_CylindricalSurface".into()
    } else if s.gp_cone().is_some() {
        "Geom_ConicalSurface".into()
    } else if s.gp_sphere().is_some() {
        "Geom_SphericalSurface".into()
    } else if s.gp_torus().is_some() {
        "Geom_ToroidalSurface".into()
    } else if s.is_surface_of_revolution() {
        "Geom_SurfaceOfRevolution".into()
    } else if s.is_surface_of_linear_extrusion() {
        "Geom_SurfaceOfLinearExtrusion".into()
    } else if s.is_bspline_surface() {
        "Geom_BSplineSurface".into()
    } else {
        "?".into()
    }
}

fn main() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..").join("..");
    let args: Vec<String> = std::env::args().skip(1).collect();
    let name = args.first().map_or("Shape-1", |s| s.as_str());
    let path = if name.contains('/') {
        root.join(name)
    } else {
        root.join("data").join(format!("{name}.step"))
    };
    let shape: TopoShape = if name.ends_with(".iges") || name.ends_with(".igs") {
        read_iges_shape(&path.to_string_lossy()).expect("read iges")
    } else {
        let model = read_step_file(&path.to_string_lossy()).expect("read step");
        model.shapes[0].shape.clone()
    };
    let mut out = String::new();
    let mut faces: Vec<TopoShape> = Vec::new();
    if shape.shape_type() == ShapeType::Face {
        faces.push(shape.clone());
    } else {
        collect_faces(&shape, &mut faces);
    }
    for (fi, f) in faces.iter().enumerate() {
        let face = Face(f.clone());
        let surf = BRepTool::face_surface(&face).expect("face surface");
        let (u1, u2) = surf.u_range();
        let (v1, v2) = surf.v_range();
        out.push_str(&format!(
            "FACE[{fi}] {} surf={} u=[{u1:.6},{u2:.6}] v=[{v1:.6},{v2:.6}] uclosed={} vclosed={}\n",
            if face.0.orientation().is_reversed() {
                "REV"
            } else {
                "FWD"
            },
            surf_kind(surf.as_ref()),
            surf.is_u_closed() as u8,
            surf.is_v_closed() as u8,
        ));
        let mut ei = 0usize;
        for (wi, w) in wires_of_face(&face).iter().enumerate() {
            let we = edges_of_wire(w);
            out.push_str(&format!("  WIRE[{wi}] edges={}\n", we.len()));
            for e in we {
                let edge = Edge(e.0.clone());
                let (e3a, e3b) = BRepTool::edge_parameters(&edge);
                let c3 = BRepTool::edge_curve(&edge);
                let mut line = format!(
                    "  EDGE[{ei}] {} c3={} c3r=[{e3a:.6},{e3b:.6}]",
                    if edge.0.orientation().is_reversed() {
                        "REV"
                    } else {
                        "FWD"
                    },
                    c3.as_ref().map_or("null".into(), |c| c3_kind(c.as_ref())),
                );
                match occt_topo::boptools_2d::curve_on_surface_range(&edge, &face) {
                    Some((pc, f0, f1)) => {
                        let p1 = pc.d0(f0);
                        let p2 = pc.d0(f1);
                        let pm = pc.d0(0.5 * (f0 + f1));
                        line.push_str(&format!(
                            " pc={} pcr=[{f0:.6},{f1:.6}] p1=({:.17e},{:.17e}) p2=({:.17e},{:.17e}) pm=({:.17e},{:.17e})",
                            c2_kind(pc.as_ref()),
                            p1.x(),
                            p1.y(),
                            p2.x(),
                            p2.y(),
                            pm.x(),
                            pm.y()
                        ));
                        if let Some((poles, _)) = pc.bspline_poles2d() {
                            line.push_str(&format!(
                                " deg={:?} npol={} nknot={:?}",
                                pc.bspline_degree(),
                                poles.len(),
                                pc.bspline_nb_knots()
                            ));
                        }
                        if let Some(l) = pc.gp_lin2d() {
                            let d = l.direction();
                            line.push_str(&format!(" dir=({:.6},{:.6})", d.x(), d.y()));
                        }
                    }
                    None => line.push_str(" pc=null"),
                }
                // Replay `ShapeFix_Edge::FixAddPCurve`'s projection call with a
                // fresh projector (empty cache) - the same replay the oracle
                // probe's `REPROJ` block does, so the two verdicts can be
                // compared per edge.
                if let Some(c3) = c3.as_ref() {
                    let prec = BRepTool::edge_tolerance(&edge);
                    let (v1, v2) = occt_topo::topo_tools::edge_vertices(&edge);
                    let tf = v1.as_ref().map_or(0.0, BRepTool::vertex_tolerance);
                    let tl = v2.as_ref().map_or(0.0, BRepTool::vertex_tolerance);
                    let mut cache = occt_topo::pcurve_full::ProjectorCache::new();
                    let rep = occt_topo::pcurve_full::project_curve_on_surface_perform(
                        c3.as_ref(),
                        surf.as_ref(),
                        e3a,
                        e3b,
                        prec,
                        tf,
                        tl,
                        &mut cache,
                    );
                    match rep {
                        Some(r) => {
                            let ra = r.first_parameter();
                            let rb = r.last_parameter();
                            let q1 = r.d0(ra.min(rb));
                            let q2 = r.d0(ra.max(rb));
                            line.push_str(&format!(
                                " reproj={} rrng=[{ra:.6},{rb:.6}] rq1=({:.6},{:.6}) rq2=({:.6},{:.6})",
                                c2_kind(r.as_ref()),
                                q1.x(),
                                q1.y(),
                                q2.x(),
                                q2.y()
                            ));
                            if let Some((poles, _)) = r.bspline_poles2d() {
                                line.push_str(&format!(" rnpol={}", poles.len()));
                            }
                        }
                        None => line.push_str(" reproj=null"),
                    }
                }
                // `PCDUMP_FULL=1`: the same facts at 17 significant digits, so a
                // sub-ulp pcurve residue (the OCCT pole-edge v offset seen in
                // `--wiredbg`) can be compared. Additive; the default line above
                // keeps the recorded `_port_sph_pc.txt` format.
                if std::env::var("PCDUMP_FULL").is_ok() {
                    if let Some((pc, f0, f1)) =
                        occt_topo::boptools_2d::curve_on_surface_range(&edge, &face)
                    {
                        let q1 = pc.d0(f0);
                        let q2 = pc.d0(f1);
                        let mut fl = format!(
                            "FULL EDGE[{ei}] kind={} sameParam={} sameRange={} deg={} rng=[{f0:.17e},{f1:.17e}] q1=({:.17e},{:.17e}) q2=({:.17e},{:.17e})",
                            c2_kind(pc.as_ref()),
                            BRepTool::same_parameter(&edge) as u8,
                            occt_topo::tgeometry::GeometryRegistry::global()
                                .same_range(&edge.0) as u8,
                            BRepTool::is_degenerated(&edge) as u8,
                            q1.x(),
                            q1.y(),
                            q2.x(),
                            q2.y()
                        );
                        if let Some(l) = pc.gp_lin2d() {
                            let d = l.direction();
                            let o = l.location();
                            fl.push_str(&format!(
                                " loc=({:.17e},{:.17e}) dir=({:.17e},{:.17e})",
                                o.x(),
                                o.y(),
                                d.x(),
                                d.y()
                            ));
                        }
                        out.push_str(&fl);
                        out.push('\n');
                    }
                }
                out.push_str(&line);
                out.push('\n');
                ei += 1;
            }
        }
    }
    match args.get(1) {
        Some(p) => {
            let mut fp = std::fs::File::create(p).expect("create dump");
            let _ = fp.write_all(out.as_bytes());
        }
        None => print!("{out}"),
    }
    let _ = Arc::new(());
}
