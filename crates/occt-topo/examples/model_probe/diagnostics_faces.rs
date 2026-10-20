//! Face-level diagnostics (--faccori, --adj, --uvf, ...).
use super::tname;
use super::Ctx;
use occt_topo::brep_bnd_lib;
use occt_topo::brep_tool::BRepTool;
use occt_topo::meshing::range_splitter::mesh_surface_type;
use occt_topo::tgeometry::GeometryRegistry;
use occt_topo::topo_tools::{edge_vertices, edges_of, edges_of_wire, faces_of, wires_of_face};

pub(crate) fn run(ctx: &Ctx) -> bool {
    let model = &ctx.model;
    let shape = ctx.shape.clone();
    // TEMP: per-face occurrence orientation + bbox, to pair the OBJ winding
    // clusters with the faces whose surface normal opposes the solid.
    if std::env::args().any(|a| a == "--faccori") {
        for (i, f) in faces_of(&shape).iter().enumerate() {
            let bb = brep_bnd_lib::shape_bnd_box(&f.0);
            let (b0, b1) = (bb.corner_min(), bb.corner_max());
            println!(
                "FACCORI face={} ori={:?} type={} bbox=({:.6},{:.6},{:.6})-({:.6},{:.6},{:.6})",
                i,
                f.0.orientation(),
                BRepTool::face_surface(f)
                    .map(|s| tname(mesh_surface_type(s.as_ref())))
                    .unwrap_or("none"),
                b0.x(),
                b0.y(),
                b0.z(),
                b1.x(),
                b1.y(),
                b1.z()
            );
        }
        return true;
    }
    if let Some(pos) = std::env::args().position(|a| a == "--adj") {
        // Which faces of the model share the edges of one face, and whether the
        // pcurve lookup finds a representation on each of them (compared with
        // OCCT `occt_probe --edev <face>`).
        let want: usize = std::env::args()
            .nth(pos + 1)
            .and_then(|s| s.parse().ok())
            .expect("--adj <faceIdx>");
        let all = faces_of(&shape);
        let f = all[want].clone();
        let reg = GeometryRegistry::global();
        for (wi, w) in wires_of_face(&f).iter().enumerate() {
            for (ei, e) in edges_of_wire(w).iter().enumerate() {
                let reps = reg.edge_pcurve_reps(&e.0);
                let kinds: Vec<&str> = reps
                    .iter()
                    .map(|(s, _, _, _)| tname(mesh_surface_type(s.as_ref())))
                    .collect();
                print!("ADJ w{wi}e{ei} nreps={} rep_types={kinds:?}", reps.len());
                for (fi, ff) in all.iter().enumerate() {
                    let shares = edges_of(&ff.0)
                        .iter()
                        .any(|q| std::sync::Arc::as_ptr(&q.0.tshape) == std::sync::Arc::as_ptr(&e.0.tshape));
                    if !shares {
                        continue;
                    }
                    let key = GeometryRegistry::shape_key(&ff.0);
                    let n = reg.edge_pcurves(&e.0, key).len();
                    let t = BRepTool::face_surface(ff)
                        .map(|s| tname(mesh_surface_type(s.as_ref())))
                        .unwrap_or("none");
                    print!(" | f{fi}:{t}:key={key}:pc={n}");
                }
                println!();
            }
        }
        return true;
    }
    if std::env::args().any(|a| a == "--pcmismatch") {
        let reg = GeometryRegistry::global();
        let mut total_faces = 0usize;
        let mut bad_faces = 0usize;
        for (ni, s) in model.shapes.iter().enumerate() {
            let mut nb = 0usize;
            let mut bf = 0usize;
            for f in faces_of(&s.shape) {
                nb += 1;
                let key = GeometryRegistry::shape_key(&f.0);
                let mut bad = false;
                for w in wires_of_face(&f) {
                    for e in edges_of_wire(&w) {
                        if !reg.edge_pcurve_reps(&e.0).is_empty()
                            && reg.edge_pcurves(&e.0, key).is_empty()
                        {
                            bad = true;
                        }
                    }
                }
                if bad {
                    bf += 1;
                }
            }
            total_faces += nb;
            bad_faces += bf;
            if bf > 0 {
                println!("PCMISMATCH n={ni} name={} bad={bf}/{nb}", s.name);
            }
        }
        println!("PCMISMATCH total bad={bad_faces}/{total_faces}");
        return true;
    }
    if let Some(pos) = std::env::args().position(|a| a == "--pcdump") {
        // Both stored 2D representations of every edge on one face, i.e. the
        // port side of `BRep_Tool::CurveOnClosedSurface(E, F, c1, c2, f, l)`.
        let want: usize = std::env::args()
            .nth(pos + 1)
            .and_then(|s| s.parse().ok())
            .expect("--pcdump <faceIdx>");
        let f = faces_of(&shape)[want].clone();
        let face_key = GeometryRegistry::shape_key(&f.0);
        println!("PCD2 face={want} key={face_key}");
        let reg = GeometryRegistry::global();
        for (wi, w) in wires_of_face(&f).iter().enumerate() {
            for (ei, e) in edges_of_wire(w).iter().enumerate() {
                let pcs = reg.edge_pcurves(&e.0, face_key);
                println!("PCD2  w{wi}e{ei} ori={:?} n={}", e.0.orientation(), pcs.len());
                for (q, c) in pcs.iter().enumerate() {
                    println!(
                        "PCD2   q={q} kind={:?} f={:.17e} l={:.17e}",
                        occt_geom2d::adaptor_curve::type2d(c.as_ref()),
                        c.first_parameter(),
                        c.last_parameter()
                    );
                    if let Some(b) = c.bspline_copy2d() {
                        println!(
                            "PCD2   q={q} deg={} np={} nk={}",
                            b.degree(),
                            b.xs.len(),
                            b.knots.len()
                        );
                        for i in 0..b.xs.len() {
                            println!("PCD2   q={q} pole[{}]={:.17e},{:.17e}", i + 1, b.xs[i], b.ys[i]);
                        }
                    }
                }
            }
        }
        return true;
    }
    if let Some(pos) = std::env::args().position(|a| a == "--uvdbg") {
        // Per-edge UV contribution of one face: which edge widens the
        // `BRepTools::UVBounds` window, and with what pcurve range/endpoints.
        // Compare with `occt_probe --typefaces <type>`'s `TF   eN ... pcpar=`.
        let want: usize = std::env::args()
            .nth(pos + 1)
            .and_then(|s| s.parse().ok())
            .expect("--uvdbg <faceIdx>");
        let doc = occt_topo::brep_exchange::document_shape(&model).expect("document shape");
        let f = faces_of(&doc)[want].clone();
        let (u0, u1, v0, v1) = occt_topo::brep_tools::uv_bounds(&f);
        println!("UVDBG face={want} u=[{u0:.6},{u1:.6}] v=[{v0:.6},{v1:.6}] ftol={:.12}",
            BRepTool::face_tolerance(&f));
        for (wi, w) in wires_of_face(&f).iter().enumerate() {
            for (ei, e) in edges_of_wire(w).iter().enumerate() {
                let r = occt_topo::boptools_2d::curve_on_surface_range(e, &f);
                match r {
                    Some((c2d, t1, t2)) => {
                        let p1 = c2d.value(t1);
                        let p2 = c2d.value(t2);
                        let (ev1, ev2) = edge_vertices(e);
                        println!(
                            "UVDBG  w{wi}e{ei} etol={:.12} v1tol={:.12} v2tol={:.12}",
                            BRepTool::edge_tolerance(e),
                            ev1.as_ref().map(BRepTool::vertex_tolerance).unwrap_or(-1.0),
                            ev2.as_ref().map(BRepTool::vertex_tolerance).unwrap_or(-1.0)
                        );
                        // Sample the pcurve between its ends so a bulge beyond
                        // the endpoints is visible without the private
                        // `BndLib_Add2dCurve` port.
                        let mut su = (f64::INFINITY, f64::NEG_INFINITY);
                        let mut sv = (f64::INFINITY, f64::NEG_INFINITY);
                        for k in 0..=64 {
                            let t = t1 + (t2 - t1) * (k as f64) / 64.0;
                            let p = c2d.value(t);
                            su.0 = su.0.min(p.x());
                            su.1 = su.1.max(p.x());
                            sv.0 = sv.0.min(p.y());
                            sv.1 = sv.1.max(p.y());
                        }
                        println!(
                            "UVDBG  w{wi}e{ei} ori={:?} deg={} kind={:?} stored_reps={} range=({t1:.9},{t2:.9}) first=({:.9},{:.9}) last=({:.9},{:.9}) samp_u=[{:.9},{:.9}] samp_v=[{:.9},{:.9}]",
                            e.0.orientation(),
                            BRepTool::is_degenerated(e),
                            occt_geom2d::adaptor_curve::type2d(c2d.as_ref()),
                            GeometryRegistry::global().edge_pcurve_reps(e).len(),
                            p1.x(), p1.y(), p2.x(), p2.y(),
                            su.0, su.1, sv.0, sv.1
                        );
                        if let Some(b) = c2d.bspline_copy2d() {
                            let (mut px0, mut px1) = (f64::INFINITY, f64::NEG_INFINITY);
                            let (mut py0, mut py1) = (f64::INFINITY, f64::NEG_INFINITY);
                            for i in 0..b.nb_poles() {
                                px0 = px0.min(b.xs[i]);
                                px1 = px1.max(b.xs[i]);
                                py0 = py0.min(b.ys[i]);
                                py1 = py1.max(b.ys[i]);
                            }
                            println!(
                                "UVDBG      poles={} pole_u=[{px0:.9},{px1:.9}] pole_v=[{py0:.9},{py1:.9}] curve_range=({:.9},{:.9})",
                                b.nb_poles(),
                                b.first_parameter(),
                                b.last_parameter()
                            );
                        }
                        // OCCT's `BRep_Tool::CurveOnPlane` ends in a concrete
                        // `Geom2d_BSplineCurve`, so its UV window is the pole
                        // hull. Project this edge's 3D poles into the plane to
                        // show the window that hull would give here.
                        if let (Some(surf), Some(c3)) = (
                            BRepTool::face_surface(&f),
                            GeometryRegistry::global().edge_curve(&e.0),
                        ) {
                            if let (Some(pln), Some(poles)) = (surf.gp_pln(), c3.bspline_poles()) {
                                let (mut qx0, mut qx1) = (f64::INFINITY, f64::NEG_INFINITY);
                                let (mut qy0, mut qy1) = (f64::INFINITY, f64::NEG_INFINITY);
                                for p in poles {
                                    let uv = occt_geom::projlib::eval_pln_pnt2d(&pln, p);
                                    qx0 = qx0.min(uv.x());
                                    qx1 = qx1.max(uv.x());
                                    qy0 = qy0.min(uv.y());
                                    qy1 = qy1.max(uv.y());
                                }
                                println!(
                                    "UVDBG      c3_poles={} proj_pole_u=[{qx0:.9},{qx1:.9}] proj_pole_v=[{qy0:.9},{qy1:.9}]",
                                    poles.len()
                                );
                            }
                        }
                    }
                    None => println!(
                        "UVDBG  w{wi}e{ei} ori={:?} deg={} no-pcurve",
                        e.0.orientation(),
                        BRepTool::is_degenerated(e)
                    ),
                }
            }
        }
        return true;
    }
    if std::env::args().any(|a| a == "--uvfaces") {
        // Per-face UV window: `BRepAdaptor_Surface(F)` in OCCT restricts the
        // surface range to `BRepTools::UVBounds(F)` (`BRepAdaptor_Surface.cxx:56-81`),
        // so this is the pcurve-shrunk domain compared with `occt_probe --typefaces`.
        // Use the first root, like the OCCT probe's `OneShape()`; `document_shape`
        // compounds every root and would traverse a wrapper compound twice.
        let doc = shape.clone();
        for (i, f) in faces_of(&doc).iter().enumerate() {
            let surf = BRepTool::face_surface(f);
            let st = surf
                .as_ref()
                .map(|s| tname(mesh_surface_type(s.as_ref())))
                .unwrap_or("none");
            let (u0, u1, v0, v1) = occt_topo::brep_tools::uv_bounds(f);
            let (su, sv) = surf
                .as_ref()
                .map(|s| (s.u_range(), s.v_range()))
                .unwrap_or(((0.0, 0.0), (0.0, 0.0)));
            let bb = brep_bnd_lib::shape_bnd_box(&f.0);
            let (b0, b1) = (bb.corner_min(), bb.corner_max());
            let trim = surf
                .as_ref()
                .and_then(|s| s.rectangular_trimmed_spec());
            println!(
                "UVF face={i} type={st} u=[{u0:.6},{u1:.6}] v=[{v0:.6},{v1:.6}] nat_u={su:?} nat_v={sv:?} trim={trim:?} bbox=({:.6},{:.6},{:.6})-({:.6},{:.6},{:.6})",
                b0.x(), b0.y(), b0.z(), b1.x(), b1.y(), b1.z()
            );
        }
        return true;
    }
    false
}
