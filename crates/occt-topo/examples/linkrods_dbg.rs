use occt_core::gp::GpPnt2d;
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_topo::brep_tools;
use occt_topo::boptools_2d::curve_on_surface;
use occt_topo::brep_exchange::brep_to_obj;
use occt_topo::brep_tool::BRepTool;
use occt_topo::meshing::data_model::MeshFace;
use occt_topo::meshing::face_discret::Classifier;
use occt_topo::meshing::parameters::MeshParameters;
use occt_topo::meshing::incremental_mesh::IncrementalMesh;
use occt_topo::meshing::range_splitter::{
    classify_surface, create_range_splitter, RangeSplitter, SurfaceType,
};
use occt_topo::shhealing::{edge_proj_aux_params, fix_shifted_wire};
use occt_topo::step::read_step_file;
use occt_topo::abs::Orientation;
use occt_topo::topo_tools_full::{edge_vertices, edges_of, edges_of_wire, faces_of, wires_of_face};

fn vf(obj: &str) -> (usize, usize) {
    let mut v = 0usize;
    let mut f = 0usize;
    for line in obj.lines() {
        if line.starts_with("v ") {
            v += 1;
        } else if line.starts_with("f ") {
            f += 1;
        }
    }
    (v, f)
}

fn type_name(t: SurfaceType) -> &'static str {
    match t {
        SurfaceType::Plane => "Plane",
        SurfaceType::Cylinder => "Cylinder",
        SurfaceType::Cone => "Cone",
        SurfaceType::Sphere => "Sphere",
        SurfaceType::Torus => "Torus",
        SurfaceType::SurfaceOfRevolution => "Rev",
        SurfaceType::SurfaceOfExtrusion => "Extr",
        SurfaceType::BezierSurface => "Bezier",
        SurfaceType::BSplineSurface => "BSpline",
        SurfaceType::OffsetSurface => "Offset",
        SurfaceType::OtherSurface => "Other",
    }
}

fn prs3d_lin(shape: &occt_topo::shape::TopoShape, maximal_chordial: f64) -> f64 {
    occt_topo::brep_exchange::prs3d_get_deflection(shape, maximal_chordial)
}

fn dump_grid(shape: &occt_topo::shape::TopoShape, stem: &str) {
    let lin = prs3d_lin(shape, 0.1);
    let mut params = MeshParameters {
        deflection: lin,
        angle: 20.0 * std::f64::consts::PI / 180.0,
        ..MeshParameters::default()
    };
    if params.deflection_interior < 1e-7 {
        params.deflection_interior = params.deflection;
    }
    if params.min_size < 1e-7 {
        params.min_size = (0.1 * params.deflection.min(params.deflection_interior)).max(1e-7);
    }
    if params.angle_interior < 1e-12 {
        params.angle_interior = 2.0 * params.angle;
    }
    println!("{stem} lin={lin:.6} min_size={:.6}", params.min_size);
    let mut none_n = 0usize;
    let mut empty_n = 0usize;
    let mut sum_nodes = 0usize;
    let mut face_i = 0usize;
    for face in faces_of(shape) {
        face_i += 1;
        let Some(surf) = BRepTool::face_surface(&face) else {
            println!("{stem} f{face_i} no-surf");
            continue;
        };
        let ty = classify_surface(surf.as_ref());
        let (su0, su1, sv0, sv1) = BRepTool::uv_bounds(&face);
        let (u0, u1, v0, v1) = brep_tools::uv_bounds(&face);
        let nu_p = surf.nb_u_poles();
        let nv_p = surf.nb_v_poles();
        let (ru, rv) = surf
            .uv_resolution(lin)
            .unwrap_or((f64::NAN, f64::NAN));
        let adapt_nb = |is_u: bool, range: (f64, f64)| -> i32 {
            let iso = if is_u {
                let vfirst = surf.v_intervals(0).first().copied().unwrap_or(sv0);
                surf.v_iso_curve(vfirst)
            } else {
                let ufirst = surf.u_intervals(0).first().copied().unwrap_or(su0);
                surf.u_iso_curve(ufirst)
            };
            match iso {
                Some(c) => {
                    if let (Some(knots), Some(degree)) = (c.bspline_knots(), c.nurbs_degree()) {
                        let iv = occt_core::bspl::adaptor_intervals(
                            knots,
                            degree,
                            c.is_periodic(),
                            6,
                            range.0,
                            range.1,
                            c.resolution(CONFUSION).min(PCONFUSION),
                        );
                        iv.len().saturating_sub(1) as i32
                    } else {
                        1
                    }
                }
                None => {
                    if is_u {
                        surf.nb_u_intervals(6)
                    } else {
                        surf.nb_v_intervals(6)
                    }
                }
            }
        };
        let au = adapt_nb(true, (u0, u1));
        let av = adapt_nb(false, (v0, v1));
        let fb_u = au == 1 && nu_p > 2;
        let fb_v = av == 1 && nv_p > 2;
        let mut mf = MeshFace::new(face.clone());
        mf.set_deflection(lin);
        let mut sp = create_range_splitter(surf.as_ref());
        sp.reset(&mf, &params);
        let mut wires_uv: Vec<Vec<GpPnt2d>> = Vec::new();
        for w in wires_of_face(&face) {
            let mut poly = Vec::new();
            for e in edges_of_wire(&w) {
                let Some(pc) = curve_on_surface(&e, &face) else {
                    continue;
                };
                let (fp, lp) = (pc.first_parameter(), pc.last_parameter());
                if !fp.is_finite() || !lp.is_finite() {
                    continue;
                }
                for k in 0..=8 {
                    let t = fp + (lp - fp) * k as f64 / 8.0;
                    let p = pc.d0(t);
                    if p.x().is_finite() && p.y().is_finite() {
                        sp.add_point(p);
                        poly.push(p);
                    }
                }
            }
            if poly.len() >= 3 {
                wires_uv.push(poly);
            }
        }
        if !u0.is_finite() {
            sp.add_point(GpPnt2d::new(0.0, 0.0));
            sp.add_point(GpPnt2d::new(1.0, 1.0));
        }
        sp.adjust_range();
        let nodes = sp.generate_surface_nodes(&params);
        let (nlen, tag) = match &nodes {
            None => {
                none_n += 1;
                (0usize, "none")
            }
            Some(v) if v.is_empty() => {
                empty_n += 1;
                (0usize, "empty")
            }
            Some(v) => {
                sum_nodes += v.len();
                (v.len(), "ok")
            }
        };
        let mut inside = 0usize;
        if let Some(v) = &nodes {
            if !wires_uv.is_empty() && u0.is_finite() {
                let mut clf = Classifier::new();
                let tol = sp.tolerance_uv();
                for w in &wires_uv {
                    clf.register_wire(w, tol, sp.range_u(), sp.range_v());
                }
                inside = v.iter().filter(|p| clf.is_inside(p)).count();
            }
        }
        let mut u_n = 0usize;
        let mut v_n = 0usize;
        let mut max_sq = f64::NAN;
        if let Some(v) = &nodes {
            let mut us: Vec<f64> = v.iter().map(|p| p.x()).collect();
            let mut vs: Vec<f64> = v.iter().map(|p| p.y()).collect();
            us.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            vs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            us.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
            vs.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
            u_n = us.len();
            v_n = vs.len();
            if ty == SurfaceType::BSplineSurface && vs.len() >= 2 {
                let um = 0.5 * (u0 + u1);
                if let Some(iso) = surf.u_iso_curve(um) {
                    max_sq = 0.0;
                    for w in vs.windows(2) {
                        let (a, b) = (w[0], w[1]);
                        let pa = iso.d0(a);
                        let pb = iso.d0(b);
                        let pm = iso.d0(0.5 * (a + b));
                        let sq = occt_topo::meshing::geom_tool::GeomTool::square_deflection_of_segment(
                            &pa, &pb, &pm,
                        );
                        if sq > max_sq {
                            max_sq = sq;
                        }
                    }
                }
            }
        }
        let mut u_max_sq = f64::NAN;
        let mut nrm_ang = f64::NAN;
        if ty == SurfaceType::BSplineSurface && u0.is_finite() && u1.is_finite() && v0.is_finite() {
            let vm = 0.5 * (v0 + v1);
            let um = 0.5 * (u0 + u1);
            if let Some(iso) = surf.v_iso_curve(vm) {
                let pa = iso.d0(u0);
                let pb = iso.d0(u1);
                let pm = iso.d0(um);
                u_max_sq = occt_topo::meshing::geom_tool::GeomTool::square_deflection_of_segment(
                    &pa, &pb, &pm,
                );
                let pmid = iso.d0(um);
                let pql = iso.d0(0.5 * (u0 + um));
                let pqr = iso.d0(0.5 * (um + u1));
                let left = occt_topo::meshing::geom_tool::GeomTool::square_deflection_of_segment(
                    &pa, &pmid, &pql,
                );
                let right = occt_topo::meshing::geom_tool::GeomTool::square_deflection_of_segment(
                    &pmid, &pb, &pqr,
                );
                u_max_sq = u_max_sq.max(left).max(right);
            }
            let (s1, n1) = occt_topo::meshing::geomlib_norm::norm_estim(
                surf.as_ref(),
                u0,
                vm,
                CONFUSION,
            );
            let (s2, n2) = occt_topo::meshing::geomlib_norm::norm_estim(
                surf.as_ref(),
                um,
                vm,
                CONFUSION,
            );
            if s1 == 0 && s2 == 0 {
                if let (Some(a), Some(b)) = (n1, n2) {
                    nrm_ang = a.angle(&b);
                }
            }
        }
        println!(
            "{stem} f{face_i} {} poles={nu_p}/{nv_p} adapt={au}/{av} fb={fb_u}/{fb_v} faceUV=[{u0:.4},{u1:.4}]x[{v0:.4},{v1:.4}] grid={nlen} uv={u_n}/{v_n} in={inside} maxsq={max_sq:.6} umaxsq={u_max_sq:.6} nrmang={nrm_ang:.4} {tag}",
            type_name(ty)
        );
    }
    println!("{stem} faces={face_i} grid_sum={sum_nodes} none={none_n} empty={empty_n}");
}

fn dump_nattr(shape: &occt_topo::shape::TopoShape, stem: &str) {
    let mut n = 0usize;
    let mut n_big = 0usize;
    for e in edges_of(shape) {
        let Some(curve) = BRepTool::edge_curve(&e) else {
            continue;
        };
        let (cf, cl) = (curve.first_parameter(), curve.last_parameter());
        if !cf.is_finite() || !cl.is_finite() {
            continue;
        }
        let (v1, v2) = edge_vertices(&e);
        let Some(v1) = v1 else { continue };
        let Some(v2) = v2 else { continue };
        let p1 = BRepTool::vertex_point(&v1);
        let p2 = BRepTool::vertex_point(&v2);
        let pf = curve.d0(cf);
        let pl = curve.d0(cl);
        let nat = pf.distance(&p1).max(pl.distance(&p2));
        let swp = pf.distance(&p2).max(pl.distance(&p1));
        let (ef, el) = BRepTool::edge_parameters(&e);
        n += 1;
        if nat > 1e-3 && swp > 1e-3 {
            n_big += 1;
            let kind = if curve.is_line() {
                "line"
            } else if curve.gp_circ().is_some() {
                "circ"
            } else if curve.bspline_poles().is_some() {
                "bspl"
            } else {
                "other"
            };
            println!(
                "{stem} e nat={nat:.4} swp={swp:.4} {kind} curve=[{cf:.4},{cl:.4}] edge=[{ef:.4},{el:.4}]"
            );
        }
    }
    println!("{stem} edges={n} mismatch={n_big}");
}

fn dump_faces(shape: &occt_topo::shape::TopoShape, stem: &str) {
    let lin = prs3d_lin(shape, 0.1);
    const PRS3D_DEV_ANGLE: f64 = 20.0 * std::f64::consts::PI / 180.0;
    let inc = IncrementalMesh::from_deflection(shape, lin, false, PRS3D_DEV_ANGLE);
    let faces = faces_of(shape);
    let mut tot_v = [0usize; 11];
    let mut tot_t = [0usize; 11];
    let mut tot_n = [0usize; 11];
    println!("{stem} lin={lin:.6} faces={}", inc.face_stats().len());
    let mut seen = std::collections::HashSet::new();
    for st in inc.face_stats() {
        seen.insert(st.index);
        let ty = type_name(st.surface);
        let bucket = match st.surface {
            SurfaceType::Plane => 0,
            SurfaceType::Cylinder => 1,
            SurfaceType::Cone => 2,
            SurfaceType::Sphere => 3,
            SurfaceType::Torus => 4,
            SurfaceType::SurfaceOfRevolution => 5,
            SurfaceType::SurfaceOfExtrusion => 6,
            SurfaceType::BezierSurface => 7,
            SurfaceType::BSplineSurface => 8,
            SurfaceType::OffsetSurface => 9,
            SurfaceType::OtherSurface => 10,
        };
        tot_v[bucket] += st.vertices;
        tot_t[bucket] += st.triangles;
        tot_n[bucket] += 1;
        let extra = faces.get(st.index).and_then(|face| {
            let surf = BRepTool::face_surface(face)?;
            let (u0, u1, v0, v1) = BRepTool::uv_bounds(face);
            Some(format!(
                " uv=[{:.4},{:.4}]x[{:.4},{:.4}] poles={}x{}",
                u0,
                u1,
                v0,
                v1,
                surf.nb_u_poles(),
                surf.nb_v_poles()
            ))
        });
        println!(
            "{stem} f{} {} {}/{}{}",
            st.index,
            ty,
            st.vertices,
            st.triangles,
            extra.unwrap_or_default()
        );
    }
    for (i, face) in faces.iter().enumerate() {
        if seen.contains(&i) {
            continue;
        }
        let ty = BRepTool::face_surface(face)
            .map(|s| type_name(classify_surface(s.as_ref())))
            .unwrap_or("no-surf");
        let (u0, u1, v0, v1) = BRepTool::uv_bounds(face);
        let poles = BRepTool::face_surface(face)
            .map(|s| format!("{}x{}", s.nb_u_poles(), s.nb_v_poles()))
            .unwrap_or_else(|| "-".to_string());
        println!(
            "{stem} f{i} {ty} SKIPPED uv=[{u0:.4},{u1:.4}]x[{v0:.4},{v1:.4}] poles={poles}"
        );
    }
    for (i, name) in [
        "Plane",
        "Cylinder",
        "Cone",
        "Sphere",
        "Torus",
        "Rev",
        "Extr",
        "Bezier",
        "BSpline",
        "Offset",
        "Other",
    ]
    .iter()
    .enumerate()
    {
        if tot_n[i] > 0 {
            println!(
                "{stem} sum {name} n={} {}/{}",
                tot_n[i], tot_v[i], tot_t[i]
            );
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let apply_shift = args.iter().any(|a| a == "--shift");
    let ranges_only = args.iter().any(|a| a == "--ranges");
    let grid_only = args.iter().any(|a| a == "--grid");
    let nattr_only = args.iter().any(|a| a == "--nattr");
    let faces_only = args.iter().any(|a| a == "--faces");
    let only = args.into_iter().find(|a| {
        a != "--shift"
            && a != "--ranges"
            && a != "--grid"
            && a != "--nattr"
            && a != "--faces"
            && a != "--sph"
    });
    let stems = [
        "OffsetPlaneHoleEdge",
        "HoledPlate",
        "Cube",
        "Cone",
        "Cylinder",
        "Torus",
        "Sphere",
        "Shape",
        "Shape-1",
        "Shape-2",
        "linkrods",
        "screw",
        "rev",
        "Offset",
        "Extrusion",
    ];
    for stem in stems {
        if let Some(ref s) = only {
            if stem != s.as_str() {
                continue;
            }
        }
        let path = format!(r"D:\source\repos\dogs\data\{stem}.step");
        match read_step_file(&path) {
            Ok(model) => {
                let shape = &model.shapes[0].shape;
                if apply_shift {
                    for face in faces_of(shape) {
                        for w in wires_of_face(&face) {
                            let _ = fix_shifted_wire(&w, &face);
                        }
                    }
                }
                if nattr_only {
                    dump_nattr(shape, stem);
                    continue;
                }
                if grid_only {
                    dump_grid(shape, stem);
                    continue;
                }
                if faces_only {
                    dump_faces(shape, stem);
                    continue;
                }
                if ranges_only {
                    for face in faces_of(shape) {
                        for w in wires_of_face(&face) {
                            for e in edges_of_wire(&w) {
                                let (a, b) = BRepTool::edge_parameters(&e);
                                let Some(pc) = curve_on_surface(&e, &face) else {
                                    continue;
                                };
                                let (fp, lp) = (pc.first_parameter(), pc.last_parameter());
                                let (v1, v2) = edge_vertices(&e);
                                let vd = match (v1, v2) {
                                    (Some(x), Some(y)) => BRepTool::vertex_point(&x)
                                        .distance(&BRepTool::vertex_point(&y)),
                                    _ => -1.0,
                                };
                                let (p1, p2) = BRepTool::edge_curve(&e)
                                    .map(|c| (c.d0(a), c.d0(b)))
                                    .unwrap_or_else(|| (occt_core::gp::GpPnt::zero(), occt_core::gp::GpPnt::zero()));
                                let ends = p1.distance(&p2);
                                let proj = edge_proj_aux_params(&e, &face, 1e-3);
                                println!(
                                    "{stem} 3d=[{a:.6},{b:.6}] 2d=[{fp:.6},{lp:.6}] d3={:.6} d2={:.6} vd={vd:.6} ends={ends:.6} line2d={} proj={proj:?} per={}",
                                    b - a,
                                    lp - fp,
                                    !pc.first_parameter().is_finite(),
                                    pc.is_periodic()
                                );
                            }
                        }
                    }
                    println!("{stem} read-ok");
                    continue;
                }
                let obj = brep_to_obj(shape, 0.1);
                let (v, f) = vf(&obj);
                println!("{stem} {v}/{f}");
                let mut internal_on_edge = 0usize;
                for e in edges_of(shape) {
                    let kids = e
                        .0
                        .tshape
                        .read()
                        .expect("poisoned TShape lock")
                        .children
                        .clone();
                    internal_on_edge += kids
                        .iter()
                        .filter(|c| c.is_vertex() && c.orientation() == Orientation::Internal)
                        .count();
                }
                let mut face_internal = 0usize;
                let mut ty_n = [0u32; 12];
                for face in faces_of(shape) {
                    if let Some(surf) = BRepTool::face_surface(&face) {
                        let ty = classify_surface(surf.as_ref());
                        let i = match ty {
                            SurfaceType::Plane => 0,
                            SurfaceType::Cylinder => 1,
                            SurfaceType::Cone => 2,
                            SurfaceType::Sphere => 3,
                            SurfaceType::Torus => 4,
                            SurfaceType::SurfaceOfRevolution => 5,
                            SurfaceType::SurfaceOfExtrusion => 6,
                            SurfaceType::BezierSurface => 7,
                            SurfaceType::BSplineSurface => 8,
                            SurfaceType::OffsetSurface => 9,
                            SurfaceType::OtherSurface => 10,
                        };
                        ty_n[i] += 1;
                    }
                    let kids = face
                        .0
                        .tshape
                        .read()
                        .expect("poisoned TShape lock")
                        .children
                        .clone();
                    face_internal += kids
                        .iter()
                        .filter(|c| c.is_vertex() && c.orientation() == Orientation::Internal)
                        .count();
                }
                println!(
                    "{stem} internal_on_edge={internal_on_edge} face_child_internal_v={face_internal} faces Pln={} Cyl={} Cone={} Sph={} Tor={} Rev={} Extr={} Bez={} BSpl={} Off={} Oth={}",
                    ty_n[0], ty_n[1], ty_n[2], ty_n[3], ty_n[4], ty_n[5], ty_n[6], ty_n[7], ty_n[8], ty_n[9], ty_n[10],
                );
            }
            Err(e) => println!("{stem} ERR {e}"),
        }
    }
}
