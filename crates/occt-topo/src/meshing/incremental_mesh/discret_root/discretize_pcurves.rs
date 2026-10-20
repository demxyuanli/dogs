//! `IncrementalMesh::discretize_pcurves` — pcurve discretization for a model.

use super::*;

impl IncrementalMesh {
    /// Discretize every edge's 3D curve once (the shared 3D polyline, OCCT
    /// `BRepMesh_EdgeDiscret::Tessellate3d`) and mirror its parameters into each
    /// pcurve (`Tessellate2d`), so the 2D and 3D boundary points stay aligned and
    /// adjacent faces meet exactly on a shared edge.
    pub(crate) fn discretize_pcurves(
        &self,
        model: &mut MeshModel,
        min_size: f64,
    ) -> Result<(), String> {
        // Collect the pcurve jobs first (an immutable pass), then populate
        // (a mutable pass) so `model` is never borrowed twice at once.
        let mut jobs: Vec<(usize, usize, Edge, Face, usize, f64, f64)> = Vec::new();
        for i in 0..model.edges_nb() {
            let edge = model.edge(i)?;
            // The pcurve is a property of the (edge, face) pair in the edge's
            // *natural* curve direction, independent of how a wire orients it;
            // the wire orientation is applied later when the boundary UV is
            // collected. The 3D parameter range is computed on the Forward edge
            // (`a -> b`), but each pcurve carries its own wire orientation — a
            // seam edge has a Forward and a Reversed pcurve on the same face,
            // and `make_pcurve_full` selects the matching side from it.
            let fwd_edge = Edge(edge.edge().0.oriented(Orientation::Forward));
            let (a, b) = BRepTool::edge_parameters(&fwd_edge);
            if !(a.is_finite() && b.is_finite() && b - a >= 1e-15) {
                continue;
            }
            for p in 0..edge.pcurves_nb() {
                let face_index = edge.pcurve(p)?.face();
                let orientation = edge.pcurve(p)?.orientation();
                let topo_face = model.face(face_index)?.face().clone();
                let oriented_edge = Edge(edge.edge().0.oriented(orientation));
                jobs.push((i, p, oriented_edge, topo_face, face_index, a, b));
            }
        }

        // Tessellate each edge's 3D curve once (cached), then reuse its parameter
        // sequence for every pcurve of that edge.
        let mut edge_params: HashMap<usize, Vec<f64>> = HashMap::new();
        for (edge_index, pcurve_index, topo_edge, topo_face, pc_face, a, b) in jobs {
            let params = match edge_params.get(&edge_index) {
                Some(p) => p.clone(),
                None => {
                    // `BRepMesh_EdgeDiscret::process` computes the edge's
                    // deflection (`BRepMesh_Deflection::ComputeDeflection`) before
                    // tessellating its 3D curve. `EdgeAmplifier` already lowered
                    // the stored deflection (`BRepMesh_ModelHealer.cxx:56-57`);
                    // do not overwrite it.
                    if !model.edge(edge_index)?.is_status(MeshStatus::REMESH) {
                        ModelPreProcessor::compute_edge_deflection(
                            model,
                            edge_index,
                            &self.parameters,
                        )?;
                    }
                    // `BRepMesh_EdgeDiscret::process` (`cxx:113-126`):
                    // CheckAndUpdateFlags on every pcurve before tessellation.
                    let n_pc_flags = model.edge(edge_index)?.pcurves_nb();
                    for p in 0..n_pc_flags {
                        let face_index = model.edge(edge_index)?.pcurve(p)?.face();
                        let topo_face = model.face(face_index)?.face().clone();
                        ShapeTool::check_and_update_flags(
                            model.edge_mut(edge_index)?,
                            &topo_face,
                        );
                    }
                    let curve = model.edge(edge_index)?.curve().ok_or_else(|| {
                        format!("IncrementalMesh::discretize_pcurves: edge {edge_index} has no 3D curve")
                    })?;
                    let (mut def, ang, degenerated, is_free, same_param, same_range, ori) = {
                        let e = model.edge(edge_index)?;
                        // `BRepMesh_CurveTessellator::init`: both linear and angular
                        // deflections are halved (`BRepMesh_CurveTessellator.cxx:73-74`).
                        (
                            (0.5 * e.deflection()).max(1e-9),
                            0.5 * e.angular_deflection(),
                            e.is_degenerated(),
                            e.is_free(),
                            e.same_param(),
                            e.same_range(),
                            e.edge().0.orientation(),
                        )
                    };
                    // INTERNAL edges take an extra half on the linear term (`cxx:75-77`).
                    if ori == Orientation::Internal {
                        def = (0.5 * def).max(1e-9);
                    }
                    // `CreateEdgeTessellator` (`BRepMesh_EdgeDiscret.cxx:57-63`):
                    // !SameParam uses `BRepAdaptor_Curve(edge, face)`.
                    let mut cos_value: Option<(
                        std::sync::Arc<dyn occt_geom2d::curve::Curve2d>,
                        std::sync::Arc<dyn occt_geom::Surface>,
                    )> = None;
                    let (tess_curve, tess_first, tess_last) = if !is_free && !same_param {
                        let (face_index, orientation) = {
                            let pc = model.edge(edge_index)?.pcurve(0)?;
                            (pc.face(), pc.orientation())
                        };
                        let topo_face = model.face(face_index)?.face().clone();
                        let oriented = Edge(
                            model.edge(edge_index)?.edge().0.oriented(orientation),
                        );
                        match (
                            crate::pcurve_full::make_pcurve_full(&oriented, &topo_face),
                            BRepTool::face_surface_world(&topo_face),
                        ) {
                            (Ok(c2), Some(s)) => {
                                // `BRepAdaptor_Curve(E,F)` loads the pcurve on
                                // the COS representation range (`BRep_GCurve`
                                // First/Last after CheckPCurves clamp).
                                let (pf, pl) = crate::boptools_2d::curve_on_surface_range(
                                    &oriented,
                                    &topo_face,
                                )
                                .map(|(_, f, l)| (f, l))
                                .unwrap_or((a, b));
                                cos_value = Some((c2.clone(), s.clone()));
                                (
                                    std::sync::Arc::new(CurveOnSurface::new(c2, s, pf, pl))
                                        as std::sync::Arc<dyn occt_geom::Curve>,
                                    pf,
                                    pl,
                                )
                            }
                            _ => (curve.clone(), a, b),
                        }
                    } else {
                        (curve.clone(), a, b)
                    };
                    let min_pts = tessellator_min_points(tess_curve.as_ref());
                    let tess_edge = model.edge(edge_index)?.edge().clone();
                    let mut tess = CurveTessellator::from_range_angular_min(
                        tess_curve,
                        tess_first,
                        tess_last,
                        def,
                        ang,
                        min_pts,
                        min_size,
                    );
                    tess.add_internal_vertices(&tess_edge);
                    let mut tess_params = tess.params().to_vec();
                    // `splitByDeflection2d` (`BRepMesh_CurveTessellator.cxx:157-188`).
                    if !is_free && same_param && same_range && tess_params.len() > 1 {
                        let mut pcs = Vec::new();
                        let n_pc = model.edge(edge_index)?.pcurves_nb();
                        for p in 0..n_pc {
                            let (face_index, orientation) = {
                                let pc = model.edge(edge_index)?.pcurve(p)?;
                                (pc.face(), pc.orientation())
                            };
                            let topo_face = model.face(face_index)?.face().clone();
                            let oriented = Edge(
                                model.edge(edge_index)?.edge().0.oriented(orientation),
                            );
                            if let (Ok(c2), Some(s)) = (
                                crate::pcurve_full::make_pcurve_full(&oriented, &topo_face),
                                BRepTool::face_surface_world(&topo_face),
                            ) {
                                pcs.push((c2, s));
                            }
                        }
                        split_by_deflection2d(
                            curve.as_ref(),
                            &mut tess_params,
                            &pcs,
                            def,
                            min_size,
                        );
                    }
                    // Store the shared 3D polyline. OCCT `Tessellate3d(theUpdateEnds=true)`
                    // replaces the two end points with the exact vertex points
                    // (`BRep_Tool::Pnt`) of the edge in its natural (Forward)
                    // curve direction — not the pcurve's wire orientation.
                    // Degenerated edges skip interior tessellator points
                    // (`BRepMesh_EdgeDiscret.cxx:250-270`).
                    let natural = Edge(model.edge(edge_index)?.edge().0.oriented(Orientation::Forward));
                    let (first_vertex, last_vertex) =
                        crate::topo_tools::edge_vertices(&natural);
                    let first_pnt = first_vertex.as_ref().map(|v| BRepTool::vertex_point_world(v));
                    let last_pnt = last_vertex.as_ref().map(|v| BRepTool::vertex_point_world(v));
                    let edge_tol = BRepTool::edge_tolerance(&model.edge(edge_index)?.edge());
                    let em = model.edge_mut(edge_index)?;
                    em.discretization_mut().clear(false);
                    let n = tess_params.len();
                    let keep: Vec<usize> = if degenerated && n >= 2 {
                        vec![0, n - 1]
                    } else {
                        let mut idxs: Vec<usize> = (0..n).collect();
                        // `Tessellate3d` (`cxx:256-259`): drop interiors that
                        // `CurveTessellator::Value` rejects.
                        if let Some((c2, s)) = cos_value.as_ref() {
                            idxs.retain(|&k| {
                                if k == 0 || k + 1 == n {
                                    return true;
                                }
                                let t = tess_params[k];
                                let p = curve.d0(t);
                                curve_tessellator_value_ok(
                                    c2.as_ref(),
                                    s.as_ref(),
                                    t,
                                    &p,
                                    edge_tol,
                                )
                            });
                        }
                        idxs
                    };
                    let mut params = Vec::with_capacity(keep.len());
                    for &k in &keep {
                        let t = tess_params[k];
                        let p = if k == 0 {
                            first_pnt.unwrap_or_else(|| curve.d0(t))
                        } else if k == n - 1 {
                            last_pnt.unwrap_or_else(|| curve.d0(t))
                        } else {
                            curve.d0(t)
                        };
                        em.discretization_mut().add_point(p, t);
                        params.push(t);
                    }
                    edge_params.insert(edge_index, params.clone());
                    params
                }
            };

            let pc = crate::pcurve_full::make_pcurve_full(&topo_edge, &topo_face)
                .map_err(|e| format!("IncrementalMesh::discretize_pcurves: {e}"))?;
            // `BRepMesh_EdgeDiscret::Tessellate2d` + `EdgeParameterProvider`: the
            // pcurve parameter is found by projecting the edge's 3D polyline point
            // onto the surface and then onto the pcurve (a non-SameParameter STEP
            // edge's 3D curve may lie off the surface). The 3D polyline was stored
            // above, so read it back for the projection.
            let surface = BRepTool::face_surface_world(&topo_face).ok_or_else(|| {
                "IncrementalMesh::discretize_pcurves: face has no surface".to_string()
            })?;
            let same_param = model.edge(edge_index)?.same_param();
            let (discretization, pcurve) =
                model.edge_discretization_and_pcurve_mut(edge_index, pcurve_index)?;
            let pts3d = discretization.points();
            pcurve.clear(false);
            let _pcurve_face = pcurve.face();
            // The pcurve's parameter range over this edge. A bounded curve is its
            // own `[first, last]`; an unbounded line's range is the arc length
            // between the edge's vertices projected onto the surface then the
            // pcurve (`BRep_Builder::UpdateEdge`).
            // `BRepAdaptor_Curve(E,F).First/Last` is COS range (`hxx:77-84`).
            // After RemovePCurves, `BRep_Tool.cxx:367-372` CurveOnPlane sets
            // First/Last to the 3D `f,l`, not `Geom2dCircle` `[0, 2pi]`.
            let (pf, pl) = if let Some((_, f, l)) =
                crate::boptools_2d::curve_on_surface_range(&topo_edge, &topo_face)
            {
                if f.is_finite() && l.is_finite() {
                    (f, l)
                } else if a.is_finite() && b.is_finite() && b > a {
                    (a, b)
                } else if pc.first_parameter().is_finite() && pc.last_parameter().is_finite() {
                    (pc.first_parameter(), pc.last_parameter())
                } else {
                    (a, b)
                }
            } else if a.is_finite() && b.is_finite() && b > a {
                // Bounded 3D curve (circle / B-spline / trimmed line) whose pcurve
                // is unbounded (a 2D line): the pcurve is parameterized over the
                // edge's own range (`SameParameter`), so use `[a, b]` directly.
                // Projecting the vertices collapses for a *closed* edge whose two
                // vertices coincide (a torus minor seam's inner point), which would
                // shrink the pcurve parameter range to a single value.
                (a, b)
            } else if pc.first_parameter().is_finite() && pc.last_parameter().is_finite() {
                (pc.first_parameter(), pc.last_parameter())
            } else {
                let proj = |v: Option<crate::shape::Vertex>| -> Option<f64> {
                    let p = BRepTool::vertex_point_world(&v?);
                    let pr =
                        occt_geom::geom_api::project_point_on_surface(surface.as_ref(), &p, 1e-7)?;
                    Some(crate::pcurve_full::project_uv_on_curve2d(
                        pc.as_ref(),
                        GpPnt2d::new(pr.u, pr.v),
                    ))
                };
                let (v1, v2) = crate::topo_tools::edge_vertices(&topo_edge);
                match (proj(v1), proj(v2)) {
                    (Some(u1), Some(u2)) => (u1.min(u2), u1.max(u2)),
                    _ => (a, b),
                }
            };
            // TEMP PC2D: the port's own pcurve definition for this edge, in the
            // same shape as the OCCT probe's `OCCTPC2D` line.
            if std::env::var("PC2D").is_ok() {
                let mid = 0.5 * (pf + pl);
                let a0 = pc.d0(pf);
                let am = pc.d0(mid);
                let a1 = pc.d0(pl);
                eprintln!(
                    "PC2D e={edge_index} pc={pcurve_index} type={} f={:.17e} l={:.17e} pf=({:.17e},{:.17e}) pm=({:.17e},{:.17e}) pl=({:.17e},{:.17e})",
                    std::any::type_name_of_val(pc.as_ref()),
                    pf,
                    pl,
                    a0.x(),
                    a0.y(),
                    am.x(),
                    am.y(),
                    a1.x(),
                    a1.y()
                );
            }
            // TEMP T-191 PCDEF: the port's own 2-D curve definition, in the same
            // format as the OCCT probe's `--pcdef` dump, so the two can be diffed
            // pole by pole.
            if std::env::var("PCDEF").is_ok() {
                eprintln!(
                    "PCDEF f={pc_face} e={edge_index} pc={pcurve_index} edge={:p} type={} pf={:.17e} pl={:.17e}",
                    std::sync::Arc::as_ptr(&topo_edge.0.tshape),
                    std::any::type_name_of_val(pc.as_ref()),
                    pc.first_parameter(),
                    pc.last_parameter()
                );
                let trimmed = pc.trimmed_basis().is_some();
                let offs = pc.offset_basis().is_some();
                let knots2d = pc.bspline_knots2d().map(|k| k.len());
                let dist = pc.bspline_distinct_knots_mults().map(|(k, m)| (k.len(), m));
                let copy_bs = pc.bspline_copy2d().is_some();
                eprintln!(
                    "PCDEF  kind trimmed={trimmed} offset={offs} deg={:?} knots2d={knots2d:?} dist={dist:?} copybs={copy_bs}",
                    pc.bspline_degree()
                );
                {
                    let f = pc.first_parameter();
                    let l = pc.last_parameter();
                    let mut s = String::new();
                    for k in 0..=4 {
                        let u = f + (l - f) * (k as f64) / 4.0;
                        let p = pc.d0(u);
                        s.push_str(&format!(" [{k}]{:.17e},{:.17e}", p.x(), p.y()));
                    }
                    let q = pc.d0(-0.19481982426506805);
                    s.push_str(&format!(" [t]{:.17e},{:.17e}", q.x(), q.y()));
                    eprintln!("PCDEF  val{s}");
                }
                if let Some(b) = pc.bspline_copy2d() {
                    eprintln!(
                        "PCDEF  degree={} npoles={} nknots={} periodic={} rational={} flatlen={}",
                        b.degree,
                        b.xs.len(),
                        b.knots.len(),
                        b.periodic as i32,
                        b.weights.is_some() as i32,
                        b.knots.len()
                    );
                    for i in 0..b.xs.len() {
                        eprintln!(
                            "PCDEF  pole[{}]={:.17e},{:.17e}",
                            i + 1,
                            b.xs[i],
                            b.ys[i]
                        );
                    }
                    if let Some((knots, mults)) = pc.bspline_distinct_knots_mults() {
                        for i in 0..knots.len() {
                            eprintln!(
                                "PCDEF  knot[{}]={:.17e} m={}",
                                i + 1,
                                knots[i],
                                mults[i]
                            );
                        }
                    }
                    let mut flat = String::new();
                    for k in &b.knots {
                        flat.push_str(&format!("{k:.17e},"));
                    }
                    eprintln!("PCDEF  flat={flat}");
                }
            }
            // TEMP EINFO: edge/pcurve provenance for the determinism hunt.
            if let Ok(spec) = std::env::var("EINFO") {
                let wanted: Vec<usize> = spec
                    .split(',')
                    .filter_map(|s| s.trim().parse::<usize>().ok())
                    .collect();
                let all = spec.trim().eq_ignore_ascii_case("ALL")
                    || spec
                        .split(',')
                        .any(|s| s.trim().eq_ignore_ascii_case("ALL"));
                if all || wanted.contains(&edge_index) {
                    let reg = crate::tgeometry::GeometryRegistry::global();
                    let face_key = crate::tgeometry::GeometryRegistry::shape_key(&topo_face.0);
                    let (ea, eb) = reg.edge_parameters(&topo_edge.0);
                    let c3 = reg.edge_curve(&topo_edge.0);
                    let c3kind = match c3.as_ref() {
                        None => "none",
                        Some(c) if c.is_line() => "line",
                        Some(c) if c.gp_circ().is_some() => "circle",
                        Some(c) if c.nurbs_degree().is_some() => "bspline",
                        Some(_) => "other",
                    };
                    let (c3f, c3l) = c3
                        .as_ref()
                        .map(|c| (c.first_parameter(), c.last_parameter()))
                        .unwrap_or((f64::NAN, f64::NAN));
                    let nstored = reg.edge_pcurves(&topo_edge.0, face_key).len();
                    let crange = reg.pcurve_range(&topo_edge.0, face_key);
                    let pckind = if pc.is_geom2d_line() {
                        "line2d"
                    } else if pc.is_bspline2d() {
                        "bspline2d"
                    } else {
                        "other2d"
                    };
                    eprintln!(
                        "EINFO f={pc_face} e={edge_index} pc={pcurve_index} ea={ea:.17e} eb={eb:.17e} c3f={c3f:.17e} c3l={c3l:.17e} pf={pf:.17e} pl={pl:.17e} pcf={:.17e} pcl={:.17e} nstored={nstored} crange={crange:?} c3={c3kind} pc2={pckind} deg={:?} skind={:?} uper={} vper={} urange={:?} vrange={:?} nparam={} t0={:.17e} t1={:.17e}",
                        pc.first_parameter(),
                        pc.last_parameter(),
                        pc.bspline_degree(),
                        crate::pcurve_full::classify_surface_kind(surface.as_ref()),
                        surface.is_u_periodic(),
                        surface.is_v_periodic(),
                        surface.u_range(),
                        surface.v_range(),
                        params.len(),
                        params.first().copied().unwrap_or(f64::NAN),
                        params.last().copied().unwrap_or(f64::NAN),
                    );
                }
            }
            let mut provider = if same_param {
                EdgeParameterProvider::new(pf, pl)
            } else {
                let old_first = params.first().copied().unwrap_or(a);
                let old_last = params.last().copied().unwrap_or(b);
                EdgeParameterProvider::with_stored(pf, pl, old_first, old_last)
            };
            let cos = CurveOnSurface::new(pc.clone(), surface.clone(), pf, pl);
            for (&t, p3d) in params.iter().zip(pts3d.iter()) {
                // `BRepMesh_EdgeParameterProvider::Parameter` (`hxx:109-139`)
                // + `Tessellate2d` (`BRepMesh_EdgeDiscret.cxx:312-317`).
                let u = provider.parameter_of(t, p3d, &cos);
                pcurve.add_point(pc.d0(u), t);
            }
            // TEMP PCW: full precision pcurve point dump, to be lined up with
            // the OCCT probe's `OCCTEPC w= e= pcId= k= u= v= t=` entries.
            if std::env::var("PCW").is_ok() {
                let n_pc = pcurve.parameters_nb();
                for k in 0..n_pc {
                    let p = pcurve.point(k)?;
                    let t = pcurve.parameter(k)?;
                    eprintln!(
                        "PCW f={pc_face} e={edge_index} pc={pcurve_index} k={k} u={:.17e} v={:.17e} t={:.17e}",
                        p.x(),                        p.y(),
                        t
                    );
                }
            }
        }
        if let Ok(spec) = std::env::var("EPROBE") {
            let parts: Vec<f64> =
                spec.split(',').filter_map(|s| s.trim().parse().ok()).collect();
            if parts.len() == 4 {
                for i in 0..model.edges_nb() {
                    let e = model.edge(i)?;
                    let pts = e.discretization().points();
                    let hit = pts.iter().any(|p| {
                        (p.x() - parts[0]).abs() <= parts[3]
                            && (p.y() - parts[1]).abs() <= parts[3]
                            && (p.z() - parts[2]).abs() <= parts[3]
                    });
                    if hit {
                        eprintln!(
                            "EPROBE edge={} np3d={} defl={:e} pcurves={} same_param={} free={}",
                            i,
                            pts.len(),
                            e.deflection(),
                            e.pcurves_nb(),
                            e.same_param(),
                            e.is_free()
                        );
                        for p in 0..e.pcurves_nb() {
                            let pc = e.pcurve(p)?;
                            eprintln!(
                                "EPROBE   pc={} face={} ori={:?} n={} interior={}",
                                p,
                                pc.face(),
                                pc.orientation(),
                                pc.parameters_nb(),
                                pc.is_internal()
                            );
                            eprint!("EPROBE   uv=");
                            for (i, uv) in pc.points().iter().enumerate() {
                                if i > 0 {
                                    eprint!(";");
                                }
                                eprint!("{:.17e},{:.17e}", uv.x(), uv.y());
                            }
                            eprintln!();
                            eprint!("EPROBE   tt=");
                            for (i, t) in pc.parameters().iter().enumerate() {
                                if i > 0 {
                                    eprint!(";");
                                }
                                eprint!("{:.17e}", t);
                            }
                            eprintln!();
                        }
                        eprintln!("EPROBE   pts={:?}", &pts[..pts.len().min(40)]);
                    }
                }
            }
        }
        Ok(())
    }

}
