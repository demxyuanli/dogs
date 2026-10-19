//! Face emitters: 142 curve-on-surface, 143 bounded surface and 144 trimmed surface, with the per-edge UV transfer. Source: `BRepToIGES_BRShell::TransferFace` and `BRepToIGES_BRWire::TransferEdge`.
use super::prelude::*;
use super::*;

impl IgesWriter {

    /// 'Geom2dToIGES_Geom2dCurve::Transfer2dCurve' +
    /// 'GeomToIGES_GeomCurve::TransferCurve': lift the 2-D curve to the plane
    /// Z=0 and reuse the 3-D emitter of 'emit_curve_range'.
    pub(super) fn emit_curve2d_range(&mut self, curve: &dyn Curve2d, first: f64, last: f64) -> Option<usize> {
        let promoted = promote_curve2d(curve)?;
        let idx = self.emit_curve_range(promoted.as_ref(), first, last);
        if crate::debug_flags::iges_trace() {
            eprintln!(
                "TRACE c2d de={:?} lin2d={} circ2d={} elips2d={} bsp2d={} -> kind={:?} a={first} b={last}",
                idx,
                curve.gp_lin2d().is_some(),
                curve.gp_circ2d().is_some(),
                curve.gp_elips2d().is_some(),
                curve.bspline_poles2d().is_some(),
                iges_curve_kind(promoted.as_ref()),
            );
        }
        idx
    }

    /// 'BRepToIGES_BRWire::TransferEdge(edge, face, originMap, length, false)'
    /// ('BRepToIGES_BRWire.cxx:340-585'): the edge's p-curve, corrected for the
    /// face's surface type and emitted as an IGES curve entity.
    ///
    /// 'theOriginMap' has no port equivalent - OCCT only uses it to look up the
    /// origin edge of a reversed face before 'SetShapeResult', whose transfer-map
    /// bookkeeping this port does not keep - so it is treated as empty.
    pub(super) fn transfer_edge_uv(
        &mut self,
        e: &Edge,
        f: &Face,
        uv_bounds: (f64, f64, f64, f64),
        reversal: Option<&FaceReversal>,
    ) -> Option<usize> {
        // 'cxx:348-352': 'GetPCurveMode() == 0' cannot happen (the default is
        // On) and 'theIsBRepMode' is false, so only the degeneracy guard applies.
        if BRepTool::is_degenerated(e) {
            return None;
        }
        // 'cxx:374-377' reads the surface of the face it is handed, which for a
        // REVERSED face is `TransferFace`'s U-reversed copy (`cxx:135-146`).
        let surf = match reversal {
            Some(r) => r.surf.clone(),
            None => BRepTool::face_surface(f)?,
        };
        // 'cxx:374-377': only a bare 'Geom_Plane' returns here. A
        // 'Geom_RectangularTrimmedSurface' over a plane is *not* a plane for
        // OCCT's 'IsKind'; the port's trimmed surface delegates 'gp_pln', so the
        // raw type is required.
        if surf.rectangular_trimmed_basis().is_none() && surf.gp_pln().is_some() {
            return None;
        }
        let (ufirst, _ulast, vfirst, vlast) = uv_bounds; // 'BRepTools::UVBounds' ('cxx:379')
        // 'cxx:380-397': peel a rectangular trim, then (from the raw surface
        // handle) an offset surface.
        let mut surf_base: Arc<dyn Surface> = match surf.rectangular_trimmed_basis() {
            Some(b) => b,
            None => surf.clone(),
        };
        if surf.is_offset_surface() {
            if let Some(b) = surf.offset_basis_surface() {
                surf_base = b;
            }
        }
        let is_cyl = surf_base.gp_cylinder().is_some();
        let is_cone = surf_base.gp_cone().is_some();
        let is_sphere = surf_base.gp_sphere().is_some();
        let is_torus = surf_base.gp_torus().is_some();
        let is_rev = surf_base.is_surface_of_revolution();
        let is_extr = surf_base.is_surface_of_linear_extrusion();
        let is_bspline = surf_base.is_bspline_surface();

        // 'cxx:359-360': the p-curve and the range of its CurveOnSurface
        // representation.
        let (curve2d, mut first, mut last) = crate::boptools_2d::curve_on_surface_range(e, f)?;
        if let Ok(path) = std::env::var("IGES_TRACE3D") {
            let (e3a, e3b) = BRepTool::edge_parameters(e);
            let c3 = BRepTool::edge_curve(e);
            let (sa, sb) = surf_base.u_range();
            let (ta, tb) = surf_base.v_range();
            if let Some(k) = c3.as_ref().and_then(|c| c.gp_circ()) {
                use std::io::Write;
                if let Ok(mut f) = std::fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(&path)
                {
                    let _ = writeln!(
                        f,
                        "e3=({e3a:.17e},{e3b:.17e}) c3kind={:?} rad={:.17e} | uv2=({first:.17e},{last:.17e}) | u_rng=({sa:.17e},{sb:.17e}) v_rng=({ta:.17e},{tb:.17e}) up={} sph={is_sphere} rev={is_rev} loc=({:.17e},{:.17e},{:.17e}) ax=({:.17e},{:.17e},{:.17e})",
                        crate::iges::iges_curve_kind(c3.as_ref().unwrap().as_ref()),
                        k.radius(),
                        surf_base.is_u_periodic(),
                        k.position().location().x(), k.position().location().y(), k.position().location().z(),
                        k.position().direction().x(), k.position().direction().y(), k.position().direction().z(),
                    );
                }
            }
        }
        if crate::debug_flags::iges_trace() {
            eprintln!(
                "TRACE raw pcurve lin2d={} circ2d={} elips2d={} bsp2d={} trimmed={} nbpol={} deg={:?} rng=({first},{last}) rev={} | surf pln={} cyl={is_cyl} cone={is_cone} sph={is_sphere} tor={is_torus} rev={is_rev} extr={is_extr} bsp={is_bspline} rect={} off={} kind={:?}",
                curve2d.gp_lin2d().is_some(),
                curve2d.gp_circ2d().is_some(),
                curve2d.gp_elips2d().is_some(),
                curve2d.bspline_poles2d().is_some(),
                curve2d.trimmed_basis().is_some(),
                curve2d.bspline_poles2d().map_or(0, |(x, _)| x.len()),
                curve2d.bspline_degree(),
                e.0.orientation() == crate::abs::Orientation::Reversed,
                surf.gp_pln().is_some(),
                surf.rectangular_trimmed_basis().is_some(),
                surf.is_offset_surface(),
                crate::pcurve_full::classify_surface_kind(surf_base.as_ref()),
            );
            eprintln!(
                "TRACE stored pcurve={}",
                crate::boptools_2d::curve_on_surface(e, f).is_some()
            );
            if let Ok(pcs) = std::env::var("IGES_TRACE_PC") {
                let fk = crate::tgeometry::GeometryRegistry::shape_key(&f.0);
                let list = crate::tgeometry::GeometryRegistry::global().edge_pcurves(&e.0, fk);
                let kinds: Vec<String> = list
                    .iter()
                    .map(|c| {
                        format!(
                            "lin{} circ{} elips{} bsp{} np{} trim{}",
                            c.gp_lin2d().is_some() as u8,
                            c.gp_circ2d().is_some() as u8,
                            c.gp_elips2d().is_some() as u8,
                            c.bspline_poles2d().is_some() as u8,
                            c.bspline_poles2d().map_or(0, |(x, _)| x.len()),
                            c.trimmed_basis().is_some() as u8
                        )
                    })
                    .collect();
                use std::io::Write;
                if let Ok(mut fp) = std::fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(&pcs)
                {
                    let _ = writeln!(
                        fp,
                        "PC face={:?} edge={:?} rev={} n={} [{}]",
                        Arc::as_ptr(&f.0.tshape) as usize,
                        Arc::as_ptr(&e.0.tshape) as usize,
                        e.0.orientation() == crate::abs::Orientation::Reversed,
                        list.len(),
                        kinds.join(" ")
                    );
                }
            }
        }

        // `BRepToIGES_BRShell.cxx:186-201`: on a REVERSED face every p-curve was
        // mirrored in `TransferFace` (`aCurve1->Transformed(T)`, `:200`) before
        // the edge was re-attached to the copied face, so the curve read here is
        // the mirrored one. `Transformed` keeps the curve's own parameter range
        // (`Geom2d_Transformed` forwards `FirstParameter`/`LastParameter`), so
        // `first`/`last` stay as they are.
        let curve2d: Arc<dyn Curve2d> = match reversal {
            Some(r) => {
                let axis = GpAx2d::new(
                    GpPnt2d::new(0.5 * r.two_center, vfirst),
                    GpDir2d::new(0.0, 1.0).expect("gp::DY2d is non-null"),
                );
                curve2d_transformed(&curve2d, &mirror_ax2d_trsf(&axis))
            }
            None => curve2d,
        };

        // 'cxx:403-422': 'analyticMode' is false, so the '!analyticMode' guard
        // holds; a surface of revolution whose (trim-unwrapped) basis curve is a
        // line also needs the shift.
        let mut need_shift = is_cyl || is_cone;
        if is_rev {
            if let Some(c) = surf_base.revolution_basis_curve() {
                let c: Arc<dyn Curve> = if c.is_geom_trimmed() {
                    c.untrimmed_basis().map(|(b, _, _)| b).unwrap_or(c)
                } else {
                    c
                };
                if c.is_line() {
                    need_shift = true;
                }
            }
        }
        let mut curve2d: Arc<dyn Curve2d> = if need_shift {
            // 'cxx:425-428': translate by '-Vfirst'.
            curve2d_transformed(&curve2d, &translate2d(0.0, -vfirst))
        } else {
            // 'cxx:431': 'Curve2d->Copy()'.
            Arc::from(curve2d.clone_dyn())
        };

        // 'cxx:435-455': a periodic B-spline surface is brought back into its
        // own period.
        if is_bspline {
            let (su0, su1) = surf_base.u_range();
            let (sv0, sv1) = surf_base.v_range();
            let mut u_shift = 0.0;
            let mut v_shift = 0.0;
            if surf_base.is_u_periodic() && (ufirst - su0).abs() > occt_core::precision::PCONFUSION {
                u_shift = crate::pcurve_full::adjust_to_period(ufirst, su0, su1);
            }
            if surf_base.is_v_periodic() && (vfirst - sv0).abs() > occt_core::precision::PCONFUSION {
                v_shift = crate::pcurve_full::adjust_to_period(vfirst, sv0, sv1);
            }
            if u_shift.abs() > occt_core::precision::PCONFUSION
                || v_shift.abs() > occt_core::precision::PCONFUSION
            {
                curve2d = curve2d_transformed(&curve2d, &translate2d(u_shift, v_shift));
            }
        }

        // 'cxx:457-474': the IGES surface of revolution inverts (u, v). The
        // cylinder / cone / sphere arm sits behind the '!analyticMode' guard,
        // which is true here; the revolution / torus arm has no such guard.
        if is_cyl || is_cone || is_sphere {
            curve2d = curve2d_transformed(&curve2d, &mirror_origin_dir11());
            curve2d = curve2d_transformed(&curve2d, &mirror_ox2d());
            curve2d = curve2d_transformed(
                &curve2d,
                &translate2d(0.0, 2.0 * std::f64::consts::PI),
            );
        }
        if is_rev || is_torus {
            curve2d = curve2d_transformed(&curve2d, &mirror_origin_dir11());
            curve2d = curve2d_transformed(&curve2d, &mirror_ox2d());
            curve2d = curve2d_transformed(
                &curve2d,
                &translate2d(0.0, 2.0 * std::f64::consts::PI),
            );
        }

        // UNPORTED ('cxx:476-502', all under 'analyticMode'): the cylinder /
        // cone 'myLen = PI/180' arm, the sphere / torus '180/PI' scaling and the
        // cone-apex translation are unreachable because every 'TransferEdge'
        // call from the face path passes 'theIsBRepMode = false'.
        // 'cxx:504-538': the 'theIsBRepMode && Surf->IsKind(Geom_Plane)' branch
        // is false too, so 'trans' only carries the extrusion scale.
        let len_val = surface_transfer_length(surf_base.as_ref(), vfirst, vlast);
        let mut trans = GpTrsf2d::identity();
        let mut u_fact = 1.0;
        if is_extr {
            // 'cxx:510-532' (emv, bug OCC22126): scale by '1/(Vlast-Vfirst)',
            // then 'uFact = (Vlast - Vfirst) / (us2 - us1)'.
            let _ = trans.set_scale(&GpPnt2d::zero(), 1.0 / (vlast - vfirst));
            let (us1, us2) = surf_base.u_range();
            let du = us2 - us1;
            u_fact = (vlast - vfirst) / du;
        }
        if is_cyl || is_cone || is_rev {
            // 'cxx:533-538': 'uFact = 1. / myLen'.
            u_fact = 1.0 / len_val;
        }
        curve2d = transform_pcurve(&curve2d, &trans, u_fact, &mut first, &mut last)?;

        // 'cxx:546-556': a second 'TransformPCurve' pass on a surface of linear
        // extrusion, shifting the p-curve range onto [0, 1] in u and v.
        if is_extr {
            let (us1, us2) = surf_base.u_range();
            let du = us2 - us1;
            let trans1 = translate2d(-us1 / du, -vfirst / (vlast - vfirst));
            curve2d = transform_pcurve(&curve2d, &trans1, 1.0, &mut first, &mut last)?;
        }

        // 'cxx:558-565': a REVERSED edge reverses the 2-D curve.
        if e.0.orientation() == crate::abs::Orientation::Reversed {
            let tmp_first = curve2d.reversed_parameter(last);
            let tmp_last = curve2d.reversed_parameter(first);
            curve2d = Arc::from(curve2d.reversed());
            first = tmp_first;
            last = tmp_last;
        }

        // 'cxx:566-568': 'Geom2dToIGES_Geom2dCurve::Transfer2dCurve'.
        if crate::debug_flags::iges_trace() {
            let a = curve2d.d0(first);
            let b = curve2d.d0(last);
            eprintln!(
                "TRACE emit pcurve lin2d={} circ2d={} elips2d={} bsp2d={} np={} rng=({first},{last}) d0=({},{}) d1=({},{})",
                curve2d.gp_lin2d().is_some(),
                curve2d.gp_circ2d().is_some(),
                curve2d.gp_elips2d().is_some(),
                curve2d.bspline_poles2d().is_some(),
                curve2d.bspline_poles2d().map_or(0, |(x, _)| x.len()),
                a.x(),
                a.y(),
                b.x(),
                b.y()
            );
        }
        self.emit_curve2d_range(curve2d.as_ref(), first, last)
    }

    /// A face as `BRepToIGES_BRShell::TransferFace` builds it
    /// (`BRepToIGES_BRShell.cxx:250-405`): the base surface entity, one
    /// `CurveOnSurface` (142) per wire, and a `TrimmedSurface` (144) that carries
    /// the outer/inner contour pointers.
    ///
    /// The 3-D curve of a wire comes from `BRepToIGES_BRWire::TransferWire`
    /// (`BRepToIGES_BRWire.cxx:662-781`): a single-edge wire is that edge's curve
    /// entity, a wire with two or more edges a `CompositeCurve` (102, form 0).
    /// The 2-D (UV) curve comes from the same `TransferWire` loop (R2-19): OCCT
    /// calls `TransferEdge(edge, face, originMap, length, false)` (`:720`) per
    /// edge and collects the results exactly as it does the 3-D curves
    /// (`:755-774`). The `142` preference is 3 when both curves transferred, 2
    /// for the 3-D-only case and 1 for UV-only (`BRepToIGES_BRShell.cxx:281-292`).
    ///
    /// **Registered remainder (T-84)**: `510` still carries this port's earlier
    /// layout instead of `IGESSolid_ToolFace::WriteOwnParams`
    /// (`510, surface, nb_loops, has_outer_loop, loop_ptrs...`) over `508` Loop
    /// entities, and `514`/`186` follow the same old shape.
    pub(super) fn emit_face(&mut self, f: &Face) -> Option<usize> {
        // `BRepToIGES_BRShell.cxx:257-261`: a face whose basis surface
        // transferred to a null entity is dropped (`AddWarning(start, "the basic
        // surface is a null entity"); return res;`).
        //
        // `cxx:126-237`: a REVERSED face first goes through the U-reversed
        // surface / mirrored p-curve branch; `face_reversal` reports whether that
        // branch applies (and is `None` for the orientations OCCT leaves alone).
        let reversal = face_reversal(f);
        let surf_idx = self.emit_face_surface(f, reversal.as_ref())?;
        let wires = wires_of_face(f);
        if crate::debug_flags::iges_trace_wire() {
            for w in &wires {
                let raw = w.0.tshape.read().expect("poisoned TShape lock").children.clone();
                eprintln!("WIRE ori={:?} n={}", w.0.orientation(), raw.len());
                for (i, c) in raw.iter().enumerate() {
                    let e = Edge(c.clone());
                    let (v1, v2) = edge_vertices(&e);
                    let f = |v: &Option<Vertex>| -> String {
                        match v {
                            Some(v) => {
                                let q = BRepTool::vertex_point(v);
                                format!("({:.4},{:.4},{:.4})", q.x(), q.y(), q.z())
                            }
                            None => "(?)".into(),
                        }
                    };
                    eprintln!("  child{} ori={:?} {} -> {}", i + 1, c.orientation(), f(&v1), f(&v2));
                }
                for (i, e) in reordered_wire_edges(w).iter().enumerate() {
                    let (v1, v2) = edge_vertices(e);
                    let f = |v: &Option<Vertex>| -> String {
                        match v {
                            Some(v) => {
                                let q = BRepTool::vertex_point(v);
                                format!("({:.4},{:.4},{:.4})", q.x(), q.y(), q.z())
                            }
                            None => "(?)".into(),
                        }
                    };
                    eprintln!("  out{} ori={:?} {} -> {}", i + 1, e.0.orientation(), f(&v1), f(&v2));
                }
                {
                    let stored = wire_data_edges(w);
                    let mut order = crate::meshing::wire_order::WireOrder::new();
                    for e in &stored {
                        let (a, b) = (wire_edge_first_vertex(e), wire_edge_last_vertex(e));
                        if let (Some(a), Some(b)) = (a, b) {
                            let (pa, pb) = (BRepTool::vertex_point(&a), BRepTool::vertex_point(&b));
                            eprintln!("    in begin=({:.3},{:.3},{:.3}) end=({:.3},{:.3},{:.3})", pa.x(), pa.y(), pa.z(), pb.x(), pb.y(), pb.z());
                            order.add_edge_xyz(pa, pb);
                        }
                    }
                    order.perform();
                    let chain: Vec<i32> = (1..=stored.len()).map(|i| order.ordered(i)).collect();
                    eprintln!("    ORDER status={:?} chain={:?}", order.status(), chain);
                }
            }
        }
        if crate::debug_flags::iges_trace_face() {
            let b = face_uv_bounds_finite(f);
            let kind = match BRepTool::face_surface(f) {
                Some(s) if s.gp_cylinder().is_some() => "cyl",
                Some(s) if s.gp_cone().is_some() => "cone",
                Some(s) if s.gp_sphere().is_some() => "sph",
                Some(s) if s.gp_torus().is_some() => "tor",
                Some(s) if s.gp_pln().is_some() => "pln",
                Some(s) if s.is_bspline_surface() => "bspl",
                Some(_) => "other",
                None => "none",
            };
            eprintln!(
                "TRACEFACE ptr={:?} de={:?} ori={:?} rev={} kind={kind} uv=({:.9},{:.9},{:.9},{:.9})",
                Arc::as_ptr(&f.0.tshape) as usize,
                surf_idx,
                f.0.orientation(),
                reversal.is_some(),
                b.0,
                b.1,
                b.2,
                b.3
            );
        }
        // `BRepTools::UVBounds(aFace, U1, U2, V1, V2)` (`BRepToIGES_BRShell.cxx:253`),
        // kept once for every `TransferEdge` call of this face (OCCT recomputes
        // them inside `TransferEdge`, `BRepToIGES_BRWire.cxx:379`).
        // A bare plane face never reaches `TransferEdge`'s UV computation
        // (`BRepToIGES_BRWire.cxx:374-377`), so the (non-trivial) bounds are only
        // derived for the faces that need them.
        let uv_bounds = match BRepTool::face_surface(f) {
            Some(s) if !(s.rectangular_trimmed_basis().is_none() && s.gp_pln().is_some()) => {
                face_uv_bounds_finite(f)
            }
            _ => (0.0, 1.0, 0.0, 1.0),
        };
        // `ShapeAlgo::AlgoContainer()->OuterWire(aFace)` (`cxx:275`) resolves to
        // `ShapeAnalysis::OuterWire` (`ShapeAnalysis.cxx`: the last wire, or the
        // first one with a non-negative `TotCross2D`).
        let outer_wire = crate::meshing::model_builder::outer_of_wires(&wires, f);
        let mut curve_refs: Vec<usize> = Vec::new();
        let mut outer_curve: Option<usize> = None;
        let mut inner_curves: Vec<usize> = Vec::new();
        // `BRepToIGES_BRShell.cxx:270`: `Iprefer` is initialised once and only
        // updated by the wires / free edges that transfer (OCCT never resets it).
        // The OCCT initial value is 0; every 142 below assigns it before reading,
        // so the initialiser is dead code but kept for fidelity.
        #[allow(unused_assignments)]
        let mut iprefer = 0i32;
        // `TransferWire` only writes `theCurve2d` when at least one edge
        // produced a UV curve (`cxx:755-774`), so a wire with none leaves the
        // caller's handle untouched (`BRepToIGES_BRShell.cxx:271`); the port
        // carries it the same way.
        #[allow(unused_assignments)]
        let mut carried_uv: Option<usize> = None;
        // The outer wire is transferred first (`cxx:276-294`), then the inner ones
        // (`:301-331`).
        let ordered: Vec<&crate::shape::Wire> = outer_wire
            .iter()
            .chain(wires.iter().filter(|w| {
                outer_wire
                    .as_ref()
                    .map_or(true, |o| Arc::as_ptr(&o.0.tshape) != Arc::as_ptr(&w.0.tshape))
            }))
            .collect();
        for w in ordered {
            let mut edge_refs: Vec<usize> = Vec::new();
            let mut uv_refs: Vec<usize> = Vec::new();
            let mut last_uv: Option<usize> = None;
            for e in reordered_wire_edges(w) {
                // `BRepToIGES_BRWire::TransferEdge` (`BRepToIGES_BRWire.cxx:266-335`)
                // transfers the edge's 3-D curve into a fresh entity on every call:
                // `Curve3d->Copy()` (`:300`) and `GeomToIGES_GeomCurve::TransferCurve`
                // (`:316`) create a new `IGESGeom_Line` / arc / 126 each time. OCCT
                // shares nothing between the wires of a face, the wires of different
                // faces, or the two occurrences of a seam edge, so the port must not
                // deduplicate here either (Cube's 6 faces x 4 edges = 24 lines).
                let idx = self.emit_edge_curve(&e);
                if let Some(idx) = idx {
                    edge_refs.push(idx);
                    curve_refs.push(idx);
                }
                // `BRepToIGES_BRWire::TransferWire` (`cxx:715-724`): every edge's
                // 2-D curve, collected in the same order as the 3-D ones.
                let uv = self.transfer_edge_uv(&e, f, uv_bounds, reversal.as_ref());
                last_uv = uv;
                if let Some(u) = uv {
                    uv_refs.push(u);
                }
            }
            if edge_refs.is_empty() && uv_refs.is_empty() {
                continue;
            }
            let curve3d = if edge_refs.is_empty() {
                None
            } else if edge_refs.len() == 1 {
                Some(edge_refs[0])
            } else {
                Some(self.emit_refs(102, 0, format!("102,{},{};", edge_refs.len(), ph_list(edge_refs.len())), &edge_refs))
            };
            // `cxx:755-774`: one 2-D entity is used directly, two or more become
            // a `CompositeCurve` (102). The single-entity arm mirrors
            // `theCurve2d = ent2d` (`cxx:760`), i.e. the last edge's transfer
            // result rather than `Seq2d(1)`; zero 2-D entities leave the caller's
            // handle from the previous wire in place.
            let curve_uv = if uv_refs.len() == 1 {
                last_uv
            } else if uv_refs.len() >= 2 {
                Some(self.emit_refs(102, 0, format!("102,{},{};", uv_refs.len(), ph_list(uv_refs.len())), &uv_refs))
            } else {
                carried_uv
            };
            carried_uv = curve_uv;
            // `BRepToIGES_BRShell.cxx:281-292`: the preference follows which of
            // the two representations transferred; both null leaves it at its
            // previous value.
            if curve3d.is_some() {
                iprefer = if curve_uv.is_some() { 3 } else { 2 };
            } else if curve_uv.is_some() {
                iprefer = 1;
            }
            // `IGESGeom_CurveOnSurface::Init` (`IGESGeom_CurveOnSurface.cxx:26-40`)
            // with `Imode = 0` (`cxx:269`);
            // `IGESGeom_ToolCurveOnSurface::WriteOwnParams` (`:144-152`) writes
            // `142, creation_mode, surface, curve_uv, curve_3d, preference_mode;`
            // - field 3 is the UV curve, field 4 the 3-D one.
            let cs = match (curve_uv, curve3d) {
                (Some(uv), Some(c3)) => {
                    self.emit_refs(142, 0, format!("142,0,#0,#1,#2,{iprefer};"), &[surf_idx, uv, c3])
                }
                (Some(uv), None) => {
                    self.emit_refs(142, 0, format!("142,0,#0,#1,0,{iprefer};"), &[surf_idx, uv])
                }
                (None, Some(c3)) => {
                    self.emit_refs(142, 0, format!("142,0,#0,0,#1,{iprefer};"), &[surf_idx, c3])
                }
                (None, None) => {
                    self.emit_refs(142, 0, format!("142,0,#0,0,0,{iprefer};"), &[surf_idx])
                }
            };
            match &outer_wire {
                Some(o) if Arc::as_ptr(&o.0.tshape) == Arc::as_ptr(&w.0.tshape) => {
                    outer_curve = Some(cs)
                }
                _ => inner_curves.push(cs),
            }
        }
        // `cxx:334-365`: edges of the face that are not part of any wire become
        // further inner contours. OCCT transfers both their 3-D curve (`:344`) and
        // their UV curve (`TransferEdge(edge, face, originMap, length, false)`,
        // `:346`).
        let wire_edge_keys: std::collections::HashSet<usize> = wires
            .iter()
            .flat_map(|w| edges_of_wire(w))
            .map(|e| Arc::as_ptr(&e.0.tshape) as usize)
            .collect();
        for e in children_of_type(&f.0, ShapeType::Edge) {
            let key = Arc::as_ptr(&e.tshape) as usize;
            if wire_edge_keys.contains(&key) {
                continue;
            }
            let edge = Edge(e);
            // Each free edge is its own `TransferEdge` result too
            // (`BRepToIGES_BRShell.cxx:344-346`), with no sharing.
            let idx = self.emit_edge_curve(&edge);
            if let Some(idx) = idx {
                curve_refs.push(idx);
            }
            let uv = self.transfer_edge_uv(&edge, f, uv_bounds, reversal.as_ref());
            // `BRepToIGES_BRShell.cxx:347-358`: the preference follows which of
            // the two representations transferred.
            if idx.is_some() {
                iprefer = if uv.is_some() { 3 } else { 2 };
            } else if uv.is_some() {
                iprefer = 1;
            }
            let cs = match (uv, idx) {
                (Some(u), Some(c3)) => {
                    self.emit_refs(142, 0, format!("142,0,#0,#1,#2,{iprefer};"), &[surf_idx, u, c3])
                }
                (Some(u), None) => {
                    self.emit_refs(142, 0, format!("142,0,#0,#1,0,{iprefer};"), &[surf_idx, u])
                }
                (None, Some(c3)) => {
                    self.emit_refs(142, 0, format!("142,0,#0,0,#1,{iprefer};"), &[surf_idx, c3])
                }
                (None, None) => {
                    self.emit_refs(142, 0, format!("142,0,#0,0,0,{iprefer};"), &[surf_idx])
                }
            };
            inner_curves.push(cs);
        }

        // `cxx:380-400`: `isWholeSurface` is `BRep_Tool::NaturalRestriction(face)`,
        // forced to false for a plane / cylinder / cone - the guard there tests the
        // `CurveOnSurface` handle, which is non-null whenever the face has a wire.
        let mut is_whole = BRepTool::natural_restriction(f);
        let surf = BRepTool::face_surface(f);
        if let Some(s) = surf.as_ref() {
            let k = classify_surface(s.as_ref());
            if (k == SurfaceKind::Plane || k == SurfaceKind::Cylinder || k == SurfaceKind::Cone)
                && outer_curve.is_some()
            {
                is_whole = false;
            }
        }
        // A face this port built without boundary wires carries no contour on
        // this side either; OCCT's `TransferWire` leaves the caller's handle
        // untouched when no edge transfers (`BRepToIGES_BRWire.cxx:715-774`).

        if curve_refs.is_empty() {
            return Some(surf_idx);
        }
        let outer_flag = u8::from(!is_whole);
        let n_inner = inner_curves.len();
        // `IGESGeom_ToolTrimmedSurface::WriteOwnParams`
        // (`IGESGeom_ToolTrimmedSurface.cxx:196-218`): `144, surface,
        // outer_boundary_type, nb_inner_contours, outer_contour, inner...;`
        // with the outer contour written as `0` when the type is false.
        // `surface` and the contours are pointers (`#k`), the boundary type and
        // the inner-contour count are plain integers.
        let has_outer = !is_whole && outer_curve.is_some();
        let mut ref_list = vec![surf_idx];
        if has_outer {
            ref_list.push(outer_curve.unwrap());
        }
        ref_list.extend(inner_curves.iter().copied());
        let mut contour_refs = if has_outer { ph(1) } else { "0".to_string() };
        for k in if has_outer { 2 } else { 1 }..ref_list.len() {
            contour_refs.push_str(&format!(",{}", ph(k)));
        }
        Some(self.emit_refs(
            144,
            0,
            format!("144,#0,{outer_flag},{n_inner},{contour_refs};"),
            &ref_list,
        ))
    }
}
