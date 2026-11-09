use super::prelude::*;
use super::statics::*;
use super::fill_contours::fill_contours;
use super::rebuild_connecting_edge::rebuild_connecting_edge;

impl BRepOffsetMakeOffset {
    /// `MakeMissingWalls(theRange)` (`cxx:3148-3687`).
    pub fn make_missing_walls(&mut self) {
        let reg = GeometryRegistry::global();
        let builder = TopoBuilder::new();
        let offset_val = self.my_offset.abs();

        // `cxx:3151-3157`.
        let mut contours: Vec<(TopoShape, Vec<TopoShape>)> = Vec::new();
        let mut map_ef: HashMap<ShapeMapKey, TopoShape> = HashMap::new();
        fill_contours(&self.my_face_comp.0, &self.my_analyse, &mut contours, &mut map_ef);

        let mut walls: Vec<TopoShape> = Vec::new();
        let mut perform_sewing = self.my_is_perform_sewing;

        // `cxx:3160-3186`.
        for (_start_key, start_vertex_shape, contour) in contours.iter().map(|(v, c)| {
            (shape_map_key(v), v.clone(), c.clone())
        }) {
            let Some(start_vertex) = Vertex::wrap(start_vertex_shape) else {
                continue;
            };
            let mut first_step = true;
            let mut prev_edge: Option<Edge> = None;
            let mut prev_vertex = start_vertex.0.clone();
            let mut start_edge: Option<Edge> = None;
            let mut build_from_scratch = false;

            for edge_shape in contour {
                let Some(an_edge) = Edge::wrap(edge_shape.clone()) else {
                    continue;
                };
                let Some(a_face_of_edge) = map_ef
                    .get(&shape_map_key(&an_edge.0))
                    .and_then(|f| Face::wrap(f.clone()))
                else {
                    continue;
                };

                // `cxx:3181-3184`.
                if !self.my_init_offset_edge.has_image(&an_edge.0) {
                    continue;
                }

                // `cxx:3187-3193`: the offset image of the edge and its
                // endpoints.
                let mut loe: Vec<TopoShape> = Vec::new();
                self.my_init_offset_edge.last_image(&an_edge.0, &mut loe);
                let Some(loe_last) = loe.last() else {
                    continue;
                };
                let mut loe2: Vec<TopoShape> = Vec::new();
                self.my_image_offset.last_image(loe_last, &mut loe2);
                let Some(oe_shape) = loe2.last() else {
                    continue;
                };
                let Some(oe) = Edge::wrap(oe_shape.clone()) else {
                    continue;
                };
                let (v4, v3) = edge_vertices(&oe);
                let (v1, v2) = edge_vertices(&an_edge);
                let (Some(v1), Some(v2)) = (v1, v2) else {
                    continue;
                };
                let a_c = reg.edge_curve(&an_edge.0);
                let (a_f, a_l) = BRepTool::edge_parameters(&an_edge);

                // `cxx:3196-3225`: an offset edge that got no vertices.
                let (v3, v4) = if v3.is_none() && v4.is_none() {
                    let oe_ori = oe.orientation();
                    let mut oe_work = oe.0.clone();
                    oe_work.set_orientation(Orientation::Forward);
                    let (an_oe_f, an_oe_l) = BRepTool::edge_parameters(&Edge(oe_work.clone()));
                    let Some(an_oe_c) = reg.edge_curve(&oe_work) else {
                        continue;
                    };
                    let p1 = an_oe_c.d0(a_f);
                    let p2 = an_oe_c.d0(a_l);
                    let a_tol = reg
                        .vertex_tolerance(&v1.0)
                        .max(reg.vertex_tolerance(&v2.0));
                    let mut new_v1 = builder.make_vertex(p1, a_tol);
                    new_v1.0.set_orientation(Orientation::Forward);
                    let mut new_v2 = builder.make_vertex(p2, a_tol);
                    new_v2.0.set_orientation(Orientation::Reversed);
                    builder.add(&mut oe_work, &new_v1.0);
                    builder.add(&mut oe_work, &new_v2.0);
                    reg.set_ranges_all(&oe_work, an_oe_f, an_oe_l);
                    oe_work.set_orientation(oe_ori);
                    edge_vertices(&Edge(oe_work))
                } else {
                    (v3, v4)
                };

                // `cxx:3226-3258`: a non-closed input edge whose endpoints
                // collapsed on one point.
                if let Some(a_c) = &a_c {
                    if !a_c.is_closed() && !a_c.is_periodic() {
                        let pnt_f = BRepTool::vertex_point(&v1);
                        let pnt_l = BRepTool::vertex_point(&v2);
                        let dist_e = pnt_f.square_distance(&pnt_l);
                        if dist_e < SQUARE_CONFUSION {
                            continue;
                        }
                        let edge_tol = reg.edge_tolerance(&an_edge.0);
                        if dist_e < edge_tol {
                            if let (Some(ov4), Some(ov3)) = (v4.as_ref(), v3.as_ref()) {
                                let pnt_of = BRepTool::vertex_point(ov4);
                                let pnt_ol = BRepTool::vertex_point(ov3);
                                if pnt_of.square_distance(&pnt_ol) > RESOLUTION {
                                    let sin_value =
                                        (2.0 * edge_tol / pnt_of.distance(&pnt_ol)).min(1.0);
                                    let max_angle = sin_value.asin().abs().min(std::f64::consts::FRAC_PI_4);
                                    let current_angle = occt_core::gp::GpVec::from_pnts(&pnt_f, &pnt_l)
                                        .angle(&occt_core::gp::GpVec::from_pnts(&pnt_of, &pnt_ol));
                                    if a_c.is_line() && current_angle.abs() > max_angle {
                                        // `cxx:3252-3254`: not collinear with the
                                        // offset edge.
                                        build_from_scratch = true;
                                        perform_sewing = true;
                                        continue;
                                    }
                                }
                            }
                        }
                    }
                }

                // `cxx:3260-3270`: the offset edge runs the other way.
                let mut to_reverse = false;
                let (mut v1, mut v2, mut v3, mut v4) = (v1, v2, v3, v4);
                if !v1.0.same_tshape(&prev_vertex) {
                    std::mem::swap(&mut v1, &mut v2);
                    std::mem::swap(&mut v3, &mut v4);
                    to_reverse = true;
                }

                // `cxx:3272`.
                let mut oe_work = oe.0.clone();
                oe_work.set_orientation(an_edge.orientation().reversed());

                // `cxx:3274`.
                let arc_on_v2 = self.my_join == occt_core::kernel::geomabs::JoinType::Arc
                    && self.my_init_offset_edge.has_image(&v2.0);

                // `cxx:3275-3294`.
                let mut e4 = if first_step || build_from_scratch {
                    let Some(v4) = v4.as_ref() else {
                        continue;
                    };
                    let Some(e) = make_edge_between(&v1, v4) else {
                        continue;
                    };
                    if first_step {
                        start_edge = Some(e.clone());
                    }
                    e
                } else {
                    let Some(e) = prev_edge.clone() else {
                        continue;
                    };
                    e
                };
                let mut e3 = if v2.0.same_tshape(&start_vertex.0) && !arc_on_v2 {
                    let Some(e) = start_edge.clone() else {
                        continue;
                    };
                    e
                } else {
                    let Some(v3) = v3.as_ref() else {
                        continue;
                    };
                    let Some(e) = make_edge_between(&v2, v3) else {
                        continue;
                    };
                    e
                };
                // `cxx:3295`.
                e4.0.set_orientation(e4.0.orientation().reversed());
                // `cxx:3297-3301`.
                if build_from_scratch {
                    e3.0.set_orientation(e3.0.orientation().reversed());
                    e4.0.set_orientation(e4.0.orientation().reversed());
                }

                // `cxx:3303-3306`.
                let mut edge_fwd = an_edge.0.clone();
                edge_fwd.set_orientation(Orientation::Forward);
                let par_v1 = BRepTool::parameter_of_vertex(&edge_fwd, &v1.0);
                let par_v2 = BRepTool::parameter_of_vertex(&edge_fwd, &v2.0);

                // `cxx:3307-3323`.
                let mut the_wire = Wire::new();
                let members: [TopoShape; 4] = if to_reverse {
                    [
                        oriented_by(&an_edge.0, Orientation::Reversed),
                        oriented_by(&e3.0, Orientation::Reversed),
                        oriented_by(&oe_work, Orientation::Reversed),
                        oriented_by(&e4.0, Orientation::Reversed),
                    ]
                } else {
                    [an_edge.0.clone(), e3.0.clone(), oe_work.clone(), e4.0.clone()]
                };
                for m in &members {
                    builder.add(&mut the_wire.0, m);
                }
                // `cxx:3325-3326`.
                let _ = build_curves3d_tol(&the_wire.0, self.my_tol, Continuity::C1, 14, 0);
                the_wire.0.set_closed(true);

                // `cxx:3327-3337`.
                let mut new_face: Option<Face> = None;
                let mut the_surf: Option<Arc<dyn Surface>> = None;
                let mut edge_line2d: Option<Geom2dLine> = None;
                let mut oe_line2d: Option<Geom2dLine> = None;
                let ba_curve = reg.edge_curve(&an_edge.0);
                let ba_curve_oe = reg.edge_curve(&oe_work);
                let Some(ba_curve) = ba_curve else {
                    continue;
                };
                let (fpar, lpar) = (ba_curve.first_parameter(), ba_curve.last_parameter());
                let pon_e = ba_curve.d0(fpar);
                let Some(ba_curve_oe) = ba_curve_oe else {
                    continue;
                };
                let pon_oe = ba_curve_oe.d0(fpar);
                let Ok(offset_dir) = GpDir::from_vec(&occt_core::gp::GpVec::from_pnts(&pon_e, &pon_oe))
                else {
                    continue;
                };
                let mut is_planar = false;

                // `cxx:3338-3459`.
                if let (Some(a_circ), Some(a_circ_oe)) = (ba_curve.gp_circ(), ba_curve_oe.gp_circ()) {
                    let an_axis_line = GpLin::new(a_circ.axis());
                    let mut circ_axis_dir = *a_circ.axis().direction();
                    let axes_parallel =
                        circ_axis_dir.is_parallel_tol(a_circ_oe.axis().direction(), CONFUSION);
                    let on_axis = an_axis_line.contains(&a_circ_oe.location(), CONFUSION);
                    if axes_parallel && on_axis {
                        if (a_circ.radius() - a_circ_oe.radius()).abs() <= CONFUSION {
                            // `cxx:3348-3350`: cylinder, `gce_MakeCylinder(gp_Circ)`
                            // (`gp_Cylinder(Circ.Position(), Circ.Radius())`).
                            if let Ok(cyl) =
                                GpCylinder::new(a_circ.position().to_ax3(), a_circ.radius())
                            {
                                the_surf = Some(Arc::new(GeomCylinder::new(cyl)));
                            }
                        } else if a_circ.location().distance(&a_circ_oe.location()) <= CONFUSION {
                            // `cxx:3351-3403`: a planar wall.
                            is_planar = true;
                            let pon_el = ba_curve.d0(lpar);
                            if pon_el.distance(&pon_e) <= PCONFUSION {
                                // `cxx:3364-3373`: the outer ring is the larger
                                // circle.
                                let (a_e1, a_e2) = if a_circ.radius() > a_circ_oe.radius() {
                                    (an_edge.0.clone(), oe_work.clone())
                                } else {
                                    (oe_work.clone(), an_edge.0.clone())
                                };
                                if let Ok(ax3) = GpAx3::from_point_dir(a_circ.location(), circ_axis_dir)
                                {
                                    let a_pl = GpPln::new(ax3);
                                    let mut wires = [Wire::new(), Wire::new()];
                                    let edges = [a_e1, a_e2];
                                    for i in 0..2 {
                                        builder.add_edge(&mut wires[i], &Edge(edges[i].clone()));
                                    }
                                    // `cxx:3381-3397`: orient each ring against
                                    // the classifier's hole answer.
                                    for i in 0..2 {
                                        let mut a_face = builder.make_face_plane(&a_pl);
                                        builder.add(&mut a_face.0, &wires[i].0);
                                        if let Ok(clsf) = FClass2d::new(&a_face, CONFUSION) {
                                            let is_hole = clsf.is_hole();
                                            if (is_hole && i == 0) || (!is_hole && i == 1) {
                                                let mut w = Wire::new();
                                                builder.add_edge(
                                                    &mut w,
                                                    &Edge(oriented_by(
                                                        &edges[i],
                                                        Orientation::Reversed,
                                                    )),
                                                );
                                                wires[i] = w;
                                            }
                                        }
                                    }
                                    // `cxx:3399-3401`.
                                    let mut face = builder.make_face_plane(&a_pl);
                                    builder.add(&mut face.0, &wires[0].0);
                                    builder.add(&mut face.0, &wires[1].0);
                                    new_face = Some(face);
                                }
                            }
                        } else {
                            // `cxx:3404-3416`: cone.
                            if let Some(mut cone) = gce_make_cone(
                                &a_circ.location(),
                                &a_circ_oe.location(),
                                a_circ.radius(),
                                a_circ_oe.radius(),
                            ) {
                                let mut the_ax3 = a_circ.position().to_ax3();
                                if circ_axis_dir.dot(cone.axis().direction()) < 0.0 {
                                    the_ax3.z_reverse();
                                    circ_axis_dir.reverse();
                                }
                                cone.set_position(&the_ax3);
                                the_surf = Some(Arc::new(GeomCone::new(cone)));
                            }
                        }
                        if !is_planar {
                            let Some(surf) = the_surf.clone() else {
                                continue;
                            };
                            // `cxx:3420-3426`.
                            let x2d = GpDir2d::new(1.0, 0.0).expect("DX2d");
                            let el = line2d_at(GpPnt2d::new(0.0, 0.0), x2d);
                            reg.update_edge_pcurve_on_surface(
                                &an_edge.0,
                                &surf,
                                Some(Arc::new(el.clone())),
                                CONFUSION,
                            );
                            let coeff = if offset_dir.dot(&circ_axis_dir) > 0.0 {
                                1.0
                            } else {
                                -1.0
                            };
                            let oel = line2d_at(GpPnt2d::new(0.0, offset_val * coeff), x2d);
                            reg.update_edge_pcurve_on_surface(
                                &oe_work,
                                &surf,
                                Some(Arc::new(oel.clone())),
                                CONFUSION,
                            );
                            let dir_y = GpDir2d::new(0.0, coeff).expect("coefficient");
                            let a_line2d = line2d_at(GpPnt2d::new(par_v2, 0.0), dir_y);
                            let a_line2d2 = line2d_at(GpPnt2d::new(par_v1, 0.0), dir_y);
                            edge_line2d = Some(el);
                            oe_line2d = Some(oel);
                            if e3.0.same_tshape(&e4.0) {
                                if coeff > 0.0 {
                                    reg.set_edge_pcurves_on_surface(
                                        &e3.0,
                                        &surf,
                                        vec![Arc::new(a_line2d), Arc::new(a_line2d2)],
                                        CONFUSION,
                                    );
                                } else {
                                    reg.set_edge_pcurves_on_surface(
                                        &e3.0,
                                        &surf,
                                        vec![Arc::new(a_line2d2), Arc::new(a_line2d)],
                                        CONFUSION,
                                    );
                                    // `cxx:3436-3442`: the wire is rebuilt.
                                    let mut w = Wire::new();
                                    builder.add(
                                        &mut w.0,
                                        &oriented_by(&an_edge.0, Orientation::Reversed),
                                    );
                                    builder.add(&mut w.0, &e4.0);
                                    builder.add(
                                        &mut w.0,
                                        &oriented_by(&oe_work, Orientation::Forward),
                                    );
                                    builder.add(&mut w.0, &e3.0);
                                    w.0.set_closed(true);
                                    the_wire = w;
                                }
                            } else {
                                // `cxx:3447-3454`.
                                reg.set_same_parameter(&e3.0, false);
                                reg.set_same_range(&e3.0, false);
                                reg.set_same_parameter(&e4.0, false);
                                reg.set_same_range(&e4.0, false);
                                reg.update_edge_pcurve_on_surface(
                                    &e3.0,
                                    &surf,
                                    Some(Arc::new(a_line2d)),
                                    CONFUSION,
                                );
                                reg.set_pcurve_range_on_surface(&e3.0, &surf, 0.0, offset_val);
                                reg.update_edge_pcurve_on_surface(
                                    &e4.0,
                                    &surf,
                                    Some(Arc::new(a_line2d2)),
                                    CONFUSION,
                                );
                                reg.set_pcurve_range_on_surface(&e4.0, &surf, 0.0, offset_val);
                            }
                            // `cxx:3456`.
                            new_face = Some(make_face_from_surface_and_wire(surf, &the_wire));
                        }
                    }
                }

                // `cxx:3460-3545`.
                if new_face.is_none() {
                    let an_edge_tol = reg.edge_tolerance(&an_edge.0);
                    let find_plane = FindSurface::new(&the_wire.0, an_edge_tol, true, false);
                    is_planar = false;
                    if find_plane.is_found() && find_plane.tolerance_reached() <= an_edge_tol {
                        let (f, l) = BRepTool::edge_parameters(&an_edge);
                        if let (Some(a_gc), Some(a_pln)) = (
                            reg.edge_curve(&an_edge.0),
                            find_plane.surface().and_then(|s| s.gp_pln()),
                        ) {
                            let a_max_dist =
                                compute_max_dist(&a_pln, a_gc.as_ref(), f, l);
                            if a_max_dist <= an_edge_tol {
                                let plane_surf = Arc::new(GeomPlane::new(a_pln.clone()));
                                let mf = make_face_from_surface_and_wire(plane_surf, &the_wire);
                                // `cxx:3478-3489`: raise the tolerance of every
                                // other edge of the wire to its plane distance.
                                for an_e in children_keep_orientation(&the_wire.0) {
                                    if an_e.shape_type() != ShapeType::Edge {
                                        continue;
                                    }
                                    if an_e.same_tshape(&an_edge.0) {
                                        continue;
                                    }
                                    let Some(gc) = reg.edge_curve(&an_e) else {
                                        continue;
                                    };
                                    let (ef, el) = BRepTool::edge_parameters(&Edge(an_e.clone()));
                                    let md = compute_max_dist(&a_pln, gc.as_ref(), ef, el);
                                    reg.raise_edge_tolerance(&an_e, md);
                                }
                                new_face = Some(mf);
                                is_planar = true;
                            }
                        }
                    }
                    // `cxx:3495-3544`: extrusion through the two curves.
                    if !is_planar {
                        let (ef, el) = BRepTool::edge_parameters(&an_edge);
                        let Some(edge_curve) = reg.edge_curve(&an_edge.0) else {
                            continue;
                        };
                        let (oef, oel) = BRepTool::edge_parameters(&Edge(oe_work.clone()));
                        let Some(offset_curve) = reg.edge_curve(&oe_work) else {
                            continue;
                        };
                        let sections: Vec<Arc<dyn occt_geom::curve::Curve>> = vec![
                            Arc::new(GeomTrimmedCurve::new(edge_curve, ef, el)),
                            Arc::new(GeomTrimmedCurve::new(offset_curve, oef, oel)),
                        ];
                        let Ok(surf) = generator_surface(&sections) else {
                            continue;
                        };
                        let (uf, ul) = surf.u_range();
                        let (vf, vl) = surf.v_range();
                        let y2d = GpDir2d::new(1.0, 0.0).expect("DX2d");
                        let el2d = line2d_at(GpPnt2d::new(0.0, vf), y2d);
                        reg.update_edge_pcurve_on_surface(
                            &an_edge.0,
                            &surf,
                            Some(Arc::new(el2d.clone())),
                            CONFUSION,
                        );
                        let oel2d = line2d_at(GpPnt2d::new(0.0, vl), y2d);
                        reg.update_edge_pcurve_on_surface(
                            &oe_work,
                            &surf,
                            Some(Arc::new(oel2d.clone())),
                            CONFUSION,
                        );
                        let uon_v1 = if to_reverse { ul } else { uf };
                        let uon_v2 = if to_reverse { uf } else { ul };
                        let dir_y = GpDir2d::new(0.0, 1.0).expect("DY2d");
                        let a_line2d = line2d_at(GpPnt2d::new(uon_v2, 0.0), dir_y);
                        let a_line2d2 = line2d_at(GpPnt2d::new(uon_v1, 0.0), dir_y);
                        edge_line2d = Some(el2d);
                        oe_line2d = Some(oel2d);
                        if e3.0.same_tshape(&e4.0) {
                            // `cxx:3521-3524`.
                            reg.set_edge_pcurves_on_surface(
                                &e3.0,
                                &surf,
                                vec![Arc::new(a_line2d), Arc::new(a_line2d2)],
                                CONFUSION,
                            );
                            if let Some(bspl_c34) = surf.u_iso_curve(uf) {
                                reg.set_edge_curve3d(&e3.0, bspl_c34, CONFUSION);
                            }
                            reg.set_ranges_all(&e3.0, vf, vl);
                        } else {
                            // `cxx:3528-3541`.
                            reg.set_same_parameter(&e3.0, false);
                            reg.set_same_range(&e3.0, false);
                            reg.set_same_parameter(&e4.0, false);
                            reg.set_same_range(&e4.0, false);
                            reg.update_edge_pcurve_on_surface(
                                &e3.0,
                                &surf,
                                Some(Arc::new(a_line2d)),
                                CONFUSION,
                            );
                            reg.set_pcurve_range_on_surface(&e3.0, &surf, vf, vl);
                            reg.update_edge_pcurve_on_surface(
                                &e4.0,
                                &surf,
                                Some(Arc::new(a_line2d2)),
                                CONFUSION,
                            );
                            reg.set_pcurve_range_on_surface(&e4.0, &surf, vf, vl);
                            if let Some(bspl_c3) = surf.u_iso_curve(uon_v2) {
                                reg.set_edge_curve3d(&e3.0, bspl_c3, CONFUSION);
                            }
                            reg.set_edge_range(&e3.0, vf, vl);
                            if let Some(bspl_c4) = surf.u_iso_curve(uon_v1) {
                                reg.set_edge_curve3d(&e4.0, bspl_c4, CONFUSION);
                            }
                            reg.set_edge_range(&e4.0, vf, vl);
                        }
                        // `cxx:3543`.
                        new_face = Some(make_face_from_surface_and_wire(surf, &the_wire));
                    }
                }

                // `cxx:3546-3601`: widen the two connecting edges in the
                // parameter range the original / offset curves disagree on.
                if !is_planar {
                    if let (Some(surf), Some(el2d), Some(oel2d)) = (
                        the_surf.clone(),
                        edge_line2d.as_ref(),
                        oe_line2d.as_ref(),
                    ) {
                        let fpar_oe = ba_curve_oe.first_parameter();
                        let lpar_oe = ba_curve_oe.last_parameter();
                        if (fpar - fpar_oe).abs() > CONFUSION {
                            let target = if to_reverse { &e3 } else { &e4 };
                            let fp2d = el2d.d0(fpar);
                            let fp2d_oe = oel2d.d0(fpar_oe);
                            rebuild_connecting_edge(target, &surf, fp2d, fp2d_oe);
                        }
                        if (lpar - lpar_oe).abs() > CONFUSION {
                            let target = if to_reverse { &e4 } else { &e3 };
                            let lp2d = el2d.d0(lpar);
                            let lp2d_oe = oel2d.d0(lpar_oe);
                            rebuild_connecting_edge(target, &surf, lp2d, lp2d_oe);
                        }
                    }
                }

                // `cxx:3603-3609`.
                if !is_planar {
                    if let Some(face) = new_face.as_ref() {
                        same_parameter_face(face, CONFUSION, false);
                    }
                }

                // `cxx:3610-3618`.
                let Some(mut face) = new_face else {
                    continue;
                };
                let an_or = orientation_of_edge_in_face(&an_edge, &a_face_of_edge);
                let or_in_new = orientation_of_edge_in_face(&an_edge, &face);
                if or_in_new != an_or.reversed() {
                    face.0.set_orientation(face.0.orientation().reversed());
                }
                walls.push(face.0.clone());

                // `cxx:3619-3683`.
                if arc_on_v2 {
                    let arc_images = self.my_init_offset_edge.image(&v2.0);
                    let Some(an_arc_shape) = arc_images.first() else {
                        continue;
                    };
                    let Some(mut an_arc) = Edge::wrap(an_arc_shape.clone()) else {
                        continue;
                    };
                    let (arc_v1, arc_v2) = edge_vertices(&an_arc);
                    let (Some(arc_v1), Some(arc_v2)) = (arc_v1, arc_v2) else {
                        continue;
                    };
                    let mut arc_reverse = false;
                    let (_arc_v1, arc_v2) = if !v3
                        .as_ref()
                        .map(|v| arc_v1.0.same_tshape(&v.0))
                        .unwrap_or(false)
                    {
                        arc_reverse = true;
                        (arc_v2, arc_v1)
                    } else {
                        (arc_v1, arc_v2)
                    };
                    // `cxx:3632-3639`.
                    let mut ea1 = e3.0.clone();
                    ea1.set_orientation(ea1.orientation().reversed());
                    if to_reverse {
                        ea1.set_orientation(ea1.orientation().reversed());
                    }
                    // `cxx:3641-3648`.
                    let mut ea2 = if v2.0.same_tshape(&start_vertex.0) {
                        let Some(e) = start_edge.clone() else {
                            continue;
                        };
                        e.0
                    } else {
                        let Some(e) = make_edge_between(&v2, &arc_v2) else {
                            continue;
                        };
                        e.0
                    };
                    // `cxx:3649-3654`.
                    an_arc
                        .0
                        .set_orientation(if arc_reverse { Orientation::Reversed } else { Orientation::Forward });
                    if ea1.orientation() == Orientation::Reversed {
                        an_arc.0.set_orientation(an_arc.0.orientation().reversed());
                    }
                    ea2.set_orientation(ea1.orientation().reversed());
                    let mut arc_wire = Wire::new();
                    builder.add(&mut arc_wire.0, &ea1);
                    builder.add(&mut arc_wire.0, &an_arc.0);
                    builder.add(&mut arc_wire.0, &ea2);
                    let _ = build_curves3d_tol(&arc_wire.0, self.my_tol, Continuity::C1, 14, 0);
                    arc_wire.0.set_closed(true);
                    // `cxx:3662-3664`.
                    let arc_face = make_face_from_wire(&arc_wire, true);
                    if let Some(arc_face) = arc_face {
                        walls.push(arc_face.0.clone());
                    }
                    // `cxx:3665-3668`.
                    prev_edge = Some(Edge(oriented_by(&ea2, Orientation::Forward)));
                    prev_vertex = v2.0.clone();
                } else if build_from_scratch {
                    // `cxx:3672-3677`.
                    prev_edge = Some(e4.clone());
                    prev_vertex = v1.0.clone();
                    build_from_scratch = false;
                } else {
                    // `cxx:3678-3682`.
                    prev_edge = Some(e3.clone());
                    prev_vertex = v2.0.clone();
                }
                // `cxx:3684`.
                first_step = false;
            }
        }

        // `cxx:3618` / `:3664`: `myWalls.Append(...)`.
        self.my_walls.extend(walls);
        self.my_is_perform_sewing = perform_sewing;
    }
}
