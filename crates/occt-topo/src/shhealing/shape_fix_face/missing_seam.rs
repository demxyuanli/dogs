//! `ShapeFix_Face::FixMissingSeam` (`ShapeFix_Face.cxx:1722-2330`).

use super::prelude::*;
use super::*;

impl ShapeFixFace {
    /// `ShapeFix_Face::FixMissingSeam` (`ShapeFix_Face.cxx:1722-2330`), ported
    /// through the setup and wire-pair selection (`cxx:1722-1898`).
    ///
    /// UNPORTED: the seam construction and the `w2 != null` merge
    /// (`cxx:1899-2330`) — port them before wiring `ShapeFix_Face::Perform`.
    pub fn fix_missing_seam(&mut self) -> bool {
        let surf = match &self.surf {
            Some(s) => s.clone(),
            None => return false, // cxx:1724-1727
        };
        let mut face = match &self.face {
            Some(f) => f.clone(),
            None => return false,
        };
        // cxx:1729-1735.
        let uclosed = crate::pcurve_full::sa_is_u_closed(surf.as_ref(), CONFUSION);
        let vclosed = crate::pcurve_full::sa_is_v_closed(surf.as_ref(), CONFUSION);
        if !uclosed && !vclosed {
            return false;
        }

        // cxx:1737-1741: `myFace = TopoDS::Face(Context()->Apply(myFace))`.
        // This is what hands the face built by `FixPeriodicDegenerated`
        // (`cxx:3108-3111` applies the context as well) to the wire-pair
        // selection below.
        let applied = self.context.apply(&face.0);
        if applied.is_face() {
            face = Face(applied);
            self.face = Some(face.clone());
        }

        // cxx:1744-1751: a BSpline surface must be U- or V-periodic.
        if surf.is_bspline_surface() && !surf.is_u_periodic() && !surf.is_v_periodic() {
            return false;
        }

        // cxx:1753-1756.
        let (mut suf, mut sul) = surf.u_range();
        let (mut svf, mut svl) = surf.v_range();
        let (f_u1, f_u2, f_v1, f_v2) = brep_tools::uv_bounds(&face);

        // cxx:1758-1802: replace infinite surface bounds with the face's.
        if Precision::is_infinite(suf) || Precision::is_infinite(sul) {
            if Precision::is_infinite(suf) {
                suf = f_u1;
            }
            if Precision::is_infinite(sul) {
                sul = f_u2;
            }
            if (sul - suf).abs() < PCONFUSION {
                if Precision::is_infinite(suf) {
                    suf -= 1000.0;
                } else {
                    sul += 1000.0;
                }
            }
        }
        if Precision::is_infinite(svf) || Precision::is_infinite(svl) {
            if Precision::is_infinite(svf) {
                svf = f_v1;
            }
            if Precision::is_infinite(svl) {
                svl = f_v2;
            }
            if (svl - svf).abs() < PCONFUSION {
                if Precision::is_infinite(svf) {
                    svf -= 1000.0;
                } else {
                    svl += 1000.0;
                }
            }
        }

        // cxx:1804-1805.
        let u_range = (sul - suf).abs().min(INFINITE);
        let v_range = (svl - svf).abs().min(INFINITE);

        // cxx:1810-1822: `TopoDS_Iterator(myFace, false)` — oriented wires go
        // to `ws`, every other direct child (non-wire or non-oriented) to
        // `aSeqNonManif`.
        let children: Vec<TopoShape> = {
            let ts = face.0.tshape.read().unwrap();
            ts.children.clone()
        };
        let mut ws: Vec<Wire> = Vec::new();
        let mut non_manifold: Vec<TopoShape> = Vec::new();
        for child in children {
            let o = child.orientation();
            if child.shape_type() == ShapeType::Wire
                && (o == Orientation::Forward || o == Orientation::Reversed)
            {
                ws.push(Wire(child));
            } else {
                non_manifold.push(child);
            }
        }

        // cxx:1824-1881: select the wire pair.
        let mut w1: Option<Wire> = None;
        let mut w2: Option<Wire> = None;
        let mut ismodeu = 0i32;
        let mut ismodev = 0i32;
        let mut isdeg1 = 0usize;
        let mut isdeg2 = 0usize;
        let mut i = 0usize;
        while i < ws.len() {
            let wire = ws[i].clone();
            let (isuopen, isvopen, isdeg) = match healing_check_wire(&wire, &face, u_range, v_range) {
                Some(v) => v,
                None => {
                    i += 1;
                    continue; // cxx:1831-1834
                }
            };
            if w1.is_none() {
                ismodeu = isuopen;
                ismodev = isvopen;
                isdeg1 = if isdeg { i + 1 } else { 0 };
                w1 = Some(wire);
            } else if w2.is_none() {
                if ismodeu == -isuopen && ismodev == -isvopen {
                    isdeg2 = if isdeg { i + 1 } else { 0 };
                    w2 = Some(wire);
                } else if ismodeu == isuopen && ismodev == isvopen {
                    // cxx:1849-1865.
                    isdeg2 = usize::from(isdeg);
                    w2 = Some(wire);
                    if isdeg1 != 0 {
                        if let Some(w) = w1.as_mut() {
                            w.0.reverse();
                        }
                        ismodeu = -ismodeu;
                        ismodev = -ismodev;
                    } else if let Some(w) = w2.as_mut() {
                        w.0.reverse();
                    }
                }
            } else if isdeg || isdeg1 != 0 || isdeg2 != 0 {
                // cxx:1868-1879.
                let rem = if isdeg {
                    i + 1
                } else if isdeg2 != 0 {
                    isdeg2
                } else {
                    isdeg1
                };
                ws.remove(rem - 1);
                w1 = None;
                w2 = None;
                i = 0;
                continue;
            }
            i += 1;
        }

        // cxx:1883-1893: a torus with MajorRadius < MinorRadius counts as
        // degenerated, unless a second wire exists.
        let torus = surf.gp_torus();
        let mut an_is_degenerated_tor = false;
        if let Some(t) = &torus {
            an_is_degenerated_tor = t.major_radius < t.minor_radius;
        }
        if an_is_degenerated_tor && w2.is_some() {
            an_is_degenerated_tor = false;
        }

        // cxx:1895-1898.
        if w1.is_none() {
            return false;
        }
        if w2.is_none() {
            // cxx:1899-1992: add a degenerated edge for spheres and BSpline
            // cone-like surfaces, then handle them as a regular pair.
            let mut builder_pt = GpPnt2d::new(0.0, 0.0);
            let dir;
            let range;
            if ismodeu != 0 && an_is_degenerated_tor {
                let t = torus.as_ref().expect("torus");
                let ra = t.major_radius;
                let ri = t.minor_radius;
                let phi = (-ra / ri).acos();
                builder_pt = GpPnt2d::new(
                    0.0,
                    if ismodeu > 0 { std::f64::consts::PI + phi } else { phi },
                );
                dir = GpDir2d::new(-(ismodeu as f64), 0.0).expect("dir");
                range = 2.0 * std::f64::consts::PI;
            } else if ismodeu != 0 && surf.gp_sphere().is_some() {
                builder_pt = GpPnt2d::new(
                    if ismodeu < 0 { 0.0 } else { 2.0 * std::f64::consts::PI },
                    ismodeu as f64 * 0.5 * std::f64::consts::PI,
                );
                dir = GpDir2d::new(-(ismodeu as f64), 0.0).expect("dir");
                range = 2.0 * std::f64::consts::PI;
            } else if ismodev != 0 && surf.is_bspline_surface() {
                let u_coord;
                if surf.d0(suf, svf).distance(&surf.d0(suf, (svf + svl) / 2.0)) < CONFUSION {
                    u_coord = suf;
                } else if surf.d0(sul, svf).distance(&surf.d0(sul, (svf + svl) / 2.0)) < CONFUSION {
                    u_coord = sul;
                } else {
                    return false; // cxx:1940-1943
                }
                builder_pt = GpPnt2d::new(u_coord, if ismodev < 0 { 0.0 } else { v_range });
                dir = GpDir2d::new(0.0, -(ismodev as f64)).expect("dir");
                range = v_range;
            } else if ismodeu != 0 && surf.is_bspline_surface() {
                let v_coord;
                if surf.d0(suf, svf).distance(&surf.d0((suf + sul) / 2.0, svf)) < CONFUSION {
                    v_coord = svf;
                } else if surf.d0(sul, svl).distance(&surf.d0((suf + sul) / 2.0, svl)) < CONFUSION {
                    v_coord = svl;
                } else {
                    return false; // cxx:1962-1965
                }
                builder_pt = GpPnt2d::new(if ismodeu < 0 { 0.0 } else { u_range }, v_coord);
                dir = GpDir2d::new(-(ismodeu as f64), 0.0).expect("dir");
                range = u_range;
            } else {
                return false; // cxx:1972-1975
            }

            // cxx:1977-1991.
            let builder = TopoBuilder::new();
            let line: Arc<dyn Curve2d> = Arc::new(Geom2dLine::from_pnt_dir(builder_pt, dir));
            let reg = GeometryRegistry::global();
            let face_key = GeometryRegistry::shape_key(&face.0);
            let mut edge = builder.make_shape(ShapeType::Edge);
            reg.set_degenerated(&edge, true);
            reg.set_edge_pcurve(&edge, face_key, line);
            reg.set_pcurve_range(&edge, face_key, 0.0, range);
            reg.set_edge_tolerance(&edge, CONFUSION);
            let mut v = builder.make_vertex(surf.d0(builder_pt.x(), builder_pt.y()), CONFUSION);
            v.0.set_orientation(Orientation::Forward);
            builder.add(&mut edge, &v.0);
            v.0.set_orientation(Orientation::Reversed);
            builder.add(&mut edge, &v.0);
            let new_w2 = builder.make_wire(&[Edge(edge)]);
            ws.push(new_w2.clone());
            w2 = Some(new_w2);
        }

        // cxx:1994-2006: UV bounds of the two wires on the face's surface.
        let wire_uv_bounds = |w: &Wire| -> (f64, f64, f64, f64) {
            let mut b = BndBox2d::new();
            add_uv_bounds_on_wire(&face, w, &mut b);
            match b.get() {
                Some((u0, v0, u1, v1)) => (u0, u1, v0, v1),
                None => (0.0, 0.0, 0.0, 0.0),
            }
        };
        let (w1, w2) = (w1.unwrap(), w2.unwrap());
        let m1 = wire_uv_bounds(&w1);
        let m2 = wire_uv_bounds(&w2);
        let coord = if ismodeu != 0 { 1usize } else { 0usize };
        let isneg = if ismodeu != 0 { ismodeu } else { -ismodev };
        let period = if ismodeu != 0 { u_range } else { v_range };
        // `m1[coord]` / `m2[coord]` in `ShapeFix_Face.cxx:1996-2127`: the arrays
        // are filled by `ShapeAnalysis::GetFaceUVBounds(F, UMin, UMax, VMin,
        // VMax)` (cxx:268-273), so `m[..][0]` holds the U bounds and `m[..][1]`
        // the V bounds, while `coord = (ismodeu ? 1 : 0)` — `coord == 1`
        // therefore selects **V**, not U. `pick(m, 1)` returns the V pair.
        let pick = |m: &(f64, f64, f64, f64), c: usize| -> (f64, f64) {
            if c == 1 { (m.2, m.3) } else { (m.0, m.1) }
        };
        let mut w1 = w1;
        let mut w2 = w2;
        if !vclosed || !uclosed || an_is_degenerated_tor {
            let (a0, a1) = pick(&m1, coord);
            let (b0, b1) = pick(&m2, coord);
            let delta_other = 0.5 * (b0 + b1) - 0.5 * (a0 + a1);
            if delta_other * (isneg as f64) < 0.0 {
                w1.0.reverse(); // cxx:2016-2020
                w2.0.reverse();
            }
        }

        // UNPORTED (cxx:2029-2032): `ShapeFix_Wire::FixReorder` — the port's
        // `fix_reorder_wire` only reports the 3D order, it does not rewrite the
        // stored edge list, so `w11`/`w21` keep the original order.
        let _ = fix_reorder_wire(&w1, &face);
        let w11 = w1.clone();
        let w21 = w2.clone();

        // cxx:2036-2066: rebuild the face, replacing w1/w2 with w11/w21 and
        // orienting the other wires as holes.
        let builder = TopoBuilder::new();
        let mut tmp_wires: Vec<Wire> = Vec::new();
        for w in &ws {
            let mut wire = w.clone();
            if crate::topo_tools::is_same(&wire.0, &w1.0) {
                wire = w11.clone();
            } else if crate::topo_tools::is_same(&wire.0, &w2.0) {
                wire = w21.clone();
            } else {
                // cxx:2052-2062.
                let cur_face = builder.make_face(surf.clone(), &[wire.clone()]);
                let mut cur_face = cur_face;
                cur_face.0.set_orientation(face.0.orientation());
                if is_outer_bound(&cur_face) {
                    wire.0.reverse();
                }
            }
            tmp_wires.push(wire);
        }
        let mut tmp_f = builder.make_face(surf.clone(), &tmp_wires);
        tmp_f.0.set_orientation(face.0.orientation());

        // cxx:2068-2136: torus-like FixShifted along the missing seam direction.
        let m1 = m1;
        let mut m1 = m1;
        let m2 = m2;
        let mut uf2 = suf;
        let mut vf2 = svf;
        if uclosed && vclosed && !an_is_degenerated_tor {
            let shiftw2 = adjust_by_period(
                0.5 * (pick(&m2, coord).0 + pick(&m2, coord).1),
                0.5 * (pick(&m1, coord).0 + pick(&m1, coord).1
                    + (isneg as f64) * (period + PCONFUSION)),
                period,
            );
            let (m2c0, m2c1) = pick(&m2, coord);
            let (m1c0, m1c1) = pick(&m1, coord);
            // cxx:2077-2078: `m1[coord][0] = min(...)`, `m1[coord][1] = max(...)`
            // — the update lands on whichever pair `coord` selects (V when
            // `ismodeu`), not on U.
            if coord == 1 {
                m1.2 = m1c0.min(m2c0 + shiftw2);
                m1.3 = m1c1.max(m2c1 + shiftw2);
            } else {
                m1.0 = m1c0.min(m2c0 + shiftw2);
                m1.1 = m1c1.max(m2c1 + shiftw2);
            }
            let tmp_face_wires = wires_of_face(&tmp_f);
            for w in &tmp_face_wires {
                if crate::topo_tools::is_same(&w.0, &w11.0) {
                    continue;
                }
                let shift = if crate::topo_tools::is_same(&w.0, &w21.0) {
                    shiftw2
                } else {
                    let mw = wire_uv_bounds(w);
                    adjust_by_period(
                        0.5 * (pick(&mw, coord).0 + pick(&mw, coord).1),
                        0.5 * (pick(&m1, coord).0 + pick(&m1, coord).1),
                        period,
                    )
                };
                if shift != 0.0 {
                    // cxx:2107 `V.SetCoord(coord + 1, shift)`: `coord == 1`
                    // writes gp_Vec2d's Y (V) component, `coord == 0` its X (U).
                    let vshift = if coord == 1 {
                        GpVec2d::new(0.0, shift)
                    } else {
                        GpVec2d::new(shift, 0.0)
                    };
                    let tkey = GeometryRegistry::shape_key(&tmp_f.0);
                    for e in edges_of_wire(w) {
                        // `sa.PCurve(E, tmpF, C, a, b)` (`ShapeFix_Face.cxx:2114`):
                        // surface overload, so the face orientation does not
                        // select the seam representation.
                        if let Some((c2d, _a, _b)) = boptools_2d::curve_on_surface_oriented(&e, &tmp_f, false) {
                            let mut nc = c2d.clone_dyn();
                            let mut tr = occt_core::gp::GpTrsf2d::identity();
                            tr.set_translation_vec(&vshift);
                            nc.transform(&tr);
                            // `C->Translate(V)` (`cxx:2118`) moves the one
                            // representation `sae.PCurve` returned; the other
                            // seam representation keeps its values.
                            let reg = GeometryRegistry::global();
                            let mut pcs = reg.edge_pcurves(&e.0, tkey);
                            let nc: Arc<dyn Curve2d> = Arc::from(nc);
                            if pcs.len() >= 2 {
                                if e.0.orientation().is_reversed() {
                                    pcs[1] = nc;
                                } else {
                                    pcs[0] = nc;
                                }
                                reg.set_edge_pcurves(&e.0, tkey, pcs);
                            } else {
                                reg.set_edge_pcurve(&e.0, tkey, nc);
                            }
                        }
                    }
                }
            }
            // cxx:2122-2135.
            if pick(&m1, coord).1 - pick(&m1, coord).0 <= period {
                let other = 0.5 * (pick(&m1, coord).0 + pick(&m1, coord).1 - period);
                if ismodeu != 0 {
                    vf2 = other;
                } else {
                    uf2 = other;
                }
            }
        }

        // cxx:2138-2234: find the best place to insert the seam.
        let mut found_u = 0i32;
        let mut found_v = 0i32;
        let wd1_edges = edges_of_wire(&w11);
        let wd2_edges = edges_of_wire(&w21);
        let nb1 = wd1_edges.len();
        let nb2 = wd2_edges.len();
        let mut i1 = 1usize;
        while i1 <= nb1 + nb2 {
            let edge1 = if i1 <= nb1 {
                wd1_edges[i1 - 1].clone()
            } else {
                wd2_edges[i1 - nb1 - 1].clone()
            };
            // `sa.PCurve(edge1, tmpF, c2d, f, l, true)` (`ShapeFix_Face.cxx:2149`).
            let (c2d, f, l) = match boptools_2d::curve_on_surface_oriented(&edge1, &tmp_f, true) {
                Some(v) => v,
                None => return false, // cxx:2149-2152
            };
            let mut pos1 = c2d.d0(l);
            let mut skip_u = !uclosed;
            if uclosed && ismodeu != 0 {
                pos1 = GpPnt2d::new(pos1.x() + adjust_by_period(pos1.x(), suf, u_range), pos1.y());
                if found_u == 2 && pos1.x().abs() > uf2.abs() {
                    skip_u = true;
                } else if found_u == 0 || (found_u == 1 && pos1.x().abs() < uf2.abs()) {
                    found_u = 1;
                    uf2 = pos1.x();
                }
            }
            let mut skip_v = !vclosed;
            if vclosed && ismodeu == 0 {
                pos1 = GpPnt2d::new(pos1.x(), pos1.y() + adjust_by_period(pos1.y(), svf, v_range));
                if found_v == 2 && pos1.y().abs() > vf2.abs() {
                    skip_v = true;
                } else if found_v == 0 || (found_v == 1 && pos1.y().abs() < vf2.abs()) {
                    found_v = 1;
                    vf2 = pos1.y();
                }
            }
            if skip_u && skip_v {
                if i1 <= nb1 {
                    i1 += 1;
                    continue;
                } else {
                    break;
                }
            }
            for i2 in 1..=nb2 {
                if i1 > nb1 {
                    break;
                }
                let edge2 = wd2_edges[i2 - 1].clone();
                // `sa.PCurve(edge2, tmpF, c2d, f, l, true)` (`ShapeFix_Face.cxx:2198`).
                let (c2d2, f2, _l2) = match boptools_2d::curve_on_surface_oriented(&edge2, &tmp_f, true) {
                    Some(v) => v,
                    None => return false, // cxx:2198-2201
                };
                let mut pos2 = c2d2.d0(f2);
                if uclosed && ismodeu != 0 {
                    pos2 = GpPnt2d::new(pos2.x() + adjust_by_period(pos2.x(), pos1.x(), u_range), pos2.y());
                    if (pos2.x() - pos1.x()).abs() < PCONFUSION
                        && (found_u != 2 || pos1.x().abs() < uf2.abs())
                    {
                        found_u = 2;
                        uf2 = pos1.x();
                    }
                }
                if vclosed && ismodeu == 0 {
                    pos2 = GpPnt2d::new(pos2.x(), pos2.y() + adjust_by_period(pos2.y(), pos1.y(), v_range));
                    if (pos2.y() - pos1.y()).abs() < PCONFUSION
                        && (found_v != 2 || pos1.y().abs() < vf2.abs())
                    {
                        found_v = 2;
                        vf2 = pos1.y();
                    }
                }
            }
            i1 += 1;
        }

        // cxx:2226-2234.
        if uf2 < suf || uf2 > sul {
            uf2 += compose_adjust_to_period(uf2, suf, suf + u_range);
        }
        if vf2 < svf || vf2 > svl {
            vf2 += compose_adjust_to_period(vf2, svf, svf + v_range);
        }

        // cxx:2236-2261: fictive grid + ComposeShell inserts the seam.
        let rts: Arc<dyn Surface> = Arc::new(
            GeomRectangularTrimmedSurface::uv(
                surf.clone(),
                uf2,
                uf2 + u_range,
                vf2,
                vf2 + v_range,
            )
            .unwrap_or_else(|e| panic!("{e}")),
        );
        let grid = CompositeSurface::with_grid(vec![vec![rts.clone()]], Parametrisation::Natural);
        // cxx:2246-2250: re-add the non-manifold children.
        for s in &non_manifold {
            builder.add(&mut tmp_f.0, s);
        }
        let mut comp = ComposeShell::new();
        comp.init(grid, &tmp_f, CONFUSION); // cxx:2252-2253
        comp.set_closed_mode(true); // cxx:2258
        comp.set_context(self.context.clone()); // cxx:2259
        comp.set_max_tolerance(self.max_tol); // cxx:2260
        comp.perform(); // cxx:2261
        // cxx:2263-2264: reset mySurf to the trimmed surface.
        self.surf = Some(rts);
        // cxx:2266 myResult = CompShell.Result().
        let mut result = comp.result().cloned();
        // cxx:2270-2322: drop the small wires / faces ComposeShell can generate.
        //  * cxx:2273-2300: every wire of every result face goes through
        //    ShapeFix_Wire::FixSmall(true, Precision()); a wire that loses all
        //    its edges is discarded and a face that loses all its wires is
        //    removed (Context()->Remove).
        //  * cxx:2303-2321: when the result keeps more than one face
        //    (nbFaces > 1), FixSmallAreaWire(true) (cxx:2331-2394) discards
        //    every wire for which ShapeAnalysis_Wire::CheckSmallArea (cxx:2004)
        //    reports a null area; a face whose wires all disappear is removed.
        // UNPORTED: Context()->Apply (cxx:2302/2322) is expressed by rebuilding
        // the result here; BRepTools::Update (cxx:2320) is a no-op.
        if let Some(res) = &result {
            if res.shape_type() != ShapeType::Face {
                let builder = TopoBuilder::new();
                // First round: FixSmall.
                let mut round1: Vec<(Face, Vec<Wire>)> = Vec::new();
                for f in crate::topo_tools::faces_of(res) {
                    let mut kept: Vec<Wire> = Vec::new();
                    for w in wires_of_face(&f) {
                        let mut w2 = w.clone();
                        crate::shhealing::fix_small_all(&mut w2, &f, CONFUSION, CONFUSION, false);
                        if !edges_of_wire(&w2).is_empty() {
                            kept.push(w2);
                        }
                    }
                    if !kept.is_empty() {
                        round1.push((f, kept));
                    }
                }
                let nb_faces = round1.len();
                let mut kept_faces: Vec<Face> = Vec::new();
                for (f, wires) in round1 {
                    let mut kept_wires = wires;
                    if nb_faces > 1 {
                        kept_wires.retain(|w| !crate::shhealing::check_small_area(w, &f));
                    }
                    if kept_wires.is_empty() {
                        continue; // cxx:2377-2385 / 2294-2298: remove the face
                    }
                    let orig_n = wires_of_face(&f).len();
                    if kept_wires.len() == orig_n {
                        kept_faces.push(f);
                    } else if let Some(s) = BRepTool::face_surface(&f) {
                        let mut nf = builder.make_face(s, &kept_wires);
                        nf.0.set_orientation(f.0.orientation());
                        kept_faces.push(nf);
                    }
                }
                result = match kept_faces.len() {
                    0 => None,
                    1 => Some(kept_faces.remove(0).0.clone()),
                    _ => Some(builder.make_shell(&kept_faces).0),
                };
            }
        }
        self.result = result;
        // cxx:2268: Context()->Replace(myFace, myResult).
        if let Some(res) = &self.result {
            self.context.replace(&face.0, res);
        }
        true
    }
}
