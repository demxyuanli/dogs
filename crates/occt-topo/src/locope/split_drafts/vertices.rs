//! Vertex transitions of the split face: the section edges between the drafted
//! walls (`LocOpe_SplitDrafts.cxx:471-866`).
//!
//! For every vertex of every edge in `MapEV` the reference decides what the two
//! walls do at that vertex:
//!
//! * **Choice 1** - the vertex joins two edges that are *not* on the original
//!   wire: the edge that carries it away from `F` is recovered from `Res`, and
//!   the section edge is the one already present, shortened;
//! * **Choice 2** - it joins two wire edges: the section edge is the line along
//!   `NormalFg` through the vertex;
//! * **Choice 3** - mixed: the section edge is the intersection of the two
//!   `MapEV` faces.
//!
//! The curve is then cut against `NewSg` and `NewSd`, and the two hits either
//! become a new edge in `MapEV` (when they differ) or collapse onto the vertex.

use std::sync::Arc;

use occt_core::gp::GpLin;
use occt_core::precision::{Precision, CONFUSION};
use occt_geom::extrema_pc::extrema_ext_pc_range_done;
use occt_geom::{Curve, GeomLine, Surface};
use occt_geom2d::curve::Curve2d;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::locope::shape_maps::ShapeValueMap;
use crate::locope::split_drafts::perform::intcs_perform;
use crate::locope::split_drafts::state::PerformState;
use crate::locope::split_drafts::statics::{contains, empty_copied_edge};
use crate::locope::split_drafts::SplitDrafts;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{explorer_shapes, is_same};

/// `MapSg.Bind(vtx, list); MapSg(vtx).Append(v)` (`cxx:817-828`): `Bind`
/// replaces whatever was bound, so the result is the single-element list.
fn bind_one(map: &mut ShapeValueMap<Vec<TopoShape>>, key: &TopoShape, v: &TopoShape) {
    map.bind(key, vec![v.clone()]);
}

impl SplitDrafts {
    /// `cxx:471-866`. `Ok(false)` is the reference's bare `return`.
    pub(crate) fn build_vertices(
        &mut self,
        f: &Face,
        st: &mut PerformState,
    ) -> Result<bool, String> {
        let b = TopoBuilder::new();
        let registry = GeometryRegistry::global();
        // `cxx:476`: `MapEV.Extent()` taken before vertices join the map.
        let nbedges = st.map_ev.extent();
        let f_faces: Vec<TopoShape> = self.my_map.find(&f.0).cloned().unwrap_or_default();

        for index in 1..=nbedges {
            let edg1 = st.map_ev.find_key(index).clone();
            for vtx in explorer_shapes(&edg1, ShapeType::Vertex) {
                // `cxx:479-482`.
                if st.map_ev.find_index(&vtx).is_some() {
                    continue;
                }

                // `cxx:489-507`: the other edge of `MapEV` that shares `vtx`.
                let mut j: Option<usize> = None;
                for jj in 1..=nbedges {
                    if jj == index {
                        continue;
                    }
                    let key = st.map_ev.find_key(jj).clone();
                    if explorer_shapes(&key, ShapeType::Vertex)
                        .iter()
                        .any(|v| is_same(&v, &vtx))
                    {
                        j = Some(jj);
                        break;
                    }
                }

                // `cxx:508-517`.
                let edg2 = match j {
                    Some(jj) => st.map_ev.find_key(jj).clone(),
                    None => edg1.clone(),
                };

                // `cxx:518-546`.
                let w1 = st.map_w.is_bound(&edg1);
                let choice = match j {
                    None => {
                        if w1 {
                            2
                        } else {
                            1
                        }
                    }
                    Some(jj) => {
                        let w2 = st.map_w.is_bound(st.map_ev.find_key(jj));
                        if w1 == w2 {
                            if w1 {
                                2
                            } else {
                                1
                            }
                        } else {
                            3
                        }
                    }
                };

                // `cxx:547-679`: the curve the section edge starts from.
                let newc: Arc<dyn Curve>;
                let mut new_cs1: Option<Arc<dyn Curve2d>> = None;
                let mut new_cs2: Option<Arc<dyn Curve2d>> = None;
                let mut knownp: f64;
                let mut ebind: Option<TopoShape> = None;

                match choice {
                    1 => {
                        // `cxx:553-583`: the edge of `Res` that carries `vtx`
                        // away from the two edges of `MapEV`.
                        let res = match &st.res {
                            Some(r) => r.clone(),
                            None => return Ok(false),
                        };
                        let mut found: Option<TopoShape> = None;
                        for cand in explorer_shapes(&res, ShapeType::Edge) {
                            if is_same(&cand, &edg1) || is_same(&cand, &edg2) {
                                continue;
                            }
                            let forward = cand.oriented(Orientation::Forward);
                            if explorer_shapes(&forward, ShapeType::Vertex)
                                .iter()
                                .any(|v| is_same(&v, &vtx))
                            {
                                found = Some(cand);
                                break;
                            }
                        }
                        // `cxx:586-592`: the `else` is a `return`.
                        let Some(found) = found else {
                            return Ok(false);
                        };
                        // `cxx:579-583`.
                        let edge = Edge(found.clone());
                        let c = match BRepTool::edge_curve(&edge) {
                            Some(c) => c,
                            None => return Ok(false),
                        };
                        newc = Arc::from(c.transformed(&edge.0.location().transformation()));
                        knownp = BRepTool::parameter_of_vertex(&found, &vtx);
                        ebind = Some(found);
                    }
                    2 => {
                        // `cxx:595-601`.
                        let the_line = GpLin::new(st.normal_fg).translated_pnts(
                            st.normal_fg.location(),
                            &BRepTool::vertex_point(&Vertex(vtx.clone())),
                        );
                        newc = Arc::new(GeomLine::new(the_line));
                        knownp = 0.0;
                    }
                    _ => {
                        // `cxx:603-677`.
                        let (Some(f1), Some(f2)) = (
                            st.map_ev.find_from_key(&edg1).cloned(),
                            st.map_ev.find_from_key(&edg2).cloned(),
                        ) else {
                            return Ok(false);
                        };
                        let Some(s1) = BRepTool::face_surface(&Face(f1)) else {
                            return Ok(false);
                        };
                        let Some(s2) = BRepTool::face_surface(&Face(f2)) else {
                            return Ok(false);
                        };
                        let app_s1 = s1.gp_pln().is_none();
                        let app_s2 = s2.gp_pln().is_none();
                        st.i2s.perform(
                            s1.clone(),
                            s2.clone(),
                            CONFUSION,
                            true,
                            app_s1,
                            app_s2,
                        );
                        if !st.i2s.is_done() || st.i2s.nb_lines() <= 0 {
                            return Ok(false);
                        }
                        let pv = BRepTool::vertex_point(&Vertex(vtx.clone()));
                        let mut glob2min = f64::MAX;
                        let mut pmin = 0.0;
                        let mut imin: i32 = 0;
                        for i in 1..=st.i2s.nb_lines() {
                            let Some(line) = st.i2s.line(i).map(|l| l.curve.clone()) else {
                                continue;
                            };
                            let (first, last) = (line.first_parameter(), line.last_parameter());
                            let (sols, done) =
                                extrema_ext_pc_range_done(line.as_ref(), &pv, first, last);
                            if !done {
                                continue;
                            }
                            // `Extrema_ExtPC::TrimmedSquareDistances` (`:639`).
                            let d_first = pv.square_distance(&line.d0(first));
                            let d_last = pv.square_distance(&line.d0(last));
                            let mut dist2min = d_first;
                            let mut thepmin = first;
                            if d_last < dist2min {
                                thepmin = last;
                            }
                            for sol in &sols {
                                if sol.sq_dist < dist2min {
                                    dist2min = sol.sq_dist;
                                    thepmin = sol.u;
                                }
                            }
                            if dist2min < glob2min {
                                glob2min = dist2min;
                                pmin = thepmin;
                                imin = i;
                            }
                        }
                        if imin == 0 {
                            return Ok(false);
                        }
                        let Some(line) = st.i2s.line(imin).map(|l| l.curve.clone()) else {
                            return Ok(false);
                        };
                        newc = line;
                        knownp = pmin;
                        if app_s1 {
                            new_cs1 = st.i2s.line_on_s1(imin).cloned();
                        }
                        if app_s2 {
                            new_cs2 = st.i2s.line_on_s2(imin).cloned();
                        }
                    }
                }

                // `cxx:680-721`: pick the plane and, for Choice 1, the single
                // intersection and the vertex it collapses onto.
                let mut nbfois = 2;
                let mut vtx2: Option<TopoShape> = None;
                let mut p2 = 0.0;
                let mut is_left = false;
                let mut surface: Arc<dyn Surface> = st.new_sg.clone();
                if choice == 1 {
                    let Some(ancestors) = st.the_map_ef.find_from_key(&edg1) else {
                        return Ok(false);
                    };
                    let mut found = false;
                    for anc in ancestors {
                        if !contains(&f_faces, anc) {
                            continue;
                        }
                        if contains(&st.the_left, anc) {
                            surface = st.new_sg.clone();
                            is_left = true;
                        } else {
                            surface = st.new_sd.clone();
                            is_left = false;
                        }
                        nbfois = 1;
                        vtx2 = Some(vtx.clone());
                        p2 = knownp;
                        found = true;
                        break;
                    }
                    // `cxx:716-720`.
                    if !found {
                        return Ok(false);
                    }
                }

                // `cxx:723-760`.
                let mut vtx1: Option<TopoShape> = None;
                let mut p1 = 0.0;
                for it in 1..=nbfois {
                    if it == 2 {
                        surface = st.new_sd.clone();
                    }
                    let Some(intcs) = intcs_perform(newc.as_ref(), surface.as_ref()) else {
                        return Ok(false);
                    };
                    if intcs.nb_points() == 0 {
                        return Ok(false);
                    }
                    let mut imin: usize = 0;
                    let mut delta = (knownp - intcs.point(0).param()).abs();
                    for i in 1..intcs.nb_points() {
                        let newdelta = (knownp - intcs.point(i).param()).abs();
                        if newdelta < delta {
                            imin = i;
                            delta = newdelta;
                        }
                    }
                    let v = b.make_vertex(intcs.point(imin).pnt(), CONFUSION);
                    if it == 1 {
                        p1 = intcs.point(imin).param();
                        knownp = p1;
                        vtx1 = Some(v.0);
                    } else {
                        p2 = intcs.point(imin).param();
                        vtx2 = Some(v.0);
                    }
                }
                let (Some(vtx1), Some(vtx2)) = (vtx1, vtx2) else {
                    return Ok(false);
                };
                let vtx1 = Vertex(vtx1);
                let vtx2 = Vertex(vtx2);

                if (p1 - p2).abs() > Precision::PCONFUSION {
                    // `cxx:761-838`.
                    let first = newc.first_parameter();
                    let last = newc.last_parameter();
                    let mut new_edge = b.make_edge(newc.clone(), first, last);
                    registry.set_edge_tolerance(&new_edge.0, CONFUSION);
                    if p1 < p2 {
                        let mut a = vtx1.0.clone();
                        a.set_orientation(Orientation::Forward);
                        let mut zv = vtx2.0.clone();
                        zv.set_orientation(Orientation::Reversed);
                        b.add(&mut new_edge.0, &a);
                        b.add(&mut new_edge.0, &zv);
                    } else {
                        let mut a = vtx1.0.clone();
                        a.set_orientation(Orientation::Reversed);
                        let mut zv = vtx2.0.clone();
                        zv.set_orientation(Orientation::Forward);
                        b.add(&mut new_edge.0, &a);
                        b.add(&mut new_edge.0, &zv);
                    }
                    registry.set_vertex_parameter_on_edge(&new_edge.0, &vtx1.0, p1);
                    registry.set_vertex_parameter_on_edge(&new_edge.0, &vtx2.0, p2);

                    // `cxx:777-791`.
                    if let Some(cs1) = &new_cs1 {
                        if let Some(f1) = st.map_ev.find_from_key(&edg1) {
                            if let Some(s1) = BRepTool::face_surface(&Face(f1.clone())) {
                                registry.update_edge_pcurve_on_surface(
                                    &new_edge.0,
                                    &s1,
                                    Some(cs1.clone()),
                                    CONFUSION,
                                );
                            }
                        }
                    }
                    if let Some(cs2) = &new_cs2 {
                        if let Some(f2) = st.map_ev.find_from_key(&edg2) {
                            if let Some(s2) = BRepTool::face_surface(&Face(f2.clone())) {
                                registry.update_edge_pcurve_on_surface(
                                    &new_edge.0,
                                    &s2,
                                    Some(cs2.clone()),
                                    CONFUSION,
                                );
                            }
                        }
                    }

                    // `cxx:793`.
                    st.map_ev.add(&vtx, new_edge.0.clone());

                    if choice == 1 {
                        // `cxx:795-828`: the recovered edge with `vtx` swapped
                        // for `vtx1`.
                        let ebind = ebind.expect("Choice 1 always binds Ebind");
                        let mut ne = empty_copied_edge(&Edge(ebind.clone()));
                        for thevtx in explorer_shapes(&ebind, ShapeType::Vertex) {
                            if is_same(&thevtx, &vtx) {
                                let mut a = vtx1.0.clone();
                                a.set_orientation(thevtx.orientation());
                                b.add(&mut ne.0, &a);
                                registry.set_vertex_parameter_on_edge(&ne.0, &vtx1.0, p1);
                            } else {
                                b.add(&mut ne.0, &thevtx);
                                let theprm = BRepTool::parameter_of_vertex(&ebind, &thevtx);
                                registry.set_vertex_parameter_on_edge(&ne.0, &thevtx, theprm);
                                let tol = BRepTool::vertex_tolerance(&Vertex(thevtx.clone()));
                                registry.set_vertex_tolerance(&thevtx, tol);
                            }
                        }
                        st.map_son_s
                            .bind(&ebind, ne.0.oriented(Orientation::Forward));
                        if is_left {
                            bind_one(&mut st.map_sg, &vtx, &vtx1.0);
                        } else {
                            bind_one(&mut st.map_sd, &vtx, &vtx1.0);
                        }
                    } else {
                        // `cxx:829-837`.
                        bind_one(&mut st.map_sg, &vtx, &vtx1.0);
                        bind_one(&mut st.map_sd, &vtx, &vtx2.0);
                    }
                } else {
                    // `cxx:839-866`.
                    st.map_ev.add(&vtx, vtx2.0.clone());
                    if choice == 1 {
                        if is_left {
                            bind_one(&mut st.map_sg, &vtx, &vtx);
                        } else {
                            bind_one(&mut st.map_sd, &vtx, &vtx);
                        }
                    } else {
                        bind_one(&mut st.map_sg, &vtx, &vtx2.0);
                        bind_one(&mut st.map_sd, &vtx, &vtx2.0);
                    }
                }
            }
        }
        let _ = ();
        Ok(true)
    }
}
