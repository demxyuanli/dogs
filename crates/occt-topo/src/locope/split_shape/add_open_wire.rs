//! `LocOpe_SplitShape::AddOpenWire` (`LocOpe_SplitShape.cxx:781-1271`).
//!
//! An open splitting wire enters a face through one of its wires (`wfirst`) and
//! leaves through another (`wlast`):
//!
//! * when both ends land on the same wire, that wire is torn along the section:
//!   the section edges are appended in both directions and the walk builds two
//!   new boundary wires, one per resulting face;
//! * otherwise one new wire replaces `wfirst` / `wlast` (the section again in
//!   both directions) and the face is rebuilt in place.
//!
//! The reference prints two diagnostics through `std::cout` (`cxx:1043-1047`,
//! `cxx:1146-1149`); the port keeps the corresponding `return false` and the
//! "wire in neither new face" branch without printing.

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d;
use crate::brep_approx::surface_tool::SurfaceTool;
use crate::brep_tool::BRepTool;
use crate::brep_tools;
use crate::builder::TopoBuilder;
use crate::locope::shape_maps::{ShapeIndexedMap, ShapeSet};
use crate::locope::split_shape::direction::{
    check_overlapping, choix_uv, get_direction, same_uv,
};
use crate::locope::split_shape::is_inside::is_inside_wire;
use crate::locope::split_shape::{empty_copied_face, SplitShape};
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
use crate::topopebrep::build::tools2d::{first_vertex, last_vertex};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;
use crate::topo_tools::{children_keep_orientation, edge_vertices, explorer_shapes, is_same, wire_vertices};

impl SplitShape {
    /// `AddOpenWire(W, F)` (`cxx:781-1271`).
    pub fn add_open_wire(&mut self, w: &Wire, f: &Face) -> bool {
        // `cxx:783-787`.
        let lf = self.map_of(&f.0);

        // `cxx:789-797`: the two ends of the forward-oriented section.
        let (v_first, v_last) =
            wire_vertices(&Wire(w.0.oriented(Orientation::Forward)));
        let (Some(mut v_last), Some(v_first)) = (v_last, v_first) else {
            return false;
        };
        let tolf = BRepTool::vertex_tolerance(&v_first);
        let mut toll = BRepTool::vertex_tolerance(&v_last);
        let mut tol1 = tolf.max(toll);

        // `cxx:800-840`: the first descendant that contains the wire, plus the
        // two face wires its ends touch.
        let mut matched: Option<usize> = None;
        let mut wfirst: Option<Wire> = None;
        let mut wlast: Option<Wire> = None;
        for (i, s) in lf.iter().enumerate() {
            let Some(mut fac) = Face::wrap(s.clone()) else {
                continue;
            };
            if !is_inside_wire(&fac, w) {
                continue;
            }
            fac.0.set_orientation(Orientation::Forward);
            let mut ffound = false;
            let mut lfound = false;
            for sub in explorer_shapes(&fac.0, ShapeType::Wire) {
                let Some(wir) = Wire::wrap(sub.clone()) else {
                    continue;
                };
                for v in explorer_shapes(&sub, ShapeType::Vertex) {
                    if !ffound && is_same(&v, &v_first.0) {
                        ffound = true;
                        wfirst = Some(wir.clone());
                    } else if !lfound && is_same(&v, &v_last.0) {
                        lfound = true;
                        wlast = Some(wir.clone());
                    }
                    if ffound && lfound {
                        break;
                    }
                }
                if ffound && lfound {
                    break;
                }
            }
            if ffound && lfound {
                matched = Some(i);
                break;
            }
        }
        let Some(found) = matched else {
            return false;
        };
        let (Some(mut wfirst), Some(mut wlast)) = (wfirst, wlast) else {
            return false;
        };

        // `cxx:847-850`: `FaceRef` leaves the descendant list.
        let Some(mut face_ref) = Face::wrap(lf[found].clone()) else {
            return false;
        };
        face_ref.0.set_orientation(Orientation::Forward);
        let mut new_lf: Vec<TopoShape> = lf[..found].to_vec();
        new_lf.extend(lf[found + 1..].iter().cloned());

        let builder = TopoBuilder::new();

        // `cxx:851-857`: the surface adaptor, its periodicity and the
        // tolerance-based resolution the walk advances by.
        let bas = SurfaceAdaptor::of_face(&face_ref);
        let is_periodic = SurfaceTool::is_u_periodic(&bas) || SurfaceTool::is_v_periodic(&bas);
        tol1 = bas.u_resolution(tol1).max(bas.v_resolution(tol1));
        // `tol1` (`cxx:790-796`, `:857`, `:1037-1038`, `:1049`) is only ever
        // rewritten by the reference; nothing reads it.
        let _ = tol1;

        if is_same(&wfirst.0, &wlast.0) {
            // `cxx:860-1172`: two faces replace `FaceRef`.
            let mut new_w1 = Wire::new();
            new_w1.0.set_orientation(Orientation::Forward);
            let mut new_w2 = Wire::new();
            new_w2.0.set_orientation(Orientation::Forward);

            // `cxx:863-871`: the edges of `wfirst`, flagging the closed ones.
            let mut wires_first: Vec<TopoShape> = Vec::new();
            for sub in explorer_shapes(&wfirst.0, ShapeType::Edge) {
                let Some(e) = Edge::wrap(sub.clone()) else {
                    continue;
                };
                // `cxx:866-869`: `BRep_Tool::IsClosed(E, FaceRef)`.
                if BRepTool::is_closed_edge_face(&e, &face_ref) {
                    self.add_double_edge(&sub);
                }
                wires_first.push(sub);
            }

            // `cxx:880-889`: the section edges, each in both directions.
            let mut nb_e = 0usize;
            for sub in explorer_shapes(&w.0.oriented(Orientation::Forward), ShapeType::Edge) {
                nb_e += 1;
                let Some(e) = Edge::wrap(sub.clone()) else {
                    continue;
                };
                let orient = e.0.orientation();
                wires_first.push(sub.clone());
                wires_first.push(e.0.oriented(orient.reversed()));
                self.add_double_edge(&sub);
            }

            // `cxx:891-919`: the edge of `W` that starts at `Vfirst`.
            let mut last_edge: Option<Edge> = None;
            for sub in explorer_shapes(&w.0.oriented(Orientation::Forward), ShapeType::Edge) {
                let Some(e) = Edge::wrap(sub.clone()) else {
                    continue;
                };
                let mut hit = false;
                for v in explorer_shapes(&sub, ShapeType::Vertex) {
                    if is_same(&v, &v_first.0) {
                        hit = true;
                        break;
                    }
                }
                if hit {
                    last_edge = Some(e);
                    break;
                }
            }
            let Some(mut last_edge) = last_edge else {
                return false;
            };

            // `cxx:922-931`: the 2D point the walk starts from.
            let Some(face_on_wfirst) = Face::wrap(face_ref.0.oriented(wfirst.0.orientation()))
            else {
                return false;
            };
            let Some((c2d, f_par, l_par)) =
                boptools_2d::curve_on_surface_range(&last_edge, &face_on_wfirst)
            else {
                return false;
            };
            let pfirst = if last_edge.0.orientation() == Orientation::Forward {
                c2d.d0(f_par)
            } else {
                c2d.d0(l_par)
            };

            // `cxx:933-953`: the edge of `W` that ends at `Vlast`.
            for sub in explorer_shapes(&w.0.oriented(Orientation::Forward), ShapeType::Edge) {
                let Some(e) = Edge::wrap(sub.clone()) else {
                    continue;
                };
                if nb_e > 1 && is_same(&e.0, &last_edge.0) {
                    continue;
                }
                let mut hit = false;
                for v in explorer_shapes(&sub, ShapeType::Vertex) {
                    if is_same(&v, &v_last.0) {
                        hit = true;
                        break;
                    }
                }
                if hit {
                    last_edge = e;
                    break;
                }
            }

            // `cxx:955-956`: the outgoing direction at that end.
            let Some(face_on_wfirst) = Face::wrap(face_ref.0.oriented(wfirst.0.orientation()))
            else {
                return false;
            };
            let Some((mut plast, mut dlast)) = get_direction(&last_edge, &face_on_wfirst, false)
            else {
                return false;
            };

            // `cxx:958-967`.
            let cond = if is_periodic {
                !(is_same(&v_first.0, &v_last.0) && same_uv(&pfirst, &plast, &bas))
            } else {
                !is_same(&v_first.0, &v_last.0)
            };

            // `cxx:969-1049`: chain the remaining edges from `Vlast` on.
            let mut map_e: ShapeSet = ShapeSet::new();
            while cond {
                // `cxx:971-988`: candidates starting at `Vlast`.
                let mut poss_e: ShapeIndexedMap<TopoShape> = ShapeIndexedMap::new();
                for shape in &wires_first {
                    let Some(e) = Edge::wrap(shape.clone()) else {
                        continue;
                    };
                    let (vdeb, vfin) = edge_vertices(&e);
                    let starts_here = match e.0.orientation() {
                        Orientation::Forward => {
                            vdeb.map(|v| is_same(&v.0, &v_last.0)).unwrap_or(false)
                        }
                        Orientation::Reversed => {
                            vfin.map(|v| is_same(&v.0, &v_last.0)).unwrap_or(false)
                        }
                        _ => false,
                    };
                    if starts_here {
                        poss_e.add(&e.0, e.0.clone());
                    }
                }

                // `cxx:989-992`.
                let nb_poss = poss_e.extent();
                if nb_poss == 0 {
                    break;
                }

                // `cxx:994-1010`.
                let next_edge = if nb_poss == 1 {
                    let Some(e) = Edge::wrap(poss_e.find_key(1).clone()) else {
                        return false;
                    };
                    let Some(face) = Face::wrap(face_ref.0.oriented(wfirst.0.orientation()))
                    else {
                        return false;
                    };
                    let Some((p, d)) = get_direction(&e, &face, false) else {
                        return false;
                    };
                    plast = p;
                    dlast = d;
                    Some(e)
                } else {
                    // `cxx:1002-1010`: `ChoixUV` updates `plast` / `dlast`.
                    let Some(face) = Face::wrap(face_ref.0.oriented(wfirst.0.orientation()))
                    else {
                        return false;
                    };
                    let Some((e, _p, _d)) =
                        choix_uv(&last_edge, &face, &poss_e, &mut plast, &mut dlast)
                    else {
                        return false;
                    };
                    Some(e)
                };

                // `cxx:1012-1039`.
                let Some(next_edge) = next_edge else {
                    return false;
                };
                if map_e.contains(&next_edge.0) {
                    break;
                }
                builder.add(&mut new_w1.0, &next_edge.0);
                map_e.add(&next_edge.0);
                last_edge = next_edge;

                let next_vlast = if last_edge.0.orientation() == Orientation::Forward {
                    last_vertex(&last_edge.0)
                } else {
                    first_vertex(&last_edge.0)
                };
                let Some(next_vlast) = next_vlast.and_then(Vertex::wrap) else {
                    return false;
                };
                v_last = next_vlast;
                toll = BRepTool::vertex_tolerance(&v_last);
                tol1 = tolf.max(toll);

                // `cxx:1049`.
                tol1 = bas.u_resolution(tol1).max(bas.v_resolution(tol1));
                // `tol1` is written but never read by the reference.
                let _ = tol1;
            }

            // `cxx:1051-1089`: the rest of the boundary for the second face.
            let mut nb_add_bound = 0usize;
            let mut an_e1: Option<TopoShape> = None;
            let mut an_e2: Option<TopoShape> = None;
            for shape in &wires_first {
                if !map_e.contains(shape) {
                    builder.add(&mut new_w2.0, shape);
                    map_e.add(shape);
                    nb_add_bound += 1;
                    if an_e1.is_none() {
                        an_e1 = Some(shape.clone());
                    } else {
                        an_e2 = Some(shape.clone());
                    }
                }
            }
            // `cxx:1072-1082`: a face bounded by fewer than two edges, or by two
            // overlapping ones, is degenerate.
            if nb_add_bound < 2 {
                return false;
            }
            if nb_add_bound == 2 {
                if let (Some(e1), Some(e2)) = (an_e1, an_e2) {
                    let (Some(e1), Some(e2)) = (Edge::wrap(e1), Edge::wrap(e2)) else {
                        return false;
                    };
                    if check_overlapping(&e1, &e2, &face_ref) {
                        return false;
                    }
                }
            }

            // `cxx:1085-1118`: the same checks on `newW1`.
            let mut nb_add_bound = 0usize;
            let mut an_e11: Option<TopoShape> = None;
            let mut an_e12: Option<TopoShape> = None;
            for child in children_keep_orientation(&new_w1.0) {
                if child.shape_type() != ShapeType::Edge {
                    continue;
                }
                nb_add_bound += 1;
                if an_e11.is_none() {
                    an_e11 = Some(child);
                } else {
                    an_e12 = Some(child);
                }
            }
            if nb_add_bound < 2 {
                return false;
            }
            if nb_add_bound == 2 {
                if let (Some(e1), Some(e2)) = (an_e11, an_e12) {
                    let (Some(e1), Some(e2)) = (Edge::wrap(e1), Edge::wrap(e2)) else {
                        return false;
                    };
                    if check_overlapping(&e1, &e2, &face_ref) {
                        return false;
                    }
                }
            }

            // `cxx:1119-1129`.
            let mut new_f1 = empty_copied_face(&face_ref);
            new_f1.0.set_orientation(Orientation::Forward);
            let mut new_f2 = empty_copied_face(&face_ref);
            new_f2.0.set_orientation(Orientation::Forward);
            builder.add(&mut new_f1.0, &new_w1.0);
            builder.add(&mut new_f2.0, &new_w2.0);

            // `cxx:1130-1150`: every other wire of the face joins the piece it
            // lies in.
            for sub in explorer_shapes(&face_ref.0.oriented(Orientation::Forward), ShapeType::Wire) {
                if is_same(&sub, &wfirst.0) {
                    continue;
                }
                let Some(wir) = Wire::wrap(sub.clone()) else {
                    continue;
                };
                if is_inside_wire(&new_f1, &wir) {
                    builder.add(&mut new_f1.0, &sub);
                } else if is_inside_wire(&new_f2, &wir) {
                    builder.add(&mut new_f2.0, &sub);
                }
                // `cxx:1144-1149`: the reference prints a warning here.
            }

            // `cxx:1151-1152`.
            new_lf.push(new_f1.0);
            new_lf.push(new_f2.0);

            // `cxx:1154-1171`: `wfirst` gives way to the two new wires.
            for sub in explorer_shapes(&f.0, ShapeType::Wire) {
                if !self.map_bound(&sub) {
                    return false;
                }
                let ls = self.map_of(&sub);
                if let Some(pos) = ls.iter().position(|s| is_same(s, &wfirst.0)) {
                    let mut ls = ls;
                    ls.remove(pos);
                    ls.push(new_w1.0.clone());
                    ls.push(new_w2.0.clone());
                    self.set_map_of(&sub, ls);
                }
            }
        } else {
            // `cxx:1173-1267`: one face replaces `FaceRef`.
            let Some(outer_w) = brep_tools::outer_wire(&face_ref) else {
                return false;
            };
            let mut new_wire = Wire::new();
            new_wire.0.set_orientation(Orientation::Forward);

            // `cxx:1184-1192`.
            let or_relat = if wfirst.0.orientation() == wlast.0.orientation() {
                Orientation::Forward
            } else {
                Orientation::Reversed
            };

            // `cxx:1194-1198`.
            if is_same(&wlast.0, &outer_w.0) {
                wlast = wfirst.clone();
                wfirst = outer_w;
            }

            // `cxx:1200-1206`: the edges of `wfirst`, kept as they are.
            for sub in explorer_shapes(&wfirst.0.oriented(Orientation::Forward), ShapeType::Edge) {
                builder.add(&mut new_wire.0, &sub);
            }

            // `cxx:1208-1213`: the edges of `wlast`, related through `orRelat`.
            for sub in explorer_shapes(&wlast.0.oriented(Orientation::Forward), ShapeType::Edge) {
                let Some(e) = Edge::wrap(sub.clone()) else {
                    continue;
                };
                let orient = Orientation::compose(or_relat, e.0.orientation());
                builder.add(&mut new_wire.0, &e.0.oriented(orient));
            }

            // `cxx:1215-1223`: the section, forwards and backwards.
            for sub in explorer_shapes(&w.0.oriented(Orientation::Forward), ShapeType::Edge) {
                let Some(e) = Edge::wrap(sub.clone()) else {
                    continue;
                };
                let orient = e.0.orientation();
                builder.add(&mut new_wire.0, &e.0.oriented(orient));
                builder.add(&mut new_wire.0, &e.0.oriented(orient.reversed()));
                self.add_double_edge(&e.0.oriented(orient));
            }

            // `cxx:1227-1229`.
            let mut new_face = empty_copied_face(&face_ref);
            new_face.0.set_orientation(Orientation::Forward);
            face_ref.0.set_orientation(Orientation::Forward);

            // `cxx:1230-1241`: `wfirst` is replaced by the rebuilt wire, `wlast`
            // is dropped, the rest is kept.
            for sub in explorer_shapes(&face_ref.0.oriented(Orientation::Forward), ShapeType::Wire) {
                let Some(wir) = Wire::wrap(sub.clone()) else {
                    continue;
                };
                if is_same(&sub, &wfirst.0) {
                    builder.add(&mut new_face.0, &new_wire.0.oriented(wir.0.orientation()));
                } else if !is_same(&sub, &wlast.0) {
                    builder.add(&mut new_face.0, &sub);
                }
            }

            // `cxx:1242`.
            new_lf.push(new_face.0);

            // `cxx:1244-1264`: the descendants of `wfirst` / `wlast` become the
            // rebuilt wire.
            for sub in explorer_shapes(&f.0, ShapeType::Wire) {
                if !self.map_bound(&sub) {
                    return false;
                }
                let ls = self.map_of(&sub);
                let mut kept: Vec<TopoShape> = Vec::new();
                let mut touch = false;
                for s in ls {
                    if is_same(&s, &wfirst.0) || is_same(&s, &wlast.0) {
                        touch = true;
                    } else {
                        kept.push(s);
                    }
                }
                if touch {
                    kept.push(new_wire.0.clone());
                    self.set_map_of(&sub, kept);
                }
            }
        }

        // `cxx:575` / `cxx:1242`: the rebuilt faces join the descendant list.
        self.map_bind(&f.0, new_lf);
        true
    }
}
