//! `ShapeFix_Face::FixWiresTwoCoincEdges` (`ShapeFix_Face.cxx:2829-2900`).
//!
//! Drops every boundary wire that holds exactly two coincident edges - the
//! pair a `ShapeExtend_WireData` exposes when one seam edge is stored twice -
//! and rebuilds the face without it. `ShapeFix_Face::Perform` runs it
//! (`cxx:676-679`) right before `FixIntersectingWires` (`cxx:680`) and
//! `FixOrientation` (`cxx:692`).
//!
//! The method lives beside `shape_fix_face.rs` because it is only a few dozen
//! lines and that file is already past the repository's split threshold; the
//! `impl ShapeFixFace` block keeps the public path `ShapeFixFace::<method>`
//! unchanged.

use std::collections::HashSet;
use std::sync::Arc;

use occt_core::bnd::BndBox2d;
use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::Surface;

use crate::abs::{Orientation, ShapeType};
use crate::bnd_lib_add2d::{add_adaptor, add_adaptor_range};
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::fclass2d::{FaceState, FClass2d};
use crate::pcurve_full::surface_value_of_uv;
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
use crate::shape_fix_compose_shell::{reverse_wire_data_on_face, ReShape};
use crate::shhealing::adjust_by_period;
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{edges_of_wire, is_same};

use super::shape_fix_face::wire_from_wire_data;
use super::ShapeFixFace;

impl ShapeFixFace {
    /// `ShapeFix_Face::FixWiresTwoCoincEdges()` (`ShapeFix_Face.cxx:2829-2900`):
    /// returns `isFixed`. On success `myFace` becomes the rebuilt face and the
    /// context records `Replace(old_face, new_face)` (`cxx:2894`).
    pub fn fix_wires_two_coinc_edges(&mut self) -> bool {
        // `cxx:2831-2835`: `myFace = TopoDS::Face(Context()->Apply(myFace))`.
        let mut face = match self.face.clone() {
            Some(f) => f,
            None => return false,
        };
        {
            let applied = self.context.apply(&face.0);
            if applied.is_face() {
                face = Face(applied);
                self.face = Some(face.clone());
            }
        }

        // `cxx:2837`: `TopAbs_Orientation ori = myFace.Orientation()`.
        let ori = face.0.orientation();

        // `cxx:2841-2856`: `TopoDS_Iterator(myFace, false)`, oriented wires only.
        let children: Vec<TopoShape> = {
            let ts = face.0.tshape.read().expect("poisoned TShape lock");
            ts.children.clone()
        };
        let mut nb_wires = 0usize;
        for child in &children {
            if child.shape_type() != ShapeType::Wire {
                continue;
            }
            let o = child.orientation();
            if o == Orientation::Forward || o == Orientation::Reversed {
                nb_wires += 1;
            }
        }
        // `cxx:2857-2860`.
        if nb_wires < 2 {
            return false;
        }

        // `cxx:2843-2846`: `emptyCopied` / `Orientation(TopAbs_FORWARD)`.
        let builder = TopoBuilder::new();
        let mut new_face = match BRepTool::face_surface(&face) {
            Some(surf) => builder.make_face(surf, &[]),
            None => Face::new(),
        };
        // The location is applied before the children so `TopoDS_Builder::Add`
        // compensates through `myFace`'s location exactly like `EmptyCopied`
        // plus `B.Add` do.
        new_face.0.set_orientation(Orientation::Forward);
        new_face.0.set_location(face.0.location());
        if let Some(geom) = GeometryRegistry::global().face_geom(&face.0) {
            GeometryRegistry::global().set_face(&new_face.0, geom);
        }

        // `cxx:2861-2888`.
        let mut is_fixed = false;
        for child in &children {
            // `cxx:2862-2867`: every non-wire or non-oriented child is copied
            // unchanged, in iteration order.
            if child.shape_type() != ShapeType::Wire {
                builder.add(&mut new_face.0, child);
                continue;
            }
            let o = child.orientation();
            if o != Orientation::Forward && o != Orientation::Reversed {
                builder.add(&mut new_face.0, child);
                continue;
            }
            let wire = Wire(child.clone());
            let ents = edges_of_wire(&wire);
            if ents.len() == 2 {
                // `cxx:2872-2878`.
                let mut e1 = ents[0].clone();
                e1.0.set_orientation(Orientation::Forward);
                let mut e2 = ents[1].clone();
                e2.0.set_orientation(Orientation::Forward);
                // `E1 == E2` is `TopoDS_Shape::operator==`, i.e. `IsEqual`:
                // same `TShape`, same location, same orientation. Both edges are
                // read off one `ShapeExtend_WireData`, so they carry the wire's
                // composed location and both were forced FORWARD just above -
                // the `TShape` handle is then the whole test.
                if !is_same(&e1.0, &e2.0) {
                    builder.add(&mut new_face.0, child);
                } else {
                    is_fixed = true;
                }
            } else {
                // `cxx:2883-2886`.
                builder.add(&mut new_face.0, child);
            }
        }

        // `cxx:2890-2897`.
        if is_fixed {
            new_face.0.set_orientation(ori);
            self.context.replace(&face.0, &new_face.0);
            self.face = Some(new_face);
        }
        is_fixed
    }
}

/// Puts `face`'s children back after a single-wire view. The view shares the
/// face `TShape`, so `CurveOnSurface` still finds pcurves stored by face key
/// (`BRep_Tool` keys them by surface).
struct FaceChildrenGuard<'a> {
    shape: &'a TopoShape,
    saved: Vec<TopoShape>,
}

impl Drop for FaceChildrenGuard<'_> {
    fn drop(&mut self) {
        if let Ok(mut ts) = self.shape.tshape.write() {
            ts.children = std::mem::take(&mut self.saved);
        }
    }
}

fn shape_ptr(s: &TopoShape) -> usize {
    Arc::as_ptr(&s.tshape) as usize
}

fn same_oriented(a: &TopoShape, b: &TopoShape) -> bool {
    is_same(a, b) && a.orientation() == b.orientation()
}

fn reversed_wire(wire: &Wire, face: &Face) -> Wire {
    let (mut data, non_manifold): (Vec<Edge>, Vec<Edge>) =
        edges_of_wire(wire).into_iter().partition(|e| {
            let o = e.0.orientation();
            o == Orientation::Forward || o == Orientation::Reversed
        });
    if wire.0.orientation().is_reversed() {
        data.reverse();
    }
    reverse_wire_data_on_face(&mut data, face);
    wire_from_wire_data(&data, &non_manifold)
}

fn wire_uv_box(wire: &Wire, face: &Face) -> BndBox2d {
    let mut a_box = BndBox2d::new();
    for edge in edges_of_wire(wire) {
        let Some((cw, cf, cl)) = curve_on_surface_range(&edge, face) else {
            continue;
        };
        let a_first = cw.first_parameter();
        let a_last = cw.last_parameter();
        // cxx:1324-1331.
        if cw.is_bspline2d() && (cf < a_first || cl > a_last) {
            add_adaptor(cw.as_ref(), CONFUSION, &mut a_box);
        } else {
            add_adaptor_range(cw.as_ref(), cf, cl, CONFUSION, &mut a_box);
        }
    }
    a_box
}

/// `FClass2d` of `wire` alone on a FORWARD copy of `face` (`cxx:1364-1371`).
fn class_one_wire(face: &Face, wire: &Wire) -> Result<FClass2d, String> {
    let mut view = face.clone();
    view.0.set_orientation(Orientation::Forward);
    let saved = view
        .0
        .tshape
        .read()
        .expect("poisoned TShape lock")
        .children
        .clone();
    {
        let mut ts = view.0.tshape.write().expect("poisoned TShape lock");
        ts.children = vec![wire.0.clone()];
    }
    let _guard = FaceChildrenGuard {
        shape: &view.0,
        saved,
    };
    FClass2d::new(&view, PCONFUSION)
}

/// `FixOrientation` when the face has several wires (`ShapeFix_Face.cxx:1276-1606`).
/// Returns how many wires were reversed.
pub(super) fn orient_several_wires(
    face: &Face,
    surf: &dyn Surface,
    ws: &mut [Wire],
    all_sub: &[TopoShape],
    map_wires: &mut Vec<(Wire, Vec<Wire>)>,
) -> usize {
    let nb = ws.len();
    if nb < 2 {
        return 0;
    }
    let uclosed = surf.is_u_closed();
    let vclosed = surf.is_v_closed();
    let (suf, sul) = surf.u_range();
    let (svf, svl) = surf.v_range();
    let u_range = sul - suf;
    let v_range = svl - svf;

    let mut boxes = Vec::with_capacity(nb);
    let mut u_middle = 0.0;
    let mut v_middle = 0.0;
    let mut is_first = true;
    for wire in ws.iter() {
        let mut a_box = wire_uv_box(wire, face);
        if let Some((xmin, ymin, xmax, ymax)) = a_box.get() {
            if is_first {
                is_first = false;
                u_middle = (xmin + xmax) * 0.5;
                v_middle = (ymin + ymax) * 0.5;
            } else {
                let x_shift = if surf.is_u_closed() {
                    adjust_by_period(0.5 * (xmin + xmax), u_middle, u_range)
                } else {
                    0.0
                };
                let y_shift = if surf.is_v_closed() {
                    adjust_by_period(0.5 * (ymin + ymax), v_middle, v_range)
                } else {
                    0.0
                };
                a_box.update(xmin + x_shift, ymin + y_shift, xmax + x_shift, ymax + y_shift);
            }
        }
        boxes.push(a_box);
    }

    let mut si = vec![0i32; nb];
    let mut inners: Vec<Vec<Wire>> = vec![Vec::new(); nb];
    let mut map_int: HashSet<usize> = HashSet::new();
    let mut nrev = 0usize;

    for i in 0..nb {
        let aw = ws[i].clone();
        let a_box1 = boxes[i];
        let Ok(clas) = class_one_wire(face, &aw) else {
            continue;
        };
        let mut check_shift = true;
        let staout = clas.perform_infinite_point_tab_orien();
        let mut sta = FaceState::Out;
        let mut int_wires: Vec<Wire> = Vec::new();
        let mut a_wire_it = 0i32;
        let mut stop_j = false;

        for a_sh2 in all_sub {
            a_wire_it += 1;
            if same_oriented(&aw.0, a_sh2) {
                continue;
            }
            let mut stb = FaceState::Unknown;
            if a_sh2.shape_type() == ShapeType::Vertex {
                a_wire_it -= 1;
                let p = BRepTool::vertex_point(&Vertex(a_sh2.clone()));
                let p2d = surface_value_of_uv(surf, &p, CONFUSION);
                stb = clas.perform_tab_orien(p2d);
                if stb == staout && (uclosed || vclosed) {
                    if uclosed {
                        stb = clas.perform_tab_orien(GpPnt2d::new(p2d.x() + u_range, p2d.y()));
                    }
                    if stb == staout && vclosed {
                        stb = clas.perform_tab_orien(GpPnt2d::new(p2d.x(), p2d.y() + v_range));
                    }
                }
            } else if a_sh2.shape_type() == ShapeType::Wire {
                check_shift = true;
                let idx = (a_wire_it - 1) as usize;
                if idx >= boxes.len() || boxes[idx].is_out_box(&a_box1) {
                    continue;
                }
                let bw = Wire(a_sh2.clone());
                for ed in edges_of_wire(&bw) {
                    let Some((cw, cf, cl)) = curve_on_surface_range(&ed, face) else {
                        continue;
                    };
                    let unp = cw.d0((cf + cl) / 2.0);
                    let ste = clas.perform_tab_orien(unp);
                    if ste == FaceState::Out || ste == FaceState::In {
                        if stb == FaceState::Unknown {
                            stb = ste;
                        } else if stb != ste {
                            sta = FaceState::Unknown;
                            si[i] = 0;
                            stop_j = true;
                            break;
                        }
                    }
                    if stb == staout && check_shift && (uclosed || vclosed) {
                        check_shift = false;
                        let mut found = false;
                        let mut unp1 = unp;
                        if uclosed {
                            unp1 = GpPnt2d::new(unp.x() + u_range, unp.y());
                            found = staout != clas.perform_tab_orien(unp1);
                            if !found {
                                unp1 = GpPnt2d::new(unp.x() - u_range, unp.y());
                                found = staout != clas.perform_tab_orien(unp1);
                            }
                        }
                        if vclosed && !found {
                            unp1 = GpPnt2d::new(unp.x(), unp.y() + v_range);
                            found = staout != clas.perform_tab_orien(unp1);
                            if !found {
                                unp1 = GpPnt2d::new(unp.x(), unp.y() - v_range);
                                found = staout != clas.perform_tab_orien(unp1);
                            }
                        }
                        if !found && uclosed && vclosed {
                            'diag: for dx in [-1.0_f64, 1.0] {
                                for dy in [-1.0_f64, 1.0] {
                                    unp1 = GpPnt2d::new(
                                        unp.x() + u_range * dx,
                                        unp.y() + v_range * dy,
                                    );
                                    found = staout != clas.perform_tab_orien(unp1);
                                    if found {
                                        break 'diag;
                                    }
                                }
                            }
                        }
                        if found {
                            stb = if stb == FaceState::In {
                                FaceState::Out
                            } else {
                                FaceState::In
                            };
                            let vec = GpVec2d::new(unp1.x() - unp.x(), unp1.y() - unp.y());
                            super::face_geom_helpers::shift_2d_wire(&bw, face, &vec, surf, false);
                        }
                    }
                }
            }
            if stb == staout {
                sta = FaceState::In;
            } else if a_sh2.shape_type() == ShapeType::Wire {
                int_wires.push(Wire(a_sh2.clone()));
                map_int.insert(shape_ptr(a_sh2));
            } else {
                map_int.insert(shape_ptr(a_sh2));
            }
            if stop_j {
                break;
            }
        }

        if sta == FaceState::Unknown {
            continue;
        }
        inners[i] = int_wires.clone();
        if sta == FaceState::Out {
            if staout == FaceState::In {
                ws[i] = reversed_wire(&aw, face);
                nrev += 1;
                si[i] = 1;
                map_wires.push((ws[i].clone(), int_wires));
            } else {
                si[i] = 1;
                map_wires.push((aw, int_wires));
            }
        } else if staout == FaceState::Out {
            si[i] = 2;
        } else {
            si[i] = 3;
        }
    }

    for i in 0..nb {
        let tmpi = si[i];
        if tmpi <= 1 {
            continue;
        }
        let key = shape_ptr(&ws[i].0);
        if !map_int.contains(&key) {
            if tmpi == 3 {
                let iw = inners[i].clone();
                ws[i] = reversed_wire(&ws[i], face);
                nrev += 1;
                map_wires.push((ws[i].clone(), iw));
            } else {
                map_wires.push((ws[i].clone(), inners[i].clone()));
            }
        } else if tmpi == 2 {
            ws[i] = reversed_wire(&ws[i], face);
            nrev += 1;
        }
    }
    nrev
}

