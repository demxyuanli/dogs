//! `TopOpeBRepBuild_CorrectFace2d` — the copy/connect half.
//!
//! Sources (OCCT 8.0.0, `TopOpeBRepBuild_CorrectFace2d.cxx`):
//! `BuildCopyData` `:116-222`, `ConnectWire` `:340-585`, `GetP2dFL` `:687-702`,
//! `MakeHeadList` `:706-744`, `CheckList` `:748-847`.
//!
//! `BuildCopyData` deep-copies the face by `EmptyCopied` on face / wire / edge
//! (`BRep_TFace::EmptyCopy`, `BRep_TEdge::EmptyCopy` copy the geometry
//! representations), remaps the avoid map onto the copies and, in the
//! `aNeedToUsePMap` arm, reuses the edges already stored in `myMap`.

use occt_core::gp::GpPnt2d;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{explorer_shapes, is_equal, vertex_position};
use crate::topopebrep::build::tools::{write_edge_pcurve, write_edge_pcurves};
use crate::topopebrep::ds::data_structure::{IndexedShapeData, IndexedShapeSet};

use super::CorrectFace2d;

const TWO_PI: f64 = 2.0 * std::f64::consts::PI;

/// `aFace.EmptyCopied()` for a `TopoDS_Face`: `BRep_TFace::EmptyCopy`
/// (`BRep_TFace.cxx:36-43`) carries the surface, the location and the
/// tolerance over.
fn empty_copied_face(f: &Face) -> Option<Face> {
    let fou = crate::topopebrep::ds::build_tool::BuildTool::new().copy_face(&f.0);
    Face::wrap(fou)
}

/// `E.EmptyCopied()` for a `TopoDS_Edge`: `BRep_TEdge::EmptyCopy`
/// (`BRep_TEdge.cxx:105-133`) copies the tolerance, every curve
/// representation and the `SameParameter` / `SameRange` / `Degenerated` flags.
///
/// `BuildTool::copy_edge` already re-registers the geometry side-table for the
/// new `TShape`; the pcurve slots are copied unconditionally because
/// `copy_edge` skips them on a pcurve-only edge (no `BRep_Curve3D`).
fn empty_copied_edge(e: &TopoShape) -> Option<TopoShape> {
    let eou = crate::topopebrep::ds::build_tool::BuildTool::new().copy_edge(e);
    let reg = GeometryRegistry::global();
    reg.copy_edge_pcurve_slots(&eou, e);
    let tol = BRepTool::edge_tolerance(&Edge::wrap(e.clone())?);
    reg.raise_edge_tolerance(&eou, tol);
    Some(eou)
}

impl CorrectFace2d {
    /// `BuildCopyData(aFace, anAvoidMap, aCopyFace, aCopyAvoidMap,
    /// aNeedToUsePMap)` (`cxx:116-222`).
    ///
    /// `aFace` / `anAvoidMap` are taken by value-shaped references because the
    /// OCCT callers pass members of the same object; only reads happen on them.
    pub(super) fn build_copy_data(
        &mut self,
        a_face: &Face,
        an_avoid_map: &IndexedShapeSet,
        a_copy_face: &mut Face,
        a_copy_avoid_map: &mut IndexedShapeSet,
        a_need_to_use_p_map: bool,
    ) {
        let builder = TopoBuilder::new();
        let mut ed_map: IndexedShapeData<TopoShape> = IndexedShapeData::new();

        // 1. Copy `aFace` =>> `aCopyFace` (`cxx:129-190`).
        let Some(copy_face) = empty_copied_face(a_face) else {
            return;
        };
        *a_copy_face = copy_face;

        for w in explorer_shapes(&a_face.0, ShapeType::Wire) {
            let mut a_wire = w.empty_copied();
            for e in explorer_shapes(&w, ShapeType::Edge) {
                if let Some(mut an_edge) = ed_map.find_from_key(&e).cloned() {
                    an_edge.set_orientation(e.orientation());
                    builder.add(&mut a_wire, &an_edge);
                    continue;
                }

                if a_need_to_use_p_map {
                    let mut hit: Option<TopoShape> = None;
                    if let Some(map) = &self.my_map {
                        let mut m = map.borrow_mut();
                        if let Some(v) = m.change_from_key(&e) {
                            v.set_orientation(e.orientation());
                            hit = Some(v.clone());
                        }
                    }
                    if let Some(an_edge) = hit {
                        ed_map.add(&e, an_edge.clone());
                        builder.add(&mut a_wire, &an_edge);
                        continue;
                    }
                }

                // Add the edge (`cxx:168-183`).
                let Some(mut an_edge) = empty_copied_edge(&e) else {
                    continue;
                };
                ed_map.add(&e, an_edge.clone());
                for v in explorer_shapes(&e, ShapeType::Vertex) {
                    builder.add(&mut an_edge, &v);
                }
                an_edge.set_orientation(e.orientation());
                builder.add(&mut a_wire, &an_edge);
            }
            // Add the wire (`cxx:184-190`).
            a_wire.set_orientation(w.orientation());
            ed_map.add(&w, a_wire.clone());
            builder.add(&mut a_copy_face.0, &a_wire);
        }

        // 2. Copy `anAvoidMap` =>> `aCopyAvoidMap` (`cxx:192-207`).
        let a_nb = an_avoid_map.extent();
        for i in 1..=a_nb {
            let Some(a_sh) = an_avoid_map.find_key(i).cloned() else {
                continue;
            };
            if let Some(mut a_copy_shape) = ed_map.find_from_key(&a_sh).cloned() {
                a_copy_shape.set_orientation(a_sh.orientation());
                a_copy_avoid_map.add(&a_copy_shape);
            }
        }

        // 3. Inversed `EdMap` (`cxx:209-221`).
        if a_need_to_use_p_map {
            self.my_ed_map_inversed.clear();
            let a_nb = ed_map.extent();
            for i in 1..=a_nb {
                let (Some(a_sh), Some(a_sh_copy)) = (
                    ed_map.find_key(i).cloned(),
                    ed_map.find_from_index(i).cloned(),
                ) else {
                    continue;
                };
                self.my_ed_map_inversed.add(&a_sh_copy, a_sh);
            }
        }
    }

    /// `ConnectWire(aCopyFace, aCopyAvoidMap, aTryFlag)` (`cxx:340-585`).
    ///
    /// Builds the connection order of `myCurrentWire` (through
    /// `myOrderedWireList`), picks the first movable edge and translates the
    /// p-curves of every following edge so that its first 2d point lands on the
    /// previous edge's last 2d point. `aTryFlag` selects `BRep_Builder::UpdateEdge`
    /// (a throw-away copy) or the member `UpdateEdge` (which also writes the
    /// original face's edges through `myMap`).
    pub(super) fn connect_wire(
        &mut self,
        a_copy_face: &mut Face,
        a_copy_avoid_map: &IndexedShapeSet,
        a_try_flag: bool,
    ) -> i32 {
        // 1. Make right order (`cxx:347-353`).
        let priz = self.make_right_wire();
        if priz != 0 {
            return priz;
        }

        // 2. Define the first edge on the wire from `aCopyAvoidMap`
        //    (`cxx:355-381`).
        let mut a_nb_edges = 0;
        let mut an_avoid_map = IndexedShapeSet::new();
        for an_edge in explorer_shapes(&self.my_current_wire, ShapeType::Edge) {
            if a_copy_avoid_map.contains(&an_edge) {
                an_avoid_map.add(&an_edge);
            }
            a_nb_edges += 1;
        }

        if a_nb_edges == 1 {
            // Nothing to do with this wire.
            return 0;
        }

        let a_nb_avoid_edges_on_wire = an_avoid_map.extent();
        if a_nb_avoid_edges_on_wire == a_nb_edges {
            // All edges of the wire are good.
            return 0;
        }

        // If at least one non-degenerated edge exists among `anAvoidMap`, take
        // it as `aFirstEdge` (`cxx:384-398`).
        let mut a_first_edge: Option<TopoShape> = None;
        for i in 1..=a_nb_avoid_edges_on_wire {
            let Some(an_edge) = an_avoid_map.find_key(i).cloned() else {
                continue;
            };
            let Some(edge) = Edge::wrap(an_edge.clone()) else {
                continue;
            };
            if !BRepTool::is_degenerated(&edge) {
                a_first_edge = Some(an_edge);
                break;
            }
        }

        if a_first_edge.is_none() {
            // Every edge of `anAvoidMap` is degenerated: take the one with the
            // largest distance between its first and last 2d points
            // (`cxx:398-418`).
            let mut a_d_max = 0.0;
            for i in 1..=a_nb_avoid_edges_on_wire {
                let Some(an_edge) = an_avoid_map.find_key(i).cloned() else {
                    continue;
                };
                let Some(edge) = Edge::wrap(an_edge.clone()) else {
                    continue;
                };
                let Some((pf, pl)) = Self::get_p2d_fl(a_copy_face, &edge) else {
                    continue;
                };
                let a_r = pf.distance(&pl);
                if (a_r - a_d_max) > 1e-7 {
                    a_d_max = a_r;
                    a_first_edge = Some(an_edge);
                }
            }
        }

        // 3. Build `HeadList` from `myOrderedWireList`, `aFirstEdge` first
        //    (`cxx:420-423`).
        let mut head_list: Vec<TopoShape> = Vec::new();
        self.make_head_list(a_first_edge.as_ref(), &mut head_list);

        // 5. Connect the edges (`cxx:430-...`).
        let mut an_ed_map = IndexedShapeSet::new();
        let a_d_tolerance = occt_core::precision::CONFUSION;
        let a_seq: Vec<TopoShape> = head_list;
        let a_nb_edges = a_seq.len() as i32;

        let Some(a_fst) = a_seq.first().cloned() else {
            // `aSeqEdges(1)` on an empty sequence: OCCT reads out of range.
            return 1;
        };
        let Some(a_fst_edge) = Edge::wrap(a_fst) else {
            return 1;
        };
        let Some((_pa, mut pb)) = Self::get_p2d_fl(a_copy_face, &a_fst_edge) else {
            return 1;
        };
        let pa1 = _pa;

        for i in 2..=a_nb_edges {
            let Some(an_edge_shape) = a_seq.get((i - 1) as usize).cloned() else {
                continue;
            };
            let Some(an_edge) = Edge::wrap(an_edge_shape.clone()) else {
                continue;
            };
            let Some((pf, pl)) = Self::get_p2d_fl(a_copy_face, &an_edge) else {
                continue;
            };

            let a_d = pf.distance(&pb);
            if a_d < a_d_tolerance {
                // Already connected; `PA` is not read again in OCCT either
                // (`cxx:461-465` only refreshes `PB`).
                pb = pl;
                continue;
            }

            // Translation vector (`cxx:468-469`).
            let a_tr_v = occt_core::gp::GpVec2d::new(pf.x() - pb.x(), pf.y() - pb.y());

            // The translation vector must be equal to `2*PI*n`, otherwise the
            // wire is wrong (NIZHNY-MZV, `cxx:470-501`).
            let u = a_tr_v.x();
            let v = a_tr_v.y();
            let (up, vp) = match BRepTool::face_surface(a_copy_face) {
                Some(s) => (s.is_u_periodic(), s.is_v_periodic()),
                None => (false, false),
            };
            let mut non_pu = u.abs() < 1e-7;
            let mut non_pv = v.abs() < 1e-7;
            if !non_pu && up {
                let d_u = (u.abs() % TWO_PI).abs();
                non_pu = d_u > 1e-7 && (TWO_PI - d_u > 1e-7);
            }
            if !non_pv && vp {
                let d_v = (v.abs() % TWO_PI).abs();
                non_pv = d_v > 1e-7 && (TWO_PI - d_v > 1e-7);
            }
            if non_pu && non_pv && !BRepTool::is_degenerated(&an_edge) {
                return 1;
            }

            let tol = self.my_face_tolerance;

            if BRepTool::is_closed_edge_face(&an_edge, a_copy_face) {
                // a. Closed edge (`cxx:502-546`).
                if an_ed_map.contains(&an_edge_shape) {
                    continue;
                }
                an_ed_map.add(&an_edge_shape);

                let an_ori = an_edge_shape.orientation();
                let (an_ef, an_er) = if an_ori == Orientation::Forward {
                    let mut rev = an_edge_shape.clone();
                    rev.set_orientation(Orientation::Reversed);
                    (an_edge_shape.clone(), rev)
                } else {
                    let mut rev = an_edge_shape.clone();
                    rev.set_orientation(Orientation::Forward);
                    (rev, an_edge_shape.clone())
                };

                let Some(a_tr_cf) = self.translate_curve2d(&an_ef, a_copy_face, &a_tr_v) else {
                    continue;
                };
                let Some(a_tr_cr) = self.translate_curve2d(&an_er, a_copy_face, &a_tr_v) else {
                    continue;
                };

                if a_try_flag {
                    // Use the builder in a trying case.
                    if let Some(e) = Edge::wrap(an_edge_shape.clone()) {
                        write_edge_pcurves(&e, a_tr_cf, a_tr_cr, a_copy_face, tol);
                    }
                } else if let Some(e) = Edge::wrap(an_edge_shape.clone()) {
                    // Use the "false builder" otherwise.
                    self.update_edge_seam(&e, a_tr_cf, a_tr_cr, a_copy_face, tol);
                }
            } else {
                // b. Usual edge (`cxx:547-557`).
                let Some(a_tr_c) = self.translate_curve2d(&an_edge, a_copy_face, &a_tr_v) else {
                    continue;
                };
                if a_try_flag {
                    if let Some(e) = Edge::wrap(an_edge_shape.clone()) {
                        write_edge_pcurve(&e, a_tr_c, a_copy_face, tol);
                    }
                } else if let Some(e) = Edge::wrap(an_edge_shape.clone()) {
                    self.update_edge(&e, a_tr_c, a_copy_face, tol);
                }
            }

            let Some((_pf, pl)) = Self::get_p2d_fl(a_copy_face, &an_edge) else {
                continue;
            };
            pb = pl;

            // In a trying case, check the first and last 2d point of the
            // contour (`cxx:562-578`).
            if a_try_flag && i == a_nb_edges {
                let a_d = pa1.distance(&pb);
                if a_d > a_d_tolerance {
                    return 1;
                }
            }
        }

        0
    }

    /// `GetP2dFL(aFace, anEdge, P2dF, P2dL)` (`cxx:687-702`).
    ///
    /// `BRep_Tool::CurveOnSurface(E, F, First, Last)` and its two endpoints,
    /// swapped when the edge is `TopAbs_REVERSED`.
    pub fn get_p2d_fl(a_face: &Face, an_edge: &Edge) -> Option<(GpPnt2d, GpPnt2d)> {
        let (c2d, first, last) = crate::boptools_2d::curve_on_surface_range(an_edge, a_face)?;
        let mut p2d_f = c2d.d0(first);
        let mut p2d_l = c2d.d0(last);
        if an_edge.0.orientation() == Orientation::Reversed {
            std::mem::swap(&mut p2d_f, &mut p2d_l);
        }
        Some((p2d_f, p2d_l))
    }

    /// `MakeHeadList(aFirstEdge, HeadList)` (`cxx:706-744`).
    ///
    /// Rotates `myOrderedWireList` so that it starts at `aFirstEdge`: the
    /// partial list from `aFirstEdge` on, then the edges before it.
    pub(super) fn make_head_list(
        &self,
        a_first_edge: Option<&TopoShape>,
        head_list: &mut Vec<TopoShape>,
    ) {
        let mut a_fe: Option<TopoShape> = a_first_edge.cloned();
        let mut a_tail_list: Vec<TopoShape> = Vec::new();
        let mut a_flag = false;

        for an_edge in &self.my_ordered_wire_list {
            // A null `aFirstEdge` falls back to the first non-degenerated edge
            // (NIZHNY-MZV, `cxx:718-721`).
            if a_fe.is_none() && !edge_is_degenerated(an_edge) {
                a_fe = Some(an_edge.clone());
            }
            if let Some(fe) = &a_fe {
                if is_equal(an_edge, fe) {
                    // Turn the switch ON.
                    a_flag = true;
                }
            }
            if a_flag {
                head_list.push(an_edge.clone());
            }
        }

        for an_edge in &self.my_ordered_wire_list {
            if let Some(fe) = &a_fe {
                if is_equal(an_edge, fe) {
                    break;
                }
            }
            a_tail_list.push(an_edge.clone());
        }
        head_list.extend(a_tail_list);
    }

    /// `CheckList(aFace, HeadList)` (`cxx:748-847`).
    ///
    /// Reverses the head list when the first edge of the contour runs away
    /// from the second one. No caller in OCCT 8.0.0: both call sites in
    /// `ConnectWire` are commented out (`cxx:427-428`).
    pub fn check_list(a_face: &Face, head_list: &mut Vec<TopoShape>) {
        let a_d_tolerance = occt_core::precision::CONFUSION;

        let a_seq: Vec<TopoShape> = head_list.clone();
        let mut a_nb = a_seq.len() as i32;
        let Some(a_seq_first) = a_seq.first().cloned() else {
            return;
        };

        let mut r1 = a_seq_first.orientation();
        let mut a_first_check = false;

        if r1 == Orientation::Reversed {
            // When the first edge is reversed and everything in 2d is fine, a
            // single FORWARD edge is enough to keep the list as it is
            // (`cxx:763-793`).
            let Some(a_fst) = Edge::wrap(a_seq_first.clone()) else {
                return;
            };
            let Some((_pa, mut pb)) = Self::get_p2d_fl(a_face, &a_fst) else {
                return;
            };
            for i in 2..=a_nb {
                let Some(an_edge_shape) = a_seq.get((i - 1) as usize).cloned() else {
                    continue;
                };
                let Some(an_edge) = Edge::wrap(an_edge_shape.clone()) else {
                    continue;
                };
                let Some((pf, pl)) = Self::get_p2d_fl(a_face, &an_edge) else {
                    continue;
                };
                if pf.distance(&pb) < a_d_tolerance {
                    pb = pl;
                    if an_edge_shape.orientation() == Orientation::Forward {
                        a_first_check = true;
                        break;
                    }
                }
            }
        }

        r1 = a_seq_first.orientation();

        // `aV1R` / `aV2F` (`cxx:795-827`).
        let mut a_v1r: Option<TopoShape> = None;
        for v in explorer_shapes(&a_seq_first, ShapeType::Vertex) {
            if v.orientation() == Orientation::Reversed {
                a_v1r = Some(v);
                break;
            }
        }

        let a_second_edge = a_seq.get(1).cloned();
        let mut a_v2f: Option<TopoShape> = None;
        if let Some(second) = &a_second_edge {
            for v in explorer_shapes(second, ShapeType::Vertex) {
                if v.orientation() == Orientation::Forward {
                    a_v2f = Some(v);
                    break;
                }
            }
        }

        let (Some(a_v1r), Some(a_v2f)) = (
            a_v1r.and_then(crate::shape::Vertex::wrap),
            a_v2f.and_then(crate::shape::Vertex::wrap),
        ) else {
            // `BRep_Tool::Pnt` on a null vertex throws in OCCT.
            return;
        };
        let p1 = vertex_position(&a_v1r);
        let p2 = vertex_position(&a_v2f);

        // NIZNHY-PKV, `cxx:829-846`.
        let dist = p1.distance(&p2);
        if (!(dist < a_d_tolerance) && r1 == Orientation::Forward)
            || (!a_first_check && r1 == Orientation::Reversed)
        {
            // The order of the list has to be reversed.
            let mut extended = a_seq.clone();
            extended.push(a_seq_first);
            head_list.clear();
            a_nb = extended.len() as i32;
            for i in (2..=a_nb).rev() {
                let Some(s) = extended.get((i - 1) as usize).cloned() else {
                    continue;
                };
                head_list.push(s);
            }
        }
    }

}

/// `BRep_Tool::Degenerated(TopoDS::Edge(anEdge))` (`cxx:719`).
fn edge_is_degenerated(an_edge: &TopoShape) -> bool {
    Edge::wrap(an_edge.clone())
        .map(|e| BRepTool::is_degenerated(&e))
        .unwrap_or(false)
}
