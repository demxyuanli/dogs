//! ShapeFix_ComposeShell::DispatchWires (ShapeFix_ComposeShell.cxx:3275-3584).

use std::sync::Arc;

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_core::precision::PCONFUSION;
use occt_geom::Surface;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape};
use crate::shhealing::fix_shifted_wire;
use crate::tgeometry::{EdgeGeom, GeometryRegistry};
use crate::shhealing::{ShapeBuildEdge, ShapeFixEdge};

use super::helpers::{compose_adjust_to_period, middle_point, TOLINT};
use super::reshape::{MapReShape, ReShape};
use super::shell::ComposeShell;


use super::wire_segment::WireSegment;

fn oriented_forward(edge: &Edge) -> Edge {
    let mut e = edge.clone();
    e.0.set_orientation(Orientation::Forward);
    e
}

/// `ShapeFix_ComposeShell.cxx:3287-3327` ("pdn: shift pcurves in the seam to
/// make OK shape w/o fixshifted"): for every REVERSED closed-on-face edge whose
/// two representations coincide, translate one of them by a period so the two
/// occurrences of the seam land on the two ends of the periodic range.
///
/// `c21` is `BRep_Tool::CurveOnSurface(E, myFace)` and `c22` is
/// `CurveOnSurface(E.Reversed(), myFace)` (`cxx:3303-3307`), i.e. the *face*
/// overload: `myFace`'s own orientation decides which representation each call
/// returns (`BRep_Tool.cxx:308-314` composes `F.Orientation()` into the edge
/// view before the closed-surface arm at `cxx:354-357` picks `PCurve1` /
/// `PCurve2`). `myFace` is the FORWARD view always (`cxx:163-167` in `Init`), so
/// for a REVERSED `E` `c21` reads `PCurve2` and `c22` reads `PCurve1`; `c22` is
/// the translated one (`cxx:3322`), i.e. the period goes on `PCurve1` of the
/// EDGE's own orientation. Picking the slot by the stored edge orientation
/// alone swaps the pair and flips the UV winding of the wire.
fn shift_closed_seam_pcurves(
    face: &Face,
    wires: &[WireSegment],
    u_closed: bool,
    v_closed: bool,
    u_period: f64,
    v_period: f64,
) {
    let reg = GeometryRegistry::global();
    let face_key = GeometryRegistry::shape_key(&face.0);
    for w in wires.iter() {
        if w.is_vertex() {
            continue; // cxx:3291-3294
        }
        for jl in 1..=w.nb_edges() {
            let Some(e) = w.edge(jl) else { continue };
            // cxx:3300: `E.Orientation() == REVERSED && IsClosed(E, myFace)`.
            if e.0.orientation() != Orientation::Reversed || !BRepTool::is_closed_edge_face(e, face) {
                continue;
            }
            let mut erev = e.clone();
            erev.0.reverse();
            let (Some((c21, f1, l1)), Some((c22, f2, l2))) = (
                crate::boptools_2d::curve_on_surface_range(e, face),
                crate::boptools_2d::curve_on_surface_range(&erev, face),
            ) else {
                continue;
            };
            let pf1 = c21.d0(f1);
            let pl1 = c21.d0(l1);
            let pf2 = c22.d0(f2);
            let pl2 = c22.d0(l2);
            let d_preci = PCONFUSION * PCONFUSION;
            if !(Arc::ptr_eq(&c21, &c22) // cxx:3312
                || pf1.square_distance(&pf2) < d_preci
                || pl1.square_distance(&pl2) < d_preci)
            {
                continue;
            }
            let mut shift = GpVec2d::new(0.0, 0.0);
            if u_closed && (pf2.x() - pl2.x()).abs() < PCONFUSION {
                shift = GpVec2d::new(u_period, shift.y()); // cxx:3316-3319
            }
            if v_closed && (pf2.y() - pl2.y()).abs() < PCONFUSION {
                shift = GpVec2d::new(shift.x(), v_period);
            }
            // cxx:3322: `c22->Translate(shift)` on the shared handle, so the
            // same-object case moves both slots.
            let mut nc = c22.clone_dyn();
            let mut tr = occt_core::gp::GpTrsf2d::identity();
            tr.set_translation_vec(&shift);
            nc.transform(&tr);
            let mut pcs = reg.edge_pcurves(&e.0, face_key);
            if Arc::ptr_eq(&c21, &c22) {
                // A single shared representation: OCCT translates the object
                // once, so every slot that held it reads the shifted curve.
                for p in pcs.iter_mut() {
                    let mut c = p.clone_dyn();
                    c.transform(&tr);
                    *p = Arc::from(c);
                }
            } else {
                match pcs.iter().position(|p| Arc::ptr_eq(p, &c22)) {
                    Some(k) => pcs[k] = Arc::from(nc),
                    None => continue,
                }
            }
            reg.set_edge_pcurves(&e.0, face_key, pcs);
        }
    }
}

impl ComposeShell {
    /// `DispatchWires(faces, wires)` (`cxx:3275-3584`).
    pub fn dispatch_wires(&mut self, faces: &mut Vec<TopoShape>, wires: &mut Vec<WireSegment>) {
        let builder = TopoBuilder::new();
        let face_val = match self.face() {
            Some(f) => f.clone(),
            None => return,
        };

        // cxx:3280-3360: closed mode fixes.
        if self.closed_mode {
            // cxx:3287-3327: shift the seam pcurves.
            shift_closed_seam_pcurves(
                &face_val,
                wires,
                self.u_closed,
                self.v_closed,
                self.u_period,
                self.v_period,
            );

            let mut i = 0usize;
            while i < wires.len() {
                if wires[i].is_vertex() {
                    i += 1;
                    continue;
                }
                // cxx:3337-3343: skip a wire with a single degenerated edge.
                if wires[i].nb_edges() == 0
                    || (wires[i].nb_edges() == 1 && BRepTool::is_degenerated(&wires[i].edges()[0]))
                {
                    wires.remove(i);
                    continue;
                }
                // cxx:3345-3346: sfw.Load(sbwd); sfw.FixShifted().
                let wire = builder.make_wire(wires[i].edges());
                let _ = fix_shifted_wire(&wire, &face_val);
                // cxx:3348-3357: ShapeBuild_Edge::RemovePCurve on degenerated edges.
                let sbe = ShapeBuildEdge;
                for j in 0..wires[i].nb_edges() {
                    if BRepTool::is_degenerated(&wires[i].edges()[j]) {
                        sbe.remove_pcurve(&wires[i].edges()[j], &face_val);
                    }
                }
                // cxx:3358: sfw.FixDegenerated(). The port's `fix_degenerated_all`
                // rewrites a `Wire`'s edge list, so rebuild the segment and
                // restore the patch indices of the surviving edges.
                let old_edges: Vec<Edge> = wires[i].edges().to_vec();
                let old_patches: Vec<Option<(i32, i32, i32, i32)>> = (1..=old_edges.len())
                    .map(|j| wires[i].patch_index(j))
                    .collect();
                let mut dw = builder.make_wire(&old_edges);
                if crate::shhealing::fix_degenerated_all(&mut dw, &face_val, self.precision) {
                    let new_edges = crate::topo_tools::edges_of_wire(&dw);
                    let manifold = wires[i].manifold_mode();
                    wires[i].clear();
                    wires[i].load_edges(&new_edges, manifold);
                    for (k, e) in new_edges.iter().enumerate() {
                        if let Some(pos) = old_edges
                            .iter()
                            .position(|o| crate::topo_tools::is_same(&o.0, &e.0))
                        {
                            if let Some((a, b, c, d)) = old_patches[pos] {
                                wires[i].set_patch_index(k + 1, a, b, c, d);
                            }
                        }
                    }
                }
                i += 1;
            }
        }

        // cxx:3362-3376: center points.
        let nb = wires.len();
        if nb == 0 {
            return;
        }
        let mut m_pnts: Vec<GpPnt2d> = Vec::with_capacity(nb);
        for i in 0..nb {
            m_pnts.push(middle_point(&wires[i], &face_val));
        }

        // cxx:3378-3382.
        let mut rs = MapReShape::new();
        let sbe = ShapeBuildEdge;
        let mut sfe = ShapeFixEdge::new();
        let reg = GeometryRegistry::global();
        let (u1, u2, v1, v2) = self.grid.bounds();

        // cxx:3362: `DispatchWires` has a single closed-mode pass - the seam
        // period shift (cxx:3287-3327) followed by one `FixShifted` /
        // `FixDegenerated` per wire (cxx:3330-3360), both run above. Nothing
        // re-shifts or re-fixes here; the wire data is already final.

        // cxx:3387-3533.
        for i in 0..nb {
            let mut pnt = m_pnts[i];
            let mut ush = 0.0f64;
            let mut vsh = 0.0f64;
            if self.u_closed {
                ush = compose_adjust_to_period(pnt.x(), u1, u2);
                pnt = GpPnt2d::new(pnt.x() + ush, pnt.y());
            }
            if self.v_closed {
                vsh = compose_adjust_to_period(pnt.y(), v1, v2);
                pnt = GpPnt2d::new(pnt.x(), pnt.y() + vsh);
            }
            m_pnts[i] = pnt;
            let ind_u = self.grid.locate_u_parameter(pnt.x());
            let ind_v = self.grid.locate_v_parameter(pnt.y());
            let (u_fact, mut trsf, mut need_t) = self.grid.global_to_local_transformation(ind_u, ind_v);
            if ush != 0.0 || vsh != 0.0 {
                // cxx:3411-3414: T.Multiply(Sh).
                let mut sh = occt_core::gp::GpTrsf2d::identity();
                sh.set_translation_vec(&occt_core::gp::GpVec2d::new(ush, vsh));
                trsf = trsf.multiplied(&sh);
                need_t = true;
            }
            if wires[i].is_vertex() {
                continue; // cxx:3416-3419
            }
            let surf: Arc<dyn Surface> = match self.grid.patch(ind_u, ind_v) {
                Some(s) => s.clone(),
                None => continue,
            };
            let face = builder.make_face(surf, &[]);
            let face_key = GeometryRegistry::shape_key(&face.0);
            for j in 1..=wires[i].nb_edges() {
                let edge = wires[i].edge(j).expect("edge").clone();
                let is_manifold = matches!(edge.0.orientation(), Orientation::Forward | Orientation::Reversed);
                let mut an_init_edge = edge.clone();
                let mut new_edge;
                if rs.is_recorded(&edge.0) {
                    // cxx:3438-3443.
                    new_edge = Edge(rs.value(&edge.0).expect("recorded"));
                } else {
                    if !is_manifold {
                        an_init_edge.0.set_orientation(Orientation::Forward);
                    }
                    new_edge = sbe.copy(&an_init_edge, false); // cxx:3451
                    if !is_manifold {
                        new_edge.0.set_orientation(edge.0.orientation()); // cxx:3453-3455
                    }
                    rs.replace(&edge.0, &new_edge.0); // cxx:3456
                    self.context_mut().replace(&edge.0, &new_edge.0); // cxx:3457
                }
                sbe.reassign_pcurve(&new_edge, &face_val, &face); // cxx:3460

                // cxx:3462-3504: transform the pcurve into the patch's
                // parametric space.
                if need_t {
                    if let Some((c2d, f, l)) = crate::boptools_2d::curve_on_surface_oriented(&new_edge, &face, false) {
                        let mut newf = f;
                        let mut newl = l;
                        let c2dnew = sbe.transform_pcurve(&c2d, &trsf, u_fact, &mut newf, &mut newl);
                        // cxx:3471-3493: a closed edge keeps both pcurves.
                        // `PCurve(newEdge)` is transformed; `PCurve(newEdge.Reversed())`
                        // stays. FORWARD stores (transformed, other); REVERSED
                        // stores (other, transformed).
                        if BRepTool::is_closed_edge_face(&new_edge, &face) {
                            let mut rev = new_edge.clone();
                            rev.0.reverse();
                            if let Some((c2d2, _, _)) =
                                crate::boptools_2d::curve_on_surface_oriented(&rev, &face, false)
                            {
                                let pcs = if new_edge.0.orientation() == Orientation::Forward {
                                    vec![c2dnew, c2d2]
                                } else {
                                    vec![c2d2, c2dnew]
                                };
                                reg.set_edge_pcurves(&new_edge.0, face_key, pcs);
                            } else {
                                reg.set_edge_pcurve(&new_edge.0, face_key, c2dnew);
                            }
                        } else {
                            reg.set_edge_pcurve(&new_edge.0, face_key, c2dnew);
                        }
                        reg.set_pcurve_range(&new_edge.0, face_key, newf, newl);
                        if (newf != f || newl != l) && !BRepTool::is_degenerated(&new_edge) {
                            reg.set_same_range(&new_edge.0, false); // cxx:3499-3502
                        }
                    }
                }

                let same_range = reg.edge_geom(&new_edge.0).map(|g| g.same_range).unwrap_or(true);
                if !same_range {
                    // cxx:3506-3526.
                    let etmp = if !is_manifold {
                        let afe = oriented_forward(&new_edge);
                        sbe.copy(&afe, false)
                    } else {
                        sbe.copy(&new_edge, false)
                    };
                    sfe.fix_add_curve3d(&etmp);
                    if let Some(c3d) = reg.edge_curve(&etmp.0) {
                        let (cf, cl) = reg.edge_parameters(&etmp.0);
                        let mut g = EdgeGeom::new(c3d, cf, cl);
                        if let Some(o) = reg.edge_geom(&new_edge.0) {
                            g.tolerance = o.tolerance;
                            g.same_parameter = o.same_parameter;
                            g.same_range = o.same_range;
                            g.degenerated = o.degenerated;
                        }
                        reg.set_edge(&new_edge.0, g);
                        sbe.set_range3d(&new_edge, cf, cl);
                    }
                } else {
                    sfe.fix_add_curve3d(&new_edge); // cxx:3529
                }
                wires[i].set_edge(j, new_edge); // cxx:3531
                let _ = face_key;
            }
        }

        // cxx:3535-3583: collect packets by surface and dispatch.
        let mut used = vec![false; nb];
        loop {
            let mut loops: Vec<TopoShape> = Vec::new();
            let mut surf: Option<Arc<dyn Surface>> = None;
            for i in 0..nb {
                if used[i] {
                    continue;
                }
                let s = self.grid.patch_pnt(&m_pnts[i]).cloned(); // cxx:3549
                // cxx:3550-3556. OCCT compares `Geom_Surface` handles, i.e.
                // pointers, and the *first* unused wire of the packet sets
                // `Surf` even when its `Patch` is null: once `Surf` is null
                // every following unused wire is taken into the packet on its
                // own surface. The packet is then dispatched on whatever `Surf`
                // ended up being (cxx:3577-3580), never on a per-wire surface.
                if surf.is_none() {
                    surf = s; // cxx:3550-3552
                } else if !same_surface(&surf, &s) {
                    continue; // cxx:3553-3556
                }
                used[i] = true;
                if wires[i].is_vertex() {
                    if let Some(v) = wires[i].vertex() {
                        if v.0.orientation() == Orientation::Internal {
                            loops.push(v.0.clone()); // cxx:3561-3567
                        }
                    }
                } else {
                    let w = builder.make_wire(wires[i].edges()); // cxx:3571-3574
                    loops.push(w.0);
                }
            }
            let Some(surf) = surf else {
                break; // cxx:3577-3580
            };
            self.make_faces_on_patch(faces, &surf, &mut loops); // cxx:3582
        }
        let _ = (ShapeType::Face, TOLINT);
    }
}

/// `S != Surf` (`ShapeFix_ComposeShell.cxx:3553`) on two `Geom_Surface`
/// handles: OCCT's `Handle` comparison is a pointer comparison, so two
/// distinct `Geom_Surface` instances never match however equal they are; a
/// null handle equals only another null handle.
fn same_surface(a: &Option<Arc<dyn Surface>>, b: &Option<Arc<dyn Surface>>) -> bool {
    match (a, b) {
        (Some(x), Some(y)) => Arc::ptr_eq(x, y),
        (None, None) => true,
        _ => false,
    }
}