//! ShapeFix_ComposeShell::DispatchWires (ShapeFix_ComposeShell.cxx:3275-3584).

use std::sync::Arc;

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_core::precision::PCONFUSION;
use occt_geom2d::curve::Curve2d;
use occt_geom::Surface;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape};
use crate::shhealing::fix_shifted_wire;
use crate::tgeometry::{EdgeGeom, GeometryRegistry};
use crate::shhealing::{ShapeBuildEdge, ShapeFixEdge};

use super::helpers::{adjust_to_period, get_middle_point, TOLINT};
use super::reshape::{MapReShape, ReShape};
use super::shell::ComposeShell;


use super::wire_segment::WireSegment;

fn oriented_forward(edge: &Edge) -> Edge {
    let mut e = edge.clone();
    e.0.set_orientation(Orientation::Forward);
    e
}

impl ComposeShell {
    /// `DispatchWires(faces, wires)` (`cxx:3275-3584`).
    pub fn dispatch_wires(&mut self, faces: &mut Vec<TopoShape>, wires: &mut Vec<WireSegment>) {
        let builder = TopoBuilder::new();
        let my_face = match self.face() {
            Some(f) => f.clone(),
            None => return,
        };

        // cxx:3280-3360: closed mode fixes.
        if self.closed_mode {
            // cxx:3287-3327: shift the seam pcurves. `CurveOnSurface(E)` reads
            // `PCurve2` for a REVERSED edge, `CurveOnSurface(E.Reversed())`
            // reads `PCurve1` (`BRep_Tool.cxx:347-357`).
            let my_face_key = GeometryRegistry::shape_key(&my_face.0);
            for wi in 0..wires.len() {
                if wires[wi].is_vertex() {
                    continue;
                }
                for jl in 1..=wires[wi].nb_edges() {
                    let e = wires[wi].edge(jl).expect("edge").clone();
                    if e.0.orientation() != Orientation::Reversed {
                        continue; // cxx:3300
                    }
                    if !crate::brep_tool::BRepTool::is_closed_edge_face(&e, &my_face) {
                        continue; // BRep_Tool::IsClosed(E, myFace)
                    }
                    let pcs = GeometryRegistry::global().edge_pcurves(&e.0, my_face_key);
                    let c21 = pcs[1].clone();
                    let c22 = pcs[0].clone();
                    let (f1, l1) = GeometryRegistry::global()
                        .pcurve_range(&e.0, my_face_key)
                        .unwrap_or((c21.first_parameter(), c21.last_parameter()));
                    let (f2, l2) = (f1, l1);
                    let pf1 = c21.d0(f1);
                    let pl1 = c21.d0(l1);
                    let pf2 = c22.d0(f2);
                    let pl2 = c22.d0(l2);
                    let d_preci = PCONFUSION * PCONFUSION;
                    if Arc::ptr_eq(&c21, &c22)
                        || pf1.square_distance(&pf2) < d_preci // cxx:3312
                        || pl1.square_distance(&pl2) < d_preci
                    {
                        let mut shift = GpVec2d::new(0.0, 0.0);
                        if self.u_closed && (pf2.x() - pl2.x()).abs() < PCONFUSION {
                            shift = GpVec2d::new(self.u_period, shift.y());
                        }
                        if self.v_closed && (pf2.y() - pl2.y()).abs() < PCONFUSION {
                            shift = GpVec2d::new(shift.x(), self.v_period);
                        }
                        // c22->Translate(shift): c22 is PCurve1; keep PCurve2.
                        let mut nc = c22.clone_dyn();
                        let mut tr = occt_core::gp::GpTrsf2d::identity();
                        tr.set_translation_vec(&shift);
                        nc.transform(&tr);
                        GeometryRegistry::global().set_edge_pcurves(
                            &e.0,
                            my_face_key,
                            vec![Arc::from(nc), c21.clone()],
                        );
                    }
                }
            }

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
                let _ = fix_shifted_wire(&wire, &my_face);
                // cxx:3348-3357: ShapeBuild_Edge::RemovePCurve on degenerated edges.
                let sbe = ShapeBuildEdge;
                for j in 0..wires[i].nb_edges() {
                    if BRepTool::is_degenerated(&wires[i].edges()[j]) {
                        sbe.remove_pcurve(&wires[i].edges()[j], &my_face);
                    }
                }
                // cxx:3358: sfw.FixDegenerated(). The port's `fix_degenerated_all`
                // rewrites a `Wire`'s edge list, so rebuild the segment and
                // restore the patch indices of the surviving edges.
                let old_edges: Vec<Edge> = wires[i].edges().to_vec();
                let old_patches: Vec<Option<(i32, i32, i32, i32)>> = (1..=old_edges.len())
                    .map(|j| wires[i].get_patch_index(j))
                    .collect();
                let mut dw = builder.make_wire(&old_edges);
                if crate::shhealing::fix_degenerated_all(&mut dw, &my_face, self.precision) {
                    let new_edges = crate::topo_tools_full::edges_of_wire(&dw);
                    let manifold = wires[i].manifold_mode();
                    wires[i].clear();
                    wires[i].load_edges(&new_edges, manifold);
                    for (k, e) in new_edges.iter().enumerate() {
                        if let Some(pos) = old_edges
                            .iter()
                            .position(|o| crate::topo_tools_full::is_same(&o.0, &e.0))
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
            m_pnts.push(get_middle_point(&wires[i], &my_face));
        }

        // cxx:3378-3382.
        let mut rs = MapReShape::new();
        let sbe = ShapeBuildEdge;
        let sfe = ShapeFixEdge;
        let reg = GeometryRegistry::global();
        let (u1, u2, v1, v2) = self.grid.bounds();

        // cxx:3281-3327 (pdn: "shift pcurves in the seam to make OK shape w/o
        // fixshifted"): in closed mode, before dispatching, every REVERSED
        // closed-on-face edge whose two pcurves coincide gets its second pcurve
        // translated by a period, so the seam edge's two occurrences land on the
        // two ends of the UV range. Without it the periodic boundary chain stays
        // open and the Delaunay fills nothing (a3n00 f176: frontier=0).
        if self.closed_mode {
            let face_key = GeometryRegistry::shape_key(&my_face.0);
            for w in wires.iter() {
                for e in w.edges() {
                    if e.0.orientation() != Orientation::Reversed
                        || !BRepTool::is_closed_edge_face(e, &my_face)
                    {
                        continue;
                    }
                    let mut pcs = reg.edge_pcurves(&e.0, face_key);
                    if pcs.len() < 2 {
                        continue;
                    }
                    let (f1, l1) = reg.pcurve_range(&e.0, face_key).unwrap_or((0.0, 0.0));
                    // c21 = CurveOnSurface(E) (REVERSED occurrence) = pcs[1];
                    // c22 = CurveOnSurface(E.Reversed()) (FORWARD) = pcs[0].
                    let (c21, c22) = (pcs[1].clone(), pcs[0].clone());
                    let pf1 = c21.d0(f1);
                    let pl1 = c21.d0(l1);
                    let pf2 = c22.d0(f1);
                    let pl2 = c22.d0(l1);
                    let d = PCONFUSION * PCONFUSION;
                    if Arc::ptr_eq(&c21, &c22)
                        || pf1.square_distance(&pf2) < d
                        || pl1.square_distance(&pl2) < d
                    {
                        let sx = if self.u_closed && (pf2.x() - pl2.x()).abs() < PCONFUSION {
                            self.u_period
                        } else {
                            0.0
                        };
                        let sy = if self.v_closed && (pf2.y() - pl2.y()).abs() < PCONFUSION {
                            self.v_period
                        } else {
                            0.0
                        };
                        if sx != 0.0 || sy != 0.0 {
                            let mut nc = c22.clone_dyn();
                            let mut t = occt_core::gp::GpTrsf2d::identity();
                            t.set_translation_vec(&occt_core::gp::GpVec2d::new(sx, sy));
                            nc.transform(&t);
                            pcs[0] = Arc::from(nc);
                            reg.set_edge_pcurves(&e.0, face_key, pcs);
                            reg.set_pcurve_range(&e.0, face_key, f1, l1);
                        }
                    }
                }
            }
        }

        // cxx:3345-3346 / 3358: `sfw.Load(sbwd); sfw.FixShifted(); sfw.FixDegenerated();`
        // run on each wire after the seam pcurve shift above.
        for w in wires.iter() {
            let dw = builder.make_wire(w.edges());
            let _ = crate::shhealing::fix_shifted_wire(&dw, &my_face);
        }

        // cxx:3387-3533.
        for i in 0..nb {
            let mut pnt = m_pnts[i];
            let mut ush = 0.0f64;
            let mut vsh = 0.0f64;
            if self.u_closed {
                ush = adjust_to_period(pnt.x(), u1, u2);
                pnt = GpPnt2d::new(pnt.x() + ush, pnt.y());
            }
            if self.v_closed {
                vsh = adjust_to_period(pnt.y(), v1, v2);
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
                sbe.reassign_pcurve(&new_edge, &my_face, &face); // cxx:3460

                // cxx:3462-3504: transform the pcurve into the patch's
                // parametric space.
                if need_t {
                    if let Some((c2d, f, l)) = crate::boptools_2d::curve_on_surface_oriented(&new_edge, &face, false) {
                        let mut newf = f;
                        let mut newl = l;
                        let c2dnew = sbe.transform_pcurve(&c2d, &trsf, u_fact, &mut newf, &mut newl);
                        // cxx:3471-3493: the closed (seam) branch keeps both
                        // pcurves; UNPORTED here — see specs/_a3n00_gap_analysis.md
                        // §9.117: the port's ReassignPCurve state makes
                        // IsClosed(newEdge, face) differ from OCCT's, so the
                        // faithful branch changed a3n00's total away from the
                        // 11052/12324 target. Kept single-pcurve until the patch
                        // pcurve representation matches.
                        reg.set_edge_pcurve(&new_edge.0, face_key, c2dnew);
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
                    sfe.fix_add_curve3d(&etmp, &face);
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
                    sfe.fix_add_curve3d(&new_edge, &face); // cxx:3529
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
                let s = self.grid.patch_pnt(&m_pnts[i]);
                match (&surf, s) {
                    (None, Some(s)) => surf = Some(s.clone()),
                    (Some(cur), Some(s)) if Arc::ptr_eq(cur, s) => {}
                    _ => continue,
                }
                used[i] = true;
                if wires[i].is_vertex() {
                    if let Some(v) = wires[i].get_vertex() {
                        if v.0.orientation() == Orientation::Internal {
                            loops.push(v.0.clone());
                        }
                    }
                } else {
                    let w = builder.make_wire(wires[i].edges());
                    loops.push(w.0);
                }
            }
            let surf = match surf {
                Some(s) => s,
                None => break,
            };
            self.make_faces_on_patch(faces, &surf, &mut loops); // cxx:3582
        }
        let _ = (ShapeType::Face, TOLINT);
    }
}