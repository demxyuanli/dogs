//! `BRepOffset_Tool::Inter2d` (`BRepOffset_Tool.cxx:2253-2532`): the 2D
//! intersection of two edges' pcurves on one face, answering the intersection
//! vertices in `LV`.
//!
//! Control flow, in the reference's order:
//! 1. `BuildPCurves` on both edges (`cxx:2266-2267`) - the pcurves of an edge
//!    extended by `ExtentEdge` are not necessarily up to date;
//! 2. `NbPC1` / `NbPC2` = 1 + `BRep_Tool::IsClosed(E, F)` (`cxx:2286-2292`):
//!    a closed edge contributes its reversed pcurve as a second candidate;
//! 3. two passes (`itry` 0 then 1, `cxx:2296`): pass 0 tests the four endpoint
//!    pairs for coincidence and then projects each endpoint of one pcurve on the
//!    other; pass 1 re-reads the ranges from the pcurves themselves before the
//!    `Geom2dInt_GInter` intersection. A solution in pass 0 still lets the
//!    remaining (i, j) pairs run, because only the *tests* are guarded by
//!    `!YaSol` (`cxx:2372`, `:2404`, `:2436`) - the loops are not.
//! 4. each hit becomes a vertex at `S->Value(P2d)` with `INTERNAL`
//!    orientation, its parameter written on both edges (`cxx:2478-2486`);
//! 5. more than one vertex is compacted to the pair with the smallest and
//!    largest parameter on `E1` (`cxx:2492-2523`).
//!
//! Port guards (no OCCT counterpart): a pcurve that neither exists nor can be
//! built skips the (i, j) pair instead of dereferencing a null handle, and a
//! face without a surface leaves `LV` empty.
//!
//! `Extrema_ExtPC2d(P, AC)` is `ExtremaExtPC2d::new_range` over the adaptor's
//! `[fl1[0], fl1[1]]`; `Geom2dInt_GInter(AC1, AC2, TolConf, TolConf)` is
//! `Geom2dIntGInter` with the two `domain_over_range` domains.

use occt_core::gp::GpPnt2d;
use occt_core::precision::{CONFUSION, INFINITE, SQUARE_CONFUSION};
use occt_geom2d::extrema2d::ExtremaExtPC2d;
use occt_geom2d::geom2d_int::Geom2dIntGInter;
use occt_geom2d::Curve2d;
use std::sync::Arc;

use crate::abs::Orientation;
use crate::brep_builder_update_vertex::update_vertex_on_edge;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, Vertex};
use crate::tgeometry::GeometryRegistry;

use super::pcurve_builder::build_pcurves;
use super::BRepOffsetTool;

impl BRepOffsetTool {
    /// `Inter2d(F, E1, E2, LV, TolConf)` (`cxx:2253-2532`).
    ///
    /// `lv` is cleared first and receives the intersection vertices
    /// (`INTERNAL` orientation).
    pub fn inter2d(f: &Face, e1: &Edge, e2: &Edge, lv: &mut Vec<Vertex>, tol_conf: f64) {
        lv.clear();

        // `cxx:2266-2267`.
        let _ = build_pcurves(e1, f);
        let _ = build_pcurves(e2, f);

        // `cxx:2286-2292`: `BRep_Tool::IsClosed(E, F)` counts the pcurves on
        // the face, so a second entry means the edge lies twice on it.
        let face_key = GeometryRegistry::shape_key(&f.0);
        let reg = GeometryRegistry::global();
        let nb_pc1 = if reg.edge_pcurves(&e1.0, face_key).len() > 1 {
            2
        } else {
            1
        };
        let nb_pc2 = if reg.edge_pcurves(&e2.0, face_key).len() > 1 {
            2
        } else {
            1
        };

        // `cxx:2291-2292`.
        let Some(s) = BRepTool::face_surface(f) else {
            return;
        };

        // `cxx:2296`.
        let mut ya_sol = false;
        let mut itry = 0;
        while !ya_sol && itry < 2 {
            for i in 1..=nb_pc1 {
                // `cxx:2301-2308`.
                let Some(c1) = curve_for(e1, f, i) else {
                    continue;
                };
                for j in 1..=nb_pc2 {
                    // `cxx:2314-2322`.
                    let Some(c2) = curve_for(e2, f, j) else {
                        continue;
                    };
                    let (c1pc, c2pc) = (&c1.1, &c2.1);
                    // `cxx:2298-2299` / `:2336-2343`: the ranges are the
                    // out-parameters of `CurveOnSurface`; pass 1 overrides them
                    // with the pcurves' own domain.
                    let (mut fl1, mut fl2) = ((c1.2, c1.3), (c2.2, c2.3));
                    if itry == 1 {
                        fl1 = (c1pc.first_parameter(), c1pc.last_parameter());
                        fl2 = (c2pc.first_parameter(), c2pc.last_parameter());
                    }

                    let mut u1 = 0.0;
                    let mut u2 = 0.0;
                    let mut p2d = GpPnt2d::zero();
                    let mut current_find = false;

                    // `cxx:2347-2432`: only when nothing has been found yet by
                    // an earlier pair.
                    if !ya_sol {
                        let found = Self::inter2d_endpoints(c1pc, fl1, c2pc, fl2, tol_conf);
                        if let Some((a, b, p)) = found {
                            ya_sol = true;
                            current_find = true;
                            u1 = a;
                            u2 = b;
                            p2d = p;
                        }
                    }

                    // `cxx:2436-2474`.
                    if !ya_sol {
                        let mut inter = Geom2dIntGInter::new();
                        let d1 = Geom2dIntGInter::domain_over_range(
                            c1pc.as_ref(),
                            fl1.0,
                            fl1.1,
                            tol_conf,
                        );
                        let d2 = Geom2dIntGInter::domain_over_range(
                            c2pc.as_ref(),
                            fl2.0,
                            fl2.1,
                            tol_conf,
                        );
                        inter.perform(c1pc.as_ref(), &d1, c2pc.as_ref(), &d2, tol_conf, tol_conf);
                        if inter.nb_points() > 0 {
                            ya_sol = true;
                            current_find = true;
                            let p = inter.point(1);
                            u1 = p.param_on_first();
                            u2 = p.param_on_second();
                            p2d = *p.value();
                        } else if inter.nb_segments() > 0 {
                            ya_sol = true;
                            current_find = true;
                            let seg = inter.segment(1);
                            let int_p1 = *seg.first_point();
                            let int_p2 = *seg.last_point();
                            let u1_on_1 = int_p1.param_on_first();
                            let u1_on_2 = int_p2.param_on_first();
                            let u2_on_1 = int_p1.param_on_second();
                            let u2_on_2 = int_p2.param_on_second();
                            u1 = (u1_on_1 + u1_on_2) / 2.0;
                            u2 = (u2_on_1 + u2_on_2) / 2.0;
                            let p2d1 = c1pc.d0(u1);
                            let p2d2 = c2pc.d0(u2);
                            p2d.set_x((p2d1.x() + p2d2.x()) / 2.0);
                            p2d.set_y((p2d1.y() + p2d2.y()) / 2.0);
                        }
                    }

                    // `cxx:2476-2487`.
                    if current_find {
                        let p = s.d0(p2d.x(), p2d.y());
                        // `BRepLib_MakeVertex(P)` (`BRepLib_MakeVertex.cxx:26-31`) uses
                        // `BRepLib::Precision()` = `Precision::Confusion()`.
                        let mut v = TopoBuilder::new().make_vertex(p, CONFUSION);
                        v.0.set_orientation(Orientation::Internal);
                        let e1f = Edge(e1.0.oriented(Orientation::Forward));
                        let e2f = Edge(e2.0.oriented(Orientation::Forward));
                        let _ = update_vertex_on_edge(&v.0, u1, &e1f.0, tol_conf);
                        let _ = update_vertex_on_edge(&v.0, u2, &e2f.0, tol_conf);
                        lv.push(v);
                    }
                }
            }
            itry += 1;
        }

        // `cxx:2492-2523`: keep the vertices closest to the start and the end.
        if lv.len() > 1 {
            let mut vf: Option<Vertex> = None;
            let mut vl: Option<Vertex> = None;
            let mut u_min = INFINITE;
            let mut u_max = -INFINITE;
            let e1f = Edge(e1.0.oriented(Orientation::Forward));
            for cv in lv.iter() {
                let u = BRepTool::parameter_of_vertex(&e1f.0, &cv.0);
                if u < u_min {
                    vf = Some(cv.clone());
                    u_min = u;
                }
                if u > u_max {
                    vl = Some(cv.clone());
                    u_max = u;
                }
            }
            lv.clear();
            if let Some(vf) = vf {
                lv.push(vf);
            }
            if let Some(vl) = vl {
                lv.push(vl);
            }
        }
    }

    /// `cxx:2347-2432`: the two endpoint tests of `Inter2d`'s pass 0.
    ///
    /// Returns `(U1, U2, P2d)` of the first hit, exactly as the reference
    /// assigns it. The endpoint-coincidence loop does not break on a hit: the
    /// last coincident pair wins (`cxx:2356-2375`). The two projection loops do
    /// break (`cxx:2393-2401`, `:2423-2431`).
    fn inter2d_endpoints(
        c1: &Arc<dyn Curve2d>,
        fl1: (f64, f64),
        c2: &Arc<dyn Curve2d>,
        fl2: (f64, f64),
        tol_conf: f64,
    ) -> Option<(f64, f64, GpPnt2d)> {
        let fl1a = [fl1.0, fl1.1];
        let fl2a = [fl2.0, fl2.1];
        let p1 = [c1.d0(fl1.0), c1.d0(fl1.1)];
        let p2 = [c2.d0(fl2.0), c2.d0(fl2.1)];

        // `cxx:2354-2376`.
        let mut result: Option<(f64, f64, GpPnt2d)> = None;
        for i1 in 0..2 {
            for i2 in 0..2 {
                if fl1a[i1].abs() < INFINITE && fl2a[i2].abs() < INFINITE {
                    if p1[i1].is_equal(&p2[i2], tol_conf) {
                        let u1 = fl1a[i1];
                        result = Some((u1, fl2a[i2], c1.d0(u1)));
                    }
                }
            }
        }
        if result.is_some() {
            return result;
        }

        // `cxx:2377-2403`: project each C1 endpoint on C2.
        for i1 in 0..2 {
            let ext = ExtremaExtPC2d::new_range(c2.as_ref(), &p1[i1], fl2a[0], fl2a[1]);
            if !ext.is_done() || ext.nb_ext() == 0 {
                continue;
            }
            let mut dist2_min = ext.square_distance(1);
            let mut index_min = 1;
            for ind in 2..=ext.nb_ext() {
                let dist2 = ext.square_distance(ind);
                if dist2 < dist2_min {
                    dist2_min = dist2;
                    index_min = ind;
                }
            }
            if dist2_min <= SQUARE_CONFUSION {
                let u2 = ext.solutions()[index_min - 1].u;
                return Some((fl1a[i1], u2, p1[i1]));
            }
        }

        // `cxx:2404-2433`: project each C2 endpoint on C1.
        for i2 in 0..2 {
            let ext = ExtremaExtPC2d::new_range(c1.as_ref(), &p2[i2], fl1a[0], fl1a[1]);
            if !ext.is_done() || ext.nb_ext() == 0 {
                continue;
            }
            let mut dist2_min = ext.square_distance(1);
            let mut index_min = 1;
            for ind in 2..=ext.nb_ext() {
                let dist2 = ext.square_distance(ind);
                if dist2 < dist2_min {
                    dist2_min = dist2;
                    index_min = ind;
                }
            }
            if dist2_min <= SQUARE_CONFUSION {
                let u1 = ext.solutions()[index_min - 1].u;
                return Some((u1, fl2a[i2], p2[i2]));
            }
        }
        None
    }
}

/// `cxx:2301-2308` / `:2314-2322`: the pcurve of the `idx`-th candidate -
/// `E` for `idx == 1`, `E.Reversed()` otherwise (`i` / `j` never exceed 2) -
/// together with the range `CurveOnSurface(E, F, f, l)` reports for it.
fn curve_for(e: &Edge, f: &Face, idx: i32) -> Option<(Edge, Arc<dyn Curve2d>, f64, f64)> {
    let cand = if idx == 1 {
        e.clone()
    } else {
        Edge(e.0.oriented(Orientation::Reversed))
    };
    let (pc, first, last) = crate::boptools_2d::curve_on_surface_range(&cand, f)?;
    Some((cand, pc, first, last))
}
