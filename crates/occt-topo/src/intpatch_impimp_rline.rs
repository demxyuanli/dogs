//! `ProcessSegments` / `ProcessRLine` for ImpImp restriction solutions.
//! Source: `IntPatch_ImpImpIntersection.cxx` ProcessSegments, ProcessRLine,
//! IsRLineGood. Transition recomputation (`IntSurf::MakeTransition`) is omitted;
//! domain flags and UV match `PutPointsOnLine`.

use occt_core::gp::GpPnt2d;
use occt_geom::Surface;

use crate::geom_int::{surface_parameters, GLineKind, GeomIntLine, RLine};
use crate::geom_int::RestrictionArc;
use crate::int_tools_wline::PatchPoint;

use super::put::square_distance_to_gline;
use super::sonb::{BoundSegment, PathPoint};

/// `ProcessSegments`: each whole-arc restriction solution becomes an `RLine`.
pub(crate) fn process_segments(
    listedg: &[BoundSegment],
    slin: &mut Vec<GeomIntLine>,
    s1: &dyn Surface,
    s2: &dyn Surface,
    on_first: bool,
    tol_arc: f64,
) {
    for seg in listedg {
        let mut rline = RLine::from_arc(seg.arc, on_first);
        let mut dofirst = seg.first.is_some();
        let mut dolast = seg.last.is_some();
        let mut procf = false;
        let mut procl = false;
        let paramf = seg
            .first
            .map(|p| p.param_on_arc)
            .unwrap_or(seg.arc.first);
        let paraml = seg.last.map(|p| p.param_on_arc).unwrap_or(seg.arc.last);
        let host = if on_first { s1 } else { s2 };
        let degenerate = edge_degenerated(host, &seg.arc, paramf, paraml);

        if (dofirst || dolast) && !degenerate {
            for line in slin.iter() {
                for ptvtx in line_vertices(line) {
                    if dofirst {
                        if let Some(pf) = seg.first {
                            if ptvtx.p.distance(&pf.p) <= tol_arc {
                                let mut newpt = *ptvtx;
                                newpt.param_on_line = paramf;
                                rline.add_vertex(newpt);
                                if !procf {
                                    procf = true;
                                    rline.has_first_point = true;
                                }
                            }
                        }
                    }
                    if dolast {
                        if let Some(pl) = seg.last {
                            if ptvtx.p.distance(&pl.p) <= tol_arc {
                                let mut newpt = *ptvtx;
                                newpt.param_on_line = paraml;
                                rline.add_vertex(newpt);
                                if !procl {
                                    procl = true;
                                    rline.has_last_point = true;
                                }
                            }
                        }
                    }
                }
                if procf {
                    dofirst = false;
                }
                if procl {
                    dolast = false;
                }
            }
        }

        if dofirst {
            if let Some(pf) = seg.first {
                rline.add_vertex(segment_point(pf, s1, s2, on_first, paramf));
                rline.has_first_point = true;
            }
        }
        if dolast {
            if let Some(pl) = seg.last {
                rline.add_vertex(segment_point(pl, s1, s2, on_first, paraml));
                rline.has_last_point = true;
            }
        }
        slin.push(GeomIntLine::Restriction(rline));
    }
}

/// `ProcessRLine`: project GLine vertices onto restriction solutions; drop
/// an RLine that coincides with a GLine unless `keep_rline`.
pub(crate) fn process_rline(
    slin: &mut Vec<GeomIntLine>,
    s1: &dyn Surface,
    s2: &dyn Surface,
    tol_arc: f64,
    keep_rline: bool,
) {
    let mut tol = 100.0 * tol_arc;
    if tol > 0.1 {
        tol = 0.1;
    }
    let mut i = 0;
    while i < slin.len() {
        let Some(r_src) = slin[i].as_rline() else {
            i += 1;
            continue;
        };
        let mut rline = r_src.clone();
        let Some(arc) = rline.uv_arc else {
            i += 1;
            continue;
        };
        let on_first = rline.arc_on_s1;
        let paramf = if rline.has_first_point && !rline.vertices.is_empty() {
            rline.vertices[0].param_on_line
        } else {
            rline.param_f
        };
        let paraml = if rline.has_last_point && !rline.vertices.is_empty() {
            rline.vertices[rline.vertices.len() - 1].param_on_line
        } else {
            rline.param_l
        };
        let host = if on_first { s1 } else { s2 };
        let degenerate = edge_degenerated(host, &arc, paramf, paraml);
        let mut seen: Vec<occt_core::gp::GpPnt> = Vec::new();
        let mut delete = false;
        let nlin = slin.len();
        for j in 0..nlin {
            let Some(g) = slin[j].as_gline() else {
                continue;
            };
            let kind = g.kind.clone();
            let verts: Vec<PatchPoint> = g.vertices.clone();
            if !degenerate {
                for pt in &verts {
                    if (on_first && pt.on_dom_s1) || (!on_first && pt.on_dom_s2) {
                        continue;
                    }
                    if seen.iter().any(|q| pt.p.distance(q) < tol_arc) {
                        continue;
                    }
                    let (u, v) = if on_first {
                        (pt.u1, pt.v1)
                    } else {
                        (pt.u2, pt.v2)
                    };
                    let (paramproj, p2d) = arc.project_uv(GpPnt2d::new(u, v));
                    let ptproj = host.d0(p2d.x(), p2d.y());
                    if pt.p.distance(&ptproj) <= 100.0 * tol
                        && paramproj >= paramf
                        && paramproj <= paraml
                    {
                        seen.push(pt.p);
                        let mut newpt = *pt;
                        newpt.param_on_line = paramproj;
                        rline.add_vertex(newpt);
                    }
                }
            }
            if !keep_rline && !is_rline_good(&kind, &rline, s1, s2, tol) {
                delete = true;
                break;
            }
        }
        if delete {
            slin.remove(i);
        } else {
            slin[i] = GeomIntLine::Restriction(rline);
            i += 1;
        }
    }
}

fn is_rline_good(
    gkind: &GLineKind,
    rline: &RLine,
    s1: &dyn Surface,
    s2: &dyn Surface,
    the_tol: f64,
) -> bool {
    let sq = the_tol * the_tol;
    let n = rline.vertices.len();
    if n < 2 {
        return false;
    }
    if n == 2 {
        if rline.vertices[0].p.square_distance(&rline.vertices[1].p) < sq {
            return false;
        }
        let Some(arc) = rline.uv_arc else {
            return false;
        };
        let uv = arc.value(0.5 * (arc.first + arc.last));
        let pmid = if rline.arc_on_s1 {
            s1.d0(uv.x(), uv.y())
        } else {
            s2.d0(uv.x(), uv.y())
        };
        return square_distance_to_gline(gkind, &pmid) > sq;
    }
    for v in rline.vertices.iter().skip(1).take(n.saturating_sub(2)) {
        if square_distance_to_gline(gkind, &v.p) > sq {
            return true;
        }
    }
    false
}

fn edge_degenerated(surf: &dyn Surface, arc: &RestrictionArc, paramf: f64, paraml: f64) -> bool {
    for edg in 0..=10 {
        let t = paramf + (paraml - paramf) * edg as f64 * 0.1;
        let uv = arc.value(t);
        let (_, d1u, _) = surf.d1(uv.x(), uv.y());
        if d1u.magnitude() > 1e-7 {
            return false;
        }
    }
    true
}

fn line_vertices(line: &GeomIntLine) -> &[PatchPoint] {
    match line {
        GeomIntLine::Geometric(g) => &g.vertices,
        GeomIntLine::Restriction(r) => &r.vertices,
        GeomIntLine::Analytic(a) => &a.vertices,
        GeomIntLine::Walking(w) => &w.vertices,
    }
}

fn segment_point(
    pt: PathPoint,
    s1: &dyn Surface,
    s2: &dyn Surface,
    on_first: bool,
    param: f64,
) -> PatchPoint {
    let other = if on_first { s2 } else { s1 };
    let (u_oth, v_oth) = surface_parameters(other, &pt.p).unwrap_or((0.0, 0.0));
    let (u1, v1, u2, v2) = if on_first {
        (pt.u, pt.v, u_oth, v_oth)
    } else {
        (u_oth, v_oth, pt.u, pt.v)
    };
    let mut vp = PatchPoint::new(pt.p, param, u1, v1, u2, v2);
    if on_first {
        vp.on_dom_s1 = true;
    } else {
        vp.on_dom_s2 = true;
    }
    vp
}
