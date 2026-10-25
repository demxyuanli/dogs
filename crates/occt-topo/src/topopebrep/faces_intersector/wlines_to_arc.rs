//! `TestWLinesToAnArc` and `MergeWLinesIfAllSegmentsAlongRestriction`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRep/TopOpeBRep_FacesIntersector.cxx:1060-1348`.
//!
//! Translation notes:
//! - OCCT identifies the queued walking lines by `Handle` identity when the
//!   chain is replaced (`cxx:1191-1200`). The port's `GeomIntLine` owns its
//!   payload, so the chain records the indices into `slin` instead; the indices
//!   are removed in descending order, which deletes the same elements the
//!   handle comparison finds.
//! - `slin.Remove(j)` shifts the sequence but the outer index `i` is not
//!   corrected in the source (`cxx:1196`), so the port keeps a plain stale
//!   index loop as well.

use occt_core::gp::GpPnt;
use occt_geom::Surface;

use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::geom_int::GeomIntLine;

use super::merged_wline::{get_arc, get_merged_wline_on_restriction};
use super::rline_build::{build_rline, test_wline_along_restriction};

/// Whether `wline` may join the chain: it either already carries the arc of
/// side `rank`, or all its samples lie on that side's restriction.
fn append_wline_on_rank(
    wline: &crate::int_tools_wline::WLine,
    rank: i32,
    surface1: &dyn Surface,
    domain1: &mut TopolTool,
    surface2: &dyn Surface,
    domain2: &mut TopolTool,
    tol_arc: f64,
) -> bool {
    if rank == 1 {
        if !wline.has_arc_on_s1()
            && !test_wline_along_restriction(wline, rank, Some(surface1), domain1, tol_arc)
        {
            return false;
        }
    }
    if rank == 2
        && !wline.has_arc_on_s2()
        && !test_wline_along_restriction(wline, rank, Some(surface2), domain2, tol_arc)
    {
        return false;
    }
    true
}

/// `TestWLinesToAnArc(Slin, Surface1, Domain1, Surface2, Domain2, TolArc)`
/// (`cxx:1064-1210`): replace each chain of walking lines that runs along one
/// restriction arc with a single restriction line.
pub fn test_wlines_to_an_arc(
    slin: &mut Vec<GeomIntLine>,
    surface1: &dyn Surface,
    domain1: &mut TopolTool,
    surface2: &dyn Surface,
    domain2: &mut TopolTool,
    tol_arc: f64,
) {
    let mut seq_of_rline: Vec<GeomIntLine> = Vec::new();

    for rank in 1..=2 {
        // Indices into `slin`, in append order (`aSeqOfWLine`).
        let mut seq_of_wline: Vec<usize> = Vec::new();

        // The index is deliberately not corrected after a removal: the source
        // loop keeps using its own `i` over the shortened sequence.
        let mut i = 0usize;
        while i < slin.len() {
            let Some(wline) = slin[i].as_wline() else {
                i += 1;
                continue;
            };

            let nbvtx = wline.nb_vertex();
            let vtx1 = wline.vertex(1);
            let vtx2 = wline.vertex(nbvtx);

            let isvertex = if rank == 1 {
                vtx1.is_vertex_on_s1()
            } else {
                vtx1.is_vertex_on_s2()
            };

            let mut wline_was_appended = false;

            if isvertex {
                let appendwline = append_wline_on_rank(
                    wline, rank, surface1, domain1, surface2, domain2, tol_arc,
                );
                wline_was_appended = appendwline;
                if appendwline {
                    seq_of_wline.push(i);
                }
            } else {
                if seq_of_wline.is_empty() {
                    i += 1;
                    continue;
                }
                let a_last_wline = slin[*seq_of_wline.last().unwrap()]
                    .as_wline()
                    .expect("chain holds a walking line");
                let a_last_point = a_last_wline.vertex(a_last_wline.nb_vertex());
                let a_tolerance = if a_last_point.tolerance > vtx1.tolerance {
                    a_last_point.tolerance
                } else {
                    vtx1.tolerance
                };
                if a_last_point.p.is_equal_with_tolerance(&vtx1.p, a_tolerance) {
                    let appendwline = append_wline_on_rank(
                        wline, rank, surface1, domain1, surface2, domain2, tol_arc,
                    );
                    wline_was_appended = appendwline;
                    if appendwline {
                        seq_of_wline.push(i);
                    }
                } else {
                    seq_of_wline.clear();
                }
            }

            let isvertex_last = if rank == 1 {
                vtx2.is_vertex_on_s1()
            } else {
                vtx2.is_vertex_on_s2()
            };

            if wline_was_appended && isvertex_last && !seq_of_wline.is_empty() {
                // `cxx:1176-1187`: build the restriction line from the chain.
                let chain: Vec<crate::int_tools_wline::WLine> = seq_of_wline
                    .iter()
                    .map(|&k| slin[k].as_wline().expect("chain entry").clone())
                    .collect();
                let an_rline = if rank == 1 {
                    build_rline(&chain, rank, Some(surface1), domain1, tol_arc)
                } else {
                    build_rline(&chain, rank, Some(surface2), domain2, tol_arc)
                };

                if let Some(an_rline) = an_rline {
                    seq_of_rline.push(GeomIntLine::Restriction(an_rline));
                    // `cxx:1190-1201`: drop the chain members from `slin`.
                    let mut doomed = seq_of_wline.clone();
                    doomed.sort_unstable();
                    for idx in doomed.into_iter().rev() {
                        if idx < slin.len() {
                            slin.remove(idx);
                        }
                    }
                }
                seq_of_wline.clear();
            }

            i += 1;
        }
    }

    // `cxx:1206-1209`
    for rline in seq_of_rline {
        slin.push(rline);
    }
}

/// `MergeWLinesIfAllSegmentsAlongRestriction(Slin, Surface1, Domain1,
/// Surface2, Domain2, TolArc)` (`cxx:1226-1344`): when every walking line lies
/// on one restriction, and no restriction line can carry them, replace the
/// whole sequence with one merged walking line sampled from that arc.
pub fn merge_wlines_if_all_segments_along_restriction(
    slin: &mut Vec<GeomIntLine>,
    surface1: &dyn Surface,
    domain1: &mut TopolTool,
    surface2: &dyn Surface,
    domain2: &mut TopolTool,
    tol_arc: f64,
) {
    let tol = 1.0e-9;

    // `cxx:1234-1277`: are all walking lines placed along one restriction?
    let mut wl_on_rs1 = 0;
    let mut wl_on_rs2 = 0;
    let mut nb_wlines = 0;
    let mut sq_vertex_points: Vec<GpPnt> = Vec::new();

    for rank in 1..=2 {
        nb_wlines = 0;
        for line in slin.iter() {
            let Some(wline) = line.as_wline() else {
                continue;
            };
            nb_wlines += 1;
            let nbvtx = wline.nb_vertex();
            let vtx1 = wline.vertex(1);
            let vtx2 = wline.vertex(nbvtx);
            if rank == 1 {
                sq_vertex_points.push(vtx1.p);
                sq_vertex_points.push(vtx2.p);
                if test_wline_along_restriction(
                    wline,
                    rank,
                    Some(surface1),
                    domain1,
                    tol_arc,
                ) {
                    wl_on_rs1 += 1;
                }
            } else if test_wline_along_restriction(wline, rank, Some(surface2), domain2, tol_arc)
            {
                wl_on_rs2 += 1;
            }
        }
        if nb_wlines == wl_on_rs1 || nb_wlines == wl_on_rs2 {
            break;
        }
    }

    // `cxx:1279-1293`: which side carries the merge, if any.
    let wline_rank = if wl_on_rs1 == nb_wlines {
        1
    } else if wl_on_rs2 == nb_wlines {
        2
    } else {
        return;
    };

    // `cxx:1295-1301`: avoid closed (degenerated) edges.
    if sq_vertex_points.len() <= 2 {
        return;
    }
    if sq_vertex_points[0]
        .is_equal_with_tolerance(&sq_vertex_points[sq_vertex_points.len() - 1], tol)
    {
        return;
    }

    let mut tol_vrtx = 1.0e-5;
    let test_point_index = if sq_vertex_points.len() > 3 {
        sq_vertex_points.len() / 2
    } else {
        2
    };
    let test_point = sq_vertex_points[test_point_index - 1];
    let mut fp = 0.0;
    let mut lp = 0.0;
    let mut line_on2_s: Vec<crate::int_tools_wline::PntOn2S> = Vec::new();

    // `cxx:1310-1334`: side 1 merges on Surface1 / Domain1 and classifies
    // against Surface2; side 2 is the mirror image.
    let arcnumber = if wline_rank == 1 {
        get_arc(
            slin,
            wline_rank,
            surface1,
            domain1,
            surface2,
            &test_point,
            &mut tol_vrtx,
            &mut line_on2_s,
            &mut fp,
            &mut lp,
        )
    } else {
        get_arc(
            slin,
            wline_rank,
            surface2,
            domain2,
            surface1,
            &test_point,
            &mut tol_vrtx,
            &mut line_on2_s,
            &mut fp,
            &mut lp,
        )
    };

    if arcnumber == 0 {
        return;
    }

    // `cxx:1337-1343`
    if let Some(an_wline) = get_merged_wline_on_restriction(slin, tol_vrtx, &line_on2_s) {
        slin.clear();
        slin.push(GeomIntLine::Walking(an_wline));
    }
}
