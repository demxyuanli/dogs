//! `IsPointOnLine` / `DublicateOfLinesProcessing` / `AddWLine`
//! (`IntPatch_PrmPrmIntersection.cxx:278-312`, `:3993-4143`).
//!
//! These three file-static helpers of `IntPatch_PrmPrmIntersection.cxx` decide
//! whether a freshly walked `IntWalk_PWalking` line duplicates an already
//! collected `IntPatch_WLine`, and remove the loser.

use occt_core::gp::{GpDir, GpPnt2d, GpVec};

use crate::geom_int::GeomIntLine;
use crate::int_tools_wline::{PntOn2S, WLine};

/// `IsPointOnLine(thePOn2S, theWLine, Deflection)` (`cxx:3993-4081`).
///
/// The three `IsOut*Box` guards (`cxx:4004-4006`) come first, then the
/// segment walk. The component-wise `HM.X() < Deflection` chain at
/// `cxx:4026-4052` is OCCT's own test, not a rewritten distance test.
pub(crate) fn is_point_on_line(p: &PntOn2S, wl: &WLine, deflection: f64) -> bool {
    let deflection2 = deflection * deflection;
    let (pu1, pv1, pu2, pv2) = p.parameters();

    if wl.is_out_surf1_box(&GpPnt2d::new(pu1, pv1)) {
        return false;
    }
    if wl.is_out_surf2_box(&GpPnt2d::new(pu2, pv2)) {
        return false;
    }
    if wl.is_out_box(&p.value()) {
        return false;
    }

    let nb = wl.nb_pnts();
    let p_start = p.value();
    let mut is_on_line = false;
    for ll in 1..nb {
        let pa = wl.point(ll).value();
        let pb = wl.point(ll + 1).value();
        let am = GpVec::from_pnts(&pa, &p_start);
        let mb = GpVec::from_pnts(&p_start, &pb);
        let ammb = am.dot(&mb);

        if ammb > 0.0 {
            // `gp_Dir ABN(Pb - Pa)` (`cxx:4020`). `AM.MB > 0` implies
            // `Pa != Pb`, so the null-vector construction cannot happen here.
            let Ok(abn) = GpDir::new(pb.x() - pa.x(), pb.y() - pa.y(), pb.z() - pa.z()) else {
                return false;
            };
            let lan = abn.x() * am.x() + abn.y() * am.y() + abn.z() * am.z();
            let ah = GpVec::from_xyz(&occt_core::gp::GpXyz::new(
                lan * abn.x(),
                lan * abn.y(),
                lan * abn.z(),
            ));
            let hm = am.subtracted(&ah);
            let mut d = 0.0;
            if hm.x() < deflection {
                d += hm.x() * hm.x();
                if hm.y() < deflection {
                    d += hm.y() * hm.y();
                    if hm.z() < deflection {
                        d += hm.z() * hm.z();
                    } else {
                        d = deflection2;
                    }
                } else {
                    d = deflection2;
                }
            } else {
                d = deflection2;
            }
            if d < deflection2 {
                is_on_line = true;
            }
        } else {
            let dab = pa.square_distance(&pb);
            let dap = pa.square_distance(&p_start);
            if dap < dab {
                is_on_line = true;
            } else {
                let dbp = pb.square_distance(&p_start);
                if dbp < dab {
                    is_on_line = true;
                }
            }
        }
    }
    is_on_line
}

/// `DublicateOfLinesProcessing(thePW, theWLID, theLines, theIsRejectReq)`
/// (`cxx:278-312`).
///
/// `theWLID` is 1-based, as in OCCT's `theLines.Value(theWLID)`.
pub(crate) fn duplicate_of_lines_processing(
    pw: &crate::intwalk::PWalking,
    wl_id: usize,
    lines: &mut Vec<GeomIntLine>,
    is_reject_req: &mut bool,
) {
    if wl_id < 1 || wl_id > lines.len() {
        return;
    }
    let Some(exist) = lines[wl_id - 1].as_wline() else {
        return;
    };
    let nb_prev = exist.nb_pnts();
    let nb_curr = pw.nb_points();

    if nb_prev < nb_curr {
        // Remove the previous line.
        lines.remove(wl_id - 1);
        *is_reject_req = false;
    } else if nb_prev == nb_curr {
        let mut l_prev = 0.0;
        let mut l_curr = 0.0;
        for nb_pp in 1..nb_prev {
            let p1prev = exist.point(nb_pp).value();
            let p2prev = exist.point(nb_pp + 1).value();
            let p1curr = pw.value(nb_pp).value();
            let p2curr = pw.value(nb_pp + 1).value();
            l_prev += p1prev.distance(&p2prev);
            l_curr += p1curr.distance(&p2curr);
        }
        if l_prev < l_curr {
            lines.remove(wl_id - 1);
            *is_reject_req = false;
        }
    }
}

/// `AddWLine(theLines, theWLine, Deflection)` (`cxx:4085-4143`).
///
/// Drops every already collected line whose middle point and vertices all lie
/// on `wline`, then appends `wline`.
pub(crate) fn add_wline(lines: &mut Vec<GeomIntLine>, wline: WLine, deflection: f64) {
    let mut i = 1usize;
    while i <= lines.len() {
        let mut is_to_remove = false;
        if let Some(wl) = lines[i - 1].as_wline() {
            let mid = (wl.nb_pnts() + 1) / 2;
            if mid > 0 {
                let pnt = wl.point(mid);
                // `IsPointOnLine(aPnt, theWLine, Deflection)` (`cxx:4109`).
                if is_point_on_line(pnt, &wline, deflection) {
                    is_to_remove = true;
                    for j in 1..=wl.nb_vertex() {
                        if !is_point_on_line(&wl.vertex(j).pnt_on_2s(), &wline, deflection) {
                            is_to_remove = false;
                            break;
                        }
                    }
                }
            }
        }
        if is_to_remove {
            lines.remove(i - 1);
        } else {
            i += 1;
        }
    }
    lines.push(GeomIntLine::Walking(wline));
}
