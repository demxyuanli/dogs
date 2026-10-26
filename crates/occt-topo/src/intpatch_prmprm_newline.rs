//! `IntPatch_PrmPrmIntersection::NewLine` (`IntPatch_PrmPrmIntersection.cxx:1470-1580`).
//!
//! Re-samples one walking line between two point indices at constant
//! curvilinear abscissa (measured on the `U1`/`V1` parameters) and snaps every
//! inserted parameter onto the true intersection with `PW.PerformFirstPoint`.

use occt_geom::Surface;

use crate::int_tools_wline::{WLine, WLineWay};
use crate::intpatch_line::IntPatchLine;
use crate::intwalk::PWalking;

/// Top-level body of `NewLine` (`cxx:1477-1579`).
///
/// `low` / `high` / `nb_pnts_to_insert` are the OCCT 1-based arguments. The
/// `Standard_OutOfRange` guard lives in the caller (`PrmPrmIntersection::new_line`).
pub(crate) fn new_line(
    s1: &dyn Surface,
    s2: &dyn Surface,
    the_line: &WLine,
    low: i32,
    high: i32,
    nb_pnts_to_insert: i32,
) -> WLine {
    let nb_pnts = nb_pnts_to_insert + high - low;

    // `U1` / `V1` / `U2` / `V2` / `AC` are `Array1(Low, High)` (`cxx:1494-1498`).
    let len = (high - low + 1) as usize;
    let mut u1 = vec![0.0f64; len];
    let mut v1 = vec![0.0f64; len];
    let mut u2 = vec![0.0f64; len];
    let mut v2 = vec![0.0f64; len];
    let mut ac = vec![0.0f64; len];
    let at = |i: i32| (i - low) as usize;

    let (a1, b1, a2, b2) = the_line.point(low).parameters();
    u1[at(low)] = a1;
    v1[at(low)] = b1;
    u2[at(low)] = a2;
    v2[at(low)] = b2;
    ac[at(low)] = 0.0;

    let mut pw = PWalking::new(s1, s2, 0.000001, 0.000001, 0.001, 0.001);

    // `cxx:1511-1526`: curvilinear abscissa on `U1` / `V1`.
    for i in (low + 1)..=high {
        let (a1, b1, a2, b2) = the_line.point(i).parameters();
        u1[at(i)] = a1;
        v1[at(i)] = b1;
        u2[at(i)] = a2;
        v2[at(i)] = b2;

        let du1 = a1 - u1[at(i - 1)];
        let dv1 = b1 - v1[at(i - 1)];
        ac[at(i)] = ac[at(i - 1)] + (du1 * du1 + dv1 * dv1).sqrt();
    }

    let mut result = Vec::new();
    result.push(*the_line.point(low));

    let ds = ac[at(high)] / f64::from(nb_pnts - 1);
    let mut indice = low;

    let dsmin = ds * 0.3;
    let smax = ac[at(high)];

    // `for (i = 2, s = ds; (i < NbPnts) && (s < smax); i++, s += ds)`
    let mut i = 2i32;
    let mut s = ds;
    while i < nb_pnts && s < smax {
        while ac[at(indice + 1)] <= s {
            result.push(*the_line.point(indice));
            indice += 1;
        }
        let a = s - ac[at(indice)];
        let b = ac[at(indice + 1)] - s;
        let nab = 1.0 / (a + b);

        if nab > ds && a > dsmin && b > dsmin {
            let start = [
                (u1[at(indice)] * b + u1[at(indice + 1)] * a) * nab,
                (v1[at(indice)] * b + v1[at(indice + 1)] * a) * nab,
                (u2[at(indice)] * b + u2[at(indice + 1)] * a) * nab,
                (v2[at(indice)] * b + v2[at(indice + 1)] * a) * nab,
            ];
            if let Some(p) = pw.perform_first_point(start) {
                result.push(p);
            }
        } else {
            s += dsmin;
        }

        i += 1;
        s += ds;
    }

    result.push(*the_line.point(high));

    let mut rwl = WLine::new();
    rwl.base = IntPatchLine::with_transitions(false, crate::intsurf_transition::IntSurfTypeTrans::Undecided, crate::intsurf_transition::IntSurfTypeTrans::Undecided);
    rwl.set_creating_way(WLineWay::PrmPrm);
    for p in result {
        rwl.add(p);
    }
    rwl
}
