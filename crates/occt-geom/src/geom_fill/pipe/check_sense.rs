//! The file-local `CheckSense` of `GeomFill_Pipe.cxx:79-222` (OCCT 8.0.0).
//!
//! It decides, for a sequence of sections, whether each section runs in the
//! same rotational sense as the first one around the first section's inertia
//! axis. When it does not (and no section is a straight line) the section is
//! reversed.
//!
//! Translation notes:
//! - `C2->Reverse()` mutates the shared `Geom_Curve` in OCCT. Curves are
//!   immutable behind `Arc` here, so the reversed curve is a new `Arc`
//!   ([`Curve::reversed`]) placed in the returned sequence; the caller's
//!   sequence is left untouched.
//! - The sample step `h` is computed once from the **first** curve
//!   (`cxx:90`) and reused for every later curve (`cxx:140-144`); that is
//!   OCCT's actual control flow and is kept as-is.

use std::sync::Arc;

use occt_core::gp::{GpAx2, GpPnt, GpVec};
use occt_core::precision::CONFUSION;

use crate::curve::Curve;
use crate::geom_lib::geom_lib_axe_of_inertia;

/// `CheckSense(Seq1, Seq2)` (`cxx:79-222`).
///
/// Returns `(no_sing, Seq2)`; `no_sing` is false when any section degenerated
/// to a straight line, and `Seq2` holds the sections with reversed copies
/// where the rotational sense disagreed.
pub fn check_sense(seq1: &[Arc<dyn Curve>]) -> (bool, Vec<Arc<dyn Curve>>) {
    let mut no_sing = true;
    let mut seq2: Vec<Arc<dyn Curve>> = Vec::new();

    let Some(c1) = seq1.first() else {
        // `Seq1.Value(1)` on an empty sequence is OCCT's `Standard_OutOfRange`.
        return (no_sing, seq2);
    };

    let f = c1.first_parameter();
    let l = c1.last_parameter();
    // `int iP, NP = 21; NCollection_Array1<gp_Pnt> Tab(1, NP)`.
    const NP: usize = 21;
    let mut tab: Vec<GpPnt> = Vec::with_capacity(NP);
    let mut u = f;
    let h = (f - l).abs() / 20.0;
    for _ in 0..NP {
        tab.push(c1.d0(u));
        u += h;
        if (u - f) * (u - l) > 0.0 {
            u = l;
        }
    }

    let (axe_ref, mut sing) = geom_lib_axe_of_inertia(&tab, CONFUSION);

    // Si la section est une droite, ca ne marche pas.
    if sing {
        no_sing = false;
    }

    let pos = axe_ref.location();
    let dir_ref = GpVec::from_xyz(axe_ref.direction().xyz());

    let mut u = (f + l - h) / 2.0 - h;
    let mut p1 = c1.d0(u);
    u += h;
    let mut p2 = c1.d0(u);
    let alpha1 = GpVec::from_pnts(&pos, &p1).angle_with_ref(
        &GpVec::from_pnts(&pos, &p2),
        &dir_ref,
    );
    p1 = p2;
    u += h;
    p2 = c1.d0(u);
    let alpha2 = GpVec::from_pnts(&pos, &p1).angle_with_ref(
        &GpVec::from_pnts(&pos, &p2),
        &dir_ref,
    );
    p1 = p2;
    u += h;
    p2 = c1.d0(u);
    let alpha3 = GpVec::from_pnts(&pos, &p1).angle_with_ref(
        &GpVec::from_pnts(&pos, &p2),
        &dir_ref,
    );
    seq2.push(c1.clone());

    for c2 in seq1.iter().skip(1) {
        let f = c2.first_parameter();
        let l = c2.last_parameter();
        let mut u = f;
        for slot in tab.iter_mut() {
            *slot = c2.d0(u);
            u += h;
            if (u - f) * (u - l) > 0.0 {
                u = l;
            }
        }
        let (axe, s): (GpAx2, bool) = geom_lib_axe_of_inertia(&tab, CONFUSION);
        sing = s;

        // Si la section est une droite, ca ne marche pas.
        if sing {
            no_sing = false;
        }

        let pos = axe.location();
        let mut u = (f + l - h) / 2.0 - h;
        let mut p1 = c2.d0(u);
        u += h;
        let mut p2 = c2.d0(u);
        let beta1 =
            GpVec::from_pnts(&pos, &p1).angle_with_ref(&GpVec::from_pnts(&pos, &p2), &dir_ref);
        p1 = p2;
        u += h;
        p2 = c2.d0(u);
        let beta2 =
            GpVec::from_pnts(&pos, &p1).angle_with_ref(&GpVec::from_pnts(&pos, &p2), &dir_ref);
        p1 = p2;
        u += h;
        p2 = c2.d0(u);
        let beta3 =
            GpVec::from_pnts(&pos, &p1).angle_with_ref(&GpVec::from_pnts(&pos, &p2), &dir_ref);

        // Meme sens ? (`cxx:169-212`).
        let pasnul1 = alpha1.abs() > CONFUSION && beta1.abs() > CONFUSION;
        let pasnul2 = alpha2.abs() > CONFUSION && beta2.abs() > CONFUSION;
        let pasnul3 = alpha3.abs() > CONFUSION && beta3.abs() > CONFUSION;
        let ok;
        if pasnul1 && pasnul2 && pasnul3 {
            if alpha1 * beta1 > 0.0 {
                ok = (alpha2 * beta2 > 0.0) || (alpha3 * beta3 > 0.0);
            } else {
                ok = (alpha2 * beta2 > 0.0) && (alpha3 * beta3 > 0.0);
            }
        } else if pasnul1 && pasnul2 && !pasnul3 {
            ok = (alpha1 * beta1 > 0.0) || (alpha2 * beta2 > 0.0);
        } else if pasnul1 && !pasnul2 && pasnul3 {
            ok = (alpha1 * beta1 > 0.0) || (alpha3 * beta3 > 0.0);
        } else if !pasnul1 && pasnul2 && pasnul3 {
            ok = (alpha2 * beta2 > 0.0) || (alpha3 * beta3 > 0.0);
        } else if pasnul1 {
            ok = alpha1 * beta1 > 0.0;
        } else if pasnul2 {
            ok = alpha2 * beta2 > 0.0;
        } else if pasnul3 {
            ok = alpha3 * beta3 > 0.0;
        } else {
            // Every candidate angle vanished: OCCT leaves `ok` true.
            ok = true;
        }

        if no_sing && !ok {
            seq2.push(Arc::from(c2.reversed()));
        } else {
            seq2.push(c2.clone());
        }
    }

    (no_sing, seq2)
}
