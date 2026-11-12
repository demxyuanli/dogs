//! `BRepCheck_Edge::Tolerance` (`BRepCheck_Edge.cxx:598-707`): the largest
//! distance between the curve representations of an edge, sampled at
//! `NCONTROL = 23` parameters, with the 5 % margin of `cxx:705`.
//!
//! Control flow of the reference:
//!
//! 1. `nbRep = TE->Curves().Extent()`; a single representation answers
//!    `Precision::Confusion()` (`cxx:601-605`).
//! 2. `First` / `Last` come from `myHCurve` when the `Minimum` pass kept one,
//!    otherwise from `BRep_Tool::Range` (`cxx:607-616`).
//! 3. Every `Curve3D` representation on a non-degenerated edge becomes the
//!    centre adaptor: the swap at `cxx:631-636` moves it to slot 1. Every
//!    `CurveOnSurface` becomes an `Adaptor3d_CurveOnSurface`, and a
//!    `CurveOnClosedSurface` adds a second one (`cxx:653-666`, `nbRep++`). Any
//!    other representation kind decrements `nbRep` (`cxx:668-670`).
//! 4. The samples of `cxx:673-704` take the maximum square distance from slot 1
//!    to every other slot; an infinite coordinate aborts with
//!    `Precision::Infinite()`.
//!
//! Translation notes:
//! - `theRep(i)` is an `Adaptor3d_Curve` or an `Adaptor3d_CurveOnSurface`;
//!   [`TolRep`] is that union, evaluating through `surface->Value(2d)` for the
//!   surface case (`Adaptor3d_CurveOnSurface::Value`, `cxx` via
//!   `GeomAdaptor_Surface::Value`).
//! - `myHCurve` is `None` here: the port's `Minimum` pass does not retain an
//!   adaptor, so `First`/`Last` always come from the edge range.
//! - The first pcurve uses `edge.Location() * representation.Location()`
//!   (`cxx:644`); the second pcurve of a closed representation uses the
//!   representation location alone (`cxx:657-658`).
//! - Representations are keyed by their `Geom_Surface` handle, which is what
//!   `BRep_CurveRepresentation::IsCurveOnSurface(S, l)` compares as well; the
//!   visiting order of distinct surfaces follows the registry's pcurve map, so
//!   when the edge has no 3D curve the role of slot 1 can differ from OCCT's
//!   list order (only then: the swap above pins slot 1 to the 3D curve).

use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_core::precision::{Precision, CONFUSION, INFINITE};
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_geom2d::Curve2d;

use crate::shape::Edge;
use crate::tgeometry::GeometryRegistry;

/// `NCONTROL` (`BRepCheck_Edge.cxx:61`).
const NCONTROL: usize = 23;

/// One entry of the reference's `theRep` array (`cxx:606`): an
/// `Adaptor3d_Curve` or an `Adaptor3d_CurveOnSurface`.
enum TolRep {
    /// `GeomAdaptor_Curve` over the `Curve3D` representation.
    Curve3d(Arc<dyn Curve>),
    /// `Adaptor3d_CurveOnSurface` over one pcurve and its surface.
    CurveOnSurface(Arc<dyn Curve2d>, Arc<dyn Surface>),
}

impl TolRep {
    /// `Adaptor3d_Curve::Value` / `Adaptor3d_CurveOnSurface::Value`.
    fn value(&self, prm: f64) -> GpPnt {
        match self {
            Self::Curve3d(c) => c.d0(prm),
            Self::CurveOnSurface(pc, s) => {
                let uv = pc.d0(prm);
                s.d0(uv.x(), uv.y())
            }
        }
    }
}

/// `BRepCheck_Edge::Tolerance` (`cxx:598-707`).
pub fn edge_tolerance(e: &Edge) -> f64 {
    let reg = GeometryRegistry::global();
    let curve3d = reg.edge_curve(&e.0);
    let degenerated = reg.is_degenerated_edge(&e.0);
    let reps = reg.edge_pcurve_reps_with_location(&e.0);

    // `TE->Curves().Extent()` (`cxx:601`): one 3D representation plus one per
    // `CurveOnSurface`. A closed-surface representation counts once, which is
    // why its two pcurves are grouped by their surface below.
    let nb_reps_initial = usize::from(curve3d.is_some()) + count_surface_reps(&reps);
    if nb_reps_initial <= 1 {
        return CONFUSION;
    }

    // `cxx:607-616`: `myHCurve` is not retained by the port's `Minimum`.
    let (first, last) = reg.edge_parameters(&e.0);

    // `cxx:617-671`: build `theRep`, applying the slot swap of the 3D curve.
    let mut nb_rep = nb_reps_initial;
    let mut the_rep: Vec<TolRep> = Vec::with_capacity(nb_reps_initial * 2);
    let loc = e.0.location();

    if let Some(c3d) = &curve3d {
        if !degenerated {
            // `cxx:626-637`.
            let transformed = if loc.is_identity() {
                c3d.clone()
            } else {
                Arc::from(c3d.transformed(&loc.transformation()))
            };
            push_centre(&mut the_rep, TolRep::Curve3d(transformed));
        } else {
            // `cxx:668-670`: a degenerate edge's 3D representation is neither.
            nb_rep = nb_rep.saturating_sub(1);
        }
    }

    // `cxx:639-666`: one `CurveOnSurface` per (surface, first pcurve), plus the
    // second pcurve of a closed-surface representation.
    let mut i = 0;
    while i < reps.len() {
        let (surf, pc1, _, _, rep_loc) = &reps[i];
        // The closed-surface second pcurve is the next entry on the same
        // surface (`GeometryRegistry::edge_pcurves` keeps forward-then-reversed
        // order).
        let has_second = reps
            .get(i + 1)
            .is_some_and(|(s2, _, _, _, l2)| same_surface(s2, surf) && l2.is_equal(rep_loc));
        let first_loc = loc.transformation().multiplied(&rep_loc.transformation());
        let sref = if first_loc.form() == occt_core::gp::TrsfForm::Identity {
            surf.clone()
        } else {
            Arc::from(surf.transformed(&first_loc))
        };
        the_rep.push(TolRep::CurveOnSurface(pc1.clone(), sref));
        if has_second {
            // `cxx:653-666`: `PCurve2` with the representation location only.
            let second_surface = if rep_loc.is_identity() {
                surf.clone()
            } else {
                Arc::from(surf.transformed(&rep_loc.transformation()))
            };
            the_rep.push(TolRep::CurveOnSurface(
                reps[i + 1].1.clone(),
                second_surface,
            ));
            nb_rep += 1;
            i += 2;
        } else {
            i += 1;
        }
    }

    // `cxx:673-704`.
    let mut tol_cal = 0.0_f64;
    for i in 0..NCONTROL {
        let prm = ((NCONTROL - 1 - i) as f64 * first + i as f64 * last) / (NCONTROL - 1) as f64;
        let Some(centre) = the_rep.first() else {
            return CONFUSION;
        };
        let centre_p = centre.value(prm);
        if is_infinite_pnt(&centre_p) {
            return INFINITE;
        }
        for rep in the_rep.iter().skip(1).take(nb_rep.saturating_sub(1)) {
            let oth_p = rep.value(prm);
            if is_infinite_pnt(&oth_p) {
                return INFINITE;
            }
            let dist2 = centre_p.square_distance(&oth_p);
            if dist2 > tol_cal {
                tol_cal = dist2;
            }
        }
    }
    // On prend 5% de marge car au dessus on controle severement (`cxx:705`).
    tol_cal.sqrt() * 1.05
}

/// `cxx:631-636`: `it = iRep; if (iRep > 1) { theRep(iRep) = theRep(1); it = 1; }
/// theRep(it) = <centre>` — the centre adaptor always ends up in slot 1 and
/// whatever was there moves to the centre's old slot.
fn push_centre(the_rep: &mut Vec<TolRep>, centre: TolRep) {
    if the_rep.is_empty() {
        the_rep.push(centre);
        return;
    }
    let old_slot1 = the_rep.remove(0);
    the_rep.push(old_slot1);
    the_rep.insert(0, centre);
}

/// Number of `CurveOnSurface` representations: one per distinct surface, two
/// pcurves on the same surface still being one closed-surface representation.
fn count_surface_reps(
    reps: &[(
        Arc<dyn Surface>,
        Arc<dyn Curve2d>,
        f64,
        f64,
        occt_core::toploc::TopLocLocation,
    )],
) -> usize {
    let mut count = 0;
    let mut i = 0;
    while i < reps.len() {
        count += 1;
        let closed = reps.get(i + 1).is_some_and(|(s2, _, _, _, location)| {
            same_surface(s2, &reps[i].0) && location.is_equal(&reps[i].4)
        });
        i += if closed { 2 } else { 1 };
    }
    count
}

/// `Handle(Geom_Surface)` identity, the key of
/// `BRep_CurveRepresentation::IsCurveOnSurface(S, l)`.
fn same_surface(a: &Arc<dyn Surface>, b: &Arc<dyn Surface>) -> bool {
    Arc::as_ptr(a) as *const () as usize == Arc::as_ptr(b) as *const () as usize
}

/// `Precision::IsInfinite` on all three coordinates (`cxx:681-684`,
/// `cxx:688-692`).
fn is_infinite_pnt(p: &GpPnt) -> bool {
    Precision::is_infinite(p.x()) || Precision::is_infinite(p.y()) || Precision::is_infinite(p.z())
}
