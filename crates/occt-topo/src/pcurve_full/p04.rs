//! `ShapeAnalysis_Surface` singularity machinery and the
//! `ShapeConstruct_ProjectCurveOnSurface` degenerate-point chain built on it:
//! `projectDegeneratedPoints` + `correctExtremity` (`cxx:1466-1497` call sites,
//! `cxx:1910-2037` implementation) and the `Handle AdjustOverDegen` block
//! (`cxx:1746-1890`, on by default: `myAdjustOverDegen = 1`, `cxx:644`).
//!
//! Both run on `thePoints2d` *before* the already-ported U/V-closed blocks
//! (`cxx:1466` / `cxx:1746` vs `cxx:1499`), so they change the input those
//! blocks consume.
//!
//! PARKED (reported, not fudged):
//! - `ShapeAnalysis_Surface::myGap`: a side effect of `IsDegenerated`
//!   (`ShapeAnalysis_Surface.cxx:362`) / `ProjectDegenerated`
//!   (`ShapeAnalysis_Surface.cxx:449`, `:498`) used only as the
//!   `Precision::Confusion() + gap` seed of `NextValueOfUV`. `NextValueOfUV`
//!   is now ported (`p01::next_value_of_uv`) and its callers pass
//!   `Precision::Confusion() + gap` / `1000 * gap` explicitly
//!   (`ShapeConstruct_ProjectCurveOnSurface.cxx:1358`, `:1413`, `:1432`); the
//!   gap here is recovered as the residue of the returned point
//!   (`p01::value_of_uv_with_gap`).
//! - The singularity arrays are recomputed per pcurve instead of being cached in
//!   `ShapeAnalysis_Surface::myNbDeg` / `myPreci` / ... (`ShapeAnalysis_Surface.hxx`).
//!   `ComputeSingularities` is a pure function of the surface and its bounds, so
//!   the values each block sees are identical.
//! - `ShapeAnalysis_Surface::ProjectDegenerated(P3d, preci, neighbour, result)`
//!   (`ShapeAnalysis_Surface.cxx:414-458`) and `DegeneratedValues`
//!   (`ShapeAnalysis_Surface.cxx:373-411`): no caller in the
//!   `ShapeConstruct_ProjectCurveOnSurface` chain (the free function
//!   `projectDegeneratedPoints` at `cxx:600-635` calls the *sequence* overload).

use super::p01::{sa_is_u_closed, sa_is_v_closed, value_of_uv};
use super::prelude::*;
use occt_core::precision::{CONFUSION, PCONFUSION};

/// One entry of `ShapeAnalysis_Surface::myP3d` / `myPreci` / `myFirstP2d` /
/// `myLastP2d` / `myFirstPar` / `myLastPar` / `myUIsoDeg`
/// (`ShapeAnalysis_Surface.hxx:296-303`). `preci` is the 3D distance from the
/// singular point within which a point counts as degenerated.
#[derive(Debug, Clone, Copy)]
pub(super) struct Singularity {
    pub preci: f64,
    pub p3d: GpPnt,
    pub first_p2d: GpPnt2d,
    pub last_p2d: GpPnt2d,
    pub first_par: f64,
    pub last_par: f64,
    pub u_iso_deg: bool,
}

/// OCCT's `Coord(1)` / `Coord(2)`: 1 is U (the `x` coordinate), 2 is V.
fn coord_of(p: &GpPnt2d, idx: usize) -> f64 {
    if idx == 1 {
        p.x()
    } else {
        p.y()
    }
}

/// `ShapeAnalysis_Surface::SortSingularities`
/// (`ShapeAnalysis_Surface.cxx:1845-1883`): selection sort of the parallel
/// arrays by ascending `myPreci`, so every later query can stop at the first
/// entry whose `preci` exceeds its tolerance.
fn sort_singularities(v: &mut [Singularity]) {
    let n = v.len();
    for i in 0..n.saturating_sub(1) {
        let mut min_preci = v[i].preci;
        let mut min_index = i;
        for j in i + 1..n {
            if min_preci > v[j].preci {
                min_preci = v[j].preci;
                min_index = j;
            }
        }
        if min_index != i {
            v.swap(i, min_index);
        }
    }
}

/// `ShapeAnalysis_Surface::ComputeSingularities`
/// (`ShapeAnalysis_Surface.cxx:181-293`). `ShapeAnalysis_Surface::Bounds`
/// (`ShapeAnalysis_Surface.lxx:54-64`) is the surface's own parameter window,
/// i.e. [`Surface::u_range`] / [`Surface::v_range`].
pub(super) fn compute_singularities(s: &dyn Surface) -> Vec<Singularity> {
    let (su1, su2) = s.u_range();
    let (sv1, sv2) = s.v_range();
    let mut out: Vec<Singularity> = Vec::new();
    if let Some(cone) = s.gp_cone() {
        // `cxx:202-212`.
        let v_apex = -cone.radius() / cone.semi_angle().sin();
        out.push(Singularity {
            preci: 0.0,
            p3d: cone.apex(),
            first_p2d: GpPnt2d::new(su1, v_apex),
            last_p2d: GpPnt2d::new(su2, v_apex),
            first_par: su1,
            last_par: su2,
            u_iso_deg: false,
        });
    } else if let Some(torus) = s.gp_torus() {
        // `cxx:214-232`.
        let minor_r = torus.minor_radius();
        let major_r = torus.major_radius();
        let ang = (major_r / minor_r).min(1.0).acos();
        let preci = (major_r - minor_r).max(0.0);
        out.push(Singularity {
            preci,
            p3d: s.d0(0.0, PI - ang),
            first_p2d: GpPnt2d::new(su1, PI - ang),
            last_p2d: GpPnt2d::new(su2, PI - ang),
            first_par: su1,
            last_par: su2,
            u_iso_deg: false,
        });
        // `myNbDeg = (majorR > minorR ? 1 : 2)`: the second entry is only kept
        // for a horn/lemon torus (`cxx:223-232`).
        if major_r <= minor_r {
            out.push(Singularity {
                preci,
                p3d: s.d0(0.0, PI + ang),
                first_p2d: GpPnt2d::new(su2, PI + ang),
                last_p2d: GpPnt2d::new(su1, PI + ang),
                first_par: su1,
                last_par: su2,
                u_iso_deg: false,
            });
        }
    } else if s.gp_sphere().is_some() {
        // `cxx:234-246`: `myP3d[0]` is the pole at `sv2`, `myP3d[1]` at `sv1`.
        out.push(Singularity {
            preci: 0.0,
            p3d: s.d0(su1, sv2),
            first_p2d: GpPnt2d::new(su2, sv2),
            last_p2d: GpPnt2d::new(su1, sv2),
            first_par: su1,
            last_par: su2,
            u_iso_deg: false,
        });
        out.push(Singularity {
            preci: 0.0,
            p3d: s.d0(su1, sv1),
            first_p2d: GpPnt2d::new(su1, sv1),
            last_p2d: GpPnt2d::new(su2, sv1),
            first_par: su1,
            last_par: su2,
            u_iso_deg: false,
        });
    } else if is_bounded_surface_kind(s) {
        // `cxx:248-291`: a surface with no known singularity gets four
        // synthetic entries at the corners of the parameter window
        // (`Geom_RectangularTrimmedSurface`, `Geom_BoundedSurface` ->
        // `Geom_BSplineSurface` / `Geom_BezierSurface`,
        // `Geom_SurfaceOfRevolution`, `Geom_OffsetSurface`). Each `preci` is the
        // largest 3D distance from the window's edge midpoint to the corners of
        // that edge, so an edge that collapses to a point yields `preci == 0`
        // and is the only case `NbSingularities(myPreci) > 0` sees.
        let mid_u = 0.5 * (su1 + su2);
        let mid_v = 0.5 * (sv1 + sv2);
        let p3d = [
            s.d0(su1, mid_v),
            s.d0(su2, mid_v),
            s.d0(mid_u, sv1),
            s.d0(mid_u, sv2),
        ];
        let first_p2d = [
            GpPnt2d::new(su1, sv2),
            GpPnt2d::new(su2, sv1),
            GpPnt2d::new(su1, sv1),
            GpPnt2d::new(su2, sv2),
        ];
        let last_p2d = [
            GpPnt2d::new(su1, sv1),
            GpPnt2d::new(su2, sv2),
            GpPnt2d::new(su2, sv1),
            GpPnt2d::new(su1, sv2),
        ];
        let c1 = s.d0(su1, sv1);
        let c2 = s.d0(su1, sv2);
        let c3 = s.d0(su2, sv1);
        let c4 = s.d0(su2, sv2);
        let preci = [
            c1.distance(&c2).max(p3d[0].distance(&c1)).max(p3d[0].distance(&c2)),
            c3.distance(&c4).max(p3d[1].distance(&c3)).max(p3d[1].distance(&c4)),
            c1.distance(&c3).max(p3d[2].distance(&c1)).max(p3d[2].distance(&c3)),
            c2.distance(&c4).max(p3d[3].distance(&c2)).max(p3d[3].distance(&c4)),
        ];
        let first_par = [sv1, sv1, su1, su1];
        let last_par = [sv2, sv2, su2, su2];
        let u_iso_deg = [true, true, false, false];
        for i in 0..4 {
            out.push(Singularity {
                preci: preci[i],
                p3d: p3d[i],
                first_p2d: first_p2d[i],
                last_p2d: last_p2d[i],
                first_par: first_par[i],
                last_par: last_par[i],
                u_iso_deg: u_iso_deg[i],
            });
        }
    }
    sort_singularities(&mut out);
    out
}

/// `mySurf->Surface()->IsKind(STANDARD_TYPE(Geom_BoundedSurface)) ||
/// IsKind(Geom_SurfaceOfRevolution) || IsKind(Geom_OffsetSurface)`
/// (`ShapeAnalysis_Surface.cxx:248-249`).
///
/// `Geom_BoundedSurface` (`Geom_BoundedSurface.hxx:50`) only has
/// `Geom_BSplineSurface` (`Geom_BSplineSurface.hxx:150`),
/// `Geom_BezierSurface` (`Geom_BezierSurface.hxx:105`) and
/// `Geom_RectangularTrimmedSurface`
/// (`Geom_RectangularTrimmedSurface.hxx:50`) as subclasses - a
/// `Geom_SurfaceOfRevolution` / `Geom_SurfaceOfLinearExtrusion` derives from
/// `Geom_SweptSurface : Geom_Surface` (`Geom_SweptSurface.hxx:33`) and a
/// `Geom_OffsetSurface` from `Geom_Surface` (`Geom_OffsetSurface.hxx:62`), so
/// both are queried on their own. `SurfaceOfLinearExtrusion` is *not* in the
/// OCCT condition and gets no singularity here either.
///
/// `osculating_bspline` is `Some` exactly for `Geom_BSplineSurface` and
/// `Geom_BezierSurface` and for nothing else in this port, so it is the
/// `Geom_BezierSurface` arm here.
fn is_bounded_surface_kind(s: &dyn Surface) -> bool {
    s.rectangular_trimmed_basis().is_some()
        || s.is_bspline_surface()
        || s.osculating_bspline().is_some()
        || s.is_surface_of_revolution()
        || s.is_offset_surface()
}

/// `ShapeAnalysis_Surface::NbSingularities(preci)`
/// (`ShapeAnalysis_Surface.cxx:305-320`).
pub(super) fn nb_singularities(sing: &[Singularity], preci: f64) -> usize {
    sing.iter().filter(|s| s.preci <= preci).count()
}

/// `ShapeAnalysis_Surface::IsDegenerated(P3d, preci)`
/// (`ShapeAnalysis_Surface.cxx:354-370`). The `myGap` side effect is parked
/// (see the module header).
pub(super) fn is_degenerated(sing: &[Singularity], p3d: &GpPnt, preci: f64) -> bool {
    for s in sing {
        if s.preci > preci {
            break;
        }
        if s.p3d.distance(p3d) <= preci {
            return true;
        }
    }
    false
}

/// The `ShapeAnalysis_Surface` singularity queries of
/// `ShapeAnalysis_Wire::CheckDegenerated` (`ShapeAnalysis_Wire.cxx:896-1116`),
/// over one `ComputeSingularities` result (OCCT caches the parallel arrays in
/// `ShapeAnalysis_Surface::myNbDeg` / `myPreci` / `myP3d` / `myFirstP2d` /
/// `myLastP2d` / `myFirstPar` / `myLastPar`; the detector is built once per
/// face, so the values each query sees are identical).
pub(crate) struct SurfaceSingularities {
    entries: Vec<Singularity>,
}

impl SurfaceSingularities {
    /// `ShapeAnalysis_Surface::ComputeSingularities`
    /// (`ShapeAnalysis_Surface.cxx:181-293`).
    pub(crate) fn compute(s: &dyn Surface) -> Self {
        Self { entries: compute_singularities(s) }
    }

    /// `ShapeAnalysis_Surface::DegeneratedValues(P3d, preci, firstP2d,
    /// lastP2d, firstPar, lastPar, forward)` (`ShapeAnalysis_Surface.cxx:373-411`)
    /// and the `NbSingularities(preci)` + `Singularity(i, ...)` minimum-gap
    /// search at `ShapeAnalysis_Wire.cxx:1013-1031`: both return the singular
    /// entry with the smallest gap to `p3d`, searched in ascending `preci`
    /// order (`SortSingularities`, `cxx:1845-1883`).
    pub(crate) fn min_gap(
        &self,
        p3d: &GpPnt,
        preci: f64,
    ) -> Option<(GpPnt2d, GpPnt2d, f64, f64)> {
        let mut best: Option<(f64, usize)> = None;
        for (i, s) in self.entries.iter().enumerate() {
            if s.preci > preci {
                break;
            }
            let gap = s.p3d.distance(p3d);
            if gap <= preci && best.is_none_or(|(g, _)| gap < g) {
                best = Some((gap, i));
            }
        }
        best.map(|(_, i)| {
            let s = &self.entries[i];
            (s.first_p2d, s.last_p2d, s.first_par, s.last_par)
        })
    }

    /// `ShapeAnalysis_Surface::IsDegenerated(P3d, preci)`
    /// (`ShapeAnalysis_Surface.cxx:354-370`).
    pub(crate) fn is_degenerated(&self, p3d: &GpPnt, preci: f64) -> bool {
        is_degenerated(&self.entries, p3d, preci)
    }
}

/// `ShapeAnalysis_Surface::ProjectDegenerated(nbrPnt, points, pnt2d, preci,
/// direct)` (`ShapeAnalysis_Surface.cxx:462-545`), the sequence overload that
/// the free function `projectDegeneratedPoints` (`cxx:600-635`) forwards to.
/// Rewrites the isoline coordinate of the pcurve samples that sit on a
/// singularity so the whole degenerate run collapses onto one parameter value.
pub(super) fn project_degenerated_points(
    s: &dyn Surface,
    sing: &[Singularity],
    preci: f64,
    points: &[GpPnt],
    pnt2d: &mut [GpPnt2d],
    direct: bool,
) -> bool {
    let nbr_pnt = pnt2d.len();
    if nbr_pnt == 0 || points.len() != nbr_pnt {
        return false;
    }
    let step: isize = if direct { 1 } else { -1 };
    let prec2 = preci * preci;
    let j0: isize = if direct { 0 } else { nbr_pnt as isize - 1 };
    // `cxx:480-492`: the singularity nearest to the first (resp. last) sample.
    let mut ind_min: Option<usize> = None;
    let mut gap_min = f64::MAX;
    for (i, sg) in sing.iter().enumerate() {
        if sg.preci > preci {
            break;
        }
        let mut gap2 = sg.p3d.square_distance(&points[j0 as usize]);
        if gap2 > prec2 {
            let q = pnt2d[j0 as usize];
            gap2 = gap2.min(sg.p3d.square_distance(&s.d0(q.x(), q.y())));
        }
        if gap2 <= prec2 && gap_min > gap2 {
            gap_min = gap2;
            ind_min = Some(i);
        }
    }
    let Some(ind) = ind_min else {
        return false;
    };
    let singular = sing[ind];
    // `cxx:503-510`: how far the degenerate run reaches.
    let mut k = j0 + step;
    while k >= 0 && k < nbr_pnt as isize {
        let pk = pnt2d[k as usize];
        let p1 = &points[k as usize];
        let p3 = s.d0(pk.x(), pk.y());
        if singular.p3d.square_distance(p1) > prec2 && singular.p3d.square_distance(&p3) > prec2 {
            break;
        }
        k += step;
    }
    // `cxx:512-530` (PRO7226 #489490): the whole pcurve is degenerate, spread
    // the samples evenly between the two ends of the window.
    if k < 0 || k >= nbr_pnt as isize {
        let (x1, x2) = if singular.u_iso_deg {
            (pnt2d[0].y(), pnt2d[nbr_pnt - 1].y())
        } else {
            (pnt2d[0].x(), pnt2d[nbr_pnt - 1].x())
        };
        for j in 0..nbr_pnt {
            let x = (x1 * (nbr_pnt - 1 - j) as f64 + x2 * j as f64) / (nbr_pnt - 1) as f64;
            if !singular.u_iso_deg {
                pnt2d[j].set_x(x);
            } else {
                pnt2d[j].set_y(x);
            }
        }
        return true;
    }
    // `cxx:532-543`: pin every sample of the run to the isolating coordinate of
    // the last sample outside it.
    let pk = pnt2d[k as usize];
    let mut j = k - step;
    while j >= 0 && j < nbr_pnt as isize {
        if !singular.u_iso_deg {
            pnt2d[j as usize].set_x(pk.x());
        } else {
            pnt2d[j as usize].set_y(pk.y());
        }
        j -= step;
    }
    true
}

/// `approxPCurve` `cxx:1466-1497`: the two `projectDegeneratedPoints` passes
/// (`direct = true` then `false`) followed by the singularity loop that calls
/// `correctExtremity` for the ends of the 3D curve that land on a singularity.
///
/// With `isAnIsoparametric` PARKed (`p01.rs`), the `if (!isoPar2d3d)` guard at
/// `cxx:1467` always holds: `isAnIsoparametric` sets `theIsoPar2d3d = false` on
/// entry (`cxx:2482`) and the port never sets it.
#[allow(clippy::too_many_arguments)]
pub(super) fn correct_degenerated_points(
    s: &dyn Surface,
    curve: &dyn Curve,
    sing: &[Singularity],
    preci: f64,
    the_tol_first: f64,
    the_tol_last: f64,
    points: &[GpPnt],
    params: &[f64],
    pnt2d: &mut [GpPnt2d],
) {
    project_degenerated_points(s, sing, preci, points, pnt2d, true);
    project_degenerated_points(s, sing, preci, points, pnt2d, false);

    // `cxx:1473-1497`. `theTolFirst` / `theTolLast` are `Perform`'s arguments;
    // a negative value means `Precision::Confusion()` (`cxx:1475-1476`).
    let a_tol_first = if the_tol_first < 0.0 { CONFUSION } else { the_tol_first };
    let a_tol_last = if the_tol_last < 0.0 { CONFUSION } else { the_tol_last };
    let a_point_first = points[0];
    let a_point_last = points[points.len() - 1];
    for sg in sing {
        if sg.preci <= CONFUSION && a_point_first.distance(&sg.p3d) <= a_tol_first {
            correct_extremity(curve, params, pnt2d, true, &sg.first_p2d, sg.u_iso_deg, preci, s);
        }
        if sg.preci <= CONFUSION && a_point_last.distance(&sg.p3d) <= a_tol_last {
            // `cxx:1495` passes `aFirstP2d` for the last point too.
            correct_extremity(curve, params, pnt2d, false, &sg.first_p2d, sg.u_iso_deg, preci, s);
        }
    }
}

/// `ShapeAnalysis_Surface::MinLOfBoundedCurve`-free helper: the 2D line/line
/// intersection `IntCurve_IntConicConic` computes for two infinite
/// `gp_Lin2d` (`IntCurve_IntConicConic_1.cxx:730-773`
/// `LineLineGeometricIntersection`, then the point the `Perform` overload
/// (`IntCurve_IntConicConic_1.cxx:1381-1942`) appends for two default
/// `IntRes2d_Domain`s).
///
/// Two parallel lines (`|D| < TOLERANCE_ANGULAIRE`, `1.e-15`,
/// `IntCurve_IntConicConic_1.cxx:42`) produce no isolated point either way:
/// coincident lines (`|D2| <= Tol`) make `LineLineGeometricIntersection`
/// (`IntCurve_IntConicConic_1.cxx:754`) report `nbsol = 2`, and with both
/// `IntRes2d_Domain`s default-constructed (`cxx:1925`) `Perform` appends an
/// infinite *segment* only (`IntCurve_IntConicConic_1.cxx:2056-2058`), leaving
/// `NbPoints() == 0` - `Intersector.Point(1)` would raise `Standard_OutOfRange`
/// there, so the port keeps the iso fallback for that case as well.
/// Otherwise the single point is `ElCLib::Value(U2, L2)`
/// (`IntCurve_IntConicConic_1.cxx:1891`), with
/// `U2 = (Uo21y * U1x - Uo21x * U1y) / D`
/// (`IntCurve_IntConicConic_1.cxx:759`). `CheckLLCoincidence`
/// (`IntCurve_IntConicConic_1.cxx:1363-1377`) returns `false` for domains
/// without first/last points, so `nbsol` stays 1.
fn intersect_lines_2d(
    l1_loc: &GpPnt2d,
    l1_dir: &GpDir2d,
    l2_loc: &GpPnt2d,
    l2_dir: &GpDir2d,
) -> Option<GpPnt2d> {
    let d = l1_dir.y * l2_dir.x - l1_dir.x * l2_dir.y;
    if d.abs() < 1.0e-15 {
        return None;
    }
    let uo21x = l2_loc.x() - l1_loc.x();
    let uo21y = l2_loc.y() - l1_loc.y();
    let u2 = (uo21y * l1_dir.x - uo21x * l1_dir.y) / d;
    Some(GpPnt2d::new(
        l2_loc.x() + u2 * l2_dir.x,
        l2_loc.y() + u2 * l2_dir.y,
    ))
}

/// `adjustSecondToFirstPoint` (`ShapeConstruct_ProjectCurveOnSurface.cxx:253-273`):
/// fold `theSecondPoint` into the half-period window around `theFirstPoint`,
/// per periodic direction.
fn adjust_second_to_first_point(surf: &dyn Surface, first: &GpPnt2d, second: &mut GpPnt2d) {
    if surf.is_u_periodic() {
        let p = surf.u_period();
        let u = occt_core::bspl::locate::in_period(
            second.x(),
            first.x() - 0.5 * p,
            first.x() + 0.5 * p,
        );
        second.set_x(u);
    }
    if surf.is_v_periodic() {
        let p = surf.v_period();
        let v = occt_core::bspl::locate::in_period(
            second.y(),
            first.y() - 0.5 * p,
            first.y() + 0.5 * p,
        );
        second.set_y(v);
    }
}

/// `ShapeConstruct_ProjectCurveOnSurface::correctExtremity`
/// (`ShapeConstruct_ProjectCurveOnSurface.cxx:1910-2037`): walk the pcurve
/// samples from the *second* sample towards the singularity at the end of the
/// 3D curve and replace the end sample by the point where the running chord
/// crosses the endpoint's isoparametric line, i.e. where the pcurve should
/// leave the deformed region.
///
/// `ShapeAnalysis_Surface::NextValueOfUV(FirstPointOfLine, aP3d, myPreci,
/// Precision::Confusion())` (`cxx:2012-2013`) is replaced by
/// [`super::p01::value_of_uv`] (see `p01.rs` on the seeding).
#[allow(clippy::too_many_arguments)]
fn correct_extremity(
    curve: &dyn Curve,
    params: &[f64],
    pnt2d: &mut [GpPnt2d],
    is_first_point: bool,
    point_on_iso_line: &GpPnt2d,
    is_u_iso: bool,
    preci: f64,
    surf: &dyn Surface,
) {
    let nb_pnt = pnt2d.len();
    // OCCT indexes `thePoints2d(3)` / `thePoints2d(2)` unconditionally
    // (`cxx:1935-1936`); `generateCurvePoints` always yields at least
    // `THE_NCONTROL` samples, so this guard is unreachable in practice.
    if nb_pnt < 3 || params.len() != nb_pnt {
        return;
    }
    let ind_coord = if is_u_iso { 2 } else { 1 };
    let singularity_coord = coord_of(point_on_iso_line, 3 - ind_coord);
    let endpoint = if is_first_point {
        pnt2d[0]
    } else {
        pnt2d[nb_pnt - 1]
    };
    let finish_coord = coord_of(&endpoint, 3 - ind_coord);

    // `cxx:1923-1927`: the isoparametric line through the endpoint, parallel to
    // the coordinate that stays free (`gp::DY2d()` for a U-iso singularity,
    // `gp::DX2d()` otherwise).
    let iso_dir = if is_u_iso {
        GpDir2d::new(0.0, 1.0).expect("(0, 1) is a unit 2d direction")
    } else {
        GpDir2d::new(1.0, 0.0).expect("(1, 0) is a unit 2d direction")
    };
    let is_periodic = if is_u_iso {
        surf.is_v_periodic()
    } else {
        surf.is_u_periodic()
    };

    // `cxx:1932-1948`.
    let mut first;
    let mut second;
    let finish_param;
    let mut first_param;
    let mut second_param;
    if is_first_point {
        first = pnt2d[2];
        second = pnt2d[1];
        finish_param = params[0];
        first_param = params[2];
        second_param = params[1];
    } else {
        first = pnt2d[nb_pnt - 3];
        second = pnt2d[nb_pnt - 2];
        finish_param = params[nb_pnt - 1];
        first_param = params[nb_pnt - 3];
        second_param = params[nb_pnt - 2];
    }

    // `cxx:1950-1958`: the singularity and the second sample must lie on the
    // same side of the endpoint, otherwise there is nothing to correct.
    if singularity_coord > finish_coord && coord_of(&second, 3 - ind_coord) > finish_coord {
        return;
    }
    if singularity_coord < finish_coord && coord_of(&second, 3 - ind_coord) < finish_coord {
        return;
    }
    {
        // `cxx:1960-1966`: the second sample must already be far enough along
        // the degenerate coordinate for the end to be worth correcting.
        let prev_dist = (coord_of(&second, ind_coord) - coord_of(&first, ind_coord)).abs();
        let cur_dist = (coord_of(&endpoint, ind_coord) - coord_of(&second, ind_coord)).abs();
        if cur_dist <= 2.0 * prev_dist {
            return;
        }
    }

    // `cxx:1968-1970`.
    let iso_fallback = |second: &GpPnt2d| {
        if is_u_iso {
            GpPnt2d::new(finish_coord, second.y())
        } else {
            GpPnt2d::new(second.x(), finish_coord)
        }
    };
    let mut finish_point = iso_fallback(&second);

    loop {
        // `cxx:1974-1978`
        if (coord_of(&second, 3 - ind_coord) - finish_coord).abs() <= 2.0 * PCONFUSION {
            break;
        }
        let a_vec = GpVec2d::new(second.x() - first.x(), second.y() - first.y());
        if a_vec.square_magnitude() <= 1.0e-32 {
            break;
        }
        // `gp_Dir2d aDir(aVec.X(), aVec.Y())` (`cxx:1986`). OCCT raises
        // `Standard_ConstructionError` when the vector is shorter than
        // `gp::Resolution()`; this port stops the walk instead (no exception
        // machinery, and the `1.e-32` guard above leaves only a degenerate
        // direction).
        let Ok(chord_dir) = GpDir2d::from_vec2d(&a_vec) else {
            break;
        };
        // `cxx:1988-1999`: intersect with the iso line and take the point on
        // the chord; on an empty intersection OCCT keeps the iso fallback.
        finish_point = match intersect_lines_2d(&endpoint, &iso_dir, &first, &chord_dir) {
            Some(p) => p,
            None => iso_fallback(&second),
        };

        // `cxx:2002-2013`.
        let prev_point = first;
        first = second;
        first_param = second_param;
        second_param = 0.5 * (first_param + finish_param);
        if (second_param - first_param).abs() <= 2.0 * PCONFUSION {
            break;
        }
        let a_p3d = curve.d0(second_param);
        second = value_of_uv(surf, &a_p3d, preci);
        if is_periodic {
            adjust_second_to_first_point(surf, &first, &mut second);
        }

        // `cxx:2017-2023`: stop once the samples stop separating.
        let prev_dist = (coord_of(&first, ind_coord) - coord_of(&prev_point, ind_coord)).abs();
        let cur_dist = (coord_of(&second, ind_coord) - coord_of(&first, ind_coord)).abs();
        if cur_dist > 2.0 * prev_dist {
            break;
        }
    }

    // `cxx:2026-2034`.
    if is_first_point {
        pnt2d[0] = finish_point;
    } else {
        pnt2d[nb_pnt - 1] = finish_point;
    }
}

/// `approxPCurve` `cxx:1746-1890` "Handle AdjustOverDegen", enabled by default
/// (`myAdjustOverDegen = 1`, `cxx:644`; the mode is
/// `ShapeConstruct_ProjectCurveOnSurface::AdjustOverDegenMode`, `cxx:692-696`).
///
/// On a closed surface that carries singularities, the point-wise projection of
/// a curve running over the singularity leaves a sample pair half a period
/// apart (`|CurX - PrevX| ~ Up/2`, resp. `Vp/2`). OCCT pins the sample on the
/// window border (`uf` when `myAdjustOverDegen` is set, `ul` otherwise) and
/// rewinds the runs before/after it into `+/-(Up/2 + PConfusion)`, which is the
/// same `Up/2` window as `ElCLib::InPeriod`.
///
/// `mySurf->IsDegenerated(gp_Pnt(0,0,0), myPreci)` (`cxx:1752`, `cxx:1822`) is
/// called for its singularity-computation side effect only; its result is
/// discarded in the C++ and the singularities are already available here.
/// `myStatus |= ShapeExtend_DONE4` (`cxx:1815`, `cxx:1884`) is a status flag
/// only - it changes no geometry and is not carried by this port. `preci` is
/// `myPreci`, the projection tolerance (`cxx:689`), used by
/// `NbSingularities(myPreci)` (`cxx:1753`, `cxx:1823`) and
/// `IsDegenerated(..., Precision::Confusion())` (`cxx:1763`, `cxx:1833`).
///
/// PORTED: the guards below are `IsUClosed(myPreci)` / `IsVClosed(myPreci)`
/// (`cxx:1750`, `cxx:1820`), the full `ShapeAnalysis_Surface` versions
/// ([`super::p01::sa_is_u_closed`] / [`super::p01::sa_is_v_closed`],
/// `ShapeAnalysis_Surface.cxx:660-863`, `:865-1077`), not the bare
/// `Geom_Surface` flags.
pub(super) fn adjust_over_degenerated(
    s: &dyn Surface,
    sing: &[Singularity],
    preci: f64,
    points: &[GpPnt],
    pnt2d: &mut [GpPnt2d],
) {
    let (uf, ul) = s.u_range();
    let (vf, vl) = s.v_range();
    let up = ul - uf;
    let vp = vl - vf;
    if points.len() != pnt2d.len() {
        return;
    }
    let n = pnt2d.len();
    // `cxx:1747`: the block is guarded by `myAdjustOverDegen != -1`; the
    // constructor sets `myAdjustOverDegen(1)` (`cxx:644`), there is no setter
    // and `approxPCurve` never writes it, so the guard always holds here.
    // `cxx:1749-1819` (U) / `cxx:1820-1889` (V). `pin` is the window border the
    // sample is pinned to: `myAdjustOverDegen ? uf : ul` (`cxx:1778`,
    // `cxx:1847`), i.e. `uf` / `vf` for the default mode 1.
    let (is_u, period, pin, win_lo, win_hi) = if sa_is_u_closed(s, preci) {
        (true, up, uf, uf, ul)
    } else if sa_is_v_closed(s, preci) {
        (false, vp, vf, vf, vl)
    } else {
        return;
    };
    if nb_singularities(sing, preci) == 0 {
        return;
    }
    let coord = |p: &GpPnt2d| if is_u { p.x() } else { p.y() };
    let set_coord = |p: &mut GpPnt2d, v: f64| {
        if is_u {
            p.set_x(v);
        } else {
            p.set_y(v);
        }
    };

    // `cxx:1759-1775` / `cxx:1829-1845`.
    let mut prev_coord = 0.0;
    let mut on_bound = false;
    let mut prev_on_bound = false;
    let mut start = true;
    let mut ind = n;
    for i in 0..n {
        let cur_coord = coord(&pnt2d[i]);
        // `cxx:1763-1766`: a sample on the singularity itself carries no
        // parameter information.
        if is_degenerated(sing, &points[i], CONFUSION) {
            continue;
        }
        on_bound =
            ((cur_coord - 0.5 * (win_lo + win_hi)).abs() - 0.5 * period).abs() <= PCONFUSION;
        if !start && ((cur_coord - prev_coord).abs() - 0.5 * period).abs() <= 0.01 * period {
            ind = i;
            break;
        }
        start = false;
        prev_coord = cur_coord;
        prev_on_bound = on_bound;
    }
    if ind >= n {
        return;
    }
    // `cxx:1777-1816` / `cxx:1848-1886`.
    let window = 0.5 * period + PCONFUSION;
    if prev_on_bound {
        set_coord(&mut pnt2d[ind - 1], pin);
        for j in (0..ind.saturating_sub(1)).rev() {
            let mut cur = coord(&pnt2d[j]);
            while cur < pin - window {
                cur += period;
                set_coord(&mut pnt2d[j], cur);
            }
            while cur > pin + window {
                cur -= period;
                set_coord(&mut pnt2d[j], cur);
            }
        }
    } else if on_bound {
        set_coord(&mut pnt2d[ind], pin);
        for j in ind + 1..n {
            let mut cur = coord(&pnt2d[j]);
            while cur < pin - window {
                cur += period;
                set_coord(&mut pnt2d[j], cur);
            }
            while cur > pin + window {
                cur -= period;
                set_coord(&mut pnt2d[j], cur);
            }
        }
    }
}
