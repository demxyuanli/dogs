//! `IntTools_EdgeEdge` general-curve branch — the parameter-box recursion.
//!
//! Faithful port of the OCCT members and file-local helpers that
//! `IntTools_EdgeEdge::Perform` uses for every pair that is not line/line
//! (`IntTools_EdgeEdge.cxx`):
//!
//! | OCCT | lines | here |
//! |---|---|---|
//! | `Prepare` (res coeff / resolution / PTol / type swap) | `:89-181` | [`EdgeEdge::prepare`] (`edge_edge.rs`) |
//! | `FindSolutions(ranges1, ranges2, bSplit2)` | `:290-349` | [`EdgeEdge::find_solutions`] |
//! | `FindSolutions(R1, B1, R2, B2, …)` (recursion) | `:353-549` | [`EdgeEdge::find_solutions_range`] |
//! | `FindParameters` | `:553-671` | [`find_parameters`] |
//! | `MergeSolutions` | `:675-779` | [`EdgeEdge::merge_solutions_boxes`] |
//! | `AddSolution` | `:780-825` | [`EdgeEdge::add_solution_occt`] |
//! | `FindBestSolution` | `:826-901` | [`EdgeEdge::find_best_solution_ranges`] |
//! | `IsIntersection` | `:1060-1146` | [`EdgeEdge::is_intersection`] |
//! | `CheckCoincidence` | `:1150-1206` | [`EdgeEdge::check_coincidence`] |
//! | `FindDistPC` / `DistPC` | `:1210-1362` | [`find_dist_pc`] / [`dist_pc`] |
//! | `SplitRangeOnSegments` | `:1366-1406` | [`split_range_on_segments`] |
//! | `BndBuildBox` | `:1410-1419` | [`bnd_build_box`] |
//! | `PointBoxDistance` | `:1423-1452` | [`point_box_distance`] |
//! | `TypeToInteger` | `:1456-1482` | [`type_to_integer`] |
//! | `ResolutionCoeff` | `:1486-1559` | [`resolution_coeff`] |
//! | `Resolution` | `:1561-1607` | [`resolution_of`] |
//! | `CurveDeflection` | `:1611-1638` | [`curve_deflection`] |
//! | `IsClosed` | `:1642-1659` | [`is_closed`] |
//!
//! Result mapping (the port's [`CommonPrt`] keeps less than
//! `IntTools_CommonPrt`): a `TopAbs_EDGE` solution becomes a
//! [`CommonPartType::Edge`] entry whose `range` is the OCCT `Range1` (already
//! expressed on the *caller's* first edge, see [`EdgeEdge::add_solution_occt`]),
//! and a `TopAbs_VERTEX` solution becomes one [`PntOn2Faces`] whose `uv1`/`uv2`
//! are the `FindBestSolution` parameters on the caller's first/second edge.
//! `IntTools_CommonPrt::{Range2, Edge1, Edge2, BoundingPoints}` are not stored
//! by that reduced struct (pre-existing data-model reduction; no consumer reads
//! them today).
//!
//! The B-rep adaptor range (`myCurve1`/`myCurve2`) is the *edge* range, so the
//! port calls `D0`/`D1` on the edge curve directly with edge parameters —
//! `IntTools_EdgeEdge.cxx:94-95,164-165` initializes `BRepAdaptor_Curve` on the
//! edge, whose `FirstParameter`/`LastParameter` are the edge's.

use super::prelude::*;
use super::*;

use occt_core::bnd::BndBox;
use occt_core::precision::{epsilon, INFINITE, RESOLUTION, SQUARE_CONFUSION};
use occt_geom::extrema_pc::extrema_ext_pc_range;

use crate::int_tools_curve_box::add_curve_to_box;

// ---------------------------------------------------------------------------
// Curve type (`GeomAdaptor_Curve::GetType`)
// ---------------------------------------------------------------------------

/// `GeomAbs_CurveType` — the arms `IntTools_EdgeEdge` distinguishes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum GeomCurveKind {
    Line,
    Circle,
    Ellipse,
    Hyperbola,
    Parabola,
    BezierCurve,
    BSplineCurve,
    OffsetCurve,
    OtherCurve,
}

/// `GeomAdaptor_Curve::GetType()` for a port curve.
///
/// The order is `GeomAdaptor_Curve::load` (`GeomAdaptor_Curve.cxx:252-311`)
/// verbatim: trimmed curves are unwrapped to their basis (the port's trimmed
/// types forward these queries), then Circle, Line, Ellipse, Parabola,
/// Hyperbola, Bezier, BSpline, Offset, Other. Only the concrete class decides —
/// a geometrically straight B-spline is `GeomAbs_BSplineCurve` in OCCT, so
/// `Curve::is_line` (a geometric test) is deliberately **not** used here.
pub(crate) fn curve_kind(c: &dyn Curve) -> GeomCurveKind {
    if c.gp_circ().is_some() {
        GeomCurveKind::Circle
    } else if c.gp_line().is_some() {
        GeomCurveKind::Line
    } else if c.gp_ellipse().is_some() {
        GeomCurveKind::Ellipse
    } else if c.gp_parabola().is_some() {
        GeomCurveKind::Parabola
    } else if c.gp_hyperbola().is_some() {
        GeomCurveKind::Hyperbola
    } else if c.bezier_poles().is_some() {
        GeomCurveKind::BezierCurve
    } else if c.bspline_knots().is_some() || c.nurbs_degree().is_some() {
        GeomCurveKind::BSplineCurve
    } else if c.offset_curve().is_some() {
        GeomCurveKind::OffsetCurve
    } else {
        GeomCurveKind::OtherCurve
    }
}

/// `TypeToInteger` (`IntTools_EdgeEdge.cxx:1456-1482`).
pub(crate) fn type_to_integer(kind: GeomCurveKind) -> i32 {
    match kind {
        GeomCurveKind::Line => 0,
        GeomCurveKind::Hyperbola | GeomCurveKind::Parabola => 1,
        GeomCurveKind::Circle | GeomCurveKind::Ellipse => 2,
        GeomCurveKind::BezierCurve | GeomCurveKind::BSplineCurve => 3,
        _ => 4,
    }
}

// ---------------------------------------------------------------------------
// Resolution helpers
// ---------------------------------------------------------------------------

/// `ResolutionCoeff` (`IntTools_EdgeEdge.cxx:1486-1559`).
pub(super) fn resolution_coeff(c: &dyn Curve, kind: GeomCurveKind, range: IntRange) -> f64 {
    match kind {
        GeomCurveKind::Circle => match c.circle_radius() {
            Some(r) => 1.0 / (2.0 * r),
            None => 0.0,
        },
        GeomCurveKind::Ellipse => match c.gp_ellipse() {
            Some(e) => 1.0 / e.major_radius(),
            None => 0.0,
        },
        GeomCurveKind::OffsetCurve => {
            // `cxx:1501-1520`: offset-of-line leaves the switch with 0 (its
            // `Resolution` arm returns `theR3D`), offset-of-circle/ellipse
            // return an analytic coefficient, every other basis falls through
            // into the sampling arm.
            let Some((basis, offset)) = c.offset_curve() else {
                return 0.0;
            };
            match curve_kind(&*basis) {
                GeomCurveKind::Line => 0.0,
                GeomCurveKind::Circle => match basis.circle_radius() {
                    Some(r) => 1.0 / (2.0 * (offset + r)),
                    None => 0.0,
                },
                GeomCurveKind::Ellipse => match basis.gp_ellipse() {
                    Some(e) => 1.0 / (offset + e.major_radius()),
                    None => 0.0,
                },
                _ => sample_resolution_coeff(c, range),
            }
        }
        GeomCurveKind::Hyperbola | GeomCurveKind::Parabola | GeomCurveKind::OtherCurve => {
            sample_resolution_coeff(c, range)
        }
        _ => 0.0,
    }
}

/// The shared sampling arm of `ResolutionCoeff` (`cxx:1524-1551`): 30 equal
/// steps over the range, the smallest `dt / chord` ratio (or the initial `10.`
/// when every step is degenerate).
fn sample_resolution_coeff(c: &dyn Curve, range: IntRange) -> f64 {
    let a_nb_p = 30;
    let (a_t1, a_t2) = (range.first, range.last);
    let a_dt = (a_t2 - a_t1) / a_nb_p as f64;
    let mut a_t = a_t1;
    let mut k_min = 10.0;
    let mut a_p1 = c.d0(a_t1);
    for _ in 1..=a_nb_p {
        a_t += a_dt;
        let a_p2 = c.d0(a_t);
        let k = a_dt / a_p1.distance(&a_p2);
        if k < k_min {
            k_min = k;
        }
        a_p1 = a_p2;
    }
    k_min
}

/// `Resolution` (`IntTools_EdgeEdge.cxx:1561-1607`).
pub(super) fn resolution_of(c: &dyn Curve, kind: GeomCurveKind, res_coeff: f64, r3d: f64) -> f64 {
    match kind {
        GeomCurveKind::Line => r3d,
        GeomCurveKind::Circle => {
            let a_dt = res_coeff * r3d;
            if a_dt <= 1.0 {
                2.0 * a_dt.asin()
            } else {
                2.0 * std::f64::consts::PI
            }
        }
        // `Geom_BezierCurve::Resolution` / `Geom_BSplineCurve::Resolution`
        // (both `BSplCLib::Resolution`), reached through the `Curve` trait.
        GeomCurveKind::BezierCurve | GeomCurveKind::BSplineCurve => c.resolution(r3d),
        GeomCurveKind::OffsetCurve => {
            let Some((basis, _)) = c.offset_curve() else {
                return res_coeff * r3d;
            };
            match curve_kind(&*basis) {
                GeomCurveKind::Line => r3d,
                GeomCurveKind::Circle => {
                    let a_dt = res_coeff * r3d;
                    if a_dt <= 1.0 {
                        2.0 * a_dt.asin()
                    } else {
                        2.0 * std::f64::consts::PI
                    }
                }
                _ => res_coeff * r3d,
            }
        }
        _ => res_coeff * r3d,
    }
}

/// `CurveDeflection` (`IntTools_EdgeEdge.cxx:1611-1638`): the total turn of the
/// tangent over 10 equal steps, skipping steps with a degenerate derivative.
pub(super) fn curve_deflection(c: &dyn Curve, range: IntRange) -> f64 {
    let a_nb_p = 10;
    let (a_t1, a_t2) = (range.first, range.last);
    let a_dt = (a_t2 - a_t1) / a_nb_p as f64;
    let mut a_t = a_t1;
    let mut a_defl = 0.0;
    let mut a_v1 = c.d1(a_t1).1;
    for _ in 1..=a_nb_p {
        a_t += a_dt;
        let a_v2 = c.d1(a_t).1;
        if a_v1.magnitude() > RESOLUTION && a_v2.magnitude() > RESOLUTION {
            a_defl += a_v1.angle(&a_v2);
        }
        a_v1 = a_v2;
    }
    a_defl
}

/// `IsClosed` (`IntTools_EdgeEdge.cxx:1642-1659`).
pub(super) fn is_closed(c: &dyn Curve, a_t1: f64, a_t2: f64, tol: f64, res: f64) -> bool {
    if (a_t1 - a_t2).abs() < res {
        return false;
    }
    let a_p1 = c.d0(a_t1);
    let a_p2 = c.d0(a_t2);
    a_p1.distance(&a_p2) < tol
}

// ---------------------------------------------------------------------------
// Range / box helpers
// ---------------------------------------------------------------------------

/// `SplitRangeOnSegments` (`IntTools_EdgeEdge.cxx:1366-1406`).
pub(super) fn split_range_on_segments(
    a_t1: f64,
    a_t2: f64,
    the_resolution: f64,
    the_nb_seg: i32,
) -> Vec<IntRange> {
    let a_diff = a_t2 - a_t1;
    if a_diff < the_resolution || the_nb_seg == 1 {
        return vec![IntRange::new_unchecked(a_t1, a_t2)];
    }
    let mut a_nb_segments = the_nb_seg;
    let mut a_dt = a_diff / a_nb_segments as f64;
    if a_dt < the_resolution {
        let a_seg = a_diff / the_resolution;
        a_nb_segments = a_seg as i32 + 1;
        a_dt = a_diff / a_nb_segments as f64;
    }
    let mut out = Vec::new();
    let mut a_t1x = a_t1;
    for _ in 1..a_nb_segments {
        let a_t2x = a_t1x + a_dt;
        out.push(IntRange::new_unchecked(a_t1x, a_t2x));
        a_t1x = a_t2x;
    }
    out.push(IntRange::new_unchecked(a_t1x, a_t2));
    out
}

/// `BndBuildBox` (`IntTools_EdgeEdge.cxx:1410-1419`) — the OCCT call is
/// `BndLib_Add3dCurve::Add(theBAC, aT1, aT2, theTol, aB)`.
pub(super) fn bnd_build_box(c: &dyn Curve, a_t1: f64, a_t2: f64, the_tol: f64) -> BndBox {
    let mut a_b = BndBox::new();
    add_curve_to_box(c, a_t1, a_t2, the_tol, &mut a_b);
    a_b
}

/// `PointBoxDistance` (`IntTools_EdgeEdge.cxx:1423-1452`).
///
/// `Bnd_Box::Get` already adds the box gap, exactly as `BndBox::distance`
/// does, so this is the same formula. A void box is `+inf` here (OCCT's `Get`
/// on a void box yields `RealLast` coordinates, i.e. a huge finite distance).
pub(super) fn point_box_distance(b: &BndBox, p: &GpPnt) -> f64 {
    b.distance(p)
}

// ---------------------------------------------------------------------------
// Projection onto a restricted range (`GeomAPI_ProjectPointOnCurve`)
// ---------------------------------------------------------------------------

/// `GeomAPI_ProjectPointOnCurve` restricted to `[t1, t2]`
/// (`GeomAPI_ProjectPointOnCurve.cxx:135-156`): `NbPoints() == 0` is the
/// "cannot project" answer, otherwise the *smallest* squared distance over all
/// extrema (no `myIsMin` filter) and its parameter.
struct ProjPc<'a> {
    curve: &'a dyn Curve,
    t1: f64,
    t2: f64,
}

impl<'a> ProjPc<'a> {
    fn new(curve: &'a dyn Curve, t1: f64, t2: f64) -> Self {
        Self { curve, t1, t2 }
    }

    /// `Perform` + `NbPoints()`: `None` when nothing projects.
    fn project(&self, p: &GpPnt) -> Option<(f64, f64)> {
        let sols = extrema_ext_pc_range(self.curve, p, self.t1, self.t2);
        let best = sols.iter().reduce(|a, b| if b.sq_dist < a.sq_dist { b } else { a })?;
        Some((best.sq_dist.sqrt(), best.u))
    }
}

/// Running extremum search state shared by `DistPC` / `FindDistPC`
/// (`IntTools_EdgeEdge.cxx:1210-1362`); `d_max` is a *maximum* for `i_c = 1`
/// and a *minimum* for `i_c = -1`, exactly as in OCCT.
#[derive(Debug, Clone, Copy)]
struct DistState {
    d_max: f64,
    t1_max: f64,
    t2_max: f64,
}

/// The 7-argument `DistPC` (`IntTools_EdgeEdge.cxx:1332-1362`): the error code
/// plus the distance and the second parameter.
fn dist_pc(
    proj: &ProjPc<'_>,
    c1: &dyn Curve,
    a_t1: f64,
    the_criteria: f64,
    i_c: i32,
) -> (i32, f64, f64) {
    let a_p1 = c1.d0(a_t1);
    let Some((a_d, a_t2)) = proj.project(&a_p1) else {
        // `iErr = 1`: the point of C1 cannot be projected on C2.
        return (1, 0.0, 0.0);
    };
    let i_err = if (i_c as f64) * (a_d - the_criteria) > 0.0 {
        2 // the distance is too big or too small
    } else {
        0
    };
    (i_err, a_d, a_t2)
}

/// The 9-argument `DistPC` (`IntTools_EdgeEdge.cxx:1301-1328`): as above plus
/// the running extremum update. Returns `(iErr, aD)`.
fn dist_pc_max(
    proj: &ProjPc<'_>,
    c1: &dyn Curve,
    a_t1: f64,
    the_criteria: f64,
    i_c: i32,
    st: &mut DistState,
) -> (i32, f64) {
    let (i_err, a_d, a_t2) = dist_pc(proj, c1, a_t1, the_criteria, i_c);
    if i_err == 1 {
        return (i_err, a_d);
    }
    if (i_c as f64) * (a_d - st.d_max) > 0.0 {
        st.d_max = a_d;
        st.t1_max = a_t1;
        st.t2_max = a_t2;
    }
    (i_err, a_d)
}

/// `FindDistPC` (`IntTools_EdgeEdge.cxx:1210-1297`) — golden-section search of
/// the extremal distance of C1 on `[a_t1a, a_t1b]` against the projector.
fn find_dist_pc(
    proj: &ProjPc<'_>,
    c1: &dyn Curve,
    a_t1a: f64,
    a_t1b: f64,
    the_criteria: f64,
    the_eps: f64,
    st: &mut DistState,
    b_max_dist: bool,
) -> i32 {
    let i_c = if b_max_dist { 1 } else { -1 };
    let a_gs = 0.618_033_988_749_894_8; // = 0.5*(1+sqrt(5))-1
    let mut a_a = a_t1a;
    let mut a_b = a_t1b;
    st.t1_max = 0.0;
    st.t2_max = 0.0;

    // check bounds
    let (i_err, _) = dist_pc_max(proj, c1, a_a, the_criteria, i_c, st);
    if i_err == 2 {
        return i_err;
    }
    let (i_err, _) = dist_pc_max(proj, c1, a_b, the_criteria, i_c, st);
    if i_err == 2 {
        return i_err;
    }
    let mut a_xp = a_a + (a_b - a_a) * a_gs;
    let mut a_xl = a_b - (a_b - a_a) * a_gs;
    let (i_err, mut a_yp) = dist_pc_max(proj, c1, a_xp, the_criteria, i_c, st);
    if i_err != 0 {
        return i_err;
    }
    let (i_err, mut a_yl) = dist_pc_max(proj, c1, a_xl, the_criteria, i_c, st);
    if i_err != 0 {
        return i_err;
    }

    let an_eps = the_eps.max(epsilon(a_a.abs().max(a_b.abs())) * 10.0);
    loop {
        let i_err;
        if (i_c as f64) * (a_yp - a_yl) > 0.0 {
            a_a = a_xl;
            a_xl = a_xp;
            a_yl = a_yp;
            a_xp = a_a + (a_b - a_a) * a_gs;
            let (e, d) = dist_pc_max(proj, c1, a_xp, the_criteria, i_c, st);
            i_err = e;
            a_yp = d;
        } else {
            a_b = a_xp;
            a_xp = a_xl;
            a_yp = a_yl;
            a_xl = a_b - (a_b - a_a) * a_gs;
            let (e, d) = dist_pc_max(proj, c1, a_xl, the_criteria, i_c, st);
            i_err = e;
            a_yl = d;
        }
        if i_err != 0 {
            if i_err == 2 && !b_max_dist {
                let a_xp_mid = (a_a + a_b) * 0.5;
                let _ = dist_pc_max(proj, c1, a_xp_mid, the_criteria, i_c, st);
            }
            return i_err;
        }
        if (a_b - a_a) < an_eps {
            break;
        }
    }
    0
}

// ---------------------------------------------------------------------------
// `FindParameters`
// ---------------------------------------------------------------------------

/// `IntTools_EdgeEdge::FindParameters` (`IntTools_EdgeEdge.cxx:553-671`).
///
/// Walks from each end of `[a_t1, a_t2]` towards the box, growing the step by
/// the local curve resolution, then bisects (golden section) for the boundary
/// parameter. `a_tb1`/`a_tb2` are the caller's *current* range values, which
/// OCCT passes by reference and only overwrites when a bounding point is
/// found. `None` reproduces OCCT's `false` return (the edge is out of the box).
#[allow(clippy::too_many_arguments)]
pub(super) fn find_parameters(
    c: &dyn Curve,
    kind: GeomCurveKind,
    a_t1: f64,
    a_t2: f64,
    the_tol: f64,
    the_res: f64,
    the_ptol: f64,
    the_res_coeff: f64,
    the_cbox: &BndBox,
    a_tb1_in: f64,
    a_tb2_in: f64,
) -> Option<(f64, f64)> {
    let a_cf = 0.618_033_988_749_894_8; // = 0.5*(1+sqrt(5))/2
    let mut a_cbx = *the_cbox;
    a_cbx.set_gap(a_cbx.gap() + the_tol);
    let a_max_dt = (a_t2 - a_t1) * 0.01;

    let mut a_tb1 = a_tb1_in;
    let mut a_tb2 = a_tb2_in;
    let mut b_ret = false;

    for i in 0..2 {
        let mut a_tb = if i == 0 { a_t1 } else { a_t2 };
        let mut a_t = if i == 0 { a_t2 } else { a_tb1 };
        let a_c = if i == 0 { 1.0 } else { -1.0 };
        let mut a_dt = the_res;
        let mut a_distp = 0.0;
        b_ret = false;
        let mut k = 1.0;
        // Looking for the point on the edge which is in the box.
        while a_c * (a_t - a_tb) >= 0.0 {
            let a_p = c.d0(a_tb);
            let a_dist = point_box_distance(the_cbox, &a_p);
            if a_dist > the_tol {
                if a_distp > 0.0 {
                    let mut to_grow = false;
                    if (a_distp - a_dist).abs() / a_distp < 0.1 {
                        a_dt = resolution_of(c, kind, the_res_coeff, k * a_dist);
                        if a_dt < a_max_dt {
                            to_grow = true;
                            k *= 2.0;
                        }
                    }
                    if !to_grow {
                        k = 1.0;
                        a_dt = resolution_of(c, kind, the_res_coeff, a_dist);
                    }
                }
                a_tb += a_c * a_dt;
            } else {
                b_ret = true;
                break;
            }
            a_distp = a_dist;
        }
        if !b_ret {
            if i == 0 {
                // Edge is out of the box.
                return None;
            }
            b_ret = true;
            a_tb = a_tb1;
            a_dt = a_t2 - a_tb1;
        }
        a_t = if i == 0 { a_t1 } else { a_t2 };
        if a_tb != a_t {
            // One point IN, one point OUT: look for the bounding point.
            let mut a_tin = a_tb;
            let mut a_tout = a_tb - a_c * a_dt;
            let mut a_diff = a_tin - a_tout;
            while a_diff.abs() > the_ptol {
                a_tb = a_tout + a_diff * a_cf;
                let a_p = c.d0(a_tb);
                if a_cbx.is_out(&a_p) {
                    a_tout = a_tb;
                } else {
                    a_tin = a_tb;
                }
                a_diff = a_tin - a_tout;
            }
            if i == 0 {
                a_tb1 = a_tb;
            } else {
                a_tb2 = a_tb;
            }
        }
    }
    if b_ret {
        Some((a_tb1, a_tb2))
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// Member helpers
// ---------------------------------------------------------------------------

impl EdgeEdge {
    /// `IntTools_EdgeEdge::FindSolutions(theRanges1, theRanges2, bSplit2)`
    /// (`IntTools_EdgeEdge.cxx:290-349`) followed by `MergeSolutions`
    /// (`Perform`, `:241-242`).
    ///
    /// Replaces the previous substitute ("the extrema of the distance within
    /// tolerance are the intersections"): OCCT recurses on the parameter boxes
    /// so that every crossing whose boxes stay in contact is reported, even
    /// when it is not a distance extremum of the whole range.
    pub(super) fn find_solutions(&mut self) {
        let (a_t11, a_t12) = (self.r1().first, self.r1().last);
        let (a_t21, a_t22) = (self.r2().first, self.r2().last);
        let c1 = self.curve1.clone().unwrap();
        let c2 = self.curve2.clone().unwrap();

        let mut b_is_closed2 = is_closed(&*c2, a_t21, a_t22, self.tol2, self.res2);
        if b_is_closed2 {
            let a_b1 = bnd_build_box(&*c1, a_t11, a_t12, self.tol1);
            let a_p = c2.d0(a_t21);
            b_is_closed2 = !a_b1.is_out(&a_p);
        }

        let mut ranges1: Vec<IntRange> = Vec::new();
        let mut ranges2: Vec<IntRange> = Vec::new();
        let b_split2;

        if !b_is_closed2 {
            let a_b1 = bnd_build_box(&*c1, a_t11, a_t12, self.tol1);
            let a_b2 = bnd_build_box(&*c2, a_t21, a_t22, self.tol2);
            self.find_solutions_range(
                IntRange::new_unchecked(a_t11, a_t12),
                &a_b1,
                IntRange::new_unchecked(a_t21, a_t22),
                &a_b2,
                &mut ranges1,
                &mut ranges2,
            );
            self.merge_solutions_boxes(&ranges1, &ranges2, false);
            return;
        }

        if self.check_coincidence(a_t11, a_t12, a_t21, a_t22, self.tol, self.res1) == 0 {
            // `cxx:320-325`: the whole range is one common part when
            // `CheckCoincidence` reports the patches as coincident.
            ranges1.push(self.r1());
            ranges2.push(self.r2());
            self.merge_solutions_boxes(&ranges1, &ranges2, false);
            return;
        }

        let a_nb1 = if is_closed(&*c1, a_t11, a_t12, self.tol1, self.res1) {
            2
        } else {
            1
        };
        let a_segments1 = split_range_on_segments(a_t11, a_t12, self.res1, a_nb1);
        let a_segments2 = split_range_on_segments(a_t21, a_t22, self.res2, 2);

        for a_r1 in &a_segments1 {
            let a_b1 = bnd_build_box(&*c1, a_r1.first, a_r1.last, self.tol1);
            for a_r2 in &a_segments2 {
                let a_b2 = bnd_build_box(&*c2, a_r2.first, a_r2.last, self.tol2);
                self.find_solutions_range(*a_r1, &a_b1, *a_r2, &a_b2, &mut ranges1, &mut ranges2);
            }
        }
        b_split2 = a_segments2.len() > 1;
        self.merge_solutions_boxes(&ranges1, &ranges2, b_split2);
    }

    /// `IntTools_EdgeEdge::FindSolutions(theR1, theBox1, theR2, theBox2, …)`
    /// (`IntTools_EdgeEdge.cxx:353-549`) — the recursion.
    #[allow(clippy::too_many_arguments)]
    fn find_solutions_range(
        &self,
        the_r1: IntRange,
        the_box1: &BndBox,
        the_r2: IntRange,
        the_box2: &BndBox,
        the_ranges1: &mut Vec<IntRange>,
        the_ranges2: &mut Vec<IntRange>,
    ) {
        let c1 = self.curve1.as_ref().unwrap();
        let c2 = self.curve2.as_ref().unwrap();
        let (mut a_t11, mut a_t12) = (the_r1.first, the_r1.last);
        let (mut a_t21, mut a_t22) = (the_r2.first, the_r2.last);
        let mut a_b1 = *the_box1;
        // `cxx:371`: `aB2` starts at the incoming box and is rebuilt at step 2;
        // the rebuilt value carries into the next iteration's step 1.
        let mut a_b2 = *the_box2;

        let mut b_out;
        let mut b_thin = false;
        let mut b_stop = false;
        let mut i_com = 1;
        loop {
            let (a_tb11, a_tb12) = (a_t11, a_t12);
            let (a_tb21, a_tb22) = (a_t21, a_t22);

            // 1. Find the parameters of the second edge inside the box of the
            //    first one.
            b_out = a_b1.is_out_box(&a_b2);
            if b_out {
                break;
            }
            b_thin = (a_t12 - a_t11) < self.res1
                || (a_b1.is_x_thin(self.tol)
                    && a_b1.is_y_thin(self.tol)
                    && a_b1.is_z_thin(self.tol));
            match find_parameters(
                &**c2,
                self.kind2,
                a_tb21,
                a_tb22,
                self.tol2,
                self.res2,
                self.ptol2,
                self.res_coeff2,
                &a_b1,
                a_t21,
                a_t22,
            ) {
                None => {
                    b_out = true;
                }
                Some((t21, t22)) => {
                    a_t21 = t21;
                    a_t22 = t22;
                    b_out = false;
                }
            }
            if b_out || b_thin {
                break;
            }

            // 2. Build the box of the second edge and find the parameters of
            //    the first one inside it.
            a_b2 = bnd_build_box(&**c2, a_t21, a_t22, self.tol2);
            b_out = a_b1.is_out_box(&a_b2);
            if b_out {
                break;
            }
            b_thin = (a_t22 - a_t21) < self.res2
                || (a_b2.is_x_thin(self.tol)
                    && a_b2.is_y_thin(self.tol)
                    && a_b2.is_z_thin(self.tol));
            match find_parameters(
                &**c1,
                self.kind1,
                a_tb11,
                a_tb12,
                self.tol1,
                self.res1,
                self.ptol1,
                self.res_coeff1,
                &a_b2,
                a_t11,
                a_t12,
            ) {
                None => {
                    b_out = true;
                }
                Some((t11, t12)) => {
                    a_t11 = t11;
                    a_t12 = t12;
                    b_out = false;
                }
            }
            if b_out || b_thin {
                break;
            }

            // 3. Check whether it makes sense to continue.
            let mut a_small_step1 = (a_tb12 - a_tb11) / 250.0;
            let mut a_small_step2 = (a_tb22 - a_tb21) / 250.0;
            if a_small_step1 < self.res1 {
                a_small_step1 = self.res1;
            }
            if a_small_step2 < self.res2 {
                a_small_step2 = self.res2;
            }
            if ((a_t11 - a_tb11) < a_small_step1)
                && ((a_tb12 - a_t12) < a_small_step1)
                && ((a_t21 - a_tb21) < a_small_step2)
                && ((a_tb22 - a_t22) < a_small_step2)
            {
                b_stop = true;
            } else {
                a_b1 = bnd_build_box(&**c1, a_t11, a_t12, self.tol1);
            }
            if b_stop {
                break;
            }
        }

        if b_out {
            // no intersection
            return;
        }
        if !b_thin {
            // check the curves for coincidence on the ranges
            i_com = self.check_coincidence(a_t11, a_t12, a_t21, a_t22, self.tol, self.res1);
            if i_com == 0 {
                b_thin = true;
            }
        }
        if b_thin {
            if i_com != 0 {
                // check the intermediate points
                let a_t1 = (a_t11 + a_t12) * 0.5;
                let a_p1 = c1.d0(a_t1);
                let proj = ProjPc::new(&**c2, a_t21, a_t22);
                let b_sol = match proj.project(&a_p1) {
                    Some((d, _)) => d <= self.tol,
                    None => {
                        let a_t2 = (a_t21 + a_t22) * 0.5;
                        let a_p2 = c2.d0(a_t2);
                        a_p1.distance(&a_p2) <= self.tol
                    }
                };
                if !b_sol {
                    return;
                }
            }
            // add the common part
            the_ranges1.push(IntRange::new_unchecked(a_t11, a_t12));
            the_ranges2.push(IntRange::new_unchecked(a_t21, a_t22));
            return;
        }

        if !self.is_intersection(a_t11, a_t12, a_t21, a_t22) {
            return;
        }

        // Split the ranges on segments and repeat.
        let a_b1_first = bnd_build_box(&**c1, a_t11, a_t12, self.tol1);
        let a_b1_sq_extent = a_b1_first.square_extent();
        let a_r2 = IntRange::new_unchecked(a_t21, a_t22);
        let a_b2_seg = bnd_build_box(&**c2, a_t21, a_t22, self.tol2);
        let a_segments1 = split_range_on_segments(a_t11, a_t12, self.res1, 3);
        let a_nb1 = a_segments1.len();
        for a_r1 in a_segments1 {
            let a_b1s = bnd_build_box(&**c1, a_r1.first, a_r1.last, self.tol1);
            if !a_b1s.is_out_box(&a_b2_seg)
                && (a_nb1 == 1 || a_b1s.square_extent() < a_b1_sq_extent)
            {
                self.find_solutions_range(a_r1, &a_b1s, a_r2, &a_b2_seg, the_ranges1, the_ranges2);
            }
        }
    }

    /// `IntTools_EdgeEdge::MergeSolutions` (`IntTools_EdgeEdge.cxx:675-779`):
    /// groups the parameter ranges found by the recursion and stores them.
    fn merge_solutions_boxes(
        &mut self,
        the_ranges1: &[IntRange],
        the_ranges2: &[IntRange],
        b_split2: bool,
    ) {
        let a_nb_cp = the_ranges1.len();
        if a_nb_cp == 0 {
            return;
        }
        let c1 = self.curve1.clone().unwrap();
        let c2 = self.curve2.clone().unwrap();
        let a_res1 = resolution_of(&*c1, self.kind1, self.res_coeff1, self.tol);
        let a_res2 = resolution_of(&*c2, self.kind2, self.res_coeff2, self.tol);
        let (a_t11, a_t12) = (self.r1().first, self.r1().last);
        let (a_t21, a_t22) = (self.r2().first, self.r2().last);
        let d_tr1 = 20.0 * a_res1;
        let d_tr2 = 20.0 * a_res2;
        let mut a_type_edge = false;
        let mut used = vec![false; a_nb_cp];

        let mut i = 0usize;
        while i < a_nb_cp {
            if used[i] {
                i += 1;
                continue;
            }
            let (mut a_ti11, mut a_ti12) = (the_ranges1[i].first, the_ranges1[i].last);
            let (mut a_ti21, mut a_ti22) = (the_ranges2[i].first, the_ranges2[i].last);
            used[i] = true;

            let mut j = i + 1;
            while j < a_nb_cp {
                if used[j] {
                    j += 1;
                    continue;
                }
                let (a_tj11, a_tj12) = (the_ranges1[j].first, the_ranges1[j].last);
                let (a_tj21, a_tj22) = (the_ranges2[j].first, the_ranges2[j].last);
                let mut b_cond = ((a_ti12 - a_tj11).abs() < d_tr1)
                    || (a_tj11 > a_ti11 && a_tj11 < a_ti12)
                    || (a_ti11 > a_tj11 && a_ti11 < a_tj12)
                    || (b_split2 && (a_tj12 - a_ti11).abs() < d_tr1);
                if b_cond && b_split2 {
                    b_cond = ((a_ti22.max(a_tj22) - a_ti21.min(a_tj21))
                        - ((a_ti22 - a_ti21) + (a_tj22 - a_tj21)))
                        .abs()
                        < d_tr2
                        || (a_tj21 > a_ti21 && a_tj21 < a_ti22)
                        || (a_ti21 > a_tj21 && a_ti21 < a_tj22);
                }
                if b_cond {
                    a_ti11 = a_ti11.min(a_tj11);
                    a_ti12 = a_ti12.max(a_tj12);
                    a_ti21 = a_ti21.min(a_tj21);
                    a_ti22 = a_ti22.max(a_tj22);
                    used[j] = true;
                    j += 1;
                } else if !b_split2 {
                    i = j;
                    break;
                } else {
                    j += 1;
                }
            }

            if ((a_t11 - a_ti11).abs() < self.res1 && (a_t12 - a_ti12).abs() < self.res1)
                || ((a_t21 - a_ti21).abs() < self.res2 && (a_t22 - a_ti22).abs() < self.res2)
            {
                a_type_edge = true;
                self.common_parts.clear();
            }

            self.add_solution_occt(a_ti11, a_ti12, a_ti21, a_ti22, a_type_edge);
            if a_type_edge {
                break;
            }
            if b_split2 {
                i += 1;
            }
        }
    }

    /// `IntTools_EdgeEdge::AddSolution` (`IntTools_EdgeEdge.cxx:780-825`).
    ///
    /// `mySwap` is honoured exactly as OCCT does: the stored range and the
    /// vertex parameters are expressed on the **caller's** first/second edge.
    fn add_solution_occt(&mut self, a_t11: f64, a_t12: f64, a_t21: f64, a_t22: f64, is_edge: bool) {
        if is_edge {
            let range = if !self.swapped {
                IntRange::new_unchecked(a_t11, a_t12)
            } else {
                IntRange::new_unchecked(a_t21, a_t22)
            };
            self.common_parts
                .push(CommonPrt::with(CommonPartType::Edge, range, None, Vec::new()));
            return;
        }
        let (a_t1, a_t2) = self.find_best_solution_ranges(a_t11, a_t12, a_t21, a_t22);
        // `cxx:810-819`: when swapped, `aT1` lives on the caller's *second*
        // edge and the vertex parameters are exchanged.
        let (e1, e2) = if !self.swapped { (a_t1, a_t2) } else { (a_t2, a_t1) };
        let Some(c1) = self.curve1.as_ref() else {
            return;
        };
        let p = c1.d0(a_t1);
        self.points.push(PntOn2Faces::new(0, 1, p, p, (e1, 0.0), (e2, 0.0)));
    }

    /// `IntTools_EdgeEdge::FindBestSolution` (`IntTools_EdgeEdge.cxx:826-901`):
    /// the parameter pair of the closest approach inside the solution range.
    fn find_best_solution_ranges(
        &self,
        a_t11: f64,
        a_t12: f64,
        a_t21: f64,
        a_t22: f64,
    ) -> (f64, f64) {
        const A_SOL_CRITERIA: f64 = 5.0e-16;
        const A_TOUCH_CRITERIA: f64 = 5.0e-13;
        let c1 = self.curve1.as_ref().unwrap();
        let c2 = self.curve2.as_ref().unwrap();
        let mut a_d_min = INFINITE;
        let mut a_t1 = 0.0;
        let mut a_t2 = 0.0;
        let mut b_touch = false;
        let mut b_touch_confirm = false;
        let (mut a_t11_touch, mut a_t12_touch) = (a_t11, a_t12);
        let (mut a_t21_touch, mut a_t22_touch) = (a_t21, a_t22);
        let mut is_sol_found = false;

        let a_res1 = resolution_of(&**c1, self.kind1, self.res_coeff1, self.tol);
        let a_nb_s = 10;
        let a_ranges = split_range_on_segments(a_t11, a_t12, 3.0 * a_res1, a_nb_s);
        let proj = ProjPc::new(&**c2, a_t21, a_t22);

        for a_r1 in a_ranges {
            let mut st = DistState { d_max: self.tol, t1_max: 0.0, t2_max: 0.0 };
            let i_err = find_dist_pc(
                &proj,
                &**c1,
                a_r1.first,
                a_r1.last,
                A_SOL_CRITERIA,
                self.ptol1,
                &mut st,
                false,
            );
            let a_d = st.d_max;
            if i_err != 1 {
                if a_d < a_d_min {
                    a_t1 = st.t1_max;
                    a_t2 = st.t2_max;
                    a_d_min = a_d;
                    is_sol_found = true;
                }
                if a_d < A_TOUCH_CRITERIA {
                    if b_touch {
                        a_t12_touch = st.t1_max;
                        a_t22_touch = st.t2_max;
                        b_touch_confirm = true;
                    } else {
                        a_t11_touch = st.t1_max;
                        a_t21_touch = st.t2_max;
                        b_touch = true;
                    }
                }
            }
        }

        if !is_sol_found || b_touch_confirm {
            a_t1 = (a_t11_touch + a_t12_touch) * 0.5;
            let (i_err, _d, t2) = dist_pc(&proj, &**c1, a_t1, A_SOL_CRITERIA, -1);
            a_t2 = t2;
            if i_err == 1 {
                a_t2 = (a_t21_touch + a_t22_touch) * 0.5;
            }
        }
        (a_t1, a_t2)
    }

    /// `IntTools_EdgeEdge::CheckCoincidence` (`IntTools_EdgeEdge.cxx:1150-1206`).
    ///
    /// `0` = the patches are coincident, `1` = a point of C1 cannot be
    /// projected, `2` = the distance is too big.
    fn check_coincidence(
        &self,
        a_t11: f64,
        a_t12: f64,
        a_t21: f64,
        a_t22: f64,
        the_criteria: f64,
        the_curve_res1: f64,
    ) -> i32 {
        let c1 = self.curve1.as_ref().unwrap();
        let c2 = self.curve2.as_ref().unwrap();
        let proj = ProjPc::new(&**c2, a_t21, a_t22);
        let a_nb = 10;
        let a_ranges = split_range_on_segments(a_t11, a_t12, the_curve_res1, a_nb);
        let a_nb1 = a_ranges.len() as i32;

        // 1. Express evaluation. The 7-argument `DistPC` writes the distance
        //    straight into `aDmax` (`cxx:1174`), so there is no extremum
        //    tracking here.
        let mut a_d_max = -1.0;
        let mut i_err = 0;
        for k in 0..a_nb1.saturating_sub(1) {
            let a_t1b = a_ranges[k as usize].last;
            let (e, d, _t2) = dist_pc(&proj, &**c1, a_t1b, the_criteria, 1);
            a_d_max = d;
            i_err = e;
            if i_err != 0 {
                return i_err;
            }
        }
        // If the ranges are shorter than `theCurveRes1`, there is no need for
        // the deep evaluation.
        if a_nb1 < a_nb {
            return i_err;
        }

        // 2. Deep evaluation.
        for k in 1..a_nb1.saturating_sub(1) {
            let (a_t1a, a_t1b) = (a_ranges[k as usize].first, a_ranges[k as usize].last);
            let mut st = DistState { d_max: a_d_max, t1_max: 0.0, t2_max: 0.0 };
            i_err = find_dist_pc(
                &proj,
                &**c1,
                a_t1a,
                a_t1b,
                the_criteria,
                the_curve_res1,
                &mut st,
                true,
            );
            a_d_max = st.d_max;
            if i_err != 0 {
                return i_err;
            }
        }
        i_err
    }

    /// `IntTools_EdgeEdge::IsIntersection` (`IntTools_EdgeEdge.cxx:1060-1146`).
    fn is_intersection(&self, a_t11: f64, a_t12: f64, a_t21: f64, a_t22: f64) -> bool {
        let c1 = self.curve1.as_ref().unwrap();
        let c2 = self.curve2.as_ref().unwrap();
        let a_coef;
        if ((a_t12 - a_t11) > 1.0e5 * self.res1) && ((a_t22 - a_t21) > 1.0e5 * self.res2) {
            a_coef = 5000.0;
        } else {
            let a_tr_min = ((a_t12 - a_t11) / self.res1).min((a_t22 - a_t21) / self.res2);
            let a_coef_tmp = a_tr_min / 100.0;
            a_coef = if a_coef_tmp < 1.0 { 1.0 } else { a_coef_tmp };
        }
        let mut a_criteria = a_coef * self.tol;
        a_criteria *= a_criteria;

        let (a_p11, a_v11) = c1.d1(a_t11);
        let (a_p12, a_v12) = c1.d1(a_t12);
        let (a_p21, a_v21) = c2.d1(a_t21);
        let (a_p22, a_v22) = c2.d1(a_t22);

        let b_small_11_21 = a_p11.square_distance(&a_p21) < a_criteria;
        let b_small_11_22 = a_p11.square_distance(&a_p22) < a_criteria;
        let b_small_12_21 = a_p12.square_distance(&a_p21) < a_criteria;
        let b_small_12_22 = a_p12.square_distance(&a_p22) < a_criteria;

        let mut b_ret = true;
        if (b_small_11_21 && b_small_12_22) || (b_small_11_22 && b_small_12_21) {
            if a_coef == 1.0 {
                return b_ret;
            }
            let an_angle_criteria = 5.0e-3;
            let mut an_angle1 = 0.0;
            let mut an_angle2 = 0.0;
            let v_ok = |v: &GpVec| v.dot(v) > SQUARE_CONFUSION;
            if v_ok(&a_v11) && v_ok(&a_v12) && v_ok(&a_v21) && v_ok(&a_v22) {
                if b_small_11_21 && b_small_12_22 {
                    an_angle1 = a_v11.angle(&a_v21);
                    an_angle2 = a_v12.angle(&a_v22);
                } else {
                    an_angle1 = a_v11.angle(&a_v22);
                    an_angle2 = a_v12.angle(&a_v21);
                }
            }
            if (an_angle1 < an_angle_criteria
                || (std::f64::consts::PI - an_angle1) < an_angle_criteria)
                || (an_angle2 < an_angle_criteria
                    || (std::f64::consts::PI - an_angle2) < an_angle_criteria)
            {
                let proj = ProjPc::new(&**c2, a_t21, a_t22);
                let mut st = DistState { d_max: INFINITE, t1_max: 0.0, t2_max: 0.0 };
                let i_err =
                    find_dist_pc(&proj, &**c1, a_t11, a_t12, self.tol, self.res1, &mut st, false);
                b_ret = i_err == 2;
            }
        }
        b_ret
    }
}
