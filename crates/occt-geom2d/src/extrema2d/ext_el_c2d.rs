//! `Extrema_ExtElC2d` -- all extremal distances between two 2D elementary
//! curves. Source: `Extrema_ExtElC2d.cxx` (542 lines) and
//! `Extrema_ExtElC2d.hxx:36-101`.
//!
//! Nine constructors, one per curve pair: `Lin2d x Lin2d` (`:46-101`),
//! `Lin2d x Circ2d` (`:105-167`), `Lin2d x Elips2d` (`:170-222`),
//! `Lin2d x Hypr2d` (`:227-269`), `Lin2d x Parab2d` (`:274-308`),
//! `Circ2d x Circ2d` (`:311-367`), `Circ2d x Elips2d` (`:370-406`),
//! `Circ2d x Hypr2d` (`:410-446`), `Circ2d x Parab2d` (`:450-486`).
//!
//! The three `Circ2d`/conic arms delegate to `Extrema_ExtPElC2d`
//! (`Extrema_ExtPElC2d.cxx`), ported as `super::analytic_solvers::{circle_all2d,
//! ellipse_all2d, hyperbola_all2d, parabola_all2d}`. OCCT drives them with
//! `Precision::Confusion()` and `[0, 2*PI]` for the ellipse, `RealFirst()/
//! RealLast()` for the hyperbola and the parabola; the port passes the same
//! ranges.
//!
//! `Extrema_POnCurv2d` is the existing port
//! `crate::geom2d_int::POnCurv2d` (`geom2d_int/gfunc_ext_pc.rs:22-41`); it is
//! reused rather than duplicated so both callers share one definition.

use super::analytic_solvers::{circle_all2d, ellipse_all2d, hyperbola_all2d, parabola_all2d, PI};
use super::prelude::*;
use crate::geom2d_int::POnCurv2d;
use occt_core::elib::clib2d;

/// `RealLast()` (`Standard_Real.hxx`): the largest finite `double`.
const REAL_LAST: f64 = f64::MAX;
/// `RealFirst()`.
const REAL_FIRST: f64 = f64::MIN;

/// Default `Extrema_POnCurv2d`: parameter `0.0`, point `(0, 0)`.
const ORIGIN: POnCurv2d = POnCurv2d {
    parameter: 0.0,
    point: GpPnt2d::new(0.0, 0.0),
};

/// `gp_Lin2d::SquareDistance(const gp_Pnt2d&)` (`gp_Lin2d.hxx:240-246`):
/// `Cross(theP - Location, Direction)^2`.
fn line_point_square_distance(l: &GpLin2d, p: &GpPnt2d) -> f64 {
    let d = l.pos.vdir;
    let cx = p.x() - l.pos.loc.x();
    let cy = p.y() - l.pos.loc.y();
    let a_d = cx * d.y - cy * d.x;
    a_d * a_d
}

/// `Extrema_ExtElC2d` (`Extrema_ExtElC2d.hxx:36-101`): `mySqDist[8]`,
/// `myPoint[8][2]`, `myNbExt`, `myDone`, `myIsPar`.
#[derive(Debug, Clone)]
pub struct ExtremaExtElC2d {
    my_done: bool,
    my_is_par: bool,
    my_nb_ext: usize,
    my_sq_dist: [f64; 8],
    my_point: [[POnCurv2d; 2]; 8],
}

impl Default for ExtremaExtElC2d {
    /// `Extrema_ExtElC2d()` (`cxx:33-42`).
    fn default() -> Self {
        Self::new()
    }
}

impl ExtremaExtElC2d {
    /// `Extrema_ExtElC2d()` (`cxx:33-42`).
    pub fn new() -> Self {
        Self {
            my_done: false,
            my_is_par: false,
            my_nb_ext: 0,
            my_sq_dist: [REAL_LAST; 8],
            my_point: [[ORIGIN; 2]; 8],
        }
    }

    /// `IsDone()` (`cxx:490-493`).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `IsParallel()` (`cxx:497-505`). OCCT throws `StdFail_NotDone` when the
    /// computation failed; the port reports `myIsPar` (hence `false`).
    pub fn is_parallel(&self) -> bool {
        self.my_done && self.my_is_par
    }

    /// `NbExt()` (`cxx:509-517`). OCCT throws when not done; the port returns 0.
    pub fn nb_ext(&self) -> usize {
        if self.my_done {
            self.my_nb_ext
        } else {
            0
        }
    }

    /// `SquareDistance(N)` (`cxx:520-527`), `N` 1-based.
    pub fn square_distance(&self, n: usize) -> f64 {
        self.my_sq_dist[n - 1]
    }

    /// `Points(N, P1, P2)` (`cxx:532-540`), `N` 1-based. `P1` is on the first
    /// curve, `P2` on the second one.
    pub fn points(&self, n: usize) -> (POnCurv2d, POnCurv2d) {
        (self.my_point[n - 1][0], self.my_point[n - 1][1])
    }

    /// `Extrema_ExtElC2d(const gp_Lin2d&, const gp_Lin2d&, const double)`
    /// (`cxx:46-101`). The `AngTol` argument is ignored by OCCT, which tests
    /// `Precision::Angular()` instead.
    pub fn new_line_line(c1: &GpLin2d, c2: &GpLin2d) -> Self {
        let mut r = Self::new();
        let d1 = GpVec2d::from_dir2d(&c1.pos.vdir);
        let d2 = GpVec2d::from_dir2d(&c2.pos.vdir);
        if d1.is_parallel(&d2, ANGULAR) {
            r.my_is_par = true;
            r.my_sq_dist[0] = line_point_square_distance(c2, &c1.pos.loc);
            r.my_nb_ext = 1;
        } else {
            // Cramer's rule for the intersection; no division by zero since the
            // directions are not parallel.
            let a_p1p2 = GpVec2d::new(c2.pos.loc.x() - c1.pos.loc.x(), c2.pos.loc.y() - c1.pos.loc.y());
            let a_delim = 1.0 / d1.crossed(&d2);
            let a_param1 = a_p1p2.crossed(&d2) * a_delim;
            let a_param2 = -(d1.crossed(&a_p1p2)) * a_delim;
            let p1 = clib2d::line_value_ax2d(a_param1, &c1.pos);
            let p2 = clib2d::line_value_ax2d(a_param2, &c2.pos);
            r.my_sq_dist[0] = 0.0;
            r.my_point[0][0] = POnCurv2d::new(a_param1, p1);
            r.my_point[0][1] = POnCurv2d::new(a_param2, p2);
            r.my_nb_ext = 1;
        }
        r.my_done = true;
        r
    }

    /// `Extrema_ExtElC2d(const gp_Lin2d&, const gp_Circ2d&, const double)`
    /// (`cxx:105-167`). The `Tol` argument is unused by OCCT.
    pub fn new_line_circle(c1: &GpLin2d, c2: &GpCirc2d) -> Self {
        let mut r = Self::new();
        let d = c1.pos.vdir;
        let dx = d.dot(&c2.pos.vxdir);
        let dy = d.dot(&c2.pos.vydir);
        let o1 = c1.pos.loc;
        let mut teta = [0.0f64; 2];
        if dy.abs() <= f64::EPSILON {
            teta[0] = PI / 2.0;
        } else {
            teta[0] = (-dx / dy).atan();
        }
        teta[1] = teta[0] + PI;
        if teta[0] < 0.0 {
            teta[0] += 2.0 * PI;
        }
        for k in 0..2 {
            let p2 = clib2d::circle_value_ax22d(teta[k], &c2.pos, c2.radius);
            let u1 = (p2.x() - o1.x()) * d.x + (p2.y() - o1.y()) * d.y;
            let p1 = clib2d::line_value_ax2d(u1, &c1.pos);
            r.my_sq_dist[k] = p1.square_distance(&p2);
            r.my_point[k][0] = POnCurv2d::new(u1, p1);
            r.my_point[k][1] = POnCurv2d::new(teta[k], p2);
        }
        r.my_nb_ext = 2;
        r.my_done = true;
        r
    }

    /// `Extrema_ExtElC2d(const gp_Lin2d&, const gp_Elips2d&)` (`cxx:170-222`).
    pub fn new_line_ellipse(c1: &GpLin2d, c2: &GpElips2d) -> Self {
        let mut r = Self::new();
        let d = c1.pos.vdir;
        let dx = d.dot(&c2.pos.vxdir);
        let dy = d.dot(&c2.pos.vydir);
        let r1 = c2.major_radius;
        let r2 = c2.minor_radius;
        let o1 = c1.pos.loc;
        let mut teta = [0.0f64; 2];
        if dy.abs() <= f64::EPSILON {
            teta[0] = PI / 2.0;
        } else {
            teta[0] = (-dx * r2 / (dy * r1)).atan();
        }
        teta[1] = teta[0] + PI;
        if teta[0] < 0.0 {
            teta[0] += 2.0 * PI;
        }
        for k in 0..2 {
            let p2 = clib2d::ellipse_value_ax22d(teta[k], &c2.pos, r1, r2);
            let u1 = (p2.x() - o1.x()) * d.x + (p2.y() - o1.y()) * d.y;
            let p1 = clib2d::line_value_ax2d(u1, &c1.pos);
            r.my_sq_dist[k] = p1.square_distance(&p2);
            r.my_point[k][0] = POnCurv2d::new(u1, p1);
            r.my_point[k][1] = POnCurv2d::new(teta[k], p2);
        }
        r.my_nb_ext = 2;
        r.my_done = true;
        r
    }

    /// `Extrema_ExtElC2d(const gp_Lin2d&, const gp_Hypr2d&)` (`cxx:227-269`).
    /// Nothing is computed (and `IsDone()` stays `false`) when the line
    /// direction has no component on the hyperbola Y axis, or when
    /// `R - r*Dx/Dy` vanishes.
    pub fn new_line_hyperbola(c1: &GpLin2d, c2: &GpHypr2d) -> Self {
        let mut r = Self::new();
        let d = c1.pos.vdir;
        let dx = d.dot(&c2.pos.vxdir);
        let dy = d.dot(&c2.pos.vydir);
        let big_r = c2.major_radius;
        let small_r = c2.minor_radius;
        if dy.abs() < f64::EPSILON {
            return r;
        }
        if (big_r - small_r * dx / dy).abs() < f64::EPSILON {
            return r;
        }
        let v2 = (big_r + small_r * dx / dy) / (big_r - small_r * dx / dy);
        let u2 = if v2 > 0.0 { v2.sqrt().ln() } else { 0.0 };
        let p2 = clib2d::hyperbola_value_ax22d(u2, &c2.pos, big_r, small_r);
        let o1 = c1.pos.loc;
        let u1 = (p2.x() - o1.x()) * d.x + (p2.y() - o1.y()) * d.y;
        let p1 = clib2d::line_value_ax2d(u1, &c1.pos);
        r.my_sq_dist[0] = p1.square_distance(&p2);
        r.my_point[0][0] = POnCurv2d::new(u1, p1);
        r.my_point[0][1] = POnCurv2d::new(u2, p2);
        r.my_nb_ext = 1;
        r.my_done = true;
        r
    }

    /// `Extrema_ExtElC2d(const gp_Lin2d&, const gp_Parab2d&)` (`cxx:274-308`).
    /// Nothing is computed when the line direction has no component on the
    /// parabola Y axis.
    pub fn new_line_parabola(c1: &GpLin2d, c2: &GpParab2d) -> Self {
        let mut r = Self::new();
        let d = c1.pos.vdir;
        // `C2.MirrorAxis().Direction()` / `C2.Axis().YAxis().Direction()`.
        let x2 = c2.mirror_axis().vdir;
        let y2 = c2.axis().vydir;
        let dx = d.dot(&x2);
        let dy = d.dot(&y2);
        if dy.abs() < f64::EPSILON {
            return r;
        }
        let p = c2.parameter();
        let u2 = dx * p / dy;
        let p2 = clib2d::parabola_value_ax22d(u2, &c2.pos, c2.focal);
        let o1 = c1.pos.loc;
        let u1 = (p2.x() - o1.x()) * d.x + (p2.y() - o1.y()) * d.y;
        let p1 = clib2d::line_value_ax2d(u1, &c1.pos);
        r.my_sq_dist[0] = p1.square_distance(&p2);
        r.my_point[0][0] = POnCurv2d::new(u1, p1);
        r.my_point[0][1] = POnCurv2d::new(u2, p2);
        r.my_nb_ext = 1;
        r.my_done = true;
        r
    }

    /// `Extrema_ExtElC2d(const gp_Circ2d&, const gp_Circ2d&)` (`cxx:311-367`).
    /// Concentric circles give the constant distance `|R1 - R2|`, otherwise the
    /// four collinear critical pairs. `myDone` is set before the work, so the
    /// four-pair branch also reports `IsDone()`.
    pub fn new_circle_circle(c1: &GpCirc2d, c2: &GpCirc2d) -> Self {
        let mut r = Self::new();
        r.my_done = true;
        let o1 = c1.pos.point;
        let o2 = c2.pos.point;
        let do1o2 = GpVec2d::new(o2.x() - o1.x(), o2.y() - o1.y());
        let a_sq_d_centers = do1o2.square_magnitude();
        if a_sq_d_centers < SQUARE_CONFUSION {
            r.my_is_par = true;
            r.my_nb_ext = 1;
            let a_dr = c1.radius - c2.radius;
            r.my_sq_dist[0] = a_dr * a_dr;
            return r;
        }
        let r1 = c1.radius;
        let r2 = c2.radius;
        let o1o2 = do1o2.multiplied_scalar(1.0 / a_sq_d_centers.sqrt());
        let p1 = [
            o1.translated_vec(&o1o2.multiplied_scalar(r1)),
            o1.translated_vec(&o1o2.multiplied_scalar(-r1)),
        ];
        let usol1 = [clib2d::parameter_circ2d(c1, &p1[0]), clib2d::parameter_circ2d(c1, &p1[1])];
        let p2 = [
            o2.translated_vec(&o1o2.multiplied_scalar(r2)),
            o2.translated_vec(&o1o2.multiplied_scalar(-r2)),
        ];
        let usol2 = [clib2d::parameter_circ2d(c2, &p2[0]), clib2d::parameter_circ2d(c2, &p2[1])];
        for no_sol in 0..2 {
            for kk in 0..2 {
                let n = r.my_nb_ext;
                r.my_sq_dist[n] = p2[kk].square_distance(&p1[no_sol]);
                r.my_point[n][0] = POnCurv2d::new(usol1[no_sol], p1[no_sol]);
                r.my_point[n][1] = POnCurv2d::new(usol2[kk], p2[kk]);
                r.my_nb_ext += 1;
            }
        }
        r
    }

    /// `Extrema_ExtElC2d(const gp_Circ2d&, const gp_Elips2d&)` (`cxx:370-406`):
    /// `Extrema_ExtPElC2d(C1.Location(), C2, Confusion, 0, 2*PI)` then, for each
    /// of its points, `Extrema_ExtPElC2d(that point, C1, Confusion, 0, 2*PI)`.
    pub fn new_circle_ellipse(c1: &GpCirc2d, c2: &GpElips2d) -> Self {
        let mut r = Self::new();
        let o1 = c1.pos.point;
        let (elip_done, elip) = ext_p_el_c2d_ellipse(&o1, c2, CONFUSION, 0.0, 2.0 * PI);
        if elip_done {
            for e_pt in &elip {
                let (circ_done, circ) =
                    ext_p_el_c2d_circle(&e_pt.p2, c1, CONFUSION, 0.0, 2.0 * PI);
                if circ_done {
                    for c_pt in &circ {
                        if r.my_nb_ext == 8 {
                            return r;
                        }
                        let n = r.my_nb_ext;
                        r.my_sq_dist[n] = e_pt.p2.square_distance(&c_pt.p2);
                        r.my_point[n][0] = POnCurv2d::new(c_pt.u2, c_pt.p2);
                        r.my_point[n][1] = POnCurv2d::new(e_pt.u2, e_pt.p2);
                        r.my_nb_ext += 1;
                    }
                }
                r.my_done = true;
            }
        }
        r
    }

    /// `Extrema_ExtElC2d(const gp_Circ2d&, const gp_Hypr2d&)` (`cxx:410-446`).
    /// OCCT drives the point/curve solver with `RealFirst()/RealLast()`.
    pub fn new_circle_hyperbola(c1: &GpCirc2d, c2: &GpHypr2d) -> Self {
        let mut r = Self::new();
        let o1 = c1.pos.point;
        let (hyp_done, hyp) = ext_p_el_c2d_hyperbola(&o1, c2, REAL_FIRST, REAL_LAST);
        if hyp_done {
            for e_pt in &hyp {
                let (circ_done, circ) =
                    ext_p_el_c2d_circle(&e_pt.p2, c1, CONFUSION, 0.0, 2.0 * PI);
                if circ_done {
                    for c_pt in &circ {
                        if r.my_nb_ext == 8 {
                            return r;
                        }
                        let n = r.my_nb_ext;
                        r.my_sq_dist[n] = e_pt.p2.square_distance(&c_pt.p2);
                        r.my_point[n][0] = POnCurv2d::new(c_pt.u2, c_pt.p2);
                        r.my_point[n][1] = POnCurv2d::new(e_pt.u2, e_pt.p2);
                        r.my_nb_ext += 1;
                    }
                }
                r.my_done = true;
            }
        }
        r
    }

    /// `Extrema_ExtElC2d(const gp_Circ2d&, const gp_Parab2d&)` (`cxx:450-486`).
    pub fn new_circle_parabola(c1: &GpCirc2d, c2: &GpParab2d) -> Self {
        let mut r = Self::new();
        let o1 = c1.pos.point;
        let (par_done, par) = ext_p_el_c2d_parabola(&o1, c2, REAL_FIRST, REAL_LAST);
        if par_done {
            for e_pt in &par {
                let (circ_done, circ) =
                    ext_p_el_c2d_circle(&e_pt.p2, c1, CONFUSION, 0.0, 2.0 * PI);
                if circ_done {
                    for c_pt in &circ {
                        if r.my_nb_ext == 8 {
                            return r;
                        }
                        let n = r.my_nb_ext;
                        r.my_sq_dist[n] = e_pt.p2.square_distance(&c_pt.p2);
                        r.my_point[n][0] = POnCurv2d::new(c_pt.u2, c_pt.p2);
                        r.my_point[n][1] = POnCurv2d::new(e_pt.u2, e_pt.p2);
                        r.my_nb_ext += 1;
                    }
                }
                r.my_done = true;
            }
        }
        r
    }
}

// ---------------------------------------------------------------------------
// `Extrema_ExtPElC2d` arms, in the shape `Extrema_ExtElC2d` consumes them:
// `(IsDone(), solutions)` with the curve point in `p2` and its parameter in
// `u2` of each `Extrema2d`.
// ---------------------------------------------------------------------------

/// `Extrema_ExtPElC2d::Perform(P, gp_Elips2d, Tol, Uinf, Usup)`
/// (`Extrema_ExtPElC2d.cxx:166-213`): not done only when `P` is the centre of a
/// circle-like ellipse (`|A - B| <= Tol`).
fn ext_p_el_c2d_ellipse(
    p: &GpPnt2d,
    e: &GpElips2d,
    tol: f64,
    uinf: f64,
    usup: f64,
) -> (bool, Vec<super::analytic_solvers::Extrema2d>) {
    let or = e.pos.point;
    if or.is_equal(p, tol) && (e.major_radius - e.minor_radius).abs() <= tol {
        return (false, Vec::new());
    }
    (true, ellipse_all2d(e, p, uinf, usup))
}

/// `Extrema_ExtPElC2d::Perform(P, gp_Circ2d, Tol, Uinf, Usup)`
/// (`Extrema_ExtPElC2d.cxx:94-154`): not done when `P` is the centre.
fn ext_p_el_c2d_circle(
    p: &GpPnt2d,
    c: &GpCirc2d,
    tol: f64,
    uinf: f64,
    usup: f64,
) -> (bool, Vec<super::analytic_solvers::Extrema2d>) {
    if c.pos.point.is_equal(p, tol) {
        return (false, Vec::new());
    }
    (true, circle_all2d(c, p, uinf, usup))
}

/// `Extrema_ExtPElC2d::Perform(P, gp_Hypr2d, Tol, Uinf, Usup)`
/// (`Extrema_ExtPElC2d.cxx:227-285`): always done once the quartic is solved.
/// `Tol` only feeds the `Tol2` duplicate filter, which is `Precision::Confusion`
/// at every `Extrema_ExtElC2d` call site and is built into `hyperbola_all2d`.
fn ext_p_el_c2d_hyperbola(
    p: &GpPnt2d,
    h: &GpHypr2d,
    uinf: f64,
    usup: f64,
) -> (bool, Vec<super::analytic_solvers::Extrema2d>) {
    (true, hyperbola_all2d(h, p, uinf, usup))
}

/// `Extrema_ExtPElC2d::Perform(P, gp_Parab2d, Tol, Uinf, Usup)`
/// (`Extrema_ExtPElC2d.cxx:297-352`): always done once the cubic is solved.
/// `Tol` only feeds the `Tol2` duplicate filter (see above).
fn ext_p_el_c2d_parabola(
    p: &GpPnt2d,
    pa: &GpParab2d,
    uinf: f64,
    usup: f64,
) -> (bool, Vec<super::analytic_solvers::Extrema2d>) {
    (true, parabola_all2d(pa, p, uinf, usup))
}
