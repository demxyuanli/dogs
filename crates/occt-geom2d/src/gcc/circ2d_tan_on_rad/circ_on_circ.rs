//! `GccAna_Circ2dTanOnRad(QualifiedCirc, gp_Circ2d, Radius, Tolerance)`
//! (`GccAna_Circ2dTanOnRad_3.cxx:40-317`): circles of the given radius tangent
//! to a circle with the centre on a circle.
//!
//! `cxx:93-233` eliminates the configurations the qualifier rules out. Three
//! arms produce a single *boundary* centre through `sign` (`cxx:120-123`,
//! `:159-170`, `:186-200`, `:215-228`); the remaining arms build `nparal`
//! parallel(s) of `C1` at distance `|R1 -/+ Radius|` and intersect them with
//! `OnCirc` (`cxx:237-277`). `cxx:279-304` turns `sign` into the single boundary
//! centre, `cxx:305-315` fills the tangency points and parameters of every
//! stored solution.
//!
//! Translation notes:
//! - `signe[5]` / `disparal[2]` of `cxx:57-65` are stack arrays in OCCT and
//!   become `[i32; 5]` / `[f64; 2]` here, `signe` zero-filled as OCCT does by
//!   hand (`cxx:58-62`), `disparal` written before any read,
//! - `Center(1, 4)` and `dir1on(1, 4)` are uninitialised in OCCT (`cxx:76-77`);
//!   the port seeds `center` with the origin and `dir1on` with `gp::DX2d()` and
//!   reads a slot only after the write OCCT performs,
//! - `cxx:287-301` stores `qualifier1(NbrSol)` while `NbrSol` is still `0`, i.e.
//!   one slot **before** the `qualifier1(1, 4)` array it indexes, and the value
//!   written there is never observed (the branch's only solution stores
//!   `qualifier1` nowhere else). The port drops that out-of-bounds store and
//!   documents it, instead of inventing a slot for it,
//! - `cxx:270` tests `Center(NbrSol).Distance(center1) <= R1`, `cxx:249`/`:284`
//!   recompute `distcc1` for the qualifier; both are kept as written,
//! - the `IsEnclosed` arm of `cxx:106-117` leaves `WellDone` **false** when
//!   `|Radius - R1| < Tol` and `OnCirc` does not pass through `C1`'s centre,
//!   while the `IsEnclosing` arm (`cxx:141-156`) leaves it **true** in the same
//!   situation; the asymmetry is preserved.
//!
//! UNPORTED: none.

use occt_core::elib::clib2d::parameter_circ2d;
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::intana2d::IntAna2dAnaIntersection;

use super::super::circ2d_2tan_on::{axis_at, dir_between, dirx, qualifier_of_circle};
use super::super::gcc_ent::GccEntQualifiedCirc;
use super::{GccAnaCirc2dTanOnRad, MAX_SOL};

pub(crate) fn new(
    qualified1: &GccEntQualifiedCirc,
    on_circ: &GpCirc2d,
    radius: f64,
    tolerance: f64,
) -> GccAnaCirc2dTanOnRad {
    let mut me = GccAnaCirc2dTanOnRad::empty();
    let tol = tolerance.abs();
    let mut signe = [0i32; 5];
    let mut disparal = [0.0f64; 2];
    let mut nparal = 0usize;
    let mut sign = 0i32;

    if !(qualified1.is_enclosed()
        || qualified1.is_enclosing()
        || qualified1.is_outside()
        || qualified1.is_unqualified())
    {
        panic!("GccEnt_BadQualifier in GccAna_Circ2dTanOnRad");
    }

    let c1 = qualified1.qualified();
    let mut center = [GpPnt2d::zero(); MAX_SOL];
    let mut dir1on = [dirx(); MAX_SOL];

    if radius < 0.0 {
        panic!("Standard_NegativeValue in GccAna_Circ2dTanOnRad");
    }

    let r1 = c1.radius();
    let r2 = on_circ.radius();
    let c1loc = c1.location();
    let (c1x, c1y) = (c1loc.x(), c1loc.y());
    let center1 = GpPnt2d::new(c1x, c1y);
    let dist = on_circ.location().distance(&center1);
    let onloc = on_circ.location();
    let (onx, ony) = (onloc.x(), onloc.y());

    if qualified1.is_enclosed() {
        if radius - r1 > tol || r1 - dist - r2 - radius > tol {
            me.well_done = true;
        } else if (radius - r1).abs() < tol {
            if on_circ.distance(&center1) < tol {
                me.nbr_sol = 1;
                me.cirsol[0] = c1;
                me.qualifier1[0] = qualified1.qualifier();
                me.the_same1[0] = 1;
                me.pntcen3[0] = center1;
                me.parcen3[0] = parameter_circ2d(on_circ, &me.pntcen3[0]);
                me.well_done = true;
            } else {
                me.well_done = false;
            }
        } else if (r1 - dist - r2 - radius).abs() <= tol {
            dir1on[0] = dir_between(&GpPnt2d::new(onx, ony), &center1);
            sign = -1;
        } else {
            nparal = 1;
            disparal[0] = (r1 - radius).abs();
        }
    } else if qualified1.is_enclosing() {
        if tol < r1 - radius || tol < r1 + dist - r2 - radius || radius - r1 - dist - r2 > tol {
            me.well_done = true;
        } else if (radius - r1).abs() < tol {
            if on_circ.distance(&center1) < tol {
                me.nbr_sol = 1;
                me.cirsol[0] = c1;
                me.qualifier1[0] = qualified1.qualifier();
                me.the_same1[0] = 1;
                me.pntcen3[0] = center1;
                me.parcen3[0] = parameter_circ2d(on_circ, &me.pntcen3[0]);
                me.well_done = true;
            } else {
                me.well_done = true;
            }
        } else if (r1 + dist - r2 - radius).abs() < tol
            || (radius - r1 - dist - r2).abs() < tol
        {
            dir1on[0] = dir_between(&GpPnt2d::new(onx, ony), &center1);
            if (r1 + dist - r2 - radius).abs() < tol {
                sign = 1;
            } else {
                sign = -1;
            }
        } else {
            nparal = 1;
            disparal[0] = (r1 - radius).abs();
        }
    } else if qualified1.is_outside() {
        if dist - r2 - r1 - radius > tol || radius - dist - r2 + r1 > tol {
            me.well_done = true;
        } else {
            dir1on[0] = dir_between(&GpPnt2d::new(onx, ony), &center1);
            if (r1 - dist - r2 + radius).abs() < tol {
                sign = -1;
            } else if (dist - r1 - r2 - radius).abs() < tol {
                sign = 1;
            } else {
                nparal = 1;
                disparal[0] = (r1 + radius).abs();
            }
        }
    } else if dist - r2 - r1 - radius > tol || radius - dist - r1 - r2 > tol {
        me.well_done = true;
    } else if r1 - dist - r2 > tol && tol < r1 - r2 - dist - radius {
        me.well_done = true;
    } else {
        dir1on[0] = dir_between(&GpPnt2d::new(onx, ony), &center1);
        if (dist - r1 - r2 - radius).abs() < tol {
            sign = 1;
        } else if dist + r1 + r2 - radius < 0.0 {
            sign = -1;
        } else if r1 - dist + r2 > 0.0 && r1 - r2 - dist - radius > 0.0 {
            sign = -1;
        } else {
            nparal = 2;
            disparal[0] = (r1 + radius).abs();
            disparal[1] = (r1 - radius).abs();
        }
    }

    for jj in 0..nparal {
        let mut intp = IntAna2dAnaIntersection::new();
        intp.perform_circ_circ(
            on_circ,
            &GpCirc2d::new(axis_at(center1), disparal[jj]),
        );
        if intp.is_done() {
            if !intp.is_empty() {
                for i in 1..=intp.nb_points() {
                    if me.nbr_sol == MAX_SOL {
                        break;
                    }
                    me.nbr_sol += 1;
                    let n = me.nbr_sol;
                    center[n - 1] = *intp.point(i).value();
                    me.cirsol[n - 1] = GpCirc2d::new(axis_at(center[n - 1]), radius);
                    let distcc1 = center[n - 1].distance(&center1);
                    me.qualifier1[n - 1] = qualifier_of_circle(qualified1, distcc1, radius, r1, tol);
                    dir1on[n - 1] = dir_between(&center[n - 1], &center1);
                    if r1 > radius && center[n - 1].distance(&center1) <= r1 {
                        signe[n - 1] = -1;
                    } else {
                        signe[n - 1] = 1;
                    }
                }
            }
            me.well_done = true;
        }
    }

    if sign != 0 {
        center[0] = GpPnt2d::new(
            onloc.x() + sign as f64 * r2 * dir1on[0].x,
            onloc.y() + sign as f64 * r2 * dir1on[0].y,
        );
        me.cirsol[0] = GpCirc2d::new(axis_at(center[0]), radius);
        // `cxx:287-301`: the qualifier is stored at `qualifier1(NbrSol)` with
        // `NbrSol == 0`, i.e. one slot before this array; OCCT never reads that
        // slot back, so the port drops the store.
        me.nbr_sol = 1;
        signe[0] = 1;
        me.well_done = true;
    }

    let mut ii = 0usize;
    while signe[ii] != 0 {
        let k = ii + 1;
        me.pnttg1sol[k - 1] = GpPnt2d::new(
            center[k - 1].x() + signe[ii] as f64 * radius * dir1on[k - 1].x,
            center[k - 1].y() + signe[ii] as f64 * radius * dir1on[k - 1].y,
        );
        me.pntcen3[k - 1] = me.cirsol[k - 1].location();
        me.pararg1[k - 1] = parameter_circ2d(&c1, &me.pnttg1sol[k - 1]);
        me.par1sol[k - 1] = parameter_circ2d(&me.cirsol[k - 1], &me.pnttg1sol[k - 1]);
        me.parcen3[k - 1] = parameter_circ2d(on_circ, &me.pntcen3[k - 1]);
        ii += 1;
    }
    me
}
