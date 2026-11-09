//! `GccAna_Circ2d3Tan(gp_Pnt2d, gp_Pnt2d, gp_Pnt2d, Tolerance)`
//! (`GccAna_Circ2d3Tan_9.cxx:40-160`): the circle passing through three points.
//!
//! Three cases (`cxx:27-39`):
//! 1. all three points coincide: a zero-radius circle at `Point1`
//!    (`cxx:83-102`),
//! 2. two points coincide: the center is the intersection of the perpendicular
//!    bisector of the two distinct points with the line through them
//!    (`cxx:105-128`),
//! 3. three distinct points: the intersection of the two bisectors
//!    (`cxx:109-115`).
//!
//! Translation notes:
//! - unlike every other constructor this one compares against `Tolerance`
//!   without `std::abs` (`cxx:83`, `cxx:108-126`); the port keeps the raw
//!   argument,
//! - `L1`/`L2` are only conditionally assigned (`cxx:105-128`); when a guard
//!   does not fire OCCT intersects the default `gp_Lin2d`, the port keeps the
//!   default line as well,
//! - the `(1, 1)` OCCT arrays (`cxx:50-64`) are sized at [`super::MAX_SOL`].

use occt_core::elib::clib2d::parameter_circ2d;
use occt_core::gp::circ2d::GpCirc2d;
use occt_core::gp::dir2d::GpDir2d;
use occt_core::gp::lin2d::GpLin2d;
use occt_core::gp::pnt2d::GpPnt2d;
use occt_core::intana2d::IntAna2dAnaIntersection;

use super::super::gcc_ent::GccEntPosition;
use super::{axis_at, GccAnaCirc2d3Tan};

pub(crate) fn new(
    point1: &GpPnt2d,
    point2: &GpPnt2d,
    point3: &GpPnt2d,
    tolerance: f64,
) -> GccAnaCirc2d3Tan {
    let mut me = GccAnaCirc2d3Tan::empty();

    let dist1 = point1.distance(point2);
    let dist2 = point1.distance(point3);
    let dist3 = point2.distance(point3);

    me.qualifier1[0] = GccEntPosition::NoQualifier;
    me.qualifier2[0] = GccEntPosition::NoQualifier;
    me.qualifier3[0] = GccEntPosition::NoQualifier;

    if dist1 < tolerance && dist2 < tolerance && dist3 < tolerance {
        me.nbr_sol += 1;
        me.well_done = true;
        me.cirsol[0] = GpCirc2d::new(axis_at(*point1), 0.0);
        me.the_same1[0] = 0;
        me.the_same2[0] = 0;
        me.the_same3[0] = 0;
        me.pnttg1sol[0] = *point1;
        me.pnttg2sol[0] = *point2;
        me.pnttg3sol[0] = *point3;
        me.par1sol[0] = 0.0;
        me.par2sol[0] = 0.0;
        me.par3sol[0] = 0.0;
        me.pararg1[0] = 0.0;
        me.pararg2[0] = 0.0;
        me.pararg3[0] = 0.0;
    } else {
        let mut l1 = GpLin2d::default();
        let mut l2 = GpLin2d::default();
        if dist1 >= tolerance {
            l1 = GpLin2d::from_pnt_dir(
                GpPnt2d::new(
                    (point1.x() + point2.x()) / 2.0,
                    (point1.y() + point2.y()) / 2.0,
                ),
                GpDir2d::new(point1.y() - point2.y(), point2.x() - point1.x())
                    .expect("gp_Dir2d(Point1.Y() - Point2.Y(), Point2.X() - Point1.X())"),
            );
        }
        if dist2 >= tolerance {
            l2 = GpLin2d::from_pnt_dir(
                GpPnt2d::new(
                    (point1.x() + point3.x()) / 2.0,
                    (point1.y() + point3.y()) / 2.0,
                ),
                GpDir2d::new(point1.y() - point3.y(), point3.x() - point1.x())
                    .expect("gp_Dir2d(Point1.Y() - Point3.Y(), Point3.X() - Point1.X())"),
            );
        }
        if dist2 <= tolerance {
            l2 = GpLin2d::from_pnt_dir(
                *point1,
                GpDir2d::new(point1.y() - point2.y(), point2.x() - point1.x())
                    .expect("gp_Dir2d(Point1.Y() - Point2.Y(), Point2.X() - Point1.X())"),
            );
        } else if dist1 <= tolerance {
            l1 = GpLin2d::from_pnt_dir(
                *point1,
                GpDir2d::new(point1.y() - point3.y(), point3.x() - point1.x())
                    .expect("gp_Dir2d(Point1.Y() - Point3.Y(), Point3.X() - Point1.X())"),
            );
        } else if dist3 <= tolerance {
            l2 = GpLin2d::from_pnt_dir(
                *point1,
                GpDir2d::new(point1.y() - point2.y(), point2.x() - point1.x())
                    .expect("gp_Dir2d(Point1.Y() - Point2.Y(), Point2.X() - Point1.X())"),
            );
        }
        // `cxx:129`: `IntAna2d_AnaIntersection Intp(L1, L2);`
        let mut intp = IntAna2dAnaIntersection::new();
        intp.perform_lin_lin(&l1, &l2);
        if intp.is_done() {
            if !intp.is_empty() {
                for i in 1..=intp.nb_points() {
                    if me.nbr_sol >= super::MAX_SOL {
                        break;
                    }
                    me.nbr_sol += 1;
                    let n = me.nbr_sol;
                    let center = *intp.point(i).value();
                    me.cirsol[n - 1] =
                        GpCirc2d::new(axis_at(center), point1.distance(&center));

                    me.the_same1[n - 1] = 0;
                    me.the_same2[n - 1] = 0;
                    me.the_same3[n - 1] = 0;
                    me.pnttg1sol[n - 1] = *point1;
                    me.pnttg2sol[n - 1] = *point2;
                    me.pnttg3sol[n - 1] = *point3;
                    me.par1sol[n - 1] =
                        parameter_circ2d(&me.cirsol[n - 1], &me.pnttg1sol[n - 1]);
                    me.par2sol[n - 1] =
                        parameter_circ2d(&me.cirsol[n - 1], &me.pnttg2sol[n - 1]);
                    me.par3sol[n - 1] =
                        parameter_circ2d(&me.cirsol[n - 1], &me.pnttg3sol[n - 1]);
                    me.pararg1[n - 1] = 0.0;
                    me.pararg2[n - 1] = 0.0;
                    me.pararg3[n - 1] = 0.0;
                }
            }
            me.well_done = true;
        }
    }
    me
}
