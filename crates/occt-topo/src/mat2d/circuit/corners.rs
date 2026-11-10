//! `MAT2d_Circuit`'s corner handling: the static `CrossProd`
//! (`MAT2d_Circuit.cxx:861-873`), `IsSharpCorner` (`cxx:229-341`),
//! `InitOpen` (`cxx:484-503`) and `InsertCorner` (`cxx:674-703`).
//!
//! Both the open-line ends and the sharp junctions are decided by
//! `IsSharpCorner`: with `GeomAbs_Arc` the sign of the tangent cross product
//! answers it, and a re-entrant (rebroussement) junction falls back to
//! intersecting the two offset curves of the incident items, which is why
//! `Adaptor2d_OffsetCurve` and `Geom2dInt_GInter` are reached from here.

use std::sync::Arc;

use occt_core::gp::GpDir2d;
use occt_core::kernel::geomabs::JoinType;
use occt_core::precision::CONFUSION;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::geom2d_int::Geom2dIntGInter;
use occt_geom2d::offset2d::OffsetCurve2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use super::Mat2dCircuit;
use crate::mat2d::item::Mat2dItem;

/// `CrossProd(Geom1, Geom2, DotProd)` (`cxx:861-873`): the cross product of
/// the tangent at the end of `Geom1` with the tangent at the start of `Geom2`,
/// plus their dot product through `dot_prod`.
pub(crate) fn cross_prod(geom1: &Mat2dItem, geom2: &Mat2dItem, dot_prod: &mut f64) -> f64 {
    let c1 = geom1.expect_curve();
    let dir1 = GpDir2d::from_vec2d(&c1.eval_dn(c1.last_parameter(), 1)).expect("gp_Dir2d: null tangent");
    let c2 = geom2.expect_curve();
    let dir2 = GpDir2d::from_vec2d(&c2.eval_dn(c2.first_parameter(), 1)).expect("gp_Dir2d: null tangent");
    *dot_prod = dir1.dot(&dir2);
    dir1.crossed(&dir2)
}

/// `InitOpen(Line)` (`cxx:484-503`): a point at each end of an open line, plus
/// a point at every salient junction.
pub(crate) fn init_open(line: &mut Vec<Mat2dItem>) {
    // `cxx:487-492`.
    let start = line
        .first()
        .expect("MAT2d_Circuit::InitOpen: empty line")
        .start_point();
    line.insert(0, Mat2dItem::Point(start));
    let end = line
        .last()
        .expect("MAT2d_Circuit::InitOpen: empty line")
        .end_point();
    line.push(Mat2dItem::Point(end));

    // `cxx:494-502`.
    let mut i = 2usize;
    while i + 2 <= line.len() {
        let mut dot_prod = 0.0;
        let pro_vec = cross_prod(&line[i - 1], &line[i], &mut dot_prod);
        if pro_vec.abs() > 1.0e-8 || dot_prod < 0.0 {
            let p = line[i - 1].end_point();
            line.insert(i, Mat2dItem::Point(p));
            i += 1;
        }
        i += 1;
    }
}

impl Mat2dCircuit {
    /// `IsSharpCorner(Geom1, Geom2, Direction)` (`cxx:229-341`): true when the
    /// junction between `geom1` and `geom2` is a salient corner for
    /// `direction`.
    pub fn is_sharp_corner(&self, geom1: &Mat2dItem, geom2: &Mat2dItem, direction: f64) -> bool {
        let mut dot_prod = 0.0;
        let mut pro_vec = cross_prod(geom1, geom2, &mut dot_prod);
        let mut nb_test = 1;
        // `cxx:234`: `DU = Precision::Confusion()`.
        let du = CONFUSION;
        let c1 = geom1.expect_curve();
        let c2 = geom2.expect_curve();
        // `cxx:242`.
        let mut tol_ang = 1.0e-8;

        if self.join_type() == JoinType::Arc {
            while nb_test <= 10 {
                if pro_vec * direction < -tol_ang {
                    return true; // Salient.
                }
                if pro_vec * direction > tol_ang {
                    return false; // Re-entrant.
                } else {
                    if dot_prod > 0.0 {
                        return false; // Flat.
                    }
                    tol_ang = 1.0e-8;
                    let u1 = c1.last_parameter() - nb_test as f64 * du;
                    let u2 = c2.first_parameter() + nb_test as f64 * du;
                    let dir1 = GpDir2d::from_vec2d(&c1.eval_dn(u1, 1)).expect("gp_Dir2d: null tangent");
                    let dir2 = GpDir2d::from_vec2d(&c2.eval_dn(u2, 1)).expect("gp_Dir2d: null tangent");
                    dot_prod = dir1.dot(&dir2);
                    pro_vec = dir1.crossed(&dir2);
                    nb_test += 1;
                }
            }

            // `cxx:274-309`: rebroussement - the two curves are offset towards
            // the computation side and intersected; no intersection means a
            // salient corner.
            let tol = CONFUSION;
            let mil_c1 = (c1.last_parameter() + c1.first_parameter()) * 0.5;
            let mil_c2 = (c2.last_parameter() + c2.first_parameter()) * 0.5;
            let p = c1.d0(c1.last_parameter());
            let p1 = c1.d0(mil_c1);
            let p2 = c2.d0(mil_c2);

            let mut d = p1.distance(&p).min(p2.distance(&p));
            d /= 10.0;
            if direction < 0.0 {
                d = -d;
            }

            // `Adaptor2d_OffsetCurve(HC1, D, MilC1, C1->LastParameter())` and
            // `Adaptor2d_OffsetCurve(HC2, D, C2->FirstParameter(), MilC2)`.
            let oc1 = offset_over(c1, d, mil_c1, c1.last_parameter());
            let oc2 = offset_over(c2, d, c2.first_parameter(), mil_c2);

            let mut intersect = Geom2dIntGInter::new();
            intersect.perform_curves(&oc1, &oc2, tol, tol);

            // `cxx:309`.
            return !intersect.is_done() || intersect.nb_points() == 0;
        } else if self.join_type() == JoinType::Intersection {
            // `cxx:311-340`.
            if pro_vec.abs() <= tol_ang && dot_prod < 0.0 {
                while nb_test <= 10 {
                    let u1 = c1.last_parameter() - nb_test as f64 * du;
                    let u2 = c2.first_parameter() + nb_test as f64 * du;
                    let dir1 = GpDir2d::from_vec2d(&c1.eval_dn(u1, 1)).expect("gp_Dir2d: null tangent");
                    let dir2 = GpDir2d::from_vec2d(&c2.eval_dn(u2, 1)).expect("gp_Dir2d: null tangent");
                    // cxx:331 -- DotProd is recomputed here but never read again in
                    // this loop (only ProVec decides), so the store is dead in OCCT too.
                    pro_vec = dir1.crossed(&dir2);
                    if pro_vec * direction < -tol_ang {
                        return true;
                    }
                    if pro_vec * direction > tol_ang {
                        return false;
                    }
                    nb_test += 1;
                }
                return false;
            } else {
                return false;
            }
        }
        false
    }

    /// `InsertCorner(Line)` (`cxx:674-703`): a point is inserted at the start of
    /// every item that begins on a salient corner.
    pub fn insert_corner(&self, line: &mut Vec<Mat2dItem>, direction: f64) {
        let mut i = 1usize;
        while i <= line.len() {
            let isuiv = if i == line.len() { 1 } else { i + 1 };
            let insert = self.is_sharp_corner(&line[i - 1], &line[isuiv - 1], direction);
            if insert {
                let p = line[isuiv - 1].start_point();
                line.insert(i, Mat2dItem::Point(p));
                i += 1;
            }
            i += 1;
        }
    }
}

/// `Adaptor2d_OffsetCurve(Basis, Offset, WFirst, WLast)`: the offset of `c` by
/// `d`, restricted to `[first, last]`.
fn offset_over(c: &Geom2dTrimmedCurve, d: f64, first: f64, last: f64) -> Geom2dTrimmedCurve {
    let basis: Arc<dyn Curve2d> = Arc::new(c.clone());
    let offset: Arc<dyn Curve2d> = Arc::new(OffsetCurve2d::new(basis, d));
    Geom2dTrimmedCurve::new(offset, first, last)
}
