//! Port of `IntCurve_Polygon2dGen` instantiated as
//! `Geom2dInt_ThePolygon2dOfTheIntPCurvePCurveOfGInter`
//! (`IntCurve_Polygon2dGen.gxx:38-383` + `.lxx:20-93`,
//! `Geom2dInt_ThePolygon2dOfTheIntPCurvePCurveOfGInter_0.cxx:27-...`).
//!
//! This is the `ThePolygon2d` / `IntCurve_ThePolygon2d` of the
//! `IntCurve_IntPolyPolyGen` engine. It discretizes a `Curve2d` over an
//! `IntRes2dDomain` into a polyline, estimates the deflection, refines the
//! point set against the other polygon's bounding box (`ComputeWithBox`) and
//! exposes the `Intf_Polygon2d` interface the `occt-core` interference engine
//! consumes.
//!
//! The curve is deliberately NOT stored: every OCCT member that needs it
//! (`cxx:122` `ComputeWithBox`) takes it as a parameter, exactly like the
//! original.

use occt_core::bnd::box2d::BndBox2d;
use occt_core::gp::{GpDir2d, GpLin2d, GpPnt2d, GpVec2d};
use occt_core::intf::IntfPolygon2d;
use occt_core::intres2d::IntRes2dDomain;

use crate::curve::Curve2d;
use super::curve_tool;

/// `MAJORATION_DEFLECTION` (`gxx:27`).
const MAJORATION_DEFLECTION: f64 = 1.5;

/// `Geom2dInt_ThePolygon2dOfTheIntPCurvePCurveOfGInter`
/// (`Geom2dInt_ThePolygon2dOfTheIntPCurvePCurveOfGInter.hxx:36-...`).
///
/// The 1-based `NCollection_Array1` members of the generator are kept as
/// `Vec`s accessed through `t*` helpers so that the transcription below stays
/// line-by-line comparable with `IntCurve_Polygon2dGen.gxx`.
#[derive(Clone, Debug)]
pub struct Geom2dIntPolygon2d {
    /// `ThePnts` (`gxx:45`).
    the_pnts: Vec<GpPnt2d>,
    /// `TheParams` (`gxx:46`).
    the_params: Vec<f64>,
    /// `TheIndex` (`gxx:47`).
    the_index: Vec<usize>,
    /// `NbPntIn` (`hxx`).
    nb_pnt_in: usize,
    /// `TheMaxNbPoints` (`hxx`).
    the_max_nb_points: usize,
    /// `myBox` (inherited from `Intf_Polygon2d`).
    my_box: BndBox2d,
    /// `TheDeflection` (`hxx`).
    the_deflection: f64,
    /// `ClosedPolygon` (`hxx`).
    closed_polygon: bool,
    /// `Binf` (`hxx`).
    binf: f64,
    /// `Bsup` (`hxx`).
    bsup: f64,
}

impl Geom2dIntPolygon2d {
    /// `IntCurve_Polygon2dGen(C, tNbPts, D, Tol)` (`gxx:38-117`).
    pub fn new(c: &dyn Curve2d, t_nb_pts: usize, d: &IntRes2dDomain, tol: f64) -> Self {
        let nb_pts = if t_nb_pts < 3 { 3 } else { t_nb_pts };
        let alloc = if t_nb_pts < 3 { 6 } else { t_nb_pts + t_nb_pts };
        let mut p = Self {
            the_pnts: vec![GpPnt2d::zero(); alloc],
            the_params: vec![0.0; alloc],
            the_index: vec![0; alloc],
            nb_pnt_in: nb_pts,
            the_max_nb_points: nb_pts + nb_pts,
            my_box: BndBox2d::new(),
            the_deflection: 0.0,
            closed_polygon: false,
            binf: 0.0,
            bsup: 0.0,
        };
        //-----------------------------------------------------
        //--- Initialization of the breaking with d_Parametre constant
        //---
        p.binf = d.first_parameter();
        p.bsup = d.last_parameter();
        //-----------------------------------------------------
        let mut u = p.binf;
        let u1 = p.bsup;
        let du = (u1 - u) / (nb_pts - 1) as f64;
        let mut i = 1usize;
        loop {
            let pnt = curve_tool::value(c, u);
            p.my_box.add_point(&pnt);
            p.set_ti(i, i);
            p.set_tp(i, pnt);
            p.set_tpar(i, u);
            u += du;
            i += 1;
            if i > nb_pts {
                break;
            }
        }

        //-----------------------------------------------------
        //--- Calculate a maximal deflection
        //---
        p.the_deflection = 0.000000001_f64.min(tol / 100.);
        i = 1;
        u = d.first_parameter();
        u += du * 0.5;
        loop {
            let pm = curve_tool::value(c, u);
            let p1 = p.tp(i);
            let p2 = p.tp(i + 1);
            u += du;
            i += 1;

            let mut t = 0.0;
            let mut dx = p1.x() - p2.x();
            if dx < 0.0 {
                dx = -dx;
            }
            let mut dy = p1.y() - p2.y();
            if dy < 0.0 {
                dy = -dy;
            }
            if dx + dy > 1e-12 {
                if let Ok(dir) = GpDir2d::from_vec2d(&GpVec2d::new(p2.x() - p1.x(), p2.y() - p1.y())) {
                    let l = GpLin2d::from_pnt_dir(p1, dir);
                    t = l.distance(&pm);
                }
                if t > p.the_deflection {
                    p.the_deflection = t;
                }
            }
            if i >= nb_pts {
                break;
            }
        }

        p.my_box.enlarge(p.the_deflection * MAJORATION_DEFLECTION);
        p.closed_polygon = false;
        p
    }

    // -----------------------------------------------------------------------
    // 1-based `NCollection_Array1` accessors (`Value` / `SetValue`).
    // -----------------------------------------------------------------------

    #[inline]
    fn tp(&self, i: usize) -> GpPnt2d {
        self.the_pnts[i - 1]
    }
    #[inline]
    fn set_tp(&mut self, i: usize, v: GpPnt2d) {
        self.the_pnts[i - 1] = v;
    }
    #[inline]
    fn tpar(&self, i: usize) -> f64 {
        self.the_params[i - 1]
    }
    #[inline]
    fn set_tpar(&mut self, i: usize, v: f64) {
        self.the_params[i - 1] = v;
    }
    #[inline]
    fn ti(&self, i: usize) -> usize {
        self.the_index[i - 1]
    }
    #[inline]
    fn set_ti(&mut self, i: usize, v: usize) {
        self.the_index[i - 1] = v;
    }

    /// `NbPntIn` (`hxx`).
    #[inline]
    pub fn nb_pnt_in(&self) -> usize {
        self.nb_pnt_in
    }

    /// `SetDeflectionOverEstimation(x)` (`lxx:25-29`).
    pub fn set_deflection_over_estimation(&mut self, x: f64) {
        self.the_deflection = x;
        self.my_box.enlarge(self.the_deflection);
    }

    /// `Closed(flag)` (`lxx:32-35`).
    pub fn closed(&mut self, flag: bool) {
        self.closed_polygon = flag;
    }

    /// `InfParameter()` (`lxx:43-46`).
    pub fn inf_parameter(&self) -> f64 {
        self.tpar(self.ti(1))
    }

    /// `SupParameter()` (`lxx:50-53`).
    pub fn sup_parameter(&self) -> f64 {
        self.tpar(self.ti(self.nb_pnt_in))
    }

    /// `CalculRegion(x, y, x1, x2, y1, y2)` (`lxx:57-93`).
    fn calcul_region(x: f64, y: f64, x1: f64, x2: f64, y1: f64, y2: f64) -> u8 {
        let mut r: u8;
        if x < x1 {
            r = 1;
        } else if x > x2 {
            r = 2;
        } else {
            r = 0;
        }
        if y < y1 {
            r |= 4;
        } else if y > y2 {
            r |= 8;
        }
        r
    }

    /// `ComputeWithBox(C, BoxOtherPolygon)` (`gxx:122-266`).
    pub fn compute_with_box(&mut self, c: &dyn Curve2d, box_other_polygon: &BndBox2d) {
        if self.my_box.is_out_box(box_other_polygon) {
            self.nb_pnt_in = 2;
            self.my_box.set_void();
        } else {
            let (mut bx0, mut by0, mut bx1, mut by1) =
                box_other_polygon.get().expect("ComputeWithBox: void box");
            bx0 -= self.the_deflection;
            by0 -= self.the_deflection;
            bx1 += self.the_deflection;
            by1 += self.the_deflection;

            let mut max_index_used = 1usize;
            let mut nbp = 0usize;

            let x = self.tp(self.ti(1)).x();
            let y = self.tp(self.ti(1)).y();
            let mut rprec = Self::calcul_region(x, y, bx0, bx1, by0, by1);
            for i in 2..=self.nb_pnt_in {
                let p2d = self.tp(self.ti(i));
                let ri = Self::calcul_region(p2d.x(), p2d.y(), bx0, bx1, by0, by1);
                if (ri & rprec) == 0 {
                    if nbp != 0 {
                        if self.ti(nbp) != self.ti(i - 1) {
                            nbp += 1;
                            let v = self.ti(i - 1);
                            self.set_ti(nbp, v);
                        }
                    } else {
                        nbp += 1;
                        let v = self.ti(i - 1);
                        self.set_ti(nbp, v);
                    }
                    nbp += 1;
                    let v = self.ti(i);
                    self.set_ti(nbp, v);
                    if self.ti(i) > max_index_used {
                        max_index_used = self.ti(i);
                    }
                    rprec = ri;
                }
                rprec = ri;
            }

            if nbp == 1 {
                self.nb_pnt_in = 2;
                self.my_box.set_void();
            } else {
                self.my_box.set_void();
                if nbp != 0 {
                    let v = self.tp(self.ti(1));
                    self.my_box.add_point(&v);
                }
                let mut nb_passage_deflection = 0;
                loop {
                    nb_passage_deflection += 1;
                    let mut new_deflection = self.the_deflection;
                    // `for (i = 2; i <= nbp; i++)` (`gxx:198`) with the `i--`
                    // at `gxx:239`; expressed as a `while` so the rewind is
                    // observable.
                    let mut i = 2usize;
                    while i <= nbp {
                        let ii = self.ti(i);
                        let iim1 = self.ti(i - 1);
                        let pi = self.tp(ii);
                        let pim1 = self.tp(iim1);
                        self.my_box.add_point(&pi);
                        let regi = Self::calcul_region(pi.x(), pi.y(), bx0, bx1, by0, by1);
                        let regim1 = Self::calcul_region(pim1.x(), pim1.y(), bx0, bx1, by0, by1);
                        if (regi & regim1) == 0 {
                            let u = 0.5 * (self.tpar(ii) + self.tpar(iim1));
                            let pm = curve_tool::value(c, u);
                            let mut t = 0.0;
                            let mut dx = pim1.x() - pi.x();
                            if dx < 0.0 {
                                dx = -dx;
                            }
                            let mut dy = pim1.y() - pi.y();
                            if dy < 0.0 {
                                dy = -dy;
                            }
                            if dx + dy > 1e-12 {
                                let p1 = pim1;
                                // `gp_Lin2d L(Pim1, gp_Dir2d(gp_Vec2d(Pim1, Pi)))`
                                // (`gxx:220`). OCCT's `gp_Dir2d` raises
                                // `Standard_ConstructionError` on a direction
                                // shorter than `gp::Resolution()`; the port
                                // leaves `t` at 0 in that degenerate case.
                                if let Ok(dir) = GpDir2d::from_vec2d(&GpVec2d::new(
                                    pi.x() - pim1.x(),
                                    pi.y() - pim1.y(),
                                )) {
                                    let l = GpLin2d::from_pnt_dir(pim1, dir);
                                    t = l.distance(&pm);
                                }
                                if (max_index_used < (self.the_max_nb_points - 1))
                                    && (t > (self.the_deflection * 0.5))
                                {
                                    nbp += 1;
                                    let mut j = nbp;
                                    while j >= i + 1 {
                                        let v = self.ti(j - 1);
                                        self.set_ti(j, v);
                                        j -= 1;
                                    }
                                    max_index_used += 1;
                                    self.set_ti(i, max_index_used);
                                    self.set_tp(max_index_used, pm);
                                    self.set_tpar(max_index_used, u);

                                    let u1m = 0.5 * (u + self.tpar(self.ti(i - 1)));
                                    let p1m = curve_tool::value(c, u1m);
                                    if let Ok(dir) = GpDir2d::from_vec2d(&GpVec2d::new(
                                        pm.x() - p1.x(),
                                        pm.y() - p1.y(),
                                    )) {
                                        let l1m = GpLin2d::from_pnt_dir(p1, dir);
                                        t = l1m.distance(&p1m);
                                    }
                                    // `i--` (`gxx:239`), cancelled by the
                                    // `i++` of the `for` (`gxx:198`).
                                    i -= 1;
                                }
                            } else if t > new_deflection {
                                new_deflection = t;
                            }
                        }
                        i += 1;
                    }
                    let ratio_deflection = if new_deflection != 0.0 {
                        self.the_deflection / new_deflection
                    } else {
                        10.0
                    };
                    self.the_deflection = new_deflection;
                    self.nb_pnt_in = nbp;
                    if !((ratio_deflection < 3.0)
                        && (nb_passage_deflection < 3)
                        && (max_index_used < (self.the_max_nb_points - 2)))
                    {
                        break;
                    }
                }
            }

            self.the_deflection *= MAJORATION_DEFLECTION;
            self.my_box.enlarge(self.the_deflection);
        }
        self.closed_polygon = false;
    }

    /// `AutoIntersectionIsPossible()` (`gxx:268-282`).
    pub fn auto_intersection_is_possible(&self) -> bool {
        let v_ref = GpVec2d::new(
            self.tp(self.ti(2)).x() - self.tp(self.ti(1)).x(),
            self.tp(self.ti(2)).y() - self.tp(self.ti(1)).y(),
        );
        for i in 3..=self.nb_pnt_in {
            let v = GpVec2d::new(
                self.tp(self.ti(i)).x() - self.tp(self.ti(i - 1)).x(),
                self.tp(self.ti(i)).y() - self.tp(self.ti(i - 1)).y(),
            );
            if v.dot(&v_ref) < 0.0 {
                return true;
            }
        }
        false
    }

    /// `ApproxParamOnCurve(Aindex, TheParamOnLine)` (`gxx:287-310`).
    pub fn approx_param_on_curve(&self, a_index: usize, the_param_on_line: f64) -> f64 {
        let mut index = a_index;
        let mut param_on_line = the_param_on_line;
        if index > self.nb_pnt_in {
            // `gxx:291-293` prints to stdout; the port carries no debug output.
        }
        if (index == self.nb_pnt_in) && (param_on_line == 0.0) {
            index -= 1;
            param_on_line = 1.0;
        }
        if index == 0 {
            index = 1;
            param_on_line = 0.0;
        }
        let indexp1 = self.ti(index + 1);
        index = self.ti(index);

        let du = self.tpar(indexp1) - self.tpar(index);
        self.tpar(index) + param_on_line * du
    }
}

impl IntfPolygon2d for Geom2dIntPolygon2d {
    /// `Bounding()` (`Intf_Polygon2d.lxx:19-22`).
    fn bounding(&self) -> &BndBox2d {
        &self.my_box
    }

    /// `Closed()` (`lxx:32-35` / `Intf_Polygon2d.cxx:20-23`).
    fn closed(&self) -> bool {
        self.closed_polygon
    }

    /// `DeflectionOverEstimation()` (`lxx:20-23`).
    fn deflection_over_estimation(&self) -> f64 {
        self.the_deflection
    }

    /// `NbSegments()` (`lxx:38-41`).
    ///
    /// OCCT returns `Standard_Integer`, so a polygon built with no crossing
    /// region (`NbPntIn() == 0`, `ComputeWithBox` `gxx:255-256` sets
    /// `NbPntIn = nbp = 0`) yields `-1`: callers' `for (i = 1; i <= N; i++)`
    /// loops then iterate zero times. Saturating keeps that behaviour without
    /// an unsigned underflow.
    fn nb_segments(&self) -> usize {
        if self.closed_polygon {
            self.nb_pnt_in
        } else {
            self.nb_pnt_in.saturating_sub(1)
        }
    }

    /// `Segment(theIndex, theBegin, theEnd)` (`gxx:370-381`).
    fn segment(&self, index: usize, the_begin: &mut GpPnt2d, the_end: &mut GpPnt2d) {
        let mut ind = index;
        *the_begin = self.the_pnts[self.ti(index) - 1];
        if index >= self.nb_pnt_in {
            if !self.closed_polygon {
                panic!("Geom2dIntPolygon2d::Segment: out of range");
            }
            ind = 0;
        }
        *the_end = self.the_pnts[self.ti(ind + 1) - 1];
    }
}
