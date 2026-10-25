//! `Approx_MCurvesToBSpCurve` — glues a sequence of `AppParCurves_MultiCurve`
//! (all Bezier of the same degree) into one `AppParCurves_MultiBSpCurve`.
//!
//! Source: `ModelingData/TKGeomBase/Approx/Approx_MCurvesToBSpCurve.cxx:80-286`
//! (`Perform`, the only non-trivial member), with the class members
//! `myCurves` / `mySpline` / `myDone` from `Approx_MCurvesToBSpCurve.hxx:28-57`.
//!
//! Translation notes:
//! - `Convert_CompBezierCurvesToBSplineCurve` and
//!   `Convert_CompBezierCurves2dToBSplineCurve2d` are
//!   `occt_core::convert::CompBezierCurves{,2d}`.
//! - `BSplCLib::IncreaseDegree(deg, ThePoles, NoWeights, Points, NoWeights)`
//!   is `occt_core::convert::increase_degree_bezier` (Bezier syntax,
//!   `BSplCLib_CurveComputation.pxx:2080-2098`).
//! - `NCollection_Sequence<AppParCurves_MultiCurve> myCurves` becomes
//!   `Vec<MultiCurve>`; OCCT's 1-based `Value(i)` is `my_curves[i - 1]`.
//! - The `TheMults(i + 1)` / `TheMults(i)` tests index the multiplicity array
//!   1-based; `the_mults[i]` / `the_mults[i - 1]` are the same elements.

use crate::app_par_curves::multi_bspline_curve::MultiBSpCurve;
use crate::app_par_curves::multi_curve::MultiCurve;
use crate::app_par_curves::multi_point::MultiPoint;
use occt_core::convert::{
    increase_degree_bezier, increase_degree_bezier_2d, CompBezierCurves, CompBezierCurves2d,
};
use occt_core::gp::{GpPnt, GpPnt2d};

/// `Approx_MCurvesToBSpCurve` (`Approx_MCurvesToBSpCurve.hxx:28-57`).
#[derive(Debug, Clone)]
pub struct MCurvesToBSpCurve {
    /// `myCurves` (`hxx:44`).
    my_curves: Vec<MultiCurve>,
    /// `mySpline` (`hxx:46`).
    my_spline: MultiBSpCurve,
    /// `myDone` (`hxx:47`).
    my_done: bool,
}

impl Default for MCurvesToBSpCurve {
    fn default() -> Self {
        Self::new()
    }
}

impl MCurvesToBSpCurve {
    /// `Approx_MCurvesToBSpCurve()` (`cxx:57-60`).
    pub fn new() -> Self {
        Self {
            my_curves: Vec::new(),
            my_spline: MultiBSpCurve::new(),
            my_done: false,
        }
    }

    /// `Reset()` (`cxx:62-66`).
    pub fn reset(&mut self) {
        self.my_done = false;
        self.my_curves.clear();
    }

    /// `Append(MC)` (`cxx:68-71`).
    pub fn append(&mut self, mc: &MultiCurve) {
        self.my_curves.push(mc.clone());
    }

    /// `Perform()` (`cxx:73-76`): `Perform(myCurves)`.
    pub fn perform(&mut self) {
        let seq = self.my_curves.clone();
        self.perform_seq(&seq);
    }

    /// `Perform(TheSeq)` (`cxx:79-286`).
    pub fn perform_seq(&mut self, the_seq: &[MultiCurve]) {
        let nbcu = the_seq.len();
        if nbcu == 1 {
            let cu = &the_seq[0];
            let deg = cu.degree() as i32;
            let knots = [0.0f64, 1.0];
            let mults = [deg + 1, deg + 1];
            self.my_spline = MultiBSpCurve::from_curve(cu, &knots, &mults);
        } else if nbcu > 1 {
            let p = the_seq[nbcu - 1].value(1);
            let nb3d = p.nb_points() as i32;
            let nb2d = p.nb_points2d() as i32;

            let mut conv = CompBezierCurves::new(1.0e-4);
            let mut conv2d = CompBezierCurves2d::new(1.0e-4);
            if nb3d != 0 {
                for cu in the_seq.iter() {
                    conv.add_curve(&cu.curve3d(1));
                }
                conv.perform();
            } else if nb2d != 0 {
                for cu in the_seq.iter() {
                    conv2d.add_curve(&cu.curve2d(1 + nb3d));
                }
                conv2d.perform();
            }

            // Recuperation (`cxx:136-163`).
            let nbpolesspl = if nb3d != 0 {
                conv.nb_poles()
            } else if nb2d != 0 {
                conv2d.nb_poles()
            } else {
                0
            };

            let mut tab_mu: Vec<MultiPoint> =
                vec![MultiPoint::with_size(nb3d.max(0) as usize, nb2d.max(0) as usize); nbpolesspl];
            let the_knots: Vec<f64>;
            let the_mults: Vec<i32>;
            let deg: i32;
            if nb3d != 0 {
                let (k, m) = conv.knots_and_mults();
                the_knots = k.to_vec();
                the_mults = m.to_vec();
                deg = conv.degree();
                let poles_spl = conv.poles();
                for j in 0..nbpolesspl {
                    tab_mu[j].set_point(1, &poles_spl[j]);
                }
            } else if nb2d != 0 {
                let (k, m) = conv2d.knots_and_mults();
                the_knots = k.to_vec();
                the_mults = m.to_vec();
                deg = conv2d.degree();
                let poles_spl2d = conv2d.poles();
                for j in 0..nbpolesspl {
                    tab_mu[j].set_point2d(1 + nb3d, &poles_spl2d[j]);
                }
            } else {
                the_knots = Vec::new();
                the_mults = Vec::new();
                deg = 0;
            }

            // `cxx:180-271`.
            let mut kpol: i32 = 1;
            let mut kpoles3d: i32 = 1;
            let mut kpoles2d: i32 = 1;
            let thefirst = if nb3d != 0 { 1 } else { 2 };

            for (ii, cu) in the_seq.iter().enumerate() {
                let i = ii as i32 + 1;
                let mydegre = cu.degree() as i32;
                // `TheMults(i + 1)` (1-based) -> `the_mults[ii + 1]`.
                let mut last = if the_mults.len() > ii + 1 && the_mults[ii + 1] == deg {
                    deg + 1 // Continuite C0
                } else {
                    deg // Continuite C1
                };
                if i == nbcu as i32 {
                    last = deg + 1;
                }
                let first = if i == 1 {
                    1
                } else if ii >= 1 && (the_mults[ii] == deg - 1 || the_mults[ii] == deg) {
                    2
                } else {
                    1
                };

                for j in 2..=nb3d {
                    kpol = kpoles3d;
                    let the_poles = cu.curve3d(j);
                    let inc = deg - mydegre;
                    let points: Vec<GpPnt> = if inc > 0 {
                        increase_degree_bezier(deg, &the_poles)
                    } else {
                        the_poles.clone()
                    };
                    let mut k = first;
                    while k <= last {
                        tab_mu[(kpol - 1) as usize].set_point(j, &points[(k - 1) as usize]);
                        kpol += 1;
                        k += 1;
                    }
                }
                kpoles3d = kpol;

                for j in thefirst..=nb2d {
                    kpol = kpoles2d;
                    let the_poles2d = cu.curve2d(j + nb3d);
                    let inc = deg - mydegre;
                    let points2d: Vec<GpPnt2d> = if inc > 0 {
                        increase_degree_bezier_2d(deg, &the_poles2d)
                    } else {
                        the_poles2d.clone()
                    };
                    let mut k = first;
                    while k <= last {
                        tab_mu[(kpol - 1) as usize]
                            .set_point2d(j + nb3d, &points2d[(k - 1) as usize]);
                        kpol += 1;
                        k += 1;
                    }
                }
                kpoles2d = kpol;
            }

            self.my_spline = MultiBSpCurve::from_tab(&tab_mu, &the_knots, &the_mults);
        }
        self.my_done = true;
    }

    /// `Value()` (`cxx:288-291`).
    pub fn value(&self) -> &MultiBSpCurve {
        &self.my_spline
    }

    /// `ChangeValue()` (`cxx:293-296`).
    pub fn change_value(&mut self) -> &mut MultiBSpCurve {
        &mut self.my_spline
    }

    /// `IsDone()`: the OCCT header exposes `myDone` directly (`hxx:52`).
    pub fn is_done(&self) -> bool {
        self.my_done
    }
}
