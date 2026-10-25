//! `Approx_ComputeLine::Compute` and `Approx_ComputeLine::ComputeCurve` —
//! the degree search and the tangent-constrained interpolator.
//!
//! Source: `ModelingData/TKGeomBase/Approx/Approx_ComputeLine.gxx:1324-1690`.
//!
//! Translation notes:
//! - `Compute`'s `TheTol3d` / `TheTol2d` out-parameters are never read by
//!   `Perform` (`gxx:1116`, `:1163`, `:1183` discard them), so they stay local
//!   and the return value is `(done, indbad)`.
//! - `AppParCurves_MultiCurve mySCU(deg + 1)` is assigned immediately from
//!   `BezierValue()` (`gxx:1357-1362`, `:1664`), so the port skips the
//!   intermediate default construction.
//! - The `nbp == 2` tangent-parallelism block is inside `#ifdef OCCT_DEBUG`
//!   (`gxx:1499-1557`); only the pole interpolation it guards is ported.
//! - `RealLast()` is `f64::MAX` (`Standard_Real.hxx:99`).

use occt_core::gp::GpPnt;
use occt_core::precision;
use occt_math::MathVector;

use crate::approx::line_tool::LineTool;
use crate::app_par_curves::constraint::Constraint;
use crate::app_par_curves::gradient::Gradient;
use crate::app_par_curves::least_square::LeastSquare;
use crate::app_par_curves::multi_curve::MultiCurve;
use crate::app_par_curves::multi_point::MultiPoint;
use crate::approx::compute_line::check_multi_curve::check_multi_curve;
use crate::approx::compute_line::ComputeLine;

impl ComputeLine {
    /// `Compute(Line, fpt, lpt, Para, TheTol3d, TheTol2d, indbad)`
    /// (`gxx:1324-1441`).
    pub(crate) fn compute<L: LineTool>(
        &mut self,
        line: &L,
        fpt: i32,
        lpt: i32,
        para: &mut MathVector,
    ) -> (bool, i32) {
        let mut indbad = 0i32;
        let nbp = lpt - fpt + 1;

        let mut par_sav = MathVector::new(para.lower(), para.upper());
        for i in para.lower()..=para.upper() {
            par_sav.set_value(i, para.value(i));
        }

        let mut mdegmax = self.my_degremax;
        if nbp < mdegmax + 5 && self.my_cut {
            mdegmax = nbp - 5;
        }
        if mdegmax < self.my_degremin {
            mdegmax = self.my_degremin;
        }

        self.current_tol3d = f64::MAX;
        self.current_tol2d = f64::MAX;

        let mut deg = (nbp - 1).min(self.my_degremin);
        while deg <= mdegmax {
            let (mydone, my_scu, the_tol3d, the_tol2d) = if self.my_squares {
                let mut sq = LeastSquare::with_parameters(
                    line,
                    fpt,
                    lpt,
                    self.my_first_c,
                    self.my_last_c,
                    para,
                    deg + 1,
                );
                let mydone = sq.is_done();
                let scu = sq.bezier_value();
                let (_, the_tol3d, the_tol2d) = sq.error();
                (mydone, scu, the_tol3d, the_tol2d)
            } else {
                let grad = Gradient::new(
                    line,
                    fpt,
                    lpt,
                    &self.my_constraints,
                    para,
                    deg,
                    self.my_tol3d,
                    self.my_tol2d,
                    self.my_iter_max,
                );
                let mydone = grad.is_done();
                let scu = grad.value();
                if scu.nb_curves() == 0 {
                    deg += 1;
                    continue;
                }
                let the_tol3d = grad.max_error3d();
                let the_tol2d = grad.max_error2d();
                (mydone, scu, the_tol3d, the_tol2d)
            };

            let mut uu1 = para.value(para.lower());
            let mut restau = false;
            for i in para.lower() + 1..=para.upper() {
                let uu2 = para.value(i);
                if uu2 <= uu1 {
                    restau = true;
                    break;
                }
                uu1 = uu2;
            }
            if restau {
                for i in para.lower()..=para.upper() {
                    para.set_value(i, par_sav.value(i));
                }
            }

            if mydone && the_tol3d <= self.my_tol3d && the_tol2d <= self.my_tol2d {
                // Stockage de la multicurve approximee.
                self.tolreached = true;
                if !check_multi_curve(&my_scu, line, fpt, lpt, &mut indbad) {
                    return (false, indbad);
                }
                self.my_multi_curves.push(my_scu);
                // Storage of parameters of the approximated MultiLine part.
                let the_par: Vec<f64> = (para.lower()..=para.upper())
                    .map(|i| para.value(i))
                    .collect();
                self.my_par.push((para.lower() as i32, the_par));
                self.tol_ers3d.push(the_tol3d);
                self.tol_ers2d.push(the_tol2d);
                return (true, indbad);
            }

            if the_tol3d <= self.current_tol3d && the_tol2d <= self.current_tol2d {
                self.the_multi_curve = my_scu;
                self.current_tol3d = the_tol3d;
                self.current_tol2d = the_tol2d;
                let the_par: Vec<f64> = (para.lower()..=para.upper())
                    .map(|i| para.value(i))
                    .collect();
                self.my_parameters = Some((para.lower() as i32, the_par));
            }

            deg += 1;
        }

        (false, indbad)
    }

    /// `ComputeCurve(Line, firstpt, lastpt)` (`gxx:1443-1690`).
    pub(crate) fn compute_curve<L: LineTool>(
        &mut self,
        line: &L,
        firstpt: i32,
        lastpt: i32,
    ) -> bool {
        let myfirstpt = firstpt;
        let mylastpt = lastpt;
        let nbp = lastpt - firstpt + 1;
        let para = self.compute_parameters(line, firstpt, lastpt);

        let nb_p3d = line.nb_p3d();
        let nb_p2d = line.nb_p2d();

        let (tab_p1, tab_p12d) = line_values(line, myfirstpt, nb_p3d, nb_p2d);
        let (tab_p2, tab_p22d) = line_values(line, mylastpt, nb_p3d, nb_p2d);

        if nbp == 2 {
            let mut my_scu = MultiCurve::with_nb_poles((self.my_degremin + 1) as usize);
            let deg = self.my_degremin;

            if nb_p3d != 0 && nb_p2d != 0 {
                my_scu.set_value(1, &MultiPoint::from_3d_2d(&tab_p1, &tab_p12d));
                my_scu.set_value(deg + 1, &MultiPoint::from_3d_2d(&tab_p2, &tab_p22d));
                for i in 2..=deg {
                    let tab_p = interpolate3d(&tab_p1, &tab_p2, i, deg);
                    let tab_p2d = interpolate2d(&tab_p12d, &tab_p22d, i, deg);
                    my_scu.set_value(i, &MultiPoint::from_3d_2d(&tab_p, &tab_p2d));
                }
            } else if nb_p3d != 0 {
                my_scu.set_value(1, &MultiPoint::from_3d(&tab_p1));
                my_scu.set_value(deg + 1, &MultiPoint::from_3d(&tab_p2));
                for i in 2..=deg {
                    let tab_p = interpolate3d(&tab_p1, &tab_p2, i, deg);
                    my_scu.set_value(i, &MultiPoint::from_3d(&tab_p));
                }
            } else if nb_p2d != 0 {
                my_scu.set_value(1, &MultiPoint::from_2d(&tab_p12d));
                my_scu.set_value(deg + 1, &MultiPoint::from_2d(&tab_p22d));
                for i in 2..=deg {
                    let tab_p2d = interpolate2d(&tab_p12d, &tab_p22d, i, deg);
                    my_scu.set_value(i, &MultiPoint::from_2d(&tab_p2d));
                }
            }

            let mydone = true;
            // Stockage de la multicurve approximee.
            self.tolreached = true;
            self.my_multi_curves.push(my_scu);
            let the_par: Vec<f64> = (para.lower()..=para.upper())
                .map(|i| para.value(i))
                .collect();
            self.my_par.push((para.lower() as i32, the_par));
            self.tol_ers3d.push(precision::CONFUSION);
            self.tol_ers2d.push(precision::PCONFUSION);
            return mydone;
        }

        // with the tangents.
        let deg = nbp + 1;
        let cons = Constraint::TangencyPoint;
        let nb_vars = (nb_p3d * 3 + nb_p2d * 2) as usize;
        let mut v1 = MathVector::new(1, nb_vars);
        let mut v2 = MathVector::new(1, nb_vars);
        self.first_tangency_vector(line, myfirstpt, &mut v1);
        let mut lambda1 = self.search_first_lambda(line, &para, &v1, myfirstpt);

        self.last_tangency_vector(line, mylastpt, &mut v2);
        let mut lambda2 = self.search_last_lambda(line, &para, &v2, mylastpt);

        let mut lsq = LeastSquare::with_parameters(
            line,
            myfirstpt,
            mylastpt,
            cons,
            cons,
            &para,
            deg + 1,
        );

        lambda1 /= deg as f64;
        lambda2 /= deg as f64;
        lsq.perform_tangents(&para, &v1, &v2, lambda1, lambda2);
        let mydone = lsq.is_done();
        let my_scu = lsq.bezier_value();

        if mydone {
            let (_, the_tol3d, the_tol2d) = lsq.error();
            // Stockage de la multicurve approximee.
            self.tolreached = true;
            self.my_multi_curves.push(my_scu);
            let the_par: Vec<f64> = (para.lower()..=para.upper())
                .map(|i| para.value(i))
                .collect();
            self.my_par.push((para.lower() as i32, the_par));
            self.tol_ers3d.push(the_tol3d);
            self.tol_ers2d.push(the_tol2d);
            return true;
        }
        mydone
    }
}

/// `LineTool::Value(Line, i, tabP, tabP2d)` dispatched on the column counts
/// (`gxx:1473-1492`).
fn line_values<L: LineTool>(
    line: &L,
    index: i32,
    nb_p3d: i32,
    nb_p2d: i32,
) -> (Vec<GpPnt>, Vec<occt_core::gp::GpPnt2d>) {
    if nb_p3d != 0 && nb_p2d != 0 {
        line.value(index)
    } else if nb_p2d != 0 {
        (Vec::new(), line.value_2d(index))
    } else if nb_p3d != 0 {
        (line.value_3d(index), Vec::new())
    } else {
        (Vec::new(), Vec::new())
    }
}

/// `tabP(j).SetXYZ(P1.XYZ() + (i - 1) * (P2.XYZ() - P1.XYZ()) / mydegremin)`
/// (`gxx:1575-1578`).
fn interpolate3d(tab_p1: &[GpPnt], tab_p2: &[GpPnt], i: i32, deg: i32) -> Vec<GpPnt> {
    (0..tab_p1.len())
        .map(|j| {
            let p1 = &tab_p1[j];
            let p2 = &tab_p2[j];
            GpPnt::new(
                p1.x() + (i - 1) as f64 * (p2.x() - p1.x()) / deg as f64,
                p1.y() + (i - 1) as f64 * (p2.y() - p1.y()) / deg as f64,
                p1.z() + (i - 1) as f64 * (p2.z() - p1.z()) / deg as f64,
            )
        })
        .collect()
}

/// `tabP2d(j).SetXY(P12d.XY() + (i - 1) * (P22d.XY() - P12d.XY()) / mydegremin)`
/// (`gxx:1580-1582`).
fn interpolate2d(
    tab_p12d: &[occt_core::gp::GpPnt2d],
    tab_p22d: &[occt_core::gp::GpPnt2d],
    i: i32,
    deg: i32,
) -> Vec<occt_core::gp::GpPnt2d> {
    (0..tab_p12d.len())
        .map(|j| {
            let p1 = &tab_p12d[j];
            let p2 = &tab_p22d[j];
            occt_core::gp::GpPnt2d::new(
                p1.x() + (i - 1) as f64 * (p2.x() - p1.x()) / deg as f64,
                p1.y() + (i - 1) as f64 * (p2.y() - p1.y()) / deg as f64,
            )
        })
        .collect()
}





