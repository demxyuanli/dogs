//! `AppDef_Variational::Approximate` — run the motor and turn the resulting
//! piecewise polynomial into B-spline poles.
//!
//! Source: `AppDef_Variational.cxx:489-662` (OCCT 8.0.0). The piecewise
//! polynomial to poles conversion is `Convert_CompPolynomialToPoles`
//! (`Convert_CompPolynomialToPoles.cxx`), already ported as
//! [`occt_core::bspl::comp_poly::CompPolynomialToPoles`].

use occt_core::bspl::comp_poly::CompPolynomialToPoles;
use occt_core::gp::{GpPnt, GpPnt2d};
use occt_math::vector::MathVector;

use super::Variational;
use crate::app_par_curves::multi_bspline_curve::MultiBSpCurve;
use crate::app_par_curves::multi_point::MultiPoint;

impl Variational {
    /// `Approximate()` (`cxx:489-662`). Returns `Err(())` where OCCT throws
    /// `StdFail_NotDone`.
    pub fn approximate(&mut self) -> Result<(), ()> {
        if !self.is_created {
            return Err(());
        }

        let (w_quadratic, w_quality) = self.smooth_criterion.weight();

        // `cxx:498`: `Ecarts` is indexed by point number.
        let mut ecarts = MathVector::new(self.first_point as usize, self.last_point as usize);

        let mut the_curve = match self.smooth_criterion.get_curve() {
            Some(c) => c,
            None => return Err(()),
        };

        // `cxx:509`.
        self.the_motor(w_quality, &mut the_curve, &mut ecarts);

        // `cxx:511-514`.
        if self.with_min_max && self.tolerance < self.max_error {
            self.adjusting(w_quadratic, w_quality, &mut the_curve, &mut ecarts);
        }

        // `cxx:518-660`: build the poles.
        let nb_elem;
        let mut coeffs;
        let mut nb_coeff;
        let mut intervals;
        let dim = self.dimension;
        {
            let cu = the_curve.borrow();
            nb_elem = cu.nb_elements() as i32;
            coeffs = vec![0.0f64; (self.max_segment * (self.max_degree + 1) * dim).max(1) as usize];
            nb_coeff = vec![0i32; self.max_segment.max(1) as usize];
            intervals = vec![0.0f64; nb_elem as usize + 2];
            intervals[1..=nb_elem as usize + 1]
                .copy_from_slice(&cu.knots()[1..=nb_elem as usize + 1]);
        }

        let mut curve = the_curve.borrow_mut();
        curve.get_polynom(&mut coeffs);
        for ii in 0..nb_elem as usize {
            nb_coeff[ii] = curve.degree(ii + 1) + 1;
        }
        drop(curve);

        // `cxx:551-556`: every span is mapped onto [-1, 1].
        let poly_intervals = vec![[-1.0f64, 1.0]; nb_elem.max(0) as usize];

        // `cxx:589`: the scalar-continuity constructor fills the per-span
        // array with `myNivCont`.
        let continuity = vec![self.niv_cont; nb_elem.max(0) as usize];

        let converter = CompPolynomialToPoles::new(
            nb_elem,
            dim,
            self.max_degree,
            &continuity,
            &nb_coeff[..nb_elem as usize],
            &coeffs,
            &poly_intervals,
            &intervals[1..],
        )
        .map_err(|_| ())?;

        if converter.done {
            // `cxx:595-658`.
            let mut tab_mu: Vec<MultiPoint> = Vec::with_capacity(converter.poles.len());
            for row in &converter.poles {
                let mut index = 0usize;
                let mut tab_p3d: Vec<GpPnt> = Vec::with_capacity(self.nb_p3d as usize);
                for _ in 0..self.nb_p3d {
                    tab_p3d.push(GpPnt::new(row[index], row[index + 1], row[index + 2]));
                    index += 3;
                }
                let mut tab_p2d: Vec<GpPnt2d> = Vec::with_capacity(self.nb_p2d as usize);
                for _ in 0..self.nb_p2d {
                    tab_p2d.push(GpPnt2d::new(row[index], row[index + 1]));
                    index += 2;
                }
                let multipoint = match (self.nb_p3d != 0, self.nb_p2d != 0) {
                    (true, true) => MultiPoint::from_3d_2d(&tab_p3d, &tab_p2d),
                    (true, false) => MultiPoint::from_3d(&tab_p3d),
                    _ => MultiPoint::from_2d(&tab_p2d),
                };
                tab_mu.push(multipoint);
            }

            // `cxx:598`: `myKnots` keeps the 1-based converter knot array.
            let mut knots = vec![0.0f64];
            knots.extend_from_slice(&converter.knots);
            self.knots = knots;

            self.mbsp_curve = MultiBSpCurve::from_tab(&tab_mu, &converter.knots, &converter.mults);
            self.is_done = true;
        }
        Ok(())
    }
}
