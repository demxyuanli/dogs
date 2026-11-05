//! `Value` / `Derivatives` / `Values` of `BlendFunc_CSConstRad`.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_CSConstRad.cxx:228-375`
//! (OCCT 8.0.0).
//!
//! OCCT writes the three bodies separately and they do not round the same way:
//! `Value` normalises the projected normal term by term
//! (`nplan.Dot(ns) / norm * nplan` and `-1. / norm * ns`) while `Derivatives` /
//! `Values` build `nplan.Dot(ns) * nplan - ns` first and divide afterwards. The
//! three bodies are kept apart here for that reason.

use occt_core::gp::{GpPnt2d, GpVec};
use occt_math::MathVector;

use super::{pnt_vec, CsConstRad};

impl CsConstRad {
    /// `Value(X, F)` (`cxx:228-247`): the three equation values, without
    /// derivatives.
    pub(super) fn evaluate_value(&mut self, x: &MathVector) -> bool {
        let (pts, d1u1, d1v1) = self.surf.d1(x.value(1), x.value(2));
        self.pts = pts;
        self.ptc = self.curv.value(x.value(3));

        self.e
            .set_value(1, self.nplan.dot(&pnt_vec(&self.pts)) + self.the_d);
        self.e
            .set_value(2, self.nplan.dot(&pnt_vec(&self.ptc)) + self.the_d);

        let ns = d1u1.crossed(&d1v1);
        let norm = self.nplan.crossed(&ns).magnitude();
        let ns = GpVec::zero().set_linear_form_2(
            self.nplan.dot(&ns) / norm,
            &self.nplan,
            -1.0 / norm,
            &ns,
        );
        // `vref.SetLinearForm(ray, ns, gp_Vec(ptc, pts))` (`cxx:240-241`).
        let vref = ns
            .multiplied_scalar(self.ray)
            .added(&GpVec::from_pnts(&self.ptc, &self.pts));

        self.e.set_value(3, vref.square_magnitude() - self.ray * self.ray);

        self.pt2d = GpPnt2d::new(x.value(1), x.value(2));
        true
    }

    /// `Derivatives(X, D)` (`cxx:251-306`): the 3x3 jacobian, without the
    /// equation values.
    pub(super) fn evaluate_derivatives(&mut self, x: &MathVector) -> bool {
        let (pts, d1u1, d1v1, d2u1, d2v1, d2uv1) = self.surf.d2(x.value(1), x.value(2));
        self.pts = pts;
        let (ptc, d1) = self.curv.d1(x.value(3));
        self.ptc = ptc;

        self.dedx.set_value(1, 1, self.nplan.dot(&d1u1));
        self.dedx.set_value(1, 2, self.nplan.dot(&d1v1));
        self.dedx.set_value(1, 3, 0.0);

        self.dedx.set_value(2, 1, 0.0);
        self.dedx.set_value(2, 2, 0.0);
        self.dedx.set_value(2, 3, self.nplan.dot(&d1));

        let ns = d1u1.crossed(&d1v1);
        let ncrossns = self.nplan.crossed(&ns);
        let norm = ncrossns.magnitude();
        let ndotns = self.nplan.dot(&ns);

        // `vref.SetLinearForm(ndotns, nplan, -1., ns); vref.Divide(norm);`
        // (`cxx:274-275`).
        let vref = GpVec::zero().set_linear_form_2(ndotns, &self.nplan, -1.0, &ns);
        let vref = vref
            .divide(norm)
            .multiplied_scalar(self.ray)
            .added(&GpVec::from_pnts(&self.ptc, &self.pts));

        // Derived in relation to u (`cxx:278-288`).
        let mut temp = d2u1.crossed(&d1v1).added(&d1u1.crossed(&d2uv1));
        let grosterme = ncrossns.dot(&self.nplan.crossed(&temp)) / norm / norm;
        let resul = GpVec::zero()
            .set_linear_form_3(
                -self.ray / norm * (grosterme * ndotns - self.nplan.dot(&temp)),
                &self.nplan,
                self.ray * grosterme / norm,
                &ns,
                -self.ray / norm,
                &temp,
            )
            .added(&d1u1);
        self.dedx.set_value(3, 1, 2.0 * resul.dot(&vref));

        // Derived in relation to v (`cxx:291-301`).
        temp = d2uv1.crossed(&d1v1).added(&d1u1.crossed(&d2v1));
        let grosterme = ncrossns.dot(&self.nplan.crossed(&temp)) / norm / norm;
        let resul = GpVec::zero()
            .set_linear_form_3(
                -self.ray / norm * (grosterme * ndotns - self.nplan.dot(&temp)),
                &self.nplan,
                self.ray * grosterme / norm,
                &ns,
                -self.ray / norm,
                &temp,
            )
            .added(&d1v1);
        self.dedx.set_value(3, 2, 2.0 * resul.dot(&vref));

        self.dedx.set_value(3, 3, -2.0 * d1.dot(&vref));

        self.pt2d = GpPnt2d::new(x.value(1), x.value(2));
        true
    }

    /// `Values(X, F, D)` (`cxx:312-375`): the equation values and the jacobian
    /// in one pass.
    pub(super) fn evaluate_values(&mut self, x: &MathVector) -> bool {
        let (pts, d1u1, d1v1, d2u1, d2v1, d2uv1) = self.surf.d2(x.value(1), x.value(2));
        self.pts = pts;
        let (ptc, d1) = self.curv.d1(x.value(3));
        self.ptc = ptc;

        self.e
            .set_value(1, self.nplan.dot(&pnt_vec(&self.pts)) + self.the_d);
        self.e
            .set_value(2, self.nplan.dot(&pnt_vec(&self.ptc)) + self.the_d);

        self.dedx.set_value(1, 1, self.nplan.dot(&d1u1));
        self.dedx.set_value(1, 2, self.nplan.dot(&d1v1));
        self.dedx.set_value(1, 3, 0.0);

        self.dedx.set_value(2, 1, 0.0);
        self.dedx.set_value(2, 2, 0.0);
        self.dedx.set_value(2, 3, self.nplan.dot(&d1));

        let ns = d1u1.crossed(&d1v1);
        let ncrossns = self.nplan.crossed(&ns);
        let norm = ncrossns.magnitude();
        let ndotns = self.nplan.dot(&ns);

        let vref = GpVec::zero().set_linear_form_2(ndotns, &self.nplan, -1.0, &ns);
        let vref = vref
            .divide(norm)
            .multiplied_scalar(self.ray)
            .added(&GpVec::from_pnts(&self.ptc, &self.pts));

        self.e.set_value(3, vref.square_magnitude() - self.ray * self.ray);

        // Derived in relation to u (`cxx:344-355`).
        let mut temp = d2u1.crossed(&d1v1).added(&d1u1.crossed(&d2uv1));
        let grosterme = ncrossns.dot(&self.nplan.crossed(&temp)) / norm / norm;
        let resul = GpVec::zero()
            .set_linear_form_3(
                -self.ray / norm * (grosterme * ndotns - self.nplan.dot(&temp)),
                &self.nplan,
                self.ray * grosterme / norm,
                &ns,
                -self.ray / norm,
                &temp,
            )
            .added(&d1u1);
        self.dedx.set_value(3, 1, 2.0 * resul.dot(&vref));

        // Derived in relation to v (`cxx:358-367`).
        temp = d2uv1.crossed(&d1v1).added(&d1u1.crossed(&d2v1));
        let grosterme = ncrossns.dot(&self.nplan.crossed(&temp)) / norm / norm;
        let resul = GpVec::zero()
            .set_linear_form_3(
                -self.ray / norm * (grosterme * ndotns - self.nplan.dot(&temp)),
                &self.nplan,
                self.ray * grosterme / norm,
                &ns,
                -self.ray / norm,
                &temp,
            )
            .added(&d1v1);
        self.dedx.set_value(3, 2, 2.0 * resul.dot(&vref));

        self.dedx.set_value(3, 3, -2.0 * d1.dot(&vref));

        self.pt2d = GpPnt2d::new(x.value(1), x.value(2));
        true
    }
}
