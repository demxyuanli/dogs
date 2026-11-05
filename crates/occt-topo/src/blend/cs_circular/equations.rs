//! `Value` / `Derivatives` / `Values` of `BlendFunc_CSCircular`.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_CSCircular.cxx:255-397`
//! (OCCT 8.0.0).
//!
//! `Derivatives` and `Values` share one body except that `Values` also writes
//! the two equation values; both differ from `Value` in the projected normal
//! (`nsov`), which `Value` recomputes instead of keeping.

use occt_core::gp::{GpPnt2d, GpVec};
use occt_math::MathVector;

use super::constants::EPS;
use super::CsCircular;

impl CsCircular {
    /// `Value(X, F)` (`cxx:255-277`): the two equation values, without
    /// derivatives.
    pub(super) fn evaluate_value(&mut self, x: &MathVector) -> bool {
        let (pts, d1u1, d1v1) = self.surf.d1(x.value(1), x.value(2));
        self.pts = pts;
        let (ptc, _d1c) = self.curv.d1(self.prmc);
        self.ptc = ptc;

        self.e
            .set_value(1, self.nplan.dot(&GpVec::from_pnts(&self.ptc, &self.pts)));

        let ns = d1u1.crossed(&d1v1);
        let mut norm = self.nplan.crossed(&ns).magnitude();
        if norm < EPS {
            norm = 1.0;
        }

        let ns = GpVec::zero().set_linear_form_2(
            self.nplan.dot(&ns) / norm,
            &self.nplan,
            -1.0 / norm,
            &ns,
        );
        let vref = ns
            .multiplied_scalar(self.ray)
            .added(&GpVec::from_pnts(&self.ptc, &self.pts));

        self.e.set_value(2, vref.square_magnitude() - self.ray * self.ray);

        self.pt2d = GpPnt2d::new(x.value(1), x.value(2));
        true
    }

    /// `Derivatives(X, D)` (`cxx:282-337`): the 2x2 jacobian, without the
    /// equation values.
    pub(super) fn evaluate_derivatives(&mut self, x: &MathVector) -> bool {
        let (pts, d1u1, d1v1, d2u1, d2v1, d2uv1) = self.surf.d2(x.value(1), x.value(2));
        self.pts = pts;
        let (ptc, _d1c) = self.curv.d1(self.prmc);
        self.ptc = ptc;

        self.dedx.set_value(1, 1, self.nplan.dot(&d1u1));
        self.dedx.set_value(1, 2, self.nplan.dot(&d1v1));

        let ns = d1u1.crossed(&d1v1);
        let ncrossns = self.nplan.crossed(&ns);
        let mut norm = ncrossns.magnitude();
        if norm < EPS {
            norm = 1.0;
        }
        let ndotns = self.nplan.dot(&ns);

        // `nsov` is the projected normal; `resul` keeps using the raw `ns`
        // (`cxx:305-312`).
        let nsov = GpVec::zero().set_linear_form_2(
            self.nplan.dot(&ns) / norm,
            &self.nplan,
            -1.0 / norm,
            &ns,
        );
        let vref = nsov
            .multiplied_scalar(self.ray)
            .added(&GpVec::from_pnts(&self.ptc, &self.pts));

        // Derived in relation to u (`cxx:309-317`).
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
        self.dedx.set_value(2, 1, 2.0 * resul.dot(&vref));

        // Derived in relation to v (`cxx:322-332`).
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
        self.dedx.set_value(2, 2, 2.0 * resul.dot(&vref));

        self.pt2d = GpPnt2d::new(x.value(1), x.value(2));
        true
    }

    /// `Values(X, F, D)` (`cxx:339-397`): the equation values and the jacobian
    /// in one pass.
    pub(super) fn evaluate_values(&mut self, x: &MathVector) -> bool {
        let (pts, d1u1, d1v1, d2u1, d2v1, d2uv1) = self.surf.d2(x.value(1), x.value(2));
        self.pts = pts;
        let (ptc, _d1c) = self.curv.d1(self.prmc);
        self.ptc = ptc;

        let ns = d1u1.crossed(&d1v1);
        let ncrossns = self.nplan.crossed(&ns);
        let mut norm = ncrossns.magnitude();
        if norm < EPS {
            norm = 1.0;
        }
        let ndotns = self.nplan.dot(&ns);

        let nsov = GpVec::zero().set_linear_form_2(ndotns / norm, &self.nplan, -1.0 / norm, &ns);
        let vref = nsov
            .multiplied_scalar(self.ray)
            .added(&GpVec::from_pnts(&self.ptc, &self.pts));

        self.e
            .set_value(1, self.nplan.dot(&GpVec::from_pnts(&self.ptc, &self.pts)));
        self.e.set_value(2, vref.square_magnitude() - self.ray * self.ray);

        self.dedx.set_value(1, 1, self.nplan.dot(&d1u1));
        self.dedx.set_value(1, 2, self.nplan.dot(&d1v1));

        // Derived in relation to u (`cxx:368-379`).
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
        self.dedx.set_value(2, 1, 2.0 * resul.dot(&vref));

        // Derived in relation to v (`cxx:382-392`).
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
        self.dedx.set_value(2, 2, 2.0 * resul.dot(&vref));

        self.pt2d = GpPnt2d::new(x.value(1), x.value(2));
        true
    }
}
