//! The four equations of `BlendFunc_Ruled` and their derivatives.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_Ruled.cxx:178-201`
//! (`Value`), `:203-294` (`Derivatives`), `:296-393` (`Values`) (OCCT 8.0.0).
//!
//! Translation notes:
//! - `Value` only needs the first derivatives of the two surfaces
//!   (`cxx:179-180`), so it is kept apart from the second-order body.
//! - `Derivatives` and `Values` repeat the same body in OCCT, `Values` adding
//!   the four residuals; the port writes it once with `f: Option<&mut
//!   MathVector>` and both entry points call it.
//! - There is no `Eps` guard on `norm1` / `norm2` (`cxx:189-194`,
//!   `:307-318`), unlike `ConstRad` / `EvolRad`: a projection of the surface
//!   normal onto the section plane that vanishes divides by zero exactly as
//!   OCCT does.

use occt_core::gp::GpVec;
use occt_math::{MathMatrix, MathVector};

use super::{pnt_vec, Ruled};

impl Ruled {
    /// `Value(X, F)` (`cxx:178-201`).
    pub(super) fn compute_value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        let (pts1, d1u1, d1v1) = self.surf1.d1(x.value(1), x.value(2));
        let (pts2, d1u2, d1v2) = self.surf2.d1(x.value(3), x.value(4));
        self.pts1 = pts1;
        self.pts2 = pts2;

        let temp = GpVec::from_pnts(&self.pts1, &self.pts2);

        let ns1 = d1u1.crossed(&d1v1);
        let ns2 = d1u2.crossed(&d1v2);

        let norm1 = self.nplan.crossed(&ns1).magnitude();
        let norm2 = self.nplan.crossed(&ns2).magnitude();
        let ns1 = GpVec::zero().set_linear_form_2(
            self.nplan.dot(&ns1) / norm1,
            &self.nplan,
            -1.0 / norm1,
            &ns1,
        );
        let ns2 = GpVec::zero().set_linear_form_2(
            self.nplan.dot(&ns2) / norm2,
            &self.nplan,
            -1.0 / norm2,
            &ns2,
        );

        f.set_value(1, self.nplan.dot(&pnt_vec(&self.pts1)) + self.the_d);
        f.set_value(2, self.nplan.dot(&pnt_vec(&self.pts2)) + self.the_d);
        f.set_value(3, temp.dot(&ns1));
        f.set_value(4, temp.dot(&ns2));

        true
    }

    /// `Derivatives(X, D)` (`cxx:203-294`) and `Values(X, F, D)`
    /// (`cxx:296-393`). `f` is `None` for `Derivatives`.
    pub(super) fn compute_values(
        &mut self,
        x: &MathVector,
        f: Option<&mut MathVector>,
        d: &mut MathMatrix,
    ) -> bool {
        let (pts1, d1u1, d1v1, d2u1, d2v1, d2uv1) =
            self.surf1.d2(x.value(1), x.value(2));
        let (pts2, d1u2, d1v2, d2u2, d2v2, d2uv2) =
            self.surf2.d2(x.value(3), x.value(4));
        self.pts1 = pts1;
        self.pts2 = pts2;

        let p1p2 = GpVec::from_pnts(&self.pts1, &self.pts2);

        let ns1 = d1u1.crossed(&d1v1);
        let ns2 = d1u2.crossed(&d1v2);
        let ncrossns1 = self.nplan.crossed(&ns1);
        let ncrossns2 = self.nplan.crossed(&ns2);
        let norm1 = ncrossns1.magnitude();
        let norm2 = ncrossns2.magnitude();

        let ndotns1 = self.nplan.dot(&ns1);
        let ndotns2 = self.nplan.dot(&ns2);

        let nor1 = GpVec::zero().set_linear_form_2(
            ndotns1 / norm1,
            &self.nplan,
            -1.0 / norm1,
            &ns1,
        );
        let nor2 = GpVec::zero().set_linear_form_2(
            ndotns2 / norm2,
            &self.nplan,
            -1.0 / norm2,
            &ns2,
        );

        // `F(1)` .. `F(4)` (`cxx:324-327`); `Values` only.
        if let Some(f) = f {
            f.set_value(1, self.nplan.dot(&pnt_vec(&self.pts1)) + self.the_d);
            f.set_value(2, self.nplan.dot(&pnt_vec(&self.pts2)) + self.the_d);
            f.set_value(3, p1p2.dot(&nor1));
            f.set_value(4, p1p2.dot(&nor2));
        }

        d.set_value(1, 1, self.nplan.dot(&d1u1));
        d.set_value(1, 2, self.nplan.dot(&d1v1));
        d.set_value(1, 3, 0.0);
        d.set_value(1, 4, 0.0);

        d.set_value(2, 1, 0.0);
        d.set_value(2, 2, 0.0);
        d.set_value(2, 3, self.nplan.dot(&d1u2));
        d.set_value(2, 4, self.nplan.dot(&d1v2));

        // Derivative of nor1 with respect to u1 (`cxx:240-250`).
        let mut temp = d2u1.crossed(&d1v1).added(&d1u1.crossed(&d2uv1));
        let mut grosterme = ncrossns1.dot(&self.nplan.crossed(&temp)) / norm1 / norm1;
        let mut resul = GpVec::zero().set_linear_form_3(
            -(grosterme * ndotns1 - self.nplan.dot(&temp)) / norm1,
            &self.nplan,
            grosterme / norm1,
            &ns1,
            -1.0 / norm1,
            &temp,
        );

        d.set_value(3, 1, -d1u1.dot(&nor1) + p1p2.dot(&resul));

        // Derivative with respect to v1 (`cxx:253-260`).
        temp = d2uv1.crossed(&d1v1).added(&d1u1.crossed(&d2v1));
        grosterme = ncrossns1.dot(&self.nplan.crossed(&temp)) / norm1 / norm1;
        resul = GpVec::zero().set_linear_form_3(
            -(grosterme * ndotns1 - self.nplan.dot(&temp)) / norm1,
            &self.nplan,
            grosterme / norm1,
            &ns1,
            -1.0 / norm1,
            &temp,
        );

        d.set_value(3, 2, -d1v1.dot(&nor1) + p1p2.dot(&resul));

        d.set_value(3, 3, d1u2.dot(&nor1));
        d.set_value(3, 4, d1v2.dot(&nor1));

        d.set_value(4, 1, -d1u1.dot(&nor2));
        d.set_value(4, 2, -d1v1.dot(&nor2));

        // Derivative of nor2 with respect to u2 (`cxx:272-280`).
        temp = d2u2.crossed(&d1v2).added(&d1u2.crossed(&d2uv2));
        grosterme = ncrossns2.dot(&self.nplan.crossed(&temp)) / norm2 / norm2;
        resul = GpVec::zero().set_linear_form_3(
            -(grosterme * ndotns2 - self.nplan.dot(&temp)) / norm2,
            &self.nplan,
            grosterme / norm2,
            &ns2,
            -1.0 / norm2,
            &temp,
        );

        d.set_value(4, 3, d1u2.dot(&nor2) + p1p2.dot(&resul));

        // Derivative with respect to v2 (`cxx:283-290`).
        temp = d2uv2.crossed(&d1v2).added(&d1u2.crossed(&d2v2));
        grosterme = ncrossns2.dot(&self.nplan.crossed(&temp)) / norm2 / norm2;
        resul = GpVec::zero().set_linear_form_3(
            -(grosterme * ndotns2 - self.nplan.dot(&temp)) / norm2,
            &self.nplan,
            grosterme / norm2,
            &ns2,
            -1.0 / norm2,
            &temp,
        );

        d.set_value(4, 4, d1v2.dot(&nor2) + p1p2.dot(&resul));

        true
    }
}
