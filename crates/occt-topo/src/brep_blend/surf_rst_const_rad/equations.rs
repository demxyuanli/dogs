//! The three equations of `BRepBlend_SurfRstConstRad` and their derivatives.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_SurfRstConstRad.cxx:35-47`
//! (`t3dto2d`), `:90-108` (`Value`), `:113-170` (`Derivatives`), `:175-241`
//! (`Values`), `:600-604` (`Set(TypeSection)`) (OCCT 8.0.0).
//!
//! The vector `X` is `(u, v, w)`: `u, v` parameterise `surf` and `w`
//! parameterises the restriction `rst` on `surfrst`.

use occt_core::gp::{GpPnt, GpVec};
use occt_geom::curve::Curve;
use occt_math::{MathMatrix, MathVector};

use super::SurfRstConstRad;

/// `X.XYZ()`, the vector from the origin to `X` (`gp_Pnt::XYZ`).
pub(super) fn pnt_vec(p: &GpPnt) -> GpVec {
    GpVec::from_xyz(p.xyz())
}

/// `t3dto2d(a, b, A, B, C)` (`cxx:35-47`): express `A` on the (not
/// orthonormal) pair `B, C`, i.e. solve `A = a*B + b*C` in the plane they span.
pub(super) fn t3dto2d(a: &mut f64, b: &mut f64, va: &GpVec, vb: &GpVec, vc: &GpVec) {
    let ab = va.dot(vb);
    let ac = va.dot(vc);
    let bc = vb.dot(vc);
    let bb = vb.dot(vb);
    let cc = vc.dot(vc);
    let deno = bb * cc - bc * bc;
    *a = (ab * cc - ac * bc) / deno;
    *b = (ac * bb - ab * bc) / deno;
}

impl SurfRstConstRad {
    /// `Value(X, F)` (`cxx:90-108`).
    pub(super) fn compute_value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        let (pts, d1u1, d1v1) = self.surf.d1(x.value(1), x.value(2));
        self.pts = pts;
        self.ptrst = self.cons.d0(x.value(3));

        f.set_value(1, self.nplan.dot(&pnt_vec(&self.pts)) + self.the_d);
        f.set_value(2, self.nplan.dot(&pnt_vec(&self.ptrst)) + self.the_d);

        let ns = d1u1.crossed(&d1v1);
        let norm = self.nplan.crossed(&ns).magnitude();
        let ns = GpVec::zero().set_linear_form_2(
            self.nplan.dot(&ns) / norm,
            &self.nplan,
            -1.0 / norm,
            &ns,
        );
        let vref = GpVec::zero().set_linear_form_2(
            self.ray,
            &ns,
            1.0,
            &GpVec::from_pnts(&self.ptrst, &self.pts),
        );
        let vref = vref.divided(self.ray);
        f.set_value(3, (vref.square_magnitude() - 1.0) * self.ray * self.ray);
        true
    }

    /// `Derivatives(X, D)` (`cxx:113-170`).
    pub(super) fn compute_derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        let (pts, d1u1, d1v1, d2u1, d2v1, d2uv1) =
            self.surf.d2(x.value(1), x.value(2));
        self.pts = pts;
        let (ptrst, d1) = self.cons.d1(x.value(3));
        self.ptrst = ptrst;

        d.set_value(1, 1, self.nplan.dot(&d1u1));
        d.set_value(1, 2, self.nplan.dot(&d1v1));
        d.set_value(1, 3, 0.0);

        d.set_value(2, 1, 0.0);
        d.set_value(2, 2, 0.0);
        d.set_value(2, 3, self.nplan.dot(&d1));

        let ns = d1u1.crossed(&d1v1);
        let ncrossns = self.nplan.crossed(&ns);
        let norm = ncrossns.magnitude();
        let ndotns = self.nplan.dot(&ns);

        // `vref.SetLinearForm(ndotns, nplan, -1., ns); vref.Divide(norm);`
        // then `vref.SetLinearForm(ray, vref, gp_Vec(ptrst, pts))` (`cxx:135-138`).
        let vref = GpVec::zero().set_linear_form_2(ndotns, &self.nplan, -1.0, &ns);
        let vref = vref.divided(norm);
        let vref = GpVec::zero().set_linear_form_2(
            self.ray,
            &vref,
            1.0,
            &GpVec::from_pnts(&self.ptrst, &self.pts),
        );

        // Derivative by u1 (`cxx:140-151`).
        let temp = d2u1.crossed(&d1v1).added(&d1u1.crossed(&d2uv1));
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
        d.set_value(3, 1, 2.0 * resul.dot(&vref));

        // Derivative by v1 (`cxx:153-164`).
        let temp = d2uv1.crossed(&d1v1).added(&d1u1.crossed(&d2v1));
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
        d.set_value(3, 2, 2.0 * resul.dot(&vref));

        d.set_value(3, 3, -2.0 * d1.dot(&vref));
        true
    }

    /// `Values(X, F, D)` (`cxx:175-241`): the same two expressions as
    /// [`Self::compute_value`] / [`Self::compute_derivatives`], inlined the way
    /// OCCT inlines them. The third value keeps the numerically stabler form
    /// `(vref / ray)^2 - 1` scaled by `ray^2` (`cxx:207-208`, replacing the
    /// commented-out `vref.SquareMagnitude() - ray * ray`), exactly like
    /// `Value` (`cxx:103-104`); the derivatives below still differentiate the
    /// undivided `|vref|^2` (`cxx:220`).
    pub(super) fn compute_values(
        &mut self,
        x: &MathVector,
        f: &mut MathVector,
        d: &mut MathMatrix,
    ) -> bool {
        let (pts, d1u1, d1v1, d2u1, d2v1, d2uv1) =
            self.surf.d2(x.value(1), x.value(2));
        self.pts = pts;
        let (ptrst, d1) = self.cons.d1(x.value(3));
        self.ptrst = ptrst;

        f.set_value(1, self.nplan.dot(&pnt_vec(&self.pts)) + self.the_d);
        f.set_value(2, self.nplan.dot(&pnt_vec(&self.ptrst)) + self.the_d);

        d.set_value(1, 1, self.nplan.dot(&d1u1));
        d.set_value(1, 2, self.nplan.dot(&d1v1));
        d.set_value(1, 3, 0.0);

        d.set_value(2, 1, 0.0);
        d.set_value(2, 2, 0.0);
        d.set_value(2, 3, self.nplan.dot(&d1));

        let ns = d1u1.crossed(&d1v1);
        let ncrossns = self.nplan.crossed(&ns);
        let norm = ncrossns.magnitude();
        let ndotns = self.nplan.dot(&ns);

        let vref = GpVec::zero().set_linear_form_2(ndotns, &self.nplan, -1.0, &ns);
        let vref = vref.divided(norm);
        let vref = GpVec::zero().set_linear_form_2(
            self.ray,
            &vref,
            1.0,
            &GpVec::from_pnts(&self.ptrst, &self.pts),
        );

        // `temp = vref / ray; F(3) = (temp.SquareMagnitude() - 1) * ray * ray;`
        // (`cxx:207-208`).
        let temp = vref.divided(self.ray);
        f.set_value(3, (temp.square_magnitude() - 1.0) * self.ray * self.ray);

        // Derivative by u1 (`cxx:211-222`).
        let temp = d2u1.crossed(&d1v1).added(&d1u1.crossed(&d2uv1));
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
        d.set_value(3, 1, 2.0 * resul.dot(&vref));

        // Derivative by v1 (`cxx:225-236`).
        let temp = d2uv1.crossed(&d1v1).added(&d1u1.crossed(&d2v1));
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
        d.set_value(3, 2, 2.0 * resul.dot(&vref));

        d.set_value(3, 3, -2.0 * d1.dot(&vref));
        true
    }
}
