//! The three equations of `BRepBlend_SurfRstEvolRad` and their derivatives.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_SurfRstEvolRad.cxx:39-49`
//! (`t3dto2d`), `:153-171` (`Value`), `:175-231` (`Derivatives`), `:235-296`
//! (`Values`) (OCCT 8.0.0).
//!
//! The vector `X` is `(u, v, w)`: `u, v` parameterise `surf` and `w`
//! parameterises the restriction `rst` on `surfrst`. `ray` / `dray` are the
//! current radius and its derivative, already signed by `sg1` in `Set(Param)`
//! (`cxx:321-324`) or in `Section` (`cxx:867-869`).

use occt_core::gp::{GpPnt, GpVec};
use occt_geom::curve::Curve;
use occt_math::{MathMatrix, MathVector};

use super::SurfRstEvolRad;

/// `X.XYZ()`, the vector from the origin to `X` (`gp_Pnt::XYZ`).
pub(super) fn pnt_vec(p: &GpPnt) -> GpVec {
    GpVec::from_xyz(p.xyz())
}

/// `t3dto2d(a, b, A, B, C)` (`cxx:39-49`): express `A` on the (not
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

impl SurfRstEvolRad {
    /// `Value(X, F)` (`cxx:153-171`).
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
        // `vref.SetLinearForm(ray, ns, gp_Vec(ptrst, pts));`
        // `F(3) = vref.SquareMagnitude() - ray * ray;` (`cxx:167-168`).
        let vref = GpVec::zero().set_linear_form_2(
            self.ray,
            &ns,
            1.0,
            &GpVec::from_pnts(&self.ptrst, &self.pts),
        );
        f.set_value(3, vref.square_magnitude() - self.ray * self.ray);
        true
    }

    /// `Derivatives(X, D)` (`cxx:175-231`).
    pub(super) fn compute_derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        let (pts, d1u1, d1v1, d2u1, d2v1, d2uv1) = self.surf.d2(x.value(1), x.value(2));
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

        let vref = GpVec::zero().set_linear_form_2(ndotns, &self.nplan, -1.0, &ns);
        let vref = vref.divided(norm);
        let vref = GpVec::zero().set_linear_form_2(
            self.ray,
            &vref,
            1.0,
            &GpVec::from_pnts(&self.ptrst, &self.pts),
        );

        // Derivative corresponding to u1 (`cxx:202-212`).
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

        // Derivative corresponding to v1 (`cxx:214-225`).
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

    /// `Values(X, F, D)` (`cxx:235-296`): the same expressions as
    /// [`Self::compute_value`] / [`Self::compute_derivatives`], inlined the way
    /// OCCT inlines them. Unlike `SurfRstConstRad`, `F(3)` is the plain
    /// `vref.SquareMagnitude() - ray * ray` here too (`cxx:265`), so `Value`
    /// and `Values` agree numerically.
    pub(super) fn compute_values(
        &mut self,
        x: &MathVector,
        f: &mut MathVector,
        d: &mut MathMatrix,
    ) -> bool {
        let (pts, d1u1, d1v1, d2u1, d2v1, d2uv1) = self.surf.d2(x.value(1), x.value(2));
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

        f.set_value(3, vref.square_magnitude() - self.ray * self.ray);

        // Derivative corresponding to u1 (`cxx:267-277`).
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

        // Derivative corresponding to v1 (`cxx:279-290`).
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
