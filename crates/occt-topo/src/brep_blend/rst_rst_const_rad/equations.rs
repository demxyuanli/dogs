//! The two equations of `BRepBlend_RstRstConstRad` and their derivatives,
//! plus the shared circle-centre construction.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_RstRstConstRad.cxx:39-49`
//! (`t3dto2d`), `:95-135` (`Value` / `Derivatives` / `Values`), `:508-543`
//! (`CenterCircleRst1Rst2`) (OCCT 8.0.0).
//!
//! The vector `X` is `(u, v)`: `u` parameterises `rst1` on `surf1`, `v`
//! parameterises `rst2` on `surf2`. Both equations say the two restriction
//! points lie in the section plane: `nplan . P + theD`.

use occt_core::gp::{GpPnt, GpVec};
use occt_geom::curve::Curve;
use occt_math::{MathMatrix, MathVector};

use super::SurfRstRstConstRad;

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

impl SurfRstRstConstRad {
    /// `Value(X, F)` (`cxx:95-107`).
    pub(super) fn compute_value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        self.ptrst1 = self.cons1.d0(x.value(1));
        self.ptrst2 = self.cons2.d0(x.value(2));

        f.set_value(1, self.nplan.dot(&pnt_vec(&self.ptrst1)) + self.the_d);
        f.set_value(2, self.nplan.dot(&pnt_vec(&self.ptrst2)) + self.the_d);
        true
    }

    /// `Derivatives(X, D)` (`cxx:111-127`).
    pub(super) fn compute_derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        let (ptrst1, d11) = self.cons1.d1(x.value(1));
        self.ptrst1 = ptrst1;
        let (ptrst2, d21) = self.cons2.d1(x.value(2));
        self.ptrst2 = ptrst2;

        d.set_value(1, 1, self.nplan.dot(&d11));
        d.set_value(1, 2, 0.0);

        d.set_value(2, 1, 0.0);
        d.set_value(2, 2, self.nplan.dot(&d21));
        true
    }

    /// `Values(X, F, D)` (`cxx:131-135`): `Value` then `Derivatives`.
    pub(super) fn compute_values(
        &mut self,
        x: &MathVector,
        f: &mut MathVector,
        d: &mut MathMatrix,
    ) -> bool {
        self.compute_value(x, f);
        self.compute_derivatives(x, d);
        true
    }

    /// `CenterCircleRst1Rst2(PtRst1, PtRst2, np, Center, VdMed)` (`cxx:508-543`):
    /// the centre of the circle of radius `|ray|` through the two restriction
    /// points, in the plane of normal `np`.
    ///
    /// `VdMed` is the un-normalised direction of the perpendicular bisector; the
    /// caller `Section` divides it by its magnitude to get the derivative of the
    /// centre, so the port returns it as OCCT leaves it.
    ///
    /// Returns `None` where OCCT returns `false` (`Dist < -1.E-07`, `cxx:522-525`):
    /// the two points are further apart than the circle diameter.
    pub(super) fn center_circle_rst1_rst2(
        &self,
        pt_rst1: &GpPnt,
        pt_rst2: &GpPnt,
        np: &GpVec,
    ) -> (Option<GpPnt>, GpVec) {
        let rst1rst2 = GpVec::from_pnts(pt_rst1, pt_rst2);
        let mut vd_med = rst1rst2.crossed(np);
        let norm2 = rst1rst2.square_magnitude();
        let dist0 = self.ray * self.ray - 0.25 * norm2;

        if self.choix > 2 {
            vd_med.reverse();
        }

        if dist0 < -super::DIST_EPS {
            return (None, vd_med);
        }

        let half = rst1rst2.multiplied_scalar(0.5);
        let base = GpPnt::from_xyz(&pt_rst1.coord.added(&half.coord));
        if dist0 > super::DIST_EPS {
            let dist = dist0.sqrt();
            let vdmed_nor = vd_med.normalized();
            let center =
                GpPnt::from_xyz(&base.coord.added(&vdmed_nor.coord.multiply(dist)));
            (Some(center), vd_med)
        } else {
            (Some(base), vd_med)
        }
    }
}
