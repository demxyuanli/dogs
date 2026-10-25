//! `IsSolution` / `Decroch` of `BRepBlend_RstRstConstRad`.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_RstRstConstRad.cxx:187-297`
//! (`IsSolution`), `:405-483` (`Decroch`) (OCCT 8.0.0).

use occt_core::gp::{GpPnt, GpVec, GpVec2d};
use occt_geom::curve::Curve;
use occt_math::{svd_padded, svd_solve_occt, Gauss, MathMatrix, MathVector};

use crate::blend::BlendDecrochStatus;

use super::equations::pnt_vec;
use super::SurfRstRstConstRad;

impl SurfRstRstConstRad {
    /// `IsSolution(Sol, Tol)` (`cxx:187-297`): both values are within
    /// tolerance, then the two parameters are solved for from the jacobian
    /// (`math_Gauss`, falling back to `math_SVD`), the circle centre is built,
    /// and `maxang` / `minang` / `distmin` are updated from the section angle.
    pub(super) fn check_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        let mut valsol = MathVector::new(1, 2);
        let mut gradsol = MathMatrix::new(1, 2, 1, 2);
        self.compute_values(sol, &mut valsol, &mut gradsol);

        if !(valsol.value(1).abs() <= tol && valsol.value(2).abs() <= tol) {
            self.istangent = true;
            return false;
        }

        // Calculation of tangents (`cxx:202-253`).
        self.prmrst1 = sol.value(1);
        self.pt2drst1 = self.rst1.d0(self.prmrst1);
        self.prmrst2 = sol.value(2);
        self.pt2drst2 = self.rst2.d0(self.prmrst2);

        let (ptrst1, d11) = self.cons1.d1(sol.value(1));
        self.ptrst1 = ptrst1;
        let (ptrst2, d21) = self.cons2.d1(sol.value(2));
        self.ptrst2 = ptrst2;

        let dnplan = GpVec::zero().set_linear_form_2(
            1.0 / self.normtg,
            &self.d2gui,
            -1.0 / self.normtg * self.nplan.dot(&self.d2gui),
            &self.nplan,
        );

        let mut secmember = MathVector::new(1, 2);
        let temp = pnt_vec(&self.ptrst1).subtracted(&pnt_vec(&self.ptgui));
        secmember.set_value(1, self.normtg - dnplan.dot(&temp));
        let temp = pnt_vec(&self.ptrst2).subtracted(&pnt_vec(&self.ptgui));
        secmember.set_value(2, self.normtg - dnplan.dot(&temp));

        let resol = Gauss::new(&gradsol, super::GAUSS_DEFAULT_EPS);
        if resol.is_done() {
            if let Some(x) = resol.solve(&secmember) {
                secmember = x;
            }
            self.istangent = false;
        } else {
            let sing_rs = svd_padded(&gradsol);
            if sing_rs.is_done {
                let dedt = secmember.clone();
                svd_solve_occt(&sing_rs, &dedt, &mut secmember, super::SVD_EPS);
                self.istangent = false;
            } else {
                self.istangent = true;
            }
        }

        if !self.istangent {
            self.tgrst1 = d11.multiplied_scalar(secmember.value(1));
            self.tgrst2 = d21.multiplied_scalar(secmember.value(2));

            let (_, d1urst1, d1vrst1) = self.surf1.d1(self.pt2drst1.x(), self.pt2drst1.y());
            let mut a = 0.0;
            let mut b = 0.0;
            super::equations::t3dto2d(&mut a, &mut b, &self.tgrst1, &d1urst1, &d1vrst1);
            self.tg2drst1 = GpVec2d::new(a, b);
            let (_, d1urst2, d1vrst2) = self.surf2.d1(self.pt2drst2.x(), self.pt2drst2.y());
            // OCCT passes `tgrst1` here, not `tgrst2` (`cxx:253`); kept as is.
            super::equations::t3dto2d(&mut a, &mut b, &self.tgrst1, &d1urst2, &d1vrst2);
            self.tg2drst2 = GpVec2d::new(a, b);
        }

        // Circle centre (`cxx:255-264`).
        let (center, _not_used) = self.center_circle_rst1_rst2(&self.ptrst1, &self.ptrst2, &self.nplan);
        let center = match center {
            Some(c) => c,
            None => return false,
        };

        let n1 = GpVec::from_pnts(&center, &self.ptrst1).normalized();
        let n2 = GpVec::from_pnts(&center, &self.ptrst2).normalized();

        let cosa = n1.dot(&n2);
        let mut sina = self.nplan.dot(&n1.crossed(&n2));
        if self.choix % 2 != 0 {
            // `// nplan is changed into -nplan` (`cxx:276`).
            sina = -sina;
        }

        let mut angle = cosa.acos();
        if sina < 0.0 {
            angle = 2.0 * std::f64::consts::PI - angle;
        }

        if angle > self.maxang {
            self.maxang = angle;
        }
        if angle < self.minang {
            self.minang = angle;
        }
        self.distmin = self.distmin.min(self.ptrst1.distance(&self.ptrst2));

        true
    }

    /// `Decroch(Sol, NRst1, TgRst1, NRst2, TgRst2)` (`cxx:405-483`): each
    /// restriction leaves when the in-plane restriction normal stops pointing
    /// against the section tangent (`dot > -1e-10`).
    ///
    /// `surfref1` / `rstref1` / `surfref2` / `rstref2` are the handles
    /// `Set(SurfRef1, RstRef1, SurfRef2, RstRef2)` stored. OCCT dereferences
    /// them without a null check (`cxx:415-422`); the port's `Option`s are
    /// unwrapped, so a missing `Set` panics exactly where OCCT would fault.
    pub(super) fn compute_decroch(
        &self,
        sol: &MathVector,
        n_rst1: &mut GpVec,
        tg_rst1: &mut GpVec,
        n_rst2: &mut GpVec,
        tg_rst2: &mut GpVec,
    ) -> BlendDecrochStatus {
        let rstref1 = self
            .rstref1
            .as_ref()
            .expect("BRepBlend_RstRstConstRad::Decroch: Set(SurfRef1, RstRef1, ..) not called");
        let surfref1 = self
            .surfref1
            .as_ref()
            .expect("BRepBlend_RstRstConstRad::Decroch: Set(SurfRef1, RstRef1, ..) not called");
        let rstref2 = self
            .rstref2
            .as_ref()
            .expect("BRepBlend_RstRstConstRad::Decroch: Set(.., SurfRef2, RstRef2) not called");
        let surfref2 = self
            .surfref2
            .as_ref()
            .expect("BRepBlend_RstRstConstRad::Decroch: Set(.., SurfRef2, RstRef2) not called");

        let uv = rstref1.d0(sol.value(1));
        let (pt_tmp1, d1u, d1v) = surfref1.d1(uv.x(), uv.y());
        *n_rst1 = d1u.crossed(&d1v);

        let uv = rstref2.d0(sol.value(2));
        let (pt_tmp2, d1u, d1v) = surfref2.d1(uv.x(), uv.y());
        *n_rst2 = d1u.crossed(&d1v);

        let (center, _not_used) = self.center_circle_rst1_rst2(&pt_tmp1, &pt_tmp2, &self.nplan);
        // OCCT ignores the failure here (`cxx:424`); `Center` keeps whatever
        // `CenterCircleRst1Rst2` left, which in the failure arm is untouched.
        let center = center.unwrap_or_else(GpPnt::zero);

        let norm = self.nplan.crossed(n_rst1).magnitude();
        let unsurnorm = 1.0 / norm;
        let mut n_rst1_in_plane = GpVec::zero().set_linear_form_2(
            self.nplan.dot(n_rst1) * unsurnorm,
            &self.nplan,
            -unsurnorm,
            n_rst1,
        );

        let centptrst = GpVec::from_pnts(&center, &pt_tmp1);
        if centptrst.dot(&n_rst1_in_plane) < 0.0 {
            n_rst1_in_plane.reverse();
        }
        *tg_rst1 = self.nplan.crossed(&centptrst);

        let norm = self.nplan.crossed(n_rst2).magnitude();
        let unsurnorm = 1.0 / norm;
        let mut n_rst2_in_plane = GpVec::zero().set_linear_form_2(
            self.nplan.dot(n_rst2) * unsurnorm,
            &self.nplan,
            -unsurnorm,
            n_rst2,
        );

        let centptrst = GpVec::from_pnts(&center, &pt_tmp2);
        if centptrst.dot(&n_rst2_in_plane) < 0.0 {
            n_rst2_in_plane.reverse();
        }
        *tg_rst2 = self.nplan.crossed(&centptrst);

        if self.choix % 2 != 0 {
            tg_rst1.reverse();
            tg_rst2.reverse();
        }

        // The vectors are returned (`cxx:459-481`).
        if n_rst1_in_plane.dot(tg_rst1) > -super::DECROCH_DOT_EPS {
            if n_rst2_in_plane.dot(tg_rst2) < super::DECROCH_DOT_EPS {
                BlendDecrochStatus::DecrochBoth
            } else {
                BlendDecrochStatus::DecrochRst1
            }
        } else if n_rst2_in_plane.dot(tg_rst2) < super::DECROCH_DOT_EPS {
            BlendDecrochStatus::DecrochRst2
        } else {
            BlendDecrochStatus::NoDecroch
        }
    }
}
