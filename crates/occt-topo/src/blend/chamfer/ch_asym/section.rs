//! The two `Section` bodies of [`ChAsym`].
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_ChAsym.cxx:580-722`
//! (OCCT 8.0.0). The second overload repeats the tangent solve of `IsSolution`
//! (`cxx:632-709`) verbatim, including the `math_Gauss` (`1.e-9`) then
//! `math_SVD` (`1.e-6`) fallback; it is kept duplicated rather than shared so
//! each body stays a direct transcription of its `.cxx` range.

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_math::{Gauss, MathFunctionSetWithDerivatives, MathMatrix, MathVector, SVD};

use crate::blend::{BlendAppFunction, BlendFunction, BlendPoint};

use super::ChAsym;

impl ChAsym {
    /// `Section(P, Poles, Poles2d, Weights)` (`cxx:580-607`): solve the four
    /// equations at `P.Parameter()` and take the two points it updated.
    pub(super) fn section_poles_impl(
        &mut self,
        p: &BlendPoint,
        poles: &mut [GpPnt],
        poles2d: &mut [GpPnt2d],
        weights: &mut [f64],
    ) {
        let prm = p.parameter();
        let low = 0;
        let upp = poles.len() - 1;
        let mut x = MathVector::new(1, 4);
        let mut f = MathVector::new(1, 4);

        let (u1, v1) = p
            .parameters_on_s1()
            .expect("BlendFunc_ChAsym::Section: parameters on S1 undefined");
        let (u2, v2) = p
            .parameters_on_s2()
            .expect("BlendFunc_ChAsym::Section: parameters on S2 undefined");
        x.set_value(1, u1);
        x.set_value(2, v1);
        x.set_value(3, u2);
        x.set_value(4, v2);
        let last2d = poles2d.len() - 1;
        poles2d[0] = GpPnt2d::new(u1, v1);
        poles2d[last2d] = GpPnt2d::new(u2, v2);

        self.set(prm);
        self.value(&x, &mut f);
        poles[low] = *self.point_on_s1();
        poles[upp] = *self.point_on_s2();
        weights[low] = 1.0;
        weights[upp] = 1.0;
    }

    /// `Section(P, Poles, DPoles, Poles2d, DPoles2d, Weights, DWeights)`
    /// (`cxx:609-722`): the first/last section with the two 3d tangents, the
    /// two 2d tangents and the zero weight derivatives.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn section_first_last(
        &mut self,
        p: &BlendPoint,
        poles: &mut [GpPnt],
        dpoles: &mut [GpVec],
        poles2d: &mut [GpPnt2d],
        dpoles2d: &mut [GpVec2d],
        weights: &mut [f64],
        dweights: &mut [f64],
    ) -> bool {
        let mut sol = MathVector::new(1, 4);
        let mut valsol = MathVector::new(1, 4);
        let mut secmember = MathVector::new(1, 4);
        let mut gradsol = MathMatrix::new(1, 4, 1, 4);
        let prm = p.parameter();
        let low = 0;
        let upp = poles.len() - 1;

        let (u1, v1) = p
            .parameters_on_s1()
            .expect("BlendFunc_ChAsym::Section: parameters on S1 undefined");
        let (u2, v2) = p
            .parameters_on_s2()
            .expect("BlendFunc_ChAsym::Section: parameters on S2 undefined");
        sol.set_value(1, u1);
        sol.set_value(2, v1);
        sol.set_value(3, u2);
        sol.set_value(4, v2);
        self.set(prm);

        let last2d = poles2d.len() - 1;
        poles2d[0] = GpPnt2d::new(sol.value(1), sol.value(2));
        poles2d[last2d] = GpPnt2d::new(sol.value(3), sol.value(4));
        poles[low] = *self.point_on_s1();
        poles[upp] = *self.point_on_s2();
        weights[low] = 1.0;
        weights[upp] = 1.0;

        let (ptgui, d1gui, d2gui) = self.curv.d2(self.param);
        let mut normg = d1gui.magnitude();
        let mut np = d1gui.normalized();
        let mut dnp = d2gui
            .subtracted(&np.multiplied_scalar(np.dot(&d2gui)))
            .divided(normg);

        if self.choix % 2 != 0 {
            np.reverse();
            dnp.reverse();
            normg = -normg;
        }

        let (p1, d1u1, d1v1) = self.surf1.d1(sol.value(1), sol.value(2));
        self.pt1 = p1;
        let nsurf1 = d1u1.crossed(&d1v1);
        self.tsurf1 = nsurf1.crossed(&np);
        let dwtsurf1 = nsurf1.crossed(&dnp);

        let (p2, d1u2, d1v2) = self.surf2.d1(sol.value(3), sol.value(4));
        self.pt2 = p2;

        let pguis1 = GpVec::from_pnts(&ptgui, &self.pt1);
        let pguis2 = GpVec::from_pnts(&ptgui, &self.pt2);
        let s1s2 = GpVec::from_pnts(&self.pt1, &self.pt2);
        let psca_inv = 1.0 / self.tsurf1.dot(&s1s2);
        let maxpiv = 1.0e-9;
        let nordu1 = d1u1.magnitude();
        let nordv1 = d1v1.magnitude();

        // `cxx:665` evaluates this into `temp` and overwrites `temp` at `cxx:675`
        // before any read; the discarded evaluation is kept.
        let _ = 2.0 * (nordu1 + nordv1) * s1s2.magnitude() + 2.0 * nordu1 * nordv1;

        self.values(&sol, &mut valsol, &mut gradsol);

        secmember.set_value(1, normg - dnp.dot(&pguis1));
        secmember.set_value(2, normg - dnp.dot(&pguis2));
        secmember.set_value(3, -2.0 * d1gui.dot(&pguis1));

        let cross_vec = self.tsurf1.crossed(&s1s2);
        let f4 = np.dot(&cross_vec) * psca_inv;
        let mut temp = dnp.dot(&cross_vec) + np.dot(&dwtsurf1.crossed(&s1s2));
        temp -= f4 * dwtsurf1.dot(&s1s2);
        secmember.set_value(4, psca_inv * temp);

        let resol = Gauss::new(&gradsol, maxpiv);

        if resol.is_done() {
            if let Some(x) = resol.solve(&secmember) {
                secmember = x;
            }
            self.istangent = false;
        } else {
            let sing_rs = SVD::new(&gradsol);
            if sing_rs.is_done {
                let dedt = secmember.clone();
                if let Ok(x) = sing_rs.solve_with_eps(&dedt, 1.0e-6) {
                    secmember = x;
                }
                self.istangent = false;
            } else {
                self.istangent = true;
            }
        }

        if !self.istangent {
            self.tg1 = self
                .tg1
                .set_linear_form_2(secmember.value(1), &d1u1, secmember.value(2), &d1v1);
            self.tg2 = self
                .tg2
                .set_linear_form_2(secmember.value(3), &d1u2, secmember.value(4), &d1v2);
            self.tg12d = GpVec2d::new(secmember.value(1), secmember.value(2));
            self.tg22d = GpVec2d::new(secmember.value(3), secmember.value(4));
        }

        let dist = self.pt1.distance(&self.pt2);
        self.distmin = self.distmin.min(dist);

        if !self.istangent {
            let t1 = self.tangent2d_on_s1();
            dpoles2d[0] = GpVec2d::new(t1.x(), t1.y());
            let t2 = self.tangent2d_on_s2();
            dpoles2d[last2d] = GpVec2d::new(t2.x(), t2.y());

            dpoles[low] = *self.tangent_on_s1();
            dpoles[upp] = *self.tangent_on_s2();
            dweights[low] = 0.0;
            dweights[upp] = 0.0;
        }

        !self.istangent
    }
}
