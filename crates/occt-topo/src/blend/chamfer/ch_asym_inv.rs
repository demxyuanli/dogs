//! `BlendFunc_ChAsymInv` — the inverse of [`ChAsym`]: the same four equations
//! with the solution restricted to a curve on one of the two supports.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_ChAsymInv.hxx:24-84`,
//! `BlendFunc_ChAsymInv.cxx:24-404` (OCCT 8.0.0).
//!
//! The vector `X` is `(t, w, U, V)`: the parameter on the curve on surface, the
//! parameter on the guide line, and the parameters of the partner point on the
//! other surface.
//!
//! Translation notes:
//! - `csurf` is an `occ::handle<Adaptor2d_Curve2d>`; the port stores it as
//!   `Option<Arc<dyn Curve2d>>` because `Set(OnFirst, C)` can be handed the
//!   empty handle. OCCT dereferences it unguarded in `Value` / `Derivatives` /
//!   `Values` / `IsSolution`; the port answers `false` instead of faulting, and
//!   `GetTolerance` / `GetBounds` use `0.` for the missing curve, as
//!   [`super::gen_chamf_inv`] does.
//! - `angle` (`hxx:76`) is written by `Set(Dist1, Angle, Choix)` and read by no
//!   other body of the `.cxx`; it stays as a dead member.
//! - The four-argument `Values` ignores the `Value` answer in OCCT and returns
//!   the `ComputeValues` one; the port keeps that.

use std::sync::Arc;

use occt_core::precision::Precision;
use occt_geom::approx_same_parameter;
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_geom2d::curve::Curve2d;
use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use crate::blend::BlendFuncInv;

/// `class BlendFunc_ChAsymInv : public Blend_FuncInv` (`hxx:24-84`).
pub struct ChAsymInv {
    /// `surf1` (`hxx:72`).
    surf1: Arc<dyn Surface>,
    /// `surf2` (`hxx:73`).
    surf2: Arc<dyn Surface>,
    /// `dist1` (`hxx:74`): the distance on the first support.
    dist1: f64,
    /// `angle` (`hxx:76`): stored, never read (`cxx:43-49`).
    #[allow(dead_code)]
    angle: f64,
    /// `tgang` (`hxx:77`): `tan(Angle)`.
    tgang: f64,
    /// `curv` (`hxx:78`): the guide line.
    curv: Arc<dyn Curve>,
    /// `csurf` (`hxx:79`): the curve on surface the solution lies on.
    csurf: Option<Arc<dyn Curve2d>>,
    /// `choix` (`hxx:80`): the "quadrant" of the chamfer.
    choix: i32,
    /// `first` (`hxx:81`): `true` when `csurf` lies on `surf1`.
    first: bool,
    /// `FX` (`hxx:82`): the values of the four equations.
    fx: MathVector,
    /// `DX` (`hxx:83`): the jacobian of the four equations.
    dx: MathMatrix,
}

impl ChAsymInv {
    /// `BlendFunc_ChAsymInv(S1, S2, C)` (`cxx:26-40`).
    pub fn new(surf1: Arc<dyn Surface>, surf2: Arc<dyn Surface>, curv: Arc<dyn Curve>) -> Self {
        Self {
            surf1,
            surf2,
            dist1: f64::MAX,
            angle: f64::MAX,
            tgang: f64::MAX,
            curv,
            csurf: None,
            choix: 0,
            first: false,
            fx: MathVector::new(1, 4),
            dx: MathMatrix::new(1, 4, 1, 4),
        }
    }

    /// `Set(Dist1, Angle, Choix)` (`cxx:43-49`).
    pub fn set_dist_angle(&mut self, dist1: f64, angle: f64, choix: i32) {
        self.dist1 = dist1.abs();
        self.angle = angle;
        self.tgang = angle.tan();
        self.choix = choix;
    }

    /// `IsSolution(Sol, Tol)` (`cxx:136-173`): the four residuals of
    /// `Value(Sol, ...)` against `Tol`, scaled by the section-plane geometry
    /// computed from `nplan` with no `choix` reversal.
    pub(super) fn check_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        let Some(csurf) = self.csurf.clone() else {
            return false;
        };

        let mut valsol = MathVector::new(1, 4);

        let (_ptgui, d1gui) = self.curv.d1(sol.value(2));
        let nplan = d1gui.normalized();

        let pt2d = csurf.value(sol.value(1));

        let (pts1, d1u1, d1v1, pts2) = if self.first {
            let (p1, u1, v1) = self.surf1.d1(pt2d.x(), pt2d.y());
            (p1, u1, v1, self.surf2.value(sol.value(3), sol.value(4)))
        } else {
            let (p1, u1, v1) = self.surf1.d1(sol.value(3), sol.value(4));
            (p1, u1, v1, self.surf2.value(pt2d.x(), pt2d.y()))
        };

        let nsurf1 = d1u1.crossed(&d1v1);
        let tsurf1 = nsurf1.crossed(&nplan);

        let s1s2 = occt_core::gp::GpVec::from_pnts(&pts1, &pts2);
        let psca_inv = 1.0 / tsurf1.dot(&s1s2);
        let nordu1 = d1u1.magnitude();
        let nordv1 = d1v1.magnitude();

        let temp = 2.0 * (nordu1 + nordv1) * s1s2.magnitude() + 2.0 * nordu1 * nordv1;

        self.value(sol, &mut valsol);

        valsol.value(1).abs() < tol
            && valsol.value(2).abs() < tol
            && valsol.value(3).abs() < 2.0 * self.dist1 * tol
            && valsol.value(4).abs() < tol * (1.0 + self.tgang) * psca_inv.abs() * temp
    }

    /// The `FX` block of `ComputeValues` (`cxx:248-256`).
    fn fill_fx(
        &mut self,
        nplan: &occt_core::gp::GpVec,
        ptgui: &occt_core::gp::GpPnt,
        pts1: &occt_core::gp::GpPnt,
        pts2: &occt_core::gp::GpPnt,
        nps1: &occt_core::gp::GpVec,
        f4: f64,
    ) {
        let dist = ptgui.xyz().dot(nplan.xyz());
        self.fx.set_value(1, pts1.xyz().dot(nplan.xyz()) - dist);
        self.fx.set_value(2, pts2.xyz().dot(nplan.xyz()) - dist);
        self.fx
            .set_value(3, self.dist1 * self.dist1 - nps1.square_magnitude());
        self.fx.set_value(4, self.tgang - f4);
    }

    /// `ComputeValues(X, DegF, DegL)` (`cxx:177-345`): fills `FX` for
    /// `DegF <= 0` and `DX` for `DegL >= 1`.
    pub fn compute_values(&mut self, x: &MathVector, deg_f: i32, deg_l: i32) -> bool {
        if deg_f > deg_l {
            return false;
        }
        let Some(csurf) = self.csurf.clone() else {
            return false;
        };

        if deg_f == 0 && deg_l == 0 {
            // `cxx:191-212`: first order scene, only `FX` is wanted.
            let (ptgui, d1gui) = self.curv.d1(x.value(2));
            let mut nplan = d1gui.normalized();

            if self.choix % 2 != 0 {
                nplan.reverse();
            }
            let pt2d = csurf.value(x.value(1));

            let (pts1, d1u1, d1v1, pts2) = if self.first {
                let (p1, u1, v1) = self.surf1.d1(pt2d.x(), pt2d.y());
                (p1, u1, v1, self.surf2.value(x.value(3), x.value(4)))
            } else {
                let (p1, u1, v1) = self.surf1.d1(x.value(3), x.value(4));
                (p1, u1, v1, self.surf2.value(pt2d.x(), pt2d.y()))
            };

            let nps1 = occt_core::gp::GpVec::from_pnts(&ptgui, &pts1);
            let s1s2 = occt_core::gp::GpVec::from_pnts(&pts1, &pts2);
            let nsurf1 = d1u1.crossed(&d1v1);
            let tsurf1 = nsurf1.crossed(&nplan);
            let psca_inv = 1.0 / s1s2.dot(&tsurf1);
            let f4 = nplan.dot(&tsurf1.crossed(&s1s2)) * psca_inv;

            self.fill_fx(&nplan, &ptgui, &pts1, &pts2, &nps1, f4);
            return true;
        }

        // `cxx:213-239`: second order scene, used by `FX` and/or `DX`.
        let (ptgui, d1gui, d2gui) = self.curv.d2(x.value(2));
        let mut nplan = d1gui.normalized();
        let mut normg = d1gui.magnitude();
        let mut dnplan = d2gui
            .subtracted(&nplan.multiplied_scalar(nplan.dot(&d2gui)))
            .divided(normg);

        if self.choix % 2 != 0 {
            nplan.reverse();
            dnplan.reverse();
            normg = -normg;
        }

        let (pt2d, v2d) = csurf.d1(x.value(1));

        let (pts1, d1u1, d1v1, d2u1, d2v1, d2uv1, pts2, d1u2, d1v2) = if self.first {
            let (p1, u1, v1, uu1, vv1, uv1) = self.surf1.d2(pt2d.x(), pt2d.y());
            let (p2, u2, v2) = self.surf2.d1(x.value(3), x.value(4));
            (p1, u1, v1, uu1, vv1, uv1, p2, u2, v2)
        } else {
            let (p1, u1, v1, uu1, vv1, uv1) = self.surf1.d2(x.value(3), x.value(4));
            let (p2, u2, v2) = self.surf2.d1(pt2d.x(), pt2d.y());
            (p1, u1, v1, uu1, vv1, uv1, p2, u2, v2)
        };

        let nps1 = occt_core::gp::GpVec::from_pnts(&ptgui, &pts1);
        let nps2 = occt_core::gp::GpVec::from_pnts(&ptgui, &pts2);
        let s1s2 = occt_core::gp::GpVec::from_pnts(&pts1, &pts2);
        let nsurf1 = d1u1.crossed(&d1v1);
        let tsurf1 = nsurf1.crossed(&nplan);
        let psca_inv = 1.0 / s1s2.dot(&tsurf1);
        let f4 = nplan.dot(&tsurf1.crossed(&s1s2)) * psca_inv;

        if deg_f == 0 {
            self.fill_fx(&nplan, &ptgui, &pts1, &pts2, &nps1, f4);
        }

        if deg_l == 1 {
            let dwtsurf1 = nsurf1.crossed(&dnplan);

            if self.first {
                let dw1pts1 = d1u1
                    .multiplied_scalar(v2d.x())
                    .added(&d1v1.multiplied_scalar(v2d.y()));
                let dw1du1 = d2u1
                    .multiplied_scalar(v2d.x())
                    .added(&d2uv1.multiplied_scalar(v2d.y()));
                let dw1dv1 = d2uv1
                    .multiplied_scalar(v2d.x())
                    .added(&d2v1.multiplied_scalar(v2d.y()));
                let dw1csurf = dw1du1
                    .crossed(&d1v1)
                    .added(&d1u1.crossed(&dw1dv1))
                    .crossed(&nplan);

                self.dx.set_value(1, 1, nplan.dot(&dw1pts1));
                self.dx.set_value(1, 2, dnplan.dot(&nps1) - normg);
                self.dx.set_value(1, 3, 0.0);
                self.dx.set_value(1, 4, 0.0);

                self.dx.set_value(2, 1, 0.0);
                self.dx.set_value(2, 2, dnplan.dot(&nps2) - normg);
                self.dx.set_value(2, 3, nplan.dot(&d1u2));
                self.dx.set_value(2, 4, nplan.dot(&d1v2));

                let temp_vec = nps1.multiplied_scalar(2.0);
                self.dx.set_value(3, 1, -dw1pts1.dot(&temp_vec));
                self.dx.set_value(3, 2, d1gui.dot(&temp_vec));
                self.dx.set_value(3, 3, 0.0);
                self.dx.set_value(3, 4, 0.0);

                let mut temp = f4 * (dw1csurf.dot(&s1s2) - tsurf1.dot(&dw1pts1));
                temp += nplan.dot(&tsurf1.crossed(&dw1pts1).subtracted(&dw1csurf.crossed(&s1s2)));
                self.dx.set_value(4, 1, psca_inv * temp);

                let mut temp = f4 * dwtsurf1.dot(&s1s2);
                temp -= dnplan.dot(&temp_vec) + nplan.dot(&dwtsurf1.crossed(&s1s2));
                self.dx.set_value(4, 2, psca_inv * temp);

                let temp = f4 * tsurf1.dot(&d1u2) - nplan.dot(&tsurf1.crossed(&d1u2));
                self.dx.set_value(4, 3, psca_inv * temp);

                let temp = f4 * tsurf1.dot(&d1v2) - nplan.dot(&tsurf1.crossed(&d1v2));
                self.dx.set_value(4, 4, psca_inv * temp);
            } else {
                let d1utsurf1 = d2u1
                    .crossed(&d1v1)
                    .added(&d1u1.crossed(&d2uv1))
                    .crossed(&nplan);
                let d1vtsurf1 = d2uv1
                    .crossed(&d1v1)
                    .added(&d1u1.crossed(&d2v1))
                    .crossed(&nplan);
                let dw2pts2 = d1u2
                    .multiplied_scalar(v2d.x())
                    .added(&d1v2.multiplied_scalar(v2d.y()));

                self.dx.set_value(1, 1, 0.0);
                self.dx.set_value(1, 2, dnplan.dot(&nps1) - normg);
                self.dx.set_value(1, 3, nplan.dot(&d1u1));
                self.dx.set_value(1, 4, nplan.dot(&d1v1));

                self.dx.set_value(2, 1, nplan.dot(&dw2pts2));
                self.dx.set_value(2, 2, dnplan.dot(&nps2) - normg);
                self.dx.set_value(2, 3, 0.0);
                self.dx.set_value(2, 4, 0.0);

                let temp_vec = nps1.multiplied_scalar(2.0);
                self.dx.set_value(3, 1, 0.0);
                self.dx.set_value(3, 2, d1gui.dot(&temp_vec));

                let temp_vec = temp_vec.reversed();
                self.dx.set_value(3, 3, d1u1.dot(&temp_vec));
                self.dx.set_value(3, 4, d1v1.dot(&temp_vec));

                let temp = f4 * tsurf1.dot(&dw2pts2) - nplan.dot(&tsurf1.crossed(&dw2pts2));
                self.dx.set_value(4, 1, psca_inv * temp);

                let mut temp = f4 * dwtsurf1.dot(&s1s2);
                temp -= dnplan.dot(&temp_vec) + nplan.dot(&dwtsurf1.crossed(&s1s2));
                self.dx.set_value(4, 2, psca_inv * temp);

                let mut temp = f4 * (d1utsurf1.dot(&s1s2) - tsurf1.dot(&d1u1));
                temp += nplan.dot(&tsurf1.crossed(&d1u1).subtracted(&d1utsurf1.crossed(&s1s2)));
                self.dx.set_value(4, 3, psca_inv * temp);

                let mut temp = f4 * (d1vtsurf1.dot(&s1s2) - tsurf1.dot(&d1v1));
                temp += nplan.dot(&tsurf1.crossed(&d1v1).subtracted(&d1vtsurf1.crossed(&s1s2)));
                self.dx.set_value(4, 4, psca_inv * temp);
            }
        }

        true
    }
}

impl MathFunctionSetWithDerivatives for ChAsymInv {
    /// `NbVariables()` = `4`, not overridden.
    fn nb_variables(&self) -> usize {
        4
    }

    /// `NbEquations()` (`cxx:53-56`).
    fn nb_equations(&self) -> usize {
        4
    }

    /// `Value(X, F)` (`cxx:349-354`): `ComputeValues(X, 0, 0)`.
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        let error = self.compute_values(x, 0, 0);
        *f = self.fx.clone();
        error
    }

    /// `Derivatives(X, D)` (`cxx:358-363`): `ComputeValues(X, 1, 1)`.
    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        let error = self.compute_values(x, 1, 1);
        *d = self.dx.clone();
        error
    }

    /// `Values(X, F, D)` (`cxx:367-404`): `ComputeValues(X, 0, 1)`.
    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        let error = self.compute_values(x, 0, 1);
        *f = self.fx.clone();
        *d = self.dx.clone();
        error
    }
}

impl BlendFuncInv for ChAsymInv {
    /// `Set(OnFirst, COnSurf)` (`cxx:60-64`).
    fn set(&mut self, on_first: bool, c_on_surf: Option<Arc<dyn Curve2d>>) {
        self.first = on_first;
        self.csurf = c_on_surf;
    }

    /// `GetTolerance(Tolerance, Tol)` (`cxx:67-83`).
    fn get_tolerance(&self, tolerance: &mut MathVector, tol: f64) {
        let res = match &self.csurf {
            Some(c) => c.resolution(tol),
            None => 0.0,
        };
        tolerance.set_value(1, res);
        tolerance.set_value(2, self.curv.resolution(tol));
        let partner = if self.first { &self.surf2 } else { &self.surf1 };
        tolerance.set_value(3, approx_same_parameter::u_resolution(partner.as_ref(), tol));
        tolerance.set_value(4, approx_same_parameter::v_resolution(partner.as_ref(), tol));
    }

    /// `GetBounds(InfBound, SupBound)` (`cxx:86-132`).
    fn get_bounds(&self, inf_bound: &mut MathVector, sup_bound: &mut MathVector) {
        let (cfirst, clast) = match &self.csurf {
            Some(c) => (c.first_parameter(), c.last_parameter()),
            None => (0.0, 0.0),
        };
        inf_bound.set_value(1, cfirst);
        inf_bound.set_value(2, self.curv.first_parameter());
        sup_bound.set_value(1, clast);
        sup_bound.set_value(2, self.curv.last_parameter());

        let partner = if self.first { &self.surf2 } else { &self.surf1 };
        let (umin, umax) = partner.u_range();
        let (vmin, vmax) = partner.v_range();
        inf_bound.set_value(3, umin);
        inf_bound.set_value(4, vmin);
        sup_bound.set_value(3, umax);
        sup_bound.set_value(4, vmax);

        for i in 3..=4 {
            if !Precision::is_infinite(inf_bound.value(i))
                && !Precision::is_infinite(sup_bound.value(i))
            {
                let range = sup_bound.value(i) - inf_bound.value(i);
                inf_bound.set_value(i, inf_bound.value(i) - range);
                sup_bound.set_value(i, sup_bound.value(i) + range);
            }
        }
    }

    /// `IsSolution(Sol, Tol)` (`cxx:136-173`).
    fn is_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        self.check_solution(sol, tol)
    }
}
