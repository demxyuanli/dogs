//! `BlendFunc_ChAsym` — an asymmetric chamfer between two surfaces: a constant
//! distance on the first support and a constant angle with the second one.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_ChAsym.hxx:34-220`,
//! `BlendFunc_ChAsym.cxx:32-766` (OCCT 8.0.0).
//!
//! The vector `X` is `(U1, V1, U2, V2)`. The four equations are
//! 1. `P1 . np - ptgui . np = 0` — the first point lies in the section plane,
//! 2. `P2 . np - ptgui . np = 0` — so does the second one,
//! 3. `dist1^2 - |P1 - ptgui|^2 = 0` — the first one is at `dist1`,
//! 4. `tgang - F4 = 0` — the section leaves the second support at `angle`.
//!
//! Translation notes:
//! - `tcurv` (`hxx:203`) is `curv->Trim(First, Last, 1.e-12)` (`cxx:73-76`).
//!   `Geom_TrimmedCurve` keeps the basis parameterisation (`D0` / `D1` / `D2`
//!   delegate to the basis at the same `U`), so every `tcurv` evaluation in
//!   this class equals the one on `curv`; the port has no `tcurv` member and
//!   `set_interval` only records nothing. `GetBounds` (`cxx:88-110`) never
//!   consults it either.
//! - `angle` (`hxx:206`) is written by `Set(Dist1, Angle, Choix)` and read by
//!   no other body of the `.cxx`; it stays as a dead member.
//! - The second `Section` overload (`cxx:609-722`) repeats the tangent solve of
//!   `IsSolution` verbatim and is in [`section`].

use std::sync::Arc;

use occt_core::elib::clib;
use occt_core::gp::{GpDir, GpLin, GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::Precision;
use occt_geom::approx_same_parameter;
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_math::{Gauss, MathFunctionSetWithDerivatives, MathMatrix, MathVector, SVD};

use crate::blend::{next_shape, BlendAppFunction, BlendFunction, BlendPoint};

mod section;

/// `class BlendFunc_ChAsym : public Blend_Function` (`hxx:34-220`).
pub struct ChAsym {
    /// `surf1` (`hxx:199`).
    surf1: Arc<dyn Surface>,
    /// `surf2` (`hxx:200`).
    surf2: Arc<dyn Surface>,
    /// `curv` (`hxx:201`): the guide line.
    curv: Arc<dyn Curve>,
    /// `param` (`hxx:204`): the parameter on the guide line of the section
    /// being solved.
    param: f64,
    /// `dist1` (`hxx:205`): the distance on the first support.
    dist1: f64,
    /// `angle` (`hxx:206`): stored, never read (`cxx:762-766`).
    #[allow(dead_code)]
    angle: f64,
    /// `tgang` (`hxx:207`): `tan(Angle)`.
    tgang: f64,
    /// `nplan` (`hxx:208`): the normal of the section plane.
    nplan: GpVec,
    /// `pt1` (`hxx:209`): the point on the first support.
    pt1: GpPnt,
    /// `tsurf1` (`hxx:210`): `Nsurf1 ^ np`.
    tsurf1: GpVec,
    /// `pt2` (`hxx:211`): the point on the second support.
    pt2: GpPnt,
    /// `FX` (`hxx:212`): the values of the four equations.
    fx: MathVector,
    /// `DX` (`hxx:213`): the jacobian of the four equations.
    dx: MathMatrix,
    /// `istangent` (`hxx:214`): `true` when the tangent solve failed.
    istangent: bool,
    /// `tg1` (`hxx:215`).
    tg1: GpVec,
    /// `tg12d` (`hxx:216`).
    tg12d: GpVec2d,
    /// `tg2` (`hxx:217`).
    tg2: GpVec,
    /// `tg22d` (`hxx:218`).
    tg22d: GpVec2d,
    /// `choix` (`hxx:219`): the "quadrant" of the chamfer.
    choix: i32,
    /// `distmin` (`hxx:220`): the smallest distance between the two points of a
    /// section seen so far, seeded to `RealLast()` (`cxx:48`).
    distmin: f64,
}

impl ChAsym {
    /// `BlendFunc_ChAsym(S1, S2, C)` (`cxx:34-50`).
    pub fn new(surf1: Arc<dyn Surface>, surf2: Arc<dyn Surface>, curv: Arc<dyn Curve>) -> Self {
        Self {
            surf1,
            surf2,
            curv,
            param: 0.0,
            dist1: f64::MAX,
            angle: f64::MAX,
            tgang: f64::MAX,
            nplan: GpVec::zero(),
            pt1: GpPnt::zero(),
            tsurf1: GpVec::zero(),
            pt2: GpPnt::zero(),
            fx: MathVector::new(1, 4),
            dx: MathMatrix::new(1, 4, 1, 4),
            istangent: true,
            tg1: GpVec::zero(),
            tg12d: GpVec2d::new(0.0, 0.0),
            tg2: GpVec::zero(),
            tg22d: GpVec2d::new(0.0, 0.0),
            choix: 0,
            distmin: f64::MAX,
        }
    }

    /// `Set(Dist1, Angle, Choix)` (`cxx:762-766`).
    pub fn set_dist_angle(&mut self, dist1: f64, angle: f64, choix: i32) {
        self.dist1 = dist1.abs();
        self.angle = angle;
        self.tgang = angle.tan();
        self.choix = choix;
    }

    /// `ComputeValues(X, DegF, DegL)` (`cxx:217-308`): fills `FX` for
    /// `DegF <= 0` and `DX` for `DegL >= 1` at the current guide parameter.
    pub fn compute_values(&mut self, x: &MathVector, deg_f: i32, deg_l: i32) -> bool {
        if deg_f > deg_l {
            return false;
        }

        let (ptgui, d1gui) = self.curv.d1(self.param);
        self.nplan = d1gui.normalized();
        let mut np = self.nplan;
        if self.choix % 2 != 0 {
            np.reverse();
        }

        let d1u1;
        let d1v1;
        let d2u1;
        let d2v1;
        let d2uv1;
        let d1u2;
        let d1v2;

        if deg_f == 0 && deg_l == 0 {
            let (p1, u1, v1) = self.surf1.d1(x.value(1), x.value(2));
            self.pt1 = p1;
            d1u1 = u1;
            d1v1 = v1;
            self.pt2 = self.surf2.value(x.value(3), x.value(4));
            d2u1 = GpVec::zero();
            d2v1 = GpVec::zero();
            d2uv1 = GpVec::zero();
            d1u2 = GpVec::zero();
            d1v2 = GpVec::zero();
        } else {
            let (p1, u1, v1, uu1, vv1, uv1) =
                self.surf1.d2(x.value(1), x.value(2));
            self.pt1 = p1;
            d1u1 = u1;
            d1v1 = v1;
            d2u1 = uu1;
            d2v1 = vv1;
            d2uv1 = uv1;
            let (p2, u2, v2) = self.surf2.d1(x.value(3), x.value(4));
            self.pt2 = p2;
            d1u2 = u2;
            d1v2 = v2;
        }

        let nsurf1 = d1u1.crossed(&d1v1);
        self.tsurf1 = nsurf1.crossed(&np);

        let nps1 = GpVec::from_pnts(&ptgui, &self.pt1);
        let s1s2 = GpVec::from_pnts(&self.pt1, &self.pt2);
        let psca_inv = 1.0 / self.tsurf1.dot(&s1s2);
        let f4 = np.dot(&self.tsurf1.crossed(&s1s2)) * psca_inv;

        if deg_f == 0 {
            let dist = ptgui.xyz().dot(np.xyz());

            self.fx.set_value(1, self.pt1.xyz().dot(np.xyz()) - dist);
            self.fx.set_value(2, self.pt2.xyz().dot(np.xyz()) - dist);
            self.fx.set_value(3, self.dist1 * self.dist1 - nps1.square_magnitude());
            self.fx.set_value(4, self.tgang - f4);
        }

        if deg_l == 1 {
            let d1utsurf1 = d2u1
                .crossed(&d1v1)
                .added(&d1u1.crossed(&d2uv1))
                .crossed(&np);
            let d1vtsurf1 = d2uv1
                .crossed(&d1v1)
                .added(&d1u1.crossed(&d2v1))
                .crossed(&np);

            self.dx.set_value(1, 1, np.dot(&d1u1));
            self.dx.set_value(1, 2, np.dot(&d1v1));
            self.dx.set_value(1, 3, 0.0);
            self.dx.set_value(1, 4, 0.0);

            self.dx.set_value(2, 1, 0.0);
            self.dx.set_value(2, 2, 0.0);
            self.dx.set_value(2, 3, np.dot(&d1u2));
            self.dx.set_value(2, 4, np.dot(&d1v2));

            let temp_vec = nps1.multiplied_scalar(-2.0);
            self.dx.set_value(3, 1, d1u1.dot(&temp_vec));
            self.dx.set_value(3, 2, d1v1.dot(&temp_vec));
            self.dx.set_value(3, 3, 0.0);
            self.dx.set_value(3, 4, 0.0);

            let mut temp = f4 * (d1utsurf1.dot(&s1s2) - self.tsurf1.dot(&d1u1));
            temp += np.dot(
                &self
                    .tsurf1
                    .crossed(&d1u1)
                    .subtracted(&d1utsurf1.crossed(&s1s2)),
            );
            self.dx.set_value(4, 1, temp * psca_inv);

            let mut temp = f4 * (d1vtsurf1.dot(&s1s2) - self.tsurf1.dot(&d1v1));
            temp += np.dot(
                &self
                    .tsurf1
                    .crossed(&d1v1)
                    .subtracted(&d1vtsurf1.crossed(&s1s2)),
            );
            self.dx.set_value(4, 2, temp * psca_inv);

            let temp = f4 * self.tsurf1.dot(&d1u2) - np.dot(&self.tsurf1.crossed(&d1u2));
            self.dx.set_value(4, 3, temp * psca_inv);

            let temp = f4 * self.tsurf1.dot(&d1v2) - np.dot(&self.tsurf1.crossed(&d1v2));
            self.dx.set_value(4, 4, temp * psca_inv);
        }

        true
    }

    /// `IsSolution(Sol, Tol)` (`cxx:114-205`): the four residuals are within
    /// `Tol` (`2 * dist1 * Tol` for the third, the scaled `temp` for the
    /// fourth), then the section tangents are solved for with `math_Gauss`
    /// (pivot `1.e-9`) and, on failure, `math_SVD` with `Eps = 1.e-6`.
    pub(super) fn check_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        let mut valsol = MathVector::new(1, 4);
        let mut secmember = MathVector::new(1, 4);
        let mut gradsol = MathMatrix::new(1, 4, 1, 4);

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

        let mut temp = 2.0 * (nordu1 + nordv1) * s1s2.magnitude() + 2.0 * nordu1 * nordv1;

        self.values(sol, &mut valsol, &mut gradsol);

        if valsol.value(1).abs() < tol
            && valsol.value(2).abs() < tol
            && valsol.value(3).abs() < 2.0 * self.dist1 * tol
            && valsol.value(4).abs() < tol * (1.0 + self.tgang) * psca_inv.abs() * temp
        {
            secmember.set_value(1, normg - dnp.dot(&pguis1));
            secmember.set_value(2, normg - dnp.dot(&pguis2));
            secmember.set_value(3, -2.0 * d1gui.dot(&pguis1));

            let cross_vec = self.tsurf1.crossed(&s1s2);
            let f4 = np.dot(&cross_vec) * psca_inv;
            temp = dnp.dot(&cross_vec) + np.dot(&dwtsurf1.crossed(&s1s2));
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

            return true;
        }

        self.istangent = true;
        false
    }

    /// `Tangent(U1, V1, U2, V2, TgF, TgL, NmF, NmL)` (`cxx:428-481`): the
    /// tangents are reversed for the quadrants `2`/`5` (both), `4`/`7` (last)
    /// and `3`/`8` (first).
    #[allow(clippy::too_many_arguments)]
    pub(super) fn compute_tangent(
        &self,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        tg_first: &mut GpVec,
        tg_last: &mut GpVec,
        norm_first: &mut GpVec,
        norm_last: &mut GpVec,
    ) {
        let (_ptgui, d1gui) = self.curv.d1(self.param);
        let np = d1gui.normalized();

        let (_pt1, d1u1, d1v1) = self.surf1.d1(u1, v1);
        *norm_first = d1u1.crossed(&d1v1);

        let (_pt2, d1u2, d1v2) = self.surf2.d1(u2, v2);
        *norm_last = d1u2.crossed(&d1v2);

        *tg_first = np.crossed(norm_first).normalized();
        *tg_last = np.crossed(norm_last).normalized();

        let mut rev_first = false;
        let mut rev_last = false;

        if self.choix == 2 || self.choix == 5 {
            rev_first = true;
            rev_last = true;
        }

        if self.choix == 4 || self.choix == 7 {
            rev_last = true;
        }

        if self.choix == 3 || self.choix == 8 {
            rev_first = true;
        }

        if rev_first {
            tg_first.reverse();
        }
        if rev_last {
            tg_last.reverse();
        }
    }

    /// `Section(Param, U1, V1, U2, V2, Pdeb, Pfin, C)` (`cxx:484-502`): the
    /// line from the point on the first support to the point on the second,
    /// parameterised from `0` to the parameter of the second point.
    pub fn analytic_section(
        &self,
        _param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) -> (f64, f64, GpLin) {
        let pt1 = self.surf1.value(u1, v1);
        let pt2 = self.surf2.value(u2, v2);
        let dir = GpDir::from_vec(&GpVec::from_pnts(&pt1, &pt2))
            .expect("BlendFunc_ChAsym::Section: degenerate segment");
        let mut lin = GpLin::from_pnt_dir(pt1, dir);

        let pdeb = 0.0;
        let pfin = clib::line_parameter(&lin, &pt2);
        lin.set_location(pt1);
        (pdeb, pfin, lin)
    }
}

impl MathFunctionSetWithDerivatives for ChAsym {
    /// `NbVariables()` = `4`, not overridden.
    fn nb_variables(&self) -> usize {
        4
    }

    /// `NbEquations()` (`cxx:54-57`).
    fn nb_equations(&self) -> usize {
        4
    }

    /// `Value(X, F)` (`cxx:310-315`): `ComputeValues(X, 0, 0)`.
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        let error = self.compute_values(x, 0, 0);
        *f = self.fx.clone();
        error
    }

    /// `Derivatives(X, D)` (`cxx:319-324`): `ComputeValues(X, 1, 1)`.
    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        let error = self.compute_values(x, 1, 1);
        *d = self.dx.clone();
        error
    }

    /// `Values(X, F, D)` (`cxx:328-334`): `ComputeValues(X, 0, 1)`.
    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        let error = self.compute_values(x, 0, 1);
        *f = self.fx.clone();
        *d = self.dx.clone();
        error
    }
}

impl BlendAppFunction for ChAsym {
    /// `Set(Param)` (`cxx:60-63`).
    fn set(&mut self, param: f64) {
        self.param = param;
    }

    /// `Set(First, Last)` (`cxx:73-76`): `tcurv = curv->Trim(First, Last,
    /// 1.e-12)`. The trimmed curve keeps the basis parameterisation, so no
    /// member changes here (see the module notes).
    fn set_interval(&mut self, _first: f64, _last: f64) {}

    /// `GetTolerance(Tolerance, Tol)` (`cxx:79-85`).
    fn get_tolerance(&self, tolerance: &mut MathVector, tol: f64) {
        tolerance.set_value(1, approx_same_parameter::u_resolution(self.surf1.as_ref(), tol));
        tolerance.set_value(2, approx_same_parameter::v_resolution(self.surf1.as_ref(), tol));
        tolerance.set_value(3, approx_same_parameter::u_resolution(self.surf2.as_ref(), tol));
        tolerance.set_value(4, approx_same_parameter::v_resolution(self.surf2.as_ref(), tol));
    }

    /// `GetBounds(InfBound, SupBound)` (`cxx:88-110`): the box of both
    /// supports, each finite range widened by its own length.
    fn get_bounds(&self, inf_bound: &mut MathVector, sup_bound: &mut MathVector) {
        let (u1min, u1max) = self.surf1.u_range();
        let (v1min, v1max) = self.surf1.v_range();
        let (u2min, u2max) = self.surf2.u_range();
        let (v2min, v2max) = self.surf2.v_range();

        inf_bound.set_value(1, u1min);
        inf_bound.set_value(2, v1min);
        inf_bound.set_value(3, u2min);
        inf_bound.set_value(4, v2min);
        sup_bound.set_value(1, u1max);
        sup_bound.set_value(2, v1max);
        sup_bound.set_value(3, u2max);
        sup_bound.set_value(4, v2max);

        for i in 1..=4 {
            if !Precision::is_infinite(inf_bound.value(i))
                && !Precision::is_infinite(sup_bound.value(i))
            {
                let range = sup_bound.value(i) - inf_bound.value(i);
                inf_bound.set_value(i, inf_bound.value(i) - range);
                sup_bound.set_value(i, sup_bound.value(i) + range);
            }
        }
    }

    /// `IsSolution(Sol, Tol)` (`cxx:114-205`).
    fn is_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        self.check_solution(sol, tol)
    }

    /// `GetMinimalDistance()` (`cxx:210-213`).
    fn get_minimal_distance(&self) -> Result<f64, String> {
        Ok(self.distmin)
    }

    /// `Pnt1()` = `PointOnS1()` (`Blend_Function.cxx:25-28`).
    fn pnt1(&self) -> &GpPnt {
        &self.pt1
    }

    /// `Pnt2()` = `PointOnS2()` (`Blend_Function.cxx:30-33`).
    fn pnt2(&self) -> &GpPnt {
        &self.pt2
    }

    /// `IsRational()` (`cxx:506-509`): `false`.
    fn is_rational(&self) -> bool {
        false
    }

    /// `GetSectionSize()` throws
    /// `Standard_NotImplemented("BlendFunc_ChAsym::GetSectionSize()")`
    /// (`cxx:513-516`).
    fn get_section_size(&self) -> f64 {
        panic!("Standard_NotImplemented: BlendFunc_ChAsym::GetSectionSize()");
    }

    /// `GetMinimalWeight(Weights)` (`cxx:519-523`): all weights are `1`.
    fn get_minimal_weight(&self, weights: &mut [f64]) {
        weights.fill(1.0);
    }

    /// `NbIntervals(S)` (`cxx:526-529`).
    fn nb_intervals(&self, s: Shape) -> usize {
        self.curv.nb_intervals(next_shape(s) as u8) as usize
    }

    /// `Intervals(T, S)` (`cxx:532-536`).
    fn intervals(&self, t: &mut [f64], s: Shape) {
        let src = self.curv.parameter_intervals(next_shape(s) as u8);
        let count = t.len().min(src.len());
        t[..count].copy_from_slice(&src[..count]);
    }

    /// `GetShape(NbPoles, NbKnots, Degree, NbPoles2d)` (`cxx:543-549`).
    fn get_shape(
        &mut self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles2d: &mut i32,
    ) {
        *nb_poles = 2;
        *nb_poles2d = 2;
        *nb_knots = 2;
        *degree = 1;
    }

    /// `GetTolerance(BoundTol, SurfTol, AngleTol, Tol3d, Tol1d)`
    /// (`cxx:555-560`).
    fn get_tolerance_3d(
        &self,
        bound_tol: f64,
        _surf_tol: f64,
        _angle_tol: f64,
        tol3d: &mut MathVector,
        _tol1d: &mut MathVector,
    ) {
        tol3d.init(bound_tol);
    }

    /// `Knots(TKnots)` (`cxx:563-567`).
    fn knots(&mut self, tknots: &mut [f64]) {
        tknots[0] = 0.0;
        tknots[1] = 1.0;
    }

    /// `Mults(TMults)` (`cxx:572-576`).
    fn mults(&mut self, tmults: &mut [i32]) {
        tmults[0] = 2;
        tmults[1] = 2;
    }

    /// `Section(P, Poles, DPoles, Poles2d, DPoles2d, Weights, DWeights)`
    /// (`cxx:609-722`).
    #[allow(clippy::too_many_arguments)]
    fn section(
        &mut self,
        p: &BlendPoint,
        poles: &mut [GpPnt],
        dpoles: &mut [GpVec],
        poles2d: &mut [GpPnt2d],
        dpoles2d: &mut [GpVec2d],
        weights: &mut [f64],
        dweights: &mut [f64],
    ) -> bool {
        self.section_first_last(p, poles, dpoles, poles2d, dpoles2d, weights, dweights)
    }

    /// `Section(P, Poles, Poles2d, Weights)` (`cxx:580-607`).
    fn section_poles(
        &mut self,
        p: &BlendPoint,
        poles: &mut [GpPnt],
        poles2d: &mut [GpPnt2d],
        weights: &mut [f64],
    ) {
        self.section_poles_impl(p, poles, poles2d, weights);
    }

    /// `Section(P, ..., D2Poles, ...)` returns `false`
    /// (`cxx:725-740`).
    #[allow(clippy::too_many_arguments)]
    fn section_with_derivatives2(
        &mut self,
        _p: &BlendPoint,
        _poles: &mut [GpPnt],
        _dpoles: &mut [GpVec],
        _d2poles: &mut [GpVec],
        _poles2d: &mut [GpPnt2d],
        _dpoles2d: &mut [GpVec2d],
        _d2poles2d: &mut [GpVec2d],
        _weights: &mut [f64],
        _dweights: &mut [f64],
        _d2weights: &mut [f64],
    ) -> bool {
        false
    }

    /// `Resolution(IC2d, Tol, TolU, TolV)` (`cxx:744-758`).
    fn resolution(&self, ic2d: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
        if ic2d == 1 {
            *tol_u = approx_same_parameter::u_resolution(self.surf1.as_ref(), tol);
            *tol_v = approx_same_parameter::v_resolution(self.surf1.as_ref(), tol);
        } else {
            *tol_u = approx_same_parameter::u_resolution(self.surf2.as_ref(), tol);
            *tol_v = approx_same_parameter::v_resolution(self.surf2.as_ref(), tol);
        }
    }
}

impl BlendFunction for ChAsym {
    /// `PointOnS1()` (`cxx:337-340`).
    fn point_on_s1(&self) -> &GpPnt {
        &self.pt1
    }

    /// `PointOnS2()` (`cxx:344-347`).
    fn point_on_s2(&self) -> &GpPnt {
        &self.pt2
    }

    /// `IsTangencyPoint()` (`cxx:351-354`).
    fn is_tangency_point(&self) -> bool {
        self.istangent
    }

    /// `TangentOnS1()` (`cxx:358-365`).
    fn tangent_on_s1(&self) -> &GpVec {
        &self.tg1
    }

    /// `Tangent2dOnS1()` (`cxx:369-376`).
    fn tangent2d_on_s1(&self) -> &GpVec2d {
        &self.tg12d
    }

    /// `TangentOnS2()` (`cxx:381-388`).
    fn tangent_on_s2(&self) -> &GpVec {
        &self.tg2
    }

    /// `Tangent2dOnS2()` (`cxx:391-398`).
    fn tangent2d_on_s2(&self) -> &GpVec2d {
        &self.tg22d
    }

    /// `Tangent(U1, V1, U2, V2, TgF, TgL, NmF, NmL)` (`cxx:428-481`).
    #[allow(clippy::too_many_arguments)]
    fn tangent(
        &self,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        tg_first: &mut GpVec,
        tg_last: &mut GpVec,
        norm_first: &mut GpVec,
        norm_last: &mut GpVec,
    ) {
        self.compute_tangent(u1, v1, u2, v2, tg_first, tg_last, norm_first, norm_last);
    }
}
