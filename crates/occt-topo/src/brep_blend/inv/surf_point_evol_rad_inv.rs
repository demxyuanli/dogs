//! `BRepBlend_SurfPointEvolRadInv` — solution of the evolving-radius
//! point/surface reframing on a given point of the curve.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_SurfPointEvolRadInv.hxx:32-84`,
//! `.cxx:17-269` (OCCT 8.0.0).
//!
//! Same three equations as [`super::surf_point_const_rad_inv`], with the radius
//! read from the evolution law at `w` instead of being a constant:
//! `ray = sg1 * tevol->Value(w)` (`cxx:58`).
//!
//! Translation notes:
//! - `Set(Choix)` (`cxx:38-56`) only sets `sg1`; unlike the constant-radius
//!   class it does not touch `ray`, which `Value` / `Derivatives` recompute on
//!   every call.
//! - `Derivatives` (`cxx:60-144`) and `Values` (`cxx:148-243`) both add the
//!   law-derivative term `dwref.SetLinearForm(1., dwref, dray * unsurnorm, temp)`
//!   (`cxx:126`, `cxx:222`) and subtract `2 * dray * ray` from `D(3,1)`
//!   (`cxx:131`, `cxx:227`) respectively. Those two lines are the only
//!   difference from the constant-radius class and are kept verbatim.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpVec};
use occt_geom::curve::Curve;
use occt_geom::law::LawFunction;
use occt_geom::surface::Surface;
use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use crate::blend::BlendSurfPointFuncInv;

use super::pnt_vec;

/// `class BRepBlend_SurfPointEvolRadInv : public Blend_SurfPointFuncInv`
/// (`hxx:32-84`).
pub struct SurfPointEvolRadInv {
    /// `surf` (`hxx:79`).
    surf: Arc<dyn Surface>,
    /// `curv` (`hxx:80`).
    curv: Arc<dyn Curve>,
    /// `tevol` (`hxx:81`).
    tevol: Arc<dyn LawFunction>,
    /// `point` (`hxx:82`).
    point: GpPnt,
    /// `ray` (`hxx:83`): written by `Value` / `Derivatives` / `Values`.
    ray: f64,
    /// `choix` (`hxx:84`).
    choix: i32,
    /// `sg1` (`hxx:85`): the sign the quadrant selects.
    sg1: f64,
}

impl SurfPointEvolRadInv {
    /// `BRepBlend_SurfPointEvolRadInv(S, C, Evol)` (`cxx:21-29`): `ray`,
    /// `choix` and `sg1` are left at `0.` / `0` / `0.` by the default
    /// constructors of the remaining members.
    pub fn new(surf: Arc<dyn Surface>, curv: Arc<dyn Curve>, evol: Arc<dyn LawFunction>) -> Self {
        Self {
            surf,
            curv,
            tevol: evol,
            point: GpPnt::zero(),
            ray: 0.0,
            choix: 0,
            sg1: 0.0,
        }
    }

    /// `Set(Choix)` (`cxx:38-56`).
    pub fn set_choix(&mut self, choix: i32) {
        self.choix = choix;
        match choix {
            1 | 2 => self.sg1 = -1.0,
            3 | 4 => self.sg1 = 1.0,
            _ => self.sg1 = -1.0,
        }
    }

    /// `ray` (`hxx:83`), read by `IsSolution`.
    pub fn ray(&self) -> f64 {
        self.ray
    }
}

impl MathFunctionSetWithDerivatives for SurfPointEvolRadInv {
    /// `Blend_SurfPointFuncInv::NbVariables()` = `3`
    /// (`Blend_SurfPointFuncInv.cxx:21-24`).
    fn nb_variables(&self) -> usize {
        3
    }

    /// `NbEquations()` (`cxx:60-63`): `3`.
    fn nb_equations(&self) -> usize {
        3
    }

    /// `Value(X, F)` (`cxx:65-87`).
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        let (ptcur, d1cur) = self.curv.d1(x.value(1));
        self.ray = self.sg1 * self.tevol.value(x.value(1));
        let nplan = d1cur.normalized();
        // `theD = -(nplan.Dot(ptcur.XYZ()))` (`cxx:71`).
        let the_d = -nplan.dot(&pnt_vec(&ptcur));

        let (pts, d1u, d1v) = self.surf.d1(x.value(2), x.value(3));
        f.set_value(1, nplan.dot(&pnt_vec(&self.point)) + the_d);
        f.set_value(2, nplan.dot(&pnt_vec(&pts)) + the_d);

        let mut ns = d1u.crossed(&d1v);
        let norm = nplan.crossed(&ns).magnitude();
        let unsurnorm = 1.0 / norm;
        ns = ns.set_linear_form_2(nplan.dot(&ns), &nplan, -1.0, &ns);
        ns = ns.multiplied_scalar(unsurnorm);
        let mut r = pnt_vec(&pts).subtracted(&pnt_vec(&self.point));
        r = r.set_linear_form_2(self.ray, &ns, 1.0, &r);
        f.set_value(3, r.square_magnitude() - self.ray * self.ray);
        true
    }

    /// `Derivatives(X, D)` (`cxx:91-146`).
    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        let (ptcur, d1cur, d2cur) = self.curv.d2(x.value(1));
        let (ray, dray) = self.tevol.d1(x.value(1));
        self.ray = self.sg1 * ray;
        let dray = self.sg1 * dray;
        let normd1cur = d1cur.magnitude();
        let unsurnormd1cur = 1.0 / normd1cur;
        let nplan = d1cur.multiplied_scalar(unsurnormd1cur);
        let ptcur_v = pnt_vec(&ptcur);

        let mut dnplan = GpVec::zero().set_linear_form_2(-nplan.dot(&d2cur), &nplan, 1.0, &d2cur);
        dnplan = dnplan.multiplied_scalar(unsurnormd1cur);
        let dthe_d = -nplan.dot(&d1cur) - dnplan.dot(&ptcur_v);
        d.set_value(1, 1, dnplan.dot(&pnt_vec(&self.point)) + dthe_d);
        d.set_value(1, 2, 0.0);
        d.set_value(1, 3, 0.0);

        let (pts, d1u, d1v, d2u, d2v, duv) = self.surf.d2(x.value(2), x.value(3));
        let pts_v = pnt_vec(&pts);
        d.set_value(2, 1, dnplan.dot(&pts_v) + dthe_d);
        d.set_value(2, 2, nplan.dot(&d1u));
        d.set_value(2, 3, nplan.dot(&d1v));

        let nsurf = d1u.crossed(&d1v);
        let dunsurf = d2u.crossed(&d1v).added(&d1u.crossed(&duv));
        let dvnsurf = d1u.crossed(&d2v).added(&duv.crossed(&d1v));

        let nplancrosnsurf = nplan.crossed(&nsurf);
        let dwnplancrosnsurf = dnplan.crossed(&nsurf);
        let dunplancrosnsurf = nplan.crossed(&dunsurf);
        let dvnplancrosnsurf = nplan.crossed(&dvnsurf);

        let norm2 = nplancrosnsurf.square_magnitude();
        let norm = norm2.sqrt();
        let unsurnorm = 1.0 / norm;
        let raysurnorm = self.ray * unsurnorm;
        let unsurnorm2 = unsurnorm * unsurnorm;
        let raysurnorm2 = self.ray * unsurnorm2;
        let dwnorm = unsurnorm * nplancrosnsurf.dot(&dwnplancrosnsurf);
        let dunorm = unsurnorm * nplancrosnsurf.dot(&dunplancrosnsurf);
        let dvnorm = unsurnorm * nplancrosnsurf.dot(&dvnplancrosnsurf);

        let nplandotnsurf = nplan.dot(&nsurf);
        let dwnplandotnsurf = dnplan.dot(&nsurf);
        let dunplandotnsurf = nplan.dot(&dunsurf);
        let dvnplandotnsurf = nplan.dot(&dvnsurf);

        let temp = GpVec::zero().set_linear_form_2(nplandotnsurf, &nplan, -1.0, &nsurf);
        let dwtemp = GpVec::zero().set_linear_form_2(nplandotnsurf, &dnplan, dwnplandotnsurf, &nplan);
        let dutemp = GpVec::zero().set_linear_form_2(dunplandotnsurf, &nplan, -1.0, &dunsurf);
        let dvtemp = GpVec::zero().set_linear_form_2(dvnplandotnsurf, &nplan, -1.0, &dvnsurf);

        let corde = GpVec::from_pnts(&self.point, &pts);
        let mut r = GpVec::zero().set_linear_form_2(raysurnorm, &temp, 1.0, &corde);
        let mut dwref =
            GpVec::zero().set_linear_form_2(raysurnorm, &dwtemp, -raysurnorm2 * dwnorm, &temp);
        // `dwref.SetLinearForm(1., dwref, dray * unsurnorm, temp)` (`cxx:126`).
        let dwcopy = dwref;
        dwref = dwcopy.set_linear_form_2(1.0, &dwcopy, dray * unsurnorm, &temp);
        let duref = GpVec::zero().set_linear_form_2plus(
            raysurnorm,
            &dutemp,
            -raysurnorm2 * dunorm,
            &temp,
            &d1u,
        );
        let dvref = GpVec::zero().set_linear_form_2plus(
            raysurnorm,
            &dvtemp,
            -raysurnorm2 * dvnorm,
            &temp,
            &d1v,
        );

        r = r.added(&r);
        d.set_value(3, 1, r.dot(&dwref) - 2.0 * dray * self.ray);
        d.set_value(3, 2, r.dot(&duref));
        d.set_value(3, 3, r.dot(&dvref));

        true
    }

    /// `Values(X, F, D)` (`cxx:148-243`): `Value` and `Derivatives` inlined.
    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        let (ptcur, d1cur, d2cur) = self.curv.d2(x.value(1));
        let (ray, dray) = self.tevol.d1(x.value(1));
        self.ray = self.sg1 * ray;
        let dray = self.sg1 * dray;
        let (pts, d1u, d1v, d2u, d2v, duv) = self.surf.d2(x.value(2), x.value(3));
        let normd1cur = d1cur.magnitude();
        let unsurnormd1cur = 1.0 / normd1cur;
        let nplan = d1cur.multiplied_scalar(unsurnormd1cur);
        let the_d = -nplan.dot(&pnt_vec(&ptcur));
        let pts_v = pnt_vec(&pts);

        f.set_value(1, nplan.dot(&pnt_vec(&self.point)) + the_d);
        f.set_value(2, nplan.dot(&pts_v) + the_d);

        let mut dnplan = GpVec::zero().set_linear_form_2(-nplan.dot(&d2cur), &nplan, 1.0, &d2cur);
        dnplan = dnplan.multiplied_scalar(unsurnormd1cur);
        let dthe_d = -nplan.dot(&d1cur) - dnplan.dot(&pnt_vec(&ptcur));
        d.set_value(1, 1, dnplan.dot(&pnt_vec(&self.point)) + dthe_d);
        d.set_value(1, 2, 0.0);
        d.set_value(1, 3, 0.0);
        d.set_value(2, 1, dnplan.dot(&pts_v) + dthe_d);
        d.set_value(2, 2, nplan.dot(&d1u));
        d.set_value(2, 3, nplan.dot(&d1v));

        let nsurf = d1u.crossed(&d1v);
        let dunsurf = d2u.crossed(&d1v).added(&d1u.crossed(&duv));
        let dvnsurf = d1u.crossed(&d2v).added(&duv.crossed(&d1v));

        let nplancrosnsurf = nplan.crossed(&nsurf);
        let dwnplancrosnsurf = dnplan.crossed(&nsurf);
        let dunplancrosnsurf = nplan.crossed(&dunsurf);
        let dvnplancrosnsurf = nplan.crossed(&dvnsurf);

        let norm2 = nplancrosnsurf.square_magnitude();
        let norm = norm2.sqrt();
        let unsurnorm = 1.0 / norm;
        let raysurnorm = self.ray * unsurnorm;
        let unsurnorm2 = unsurnorm * unsurnorm;
        let raysurnorm2 = self.ray * unsurnorm2;
        let dwnorm = unsurnorm * nplancrosnsurf.dot(&dwnplancrosnsurf);
        let dunorm = unsurnorm * nplancrosnsurf.dot(&dunplancrosnsurf);
        let dvnorm = unsurnorm * nplancrosnsurf.dot(&dvnplancrosnsurf);

        let nplandotnsurf = nplan.dot(&nsurf);
        let dwnplandotnsurf = dnplan.dot(&nsurf);
        let dunplandotnsurf = nplan.dot(&dunsurf);
        let dvnplandotnsurf = nplan.dot(&dvnsurf);

        let temp = GpVec::zero().set_linear_form_2(nplandotnsurf, &nplan, -1.0, &nsurf);
        let dwtemp = GpVec::zero().set_linear_form_2(nplandotnsurf, &dnplan, dwnplandotnsurf, &nplan);
        let dutemp = GpVec::zero().set_linear_form_2(dunplandotnsurf, &nplan, -1.0, &dunsurf);
        let dvtemp = GpVec::zero().set_linear_form_2(dvnplandotnsurf, &nplan, -1.0, &dvnsurf);

        let corde = GpVec::from_pnts(&self.point, &pts);
        let mut r = GpVec::zero().set_linear_form_2(raysurnorm, &temp, 1.0, &corde);
        f.set_value(3, r.square_magnitude() - self.ray * self.ray);
        let mut dwref =
            GpVec::zero().set_linear_form_2(raysurnorm, &dwtemp, -raysurnorm2 * dwnorm, &temp);
        // `dwref.SetLinearForm(1., dwref, dray * unsurnorm, temp)` (`cxx:222`).
        let dwcopy = dwref;
        dwref = dwcopy.set_linear_form_2(1.0, &dwcopy, dray * unsurnorm, &temp);
        let duref = GpVec::zero().set_linear_form_2plus(
            raysurnorm,
            &dutemp,
            -raysurnorm2 * dunorm,
            &temp,
            &d1u,
        );
        let dvref = GpVec::zero().set_linear_form_2plus(
            raysurnorm,
            &dvtemp,
            -raysurnorm2 * dvnorm,
            &temp,
            &d1v,
        );

        r = r.added(&r);
        d.set_value(3, 1, r.dot(&dwref) - 2.0 * dray * self.ray);
        d.set_value(3, 2, r.dot(&duref));
        d.set_value(3, 3, r.dot(&dvref));

        true
    }
}

impl BlendSurfPointFuncInv for SurfPointEvolRadInv {
    /// `Set(P)` (`cxx:247-250`).
    fn set(&mut self, p: GpPnt) {
        self.point = p;
    }

    /// `GetTolerance(Tolerance, Tol)` (`cxx:252-258`).
    fn get_tolerance(&self, tolerance: &mut MathVector, tol: f64) {
        tolerance.set_value(1, self.curv.resolution(tol));
        tolerance.set_value(
            2,
            occt_geom::approx_same_parameter::u_resolution(self.surf.as_ref(), tol),
        );
        tolerance.set_value(
            3,
            occt_geom::approx_same_parameter::v_resolution(self.surf.as_ref(), tol),
        );
    }

    /// `GetBounds(InfBound, SupBound)` (`cxx:260-273`).
    fn get_bounds(&self, inf_bound: &mut MathVector, sup_bound: &mut MathVector) {
        inf_bound.set_value(1, self.curv.first_parameter());
        sup_bound.set_value(1, self.curv.last_parameter());
        let (u1, u2) = self.surf.u_range();
        inf_bound.set_value(2, u1);
        sup_bound.set_value(2, u2);
        let (v1, v2) = self.surf.v_range();
        inf_bound.set_value(3, v1);
        sup_bound.set_value(3, v2);
    }

    /// `IsSolution(Sol, Tol)` (`cxx:275-282`).
    fn is_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        let mut valsol = MathVector::new(1, 3);
        self.value(sol, &mut valsol);
        valsol.value(1).abs() <= tol
            && valsol.value(2).abs() <= tol
            && valsol.value(3).abs() <= 2.0 * tol * self.ray.abs()
    }
}
