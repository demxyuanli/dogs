//! `BRepBlend_SurfPointConstRadInv` — solution of the constant-radius
//! point/surface reframing on a given point of the curve.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_SurfPointConstRadInv.hxx:32-84`,
//! `.cxx:17-276` (OCCT 8.0.0).
//!
//! The vector `X` is `(w, U, V)`: `w` the parameter on the guide curve, `U, V`
//! the parameters on the partner surface (`hxx:42-46`).
//!
//! The three equations are, with `nplan` the unit tangent of the guide curve at
//! `w` and `theD = -(nplan . ptcur)`:
//! 1. `nplan . point + theD` — the given point lies in the section plane,
//! 2. `nplan . pts + theD` — the surface point lies in the same plane,
//! 3. `|corde + ray * ns|^2 - ray^2` — the surface point is at distance
//!    `|ray|` from the point, `ns` being the surface normal projected into the
//!    section plane.
//!
//! Translation notes:
//! - `Set(R, Choix)` picks the sign of `ray` from the fillet quadrant
//!   (`cxx:33-51`): `1`/`2` and the `default:` arm give `-|R|`, `3`/`4` give
//!   `+|R|`. The port keeps the switch, not a sign formula.
//! - `Value` (`cxx:60-82`) and `Derivatives` (`cxx:86-176`) are written
//!   separately in OCCT, with `Values` (`cxx:180-262`) inlining both; the port
//!   keeps the three methods, so `values` reproduces OCCT's inlined form (which
//!   expects `d2` from both supports) rather than calling `value` +
//!   `derivatives`.
//! - `ns.SetLinearForm(nplan.Dot(ns), nplan, -1., ns)` (`cxx:74`) becomes
//!   `set_linear_form_2`; `ref.SetLinearForm(ray, ns, ref)` (`cxx:79`) is the
//!   same call with coefficient `1.0` on `ref`.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpVec};
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use crate::blend::BlendSurfPointFuncInv;

use super::pnt_vec;

/// `class BRepBlend_SurfPointConstRadInv : public Blend_SurfPointFuncInv`
/// (`hxx:32-84`).
pub struct SurfPointConstRadInv {
    /// `surf` (`hxx:79`).
    surf: Arc<dyn Surface>,
    /// `curv` (`hxx:80`).
    curv: Arc<dyn Curve>,
    /// `point` (`hxx:81`).
    point: GpPnt,
    /// `ray` (`hxx:82`).
    ray: f64,
    /// `choix` (`hxx:83`). Only `Set(R, Choix)` writes it; no method reads it,
    /// exactly as in OCCT.
    choix: i32,
}

impl SurfPointConstRadInv {
    /// `BRepBlend_SurfPointConstRadInv(S, C)` (`cxx:21-29`): `ray = 0.`,
    /// `choix = 0`.
    pub fn new(surf: Arc<dyn Surface>, curv: Arc<dyn Curve>) -> Self {
        Self {
            surf,
            curv,
            point: GpPnt::zero(),
            ray: 0.0,
            choix: 0,
        }
    }

    /// `Set(R, Choix)` (`cxx:33-51`).
    pub fn set_radius(&mut self, r: f64, choix: i32) {
        self.choix = choix;
        match choix {
            1 | 2 => self.ray = -r.abs(),
            3 | 4 => self.ray = r.abs(),
            _ => self.ray = -r.abs(),
        }
    }

    /// `ray` (`hxx:82`), read by `IsSolution`.
    pub fn ray(&self) -> f64 {
        self.ray
    }
}

impl MathFunctionSetWithDerivatives for SurfPointConstRadInv {
    /// `Blend_SurfPointFuncInv::NbVariables()` = `3`
    /// (`Blend_SurfPointFuncInv.cxx:21-24`).
    fn nb_variables(&self) -> usize {
        3
    }

    /// `NbEquations()` (`cxx:55-58`): `3`.
    fn nb_equations(&self) -> usize {
        3
    }

    /// `Value(X, F)` (`cxx:60-82`).
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        let (ptcur, d1cur) = self.curv.d1(x.value(1));
        let nplan = d1cur.normalized();
        // `theD = nplan.Dot(ptcur.XYZ()); theD = theD * (-1.);` (`cxx:67-70`).
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

    /// `Derivatives(X, D)` (`cxx:86-176`).
    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        let (ptcur, d1cur, d2cur) = self.curv.d2(x.value(1));
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
        let dwref = GpVec::zero().set_linear_form_2(raysurnorm, &dwtemp, -raysurnorm2 * dwnorm, &temp);
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

        // `ref.Add(ref)` (`cxx:171`): `ref` becomes twice itself before the
        // three dot products, i.e. the gradient of `|ref|^2`.
        r = r.added(&r);
        d.set_value(3, 1, r.dot(&dwref));
        d.set_value(3, 2, r.dot(&duref));
        d.set_value(3, 3, r.dot(&dvref));

        true
    }

    /// `Values(X, F, D)` (`cxx:180-262`): `Value` and `Derivatives` inlined,
    /// with `d2` taken from both supports.
    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        let (ptcur, d1cur, d2cur) = self.curv.d2(x.value(1));
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
        let dwref = GpVec::zero().set_linear_form_2(raysurnorm, &dwtemp, -raysurnorm2 * dwnorm, &temp);
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
        d.set_value(3, 1, r.dot(&dwref));
        d.set_value(3, 2, r.dot(&duref));
        d.set_value(3, 3, r.dot(&dvref));

        true
    }
}

impl BlendSurfPointFuncInv for SurfPointConstRadInv {
    /// `Set(P)` (`cxx:266-269`).
    fn set(&mut self, p: GpPnt) {
        self.point = p;
    }

    /// `GetTolerance(Tolerance, Tol)` (`cxx:271-276`).
    fn get_tolerance(&self, tolerance: &mut MathVector, tol: f64) {
        tolerance.set_value(1, self.curv.resolution(tol));
        let (ru, rv) = uv_resolution(self.surf.as_ref(), tol);
        tolerance.set_value(2, ru);
        tolerance.set_value(3, rv);
    }

    /// `GetBounds(InfBound, SupBound)` (`cxx:250-262`).
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

    /// `IsSolution(Sol, Tol)` (`cxx:246-252`): the third equation is scaled by
    /// `2 * |ray|` because it is a squared distance.
    fn is_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        let mut valsol = MathVector::new(1, 3);
        self.value(sol, &mut valsol);
        valsol.value(1).abs() <= tol
            && valsol.value(2).abs() <= tol
            && valsol.value(3).abs() <= 2.0 * tol * self.ray.abs()
    }
}

/// `surf->UResolution(Tol)` / `surf->VResolution(Tol)`, the pair OCCT's
/// `GetTolerance` writes into slots 2 and 3.
fn uv_resolution(s: &dyn Surface, tol: f64) -> (f64, f64) {
    (
        occt_geom::approx_same_parameter::u_resolution(s, tol),
        occt_geom::approx_same_parameter::v_resolution(s, tol),
    )
}
