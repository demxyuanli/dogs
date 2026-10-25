//! `BRepBlend_SurfCurvEvolRadInv` — solution of the evolving-radius
//! surface/curve reframing on a given restriction of the surface.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_SurfCurvEvolRadInv.hxx:33-90`,
//! `.cxx:17-303` (OCCT 8.0.0).
//!
//! Same three equations as [`super::surf_curv_const_rad_inv`], with
//! `ray = sg1 * tevol->Value(wguide)` (`cxx:63`).
//!
//! Translation notes:
//! - `Set(Choix)` (`cxx:35-53`) only sets `sg1`; `ray` is recomputed on every
//!   evaluation.
//! - `Derivatives` (`cxx:57-152`) and `Values` (`cxx:156-262`) add
//!   `dwguiref.SetLinearForm(1., dwguiref, dray * unsurnorm, temp)`
//!   (`cxx:136`, `cxx:250`) and subtract `2 * dray * ray` from `D(3,1)`
//!   (`cxx:141`, `cxx:255`): the law-derivative terms that the constant-radius
//!   class does not have.
//! - `D(3, 2) = -ref.Dot(d1cur)` keeps the sign of `corde =
//!   gp_Vec(ptcur, pts)`.
//! - `rst` is stored as `Option` for the empty handle, as in
//!   [`super::surf_curv_const_rad_inv`].

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpVec};
use occt_geom::curve::Curve;
use occt_geom::law::LawFunction;
use occt_geom::surface::Surface;
use occt_geom2d::curve::Curve2d;
use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use crate::blend::BlendSurfCurvFuncInv;

use super::pnt_vec;

/// `class BRepBlend_SurfCurvEvolRadInv : public Blend_SurfCurvFuncInv`
/// (`hxx:33-90`).
pub struct SurfCurvEvolRadInv {
    /// `surf` (`hxx:79`).
    surf: Arc<dyn Surface>,
    /// `curv` (`hxx:80`).
    curv: Arc<dyn Curve>,
    /// `guide` (`hxx:81`).
    guide: Arc<dyn Curve>,
    /// `rst` (`hxx:82`): `occ::handle<Adaptor2d_Curve2d>`.
    rst: Option<Arc<dyn Curve2d>>,
    /// `tevol` (`hxx:83`).
    tevol: Arc<dyn LawFunction>,
    /// `ray` (`hxx:84`).
    ray: f64,
    /// `choix` (`hxx:85`).
    choix: i32,
    /// `sg1` (`hxx:86`).
    sg1: f64,
}

impl SurfCurvEvolRadInv {
    /// `BRepBlend_SurfCurvEvolRadInv(S, C, Cg, Evol)` (`cxx:21-32`): `tevol`
    /// assigned, `rst` null, `ray` / `choix` / `sg1` left at their member
    /// initialisers.
    pub fn new(
        surf: Arc<dyn Surface>,
        curv: Arc<dyn Curve>,
        guide: Arc<dyn Curve>,
        evol: Arc<dyn LawFunction>,
    ) -> Self {
        Self {
            surf,
            curv,
            guide,
            rst: None,
            tevol: evol,
            ray: 0.0,
            choix: 0,
            sg1: 0.0,
        }
    }

    /// `Set(Choix)` (`cxx:35-53`).
    pub fn set_choix(&mut self, choix: i32) {
        self.choix = choix;
        match choix {
            1 | 2 => self.sg1 = -1.0,
            3 | 4 => self.sg1 = 1.0,
            _ => self.sg1 = -1.0,
        }
    }

    /// `ray` (`hxx:84`), read by `IsSolution`.
    pub fn ray(&self) -> f64 {
        self.ray
    }
}

impl MathFunctionSetWithDerivatives for SurfCurvEvolRadInv {
    /// `Blend_SurfCurvFuncInv::NbVariables()` = `3`
    /// (`Blend_SurfCurvFuncInv.cxx:22-25`).
    fn nb_variables(&self) -> usize {
        3
    }

    /// `NbEquations()` (`cxx:57-60`): `3`.
    fn nb_equations(&self) -> usize {
        3
    }

    /// `Value(X, F)` (`cxx:62-90`).
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        let Some(rst) = self.rst.clone() else {
            return false;
        };
        let (ptgui, d1gui) = self.guide.d1(x.value(1));
        self.ray = self.sg1 * self.tevol.value(x.value(1));
        let nplan = d1gui.normalized();
        // `theD = nplan.Dot(ptgui.XYZ()); theD = theD * (-1.)` (`cxx:68-72`).
        let the_d = -nplan.dot(&pnt_vec(&ptgui));

        let ptcur = self.curv.value(x.value(2));
        f.set_value(1, nplan.dot(&pnt_vec(&ptcur)) + the_d);

        let p2drst = rst.value(x.value(3));
        let (pts, du, dv) = self.surf.d1(p2drst.x(), p2drst.y());
        f.set_value(2, nplan.dot(&pnt_vec(&pts)) + the_d);

        let mut ns = du.crossed(&dv);
        let norm = nplan.crossed(&ns).magnitude();
        let unsurnorm = 1.0 / norm;
        ns = ns.set_linear_form_2(nplan.dot(&ns), &nplan, -1.0, &ns);
        ns = ns.multiplied_scalar(unsurnorm);
        let mut r = GpVec::from_pnts(&ptcur, &pts);
        r = r.set_linear_form_2(self.ray, &ns, 1.0, &r);
        f.set_value(3, r.square_magnitude() - self.ray * self.ray);
        true
    }

    /// `Derivatives(X, D)` (`cxx:94-153`).
    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        let Some(rst) = self.rst.clone() else {
            return false;
        };
        let (ptgui, d1gui, d2gui) = self.guide.d2(x.value(1));
        let normd1gui = d1gui.magnitude();
        let unsurnormd1gui = 1.0 / normd1gui;
        let (ray, dray) = self.tevol.d1(x.value(1));
        self.ray = self.sg1 * ray;
        let dray = self.sg1 * dray;
        let nplan = d1gui.multiplied_scalar(unsurnormd1gui);
        let mut dnplan = GpVec::zero().set_linear_form_2(-nplan.dot(&d2gui), &nplan, 1.0, &d2gui);
        dnplan = dnplan.multiplied_scalar(unsurnormd1gui);
        let dthe_d = -nplan.dot(&d1gui) - dnplan.dot(&pnt_vec(&ptgui));

        let (ptcur, d1cur) = self.curv.d1(x.value(2));
        d.set_value(1, 1, dnplan.dot(&pnt_vec(&ptcur)) + dthe_d);
        d.set_value(1, 2, nplan.dot(&d1cur));
        d.set_value(1, 3, 0.0);

        let (p2drst, d1rst) = rst.d1(x.value(3));
        let (pts, d1u, d1v, d2u, d2v, duv) = self.surf.d2(p2drst.x(), p2drst.y());
        let pts_v = pnt_vec(&pts);
        d.set_value(2, 1, dnplan.dot(&pts_v) + dthe_d);
        d.set_value(2, 2, 0.0);
        let dwrstpts = GpVec::zero().set_linear_form_2(d1rst.x(), &d1u, d1rst.y(), &d1v);
        d.set_value(2, 3, nplan.dot(&dwrstpts));

        let nsurf = d1u.crossed(&d1v);
        let dunsurf = d2u.crossed(&d1v).added(&d1u.crossed(&duv));
        let dvnsurf = d1u.crossed(&d2v).added(&duv.crossed(&d1v));
        let dwrstnsurf = GpVec::zero().set_linear_form_2(d1rst.x(), &dunsurf, d1rst.y(), &dvnsurf);

        let nplancrosnsurf = nplan.crossed(&nsurf);
        let dwguinplancrosnsurf = dnplan.crossed(&nsurf);
        let dwrstnplancrosnsurf = nplan.crossed(&dwrstnsurf);

        let norm2 = nplancrosnsurf.square_magnitude();
        let norm = norm2.sqrt();
        let unsurnorm = 1.0 / norm;
        let raysurnorm = self.ray * unsurnorm;
        let unsurnorm2 = unsurnorm * unsurnorm;
        let raysurnorm2 = self.ray * unsurnorm2;
        let dwguinorm = unsurnorm * nplancrosnsurf.dot(&dwguinplancrosnsurf);
        let dwrstnorm = unsurnorm * nplancrosnsurf.dot(&dwrstnplancrosnsurf);

        let nplandotnsurf = nplan.dot(&nsurf);
        let dwguinplandotnsurf = dnplan.dot(&nsurf);
        let dwrstnplandotnsurf = nplan.dot(&dwrstnsurf);

        let temp = GpVec::zero().set_linear_form_2(nplandotnsurf, &nplan, -1.0, &nsurf);
        let dwguitemp =
            GpVec::zero().set_linear_form_2(nplandotnsurf, &dnplan, dwguinplandotnsurf, &nplan);
        let dwrsttemp =
            GpVec::zero().set_linear_form_2(dwrstnplandotnsurf, &nplan, -1.0, &dwrstnsurf);

        let corde = GpVec::from_pnts(&ptcur, &pts);
        let mut r = GpVec::zero().set_linear_form_2(raysurnorm, &temp, 1.0, &corde);
        let mut dwguiref =
            GpVec::zero().set_linear_form_2(raysurnorm, &dwguitemp, -raysurnorm2 * dwguinorm, &temp);
        // `dwguiref.SetLinearForm(1., dwguiref, dray * unsurnorm, temp)` (`cxx:136`).
        let dwcopy = dwguiref;
        dwguiref = dwcopy.set_linear_form_2(1.0, &dwcopy, dray * unsurnorm, &temp);
        let dwrstref = GpVec::zero().set_linear_form_2plus(
            raysurnorm,
            &dwrsttemp,
            -raysurnorm2 * dwrstnorm,
            &temp,
            &dwrstpts,
        );

        r = r.added(&r);
        d.set_value(3, 1, r.dot(&dwguiref) - 2.0 * dray * self.ray);
        d.set_value(3, 2, -r.dot(&d1cur));
        d.set_value(3, 3, r.dot(&dwrstref));

        true
    }

    /// `Values(X, F, D)` (`cxx:156-262`): `Value` and `Derivatives` inlined.
    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        let Some(rst) = self.rst.clone() else {
            return false;
        };
        let (ptgui, d1gui, d2gui) = self.guide.d2(x.value(1));
        let (ray, dray) = self.tevol.d1(x.value(1));
        self.ray = self.sg1 * ray;
        let dray = self.sg1 * dray;
        let normd1gui = d1gui.magnitude();
        let unsurnormd1gui = 1.0 / normd1gui;
        let nplan = d1gui.multiplied_scalar(unsurnormd1gui);
        let the_d = -nplan.dot(&pnt_vec(&ptgui));
        let mut dnplan = GpVec::zero().set_linear_form_2(-nplan.dot(&d2gui), &nplan, 1.0, &d2gui);
        dnplan = dnplan.multiplied_scalar(unsurnormd1gui);
        let dthe_d = -nplan.dot(&d1gui) - dnplan.dot(&pnt_vec(&ptgui));

        let (ptcur, d1cur) = self.curv.d1(x.value(2));
        f.set_value(1, nplan.dot(&pnt_vec(&ptcur)) + the_d);
        d.set_value(1, 1, dnplan.dot(&pnt_vec(&ptcur)) + dthe_d);
        d.set_value(1, 2, nplan.dot(&d1cur));
        d.set_value(1, 3, 0.0);

        let (p2drst, d1rst) = rst.d1(x.value(3));
        let (pts, d1u, d1v, d2u, d2v, duv) = self.surf.d2(p2drst.x(), p2drst.y());
        let pts_v = pnt_vec(&pts);
        f.set_value(2, nplan.dot(&pts_v) + the_d);
        d.set_value(2, 1, dnplan.dot(&pts_v) + dthe_d);
        d.set_value(2, 2, 0.0);
        let dwrstpts = GpVec::zero().set_linear_form_2(d1rst.x(), &d1u, d1rst.y(), &d1v);
        d.set_value(2, 3, nplan.dot(&dwrstpts));

        let nsurf = d1u.crossed(&d1v);
        let dunsurf = d2u.crossed(&d1v).added(&d1u.crossed(&duv));
        let dvnsurf = d1u.crossed(&d2v).added(&duv.crossed(&d1v));
        let dwrstnsurf = GpVec::zero().set_linear_form_2(d1rst.x(), &dunsurf, d1rst.y(), &dvnsurf);

        let nplancrosnsurf = nplan.crossed(&nsurf);
        let dwguinplancrosnsurf = dnplan.crossed(&nsurf);
        let dwrstnplancrosnsurf = nplan.crossed(&dwrstnsurf);

        let norm2 = nplancrosnsurf.square_magnitude();
        let norm = norm2.sqrt();
        let unsurnorm = 1.0 / norm;
        let raysurnorm = self.ray * unsurnorm;
        let unsurnorm2 = unsurnorm * unsurnorm;
        let raysurnorm2 = self.ray * unsurnorm2;
        let dwguinorm = unsurnorm * nplancrosnsurf.dot(&dwguinplancrosnsurf);
        let dwrstnorm = unsurnorm * nplancrosnsurf.dot(&dwrstnplancrosnsurf);

        let nplandotnsurf = nplan.dot(&nsurf);
        let dwguinplandotnsurf = dnplan.dot(&nsurf);
        let dwrstnplandotnsurf = nplan.dot(&dwrstnsurf);

        let temp = GpVec::zero().set_linear_form_2(nplandotnsurf, &nplan, -1.0, &nsurf);
        let dwguitemp =
            GpVec::zero().set_linear_form_2(nplandotnsurf, &dnplan, dwguinplandotnsurf, &nplan);
        let dwrsttemp =
            GpVec::zero().set_linear_form_2(dwrstnplandotnsurf, &nplan, -1.0, &dwrstnsurf);

        let corde = GpVec::from_pnts(&ptcur, &pts);
        let mut r = GpVec::zero().set_linear_form_2(raysurnorm, &temp, 1.0, &corde);
        f.set_value(3, r.square_magnitude() - self.ray * self.ray);
        let mut dwguiref =
            GpVec::zero().set_linear_form_2(raysurnorm, &dwguitemp, -raysurnorm2 * dwguinorm, &temp);
        // `dwguiref.SetLinearForm(1., dwguiref, dray * unsurnorm, temp)` (`cxx:250`).
        let dwcopy = dwguiref;
        dwguiref = dwcopy.set_linear_form_2(1.0, &dwcopy, dray * unsurnorm, &temp);
        let dwrstref = GpVec::zero().set_linear_form_2plus(
            raysurnorm,
            &dwrsttemp,
            -raysurnorm2 * dwrstnorm,
            &temp,
            &dwrstpts,
        );

        r = r.added(&r);
        d.set_value(3, 1, r.dot(&dwguiref) - 2.0 * dray * self.ray);
        d.set_value(3, 2, -r.dot(&d1cur));
        d.set_value(3, 3, r.dot(&dwrstref));
        true
    }
}

impl BlendSurfCurvFuncInv for SurfCurvEvolRadInv {
    /// `Set(Rst)` (`cxx:266-269`).
    fn set(&mut self, rst: Option<Arc<dyn Curve2d>>) {
        self.rst = rst;
    }

    /// `GetTolerance(Tolerance, Tol)` (`cxx:271-280`).
    fn get_tolerance(&self, tolerance: &mut MathVector, tol: f64) {
        tolerance.set_value(1, self.guide.resolution(tol));
        tolerance.set_value(2, self.curv.resolution(tol));
        let ru = occt_geom::approx_same_parameter::u_resolution(self.surf.as_ref(), tol);
        let rv = occt_geom::approx_same_parameter::v_resolution(self.surf.as_ref(), tol);
        let r = match &self.rst {
            Some(rst) => rst.resolution(ru.min(rv)),
            None => 0.0,
        };
        tolerance.set_value(3, r);
    }

    /// `GetBounds(InfBound, SupBound)` (`cxx:282-294`).
    fn get_bounds(&self, inf_bound: &mut MathVector, sup_bound: &mut MathVector) {
        inf_bound.set_value(1, self.guide.first_parameter());
        sup_bound.set_value(1, self.guide.last_parameter());
        inf_bound.set_value(2, self.curv.first_parameter());
        sup_bound.set_value(2, self.curv.last_parameter());
        let (f, l) = match &self.rst {
            Some(rst) => (rst.first_parameter(), rst.last_parameter()),
            None => (0.0, 0.0),
        };
        inf_bound.set_value(3, f);
        sup_bound.set_value(3, l);
    }

    /// `IsSolution(Sol, Tol)` (`cxx:296-303`).
    fn is_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        let mut valsol = MathVector::new(1, 3);
        self.value(sol, &mut valsol);
        valsol.value(1).abs() <= tol
            && valsol.value(2).abs() <= tol
            && valsol.value(3).abs() <= 2.0 * tol * self.ray.abs()
    }
}
