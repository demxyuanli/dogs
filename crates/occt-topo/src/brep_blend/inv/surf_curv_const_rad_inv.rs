//! `BRepBlend_SurfCurvConstRadInv` — solution of the constant-radius
//! surface/curve reframing on a given restriction of the surface.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_SurfCurvConstRadInv.hxx:33-90`,
//! `.cxx:17-284` (OCCT 8.0.0).
//!
//! The vector `X` is `(wguide, wcurv, wrst)`: the parameter on the guide line,
//! the parameter on the curve and the parameter on the restriction
//! (`hxx:43-48`).
//!
//! The three equations mirror [`super::surf_point_const_rad_inv`], with the
//! "point" replaced by a curve point `curv->Value(wcurv)` and the "surface
//! point" taken at the restriction `rst->Value(wrst)`:
//! 1. `nplan . ptcur + theD`,
//! 2. `nplan . pts + theD`,
//! 3. `|corde + ray * ns|^2 - ray^2`.
//!
//! Translation notes:
//! - `rst` is an `occ::handle<Adaptor2d_Curve2d>`; the port stores it as
//!   `Option<Arc<dyn Curve2d>>` because [`BlendSurfCurvFuncInv::set`] admits
//!   the empty handle. OCCT dereferences it unguarded; the port maps the empty
//!   handle to "not computed" (`false`) instead of undefined behaviour. The
//!   walking builders always set it before a solve.
//! - `D(3, 2) = -ref.Dot(d1cur)` (`cxx:216`, `cxx:257`) carries the minus sign
//!   of `corde = gp_Vec(ptcur, pts)`; the port keeps it.
//! - `dwrstpts` / `dwrstnsurf` are `SetLinearForm(x, d1u, y, d1v)` on the
//!   2d tangent of the restriction (`cxx:190-191`).

use std::sync::Arc;

use occt_core::gp::GpVec;
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_geom2d::curve::Curve2d;
use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use crate::blend::BlendSurfCurvFuncInv;

use super::pnt_vec;

/// `class BRepBlend_SurfCurvConstRadInv : public Blend_SurfCurvFuncInv`
/// (`hxx:33-90`).
pub struct SurfCurvConstRadInv {
    /// `surf` (`hxx:79`).
    surf: Arc<dyn Surface>,
    /// `curv` (`hxx:80`).
    curv: Arc<dyn Curve>,
    /// `guide` (`hxx:81`).
    guide: Arc<dyn Curve>,
    /// `rst` (`hxx:82`): `occ::handle<Adaptor2d_Curve2d>`.
    rst: Option<Arc<dyn Curve2d>>,
    /// `ray` (`hxx:83`).
    ray: f64,
    /// `choix` (`hxx:84`).
    choix: i32,
}

impl SurfCurvConstRadInv {
    /// `BRepBlend_SurfCurvConstRadInv(S, C, Cg)` (`cxx:21-32`): `rst` stays
    /// null, `ray = 0.`, `choix = 0`.
    pub fn new(surf: Arc<dyn Surface>, curv: Arc<dyn Curve>, guide: Arc<dyn Curve>) -> Self {
        Self {
            surf,
            curv,
            guide,
            rst: None,
            ray: 0.0,
            choix: 0,
        }
    }

    /// `Set(R, Choix)` (`cxx:36-54`).
    pub fn set_radius(&mut self, r: f64, choix: i32) {
        self.choix = choix;
        match choix {
            1 | 2 => self.ray = -r.abs(),
            3 | 4 => self.ray = r.abs(),
            _ => self.ray = -r.abs(),
        }
    }

    /// `ray` (`hxx:83`), read by `IsSolution`.
    pub fn ray(&self) -> f64 {
        self.ray
    }
}

impl MathFunctionSetWithDerivatives for SurfCurvConstRadInv {
    /// `Blend_SurfCurvFuncInv::NbVariables()` = `3`
    /// (`Blend_SurfCurvFuncInv.cxx:22-25`).
    fn nb_variables(&self) -> usize {
        3
    }

    /// `NbEquations()` (`cxx:58-61`): `3`.
    fn nb_equations(&self) -> usize {
        3
    }

    /// `Value(X, F)` (`cxx:63-88`). Empty `rst` gives `false`.
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        let Some(rst) = self.rst.clone() else {
            return false;
        };
        let (ptgui, d1gui) = self.guide.d1(x.value(1));
        let nplan = d1gui.normalized();
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

    /// `Derivatives(X, D)` (`cxx:92-184`).
    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        let Some(rst) = self.rst.clone() else {
            return false;
        };
        let (ptgui, d1gui, d2gui) = self.guide.d2(x.value(1));
        let normd1gui = d1gui.magnitude();
        let unsurnormd1gui = 1.0 / normd1gui;
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
        let dwguiref =
            GpVec::zero().set_linear_form_2(raysurnorm, &dwguitemp, -raysurnorm2 * dwguinorm, &temp);
        let dwrstref = GpVec::zero().set_linear_form_2plus(
            raysurnorm,
            &dwrsttemp,
            -raysurnorm2 * dwrstnorm,
            &temp,
            &dwrstpts,
        );

        r = r.added(&r);
        d.set_value(3, 1, r.dot(&dwguiref));
        d.set_value(3, 2, -r.dot(&d1cur));
        d.set_value(3, 3, r.dot(&dwrstref));

        true
    }

    /// `Values(X, F, D)` (`cxx:188-270`): `Value` and `Derivatives` inlined.
    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        let Some(rst) = self.rst.clone() else {
            return false;
        };
        let (ptgui, d1gui, d2gui) = self.guide.d2(x.value(1));
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
        let dwguiref =
            GpVec::zero().set_linear_form_2(raysurnorm, &dwguitemp, -raysurnorm2 * dwguinorm, &temp);
        let dwrstref = GpVec::zero().set_linear_form_2plus(
            raysurnorm,
            &dwrsttemp,
            -raysurnorm2 * dwrstnorm,
            &temp,
            &dwrstpts,
        );

        r = r.added(&r);
        d.set_value(3, 1, r.dot(&dwguiref));
        d.set_value(3, 2, -r.dot(&d1cur));
        d.set_value(3, 3, r.dot(&dwrstref));
        true
    }
}

impl BlendSurfCurvFuncInv for SurfCurvConstRadInv {
    /// `Set(Rst)` (`cxx:274-277`).
    fn set(&mut self, rst: Option<Arc<dyn Curve2d>>) {
        self.rst = rst;
    }

    /// `GetTolerance(Tolerance, Tol)` (`cxx:279-288`): `Tolerance(3)` is the
    /// restriction resolution at the smaller of the surface's two parametric
    /// resolutions.
    fn get_tolerance(&self, tolerance: &mut MathVector, tol: f64) {
        tolerance.set_value(1, self.guide.resolution(tol));
        tolerance.set_value(2, self.curv.resolution(tol));
        let ru = occt_geom::approx_same_parameter::u_resolution(self.surf.as_ref(), tol);
        let rv = occt_geom::approx_same_parameter::v_resolution(self.surf.as_ref(), tol);
        let r = match &self.rst {
            // `rst->Resolution(std::min(ru, rv))`: an empty handle has no
            // resolution; OCCT never reaches this with a null `rst`.
            Some(rst) => rst.resolution(ru.min(rv)),
            None => 0.0,
        };
        tolerance.set_value(3, r);
    }

    /// `GetBounds(InfBound, SupBound)` (`cxx:290-302`).
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

    /// `IsSolution(Sol, Tol)` (`cxx:304-311`).
    fn is_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        let mut valsol = MathVector::new(1, 3);
        self.value(sol, &mut valsol);
        valsol.value(1).abs() <= tol
            && valsol.value(2).abs() <= tol
            && valsol.value(3).abs() <= 2.0 * tol * self.ray.abs()
    }
}
