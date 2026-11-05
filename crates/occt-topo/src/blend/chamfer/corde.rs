//! `BlendFunc_Corde` — the point of a surface at a given distance from a guide
//! line, inside the section plane normal to that line.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_Corde.hxx:28-100`,
//! `BlendFunc_Corde.cxx:30-201` (OCCT 8.0.0).
//!
//! The vector `X` is `(U, V)` on the surface. The two equations are
//! 1. `nplan . pts + theD` — the point lies in the section plane,
//! 2. `|pts - ptgui|^2 - dis^2` — the point is at the distance `dis` from the
//!    guide-line point.
//!
//! This is a helper of [`super::Chamfer`] / [`super::ChamfInv`], not a
//! `Blend_Function`.
//!
//! Translation notes:
//! - `pt2d` (`hxx:89`) is declared but never written or read by any body of
//!   the `.cxx`; it stays as a `#[allow(dead_code)]` member.
//! - `TangentOnS` / `Tangent2dOnS` throw `Standard_DomainError` when
//!   `istangent` (`cxx:124-141`); the accessors return the stored value (every
//!   OCCT caller guards the call with `!IsTangencyPoint()`).

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_math::{Gauss, MathMatrix, MathVector};

use super::shared::{GAUSS_MIN_PIVOT, pnt_vec};

/// `class BlendFunc_Corde` (`hxx:28-100`).
pub struct Corde {
    /// `surf` (`hxx:84`).
    surf: Arc<dyn Surface>,
    /// `guide` (`hxx:85`).
    guide: Arc<dyn Curve>,
    /// `pts` (`hxx:86`).
    pts: GpPnt,
    /// `pt2d` (`hxx:89`): never used by the `.cxx`.
    #[allow(dead_code)]
    pt2d: GpPnt2d,
    /// `dis` (`hxx:90`).
    dis: f64,
    /// `normtg` (`hxx:91`).
    normtg: f64,
    /// `theD` (`hxx:92`).
    the_d: f64,
    /// `ptgui` (`hxx:93`).
    ptgui: GpPnt,
    /// `nplan` (`hxx:94`).
    nplan: GpVec,
    /// `d1gui` (`hxx:95`).
    d1gui: GpVec,
    /// `d2gui` (`hxx:96`).
    d2gui: GpVec,
    /// `tgs` (`hxx:97`).
    tgs: GpVec,
    /// `tg2d` (`hxx:98`).
    tg2d: GpVec2d,
    /// `istangent` (`hxx:99`).
    istangent: bool,
}

impl Corde {
    /// `BlendFunc_Corde(S, CG)` (`cxx:30-39`).
    pub fn new(surf: Arc<dyn Surface>, guide: Arc<dyn Curve>) -> Self {
        Self {
            surf,
            guide,
            pts: GpPnt::zero(),
            pt2d: GpPnt2d::new(0.0, 0.0),
            dis: 0.0,
            normtg: 0.0,
            the_d: 0.0,
            ptgui: GpPnt::zero(),
            nplan: GpVec::zero(),
            d1gui: GpVec::zero(),
            d2gui: GpVec::zero(),
            tgs: GpVec::zero(),
            tg2d: GpVec2d::new(0.0, 0.0),
            istangent: false,
        }
    }

    /// `SetDist(Dist)` (`cxx:43-46`).
    pub fn set_dist(&mut self, dist: f64) {
        self.dis = dist;
    }

    /// `SetParam(Param)` (`cxx:50-57`): the guide-line frame.
    pub fn set_param(&mut self, param: f64) {
        let (ptgui, d1gui, d2gui) = self.guide.d2(param);
        self.ptgui = ptgui;
        self.d1gui = d1gui;
        self.d2gui = d2gui;
        self.normtg = self.d1gui.magnitude();
        self.nplan = self.d1gui.normalized();
        self.the_d = -self.nplan.dot(&pnt_vec(&self.ptgui));
    }

    /// `Value(X, F)` (`cxx:63-73`): always `true`.
    pub fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        let (pts, _du, _dv) = self.surf.d1(x.value(1), x.value(2));
        self.pts = pts;

        f.set_value(1, self.nplan.dot(&pnt_vec(&self.pts)) + self.the_d);
        let vref = GpVec::from_pnts(&self.ptgui, &self.pts);
        f.set_value(2, vref.square_magnitude() - self.dis * self.dis);

        true
    }

    /// `Derivatives(X, D)` (`cxx:81-91`): always `true`.
    pub fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        let (pts, d1u, d1v) = self.surf.d1(x.value(1), x.value(2));
        self.pts = pts;

        let vref = GpVec::from_pnts(&self.ptgui, &self.pts);
        d.set_value(1, 1, self.nplan.dot(&d1u));
        d.set_value(1, 2, self.nplan.dot(&d1v));
        d.set_value(2, 1, 2.0 * vref.dot(&d1u));
        d.set_value(2, 2, 2.0 * vref.dot(&d1v));

        true
    }

    /// `PointOnS()` (`cxx:95-98`).
    pub fn point_on_s(&self) -> &GpPnt {
        &self.pts
    }

    /// `PointOnGuide()` (`cxx:102-105`).
    pub fn point_on_guide(&self) -> &GpPnt {
        &self.ptgui
    }

    /// `NPlan()` (`cxx:109-112`).
    pub fn n_plan(&self) -> &GpVec {
        &self.nplan
    }

    /// `IsTangencyPoint()` (`cxx:116-119`).
    pub fn is_tangency_point(&self) -> bool {
        self.istangent
    }

    /// `TangentOnS()` (`cxx:124-131`).
    pub fn tangent_on_s(&self) -> &GpVec {
        &self.tgs
    }

    /// `Tangent2dOnS()` (`cxx:135-142`).
    pub fn tangent2d_on_s(&self) -> &GpVec2d {
        &self.tg2d
    }

    /// `DerFguide(Sol, DerF)` (`cxx:146-157`): the derivative of the two
    /// equations with respect to the guide-line parameter.
    pub fn der_fguide(&mut self, sol: &MathVector, der_f: &mut GpVec2d) {
        let (pts, _du, _dv) = self.surf.d1(sol.value(1), sol.value(2));
        self.pts = pts;

        let dnplan = GpVec::zero().set_linear_form_2(
            1.0 / self.normtg,
            &self.d2gui,
            -1.0 / self.normtg * self.nplan.dot(&self.d2gui),
            &self.nplan,
        );

        let temp = GpVec::from_pnts(&self.ptgui, &self.pts);

        der_f.set_coord(
            dnplan.dot(&temp) - self.nplan.dot(&self.d1gui),
            -2.0 * self.d1gui.dot(&temp),
        );
    }

    /// `IsSolution(Sol, Tol)` (`cxx:161-200`): the two values are within `Tol`
    /// (`Tol^2` for the second), then the section tangent is solved for from
    /// the jacobian with `math_Gauss` (pivot `1.e-20`, no fallback).
    pub fn is_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        let mut secmember = MathVector::new(1, 2);
        let mut valsol = MathVector::new(1, 2);
        let mut gradsol = MathMatrix::new(1, 2, 1, 2);

        self.value(sol, &mut valsol);
        self.derivatives(sol, &mut gradsol);

        if valsol.value(1).abs() <= tol && valsol.value(2).abs() <= tol * tol {
            let (pts, d1u, d1v) = self.surf.d1(sol.value(1), sol.value(2));
            self.pts = pts;

            let dnplan = GpVec::zero().set_linear_form_2(
                1.0 / self.normtg,
                &self.d2gui,
                -1.0 / self.normtg * self.nplan.dot(&self.d2gui),
                &self.nplan,
            );

            let temp = GpVec::from_pnts(&self.ptgui, &self.pts);

            secmember.set_value(1, self.nplan.dot(&self.d1gui) - dnplan.dot(&temp));
            secmember.set_value(2, 2.0 * self.d1gui.dot(&temp));

            let resol = Gauss::new(&gradsol, GAUSS_MIN_PIVOT);
            if resol.is_done() {
                if let Some(x) = resol.solve(&secmember) {
                    secmember = x;
                }
                self.tgs = GpVec::zero().set_linear_form_2(
                    secmember.value(1),
                    &d1u,
                    secmember.value(2),
                    &d1v,
                );
                self.tg2d = GpVec2d::new(secmember.value(1), secmember.value(2));
                self.istangent = false;
            } else {
                self.istangent = true;
            }
            return true;
        }

        false
    }
}
