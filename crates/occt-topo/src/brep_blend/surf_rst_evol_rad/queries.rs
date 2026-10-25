//! The query half of `BRepBlend_SurfRstEvolRad`: tolerances, bounds, section
//! shape, weights, intervals and knots.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_SurfRstEvolRad.cxx:52-110`
//! (`FusionneIntervalles`), `:337-368` (`GetTolerance` / `GetBounds`),
//! `:495-497` (`GetMinimalDistance`), `:726-836` (`IsRational`, `GetSectionSize`,
//! `GetMinimalWeight`, `NbIntervals`, `Intervals`, `GetShape`,
//! `GetTolerance(BoundTol, ...)`, `Knots`, `Mults`) (OCCT 8.0.0).

use occt_core::kernel::geomabs::Shape;
use occt_core::precision::{Precision, PCONFUSION};
use occt_geom::curve::Curve;
use occt_geom::geom_fill::section_convert;
use occt_math::MathVector;

use crate::blend::{
    get_minimal_weights, get_shape as blend_get_shape, next_shape, BlendFuncSectionShape,
};

use super::SurfRstEvolRad;

/// `FusionneIntervalles(I1, I2, Seq)` (`cxx:52-110`): merge two ascending
/// interval tables, coalescing parameters closer than
/// `PConfusion() * 0.99` into their midpoint.
///
/// OCCT fills a `NCollection_Sequence` 1-based; the port appends to a `Vec`.
pub(super) fn fusionne_intervalles(i1: &[f64], i2: &[f64], seq: &mut Vec<f64>) {
    let epspar = PCONFUSION * 0.99;
    let mut ind1 = 0usize;
    let mut ind2 = 0usize;
    while ind1 < i1.len() && ind2 < i2.len() {
        let v1 = i1[ind1];
        let v2 = i2[ind2];
        if (v1 - v2).abs() <= epspar {
            // Here the elements of I1 and I2 fit (`cxx:74-81`).
            seq.push((v1 + v2) / 2.0);
            ind1 += 1;
            ind2 += 1;
        } else if v1 < v2 {
            seq.push(v1);
            ind1 += 1;
        } else {
            seq.push(v2);
            ind2 += 1;
        }
    }
    while ind1 < i1.len() {
        seq.push(i1[ind1]);
        ind1 += 1;
    }
    while ind2 < i2.len() {
        seq.push(i2[ind2]);
        ind2 += 1;
    }
}

impl SurfRstEvolRad {
    /// `GetTolerance(Tolerance, Tol)` (`cxx:337-343`).
    pub(super) fn fill_parametric_tolerance(&self, tolerance: &mut MathVector, tol: f64) {
        tolerance.set_value(
            1,
            occt_geom::approx_same_parameter::u_resolution(self.surf.as_ref(), tol),
        );
        tolerance.set_value(
            2,
            occt_geom::approx_same_parameter::v_resolution(self.surf.as_ref(), tol),
        );
        tolerance.set_value(3, self.cons.resolution(tol));
    }

    /// `GetBounds(InfBound, SupBound)` (`cxx:347-368`).
    pub(super) fn fill_bounds(&self, inf_bound: &mut MathVector, sup_bound: &mut MathVector) {
        let (u1, u2) = self.surf.u_range();
        let (v1, v2) = self.surf.v_range();
        inf_bound.set_value(1, u1);
        inf_bound.set_value(2, v1);
        inf_bound.set_value(3, self.cons.first_parameter());
        sup_bound.set_value(1, u2);
        sup_bound.set_value(2, v2);
        sup_bound.set_value(3, self.cons.last_parameter());

        if !Precision::is_infinite(inf_bound.value(1)) && !Precision::is_infinite(sup_bound.value(1))
        {
            let range = sup_bound.value(1) - inf_bound.value(1);
            inf_bound.set_value(1, inf_bound.value(1) - range);
            sup_bound.set_value(1, sup_bound.value(1) + range);
        }
        if !Precision::is_infinite(inf_bound.value(2)) && !Precision::is_infinite(sup_bound.value(2))
        {
            let range = sup_bound.value(2) - inf_bound.value(2);
            inf_bound.set_value(2, inf_bound.value(2) - range);
            sup_bound.set_value(2, sup_bound.value(2) + range);
        }
    }

    /// `IsRational()` (`cxx:728-731`).
    pub(super) fn shape_is_rational(&self) -> bool {
        self.my_sshape == BlendFuncSectionShape::Rational
            || self.my_sshape == BlendFuncSectionShape::QuasiAngular
    }

    /// `GetSectionSize()` (`cxx:735-738`): the length of the maximum section,
    /// using the radius the law gave at the last `Set` / `Section`.
    pub(super) fn section_size(&self) -> f64 {
        self.maxang * self.ray.abs()
    }

    /// `GetMinimalWeight(Weights)` (`cxx:742-747`). OCCT raises out of
    /// `GeomConvert` for an angle out of range (i.e. before any `IsSolution`);
    /// the port panics at the same place.
    pub(super) fn minimal_weights(&self, weights: &mut [f64]) {
        get_minimal_weights(
            self.my_sshape,
            self.my_tconv,
            self.minang,
            self.maxang,
            weights,
        )
        .expect("BRepBlend_SurfRstEvolRad::GetMinimalWeight: section angle out of range");
    }

    /// `NbIntervals(S)` (`cxx:751-769`): `guide->NbIntervals(NextShape(S))`
    /// merged with `fevol->NbIntervals(S)`; when the law has a single interval
    /// the guide's count is returned unchanged.
    pub(super) fn nb_intervals_impl(&self, s: Shape) -> i32 {
        let nb_curve = self.guide.nb_intervals(next_shape(s) as u8);
        let nb_law = self
            .fevol
            .nb_intervals(s)
            .expect("BRepBlend_SurfRstEvolRad::NbIntervals: law NbIntervals failed");

        if nb_law == 1 {
            return nb_curve;
        }

        let int_c = self.guide.parameter_intervals(next_shape(s) as u8);
        let int_l = self
            .fevol
            .intervals(s)
            .expect("BRepBlend_SurfRstEvolRad::NbIntervals: law Intervals failed");
        let mut inter = Vec::new();
        fusionne_intervalles(&int_c, &int_l, &mut inter);
        inter.len() as i32 - 1
    }

    /// `Intervals(T, S)` (`cxx:773-797`).
    pub(super) fn intervals_impl(&self, t: &mut [f64], s: Shape) {
        let nb_law = self
            .fevol
            .nb_intervals(s)
            .expect("BRepBlend_SurfRstEvolRad::Intervals: law NbIntervals failed");

        if nb_law == 1 {
            let src = self.guide.parameter_intervals(next_shape(s) as u8);
            let count = t.len().min(src.len());
            t[..count].copy_from_slice(&src[..count]);
        } else {
            let int_c = self.guide.parameter_intervals(next_shape(s) as u8);
            let int_l = self
                .fevol
                .intervals(s)
                .expect("BRepBlend_SurfRstEvolRad::Intervals: law Intervals failed");
            let mut inter = Vec::new();
            fusionne_intervalles(&int_c, &int_l, &mut inter);
            let count = t.len().min(inter.len());
            t[..count].copy_from_slice(&inter[..count]);
        }
    }

    /// `GetShape(NbPoles, NbKnots, Degree, NbPoles2d)` (`cxx:801-805`).
    pub(super) fn section_shape(
        &mut self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles2d: &mut i32,
    ) {
        *nb_poles2d = 2;
        blend_get_shape(
            self.my_sshape,
            self.maxang,
            nb_poles,
            nb_knots,
            degree,
            &mut self.my_tconv,
        );
    }

    /// `GetTolerance(BoundTol, SurfTol, AngleTol, Tol3d, Tol1d)`
    /// (`cxx:807-821`). `GeomFill::GetTolerance` raises for an angle out of
    /// range in OCCT; the port panics at the same place.
    pub(super) fn fill_section_tolerance(
        &self,
        bound_tol: f64,
        surf_tol: f64,
        angle_tol: f64,
        tol3d: &mut MathVector,
        tol1d: &mut MathVector,
    ) {
        let tol = section_convert::get_tolerance(
            self.my_tconv,
            self.minang,
            self.ray.abs(),
            angle_tol,
            surf_tol,
        )
        .expect("BRepBlend_SurfRstEvolRad::GetTolerance: section angle out of range");
        tol1d.init(surf_tol);
        tol3d.init(surf_tol);
        let low = tol3d.lower();
        let up = tol3d.upper();
        let inner = tol.min(surf_tol);
        tol3d.set_value(low + 1, inner);
        tol3d.set_value(up - 1, inner);
        let outer = tol.min(bound_tol);
        tol3d.set_value(low, outer);
        tol3d.set_value(up, outer);
    }

    /// `Knots(TKnots)` (`cxx:825-827`).
    pub(super) fn section_knots(&self, tknots: &mut [f64]) {
        section_convert::knots(self.my_tconv, tknots);
    }

    /// `Mults(TMults)` (`cxx:831-833`).
    pub(super) fn section_mults(&self, tmults: &mut [i32]) {
        section_convert::mults(self.my_tconv, tmults);
    }
}
