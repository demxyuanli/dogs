//! `BlendFunc_GenChamfer` — the deferred base of the chamfer solvers between
//! two surfaces, driven by a guide line.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_GenChamfer.hxx:31-149`,
//! `BlendFunc_GenChamfer.cxx:29-300` (OCCT 8.0.0).
//!
//! The vector `X` is `(U1, V1, U2, V2)`: the parameters of the two points, one
//! on each support. The four equations themselves are written by the concrete
//! classes (`Value` / `Derivatives` stay pure here); this class supplies the
//! section shape (a degree 1 segment, `2` poles), the tolerance / bound
//! queries, the analytic section used by the visualization code and the two
//! `Section` bodies (see [`super::shared`]).
//!
//! Translation notes:
//! - `Set(First, Last)` (`cxx:50`) is empty; the port's `set_interval` of
//!   every concrete class is the same no-op.
//! - `Section(Param, U1, V1, U2, V2, Pdeb, Pfin, C)` (`cxx:104-139`) is not a
//!   virtual of `Blend_Function`; it stays an inherent method here and answers
//!   `(Pdeb, Pfin, C)`.
//! - `GetSectionSize` is pure at this level; the concrete classes throw
//!   `Standard_NotImplemented` (`ConstThroat.cxx:276-280`).

use std::sync::Arc;

use occt_core::elib::clib;
use occt_core::gp::{GpDir, GpLin, GpVec};
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::Precision;
use occt_geom::approx_same_parameter;
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_math::MathVector;

use crate::blend::next_shape;

/// `class BlendFunc_GenChamfer : public Blend_Function` (`hxx:31-149`).
pub struct GenChamfer {
    /// `surf1` (`hxx`): the first support.
    pub(super) surf1: Arc<dyn Surface>,
    /// `surf2` (`hxx`): the second support.
    pub(super) surf2: Arc<dyn Surface>,
    /// `curv` (`hxx`): the guide line.
    pub(super) curv: Arc<dyn Curve>,
    /// `choix` (`hxx`): the "quadrant" of the chamfer.
    pub(super) choix: i32,
    /// `tol` (`hxx`): the 3d tolerance last given to `IsSolution`.
    pub(super) tol: f64,
    /// `distmin` (`hxx`): the smallest distance between the two extremities
    /// seen so far, seeded to `RealLast()` (`cxx:38`).
    pub(super) distmin: f64,
}

impl GenChamfer {
    /// `BlendFunc_GenChamfer(S1, S2, CG)` (`cxx:29-40`).
    pub fn new(surf1: Arc<dyn Surface>, surf2: Arc<dyn Surface>, curv: Arc<dyn Curve>) -> Self {
        Self {
            surf1,
            surf2,
            curv,
            choix: 0,
            tol: 0.0,
            distmin: f64::MAX,
        }
    }

    /// `GetTolerance(Tolerance, Tol)` (`cxx:54-61`): the parametric resolution
    /// of each support.
    pub(super) fn fill_tolerance(&self, tolerance: &mut MathVector, tol: f64) {
        tolerance.set_value(1, approx_same_parameter::u_resolution(self.surf1.as_ref(), tol));
        tolerance.set_value(2, approx_same_parameter::v_resolution(self.surf1.as_ref(), tol));
        tolerance.set_value(3, approx_same_parameter::u_resolution(self.surf2.as_ref(), tol));
        tolerance.set_value(4, approx_same_parameter::v_resolution(self.surf2.as_ref(), tol));
    }

    /// `GetBounds(InfBound, SupBound)` (`cxx:65-84`): the box of both supports,
    /// each finite range widened by its own length.
    pub(super) fn fill_bounds(&self, inf_bound: &mut MathVector, sup_bound: &mut MathVector) {
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

    /// `GetMinimalDistance()` (`cxx:88-91`).
    pub(super) fn minimal_distance(&self) -> f64 {
        self.distmin
    }

    /// `Section(Param, U1, V1, U2, V2, Pdeb, Pfin, C)` (`cxx:104-139`): the
    /// line from the point on the first support to the point on the second,
    /// parameterised from `0` to the parameter of the second point.
    ///
    /// The body only reads the two surfaces, so it is written once in
    /// [`gen_chamfer_section`] and reached by name below; `Param` is unused in
    /// OCCT too (`cxx:105`).
    pub fn analytic_section(
        &self,
        _param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) -> (f64, f64, GpLin) {
        gen_chamfer_section(&*self.surf1, &*self.surf2, u1, v1, u2, v2)
            .expect("BlendFunc_GenChamfer::Section: degenerate segment")
    }

    /// `IsRational()` (`cxx:125-128`): `false`.
    pub(super) fn is_rational(&self) -> bool {
        false
    }

    /// `GetMinimalWeight(Weights)` (`cxx:132-135`): all weights are `1`.
    pub(super) fn fill_minimal_weights(&self, weights: &mut [f64]) {
        weights.fill(1.0);
    }

    /// `NbIntervals(S)` (`cxx:139-142`): `curv->NbIntervals(NextShape(S))`.
    pub(super) fn nb_intervals(&self, s: Shape) -> i32 {
        self.curv.nb_intervals(next_shape(s) as u8)
    }

    /// `Intervals(T, S)` (`cxx:146-149`): `curv->Intervals(T, NextShape(S))`.
    pub(super) fn guide_intervals(&self, t: &mut [f64], s: Shape) {
        let src = self.curv.parameter_intervals(next_shape(s) as u8);
        let count = t.len().min(src.len());
        t[..count].copy_from_slice(&src[..count]);
    }

    /// `GetShape(NbPoles, NbKnots, Degree, NbPoles2d)` (`cxx:153-159`): a
    /// degree 1 segment, `2` poles and `2` poles 2d.
    pub(super) fn section_shape(
        &self,
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

    /// `GetTolerance(BoundTol, SurfTol, AngleTol, Tol3d, Tol1d)` (`cxx:165-172`):
    /// every 3d tolerance is `BoundTol`; `Tol1d` is left untouched.
    pub(super) fn section_tolerance_3d(&self, bound_tol: f64, tol3d: &mut MathVector) {
        tol3d.init(bound_tol);
    }

    /// `Knots(TKnots)` (`cxx:176-180`).
    pub(super) fn section_knots(&self, tknots: &mut [f64]) {
        tknots[0] = 0.0;
        tknots[1] = 1.0;
    }

    /// `Mults(TMults)` (`cxx:184-188`).
    pub(super) fn section_mults(&self, tmults: &mut [i32]) {
        tmults[0] = 2;
        tmults[1] = 2;
    }

    /// `Resolution(IC2d, Tol, TolU, TolV)` (`cxx:286-300`): the resolution of
    /// the first support on `IC2d == 1`, of the second one otherwise.
    pub(super) fn resolution(&self, ic2d: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
        if ic2d == 1 {
            *tol_u = approx_same_parameter::u_resolution(self.surf1.as_ref(), tol);
            *tol_v = approx_same_parameter::v_resolution(self.surf1.as_ref(), tol);
        } else {
            *tol_u = approx_same_parameter::u_resolution(self.surf2.as_ref(), tol);
            *tol_v = approx_same_parameter::v_resolution(self.surf2.as_ref(), tol);
        }
    }
}

/// `Section(Param, U1, V1, U2, V2, Pdeb, Pfin, C)` (`cxx:104-122`).
///
/// The line through the point of the first support and the point of the
/// second, with `Pdeb = 0` and `Pfin` the line parameter of the second point.
/// `Param` is commented out in OCCT's signature (`cxx:105`) and is not taken.
///
/// `BlendFunc_ChAsym::Section` (`BlendFunc_ChAsym.cxx:484-503`) has the same
/// body, so the `ChFi3d_ChBuilder` call sites reach this one function for both
/// families.
///
/// `gp_Dir(gp_Vec(Pt1, Pt2))` raises `Standard_ConstructionError` on a
/// degenerate segment (`cxx:113`); the port answers `Err`.
pub(crate) fn gen_chamfer_section(
    surf1: &dyn Surface,
    surf2: &dyn Surface,
    u1: f64,
    v1: f64,
    u2: f64,
    v2: f64,
) -> Result<(f64, f64, GpLin), String> {
    let pts1 = surf1.value(u1, v1);
    let pts2 = surf2.value(u2, v2);
    let dir = GpDir::from_vec(&GpVec::from_pnts(&pts1, &pts2))
        .map_err(|e| format!("BlendFunc_GenChamfer::Section: {e}"))?;
    let lin = GpLin::from_pnt_dir(pts1, dir);
    let pdeb = 0.0;
    let pfin = clib::line_parameter(&lin, &pts2);
    Ok((pdeb, pfin, lin))
}
