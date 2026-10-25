//! `Adaptor3d_CurveOnSurface` — the 3D curve carried by a pcurve on a surface.
//!
//! Source: `ModelingData/TKG3d/Adaptor3d/Adaptor3d_CurveOnSurface.hxx`,
//! `Adaptor3d_CurveOnSurface.cxx` (OCCT 8.0.0).
//!
//! Ported: the evaluation core the blending functions need — `FirstParameter` /
//! `LastParameter` (`cxx:977-987`), `Continuity` (`cxx:991-1006`),
//! `EvalD0` (`cxx:1176-1197`), `EvalD1` (`cxx:1202-1239`), `EvalD2`
//! (`cxx:1244-1289`), `Resolution` (`cxx:1364-1371`), `IsPeriodic`
//! (`cxx:1152-1161`) and `Period` (`cxx:1164-1172`) for the generic (non
//! `GeomAbs_Line` / `GeomAbs_Circle`) path.
//!
//! `EvalKPart` (`cxx:1552-1735`) builds the type table that the `GeomAbs_Line`
//! and `GeomAbs_Circle` fast paths and the `GetType` / `Line` / `Circle` /
//! `Ellipse` / `Hyperbola` / `Parabola` queries read. Those fast paths evaluate
//! `ElCLib::*` on the exact 3D image of the pcurve, which is the same point and
//! the same derivatives as the chain rule below, so only the chain rule is
//! written here. The type table itself is UNPORTED, which is why
//! `is_periodic` / `period` fall back to the 2D curve instead of reporting
//! `2*pi` for a pcurve that `EvalKPart` would classify as a circle/ellipse.
//!
//! `myFirstSurf` / `myLastSurf` (`EvalFirstLastSurf`, `cxx:1736-1832`) switch the
//! evaluation to the surface of the adjacent face when the pcurve reaches the
//! parameter-space boundary of a B-spline / extrusion / revolution surface. That
//! is UNPORTED: this port always evaluates on `mySurface`.
//!
//! `NbIntervals` / `Intervals` (`cxx:1045-1130`) add the roots of the pcurve
//! crossing each surface discontinuity to the 2D curve's own break points
//! (`Adaptor3d_InterFunc` + `math_FunctionRoots` + `AddIntervals`). UNPORTED:
//! the port forwards to the 2D curve's break points only.
//!
//! `Trim` (`cxx:1133-1141`) rebuilds the adaptor from `myCurve->Trim(First,
//! Last, Tol)`; `Geom2d_Curve` has no `Trim` in this port, so it is UNPORTED
//! too. `Transform` / `Reverse` do not exist on OCCT's `Adaptor3d_Curve`; they
//! are no-ops here (`Transform`) and pcurve-forwarding (`Reverse`).

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpTrsf, GpVec};
use occt_geom2d::curve::Curve2d;

use crate::approx_same_parameter::{u_resolution, v_resolution};
use crate::curve::Curve;
use crate::surface::Surface;

/// `class Adaptor3d_CurveOnSurface : public Adaptor3d_Curve` (`hxx:29-...`).
#[derive(Clone)]
pub struct CurveOnSurface {
    /// `myCurve` (`hxx:...`): the pcurve.
    curve: Arc<dyn Curve2d>,
    /// `mySurface` (`hxx:...`): the surface the pcurve lives on.
    surface: Arc<dyn Surface>,
}

impl CurveOnSurface {
    /// `Adaptor3d_CurveOnSurface()` followed by
    /// `Load(C, S)` (`cxx:968-972`).
    pub fn new(curve: Arc<dyn Curve2d>, surface: Arc<dyn Surface>) -> Self {
        Self { curve, surface }
    }

    /// `GetCurve()` (`cxx:1524-1527`).
    pub fn curve(&self) -> &Arc<dyn Curve2d> {
        &self.curve
    }

    /// `GetSurface()` (`cxx:1531-1534`).
    pub fn surface(&self) -> &Arc<dyn Surface> {
        &self.surface
    }
}

impl Curve for CurveOnSurface {
    /// `EvalD0` (`cxx:1176-1197`), generic arm (`cxx:1189-1195`):
    /// `myCurve->D0(U, Puv); mySurface->D0(Puv.X(), Puv.Y(), P)`.
    fn d0(&self, u: f64) -> GpPnt {
        let uv = self.curve.d0(u);
        self.surface.d0(uv.x(), uv.y())
    }

    /// `EvalD1` (`cxx:1202-1239`), generic arm (`cxx:1233-1237`):
    /// `D1 = Duv.X() * D1U + Duv.Y() * D1V`.
    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        let (uv, duv) = self.curve.d1(u);
        let (p, d1u, d1v) = self.surface.d1(uv.x(), uv.y());
        (p, GpVec::zero().set_linear_form_2(duv.x(), &d1u, duv.y(), &d1v))
    }

    /// `EvalD2` (`cxx:1244-1289`), generic arm (`cxx:1282-1287`):
    /// `D2 = D2W.X() * D1U + D2W.Y() * D1V + 2 * DW.X() * DW.Y() * D2UV`,
    /// then `+= DW.X()^2 * D2U + DW.Y()^2 * D2V`.
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        let (uv, dw, d2w) = self.curve.d2(u);
        let (p, d1u, d1v, d2u, d2v, d2uv) = self.surface.d2(uv.x(), uv.y());
        let d1 = GpVec::zero().set_linear_form_2(dw.x(), &d1u, dw.y(), &d1v);
        // `aRes.D2.SetLinearForm(D2W.X(), D1U, D2W.Y(), D1V, 2*DW.X*DW.Y, D2UV)`
        // (`cxx:1285`), then
        // `aRes.D2.SetLinearForm(DW.X()*DW.X(), D2U, DW.Y()*DW.Y(), D2V, aRes.D2)`
        // (`cxx:1286`): the previous value enters as the trailing vector.
        let d2 = GpVec::zero().set_linear_form_3(
            d2w.x(),
            &d1u,
            d2w.y(),
            &d1v,
            2.0 * dw.x() * dw.y(),
            &d2uv,
        );
        let d2 =
            GpVec::zero().set_linear_form_2plus(dw.x() * dw.x(), &d2u, dw.y() * dw.y(), &d2v, &d2);
        (p, d1, d2)
    }

    /// `FirstParameter()` (`cxx:977-980`): `myCurve->FirstParameter()`.
    fn first_parameter(&self) -> f64 {
        self.curve.first_parameter()
    }

    /// `LastParameter()` (`cxx:984-987`): `myCurve->LastParameter()`.
    fn last_parameter(&self) -> f64 {
        self.curve.last_parameter()
    }

    /// `Continuity()` (`cxx:991-1006`): the minimum of the pcurve's continuity
    /// and the surface's `UContinuity` / `VContinuity`. The port's
    /// [`Surface::continuity`] carries the single value OCCT's two queries
    /// share for every concrete surface, so the minimum is taken over the pair.
    fn continuity(&self) -> u8 {
        self.curve.continuity().min(self.surface.continuity())
    }

    /// `Resolution(R3d)` (`cxx:1364-1371`):
    /// `myCurve->Resolution(min(surf.UResolution(R3d), surf.VResolution(R3d)))`.
    fn resolution(&self, r3d: f64) -> f64 {
        let ru = u_resolution(self.surface.as_ref(), r3d);
        let rv = v_resolution(self.surface.as_ref(), r3d);
        self.curve.resolution(ru.min(rv))
    }

    /// `IsPeriodic()` (`cxx:1152-1161`). UNPORTED: the `GeomAbs_Circle` /
    /// `GeomAbs_Ellipse` arm of the type table, see the module note.
    fn is_periodic(&self) -> bool {
        self.curve.is_periodic()
    }

    /// `Period()` (`cxx:1164-1172`). UNPORTED: the `GeomAbs_Circle` /
    /// `GeomAbs_Ellipse` arm of the type table, see the module note.
    fn period(&self) -> f64 {
        self.curve.period()
    }

    /// `Trim` (`cxx:1133-1141`) is UNPORTED (see the module note), and
    /// `Adaptor3d_Curve` has no `Transform`; this is a no-op.
    fn transform(&mut self, _t: &GpTrsf) {}

    /// `Adaptor3d_Curve` has no `Reverse`; OCCT offers no operation to forward
    /// to, so this is a no-op (a shared `Arc<dyn Curve2d>` cannot be reversed in
    /// place either).
    fn reverse(&mut self) {}

    fn clone_dyn(&self) -> Box<dyn Curve> {
        Box::new(self.clone())
    }

    /// `NbIntervals(S)` (`cxx:1045-1111`) is UNPORTED (see the module note); the
    /// 2D curve's own break points are the only ones reported.
    fn nb_intervals(&self, s: u8) -> i32 {
        self.curve.nb_intervals(s)
    }

    /// `Intervals(T, S)` (`cxx:1119-1130`) is UNPORTED (see the module note).
    fn parameter_intervals(&self, s: u8) -> Vec<f64> {
        self.curve.parameter_intervals(s)
    }
}
