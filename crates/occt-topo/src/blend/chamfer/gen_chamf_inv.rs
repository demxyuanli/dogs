//! `BlendFunc_GenChamfInv` — the deferred base of the inverse chamfer solvers:
//! a solution on a restriction (curve on surface) of one of the two supports.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_GenChamfInv.hxx:31-107`,
//! `BlendFunc_GenChamfInv.cxx:23-122` (OCCT 8.0.0).
//!
//! The vector `X` is `(t, w, U, V)`: the parameter on the curve on surface, the
//! parameter on the guide line, and the parameters of the partner point on the
//! other surface.
//!
//! Translation notes:
//! - `csurf` is an `occ::handle<Adaptor2d_Curve2d>`; the port stores it as
//!   `Option<Arc<dyn Curve2d>>` because `Set(OnFirst, C)` can be handed the
//!   empty handle. OCCT dereferences it unguarded in `Value` / `Derivatives` /
//!   `Values`; the port answers `false` instead of faulting. `GetTolerance`
//!   keeps the null case at `0.` and `GetBounds` at `(0., 0.)`.
//! - `Values` (`cxx:119-122`) ignores the boolean `Value` answers and returns
//!   `true`; the port keeps that.

use std::sync::Arc;

use occt_core::precision::Precision;
use occt_geom::approx_same_parameter;
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_geom2d::curve::Curve2d;
use occt_math::MathVector;

/// `class BlendFunc_GenChamfInv : public Blend_FuncInv` (`hxx:31-107`).
pub struct GenChamfInv {
    /// `surf1` (`hxx`): the first support.
    pub(super) surf1: Arc<dyn Surface>,
    /// `surf2` (`hxx`): the second support.
    pub(super) surf2: Arc<dyn Surface>,
    /// `curv` (`hxx`): the guide line.
    pub(super) curv: Arc<dyn Curve>,
    /// `choix` (`hxx`): the "quadrant" of the chamfer.
    pub(super) choix: i32,
    /// `first` (`hxx`): `true` when `csurf` lies on `surf1`.
    pub(super) first: bool,
    /// `csurf` (`hxx`): the curve on surface; `Set(OnFirst, C)` writes it.
    pub(super) csurf: Option<Arc<dyn Curve2d>>,
}

impl GenChamfInv {
    /// `BlendFunc_GenChamfInv(S1, S2, C)` (`cxx:23-32`): `csurf` stays null,
    /// `choix = 0`, `first = false`.
    pub fn new(surf1: Arc<dyn Surface>, surf2: Arc<dyn Surface>, curv: Arc<dyn Curve>) -> Self {
        Self {
            surf1,
            surf2,
            curv,
            choix: 0,
            first: false,
            csurf: None,
        }
    }

    /// `Set(OnFirst, C)` (`cxx:45-49`).
    pub(super) fn set_first(&mut self, on_first: bool, c_on_surf: Option<Arc<dyn Curve2d>>) {
        self.first = on_first;
        self.csurf = c_on_surf;
    }

    /// `GetTolerance(Tolerance, Tol)` (`cxx:52-65`).
    pub(super) fn fill_tolerance(&self, tolerance: &mut MathVector, tol: f64) {
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

    /// `GetBounds(InfBound, SupBound)` (`cxx:69-114`): the box of `csurf` and
    /// the guide line, plus the box of the partner surface widened by its own
    /// range when both ends are finite.
    pub(super) fn fill_bounds(&self, inf_bound: &mut MathVector, sup_bound: &mut MathVector) {
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

        if !Precision::is_infinite(inf_bound.value(3)) && !Precision::is_infinite(sup_bound.value(3))
        {
            let range = sup_bound.value(3) - inf_bound.value(3);
            inf_bound.set_value(3, inf_bound.value(3) - range);
            sup_bound.set_value(3, sup_bound.value(3) + range);
        }
        if !Precision::is_infinite(inf_bound.value(4)) && !Precision::is_infinite(sup_bound.value(4))
        {
            let range = sup_bound.value(4) - inf_bound.value(4);
            inf_bound.set_value(4, inf_bound.value(4) - range);
            sup_bound.set_value(4, sup_bound.value(4) + range);
        }
    }

    /// The point `(U, V)` of the partner surface, in the order `Value` wants:
    /// `csurf` coordinates first when `first`, last otherwise (`cxx:100-117`).
    pub(super) fn partner_parameters(
        &self,
        p2d: &occt_core::gp::GpPnt2d,
        x: &MathVector,
    ) -> (f64, f64, f64, f64) {
        if self.first {
            (p2d.x(), p2d.y(), x.value(3), x.value(4))
        } else {
            (x.value(3), x.value(4), p2d.x(), p2d.y())
        }
    }
}
