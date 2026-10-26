//! `GeomConvert_ApproxSurface` (`GeomConvert_ApproxSurface.cxx:312-460`).
//!
//! The surface is approximated by `AdvApp2Var_ApproxAFunc2Var` with one 3D
//! subspace (`nb1 = nb2 = 0`, `nb3 = 1`), `FavorIso = GeomAbs_IsoV`, the
//! requested tolerances and cuts, and the preferential C2/C3 breaks of
//! `AdvApprox_PrefAndRec` (`cxx:381-398`).
//!
//! The `GeomConvert_ApproxSurface_Eval` callback (`cxx:55-310`) is
//! `adv_app2var::eval_surface`.

use occt_core::kernel::geomabs::Shape;

use crate::adv_app2var::{approximate_until_cdi, context_until_roots};
use crate::bspline_surface::GeomBSplineSurface;
use crate::surface::Surface;

/// `GeomConvert_ApproxSurface` (`GeomConvert_ApproxSurface.hxx:32-108`).
pub struct GeomConvertApproxSurface {
    /// `myIsDone` (`hxx:106`).
    is_done: bool,
    /// `myHasResult` (`hxx:107`).
    has_result: bool,
    /// `myBSplSurf` (`hxx:108`).
    bspl_surf: Option<GeomBSplineSurface>,
    /// `myMaxError` (`hxx:109`).
    max_error: f64,
}

impl GeomConvertApproxSurface {
    /// `GeomConvert_ApproxSurface(Adaptor3d_Surface, ...)`
    /// (`GeomConvert_ApproxSurface.cxx:325-343`). The `Geom_Surface`
    /// constructor (`:312-331`) first wraps the surface into a
    /// `GeomAdaptor_Surface`, which is already the `Surface` trait object here.
    ///
    /// The C++ constructor throws through `AdvApp2Var_ApproxAFunc2Var::Init`
    /// (`AdvApp2Var_ApproxAFunc2Var.cxx:171-205`, `:209-215`) when the
    /// continuity is not C0/C1/C2 or the maximal degree is below
    /// `2 * order + 1`. Those cases leave `myHasResult = false` here instead
    /// of unwinding, because the single caller only passes C1.
    pub fn new(
        surf: &dyn Surface,
        tol3d: f64,
        u_continuity: Shape,
        v_continuity: Shape,
        max_deg_u: i32,
        max_deg_v: i32,
        max_segments: i32,
        precis_code: i32,
    ) -> Self {
        let mut this = Self {
            is_done: false,
            has_result: false,
            bspl_surf: None,
            max_error: 0.0,
        };
        this.approximate(
            surf,
            tol3d,
            u_continuity,
            v_continuity,
            max_deg_u,
            max_deg_v,
            max_segments,
            precis_code,
        );
        this
    }

    /// `GeomConvert_ApproxSurface::Approximate` (`cxx:345-430`).
    fn approximate(
        &mut self,
        surf: &dyn Surface,
        tol3d: f64,
        u_continuity: Shape,
        v_continuity: Shape,
        max_deg_u: i32,
        max_deg_v: i32,
        max_segments: i32,
        precis_code: i32,
    ) {
        let (u0, u1) = surf.u_range();
        let (v0, v1) = surf.v_range();

        // `cxx:361-376`: one 3D subspace; the 1D / 2D tolerances are the
        // all-zero arrays `nul1` / `nul2` and are dropped with
        // `nb1 = nb2 = 0`.
        let three_d_tol = [tol3d];
        let three_d_tol_front = [[tol3d; 4]];

        // `cxx:378`: `GeomAbs_IsoV` is `ifav = 2` (`AdvApp2Var_ApproxAFunc2Var.cxx:173`).
        let favor_iso = 2;
        let Ok(ord_u) = continuity_order(u_continuity) else {
            return;
        };
        let Ok(ord_v) = continuity_order(v_continuity) else {
            return;
        };
        // `AdvApp2Var_ApproxAFunc2Var.cxx:207-215`.
        let nd_u = (max_deg_u + 1).max(2 * ord_u + 2);
        let nd_v = (max_deg_v + 1).max(2 * ord_v + 2);
        if nd_u < 2 * ord_u + 2 || nd_v < 2 * ord_v + 2 {
            return;
        }
        // `AdvApp2Var_ApproxAFunc2Var.cxx:217`.
        let precis_code = precis_code.clamp(0, 3);

        // `cxx:400-428`: the surfaces are cut by the first `MaxPatch` squares.
        // The C2 / C3 breaks of `AdvApprox_PrefAndRec` (`cxx:381-398`) are
        // rebuilt by `approximate_until_cdi` from the same surface
        // (`adv_app2var/iso_fnc.rs:39-41`), so they are not duplicated here.
        let Ok(ctx) = context_until_roots(
            favor_iso,
            ord_u,
            ord_v,
            nd_u,
            nd_v,
            precis_code,
            0,
            0,
            1,
            &[],
            &[],
            &three_d_tol,
            &[],
            &[],
            &three_d_tol_front,
        ) else {
            return;
        };
        let Ok(result) = approximate_until_cdi(surf, u0, u1, v0, v1, &ctx, max_segments) else {
            return;
        };
        self.max_error = result.max_error;
        self.is_done = result.is_done;
        self.has_result = true;
        self.bspl_surf = Some(result.surface);
    }

    /// `Surface()` (`cxx:434-437`).
    pub fn surface(&self) -> Option<&GeomBSplineSurface> {
        self.bspl_surf.as_ref()
    }

    /// `Surface()` by value (`cxx:434-437`).
    pub fn into_surface(self) -> Option<GeomBSplineSurface> {
        self.bspl_surf
    }

    /// `IsDone()` (`cxx:441-444`).
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// `HasResult()` (`cxx:447-451`).
    pub fn has_result(&self) -> bool {
        self.has_result
    }

    /// `MaxError()` (`cxx:455-458`).
    pub fn max_error(&self) -> f64 {
        self.max_error
    }
}

/// `AdvApp2Var_ApproxAFunc2Var::Init` continuity switch
/// (`AdvApp2Var_ApproxAFunc2Var.cxx:178-205`).
fn continuity_order(continuity: Shape) -> Result<i32, ()> {
    match continuity {
        Shape::C0 => Ok(0),
        Shape::C1 => Ok(1),
        Shape::C2 => Ok(2),
        _ => Err(()),
    }
}
