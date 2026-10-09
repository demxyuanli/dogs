//! `IntCurveSurface_TheExactHInter`.
//!
//! Source: `IntImp_IntCS.gxx` as instantiated by
//! `IntCurveSurface_TheExactHInter_0.cxx` (`TheFunction =
//! IntCurveSurface_TheCSFunctionOfHInter`, `ThePSurfaceTool =
//! Adaptor3d_HSurfaceTool`, `TheCurveTool = IntCurveSurface_TheHCurveTool`).
//! Given a start point `(u, v, w)` coming out of the polygon / polyhedron
//! interference, it runs `math_FunctionSetRoot` on the 3x3 system
//! `S(u,v) = C(w)` over the supplied parameter box and accepts the root only
//! when the squared residual `Function().Root()` is within `TolTangency^2`.

use occt_core::gp::GpPnt;
use occt_geom::Surface;
use occt_math::{MathFunctionSetRoot, MathVector};

use super::cs_function::TheCSFunctionOfHInter;

/// `THE_TOLTANGENCY` (`IntCurveSurface_Inter.pxx:12`, `InterUtils.pxx`).
pub const THE_TOLTANGENCY: f64 = 1e-8;

/// `Adaptor3d_HSurfaceTool::UResolution/VResolution(S, Tol)`. The trait
/// `Surface::uv_resolution` returns the exact values where the port has them
/// (B-spline / offset / trimmed surfaces); the analytic family leaves it
/// `None`, in which case the linear estimate `Tol / |D1|` used by the
/// adaptor's `Geom_Surface::UResolution` is applied.
fn surface_uv_resolution(surface: &dyn Surface, u: f64, v: f64, tol: f64) -> (f64, f64) {
    if let Some((ru, rv)) = surface.uv_resolution(tol) {
        return (ru, rv);
    }
    let (_, du, dv) = surface.d1(u, v);
    let ur = if du.magnitude() > 1.0e-12 {
        tol / du.magnitude()
    } else {
        tol
    };
    let vr = if dv.magnitude() > 1.0e-12 {
        tol / dv.magnitude()
    } else {
        tol
    };
    (ur, vr)
}

/// `IntCurveSurface_TheExactHInter` (`...hxx:28-107`).
pub struct TheExactHInter<'a> {
    done: bool,
    empty: bool,
    my_function: TheCSFunctionOfHInter<'a>,
    w: f64,
    u: f64,
    v: f64,
    tol: f64,
}

impl<'a> TheExactHInter<'a> {
    /// `IntImp_IntCS(F, TolTangency)` (`...gxx:80-90`).
    pub fn new(function: TheCSFunctionOfHInter<'a>, tol_tangency: f64) -> Self {
        let mut tol = tol_tangency * tol_tangency;
        if tol < occt_core::precision::SQUARE_CONFUSION {
            tol = occt_core::precision::SQUARE_CONFUSION;
        }
        Self {
            done: true,
            empty: true,
            my_function: function,
            w: 0.0,
            u: 0.0,
            v: 0.0,
            tol,
        }
    }

    /// `IntImp_IntCS(U, V, W, F, TolTangency, MarginCoef = 0.0)`
    /// (`...gxx:24-79`).
    pub fn with_start(
        u: f64,
        v: f64,
        w: f64,
        function: TheCSFunctionOfHInter<'a>,
        tol_tangency: f64,
        margin_coef: f64,
    ) -> Self {
        let mut me = Self::new(function, tol_tangency);
        let s = me.my_function.auxillar_surface();
        let c = me.my_function.auxillar_curve();

        let w0 = c.first_parameter();
        let w1 = c.last_parameter();

        let (mut u0, mut u1) = s.u_range();
        let (mut v0, mut v1) = s.v_range();

        if margin_coef > 0.0 {
            if !occt_core::precision::Precision::is_infinite(u0) && !occt_core::precision::Precision::is_infinite(u1) {
                let mut marg = (u1 - u0) * margin_coef;
                if u0 > u1 {
                    marg = -marg;
                }
                u0 -= marg;
                u1 += marg;
            }
            if !occt_core::precision::Precision::is_infinite(v0) && !occt_core::precision::Precision::is_infinite(v1) {
                let mut marg = (v1 - v0) * margin_coef;
                if v0 > v1 {
                    marg = -marg;
                }
                v0 -= marg;
                v1 += marg;
            }
        }

        let mut rsnld = MathFunctionSetRoot::with_iterations(&me.my_function, 100);
        me.perform(u, v, w, &mut rsnld, u0, u1, v0, v1, w0, w1);
        me
    }

    /// `IntImp_IntCS::Perform` (`...gxx:92-152`). The parameter order follows
    /// the definition (`u0, u1, v0, v1, w0, w1`), not the header comment.
    #[allow(clippy::too_many_arguments)]
    pub fn perform(
        &mut self,
        u: f64,
        v: f64,
        w: f64,
        rsnld: &mut MathFunctionSetRoot,
        u0: f64,
        u1: f64,
        v0: f64,
        v1: f64,
        w0: f64,
        w1: f64,
    ) {
        self.done = true;
        let mut uvap = MathVector::new(1, 3);
        uvap.set_value(1, u);
        uvap.set_value(2, v);
        uvap.set_value(3, w);

        let mut born_inf = MathVector::new(1, 3);
        let mut born_sup = MathVector::new(1, 3);
        born_inf.set_value(1, u0);
        born_inf.set_value(2, v0);
        born_sup.set_value(1, u1);
        born_sup.set_value(2, v1);
        born_inf.set_value(3, w0);
        born_sup.set_value(3, w1);

        let s = self.my_function.auxillar_surface();
        let c = self.my_function.auxillar_curve();
        let (ru, rv) = surface_uv_resolution(s, u, v, occt_core::precision::CONFUSION);
        let mut tolerance = MathVector::new(1, 3);
        tolerance.set_value(1, ru);
        tolerance.set_value(2, rv);
        tolerance.set_value(3, c.resolution(occt_core::precision::CONFUSION));
        rsnld.set_tolerance(&tolerance);

        let mut autretentative = 0;
        self.done = false;
        while !self.done && autretentative < 3 {
            if autretentative == 1 {
                uvap.set_value(3, w0);
            } else if autretentative == 2 {
                uvap.set_value(3, w1);
            }
            autretentative += 1;
            rsnld.perform_with_bounds(&mut self.my_function, &uvap, &born_inf, &born_sup, false);
            if rsnld.is_done() {
                if self.my_function.root().abs() <= self.tol {
                    let root = rsnld.root();
                    self.u = root.value(1);
                    self.v = root.value(2);
                    self.w = root.value(3);
                    self.empty = false;
                    self.done = true;
                }
            }
        }
    }

    /// `IntImp_IntCS::IsDone` (`...gxx:154-157`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `IntImp_IntCS::IsEmpty` (`...gxx:159-164`).
    pub fn is_empty(&self) -> bool {
        assert!(self.done, "IntCurveSurface_TheExactHInter: NotDone");
        self.empty
    }

    /// `IntImp_IntCS::Point` (`...gxx:166-175`).
    pub fn point(&self) -> &GpPnt {
        assert!(self.done, "IntCurveSurface_TheExactHInter: NotDone");
        assert!(!self.empty, "IntCurveSurface_TheExactHInter: DomainError");
        self.my_function.point()
    }

    /// `IntImp_IntCS::ParameterOnSurface` (`...gxx:177-185`).
    pub fn parameter_on_surface(&self) -> (f64, f64) {
        assert!(self.done, "IntCurveSurface_TheExactHInter: NotDone");
        assert!(!self.empty, "IntCurveSurface_TheExactHInter: DomainError");
        (self.u, self.v)
    }

    /// `IntImp_IntCS::ParameterOnCurve` (`...gxx:187-193`).
    pub fn parameter_on_curve(&self) -> f64 {
        assert!(self.done, "IntCurveSurface_TheExactHInter: NotDone");
        assert!(!self.empty, "IntCurveSurface_TheExactHInter: DomainError");
        self.w
    }

    /// `IntImp_IntCS::Function` (`...gxx:195-198`).
    pub fn function(&self) -> &TheCSFunctionOfHInter<'a> {
        &self.my_function
    }
}
