//! Bodies shared by the `BlendFunc_GenChamfer` family.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_GenChamfer.cxx:205-284`
//! (the two `Section` overloads) and the `math_Gauss` pivot the throat
//! solvers use (`math_Gauss` defaults to `1.e-20`).
//!
//! OCCT writes these bodies once in `BlendFunc_GenChamfer` and reaches the
//! derived class through the virtual `Value` / `Derivatives` / `Values` /
//! `IsSolution` / `PointOnS*` / `Tangent*` methods. Rust has no inheritance, so
//! the same body is a generic function over the interface traits; the caller
//! passes the base-class `tol` the C++ body reads directly.

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_math::{MathMatrix, MathVector};

use crate::blend::{BlendFunction, BlendPoint};

/// The pivot `math_Gauss` uses when constructed from a matrix alone
/// (`math_Gauss::math_Gauss` defaults `Eps` to `1.e-20`).
pub(super) const GAUSS_MIN_PIVOT: f64 = 1.0e-20;

/// `p.XYZ()` as a `gp_Vec`.
pub(super) fn pnt_vec(p: &GpPnt) -> GpVec {
    GpVec::from_xyz(p.xyz())
}

/// `gp_Pnt((a.XYZ() + b.XYZ()) / 2)`.
pub(super) fn mid_pnt(a: &GpPnt, b: &GpPnt) -> GpPnt {
    GpPnt::from_xyz(&a.xyz().added(b.xyz()).multiplied(0.5))
}

/// `BlendFunc_GenChamfer::Section(P, Poles, DPoles, Poles2d, DPoles2d, Weights,
/// DWeights)` (`BlendFunc_GenChamfer.cxx:205-254`).
#[allow(clippy::too_many_arguments)]
pub(super) fn gen_section_with_derivatives<T>(
    f: &mut T,
    tol: f64,
    p: &BlendPoint,
    poles: &mut [GpPnt],
    dpoles: &mut [GpVec],
    poles2d: &mut [GpPnt2d],
    dpoles2d: &mut [GpVec2d],
    weights: &mut [f64],
    dweights: &mut [f64],
) -> bool
where
    T: BlendFunction + ?Sized,
{
    let mut sol = MathVector::new(1, 4);
    let mut valsol = MathVector::new(1, 4);
    let mut gradsol = MathMatrix::new(1, 4, 1, 4);

    let prm = p.parameter();
    let low = 0usize;
    let upp = poles.len() - 1;
    let low2d = 0usize;
    let upp2d = poles2d.len() - 1;

    let (u1, v1) = p
        .parameters_on_s1()
        .expect("BlendFunc_GenChamfer::Section: point is not on the first surface");
    let (u2, v2) = p
        .parameters_on_s2()
        .expect("BlendFunc_GenChamfer::Section: point is not on the second surface");
    sol.set_value(1, u1);
    sol.set_value(2, v1);
    sol.set_value(3, u2);
    sol.set_value(4, v2);

    f.set(prm);

    f.values(&sol, &mut valsol, &mut gradsol);
    f.is_solution(&sol, tol);

    let istgt = f.is_tangency_point();

    poles2d[low2d] = GpPnt2d::new(sol.value(1), sol.value(2));
    poles2d[upp2d] = GpPnt2d::new(sol.value(3), sol.value(4));
    if !istgt {
        let t = *f.tangent2d_on_s1();
        dpoles2d[low2d] = GpVec2d::new(t.x(), t.y());
        let t = *f.tangent2d_on_s2();
        dpoles2d[upp2d] = GpVec2d::new(t.x(), t.y());
    }
    poles[low] = *f.point_on_s1();
    poles[upp] = *f.point_on_s2();
    weights[low] = 1.0;
    weights[upp] = 1.0;
    if !istgt {
        dpoles[low] = *f.tangent_on_s1();
        dpoles[upp] = *f.tangent_on_s2();
        dweights[low] = 0.0;
        dweights[upp] = 0.0;
    }

    !istgt
}

/// `BlendFunc_GenChamfer::Section(P, Poles, Poles2d, Weights)`
/// (`BlendFunc_GenChamfer.cxx:256-284`).
pub(super) fn gen_section_poles<T>(
    f: &mut T,
    p: &BlendPoint,
    poles: &mut [GpPnt],
    poles2d: &mut [GpPnt2d],
    weights: &mut [f64],
) where
    T: BlendFunction + ?Sized,
{
    let prm = p.parameter();
    let low = 0usize;
    let upp = poles.len() - 1;
    let low2d = 0usize;
    let upp2d = poles2d.len() - 1;

    let mut x = MathVector::new(1, 4);
    let mut fv = MathVector::new(1, 4);

    let (u1, v1) = p
        .parameters_on_s1()
        .expect("BlendFunc_GenChamfer::Section: point is not on the first surface");
    let (u2, v2) = p
        .parameters_on_s2()
        .expect("BlendFunc_GenChamfer::Section: point is not on the second surface");
    x.set_value(1, u1);
    x.set_value(2, v1);
    x.set_value(3, u2);
    x.set_value(4, v2);
    poles2d[low2d] = GpPnt2d::new(u1, v1);
    poles2d[upp2d] = GpPnt2d::new(u2, v2);

    f.set(prm);
    f.value(&x, &mut fv);
    poles[low] = *f.point_on_s1();
    poles[upp] = *f.point_on_s2();
    weights[low] = 1.0;
    weights[upp] = 1.0;
}
