//! The `Point` / `Vec` hooks of the three `BRepBlend_AppFunc*` subclasses.
//!
//! Source: `BRepBlend_AppFunc.cxx:33-44`, `BRepBlend_AppFuncRst.cxx:32-44`
//! and `BRepBlend_AppFuncRstRst.cxx:33-44` (OCCT 8.0.0); declared pure virtual
//! in `BRepBlend_AppFuncRoot.hxx:150-156`.
//!
//! `Point` packs the solver solution `Sol` back into a `Blend_Point` with the
//! `SetValue` overload of the contact kind, and `Vec` unpacks an existing
//! point into the initial solution `XInit`. Both are `const` in OCCT, so the
//! trait takes `&self`.

use occt_math::MathVector;

use crate::blend::{BlendAppFunction, BlendPoint};

/// The `Point` / `Vec` vtable of `BRepBlend_AppFuncRoot.hxx:150-156`.
pub trait AppFuncContact {
    /// `Point(Func, Param, Sol, Pnt)` (`hxx:150-152`): write the point of
    /// parameter `param` for the solution `sol`. `Err` is OCCT's
    /// `Standard_DomainError` from the `Blend_Point` accessors.
    fn point(
        &self,
        func: &dyn BlendAppFunction,
        param: f64,
        sol: &MathVector,
    ) -> Result<BlendPoint, String>;

    /// `Vec(Sol, Pnt)` (`hxx:154`): expand `pnt` into the initial solution.
    fn vec(&self, sol: &mut MathVector, pnt: &BlendPoint) -> Result<(), String>;
}

/// `BRepBlend_AppFunc::Point` / `Vec` (`cxx:33-44`): surface / surface, the
/// solution is `U1, V1, U2, V2`.
pub struct SurfSurfContact;

impl AppFuncContact for SurfSurfContact {
    /// `Pnt.SetValue(Func.Pnt1(), Func.Pnt2(), Param, Sol(1), Sol(2), Sol(3),
    /// Sol(4))` (`cxx:33-38`), i.e. the
    /// `SetValue(P1, P2, Param, U1, V1, U2, V2)` overload
    /// (`Blend_Point.cxx:111-131`).
    fn point(
        &self,
        func: &dyn BlendAppFunction,
        param: f64,
        sol: &MathVector,
    ) -> Result<BlendPoint, String> {
        let mut pnt = BlendPoint::new();
        pnt.set_on_two_surfaces(
            *func.pnt1(),
            *func.pnt2(),
            param,
            sol.value(1),
            sol.value(2),
            sol.value(3),
            sol.value(4),
        );
        Ok(pnt)
    }

    /// `Pnt.ParametersOnS1(Sol(1), Sol(2)); Pnt.ParametersOnS2(Sol(3),
    /// Sol(4))` (`cxx:40-44`).
    fn vec(&self, sol: &mut MathVector, pnt: &BlendPoint) -> Result<(), String> {
        let (u1, v1) = pnt.parameters_on_s1()?;
        let (u2, v2) = pnt.parameters_on_s2()?;
        sol.set_value(1, u1);
        sol.set_value(2, v1);
        sol.set_value(3, u2);
        sol.set_value(4, v2);
        Ok(())
    }
}

/// `BRepBlend_AppFuncRst::Point` / `Vec` (`cxx:32-44`): surface / restriction,
/// the solution is `U, V, W`.
pub struct SurfRstContact;

impl AppFuncContact for SurfRstContact {
    /// `Pnt.SetValue(Func.Pnt1(), Func.Pnt2(), Param, Sol(1), Sol(2),
    /// Sol(3))` (`cxx:32-38`), i.e. the `SetValue(Ps, Pc, Param, U, V, W)`
    /// overload (`Blend_Point.cxx:207-225`).
    fn point(
        &self,
        func: &dyn BlendAppFunction,
        param: f64,
        sol: &MathVector,
    ) -> Result<BlendPoint, String> {
        let mut pnt = BlendPoint::new();
        pnt.set_on_surface_and_curve(
            *func.pnt1(),
            *func.pnt2(),
            param,
            sol.value(1),
            sol.value(2),
            sol.value(3),
        );
        Ok(pnt)
    }

    /// `Pnt.ParametersOnS(Sol(1), Sol(2)); Sol(3) = Pnt.ParameterOnC()`
    /// (`cxx:40-44`).
    fn vec(&self, sol: &mut MathVector, pnt: &BlendPoint) -> Result<(), String> {
        let (u, v) = pnt.parameters_on_s()?;
        sol.set_value(1, u);
        sol.set_value(2, v);
        sol.set_value(3, pnt.parameter_on_c()?);
        Ok(())
    }
}

/// `BRepBlend_AppFuncRstRst::Point` / `Vec` (`cxx:33-44`): restriction /
/// restriction, the solution is `W1, W2`.
pub struct RstRstContact;

impl AppFuncContact for RstRstContact {
    /// `Pnt.SetValue(Func.Pnt1(), Func.Pnt2(), Param, Sol(1), Sol(2))`
    /// (`cxx:33-38`), i.e. the `SetValue(P1, P2, Param, PC1, PC2)` overload
    /// (`Blend_Point.cxx` `SetValue(..., PC1, PC2)`).
    fn point(
        &self,
        func: &dyn BlendAppFunction,
        param: f64,
        sol: &MathVector,
    ) -> Result<BlendPoint, String> {
        let mut pnt = BlendPoint::new();
        pnt.set_on_two_curves(*func.pnt1(), *func.pnt2(), param, sol.value(1), sol.value(2));
        Ok(pnt)
    }

    /// `Sol(1) = Pnt.ParameterOnC1(); Sol(2) = Pnt.ParameterOnC2()`
    /// (`cxx:40-44`).
    fn vec(&self, sol: &mut MathVector, pnt: &BlendPoint) -> Result<(), String> {
        sol.set_value(1, pnt.parameter_on_c1()?);
        sol.set_value(2, pnt.parameter_on_c2()?);
        Ok(())
    }
}
