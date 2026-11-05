//! Curve-surface extrema. Port of `Extrema_ExtCS`, `Extrema_GenExtCS`,
//! `Extrema_GenLocateExtCS`, `Extrema_FuncExtCS`, `Extrema_ExtElCS` and the
//! three `Extrema_GlobOptFunc*` objectives (TKGeomBase).
//!
//! | OCCT | source | port |
//! |---|---|---|
//! | `Extrema_FuncExtCS` | `Extrema_FuncExtCS.cxx:45-238`, `.hxx:35-100` | [`func_ext_cs`] |
//! | `Extrema_GlobOptFuncCS` | `Extrema_GlobOptFuncCS.cxx:28-174`, `.hxx:27-58` | [`glob_opt_func_cs`] |
//! | `Extrema_GlobOptFuncConicS` | `Extrema_GlobOptFuncConicS.cxx:22-239`, `.hxx:24-74` | [`glob_opt_func_conic_s`] |
//! | `Extrema_GlobOptFuncCQuadric` | `Extrema_GlobOptFuncCQuadric.cxx:22-274`, `.hxx:23-73` | [`glob_opt_func_cquadric`] |
//! | `Extrema_GenExtCS` | `Extrema_GenExtCS.cxx:136-877`, `.hxx:31-163` | [`gen_ext_cs`] |
//! | `Extrema_GenLocateExtCS` | `Extrema_GenLocateExtCS.cxx:28-139`, `.hxx:31-77` | [`gen_locate_ext_cs`] |
//! | `Extrema_ExtElCS` | `Extrema_ExtElCS.cxx:45-806`, `.hxx:45-140` | [`ext_el_cs`] |
//! | `Extrema_ExtCS` | `Extrema_ExtCS.cxx:37-630`, `.hxx:53-130` | [`ext_cs`] |
//!
//! The engine chain matches OCCT exactly: `Extrema_ExtCS::Perform` first tries
//! the analytic `Extrema_ExtElCS` arms for the elementary curve/surface pairs
//! (`cxx:127-297`), and falls back to `Extrema_GenExtCS` (`cxx:418-449`) with a
//! `NbT = 12` / `NbU = NbV = 10` sampling grid unless the exact arm succeeded.
//! `Extrema_GenExtCS` itself reduces the dimension of the search by surface and
//! curve type (`Extrema_GenExtCS.cxx:325-336`): 1 variable for a conic against
//! a quadric (`GlobMinCQuadric`), 2 for a conic against any other surface
//! (`GlobMinConicS`), 3 otherwise (`GlobMinGenCS`). Every arm seeds
//! `math_PSO` from pre-computed sample points and then polishes with
//! `math_FunctionSetRoot` over the [`FuncExtCS`] orthogonality system.
//!
//! `GeomGridEval_Surface::EvaluateGrid` (`GeomGridEval_Surface.cxx:340-364`,
//! `GeomGridEval_OtherSurface.cxx:112-138`) is a batched `S.Value(u, v)` over
//! the two parameter arrays; the port evaluates `d0` at the same `(u, v)`.
//!
//! **All** `Extrema_ExtElCS` arms listed above are ported; the
//! `Standard_NotImplemented` ones (`Cone`, `Torus`) are the OCCT ones.

mod prelude {

    pub(crate) use occt_core::elib::{clib, slib};
    pub(crate) use occt_core::gp::{
        GpAx3, GpCone, GpCylinder, GpDir, GpElips, GpHypr, GpLin, GpParab, GpPln, GpPnt, GpSphere,
        GpTorus, GpVec,
    };
    pub(crate) use occt_core::precision::{
        CONFUSION, PCONFUSION, SQUARE_CONFUSION, ANGULAR,
    };

    pub(crate) use crate::curve::Curve;
    pub(crate) use crate::extrema_el_c::ExtremaPOnCurv;
    pub(crate) use crate::extrema_surf::ExtremaPOnSurf;
    pub(crate) use crate::surface::Surface;

}

mod ext_cs;
mod ext_el_cs;
mod func_ext_cs;
mod gen_ext_cs;
mod gen_locate_ext_cs;
mod glob_opt_func_conic_s;
mod glob_opt_func_cquadric;
mod glob_opt_func_cs;

pub use ext_cs::*;
pub use ext_el_cs::*;
pub use func_ext_cs::*;
pub use gen_ext_cs::*;
pub use gen_locate_ext_cs::*;
pub use glob_opt_func_conic_s::*;
pub use glob_opt_func_cquadric::*;
pub use glob_opt_func_cs::*;
