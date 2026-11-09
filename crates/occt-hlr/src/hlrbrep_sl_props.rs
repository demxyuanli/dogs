//! `HLRBRep_SLProps` (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep_SLProps.hxx:29-33`)
//! and `HLRBRep_SLPropsATool` (`HLRBRep_SLPropsATool.hxx:27-70`,
//! `HLRBRep_SLPropsATool.lxx:19-99`).
//!
//! `HLRBRep_SLProps` is the alias
//! `GeomLProp_SLPropsBase<HLRBRep_SurfacePtr, LProp_SurfaceUtils::ToolAccess<
//! HLRBRep_SLPropsATool>>` (`hxx:30-33`), i.e. the surface local-properties
//! engine of the exact-HLR layer. Only the subset `HLRBRep_Data` uses is
//! ported: `SetSurface`, `SetParameters`, `IsNormalDefined` and `Normal`
//! (`HLRBRep_Data.cxx:791`, `:871-877`, `:1592-1595`, `:1802-1819`,
//! `:1909-1912`).
//!
//! ## Ported
//!
//! * [`HlrBRepSlProps::set_surface`] - `SetSurface` (`GeomLProp_SLProps.hxx:120-124`):
//!   stores the surface and resets `myCN` to 4.
//! * [`HlrBRepSlProps::set_parameters`] - `SetParameters`
//!   (`GeomLProp_SLProps.hxx:128-145` + `GeomLProp_SurfaceUtils.hxx:486-511`):
//!   stores `(U, V)`, evaluates the derivatives through
//!   `EvalSurfDerivatives` (`GeomLProp_SurfaceUtils.hxx:198-217`) and resets
//!   the four `LProp_Status` fields.
//! * [`HlrBRepSlProps::is_normal_defined`] - `IsNormalDefined`
//!   (`GeomLProp_SLProps.hxx:301-305` + `GeomLProp_SurfaceUtils.hxx:609-629`):
//!   the `LProp_Undefined` short-circuit, the `>= LProp_Defined` cache hit and
//!   the `ComputeSurfNormal` (`GeomLProp_SurfaceUtils.hxx:350-358`) fallback,
//!   which is `CSLib::Normal(D1U, D1V, SinTol, ...)` (`CSLib.cxx:45-80`, ported
//!   as [`occt_core::cslib::normal_d1`]).
//! * [`HlrBRepSlProps::normal`] - `Normal` (`GeomLProp_SurfaceUtils.hxx:632-636`):
//!   throws `LProp_NotDefined` when the normal is not defined -> `None`.
//!
//! ## UNPORTED
//!
//! * `IsTangentUDefined` / `IsTangentVDefined` / `TangentU` / `TangentV`
//!   (`GeomLProp_SLProps.hxx:231-299`) and `IsCurvatureDefined` / `IsUmbilic` /
//!   `MaxCurvature` / `MinCurvature` / `MeanCurvature` / `GaussianCurvature` /
//!   `CurvatureDirections` (`hxx:307-369`): no `HLRBRep_Data` path reaches them
//!   (`mySLProps` is only ever asked for the first derivatives and the normal).
//! * `HLRBRep_SLPropsATool::DN` (`lxx:72-80`) and `Continuity` / `Bounds`
//!   (`lxx:84-99`): `DN` has no caller in the ported set, and the constant `2`
//!   / `RealFirst..RealLast` bounds are only read by `IsCurvatureDefined`.
//! * The second derivatives evaluated when `der_order == 2` (which is what
//!   `HLRBRep_Data` constructs, `cxx:532`) are dropped: no HLR caller reads
//!   `D2U`/`D2V`/`DUV`. When the adaptor cannot answer `D2` the evaluation
//!   falls back to `D1`; OCCT would let the `Standard_ConstructionError` out of
//!   `HLRBRep_Surface::D2` (`HLRBRep_Surface.lxx:151-163`).
//!
//! `myLinTol` is `Epsilon(1.)` for the `HLRBRep_Data` instance (`cxx:532`), the
//! value `ComputeSurfNormal` passes to `CSLib::Normal` as `SinTol`.

use occt_core::cslib::{normal_d1, CSLibDerivativeStatus};
use occt_core::gp::{GpDir, GpPnt, GpVec};
use occt_core::precision::REAL_EPSILON;
use occt_geom::geom_lprop::LPropStatus;

use crate::hlrbrep_surface::HlrBRepSurface;

/// `HLRBRep_SLProps` = `GeomLProp_SLPropsBase<HLRBRep_SurfacePtr, ToolAccess>`.
#[derive(Clone)]
pub struct HlrBRepSlProps {
    /// `SurfaceType mySurf` (`GeomLProp_SLProps.hxx:376`): `HLRBRep_SurfacePtr`,
    /// a non-owning surface pointer in OCCT. The port owns a clone.
    surface: Option<HlrBRepSurface>,
    /// `double myU` (`hxx:377`), initialised to `RealLast()`.
    u: f64,
    /// `double myV` (`hxx:378`), initialised to `RealLast()`.
    v: f64,
    /// `int myDerOrder` (`hxx:379`).
    der_order: i32,
    /// `int myCN` (`hxx:380`): `HLRBRep_SLPropsATool::Continuity` answers `2`,
    /// `SetSurface` overwrites it with 4.
    cn: i32,
    /// `double myLinTol` (`hxx:381`).
    lin_tol: f64,
    /// `gp_Pnt myPnt` (`hxx:382`).
    pnt: GpPnt,
    /// `gp_Vec myD1u` (`hxx:383`).
    d1u: GpVec,
    /// `gp_Vec myD1v` (`hxx:384`).
    d1v: GpVec,
    /// `gp_Dir myNormal` (`hxx:389`).
    normal: Option<GpDir>,
    /// `LProp_Status myNormalStatus` (`hxx:397`).
    normal_status: LPropStatus,
    /// `LProp_Status myUTangentStatus` / `myVTangentStatus` (`hxx:395-396`),
    /// reset by `SetParameters` and otherwise unread here.
    u_tangent_status: LPropStatus,
    v_tangent_status: LPropStatus,
    /// `LProp_Status myCurvatureStatus` (`hxx:398`), reset by `SetParameters`.
    curvature_status: LPropStatus,
}

impl Default for HlrBRepSlProps {
    /// `GeomLProp_SLPropsBase(const int N, const double Resolution)`
    /// (`GeomLProp_SLProps.hxx:100-115`): `myCN(0)`, `myU = myV = RealLast()`,
    /// every status `LProp_Undecided`, no surface.
    fn default() -> Self {
        Self {
            surface: None,
            u: f64::MAX,
            v: f64::MAX,
            der_order: 0,
            cn: 0,
            lin_tol: 0.0,
            pnt: GpPnt::default(),
            d1u: GpVec::default(),
            d1v: GpVec::default(),
            normal: None,
            normal_status: LPropStatus::Undecided,
            u_tangent_status: LPropStatus::Undecided,
            v_tangent_status: LPropStatus::Undecided,
            curvature_status: LPropStatus::Undecided,
        }
    }
}

impl HlrBRepSlProps {
    /// `HLRBRep_SLProps(const int N, const double Resolution)`
    /// (`GeomLProp_SLProps.hxx:100-115`): the constructor `HLRBRep_Data` uses
    /// (`cxx:532`, `N = 2`, `Resolution = Epsilon(1.)`).
    pub fn new(der_order: i32, resolution: f64) -> Self {
        Self {
            der_order,
            lin_tol: resolution,
            ..Self::default()
        }
    }

    /// The default `HLRBRep_Data` instance: `mySLProps(2, Epsilon(1.))`
    /// (`HLRBRep_Data.cxx:532`).
    pub fn for_hlr_data() -> Self {
        Self::new(2, REAL_EPSILON)
    }

    /// `SetSurface(S)` (`GeomLProp_SLProps.hxx:120-124`): `mySurf = S;
    /// myCN = 4;`.
    pub fn set_surface(&mut self, s: &HlrBRepSurface) {
        self.surface = Some(s.clone());
        self.cn = 4;
    }

    /// `SetParameters(U, V)` (`GeomLProp_SLProps.hxx:128-145`)
    /// + `LProp_SurfaceUtils::SetParameters` (`GeomLProp_SurfaceUtils.hxx:486-511`).
    pub fn set_parameters(&mut self, u: f64, v: f64) {
        self.u = u;
        self.v = v;
        // `EvalSurfDerivatives<Access>(theSurf, theU, theV, theDerOrder, ...)`
        // (`GeomLProp_SurfaceUtils.hxx:198-217`).
        if let Some(s) = &self.surface {
            match self.der_order {
                0 => {
                    if let Some(p) = s.d0(u, v) {
                        self.pnt = p;
                    }
                }
                1 => {
                    if let Some((p, d1u, d1v)) = s.d1(u, v) {
                        self.pnt = p;
                        self.d1u = d1u;
                        self.d1v = d1v;
                    }
                }
                _ => match s.d2(u, v) {
                    Some((p, d1u, d1v, _d2u, _d2v, _duv)) => {
                        self.pnt = p;
                        self.d1u = d1u;
                        self.d1v = d1v;
                    }
                    // `HLRBRep_Surface::D2` throws in OCCT when the adaptor
                    // cannot answer; the port keeps the first derivatives.
                    None => {
                        if let Some((p, d1u, d1v)) = s.d1(u, v) {
                            self.pnt = p;
                            self.d1u = d1u;
                            self.d1v = d1v;
                        }
                    }
                },
            }
        }
        // `theUTanSt = theVTanSt = theNormSt = theCurvSt = LProp_Undecided;`
        // (`GeomLProp_SurfaceUtils.hxx:507-510`).
        self.u_tangent_status = LPropStatus::Undecided;
        self.v_tangent_status = LPropStatus::Undecided;
        self.normal_status = LPropStatus::Undecided;
        self.curvature_status = LPropStatus::Undecided;
    }

    /// `Value()` (`GeomLProp_SLProps.hxx:148`): the point evaluated by
    /// `SetParameters`.
    pub fn value(&self) -> &GpPnt {
        &self.pnt
    }

    /// `D1U()` / `D1V()` (`GeomLProp_SLProps.hxx:152-176`). OCCT recomputes
    /// them through `EnsureSurfDeriv` when `myDerOrder < 1`; here they are the
    /// values `SetParameters` stored.
    pub fn d1(&self) -> (&GpVec, &GpVec) {
        (&self.d1u, &self.d1v)
    }

    /// `IsNormalDefined()` (`GeomLProp_SLProps.hxx:301-305` +
    /// `GeomLProp_SurfaceUtils.hxx:609-629`).
    pub fn is_normal_defined(&mut self) -> bool {
        if self.normal_status == LPropStatus::Undefined {
            return false;
        }
        if self.normal_status >= LPropStatus::Defined {
            return true;
        }
        // `ComputeSurfNormal(theD1u, theD1v, theLinTol, theNormal)`
        // (`GeomLProp_SurfaceUtils.hxx:350-358`).
        let (status, n) = normal_d1(&self.d1u, &self.d1v, self.lin_tol);
        if status == CSLibDerivativeStatus::Done {
            self.normal = n;
            self.normal_status = LPropStatus::Computed;
            return true;
        }
        self.normal_status = LPropStatus::Undefined;
        false
    }

    /// `Normal()` (`GeomLProp_SurfaceUtils.hxx:632-636`): throws
    /// `LProp_NotDefined` unless `IsNormalDefined()`.
    pub fn normal(&self) -> Option<&GpDir> {
        if self.normal_status >= LPropStatus::Defined {
            self.normal.as_ref()
        } else {
            None
        }
    }
}
