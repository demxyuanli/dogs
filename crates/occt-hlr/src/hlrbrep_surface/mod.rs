//! `HLRBRep_Surface` (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep_Surface.hxx:42-198`,
//! `.cxx:27-336`, `.lxx:21-274`).
//!
//! The surface side of the exact-HLR geometry layer: a `BRepAdaptor_Surface`
//! plus the cached classified type (`myType`) and the projector pointer. The
//! adaptor itself is [`SurfaceAdaptor`] (`BRepAdaptor_Surface`); every
//! `HLRBRep_BSurfaceTool::*` call in the `.lxx` is `GeomAdaptor_Surface::*` on
//! that adaptor, so the delegation is direct.
//!
//! UNPORTED / deviations (each with its OCCT line):
//! * `UIntervalContinuity()` / `VIntervalContinuity()` (`hxx:87-89`) are
//!   declared in the header and never defined in `HLRBRep_Surface.cxx` /
//!   `HLRBRep_Surface.lxx` - dead declarations, nothing to port.
//! * `UContinuity()` / `VContinuity()`: the BSpline arm of
//!   `GeomAdaptor_Surface::UContinuity` (`GeomAdaptor_Surface.cxx:543-554`)
//!   changes the port's answer to `None`. It needs
//!   `GeomAdaptor_Surface::LocalContinuity` (`cxx:110-165`: distinct knots plus
//!   multiplicities, the adaptor's restricted `myUFirst`/`myULast`, and
//!   `BSplCLib::LocateParameter`). `occt_core::bspl::intervals::local_continuity`
//!   is the *curve* function (`GeomAdaptor_Curve.cxx:143-225`) with its own
//!   index walk, so it is not the same branch. The remaining arms are ported:
//!   analytic and Bezier surfaces answer `GeomAbs_CN` (`cxx:577-586`), the
//!   offset arm folds the basis continuity (`cxx:555-570`), the revolution /
//!   extrusion arm answers the basis curve's continuity (`cxx:572-574`,
//!   `:624-626`), and `GeomAbs_OtherSurface` answers `None` (OCCT raises
//!   `Standard_NoSuchObject`, `cxx:575-576`).
//! * Offset arm, basis query: OCCT asks the *basis adaptor's* direction-specific
//!   `UContinuity()` / `VContinuity()` (`cxx:556`, `:609`). The port has only
//!   `Geom_Surface::Continuity()` (`cxx` `Geom_BSplineSurface::Continuity` is
//!   `Max(U, V)`), so a basis whose two directions differ yields the smoother
//!   of the two. `DirectionAdaptor`-level continuity is not ported.
//! * `SideRowsOfPoles` needs `GProp_PEquation`
//!   (`GProp_PEquation.cxx:28-177`), ported as
//!   [`occt_geom::gprop_pequation::GPropPEquation`].
//! * `D3` (`hxx:132-144`) has no `TryD3` arm in this port's `Surface` trait
//!   either: [`HLRBRepSurface::d3`] returns `None` where the trait's
//!   `Surface::d3` cannot answer (documented in `occt-geom/src/surface.rs`).

use std::sync::Arc;

use occt_core::gp::{GpAx1, GpAx3, GpCone, GpCylinder, GpDir, GpPln, GpPnt, GpSphere, GpTorus, GpVec};
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_geom::surface_d3::SurfD3;

use occt_topo::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::hlr_algo::HlrAlgoProjector;
use occt_topo::shape::Face;
use occt_topo::topopebrep::tool::adaptor::SurfaceAdaptor;

mod side;
mod adaptor;

pub use adaptor::HlrBRepSurfaceView;

/// `HLRBRep_Surface`.
#[derive(Clone)]
pub struct HlrBRepSurface {
    /// `BRepAdaptor_Surface mySurf` (`hxx:190`).
    my_surf: SurfaceAdaptor,
    /// `GeomAbs_SurfaceType myType` (`hxx:191`): the classified type.
    type_: GeomAbsSurfaceType,
    /// `const HLRAlgo_Projector* myProj` (`hxx:192`).
    proj: Option<Arc<HlrAlgoProjector>>,
}

impl Default for HlrBRepSurface {
    /// `HLRBRep_Surface()` (`cxx:28-31`): `myType = GeomAbs_OtherSurface`,
    /// `myProj = nullptr`, `mySurf` uninitialized.
    fn default() -> Self {
        Self {
            my_surf: SurfaceAdaptor::new(),
            type_: GeomAbsSurfaceType::OtherSurface,
            proj: None,
        }
    }
}

impl HlrBRepSurface {
    /// `HLRBRep_Surface()` (`cxx:28-31`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `Projector(const HLRAlgo_Projector*)` (`hxx:49`).
    pub fn set_projector(&mut self, proj: Option<Arc<HlrAlgoProjector>>) {
        self.proj = proj;
    }

    /// The stored `myProj` pointer.
    pub fn projector(&self) -> Option<&Arc<HlrAlgoProjector>> {
        self.proj.as_ref()
    }

    /// `Surface()` / `Surface() const` (`lxx:23-33`): the adaptor itself.
    pub fn surface_adaptor(&self) -> &SurfaceAdaptor {
        &self.my_surf
    }

    /// The mutable `Surface()` overload.
    pub fn surface_adaptor_mut(&mut self) -> &mut SurfaceAdaptor {
        &mut self.my_surf
    }

    /// `Surface(const TopoDS_Face&)` (`cxx:34-67`): `mySurf.Initialize(F, true)`
    /// (`BRepAdaptor_Surface.cxx:70-76`, the restricted-range overload), then
    /// the type classification.
    pub fn set_surface(&mut self, f: &Face) {
        self.my_surf = SurfaceAdaptor::of_face_restricted(f);
        let kind = self.my_surf.get_type();
        self.type_ = match kind {
            GeomAbsSurfaceType::Plane
            | GeomAbsSurfaceType::Cylinder
            | GeomAbsSurfaceType::Cone
            | GeomAbsSurfaceType::Sphere
            | GeomAbsSurfaceType::Torus => kind,
            // `cxx:53-59`: a Bezier patch of degree 1 in both directions is a
            // plane, anything else keeps its Bezier type.
            GeomAbsSurfaceType::BezierSurface => {
                if self.surface().map(|s| (s.u_degree(), s.v_degree())) == Some((1, 1)) {
                    GeomAbsSurfaceType::Plane
                } else {
                    kind
                }
            }
            // `cxx:63-66`.
            _ => GeomAbsSurfaceType::OtherSurface,
        };
    }

    /// The underlying `Handle(Geom_Surface)` (`SurfaceAdaptor::surface`).
    pub fn surface(&self) -> Option<&Arc<dyn Surface>> {
        self.my_surf.surface()
    }

    /// `GetType()` (`lxx:191-194`): the cached `myType`.
    pub fn get_type(&self) -> GeomAbsSurfaceType {
        self.type_
    }

    /// The adaptor's own type (`HLRBRep_BSurfaceTool::GetType(mySurf)`), which
    /// `Plane()` and `IsAbove` switch on.
    pub fn adaptor_type(&self) -> GeomAbsSurfaceType {
        self.my_surf.get_type()
    }

    /// `FirstUParameter()` (`lxx:37-40`).
    pub fn first_u_parameter(&self) -> f64 {
        self.my_surf.first_u_parameter()
    }

    /// `LastUParameter()` (`lxx:44-47`).
    pub fn last_u_parameter(&self) -> f64 {
        self.my_surf.last_u_parameter()
    }

    /// `FirstVParameter()` (`lxx:51-54`).
    pub fn first_v_parameter(&self) -> f64 {
        self.my_surf.first_v_parameter()
    }

    /// `LastVParameter()` (`lxx:58-61`).
    pub fn last_v_parameter(&self) -> f64 {
        self.my_surf.last_v_parameter()
    }

    /// `UContinuity()` (`lxx:65-68`). See the module notes for the BSpline arm.
    pub fn u_continuity(&self) -> Option<u8> {
        self.direction_continuity(true)
    }

    /// `VContinuity()` (`lxx:72-75`).
    pub fn v_continuity(&self) -> Option<u8> {
        self.direction_continuity(false)
    }

    fn direction_continuity(&self, along_u: bool) -> Option<u8> {
        /// `GeomAbs_CN` (`GeomAbs_Shape.hxx:31`).
        const GEOM_ABS_CN: u8 = 6;
        let s = self.surface()?;
        match self.my_surf.get_type() {
            GeomAbsSurfaceType::BSplineSurface => None,
            GeomAbsSurfaceType::OffsetSurface => {
                // `GeomAdaptor_Surface.cxx:555-570`: the basis continuity folds
                // to `GeomAbs_CN` for `CN`/`C3` and to `GeomAbs_C1` for
                // `G2`/`C2`; `G1`, `C1` and `C0` all answer `GeomAbs_C0`.
                const GEOM_ABS_G2: u8 = 3;
                const GEOM_ABS_C2: u8 = 4;
                const GEOM_ABS_C3: u8 = 5;
                const GEOM_ABS_C1: u8 = 2;
                const GEOM_ABS_C0: u8 = 0;
                let basis = s.offset_basis_surface()?;
                let c = basis.continuity();
                Some(if c == GEOM_ABS_CN || c == GEOM_ABS_C3 {
                    GEOM_ABS_CN
                } else if c == GEOM_ABS_G2 || c == GEOM_ABS_C2 {
                    GEOM_ABS_C1
                } else {
                    GEOM_ABS_C0
                })
            }
            GeomAbsSurfaceType::SurfaceOfExtrusion if along_u => {
                s.extrusion_basis_curve().map(|c| c.continuity())
            }
            GeomAbsSurfaceType::SurfaceOfRevolution if !along_u => {
                s.revolution_basis_curve().map(|c| c.continuity())
            }
            GeomAbsSurfaceType::OtherSurface => None,
            _ => Some(GEOM_ABS_CN),
        }
    }

    /// `NbUIntervals(S)` (`lxx:79-82`).
    pub fn nb_u_intervals(&self, s: u8) -> i32 {
        match self.surface() {
            Some(surf) => surf.nb_u_intervals(s),
            None => 0,
        }
    }

    /// `NbVIntervals(S)` (`lxx:85-88`).
    pub fn nb_v_intervals(&self, s: u8) -> i32 {
        match self.surface() {
            Some(surf) => surf.nb_v_intervals(s),
            None => 0,
        }
    }

    /// `IsUClosed()` (`lxx:92-95`).
    pub fn is_u_closed(&self) -> bool {
        self.surface().map(|s| s.is_u_closed()).unwrap_or(false)
    }

    /// `IsVClosed()` (`lxx:97-100`).
    pub fn is_v_closed(&self) -> bool {
        self.surface().map(|s| s.is_v_closed()).unwrap_or(false)
    }

    /// `IsUPeriodic()` (`lxx:104-107`).
    pub fn is_u_periodic(&self) -> bool {
        self.surface().map(|s| s.is_u_periodic()).unwrap_or(false)
    }

    /// `UPeriod()` (`lxx:111-114`).
    pub fn u_period(&self) -> f64 {
        self.surface().map(|s| s.u_period()).unwrap_or(0.0)
    }

    /// `IsVPeriodic()` (`lxx:118-121`).
    pub fn is_v_periodic(&self) -> bool {
        self.surface().map(|s| s.is_v_periodic()).unwrap_or(false)
    }

    /// `VPeriod()` (`lxx:125-128`).
    pub fn v_period(&self) -> f64 {
        self.surface().map(|s| s.v_period()).unwrap_or(0.0)
    }

    /// `D0(U, V, P)` (`lxx:134-137`).
    pub fn d0(&self, u: f64, v: f64) -> Option<GpPnt> {
        self.my_surf.d0(u, v)
    }

    /// `D1(U, V, P, D1U, D1V)` (`lxx:141-148`).
    pub fn d1(&self, u: f64, v: f64) -> Option<(GpPnt, GpVec, GpVec)> {
        self.my_surf.d1(u, v)
    }

    /// `D2(U, V, P, D1U, D1V, D2U, D2V, D2UV)` (`lxx:151-163`).
    #[allow(clippy::type_complexity)]
    pub fn d2(&self, u: f64, v: f64) -> Option<(GpPnt, GpVec, GpVec, GpVec, GpVec, GpVec)> {
        self.my_surf.d2(u, v)
    }

    /// `D3(U, V, ...)` (`lxx:165-180`).
    pub fn d3(&self, u: f64, v: f64) -> Option<SurfD3> {
        self.surface()?.d3(u, v)
    }

    /// `DN(U, V, Nu, Nv)` (`lxx:184-187`): `GeomAdaptor_Surface::DN`
    /// (`GeomAdaptor_Surface.cxx:1745-1812`), which for the analytic types is
    /// `ElSLib::*DN` (`ElSLib.cxx:169-557`) and otherwise
    /// `Geom_Surface::EvalDN`. Both are the same calls this port's `Surface`
    /// implementations make in `try_eval_dn`.
    pub fn dn(&self, u: f64, v: f64, nu: i32, nv: i32) -> Result<GpVec, String> {
        match self.surface() {
            Some(s) => s.try_eval_dn(u, v, nu, nv),
            None => Err("HLRBRep_Surface::DN: no surface".into()),
        }
    }

    /// `Value(U, V)` (`cxx:310-315`): `D0` and return the point.
    pub fn value(&self, u: f64, v: f64) -> Option<GpPnt> {
        self.d0(u, v)
    }

    /// `Plane()` (`cxx:319-335`): a Bezier patch answers the tangent plane at
    /// (0.5, 0.5), every other type answers `GeomAdaptor_Surface::Plane`
    /// (`GeomAdaptor_Surface.cxx:1962-1969`), which throws unless the adaptor
    /// type is `GeomAbs_Plane`. `None` is that throw.
    pub fn plane(&self) -> Option<GpPln> {
        match self.my_surf.get_type() {
            GeomAbsSurfaceType::BezierSurface => {
                let (p, d1u, d1v) = self.d1(0.5, 0.5)?;
                let n = d1u.crossed(&d1v);
                let dir = GpDir::from_xyz(&n.coord).ok()?;
                Some(GpPln::new(GpAx3::from_ax1(&GpAx1::new(p, dir))))
            }
            GeomAbsSurfaceType::Plane => self.surface()?.gp_pln(),
            _ => None,
        }
    }

    /// `Cylinder()` (`lxx:198-201`): throws in OCCT unless the adaptor type
    /// matches, so `None` here.
    pub fn cylinder(&self) -> Option<GpCylinder> {
        self.surface()?.gp_cylinder()
    }

    /// `Cone()` (`lxx:205-208`).
    pub fn cone(&self) -> Option<GpCone> {
        self.surface()?.gp_cone()
    }

    /// `Sphere()` (`lxx:212-215`).
    pub fn sphere(&self) -> Option<GpSphere> {
        self.surface()?.gp_sphere()
    }

    /// `Torus()` (`lxx:219-222`).
    pub fn torus(&self) -> Option<GpTorus> {
        self.surface()?.gp_torus()
    }

    /// `Axis()` (`lxx:268-271`): `GeomAdaptor_Surface::AxeOfRevolution`, which
    /// throws unless the surface is a surface of revolution.
    pub fn axis(&self) -> Option<GpAx1> {
        self.surface()?.revolution_axis()
    }

    /// `UDegree()` (`lxx:227-230`), `GeomAdaptor_Surface::UDegree`
    /// (`GeomAdaptor_Surface.cxx:2017-2032`): BSpline / Bezier patches answer
    /// their degree, a surface of extrusion answers the basis curve degree,
    /// every other type throws -> `None`.
    pub fn u_degree(&self) -> Option<i32> {
        let s = self.surface()?;
        match self.my_surf.get_type() {
            GeomAbsSurfaceType::BSplineSurface | GeomAbsSurfaceType::BezierSurface => {
                Some(s.u_degree())
            }
            GeomAbsSurfaceType::SurfaceOfExtrusion => {
                s.extrusion_basis_curve().and_then(|c| curve_degree(c.as_ref()))
            }
            _ => None,
        }
    }

    /// `VDegree()` (`lxx:239-242`), `GeomAdaptor_Surface.cxx:2055-2075`.
    pub fn v_degree(&self) -> Option<i32> {
        let s = self.surface()?;
        match self.my_surf.get_type() {
            GeomAbsSurfaceType::BSplineSurface | GeomAbsSurfaceType::BezierSurface => {
                Some(s.v_degree())
            }
            GeomAbsSurfaceType::SurfaceOfRevolution => {
                s.revolution_basis_curve().and_then(|c| curve_degree(c.as_ref()))
            }
            _ => None,
        }
    }

    /// `NbUPoles()` (`lxx:232-235`), `GeomAdaptor_Surface.cxx:2036-2051`.
    pub fn nb_u_poles(&self) -> Option<i32> {
        let s = self.surface()?;
        match self.my_surf.get_type() {
            GeomAbsSurfaceType::BSplineSurface | GeomAbsSurfaceType::BezierSurface => {
                Some(s.nb_u_poles())
            }
            GeomAbsSurfaceType::SurfaceOfExtrusion => {
                s.extrusion_basis_curve().and_then(|c| curve_nb_poles(c.as_ref()))
            }
            _ => None,
        }
    }

    /// `NbVPoles()` (`lxx:247-250`), `GeomAdaptor_Surface.cxx:2082-2091`.
    pub fn nb_v_poles(&self) -> Option<i32> {
        let s = self.surface()?;
        match self.my_surf.get_type() {
            GeomAbsSurfaceType::BSplineSurface | GeomAbsSurfaceType::BezierSurface => {
                Some(s.nb_v_poles())
            }
            GeomAbsSurfaceType::SurfaceOfRevolution => {
                s.revolution_basis_curve().and_then(|c| curve_nb_poles(c.as_ref()))
            }
            _ => None,
        }
    }

    /// `NbUKnots()` (`lxx:254-257`), `GeomAdaptor_Surface.cxx:2093-2104`.
    pub fn nb_u_knots(&self) -> Option<i32> {
        let s = self.surface()?;
        match self.my_surf.get_type() {
            GeomAbsSurfaceType::BSplineSurface => Self::distinct_count(s.bspline_surface_uknots()?),
            GeomAbsSurfaceType::SurfaceOfExtrusion => {
                let c = s.extrusion_basis_curve()?;
                curve_nb_knots(c.as_ref())
            }
            _ => None,
        }
    }

    /// `NbVKnots()` (`lxx:260-263`), `GeomAdaptor_Surface.cxx:2108-2115`.
    pub fn nb_v_knots(&self) -> Option<i32> {
        let s = self.surface()?;
        match self.my_surf.get_type() {
            GeomAbsSurfaceType::BSplineSurface => Self::distinct_count(s.bspline_surface_vknots()?),
            _ => None,
        }
    }

    /// `Geom_BSplineSurface::NbUKnots()` / `NbVKnots()`: the number of
    /// *distinct* knots of the stored knot vector.
    fn distinct_count(knots: &[f64]) -> Option<i32> {
        let mut n = 0i32;
        let mut previous: Option<f64> = None;
        for k in knots {
            if previous != Some(*k) {
                n += 1;
                previous = Some(*k);
            }
        }
        Some(n)
    }
}

/// `GeomAdaptor_Curve::Degree()` (`GeomAdaptor_Curve.cxx:1203-1218`): a Bezier
/// curve answers its degree, a BSpline curve its `Degree()`, every other type
/// throws `Standard_NoSuchObject` -> `None`.
fn curve_degree(c: &dyn Curve) -> Option<i32> {
    match c.bezier_poles() {
        Some(poles) => Some(poles.len() as i32 - 1),
        None => c.nurbs_degree().map(|d| d as i32),
    }
}

/// `GeomAdaptor_Curve::NbPoles()` (`GeomAdaptor_Curve.cxx:1236-1251`).
fn curve_nb_poles(c: &dyn Curve) -> Option<i32> {
    match c.bezier_poles() {
        Some(poles) => Some(poles.len() as i32),
        None => c.bspline_poles().map(|poles| poles.len() as i32),
    }
}

/// `GeomAdaptor_Curve::NbKnots()` (`GeomAdaptor_Curve.cxx:1254-1261`): only a
/// BSpline curve answers, every other type throws.
fn curve_nb_knots(c: &dyn Curve) -> Option<i32> {
    HlrBRepSurface::distinct_count(c.bspline_knots()?)
}
