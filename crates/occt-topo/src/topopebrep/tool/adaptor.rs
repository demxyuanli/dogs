//! The slice of `BRepAdaptor_Surface` / `BRepAdaptor_Curve` that
//! `TopOpeBRepTool_ShapeTool` reads: `GetType()`, `FirstUParameter()` /
//! `FirstVParameter()` / `FirstParameter()`, and `D1`.
//!
//! Sources:
//! - `TKBRep/BRepAdaptor/BRepAdaptor_Surface.cxx:56-82` (`Initialize`)
//! - `TKBRep/BRepAdaptor/BRepAdaptor_Surface.hxx:28-40` (`Load(S, Trsf)`)
//! - `TKG3d/GeomAdaptor/GeomAdaptor_Surface.cxx:402-431` (`load`: type chain)
//! - `TKG3d/GeomAdaptor/GeomAdaptor_Surface.hxx:150-190`
//! - `TKBRep/BRepAdaptor/BRepAdaptor_Curve.hxx:36-98`
//!
//! The adaptor applies the shape's `TopLoc_Location` to every evaluated point
//! and derivative (`GeomAdaptor_Surface::D1` / `GeomAdaptor_Curve::D1` return
//! `myTrsf`-transformed values), so these structs hold the world-coordinate
//! geometry, matching [`crate::brep_tool::BRepTool::face_surface_world`] and
//! [`crate::brep_tool::BRepTool::edge_curve_adaptor_world`].

use std::sync::Arc;

use occt_core::gp::{
    GpCirc, GpCone, GpCylinder, GpElips, GpHypr, GpLin, GpParab, GpPln, GpPnt, GpTrsf, GpVec,
};
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;

use crate::brep_tool::BRepTool;
use crate::geom_abs_curve_type::{curve_kind, GeomAbsCurveType};
use crate::geom_bnd_lib_surface3d::{surface_kind, GeomAbsSurfaceType};
use crate::shape::{Edge, Face};

/// `BRepAdaptor_Surface(F, Restriction)` (`BRepAdaptor_Surface.cxx:32-35`,
/// `:56-82`).
#[derive(Clone)]
pub struct SurfaceAdaptor {
    /// `GeomAdaptor_Surface::Surface()` after `Load(aSurface, L.Transformation())`
    /// — i.e. the transformed `BRep_Tool::Surface(F, L)`.
    surface: Option<Arc<dyn Surface>>,
    /// `myUFirst` = `Bounds(U1)`.
    first_u: f64,
    /// `myVFirst` = `Bounds(V1)`.
    first_v: f64,
    /// `myULast` = `Bounds(U2)` — the surface's own upper U bound for the
    /// unrestricted `Load`, the face's `UVBounds` UMax for the restricted one
    /// (`BRepAdaptor_Surface.cxx:70-80`).
    last_u: f64,
    /// `myVLast` = `Bounds(V2)`, see [`Self::last_u`].
    last_v: f64,
    /// `mySurfaceType` (`GeomAdaptor_Surface.cxx:402-431`).
    kind: GeomAbsSurfaceType,
    /// `BRepAdaptor_Surface::Face()` (`BRepAdaptor_Surface.hxx:44`): the face the
    /// adaptor was initialized with. `BRepAdaptor_Surface::Initialize(F, ...)`
    /// stores `myFace(F)` (`BRepAdaptor_Surface.cxx:56-82`).
    face: Option<Face>,
}

impl Default for SurfaceAdaptor {
    /// `BRepAdaptor_Surface()` (`BRepAdaptor_Surface.cxx:28`).
    fn default() -> Self {
        Self::new()
    }
}

impl SurfaceAdaptor {
    /// `BRepAdaptor_Surface()` (`BRepAdaptor_Surface.cxx:28`, declared at
    /// `hxx:52-53` as "an undefined surface with no face loaded"): the
    /// `= default` constructor leaves `myFace` null and the inner
    /// `GeomAdaptor_Surface` unloaded.
    pub fn new() -> Self {
        Self {
            surface: None,
            first_u: 0.0,
            first_v: 0.0,
            last_u: 0.0,
            last_v: 0.0,
            kind: GeomAbsSurfaceType::OtherSurface,
            face: None,
        }
    }

    /// `BRepAdaptor_Surface(F, Restriction)`.
    ///
    /// `of_face` is the `Restriction == false` arm
    /// (`BRepAdaptor_Surface.cxx:77-80`): `Load(aSurface, L.Transformation())`,
    /// i.e. the surface's own parameter bounds.
    pub fn of_face(f: &Face) -> Self {
        let surface = BRepTool::face_surface_world(f);
        match surface {
            Some(s) => {
                let (u1, u2) = s.u_range();
                let (v1, v2) = s.v_range();
                let kind = surface_kind(s.as_ref());
                Self {
                    surface: Some(s),
                    first_u: u1,
                    first_v: v1,
                    last_u: u2,
                    last_v: v2,
                    kind,
                    face: Some(f.clone()),
                }
            }
            None => Self {
                surface: None,
                first_u: 0.0,
                first_v: 0.0,
                last_u: 0.0,
                last_v: 0.0,
                kind: GeomAbsSurfaceType::OtherSurface,
                face: Some(f.clone()),
            },
        }
    }

    /// `BRepAdaptor_Surface(F, Restriction = true)`
    /// (`BRepAdaptor_Surface.cxx:70-76`): `BRepTools::UVBounds(F, umin, umax,
    /// vmin, vmax)` then `Load(aSurface, umin, umax, vmin, vmax,
    /// L.Transformation())`.
    ///
    /// The bound overload only replaces the four parameters; the surface
    /// itself stays untrimmed (`GeomAdaptor_Surface::Load(S, U1, U2, V1, V2,
    /// Trsf)` `GeomAdaptor_Surface.hxx`), so the geometry handed out by
    /// [`Self::surface`] is the same as [`Self::of_face`]'s.
    pub fn of_face_restricted(f: &Face) -> Self {
        let surface = BRepTool::face_surface_world(f);
        let (umin, umax, vmin, vmax) = BRepTool::uv_bounds(f);
        match surface {
            Some(s) => {
                let kind = surface_kind(s.as_ref());
                Self {
                    surface: Some(s),
                    first_u: umin,
                    first_v: vmin,
                    last_u: umax,
                    last_v: vmax,
                    kind,
                    face: Some(f.clone()),
                }
            }
            None => Self {
                surface: None,
                first_u: umin,
                first_v: vmin,
                last_u: umax,
                last_v: vmax,
                kind: GeomAbsSurfaceType::OtherSurface,
                face: Some(f.clone()),
            },
        }
    }

    /// `GeomAdaptor_Surface::Load(S)` (`GeomAdaptor_Surface.cxx:402-431`, the
    /// `Load(theSurface)` arm): the surface with its own parameter bounds.
    pub fn of_surface(s: Arc<dyn Surface>) -> Self {
        let (u1, u2) = s.u_range();
        let (v1, v2) = s.v_range();
        let kind = surface_kind(s.as_ref());
        Self {
            surface: Some(s),
            first_u: u1,
            first_v: v1,
            last_u: u2,
            last_v: v2,
            kind,
            face: None,
        }
    }

    /// `GeomAdaptor_Surface::Load(S, U1, U2, V1, V2)` (`GeomAdaptor_Surface.hxx`,
    /// the four-parameter `Load`): the same surface, restricted to the given
    /// parameter box. `ChFi3d_BoundSurf` builds one of these to keep the
    /// intersection inside a range around the fillet.
    pub fn of_surface_with_range(
        s: Arc<dyn Surface>,
        first_u: f64,
        last_u: f64,
        first_v: f64,
        last_v: f64,
    ) -> Self {
        let kind = surface_kind(s.as_ref());
        Self {
            surface: Some(s),
            first_u,
            first_v,
            last_u,
            last_v,
            kind,
            face: None,
        }
    }

    /// `Adaptor3d_Surface::Value(U, V)` — `GeomAdaptor_Surface::Value` forwards
    /// to the loaded surface (`GeomAdaptor_Surface.cxx`).
    pub fn value(&self, u: f64, v: f64) -> Option<GpPnt> {
        self.surface.as_ref().map(|s| s.d0(u, v))
    }

    /// `BRepAdaptor_Surface::Face()` (`BRepAdaptor_Surface.hxx:44`).
    pub fn face(&self) -> Option<&Face> {
        self.face.as_ref()
    }

    /// `Adaptor3d_Surface::Surface()`.
    pub fn surface(&self) -> Option<&Arc<dyn Surface>> {
        self.surface.as_ref()
    }

    /// `GeomAdaptor_Surface::FirstUParameter()` (`GeomAdaptor_Surface.hxx:185`).
    pub fn first_u_parameter(&self) -> f64 {
        self.first_u
    }

    /// `GeomAdaptor_Surface::FirstVParameter()` (`GeomAdaptor_Surface.hxx:189`).
    pub fn first_v_parameter(&self) -> f64 {
        self.first_v
    }

    /// `GeomAdaptor_Surface::LastUParameter()` (`GeomAdaptor_Surface.hxx:187`):
    /// `myULast`, i.e. the surface's own upper U bound for [`Self::of_face`] and
    /// the face's `UVBounds` UMax for [`Self::of_face_restricted`].
    pub fn last_u_parameter(&self) -> f64 {
        self.last_u
    }

    /// `GeomAdaptor_Surface::LastVParameter()` (`GeomAdaptor_Surface.hxx:191`).
    pub fn last_v_parameter(&self) -> f64 {
        self.last_v
    }

    /// `Adaptor3d_Surface::GetType()`.
    pub fn get_type(&self) -> GeomAbsSurfaceType {
        self.kind
    }

    /// `Adaptor3d_Surface::Plane()` (`GeomAdaptor_Surface::Plane`,
    /// `GeomAdaptor_Surface.cxx:172-176`): the `gp_Pln` of a loaded
    /// `Geom_Plane`. OCCT raises `Standard_NoSuchObject` when the type is not
    /// `GeomAbs_Plane`; the port answers `None`.
    pub fn plane(&self) -> Option<GpPln> {
        self.surface.as_ref().and_then(|s| s.gp_pln())
    }

    /// `Adaptor3d_Surface::Cylinder()` (`GeomAdaptor_Surface.cxx`): the
    /// `gp_Cylinder` of a loaded `Geom_CylindricalSurface`; `None` otherwise,
    /// where OCCT raises `Standard_NoSuchObject`.
    pub fn cylinder(&self) -> Option<GpCylinder> {
        self.surface.as_ref().and_then(|s| s.gp_cylinder())
    }

    /// `Adaptor3d_Surface::Cone()` (`GeomAdaptor_Surface.cxx`): the `gp_Cone`
    /// of a loaded `Geom_ConicalSurface`; `None` otherwise, where OCCT raises
    /// `Standard_NoSuchObject`.
    pub fn cone(&self) -> Option<GpCone> {
        self.surface.as_ref().and_then(|s| s.gp_cone())
    }

    /// `Adaptor3d_Surface::D1(U, V, P, D1U, D1V)`.
    pub fn d1(&self, u: f64, v: f64) -> Option<(GpPnt, GpVec, GpVec)> {
        self.surface.as_ref().map(|s| s.d1(u, v))
    }

    /// `Adaptor3d_Surface::Value(U, V, P)` = `D0`, which
    /// `GeomAdaptor_Surface::D0` forwards to `mySurf->D0`.
    pub fn d0(&self, u: f64, v: f64) -> Option<GpPnt> {
        self.surface.as_ref().map(|s| s.d0(u, v))
    }

    /// `Adaptor3d_Surface::D2(U, V, P, D1U, D1V, D2U, D2V, D2UV)`:
    /// `GeomAdaptor_Surface::D2` forwards to `mySurf->D2` for every untrimmed
    /// surface, so this is [`occt_geom::Surface::d2`] on the loaded surface.
    pub fn d2(&self, u: f64, v: f64) -> Option<(GpPnt, GpVec, GpVec, GpVec, GpVec, GpVec)> {
        self.surface.as_ref().map(|s| s.d2(u, v))
    }

    /// `BRepAdaptor_Surface::Initialize(F)` (`BRepAdaptor_Surface.cxx:56-82`):
    /// `Load(BRep_Tool::Surface(F, L), L.Transformation())`, keeping the new
    /// face. `of_face` is the constructor form of the same state.
    pub fn initialize(&mut self, f: &Face) {
        *self = Self::of_face(f);
    }
}

/// `BRepAdaptor_Curve(E)` (`BRepAdaptor_Curve.hxx:36-98`).
#[derive(Clone)]
pub struct CurveAdaptor {
    /// `GeomAdaptor_Curve::Curve()` after `Load(BRep_Tool::Curve(E, L))`.
    curve: Option<Arc<dyn Curve>>,
    /// `myFirst` (`GeomAdaptor_Curve.cxx:243`).
    first: f64,
    /// `myLast`.
    last: f64,
    /// `myTypeCurve` (`GeomAdaptor_Curve.cxx:238-313`).
    kind: GeomAbsCurveType,
    /// `BRepAdaptor_Curve::Tolerance()` = `BRep_Tool::Tolerance(E)`.
    tolerance: f64,
}

/// `RwLock`-free view of `BRepAdaptor_Curve`'s members: the inner `dyn Curve`
/// carries no `Debug`, so it prints as loaded / not loaded.
impl std::fmt::Debug for CurveAdaptor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CurveAdaptor")
            .field("has_curve", &self.curve.is_some())
            .field("first", &self.first)
            .field("last", &self.last)
            .field("kind", &self.kind)
            .field("tolerance", &self.tolerance)
            .finish()
    }
}

impl CurveAdaptor {
    /// `BRepAdaptor_Curve(const TopoDS_Edge&)` — `Initialize(E)` with no face,
    /// so the 3D curve is used when the edge has one and the pcurve adaptor
    /// otherwise.
    pub fn of_edge(e: &Edge) -> Self {
        let curve = BRepTool::edge_curve_adaptor_world(e);
        let (first, last) = BRepTool::edge_parameters(e);
        let kind = match curve.as_ref() {
            Some(c) => curve_kind(c.as_ref()),
            None => GeomAbsCurveType::OtherCurve,
        };
        let tolerance = BRepTool::edge_tolerance(e);
        Self {
            curve,
            first,
            last,
            kind,
            tolerance,
        }
    }

    /// `Adaptor3d_Curve::Curve()`.
    pub fn curve(&self) -> Option<&Arc<dyn Curve>> {
        self.curve.as_ref()
    }

    /// `Adaptor3d_Curve::FirstParameter()`.
    pub fn first_parameter(&self) -> f64 {
        self.first
    }

    /// `Adaptor3d_Curve::LastParameter()`.
    pub fn last_parameter(&self) -> f64 {
        self.last
    }

    /// `Adaptor3d_Curve::GetType()`.
    pub fn get_type(&self) -> GeomAbsCurveType {
        self.kind
    }

    /// `Adaptor3d_Curve::Continuity()` — needed by
    /// `LProp_CurveUtils::IsTangentDefined` (`LProp_CurveUtils.hxx:332-378`).
    pub fn continuity(&self) -> u8 {
        match self.curve.as_ref() {
            Some(c) => c.continuity(),
            None => 0,
        }
    }

    /// `Adaptor3d_Curve::D0(U, P)`.
    pub fn d0(&self, u: f64) -> Option<GpPnt> {
        self.curve.as_ref().map(|c| c.d0(u))
    }

    /// `Adaptor3d_Curve::D1(U, P, D1)`.
    pub fn d1(&self, u: f64) -> Option<(GpPnt, GpVec)> {
        self.curve.as_ref().map(|c| c.d1(u))
    }

    /// `Adaptor3d_Curve::D2(U, P, D1, D2)`.
    pub fn d2(&self, u: f64) -> Option<(GpPnt, GpVec, GpVec)> {
        self.curve.as_ref().map(|c| c.d2(u))
    }

    /// `Adaptor3d_Curve::D3(U, P, D1, D2, D3)`.
    pub fn d3(&self, u: f64) -> Option<(GpPnt, GpVec, GpVec, GpVec)> {
        self.curve.as_ref().map(|c| c.d3(u))
    }

    /// `Adaptor3d_Curve::Resolution(R3d)` — `GeomAdaptor_Curve::Resolution`
    /// (`GeomAdaptor_Curve.cxx:1116-1148`), i.e. the underlying curve's.
    pub fn resolution(&self, r3d: f64) -> f64 {
        match self.curve.as_ref() {
            Some(c) => c.resolution(r3d),
            None => r3d * occt_core::precision::PCONFUSION / occt_core::precision::CONFUSION,
        }
    }

    /// `BRepAdaptor_Curve::Tolerance()` = `BRep_Tool::Tolerance(E)`.
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// `BRepAdaptor_Curve::Trim(First, Last, Tol)` (`BRepAdaptor_Curve.cxx:464-471`)
    /// clones the adaptor and calls
    /// `GeomAdaptor_Curve::Trim(First, Last, Tol)` (`GeomAdaptor_Curve.cxx:302-315`),
    /// which only overwrites `myFirst` / `myLast`: evaluation keeps the basis
    /// parameters, so the inner curve is untouched.
    pub fn trim(&self, first: f64, last: f64) -> CurveAdaptor {
        let mut trimmed = self.clone();
        trimmed.first = first;
        trimmed.last = last;
        trimmed
    }

    /// `Adaptor3d_Curve::IsClosed()` — `GeomAdaptor_Curve::IsClosed`
    /// (`GeomAdaptor_Curve.cxx:576-585`): an infinite range is never closed,
    /// otherwise the two ends must coincide within `Precision::Confusion()`.
    pub fn is_closed(&self) -> bool {
        let Some(curve) = self.curve.as_ref() else {
            return false;
        };
        if occt_core::precision::Precision::is_infinite(self.last)
            || occt_core::precision::Precision::is_infinite(self.first)
        {
            return false;
        }
        curve.d0(self.first).distance(&curve.d0(self.last)) <= occt_core::precision::CONFUSION
    }
}

/// `BRepAdaptor_Curve` *is* an `Adaptor3d_Curve` (`BRepAdaptor_Curve.hxx:36-98`),
/// so this is the trait view for algorithms that take one. `ChFiDS_Spine`
/// needs it for its `GCPnts_AbscissaPoint` calls on `myCurve`
/// (`ChFiDS_Spine.cxx:556-561`, `:609`, `:737`, `:786`, `:836`): those read the
/// *edge* range from the adaptor (`FirstParameter()` / `LastParameter()` store
/// `myFirst` / `myLast`, `GeomAdaptor_Curve.cxx:243-245`), not the inner
/// curve's own range, which is why the range is forwarded from `self.first` /
/// `self.last` here.
///
/// A null inner curve makes every geometric query raise `Standard_NullObject`
/// in OCCT; the trait cannot, so those unwrap with that message.
impl Curve for CurveAdaptor {
    fn d0(&self, u: f64) -> GpPnt {
        self.geom_curve().d0(u)
    }

    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        self.geom_curve().d1(u)
    }

    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        self.geom_curve().d2(u)
    }

    fn d3(&self, u: f64) -> (GpPnt, GpVec, GpVec, GpVec) {
        self.geom_curve().d3(u)
    }

    fn first_parameter(&self) -> f64 {
        self.first
    }

    fn last_parameter(&self) -> f64 {
        self.last
    }

    fn continuity(&self) -> u8 {
        self.geom_curve().continuity()
    }

    fn nb_intervals(&self, continuity: u8) -> i32 {
        match self.curve.as_ref() {
            Some(c) => c.nb_intervals(continuity),
            None => 1,
        }
    }

    fn parameter_intervals(&self, continuity: u8) -> Vec<f64> {
        match self.curve.as_ref() {
            Some(c) => c.parameter_intervals(continuity),
            None => vec![self.first, self.last],
        }
    }

    fn resolution(&self, r3d: f64) -> f64 {
        CurveAdaptor::resolution(self, r3d)
    }

    fn is_periodic(&self) -> bool {
        match self.curve.as_ref() {
            Some(c) => c.is_periodic(),
            None => false,
        }
    }

    fn period(&self) -> f64 {
        match self.curve.as_ref() {
            Some(c) => c.period(),
            None => 0.0,
        }
    }

    fn circle_radius(&self) -> Option<f64> {
        self.curve.as_ref().and_then(|c| c.circle_radius())
    }

    fn gp_line(&self) -> Option<GpLin> {
        self.curve.as_ref().and_then(|c| c.gp_line())
    }

    fn gp_circ(&self) -> Option<GpCirc> {
        self.curve.as_ref().and_then(|c| c.gp_circ())
    }

    fn gp_ellipse(&self) -> Option<GpElips> {
        self.curve.as_ref().and_then(|c| c.gp_ellipse())
    }

    fn gp_hyperbola(&self) -> Option<GpHypr> {
        self.curve.as_ref().and_then(|c| c.gp_hyperbola())
    }

    fn gp_parabola(&self) -> Option<GpParab> {
        self.curve.as_ref().and_then(|c| c.gp_parabola())
    }

    fn is_line(&self) -> bool {
        self.curve.as_ref().is_some_and(|c| c.is_line())
    }

    fn nurbs_degree(&self) -> Option<usize> {
        self.curve.as_ref().and_then(|c| c.nurbs_degree())
    }

    fn bspline_poles(&self) -> Option<&[GpPnt]> {
        self.curve.as_ref().and_then(|c| c.bspline_poles())
    }

    fn bspline_weights(&self) -> Option<&[f64]> {
        self.curve.as_ref().and_then(|c| c.bspline_weights())
    }

    fn bezier_poles(&self) -> Option<&[GpPnt]> {
        self.curve.as_ref().and_then(|c| c.bezier_poles())
    }

    fn bezier_weights(&self) -> Option<&[f64]> {
        self.curve.as_ref().and_then(|c| c.bezier_weights())
    }

    /// `Adaptor3d_Curve` has no `Transform`; only `Geom_Curve` does.
    fn transform(&mut self, _t: &GpTrsf) {}

    /// `Adaptor3d_Curve` has no `Reverse`; only `Geom_Curve` does.
    fn reverse(&mut self) {}

    fn clone_dyn(&self) -> Box<dyn Curve> {
        Box::new(self.clone())
    }
}

impl CurveAdaptor {
    /// The inner curve, or the `Standard_NullObject` OCCT raises when the
    /// adaptor was never initialized.
    fn geom_curve(&self) -> &dyn Curve {
        self.curve
            .as_ref()
            .map(|c| c.as_ref())
            .expect("BRepAdaptor_Curve: no curve loaded (Standard_NullObject)")
    }
}

/// `GeomAdaptor_Curve` — a `Geom_Curve` plus the range it is evaluated on.
///
/// Source: `TKG3d/GeomAdaptor/GeomAdaptor_Curve.cxx:262-271` (`Load(C, U1,
/// U2)`), `:243-245` (`myFirst` / `myLast`), `:576-585` (`IsClosed`).
///
/// `Load` keeps the *basis* parameters: `Value(U)` forwards to
/// `myCurve->Value(U)` and only `FirstParameter` / `LastParameter` report the
/// trimmed range. That is what the two `ChFi3d` callers below rely on —
/// `ChFi3d_ComputeArete` loads `GeomAdaptor_Curve(C3d, Pardeb, Parfin)` and
/// hands it to `ChFi3d_ComputePCurv`, whose p-curve is trimmed to the same
/// basis parameters.
#[derive(Clone)]
pub struct GeomAdaptorCurve {
    /// `myCurve`.
    curve: Option<Arc<dyn Curve>>,
    /// `myFirst`.
    first: f64,
    /// `myLast`.
    last: f64,
}

impl GeomAdaptorCurve {
    /// `GeomAdaptor_Curve()` (`cxx:28-33`).
    pub fn new() -> Self {
        Self {
            curve: None,
            first: 0.0,
            last: 0.0,
        }
    }

    /// `Load(C)` (`GeomAdaptor_Curve.cxx:225-231`): the curve's own range.
    pub fn of_curve(curve: Arc<dyn Curve>) -> Self {
        Self {
            first: curve.first_parameter(),
            last: curve.last_parameter(),
            curve: Some(curve),
        }
    }

    /// `Load(C, U1, U2)` (`cxx:262-271`).
    pub fn with_range(curve: Arc<dyn Curve>, first: f64, last: f64) -> Self {
        Self {
            curve: Some(curve),
            first,
            last,
        }
    }

    /// `Curve()` (`hxx`).
    pub fn curve(&self) -> Option<&Arc<dyn Curve>> {
        self.curve.as_ref()
    }

    /// `FirstParameter()` (`cxx:243-245`).
    pub fn first_parameter(&self) -> f64 {
        self.first
    }

    /// `LastParameter()` (`cxx:243-245`).
    pub fn last_parameter(&self) -> f64 {
        self.last
    }

    /// `Value(U)`.
    pub fn value(&self, u: f64) -> GpPnt {
        self.geom_curve().d0(u)
    }

    /// See [`CurveAdaptor::geom_curve`]: a null `myCurve` is
    /// `Standard_NullObject`.
    fn geom_curve(&self) -> &dyn Curve {
        self.curve
            .as_ref()
            .map(|c| c.as_ref())
            .expect("GeomAdaptor_Curve: no curve loaded (Standard_NullObject)")
    }
}

impl Default for GeomAdaptorCurve {
    fn default() -> Self {
        Self::new()
    }
}

/// `Adaptor3d_Curve` view of [`GeomAdaptorCurve`], the same way
/// [`CurveAdaptor`] is the view of `BRepAdaptor_Curve`.
impl Curve for GeomAdaptorCurve {
    fn d0(&self, u: f64) -> GpPnt {
        self.geom_curve().d0(u)
    }

    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        self.geom_curve().d1(u)
    }

    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        self.geom_curve().d2(u)
    }

    fn d3(&self, u: f64) -> (GpPnt, GpVec, GpVec, GpVec) {
        self.geom_curve().d3(u)
    }

    fn first_parameter(&self) -> f64 {
        self.first
    }

    fn last_parameter(&self) -> f64 {
        self.last
    }

    fn continuity(&self) -> u8 {
        self.geom_curve().continuity()
    }

    fn nb_intervals(&self, continuity: u8) -> i32 {
        match self.curve.as_ref() {
            Some(c) => c.nb_intervals(continuity),
            None => 1,
        }
    }

    fn parameter_intervals(&self, continuity: u8) -> Vec<f64> {
        match self.curve.as_ref() {
            Some(c) => c.parameter_intervals(continuity),
            None => vec![self.first, self.last],
        }
    }

    fn resolution(&self, r3d: f64) -> f64 {
        match self.curve.as_ref() {
            Some(c) => c.resolution(r3d),
            None => r3d * occt_core::precision::PCONFUSION / occt_core::precision::CONFUSION,
        }
    }

    fn is_periodic(&self) -> bool {
        self.curve.as_ref().is_some_and(|c| c.is_periodic())
    }

    fn period(&self) -> f64 {
        self.curve.as_ref().map(|c| c.period()).unwrap_or(0.0)
    }

    fn circle_radius(&self) -> Option<f64> {
        self.curve.as_ref().and_then(|c| c.circle_radius())
    }

    fn is_line(&self) -> bool {
        self.curve.as_ref().is_some_and(|c| c.is_line())
    }

    fn gp_line(&self) -> Option<GpLin> {
        self.curve.as_ref().and_then(|c| c.gp_line())
    }

    fn gp_circ(&self) -> Option<GpCirc> {
        self.curve.as_ref().and_then(|c| c.gp_circ())
    }

    fn nurbs_degree(&self) -> Option<usize> {
        self.curve.as_ref().and_then(|c| c.nurbs_degree())
    }

    fn bspline_poles(&self) -> Option<&[GpPnt]> {
        self.curve.as_ref().and_then(|c| c.bspline_poles())
    }

    fn bspline_weights(&self) -> Option<&[f64]> {
        self.curve.as_ref().and_then(|c| c.bspline_weights())
    }

    fn bezier_poles(&self) -> Option<&[GpPnt]> {
        self.curve.as_ref().and_then(|c| c.bezier_poles())
    }

    /// `Adaptor3d_Curve` has no `Transform`; only `Geom_Curve` does.
    fn transform(&mut self, _t: &GpTrsf) {}

    /// `Adaptor3d_Curve` has no `Reverse`; only `Geom_Curve` does.
    fn reverse(&mut self) {}

    fn clone_dyn(&self) -> Box<dyn Curve> {
        Box::new(self.clone())
    }
}
