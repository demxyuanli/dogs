use std::sync::Arc;

use occt_core::gp::{GpCone, GpCylinder, GpPln, GpPnt, GpSphere, GpTorus, GpVec};
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;

use crate::brep_tool::BRepTool;
use crate::shape::Face;
use occt_geom::geom_bnd_lib::surface3d::{surface_kind, GeomAbsSurfaceType};

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

    /// `Adaptor3d_Surface::Sphere()` (`GeomAdaptor_Surface::Sphere`): the
    /// `gp_Sphere` of a loaded `Geom_SphericalSurface`; `None` otherwise.
    pub fn sphere(&self) -> Option<GpSphere> {
        self.surface.as_ref().and_then(|s| s.gp_sphere())
    }

    /// `Adaptor3d_Surface::Torus()` (`GeomAdaptor_Surface::Torus`): the
    /// `gp_Torus` of a loaded `Geom_ToroidalSurface`; `None` otherwise.
    pub fn torus(&self) -> Option<GpTorus> {
        self.surface.as_ref().and_then(|s| s.gp_torus())
    }

    /// `Adaptor3d_Surface::UResolution(Ruv)` (`GeomAdaptor_Surface.cxx:1818-1895`),
    /// the same helper `IntTools_WLineTool` uses.
    pub fn u_resolution(&self, ruv: f64) -> f64 {
        match self.surface.as_ref() {
            Some(s) => occt_geom::approx_same_parameter::u_resolution(s.as_ref(), ruv),
            None => occt_core::precision::CONFUSION,
        }
    }

    /// `Adaptor3d_Surface::VResolution(Ruv)` (`GeomAdaptor_Surface.cxx:1900-1957`).
    pub fn v_resolution(&self, ruv: f64) -> f64 {
        match self.surface.as_ref() {
            Some(s) => occt_geom::approx_same_parameter::v_resolution(s.as_ref(), ruv),
            None => occt_core::precision::CONFUSION,
        }
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

    /// `Adaptor3d_Surface::UDegree()` (`GeomAdaptor_Surface.cxx:2015-2032`):
    /// BSpline / Bezier take the surface degree, a surface of extrusion takes
    /// its basis curve's degree, everything else throws
    /// `Standard_NoSuchObject` -> `None`.
    pub fn u_degree(&self) -> Option<i32> {
        let s = self.surface.as_ref()?;
        match self.kind {
            GeomAbsSurfaceType::BSplineSurface | GeomAbsSurfaceType::BezierSurface => {
                Some(s.u_degree())
            }
            GeomAbsSurfaceType::SurfaceOfExtrusion => {
                s.extrusion_basis_curve().and_then(|c| curve_degree(c.as_ref()))
            }
            _ => None,
        }
    }

    /// `Adaptor3d_Surface::VDegree()` (`GeomAdaptor_Surface.cxx:2055-2070`).
    pub fn v_degree(&self) -> Option<i32> {
        let s = self.surface.as_ref()?;
        match self.kind {
            GeomAbsSurfaceType::BSplineSurface | GeomAbsSurfaceType::BezierSurface => {
                Some(s.v_degree())
            }
            GeomAbsSurfaceType::SurfaceOfRevolution => {
                s.revolution_basis_curve().and_then(|c| curve_degree(c.as_ref()))
            }
            _ => None,
        }
    }

    /// `Adaptor3d_Surface::NbUPoles()` (`GeomAdaptor_Surface.cxx:2036-2051`).
    pub fn nb_u_poles(&self) -> Option<i32> {
        let s = self.surface.as_ref()?;
        match self.kind {
            GeomAbsSurfaceType::BSplineSurface | GeomAbsSurfaceType::BezierSurface => {
                Some(s.nb_u_poles())
            }
            GeomAbsSurfaceType::SurfaceOfExtrusion => {
                s.extrusion_basis_curve().and_then(|c| curve_nb_poles(c.as_ref()))
            }
            _ => None,
        }
    }

    /// `Adaptor3d_Surface::NbVPoles()` (`GeomAdaptor_Surface.cxx:2076-2091`).
    pub fn nb_v_poles(&self) -> Option<i32> {
        let s = self.surface.as_ref()?;
        match self.kind {
            GeomAbsSurfaceType::BSplineSurface | GeomAbsSurfaceType::BezierSurface => {
                Some(s.nb_v_poles())
            }
            GeomAbsSurfaceType::SurfaceOfRevolution => {
                s.revolution_basis_curve().and_then(|c| curve_nb_poles(c.as_ref()))
            }
            _ => None,
        }
    }

    /// `Adaptor3d_Surface::NbUKnots()` (`GeomAdaptor_Surface.cxx:2095-2105`).
    pub fn nb_u_knots(&self) -> Option<i32> {
        let s = self.surface.as_ref()?;
        match self.kind {
            GeomAbsSurfaceType::BSplineSurface => {
                distinct_knot_count(s.bspline_surface_uknots()?)
            }
            GeomAbsSurfaceType::SurfaceOfExtrusion => {
                s.extrusion_basis_curve().and_then(|c| curve_nb_knots(c.as_ref()))
            }
            _ => None,
        }
    }

    /// `Adaptor3d_Surface::NbVKnots()` (`GeomAdaptor_Surface.cxx:2111-2118`).
    pub fn nb_v_knots(&self) -> Option<i32> {
        let s = self.surface.as_ref()?;
        match self.kind {
            GeomAbsSurfaceType::BSplineSurface => {
                distinct_knot_count(s.bspline_surface_vknots()?)
            }
            _ => None,
        }
    }

    /// `BRepAdaptor_Surface::Initialize(F)` (`BRepAdaptor_Surface.cxx:56-82`):
    /// `Load(BRep_Tool::Surface(F, L), L.Transformation())`, keeping the new
    /// face. `of_face` is the constructor form of the same state.
    pub fn initialize(&mut self, f: &Face) {
        *self = Self::of_face(f);
    }
}

/// `GeomAdaptor_Curve::Degree()` (`GeomAdaptor_Curve.cxx:1203-1218`).
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
/// BSpline curve answers.
fn curve_nb_knots(c: &dyn Curve) -> Option<i32> {
    distinct_knot_count(c.bspline_knots()?)
}

/// `Geom_BSplineSurface::NbUKnots()` / a flat knot vector's distinct values.
fn distinct_knot_count(knots: &[f64]) -> Option<i32> {
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
