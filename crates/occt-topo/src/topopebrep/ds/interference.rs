//! `TopOpeBRepDS_Interference` and its subclasses — the attachment of a new
//! geometry on a geometry (an intersection point on an edge or on a curve).
//!
//! Sources:
//! - `TopOpeBRepDS_Interference.hxx:29-104`, `.cxx:20-165`
//! - `TopOpeBRepDS_SurfaceCurveInterference.hxx:28-56`, `.cxx:22-55`
//! - `TopOpeBRepDS_CurvePointInterference.hxx:28-45`, `.cxx:22-47`
//! - `TopOpeBRepDS_SolidSurfaceInterference.hxx:26-40`, `.cxx:22-33`
//! - `TopOpeBRepDS_ShapeShapeInterference.hxx:29-70`, `.cxx:24-54`
//! - `TopOpeBRepDS_EdgeVertexInterference.hxx:28-84`, `.cxx:24-60`
//! - `TopOpeBRepDS_FaceEdgeInterference.hxx:26-45`, `.cxx:22-37`
//!
//! OCCT spells the hierarchy with `Standard_Transient` RTTI and keeps
//! `NCollection_List<occ::handle<TopOpeBRepDS_Interference>>`, i.e. one
//! homogeneous list of base handles. The Rust translation keeps one
//! [`Interference`] value with an [`InterferenceData`] payload that carries the
//! subclass fields, so a list stays homogeneous too. `IsKind(STANDARD_TYPE(X))`
//! is [`Interference::is_shape_shape`] / `as_*`, and `DynamicType()` is
//! [`Interference::dynamic_type`].

use std::sync::{Arc, RwLock};

use occt_geom2d::curve::Curve2d;

use super::kind::{Config, Kind};
use super::transition::Transition;

/// `occ::handle<TopOpeBRepDS_Interference>`.
///
/// OCCT hands out the same object by reference (`ChangeTransition()` mutates
/// what every holder sees), so the Rust handle is a shared, lockable cell and
/// not an `Arc` of an immutable value.
pub type InterferenceHandle = Arc<RwLock<Interference>>;

/// The subclass fields of `TopOpeBRepDS_ShapeShapeInterference`
/// (`ShapeShapeInterference.hxx:68-69`).
#[derive(Debug, Clone, Copy, Default)]
pub struct ShapeShapeData {
    /// `myGBound`
    pub g_bound: bool,
    /// `myC`
    pub config: Config,
}

/// `TopOpeBRepDS_SurfaceCurveInterference` (`SurfaceCurveInterference.hxx:53`).
#[derive(Clone, Default)]
pub struct SurfaceCurveData {
    /// `myPCurve`
    pub p_curve: Option<Arc<dyn Curve2d>>,
}

impl std::fmt::Debug for SurfaceCurveData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SurfaceCurveData")
            .field("is_null", &self.p_curve.is_none())
            .finish()
    }
}

/// `TopOpeBRepDS_EdgeVertexInterference` (`EdgeVertexInterference.hxx:81`):
/// a `ShapeShapeInterference` plus the parameter of the vertex on the edge.
#[derive(Debug, Clone, Copy, Default)]
pub struct EdgeVertexData {
    /// The `TopOpeBRepDS_ShapeShapeInterference` part.
    pub shape_shape: ShapeShapeData,
    /// `myParam`
    pub param: f64,
}

/// The subclass payload of an [`Interference`]. One variant per
/// `TopOpeBRepDS_Interference` subclass, plus [`InterferenceData::Base`] for
/// the class itself.
#[derive(Debug, Clone, Default)]
pub enum InterferenceData {
    /// `TopOpeBRepDS_Interference`
    #[default]
    Base,
    /// `TopOpeBRepDS_SurfaceCurveInterference`
    SurfaceCurve(SurfaceCurveData),
    /// `TopOpeBRepDS_CurvePointInterference` (`myParam`)
    CurvePoint { param: f64 },
    /// `TopOpeBRepDS_SolidSurfaceInterference`
    SolidSurface,
    /// `TopOpeBRepDS_ShapeShapeInterference`
    ShapeShape(ShapeShapeData),
    /// `TopOpeBRepDS_EdgeVertexInterference`
    EdgeVertex(EdgeVertexData),
    /// `TopOpeBRepDS_FaceEdgeInterference`
    FaceEdge(ShapeShapeData),
}

/// Source: `TopOpeBRepDS_Interference.hxx:96-101`
#[derive(Debug, Clone)]
pub struct Interference {
    transition: Transition,
    support: i32,
    geometry: i32,
    support_type: Kind,
    geometry_type: Kind,
    data: InterferenceData,
}

impl Default for Interference {
    /// `TopOpeBRepDS_Interference()` (`cxx:26`): the OCCT constructor is
    /// `= default`, leaving the scalar members indeterminate. The Rust default
    /// picks the zero values of the same fields so a constructed base is
    /// deterministic.
    fn default() -> Self {
        Self {
            transition: Transition::default(),
            support: 0,
            geometry: 0,
            support_type: Kind::Point,
            geometry_type: Kind::Point,
            data: InterferenceData::Base,
        }
    }
}

impl Interference {
    /// `TopOpeBRepDS_Interference(T, ST, S, GT, G)` (`cxx:28-39`).
    pub fn new(
        transition: Transition,
        support_type: Kind,
        support: i32,
        geometry_type: Kind,
        geometry: i32,
    ) -> Self {
        Self {
            transition,
            support,
            geometry,
            support_type,
            geometry_type,
            data: InterferenceData::Base,
        }
    }

    /// The base class part of
    /// `TopOpeBRepDS_SurfaceCurveInterference(T, ST, S, GT, G, PC)`
    /// (`SurfaceCurveInterference.cxx:28-37`).
    pub fn surface_curve(
        transition: Transition,
        support_type: Kind,
        support: i32,
        geometry_type: Kind,
        geometry: i32,
        p_curve: Option<Arc<dyn Curve2d>>,
    ) -> Self {
        Self {
            transition,
            support,
            geometry,
            support_type,
            geometry_type,
            data: InterferenceData::SurfaceCurve(SurfaceCurveData { p_curve }),
        }
    }

    /// `TopOpeBRepDS_SurfaceCurveInterference(const handle<Interference>& I)`
    /// (`SurfaceCurveInterference.cxx:42-46`): the base part is copied, the
    /// 2d curve is the null handle.
    pub fn surface_curve_from_interference(i: &Interference) -> Self {
        Self {
            data: InterferenceData::SurfaceCurve(SurfaceCurveData::default()),
            ..i.clone()
        }
    }

    /// `TopOpeBRepDS_CurvePointInterference(T, ST, S, GT, G, P)`
    /// (`CurvePointInterference.cxx:24-35`).
    pub fn curve_point(
        transition: Transition,
        support_type: Kind,
        support: i32,
        geometry_type: Kind,
        geometry: i32,
        param: f64,
    ) -> Self {
        Self {
            transition,
            support,
            geometry,
            support_type,
            geometry_type,
            data: InterferenceData::CurvePoint { param },
        }
    }

    /// `TopOpeBRepDS_SolidSurfaceInterference(T, ST, S, GT, G)`
    /// (`SolidSurfaceInterference.cxx:24-31`): no subclass field.
    pub fn solid_surface(
        transition: Transition,
        support_type: Kind,
        support: i32,
        geometry_type: Kind,
        geometry: i32,
    ) -> Self {
        Self {
            transition,
            support,
            geometry,
            support_type,
            geometry_type,
            data: InterferenceData::SolidSurface,
        }
    }

    /// `TopOpeBRepDS_ShapeShapeInterference(T, ST, S, GT, G, GBound, C)`
    /// (`ShapeShapeInterference.cxx:28-39`).
    pub fn shape_shape(
        transition: Transition,
        support_type: Kind,
        support: i32,
        geometry_type: Kind,
        geometry: i32,
        g_bound: bool,
        config: Config,
    ) -> Self {
        Self {
            transition,
            support,
            geometry,
            support_type,
            geometry_type,
            data: InterferenceData::ShapeShape(ShapeShapeData { g_bound, config }),
        }
    }

    /// `TopOpeBRepDS_EdgeVertexInterference(T, ST, S, G, GIsBound, C, P)`
    /// (`EdgeVertexInterference.cxx:28-40`): geometry type is `VERTEX`.
    pub fn edge_vertex(
        transition: Transition,
        support_type: Kind,
        support: i32,
        geometry: i32,
        g_is_bound: bool,
        config: Config,
        param: f64,
    ) -> Self {
        Self {
            transition,
            support,
            geometry,
            support_type,
            geometry_type: Kind::Vertex,
            data: InterferenceData::EdgeVertex(EdgeVertexData {
                shape_shape: ShapeShapeData {
                    g_bound: g_is_bound,
                    config,
                },
                param,
            }),
        }
    }

    /// `TopOpeBRepDS_EdgeVertexInterference(T, S, G, GIsBound, C, P)`
    /// (`EdgeVertexInterference.cxx:45-55`): support type is `EDGE`.
    pub fn edge_vertex_on_edge(
        transition: Transition,
        support: i32,
        geometry: i32,
        g_is_bound: bool,
        config: Config,
        param: f64,
    ) -> Self {
        Self::edge_vertex(
            transition,
            Kind::Edge,
            support,
            geometry,
            g_is_bound,
            config,
            param,
        )
    }

    /// `TopOpeBRepDS_FaceEdgeInterference(T, S, G, GIsBound, C)`
    /// (`FaceEdgeInterference.cxx:22-37`): support type `FACE`, geometry type
    /// `EDGE`.
    pub fn face_edge(
        transition: Transition,
        support: i32,
        geometry: i32,
        g_is_bound: bool,
        config: Config,
    ) -> Self {
        Self {
            transition,
            support,
            geometry,
            support_type: Kind::Face,
            geometry_type: Kind::Edge,
            data: InterferenceData::FaceEdge(ShapeShapeData {
                g_bound: g_is_bound,
                config,
            }),
        }
    }

    /// `TopOpeBRepDS_Interference(const handle<Interference>& I)` (`cxx:44-51`),
    /// used by the subclass copy constructors.
    pub fn from_interference(i: &Interference) -> Self {
        Self {
            transition: i.transition(),
            support: i.support(),
            geometry: i.geometry(),
            support_type: i.support_type(),
            geometry_type: i.geometry_type(),
            data: i.data.clone(),
        }
    }

    /// `DynamicType().Name()` of the represented class.
    pub fn dynamic_type(&self) -> &'static str {
        match self.data {
            InterferenceData::Base => "TopOpeBRepDS_Interference",
            InterferenceData::SurfaceCurve(_) => "TopOpeBRepDS_SurfaceCurveInterference",
            InterferenceData::CurvePoint { .. } => "TopOpeBRepDS_CurvePointInterference",
            InterferenceData::SolidSurface => "TopOpeBRepDS_SolidSurfaceInterference",
            InterferenceData::ShapeShape(_) => "TopOpeBRepDS_ShapeShapeInterference",
            InterferenceData::EdgeVertex(_) => "TopOpeBRepDS_EdgeVertexInterference",
            InterferenceData::FaceEdge(_) => "TopOpeBRepDS_FaceEdgeInterference",
        }
    }

    /// `IsKind(STANDARD_TYPE(TopOpeBRepDS_ShapeShapeInterference))`: true for
    /// `ShapeShape`, `EdgeVertex` and `FaceEdge`.
    pub fn is_shape_shape(&self) -> bool {
        matches!(
            self.data,
            InterferenceData::ShapeShape(_)
                | InterferenceData::EdgeVertex(_)
                | InterferenceData::FaceEdge(_)
        )
    }

    /// `down_cast<TopOpeBRepDS_ShapeShapeInterference>` data, i.e. the
    /// `ShapeShapeData` of any of the three subclass variants.
    pub fn shape_shape_data(&self) -> Option<&ShapeShapeData> {
        match &self.data {
            InterferenceData::ShapeShape(d) | InterferenceData::FaceEdge(d) => Some(d),
            InterferenceData::EdgeVertex(d) => Some(&d.shape_shape),
            _ => None,
        }
    }

    /// Mutable `ShapeShapeData` of any shape-shape subclass.
    pub fn shape_shape_data_mut(&mut self) -> Option<&mut ShapeShapeData> {
        match &mut self.data {
            InterferenceData::ShapeShape(d) | InterferenceData::FaceEdge(d) => Some(d),
            InterferenceData::EdgeVertex(d) => Some(&mut d.shape_shape),
            _ => None,
        }
    }

    /// `down_cast<TopOpeBRepDS_SurfaceCurveInterference>`.
    pub fn as_surface_curve(&self) -> Option<&SurfaceCurveData> {
        match &self.data {
            InterferenceData::SurfaceCurve(d) => Some(d),
            _ => None,
        }
    }

    /// Mutable `down_cast<TopOpeBRepDS_SurfaceCurveInterference>`.
    pub fn as_surface_curve_mut(&mut self) -> Option<&mut SurfaceCurveData> {
        match &mut self.data {
            InterferenceData::SurfaceCurve(d) => Some(d),
            _ => None,
        }
    }

    /// `down_cast<TopOpeBRepDS_EdgeVertexInterference>`.
    pub fn as_edge_vertex(&self) -> Option<&EdgeVertexData> {
        match &self.data {
            InterferenceData::EdgeVertex(d) => Some(d),
            _ => None,
        }
    }

    /// The subclass payload.
    pub fn data(&self) -> &InterferenceData {
        &self.data
    }

    /// `ShapeShapeInterference::Parameter()` (`EdgeVertexInterference.cxx:56-59`)
    /// / `CurvePointInterference::Parameter()` (`CurvePointInterference.cxx:37-40`):
    /// the parameter of the interference on its support. `None` for the
    /// subclasses that carry no parameter.
    pub fn parameter(&self) -> Option<f64> {
        match self.data {
            InterferenceData::CurvePoint { param } => Some(param),
            InterferenceData::EdgeVertex(d) => Some(d.param),
            _ => None,
        }
    }

    /// `occ::down_cast<TopOpeBRepDS_CurvePointInterference>(I).IsNull() == false`
    /// (`CurvePointInterference.hxx:28-45`): true only for that subclass, not
    /// for `TopOpeBRepDS_EdgeVertexInterference`, which also carries a
    /// parameter.
    pub fn is_curve_point(&self) -> bool {
        matches!(self.data, InterferenceData::CurvePoint { .. })
    }

    /// `IsKind(STANDARD_TYPE(TopOpeBRepDS_EdgeVertexInterference))`.
    pub fn is_edge_vertex(&self) -> bool {
        matches!(self.data, InterferenceData::EdgeVertex(_))
    }

    /// `IsKind(STANDARD_TYPE(TopOpeBRepDS_FaceEdgeInterference))`.
    pub fn is_face_edge(&self) -> bool {
        matches!(self.data, InterferenceData::FaceEdge(_))
    }

    /// `Parameter(P)` of the parameter-carrying subclasses
    /// (`EdgeVertexInterference.cxx:62-65`, `CurvePointInterference.cxx:42-45`).
    pub fn set_parameter(&mut self, p: f64) {
        match &mut self.data {
            InterferenceData::CurvePoint { param } => *param = p,
            InterferenceData::EdgeVertex(d) => d.param = p,
            _ => {}
        }
    }

    /// `ShapeShapeInterference::Config()` (`ShapeShapeInterference.cxx:43-46`).
    pub fn config(&self) -> Option<Config> {
        self.shape_shape_data().map(|d| d.config)
    }

    /// `ShapeShapeInterference::GBound()` (`cxx:50-53`).
    pub fn g_bound(&self) -> Option<bool> {
        self.shape_shape_data().map(|d| d.g_bound)
    }

    /// `ShapeShapeInterference::SetGBound(b)` (`cxx:55-58`).
    pub fn set_g_bound(&mut self, b: bool) {
        if let Some(d) = self.shape_shape_data_mut() {
            d.g_bound = b;
        }
    }

    /// `SurfaceCurveInterference::PCurve()` (`cxx:50-53`).
    pub fn p_curve(&self) -> Option<&Arc<dyn Curve2d>> {
        self.as_surface_curve().and_then(|d| d.p_curve.as_ref())
    }

    /// `SurfaceCurveInterference::PCurve(PC)` (`cxx:57-60`).
    pub fn set_p_curve(&mut self, pc: Option<Arc<dyn Curve2d>>) {
        if let Some(d) = self.as_surface_curve_mut() {
            d.p_curve = pc;
        }
    }

    /// `Transition()` (`cxx:55-58`). OCCT returns a const reference; the value
    /// is returned by copy here.
    pub fn transition(&self) -> Transition {
        self.transition.clone()
    }

    /// `ChangeTransition()` (`cxx:62-65`).
    pub fn change_transition(&mut self) -> &mut Transition {
        &mut self.transition
    }

    /// `Transition(T)` (`cxx:70-73`).
    pub fn set_transition(&mut self, t: Transition) {
        self.transition = t;
    }

    /// `GKGSKS(GK, G, SK, S)` (`cxx:77-86`): geometry kind + index, support
    /// kind + index.
    pub fn gkgsks(&self) -> (Kind, i32, Kind, i32) {
        (
            self.geometry_type,
            self.geometry,
            self.support_type,
            self.support,
        )
    }

    /// `SupportType()` (`cxx:90-93`).
    pub fn support_type(&self) -> Kind {
        self.support_type
    }

    /// `Support()` (`cxx:95-98`).
    pub fn support(&self) -> i32 {
        self.support
    }

    /// `GeometryType()` (`cxx:102-105`).
    pub fn geometry_type(&self) -> Kind {
        self.geometry_type
    }

    /// `Geometry()` (`cxx:111-114`).
    pub fn geometry(&self) -> i32 {
        self.geometry
    }

    /// `SetGeometry(GI)` (`cxx:118-121`).
    pub fn set_geometry(&mut self, gi: i32) {
        self.geometry = gi;
    }

    /// `SupportType(ST)` (`cxx:125-128`).
    pub fn set_support_type(&mut self, st: Kind) {
        self.support_type = st;
    }

    /// `Support(S)` (`cxx:132-135`).
    pub fn set_support(&mut self, s: i32) {
        self.support = s;
    }

    /// `GeometryType(GT)` (`cxx:139-142`).
    pub fn set_geometry_type(&mut self, gt: Kind) {
        self.geometry_type = gt;
    }

    /// `Geometry(G)` (`cxx:146-149`). Same field as
    /// [`Interference::set_geometry`] (`SetGeometry`); OCCT exposes both
    /// spellings, the Rust API keeps one.
    pub fn set_geometry_index(&mut self, g: i32) {
        self.geometry = g;
    }

    /// `HasSameSupport(I)` (`cxx:153-157`).
    pub fn has_same_support(&self, other: &Interference) -> bool {
        self.support_type == other.support_type && self.support == other.support
    }

    /// `HasSameGeometry(I)` (`cxx:161-165`).
    pub fn has_same_geometry(&self, other: &Interference) -> bool {
        self.geometry_type == other.geometry_type && self.geometry == other.geometry
    }
}
