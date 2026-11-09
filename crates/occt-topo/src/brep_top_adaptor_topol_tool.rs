//! `BRepTopAdaptor_TopolTool` — a face's restriction arcs and its 2D
//! classifier, seen through the `Adaptor3d_TopolTool` interface.
//!
//! Source: `ModelingAlgorithms/TKTopAlgo/BRepTopAdaptor/BRepTopAdaptor_TopolTool.hxx:36-137`,
//! `.cxx:43-653`.
//!
//! Translation notes:
//! - `Initialize(S)` (`.cxx:71-95`) requires an `Adaptor3d_Surface` that really
//!   is a `BRepAdaptor_Surface`, because it reads the face out of it
//!   (`.cxx:76-80`); the port's `Surface` trait carries no face, so the face is
//!   passed explicitly to [`TopolTool::initialize`]. A null face mirrors OCCT's
//!   `Standard_ConstructionError`.
//! - `myCurves` (`.cxx:87-92`) is `TopExp_Explorer(myFace, TopAbs_EDGE)` order;
//!   [`crate::topexp::Explorer`] is the port of that walk.
//! - `myFClass2d` (`.cxx:176-202`) is a lazily built `BRepTopAdaptor_FClass2d`
//!   whose tolerance is frozen on the first `Classify` / `IsThePointOn` call.
//!   The port's [`FClass2d`] is the port of `IntTools_FClass2d`, which samples
//!   the same face wires; the lazy, first-tolerance-wins storage is kept.
//! - `TestOnRestriction(P, Tol, RecadreOnPeriodic)` has no port counterpart
//!   taking `RecadreOnPeriodic`; [`FClass2d::test_on_restriction`] always folds
//!   periodic images, so an `IsThePointOn(..., false)` call also folds here.
//! - `ComputeSamplePoints` (`.cxx:345-524`) distinguishes Bezier / BSpline
//!   surfaces, surface of revolution / extrusion and the analytic kinds. The
//!   port's `SurfaceKind` (`brep_surface.rs`) only reports the analytic family,
//!   so the remaining families are detected through the `Surface` trait
//!   predicates (`is_bezier_surface`, `is_bspline_surface`,
//!   `is_surface_of_revolution`, `is_surface_of_linear_extrusion`).
//! - `Analyse` (`.cxx:232-342`) reads the poles of a Bezier / BSpline surface
//!   as a 1-based `(i = U, j = V)` array; `Surface::bezier_control_grid` is
//!   row-major `poles[i * n_v + j]` and `GeomBSplineSurface::pole(u, v)` is
//!   0-based, both mapped onto OCCT's 1-based `Value(i, j)`.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::precision::Precision;
use occt_geom::Surface;

use crate::abs::{Orientation, ShapeType};
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_surface::classify_surface;
use crate::brep_tool::BRepTool;
use crate::brep_top_adaptor_hvertex::HVertex;
use crate::fclass2d::{FClass2d, FaceState};
use crate::shape::{Edge, Face, Vertex};
use crate::topexp::Explorer;

/// `class BRepTopAdaptor_TopolTool : public Adaptor3d_TopolTool`
/// (`hxx:36-137`).
#[derive(Clone)]
pub struct TopolTool {
    /// `myS` (`Adaptor3d_TopolTool`): the surface the arcs belong to.
    surface: Option<Arc<dyn Surface>>,
    /// `myFace` (`.cxx:80`).
    face: Option<Face>,
    /// `myFClass2d` (`.cxx:180-187`): built on the first classification call,
    /// with that call's tolerance.
    fclass2d: Option<FClass2d>,
    /// `myCurves` (`.cxx:87-92`): one `BRepAdaptor_Curve2d` per face edge, in
    /// `TopExp_Explorer` order.
    curves: Vec<Curve2dAdaptor>,
    /// `myCIterator` position; `curves.len()` means `More()` is false.
    curve_index: usize,
    /// `myCurve` (`.cxx:99-106`): the arc handed to `InitVertexIterator`.
    current_curve: Option<Curve2dAdaptor>,
    /// `myVIterator`: the vertexes of `myCurve`'s edge.
    vertexes: Vec<Vertex>,
    /// `myVIterator` position; `vertexes.len()` means `MoreVertex()` is false.
    vertex_index: usize,
    /// `myNbSamplesU` / `myNbSamplesV` (`.cxx:360`, `-1` = not computed).
    nb_samples_u: i32,
    nb_samples_v: i32,
    /// `myU0` / `myV0` / `myDU` / `myDV` (`.cxx:517-521`).
    u0: f64,
    v0: f64,
    du: f64,
    dv: f64,
}

impl Default for TopolTool {
    /// `BRepTopAdaptor_TopolTool()` (`.cxx:43-52`): no face, no curves,
    /// `myNbSamplesU = -1`.
    fn default() -> Self {
        Self::new()
    }
}

impl TopolTool {
    /// `BRepTopAdaptor_TopolTool()` (`.cxx:43-52`).
    pub fn new() -> Self {
        Self {
            surface: None,
            face: None,
            fclass2d: None,
            curves: Vec::new(),
            curve_index: 0,
            current_curve: None,
            vertexes: Vec::new(),
            vertex_index: 0,
            nb_samples_u: -1,
            nb_samples_v: 0,
            u0: 0.0,
            v0: 0.0,
            du: 0.0,
            dv: 0.0,
        }
    }

    /// `Initialize(S)` (`.cxx:71-95`): remember the surface, take `S`'s face,
    /// drop the classifier and the sample counts, and rebuild `myCurves` from
    /// every edge of the face.
    ///
    /// `None` mirrors OCCT's `Standard_ConstructionError` for a surface that is
    /// not a `BRepAdaptor_Surface` (`cxx:76-80`) and for a null face.
    pub fn initialize(&mut self, surface: Arc<dyn Surface>, face: &Face) -> Result<(), String> {
        self.face = Some(face.clone());
        self.fclass2d = None;
        self.nb_samples_u = -1;
        self.surface = Some(surface);
        self.curves.clear();
        let mut ex = Explorer::new(face.shape(), ShapeType::Edge);
        while ex.more() {
            let e = Edge::wrap(ex.current().clone())
                .ok_or_else(|| "BRepTopAdaptor_TopolTool::Initialize: not an edge".to_string())?;
            self.curves.push(Curve2dAdaptor::of_edge_face(&e, face));
            ex.next();
        }
        self.curve_index = 0;
        self.vertexes.clear();
        self.vertex_index = 0;
        Ok(())
    }

    /// `Initialize(C)` (`.cxx:98-106`): the arc whose vertexes the iterator
    /// walks.
    pub fn initialize_curve(&mut self, c: &Curve2dAdaptor) {
        self.current_curve = Some(c.clone());
    }

    /// The `Initialize(S)` state for a surface that has no topological face.
    ///
    /// OCCT's `Initialize(S)` (`.cxx:71-95`) takes the face out of the
    /// `BRepAdaptor_Surface` and throws `Standard_ConstructionError` when
    /// there is none; a `GeomAdaptor_Surface` domain therefore never reaches
    /// `BRepTopAdaptor_TopolTool`. The port keeps one `geom_int::TopolTool`
    /// view that also covers the surface-only case
    /// (`geom_int::TopolTool::from_surface`), and this is its counterpart
    /// here: the surface is remembered (so `NbSamplesU/V` still work) and
    /// `myCurves` stays empty, i.e. the domain has no restriction arcs.
    pub fn initialize_surface_only(&mut self, surface: Arc<dyn Surface>) {
        self.face = None;
        self.fclass2d = None;
        self.nb_samples_u = -1;
        self.surface = Some(surface);
        self.curves.clear();
        self.curve_index = 0;
        self.vertexes.clear();
        self.vertex_index = 0;
    }

    /// `Initialize()` (`.cxx:64-67`): `Standard_NotImplemented`.
    pub fn initialize_not_implemented() -> Result<(), String> {
        Err("BRepTopAdaptor_TopolTool::Initialize() is not implemented".to_string())
    }

    /// `Init()` (`.cxx:109-112`).
    pub fn init(&mut self) {
        self.curve_index = 0;
    }

    /// `More()` (`.cxx:116-119`).
    pub fn more(&self) -> bool {
        self.curve_index < self.curves.len()
    }

    /// `Next()` (`.cxx:123-126`).
    pub fn next(&mut self) {
        self.curve_index += 1;
    }

    /// `Value()` (`.cxx:130-135`): the current arc.
    pub fn value(&self) -> Option<&Curve2dAdaptor> {
        self.curves.get(self.curve_index)
    }

    /// `Edge()` (`.cxx:138-143`): the current arc's edge.
    pub fn edge(&self) -> Option<&Edge> {
        self.curves.get(self.curve_index).map(|c| c.edge())
    }

    /// `InitVertexIterator()` (`.cxx:148-151`):
    /// `myVIterator.Init(myCurve->Edge(), TopAbs_VERTEX)`.
    ///
    /// A missing `myCurve` mirrors a null `TopoDS_Edge` handed to
    /// `TopExp_Explorer`, which yields nothing.
    pub fn init_vertex_iterator(&mut self) {
        self.vertexes.clear();
        self.vertex_index = 0;
        let Some(c) = self.current_curve.as_ref() else {
            return;
        };
        let mut ex = Explorer::new(c.edge().shape(), ShapeType::Vertex);
        while ex.more() {
            if let Some(v) = Vertex::wrap(ex.current().clone()) {
                self.vertexes.push(v);
            }
            ex.next();
        }
    }

    /// `MoreVertex()` (`.cxx:155-158`).
    pub fn more_vertex(&self) -> bool {
        self.vertex_index < self.vertexes.len()
    }

    /// `NextVertex()` (`.cxx:160-163`).
    pub fn next_vertex(&mut self) {
        self.vertex_index += 1;
    }

    /// `Vertex()` (`.cxx:167-170`): `BRepTopAdaptor_HVertex(current, myCurve)`.
    pub fn vertex(&self) -> Option<HVertex> {
        self.vertexes
            .get(self.vertex_index)
            .map(|v| HVertex::new(v))
    }

    /// `Classify(P, Tol, RecadreOnPeriodic = true)` (`.cxx:174-188`).
    ///
    /// `myFace.IsNull()` answers `TopAbs_UNKNOWN` (`.cxx:178-181`); the
    /// classifier is built once, with the first call's tolerance.
    pub fn classify(&mut self, p: GpPnt2d, tol: f64, recadre_on_periodic: bool) -> FaceState {
        let Some(face) = self.face.as_ref() else {
            return FaceState::Unknown;
        };
        if self.fclass2d.is_none() {
            match FClass2d::new(face, tol) {
                Ok(c) => self.fclass2d = Some(c),
                // OCCT's `BRepTopAdaptor_FClass2d` constructor has no failure
                // path; a port failure leaves `myFClass2d` null, which OCCT
                // would dereference. Report `UNKNOWN` instead.
                Err(_) => return FaceState::Unknown,
            }
        }
        match self.fclass2d.as_ref() {
            Some(c) => c.perform_recadre(p, recadre_on_periodic),
            None => FaceState::Unknown,
        }
    }

    /// `Identical(V1, V2)` (`.cxx:616-620`): `V1->IsSame(V2)`.
    pub fn identical(v1: &HVertex, v2: &HVertex) -> bool {
        v1.is_same(v2)
    }

    /// `IsThePointOn(P, Tol, RecadreOnPeriodic = true)` (`.cxx:191-202`).
    pub fn is_the_point_on(&mut self, p: GpPnt2d, tol: f64, _recadre_on_periodic: bool) -> bool {
        let Some(face) = self.face.as_ref() else {
            return false;
        };
        if self.fclass2d.is_none() {
            match FClass2d::new(face, tol) {
                Ok(c) => self.fclass2d = Some(c),
                Err(_) => return false,
            }
        }
        match self.fclass2d.as_ref() {
            Some(c) => c.test_on_restriction(p, tol) == FaceState::On,
            None => false,
        }
    }

    /// `Orientation(C)` (`.cxx:214-218`): the arc's edge orientation.
    pub fn orientation_of_curve(c: &Curve2dAdaptor) -> Orientation {
        c.edge().orientation()
    }

    /// `Has3d()` (`.cxx:598-601`).
    pub fn has3d(&self) -> bool {
        true
    }

    /// `Tol3d(C)` (`.cxx:605-620`): `BRep_Tool::Tolerance(edge)`.
    pub fn tol3d_of_curve(c: &Curve2dAdaptor) -> f64 {
        BRepTool::edge_tolerance(c.edge())
    }

    /// `Tol3d(V)` (`.cxx:623-637`): `BRep_Tool::Tolerance(vertex)`.
    pub fn tol3d_of_vertex(v: &HVertex) -> f64 {
        BRepTool::vertex_tolerance(v.vertex())
    }

    /// `Pnt(V)` (`.cxx:640-653`): `BRep_Tool::Pnt(vertex)`.
    pub fn pnt(v: &HVertex) -> GpPnt {
        BRepTool::vertex_point_world(v.vertex())
    }

    /// `NbSamplesU()` (`.cxx:525-533`).
    pub fn nb_samples_u(&mut self) -> i32 {
        if self.nb_samples_u < 0 {
            self.compute_sample_points();
        }
        self.nb_samples_u
    }

    /// `NbSamplesV()` (`.cxx:536-544`).
    pub fn nb_samples_v(&mut self) -> i32 {
        if self.nb_samples_u < 0 {
            self.compute_sample_points();
        }
        self.nb_samples_v
    }

    /// `NbSamples()` (`.cxx:547-555`).
    pub fn nb_samples(&mut self) -> i32 {
        if self.nb_samples_u < 0 {
            self.compute_sample_points();
        }
        self.nb_samples_u * self.nb_samples_v
    }

    /// `SamplePoint(i, P2d, P3d)` (`.cxx:558-566`).
    pub fn sample_point(&self, i: i32) -> Option<(GpPnt2d, GpPnt)> {
        let s = self.surface.as_ref()?;
        let iv = 1 + i / self.nb_samples_u;
        let iu = 1 + i - (iv - 1) * self.nb_samples_u;
        let u = self.u0 + iu as f64 * self.du;
        let v = self.v0 + iv as f64 * self.dv;
        Some((GpPnt2d::new(u, v), s.value(u, v)))
    }

    /// `DomainIsInfinite()` (`.cxx:570-596`).
    pub fn domain_is_infinite(&self) -> bool {
        let Some(s) = self.surface.as_ref() else {
            return false;
        };
        let (uinf, usup) = s.u_range();
        let (vinf, vsup) = s.v_range();
        if Precision::is_negative_infinite(uinf) {
            return true;
        }
        if Precision::is_positive_infinite(usup) {
            return true;
        }
        if Precision::is_negative_infinite(vinf) {
            return true;
        }
        if Precision::is_positive_infinite(vsup) {
            return true;
        }
        false
    }

    /// `ComputeSamplePoints()` (`.cxx:345-524`).
    fn compute_sample_points(&mut self) {
        let Some(s) = self.surface.clone() else {
            self.nb_samples_u = 0;
            self.nb_samples_v = 0;
            return;
        };
        let (mut uinf, mut usup) = s.u_range();
        let (mut vinf, mut vsup) = s.v_range();
        if usup < uinf {
            std::mem::swap(&mut uinf, &mut usup);
        }
        if vsup < vinf {
            std::mem::swap(&mut vinf, &mut vsup);
        }
        // `.cxx:366-392`: OCCT compares against the exact `RealFirst()` /
        // `RealLast()` sentinels its adaptors report for an infinite side; the
        // port's adaptors report `±Precision::Infinite()`, so the same test is
        // written with the named predicates.
        let is_first = Precision::is_negative_infinite;
        let is_last = Precision::is_positive_infinite;
        if is_first(uinf) && is_last(usup) {
            uinf = -1.0e5;
            usup = 1.0e5;
        } else if is_first(uinf) {
            uinf = usup - 2.0e5;
        } else if is_last(usup) {
            usup = uinf + 2.0e5;
        }
        if is_first(vinf) && is_last(vsup) {
            vinf = -1.0e5;
            vsup = 1.0e5;
        } else if is_first(vinf) {
            vinf = vsup - 2.0e5;
        } else if is_last(vsup) {
            vsup = vinf + 2.0e5;
        }

        let bspline = if s.is_bspline_surface() {
            s.osculating_bspline()
        } else {
            None
        };
        let bezier_grid = if s.is_bezier_surface() {
            s.bezier_control_grid()
        } else {
            None
        };

        let (mut nbsu, mut nbsv) = match classify_surface(s.as_ref()) {
            crate::brep_surface::SurfaceKind::Plane => (2, 2),
            _ if bspline.is_some() => {
                let b = bspline.as_ref().expect("checked");
                let mut nbsv = (b.knots_v().len() as i32) * (b.deg_v() as i32);
                if nbsv < 4 {
                    nbsv = 4;
                }
                let mut nbsu = (b.knots_u().len() as i32) * (b.deg_u() as i32);
                if nbsu < 4 {
                    nbsu = 4;
                }
                (nbsu, nbsv)
            }
            _ if bezier_grid.is_some() => {
                let (_, nu, nv) = bezier_grid.as_ref().expect("checked");
                (3 + *nu as i32, 3 + *nv as i32)
            }
            crate::brep_surface::SurfaceKind::Cylinder
            | crate::brep_surface::SurfaceKind::Cone
            | crate::brep_surface::SurfaceKind::Sphere
            | crate::brep_surface::SurfaceKind::Torus => {
                let mut nbsu = (8.0 * (usup - uinf)) as i32;
                let mut nbsv = (7.0 * (vsup - vinf)) as i32;
                if nbsu < 5 {
                    nbsu = 5;
                }
                if nbsv < 5 {
                    nbsv = 5;
                }
                if nbsu > 30 {
                    nbsu = 30;
                }
                if nbsv > 15 {
                    nbsv = 15;
                }
                (nbsu, nbsv)
            }
            _ if s.is_surface_of_revolution() || s.is_surface_of_linear_extrusion() => (25, 15),
            _ => (10, 10),
        };

        // `.cxx:463-473`.
        if nbsu < 10 {
            nbsu = 10;
        }
        if nbsv < 10 {
            nbsv = 10;
        }
        self.nb_samples_u = nbsu;
        self.nb_samples_v = nbsv;

        // `.cxx:478-503`: only a Bezier / BSpline whose counts exceed 10 is
        // re-analysed from its pole grid.
        if (nbsu > 10 || nbsv > 10) && bspline.is_some() {
            let b = bspline.as_ref().expect("checked");
            let nbup = b.nb_poles_u();
            let nbvp = b.nb_poles_v();
            let (u, v) = analyse_poles(nbup, nbvp, |i, j| b.pole(i - 1, j - 1));
            self.nb_samples_u = u;
            self.nb_samples_v = v;
            nbsu = u;
            nbsv = v;
        } else if (nbsu > 10 || nbsv > 10) && bezier_grid.is_some() {
            let (poles, nu, nv) = bezier_grid.as_ref().expect("checked");
            let (u, v) = analyse_poles(*nu, *nv, |i, j| poles[(i - 1) * nv + (j - 1)]);
            self.nb_samples_u = u;
            self.nb_samples_v = v;
            nbsu = u;
            nbsv = v;
        }

        // `.cxx:505-514`.
        if nbsu < 10 {
            nbsu = 10;
        }
        if nbsv < 10 {
            nbsv = 10;
        }
        self.nb_samples_u = nbsu;
        self.nb_samples_v = nbsv;

        self.u0 = uinf;
        self.v0 = vinf;
        self.du = (usup - uinf) / (self.nb_samples_u as f64 + 1.0);
        self.dv = (vsup - vinf) / (self.nb_samples_v as f64 + 1.0);
    }
}

/// `Analyse` (`.cxx:232-342`): count the pole-grid direction changes of a
/// Bezier / BSpline surface and size the sample grid from them.
///
/// `pole(i, j)` is OCCT's 1-based `array2.Value(i, j)` with `i` the U index and
/// `j` the V index. `sh` is deliberately *not* reset between the two loops
/// (`.cxx:276-284` and `.cxx:316-330` share it), matching OCCT.
fn analyse_poles<F>(nbup: usize, nbvp: usize, pole: F) -> (i32, i32)
where
    F: Fn(usize, usize) -> GpPnt,
{
    let nbup = nbup as i32;
    let nbvp = nbvp as i32;
    let mut sh = 1;
    let mut nbch = 0;
    if nbvp > 2 {
        for i in 2..nbup {
            let a = pole(i as usize, 1);
            let b = pole(i as usize, 2);
            let c = pole(i as usize, 3);
            let mut vi = [
                c.x() - b.x() - b.x() + a.x(),
                c.y() - b.y() - b.y() + a.y(),
                c.z() - b.z() - b.z() + a.z(),
            ];
            let mut locnbch = 0;
            for j in 3..nbvp {
                let ax = pole(i as usize, j as usize - 1);
                let bx = pole(i as usize, j as usize);
                let cx = pole(i as usize, j as usize + 1);
                let vip1 = [
                    cx.x() - bx.x() - bx.x() + ax.x(),
                    cx.y() - bx.y() - bx.y() + ax.y(),
                    cx.z() - bx.z() - bx.z() + ax.z(),
                ];
                let pd = vi[0] * vip1[0] + vi[1] * vip1[1] + vi[2] * vip1[2];
                vi = vip1;
                if pd > 1.0e-7 || pd < -1.0e-7 {
                    if pd > 0.0 {
                        if sh == -1 {
                            sh = 1;
                            locnbch += 1;
                        }
                    } else if sh == 1 {
                        sh = -1;
                        locnbch += 1;
                    }
                }
            }
            if locnbch > nbch {
                nbch = locnbch;
            }
        }
    }
    let samples_v = nbch + 5;

    nbch = 0;
    if nbup > 2 {
        for j in 2..nbvp {
            let a = pole(1, j as usize);
            let b = pole(2, j as usize);
            let c = pole(3, j as usize);
            let mut vi = [
                c.x() - b.x() - b.x() + a.x(),
                c.y() - b.y() - b.y() + a.y(),
                c.z() - b.z() - b.z() + a.z(),
            ];
            let mut locnbch = 0;
            for i in 3..nbup {
                let ax = pole(i as usize - 1, j as usize);
                let bx = pole(i as usize, j as usize);
                let cx = pole(i as usize + 1, j as usize);
                let vip1 = [
                    cx.x() - bx.x() - bx.x() + ax.x(),
                    cx.y() - bx.y() - bx.y() + ax.y(),
                    cx.z() - bx.z() - bx.z() + ax.z(),
                ];
                let pd = vi[0] * vip1[0] + vi[1] * vip1[1] + vi[2] * vip1[2];
                vi = vip1;
                if pd > 1.0e-7 || pd < -1.0e-7 {
                    if pd > 0.0 {
                        if sh == -1 {
                            sh = 1;
                            locnbch += 1;
                        }
                    } else if sh == 1 {
                        sh = -1;
                        locnbch += 1;
                    }
                }
            }
            if locnbch > nbch {
                nbch = locnbch;
            }
        }
    }
    (nbch + 5, samples_v)
}
