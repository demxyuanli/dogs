//! `Adaptor3d_TopolTool` stand-in used by GeomInt LineConstructor / IntSS.
//!
//! With a face, classification is `IntTools_FClass2d`. Without a face, the
//! surface UV box is the domain (`Adaptor3d_TopolTool(GeomAdaptor_Surface)`).

use std::sync::Arc;

use occt_core::gp::{GpDir2d, GpPnt2d, GpVec2d};
use occt_core::precision::{Precision, CONFUSION, INFINITE};
use occt_geom::Surface;
use occt_geom2d::{curve::Curve2d, Geom2dLine};

use crate::abs::Orientation;
use crate::brep_surface::{classify_surface, face_uv_bounds, SurfaceKind};
use crate::fclass2d::{FClass2d, FaceState};
use crate::shape::Face;

/// Resolution stored on the restriction vertexes built by
/// `Adaptor3d_TopolTool::Initialize(C)` (`Adaptor3d_TopolTool.cxx:245`,
/// `:251`).
const HVERTEX_RESOLUTION: f64 = 1.0e-8;

/// `Adaptor3d_HVertex` (`Adaptor3d_HVertex.hxx:24-51`) — a 2D point on a
/// restriction arc, its orientation and its parametric resolution.
#[derive(Debug, Clone, Copy, PartialEq)]
struct HVertex {
    /// `Value()` (`Adaptor3d_HVertex.cxx:38-41`).
    pnt: GpPnt2d,
    /// `Orientation()` (`Adaptor3d_HVertex.cxx:53-56`).
    orientation: Orientation,
    /// `Resolution(C)` (`Adaptor3d_HVertex.cxx:48-51`).
    resolution: f64,
    /// `Parameter(C)` (`Adaptor3d_HVertex.cxx:43-46`). The port's restriction
    /// arcs are straight 2D isolines, so `ElCLib::Parameter(Line, P)` is the
    /// arc parameter the vertex was built at.
    param: f64,
}

/// One entry of `Adaptor3d_TopolTool::myVtx` (`Adaptor3d_TopolTool.hxx`) plus
/// the restriction arc it was built on.
#[derive(Debug, Clone, Copy, PartialEq)]
struct ArcVertex {
    arc: RestrictionArc,
    hvertex: HVertex,
}

/// Domain classifier for one surface (`Adaptor3d_TopolTool`).
#[derive(Clone)]
pub struct TopolTool {
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
    class2d: Option<FClass2d>,
    /// `Adaptor3d_TopolTool::myVtx` — the vertexes of the restriction arcs.
    /// OCCT rebuilds them per arc in `Initialize(C)`
    /// (`Adaptor3d_TopolTool.cxx:235-254`); the port holds every arc's
    /// vertexes in one table and identifies them by index.
    vertices: Vec<ArcVertex>,
}

impl TopolTool {
    /// UV-box domain of `s` (no face wires).
    pub fn from_surface(s: &dyn Surface) -> Self {
        let (umin, umax) = s.u_range();
        let (vmin, vmax) = s.v_range();
        Self {
            umin,
            umax,
            vmin,
            vmax,
            class2d: None,
            vertices: Vec::new(),
        }
        .with_arc_vertices()
    }

    /// Face domain via `FClass2d`, falling back to UV bounds.
    pub fn from_face(face: &Face, s: &dyn Surface) -> Self {
        let (umin, umax, vmin, vmax) = face_uv_bounds(face);
        let class2d = FClass2d::new(face, crate::brep_tool::BRepTool::face_tolerance(face)).ok();
        let (umin, umax, vmin, vmax) = if umin.is_finite() {
            (umin, umax, vmin, vmax)
        } else {
            let (u0, u1) = s.u_range();
            let (v0, v1) = s.v_range();
            (u0, u1, v0, v1)
        };
        Self {
            umin,
            umax,
            vmin,
            vmax,
            class2d,
            vertices: Vec::new(),
        }
        .with_arc_vertices()
    }

    pub fn u_bounds(&self) -> (f64, f64) {
        (self.umin, self.umax)
    }

    pub fn v_bounds(&self) -> (f64, f64) {
        (self.vmin, self.vmax)
    }

    /// Fill `myVtx` from the restriction arcs, i.e. run
    /// `Adaptor3d_TopolTool::Initialize(C)` (`Adaptor3d_TopolTool.cxx:235-254`)
    /// once per arc of the UV box.
    fn with_arc_vertices(mut self) -> Self {
        let mut vertices = Vec::new();
        for arc in self.restriction_arcs() {
            // `theUinf > -myInfinite` (`Adaptor3d_TopolTool.cxx:243`).
            if arc.first > -INFINITE {
                vertices.push(ArcVertex {
                    arc,
                    hvertex: HVertex {
                        pnt: arc.value(arc.first),
                        orientation: Orientation::Forward,
                        resolution: HVERTEX_RESOLUTION,
                        param: arc.first,
                    },
                });
            }
            // `theUsup < myInfinite` (`Adaptor3d_TopolTool.cxx:249`).
            if arc.last < INFINITE {
                vertices.push(ArcVertex {
                    arc,
                    hvertex: HVertex {
                        pnt: arc.value(arc.last),
                        orientation: Orientation::Reversed,
                        resolution: HVERTEX_RESOLUTION,
                        param: arc.last,
                    },
                });
            }
        }
        self.vertices = vertices;
        self
    }

    /// Table indexes of the vertexes `Initialize(arc)` builds for `arc`
    /// (`Adaptor3d_TopolTool.cxx:235-254`), in table order.
    pub fn arc_vertexes(&self, arc: &RestrictionArc) -> Vec<usize> {
        self.vertices
            .iter()
            .enumerate()
            .filter(|(_, entry)| entry.arc == *arc)
            .map(|(index, _)| index)
            .collect()
    }

    /// `IntPatch_HInterTool::Parameter(V, C)` (`IntPatch_HInterTool.cxx:312-316`)
    /// = `Adaptor3d_HVertex::Parameter` (`Adaptor3d_HVertex.cxx:43-46`).
    pub fn vertex_parameter(&self, index: usize) -> f64 {
        self.vertices[index].hvertex.param
    }

    /// `IntPatch_HInterTool::Tolerance(V, C)` (`IntPatch_HInterTool.cxx:306-310`)
    /// = `Adaptor3d_HVertex::Resolution` (`Adaptor3d_HVertex.cxx:48-51`).
    pub fn vertex_resolution(&self, index: usize) -> f64 {
        self.vertices[index].hvertex.resolution
    }

    /// `Adaptor3d_TopolTool::Identical(V1, V2)`
    /// (`Adaptor3d_TopolTool.cxx:616-620`) = `Adaptor3d_HVertex::IsSame`
    /// (`Adaptor3d_HVertex.cxx:58-61`).
    pub fn identical(&self, first: usize, second: usize) -> bool {
        self.vertices[first]
            .hvertex
            .pnt
            .distance(&self.vertices[second].hvertex.pnt)
            <= CONFUSION
    }

    /// `Adaptor3d_TopolTool::Orientation(HVertex)`
    /// (`Adaptor3d_TopolTool.cxx:611-614`).
    pub fn vertex_orientation(&self, index: usize) -> Orientation {
        self.vertices[index].hvertex.orientation
    }

    /// `Adaptor3d_TopolTool::Classify(P, Tol)`.
    pub fn classify(&self, p: GpPnt2d, tol: f64) -> FaceState {
        if let Some(c) = &self.class2d {
            return c.perform(p);
        }
        classify_uv_box(p, self.umin, self.umax, self.vmin, self.vmax, tol)
    }

    /// `Adaptor3d_TopolTool::NbSamplesU`.
    pub fn nb_samples_u(&self, s: &dyn Surface) -> i32 {
        sample_counts(s).0
    }

    /// `Adaptor3d_TopolTool::NbSamplesV`.
    pub fn nb_samples_v(&self, s: &dyn Surface) -> i32 {
        sample_counts(s).1
    }

    /// `Adaptor3d_TopolTool::NbSamples`.
    pub fn nb_samples(&self, s: &dyn Surface) -> i32 {
        let (nu, nv) = sample_counts(s);
        nu * nv
    }

    /// `Adaptor3d_TopolTool::SamplePoint` (1-based, `myUPars` null branch).
    pub fn sample_point(&self, s: &dyn Surface, i: i32) -> (GpPnt2d, occt_core::gp::GpPnt) {
        let (nu, nv) = sample_counts(s);
        let (uinf, usup) = finite_uv(self.umin, self.umax);
        let (vinf, vsup) = finite_uv(self.vmin, self.vmax);
        let du = (usup - uinf) / (nu as f64 + 1.0);
        let dv = (vsup - vinf) / (nv as f64 + 1.0);
        let iv = 1 + i / nu;
        let iu = 1 + i - (iv - 1) * nu;
        let u = uinf + iu as f64 * du;
        let v = vinf + iv as f64 * dv;
        let p2 = GpPnt2d::new(u, v);
        (p2, s.d0(u, v))
    }

    /// Finite sides of the UV box, as 2D isolines (`Adaptor3d_TopolTool::Init`
    /// on a `GeomAdaptor_Surface` without a topological face). Infinite
    /// bounds contribute no arc.
    pub fn restriction_arcs(&self) -> Vec<RestrictionArc> {
        let mut arcs = Vec::new();
        let u_fin = !Precision::is_infinite(self.umin) && !Precision::is_infinite(self.umax);
        let v_fin = !Precision::is_infinite(self.vmin) && !Precision::is_infinite(self.vmax);
        if v_fin {
            if !Precision::is_infinite(self.umin) {
                arcs.push(RestrictionArc::v_iso(self.umin, self.vmin, self.vmax));
            }
            if !Precision::is_infinite(self.umax) && (self.umax - self.umin).abs() > 1e-14 {
                arcs.push(RestrictionArc::v_iso(self.umax, self.vmin, self.vmax));
            }
        }
        if u_fin {
            if !Precision::is_infinite(self.vmin) {
                arcs.push(RestrictionArc::u_iso(self.vmin, self.umin, self.umax));
            }
            if !Precision::is_infinite(self.vmax) && (self.vmax - self.vmin).abs() > 1e-14 {
                arcs.push(RestrictionArc::u_iso(self.vmax, self.umin, self.umax));
            }
        }
        arcs
    }
}

/// One UV-box side: a 2D line segment `t ∈ [first, last]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RestrictionArc {
    pub first: f64,
    pub last: f64,
    u0: f64,
    du: f64,
    v0: f64,
    dv: f64,
}

impl RestrictionArc {
    fn v_iso(u: f64, vmin: f64, vmax: f64) -> Self {
        Self {
            first: vmin,
            last: vmax,
            u0: u,
            du: 0.0,
            v0: 0.0,
            dv: 1.0,
        }
    }

    fn u_iso(v: f64, umin: f64, umax: f64) -> Self {
        Self {
            first: umin,
            last: umax,
            u0: 0.0,
            du: 1.0,
            v0: v,
            dv: 0.0,
        }
    }

    pub fn value(&self, t: f64) -> GpPnt2d {
        GpPnt2d::new(self.u0 + self.du * t, self.v0 + self.dv * t)
    }

    /// Isoline derivative (`Adaptor2d_Curve2d::D1`).
    pub fn d1(&self, t: f64) -> (GpPnt2d, GpVec2d) {
        (self.value(t), GpVec2d::new(self.du, self.dv))
    }

    /// Isoline as `Geom2d_Line` with the same parameter as `value`.
    pub fn to_curve2d(&self) -> Arc<dyn Curve2d> {
        let loc = self.value(0.0);
        let dir = if self.du.abs() <= self.dv.abs() {
            GpDir2d::new(0.0, 1.0).expect("v iso")
        } else {
            GpDir2d::new(1.0, 0.0).expect("u iso")
        };
        Arc::new(Geom2dLine::from_pnt_dir(loc, dir))
    }

    /// Orthogonal parameter of `uv` on this isoline (`HInterTool::Project`).
    pub fn project_uv(&self, uv: GpPnt2d) -> (f64, GpPnt2d) {
        let t = if self.du.abs() <= self.dv.abs() {
            uv.y()
        } else {
            uv.x()
        };
        (t, self.value(t))
    }
}

fn classify_uv_box(p: GpPnt2d, umin: f64, umax: f64, vmin: f64, vmax: f64, tol: f64) -> FaceState {
    let u = p.x();
    let v = p.y();
    let u_lo = if Precision::is_infinite(umin) {
        false
    } else {
        u < umin - tol
    };
    let u_hi = if Precision::is_infinite(umax) {
        false
    } else {
        u > umax + tol
    };
    let v_lo = if Precision::is_infinite(vmin) {
        false
    } else {
        v < vmin - tol
    };
    let v_hi = if Precision::is_infinite(vmax) {
        false
    } else {
        v > vmax + tol
    };
    if u_lo || u_hi || v_lo || v_hi {
        return FaceState::Out;
    }
    let on_u = (!Precision::is_infinite(umin) && (u - umin).abs() <= tol)
        || (!Precision::is_infinite(umax) && (u - umax).abs() <= tol);
    let on_v = (!Precision::is_infinite(vmin) && (v - vmin).abs() <= tol)
        || (!Precision::is_infinite(vmax) && (v - vmax).abs() <= tol);
    if on_u || on_v {
        FaceState::On
    } else {
        FaceState::In
    }
}

/// `Adaptor3d_TopolTool::ComputeSamplePoints` counts (minimum 6).
fn sample_counts(s: &dyn Surface) -> (i32, i32) {
    let (mut nbsu, mut nbsv) = match classify_surface(s) {
        SurfaceKind::Plane => (2, 2),
        SurfaceKind::Cylinder | SurfaceKind::Cone | SurfaceKind::Sphere | SurfaceKind::Torus => {
            (15, 15)
        }
        SurfaceKind::Other => (10, 10),
    };
    if nbsu < 6 {
        nbsu = 6;
    }
    if nbsv < 6 {
        nbsv = 6;
    }
    (nbsu, nbsv)
}

fn finite_uv(a: f64, b: f64) -> (f64, f64) {
    let mut lo = a;
    let mut hi = b;
    if hi < lo {
        std::mem::swap(&mut lo, &mut hi);
    }
    if Precision::is_infinite(lo) && Precision::is_infinite(hi) {
        return (-1.0e5, 1.0e5);
    }
    if Precision::is_infinite(lo) {
        lo = hi - 2.0e5;
    }
    if Precision::is_infinite(hi) {
        hi = lo + 2.0e5;
    }
    (lo, hi)
}
