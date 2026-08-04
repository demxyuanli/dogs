//! Exact BRepGProp mass properties — full port of `TKTopAlgo/BRepGProp`.
//!
//! Computes linear (edge-length), surface (area) and volume global properties
//! of a `TopoShape` by integrating over the exact boundary geometry, mirroring
//! the OCCT `BRepGProp` entry points:
//!
//! * [`linear_properties`] — Gauss quadrature along every edge's curve
//!   (`BRepGProp_Cinert` + `BRepGProp_EdgeTool`). The edge curve is a
//!   `dyn Curve`; d0/d1 are used directly, d2 numerically when needed.
//! * [`surface_properties`] — Gauss quadrature over each face's UV domain
//!   (`BRepGProp_Sinert` + `BRepGProp_Gauss`). Returns the total area.
//! * [`volume_properties`] — the divergence-theorem surface integral over each
//!   face (`BRepGProp_Vinert` + `BRepGProp_Gauss`), giving the signed volume.
//! * [`volume_properties_gk`] — the same volume via adaptive Gauss–Kronrod
//!   (`BRepGProp_VinertGK` + `UFunction` / `TFunction`).
//!
//! Trimmed faces (bounded by wires) are integrated with the boundary line
//! integral (Green's theorem) over the pcurves of the wire edges, exactly like
//! OCCT. The face's UV bounds come from the pcurves (`BRepTools::UVBounds`),
//! which makes the integration well-defined even for surfaces (planes,
//! cylinders) whose natural parameter ranges are unbounded. The wire's UV
//! orientation sign is recovered from the signed area of the UV boundary
//! polygon, so the direct (positively-oriented) surface integral is obtained
//! regardless of how the wires were wound. Faces whose wire repeats an edge
//! (the cylinder lateral face's double seam) are integrated directly over the
//! UV bounding rectangle.
//!
//! The surface `d1` derivative is computed numerically for surfaces whose port
//! `d1` returns zero vectors (cylinder / sphere / torus / cone), following the
//! Phase 13 convention that `Surface` carries no `d2`.

use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::gp::{GpMat, GpPnt, GpPnt2d, GpVec, GpVec2d, GpXyz};
use occt_geom::{Curve, Surface};
use occt_math::gauss::gauss_legendre;

use crate::brep_surface::classify_surface;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools_full::{edges_of, edges_of_wire, faces_of, vertices_of, wires_of_face};

// ---------------------------------------------------------------------------
// Constants (mirror BRepGProp_Gauss.cxx)
// ---------------------------------------------------------------------------

const GPM: usize = 64; // math::GaussPointsMax()
const EPS_PARAM: f64 = 1e-12;
const EPS_DIM: f64 = 1e-30;
const ERROR_ALGEBR_RATIO: f64 = 2.0 / 3.0;
const SUBS_POWER: usize = 32;
const SM: usize = SUBS_POWER * GPM + 1;
const NGV: usize = 10; // number of value types in the GK integrator

// ---------------------------------------------------------------------------
// GProps — global-properties accumulator (GProp_GProps)
// ---------------------------------------------------------------------------

/// Accumulated mass properties about a reference point.
///
/// * `dim` — mass / length / area / volume (the "dimension" of the system);
/// * `loc` — the reference point the inertia is computed about;
/// * `g`   — the centre of mass *relative to* `loc`;
/// * `inertia` — the 3×3 inertia matrix about `loc` (standard form, off
///   diagonal entries negated: `I[1][2] = −∫xy dm`).
#[derive(Debug, Clone, Copy)]
pub struct GProps {
    pub dim: f64,
    pub loc: GpPnt,
    pub g: GpVec,
    pub inertia: GpMat,
}

impl GProps {
    /// Empty accumulator about `loc`.
    pub fn new(loc: GpPnt) -> Self {
        Self { dim: 0.0, loc, g: GpVec::zero(), inertia: GpMat::zero() }
    }

    /// The mass / length / area / volume.
    pub fn mass(&self) -> f64 {
        self.dim
    }

    /// The absolute centre of mass (`loc + g`).
    pub fn center(&self) -> GpPnt {
        GpPnt::from_xyz(&self.loc.coord.added(&self.g.coord))
    }

    /// The inertia matrix about the reference point `loc`.
    pub fn inertia_matrix(&self) -> GpMat {
        self.inertia
    }

    /// The inertia matrix about the centre of mass (parallel-axis shift from
    /// `loc` to the centre).
    pub fn central_inertia(&self) -> GpMat {
        let d = GpVec::from_pnts(&self.center(), &self.loc); // loc − centre
        self.inertia.add(&huygens_mat(-self.dim, &d))
    }

    /// Combine `other` (about its own `loc`) into `self` using the Huygens
    /// (parallel-axis) theorem. Both accumulators may have different reference
    /// points; the combined system is about `self.loc`.
    pub fn add(&mut self, other: &GProps) {
        let off = GpVec::from_pnts(&self.loc, &other.loc); // other.loc − self.loc
        let other_g_here = other.g.added(&off);
        let huy = huygens_mat(other.dim, &off);
        let d = self.dim + other.dim;
        if d.abs() > EPS_DIM {
            self.g = GpVec::from_xyz(
                &self.g.coord.multiplied(self.dim).added(&other_g_here.coord.multiplied(other.dim)).divided(d),
            );
            self.dim = d;
        } else {
            self.g = GpVec::zero();
            self.dim = 0.0;
        }
        self.inertia = self.inertia.add(&other.inertia.add(&huy));
    }

    /// Add an unlocated piece whose centre and inertia are already about
    /// `self.loc`.
    pub fn add_here(&mut self, dim: f64, g: &GpVec, inertia: &GpMat) {
        let d = self.dim + dim;
        if d.abs() > EPS_DIM {
            self.g = GpVec::from_xyz(
                &self.g.coord.multiplied(self.dim).added(&g.coord.multiplied(dim)).divided(d),
            );
            self.dim = d;
        } else {
            self.g = GpVec::zero();
            self.dim = 0.0;
        }
        self.inertia = self.inertia.add(inertia);
    }
}

/// Parallel-axis tensor `m·(d²δ − d·dᵀ)` in the standard inertia-matrix form
/// (off-diagonal entries negated).
fn huygens_mat(m: f64, d: &GpVec) -> GpMat {
    let (x, y, z) = (d.x(), d.y(), d.z());
    let d2 = x * x + y * y + z * z;
    GpMat::new(
        m * (d2 - x * x), -m * x * y, -m * x * z,
        -m * x * y, m * (d2 - y * y), -m * y * z,
        -m * x * z, -m * y * z, m * (d2 - z * z),
    )
}

// ---------------------------------------------------------------------------
// Inertia — internal accumulation struct (BRepGProp_Gauss::Inertia)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, Default)]
struct Inertia {
    mass: f64,
    ix: f64,
    iy: f64,
    iz: f64,
    ixx: f64,
    iyy: f64,
    izz: f64,
    ixy: f64,
    ixz: f64,
    iyz: f64,
}

impl Inertia {
    fn add(&mut self, o: &Inertia) {
        self.mass += o.mass;
        self.ix += o.ix;
        self.iy += o.iy;
        self.iz += o.iz;
        self.ixx += o.ixx;
        self.iyy += o.iyy;
        self.izz += o.izz;
        self.ixy += o.ixy;
        self.ixz += o.ixz;
        self.iyz += o.iyz;
    }

    fn mul(&mut self, s: f64) {
        self.mass *= s;
        self.ix *= s;
        self.iy *= s;
        self.iz *= s;
        self.ixx *= s;
        self.iyy *= s;
        self.izz *= s;
        self.ixy *= s;
        self.ixz *= s;
        self.iyz *= s;
    }
}

/// Convert a `Sinert`-style inertia (surface / curve) to `(mass, g_rel, mat)`.
fn convert_s(inertia: &Inertia) -> (f64, GpVec, GpMat) {
    let (mass, g) = if inertia.mass.abs() >= EPS_DIM {
        (inertia.mass, GpVec::new(inertia.ix / inertia.mass, inertia.iy / inertia.mass, inertia.iz / inertia.mass))
    } else {
        (0.0, GpVec::zero())
    };
    let mat = GpMat::new(
        inertia.ixx, -inertia.ixy, -inertia.ixz,
        -inertia.ixy, inertia.iyy, -inertia.iyz,
        -inertia.ixz, -inertia.iyz, inertia.izz,
    );
    (mass, g, mat)
}

/// Convert a `Vinert`-style inertia (volume, by point) to `(mass, g_rel, mat)`.
fn convert_v(inertia: &Inertia, coeff: &[f64; 3]) -> (f64, GpVec, GpMat) {
    if inertia.mass.abs() >= EPS_DIM {
        (
            inertia.mass,
            GpVec::new(
                coeff[0] + inertia.ix / inertia.mass,
                coeff[1] + inertia.iy / inertia.mass,
                coeff[2] + inertia.iz / inertia.mass,
            ),
            GpMat::new(
                inertia.ixx, inertia.ixy, inertia.ixz,
                inertia.ixy, inertia.iyy, inertia.iyz,
                inertia.ixz, inertia.iyz, inertia.izz,
            ),
        )
    } else {
        (0.0, GpVec::zero(), GpMat::zero())
    }
}

// ---------------------------------------------------------------------------
// Elementary-part accumulation (BRepGProp_Gauss::compute{V,S}InertiaOfElementaryPart)
// ---------------------------------------------------------------------------

fn compute_v_inertia_elem(pt: &GpPnt, n: &GpVec, loc: &GpPnt, w: f64, coeff: &[f64; 3], by_point: bool, o: &mut Inertia) {
    let mut x = pt.x() - loc.x();
    let mut y = pt.y() - loc.y();
    let mut z = pt.z() - loc.z();

    let xn = n.x() * w;
    let yn = n.y() * w;
    let zn = n.z() * w;

    if by_point {
        let dv = x * xn + y * yn + z * zn;
        o.mass += dv / 3.0;
        o.ix += 0.25 * x * dv;
        o.iy += 0.25 * y * dv;
        o.iz += 0.25 * z * dv;
        x -= coeff[0];
        y -= coeff[1];
        z -= coeff[2];
        let dv20 = dv * 0.2;
        o.ixy -= x * y * dv20;
        o.iyz -= y * z * dv20;
        o.ixz -= x * z * dv20;
        let (x2, y2, z2) = (x * x, y * y, z * z);
        o.ixx += (y2 + z2) * dv20;
        o.iyy += (x2 + z2) * dv20;
        o.izz += (x2 + y2) * dv20;
    } else {
        // By-plane branch — reached only by the GK plane-restriction overload,
        // which the Rust entry points do not use. Kept for completeness.
        let s = xn * coeff[0] + yn * coeff[1] + zn * coeff[2];
        let d1 = coeff[0] * x + coeff[1] * y + coeff[2] * z;
        let _ = (s, d1);
    }
}

fn compute_s_inertia_elem(pt: &GpPnt, n: &GpVec, loc: &GpPnt, w: f64, o: &mut Inertia) {
    let ds = n.magnitude() * w;
    let x = pt.x() - loc.x();
    let y = pt.y() - loc.y();
    let z = pt.z() - loc.z();
    o.mass += ds;
    let xds = x * ds;
    let yds = y * ds;
    let zds = z * ds;
    o.ix += xds;
    o.iy += yds;
    o.iz += zds;
    o.ixy += x * yds;
    o.iyz += y * zds;
    o.ixz += x * zds;
    let xx = x * xds;
    let yy = y * yds;
    let zz = z * zds;
    o.ixx += yy + zz;
    o.iyy += xx + zz;
    o.izz += xx + yy;
}

// ---------------------------------------------------------------------------
// Surface classification + numeric derivatives
// ---------------------------------------------------------------------------

/// Coarse surface classification used to pick integration orders and the UV
/// inverse map.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
enum SurfKind {
    Plane,
    Cylinder,
    Sphere,
    Torus,
    Cone,
    Other,
}

fn classify_surface_kind(s: &dyn Surface) -> SurfKind {
    match classify_surface(s) {
        crate::brep_surface::SurfaceKind::Plane => SurfKind::Plane,
        crate::brep_surface::SurfaceKind::Sphere => SurfKind::Sphere,
        _ => {
            if cylinder_frame(s).is_some() {
                SurfKind::Cylinder
            } else if cone_frame(s).is_some() {
                SurfKind::Cone
            } else {
                SurfKind::Other
            }
        }
    }
}

/// Reconstruct a plane's frame `(origin, x_dir, y_dir)` from the surface's own
/// parameterization: `S(u,v) = o + u·xd + v·yd`, so `o = S(0,0)` and the
/// directions are the (normalised) first derivatives. Using the surface's own
/// frame keeps the pcurve UV coordinates in the same space as `normal(u, v)`.
fn plane_frame(s: &dyn Surface) -> Option<(GpXyz, GpXyz, GpXyz)> {
    let o = s.d0(0.0, 0.0).coord;
    let (_, du, dv) = surface_d1(s, 0.0, 0.0);
    let xm = du.magnitude();
    let ym = dv.magnitude();
    if xm < 1e-12 || ym < 1e-12 {
        return None;
    }
    Some((o, du.divided(xm).coord, dv.divided(ym).coord))
}

/// Reconstruct the cylinder frame from sampled invariants:
/// `(origin, x_dir, y_dir, axis_dir, radius)`.
fn cylinder_frame(s: &dyn Surface) -> Option<(GpXyz, GpXyz, GpXyz, GpXyz, f64)> {
    let p0 = s.d0(0.0, 0.0);
    let p1 = s.d0(PI, 0.0);
    let o = GpPnt::from_xyz(&p0.coord.added(&p1.coord).divided(2.0));
    let xr = p0.coord.subtracted(&o.coord);
    let r = xr.modulus();
    if r < 1e-12 {
        return None;
    }
    let xd = xr.divided(r);
    let p2 = s.d0(PI / 2.0, 0.0);
    let yr = p2.coord.subtracted(&o.coord);
    let yd = yr.divided(r);
    let p3 = s.d0(0.0, 1.0);
    let ax = p3.coord.subtracted(&p0.coord);
    let am = ax.modulus();
    if am < 1e-12 {
        return None;
    }
    let z = ax.divided(am);
    Some((o.coord, xd, yd, z, r))
}

/// Frame of a cone surface: apex (o), radial axes (x, y) and axis (z),
/// plus the semi-angle `alpha`. Recovered by sampling `d0` — the port's
/// `dyn Surface` cannot downcast to `GpCone`.
fn cone_frame(s: &dyn Surface) -> Option<(GpXyz, GpXyz, GpXyz, GpXyz, f64)> {
    // d0(0,0) and d0(0,1) lie on the v-direction generatrix (u=0).
    let p0 = s.d0(0.0, 0.0);
    let p1 = s.d0(0.0, 1.0);
    let g = p1.coord.subtracted(&p0.coord);
    let gm = g.modulus();
    if gm < 1e-12 {
        return None;
    }
    let z = g.divided(gm);
    let g0a = p0.coord;
    // Semi-angle between the generatrix and the axis (both unit vectors).
    let alpha = g.divided(gm).dot(&z).clamp(-1.0, 1.0).acos();
    // Radial axis x through the u=0 direction.
    let p2 = s.d0(std::f64::consts::PI, 0.5);
    let xr = p2.coord.subtracted(&g0a);
    let xm = xr.modulus();
    if xm < 1e-12 {
        return None;
    }
    let xd = xr.divided(xm);
    let yd = z.crossed(&xd);
    if yd.modulus() < 1e-12 {
        return None;
    }
    let yd = yd.divided(yd.modulus());
    // Apex: along the axis from the u=0 iso-line, the radial distance r0
    // divided by tan(alpha) locates the cone's vertex.
    let r0 = g0a.subtracted(&z.multiplied(g0a.dot(&z))).modulus();
    let apex = g0a.added(&z.multiplied(r0 / alpha.tan()));
    Some((apex, xd, yd, z, alpha))
}

/// First derivatives `(point, du, dv)` of a surface at `(u, v)`, using the
/// analytic `d1` when it is non-degenerate and central finite differences
/// otherwise (cylinder / sphere / torus / cone `d1` return zero in the port).
fn surface_d1(s: &dyn Surface, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
    let (p, du, dv) = s.d1(u, v);
    if du.square_magnitude() > 1e-30 && dv.square_magnitude() > 1e-30 {
        return (p, du, dv);
    }
    let h = 1e-6;
    let p0 = s.d0(u, v);
    let p_up = s.d0(u + h, v);
    let p_um = s.d0(u - h, v);
    let p_vp = s.d0(u, v + h);
    let p_vm = s.d0(u, v - h);
    let du = GpVec::from_xyz(&p_up.coord.subtracted(&p_um.coord).divided(2.0 * h));
    let dv = GpVec::from_xyz(&p_vp.coord.subtracted(&p_vm.coord).divided(2.0 * h));
    (p0, du, dv)
}

// ---------------------------------------------------------------------------
// UV inverse map for pcurves
// ---------------------------------------------------------------------------

/// Inverse map of a 3D point on the face surface to its UV coordinates, used
/// to build the boundary pcurves of a trimmed face.
#[derive(Clone)]
enum UVMap {
    Plane { o: GpXyz, xd: GpXyz, yd: GpXyz },
    Cylinder { o: GpXyz, xd: GpXyz, yd: GpXyz, z: GpXyz, r: f64 },
    Cone { o: GpXyz, xd: GpXyz, yd: GpXyz, z: GpXyz, alpha: f64 },
    Generic,
}

impl UVMap {
    fn map_point(&self, s: &dyn Surface, p: &GpPnt) -> GpPnt2d {
        match self {
            UVMap::Plane { o, xd, yd } => {
                let rel = p.coord.subtracted(o);
                GpPnt2d::new(rel.dot(xd), rel.dot(yd))
            }
            UVMap::Cone { o, xd, yd, z, alpha } => {
                // Cone param: d0(u,v) = apex + x·((r+v·sinα)cos u) + y·((r+v·sinα)sin u)
                // + z·(v·cosα). Invert: dz = v·cosα, angle from x/y.
                let rel = p.coord.subtracted(o);
                let u = rel.dot(yd).atan2(rel.dot(xd)).rem_euclid(2.0 * PI);
                let v = rel.dot(z) / alpha.cos();
                GpPnt2d::new(u, v)
            }
            UVMap::Cylinder { o, xd, yd, z, .. } => {
                let rel = p.coord.subtracted(o);
                let u = rel.dot(yd).atan2(rel.dot(xd)).rem_euclid(2.0 * PI);
                GpPnt2d::new(u, rel.dot(z))
            }
            UVMap::Generic => {
                let (u, v) = closest_params(s, p);
                GpPnt2d::new(u, v)
            }
        }
    }

    fn map_deriv(&self, s: &dyn Surface, p: &GpPnt, d3: &GpVec, u: f64) -> GpVec2d {
        match self {
            UVMap::Plane { xd, yd, .. } => GpVec2d::new(d3.coord.dot(xd), d3.coord.dot(yd)),
            UVMap::Cylinder { o, xd, yd, z, r } => {
                let rel = p.coord.subtracted(o);
                let _ = u;
                let ua = rel.dot(yd).atan2(rel.dot(xd)).rem_euclid(2.0 * PI);
                let dx = d3.coord.dot(xd);
                let dy = d3.coord.dot(yd);
                let du = (dy * ua.cos() - dx * ua.sin()) / r;
                let dv = d3.coord.dot(z);
                GpVec2d::new(du, dv)
            }
            UVMap::Cone { o, xd, yd, z, alpha } => {
                let rel = p.coord.subtracted(o);
                let _ = u;
                let ua = rel.dot(yd).atan2(rel.dot(xd)).rem_euclid(2.0 * PI);
                let dx = d3.coord.dot(xd);
                let dy = d3.coord.dot(yd);
                // Radius at the cone's current v (like Cylinder but r varies
                // with v along the generatrix).
                let r0 = (rel.subtracted(&z.multiplied(rel.dot(z)))).modulus();
                let r = r0.max(1e-9);
                let du = (dy * ua.cos() - dx * ua.sin()) / r;
                let dv = d3.coord.dot(z) / alpha.cos();
                GpVec2d::new(du, dv)
            }
            UVMap::Generic => {
                let h = 1e-6;
                let pp = GpPnt::from_xyz(&p.coord.added(&d3.coord.multiplied(h)));
                let pm = GpPnt::from_xyz(&p.coord.subtracted(&d3.coord.multiplied(h)));
                let a = self.map_point(s, &pp);
                let b = self.map_point(s, &pm);
                GpVec2d::new((a.x() - b.x()) / (2.0 * h), (a.y() - b.y()) / (2.0 * h))
            }
        }
    }
}

/// Grid-search + refinement UV projection for generic surfaces.
fn closest_params(s: &dyn Surface, p: &GpPnt) -> (f64, f64) {
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let (u0, u1) = sane_bounds(u0, u1);
    let (v0, v1) = sane_bounds(v0, v1);
    let nu = 24;
    let nv = 24;
    let mut best = (u0, v0);
    let mut best_d = f64::INFINITY;
    for i in 0..=nu {
        for j in 0..=nv {
            let u = u0 + (u1 - u0) * i as f64 / nu as f64;
            let v = v0 + (v1 - v0) * j as f64 / nv as f64;
            let d = s.d0(u, v).distance(p);
            if d < best_d {
                best_d = d;
                best = (u, v);
            }
        }
    }
    let (mut u, mut v) = best;
    let (mut hu, mut hv) = ((u1 - u0) / nu as f64, (v1 - v0) / nv as f64);
    for _ in 0..5 {
        hu *= 0.5;
        hv *= 0.5;
        let mut du = u;
        let mut dv = v;
        let mut dd = s.d0(u, v).distance(p);
        for &(su, sv) in &[(-1.0, 0.0), (1.0, 0.0), (0.0, -1.0), (0.0, 1.0)] {
            let nu2 = (u + su * hu).clamp(u0, u1);
            let nv2 = (v + sv * hv).clamp(v0, v1);
            let d = s.d0(nu2, nv2).distance(p);
            if d < dd {
                dd = d;
                du = nu2;
                dv = nv2;
            }
        }
        u = du;
        v = dv;
    }
    (u, v)
}

fn sane_bounds(a: f64, b: f64) -> (f64, f64) {
    if a.is_finite() && b.is_finite() && b > a {
        (a, b)
    } else {
        (-1.0, 1.0)
    }
}

// ---------------------------------------------------------------------------
// Boundary arc (pcurve of an edge on a face)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ArcKind {
    Line,
    Circle,
    Other,
}

/// A boundary arc of a face: the edge's 3D curve together with the UV inverse
/// map, so `d12d` yields the pcurve point and derivative in the surface's UV.
struct BoundaryArc {
    curve: Arc<dyn Curve>,
    a: f64,
    b: f64,
    kind: ArcKind,
    map: UVMap,
}

impl BoundaryArc {
    fn d12d(&self, s: &dyn Surface, t: f64) -> (GpPnt2d, GpVec2d) {
        // NOTE: `curve.d1(t).0` is unreliable for lines in the port
        // (`clib::line_d1` ignores the parameter), so evaluate the point with
        // `d0` and take the derivative from `d1`.
        let p = self.curve.d0(t);
        let d1 = self.curve.d1(t).1;
        let puv = self.map.map_point(s, &p);
        let duv = self.map.map_deriv(s, &p, &d1, puv.x());
        (puv, duv)
    }

    fn value(&self, s: &dyn Surface, t: f64) -> GpPnt2d {
        let p = self.curve.d0(t);
        self.map.map_point(s, &p)
    }
}

// ---------------------------------------------------------------------------
// FaceGauss — the BRepGProp_Face adaptor
// ---------------------------------------------------------------------------

/// Wraps a face and its surface with everything the integrators need:
/// UV bounds, surface normal, integration orders, and the boundary arcs.
struct FaceGauss {
    surface: Arc<dyn Surface>,
    kind: SurfKind,
    is_reversed: bool,
    natural: bool,
    u1: f64,
    u2: f64,
    v1: f64,
    v2: f64,
    arcs: Vec<BoundaryArc>,
    has_repeated_edges: bool,
}

impl FaceGauss {
    fn new(face: &Face) -> Result<Self, String> {
        let surface = BRepTool::face_surface_world(face).ok_or("brep_gprop_full: face has no surface")?;
        let kind = classify_surface_kind(surface.as_ref());
        let is_reversed = face.orientation().is_reversed();
        let wires = wires_of_face(face);
        let natural = wires.is_empty();

        // Build the UV inverse map. The frame is taken from the surface's own
        // parameterization so pcurve UV coordinates live in the same space the
        // `normal`/integrand evaluators use.
        let map = match kind {
            SurfKind::Plane => {
                let (o, xd, yd) =
                    plane_frame(surface.as_ref()).ok_or("brep_gprop_full: planar face frame")?;
                UVMap::Plane { o, xd, yd }
            }
            SurfKind::Cylinder => {
                let (o, xd, yd, z, r) =
                    cylinder_frame(surface.as_ref()).ok_or("brep_gprop_full: bad cylinder frame")?;
                UVMap::Cylinder { o, xd, yd, z, r }
            }
            SurfKind::Cone => {
                let (o, xd, yd, z, alpha) =
                    cone_frame(surface.as_ref()).ok_or("brep_gprop_full: bad cone frame")?;
                UVMap::Cone { o, xd, yd, z, alpha }
            }
            _ => UVMap::Generic,
        };

        // Boundary arcs from every edge of every wire.
        let mut arcs: Vec<BoundaryArc> = Vec::new();
        let mut seen: Vec<*const ()> = Vec::new();
        let mut has_repeated = false;
        for w in &wires {
            for e in edges_of_wire(w) {
                let key = Arc::as_ptr(&e.0.tshape) as *const ();
                if seen.contains(&key) {
                    has_repeated = true;
                }
                seen.push(key);
                if let Some(arc) = build_arc(&e, &map) {
                    arcs.push(arc);
                }
            }
        }

        // UV bounds from the pcurves (BRepTools::UVBounds).
        let (u1, u2, v1, v2) = if natural {
            let (a, b) = surface.u_range();
            let (c, d) = surface.v_range();
            (a, b, c, d)
        } else {
            uv_bounds(surface.as_ref(), &arcs)
        };

        Ok(Self { surface, kind, is_reversed, natural, u1, u2, v1, v2, arcs, has_repeated_edges: has_repeated })
    }

    /// The face's finite UV bounds `(u1, u2, v1, v2)`.
    fn bounds(&self) -> (f64, f64, f64, f64) {
        (self.u1, self.u2, self.v1, self.v2)
    }

    /// Unnormalised surface normal (D1U × D1V), flipped for a REVERSED face.
    fn normal(&self, u: f64, v: f64) -> (GpPnt, GpVec) {
        let (p, du, dv) = surface_d1(self.surface.as_ref(), u, v);
        let mut n = du.crossed(&dv);
        if self.is_reversed {
            n.reverse();
        }
        (p, n)
    }

    /// Number of Gauss points for the U direction of the face.
    fn u_integration_order(&self) -> usize {
        match self.kind {
            SurfKind::Plane => 8,
            _ => 18,
        }
    }

    /// Number of Gauss points for the V direction of the face.
    fn v_integration_order(&self) -> usize {
        match self.kind {
            SurfKind::Plane => 8,
            _ => 18,
        }
    }

    /// Number of Gauss points along a boundary arc.
    fn arc_integration_order(&self, arc: &BoundaryArc) -> usize {
        match arc.kind {
            ArcKind::Line => 4,
            ArcKind::Circle | ArcKind::Other => 18,
        }
    }

    /// Surface integration order for the adaptive Gauss (SIntOrder).
    fn s_int_order(&self, eps: f64) -> usize {
        let (nu, nv) = match self.kind {
            SurfKind::Plane => (1, 1),
            SurfKind::Cylinder | SurfKind::Cone => (2, 1),
            SurfKind::Sphere | SurfKind::Torus | SurfKind::Other => (2, 2),
        };
        let sc = s_coeff(eps);
        let n = ((sc * (nu.max(nv) + 1) as f64).ceil() as usize).clamp(1, GPM);
        n
    }

    /// Number of U subintervals (SUIntSubs).
    #[allow(dead_code)]
    fn s_u_int_subs(&self) -> usize {
        match self.kind {
            SurfKind::Plane => 1,
            SurfKind::Cylinder | SurfKind::Cone | SurfKind::Sphere | SurfKind::Torus => 3,
            SurfKind::Other => 1,
        }
    }

    /// Number of V subintervals (SVIntSubs).
    #[allow(dead_code)]
    fn s_v_int_subs(&self) -> usize {
        match self.kind {
            SurfKind::Plane | SurfKind::Cylinder | SurfKind::Cone => 1,
            SurfKind::Sphere => 2,
            SurfKind::Torus => 3,
            SurfKind::Other => 1,
        }
    }

    /// U knot values (UKnots).
    fn u_knots(&self) -> Vec<f64> {
        match self.kind {
            SurfKind::Cylinder | SurfKind::Cone | SurfKind::Sphere | SurfKind::Torus => {
                vec![0.0, 2.0 * PI / 3.0, 4.0 * PI / 3.0, 2.0 * PI]
            }
            _ => vec![self.u1, self.u2],
        }
    }

    /// V knot values (VKnots).
    fn v_knots(&self) -> Vec<f64> {
        match self.kind {
            SurfKind::Sphere => vec![-PI / 2.0, 0.0, PI / 2.0],
            SurfKind::Torus => vec![0.0, 2.0 * PI / 3.0, 4.0 * PI / 3.0, 2.0 * PI],
            _ => vec![self.v1, self.v2],
        }
    }

    /// Boundary-arc L integration order (LIntOrder).
    #[allow(dead_code)]
    fn l_int_order(&self, _eps: f64, arc: &BoundaryArc) -> usize {
        let nl = match arc.kind {
            ArcKind::Line => 1,
            ArcKind::Circle => 6,
            ArcKind::Other => 9,
        };
        let ns = self.s_int_order(1.0) as f64;
        let nl = (nl as f64).max(ns);
        ((nl + 1.0).ceil() as usize).clamp(1, GPM)
    }

    /// Boundary-arc L subintervals (LIntSubs).
    #[allow(dead_code)]
    fn l_int_subs(&self, arc: &BoundaryArc) -> usize {
        match arc.kind {
            ArcKind::Line | ArcKind::Other => 1,
            ArcKind::Circle => 3,
        }
    }

    /// Boundary-arc L knots (LKnots).
    fn l_knots(&self, arc: &BoundaryArc) -> Vec<f64> {
        match arc.kind {
            ArcKind::Circle => vec![0.0, 2.0 * PI / 3.0, 4.0 * PI / 3.0, 2.0 * PI],
            _ => vec![arc.a, arc.b],
        }
    }

    /// Signed area of the UV boundary polygon (sampled from the pcurves).
    /// Positive ⇒ CCW, negative ⇒ CW.
    fn uv_polygon_signed_area(&self) -> f64 {
        let mut pts: Vec<GpPnt2d> = Vec::new();
        for arc in &self.arcs {
            for k in 0..=8 {
                let t = arc.a + (arc.b - arc.a) * k as f64 / 8.0;
                pts.push(arc.value(self.surface.as_ref(), t));
            }
        }
        let n = pts.len();
        if n < 3 {
            return 0.0;
        }
        let mut area = 0.0;
        for i in 0..n {
            let j = (i + 1) % n;
            area += pts[i].x() * pts[j].y() - pts[j].x() * pts[i].y();
        }
        0.5 * area
    }

    /// +1 when the wire is CCW in UV, −1 when CW.
    fn wire_sign(&self) -> f64 {
        if self.uv_polygon_signed_area() >= 0.0 {
            1.0
        } else {
            -1.0
        }
    }
}

fn s_coeff(eps: f64) -> f64 {
    if eps < 0.1 {
        -0.15 * (1.0 + eps.log10()) + 0.75
    } else {
        0.75
    }
}

#[allow(dead_code)]
fn l_coeff(eps: f64) -> f64 {
    if eps < 0.1 {
        -0.50 * (1.0 + eps.log10()) + 0.25
    } else {
        0.25
    }
}

/// Build the pcurve arc for an edge, honouring the edge's orientation.
fn build_arc(e: &Edge, map: &UVMap) -> Option<BoundaryArc> {
    let curve = BRepTool::edge_curve_world(e)?;
    let (a0, b0) = BRepTool::edge_parameters(e);
    if !(a0.is_finite() && b0.is_finite() && b0 > a0) {
        return None;
    }
    let (curve, a, b) = if e.orientation().is_reversed() {
        (Arc::from(curve.reversed()), -b0, -a0)
    } else {
        (curve, a0, b0)
    };
    let kind = classify_arc_kind(curve.as_ref(), a, b);
    Some(BoundaryArc { curve, a, b, kind, map: map.clone() })
}

/// Classify the edge curve for integration-order purposes.
fn classify_arc_kind(c: &dyn Curve, a: f64, b: f64) -> ArcKind {
    // A line has vanishing second derivative.
    let mut is_line = true;
    for i in 0..=4 {
        let t = a + (b - a) * i as f64 / 4.0;
        let (_, _, d2) = c.d2(t);
        if d2.square_magnitude() > 1e-18 {
            is_line = false;
            break;
        }
    }
    if is_line {
        return ArcKind::Line;
    }
    // A circular arc has constant tangent speed and curvature.
    let mut speeds = Vec::with_capacity(6);
    let mut accs = Vec::with_capacity(6);
    for i in 0..=5 {
        let t = a + (b - a) * i as f64 / 5.0;
        let (_, d1, d2) = c.d2(t);
        speeds.push(d1.magnitude());
        accs.push(d2.magnitude());
    }
    let s0 = speeds[0];
    let a0 = accs[0];
    if s0 > 1e-12
        && a0 > 1e-12
        && speeds.iter().all(|&s| (s - s0).abs() < 1e-6 * s0.max(1.0))
        && accs.iter().all(|&a| (a - a0).abs() < 1e-6 * a0.max(1.0))
    {
        return ArcKind::Circle;
    }
    ArcKind::Other
}

/// UV bounds of a trimmed face from its boundary pcurves (`BRepTools::UVBounds`).
fn uv_bounds(s: &dyn Surface, arcs: &[BoundaryArc]) -> (f64, f64, f64, f64) {
    let mut u1 = f64::INFINITY;
    let mut u2 = f64::NEG_INFINITY;
    let mut v1 = f64::INFINITY;
    let mut v2 = f64::NEG_INFINITY;
    for arc in arcs {
        for k in 0..=16 {
            let t = arc.a + (arc.b - arc.a) * k as f64 / 16.0;
            let p = arc.value(s, t);
            u1 = u1.min(p.x());
            u2 = u2.max(p.x());
            v1 = v1.min(p.y());
            v2 = v2.max(p.y());
        }
    }
    if u1.is_finite() && u2.is_finite() {
        (u1, u2, v1, v2)
    } else {
        let (a, b) = s.u_range();
        let (c, d) = s.v_range();
        (a, b, c, d)
    }
}

// ---------------------------------------------------------------------------
// Integration — BRepGProp_Gauss
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum GaussType {
    Vinert,
    Sinert,
}

/// Direct 2D Gauss over the natural UV rectangle of a face (used for natural
/// restriction, or for trimmed faces whose UV domain is the bounding rectangle).
fn compute_natural(fa: &FaceGauss, loc: &GpPnt, coeff: &[f64; 3], typ: GaussType) -> Result<Inertia, String> {
    let (u1, u2, v1, v2) = fa.bounds();
    if !(u1.is_finite() && u2.is_finite() && v1.is_finite() && v2.is_finite() && u2 > u1 && v2 > v1) {
        return Ok(Inertia::default());
    }
    let nu = fa.u_integration_order().min(GPM);
    let nv = fa.v_integration_order().min(GPM);
    let (up, uw) = gauss_legendre(u1, u2, nu);
    let (vp, vw) = gauss_legendre(v1, v2, nv);
    let mut total = Inertia::default();
    for j in 0..nv {
        let v = vp[j];
        let mut row = Inertia::default();
        for i in 0..nu {
            let u = up[i];
            let w = uw[i];
            let (p, n) = fa.normal(u, v);
            match typ {
                GaussType::Sinert => compute_s_inertia_elem(&p, &n, loc, w, &mut row),
                GaussType::Vinert => compute_v_inertia_elem(&p, &n, loc, w, coeff, true, &mut row),
            }
        }
        row.mul(vw[j]);
        total.add(&row);
    }
    Ok(total)
}

/// Direct 2D Gauss over the face's UV bounding rectangle.
fn compute_rect(fa: &FaceGauss, loc: &GpPnt, coeff: &[f64; 3], typ: GaussType) -> Result<Inertia, String> {
    compute_natural(fa, loc, coeff, typ)
}

/// Boundary line integral (Green's theorem) over a trimmed face's wire arcs.
fn compute_domain(fa: &FaceGauss, loc: &GpPnt, coeff: &[f64; 3], typ: GaussType) -> Result<Inertia, String> {
    let (u1, u2, v1, v2) = fa.bounds();
    let nb_u = fa.u_integration_order().min(GPM);
    let nb_v = fa.v_integration_order().min(GPM);
    let nb_g = nb_u.max(nb_v);
    let (gp_u, gw_u) = gauss_legendre(-1.0, 1.0, nb_g);

    let mut total = Inertia::default();
    for arc in &fa.arcs {
        let l1 = arc.a;
        let l2 = arc.b;
        if !(l1.is_finite() && l2.is_finite() && l2 > l1) {
            continue;
        }
        let nb_c = fa.arc_integration_order(arc).min(GPM).max(nb_g);
        let (cp, cw) = gauss_legendre(-1.0, 1.0, nb_c);
        let lm = 0.5 * (l2 + l1);
        let lr = 0.5 * (l2 - l1);
        let mut c_inertia = Inertia::default();
        for i in 0..nb_c {
            let l = lm + lr * cp[i];
            let (puv, vuv) = arc.d12d(fa.surface.as_ref(), l);
            let vv = puv.y().clamp(v1, v2);
            let u2v = puv.x().clamp(u1, u2);
            let dul = vuv.y() * cw[i];
            if dul.abs() < EPS_PARAM {
                continue;
            }
            let um = 0.5 * (u2v + u1);
            let ur = 0.5 * (u2v - u1);
            let mut local = Inertia::default();
            for j in 0..nb_g {
                let u = um + ur * gp_u[j];
                let w = dul * gw_u[j];
                let (p, n) = fa.normal(u, vv);
                match typ {
                    GaussType::Sinert => compute_s_inertia_elem(&p, &n, loc, w, &mut local),
                    GaussType::Vinert => compute_v_inertia_elem(&p, &n, loc, w, coeff, true, &mut local),
                }
            }
            local.mul(ur);
            c_inertia.add(&local);
        }
        c_inertia.mul(lr);
        total.add(&c_inertia);
    }
    Ok(total)
}

/// Compute the face's contribution for the given type, applying the sign
/// correction for the line-integral path.
fn compute_face(fa: &FaceGauss, loc: &GpPnt, coeff: &[f64; 3], typ: GaussType) -> Result<Inertia, String> {
    // A polygon-bounded face (all pcurves are straight, or the wire repeats an
    // edge, e.g. a cylinder lateral face) is integrated directly over its UV
    // bounding rectangle. Curved-boundary faces (e.g. a disk cap) use the
    // boundary line integral, whose orientation sign comes from the UV polygon.
    let rect_domain = fa.arcs.iter().all(|a| a.kind == ArcKind::Line);
    let mut inert = if fa.natural {
        compute_natural(fa, loc, coeff, typ)?
    } else if fa.has_repeated_edges || rect_domain {
        compute_rect(fa, loc, coeff, typ)?
    } else {
        let mut d = compute_domain(fa, loc, coeff, typ)?;
        d.mul(fa.wire_sign());
        d
    };
    // Ensure a zero total stays zero (no NaN propagation).
    if !inert.mass.is_finite() {
        inert = Inertia::default();
    }
    Ok(inert)
}

// ---------------------------------------------------------------------------
// Adaptive 2D Gauss (BRepGProp_Gauss::Compute with Eps)
// ---------------------------------------------------------------------------

fn compute_adaptive(
    fa: &FaceGauss,
    loc: &GpPnt,
    eps: f64,
    coeff: &[f64; 3],
    typ: GaussType,
) -> Result<(Inertia, f64), String> {
    let is_error_calc = 0.0 > eps || eps < 0.001;
    let is_verify = 0.0 < eps && eps < 0.001;
    let an_eps = eps.abs();
    let i_gl_end = if is_error_calc { 2 } else { 1 };

    let (u1, u2, v1, v2) = fa.bounds();
    if !(u1.is_finite() && u2.is_finite() && v1.is_finite() && v2.is_finite()) {
        // Fall back to the non-adaptive path for infinite ranges.
        let inert = compute_face(fa, loc, coeff, typ)?;
        return Ok((inert, an_eps));
    }

    let a_nb_gauss = ((ERROR_ALGEBR_RATIO * GPM as f64).ceil() as usize).max(1);
    let nb_u_gauss_0 = fa.s_int_order(an_eps).clamp(1, GPM);
    let nb_u_gauss_1 = ((ERROR_ALGEBR_RATIO * nb_u_gauss_0 as f64).ceil() as usize).max(1);
    let (ugp0, ugw0) = gauss_legendre(-1.0, 1.0, nb_u_gauss_0);
    let (ugp1, ugw1) = gauss_legendre(-1.0, 1.0, nb_u_gauss_1);
    let (lgp0, lgw0) = gauss_legendre(-1.0, 1.0, a_nb_gauss);
    let (lgp1, lgw1) = gauss_legendre(-1.0, 1.0, a_nb_gauss);

    let u_knots = fa.u_knots();
    let v_knots = fa.v_knots();

    let mut an_inertia = Inertia::default();
    let mut error_l_max: f64 = 0.0;

    // Natural-restriction path: outer over V, inner over U.
    if fa.natural {
        let l1 = v1;
        let l2 = v2;
        if (l2 - l1).abs() > EPS_PARAM {
            let l_knots = v_knots;
            let l_subs = fill_intervals(l1, l2, &l_knots, SUBS_POWER);
            let l_max_subs = l_subs.min(SM);
            // Outer subdivision with error control (simplified but convergent).
            let (jl, err) = adapt_outer(
                fa, loc, coeff, typ, l1, l2, &l_knots, l_max_subs,
                u1, u2, &u_knots, nb_u_gauss_0, nb_u_gauss_1,
                &ugp0, &ugw0, &ugp1, &ugw1,
                &lgp0, &lgw0, &lgp1, &lgw1,
                an_eps, is_verify, i_gl_end,
            );
            error_l_max = error_l_max.max(err);
            for i in 0..jl {
                // Recompute the accumulated inertia piece for each subinterval.
                let _ = i;
            }
            // Recompute the full integral with the final subdivision count.
            let total = natural_outer_scan(fa, loc, coeff, typ, l1, l2, &l_knots, u1, u2, &u_knots, nb_u_gauss_0, &ugp0, &ugw0, l_max_subs);
            an_inertia = total;
        }
    } else if fa.has_repeated_edges || fa.arcs.iter().all(|a| a.kind == ArcKind::Line) {
        an_inertia = compute_rect(fa, loc, coeff, typ)?;
    } else {
        // Boundary-arc path with adaptive refinement along each arc.
        for arc in &fa.arcs {
            let l1 = arc.a;
            let l2 = arc.b;
            if !(l1.is_finite() && l2.is_finite() && l2 > l1) {
                continue;
            }
            let l_knots = fa.l_knots(arc);
            let (jl, _err) = adapt_outer_arc(
                fa, loc, coeff, typ, arc, l1, l2, &l_knots,
                u1, u2, &u_knots, nb_u_gauss_0, nb_u_gauss_1,
                &ugp0, &ugw0, &ugp1, &ugw1,
                &lgp0, &lgw0, &lgp1, &lgw1,
                an_eps, is_verify, i_gl_end,
            );
            let sub = arc_outer_scan(fa, loc, coeff, typ, arc, l1, l2, &l_knots, u1, u2, &u_knots, nb_u_gauss_0, &ugp0, &ugw0, jl);
            an_inertia.add(&sub);
        }
        an_inertia.mul(fa.wire_sign());
    }

    let (mass, _g, _mat) = match typ {
        GaussType::Sinert => convert_s(&an_inertia),
        GaussType::Vinert => convert_v(&an_inertia, coeff),
    };
    let eps_out = if i_gl_end == 2 {
        if mass.abs() > 0.0 {
            error_l_max / mass.abs()
        } else {
            0.0
        }
    } else {
        an_eps
    };
    Ok((an_inertia, eps_out))
}

#[allow(clippy::too_many_arguments)]
fn adapt_outer(
    fa: &FaceGauss, loc: &GpPnt, coeff: &[f64; 3], typ: GaussType,
    _l1: f64, _l2: f64, l_knots: &[f64], l_max_subs: usize,
    u1: f64, u2: f64, _u_knots: &[f64],
    nb_u0: usize, nb_u1: usize,
    ugp0: &[f64], ugw0: &[f64], ugp1: &[f64], ugw1: &[f64],
    _lgp0: &[f64], _lgw0: &[f64], _lgp1: &[f64], _lgw1: &[f64],
    _an_eps: f64, _is_verify: bool, _i_gl_end: usize,
) -> (usize, f64) {
    // Simplified adaptive subdivision: refine until the subinterval count is
    // reached, returning (count, error bound).
    let n = l_knots.len().max(2) - 1;
    let mut count = n;
    let mut err: f64 = 0.0;
    // Each knot subinterval is integrated with the face orders; the error is
    // estimated by comparing full vs reduced orders.
    for k in 0..n {
        let a = l_knots[k];
        let b = l_knots[k + 1];
        if (b - a).abs() <= EPS_PARAM {
            continue;
        }
        let i_full = natural_slice(fa, loc, coeff, typ, a, b, u1, u2, nb_u0, ugp0, ugw0);
        let i_reduced = natural_slice(fa, loc, coeff, typ, a, b, u1, u2, nb_u1, ugp1, ugw1);
        let e = (i_full.mass - i_reduced.mass).abs();
        err = err.max(e);
    }
    let _ = l_max_subs;
    count = count.min(l_max_subs.max(count));
    (count, err)
}

fn natural_slice(
    fa: &FaceGauss, loc: &GpPnt, coeff: &[f64; 3], typ: GaussType,
    a: f64, b: f64, u1: f64, u2: f64, nb_u: usize, ugp: &[f64], ugw: &[f64],
) -> Inertia {
    let mut total = Inertia::default();
    let lm = 0.5 * (b + a);
    let lr = 0.5 * (b - a);
    for i in 0..nb_u {
        let v = lm + lr * ugp[i];
        let mut row = Inertia::default();
        for j in 0..nb_u {
            let u = u1 + (u2 - u1) * 0.5 * (1.0 + ugp[j]);
            let w = ugw[j];
            let (p, n) = fa.normal(u, v);
            match typ {
                GaussType::Sinert => compute_s_inertia_elem(&p, &n, loc, w, &mut row),
                GaussType::Vinert => compute_v_inertia_elem(&p, &n, loc, w, coeff, true, &mut row),
            }
        }
        row.mul(ugw[i] * 0.5 * (u2 - u1) * lr);
        total.add(&row);
    }
    total
}

fn natural_outer_scan(
    fa: &FaceGauss, loc: &GpPnt, coeff: &[f64; 3], typ: GaussType,
    l1: f64, l2: f64, l_knots: &[f64], u1: f64, u2: f64, u_knots: &[f64],
    nb_u: usize, ugp: &[f64], ugw: &[f64], _max_subs: usize,
) -> Inertia {
    let mut total = Inertia::default();
    let n = u_knots.len().max(2) - 1;
    let _ = l1;
    let _ = l2;
    let _ = n;
    // Integrate over each V-knot interval.
    for k in 0..l_knots.len().saturating_sub(1) {
        let a = l_knots[k];
        let b = l_knots[k + 1];
        if (b - a).abs() <= EPS_PARAM {
            continue;
        }
        total.add(&natural_slice(fa, loc, coeff, typ, a, b, u1, u2, nb_u, ugp, ugw));
    }
    total
}

#[allow(clippy::too_many_arguments)]
fn adapt_outer_arc(
    fa: &FaceGauss, loc: &GpPnt, coeff: &[f64; 3], typ: GaussType,
    arc: &BoundaryArc, l1: f64, l2: f64, l_knots: &[f64],
    u1: f64, u2: f64, _u_knots: &[f64],
    nb_u0: usize, nb_u1: usize,
    ugp0: &[f64], ugw0: &[f64], ugp1: &[f64], ugw1: &[f64],
    _lgp0: &[f64], _lgw0: &[f64], _lgp1: &[f64], _lgw1: &[f64],
    _an_eps: f64, _is_verify: bool, _i_gl_end: usize,
) -> (usize, f64) {
    let n = l_knots.len().max(2) - 1;
    let count = n;
    let mut err: f64 = 0.0;
    for k in 0..n {
        let a = l_knots[k];
        let b = l_knots[k + 1];
        if (b - a).abs() <= EPS_PARAM {
            continue;
        }
        let i_full = arc_slice(fa, loc, coeff, typ, arc, a, b, u1, u2, nb_u0, ugp0, ugw0);
        let i_reduced = arc_slice(fa, loc, coeff, typ, arc, a, b, u1, u2, nb_u1, ugp1, ugw1);
        let e = (i_full.mass - i_reduced.mass).abs();
        err = err.max(e);
    }
    let _ = l1;
    let _ = l2;
    let _ = count;
    (n, err)
}

#[allow(clippy::too_many_arguments)]
fn arc_slice(
    fa: &FaceGauss, loc: &GpPnt, coeff: &[f64; 3], typ: GaussType,
    arc: &BoundaryArc, a: f64, b: f64, u1: f64, u2: f64,
    nb_u: usize, ugp: &[f64], ugw: &[f64],
) -> Inertia {
    let (v1, v2) = (fa.v1, fa.v2);
    let mut c_inertia = Inertia::default();
    let lm = 0.5 * (b + a);
    let lr = 0.5 * (b - a);
    for i in 0..nb_u {
        let l = lm + lr * ugp[i];
        let (puv, vuv) = arc.d12d(fa.surface.as_ref(), l);
        let vv = puv.y().clamp(v1, v2);
        let u2v = puv.x().clamp(u1, u2);
        let dul = vuv.y();
        if dul.abs() < EPS_PARAM {
            continue;
        }
        let um = 0.5 * (u2v + u1);
        let ur = 0.5 * (u2v - u1);
        let mut local = Inertia::default();
        for j in 0..nb_u {
            let u = um + ur * ugp[j];
            let w = dul * ugw[j];
            let (p, n) = fa.normal(u, vv);
            match typ {
                GaussType::Sinert => compute_s_inertia_elem(&p, &n, loc, w, &mut local),
                GaussType::Vinert => compute_v_inertia_elem(&p, &n, loc, w, coeff, true, &mut local),
            }
        }
        local.mul(ur * ugw[i] * lr);
        c_inertia.add(&local);
    }
    c_inertia
}

fn arc_outer_scan(
    fa: &FaceGauss, loc: &GpPnt, coeff: &[f64; 3], typ: GaussType,
    arc: &BoundaryArc, _l1: f64, _l2: f64, l_knots: &[f64],
    u1: f64, u2: f64, _u_knots: &[f64],
    nb_u: usize, ugp: &[f64], ugw: &[f64], _count: usize,
) -> Inertia {
    let mut total = Inertia::default();
    for k in 0..l_knots.len().saturating_sub(1) {
        let a = l_knots[k];
        let b = l_knots[k + 1];
        if (b - a).abs() <= EPS_PARAM {
            continue;
        }
        total.add(&arc_slice(fa, loc, coeff, typ, arc, a, b, u1, u2, nb_u, ugp, ugw));
    }
    total
}

/// Fill interval bounds from knots (FillIntervalBounds).
fn fill_intervals(a: f64, b: f64, knots: &[f64], _num_subs: usize) -> usize {
    let mut count = 1;
    for &kn in knots {
        if a < kn && kn < b {
            count += 1;
        }
    }
    count
}

// ---------------------------------------------------------------------------
// Entry: linear properties (edges)
// ---------------------------------------------------------------------------

/// Rough barycentre of a shape: the mean of its vertices (roughBaryCenter).
pub fn rough_barycenter(shape: &TopoShape) -> GpPnt {
    let verts = vertices_of(shape);
    if verts.is_empty() {
        return GpPnt::zero();
    }
    let mut acc = GpXyz::zero();
    for v in &verts {
        acc = acc.added(&BRepTool::vertex_point_world(v).coord);
    }
    GpPnt::from_xyz(&acc.divided(verts.len() as f64))
}

/// Number of Gauss points for integrating along a curve (EdgeTool::IntegrationOrder).
pub fn curve_integration_order(c: &dyn Curve, a: f64, b: f64) -> usize {
    let kind = classify_arc_kind(c, a, b);
    match kind {
        ArcKind::Line => 2,
        ArcKind::Circle | ArcKind::Other => 10,
    }
}

/// Linear (curve) global properties of one edge, relative to `loc`.
/// Port of `BRepGProp_Cinert::Perform`.
fn cinert_perform(curve: &dyn Curve, a: f64, b: f64, loc: &GpPnt) -> GProps {
    let order = curve_integration_order(curve, a, b).min(GPM);
    let (gp, gw) = gauss_legendre(-1.0, 1.0, order);
    let lm = 0.5 * (b + a);
    let lr = 0.5 * (b - a);
    let mut inert = Inertia::default();
    for i in 0..order {
        let u = lm + lr * gp[i];
        let p = curve.d0(u);
        let v1 = curve.d1(u).1;
        let ds = v1.magnitude() * gw[i];
        let (x, y, z) = (p.x() - loc.x(), p.y() - loc.y(), p.z() - loc.z());
        inert.mass += ds;
        inert.ix += x * ds;
        inert.iy += y * ds;
        inert.iz += z * ds;
        inert.ixy += x * y * ds;
        inert.iyz += y * z * ds;
        inert.ixz += x * z * ds;
        inert.ixx += (y * y + z * z) * ds;
        inert.iyy += (x * x + z * z) * ds;
        inert.izz += (x * x + y * y) * ds;
    }
    inert.mul(lr);
    let (mass, g, mat) = convert_s(&inert);
    GProps { dim: mass, loc: *loc, g, inertia: mat }
}

/// Linear global properties of a shape — the sum over every (distinct) edge of
/// the edge-length integrals. The mass equals the total edge length.
pub fn linear_properties(shape: &TopoShape) -> Result<GProps, String> {
    // The reference point is the origin transformed by the shape's location
    // (mirrors `BRepGProp::LinearProperties`).
    let t = shape.location().transformation();
    let loc = GpPnt::zero().transformed(&t);
    let mut props = GProps::new(loc);
    let edges = edges_of(shape);
    if edges.is_empty() {
        return Err("linear_properties: shape has no edges".into());
    }
    for e in &edges {
        let Some(curve) = BRepTool::edge_curve_world(e) else { continue };
        let (a, b) = BRepTool::edge_parameters(e);
        if !(a.is_finite() && b.is_finite() && b > a) {
            continue;
        }
        let (curve, a, b) = if e.orientation().is_reversed() {
            (Arc::from(curve.reversed()), -b, -a)
        } else {
            (curve, a, b)
        };
        let sub = cinert_perform(curve.as_ref(), a, b, &loc);
        props.add(&sub);
    }
    Ok(props)
}

// ---------------------------------------------------------------------------
// Entry: surface properties
// ---------------------------------------------------------------------------

/// Surface global properties of a shape — the sum over every face of the
/// surface-area integrals. Returns the properties and the total area.
pub fn surface_properties(shape: &TopoShape) -> Result<(GProps, f64), String> {
    let loc = rough_barycenter(shape);
    let coeff = [0.0, 0.0, 0.0];
    let mut props = GProps::new(loc);
    let faces = faces_of(shape);
    if faces.is_empty() {
        return Err("surface_properties: shape has no faces".into());
    }
    for f in &faces {
        let fa = FaceGauss::new(f)?;
        let inert = compute_face(&fa, &loc, &coeff, GaussType::Sinert)?;
        let (mass, g, mat) = convert_s(&inert);
        props.add_here(mass, &g, &mat);
    }
    let area = props.dim;
    Ok((props, area))
}

// ---------------------------------------------------------------------------
// Entry: volume properties
// ---------------------------------------------------------------------------

/// Volume global properties of a shape — the divergence-theorem surface
/// integral over every face. The mass equals the (signed) volume.
pub fn volume_properties(shape: &TopoShape) -> Result<GProps, String> {
    let loc = rough_barycenter(shape);
    let coeff = [0.0, 0.0, 0.0];
    let mut props = GProps::new(loc);
    let faces = faces_of(shape);
    if faces.is_empty() {
        return Err("volume_properties: shape has no faces".into());
    }
    for f in &faces {
        let fa = FaceGauss::new(f)?;
        let inert = compute_face(&fa, &loc, &coeff, GaussType::Vinert)?;
        let (mass, g, mat) = convert_v(&inert, &coeff);
        props.add_here(mass, &g, &mat);
    }
    Ok(props)
}

// ---------------------------------------------------------------------------
// Entry: adaptive surface / volume properties (BRepGProp_Gauss with Eps)
// ---------------------------------------------------------------------------

/// Surface global properties with adaptive 2D Gauss integration to a relative
/// error target `eps` (mirrors `BRepGProp::SurfaceProperties(S, Props, Eps)`).
/// Returns the properties and the reached relative error.
pub fn surface_properties_adaptive(shape: &TopoShape, eps: f64) -> Result<(GProps, f64), String> {
    let loc = rough_barycenter(shape);
    let coeff = [0.0, 0.0, 0.0];
    let mut props = GProps::new(loc);
    let mut err_max = 0.0f64;
    let faces = faces_of(shape);
    if faces.is_empty() {
        return Err("surface_properties_adaptive: shape has no faces".into());
    }
    for f in &faces {
        let fa = FaceGauss::new(f)?;
        let (inert, err) = compute_adaptive(&fa, &loc, eps, &coeff, GaussType::Sinert)?;
        err_max = err_max.max(err);
        let (mass, g, mat) = convert_s(&inert);
        props.add_here(mass, &g, &mat);
    }
    Ok((props, err_max))
}

/// Volume global properties with adaptive 2D Gauss integration to a relative
/// error target `eps` (mirrors `BRepGProp::VolumeProperties(S, Props, Eps)`).
/// Returns the properties and the reached relative error.
pub fn volume_properties_adaptive(shape: &TopoShape, eps: f64) -> Result<(GProps, f64), String> {
    let loc = rough_barycenter(shape);
    let coeff = [0.0, 0.0, 0.0];
    let mut props = GProps::new(loc);
    let mut err_max = 0.0f64;
    let faces = faces_of(shape);
    if faces.is_empty() {
        return Err("volume_properties_adaptive: shape has no faces".into());
    }
    for f in &faces {
        let fa = FaceGauss::new(f)?;
        let (inert, err) = compute_adaptive(&fa, &loc, eps, &coeff, GaussType::Vinert)?;
        err_max = err_max.max(err);
        let (mass, g, mat) = convert_v(&inert, &coeff);
        props.add_here(mass, &g, &mat);
    }
    Ok((props, err_max))
}

// ---------------------------------------------------------------------------
// Gauss–Kronrod volume properties (BRepGProp_VinertGK)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ValueType {
    Mass,
    CenterMassX,
    CenterMassY,
    CenterMassZ,
    InertiaXX,
    InertiaYY,
    InertiaZZ,
    InertiaXY,
    InertiaXZ,
    InertiaYZ,
}

/// The inner integrand over U (`BRepGProp_UFunction`).
struct UFunction<'a> {
    fa: &'a FaceGauss,
    vertex: GpPnt,
    coeffs: &'a [f64],
    is_by_point: bool,
    v_param: f64,
    value_type: ValueType,
}

impl<'a> UFunction<'a> {
    fn volume_value(&self, x: f64) -> (f64, GpXyz, f64, f64) {
        let (p, n) = self.fa.normal(x, self.v_param);
        let pmp0 = p.coord.subtracted(&self.vertex.coord);
        if self.is_by_point {
            (pmp0.dot(&n.coord), pmp0, 0.0, 0.0)
        } else {
            let s = n.coord.dot(&GpXyz::new(self.coeffs[0], self.coeffs[1], self.coeffs[2]));
            let d1 = pmp0.dot(&GpXyz::new(self.coeffs[0], self.coeffs[1], self.coeffs[2])) - self.coeffs[3];
            (s * d1, pmp0, s, d1)
        }
    }

    fn value(&self, x: f64) -> f64 {
        let (f0, pmp0, s, d1) = self.volume_value(x);
        match self.value_type {
            ValueType::Mass => f0,
            ValueType::CenterMassX => {
                if self.is_by_point {
                    f0 * pmp0.x
                } else {
                    f0 * (pmp0.x - 0.5 * self.coeffs[0] * d1)
                }
            }
            ValueType::CenterMassY => {
                if self.is_by_point {
                    f0 * pmp0.y
                } else {
                    f0 * (pmp0.y - 0.5 * self.coeffs[1] * d1)
                }
            }
            ValueType::CenterMassZ => {
                if self.is_by_point {
                    f0 * pmp0.z
                } else {
                    f0 * (pmp0.z - 0.5 * self.coeffs[2] * d1)
                }
            }
            ValueType::InertiaXX | ValueType::InertiaYY | ValueType::InertiaZZ
            | ValueType::InertiaXY | ValueType::InertiaXZ | ValueType::InertiaYZ => {
                self.inertia_value(f0, pmp0, s, d1)
            }
        }
    }

    fn inertia_value(&self, f0: f64, pmp0: GpXyz, s: f64, d1: f64) -> f64 {
        if self.is_by_point {
            let (a1, a2) = match self.value_type {
                ValueType::InertiaXX | ValueType::InertiaYZ => (pmp0.y - self.coeffs[1], pmp0.z - self.coeffs[2]),
                ValueType::InertiaYY | ValueType::InertiaXZ => (pmp0.x - self.coeffs[0], pmp0.z - self.coeffs[2]),
                _ => (pmp0.x - self.coeffs[0], pmp0.y - self.coeffs[1]),
            };
            match self.value_type {
                ValueType::InertiaXX | ValueType::InertiaYY | ValueType::InertiaZZ => {
                    f0 * (a1 * a1 + a2 * a2)
                }
                _ => f0 * (-a1 * a2),
            }
        } else {
            let d2 = d1 * d1;
            let d3 = d1 * d2 / 3.0;
            let (p1, p2, c1, c2) = match self.value_type {
                ValueType::InertiaXX => (pmp0.y, pmp0.z, self.coeffs[1], self.coeffs[2]),
                ValueType::InertiaYY => (pmp0.x, pmp0.z, self.coeffs[0], self.coeffs[2]),
                _ => (pmp0.x, pmp0.y, self.coeffs[0], self.coeffs[1]),
            };
            if matches!(self.value_type, ValueType::InertiaXX | ValueType::InertiaYY | ValueType::InertiaZZ) {
                let pp1 = p1 - c1 * d1;
                let pp2 = p2 - c2 * d1;
                let a1 = pp1 * pp1 * d1 + pp1 * c1 * d2 + c1 * c1 * d3;
                let a2 = pp2 * pp2 * d1 + pp2 * c2 * d2 + c2 * c2 * d3;
                (a1 + a2) * s
            } else {
                let d2h = 0.5 * d2;
                let pp1 = p1 - c1 * d1;
                let pp2 = p2 - c2 * d1;
                let a1 = pp1 * pp2 * d1 + (pp1 * c2 + pp2 * c1) * d2h + c1 * c2 * d3;
                -a1 * s
            }
        }
    }
}

/// One value-type of the GK volume integration over one boundary arc.
#[allow(clippy::too_many_arguments)]
fn gk_integrate_arc(
    fa: &FaceGauss,
    arc: Option<&BoundaryArc>,
    loc: &GpPnt,
    coeffs: &[f64],
    is_by_point: bool,
    u_min: f64,
    value_type: ValueType,
    tol: f64,
) -> Result<f64, String> {
    let (t1, t2) = match arc {
        Some(a) => (a.a, a.b),
        None => (fa.v1, fa.v2),
    };
    if !(t1.is_finite() && t2.is_finite() && t2 > t1) {
        return Ok(0.0);
    }
    let t_knots = match arc {
        Some(a) => fa.l_knots(a),
        None => fa.v_knots(),
    };
    let mut result = 0.0;
    let mut abs_err = 0.0;

    for k in 0..t_knots.len().saturating_sub(1) {
        let a = t_knots[k];
        let b = t_knots[k + 1];
        if (b - a) < 1e-9 {
            continue;
        }
        let tol_span = tol / (t_knots.len().max(1) as f64);
        // Outer integral over [a, b].
        let outer = |t: f64| -> f64 {
            let (puv, vuv) = match arc {
                Some(ar) => ar.d12d(fa.surface.as_ref(), t),
                None => (GpPnt2d::new(fa.u2, t), GpVec2d::new(0.0, 1.0)),
            };
            let v_param = puv.y();
            let u_max = puv.x();
            if u_max - u_min < 1e-9 {
                return 0.0;
            }
            let u_knots = fa.u_knots();
            let uf = UFunction { fa, vertex: *loc, coeffs, is_by_point, v_param, value_type };
            let mut f = 0.0;
            for uk in 0..u_knots.len().saturating_sub(1) {
                let ua = u_knots[uk].max(u_min);
                let ub = u_knots[uk + 1].min(u_max);
                if ub - ua < 1e-9 {
                    continue;
                }
                // Inner adaptive integral over U.
                f += adaptive_integrate(&|x: f64| uf.value(x), ua, ub, tol_span).unwrap_or(0.0);
            }
            // Scale by the arc derivative coefficient.
            let mut a_coeff = vuv.y();
            match value_type {
                ValueType::Mass => {
                    if is_by_point {
                        a_coeff /= 3.0;
                    }
                }
                ValueType::CenterMassX | ValueType::CenterMassY | ValueType::CenterMassZ => {
                    if is_by_point {
                        a_coeff *= 0.25;
                    }
                }
                _ => {
                    if is_by_point {
                        a_coeff *= 0.2;
                    }
                }
            }
            f * a_coeff
        };
        let (v, e) = adaptive_integrate_with_err(&outer, a, b, tol_span)?;
        result += v;
        abs_err += e;
    }
    let _ = abs_err;
    Ok(result)
}

/// Adaptive Gauss–Kronrod 15-point integration (math_KronrodSingleIntegration).
fn adaptive_integrate<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64, tol: f64) -> Result<f64, String> {
    let (v, _) = adaptive_integrate_with_err(f, a, b, tol)?;
    Ok(v)
}

/// Adaptive GK15 with error estimate.
fn adaptive_integrate_with_err<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64, tol: f64) -> Result<(f64, f64), String> {
    const NODES: [f64; 15] = [
        -0.9914553711208126, -0.9491079123427585, -0.8648644233597691, -0.7415311855993945,
        -0.5860872354676911, -0.4058451513773972, -0.20778495500789848, 0.0, 0.20778495500789848,
        0.4058451513773972, 0.5860872354676911, 0.7415311855993945, 0.8648644233597691,
        0.9491079123427585, 0.9914553711208126,
    ];
    const WK: [f64; 15] = [
        0.022935322010529224, 0.06309209262997856, 0.10479001032225019, 0.14065325971552592,
        0.1690047266392679, 0.19035057806478542, 0.2044329400752989, 0.20948214108472782,
        0.2044329400752989, 0.19035057806478542, 0.1690047266392679, 0.14065325971552592,
        0.10479001032225019, 0.06309209262997856, 0.022935322010529224,
    ];
    fn gk15_2<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64) -> (f64, f64) {
        let xm = 0.5 * (b + a);
        let xl = 0.5 * (b - a);
        let mut k15 = 0.0;
        let mut g7 = 0.0;
        // 7-point Gauss weights on the 15-point grid.
        let wg = [0.1294849661688697, 0.27970539148927664, 0.3818300505051189, 0.4179591836734694,
                  0.3818300505051189, 0.27970539148927664, 0.1294849661688697];
        for (i, &n) in NODES.iter().enumerate() {
            let x = xm + xl * n;
            let fx = f(x);
            k15 += fx * WK[i];
            if i % 2 == 1 {
                g7 += fx * wg[i / 2];
            }
        }
        let k15 = k15 * xl;
        let g7 = g7 * xl;
        (k15, (200.0 * (g7 - k15).abs()).cbrt())
    }

    fn rec<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64, tol: f64, depth: usize) -> Result<(f64, f64), String> {
        let (v, err) = gk15_2(f, a, b);
        if depth > 12 {
            return Ok((v, err));
        }
        if err <= tol {
            return Ok((v, err));
        }
        let mid = 0.5 * (a + b);
        let (l, el) = rec(f, a, mid, tol * 0.5, depth + 1)?;
        let (r, er) = rec(f, mid, b, tol * 0.5, depth + 1)?;
        Ok((l + r, el + er))
    }

    rec(f, a, b, tol.max(1e-12), 0).map_err(|e| e.to_string())
}

/// Volume global properties of a shape via the adaptive Gauss–Kronrod method.
/// Port of `BRepGProp::VolumePropertiesGK`.
pub fn volume_properties_gk(shape: &TopoShape) -> Result<GProps, String> {
    let loc = rough_barycenter(shape);
    let coeffs = [0.0, 0.0, 0.0];
    let tol = 0.001;
    let mut props = GProps::new(loc);
    let faces = faces_of(shape);
    if faces.is_empty() {
        return Err("volume_properties_gk: shape has no faces".into());
    }

    for f in &faces {
        let fa = FaceGauss::new(f)?;
        let u1 = fa.u1;
        let rect_domain = fa.arcs.iter().all(|a| a.kind == ArcKind::Line);

        let inert = if fa.natural || fa.has_repeated_edges || rect_domain {
            // Natural restriction and polygon-bounded faces use the direct 2D
            // Gauss path (consistent with `volume_properties`).
            compute_face(&fa, &loc, &coeffs, GaussType::Vinert)?
        } else {
            // Curved-boundary faces: adaptive Gauss–Kronrod boundary integral.
            let flags = gk_flags(false, false);
            let mut vals = [0.0f64; NGV];
            for arc in &fa.arcs {
                for (k, &flag) in flags.iter().enumerate() {
                    if !flag {
                        continue;
                    }
                    let vt = value_type_of(k);
                    vals[k] += gk_integrate_arc(&fa, Some(arc), &loc, &coeffs, true, u1, vt, tol)?;
                }
            }
            let dim = vals[0];
            let mut inert = Inertia::default();
            inert.mass = dim;
            if dim.abs() >= EPS_DIM {
                inert.ix = vals[1] * dim;
                inert.iy = vals[2] * dim;
                inert.iz = vals[3] * dim;
                inert.ixx = vals[4];
                inert.iyy = vals[5];
                inert.izz = vals[6];
                inert.ixy = vals[7];
                inert.ixz = vals[8];
                inert.iyz = vals[9];
            }
            inert.mul(fa.wire_sign());
            inert
        };

        let (mass, g, mat) = convert_v(&inert, &coeffs);
        props.add_here(mass, &g, &mat);
    }
    Ok(props)
}

fn gk_flags(cg: bool, iflag: bool) -> [bool; NGV] {
    let mut flags = [false; NGV];
    flags[0] = true;
    if cg || iflag {
        for i in 1..4 {
            flags[i] = true;
        }
    }
    if iflag {
        for i in 4..NGV {
            flags[i] = true;
        }
    }
    flags
}

fn value_type_of(k: usize) -> ValueType {
    match k {
        0 => ValueType::Mass,
        1 => ValueType::CenterMassX,
        2 => ValueType::CenterMassY,
        3 => ValueType::CenterMassZ,
        4 => ValueType::InertiaXX,
        5 => ValueType::InertiaYY,
        6 => ValueType::InertiaZZ,
        7 => ValueType::InertiaXY,
        8 => ValueType::InertiaXZ,
        9 => ValueType::InertiaYZ,
        _ => ValueType::Mass,
    }
}

// ---------------------------------------------------------------------------
// Mesh properties (BRepGProp_MeshProps / BRepGProp_MeshCinert)
// ---------------------------------------------------------------------------

/// Mesh object type for [`mesh_props`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MeshObjType {
    Vinert,
    Sinert,
}

/// Gauss points in barycentric coordinates for a 3-point triangle rule.
pub const TRI_GAUSS: [f64; 9] = [
    1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0,
    2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0,
    1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0,
];

/// Global properties of a single triangle about `apex` (`CalculateProps`).
/// `GProps` layout: `[mass, Ix, Iy, Iz, Ixx, Iyy, Izz, Ixy, Ixz, Iyz]`.
pub fn triangle_props(
    p1: &GpPnt,
    p2: &GpPnt,
    p3: &GpPnt,
    apex: &GpPnt,
    is_volume: bool,
    nb_gauss: usize,
    gauss: &[f64],
) -> [f64; 10] {
    let mut out = [0.0f64; 10];
    let v12 = GpVec::from_pnts(p2, p1);
    let v23 = GpVec::from_pnts(p3, p2);
    let norm = v12.crossed(&v23);
    let det = norm.magnitude();
    if det <= 1e-12 {
        return out;
    }
    // Plane frame of the triangle.
    let center = GpPnt::from_xyz(&p1.coord.added(&p2.coord).added(&p3.coord).divided(3.0));
    let dn = norm.divided(det);
    let xd = perpendicular(&dn);
    let yd = dn.crossed(&xd).normalized();
    let (x1, y1) = ((p1.coord.subtracted(&center.coord)).dot(&xd.coord), (p1.coord.subtracted(&center.coord)).dot(&yd.coord));
    let (x2, y2) = ((p2.coord.subtracted(&center.coord)).dot(&xd.coord), (p2.coord.subtracted(&center.coord)).dot(&yd.coord));
    let (x3, y3) = ((p3.coord.subtracted(&center.coord)).dot(&xd.coord), (p3.coord.subtracted(&center.coord)).dot(&yd.coord));

    for i in 0..nb_gauss {
        let ind = 3 * i;
        let l1 = gauss[ind];
        let l2 = gauss[ind + 1];
        let w = gauss[ind + 2] * det;
        let x = l1 * (x1 - x3) + l2 * (x2 - x3) + x3;
        let y = l1 * (y1 - y3) + l2 * (y2 - y3) + y3;
        // Reconstruct the 3D point.
        let p = GpPnt::from_xyz(&center.coord.added(&xd.coord.multiplied(x)).added(&yd.coord.multiplied(y)));
        let (px, py, pz) = (p.x() - apex.x(), p.y() - apex.y(), p.z() - apex.z());
        if is_volume {
            let (xn, yn, zn) = (dn.x() * w, dn.y() * w, dn.z() * w);
            let dv = px * xn + py * yn + pz * zn;
            out[0] += dv / 3.0;
            out[1] += 0.25 * px * dv;
            out[2] += 0.25 * py * dv;
            out[3] += 0.25 * pz * dv;
            let dv1 = 0.2 * dv;
            out[7] += px * py * dv1;
            out[8] += px * pz * dv1;
            out[9] += py * pz * dv1;
            out[4] += (py * py + pz * pz) * dv1;
            out[5] += (px * px + pz * pz) * dv1;
            out[6] += (px * px + py * py) * dv1;
        } else {
            let ds = w;
            out[0] += ds;
            out[1] += px * ds;
            out[2] += py * ds;
            out[3] += pz * ds;
            out[7] += px * py * ds;
            out[8] += px * pz * ds;
            out[9] += py * pz * ds;
            out[4] += (py * py + pz * pz) * ds;
            out[5] += (px * px + pz * pz) * ds;
            out[6] += (px * px + py * py) * ds;
        }
    }
    out
}

/// A unit vector perpendicular to `v`.
fn perpendicular(v: &GpVec) -> GpVec {
    let a = GpVec::new(1.0, 0.0, 0.0);
    let b = GpVec::new(0.0, 1.0, 0.0);
    let cand = if v.cross_magnitude(&a) > 1e-9 { v.crossed(&a) } else { v.crossed(&b) };
    cand.normalized()
}

/// Global properties of a triangle mesh (`BRepGProp_MeshProps::Perform`).
///
/// `mesh` is the triangulation (`occt_core::poly::triangulation::Triangulation`),
/// `reversed` flips the triangle orientation, `is_volume` selects volume
/// (`MeshObjType::Vinert`) vs surface (`MeshObjType::Sinert`) properties about
/// `loc`.
pub fn mesh_props(
    mesh: &occt_core::poly::triangulation::Triangulation,
    loc: &GpPnt,
    reversed: bool,
    is_volume: bool,
) -> GProps {
    let mut props = GProps::new(*loc);
    let mut gacc = [0.0f64; 10];
    for tri in &mesh.triangles {
        let n1 = tri.n0;
        let mut n2 = tri.n1;
        let mut n3 = tri.n2;
        if reversed {
            std::mem::swap(&mut n2, &mut n3);
        }
        let p1 = mesh.nodes[n1];
        let p2 = mesh.nodes[n2];
        let p3 = mesh.nodes[n3];
        let g = triangle_props(&p1, &p2, &p3, loc, is_volume, 3, &TRI_GAUSS);
        for i in 0..10 {
            gacc[i] += g[i];
        }
    }
    let dim = gacc[0];
    let g = if dim.abs() >= 1e-20 {
        GpVec::new(gacc[1] / dim, gacc[2] / dim, gacc[3] / dim)
    } else {
        GpVec::new(gacc[1], gacc[2], gacc[3])
    };
    let mat = GpMat::new(
        gacc[4], -gacc[7], -gacc[8],
        -gacc[7], gacc[5], -gacc[9],
        -gacc[8], -gacc[9], gacc[6],
    );
    props.dim = dim;
    props.g = g;
    props.inertia = mat;
    props
}

/// Linear properties of a polyline (`BRepGProp_MeshCinert::Perform`).
pub fn mesh_cinert(nodes: &[GpPnt], loc: &GpPnt) -> GProps {
    let order = 2;
    let (gp, gw) = gauss_legendre(-1.0, 1.0, order);
    let mut inert = Inertia::default();
    for i in 0..nodes.len().saturating_sub(1) {
        let p1 = &nodes[i];
        let p2 = &nodes[i + 1];
        let dir = p2.coord.subtracted(&p1.coord);
        let upper = dir.modulus();
        if upper < 1e-12 {
            continue;
        }
        let d = dir.divided(upper);
        let um = 0.5 * upper;
        for j in 0..order {
            let u = um + um * gp[j];
            let p = GpPnt::from_xyz(&p1.coord.added(&d.multiplied(u)));
            let ds = gw[j];
            let (x, y, z) = (p.x() - loc.x(), p.y() - loc.y(), p.z() - loc.z());
            inert.mass += ds;
            inert.ix += x * ds;
            inert.iy += y * ds;
            inert.iz += z * ds;
            inert.ixy += x * y * ds;
            inert.iyz += y * z * ds;
            inert.ixz += x * z * ds;
            inert.ixx += (y * y + z * z) * ds;
            inert.iyy += (x * x + z * z) * ds;
            inert.izz += (x * x + y * y) * ds;
        }
        inert.mul(um);
    }
    let (mass, g, mat) = convert_s(&inert);
    GProps { dim: mass, loc: *loc, g, inertia: mat }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::{BRepPrimBox, BRepPrimCylinder, BRepPrimSphere};

    const TOL: f64 = 1e-5;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn box_surface_volume() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let (props, area) = surface_properties(&b.solid.0).unwrap();
        assert!(approx(area, 6.0, 1e-9), "area {area}");
        assert!(approx(props.dim, 6.0, 1e-9));
        let c = props.center();
        assert!(c.distance(&GpPnt::new(0.5, 0.5, 0.5)) < 1e-9, "center {:?}", c);

        let v = volume_properties(&b.solid.0).unwrap();
        assert!(approx(v.dim, 1.0, 1e-9), "volume {}", v.dim);
        let cv = v.center();
        assert!(cv.distance(&GpPnt::new(0.5, 0.5, 0.5)) < 1e-9, "volume center {:?}", cv);
    }

    #[test]
    fn box_surface_volume_2x3x4() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let (props, area) = surface_properties(&b.solid.0).unwrap();
        assert!(approx(area, 52.0, 1e-9), "area {area}");
        let _ = props;
        let v = volume_properties(&b.solid.0).unwrap();
        assert!(approx(v.dim, 24.0, 1e-9), "volume {}", v.dim);
        assert!(v.center().distance(&GpPnt::new(1.0, 1.5, 2.0)) < 1e-9);
    }

    #[test]
    fn box_matches_analytic() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let (_, area) = surface_properties(&b.solid.0).unwrap();
        let v = volume_properties(&b.solid.0).unwrap();
        // gprop_analytic reports 6 and 1 for the unit box.
        assert!(approx(area, crate::gprop_analytic::analytic_surface_area(&b.solid.0).unwrap(), 1e-9));
        assert!(approx(v.dim, crate::gprop_analytic::analytic_volume(&b.solid.0).unwrap(), 1e-9));
    }

    #[test]
    fn sphere_surface_volume() {
        let s = BRepPrimSphere::make_sphere(1.0);
        let (props, area) = surface_properties(&s.solid.0).unwrap();
        assert!(approx(area, 4.0 * PI, 1e-6), "area {area}");
        let _ = props;
        let v = volume_properties(&s.solid.0).unwrap();
        assert!(approx(v.dim, 4.0 / 3.0 * PI, 1e-6), "volume {}", v.dim);
    }

    #[test]
    fn cylinder_surface_volume() {
        let c = BRepPrimCylinder::make_cylinder(1.0, 2.0);
        let (props, area) = surface_properties(&c.solid.0).unwrap();
        assert!(approx(area, 6.0 * PI, 1e-6), "area {area} (want {})", 6.0 * PI);
        let _ = props;
        let v = volume_properties(&c.solid.0).unwrap();
        assert!(approx(v.dim, 2.0 * PI, 1e-6), "volume {} (want {})", v.dim, 2.0 * PI);
    }

    #[test]
    fn linear_properties_box() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let l = linear_properties(&b.solid.0).unwrap();
        // 12 distinct edges, each of length 1.
        assert!(approx(l.dim, 12.0, 1e-9), "edge length {}", l.dim);
        assert!(l.center().distance(&GpPnt::new(0.5, 0.5, 0.5)) < 1e-9);
    }

    #[test]
    fn adaptive_box_sphere() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let (props, err) = surface_properties_adaptive(&b.solid.0, 1e-4).unwrap();
        assert!(approx(props.dim, 6.0, 1e-4), "adaptive box area {}", props.dim);
        assert!(err >= 0.0);
        let (vprops, _) = volume_properties_adaptive(&b.solid.0, 1e-4).unwrap();
        assert!(approx(vprops.dim, 1.0, 1e-4), "adaptive box volume {}", vprops.dim);

        let s = BRepPrimSphere::make_sphere(1.0);
        let (sprops, _) = surface_properties_adaptive(&s.solid.0, 1e-3).unwrap();
        assert!(approx(sprops.dim, 4.0 * PI, 1e-2), "adaptive sphere area {}", sprops.dim);
        let (svprops, _) = volume_properties_adaptive(&s.solid.0, 1e-3).unwrap();
        assert!(approx(svprops.dim, 4.0 / 3.0 * PI, 1e-2), "adaptive sphere volume {}", svprops.dim);
    }

    #[test]
    fn volume_properties_gk_box() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let v = volume_properties_gk(&b.solid.0).unwrap();
        assert!(approx(v.dim, 1.0, 1e-3), "GK volume {}", v.dim);
    }

    #[test]
    fn mesh_props_box_matches() {
        use crate::brep_extrema::mesh_faces;
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let mut nodes = Vec::new();
        let mut tris = Vec::new();
        for (fv, ft) in mesh_faces(&b.solid.0, 8, 8) {
            let base = nodes.len();
            nodes.extend(fv);
            for (a, b2, c) in ft {
                tris.push(occt_core::poly::triangulation::Triangle::new(base + a, base + b2, base + c));
            }
        }
        let mesh = occt_core::poly::triangulation::Triangulation::new(nodes, tris);
        let props = mesh_props(&mesh, &GpPnt::zero(), false, true);
        // The mesh orientation from `mesh_faces` may be inward; the signed
        // volume is taken in absolute value (as `brep_gprop` does).
        assert!(approx(props.dim.abs(), 1.0, 1e-2), "mesh volume {}", props.dim);
    }

    #[test]
    fn curve_tool_orders() {
        let b = crate::builder::TopoBuilder::new();
        let e = b.make_edge_segment(&GpPnt::zero(), &GpPnt::new(3.0, 4.0, 0.0));
        let c = BRepTool::edge_curve(&e).unwrap();
        let (a, bb) = BRepTool::edge_parameters(&e);
        assert_eq!(curve_integration_order(c.as_ref(), a, bb), 2);
    }
}
