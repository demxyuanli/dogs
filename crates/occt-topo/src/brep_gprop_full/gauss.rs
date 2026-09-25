use super::prelude::*;
use super::*;

// ---------------------------------------------------------------------------
// Constants (mirror BRepGProp_Gauss.cxx)
// ---------------------------------------------------------------------------

pub(super) const GPM: usize = 64; // math::GaussPointsMax()
pub(super) const EPS_PARAM: f64 = 1e-12;
pub(super) const EPS_DIM: f64 = 1e-30;
pub(super) const ERROR_ALGEBR_RATIO: f64 = 2.0 / 3.0;
pub(super) const SUBS_POWER: usize = 32;
pub(super) const SM: usize = SUBS_POWER * GPM + 1;
pub(super) const NGV: usize = 10; // number of value types in the GK integrator

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
pub(super) fn huygens_mat(m: f64, d: &GpVec) -> GpMat {
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
pub(super) struct Inertia {
    pub(super) mass: f64,
    pub(super) ix: f64,
    pub(super) iy: f64,
    pub(super) iz: f64,
    pub(super) ixx: f64,
    pub(super) iyy: f64,
    pub(super) izz: f64,
    pub(super) ixy: f64,
    pub(super) ixz: f64,
    pub(super) iyz: f64,
}

impl Inertia {
    pub(super) fn add(&mut self, o: &Inertia) {
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

    pub(super) fn mul(&mut self, s: f64) {
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
pub(super) fn convert_s(inertia: &Inertia) -> (f64, GpVec, GpMat) {
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
pub(super) fn convert_v(inertia: &Inertia, coeff: &[f64; 3]) -> (f64, GpVec, GpMat) {
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

pub(super) fn compute_v_inertia_elem(pt: &GpPnt, n: &GpVec, loc: &GpPnt, w: f64, coeff: &[f64; 3], by_point: bool, o: &mut Inertia) {
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

pub(super) fn compute_s_inertia_elem(pt: &GpPnt, n: &GpVec, loc: &GpPnt, w: f64, o: &mut Inertia) {
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
pub(super) enum SurfKind {
    Plane,
    Cylinder,
    Sphere,
    Torus,
    Cone,
    Other,
}

pub(super) fn classify_surface_kind(s: &dyn Surface) -> SurfKind {
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
pub(super) fn plane_frame(s: &dyn Surface) -> Option<(GpXyz, GpXyz, GpXyz)> {
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
pub(super) fn cylinder_frame(s: &dyn Surface) -> Option<(GpXyz, GpXyz, GpXyz, GpXyz, f64)> {
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
pub(super) fn cone_frame(s: &dyn Surface) -> Option<(GpXyz, GpXyz, GpXyz, GpXyz, f64)> {
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
pub(super) fn surface_d1(s: &dyn Surface, u: f64, v: f64) -> (GpPnt, GpVec, GpVec) {
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
pub(super) enum UVMap {
    Plane { o: GpXyz, xd: GpXyz, yd: GpXyz },
    Cylinder { o: GpXyz, xd: GpXyz, yd: GpXyz, z: GpXyz, r: f64 },
    Cone { o: GpXyz, xd: GpXyz, yd: GpXyz, z: GpXyz, alpha: f64 },
    Generic,
}

impl UVMap {
    pub(super) fn map_point(&self, s: &dyn Surface, p: &GpPnt) -> GpPnt2d {
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

    pub(super) fn map_deriv(&self, s: &dyn Surface, p: &GpPnt, d3: &GpVec, u: f64) -> GpVec2d {
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
pub(super) fn closest_params(s: &dyn Surface, p: &GpPnt) -> (f64, f64) {
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

pub(super) fn sane_bounds(a: f64, b: f64) -> (f64, f64) {
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
pub(super) enum ArcKind {
    Line,
    Circle,
    Other,
}

/// A boundary arc of a face — how it yields UV.
///
/// OCCT's `BRepGProp_Face` trims on the face's **pcurve**
/// (`BRep_Tool::CurveOnSurface`): only a pcurve can distinguish the two sides of
/// a seam (`u = 0` vs `u = 2π`) and a degenerate pole edge has no 3D curve at
/// all. The 3D-curve + UV-inverse form is kept as the fallback for faces whose
/// pcurves are absent (a port shortcut, documented as such).
pub(super) enum ArcGeom {
    Pcurve(Arc<dyn Curve2d>),
    Curve3d(Arc<dyn Curve>, UVMap),
}

/// A boundary arc of a face: its UV representation plus the edge's parameter
/// range, so `d12d` yields the pcurve point and derivative in the surface's UV.
pub(super) struct BoundaryArc {
    pub(super) geom: ArcGeom,
    pub(super) a: f64,
    pub(super) b: f64,
    pub(super) kind: ArcKind,
    /// The edge is REVERSED, so the arc is walked from `b` to `a`.
    pub(super) reversed: bool,
}

impl BoundaryArc {
    /// `BRepGProp_Face::Load(const TopoDS_Edge&)` (`BRepGProp_Face.cxx:173-179`)
    /// integrates the *reversed* pcurve over the reversed range when the edge is
    /// REVERSED. Mirroring the parameter inside `[a, b]` walks the identical UV
    /// points in the opposite order, which is what that replacement does
    /// geometrically (and it keeps `a < b`, as the Gauss loops require).
    fn param(&self, t: f64) -> f64 {
        if self.reversed {
            self.a + self.b - t
        } else {
            t
        }
    }

    pub(super) fn d12d(&self, s: &dyn Surface, t: f64) -> (GpPnt2d, GpVec2d) {
        let t = self.param(t);
        let (p, mut d) = match &self.geom {
            ArcGeom::Pcurve(pc) => pc.d1(t),
            ArcGeom::Curve3d(c, map) => {
                let p = c.d0(t);
                let d1 = c.d1(t).1;
                let puv = map.map_point(s, &p);
                let duv = map.map_deriv(s, &p, &d1, puv.x());
                (puv, duv)
            }
        };
        if self.reversed {
            d = GpVec2d::new(-d.x(), -d.y());
        }
        (p, d)
    }

    pub(super) fn value(&self, s: &dyn Surface, t: f64) -> GpPnt2d {
        let t = self.param(t);
        match &self.geom {
            ArcGeom::Pcurve(pc) => pc.d0(t),
            ArcGeom::Curve3d(c, map) => map.map_point(s, &c.d0(t)),
        }
    }
}

// ---------------------------------------------------------------------------
// FaceGauss — the BRepGProp_Face adaptor
// ---------------------------------------------------------------------------

/// Wraps a face and its surface with everything the integrators need:
/// UV bounds, surface normal, integration orders, and the boundary arcs.
pub(super) struct FaceGauss {
    pub(super) surface: Arc<dyn Surface>,
    pub(super) kind: SurfKind,
    pub(super) is_reversed: bool,
    pub(super) natural: bool,
    pub(super) u1: f64,
    pub(super) u2: f64,
    pub(super) v1: f64,
    pub(super) v2: f64,
    pub(super) arcs: Vec<BoundaryArc>,
    pub(super) has_repeated_edges: bool,
}
