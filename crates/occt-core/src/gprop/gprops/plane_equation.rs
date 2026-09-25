use super::prelude::*;

/// Plane equation `A·x + B·y + C·z + D = 0`.
///
/// Constructed from three points or from a point and a normal vector; supports
/// signed distance, orthogonal projection (nearest point) and conversion back
/// to a [`GpPln`].
///
/// *Note:* OCCT's `GProp_PEquation` is a principal-axis point-cloud fitter;
/// this port exposes the explicit plane-equation form requested for the
/// migration (three-point / point+normal construction).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PEquation {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub d: f64,
}

impl PEquation {
    /// Build the plane through `p1`, `p2`, `p3`.
    /// Returns `Err` when the points are collinear.
    pub fn from_points(p1: &GpPnt, p2: &GpPnt, p3: &GpPnt) -> Result<Self, String> {
        let u = p2.coord.subtracted(&p1.coord);
        let v = p3.coord.subtracted(&p1.coord);
        let n = u.crossed(&v);
        if n.modulus() <= CONFUSION {
            return Err("PEquation::from_points: collinear points".into());
        }
        Self::from_point_normal(p1, &GpVec::from_xyz(&n))
    }

    /// Build the plane with unit-consistent normal `n` passing through `p`.
    /// Returns `Err` when `n` is (near-)zero.
    pub fn from_point_normal(p: &GpPnt, n: &GpVec) -> Result<Self, String> {
        let nx = n.x();
        let ny = n.y();
        let nz = n.z();
        if !(nx.is_finite() && ny.is_finite() && nz.is_finite())
            || n.square_magnitude() <= CONFUSION * CONFUSION
        {
            return Err("PEquation::from_point_normal: zero normal".into());
        }
        let (a, b, c) = (n.x(), n.y(), n.z());
        let d = -(a * p.x() + b * p.y() + c * p.z());
        Ok(Self { a, b, c, d })
    }

    /// Build from a plane: its normal and location define the equation.
    pub fn from_plane(pln: &GpPln) -> Self {
        let axis = pln.axis();
        let n = axis.direction();
        let p = pln.location();
        Self::from_point_normal(&p, &GpVec::from_xyz(n.xyz())).expect("plane normal is non-zero")
    }

    /// The `(A, B, C, D)` coefficients of `A·x + B·y + C·z + D = 0`.
    pub fn coefficients(&self) -> (f64, f64, f64, f64) {
        (self.a, self.b, self.c, self.d)
    }

    /// The un-normalised normal vector `(A, B, C)`.
    pub fn normal(&self) -> GpVec {
        GpVec::new(self.a, self.b, self.c)
    }

    /// The unit normal vector `(A, B, C) / |(A, B, C)|`.
    pub fn normal_unit(&self) -> GpVec {
        let n = self.normal();
        let m = n.magnitude();
        if m > CONFUSION {
            n.divided(m)
        } else {
            n
        }
    }

    /// Signed distance from `p` to the plane (positive on the `(A,B,C)` side).
    /// `NaN` when the plane is degenerate.
    pub fn signed_distance(&self, p: &GpPnt) -> f64 {
        let denom = (self.a * self.a + self.b * self.b + self.c * self.c).sqrt();
        if denom < CONFUSION {
            return f64::NAN;
        }
        (self.a * p.x() + self.b * p.y() + self.c * p.z() + self.d) / denom
    }

    /// Unsigned distance from `p` to the plane.
    pub fn distance_to(&self, p: &GpPnt) -> f64 {
        self.signed_distance(p).abs()
    }

    /// Orthogonal projection of `p` onto the plane.
    pub fn project(&self, p: &GpPnt) -> GpPnt {
        let sd = self.signed_distance(p);
        let n = self.normal_unit();
        GpPnt::new(p.x() - sd * n.x(), p.y() - sd * n.y(), p.z() - sd * n.z())
    }

    /// Nearest point on the plane to `p` (same as [`Self::project`]).
    pub fn nearest_point(&self, p: &GpPnt) -> GpPnt {
        self.project(p)
    }

    /// Any point belonging to the plane (component with the largest
    /// coefficient set to zero, avoiding division by a near-zero value).
    pub fn point_on_plane(&self) -> GpPnt {
        let (a, b, c) = (self.a, self.b, self.c);
        if c.abs() >= a.abs() && c.abs() >= b.abs() {
            GpPnt::new(0.0, 0.0, -self.d / c)
        } else if b.abs() >= a.abs() {
            GpPnt::new(0.0, -self.d / b, 0.0)
        } else {
            GpPnt::new(-self.d / a, 0.0, 0.0)
        }
    }

    /// The [`GpPln`] equivalent of this plane equation.
    pub fn to_plane(&self) -> GpPln {
        let p = self.point_on_plane();
        let n = self.normal();
        let dir = GpDir::from_xyz(&n.coord).unwrap_or(GpDir::default_dir());
        GpPln::new(GpAx3::from_ax1(&GpAx1::new(p, dir)))
    }
}

// ---------------------------------------------------------------------------
// Element contribution helpers
// ---------------------------------------------------------------------------

/// Component `axis` (0=x, 1=y, 2=z) of a coordinate.
#[inline]
pub(super) fn comp(v: &GpXyz, axis: usize) -> f64 {
    match axis {
        0 => v.x,
        1 => v.y,
        _ => v.z,
    }
}

/// Build the inertia tensor `I = (tr S)·δ − S` from the second-moment matrix
/// `S` with `S[i][j] = ∫ r_i r_j dμ`.
#[inline]
pub(super) fn second_moments_to_inertia(s: &[[f64; 3]; 3]) -> GpMat {
    GpMat::new(
        s[1][1] + s[2][2],
        -s[0][1],
        -s[0][2],
        -s[1][0],
        s[0][0] + s[2][2],
        -s[1][2],
        -s[2][0],
        -s[2][1],
        s[0][0] + s[1][1],
    )
}

/// Inertia tensor of a unit point mass at `p`, about the origin:
/// `|p|²·δ − p·pᵀ`.
pub(super) fn point_inertia_origin(p: &GpPnt) -> GpMat {
    let (x, y, z) = (p.x(), p.y(), p.z());
    GpMat::new(
        y * y + z * z,
        -x * y,
        -x * z,
        -x * y,
        x * x + z * z,
        -y * z,
        -x * z,
        -y * z,
        x * x + y * y,
    )
}

/// Inertia tensor of a straight segment `p1 → p2` of length `len`, about the
/// origin. Integrates `r⊗r` along the segment in closed form:
/// `∫ r_i r_j ds = len·(a_i a_j + (a_i d_j + d_i a_j)/2 + d_i d_j/3)`
/// with `a = p1`, `d = p2 − p1`.
pub(super) fn segment_inertia_origin(p1: &GpPnt, p2: &GpPnt, len: f64) -> GpMat {
    let a = p1.coord;
    let d = p2.coord.subtracted(&a);
    let mut s = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let ai = comp(&a, i);
            let aj = comp(&a, j);
            let di = comp(&d, i);
            let dj = comp(&d, j);
            s[i][j] = len * (ai * aj + (ai * dj + di * aj) / 2.0 + di * dj / 3.0);
        }
    }
    second_moments_to_inertia(&s)
}

/// Area of triangle `abc`.
pub(super) fn triangle_area(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> f64 {
    let ab = b.coord.subtracted(&a.coord);
    let ac = c.coord.subtracted(&a.coord);
    0.5 * ab.crossed(&ac).modulus()
}

/// Inertia tensor of a flat triangle of area `area`, about the origin.
/// Barycentric integration over the triangle gives
/// `∫ r_i r_j dA = (A/6)·Σ_k r_k_i r_k_j + (A/12)·Σ_{k≠l} r_k_i r_l_j`.
pub(super) fn triangle_inertia_origin(a: &GpPnt, b: &GpPnt, c: &GpPnt, area: f64) -> GpMat {
    let v = [a.coord, b.coord, c.coord];
    let mut s = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let mut diag = 0.0;
            let mut cross = 0.0;
            for k in 0..3 {
                diag += comp(&v[k], i) * comp(&v[k], j);
                for l in 0..3 {
                    if k != l {
                        cross += comp(&v[k], i) * comp(&v[l], j);
                    }
                }
            }
            s[i][j] = area / 6.0 * diag + area / 12.0 * cross;
        }
    }
    second_moments_to_inertia(&s)
}

/// Unsigned volume of tetrahedron `abcd`.
pub(super) fn tetra_volume(a: &GpPnt, b: &GpPnt, c: &GpPnt, d: &GpPnt) -> f64 {
    let ab = b.coord.subtracted(&a.coord);
    let ac = c.coord.subtracted(&a.coord);
    let ad = d.coord.subtracted(&a.coord);
    ab.dot_cross(&ac, &ad).abs() / 6.0
}

/// Inertia tensor of a tetrahedron of volume `vol`, about the origin.
/// Barycentric integration gives
/// `∫ r_i r_j dV = (V/10)·Σ_k r_k_i r_k_j + (V/20)·Σ_{k≠l} r_k_i r_l_j`.
pub(super) fn tetra_inertia_origin(a: &GpPnt, b: &GpPnt, c: &GpPnt, d: &GpPnt, vol: f64) -> GpMat {
    let v = [a.coord, b.coord, c.coord, d.coord];
    let mut s = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let mut diag = 0.0;
            let mut cross = 0.0;
            for k in 0..4 {
                diag += comp(&v[k], i) * comp(&v[k], j);
                for l in 0..4 {
                    if k != l {
                        cross += comp(&v[k], i) * comp(&v[l], j);
                    }
                }
            }
            s[i][j] = vol / 10.0 * diag + vol / 20.0 * cross;
        }
    }
    second_moments_to_inertia(&s)
}

/// Split an axis-aligned box into 12 tetrahedra (each of the 6 faces, fan-split
/// into 2 triangles, connected to the box centre). Exactly tiles the box.
pub(super) fn box_tetrahedra(min: &GpPnt, max: &GpPnt) -> Vec<[GpPnt; 4]> {
    let center = GpPnt::new(
        (min.x() + max.x()) / 2.0,
        (min.y() + max.y()) / 2.0,
        (min.z() + max.z()) / 2.0,
    );
    let (x0, x1) = (min.x(), max.x());
    let (y0, y1) = (min.y(), max.y());
    let (z0, z1) = (min.z(), max.z());
    let v = [
        GpPnt::new(x0, y0, z0),
        GpPnt::new(x1, y0, z0),
        GpPnt::new(x1, y1, z0),
        GpPnt::new(x0, y1, z0),
        GpPnt::new(x0, y0, z1),
        GpPnt::new(x1, y0, z1),
        GpPnt::new(x1, y1, z1),
        GpPnt::new(x0, y1, z1),
    ];
    let faces: [[usize; 4]; 6] = [
        [0, 1, 2, 3], // z = z0
        [4, 5, 6, 7], // z = z1
        [0, 1, 5, 4], // y = y0
        [2, 3, 7, 6], // y = y1
        [0, 3, 7, 4], // x = x0
        [1, 2, 6, 5], // x = x1
    ];
    let mut out = Vec::with_capacity(12);
    for f in faces {
        out.push([v[f[0]], v[f[1]], v[f[2]], center]);
        out.push([v[f[0]], v[f[2]], v[f[3]], center]);
    }
    out
}

// ---------------------------------------------------------------------------
// Symmetric 3×3 Jacobi eigensolver
// ---------------------------------------------------------------------------

/// Eigen-decomposition of a symmetric 3×3 matrix by cyclic Jacobi rotations.
///
/// Returns `(eigenvalues descending, matching unit eigenvectors)`. Only valid
/// for symmetric input; the off-diagonal annihilation tolerance is absolute.
pub(super) fn jacobi_symmetric3(a: [[f64; 3]; 3]) -> Result<([f64; 3], [GpVec; 3]), String> {
    let mut m = a;
    let mut v = [[0.0f64; 3]; 3];
    for i in 0..3 {
        v[i][i] = 1.0;
    }
    let mut converged = false;
    for _ in 0..64 {
        let mut p = 0usize;
        let mut q = 1usize;
        let mut mx = m[0][1].abs();
        for i in 0..3 {
            for j in (i + 1)..3 {
                if m[i][j].abs() > mx {
                    mx = m[i][j].abs();
                    p = i;
                    q = j;
                }
            }
        }
        if mx < 1e-13 {
            converged = true;
            break;
        }
        let app = m[p][p];
        let aqq = m[q][q];
        let apq = m[p][q];
        let tau = (aqq - app) / (2.0 * apq);
        let t = tau.signum() / (tau.abs() + (1.0 + tau * tau).sqrt());
        let c = 1.0 / (1.0 + t * t).sqrt();
        let s = t * c;
        for k in 0..3 {
            if k == p || k == q {
                continue;
            }
            let akp = m[k][p];
            let akq = m[k][q];
            m[k][p] = c * akp - s * akq;
            m[p][k] = m[k][p];
            m[k][q] = s * akp + c * akq;
            m[q][k] = m[k][q];
        }
        m[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        m[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
        m[p][q] = 0.0;
        m[q][p] = 0.0;
        for k in 0..3 {
            let vkp = v[k][p];
            let vkq = v[k][q];
            v[k][p] = c * vkp - s * vkq;
            v[k][q] = s * vkp + c * vkq;
        }
    }
    if !converged {
        return Err("jacobi_symmetric3: did not converge".into());
    }
    let mut pairs: Vec<(f64, [f64; 3])> =
        (0..3).map(|i| (m[i][i], [v[0][i], v[1][i], v[2][i]])).collect();
    pairs.sort_by(|x, y| y.0.total_cmp(&x.0)); // descending
    let vals = [pairs[0].0, pairs[1].0, pairs[2].0];
    let mut vecs = [GpVec::new(0.0, 0.0, 0.0); 3];
    for (idx, (_, e)) in pairs.iter().enumerate() {
        let n = (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]).sqrt();
        let n = if n > 1e-30 { n } else { 1.0 };
        vecs[idx] = GpVec::new(e[0] / n, e[1] / n, e[2] / n);
    }
    Ok((vals, vecs))
}
