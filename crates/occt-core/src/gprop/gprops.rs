//! GProp — global properties framework. Source: `GProp/`.
//!
//! Ports the OCCT `GProp` package: the composite global-properties
//! accumulator [`GProps`] (`GProp_GProps`), the principal-properties
//! presentation [`PrincipalProps`] (`GProp_PrincipalProps`), the element
//! frameworks `PGProps` / `SelGProps` / `VelGProps` / `CelGProps`, the plane
//! equation [`PEquation`], the Huygens operator [`h_operator`] and the
//! [`ValueType`] enumeration.
//!
//! `GProps` composes the global properties (mass / length / area / volume,
//! centre of mass, quadratic inertia matrix) of a *compound geometric system*.
//! Elementary pieces are added either directly ([`GProps::add_point_mass`]) or
//! through the element frameworks, each of which accumulates one family of
//! geometric cells into a `GProps`:
//!
//! * [`PGProps`] — point set (`add_point`);
//! * [`SelGProps`] — surface (triangles);
//! * [`VelGProps`] — solid (tetrahedra / boxes / octahedra);
//! * [`CelGProps`] — curve (segments).
//!
//! All inertia is accumulated about the framework's reference point `loc`
//! (origin by default). [`GProps::matrix_of_inertia`] shifts it to the centre
//! of mass by the parallel-axis (Huygens) theorem, and [`h_operator`] provides
//! the parallel-axis term `m·(d²δ − d·dᵀ)` used everywhere.
//!
//! Each element contribution is accumulated exactly: a point mass uses
//! `r²δ − r·rᵀ`, a segment / triangle / tetrahedron uses the closed-form
//! integral of `r⊗r` over the cell (barycentric formulas), so a box
//! decomposed into tetrahedra reproduces the analytic box inertia tensor
//! exactly (up to floating-point rounding).

use crate::gp::{GpAx1, GpAx3, GpDir, GpMat, GpPln, GpPnt, GpVec, GpXyz};
use crate::precision::CONFUSION;
use std::ops::{Deref, DerefMut};

/// Kind of global-property value a `GProp` computation can deliver.
/// Source: `GProp_ValueType.hxx`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValueType {
    Mass,
    CenterMassX,
    CenterMassY,
    CenterMassZ,
    InertiaXx,
    InertiaYy,
    InertiaZz,
    InertiaXy,
    InertiaXz,
    InertiaYz,
    Unknown,
}

/// Huygens operator (parallel-axis term) of a system of mass `mass` whose
/// centre of mass is `g`, expressed at the point `q`:
///
/// ```text
/// H = mass · (|g − q|² · δ − (g − q)(g − q)ᵀ)
/// ```
///
/// Huygens' theorem shifts an inertia matrix about one point to another:
/// `I₀/Q = I₀/G + H`. Source: `GProp::HOperator`.
pub fn h_operator(g: &GpPnt, q: &GpPnt, mass: f64) -> GpMat {
    let d = g.coord.subtracted(&q.coord);
    let x = d.x;
    let y = d.y;
    let z = d.z;
    let ixx = y * y + z * z;
    let iyy = x * x + z * z;
    let izz = x * x + y * y;
    let ixy = -x * y;
    let iyz = -y * z;
    let ixz = -x * z;
    GpMat::new(ixx, ixy, ixz, ixy, iyy, iyz, ixz, iyz, izz).multiply_scalar(mass)
}

/// General mechanism to compute the global properties of a compound geometric
/// system in 3D by composition of the global properties of elementary pieces
/// (points, curves, surfaces, solids). Source: `GProp_GProps`.
///
/// Internal state:
/// * `loc` — reference point used for inertia accumulation;
/// * `g` — centre of mass *relative to* `loc` (absolute COM = `loc + g`);
/// * `dim` — total mass / length / area / volume;
/// * `inertia` — quadratic-moment (inertia) matrix about `loc`.
#[derive(Debug, Clone)]
pub struct GProps {
    g: GpPnt,
    loc: GpPnt,
    dim: f64,
    inertia: GpMat,
}

impl Default for GProps {
    fn default() -> Self {
        Self::new()
    }
}

impl GProps {
    /// Empty system whose reference point is the absolute origin.
    pub fn new() -> Self {
        Self {
            g: GpPnt::zero(),
            loc: GpPnt::zero(),
            dim: 0.0,
            inertia: GpMat::zero(),
        }
    }

    /// Empty system whose inertia is accumulated about `system_location`.
    /// For greater accuracy choose a point close to the system (e.g. near its
    /// expected centre of mass).
    pub fn new_at(system_location: GpPnt) -> Self {
        Self {
            g: GpPnt::zero(),
            loc: system_location,
            dim: 0.0,
            inertia: GpMat::zero(),
        }
    }

    /// Reference point used for inertia accumulation.
    pub fn reference_point(&self) -> GpPnt {
        self.loc
    }

    /// Override the reference point. Should be called before any `add*` call,
    /// otherwise the accumulated inertia and centre-of-mass are left relative
    /// to the previous reference point (matching the OCCT usage contract).
    pub fn set_location(&mut self, loc: GpPnt) {
        self.loc = loc;
    }

    /// Total mass (or length / area / volume when the pieces carry unit
    /// density) of the current system.
    pub fn mass(&self) -> f64 {
        self.dim
    }

    /// Centre of mass of the current system in the absolute frame.
    pub fn centre_of_mass(&self) -> GpPnt {
        GpPnt::from_xyz(&self.loc.coord.added(&self.g.coord))
    }

    /// The symmetric 3×3 matrix of inertia in the central coordinate system
    /// `(G, Gx, Gy, Gz)` whose origin `G` is the centre of mass. To compute it
    /// at another location use [`h_operator`] (Huygens' theorem).
    pub fn matrix_of_inertia(&self) -> GpMat {
        self.inertia.subtract(&h_operator(&self.g, &GpPnt::zero(), self.dim))
    }

    /// Static moments of inertia about the three axes of the absolute
    /// Cartesian frame: `(m·x_G, m·y_G, m·z_G)`.
    pub fn static_moments(&self) -> (f64, f64, f64) {
        let c = self.centre_of_mass();
        (c.x() * self.dim, c.y() * self.dim, c.z() * self.dim)
    }

    /// Moment of inertia of the system about the axis `a`.
    ///
    /// The inertia matrix is evaluated about `a.location()` (moving from the
    /// stored reference point or the centre of mass by Huygens' theorem) and
    /// contracted with the axis direction: `I = uᵀ·I·u`.
    pub fn moment_of_inertia(&self, a: &GpAx1) -> f64 {
        let dir = *a.direction().xyz();
        if self.loc.distance(a.location()) <= CONFUSION {
            let iv = self.inertia.multiplied(&dir);
            dir.dot(&iv)
        } else {
            let com = self.centre_of_mass();
            let axis_inertia =
                self.matrix_of_inertia().add(&h_operator(&com, a.location(), self.dim));
            let iv = axis_inertia.multiplied(&dir);
            dir.dot(&iv)
        }
    }

    /// Radius of gyration of the system about the axis `a`:
    /// `sqrt(moment_of_inertia(a) / mass)`. Returns 0 for a massless system.
    pub fn radius_of_gyration(&self, a: &GpAx1) -> f64 {
        let m = self.moment_of_inertia(a);
        if self.dim.abs() < 1e-30 {
            0.0
        } else {
            (m / self.dim).abs().sqrt()
        }
    }

    /// Principal properties of inertia of the system: the eigenvalues
    /// (principal moments) and eigenvectors (principal axes) of the central
    /// inertia matrix, packaged in a [`PrincipalProps`].
    pub fn principal_properties(&self) -> Result<PrincipalProps, String> {
        let mat = self.matrix_of_inertia();
        let a = [
            [mat.m[0][0], mat.m[0][1], mat.m[0][2]],
            [mat.m[1][0], mat.m[1][1], mat.m[1][2]],
            [mat.m[2][0], mat.m[2][1], mat.m[2][2]],
        ];
        let (vals, vecs) = jacobi_symmetric3(a)?;
        let (i1, i2, i3) = (vals[0], vals[1], vals[2]);
        let (v1, v2, v3) = (vecs[0], vecs[1], vecs[2]);
        let (r1, r2, r3) = if self.dim.abs() < 1e-30 {
            (0.0, 0.0, 0.0)
        } else {
            (
                (i1 / self.dim).abs().sqrt(),
                (i2 / self.dim).abs().sqrt(),
                (i3 / self.dim).abs().sqrt(),
            )
        };
        Ok(PrincipalProps::new_full(i1, i2, i3, r1, r2, r3, v1, v2, v3, self.centre_of_mass()))
    }

    /// Compose the global properties of `item` into the current system.
    /// `density` (default 1.0 in OCCT) scales the mass and inertia of `item`.
    ///
    /// The reference point of `item` may differ from `self.loc`; Huygens'
    /// theorem is applied automatically to transfer the inertia values.
    /// Returns `Err` when `density <= CONFUSION` (OCCT `Standard_DomainError`).
    pub fn add(&mut self, item: &GProps, density: f64) -> Result<(), String> {
        if !density.is_finite() || density <= CONFUSION {
            return Err("density <= resolution".into());
        }
        if self.loc.distance(&item.loc) <= CONFUSION {
            // Same reference point: weight the centres of mass directly.
            let mut gxyz = item.g.coord.multiplied(item.dim * density);
            gxyz = gxyz.added(&self.g.coord.multiplied(self.dim));
            self.dim += item.dim * density;
            self.g = if self.dim.abs() >= 1e-20 {
                GpPnt::from_xyz(&gxyz.divided(self.dim))
            } else {
                GpPnt::zero()
            };
            self.inertia = self.inertia.add(&item.inertia.multiply_scalar(density));
        } else {
            // Different reference points: shift item's inertia to self.loc.
            let itemloc = self.loc.coord.subtracted(&item.loc.coord);
            let itemg = item.loc.coord.added(&item.g.coord);
            let mut gxyz = item.g.coord.subtracted(&itemloc);
            gxyz = gxyz.multiplied(item.dim * density);
            gxyz = gxyz.added(&self.g.coord.multiplied(self.dim));
            self.dim += item.dim * density;
            self.g = if self.dim.abs() >= 1e-20 {
                GpPnt::from_xyz(&gxyz.divided(self.dim))
            } else {
                GpPnt::zero()
            };
            // Item inertia about item's own centre of mass, then about self.loc.
            let mut item_inertia = item.inertia;
            if item.g.coord.modulus() > CONFUSION {
                let h = h_operator(&GpPnt::from_xyz(&itemg), &item.loc, item.dim);
                item_inertia = item_inertia.subtract(&h);
            }
            let h = h_operator(&GpPnt::from_xyz(&itemg), &self.loc, item.dim);
            item_inertia = item_inertia.add(&h);
            self.inertia = self.inertia.add(&item_inertia.multiply_scalar(density));
        }
        Ok(())
    }

    /// Add a single point mass of magnitude `density` at `p`.
    pub fn add_point_mass(&mut self, p: &GpPnt, density: f64) -> Result<(), String> {
        self.add_element(*p, 1.0, point_inertia_origin(p), density)
    }

    /// Add a point set of unit-mass points (shortcut for building a
    /// [`PGProps`] and composing it with unit density).
    pub fn add_points(&mut self, points: &[GpPnt]) -> Result<(), String> {
        for p in points {
            self.add_point_mass(p, 1.0)?;
        }
        Ok(())
    }

    /// Internal helper: compose one element whose inertia is known about the
    /// origin and whose centre is `centroid`.
    fn add_element(
        &mut self,
        centroid: GpPnt,
        dim: f64,
        inertia_origin: GpMat,
        density: f64,
    ) -> Result<(), String> {
        let el = GProps {
            g: centroid,
            loc: GpPnt::zero(),
            dim,
            inertia: inertia_origin,
        };
        self.add(&el, density)
    }
}

/// Presentation of the principal properties of inertia of a system.
/// Source: `GProp_PrincipalProps`.
///
/// The principal axes pass through the centre of mass and are parallel to the
/// eigenvectors of the central inertia matrix; the principal moments are the
/// associated eigenvalues (reported in descending order), and the principal
/// radii of gyration are `sqrt(|I|/m)`.
#[derive(Debug, Clone)]
pub struct PrincipalProps {
    i1: f64,
    i2: f64,
    i3: f64,
    r1: f64,
    r2: f64,
    r3: f64,
    v1: GpVec,
    v2: GpVec,
    v3: GpVec,
    g: GpPnt,
}

impl PrincipalProps {
    /// Creates an undefined `PrincipalProps` (all fields sentinel `f64::MAX`).
    pub fn new() -> Self {
        let mx = f64::MAX;
        Self {
            i1: mx,
            i2: mx,
            i3: mx,
            r1: mx,
            r2: mx,
            r3: mx,
            v1: GpVec::new(1.0, 0.0, 0.0),
            v2: GpVec::new(0.0, 1.0, 0.0),
            v3: GpVec::new(0.0, 0.0, 1.0),
            g: GpPnt::new(mx, mx, mx),
        }
    }

    /// Full constructor with explicit moments, radii, axes and centre of mass.
    pub fn new_full(
        i1: f64,
        i2: f64,
        i3: f64,
        r1: f64,
        r2: f64,
        r3: f64,
        v1: GpVec,
        v2: GpVec,
        v3: GpVec,
        g: GpPnt,
    ) -> Self {
        Self { i1, i2, i3, r1, r2, r3, v1, v2, v3, g }
    }

    /// Principal moments of inertia `(I1, I2, I3)` in descending order.
    pub fn moments(&self) -> (f64, f64, f64) {
        (self.i1, self.i2, self.i3)
    }

    /// First (largest-moment) principal axis of inertia.
    pub fn first_axis_of_inertia(&self) -> GpVec {
        self.v1
    }

    /// Second principal axis of inertia.
    pub fn second_axis_of_inertia(&self) -> GpVec {
        self.v2
    }

    /// Third (smallest-moment) principal axis of inertia.
    pub fn third_axis_of_inertia(&self) -> GpVec {
        self.v3
    }

    /// Principal radii of gyration `(R1, R2, R3)` about the principal axes.
    pub fn radius_of_gyration(&self) -> (f64, f64, f64) {
        (self.r1, self.r2, self.r3)
    }

    /// Centre of mass of the system the principal properties were computed for.
    pub fn centre_of_mass(&self) -> GpPnt {
        self.g
    }

    /// True if the system has an axis of symmetry (two principal moments
    /// equal within a relative tolerance of `1e-10`).
    pub fn has_symmetry_axis(&self) -> bool {
        self.has_symmetry_axis_tol(1e-10)
    }

    /// True if the system has an axis of symmetry, comparing moments with a
    /// caller-provided relative tolerance (plus machine epsilon when `tol` is 0).
    pub fn has_symmetry_axis_tol(&self, tol: f64) -> bool {
        let eps1 = self.i1.abs() * tol + f64::EPSILON * self.i1.abs();
        let eps2 = self.i2.abs() * tol + f64::EPSILON * self.i2.abs();
        (self.i1 - self.i2).abs() <= eps1
            || (self.i1 - self.i3).abs() <= eps1
            || (self.i2 - self.i3).abs() <= eps2
    }

    /// True if the system has a point of symmetry (all three principal moments
    /// equal within a relative tolerance of `1e-10`).
    pub fn has_symmetry_point(&self) -> bool {
        self.has_symmetry_point_tol(1e-10)
    }

    /// True if the system has a point of symmetry, comparing moments with a
    /// caller-provided relative tolerance (plus machine epsilon when `tol` is 0).
    pub fn has_symmetry_point_tol(&self, tol: f64) -> bool {
        let eps1 = self.i1.abs() * tol + f64::EPSILON * self.i1.abs();
        (self.i1 - self.i2).abs() <= eps1 && (self.i1 - self.i3).abs() <= eps1
    }
}

impl Default for PrincipalProps {
    fn default() -> Self {
        Self::new()
    }
}

/// Point-set framework: accumulates the global properties of a weighted set of
/// 3D points. Source: `GProp_PGProps`.
///
/// Inertia is accumulated at the absolute origin. The reference point of the
/// underlying `GProps` is the origin, so a `PGProps` can be composed into a
/// larger system via [`GProps::add`].
#[derive(Debug, Clone)]
pub struct PGProps {
    props: GProps,
}

impl PGProps {
    /// Empty point set at the origin with zero mass.
    pub fn new() -> Self {
        Self { props: GProps::new() }
    }

    /// Point set from a slice of unit-mass points.
    pub fn from_points(points: &[GpPnt]) -> Self {
        let mut p = Self::new();
        for pt in points {
            p.add_point(pt).expect("unit density is valid");
        }
        p
    }

    /// Add a point of unit mass.
    pub fn add_point(&mut self, p: &GpPnt) -> Result<(), String> {
        self.add_point_density(p, 1.0)
    }

    /// Add a point with a given mass.
    /// Returns `Err` when `density <= CONFUSION` (OCCT `Standard_DomainError`).
    pub fn add_point_density(&mut self, p: &GpPnt, density: f64) -> Result<(), String> {
        self.add_element(*p, 1.0, point_inertia_origin(p), density)
    }

    /// Barycentre of a set of unit-mass points.
    pub fn barycentre(points: &[GpPnt]) -> GpPnt {
        let n = points.len();
        if n == 0 {
            return GpPnt::zero();
        }
        let mut sum = GpXyz::zero();
        for p in points {
            sum = sum.added(&p.coord);
        }
        GpPnt::from_xyz(&sum.divided(n as f64))
    }

    /// Weighted barycentre and total mass of a point set.
    /// Returns `Err` on length mismatch or a non-positive density.
    pub fn weighted_barycentre(
        points: &[GpPnt],
        densities: &[f64],
    ) -> Result<(f64, GpPnt), String> {
        if points.len() != densities.len() {
            return Err("PGProps::weighted_barycentre: length mismatch".into());
        }
        let mut mass = 0.0;
        let mut sum = GpXyz::zero();
        for (p, d) in points.iter().zip(densities) {
            if !d.is_finite() || *d <= CONFUSION {
                return Err("density <= resolution".into());
            }
            mass += *d;
            sum = sum.added(&p.coord.multiplied(*d));
        }
        if mass.abs() < 1e-30 {
            return Err("PGProps::weighted_barycentre: zero total mass".into());
        }
        Ok((mass, GpPnt::from_xyz(&sum.divided(mass))))
    }
}

impl Default for PGProps {
    fn default() -> Self {
        Self::new()
    }
}

impl Deref for PGProps {
    type Target = GProps;
    fn deref(&self) -> &GProps {
        &self.props
    }
}

impl DerefMut for PGProps {
    fn deref_mut(&mut self) -> &mut GProps {
        &mut self.props
    }
}

/// Surface framework: accumulates the global properties of a bounded surface
/// approximated by planar triangular facets. Source: `GProp_SelGProps`.
#[derive(Debug, Clone)]
pub struct SelGProps {
    props: GProps,
}

impl SelGProps {
    /// Empty surface system.
    pub fn new() -> Self {
        Self { props: GProps::new() }
    }

    /// Set the reference point used for inertia accumulation (call before
    /// adding triangles, matching the OCCT `SetLocation` contract).
    pub fn set_location(&mut self, loc: GpPnt) {
        self.props.loc = loc;
    }

    /// Add a triangular facet with the given `density` (the mass contribution
    /// is the facet area times `density`).
    pub fn add_triangle(
        &mut self,
        a: &GpPnt,
        b: &GpPnt,
        c: &GpPnt,
        density: f64,
    ) -> Result<(), String> {
        let area = triangle_area(a, b, c);
        if area <= 1e-30 {
            return Ok(());
        }
        let centroid = GpPnt::new(
            (a.x() + b.x() + c.x()) / 3.0,
            (a.y() + b.y() + c.y()) / 3.0,
            (a.z() + b.z() + c.z()) / 3.0,
        );
        self.add_element(centroid, area, triangle_inertia_origin(a, b, c, area), density)
    }
}

impl Default for SelGProps {
    fn default() -> Self {
        Self::new()
    }
}

impl Deref for SelGProps {
    type Target = GProps;
    fn deref(&self) -> &GProps {
        &self.props
    }
}

impl DerefMut for SelGProps {
    fn deref_mut(&mut self) -> &mut GProps {
        &mut self.props
    }
}

/// Solid framework: accumulates the global properties of a closed 3D region
/// decomposed into tetrahedra. Source: `GProp_VelGProps`.
#[derive(Debug, Clone)]
pub struct VelGProps {
    props: GProps,
}

impl VelGProps {
    /// Empty solid system.
    pub fn new() -> Self {
        Self { props: GProps::new() }
    }

    /// Set the reference point used for inertia accumulation (call before
    /// adding cells, matching the OCCT `SetLocation` contract).
    pub fn set_location(&mut self, loc: GpPnt) {
        self.props.loc = loc;
    }

    /// Add a tetrahedron cell with the given `density` (the mass contribution
    /// is the unsigned cell volume times `density`).
    pub fn add_tetrahedron(
        &mut self,
        a: &GpPnt,
        b: &GpPnt,
        c: &GpPnt,
        d: &GpPnt,
        density: f64,
    ) -> Result<(), String> {
        let vol = tetra_volume(a, b, c, d);
        if vol <= 1e-30 {
            return Ok(());
        }
        let centroid = GpPnt::new(
            (a.x() + b.x() + c.x() + d.x()) / 4.0,
            (a.y() + b.y() + c.y() + d.y()) / 4.0,
            (a.z() + b.z() + c.z() + d.z()) / 4.0,
        );
        self.add_element(centroid, vol, tetra_inertia_origin(a, b, c, d, vol), density)
    }

    /// Add an axis-aligned box `[min, max]` by decomposing it into 12
    /// tetrahedra (6 face pyramids, each split into 2 cells).
    pub fn add_box(&mut self, min: &GpPnt, max: &GpPnt, density: f64) -> Result<(), String> {
        for t in box_tetrahedra(min, max) {
            self.add_tetrahedron(&t[0], &t[1], &t[2], &t[3], density)?;
        }
        Ok(())
    }

    /// Add an octahedron with apexes `top` / `bottom` and an equatorial
    /// quadrilateral `equator` (vertices in cyclic order), split into 4
    /// tetrahedra (2 per pyramid half).
    pub fn add_octahedron(
        &mut self,
        top: &GpPnt,
        bottom: &GpPnt,
        equator: &[GpPnt; 4],
        density: f64,
    ) -> Result<(), String> {
        self.add_tetrahedron(top, &equator[0], &equator[1], &equator[2], density)?;
        self.add_tetrahedron(top, &equator[0], &equator[2], &equator[3], density)?;
        self.add_tetrahedron(bottom, &equator[0], &equator[1], &equator[2], density)?;
        self.add_tetrahedron(bottom, &equator[0], &equator[2], &equator[3], density)?;
        Ok(())
    }
}

impl Default for VelGProps {
    fn default() -> Self {
        Self::new()
    }
}

impl Deref for VelGProps {
    type Target = GProps;
    fn deref(&self) -> &GProps {
        &self.props
    }
}

impl DerefMut for VelGProps {
    fn deref_mut(&mut self) -> &mut GProps {
        &mut self.props
    }
}

/// Curve framework: accumulates the global properties of a bounded curve
/// approximated by straight segments. Source: `GProp_CelGProps`.
#[derive(Debug, Clone)]
pub struct CelGProps {
    props: GProps,
}

impl CelGProps {
    /// Empty curve system.
    pub fn new() -> Self {
        Self { props: GProps::new() }
    }

    /// Set the reference point used for inertia accumulation (call before
    /// adding segments, matching the OCCT `SetLocation` contract).
    pub fn set_location(&mut self, loc: GpPnt) {
        self.props.loc = loc;
    }

    /// Add a straight segment `p1 → p2` with the given `density` (the mass
    /// contribution is the segment length times `density`).
    pub fn add_segment(
        &mut self,
        p1: &GpPnt,
        p2: &GpPnt,
        density: f64,
    ) -> Result<(), String> {
        let len = p1.distance(p2);
        if len <= 1e-30 {
            return Ok(());
        }
        let centroid = GpPnt::new(
            (p1.x() + p2.x()) / 2.0,
            (p1.y() + p2.y()) / 2.0,
            (p1.z() + p2.z()) / 2.0,
        );
        self.add_element(centroid, len, segment_inertia_origin(p1, p2, len), density)
    }
}

impl Default for CelGProps {
    fn default() -> Self {
        Self::new()
    }
}

impl Deref for CelGProps {
    type Target = GProps;
    fn deref(&self) -> &GProps {
        &self.props
    }
}

impl DerefMut for CelGProps {
    fn deref_mut(&mut self) -> &mut GProps {
        &mut self.props
    }
}

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
fn comp(v: &GpXyz, axis: usize) -> f64 {
    match axis {
        0 => v.x,
        1 => v.y,
        _ => v.z,
    }
}

/// Build the inertia tensor `I = (tr S)·δ − S` from the second-moment matrix
/// `S` with `S[i][j] = ∫ r_i r_j dμ`.
#[inline]
fn second_moments_to_inertia(s: &[[f64; 3]; 3]) -> GpMat {
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
fn point_inertia_origin(p: &GpPnt) -> GpMat {
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
fn segment_inertia_origin(p1: &GpPnt, p2: &GpPnt, len: f64) -> GpMat {
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
fn triangle_area(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> f64 {
    let ab = b.coord.subtracted(&a.coord);
    let ac = c.coord.subtracted(&a.coord);
    0.5 * ab.crossed(&ac).modulus()
}

/// Inertia tensor of a flat triangle of area `area`, about the origin.
/// Barycentric integration over the triangle gives
/// `∫ r_i r_j dA = (A/6)·Σ_k r_k_i r_k_j + (A/12)·Σ_{k≠l} r_k_i r_l_j`.
fn triangle_inertia_origin(a: &GpPnt, b: &GpPnt, c: &GpPnt, area: f64) -> GpMat {
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
fn tetra_volume(a: &GpPnt, b: &GpPnt, c: &GpPnt, d: &GpPnt) -> f64 {
    let ab = b.coord.subtracted(&a.coord);
    let ac = c.coord.subtracted(&a.coord);
    let ad = d.coord.subtracted(&a.coord);
    ab.dot_cross(&ac, &ad).abs() / 6.0
}

/// Inertia tensor of a tetrahedron of volume `vol`, about the origin.
/// Barycentric integration gives
/// `∫ r_i r_j dV = (V/10)·Σ_k r_k_i r_k_j + (V/20)·Σ_{k≠l} r_k_i r_l_j`.
fn tetra_inertia_origin(a: &GpPnt, b: &GpPnt, c: &GpPnt, d: &GpPnt, vol: f64) -> GpMat {
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
fn box_tetrahedra(min: &GpPnt, max: &GpPnt) -> Vec<[GpPnt; 4]> {
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
fn jacobi_symmetric3(a: [[f64; 3]; 3]) -> Result<([f64; 3], [GpVec; 3]), String> {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gp::dir::DirAxis;
    use std::f64::consts::PI;

    fn assert_approx(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() < tol, "{a} != {b} within {tol}");
    }

    fn z_axis() -> GpAx1 {
        GpAx1::new(GpPnt::zero(), GpDir::from_axis(DirAxis::Z))
    }

    #[test]
    fn h_operator_point_above_origin() {
        let g = GpPnt::new(0.0, 0.0, 2.0);
        let m = h_operator(&g, &GpPnt::zero(), 1.0);
        assert_approx(m.m[0][0], 4.0, 1e-12); // y²+z²
        assert_approx(m.m[1][1], 4.0, 1e-12); // x²+z²
        assert_approx(m.m[2][2], 0.0, 1e-12); // x²+y²
        assert_approx(m.m[0][1], 0.0, 1e-12);
        assert_approx(m.m[0][2], 0.0, 1e-12);
    }

    #[test]
    fn empty_system_is_zero() {
        let s = GProps::new();
        assert_approx(s.mass(), 0.0, 1e-30);
        assert_approx(s.centre_of_mass().x(), 0.0, 1e-30);
        assert_approx(s.centre_of_mass().y(), 0.0, 1e-30);
        assert_approx(s.centre_of_mass().z(), 0.0, 1e-30);
        for r in 0..3 {
            for c in 0..3 {
                assert_approx(s.matrix_of_inertia().m[r][c], 0.0, 1e-30);
            }
        }
    }

    #[test]
    fn density_must_be_positive() {
        let mut s = GProps::new();
        let p = GProps::new();
        assert!(s.add(&p, 0.0).is_err());
        assert!(s.add(&p, -1.0).is_err());
        assert!(s.add(&p, 1e-8).is_err());
        assert!(s.add(&p, f64::NAN).is_err());
        assert!(s.add(&p, 1.0).is_ok());
    }

    #[test]
    fn point_mass_props() {
        let mut s = GProps::new();
        s.add_point_mass(&GpPnt::new(1.0, 0.0, 0.0), 2.0).unwrap();
        assert_approx(s.mass(), 2.0, 1e-12);
        assert_approx(s.centre_of_mass().x(), 1.0, 1e-12);
        // Moment about the z axis through the origin: m·1² = 2.
        assert_approx(s.moment_of_inertia(&z_axis()), 2.0, 1e-12);
        assert_approx(s.radius_of_gyration(&z_axis()), 1.0, 1e-12);
        // Static moments: COM·mass = (2, 0, 0).
        let (ix, iy, iz) = s.static_moments();
        assert_approx(ix, 2.0, 1e-12);
        assert_approx(iy, 0.0, 1e-12);
        assert_approx(iz, 0.0, 1e-12);
    }

    #[test]
    fn point_cloud_cube_corners() {
        // 8 unit masses at the corners of [0,1]³: mass 8, COM at the centre,
        // each principal moment about the COM = Σ(y²+z²) = 8·(0.25+0.25) = 4.
        let pts = [
            GpPnt::new(0., 0., 0.),
            GpPnt::new(1., 0., 0.),
            GpPnt::new(1., 1., 0.),
            GpPnt::new(0., 1., 0.),
            GpPnt::new(0., 0., 1.),
            GpPnt::new(1., 0., 1.),
            GpPnt::new(1., 1., 1.),
            GpPnt::new(0., 1., 1.),
        ];
        let mut pg = PGProps::new();
        for p in &pts {
            pg.add_point(p).unwrap();
        }
        assert_approx(pg.mass(), 8.0, 1e-12);
        let c = pg.centre_of_mass();
        assert_approx(c.x(), 0.5, 1e-12);
        assert_approx(c.y(), 0.5, 1e-12);
        assert_approx(c.z(), 0.5, 1e-12);
        let pp = pg.principal_properties().unwrap();
        let (i1, i2, i3) = pp.moments();
        assert_approx(i1, 4.0, 1e-9);
        assert_approx(i2, 4.0, 1e-9);
        assert_approx(i3, 4.0, 1e-9);
        assert!(pp.has_symmetry_point());
        assert!(pp.has_symmetry_axis());
        let (r1, _, _) = pp.radius_of_gyration();
        assert_approx(r1, 0.5f64.sqrt(), 1e-9);
    }

    #[test]
    fn solid_unit_cube() {
        let mut v = VelGProps::new();
        v.add_box(&GpPnt::new(0., 0., 0.), &GpPnt::new(1., 1., 1.), 1.0).unwrap();
        assert_approx(v.mass(), 1.0, 1e-12);
        let c = v.centre_of_mass();
        assert_approx(c.x(), 0.5, 1e-9);
        assert_approx(c.y(), 0.5, 1e-9);
        assert_approx(c.z(), 0.5, 1e-9);
        // Solid cube: I = m/12·(1+1) = 1/6 along each axis.
        let pp = v.principal_properties().unwrap();
        let (i1, i2, i3) = pp.moments();
        assert_approx(i1, 1.0 / 6.0, 1e-9);
        assert_approx(i2, 1.0 / 6.0, 1e-9);
        assert_approx(i3, 1.0 / 6.0, 1e-9);
        assert!(pp.has_symmetry_point());
        let (r1, _, _) = pp.radius_of_gyration();
        assert_approx(r1, (1.0f64 / 6.0).sqrt(), 1e-9);
    }

    #[test]
    fn solid_box_2x3x4() {
        // Density 1, box 2×3×4: volume 24, COM at the centre.
        let mut v = VelGProps::new();
        v.add_box(&GpPnt::new(0., 0., 0.), &GpPnt::new(2., 3., 4.), 1.0).unwrap();
        assert_approx(v.mass(), 24.0, 1e-9);
        let c = v.centre_of_mass();
        assert_approx(c.x(), 1.0, 1e-9);
        assert_approx(c.y(), 1.5, 1e-9);
        assert_approx(c.z(), 2.0, 1e-9);
        // Analytic box inertia about COM:
        //   Ixx = m/12·(dy²+dz²) = 24/12·(9+16) = 50
        //   Iyy = m/12·(dx²+dz²) = 24/12·(4+16) = 40
        //   Izz = m/12·(dx²+dy²) = 24/12·(4+9)  = 26
        let m = v.matrix_of_inertia();
        assert_approx(m.m[0][0], 50.0, 1e-9);
        assert_approx(m.m[1][1], 40.0, 1e-9);
        assert_approx(m.m[2][2], 26.0, 1e-9);
        assert_approx(m.m[0][1], 0.0, 1e-9);
        assert_approx(m.m[0][2], 0.0, 1e-9);
        assert_approx(m.m[1][2], 0.0, 1e-9);
        let pp = v.principal_properties().unwrap();
        let (i1, i2, i3) = pp.moments();
        assert_approx(i1, 50.0, 1e-9);
        assert_approx(i2, 40.0, 1e-9);
        assert_approx(i3, 26.0, 1e-9);
        assert!(!pp.has_symmetry_point());
        assert!(!pp.has_symmetry_axis());
        // Principal axes of the axis-aligned box are the coordinate axes.
        let a = pp.first_axis_of_inertia();
        assert!(a.x().abs() > 0.9);
    }

    #[test]
    fn segment_properties() {
        let mut c = CelGProps::new();
        c.add_segment(&GpPnt::new(0., 0., 0.), &GpPnt::new(2., 0., 0.), 1.0).unwrap();
        assert_approx(c.mass(), 2.0, 1e-12);
        assert_approx(c.centre_of_mass().x(), 1.0, 1e-12);
        // About COM: Ixx = 0, Iyy = Izz = m·L²/12 = 2·4/12 = 2/3.
        let m = c.matrix_of_inertia();
        assert_approx(m.m[0][0], 0.0, 1e-12);
        assert_approx(m.m[1][1], 2.0 / 3.0, 1e-12);
        assert_approx(m.m[2][2], 2.0 / 3.0, 1e-12);
    }

    #[test]
    fn surface_disk() {
        // Triangulate a radius-1 disk from its centre: area → π, Izz → π/2.
        let r = 1.0f64;
        let n = 256usize;
        let mut s = SelGProps::new();
        let center = GpPnt::new(0., 0., 0.);
        for i in 0..n {
            let t1 = 2.0 * PI * i as f64 / n as f64;
            let t2 = 2.0 * PI * (i + 1) as f64 / n as f64;
            let p1 = GpPnt::new(r * t1.cos(), r * t1.sin(), 0.);
            let p2 = GpPnt::new(r * t2.cos(), r * t2.sin(), 0.);
            s.add_triangle(&center, &p1, &p2, 1.0).unwrap();
        }
        assert_approx(s.mass(), PI * r * r, 0.01);
        let i = s.matrix_of_inertia();
        assert_approx(i.m[2][2], PI * r * r * r * r / 2.0, 0.02);
        assert_approx(s.centre_of_mass().x(), 0.0, 1e-9);
        assert_approx(s.centre_of_mass().y(), 0.0, 1e-9);
    }

    #[test]
    fn surface_sphere_shell() {
        // UV-grid triangulation of a unit sphere surface: area → 4π,
        // Izz → ∫(x²+y²)dA = (8/3)π (isotropic shell, about the centre).
        let r = 1.0f64;
        let (nu, nv) = (48usize, 32usize);
        let mut s = SelGProps::new();
        let pt = |u: f64, v: f64| GpPnt::new(r * u.cos() * v.sin(), r * u.sin() * v.sin(), r * v.cos());
        for i in 0..nu {
            let u0 = 2.0 * PI * i as f64 / nu as f64;
            let u1 = 2.0 * PI * (i + 1) as f64 / nu as f64;
            for j in 0..nv {
                let v0 = PI * j as f64 / nv as f64;
                let v1 = PI * (j + 1) as f64 / nv as f64;
                s.add_triangle(&pt(u0, v0), &pt(u1, v0), &pt(u0, v1), 1.0).unwrap();
                s.add_triangle(&pt(u1, v0), &pt(u1, v1), &pt(u0, v1), 1.0).unwrap();
            }
        }
        assert_approx(s.mass(), 4.0 * PI * r * r, 0.5);
        let i = s.matrix_of_inertia();
        assert_approx(i.m[2][2], (8.0 / 3.0) * PI * r.powi(4), 0.5);
    }

    #[test]
    fn solid_cylinder() {
        // Decompose a radius-1, height-2 cylinder into 48 triangular prisms,
        // each split into 3 tetrahedra. Volume → πr²h, Izz → ½mr².
        let r = 1.0f64;
        let h = 2.0f64;
        let n = 48usize;
        let z0 = -h / 2.0;
        let z1 = h / 2.0;
        let mut v = VelGProps::new();
        for i in 0..n {
            let t1 = 2.0 * PI * i as f64 / n as f64;
            let t2 = 2.0 * PI * (i + 1) as f64 / n as f64;
            let ob = GpPnt::new(0.0, 0.0, z0);
            let ot = GpPnt::new(0.0, 0.0, z1);
            let b1 = GpPnt::new(r * t1.cos(), r * t1.sin(), z0);
            let b2 = GpPnt::new(r * t2.cos(), r * t2.sin(), z0);
            let tp1 = GpPnt::new(r * t1.cos(), r * t1.sin(), z1);
            let tp2 = GpPnt::new(r * t2.cos(), r * t2.sin(), z1);
            v.add_tetrahedron(&ob, &b1, &b2, &tp2, 1.0).unwrap();
            v.add_tetrahedron(&ob, &b1, &tp2, &tp1, 1.0).unwrap();
            v.add_tetrahedron(&ob, &ot, &tp1, &tp2, 1.0).unwrap();
        }
        let m = v.mass();
        assert_approx(m, PI * r * r * h, 0.05 * PI * r * r * h);
        let i = v.matrix_of_inertia();
        assert_approx(i.m[2][2], 0.5 * m * r * r, 0.05 * (0.5 * m * r * r));
        // A cylinder has an axis of symmetry.
        let pp = v.principal_properties().unwrap();
        assert!(pp.has_symmetry_axis());
        assert!(!pp.has_symmetry_point());
    }

    #[test]
    fn compose_with_density() {
        let mut pg = PGProps::new();
        pg.add_point(&GpPnt::new(1., 0., 0.)).unwrap();
        let mut sys = GProps::new();
        sys.add(&pg, 3.0).unwrap(); // density 3
        assert_approx(sys.mass(), 3.0, 1e-12);
        assert_approx(sys.centre_of_mass().x(), 1.0, 1e-12);
    }

    #[test]
    fn composite_two_locations() {
        // A unit mass at the origin and another at (10,0,0): COM at (5,0,0),
        // inertia about COM: each contributes m·5² to Iyy and Izz → 50 each.
        let mut a = GProps::new();
        a.add_point_mass(&GpPnt::new(0., 0., 0.), 1.0).unwrap();
        let mut b = GProps::new_at(GpPnt::new(10., 0., 0.));
        b.add_point_mass(&GpPnt::new(10., 0., 0.), 1.0).unwrap();
        let mut sys = GProps::new();
        sys.add(&a, 1.0).unwrap();
        sys.add(&b, 1.0).unwrap();
        assert_approx(sys.mass(), 2.0, 1e-12);
        assert_approx(sys.centre_of_mass().x(), 5.0, 1e-9);
        assert_approx(sys.centre_of_mass().y(), 0.0, 1e-9);
        let m = sys.matrix_of_inertia();
        assert_approx(m.m[1][1], 50.0, 1e-9);
        assert_approx(m.m[2][2], 50.0, 1e-9);
        assert_approx(m.m[0][0], 0.0, 1e-9);
    }

    #[test]
    fn principal_axes_orthonormal() {
        let mut v = VelGProps::new();
        v.add_box(&GpPnt::new(0., 0., 0.), &GpPnt::new(2., 3., 4.), 1.0).unwrap();
        let pp = v.principal_properties().unwrap();
        let (a, b, c) = (
            pp.first_axis_of_inertia(),
            pp.second_axis_of_inertia(),
            pp.third_axis_of_inertia(),
        );
        assert_approx(a.dot(&b), 0.0, 1e-9);
        assert_approx(a.dot(&c), 0.0, 1e-9);
        assert_approx(b.dot(&c), 0.0, 1e-9);
        assert_approx(a.magnitude(), 1.0, 1e-9);
        assert_approx(b.magnitude(), 1.0, 1e-9);
        assert_approx(c.magnitude(), 1.0, 1e-9);
    }

    #[test]
    fn barycentre_helpers() {
        let pts = [GpPnt::new(0., 0., 0.), GpPnt::new(2., 0., 0.), GpPnt::new(0., 2., 0.)];
        let b = PGProps::barycentre(&pts);
        assert_approx(b.x(), 2.0 / 3.0, 1e-12);
        assert_approx(b.y(), 2.0 / 3.0, 1e-12);
        let (mass, g) = PGProps::weighted_barycentre(&pts, &[1.0, 2.0, 1.0]).unwrap();
        assert_approx(mass, 4.0, 1e-12);
        assert_approx(g.x(), 1.0, 1e-12); // (0·1 + 2·2 + 0·1)/4
        assert_approx(g.y(), 0.5, 1e-12); // (0·1 + 0·2 + 2·1)/4
        assert!(PGProps::weighted_barycentre(&pts, &[1.0, 0.0, 1.0]).is_err());
    }

    #[test]
    fn plane_equation_from_points() {
        let eq = PEquation::from_points(
            &GpPnt::new(0., 0., 0.),
            &GpPnt::new(1., 0., 0.),
            &GpPnt::new(0., 1., 0.),
        )
        .unwrap();
        assert_approx(eq.a, 0.0, 1e-12);
        assert_approx(eq.b, 0.0, 1e-12);
        assert!(eq.c.abs() > 0.9);
        let p = GpPnt::new(1.0, 2.0, 5.0);
        assert_approx(eq.distance_to(&p), 5.0, 1e-12);
        let proj = eq.project(&p);
        assert_approx(proj.x(), 1.0, 1e-9);
        assert_approx(proj.y(), 2.0, 1e-9);
        assert_approx(proj.z(), 0.0, 1e-9);
        // Collinear points are rejected.
        assert!(
            PEquation::from_points(
                &GpPnt::new(0., 0., 0.),
                &GpPnt::new(1., 0., 0.),
                &GpPnt::new(2., 0., 0.),
            )
            .is_err()
        );
    }

    #[test]
    fn plane_equation_from_point_normal() {
        let eq =
            PEquation::from_point_normal(&GpPnt::new(0., 0., 3.), &GpVec::new(0., 0., 2.)).unwrap();
        assert_approx(eq.signed_distance(&GpPnt::new(0., 0., 5.)), 2.0, 1e-12);
        assert_approx(eq.distance_to(&GpPnt::new(1., 1., 3.)), 0.0, 1e-12);
        let n = eq.normal_unit();
        assert_approx(n.z(), 1.0, 1e-12);
        // Zero normal is rejected.
        assert!(PEquation::from_point_normal(&GpPnt::zero(), &GpVec::zero()).is_err());
    }

    #[test]
    fn plane_equation_to_plane() {
        let eq =
            PEquation::from_point_normal(&GpPnt::new(1., 2., 3.), &GpVec::new(1., 1., 1.)).unwrap();
        let pln = eq.to_plane();
        // Every point of the plane satisfies the equation.
        let loc = pln.location();
        assert_approx(eq.distance_to(&loc), 0.0, 1e-9);
        let axis = pln.axis();
        let nu = eq.normal_unit();
        assert!(axis.direction().xyz().dot(&nu.coord) > 0.999);
    }

    #[test]
    fn value_type_enum() {
        // Sanity: the enum carries the OCCT GProp_ValueType variants.
        let _ = [
            ValueType::Mass,
            ValueType::CenterMassX,
            ValueType::CenterMassY,
            ValueType::CenterMassZ,
            ValueType::InertiaXx,
            ValueType::InertiaYy,
            ValueType::InertiaZz,
            ValueType::InertiaXy,
            ValueType::InertiaXz,
            ValueType::InertiaYz,
            ValueType::Unknown,
        ];
        assert_eq!(ValueType::Mass as u8, 0);
        assert_eq!(ValueType::Unknown as u8, 10);
    }

    #[test]
    fn octahedron_volume() {
        // Regular octahedron with vertices (±1,0,0),(0,±1,0),(0,0,±1):
        // volume = 4/3, centroid at origin.
        let top = GpPnt::new(0., 0., 1.);
        let bottom = GpPnt::new(0., 0., -1.);
        let equator = [
            GpPnt::new(1., 0., 0.),
            GpPnt::new(0., 1., 0.),
            GpPnt::new(-1., 0., 0.),
            GpPnt::new(0., -1., 0.),
        ];
        let mut v = VelGProps::new();
        v.add_octahedron(&top, &bottom, &equator, 1.0).unwrap();
        assert_approx(v.mass(), 4.0 / 3.0, 1e-9);
        let c = v.centre_of_mass();
        assert_approx(c.x(), 0.0, 1e-9);
        assert_approx(c.y(), 0.0, 1e-9);
        assert_approx(c.z(), 0.0, 1e-9);
        // The regular octahedron has a point of symmetry.
        assert!(v.principal_properties().unwrap().has_symmetry_point());
    }
}
