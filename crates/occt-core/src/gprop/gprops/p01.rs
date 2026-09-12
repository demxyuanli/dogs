use super::prelude::*;
use super::*;

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
    pub(super) g: GpPnt,
    pub(super) loc: GpPnt,
    pub(super) dim: f64,
    pub(super) inertia: GpMat,
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
    pub(super) fn add_element(
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
    pub(super) i1: f64,
    pub(super) i2: f64,
    pub(super) i3: f64,
    pub(super) r1: f64,
    pub(super) r2: f64,
    pub(super) r3: f64,
    pub(super) v1: GpVec,
    pub(super) v2: GpVec,
    pub(super) v3: GpVec,
    pub(super) g: GpPnt,
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
    pub(super) props: GProps,
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
    pub(super) props: GProps,
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
    pub(super) props: GProps,
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
    pub(super) props: GProps,
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
