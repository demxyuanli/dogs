use super::prelude::*;
use super::*;

// ---------------------------------------------------------------------------
// Centroid
// ---------------------------------------------------------------------------

/// Exact centroid of an analytic solid: bbox center for boxes, surface-center
/// for spheres/toruses, mid-axis point for cylinders, ¾-height for cones.
/// Falls back to the mesh centroid otherwise.
pub fn analytic_centroid(shape: &TopoShape) -> Result<GpPnt, String> {
    let faces = faces_of(shape);
    if faces.is_empty() {
        return Err("analytic_centroid: shape has no faces".into());
    }
    // All-planar solid (a box): center of the vertex bounding box.
    if faces.iter().all(face_is_planar) {
        let verts = vertices_of(shape);
        if verts.len() >= 4 {
            let mut mn = [f64::INFINITY; 3];
            let mut mx = [f64::NEG_INFINITY; 3];
            for v in &verts {
                let p = vertex_position(v);
                mn[0] = mn[0].min(p.x());
                mn[1] = mn[1].min(p.y());
                mn[2] = mn[2].min(p.z());
                mx[0] = mx[0].max(p.x());
                mx[1] = mx[1].max(p.y());
                mx[2] = mx[2].max(p.z());
            }
            if mx[0] > mn[0] && mx[1] > mn[1] && mx[2] > mn[2] {
                return Ok(GpPnt::new(
                    (mn[0] + mx[0]) / 2.0,
                    (mn[1] + mx[1]) / 2.0,
                    (mn[2] + mx[2]) / 2.0,
                ));
            }
        }
    }
    // Sphere → center.
    for f in &faces {
        if let Some(surf) = BRepTool::face_surface(f) {
            if classify_surface_full(surf.as_ref()) == SurfaceKind::Sphere {
                if let Some(c) = sphere_center(surf.as_ref()) {
                    return Ok(c);
                }
            }
        }
    }
    // Cylinder → mid-axis point.
    for f in &faces {
        if let Some(surf) = BRepTool::face_surface(f) {
            if let Some((center, ax, _)) = cylinder_params(surf.as_ref()) {
                let h = shape_axial_extent(shape, &ax);
                return Ok(pnt_add_vec(center, &ax.multiplied_scalar(h / 2.0)));
            }
        }
    }
    // Cone → ¾ height from the base toward the apex.
    if let Some((r, _h)) = cone_radius_height(shape) {
        if let Some(pln) = faces.iter().find(|f| face_is_planar(f)).and_then(face_plane) {
            let ax1 = pln.axis();
            let ax = GpVec::from_xyz(ax1.direction().xyz());
            let base_center = pln.location();
            // The apex is the vertex farthest from the base center along the axis.
            let mut t_extreme = 0.0f64;
            for v in vertices_of(shape) {
                let p = vertex_position(&v);
                let t = (p.x() - base_center.x()) * ax.x()
                    + (p.y() - base_center.y()) * ax.y()
                    + (p.z() - base_center.z()) * ax.z();
                if t.abs() > t_extreme.abs() {
                    t_extreme = t;
                }
            }
            if t_extreme.abs() > 1e-12 {
                let _ = r;
                return Ok(pnt_add_vec(base_center, &ax.multiplied_scalar(t_extreme * 0.75)));
            }
        }
    }
    // Torus → center.
    for f in &faces {
        if let Some(surf) = BRepTool::face_surface(f) {
            if let Some((center, _, _, _)) = torus_params(surf.as_ref()) {
                return Ok(center);
            }
        }
    }
    // General fallback: mesh centroid.
    Ok(crate::brep_gprop::centroid(shape, 0.05).unwrap_or(GpPnt::zero()))
}

/// Combined exact mass properties.
pub fn analytic_properties(shape: &TopoShape) -> Result<AnalyticProps, String> {
    let surface_area = analytic_surface_area(shape)?;
    let volume = analytic_volume(shape)?;
    let centroid = analytic_centroid(shape)?;
    let exact = is_analytic(shape) && volume > 0.0;
    Ok(AnalyticProps { surface_area, volume, centroid, exact })
}

/// Whether every face of `shape` classifies as an analytic surface
/// (plane / sphere / cylinder / cone / torus).
pub fn is_analytic(shape: &TopoShape) -> bool {
    let faces = faces_of(shape);
    if faces.is_empty() {
        return false;
    }
    faces.iter().all(|f| match BRepTool::face_surface(f) {
        Some(surf) => classify_surface_full(surf.as_ref()) != SurfaceKind::Other,
        None => false,
    })
}

// ---------------------------------------------------------------------------
// Inertia tensor + principal axes
// ---------------------------------------------------------------------------

/// Symmetric 3×3 inertia tensor `(Ixx, Iyy, Izz, Ixy, Ixz, Iyz)` about the
/// centroid (off-diagonal entries are symmetric: `Ixy = Iyx`, etc.).
#[derive(Debug, Clone, Copy)]
pub struct InertiaTensor {
    pub ixx: f64,
    pub iyy: f64,
    pub izz: f64,
    pub ixy: f64,
    pub ixz: f64,
    pub iyz: f64,
}

impl InertiaTensor {
    /// The 3×3 symmetric matrix form `[[Ixx,Ixy,Ixz],[Ixy,Iyy,Iyz],[Ixz,Iyz,Izz]]`.
    pub fn matrix(&self) -> [[f64; 3]; 3] {
        [
            [self.ixx, self.ixy, self.ixz],
            [self.ixy, self.iyy, self.iyz],
            [self.ixz, self.iyz, self.izz],
        ]
    }

    /// The trace `Ixx + Iyy + Izz` (twice the sum of the squared radii of
    /// gyration about the coordinate axes).
    pub fn trace(&self) -> f64 {
        self.ixx + self.iyy + self.izz
    }

    /// The zero tensor.
    pub fn zero() -> InertiaTensor {
        InertiaTensor { ixx: 0.0, iyy: 0.0, izz: 0.0, ixy: 0.0, ixz: 0.0, iyz: 0.0 }
    }

    /// A pure diagonal tensor (a solid whose coordinate axes are already
    /// principal).
    pub fn diagonal(ixx: f64, iyy: f64, izz: f64) -> InertiaTensor {
        InertiaTensor { ixx, iyy, izz, ixy: 0.0, ixz: 0.0, iyz: 0.0 }
    }

    /// Build from a 3×3 matrix, symmetrizing the off-diagonal entries
    /// (`M[i][j]` and `M[j][i]` are averaged).
    pub fn from_matrix(m: &[[f64; 3]; 3]) -> InertiaTensor {
        InertiaTensor {
            ixx: m[0][0],
            iyy: m[1][1],
            izz: m[2][2],
            ixy: 0.5 * (m[0][1] + m[1][0]),
            ixz: 0.5 * (m[0][2] + m[2][0]),
            iyz: 0.5 * (m[1][2] + m[2][1]),
        }
    }

    /// Element-wise sum — combine the inertia tensors of disjoint pieces
    /// (both must be about the same point).
    pub fn sum(&self, other: &InertiaTensor) -> InertiaTensor {
        InertiaTensor {
            ixx: self.ixx + other.ixx,
            iyy: self.iyy + other.iyy,
            izz: self.izz + other.izz,
            ixy: self.ixy + other.ixy,
            ixz: self.ixz + other.ixz,
            iyz: self.iyz + other.iyz,
        }
    }

    /// Scale all moments by `factor` (e.g. a density ratio).
    pub fn scaled(&self, factor: f64) -> InertiaTensor {
        InertiaTensor {
            ixx: self.ixx * factor,
            iyy: self.iyy * factor,
            izz: self.izz * factor,
            ixy: self.ixy * factor,
            ixz: self.ixz * factor,
            iyz: self.iyz * factor,
        }
    }

    /// Shift the tensor to a parallel point displaced by `offset` via the
    /// parallel-axis theorem: `I' = I + m·(d²δ − d·dᵀ)`.
    pub fn translated(&self, mass: f64, offset: &GpVec) -> InertiaTensor {
        let (x, y, z) = (offset.x(), offset.y(), offset.z());
        let d2 = x * x + y * y + z * z;
        InertiaTensor {
            ixx: self.ixx + mass * (d2 - x * x),
            iyy: self.iyy + mass * (d2 - y * y),
            izz: self.izz + mass * (d2 - z * z),
            ixy: self.ixy + mass * (-x * y),
            ixz: self.ixz + mass * (-x * z),
            iyz: self.iyz + mass * (-y * z),
        }
    }

    /// Principal moments (eigenvalues, descending) and principal axes
    /// (eigenvectors) of this tensor, via the symmetric Jacobi solver.
    pub fn principal(&self) -> Result<(Vec<f64>, Vec<GpVec>), String> {
        let mat = self.matrix();
        let mut m = MathMatrix::new(1, 3, 1, 3);
        for i in 0..3 {
            for j in 0..3 {
                m.set_value(i + 1, j + 1, mat[i][j]);
            }
        }
        let (vals, vecs) = jacobi_eigen_symmetric(&m, 1e-9)?;
        let mut pairs: Vec<(f64, GpVec)> = (0..3)
            .map(|i| (vals[i], GpVec::new(vecs[i].value(1), vecs[i].value(2), vecs[i].value(3))))
            .collect();
        pairs.sort_by(|a, b| b.0.total_cmp(&a.0));
        Ok((
            pairs.iter().map(|(v, _)| *v).collect(),
            pairs.iter().map(|(_, v)| *v).collect(),
        ))
    }

    /// Radii of gyration about the three coordinate axes: `sqrt(Ixx/m)`,
    /// `sqrt(Iyy/m)`, `sqrt(Izz/m)`.
    pub fn radius_of_gyration(&self, mass: f64) -> (f64, f64, f64) {
        let m = mass.max(1e-30);
        ((self.ixx / m).sqrt(), (self.iyy / m).sqrt(), (self.izz / m).sqrt())
    }

    /// The moment of inertia about a unit direction `u`: `I = u·(I·u)`.
    pub fn moment_about_axis(&self, u: &GpVec) -> f64 {
        let (x, y, z) = (u.x(), u.y(), u.z());
        self.ixx * x * x + self.iyy * y * y + self.izz * z * z
            + 2.0 * (self.ixy * x * y + self.ixz * x * z + self.iyz * y * z)
    }

    /// The inertia matrix applied to a vector (matrix-vector product `I·v`).
    pub fn apply(&self, v: &GpVec) -> GpVec {
        let (x, y, z) = (v.x(), v.y(), v.z());
        GpVec::new(
            self.ixx * x + self.ixy * y + self.ixz * z,
            self.ixy * x + self.iyy * y + self.iyz * z,
            self.ixz * x + self.iyz * y + self.izz * z,
        )
    }

    /// Angular momentum of a rigid body with this inertia tensor and angular
    /// velocity `omega`: `L = I·ω`.
    pub fn angular_momentum(&self, omega: &GpVec) -> GpVec {
        self.apply(omega)
    }
}

/// Axis-aligned box extents `(dx, dy, dz)` from the vertex bounding box, for
/// a solid with the box's structural signature (8 vertices / 12 edges / 6 faces).
pub(super) fn box_dimensions(shape: &TopoShape) -> Option<(f64, f64, f64)> {
    let verts = vertices_of(shape);
    if verts.is_empty() {
        return None;
    }
    let mut mn = [f64::INFINITY; 3];
    let mut mx = [f64::NEG_INFINITY; 3];
    for v in &verts {
        let p = vertex_position(v);
        mn[0] = mn[0].min(p.x());
        mn[1] = mn[1].min(p.y());
        mn[2] = mn[2].min(p.z());
        mx[0] = mx[0].max(p.x());
        mx[1] = mx[1].max(p.y());
        mx[2] = mx[2].max(p.z());
    }
    if mx[0] > mn[0] && mx[1] > mn[1] && mx[2] > mn[2] {
        Some((mx[0] - mn[0], mx[1] - mn[1], mx[2] - mn[2]))
    } else {
        None
    }
}

/// Radius of a sphere primitive, from any spherical face.
pub(super) fn sphere_radius_from_shape(shape: &TopoShape) -> Option<f64> {
    for f in faces_of(shape) {
        if let Some(surf) = BRepTool::face_surface(&f) {
            if classify_surface_full(surf.as_ref()) == SurfaceKind::Sphere {
                let center = sphere_center(surf.as_ref())?;
                let r = surf.d0(0.0, 0.0).distance(&center);
                if r > 0.0 {
                    return Some(r);
                }
            }
        }
    }
    None
}

/// Cylinder `(axis, radius, height)` from the lateral surface and the axial
/// extent of the whole shape.
pub(super) fn cylinder_axis_radius_height(shape: &TopoShape) -> Option<(GpVec, f64, f64)> {
    for f in faces_of(shape) {
        if let Some(surf) = BRepTool::face_surface(&f) {
            if let Some((_, ax, r)) = cylinder_params(surf.as_ref()) {
                let h = shape_axial_extent(shape, &ax);
                if h > 0.0 {
                    return Some((ax, r, h));
                }
            }
        }
    }
    None
}

/// Cone `(axis, base radius, height)`; the axis is the planar base face normal.
pub(super) fn cone_axis_radius_height(shape: &TopoShape) -> Option<(GpVec, f64, f64)> {
    let (r, h) = cone_radius_height(shape)?;
    let faces = faces_of(shape);
    let base = faces.iter().find(|f| face_is_planar(f))?;
    let pln = face_plane(base)?;
    let ax = GpVec::from_xyz(pln.axis().direction().xyz());
    Some((ax, r, h))
}

/// Torus `(axis, major radius, minor radius)`.
pub(super) fn torus_axis_radii(shape: &TopoShape) -> Option<(GpVec, f64, f64)> {
    for f in faces_of(shape) {
        if let Some(surf) = BRepTool::face_surface(&f) {
            if let Some((_, ax, big_r, small_r)) = torus_params(surf.as_ref()) {
                return Some((ax, big_r, small_r));
            }
        }
    }
    None
}

/// Axis-aligned box extents `(dx, dy, dz)` of an analytic box solid.
pub fn analytic_box_dimensions(shape: &TopoShape) -> Option<(f64, f64, f64)> {
    box_dimensions(shape)
}

/// Radius of an analytic sphere solid.
pub fn analytic_sphere_radius(shape: &TopoShape) -> Option<f64> {
    sphere_radius_from_shape(shape)
}

/// `(center, axis, radius, height)` of an analytic cylinder solid.
pub fn analytic_cylinder_parameters(shape: &TopoShape) -> Option<(GpPnt, GpVec, f64, f64)> {
    for f in faces_of(shape) {
        if let Some(surf) = BRepTool::face_surface(&f) {
            if let Some((center, ax, r)) = cylinder_params(surf.as_ref()) {
                let h = shape_axial_extent(shape, &ax);
                if h > 0.0 {
                    return Some((center, ax, r, h));
                }
            }
        }
    }
    None
}

/// `(axis, base radius, height)` of an analytic cone solid.
pub fn analytic_cone_parameters(shape: &TopoShape) -> Option<(GpVec, f64, f64)> {
    cone_axis_radius_height(shape)
}

/// `(axis, major radius, minor radius)` of an analytic torus solid.
pub fn analytic_torus_parameters(shape: &TopoShape) -> Option<(GpVec, f64, f64)> {
    torus_axis_radii(shape)
}

/// Build the global-frame inertia tensor from the moments about the symmetry
/// axis (`i_axial`) and about any axis perpendicular to it (`i_radial`):
/// `I = I_rad·δ + (I_ax − I_rad)·(a⊗a)` for a unit axis `a`.
pub(super) fn tensor_from_axis(ax: &GpVec, i_radial: f64, i_axial: f64) -> InertiaTensor {
    let (x, y, z) = (ax.x(), ax.y(), ax.z());
    let d = i_axial - i_radial;
    InertiaTensor {
        ixx: i_radial + d * x * x,
        iyy: i_radial + d * y * y,
        izz: i_radial + d * z * z,
        ixy: d * x * y,
        ixz: d * x * z,
        iyz: d * y * z,
    }
}

/// Fallback inertia tensor: treat each vertex as carrying an equal share of the
/// total mass and accumulate `Σ mᵢ·(rᵢ²δ − rᵢ rᵢᵀ)` about the origin, then shift
/// to the centroid by the parallel-axis theorem. This is an approximation (the
/// exact analytic forms above are used whenever the shape classifies) and is
/// only a coarse estimate for non-analytic solids.
pub(super) fn vertex_mass_tensor(shape: &TopoShape, density: f64, volume: f64) -> Option<InertiaTensor> {
    let verts = vertices_of(shape);
    if verts.is_empty() {
        return None;
    }
    let m = volume * density;
    let m_i = m / verts.len() as f64;
    let mut ixx = 0.0;
    let mut iyy = 0.0;
    let mut izz = 0.0;
    let mut ixy = 0.0;
    let mut ixz = 0.0;
    let mut iyz = 0.0;
    for v in &verts {
        let p = vertex_position(v);
        let (x, y, z) = (p.x(), p.y(), p.z());
        let r2 = x * x + y * y + z * z;
        ixx += m_i * (r2 - x * x);
        iyy += m_i * (r2 - y * y);
        izz += m_i * (r2 - z * z);
        ixy += m_i * (-x * y);
        ixz += m_i * (-x * z);
        iyz += m_i * (-y * z);
    }
    // Parallel-axis theorem: I_cm = I_origin − m·(d²δ − d·dᵀ) with d = centroid.
    let c = analytic_centroid(shape).unwrap_or(GpPnt::zero());
    let (cx, cy, cz) = (c.x(), c.y(), c.z());
    let d2 = cx * cx + cy * cy + cz * cz;
    ixx -= m * (d2 - cx * cx);
    iyy -= m * (d2 - cy * cy);
    izz -= m * (d2 - cz * cz);
    ixy -= m * (-cx * cy);
    ixz -= m * (-cx * cz);
    iyz -= m * (-cy * cz);
    Some(InertiaTensor { ixx, iyy, izz, ixy, ixz, iyz })
}

/// Whether `shape` is an axis-aligned box solid: all-planar with the box's
/// structural signature (8 vertices / 12 edges / 6 faces).
pub fn shape_is_box(shape: &TopoShape) -> bool {
    let faces = faces_of(shape);
    faces.len() == 6
        && vertices_of(shape).len() == 8
        && edges_of(shape).len() == 12
        && faces.iter().all(face_is_planar)
}

/// Analytic kind of a closed solid, matching the inertia integration paths.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolidKind {
    Box,
    Sphere,
    Cylinder,
    Cone,
    Torus,
    Other,
}

/// Classify a closed analytic solid by the analytic surface it is built from.
pub fn classify_solid(shape: &TopoShape) -> SolidKind {
    if shape_is_box(shape) {
        return SolidKind::Box;
    }
    if sphere_radius_from_shape(shape).is_some() {
        return SolidKind::Sphere;
    }
    if cylinder_axis_radius_height(shape).is_some() {
        return SolidKind::Cylinder;
    }
    if cone_axis_radius_height(shape).is_some() {
        return SolidKind::Cone;
    }
    if torus_axis_radii(shape).is_some() {
        return SolidKind::Torus;
    }
    SolidKind::Other
}

/// 3×3 inertia tensor of `shape` about its centroid, with uniform `density`.
///
/// Analytic solids are integrated in closed form:
///
/// * box `dx×dy×dz` → `Ixx = m/12·(dy²+dz²)`, …;
/// * sphere radius `r` → `Ixx=Iyy=Izz = (2/5)·m·r²`;
/// * cylinder `r×h` → `Izz = (1/2)·m·r²`, `Ixx=Iyy = m·(3r²+h²)/12`;
/// * cone `r×h` → `Izz = (3/10)·m·r²`, `Ixx=Iyy = (3/20)·m·r² + (3/80)·m·h²`;
/// * torus `R×r` → `Izz = m·(R² + 3r²/4)`, `Ixx=Iyy = m·(R²/2 + 5r²/8)`.
///
/// Everything else falls back to the vertex-mass approximation. All results are
/// about the shape's centroid (no further offset is applied).
pub fn inertia_tensor(shape: &TopoShape, density: f64) -> Result<InertiaTensor, String> {
    let faces = faces_of(shape);
    if faces.is_empty() {
        return Err("inertia_tensor: shape has no faces".into());
    }
    if !density.is_finite() {
        return Err("inertia_tensor: density must be finite".into());
    }
    let volume = analytic_volume(shape)?;
    if volume <= 0.0 {
        return Err("inertia_tensor: shape has no volume".into());
    }
    let m = volume * density;

    // Box: axis-aligned rectangular solid.
    if shape_is_box(shape) {
        if let Some((dx, dy, dz)) = box_dimensions(shape) {
            return Ok(InertiaTensor {
                ixx: m / 12.0 * (dy * dy + dz * dz),
                iyy: m / 12.0 * (dx * dx + dz * dz),
                izz: m / 12.0 * (dx * dx + dy * dy),
                ixy: 0.0,
                ixz: 0.0,
                iyz: 0.0,
            });
        }
    }
    // Sphere.
    if let Some(r) = sphere_radius_from_shape(shape) {
        let s = (2.0 / 5.0) * m * r * r;
        return Ok(InertiaTensor { ixx: s, iyy: s, izz: s, ixy: 0.0, ixz: 0.0, iyz: 0.0 });
    }
    // Cylinder.
    if let Some((ax, r, h)) = cylinder_axis_radius_height(shape) {
        let i_axial = 0.5 * m * r * r;
        let i_radial = m * (3.0 * r * r + h * h) / 12.0;
        return Ok(tensor_from_axis(&ax, i_radial, i_axial));
    }
    // Cone.
    if let Some((ax, r, h)) = cone_axis_radius_height(shape) {
        let i_axial = (3.0 / 10.0) * m * r * r;
        let i_radial = (3.0 / 20.0) * m * r * r + (3.0 / 80.0) * m * h * h;
        return Ok(tensor_from_axis(&ax, i_radial, i_axial));
    }
    // Torus.
    if let Some((ax, big_r, small_r)) = torus_axis_radii(shape) {
        let i_axial = m * (big_r * big_r + 0.75 * small_r * small_r);
        let i_radial = m * (0.5 * big_r * big_r + (5.0 / 8.0) * small_r * small_r);
        return Ok(tensor_from_axis(&ax, i_radial, i_axial));
    }

    // Fallback: vertex-mass approximation.
    if let Some(t) = vertex_mass_tensor(shape, density, volume) {
        return Ok(t);
    }
    Err("inertia_tensor: unsupported shape".into())
}

/// The 3×3 symmetric inertia matrix of `shape` about its centroid:
/// `[[Ixx,Ixy,Ixz],[Ixy,Iyy,Iyz],[Ixz,Iyz,Izz]]`.
pub fn inertia_matrix(shape: &TopoShape, density: f64) -> Result<[[f64; 3]; 3], String> {
    Ok(inertia_tensor(shape, density)?.matrix())
}

/// Principal moments of inertia (eigenvalues, sorted descending) and the
/// corresponding principal axes (eigenvectors) of the inertia matrix.
///
/// The eigensolver is the symmetric Jacobi rotation method.
pub fn principal_inertia(shape: &TopoShape, density: f64) -> Result<(Vec<f64>, Vec<GpVec>), String> {
    let mat = inertia_matrix(shape, density)?;
    let mut m = MathMatrix::new(1, 3, 1, 3);
    for i in 0..3 {
        for j in 0..3 {
            m.set_value(i + 1, j + 1, mat[i][j]);
        }
    }
    let (vals, vecs) = jacobi_eigen_symmetric(&m, 1e-9)?;
    let mut pairs: Vec<(f64, GpVec)> = (0..3)
        .map(|i| {
            let v = GpVec::new(vecs[i].value(1), vecs[i].value(2), vecs[i].value(3));
            (vals[i], v)
        })
        .collect();
    pairs.sort_by(|a, b| b.0.total_cmp(&a.0));
    let moments: Vec<f64> = pairs.iter().map(|(v, _)| *v).collect();
    let axes: Vec<GpVec> = pairs.iter().map(|(_, v)| *v).collect();
    Ok((moments, axes))
}

/// The 3×3 rotation matrix whose columns are the shape's principal axes
/// (an orthonormal frame for the principal coordinate system).
pub fn principal_axes_matrix(shape: &TopoShape, density: f64) -> Result<[[f64; 3]; 3], String> {
    let (_, axes) = principal_inertia(shape, density)?;
    Ok([
        [axes[0].x(), axes[1].x(), axes[2].x()],
        [axes[0].y(), axes[1].y(), axes[2].y()],
        [axes[0].z(), axes[1].z(), axes[2].z()],
    ])
}

/// Inertia tensor about an arbitrary point `p`, obtained from the centroid
/// tensor via the parallel-axis theorem: `I_p = I_cm + m·(d²δ − d·dᵀ)` with
/// `d = p − centroid`.
pub fn inertia_tensor_at(shape: &TopoShape, density: f64, p: &GpPnt) -> Result<InertiaTensor, String> {
    let t = inertia_tensor(shape, density)?;
    let volume = analytic_volume(shape)?;
    let mass = volume * density;
    let c = analytic_centroid(shape)?;
    let offset = GpVec::from_pnts(&c, p);
    Ok(t.translated(mass, &offset))
}

/// Inertia tensor of a set of disjoint solids about their combined centre of
/// mass.
///
/// Each solid's tensor is computed about its own centroid (via
/// [`inertia_tensor`]), shifted to the combined centre of mass by the
/// parallel-axis theorem, and summed element-wise. Returns the combined tensor
/// and the combined centroid.
pub fn inertia_tensor_composite(
    shapes: &[&TopoShape],
    density: f64,
) -> Result<(InertiaTensor, GpPnt), String> {
    if shapes.is_empty() {
        return Err("inertia_tensor_composite: no shapes".into());
    }
    let mut total_mass = 0.0;
    let mut sum_mc = occt_core::gp::GpXyz::zero();
    let mut entries: Vec<(InertiaTensor, f64, GpPnt)> = Vec::with_capacity(shapes.len());
    for s in shapes {
        let volume = analytic_volume(s)?;
        let mass = volume * density;
        let c = analytic_centroid(s)?;
        sum_mc = sum_mc.added(&c.coord.multiplied(mass));
        total_mass += mass;
        entries.push((inertia_tensor(s, density)?, mass, c));
    }
    if total_mass <= 0.0 {
        return Err("inertia_tensor_composite: no mass".into());
    }
    let center = GpPnt::from_xyz(&sum_mc.divided(total_mass));
    let mut acc = InertiaTensor::zero();
    for (t, mass, c) in &entries {
        let offset = GpVec::from_pnts(c, &center);
        acc = acc.sum(&t.translated(*mass, &offset));
    }
    Ok((acc, center))
}

/// Density-weighted mass of an analytic solid: `density × volume`.
pub fn analytic_mass(shape: &TopoShape, density: f64) -> Result<f64, String> {
    if !density.is_finite() || density < 0.0 {
        return Err("analytic_mass: density must be finite and non-negative".into());
    }
    Ok(analytic_volume(shape)? * density)
}

// ---------------------------------------------------------------------------
// Analytic curve length
// ---------------------------------------------------------------------------

/// Analytic kind of a 3D curve, reconstructed from sampled geometry invariants
/// (the port keeps no `Any`-typed geometry, so the concrete `Geom*` type cannot
/// be downcast).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CurveKind {
    /// Straight segment: constant tangent direction, zero second derivative.
    Line,
    /// Circular arc: constant tangent speed and curvature, planar binormal.
    Circle,
    /// Bounded planar curve that is not a circle (ellipse, etc.).
    Ellipse,
    /// Anything else (B-spline, helix, …), handled by chord-length sampling.
    Other,
}

/// Classify a curve over the parameter window `[a, b]`:
///
/// * `Line` — the second derivative vanishes at every sample;
/// * `Circle` — constant tangent speed `r` and constant nonzero curvature with
///   a planar binormal direction, so the arc length is exactly `r·|b−a|`;
/// * `Ellipse` — the sampled points are coplanar but the tangent speed is not
///   constant (an ellipse's arc length needs an elliptic integral);
/// * `Other` — non-planar or unbounded curves.
pub fn classify_curve(curve: &dyn occt_geom::Curve, a: f64, b: f64) -> CurveKind {
    if is_line_curve(curve, a, b) {
        CurveKind::Line
    } else if circle_radius(curve, a, b).is_some() {
        CurveKind::Circle
    } else if curve_is_planar(curve, a, b) {
        CurveKind::Ellipse
    } else {
        CurveKind::Other
    }
}
