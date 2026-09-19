use super::prelude::*;
use super::*;

// ===========================================================================
// Meridional (axisymmetric) blend geometry
// ===========================================================================

/// A surface of revolution's cross-section (meridional) radius as a function of
/// height along a common axis. `NaN` marks heights off the surface.
#[derive(Debug, Clone, Copy)]
pub(super) enum MeridionalShape {
    /// Constant radius (cylinder).
    Cylinder { rho: f64 },
    /// Cone: radius = (z − apex_z)·d·tan on the nappe where that is positive.
    Cone { apex_z: f64, d: f64, tan: f64, sin_alpha: f64 },
    /// Sphere: radius = √(rho² − (z − center_z)²).
    Sphere { center_z: f64, rho: f64 },
}

impl MeridionalShape {
    /// Cross-section radius of the *original* surface at height `z`.
    pub(super) fn radius(&self, z: f64) -> f64 {
        match *self {
            MeridionalShape::Cylinder { rho } => rho,
            MeridionalShape::Cone { apex_z, d, tan, .. } => {
                let r = (z - apex_z) * d * tan;
                if r > 0.0 {
                    r
                } else {
                    f64::NAN
                }
            }
            MeridionalShape::Sphere { center_z, rho } => {
                let d2 = rho * rho - (z - center_z) * (z - center_z);
                if d2 >= 0.0 {
                    d2.sqrt()
                } else {
                    f64::NAN
                }
            }
        }
    }

    /// Cross-section radius of the surface offset outward by `r`.
    pub(super) fn offset_radius(&self, z: f64, r: f64) -> f64 {
        match *self {
            MeridionalShape::Cylinder { rho } => rho + r,
            MeridionalShape::Cone { apex_z, d, tan, sin_alpha } => {
                // The outward normal of a cone has a +r·sin(α) component
                // opposite the base, so the offset apex moves by r/sin(α)
                // opposite the base direction.
                let apex_z_off = apex_z - (r / sin_alpha) * d;
                let ro = (z - apex_z_off) * d * tan;
                if ro > 0.0 {
                    ro
                } else {
                    f64::NAN
                }
            }
            MeridionalShape::Sphere { center_z, rho } => {
                let d2 = (rho + r) * (rho + r) - (z - center_z) * (z - center_z);
                if d2 >= 0.0 {
                    d2.sqrt()
                } else {
                    f64::NAN
                }
            }
        }
    }

    /// z-interval on which `offset_radius(z, r)` is finite (the offset nappe).
    /// `None` marks an unbounded end. Used to bound the offset-centreline search.
    pub(super) fn offset_z_range(&self, r: f64) -> (Option<f64>, Option<f64>) {
        match *self {
            MeridionalShape::Cylinder { .. } => (None, None),
            MeridionalShape::Cone { apex_z, d, sin_alpha, .. } => {
                let z0 = apex_z - (r / sin_alpha) * d;
                if d > 0.0 {
                    (Some(z0), None)
                } else {
                    (None, Some(z0))
                }
            }
            MeridionalShape::Sphere { center_z, rho } => {
                (Some(center_z - (rho + r)), Some(center_z + (rho + r)))
            }
        }
    }
}

/// Closest point on a meridional curve `rho(z)` to the ball centre
/// `(rho_c, z_c)`, by coarse grid + ternary refinement. Returns `(rho, z)`.
pub(super) fn closest_point_meridional(
    rho: &dyn Fn(f64) -> f64,
    rho_c: f64,
    z_c: f64,
    z_lo: f64,
    z_hi: f64,
) -> Option<(f64, f64)> {
    let n = 256;
    let mut best = (f64::NAN, f64::NAN, f64::INFINITY);
    for i in 0..=n {
        let z = z_lo + (z_hi - z_lo) * i as f64 / n as f64;
        let r = rho(z);
        if !r.is_finite() || r < 0.0 {
            continue;
        }
        let d2 = (r - rho_c) * (r - rho_c) + (z - z_c) * (z - z_c);
        if d2 < best.2 {
            best = (r, z, d2);
        }
    }
    if !best.0.is_finite() {
        return None;
    }
    let mut lo = (best.1 - (z_hi - z_lo) / n as f64).max(z_lo);
    let mut hi = (best.1 + (z_hi - z_lo) / n as f64).min(z_hi);
    for _ in 0..96 {
        let m1 = lo + (hi - lo) / 3.0;
        let m2 = hi - (hi - lo) / 3.0;
        let d1 = (rho(m1) - rho_c) * (rho(m1) - rho_c) + (m1 - z_c) * (m1 - z_c);
        let d2v = (rho(m2) - rho_c) * (rho(m2) - rho_c) + (m2 - z_c) * (m2 - z_c);
        if d1 < d2v {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    let z = 0.5 * (lo + hi);
    Some((rho(z), z))
}

/// Compute the tangency circle of the rolling ball against a surface of
/// revolution described by its meridional radius function `rho`.
///
/// `center`/`axis` are the centreline circle (on the axis) and its unit
/// normal; the ball centre is at radius `r_c` in the plane of the centreline
/// (height 0). The search for the contact height is bounded by `z_lo`/`z_hi`
/// (relative to the centreline plane).
pub(super) fn contact_circle_meridional(
    center: &GpPnt,
    axis: &GpVec,
    r_c: f64,
    r: f64,
    rho: &dyn Fn(f64) -> f64,
    z_lo: f64,
    z_hi: f64,
) -> Result<ContactCircleGeom, String> {
    let (rho_t, z_t) = closest_point_meridional(rho, r_c, 0.0, z_lo, z_hi)
        .ok_or("fillet_curved: cannot locate the ball contact on a curved face")?;
    let dist = ((rho_t - r_c) * (rho_t - r_c) + z_t * z_t).sqrt();
    if (dist - r).abs() > 1e-2 * r.max(1.0) {
        return Err(format!(
            "fillet_curved: rolling ball of radius {r} does not reach this face \
             (contact distance {dist})"
        ));
    }
    let cos_v = ((rho_t - r_c) / r).clamp(-1.0, 1.0);
    let sin_v = (z_t / r).clamp(-1.0, 1.0);
    let v = sin_v.atan2(cos_v);
    Ok(ContactCircleGeom {
        center: center.translated_vec(&axis.multiplied_scalar(r * v.sin())),
        normal: *axis,
        radius: r_c + r * v.cos(),
    })
}

/// Tangency circle of the ball against a plane perpendicular to the axis. The
/// plane is at height `plane_z` relative to the centreline plane.
pub(super) fn contact_circle_plane(
    center: &GpPnt,
    axis: &GpVec,
    r_c: f64,
    r: f64,
    plane_z: f64,
) -> Result<ContactCircleGeom, String> {
    let t = (plane_z / r).clamp(-1.0, 1.0);
    let v = t.asin();
    Ok(ContactCircleGeom {
        center: center.translated_vec(&axis.multiplied_scalar(r * v.sin())),
        normal: *axis,
        radius: r_c + r * v.cos(),
    })
}

/// Find the height where two offset meridional radius functions are equal
/// (the rolling-ball centreline circle). Returns `(z, radius)`.
pub(super) fn offset_centerline_height(
    rho1: &dyn Fn(f64) -> f64,
    rho2: &dyn Fn(f64) -> f64,
    z_lo: f64,
    z_hi: f64,
) -> Result<(f64, f64), String> {
    let n = 1024;
    let mut prev = rho1(z_lo) - rho2(z_lo);
    for i in 1..=n {
        let z = z_lo + (z_hi - z_lo) * i as f64 / n as f64;
        let cur = rho1(z) - rho2(z);
        if prev.is_finite() && cur.is_finite() {
            if (prev < 0.0) != (cur < 0.0) {
                let mut lo = z_lo + (z_hi - z_lo) * (i - 1) as f64 / n as f64;
                let mut hi = z;
                let mut flo = prev;
                for _ in 0..80 {
                    let mid = 0.5 * (lo + hi);
                    let fm = rho1(mid) - rho2(mid);
                    if (flo < 0.0) == (fm < 0.0) {
                        lo = mid;
                        flo = fm;
                    } else {
                        hi = mid;
                    }
                }
                let zc = 0.5 * (lo + hi);
                let rc = rho1(zc);
                if rc.is_finite() && rc > 1e-9 {
                    return Ok((zc, rc));
                }
            }
        }
        prev = cur;
    }
    Err("fillet_curved: the two offset surfaces do not intersect in a circle".to_string())
}

/// Finite z-window bracketing the offset-centreline crossing of two coaxial
/// surfaces of revolution. Intersects the two offset nappes; an unbounded end
/// (cylinder, or a cone whose nappe is open on that side) is clamped by walking
/// outward until the faster-growing offset radius exceeds the other by a factor
/// of 8 — the crossing (unique for these coaxial pairs) lies inside.
pub(super) fn overlap_window(m1: &MeridionalShape, m2: &MeridionalShape, r: f64) -> Result<(f64, f64), String> {
    let (l1, h1) = m1.offset_z_range(r);
    let (l2, h2) = m2.offset_z_range(r);
    let lo = match (l1, l2) {
        (Some(a), Some(b)) => Some(a.max(b)),
        (Some(a), None) => Some(a),
        (None, Some(b)) => Some(b),
        (None, None) => None,
    };
    let hi = match (h1, h2) {
        (Some(a), Some(b)) => Some(a.min(b)),
        (Some(a), None) => Some(a),
        (None, Some(b)) => Some(b),
        (None, None) => None,
    };
    if let (Some(a), Some(b)) = (lo, hi) {
        if a >= b {
            return Err("fillet_curved: the two offset nappes do not overlap along the axis".to_string());
        }
    }
    let (mut lo, mut hi) = (lo, hi);
    if lo.is_some() && hi.is_some() {
        return Ok((lo.unwrap(), hi.unwrap()));
    }
    let (anchor, unbounded_below) = match (lo, hi) {
        (Some(a), _) => (a, false),
        (None, Some(b)) => (b, true),
        (None, None) => {
            return Err("fillet_curved: both offset surfaces are unbounded along the axis".to_string())
        }
    };
    let rho_anchor = m1.offset_radius(anchor, r).max(m2.offset_radius(anchor, r)).max(1e-9);
    let e = rho_anchor.max(1.0) * 1e-4;
    // Probe strictly inside the overlap (the boundary itself may be a nappe
    // apex, where `offset_radius` is NaN) to measure the growth rate.
    let sgn = if unbounded_below { -1.0 } else { 1.0 };
    let ra1 = m1.offset_radius(anchor + sgn * e, r);
    let ra2 = m2.offset_radius(anchor + sgn * e, r);
    let rb1 = m1.offset_radius(anchor + sgn * 2.0 * e, r);
    let rb2 = m2.offset_radius(anchor + sgn * 2.0 * e, r);
    let slope = {
        let s1 = ((rb1 - ra1) / e).abs();
        let s2 = ((rb2 - ra2) / e).abs();
        (if s1.is_finite() { s1 } else { 1e-9 })
            .max(if s2.is_finite() { s2 } else { 1e-9 })
            .max(1e-9)
    };
    let span = 8.0 * rho_anchor / slope;
    if lo.is_none() {
        lo = Some(anchor - span);
    }
    if hi.is_none() {
        hi = Some(anchor + span);
    }
    let (a, b) = (lo.unwrap(), hi.unwrap());
    if a >= b || !a.is_finite() || !b.is_finite() {
        return Err("fillet_curved: the two offset surfaces have no reachable centerline".to_string());
    }
    Ok((a, b))
}

/// Build a `MeridionalShape::Cone` from a `ConeInfo` in the frame of `axis`.
pub(super) fn cone_meridional(cone: &ConeInfo, axis: &GpVec, ref_pt: &GpPnt) -> MeridionalShape {
    let d = cone.axis.dot(axis);
    let apex_z = GpVec::from_pnts(ref_pt, &cone.apex).dot(axis);
    MeridionalShape::Cone {
        apex_z,
        d,
        tan: cone.semi_angle.tan(),
        sin_alpha: cone.semi_angle.sin(),
    }
}

// ===========================================================================
// Blend geometry per supported pair
// ===========================================================================

/// Build the torus band blend for a plane + sphere pair.
pub(super) fn plane_sphere_blend(pln: &PlaneInfo, sph: &SphereInfo, r: f64) -> Result<TorusBlendGeom, String> {
    let circle = plane_sphere_centerline(pln, sph, r)?;
    let n = circle.normal;
    let delta = n.dot(&GpVec::from_pnts(&pln.origin, &sph.center)) - r;
    // Tube angle of the sphere-tangency circle: cos v = −rho/(R_s+r), sin v = δ/(R_s+r).
    let v_sphere = delta.atan2(-circle.radius);
    let v_plane = 3.0 * PI / 2.0; // torus bottom: touches the plane
    let v_sphere = if v_sphere < 0.0 { v_sphere + 2.0 * PI } else { v_sphere };
    let plane_contact = ContactCircleGeom {
        center: circle.center.translated_vec(&n.multiplied_scalar(r * (v_plane.sin()))),
        normal: n,
        radius: circle.radius + r * v_plane.cos(),
    };
    let sphere_contact = ContactCircleGeom {
        center: circle.center.translated_vec(&n.multiplied_scalar(r * v_sphere.sin())),
        normal: n,
        radius: circle.radius + r * v_sphere.cos(),
    };
    Ok(TorusBlendGeom {
        center: circle.center,
        axis: n,
        major: circle.radius,
        minor: r,
        contact_a: plane_contact,
        contact_b: sphere_contact,
    })
}

/// Build the torus band blend for a plane + cylinder pair. The cylinder axis
/// must be perpendicular to the plane (a "standing" cylinder); otherwise the
/// offset-plane / offset-cylinder intersection is an ellipse and the exact
/// torus construction does not apply.
pub(super) fn plane_cylinder_blend(pln: &PlaneInfo, cyl: &CylinderInfo, r: f64, tol: f64) -> Result<TorusBlendGeom, String> {
    let n = pln.normal.normalized();
    let a = cyl.axis.normalized();
    if n.xyz().crossed(&a.xyz()).modulus() > tol.max(1e-6) {
        return Err(
            "fillet_curved: plane+cylinder blend requires the cylinder axis to be perpendicular to the plane"
                .to_string(),
        );
    }
    let rho = cyl.radius + r;
    // Project the axis point onto the offset plane n·x = n·O + r.
    let n_origin = n.dot(&GpVec::from_pnts(&GpPnt::zero(), &pln.origin));
    let n_axis = n.dot(&GpVec::from_pnts(&GpPnt::zero(), &cyl.axis_origin));
    let center = cyl.axis_origin.translated_vec(&n.multiplied_scalar(n_origin + r - n_axis));
    let v_cyl = PI; // inner equator: touches the cylinder
    let v_plane = 3.0 * PI / 2.0; // bottom: touches the plane
    let plane_contact = ContactCircleGeom {
        center: center.translated_vec(&n.multiplied_scalar(r * v_plane.sin())),
        normal: n,
        radius: rho + r * v_plane.cos(),
    };
    let cyl_contact = ContactCircleGeom {
        center: center.translated_vec(&n.multiplied_scalar(r * v_cyl.sin())),
        normal: n,
        radius: rho + r * v_cyl.cos(),
    };
    Ok(TorusBlendGeom {
        center,
        axis: n,
        major: rho,
        minor: r,
        contact_a: plane_contact,
        contact_b: cyl_contact,
    })
}

/// Build the torus band blend for a sphere + sphere pair.
pub(super) fn sphere_sphere_blend(s1: &SphereInfo, s2: &SphereInfo, r: f64) -> Result<TorusBlendGeom, String> {
    let circle = sphere_sphere_centerline(s1, s2, r)?;
    let dv = GpVec::from_pnts(&s1.center, &s2.center);
    let d = dv.magnitude();
    let dhat = dv.normalized();
    let r1 = s1.radius + r;
    let r2 = s2.radius + r;
    let a = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
    // Tube angles: cos v_i = −rho / (R_i + r); sin v_1 = −a / (R_1 + r),
    // sin v_2 = (d − a) / (R_2 + r).
    let v1 = (-a).atan2(-circle.radius);
    let v2 = (d - a).atan2(-circle.radius);
    let contact1 = ContactCircleGeom {
        center: circle.center.translated_vec(&dhat.multiplied_scalar(r * v1.sin())),
        normal: dhat,
        radius: circle.radius + r * v1.cos(),
    };
    let contact2 = ContactCircleGeom {
        center: circle.center.translated_vec(&dhat.multiplied_scalar(r * v2.sin())),
        normal: dhat,
        radius: circle.radius + r * v2.cos(),
    };
    Ok(TorusBlendGeom {
        center: circle.center,
        axis: dhat,
        major: circle.radius,
        minor: r,
        contact_a: contact1,
        contact_b: contact2,
    })
}

/// Build the torus band blend for a plane + cone pair (concave base crease).
pub(super) fn plane_cone_blend(pln: &PlaneInfo, cone: &ConeInfo, r: f64, tol: f64) -> Result<TorusBlendGeom, String> {
    let n = pln.normal.normalized();
    let a = cone.axis.normalized();
    if n.xyz().crossed(&a.xyz()).modulus() > tol.max(1e-6) {
        return Err(
            "fillet_curved: plane+cone blend requires the cone axis perpendicular to the plane"
                .to_string(),
        );
    }
    let axis = n;
    let ref_pt = cone.apex;
    let m_cone = cone_meridional(cone, &axis, &ref_pt);
    let z_p = GpVec::from_pnts(&ref_pt, &pln.origin).dot(&axis);
    let z_c = z_p + r;
    let rc = m_cone.offset_radius(z_c, r);
    if !rc.is_finite() || rc <= 1e-9 {
        return Err("fillet_curved: plane+cone offset surfaces do not yield a centerline".to_string());
    }
    let center = ref_pt.translated_vec(&axis.multiplied_scalar(z_c));
    let contact_plane = contact_circle_plane(&center, &axis, rc, r, z_p - z_c)?;
    let rho_cone = |z: f64| m_cone.radius(z_c + z);
    let contact_cone = contact_circle_meridional(&center, &axis, rc, r, &rho_cone, -1e3, 1e3)?;
    Ok(TorusBlendGeom {
        center,
        axis,
        major: rc,
        minor: r,
        contact_a: contact_plane,
        contact_b: contact_cone,
    })
}

/// Build the torus band blend for a coaxial cylinder + cone pair.
pub(super) fn cylinder_cone_blend(cyl: &CylinderInfo, cone: &ConeInfo, r: f64, tol: f64) -> Result<TorusBlendGeom, String> {
    let a_cyl = cyl.axis.normalized();
    let a_cone = cone.axis.normalized();
    if a_cyl.xyz().crossed(&a_cone.xyz()).modulus() > tol.max(1e-6) {
        return Err("fillet_curved: cylinder+cone blend requires coaxial surfaces".to_string());
    }
    let axis = a_cone;
    let ref_pt = cone.apex;
    let m_cyl = MeridionalShape::Cylinder { rho: cyl.radius };
    let m_cone = cone_meridional(cone, &axis, &ref_pt);
    let rho_cyl = |z: f64| m_cyl.offset_radius(z, r);
    let rho_cone = |z: f64| m_cone.offset_radius(z, r);
    let (z_lo, z_hi) = overlap_window(&m_cyl, &m_cone, r)?;
    let (zc, rc) = offset_centerline_height(&rho_cyl, &rho_cone, z_lo, z_hi)?;
    let center = ref_pt.translated_vec(&axis.multiplied_scalar(zc));
    let rho_cyl_s = |z: f64| m_cyl.radius(zc + z);
    let contact_cyl = contact_circle_meridional(&center, &axis, rc, r, &rho_cyl_s, -1e3, 1e3)?;
    let rho_cone_s = |z: f64| m_cone.radius(zc + z);
    let contact_cone = contact_circle_meridional(&center, &axis, rc, r, &rho_cone_s, -1e3, 1e3)?;
    Ok(TorusBlendGeom {
        center,
        axis,
        major: rc,
        minor: r,
        contact_a: contact_cyl,
        contact_b: contact_cone,
    })
}

/// Build the torus band blend for a sphere + coaxial cone pair.
pub(super) fn cone_sphere_blend(cone: &ConeInfo, sph: &SphereInfo, r: f64, tol: f64) -> Result<TorusBlendGeom, String> {
    let a_cone = cone.axis.normalized();
    let cs = GpVec::from_pnts(&cone.apex, &sph.center);
    let radial = cs.subtracted(&a_cone.multiplied_scalar(cs.dot(&a_cone)));
    if radial.magnitude() > tol.max(1e-6) {
        return Err(
            "fillet_curved: cone+sphere blend requires the sphere centre on the cone axis".to_string(),
        );
    }
    let axis = a_cone;
    let ref_pt = cone.apex;
    let m_cone = cone_meridional(cone, &axis, &ref_pt);
    let m_sph = MeridionalShape::Sphere { center_z: cs.dot(&axis), rho: sph.radius };
    let rho_cone = |z: f64| m_cone.offset_radius(z, r);
    let rho_sph = |z: f64| m_sph.offset_radius(z, r);
    let (z_lo, z_hi) = overlap_window(&m_cone, &m_sph, r)?;
    let (zc, rc) = offset_centerline_height(&rho_cone, &rho_sph, z_lo, z_hi)?;
    let center = ref_pt.translated_vec(&axis.multiplied_scalar(zc));
    let rho_cone_s = |z: f64| m_cone.radius(zc + z);
    let contact_cone = contact_circle_meridional(&center, &axis, rc, r, &rho_cone_s, -1e3, 1e3)?;
    let rho_sph_s = |z: f64| m_sph.radius(zc + z);
    let contact_sph = contact_circle_meridional(&center, &axis, rc, r, &rho_sph_s, -1e3, 1e3)?;
    Ok(TorusBlendGeom {
        center,
        axis,
        major: rc,
        minor: r,
        contact_a: contact_cone,
        contact_b: contact_sph,
    })
}

/// Build the torus band blend for a pair of coaxial cones.
pub(super) fn cone_cone_blend(c1: &ConeInfo, c2: &ConeInfo, r: f64, tol: f64) -> Result<TorusBlendGeom, String> {
    let a1 = c1.axis.normalized();
    let a2 = c2.axis.normalized();
    if a1.xyz().crossed(&a2.xyz()).modulus() > tol.max(1e-6) {
        return Err("fillet_curved: cone+cone blend requires coaxial cones".to_string());
    }
    let axis = a1;
    let ref_pt = c1.apex;
    let m1 = cone_meridional(c1, &axis, &ref_pt);
    let m2 = cone_meridional(c2, &axis, &ref_pt);
    let rho1 = |z: f64| m1.offset_radius(z, r);
    let rho2 = |z: f64| m2.offset_radius(z, r);
    let (z_lo, z_hi) = overlap_window(&m1, &m2, r)?;
    let (zc, rc) = offset_centerline_height(&rho1, &rho2, z_lo, z_hi)?;
    let center = ref_pt.translated_vec(&axis.multiplied_scalar(zc));
    let rho1s = |z: f64| m1.radius(zc + z);
    let ca = contact_circle_meridional(&center, &axis, rc, r, &rho1s, -1e3, 1e3)?;
    let rho2s = |z: f64| m2.radius(zc + z);
    let cb = contact_circle_meridional(&center, &axis, rc, r, &rho2s, -1e3, 1e3)?;
    Ok(TorusBlendGeom {
        center,
        axis,
        major: rc,
        minor: r,
        contact_a: ca,
        contact_b: cb,
    })
}

// ===========================================================================
// Topology rebuild
// ===========================================================================

/// How to rebuild an adjacent face.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum AdjacentKind {
    Plane,
    Sphere,
    Cylinder,
    Cone,
}

/// Build a full circular edge on the circle `c`, with its seam at `seam` (a
/// point on the circle, angle 0).
pub(super) fn build_circle_edge(b: &TopoBuilder, c: &ContactCircleGeom, seam: &GpPnt) -> Result<Edge, String> {
    let nd = GpDir::from_xyz(&c.normal.xyz()).map_err(|e| e.to_string())?;
    let xd = GpDir::from_vec(&GpVec::from_pnts(&c.center, seam))
        .map_err(|_| "fillet_curved: tangency circle seam coincides with its centre".to_string())?;
    let ax2 = GpAx2::new(c.center, nd, xd).map_err(|e| format!("fillet_curved: circle frame: {e}"))?;
    let mut e = b.make_edge_circle(&ax2, c.radius, 0.0, 2.0 * PI);
    let v = b.make_vertex(*seam, 0.0);
    b.add_edge_vertices(&mut e, &v, &v);
    Ok(e)
}

/// A default seam point for a tangency circle (angle 0 along a fixed
/// perpendicular). Used where the seam angle is not constrained by a neighbour.
pub(super) fn default_seam(c: &ContactCircleGeom) -> GpPnt {
    let e1 = perp_dir(&c.normal);
    c.center.translated_vec(&e1.multiplied_scalar(c.radius))
}

/// The closed (circular) edge of a cylinder lateral face that is not the
/// filleted edge — i.e. the top rim circle.
pub(super) fn find_top_circle(face: &Face, fillet_edge: &Edge) -> Result<Edge, String> {
    for w in wires_of_face(face) {
        for e in edges_of_wire(&w) {
            if is_same(&e.0, &fillet_edge.0) {
                continue;
            }
            if let Some((a, b)) = BRepTool::edge_vertices(&e) {
                if a.distance(&b) < 1e-9 {
                    return Ok(e.clone());
                }
            }
        }
    }
    Err("fillet_curved: cylinder lateral face has no top rim circle".to_string())
}

/// The point on a cylinder-tangency circle that shares a generator with the
/// cylinder's top rim seam, so the rebuilt lateral face's seam sub-edge stays a
/// straight generator.
pub(super) fn cylinder_contact_seam(face: &Face, fillet_edge: &Edge, contact: &ContactCircleGeom) -> Result<GpPnt, String> {
    let top = find_top_circle(face, fillet_edge)?;
    let (s_top, _) = BRepTool::edge_vertices(&top).ok_or("fillet_curved: top rim has no vertices")?;
    let v = GpVec::from_pnts(&contact.center, &s_top);
    let along = v.dot(&contact.normal);
    let radial = v.subtracted(&contact.normal.multiplied_scalar(along));
    let m = radial.magnitude();
    if m < 1e-9 {
        return Err("fillet_curved: cylinder top seam lies on the axis".to_string());
    }
    Ok(contact.center.translated_vec(&radial.multiplied_scalar(contact.radius / m)))
}

/// Rebuild a planar adjacent face: keep every original wire except the one
/// containing the filleted edge, and add a new hole wire made of the enlarged
/// plane-tangency circle.
pub(super) fn rebuild_plane_face(face: &Face, fillet_edge: &Edge, new_circle: &Edge, b: &TopoBuilder) -> Result<Face, String> {
    let surf = BRepTool::face_surface(face).ok_or("fillet_curved: planar face has no surface")?;
    let mut wires: Vec<Wire> = Vec::new();
    for w in wires_of_face(face) {
        let edges = edges_of_wire(&w);
        let orig_len = edges.len();
        let kept: Vec<Edge> = edges.into_iter().filter(|e| !is_same(&e.0, &fillet_edge.0)).collect();
        if kept.is_empty() {
            continue;
        }
        // A wire may hold the filleted edge next to other boundary edges (e.g. a
        // single-wire annulus); keep the remainder rather than dropping the wire.
        wires.push(if kept.len() < orig_len {
            b.make_wire(&kept)
        } else {
            w
        });
    }
    if wires.is_empty() {
        return Err(
            "fillet_curved: the planar face is bounded only by the filleted edge; its trim would be unbounded"
                .to_string(),
        );
    }
    wires.push(b.make_wire(&[new_circle.clone()]));
    Ok(b.make_face(surf, &wires))
}

/// Rebuild a spherical adjacent face: its new boundary is the single
/// surface-tangency circle.
pub(super) fn rebuild_single_circle_face(face: &Face, new_circle: &Edge, b: &TopoBuilder) -> Result<Face, String> {
    let surf = BRepTool::face_surface(face).ok_or("fillet_curved: face has no surface")?;
    let wire = b.make_wire(&[new_circle.clone()]);
    Ok(b.make_face(surf, &[wire]))
}

/// Rebuild a cylinder lateral face: the base (filleted) circle is replaced by
/// the tangency circle, the seam is re-generated from the tangency circle up to
/// the original top rim, and the top rim circle is reused (shared with the cap).
pub(super) fn rebuild_cylinder_face(
    face: &Face,
    fillet_edge: &Edge,
    contact_edge: &Edge,
    contact_seam: &GpPnt,
    b: &TopoBuilder,
) -> Result<Face, String> {
    let surf = BRepTool::face_surface(face).ok_or("fillet_curved: cylinder face has no surface")?;
    let top = find_top_circle(face, fillet_edge)?;
    let (s_top, _) = BRepTool::edge_vertices(&top).ok_or("fillet_curved: top rim has no vertices")?;
    let seam_up = b.make_edge_segment(contact_seam, &s_top);
    let wire = b.make_wire(&[contact_edge.clone(), seam_up.clone(), top, seam_up]);
    Ok(b.make_face(surf, &[wire]))
}

/// The apex vertex point of a cone lateral face: the face vertex that does not
/// lie on the (filleted) base circle.
pub(super) fn find_cone_apex(face: &Face, fillet_edge: &Edge) -> Result<GpPnt, String> {
    let (se, _) = BRepTool::edge_vertices(fillet_edge)
        .ok_or("fillet_curved: cone base circle has no vertices")?;
    let seam_pt = se;
    for w in wires_of_face(face) {
        for e in edges_of_wire(&w) {
            if is_same(&e.0, &fillet_edge.0) {
                continue;
            }
            if let Some((a, b)) = BRepTool::edge_vertices(&e) {
                if a.distance(&seam_pt) > 1e-6 {
                    return Ok(a);
                }
                if b.distance(&seam_pt) > 1e-6 {
                    return Ok(b);
                }
            }
        }
    }
    Err("fillet_curved: cone lateral face has no apex vertex".to_string())
}
