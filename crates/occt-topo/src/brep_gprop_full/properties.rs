use super::prelude::*;
use super::*;


/// Fill interval bounds from knots (FillIntervalBounds).
pub(super) fn fill_intervals(a: f64, b: f64, knots: &[f64], _num_subs: usize) -> usize {
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
pub(super) fn cinert_perform(curve: &dyn Curve, a: f64, b: f64, loc: &GpPnt) -> GProps {
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
pub(super) enum ValueType {
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
pub(super) struct UFunction<'a> {
    pub(super) fa: &'a FaceGauss,
    pub(super) vertex: GpPnt,
    pub(super) coeffs: &'a [f64],
    pub(super) is_by_point: bool,
    pub(super) v_param: f64,
    pub(super) value_type: ValueType,
}

impl<'a> UFunction<'a> {
    pub(super) fn volume_value(&self, x: f64) -> (f64, GpXyz, f64, f64) {
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

    pub(super) fn value(&self, x: f64) -> f64 {
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

    pub(super) fn inertia_value(&self, f0: f64, pmp0: GpXyz, s: f64, d1: f64) -> f64 {
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
pub(super) fn gk_integrate_arc(
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
pub(super) fn adaptive_integrate<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64, tol: f64) -> Result<f64, String> {
    let (v, _) = adaptive_integrate_with_err(f, a, b, tol)?;
    Ok(v)
}

/// Adaptive GK15 with error estimate.
pub(super) fn adaptive_integrate_with_err<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64, tol: f64) -> Result<(f64, f64), String> {
    pub(super) const NODES: [f64; 15] = [
        -0.9914553711208126, -0.9491079123427585, -0.8648644233597691, -0.7415311855993945,
        -0.5860872354676911, -0.4058451513773972, -0.20778495500789848, 0.0, 0.20778495500789848,
        0.4058451513773972, 0.5860872354676911, 0.7415311855993945, 0.8648644233597691,
        0.9491079123427585, 0.9914553711208126,
    ];
    pub(super) const WK: [f64; 15] = [
        0.022935322010529224, 0.06309209262997856, 0.10479001032225019, 0.14065325971552592,
        0.1690047266392679, 0.19035057806478542, 0.2044329400752989, 0.20948214108472782,
        0.2044329400752989, 0.19035057806478542, 0.1690047266392679, 0.14065325971552592,
        0.10479001032225019, 0.06309209262997856, 0.022935322010529224,
    ];
    pub(super) fn gk15_2<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64) -> (f64, f64) {
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

    pub(super) fn rec<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64, tol: f64, depth: usize) -> Result<(f64, f64), String> {
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

pub(super) fn gk_flags(cg: bool, iflag: bool) -> [bool; NGV] {
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

pub(super) fn value_type_of(k: usize) -> ValueType {
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
pub(super) fn perpendicular(v: &GpVec) -> GpVec {
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
