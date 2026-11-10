//! Periodic-domain checks, `ParametricResolution` and `U/VResolution`.

use super::*;

pub(crate) fn check_reparam(
    flat: &[f64],
    poles: &[f64],
    nbp: usize,
    pc3d: &[f64],
    c3d: &dyn Curve,
    c2d: &dyn Curve2d,
    s: &dyn Surface,
    tol: &mut f64,
    oldtol: f64,
) -> bool {
    let mut param_first = 3.0 * pc3d[0] - 2.0 * pc3d[nbp - 1];
    let mut param_last = 3.0 * pc3d[nbp - 1] - 2.0 * pc3d[0];
    let first_par = c2d.first_parameter();
    let last_par = c2d.last_parameter();
    if param_first < first_par {
        param_first = first_par;
    }
    if param_last > last_par {
        param_last = last_par;
    }
    let d = *tol;
    let nn = 2 * nbp;
    let unsurnn = 1.0 / nn as f64;
    let mut tprev = param_first;
    let mut d2: f64 = 0.0;
    for i in 0..=nn {
        let t = unsurnn * i as f64;
        let tc3d = pc3d[0] * (1.0 - t) + pc3d[nbp - 1] * t;
        let pc3d_pt = c3d.d0(tc3d);
        let tcons = eval_1d_d0(flat, poles, tc3d);
        if tcons < tprev || tcons > param_last {
            *tol = INFINITE;
            return false;
        }
        tprev = tcons;
        let pcons = cons_value(c2d, s, tcons);
        d2 = d2.max(pc3d_pt.square_distance(&pcons));
    }
    *tol = d2.sqrt();
    for i in 1..poles.len() {
        if poles[i - 1] > poles[i] {
            return false;
        }
    }
    *tol <= d || *tol > 0.8 * oldtol
}

/// `Precision::Parametric(R3d)` (`Precision.hxx:328`, `= R3d * 0.01`), the
/// `default` arm of `GeomAdaptor_Surface::UResolution` / `VResolution`.
#[inline]
pub(crate) fn parametric_resolution(r3d: f64) -> f64 {
    r3d * PCONFUSION / CONFUSION
}

/// `GeomAdaptor_Surface::UResolution` (`GeomAdaptor_Surface.cxx:1818-1896`).
/// `r3d` is OCCT's `R3d` argument.
///
/// Shared by every OCCT call site that asks a surface for a parametric
/// resolution: `Approx_SameParameter.cxx:429` / `:504`
/// (`mySurf->UResolution(besttol)`), `BRepLib.cxx:1085`
/// (`DSdu = 1./surf->UResolution(1.)`) and `IntTools_WLineTool.cxx:153`
/// (`aGAS1.UResolution(aDelta)`).
pub fn u_resolution(s: &dyn Surface, r3d: f64) -> f64 {
    // `GeomAdaptor_Surface::load` peels a `Geom_RectangularTrimmedSurface`
    // down to its basis surface (`cxx:423-430`), so the switch below sees
    // the basis type. NOTE: OCCT keeps the trimmed range as `myUFirst..myVLast`;
    // this recursion uses the basis range, which only matters for the
    // unbounded-V Cone guard (`cxx:1855-1858`).
    if let Some(b) = s.rectangular_trimmed_basis() {
        return u_resolution(b.as_ref(), r3d);
    }
    // `case GeomAbs_SurfaceOfExtrusion: BasisCurve->Resolution(R3d)`
    // (`cxx:1824-1827`), live through `GeomSurfaceOfLinearExtrusion` and the
    // trimmed/offset wrappers that forward `extrusion_basis_curve`.
    if let Some(c) = s.extrusion_basis_curve() {
        return c.resolution(r3d);
    }
    // `case GeomAbs_OffsetSurface: BasisAdaptor->UResolution(R3d)`
    // (`cxx:1882-1885`).
    if let Some(b) = s.offset_basis_surface() {
        return u_resolution(b.as_ref(), r3d);
    }
    // `case GeomAbs_Cone` (`cxx:1853-1866`). The guard is the infinite-domain
    // one: `myVLast - myVFirst > 1.e10` means "not truly bounded".
    if let Some((radius, semi_angle)) = s.cone_ref() {
        let (v_first, v_last) = s.v_range();
        if v_last - v_first > 1.0e10 {
            return parametric_resolution(r3d);
        }
        // `S->VIso(V)` radius is `Radius + V*sin(Angle)`, taken absolute
        // (`ElSLib.cxx:1793-1811`).
        let rayon1 = (radius + v_last * semi_angle.sin()).abs();
        let rayon2 = (radius + v_first * semi_angle.sin()).abs();
        let r = if rayon1 > rayon2 { rayon1 } else { rayon2 };
        return if r > CONFUSION { r3d / r } else { 0.0 };
    }
    // `case GeomAbs_Plane: return R3d` (`cxx:1867-1869`).
    if s.gp_pln().is_some() {
        return r3d;
    }
    // `GeomAbs_Torus` / `GeomAbs_Sphere` / `GeomAbs_Cylinder` compute
    // `Res = R3d / (2. * R)` and then fall into the `2.*asin(Res)` tail
    // (`cxx:1828-1852`, `:1890-1895`).
    let res = if let Some(t) = s.gp_torus() {
        let r = t.major_radius + t.minor_radius;
        if r > CONFUSION { r3d / (2.0 * r) } else { 0.0 }
    } else if let Some(sp) = s.gp_sphere() {
        if sp.radius > CONFUSION { r3d / (2.0 * sp.radius) } else { 0.0 }
    } else if let Some(cy) = s.gp_cylinder() {
        if cy.radius > CONFUSION { r3d / (2.0 * cy.radius) } else { 0.0 }
    } else if let Some((ures, _)) = s.uv_resolution(r3d) {
        // `GeomAbs_BSplineSurface` (`cxx:1875-1881`) and the offset/revolution
        // wrappers. UNPORTED: `Geom_BezierSurface::Resolution`
        // (`GeomAdaptor_Surface.cxx:1870-1874`) has no Rust implementation, so
        // a Bezier surface falls through to the `default` arm below.
        return ures;
    } else {
        // `default: return Precision::Parametric(R3d)` (`cxx:1886-1887`).
        return parametric_resolution(r3d);
    };
    if res <= 1.0 {
        2.0 * res.asin()
    } else {
        2.0 * std::f64::consts::PI
    }
}

/// `GeomAdaptor_Surface::VResolution` (`GeomAdaptor_Surface.cxx:1900-1958`).
/// `r3d` is OCCT's `R3d` argument.
pub fn v_resolution(s: &dyn Surface, r3d: f64) -> f64 {
    // Trimmed-surface unwrap, as in `u_resolution` above.
    if let Some(b) = s.rectangular_trimmed_basis() {
        return v_resolution(b.as_ref(), r3d);
    }
    // `case GeomAbs_SurfaceOfRevolution: BasisCurve->Resolution(R3d)`
    // (`cxx:1905-1908`).
    if let Some(c) = s.revolution_basis_curve() {
        return c.resolution(r3d);
    }
    // `case GeomAbs_OffsetSurface: BasisAdaptor->VResolution(R3d)`
    // (`cxx:1944-1947`).
    if let Some(b) = s.offset_basis_surface() {
        return v_resolution(b.as_ref(), r3d);
    }
    // `GeomAbs_Torus` (MinorRadius, `cxx:1910-1918`) and `GeomAbs_Sphere`
    // (`cxx:1919-1927`) share the `2.*asin(Res)` tail (`cxx:1952-1957`).
    let res = if let Some(t) = s.gp_torus() {
        let r = t.minor_radius;
        if r > CONFUSION { r3d / (2.0 * r) } else { 0.0 }
    } else if let Some(sp) = s.gp_sphere() {
        if sp.radius > CONFUSION { r3d / (2.0 * sp.radius) } else { 0.0 }
    } else if s.gp_pln().is_some()
        || s.gp_cylinder().is_some()
        || s.cone_ref().is_some()
        || s.extrusion_basis_curve().is_some()
    {
        // `case GeomAbs_SurfaceOfExtrusion / GeomAbs_Cylinder / GeomAbs_Cone /
        // GeomAbs_Plane: return R3d` (`cxx:1928-1933`).
        return r3d;
    } else if let Some((_, vres)) = s.uv_resolution(r3d) {
        return vres;
    } else {
        // `default: return Precision::Parametric(R3d)` (`cxx:1948-1949`).
        return parametric_resolution(r3d);
    };
    if res <= 1.0 {
        2.0 * res.asin()
    } else {
        2.0 * std::f64::consts::PI
    }
}

pub(crate) fn curve2d_from_approx(approx: &ApproxAFunction1dPair) -> Option<Arc<dyn Curve2d>> {
    if !approx.has_result || approx.poles_u.is_empty() {
        return None;
    }
    let deg = approx.degree.max(1) as usize;
    let flat = knot_sequence(&approx.knots, &approx.mults, approx.degree);
    let c = Geom2dBSplineCurve::new(approx.poles_u.clone(), approx.poles_v.clone(), flat, deg).ok()?;
    Some(Arc::new(c))
}

