use super::prelude::*;
use super::*;

impl TopolTool {
    /// Empty tool — no surface. Call [`initialize`](Self::initialize) before
    /// sampling.
    pub fn new() -> Self {
        Self {
            surface: None,
            nb_samples_u: 0,
            nb_samples_v: 0,
            u0: 0.0,
            v0: 0.0,
            du: 1.0,
            dv: 1.0,
            u_range: (0.0, 1.0),
            v_range: (0.0, 1.0),
            u_pars: None,
            v_pars: None,
        }
    }

    /// Whether a surface has been attached.
    pub fn is_initialized(&self) -> bool {
        self.surface.is_some()
    }

    /// Binds `surface` (cloned) and computes the sample grid. Any adaptive
    /// parameters from a previous [`sample_pnts`](Self::sample_pnts) call are
    /// cleared. Source: `IntTools_TopolTool::Initialize`.
    pub fn initialize(&mut self, surface: &dyn Surface) {
        self.surface = Some(surface.clone_dyn());
        self.nb_samples_u = 0;
        self.nb_samples_v = 0;
        self.u0 = 0.0;
        self.v0 = 0.0;
        self.du = 1.0;
        self.dv = 1.0;
        self.u_pars = None;
        self.v_pars = None;
        self.compute_sample_points();
    }

    /// Computes `nb_samples_u`/`nb_samples_v` (and the `du`/`dv` steps) from
    /// the surface type and its parameter domain. Source:
    /// `IntTools_TopolTool::ComputeSamplePoints`.
    ///
    /// `ponytail:` analytic type is recovered via [`classify_surface`]; radii
    /// for the angular-step count are measured by sampling the surface, and
    /// BSpline/Bezier pole/knot counts (not reachable through `dyn Surface`)
    /// fall back to a 10×10 base grid.
    pub fn compute_sample_points(&mut self) {
        let Some(surface) = self.surface.as_ref() else { return };
        let s: &dyn Surface = surface.as_ref();
        let (mut uinf, mut usup) = s.u_range();
        let (mut vinf, mut vsup) = s.v_range();

        if usup < uinf {
            std::mem::swap(&mut uinf, &mut usup);
        }
        if vsup < vinf {
            std::mem::swap(&mut vinf, &mut vsup);
        }

        // Clamp unbounded directions to a big-but-finite span (OCCT sentinel).
        let is_big_uinf = !uinf.is_finite() && uinf < 0.0;
        let is_big_usup = !usup.is_finite() && usup > 0.0;
        let is_big_vinf = !vinf.is_finite() && vinf < 0.0;
        let is_big_vsup = !vsup.is_finite() && vsup > 0.0;
        if is_big_uinf && is_big_usup {
            uinf = -BIG_RANGE;
            usup = BIG_RANGE;
        } else if is_big_uinf {
            uinf = usup - 2.0 * BIG_RANGE;
        } else if is_big_usup {
            usup = uinf + 2.0 * BIG_RANGE;
        }
        if is_big_vinf && is_big_vsup {
            vinf = -BIG_RANGE;
            vsup = BIG_RANGE;
        } else if is_big_vinf {
            vinf = vsup - 2.0 * BIG_RANGE;
        } else if is_big_vsup {
            vsup = vinf + 2.0 * BIG_RANGE;
        }

        self.u0 = uinf;
        self.v0 = vinf;
        self.u_range = (uinf, usup);
        self.v_range = (vinf, vsup);
        self.u_pars = None;
        self.v_pars = None;

        let typ = classify_surface(s);
        let (mut nbsu, mut nbsv): (usize, usize) = match typ {
            SurfaceType::Plane => (10, 10),
            SurfaceType::Cylinder => {
                let radius = cylinder_radius(s);
                let max_angle = max_angle_for_radius(radius);
                let nbsu = if max_angle > ANGULAR {
                    ((usup - uinf) / max_angle) as usize
                } else {
                    0
                };
                let nbsv = ((vsup - vinf) / 10.0) as usize;
                (nbsu.max(2).min(MAX_NB_SAMPLE), nbsv.max(2).min(MAX_NB_SAMPLE))
            }
            SurfaceType::Cone => {
                let radius = cone_radius_at(s, vinf).max(cone_radius_at(s, vsup));
                let max_angle = max_angle_for_radius(radius);
                let nbsu = if max_angle > ANGULAR {
                    ((usup - uinf) / max_angle) as usize
                } else {
                    0
                };
                let nbsv = ((vsup - vinf) / 10.0) as usize;
                (nbsu.max(10).min(MAX_NB_SAMPLE), nbsv.max(10).min(MAX_NB_SAMPLE))
            }
            SurfaceType::Sphere => {
                let radius = sphere_radius(s);
                let max_angle = max_angle_for_radius(radius);
                let nbsu = if max_angle > ANGULAR {
                    ((usup - uinf) / max_angle) as usize
                } else {
                    0
                };
                let nbsv = if max_angle > ANGULAR {
                    ((vsup - vinf) / max_angle) as usize
                } else {
                    0
                };
                (nbsu.max(10).min(MAX_NB_SAMPLE), nbsv.max(10).min(MAX_NB_SAMPLE))
            }
            SurfaceType::Torus => {
                let (_major, minor) = torus_radii(s);
                let max_angle = max_angle_for_radius(minor);
                let nbsu = if max_angle > ANGULAR {
                    ((usup - uinf) / max_angle) as usize
                } else {
                    0
                };
                let nbsv = if max_angle > ANGULAR {
                    ((vsup - vinf) / max_angle) as usize
                } else {
                    0
                };
                (nbsu.max(10).min(MAX_NB_SAMPLE), nbsv.max(10).min(MAX_NB_SAMPLE))
            }
            SurfaceType::BezierSurface | SurfaceType::BSplineSurface => {
                // Pole/knot counts are unreachable through `dyn Surface`; the
                // base grid is refined adaptively by `sample_pnts`.
                (10, 10)
            }
            SurfaceType::SurfaceOfExtrusion => {
                let nbsv = ((vsup - vinf) / 10.0) as usize;
                (15, nbsv.max(15).min(MAX_NB_SAMPLE))
            }
            SurfaceType::SurfaceOfRevolution => (15, 15),
            SurfaceType::OffsetSurface | SurfaceType::OtherSurface => (10, 10),
        };
        if nbsu == 0 {
            nbsu = 10;
        }
        if nbsv == 0 {
            nbsv = 10;
        }

        self.nb_samples_u = nbsu;
        self.nb_samples_v = nbsv;
        self.du = (usup - uinf) / (nbsu + 1) as f64;
        self.dv = (vsup - vinf) / (nbsv + 1) as f64;
    }

    /// Number of sample points along U. Returns 0 when no surface is attached.
    pub fn nb_samples_u(&self) -> usize {
        self.nb_samples_u
    }

    /// Number of sample points along V. Returns 0 when no surface is attached.
    pub fn nb_samples_v(&self) -> usize {
        self.nb_samples_v
    }

    /// Total number of sample points: `nb_samples_u * nb_samples_v`.
    pub fn nb_samples(&self) -> usize {
        self.nb_samples_u * self.nb_samples_v
    }

    /// The clamped `[u_first, u_last]` parameter domain.
    pub fn u_range(&self) -> (f64, f64) {
        self.u_range
    }

    /// The clamped `[v_first, v_last]` parameter domain.
    pub fn v_range(&self) -> (f64, f64) {
        self.v_range
    }

    /// The U step between interior samples.
    pub fn u_step(&self) -> f64 {
        self.du
    }

    /// The V step between interior samples.
    pub fn v_step(&self) -> f64 {
        self.dv
    }

    /// Returns the `index`-th sample: the `(u, v)` parameter pair plus the 3D
    /// surface point `surface.d0(u, v)`.
    ///
    /// `index` is 1-based, from `1` to [`nb_samples`](Self::nb_samples), in
    /// row-major order (`U` fastest). When [`sample_pnts`](Self::sample_pnts)
    /// has been called, the adaptive parameter grid is used instead of the
    /// uniform one. Source: `IntTools_TopolTool::SamplePoint`.
    pub fn sample_point(&self, index: usize) -> Result<(GpPnt2d, GpPnt), String> {
        let s = self
            .surface
            .as_ref()
            .ok_or("TopolTool::sample_point: no surface initialized")?;
        let s: &dyn Surface = s.as_ref();
        if self.nb_samples_u == 0 || self.nb_samples_v == 0 {
            return Err("TopolTool::sample_point: sample grid not computed".into());
        }
        let (nu, nv) = match (&self.u_pars, &self.v_pars) {
            (Some(u), Some(v)) => (u.len(), v.len()),
            _ => (self.nb_samples_u, self.nb_samples_v),
        };
        if index == 0 || index > nu * nv {
            return Err(format!(
                "TopolTool::sample_point: index {index} out of range [1, {}]",
                nu * nv
            ));
        }
        let (iu, iv) = ((index - 1) % nu, (index - 1) / nu);
        let (u, v) = match (&self.u_pars, &self.v_pars) {
            (Some(u), Some(v)) => (u[iu], v[iv]),
            _ => (
                self.u0 + (iu + 1) as f64 * self.du,
                self.v0 + (iv + 1) as f64 * self.dv,
            ),
        };
        let p3d = s.d0(u, v);
        Ok((GpPnt2d::new(u, v), p3d))
    }

    /// Builds a full grid of sample points.
    ///
    /// For BSpline/Bezier surfaces the `U` and `V` parameter sequences are
    /// refined adaptively until the chord deviation of each interval stays
    /// under `deflection`; every other surface uses a uniform grid of
    /// `max(nb_samples, nu_min) × max(nb_samples, nv_min)` interior samples.
    /// The `u_pars`/`v_pars` are stored so subsequent
    /// [`sample_point`](Self::sample_point) calls walk the refined grid.
    ///
    /// Returns `(u, v, surface.d0(u, v))` for every parameter pair, in
    /// row-major order. Source: `IntTools_TopolTool::SamplePnts` +
    /// `Adaptor3d_TopolTool::SamplePnts`.
    pub fn sample_pnts(
        &mut self,
        deflection: f64,
        nu_min: usize,
        nv_min: usize,
    ) -> Result<Vec<(GpPnt2d, GpPnt)>, String> {
        if self.surface.is_none() {
            return Err("TopolTool::sample_pnts: no surface initialized".into());
        }
        // Recompute the analytic base grid (idempotent; also clears any stale
        // adaptive parameters from a previous call).
        self.compute_sample_points();
        let s = self.surface.as_ref().expect("checked above");
        let s: &dyn Surface = s.as_ref();
        let (u0, u1) = self.u_range;
        let (v0, v1) = self.v_range;
        let nbsu = self.nb_samples_u.max(nu_min).max(1);
        let nbsv = self.nb_samples_v.max(nv_min).max(1);

        let uniform_u = || -> Vec<f64> {
            let du = (u1 - u0) / (nbsu + 1) as f64;
            (1..=nbsu).map(|i| u0 + i as f64 * du).collect()
        };
        let uniform_v = || -> Vec<f64> {
            let dv = (v1 - v0) / (nbsv + 1) as f64;
            (1..=nbsv).map(|i| v0 + i as f64 * dv).collect()
        };

        let typ = classify_surface(s);
        let (u_pars, v_pars) = if matches!(
            typ,
            SurfaceType::BSplineSurface | SurfaceType::BezierSurface
        ) {
            let u_base = uniform_u();
            let v_base = uniform_v();
            let v_fixed = 0.5 * (v0 + v1);
            let u_fixed = 0.5 * (u0 + u1);
            (
                refine_params(s, deflection, u_base, v_fixed, true, (u0, u1)),
                refine_params(s, deflection, v_base, u_fixed, false, (v0, v1)),
            )
        } else {
            (uniform_u(), uniform_v())
        };

        self.u_pars = Some(u_pars.clone());
        self.v_pars = Some(v_pars.clone());
        self.nb_samples_u = u_pars.len();
        self.nb_samples_v = v_pars.len();

        let mut out = Vec::with_capacity(u_pars.len() * v_pars.len());
        for &u in &u_pars {
            for &v in &v_pars {
                out.push((GpPnt2d::new(u, v), s.d0(u, v)));
            }
        }
        Ok(out)
    }
}

impl Default for TopolTool {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Surface-classification helpers (measured through `dyn Surface`)
// ---------------------------------------------------------------------------

/// Angular step for a circle of `radius` under the 1e-2 sampling deflection.
/// Source: `IntTools_TopolTool::ComputeSamplePoints`
/// (`acos(1 - deflection / radius) * 2`, floored by `π/2`).
pub(super) fn max_angle_for_radius(radius: f64) -> f64 {
    let mut max_angle = std::f64::consts::PI * 0.5;
    if radius > SAMPLE_DEFLECTION {
        max_angle = (1.0 - SAMPLE_DEFLECTION / radius).acos() * 2.0;
    }
    max_angle
}

/// Radius of a cylinder measured from the surface: half the distance between
/// the diametrically opposite points `(0, v)` and `(π, v)`.
pub(super) fn cylinder_radius(s: &dyn Surface) -> f64 {
    s.d0(0.0, 0.0).distance(&s.d0(std::f64::consts::PI, 0.0)) * 0.5
}

/// Radius of a sphere measured from the surface.
pub(super) fn sphere_radius(s: &dyn Surface) -> f64 {
    s.d0(0.0, 0.0).distance(&s.d0(std::f64::consts::PI, 0.0)) * 0.5
}

/// Radius of the U-circle of a cone at parameter `v`, measured from the
/// surface.
pub(super) fn cone_radius_at(s: &dyn Surface, v: f64) -> f64 {
    s.d0(0.0, v).distance(&s.d0(std::f64::consts::PI, v)) * 0.5
}

/// `(major, minor)` radii of a torus measured from the surface.
pub(super) fn torus_radii(s: &dyn Surface) -> (f64, f64) {
    let minor = s.d0(0.0, 0.0).distance(&s.d0(0.0, std::f64::consts::PI)) * 0.5;
    let major_plus_minor = s.d0(0.0, 0.0).distance(&s.d0(std::f64::consts::PI, 0.0)) * 0.5;
    (major_plus_minor - minor, minor)
}

/// Chord deviation of the surface point at parameter `m` against the linear
/// interpolation of the points at `a` and `b`, along the `is_u` iso-line
/// (the orthogonal parameter held at `fixed`).
pub(super) fn slice_deviation(s: &dyn Surface, a: f64, b: f64, m: f64, fixed: f64, is_u: bool) -> f64 {
    let pa = if is_u { s.d0(a, fixed) } else { s.d0(fixed, a) };
    let pb = if is_u { s.d0(b, fixed) } else { s.d0(fixed, b) };
    let pm = if is_u { s.d0(m, fixed) } else { s.d0(fixed, m) };
    let t = (m - a) / (b - a);
    let interp = GpPnt::new(
        pa.x() + t * (pb.x() - pa.x()),
        pa.y() + t * (pb.y() - pa.y()),
        pa.z() + t * (pb.z() - pa.z()),
    );
    pm.distance(&interp)
}

/// Deflection-adaptive 1D refinement of a parameter sequence.
///
/// Starts from the interior `base` samples plus the two domain endpoints;
/// every interval whose midpoint deviates from its chord by more than
/// `deflection` is split at the midpoint. Repeats until no interval needs
/// splitting (bounded by a safety cap). The returned sequence drops the domain
/// endpoints, keeping the interior-only semantics of the uniform grid (and of
/// OCCT's `SamplePoint` `iu`/`iv` in `1..=nb`).
pub(super) fn refine_params(
    s: &dyn Surface,
    deflection: f64,
    base: Vec<f64>,
    fixed: f64,
    is_u: bool,
    domain: (f64, f64),
) -> Vec<f64> {
    pub(super) const MAX_PARAMS: usize = 4096;
    let deflection = deflection.max(1e-9);
    let mut params: Vec<f64> = base;
    params.push(domain.0);
    params.push(domain.1);
    params.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    params.dedup_by(|a, b| (*a - *b).abs() < 1e-12);

    loop {
        let mut changed = false;
        let mut next: Vec<f64> = Vec::with_capacity(params.len().min(MAX_PARAMS) + 16);
        for w in params.windows(2) {
            let (a, b) = (w[0], w[1]);
            if b - a <= 1e-12 {
                next.push(a);
                continue;
            }
            let m = 0.5 * (a + b);
            let dev = slice_deviation(s, a, b, m, fixed, is_u);
            if dev > deflection && next.len() < MAX_PARAMS {
                next.push(a);
                next.push(m);
                changed = true;
            } else {
                next.push(a);
            }
        }
        if let Some(&last) = params.last() {
            next.push(last);
        }
        next.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        next.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
        params = next;
        if !changed || params.len() >= MAX_PARAMS {
            break;
        }
    }

    params
        .into_iter()
        .filter(|&p| p > domain.0 + 1e-9 && p < domain.1 - 1e-9)
        .collect()
}
