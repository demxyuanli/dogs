use super::prelude::*;
use super::*;

impl DefaultRangeSplitter {
    /// Creates an empty splitter.
    pub fn new() -> Self {
        Self {
            dface: None,
            surface: None,
            deflection: 0.0,
            range_u: (1e100, -1e100),
            range_v: (1e100, -1e100),
            delta: (1.0, 1.0),
            tolerance: (CONFUSION, CONFUSION),
            is_valid: true,
        }
    }

    /// Resets the splitter. Source: `Reset`.
    pub fn reset(&mut self, dface: &MeshFace, _params: &MeshParameters) {
        self.dface = Some(dface.clone());
        self.surface = dface.surface();
        self.deflection = dface.deflection();
        self.range_u = (1e100, -1e100);
        self.range_v = (1e100, -1e100);
        self.delta = (1.0, 1.0);
        self.tolerance = (CONFUSION, CONFUSION);
        self.is_valid = true;
    }

    /// Registers a boundary point. Source: `AddPoint`.
    pub fn add_point(&mut self, point: GpPnt2d) {
        self.range_u.0 = self.range_u.0.min(point.x());
        self.range_u.1 = self.range_u.1.max(point.x());
        self.range_v.0 = self.range_v.0.min(point.y());
        self.range_v.1 = self.range_v.1.max(point.y());
    }

    /// True when the computed range is valid. Source: `IsValid`.
    pub fn is_valid(&self) -> bool {
        self.is_valid
    }

    /// Scales a point between real parametric space and the face basis.
    /// Source: `Scale`.
    pub fn scale(&self, point: GpPnt2d, to_face_basis: bool) -> GpPnt2d {
        if to_face_basis {
            GpPnt2d::new(
                (point.x() - self.range_u.0) / self.delta.0,
                (point.y() - self.range_v.0) / self.delta.1,
            )
        } else {
            GpPnt2d::new(
                point.x() * self.delta.0 + self.range_u.0,
                point.y() * self.delta.1 + self.range_v.0,
            )
        }
    }

    /// The base splitter generates no interior nodes (null list in OCCT).
    pub fn generate_surface_nodes(&self, _params: &MeshParameters) -> Option<Vec<GpPnt2d>> {
        None
    }

    /// 3D point at the given UV parameter. Source: `Point`.
    pub fn point(&self, point2d: GpPnt2d) -> GpPnt {
        self.surface
            .as_ref()
            .map(|s| s.d0(point2d.x(), point2d.y()))
            .unwrap_or_else(GpPnt::zero)
    }

    /// Discrete U range.
    pub fn range_u(&self) -> (f64, f64) {
        self.range_u
    }
    /// Discrete V range.
    pub fn range_v(&self) -> (f64, f64) {
        self.range_v
    }
    /// Scale factors per direction.
    pub fn delta(&self) -> (f64, f64) {
        self.delta
    }
    /// Parametric tolerances per direction.
    pub fn tolerance_uv(&self) -> (f64, f64) {
        self.tolerance
    }
    /// The discrete face model.
    pub fn dface(&self) -> Option<&MeshFace> {
        self.dface.as_ref()
    }
    /// The face surface.
    pub fn surface(&self) -> Option<&Arc<dyn Surface>> {
        self.surface.as_ref()
    }
    /// The face deflection.
    pub fn deflection(&self) -> f64 {
        self.deflection
    }

    /// Length along U of the discrete range, sampled on a 20-segment grid.
    /// Source: `computeLengthU`.
    pub fn compute_length_u(&self, s: &dyn Surface) -> f64 {
        let (u0, u1) = self.range_u;
        let (v0, v1) = self.range_v;
        let mut long = 0.0;
        let du = 0.05 * (u1 - u0);
        let v_ave = 0.5 * (v1 + v0);
        let mut p11 = s.d0(u0, v0);
        let mut p21 = s.d0(u0, v_ave);
        let mut p31 = s.d0(u0, v1);
        let mut u = u0 + du;
        for _ in 1..=20 {
            let p12 = s.d0(u, v0);
            let p22 = s.d0(u, v_ave);
            let p32 = s.d0(u, v1);
            long += p11.distance(&p12) + p21.distance(&p22) + p31.distance(&p32);
            p11 = p12;
            p21 = p22;
            p31 = p32;
            u += du;
        }
        long / 3.0
    }

    /// Length along V of the discrete range, sampled on a 20-segment grid.
    /// Source: `computeLengthV`.
    pub fn compute_length_v(&self, s: &dyn Surface) -> f64 {
        let (u0, u1) = self.range_u;
        let (v0, v1) = self.range_v;
        let mut long = 0.0;
        let dv = 0.05 * (v1 - v0);
        let u_ave = 0.5 * (u1 + u0);
        let mut p11 = s.d0(u0, v0);
        let mut p21 = s.d0(u_ave, v0);
        let mut p31 = s.d0(u1, v0);
        let mut v = v0 + dv;
        for _ in 1..=20 {
            let p12 = s.d0(u0, v);
            let p22 = s.d0(u_ave, v);
            let p32 = s.d0(u1, v);
            long += p11.distance(&p12) + p21.distance(&p22) + p31.distance(&p32);
            p11 = p12;
            p21 = p22;
            p31 = p32;
            v += dv;
        }
        long / 3.0
    }

    /// Computes the parametric tolerances. Source: `computeTolerance`.
    pub fn compute_tolerance(&mut self, _len_u: f64, _len_v: f64) {
        let diff_u = self.range_u.1 - self.range_u.0;
        let diff_v = self.range_v.1 - self.range_v.0;
        let face_tol = self
            .dface
            .as_ref()
            .map(|f| BRepTool::face_tolerance(f.face()))
            .unwrap_or(CONFUSION);
        let res_u = self
            .surface
            .as_ref()
            .map(|s| param_resolution(s.as_ref(), face_tol, true) * 1.1)
            .unwrap_or(face_tol);
        let res_v = self
            .surface
            .as_ref()
            .map(|s| param_resolution(s.as_ref(), face_tol, false) * 1.1)
            .unwrap_or(face_tol);
        pub(super) const DEFLECTION_UV: f64 = 1e-5;
        self.tolerance.0 = (DEFLECTION_UV.min(res_u)).max(1e-7 * diff_u);
        self.tolerance.1 = (DEFLECTION_UV.min(res_v)).max(1e-7 * diff_v);
    }

    /// Computes the scale factors. Source: `computeDelta`.
    pub fn compute_delta(&mut self, len_u: f64, len_v: f64) {
        let diff_u = self.range_u.1 - self.range_u.0;
        let diff_v = self.range_v.1 - self.range_v.0;
        self.delta.0 = diff_u / (if len_u < self.tolerance.0 { 1.0 } else { len_u });
        self.delta.1 = diff_v / (if len_v < self.tolerance.1 { 1.0 } else { len_v });
    }

    pub(crate) fn set_range_u(&mut self, r: (f64, f64)) {
        self.range_u = r;
    }
    pub(crate) fn set_range_v(&mut self, r: (f64, f64)) {
        self.range_v = r;
    }
    pub(crate) fn set_delta(&mut self, d: (f64, f64)) {
        self.delta = d;
    }
    pub(crate) fn set_valid(&mut self, v: bool) {
        self.is_valid = v;
    }
}

impl RangeSplitter for DefaultRangeSplitter {
    fn base(&self) -> &DefaultRangeSplitter {
        self
    }
    fn base_mut(&mut self) -> &mut DefaultRangeSplitter {
        self
    }
}

/// Cylindrical surface splitter — U-periodic seam, interior nodes along the
/// parametric grid. Source: `BRepMesh_CylinderRangeSplitter`.
pub struct CylinderRangeSplitter {
    pub(super) inner: DefaultRangeSplitter,
    pub(super) du: f64,
}

impl CylinderRangeSplitter {
    /// Creates an empty splitter.
    pub fn new() -> Self {
        Self {
            inner: DefaultRangeSplitter::new(),
            du: 1.0,
        }
    }
}

impl RangeSplitter for CylinderRangeSplitter {
    fn base(&self) -> &DefaultRangeSplitter {
        &self.inner
    }
    fn base_mut(&mut self) -> &mut DefaultRangeSplitter {
        &mut self.inner
    }

    fn reset(&mut self, dface: &MeshFace, params: &MeshParameters) {
        self.reset_base(dface, params);
        let r = self.surface().map(|s| cylinder_radius(s.as_ref())).unwrap_or(1.0);
        let defl = self.deflection();
        self.du = arc_angular_step(r, defl, params.angle, params.min_size);
    }

    fn compute_delta(&mut self, _len_u: f64, len_v: f64) {
        let range_v = self.base().range_v();
        self.inner
            .set_delta((self.du / len_v.max(range_v.1 - range_v.0), 1.0));
    }

    fn generate_surface_nodes(&self, _params: &MeshParameters) -> Option<Vec<GpPnt2d>> {
        let range_u = self.base().range_u();
        let range_v = self.base().range_v();
        let radius = self.surface().map(|s| cylinder_radius(s.as_ref())).unwrap_or(0.0);
        let deflection = self.deflection();

        let su = range_u.1 - range_u.0;
        let sv = range_v.1 - range_v.0;
        let a_arc_len = su * radius;
        let mut nb_u = 0i32;
        let mut nb_v = 0i32;
        if a_arc_len > deflection {
            nb_u = (su / self.du) as i32;
            // ponytail: the OCCT V-step computation is commented out, so nbV stays 0
            // and no interior rows are produced.
        }
        let du = su / (nb_u + 1) as f64;
        let dv = sv / (nb_v + 1) as f64;

        let pas_max_v = range_v.1 - dv * 0.5;
        let pas_max_u = range_u.1 - du * 0.5;
        let mut nodes = Vec::new();
        let mut pas_v = range_v.0 + dv;
        while pas_v < pas_max_v {
            let mut pas_u = range_u.0 + du;
            while pas_u < pas_max_u {
                nodes.push(GpPnt2d::new(pas_u, pas_v));
                pas_u += du;
            }
            pas_v += dv;
        }
        Some(nodes)
    }
}

/// Conical surface splitter. Source: `BRepMesh_ConeRangeSplitter`.
pub struct ConeRangeSplitter {
    pub(super) inner: DefaultRangeSplitter,
}

impl ConeRangeSplitter {
    /// Creates an empty splitter.
    pub fn new() -> Self {
        Self {
            inner: DefaultRangeSplitter::new(),
        }
    }

    /// Returns the split steps along U and V and the number of steps.
    /// Source: `GetSplitSteps`.
    pub fn get_split_steps(
        &self,
        params: &MeshParameters,
        steps_nb: &mut (i32, i32),
    ) -> (f64, f64) {
        let range_u = self.base().range_u();
        let range_v = self.base().range_v();
        let surface = self.surface().unwrap();
        let deflection = self.deflection();
        // `BRepMesh_ConeRangeSplitter::GetSplitSteps`:
        // `aRadius = max(|RefR + V0*sin(ang)|, |RefR + V1*sin(ang)|)`.
        let a_radius = if let Some((ref_r, sang)) = surface.cone_ref() {
            (ref_r + range_v.0 * sang.sin())
                .abs()
                .max((ref_r + range_v.1 * sang.sin()).abs())
        } else {
            cone_radius_at(surface.as_ref(), range_v.0)
                .max(cone_radius_at(surface.as_ref(), range_v.1))
        };

        let mut du = arc_angular_step(a_radius, deflection, params.angle, params.min_size);

        let a_diff_u = range_u.1 - range_u.0;
        let a_diff_v = range_v.1 - range_v.0;
        let a_scale = du * a_radius;
        let a_ratio = (a_diff_v / a_scale).ln().max(1.0);
        let nb_u = (a_diff_u / du) as i32;
        let nb_v = (a_diff_v / a_scale / a_ratio) as i32;

        du = a_diff_u / (nb_u + 1) as f64;
        let dv = a_diff_v / (nb_v + a_ratio as i32) as f64;

        steps_nb.0 = nb_u;
        steps_nb.1 = nb_v;
        (du, dv)
    }
}

impl RangeSplitter for ConeRangeSplitter {
    fn base(&self) -> &DefaultRangeSplitter {
        &self.inner
    }
    fn base_mut(&mut self) -> &mut DefaultRangeSplitter {
        &mut self.inner
    }

    fn generate_surface_nodes(&self, params: &MeshParameters) -> Option<Vec<GpPnt2d>> {
        let range_u = self.base().range_u();
        let range_v = self.base().range_v();
        let mut steps_nb = (0i32, 0i32);
        let (du, dv) = self.get_split_steps(params, &mut steps_nb);

        let pas_max_v = range_v.1 - dv * 0.5;
        let pas_max_u = range_u.1 - du * 0.5;
        let mut nodes = Vec::new();
        let mut pas_v = range_v.0 + dv;
        while pas_v < pas_max_v {
            let mut pas_u = range_u.0 + du;
            while pas_u < pas_max_u {
                nodes.push(GpPnt2d::new(pas_u, pas_v));
                pas_u += du;
            }
            pas_v += dv;
        }
        Some(nodes)
    }
}

/// Spherical surface splitter — staggered U/V grid. Source:
/// `BRepMesh_SphereRangeSplitter`.
pub struct SphereRangeSplitter {
    pub(super) inner: DefaultRangeSplitter,
}

impl SphereRangeSplitter {
    /// Creates an empty splitter.
    pub fn new() -> Self {
        Self {
            inner: DefaultRangeSplitter::new(),
        }
    }

    /// Computes the step and upper bound for a range. Source: `computeStep`.
    pub(super) fn compute_step(&self, range: (f64, f64), default_step: f64) -> (f64, f64) {
        let diff = range.1 - range.0;
        let step = diff / ((diff / default_step) as i32 + 1) as f64;
        (step, range.1 - PCONFUSION)
    }
}

impl RangeSplitter for SphereRangeSplitter {
    fn base(&self) -> &DefaultRangeSplitter {
        &self.inner
    }
    fn base_mut(&mut self) -> &mut DefaultRangeSplitter {
        &mut self.inner
    }

    fn generate_surface_nodes(&self, params: &MeshParameters) -> Option<Vec<GpPnt2d>> {
        let range_v = self.base().range_v();
        let range_u = self.base().range_u();
        let radius = self.surface().map(|s| sphere_radius(s.as_ref())).unwrap_or(0.0);
        let deflection = self.deflection();
        let a_step = 0.7 * arc_angular_step(radius, deflection, params.angle, params.min_size);

        let (step_v, max_v) = self.compute_step(range_v, a_step);
        let (step_u, max_u) = self.compute_step(range_u, a_step);
        let half_du = step_u * 0.5;

        let mut shift = false;
        let mut nodes = Vec::new();
        let mut pas_v = range_v.0 + step_v;
        while pas_v < max_v {
            shift = !shift;
            let d = if shift { half_du } else { 0.0 };
            let mut pas_u = range_u.0 + d;
            while pas_u < max_u {
                nodes.push(GpPnt2d::new(pas_u, pas_v));
                pas_u += step_u;
            }
            pas_v += step_v;
        }
        Some(nodes)
    }
}

/// UV range splitter — tracks the U/V parameters of boundary points. Source:
/// `BRepMesh_UVParamRangeSplitter`.
pub struct UVParamRangeSplitter {
    pub(super) inner: DefaultRangeSplitter,
    pub(super) u_params: ParamSet,
    pub(super) v_params: ParamSet,
}

impl UVParamRangeSplitter {
    /// Creates an empty splitter.
    pub fn new() -> Self {
        Self {
            inner: DefaultRangeSplitter::new(),
            u_params: ParamSet::new(),
            v_params: ParamSet::new(),
        }
    }
}

impl RangeSplitter for UVParamRangeSplitter {
    fn base(&self) -> &DefaultRangeSplitter {
        &self.inner
    }
    fn base_mut(&mut self) -> &mut DefaultRangeSplitter {
        &mut self.inner
    }

    fn reset(&mut self, dface: &MeshFace, params: &MeshParameters) {
        self.reset_base(dface, params);
        self.u_params.clear();
        self.v_params.clear();
    }

    fn parameters_u(&self) -> Option<&ParamSet> {
        Some(&self.u_params)
    }
    fn parameters_v(&self) -> Option<&ParamSet> {
        Some(&self.v_params)
    }
    fn parameters_u_mut(&mut self) -> Option<&mut ParamSet> {
        Some(&mut self.u_params)
    }
    fn parameters_v_mut(&mut self) -> Option<&mut ParamSet> {
        Some(&mut self.v_params)
    }
}

/// Torus surface splitter — U/V periodic, boundary-parameter aware grid.
/// Source: `BRepMesh_TorusRangeSplitter`.
pub struct TorusRangeSplitter {
    pub(super) inner: UVParamRangeSplitter,
}

impl TorusRangeSplitter {
    /// Creates an empty splitter.
    pub fn new() -> Self {
        Self {
            inner: UVParamRangeSplitter::new(),
        }
    }

    /// Fills a parameter sequence from the collected params, spaced at least
    /// `aStdStep` apart. Source: `fillParams`.
    pub(super) fn fill_params(
        &self,
        params: &ParamSet,
        range: (f64, f64),
        steps_nb: i32,
        scale: f64,
    ) -> Vec<f64> {
        let mut arr: Vec<f64> = params.iter().copied().collect();
        let diff = (range.1 - range.0).abs();
        let mut step = calc_average_duv(&mut arr);
        step = step.max(diff / steps_nb as f64 / 2.0);

        let mut std_step = if arr.is_empty() { 0.0 } else { diff / arr.len() as f64 };
        if step > std_step {
            std_step = step;
        }
        std_step *= scale;

        let mut result = Vec::new();
        for &pp in &arr {
            let is_to_insert = result.iter().all(|&v: &f64| (v - pp).abs() > std_step);
            if is_to_insert {
                result.push(pp);
            }
        }
        result
    }
}

/// Average gap between consecutive sorted parameters. Source:
/// `FUN_CalcAverageDUV`.
pub(super) fn calc_average_duv(p: &mut [f64]) -> f64 {
    p.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut n = 0;
    let mut result = 0.0;
    for i in 1..p.len() {
        let d = (p[i] - p[i - 1]).abs();
        if d > 1e-7 {
            result += d;
            n += 1;
        }
    }
    if n > 0 {
        result / n as f64
    } else {
        -1.0
    }
}

impl RangeSplitter for TorusRangeSplitter {
    fn base(&self) -> &DefaultRangeSplitter {
        &self.inner.inner
    }
    fn base_mut(&mut self) -> &mut DefaultRangeSplitter {
        &mut self.inner.inner
    }

    fn reset(&mut self, dface: &MeshFace, params: &MeshParameters) {
        self.inner.reset(dface, params);
    }

    fn add_point(&mut self, point: GpPnt2d) {
        self.add_point_base(point);
        self.inner.u_params.insert(point.x());
        self.inner.v_params.insert(point.y());
    }

    fn parameters_u(&self) -> Option<&ParamSet> {
        Some(&self.inner.u_params)
    }
    fn parameters_v(&self) -> Option<&ParamSet> {
        Some(&self.inner.v_params)
    }

    fn generate_surface_nodes(&self, params: &MeshParameters) -> Option<Vec<GpPnt2d>> {
        let range_u = self.base().range_u();
        let range_v = self.base().range_v();
        let diff_u = range_u.1 - range_u.0;
        let diff_v = range_v.1 - range_v.0;

        let (r_major, r_minor) = self
            .surface()
            .map(|s| torus_radii(s.as_ref()))
            .unwrap_or((0.0, 0.0));
        let deflection = self.deflection();
        let r = r_minor;
        let R = r_major;

        let old_dv = arc_angular_step(r, deflection, params.angle, params.min_size);
        let dv = old_dv;

        let nb_v = (diff_v / dv) as i32;
        let nb_v = nb_v.max(2);
        let dv = diff_v / (nb_v + 1) as f64;

        let ru = R + r;
        let du = if ru > 1e-16 {
            let du0 = arc_angular_step(ru, deflection, params.angle, params.min_size);
            let aa = (du0 * du0 + old_dv * old_dv).sqrt();
            if aa < RESOLUTION {
                return None;
            }
            du0 * old_dv.min(du0) / aa
        } else {
            dv
        };

        let mut nb_u = (diff_u / du) as i32;
        nb_u = nb_u.max(2);
        let ratio_terms = if diff_v * r != 0.0 {
            nb_v as f64 * diff_u * R / (diff_v * r) / 5.0
        } else {
            0.0
        };
        nb_u = nb_u.max(ratio_terms as i32);
        let du = diff_u / (nb_u + 1) as f64;

        let param_u = if R < r {
            (0..=nb_u).map(|i| range_u.0 + i as f64 * du).collect()
        } else {
            self.fill_params(&self.inner.u_params, range_u, nb_u, 0.5)
        };
        let param_v = self.fill_params(&self.inner.v_params, range_v, nb_v, 2.0 / 3.0);

        let new_range_u = (range_u.0 + du * 0.1, range_u.1 - du * 0.1);
        let new_range_v = (range_v.0 + dv * 0.1, range_v.1 - dv * 0.1);

        let mut nodes = Vec::new();
        for &pas_u in &param_u {
            if pas_u >= new_range_u.0 && pas_u < new_range_u.1 {
                for &pas_v in &param_v {
                    if pas_v >= new_range_v.0 && pas_v < new_range_v.1 {
                        nodes.push(GpPnt2d::new(pas_u, pas_v));
                    }
                }
            }
        }
        Some(nodes)
    }
}

/// NURBS / Bezier / BSpline surface splitter — interval-based interior grid.
/// Source: `BRepMesh_NURBSRangeSplitter`.
pub struct NURBSRangeSplitter {
    pub(super) inner: UVParamRangeSplitter,
    pub(super) surface_type: SurfaceType,
}

impl NURBSRangeSplitter {
    /// Creates an empty splitter.
    pub fn new() -> Self {
        Self {
            inner: UVParamRangeSplitter::new(),
            surface_type: SurfaceType::OtherSurface,
        }
    }
}

impl RangeSplitter for NURBSRangeSplitter {
    fn base(&self) -> &DefaultRangeSplitter {
        &self.inner.inner
    }
    fn base_mut(&mut self) -> &mut DefaultRangeSplitter {
        &mut self.inner.inner
    }

    fn adjust_range(&mut self) {
        self.adjust_range_base();
        self.surface_type = self
            .surface()
            .map(|s| classify_surface(s.as_ref()))
            .unwrap_or(SurfaceType::OtherSurface);
        if self.surface_type == SurfaceType::BezierSurface {
            let (ru0, ru1) = self.base().range_u();
            let (rv0, rv1) = self.base().range_v();
            self.base_mut()
                .set_valid(ru0 >= -0.5 && ru1 <= 1.5 && rv0 >= -0.5 && rv1 <= 1.5);
        }
    }

    fn generate_surface_nodes(&self, params: &MeshParameters) -> Option<Vec<GpPnt2d>> {
        generate_nurbs_grid(self, params)
    }

    fn parameters_u(&self) -> Option<&ParamSet> {
        Some(&self.inner.u_params)
    }
    fn parameters_v(&self) -> Option<&ParamSet> {
        Some(&self.inner.v_params)
    }
    fn parameters_u_mut(&mut self) -> Option<&mut ParamSet> {
        Some(&mut self.inner.u_params)
    }
    fn parameters_v_mut(&mut self) -> Option<&mut ParamSet> {
        Some(&mut self.inner.v_params)
    }

    fn get_undefined_interval_nb(&self, is_u: bool, _continuity: u8) -> i32 {
        let Some(surf) = self.surface() else {
            return 1;
        };
        let n = if is_u { surf.nb_u_poles() } else { surf.nb_v_poles() };
        (n - 1).max(1)
    }
}

/// Splitter for surfaces that look like NURBS but expose no poles or other
/// interval characteristics — a single interval per direction. Source:
/// `BRepMesh_UndefinedRangeSplitter`.
pub struct UndefinedRangeSplitter {
    pub(super) inner: NURBSRangeSplitter,
}

impl UndefinedRangeSplitter {
    /// Creates an empty splitter.
    pub fn new() -> Self {
        Self {
            inner: NURBSRangeSplitter::new(),
        }
    }
}
