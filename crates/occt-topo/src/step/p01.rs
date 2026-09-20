use super::prelude::*;
use super::*;

// ---------------------------------------------------------------------------
// Formatting helpers
// ---------------------------------------------------------------------------

/// Render an f64 in STEP real syntax (always carries a decimal point in the
/// mantissa; the shortest round-tripping form).

pub(super) fn step_real(x: f64) -> String {
    if !x.is_finite() {
        return "0.".into();
    }
    let s = format!("{:?}", x);
    // Rust's `{:?}` gives e.g. "0.0", "2.5", "1e-9", "1e20", "inf".
    if let Some(pos) = s.find(['e', 'E']) {
        let mant = &s[..pos];
        if !mant.contains('.') {
            return format!("{mant}.{}", &s[pos..]);
        }
        s
    } else if s.contains('.') {
        s
    } else {
        format!("{s}.")
    }
}

/// Escape a label for a STEP string literal (single quotes doubled).
pub(super) fn esc_str(s: &str) -> String {
    s.replace('\'', "''")
}

/// `#1,#2,#3` from a slice of ids.
pub(super) fn join_refs(ids: &[usize]) -> String {
    ids.iter().map(|i| format!("#{i}")).collect::<Vec<_>>().join(",")
}

pub(super) fn dir_x() -> GpDir {
    GpDir::new(1.0, 0.0, 0.0).unwrap()
}
pub(super) fn dir_y() -> GpDir {
    GpDir::new(0.0, 1.0, 0.0).unwrap()
}
pub(super) fn dir_z() -> GpDir {
    GpDir::new(0.0, 0.0, 1.0).unwrap()
}

/// 3×3 determinant.
pub(super) fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Circumcenter of three non-collinear 3D points (perpendicular-bisector
/// system solved by Cramer's rule).
pub(super) fn circle_center3(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> Option<GpPnt> {
    let d1 = GpVec::from_pnts(a, b);
    let d2 = GpVec::from_pnts(a, c);
    let n = d1.xyz().crossed(d2.xyz());
    if n.square_modulus() < 1e-30 {
        return None;
    }
    let n2 = |p: &GpPnt| p.x() * p.x() + p.y() * p.y() + p.z() * p.z();
    let mat = [
        [d1.x(), d1.y(), d1.z()],
        [d2.x(), d2.y(), d2.z()],
        [n.x, n.y, n.z],
    ];
    let rhs = [
        0.5 * (n2(b) - n2(a)),
        0.5 * (n2(c) - n2(a)),
        a.x() * n.x + a.y() * n.y + a.z() * n.z,
    ];
    let d = det3(&mat);
    if d.abs() < 1e-30 {
        return None;
    }
    let mut o = [0.0; 3];
    for k in 0..3 {
        let mut m = mat;
        for i in 0..3 {
            m[i][k] = rhs[i];
        }
        o[k] = det3(&m) / d;
    }
    Some(GpPnt::new(o[0], o[1], o[2]))
}

// ---------------------------------------------------------------------------
// Geometry classification
// ---------------------------------------------------------------------------

/// Analytic surface family. Plane/sphere reuse `brep_surface`'s tested
/// classifiers; cylinder/cone/torus are detected from constant-v rings.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum SurfKind {
    Plane,
    Sphere,
    Cylinder,
    Cone,
    Torus,
    Other,
}

pub(super) fn surf_bounds(s: &dyn Surface) -> (f64, f64, f64, f64) {
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let clamp = |a: f64, b: f64| {
        if a.is_finite() && b.is_finite() && b > a {
            (a, b)
        } else {
            (-1.0, 1.0)
        }
    };
    let (u0, u1) = clamp(u0, u1);
    let (v0, v1) = clamp(v0, v1);
    (u0, u1, v0, v1)
}

/// Center + radius of the ring obtained by sweeping `u` at fixed `v`.
pub(super) fn ring_center_radius(s: &dyn Surface, u0: f64, u1: f64, v: f64) -> Option<(GpPnt, f64)> {
    let span = u1 - u0;
    let p0 = s.d0(u0, v);
    let p1 = s.d0(u0 + 0.25 * span, v);
    let p2 = s.d0(u0 + 0.5 * span, v);
    let c = circle_center3(&p0, &p1, &p2)?;
    Some((c, c.distance(&p0)))
}

/// Axis direction of a ruled-of-revolution surface from a constant-v ring's
/// plane normal.
pub(super) fn ring_axis(s: &dyn Surface, u0: f64, u1: f64, v: f64) -> Option<GpDir> {
    let span = u1 - u0;
    let p0 = s.d0(u0, v);
    let p1 = s.d0(u0 + 0.25 * span, v);
    let p2 = s.d0(u0 + 0.5 * span, v);
    let n = GpVec::from_pnts(&p0, &p1)
        .xyz()
        .crossed(GpVec::from_pnts(&p0, &p2).xyz());
    GpDir::from_vec(&GpVec::from_xyz(&n)).ok()
}

/// Whether the v-lines (constant u) are straight — cone generators are,
/// torus minor circles are not.
pub(super) fn v_line_straight(s: &dyn Surface, u: f64, v0: f64, v1: f64) -> bool {
    let span = v1 - v0;
    let p0 = s.d0(u, v0);
    let p1 = s.d0(u, v0 + 0.25 * span);
    let p2 = s.d0(u, v0 + 0.5 * span);
    let a = GpVec::from_pnts(&p0, &p1);
    let b = GpVec::from_pnts(&p0, &p2);
    let cross = a.xyz().crossed(b.xyz()).modulus();
    cross < 1e-6 * (a.magnitude() * b.magnitude()).max(1e-12)
}

pub(super) fn surface_kind(s: &dyn Surface) -> SurfKind {
    match classify_surface(s) {
        SurfaceKind::Plane => return SurfKind::Plane,
        SurfaceKind::Sphere => return SurfKind::Sphere,
        _ => {}
    }
    let (u0, u1, v0, v1) = surf_bounds(s);
    let vm = 0.5 * (v0 + v1);
    let v2 = vm + 0.25 * (v1 - v0);
    if let (Some((c0, r0)), Some((c1, r1))) = (
        ring_center_radius(s, u0, u1, vm),
        ring_center_radius(s, u0, u1, v2),
    ) {
        if (r1 - r0).abs() <= 1e-6 * r0.abs().max(1.0) && c0.distance(&c1) > 1e-9 {
            return SurfKind::Cylinder;
        }
        if v_line_straight(s, u0, v0, v1) {
            return SurfKind::Cone;
        }
        return SurfKind::Torus;
    }
    SurfKind::Other
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

/// STEP physical-file header fields (`FILE_DESCRIPTION` / `FILE_NAME` /
/// `FILE_SCHEMA`), used by [`StepWriter::set_header`].
#[derive(Debug, Clone)]
pub struct StepHeader {
    pub description: String,
    pub name: String,
    pub timestamp: String,
    pub author: String,
    pub organization: String,
    pub preprocessor: String,
    pub originator: String,
    pub schema: String,
}

impl Default for StepHeader {
    fn default() -> Self {
        Self {
            description: "BRep".into(),
            name: "model.step".into(),
            timestamp: "2026-07-31T00:00:00".into(),
            author: String::new(),
            organization: String::new(),
            preprocessor: "rust".into(),
            originator: String::new(),
            schema: "AUTOMOTIVE_DESIGN".into(),
        }
    }
}

/// Incremental STEP writer — emits `#N=...;` records with an internal counter.
pub struct StepWriter {
    pub(super) next_id: usize,
    pub(super) lines: Vec<String>,
    pub(super) header: StepHeader,
}

impl StepWriter {
    pub fn new() -> Self {
        Self {
            next_id: 1,
            lines: Vec::new(),
            header: StepHeader::default(),
        }
    }

    /// Replace the physical-file header (FILE_DESCRIPTION / FILE_NAME /
    /// FILE_SCHEMA). By default a minimal header is emitted.
    pub fn set_header(&mut self, header: &StepHeader) {
        self.header = header.clone();
    }

    /// Append a STEP comment (`/* ... */`) to the DATA section. Comments are
    /// ignored by parsers but are useful for provenance and debugging.
    pub fn comment(&mut self, text: &str) {
        let clean = text.replace("*/", "* /");
        self.lines.push(format!("/* {clean} */"));
    }

    /// Emit one entity body (`TYPE(a,b,...)`) and return its record id.
    pub fn emit(&mut self, entity: String) -> usize {
        let id = self.next_id;
        self.next_id += 1;
        self.lines.push(format!("#{id}={entity};"));
        id
    }

    pub fn add_cartesian_point(&mut self, p: &GpPnt) -> usize {
        self.emit(format!(
            "CARTESIAN_POINT('',({},{},{}))",
            step_real(p.x()),
            step_real(p.y()),
            step_real(p.z())
        ))
    }

    pub fn add_direction(&mut self, d: &GpDir) -> usize {
        self.emit(format!(
            "DIRECTION('',({},{},{}))",
            step_real(d.x()),
            step_real(d.y()),
            step_real(d.z())
        ))
    }

    /// VECTOR(name, direction, magnitude).
    pub fn add_vector(&mut self, d: &GpDir, magnitude: f64) -> usize {
        let did = self.add_direction(d);
        self.emit(format!("VECTOR('',#{did},{})", step_real(magnitude)))
    }

    /// AXIS2_PLACEMENT_3D(name, location, axis, ref_direction).
    pub fn add_axis2_placement_3d(&mut self, loc: &GpPnt, axis: &GpDir, ref_dir: &GpDir) -> usize {
        let lid = self.add_cartesian_point(loc);
        let aid = self.add_direction(axis);
        let rid = self.add_direction(ref_dir);
        self.emit(format!("AXIS2_PLACEMENT_3D('',#{lid},#{aid},#{rid})"))
    }

    /// Assemble the complete physical file.
    pub fn finish(self) -> String {
        let h = self.header;
        let desc = esc_str(&h.description);
        let name = esc_str(&h.name);
        let ts = esc_str(&h.timestamp);
        let author = esc_str(&h.author);
        let org = esc_str(&h.organization);
        let pre = esc_str(&h.preprocessor);
        let orig = esc_str(&h.originator);
        let schema = esc_str(&h.schema);
        let mut out = String::new();
        out.push_str("ISO-10303-21;\n");
        out.push_str("HEADER;\n");
        out.push_str(&format!("FILE_DESCRIPTION(('{desc}'),'2;1');\n"));
        out.push_str(&format!(
            "FILE_NAME('{name}','{ts}',('{author}'),('{org}'),'{pre}','{orig}','');\n"
        ));
        out.push_str(&format!("FILE_SCHEMA(('{schema}'));\n"));
        out.push_str("ENDSEC;\n");
        out.push_str("DATA;\n");
        for l in &self.lines {
            out.push_str(l);
            out.push('\n');
        }
        out.push_str("ENDSEC;\n");
        out.push_str("END-ISO-10303-21;\n");
        out
    }
}

/// Serialize a single shape with a custom physical-file header.
///
/// The shape representation and product scaffolding are written as in
/// [`write_step_with_name`], but the `HEADER` section (`FILE_DESCRIPTION`,
/// `FILE_NAME`, `FILE_SCHEMA`) comes from `header` instead of the defaults.
pub fn write_step_with_header(shape: &TopoShape, name: &str, header: &StepHeader) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    ctx.write_shape_named(name, shape)
        .ok_or_else(|| "write_step_with_header: nothing to write".to_string())?;
    ctx.w.set_header(header);
    Ok(ctx.finish())
}

// ---------------------------------------------------------------------------
// B-spline knot helpers
// ---------------------------------------------------------------------------

/// Split an expanded knot vector into its distinct values, in order.
///
/// STEP's `b_spline_curve_with_knots` / `b_spline_surface_with_knots` store
/// knots as a list of distinct values plus a parallel list of multiplicities;
/// this recovers the values half of that pair from the internal expanded form.
pub(super) fn unique_knots(knots: &[f64]) -> Vec<f64> {
    let mut out: Vec<f64> = Vec::new();
    for &k in knots {
        if out.last() != Some(&k) {
            out.push(k);
        }
    }
    out
}

/// The multiplicity of each distinct knot in an expanded knot vector.
pub(super) fn knot_multiplicities(knots: &[f64]) -> Vec<usize> {
    let mut mults: Vec<usize> = Vec::new();
    let mut i = 0;
    while i < knots.len() {
        let mut j = i + 1;
        while j < knots.len() && knots[j] == knots[i] {
            j += 1;
        }
        mults.push(j - i);
        i = j;
    }
    mults
}

/// Reconstruct an expanded knot vector from STEP's mult/knot lists.
pub(super) fn expand_knots(mults: &[usize], knots: &[f64]) -> Vec<f64> {
    let mut out = Vec::new();
    for (i, &m) in mults.iter().enumerate() {
        let k = knots.get(i).copied().unwrap_or(0.0);
        for _ in 0..m {
            out.push(k);
        }
    }
    out
}

/// A clamped uniform knot vector for `n` poles of `degree` — the default used
/// by STEP `B_SPLINE_CURVE` / `B_SPLINE_SURFACE` records that omit explicit
/// knots (multiplicity `degree + 1` at each end, interior knots evenly spaced).
pub(super) fn uniform_knots_for(n: usize, degree: usize) -> Vec<f64> {
    let len = n + degree + 1;
    let interior = if n >= degree + 1 { n - degree - 1 } else { 0 };
    let mut k = Vec::with_capacity(len);
    for _ in 0..=degree {
        k.push(0.0);
    }
    for i in 1..=interior {
        k.push(i as f64 / (interior + 1) as f64);
    }
    while k.len() < len {
        k.push(1.0);
    }
    k
}

// ---------------------------------------------------------------------------
// B-spline entity writers
// ---------------------------------------------------------------------------

/// Emit a `B_SPLINE_CURVE_WITH_KNOTS` entity and return its record id.
///
/// Rational curves carry a weights list in argument 3; polynomial curves use
/// the `SELF` keyword (weights default to 1.0 on read). The curve form,
/// closedness and self-intersection flags are written as unsensed defaults so
/// the record stays minimal while remaining schema-valid. The knot vector is
/// written as its distinct values plus multiplicities.
pub fn write_bspline_curve(step: &mut StepWriter, curve: &GeomBSplineCurve) -> Result<usize, String> {
    let pole_refs: Vec<usize> = curve
        .poles
        .iter()
        .map(|p| step.add_cartesian_point(p))
        .collect();
    let mults = knot_multiplicities(&curve.knots);
    let knots = unique_knots(&curve.knots);
    let weights = match &curve.weights {
        Some(w) => format!(
            "({})",
            w.iter().map(|wi| step_real(*wi)).collect::<Vec<_>>().join(",")
        ),
        None => "SELF".to_string(),
    };
    Ok(step.emit(format!(
        "B_SPLINE_CURVE_WITH_KNOTS('',{},({}),{},UNSPECIFIED,.F.,.F.,({}),({}),UNSPECIFIED)",
        curve.degree,
        join_refs(&pole_refs),
        weights,
        knots.iter().map(|k| step_real(*k)).collect::<Vec<_>>().join(","),
        mults.iter().map(|m| m.to_string()).collect::<Vec<_>>().join(","),
    )))
}

/// Emit a `B_SPLINE_SURFACE_WITH_KNOTS` entity and return its record id.
///
/// The pole grid is written row-major in `u` (each row is a fixed-`u` strip of
/// `v`-varying control points). Weights (when the surface is rational) mirror
/// the grid layout; polynomial surfaces use `SELF`. Each knot direction is
/// split into distinct values plus multiplicities, and the surface form /
/// closedness / self-intersection flags are written as unsensed defaults.
pub fn write_bspline_surface(step: &mut StepWriter, s: &GeomBSplineSurface) -> Result<usize, String> {
    let mut rows = Vec::with_capacity(s.poles.len());
    for row in &s.poles {
        let refs: Vec<usize> = row.iter().map(|p| step.add_cartesian_point(p)).collect();
        rows.push(format!("({})", join_refs(&refs)));
    }
    let poles_grid = format!("({})", rows.join(","));
    let weights = match &s.weights {
        Some(w) => format!(
            "({})",
            w.iter()
                .map(|row| format!(
                    "({})",
                    row.iter().map(|wi| step_real(*wi)).collect::<Vec<_>>().join(",")
                ))
                .collect::<Vec<_>>()
                .join(",")
        ),
        None => "SELF".to_string(),
    };
    let u_knots = unique_knots(&s.knots_u);
    let v_knots = unique_knots(&s.knots_v);
    let u_mult = knot_multiplicities(&s.knots_u);
    let v_mult = knot_multiplicities(&s.knots_v);
    let real_list = |v: &[f64]| v.iter().map(|k| step_real(*k)).collect::<Vec<_>>().join(",");
    let int_list = |v: &[usize]| v.iter().map(|m| m.to_string()).collect::<Vec<_>>().join(",");
    Ok(step.emit(format!(
        "B_SPLINE_SURFACE_WITH_KNOTS('',{},{},{},{},UNSPECIFIED,.F.,.F.,.F.,({}),({}),({}),({}),UNSPECIFIED)",
        s.deg_u,
        s.deg_v,
        poles_grid,
        weights,
        real_list(&u_knots),
        int_list(&u_mult),
        real_list(&v_knots),
        int_list(&v_mult),
    )))
}

/// Emit a `TRIMMED_CURVE` entity restricting a basis curve to `[a, b]`.
///
/// `master_representation` is `PARAMETER`, so the trim bounds are parameter
/// values of the basis curve (rather than points lying on it). The two bounds
/// are normalised to ascending order.
pub fn write_trimmed_curve(step: &mut StepWriter, curve_ref: usize, a: f64, b: f64) -> Result<usize, String> {
    let (a, b) = (a.min(b), a.max(b));
    Ok(step.emit(format!(
        "TRIMMED_CURVE('',#{curve_ref},1,{},{},PARAMETER)",
        step_real(a),
        step_real(b)
    )))
}

/// Emit an `OFFSET_CURVE_3D` entity: `curve_ref` shifted by `offset` along the
/// constant direction `dir_ref`.
///
/// `self_intersect` and `curve_form` are left `UNSPECIFIED` (a caller that
/// knows the basis's analytic form may override them).
pub fn write_offset_curve(
    step: &mut StepWriter,
    curve_ref: usize,
    offset: f64,
    dir_ref: usize,
) -> Result<usize, String> {
    Ok(step.emit(format!(
        "OFFSET_CURVE_3D('',#{curve_ref},#{dir_ref},{},UNSPECIFIED,UNSPECIFIED)",
        step_real(offset)
    )))
}

/// Emit a `POLYLINE` entity through `points` and return its record id.
///
/// A STEP `POLYLINE` is a connected sequence of `CARTESIAN_POINT`s; the reader
/// reconstructs it as a degree-1 B-spline curve, so an edge built from a
/// polyline round-trips as a piecewise-linear curve.
pub fn write_polyline(step: &mut StepWriter, points: &[GpPnt]) -> Result<usize, String> {
    if points.len() < 2 {
        return Err("write_polyline: need at least 2 points".into());
    }
    let refs: Vec<usize> = points
        .iter()
        .map(|p| step.add_cartesian_point(p))
        .collect();
    Ok(step.emit(format!("POLYLINE('',({}))", join_refs(&refs))))
}

// ---------------------------------------------------------------------------
// Conic entity writers (exact `GeomToStep_Make*` transcriptions)
// ---------------------------------------------------------------------------

/// `GeomToStep_MakeCircle` (`GeomToStep_MakeCircle.cxx`): the
/// `AXIS2_PLACEMENT_3D` comes from `C.Position()` (Location / Direction /
/// XDirection) and the radius from `C.Radius()`. No sampling: the previous body
/// recovered the center from three curve samples and the radius from their
/// distance.
pub(super) fn emit_circle_entity(step: &mut StepWriter, c: &GpCirc) -> usize {
    let pos = c.position();
    let ax = step.add_axis2_placement_3d(&pos.location(), &pos.direction(), &pos.x_direction());
    step.emit(format!("CIRCLE('',#{ax},{})", step_real(c.radius())))
}

/// `GeomToStep_MakeEllipse`: placement from `E.Position()`, semi-axes from
/// `E.MajorRadius()` / `E.MinorRadius()`.
pub(super) fn emit_ellipse_entity(step: &mut StepWriter, e: &GpElips) -> usize {
    let pos = e.pos;
    let ax = step.add_axis2_placement_3d(&pos.location(), &pos.direction(), &pos.x_direction());
    step.emit(format!(
        "ELLIPSE('',#{ax},{},{})",
        step_real(e.major_radius),
        step_real(e.minor_radius)
    ))
}

/// `GeomToStep_MakeParabola`: placement from `P.Position()`, focal length from
/// `P.Focal()`.
pub(super) fn emit_parabola_entity(step: &mut StepWriter, p: &GpParab) -> usize {
    let pos = p.pos;
    let ax = step.add_axis2_placement_3d(&pos.location(), &pos.direction(), &pos.x_direction());
    step.emit(format!("PARABOLA('',#{ax},{})", step_real(p.focal)))
}

/// `GeomToStep_MakeHyperbola`: placement from `H.Position()`, semi-axes from
/// `H.MajorRadius()` / `H.MinorRadius()`.
pub(super) fn emit_hyperbola_entity(step: &mut StepWriter, h: &GpHypr) -> usize {
    let pos = h.pos;
    let ax = step.add_axis2_placement_3d(&pos.location(), &pos.direction(), &pos.x_direction());
    step.emit(format!(
        "HYPERBOLA('',#{ax},{},{})",
        step_real(h.major_radius),
        step_real(h.minor_radius)
    ))
}

/// Write the full analytic (conic) parameterisation of `curve` if it is a
/// circle, ellipse, hyperbola or parabola. Returns `Ok(None)` for non-conic
/// curves (lines, B-splines, generic trimmed/offset geometry).
///
/// Type dispatch follows `GeomToStep_MakeCurve.cxx:60-65` (`Geom_Conic` →
/// `GeomToStep_MakeConic`) and `GeomToStep_MakeConic.cxx`: `IsKind` order
/// Circle → Ellipse → Hyperbola → Parabola. `Geom_Line` and the bounded curves
/// are handled by the caller (`GeomToStep_MakeLine` / `MakeBoundedCurve`,
/// `MakeCurve.cxx:54-59` and `:94-99`), which is why a line yields `None` here.
pub fn write_conic_params(step: &mut StepWriter, curve: &dyn Curve) -> Result<Option<usize>, String> {
    Ok(if let Some(c) = curve.gp_circ() {
        Some(emit_circle_entity(step, &c))
    } else if let Some(e) = curve.gp_ellipse() {
        Some(emit_ellipse_entity(step, &e))
    } else if let Some(h) = curve.gp_hyperbola() {
        Some(emit_hyperbola_entity(step, &h))
    } else if let Some(p) = curve.gp_parabola() {
        Some(emit_parabola_entity(step, &p))
    } else {
        None
    })
}
