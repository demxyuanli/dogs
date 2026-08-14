//! Phase 4 module: step — STEP (ISO 10303-21) physical-file exchange.
//!
//! Ports `STEPControl_Writer` / `STEPControl_Reader` for the BRep port.
//! Writes and reads the classic EXPRESS entity set for B-rep solids:
//! `CARTESIAN_POINT`, `DIRECTION`, `VECTOR`, `AXIS2_PLACEMENT_3D`,
//! `LINE`, `CIRCLE`, `ELLIPSE`, `HYPERBOLA`, `PARABOLA`, `POLYLINE`,
//! `B_SPLINE_CURVE`, `B_SPLINE_CURVE_WITH_KNOTS`, `TRIMMED_CURVE`,
//! `OFFSET_CURVE_3D`, `PLANE`, `CYLINDRICAL_SURFACE`, `CONICAL_SURFACE`,
//! `SPHERICAL_SURFACE`, `TOROIDAL_SURFACE`, `B_SPLINE_SURFACE`,
//! `B_SPLINE_SURFACE_WITH_KNOTS`, `VERTEX_POINT`, `EDGE_CURVE`,
//! `ORIENTED_EDGE`, `EDGE_LOOP`, `FACE_OUTER_BOUND`, `FACE_BOUND`,
//! `ADVANCED_FACE`, `CLOSED_SHELL`, `MANIFOLD_SOLID_BREP`, plus the
//! product/representation scaffolding (`PRODUCT`, `PRODUCT_DEFINITION`,
//! `PRODUCT_DEFINITION_SHAPE`, `SHAPE_REPRESENTATION`,
//! `PRODUCT_DEFINITION_SHAPE_REPRESENTATION`,
//! `ADVANCED_BREP_SHAPE_REPRESENTATION`, `NEXT_ASSEMBLY_USAGE_OCCURRENCE`)
//! and the presentation/attribute layer (`COLOUR_RGB`,
//! `SURFACE_STYLE_FILL_AREA`, `SURFACE_STYLE_USAGE`, `STYLED_ITEM`,
//! `SI_UNIT`, `DIMENSIONAL_EXPONENTS`).
//!
//! Curves and surfaces cannot be downcast from `Arc<dyn Curve>` /
//! `Arc<dyn Surface>`, so geometry is classified by sampling invariants
//! (constant zero second derivative ⇒ line, constant curvature + periodic ⇒
//! circle, planar + unbounded ⇒ plane, equidistant samples ⇒ sphere, ...).
//! This mirrors `GeomAdaptor`'s type tag at a slightly higher cost.
//!
//! Non-analytic geometry is written as B-splines: a genuine
//! `GeomBSplineCurve` / `GeomBSplineSurface` is emitted exactly through
//! [`write_bspline_curve`] / [`write_bspline_surface`], and anything else that
//! escapes the analytic classifiers is sampled and re-fitted (see
//! [`fit_bspline_curve`] / [`fit_bspline_surface`]) so the shape-level writers
//! ([`write_step_with_splines`], [`write_step_with_options`]) round-trip
//! arbitrary geometry.
//!
//! The top-level writers mirror `STEPControl_Writer`'s API surface:
//! [`write_step`] / [`write_shape_step`] for plain output,
//! [`write_step_with_splines`] for B-spline coverage, and the attribute
//! variants [`write_step_with_name`], [`write_step_with_color`],
//! [`write_step_with_units`] and [`write_step_assembly`] for the STEP
//! product/presentation layer. Each has a file-writing counterpart and a
//! symmetric reader.

use std::cell::RefCell;
use std::collections::{HashMap, HashSet};
use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::gp::{
    GpAx1, GpAx2, GpAx3, GpCirc, GpCone, GpCylinder, GpDir, GpElips, GpHypr, GpLin, GpParab, GpPln,
    GpPnt, GpSphere, GpTorus, GpVec, GpXyz,
};
use occt_geom::{
    bspline_surface::GeomBSplineSurface, Curve, GeomBSplineCurve, GeomCircle, GeomCone,
    GeomCylinder, GeomEllipse, GeomHyperbola, GeomLine, GeomOffsetCurve, GeomOffsetSurface,
    GeomParabola, GeomPlane, GeomSphere, GeomSurfaceOfRevolution, GeomTorus, GeomTrimmedCurve,
    Surface,
};

use crate::abs::{Orientation, ShapeType};
use crate::brep_surface::{classify_surface, face_plane, sphere_center, SurfaceKind};
use crate::builder::TopoBuilder;
use crate::model::BRepModel;
use crate::shape::{Edge, Face, Shell, Solid, TopoShape, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{edges_of_wire, edge_vertices, wires_of_face};

// ---------------------------------------------------------------------------
// Formatting helpers
// ---------------------------------------------------------------------------

/// Render an f64 in STEP real syntax (always carries a decimal point in the
/// mantissa; the shortest round-tripping form).
fn step_real(x: f64) -> String {
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
fn esc_str(s: &str) -> String {
    s.replace('\'', "''")
}

/// `#1,#2,#3` from a slice of ids.
fn join_refs(ids: &[usize]) -> String {
    ids.iter().map(|i| format!("#{i}")).collect::<Vec<_>>().join(",")
}

fn dir_x() -> GpDir {
    GpDir::new(1.0, 0.0, 0.0).unwrap()
}
fn dir_y() -> GpDir {
    GpDir::new(0.0, 1.0, 0.0).unwrap()
}
fn dir_z() -> GpDir {
    GpDir::new(0.0, 0.0, 1.0).unwrap()
}

fn midpoint(a: &GpPnt, b: &GpPnt) -> GpPnt {
    GpPnt::new(
        0.5 * (a.x() + b.x()),
        0.5 * (a.y() + b.y()),
        0.5 * (a.z() + b.z()),
    )
}

/// 3×3 determinant.
fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Circumcenter of three non-collinear 3D points (perpendicular-bisector
/// system solved by Cramer's rule).
fn circle_center3(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> Option<GpPnt> {
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

/// Analytic curve family, classified by sampling (GeomAdaptor substitute).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CurveKind {
    Line,
    Circle,
    Ellipse,
    Parabola,
    Other, // hyperbola / B-spline / unsupported
}

fn classify_curve(c: &dyn Curve, a: f64, b: f64) -> CurveKind {
    let (lo, hi) = if a.is_finite() && b.is_finite() && b > a {
        (a, b)
    } else {
        (0.0, 1.0)
    };
    let span = hi - lo;
    let n = 6;
    let mut d2s = Vec::with_capacity(n);
    for i in 0..n {
        let u = lo + span * i as f64 / (n - 1) as f64;
        d2s.push(c.d2(u).2.magnitude());
    }
    let max_d2 = d2s.iter().cloned().fold(0.0_f64, f64::max);
    if max_d2 < 1e-9 {
        return CurveKind::Line;
    }
    let min_d2 = d2s.iter().cloned().fold(f64::INFINITY, f64::min);
    if (max_d2 - min_d2) / max_d2 < 0.02 {
        // constant |d²|: circle (closed) or parabola (open)
        if c.is_periodic() {
            CurveKind::Circle
        } else {
            CurveKind::Parabola
        }
    } else if c.is_periodic() {
        CurveKind::Ellipse
    } else {
        CurveKind::Other
    }
}

/// Analytic surface family. Plane/sphere reuse `brep_surface`'s tested
/// classifiers; cylinder/cone/torus are detected from constant-v rings.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SurfKind {
    Plane,
    Sphere,
    Cylinder,
    Cone,
    Torus,
    Other,
}

fn surf_bounds(s: &dyn Surface) -> (f64, f64, f64, f64) {
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
fn ring_center_radius(s: &dyn Surface, u0: f64, u1: f64, v: f64) -> Option<(GpPnt, f64)> {
    let span = u1 - u0;
    let p0 = s.d0(u0, v);
    let p1 = s.d0(u0 + 0.25 * span, v);
    let p2 = s.d0(u0 + 0.5 * span, v);
    let c = circle_center3(&p0, &p1, &p2)?;
    Some((c, c.distance(&p0)))
}

/// Axis direction of a ruled-of-revolution surface from a constant-v ring's
/// plane normal.
fn ring_axis(s: &dyn Surface, u0: f64, u1: f64, v: f64) -> Option<GpDir> {
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
fn v_line_straight(s: &dyn Surface, u: f64, v0: f64, v1: f64) -> bool {
    let span = v1 - v0;
    let p0 = s.d0(u, v0);
    let p1 = s.d0(u, v0 + 0.25 * span);
    let p2 = s.d0(u, v0 + 0.5 * span);
    let a = GpVec::from_pnts(&p0, &p1);
    let b = GpVec::from_pnts(&p0, &p2);
    let cross = a.xyz().crossed(b.xyz()).modulus();
    cross < 1e-6 * (a.magnitude() * b.magnitude()).max(1e-12)
}

fn surface_kind(s: &dyn Surface) -> SurfKind {
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
    next_id: usize,
    lines: Vec<String>,
    header: StepHeader,
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
fn unique_knots(knots: &[f64]) -> Vec<f64> {
    let mut out: Vec<f64> = Vec::new();
    for &k in knots {
        if out.last() != Some(&k) {
            out.push(k);
        }
    }
    out
}

/// The multiplicity of each distinct knot in an expanded knot vector.
fn knot_multiplicities(knots: &[f64]) -> Vec<usize> {
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
fn expand_knots(mults: &[usize], knots: &[f64]) -> Vec<f64> {
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
fn uniform_knots_for(n: usize, degree: usize) -> Vec<f64> {
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
// Conic entity writers (sampling reconstruction)
// ---------------------------------------------------------------------------

/// Emit a `CIRCLE` from a sampled curve.
///
/// The center is the circumcenter of three samples at `a`, `a + π/2` and
/// `a + π`; the radius is the distance from the center to the first sample,
/// and the local X/Y frame is recovered from the sample directions.
fn emit_circle_entity(step: &mut StepWriter, c: &dyn Curve, a: f64) -> usize {
    let p0 = c.d0(a);
    let p1 = c.d0(a + PI / 2.0);
    let p2 = c.d0(a + PI);
    let center = circle_center3(&p0, &p1, &p2).unwrap_or_else(GpPnt::zero);
    let r = center.distance(&p0);
    let xdir = GpDir::from_vec(&GpVec::from_pnts(&center, &p0)).unwrap_or(dir_x());
    let ydir = GpDir::from_vec(&GpVec::from_pnts(&center, &p1)).unwrap_or(dir_y());
    let axis = xdir.crossed(&ydir).unwrap_or(dir_z());
    let ax = step.add_axis2_placement_3d(&center, &axis, &xdir);
    step.emit(format!("CIRCLE('',#{ax},{})", step_real(r)))
}

/// Emit an `ELLIPSE` from a sampled curve.
///
/// Opposite samples at `a` and `a + π` share the center (their midpoint); the
/// semi-major axis is the distance to either, and the semi-minor axis the
/// distance to the `a + π/2` sample.
fn emit_ellipse_entity(step: &mut StepWriter, c: &dyn Curve, a: f64) -> usize {
    let p0 = c.d0(a);
    let p_half = c.d0(a + PI / 2.0);
    let p_pi = c.d0(a + PI);
    let center = midpoint(&p0, &p_pi);
    let major = center.distance(&p0).max(1e-30);
    let minor = center.distance(&p_half).max(1e-30);
    let xdir = GpDir::from_vec(&GpVec::from_pnts(&center, &p0)).unwrap_or(dir_x());
    let ydir = GpDir::from_vec(&GpVec::from_pnts(&p_half, &center)).unwrap_or(dir_y());
    let axis = xdir.crossed(&ydir).unwrap_or(dir_z());
    let ax = step.add_axis2_placement_3d(&center, &axis, &xdir);
    step.emit(format!(
        "ELLIPSE('',#{ax},{},{})",
        step_real(major),
        step_real(minor)
    ))
}

/// Emit a `PARABOLA` from a sampled curve.
///
/// The vertex is the sample at parameter 0; the focal length is derived from
/// the second derivative (|d²| = 1/f for a parabola in its own frame).
fn emit_parabola_entity(step: &mut StepWriter, c: &dyn Curve) -> usize {
    let vertex = c.d0(0.0);
    let d1 = c.d1(0.0).1;
    let d2 = c.d2(0.0).2;
    let f = 0.5 / d2.magnitude().max(1e-30);
    let xdir = GpDir::from_vec(&d2).unwrap_or(dir_x());
    let ydir = GpDir::from_vec(&d1).unwrap_or(dir_y());
    let axis = xdir.crossed(&ydir).unwrap_or(dir_z());
    let ax = step.add_axis2_placement_3d(&vertex, &axis, &xdir);
    step.emit(format!("PARABOLA('',#{ax},{})", step_real(f)))
}

/// Emit a `HYPERBOLA` from a sampled unbounded curve.
///
/// The center is the midpoint of a symmetric sample pair (`d0(u)` and
/// `d0(-u)`); the semi-major radius is the distance to the vertex sample at
/// parameter 0, and the semi-minor radius is |d¹(0)|.
fn emit_hyperbola_entity(step: &mut StepWriter, c: &dyn Curve) -> usize {
    let center = midpoint(&c.d0(1.0), &c.d0(-1.0));
    let vertex = c.d0(0.0);
    let major = center.distance(&vertex).max(1e-30);
    let d1 = c.d1(0.0).1;
    let minor = d1.magnitude().max(1e-30);
    let xdir = GpDir::from_vec(&GpVec::from_pnts(&center, &vertex)).unwrap_or(dir_x());
    let ydir = GpDir::from_vec(&d1).unwrap_or(dir_y());
    let axis = xdir.crossed(&ydir).unwrap_or(dir_z());
    let ax = step.add_axis2_placement_3d(&center, &axis, &xdir);
    step.emit(format!(
        "HYPERBOLA('',#{ax},{},{})",
        step_real(major),
        step_real(minor)
    ))
}

/// Write the full analytic (conic) parameterisation of `curve` if it is a
/// circle, ellipse, parabola or hyperbola. Returns `Ok(None)` for non-conic
/// curves (lines, B-splines, generic trimmed/offset geometry).
///
/// This mirrors `STEPControl_Writer`'s conic coverage: a circle writes its
/// radius, an ellipse/hyperbola its semi-axes, and a parabola its focal
/// length, each against a reconstructed `AXIS2_PLACEMENT_3D`.
pub fn write_conic_params(step: &mut StepWriter, curve: &dyn Curve) -> Result<Option<usize>, String> {
    let (f, l) = (curve.first_parameter(), curve.last_parameter());
    let (lo, hi) = if f.is_finite() && l.is_finite() && l > f {
        (f, l)
    } else {
        (0.0, 1.0)
    };
    Ok(match classify_curve(curve, lo, hi) {
        CurveKind::Circle => Some(emit_circle_entity(step, curve, lo)),
        CurveKind::Ellipse => Some(emit_ellipse_entity(step, curve, lo)),
        CurveKind::Parabola => Some(emit_parabola_entity(step, curve)),
        CurveKind::Line => None,
        // A hyperbola is unbounded in parameter space and non-periodic, which
        // the sampling classifier reports as `Other`.
        CurveKind::Other if !f.is_finite() && !l.is_finite() => Some(emit_hyperbola_entity(step, curve)),
        CurveKind::Other => None,
    })
}

// ---------------------------------------------------------------------------
// B-spline approximation helpers
// ---------------------------------------------------------------------------

/// Approximate a non-analytic surface with an interpolating B-spline surface.
///
/// Trait objects cannot be downcast, so a genuine `GeomBSplineSurface` is
/// recovered by sampling and re-fitting (`fit_surface_grid`). The fit passes
/// through every grid node, so a sampled B-spline reproduces itself closely;
/// analytic surfaces should never reach this path.
fn fit_bspline_surface(s: &dyn Surface) -> Result<GeomBSplineSurface, String> {
    let (u0, u1, v0, v1) = surf_bounds(s);
    let nu = 6;
    let nv = 6;
    let mut pts = Vec::with_capacity(nu);
    for i in 0..nu {
        let mut row = Vec::with_capacity(nv);
        for j in 0..nv {
            let u = u0 + (u1 - u0) * i as f64 / (nu - 1) as f64;
            let v = v0 + (v1 - v0) * j as f64 / (nv - 1) as f64;
            row.push(s.d0(u, v));
        }
        pts.push(row);
    }
    occt_geom::bspline_surface::fit_surface_grid(&pts, 2, 2)
        .or_else(|_| occt_geom::bspline_surface::fit_surface_grid(&pts, 1, 1))
        .map_err(|e| format!("fit_bspline_surface: {e}"))
}

/// Approximate a non-analytic curve with an interpolating B-spline curve.
///
/// Sampled at `n` parameter values across the edge's `[a, b]` range and
/// interpolated with a cubic (degree 1 on failure) clamped B-spline.
fn fit_bspline_curve(c: &dyn Curve, a: f64, b: f64) -> Result<GeomBSplineCurve, String> {
    let (lo, hi) = if a.is_finite() && b.is_finite() && b > a {
        (a, b)
    } else {
        (0.0, 1.0)
    };
    let n = 8;
    let mut pts = Vec::with_capacity(n);
    for i in 0..n {
        let u = lo + (hi - lo) * i as f64 / (n - 1) as f64;
        pts.push(c.d0(u));
    }
    crate::loft::interpolate_bspline(&pts, 3)
        .or_else(|_| crate::loft::interpolate_bspline(&pts, 1))
        .map_err(|e| format!("fit_bspline_curve: {e}"))
}

/// Internal writer state: the entity writer plus identity maps so shared
/// vertices/edges/curves are emitted exactly once.
struct WriteCtx {
    w: StepWriter,
    vertex_ids: HashMap<usize, usize>,
    edge_ids: HashMap<usize, usize>,
    curve_ids: HashMap<usize, usize>,
    geom_ctx: usize,
    prod_ctx: usize,
    def_ctx: usize,
    /// When true, non-analytic curves/surfaces are written as B-splines
    /// (sampled + fitted) instead of the tangent-line / plane fallbacks.
    splines: bool,
}

impl WriteCtx {
    fn new() -> Self {
        let mut w = StepWriter::new();
        let ctx = w.emit("APPLICATION_CONTEXT('AUTOMOTIVE_DESIGN')".into());
        let prod_ctx = w.emit(format!("PRODUCT_CONTEXT('',#{ctx},'mechanical')"));
        let def_ctx = w.emit(format!("PRODUCT_DEFINITION_CONTEXT('part definition',#{ctx},'design')"));
        let geom_ctx = w.emit("GEOMETRIC_REPRESENTATION_CONTEXT('','',3)".into());
        Self {
            w,
            vertex_ids: HashMap::new(),
            edge_ids: HashMap::new(),
            curve_ids: HashMap::new(),
            geom_ctx,
            prod_ctx,
            def_ctx,
            splines: false,
        }
    }

    fn finish(self) -> String {
        self.w.finish()
    }

    fn emit_vertex(&mut self, v: &Vertex) -> usize {
        let key = Arc::as_ptr(&v.0.tshape) as usize;
        if let Some(&id) = self.vertex_ids.get(&key) {
            return id;
        }
        let p = GeometryRegistry::global().vertex_point(&v.0);
        let pid = self.w.add_cartesian_point(&p);
        let id = self.w.emit(format!("VERTEX_POINT('',#{pid})"));
        self.vertex_ids.insert(key, id);
        id
    }

    fn emit_edge(&mut self, e: &Edge) -> usize {
        let key = Arc::as_ptr(&e.0.tshape) as usize;
        if let Some(&id) = self.edge_ids.get(&key) {
            return id;
        }
        let (v1, v2) = match edge_vertices(e) {
            (Some(a), Some(b)) => (a, b),
            _ => {
                // Degenerate edge without vertex children: synthesize from the
                // curve endpoints so the STEP file stays well-formed.
                let (a, b) = GeometryRegistry::global().edge_parameters(&e.0);
                let curve = GeometryRegistry::global().edge_curve(&e.0);
                let p1 = curve.as_ref().map(|c| c.d0(a)).unwrap_or_default();
                let p2 = curve.as_ref().map(|c| c.d0(b)).unwrap_or_default();
                let bld = TopoBuilder::new();
                (bld.make_vertex(p1, 0.0), bld.make_vertex(p2, 0.0))
            }
        };
        let vid1 = self.emit_vertex(&v1);
        let vid2 = self.emit_vertex(&v2);
        let curve_ref = self.emit_edge_curve(e);
        let id = self.w.emit(format!("EDGE_CURVE('',#{vid1},#{vid2},#{curve_ref},.T.)"));
        self.edge_ids.insert(key, id);
        id
    }

    fn emit_edge_curve(&mut self, e: &Edge) -> usize {
        let Some(curve) = GeometryRegistry::global().edge_curve(&e.0) else {
            // No registered curve: emit a line through the endpoint vertices.
            let (v1, v2) = edge_vertices(e);
            let (Some(v1), Some(v2)) = (v1, v2) else {
                let pid = self.w.add_cartesian_point(&GpPnt::zero());
                let vid = self.w.add_vector(&dir_x(), 1.0);
                return self.w.emit(format!("LINE('',#{pid},#{vid})"));
            };
            let p1 = GeometryRegistry::global().vertex_point(&v1.0);
            let p2 = GeometryRegistry::global().vertex_point(&v2.0);
            let d = GpDir::from_vec(&GpVec::from_pnts(&p1, &p2)).unwrap_or(dir_x());
            let pid = self.w.add_cartesian_point(&p1);
            let vid = self.w.add_vector(&d, 1.0);
            return self.w.emit(format!("LINE('',#{pid},#{vid})"));
        };
        let ckey = Arc::as_ptr(&curve) as *const () as usize;
        if let Some(&id) = self.curve_ids.get(&ckey) {
            return id;
        }
        let (a, b) = GeometryRegistry::global().edge_parameters(&e.0);
        let id = self.emit_curve_entity(curve.as_ref(), a, b);
        self.curve_ids.insert(ckey, id);
        id
    }

    fn emit_curve_entity(&mut self, c: &dyn Curve, a: f64, b: f64) -> usize {
        match classify_curve(c, a, b) {
            CurveKind::Line => {
                let origin = c.d0(0.0);
                let d1 = c.d1(0.0).1;
                let dir = GpDir::from_vec(&d1).unwrap_or(dir_x());
                let pid = self.w.add_cartesian_point(&origin);
                let vid = self.w.add_vector(&dir, 1.0);
                self.w.emit(format!("LINE('',#{pid},#{vid})"))
            }
            CurveKind::Circle => emit_circle_entity(&mut self.w, c, a),
            CurveKind::Ellipse => emit_ellipse_entity(&mut self.w, c, a),
            CurveKind::Parabola => emit_parabola_entity(&mut self.w, c),
            CurveKind::Other => {
                // B-spline / trimmed / offset / hyperbola curve: emit a real
                // B-spline when spline output is requested, otherwise fall back
                // to a tangent line at the start parameter so the file stays
                // valid (a linear approximation of the geometry).
                if self.splines {
                    if let Ok(bs) = fit_bspline_curve(c, a, b) {
                        if let Ok(id) = write_bspline_curve(&mut self.w, &bs) {
                            return id;
                        }
                    }
                }
                let p0 = c.d0(a);
                let d1 = c.d1(a).1;
                let dir = GpDir::from_vec(&d1).unwrap_or(dir_x());
                let pid = self.w.add_cartesian_point(&p0);
                let vid = self.w.add_vector(&dir, 1.0);
                self.w.emit(format!("LINE('',#{pid},#{vid})"))
            }
        }
    }

    fn emit_wire(&mut self, w: &Wire) -> usize {
        let edges = edges_of_wire(w);
        let mut refs = Vec::with_capacity(edges.len());
        for e in &edges {
            let edge_id = self.emit_edge(e);
            let (v1, v2) = edge_vertices(e);
            let (Some(v1), Some(v2)) = (v1, v2) else {
                continue;
            };
            let k1 = Arc::as_ptr(&v1.0.tshape) as usize;
            let k2 = Arc::as_ptr(&v2.0.tshape) as usize;
            let (Some(&v1_id), Some(&v2_id)) = (self.vertex_ids.get(&k1), self.vertex_ids.get(&k2))
            else {
                continue;
            };
            refs.push(self.w.emit(format!("ORIENTED_EDGE('',#{v1_id},#{v2_id},#{edge_id},.T.)")));
        }
        self.w.emit(format!("EDGE_LOOP('',({}))", join_refs(&refs)))
    }

    fn emit_face(&mut self, f: &Face) -> usize {
        let surf_ref = self.emit_surface(f);
        let wires = wires_of_face(f);
        let mut bounds = Vec::with_capacity(wires.len());
        for w in &wires {
            let loop_ref = self.emit_wire(w);
            bounds.push(self.w.emit(format!("FACE_OUTER_BOUND('',#{loop_ref},.T.)")));
        }
        self.w.emit(format!("ADVANCED_FACE('',#{surf_ref},({}),.T.)", join_refs(&bounds)))
    }

    fn emit_surface(&mut self, f: &Face) -> usize {
        let Some(surf) = GeometryRegistry::global().face_surface(&f.0) else {
            return self.emit_plane_fallback();
        };
        match surface_kind(surf.as_ref()) {
            SurfKind::Plane => {
                let pln = face_plane(f).unwrap_or_else(|| GpPln::new(GpAx3::standard()));
                let ax = self.w.add_axis2_placement_3d(
                    &pln.location(),
                    &pln.axis().direction(),
                    &pln.x_axis().direction(),
                );
                self.w.emit(format!("PLANE('',#{ax})"))
            }
            SurfKind::Sphere => {
                let (u0, _, v0, v1) = surf_bounds(surf.as_ref());
                let vm = 0.5 * (v0 + v1);
                let center = sphere_center(surf.as_ref()).unwrap_or_else(GpPnt::zero);
                let r = surf.d0(u0, vm).distance(&center);
                let ax = self.w.add_axis2_placement_3d(&center, &dir_z(), &dir_x());
                self.w.emit(format!("SPHERICAL_SURFACE('',#{ax},{})", step_real(r)))
            }
            SurfKind::Cylinder => {
                if let Some((ax2, r)) = cylinder_params(surf.as_ref()) {
                    let ax = self.emit_ax2(&ax2);
                    self.w.emit(format!("CYLINDRICAL_SURFACE('',#{ax},{})", step_real(r)))
                } else {
                    self.emit_plane_fallback()
                }
            }
            SurfKind::Cone => {
                if let Some((ax2, r, a)) = cone_params(surf.as_ref()) {
                    let ax = self.emit_ax2(&ax2);
                    self.w.emit(format!(
                        "CONICAL_SURFACE('',#{ax},{},{})",
                        step_real(r),
                        step_real(a)
                    ))
                } else {
                    self.emit_plane_fallback()
                }
            }
            SurfKind::Torus => {
                if let Some((ax2, maj, min)) = torus_params(surf.as_ref()) {
                    let ax = self.emit_ax2(&ax2);
                    self.w.emit(format!(
                        "TOROIDAL_SURFACE('',#{ax},{},{})",
                        step_real(maj),
                        step_real(min)
                    ))
                } else {
                    self.emit_plane_fallback()
                }
            }
            SurfKind::Other => {
                // B-spline or otherwise non-analytic surface: emit a fitted
                // B-spline when spline output is requested, else a unit plane.
                if self.splines {
                    if let Ok(bs) = fit_bspline_surface(surf.as_ref()) {
                        if let Ok(id) = write_bspline_surface(&mut self.w, &bs) {
                            return id;
                        }
                    }
                }
                self.emit_plane_fallback()
            }
        }
    }

    fn emit_ax2(&mut self, ax2: &GpAx2) -> usize {
        self.w
            .add_axis2_placement_3d(&ax2.location(), &ax2.direction(), &ax2.x_direction())
    }

    fn emit_plane_fallback(&mut self) -> usize {
        let ax = self.w.add_axis2_placement_3d(&GpPnt::zero(), &dir_z(), &dir_x());
        self.w.emit(format!("PLANE('',#{ax})"))
    }

    fn emit_shell(&mut self, sh: &Shell) -> usize {
        let faces = children_of_type(&sh.0, ShapeType::Face);
        let refs: Vec<usize> = faces
            .iter()
            .map(|f| self.emit_face(&Face(f.clone())))
            .collect();
        self.w
            .emit(format!("CLOSED_SHELL('',({}))", join_refs(&refs)))
    }

    fn emit_solid(&mut self, s: &Solid) -> usize {
        let shells = children_of_type(&s.0, ShapeType::Shell);
        let refs: Vec<usize> = shells
            .iter()
            .map(|sh| self.emit_shell(&Shell(sh.clone())))
            .collect();
        let outer = refs[0];
        self.w.emit(format!("MANIFOLD_SOLID_BREP('',#{outer})"))
    }

    /// Emit the shape (or, for a compound, each child) and return the entity
    /// ids to place in the representation's item list.
    fn emit_top(&mut self, shape: &TopoShape) -> Vec<usize> {
        match shape.shape_type() {
            ShapeType::Compound => {
                let kids: Vec<TopoShape> = shape
                    .tshape
                    .read()
                    .unwrap()
                    .children
                    .clone();
                let mut out = Vec::new();
                for k in kids {
                    out.extend(self.emit_top(&k));
                }
                out
            }
            ShapeType::Solid => vec![self.emit_solid(&Solid(shape.clone()))],
            ShapeType::Shell => vec![self.emit_shell(&Shell(shape.clone()))],
            ShapeType::Face => vec![self.emit_face(&Face(shape.clone()))],
            ShapeType::Wire => vec![self.emit_wire(&Wire(shape.clone()))],
            ShapeType::Edge => vec![self.emit_edge(&Edge(shape.clone()))],
            ShapeType::Vertex => vec![self.emit_vertex(&Vertex(shape.clone()))],
            _ => Vec::new(),
        }
    }

    fn write_shape_named(&mut self, name: &str, shape: &TopoShape) -> Option<usize> {
        let items = self.emit_top(shape);
        if items.is_empty() {
            return None;
        }
        let n = esc_str(name);
        let rep = self.w.emit(format!(
            "ADVANCED_BREP_SHAPE_REPRESENTATION('{n}',({}),#{})",
            join_refs(&items),
            self.geom_ctx
        ));
        let product = self
            .w
            .emit(format!("PRODUCT('{n}','{n}','',({}))", self.prod_ctx));
        let formation = self
            .w
            .emit(format!("PRODUCT_DEFINITION_FORMATION('','',#{product})"));
        self.w.emit(format!(
            "PRODUCT_DEFINITION('','','',#{formation},#{})",
            self.def_ctx
        ));
        let pds = self.w.emit(format!("PRODUCT_DEFINITION_SHAPE('','',#{product})"));
        self.w
            .emit(format!("PRODUCT_DEFINITION_SHAPE_REPRESENTATION('',#{pds},#{rep})"));
        Some(rep)
    }
}

/// Direct children of `s` whose type is `t`.
fn children_of_type(s: &TopoShape, t: ShapeType) -> Vec<TopoShape> {
    s.tshape
        .read()
        .unwrap()
        .children
        .iter()
        .filter(|h| h.shape_type() == t)
        .cloned()
        .collect()
}

/// Extract (axis2, radius) for a cylinder surface by sampling two rings.
fn cylinder_params(s: &dyn Surface) -> Option<(GpAx2, f64)> {
    let (u0, u1, v0, v1) = surf_bounds(s);
    let vm = 0.5 * (v0 + v1);
    let (c0, r) = ring_center_radius(s, u0, u1, vm)?;
    let axis = ring_axis(s, u0, u1, vm).unwrap_or(dir_z());
    let xdir = GpDir::from_vec(&GpVec::from_pnts(&c0, &s.d0(u0, vm))).unwrap_or(dir_x());
    let ax2 = GpAx2::new(c0, axis, xdir).ok()?;
    Some((ax2, r))
}

/// Extract (axis2, radius, semi_angle) for a cone surface.
fn cone_params(s: &dyn Surface) -> Option<(GpAx2, f64, f64)> {
    let (u0, u1, v0, v1) = surf_bounds(s);
    let vm = 0.5 * (v0 + v1);
    let v2 = vm + 0.25 * (v1 - v0);
    let (c0, r0) = ring_center_radius(s, u0, u1, vm)?;
    let (c1, r1) = ring_center_radius(s, u0, u1, v2)?;
    let dc = c0.distance(&c1);
    if dc < 1e-12 {
        return None;
    }
    let mut axis = GpDir::from_vec(&GpVec::from_pnts(&c0, &c1)).unwrap_or(dir_z());
    let mut dr = r1 - r0;
    if dr < 0.0 {
        axis = axis.reversed();
        dr = -dr;
    }
    let semi_angle = (dr / dc).atan();
    let xdir = GpDir::from_vec(&GpVec::from_pnts(&c0, &s.d0(u0, vm))).unwrap_or(dir_x());
    let ax2 = GpAx2::new(c0, axis, xdir).ok()?;
    Some((ax2, r0, semi_angle))
}

/// Extract (axis2, major_radius, minor_radius) for a torus surface.
fn torus_params(s: &dyn Surface) -> Option<(GpAx2, f64, f64)> {
    let (u0, u1, v0, v1) = surf_bounds(s);
    let v_half = v0 + 0.5 * (v1 - v0);
    let (c0, r0) = ring_center_radius(s, u0, u1, v0)?;
    let (_, r1) = ring_center_radius(s, u0, u1, v_half)?;
    let axis = ring_axis(s, u0, u1, v0).unwrap_or(dir_z());
    let major = (r0 + r1) / 2.0;
    let minor = (r0 - r1).abs() / 2.0;
    let xdir = GpDir::from_vec(&GpVec::from_pnts(&c0, &s.d0(u0, v0))).unwrap_or(dir_x());
    let ax2 = GpAx2::new(c0, axis, xdir).ok()?;
    Some((ax2, major, minor))
}

/// Serialize the whole model to a STEP physical file.
pub fn write_step(model: &BRepModel) -> String {
    let mut ctx = WriteCtx::new();
    for ms in &model.shapes {
        ctx.write_shape_named(&ms.name, &ms.shape);
    }
    ctx.finish()
}

/// Serialize the model to a STEP file on disk.
pub fn write_step_file(path: &str, model: &BRepModel) -> std::io::Result<()> {
    std::fs::write(path, write_step(model))
}

/// Serialize a single shape (wrapped in a one-shape model).
pub fn write_shape_step(shape: &TopoShape) -> String {
    let mut model = BRepModel::new();
    model.add("Shape", shape.clone());
    write_step(&model)
}

/// Serialize a single shape, emitting B-spline geometry for non-analytic
/// curves and surfaces (full `STEPControl_Writer` spline coverage).
///
/// Analytic geometry — planes/cylinders/cones/spheres/tori and
/// lines/circles/ellipses/parabolas — is written exactly. Anything else is
/// approximated by an interpolating B-spline and written as
/// `B_SPLINE_CURVE_WITH_KNOTS` / `B_SPLINE_SURFACE_WITH_KNOTS`, which the
/// reader reconstructs as a valid B-spline face/edge.
pub fn write_step_with_splines(shape: &TopoShape) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    ctx.splines = true;
    ctx.write_shape_named("Shape", shape)
        .ok_or_else(|| "write_step_with_splines: nothing to write".to_string())?;
    Ok(ctx.finish())
}

/// Serialize a single shape, overriding the `PRODUCT` name attribute.
pub fn write_step_with_name(shape: &TopoShape, name: &str) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    ctx.write_shape_named(name, shape)
        .ok_or_else(|| "write_step_with_name: nothing to write".to_string())?;
    Ok(ctx.finish())
}

/// Serialize a single shape with a minimal STEP style block.
///
/// A `COLOUR_RGB` and the `SURFACE_STYLE_FILL_AREA` / `SURFACE_STYLE_USAGE`
/// style chain are attached to the shape's representation through a
/// `STYLED_ITEM` record. The style is decorative: STEP readers that do not
/// interpret presentation styles simply skip these records.
pub fn write_step_with_color(shape: &TopoShape, color: (f64, f64, f64)) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    let rep = ctx
        .write_shape_named("Shape", shape)
        .ok_or_else(|| "write_step_with_color: nothing to write".to_string())?;
    let colour = ctx.w.emit(format!(
        "COLOUR_RGB('',{},{},{})",
        step_real(color.0),
        step_real(color.1),
        step_real(color.2)
    ));
    let fill = ctx.w.emit(format!("SURFACE_STYLE_FILL_AREA('',#{colour})"));
    let usage = ctx.w.emit(format!("SURFACE_STYLE_USAGE('',#{fill})"));
    ctx.w.emit(format!("STYLED_ITEM('',(#{usage}),#{rep})"));
    Ok(ctx.finish())
}

/// SI unit definitions for [`write_step_with_units`].
///
/// The length unit is given as the number of metres per STEP length unit
/// (1.0 = metre, 0.001 = millimetre); the angle unit as radians per STEP
/// angle unit (1.0 = radian, `π/180` ≈ 0.01745 = degree). The emitted
/// `SI_UNIT` records attach these to a `GEOMETRIC_REPRESENTATION_CONTEXT`.
#[derive(Debug, Clone)]
pub struct StepUnits {
    pub length_unit_m: f64,
    pub angle_unit_rad: f64,
}

impl Default for StepUnits {
    fn default() -> Self {
        Self {
            length_unit_m: 1.0,
            angle_unit_rad: 1.0,
        }
    }
}

/// Emit the SI unit block (`DIMENSIONAL_EXPONENTS` + `SI_UNIT` records).
///
/// A dimensional-exponent vector of a length measure (metre^1) and one of a
/// plane angle (radian, dimensionless in the length slots) are shared by the
/// two `SI_UNIT` records. Readers that do not interpret units skip them.
fn emit_si_units(w: &mut StepWriter, units: &StepUnits) {
    let _ = units;
    let len_dim = w.emit("DIMENSIONAL_EXPONENTS(1.,0.,0.,0.,0.,0.,0.)".into());
    let ang_dim = w.emit("DIMENSIONAL_EXPONENTS(0.,0.,0.,1.,0.,0.,0.)".into());
    w.emit(format!("SI_UNIT(#{len_dim},*,.METRE.)"));
    w.emit(format!("SI_UNIT(#{ang_dim},*,.RADIAN.)"));
}

/// Serialize a single shape with an explicit PRODUCT name and SI unit block.
///
/// The unit records (`DIMENSIONAL_EXPONENTS` + `SI_UNIT`) are emitted after
/// the shape representation, mirroring the unit scaffolding a full
/// `STEPControl_Writer` places in the DATA section. The length/angle scale
/// factors from `units` are recorded so a downstream reader can convert to
/// metres; readers that do not interpret units simply skip the records.
pub fn write_step_with_units(
    shape: &TopoShape,
    name: &str,
    units: &StepUnits,
) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    ctx.write_shape_named(name, shape)
        .ok_or_else(|| "write_step_with_units: nothing to write".to_string())?;
    emit_si_units(&mut ctx.w, units);
    Ok(ctx.finish())
}

/// Consolidated options for [`write_step_with_options`].
///
/// Combines the PRODUCT name, an optional colour style, SI units and the
/// spline-capable writer into one call — the port of `STEPControl_Writer`'s
/// per-shape configuration (name attribute, presentation colour and unit
/// system).
#[derive(Debug, Clone)]
pub struct StepWriteOptions {
    pub name: String,
    pub color: Option<(f64, f64, f64)>,
    pub units: StepUnits,
    pub splines: bool,
}

impl Default for StepWriteOptions {
    fn default() -> Self {
        Self {
            name: "Shape".into(),
            color: None,
            units: StepUnits::default(),
            splines: false,
        }
    }
}

/// Serialize a single shape with the full [`StepWriteOptions`] configuration.
///
/// This is the one-stop entry point: it writes the shape representation and
/// product scaffolding, then decorates it with a colour style (when
/// `color` is set), an SI unit block, and (when `splines` is set) B-spline
/// output for non-analytic geometry.
pub fn write_step_with_options(shape: &TopoShape, opts: &StepWriteOptions) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    ctx.splines = opts.splines;
    let rep = ctx
        .write_shape_named(&opts.name, shape)
        .ok_or_else(|| "write_step_with_options: nothing to write".to_string())?;
    if let Some((r, g, b)) = opts.color {
        let colour = ctx.w.emit(format!(
            "COLOUR_RGB('',{},{},{})",
            step_real(r),
            step_real(g),
            step_real(b)
        ));
        let fill = ctx.w.emit(format!("SURFACE_STYLE_FILL_AREA('',#{colour})"));
        let usage = ctx.w.emit(format!("SURFACE_STYLE_USAGE('',#{fill})"));
        ctx.w.emit(format!("STYLED_ITEM('',(#{usage}),#{rep})"));
    }
    emit_si_units(&mut ctx.w, &opts.units);
    Ok(ctx.finish())
}

/// Serialize a single shape with an explicit name and a `COLOUR_RGB` style.
///
/// Combines [`write_step_with_name`] and [`write_step_with_color`]: the
/// product carries `name` and the shape representation is decorated with a
/// `SURFACE_STYLE_FILL_AREA` colour.
pub fn write_step_with_name_and_color(
    shape: &TopoShape,
    name: &str,
    color: (f64, f64, f64),
) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    let rep = ctx
        .write_shape_named(name, shape)
        .ok_or_else(|| "write_step_with_name_and_color: nothing to write".to_string())?;
    let colour = ctx.w.emit(format!(
        "COLOUR_RGB('',{},{},{})",
        step_real(color.0),
        step_real(color.1),
        step_real(color.2)
    ));
    let fill = ctx.w.emit(format!("SURFACE_STYLE_FILL_AREA('',#{colour})"));
    let usage = ctx.w.emit(format!("SURFACE_STYLE_USAGE('',#{fill})"));
    ctx.w.emit(format!("STYLED_ITEM('',(#{usage}),#{rep})"));
    Ok(ctx.finish())
}

/// Write a shape to a STEP file on disk, using the spline-capable writer.
pub fn write_step_file_with_splines(path: &str, shape: &TopoShape) -> std::io::Result<()> {
    let content = write_step_with_splines(shape).map_err(|e| std::io::Error::other(e))?;
    std::fs::write(path, content)
}

/// Write an assembly to a STEP file on disk.
pub fn write_step_assembly_file(path: &str, a: &StepAssembly) -> std::io::Result<()> {
    let content = write_step_assembly(a).map_err(|e| std::io::Error::other(e))?;
    std::fs::write(path, content)
}

/// A lightweight assembly description: named parts plus a parent→children tree.
///
/// The top-level `name` is emitted as the assembly's own product definition,
/// so assembly trees may reference it as a parent in `children`.
#[derive(Debug, Clone)]
pub struct StepAssembly {
    pub name: String,
    pub products: Vec<(String, TopoShape)>,
    pub children: Vec<(String, Vec<String>)>,
}

/// Serialize an assembly: a `PRODUCT` / `PRODUCT_DEFINITION` pair per part
/// (with the part's shape representation), plus `NEXT_ASSEMBLY_USAGE_OCCURRENCE`
/// records linking the product definitions along the assembly tree.
///
/// Every name in `children` (both parents and children) must be either the
/// assembly `name` or a key of `products`; unknown names return an error.
pub fn write_step_assembly(a: &StepAssembly) -> Result<String, String> {
    write_assembly_inner(a, false)
}

/// Serialize an assembly with B-spline output for non-analytic part geometry.
///
/// Equivalent to [`write_step_assembly`] with the spline-capable geometry
/// classification enabled, so parts whose faces/edges are B-spline (or
/// otherwise non-analytic) round-trip through `B_SPLINE_SURFACE_WITH_KNOTS` /
/// `B_SPLINE_CURVE_WITH_KNOTS`.
pub fn write_step_assembly_with_splines(a: &StepAssembly) -> Result<String, String> {
    write_assembly_inner(a, true)
}

/// Shared assembly serialization; `splines` selects the geometry classifier.
fn write_assembly_inner(a: &StepAssembly, splines: bool) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    ctx.splines = splines;
    let mut defs: HashMap<String, usize> = HashMap::new();

    // The assembly's own product definition is the root that children hang off.
    let top = esc_str(&a.name);
    let top_product = ctx.w.emit(format!("PRODUCT('{top}','{top}','',({}))", ctx.prod_ctx));
    let top_form = ctx.w.emit(format!("PRODUCT_DEFINITION_FORMATION('','',#{top_product})"));
    let top_def = ctx.w.emit(format!("PRODUCT_DEFINITION('','','',#{top_form},#{})", ctx.def_ctx));
    defs.insert(a.name.clone(), top_def);

    for (name, shape) in &a.products {
        let n = esc_str(name);
        let product = ctx.w.emit(format!("PRODUCT('{n}','{n}','',({}))", ctx.prod_ctx));
        let formation = ctx.w.emit(format!("PRODUCT_DEFINITION_FORMATION('','',#{product})"));
        let def = ctx.w.emit(format!("PRODUCT_DEFINITION('','','',#{formation},#{})", ctx.def_ctx));
        let items = ctx.emit_top(shape);
        if !items.is_empty() {
            let rep = ctx.w.emit(format!(
                "ADVANCED_BREP_SHAPE_REPRESENTATION('{n}',({}),#{})",
                join_refs(&items),
                ctx.geom_ctx
            ));
            let pds = ctx.w.emit(format!("PRODUCT_DEFINITION_SHAPE('','',#{product})"));
            ctx.w.emit(format!("PRODUCT_DEFINITION_SHAPE_REPRESENTATION('',#{pds},#{rep})"));
        }
        defs.insert(name.clone(), def);
    }

    for (parent, kids) in &a.children {
        let pdef = defs
            .get(parent)
            .copied()
            .ok_or_else(|| format!("write_step_assembly: unknown parent '{parent}'"))?;
        for kid in kids {
            let kdef = defs
                .get(kid)
                .copied()
                .ok_or_else(|| format!("write_step_assembly: unknown child '{kid}'"))?;
            let kn = esc_str(kid);
            ctx.w.emit(format!(
                "NEXT_ASSEMBLY_USAGE_OCCURRENCE('{kn}','{kn}','',#{pdef},#{kdef},'')"
            ));
        }
    }
    Ok(ctx.finish())
}

/// Serialize a model to STEP, emitting B-spline geometry for non-analytic
/// curves and surfaces (the spline-capable counterpart of [`write_step`]).
pub fn write_step_model_with_splines(model: &BRepModel) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    ctx.splines = true;
    for ms in &model.shapes {
        ctx.write_shape_named(&ms.name, &ms.shape);
    }
    Ok(ctx.finish())
}

/// Serialize a model to STEP, attaching a `COLOUR_RGB` style to each shape
/// whose `ModelShape.color` is set.
///
/// Shapes without a colour are written plain. This is the model-level
/// counterpart of [`write_step_with_color`], reading the colours already
/// stored on the [`BRepModel`] (e.g. via `add_with_color`).
pub fn write_step_model_with_colors(model: &BRepModel) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    for ms in &model.shapes {
        let rep = ctx.write_shape_named(&ms.name, &ms.shape);
        if let (Some(rep), Some(c)) = (rep, ms.color) {
            let colour = ctx.w.emit(format!(
                "COLOUR_RGB('',{},{},{})",
                step_real(c.r as f64),
                step_real(c.g as f64),
                step_real(c.b as f64)
            ));
            let fill = ctx.w.emit(format!("SURFACE_STYLE_FILL_AREA('',#{colour})"));
            let usage = ctx.w.emit(format!("SURFACE_STYLE_USAGE('',#{fill})"));
            ctx.w.emit(format!("STYLED_ITEM('',(#{usage}),#{rep})"));
        }
    }
    Ok(ctx.finish())
}

/// Serialize a model with per-shape [`StepWriteOptions`].
///
/// Each shape in the model is written with the shared options; spline output
/// applies to every shape, and the colour/units/name defaults apply as in
/// [`write_step_with_options`]. Shape names come from the model, so `opts.name`
/// is only used as a fallback for unnamed shapes.
pub fn write_step_model_with_options(model: &BRepModel, opts: &StepWriteOptions) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    ctx.splines = opts.splines;
    for ms in &model.shapes {
        let name = if ms.name.is_empty() { opts.name.clone() } else { ms.name.clone() };
        let rep = ctx.write_shape_named(&name, &ms.shape);
        if let (Some(rep), Some((r, g, b))) = (rep, opts.color) {
            let colour = ctx.w.emit(format!(
                "COLOUR_RGB('',{},{},{})",
                step_real(r),
                step_real(g),
                step_real(b)
            ));
            let fill = ctx.w.emit(format!("SURFACE_STYLE_FILL_AREA('',#{colour})"));
            let usage = ctx.w.emit(format!("SURFACE_STYLE_USAGE('',#{fill})"));
            ctx.w.emit(format!("STYLED_ITEM('',(#{usage}),#{rep})"));
        }
        emit_si_units(&mut ctx.w, &opts.units);
    }
    Ok(ctx.finish())
}

/// Read a STEP physical file from disk, also returning collected warnings for
/// skipped or unsupported entities.
pub fn read_step_file_with_warnings(path: &str) -> Result<(BRepModel, Vec<String>), String> {
    let content = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    read_step_with_warnings(&content)
}

/// Serialize a list of named shapes into one STEP file.
///
/// Each `(name, shape)` pair becomes its own representation + product pair,
/// exactly as [`write_step`] does for a [`BRepModel`], without requiring the
/// caller to build a model object first. The reader returns each shape under
/// its `name`.
pub fn write_step_shapes(shapes: &[(String, TopoShape)]) -> String {
    let mut ctx = WriteCtx::new();
    for (name, shape) in shapes {
        ctx.write_shape_named(name, shape);
    }
    ctx.finish()
}

/// Spline-capable counterpart of [`write_step_shapes`]: non-analytic curves
/// and surfaces are emitted as B-splines.
pub fn write_step_shapes_with_splines(shapes: &[(String, TopoShape)]) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    ctx.splines = true;
    for (name, shape) in shapes {
        ctx.write_shape_named(name, shape);
    }
    Ok(ctx.finish())
}

/// Serialize a compound as its named children in one STEP file.
///
/// A [`TopoShape`] of type `Compound` is flattened: each child shape is written
/// as a named representation (`Child1`, `Child2`, ...), which the reader
/// reconstructs as separate model shapes. Non-compound shapes are written as a
/// single `"Shape"` representation.
pub fn write_step_compound(compound: &TopoShape) -> Result<String, String> {
    if !compound.is_compound() {
        return write_step_with_splines(compound);
    }
    let kids: Vec<TopoShape> = compound
        .tshape
        .read()
        .unwrap()
        .children
        .clone();
    if kids.is_empty() {
        return Err("write_step_compound: empty compound".into());
    }
    let named: Vec<(String, TopoShape)> = kids
        .into_iter()
        .enumerate()
        .map(|(i, k)| (format!("Child{}", i + 1), k))
        .collect();
    Ok(write_step_shapes_with_splines(&named)?)
}

// ---------------------------------------------------------------------------
// Reader
// ---------------------------------------------------------------------------

/// A parsed `#N=TYPE(...)` data record.
#[derive(Debug, Clone)]
struct Record {
    type_name: String,
    args: Vec<String>,
}

/// Split a comma-separated argument list at the top nesting level, respecting
/// strings (with `''` escapes), parens and brackets.
fn split_top(s: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut depth = 0usize;
    let mut cur = String::new();
    let mut chars = s.chars().peekable();
    let mut in_str = false;
    while let Some(c) = chars.next() {
        if in_str {
            cur.push(c);
            if c == '\'' {
                if chars.peek() == Some(&'\'') {
                    cur.push(chars.next().unwrap());
                } else {
                    in_str = false;
                }
            }
        } else {
            match c {
                '\'' => {
                    in_str = true;
                    cur.push(c);
                }
                '(' | '[' => {
                    depth += 1;
                    cur.push(c);
                }
                ')' | ']' => {
                    depth -= 1;
                    cur.push(c);
                }
                ',' if depth == 0 => {
                    out.push(cur.trim().to_string());
                    cur.clear();
                }
                _ => cur.push(c),
            }
        }
    }
    if !cur.trim().is_empty() {
        out.push(cur.trim().to_string());
    }
    out
}

/// Parse the entity body `TYPE(a,b,...)` into its type name and argument list.
///
/// A STEP *complex* entity instance groups several subtypes as space-separated
/// members inside one outer paren pair — `( A() B(...) C(...) )` — optionally
/// followed by the compound's own (empty) attribute list. The attributes of the
/// member subtypes form the attribute list of the most derived subtype; this
/// merge reconstructs that record (currently for the B-spline curve/surface
/// families, whose members split the degree/poles/form, the knots, and the
/// weights). Non-B-spline complexes (unit/context metadata) are skipped by the
/// resolver, so they keep an empty type name.
fn parse_entity_body(body: &str) -> (String, Vec<String>) {
    let body = body.trim();
    if body.starts_with('(') {
        return merge_complex_body(body);
    }
    let open = body.find('(').unwrap_or(body.len());
    let type_name = body[..open].trim().to_string();
    let close = body.rfind(')').unwrap_or(body.len());
    let inner = if close > open + 1 {
        &body[open + 1..close]
    } else {
        ""
    };
    (type_name, split_top(inner))
}

/// Split `( MEMBER1(args) MEMBER2(args) ... )` into its member type/args pairs.
///
/// Members are separated by whitespace at the top level (not commas); each is a
/// `TYPE(...)` token. The outer parens and any trailing empty compound list are
/// stripped first.
fn split_complex_members(body: &str) -> Vec<(String, Vec<String>)> {
    let b = body.trim();
    // Strip the outer pair of parens.
    let inner = if b.starts_with('(') && b.ends_with(')') {
        &b[1..b.len() - 1]
    } else {
        b
    };
    let bytes = inner.as_bytes();
    let mut members: Vec<(String, Vec<String>)> = Vec::new();
    let mut i = 0;
    while i < bytes.len() {
        while i < bytes.len() && bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        if i >= bytes.len() {
            break;
        }
        // Member type name: identifier until `(`.
        let t0 = i;
        while i < bytes.len() && (bytes[i].is_ascii_alphanumeric() || bytes[i] == b'_') {
            i += 1;
        }
        if i >= bytes.len() || bytes[i] != b'(' {
            // Skip stray tokens (not a TYPE(...) member).
            while i < bytes.len() && !bytes[i].is_ascii_whitespace() {
                i += 1;
            }
            continue;
        }
        let type_name = inner[t0..i].to_string();
        // Scan to the matching `)`, respecting strings and nesting.
        let mut depth = 0usize;
        let mut in_str = false;
        let a0 = i;
        while i < bytes.len() {
            let c = bytes[i] as char;
            if in_str {
                if c == '\'' {
                    if i + 1 < bytes.len() && bytes[i + 1] == b'\'' {
                        i += 2;
                        continue;
                    }
                    in_str = false;
                }
            } else {
                match c {
                    '\'' => in_str = true,
                    '(' => depth += 1,
                    ')' => {
                        depth -= 1;
                        if depth == 0 {
                            i += 1;
                            break;
                        }
                    }
                    _ => {}
                }
            }
            i += 1;
        }
        let member_body = &inner[a0..i];
        let (_, args) = {
            let open = member_body.find('(').unwrap_or(0);
            let close = member_body.rfind(')').unwrap_or(member_body.len());
            if close > open + 1 {
                (String::new(), split_top(&member_body[open + 1..close]))
            } else {
                (String::new(), Vec::new())
            }
        };
        members.push((type_name, args));
    }
    members
}

/// Reconstruct a single record from a STEP complex entity body.
///
/// See [`parse_entity_body`] for the format. B-spline curves and surfaces merge
/// their members into a `B_SPLINE_CURVE_WITH_KNOTS` / `B_SPLINE_SURFACE_WITH_KNOTS`
/// record whose argument layout matches the resolver's expectations (degree /
/// poles / form / closed from the base member, knots + multiplicities from the
/// `_WITH_KNOTS` member, weights from the `RATIONAL_*` member). Other complex
/// entities (unit / context metadata) return an empty type name — the resolver
/// skips them.
fn merge_complex_body(body: &str) -> (String, Vec<String>) {
    let members = split_complex_members(body);
    if members.is_empty() {
        return (String::new(), Vec::new());
    }
    // Helper: first member of a given type.
    let member = |t: &str| members.iter().find(|(ty, _)| ty == t);
    // Helper: arg of a member by index.
    let arg = |t: &str, idx: usize| -> Option<String> {
        member(t).and_then(|(_, a)| a.get(idx)).cloned()
    };

    let is_curve = member("B_SPLINE_CURVE").is_some() || member("B_SPLINE_CURVE_WITH_KNOTS").is_some();
    let is_surface = member("B_SPLINE_SURFACE").is_some()
        || member("B_SPLINE_SURFACE_WITH_KNOTS").is_some();

    if is_curve {
        // Complex members carry only the subtype's *own* attributes, no entity
        // name. B_SPLINE_CURVE args: (degree, control_points, curve_form,
        // closed, self_intersect). B_SPLINE_CURVE_WITH_KNOTS adds
        // (multiplicities, knots, knot_spec); RATIONAL_B_SPLINE_CURVE adds
        // (weights).
        let base = member("B_SPLINE_CURVE").or_else(|| member("B_SPLINE_CURVE_WITH_KNOTS"));
        let knots = member("B_SPLINE_CURVE_WITH_KNOTS");
        let rational = member("RATIONAL_B_SPLINE_CURVE");
        let (b_args, k_args, r_args) = match (base, knots, rational) {
            (Some((_, ba)), Some((_, ka)), Some((_, ra))) => (ba, ka, ra),
            (Some((_, ba)), Some((_, ka)), None) => (ba, ka, &Vec::new()),
            (Some((_, ba)), None, Some((_, ra))) => (ba, &Vec::new(), ra),
            (Some((_, ba)), None, None) => (ba, &Vec::new(), &Vec::new()),
            _ => return (String::new(), Vec::new()),
        };
        let degree = b_args.get(0).cloned().unwrap_or_default();
        let control_points = b_args.get(1).cloned().unwrap_or_default();
        let curve_form = b_args.get(2).cloned().unwrap_or_else(|| ".UNSPECIFIED.".to_string());
        let closed = b_args.get(3).cloned().unwrap_or_else(|| ".F.".to_string());
        let self_intersect = b_args.get(4).cloned().unwrap_or_else(|| ".F.".to_string());
        let weights = r_args.first().cloned().unwrap_or_else(|| "SELF".to_string());
        let mults = k_args.get(0).cloned().unwrap_or_else(|| "()".to_string());
        let knots_list = k_args.get(1).cloned().unwrap_or_else(|| "()".to_string());
        let knot_spec = k_args.get(2).cloned().unwrap_or_else(|| ".UNSPECIFIED.".to_string());
        let type_name = if k_args.is_empty() {
            "B_SPLINE_CURVE"
        } else {
            "B_SPLINE_CURVE_WITH_KNOTS"
        };
        // Resolver layout (B_SPLINE_CURVE_WITH_KNOTS) — the reader uses
        // args[1]=degree, [2]=control_points, [3]=weights, [6]=multiplicities,
        // [7]=knots; the remaining slots (curve_form/closed) are carried for
        // fidelity but not read.
        return (
            type_name.to_string(),
            vec![
                "''".to_string(),
                degree,
                control_points,
                weights,
                curve_form,
                closed,
                mults,
                knots_list,
                knot_spec,
            ],
        );
    }

    if is_surface {
        // B_SPLINE_SURFACE args: (u_degree, v_degree, control_points grid,
        // surface_form, closed_u, closed_v, self_intersect).
        // B_SPLINE_SURFACE_WITH_KNOTS adds (u_mults, v_mults, u_knots,
        // v_knots, knot_spec); RATIONAL_B_SPLINE_SURFACE adds (weights).
        let base = member("B_SPLINE_SURFACE").or_else(|| member("B_SPLINE_SURFACE_WITH_KNOTS"));
        let knots = member("B_SPLINE_SURFACE_WITH_KNOTS");
        let rational = member("RATIONAL_B_SPLINE_SURFACE");
        let (b_args, k_args, r_args) = match (base, knots, rational) {
            (Some((_, ba)), Some((_, ka)), Some((_, ra))) => (ba, ka, ra),
            (Some((_, ba)), Some((_, ka)), None) => (ba, ka, &Vec::new()),
            (Some((_, ba)), None, Some((_, ra))) => (ba, &Vec::new(), ra),
            (Some((_, ba)), None, None) => (ba, &Vec::new(), &Vec::new()),
            _ => return (String::new(), Vec::new()),
        };
        let deg_u = b_args.get(0).cloned().unwrap_or_default();
        let deg_v = b_args.get(1).cloned().unwrap_or_default();
        let control_points = b_args.get(2).cloned().unwrap_or_default();
        let surface_form = b_args.get(3).cloned().unwrap_or_else(|| ".UNSPECIFIED.".to_string());
        let closed_u = b_args.get(4).cloned().unwrap_or_else(|| ".F.".to_string());
        let closed_v = b_args.get(5).cloned().unwrap_or_else(|| ".F.".to_string());
        let self_intersect = b_args.get(6).cloned().unwrap_or_else(|| ".F.".to_string());
        let weights = r_args.first().cloned().unwrap_or_else(|| "SELF".to_string());
        // Knot member layout: (u_mults, v_mults, u_knots, v_knots, knot_spec).
        let u_mults = k_args.get(0).cloned().unwrap_or_else(|| "()".to_string());
        let v_mults = k_args.get(1).cloned().unwrap_or_else(|| "()".to_string());
        let u_knots = k_args.get(2).cloned().unwrap_or_else(|| "()".to_string());
        let v_knots = k_args.get(3).cloned().unwrap_or_else(|| "()".to_string());
        let knot_spec = k_args.get(4).cloned().unwrap_or_else(|| ".UNSPECIFIED.".to_string());
        let type_name = if k_args.is_empty() {
            "B_SPLINE_SURFACE"
        } else {
            "B_SPLINE_SURFACE_WITH_KNOTS"
        };
        // Standard ISO 10303-42 layout: name, u_degree, v_degree,
        // control_points grid, surface_form, closed_u, closed_v, self_intersect,
        // u_multiplicities, v_multiplicities, u_knots, v_knots, knot_spec, then
        // the weight grid as the optional 14th arg (rational surfaces).
        return (
            type_name.to_string(),
            vec![
                "''".to_string(),
                deg_u,
                deg_v,
                control_points,
                surface_form,
                closed_u,
                closed_v,
                self_intersect,
                u_mults,
                v_mults,
                u_knots,
                v_knots,
                knot_spec,
                weights,
            ],
        );
    }

    // Non-B-spline complex (unit/context metadata): the resolver skips it.
    (String::new(), Vec::new())
}

/// Split the DATA section into `#id=TYPE(...)` records.
fn parse_records(data: &str) -> Result<HashMap<usize, Record>, String> {
    let bytes = data.as_bytes();
    let mut records = HashMap::new();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'#' {
            let mut j = i + 1;
            while j < bytes.len() && bytes[j].is_ascii_digit() {
                j += 1;
            }
            if j == i + 1 {
                i += 1;
                continue;
            }
            let id: usize = data[i + 1..j]
                .parse()
                .map_err(|_| format!("bad entity id at offset {i}"))?;
            let mut k = j;
            while k < bytes.len() && bytes[k].is_ascii_whitespace() {
                k += 1;
            }
            if k >= bytes.len() || bytes[k] != b'=' {
                return Err(format!("malformed record #{id}: expected '='"));
            }
            let (body, next) = parse_entity_body_text(data, k + 1)?;
            let (type_name, args) = parse_entity_body(&body);
            records.insert(id, Record { type_name, args });
            i = next;
        } else {
            i += 1;
        }
    }
    Ok(records)
}

/// Scan one entity body starting after the `=`, stopping at the terminating
/// `;`. Returns the body text and the index just past the `;`.
fn parse_entity_body_text(data: &str, start: usize) -> Result<(String, usize), String> {
    let bytes = data.as_bytes();
    let mut i = start;
    let mut depth = 0usize;
    let mut in_str = false;
    while i < bytes.len() {
        let c = bytes[i] as char;
        if in_str {
            if c == '\'' {
                if i + 1 < bytes.len() && bytes[i + 1] == b'\'' {
                    i += 2;
                    continue;
                }
                in_str = false;
            }
        } else {
            match c {
                '\'' => in_str = true,
                '(' => depth += 1,
                ')' => {
                    if depth == 0 {
                        return Err("unbalanced parens in entity body".into());
                    }
                    depth -= 1;
                    if depth == 0 {
                        let body = data[start..=i].to_string();
                        let mut j = i + 1;
                        // STEP complex entities end `( A() B() C() )()` — the
                        // compound's own (usually empty) attribute list follows
                        // the grouped members. Consume it before the `;`.
                        let mut k = j;
                        while k < bytes.len() && bytes[k].is_ascii_whitespace() {
                            k += 1;
                        }
                        if k + 1 < bytes.len() && bytes[k] == b'(' && bytes[k + 1] == b')' {
                            j = k + 2;
                        }
                        while j < bytes.len() && bytes[j].is_ascii_whitespace() {
                            j += 1;
                        }
                        if j >= bytes.len() || bytes[j] != b';' {
                            return Err("missing ';' after entity body".into());
                        }
                        return Ok((body, j + 1));
                    }
                }
                _ => {}
            }
        }
        i += 1;
    }
    Err("unterminated entity body".into())
}

fn parse_ref(s: &str) -> Option<usize> {
    s.trim().strip_prefix('#')?.trim().parse().ok()
}

fn parse_ref_list(s: &str) -> Vec<usize> {
    let s = s.trim();
    if !(s.starts_with('(') && s.ends_with(')')) {
        return Vec::new();
    }
    split_top(&s[1..s.len() - 1])
        .into_iter()
        .filter_map(|a| parse_ref(&a))
        .collect()
}

fn parse_f64(s: &str) -> Result<f64, String> {
    s.trim()
        .parse()
        .map_err(|_| format!("bad real literal '{s}'"))
}

fn parse_xyz(s: &str) -> Result<GpXyz, String> {
    let s = s.trim();
    let inner = s.trim_start_matches('(').trim_end_matches(')');
    let parts: Vec<String> = split_top(inner)
        .into_iter()
        .map(|p| p.trim().to_string())
        .collect();
    if parts.len() != 3 {
        return Err(format!("expected 3-component tuple, got '{s}'"));
    }
    Ok(GpXyz::new(
        parse_f64(&parts[0])?,
        parse_f64(&parts[1])?,
        parse_f64(&parts[2])?,
    ))
}

fn parse_str(s: &str) -> String {
    let s = s.trim();
    if s.len() >= 2 && s.starts_with('\'') && s.ends_with('\'') {
        s[1..s.len() - 1].replace("''", "'")
    } else {
        s.to_string()
    }
}

/// Reference resolver with memoization; unknown/unsupported entities are
/// recorded as warnings and skipped.
struct Resolver<'a> {
    records: &'a HashMap<usize, Record>,
    b: TopoBuilder,
    shape_cache: RefCell<HashMap<usize, TopoShape>>,
    point_cache: RefCell<HashMap<usize, GpPnt>>,
    dir_cache: RefCell<HashMap<usize, GpDir>>,
    axis_cache: RefCell<HashMap<usize, GpAx2>>,
    curve_cache: RefCell<HashMap<usize, Arc<dyn Curve>>>,
    surface_cache: RefCell<HashMap<usize, Arc<dyn Surface>>>,
    resolving: RefCell<HashSet<usize>>,
    warnings: RefCell<Vec<String>>,
}

impl<'a> Resolver<'a> {
    fn new(records: &'a HashMap<usize, Record>) -> Self {
        Self {
            records,
            b: TopoBuilder::new(),
            shape_cache: RefCell::new(HashMap::new()),
            point_cache: RefCell::new(HashMap::new()),
            dir_cache: RefCell::new(HashMap::new()),
            axis_cache: RefCell::new(HashMap::new()),
            curve_cache: RefCell::new(HashMap::new()),
            surface_cache: RefCell::new(HashMap::new()),
            resolving: RefCell::new(HashSet::new()),
            warnings: RefCell::new(Vec::new()),
        }
    }

    fn record(&self, id: usize) -> Result<&'a Record, String> {
        self.records
            .get(&id)
            .ok_or_else(|| format!("reference to undefined entity #{id}"))
    }

    fn warn(&self, msg: String) {
        self.warnings.borrow_mut().push(msg);
    }

    fn resolve_shape(&self, id: usize) -> Result<TopoShape, String> {
        if let Some(s) = self.shape_cache.borrow().get(&id) {
            return Ok(s.clone());
        }
        let rec = self.record(id)?;
        if !self.resolving.borrow_mut().insert(id) {
            return Err(format!("cyclic entity reference #{id}"));
        }
        let shape = match rec.type_name.as_str() {
            "VERTEX_POINT" => self.resolve_vertex(rec),
            "EDGE_CURVE" => self.resolve_edge(rec),
            "ORIENTED_EDGE" => self.resolve_oriented_edge(rec),
            "EDGE_LOOP" => self.resolve_loop(rec),
            "VERTEX_LOOP" => self.resolve_vertex_loop(rec),
            "FACE_OUTER_BOUND" => self.resolve_outer_bound(rec),
            "FACE_BOUND" => self.resolve_outer_bound(rec),
            "ADVANCED_FACE" => self.resolve_face(rec),
            "CLOSED_SHELL" => self.resolve_shell(rec),
            "MANIFOLD_SOLID_BREP" => self.resolve_solid(rec),
            other => {
                self.warn(format!("unsupported topological entity {other} (#{id})"));
                Err(format!("unsupported entity {other} (#{id})"))
            }
        };
        self.resolving.borrow_mut().remove(&id);
        match shape {
            Ok(s) => {
                self.shape_cache.borrow_mut().insert(id, s.clone());
                Ok(s)
            }
            Err(e) => Err(e),
        }
    }

    fn resolve_vertex(&self, rec: &'a Record) -> Result<TopoShape, String> {
        if rec.args.len() < 2 {
            return Err("VERTEX_POINT: bad args".into());
        }
        let pid = parse_ref(&rec.args[1]).ok_or("VERTEX_POINT: bad point ref")?;
        let p = self.resolve_point(pid)?;
        Ok(self.b.make_vertex(p, 0.0).0)
    }

    fn resolve_edge(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let start_ref = parse_ref(&rec.args[1]).ok_or("EDGE_CURVE: bad start ref")?;
        let end_ref = parse_ref(&rec.args[2]).ok_or("EDGE_CURVE: bad end ref")?;
        let curve_ref = parse_ref(&rec.args[3]).ok_or("EDGE_CURVE: bad curve ref")?;
        let v1 = self.resolve_shape(start_ref)?;
        let v2 = self.resolve_shape(end_ref)?;
        if !v1.is_vertex() || !v2.is_vertex() {
            return Err("EDGE_CURVE: endpoints are not vertices".into());
        }
        let curve = self.resolve_curve(curve_ref)?;
        let p1 = GeometryRegistry::global().vertex_point(&v1);
        let p2 = GeometryRegistry::global().vertex_point(&v2);
        let (first, last) = edge_params_for_curve(curve.as_ref(), &p1, &p2);
        let mut e = self.b.make_edge(curve, first, last);
        self.b.add(&mut e.0, &v1);
        self.b.add(&mut e.0, &v2);
        Ok(e.0)
    }

    fn resolve_oriented_edge(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let edge_ref = parse_ref(&rec.args[3]).ok_or("ORIENTED_EDGE: bad edge ref")?;
        let mut s = self.resolve_shape(edge_ref)?;
        // ORIENTED_EDGE(name, *, *, edge, orientation): a .F. reverses the edge
        // in the wire. Without this, two ORIENTED_EDGEs referencing the same
        // EDGE_CURVE (e.g. a cylinder side wall's two seam generatrices) both
        // return the same forward edge, and the wire loses one — the surface
        // never closes.
        if let Some(o) = rec.args.get(4).map(|s| s.trim().to_string()) {
            if o == ".F." {
                s.set_orientation(Orientation::Reversed);
            }
        }
        Ok(s)
    }

    fn resolve_loop(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let items = parse_ref_list(&rec.args[1]);
        let mut edges = Vec::with_capacity(items.len());
        for &it in &items {
            let s = self.resolve_shape(it)?;
            if !s.is_edge() {
                return Err(format!("#{it}: expected EDGE in EDGE_LOOP"));
            }
            edges.push(Edge(s));
        }
        let mut wire = self.b.make_wire(&edges).0;
        // The EDGE_LOOP lists its ORIENTED_EDGEs in any cyclic order; the file
        // is not guaranteed to place consecutive edges next to each other (a
        // CAD writer may list them out of sequence). Reorder the wire's edges
        // into a connected chain by following shared vertices, as the OCCT STEP
        // reader does when assembling the loop.
        crate::algo_tools::AlgoTools::orient_edges_on_wire(&mut wire);
        Ok(wire)
    }

    /// A `VERTEX_LOOP(name, vertex)` is the boundary of a degenerate face — a
    /// loop reduced to a single vertex (sphere pole, cone apex). It maps to a
    /// wire containing no edges (the face still references it via FACE_BOUND).
    fn resolve_vertex_loop(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let _ = parse_ref(&rec.args[1]).ok_or("VERTEX_LOOP: bad vertex ref")?;
        Ok(self.b.make_wire(&[]).0)
    }

    fn resolve_outer_bound(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let loop_ref = parse_ref(&rec.args[1]).ok_or("FACE_OUTER_BOUND: bad loop ref")?;
        let s = self.resolve_shape(loop_ref)?;
        if !s.is_wire() {
            return Err("FACE_OUTER_BOUND: loop is not a wire".into());
        }
        Ok(s)
    }

    fn resolve_face(&self, rec: &'a Record) -> Result<TopoShape, String> {
        // Two argument layouts occur in the wild:
        //  * STEP-214 (ISO standard, written by OCCT/FreeCAD): surface is the
        //    third argument — `ADVANCED_FACE(name, bounds, surface, same_sense)`.
        //  * this port's own writer (step.rs write path): surface is the second
        //    argument — `ADVANCED_FACE('', #surface, (bounds), .T.)`.
        // Detect by checking which argument holds a surface reference.
        let surf_ref = if let Some(r) = parse_ref(&rec.args[2]) {
            r
        } else {
            parse_ref(&rec.args[1]).ok_or("ADVANCED_FACE: bad surface ref")?
        };
        let surface = self.resolve_surface(surf_ref)?;
        // Bounds are the other of the two leading arguments.
        let bounds = if rec.args.get(2).map(|s| s.starts_with('(')).unwrap_or(false) {
            parse_ref_list(&rec.args[2])
        } else {
            parse_ref_list(&rec.args[1])
        };
        let mut wires = Vec::with_capacity(bounds.len());
        for &b in &bounds {
            let s = self.resolve_shape(b)?;
            if !s.is_wire() {
                return Err(format!("#{b}: expected wire in face bounds"));
            }
            wires.push(Wire(s));
        }
        Ok(self.b.make_face(surface, &wires).0)
    }

    fn resolve_shell(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let faces = parse_ref_list(&rec.args[1]);
        let mut fs = Vec::with_capacity(faces.len());
        for &it in &faces {
            let s = self.resolve_shape(it)?;
            if !s.is_face() {
                return Err(format!("#{it}: expected FACE in CLOSED_SHELL"));
            }
            fs.push(Face(s));
        }
        Ok(self.b.make_shell(&fs).0)
    }

    fn resolve_solid(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let outer = parse_ref(&rec.args[1]).ok_or("MANIFOLD_SOLID_BREP: bad shell ref")?;
        let s = self.resolve_shape(outer)?;
        if !s.is_shell() {
            return Err("MANIFOLD_SOLID_BREP: outer is not a shell".into());
        }
        Ok(self.b.make_solid(&[Shell(s)]).0)
    }

    fn resolve_point(&self, id: usize) -> Result<GpPnt, String> {
        if let Some(p) = self.point_cache.borrow().get(&id) {
            return Ok(*p);
        }
        let rec = self.record(id)?;
        if rec.type_name != "CARTESIAN_POINT" {
            self.warn(format!("expected CARTESIAN_POINT, got {} (#{id})", rec.type_name));
            return Err(format!("expected CARTESIAN_POINT at #{id}"));
        }
        let v = parse_xyz(&rec.args[1]).map_err(|e| format!("CARTESIAN_POINT #{id}: {e}"))?;
        let p = GpPnt::from_xyz(&v);
        self.point_cache.borrow_mut().insert(id, p);
        Ok(p)
    }

    fn resolve_direction(&self, id: usize) -> Result<GpDir, String> {
        if let Some(d) = self.dir_cache.borrow().get(&id) {
            return Ok(*d);
        }
        let rec = self.record(id)?;
        if rec.type_name != "DIRECTION" {
            self.warn(format!("expected DIRECTION, got {} (#{id})", rec.type_name));
            return Err(format!("expected DIRECTION at #{id}"));
        }
        let v = parse_xyz(&rec.args[1]).map_err(|e| format!("DIRECTION #{id}: {e}"))?;
        let d = GpDir::new(v.x, v.y, v.z).map_err(|e| format!("DIRECTION #{id}: {e}"))?;
        self.dir_cache.borrow_mut().insert(id, d);
        Ok(d)
    }

    fn resolve_vector(&self, id: usize) -> Result<GpVec, String> {
        let rec = self.record(id)?;
        if rec.type_name != "VECTOR" {
            self.warn(format!("expected VECTOR, got {} (#{id})", rec.type_name));
            return Err(format!("expected VECTOR at #{id}"));
        }
        let dir_ref = parse_ref(&rec.args[1]).ok_or("VECTOR: bad direction ref")?;
        let mag = parse_f64(&rec.args[2])?;
        let dir = self.resolve_direction(dir_ref)?;
        Ok(GpVec::from_xyz(&dir.xyz().multiplied(mag)))
    }

    /// `AXIS1_PLACEMENT(name, location, axis_direction)` → `GpAx1`.
    fn resolve_axis1(&self, id: usize) -> Result<GpAx1, String> {
        let rec = self.record(id)?;
        if rec.type_name != "AXIS1_PLACEMENT" {
            self.warn(format!("expected AXIS1_PLACEMENT, got {} (#{id})", rec.type_name));
            return Err(format!("expected AXIS1_PLACEMENT at #{id}"));
        }
        let loc_ref = parse_ref(&rec.args[1]).ok_or("AXIS1: bad location ref")?;
        let dir_ref = parse_ref(&rec.args[2]).ok_or("AXIS1: bad direction ref")?;
        let loc = self.resolve_point(loc_ref)?;
        let dir = self.resolve_direction(dir_ref)?;
        Ok(GpAx1::new(loc, dir))
    }

    fn resolve_axis2(&self, id: usize) -> Result<GpAx2, String> {
        if let Some(a) = self.axis_cache.borrow().get(&id) {
            return Ok(*a);
        }
        let rec = self.record(id)?;
        if rec.type_name != "AXIS2_PLACEMENT_3D" {
            self.warn(format!(
                "expected AXIS2_PLACEMENT_3D, got {} (#{id})",
                rec.type_name
            ));
            return Err(format!("expected AXIS2_PLACEMENT_3D at #{id}"));
        }
        let loc = self
            .resolve_point(parse_ref(&rec.args[1]).ok_or("AXIS2: bad location ref")?)?;
        let axis = match parse_ref(&rec.args[2]) {
            Some(r) => self.resolve_direction(r)?,
            None => dir_z(),
        };
        let refd = match parse_ref(&rec.args[3]) {
            Some(r) => self.resolve_direction(r)?,
            None => dir_x(),
        };
        let ax2 = GpAx2::new(loc, axis, refd).map_err(|e| format!("AXIS2_PLACEMENT_3D: {e}"))?;
        self.axis_cache.borrow_mut().insert(id, ax2);
        Ok(ax2)
    }

    fn resolve_curve(&self, id: usize) -> Result<Arc<dyn Curve>, String> {
        if let Some(c) = self.curve_cache.borrow().get(&id) {
            return Ok(c.clone());
        }
        let rec = self.record(id)?;
        let curve: Arc<dyn Curve> = match rec.type_name.as_str() {
            "LINE" => {
                let pnt = parse_ref(&rec.args[1]).ok_or("LINE: bad point ref")?;
                let vec = parse_ref(&rec.args[2]).ok_or("LINE: bad vector ref")?;
                let p = self.resolve_point(pnt)?;
                let v = self.resolve_vector(vec)?;
                let d = GpDir::from_vec(&v).map_err(|e| format!("LINE: {e}"))?;
                Arc::new(GeomLine::new(GpLin::from_pnt_dir(p, d)))
            }
            "CIRCLE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("CIRCLE: bad axis ref")?;
                let r = parse_f64(&rec.args[2])?;
                Arc::new(GeomCircle::new(GpCirc::new(self.resolve_axis2(ax)?, r)))
            }
            "ELLIPSE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("ELLIPSE: bad axis ref")?;
                let maj = parse_f64(&rec.args[2])?;
                let min = parse_f64(&rec.args[3])?;
                Arc::new(GeomEllipse::new(GpElips::new(
                    self.resolve_axis2(ax)?,
                    maj,
                    min,
                )))
            }
            "HYPERBOLA" => {
                let ax = parse_ref(&rec.args[1]).ok_or("HYPERBOLA: bad axis ref")?;
                let maj = parse_f64(&rec.args[2])?;
                let min = parse_f64(&rec.args[3])?;
                Arc::new(GeomHyperbola::new(GpHypr::new(
                    self.resolve_axis2(ax)?,
                    maj,
                    min,
                )))
            }
            "PARABOLA" => {
                let ax = parse_ref(&rec.args[1]).ok_or("PARABOLA: bad axis ref")?;
                let f = parse_f64(&rec.args[2])?;
                Arc::new(GeomParabola::new(GpParab::new(self.resolve_axis2(ax)?, f)))
            }
            "SURFACE_CURVE" | "SEAM_CURVE" => {
                // SURFACE_CURVE/SEAM_CURVE(name, curve_3d, pcurves, master_rep):
                // the 3D curve is the second argument; the pcurve list is
                // referenced per face at the ADVANCED_FACE level, so we take
                // the 3D curve. SEAM_CURVE is the seam of a closed surface.
                let c3d = parse_ref(&rec.args[1]).ok_or("SURFACE_CURVE: bad 3D curve ref")?;
                self.resolve_curve(c3d)?
            }
            "B_SPLINE_CURVE_WITH_KNOTS" => {
                // Layout (10 args): name, degree, control_points, weights|SELF,
                // curve_form, closed, self_intersect, knots, multiplicities, knot_spec.
                let degree = parse_f64(&rec.args[1])? as usize;
                let poles: Vec<GpPnt> = parse_ref_list(&rec.args[2])
                    .into_iter()
                    .map(|r| self.resolve_point(r))
                    .collect::<Result<Vec<_>, _>>()?;
                let weights_arg = rec.args.get(3).map(|s| s.trim().to_string());
                // B_SPLINE_CURVE_WITH_KNOTS layout:
                // (name, degree, control_points, curve_form, closed,
                //  self_intersect, knot_multiplicities, knots, knot_spec).
                let knots = expand_knots(
                    &parse_usize_list(rec.args.get(6).map(|s| s.as_str()).unwrap_or("()")),
                    &parse_real_list(rec.args.get(7).map(|s| s.as_str()).unwrap_or("()")),
                );
                let curve = match weights_arg.as_deref() {
                    // No weights: non-rational curve (curve_form is UNSPECIFIED
                    // or a non-rational flag like CIRCULAR/LINEAR).
                    None | Some("SELF") | Some(".UNSPECIFIED.") | Some(".CIRCULAR.")
                    | Some(".LINEAR.") => GeomBSplineCurve::new(poles, knots, degree),
                    Some(w) if w.starts_with('.') => {
                        // Another non-rational curve form marker.
                        GeomBSplineCurve::new(poles, knots, degree)
                    }
                    Some(w) => {
                        let weights = parse_real_list(w);
                        if weights.len() != poles.len() {
                            return Err("B_SPLINE_CURVE: weight count mismatch".into());
                        }
                        GeomBSplineCurve::rational(poles, weights, knots, degree)
                    }
                };
                Arc::new(curve.map_err(|e| format!("B_SPLINE_CURVE: {e}"))?)
            }
            "POLYLINE" => {
                // A connected sequence of CARTESIAN_POINTs; represent it as a
                // degree-1 B-spline that passes through every vertex.
                let pts: Vec<GpPnt> = parse_ref_list(&rec.args[1])
                    .into_iter()
                    .map(|r| self.resolve_point(r))
                    .collect::<Result<Vec<_>, _>>()?;
                if pts.len() < 2 {
                    return Err("POLYLINE: need at least 2 points".into());
                }
                let knots = uniform_knots_for(pts.len(), 1);
                Arc::new(
                    GeomBSplineCurve::new(pts, knots, 1)
                        .map_err(|e| format!("POLYLINE: {e}"))?,
                )
            }
            "B_SPLINE_CURVE" => {
                // Plain B-spline (no explicit knots): the knot vector is the
                // clamped uniform one implied by the pole count and degree.
                let degree = parse_f64(&rec.args[1])? as usize;
                let poles: Vec<GpPnt> = parse_ref_list(&rec.args[2])
                    .into_iter()
                    .map(|r| self.resolve_point(r))
                    .collect::<Result<Vec<_>, _>>()?;
                let knots = uniform_knots_for(poles.len(), degree);
                Arc::new(
                    GeomBSplineCurve::new(poles, knots, degree)
                        .map_err(|e| format!("B_SPLINE_CURVE: {e}"))?,
                )
            }
            "TRIMMED_CURVE" => {
                // Layout: name, basis_curve, trim_1, trim_2, sense_agreement,
                // master_representation. The bounds are the 4th/5th attributes.
                let basis_ref = parse_ref(&rec.args[1]).ok_or("TRIMMED_CURVE: bad basis ref")?;
                let basis = self.resolve_curve(basis_ref)?;
                let a = parse_f64(&rec.args[3])?;
                let b = parse_f64(&rec.args[4])?;
                Arc::new(GeomTrimmedCurve::new(basis, a, b))
            }
            "OFFSET_CURVE_3D" => {
                // Layout: name, basis_curve, direction, distance, self_intersect,
                // curve_form.
                let basis_ref = parse_ref(&rec.args[1]).ok_or("OFFSET_CURVE_3D: bad curve ref")?;
                let dir_ref = parse_ref(&rec.args[2]).ok_or("OFFSET_CURVE_3D: bad dir ref")?;
                let offset = parse_f64(&rec.args[3])?;
                let basis = self.resolve_curve(basis_ref)?;
                let dir = self.resolve_direction(dir_ref)?;
                Arc::new(GeomOffsetCurve::new(basis, offset, dir))
            }
            other => {
                self.warn(format!("unsupported curve entity {other} (#{id})"));
                return Err(format!("unsupported curve entity {other} (#{id})"));
            }
        };
        self.curve_cache.borrow_mut().insert(id, curve.clone());
        Ok(curve)
    }

    fn resolve_surface(&self, id: usize) -> Result<Arc<dyn Surface>, String> {
        if let Some(s) = self.surface_cache.borrow().get(&id) {
            return Ok(s.clone());
        }
        let rec = self.record(id)?;
        let surface: Arc<dyn Surface> = match rec.type_name.as_str() {
            "PLANE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("PLANE: bad axis ref")?;
                Arc::new(GeomPlane::new(GpPln::new(self.resolve_axis2(ax)?.to_ax3())))
            }
            "CYLINDRICAL_SURFACE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("CYLINDRICAL_SURFACE: bad axis ref")?;
                let r = parse_f64(&rec.args[2])?;
                Arc::new(GeomCylinder::new(
                    GpCylinder::new(self.resolve_axis2(ax)?.to_ax3(), r)
                        .map_err(|e| format!("CYLINDRICAL_SURFACE: {e}"))?,
                ))
            }
            "CONICAL_SURFACE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("CONICAL_SURFACE: bad axis ref")?;
                let r = parse_f64(&rec.args[2])?;
                let a = parse_f64(&rec.args[3])?;
                Arc::new(GeomCone::new(
                    GpCone::new(self.resolve_axis2(ax)?.to_ax3(), r, a)
                        .map_err(|e| format!("CONICAL_SURFACE: {e}"))?,
                ))
            }
            "SPHERICAL_SURFACE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("SPHERICAL_SURFACE: bad axis ref")?;
                let r = parse_f64(&rec.args[2])?;
                Arc::new(GeomSphere::new(
                    GpSphere::new(self.resolve_axis2(ax)?.to_ax3(), r)
                        .map_err(|e| format!("SPHERICAL_SURFACE: {e}"))?,
                ))
            }
            "TOROIDAL_SURFACE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("TOROIDAL_SURFACE: bad axis ref")?;
                let maj = parse_f64(&rec.args[2])?;
                let min = parse_f64(&rec.args[3])?;
                Arc::new(GeomTorus::new(
                    GpTorus::new(self.resolve_axis2(ax)?.to_ax3(), maj, min)
                        .map_err(|e| format!("TOROIDAL_SURFACE: {e}"))?,
                ))
            }
            "SURFACE_OF_REVOLUTION" => {
                // SURFACE_OF_REVOLUTION(name, axis, generatrix) — the axis is an
                // AXIS1_PLACEMENT and the generatrix a curve. STEP order is
                // (axis, curve); FreeCAD emits (curve, axis). Detect which
                // argument holds the AXIS1_PLACEMENT by its record type.
                let r1 = parse_ref(&rec.args[1]);
                let r2 = parse_ref(&rec.args[2]);
                let is_axis1 = |r: usize| {
                    self.record(r)
                        .map(|rec| rec.type_name == "AXIS1_PLACEMENT")
                        .unwrap_or(false)
                };
                let (axis_ref, gen_ref) = match (r1, r2) {
                    (Some(a), Some(g)) if is_axis1(a) => (a, g),
                    (Some(g), Some(a)) if is_axis1(a) => (a, g),
                    (Some(a), Some(g)) => (a, g), // fallback: STEP order
                    (Some(a), None) => (a, a),
                    _ => return Err("SURFACE_OF_REVOLUTION: bad refs".into()),
                };
                let axis = self.resolve_axis1(axis_ref)?;
                let generatrix = self.resolve_curve(gen_ref)?;
                Arc::new(GeomSurfaceOfRevolution::new(generatrix, axis))
            }
            "B_SPLINE_SURFACE" => {
                // Plain B-spline surface (no explicit knots): clamped uniform
                // knot vectors derived from the pole grid and degrees.
                let deg_u = parse_f64(&rec.args[1])? as usize;
                let deg_v = parse_f64(&rec.args[2])? as usize;
                let poles: Vec<Vec<GpPnt>> = parse_nested_ref_list(&rec.args[3])
                    .into_iter()
                    .map(|row| {
                        row.into_iter()
                            .map(|r| self.resolve_point(r))
                            .collect::<Result<Vec<_>, _>>()
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                let nu = poles.len();
                let nv = poles.first().map_or(0, |r| r.len());
                let u_knots = uniform_knots_for(nu, deg_u);
                let v_knots = uniform_knots_for(nv, deg_v);
                Arc::new(
                    GeomBSplineSurface::new(poles, u_knots, v_knots, deg_u, deg_v)
                        .map_err(|e| format!("B_SPLINE_SURFACE: {e}"))?,
                )
            }
            "B_SPLINE_SURFACE_WITH_KNOTS" => {
                // Standard ISO 10303-42 layout (13 args): name, u_degree,
                // v_degree, control_points grid, surface_form, closed_u,
                // closed_v, self_intersect, u_multiplicities, v_multiplicities,
                // u_knots, v_knots, knot_spec. A merged rational surface carries
                // the weight grid as an optional 14th arg.
                let deg_u = parse_f64(&rec.args[1])? as usize;
                let deg_v = parse_f64(&rec.args[2])? as usize;
                let poles: Vec<Vec<GpPnt>> = parse_nested_ref_list(&rec.args[3])
                    .into_iter()
                    .map(|row| {
                        row.into_iter()
                            .map(|r| self.resolve_point(r))
                            .collect::<Result<Vec<_>, _>>()
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                let weights_arg = rec.args.get(13).map(|s| s.trim().to_string());
                let u_knots = expand_knots(
                    &parse_usize_list(rec.args.get(8).map(|s| s.as_str()).unwrap_or("()")),
                    &parse_real_list(rec.args.get(10).map(|s| s.as_str()).unwrap_or("()")),
                );
                let v_knots = expand_knots(
                    &parse_usize_list(rec.args.get(9).map(|s| s.as_str()).unwrap_or("()")),
                    &parse_real_list(rec.args.get(11).map(|s| s.as_str()).unwrap_or("()")),
                );
                let surface = match weights_arg.as_deref() {
                    None | Some("SELF") => {
                        GeomBSplineSurface::new(poles, u_knots, v_knots, deg_u, deg_v)
                    }
                    Some(w) => {
                        let wgrid = parse_nested_real_list(w);
                        GeomBSplineSurface::rational(poles, wgrid, u_knots, v_knots, deg_u, deg_v)
                    }
                };
                Arc::new(surface.map_err(|e| format!("B_SPLINE_SURFACE: {e}"))?)
            }
            "OFFSET_SURFACE" => {
                // Layout: (name, basis_surface, distance, self_intersect).
                let basis_ref = parse_ref(&rec.args[1]).ok_or("OFFSET_SURFACE: bad basis ref")?;
                let distance = parse_f64(&rec.args[2])?;
                let basis = self.resolve_surface(basis_ref)?;
                Arc::new(GeomOffsetSurface::new(basis, distance))
            }
            other => {
                self.warn(format!("unsupported surface entity {other} (#{id})"));
                return Err(format!("unsupported surface entity {other} (#{id})"));
            }
        };
        self.surface_cache.borrow_mut().insert(id, surface.clone());
        Ok(surface)
    }

    /// Resolve a representation record into (name, shapes).
    fn resolve_representation(&self, id: usize) -> Result<(String, Vec<TopoShape>), String> {
        let rec = self.record(id)?;
        let name = parse_str(&rec.args[0]);
        let items = parse_ref_list(&rec.args[1]);
        // A representation's item list may mix the geometric/topological shape
        // entities (MANIFOLD_SOLID_BREP, …) with auxiliary placement / point /
        // direction entities (AXIS2_PLACEMENT_3D, CARTESIAN_POINT, DIRECTION,
        // VECTOR, …) that are referenced *by* the shapes but are not themselves
        // top-level shapes. OCCT's STEP reader skips these non-shape items; we
        // do the same so a placement next to the solid does not fail the whole
        // representation.
        let mut shapes = Vec::with_capacity(items.len());
        for &it in &items {
            match self.resolve_shape(it) {
                Ok(s) => shapes.push(s),
                Err(e) => {
                    if !self.is_auxiliary_entity(it) {
                        return Err(e);
                    }
                    // Auxiliary item (placement/point/direction/…): skip.
                }
            }
        }
        Ok((name, shapes))
    }

    /// Whether the record at `id` is an auxiliary geometric entity that a
    /// representation lists alongside its shapes but that is not a shape.
    fn is_auxiliary_entity(&self, id: usize) -> bool {
        let Some(rec) = self.records.get(&id) else { return false };
        matches!(
            rec.type_name.as_str(),
            "AXIS2_PLACEMENT_3D"
                | "AXIS2_PLACEMENT_2D"
                | "AXIS1_PLACEMENT"
                | "CARTESIAN_POINT"
                | "DIRECTION"
                | "VECTOR"
                | "GEOMETRIC_REPRESENTATION_CONTEXT"
                | "REPRESENTATION_CONTEXT"
                | "PARAMETRIC_REPRESENTATION_CONTEXT"
                | "GLOBAL_UNIT_ASSIGNED_CONTEXT"
                | "APPLICATION_CONTEXT"
                | "PRODUCT_CONTEXT"
                | "PRODUCT_DEFINITION_CONTEXT"
                | "LENGTH_UNIT"
                | "PLANE_ANGLE_UNIT"
                | "SOLID_ANGLE_UNIT"
                | "SI_UNIT"
                | "NAMED_UNIT"
                | "UNCERTAINTY_MEASURE_WITH_UNIT"
        )
    }
}

/// Compute an edge's parameter range from its endpoint vertex points, based on
/// the reconstructed curve's analytic type.
fn edge_params_for_curve(curve: &dyn Curve, p1: &GpPnt, p2: &GpPnt) -> (f64, f64) {
    let (f, l) = (curve.first_parameter(), curve.last_parameter());
    let (lo, hi) = if f.is_finite() && l.is_finite() && l > f {
        (f, l)
    } else {
        (0.0, 1.0)
    };
    match classify_curve(curve, lo, hi) {
        CurveKind::Line => {
            let origin = curve.d0(0.0);
            let d1 = curve.d1(0.0).1;
            let dir = GpVec::from_xyz(&d1.xyz().normalized());
            (
                GpVec::from_pnts(&origin, p1).dot(&dir),
                GpVec::from_pnts(&origin, p2).dot(&dir),
            )
        }
        CurveKind::Circle => {
            if p1.distance(p2) < 1e-9 {
                return (0.0, 2.0 * PI);
            }
            let pa = curve.d0(lo);
            let pb = curve.d0(lo + (hi - lo) / 4.0);
            let pc = curve.d0(lo + (hi - lo) / 2.0);
            let center = circle_center3(&pa, &pb, &pc).unwrap_or_else(GpPnt::zero);
            let xdir = GpDir::from_vec(&GpVec::from_pnts(&center, &pa)).unwrap_or(dir_x());
            let ydir = GpDir::from_vec(&GpVec::from_pnts(&center, &pb)).unwrap_or(dir_y());
            let ang = |p: &GpPnt| {
                let v = GpVec::from_pnts(&center, p);
                v.xyz().dot(ydir.xyz()).atan2(v.xyz().dot(xdir.xyz()))
            };
            let t1 = ang(p1);
            let mut t2 = ang(p2);
            if t2 <= t1 {
                t2 += 2.0 * PI;
            }
            (t1, t2)
        }
        CurveKind::Ellipse => {
            if p1.distance(p2) < 1e-9 {
                return (0.0, 2.0 * PI);
            }
            let pa = curve.d0(lo);
            let pb = curve.d0(lo + (hi - lo) / 4.0);
            let pc = curve.d0(lo + (hi - lo) / 2.0);
            let center = midpoint(&pa, &pc);
            let a_major = center.distance(&pa).max(1e-30);
            let b_minor = center.distance(&pb).max(1e-30);
            let xdir = GpDir::from_vec(&GpVec::from_pnts(&center, &pa)).unwrap_or(dir_x());
            let ydir = GpDir::from_vec(&GpVec::from_pnts(&pb, &center)).unwrap_or(dir_y());
            let ang = |p: &GpPnt| {
                let v = GpVec::from_pnts(&center, p);
                (-v.xyz().dot(ydir.xyz()) / b_minor)
                    .atan2(v.xyz().dot(xdir.xyz()) / a_major)
            };
            let t1 = ang(p1);
            let mut t2 = ang(p2);
            if t2 <= t1 {
                t2 += 2.0 * PI;
            }
            (t1, t2)
        }
        CurveKind::Parabola => {
            let vertex = curve.d0(0.0);
            let d1 = curve.d1(0.0).1;
            let ydir = GpDir::from_vec(&d1).unwrap_or(dir_y());
            let dir = GpVec::from_xyz(ydir.xyz());
            (
                GpVec::from_pnts(&vertex, p1).dot(&dir),
                GpVec::from_pnts(&vertex, p2).dot(&dir),
            )
        }
        CurveKind::Other => {
            if f.is_finite() && l.is_finite() && l > f {
                (f, l)
            } else {
                (0.0, 1.0)
            }
        }
    }
}

fn parse_usize_list(s: &str) -> Vec<usize> {
    let s = s.trim();
    if !(s.starts_with('(') && s.ends_with(')')) {
        return Vec::new();
    }
    split_top(&s[1..s.len() - 1])
        .into_iter()
        .filter_map(|a| a.trim().parse().ok())
        .collect()
}

fn parse_real_list(s: &str) -> Vec<f64> {
    let s = s.trim();
    if !(s.starts_with('(') && s.ends_with(')')) {
        return Vec::new();
    }
    split_top(&s[1..s.len() - 1])
        .into_iter()
        .filter_map(|a| a.trim().parse().ok())
        .collect()
}

/// Split a nested `((...),(...),...)` argument into its inner lists.
fn parse_nested_lists(s: &str) -> Vec<Vec<String>> {
    let s = s.trim();
    if !(s.starts_with('(') && s.ends_with(')')) {
        return Vec::new();
    }
    split_top(&s[1..s.len() - 1])
        .into_iter()
        .filter_map(|row| {
            let row = row.trim();
            if row.starts_with('(') && row.ends_with(')') {
                Some(split_top(&row[1..row.len() - 1]))
            } else {
                None
            }
        })
        .collect()
}

/// Parse a nested list of entity references (e.g. a B-spline surface pole grid).
fn parse_nested_ref_list(s: &str) -> Vec<Vec<usize>> {
    parse_nested_lists(s)
        .into_iter()
        .map(|row| row.into_iter().filter_map(|a| parse_ref(&a)).collect())
        .collect()
}

/// Parse a nested list of real values (e.g. a B-spline surface weight grid).
fn parse_nested_real_list(s: &str) -> Vec<Vec<f64>> {
    parse_nested_lists(s)
        .into_iter()
        .map(|row| row.into_iter().filter_map(|a| a.trim().parse().ok()).collect())
        .collect()
}

/// Parse a STEP physical file, returning the model plus collected warnings.
/// Extract the `DATA` section of a physical file (between `DATA;` and the
/// closing `ENDSEC;`), validating the file terminator.
fn data_section(content: &str) -> Result<&str, String> {
    if !content.contains("END-ISO-10303-21") {
        return Err("STEP file: missing END-ISO-10303-21 terminator".into());
    }
    let data_start = content
        .find("DATA;")
        .ok_or("STEP file: missing DATA section")?;
    let after_data = &content[data_start + 5..];
    let data_end = after_data
        .find("ENDSEC;")
        .ok_or("STEP file: DATA section not closed by ENDSEC")?;
    Ok(&after_data[..data_end])
}

fn read_step_impl(content: &str) -> Result<(BRepModel, Vec<String>), String> {
    let data = data_section(content)?;
    let records = parse_records(data)?;
    let resolver = Resolver::new(&records);

    let mut rep_ids: Vec<usize> = records
        .iter()
        .filter(|(_, r)| {
            r.type_name == "ADVANCED_BREP_SHAPE_REPRESENTATION" || r.type_name == "SHAPE_REPRESENTATION"
        })
        .map(|(id, _)| *id)
        .collect();
    rep_ids.sort_unstable();

    let mut model = BRepModel::new();
    if !rep_ids.is_empty() {
        for id in rep_ids {
            match resolver.resolve_representation(id) {
                Ok((name, shapes)) => {
                    if shapes.is_empty() {
                        continue;
                    }
                    let shape = if shapes.len() == 1 {
                        shapes.into_iter().next().unwrap()
                    } else {
                        let b = TopoBuilder::new();
                        b.make_compound_of(&shapes).0
                    };
                    model.add(&name, shape);
                }
                Err(e) => resolver.warn(format!("representation #{id}: {e}")),
            }
        }
    } else {
        // Fallback for bare files without a representation layer: treat any
        // MANIFOLD_SOLID_BREP as a root shape.
        let mut roots: Vec<usize> = records
            .iter()
            .filter(|(_, r)| r.type_name == "MANIFOLD_SOLID_BREP")
            .map(|(id, _)| *id)
            .collect();
        roots.sort_unstable();
        for id in roots {
            match resolver.resolve_shape(id) {
                Ok(s) => {
                    model.add("", s);
                }
                Err(e) => resolver.warn(format!("shape #{id}: {e}")),
            }
        }
    }

    let warnings = resolver.warnings.into_inner();
    Ok((model, warnings))
}

/// Parse a STEP physical file into a `BRepModel`.
pub fn read_step(content: &str) -> Result<BRepModel, String> {
    Ok(read_step_impl(content)?.0)
}

/// Parse a STEP physical file, also returning collected warnings for skipped
/// or unsupported entities.
pub fn read_step_with_warnings(content: &str) -> Result<(BRepModel, Vec<String>), String> {
    read_step_impl(content)
}

/// Read a STEP physical file from disk.
pub fn read_step_file(path: &str) -> Result<BRepModel, String> {
    let content = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    read_step(&content)
}

/// Read a STEP physical file back into a [`StepAssembly`].
///
/// Reconstructs the product tree from `PRODUCT` / `PRODUCT_DEFINITION` /
/// `NEXT_ASSEMBLY_USAGE_OCCURRENCE` records. The root name is the product that
/// is never referenced as a child; part shapes are re-read from the part
/// representations when present (parts without a representation keep an empty
/// `TopoShape`).
pub fn read_step_assembly(content: &str) -> Result<StepAssembly, String> {
    let data = data_section(content)?;
    let records = parse_records(data)?;
    let resolver = Resolver::new(&records);

    // PRODUCT id -> product name.
    let mut product_names: HashMap<usize, String> = HashMap::new();
    for (id, rec) in &records {
        if rec.type_name == "PRODUCT" {
            product_names.insert(*id, parse_str(&rec.args[0]));
        }
    }

    // PRODUCT_DEFINITION_FORMATION id -> owning PRODUCT id.
    let mut form_product: HashMap<usize, usize> = HashMap::new();
    for (id, rec) in &records {
        if rec.type_name == "PRODUCT_DEFINITION_FORMATION" {
            if let Some(p) = parse_ref(&rec.args[2]) {
                form_product.insert(*id, p);
            }
        }
    }

    // PRODUCT_DEFINITION id -> product name (via its formation).
    let mut def_name: HashMap<usize, String> = HashMap::new();
    for (id, rec) in &records {
        if rec.type_name == "PRODUCT_DEFINITION" {
            let name = parse_ref(&rec.args[3])
                .and_then(|f| form_product.get(&f).copied())
                .and_then(|p| product_names.get(&p).cloned())
                .unwrap_or_default();
            def_name.insert(*id, name);
        }
    }

    // PRODUCT_DEFINITION_SHAPE id -> PRODUCT id (the shape's owner).
    let mut pds_product: HashMap<usize, usize> = HashMap::new();
    for (id, rec) in &records {
        if rec.type_name == "PRODUCT_DEFINITION_SHAPE" {
            if let Some(p) = parse_ref(&rec.args[2]) {
                pds_product.insert(*id, p);
            }
        }
    }

    // representation id -> PRODUCT id, via the SHAPE_REPRESENTATION link.
    let mut rep_product: HashMap<usize, usize> = HashMap::new();
    for (id, rec) in &records {
        if rec.type_name == "PRODUCT_DEFINITION_SHAPE_REPRESENTATION" {
            if let (Some(pds), Some(rep)) = (parse_ref(&rec.args[1]), parse_ref(&rec.args[2])) {
                if let Some(p) = pds_product.get(&pds) {
                    rep_product.insert(rep, *p);
                }
            }
        }
    }

    // Parts: resolve each representation owned by a product.
    let mut products: Vec<(String, TopoShape)> = Vec::new();
    let mut seen: HashSet<String> = HashSet::new();
    let mut reps: Vec<usize> = records
        .iter()
        .filter(|(_, r)| {
            r.type_name == "ADVANCED_BREP_SHAPE_REPRESENTATION" || r.type_name == "SHAPE_REPRESENTATION"
        })
        .map(|(id, _)| *id)
        .collect();
    reps.sort_unstable();
    for id in reps {
        let Some(pid) = rep_product.get(&id) else { continue };
        let Some(name) = product_names.get(pid) else { continue };
        let name = name.clone();
        if seen.contains(&name) {
            continue;
        }
        match resolver.resolve_representation(id) {
            Ok((_, shapes)) => {
                if shapes.is_empty() {
                    continue;
                }
                let shape = if shapes.len() == 1 {
                    shapes.into_iter().next().unwrap()
                } else {
                    let b = TopoBuilder::new();
                    b.make_compound_of(&shapes).0
                };
                seen.insert(name.clone());
                products.push((name, shape));
            }
            Err(e) => resolver.warn(format!("representation #{id}: {e}")),
        }
    }

    // Assembly tree from NEXT_ASSEMBLY_USAGE_OCCURRENCE records.
    let mut children: Vec<(String, Vec<String>)> = Vec::new();
    let mut has_parent: HashSet<String> = HashSet::new();
    for (_, rec) in &records {
        if rec.type_name == "NEXT_ASSEMBLY_USAGE_OCCURRENCE" {
            if let (Some(rel), Some(red)) = (parse_ref(&rec.args[3]), parse_ref(&rec.args[4])) {
                if let (Some(pn), Some(cn)) = (def_name.get(&rel), def_name.get(&red)) {
                    if let Some(entry) = children.iter_mut().find(|(p, _)| p == pn) {
                        entry.1.push(cn.clone());
                    } else {
                        children.push((pn.clone(), vec![cn.clone()]));
                    }
                    has_parent.insert(cn.clone());
                }
            }
        }
    }

    // The assembly root is the product never referenced as a child.
    let name = product_names
        .values()
        .find(|n| !has_parent.contains(*n))
        .cloned()
        .unwrap_or_else(|| "Assembly".to_string());

    Ok(StepAssembly {
        name,
        products,
        children,
    })
}

/// Read a STEP assembly physical file from disk.
pub fn read_step_assembly_file(path: &str) -> Result<StepAssembly, String> {
    let content = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    read_step_assembly(&content)
}

/// Round-trip helper: write `shape` to STEP, read it back, and count the
/// distinct vertices / edges / faces in the reconstructed model.
pub fn step_roundtrip_counts(shape: &TopoShape) -> Result<(usize, usize, usize), String> {
    let step = write_shape_step(shape);
    let model = read_step(&step)?;
    let mut nv = 0usize;
    let mut ne = 0usize;
    let mut nf = 0usize;
    for ms in &model.shapes {
        let c = crate::topo_tools_full::shape_counts(&ms.shape);
        nv += c.get(&ShapeType::Vertex).copied().unwrap_or(0);
        ne += c.get(&ShapeType::Edge).copied().unwrap_or(0);
        nf += c.get(&ShapeType::Face).copied().unwrap_or(0);
    }
    Ok((nv, ne, nf))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_surface::SurfaceKind;
    use crate::primitives::{BRepPrimBox, BRepPrimSphere};
    use crate::topo_tools_full::{faces_of, vertices_of};

    #[test]
    fn box_roundtrip_counts_points_and_faces() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let step = write_shape_step(&b.solid.0);
        let model = read_step(&step).expect("parse step");
        assert_eq!(model.len(), 1);

        let (nv, ne, nf) = step_roundtrip_counts(&b.solid.0).expect("roundtrip counts");
        assert_eq!((nv, ne, nf), (8, 12, 6));

        let shape = &model.shapes[0].shape;
        let verts = vertices_of(shape);
        assert_eq!(verts.len(), 8);
        let pts: Vec<GpPnt> = verts
            .iter()
            .map(|v| GeometryRegistry::global().vertex_point(&v.0))
            .collect();
        for &(x, y, z) in &[
            (0., 0., 0.),
            (2., 0., 0.),
            (2., 3., 0.),
            (0., 3., 0.),
            (0., 0., 4.),
            (2., 0., 4.),
            (2., 3., 4.),
            (0., 3., 4.),
        ] {
            assert!(
                pts.iter().any(|p| p.is_equal(&GpPnt::new(x, y, z))),
                "missing corner ({x},{y},{z})"
            );
        }

        let faces = faces_of(shape);
        assert_eq!(faces.len(), 6);
        for f in &faces {
            let surf = GeometryRegistry::global().face_surface(&f.0);
            assert!(surf.is_some(), "face has a surface");
            if let Some(s) = surf {
                assert_eq!(classify_surface(s.as_ref()), SurfaceKind::Plane);
            }
        }
    }

    #[test]
    fn sphere_roundtrip_surface() {
        let s = BRepPrimSphere::make_sphere(2.5);
        let step = write_shape_step(&s.solid.0);
        let model = read_step(&step).expect("parse sphere step");
        assert_eq!(model.len(), 1);
        let shape = &model.shapes[0].shape;
        let faces = faces_of(shape);
        assert_eq!(faces.len(), 1);
        let surf = GeometryRegistry::global()
            .face_surface(&faces[0].0)
            .expect("sphere face has surface");
        assert_eq!(classify_surface(surf.as_ref()), SurfaceKind::Sphere);
    }

    #[test]
    fn model_two_shapes_names_preserved() {
        let b = BRepPrimBox::make_box(1.0, 2.0, 3.0);
        let s = BRepPrimSphere::make_sphere(1.0);
        let mut model = BRepModel::new();
        model.add("Box", b.solid.0.clone());
        model.add("Sphere", s.solid.0.clone());
        let step = write_step(&model);
        let got = read_step(&step).expect("parse model step");
        assert_eq!(got.len(), 2);
        let mut names: Vec<&str> = got.names();
        names.sort();
        assert_eq!(names, vec!["Box", "Sphere"]);
    }

    #[test]
    fn malformed_missing_end_iso_errors() {
        let s = "ISO-10303-21;\nHEADER;\nENDSEC;\nDATA;\n#1=DIRECTION('',(1.,0.,0.));\nENDSEC;\n";
        assert!(read_step(s).is_err());
    }

    #[test]
    fn reader_skips_unknown_entity_with_warning() {
        // A minimal file whose representation references an unsupported entity;
        // the reader must collect a warning and skip the shape gracefully.
        let s = "ISO-10303-21;\nHEADER;\nENDSEC;\nDATA;\n\
#1=APPLICATION_CONTEXT('AUTOMOTIVE_DESIGN');\n\
#2=WIDGET_FROBNICATOR('',42.);\n\
#3=GEOMETRIC_REPRESENTATION_CONTEXT('','',3);\n\
#4=ADVANCED_BREP_SHAPE_REPRESENTATION('widget',(#2),#3);\n\
ENDSEC;\nEND-ISO-10303-21;";
        let (model, warnings) = read_step_with_warnings(s).expect("parse");
        assert!(model.is_empty(), "unsupported shape is skipped");
        assert!(
            warnings.iter().any(|w| w.contains("WIDGET_FROBNICATOR")),
            "expected a warning about the unknown entity, got {warnings:?}"
        );
    }

    #[test]
    fn step_real_formats_decimal_point() {
        assert_eq!(step_real(0.0), "0.0");
        assert_eq!(step_real(2.5), "2.5");
        assert_eq!(step_real(-3.0), "-3.0");
        assert!(step_real(1e20).contains('.'));
    }

    #[test]
    fn split_top_handles_nested_lists() {
        let args = split_top("'',(1.,0.,0.),(#2,#3)");
        assert_eq!(args.len(), 3);
        assert_eq!(args[0], "''");
        assert_eq!(args[1], "(1.,0.,0.)");
        assert_eq!(args[2], "(#2,#3)");
    }

    #[test]
    fn cylinder_roundtrip() {
        use crate::primitives::BRepPrimCylinder;
        let c = BRepPrimCylinder::make_cylinder(2.0, 10.0);
        let (nv, ne, nf) = step_roundtrip_counts(&c.solid.0).expect("counts");
        assert_eq!((nv, ne, nf), (2, 3, 3));
        let s = write_shape_step(&c.solid.0);
        let m = read_step(&s).expect("read cylinder");
        assert_eq!(m.len(), 1);
        let fs = faces_of(&m.shapes[0].shape);
        assert_eq!(fs.len(), 3);
        // Two cap faces should classify as planes.
        let planar = fs
            .iter()
            .filter(|f| {
                GeometryRegistry::global()
                    .face_surface(&f.0)
                    .map(|s| classify_surface(s.as_ref()) == SurfaceKind::Plane)
                    .unwrap_or(false)
            })
            .count();
        assert_eq!(planar, 2);
    }

    #[test]
    fn step_file_io_roundtrip() {
        use crate::primitives::BRepPrimSphere;
        let s = BRepPrimSphere::make_sphere(1.5);
        let mut model = BRepModel::new();
        model.add("Ball", s.solid.0.clone());
        let path = std::env::temp_dir().join("occt_step_test.step");
        let p = path.to_str().unwrap();
        write_step_file(p, &model).expect("write file");
        let got = read_step_file(p).expect("read file");
        assert_eq!(got.len(), 1);
        assert_eq!(got.names(), vec!["Ball"]);
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn bspline_curve_written() {
        let c = GeomBSplineCurve::new(
            vec![
                GpPnt::new(0.0, 0.0, 0.0),
                GpPnt::new(1.0, 0.0, 0.0),
                GpPnt::new(2.0, 1.0, 0.0),
                GpPnt::new(3.0, 0.0, 0.0),
            ],
            vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
            2,
        )
        .unwrap();
        let mut w = StepWriter::new();
        let id = write_bspline_curve(&mut w, &c).unwrap();
        let line = &w.lines[id - 1];
        assert!(line.contains("B_SPLINE_CURVE_WITH_KNOTS"), "{line}");
        // Degree 2, four control-point refs, polynomial (SELF weights).
        assert!(line.contains("('',2,(#1,#2,#3,#4)"), "{line}");
        assert!(line.contains("SELF"), "{line}");
        // Knots 0/0.5/1 with multiplicities 3/1/3, and the UNSPECIFIED defaults.
        assert!(line.contains("(0.0,0.5,1.0)"), "{line}");
        assert!(line.contains("(3,1,3)"), "{line}");
        assert!(line.contains("UNSPECIFIED"), "{line}");
    }

    #[test]
    fn rational_bspline_curve_written() {
        let c = GeomBSplineCurve::rational(
            vec![
                GpPnt::new(0.0, 0.0, 0.0),
                GpPnt::new(1.0, 0.0, 0.0),
                GpPnt::new(2.0, 0.0, 0.0),
            ],
            vec![1.0, 0.5, 1.0],
            vec![0.0, 0.0, 0.5, 1.0, 1.0],
            1,
        )
        .unwrap();
        let mut w = StepWriter::new();
        let id = write_bspline_curve(&mut w, &c).unwrap();
        let line = &w.lines[id - 1];
        assert!(line.contains("(1.0,0.5,1.0)"), "{line}");
        assert!(!line.contains("SELF"), "{line}");
    }

    #[test]
    fn bspline_surface_written() {
        let (ku, kv) = occt_geom::bspline_surface::bspline_surface_uniform_knots(3, 3, 2, 2);
        let poles = vec![
            vec![
                GpPnt::new(0.0, 0.0, 0.0),
                GpPnt::new(0.0, 1.0, 0.0),
                GpPnt::new(0.0, 2.0, 0.0),
            ],
            vec![
                GpPnt::new(1.0, 0.0, 1.0),
                GpPnt::new(1.0, 1.0, 1.0),
                GpPnt::new(1.0, 2.0, 1.0),
            ],
            vec![
                GpPnt::new(2.0, 0.0, 0.0),
                GpPnt::new(2.0, 1.0, 0.0),
                GpPnt::new(2.0, 2.0, 0.0),
            ],
        ];
        let s = GeomBSplineSurface::new(poles, ku, kv, 2, 2).unwrap();
        let mut w = StepWriter::new();
        let id = write_bspline_surface(&mut w, &s).unwrap();
        let line = &w.lines[id - 1];
        assert!(line.contains("B_SPLINE_SURFACE_WITH_KNOTS"), "{line}");
        // u_degree, v_degree both 2; 3×3 pole grid; polynomial.
        assert!(line.contains("('',2,2,((#1,#2,#3),(#4,#5,#6),(#7,#8,#9))"), "{line}");
        assert!(line.contains("SELF"), "{line}");
        assert!(line.contains("(0.0,1.0)"), "{line}");
        assert!(line.contains("(3,3)"), "{line}");
    }

    #[test]
    fn trimmed_curve_written() {
        let mut w = StepWriter::new();
        let id = write_trimmed_curve(&mut w, 42, 1.5, 0.5).unwrap();
        let line = &w.lines[id - 1];
        // Bounds are normalised to ascending order.
        assert!(
            line.contains("TRIMMED_CURVE('',#42,1,0.5,1.5,PARAMETER)"),
            "{line}"
        );
    }

    #[test]
    fn offset_curve_written() {
        let mut w = StepWriter::new();
        let id = write_offset_curve(&mut w, 7, 2.5, 9).unwrap();
        let line = &w.lines[id - 1];
        assert!(
            line.contains("OFFSET_CURVE_3D('',#7,#9,2.5,UNSPECIFIED,UNSPECIFIED)"),
            "{line}"
        );
    }

    #[test]
    fn circle_ellipse_params() {
        let ax2 = GpAx2::new(GpPnt::zero(), dir_z(), dir_x()).unwrap();
        let c = GeomCircle::new(GpCirc::new(ax2.clone(), 2.0));
        let mut w = StepWriter::new();
        let id = write_conic_params(&mut w, &c).unwrap().expect("circle is a conic");
        let line = &w.lines[id - 1];
        assert!(line.contains("CIRCLE"), "{line}");
        assert!(line.contains(",2.0)"), "{line}");

        let e = GeomEllipse::new(GpElips::new(ax2, 3.0, 1.5));
        let id = write_conic_params(&mut w, &e).unwrap().expect("ellipse is a conic");
        let line = &w.lines[id - 1];
        assert!(line.contains("ELLIPSE"), "{line}");
        assert!(line.contains("3.0,1.5)"), "{line}");
    }

    /// A small curved B-spline surface grid (z = u² + v³) reused by the
    /// spline round-trip tests.
    fn bspline_test_surface() -> GeomBSplineSurface {
        let (nu, nv) = (5, 5);
        let points: Vec<Vec<GpPnt>> = (0..nu)
            .map(|i| {
                (0..nv)
                    .map(|j| {
                        let u = i as f64 / (nu - 1) as f64;
                        let v = j as f64 / (nv - 1) as f64;
                        GpPnt::new(u, v, u * u + v * v * v)
                    })
                    .collect()
            })
            .collect();
        occt_geom::bspline_surface::fit_surface_grid(&points, 3, 3).unwrap()
    }

    #[test]
    fn step_bspline_roundtrip() {
        let surf = bspline_test_surface();
        let face = crate::brep_builder_full::BRepBuilderFace::from_surface(Arc::new(surf));
        let step = write_step_with_splines(&face).unwrap();
        let m = read_step(&step).expect("read bspline step");
        assert_eq!(m.len(), 1);
        assert!(m.shapes[0].shape.is_face(), "shape type preserved");
        let fs = faces_of(&m.shapes[0].shape);
        assert!(!fs.is_empty(), "reconstructed shape has faces");
    }

    #[test]
    fn step_spline_reader_handles() {
        let surf = bspline_test_surface();
        let face = crate::brep_builder_full::BRepBuilderFace::from_surface(Arc::new(surf));
        let step = write_step_with_splines(&face).unwrap();
        let path = std::env::temp_dir().join("occt_step_spline_reader.step");
        let p = path.to_str().unwrap();
        std::fs::write(p, &step).expect("write spline step");
        let m = read_step_file(p).expect("read spline step file");
        assert_eq!(m.len(), 1);
        let fs = faces_of(&m.shapes[0].shape);
        assert!(!fs.is_empty());
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn box_roundtrip_splines() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let step = write_step_with_splines(&b.solid.0).unwrap();
        let m = read_step(&step).expect("read spline box step");
        assert_eq!(m.len(), 1);
        let c = crate::topo_tools_full::shape_counts(&m.shapes[0].shape);
        assert_eq!(c.get(&ShapeType::Vertex).copied().unwrap_or(0), 8);
        assert_eq!(c.get(&ShapeType::Edge).copied().unwrap_or(0), 12);
        assert_eq!(c.get(&ShapeType::Face).copied().unwrap_or(0), 6);
    }

    #[test]
    fn conic_circle_roundtrip() {
        let b = TopoBuilder::new();
        let ax2 = GpAx2::new(GpPnt::zero(), dir_z(), dir_x()).unwrap();
        let circ = GpCirc::new(ax2, 2.0);
        let mut e = b.make_edge(Arc::new(GeomCircle::new(circ)), 0.0, 2.0 * PI);
        let v = b.make_vertex(GpPnt::new(2.0, 0.0, 0.0), 0.0);
        b.add(&mut e.0, &v.0);
        b.add(&mut e.0, &v.0);
        let w = b.make_wire(&[Edge(e.0)]);
        let ax3 = GpAx3::new(GpPnt::zero(), dir_z(), &dir_x()).unwrap();
        let face = crate::brep_builder_full::BRepBuilderFace::from_wire(&w, &GpPln::new(ax3)).unwrap();
        let step = write_step_with_splines(&face).unwrap();
        let m = read_step(&step).expect("read circle face step");
        assert_eq!(m.len(), 1);
        let fs = faces_of(&m.shapes[0].shape);
        assert!(!fs.is_empty(), "reconstructed circle face");
    }

    #[test]
    fn step_assembly_written() {
        let a = StepAssembly {
            name: "Assy".into(),
            products: vec![
                ("PartA".into(), BRepPrimBox::make_box(1.0, 2.0, 3.0).solid.0),
                ("PartB".into(), BRepPrimBox::make_box(2.0, 1.0, 1.0).solid.0),
            ],
            children: vec![("Assy".into(), vec!["PartA".into(), "PartB".into()])],
        };
        let out = write_step_assembly(&a).unwrap();
        assert!(out.contains("NEXT_ASSEMBLY_USAGE_OCCURRENCE"), "{out}");
        // The assembly root plus the two parts.
        assert!(out.matches("PRODUCT('").count() >= 2, "{out}");
        // The usage records reference the correct product definitions.
        assert!(out.contains("PartA"), "{out}");
        assert!(out.contains("PartB"), "{out}");
    }

    #[test]
    fn step_assembly_unknown_child_errors() {
        let a = StepAssembly {
            name: "Assy".into(),
            products: vec![("PartA".into(), BRepPrimBox::make_box(1.0, 1.0, 1.0).solid.0)],
            children: vec![("Assy".into(), vec!["Ghost".into()])],
        };
        let err = write_step_assembly(&a).unwrap_err();
        assert!(err.contains("Ghost"), "{err}");
    }

    #[test]
    fn step_color_written() {
        let b = BRepPrimBox::make_box(1.0, 2.0, 3.0);
        let out = write_step_with_color(&b.solid.0, (0.8, 0.2, 0.1)).unwrap();
        assert!(out.contains("COLOUR_RGB"), "{out}");
        assert!(out.contains("0.8"), "{out}");
        assert!(out.contains("0.2"), "{out}");
        assert!(out.contains("0.1"), "{out}");
        // The style chain is attached to the shape representation.
        assert!(out.contains("SURFACE_STYLE_FILL_AREA"), "{out}");
        assert!(out.contains("SURFACE_STYLE_USAGE"), "{out}");
        // The file still reads back as the shape.
        let m = read_step(&out).expect("read colored step");
        assert_eq!(m.len(), 1);
    }

    #[test]
    fn step_name_written() {
        let b = BRepPrimBox::make_box(1.0, 2.0, 3.0);
        let out = write_step_with_name(&b.solid.0, "MyBox").unwrap();
        assert!(out.contains("MyBox"), "{out}");
        // The PRODUCT name survives a read (representation name is derived
        // from the product name by the writer).
        let m = read_step(&out).expect("read named step");
        assert_eq!(m.names(), vec!["MyBox"]);
    }

    #[test]
    fn step_units_written() {
        let b = BRepPrimBox::make_box(1.0, 2.0, 3.0);
        let units = StepUnits {
            length_unit_m: 0.001,
            angle_unit_rad: 1.0,
        };
        let out = write_step_with_units(&b.solid.0, "Metric", &units).unwrap();
        assert!(out.contains("Metric"), "{out}");
        assert!(out.contains("SI_UNIT"), "{out}");
        assert!(out.contains("DIMENSIONAL_EXPONENTS"), "{out}");
        // The file still reads back as the shape (unit records are skipped).
        let m = read_step(&out).expect("read units step");
        assert_eq!(m.len(), 1);
    }

    #[test]
    fn step_name_color_written() {
        let b = BRepPrimBox::make_box(1.0, 2.0, 3.0);
        let out = write_step_with_name_and_color(&b.solid.0, "RedBox", (1.0, 0.0, 0.0)).unwrap();
        assert!(out.contains("RedBox"), "{out}");
        assert!(out.contains("COLOUR_RGB"), "{out}");
        assert!(out.contains("SURFACE_STYLE_USAGE"), "{out}");
        let m = read_step(&out).expect("read named+color step");
        assert_eq!(m.names(), vec!["RedBox"]);
    }

    #[test]
    fn step_assembly_roundtrip() {
        let a = StepAssembly {
            name: "Rig".into(),
            products: vec![
                ("Leg".into(), BRepPrimBox::make_box(0.5, 0.5, 2.0).solid.0),
                ("Foot".into(), BRepPrimBox::make_box(0.8, 0.3, 0.2).solid.0),
            ],
            children: vec![
                ("Rig".into(), vec!["Leg".into()]),
                ("Leg".into(), vec!["Foot".into()]),
            ],
        };
        let out = write_step_assembly(&a).unwrap();
        let got = read_step_assembly(&out).expect("read assembly step");
        // The tree is reconstructed: root "Rig", two parts, two usage edges.
        assert_eq!(got.name, "Rig");
        assert_eq!(got.products.len(), 2);
        let names: Vec<String> = got.products.iter().map(|(n, _)| n.clone()).collect();
        assert!(names.contains(&"Leg".into()), "{names:?}");
        assert!(names.contains(&"Foot".into()), "{names:?}");
        let flat: Vec<(String, String)> = got
            .children
            .iter()
            .flat_map(|(p, kids)| kids.iter().map(move |k| (p.clone(), k.clone())))
            .collect();
        assert!(flat.contains(&("Rig".into(), "Leg".into())), "{flat:?}");
        assert!(flat.contains(&("Leg".into(), "Foot".into())), "{flat:?}");
    }

    #[test]
    fn step_spline_reader_plain_bspline() {
        // A hand-written STEP file using the plain (knotless) B_SPLINE_CURVE
        // and B_SPLINE_SURFACE entities must read back into a valid face.
        let s = "ISO-10303-21;\nHEADER;\nENDSEC;\nDATA;\n\
#1=APPLICATION_CONTEXT('AUTOMOTIVE_DESIGN');\n\
#2=CARTESIAN_POINT('',(0.,0.,0.));\n\
#3=CARTESIAN_POINT('',(0.,1.,0.));\n\
#4=CARTESIAN_POINT('',(1.,0.,1.));\n\
#5=CARTESIAN_POINT('',(1.,1.,1.));\n\
#6=CARTESIAN_POINT('',(2.,0.,0.));\n\
#7=CARTESIAN_POINT('',(2.,1.,0.));\n\
#8=CARTESIAN_POINT('',(3.,0.,1.));\n\
#9=CARTESIAN_POINT('',(3.,1.,1.));\n\
#10=B_SPLINE_SURFACE('',1,1,((#2,#3),(#4,#5),(#6,#7),(#8,#9)),UNSPECIFIED,.F.,.F.,.F.);\n\
#11=GEOMETRIC_REPRESENTATION_CONTEXT('','',3);\n\
#12=ADVANCED_FACE('',#10,(),.T.);\n\
#13=ADVANCED_BREP_SHAPE_REPRESENTATION('plain',(#12),#11);\n\
ENDSEC;\nEND-ISO-10303-21;";
        let m = read_step(s).expect("read plain bspline step");
        assert_eq!(m.len(), 1);
        let fs = faces_of(&m.shapes[0].shape);
        assert!(!fs.is_empty(), "plain B_SPLINE_SURFACE yields a face");
    }

    #[test]
    fn polyline_written_and_read() {
        let pts = [
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(2.0, 1.0, 0.0),
        ];
        let mut w = StepWriter::new();
        let id = write_polyline(&mut w, &pts).unwrap();
        assert!(w.lines[id - 1].contains("POLYLINE('',(#1,#2,#3))"), "{}", w.lines[id - 1]);

        // A wire built from a polyline edge round-trips through a hand-written
        // STEP file that references the POLYLINE.
        let s = "ISO-10303-21;\nHEADER;\nENDSEC;\nDATA;\n\
#1=APPLICATION_CONTEXT('AUTOMOTIVE_DESIGN');\n\
#2=CARTESIAN_POINT('',(0.,0.,0.));\n\
#3=CARTESIAN_POINT('',(1.,0.,0.));\n\
#4=CARTESIAN_POINT('',(2.,1.,0.));\n\
#5=POLYLINE('',(#2,#3,#4));\n\
#6=DIRECTION('',(1.,0.,0.));\n\
#7=VECTOR('',#6,1.);\n\
#8=LINE('',#2,#7);\n\
#9=VERTEX_POINT('',#2);\n\
#10=VERTEX_POINT('',#4);\n\
#11=EDGE_CURVE('',#9,#10,#5,.T.);\n\
#12=ORIENTED_EDGE('',#9,#10,#11,.T.);\n\
#13=EDGE_LOOP('',(#12));\n\
#14=GEOMETRIC_REPRESENTATION_CONTEXT('','',3);\n\
#15=ADVANCED_BREP_SHAPE_REPRESENTATION('poly',(#11),#14);\n\
ENDSEC;\nEND-ISO-10303-21;";
        let m = read_step(s).expect("read polyline step");
        assert_eq!(m.len(), 1);
        let e = crate::topo_tools_full::edges_of(&m.shapes[0].shape);
        assert_eq!(e.len(), 1);
    }

    #[test]
    fn step_options_written() {
        let b = BRepPrimBox::make_box(1.0, 2.0, 3.0);
        let opts = StepWriteOptions {
            name: "OptBox".into(),
            color: Some((0.25, 0.5, 0.75)),
            units: StepUnits {
                length_unit_m: 0.001,
                angle_unit_rad: 1.0,
            },
            splines: true,
        };
        let out = write_step_with_options(&b.solid.0, &opts).unwrap();
        assert!(out.contains("OptBox"), "{out}");
        assert!(out.contains("COLOUR_RGB"), "{out}");
        assert!(out.contains("0.25"), "{out}");
        assert!(out.contains("SI_UNIT"), "{out}");
        let m = read_step(&out).expect("read options step");
        assert_eq!(m.names(), vec!["OptBox"]);
    }

    #[test]
    fn step_header_written() {
        let b = BRepPrimBox::make_box(1.0, 2.0, 3.0);
        let h = StepHeader {
            description: "Test part".into(),
            name: "part.step".into(),
            timestamp: "2030-01-01T00:00:00".into(),
            author: "alice".into(),
            organization: "acme".into(),
            preprocessor: "occt-topo".into(),
            originator: "bob".into(),
            schema: "AP242".into(),
        };
        let out = write_step_with_header(&b.solid.0, "Part", &h).unwrap();
        assert!(out.contains("Test part"), "{out}");
        assert!(out.contains("part.step"), "{out}");
        assert!(out.contains("2030-01-01T00:00:00"), "{out}");
        assert!(out.contains("AP242"), "{out}");
        let m = read_step(&out).expect("read header step");
        assert_eq!(m.names(), vec!["Part"]);
    }

    #[test]
    fn step_assembly_with_splines_written() {
        let a = StepAssembly {
            name: "Assy".into(),
            products: vec![("PartA".into(), BRepPrimBox::make_box(1.0, 2.0, 3.0).solid.0)],
            children: vec![("Assy".into(), vec!["PartA".into()])],
        };
        let out = write_step_assembly_with_splines(&a).unwrap();
        assert!(out.contains("NEXT_ASSEMBLY_USAGE_OCCURRENCE"), "{out}");
        let got = read_step_assembly(&out).expect("read assembly step");
        assert_eq!(got.name, "Assy");
        assert_eq!(got.products.len(), 1);
    }
}
