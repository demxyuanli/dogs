//! IGES (ANSI Y14.26M) B-Rep writer — a port of `IGESControl_Writer`.
//!
//! Emits an 80-column fixed-width IGES file with the classic entity set for
//! B-rep solids: 116 POINT, 110 LINE, 100 CIRCULAR ARC, 108 PLANE,
//! 120 SURFACE OF REVOLUTION, 144 TRIMMED SURFACE, 510 FACE, 514 SHELL and
//! 186 MANIFOLD SOLID BREP OBJECT. Curves and surfaces cannot be downcast from
//! `Arc<dyn Curve>` / `Arc<dyn Surface>`, so geometry is classified by sampling
//! invariants (constant zero second derivative ⇒ line, constant curvature ⇒
//! circular arc, equidistant samples ⇒ sphere). Faces without boundary wires
//! (e.g. an untrimmed sphere) get a synthesized revolution-surface and
//! meridian-arc boundary.
//!
//! The physical layout is the standard sectioned format (S/G/D/P/T) with each
//! line exactly 80 characters. Section lines carry their section letter and a
//! 7-digit sequence number in the first 8 columns, followed by the data.

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::gp::{GpAx3, GpDir, GpPln, GpPnt, GpVec};
use occt_geom::{Curve, Surface};

use crate::abs::ShapeType;
use crate::brep_surface::{classify_surface, face_is_planar, face_plane, sphere_center, SurfaceKind};
use crate::brep_tool::BRepTool;
use crate::model::BRepModel;
use crate::shape::{Edge, Face, Shell, Solid, TopoShape};
use crate::topo_tools_full::{edge_vertices, edges_of_wire, vertices_of, wires_of_face};

// ---------------------------------------------------------------------------
// Formatting helpers
// ---------------------------------------------------------------------------

/// Pad or truncate `s` to exactly 80 columns.
pub fn rec80(s: &str) -> String {
    let mut line = String::with_capacity(80);
    line.push_str(s);
    if line.len() > 80 {
        line.truncate(80);
    } else {
        line.push_str(&" ".repeat(80 - line.len()));
    }
    line
}

/// Data columns of a Start/Global/Directory card
/// (`IGESData_IGESWriter.cxx:38-39`: `MaxcarsG = 72`).
const MAXCARS_G: usize = 72;
/// Data columns of a Parameter card (`cxx:39`: `MaxcarsP = 64`).
const MAXCARS_P: usize = 64;

/// One IGES card: the data field left-justified in `width` columns, then
/// `section` and the card's sequence number in the last seven columns
/// (`IGESData_IGESWriter.cxx:836-880` for D, `:769-794` for G/S).
fn sec_line(section: char, seq: usize, content: &str, width: usize) -> String {
    let mut s = String::with_capacity(80);
    for ch in content.chars().take(width) {
        s.push(ch);
    }
    s.push_str(&" ".repeat(width - s.len()));
    s.push(section);
    s.push_str(&format!("{seq:>7}"));
    s
}

/// One Parameter card (`IGESData_IGESWriter.cxx:902-925`):
/// `data` (64 columns) + blank + the owning entity's directory-entry pointer
/// (`2*i - 1`) + `P` + the card's sequence number.
fn param_line(content: &str, de_pointer: usize, seq: usize) -> String {
    let mut s = String::with_capacity(80);
    for ch in content.chars().take(MAXCARS_P) {
        s.push(ch);
    }
    s.push_str(&" ".repeat(MAXCARS_P - s.len()));
    s.push_str(&format!(" {de_pointer:>7}P{seq:>7}"));
    s
}

/// Right-justify a value in an 8-column IGES field.
fn field8(v: &str) -> String {
    format!("{v:>8}")
}

/// Render an f64 in IGES parameter syntax (always carries a decimal point or
/// exponent marker).
fn num(x: f64) -> String {
    if !x.is_finite() {
        return "0.0".into();
    }
    let s = format!("{:?}", x);
    if s.contains('.') || s.contains('e') || s.contains('E') {
        s
    } else {
        format!("{s}.")
    }
}

/// Finite sampling bounds of a surface (unbounded ranges clamp to ±1).
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

/// 3×3 determinant.
fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Circumcenter of three non-collinear 3D points (Cramer's rule on the
/// perpendicular-bisector system).
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
// Curve classification
// ---------------------------------------------------------------------------

/// Curve kind for the IGES curve dispatch.
///
/// `GeomToIGES_GeomCurve::TransferCurve` (`GeomToIGES_GeomCurve.cxx:75-116`)
/// dispatches on the curve's **exact type** (`IsKind`), never on sampled
/// geometry: `Geom_BoundedCurve` (Bezier / B-spline / trimmed) → the 126/100
/// family, `Geom_Conic` → 104 (circle 100), `Geom_OffsetCurve`, `Geom_Line` →
/// 110. The port has the equivalent type queries on [`Curve`], so the previous
/// "6 samples of |C''| with a 2% spread" classifier (audit A26) is replaced by
/// this dispatch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum IgCurveKind {
    Line,
    Circle,
    /// Ellipse / hyperbola / parabola: OCCT emits conic-arc entity 104
    /// (`GeomToIGES_GeomCurve::TransferConic`); this port has no 104 emitter yet
    /// (T-78).
    Conic,
    /// Bezier / B-spline: OCCT emits entity 126
    /// (`GeomToIGES_GeomCurve::TransferBSplineCurve`), see
    /// [`IgesWriter::emit_bspline_curve`].
    Bounded,
    Other,
}

fn iges_curve_kind(c: &dyn Curve) -> IgCurveKind {
    if c.is_line() {
        IgCurveKind::Line
    } else if c.gp_circ().is_some() {
        IgCurveKind::Circle
    } else if c.gp_ellipse().is_some() || c.gp_hyperbola().is_some() || c.gp_parabola().is_some() {
        IgCurveKind::Conic
    } else if c.bspline_poles().is_some() || c.bezier_poles().is_some() {
        IgCurveKind::Bounded
    } else {
        IgCurveKind::Other
    }
}

/// `ArePolesPlanar` (`GeomToIGES_GeomCurve.cxx:170-199`): the area vector
/// `P(n) x P(1) + sum_{i<n} P(i) x P(i+1)`, normalised; every pole must then sit
/// at the same distance from the plane through `P(1)`. Fewer than three poles is
/// planar by definition and gets an arbitrary normal perpendicular to the first
/// segment (`GetAnyNormal`).
fn poles_planar_and_normal(poles: &[GpPnt]) -> (bool, GpVec) {
    let xyz = |p: &GpPnt| p.coord;
    if poles.len() < 3 {
        let d = xyz(&poles[0]).subtracted(&xyz(&poles[1]));
        return (true, any_normal(&d));
    }
    let mut n = xyz(&poles[poles.len() - 1]).crossed(&xyz(&poles[0]));
    for i in 0..poles.len() - 1 {
        n = n.added(&xyz(&poles[i]).crossed(&xyz(&poles[i + 1])));
    }
    let modulus = n.modulus();
    if modulus < occt_core::precision::CONFUSION {
        return (false, GpVec::new(0.0, 0.0, 1.0));
    }
    let normal = GpVec::from_xyz(&n.divided(modulus));
    let scl = xyz(&poles[0]).dot(&normal.xyz());
    for p in &poles[1..] {
        if (xyz(p).dot(&normal.xyz()) - scl).abs() > occt_core::precision::CONFUSION {
            return (false, normal);
        }
    }
    (true, normal)
}

/// `GetAnyNormal` (`GeomToIGES_GeomCurve.cxx:120-140`): a unit vector normal to
/// `v`, built from the axis least aligned with it.
fn any_normal(v: &occt_core::gp::GpXyz) -> GpVec {
    let ax = v.x().abs();
    let ay = v.y().abs();
    let az = v.z().abs();
    let out = if ax <= ay && ax <= az {
        GpVec::new(0.0, -v.z(), v.y())
    } else if ay <= az {
        GpVec::new(v.z(), 0.0, -v.x())
    } else {
        GpVec::new(-v.y(), v.x(), 0.0)
    };
    let m = out.xyz().modulus();
    if m < 1e-30 {
        GpVec::new(0.0, 0.0, 1.0)
    } else {
        GpVec::new(out.x() / m, out.y() / m, out.z() / m)
    }
}

/// The `Geom_SweptSurface` family of `GeomToIGES_GeomSurface::TransferSurface`
/// (`GeomToIGES_GeomSurface.cxx:1000-1025`): a linear extrusion becomes entity
/// 122, a revolution entity 120.
#[derive(Clone, Copy, PartialEq, Eq)]
enum SweptKind {
    Extrusion,
    Revolution,
}

/// Resolve the swept surface behind a face's surface the way
/// `GeomToIGES_GeomSurface::TransferSurface` dispatches: a
/// `Geom_RectangularTrimmedSurface` is peeled first (`cxx:492-515`, the
/// `Bounded` branch at `:177-181` recurses on `BasisSurface`), then the exact
/// `Geom_SweptSurface` type decides (`cxx:1013-1022`, extrusion before
/// revolution). Anything else - including a `Geom_OffsetSurface`, which is not a
/// `Geom_SweptSurface` - has no swept branch.
fn swept_surface_kind(s: &dyn Surface) -> Option<SweptKind> {
    if let Some(b) = s.rectangular_trimmed_basis() {
        return swept_surface_kind(b.as_ref());
    }
    if s.is_surface_of_linear_extrusion() {
        Some(SweptKind::Extrusion)
    } else if s.is_surface_of_revolution() {
        Some(SweptKind::Revolution)
    } else {
        None
    }
}

/// `Geom_BSplineCurve::IsEqual` (`Geom_BSplineCurve_1.cxx:662-...`) between the
/// two iso curves that `Geom_BSplineSurface::IsUClosed`/`IsVClosed`
/// (`Geom_BSplineSurface_1.cxx:1350-1389`) compares: the pole lists must agree
/// **component-wise** within `Precision::Confusion()` (not by point distance), and
/// every other quantity `IsEqual` looks at (degree, knot count and values,
/// multiplicity list, rational flag, pole count) is shared by construction,
/// because both iso curves are built from the same surface
/// (`Geom_BSplineSurface::UIso`/`VIso`) and therefore use the same V/U knots,
/// multiplicities, degree and periodicity. Differing lengths fail the comparison,
/// as `IsEqual` does on a count mismatch.
fn iso_rows_equal(a: &[GpPnt], b: &[GpPnt]) -> bool {
    let tol = occt_core::precision::CONFUSION;
    a.len() == b.len()
        && a.iter().zip(b).all(|(p, q)| {
            (p.x() - q.x()).abs() <= tol && (p.y() - q.y()).abs() <= tol && (p.z() - q.z()).abs() <= tol
        })
}

/// The rational arm of `Geom_BSplineCurve::IsEqual`
/// (`Geom_BSplineCurve_1.cxx`, `fabs(w1 - w2) > Epsilon(w1)`): weights are
/// compared against `Epsilon` of the **first** curve's weight
/// (`Standard_Real.hxx:242-246`).
fn weight_equal(a: f64, b: f64) -> bool {
    (a - b).abs() <= occt_core::precision::epsilon(a)
}

/// The unit of the model space the Global section declares
/// (`IGESData_IGESWriter.cxx:38-45` writes `2Hmm` with unit flag 2 = millimetre,
/// so `GeomToIGES_GeomEntity::GetUnit()` is 1: distances are written as they are).
const IGES_UNIT: f64 = 1.0;

/// `IGESConvGeom_GeomBuilder::IsIdentity` (`IGESConvGeom_GeomBuilder.cxx:159-175`,
/// `epsl = epsa = 1.E-10` at `:29-30`): the frame's 3x3 part is the identity
/// within `epsa` and its translation part vanishes within `epsl`. In OCCT `thepos`
/// is the frame's **local → global** transformation (`SetPosition(pos)` stores
/// `Trsf(pos, gp::XOY())`, `:137-143`) and `EvalXYZ` applies its inverse
/// (`:212-216`), so `thepos` is the identity exactly when the frame is the
/// absolute one.
fn frame_is_identity(frame: &GpAx3) -> bool {
    let eps = 1e-10;
    let (x, y, z) = (
        *frame.x_direction().xyz(),
        *frame.y_direction().xyz(),
        *frame.direction().xyz(),
    );
    let o = frame.location().coord;
    let m = [
        [x.x(), y.x(), z.x()],
        [x.y(), y.y(), z.y()],
        [x.z(), y.z(), z.z()],
    ];
    for (i, row) in m.iter().enumerate() {
        for (j, v) in row.iter().enumerate() {
            let cons = if i == j { 1.0 } else { 0.0 };
            if (*v - cons).abs() > eps {
                return false;
            }
        }
    }
    o.x().abs() <= eps && o.y().abs() <= eps && o.z().abs() <= eps
}

/// `IGESConvGeom_GeomBuilder::EvalXYZ` (`IGESConvGeom_GeomBuilder.cxx:212-216`):
/// the point expressed in the frame's local coordinates (the frame's axes are
/// orthonormal, so this is `thepos.Inverted()` applied to the point).
fn frame_local_point(frame: &GpAx3, p: &GpPnt) -> (f64, f64, f64) {
    let d = p.coord.subtracted(&frame.location().coord);
    (
        d.dot(frame.x_direction().xyz()),
        d.dot(frame.y_direction().xyz()),
        d.dot(&frame.direction().xyz()),
    )
}

/// `gp_Elips2d::Coefficients` (`gp_Elips2d.cxx:25-62`) evaluated on the frame
/// `gp_Ax22d(gp::Origin2d(), gp::DX2d(), gp::DY2d())` that
/// `GeomToIGES_GeomCurve.cxx:671` builds. With that frame the `gp_Trsf2d T` there
/// is the identity, so the general sums collapse to
/// `A = 1/DMaj, B = 1/DMin, C = D = E = 0, F = -1`; the degenerate arm
/// (`DMin <= gp::Resolution()`, `:35-45`) keeps its own form. The result is in
/// `gp_Elips2d::Coefficients`' **own** output order — the caller in
/// `GeomToIGES_GeomCurve.cxx:675` passes it in the shuffled order `(A, C, B, D,
/// E, F)`.
fn gp_elips2d_coefficients(major: f64, minor: f64) -> (f64, f64, f64, f64, f64, f64) {
    let dmin = minor * minor;
    let dmaj = major * major;
    if dmin <= occt_core::precision::REAL_SMALL && dmaj <= occt_core::precision::REAL_SMALL {
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    } else if dmin <= occt_core::precision::REAL_SMALL {
        (1.0, 0.0, 0.0, 0.0, 0.0, -dmaj)
    } else {
        (1.0 / dmaj, 1.0 / dmin, 0.0, 0.0, 0.0, -1.0)
    }
}

/// `gp_Hypr2d::Coefficients` (`gp_Hypr2d.cxx:25-61`) on the same identity frame:
/// `A = 1/DMaj, B = -1/DMin, C = D = E = 0, F = -1`, with the same degenerate arm.
fn gp_hypr2d_coefficients(major: f64, minor: f64) -> (f64, f64, f64, f64, f64, f64) {
    let dmin = minor * minor;
    let dmaj = major * major;
    if dmin <= occt_core::precision::REAL_SMALL && dmaj <= occt_core::precision::REAL_SMALL {
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    } else if dmin <= occt_core::precision::REAL_SMALL {
        (1.0, 0.0, 0.0, 0.0, 0.0, -dmaj)
    } else {
        (1.0 / dmaj, -1.0 / dmin, 0.0, 0.0, 0.0, -1.0)
    }
}

/// `gp_Parab2d::Coefficients` (`gp_Parab2d.cxx:44-60`) on the same identity frame,
/// with `P = 2 * focalLength` and `focalLength = start->Focal() * 2`
/// (`GeomToIGES_GeomCurve.cxx:818-819`): `A = C = 0`, `B = 1`, `D = -P`,
/// `E = F = 0`.
fn gp_parab2d_coefficients(focal: f64) -> (f64, f64, f64, f64, f64, f64) {
    let p = 2.0 * (focal * 2.0);
    (0.0, 1.0, 0.0, -p, 0.0, 0.0)
}

/// `GeomToIGES_GeomCurve.cxx:722-729` (hyperbola) and `:795-802` (parabola): a
/// range end past `Precision::Infinite()` becomes ±`Precision::Infinite()`. The
/// ellipse branch has no such arm.
fn occt_infinite_range(a: f64, b: f64) -> (f64, f64) {
    let u1 = if occt_core::precision::Precision::is_negative_infinite(a) {
        -occt_core::precision::INFINITE
    } else {
        a
    };
    let u2 = if occt_core::precision::Precision::is_positive_infinite(b) {
        occt_core::precision::INFINITE
    } else {
        b
    };
    (u1, u2)
}

/// `IGESGeom_ConicArc::ComputedFormNumber` (`IGESGeom_ConicArc.cxx:97-124`):
/// form 1 ellipse, 2 hyperbola, 3 parabola, 0 when the coefficients identify no
/// conic. `eps = 1.E-08` there, `eps2 = eps*eps`, `eps4 = eps2*eps2`.
fn conic_form_number(a: f64, b: f64, c: f64, d: f64, e: f64, f: f64) -> i32 {
    let eps = 1e-8f64;
    let eps4 = (eps * eps) * (eps * eps);
    let q1 = a * (c * f - e * e / 4.0) + b / 2.0 * (e * d / 4.0 - b * f / 2.0)
        + d / 2.0 * (b * e / 4.0 - c * d / 2.0);
    let q2 = a * c - b * b / 4.0;
    let q3 = a + c;
    if q2 > eps4 && q1 * q3 < 0.0 {
        1
    } else if q2 < -eps4 && q1.abs() > eps4 {
        2
    } else if q2.abs() <= eps4 && q1.abs() > eps4 {
        3
    } else {
        0
    }
}

// ---------------------------------------------------------------------------
// Entity bookkeeping
// ---------------------------------------------------------------------------

/// One IGES entity: directory-entry metadata + parameter-data string.
struct Ent {
    ty: i32,
    form: i32,
    /// Pointer to the entity's transformation matrix (`#124`), `HasTransf()` in
    /// `IGESData_IGESEntity`; written to DE card 1 field 7
    /// (`IGESData_IGESWriter.cxx:324-331`, `v[6] = themodel->DNum(...)`).
    trsf: Option<usize>,
    /// Parameter data, always ending with the record delimiter `;`.
    params: String,
}

impl Ent {
    fn param_line_count(&self) -> usize {
        self.params.len().div_ceil(64).max(1)
    }

    fn param_chunks(&self) -> Vec<String> {
        if self.params.is_empty() {
            return vec![String::new()];
        }
        let mut out = Vec::new();
        let mut start = 0;
        while start < self.params.len() {
            let end = (start + 64).min(self.params.len());
            out.push(self.params[start..end].to_string());
            start = end;
        }
        out
    }

    /// The two Directory Entry cards of the entity
    /// (`IGESData_IGESWriter.cxx:276-365` filling `v[0..16]`, `:809-880` laying
    /// the cards out): eight 8-column fields plus four 2-column status fields on
    /// the first card and nine 8-column fields on the second - 72 data columns
    /// each. Card 1 field 2 is the sequence number of the entity's **first
    /// parameter card** (`v[1] = thepnum.Value(i)`, `cxx:834`), field 7 the
    /// transformation matrix pointer (`cxx:324-331`); card 2 field 4 is the number
    /// of parameter cards (`v[15]`, `cxx:835`), field 5 the form number
    /// (`v[16] = anent->FormNumber()`, `cxx:364`).
    ///
    /// The remaining fields carry `IGESData_IGESEntity`'s defaults: line weight
    /// `theLWeightNum = 0` and color `DefColor() == DefVoid` ⇒ 0
    /// (`IGESData_IGESEntity.cxx:53-63`, `IGESData_IGESWriter.cxx:347-361`), a
    /// blank label (`theShortLabel` is null unless `SetShortLabel` is called,
    /// `cxx:366-379`) and the subscript number `theSubScriptN = 0`
    /// (`IGESData_IGESEntity.cxx:60`, written right-justified by `cxx:380-391`).
    fn directory_lines(&self, p_start: usize) -> (String, String) {
        let ty = field8(&self.ty.to_string());
        let pstart = field8(&p_start.to_string());
        let pcount = field8(&self.param_line_count().to_string());
        let form = field8(&self.form.to_string());
        let trsf = field8(&self.trsf.map_or("0".to_string(), |t| t.to_string()));
        let l1 = format!(
            "{ty}{pstart}{}{}{}{}{trsf}{}{:>2}{:>2}{:>2}{:>2}",
            field8("0"),
            field8("0"),
            field8("0"),
            field8("0"),
            field8("0"),
            0,
            0,
            0,
            0
        );
        let l2 = format!(
            "{ty}{}{}{pcount}{form}{}{}{}{}",
            field8("0"),
            field8("0"),
            field8(""),
            field8(""),
            field8(""),
            field8("0")
        );
        (l1, l2)
    }
}

/// Incremental IGES writer holding the entity list and dedup maps.
struct IgesWriter {
    entities: Vec<Ent>,
    point_entities: HashMap<usize, usize>,
    edge_curve_entities: HashMap<usize, usize>,
}

impl IgesWriter {
    fn new() -> Self {
        Self {
            entities: Vec::new(),
            point_entities: HashMap::new(),
            edge_curve_entities: HashMap::new(),
        }
    }

    /// Register one entity. The DE label, line weight, color and status fields are
    /// not parameters here: OCCT's writer takes them from the entity's own
    /// `IGESData_IGESEntity` state, which is at its defaults for generated
    /// geometry (`IGESData_IGESWriter.cxx:276-365`), see
    /// [`Ent::directory_lines`].
    fn emit(&mut self, ty: i32, form: i32, params: String) -> usize {
        self.entities.push(Ent {
            ty,
            form,
            trsf: None,
            params,
        });
        self.entities.len()
    }

    /// `IGESData_IGESEntity::InitTransf` — record the entity's transformation
    /// matrix pointer, written to DE card 1 field 7 (`IGESData_IGESWriter.cxx:324-331`).
    fn set_trsf(&mut self, de: usize, trsf: usize) {
        self.entities[de - 1].trsf = Some(trsf);
    }

    // ---- geometry entities ----

    fn emit_point(&mut self, p: &GpPnt) -> usize {
        self.emit(
            116,
            0,
            format!("116,{},{},{};", num(p.x()), num(p.y()), num(p.z())),
        )
    }

    fn emit_line(&mut self, a: &GpPnt, b: &GpPnt) -> usize {
        self.emit(
            110,
            0,
            format!(
                "110,{},{},{},{},{},{};",
                num(a.x()),
                num(a.y()),
                num(a.z()),
                num(b.x()),
                num(b.y()),
                num(b.z())
            ),
        )
    }

    fn emit_plane(&mut self, pln: &GpPln) -> usize {
        let n = *pln.axis().direction().xyz();
        let loc = pln.location().coord;
        let d = n.dot(&loc);
        self.emit(
            108,
            1,
            format!("108,{},{},{},{};", num(n.x), num(n.y), num(n.z), num(d)),
        )
    }

    fn emit_circular_arc(&mut self, center: &GpPnt, start: &GpPnt, end: &GpPnt, plane_pt: &GpPnt) -> usize {
        self.emit(
            100,
            0,
            format!(
                "100,0.,{},{},{},{},{},{},{},{},{},{},{},{};",
                num(center.x()),
                num(center.y()),
                num(center.z()),
                num(start.x()),
                num(start.y()),
                num(start.z()),
                num(end.x()),
                num(end.y()),
                num(end.z()),
                num(plane_pt.x()),
                num(plane_pt.y()),
                num(plane_pt.z())
            ),
        )
    }

    /// Arc through three points: center from the circumcircle, start/end from
    /// the first/last, and a fourth point on the arc plane.
    fn emit_arc_3p(&mut self, p1: &GpPnt, p2: &GpPnt, p3: &GpPnt) -> usize {
        let center = circle_center3(p1, p2, p3).unwrap_or_else(GpPnt::zero);
        let n = GpVec::from_pnts(p1, p2).xyz().crossed(GpVec::from_pnts(p1, p3).xyz());
        let plane_pt = if n.square_modulus() < 1e-30 {
            GpPnt::new(center.x(), center.y(), center.z() + 1.0)
        } else {
            GpPnt::from_xyz(&center.coord.added(&n.divided(n.modulus())))
        };
        self.emit_circular_arc(&center, p1, p3, &plane_pt)
    }

    // ---- topology entities ----

    fn emit_vertices(&mut self, shape: &TopoShape) {
        for v in vertices_of(shape) {
            let key = Arc::as_ptr(&v.0.tshape) as usize;
            if self.point_entities.contains_key(&key) {
                continue;
            }
            let p = BRepTool::vertex_point(&v);
            let idx = self.emit_point(&p);
            self.point_entities.insert(key, idx);
        }
    }

    /// `GeomToIGES_GeomCurve::TransferCurve(Geom_Curve, Udeb, Ufin)`
    /// (`GeomToIGES_GeomCurve.cxx:94-126` dispatching, `:133-161` on the bounded
    /// family, `:481-528` on the conics): a `Geom_Line` becomes entity 110, a
    /// `Geom_Circle` entity 100, a B-spline/Bezier entity 126. Returns `None`
    /// where OCCT's transfer returns a null handle (conics: the 104 writer is not
    /// ported yet; any other curve type has no branch at all).
    fn emit_curve_range(&mut self, curve: &dyn Curve, a: f64, b: f64) -> Option<usize> {
        match iges_curve_kind(curve) {
            IgCurveKind::Line => {
                let p1 = curve.d0(a);
                let p2 = curve.d0(b);
                Some(self.emit_line(&p1, &p2))
            }
            IgCurveKind::Circle => {
                let p1 = curve.d0(a);
                let p2 = curve.d0(b);
                // `GeomToIGES_GeomCurve::TransferCircle`
                // (`GeomToIGES_GeomCurve.cxx:292-329`) takes the centre, the axis
                // and the radius from the `Geom_Circle` itself - never from a
                // three-point fit on sampled points. A closed edge (start == end)
                // is written with the arc's start/end coincident (IGES closed
                // circular arc form).
                let (center, plane_pt) = match curve.gp_circ() {
                    Some(c) => {
                        let pos = c.position();
                        let loc = pos.location();
                        let n = pos.direction();
                        (
                            loc,
                            GpPnt::new(loc.x() + n.x(), loc.y() + n.y(), loc.z() + n.z()),
                        )
                    }
                    None => (GpPnt::zero(), GpPnt::new(0.0, 0.0, 1.0)),
                };
                Some(self.emit_circular_arc(&center, &p1, &p2, &plane_pt))
            }
            // UNPORTED (audit A26 / task T-78): `TransferCurve(Geom_BSplineCurve)`
            // (`cxx:279-423`) first makes a periodic curve non-periodic
            // (`SetNotPeriodic`) and calls `Segment` when the requested range is
            // narrower than the curve's own; this port has neither, so those cases
            // return `None`.
            IgCurveKind::Bounded => self.emit_bspline_curve(curve, a, b),
            // `TransferConic` (`GeomToIGES_GeomCurve.cxx:533-603` is the circle;
            // the ellipse/hyperbola/parabola transfers are `:608-700`, `:707-773`,
            // `:780-845`) writes entity 104. Any other curve type has no branch in
            // `TransferCurve` and yields a null handle (UNPORTED, audit A26 /
            // task T-78).
            IgCurveKind::Conic => self.emit_conic_arc(curve, a, b),
            IgCurveKind::Other => None,
        }
    }

    fn emit_edge_curve(&mut self, e: &Edge) -> usize {
        let Some(curve) = BRepTool::edge_curve(e) else {
            let (v1, v2) = edge_vertices(e);
            let p1 = v1.map(|v| BRepTool::vertex_point(&v)).unwrap_or_default();
            let p2 = v2.map(|v| BRepTool::vertex_point(&v)).unwrap_or_default();
            return self.emit_line(&p1, &p2);
        };
        let (a, b) = BRepTool::edge_parameters(e);
        match self.emit_curve_range(curve.as_ref(), a, b) {
            Some(idx) => idx,
            None => {
                // The transfer above returned a null handle: OCCT writes no curve
                // at all for those types, the port keeps the chord-line stand-in
                // (UNPORTED, audit A26 / task T-78, see `emit_curve_range`).
                let p1 = curve.d0(a);
                let p2 = curve.d0(b);
                self.emit_line(&p1, &p2)
            }
        }
    }

    /// Entity 126 (`GeomToIGES_GeomCurve::TransferCurve(Geom_BSplineCurve)`,
    /// `GeomToIGES_GeomCurve.cxx:279-423`, written by
    /// `IGESGeom_ToolBSplineCurve::WriteOwnParams`,
    /// `IGESGeom_ToolBSplineCurve.cxx:64-105`):
    ///
    /// `126, index, degree, planar, closed, polynomial, periodic,
    ///   knots[-degree .. index+1], weights[0 .. index], poles[0 .. index],
    ///   UMin, UMax, Normal;`
    ///
    /// with `index = nb_poles - 1` and the knot array the curve's flattened knot
    /// sequence. The planar flag and the normal come from `ArePolesPlanar`
    /// (`cxx:170-199`): `P(n) x P(1) + sum P(i) x P(i+1)`, normalised (or
    /// `(0,0,1)` with `planar = false` when the sum is shorter than
    /// `Precision::Confusion()`), then every pole must sit at the same distance
    /// from the plane through `P(1)`. Returns `None` for the cases this port
    /// cannot express.
    fn emit_bspline_curve(&mut self, curve: &dyn Curve, first: f64, last: f64) -> Option<usize> {
        // `cxx:294-307`: a periodic curve is converted to a non-periodic copy
        // before writing.
        if curve.is_periodic() {
            return None;
        }
        let (curve_first, curve_last) = (curve.first_parameter(), curve.last_parameter());
        // `cxx:320-356`: a narrower range is obtained with
        // `Geom_BSplineCurve::Segment`.
        if first > curve_first + occt_core::precision::PCONFUSION
            || last < curve_last - occt_core::precision::PCONFUSION
        {
            return None;
        }

        let (poles, knots, degree) = if let (Some(p), Some(k), Some(d)) = (
            curve.bspline_poles(),
            curve.bspline_knots(),
            curve.nurbs_degree(),
        ) {
            (p.to_vec(), k.to_vec(), d)
        } else if let Some(p) = curve.bezier_poles() {
            // `TransferCurve(Geom_BezierCurve)` (`cxx:430-450`) goes through
            // `GeomConvert::CurveToBSplineCurve`: one span whose end knots carry
            // multiplicity degree+1.
            let d = p.len().checked_sub(1)?;
            let mut k = vec![curve_first; d + 1];
            k.extend(std::iter::repeat(curve_last).take(d + 1));
            (p.to_vec(), k, d)
        } else {
            return None;
        };
        if poles.len() < 2 || degree == 0 || knots.len() != poles.len() + degree + 1 {
            return None;
        }
        let index = poles.len() - 1;

        let (planar, normal) = poles_planar_and_normal(&poles);
        // `cxx:415-418`: the normal is flipped when it points down.
        let normal = if normal.z() < 0.0 {
            GpVec::new(-normal.x(), -normal.y(), -normal.z())
        } else {
            normal
        };

        let weights: Vec<f64> = match curve.bspline_weights() {
            Some(w) if w.len() == poles.len() => w.to_vec(),
            _ => vec![1.0; poles.len()],
        };
        // `cxx:360`: `IPolyn = !IsRational()`.
        let polynomial = curve.bspline_weights().is_none();
        let closed = poles
            .first()
            .zip(poles.last())
            .map(|(f, l)| f.distance(l) <= occt_core::precision::CONFUSION)
            .unwrap_or(false);

        let mut s = format!(
            "126,{index},{degree},{},{},{},0",
            u8::from(planar),
            u8::from(closed),
            u8::from(polynomial)
        );
        for k in &knots {
            s.push(',');
            s.push_str(&num(*k));
        }
        for w in &weights {
            s.push(',');
            s.push_str(&num(*w));
        }
        for p in &poles {
            s.push_str(&format!(",{},{},{}", num(p.x()), num(p.y()), num(p.z())));
        }
        s.push_str(&format!(
            ",{},{},{},{},{};",
            num(first),
            num(last),
            num(normal.x()),
            num(normal.y()),
            num(normal.z())
        ));
        Some(self.emit(126, 0, s))
    }

    /// Entity 123 Direction (`IGESGeom_ToolDirection::WriteOwnParams`,
    /// `IGESGeom_ToolDirection.cxx:64-74`): the three components of the unit
    /// vector.
    fn emit_direction(&mut self, d: &GpDir) -> usize {
        self.emit(
            123,
            0,
            format!("123,{},{},{};", num(d.x()), num(d.y()), num(d.z())),
        )
    }

    /// Entity 196 (`GeomToIGES_GeomSurface::TransferSphericalSurface`,
    /// `GeomToIGES_GeomSurface.cxx:1366-1400`, written by
    /// `IGESSolid_ToolSphericalSurface::WriteOwnParams`,
    /// `IGESSolid_ToolSphericalSurface.cxx:70-80`):
    /// `196, centre_point, radius, axis_direction, reference_direction;` — the
    /// parametrised form, which is the one OCCT's transfer produces.
    fn emit_spherical_surface(
        &mut self,
        center: &GpPnt,
        radius: f64,
        axis: &GpDir,
        x_dir: &GpDir,
    ) -> usize {
        let c = self.emit_point(center);
        let a = self.emit_direction(axis);
        let r = self.emit_direction(x_dir);
        self.emit(
            196,
            0,
            format!("196,{c},{},{a},{r};", num(radius)),
        )
    }

    /// Entity 192 (`GeomToIGES_GeomSurface::TransferCylindricalSurface`,
    /// `GeomToIGES_GeomSurface.cxx:1280-1314`, written by
    /// `IGESSolid_ToolCylindricalSurface::WriteOwnParams`):
    /// `192, location_point, axis_direction, radius, reference_direction;`.
    fn emit_cylindrical_surface(
        &mut self,
        location: &GpPnt,
        axis: &GpDir,
        radius: f64,
        x_dir: &GpDir,
    ) -> usize {
        let l = self.emit_point(location);
        let a = self.emit_direction(axis);
        let r = self.emit_direction(x_dir);
        self.emit(
            192,
            0,
            format!("192,{l},{a},{},{r};", num(radius)),
        )
    }

    /// Entity 194 (`TransferConicalSurface`, `cxx:1318-1362`, written by
    /// `IGESSolid_ToolConicalSurface::WriteOwnParams`):
    /// `194, location_point, axis_direction, ref_radius, semi_angle_deg,
    /// reference_direction;`. A negative semi-angle is written by mirroring the
    /// reference point through the apex, negating the angle and reversing the
    /// reference direction (`cxx:1344-1350`).
    fn emit_conical_surface(
        &mut self,
        location: &GpPnt,
        apex: &GpPnt,
        axis: &GpDir,
        ref_radius: f64,
        semi_angle: f64,
        x_dir: &GpDir,
    ) -> usize {
        let (loc, angle, xd) = if semi_angle < 0.0 {
            (
                GpPnt::new(
                    2.0 * apex.x() - location.x(),
                    2.0 * apex.y() - location.y(),
                    2.0 * apex.z() - location.z(),
                ),
                -semi_angle,
                GpDir::new(-x_dir.x(), -x_dir.y(), -x_dir.z())
                    .unwrap_or_else(|_| *x_dir),
            )
        } else {
            (*location, semi_angle, *x_dir)
        };
        let l = self.emit_point(&loc);
        let a = self.emit_direction(axis);
        let r = self.emit_direction(&xd);
        self.emit(
            194,
            0,
            format!(
                "194,{l},{a},{},{},{r};",
                num(ref_radius),
                num(angle * 180.0 / std::f64::consts::PI)
            ),
        )
    }

    /// Entity 198 (`TransferToroidalSurface`, `cxx:1402-1435`, written by
    /// `IGESSolid_ToolToroidalSurface::WriteOwnParams`):
    /// `198, centre_point, axis_direction, major_radius, minor_radius,
    /// reference_direction;`.
    fn emit_toroidal_surface(
        &mut self,
        center: &GpPnt,
        axis: &GpDir,
        major: f64,
        minor: f64,
        x_dir: &GpDir,
    ) -> usize {
        let c = self.emit_point(center);
        let a = self.emit_direction(axis);
        let r = self.emit_direction(x_dir);
        self.emit(
            198,
            0,
            format!("198,{c},{a},{},{},{r};", num(major), num(minor)),
        )
    }

    /// Entity 128 (`GeomToIGES_GeomSurface::TransferBSplineSurface`,
    /// `GeomToIGES_GeomSurface.cxx:191-460`, written by
    /// `IGESGeom_ToolBSplineSurface::WriteOwnParams`,
    /// `IGESGeom_ToolBSplineSurface.cxx:64-125`):
    ///
    /// `128, indU, indV, degU, degV, closedU, closedV, polynomial, periodicU,
    ///   periodicV, knotU[-degU .. indU+1], knotV[-degV .. indV+1],
    ///   weights[0..indU][0..indV], poles[0..indU][0..indV],
    ///   UMin, UMax, VMin, VMax;`
    ///
    /// with `indU = nb_poles_u - 1`, `indV = nb_poles_v - 1` and the flattened
    /// knot vectors. A Bezier surface becomes the equivalent single-span
    /// B-spline (end knots with multiplicity degree+1) as
    /// `GeomConvert::SurfaceToBSplineSurface` would produce.
    ///
    /// `closedU/V` come from `IsUClosed`/`IsVClosed` (`cxx:347-348`, see
    /// [`iso_rows_equal`]), `polynomial` from `Polynom = !(RationU || RationV)`
    /// (`cxx:355`; the port's B-spline surface carries one rational flag for both
    /// directions) and the written range from the non-periodic arm of the bounds
    /// fix (`cxx:245-255`, `:274-284`). The periodic arm (`cxx:256-273`,
    /// `:285-302`, `:304-343`: `ShapeAnalysis::AdjustToPeriod`,
    /// `SetUOrigin`/`SetVOrigin`, `SetUNotPeriodic`/`SetVNotPeriodic` →
    /// `BSplSLib::Unperiodize`) is UNPORTED (audit A26 / task T-78): OCCT re-origins
    /// and unperiodizes the surface before reading its knots and poles and writes
    /// `periodicU/V` from the **original** `IsUPeriodic`/`IsVPeriodic`
    /// (`cxx:235-236`, `:448-449`), while this port keeps the periodic
    /// representation and writes `periodicU/V = 0`.
    fn emit_bspline_surface(
        &mut self,
        surf: &dyn Surface,
        u0: f64,
        u1: f64,
        v0: f64,
        v1: f64,
    ) -> Option<usize> {
        let (poles, knots_u, knots_v, deg_u, deg_v) = match (
            surf.bspline_surface_poles(),
            surf.bspline_surface_uknots(),
            surf.bspline_surface_vknots(),
        ) {
            (Some(p), Some(ku), Some(kv)) => (
                p.to_vec(),
                ku.to_vec(),
                kv.to_vec(),
                surf.u_degree().max(1) as usize,
                surf.v_degree().max(1) as usize,
            ),
            _ => return None,
        };
        if poles.is_empty() || poles[0].is_empty() {
            return None;
        }
        let (nu, nv) = (poles.len(), poles[0].len());
        if knots_u.len() != nu + deg_u + 1 || knots_v.len() != nv + deg_v + 1 {
            return None;
        }
        let (ind_u, ind_v) = (nu - 1, nv - 1);

        // `TransferSurface(Geom_RectangularTrimmedSurface)` recurses on
        // `BasisSurface()` before this branch (`cxx:492-515`), so every quantity
        // read below - poles, knots, degree, periodicity and `Bounds` - belongs to
        // the untrimmed basis (`st` there plays the role of `start` here).
        let basis = surf.rectangular_trimmed_basis();
        let bs: &dyn Surface = match basis.as_deref() {
            Some(b) => b,
            None => surf,
        };

        let stored_weights = bs.bspline_surface_weights();
        let rational = stored_weights.is_some();
        let weights: Vec<Vec<f64>> = match stored_weights {
            Some(w) if w.len() == nu && w.iter().all(|r| r.len() == nv) => w.to_vec(),
            _ => vec![vec![1.0; nv]; nu],
        };
        let polynomial = !rational;
        // `CloseU = mysurface->IsUClosed()` / `CloseV` (`cxx:347-348`).
        let closed_u = bs.is_u_periodic()
            || (iso_rows_equal(&poles[0], &poles[nu - 1])
                && (!rational || weights[0].iter().zip(&weights[nu - 1]).all(|(a, b)| weight_equal(*a, *b))));
        let closed_v = bs.is_v_periodic()
            || (0..nu).all(|i| {
                iso_rows_equal(&poles[i][0..1], &poles[i][nv - 1..nv])
                    && (!rational || weight_equal(weights[i][0], weights[i][nv - 1]))
            });

        // `cxx:244-255` / `:274-284`: a non-periodic surface clamps the written
        // range to its own bounds (the periodic arm is UNPORTED, see above).
        let (su0, su1) = bs.u_range();
        let (sv0, sv1) = bs.v_range();
        let (u0, u1) = if bs.is_u_periodic() {
            (u0, u1)
        } else {
            (u0.max(su0), u1.min(su1))
        };
        let (v0, v1) = if bs.is_v_periodic() {
            (v0, v1)
        } else {
            (v0.max(sv0), v1.min(sv1))
        };

        let mut s = format!(
            "128,{ind_u},{ind_v},{deg_u},{deg_v},{},{},{},0,0",
            u8::from(closed_u),
            u8::from(closed_v),
            u8::from(polynomial)
        );
        for k in &knots_u {
            s.push(',');
            s.push_str(&num(*k));
        }
        for k in &knots_v {
            s.push(',');
            s.push_str(&num(*k));
        }
        for j in 0..nv {
            for i in 0..nu {
                s.push(',');
                s.push_str(&num(weights[i][j]));
            }
        }
        for j in 0..nv {
            for i in 0..nu {
                let p = &poles[i][j];
                s.push_str(&format!(",{},{},{}", num(p.x()), num(p.y()), num(p.z())));
            }
        }
        s.push_str(&format!(
            ",{},{},{},{};",
            num(u0),
            num(u1),
            num(v0),
            num(v1)
        ));
        Some(self.emit(128, 0, s))
    }

    /// Entity 120 (`GeomToIGES_GeomSurface::TransferSurface(
    /// Geom_SurfaceOfRevolution)`, `GeomToIGES_GeomSurface.cxx:1111-1188`, written
    /// by `IGESGeom_ToolSurfaceOfRevolution::WriteOwnParams`, `:119-127`):
    ///
    /// `120, axis_line, generatrix, start_angle, end_angle;`
    ///
    /// - the generatrix is `GC.TransferCurve(BasisCurve, V1, V2)` (`cxx:1155`);
    /// - the axis is an `IGESGeom_Line` (`cxx:1163-1183`) whose start point is
    ///   the axis location and whose second point is `Location - Direction` -
    ///   the CAS.CADE axis reversed (`#30 rln`, `#36 BUC60328 face 7`), the same
    ///   convention the reader inverts (`IGESToBRep_TopoSurface.cxx:768`);
    /// - the angles come from `Surf->Init(Axis, Generatrix, 2*M_PI - U2,
    ///   2*M_PI - U1)` (`cxx:1185`), the CAS.CADE/IGES phase convention.
    fn emit_surface_of_revolution(
        &mut self,
        surf: &dyn Surface,
        u0: f64,
        u1: f64,
        v0: f64,
        v1: f64,
    ) -> Option<usize> {
        let basis = surf.revolution_basis_curve()?;
        let axis = surf.revolution_axis()?;
        // `GC.TransferCurve(Curve, V1, V2)` (`cxx:1155`); a generatrix type with
        // no IGES writer keeps the port's chord-line stand-in (UNPORTED, audit
        // A26 / task T-78: `TransferConic` entity 104).
        let generatrix = match self.emit_curve_range(basis.as_ref(), v0, v1) {
            Some(i) => i,
            None => {
                let p0 = basis.d0(v0);
                let p1 = basis.d0(v1);
                self.emit_line(&p0, &p1)
            }
        };
        let loc = axis.location();
        let d = axis.direction();
        let axis_line = self.emit_line(
            &loc,
            &GpPnt::new(loc.x() - d.x(), loc.y() - d.y(), loc.z() - d.z()),
        );
        let tau = 2.0 * std::f64::consts::PI;
        Some(self.emit(
            120,
            0,
            format!(
                "120,{axis_line},{generatrix},{},{};",
                num(tau - u1),
                num(tau - u0)
            ),
        ))
    }

    /// Entity 122 (`GeomToIGES_GeomSurface::TransferSurface(
    /// Geom_SurfaceOfLinearExtrusion)`, `GeomToIGES_GeomSurface.cxx:1032-1104`,
    /// written by `IGESGeom_ToolTabulatedCylinder::WriteOwnParams`, `:90-98`):
    ///
    /// `122, directrix, end_point.x, end_point.y, end_point.z;`
    ///
    /// - the U range is re-read from the surface's own `Bounds`
    ///   (`cxx:1067-1071`, the `OCC9490` fix), not from the caller's range;
    /// - the V range keeps OCCT's infinite handling (`cxx:1058-1065`);
    /// - the directrix is the basis curve translated so that its origin coincides
    ///   with the directrix origin (`cxx:1075-1096`): when
    ///   `|V1| > Precision::Confusion()` it is shifted by
    ///   `Value(U1,V1) - Value(U1,0)`, i.e. by `V1 * Direction`;
    /// - the terminate point is `start->Value(U1, V2)` (`cxx:1078`).
    fn emit_tabulated_cylinder(&mut self, surf: &dyn Surface, v0: f64, v1: f64) -> Option<usize> {
        let basis = surf.extrusion_basis_curve()?;
        let (u1, u2) = surf.u_range();
        let vv1 = if occt_core::precision::Precision::is_negative_infinite(v0) {
            -occt_core::precision::INFINITE
        } else {
            v0
        };
        let vv2 = if occt_core::precision::Precision::is_positive_infinite(v1) {
            occt_core::precision::INFINITE
        } else {
            v1
        };
        let end = surf.value(u1, vv2);
        let copy = if vv1.abs() > occt_core::precision::CONFUSION {
            let shift = GpVec::from_pnts(&surf.value(u1, 0.0), &surf.value(u1, vv1));
            basis.translated(&shift)
        } else {
            basis.clone_dyn()
        };
        let directrix = match self.emit_curve_range(copy.as_ref(), u1, u2) {
            Some(i) => i,
            None => {
                let p0 = copy.d0(u1);
                let p1 = copy.d0(u2);
                self.emit_line(&p0, &p1)
            }
        };
        Some(self.emit(
            122,
            0,
            format!(
                "122,{directrix},{},{},{};",
                num(end.x()),
                num(end.y()),
                num(end.z())
            ),
        ))
    }

    /// Entity 124 (`IGESConvGeom_GeomBuilder::MakeTransformation`,
    /// `IGESConvGeom_GeomBuilder.cxx:218-237`, written by
    /// `IGESGeom_ToolTransformationMatrix::WriteOwnParams`,
    /// `IGESGeom_ToolTransformationMatrix.cxx:90-104`): the frame's 3x4 matrix
    /// `R11 R12 R13 T1 R21 R22 R23 T2 R31 R32 R33 T3`, the translation column
    /// divided by the file unit; form 1 when the frame is left-handed
    /// (`rs->SetFormNumber(1)` when `thepos.IsNegative()`).
    fn emit_transformation_matrix(&mut self, frame: &GpAx3, unit: f64) -> usize {
        let x = *frame.x_direction().xyz();
        let y = *frame.y_direction().xyz();
        let z = *frame.direction().xyz();
        let o = frame.location().coord;
        let rows = [
            [x.x(), y.x(), z.x()],
            [x.y(), y.y(), z.y()],
            [x.z(), y.z(), z.z()],
        ];
        let form = if det3(&rows) < 0.0 { 1 } else { 0 };
        let mut s = String::from("124");
        for (i, row) in rows.iter().enumerate() {
            for v in row {
                s.push(',');
                s.push_str(&num(*v));
            }
            s.push(',');
            let t = match i {
                0 => o.x(),
                1 => o.y(),
                _ => o.z(),
            };
            s.push_str(&num(t / unit));
        }
        s.push(';');
        self.emit(124, form, s)
    }

    /// Entity 104 (`GeomToIGES_GeomCurve::TransferCurve(Geom_Ellipse)`,
    /// `GeomToIGES_GeomCurve.cxx:608-700`; `(Geom_Hyperbola)`, `:707-773`;
    /// `(Geom_Parabola)`, `:780-845`; written by
    /// `IGESGeom_ToolConicArc::WriteOwnParams`, `IGESGeom_ToolConicArc.cxx:103-125`):
    ///
    /// `104, A, B, C, D, E, F, ZT, start.x, start.y, end.x, end.y;`
    ///
    /// The coefficients come from the 2d conic built on the identity frame (see
    /// [`gp_elips2d_coefficients`]) with the radii divided by the file unit, the
    /// DE form number from [`conic_form_number`] (`IGESGeom_ConicArc::Init`,
    /// `IGESGeom_ConicArc.cxx:33-53`), `ZT = 0` and the arc's end points in the
    /// conic's own frame (`Build.EvalXYZ`, `:669-670`). A frame other than the
    /// absolute one is recorded as entity 124 on the DE card
    /// (`:692-697`). Returns `None` for the full-period ellipse, which OCCT routes
    /// through `GeomConvert_ApproxCurve` instead (`:620-645`).
    fn emit_conic_arc(&mut self, curve: &dyn Curve, a: f64, b: f64) -> Option<usize> {
        use std::f64::consts::PI;
        let unit = IGES_UNIT;
        let (abc, pos, u1, u2) = if let Some(e) = curve.gp_ellipse() {
            if (b - a - 2.0 * PI).abs() <= occt_core::precision::PCONFUSION {
                // UNPORTED (audit A26 / task T-78): `cxx:620-645` converts the
                // full-period ellipse with `GeomConvert_ApproxCurve` (then
                // `GeomConvert::CurveToBSplineCurve` + `Reparametrize`) and
                // transfers that B-spline; this port has no `GeomConvert_ApproxCurve`.
                return None;
            }
            // `cxx:649-654`: `|Udeb| <= gp::Resolution()` is snapped to 0.
            let u1 = if a.abs() <= occt_core::precision::REAL_SMALL {
                0.0
            } else {
                a
            };
            let (fa, fb, fc, fd, fe, ff) =
                gp_elips2d_coefficients(e.major_radius() / unit, e.minor_radius() / unit);
            // `E2d.Coefficients(A, C, B, D, E, F)` (`cxx:675`) then
            // `Init(A, 2*B, C, 2*D, 2*E, F)` (`cxx:678-686`).
            ([fa, 2.0 * fc, fb, 2.0 * fd, 2.0 * fe, ff], *e.position(), u1, b)
        } else if let Some(h) = curve.gp_hyperbola() {
            let (u1, u2) = occt_infinite_range(a, b);
            let (fa, fb, fc, fd, fe, ff) =
                gp_hypr2d_coefficients(h.major_radius / unit, h.minor_radius / unit);
            // `H2d.Coefficients(A, C, B, D, E, F)` (`cxx:749`) then
            // `Init(A, B, C, D, E, F)` (`cxx:751-759`; no doubling here).
            ([fa, fc, fb, fd, fe, ff], *h.position(), u1, u2)
        } else if let Some(p) = curve.gp_parabola() {
            let (u1, u2) = occt_infinite_range(a, b);
            let (fa, fb, fc, fd, fe, ff) = gp_parab2d_coefficients(p.focal / unit);
            // `P2d.Coefficients(A, C, B, D, E, F)` (`cxx:821`) then
            // `Init(A, B, C, D, E, F)` (`cxx:823-831`; no doubling here).
            ([fa, fc, fb, fd, fe, ff], *p.position(), u1, u2)
        } else {
            return None;
        };
        let frame = pos.to_ax3();
        let (sx, sy, _) = frame_local_point(&frame, &curve.d0(u1));
        let (ex, ey, _) = frame_local_point(&frame, &curve.d0(u2));
        let [fa, fb, fc, fd, fe, ff] = abc;
        let form = conic_form_number(fa, fb, fc, fd, fe, ff);
        let de = self.emit(
            104,
            form,
            format!(
                "104,{},{},{},{},{},{},{},{},{},{},{};",
                num(fa),
                num(fb),
                num(fc),
                num(fd),
                num(fe),
                num(ff),
                num(0.0),
                num(sx),
                num(sy),
                num(ex),
                num(ey)
            ),
        );
        // `cxx:688-697`: the transformation matrix is created **after** the conic
        // arc and recorded on it (`Conic->InitTransf(TMat)`).
        if !frame_is_identity(&frame) {
            let t = self.emit_transformation_matrix(&frame, unit);
            self.set_trsf(de, t);
        }
        Some(de)
    }

    /// Base surface entity for a face, plus any synthesized boundary curves
    /// (used when the face carries no boundary wires, e.g. a sphere).
    fn emit_face_surface(&mut self, f: &Face) -> (usize, Vec<usize>) {
        let Some(surf) = BRepTool::face_surface(f) else {
            return (self.emit_plane(&GpPln::default()), Vec::new());
        };
        if face_is_planar(f) {
            let pln = face_plane(f).unwrap_or_else(GpPln::default);
            return (self.emit_plane(&pln), Vec::new());
        }
        // `GeomToIGES_GeomSurface::TransferSurface` (`cxx:520-600`) dispatches on
        // the surface's exact type; the elementary surfaces become the IGES
        // surface entities 192/194/196/198 with the location point, the axis and
        // the reference direction taken from the `gp_*` surface itself.
        match classify_surface(surf.as_ref()) {
            SurfaceKind::Cylinder => {
                if let Some(cy) = surf.gp_cylinder() {
                    let pos = cy.position();
                    return (
                        self.emit_cylindrical_surface(
                            &pos.location(),
                            pos.axis().direction(),
                            cy.radius(),
                            pos.x_direction(),
                        ),
                        Vec::new(),
                    );
                }
            }
            SurfaceKind::Cone => {
                if let Some(co) = surf.gp_cone() {
                    let pos = co.position();
                    return (
                        self.emit_conical_surface(
                            &pos.location(),
                            &co.apex(),
                            pos.axis().direction(),
                            co.radius(),
                            co.semi_angle(),
                            pos.x_direction(),
                        ),
                        Vec::new(),
                    );
                }
            }
            SurfaceKind::Torus => {
                if let Some(t) = surf.gp_torus() {
                    let pos = t.position();
                    return (
                        self.emit_toroidal_surface(
                            &pos.location(),
                            pos.axis().direction(),
                            t.major_radius(),
                            t.minor_radius(),
                            pos.x_direction(),
                        ),
                        Vec::new(),
                    );
                }
            }
            _ => {}
        }
        if classify_surface(surf.as_ref()) == SurfaceKind::Sphere {
            if let Some(center) = sphere_center(surf.as_ref()) {
                let (u0, _, v0, v1) = surf_bounds(surf.as_ref());
                let vm = 0.5 * (v0 + v1);
                let r = surf.d0(u0, vm).distance(&center);
                if r > 1e-9 {
                    // `GeomToIGES_GeomSurface::TransferSphericalSurface`
                    // (`GeomToIGES_GeomSurface.cxx:1366-1400`): entity 196 is a
                    // centre **point** entity (#116), the radius, the axis
                    // direction (#123) and the reference direction (#123, the
                    // sphere's X axis), written by
                    // `IGESSolid_ToolSphericalSurface::WriteOwnParams` as
                    // `196, centre, radius, axis, refdir;`.
                    // `GeomToIGES_GeomSurface.cxx:1384-1393`: the axis and the
                    // reference direction come from the `gp_Sphere`'s position
                    // (its main axis and its X axis).
                    let (axis, x_dir) = match surf.gp_sphere() {
                        Some(sp) => {
                            let axis = *sp.position().axis().direction();
                            let x_dir = *sp.position().x_direction();
                            (axis, x_dir)
                        }
                        None => (
                            GpDir::new(0.0, 0.0, 1.0).expect("z"),
                            GpDir::new(1.0, 0.0, 0.0).expect("x"),
                        ),
                    };
                    let sph = self.emit_spherical_surface(&center, r, &axis, &x_dir);
                    // The two meridian arcs still serve as the face's seam
                    // boundary curves (#144 needs a boundary), matching the
                    // sphere the port's primitives build.
                    let south = GpPnt::new(center.x(), center.y(), center.z() - r);
                    let north = GpPnt::new(center.x(), center.y(), center.z() + r);
                    let gen = self.emit_arc_3p(
                        &south,
                        &GpPnt::new(center.x() + r, center.y(), center.z()),
                        &north,
                    );
                    let mer2 = self.emit_arc_3p(
                        &north,
                        &GpPnt::new(center.x() - r, center.y(), center.z()),
                        &south,
                    );
                    return (sph, vec![gen, mer2]);
                }
            }
        }
        // `GeomToIGES_GeomSurface::TransferSurface` (`cxx:520-600`) also has the
        // B-spline branch (`TransferBSplineSurface`), which the port now writes
        // as entity 128.
        let (u0, u1, v0, v1) = face_uv_bounds_finite(f);
        if let Some(idx) = self.emit_bspline_surface(surf.as_ref(), u0, u1, v0, v1) {
            return (idx, Vec::new());
        }
        // `GeomToIGES_GeomSurface::TransferSurface` (`cxx:1000-1025`) on the swept
        // family - reached after the bounded family (B-spline / Bezier / trimmed,
        // `cxx:125-139`), which is the order `cxx:112-147` dispatches in.
        match swept_surface_kind(surf.as_ref()) {
            Some(SweptKind::Extrusion) => {
                if let Some(idx) = self.emit_tabulated_cylinder(surf.as_ref(), v0, v1) {
                    return (idx, Vec::new());
                }
            }
            Some(SweptKind::Revolution) => {
                if let Some(idx) = self.emit_surface_of_revolution(surf.as_ref(), u0, u1, v0, v1) {
                    return (idx, Vec::new());
                }
            }
            None => {}
        }
        // Unclassified curved face: fall back to a plane.
        // ponytail: covers the surfaces the IGES writer has no emitter for.
        let pln = face_plane(f).unwrap_or_else(GpPln::default);
        (self.emit_plane(&pln), Vec::new())
    }

    /// A face as `BRepToIGES_BRShell::TransferFace` builds it
    /// (`BRepToIGES_BRShell.cxx:250-405`): the base surface entity, one
    /// `CurveOnSurface` (142) per wire, and a `TrimmedSurface` (144) that carries
    /// the outer/inner contour pointers.
    ///
    /// The 3-D curve of a wire comes from `BRepToIGES_BRWire::TransferWire`
    /// (`BRepToIGES_BRWire.cxx:662-781`): a single-edge wire is that edge's curve
    /// entity, a wire with two or more edges a `CompositeCurve` (102, form 0).
    /// The 2-D (UV) curve, which OCCT builds with
    /// `TransferEdge(edge, face, originMap, length, false)` (`:720`), is UNPORTED
    /// (audit A26 / task T-78) - the transfer therefore takes the "3-D only" arm
    /// of `BRepToIGES_BRShell.cxx:285-288`, writing a null `CurveUV` and
    /// `PreferenceMode = 2`. Edges of the face that belong to no wire
    /// (`:334-365`) are likewise UNPORTED.
    ///
    /// **Registered remainder (T-84)**: `510` still carries this port's earlier
    /// layout instead of `IGESSolid_ToolFace::WriteOwnParams`
    /// (`510, surface, nb_loops, has_outer_loop, loop_ptrs...`) over `508` Loop
    /// entities, and `514`/`186` follow the same old shape.
    fn emit_face(&mut self, f: &Face) -> usize {
        let (surf_idx, mut synth) = self.emit_face_surface(f);
        let wires = wires_of_face(f);
        // `ShapeAlgo::AlgoContainer()->OuterWire(aFace)` (`cxx:275`) resolves to
        // `ShapeAnalysis::OuterWire` (`ShapeAnalysis.cxx`: the last wire, or the
        // first one with a non-negative `TotCross2D`).
        let outer_wire = crate::meshing::model_builder::outer_of_wires(&wires, f);
        let mut curve_refs: Vec<usize> = Vec::new();
        let mut outer_curve: Option<usize> = None;
        let mut inner_curves: Vec<usize> = Vec::new();
        // The outer wire is transferred first (`cxx:276-294`), then the inner ones
        // (`:301-331`).
        let ordered: Vec<&crate::shape::Wire> = outer_wire
            .iter()
            .chain(wires.iter().filter(|w| {
                outer_wire
                    .as_ref()
                    .map_or(true, |o| Arc::as_ptr(&o.0.tshape) != Arc::as_ptr(&w.0.tshape))
            }))
            .collect();
        for w in ordered {
            let mut edge_refs: Vec<usize> = Vec::new();
            for e in edges_of_wire(w) {
                let key = Arc::as_ptr(&e.0.tshape) as usize;
                let idx = match self.edge_curve_entities.get(&key) {
                    Some(&i) => i,
                    None => {
                        let i = self.emit_edge_curve(&e);
                        self.edge_curve_entities.insert(key, i);
                        i
                    }
                };
                edge_refs.push(idx);
                curve_refs.push(idx);
            }
            if edge_refs.is_empty() {
                continue;
            }
            let curve3d = if edge_refs.len() == 1 {
                edge_refs[0]
            } else {
                let refs = edge_refs
                    .iter()
                    .map(|i| i.to_string())
                    .collect::<Vec<_>>()
                    .join(",");
                self.emit(102, 0, format!("102,{},{};", edge_refs.len(), refs))
            };
            // `IGESGeom_CurveOnSurface::Init` (`IGESGeom_CurveOnSurface.cxx:26-40`)
            // with `Imode = 0` (`cxx:269`) and the "3-D only" preference above;
            // `IGESGeom_ToolCurveOnSurface::WriteOwnParams` (`:104-115`) writes
            // `142, creation_mode, surface, curve_uv, curve_3d, preference_mode;`.
            let cs = self.emit(142, 0, format!("142,0,{surf_idx},0,{curve3d},2;"));
            match &outer_wire {
                Some(o) if Arc::as_ptr(&o.0.tshape) == Arc::as_ptr(&w.0.tshape) => {
                    outer_curve = Some(cs)
                }
                _ => inner_curves.push(cs),
            }
        }
        curve_refs.append(&mut synth);

        // `cxx:380-400`: `isWholeSurface` is `BRep_Tool::NaturalRestriction(face)`,
        // forced to false for a plane / cylinder / cone - the guard there tests the
        // `CurveOnSurface` handle, which is non-null whenever the face has a wire.
        let mut is_whole = BRepTool::natural_restriction(f);
        let surf = BRepTool::face_surface(f);
        if let Some(s) = surf.as_ref() {
            let k = classify_surface(s.as_ref());
            if (k == SurfaceKind::Plane || k == SurfaceKind::Cylinder || k == SurfaceKind::Cone)
                && outer_curve.is_some()
            {
                is_whole = false;
            }
        }
        let trim_idx = if curve_refs.is_empty() {
            surf_idx
        } else {
            let outer_flag = u8::from(!is_whole);
            let n_inner = inner_curves.len();
            // `IGESGeom_ToolTrimmedSurface::WriteOwnParams`
            // (`IGESGeom_ToolTrimmedSurface.cxx:196-218`): `144, surface,
            // outer_boundary_type, nb_inner_contours, outer_contour, inner...;`
            // with the outer contour written as `0` when the type is false.
            let mut refs = match (is_whole, outer_curve) {
                (false, Some(i)) => i.to_string(),
                _ => "0".to_string(),
            };
            for c in &inner_curves {
                refs.push_str(&format!(",{c}"));
            }
            self.emit(
                144,
                0,
                format!("144,{surf_idx},{outer_flag},{n_inner},{refs};"),
            )
        };

        let mut params = format!("510,{trim_idx}");
        for c in &curve_refs {
            params.push_str(&format!(",{c}"));
        }
        params.push(';');
        self.emit(510, 0, params)
    }

    fn emit_shell(&mut self, sh: &Shell) -> usize {
        let faces = children_of_type(&sh.0, ShapeType::Face);
        let refs: Vec<String> = faces
            .iter()
            .map(|f| self.emit_face(&Face(f.clone())).to_string())
            .collect();
        self.emit(514, 0, format!("514,{};", refs.join(",")))
    }

    fn emit_solid(&mut self, s: &Solid) -> usize {
        let shells = children_of_type(&s.0, ShapeType::Shell);
        let refs: Vec<String> = shells
            .iter()
            .map(|sh| self.emit_shell(&Shell(sh.clone())).to_string())
            .collect();
        let outer = refs.first().cloned().unwrap_or_else(|| "0".to_string());
        self.emit(186, 0, format!("186,{outer};"))
    }

    fn emit_shape(&mut self, shape: &TopoShape) {
        match shape.shape_type() {
            ShapeType::Compound => {
                let kids: Vec<TopoShape> = shape
                    .tshape
                    .read()
                    .unwrap()
                    .children
                    .clone();
                for k in kids {
                    self.emit_shape(&k);
                }
            }
            ShapeType::Solid => {
                self.emit_vertices(shape);
                self.emit_solid(&Solid(shape.clone()));
            }
            ShapeType::Shell => {
                self.emit_vertices(shape);
                self.emit_shell(&Shell(shape.clone()));
            }
            ShapeType::Face => {
                self.emit_vertices(shape);
                self.emit_face(&Face(shape.clone()));
            }
            _ => {}
        }
    }

    // ---- file assembly ----

    fn global_lines(&self) -> Vec<String> {
        let data = format!(
            "1H,,1H;,4HIGES,5Hmodel,7Hocctp2,8H20260731,32,38,6,308,15,14,1,2,2Hmm,1,1,8H20260731.01,1.0E-6,1.0E6,8Hocct-rs,4Hrust,3,1,8H20260731.01,0;"
        );
        let mut out = Vec::new();
        let mut s = data.as_str();
        while !s.is_empty() {
            let end = s.len().min(64);
            out.push(s[..end].to_string());
            s = &s[end..];
        }
        out
    }

    fn finish(self) -> String {
        let mut out = String::new();
        let mut s_seq = 1usize;
        let mut g_seq = 1usize;
        let mut d_seq = 1usize;
        let mut p_seq = 1usize;

        out.push_str(&sec_line(
            'S',
            s_seq,
            "IGES B-REP MODEL GENERATED BY RUST OCCT PORT",
            MAXCARS_G,
        ));
        out.push('\n');
        s_seq += 1;

        for gl in self.global_lines() {
            out.push_str(&sec_line('G', g_seq, &gl, MAXCARS_G));
            out.push('\n');
            g_seq += 1;
        }

        // Directory entries: entity `i` owns the two DE pointers `2i-1` and `2i`.
        let mut p_start = 1usize;
        for (i, ent) in self.entities.iter().enumerate() {
            let (l1, l2) = ent.directory_lines(p_start);
            p_start += ent.param_line_count();
            out.push_str(&sec_line('D', 2 * i + 1, &l1, MAXCARS_G));
            out.push('\n');
            out.push_str(&sec_line('D', 2 * i + 2, &l2, MAXCARS_G));
            out.push('\n');
            d_seq += 2;
        }

        for (i, ent) in self.entities.iter().enumerate() {
            for chunk in ent.param_chunks() {
                out.push_str(&param_line(&chunk, 2 * i + 1, p_seq));
                out.push('\n');
                p_seq += 1;
            }
        }

        // Terminate card (`IGESData_IGESWriter.cxx:942-943`): the last sequence
        // number of each section, blank-filled, then 40 blanks and `T0000001`.
        let t = format!(
            "S{:>7}G{:>7}D{:>7}P{:>7}{}T0000001",
            s_seq - 1,
            g_seq - 1,
            d_seq,
            p_seq - 1,
            " ".repeat(40)
        );
        out.push_str(&t);
        out.push('\n');
        out
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

/// Finite UV bounds of a face (derived from the boundary curves for planes).
fn face_uv_bounds_finite(f: &Face) -> (f64, f64, f64, f64) {
    if let Some(s) = BRepTool::face_surface(f) {
        crate::wireframe::face_uv_bounds(f, s.as_ref())
    } else {
        (0.0, 1.0, 0.0, 1.0)
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Serialize the whole model to an IGES file (ANSI Y14.26M).
pub fn write_iges(model: &BRepModel) -> String {
    let mut w = IgesWriter::new();
    for ms in &model.shapes {
        w.emit_shape(&ms.shape);
    }
    w.finish()
}

/// Serialize the model to an IGES file on disk.
pub fn write_iges_file(path: &str, model: &BRepModel) -> std::io::Result<()> {
    std::fs::write(path, write_iges(model))
}

/// Serialize a single shape (wrapped in a one-shape model).
pub fn write_shape_iges(shape: &TopoShape) -> String {
    let mut model = BRepModel::new();
    model.add("Shape", shape.clone());
    write_iges(&model)
}

/// Count `(geometry entities, total records)` for testing.
///
/// Geometry entities are directory entries whose type number falls in the
/// geometric range (100–199: point/line/arc/plane/revolution/trimmed/MSBO, or
/// 500–599: B-rep face/shell). Total records is the number of directory
/// entries (D-section lines / 2).
pub fn iges_entity_counts(text: &str) -> (usize, usize) {
    // The section letter sits in column 73 and the entity type in columns 1-8
    // (`IGESData_IGESWriter.cxx:836-859`).
    let d_lines: Vec<&str> = text
        .lines()
        .filter(|l| l.as_bytes().get(72) == Some(&b'D'))
        .collect();
    let total = d_lines.len() / 2;
    let mut geo = 0usize;
    for (i, l) in d_lines.iter().enumerate() {
        if i % 2 == 0 {
            let ty: i32 = l.get(0..8).unwrap_or("0").trim().parse().unwrap_or(0);
            if (100..200).contains(&ty) || (500..600).contains(&ty) {
                geo += 1;
            }
        }
    }
    (geo, total)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::{BRepPrimBox, BRepPrimSphere};

    #[test]
    fn box_iges_contains_expected_entities() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let iges = write_shape_iges(&b.solid.0);
        for needle in ["116", "110", "108", "510", "514", "186"] {
            assert!(iges.contains(needle), "missing entity {needle}");
        }
        for line in iges.lines() {
            assert_eq!(line.len(), 80, "line length {}: {:?}", line.len(), line);
        }
        let (geo, total) = iges_entity_counts(&iges);
        assert!(geo > 0 && total > 0, "entity counts ({geo}, {total})");
        let last = iges.lines().last().unwrap();
        assert!(last.contains("T0000001"), "terminate record: {last}");
    }

    #[test]
    fn sphere_iges_has_arc_and_solid() {
        let s = BRepPrimSphere::make_sphere(2.0);
        let iges = write_shape_iges(&s.solid.0);
        assert!(
            iges.contains("100") || iges.contains("128"),
            "expected a circular arc or NURBS curve"
        );
        assert!(iges.contains("186"), "missing 186");
        for line in iges.lines() {
            assert_eq!(line.len(), 80, "line length {}", line.len());
        }
    }

    #[test]
    fn iges_file_written() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let mut model = BRepModel::new();
        model.add("Box", b.solid.0.clone());
        let path = std::env::temp_dir().join("occt_iges_test.igs");
        let p = path.to_str().unwrap();
        write_iges_file(p, &model).expect("write iges file");
        let content = std::fs::read_to_string(p).expect("read iges file");
        // The card format puts the data in columns 1-72, the section letter in
        // column 73 and the sequence number in columns 74-80
        // (`IGESData_IGESWriter.cxx:763-793` for the Start section).
        let first = content.lines().next().expect("first card");
        assert_eq!(first.len(), 80, "card length: {first:?}");
        assert_eq!(
            first.as_bytes().get(72),
            Some(&b'S'),
            "first card must be a Start card: {first:?}"
        );
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn section_terminators_present() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let iges = write_shape_iges(&b.solid.0);
        let last = iges.lines().last().unwrap();
        assert!(
            last.contains('S')
                && last.contains('G')
                && last.contains('D')
                && last.contains('P')
                && last.contains('T'),
            "terminate record missing section markers: {last}"
        );
    }

    #[test]
    fn rec80_pads_and_truncates() {
        assert_eq!(rec80("abc").len(), 80);
        let long = "x".repeat(100);
        assert_eq!(rec80(&long).len(), 80);
    }
}
