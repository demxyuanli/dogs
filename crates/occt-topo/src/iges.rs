//! IGES (ANSI Y14.26M) B-Rep writer — a port of `IGESControl_Writer`.
//!
//! Emits an 80-column fixed-width IGES file with the entity set OCCT's default
//! Faces mode (`write.iges.brep.mode = 0`, `IGESData.cxx:90-94`) produces for a
//! B-rep: 116 POINT, 110 LINE, 100 CIRCULAR ARC, 108 PLANE, 120 SURFACE OF
//! REVOLUTION, 122 TABULATED CYLINDER, 124 TRANSFORMATION MATRIX, 126 B-SPLINE
//! CURVE, 128 B-SPLINE SURFACE, 142 CURVE ON SURFACE, 144 TRIMMED SURFACE and
//! 402 GROUP. Elementary surfaces (cylinder, cone, sphere, torus) become 120,
//! **not** the 192/194/196/198 solids — those are the BRep-mode arm
//! (`myBRepMode && myAnalytic`) and are not selected here.
//!
//! Curves and surfaces cannot be downcast from `Arc<dyn Curve>` /
//! `Arc<dyn Surface>`, so the `Surface` / `Curve` traits expose the exact
//! analytic type (`gp_pln`, `gp_cylinder`, … / `gp_circ`, …) and every dispatch
//! is a type query, mirroring `GeomToIGES_GeomSurface::TransferSurface`'s
//! `IsKind` chain. A face whose basis surface has no branch in that chain is
//! dropped (`BRepToIGES_BRShell.cxx:257-261`), never replaced by a plane.
//!
//! The physical layout is the standard sectioned format (S/G/D/P/T) with each
//! line exactly 80 characters. Section lines carry their section letter and a
//! 7-digit sequence number in the first 8 columns, followed by the data.
//!
//! **Model contents (T-85 step 2)**: only the entities reachable from the
//! shapes handed to the writer are written, numbered in
//! `Interface_InterfaceModel::AddWithRefs` order (`Interface_InterfaceModel.cxx:652-692`),
//! exactly as `IGESControl_Writer::AddEntity` → `myModel->AddWithRefs`
//! (`IGESControl_Writer.cxx:243-252`) does. Parameter pointers therefore cannot
//! be written as they are built: the emitters record them as `#k` placeholders
//! ([`emit_refs`](IgesWriter::emit_refs)) and
//! [`final_entities`](IgesWriter::final_entities) resolves them after the
//! renumbering.
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use occt_core::gp::{
    GpAx2, GpAx2d, GpAx22d, GpAx3, GpCirc, GpCirc2d, GpDir, GpDir2d, GpElips, GpElips2d, GpHypr,
    GpHypr2d, GpLin, GpLin2d, GpMat2d, GpParab, GpParab2d, GpPln, GpPnt, GpPnt2d, GpTrsf,
    GpTrsf2d, GpVec, GpVec2d, GpXY, GpXyz, TrsfForm,
};
use occt_geom::{
    Curve, GeomBSplineCurve, GeomBezierCurve, GeomCircle, GeomEllipse, GeomHyperbola, GeomLine,
    GeomParabola, Surface,
};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;
use occt_geom2d::{
    Geom2dBSplineCurve, Geom2dBezierCurve, Geom2dCircle, Geom2dEllipse, Geom2dHyperbola,
    Geom2dLine, Geom2dParabola,
};

use crate::abs::ShapeType;
use crate::brep_surface::{classify_surface, SurfaceKind};
use crate::brep_tool::BRepTool;
use crate::model::BRepModel;
use crate::shape::{Edge, Face, Shell, Solid, TopoShape, Vertex};
use crate::topo_tools_full::{edge_vertices, edges_of_wire, wires_of_face};

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
/// The sequence is `%7.7d` (`Sprintf(finlin, "G%7.7d", i)` and `D%7.7d`),
/// zero-filled, not blank-filled.
fn sec_line(section: char, seq: usize, content: &str, width: usize) -> String {
    let mut s = String::with_capacity(80);
    for ch in content.chars().take(width) {
        s.push(ch);
    }
    s.push_str(&" ".repeat(width - s.len()));
    s.push(section);
    s.push_str(&format!("{seq:07}"));
    s
}

/// One Parameter card (`IGESData_IGESWriter.cxx:902-925`):
/// `data` (64 columns) + blank + the owning entity's directory-entry pointer
/// (`2*i - 1`) + `P` + the card's sequence number.
/// `Sprintf(finlin, " %7.7dP%7.7d", 2*i-1, j)` (`cxx:903`).
fn param_line(content: &str, de_pointer: usize, seq: usize) -> String {
    let mut s = String::with_capacity(80);
    for ch in content.chars().take(MAXCARS_P) {
        s.push(ch);
    }
    s.push_str(&" ".repeat(MAXCARS_P - s.len()));
    s.push_str(&format!(" {de_pointer:07}P{seq:07}"));
    s
}

/// The `Interface_LineBuffer` a section writer fills
/// (`IGESData_IGESWriter.cxx:49-50`: `thecurr(MaxcarsG + 1)`, `:248`:
/// `thecurr.SetMax(MaxcarsP)` at the start of the parameter section): text is
/// appended while it fits on the current card, and a card that overflows is
/// dumped **before** the text that does not fit is written
/// (`Interface_LineBuffer.cxx:68-80` `CanGet`, `:150-157` `Keep`).
///
/// This is what decides where a Parameter card ends: not a fixed 64-column cut,
/// but the first field boundary that does not fit. `Send(v)`
/// (`IGESData_IGESWriter.cxx:560-565`, `:586-589`) writes the separator and the
/// value as two separate additions, so a card can also end on a separator whose
/// value went to the next card.
struct CardBuffer {
    /// Cards already dumped, i.e. `thepars` of `IGESData_IGESWriter`.
    cards: Vec<String>,
    /// The card being filled, i.e. `thecurr`.
    line: String,
    /// `thecurr`'s maximum length: `MaxcarsP` (64) in the parameter section,
    /// `MaxcarsG` (72) in the global one (`IGESData_IGESWriter.cxx:39-40`).
    max: usize,
}

impl CardBuffer {
    fn new(max: usize) -> Self {
        Self {
            cards: Vec::new(),
            line: String::new(),
            max,
        }
    }

    /// `Interface_LineBuffer::Moved` through `IGESData_IGESWriter`'s flush
    /// (`IGESData_IGESWriter.cxx:508-512`, `:541-550`): the current card is
    /// appended to the section **as is** — trailing blanks are added at printing
    /// time (`cxx:913`), and a flush of a short card keeps it short.
    fn flush(&mut self) {
        self.cards.push(std::mem::take(&mut self.line));
    }

    /// `IGESData_IGESWriter::AddString` (`cxx:497-533`). The `+ 1` of the
    /// `CanGet` call is OCCT's own reservation (comment at `cxx:504`): one
    /// character must remain free, so that the *next* separator is never the
    /// first character of a line.
    fn add_string(&mut self, text: &str) {
        let bytes = text.as_bytes();
        if self.line.len() + bytes.len() + 1 > self.max {
            self.flush();
        }
        // Text longer than a card fills whole cards (`cxx:516-531`); the buffer
        // is empty here, so each iteration dumps exactly `max` characters.
        let mut rest = bytes;
        while rest.len() > self.max {
            let n = self.max - self.line.len();
            self.line
                .push_str(&String::from_utf8_lossy(&rest[..n.min(rest.len())]));
            self.flush();
            rest = &rest[self.max..];
        }
        let n = (self.max - self.line.len()).min(rest.len());
        self.line.push_str(&String::from_utf8_lossy(&rest[..n]));
    }

    /// `IGESData_IGESWriter::AddChar` (`cxx:535-553`), used for the parameter
    /// separators and the record delimiter.
    fn add_char(&mut self, c: u8) {
        if self.line.len() + 1 > self.max {
            self.flush();
        }
        self.line.push(c as char);
    }
}

/// Right-justify a value in an 8-column IGES field.
fn field8(v: &str) -> String {
    format!("{v:>8}")
}

/// Every real parameter of every entity goes through
/// `IGESData_IGESWriter::Send(const double)` (`IGESData_IGESWriter.cxx:582-589`),
/// which is `thefloatw.Write(val, lval)`, and the writer constructs `thefloatw`
/// as `Interface_FloatWriter(9)` (`IGESData_IGESWriter.cxx:49-50`, copy `:67`).
/// `Interface_FloatWriter(9)` runs `SetDefaults(9)`
/// (`Interface_FloatWriter.cxx:23-27` -> `:50-70`), which fixes the two `Sprintf`
/// formats and the range they are chosen by:
///
/// - `"%11.9f"` while `|val|` is in `[0.1, 1000)`, `"%11.9E"` outside it
///   (`cxx:139-147`);
/// - `thezerosup = true`, so `Convert` then strips trailing zeros from the
///   mantissa and drops an `E+00` exponent (`cxx:148-174`).
///
/// `%11.9E`, unlike Rust's `{:.9E}`, always writes the exponent sign and at
/// least two exponent digits (`1.5E+3` against `1.500000000E+03`), and the zero
/// suppression then leaves the shortest form (`150.`, `-0.`, `2.842170943E-14`).
const FLOAT_RANGE_MIN: f64 = 0.1;
const FLOAT_RANGE_MAX: f64 = 1000.0;

/// `Sprintf(text, "%11.9f", val)` (`cxx:143`). The field width is a minimum and
/// `"0.000000000"` already fills 11 columns, so no padding can ever appear.
fn float_fixed(x: f64) -> String {
    format!("{x:.9}")
}

/// `Sprintf(text, "%11.9E", val)` (`cxx:145`), re-normalised to C's `%E`
/// layout: sign always present, exponent padded to two digits.
fn float_sci(x: f64) -> String {
    let s = format!("{x:.9E}");
    let (mantissa, exp) = match s.split_once('E') {
        Some((m, e)) => (m, e),
        None => (s.as_str(), "0"),
    };
    let (sign, digits) = match exp.strip_prefix('-') {
        Some(d) => ('-', d),
        None => ('+', exp),
    };
    let pad = if digits.len() < 2 { "0" } else { "" };
    format!("{mantissa}E{sign}{pad}{digits}")
}

/// `Interface_FloatWriter::Convert`'s zero suppression (`cxx:148-174`): scan the
/// first 16 characters for the exponent marker, remember the exponent, cut the
/// mantissa there, drop its trailing zeros, and put the exponent back unless it
/// was exactly `+00`. Without a marker the cut point is the terminator, so a
/// fixed-format value only loses its trailing zeros (`150.000000000` -> `150.`).
fn float_zero_suppress(text: &str) -> String {
    let chars: Vec<char> = text.chars().collect();
    let mut cut = chars.len();
    let mut exp = String::new();
    let mut i = 0;
    while i < 16 && i < chars.len() {
        if chars[i] == 'e' || chars[i] == 'E' {
            let mut e = String::from("E");
            e.extend(chars[i + 1..(i + 5).min(chars.len())].iter());
            // `cxx:159-162`: only the exponent 0 loses its marker
            // (`lxp[1] == '+' && lxp[2] == '0' && lxp[3] == '0' && lxp[4] == '\0'`).
            if e != "E+00" {
                exp = e;
            }
            cut = i;
            break;
        }
        i += 1;
    }
    let mut mantissa: String = chars[..cut].iter().collect();
    while mantissa.ends_with('0') {
        mantissa.pop();
    }
    mantissa + &exp
}

/// Render an f64 in IGES parameter syntax (always carries a decimal point or
/// exponent marker).
fn num(x: f64) -> String {
    // `Sprintf` of a non-finite value yields `inf`/`nan`, which is not an IGES
    // number; no OCCT transfer produces one, and the port keeps the file
    // well-formed by writing the zero both formats would zero-suppress to.
    if !x.is_finite() {
        return "0.".into();
    }
    let in_range = (x >= FLOAT_RANGE_MIN && x < FLOAT_RANGE_MAX)
        || (x <= -FLOAT_RANGE_MIN && x > -FLOAT_RANGE_MAX);
    let text = if in_range { float_fixed(x) } else { float_sci(x) };
    float_zero_suppress(&text)
}

/// `nH` text parameter (`IGESData_GlobalSection::MakeHollerith`,
/// `IGESData_GlobalSection.cxx:47-69`). An empty string is an empty field.
fn hollerith(text: &str) -> String {
    if text.is_empty() {
        return String::new();
    }
    format!("{}H{text}", text.len())
}

/// `IGESData_GlobalSection::NewDateString(..., year < 0)` (`IGESData.cxx:247`):
/// `YYYYMMDD.HHMMSS`, always 15 characters.
fn iges_date_now() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0) as i64;
    let tod = secs.rem_euclid(86400);
    let days = secs.div_euclid(86400);
    let z = days + 719468;
    let era = (if z >= 0 { z } else { z - 146096 }) / 146097;
    let doe = (z - era * 146097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let mut y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = ((5 * doy + 2) / 153) as i64;
    let d = doy as i64 - (153 * mp + 2) / 5 + 1;
    let m = mp + if mp < 10 { 3 } else { -9 };
    if m <= 2 {
        y += 1;
    }
    let hh = tod / 3600;
    let mm = (tod % 3600) / 60;
    let ss = tod % 60;
    format!("{y:04}{m:02}{d:02}.{hh:02}{mm:02}{ss:02}")
}

/// Pointer placeholder for the `k`-th entry of an entity's reference list.
///
/// `#` never appears in IGES parameter data (strings use the `nH…` form), so a
/// `#k` token unambiguously marks a DE pointer position until
/// [`resolve_placeholders`] turns it into the final number.
fn ph(k: usize) -> String {
    format!("#{k}")
}

/// `n` comma-separated placeholders: `#0,#1,…`.
fn ph_list(n: usize) -> String {
    (0..n).map(ph).collect::<Vec<_>>().join(",")
}

/// True when every `#k` in `params` addresses an entry of the reference list.
fn pointers_all_known(params: &str, n_refs: usize) -> bool {
    placeholder_indices(params).all(|k| k < n_refs)
}

/// The `k` of every `#k` token, in text order.
fn placeholder_indices(params: &str) -> impl Iterator<Item = usize> + '_ {
    let bytes = params.as_bytes();
    let mut i = 0usize;
    let mut found: Vec<usize> = Vec::new();
    while i < bytes.len() {
        if bytes[i] == b'#' {
            let start = i + 1;
            let mut j = start;
            while j < bytes.len() && bytes[j].is_ascii_digit() {
                j += 1;
            }
            if j > start {
                if let Ok(k) = params[start..j].parse::<usize>() {
                    found.push(k);
                }
                i = j;
                continue;
            }
        }
        i += 1;
    }
    found.into_iter()
}

/// Replace every `#k` with the final **pointer** of the `k`-th referenced entity,
/// i.e. `Interface_InterfaceModel::DNum` = `2 * Number(entity) - 1`
/// (`Interface_InterfaceModel.cxx`), which is the sequence number of the
/// entity's **first Directory card**. `IGESData_IGESWriter::Send(handle, …)`
/// writes through `DNum` (`IGESData_IGESWriter.cxx:609-621`), so the D section
/// lays out two cards per entity and every in-model reference names an odd
/// D-section line. `refs` holds the referenced entities in parameter order and
/// `new_of` the reachability renumbering; a `#k` without a target keeps its
/// literal text (only reachable when an emitter records fewer refs than it
/// writes placeholders, which `emit_refs` already checks in debug builds).
fn resolve_placeholders(params: &str, refs: &[usize], new_of: &HashMap<usize, usize>) -> String {
    let bytes = params.as_bytes();
    let mut out = String::with_capacity(params.len());
    let mut i = 0usize;
    while i < bytes.len() {
        if bytes[i] == b'#' {
            let start = i + 1;
            let mut j = start;
            while j < bytes.len() && bytes[j].is_ascii_digit() {
                j += 1;
            }
            if j > start {
                if let Ok(k) = params[start..j].parse::<usize>() {
                    if let Some(n) = refs.get(k).and_then(|r| new_of.get(r)) {
                        out.push_str(&(2 * n - 1).to_string());
                        i = j;
                        continue;
                    }
                }
            }
        }
        out.push(bytes[i] as char);
        i += 1;
    }
    out
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

/// `Precision::IsNegativeInfinite` / `IsPositiveInfinite` clamp of a parameter
/// range before it is used as a curve parameter
/// (`GeomToIGES_GeomSurface.cxx:806-813` cone, `:1058-1065` tabulated
/// cylinder, `:1144-1150` revolution).
fn clamp_infinite_range(v0: f64, v1: f64) -> (f64, f64) {
    let a = if occt_core::precision::Precision::is_negative_infinite(v0) {
        -occt_core::precision::INFINITE
    } else {
        v0
    };
    let b = if occt_core::precision::Precision::is_positive_infinite(v1) {
        occt_core::precision::INFINITE
    } else {
        v1
    };
    (a, b)
}

/// The `gp_Ax3(gp_Ax2(...))` frames the elementary-surface transfers build their
/// generatrix circles on:
///
/// - sphere: `gp_Ax2(gp::Origin(), -gp::DY(), gp::DX())`
///   (`GeomToIGES_GeomSurface.cxx:893`);
/// - torus: `gp_Ax2(gp_Pnt(MajorRadius, 0., 0.), -gp::DY(), gp::DX())`
///   (`cxx:960`), i.e. the same frame shifted to the tube centre.
///
/// The circle is written in this frame (entity 100 + its own 124, see
/// [`IgesWriter::emit_circular_arc`]) and referenced by the 120.
fn revolution_generatrix_frame(location: GpPnt) -> GpAx3 {
    // `-gp::DY()` is `gp::DY()` reversed (`gp_Dir::Reverse`), so its zero
    // components carry a negative sign; building `(0, -1, 0)` from literals
    // loses them and the 124's signed zeros no longer match OCCT. The 3-arg
    // `gp_Ax2(N, Vx)` ctor (`gp_Ax2.hxx:73-80`) sets `vxdir = N ^ (Vx ^ N)`
    // (`gp_Dir::CrossCross`) and `vydir = N ^ vxdir`; `GpAx3::new` computes the
    // same `vydir` from the orthogonalised `vxdir`.
    let n = GpDir::from_axis(occt_core::gp::dir::DirAxis::Y).reversed();
    let vx = GpDir::from_axis(occt_core::gp::dir::DirAxis::X);
    let inner = vx.xyz().crossed(n.xyz());
    let xdir = GpDir::from_xyz(&n.xyz().crossed(&inner)).unwrap_or(vx);
    GpAx3::new(location, n, &xdir).expect("axis2 placement (Origin, -DY, DX)")
}

/// `ElCLib::CircleValue` on `frame` (`clib::circle_value`).
fn circle_point(frame: &GpAx3, radius: f64, u: f64) -> GpPnt {
    let p = frame
        .location()
        .coord
        .added(&frame.x_direction().xyz().multiplied(radius * u.cos()))
        .added(&frame.y_direction().xyz().multiplied(radius * u.sin()));
    GpPnt::from_xyz(&p)
}

/// Inverse of `IGESConvGeom_GeomBuilder::SetPosition(gp_Ax3)`
/// (`GeomBuilder.cxx:137-143`, `gp_Trsf::SetTransformation` then
/// `gp_Trsf::Invert` at `gp_Trsf.cxx:418-425`).
fn builder_local_trsf(frame: &GpAx3) -> GpTrsf {
    let mut ps = GpTrsf::identity();
    ps.set_transformation_from_to(frame, &xoy_frame());
    ps.scale = 1.0 / ps.scale;
    ps.matrix = ps.matrix.transpose();
    ps.loc = ps.matrix.multiplied(&ps.loc).multiply(-ps.scale);
    ps
}

fn eval_xyz(local: &GpTrsf, p: &GpPnt) -> (f64, f64, f64) {
    let q = p.transformed(local);
    (q.x(), q.y(), q.z())
}

/// `gp::XOY()` — the standard `gp_Ax3` frame (origin, `+Z` main, `+X`).
fn xoy_frame() -> GpAx3 {
    GpAx3::new(
        GpPnt::zero(),
        GpDir::from_axis(occt_core::gp::dir::DirAxis::Z),
        &GpDir::from_axis(occt_core::gp::dir::DirAxis::X),
    )
    .expect("gp::XOY()")
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
/// family, `Geom_Conic` → 104 (circle 100), `Geom_Line` → 110. A
/// `Geom_OffsetCurve` is not classified here: with the default
/// `write.iges.offset.mode = 0` it is approximated to a B-spline first
/// ([`IgesWriter::emit_offset_curve_default`]). The port has the equivalent
/// type queries on [`Curve`], so the previous
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

/// The curve `TransferCurve` approximates when `write.iges.offset.mode` is 0,
/// plus the range that approximation is then segmented to.
///
/// A trimmed curve is forwarded to its basis (`GeomToIGES_GeomCurve.cxx:467-475`).
/// `GeomTrimmedCurve` reports `FirstParameter() == 0` and `LastParameter() == 1`
/// while `untrimmed_basis` carries the basis interval, so `(a, b)` in that
/// unit interval is mapped back. `GeomTrimmedCurveBasis` already evaluates in
/// basis parameters (`first_parameter` is the trim start), and a bare
/// `Geom_OffsetCurve` is approximated as received.
fn offset_curve_request(
    curve: &dyn Curve,
    a: f64,
    b: f64,
) -> Option<(Arc<dyn Curve>, f64, f64)> {
    if curve.is_geom_trimmed() {
        let (basis, t0, t1) = curve.untrimmed_basis()?;
        if basis.offset_curve().is_none() {
            return None;
        }
        let unit = curve.first_parameter().abs() <= f64::EPSILON
            && (curve.last_parameter() - 1.0).abs() <= f64::EPSILON;
        let basis_is_unit = t0.abs() <= 1e-12 && (t1 - 1.0).abs() <= 1e-12;
        let (u0, u1) = if unit && !basis_is_unit {
            let span = t1 - t0;
            (t0 + a * span, t0 + b * span)
        } else {
            (a, b)
        };
        return Some((basis, u0, u1));
    }
    if curve.offset_curve().is_some() {
        return Some((Arc::from(curve.clone_dyn()), a, b));
    }
    None
}

/// `mycurve->Copy()` for the `Geom_BSplineCurve` transfer
/// (`GeomToIGES_GeomCurve.cxx:301`, `:337`): the port cannot downcast
/// `&dyn Curve`, so the copy is rebuilt from the same
/// `bspline_poles`/`bspline_weights`/`bspline_knots`/`nurbs_degree`/
/// `is_periodic` queries (equal data, no re-approximation). `None` where OCCT's
/// `occ::down_cast` would yield a null handle.
fn bspline_copy_of(curve: &dyn Curve) -> Option<occt_geom::bspline_curve::GeomBSplineCurve> {
    let poles = curve.bspline_poles()?.to_vec();
    let knots = curve.bspline_knots()?.to_vec();
    let degree = curve.nurbs_degree()?;
    Some(occt_geom::bspline_curve::GeomBSplineCurve {
        poles,
        weights: curve.bspline_weights().map(|w| w.to_vec()),
        knots,
        degree,
        periodic: curve.is_periodic(),
    })
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

/// `IsPlanar` (`GeomToIGES_GeomCurve.cxx:236-272`): the planar flag and plane
/// normal OCCT computes on the curve handed to `TransferCurve` (`cxx:414`,
/// `IsPlanar(start, Norm)`), i.e. before it is unperiodized or `Segment`-ed into
/// the copy whose poles are written out. The dispatch is OCCT's: `Geom_Line`
/// takes `GetAnyNormal` of its direction, the four conics their `Axis()`
/// direction, a `Geom_TrimmedCurve`/`Geom_OffsetCurve` recurses into the basis,
/// and `Geom_BSplineCurve`/`Geom_BezierCurve` fall through to
/// `ArePolesPlanar`. Anything else returns `false` with an untouched normal.
fn curve_planar_and_normal(curve: &dyn Curve) -> (bool, GpVec) {
    if let Some(l) = curve.gp_line() {
        return (true, any_normal(&l.direction().xyz()));
    }
    let conic = curve
        .gp_circ()
        .map(|c| c.pos)
        .or_else(|| curve.gp_ellipse().map(|e| e.pos))
        .or_else(|| curve.gp_hyperbola().map(|h| h.pos))
        .or_else(|| curve.gp_parabola().map(|p| p.pos));
    if let Some(pos) = conic {
        return (true, GpVec::from_xyz(&pos.direction().xyz()));
    }
    if let Some((base, _, _)) = curve.untrimmed_basis() {
        return curve_planar_and_normal(base.as_ref());
    }
    if let Some((base, _)) = curve.offset_curve() {
        return curve_planar_and_normal(base.as_ref());
    }
    if let Some(p) = curve.bspline_poles() {
        return poles_planar_and_normal(p);
    }
    if let Some(p) = curve.bezier_poles() {
        return poles_planar_and_normal(p);
    }
    (false, GpVec::new(0.0, 0.0, 0.0))
}

/// `GetAnyNormal` (`GeomToIGES_GeomCurve.cxx:168-189`): `(0,0,1)` when `|Z|` is
/// below `Precision::Confusion()`, otherwise the unit vector `(Z, 0, -X)`.
fn any_normal(v: &occt_core::gp::GpXyz) -> GpVec {
    if v.z().abs() < occt_core::precision::CONFUSION {
        return GpVec::new(0.0, 0.0, 1.0);
    }
    let n = occt_core::gp::GpXyz::new(v.z(), 0.0, -v.x());
    let m = n.modulus();
    if m < occt_core::precision::CONFUSION {
        GpVec::new(0.0, 0.0, 1.0)
    } else {
        GpVec::from_xyz(&n.divided(m))
    }
}

/// `Rational` (`Geom_BSplineCurve.cxx:98-108`): true when two consecutive
/// weights differ by more than `gp::Resolution()`.
fn weights_are_rational(weights: &[f64]) -> bool {
    weights.windows(2).any(|pair| {
        (pair[0] - pair[1]).abs() > occt_core::precision::REAL_SMALL
    })
}

/// `(U1 + U2) / 2` over a surface's bounds, the mid point
/// `GeomToIGES_GeomSurface.cxx:1212-1213` evaluates the offset indicator at.
///
/// `Geom_Surface::Bounds` reports an unbounded direction as
/// `Precision::Infinite()` = `RealLast()`, so OCCT's mid point of a
/// doubly-infinite range is 0 and of a half-infinite range is `RealLast()/2`;
/// the port stores `f64::INFINITY` for those, and `INFINITY + -INFINITY` would be
/// NaN - a NaN parameter poisons the returned normal - so a non-finite range
/// collapses to 0 here. Elementary surfaces have constant derivatives in the
/// unbounded direction, so the indicator is unaffected either way.
fn bounds_mid_point(a: f64, b: f64) -> f64 {
    if a.is_finite() && b.is_finite() {
        (a + b) / 2.0
    } else {
        0.0
    }
}

/// `CSLib::Normal(D1U, D1V, theSinTol, theStatus, theNormal)`
/// (`CSLib.cxx:45-80`), the overload reached from `GeomLProp_SLProps::Normal()`
/// through `LProp_SurfaceUtils::ComputeSurfNormal`
/// (`GeomLProp_SurfaceUtils.hxx:348-358`).
///
/// Returns the unit normal when the status is `CSLib_Done`. Every other status
/// (`D1IsNull`, `D1uIsNull`, `D1vIsNull`, `D1uIsParallelD1v`) makes OCCT's
/// `Normal()` throw `LProp_NotDefined`, and is reported here as `None`.
fn surface_normal_from_d1(d1u: &GpVec, d1v: &GpVec, sin_tol: f64) -> Option<GpVec> {
    // `gp::Resolution()` is `RealSmall()` (`gp.hxx:60`), not `Precision`.
    const GP_RESOLUTION: f64 = occt_core::precision::REAL_SMALL;
    let mag_u = d1u.square_magnitude();
    let mag_v = d1v.square_magnitude();
    if mag_u <= GP_RESOLUTION || mag_v <= GP_RESOLUTION {
        return None;
    }
    let cross = d1u.crossed(d1v);
    let sin2 = cross.square_magnitude() / (mag_u * mag_v);
    if sin2 < sin_tol * sin_tol {
        return None;
    }
    // `gp_Dir(theD1UxD1V)`: OCCT normalises unconditionally here (the magnitude
    // is above `gp::Resolution()` because `sin2` is strictly positive), so this
    // divides explicitly instead of going through `GpVec::normalized`, which
    // keeps an unnormalised vector below `Precision::Resolution()`.
    let m = cross.magnitude();
    if m <= 0.0 {
        return None;
    }
    Some(cross.divided(m))
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
// 2-D (UV) curve transfer
//
// 'BRepToIGES_BRWire::TransferEdge(edge, face, originMap, length, false)'
// ('BRepToIGES_BRWire.cxx:340-585') writes the edge's p-curve as an IGES curve
// entity: the p-curve is corrected per surface type, transformed by
// 'ShapeBuild_Edge::TransformPCurve', reversed for a REVERSED edge, and handed
// to 'Geom2dToIGES_Geom2dCurve::Transfer2dCurve', which lifts it into the plane
// Z=0 and reuses the 3-D 'GeomToIGES_GeomCurve::TransferCurve' dispatch.
//
// The port has no 'write.surfacecurve.mode' switch, so 'GetPCurveMode()' keeps
// OCCT's default 1 (On, 'Interface_StaticStandards.cxx:85'); every call here
// passes 'theIsBRepMode = false' ('BRepToIGES_BRShell.cxx:280', ':312', ':346'),
// which makes 'analyticMode = (GetConvertSurfaceMode() == 0 && false)' false
// ('BRepToIGES_BRWire.cxx:356').
// ---------------------------------------------------------------------------

/// 'gp_Dir2d::Transform' ('gp_Dir2d.cxx:80-105'). The port's
/// 'GpTrsf2d::transforms_xy' handles points; directions follow OCCT's own rule
/// (raw matrix, normalise, reverse for a negative scale factor).
fn dir2d_transform(t: &GpTrsf2d, d: &GpDir2d) -> GpDir2d {
    match t.form() {
        TrsfForm::Identity | TrsfForm::Translation => *d,
        TrsfForm::PntMirror => GpDir2d::new(-d.x, -d.y).unwrap_or(*d),
        TrsfForm::Scale => {
            if t.scale_factor() < 0.0 {
                GpDir2d::new(-d.x, -d.y).unwrap_or(*d)
            } else {
                *d
            }
        }
        _ => {
            let mut xy = GpXY::new(d.x, d.y);
            xy.multiply_mat2d(t.h_vectorial_part());
            let _ = xy.normalize();
            if t.scale_factor() < 0.0 {
                xy.reverse();
            }
            GpDir2d::new(xy.x, xy.y).unwrap_or(*d)
        }
    }
}

/// 'gp_Pnt2d::Transform' ('gp_Pnt2d.hxx'), which is exactly what
/// 'GpTrsf2d::transforms_xy' applies.
fn pnt2d_transform(t: &GpTrsf2d, p: &GpPnt2d) -> GpPnt2d {
    let mut c = p.coord;
    t.transforms_xy(&mut c);
    GpPnt2d::from_xy(c)
}

/// 'gp_Ax22d::Transform' ('gp_Ax22d.hxx:360-367'): location, X and Y directions
/// are transformed independently (a mirroring transform may leave the axis
/// indirect, exactly as OCCT's flag-free 'gp_Ax22d' does).
fn ax22d_transform(t: &GpTrsf2d, a: &GpAx22d) -> GpAx22d {
    GpAx22d {
        point: pnt2d_transform(t, &a.point),
        vxdir: dir2d_transform(t, &a.vxdir),
        vydir: dir2d_transform(t, &a.vydir),
    }
}

/// 'gp_Trsf2d::SetMirror(const gp_Ax2d&)' ('gp_Trsf2d.cxx:31-46').
///
/// The port's 'GpTrsf2d::set_mirror_ax2d' stores the *positive* reflection
/// matrix together with 'scale = -1', so 'transforms_xy' applies '-R'; this
/// builds OCCT's own matrix (the negative reflection, 'scale = -1') so the
/// composite is the reflection R.
fn mirror_ax2d_trsf(a: &GpAx2d) -> GpTrsf2d {
    let (vx, vy) = (a.direction().x, a.direction().y);
    let (x0, y0) = (a.location().x(), a.location().y());
    let matrix = GpMat2d::new(
        1.0 - 2.0 * vx * vx,
        -2.0 * vx * vy,
        -2.0 * vx * vy,
        1.0 - 2.0 * vy * vy,
    );
    let loc = GpXY::new(
        -2.0 * ((vx * vx - 1.0) * x0 + vx * vy * y0),
        -2.0 * (vx * vy * x0 + (vy * vy - 1.0) * y0),
    );
    GpTrsf2d {
        scale: -1.0,
        shape: TrsfForm::Ax1Mirror,
        matrix,
        loc,
    }
}

/// 'Mirror(gp_Ax2d(gp::Origin2d(), gp::Dir2d(1., 1.)))'
/// ('BRepToIGES_BRWire.cxx:463', ':471').
fn mirror_origin_dir11() -> GpTrsf2d {
    let ax = GpAx2d::new(
        GpPnt2d::zero(),
        GpDir2d::new(1.0, 1.0).expect("gp::Dir2d(1,1) is non-null"),
    );
    mirror_ax2d_trsf(&ax)
}

/// 'Mirror(gp::OX2d())' ('BRepToIGES_BRWire.cxx:464', ':472').
fn mirror_ox2d() -> GpTrsf2d {
    let ax = GpAx2d::new(
        GpPnt2d::zero(),
        GpDir2d::new(1.0, 0.0).expect("gp::DX2d is non-null"),
    );
    mirror_ax2d_trsf(&ax)
}

/// 'Curve2d->Translate(gp_Vec2d(dx, dy))'.
fn translate2d(dx: f64, dy: f64) -> GpTrsf2d {
    let mut t = GpTrsf2d::identity();
    t.set_translation_vec(&GpVec2d::new(dx, dy));
    t
}

/// 'Geom2d_Curve::TransformedParameter' ('Geom2d_Curve.cxx:41-44' default;
/// 'Geom2d_Line.cxx:246-253', 'Geom2d_Parabola.cxx:249-256',
/// 'Geom2d_TrimmedCurve.cxx:297-300', 'Geom2d_OffsetCurve.cxx:423-426'): only a
/// line, a parabola and the trimmed/offset wrappers over them rescale the
/// parameter; every other type keeps it.
fn transformed_parameter(curve: &dyn Curve2d, u: f64, t: &GpTrsf2d) -> f64 {
    if curve.gp_lin2d().is_some() || curve.gp_parab2d().is_some() {
        if occt_core::precision::Precision::is_infinite(u) {
            u
        } else {
            u * t.scale_factor().abs()
        }
    } else if let Some(b) = curve.trimmed_basis() {
        transformed_parameter(b, u, t)
    } else if let Some(b) = curve.offset_basis() {
        transformed_parameter(b, u, t)
    } else {
        u
    }
}

/// 'Geom2d_*::Transform' for the concrete types this port stores as p-curves -
/// 'Geom2d_Line', 'Geom2d_Circle', 'Geom2d_Ellipse', 'Geom2d_Hyperbola',
/// 'Geom2d_Parabola', 'Geom2d_BSplineCurve', 'Geom2d_TrimmedCurve'.
///
/// The port's 'Curve2d::transform' only moves the *location* of a line / circle /
/// conic ('GpLin2d::transform' etc. leave the axis directions untouched), so the
/// transform is rebuilt from the exact-type queries instead of delegating to it.
/// Any other curve (Bezier / Offset) keeps the port's own 'transform'.
fn curve2d_transformed(curve: &Arc<dyn Curve2d>, t: &GpTrsf2d) -> Arc<dyn Curve2d> {
    // 'Geom2d_TrimmedCurve::Transform' ('Geom2d_TrimmedCurve.cxx:287-293'):
    // transform the basis, then set the trim to the transformed parameters.
    if let Some(b) = curve.trimmed_basis() {
        let nb = curve2d_transformed(&Arc::from(b.clone_dyn()), t);
        let u1 = transformed_parameter(b, curve.first_parameter(), t);
        let u2 = transformed_parameter(b, curve.last_parameter(), t);
        return Arc::new(Geom2dTrimmedCurve::new_sense(nb, u1, u2, true, false));
    }
    if let Some(l) = curve.gp_lin2d() {
        // 'Geom2d_Line::Transform' ('Geom2d_Line.cxx:239-242') ->
        // 'gp_Lin2d::Position().Transform' (location and direction).
        return Arc::new(Geom2dLine::new(GpAx2d::new(
            pnt2d_transform(t, &l.location()),
            dir2d_transform(t, l.direction()),
        )));
    }
    if let Some(c) = curve.gp_circ2d() {
        // 'gp_Circ2d::Transform': axis transformed, 'radius *= |ScaleFactor|'.
        return Arc::new(Geom2dCircle::new(GpCirc2d::new(
            ax22d_transform(t, c.position()),
            c.radius() * t.scale_factor().abs(),
        )));
    }
    if let Some(e) = curve.gp_elips2d() {
        // 'gp_Elips2d::Transform': axis transformed, both radii scaled.
        let s = t.scale_factor().abs();
        return Arc::new(Geom2dEllipse::new(GpElips2d::new(
            ax22d_transform(t, e.axis()),
            e.major_radius * s,
            e.minor_radius * s,
        )));
    }
    if let Some(h) = curve.gp_hypr2d() {
        let s = t.scale_factor().abs();
        return Arc::new(Geom2dHyperbola::new(GpHypr2d::new(
            ax22d_transform(t, h.axis()),
            h.major_radius * s,
            h.minor_radius * s,
        )));
    }
    if let Some(p) = curve.gp_parab2d() {
        return Arc::new(Geom2dParabola::new(GpParab2d::new(
            ax22d_transform(t, p.axis()),
            p.focal * t.scale_factor().abs(),
        )));
    }
    if let (Some((xs, ys)), Some(knots), Some(degree)) = (
        curve.bspline_poles2d(),
        curve.bspline_knots2d(),
        curve.bspline_degree(),
    ) {
        // 'Geom2d_BSplineCurve::Transform': the poles carry the affine map. The
        // port's 2-D B-spline is non-rational (no weights).
        let poles: Vec<GpPnt2d> = xs
            .iter()
            .zip(ys.iter())
            .map(|(x, y)| pnt2d_transform(t, &GpPnt2d::new(*x, *y)))
            .collect();
        let (out_xs, out_ys): (Vec<f64>, Vec<f64>) =
            poles.iter().map(|p| (p.x(), p.y())).unzip();
        return Arc::new(Geom2dBSplineCurve {
            xs: out_xs,
            ys: out_ys,
            weights: None,
            knots: knots.to_vec(),
            degree,
            periodic: curve.is_periodic(),
        });
    }
    Arc::from(curve.transformed(t))
}

/// 'ShapeBuild_Edge::TransformPCurve(pcurve, trans, uFact, aFirst, aLast)'
/// ('ShapeBuild_Edge.cxx:596-699').
///
/// Returns 'None' for the arms this port cannot reproduce - the caller then
/// writes no UV curve for that edge (the null handle of the OCCT transfer).
fn transform_pcurve(
    pcurve: &Arc<dyn Curve2d>,
    trans: &GpTrsf2d,
    u_fact: f64,
    first: &mut f64,
    last: &mut f64,
) -> Option<Arc<dyn Curve2d>> {
    // 'cxx:602-608': a non-identity trans is applied and the range ends are
    // mapped through 'TransformedParameter'.
    let mut result: Arc<dyn Curve2d> = if trans.form() != TrsfForm::Identity {
        let r = curve2d_transformed(pcurve, trans);
        *first = transformed_parameter(r.as_ref(), *first, trans);
        *last = transformed_parameter(r.as_ref(), *last, trans);
        r
    } else {
        Arc::from(pcurve.clone_dyn())
    };
    if u_fact == 1.0 {
        return Some(result); // 'cxx:609-612'
    }
    // 'cxx:614-618': a trimmed curve is replaced by its basis.
    if let Some(b) = result.trimmed_basis() {
        result = Arc::from(b.clone_dyn());
    }
    // 'tMatu.SetAffinity(gp::OY2d(), uFact)' ('cxx:620-621') scales the X (U)
    // coordinate by uFact ('gp_GTrsf2d::SetAffinity', 'gp_GTrsf2d.cxx:24-38',
    // with the Y axis: 'matrix = (uFact, 0; 0, 1)').
    if result.gp_lin2d().is_some() {
        // 'cxx:624-641': scale the two range points, rebuild the line through
        // them and take the range from 'ElCLib::Parameter'.
        let pf = result.d0(*first);
        let pl = result.d0(*last);
        let pf = GpPnt2d::new(pf.x() * u_fact, pf.y());
        let pl = GpPnt2d::new(pl.x() * u_fact, pl.y());
        let v = GpVec2d::new(pl.x() - pf.x(), pl.y() - pf.y());
        let dir = GpDir2d::from_vec2d(&v).ok()?;
        let line = GpLin2d::from_pnt_dir(pf, dir);
        *first = occt_core::elib::clib2d::parameter_lin2d(&line, &pf);
        *last = occt_core::elib::clib2d::parameter_lin2d(&line, &pl);
        return Some(Arc::new(Geom2dLine::new(GpAx2d::new(pf, dir))));
    }
    // 'cxx:642-656': a 2-D Bezier scales every pole on X and keeps its own
    // range; the weights (rational Bezier) are carried over unchanged.
    if result.is_bezier2d() {
        let poles: Vec<GpPnt2d> = result
            .poles2d()?
            .iter()
            .map(|p| GpPnt2d::new(p.x() * u_fact, p.y()))
            .collect();
        let bez = match result.bezier_weights2d() {
            Some(w) => Geom2dBezierCurve::rational(poles, w.to_vec()),
            None => Geom2dBezierCurve::new(poles),
        }
        .ok()?;
        return Some(Arc::new(bez));
    }
    let is_conic = result.gp_circ2d().is_some()
        || result.gp_elips2d().is_some()
        || result.gp_hypr2d().is_some()
        || result.gp_parab2d().is_some();
    // 'cxx:658-698': the conic arm trims the curve to 'aFirst'/'aLast' and
    // re-runs it through 'Geom2dConvert_ApproxCurve' (falling back to
    // 'Geom2dConvert::CurveToBSplineCurve(thecurve, Convert_QuasiAngular)');
    // any remaining 2-D type goes through 'CurveToBSplineCurve(result,
    // Convert_QuasiAngular)'; a B-spline is used as it is. All three then have
    // their poles scaled on X.
    let mut bs = if is_conic {
        let tcurve: Arc<dyn Curve2d> =
            Arc::new(Geom2dTrimmedCurve::new(result.clone(), *first, *last));
        let approx = occt_geom2d::Geom2dConvertApproxCurve::new(
            tcurve.as_ref(),
            occt_core::precision::APPROXIMATION,
            occt_core::kernel::geomabs::Shape::C1,
            100,
            6,
        );
        let bs = match approx.curve() {
            Some(curve) => curve.clone(),
            None => occt_geom2d::geom2d_convert::curve_to_bspline_curve_bspl(
                tcurve.as_ref(),
                occt_core::convert::ParameterisationType::QuasiAngular,
            )?,
        };
        // 'cxx:673-674': 'aFirst'/'aLast' follow the approximation's range.
        *first = bs.first_parameter();
        *last = bs.last_parameter();
        bs
    } else if !result.is_bspline2d() {
        occt_geom2d::geom2d_convert::curve_to_bspline_curve_bspl(
            result.as_ref(),
            occt_core::convert::ParameterisationType::QuasiAngular,
        )?
    } else {
        result.bspline_copy2d()?
    };
    // 'cxx:688-698': transform the poles. 'aFirst'/'aLast' are only rewritten in
    // the conic arm above, so a plain B-spline keeps its range.
    bs.xs = bs.xs.iter().map(|x| x * u_fact).collect();
    Some(Arc::new(bs))
}

/// 'Adaptor3d_CurveOnSurface.cxx:72-78' (to3d of a 'gp_Ax22d' on
/// 'gp_Pln(0,0,1,0)'): 'gp_Ax2(P, VX.Crossed(VY), VX)' with 'P = (x, y, 0)'.
fn promote_ax22d(a: &GpAx22d) -> Option<GpAx2> {
    let p = GpPnt::new(a.point.x(), a.point.y(), 0.0);
    let vx = GpDir::from_xyz(&GpXyz::new(a.vxdir.x, a.vxdir.y, 0.0)).ok()?;
    let vy = GpDir::from_xyz(&GpXyz::new(a.vydir.x, a.vydir.y, 0.0)).ok()?;
    let n = vx.crossed(&vy).ok()?;
    Some(GpAx2::new(p, n, vx).unwrap_or_else(|_| GpAx2::from_axis(p, n)))
}

/// 'Geom2dToIGES_Geom2dCurve::Transfer2dCurve' ('Geom2dToIGES_Geom2dCurve.cxx:52-68'):
/// 'GC.TransferCurve(GeomAPI::To3d(start, gp_Pln(0,0,1,0)), Udeb, Ufin)'.
///
/// 'GeomAPI::To3d' ('GeomAPI.cxx:56-65') puts the 2-D curve in an
/// 'Adaptor3d_CurveOnSurface' over the plane Z=0 and 'GeomAdaptor::MakeCurve'
/// ('GeomAdaptor.cxx:45-92') builds the concrete 'Geom_*' curve of the same
/// type. 'Geom2dAdaptor_Curve::load' unwraps a 'Geom2d_TrimmedCurve' and keeps
/// the basis, and 'GeomToIGES_GeomCurve::TransferCurve(Geom_TrimmedCurve)'
/// ('GeomToIGES_GeomCurve.cxx:456-477') transfers that basis with the requested
/// range - so promoting the basis with '[First, Last]' is equivalent.
fn promote_curve2d(curve: &dyn Curve2d) -> Option<Box<dyn Curve>> {
    let basis: &dyn Curve2d = curve.trimmed_basis().unwrap_or(curve);
    if let Some(l) = basis.gp_lin2d() {
        let p = GpPnt::new(l.location().x(), l.location().y(), 0.0);
        let d = GpDir::from_xyz(&GpXyz::new(l.direction().x, l.direction().y, 0.0)).ok()?;
        return Some(Box::new(GeomLine::new(GpLin::from_pnt_dir(p, d))));
    }
    if let Some(c) = basis.gp_circ2d() {
        return Some(Box::new(GeomCircle::new(GpCirc::new(
            promote_ax22d(c.position())?,
            c.radius(),
        ))));
    }
    if let Some(e) = basis.gp_elips2d() {
        return Some(Box::new(GeomEllipse::new(GpElips::new(
            promote_ax22d(e.axis())?,
            e.major_radius,
            e.minor_radius,
        ))));
    }
    if let Some(h) = basis.gp_hypr2d() {
        return Some(Box::new(GeomHyperbola::new(GpHypr::new(
            promote_ax22d(h.axis())?,
            h.major_radius,
            h.minor_radius,
        ))));
    }
    if let Some(p) = basis.gp_parab2d() {
        return Some(Box::new(GeomParabola::new(GpParab::new(
            promote_ax22d(p.axis())?,
            p.focal,
        ))));
    }
    if let (Some((xs, ys)), Some(knots), Some(degree)) = (
        basis.bspline_poles2d(),
        basis.bspline_knots2d(),
        basis.bspline_degree(),
    ) {
        // `Adaptor3d_CurveOnSurface::BSpline()` (`cxx:1489-1520`) lifts every
        // pole with `to3d` (`cxx:57-60`): `ElSLib::Value` on the plane, not
        // `(x, y, 0)`. `GeomAPI::To3d` uses `gp_Pln(0, 0, 1, 0)`
        // (`Geom2dToIGES_Geom2dCurve.cxx:52-68`). The port's 2-D B-spline is
        // non-rational, so the `IsRational()` arm (weights array) has no
        // ported source.
        let plane = occt_core::gp::GpPln::from_coefficients(0.0, 0.0, 1.0, 0.0).ok()?;
        let poles: Vec<GpPnt> = xs
            .iter()
            .zip(ys.iter())
            .map(|(x, y)| occt_core::elib::slib::plane_value(&plane, *x, *y))
            .collect();
        return Some(Box::new(GeomBSplineCurve {
            poles,
            weights: None,
            knots: knots.to_vec(),
            degree,
            periodic: basis.is_periodic(),
        }));
    }
    // `Adaptor3d_CurveOnSurface::Bezier()` (`cxx:1458-1485`) lifts the 2-D
    // Bezier's poles to Z=0 and `GeomAdaptor::MakeCurve` (`GeomAdaptor.cxx:60-77`)
    // builds the matching `Geom_BezierCurve`.
    //
    // UNPORTED (`cxx:1474-1481`): the rational sub-arm
    // (`Bez2d->IsRational()` -> `Geom_BezierCurve(Poles, Weights)`) has no port
    // equivalent - `occt_geom::GeomBezierCurve` is non-rational, the same gap
    // noted in `occt_geom::convert_bspl::curve_to_bspline_curve` (`cxx:313-321`).
    if basis.is_bezier2d() {
        let plane = occt_core::gp::GpPln::from_coefficients(0.0, 0.0, 1.0, 0.0).ok()?;
        let poles: Vec<GpPnt> = basis
            .poles2d()?
            .iter()
            .map(|p| occt_core::elib::slib::plane_value(&plane, p.x(), p.y()))
            .collect();
        return Some(Box::new(GeomBezierCurve::new(poles).ok()?));
    }
    // UNPORTED: any other type has no 'GeomAdaptor::MakeCurve' arm either
    // ('GeomAdaptor.cxx:79-80' throws OtherCurve).
    None
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
    /// The other entities this one references, in record order - `OwnShared` in
    /// OCCT terms (`IGESData_GeneralModule::OwnShared` feeds the writer's pointer
    /// resolution). Recorded so the pointers can be checked (and, for T-85, later
    /// remapped) without parsing the parameter text back into numbers.
    refs: Vec<usize>,
    /// Parameter data, always ending with the record delimiter `;`.
    params: String,
    /// The entity's status number, DE card 1 columns 65-72: BlankStatus,
    /// SubordinateStatus, UseFlag and HierarchyStatus, two columns each.
    /// Filled by [`compute_status`] from the finished model;
    /// `IGESData_IGESEntity` starts at `InitStatus(0, 0, 0, 0)`
    /// (`IGESData_IGESEntity.cxx:53-63`).
    status: [i32; 4],
}

impl Ent {
    /// Number of Parameter cards this entity occupies, which card 2 of its
    /// directory entry reports (`IGESData_IGESWriter.cxx:835`:
    /// `v[15] = thepnum.Value(i + 1) - thepnum.Value(i)`).
    fn param_line_count(&self) -> usize {
        self.param_chunks().len()
    }

    /// The entity's Parameter cards, as `IGESData_IGESWriter` builds them.
    ///
    /// `OwnParams` clears the card buffer and writes the type number, without a
    /// separator (`IGESData_IGESWriter.cxx:418-434`); every own parameter is
    /// then written by a `Send`, i.e. the separator followed by the value
    /// (`cxx:560-565` for integers, `:586-589` for reals, `:610-620` for
    /// pointers — a *void* parameter writes the separator alone, `:557-560`);
    /// `EndEntity` closes the record with `theendm` and dumps the last card
    /// (`:471-483`). Splitting the parameter text on the separator and feeding
    /// the pieces to [`CardBuffer`] is that same sequence, so the cards break
    /// where OCCT's break instead of every 64 columns.
    fn param_chunks(&self) -> Vec<String> {
        let body = self.params.strip_suffix(';').unwrap_or(self.params.as_str());
        let mut buf = CardBuffer::new(MAXCARS_P);
        for (i, field) in body.split(',').enumerate() {
            if i > 0 {
                buf.add_char(b',');
            }
            buf.add_string(field);
        }
        buf.add_char(b';');
        if !buf.line.is_empty() {
            buf.flush();
        }
        buf.cards
    }

    /// The two Directory Entry cards of the entity
    /// (`IGESData_IGESWriter.cxx:276-365` filling `v[0..16]`, `:809-880` laying
    /// the cards out): eight 8-column fields plus four 2-column status fields on
    /// the first card and nine 8-column fields on the second - 72 data columns
    /// each. The four status fields are printed `%2.2d`
    /// (`IGESData_IGESWriter.cxx:836-840`), i.e. zero filled, not blank filled:
    /// a status of 0 is `00`, so card 1 of a default entity ends on `00000000`.
    /// Card 1 field 2 is the sequence number of the entity's **first
    /// parameter card** (`v[1] = thepnum.Value(i)`, `cxx:834`), field 7 the
    /// transformation matrix pointer (`cxx:324-331`); card 2 field 4 is the number
    /// of parameter cards (`v[15]`, `cxx:835`), field 5 the form number
    /// (`v[16] = anent->FormNumber()`, `cxx:364`).
    ///
    /// The status number (`v[8..11]`) is [`Ent::status`], computed on the
    /// finished model by [`compute_status`].
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
            "{ty}{pstart}{}{}{}{}{trsf}{}{:02}{:02}{:02}{:02}",
            field8("0"),
            field8("0"),
            field8("0"),
            field8("0"),
            field8("0"),
            self.status[0],
            self.status[1],
            self.status[2],
            self.status[3]
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

/// `IGESData_BasicEditor::ComputeStatus` (`IGESData_BasicEditor.cxx:216-328`)
/// followed by `AutoCorrectModel` -> `AutoCorrect` (`:330-433`), which
/// `IGESControl_Writer::ComputeModel` runs before printing
/// (`IGESControl_Writer.cxx:256-262`). Both work on the numbered model, i.e.
/// after [`IgesWriter::final_entities`], and fill DE card 1's status number.
///
/// ComputeStatus has two phases, documented as such at `cxx:231-255`:
///
/// * **Subordinate status** (`cxx:257-281`): every entity ORs `1` into the
///   status of each entity its *own* parameters reference
///   (`gmod->OwnSharedCase` = [`Ent::refs`]), and `2` instead of `1` when the
///   referring entity is an associativity (type 402 or 404). DE-part
///   references are not `OwnShared`, so the transformation matrix
///   ([`Ent::trsf`], card 1 field 7) contributes nothing - the oracle confirms
///   it: all 347 type-124 entries are `00000000`, while every entity
///   referenced from a parameter record is `00010000` or `00020000`.
/// * **UseFlag** (`cxx:283-309`): not ported. It propagates only for types
///   200-299 (`G.GetFromEntity(ent, true, 1)`), 134, 116 and 132
///   (`G.GetFromEntity(ent, true, 4)`); the port's transfer classes create
///   none of those, and the phase reads an `Interface_Graph` status array
///   whose shared sets include the DE-part pointers this writer does not
///   model. `G.Status(i)` is therefore always 0, so the second phase
///   (`cxx:315-327`) keeps each entity's own UseFlag when it is non-zero and
///   takes 0 otherwise.
///
/// `AutoCorrect` then applies the general module's `DirChecker(CN, ent)`:
/// `Correct` forces a status field whenever it differs from the
/// `...Required(val)` the tool declared (`IGESData_DirChecker.cxx:422-451`),
/// while a field left at its `-100` "do not test" default
/// (`IGESData_DirChecker.cxx:33-36`) is untouched (`thegraphier == -100` also
/// disables the DE graphics cleanup, `:400-420`). Of the types this writer
/// emits (100/102/108/110/120/124/126/128/140/142/144/402) exactly two declare
/// a status:
///
/// * 142 - `UseFlagRequired(5)` (`IGESGeom_ToolCurveOnSurface.cxx:200-211`),
///   the source of the oracle's `00010500` on every type 142;
/// * 144 - `UseFlagRequired(0)` (`IGESGeom_ToolTrimmedSurface.cxx:264-276`),
///   which is a no-op for an entity whose UseFlag is already 0.
///
/// `AutoCorrect`'s specific-module call `smod->OwnCorrect(CN, ent)`
/// (`cxx:413-418`) is deliberately *not* ported: in OCCT 8.0.0 the 142
/// correction is registered under case 9 but casts to `IGESGeom_Boundary`
/// (`IGESGeom_SpecificModule.cxx:337-345`), so `anent.IsNull()` and the
/// CurveUV/Curve3D of a 142 keep UseFlag 0 - which is what the oracle shows.
fn compute_status(entities: &mut [Ent]) {
    let nb = entities.len();
    // NCollection_Array1<int> subs(0, nb); subs.Init(0) (`cxx:226-228`).
    let mut subs = vec![0i32; nb + 1];
    for e in entities.iter() {
        let bit = if e.ty == 402 || e.ty == 404 { 2 } else { 1 };
        for r in &e.refs {
            // themodel->Number(sh.Value()) is 0 for anything outside the model;
            // index 0 of `subs` absorbs those (`cxx:269-280`).
            if let Some(s) = subs.get_mut(*r) {
                *s |= bit;
            }
        }
    }

    for (i, e) in entities.iter_mut().enumerate() {
        // InitStatus(bl, subs.Value(i), uf, hy), `uf == 0` taking G.Status(i)
        // (`cxx:315-327`). Blank and hierarchy status stay as the entity set
        // them - the port never changes them from 0.
        e.status = [e.status[0], subs[i + 1], e.status[2], e.status[3]];

        // DirChecker::Correct (`IGESData_DirChecker.cxx:422-451`).
        let required = match e.ty {
            142 => Some(5), // IGESGeom_ToolCurveOnSurface.cxx:209
            144 => Some(0), // IGESGeom_ToolTrimmedSurface.cxx:272
            _ => None,
        };
        if let Some(val) = required {
            if val != e.status[2] {
                e.status[2] = val;
            }
        }
    }
}

/// Incremental IGES writer holding the entity list.
struct IgesWriter {
    entities: Vec<Ent>,
    /// Top-level entity of every shape handed to [`IgesWriter::emit_shape`], in
    /// that order. `IGESControl_Writer::AddEntity` feeds each one to
    /// `Interface_InterfaceModel::AddWithRefs` (`IGESControl_Writer.cxx:243-252`),
    /// so these are the DFS roots of the written model (T-85 step 2).
    roots: Vec<usize>,
    point_entities: HashMap<usize, usize>,
    /// Largest absolute model coordinate (`IGESControl_Writer.cxx:180-189`,
    /// `MaxMaxCoords` of the shape box). `0` leaves global field 20 empty,
    /// which is `hasMaxCoord == false` (`IGESData_GlobalSection.cxx:519-527`).
    max_coord: f64,
}

impl IgesWriter {
    fn new() -> Self {
        Self {
            entities: Vec::new(),
            roots: Vec::new(),
            point_entities: HashMap::new(),
            max_coord: 0.0,
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
            refs: Vec::new(),
            params,
            status: [0; 4],
        });
        self.entities.len()
    }

    /// Register one entity that references others, recording those pointers in
    /// record order (T-85's prerequisite for OCCT's reachability-based writing).
    ///
    /// The parameter text carries a `#k` placeholder for the `k`-th entry of
    /// `refs` instead of a DE number; [`IgesWriter::final_entities`] substitutes
    /// the final number after the reachability renumbering. OCCT can write the
    /// numbers directly because its model is numbered while it is built
    /// (`Interface_InterfaceModel::AddWithRefs`, `Interface_InterfaceModel.cxx:652-692`);
    /// the port numbers at write time, so the reference *index* is what the
    /// parameter text can hold.
    fn emit_refs(&mut self, ty: i32, form: i32, params: String, refs: &[usize]) -> usize {
        debug_assert!(
            pointers_all_known(&params, refs.len()),
            "entity type {ty}: parameter text has placeholders outside 0..{}",
            refs.len()
        );
        let de = self.emit(ty, form, params);
        self.entities[de - 1].refs = refs.to_vec();
        de
    }

    /// Every recorded pointer must address an entity that exists and precede the
    /// referencing entity (`IGESData_IGESWriter::Send` resolves pointers through
    /// the model, so a dangling number would be written verbatim). Debug-only:
    /// the emitters below build the pointers themselves.
    fn check_refs(&self) {
        for (i, e) in self.entities.iter().enumerate() {
            for r in &e.refs {
                debug_assert!(
                    *r >= 1 && *r <= self.entities.len() && *r != i + 1,
                    "entity {} (type {}) references {} - out of range or self",
                    i + 1,
                    e.ty,
                    r
                );
            }
        }
    }

    /// The entities OCCT's `AddWithRefs` would put in the model, in the order it
    /// adds them, renumbered and with the `#k` placeholders resolved.
    ///
    /// `Interface_InterfaceModel::AddWithRefs` (`Interface_InterfaceModel.cxx:652-692`)
    /// adds an entity and then recurses into its shared entities **in order**,
    /// skipping anything already added: a pre-order DFS from the roots
    /// (`IGESControl_Writer::AddEntity`, `IGESControl_Writer.cxx:243-252`). The
    /// shared entities are `IGESData_GeneralModule::FillSharedCase`'s list —
    /// the DE-part entities first (only field 7, the transformation matrix, is
    /// ever set by this writer) and then the own-parameter references, which is
    /// the record order of [`Ent::refs`].
    fn final_entities(&self) -> Vec<Ent> {
        let mut order: Vec<usize> = Vec::new();
        let mut seen: HashSet<usize> = HashSet::new();
        let mut stack: Vec<usize> = self.roots.iter().rev().copied().collect();
        while let Some(old) = stack.pop() {
            if old < 1 || old > self.entities.len() || !seen.insert(old) {
                continue;
            }
            order.push(old);
            let e = &self.entities[old - 1];
            let mut kids: Vec<usize> = Vec::new();
            if let Some(t) = e.trsf {
                kids.push(t);
            }
            kids.extend(e.refs.iter().copied());
            for k in kids.into_iter().rev() {
                stack.push(k);
            }
        }

        let new_of: HashMap<usize, usize> =
            order.iter().enumerate().map(|(i, old)| (*old, i + 1)).collect();
        order
            .iter()
            .map(|old| {
                let e = &self.entities[old - 1];
                Ent {
                    ty: e.ty,
                    form: e.form,
                    // `v[6] = themodel->DNum(trsf)` (`IGESData_IGESWriter.cxx:324-331`):
                    // the transformation matrix field is an entity pointer like
                    // any other, so it carries the first Directory card number.
                    // A pointer that is not in `new_of` cannot happen (every
                    // child was pushed onto the DFS stack); keep the old value
                    // rather than emitting a hole.
                    trsf: e.trsf.map(|t| {
                        new_of
                            .get(&t)
                            .map(|n| 2 * n - 1)
                            .unwrap_or(t)
                    }),
                    refs: e
                        .refs
                        .iter()
                        .map(|r| new_of.get(r).copied().unwrap_or(*r))
                        .collect(),
                    params: resolve_placeholders(&e.params, &e.refs, &new_of),
                    status: [0; 4],
                }
            })
            .collect()
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

    /// Entity 108 (`GeomToIGES_GeomSurface::TransferSurface(Geom_Plane)`,
    /// `GeomToIGES_GeomSurface.cxx:609-688`, the `write.iges.plane.mode == 0`
    /// arm - that static is initialised to 0, `IGESData.cxx:182-186`; written by
    /// `IGESGeom_ToolPlane::WriteOwnParams`, `IGESGeom_ToolPlane.cxx:129-145`):
    ///
    /// `108, A, B, C, D, bounding_curve, attach_x, attach_y, attach_z, symbol_size;`
    ///
    /// - `A,B,C,D` are `start->Coefficients(...)` with `D` negated (`cxx:625-627`)
    ///   - the plane as `A*x + B*y + C*z = D`, all divided by the file unit,
    ///   which is what the port's `d = n . location` already is;
    /// - the bounding curve is null and the symbol size is 0 (`cxx:630` passes
    ///   `nullptr` and `0`), so both write `0`;
    /// - the attachment point is the plane's location / unit (`cxx:629`);
    /// - the form stays 0: `TransferPlaneSurface` never calls `SetFormNumber`,
    ///   and `IGESGeom_Plane::Init` only re-applies the form already stored
    ///   (`IGESGeom_Plane.cxx:44`), which is 0 on a fresh entity.
    fn emit_plane(&mut self, pln: &GpPln) -> usize {
        // `GeomToIGES_GeomSurface.cxx:625-629`: `start->Coefficients(A,B,C,D)`
        // then `D = -D` "because of difference in Geom_Plane class and Type 108".
        // `gp_Pln::Coefficients` (`gp_Pln.hxx`) negates the stored main direction
        // when the placement is indirect, which is exactly what
        // `Geom_Plane::UReverse` (`gp_Pln.hxx:101`, `XReverse`) produces for a
        // REVERSED face.
        let (a, b, c, cte) = pln.coefficients();
        let d = -cte;
        let loc = pln.location().coord;
        self.emit(
            108,
            0,
            format!(
                "108,{},{},{},{},0,{},{},{},0.;",
                num(a),
                num(b),
                num(c),
                num(d),
                num(loc.x),
                num(loc.y),
                num(loc.z)
            ),
        )
    }

    /// Entity 100 (`GeomToIGES_GeomCurve::TransferCurve(Geom_Circle)`,
    /// `GeomToIGES_GeomCurve.cxx:533-603`, written by
    /// `IGESGeom_ToolCircularArc::WriteOwnParams`,
    /// `IGESGeom_ToolCircularArc.cxx:79-89`):
    ///
    /// `100, ZT, Xc, Yc, Xs, Ys, Xe, Ye;`
    ///
    /// Every coordinate is **2-D, in the arc's own plane**: `TransferCircle`
    /// places the builder on `gp_Ax3(arc.Position())` (`cxx:560-563`) and writes
    /// each point through `Build.EvalXYZ` (`cxx:580-583`). The placing frame
    /// itself rides on the entity's transformation matrix (`cxx:591-597`).
    fn emit_circular_arc(&mut self, frame: &GpAx3, radius: f64, u1: f64, u2: f64) -> usize {
        // `TransferCurve` (`GeomToIGES_GeomCurve.cxx:569-586`) evaluates the
        // centre and the two ends with `Geom_Circle::D0`, then
        // `IGESConvGeom_GeomBuilder::EvalXYZ` (`GeomBuilder.cxx:212-216`):
        // `thepos.Inverted().Transforms` after `SetPosition(gp_Ax3)`
        // (`:137-143`). `gp_Trsf::Invert` transposes (`gp_Trsf.cxx:418-425`).
        // Writing `R*(cos U, sin U)` skips that round trip; a frame whose axes
        // carry a sub-resolution component cancels the libm residual to a real
        // 0 only after the inverse is applied. `GetUnit()` is 1 for this writer
        // (the 124 below is emitted with unit 1).
        let local = builder_local_trsf(frame);
        let (xc, yc, zc) = eval_xyz(&local, &frame.location());
        let (xs, ys, _) = eval_xyz(&local, &circle_point(frame, radius, u1));
        let (xe, ye, _) = eval_xyz(&local, &circle_point(frame, radius, u2));
        let de = self.emit(
            100,
            0,
            format!(
                "100,{},{},{},{},{},{},{};",
                num(zc),
                num(xc),
                num(yc),
                num(xs),
                num(ys),
                num(xe),
                num(ye)
            ),
        );
        if !frame_is_identity(frame) {
            let trsf = self.emit_transformation_matrix(frame, 1.0);
            self.set_trsf(de, trsf);
        }
        de
    }

    /// `GeomToIGES_GeomCurve::TransferCurve(Geom_Circle, Udeb, Ufin)`
    /// (`GeomToIGES_GeomCurve.cxx:533-603`): before the arc is built a `Udeb`
    /// within `gp::Resolution()` snaps to 0 (`:549-552`), and a request spanning
    /// a whole period writes the end point coincident with the start
    /// (`:571-576`, `gka BUG 6542`).
    fn emit_circle_transfer(&mut self, frame: &GpAx3, radius: f64, u1: f64, u2: f64) -> usize {
        let first = if u1.abs() <= occt_core::precision::REAL_SMALL {
            0.0
        } else {
            u1
        };
        let last = if (u2 - u1 - 2.0 * std::f64::consts::PI).abs()
            <= occt_core::precision::PCONFUSION
        {
            first
        } else {
            u2
        };
        self.emit_circular_arc(frame, radius, first, last)
    }

    // ---- topology entities ----

    /// `TransferCurve(Geom_OffsetCurve)` under the default
    /// `write.iges.offset.mode = 0` (`IGESData.cxx:193`,
    /// `GeomToIGES_GeomCurve.cxx:916-919`):
    /// `TransferCurve(GeomConvert::CurveToBSplineCurve(start), U1, U2)`.
    ///
    /// `TransferCurve(Geom_TrimmedCurve)` forwards to the immediate basis with
    /// the same range (`cxx:467-475`; the nested-trim assignment at `:472` is
    /// overwritten, so only the basis that recursion reaches is written). The
    /// port's `GeomTrimmedCurve` evaluates on `[0, 1]` while
    /// `untrimmed_basis` reports the basis interval (`trimmed.rs`), so a
    /// remapped trim's request is mapped back before the spline is segmented.
    /// The mode `1` arm (entity 130) is not the default and is not taken.
    /// A failed approximation is the `Standard_ConstructionError` OCCT throws
    /// (`GeomConvert.cxx:350-352`); the caller then keeps its null-handle path.
    fn emit_offset_curve_default(
        &mut self,
        curve: &dyn Curve,
        a: f64,
        b: f64,
    ) -> Option<usize> {
        let (target, u0, u1) = offset_curve_request(curve, a, b)?;
        let bs = occt_geom::convert_bspl::curve_to_bspline_curve(
            target.as_ref(),
            occt_core::convert::ParameterisationType::TgtThetaOver2,
        )
        .ok()?;
        self.emit_bspline_curve(&bs, u0, u1)
    }

    /// `GeomToIGES_GeomCurve::TransferCurve(Geom_Curve, Udeb, Ufin)`
    /// (`GeomToIGES_GeomCurve.cxx:94-126` dispatching, `:133-161` on the bounded
    /// family, `:481-528` on the conics): a `Geom_Line` becomes entity 110, a
    /// `Geom_Circle` entity 100, a B-spline/Bezier entity 126. Returns `None`
    /// where OCCT's transfer returns a null handle (conics: the 104 writer is not
    /// ported yet; any other curve type has no branch at all).
    ///
    /// A `Geom_OffsetCurve` (and a trimmed curve whose basis is one,
    /// `cxx:467-475`) is taken first, matching the default
    /// `write.iges.offset.mode = 0` arm (`cxx:916-919`).
    fn emit_curve_range(&mut self, curve: &dyn Curve, a: f64, b: f64) -> Option<usize> {
        if let Some(idx) = self.emit_offset_curve_default(curve, a, b) {
            return Some(idx);
        }
        match iges_curve_kind(curve) {
            IgCurveKind::Line => {
                let p1 = curve.d0(a);
                let p2 = curve.d0(b);
                Some(self.emit_line(&p1, &p2))
            }
            // `GeomToIGES_GeomCurve::TransferCircle`
            // (`GeomToIGES_GeomCurve.cxx:533-603`): the centre, the axis and the
            // radius come from the `Geom_Circle`'s own `gp_Circ` - never from a
            // three-point fit on sampled points - and the start/end parameters are
            // used as they are, apart from the two special cases in
            // `emit_circle_transfer`.
            IgCurveKind::Circle => {
                let c = curve.gp_circ()?;
                Some(self.emit_circle_transfer(&c.pos.to_ax3(), c.radius, a, b))
            }
            // `TransferCurve(Geom_BSplineCurve)` (`cxx:279-423`): a periodic curve
            // is written through a `SetNotPeriodic` copy (`emit_bspline_curve`).
            // UNPORTED: a requested range narrower than the curve's own is
            // obtained with `Geom_BSplineCurve::Segment`, so those cases return
            // `None`.
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

    /// `BRepToIGES_BRWire::TransferEdge(edge, originMap, false)`
    /// (`BRepToIGES_BRWire.cxx:266-335`): the edge's 3-D curve entity.
    ///
    /// Returns `None` where OCCT leaves `ICurve` null and the caller skips the
    /// edge (`TransferWire`, `cxx:715-723`): the edge has no 3-D curve at all
    /// (a degenerated edge - `BRep_Builder::Degenerated` drops the curve,
    /// `BRep_Builder.cxx:1073-1085`, so `BRep_Tool::Curve` at `cxx:282` returns
    /// null; the port's `EdgeGeom::curve` is non-nullable and carries the
    /// `degenerated` flag instead, see `builder.rs:326-331`).
    fn emit_edge_curve(&mut self, e: &Edge) -> Option<usize> {
        let Some(curve) = BRepTool::edge_curve(e) else {
            let (v1, v2) = edge_vertices(e);
            let p1 = v1.map(|v| BRepTool::vertex_point(&v)).unwrap_or_default();
            let p2 = v2.map(|v| BRepTool::vertex_point(&v)).unwrap_or_default();
            return Some(self.emit_line(&p1, &p2));
        };
        if BRepTool::is_degenerated(e) {
            return None;
        }
        let (first, last) = BRepTool::edge_parameters(e);
        // `BRepToIGES_BRWire::TransferEdge` (`cxx:302-314`): every call from the
        // face path passes `theIsBRepMode = false`, so a REVERSED edge is written
        // through `Curve3d->Reverse()` with the range mapped by
        // `ReversedParameter(Last)` / `ReversedParameter(First)`.
        let (curve, a, b): (Arc<dyn Curve>, f64, f64) =
            if e.0.orientation() == crate::abs::Orientation::Reversed {
                let u1 = curve.reversed_parameter(last);
                let u2 = curve.reversed_parameter(first);
                (Arc::from(curve.reversed()), u1, u2)
            } else {
                (curve, first, last)
            };
        match self.emit_curve_range(curve.as_ref(), a, b) {
            Some(idx) => Some(idx),
            None => {
                // The transfer above returned a null handle: OCCT writes no curve
                // at all for those types, the port keeps the chord-line stand-in
                // (UNPORTED, audit A26 / task T-78, see `emit_curve_range`).
                let p1 = curve.d0(a);
                let p2 = curve.d0(b);
                Some(self.emit_line(&p1, &p2))
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
    /// sequence. The planar flag and the normal come from `IsPlanar(start, Norm)`
    /// (`cxx:414` -> `:236-272`), which runs on the curve **as received**, not on
    /// the unperiodized / `Segment`-ed copy whose poles are written out: for a
    /// B-spline that is `ArePolesPlanar(Poles())` (`cxx:191-228`),
    /// `P(n) x P(1) + sum P(i) x P(i+1)`, normalised (or `(0,0,1)` with
    /// `planar = false` when the sum is shorter than `Precision::Confusion()`),
    /// then every pole must sit at the same distance from the plane through
    /// `P(1)`. Returns `None` for the cases this port cannot express.
    fn emit_bspline_curve(&mut self, start: &dyn Curve, u_deb: f64, u_fin: f64) -> Option<usize> {
        // `cxx:294-307`: a periodic curve is written through a non-periodic copy
        // (`SetNotPeriodic`). `IPerio` is captured from the curve as received
        // (`cxx:297`) and passed to `Init` (`cxx:419`) unchanged, so the 126
        // parameter stays 1 after the copy is opened.
        let unperiodized;
        let curve: &dyn Curve = if start.is_periodic() {
            let poles = start.bspline_poles()?.to_vec();
            let knots = start.bspline_knots()?.to_vec();
            let degree = start.nurbs_degree()?;
            let mut c = occt_geom::bspline_curve::GeomBSplineCurve {
                poles,
                weights: start.bspline_weights().map(|w| w.to_vec()),
                knots,
                degree,
                periodic: true,
            };
            c.set_not_periodic();
            unperiodized = c;
            &unperiodized
        } else {
            start
        };

        // `cxx:309-318`: an infinite bound is replaced by `±Precision::Infinite()`.
        let mut umin = if occt_core::precision::Precision::is_negative_infinite(u_deb) {
            -occt_core::precision::INFINITE
        } else {
            u_deb
        };
        let mut umax = if occt_core::precision::Precision::is_positive_infinite(u_fin) {
            occt_core::precision::INFINITE
        } else {
            u_fin
        };

        // `cxx:320-331`: protect against exceptions in `Segment()`; the range is
        // pulled onto the curve's own bounds (one-sided, as OCCT writes the
        // clamped values into the entity's `UMin`/`UMax`).
        let (curve_first, curve_last) = (curve.first_parameter(), curve.last_parameter());
        if umin - curve_first < occt_core::precision::PCONFUSION {
            umin = curve_first;
        }
        if curve_last - umax < occt_core::precision::PCONFUSION {
            umax = curve_last;
        }

        // `cxx:332-356`: cut the curve for E3 — a copy is `Segment`-ed and
        // replaces the original. OCCT catches the failure and keeps the
        // untrimmed copy; a null copy likewise leaves the curve unchanged.
        let trimmed;
        let curve: &dyn Curve = if umin - curve_first > occt_core::precision::PCONFUSION
            || curve_last - umax > occt_core::precision::PCONFUSION
        {
            // UNPORTED: `TransferCurve(Geom_BezierCurve)` (`cxx:441-448`) puts
            // the Bezier in a `Geom_TrimmedCurve` and converts it through
            // `GeomConvert::CurveToBSplineCurve`, which needs
            // `Geom_BezierCurve::Segment` (not ported) — a narrow request on a
            // Bezier keeps returning `None` (the caller writes a chord line).
            let bs = bspline_copy_of(curve)?;
            let mut bs = bs;
            if (umax - umin).abs() > occt_core::precision::PCONFUSION {
                // default `theTolerance = Precision::PConfusion()`.
                let _ = bs.segment(umin, umax, occt_core::precision::PCONFUSION);
            }
            trimmed = bs;
            &trimmed
        } else {
            curve
        };

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


        let (planar, normal) = curve_planar_and_normal(start);
        // `cxx:415-418`: the normal is flipped when it points down.
        let normal = if normal.z() < 0.0 {
            GpVec::new(-normal.x(), -normal.y(), -normal.z())
        } else {
            normal
        };

        // `cxx:360`: `IPolyn = !mycurve->IsRational()`. `IsRational` is the
        // `myRational` flag (`Geom_BSplineCurve_1.cxx:683-686`), set by
        // `Rational(Weights)` (`Geom_BSplineCurve.cxx:98-108`, `:206-208`).
        // A stored weight array whose values are all equal is not rational, and
        // the entity then carries `BSplCLib::UnitWeights` (`:220-222`).
        //
        // `WriteOwnParams` does not send that flag. It sends
        // `IGESGeom_BSplineCurve::IsPolynomial()` (`IGESGeom_BSplineCurve.cxx:105-124`,
        // called from `IGESGeom_ToolBSplineCurve.cxx:242`), which reports
        // polynomial when every weight is within `1e-10` of the first one.
        let stored = curve.bspline_weights().filter(|w| w.len() == poles.len());
        let geom_rational = stored.as_ref().is_some_and(|w| weights_are_rational(w));
        let weights: Vec<f64> = if geom_rational {
            stored.unwrap().to_vec()
        } else {
            vec![1.0; poles.len()]
        };
        let polynomial = !geom_rational
            || weights
                .first()
                .is_none_or(|w0| weights.iter().all(|wi| (wi - w0).abs() <= 1.0e-10));
        let closed = poles
            .first()
            .zip(poles.last())
            .map(|(f, l)| f.distance(l) <= occt_core::precision::CONFUSION)
            .unwrap_or(false);

        let mut s = format!(
            "126,{index},{degree},{},{},{},{}",
            u8::from(planar),
            u8::from(closed),
            u8::from(polynomial),
            u8::from(start.is_periodic())
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
        // `cxx:419`: the entity carries the clamped/segmented `Umin`/`Umax`,
        // not the raw requested range.
        s.push_str(&format!(
            ",{},{},{},{},{};",
            num(umin),
            num(umax),
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
    ///
    /// **Not wired**: this is the BRep-mode arm (`myBRepMode && myAnalytic`,
    /// `cxx:579-586`). OCCT's default `write.iges.brep.mode = 0`
    /// (`IGESData.cxx:90-94`) keeps the writer in Faces mode, where a sphere
    /// becomes a 120 revolution surface instead (see
    /// [`IgesWriter::emit_local_revolution_surface`]).
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
        self.emit_refs(
            196,
            0,
            format!("196,#0,{},#1,#2;", num(radius)),
            &[c, a, r],
        )
    }

    /// Entity 192 (`GeomToIGES_GeomSurface::TransferCylindricalSurface`,
    /// `GeomToIGES_GeomSurface.cxx:1280-1314`, written by
    /// `IGESSolid_ToolCylindricalSurface::WriteOwnParams`):
    /// `192, location_point, axis_direction, radius, reference_direction;`.
    ///
    /// **Not wired**: BRep-mode only (`myBRepMode && myAnalytic`, `cxx:555-562`);
    /// in Faces mode a cylinder is a 120.
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
        self.emit_refs(
            192,
            0,
            format!("192,#0,#1,{},#2;", num(radius)),
            &[l, a, r],
        )
    }

    /// Entity 194 (`TransferConicalSurface`, `cxx:1318-1362`, written by
    /// `IGESSolid_ToolConicalSurface::WriteOwnParams`):
    /// `194, location_point, axis_direction, ref_radius, semi_angle_deg,
    /// reference_direction;`. A negative semi-angle is written by mirroring the
    /// reference point through the apex, negating the angle and reversing the
    /// reference direction (`cxx:1344-1350`).
    ///
    /// **Not wired**: BRep-mode only (`myBRepMode && myAnalytic`, `cxx:567-574`);
    /// in Faces mode a cone is a 120.
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
        self.emit_refs(
            194,
            0,
            format!(
                "194,#0,#1,{},{},#2;",
                num(ref_radius),
                num(angle * 180.0 / std::f64::consts::PI)
            ),
            &[l, a, r],
        )
    }

    /// Entity 198 (`TransferToroidalSurface`, `cxx:1402-1435`, written by
    /// `IGESSolid_ToolToroidalSurface::WriteOwnParams`):
    /// `198, centre_point, axis_direction, major_radius, minor_radius,
    /// reference_direction;`.
    ///
    /// **Not wired**: BRep-mode only (`myBRepMode && myAnalytic`, `cxx:591-598`);
    /// in Faces mode a torus is a 120.
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
        self.emit_refs(
            198,
            0,
            format!("198,#0,#1,{},{},#2;", num(major), num(minor)),
            &[c, a, r],
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
    /// [`iso_rows_equal`]) of the **unperiodized** surface, `polynomial` from
    /// `Polynom = !(RationU || RationV)` (`cxx:355`; the port's B-spline surface
    /// carries one rational flag for both directions), and `periodicU/V` from the
    /// original `IsUPeriodic`/`IsVPeriodic` (`cxx:235-236`, `:448-449`). The
    /// bounds fix implements both arms (`cxx:244-302`) and a periodic direction is
    /// unperiodized before its knots and poles are read (`cxx:303-343`,
    /// `SetUNotPeriodic`/`SetVNotPeriodic` → `BSplSLib::Unperiodize`,
    /// [`occt_core::bspl::unperiodize_direction`]). The only UNPORTED piece is the
    /// `SetUOrigin`/`SetVOrigin` re-origin (`cxx:310-320`, `:330-340`), see the
    /// note inside.
    fn emit_bspline_surface(
        &mut self,
        surf: &dyn Surface,
        u0: f64,
        u1: f64,
        v0: f64,
        v1: f64,
    ) -> Option<usize> {
        let (mut poles, mut knots_u, mut knots_v, deg_u, deg_v) = match (
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
            _ => {
                if std::env::var("IGES_TRACE").is_ok() {
                    eprintln!("TRACE bspline missing data: is_bspline={}", surf.is_bspline_surface());
                }
                return None;
            }
        };
        if poles.is_empty() || poles[0].is_empty() {
            return None;
        }
        let (nu0, nv0) = (poles.len(), poles[0].len());

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
        let mut weights: Vec<Vec<f64>> = match stored_weights {
            Some(w) if w.len() == nu0 && w.iter().all(|r| r.len() == nv0) => w.to_vec(),
            _ => vec![vec![1.0; nv0]; nu0],
        };
        let polynomial = !rational;

        // `GeomToIGES_GeomSurface.cxx:235-236`: `PeriodicU/V` are read from the
        // **original** surface and written verbatim (`:448-449`).
        let period_u = bs.is_u_periodic();
        let period_v = bs.is_v_periodic();

        // An **open** direction's flat knot vector and pole grid must describe
        // the same representation: `nb_poles + degree + 1`. A periodic direction
        // is *not* trimmed to the open length -
        // `Geom_BSplineSurface::SetUPeriodic` (`Geom_BSplineSurface_1.cxx:940-981`)
        // trims the poles to `BSplCLib::NbPoles(degree, true, mults)` while
        // `updateUKnots` rebuilds the flat sequence with the periodic extension
        // (`BSplCLib::KnotSequenceLength`, `BSplCLib.cxx:455-474`) - and the base
        // multiplicities that would let the length be predicted are no longer
        // recoverable from the extended sequence (the leading `degree + 1 -
        // mults(1)` entries are the wrapped tail, shifted by one period). Such a
        // direction is therefore checked after the unperiodization below, on the
        // open representation OCCT's writer actually emits.
        let open_flat_ok = |flat: &[f64], nb: usize, deg: usize, periodic: bool| -> bool {
            periodic || flat.len() == nb + deg + 1
        };
        if !open_flat_ok(&knots_u, nu0, deg_u, period_u)
            || !open_flat_ok(&knots_v, nv0, deg_v, period_v)
        {
            if std::env::var("IGES_TRACE").is_ok() {
                eprintln!(
                    "TRACE bspline layout mismatch: nu={nu0} nv={nv0} deg=({deg_u},{deg_v}) ku={} kv={} | nbuv=({},{}) uvdeg=({},{}) pu={period_u} pv={period_v} trim={} bs={}",
                    knots_u.len(),
                    knots_v.len(),
                    surf.nb_u_poles(),
                    surf.nb_v_poles(),
                    surf.u_degree(),
                    surf.v_degree(),
                    surf.rectangular_trimmed_basis().is_some(),
                    surf.is_bspline_surface(),
                );
            }
            return None;
        }

        // `cxx:239-301`: the written range is clamped to the surface's own
        // bounds; the periodic arm snaps an end that already sits on a bound and
        // otherwise shifts the range into the period
        // (`ShapeAnalysis::AdjustToPeriod`), truncating it to one period. This
        // runs **before** the re-origin below, which needs `uShift`/`vShift`.
        let (su0, su1) = bs.u_range();
        let (sv0, sv1) = bs.v_range();
        let (u1_arg, v1_arg) = (u1, v1);
        let (mut u0, mut u1) = (u0, u1);
        let mut u_shift = 0.0;
        if period_u {
            if (u0 - su0).abs() < occt_core::precision::PCONFUSION {
                u0 = su0;
            }
            if (u1 - su1).abs() < occt_core::precision::PCONFUSION {
                u1 = su1;
            }
            u_shift = crate::pcurve_full::adjust_to_period(u0, su0, su1);
            u0 += u_shift;
            u1 += u_shift;
            if u1 - u0 > su1 - su0 {
                u1 = u0 + (su1 - su0);
            }
        } else {
            u0 = u0.max(su0);
            u1 = u1.min(su1);
        }
        let (mut v0, mut v1) = (v0, v1);
        let mut v_shift = 0.0;
        if period_v {
            if (v0 - sv0).abs() < occt_core::precision::PCONFUSION {
                v0 = sv0;
            }
            if (v1 - sv1).abs() < occt_core::precision::PCONFUSION {
                v1 = sv1;
            }
            v_shift = crate::pcurve_full::adjust_to_period(v0, sv0, sv1);
            v0 += v_shift;
            v1 += v_shift;
            if v1 - v0 > sv1 - sv0 {
                v1 = v0 + (sv1 - sv0);
            }
        } else {
            v0 = v0.max(sv0);
            v1 = v1.min(sv1);
        }

        // `cxx:303-343`: `SetUOrigin`/`SetVOrigin` re-origin a periodic
        // B-spline direction when the written range's own period shift differs
        // from the range's start one (issue 26138: synchronize the p-curve
        // ranges with the surface bounds). `LocateU(Vmin)` gives the window
        // index the new origin starts at; the re-origin only permutes the
        // periodic representation, so the geometry is unchanged.
        let u_reorigin = period_u
            && (u_shift - crate::pcurve_full::adjust_to_period(u1_arg, su0, su1)).abs()
                > occt_core::precision::PCONFUSION;
        let v_reorigin = period_v
            && (v_shift - crate::pcurve_full::adjust_to_period(v1_arg, sv0, sv1)).abs()
                > occt_core::precision::PCONFUSION;
        if u_reorigin || v_reorigin {
            let mut ms = occt_geom::bspline_surface::GeomBSplineSurface {
                poles: poles.clone(),
                knots_u: knots_u.clone(),
                knots_v: knots_v.clone(),
                deg_u,
                deg_v,
                weights: if rational { Some(weights.clone()) } else { None },
                u_periodic: period_u,
                v_periodic: period_v,
            };
            let mut changed = false;
            if u_reorigin {
                let (left, _right) = ms.locate_u(u0, occt_core::precision::PCONFUSION);
                if ms.set_u_origin(left).is_ok() {
                    changed = true;
                }
            }
            if v_reorigin {
                let (left, _right) = ms.locate_v(v0, occt_core::precision::PCONFUSION);
                if ms.set_v_origin(left).is_ok() {
                    changed = true;
                }
            }
            if changed {
                poles = ms.poles;
                if let Some(w) = ms.weights {
                    weights = w;
                }
                knots_u = ms.knots_u;
                knots_v = ms.knots_v;
            }
        }

        // `cxx:303-343`: a periodic B-spline surface is **unperiodized** before
        // its knots and poles are read (`SetUNotPeriodic` / `SetVNotPeriodic` →
        // `BSplSLib::Unperiodize`, `Geom_BSplineSurface_1.cxx:1238-1300`), so
        // the written 128 has the open knot vector IGES expects. The pole
        // extension is the cyclic wrap of `BSplCLib::Unperiodize`
        // (`BSplCLib.cxx:3076-3079`, modulus `Poles.Length()` =
        // `BSplCLib::NbPoles(degree, true, mults)`); weights follow the same wrap
        // because `BSplSLib::Unperiodize` works on homogeneous poles.
        if period_u {
            let (nf, map) =
                occt_core::bspl::unperiodize_direction(deg_u as i32, &knots_u, nu0);
            knots_u = nf;
            poles = map.iter().map(|k| poles[*k].clone()).collect();
            weights = map.iter().map(|k| weights[*k].clone()).collect();
        }
        if period_v {
            let (nf, map) =
                occt_core::bspl::unperiodize_direction(deg_v as i32, &knots_v, nv0);
            knots_v = nf;
            for row in poles.iter_mut() {
                *row = map.iter().map(|k| row[*k].clone()).collect();
            }
            for row in weights.iter_mut() {
                *row = map.iter().map(|k| row[*k]).collect();
            }
        }

        let (nu, nv) = (poles.len(), poles[0].len());
        if knots_u.len() != nu + deg_u + 1 || knots_v.len() != nv + deg_v + 1 {
            if std::env::var("IGES_TRACE").is_ok() {
                eprintln!(
                    "TRACE bspline post-unperiodize mismatch: nu={nu} nv={nv} deg=({deg_u},{deg_v}) ku={} kv={} pu={period_u} pv={period_v}",
                    knots_u.len(),
                    knots_v.len()
                );
            }
            return None;
        }
        let (ind_u, ind_v) = (nu - 1, nv - 1);

        // `Geom_BSplineSurface_1.cxx:1026-1130` `SetUOrigin`/`SetVOrigin` — the
        // re-origin OCCT performs when the written range crosses the period
        // origin (`cxx:310-320`, `:330-340`) — is UNPORTED (audit A26 / T-78):
        // it only differs when `AdjustToPeriod(Umin, U0, U1)` differs from
        // `AdjustToPeriod(Ufin, U0, U1)`, and the port keeps the surface's own
        // knot origin in that case.
        //
        // `CloseU = mysurface->IsUClosed()` / `CloseV` (`cxx:347-348`), read from
        // the **unperiodized** surface: a periodic direction whose first and last
        // pole rows coincide is closed, which the pole comparison below detects
        // after the cyclic extension.
        let closed_u = iso_rows_equal(&poles[0], &poles[nu - 1])
            && (!rational || weights[0].iter().zip(&weights[nu - 1]).all(|(a, b)| weight_equal(*a, *b)));
        let closed_v = (0..nu).all(|i| {
            iso_rows_equal(&poles[i][0..1], &poles[i][nv - 1..nv])
                && (!rational || weight_equal(weights[i][0], weights[i][nv - 1]))
        });

        let mut s = format!(
            "128,{ind_u},{ind_v},{deg_u},{deg_v},{},{},{},{},{}",
            u8::from(closed_u),
            u8::from(closed_v),
            u8::from(polynomial),
            u8::from(period_u),
            u8::from(period_v)
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
        Some(self.emit_refs(
            120,
            0,
            format!(
                "120,#0,#1,{},{};",
                num(tau - u1),
                num(tau - u0)
            ),
            &[axis_line, generatrix],
        ))
    }

    /// `GeomToIGES_GeomSurface::TransferSurface(Geom_CylindricalSurface /
    /// Geom_ConicalSurface / Geom_SphericalSurface / Geom_ToroidalSurface)`
    /// (`GeomToIGES_GeomSurface.cxx:696-768`, `:777-852`, `:859-925`,
    /// `:933-993`).
    ///
    /// These four overloads are reached only from the `myBRepMode && myAnalytic`
    /// arm of the `Geom_ElementarySurface` branch (`cxx:552-598`). Under OCCT's
    /// default `write.iges.brep.mode = 0` (`IGESData.cxx:90-94`) the writer runs
    /// in Faces mode, so they are **not** taken and the elementary surfaces are
    /// written through the plain `TransferSurface(Geom_Surface)` overload as
    /// `IGESGeom_SurfaceOfRevolution` (120). This helper is that Faces-mode 120:
    /// the generatrix and the axis are built in the surface's **local** frame and
    /// the frame itself rides on a 124.
    ///
    /// - the axis is always the local Z axis, written as an `IGESGeom_Line` from
    ///   `(0,0,1)` to `(0,0,0)` (`cxx:744`, `:830`, `:902`, `:969`; the
    ///   `#30 rln 19.10.98` "IGES axis = reversed CAS.CADE axis");
    /// - `generatrix` is `GC.TransferCurve(<local generatrix>, V1, V2)` already
    ///   written by the caller (`cxx:819`, `:820`, `:896`, `:964`);
    /// - the angles are `Init(Axis, Generatrix, 2*PI - U2, 2*PI - U1)`
    ///   (`cxx:831`, `:907`, `:971`);
    /// - `IGESConvGeom_GeomBuilder::SetPosition(<surface>.Position())` plus
    ///   `IsIdentity()` decide whether the 124 is written
    ///   (`cxx:841-849`, `:915-921`, `:981-988`;
    ///   `IGESConvGeom_GeomBuilder.cxx:137-198`).
    fn emit_local_revolution_surface(
        &mut self,
        generatrix: usize,
        frame: &GpAx3,
        u0: f64,
        u1: f64,
    ) -> usize {
        let axis = self.emit_line(&GpPnt::new(0.0, 0.0, 1.0), &GpPnt::new(0.0, 0.0, 0.0));
        let tau = 2.0 * std::f64::consts::PI;
        let de = self.emit_refs(
            120,
            0,
            format!("120,#0,#1,{},{};", num(tau - u1), num(tau - u0)),
            &[axis, generatrix],
        );
        if !frame_is_identity(frame) {
            let trsf = self.emit_transformation_matrix(frame, 1.0);
            self.set_trsf(de, trsf);
        }
        de
    }

    /// Entity 122 (`GeomToIGES_GeomSurface::TransferSurface(Geom_SurfaceOfLinearExtrusion)`,
    /// `GeomToIGES_GeomSurface.cxx:1032-1104`, written by
    /// `IGESGeom_ToolTabulatedCylinder::WriteOwnParams`, `:90-98`):
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
        Some(self.emit_refs(
            122,
            0,
            format!(
                "122,#0,{},{},{};",
                num(end.x()),
                num(end.y()),
                num(end.z())
            ),
            &[directrix],
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
        // `IGESConvGeom_GeomBuilder::SetPosition(const gp_Ax3&)`
        // (`GeomBuilder.cxx:137-143`) runs `ps.SetTransformation(pos, gp::XOY())`,
        // then `MakeTransformation` (`GeomBuilder.cxx:218-237`) reads
        // `thepos.Value(i, j)` and divides column 4 by the file unit. Going
        // through `gp_Trsf::SetTransformation` (`gp_Trsf.cxx:172-194`) matters:
        // the trailing `matrix.Multiply(MA1)` (identity times the frame's basis
        // columns) normalises the signed zeros of `Value(i, j)`, and `IsNegative`
        // (`gp_Trsf::IsNegative`, `gp_Trsf.cxx`) is the form-number test.
        let mut ps = GpTrsf::identity();
        ps.set_transformation_from_to(frame, &xoy_frame());
        let form = if ps.is_negative() { 1 } else { 0 };
        let mut s = String::from("124");
        for i in 1..=3 {
            for j in 1..=3 {
                s.push(',');
                s.push_str(&num(ps.value(i, j)));
            }
            s.push(',');
            s.push_str(&num(ps.value(i, 4) / unit));
        }
        s.push(';');
        self.emit(124, form, s)
    }

    /// `GeomToIGES_GeomCurve::TransferCurve(Geom_Ellipse)` full-period arm
    /// (`GeomToIGES_GeomCurve.cxx:620-645`): the ellipse copy is rotated by
    /// `Udeb` (`:626-628`, direction chosen by `gp_Ax3(pos).Direct()`), converted
    /// to a B-spline, re-parameterised onto `[Udeb, Udeb + 2*PI]` (`:641-643`)
    /// and transferred as entity 126.
    ///
    /// `GeomConvert_ApproxCurve(aCopy, Precision::Approximation(), GeomAbs_C1,
    /// 100, 6)` is tried first (`:632-636`); `GeomConvert::CurveToBSplineCurve
    /// (copystart, Convert_QuasiAngular)` (`:637-640`) is only the fallback when
    /// that produces no curve.
    fn emit_whole_period_ellipse(&mut self, e: &occt_core::gp::GpElips, a: f64, b: f64) -> Option<usize> {
        use std::f64::consts::PI;
        let pos = *e.position();
        // `copystart->SetPosition(pos.Rotated(pos.Axis(), gp_Ax3(pos).Direct() ? Udeb : 2*PI - Udeb))`
        // (`cxx:626-628`).
        let angle = if pos.to_ax3().is_direct() { a } else { 2.0 * PI - a };
        let mut copy = *e;
        copy.set_position(pos.rotated(&pos.axis().clone(), angle));
        let rotated = occt_geom::ellipse::GeomEllipse::new(copy);
        // `if (approx.HasResult()) Bspline = approx.Curve(); if (Bspline.IsNull())
        //  Bspline = GeomConvert::CurveToBSplineCurve(copystart, Convert_QuasiAngular);`
        // (`cxx:632-640`).
        let approx = occt_geom::GeomConvertApproxCurve::new(
            &rotated,
            occt_core::precision::APPROXIMATION,
            occt_core::kernel::geomabs::Shape::C1,
            100,
            6,
        );
        let mut bs = match approx.curve() {
            Some(curve) => curve.clone(),
            // `GeomConvert::CurveToBSplineCurve(copystart, Convert_QuasiAngular)` (`cxx:639`).
            None => occt_geom::convert_bspl::curve_to_bspline_curve(
                &rotated,
                occt_core::convert::ParameterisationType::QuasiAngular,
            )
            .ok()?,
        };
        // `Knots = Bspline->Knots(); BSplCLib::Reparametrize(Udeb, Udeb + 2*PI, Knots);
        //  Bspline->SetKnots(Knots);` (`cxx:641-643`).
        let (mut uknots, _) = bs.distinct_knots_and_mults();
        occt_core::bspl::knots::reparametrize(a, a + 2.0 * PI, &mut uknots);
        bs.set_knots(&uknots).ok()?;
        // `return TransferCurve(Bspline, Udeb, Ufin)` (`cxx:644`). The periodic
        // flag is read inside that transfer, before `SetNotPeriodic`.
        self.emit_bspline_curve(&bs, a, b)
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
    /// (`:692-697`). A full-period ellipse is routed to
    /// [`emit_whole_period_ellipse`](IgesWriter::emit_whole_period_ellipse)
    /// instead, exactly as OCCT does at `:620-645`.
    fn emit_conic_arc(&mut self, curve: &dyn Curve, a: f64, b: f64) -> Option<usize> {
        use std::f64::consts::PI;
        let unit = IGES_UNIT;
        let (abc, pos, u1, u2) = if let Some(e) = curve.gp_ellipse() {
            if (b - a - 2.0 * PI).abs() <= occt_core::precision::PCONFUSION {
                // `GeomToIGES_GeomCurve.cxx:620-645`: a trimmed **full-period**
                // ellipse (IGES 104 cannot carry the whole-period semantics) is
                // re-parameterised onto `[Udeb, Udeb + 2*PI]` and written as a
                // B-spline (entity 126) instead.
                return self.emit_whole_period_ellipse(&e, a, b);
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

    /// Base surface entity for a face.
    ///
    /// `BRepToIGES_BRShell::TransferFace` (`BRepToIGES_BRShell.cxx:246-262`)
    /// transfers the face's surface with
    /// `GeomToIGES_GeomSurface::TransferSurface(Surf, U1, U2, V1, V2)` and, when
    /// that returns a null handle, logs `"the basic surface is a null entity"`
    /// and **drops the face** (`:257-261`). There is no plane stand-in.
    ///
    /// `reversal` is `TransferFace`'s REVERSED branch: OCCT already replaced the
    /// face surface with the U-reversed, trim-stripped one (`:126-146`), mirrored
    /// every p-curve about `U = aCenter` and recomputed the UV bounds on the new
    /// face (`:253`). The port transfers the original face and mirrors the
    /// p-curves on the fly (`transfer_edge_uv`), so the box it hands to
    /// `TransferSurface` has to be the mirror of the original one -
    /// `[2*aCenter - U1, 2*aCenter - U0]`, V untouched - exactly the box OCCT
    /// reads off the mirrored copy. Handing the unmirrored box over instead turns
    /// every elementary arm's `2*PI - U2, 2*PI - U1` into a whole-period shift of
    /// the revolution range (a quarter cylinder at `[0, PI/2]` came out as
    /// `[3*PI/2, 2*PI]`).
    fn emit_face_surface(&mut self, f: &Face, reversal: Option<&FaceReversal>) -> Option<usize> {
        let surf = match reversal {
            Some(r) => r.surf.clone(),
            None => BRepTool::face_surface(f)?,
        };
        let (u0, u1, v0, v1) = face_uv_bounds_finite(f);
        let (u0, u1) = match reversal {
            Some(r) => (r.two_center - u1, r.two_center - u0),
            None => (u0, u1),
        };
        self.emit_surface_entity(surf.as_ref(), u0, u1, v0, v1)
    }

    /// `GeomToIGES_GeomSurface::TransferSurface(const Geom_Surface&, Udeb, Ufin,
    /// Vdeb, Vfin)` (`GeomToIGES_GeomSurface.cxx:112-147`): the surface type
    /// dispatcher.
    ///
    /// The dispatch order is OCCT's: the `Geom_RectangularTrimmedSurface` wrapper
    /// is stripped first (`Geom_BoundedSurface` family, `:125-128` -> `:492-515`,
    /// which recurses on `BasisSurface()`), then the exact `Geom_BoundedSurface` /
    /// `Geom_ElementarySurface` / `Geom_SweptSurface` / `Geom_OffsetSurface` types
    /// are tried in that order. No sampling is involved; a surface with no
    /// matching arm has no branch in OCCT either and returns a null handle there
    /// and `None` here.
    fn emit_surface_entity(
        &mut self,
        s: &dyn Surface,
        u0: f64,
        u1: f64,
        v0: f64,
        v1: f64,
    ) -> Option<usize> {
        // `TransferSurface(Geom_BoundedSurface)` (`cxx:492-515`): a rectangular
        // trim restarts the whole dispatch on its basis.
        if let Some(b) = s.rectangular_trimmed_basis() {
            return self.emit_surface_entity(b.as_ref(), u0, u1, v0, v1);
        }
        if let Some(pln) = s.gp_pln() {
            return Some(self.emit_plane(&pln));
        }
        if let Some(cy) = s.gp_cylinder() {
            // `cxx:815-817`: the generatrix is `Geom_Line(gp_Pnt(R,0,0), +Z)`, so
            // `Value(t) = (R, 0, t)`; `V1`/`V2` are the clamped `UVBounds` of the
            // face (`cxx:806-813`).
            let (v0, v1) = clamp_infinite_range(v0, v1);
            let r = cy.radius();
            let generatrix = self.emit_line(&GpPnt::new(r, 0.0, v0), &GpPnt::new(r, 0.0, v1));
            return Some(self.emit_local_revolution_surface(
                generatrix,
                &cy.position(),
                u0,
                u1,
            ));
        }
        if let Some(co) = s.gp_cone() {
            // `cxx:816-818`: `Geom_Line(gp_Pnt(RefRadius,0,0),
            // gp_Dir(sin a, 0, cos a))`, so `Value(t) = (RefRadius + t*sin a, 0,
            // t*cos a)`. The semi-angle is used as it is stored, sign included.
            let (v0, v1) = clamp_infinite_range(v0, v1);
            let (semi_sin, semi_cos) = co.semi_angle().sin_cos();
            let ref_radius = co.radius();
            let generatrix = self.emit_line(
                &GpPnt::new(ref_radius + v0 * semi_sin, 0.0, v0 * semi_cos),
                &GpPnt::new(ref_radius + v1 * semi_sin, 0.0, v1 * semi_cos),
            );
            return Some(self.emit_local_revolution_surface(
                generatrix,
                &co.position(),
                u0,
                u1,
            ));
        }
        if let Some(sp) = s.gp_sphere() {
            // `cxx:893-894`: the generatrix is the half circle
            // `Geom_Circle(gp_Ax2(gp::Origin(), -gp::DY(), gp::DX()), R)`
            // travelled over the face's own V range (`Su(-PI/2) .. Sv(+PI/2)`).
            let generatrix = self.emit_circle_transfer(
                &revolution_generatrix_frame(GpPnt::zero()),
                sp.radius(),
                v0,
                v1,
            );
            return Some(self.emit_local_revolution_surface(
                generatrix,
                &sp.position(),
                u0,
                u1,
            ));
        }
        if let Some(t) = s.gp_torus() {
            // `cxx:960-962`: `Geom_Circle(gp_Ax2(gp_Pnt(Major,0,0), -gp::DY(),
            // gp::DX()), Minor)`.
            let generatrix = self.emit_circle_transfer(
                &revolution_generatrix_frame(GpPnt::new(t.major_radius(), 0.0, 0.0)),
                t.minor_radius(),
                v0,
                v1,
            );
            return Some(self.emit_local_revolution_surface(
                generatrix,
                &t.position(),
                u0,
                u1,
            ));
        }
        // `GeomToIGES_GeomSurface::TransferSurface` (`cxx:154-188`) on the
        // Bounded family: the port writes entity 128.
        if let Some(idx) = self.emit_bspline_surface(s, u0, u1, v0, v1) {
            return Some(idx);
        }
        // `cxx:1000-1025` on the swept family - reached after the bounded family
        // (B-spline / Bezier / trimmed), which is the order `cxx:112-147`
        // dispatches in.
        match swept_surface_kind(s) {
            Some(SweptKind::Extrusion) => {
                if let Some(idx) = self.emit_tabulated_cylinder(s, v0, v1) {
                    return Some(idx);
                }
            }
            Some(SweptKind::Revolution) => {
                if let Some(idx) = self.emit_surface_of_revolution(s, u0, u1, v0, v1) {
                    return Some(idx);
                }
            }
            None => {}
        }
        // `cxx:1195-1237` on the offset family (entity 140 `OffsetSurface`). The
        // basis surface is put through this same dispatcher, so the 140 nests
        // whatever entity the basis maps to.
        if let (Some(basis), Some(distance)) = (s.offset_basis_surface(), s.offset_distance()) {
            // `cxx:1208-1214`: the indicator is read at the mid point of the
            // offset's own `Bounds`, which `Geom_OffsetSurface::Bounds` delegates
            // to its basis.
            let (bu0, bu1) = s.u_range();
            let (bv0, bv1) = s.v_range();
            let (um, vm) = (bounds_mid_point(bu0, bu1), bounds_mid_point(bv0, bv1));
            let (_, d1u, d1v) = basis.d1(um, vm);
            // `cxx:1217-1219`: `GeomLProp_SLProps(TheSurf, Um, Vm, 1,
            // Precision::Confusion()).Normal()`. OCCT's `Normal()` throws
            // `LProp_NotDefined` when `CSLib::Normal` does not report
            // `CSLib_Done`; the port reports the same condition as an unhandled
            // surface and lets the face be dropped.
            let normal = surface_normal_from_d1(&d1u, &d1v, occt_core::precision::CONFUSION)?;
            // `cxx:1215`: the basis is transferred with the *caller's* parameter
            // range, not with its own bounds.
            let surface = self.emit_surface_entity(basis.as_ref(), u0, u1, v0, v1);
            return Some(self.emit_offset_surface(&normal, distance / IGES_UNIT, surface));
        }
        // Everything else has no branch in `TransferSurface` and returns a null
        // handle: the face is dropped, never replaced by a plane.
        None
    }

    /// Entity 140 (`GeomToIGES_GeomSurface::TransferSurface(Geom_OffsetSurface)`,
    /// `GeomToIGES_GeomSurface.cxx:1195-1237`, written by
    /// `IGESGeom_ToolOffsetSurface::WriteOwnParams`,
    /// `IGESGeom_ToolOffsetSurface.cxx:105-113`):
    ///
    /// `140, indicator.x, indicator.y, indicator.z, distance, surface;`
    ///
    /// - the indicator is the basis surface's unit normal at the mid point of the
    ///   offset's bounds, every component divided by the file unit
    ///   (`cxx:1217-1226`);
    /// - the distance is the signed `Offset()` divided by the file unit
    ///   (`cxx:1216`);
    /// - the surface is the recursively transferred basis (`cxx:1215`). When that
    ///   transfer yields a null handle OCCT still builds the 140 and
    ///   `IGESData_IGESWriter::Send` writes the reference as `0`, so the port
    ///   writes `0` too instead of dropping the entity.
    fn emit_offset_surface(
        &mut self,
        indicator: &GpVec,
        distance: f64,
        surface: Option<usize>,
    ) -> usize {
        let params = format!(
            "140,{},{},{},{},",
            num(indicator.x() / IGES_UNIT),
            num(indicator.y() / IGES_UNIT),
            num(indicator.z() / IGES_UNIT),
            num(distance)
        );
        match surface {
            Some(idx) => self.emit_refs(140, 0, format!("{params}#0;"), &[idx]),
            None => self.emit_refs(140, 0, format!("{params}0;"), &[]),
        }
    }

    /// 'Geom2dToIGES_Geom2dCurve::Transfer2dCurve' +
    /// 'GeomToIGES_GeomCurve::TransferCurve': lift the 2-D curve to the plane
    /// Z=0 and reuse the 3-D emitter of 'emit_curve_range'.
    fn emit_curve2d_range(&mut self, curve: &dyn Curve2d, first: f64, last: f64) -> Option<usize> {
        let promoted = promote_curve2d(curve)?;
        let idx = self.emit_curve_range(promoted.as_ref(), first, last);
        if std::env::var("IGES_TRACE").is_ok() {
            eprintln!(
                "TRACE c2d de={:?} lin2d={} circ2d={} elips2d={} bsp2d={} -> kind={:?} a={first} b={last}",
                idx,
                curve.gp_lin2d().is_some(),
                curve.gp_circ2d().is_some(),
                curve.gp_elips2d().is_some(),
                curve.bspline_poles2d().is_some(),
                iges_curve_kind(promoted.as_ref()),
            );
        }
        idx
    }

    /// 'BRepToIGES_BRWire::TransferEdge(edge, face, originMap, length, false)'
    /// ('BRepToIGES_BRWire.cxx:340-585'): the edge's p-curve, corrected for the
    /// face's surface type and emitted as an IGES curve entity.
    ///
    /// 'theOriginMap' has no port equivalent - OCCT only uses it to look up the
    /// origin edge of a reversed face before 'SetShapeResult', whose transfer-map
    /// bookkeeping this port does not keep - so it is treated as empty.
    fn transfer_edge_uv(
        &mut self,
        e: &Edge,
        f: &Face,
        uv_bounds: (f64, f64, f64, f64),
        reversal: Option<&FaceReversal>,
    ) -> Option<usize> {
        // 'cxx:348-352': 'GetPCurveMode() == 0' cannot happen (the default is
        // On) and 'theIsBRepMode' is false, so only the degeneracy guard applies.
        if BRepTool::is_degenerated(e) {
            return None;
        }
        // 'cxx:374-377' reads the surface of the face it is handed, which for a
        // REVERSED face is `TransferFace`'s U-reversed copy (`cxx:135-146`).
        let surf = match reversal {
            Some(r) => r.surf.clone(),
            None => BRepTool::face_surface(f)?,
        };
        // 'cxx:374-377': only a bare 'Geom_Plane' returns here. A
        // 'Geom_RectangularTrimmedSurface' over a plane is *not* a plane for
        // OCCT's 'IsKind'; the port's trimmed surface delegates 'gp_pln', so the
        // raw type is required.
        if surf.rectangular_trimmed_basis().is_none() && surf.gp_pln().is_some() {
            return None;
        }
        let (ufirst, _ulast, vfirst, vlast) = uv_bounds; // 'BRepTools::UVBounds' ('cxx:379')
        // 'cxx:380-397': peel a rectangular trim, then (from the raw surface
        // handle) an offset surface.
        let mut surf_base: Arc<dyn Surface> = match surf.rectangular_trimmed_basis() {
            Some(b) => b,
            None => surf.clone(),
        };
        if surf.is_offset_surface() {
            if let Some(b) = surf.offset_basis_surface() {
                surf_base = b;
            }
        }
        let is_cyl = surf_base.gp_cylinder().is_some();
        let is_cone = surf_base.gp_cone().is_some();
        let is_sphere = surf_base.gp_sphere().is_some();
        let is_torus = surf_base.gp_torus().is_some();
        let is_rev = surf_base.is_surface_of_revolution();
        let is_extr = surf_base.is_surface_of_linear_extrusion();
        let is_bspline = surf_base.is_bspline_surface();

        // 'cxx:359-360': the p-curve and the range of its CurveOnSurface
        // representation.
        let (curve2d, mut first, mut last) = crate::boptools_2d::curve_on_surface_range(e, f)?;
        if let Ok(path) = std::env::var("IGES_TRACE3D") {
            let (e3a, e3b) = BRepTool::edge_parameters(e);
            let c3 = BRepTool::edge_curve(e);
            let (sa, sb) = surf_base.u_range();
            let (ta, tb) = surf_base.v_range();
            if let Some(k) = c3.as_ref().and_then(|c| c.gp_circ()) {
                use std::io::Write;
                if let Ok(mut f) = std::fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(&path)
                {
                    let _ = writeln!(
                        f,
                        "e3=({e3a:.17e},{e3b:.17e}) c3kind={:?} rad={:.17e} | uv2=({first:.17e},{last:.17e}) | u_rng=({sa:.17e},{sb:.17e}) v_rng=({ta:.17e},{tb:.17e}) up={} sph={is_sphere} rev={is_rev} loc=({:.17e},{:.17e},{:.17e}) ax=({:.17e},{:.17e},{:.17e})",
                        crate::iges::iges_curve_kind(c3.as_ref().unwrap().as_ref()),
                        k.radius(),
                        surf_base.is_u_periodic(),
                        k.position().location().x(), k.position().location().y(), k.position().location().z(),
                        k.position().direction().x(), k.position().direction().y(), k.position().direction().z(),
                    );
                }
            }
        }
        if std::env::var("IGES_TRACE").is_ok() {
            eprintln!(
                "TRACE raw pcurve lin2d={} circ2d={} elips2d={} bsp2d={} trimmed={} nbpol={} deg={:?} rng=({first},{last}) rev={} | surf pln={} cyl={is_cyl} cone={is_cone} sph={is_sphere} tor={is_torus} rev={is_rev} extr={is_extr} bsp={is_bspline} rect={} off={} kind={:?}",
                curve2d.gp_lin2d().is_some(),
                curve2d.gp_circ2d().is_some(),
                curve2d.gp_elips2d().is_some(),
                curve2d.bspline_poles2d().is_some(),
                curve2d.trimmed_basis().is_some(),
                curve2d.bspline_poles2d().map_or(0, |(x, _)| x.len()),
                curve2d.bspline_degree(),
                e.0.orientation() == crate::abs::Orientation::Reversed,
                surf.gp_pln().is_some(),
                surf.rectangular_trimmed_basis().is_some(),
                surf.is_offset_surface(),
                crate::pcurve_full::classify_surface_kind(surf_base.as_ref()),
            );
            eprintln!(
                "TRACE stored pcurve={}",
                crate::boptools_2d::curve_on_surface(e, f).is_some()
            );
            if let Ok(pcs) = std::env::var("IGES_TRACE_PC") {
                let fk = crate::tgeometry::GeometryRegistry::shape_key(&f.0);
                let list = crate::tgeometry::GeometryRegistry::global().edge_pcurves(&e.0, fk);
                let kinds: Vec<String> = list
                    .iter()
                    .map(|c| {
                        format!(
                            "lin{} circ{} elips{} bsp{} np{} trim{}",
                            c.gp_lin2d().is_some() as u8,
                            c.gp_circ2d().is_some() as u8,
                            c.gp_elips2d().is_some() as u8,
                            c.bspline_poles2d().is_some() as u8,
                            c.bspline_poles2d().map_or(0, |(x, _)| x.len()),
                            c.trimmed_basis().is_some() as u8
                        )
                    })
                    .collect();
                use std::io::Write;
                if let Ok(mut fp) = std::fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(&pcs)
                {
                    let _ = writeln!(
                        fp,
                        "PC face={:?} edge={:?} rev={} n={} [{}]",
                        Arc::as_ptr(&f.0.tshape) as usize,
                        Arc::as_ptr(&e.0.tshape) as usize,
                        e.0.orientation() == crate::abs::Orientation::Reversed,
                        list.len(),
                        kinds.join(" ")
                    );
                }
            }
        }

        // `BRepToIGES_BRShell.cxx:186-201`: on a REVERSED face every p-curve was
        // mirrored in `TransferFace` (`aCurve1->Transformed(T)`, `:200`) before
        // the edge was re-attached to the copied face, so the curve read here is
        // the mirrored one. `Transformed` keeps the curve's own parameter range
        // (`Geom2d_Transformed` forwards `FirstParameter`/`LastParameter`), so
        // `first`/`last` stay as they are.
        let curve2d: Arc<dyn Curve2d> = match reversal {
            Some(r) => {
                let axis = GpAx2d::new(
                    GpPnt2d::new(0.5 * r.two_center, vfirst),
                    GpDir2d::new(0.0, 1.0).expect("gp::DY2d is non-null"),
                );
                curve2d_transformed(&curve2d, &mirror_ax2d_trsf(&axis))
            }
            None => curve2d,
        };

        // 'cxx:403-422': 'analyticMode' is false, so the '!analyticMode' guard
        // holds; a surface of revolution whose (trim-unwrapped) basis curve is a
        // line also needs the shift.
        let mut need_shift = is_cyl || is_cone;
        if is_rev {
            if let Some(c) = surf_base.revolution_basis_curve() {
                let c: Arc<dyn Curve> = if c.is_geom_trimmed() {
                    c.untrimmed_basis().map(|(b, _, _)| b).unwrap_or(c)
                } else {
                    c
                };
                if c.is_line() {
                    need_shift = true;
                }
            }
        }
        let mut curve2d: Arc<dyn Curve2d> = if need_shift {
            // 'cxx:425-428': translate by '-Vfirst'.
            curve2d_transformed(&curve2d, &translate2d(0.0, -vfirst))
        } else {
            // 'cxx:431': 'Curve2d->Copy()'.
            Arc::from(curve2d.clone_dyn())
        };

        // 'cxx:435-455': a periodic B-spline surface is brought back into its
        // own period.
        if is_bspline {
            let (su0, su1) = surf_base.u_range();
            let (sv0, sv1) = surf_base.v_range();
            let mut u_shift = 0.0;
            let mut v_shift = 0.0;
            if surf_base.is_u_periodic() && (ufirst - su0).abs() > occt_core::precision::PCONFUSION {
                u_shift = crate::pcurve_full::adjust_to_period(ufirst, su0, su1);
            }
            if surf_base.is_v_periodic() && (vfirst - sv0).abs() > occt_core::precision::PCONFUSION {
                v_shift = crate::pcurve_full::adjust_to_period(vfirst, sv0, sv1);
            }
            if u_shift.abs() > occt_core::precision::PCONFUSION
                || v_shift.abs() > occt_core::precision::PCONFUSION
            {
                curve2d = curve2d_transformed(&curve2d, &translate2d(u_shift, v_shift));
            }
        }

        // 'cxx:457-474': the IGES surface of revolution inverts (u, v). The
        // cylinder / cone / sphere arm sits behind the '!analyticMode' guard,
        // which is true here; the revolution / torus arm has no such guard.
        if is_cyl || is_cone || is_sphere {
            curve2d = curve2d_transformed(&curve2d, &mirror_origin_dir11());
            curve2d = curve2d_transformed(&curve2d, &mirror_ox2d());
            curve2d = curve2d_transformed(
                &curve2d,
                &translate2d(0.0, 2.0 * std::f64::consts::PI),
            );
        }
        if is_rev || is_torus {
            curve2d = curve2d_transformed(&curve2d, &mirror_origin_dir11());
            curve2d = curve2d_transformed(&curve2d, &mirror_ox2d());
            curve2d = curve2d_transformed(
                &curve2d,
                &translate2d(0.0, 2.0 * std::f64::consts::PI),
            );
        }

        // UNPORTED ('cxx:476-502', all under 'analyticMode'): the cylinder /
        // cone 'myLen = PI/180' arm, the sphere / torus '180/PI' scaling and the
        // cone-apex translation are unreachable because every 'TransferEdge'
        // call from the face path passes 'theIsBRepMode = false'.
        // 'cxx:504-538': the 'theIsBRepMode && Surf->IsKind(Geom_Plane)' branch
        // is false too, so 'trans' only carries the extrusion scale.
        let my_len = surface_transfer_length(surf_base.as_ref(), vfirst, vlast);
        let mut trans = GpTrsf2d::identity();
        let mut u_fact = 1.0;
        if is_extr {
            // 'cxx:510-532' (emv, bug OCC22126): scale by '1/(Vlast-Vfirst)',
            // then 'uFact = (Vlast - Vfirst) / (us2 - us1)'.
            let _ = trans.set_scale(&GpPnt2d::zero(), 1.0 / (vlast - vfirst));
            let (us1, us2) = surf_base.u_range();
            let du = us2 - us1;
            u_fact = (vlast - vfirst) / du;
        }
        if is_cyl || is_cone || is_rev {
            // 'cxx:533-538': 'uFact = 1. / myLen'.
            u_fact = 1.0 / my_len;
        }
        curve2d = transform_pcurve(&curve2d, &trans, u_fact, &mut first, &mut last)?;

        // 'cxx:546-556': a second 'TransformPCurve' pass on a surface of linear
        // extrusion, shifting the p-curve range onto [0, 1] in u and v.
        if is_extr {
            let (us1, us2) = surf_base.u_range();
            let du = us2 - us1;
            let trans1 = translate2d(-us1 / du, -vfirst / (vlast - vfirst));
            curve2d = transform_pcurve(&curve2d, &trans1, 1.0, &mut first, &mut last)?;
        }

        // 'cxx:558-565': a REVERSED edge reverses the 2-D curve.
        if e.0.orientation() == crate::abs::Orientation::Reversed {
            let tmp_first = curve2d.reversed_parameter(last);
            let tmp_last = curve2d.reversed_parameter(first);
            curve2d = Arc::from(curve2d.reversed());
            first = tmp_first;
            last = tmp_last;
        }

        // 'cxx:566-568': 'Geom2dToIGES_Geom2dCurve::Transfer2dCurve'.
        if std::env::var("IGES_TRACE").is_ok() {
            let a = curve2d.d0(first);
            let b = curve2d.d0(last);
            eprintln!(
                "TRACE emit pcurve lin2d={} circ2d={} elips2d={} bsp2d={} np={} rng=({first},{last}) d0=({},{}) d1=({},{})",
                curve2d.gp_lin2d().is_some(),
                curve2d.gp_circ2d().is_some(),
                curve2d.gp_elips2d().is_some(),
                curve2d.bspline_poles2d().is_some(),
                curve2d.bspline_poles2d().map_or(0, |(x, _)| x.len()),
                a.x(),
                a.y(),
                b.x(),
                b.y()
            );
        }
        self.emit_curve2d_range(curve2d.as_ref(), first, last)
    }

    /// A face as `BRepToIGES_BRShell::TransferFace` builds it
    /// (`BRepToIGES_BRShell.cxx:250-405`): the base surface entity, one
    /// `CurveOnSurface` (142) per wire, and a `TrimmedSurface` (144) that carries
    /// the outer/inner contour pointers.
    ///
    /// The 3-D curve of a wire comes from `BRepToIGES_BRWire::TransferWire`
    /// (`BRepToIGES_BRWire.cxx:662-781`): a single-edge wire is that edge's curve
    /// entity, a wire with two or more edges a `CompositeCurve` (102, form 0).
    /// The 2-D (UV) curve comes from the same `TransferWire` loop (R2-19): OCCT
    /// calls `TransferEdge(edge, face, originMap, length, false)` (`:720`) per
    /// edge and collects the results exactly as it does the 3-D curves
    /// (`:755-774`). The `142` preference is 3 when both curves transferred, 2
    /// for the 3-D-only case and 1 for UV-only (`BRepToIGES_BRShell.cxx:281-292`).
    ///
    /// **Registered remainder (T-84)**: `510` still carries this port's earlier
    /// layout instead of `IGESSolid_ToolFace::WriteOwnParams`
    /// (`510, surface, nb_loops, has_outer_loop, loop_ptrs...`) over `508` Loop
    /// entities, and `514`/`186` follow the same old shape.
    fn emit_face(&mut self, f: &Face) -> Option<usize> {
        // `BRepToIGES_BRShell.cxx:257-261`: a face whose basis surface
        // transferred to a null entity is dropped (`AddWarning(start, "the basic
        // surface is a null entity"); return res;`).
        //
        // `cxx:126-237`: a REVERSED face first goes through the U-reversed
        // surface / mirrored p-curve branch; `face_reversal` reports whether that
        // branch applies (and is `None` for the orientations OCCT leaves alone).
        let reversal = face_reversal(f);
        let surf_idx = self.emit_face_surface(f, reversal.as_ref())?;
        let wires = wires_of_face(f);
        if std::env::var("IGES_TRACE_WIRE").is_ok() {
            for w in &wires {
                let raw = w.0.tshape.read().expect("poisoned TShape lock").children.clone();
                eprintln!("WIRE ori={:?} n={}", w.0.orientation(), raw.len());
                for (i, c) in raw.iter().enumerate() {
                    let e = Edge(c.clone());
                    let (v1, v2) = edge_vertices(&e);
                    let f = |v: &Option<Vertex>| -> String {
                        match v {
                            Some(v) => {
                                let q = BRepTool::vertex_point(v);
                                format!("({:.4},{:.4},{:.4})", q.x(), q.y(), q.z())
                            }
                            None => "(?)".into(),
                        }
                    };
                    eprintln!("  child{} ori={:?} {} -> {}", i + 1, c.orientation(), f(&v1), f(&v2));
                }
                for (i, e) in reordered_wire_edges(w).iter().enumerate() {
                    let (v1, v2) = edge_vertices(e);
                    let f = |v: &Option<Vertex>| -> String {
                        match v {
                            Some(v) => {
                                let q = BRepTool::vertex_point(v);
                                format!("({:.4},{:.4},{:.4})", q.x(), q.y(), q.z())
                            }
                            None => "(?)".into(),
                        }
                    };
                    eprintln!("  out{} ori={:?} {} -> {}", i + 1, e.0.orientation(), f(&v1), f(&v2));
                }
                {
                    let stored = wire_data_edges(w);
                    let mut order = crate::meshing::wire_order::WireOrder::new();
                    for e in &stored {
                        let (a, b) = (wire_edge_first_vertex(e), wire_edge_last_vertex(e));
                        if let (Some(a), Some(b)) = (a, b) {
                            let (pa, pb) = (BRepTool::vertex_point(&a), BRepTool::vertex_point(&b));
                            eprintln!("    in begin=({:.3},{:.3},{:.3}) end=({:.3},{:.3},{:.3})", pa.x(), pa.y(), pa.z(), pb.x(), pb.y(), pb.z());
                            order.add_edge_xyz(pa, pb);
                        }
                    }
                    order.perform();
                    let chain: Vec<i32> = (1..=stored.len()).map(|i| order.ordered(i)).collect();
                    eprintln!("    ORDER status={:?} chain={:?}", order.status(), chain);
                }
            }
        }
        if std::env::var("IGES_TRACE_FACE").is_ok() {
            let b = face_uv_bounds_finite(f);
            let kind = match BRepTool::face_surface(f) {
                Some(s) if s.gp_cylinder().is_some() => "cyl",
                Some(s) if s.gp_cone().is_some() => "cone",
                Some(s) if s.gp_sphere().is_some() => "sph",
                Some(s) if s.gp_torus().is_some() => "tor",
                Some(s) if s.gp_pln().is_some() => "pln",
                Some(s) if s.is_bspline_surface() => "bspl",
                Some(_) => "other",
                None => "none",
            };
            eprintln!(
                "TRACEFACE ptr={:?} de={:?} ori={:?} rev={} kind={kind} uv=({:.9},{:.9},{:.9},{:.9})",
                Arc::as_ptr(&f.0.tshape) as usize,
                surf_idx,
                f.0.orientation(),
                reversal.is_some(),
                b.0,
                b.1,
                b.2,
                b.3
            );
        }
        // `BRepTools::UVBounds(aFace, U1, U2, V1, V2)` (`BRepToIGES_BRShell.cxx:253`),
        // kept once for every `TransferEdge` call of this face (OCCT recomputes
        // them inside `TransferEdge`, `BRepToIGES_BRWire.cxx:379`).
        // A bare plane face never reaches `TransferEdge`'s UV computation
        // (`BRepToIGES_BRWire.cxx:374-377`), so the (non-trivial) bounds are only
        // derived for the faces that need them.
        let uv_bounds = match BRepTool::face_surface(f) {
            Some(s) if !(s.rectangular_trimmed_basis().is_none() && s.gp_pln().is_some()) => {
                face_uv_bounds_finite(f)
            }
            _ => (0.0, 1.0, 0.0, 1.0),
        };
        // `ShapeAlgo::AlgoContainer()->OuterWire(aFace)` (`cxx:275`) resolves to
        // `ShapeAnalysis::OuterWire` (`ShapeAnalysis.cxx`: the last wire, or the
        // first one with a non-negative `TotCross2D`).
        let outer_wire = crate::meshing::model_builder::outer_of_wires(&wires, f);
        let mut curve_refs: Vec<usize> = Vec::new();
        let mut outer_curve: Option<usize> = None;
        let mut inner_curves: Vec<usize> = Vec::new();
        // `BRepToIGES_BRShell.cxx:270`: `Iprefer` is initialised once and only
        // updated by the wires / free edges that transfer (OCCT never resets it).
        // The OCCT initial value is 0; every 142 below assigns it before reading,
        // so the initialiser is dead code but kept for fidelity.
        #[allow(unused_assignments)]
        let mut iprefer = 0i32;
        // `TransferWire` only writes `theCurve2d` when at least one edge
        // produced a UV curve (`cxx:755-774`), so a wire with none leaves the
        // caller's handle untouched (`BRepToIGES_BRShell.cxx:271`); the port
        // carries it the same way.
        #[allow(unused_assignments)]
        let mut carried_uv: Option<usize> = None;
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
            let mut uv_refs: Vec<usize> = Vec::new();
            let mut last_uv: Option<usize> = None;
            for e in reordered_wire_edges(w) {
                // `BRepToIGES_BRWire::TransferEdge` (`BRepToIGES_BRWire.cxx:266-335`)
                // transfers the edge's 3-D curve into a fresh entity on every call:
                // `Curve3d->Copy()` (`:300`) and `GeomToIGES_GeomCurve::TransferCurve`
                // (`:316`) create a new `IGESGeom_Line` / arc / 126 each time. OCCT
                // shares nothing between the wires of a face, the wires of different
                // faces, or the two occurrences of a seam edge, so the port must not
                // deduplicate here either (Cube's 6 faces x 4 edges = 24 lines).
                let idx = self.emit_edge_curve(&e);
                if let Some(idx) = idx {
                    edge_refs.push(idx);
                    curve_refs.push(idx);
                }
                // `BRepToIGES_BRWire::TransferWire` (`cxx:715-724`): every edge's
                // 2-D curve, collected in the same order as the 3-D ones.
                let uv = self.transfer_edge_uv(&e, f, uv_bounds, reversal.as_ref());
                last_uv = uv;
                if let Some(u) = uv {
                    uv_refs.push(u);
                }
            }
            if edge_refs.is_empty() && uv_refs.is_empty() {
                continue;
            }
            let curve3d = if edge_refs.is_empty() {
                None
            } else if edge_refs.len() == 1 {
                Some(edge_refs[0])
            } else {
                Some(self.emit_refs(102, 0, format!("102,{},{};", edge_refs.len(), ph_list(edge_refs.len())), &edge_refs))
            };
            // `cxx:755-774`: one 2-D entity is used directly, two or more become
            // a `CompositeCurve` (102). The single-entity arm mirrors
            // `theCurve2d = ent2d` (`cxx:760`), i.e. the last edge's transfer
            // result rather than `Seq2d(1)`; zero 2-D entities leave the caller's
            // handle from the previous wire in place.
            let curve_uv = if uv_refs.len() == 1 {
                last_uv
            } else if uv_refs.len() >= 2 {
                Some(self.emit_refs(102, 0, format!("102,{},{};", uv_refs.len(), ph_list(uv_refs.len())), &uv_refs))
            } else {
                carried_uv
            };
            carried_uv = curve_uv;
            // `BRepToIGES_BRShell.cxx:281-292`: the preference follows which of
            // the two representations transferred; both null leaves it at its
            // previous value.
            if curve3d.is_some() {
                iprefer = if curve_uv.is_some() { 3 } else { 2 };
            } else if curve_uv.is_some() {
                iprefer = 1;
            }
            // `IGESGeom_CurveOnSurface::Init` (`IGESGeom_CurveOnSurface.cxx:26-40`)
            // with `Imode = 0` (`cxx:269`);
            // `IGESGeom_ToolCurveOnSurface::WriteOwnParams` (`:144-152`) writes
            // `142, creation_mode, surface, curve_uv, curve_3d, preference_mode;`
            // - field 3 is the UV curve, field 4 the 3-D one.
            let cs = match (curve_uv, curve3d) {
                (Some(uv), Some(c3)) => {
                    self.emit_refs(142, 0, format!("142,0,#0,#1,#2,{iprefer};"), &[surf_idx, uv, c3])
                }
                (Some(uv), None) => {
                    self.emit_refs(142, 0, format!("142,0,#0,#1,0,{iprefer};"), &[surf_idx, uv])
                }
                (None, Some(c3)) => {
                    self.emit_refs(142, 0, format!("142,0,#0,0,#1,{iprefer};"), &[surf_idx, c3])
                }
                (None, None) => {
                    self.emit_refs(142, 0, format!("142,0,#0,0,0,{iprefer};"), &[surf_idx])
                }
            };
            match &outer_wire {
                Some(o) if Arc::as_ptr(&o.0.tshape) == Arc::as_ptr(&w.0.tshape) => {
                    outer_curve = Some(cs)
                }
                _ => inner_curves.push(cs),
            }
        }
        // `cxx:334-365`: edges of the face that are not part of any wire become
        // further inner contours. OCCT transfers both their 3-D curve (`:344`) and
        // their UV curve (`TransferEdge(edge, face, originMap, length, false)`,
        // `:346`).
        let wire_edge_keys: std::collections::HashSet<usize> = wires
            .iter()
            .flat_map(|w| edges_of_wire(w))
            .map(|e| Arc::as_ptr(&e.0.tshape) as usize)
            .collect();
        for e in children_of_type(&f.0, ShapeType::Edge) {
            let key = Arc::as_ptr(&e.tshape) as usize;
            if wire_edge_keys.contains(&key) {
                continue;
            }
            let edge = Edge(e);
            // Each free edge is its own `TransferEdge` result too
            // (`BRepToIGES_BRShell.cxx:344-346`), with no sharing.
            let idx = self.emit_edge_curve(&edge);
            if let Some(idx) = idx {
                curve_refs.push(idx);
            }
            let uv = self.transfer_edge_uv(&edge, f, uv_bounds, reversal.as_ref());
            // `BRepToIGES_BRShell.cxx:347-358`: the preference follows which of
            // the two representations transferred.
            if idx.is_some() {
                iprefer = if uv.is_some() { 3 } else { 2 };
            } else if uv.is_some() {
                iprefer = 1;
            }
            let cs = match (uv, idx) {
                (Some(u), Some(c3)) => {
                    self.emit_refs(142, 0, format!("142,0,#0,#1,#2,{iprefer};"), &[surf_idx, u, c3])
                }
                (Some(u), None) => {
                    self.emit_refs(142, 0, format!("142,0,#0,#1,0,{iprefer};"), &[surf_idx, u])
                }
                (None, Some(c3)) => {
                    self.emit_refs(142, 0, format!("142,0,#0,0,#1,{iprefer};"), &[surf_idx, c3])
                }
                (None, None) => {
                    self.emit_refs(142, 0, format!("142,0,#0,0,0,{iprefer};"), &[surf_idx])
                }
            };
            inner_curves.push(cs);
        }

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
        // A face this port built without boundary wires carries no contour on
        // this side either; OCCT's `TransferWire` leaves the caller's handle
        // untouched when no edge transfers (`BRepToIGES_BRWire.cxx:715-774`).

        if curve_refs.is_empty() {
            return Some(surf_idx);
        }
        let outer_flag = u8::from(!is_whole);
        let n_inner = inner_curves.len();
        // `IGESGeom_ToolTrimmedSurface::WriteOwnParams`
        // (`IGESGeom_ToolTrimmedSurface.cxx:196-218`): `144, surface,
        // outer_boundary_type, nb_inner_contours, outer_contour, inner...;`
        // with the outer contour written as `0` when the type is false.
        // `surface` and the contours are pointers (`#k`), the boundary type and
        // the inner-contour count are plain integers.
        let has_outer = !is_whole && outer_curve.is_some();
        let mut ref_list = vec![surf_idx];
        if has_outer {
            ref_list.push(outer_curve.unwrap());
        }
        ref_list.extend(inner_curves.iter().copied());
        let mut contour_refs = if has_outer { ph(1) } else { "0".to_string() };
        for k in if has_outer { 2 } else { 1 }..ref_list.len() {
            contour_refs.push_str(&format!(",{}", ph(k)));
        }
        Some(self.emit_refs(
            144,
            0,
            format!("144,#0,{outer_flag},{n_inner},{contour_refs};"),
            &ref_list,
        ))
    }

    /// `Group` (402) when there is more than one item, the item itself otherwise -
    /// `IGESBasic_ToolGroup::WriteOwnParams` (`IGESBasic_ToolGroup.cxx:104-115`)
    /// writes `402, n, entity...;` and both `BRepToIGES_BRShell::TransferShell`
    /// (`BRepToIGES_BRShell.cxx:463-471`) and
    /// `BRepToIGES_BRSolid::TransferSolid` (`BRepToIGES_BRSolid.cxx:154-163`) use
    /// exactly this rule. The form is 1: `IGESBasic_Group::IGESBasic_Group`
    /// calls `InitTypeAndForm(402, 1)` (`IGESBasic_Group.cxx:29-31`), and
    /// `DirChecker` requires that form (`IGESBasic_ToolGroup.cxx:161`).
    fn group_or_single(&mut self, items: Vec<usize>) -> Option<usize> {
        match items.len() {
            0 => None,
            1 => Some(items[0]),
            n => {
                Some(self.emit_refs(402, 1, format!("402,{n},{};", ph_list(n)), &items))
            }
        }
    }

    /// `BRepToIGES_BRShell::TransferShell` (`BRepToIGES_BRShell.cxx:411-476`): the
    /// shell's faces are transferred as trimmed surfaces (144) and grouped.
    fn emit_shell(&mut self, sh: &Shell) -> Option<usize> {
        if std::env::var("IGES_TRACE_FACE").is_ok() {
            eprintln!("TRACESHELL ori={:?}", sh.0.orientation());
        }
        // `TopExp_Explorer(shell, FACE)` composes the shell orientation
        // (`BRepToIGES_BRShell.cxx:433`). A shell `SolidFromShell` reversed
        // (`ShapeFix_Solid.cxx:686`) therefore flips every face.
        let faces: Vec<TopoShape> = crate::iterator::cumulated_children(&sh.0)
            .into_iter()
            .filter(|h| h.shape_type() == ShapeType::Face)
            .collect();
        let refs: Vec<usize> = faces
            .iter()
            .filter_map(|f| self.emit_face(&Face(f.clone())))
            .collect();
        self.group_or_single(refs)
    }

    /// `BRepToIGES_BRSolid::TransferSolid` (`BRepToIGES_BRSolid.cxx:100-168`): the
    /// solid's shells are transferred and grouped the same way.
    fn emit_solid(&mut self, s: &Solid) -> Option<usize> {
        // `TopExp_Explorer(solid, SHELL)` (`BRepToIGES_BRSolid.cxx:124`).
        let shells: Vec<TopoShape> = crate::iterator::cumulated_children(&s.0)
            .into_iter()
            .filter(|h| h.shape_type() == ShapeType::Shell)
            .collect();
        let refs: Vec<usize> = shells
            .iter()
            .filter_map(|sh| self.emit_shell(&Shell(sh.clone())))
            .collect();
        self.group_or_single(refs)
    }

    fn note_max_coord(&mut self, shape: &TopoShape) {
        // `BRepBndLib::Add` then skip a void or open box
        // (`IGESControl_Writer.cxx:181-189`).
        let box_ = crate::bopalgo_tools_wires::shape_box_of(shape);
        if box_.is_void()
            || box_.is_open_xmin()
            || box_.is_open_xmax()
            || box_.is_open_ymin()
            || box_.is_open_ymax()
            || box_.is_open_zmin()
            || box_.is_open_zmax()
        {
            return;
        }
        let Some((xmin, xmax, ymin, ymax, zmin, zmax)) = box_.get() else {
            return;
        };
        for v in [xmin, xmax, ymin, ymax, zmin, zmax] {
            let a = v.abs();
            if a.is_finite() && a > self.max_coord {
                self.max_coord = a;
            }
        }
    }

    fn emit_shape(&mut self, shape: &TopoShape) {
        self.note_max_coord(shape);
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
                if let Some(r) = self.emit_solid(&Solid(shape.clone())) {
                    self.roots.push(r);
                }
            }
            ShapeType::Shell => {
                if let Some(r) = self.emit_shell(&Shell(shape.clone())) {
                    self.roots.push(r);
                }
            }
            ShapeType::Face => {
                // `BRepToIGES_BRShell.cxx:257-261`: a face with a null basis
                // surface entity contributes nothing to the model.
                if let Some(r) = self.emit_face(&Face(shape.clone())) {
                    self.roots.push(r);
                }
            }
            _ => {}
        }
    }

    // ---- file assembly ----

    fn global_lines(&self) -> Vec<String> {
        // `IGESData::Init` (`IGESData.cxx:249-279`) plus the per-shape update
        // in `IGESControl_Writer::AddShape` (`IGESControl_Writer.cxx:177-189`).
        // Default comma and semicolon are written as two empty fields
        // (`IGESData_GlobalSection::Params`, `cxx:445-462`), not as `1H,`.
        // Single and double precision are both `RealLast10Exp()` / `RealDigits()`
        // (`cxx:258-261`), which is `308` / `15`, and the version is
        // `SetIGESVersion(11)` (`cxx:276`). The date is `YYYYMMDD.HHMMSS`
        // (`NewDateString`, `cxx:247`); its Hollerith count has to be the
        // character length or every later field, including the version, shifts.
        let date = iges_date_now();
        let product = "Open CASCADE IGES processor 8.0";
        let system = "Open CASCADE 8.0";
        let max_coord = if self.max_coord > 0.0 {
            num(self.max_coord)
        } else {
            String::new()
        };
        let data = format!(
            ",,{product_h},{file_h},{system_h},{product_h},32,308,15,308,15,,{scale},2,{unit},1,{weight},{date_h},{resolution},{max_coord},,,11,0,{date_h},;",
            product_h = hollerith(product),
            file_h = hollerith("Filename.iges"),
            system_h = hollerith(system),
            scale = num(1.0),
            unit = hollerith("MM"),
            weight = num(0.01),
            date_h = hollerith(&date),
            resolution = num(occt_core::CONFUSION),
        );
        let mut out = Vec::new();
        let mut s = data.as_str();
        while !s.is_empty() {
            let end = s.len().min(MAXCARS_G);
            out.push(s[..end].to_string());
            s = &s[end..];
        }
        out
    }

    fn finish(self) -> String {
        self.check_refs();
        // T-85 step 2: OCCT's model holds only the entities reachable from the
        // shapes handed to `AddEntity` (`AddWithRefs`), numbered in the order it
        // adds them; the placeholders in the parameter text become the final
        // numbers here.
        let mut entities = self.final_entities();
        // `IGESControl_Writer::ComputeModel` (`IGESControl_Writer.cxx:256-262`).
        compute_status(&mut entities);
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
        for (i, ent) in entities.iter().enumerate() {
            let (l1, l2) = ent.directory_lines(p_start);
            p_start += ent.param_line_count();
            out.push_str(&sec_line('D', 2 * i + 1, &l1, MAXCARS_G));
            out.push('\n');
            out.push_str(&sec_line('D', 2 * i + 2, &l2, MAXCARS_G));
            out.push('\n');
            d_seq += 2;
        }

        for (i, ent) in entities.iter().enumerate() {
            for chunk in ent.param_chunks() {
                out.push_str(&param_line(&chunk, 2 * i + 1, p_seq));
                out.push('\n');
                p_seq += 1;
            }
        }

        // Terminate card (`IGESData_IGESWriter.cxx:942-947`): the **last sequence
        // number** of each section - `nbs`, `nbg`, `nbd * 2` and
        // `thepnum.Value(thepnum.Length()) - 1` - blank-filled, then 40 blanks and
        // `T0000001`. The Directory section therefore reports `2 * entities`, not
        // `2 * entities + 1` (as of batch 70).
        let t = format!(
            "S{:>7}G{:>7}D{:>7}P{:>7}{}T0000001",
            s_seq - 1,
            g_seq - 1,
            d_seq - 1,
            p_seq - 1,
            " ".repeat(40)
        );
        out.push_str(&t);
        out.push('\n');
        out
    }
}

/// The `ShapeExtend_WireData` edge list of a wire (`ShapeExtend_WireData.cxx:72-148`
/// with the default `chained = true`, `theManifold = true`): `TopoDS_Iterator`'s
/// composed orientation, and — when the wire itself is REVERSED — every edge
/// `Prepend`ed instead of appended (`cxx:114-121`), so the list comes out in the
/// wire's own traversal order. `INTERNAL`/`EXTERNAL` edges go to
/// `myNonmanifoldEdges` (`cxx:84-89`), which `TransferWire` never reads.
fn wire_data_edges(wire: &crate::shape::Wire) -> Vec<Edge> {
    let mut stored: Vec<Edge> = edges_of_wire(wire)
        .into_iter()
        .filter(|e| {
            let o = e.0.orientation();
            o == crate::abs::Orientation::Forward || o == crate::abs::Orientation::Reversed
        })
        .collect();
    if wire.0.orientation() == crate::abs::Orientation::Reversed {
        stored.reverse();
    }
    stored
}

/// `ShapeAnalysis_Edge::FirstVertex(edge, CumOri = true)`
/// (`ShapeAnalysis_Edge.cxx:230-238`): the FORWARD vertex of the edge *as
/// traversed*, i.e. the stored REVERSED one on a REVERSED edge.
fn wire_edge_first_vertex(e: &Edge) -> Option<Vertex> {
    let (f, l) = edge_vertices(e);
    if e.0.orientation() == crate::abs::Orientation::Reversed {
        l
    } else {
        f
    }
}

/// `ShapeAnalysis_Edge::LastVertex(edge, CumOri = true)`.
fn wire_edge_last_vertex(e: &Edge) -> Option<Vertex> {
    let (f, l) = edge_vertices(e);
    if e.0.orientation() == crate::abs::Orientation::Reversed {
        f
    } else {
        l
    }
}

/// `BRepToIGES_BRWire::TransferWire` (`BRepToIGES_BRWire.cxx:704-726`) orders the
/// edges of a wire through `ShapeFix_Wire::FixReorder()` before transferring them
/// (`cxx:706`). That is the 3D `ShapeAnalysis_Wire::CheckOrder(sawo, myClosedMode,
/// true, false)` of `ShapeFix_Wire::FixReorder()` (`ShapeFix_Wire.cxx:487-534`)
/// followed by `FixReorder(sawo)` (`cxx:1351-1400`). `Perform` ignores its
/// `closed` argument (`ShapeAnalysis_WireOrder.cxx:205`), so only the edge
/// connectivity matters.
///
/// Returns the same list as `wire_data_edges` when the order is unchanged
/// (`status == 0`, `cxx:1360-1363`), when an edge has no vertex (FAIL2,
/// `ShapeAnalysis_Wire.cxx:617-621`) or when the order could not be applied
/// (`cxx:1372-1385`). `Ordered(i) < 0` is `ShapeExtend_WireData::Edge(signed)`
/// = the edge reversed (`ShapeExtend_WireData.cxx:583-591`).
fn reordered_wire_edges(wire: &crate::shape::Wire) -> Vec<Edge> {
    let stored = wire_data_edges(wire);
    if stored.len() < 2 {
        return stored;
    }
    let mut order = crate::meshing::wire_order::WireOrder::new();
    for e in &stored {
        let (Some(v1), Some(v2)) = (wire_edge_first_vertex(e), wire_edge_last_vertex(e)) else {
            return stored;
        };
        order.add_edge_xyz(BRepTool::vertex_point(&v1), BRepTool::vertex_point(&v2));
    }
    order.perform();
    if order.status() == crate::meshing::wire_order::WireOrderStatus::Same
        || order.nb_edges() != stored.len()
    {
        return stored;
    }
    let mut out = Vec::with_capacity(stored.len());
    for i in 1..=stored.len() {
        let signed = order.ordered(i);
        if signed == 0 {
            return stored;
        }
        let mut e = stored[signed.unsigned_abs() as usize - 1].clone();
        if signed < 0 {
            e.0.reverse();
        }
        out.push(e);
    }
    out
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

/// The UV box handed to `GeomToIGES_GeomSurface::TransferSurface`.
///
/// OCCT takes it straight from the face (`BRep_Tool::UVBounds`, i.e.
/// `BRepTools::AddUVBounds` -> the edge p-curve boxes, see
/// [`crate::brep_tools::add_uv_bounds`]). Measured on the 14-model parity set
/// this is strictly better than the sampled boundary box
/// (`wireframe::face_uv_bounds`, which unwraps a periodic span into
/// `[first, first + period]`): GRAND TOTAL 1396 -> 1308, `Cone` to `SAME`,
/// `screw` 144 -> 130, `Shape-1` 654 -> 634, `Shape-2` 157 -> 105, no model
/// worse. The sampled box stays as the fallback for a void or non-finite
/// result, which `BRepTools::UVBounds` reports as `0,0,0,0`
/// (`BRepTools.cxx:70-80`) and the revolution arms here cannot use.
fn face_uv_bounds_finite(f: &Face) -> (f64, f64, f64, f64) {
    if let Some(s) = BRepTool::face_surface(f) {
        let mut b = occt_core::bnd::BndBox2d::new();
        crate::brep_tools::add_uv_bounds(f, &mut b);
        if let Some((u0, v0, u1, v1)) = b.get() {
            if u0.is_finite() && u1.is_finite() && v0.is_finite() && v1.is_finite() {
                return (u0, u1, v0, v1);
            }
        }
        crate::wireframe::face_uv_bounds(f, s.as_ref())
    } else {
        (0.0, 1.0, 0.0, 1.0)
    }
}

/// `Geom_ElementarySurface::Bounds` (`Geom_ElementarySurface.cxx:82-89`): every
/// elementary surface reports the same "infinite" range, `U [0, 2*PI],
/// V [0, 2*PI]`. `TransferFace` mirrors the p-curves of a REVERSED face about
/// `U = 0.5*(U1+U2)` of the *reversed* surface's own bounds
/// (`BRepToIGES_BRShell.cxx:176-183`), so this is `2*aCenter` for every surface
/// class that has an `UReversed` branch.
const ELEMENTARY_REVERSED_TWO_CENTER: f64 = 2.0 * std::f64::consts::PI;

/// `BRepToIGES_BRShell::TransferFace` (`BRepToIGES_BRShell.cxx:126-237`): a face
/// whose orientation is `TopAbs_REVERSED` is not transferred as it is. OCCT
/// copies the face, replaces its surface with `aSurf->UReversed()` after peeling
/// every `Geom_RectangularTrimmedSurface` (`:135-141`) and mirrors each p-curve
/// about `gp_Ax2d(gp_Pnt2d(0.5*(U1+U2), V1), gp_Dir2d(0,1))` (`:176-201`), then
/// transfers the copy.
///
/// The port transfers the original face (the copied wires only get re-oriented
/// to the originals, `:150-174`) and mirrors the p-curves it reads on the fly
/// (`transfer_edge_uv`), so this carries what both call sites need: the
/// U-reversed surface and the mirror's `2*aCenter`.
///
/// `None` means the REVERSED branch is not taken: either the face is not
/// REVERSED, or its surface class has no `Geom_Surface::UReversed` branch yet
/// (elementary surfaces, a non-periodic-U B-spline, a surface of revolution or
/// linear extrusion, and an offset whose basis can be U-reversed — see
/// `Surface::u_reversed`), in which case the face is transferred as if it were
/// FORWARD, exactly as an unimplemented OCCT branch would leave it.
struct FaceReversal {
    /// `aSurf->UReversed()` (`cxx:135-142`), already trim-stripped.
    surf: Arc<dyn Surface>,
    /// `2*aCenter` of `gp_Ax2d(gp_Pnt2d(aCenter, V1), gp_Dir2d(0,1))`: the mirror
    /// maps `U` to `2*aCenter - U` (`cxx:178-183`).
    two_center: f64,
}

fn face_reversal(f: &Face) -> Option<FaceReversal> {
    if f.0.orientation() != crate::abs::Orientation::Reversed {
        return None;
    }
    let surf = BRepTool::face_surface(f)?;
    // `cxx:135-141`: peel the rectangular trims first, "because pcurves will be
    // transformed, so trim will be shifted, accorded to new face bounds".
    let mut base: Arc<dyn Surface> = surf.clone();
    while let Some(b) = base.rectangular_trimmed_basis() {
        base = b;
    }
    let rev = base.u_reversed()?;
    // `cxx:176-182`: `aSurf->Bounds(U1, U2, V1, V2); aCenter = 0.5*(U1+U2)` on the
    // REVERSED surface, so the mirror is `U -> 2*aCenter - U`. Finite bounds
    // (B-spline, revolution, offset of those) use that sum; an unbounded
    // elementary plane keeps the `[0, 2*PI]` constant the other elementary
    // surfaces report (`Geom_ElementarySurface::Bounds`).
    let (ru0, ru1) = rev.u_range();
    let two_center = if ru0.is_finite() && ru1.is_finite() {
        ru0 + ru1
    } else {
        ELEMENTARY_REVERSED_TWO_CENTER
    };
    if std::env::var("IGES_TRACE_FACE").is_ok() {
        eprintln!(
            "TRACEREV ptr={:?} surf={:?} two_center={two_center}",
            Arc::as_ptr(&f.0.tshape) as usize,
            crate::pcurve_full::classify_surface_kind(rev.as_ref()),
        );
    }
    Some(FaceReversal {
        surf: rev,
        two_center,
    })
}

/// 'GeomToIGES_GeomSurface::Length()' ('GeomToIGES_GeomSurface.cxx:1441-1443'):
/// the 'TheLength' the matching 'TransferSurface' overload stored for the
/// face's surface. Only two arms are not 1:
///
/// * 'TransferSurface(Geom_ConicalSurface)' ('cxx:777-823') sets it to the
///   generatrix span 'gen1.Distance(gen2)', and the generatrix is a unit
///   direction, so the value is '|V2 - V1|' over the (infinite-substituted) V
///   bounds;
/// * 'TransferSurface(Geom_SurfaceOfRevolution)' ('cxx:1111-1169') sets it the
///   same way only when the trim-unwrapped basis curve is a line.
///
/// Every other overload ('cxx:618', '714', '879', '947', '1047', '1132', '1253')
/// leaves it at 1.
fn surface_transfer_length(surf: &dyn Surface, v0: f64, v1: f64) -> f64 {
    let bound = |v: f64, positive: bool| -> f64 {
        if positive {
            if occt_core::precision::Precision::is_positive_infinite(v) {
                occt_core::precision::INFINITE
            } else {
                v
            }
        } else if occt_core::precision::Precision::is_negative_infinite(v) {
            -occt_core::precision::INFINITE
        } else {
            v
        }
    };
    let (vv0, vv1) = (bound(v0, false), bound(v1, true));
    // 'cxx:737-741' (GeomToIGES_GeomSurface::TransferSurface(Geom_Cylinder)):
    // the generatrix is 'Geom_Line(gp_Pnt(R,0,0), gp::DZ())', so
    // 'TheLength = gen1.Distance(gen2) = |V2 - V1|'.
    if surf.gp_cylinder().is_some() {
        return (vv1 - vv0).abs();
    }
    // 'cxx:817-825' (GeomToIGES_GeomSurface::TransferSurface(Geom_ConicalSurface)):
    // the generatrix is 'Geom_Line(gp_Pnt(RefRadius,0,0), gp_Dir(sin a,0,cos a))',
    // a unit direction, so the same expression holds.
    if surf.gp_cone().is_some() {
        return (vv1 - vv0).abs();
    }
    if surf.is_surface_of_revolution() {
        if let Some(c) = surf.revolution_basis_curve() {
            let c: Arc<dyn Curve> = if c.is_geom_trimmed() {
                c.untrimmed_basis().map(|(b, _, _)| b).unwrap_or(c)
            } else {
                c
            };
            if c.is_line() {
                return (vv1 - vv0).abs();
            }
        }
    }
    1.0
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
        // The default `write.iges.brep.mode = 0` (Faces mode) writes the faces as
        // trimmed surfaces (144 over 142) grouped by 402 - not the 510/514/186
        // BRep-mode tree this port used to emit. Point entities (116) appear only
        // where OCCT references them (the location point of 192/194/196/198); the
        // port's earlier unconditional per-vertex 116 block was unreferenced output
        // and was removed in batch 71, so a box legitimately has none.
        for needle in ["110", "108", "142", "144", "402"] {
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
        assert!(iges.contains("144"), "missing 144");
        for line in iges.lines() {
            assert_eq!(line.len(), 80, "line length {}", line.len());
        }
    }

    #[test]
    fn iges_file_written() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let mut model = BRepModel::new();
        model.add("Box", b.solid.0.clone());
        let path = { let d = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("target").join("test-out"); let _ = std::fs::create_dir_all(&d); d }.join("occt_iges_test.igs");
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
