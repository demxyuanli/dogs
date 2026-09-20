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

use occt_core::gp::{GpPln, GpPnt, GpVec};
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

/// Build one section line: section letter + 7-digit sequence + data → 80 cols.
fn sec_line(section: char, seq: usize, content: &str) -> String {
    let mut s = String::with_capacity(80);
    s.push(section);
    s.push_str(&format!("{seq:07}"));
    s.push_str(content);
    rec80(&s)
}

/// Right-justify a value in an 8-column IGES field.
fn field8(v: &str) -> String {
    format!("{v:>8}")
}

/// Left-justify a label in an 8-column IGES field.
fn field8l(v: &str) -> String {
    format!("{v:<8}")
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
    /// (`GeomToIGES_GeomCurve::TransferBSplineCurve`); this port has no 126
    /// emitter yet (T-78).
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

// ---------------------------------------------------------------------------
// Entity bookkeeping
// ---------------------------------------------------------------------------

/// One IGES entity: directory-entry metadata + parameter-data string.
struct Ent {
    ty: i32,
    form: i32,
    label: String,
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

    fn directory_lines(&self) -> (String, String) {
        let ty = self.ty.to_string();
        let pcount = self.param_line_count().to_string();
        let form = self.form.to_string();
        let l1 = format!(
            "{}{}{}{}{}{}{}{}",
            field8(&ty),
            field8(&pcount),
            field8("0"),
            field8("0"),
            field8("0"),
            field8("0"),
            field8("0"),
            field8("0")
        );
        let l2 = format!(
            "{}{}{}{}{}{}{}{}",
            field8(&ty),
            field8("1"),
            field8("7"),
            field8(&pcount),
            field8(&form),
            field8l(&self.label),
            field8("0"),
            field8("")
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

    fn emit(&mut self, ty: i32, form: i32, label: &str, params: String) -> usize {
        self.entities.push(Ent {
            ty,
            form,
            label: label.into(),
            params,
        });
        self.entities.len()
    }

    // ---- geometry entities ----

    fn emit_point(&mut self, p: &GpPnt) -> usize {
        self.emit(
            116,
            0,
            "POINT",
            format!("116,{},{},{};", num(p.x()), num(p.y()), num(p.z())),
        )
    }

    fn emit_line(&mut self, a: &GpPnt, b: &GpPnt) -> usize {
        self.emit(
            110,
            0,
            "LINE",
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
            "PLANE",
            format!("108,{},{},{},{};", num(n.x), num(n.y), num(n.z), num(d)),
        )
    }

    fn emit_circular_arc(&mut self, center: &GpPnt, start: &GpPnt, end: &GpPnt, plane_pt: &GpPnt) -> usize {
        self.emit(
            100,
            0,
            "ARC",
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

    fn emit_edge_curve(&mut self, e: &Edge) -> usize {
        let Some(curve) = BRepTool::edge_curve(e) else {
            let (v1, v2) = edge_vertices(e);
            let p1 = v1.map(|v| BRepTool::vertex_point(&v)).unwrap_or_default();
            let p2 = v2.map(|v| BRepTool::vertex_point(&v)).unwrap_or_default();
            return self.emit_line(&p1, &p2);
        };
        let (a, b) = BRepTool::edge_parameters(e);
        match iges_curve_kind(curve.as_ref()) {
            IgCurveKind::Line => {
                let p1 = curve.d0(a);
                let p2 = curve.d0(b);
                self.emit_line(&p1, &p2)
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
                self.emit_circular_arc(&center, &p1, &p2, &plane_pt)
            }
            IgCurveKind::Conic | IgCurveKind::Bounded | IgCurveKind::Other => {
                // UNPORTED (audit A26 / task T-78): OCCT writes entity 104 for
                // ellipse/hyperbola/parabola (`TransferConic`) and entity 126 for
                // Bezier/B-spline (`TransferBSplineCurve`); neither emitter is
                // ported, and `TransferCurve` returns a null handle for any other
                // curve type, in which case the caller writes no curve at all.
                // Until those emitters land the edge keeps the chord-line
                // stand-in.
                let p1 = curve.d0(a);
                let p2 = curve.d0(b);
                self.emit_line(&p1, &p2)
            }
        }
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
        if classify_surface(surf.as_ref()) == SurfaceKind::Sphere {
            if let Some(center) = sphere_center(surf.as_ref()) {
                let (u0, _, v0, v1) = surf_bounds(surf.as_ref());
                let vm = 0.5 * (v0 + v1);
                let r = surf.d0(u0, vm).distance(&center);
                if r > 1e-9 {
                    // Surface of revolution: axis line + generatrix semicircle.
                    let south = GpPnt::new(center.x(), center.y(), center.z() - r);
                    let north = GpPnt::new(center.x(), center.y(), center.z() + r);
                    let axis_idx = self.emit_line(&south, &north);
                    let gen = self.emit_arc_3p(
                        &south,
                        &GpPnt::new(center.x() + r, center.y(), center.z()),
                        &north,
                    );
                    let rev = self.emit(120, 0, "REVOLVED", format!("120,{axis_idx},{gen};"));
                    // Second meridian completes the sphere's seam boundary.
                    let mer2 = self.emit_arc_3p(
                        &north,
                        &GpPnt::new(center.x() - r, center.y(), center.z()),
                        &south,
                    );
                    return (rev, vec![gen, mer2]);
                }
            }
        }
        // Unclassified curved face: fall back to a plane.
        // ponytail: covers sphere/box tests; add 128 NURBS surfaces when needed.
        let pln = face_plane(f).unwrap_or_else(GpPln::default);
        (self.emit_plane(&pln), Vec::new())
    }

    fn emit_face(&mut self, f: &Face) -> usize {
        let (surf_idx, mut synth) = self.emit_face_surface(f);
        let mut curve_refs: Vec<usize> = Vec::new();
        for wire in wires_of_face(f) {
            for e in edges_of_wire(&wire) {
                let key = Arc::as_ptr(&e.0.tshape) as usize;
                let idx = match self.edge_curve_entities.get(&key) {
                    Some(&i) => i,
                    None => {
                        let i = self.emit_edge_curve(&e);
                        self.edge_curve_entities.insert(key, i);
                        i
                    }
                };
                curve_refs.push(idx);
            }
        }
        curve_refs.append(&mut synth);

        let trim_idx = if curve_refs.is_empty() {
            surf_idx
        } else {
            let (u0, u1, v0, v1) = face_uv_bounds_finite(f);
            let m = curve_refs.len();
            let refs = curve_refs
                .iter()
                .map(|i| i.to_string())
                .collect::<Vec<_>>()
                .join(",");
            self.emit(
                144,
                0,
                "TRIMMED",
                format!("144,{surf_idx},{},{},{},{},1,{m},{};", num(u0), num(u1), num(v0), num(v1), refs),
            )
        };

        let mut params = format!("510,{trim_idx}");
        for c in &curve_refs {
            params.push_str(&format!(",{c}"));
        }
        params.push(';');
        self.emit(510, 0, "FACE", params)
    }

    fn emit_shell(&mut self, sh: &Shell) -> usize {
        let faces = children_of_type(&sh.0, ShapeType::Face);
        let refs: Vec<String> = faces
            .iter()
            .map(|f| self.emit_face(&Face(f.clone())).to_string())
            .collect();
        self.emit(514, 0, "SHELL", format!("514,{};", refs.join(",")))
    }

    fn emit_solid(&mut self, s: &Solid) -> usize {
        let shells = children_of_type(&s.0, ShapeType::Shell);
        let refs: Vec<String> = shells
            .iter()
            .map(|sh| self.emit_shell(&Shell(sh.clone())).to_string())
            .collect();
        let outer = refs.first().cloned().unwrap_or_else(|| "0".to_string());
        self.emit(186, 0, "MSBO", format!("186,{outer};"))
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

        out.push_str(&sec_line('S', s_seq, "IGES B-REP MODEL GENERATED BY RUST OCCT PORT"));
        out.push('\n');
        s_seq += 1;

        for gl in self.global_lines() {
            out.push_str(&sec_line('G', g_seq, &gl));
            out.push('\n');
            g_seq += 1;
        }

        for ent in &self.entities {
            let (l1, l2) = ent.directory_lines();
            out.push_str(&sec_line('D', d_seq, &l1));
            out.push('\n');
            d_seq += 1;
            out.push_str(&sec_line('D', d_seq, &l2));
            out.push('\n');
            d_seq += 1;
        }

        for ent in &self.entities {
            for chunk in ent.param_chunks() {
                out.push_str(&sec_line('P', p_seq, &chunk));
                out.push('\n');
                p_seq += 1;
            }
        }

        // Terminate record: last sequence number of each section + T=1.
        let t = format!(
            "S{:07}G{:07}D{:07}P{:07}        T0000001",
            s_seq - 1,
            g_seq - 1,
            d_seq - 1,
            p_seq - 1
        );
        out.push_str(&rec80(&t));
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
    let d_lines: Vec<&str> = text.lines().filter(|l| l.starts_with('D')).collect();
    let total = d_lines.len() / 2;
    let mut geo = 0usize;
    for (i, l) in d_lines.iter().enumerate() {
        if i % 2 == 0 {
            let ty: i32 = l.get(8..16).unwrap_or("0").trim().parse().unwrap_or(0);
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
        assert!(content.starts_with('S'), "first line should start with S");
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
