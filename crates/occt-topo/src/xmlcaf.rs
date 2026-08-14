//! XML XCAF container — an XML 1.0 assembly/attributes document.
//! Source: `XmlXCAF` (XML `XCAFDoc_*` document driver).
//!
//! A `XmlXcafDoc` is a tree of entries (`XmlEntry`), each carrying an
//! optional `TopoShape`, a list of string attributes (name / color / layer /
//! material …) and child entries — the text-format complement to the binary
//! `bincaf` container. The shape is stored topologically with the same scheme
//! `bincaf` uses, but as XML elements: type, then per-type geometry snapshots
//! (vertex points, edge curves, face surfaces with parameter ranges, and the
//! child sub-shape lists). Curves and surfaces cannot be downcast from
//! `Arc<dyn Curve>` / `Arc<dyn Surface>`, so they are classified by sampling
//! invariants (zero second derivative ⇒ line, periodic ⇒ circle, planar ⇒
//! plane, equidistant samples ⇒ sphere), exactly like `bincaf` / STEP / IGES.
//!
//! Document layout:
//! ```text
//! <?xml version="1.0" encoding="UTF-8"?>
//! <xcaf version="1.0">
//!   <entry name="...">
//!     <attribute kind="name" value="..."/>...
//!     <shape type="Solid">
//!       <shape type="Shell">
//!         <shape type="Face"><surface .../><shape type="Wire">...</shape></shape>
//!       </shape>
//!     </shape>
//!     <entry name="...">...</entry>
//!   </entry>
//! </xcaf>
//! ```
//! Attribute values are XML-escaped (`&` `<` `>` `"` `'`); the hand-rolled
//! parser unescapes them on read and rejects malformed input with `Err`.

use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::gp::{GpAx1, GpAx2, GpAx3, GpCirc, GpDir, GpPln, GpPnt, GpSphere, GpVec};
use occt_geom::{Curve, GeomCircle, GeomLine, GeomPlane, GeomSphere, Surface};

use crate::abs::ShapeType;
use crate::brep_surface::{classify_surface, face_plane, sphere_center, SurfaceKind};
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, Shell, TopoShape, Wire};
use crate::tgeometry::GeometryRegistry;

/// A named attribute attached to an assembly entry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct XmlAttribute {
    /// e.g. `"name"`, `"color"`, `"layer"`, `"material"`.
    pub kind: String,
    /// Attribute payload (free-form string, e.g. `"1.0,0.0,0.0"`).
    pub value: String,
}

/// One node of the assembly tree: a name, an optional shape plus
/// attributes/children.
#[derive(Debug, Clone, Default)]
pub struct XmlEntry {
    pub name: String,
    pub shape: Option<TopoShape>,
    pub attributes: Vec<XmlAttribute>,
    pub children: Vec<XmlEntry>,
}

/// XML XCAF document — a root entry wrapping the whole assembly tree.
#[derive(Debug, Clone)]
pub struct XmlXcafDoc {
    pub root: XmlEntry,
    pub version: String,
}

impl Default for XmlXcafDoc {
    fn default() -> Self {
        Self {
            root: XmlEntry::default(),
            version: "1.0".to_string(),
        }
    }
}

/// Serialize a `XmlXcafDoc` document to well-formed XML 1.0 text.
pub fn to_xml(doc: &XmlXcafDoc) -> Result<String, String> {
    let mut out = String::new();
    out.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    out.push_str(&format!("<xcaf version=\"{}\">\n", escape_xml(&doc.version)));
    write_entry(&mut out, &doc.root, 1);
    out.push_str("</xcaf>\n");
    Ok(out)
}

/// Parse an XML XCAF document back into a [`XmlXcafDoc`].
///
/// Accepts an optional `<?xml …?>` declaration, validates the root element
/// and tag balance, and rejects unknown elements or unterminated input with
/// `Err`.
pub fn from_xml(xml: &str) -> Result<XmlXcafDoc, String> {
    let root = parse_xml_doc(xml)?;
    if root.name != "xcaf" {
        return Err(format!("xmlcaf: expected root <xcaf>, found <{}>", root.name));
    }
    let version = root.attr("version").unwrap_or("1.0").to_string();
    let entries: Vec<&XmlEl> = root.children.iter().filter(|c| c.name == "entry").collect();
    if entries.len() != 1 {
        return Err("xmlcaf: <xcaf> must contain exactly one root <entry>".into());
    }
    let root_entry = entry_from_element(entries[0])?;
    Ok(XmlXcafDoc {
        root: root_entry,
        version,
    })
}

/// Write a `XmlXcafDoc` document to `path`.
pub fn write_xml_file(doc: &XmlXcafDoc, path: &str) -> Result<(), String> {
    let xml = to_xml(doc)?;
    std::fs::write(path, &xml).map_err(|e| format!("xmlcaf: write {path}: {e}"))
}

/// Read a `XmlXcafDoc` document from `path`.
pub fn read_xml_file(path: &str) -> Result<XmlXcafDoc, String> {
    let xml = std::fs::read_to_string(path).map_err(|e| format!("xmlcaf: read {path}: {e}"))?;
    from_xml(&xml)
}

/// XML-escape a string for use in element text or quoted attribute values.
/// Escapes `&`, `<`, `>`, `"` and `'`.
pub fn escape_xml(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for ch in s.chars() {
        match ch {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            _ => out.push(ch),
        }
    }
    out
}

/// Invert [`escape_xml`], decoding the five XML entity references. An unknown
/// entity or a stray `&` without a terminating `;` yields `Err`.
pub fn unescape_xml(s: &str) -> Result<String, String> {
    let b = s.as_bytes();
    let mut out = String::with_capacity(s.len());
    let mut i = 0;
    while i < b.len() {
        if b[i] == b'&' {
            let semi = b[i..]
                .iter()
                .position(|&c| c == b';')
                .ok_or("xmlcaf: '&' without terminating ';'")?;
            let entity = std::str::from_utf8(&b[i..i + semi + 1])
                .map_err(|e| format!("xmlcaf: bad entity bytes: {e}"))?;
            match entity {
                "&amp;" => out.push('&'),
                "&lt;" => out.push('<'),
                "&gt;" => out.push('>'),
                "&quot;" => out.push('"'),
                "&apos;" => out.push('\''),
                other => return Err(format!("xmlcaf: unknown entity '{other}'")),
            }
            i += semi + 1;
        } else {
            let n = utf8_len(b[i]);
            let chunk = std::str::from_utf8(&b[i..i + n])
                .map_err(|e| format!("xmlcaf: invalid utf-8: {e}"))?;
            out.push_str(chunk);
            i += n;
        }
    }
    Ok(out)
}

/// Serialize a `TopoShape` to a (possibly multi-line) `<shape …>` XML element.
///
/// Mirrors `bincaf`'s topological shape scheme: vertex points become
/// attributes of a self-closing `<shape type="Vertex" …/>`, edges carry a
/// `<curve>` plus vertex children, faces carry a `<surface>` plus wire
/// children, and shells/solids nest their sub-shapes.
pub fn shape_to_xml_element(shape: &TopoShape) -> String {
    let mut out = String::new();
    write_shape_xml(&mut out, shape, 0);
    out
}

/// Rebuild a `TopoShape` from a `<shape …>` opening tag and its inner
/// content. `tag` is the opening tag (`<shape type="Solid">`), `inner` is the
/// text between it and the matching `</shape>`.
pub fn shape_from_xml_element(tag: &str, inner: &str) -> Result<TopoShape, String> {
    let (name, attrs) = parse_attrs(tag)?;
    if name != "shape" {
        return Err(format!("xmlcaf: expected <shape>, found <{name}>"));
    }
    let children = parse_element_list(inner)?;
    let el = XmlEl {
        name,
        attrs,
        children,
    };
    shape_from_element(&el)
}

// ---------------------------------------------------------------------------
// Entry writer
// ---------------------------------------------------------------------------

fn write_entry(out: &mut String, e: &XmlEntry, depth: usize) {
    indent(out, depth);
    out.push_str(&format!("<entry name=\"{}\">\n", escape_xml(&e.name)));
    for a in &e.attributes {
        indent(out, depth + 1);
        out.push_str(&format!(
            "<attribute kind=\"{}\" value=\"{}\"/>\n",
            escape_xml(&a.kind),
            escape_xml(&a.value)
        ));
    }
    if let Some(shape) = &e.shape {
        out.push_str(&reindent(&shape_to_xml_element(shape), depth + 1));
    }
    for c in &e.children {
        write_entry(out, c, depth + 1);
    }
    indent(out, depth);
    out.push_str("</entry>\n");
}

/// Shift a multi-line block right by `depth` levels of two-space indent,
/// preserving its internal relative indentation.
fn reindent(block: &str, depth: usize) -> String {
    let pad = "  ".repeat(depth);
    block
        .lines()
        .map(|l| if l.trim().is_empty() { String::new() } else { format!("{pad}{l}") })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n"
}

fn indent(out: &mut String, depth: usize) {
    for _ in 0..depth {
        out.push_str("  ");
    }
}

// ---------------------------------------------------------------------------
// Shape writer
// ---------------------------------------------------------------------------

fn write_shape_xml(out: &mut String, s: &TopoShape, depth: usize) {
    let reg = GeometryRegistry::global();
    match s.shape_type() {
        ShapeType::Vertex => {
            let p = reg.vertex_point(s);
            let tol = reg.vertex_tolerance(s);
            indent(out, depth);
            out.push_str(&format!(
                "<shape type=\"Vertex\" x=\"{}\" y=\"{}\" z=\"{}\" tol=\"{}\"/>\n",
                fmt(p.x()),
                fmt(p.y()),
                fmt(p.z()),
                fmt(tol)
            ));
        }
        ShapeType::Edge => {
            indent(out, depth);
            out.push_str("<shape type=\"Edge\">\n");
            write_curve_xml(out, s, depth + 1);
            for v in children_of_type(s, ShapeType::Vertex) {
                write_shape_xml(out, &v, depth + 1);
            }
            indent(out, depth);
            out.push_str("</shape>\n");
        }
        ShapeType::Wire => {
            indent(out, depth);
            out.push_str("<shape type=\"Wire\">\n");
            for e in children_of_type(s, ShapeType::Edge) {
                write_shape_xml(out, &e, depth + 1);
            }
            indent(out, depth);
            out.push_str("</shape>\n");
        }
        ShapeType::Face => {
            indent(out, depth);
            out.push_str("<shape type=\"Face\">\n");
            write_surface_xml(out, s, depth + 1);
            for w in children_of_type(s, ShapeType::Wire) {
                write_shape_xml(out, &w, depth + 1);
            }
            indent(out, depth);
            out.push_str("</shape>\n");
        }
        ShapeType::Shell => {
            indent(out, depth);
            out.push_str("<shape type=\"Shell\">\n");
            for f in children_of_type(s, ShapeType::Face) {
                write_shape_xml(out, &f, depth + 1);
            }
            indent(out, depth);
            out.push_str("</shape>\n");
        }
        ShapeType::Solid => {
            indent(out, depth);
            out.push_str("<shape type=\"Solid\">\n");
            for sh in children_of_type(s, ShapeType::Shell) {
                write_shape_xml(out, &sh, depth + 1);
            }
            indent(out, depth);
            out.push_str("</shape>\n");
        }
        // Compound, CompSolid, Shape: arbitrary child list.
        _ => {
            indent(out, depth);
            out.push_str(&format!("<shape type=\"{}\">\n", s.shape_type().to_str()));
            let kids: Vec<TopoShape> = s
                .tshape
                .read()
                .unwrap()
                .children
                .clone();
            for k in &kids {
                write_shape_xml(out, k, depth + 1);
            }
            indent(out, depth);
            out.push_str("</shape>\n");
        }
    }
}

fn write_curve_xml(out: &mut String, s: &TopoShape, depth: usize) {
    let reg = GeometryRegistry::global();
    let (lo, hi) = finite_range(reg.edge_parameters(s));
    indent(out, depth);
    let Some(curve) = reg.edge_curve(s) else {
        out.push_str(&format!(
            "<curve kind=\"line\" x=\"0\" y=\"0\" z=\"0\" dx=\"1\" dy=\"0\" dz=\"0\" lo=\"{}\" hi=\"{}\"/>\n",
            fmt(lo),
            fmt(hi)
        ));
        return;
    };
    match classify_curve(curve.as_ref(), lo, hi) {
        CurveKind::Line => {
            let origin = curve.d0(0.0);
            let d = dir_of(&curve.d1(0.0).1);
            out.push_str(&format!(
                "<curve kind=\"line\" x=\"{}\" y=\"{}\" z=\"{}\" dx=\"{}\" dy=\"{}\" dz=\"{}\" lo=\"{}\" hi=\"{}\"/>\n",
                fmt(origin.x()),
                fmt(origin.y()),
                fmt(origin.z()),
                fmt(d.x()),
                fmt(d.y()),
                fmt(d.z()),
                fmt(lo),
                fmt(hi)
            ));
        }
        CurveKind::Circle => {
            let p0 = curve.d0(lo);
            let p1 = curve.d0(lo + PI / 2.0);
            let p2 = curve.d0(lo + PI);
            let center = circle_center3(&p0, &p1, &p2).unwrap_or_else(GpPnt::zero);
            let r = center.distance(&p0);
            let n = GpVec::from_pnts(&p0, &p1).xyz().crossed(GpVec::from_pnts(&p0, &p2).xyz());
            let normal = GpDir::from_xyz(&n).unwrap_or(dir_z());
            let xdir = GpDir::from_vec(&GpVec::from_pnts(&center, &p0)).unwrap_or(dir_x());
            out.push_str(&format!(
                "<curve kind=\"circle\" cx=\"{}\" cy=\"{}\" cz=\"{}\" nx=\"{}\" ny=\"{}\" nz=\"{}\" xdx=\"{}\" xdy=\"{}\" xdz=\"{}\" r=\"{}\" lo=\"{}\" hi=\"{}\"/>\n",
                fmt(center.x()),
                fmt(center.y()),
                fmt(center.z()),
                fmt(normal.x()),
                fmt(normal.y()),
                fmt(normal.z()),
                fmt(xdir.x()),
                fmt(xdir.y()),
                fmt(xdir.z()),
                fmt(r),
                fmt(lo),
                fmt(hi)
            ));
        }
        CurveKind::Other => {
            // ponytail: B-spline curves serialize as their tangent line.
            let origin = curve.d0(lo);
            let d = dir_of(&curve.d1(lo).1);
            out.push_str(&format!(
                "<curve kind=\"line\" x=\"{}\" y=\"{}\" z=\"{}\" dx=\"{}\" dy=\"{}\" dz=\"{}\" lo=\"{}\" hi=\"{}\"/>\n",
                fmt(origin.x()),
                fmt(origin.y()),
                fmt(origin.z()),
                fmt(d.x()),
                fmt(d.y()),
                fmt(d.z()),
                fmt(lo),
                fmt(hi)
            ));
        }
    }
}

fn write_surface_xml(out: &mut String, s: &TopoShape, depth: usize) {
    let reg = GeometryRegistry::global();
    indent(out, depth);
    let Some(surf) = reg.face_surface(s) else {
        out.push_str(
            "<surface kind=\"plane\" x=\"0\" y=\"0\" z=\"0\" nx=\"0\" ny=\"0\" nz=\"1\" u0=\"0\" u1=\"1\" v0=\"0\" v1=\"1\"/>\n",
        );
        return;
    };
    let (u0, u1, v0, v1) = clamp_uv(surf.u_range().0, surf.u_range().1, surf.v_range().0, surf.v_range().1);
    match classify_surface(surf.as_ref()) {
        SurfaceKind::Plane => {
            let pln = face_plane(&Face(s.clone())).unwrap_or_default();
            out.push_str(&format!(
                "<surface kind=\"plane\" x=\"{}\" y=\"{}\" z=\"{}\" nx=\"{}\" ny=\"{}\" nz=\"{}\" u0=\"{}\" u1=\"{}\" v0=\"{}\" v1=\"{}\"/>\n",
                fmt(pln.location().x()),
                fmt(pln.location().y()),
                fmt(pln.location().z()),
                fmt(pln.axis().direction().x()),
                fmt(pln.axis().direction().y()),
                fmt(pln.axis().direction().z()),
                fmt(u0),
                fmt(u1),
                fmt(v0),
                fmt(v1)
            ));
        }
        SurfaceKind::Sphere => {
            let center = sphere_center(surf.as_ref()).unwrap_or_else(GpPnt::zero);
            let r = surf.d0(u0, 0.5 * (v0 + v1)).distance(&center);
            out.push_str(&format!(
                "<surface kind=\"sphere\" cx=\"{}\" cy=\"{}\" cz=\"{}\" r=\"{}\" u0=\"{}\" u1=\"{}\" v0=\"{}\" v1=\"{}\"/>\n",
                fmt(center.x()),
                fmt(center.y()),
                fmt(center.z()),
                fmt(r),
                fmt(u0),
                fmt(u1),
                fmt(v0),
                fmt(v1)
            ));
        }
        _ => {
            // ponytail: cylinder/cone/torus surfaces store a plane fallback;
            // add ring-sampled parameters when a curved-face export needs them.
            let pln = face_plane(&Face(s.clone())).unwrap_or_default();
            out.push_str(&format!(
                "<surface kind=\"plane\" x=\"{}\" y=\"{}\" z=\"{}\" nx=\"{}\" ny=\"{}\" nz=\"{}\" u0=\"{}\" u1=\"{}\" v0=\"{}\" v1=\"{}\"/>\n",
                fmt(pln.location().x()),
                fmt(pln.location().y()),
                fmt(pln.location().z()),
                fmt(pln.axis().direction().x()),
                fmt(pln.axis().direction().y()),
                fmt(pln.axis().direction().z()),
                fmt(u0),
                fmt(u1),
                fmt(v0),
                fmt(v1)
            ));
        }
    }
}

/// Shortest round-trippable decimal form of a float.
fn fmt(v: f64) -> String {
    format!("{}", v)
}

// ---------------------------------------------------------------------------
// Shape reader
// ---------------------------------------------------------------------------

fn shape_from_element(el: &XmlEl) -> Result<TopoShape, String> {
    let ty_str = el.attr("type").ok_or("xmlcaf: <shape> missing 'type' attribute")?;
    let ty = ShapeType::from_str(ty_str).ok_or_else(|| format!("xmlcaf: unknown shape type '{ty_str}'"))?;
    let b = TopoBuilder::new();
    match ty {
        ShapeType::Vertex => {
            let x = num_attr(el, "x")?;
            let y = num_attr(el, "y")?;
            let z = num_attr(el, "z")?;
            let tol = num_attr(el, "tol")?;
            Ok(b.make_vertex(GpPnt::new(x, y, z), tol).0)
        }
        ShapeType::Edge => {
            let curve_el = el.child("curve").ok_or("xmlcaf: <shape type=\"Edge\"> missing <curve>")?;
            let (curve, lo, hi) = curve_from_element(curve_el)?;
            let mut e = b.make_edge(curve, lo, hi);
            for v in el.child_shapes() {
                let v = shape_from_element(v)?;
                b.add(&mut e.0, &v);
            }
            Ok(e.0)
        }
        ShapeType::Wire => {
            let mut edges = Vec::new();
            for s in el.child_shapes() {
                edges.push(Edge(shape_from_element(s)?));
            }
            Ok(b.make_wire(&edges).0)
        }
        ShapeType::Face => {
            let surf_el = el.child("surface").ok_or("xmlcaf: <shape type=\"Face\"> missing <surface>")?;
            let (surface, _uv) = surface_from_element(surf_el)?;
            let mut wires = Vec::new();
            for s in el.child_shapes() {
                wires.push(Wire(shape_from_element(s)?));
            }
            Ok(b.make_face(surface, &wires).0)
        }
        ShapeType::Shell => {
            let mut faces = Vec::new();
            for s in el.child_shapes() {
                faces.push(Face(shape_from_element(s)?));
            }
            Ok(b.make_shell(&faces).0)
        }
        ShapeType::Solid => {
            let mut shells = Vec::new();
            for s in el.child_shapes() {
                shells.push(Shell(shape_from_element(s)?));
            }
            Ok(b.make_solid(&shells).0)
        }
        // Compound, CompSolid, Shape: arbitrary child list.
        _ => {
            let mut comp = TopoShape::new(ty);
            for s in el.child_shapes() {
                let k = shape_from_element(s)?;
                b.add(&mut comp, &k);
            }
            Ok(comp)
        }
    }
}

fn curve_from_element(el: &XmlEl) -> Result<(Arc<dyn Curve>, f64, f64), String> {
    let kind = el.attr("kind").unwrap_or("line").to_string();
    let lo = num_attr(el, "lo")?;
    let hi = num_attr(el, "hi")?;
    let curve: Arc<dyn Curve> = match kind.as_str() {
        "circle" => {
            let cx = num_attr(el, "cx")?;
            let cy = num_attr(el, "cy")?;
            let cz = num_attr(el, "cz")?;
            let nx = num_attr(el, "nx")?;
            let ny = num_attr(el, "ny")?;
            let nz = num_attr(el, "nz")?;
            let xdx = num_attr(el, "xdx")?;
            let xdy = num_attr(el, "xdy")?;
            let xdz = num_attr(el, "xdz")?;
            let r = num_attr(el, "r")?;
            let ax2 = GpAx2::new(
                GpPnt::new(cx, cy, cz),
                dir3(nx, ny, nz)?,
                dir3(xdx, xdy, xdz)?,
            )
            .map_err(|e| format!("xmlcaf: bad circle axis: {e}"))?;
            Arc::new(GeomCircle::new(GpCirc::new(ax2, r)))
        }
        // Any other kind (or a missing one) decodes as a line.
        _ => {
            let x = num_attr(el, "x")?;
            let y = num_attr(el, "y")?;
            let z = num_attr(el, "z")?;
            let dx = num_attr(el, "dx")?;
            let dy = num_attr(el, "dy")?;
            let dz = num_attr(el, "dz")?;
            Arc::new(GeomLine::from_pnt_dir(GpPnt::new(x, y, z), dir3(dx, dy, dz)?))
        }
    };
    Ok((curve, lo, hi))
}

fn surface_from_element(el: &XmlEl) -> Result<(Arc<dyn Surface>, (f64, f64, f64, f64)), String> {
    let kind = el.attr("kind").unwrap_or("plane").to_string();
    let surf: Arc<dyn Surface> = match kind.as_str() {
        "sphere" => {
            let cx = num_attr(el, "cx")?;
            let cy = num_attr(el, "cy")?;
            let cz = num_attr(el, "cz")?;
            let r = num_attr(el, "r")?;
            let ax3 = GpAx3::from_ax1(&GpAx1::new(GpPnt::new(cx, cy, cz), dir_z()));
            let sphere = GpSphere::new(ax3, r).map_err(|e| format!("xmlcaf: {e}"))?;
            Arc::new(GeomSphere::new(sphere))
        }
        // Any other kind (or a missing one) decodes as a plane.
        _ => {
            let x = num_attr(el, "x")?;
            let y = num_attr(el, "y")?;
            let z = num_attr(el, "z")?;
            let nx = num_attr(el, "nx")?;
            let ny = num_attr(el, "ny")?;
            let nz = num_attr(el, "nz")?;
            let origin = GpPnt::new(x, y, z);
            let normal = dir3(nx, ny, nz)?;
            let ax3 = GpAx3::from_ax1(&GpAx1::new(origin, normal));
            Arc::new(GeomPlane::new(GpPln::new(ax3)))
        }
    };
    let u0 = num_attr(el, "u0").unwrap_or(0.0);
    let u1 = num_attr(el, "u1").unwrap_or(1.0);
    let v0 = num_attr(el, "v0").unwrap_or(0.0);
    let v1 = num_attr(el, "v1").unwrap_or(1.0);
    Ok((surf, (u0, u1, v0, v1)))
}

fn num_attr(el: &XmlEl, k: &str) -> Result<f64, String> {
    let v = el.attr(k).ok_or_else(|| format!("xmlcaf: <{}> missing attribute '{k}'", el.name))?;
    v.trim()
        .parse::<f64>()
        .map_err(|e| format!("xmlcaf: bad number '{v}' for '{k}': {e}"))
}

fn dir3(x: f64, y: f64, z: f64) -> Result<GpDir, String> {
    GpDir::new(x, y, z).map_err(|e| format!("xmlcaf: bad direction: {e}"))
}

// ---------------------------------------------------------------------------
// Minimal XML parser
// ---------------------------------------------------------------------------

/// A parsed XML element: name, unescaped attributes, and child elements.
#[derive(Debug, Clone)]
struct XmlEl {
    name: String,
    attrs: Vec<(String, String)>,
    children: Vec<XmlEl>,
}

impl XmlEl {
    fn attr(&self, k: &str) -> Option<&str> {
        self.attrs.iter().find(|(a, _)| a == k).map(|(_, v)| v.as_str())
    }
    fn child(&self, name: &str) -> Option<&XmlEl> {
        self.children.iter().find(|c| c.name == name)
    }
    fn child_shapes(&self) -> Vec<&XmlEl> {
        self.children.iter().filter(|c| c.name == "shape").collect()
    }
}

struct Cursor<'a> {
    s: &'a str,
    pos: usize,
}

impl<'a> Cursor<'a> {
    fn rest(&self) -> &'a str {
        &self.s[self.pos..]
    }
    fn skip_ws(&mut self) {
        while let Some(b) = self.s.as_bytes().get(self.pos) {
            if b.is_ascii_whitespace() {
                self.pos += 1;
            } else {
                break;
            }
        }
    }
    fn parse_name(&mut self) -> Result<String, String> {
        let b = self.s.as_bytes();
        let start = self.pos;
        while let Some(&c) = b.get(self.pos) {
            if c.is_ascii_alphanumeric() || c == b'_' || c == b'-' || c == b'.' {
                self.pos += 1;
            } else {
                break;
            }
        }
        if self.pos == start {
            return Err("xmlcaf: expected a name".into());
        }
        Ok(self.s[start..self.pos].to_string())
    }
    fn parse_element(&mut self) -> Result<XmlEl, String> {
        self.skip_ws();
        if !self.rest().starts_with('<') {
            return Err("xmlcaf: expected '<'".into());
        }
        self.pos += 1;
        let name = self.parse_name()?;
        let mut attrs = Vec::new();
        let mut self_closing = false;
        loop {
            self.skip_ws();
            let rest = self.rest();
            if rest.starts_with("/>") {
                self.pos += 2;
                self_closing = true;
                break;
            }
            if rest.starts_with('>') {
                self.pos += 1;
                break;
            }
            let aname = self.parse_name()?;
            self.skip_ws();
            if !self.rest().starts_with('=') {
                return Err("xmlcaf: expected '=' after attribute name".into());
            }
            self.pos += 1;
            self.skip_ws();
            if !self.rest().starts_with('"') {
                return Err("xmlcaf: expected '\"' before attribute value".into());
            }
            self.pos += 1;
            let end = self.s[self.pos..]
                .find('"')
                .ok_or("xmlcaf: unterminated attribute value")?;
            let raw = &self.s[self.pos..self.pos + end];
            let value = unescape_xml(raw)?;
            self.pos += end + 1;
            attrs.push((aname, value));
        }
        if self_closing {
            return Ok(XmlEl { name, attrs, children: vec![] });
        }
        let mut children = Vec::new();
        loop {
            self.skip_ws();
            let rest = self.rest();
            if rest.starts_with("</") {
                self.pos += 2;
                let cname = self.parse_name()?;
                self.skip_ws();
                if !self.rest().starts_with('>') {
                    return Err("xmlcaf: expected '>' in closing tag".into());
                }
                self.pos += 1;
                if cname != name {
                    return Err(format!("xmlcaf: mismatched closing tag </{cname}> for <{name}>"));
                }
                return Ok(XmlEl { name, attrs, children });
            }
            if rest.is_empty() {
                return Err(format!("xmlcaf: unexpected end of input inside <{name}>"));
            }
            children.push(self.parse_element()?);
        }
    }
}

/// Parse an optional `<?xml …?>` declaration followed by the single root
/// element. Trailing content after the root is rejected.
fn parse_xml_doc(s: &str) -> Result<XmlEl, String> {
    let mut c = Cursor { s, pos: 0 };
    c.skip_ws();
    if c.rest().starts_with("<?xml") {
        let end = c.rest().find("?>").ok_or("xmlcaf: unterminated XML declaration")?;
        c.pos += end + 2;
    }
    c.skip_ws();
    let root = c.parse_element()?;
    c.skip_ws();
    if !c.rest().is_empty() {
        return Err("xmlcaf: trailing content after root element".into());
    }
    Ok(root)
}

/// Parse a sequence of sibling elements (used for the inner content of a
/// `<shape>` element).
fn parse_element_list(s: &str) -> Result<Vec<XmlEl>, String> {
    let mut c = Cursor { s, pos: 0 };
    let mut els = Vec::new();
    loop {
        c.skip_ws();
        if c.rest().is_empty() {
            return Ok(els);
        }
        els.push(c.parse_element()?);
    }
}

/// Parse the name and attributes out of a single opening tag (which may end
/// in `>` or `/>`). Child content is not consumed.
fn parse_attrs(s: &str) -> Result<(String, Vec<(String, String)>), String> {
    let mut c = Cursor { s, pos: 0 };
    c.skip_ws();
    if !c.rest().starts_with('<') {
        return Err("xmlcaf: expected '<'".into());
    }
    c.pos += 1;
    let name = c.parse_name()?;
    let mut attrs = Vec::new();
    loop {
        c.skip_ws();
        let rest = c.rest();
        if rest.starts_with("/>") || rest.starts_with('>') {
            break;
        }
        let aname = c.parse_name()?;
        c.skip_ws();
        if !c.rest().starts_with('=') {
            return Err("xmlcaf: expected '=' after attribute name".into());
        }
        c.pos += 1;
        c.skip_ws();
        if !c.rest().starts_with('"') {
            return Err("xmlcaf: expected '\"' before attribute value".into());
        }
        c.pos += 1;
        let end = c.s[c.pos..]
            .find('"')
            .ok_or("xmlcaf: unterminated attribute value")?;
        let value = unescape_xml(&c.s[c.pos..c.pos + end])?;
        c.pos += end + 1;
        attrs.push((aname, value));
    }
    Ok((name, attrs))
}

/// Length in bytes of the UTF-8 character whose leading byte is `lead`.
fn utf8_len(lead: u8) -> usize {
    if lead < 0x80 {
        1
    } else if lead < 0xE0 {
        2
    } else if lead < 0xF0 {
        3
    } else {
        4
    }
}

fn entry_from_element(el: &XmlEl) -> Result<XmlEntry, String> {
    if el.name != "entry" {
        return Err(format!("xmlcaf: expected <entry>, found <{}>", el.name));
    }
    let name = el.attr("name").unwrap_or("").to_string();
    let mut attributes = Vec::new();
    let mut shape = None;
    let mut children = Vec::new();
    for c in &el.children {
        match c.name.as_str() {
            "attribute" => {
                let kind = c.attr("kind").unwrap_or("").to_string();
                let value = c.attr("value").unwrap_or("").to_string();
                attributes.push(XmlAttribute { kind, value });
            }
            "shape" => {
                shape = Some(shape_from_element(c)?);
            }
            "entry" => {
                children.push(entry_from_element(c)?);
            }
            other => return Err(format!("xmlcaf: unexpected <{other}> inside <entry>")),
        }
    }
    Ok(XmlEntry {
        name,
        shape,
        attributes,
        children,
    })
}

// ---------------------------------------------------------------------------
// Geometry classification / helpers (mirrors bincaf.rs)
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq, Eq)]
enum CurveKind {
    Line,
    Circle,
    Other,
}

fn classify_curve(c: &dyn Curve, a: f64, b: f64) -> CurveKind {
    let (lo, hi) = if a.is_finite() && b.is_finite() && b > a {
        (a, b)
    } else {
        (0.0, 1.0)
    };
    let span = hi - lo;
    if span < 1e-12 {
        return CurveKind::Other;
    }
    let mut d2s = Vec::with_capacity(6);
    for i in 0..6 {
        let u = lo + span * i as f64 / 5.0;
        d2s.push(c.d2(u).2.magnitude());
    }
    let max_d2 = d2s.iter().cloned().fold(0.0_f64, f64::max);
    if max_d2 < 1e-9 {
        return CurveKind::Line;
    }
    let min_d2 = d2s.iter().cloned().fold(f64::INFINITY, f64::min);
    if c.is_periodic() || (max_d2 - min_d2) / max_d2 < 0.02 {
        CurveKind::Circle
    } else {
        CurveKind::Other
    }
}

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

fn finite_range((a, b): (f64, f64)) -> (f64, f64) {
    if a.is_finite() && b.is_finite() && b > a {
        (a, b)
    } else {
        (0.0, 1.0)
    }
}

fn clamp_uv(u0: f64, u1: f64, v0: f64, v1: f64) -> (f64, f64, f64, f64) {
    let c = |a: f64, b: f64| if a.is_finite() && b.is_finite() && b > a { (a, b) } else { (-1.0, 1.0) };
    let (u0, u1) = c(u0, u1);
    let (v0, v1) = c(v0, v1);
    (u0, u1, v0, v1)
}

fn dir_x() -> GpDir {
    GpDir::new(1.0, 0.0, 0.0).expect("unit x")
}
fn dir_z() -> GpDir {
    GpDir::new(0.0, 0.0, 1.0).expect("unit z")
}

fn dir_of(v: &GpVec) -> GpDir {
    GpDir::from_vec(v).unwrap_or(dir_x())
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::{BRepPrimBox, BRepPrimSphere};
    use crate::topo_tools_full::faces_of;

    fn attr(kind: &str, value: &str) -> XmlAttribute {
        XmlAttribute {
            kind: kind.into(),
            value: value.into(),
        }
    }

    fn doc_with_attrs() -> XmlXcafDoc {
        XmlXcafDoc {
            version: "1.0".into(),
            root: XmlEntry {
                name: "Root".into(),
                shape: None,
                attributes: vec![
                    attr("name", "Root"),
                    attr("color", "1.0,0.0,0.0"),
                    attr("layer", "L1"),
                ],
                children: vec![XmlEntry {
                    name: "Child".into(),
                    shape: None,
                    attributes: vec![attr("name", "Child")],
                    children: vec![],
                }],
            },
        }
    }

    #[test]
    fn xml_header_and_escape() {
        let s = "a&b<c>d\"e'f";
        let esc = escape_xml(s);
        assert_eq!(esc, "a&amp;b&lt;c&gt;d&quot;e&apos;f");
        assert_eq!(unescape_xml(&esc).unwrap(), s);
        assert_eq!(escape_xml("1 < 2 && 3 > 0"), "1 &lt; 2 &amp;&amp; 3 &gt; 0");
    }

    #[test]
    fn to_xml_wellformed() {
        let doc = doc_with_attrs();
        let xml = to_xml(&doc).unwrap();
        assert!(xml.starts_with("<?xml version=\"1.0\" encoding=\"UTF-8\"?>"), "header:\n{xml}");
        assert!(xml.contains("<xcaf version=\"1.0\">"));
        // <entry> and <shape> are written with matching open/close tags.
        for tag in ["entry", "shape"] {
            let open = xml.matches(&format!("<{tag}")).count();
            let close = xml.matches(&format!("</{tag}>")).count();
            assert_eq!(open, close, "unbalanced <{tag}>");
        }
        // <attribute> is self-closing: it must never appear with a close tag.
        assert!(xml.contains("<attribute"));
        assert_eq!(xml.matches("</attribute>").count(), 0);
    }

    #[test]
    fn roundtrip_attributes() {
        let doc = doc_with_attrs();
        let got = from_xml(&to_xml(&doc).unwrap()).unwrap();
        assert_eq!(got.version, "1.0");
        assert_eq!(got.root.name, "Root");
        assert_eq!(got.root.attributes, doc.root.attributes);
        assert_eq!(got.root.children.len(), 1);
        assert_eq!(got.root.children[0].attributes, vec![attr("name", "Child")]);
        assert!(got.root.shape.is_none());
    }

    #[test]
    fn roundtrip_nested_tree() {
        let doc = XmlXcafDoc {
            version: "1.0".into(),
            root: XmlEntry {
                name: "R".into(),
                shape: None,
                attributes: vec![],
                children: vec![
                    XmlEntry {
                        name: "C1".into(),
                        shape: None,
                        attributes: vec![attr("name", "c1")],
                        children: vec![XmlEntry {
                            name: "GC".into(),
                            shape: None,
                            attributes: vec![attr("name", "gc")],
                            children: vec![],
                        }],
                    },
                    XmlEntry {
                        name: "C2".into(),
                        shape: None,
                        attributes: vec![],
                        children: vec![],
                    },
                ],
            },
        };
        let got = from_xml(&to_xml(&doc).unwrap()).unwrap();
        assert_eq!(got.root.children.len(), 2);
        assert_eq!(got.root.children[0].name, "C1");
        assert_eq!(got.root.children[0].children.len(), 1);
        assert_eq!(got.root.children[0].children[0].name, "GC");
        assert_eq!(got.root.children[0].children[0].attributes[0].value, "gc");
        assert_eq!(got.root.children[1].name, "C2");
    }

    #[test]
    fn roundtrip_box_shape() {
        let b = BRepPrimBox::make_box(2.0, 3.0, 4.0);
        let doc = XmlXcafDoc {
            version: "1.0".into(),
            root: XmlEntry {
                name: "Box".into(),
                shape: Some(b.solid.0),
                attributes: vec![],
                children: vec![],
            },
        };
        let xml = to_xml(&doc).unwrap();
        assert!(xml.contains("<shape type=\"Solid\">"), "solid in:\n{xml}");
        let got = from_xml(&xml).unwrap();
        let shape = got.root.shape.expect("shape");
        assert_eq!(shape.shape_type(), ShapeType::Solid);
        assert_eq!(faces_of(&shape).len(), 6, "box face count");
    }

    #[test]
    fn roundtrip_sphere_shape() {
        let s = BRepPrimSphere::make_sphere(2.5);
        let doc = XmlXcafDoc {
            version: "1.0".into(),
            root: XmlEntry {
                name: "Sphere".into(),
                shape: Some(s.solid.0),
                attributes: vec![],
                children: vec![],
            },
        };
        let xml = to_xml(&doc).unwrap();
        assert!(xml.contains("kind=\"sphere\""), "sphere surface in:\n{xml}");
        let got = from_xml(&xml).unwrap();
        let shape = got.root.shape.expect("shape");
        assert_eq!(shape.shape_type(), ShapeType::Solid);
        assert_eq!(faces_of(&shape).len(), 1, "sphere face count");
    }

    #[test]
    fn malformed_xml_rejected() {
        assert!(from_xml("not xml at all").is_err());
        assert!(from_xml("<?xml version=\"1.0\"?><xcaf>").is_err());
        assert!(from_xml("<xcaf><entry name=\"a\"></entry>").is_err());
        assert!(from_xml("<xcaf><entry name=\"a\"><foo/></entry></xcaf>").is_err());
        assert!(from_xml("<xcaf><entry name=\"a\"/><entry name=\"b\"/></xcaf>").is_err());
        assert!(from_xml("<xcaf><entry name=\"a\"></entry></xcaf> junk").is_err());
        assert!(from_xml("<xcaf><entry name=\"a\"><entry></entry></xcf></xcaf>").is_err());
        assert!(unescape_xml("&bogus;").is_err());
        assert!(unescape_xml("a & b").is_err());
    }

    #[test]
    fn xml_escape_roundtrip() {
        let tricky = "Widget <A> & \"B\" 'C' > 5";
        let doc = XmlXcafDoc {
            version: "1.0".into(),
            root: XmlEntry {
                name: "Root".into(),
                shape: None,
                attributes: vec![attr("name", tricky)],
                children: vec![],
            },
        };
        let xml = to_xml(&doc).unwrap();
        assert!(xml.contains("&lt;A&gt; &amp; &quot;B&quot; &apos;C&apos; &gt; 5"), "escaped:\n{xml}");
        let got = from_xml(&xml).unwrap();
        assert_eq!(got.root.attributes[0].value, tricky);
    }

    #[test]
    fn file_roundtrip() {
        let b = BRepPrimBox::make_box(1.0, 2.0, 3.0);
        let doc = XmlXcafDoc {
            version: "1.0".into(),
            root: XmlEntry {
                name: "Box".into(),
                shape: Some(b.solid.0),
                attributes: vec![attr("name", "Box")],
                children: vec![],
            },
        };
        let path = std::env::temp_dir().join("occt_xmlcaf_test.xml");
        let p = path.to_str().unwrap();
        write_xml_file(&doc, p).unwrap();
        let got = read_xml_file(p).unwrap();
        assert_eq!(got.root.attributes[0].value, "Box");
        let shape = got.root.shape.expect("shape");
        assert_eq!(shape.shape_type(), ShapeType::Solid);
        assert_eq!(faces_of(&shape).len(), 6);
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn empty_doc_ok() {
        let doc = XmlXcafDoc::default();
        let got = from_xml(&to_xml(&doc).unwrap()).unwrap();
        assert_eq!(got.version, "1.0");
        assert_eq!(got.root.name, "");
        assert!(got.root.shape.is_none());
        assert!(got.root.attributes.is_empty());
        assert!(got.root.children.is_empty());
    }

    #[test]
    fn shape_element_helpers_roundtrip() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let shape = b.solid.0.clone();
        let xml = shape_to_xml_element(&shape);
        // `tag` is the opening tag; `inner` is everything up to the final
        // closing `</shape>` (the root shape's own close tag).
        let tag = xml.lines().next().unwrap_or("").to_string();
        let open_end = xml.find('>').unwrap();
        let close_start = xml.rfind("</shape>").unwrap();
        let inner = &xml[open_end + 1..close_start];
        let got = shape_from_xml_element(&tag, inner).unwrap();
        assert_eq!(got.shape_type(), ShapeType::Solid);
        assert_eq!(faces_of(&got).len(), 6);
    }
}
