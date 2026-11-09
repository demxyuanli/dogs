use super::prelude::*;
use super::*;

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

pub(super) fn write_entry(out: &mut String, e: &XmlEntry, depth: usize) {
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
pub(super) fn reindent(block: &str, depth: usize) -> String {
    let pad = "  ".repeat(depth);
    block
        .lines()
        .map(|l| if l.trim().is_empty() { String::new() } else { format!("{pad}{l}") })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n"
}

pub(super) fn indent(out: &mut String, depth: usize) {
    for _ in 0..depth {
        out.push_str("  ");
    }
}

// ---------------------------------------------------------------------------
// Shape writer
// ---------------------------------------------------------------------------

pub(super) fn write_shape_xml(out: &mut String, s: &TopoShape, depth: usize) {
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

pub(super) fn write_curve_xml(out: &mut String, s: &TopoShape, depth: usize) {
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

pub(super) fn write_surface_xml(out: &mut String, s: &TopoShape, depth: usize) {
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
pub(super) fn fmt(v: f64) -> String {
    format!("{}", v)
}

// ---------------------------------------------------------------------------
// Shape reader
// ---------------------------------------------------------------------------

pub(super) fn shape_from_element(el: &XmlEl) -> Result<TopoShape, String> {
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

pub(super) fn curve_from_element(el: &XmlEl) -> Result<(Arc<dyn Curve>, f64, f64), String> {
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

pub(super) fn surface_from_element(el: &XmlEl) -> Result<(Arc<dyn Surface>, (f64, f64, f64, f64)), String> {
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

pub(super) fn num_attr(el: &XmlEl, k: &str) -> Result<f64, String> {
    let v = el.attr(k).ok_or_else(|| format!("xmlcaf: <{}> missing attribute '{k}'", el.name))?;
    v.trim()
        .parse::<f64>()
        .map_err(|e| format!("xmlcaf: bad number '{v}' for '{k}': {e}"))
}

pub(super) fn dir3(x: f64, y: f64, z: f64) -> Result<GpDir, String> {
    GpDir::new(x, y, z).map_err(|e| format!("xmlcaf: bad direction: {e}"))
}

// ---------------------------------------------------------------------------
// Minimal XML parser
// ---------------------------------------------------------------------------

/// A parsed XML element: name, unescaped attributes, and child elements.
#[derive(Debug, Clone)]
pub(super) struct XmlEl {
    pub(super) name: String,
    pub(super) attrs: Vec<(String, String)>,
    pub(super) children: Vec<XmlEl>,
}

impl XmlEl {
    pub(super) fn attr(&self, k: &str) -> Option<&str> {
        self.attrs.iter().find(|(a, _)| a == k).map(|(_, v)| v.as_str())
    }
    pub(super) fn child(&self, name: &str) -> Option<&XmlEl> {
        self.children.iter().find(|c| c.name == name)
    }
    pub(super) fn child_shapes(&self) -> Vec<&XmlEl> {
        self.children.iter().filter(|c| c.name == "shape").collect()
    }
}

pub(super) struct Cursor<'a> {
    pub(super) s: &'a str,
    pub(super) pos: usize,
}
