//! Phase 5 module: xcaf — STEP assembly metadata (`XCAFDoc_ShapeTool`-lite).
//!
//! A minimal stand-in for OCCT's XCAF document: tracks the product list of an
//! assembly and per-product name / color / layer attributes, serializes them
//! alongside a STEP physical file, and reads them back. Because a fully valid
//! XCAF (STEP 214 `APPLICATION_PROTOCOL`) layer is out of scope, the metadata
//! is embedded as a comment block before the `DATA` section using a clear
//! marker, which keeps the STEP file valid and round-trippable.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::f64::consts::PI;

use crate::bincaf::{BinXcaf, BinXcafEntry, XcafAttribute};
use crate::model::BRepModel;
use crate::shape::TopoShape;
use crate::xmlcaf::{XmlEntry, XmlXcafDoc};

/// Assembly metadata: product list plus per-product name/color/layer.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct StepAssembly {
    pub products: Vec<String>,
    pub names: HashMap<String, String>,
    pub colors: HashMap<String, (f64, f64, f64)>,
    pub layers: HashMap<String, String>,
}

impl StepAssembly {
    /// Register a product (idempotent).
    pub fn add_product(&mut self, name: &str) {
        if !self.products.iter().any(|p| p == name) {
            self.products.push(name.to_string());
        }
        self.names.entry(name.to_string()).or_insert_with(|| name.to_string());
    }

    /// Set the RGB color of a product.
    pub fn set_color(&mut self, name: &str, rgb: (f64, f64, f64)) {
        self.add_product(name);
        self.colors.insert(name.to_string(), rgb);
    }

    /// Set the layer of a product.
    pub fn set_layer(&mut self, name: &str, layer: &str) {
        self.add_product(name);
        self.layers.insert(name.to_string(), layer.to_string());
    }

    /// Look up a product by name (returns its display name).
    pub fn find(&self, name: &str) -> Option<&str> {
        self.names.get(name).map(|s| s.as_str())
    }
}

/// Extract the assembly metadata (name/color/layer) from each model shape.
pub fn model_to_step_assembly(model: &BRepModel) -> StepAssembly {
    let mut asm = StepAssembly::default();
    for ms in &model.shapes {
        asm.add_product(&ms.name);
        if let Some(c) = ms.color {
            asm.set_color(&ms.name, (c.r as f64, c.g as f64, c.b as f64));
        }
        if !ms.layer.is_empty() {
            asm.set_layer(&ms.name, &ms.layer);
        }
    }
    asm
}

/// Marker lines delimiting the embedded assembly block.
const BLOCK_BEGIN: &str = "/* ASSEMBLY-BEGIN";
const BLOCK_END: &str = "ASSEMBLY-END";

/// The STEP text for `model` with the assembly metadata embedded as a comment
/// block before the `DATA` section.
///
/// The STEP body itself is produced by `crate::step::write_step`, which
/// already writes the product names into `ADVANCED_BREP_SHAPE_REPRESENTATION`
/// records. The comment block adds the color/layer attributes that plain STEP
/// does not carry. Lines are formatted `name|r,g,b|layer`; names or layers
/// containing a `|` are sanitized (replaced with `/`).
pub fn write_step_with_metadata(model: &BRepModel) -> String {
    let step = crate::step::write_step(model);
    let mut block = String::new();
    block.push_str(BLOCK_BEGIN);
    block.push('\n');
    for ms in &model.shapes {
        let name = sanitize(&ms.name);
        let color = match ms.color {
            Some(c) => format2(c),
            None => String::new(),
        };
        let layer = sanitize(&ms.layer);
        block.push_str(&format!("{name}|{color}|{layer}\n"));
    }
    block.push_str(BLOCK_END);
    block.push_str("\n*/\n");
    step.replacen("DATA;", &format!("{block}DATA;"), 1)
}

/// Format a color as `r,g,b` (rounded to 3 decimals, stable round-trip).
fn format2(c: occt_core::quantity::Color) -> String {
    let r = (c.r * 1000.0).round() / 1000.0;
    let g = (c.g * 1000.0).round() / 1000.0;
    let b = (c.b * 1000.0).round() / 1000.0;
    format!("{r},{g},{b}")
}

/// Replace characters that would break the `|`-delimited comment format.
fn sanitize(s: &str) -> String {
    s.replace(['|', '\n', '\r'], "/")
}

/// Parse the embedded assembly comment block back into a [`StepAssembly`].
///
/// Lines are `name|r,g,b|layer`; empty color/layer fields are skipped.
pub fn extract_metadata_from_step(step_text: &str) -> StepAssembly {
    let mut asm = StepAssembly::default();
    let Some(begin) = step_text.find(BLOCK_BEGIN) else {
        return asm;
    };
    let rest = &step_text[begin + BLOCK_BEGIN.len()..];
    let body = &rest[..rest.find(BLOCK_END).unwrap_or(rest.len())];
    for line in body.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('*') {
            continue;
        }
        let parts: Vec<&str> = line.split('|').collect();
        let name = parts[0].trim();
        if name.is_empty() {
            continue;
        }
        asm.add_product(name);
        if parts.len() >= 2 {
            let rgb: Vec<f64> = parts[1]
                .split(',')
                .filter_map(|s| s.trim().parse().ok())
                .collect();
            if rgb.len() == 3 {
                asm.colors.insert(name.to_string(), (rgb[0], rgb[1], rgb[2]));
            }
        }
        if parts.len() >= 3 && !parts[2].trim().is_empty() {
            asm.layers.insert(name.to_string(), parts[2].trim().to_string());
        }
    }
    asm
}

/// Write `model` with metadata and read it back, verifying the shape names
/// survive both the embedded block and the STEP representation records.
pub fn model_with_metadata_step_roundtrip(model: &BRepModel) -> Result<(), String> {
    let step = write_step_with_metadata(model);
    let asm = extract_metadata_from_step(&step);
    let got = crate::step::read_step(&step)?;
    for name in model.names() {
        if asm.find(name).is_none() {
            return Err(format!("assembly metadata lost product '{name}'"));
        }
        if !got.names().iter().any(|n| *n == name) {
            return Err(format!("STEP read-back lost shape name '{name}'"));
        }
    }
    Ok(())
}

// ===========================================================================
// XCAF full document schema (XCAFDoc)
// ===========================================================================

/// The four attribute kinds the XCAF schema carries on a shape reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum XcafAttrKind {
    Name,
    Color,
    Layer,
    Material,
}

impl XcafAttrKind {
    /// Stable string tag (`"name"`, `"color"`, `"layer"`, `"material"`).
    pub fn as_str(&self) -> &'static str {
        match self {
            XcafAttrKind::Name => "name",
            XcafAttrKind::Color => "color",
            XcafAttrKind::Layer => "layer",
            XcafAttrKind::Material => "material",
        }
    }

    /// Inverse of [`XcafAttrKind::as_str`].
    pub fn from_str(s: &str) -> Option<Self> {
        match s {
            "name" => Some(XcafAttrKind::Name),
            "color" => Some(XcafAttrKind::Color),
            "layer" => Some(XcafAttrKind::Layer),
            "material" => Some(XcafAttrKind::Material),
            _ => None,
        }
    }
}

/// The attribute set carried by one XCAF shape reference.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct XcafAttrs {
    pub name: Option<String>,
    pub color: Option<(f64, f64, f64)>,
    pub layer: Option<String>,
    pub material: Option<String>,
}

impl XcafAttrs {
    /// Empty attribute set.
    pub fn new() -> Self {
        Self::default()
    }

    /// Attribute set with just a name.
    pub fn named(name: &str) -> Self {
        Self {
            name: Some(name.to_string()),
            ..Default::default()
        }
    }

    /// Get one attribute's value (color formatted as `r,g,b`).
    pub fn get(&self, k: XcafAttrKind) -> Option<String> {
        match k {
            XcafAttrKind::Name => self.name.clone(),
            XcafAttrKind::Color => self.color.map(|(r, g, b)| format!("{r},{g},{b}")),
            XcafAttrKind::Layer => self.layer.clone(),
            XcafAttrKind::Material => self.material.clone(),
        }
    }

    /// Set one attribute from a string value. Color parses `r,g,b`.
    pub fn set(&mut self, k: XcafAttrKind, value: &str) -> Result<(), String> {
        match k {
            XcafAttrKind::Name => self.name = Some(value.to_string()),
            XcafAttrKind::Color => {
                let parts: Vec<f64> = value
                    .split(',')
                    .map(|s| s.trim())
                    .filter_map(|s| s.parse().ok())
                    .collect();
                if parts.len() != 3 {
                    return Err(format!("xcaf: bad color '{value}' (need r,g,b)"));
                }
                self.color = Some((parts[0], parts[1], parts[2]));
            }
            XcafAttrKind::Layer => self.layer = Some(value.to_string()),
            XcafAttrKind::Material => self.material = Some(value.to_string()),
        }
        Ok(())
    }

    /// Whether any attribute is set.
    pub fn is_empty(&self) -> bool {
        self.name.is_none()
            && self.color.is_none()
            && self.layer.is_none()
            && self.material.is_none()
    }
}

/// One node of the XCAF assembly tree: an optional shape, its attributes, and
/// child occurrences. `instance_name` mirrors the STEP
/// `NEXT_ASSEMBLY_USAGE_OCCURRENCE` occurrence name (distinct from the product
/// `name` attribute).
#[derive(Debug, Clone, Default)]
pub struct XcafDocNode {
    pub shape: Option<TopoShape>,
    pub attrs: XcafAttrs,
    pub children: Vec<XcafDocNode>,
    pub instance_name: Option<String>,
}

impl XcafDocNode {
    /// Empty node.
    pub fn new() -> Self {
        Self::default()
    }

    /// Leaf node wrapping `shape`.
    pub fn with_shape(shape: TopoShape) -> Self {
        Self {
            shape: Some(shape),
            ..Default::default()
        }
    }

    /// Append a child occurrence.
    pub fn add_child(&mut self, child: XcafDocNode) -> &mut Self {
        self.children.push(child);
        self
    }

    /// Effective display name: the occurrence name, else the name attribute,
    /// else `""`.
    pub fn name(&self) -> &str {
        self.instance_name
            .as_deref()
            .or(self.attrs.name.as_deref())
            .unwrap_or("")
    }
}

/// A full XCAF document: an assembly tree plus a format version tag.
#[derive(Debug, Clone)]
pub struct XcafDocument {
    pub root: XcafDocNode,
    pub version: String,
}

impl Default for XcafDocument {
    fn default() -> Self {
        Self {
            root: XcafDocNode::default(),
            version: "1.0".to_string(),
        }
    }
}

impl XcafDocument {
    /// An empty document.
    pub fn new() -> Self {
        Self::default()
    }

    /// Total number of nodes in the assembly tree.
    pub fn node_count(&self) -> usize {
        count_nodes(&self.root)
    }
}

/// Total node count of a subtree (including `node` itself).
pub fn count_nodes(node: &XcafDocNode) -> usize {
    1 + node.children.iter().map(count_nodes).sum::<usize>()
}

/// Human-readable tree dump: one line per node with depth indentation and a
/// `[shape]` marker when the node carries geometry.
pub fn document_tree_lines(doc: &XcafDocument) -> Vec<String> {
    let mut out = Vec::new();
    tree_lines(&doc.root, 0, &mut out);
    out
}

fn tree_lines(n: &XcafDocNode, depth: usize, out: &mut Vec<String>) {
    let indent = "  ".repeat(depth);
    let shape = if n.shape.is_some() { " [shape]" } else { "" };
    out.push(format!("{indent}{}{shape}", n.name()));
    for c in &n.children {
        tree_lines(c, depth + 1, out);
    }
}

// ---------------------------------------------------------------------------
// Attribute <-> string-pair conversion
// ---------------------------------------------------------------------------

/// Encode the four attributes as `(kind, value)` string pairs.
///
/// `color` becomes the `"color"` pair with value `r,g,b`; the other three map
/// straight onto their string payloads.
pub fn attrs_to_strings(a: &XcafAttrs) -> Vec<XcafAttribute> {
    let mut out = Vec::new();
    if let Some(n) = &a.name {
        out.push(XcafAttribute {
            kind: "name".into(),
            value: n.clone(),
        });
    }
    if let Some((r, g, b)) = a.color {
        out.push(XcafAttribute {
            kind: "color".into(),
            value: format!("{r},{g},{b}"),
        });
    }
    if let Some(l) = &a.layer {
        out.push(XcafAttribute {
            kind: "layer".into(),
            value: l.clone(),
        });
    }
    if let Some(m) = &a.material {
        out.push(XcafAttribute {
            kind: "material".into(),
            value: m.clone(),
        });
    }
    out
}

/// Decode `(kind, value)` string pairs back into [`XcafAttrs`]. Unknown kinds
/// are ignored; a malformed `color` value is skipped.
pub fn strings_to_attrs(v: &[XcafAttribute]) -> XcafAttrs {
    let mut a = XcafAttrs::default();
    for attr in v {
        match attr.kind.as_str() {
            "name" => a.name = Some(attr.value.clone()),
            "color" => {
                let parts: Vec<f64> = attr
                    .value
                    .split(',')
                    .map(|s| s.trim())
                    .filter_map(|s| s.parse().ok())
                    .collect();
                if parts.len() == 3 {
                    a.color = Some((parts[0], parts[1], parts[2]));
                }
            }
            "layer" => a.layer = Some(attr.value.clone()),
            "material" => a.material = Some(attr.value.clone()),
            _ => {}
        }
    }
    a
}

/// Parse one `kind:value` string (e.g. `"color:1,0,0"`) into an attribute.
pub fn parse_attr_string(s: &str) -> Option<XcafAttribute> {
    let (kind, value) = s.split_once(':')?;
    XcafAttrKind::from_str(kind.trim())?;
    Some(XcafAttribute {
        kind: kind.trim().to_string(),
        value: value.trim().to_string(),
    })
}

/// The reverse of [`parse_attr_string`]: `kind:value` text.
pub fn format_attr_string(a: &XcafAttribute) -> String {
    format!("{}:{}", a.kind, a.value)
}

// ---------------------------------------------------------------------------
// BinXcaf conversion (Phase 6 binary container)
// ---------------------------------------------------------------------------

/// Map a full [`XcafDocument`] onto the binary container. The `instance_name`
/// is carried as a synthetic `"instance_name"` attribute so it survives the
/// binary round-trip.
pub fn xcaf_to_bincaf(doc: &XcafDocument) -> BinXcaf {
    BinXcaf {
        root: node_to_bin_entry(&doc.root),
    }
}

fn node_to_bin_entry(n: &XcafDocNode) -> BinXcafEntry {
    let mut attributes = attrs_to_strings(&n.attrs);
    if let Some(inst) = &n.instance_name {
        attributes.push(XcafAttribute {
            kind: "instance_name".into(),
            value: inst.clone(),
        });
    }
    BinXcafEntry {
        shape: n.shape.clone(),
        attributes,
        children: n.children.iter().map(node_to_bin_entry).collect(),
    }
}

/// Rebuild a [`XcafDocument`] from the binary container.
pub fn bincaf_to_xcaf(b: &BinXcaf) -> XcafDocument {
    XcafDocument {
        root: bin_entry_to_node(&b.root),
        version: "1.0".to_string(),
    }
}

fn bin_entry_to_node(e: &BinXcafEntry) -> XcafDocNode {
    let mut instance_name = None;
    let mut attrs = Vec::new();
    for a in &e.attributes {
        if a.kind == "instance_name" {
            instance_name = Some(a.value.clone());
        } else {
            attrs.push(a.clone());
        }
    }
    XcafDocNode {
        shape: e.shape.clone(),
        attrs: strings_to_attrs(&attrs),
        children: e.children.iter().map(bin_entry_to_node).collect(),
        instance_name,
    }
}

// ---------------------------------------------------------------------------
// XmlXcaf conversion (Phase 7 XML container)
// ---------------------------------------------------------------------------

/// Serialize a document to the XML container format.
pub fn xcaf_to_xml(doc: &XcafDocument) -> String {
    let xdoc = XmlXcafDoc {
        version: doc.version.clone(),
        root: node_to_xml_entry(&doc.root),
    };
    crate::xmlcaf::to_xml(&xdoc).unwrap_or_default()
}

fn node_to_xml_entry(n: &XcafDocNode) -> XmlEntry {
    XmlEntry {
        name: n.instance_name.clone().unwrap_or_default(),
        shape: n.shape.clone(),
        attributes: attrs_to_strings(&n.attrs)
            .iter()
            .map(xcaf_attr_to_xml)
            .collect(),
        children: n.children.iter().map(node_to_xml_entry).collect(),
    }
}

fn xcaf_attr_to_xml(a: &XcafAttribute) -> crate::xmlcaf::XmlAttribute {
    crate::xmlcaf::XmlAttribute {
        kind: a.kind.clone(),
        value: a.value.clone(),
    }
}

fn xml_attr_to_xcaf(a: &crate::xmlcaf::XmlAttribute) -> XcafAttribute {
    XcafAttribute {
        kind: a.kind.clone(),
        value: a.value.clone(),
    }
}

/// Parse a document back from the XML container format.
pub fn xml_to_xcaf(xml: &str) -> Result<XcafDocument, String> {
    let xdoc = crate::xmlcaf::from_xml(xml)?;
    Ok(XcafDocument {
        version: xdoc.version,
        root: xml_entry_to_node(&xdoc.root),
    })
}

fn xml_entry_to_node(e: &XmlEntry) -> XcafDocNode {
    let attrs: Vec<XcafAttribute> = e.attributes.iter().map(xml_attr_to_xcaf).collect();
    XcafDocNode {
        shape: e.shape.clone(),
        attrs: strings_to_attrs(&attrs),
        children: e.children.iter().map(xml_entry_to_node).collect(),
        instance_name: if e.name.is_empty() {
            None
        } else {
            Some(e.name.clone())
        },
    }
}

// ---------------------------------------------------------------------------
// Tree queries
// ---------------------------------------------------------------------------

/// Depth-first search for a node whose name attribute or instance name equals
/// `name`.
pub fn find_by_name<'a>(doc: &'a XcafDocument, name: &str) -> Option<&'a XcafDocNode> {
    find_node(&doc.root, name)
}

fn find_node<'a>(n: &'a XcafDocNode, name: &str) -> Option<&'a XcafDocNode> {
    if n.attrs.name.as_deref() == Some(name) || n.instance_name.as_deref() == Some(name) {
        return Some(n);
    }
    for c in &n.children {
        if let Some(found) = find_node(c, name) {
            return Some(found);
        }
    }
    None
}

/// Flatten the tree into `(shape, attrs)` pairs for every node that carries a
/// shape, in depth-first order.
pub fn apply_attrs_to_shapes(doc: &XcafDocument) -> Vec<(TopoShape, XcafAttrs)> {
    let mut out = Vec::new();
    collect_shapes(&doc.root, &mut out);
    out
}

fn collect_shapes(n: &XcafDocNode, out: &mut Vec<(TopoShape, XcafAttrs)>) {
    if let Some(s) = &n.shape {
        out.push((s.clone(), n.attrs.clone()));
    }
    for c in &n.children {
        collect_shapes(c, out);
    }
}

/// Build a flat [`BRepModel`] from the shapes in a document (one product per
/// shape, named by the `name` attribute).
pub fn model_from_document(doc: &XcafDocument) -> BRepModel {
    let mut model = BRepModel::new();
    let mut i = 0;
    for (shape, attrs) in apply_attrs_to_shapes(doc) {
        let name = attrs.name.clone().unwrap_or_else(|| format!("Part{i}"));
        let idx = model.add(&name, shape);
        if let Some((r, g, b)) = attrs.color {
            model.shapes[idx].color = Some(occt_core::quantity::Color {
                r: r as f32,
                g: g as f32,
                b: b as f32,
            });
        }
        if let Some(l) = &attrs.layer {
            model.shapes[idx].layer = l.clone();
        }
        i += 1;
    }
    model
}

// ---------------------------------------------------------------------------
// Deep STEP read: entity classification, summary, curve lengths
// ---------------------------------------------------------------------------

/// One classified STEP DATA record. Variants carry `(id, …)` where `id` is the
/// `#N` record number as a string.
#[derive(Debug, Clone)]
pub enum StepEntity {
    /// `CARTESIAN_POINT`: (id, [x, y, z]).
    Point((String, [f64; 3])),
    /// `LINE`: (id, name, origin, direction-vector). The direction is the
    /// VECTOR's unit direction scaled by its magnitude.
    Line((String, String, [f64; 3], [f64; 3])),
    /// `CIRCLE`: (id, name, radius).
    Circle((String, String, f64)),
    /// `MANIFOLD_SOLID_BREP`: (id, resolved shape).
    Solid((String, TopoShape)),
    /// `NEXT_ASSEMBLY_USAGE_OCCURRENCE`: (id, referenced entity ids).
    Assembly((String, Vec<String>)),
}

/// A parsed `#N=TYPE(...)` data record (lightweight, STEP-text local copy).
#[derive(Debug, Clone)]
struct StepRecord {
    type_name: String,
    args: Vec<String>,
}

impl StepRecord {
    /// Reconstruct the original entity body `TYPE(a,b,...)` from the args.
    fn body(&self) -> String {
        format!("{}({})", self.type_name, self.args.join(","))
    }
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
fn parse_entity_body(body: &str) -> (String, Vec<String>) {
    let body = body.trim();
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

/// Scan one entity body starting just after `=`, stopping at the terminating
/// `;`. Returns the body text and the index just past the `;`.
fn scan_entity_body(data: &str, start: usize) -> Result<(String, usize), String> {
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
                        return Err("STEP: unbalanced parens in entity body".into());
                    }
                    depth -= 1;
                    if depth == 0 {
                        let body = data[start..=i].to_string();
                        let mut j = i + 1;
                        while j < bytes.len() && bytes[j].is_ascii_whitespace() {
                            j += 1;
                        }
                        if j >= bytes.len() || bytes[j] != b';' {
                            return Err("STEP: missing ';' after entity body".into());
                        }
                        return Ok((body, j + 1));
                    }
                }
                _ => {}
            }
        }
        i += 1;
    }
    Err("STEP: unterminated entity body".into())
}

/// Split the DATA section into `#id=TYPE(...)` records.
fn scan_data_records(data: &str) -> Result<Vec<(usize, StepRecord)>, String> {
    let bytes = data.as_bytes();
    let mut out = Vec::new();
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
                .map_err(|_| format!("STEP: bad entity id at offset {i}"))?;
            let mut k = j;
            while k < bytes.len() && bytes[k].is_ascii_whitespace() {
                k += 1;
            }
            if k >= bytes.len() || bytes[k] != b'=' {
                return Err(format!("STEP: malformed record #{id}: expected '='"));
            }
            let (body, next) = scan_entity_body(data, k + 1)?;
            let (type_name, args) = parse_entity_body(&body);
            out.push((id, StepRecord { type_name, args }));
            i = next;
        } else {
            i += 1;
        }
    }
    Ok(out)
}

/// Parse a whole STEP physical file's DATA section into an id→record map.
fn parse_step_data(content: &str) -> Result<HashMap<usize, StepRecord>, String> {
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
    let data = &after_data[..data_end];
    let mut map = HashMap::new();
    for (id, rec) in scan_data_records(data)? {
        map.insert(id, rec);
    }
    Ok(map)
}

fn parse_ref(s: &str) -> Option<usize> {
    s.trim().strip_prefix('#')?.trim().parse().ok()
}

fn parse_str(s: &str) -> String {
    let s = s.trim();
    if s.len() >= 2 && s.starts_with('\'') && s.ends_with('\'') {
        s[1..s.len() - 1].replace("''", "'")
    } else {
        s.to_string()
    }
}

/// Parse a `(x,y,z)` tuple into `[x, y, z]`.
fn parse_xyz3(s: &str) -> Option<[f64; 3]> {
    let s = s.trim();
    let inner = s.trim_start_matches('(').trim_end_matches(')');
    let parts: Vec<String> = split_top(inner).into_iter().map(|p| p.trim().to_string()).collect();
    if parts.len() != 3 {
        return None;
    }
    let x: f64 = parts[0].parse().ok()?;
    let y: f64 = parts[1].parse().ok()?;
    let z: f64 = parts[2].parse().ok()?;
    Some([x, y, z])
}

/// All `#N` references appearing anywhere in `s`.
fn extract_refs(s: &str) -> Vec<usize> {
    let bytes = s.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'#' {
            let mut j = i + 1;
            while j < bytes.len() && bytes[j].is_ascii_digit() {
                j += 1;
            }
            if j > i + 1 {
                if let Ok(n) = s[i + 1..j].parse() {
                    out.push(n);
                }
            }
            i = j;
        } else {
            i += 1;
        }
    }
    out
}

/// Classify the DATA records of a STEP text into [`StepEntity`] values.
///
/// Points / lines / circles are decoded directly from their records; each
/// `MANIFOLD_SOLID_BREP` is resolved into a real [`TopoShape`] by re-feeding
/// its transitive record closure to the existing STEP reader; each
/// `NEXT_ASSEMBLY_USAGE_OCCURRENCE` is kept as an (id, refs) pair.
pub fn step_entities_from_text(content: &str) -> Result<Vec<StepEntity>, String> {
    let records = parse_step_data(content)?;

    // Pass 1: geometry lookups shared by LINE / VECTOR / DIRECTION.
    let mut points: HashMap<usize, [f64; 3]> = HashMap::new();
    let mut dirs: HashMap<usize, [f64; 3]> = HashMap::new();
    let mut vecs: HashMap<usize, (usize, f64)> = HashMap::new(); // (dir_ref, magnitude)
    for (&id, rec) in &records {
        match rec.type_name.as_str() {
            "CARTESIAN_POINT" => {
                if let Some(xyz) = rec.args.get(1).and_then(|a| parse_xyz3(a)) {
                    points.insert(id, xyz);
                }
            }
            "DIRECTION" => {
                if let Some(xyz) = rec.args.get(1).and_then(|a| parse_xyz3(a)) {
                    dirs.insert(id, xyz);
                }
            }
            "VECTOR" => {
                let dr = rec.args.get(1).and_then(|a| parse_ref(a));
                let mag = rec.args.get(2).and_then(|a| a.trim().parse::<f64>().ok());
                if let (Some(dr), Some(mag)) = (dr, mag) {
                    vecs.insert(id, (dr, mag));
                }
            }
            _ => {}
        }
    }

    let mut entities = Vec::new();
    let mut ids: Vec<usize> = records.keys().cloned().collect();
    ids.sort_unstable();
    for id in ids {
        let rec = &records[&id];
        match rec.type_name.as_str() {
            "CARTESIAN_POINT" => {
                if let Some(p) = points.get(&id) {
                    entities.push(StepEntity::Point((id.to_string(), *p)));
                }
            }
            "LINE" => {
                let name = rec.args.first().map(|a| parse_str(a)).unwrap_or_default();
                let pnt_ref = rec.args.get(1).and_then(|a| parse_ref(a));
                let vec_ref = rec.args.get(2).and_then(|a| parse_ref(a));
                if let (Some(pr), Some(vr)) = (pnt_ref, vec_ref) {
                    if let (Some(p), Some((dr, mag))) = (points.get(&pr), vecs.get(&vr)) {
                        let d = dirs.get(dr).cloned().unwrap_or([0.0, 0.0, 1.0]);
                        let dir = [d[0] * mag, d[1] * mag, d[2] * mag];
                        entities.push(StepEntity::Line((id.to_string(), name, *p, dir)));
                    }
                }
            }
            "CIRCLE" => {
                let name = rec.args.first().map(|a| parse_str(a)).unwrap_or_default();
                if let Some(r) = rec.args.get(2).and_then(|a| a.trim().parse::<f64>().ok()) {
                    entities.push(StepEntity::Circle((id.to_string(), name, r)));
                }
            }
            "MANIFOLD_SOLID_BREP" => {
                if let Some(shape) = resolve_solid_from_records(&records, id) {
                    entities.push(StepEntity::Solid((id.to_string(), shape)));
                }
            }
            "NEXT_ASSEMBLY_USAGE_OCCURRENCE" => {
                let refs: Vec<String> = rec
                    .args
                    .iter()
                    .filter_map(|a| parse_ref(a).map(|r| format!("#{r}")))
                    .collect();
                entities.push(StepEntity::Assembly((id.to_string(), refs)));
            }
            _ => {}
        }
    }
    Ok(entities)
}

/// Resolve one `MANIFOLD_SOLID_BREP` record into a real [`TopoShape`] by
/// collecting its transitive record closure into a synthetic STEP file and
/// reusing the existing STEP reader's full topological resolver.
fn resolve_solid_from_records(
    records: &HashMap<usize, StepRecord>,
    seed: usize,
) -> Option<TopoShape> {
    let mut reachable: HashSet<usize> = HashSet::new();
    let mut stack = vec![seed];
    while let Some(id) = stack.pop() {
        if !reachable.insert(id) {
            continue;
        }
        if let Some(rec) = records.get(&id) {
            for r in extract_refs(&rec.body()) {
                if !reachable.contains(&r) {
                    stack.push(r);
                }
            }
        }
    }
    let mut sorted: Vec<usize> = reachable.into_iter().collect();
    sorted.sort_unstable();
    let mut lines = String::new();
    for id in sorted {
        if let Some(rec) = records.get(&id) {
            lines.push_str(&format!("#{id}={};\n", rec.body()));
        }
    }
    let synthetic = format!("ISO-10303-21;\nDATA;\n{lines}ENDSEC;\nEND-ISO-10303-21;\n");
    crate::step::read_step(&synthetic)
        .ok()
        .and_then(|m| m.shapes.into_iter().next().map(|s| s.shape))
}

/// Read a STEP physical file from `path` and classify its DATA entities.
pub fn read_step_entities(path: &str) -> Result<Vec<StepEntity>, String> {
    let content = std::fs::read_to_string(path).map_err(|e| format!("STEP read {path}: {e}"))?;
    step_entities_from_text(&content)
}

/// A text summary of a STEP file: total entity count and per-type counts.
pub fn read_step_summary(path: &str) -> Result<String, String> {
    let content = std::fs::read_to_string(path).map_err(|e| format!("STEP read {path}: {e}"))?;
    let records = parse_step_data(&content)?;
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    for rec in records.values() {
        *counts.entry(rec.type_name.clone()).or_insert(0) += 1;
    }
    let total: usize = counts.values().sum();
    let mut lines = vec![format!("STEP entities: {total}")];
    for (t, c) in counts {
        lines.push(format!("{t}: {c}"));
    }
    Ok(lines.join("\n"))
}

/// For every `LINE` / `CIRCLE` entity in a STEP file, its length: the VECTOR
/// magnitude for a line and `2πr` for a circle. Returns `(id, length)` pairs.
pub fn step_curve_lengths(path: &str) -> Result<Vec<(String, f64)>, String> {
    let content = std::fs::read_to_string(path).map_err(|e| format!("STEP read {path}: {e}"))?;
    let records = parse_step_data(&content)?;

    let mut vecs: HashMap<usize, (usize, f64)> = HashMap::new(); // (dir_ref, magnitude)
    for (&id, rec) in &records {
        if rec.type_name == "VECTOR" {
            let dr = rec.args.get(1).and_then(|a| parse_ref(a));
            let mag = rec.args.get(2).and_then(|a| a.trim().parse::<f64>().ok());
            if let (Some(dr), Some(mag)) = (dr, mag) {
                vecs.insert(id, (dr, mag));
            }
        }
    }

    let mut out = Vec::new();
    let mut ids: Vec<usize> = records.keys().cloned().collect();
    ids.sort_unstable();
    for id in ids {
        let rec = &records[&id];
        match rec.type_name.as_str() {
            "LINE" => {
                if let Some(vr) = rec.args.get(2).and_then(|a| parse_ref(a)) {
                    if let Some((_, mag)) = vecs.get(&vr) {
                        out.push((id.to_string(), *mag));
                    }
                }
            }
            "CIRCLE" => {
                if let Some(r) = rec.args.get(2).and_then(|a| a.trim().parse::<f64>().ok()) {
                    out.push((id.to_string(), 2.0 * PI * r));
                }
            }
            _ => {}
        }
    }
    Ok(out)
}

// ---------------------------------------------------------------------------
// Document → STEP assembly writer
// ---------------------------------------------------------------------------

/// Serialize a [`XcafDocument`] to a STEP physical file with a simple assembly
/// structure: each shape becomes a named `ADVANCED_BREP_SHAPE_REPRESENTATION`
/// product (via the existing STEP writer) and every parent→child occurrence in
/// the tree becomes a `NEXT_ASSEMBLY_USAGE_OCCURRENCE`. Attributes that STEP
/// has no native slot for (layer, material, instance names) are embedded as a
/// `/* XCAFDOC-BEGIN … */` comment block before `DATA`, keeping the file valid.
pub fn document_to_step(doc: &XcafDocument) -> String {
    let step = crate::step::write_step(&model_from_document(doc));

    // Assembly edges parent→child, from the document tree.
    let mut edges = Vec::new();
    collect_edges(&doc.root, None, &mut edges);
    let nauo = assembly_nauo_records(&step, &edges);

    let mut block = String::new();
    block.push_str("/* XCAFDOC-BEGIN\n");
    write_tree_block(&mut block, &doc.root, 0);
    block.push_str("XCAFDOC-END\n*/\n");

    let with_block = step.replacen("DATA;", &format!("{block}DATA;"), 1);
    if nauo.is_empty() {
        return with_block;
    }
    match with_block.find("DATA;") {
        Some(pos) => {
            let (head, tail) = with_block.split_at(pos + 5);
            format!("{head}\n{nauo}{tail}")
        }
        None => with_block,
    }
}

/// Emit `NEXT_ASSEMBLY_USAGE_OCCURRENCE` records for every tree edge whose
/// endpoints both have product definitions in `step_text`.
fn assembly_nauo_records(step_text: &str, edges: &[(String, String)]) -> String {
    let Ok(records) = parse_step_data(step_text) else {
        return String::new();
    };
    let pd = product_definition_by_name(&records);
    let max_id = records.keys().max().copied().unwrap_or(0);
    let mut out = String::new();
    let mut id = max_id + 1;
    for (parent, child) in edges {
        let (Some(&prel), Some(&crel)) = (pd.get(parent), pd.get(child)) else {
            continue;
        };
        let name = child.replace('\'', "''");
        out.push_str(&format!(
            "#{id}=NEXT_ASSEMBLY_USAGE_OCCURRENCE('{name}',#{prel},#{crel},$);\n"
        ));
        id += 1;
    }
    out
}

/// Map product name → product-definition record id from a parsed STEP file.
fn product_definition_by_name(
    records: &HashMap<usize, StepRecord>,
) -> HashMap<String, usize> {
    let mut prod_names: HashMap<usize, String> = HashMap::new();
    for (&id, rec) in records {
        if rec.type_name == "PRODUCT" {
            prod_names.insert(id, parse_str(&rec.args.first().cloned().unwrap_or_default()));
        }
    }
    let mut formation_prod: HashMap<usize, usize> = HashMap::new();
    for (&id, rec) in records {
        if rec.type_name == "PRODUCT_DEFINITION_FORMATION" {
            if let Some(p) = rec.args.get(2).and_then(|a| parse_ref(a)) {
                formation_prod.insert(id, p);
            }
        }
    }
    let mut by_name = HashMap::new();
    for (&id, rec) in records {
        if rec.type_name == "PRODUCT_DEFINITION" {
            if let Some(f) = rec.args.get(3).and_then(|a| parse_ref(a)) {
                if let Some(p) = formation_prod.get(&f) {
                    if let Some(name) = prod_names.get(p) {
                        by_name.insert(name.clone(), id);
                    }
                }
            }
        }
    }
    by_name
}

/// Depth-first (parent, child) name pairs for every tree edge.
fn collect_edges(n: &XcafDocNode, parent: Option<&str>, out: &mut Vec<(String, String)>) {
    let name = n.name().to_string();
    if let Some(p) = parent {
        out.push((p.to_string(), name.clone()));
    }
    for c in &n.children {
        collect_edges(c, Some(&name), out);
    }
}

/// One line per node of the assembly tree, indented by depth:
/// `instance|name|color|layer|material`.
fn write_tree_block(out: &mut String, n: &XcafDocNode, depth: usize) {
    let inst = sanitize(n.instance_name.as_deref().unwrap_or(""));
    let name = sanitize(n.attrs.name.as_deref().unwrap_or(""));
    let color = match n.attrs.color {
        Some((r, g, b)) => format!("{r},{g},{b}"),
        None => String::new(),
    };
    let layer = sanitize(n.attrs.layer.as_deref().unwrap_or(""));
    let material = sanitize(n.attrs.material.as_deref().unwrap_or(""));
    let indent = "  ".repeat(depth);
    out.push_str(&format!(
        "{indent}{inst}|{name}|{color}|{layer}|{material}\n"
    ));
    for c in &n.children {
        write_tree_block(out, c, depth + 1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::{BRepPrimBox, BRepPrimSphere};
    use occt_core::quantity::Color;

    fn two_shape_model() -> BRepModel {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let s = BRepPrimSphere::make_sphere(0.5);
        let mut model = BRepModel::new();
        model.add_with_color("RedBox", b.solid.0.clone(), Color::RED);
        model.add("GreyBall", s.solid.0.clone());
        model.find_mut("GreyBall").unwrap().layer = "L1".into();
        model
    }

    #[test]
    fn model_to_assembly_has_both_products() {
        let model = two_shape_model();
        let asm = model_to_step_assembly(&model);
        assert!(asm.find("RedBox").is_some());
        assert!(asm.find("GreyBall").is_some());
        assert_eq!(asm.colors.get("RedBox"), Some(&(1.0, 0.0, 0.0)));
        assert_eq!(asm.layers.get("GreyBall").map(|s| s.as_str()), Some("L1"));
        assert_eq!(asm.products.len(), 2);
    }

    #[test]
    fn write_contains_names() {
        let model = two_shape_model();
        let step = write_step_with_metadata(&model);
        assert!(step.contains("RedBox"), "step:\n{step}");
        assert!(step.contains("GreyBall"));
        assert!(step.contains(BLOCK_BEGIN));
    }

    #[test]
    fn extract_roundtrips_names_and_colors() {
        let model = two_shape_model();
        let step = write_step_with_metadata(&model);
        let asm = extract_metadata_from_step(&step);
        assert!(asm.find("RedBox").is_some());
        assert!(asm.find("GreyBall").is_some());
        assert_eq!(asm.colors.get("RedBox"), Some(&(1.0, 0.0, 0.0)));
        assert_eq!(asm.layers.get("GreyBall").map(|s| s.as_str()), Some("L1"));
        assert_eq!(asm.products.len(), 2);
    }

    #[test]
    fn full_roundtrip_succeeds() {
        let model = two_shape_model();
        model_with_metadata_step_roundtrip(&model).expect("roundtrip ok");
    }
}

#[cfg(test)]
mod xcaf_doc_tests {
    use super::*;
    use crate::abs::ShapeType;
    use crate::primitives::{BRepPrimBox, BRepPrimSphere};
    use crate::topo_tools_full::faces_of;
    use std::sync::atomic::{AtomicUsize, Ordering};

    static TMP_SEQ: AtomicUsize = AtomicUsize::new(0);

    fn write_temp_step(content: &str) -> String {
        let n = TMP_SEQ.fetch_add(1, Ordering::SeqCst);
        let path = std::env::temp_dir().join(format!(
            "occt_xcaf_doc_{}_{n}.step",
            std::process::id()
        ));
        std::fs::write(&path, content).expect("write temp step");
        path.to_str().unwrap().to_string()
    }

    fn sample_doc() -> XcafDocument {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let mut leaf = XcafDocNode::with_shape(b.solid.0);
        leaf.attrs = XcafAttrs {
            name: Some("Leaf".into()),
            color: Some((0.0, 1.0, 0.0)),
            layer: Some("L1".into()),
            material: Some("steel".into()),
        };
        leaf.instance_name = Some("leaf-01".into());
        let mut mid = XcafDocNode::new();
        mid.attrs = XcafAttrs::named("Mid");
        mid.children.push(leaf);
        let mut root = XcafDocNode::new();
        root.attrs = XcafAttrs::named("Root");
        root.children.push(mid);
        XcafDocument {
            root,
            version: "1.0".into(),
        }
    }

    #[test]
    fn xcaf_to_bincaf_roundtrip() {
        let doc = sample_doc();
        let bin = xcaf_to_bincaf(&doc);
        let bytes = crate::bincaf::serialize_bincaf(&bin).expect("serialize");
        let bin2 = crate::bincaf::deserialize_bincaf(&bytes).expect("deserialize");
        let got = bincaf_to_xcaf(&bin2);
        assert_eq!(got.root.attrs.name.as_deref(), Some("Root"));
        assert_eq!(got.root.children[0].attrs.name.as_deref(), Some("Mid"));
        let leaf = &got.root.children[0].children[0];
        assert_eq!(leaf.attrs.name.as_deref(), Some("Leaf"));
        assert_eq!(leaf.attrs.color, Some((0.0, 1.0, 0.0)));
        assert_eq!(leaf.attrs.layer.as_deref(), Some("L1"));
        assert_eq!(leaf.attrs.material.as_deref(), Some("steel"));
        assert_eq!(leaf.instance_name.as_deref(), Some("leaf-01"));
        assert_eq!(
            leaf.shape.as_ref().map(|s| s.shape_type()),
            Some(ShapeType::Solid)
        );
    }

    #[test]
    fn attrs_strings_roundtrip() {
        let a = XcafAttrs {
            name: Some("Widget".into()),
            color: Some((1.0, 0.0, 0.0)),
            layer: Some("L1".into()),
            material: Some("steel".into()),
        };
        let s = attrs_to_strings(&a);
        assert_eq!(s.len(), 4);
        let b = strings_to_attrs(&s);
        assert_eq!(a, b);
        // kind:value text form.
        let p = parse_attr_string("color:1,0,0").expect("parse attr");
        assert_eq!((p.kind.as_str(), p.value.as_str()), ("color", "1,0,0"));
        assert_eq!(format_attr_string(&p), "color:1,0,0");
    }

    #[test]
    fn xcaf_xml_roundtrip() {
        let doc = sample_doc();
        let xml = xcaf_to_xml(&doc);
        let got = xml_to_xcaf(&xml).expect("parse xml");
        assert_eq!(got.version, "1.0");
        assert_eq!(got.root.attrs.name.as_deref(), Some("Root"));
        assert_eq!(got.root.children.len(), 1);
        let leaf = &got.root.children[0].children[0];
        assert_eq!(leaf.attrs.name.as_deref(), Some("Leaf"));
        assert_eq!(leaf.attrs.material.as_deref(), Some("steel"));
        assert_eq!(leaf.instance_name.as_deref(), Some("leaf-01"));
        assert_eq!(
            leaf.shape.as_ref().map(|s| s.shape_type()),
            Some(ShapeType::Solid)
        );
    }

    #[test]
    fn find_by_name_nested() {
        let doc = sample_doc();
        assert!(find_by_name(&doc, "Root").is_some());
        assert!(find_by_name(&doc, "Mid").is_some());
        let leaf = find_by_name(&doc, "Leaf").expect("leaf found");
        assert_eq!(leaf.attrs.name.as_deref(), Some("Leaf"));
        assert!(find_by_name(&doc, "Nope").is_none());
        assert_eq!(count_nodes(&doc.root), 3);
        let lines = document_tree_lines(&doc);
        assert_eq!(lines.len(), 3);
    }

    #[test]
    fn apply_attrs_flattens() {
        let doc = sample_doc();
        let flat = apply_attrs_to_shapes(&doc);
        assert_eq!(flat.len(), 1, "only the leaf carries a shape");
        assert_eq!(flat[0].1.name.as_deref(), Some("Leaf"));
        assert_eq!(flat[0].1.color, Some((0.0, 1.0, 0.0)));
        assert_eq!(flat[0].0.shape_type(), ShapeType::Solid);
        assert_eq!(faces_of(&flat[0].0).len(), 6);
    }

    #[test]
    fn read_step_points_lines() {
        let step = "ISO-10303-21;\nDATA;\n\
#1=CARTESIAN_POINT('P1',(1.,2.,3.));\n\
#2=DIRECTION('',(0.,0.,1.));\n\
#3=VECTOR('',#2,5.);\n\
#4=LINE('LineA',#1,#3);\n\
ENDSEC;\nEND-ISO-10303-21;";
        let path = write_temp_step(step);
        let entities = read_step_entities(&path).expect("read entities");
        let lengths = step_curve_lengths(&path).expect("curve lengths");
        std::fs::remove_file(&path).ok();

        let mut n_points = 0;
        let mut n_lines = 0;
        for e in entities {
            match e {
                StepEntity::Point((id, p)) => {
                    assert_eq!(id, "1");
                    assert_eq!(p, [1., 2., 3.]);
                    n_points += 1;
                }
                StepEntity::Line((id, name, origin, dir)) => {
                    assert_eq!(id, "4");
                    assert_eq!(name, "LineA");
                    assert_eq!(origin, [1., 2., 3.]);
                    assert_eq!(dir, [0., 0., 5.]);
                    n_lines += 1;
                }
                _ => {}
            }
        }
        assert_eq!(n_points, 1);
        assert_eq!(n_lines, 1);
        // LINE length = VECTOR magnitude.
        assert_eq!(lengths, vec![("4".to_string(), 5.0)]);
    }

    #[test]
    fn read_step_circle_length() {
        let step = "ISO-10303-21;\nDATA;\n\
#1=CARTESIAN_POINT('',(0.,0.,0.));\n\
#2=DIRECTION('',(0.,0.,1.));\n\
#3=DIRECTION('',(1.,0.,0.));\n\
#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);\n\
#5=CIRCLE('',#4,2.);\n\
ENDSEC;\nEND-ISO-10303-21;";
        let path = write_temp_step(step);
        let lengths = step_curve_lengths(&path).expect("lengths");
        std::fs::remove_file(&path).ok();
        assert_eq!(lengths.len(), 1);
        assert_eq!(lengths[0].0, "5");
        assert!(
            (lengths[0].1 - 2.0 * PI * 2.0).abs() < 1e-9,
            "length {}",
            lengths[0].1
        );
    }

    #[test]
    fn read_step_solid_present() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let step = crate::step::write_shape_step(&b.solid.0);
        let path = write_temp_step(&step);
        let entities = read_step_entities(&path).expect("read entities");
        std::fs::remove_file(&path).ok();
        let mut solids = 0;
        for e in &entities {
            if let StepEntity::Solid((_, shape)) = e {
                assert_eq!(shape.shape_type(), ShapeType::Solid);
                assert_eq!(faces_of(shape).len(), 6);
                solids += 1;
            }
        }
        assert_eq!(solids, 1, "one solid expected");
    }

    #[test]
    fn read_step_summary_types() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let step = crate::step::write_shape_step(&b.solid.0);
        let path = write_temp_step(&step);
        let summary = read_step_summary(&path).expect("summary");
        std::fs::remove_file(&path).ok();
        assert!(summary.contains("MANIFOLD_SOLID_BREP: 1"), "{summary}");
        assert!(summary.contains("ADVANCED_FACE"), "{summary}");
        assert!(summary.contains("CARTESIAN_POINT"), "{summary}");
    }

    #[test]
    fn document_to_step_writes() {
        let b = BRepPrimBox::make_box(1.0, 1.0, 1.0);
        let s = BRepPrimSphere::make_sphere(0.5);
        let mut part = XcafDocNode::with_shape(b.solid.0);
        part.attrs = XcafAttrs {
            name: Some("CubePart".into()),
            layer: Some("L2".into()),
            ..Default::default()
        };
        part.instance_name = Some("cube-01".into());
        let mut root = XcafDocNode::with_shape(s.solid.0);
        root.attrs = XcafAttrs::named("AssemblyRoot");
        root.children.push(part);
        let doc = XcafDocument {
            root,
            version: "1.0".into(),
        };
        let step = document_to_step(&doc);
        assert!(step.contains("CubePart"), "product name:\n{step}");
        assert!(step.contains("AssemblyRoot"), "assembly name");
        assert!(step.contains("L2"), "layer");
        assert!(step.contains("cube-01"), "instance name");
        assert!(step.contains("XCAFDOC-BEGIN"), "tree block");
    }

    #[test]
    fn document_to_step_roundtrip() {
        let b1 = BRepPrimBox::make_box(2.0, 2.0, 2.0);
        let b2 = BRepPrimSphere::make_sphere(1.0);
        let mut part = XcafDocNode::with_shape(b1.solid.0);
        part.attrs = XcafAttrs::named("BoxPart");
        let mut root = XcafDocNode::with_shape(b2.solid.0);
        root.attrs = XcafAttrs::named("Assembly");
        root.children.push(part);
        let doc = XcafDocument {
            root,
            version: "1.0".into(),
        };
        let step = document_to_step(&doc);
        assert!(
            step.contains("NEXT_ASSEMBLY_USAGE_OCCURRENCE"),
            "assembly:\n{step}"
        );
        let model = crate::step::read_step(&step).expect("read back");
        assert_eq!(model.len(), 2, "two solids");
        let names: Vec<&str> = model.names();
        assert!(names.contains(&"Assembly"), "names: {names:?}");
        assert!(names.contains(&"BoxPart"), "names: {names:?}");
        for ms in &model.shapes {
            assert_eq!(ms.shape.shape_type(), ShapeType::Solid);
        }
    }

    #[test]
    fn empty_doc_ok() {
        let doc = XcafDocument::new();
        let bin = xcaf_to_bincaf(&doc);
        let got = bincaf_to_xcaf(&bin);
        assert!(got.root.attrs.name.is_none());
        assert!(got.root.children.is_empty());
        let xml = xcaf_to_xml(&doc);
        let got2 = xml_to_xcaf(&xml).expect("xml parse");
        assert!(got2.root.children.is_empty());
        let step = document_to_step(&doc);
        assert!(step.contains("DATA;"), "{step}");
        assert!(strings_to_attrs(&[]).is_empty());
    }
}
