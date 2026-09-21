use super::prelude::*;
use super::*;

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
pub(super) const BLOCK_BEGIN: &str = "/* ASSEMBLY-BEGIN";
pub(super) const BLOCK_END: &str = "ASSEMBLY-END";

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
pub(super) fn format2(c: occt_core::quantity::Color) -> String {
    let r = (c.r * 1000.0).round() / 1000.0;
    let g = (c.g * 1000.0).round() / 1000.0;
    let b = (c.b * 1000.0).round() / 1000.0;
    format!("{r},{g},{b}")
}

/// Replace characters that would break the `|`-delimited comment format.
pub(super) fn sanitize(s: &str) -> String {
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
/// `name` attribute). `location` is the occurrence placement relative to its
/// parent (the `TopLoc_Location` of a `XCAFDoc` component).
#[derive(Debug, Clone, Default)]
pub struct XcafDocNode {
    pub shape: Option<TopoShape>,
    pub attrs: XcafAttrs,
    pub children: Vec<XcafDocNode>,
    pub instance_name: Option<String>,
    /// Placement of this occurrence in its parent's frame; `None` = identity.
    pub location: Option<GpTrsf>,
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

/// A full XCAF document: an assembly tree, named views, dimension/note
/// annotations, plus a format version tag.
#[derive(Debug, Clone)]
pub struct XcafDocument {
    pub root: XcafDocNode,
    pub version: String,
    /// Named camera views (`XCAFDoc_DocumentTool` view layer).
    pub views: Vec<XcafView>,
    /// Dimension / note / callout annotations.
    pub annotations: Vec<XcafAnnotation>,
}

impl Default for XcafDocument {
    fn default() -> Self {
        Self {
            root: XcafDocNode::default(),
            version: "1.0".to_string(),
            views: Vec::new(),
            annotations: Vec::new(),
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

// ---------------------------------------------------------------------------
// Named views + annotations (Phase 12 XCAF schema)
// ---------------------------------------------------------------------------

/// Projection type of a named camera view.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum XcafProjection {
    /// Perspective camera (foreshortened).
    Perspective,
    /// Orthographic camera (parallel projection).
    Orthographic,
}

/// A named camera view of the document. Mirrors the `XCAFDoc_*` view layer
/// (camera position / look-at target / up direction / projection).
#[derive(Debug, Clone, PartialEq)]
pub struct XcafView {
    pub name: String,
    /// Camera (eye) position in world coordinates.
    pub eye: GpPnt,
    /// Look-at target point.
    pub target: GpPnt,
    /// Up direction (unit).
    pub up: GpDir,
    pub projection: XcafProjection,
}

/// The three annotation kinds carried by [`XcafAnnotation`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum XcafAnnotationKind {
    Dimension,
    Note,
    Callout,
}

impl XcafAnnotationKind {
    /// Stable string tag (`"dimension"`, `"note"`, `"callout"`).
    pub fn as_str(&self) -> &'static str {
        match self {
            XcafAnnotationKind::Dimension => "dimension",
            XcafAnnotationKind::Note => "note",
            XcafAnnotationKind::Callout => "callout",
        }
    }

    /// Inverse of [`XcafAnnotationKind::as_str`].
    pub fn from_str(s: &str) -> Option<Self> {
        match s {
            "dimension" => Some(XcafAnnotationKind::Dimension),
            "note" => Some(XcafAnnotationKind::Note),
            "callout" => Some(XcafAnnotationKind::Callout),
            _ => None,
        }
    }
}

/// A dimension / note / callout annotation anchored at a world point.
#[derive(Debug, Clone, PartialEq)]
pub struct XcafAnnotation {
    /// Stable identifier used for lookup / removal.
    pub id: String,
    pub kind: XcafAnnotationKind,
    /// Anchor point of the annotation in world coordinates.
    pub anchor: GpPnt,
    /// Display text.
    pub text: String,
    /// Measured value (e.g. a dimension length).
    pub value: f64,
}

/// A read view of one placed occurrence in the assembly tree: its name, the
/// referenced sub-tree and the placement relative to its parent.
#[derive(Debug, Clone)]
pub struct XcafInstance {
    pub name: String,
    pub node: XcafDocNode,
    pub transform: GpTrsf,
}

// ---------------------------------------------------------------------------
// View / annotation / transform string encoding
//
// Views, annotations and node placements are carried through the binary / XML
// containers as synthetic string attributes (`__view`, `__annotation`,
// `location`), which the generic `strings_to_attrs` decoder already ignores.
// ---------------------------------------------------------------------------

/// Shortest round-trippable decimal form of a float.
pub(super) fn fmt(v: f64) -> String {
    format!("{v}")
}

/// `name|eye.xyz|target.xyz|up.xyz|projection` (name sanitized of `|`).
pub(super) fn view_to_string(v: &XcafView) -> String {
    let proj = match v.projection {
        XcafProjection::Perspective => "perspective",
        XcafProjection::Orthographic => "orthographic",
    };
    format!(
        "{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}",
        sanitize(&v.name),
        fmt(v.eye.x()),
        fmt(v.eye.y()),
        fmt(v.eye.z()),
        fmt(v.target.x()),
        fmt(v.target.y()),
        fmt(v.target.z()),
        fmt(v.up.x()),
        fmt(v.up.y()),
        fmt(v.up.z()),
        proj
    )
}

pub(super) fn view_from_string(s: &str) -> Option<XcafView> {
    let parts: Vec<&str> = s.split('|').collect();
    if parts.len() != 11 {
        return None;
    }
    let n: Vec<f64> = parts[1..10].iter().filter_map(|p| p.parse().ok()).collect();
    if n.len() != 9 {
        return None;
    }
    Some(XcafView {
        name: parts[0].to_string(),
        eye: GpPnt::new(n[0], n[1], n[2]),
        target: GpPnt::new(n[3], n[4], n[5]),
        up: GpDir::new(n[6], n[7], n[8]).ok()?,
        projection: if parts[10] == "orthographic" {
            XcafProjection::Orthographic
        } else {
            XcafProjection::Perspective
        },
    })
}

/// `id|kind|anchor.xyz|value|text` (id/text sanitized of `|`).
pub(super) fn annotation_to_string(a: &XcafAnnotation) -> String {
    format!(
        "{}|{}|{}|{}|{}|{}|{}",
        sanitize(&a.id),
        a.kind.as_str(),
        fmt(a.anchor.x()),
        fmt(a.anchor.y()),
        fmt(a.anchor.z()),
        fmt(a.value),
        sanitize(&a.text)
    )
}

pub(super) fn annotation_from_string(s: &str) -> Option<XcafAnnotation> {
    let parts: Vec<&str> = s.split('|').collect();
    if parts.len() != 7 {
        return None;
    }
    let kind = XcafAnnotationKind::from_str(parts[1])?;
    let x: f64 = parts[2].parse().ok()?;
    let y: f64 = parts[3].parse().ok()?;
    let z: f64 = parts[4].parse().ok()?;
    let value: f64 = parts[5].parse().ok()?;
    Some(XcafAnnotation {
        id: parts[0].to_string(),
        kind,
        anchor: GpPnt::new(x, y, z),
        text: parts[6].to_string(),
        value,
    })
}

/// `scale|m00|m01|...|m22|lx|ly|lz` — the full `gp_Trsf` linear+translation
/// part, so any placement round-trips exactly.
pub(super) fn transform_to_string(t: &GpTrsf) -> String {
    let m = t.matrix.m;
    format!(
        "{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}",
        fmt(t.scale),
        fmt(m[0][0]),
        fmt(m[0][1]),
        fmt(m[0][2]),
        fmt(m[1][0]),
        fmt(m[1][1]),
        fmt(m[1][2]),
        fmt(m[2][0]),
        fmt(m[2][1]),
        fmt(m[2][2]),
        fmt(t.loc.x),
        fmt(t.loc.y),
        fmt(t.loc.z)
    )
}

pub(super) fn transform_from_string(s: &str) -> Option<GpTrsf> {
    let parts: Vec<&str> = s.split('|').collect();
    if parts.len() != 13 {
        return None;
    }
    let n: Vec<f64> = parts.iter().filter_map(|p| p.parse().ok()).collect();
    if n.len() != 13 {
        return None;
    }
    let mut t = GpTrsf::identity();
    t.scale = n[0];
    t.matrix = GpMat::new(n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9]);
    t.loc = GpXyz::new(n[10], n[11], n[12]);
    t.shape = TrsfForm::CompoundTrsf;
    Some(t)
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

pub(super) fn tree_lines(n: &XcafDocNode, depth: usize, out: &mut Vec<String>) {
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
/// and per-node `location` are carried as synthetic string attributes so they
/// survive the binary round-trip; views and annotations are appended to the
/// root entry's attribute list.
pub fn xcaf_to_bincaf(doc: &XcafDocument) -> BinXcaf {
    let mut root = node_to_bin_entry(&doc.root);
    for v in &doc.views {
        root.attributes.push(XcafAttribute {
            kind: "__view".into(),
            value: view_to_string(v),
        });
    }
    for a in &doc.annotations {
        root.attributes.push(XcafAttribute {
            kind: "__annotation".into(),
            value: annotation_to_string(a),
        });
    }
    BinXcaf { root }
}
