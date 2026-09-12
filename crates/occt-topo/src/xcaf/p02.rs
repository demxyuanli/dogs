use super::prelude::*;
use super::*;


pub(super) fn node_to_bin_entry(n: &XcafDocNode) -> BinXcafEntry {
    let mut attributes = attrs_to_strings(&n.attrs);
    if let Some(inst) = &n.instance_name {
        attributes.push(XcafAttribute {
            kind: "instance_name".into(),
            value: inst.clone(),
        });
    }
    if let Some(t) = &n.location {
        attributes.push(XcafAttribute {
            kind: "location".into(),
            value: transform_to_string(t),
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
    let mut views = Vec::new();
    let mut annotations = Vec::new();
    let mut normal = Vec::new();
    for a in &b.root.attributes {
        match a.kind.as_str() {
            "__view" => {
                if let Some(v) = view_from_string(&a.value) {
                    views.push(v);
                }
            }
            "__annotation" => {
                if let Some(x) = annotation_from_string(&a.value) {
                    annotations.push(x);
                }
            }
            _ => normal.push(a.clone()),
        }
    }
    let root = BinXcafEntry {
        shape: b.root.shape.clone(),
        attributes: normal,
        children: b.root.children.clone(),
    };
    XcafDocument {
        root: bin_entry_to_node(&root),
        version: "1.0".to_string(),
        views,
        annotations,
    }
}

pub(super) fn bin_entry_to_node(e: &BinXcafEntry) -> XcafDocNode {
    let mut instance_name = None;
    let mut location = None;
    let mut attrs = Vec::new();
    for a in &e.attributes {
        match a.kind.as_str() {
            "instance_name" => instance_name = Some(a.value.clone()),
            "location" => location = transform_from_string(&a.value),
            _ => attrs.push(a.clone()),
        }
    }
    XcafDocNode {
        shape: e.shape.clone(),
        attrs: strings_to_attrs(&attrs),
        children: e.children.iter().map(bin_entry_to_node).collect(),
        instance_name,
        location,
    }
}

// ---------------------------------------------------------------------------
// XmlXcaf conversion (Phase 7 XML container)
// ---------------------------------------------------------------------------

/// Serialize a document to the XML container format.
///
/// Views and annotations are appended to the root entry's attribute list as
/// synthetic `__view` / `__annotation` attributes; per-node placements are
/// carried as `location` attributes on their own entries.
pub fn xcaf_to_xml(doc: &XcafDocument) -> String {
    let mut root = node_to_xml_entry(&doc.root);
    for v in &doc.views {
        root.attributes.push(crate::xmlcaf::XmlAttribute {
            kind: "__view".into(),
            value: view_to_string(v),
        });
    }
    for a in &doc.annotations {
        root.attributes.push(crate::xmlcaf::XmlAttribute {
            kind: "__annotation".into(),
            value: annotation_to_string(a),
        });
    }
    let xdoc = XmlXcafDoc {
        version: doc.version.clone(),
        root,
    };
    crate::xmlcaf::to_xml(&xdoc).unwrap_or_default()
}

pub(super) fn node_to_xml_entry(n: &XcafDocNode) -> XmlEntry {
    let mut attributes: Vec<crate::xmlcaf::XmlAttribute> = attrs_to_strings(&n.attrs)
        .iter()
        .map(xcaf_attr_to_xml)
        .collect();
    if let Some(t) = &n.location {
        attributes.push(crate::xmlcaf::XmlAttribute {
            kind: "location".into(),
            value: transform_to_string(t),
        });
    }
    XmlEntry {
        name: n.instance_name.clone().unwrap_or_default(),
        shape: n.shape.clone(),
        attributes,
        children: n.children.iter().map(node_to_xml_entry).collect(),
    }
}

pub(super) fn xcaf_attr_to_xml(a: &XcafAttribute) -> crate::xmlcaf::XmlAttribute {
    crate::xmlcaf::XmlAttribute {
        kind: a.kind.clone(),
        value: a.value.clone(),
    }
}

pub(super) fn xml_attr_to_xcaf(a: &crate::xmlcaf::XmlAttribute) -> XcafAttribute {
    XcafAttribute {
        kind: a.kind.clone(),
        value: a.value.clone(),
    }
}

/// Parse a document back from the XML container format.
pub fn xml_to_xcaf(xml: &str) -> Result<XcafDocument, String> {
    let xdoc = crate::xmlcaf::from_xml(xml)?;
    let mut views = Vec::new();
    let mut annotations = Vec::new();
    let mut normal = Vec::new();
    for a in &xdoc.root.attributes {
        match a.kind.as_str() {
            "__view" => {
                if let Some(v) = view_from_string(&a.value) {
                    views.push(v);
                }
            }
            "__annotation" => {
                if let Some(x) = annotation_from_string(&a.value) {
                    annotations.push(x);
                }
            }
            _ => normal.push(a.clone()),
        }
    }
    let root_entry = XmlEntry {
        name: xdoc.root.name.clone(),
        shape: xdoc.root.shape.clone(),
        attributes: normal,
        children: xdoc.root.children.clone(),
    };
    Ok(XcafDocument {
        version: xdoc.version,
        root: xml_entry_to_node(&root_entry),
        views,
        annotations,
    })
}

pub(super) fn xml_entry_to_node(e: &XmlEntry) -> XcafDocNode {
    let mut instance_name = None;
    let mut location = None;
    let mut attrs = Vec::new();
    for a in &e.attributes {
        match a.kind.as_str() {
            "instance_name" => instance_name = Some(a.value.clone()),
            "location" => location = transform_from_string(&a.value),
            _ => attrs.push(xml_attr_to_xcaf(a)),
        }
    }
    XcafDocNode {
        shape: e.shape.clone(),
        attrs: strings_to_attrs(&attrs),
        children: e.children.iter().map(xml_entry_to_node).collect(),
        instance_name: if e.name.is_empty() {
            instance_name
        } else {
            Some(e.name.clone())
        },
        location,
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

pub(super) fn find_node<'a>(n: &'a XcafDocNode, name: &str) -> Option<&'a XcafDocNode> {
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

pub(super) fn collect_shapes(n: &XcafDocNode, out: &mut Vec<(TopoShape, XcafAttrs)>) {
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
// Views / annotations / instances API (Phase 12 XCAF schema)
// ---------------------------------------------------------------------------

/// Register a named camera view on the document.
pub fn xcaf_add_view(doc: &mut XcafDocument, view: XcafView) {
    doc.views.push(view);
}

/// The named views registered on the document, in insertion order.
pub fn xcaf_list_views(doc: &XcafDocument) -> &[XcafView] {
    &doc.views
}

/// Register an annotation (dimension / note / callout) on the document.
pub fn xcaf_add_annotation(doc: &mut XcafDocument, annotation: XcafAnnotation) {
    doc.annotations.push(annotation);
}

/// Remove an annotation by its `id`; returns whether one was removed.
pub fn xcaf_remove_annotation(doc: &mut XcafDocument, id: &str) -> bool {
    let before = doc.annotations.len();
    doc.annotations.retain(|a| a.id != id);
    doc.annotations.len() != before
}

/// Look up an annotation by its `id`.
pub fn xcaf_find_annotation<'a>(doc: &'a XcafDocument, id: &str) -> Option<&'a XcafAnnotation> {
    doc.annotations.iter().find(|a| a.id == id)
}

/// All annotations on the document, in insertion order.
pub fn xcaf_annotations(doc: &XcafDocument) -> &[XcafAnnotation] {
    &doc.annotations
}

/// Add `child` as a placed occurrence of the root assembly. The placement is
/// stored as the child's `location`; the returned name is the occurrence name
/// (the child's `instance_name`, else its `name` attribute, else a generated
/// one).
pub fn xcaf_add_instance(doc: &mut XcafDocument, mut child: XcafDocNode, transform: GpTrsf) -> String {
    let name = child
        .instance_name
        .clone()
        .or_else(|| child.attrs.name.clone())
        .unwrap_or_else(|| format!("Instance{}", doc.node_count() + 1));
    child.instance_name = Some(name.clone());
    child.location = Some(transform);
    doc.root.children.push(child);
    name
}

/// Number of placed occurrences in the tree (nodes carrying a `location`).
pub fn xcaf_instance_count(doc: &XcafDocument) -> usize {
    count_located(&doc.root)
}

pub(super) fn count_located(n: &XcafDocNode) -> usize {
    let self_count = usize::from(n.location.is_some());
    self_count + n.children.iter().map(count_located).sum::<usize>()
}

/// Look up a placed occurrence by its occurrence name (or `name` attribute),
/// returning its referenced sub-tree and placement.
pub fn xcaf_find_instance(doc: &XcafDocument, name: &str) -> Option<XcafInstance> {
    let node = find_node(&doc.root, name)?;
    Some(XcafInstance {
        name: node.name().to_string(),
        node: node.clone(),
        transform: node.location.clone().unwrap_or_else(GpTrsf::identity),
    })
}

/// Expand the assembly tree into world-coordinate shapes.
///
/// Each node that carries a shape is returned once, in depth-first order, as
/// `(name, shape)` where `shape`'s geometry has been baked into world
/// coordinates (the accumulated product of ancestor placements). Nested
/// instances therefore appear at their true world position. Source:
/// `XCAFDoc_ShapeTool::GetShape` + `TopoDS::Transformed`.
pub fn xcaf_expand_instances(doc: &XcafDocument) -> Vec<(String, TopoShape)> {
    let mut out = Vec::new();
    expand_node(&mut out, &doc.root, &GpTrsf::identity());
    out
}

pub(super) fn expand_node(out: &mut Vec<(String, TopoShape)>, n: &XcafDocNode, parent: &GpTrsf) {
    let local = n.location.clone().unwrap_or_else(GpTrsf::identity);
    let world = parent.multiplied(&local);
    if let Some(shape) = &n.shape {
        if let Ok(world_shape) = crate::shape_ops::transformed_copy(shape, &world) {
            out.push((n.name().to_string(), world_shape));
        }
    }
    for c in &n.children {
        expand_node(out, c, &world);
    }
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
pub(super) struct StepRecord {
    pub(super) type_name: String,
    pub(super) args: Vec<String>,
}

impl StepRecord {
    /// Reconstruct the original entity body `TYPE(a,b,...)` from the args.
    pub(super) fn body(&self) -> String {
        format!("{}({})", self.type_name, self.args.join(","))
    }
}

/// Split a comma-separated argument list at the top nesting level, respecting
/// strings (with `''` escapes), parens and brackets.
pub(super) fn split_top(s: &str) -> Vec<String> {
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
pub(super) fn parse_entity_body(body: &str) -> (String, Vec<String>) {
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
pub(super) fn scan_entity_body(data: &str, start: usize) -> Result<(String, usize), String> {
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
pub(super) fn scan_data_records(data: &str) -> Result<Vec<(usize, StepRecord)>, String> {
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
pub(super) fn parse_step_data(content: &str) -> Result<HashMap<usize, StepRecord>, String> {
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

pub(super) fn parse_ref(s: &str) -> Option<usize> {
    s.trim().strip_prefix('#')?.trim().parse().ok()
}

pub(super) fn parse_str(s: &str) -> String {
    let s = s.trim();
    if s.len() >= 2 && s.starts_with('\'') && s.ends_with('\'') {
        s[1..s.len() - 1].replace("''", "'")
    } else {
        s.to_string()
    }
}

/// Parse a `(x,y,z)` tuple into `[x, y, z]`.
pub(super) fn parse_xyz3(s: &str) -> Option<[f64; 3]> {
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
pub(super) fn extract_refs(s: &str) -> Vec<usize> {
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
