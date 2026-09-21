use super::prelude::*;
use super::*;


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
pub(super) fn resolve_solid_from_records(
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
pub(super) fn assembly_nauo_records(step_text: &str, edges: &[(String, String)]) -> String {
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
pub(super) fn product_definition_by_name(
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
pub(super) fn collect_edges(n: &XcafDocNode, parent: Option<&str>, out: &mut Vec<(String, String)>) {
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
pub(super) fn write_tree_block(out: &mut String, n: &XcafDocNode, depth: usize) {
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

/// Serialize a placed [`XcafDocument`] assembly to a STEP physical file.
///
/// Unlike [`document_to_step`], each shape is first expanded to world
/// coordinates (its accumulated placement baked into the geometry), so the
/// STEP file carries the instances at their true position. The parent→child
/// occurrences are emitted as `NEXT_ASSEMBLY_USAGE_OCCURRENCE` records and the
/// tree / views / annotations are embedded as a `/* XCAFDOC-BEGIN … */` comment
/// block before `DATA`, keeping the file valid.
pub fn write_step_assembly_with_placements(doc: &XcafDocument) -> String {
    let expanded = xcaf_expand_instances(doc);
    let step = crate::step::write_step_shapes(&expanded);

    let mut edges = Vec::new();
    collect_edges(&doc.root, None, &mut edges);
    let nauo = assembly_nauo_records(&step, &edges);

    let mut block = String::new();
    block.push_str("/* XCAFDOC-BEGIN\n");
    write_tree_block(&mut block, &doc.root, 0);
    for v in &doc.views {
        block.push_str(&format!("VIEW|{}\n", view_to_string(v)));
    }
    for a in &doc.annotations {
        block.push_str(&format!("ANN|{}\n", annotation_to_string(a)));
    }
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
