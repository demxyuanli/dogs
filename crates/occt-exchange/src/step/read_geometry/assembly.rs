//! Assembly (product structure) reading entry points.

use super::super::prelude::*;
use super::*;

/// Product-name tables derived from the `PRODUCT` graph:
/// `product_names` (PRODUCT id -> name), `def_name` (PRODUCT_DEFINITION id ->
/// product name, the NAUO endpoint key) and `rep_name` (representation id ->
/// product name, the label name `STEPCAFControl_Reader::SetName` writes,
/// `cxx:1725-1760`).
///
/// Attribute positions follow `RWStepBasic_RWProductDefinition`
/// (`ReadStep` `:36-77`): `product_definition(id, description, formation,
/// frame_of_reference)`, so the formation is argument 3 (index 2). The name
/// falls back to the product id when the name itself is empty
/// (`STEPCAFControl_Reader.cxx:1745-1750`).
pub(super) struct ProductNameMaps {
    product_names: HashMap<usize, String>,
    pub(super) def_name: HashMap<usize, String>,
    pub(super) rep_name: HashMap<usize, String>,
}

pub(super) fn product_name_maps(records: &HashMap<usize, Record>) -> ProductNameMaps {
    // PRODUCT id -> name shown on the label: `Prod->Name()` when it has content,
    // else `Prod->Id()` (`STEPCAFControl_Reader.cxx:1741-1753`).
    let mut product_names: HashMap<usize, String> = HashMap::new();
    for (id, rec) in records {
        if rec.type_name == "PRODUCT" {
            let name = parse_str(&rec.args[0]);
            let label = if name.is_empty() { parse_str(&rec.args[1]) } else { name };
            product_names.insert(*id, label);
        }
    }

    // PRODUCT_DEFINITION_FORMATION id -> owning PRODUCT id.
    let mut form_product: HashMap<usize, usize> = HashMap::new();
    for (id, rec) in records {
        if rec.type_name == "PRODUCT_DEFINITION_FORMATION" {
            if let Some(p) = parse_ref(&rec.args[2]) {
                form_product.insert(*id, p);
            }
        }
    }

    // PRODUCT_DEFINITION id -> product name (via its formation).
    let mut def_name: HashMap<usize, String> = HashMap::new();
    for (id, rec) in records {
        if rec.type_name == "PRODUCT_DEFINITION" {
            let name = parse_ref(&rec.args[2])
                .and_then(|f| form_product.get(&f).copied())
                .and_then(|p| product_names.get(&p).cloned())
                .unwrap_or_default();
            def_name.insert(*id, name);
        }
    }

    // PRODUCT_DEFINITION_SHAPE id -> product name. Its `definition` field
    // (`characterized_definition`) is a PRODUCT_DEFINITION in a standard file
    // (`SHAPE_DEFINITION_REPRESENTATION`, as OCCT writes) and was a PRODUCT id
    // in this port's earlier files; accept both.
    let mut pds_name: HashMap<usize, String> = HashMap::new();
    for (id, rec) in records {
        if rec.type_name == "PRODUCT_DEFINITION_SHAPE" {
            if let Some(d) = parse_ref(&rec.args[2]) {
                if let Some(n) = def_name.get(&d).or_else(|| product_names.get(&d)) {
                    pds_name.insert(*id, n.clone());
                }
            }
        }
    }

    // representation id -> product name, via the shape-definition link.
    //  * `SHAPE_DEFINITION_REPRESENTATION(definition, used_representation)` is
    //    the ISO 10303-42 entity OCCT's writer emits.
    //  * `PRODUCT_DEFINITION_SHAPE_REPRESENTATION` is this port's legacy record.
    let mut rep_name: HashMap<usize, String> = HashMap::new();
    for (_, rec) in records {
        let (pds, rep) = match rec.type_name.as_str() {
            "SHAPE_DEFINITION_REPRESENTATION" => {
                (parse_ref(&rec.args[0]), parse_ref(&rec.args[1]))
            }
            "PRODUCT_DEFINITION_SHAPE_REPRESENTATION" => {
                (parse_ref(&rec.args[1]), parse_ref(&rec.args[2]))
            }
            _ => (None, None),
        };
        if let (Some(pds), Some(rep)) = (pds, rep) {
            if let Some(n) = pds_name.get(&pds) {
                rep_name.insert(rep, n.clone());
            }
        }
    }
    // The label name `SetName` writes is the product name, else the product id
    // (`STEPCAFControl_Reader.cxx:1741-1753`), which `product_names` already
    // holds.
    ProductNameMaps { product_names, def_name, rep_name }
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

    // PRODUCT / PRODUCT_DEFINITION / representation -> product-name tables.
    let ProductNameMaps { product_names, def_name, rep_name } = product_name_maps(&records);

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
        let Some(name) = rep_name.get(&id) else { continue };
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
