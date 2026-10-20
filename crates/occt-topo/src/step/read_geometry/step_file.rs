//! `read_step_impl` and the public reader entry points.

use super::super::prelude::*;
use super::*;

pub(in crate::step) fn read_step_impl(content: &str) -> Result<(BRepModel, Vec<String>), String> {
    let data = data_section(content)?;
    let records = parse_records(data)?;
    let resolver = Resolver::new(&records);

    let mut rep_ids: Vec<usize> = records
        .iter()
        .filter(|(_, r)| {
            r.type_name == "ADVANCED_BREP_SHAPE_REPRESENTATION" || r.type_name == "SHAPE_REPRESENTATION"
        })
        .map(|(id, _)| *id)
        .collect();
    rep_ids.sort_unstable();

    let mut model = BRepModel::new();
    // `STEPCAFControl_Reader::prepareUnits` (`cxx:508-522`):
    // `SetLengthUnit(doc, aScaleFactorMM, Millimeter)` stores
    // `aScaleFactorMM * 1e-3` metres per model unit, and `aScaleFactorMM` is
    // `UnitsMethods::GetCasCadeLengthUnit()` — the representation's length
    // factor, which `context_length_factor` computes in millimetres
    // (`transfer.rs`, `cxx:272-302`).
    model.length_unit = Some(super::transfer::context_length_factor(&records) / 1000.0);
    // `readNameAttribute(RefLabel)` (`RWObj_CafWriter.cxx:197`): the label name
    // `STEPCAFControl_Reader::SetName` wrote, i.e. the PRODUCT name. The
    // representation's own name is the fallback for a file without a product
    // layer.
    let names = product_name_maps(&records);
    // The port's `AddShape` label registry: the entity that produced each model
    // shape, so `ReadColors` can resolve a styled item's target back to it
    // (`STEPCAFControl_Reader.cxx:750` / `:831-861`).
    let mut rep_index: HashMap<usize, usize> = HashMap::new();

    // `STEPCAFControl_Reader::Transfer` turns the transferred roots into
    // document shapes: `num = reader.NbRootsForTransfer()` and one
    // `AddShape(reader.Shape(i))` per root (`cxx:540-557`, `:756-758`). A root is
    // a PRODUCT_DEFINITION (`STEPControl_Reader::NbRootsForTransfer`,
    // `cxx:327-458`) and `reader.Shape(i)` is that product's
    // `TransferEntity(PD)` result (`STEPControl_ActorRead.cxx:617-854`). Iterating
    // representations instead would make a *related* representation
    // (`SHAPE_REPRESENTATION_RELATIONSHIP`, e.g. the `ADVANCED_BREP_...` behind a
    // product's named `SHAPE_REPRESENTATION`) a second node of its own: it would
    // claim the shared faces of the product's node and carry neither its name nor
    // its styles.
    let roots = root_product_definitions(&records);
    if crate::debug_flags::zzstep_dbg() {
        eprintln!(
            "DBG records={} roots={} rep_ids={}",
            records.len(),
            roots.len(),
            rep_ids.len()
        );
    }
    let composer = SsrComposer::new(&resolver, &records);
    // One reverse-relation pass over the record table; the assembly walk below
    // queries it per node instead of rescanning the file for every NAUO
    // (`STEPConstruct_Assembly.cxx` reads the same relations through the
    // transfer processor's bindings).
    let assembly_index = super::read_assembly::AssemblyIndex::new(&records);
    if !roots.is_empty() {
        for pd in roots {
            // `TransferEntity(PD)`'s assembly arm (`cxx:713-758`): the product's
            // shape is the composition of its `NEXT_ASSEMBLY_USAGE_OCCURRENCE`
            // components, and `STEPCAFControl_Reader` gives each component that
            // owns geometry its own document label (`cxx:540-557`). Emitting the
            // per-product nodes keeps the group names an assembly file carries.
            if assembly_index.has_sub_assemblies(pd) {
                let mut nodes = Vec::new();
                super::read_assembly::collect_product_nodes(
                    &resolver,
                    &assembly_index,
                    &records,
                    &composer,
                    pd,
                    None,
                    &mut nodes,
                    0,
                );
                for node in nodes {
                    let name = names
                        .rep_name
                        .get(&node.rep)
                        .filter(|n| !n.is_empty())
                        .cloned()
                        .or_else(|| {
                            names.def_name.get(&node.pd).filter(|n| !n.is_empty()).cloned()
                        })
                        .unwrap_or_else(|| composer.name(node.rep));
                    if crate::debug_flags::zzstep_dbg() {
                        let n = crate::topo_tools::faces_of(&node.shape).len();
                        let b = crate::brep_bnd_lib::shape_bnd_box(&node.shape).get();
                        eprintln!(
                            "DBG pd={} rep={} faces={n} name={name} bbox={b:?}",
                            node.pd, node.rep
                        );
                    }
                    let index = model.add(&name, node.shape);
                    rep_index.insert(node.rep, index);
                    rep_index.insert(node.pd, index);
                }
                continue;
            }
            // The product's own representation: `GetLabelFromPD`'s
            // PDS -> SDR -> `used_representation` chain
            // (`STEPCAFControl_Reader.cxx:1552-1770`).
            let Some(rep) = product_representation(&records, pd) else {
                // A product with neither a representation nor sub-assemblies
                // carries no shape at all (`STEPControl_ActorRead.cxx:713-716`).
                continue;
            };
            // `TransferEntity(PD)` with one SDR and no NAUO binds the
            // representation's own shape as it is (`cxx:842-848`).
            let Some(shape) = composer.rep_shape(rep) else { continue };
            if crate::debug_flags::zzstep_dbg() {
                let n = crate::topo_tools::faces_of(&shape).len();
                eprintln!("DBG pd={pd} rep={rep} faces={n}");
            }
            let name = names
                .rep_name
                .get(&rep)
                .filter(|n| !n.is_empty())
                .cloned()
                .unwrap_or_else(|| composer.name(rep));
            let index = model.add(&name, shape);
            rep_index.insert(rep, index);
            rep_index.insert(pd, index);
        }
    } else if !rep_ids.is_empty() {
        // Fallback for a file without a product layer: keep one node per
        // representation.
        for id in &rep_ids {
            if let Some(shape) = composer.rep_shape(*id) {
                let name = names
                    .rep_name
                    .get(id)
                    .filter(|n| !n.is_empty())
                    .cloned()
                    .unwrap_or_else(|| composer.name(*id));
                rep_index.insert(*id, model.add(&name, shape));
            }
        }
    } else {
        // Fallback for bare files without a representation layer: treat any
        // MANIFOLD_SOLID_BREP as a root shape.
        let mut roots: Vec<usize> = records
            .iter()
            .filter(|(_, r)| r.type_name == "MANIFOLD_SOLID_BREP")
            .map(|(id, _)| *id)
            .collect();
        roots.sort_unstable();
        for id in roots {
            match resolver.resolve_shape(id) {
                Ok(s) => {
                    rep_index.insert(id, model.add("", s));
                }
                Err(e) => resolver.warn(format!("shape #{id}: {e}")),
            }
        }
    }

    // `STEPCAFControl_Reader::Transfer`'s color pass (`cxx:754-760`):
    // `ReadColors(WS, Doc, ...)` after the shapes have been transferred.
    super::read_colors::read_colors(&mut model, &resolver, &rep_index);

    let warnings = resolver.warnings.into_inner();
    Ok((model, warnings))
}

/// `STEPControl_Reader::NbRootsForTransfer` (`cxx:327-458`): every
/// `PRODUCT_DEFINITION` that no `NEXT_ASSEMBLY_USAGE_OCCURRENCE` uses as its
/// `related_product_definition`, in entity order.
///
/// **UNPORTED**: the `ReadAllShapes` special case for unshared
/// `MANIFOLD_SOLID_BREP` / `SHELL_BASED_SURFACE_MODEL` (`cxx:341-363`), gated by
/// `InternalParameters.ReadAllShapes == 1`, and the
/// `PRODUCT_DEFINITION_WITH_ASSOCIATED_DOCUMENTS` /
/// `DOCUMENT_PRODUCT_EQUIVALENCE` exclusion (`cxx:365-404`).
fn root_product_definitions(records: &HashMap<usize, Record>) -> Vec<usize> {
    let mut related: HashSet<usize> = HashSet::new();
    for (_, rec) in records {
        if rec.type_name == "NEXT_ASSEMBLY_USAGE_OCCURRENCE" {
            if let Some(red) = rec.args.get(4).and_then(|a| parse_ref(a)) {
                related.insert(red);
            }
        }
    }
    let mut roots: Vec<usize> = records
        .iter()
        .filter(|(_, r)| r.type_name == "PRODUCT_DEFINITION")
        .map(|(id, _)| *id)
        .filter(|id| !related.contains(id))
        .collect();
    roots.sort_unstable();
    roots
}

/// The representation a product definition owns:
/// `PRODUCT_DEFINITION_SHAPE.definition` (attribute 3) -> `SHAPE_DEFINITION_
/// REPRESENTATION.used_representation` (attribute 2), the chain
/// `STEPCAFControl_Reader::GetLabelFromPD` walks (`cxx:1552-1770`).
fn product_representation(records: &HashMap<usize, Record>, pd: usize) -> Option<usize> {
    let mut pds: Option<usize> = None;
    for (id, rec) in records {
        if rec.type_name == "PRODUCT_DEFINITION_SHAPE"
            && rec.args.get(2).and_then(|a| parse_ref(a)) == Some(pd)
        {
            pds = Some(*id);
        }
    }
    let pds = pds?;
    let mut rep: Option<usize> = None;
    for (_, rec) in records {
        if rec.type_name == "SHAPE_DEFINITION_REPRESENTATION"
            && rec.args.first().and_then(|a| parse_ref(a)) == Some(pds)
        {
            // `NCollection_DataMap::Bind` keeps the last binding.
            rep = rec.args.get(1).and_then(|a| parse_ref(a));
        }
    }
    rep
}

/// Parse a STEP physical file into a `BRepModel`.
pub fn read_step(content: &str) -> Result<BRepModel, String> {
    Ok(read_step_impl(content)?.0)
}
/// Parse a STEP physical file, also returning collected warnings for skipped
/// or unsupported entities.
pub fn read_step_with_warnings(content: &str) -> Result<(BRepModel, Vec<String>), String> {
    read_step_impl(content)
}

/// Read a STEP physical file from disk.
pub fn read_step_file(path: &str) -> Result<BRepModel, String> {
    let content = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    read_step(&content)
}
