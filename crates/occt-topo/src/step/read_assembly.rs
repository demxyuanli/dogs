//! `STEPControl_ActorRead`'s assembly arms: `TransferEntity(NAUO)`
//! (`STEPControl_ActorRead.cxx:859-981`) and
//! `STEPConstruct_Assembly::CheckSRRReversesNAUO`
//! (`STEPConstruct_Assembly.cxx:180-276`).
//!
//! `STEPControl_Reader::NbRootsForTransfer` (`cxx:327-458`) hands back only the
//! `PRODUCT_DEFINITION`s no `NEXT_ASSEMBLY_USAGE_OCCURRENCE` relates, so an
//! assembly file's whole content hangs off the single root product: the root's
//! own shape representation is usually empty (its items are the assembly's axis
//! placements) and every part arrives through a NAUO. Without this arm such a
//! file transfers to no shapes at all.
//!
//! `STEPCAFControl_Reader` keeps one document label per product, so the tree
//! this module walks is emitted as one node per product that owns geometry,
//! each placed by the composed NAUO transformations - the node set
//! `XCAFPrs_DocumentExplorer` gives `RWObj_CafWriter`.
//!
//! OCCT walks the entity graph (`Interface_Graph::Sharings` / `Shareds`, which
//! are indexed); the port's records live in a map, so the reverse edges the
//! assembly walk needs are built once in [`AssemblyIndex`].
use super::prelude::*;
use super::*;

/// One product's document node: the geometry it owns, the product definition it
/// belongs to and the representation that geometry came from.
pub(super) struct ProductNode {
    pub pd: usize,
    pub rep: usize,
    pub shape: TopoShape,
}

/// The reverse edges of the assembly graph, built in one pass:
/// * `nauos_by_relating` - `relating_product_definition` -> NAUOs
///   (`STEPControl_ActorRead.cxx:653-663`),
/// * `pds_by_definition` - `PRODUCT_DEFINITION_SHAPE.definition` -> the PDS
///   entities that characterize it,
/// * `sdrs_by_pds` - PDS -> `(sdr, used_representation)` (`getSDR`,
///   `STEPControl_ActorRead.cxx:555-589`),
/// * `pd_by_rep` - representation -> the `PRODUCT_DEFINITION` that owns it
///   (`STEPConstruct_Assembly.cxx:205-240`),
/// * `cdsrs_by_pds` - PDS -> `CONTEXT_DEPENDENT_SHAPE_REPRESENTATION`s
///   (`STEPControl_ActorRead.cxx:881-889`).
pub(super) struct AssemblyIndex {
    nauos_by_relating: HashMap<usize, Vec<usize>>,
    pds_by_definition: HashMap<usize, Vec<usize>>,
    sdrs_by_pds: HashMap<usize, Vec<(usize, usize)>>,
    pd_by_rep: HashMap<usize, usize>,
    cdsrs_by_pds: HashMap<usize, Vec<usize>>,
}

/// Sorted push, so every lookup answers in entity order the way
/// `Interface_Graph` iterators do.
fn push_sorted(map: &mut HashMap<usize, Vec<usize>>, key: usize, value: usize) {
    let list = map.entry(key).or_default();
    match list.binary_search(&value) {
        Ok(_) => {}
        Err(pos) => list.insert(pos, value),
    }
}

impl AssemblyIndex {
    pub(super) fn new(records: &HashMap<usize, Record>) -> Self {
        let mut index = Self {
            nauos_by_relating: HashMap::new(),
            pds_by_definition: HashMap::new(),
            sdrs_by_pds: HashMap::new(),
            pd_by_rep: HashMap::new(),
            cdsrs_by_pds: HashMap::new(),
        };
        for (id, rec) in records {
            match rec.type_name.as_str() {
                // `next_assembly_usage_occurrence(id, name, description,
                // relating_product_definition, related_product_definition,
                // reference_designator)`.
                "NEXT_ASSEMBLY_USAGE_OCCURRENCE" => {
                    if let Some(relating) = rec.args.get(3).and_then(|a| parse_ref(a)) {
                        push_sorted(&mut index.nauos_by_relating, relating, *id);
                    }
                }
                // `product_definition_shape(name, description,
                // definition)`.
                "PRODUCT_DEFINITION_SHAPE" => {
                    if let Some(def) = rec.args.get(2).and_then(|a| parse_ref(a)) {
                        push_sorted(&mut index.pds_by_definition, def, *id);
                    }
                }
                // `shape_definition_representation(definition,
                // used_representation)`.
                "SHAPE_DEFINITION_REPRESENTATION" => {
                    let (Some(pds), Some(rep)) = (
                        rec.args.first().and_then(|a| parse_ref(a)),
                        rec.args.get(1).and_then(|a| parse_ref(a)),
                    ) else {
                        continue;
                    };
                    index
                        .sdrs_by_pds
                        .entry(pds)
                        .or_default()
                        .push((*id, rep));
                    // SDR -> PDS -> `PRODUCT_DEFINITION`
                    // (`STEPConstruct_Assembly.cxx:205-240`): the first SDR of a
                    // representation wins, as `pd1`/`pd2` there are only filled
                    // while null.
                    if !index.pd_by_rep.contains_key(&rep) {
                        let pd = records
                            .get(&pds)
                            .filter(|r| r.type_name == "PRODUCT_DEFINITION_SHAPE")
                            .and_then(|r| r.args.get(2))
                            .and_then(|a| parse_ref(a));
                        if let Some(pd) = pd {
                            if records
                                .get(&pd)
                                .map(|r| r.type_name == "PRODUCT_DEFINITION")
                                .unwrap_or(false)
                            {
                                index.pd_by_rep.insert(rep, pd);
                            }
                        }
                    }
                }
                // `context_dependent_shape_representation(
                // representation_relation, represented_product_relation)`.
                "CONTEXT_DEPENDENT_SHAPE_REPRESENTATION" => {
                    if let Some(pds) = rec.args.get(1).and_then(|a| parse_ref(a)) {
                        push_sorted(&mut index.cdsrs_by_pds, pds, *id);
                    }
                }
                _ => {}
            }
        }
        for list in index.sdrs_by_pds.values_mut() {
            list.sort_unstable();
        }
        index
    }

    /// The `NEXT_ASSEMBLY_USAGE_OCCURRENCE`s whose `relating_product_definition`
    /// is `pd`.
    fn sub_assemblies(&self, pd: usize) -> &[usize] {
        self.nauos_by_relating
            .get(&pd)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// Whether `pd` composes sub-assemblies, i.e. whether the assembly arm has
    /// to run for it (`STEPControl_ActorRead.cxx:713-758`).
    pub(super) fn has_sub_assemblies(&self, pd: usize) -> bool {
        !self.sub_assemblies(pd).is_empty()
    }

    /// `PRODUCT_DEFINITION_SHAPE`s whose `definition` is `pd`.
    fn pds_of(&self, pd: usize) -> &[usize] {
        self.pds_by_definition
            .get(&pd)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    fn sdrs_of(&self, pds: usize) -> &[(usize, usize)] {
        self.sdrs_by_pds
            .get(&pds)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }
}

/// `relating_product_definition` / `related_product_definition` of a NAUO.
fn nauo_definitions(rec: &Record) -> Option<(usize, usize)> {
    Some((parse_ref(rec.args.get(3)?)?, parse_ref(rec.args.get(4)?)?))
}

/// `rep_1` / `rep_2` and the transformation operator of a representation
/// relationship. The `SRRWT` a `CDSR` points at is a complex instance
/// (`(REPRESENTATION_RELATIONSHIP() REPRESENTATION_RELATIONSHIP_WITH_
/// TRANSFORMATION(op) SHAPE_REPRESENTATION_RELATIONSHIP())`), which
/// `parse_entity_body` keeps in `members`.
fn representation_relationship(rec: &Record) -> Option<(usize, usize, Option<usize>)> {
    let mut rep1 = None;
    let mut rep2 = None;
    let mut op = None;
    if rec.members.is_empty() {
        // Simple `SHAPE_REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION`:
        // `(name, description, rep_1, rep_2, transformation_operator)`.
        rep1 = parse_ref(rec.args.get(2)?);
        rep2 = parse_ref(rec.args.get(3)?);
        op = rec.args.get(4).and_then(|a| parse_ref(a));
    } else {
        for (name, args) in &rec.members {
            match name.as_str() {
                // `representation_relationship(name, description, rep_1, rep_2)`.
                "REPRESENTATION_RELATIONSHIP" => {
                    rep1 = args.get(2).and_then(|a| parse_ref(a));
                    rep2 = args.get(3).and_then(|a| parse_ref(a));
                }
                // `representation_relationship_with_transformation(
                // transformation_operator)`.
                "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION" => {
                    op = args.first().and_then(|a| parse_ref(a));
                }
                _ => {}
            }
        }
    }
    Some((rep1?, rep2?, op))
}

/// `STEPConstruct_Assembly::CheckSRRReversesNAUO`
/// (`STEPConstruct_Assembly.cxx:180-276`): `true` when the CDSR's relationship
/// runs parent -> child instead of child -> parent, i.e. the transformation has
/// to be inverted before it places the child.
fn srr_reverses_nauo(
    index: &AssemblyIndex,
    records: &HashMap<usize, Record>,
    nauo: usize,
    rep1: usize,
    rep2: usize,
) -> bool {
    let Some(nauo_rec) = records.get(&nauo) else { return false };
    let Some((relating, related)) = nauo_definitions(nauo_rec) else { return false };
    // `SDR->Definition().PropertyDefinition()` of the two representations.
    let pd1 = index.pd_by_rep.get(&rep1).copied();
    let pd2 = index.pd_by_rep.get(&rep2).copied();
    if pd1 == Some(related) && pd2 == Some(relating) {
        return false;
    }
    if pd2 == Some(related) && pd1 == Some(relating) {
        return true;
    }
    // `cxx:270-273`: a relationship whose ends do not match the NAUO at all is
    // still treated as reversed when one end lines up.
    pd2 == Some(related) || pd1 == Some(relating)
}

/// What `TransferEntity(NAUO)` (`STEPControl_ActorRead.cxx:859-981`) resolves
/// for one component: the child product definition, the representation its
/// geometry comes from, and the placement (already inverted when the CDSR runs
/// backwards).
fn nauo_target(
    resolver: &Resolver,
    index: &AssemblyIndex,
    records: &HashMap<usize, Record>,
    nauo: usize,
) -> Option<(usize, usize, Option<occt_core::gp::GpTrsf>)> {
    // `graph.Sharings(NAUO)` -> the `PRODUCT_DEFINITION_SHAPE` carrying the
    // `CONTEXT_DEPENDENT_SHAPE_REPRESENTATION` (`cxx:871-886`).
    for &pds in index.pds_of(nauo) {
        let Some(cdsrs) = index.cdsrs_by_pds.get(&pds) else { continue };
        for &cdsr in cdsrs {
            let Some(cdsr_rec) = records.get(&cdsr) else { continue };
            let Some(rr_id) = cdsr_rec.args.first().and_then(|a| parse_ref(a)) else { continue };
            let Some(rr_rec) = records.get(&rr_id) else { continue };
            let Some((rep1, rep2, op)) = representation_relationship(rr_rec) else { continue };
            let reversed = srr_reverses_nauo(index, records, nauo, rep1, rep2);
            // `cxx:896`: `rep = SRRReversed ? RR->Rep2() : RR->Rep1()`.
            let rep = if reversed { rep2 } else { rep1 };
            // `iatrsf = ComputeSRRWT(RR, TP, Trsf, ...)` (`cxx:901`); the caller
            // inverts it when the relationship is reversed (`cxx:953-959`).
            let mut trsf = op.and_then(|op| compute_srrwt_from_op(resolver, records, op));
            if reversed {
                trsf = trsf.and_then(|t| t.inverted().ok());
            }
            // "find real ProductDefinition used rep" (`cxx:903-929`).
            let Some(&pd) = index.pd_by_rep.get(&rep) else { continue };
            return Some((pd, rep, trsf));
        }
    }
    None
}

/// Applies a placement to a transferred child, baking it into a copy as
/// `STEPControl_ActorRead::ApplyTransformation` does (`cxx:950-961`). The
/// geometry registry is location-blind, so the placement cannot ride on
/// `TopoShape::location` (same mechanism as `MAPPED_ITEM`).
fn place(shape: &TopoShape, trsf: Option<&occt_core::gp::GpTrsf>) -> TopoShape {
    match trsf {
        Some(t) if t.form() != occt_core::gp::TrsfForm::Identity => {
            crate::shape_ops::transformed_copy(shape, t).unwrap_or_else(|_| shape.clone())
        }
        _ => shape.clone(),
    }
}

/// Emits one document node per product that owns geometry, walking the tree
/// `TransferEntity(PD)` (`STEPControl_ActorRead.cxx:617-855`) composes.
/// `outer` is the composed placement of the enclosing sub-assemblies.
pub(super) fn collect_product_nodes(
    resolver: &Resolver,
    index: &AssemblyIndex,
    records: &HashMap<usize, Record>,
    composer: &SsrComposer,
    pd: usize,
    outer: Option<&occt_core::gp::GpTrsf>,
    out: &mut Vec<ProductNode>,
    depth: u32,
) {
    // Guards a cyclic NAUO graph (`A contains B contains A`), which OCCT's
    // `TP->IsBound` bookkeeping also stops.
    if depth > 64 {
        return;
    }
    for &pds in index.pds_of(pd) {
        for &(_, rep) in index.sdrs_of(pds) {
            if let Some(shape) = composer.rep_shape(rep) {
                out.push(ProductNode { pd, rep, shape: place(&shape, outer) });
            }
        }
    }
    for &nauo in index.sub_assemblies(pd) {
        if let Some((child_pd, _, trsf)) = nauo_target(resolver, index, records, nauo) {
            // The child's own placement applies first, the enclosing ones
            // after it (`GpTrsf::multiplied` is `self ∘ other`).
            let combined = match (outer, &trsf) {
                (None, None) => None,
                (Some(t), None) => Some(t.clone()),
                (None, Some(u)) => Some(u.clone()),
                (Some(t), Some(u)) => Some(t.multiplied(u)),
            };
            collect_product_nodes(
                resolver,
                index,
                records,
                composer,
                child_pd,
                combined.as_ref(),
                out,
                depth + 1,
            );
        }
    }
}
