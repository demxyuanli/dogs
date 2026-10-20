//! `SsrComposer` — `CONTEXT_DEPENDENT_SHAPE_REPRESENTATION` transfer.

// ---------------------------------------------------------------------------
// SHAPE_REPRESENTATION_RELATIONSHIP composition
// (`STEPControl_ActorRead::TransferRelatedSRR` and friends)
// ---------------------------------------------------------------------------

use super::super::prelude::*;
use super::*;

/// Composes shape representations the way `STEPControl_ActorRead` does.
///
/// * `TransferEntity(ShapeRepresentation, …)` (`STEPControl_ActorRead.cxx:2048-2092`):
///   the representation's own shape, then — when `read.step.shape.relationship`
///   is set, which defaults to true (`DESTEP_Parameters.hxx:170`) — every shape
///   `TransferRelatedSRR` composes for the relations sharing it. A single child
///   binds that shape, more than one binds a compound; a non-null related result
///   replaces the single result (`cxx:2077-2089`).
/// * `TransferRelatedSRR` (`cxx:2679-2729`): `theCund` grows with each related
///   shape and the **last** one is returned (`cxx:2724-2728`).
/// * `TransferEntity(SRR, …)` (`cxx:1337-1433`): for `nbrep` 1/2 take `Rep1` or
///   `Rep2` (0 takes both), compose them, and apply the `ComputeSRRWT`
///   transformation to the single result or to the compound (`cxx:1408-1417`).
///
/// **UNPORTED**: the `MECHANICAL_DESIGN_AND_DRAUGHTING_RELATIONSHIP` and
/// constructive-geometry arms of `TransferRelatedSRR` (`cxx:2707-2721`).
/// `TransferEntity(NAUO, …)` (`cxx:859-981`), the assembly-level path that
/// resolves `CONTEXT_DEPENDENT_SHAPE_REPRESENTATION` / `SRRReversed`, lives in
/// [`crate::step::read_assembly`]; the relation entities themselves are read
/// directly here, as `TransferRelatedSRR` does.
pub(in crate::step) struct SsrComposer<'a, 'r> {
    resolver: &'a Resolver<'r>,
    records: &'a HashMap<usize, Record>,
    names: RefCell<HashMap<usize, String>>,
    cache: RefCell<HashMap<usize, Option<TopoShape>>>,
    in_progress: RefCell<HashSet<usize>>,
}
impl<'a, 'r> SsrComposer<'a, 'r> {
    pub(in crate::step) fn new(resolver: &'a Resolver<'r>, records: &'a HashMap<usize, Record>) -> Self {
        Self {
            resolver,
            records,
            names: RefCell::new(HashMap::new()),
            cache: RefCell::new(HashMap::new()),
            in_progress: RefCell::new(HashSet::new()),
        }
    }

    /// The representation name recorded by `resolve_representation`.
    pub(in crate::step) fn name(&self, rep_id: usize) -> String {
        if let Some(n) = self.names.borrow().get(&rep_id) {
            return n.clone();
        }
        match self.resolver.resolve_representation(rep_id) {
            Ok((name, _)) => {
                self.names.borrow_mut().insert(rep_id, name.clone());
                name
            }
            Err(_) => format!("shape_{rep_id}"),
        }
    }

    /// `TransferEntity(ShapeRepresentation, …)` (`cxx:2048-2092`).
    pub(in crate::step) fn rep_shape(&self, rep_id: usize) -> Option<TopoShape> {
        if let Some(hit) = self.cache.borrow().get(&rep_id) {
            return hit.clone();
        }
        // Rust-only guard: a cyclic SRR graph would make OCCT's transfer
        // recurse forever (`TP->Bind` happens after the transfer, `cxx:2091`);
        // fall back to the representation's own shape instead of recursing.
        if !self.in_progress.borrow_mut().insert(rep_id) {
            return self.own_shape(rep_id);
        }
        let own = self.own_shape(rep_id);
        let mut cund: Vec<TopoShape> = own.into_iter().collect();
        let new_result = self.transfer_related_srr(rep_id, &mut cund);
        let result = if cund.is_empty() {
            None
        } else if cund.len() == 1 {
            Some(new_result.unwrap_or_else(|| cund[0].clone()))
        } else {
            Some(TopoBuilder::new().make_compound_of(&cund).0)
        };
        self.in_progress.borrow_mut().remove(&rep_id);
        self.cache.borrow_mut().insert(rep_id, result.clone());
        result
    }

    /// The representation's own items (`TransferEntity`'s `aCund` seed,
    /// `cxx:2048-2056`).
    fn own_shape(&self, rep_id: usize) -> Option<TopoShape> {
        match self.resolver.resolve_representation(rep_id) {
            Ok((name, shapes)) => {
                self.names.borrow_mut().insert(rep_id, name);
                match shapes.len() {
                    0 => None,
                    1 => Some(shapes.into_iter().next().unwrap()),
                    _ => Some(TopoBuilder::new().make_compound_of(&shapes).0),
                }
            }
            Err(e) => {
                self.resolver.warn(format!("representation #{rep_id}: {e}"));
                None
            }
        }
    }

    /// `TransferRelatedSRR` (`cxx:2679-2729`).
    fn transfer_related_srr(&self, rep_id: usize, cund: &mut Vec<TopoShape>) -> Option<TopoShape> {
        // `aGraph.Sharings(theRep)` walks the graph's entity order; the parsed
        // records live in a map, so file order is restored by entity number.
        let mut srr_ids: Vec<usize> = self
            .records
            .iter()
            .filter(|(_, r)| {
                r.type_name == "SHAPE_REPRESENTATION_RELATIONSHIP"
                    || r.type_name == "SHAPE_REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION"
            })
            .filter(|(_, r)| {
                let rep1 = r.args.get(2).and_then(|a| parse_ref(a));
                let rep2 = r.args.get(3).and_then(|a| parse_ref(a));
                rep1 == Some(rep_id) || rep2 == Some(rep_id)
            })
            .map(|(id, _)| *id)
            .collect();
        srr_ids.sort_unstable();

        let mut last: Option<TopoShape> = None;
        for srr_id in srr_ids {
            let Some(rec) = self.records.get(&srr_id) else { continue };
            let rep1 = rec.args.get(2).and_then(|a| parse_ref(a));
            let nb_rep = if rep1 == Some(rep_id) { 2 } else { 1 };
            if let Some(shape) = self.transfer_srr(srr_id, nb_rep) {
                cund.push(shape.clone());
                last = Some(shape);
            }
        }
        last
    }

    /// `TransferEntity(ShapeRepresentationRelationship, …)` (`cxx:1337-1433`).
    fn transfer_srr(&self, srr_id: usize, nb_rep: usize) -> Option<TopoShape> {
        let rec = self.records.get(&srr_id)?;
        let rep1 = rec.args.get(2).and_then(|a| parse_ref(a))?;
        let rep2 = rec.args.get(3).and_then(|a| parse_ref(a))?;
        let iatrsf = compute_srrwt(self.resolver, self.records, srr_id);

        let mut one_result: Option<TopoShape> = None;
        let mut children: Vec<TopoShape> = Vec::new();
        for i in 1..=2usize {
            if nb_rep != 0 && nb_rep != i {
                continue;
            }
            let rep = if i == 1 { rep1 } else { rep2 };
            if let Some(shape) = self.rep_shape(rep) {
                one_result = Some(shape.clone());
                children.push(shape);
            }
        }
        if children.is_empty() {
            return None;
        }
        let mut result = if children.len() == 1 {
            one_result.unwrap()
        } else {
            TopoBuilder::new().make_compound_of(&children).0
        };
        if let Some(t) = iatrsf {
            // `ApplyTransformation(OneResult|Cund, Trsf)` (`cxx:1408-1417`); the
            // `SRRReversed` inverse (`cxx:1316-1323`) belongs to the NAUO path,
            // which is UNPORTED.
            result = crate::transform::transformed(&result, &t);
        }
        Some(result)
    }
}
