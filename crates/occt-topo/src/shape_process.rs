//! `ShapeProcess_ShapeContext`'s modification map and the walk that fills it.
//! Source: `ShapeProcess_ShapeContext.hxx:112-144` (`myUntil`, `myMap`),
//! `ShapeProcess_ShapeContext.cxx:266-359` (`RecModif`,
//! `RecordModification`), and `XSAlgo_ShapeProcessor::MergeShapeTransferInfo`
//! (`XSAlgo_ShapeProcessor.cxx:133-197`).
//!
//! `ShapeProcess_OperLibrary`'s `fixshape` operation runs `ShapeFix_Shape` and
//! then keeps the result in the context:
//!
//! ```text
//! sfs->Init(ctx->Result());
//! sfs->Perform(aPS.Next());
//! TopoDS_Shape result = sfs->Shape();
//! if (result != ctx->Result() || ...) {
//!   ctx->RecordModification(sfs->Context(), msg);
//!   ctx->SetResult(result);
//! }
//! ```
//! (`ShapeProcess_OperLibrary.cxx:885-898`).
//!
//! `RecordModification` walks the *pre-heal* root (`myShape`, i.e. `ctx->Result()`
//! before `SetResult`) down to `myUntil`, and for every shape on the way stores
//! the shape's final representative in `myMap`. `MergeTransferInfo` then replays
//! that map onto the transfer binders (`XSAlgo_ShapeProcessor.cxx:200-214`), so a
//! binder that held a face which `ShapeFix_Shape` rebuilt now holds the rebuilt
//! face — the face the *model* is eventually made of.
//!
//! Without this step a healed face keeps two identities: the binder's pre-heal
//! `TopoDS_Face` and the face inside the repaired solid. `XCAFDoc_ShapeTool::
//! SearchUsingMap` (`XCAFDoc_ShapeTool.cxx:134-195`) hashes by `IsSame`, so the
//! pre-heal face finds no label, the styled item that named it loses its color
//! and the face falls back to the node color: ATU01038 exported 70 `usemtl`
//! runs instead of OCCT's 38, 798 triangles of `mat_8`/`mat_6` painted with the
//! solid's `mat_2`.
//!
//! The walk starts at the root and stops at `TopAbs_FACE`: `myUntil` is
//! `ShapeProcess_ShapeContext`'s ctor default (`cxx:29-37`), because
//! `XSAlgo_ShapeProcessor::initializeContext` only lowers it for an explicit
//! `DetalizationLevel` parameter (`XSAlgo_ShapeProcessor.cxx:88-97`) and
//! `FromSTEP` never sets one (`STEPControl_Controller.cxx:199-259`).

use std::collections::HashMap;

use crate::abs::ShapeType;
use crate::shape::TopoShape;
use crate::shape_fix_compose_shell::{MapReShape, ReShape, SharedReShape};
use occt_core::toploc::TopLocLocation;

/// `ShapeProcess_ShapeContext::myMap` (`ShapeProcess_ShapeContext.hxx:144`): an
/// `NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher>`,
/// so both storing and reading go through `IsSame` (the TShape address) —
/// orientation and, as everywhere in this port's `ReShape`, location are not
/// part of the key. The key shape is kept alive beside its replacement for the
/// same reason as `MapReShape::map`: a dangling TShape address would let a
/// later shape be read as the replaced one.
#[derive(Default)]
pub struct ShapeProcessMap {
    map: HashMap<usize, (TopoShape, TopoShape)>,
}

impl ShapeProcessMap {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }

    /// `NCollection_DataMap::Bind` (`ShapeProcess_ShapeContext.cxx:292`,
    /// `:301`).
    pub fn bind(&mut self, old: &TopoShape, new: &TopoShape) {
        self.map.insert(key(old), (old.clone(), new.clone()));
    }

    /// `NCollection_DataMap::Find`; `None` is `!IsBound`
    /// (`ShapeProcess_ShapeContext.cxx:284-286`).
    pub fn find(&self, s: &TopoShape) -> Option<TopoShape> {
        self.map.get(&key(s)).map(|(_, new)| new.clone())
    }

    /// `XSAlgo_ShapeProcessor::MergeShapeTransferInfo`
    /// (`XSAlgo_ShapeProcessor.cxx:133-197`): what the binder bound to
    /// `anOriginalShape` holds after the shape healing, or `None` when the map
    /// leaves the binder alone.
    pub fn transferred(&self, original: &TopoShape) -> Option<TopoShape> {
        // `if (theModifiedShapesMap.IsEmpty()) return;` (`cxx:140-143`).
        if self.map.is_empty() {
            return None;
        }
        // `if (theModifiedShapesMap.IsBound(anOriginalShape))
        // aShapeBinder->SetResult(theModifiedShapesMap.Find(anOriginalShape));`
        // (`cxx:155-158`).
        if let Some(found) = self.find(original) {
            return Some(found);
        }
        // `else if (!anOriginalShape.Location().IsIdentity()) { ... IsBound
        // (aTemporaryShape) ... }` (`cxx:160-168`).
        if !original.location().is_identity() {
            let temporary = original.located(&TopLocLocation::identity());
            return self.find(&temporary);
        }
        // `else { // Some of edges may be modified. ... }` (`cxx:169-190`): the
        // binder holds a composite whose own key was not rebuilt, but whose
        // edges were. A fresh `BRepTools_ReShape` replays exactly those edge
        // substitutions and the resulting shape goes to the binder.
        let mut reshaper = MapReShape::new();
        let mut has_modified_edges = false;
        for edge in crate::topo_tools::edges_of(original) {
            if let Some(modified) = self.find(&edge.0) {
                has_modified_edges = true;
                reshaper.replace(&edge.0, &modified);
            }
        }
        if has_modified_edges {
            return Some(reshaper.apply(original));
        }
        None
    }
}

fn key(s: &TopoShape) -> usize {
    std::sync::Arc::as_ptr(&s.tshape) as usize
}

/// `ShapeProcess_ShapeContext::RecordModification(ShapeBuild_ReShape)`
/// (`ShapeProcess_ShapeContext.cxx:341-357`) followed by `SetResult`
/// (`cxx:113-117`): the `myMap` a `ShapeFix_Shape` pass leaves behind, built
/// from the pass's `ShapeBuild_ReShape`, the root shape the pass started from
/// and the root it produced.
pub fn record_modification(
    root: &TopoShape,
    healed: &TopoShape,
    ctx: &mut SharedReShape,
) -> ShapeProcessMap {
    let mut map = ShapeProcessMap::new();
    rec_modif(root, ctx, &mut map, ShapeType::Face);
    // `SetResult`: `myMap.Bind(myShape, myResult)` (`cxx:116`), the root of the
    // pass bound to the root the pass produced.
    map.bind(root, healed);
    map
}

/// `ShapeProcess_ShapeContext::RecModif`
/// (`ShapeProcess_ShapeContext.cxx:266-340`).
fn rec_modif(
    s: &TopoShape,
    ctx: &mut SharedReShape,
    map: &mut ShapeProcessMap,
    until: ShapeType,
) {
    if s.is_null() {
        return; // `cxx:274-277`.
    }
    // `cxx:279-287`: the stored key is the shape with its location removed
    // ("gka - modification to keep history for shape with location", OCC21617);
    // `r` starts at that key, or at the entry already stored for it.
    let s_val = s.located(&TopLocLocation::identity());
    let mut r = map.find(&s_val).unwrap_or_else(|| s_val.clone());
    if !r.is_null() {
        // `if (repl->Status(r, res, true) && res != r) map.Bind(aS, res);`
        // (`cxx:288-294`).
        if let Some(res) = ctx.status_last(&r) {
            if !crate::topo_tools::is_equal(&res, &r) {
                map.bind(&s_val, &res);
            }
        }
        // `else if (r.ShapeType() < S.ShapeType()) { res = repl->Apply(r,
        // (TopAbs_ShapeEnum)(S.ShapeType() + 1)); if (res != r) map.Bind(aS,
        // res); }` (`cxx:296-303`): the current representative is *lower* in the
        // type order than the shape it stands for (a shell replaced by a face,
        // say), so it is rebuilt from the substitutions below it.
        else if (r.shape_type() as u8) < (s.shape_type() as u8) {
            let res = ctx.apply(&r);
            if !crate::topo_tools::is_equal(&res, &r) {
                map.bind(&s_val, &res);
            }
        }
    }
    // `if (until == TopAbs_SHAPE || S.ShapeType() >= until) return;`
    // (`cxx:332-335`).
    if until == ShapeType::Shape || (s.shape_type() as u8) >= (until as u8) {
        return;
    }
    // `for (TopoDS_Iterator it(S, false); ...) RecModif(it.Value(), ...)`
    // (`cxx:337-340`): direct children, `cumOri` false.
    let children: Vec<TopoShape> = {
        let ts = s.tshape.read().expect("poisoned TShape lock");
        ts.children.clone()
    };
    for child in children {
        rec_modif(&child, ctx, map, until);
    }
}
