//! `BOPTools_Set` — a hashed set of sub-shapes used as a same-domain key.
//!
//! Source: `BOPTools_Set.hxx/.cxx`. `BOPAlgo_Builder::BuildSplitSolids`
//! (`_3.cxx:456-459`, `:593-613`) builds a set of FACE sub-shapes for each
//! resulting solid. Two solids whose face sets compare equal are the same
//! domain: the later solid is bound in `myShapesSD` to the earlier
//! representative, and the image list stores that representative.
//!
//! `BOPTools_Set::Add` (`BOPTools_Set.cxx`) explores `theType` with the
//! equivalent of `TopExp_Explorer`: degenerated edges are skipped; an
//! INTERNAL sub-shape is stored twice (FORWARD and REVERSED, same TShape).
//! `IsEqual` requires equal `NbShapes` and `Contains` (IsSame) of every
//! member. This port stores the explorer list of [`shape_key`] values plus
//! the representative (`Shape()`).

use std::collections::HashSet;
use std::hash::{Hash, Hasher};

use crate::abs::{Orientation, ShapeType};
use crate::bop_occt_util::{iter_children, shape_key};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, TopoShape};

/// A `BOPTools_Set`: explorer list of sub-shape TShape keys plus the
/// representative (`myShape`).
#[derive(Debug, Clone)]
pub struct BopToolsSet {
    /// Representative shape (`BOPTools_Set::Shape`).
    shape: TopoShape,
    /// Explorer-order TShape keys (`myShapes`). Degenerated edges omitted;
    /// INTERNAL sub-shapes appear twice.
    keys: Vec<usize>,
}

impl PartialEq for BopToolsSet {
    fn eq(&self, other: &Self) -> bool {
        if self.keys.len() != other.keys.len() {
            return false;
        }
        let mine: HashSet<usize> = self.keys.iter().copied().collect();
        other.keys.iter().all(|k| mine.contains(k))
    }
}

impl Eq for BopToolsSet {}

impl Hash for BopToolsSet {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.keys.len().hash(state);
        let mut uniq: Vec<usize> = {
            let s: HashSet<usize> = self.keys.iter().copied().collect();
            s.into_iter().collect()
        };
        uniq.sort_unstable();
        uniq.hash(state);
    }
}

impl BopToolsSet {
    /// Empty set.
    pub fn new() -> Self {
        Self {
            shape: TopoShape::new(ShapeType::Compound),
            keys: Vec::new(),
        }
    }

    /// `BOPTools_Set::Add(theS, theType)` — explore `the_s` for `the_type`
    /// and remember `the_s` as the representative.
    pub fn add(&mut self, the_s: &TopoShape, the_type: ShapeType) {
        self.shape = the_s.clone();
        let mut keys: Vec<usize> = Vec::new();
        collect_set_items(the_s, the_type, &mut keys);
        self.keys = keys;
    }

    /// Construct and fill in one call.
    pub fn from_shape(the_s: &TopoShape, the_type: ShapeType) -> Self {
        let mut s = Self::new();
        s.add(the_s, the_type);
        s
    }

    /// The representative shape (`Shape()`).
    pub fn shape(&self) -> &TopoShape {
        &self.shape
    }

    /// Explorer-order TShape keys (`myShapes`).
    pub fn keys(&self) -> &[usize] {
        &self.keys
    }

    /// `NbShapes`.
    pub fn extent(&self) -> usize {
        self.keys.len()
    }

    /// True when the set contains the TShape of `s` (`IsSame`).
    pub fn contains_shape(&self, s: &TopoShape) -> bool {
        let k = shape_key(s);
        self.keys.contains(&k)
    }
}

/// `BOPTools_Set::Add` explorer body (`BOPTools_Set.cxx`).
fn collect_set_items(s: &TopoShape, ty: ShapeType, out: &mut Vec<usize>) {
    if s.shape_type() == ty {
        let skip_degen = ty == ShapeType::Edge && BRepTool::is_degenerated(&Edge(s.clone()));
        if !skip_degen {
            let k = shape_key(s);
            if s.orientation() == Orientation::Internal {
                out.push(k);
                out.push(k);
            } else {
                out.push(k);
            }
        }
    }
    for c in iter_children(s) {
        collect_set_items(&c, ty, out);
    }
}

impl Default for BopToolsSet {
    fn default() -> Self {
        Self::new()
    }
}

/// A set-of-sets with interned representatives (`NCollection_Map<BOPTools_Set>`
/// plus `Added` returning the stored copy).
///
/// `contains` / `added` mirror the BuildSplitSolids same-domain check:
/// ```text
/// bFlagSD = aMST.Contains(aST);
/// const BOPTools_Set& aSTx = aMST.Added(aST);
/// const TopoDS_Shape& aSx  = aSTx.Shape();
/// ```
#[derive(Debug, Clone, Default)]
pub struct BopToolsSetMap {
    sets: Vec<BopToolsSet>,
}

impl BopToolsSetMap {
    pub fn new() -> Self {
        Self { sets: Vec::new() }
    }

    /// `Contains`.
    pub fn contains(&self, set: &BopToolsSet) -> bool {
        self.sets.iter().any(|s| s == set)
    }

    /// `Add` — insert when new, return true if inserted.
    pub fn add(&mut self, set: BopToolsSet) -> bool {
        if self.contains(&set) {
            false
        } else {
            self.sets.push(set);
            true
        }
    }

    /// `Added` — insert when new and return a reference to the stored set
    /// (the interned representative).
    pub fn added(&mut self, set: BopToolsSet) -> &BopToolsSet {
        if let Some(i) = self.sets.iter().position(|s| s == &set) {
            return &self.sets[i];
        }
        self.sets.push(set);
        self.sets.last().unwrap()
    }

    /// Number of interned sets.
    pub fn extent(&self) -> usize {
        self.sets.len()
    }

    /// All interned sets.
    pub fn sets(&self) -> &[BopToolsSet] {
        &self.sets
    }
}

/// Face-set of a solid, the key BuildSplitSolids uses for same-domain solids.
pub fn solid_face_set(solid: &TopoShape) -> BopToolsSet {
    BopToolsSet::from_shape(solid, ShapeType::Face)
}

/// Resolve `solid` against `map`: if a same-domain representative already
/// exists, return `(representative, true)`; otherwise intern `solid` and
/// return `(solid, false)`.
pub fn intern_solid_sd(map: &mut BopToolsSetMap, solid: &TopoShape) -> (TopoShape, bool) {
    let set = solid_face_set(solid);
    let flag_sd = map.contains(&set);
    let stored = map.added(set).shape().clone();
    (stored, flag_sd)
}
