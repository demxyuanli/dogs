use super::prelude::*;
use super::*;

/// Stable key of a shape inside the history: the address of its shared TShape.
///
/// Reuses [`ShapeId`] so callers can correlate history entries with
/// `ShapeId::of(shape)`. The key is stable for the lifetime of the shape's
/// `Arc`; cloned views of the same shape share the key.

pub type ShapeKey = ShapeId;

/// Returns true if the shape participates in the history.
///
/// Mirrors `BRepTools_History::IsSupportedType` — only vertices, edges, faces
/// and solids are tracked.
pub fn is_supported_type(s: &TopoShape) -> bool {
    matches!(
        s.shape_type(),
        ShapeType::Vertex | ShapeType::Edge | ShapeType::Face | ShapeType::Solid
    )
}

/// Shape-identity test matching the `ShapeKey` domain (same TShape data).
pub(super) fn same_shape(a: &TopoShape, b: &TopoShape) -> bool {
    a.same_tshape(b)
}

/// Appends the source shapes to `dst`, skipping entries already present.
///
/// Mirrors the `add(list, collection)` helper of `BRepTools_History.cxx`.
pub(super) fn append_unique(dst: &mut Vec<TopoShape>, src: impl IntoIterator<Item = TopoShape>) {
    for s in src {
        if !dst.iter().any(|d| same_shape(d, &s)) {
            dst.push(s);
        }
    }
}

/// History of modifications/generations/removals plus the images naming table.
///
/// Public surface is split into per-shape queries (`image`, `modified`,
/// `generated`, `is_deleted`, `has_*`) and whole-table accessors (`images`,
/// `modified_map`, `generated_map`, `removed`).
#[derive(Debug, Clone, Default)]
pub struct BopHistory {
    /// Naming side table: each old shape → the pieces it was split/expanded
    /// into during intersection (`BOPAlgo_Builder::myImages`).
    pub(super) images: HashMap<ShapeKey, Vec<TopoShape>>,
    /// Each old shape → shapes modified from it (`BRepTools_History`).
    pub(super) modified: HashMap<ShapeKey, Vec<TopoShape>>,
    /// Each old shape → shapes generated from it (`BRepTools_History`).
    pub(super) generated: HashMap<ShapeKey, Vec<TopoShape>>,
    /// Shapes completely removed from the result (`BRepTools_History`).
    pub(super) removed: Vec<TopoShape>,
    /// Key mirror of `removed`. `merge` records removals for intermediate
    /// shapes it no longer holds a handle to, which cannot go into the
    /// `Vec<TopoShape>`; the set keeps `is_deleted` correct in every case.
    pub(super) removed_keys: HashSet<ShapeKey>,
}

impl BopHistory {
    /// Empty constructor.
    pub fn new() -> Self {
        Self::default()
    }

    /// Clears every relation.
    pub fn clear(&mut self) {
        self.images.clear();
        self.modified.clear();
        self.generated.clear();
        self.removed.clear();
        self.removed_keys.clear();
    }

    /// True when no relation is recorded at all.
    pub fn is_empty(&self) -> bool {
        self.images.is_empty()
            && self.modified.is_empty()
            && self.generated.is_empty()
            && self.removed.is_empty()
    }

    // -----------------------------------------------------------------------
    // Writing
    // -----------------------------------------------------------------------

    /// Records that `new` is an image (split piece) of `old`.
    ///
    /// The images table is the builder's naming side table, so it accepts any
    /// shape type (no `BRepTools_History` restriction).
    pub fn add_image(&mut self, old: &TopoShape, new: TopoShape) {
        self.images.entry(ShapeId::of(old)).or_default().push(new);
    }

    /// Records that `new` is modified from `old`.
    ///
    /// A shape cannot be both modified and generated from the same initial: if
    /// `new` was already recorded as a generation of `old`, it is moved to the
    /// modified list (OCCT `BRepTools_History::prepareModified`). A relation
    /// whose list became empty is unbound, keeping the "no empty lists"
    /// invariant of `BRepTools_History`.
    pub fn add_modified(&mut self, old: &TopoShape, new: TopoShape) -> Result<(), String> {
        if !is_supported_type(old) {
            return Err(unsupported_type(old));
        }
        let id = ShapeId::of(old);
        // prepareModified: a generation of `new` cannot coexist with a
        // modification of `new`.
        if let Some(list) = self.generated.get_mut(&id) {
            list.retain(|g| !same_shape(g, &new));
            if list.is_empty() {
                self.generated.remove(&id);
            }
        }
        self.modified.entry(id).or_default().push(new);
        Ok(())
    }

    /// Records that `new` is generated from `old`.
    ///
    /// If `new` was already recorded as a modification of `old`, it is moved
    /// to the generated list (OCCT `BRepTools_History::prepareGenerated`). A
    /// relation whose list became empty is unbound.
    pub fn add_generated(&mut self, old: &TopoShape, new: TopoShape) -> Result<(), String> {
        if !is_supported_type(old) {
            return Err(unsupported_type(old));
        }
        let id = ShapeId::of(old);
        if let Some(list) = self.modified.get_mut(&id) {
            list.retain(|m| !same_shape(m, &new));
            if list.is_empty() {
                self.modified.remove(&id);
            }
        }
        self.generated.entry(id).or_default().push(new);
        Ok(())
    }

    /// Records that `s` is completely removed from the result.
    ///
    /// A removed shape cannot also be modified (OCCT `BRepTools_History::Remove`
    /// unbinds the modified entry); it may still have generated shapes.
    /// Repeated removal of the same shape is a no-op.
    pub fn add_removed(&mut self, s: TopoShape) -> Result<(), String> {
        if !is_supported_type(&s) {
            return Err(unsupported_type(&s));
        }
        let id = ShapeId::of(&s);
        self.modified.remove(&id);
        if self.removed_keys.insert(id) {
            self.removed.push(s);
        }
        Ok(())
    }

    // -----------------------------------------------------------------------
    // Per-shape reading
    // -----------------------------------------------------------------------

    /// Returns the image pieces of `old`, when recorded.
    pub fn image(&self, old: &TopoShape) -> Option<&[TopoShape]> {
        self.images.get(&ShapeId::of(old)).map(|v| v.as_slice())
    }

    /// Returns the shapes modified from `old`, when recorded.
    pub fn modified(&self, old: &TopoShape) -> Option<&[TopoShape]> {
        self.modified.get(&ShapeId::of(old)).map(|v| v.as_slice())
    }

    /// Returns the shapes generated from `old`, when recorded.
    pub fn generated(&self, old: &TopoShape) -> Option<&[TopoShape]> {
        self.generated.get(&ShapeId::of(old)).map(|v| v.as_slice())
    }

    /// Returns true if `old` has been completely removed from the result.
    ///
    /// Mirrors `BRepTools_History::IsRemoved` / `BOPAlgo_BuilderShape::IsDeleted`.
    pub fn is_deleted(&self, old: &TopoShape) -> bool {
        self.removed_keys.contains(&ShapeId::of(old))
    }

    /// BuilderShape-style: true if `old` has recorded images.
    pub fn has_image(&self, old: &TopoShape) -> bool {
        self.images.contains_key(&ShapeId::of(old))
    }

    /// Mutable image list of `old` (`myImages.ChangeSeek`).
    pub fn image_list_mut(&mut self, old: &TopoShape) -> Option<&mut Vec<TopoShape>> {
        self.images.get_mut(&ShapeId::of(old))
    }

    /// BuilderShape-style: true if `old` has modified shapes.
    pub fn has_modified(&self, old: &TopoShape) -> bool {
        self.modified.contains_key(&ShapeId::of(old))
    }

    /// BuilderShape-style: true if `old` has generated shapes.
    pub fn has_generated(&self, old: &TopoShape) -> bool {
        self.generated.contains_key(&ShapeId::of(old))
    }

    // -----------------------------------------------------------------------
    // Whole-table accessors
    // -----------------------------------------------------------------------

    /// Returns the whole images naming table.
    pub fn images(&self) -> &HashMap<ShapeKey, Vec<TopoShape>> {
        &self.images
    }

    /// Returns the whole modified table.
    ///
    /// Named `modified_map` because `modified(old)` is the per-shape query and
    /// Rust does not allow overloading on arity.
    pub fn modified_map(&self) -> &HashMap<ShapeKey, Vec<TopoShape>> {
        &self.modified
    }

    /// Returns the whole generated table.
    ///
    /// Named `generated_map` for the same reason as `modified_map`.
    pub fn generated_map(&self) -> &HashMap<ShapeKey, Vec<TopoShape>> {
        &self.generated
    }

    /// Returns the removed shapes (those the history holds a handle to).
    pub fn removed(&self) -> &[TopoShape] {
        &self.removed
    }

    // -----------------------------------------------------------------------
    // Global predicates (BuilderShape `HasModified` / `HasGenerated` /
    // `HasDeleted`)
    // -----------------------------------------------------------------------

    /// True if any initial shape has modified shapes.
    pub fn has_any_modified(&self) -> bool {
        !self.modified.is_empty()
    }

    /// True if any initial shape has generated shapes.
    pub fn has_any_generated(&self) -> bool {
        !self.generated.is_empty()
    }

    /// True if any shape has been deleted.
    pub fn has_any_deleted(&self) -> bool {
        !self.removed_keys.is_empty()
    }

    /// True if the images table is non-empty.
    pub fn has_any_images(&self) -> bool {
        !self.images.is_empty()
    }

    // -----------------------------------------------------------------------
    // Merge
    // -----------------------------------------------------------------------

    /// Merges `other` (H23, an intermediate T→Q history) into this history
    /// (H12, S→T), producing H13 (S→Q).
    ///
    /// The images naming table is composed separately: each image piece of S1
    /// is mapped forward through `other.images` when present, dropped when the
    /// piece was removed in `other`, and kept otherwise. Keys that exist only
    /// in `other.images` are bound as-is.
    ///
    /// The modified/generated/removed tables follow
    /// `BRepTools_History::Merge`:
    ///
    /// 1. every Tj in the S1 lists is propagated: if removed in `other` it
    ///    disappears from the list; its generations become generations of S1;
    ///    its modifications replace it (as modifications — or as generations,
    ///    when the S1 relation was itself a generation);
    /// 2. T-shapes that carry relations in `other` but were not reached from
    ///    any S1 list are exposed as keys of the merged history (their removed
    ///    markers are cleared);
    /// 3. S1 entries whose lists became empty are unbound and marked removed;
    /// 4. removals of `other` not otherwise accounted for are carried over.
    pub fn merge(&mut self, other: &BopHistory) {
        self.merge_images(other);

        if !(other.has_any_modified() || other.has_any_generated() || other.has_any_deleted()) {
            return;
        }

        // Shapes removed out of an S1 list during phase 1, and T-shapes that
        // got their relations propagated, so phase 4 does not re-process them.
        let mut removed_propagated: HashSet<ShapeKey> = HashSet::new();
        let mut mg_propagated: HashSet<ShapeKey> = HashSet::new();

        // Phase 1 — propagate R23, M23 and G23 into the existing M12/G12
        // tables. `table_is_generated` selects the generated (true) or the
        // modified (false) table.
        for table_is_generated in [true, false] {
            let s1_keys: Vec<ShapeKey> = if table_is_generated {
                self.generated.keys().copied().collect()
            } else {
                self.modified.keys().copied().collect()
            };
            for s1 in s1_keys {
                let mut additions_g: Vec<TopoShape> = Vec::new();
                let mut additions_m: Vec<TopoShape> = Vec::new();
                let mut kept: Vec<TopoShape> = Vec::new();

                let old_list = if table_is_generated {
                    self.generated.get(&s1).cloned().unwrap_or_default()
                } else {
                    self.modified.get(&s1).cloned().unwrap_or_default()
                };

                for s2 in old_list {
                    if other.is_deleted(&s2) {
                        // R23: the intermediate disappeared, drop it from the
                        // S1 list.
                        removed_propagated.insert(ShapeId::of(&s2));
                    } else {
                        let id2 = ShapeId::of(&s2);
                        if let Some(g23) = other.generated.get(&id2) {
                            // G23: these become generations of S1.
                            append_unique(&mut additions_g, g23.iter().cloned());
                            mg_propagated.insert(id2);
                        }
                        if let Some(m23) = other.modified.get(&id2) {
                            // M23: these replace S2 — as modifications of S1,
                            // or as generations when the S1 relation itself
                            // was a generation (aI == 0 in OCCT).
                            if table_is_generated {
                                append_unique(&mut additions_g, m23.iter().cloned());
                            } else {
                                append_unique(&mut additions_m, m23.iter().cloned());
                            }
                            mg_propagated.insert(id2);
                            // S2 is replaced by its modifications: dropped.
                        } else {
                            kept.push(s2);
                        }
                    }
                }

                if table_is_generated {
                    append_unique(&mut kept, additions_g);
                    self.generated.insert(s1, kept);
                } else {
                    append_unique(&mut kept, additions_m);
                    self.modified.insert(s1, kept);
                    // Rule 2: generations of a modified intermediate become
                    // generations of the initial shape.
                    if !additions_g.is_empty() {
                        let gen = self.generated.entry(s1).or_default();
                        append_unique(gen, additions_g);
                    }
                }
            }
        }

        // Phase 2 — expose the H23 keys not reached from any S1 list as keys
        // of the merged history, and clear their removed markers.
        for table_is_generated in [true, false] {
            let other_table = if table_is_generated {
                &other.generated
            } else {
                &other.modified
            };
            let s2_keys: Vec<ShapeKey> = other_table
                .keys()
                .copied()
                .filter(|k| !mg_propagated.contains(k))
                .collect();
            for s2 in s2_keys {
                let values = other_table.get(&s2).cloned().unwrap_or_default();
                let self_table = if table_is_generated {
                    &mut self.generated
                } else {
                    &mut self.modified
                };
                let entry = self_table.entry(s2).or_default();
                append_unique(entry, values);
                self.removed_keys.remove(&s2);
                self.removed.retain(|s| ShapeId::of(s) != s2);
            }
        }

        // Phase 3 — S1 entries whose lists became empty are unbound and the
        // initial shape is marked removed.
        for table_is_generated in [true, false] {
            let keys: Vec<ShapeKey> = if table_is_generated {
                self.generated.keys().copied().collect()
            } else {
                self.modified.keys().copied().collect()
            };
            for s1 in keys {
                let empty = if table_is_generated {
                    self.generated.get(&s1).map(|v| v.is_empty()).unwrap_or(true)
                } else {
                    self.modified.get(&s1).map(|v| v.is_empty()).unwrap_or(true)
                };
                if empty {
                    if table_is_generated {
                        self.generated.remove(&s1);
                    } else {
                        self.modified.remove(&s1);
                    }
                    self.removed_keys.insert(s1);
                }
            }
        }

        // Phase 4 — carry over the removals of `other` that were not already
        // accounted for (not dropped out of an S1 list, not re-exposed as a
        // key in phase 2, not unbound in phase 3).
        for s2 in other.removed.iter() {
            let id = ShapeId::of(s2);
            if !removed_propagated.contains(&id)
                && !self.modified.contains_key(&id)
                && !self.generated.contains_key(&id)
            {
                self.removed_keys.insert(id);
                self.removed.push(s2.clone());
            }
        }
    }

    /// Composes the images naming table with `other`.
    pub(super) fn merge_images(&mut self, other: &BopHistory) {
        let mut composed: HashMap<ShapeKey, Vec<TopoShape>> = HashMap::with_capacity(
            self.images.len() + other.images.len(),
        );
        for (&s1, pieces) in &self.images {
            let mut out: Vec<TopoShape> = Vec::new();
            for s2 in pieces {
                if other.is_deleted(s2) {
                    continue; // the piece no longer exists in H23
                }
                match other.images.get(&ShapeId::of(s2)) {
                    Some(sub) => append_unique(&mut out, sub.iter().cloned()),
                    None => {
                        if !out.iter().any(|d| same_shape(d, s2)) {
                            out.push(s2.clone());
                        }
                    }
                }
            }
            composed.insert(s1, out);
        }
        // Keys that only appear in `other` are bound as-is.
        for (&s2, pieces) in &other.images {
            if !self.images.contains_key(&s2) {
                composed.insert(s2, pieces.clone());
            }
        }
        self.images = composed;
    }

    /// Dumps a brief description of the history, mirroring
    /// `BRepTools_History::Dump`.
    pub fn dump(&self) -> String {
        format!(
            "History contains:\n - {} Deleted shapes;\n - {} Modified shapes;\n - {} Generated shapes.",
            self.removed.len(),
            self.modified.len(),
            self.generated.len()
        )
    }
}

pub(super) fn unsupported_type(s: &TopoShape) -> String {
    format!("history does not support shape type {:?}", s.shape_type())
}
