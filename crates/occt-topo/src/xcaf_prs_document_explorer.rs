//! `XCAFPrs_DocumentExplorer` over the [`BRepModel`] label tree.
//!
//! Source: `XCAFPrs_DocumentExplorer.hxx:29-183`,
//! `XCAFPrs_DocumentExplorer.cxx:28-437`, `XCAFPrs_DocumentNode.hxx:25-41`.
//!
//! OCCT walks a `TDocStd_Document`: every node is a `TDF_Label`, its Product is
//! the label `GetReferredShape` resolves to, and the children are the labels
//! `TDF_ChildIterator(RefLabel)` visits. This crate has no TDF document; a
//! label is an index into [`BRepModel::labels`] and the child iterator is the
//! label's `children` vector, so the iteration order and the node contents are
//! the same.
//!
//! Locations are not modelled: a shape's placement is baked into its geometry
//! when the document is filled (`rwmesh::gltf_document`), so `Location` and
//! `LocalTrsf` are always the identity transformation. `DefineChildId` uses the
//! label index in place of `TDF_Tool::Entry` (`cxx:88-93`), keeping the
//! `parent/child.` shape of the identifier the OBJ label filter matches on.
//!
//! UNPORTED: `FindLabelFromPathId` / `FindShapeFromPathId` (`cxx:97-160`, no
//! caller in this crate) and the `XCAFDoc_ColorTool` `ColorGen` role — a label
//! stores the `ColorSurf` / `ColorCurv` roles only (`model.rs`).

use occt_core::gp::GpTrsf;
use occt_core::quantity::ColorRGBA;

use crate::model::BRepModel;
use crate::xcaf_prs::XcafPrsStyle;

/// `XCAFPrs_DocumentExplorerFlags` (`hxx:32-37`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct DocumentExplorerFlags(u32);

impl DocumentExplorerFlags {
    pub const NONE: Self = Self(0x00);
    /// `XCAFPrs_DocumentExplorerFlags_OnlyLeafNodes`: skip assembly nodes.
    pub const ONLY_LEAF_NODES: Self = Self(0x01);
    /// `XCAFPrs_DocumentExplorerFlags_NoStyle`: do not fetch styles.
    pub const NO_STYLE: Self = Self(0x02);

    pub fn has(self, flag: Self) -> bool {
        self.0 & flag.0 != 0
    }
}

/// `XCAFPrs_DocumentNode` (`hxx:25-41`): one explored document node.
#[derive(Debug, Clone, Default)]
pub struct DocumentNode {
    /// `Id` — the `DefineChildId` path.
    pub id: String,
    /// `Label` — the label in the document.
    pub label: usize,
    /// `RefLabel` — the Product label `GetReferredShape` resolves to.
    pub ref_label: usize,
    /// `Style` — the merged `XCAFPrs_Style`.
    pub style: XcafPrsStyle,
    /// `Location` — the node's global transformation (identity here).
    pub location: GpTrsf,
    /// `LocalTrsf` — the node's transformation relative to its parent (identity
    /// here).
    pub local_trsf: GpTrsf,
    /// `IsAssembly` — the label is an assembly.
    pub is_assembly: bool,
}

/// `XCAFPrs_DocumentNode::ChildIter` plus the node it belongs to.
struct StackNode {
    node: DocumentNode,
    /// Position in `labels[node.ref_label].children`.
    child: usize,
}

/// `XCAFPrs_DocumentExplorer` (`hxx:39-183`): depth-first iteration over the
/// document's root labels.
pub struct DocumentExplorer<'a> {
    model: &'a BRepModel,
    roots: Vec<usize>,
    root_index: usize,
    /// `myNodeStack`; `myTop` is `stack.len() - 1`.
    stack: Vec<StackNode>,
    current: DocumentNode,
    has_more: bool,
    flags: DocumentExplorerFlags,
    def_style: XcafPrsStyle,
}

impl<'a> DocumentExplorer<'a> {
    /// `XCAFPrs_DocumentExplorer(document, theFlags, theDefStyle)`
    /// (`cxx:173-185`): explore the document's free shapes (`GetFreeShapes`).
    pub fn new(model: &'a BRepModel, flags: DocumentExplorerFlags, def_style: XcafPrsStyle) -> Self {
        let roots = model.roots.clone();
        Self::with_roots(model, &roots, flags, def_style)
    }

    /// `XCAFPrs_DocumentExplorer(document, theRoots, ...)` (`cxx:188-197`).
    pub fn with_roots(
        model: &'a BRepModel,
        roots: &[usize],
        flags: DocumentExplorerFlags,
        def_style: XcafPrsStyle,
    ) -> Self {
        let mut explorer = Self {
            model,
            roots: roots.to_vec(),
            root_index: 0,
            stack: Vec::new(),
            current: DocumentNode::default(),
            has_more: false,
            flags,
            def_style,
        };
        explorer.init_root();
        explorer
    }

    /// `More()` (`hxx:127`).
    pub fn more(&self) -> bool {
        self.has_more
    }

    /// `Current()` (`hxx:130`).
    pub fn current(&self) -> &DocumentNode {
        &self.current
    }

    /// `CurrentDepth()` (`hxx:150-152`): 0 is a root label.
    pub fn current_depth(&self) -> usize {
        if self.current.is_assembly {
            self.stack.len().saturating_sub(1)
        } else {
            self.stack.len()
        }
    }

    /// `Current(theDepth)` (`hxx:133-148`): the node at `depth`, the current one
    /// included.
    pub fn current_at(&self, depth: usize) -> Option<&DocumentNode> {
        if depth == self.current_depth() {
            return Some(&self.current);
        }
        self.stack.get(depth).map(|s| &s.node)
    }

    /// `Init` -> `initRoot` (`cxx:214-278`).
    fn init_root(&mut self) {
        loop {
            // "reset the stack" (`cxx:242-248`).
            self.stack.clear();
            if self.root_index >= self.roots.len() {
                self.has_more = false;
                self.init_current(false);
                return;
            }
            let root = self.roots[self.root_index];
            if root >= self.model.labels.len() {
                self.root_index += 1;
                continue;
            }
            self.has_more = true;
            let ref_label = self.referred_or_self(root);
            if self.model.is_assembly(ref_label) {
                self.next();
            } else {
                self.init_current(false);
            }
            return;
        }
    }

    /// `initCurrent(theIsAssembly)` (`cxx:284-322`).
    fn init_current(&mut self, is_assembly: bool) {
        self.current = DocumentNode::default();
        if is_assembly {
            if let Some(top) = self.stack.last() {
                self.current = top.node.clone();
            }
        } else if let Some(parent) = self.stack.last() {
            let parent_node = parent.node.clone();
            let Some(&label) = self.model.labels[parent_node.ref_label].children.get(parent.child) else {
                return;
            };
            let ref_label = self.referred_or_self(label);
            self.current = DocumentNode {
                label,
                ref_label,
                style: self.merged_style(&parent_node.style, label, ref_label),
                location: GpTrsf::identity(),
                local_trsf: GpTrsf::identity(),
                id: define_child_id(label, &parent_node.id),
                is_assembly: false,
            };
        } else {
            if self.root_index >= self.roots.len() {
                return;
            }
            let label = self.roots[self.root_index];
            let ref_label = self.referred_or_self(label);
            self.current = DocumentNode {
                label,
                ref_label,
                style: self.merged_style(&self.def_style, label, ref_label),
                location: GpTrsf::identity(),
                local_trsf: GpTrsf::identity(),
                id: define_child_id(label, ""),
                is_assembly: false,
            };
        }
    }

    /// `Next()` (`cxx:328-437`).
    pub fn next(&mut self) {
        if !self.has_more {
            return;
        }
        if self.stack.is_empty() {
            let root = self.roots[self.root_index];
            let ref_label = self.referred_or_self(root);
            if !self.model.is_assembly(ref_label) {
                // "already visited once" (`cxx:338-344`).
                self.root_index += 1;
                self.init_root();
                return;
            }
            // "push and try to find" (`cxx:346-364`).
            let node = self.make_node(root, ref_label, &self.def_style, "", true);
            self.stack.push(StackNode { node, child: 0 });
            if !self.flags.has(DocumentExplorerFlags::ONLY_LEAF_NODES) {
                self.init_current(true);
                return;
            }
        } else if !self.current.is_assembly {
            if let Some(top) = self.stack.last_mut() {
                top.child += 1;
            }
        }

        loop {
            let Some(top) = self.stack.len().checked_sub(1) else {
                self.root_index += 1;
                self.init_root();
                return;
            };
            let ref_label = self.stack[top].node.ref_label;
            let children_len = self.model.labels.get(ref_label).map(|l| l.children.len()).unwrap_or(0);
            if self.stack[top].child < children_len {
                let label = self.model.labels[ref_label].children[self.stack[top].child];
                let child_ref = self.referred_or_self(label);
                if !self.model.is_assembly(child_ref) {
                    self.has_more = true;
                    self.init_current(false);
                    return;
                }
                // `cxx:394-418`: descend into the assembly.
                let parent = self.stack[top].node.clone();
                let node = self.make_node(label, child_ref, &parent.style, &parent.id, true);
                self.stack.push(StackNode { node, child: 0 });
                if !self.flags.has(DocumentExplorerFlags::ONLY_LEAF_NODES) {
                    self.has_more = true;
                    self.init_current(true);
                    return;
                }
            } else {
                // `cxx:423-435`: the child iterator is exhausted.
                self.stack.pop();
                if self.stack.is_empty() {
                    self.root_index += 1;
                    self.init_root();
                    return;
                }
                if let Some(new_top) = self.stack.last_mut() {
                    new_top.child += 1;
                }
            }
        }
    }

    fn make_node(
        &self,
        label: usize,
        ref_label: usize,
        parent_style: &XcafPrsStyle,
        parent_id: &str,
        is_assembly: bool,
    ) -> DocumentNode {
        DocumentNode {
            label,
            ref_label,
            style: self.merged_style(parent_style, label, ref_label),
            location: GpTrsf::identity(),
            local_trsf: GpTrsf::identity(),
            id: define_child_id(label, parent_id),
            is_assembly,
        }
    }

    /// `GetReferredShape(label).unwrap_or(label)`: the Product label of a
    /// reference, the label itself otherwise (`XCAFDoc_ShapeTool.cxx:939-950`).
    fn referred_or_self(&self, label: usize) -> usize {
        self.model.referred_label(label).unwrap_or(label)
    }

    /// `mergedStyle` (`cxx:28-81`): the parent style, then the Product label's
    /// style, then the Instance label's style.
    fn merged_style(&self, parent: &XcafPrsStyle, label: usize, ref_label: usize) -> XcafPrsStyle {
        if self.flags.has(DocumentExplorerFlags::NO_STYLE) {
            return parent.clone();
        }
        let mut style = parent.clone();
        apply_label_style(self.model, ref_label, &mut style);
        if label != ref_label {
            apply_label_style(self.model, label, &mut style);
        }
        style
    }
}

/// `DefineChildId` (`cxx:88-93`): `<parent>/<entry>.`, the entry being the
/// label index here.
pub fn define_child_id(label: usize, parent_id: &str) -> String {
    if parent_id.is_empty() {
        format!("{label}.")
    } else {
        format!("{parent_id}{label}.")
    }
}

/// The `theRefLabel` half of `mergedStyle` (`cxx:40-57`): material, then the
/// `ColorGen` / `ColorSurf` / `ColorCurv` roles. This crate's labels carry the
/// `ColorSurf` role (`ModelLabel::color`) and `ColorCurv`
/// (`ModelLabel::color_curv`) only, so the `ColorGen` arm has no source.
fn apply_label_style(model: &BRepModel, label: usize, style: &mut XcafPrsStyle) {
    let Some(l) = model.label(label) else { return };
    if let Some(mat) = l.material.and_then(|i| model.get_material(i)) {
        style.set_material(Some(mat.clone()));
    }
    if let Some(c) = l.color {
        style.set_color_surf(ColorRGBA { r: c.r, g: c.g, b: c.b, a: 1.0 });
    }
    if let Some(c) = l.color_curv {
        style.set_color_curv(c);
    }
}
