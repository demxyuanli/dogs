//! `BRepFeat_Builder` -- the local Boolean of the form features.
//!
//! Source: `ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_Builder.hxx:1-131`,
//! `.cxx:1-827`.
//!
//! The class is a `BOPAlgo_BOP` whose result is rebuilt from the kept parts of
//! the tool: `KeepPart` fills the `myShapes` map, `Prepare` + `RebuildFaces`
//! rewrite the face images so that only the retained edges bound them, and
//! `PerformResult` runs the remaining building stages.
//!
//! Module table (`BRepFeat` = `src/ModelingAlgorithms/TKFeat/BRepFeat`):
//!
//! | Reference | Rust |
//! |---|---|
//! | ctor / `Clear` (`cxx:45-61`) | [`BRepFeatBuilder::new`] / `clear` |
//! | both `Init` (`cxx:65-81`) | [`BRepFeatBuilder::init`] / [`BRepFeatBuilder::init_with_tool`] |
//! | both `SetOperation` (`cxx:85-103`) | [`BRepFeatBuilder::set_operation`] / [`BRepFeatBuilder::set_operation_with_flag`] |
//! | inherited `AddTool` (`BOPAlgo_ToolsProvider.cxx:47-53`) | [`BRepFeatBuilder::add_tool`] |
//! | `PartsOfTool` (`cxx:107-118`) | [`BRepFeatBuilder::parts_of_tool`] |
//! | `KeepParts` / `KeepPart` (`cxx:122-142`) | [`BRepFeatBuilder::keep_parts`] / [`BRepFeatBuilder::keep_part`] |
//! | `Prepare` (`cxx:146-156`) | `prepare` |
//! | `FillRemoved` both overloads (`cxx:158-187`, `cxx:781-797`) | `fill_removed` / `fill_removed_rec` |
//! | `PerformResult` (`cxx:192-267`) | [`BRepFeatBuilder::perform_result`] |
//! | `RebuildFaces` (`cxx:271-558`) | `rebuild_faces.rs` |
//! | `RebuildEdge` (`cxx:562-722`) | `rebuild_edge.rs` |
//! | `CheckSolidImages` (`cxx:726-777`) | `check_solid_images` |
//! | `FillIn3DParts` (`cxx:801-826`) | not ported, see the deviations |
//!
//! Deviations:
//! - the reference inherits `BOPAlgo_BOP`; the port owns a
//!   [`crate::bop_builder2::BopBuilder`] and forwards the argument / tool
//!   groups to it. `BOPAlgo_BOP::Perform` is `perform`, and
//!   `BOPAlgo_BOP::BuildShape` is `perform_result`.
//! - `BRepFeat_Builder::FillIn3DParts` (`cxx:801-826`) drops the removed faces
//!   from `myInParts` *between* `BOPAlgo_Builder::FillIn3DParts` and
//!   `BOPAlgo_BOP::BuildShape`. The port's `BopBuilder::build_bop` runs both
//!   stages inside one call and exposes no hook between them, so the override
//!   is **UNPORTED** (`BRepFeat_Builder.cxx:801-826`); the removed faces stay
//!   in `myInParts`.
//! - `CheckArgsForOpenSolid -> false` (`hxx:119-121`) is applied through
//!   `BopBuilder::set_check_args_for_open_solid` when the object is made.

mod rebuild_edge;
mod rebuild_faces;
mod shape_set;

use crate::abs::ShapeType;
use crate::bop_build_common::fill_images_containers;
use crate::bop_builder2::{BopBuilder, BoolOp2};
use crate::shape::TopoShape;
use crate::topo_tools::shapes_of;

pub use shape_set::ShapeSet;

use rebuild_edge::rebuild_edge;

/// `BRepFeat_Builder` (`hxx:47-127`).
#[derive(Debug, Clone)]
pub struct BRepFeatBuilder {
    /// The inherited `BOPAlgo_BOP`.
    b: BopBuilder,
    /// `myShapes` (`hxx:123`): the kept sub-shapes of the tool.
    pub(super) my_shapes: ShapeSet,
    /// `myRemoved` (`hxx:124`): the tool sub-shapes that were dropped.
    pub(super) my_removed: ShapeSet,
    /// `myFuse` (`hxx:125`).
    my_fuse: i32,
    /// `BOPAlgo_BOP::myOperation`, as the port carries it.
    my_operation: BoolOp2,
    /// `myArguments` of the object group (`BOPAlgo_BOP::AddArgument`).
    my_arguments: Vec<TopoShape>,
    /// `myTools` (`BOPAlgo_BOP::AddTool`).
    my_tools: Vec<TopoShape>,
}

impl Default for BRepFeatBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepFeatBuilder {
    /// `BRepFeat_Builder()` (`cxx:45-48`): `Clear()`.
    pub fn new() -> Self {
        let mut b = BopBuilder::new();
        // `hxx:119-121`: `CheckArgsForOpenSolid()` always answers false.
        b.set_check_args_for_open_solid(false);
        Self {
            b,
            my_shapes: ShapeSet::new(),
            my_removed: ShapeSet::new(),
            my_fuse: 0,
            my_operation: BoolOp2::Cut,
            my_arguments: Vec::new(),
            my_tools: Vec::new(),
        }
    }

    /// `Clear()` (`cxx:56-61`).
    pub fn clear(&mut self) {
        self.my_shapes.clear();
        self.my_removed.clear();
        self.my_arguments.clear();
        self.my_tools.clear();
        self.b.clear();
        self.b.set_check_args_for_open_solid(false);
    }

    /// `Init(theShape)` (`cxx:65-72`).
    pub fn init(&mut self, the_shape: &TopoShape) {
        self.clear();
        self.my_arguments.push(the_shape.clone());
        self.sync_arguments();
    }

    /// `Init(theShape, theTool)` (`cxx:76-81`).
    pub fn init_with_tool(&mut self, the_shape: &TopoShape, the_tool: &TopoShape) {
        self.clear();
        self.my_arguments.push(the_shape.clone());
        self.my_tools.push(the_tool.clone());
        self.sync_arguments();
    }

    /// `BOPAlgo_ToolsProvider::AddTool(theShape)` (`BOPAlgo_ToolsProvider.cxx:47-53`):
    /// append to `myTools` when `myMapTools` does not hold it yet. The map
    /// hasher of `TopTools_ShapeMapHasher` keys on TShape + Location, which is
    /// [`crate::topo_tools::is_same_occurrence`].
    pub fn add_tool(&mut self, the_tool: &TopoShape) {
        if !self
            .my_tools
            .iter()
            .any(|t| crate::topo_tools::is_same_occurrence(t, the_tool))
        {
            self.my_tools.push(the_tool.clone());
        }
        self.sync_arguments();
    }

    /// `myArguments.First()` (`NCollection_List::First`) of the object group.
    pub fn first_argument(&self) -> Option<&TopoShape> {
        self.my_arguments.first()
    }

    /// `myArguments` (`BOPAlgo_Builder::Arguments`): the object group only.
    pub fn arguments(&self) -> &[TopoShape] {
        &self.my_arguments
    }

    /// `SetOperation(theFuse)` (`cxx:85-89`): `0` is CUT, anything else FUSE.
    pub fn set_operation(&mut self, the_fuse: i32) {
        self.my_fuse = the_fuse;
        self.my_operation = if self.my_fuse != 0 {
            BoolOp2::Fuse
        } else {
            BoolOp2::Cut
        };
    }

    /// `SetOperation(theFuse, theFlag)` (`cxx:93-103`): with `theFlag` the
    /// operation is COMMON (cut) or CUT21 (fuse), since no part of the tool is
    /// selected.
    pub fn set_operation_with_flag(&mut self, the_fuse: i32, the_flag: bool) {
        self.my_fuse = the_fuse;
        self.my_operation = if !the_flag {
            if self.my_fuse != 0 {
                BoolOp2::Fuse
            } else {
                BoolOp2::Cut
            }
        } else if self.my_fuse != 0 {
            BoolOp2::Cut21
        } else {
            BoolOp2::Common
        };
    }

    /// `PartsOfTool(aLT)` (`cxx:107-118`): the solids of the built result.
    pub fn parts_of_tool(&self) -> Vec<TopoShape> {
        shapes_of(self.b.result(), ShapeType::Solid)
    }

    /// `KeepParts(theIm)` (`cxx:122-132`).
    pub fn keep_parts(&mut self, the_im: &[TopoShape]) {
        for t in the_im {
            self.keep_part(t);
        }
    }

    /// `KeepPart(thePart)` (`cxx:136-142`): `TopExp::MapShapes(thePart,
    /// myShapes)`.
    pub fn keep_part(&mut self, the_part: &TopoShape) {
        for s in explorer_all(the_part) {
            self.my_shapes.add(&s);
        }
    }

    /// `Perform()` -- the inherited `BOPAlgo_BOP::Perform` (`cxx` of
    /// `BOPAlgo_BOP.cxx:364-409`): the intersection phase, then the building
    /// phase with `myOperation`.
    pub fn perform(&mut self) -> Result<(), String> {
        self.sync_arguments();
        // Intersection + face images (`PerformInternal1` up to BuildResult).
        self.b.perform()?;
        // `BOPAlgo_BOP::BuildShape`.
        self.b
            .build_bop_operation(&self.my_arguments, &self.my_tools, self.my_operation)
    }

    /// `Prepare()` (`cxx:146-156`): the report is cleared, `myShape` starts as
    /// an empty compound and the removed tool parts are collected.
    pub fn prepare(&mut self) {
        self.fill_removed();
    }

    /// `FillRemoved()` (`cxx:158-187`).
    ///
    /// Every solid of the object group loses its images, then the images of the
    /// tool group are walked and every sub-shape that is not kept is collected
    /// in `myRemoved`.
    pub fn fill_removed(&mut self) {
        let args0 = match self.my_arguments.first() {
            Some(s) => s.clone(),
            None => return,
        };
        let args1 = match self.my_tools.first() {
            Some(s) => s.clone(),
            None => return,
        };
        for solid in shapes_of(&args0, ShapeType::Solid) {
            self.b.history_mut().unbind_image(&solid);
        }
        let list = match self.b.history().image(&args1) {
            Some(v) => v.to_vec(),
            None => return,
        };
        for s in list {
            fill_removed_rec(&s, &self.my_shapes, &mut self.my_removed);
        }
    }

    /// `PerformResult()` (`cxx:192-267`): the kept parts of the tool drive the
    /// face rebuild, then the Boolean result is assembled.
    ///
    /// `BuildResult(SOLID)` / `FillImagesCompounds` / `BuildResult(COMPOUND)` /
    /// `BuildShape` of the reference are the single
    /// [`BopBuilder::build_bop_operation`] call here: the port's `BuildBOP`
    /// runs the same assembly in one entry point.
    pub fn perform_result(&mut self) -> Result<(), String> {
        // `cxx:194`.
        self.my_operation = self.fuse_operation();
        if self.my_shapes.is_empty() {
            return self
                .b
                .build_bop_operation(&self.my_arguments, &self.my_tools, self.my_operation);
        }

        // `cxx:230`.
        self.prepare();
        // `cxx:233`.
        self.rebuild_faces()?;
        // `cxx:236`.
        fill_images_containers(&mut self.b, ShapeType::Shell)?;
        if self.b.has_errors() {
            return Err(self.first_error());
        }
        // `cxx:238`.
        self.b.fill_images_solids()?;
        if self.b.has_errors() {
            return Err(self.first_error());
        }
        // `cxx:245`.
        self.check_solid_images();
        // `cxx:249-267`: `BuildResult(SOLID)` and the compound stages.
        self.b
            .build_bop_operation(&self.my_arguments, &self.my_tools, self.my_operation)
    }

    /// `RebuildFaces()` (`cxx:271-558`), see `rebuild_faces.rs`.
    pub fn rebuild_faces(&mut self) -> Result<(), String> {
        rebuild_faces::rebuild_faces(self)
    }

    /// `RebuildEdge(theE, theF, aME, aLIm)` (`cxx:562-722`).
    pub fn rebuild_edge(
        &mut self,
        the_e: &TopoShape,
        the_f: &crate::shape::Face,
        a_me: &ShapeSet,
        a_lim: &mut Vec<TopoShape>,
    ) -> Result<(), String> {
        rebuild_edge(&mut self.b, the_e, the_f, a_me, a_lim)
    }

    /// `CheckSolidImages()` (`cxx:726-777`).
    ///
    /// Every image of an object solid is replaced by the same-domain
    /// representative among the tool images. The reference keeps `aLSImNew`
    /// across the outer loop without clearing it, so later solids inherit the
    /// images of the earlier ones; that is reproduced as written.
    pub fn check_solid_images(&mut self) {
        let args0 = match self.my_arguments.first() {
            Some(s) => s.clone(),
            None => return,
        };
        let args1 = match self.my_tools.first() {
            Some(s) => s.clone(),
            None => return,
        };
        let mut a_mst = crate::bop_tools_set::BopToolsSetMap::new();
        if let Some(list) = self.b.history().image(&args1) {
            for sol_im in list {
                a_mst.add(crate::bop_tools_set::BopToolsSet::from_shape(
                    sol_im,
                    ShapeType::Face,
                ));
            }
        }
        let mut new_images: Vec<TopoShape> = Vec::new();
        for solid in shapes_of(&args0, ShapeType::Solid) {
            let images = match self.b.history().image(&solid) {
                Some(v) => v.to_vec(),
                None => continue,
            };
            for sol_im in images {
                let set =
                    crate::bop_tools_set::BopToolsSet::from_shape(&sol_im, ShapeType::Face);
                let flag_sd = a_mst.contains(&set);
                let sx = a_mst.added(set).shape().clone();
                new_images.push(sx.clone());
                if flag_sd {
                    self.b.bind_shapes_sd(sol_im, sx);
                }
            }
            self.b.history_mut().set_images(&solid, new_images.clone());
        }
    }

    /// `Modified(S)` of `BOPAlgo_BuilderShape`: the images of `S`.
    pub fn modified(&self, s: &TopoShape) -> Vec<TopoShape> {
        self.b
            .history()
            .modified(s)
            .map(|v| v.to_vec())
            .unwrap_or_default()
    }

    /// `IsDeleted(S)` of `BOPAlgo_BuilderShape`.
    pub fn is_deleted(&self, s: &TopoShape) -> bool {
        self.b.history().is_deleted(s)
    }

    /// `Shape()` of the last run.
    pub fn shape(&self) -> &TopoShape {
        self.b.result()
    }

    /// `HasErrors` of the last run.
    pub fn has_errors(&self) -> bool {
        self.b.has_errors()
    }

    /// `myFuse ? FUSE : CUT` (`cxx:194`).
    fn fuse_operation(&self) -> BoolOp2 {
        if self.my_fuse != 0 {
            BoolOp2::Fuse
        } else {
            BoolOp2::Cut
        }
    }

    /// `BOPAlgo_BOP::myArguments` / `myTools` written back to the builder.
    fn sync_arguments(&mut self) {
        let mut all = self.my_arguments.clone();
        for t in &self.my_tools {
            if !all.iter().any(|a| a.same_tshape(t)) {
                all.push(t.clone());
            }
        }
        self.b.set_arguments(&all);
    }

    /// The builder of the inherited half.
    pub(super) fn bop(&mut self) -> &mut BopBuilder {
        &mut self.b
    }

    /// The first error of the last run (`Standard_OutOfRange` on an empty
    /// report in the reference).
    fn first_error(&self) -> String {
        self.b
            .errors()
            .first()
            .cloned()
            .unwrap_or_else(|| "BOPAlgo_AlertBuilderFailed".to_string())
    }
}

/// `TopExp::MapShapes(S)` equivalent: `S` and every sub-shape, in explorer
/// order, without duplicates.
fn explorer_all(s: &TopoShape) -> Vec<TopoShape> {
    let mut out: Vec<TopoShape> = Vec::new();
    collect_shapes(s, &mut out);
    out
}

fn collect_shapes(s: &TopoShape, out: &mut Vec<TopoShape>) {
    if !out.iter().any(|x| x.same_tshape(s)) {
        out.push(s.clone());
    }
    for c in children_of(s) {
        collect_shapes(&c, out);
    }
}

/// `TopoDS_Iterator` children of `S` in stored order.
fn children_of(s: &TopoShape) -> Vec<TopoShape> {
    crate::bop_occt_util::iter_children(s)
}

/// `FillRemoved(S, M)` (`cxx:783-796`): a kept shape stops the walk, every
/// other shape goes into `M` before its children are visited.
fn fill_removed_rec(s: &TopoShape, my_shapes: &ShapeSet, map: &mut ShapeSet) {
    if my_shapes.contains(s) {
        return;
    }
    map.add(s);
    for c in children_of(s) {
        fill_removed_rec(&c, my_shapes, map);
    }
}
