//! `BOPAlgo_MakerVolume` (`BOPAlgo_MakerVolume.cxx:1-414`,
//! `BOPAlgo_MakerVolume.hxx:1-205`, `BOPAlgo_MakerVolume.lxx:1-78`).
//!
//! Builds solids from a set of shapes that do not bound a closed volume by
//! themselves: the faces of the arguments (split by the intersection phase) are
//! collected, a box that covers the whole set is added, and
//! [`BuilderSolid`](crate::builder_solid::BuilderSolid) builds the areas. The
//! area that touches the covering box is the outside of the arguments and is
//! dropped, the remaining areas are the result.
//!
//! The class derives from `BOPAlgo_Builder` in OCCT, so [`BopBuilder`] is
//! embedded instead of inherited: it owns the [`PaveFiller`], the images table
//! (`myImages`), the history and the report. The `BOPAlgo_Builder` phases the
//! reference calls from `PerformInternal1` (`cxx:102-194`) are the port's
//! [`BopBuilder::fill_images_edges`] family, [`BopBuilder::prepare_history`] and
//! [`BopBuilder::post_treat`].
//!
//! Translation notes:
//! - `TopExp_Explorer` walks become [`explorer_shapes`], `TopoDS_Iterator`
//!   walks [`iter_children`] (`cxx:377-400`, `cxx:484-497`).
//! - `NCollection_List::Remove(iterator)` (`cxx:324`) only deletes the first
//!   element that matches, which is what [`Self::remove_box`] does.

use std::collections::{HashMap, HashSet};

use occt_core::bnd::BndBox;
use occt_core::gp::GpPnt;

use crate::abs::{Orientation, ShapeType};
use crate::bop_builder2::BopBuilder;
use crate::bop_occt_util::{iter_children, shape_key, treat_compound};
use crate::bopds::BopdsDS;
use crate::builder::TopoBuilder;
use crate::builder_solid::BuilderSolid;
use crate::int_tools_full::IntToolsContext;
use crate::pave_filler::{GlueEnum, PaveFiller};
use crate::primitives::BRepPrimBox;
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, ShapeMapKey};

/// `BOPAlgo_MakerVolume` (`BOPAlgo_MakerVolume.hxx:108-205`).
#[derive(Clone)]
pub struct BopMakerVolume {
    /// The `BOPAlgo_Builder` base: filler, images (`myImages`), report.
    builder: BopBuilder,
    /// `myIntersect` (`lxx:18`): `true` by default.
    intersect: bool,
    /// `myBBox` (`hxx:201`).
    bbox: BndBox,
    /// `mySBox` (`hxx:202`).
    sbox: Option<TopoShape>,
    /// `myFaces` (`hxx:203`).
    faces: Vec<TopoShape>,
    /// `myAvoidInternalShapes` (`hxx:204`).
    avoid_internal_shapes: bool,
}

impl Default for BopMakerVolume {
    fn default() -> Self {
        Self::new()
    }
}

impl BopMakerVolume {
    /// `BOPAlgo_MakerVolume()` (`lxx:17-21`).
    pub fn new() -> Self {
        Self {
            builder: BopBuilder::new(),
            intersect: true,
            bbox: BndBox::new(),
            sbox: None,
            faces: Vec::new(),
            avoid_internal_shapes: false,
        }
    }

    /// `Clear()` (`lxx:42-50`).
    pub fn clear(&mut self) {
        self.builder.clear();
        self.intersect = true;
        self.bbox = BndBox::new();
        self.sbox = None;
        self.faces.clear();
        self.avoid_internal_shapes = false;
    }

    /// `SetArguments` (`BOPAlgo_Builder::SetArguments`).
    pub fn set_arguments(&mut self, shapes: &[TopoShape]) {
        self.builder.set_arguments(shapes);
    }

    /// `Arguments()`.
    pub fn arguments(&self) -> &[TopoShape] {
        self.builder.arguments()
    }

    /// `SetFuzzyValue` (`BOPAlgo_Options::SetFuzzyValue`).
    pub fn set_fuzzy_value(&mut self, v: f64) {
        self.builder.set_fuzzy_value(v);
    }

    /// `SetRunParallel` (`BOPAlgo_Options::SetRunParallel`). The port has no
    /// parallel execution, so the flag is not carried.
    pub fn set_run_parallel(&mut self, _b: bool) {}

    /// `SetGlue` (`BOPAlgo_Options::SetGlue`).
    pub fn set_glue(&mut self, glue: GlueEnum) {
        self.builder.set_glue(glue);
    }

    /// `SetNonDestructive` (`BOPAlgo_Options::SetNonDestructive`).
    pub fn set_non_destructive(&mut self, b: bool) {
        self.builder.set_non_destructive(b);
    }

    /// `SetIntersect` (`lxx:53-56`).
    pub fn set_intersect(&mut self, b: bool) {
        self.intersect = b;
    }

    /// `IsIntersect()` (`lxx:63-66`).
    pub fn is_intersect(&self) -> bool {
        self.intersect
    }

    /// `SetAvoidInternalShapes` (`hxx:139-142`).
    pub fn set_avoid_internal_shapes(&mut self, the_avoid_internal: bool) {
        self.avoid_internal_shapes = the_avoid_internal;
    }

    /// `IsAvoidInternalShapes()` (`hxx:145`).
    pub fn is_avoid_internal_shapes(&self) -> bool {
        self.avoid_internal_shapes
    }

    /// `Box()` (`lxx:68-71`): the covering solid box.
    pub fn box_shape(&self) -> Option<&TopoShape> {
        self.sbox.as_ref()
    }

    /// `Faces()` (`lxx:75-78`).
    pub fn faces(&self) -> &[TopoShape] {
        &self.faces
    }

    /// `Shape()` (`BOPAlgo_BuilderShape::Shape`).
    pub fn shape(&self) -> &TopoShape {
        self.builder.result()
    }

    /// `HasErrors()` (`BOPAlgo_Options`).
    pub fn has_errors(&self) -> bool {
        self.builder.has_errors()
    }

    /// The collected fatal alerts.
    pub fn errors(&self) -> &[String] {
        self.builder.errors()
    }

    /// `BOPAlgo_BuilderShape::Modified(theS)` — the history relations the last
    /// run filled in `PrepareHistory`.
    pub fn modified(&self, shape: &TopoShape) -> Vec<TopoShape> {
        self.builder
            .history()
            .modified(shape)
            .map(|s| s.to_vec())
            .unwrap_or_default()
    }

    /// `BOPAlgo_BuilderShape::Images()` - `myHistory->Images()`.
    pub fn images(&self) -> &HashMap<ShapeMapKey, Vec<TopoShape>> {
        self.builder.history().images()
    }

    /// `BOPAlgo_BuilderShape::Origins()` - `myHistory->Origins()`.
    pub fn origins(&self) -> &HashMap<ShapeMapKey, Vec<TopoShape>> {
        self.builder.origins()
    }

    /// `BOPAlgo_BuilderShape::HasDeleted()` - `myHistory->HasDeleted()`.
    pub fn has_deleted(&self) -> bool {
        self.builder.history().has_any_deleted()
    }

    /// `BOPAlgo_BuilderShape::IsDeleted(theS)` - `myHistory->IsDeleted`.
    pub fn is_deleted(&self, shape: &TopoShape) -> bool {
        self.builder.history().is_deleted(shape)
    }

    /// The data structure of the last run (through the filler).
    pub fn ds(&self) -> &BopdsDS {
        self.builder.ds()
    }

    /// `Perform(theRange)` (`cxx:46-98`).
    ///
    /// The intersection phase is set up here (a compound of the arguments when
    /// `myIntersect` is off, `cxx:64-87`), then the building phase runs in
    /// [`Self::perform_internal1`] (`cxx:102-194`).
    pub fn perform(&mut self) -> Result<(), String> {
        // `GetReport()->Clear()` (`cxx:52`).
        let args = self.arguments().to_vec();
        let fuzzy = self.builder.fuzzy_value();
        let non_destructive = self.builder.non_destructive();
        let glue = self.builder.glue();
        let mut filler = PaveFiller::new();
        if !self.intersect {
            // `cxx:64-83`: one compound argument built from every entry.
            let b = TopoBuilder::new();
            let mut comp = b.make_compound_of(&[]);
            for a in &args {
                b.add(&mut comp.0, a);
            }
            filler.set_arguments(&[comp.0]);
        } else {
            filler.set_arguments(&args);
        }
        filler.set_fuzzy_value(fuzzy);
        filler.set_non_destructive(non_destructive);
        filler.set_glue(glue);
        // `pPF->Perform(aPS.Next(anInterPart))` (`cxx:94`).
        let performed = filler.perform();
        // `myEntryPoint = 1; PerformInternal(*pPF, ...)` (`cxx:96-97`).
        self.builder = BopBuilder::new();
        self.builder.set_arguments(&args);
        self.builder.set_filler(filler);
        if let Err(e) = performed {
            for m in self.builder.filler().errors().to_vec() {
                self.builder.add_error(m);
            }
            if !self.builder.has_errors() {
                self.builder.add_error(e.clone());
            }
            return Err(e);
        }
        self.perform_internal1()
    }

    /// `PerformInternal1(thePF, theRange)` (`cxx:102-194`).
    fn perform_internal1(&mut self) -> Result<(), String> {
        // `1. CheckData` (`cxx:110-115`).
        self.check_data()?;

        // `2. Prepare` (`cxx:117-122`): `myShape` is an empty compound
        // (`BOPAlgo_Builder::Prepare`, `BOPAlgo_Builder.cxx:155-162`).
        self.builder
            .set_result_shape(TopoShape::new(ShapeType::Compound));

        // `3. FillImages` (`cxx:127-154`): vertices, edges, wires, faces. The
        // reference does not call the shell / solid container passes here.
        if self.intersect {
            self.builder.fill_images_vertices()?;
            self.builder.fill_images_edges()?;
            crate::bop_build_common::fill_images_containers(&mut self.builder, ShapeType::Wire)?;
            self.builder.fill_images_faces()?;
        }

        // `4. CollectFaces` (`cxx:156-161`).
        self.collect_faces();

        // `5. MakeBox` (`cxx:163-167`).
        let mut box_faces: HashSet<ShapeMapKey> = HashSet::new();
        self.make_box(&mut box_faces)?;

        // `6. BuildSolids` (`cxx:169-174`).
        let mut areas = self.build_solids()?;

        // `7. RemoveBox` (`cxx:176-177`).
        Self::remove_box(&mut areas, &box_faces);

        // `8. FillInternalShapes` (`cxx:179-180`).
        self.fill_internal_shapes(&mut areas);

        // `9. BuildShape` (`cxx:182-183`).
        self.build_shape(&areas);

        // `10. PrepareHistory` (`cxx:185-190`).
        self.builder.prepare_history();
        if self.builder.has_errors() {
            return Err(self.builder.errors().first().cloned().unwrap_or_default());
        }

        // `11. PostTreat` (`cxx:192-193`).
        self.builder.post_treat()
    }

    /// `CheckData()` (`cxx:33-42`): at least one argument, then `CheckFiller`
    /// (`BOPAlgo_Builder.cxx:143-151` — merge the filler report).
    fn check_data(&mut self) -> Result<(), String> {
        if self.arguments().is_empty() {
            let msg = "BOPAlgo_AlertTooFewArguments".to_string();
            self.builder.add_error(msg.clone());
            return Err(msg);
        }
        if self.builder.filler().has_errors() {
            for m in self.builder.filler().errors().to_vec() {
                self.builder.add_error(m);
            }
            return Err(self.builder.errors().first().cloned().unwrap_or_default());
        }
        Ok(())
    }

    /// `CollectFaces()` (`cxx:213-251`).
    fn collect_faces(&mut self) {
        let mut fence: HashSet<ShapeMapKey> = HashSet::new();
        let n = self.builder.ds().nb_source_shapes();
        for i in 0..n {
            let Some(si) = self.builder.ds().shape_info(i) else {
                continue;
            };
            if si.shape_type() != ShapeType::Face {
                continue;
            }
            // `myBBox.Add(aSI.Box())` (`cxx:229-230`).
            if let Some(b) = self.builder.ds().box_of(i) {
                self.bbox.add_box(b);
            }
            let f = si.shape().clone();
            if let Some(imgs) = self.builder.history().image(&f).map(|s| s.to_vec()) {
                // `cxx:233-245`.
                for f_im in imgs {
                    if fence.insert(shape_key(&f_im)) {
                        add_face(&f_im, &mut self.faces);
                    }
                }
            } else {
                // `cxx:246-249`.
                add_face(&f, &mut self.faces);
            }
        }
    }

    /// `MakeBox(theBoxFaces)` (`cxx:255-276`).
    fn make_box(&mut self, box_faces: &mut HashSet<ShapeMapKey>) -> Result<(), String> {
        // `anExt = sqrt(myBBox.SquareExtent()) * 0.5; myBBox.Enlarge(anExt)`.
        let an_ext = self.bbox.square_extent().sqrt() * 0.5;
        self.bbox.enlarge(an_ext);
        // `myBBox.Get(aXmin, ...)` (`cxx:263`).
        let Some((xmin, ymin, zmin, xmax, ymax, zmax)) = self.bbox.get() else {
            return Err("BOPAlgo_MakerVolume: the arguments have no bounding box".to_string());
        };
        // The port's box primitive rejects a degenerate extent; OCCT would
        // build a flat box here. Report the empty box instead of panicking.
        if !(xmax > xmin && ymax > ymin && zmax > zmin) {
            return Err("BOPAlgo_MakerVolume: the bounding box is degenerate".to_string());
        }
        let pmin = GpPnt::new(xmin, ymin, zmin);
        let pmax = GpPnt::new(xmax, ymax, zmax);
        // `mySBox = BRepPrimAPI_MakeBox(aPMin, aPMax).Solid()` (`cxx:267`).
        let sbox = BRepPrimBox::make_box_corner(&pmin, &pmax);
        let sbox_shape = sbox.solid.0.clone();
        // `cxx:269-275`.
        for f in explorer_shapes(&sbox_shape, ShapeType::Face) {
            self.faces.push(f.clone());
            box_faces.insert(shape_key(&f));
        }
        self.sbox = Some(sbox_shape);
        Ok(())
    }

    /// `BuildSolids(theLSR, theRange)` (`cxx:280-298`).
    fn build_solids(&mut self) -> Result<Vec<TopoShape>, String> {
        let mut bs = BuilderSolid::new();
        // `aBS.SetShapes(myFaces); aBS.Perform()` (`cxx:285-288`).
        bs.set_shapes(self.faces.clone());
        bs.set_avoid_internal_shapes(self.avoid_internal_shapes);
        if let Err(e) = bs.perform() {
            // `AddError(new BOPAlgo_AlertSolidBuilderFailed)` (`cxx:289-293`).
            let msg = "BOPAlgo_AlertSolidBuilderFailed".to_string();
            self.builder.add_error(msg.clone());
            let _ = e;
            return Err(msg);
        }
        // `theLSR = aBS.Areas()` (`cxx:297`).
        Ok(bs.areas().to_vec())
    }

    /// `RemoveBox(theLSR, theBoxFaces)` (`cxx:302-333`).
    ///
    /// Removes the first area that contains a face of the covering box, then
    /// stops — the reference's `break` leaves the other areas (`cxx:328-331`).
    fn remove_box(areas: &mut Vec<TopoShape>, box_faces: &HashSet<ShapeMapKey>) {
        for i in 0..areas.len() {
            let mut found = false;
            for f in explorer_shapes(&areas[i], ShapeType::Face) {
                if box_faces.contains(&shape_key(&f)) {
                    found = true;
                    break;
                }
            }
            if found {
                areas.remove(i);
                break;
            }
        }
    }

    /// `FillInternalShapes(theLSR)` (`cxx:359-403`).
    fn fill_internal_shapes(&mut self, areas: &mut [TopoShape]) {
        if self.avoid_internal_shapes {
            return;
        }

        // `TreatCompound` over the arguments (`cxx:366-375`).
        let mut lsc: Vec<TopoShape> = Vec::new();
        let mut fence: HashSet<ShapeMapKey> = HashSet::new();
        for a in self.builder.ds().arguments().to_vec() {
            treat_compound(&a, &mut lsc, &mut fence);
        }

        // `cxx:377-400`: only edges and vertices, a wire contributing its edges.
        let mut lve: Vec<TopoShape> = Vec::new();
        for s in &lsc {
            match s.shape_type() {
                ShapeType::Wire => {
                    for ss in iter_children(s) {
                        if fence.insert(shape_key(&ss)) {
                            lve.push(ss);
                        }
                    }
                }
                ShapeType::Vertex | ShapeType::Edge => lve.push(s.clone()),
                _ => {}
            }
        }

        // `BOPAlgo_Tools::FillInternals(theLSR, aLVE, myImages, myContext)`
        // (`cxx:402`).
        let images = self.builder.history().images().clone();
        let ctx = IntToolsContext::new();
        crate::bopalgo_tools_class::fill_internals(areas, &lve, &images, &ctx);
    }

    /// `BuildShape(theLSR)` (`cxx:337-355`).
    fn build_shape(&mut self, areas: &[TopoShape]) {
        if areas.len() == 1 {
            self.builder.set_result_shape(areas[0].clone());
            return;
        }
        let b = TopoBuilder::new();
        let mut result = self.builder.result().clone();
        for sol in areas {
            b.add(&mut result, sol);
        }
        self.builder.set_result_shape(result);
    }
}

/// The static `AddFace` (`cxx:407-414`): a face enters `myFaces` twice, once
/// FORWARD and once REVERSED.
fn add_face(f: &TopoShape, list: &mut Vec<TopoShape>) {
    let mut ff = f.clone();
    ff.set_orientation(Orientation::Forward);
    list.push(ff.clone());
    ff.set_orientation(Orientation::Reversed);
    list.push(ff);
}
