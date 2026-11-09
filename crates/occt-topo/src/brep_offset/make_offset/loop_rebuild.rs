//! `BRepOffset_MakeOffset::MakeLoops` (`BRepOffset_MakeOffset.cxx:2991-3057`)
//! and `MakeFaces` (`:3059-3146`): the 2D unwinding and the reconstruction of
//! the untouched faces that share a rebuilt edge.
//!
//! `MakeLoops` (`:2991-3057`) splits `Modif` into the offset faces (`LF`, the
//! ones that are **not** caps) and the caps (`LC`), rebuilds the loops of `LF`
//! with `myMakeLoops.Build` (`:3024`) - or with `BuildSplitsOfTrimmedFaces`
//! (`:3020`) on the planar `Complete` arm - and then the loops of the caps on
//! their context with `myMakeLoops.BuildOnContext(LC, myAnalyse, myAsDes,
//! myImageOffset, InSide)` (`:3048`), where `InSide` is `false` for a
//! non-negative offset (`:3041-3045`).
//!
//! `MakeFaces` (`:3059-3146`) collects the first image of every
//! `myInitOffsetFace` root whose image has no entry in `myImageOffset`
//! (`:3078-3084`) and hands that list to `myMakeLoops.BuildFaces` (`:3093`) -
//! or to `BuildSplitsOfTrimmedFaces` (`:3089`) on the same planar arm.
//!
//! UNPORTED:
//!
//! * `Message_ProgressScope` / `!More()` -> `BRepOffset_UserBreak`
//!   (`:3027-3031`, `:3095-3099`).
//! * the commented-out `NCollection_Map` iterators (`:2957`, `:3004`) and the
//!   `OCCT_DEBUG` blocks (`:2995-3001`, `:3049-3051`, `:3063-3071`,
//!   `:3101-3103`).

use occt_core::kernel::geomabs::JoinType;

use crate::brep_offset::make_offset::build_offset_faces::build_splits_of_trimmed_faces;
use crate::shape::{Face, TopoShape};

use super::BRepOffsetMakeOffset;

impl BRepOffsetMakeOffset {
    /// `MakeLoops(Modif, theRange)` (`cxx:2991-3057`).
    pub fn make_loops(&mut self, modif: &crate::brep_offset::inter3d::ShapeIndexedMap) {
        // `cxx:3005-3014`: the offset faces, caps excluded.
        let mut lf: Vec<TopoShape> = Vec::new();
        for i in 1..=modif.extent() {
            let Some(shape) = modif.value(i) else {
                continue;
            };
            if !self.my_faces.contains(shape) {
                lf.push(shape.clone());
            }
        }

        // `cxx:3018-3025`.
        if self.my_join == JoinType::Intersection && self.my_inter && self.my_is_planar {
            build_splits_of_trimmed_faces(&lf, &self.my_as_des, &mut self.my_image_offset);
        } else {
            self.my_make_loops.build(
                &lf,
                &mut self.my_as_des,
                &mut self.my_image_offset,
                &self.my_image_vv,
            );
        }

        // `cxx:3036-3038`: the caps.
        let mut lc: Vec<TopoShape> = Vec::new();
        for i in 1..=self.my_faces.extent() {
            if let Some(shape) = self.my_faces.value(i) {
                lc.push(shape.clone());
            }
        }

        // `cxx:3040-3045`.
        let mut in_side = true;
        if self.my_offset > 0.0 {
            in_side = false;
        }

        // `cxx:3048`.
        self.my_make_loops.build_on_context(
            &lc,
            &self.my_analyse,
            &mut self.my_as_des,
            &mut self.my_image_offset,
            in_side,
        );
    }

    /// `MakeFaces(Modif, theRange)` (`cxx:3059-3146`).
    pub fn make_faces(&mut self) {
        // `cxx:3072-3085`: the roots whose first image is still unmodified.
        let mut lof: Vec<TopoShape> = Vec::new();
        for root in self.my_init_offset_face.roots() {
            let Some(first) = self.my_init_offset_face.image(root).into_iter().next() else {
                continue;
            };
            let Some(f) = Face::wrap(first) else {
                continue;
            };
            if !self.my_image_offset.has_image(&f.0) {
                lof.push(f.0);
            }
        }

        // `cxx:3086-3094`.
        if self.my_join == JoinType::Intersection && self.my_inter && self.my_is_planar {
            build_splits_of_trimmed_faces(&lof, &self.my_as_des, &mut self.my_image_offset);
        } else {
            self.my_make_loops
                .build_faces(&lof, &mut self.my_as_des, &mut self.my_image_offset);
        }
    }
}
