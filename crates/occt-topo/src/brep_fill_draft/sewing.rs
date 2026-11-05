//! `BRepFill_Draft::Sewing` (`BRepFill_Draft.cxx:828-916`): assemble the skin
//! with the input shape and hand back the solid.

use occt_core::precision::Precision;

use crate::abs::ShapeType;
use crate::brep_builder_api_sewing::Sewing;
use crate::brep_class3d::SolidClassifier;
use crate::builder::TopoBuilder;
use crate::fclass2d::FaceState;
use crate::shape::Shell;

use super::BrepFillDraft;

impl BrepFillDraft {
    /// `Sewing()` (`cxx:828-916`).
    ///
    /// `Ok(false)` is the reference's `return false` (no assembly performed or
    /// no contiguous edge found); the `Err` arm is the
    /// `TopoDS::Shell(myShape)` cast of the solid branch.
    pub(super) fn sewing(&mut self) -> Result<bool, String> {
        // `cxx:832-838`.
        let to_ass = self.my_top.shape_type() != ShapeType::Wire;
        if !to_ass || !self.my_done {
            return Ok(false);
        }

        // `cxx:841-845`: `new BRepBuilderAPI_Sewing(5 * myTol, true, true,
        // false)` — tolerance, then sewing / analysis on, cutting off.
        let mut ass = Sewing::new(5.0 * self.my_tol, true, true, false, false);
        if let Some(shape) = self.my_shape.clone() {
            ass.add(&shape);
        }
        ass.add(&self.my_top);

        // `cxx:848-850`.
        ass.perform()?;
        let nb_ce = ass.nb_contigous_edges();

        // `cxx:852-871`.
        let mut ok = false;
        if nb_ce > 0 {
            if let Some(res) = ass.sewed_shape().cloned() {
                match res.shape_type() {
                    ShapeType::Shell | ShapeType::Solid => {
                        self.my_shape = Some(res);
                        ok = true;
                    }
                    ShapeType::Compound => {
                        let children = crate::bop_occt_util::iter_children(&res);
                        // "Only one part => this is correct".
                        if children.len() == 1 {
                            self.my_shape = Some(children[0].clone());
                            ok = true;
                        }
                    }
                    _ => {}
                }
            }
        }

        // `cxx:874-909`.
        if ok {
            // "Update the History".
            // `Ass->IsModified` is read first and the stored handle cloned, so
            // the table is not written while it is being read.
            let nb_law = self.nb_law();
            for ii in 1..=nb_law {
                if let Some(old) = self.my_faces.value(1, ii).cloned() {
                    if ass.is_modified(&old) {
                        let new = ass.modified(&old).clone();
                        self.my_faces.set_value(1, ii, Some(new));
                    }
                }
            }
            for ii in 1..=nb_law + 1 {
                if let Some(old) = self.my_sections.value(1, ii).cloned() {
                    if ass.is_modified(&old) {
                        let new = ass.modified(&old).clone();
                        self.my_sections.set_value(1, ii, Some(new));
                    }
                }
            }

            // `cxx:895-909`: "Make a Solid".
            let closed = self.my_shape.as_ref().map(|s| s.is_closed()).unwrap_or(false);
            if closed {
                let b = TopoBuilder::new();
                let shell = self
                    .my_shape
                    .clone()
                    .and_then(Shell::wrap)
                    .ok_or("Standard_TypeMismatch: TopoDS::Shell")?;
                let mut solid = b.make_solid(&[shell]);

                let mut sc = SolidClassifier::new();
                sc.load(solid.0.clone());
                sc.perform_infinite_point(Precision::CONFUSION);
                if sc.state() == FaceState::In {
                    // `cxx:904-907`: rebuild the solid around the reversed
                    // shell.
                    if let Some(shape) = self.my_shape.as_mut() {
                        shape.reverse();
                    }
                    let shell = self
                        .my_shape
                        .clone()
                        .and_then(Shell::wrap)
                        .ok_or("Standard_TypeMismatch: TopoDS::Shell")?;
                    solid = b.make_solid(&[shell]);
                }
                self.my_shape = Some(solid.0);
            }
        }

        Ok(ok)
    }
}
