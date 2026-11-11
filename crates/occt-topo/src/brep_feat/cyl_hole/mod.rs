//! `BRepFeat_MakeCylindricalHole` — cylindrical holes along an axis.
//!
//! Source: `ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_MakeCylindricalHole.hxx:1-116`,
//! `.cxx:1-810`, `.lxx:1-54`.
//!
//! The class derives from `BRepFeat_Builder` and drives it exactly the way the
//! reference does: it builds an elongated cylinder as the tool, runs the
//! inherited `BOPAlgo_BOP::Perform`, then keeps the parts of the tool whose
//! parameter along the axis falls in the requested range. `Build()` only runs
//! `PerformResult` and `Validate` (`cxx:630-648`).
//!
//! Module table:
//!
//! | Reference | Rust |
//! |---|---|
//! | ctor / `Init` / `Status` (`lxx:21-53`) | [`MakeCylindricalHole::new`] / `init` / `init_axis` / `status` |
//! | `Perform(Radius)` (`cxx:56-100`) | `perform` |
//! | `PerformThruNext` (`cxx:106-253`) | `perform_thru_next` |
//! | `PerformUntilEnd` (`cxx:257-353`) | `perform_until_end` |
//! | `Perform(Radius, PFrom, PTo, Cont)` (`cxx:357-470`) | `perform_range` |
//! | `PerformBlind` (`cxx:474-626`) | `perform_blind` |
//! | `Build` (`cxx:630-648`) | [`MakeCylindricalHole::build`] |
//! | `Validate` (`cxx:652-688`) | `validate` |
//! | `Baryc` / `BoxParameters` / `GetOffset` / `CreateCyl` | `statics` |
//!
//! Entry points: [`crate::brepfeat::through_hole`] and
//! [`crate::brepfeat::blind_hole`] (`Perform(Radius)` and `PerformBlind` plus
//! `Build`).
//!
//! Deviations (every one of them a reference dereference the port answers with
//! a value instead):
//! - `myArguments.First()` on an empty list and a null object raise
//!   `Standard_ConstructionError` (`cxx:58-62`); the port answers `Err`.
//! - `BoxParameters` raises on an empty `Bnd_Box` (`B.Get`, `cxx:725`); the
//!   port answers `None` and the caller reports `InvalidPlacement`.
//! - `BRepPrim_Cylinder(a2, R, H)` builds its faces in the given `gp_Ax2`
//!   (`BRepPrim_Cylinder.cxx:34-64`); the port's `primitives::BRepPrimCylinder`
//!   only takes the identity frame, so `statics::cylinder` transforms a copy of
//!   it and bakes the frame into the geometry, which is what the reference
//!   constructor answers (identity `TopLoc_Location`).
//! - `GetOffset` dereferences the `BRepAdaptor_Surface` handle of a face that
//!   may be null (`cxx:757`); the port answers `false` (`None`), i.e. falls
//!   back to `offF = Radius` (`cxx:792-795`).
//! - `PerformBlind` fills `theList` with the consumed faces (`cxx:533-537`) and
//!   never reads it; the port omits the dead fill.

mod perform;
mod perform_blind;
mod statics;

pub use statics::{baryc, box_parameters, create_cyl, cylinder, get_offset};

use occt_core::gp::GpAx1;

use crate::abs::ShapeType;
use crate::brep_feat::builder::BRepFeatBuilder;
use crate::brep_feat::status::Status;
use crate::shape::{Face, TopoShape};
use crate::topo_tools::explorer_shapes;
/// `BRepFeat_MakeCylindricalHole` (`hxx:37-110`).
pub struct MakeCylindricalHole {
    /// The inherited `BRepFeat_Builder`.
    base: BRepFeatBuilder,
    /// `myAxis` (`hxx:105`).
    my_axis: GpAx1,
    /// `myAxDef` (`hxx:106`).
    my_ax_def: bool,
    /// `myStatus` (`hxx:107`).
    my_status: Status,
    /// `myIsBlind` (`hxx:108`).
    my_is_blind: bool,
    /// `myValidate` (`hxx:109`).
    my_validate: bool,
    /// `myTopFace` (`hxx:110`).
    my_top_face: Option<Face>,
    /// `myBotFace` (`hxx:111`).
    my_bot_face: Option<Face>,
}

impl Default for MakeCylindricalHole {
    fn default() -> Self {
        Self::new()
    }
}

impl MakeCylindricalHole {
    /// `BRepFeat_MakeCylindricalHole()` (`lxx:21-27`): `myAxDef` false,
    /// `myStatus` `NoError`, `myIsBlind` and `myValidate` false.
    pub fn new() -> Self {
        Self {
            base: BRepFeatBuilder::new(),
            my_axis: GpAx1::default(),
            my_ax_def: false,
            my_status: Status::NoError,
            my_is_blind: false,
            my_validate: false,
            my_top_face: None,
            my_bot_face: None,
        }
    }

    /// `Init(Axis)` (`lxx:32-36`).
    pub fn init_axis(&mut self, axis: &GpAx1) {
        self.my_axis = *axis;
        self.my_ax_def = true;
    }

    /// `Init(S, Axis)` (`lxx:40-45`): `BRepFeat_Builder::Init(S)` then the
    /// axis.
    pub fn init(&mut self, s: &TopoShape, axis: &GpAx1) {
        self.base.init(s);
        self.my_axis = *axis;
        self.my_ax_def = true;
    }

    /// `Status()` (`lxx:49-52`).
    pub fn status(&self) -> Status {
        self.my_status
    }

    /// `Shape()` of the inherited `BOPAlgo_BuilderShape`.
    pub fn shape(&self) -> &TopoShape {
        self.base.shape()
    }

    /// `Modified(S)` of the inherited `BOPAlgo_BuilderShape`.
    pub fn modified(&self, s: &TopoShape) -> Vec<TopoShape> {
        self.base.modified(s)
    }

    /// `IsDeleted(S)` of the inherited `BOPAlgo_BuilderShape`.
    pub fn is_deleted(&self, s: &TopoShape) -> bool {
        self.base.is_deleted(s)
    }

    /// `HasErrors()` of the inherited `BOPAlgo_Builder`.
    pub fn has_errors(&self) -> bool {
        self.base.has_errors()
    }

    /// `Build()` (`cxx:630-648`): `PerformResult`, then `Validate` when the
    /// caller asked for it, then `myShape = Shape()`.
    ///
    /// The reference's `PerformResult` may raise; the port turns that into
    /// `InvalidPlacement`, the same status the reference writes when
    /// `HasErrors()` answers true.
    pub fn build(&mut self) {
        if self.my_status != Status::NoError {
            return;
        }
        if self.base.perform_result().is_err() {
            self.my_status = Status::InvalidPlacement;
            return;
        }
        if self.base.has_errors() {
            self.my_status = Status::InvalidPlacement;
            return;
        }
        self.my_status = if self.my_validate {
            self.validate()
        } else {
            Status::NoError
        };
        // `cxx:641`: `myShape = Shape()` — the builder already holds that
        // shape, so the assignment has nothing to do.
    }

    /// `Validate()` (`cxx:652-688`): a blind hole must keep its top face in the
    /// result, a through hole must not keep either cap.
    fn validate(&self) -> Status {
        let faces = explorer_shapes(self.base.shape(), ShapeType::Face);
        if self.my_is_blind {
            for f in &faces {
                if self.top_face_is(f) {
                    return Status::NoError;
                }
            }
            return Status::HoleTooLong;
        }
        for f in &faces {
            if self.top_face_is(f) {
                return Status::InvalidPlacement;
            }
        }
        for f in &faces {
            if self.bot_face_is(f) {
                return Status::InvalidPlacement;
            }
        }
        Status::NoError
    }

    /// `ex.Current().IsSame(myTopFace)` (`cxx:660`, `:675`): the same `TShape`
    /// **and** the same `TopLoc_Location` (`TopoDS_Shape.hxx:266-272`); a null
    /// `myTopFace` matches nothing.
    fn top_face_is(&self, f: &TopoShape) -> bool {
        self.my_top_face
            .as_ref()
            .is_some_and(|t| crate::topo_tools::is_same_occurrence(f, &t.0))
    }

    /// `ex.Current().IsSame(myBotFace)` (`cxx:681`), same comparison as
    /// [`Self::top_face_is`].
    fn bot_face_is(&self, f: &TopoShape) -> bool {
        self.my_bot_face
            .as_ref()
            .is_some_and(|b| crate::topo_tools::is_same_occurrence(f, &b.0))
    }

    /// `myArguments.First()` plus the two `IsNull` guards of every entry point
    /// (`cxx:58-62`, `:108-112`, `:258-262`, `:361-365`, `:478-482`).
    fn checked_object(&self) -> Result<TopoShape, String> {
        let object = self
            .base
            .first_argument()
            .cloned()
            .ok_or_else(|| "Standard_ConstructionError".to_string())?;
        if object.is_null() || !self.my_ax_def {
            return Err("Standard_ConstructionError".to_string());
        }
        Ok(object)
    }

    /// `AddTool(theTool); SetOperation(Fuse); BOPAlgo_BOP::Perform();`
    /// (`cxx:97-100`, `:178-181`, `:318-321`, `:435-438`, `:567-570`).
    fn add_tool_and_perform(&mut self, tool: &TopoShape) -> Result<(), String> {
        self.base.add_tool(tool);
        // `bool Fuse = false;`
        self.base.set_operation(0);
        self.base.perform()
    }
}
