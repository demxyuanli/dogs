//! Mode fields of `ShapeFix_Shape` / `ShapeFix_Solid` / `ShapeFix_Shell` and
//! the `ShapeFix_Root::NeedFix` gates they feed.
//!
//! Sources:
//! * `ShapeFix_Shape.hxx:143-149` (fields), `ShapeFix_Shape.cxx:40-48`
//!   (ctor defaults), `ShapeFix_Shape.cxx:83-302` (`Perform`).
//! * `ShapeFix_Solid.hxx:109-111` (fields), `ShapeFix_Solid.cxx:52-70` (ctor
//!   defaults), `ShapeFix_Solid.cxx:460-643` (`Perform`).
//! * `ShapeFix_Shell.hxx:125-128` (fields), `ShapeFix_Shell.cxx:68-88` (ctor
//!   defaults), `ShapeFix_Shell.cxx:101-173` (`Perform`).
//! * `ShapeProcess_OperLibrary.cxx:810-822`: the `fixshape` operator copies the
//!   keys into the tools.
//! * `STEPControl_Controller.cxx:208-218`, `:258-259` and `IGESData.cxx:113-123`,
//!   `:163-164`: the per-reader key values.
//!
//! Every key both readers register is `-1` (`FixMode::FixOrNot`) except
//! `FixVertexPositionMode`, which both register as `0` (`FixMode::NotFix`), so
//! `NeedFix` turns every stage on except `ShapeFix::FixVertexPosition`. The one
//! value that differs between the readers is `CreateOpenSolidMode`: `0` on STEP
//! (`STEPControl_Controller.cxx:214`) and `1` on IGES (`IGESData.cxx:119`).
//! `ShapeFix_Shell`'s `myNonManifold` has no resource key of its own: it is
//! `ShapeProcess_ShapeContext::IsNonManifold` (`XSAlgo_ShapeProcessor.cxx:99-103`),
//! set only from an explicit `NonManifold` parameter, so it is false on both
//! reader paths.
//!
//! Call sites of the gates:
//! * [`ShapeFixShapeModes::need_fix_same_parameter`] / `need_fix_vertex_tol`:
//!   `step/read_topology/resolver_context.rs::close_heal_context` and the tail
//!   of `iges_read/iges_process_shape.rs::process_iges_shape`.
//! * [`ShapeFixShapeModes::need_fix_solid`] / `need_fix_shell` /
//!   `need_fix_face` / `need_fix_free_wire`: the shape-type switch of
//!   `iges_process_shape.rs::process_node` (the STEP reader dispatches by STEP
//!   entity type before any heal pass exists, so its equivalent is the entity
//!   resolver in `step/read_topology`).
//! * [`ShapeFixShellModes::need_fix_face`] / `need_fix_orientation`:
//!   `iges_process_shape.rs::perform_iges_shell` and
//!   `step/read_topology/translate_shell_solid.rs::resolve_shell`.
//! * [`ShapeFixSolidModes::need_fix_shell`] / `need_fix_shell_orientation` /
//!   `create_open_solid_mode`: `iges_process_shape.rs::perform_iges_solid` /
//!   `solid_from_one_shell` and
//!   `step/read_topology/translate_shell_solid.rs::resolve_solid_shells`.
//!
//! UNPORTED, and inert on both reader paths:
//! * `ShapeFix::FixVertexPosition` (`ShapeFix.cxx:290-350`): its gate is
//!   `NeedFix(myFixVertexPositionMode)` and both readers register the key as
//!   `0` (`STEPControl_Controller.cxx:258`, `IGESData.cxx:163`), i.e.
//!   `FixMode::NotFix`.
//! * `ShapeFix_Shell`'s `NonManifold` arm: only `FixFaceOrientation`'s
//!   multi-connex branch reads it (`ShapeFix_Shell.cxx:1487-1496`), and the
//!   port's `orient_shell_faces` already declines that case
//!   (`shhealing/shape_fix_shell.rs`).
//! * `ShapeFix_Solid`'s multi-shell arm (`CreateSolids`, `cxx:577-640`) on the
//!   STEP reader: the port only reaches `CreateSolids` from the IGES path
//!   (`iges_read/create_solids.rs`); no gate is missing there, the arm itself
//!   is not wired into `resolve_solid_shells`.

use crate::shhealing::wire_fix::{need_fix, FIX_OR_NOT, NOT_FIX};

/// `ShapeFix_Shape`'s mode fields (`ShapeFix_Shape.hxx:143-149`).
#[derive(Debug, Clone, PartialEq)]
pub struct ShapeFixShapeModes {
    /// `myFixSolidMode`: the `TopAbs_SOLID` arm (`ShapeFix_Shape.cxx:162`).
    pub fix_solid_mode: i32,
    /// `myFixShellMode`: the `TopAbs_SHELL` arm (`cxx:183`).
    pub fix_shell_mode: i32,
    /// `myFixFaceMode`: the `TopAbs_FACE` arm (`cxx:196`).
    pub fix_face_mode: i32,
    /// `myFixWireMode`: the `TopAbs_WIRE` arm (`cxx:213`).
    pub fix_wire_mode: i32,
    /// `myFixSameParameterMode`: the `SameParameter` tail (`cxx:259`), cleared
    /// for the recursing `COMPOUND` / `COMPSOLID` arm (`cxx:136`).
    pub fix_same_parameter_mode: i32,
    /// `myFixVertexPositionMode`: `ShapeFix::FixVertexPosition` before the type
    /// switch (`cxx:118`). Both readers register `0`, so it is off.
    pub fix_vertex_position_mode: i32,
    /// `myFixVertexTolMode`: the multi-face `FixVertexTolerance` tail
    /// (`cxx:267`), saved/restored around the `COMPOUND` arm (`cxx:87`,
    /// `:156`).
    pub fix_vertex_tol_mode: i32,
}

impl Default for ShapeFixShapeModes {
    /// `ShapeFix_Shape::ShapeFix_Shape()` (`ShapeFix_Shape.cxx:40-48`).
    fn default() -> Self {
        ShapeFixShapeModes {
            fix_solid_mode: FIX_OR_NOT,
            fix_shell_mode: FIX_OR_NOT,
            fix_face_mode: FIX_OR_NOT,
            fix_wire_mode: FIX_OR_NOT,
            fix_same_parameter_mode: FIX_OR_NOT,
            fix_vertex_position_mode: NOT_FIX,
            fix_vertex_tol_mode: FIX_OR_NOT,
        }
    }
}

impl ShapeFixShapeModes {
    /// `FromSTEP.FixShape`'s values (`STEPControl_Controller.cxx:208-218`,
    /// `:258-259`) as `ShapeProcess_OperLibrary.cxx:810-816` copies them.
    pub fn step_reader() -> Self {
        Self::default()
    }

    /// `FromIGES.FixShape`'s values (`IGESData.cxx:113-123`, `:163-164`). Same
    /// integers as [`Self::step_reader`]; the readers differ only in
    /// `ShapeFix_Solid`'s `CreateOpenSolidMode`.
    pub fn iges_reader() -> Self {
        Self::default()
    }

    /// `NeedFix(myFixSolidMode)` (`cxx:162`).
    pub fn need_fix_solid(&self) -> bool {
        need_fix(self.fix_solid_mode, true)
    }
    /// `NeedFix(myFixShellMode)` (`cxx:183`).
    pub fn need_fix_shell(&self) -> bool {
        need_fix(self.fix_shell_mode, true)
    }
    /// `NeedFix(myFixFaceMode)` (`cxx:196`).
    pub fn need_fix_face(&self) -> bool {
        need_fix(self.fix_face_mode, true)
    }
    /// `NeedFix(myFixWireMode)` (`cxx:213`).
    pub fn need_fix_free_wire(&self) -> bool {
        need_fix(self.fix_wire_mode, true)
    }
    /// `NeedFix(myFixSameParameterMode)` (`cxx:259`).
    pub fn need_fix_same_parameter(&self) -> bool {
        need_fix(self.fix_same_parameter_mode, true)
    }
    /// `NeedFix(myFixVertexPositionMode)` (`cxx:118`).
    pub fn need_fix_vertex_position(&self) -> bool {
        need_fix(self.fix_vertex_position_mode, true)
    }
    /// `NeedFix(myFixVertexTolMode)` (`cxx:267`).
    pub fn need_fix_vertex_tol(&self) -> bool {
        need_fix(self.fix_vertex_tol_mode, true)
    }

    /// The `COMPOUND` / `COMPSOLID` arm's overrides (`cxx:135-137`, restored at
    /// `cxx:155-156`): `FixSameParameter` and `FixVertexTol` are off for the
    /// recursive per-child `Perform` calls, so the outermost call is the only one
    /// that runs either tail.
    pub fn sub_shape_modes(&self) -> Self {
        let mut m = self.clone();
        m.fix_same_parameter_mode = NOT_FIX;
        m.fix_vertex_tol_mode = NOT_FIX;
        m
    }
}

/// `ShapeFix_Solid`'s mode fields (`ShapeFix_Solid.hxx:109-111`).
#[derive(Debug, Clone, PartialEq)]
pub struct ShapeFixSolidModes {
    /// `myFixShellMode`: the per-shell `ShapeFix_Shell::Perform` loop
    /// (`ShapeFix_Solid.cxx:483`).
    pub fix_shell_mode: i32,
    /// `myFixShellOrientationMode`: the `SolidFromShell` / `CreateSolids` tail
    /// (`cxx:514`).
    pub fix_shell_orientation_mode: i32,
    /// `myCreateOpenSolidMode`: `SolidFromShell` is entered for a shell with
    /// free edges too (`cxx:545`, `:587`, `:617`).
    pub create_open_solid_mode: bool,
}

impl Default for ShapeFixSolidModes {
    /// `ShapeFix_Solid::ShapeFix_Solid()` (`ShapeFix_Solid.cxx:52-70`):
    /// `myCreateOpenSolidMode = false`.
    fn default() -> Self {
        ShapeFixSolidModes {
            fix_shell_mode: FIX_OR_NOT,
            fix_shell_orientation_mode: FIX_OR_NOT,
            create_open_solid_mode: false,
        }
    }
}

impl ShapeFixSolidModes {
    /// `FromSTEP.FixShape` (`STEPControl_Controller.cxx:213-215`): the two modes
    /// are `-1` and `CreateOpenSolidMode` is `0`
    /// (`ShapeProcess_OperLibrary.cxx:820-822`).
    pub fn step_reader() -> Self {
        ShapeFixSolidModes {
            create_open_solid_mode: false,
            ..Self::default()
        }
    }

    /// `FromIGES.FixShape` (`IGESData.cxx:118-120`): same modes,
    /// `CreateOpenSolidMode` is `1`.
    pub fn iges_reader() -> Self {
        ShapeFixSolidModes {
            create_open_solid_mode: true,
            ..Self::default()
        }
    }

    /// `NeedFix(myFixShellMode)` (`ShapeFix_Solid.cxx:483`).
    pub fn need_fix_shell(&self) -> bool {
        need_fix(self.fix_shell_mode, true)
    }
    /// `NeedFix(myFixShellOrientationMode)` (`cxx:514`).
    pub fn need_fix_shell_orientation(&self) -> bool {
        need_fix(self.fix_shell_orientation_mode, true)
    }
}

/// `ShapeFix_Shell`'s mode fields (`ShapeFix_Shell.hxx:125-128`).
#[derive(Debug, Clone, PartialEq)]
pub struct ShapeFixShellModes {
    /// `myFixFaceMode`: the per-face `ShapeFix_Face::Perform` loop
    /// (`ShapeFix_Shell.cxx:110`).
    pub fix_face_mode: i32,
    /// `myFixOrientationMode`: `FixFaceOrientation` (`cxx:140`).
    pub fix_orientation_mode: i32,
    /// `myNonManifold`: passed to `FixFaceOrientation` (`cxx:142`).
    pub non_manifold: bool,
}

impl Default for ShapeFixShellModes {
    /// `ShapeFix_Shell::ShapeFix_Shell()` (`ShapeFix_Shell.cxx:68-88`).
    fn default() -> Self {
        ShapeFixShellModes {
            fix_face_mode: FIX_OR_NOT,
            fix_orientation_mode: FIX_OR_NOT,
            non_manifold: false,
        }
    }
}

impl ShapeFixShellModes {
    /// The reader value set (`ShapeProcess_OperLibrary.cxx:818-820`): no key
    /// changes either mode away from `-1`, and `SetNonManifold` is only called
    /// for an explicit `NonManifold` parameter
    /// (`XSAlgo_ShapeProcessor.cxx:99-103`).
    pub fn reader() -> Self {
        Self::default()
    }

    /// `NeedFix(myFixFaceMode)` (`ShapeFix_Shell.cxx:110`).
    pub fn need_fix_face(&self) -> bool {
        need_fix(self.fix_face_mode, true)
    }
    /// `NeedFix(myFixOrientationMode)` (`cxx:140`).
    pub fn need_fix_orientation(&self) -> bool {
        need_fix(self.fix_orientation_mode, true)
    }
}
