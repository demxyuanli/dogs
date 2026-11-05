//! `BRepOffsetAPI_MakeEvolved` — the evolved-shape API wrapper.
//!
//! Source: `ModelingAlgorithms/TKOffset/BRepOffsetAPI/BRepOffsetAPI_MakeEvolved.cxx:29-149`
//! and `.hxx:78-122` (OCCT 8.0.0).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `()` | `cxx:29` | [`BRepOffsetApiMakeEvolved::new`] |
//! | `(Spine, Profil, Join, AxeProf, Solid, ProfOnSpine, Tol, IsVolume, RunInParallel)` | `cxx:31-90` | [`BRepOffsetApiMakeEvolved::with_options`] |
//! | `Evolved()` | `cxx:96-105` | [`BRepOffsetApiMakeEvolved::evolved`] |
//! | `Build()` | `cxx:109-121` | [`BRepOffsetApiMakeEvolved::build`] |
//! | `Top()` | `cxx:125-128` | [`BRepOffsetApiMakeEvolved::top`] |
//! | `Bottom()` | `cxx:132-135` | [`BRepOffsetApiMakeEvolved::bottom`] |
//! | `GeneratedShapes(SpineShape, ProfShape)` | `cxx:137-148` | [`BRepOffsetApiMakeEvolved::generated_shapes`] |
//!
//! # Translation notes
//!
//! * The wrapper holds the two performers verbatim (`hxx:119-121`): a
//!   [`BRepFillEvolved`] and an [`AdvancedEvolved`] plus `myIsVolume`. Only the
//!   one selected by the constructor performs, and `Build` copies whichever of
//!   the two answered `IsDone()` (`cxx:111-117`).
//! * A spine that is neither a wire nor a face throws
//!   `Standard_TypeMismatch` (`cxx:44-48`) before either performer is touched;
//!   the port answers `Err` from the constructor.
//! * The volume arm takes the spine wire from the spine itself, or — for a
//!   face — from its first `TopoDS_Iterator` child (`cxx:51-59`), which is the
//!   face's outer wire.
//! * Both early exits (`cxx:61-64`, `cxx:75-78`) leave `myShape` null and never
//!   call `Build`, so the object stays not done; the port returns the `Self`
//!   built so far, exactly the state those returns leave behind.
//! * `Evolved()` throws `Standard_TypeMismatch` in volume mode (`cxx:98-103`);
//!   the port answers `Err`.
//! * `Top()` / `Bottom()` return the performer's null `TopoDS_Shape` before
//!   `Build` ran; that null is `None` here.
//! * `GeneratedShapes` answers the static empty list when `myEvolved` is not
//!   done (`cxx:144-147`); the port answers an empty `Vec`.
//! * `myIsVolume` is left uninitialized by the defaulted default constructor
//!   (`cxx:29`, `hxx:121`); the port initializes it to `false`, the value the
//!   `.hxx` defaulted argument gives the other constructor.

use occt_core::gp::{GpAx3, GpDir, GpPnt};
use occt_core::kernel::geomabs::JoinType;
use occt_core::precision;

use crate::abs::ShapeType;
use crate::brep_fill;
use crate::brep_fill_advanced_evolved::AdvancedEvolved;
use crate::brep_fill_evolved::BRepFillEvolved;
use crate::iterator::ShapeIterator;
use crate::shape::{Face, TopoShape, Wire};

/// `BRepOffsetAPI_MakeEvolved` (`BRepOffsetAPI_MakeEvolved.hxx:78-122`).
pub struct BRepOffsetApiMakeEvolved {
    /// `myEvolved` (`hxx:119`).
    my_evolved: BRepFillEvolved,
    /// `myVolume` (`hxx:120`).
    my_volume: AdvancedEvolved,
    /// `myIsVolume` (`hxx:121`).
    my_is_volume: bool,
    /// `myShape` (`BRepBuilderAPI_MakeShape::myShape`).
    my_shape: Option<TopoShape>,
    /// `myDone` (`BRepBuilderAPI_Command::myDone`), set by `Done()`.
    my_done: bool,
}

impl Default for BRepOffsetApiMakeEvolved {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepOffsetApiMakeEvolved {
    /// `BRepOffsetAPI_MakeEvolved()` (`cxx:29`): every member defaulted.
    pub fn new() -> Self {
        Self {
            my_evolved: BRepFillEvolved::new(),
            my_volume: AdvancedEvolved::new(),
            // `hxx:121` is uninitialized in the reference; see the module
            // note.
            my_is_volume: false,
            my_shape: None,
            my_done: false,
        }
    }

    /// `BRepOffsetAPI_MakeEvolved(theSpine, theProfile, theJoinType,
    /// theIsAxeProf, theIsSolid, theIsProfOnSpine, theTol, theIsVolume,
    /// theRunInParallel)` (`cxx:31-90`).
    ///
    /// The reference defaults are `theJoinType = GeomAbs_Arc`,
    /// `theIsAxeProf = true`, `theIsSolid = false`, `theIsProfOnSpine = false`,
    /// `theTol = 1e-7`, `theIsVolume = false`, `theRunInParallel = false`
    /// (`hxx:92-97`).
    #[allow(clippy::too_many_arguments)]
    pub fn with_options(
        spine: &TopoShape,
        profile: &Wire,
        join: JoinType,
        is_axe_prof: bool,
        is_solid: bool,
        is_prof_on_spine: bool,
        tol: f64,
        is_volume: bool,
        run_in_parallel: bool,
    ) -> Result<Self, String> {
        let mut this = Self::new();
        // `cxx:42`.
        this.my_is_volume = is_volume;

        // `cxx:44-48`.
        let spine_type = spine.shape_type();
        if spine_type != ShapeType::Wire && spine_type != ShapeType::Face {
            return Err(
                "Standard_TypeMismatch: BRepOffsetAPI_MakeEvolved: face or wire is expected as a spine"
                    .to_string(),
            );
        }

        if is_volume {
            // `cxx:50`.
            this.my_volume.set_parallel_mode(run_in_parallel);

            // `cxx:51-59`: `TopoDS::Wire(Spine)` or
            // `TopoDS::Wire(TopoDS_Iterator(Spine).Value())`.
            let a_spine = if spine_type == ShapeType::Wire {
                Wire::wrap(spine.clone())
                    .ok_or_else(|| "Standard_TypeMismatch: TopoDS::Wire".to_string())?
            } else {
                let first = ShapeIterator::of_shape(spine).next().ok_or_else(|| {
                    "Standard_NullObject: TopoDS_Iterator on a spine face".to_string()
                })?;
                Wire::wrap(first)
                    .ok_or_else(|| "Standard_TypeMismatch: TopoDS::Wire".to_string())?
            };

            // `cxx:60-64`.
            this.my_volume.perform(&a_spine, profile, tol, is_solid);
            if !this.my_volume.is_done() {
                return Ok(this);
            }
        } else {
            // `cxx:68`: `gp_Ax3 Axis(gp::Origin(), gp::DZ(), gp::DX())`.
            let mut axis = GpAx3::new(
                GpPnt::new(0.0, 0.0, 0.0),
                GpDir::new(0.0, 0.0, 1.0).expect("gp::DZ"),
                &GpDir::new(1.0, 0.0, 0.0).expect("gp::DX"),
            )
            .expect("gp_Ax3: Z and X are perpendicular");

            // `cxx:70-78`.
            if !is_axe_prof {
                let (a_axis, pos) = brep_fill::axe(spine, profile, tol.max(precision::CONFUSION))?;
                axis = a_axis;
                if is_prof_on_spine && !pos {
                    return Ok(this);
                }
            }

            // `cxx:79-86`.
            if spine_type == ShapeType::Wire {
                let a_spine = Wire::wrap(spine.clone())
                    .ok_or_else(|| "Standard_TypeMismatch: TopoDS::Wire".to_string())?;
                this.my_evolved
                    .perform_wire(&a_spine, profile, &axis, join, is_solid)?;
            } else {
                let a_spine = Face::wrap(spine.clone())
                    .ok_or_else(|| "Standard_TypeMismatch: TopoDS::Face".to_string())?;
                this.my_evolved
                    .perform_face(&a_spine, profile, &axis, join, is_solid)?;
            }
        }

        // `cxx:89`.
        this.build();
        Ok(this)
    }

    /// `Evolved()` (`cxx:96-105`).
    pub fn evolved(&self) -> Result<&BRepFillEvolved, String> {
        // `cxx:98-103`.
        if self.my_is_volume {
            return Err(
                "Standard_TypeMismatch: BRepOffsetAPI_MakeEvolved: myEvolved is accessed while in volume mode"
                    .to_string(),
            );
        }
        Ok(&self.my_evolved)
    }

    /// `Build(theRange)` (`cxx:109-121`): copy the shape of the performer that
    /// answered `IsDone()`, then `Done()`.
    pub fn build(&mut self) {
        // `cxx:111-117`.
        if self.my_evolved.is_done() {
            self.my_shape = self.my_evolved.shape().cloned();
        } else if self.my_volume.is_done() {
            self.my_shape = self.my_volume.shape().cloned();
        }
        // `cxx:119`.
        self.my_done = true;
    }

    /// `Top()` (`cxx:125-128`): the face top of `myEvolved`.
    pub fn top(&self) -> Option<&TopoShape> {
        self.my_evolved.top()
    }

    /// `Bottom()` (`cxx:132-135`): the face bottom of `myEvolved`.
    pub fn bottom(&self) -> Option<&TopoShape> {
        self.my_evolved.bottom()
    }

    /// `GeneratedShapes(SpineShape, ProfShape)` (`cxx:137-148`).
    pub fn generated_shapes(
        &self,
        spine_shape: &TopoShape,
        prof_shape: &TopoShape,
    ) -> Vec<TopoShape> {
        // `cxx:144-147`.
        if !self.my_evolved.is_done() {
            return Vec::new();
        }
        self.my_evolved.generated_shapes(spine_shape, prof_shape)
    }

    /// `IsDone()` — the `BRepBuilderAPI_Command` flag set by `Build`.
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `Shape()` (`BRepBuilderAPI_MakeShape`): the null shape is `None` until
    /// `Build` copied a performer result.
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_shape.as_ref()
    }

    /// `myIsVolume` as set by the constructor (`cxx:42`).
    pub fn is_volume(&self) -> bool {
        self.my_is_volume
    }

    /// The `AdvancedEvolved` used by the volume arm.
    pub fn volume(&self) -> &AdvancedEvolved {
        &self.my_volume
    }
}
