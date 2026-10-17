//! `RWMesh_CoordinateSystemConverter` — length-unit and axis-system conversion
//! for mesh import/export.
//!
//! Source: `RWMesh_CoordinateSystemConverter.hxx` / `.cxx`,
//! `RWMesh_CoordinateSystem.hxx`.
//!
//! The converter is a pure post-process over three element kinds — positions,
//! normals and locations — and it is *empty* unless the application defines both
//! the input and the output side:
//!
//! * unit scaling needs both length units to be positive (`cxx:39-57`);
//! * the axis change needs both `gp_Ax3`s to be set (`cxx:58-66`).
//!
//! `RWGltf_CafReader`'s constructor sets only the input side — unit 1.0 and
//! `RWMesh_CoordinateSystem_glTF` (`RWGltf_CafReader.cxx:182-183`) — so a plain
//! read performs no conversion at all and the file's Y-up metres reach the
//! document untouched. An application that wants Z-up output sets the output
//! system on the converter and the same read applies it.

use occt_core::gp::{GpAx3, GpDir, GpPnt, GpTrsf, GpXyz, TrsfForm};
use occt_core::gp::dir::DirAxis;

/// `RWMesh_CoordinateSystem` (`RWMesh_CoordinateSystem.hxx:21-40`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RwMeshCoordinateSystem {
    /// `RWMesh_CoordinateSystem_Undefined` (-1).
    Undefined,
    /// `RWMesh_CoordinateSystem_posYfwd_posZup` (Blender, `+YForward+Zup+Xright`).
    PosYfwdPosZup,
    /// `RWMesh_CoordinateSystem_negZfwd_posYup` (glTF, `-ZForward+Yup+Xright`).
    NegZfwdPosYup,
}

impl RwMeshCoordinateSystem {
    /// `RWMesh_CoordinateSystem_glTF` / `_Yup`: the name glTF uses for the
    /// `-ZForward+Yup+Xright` convention (`hxx:38-40`).
    pub const GLTF: Self = Self::NegZfwdPosYup;
    /// `RWMesh_CoordinateSystem_Zup` / `_Blender` (`hxx:33-37`).
    pub const ZUP: Self = Self::PosYfwdPosZup;
}

/// `RWMesh_CoordinateSystemConverter::StandardCoordinateSystem`
/// (`RWMesh_CoordinateSystemConverter.hxx:56-69`).
pub fn standard_coordinate_system(sys: RwMeshCoordinateSystem) -> GpAx3 {
    match sys {
        RwMeshCoordinateSystem::PosYfwdPosZup => {
            // `gp_Ax3(gp::Origin(), gp::DZ(), gp::DX())`
            GpAx3::new(
                GpPnt::default(),
                GpDir::from_axis(DirAxis::Z),
                &GpDir::from_axis(DirAxis::X),
            )
            .unwrap_or_else(|_| GpAx3::standard())
        }
        RwMeshCoordinateSystem::NegZfwdPosYup => {
            // `gp_Ax3(gp::Origin(), gp::DY(), gp::DX())`
            GpAx3::new(
                GpPnt::default(),
                GpDir::from_axis(DirAxis::Y),
                &GpDir::from_axis(DirAxis::X),
            )
            .unwrap_or_else(|_| GpAx3::standard())
        }
        RwMeshCoordinateSystem::Undefined => GpAx3::default(),
    }
}

/// `RWMesh_CoordinateSystemConverter` (`RWMesh_CoordinateSystemConverter.hxx:55-203`).
#[derive(Debug, Clone)]
pub struct RwMeshCoordinateSystemConverter {
    input_ax3: GpAx3,
    output_ax3: GpAx3,
    input_length_unit: f64,
    output_length_unit: f64,
    has_input_ax3: bool,
    has_output_ax3: bool,
    trsf: GpTrsf,
    trsf_inv: GpTrsf,
    norm_trsf: [[f64; 4]; 4],
    unit_factor: f64,
    has_scale: bool,
    is_empty: bool,
}

impl Default for RwMeshCoordinateSystemConverter {
    /// `RWMesh_CoordinateSystemConverter::RWMesh_CoordinateSystemConverter()`
    /// (`cxx:21-31`): units -1.0, no systems, identity transform, empty.
    fn default() -> Self {
        Self {
            input_ax3: GpAx3::default(),
            output_ax3: GpAx3::default(),
            input_length_unit: -1.0,
            output_length_unit: -1.0,
            has_input_ax3: false,
            has_output_ax3: false,
            trsf: GpTrsf::identity(),
            trsf_inv: GpTrsf::identity(),
            norm_trsf: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
            unit_factor: 1.0,
            has_scale: false,
            is_empty: true,
        }
    }
}

impl RwMeshCoordinateSystemConverter {
    /// Empty constructor (`cxx:21-31`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `IsEmpty()` (`hxx:77`): TRUE when no transformation is applied.
    pub fn is_empty(&self) -> bool {
        self.is_empty
    }

    /// `InputLengthUnit()` (`hxx:82`).
    pub fn input_length_unit(&self) -> f64 {
        self.input_length_unit
    }

    /// `SetInputLengthUnit()` (`hxx:84-87`).
    pub fn set_input_length_unit(&mut self, unit: f64) {
        self.init(self.input_ax3.clone(), unit, self.output_ax3.clone(), self.output_length_unit);
    }

    /// `OutputLengthUnit()` (`hxx:92`).
    pub fn output_length_unit(&self) -> f64 {
        self.output_length_unit
    }

    /// `SetOutputLengthUnit()` (`hxx:94-97`).
    pub fn set_output_length_unit(&mut self, unit: f64) {
        self.init(self.input_ax3.clone(), self.input_length_unit, self.output_ax3.clone(), unit);
    }

    /// `HasInputCoordinateSystem()` (`hxx:102`).
    pub fn has_input_coordinate_system(&self) -> bool {
        self.has_input_ax3
    }

    /// `InputCoordinateSystem()` (`hxx:105`).
    pub fn input_coordinate_system(&self) -> &GpAx3 {
        &self.input_ax3
    }

    /// `SetInputCoordinateSystem(const gp_Ax3&)` (`hxx:107-112`).
    pub fn set_input_coordinate_system_ax3(&mut self, sys: GpAx3) {
        self.has_input_ax3 = true;
        self.init(sys, self.input_length_unit, self.output_ax3.clone(), self.output_length_unit);
    }

    /// `SetInputCoordinateSystem(RWMesh_CoordinateSystem)` (`hxx:114-119`): an
    /// `Undefined` enum clears the flag but still installs the default `gp_Ax3`.
    pub fn set_input_coordinate_system(&mut self, sys: RwMeshCoordinateSystem) {
        self.has_input_ax3 = sys != RwMeshCoordinateSystem::Undefined;
        self.init(
            standard_coordinate_system(sys),
            self.input_length_unit,
            self.output_ax3.clone(),
            self.output_length_unit,
        );
    }

    /// `HasOutputCoordinateSystem()` (`hxx:122`).
    pub fn has_output_coordinate_system(&self) -> bool {
        self.has_output_ax3
    }

    /// `OutputCoordinateSystem()` (`hxx:125`).
    pub fn output_coordinate_system(&self) -> &GpAx3 {
        &self.output_ax3
    }

    /// `SetOutputCoordinateSystem(const gp_Ax3&)` (`hxx:127-132`).
    pub fn set_output_coordinate_system_ax3(&mut self, sys: GpAx3) {
        self.has_output_ax3 = true;
        self.init(self.input_ax3.clone(), self.input_length_unit, sys, self.output_length_unit);
    }

    /// `SetOutputCoordinateSystem(RWMesh_CoordinateSystem)` (`hxx:134-139`).
    pub fn set_output_coordinate_system(&mut self, sys: RwMeshCoordinateSystem) {
        self.has_output_ax3 = sys != RwMeshCoordinateSystem::Undefined;
        self.init(
            self.input_ax3.clone(),
            self.input_length_unit,
            standard_coordinate_system(sys),
            self.output_length_unit,
        );
    }

    /// `Init(theInputSystem, theInputLengthUnit, theOutputSystem,
    /// theOutputLengthUnit)` (`cxx:37-70`).
    pub fn init(
        &mut self,
        input_system: GpAx3,
        input_length_unit: f64,
        output_system: GpAx3,
        output_length_unit: f64,
    ) {
        self.input_length_unit = input_length_unit;
        self.output_length_unit = output_length_unit;
        self.input_ax3 = input_system;
        self.output_ax3 = output_system;
        if input_length_unit > 0.0 && output_length_unit > 0.0 {
            self.unit_factor = input_length_unit / output_length_unit;
            // `std::abs(myUnitFactor - 1.0) > gp::Resolution()`
            self.has_scale =
                (self.unit_factor - 1.0).abs() > occt_core::precision::RESOLUTION;
        } else {
            self.unit_factor = 1.0;
            self.has_scale = false;
        }

        let mut trsf = GpTrsf::identity();
        if self.has_input_ax3 && self.has_output_ax3 {
            // `aTrsf.SetTransformation(theOutputSystem, theInputSystem)`.
            trsf.set_transformation_from_to(&output_system, &input_system);
            // `cxx:58-63`: an identity change of basis (equal systems) collapses
            // to the default `gp_Trsf`.
            if trsf.translation_part().modulus() <= occt_core::precision::RESOLUTION
                && rotation_is_identity(&trsf)
            {
                trsf = GpTrsf::identity();
            }
        }

        self.trsf = trsf.clone();
        self.trsf_inv = trsf.inverted().unwrap_or_else(|_| GpTrsf::identity());
        self.norm_trsf = trsf.get_mat4();
        // `myIsEmpty = !myHasScale && myTrsf.Form() == gp_Identity` (`cxx:69`).
        self.is_empty = !self.has_scale && trsf.form() == TrsfForm::Identity;
    }

    /// `TransformTransformation(gp_Trsf&)` (`hxx:148-161`).
    pub fn transform_transformation(&self, trsf: &mut GpTrsf) {
        if self.has_scale {
            let t = trsf.translation_part().multiply(self.unit_factor);
            trsf.set_translation_part(&occt_core::gp::GpVec::from_xyz(&t));
        }
        if self.trsf.form() != TrsfForm::Identity {
            // `theTrsf = myTrsf * theTrsf * myTrsfInv`
            *trsf = self.trsf.multiplied(trsf).multiplied(&self.trsf_inv);
        }
    }

    /// `TransformPosition(gp_XYZ&)` (`hxx:163-173`).
    pub fn transform_position(&self, pos: &mut GpXyz) {
        if self.has_scale {
            // `thePosition.Multiply(myUnitFactor)` — `gp_XYZ::Multiply` scales
            // in place; `GpXyz::multiply` returns the scaled value instead.
            *pos = pos.multiply(self.unit_factor);
        }
        if self.trsf.form() != TrsfForm::Identity {
            self.trsf.transforms_xyz(pos);
        }
    }

    /// `TransformNormal(NCollection_Vec3<float>&)` (`hxx:175-184`): the
    /// translation is excluded by appending a 0 as the homogeneous coordinate.
    ///
    /// The C++ multiplies `NCollection_Mat4<float>` by `NCollection_Vec4<float>`,
    /// so the whole sum runs in single precision *including* the
    /// `translation * 0.0f` term. Dropping that term changes the sign of a zero
    /// component (`-0.0f + 0.0f` is `+0.0f`), which the byte-level glTF baseline
    /// shows, so the port keeps the fourth product.
    pub fn transform_normal(&self, norm: &mut GpXyz) {
        if self.trsf.form() == TrsfForm::Identity {
            return;
        }
        let v = [norm.x as f32, norm.y as f32, norm.z as f32];
        let mut out = [0.0f32; 3];
        for (i, o) in out.iter_mut().enumerate() {
            *o = self.norm_trsf[i][0] as f32 * v[0]
                + self.norm_trsf[i][1] as f32 * v[1]
                + self.norm_trsf[i][2] as f32 * v[2]
                + self.norm_trsf[i][3] as f32 * 0.0;
        }
        *norm = GpXyz::new(out[0] as f64, out[1] as f64, out[2] as f64);
    }
}

/// `gp_Quaternion::IsEqual(gp_Quaternion())` — the rotation half of the
/// identity check in `cxx:59-60`.
fn rotation_is_identity(trsf: &GpTrsf) -> bool {
    let q = trsf.get_rotation();
    occt_core::precision::RESOLUTION >= q.x().abs()
        && occt_core::precision::RESOLUTION >= q.y().abs()
        && occt_core::precision::RESOLUTION >= q.z().abs()
        && occt_core::precision::RESOLUTION >= (q.w() - 1.0).abs()
}
