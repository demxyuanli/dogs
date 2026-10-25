//! `GeomFill_SectionPlacement` — places a section on a sweep spine.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_SectionPlacement.cxx:157-983`
//! (OCCT 8.0.0). The class looks for the parameter on the path whose frame
//! best matches the section: it first tries the section's mean plane against
//! the path (plane/curve intersection), then curve/curve extrema, then a
//! point/path fallback, and scores every candidate with `Choix` (distance
//! gain against `Gabarit`, then a penalty on the angle).
//!
//! Translation notes:
//! - `occ::handle<Geom_Geometry>` is the sum type [`SectionGeometry`]: OCCT
//!   branches on `IsInstance(Geom_CartesianPoint)` in the constructor
//!   (`cxx:169-180`); every other `Geom_Geometry` is `occ::down_cast` to
//!   `Geom_Curve` (`cxx:177`), i.e. a curve here.
//! - `GeomAdaptor_Curve myAdpSection` is the `Arc<dyn Curve>` itself: the
//!   adaptor only forwards `GetType` / `FirstParameter` / `Value` / `NbPoles`
//!   / `IsClosed` / `LocateU` to that curve, all of which
//!   [`crate::curve::Curve`] answers (see `curve.rs` and `trimmed.rs` for the
//!   trimmed-curve forwards that make `GetType` match
//!   `GeomAdaptor_Curve::load`).
//! - `Extrema_ExtPC myExt` (initialized once in the constructor,
//!   `cxx:373-376`) is the free function
//!   [`crate::extrema_pc::extrema_ext_pc_range_done`] over the section's range:
//!   `Perform` is re-entrant, the "initialize" carries no other state.
//! - `IntCurveSurface_HInter` against the plane surface is
//!   [`crate::geom_fill::plane_curve_inter::plane_intersections`], the same
//!   entry `GeomFill_GuideTrihedronPlan` uses.
//! - `BndLib_Add3dCurve::Add` is `crate::geom_bnd_lib::box_curve`.
//! - `Extrema_ExtCC` is
//!   [`crate::extrema_cc::curve_curve_extrema_all_range`]. `IsParallel` has no
//!   port-side counterpart (the parallel arms report no solution), so the loop
//!   sees an empty sequence exactly when OCCT would have skipped it.
//! - `SectionAxis` (`cxx:908-949`) is **UNPORTED**: it is private, never called
//!   anywhere in OCCT, and needs `GeomLProp_CLProps` (3D), which this crate
//!   does not have yet (only `GeomLProp_CLProps2d` exists,
//!   `occt-geom2d::cl_props2d`).
//! - OCCT throws `Standard_Failure` from `Transformation` for a non-planar
//!   profile (`cxx:830`) and from `gp_Dir` on a null direction; those become
//!   `Err(String)`.

use std::sync::Arc;

use occt_core::gp::dir::DirAxis;
use occt_core::gp::{GpAx1, GpAx3, GpDir, GpMat, GpPnt, GpTrsf, GpVec};
use occt_core::precision::ANGULAR;

use crate::convert_geom::curve_as_handle;
use crate::curve::Curve;
use crate::geom_fill::location_law::LocationLaw;

mod init_section;
mod math_helpers;
mod search;

/// `Section->IsInstance(STANDARD_TYPE(Geom_CartesianPoint))` versus
/// `occ::down_cast<Geom_Curve>(Section)` (`GeomFill_SectionPlacement.cxx:169-180`).
pub enum SectionGeometry {
    Point(GpPnt),
    Curve(Arc<dyn Curve>),
}

/// `GeomFill_SectionPlacement` (`GeomFill_SectionPlacement.hxx:31-100`).
pub struct SectionPlacement {
    /// `done` (`hxx:82`).
    done: bool,
    /// `isplan` (`hxx:83`): the section has a mean plane, `TheAxe` is its axis.
    is_plan: bool,
    /// `TheAxe` (`hxx:84`).
    the_axe: GpAx1,
    /// `Gabarit` (`hxx:85`): half the section box diagonal.
    gabarit: f64,
    /// `myLaw` (`hxx:86`).
    law: Arc<dyn LocationLaw>,
    /// `myAdpSection` / `mySection` (`hxx:87-88`): `None` when the section is a
    /// point.
    section_curve: Option<Arc<dyn Curve>>,
    /// `SecParam` (`hxx:89`).
    sec_param: f64,
    /// `PathParam` (`hxx:90`).
    path_param: f64,
    /// `Dist` (`hxx:91`).
    dist: f64,
    /// `AngleMax` (`hxx:92`).
    angle_max: f64,
    /// `myIsPoint` (`hxx:94`).
    is_point: bool,
    /// `myPoint` (`hxx:95`).
    point: GpPnt,
}

impl SectionPlacement {
    /// `GeomFill_SectionPlacement(L, Section)` (`cxx:157-380`).
    pub fn new(law: Arc<dyn LocationLaw>, section: SectionGeometry) -> Self {
        let mut placement = SectionPlacement {
            done: false,
            is_plan: false,
            // `gp_Ax1` default constructor: origin, `gp::DZ()` (`gp_Ax1.hxx:53-56`).
            the_axe: GpAx1::new(
                GpPnt::new(0.0, 0.0, 0.0),
                GpDir::from_axis(DirAxis::Z),
            ),
            gabarit: 0.0,
            law,
            section_curve: None,
            sec_param: 0.0,
            path_param: 0.0,
            // `Dist(RealLast())` (`cxx:161`).
            dist: f64::MAX,
            angle_max: 0.0,
            is_point: false,
            point: GpPnt::new(0.0, 0.0, 0.0),
        };
        match section {
            SectionGeometry::Point(p) => {
                placement.is_point = true;
                placement.point = p;
                placement.is_plan = true;
            }
            SectionGeometry::Curve(c) => {
                placement.section_curve = Some(c);
            }
        }
        placement.init_section();
        placement
    }

    /// `GeomFill_SectionPlacement::SetLocation` (`cxx:383-387`).
    pub fn set_location(&mut self, law: Arc<dyn LocationLaw>) {
        self.law = law;
    }

    /// `GeomFill_SectionPlacement::Perform(Tol)` (`cxx:390-396`).
    pub fn perform(&mut self, tol: f64) {
        let path = self
            .law
            .get_curve()
            .expect("GeomFill_SectionPlacement::Perform: null LocationLaw curve");
        self.perform_path(path.as_ref(), tol);
    }

    /// `GeomFill_SectionPlacement::Perform(Path, Tol)` (`cxx:400-704`).
    pub fn perform_on(&mut self, path: &dyn Curve, tol: f64) {
        self.perform_path(path, tol);
    }

    /// `GeomFill_SectionPlacement::Perform(ParamOnPath, Tol)` (`cxx:712-753`).
    pub fn perform_at(&mut self, param: f64, tol: f64) {
        self.done = true;
        let path = self
            .law
            .get_curve()
            .expect("GeomFill_SectionPlacement::Perform: null LocationLaw curve");

        self.path_param = param;
        if self.is_point {
            let pon_path = path.d0(self.path_param);
            self.dist = pon_path.distance(&self.point);
            self.angle_max = std::f64::consts::FRAC_PI_2;
        } else {
            let section = self
                .section_curve
                .clone()
                .expect("GeomFill_SectionPlacement: null section");
            self.sec_param = section.first_parameter();

            let v_ref = GpVec::from_xyz(self.the_axe.direction().xyz());
            let (pon_path, dp1) = math_helpers::tangente(path.as_ref(), self.path_param);
            let pon_sec = section.d0(self.sec_param);
            self.dist = pon_path.distance(&pon_sec);
            if self.dist > tol {
                // On Cherche un meilleur point sur la section
                let (sols, ext_done) = init_section::extrema_on_section(
                    section.as_ref(),
                    &pon_path,
                );
                if ext_done {
                    search::dist_mini(
                        section.as_ref(),
                        &pon_path,
                        &sols,
                        ext_done,
                        &mut self.dist,
                        &mut self.sec_param,
                    );
                }
            }
            self.angle_max = math_helpers::eval_angle(&v_ref, &dp1);
            if self.is_plan {
                self.angle_max = std::f64::consts::FRAC_PI_2 - self.angle_max;
            }
        }

        self.done = true;
    }

    /// `GeomFill_SectionPlacement::IsDone` (`cxx:757-760`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `GeomFill_SectionPlacement::ParameterOnPath` (`cxx:764-767`).
    pub fn parameter_on_path(&self) -> f64 {
        self.path_param
    }

    /// `GeomFill_SectionPlacement::ParameterOnSection` (`cxx:772-775`).
    pub fn parameter_on_section(&self) -> f64 {
        self.sec_param
    }

    /// `GeomFill_SectionPlacement::Distance` (`cxx:780-783`).
    pub fn distance(&self) -> f64 {
        self.dist
    }

    /// `GeomFill_SectionPlacement::Angle` (`cxx:788-791`).
    pub fn angle(&self) -> f64 {
        self.angle_max
    }

    /// `GeomFill_SectionPlacement::Transformation(WithTranslation,
    /// WithCorrection)` (`cxx:795-887`).
    ///
    /// `Err` is the OCCT `Standard_Failure("Illegal usage: can't rotate
    /// non-planar profile")` (`cxx:830`) and the `gp_Dir` null-magnitude
    /// throws.
    pub fn transformation(
        &self,
        with_translation: bool,
        with_correction: bool,
    ) -> Result<GpTrsf, String> {
        let (m, v) = self.law.d0(self.path_param)?;

        // Calcul des reperes
        let mut p = GpPnt::from_xyz(v.xyz());
        let d = dir_from_column(&m, 3)?;
        let dn = dir_from_column(&m, 1)?;
        let paxe = GpAx3::new(p, d, &dn)?;

        let mut p_section = GpPnt::new(0.0, 0.0, 0.0);
        if with_translation || with_correction {
            p_section = if self.is_point {
                self.point
            } else {
                self.section_curve
                    .as_ref()
                    .expect("GeomFill_SectionPlacement: null section")
                    .value(self.sec_param)
            };
        }

        let mut rot = GpTrsf::identity();
        if with_correction && !self.is_point {
            if !self.is_plan {
                return Err("Illegal usage: can't rotate non-planar profile".into());
            }

            let profile_normal = self.the_axe.direction();
            let spine_start_dir = paxe.direction();

            // `!ProfileNormal.IsParallel(SpineStartDir, Precision::Angular())`
            if !profile_normal.is_parallel_tol(&spine_start_dir, ANGULAR) {
                let dir_axe_of_rotation = profile_normal.crossed(&spine_start_dir)?;
                let angle =
                    profile_normal.angle_with_ref(&spine_start_dir, &dir_axe_of_rotation);
                let axe_of_rotation = GpAx1::new(*self.the_axe.location(), dir_axe_of_rotation);
                rot = GpTrsf::identity();
                rot.set_rotation_ax1(&axe_of_rotation, angle)?;
            }
            p_section.transform(&rot);
        }

        if with_translation {
            // `P.ChangeCoord().SetLinearForm(-1, PSection.XYZ(), V.XYZ())`
            let xyz = GpVec::zero().set_linear_form_2(
                -1.0,
                &GpVec::from_xyz(p_section.xyz()),
                1.0,
                &v,
            );
            p = GpPnt::from_xyz(xyz.xyz());
        } else {
            p = GpPnt::new(0.0, 0.0, 0.0);
        }

        let saxe = GpAx3::new(
            p,
            GpDir::from_axis(DirAxis::Z),
            &GpDir::from_axis(DirAxis::X),
        )?;

        // Transfo
        let mut tf = GpTrsf::identity();
        tf.set_transformation_from_to(&saxe, &paxe);

        if with_correction {
            tf.multiply(&rot);
        }

        Ok(tf)
    }

    /// `GeomFill_SectionPlacement::Section(WithTranslation)` (`cxx:890-895`).
    ///
    /// OCCT dereferences the null `mySection` handle for a point section, so
    /// this returns `None` there.
    pub fn section(&self, with_translation: bool) -> Result<Option<Arc<dyn Curve>>, String> {
        let Some(section) = self.section_curve.as_ref() else {
            return Ok(None);
        };
        let tf = self.transformation(with_translation, false)?;
        let mut copy = section.clone_dyn();
        copy.transform(&tf);
        Ok(Some(curve_as_handle(copy)))
    }

    /// `GeomFill_SectionPlacement::ModifiedSection(WithTranslation)`
    /// (`cxx:899-904`). `None` for a point section, as [`Self::section`].
    pub fn modified_section(
        &self,
        with_translation: bool,
    ) -> Result<Option<Arc<dyn Curve>>, String> {
        let Some(section) = self.section_curve.as_ref() else {
            return Ok(None);
        };
        let tf = self.transformation(with_translation, true)?;
        let mut copy = section.clone_dyn();
        copy.transform(&tf);
        Ok(Some(curve_as_handle(copy)))
    }
}

/// `gp_Dir::SetXYZ(M.Column())` (`cxx:807-808`): a zero column is the OCCT
/// `gp_VectorWithNullMagnitude` throw.
fn dir_from_column(m: &GpMat, col: usize) -> Result<GpDir, String> {
    let xyz = m.column(col);
    GpDir::new(xyz.x, xyz.y, xyz.z).map_err(|e| e.to_string())
}
