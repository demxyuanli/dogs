//! `ProjLib_Cone` — project elementary curves onto a `gp_Cone`.
//!
//! Source: `ModelingData/TKGeomBase/ProjLib/ProjLib_Cone.hxx:27-85`,
//! `.cxx:30-180` (OCCT 8.0.0).
//!
//! `Project(gp_Lin)` (`cxx:69-104`) only answers when the line is parallel to
//! the cone's V-isoline through its own first parameter — the case of a
//! ruling; the coordinates come straight from `ElSLib::ConeParameters` /
//! `ElSLib::ConeD1`, with the apex shifted onto the line first so the
//! parameters are defined. `Project(gp_Circ)` (`cxx:108-165`) carries the
//! "U V without the Trsf round trip" code of `ElSLib` and bails out (`isDone =
//! false`) as soon as the circle's normal leaves the cone axis.

use occt_core::elib::{slib, surface_eval};
use occt_core::gp::{GpCirc, GpCone, GpDir, GpDir2d, GpElips, GpHypr, GpLin, GpLin2d, GpPnt2d, GpVec};
use occt_core::kernel::geomabs::CurveType;
use occt_core::precision;

use super::projector::dir_from_vec;
use super::projector::{Projector, ProjectorState};

/// `class ProjLib_Cone` (`hxx:27-85`).
pub struct Cone {
    base: ProjectorState,
    /// `myCone` (`hxx:82`).
    cone: GpCone,
}

impl Cone {
    /// `ProjLib_Cone()` (`cxx:32`).
    pub fn new() -> Self {
        Self {
            base: ProjectorState::new(),
            cone: GpCone::default(),
        }
    }

    /// `ProjLib_Cone(const gp_Cone&)` (`cxx:36-39`).
    pub fn from_cone(co: &GpCone) -> Self {
        let mut p = Self::new();
        p.init(co);
        p
    }

    /// `ProjLib_Cone(Co, L)` (`cxx:43-47`).
    pub fn with_line(co: &GpCone, l: &GpLin) -> Self {
        let mut p = Self::new();
        p.init(co);
        p.project_lin(l);
        p
    }

    /// `ProjLib_Cone(Co, C)` (`cxx:51-55`).
    pub fn with_circ(co: &GpCone, c: &GpCirc) -> Self {
        let mut p = Self::new();
        p.init(co);
        p.project_circ(c);
        p
    }

    /// `Init(Co)` (`cxx:59-65`).
    pub fn init(&mut self, co: &GpCone) {
        self.base.curve_type = CurveType::OtherCurve;
        self.cone = co.clone();
        self.base.periodic = false;
        self.base.done = false;
    }

    /// `myCone` (`hxx:82`).
    pub fn cone(&self) -> &GpCone {
        &self.cone
    }
}

impl Default for Cone {
    fn default() -> Self {
        Self::new()
    }
}

/// `gp_Dir2d(0., Signe)` with `Signe = +-1`.
fn dir2d_v(sign: f64) -> GpDir2d {
    GpDir2d::new(0.0, sign).expect("Standard_ConstructionError: gp_Dir2d")
}

/// `gp_Dir2d(+-1., 0.)` chosen from the sign of `z_cone . z_circ`.
fn dir2d_u(z_cone: &GpDir, z_circ: &GpDir) -> GpDir2d {
    if z_cone.dot(z_circ) > 0.0 {
        GpDir2d::new(1.0, 0.0).expect("Standard_ConstructionError: gp_Dir2d")
    } else {
        GpDir2d::new(-1.0, 0.0).expect("Standard_ConstructionError: gp_Dir2d")
    }
}

impl Projector for Cone {
    fn state(&self) -> &ProjectorState {
        &self.base
    }

    fn state_mut(&mut self) -> &mut ProjectorState {
        &mut self.base
    }

    /// `Project(const gp_Lin&)` (`cxx:69-104`).
    fn project_lin(&mut self, l: &GpLin) {
        let apex = self.cone.apex();
        let mut a_pnt = l.location();

        let mut a_delta_v = 0.0;

        // "Take another point in the line L, which does not coincide with the
        // cone apex." (`cxx:77-81`).
        if a_pnt.is_equal_with_tolerance(&apex, precision::CONFUSION) {
            a_pnt.translate_vec(&GpVec::from_xyz(l.direction().xyz()));
            a_delta_v = 1.0; // == ||L.Direction()|| == 1.0
        }

        let pos = self.cone.position();
        let (u, v) = slib::cone_parameters(&pos, self.cone.radius(), self.cone.semi_angle(), &a_pnt);

        let (_p, _vu, vv) = surface_eval::cone_d1(&self.cone, u, v);

        let Some(dv) = dir_from_vec(&vv) else {
            return;
        };
        if dv.is_parallel_tol(&l.direction(), precision::ANGULAR) {
            // "L is parallel to U-isoline of the cone."
            self.base.curve_type = CurveType::Line;

            let a_sign = if l.direction().dot(&dv) >= 0.0 { 1.0 } else { -1.0 };
            let p2d = GpPnt2d::new(u, v - a_delta_v * a_sign);
            let d2d = dir2d_v(a_sign);

            self.base.lin = GpLin2d::from_pnt_dir(p2d, d2d);

            self.base.done = true;
        }
    }

    /// `Project(const gp_Circ&)` (`cxx:108-165`).
    fn project_circ(&mut self, c: &GpCirc) {
        self.base.curve_type = CurveType::Line;

        let cone_pos = self.cone.position();
        let circ_pos = c.position();
        if !cone_pos
            .direction()
            .is_parallel_tol(&circ_pos.direction(), precision::ANGULAR)
        {
            self.base.done = false;
            return;
        }

        let z_cone = cone_pos
            .x_direction()
            .crossed(cone_pos.y_direction())
            .expect("Standard_ConstructionError: gp_Dir::Crossed");
        let z_cir = circ_pos
            .x_direction()
            .crossed(circ_pos.y_direction())
            .expect("Standard_ConstructionError: gp_Dir::Crossed");

        // "to find point U V, we use the code from ElSLib without applying the
        // Trsf to the point (unnecessary round trip)." (`cxx:129-131`).
        let x = cone_pos.x_direction().dot(circ_pos.x_direction());
        let y = cone_pos.y_direction().dot(circ_pos.x_direction());
        let z = GpVec::from_pnts(&self.cone.location(), &c.location())
            .xyz()
            .dot(cone_pos.direction().xyz());

        let mut u = if x.abs() <= precision::ANGULAR && y.abs() <= precision::ANGULAR {
            0.0
        } else if -self.cone.radius() > z * self.cone.semi_angle().tan() {
            (-y).atan2(-x)
        } else {
            y.atan2(x)
        };
        if u < 0.0 {
            u += 2.0 * std::f64::consts::PI;
        }

        let v = z / self.cone.semi_angle().cos();

        let p2d1 = GpPnt2d::new(u, v);
        let d2d = dir2d_u(&z_cone, &z_cir);

        self.base.lin = GpLin2d::from_pnt_dir(p2d1, d2d);
        self.base.done = true;
    }

    /// `Project(const gp_Elips&)` (`cxx:167-170`): the base body.
    fn project_elips(&mut self, e: &GpElips) {
        self.base.curve_type = CurveType::OtherCurve;
        let _ = e;
    }

    /// `Project(const gp_Parab&)` (`cxx:172-175`): the base body.
    fn project_parab(&mut self, p: &occt_core::gp::GpParab) {
        self.base.curve_type = CurveType::OtherCurve;
        let _ = p;
    }

    /// `Project(const gp_Hypr&)` (`cxx:177-180`): the base body.
    fn project_hypr(&mut self, h: &GpHypr) {
        self.base.curve_type = CurveType::OtherCurve;
        let _ = h;
    }
}
