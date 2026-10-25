//! `ProjLib_Cylinder` — project elementary curves onto a `gp_Cylinder`.
//!
//! Source: `ModelingData/TKGeomBase/ProjLib/ProjLib_Cylinder.hxx:27-84`,
//! `.cxx:28-175` (OCCT 8.0.0).
//!
//! Only curves that lie on the cylinder in the sense the file checks get an
//! analytic 2D form: the `Project` bodies return early — leaving `myType` at
//! `GeomAbs_OtherCurve` and `isDone == false` — when the curve is not parallel
//! to the axis (`cxx:99-104`, `cxx:126-131`), and `Project(gp_Elips)`
//! deliberately projects nothing so that the caller approximates
//! (`cxx:159-164`).

use occt_core::gp::{GpCirc, GpCylinder, GpDir2d, GpElips, GpHypr, GpLin, GpLin2d, GpParab, GpPnt2d};
use occt_core::kernel::geomabs::CurveType;
use occt_core::precision;

use super::projector::{Projector, ProjectorState};

/// `EvalPnt2d(P, Cyl)` (`cxx:72-90`): the parameters of `P` on the cylinder,
/// `U = atan2(Y, X)` measured from the cylinder's own X axis.
fn eval_pnt2d(p: &occt_core::gp::GpPnt, cy: &GpCylinder) -> GpPnt2d {
    let pos = cy.position();
    let op = occt_core::gp::GpVec::from_pnts(&cy.location(), p);
    let x = op.xyz().dot(pos.x_direction().xyz());
    let y = op.xyz().dot(pos.y_direction().xyz());
    let z = op.xyz().dot(pos.direction().xyz());
    let u = if x.abs() > precision::PCONFUSION || y.abs() > precision::PCONFUSION {
        y.atan2(x)
    } else {
        0.0
    };
    GpPnt2d::new(u, z)
}

/// `class ProjLib_Cylinder` (`hxx:27-84`).
pub struct Cylinder {
    base: ProjectorState,
    /// `myCylinder` (`hxx:81`).
    cylinder: GpCylinder,
}

impl Cylinder {
    /// `ProjLib_Cylinder()` (`cxx:28`).
    pub fn new() -> Self {
        Self {
            base: ProjectorState::new(),
            cylinder: GpCylinder::default(),
        }
    }

    /// `ProjLib_Cylinder(const gp_Cylinder&)` (`cxx:34-37`).
    pub fn from_cylinder(cy: &GpCylinder) -> Self {
        let mut p = Self::new();
        p.init(cy);
        p
    }

    /// `ProjLib_Cylinder(Cyl, L)` (`cxx:40-44`).
    pub fn with_line(cy: &GpCylinder, l: &GpLin) -> Self {
        let mut p = Self::new();
        p.init(cy);
        p.project_lin(l);
        p
    }

    /// `ProjLib_Cylinder(Cyl, C)` (`cxx:48-52`).
    pub fn with_circ(cy: &GpCylinder, c: &GpCirc) -> Self {
        let mut p = Self::new();
        p.init(cy);
        p.project_circ(c);
        p
    }

    /// `ProjLib_Cylinder(Cyl, E)` (`cxx:56-60`).
    pub fn with_elips(cy: &GpCylinder, e: &GpElips) -> Self {
        let mut p = Self::new();
        p.init(cy);
        p.project_elips(e);
        p
    }

    /// `Init(Cyl)` (`cxx:64-70`).
    pub fn init(&mut self, cy: &GpCylinder) {
        self.base.curve_type = CurveType::OtherCurve;
        self.cylinder = cy.clone();
        self.base.periodic = false;
        self.base.done = false;
    }

    /// `myCylinder` (`hxx:81`).
    pub fn cylinder(&self) -> &GpCylinder {
        &self.cylinder
    }
}

impl Default for Cylinder {
    fn default() -> Self {
        Self::new()
    }
}

impl Projector for Cylinder {
    fn state(&self) -> &ProjectorState {
        &self.base
    }

    fn state_mut(&mut self) -> &mut ProjectorState {
        &mut self.base
    }

    /// `Project(const gp_Lin&)` (`cxx:96-118`).
    fn project_lin(&mut self, l: &GpLin) {
        // "Check the line is parallel to the axis of cylinder. In other cases,
        // the projection is wrong." (`cxx:98-104`).
        let axis = self.cylinder.position().direction();
        if l.direction()
            .xyz()
            .crossed(axis.xyz())
            .square_modulus()
            > precision::ANGULAR * precision::ANGULAR
        {
            return;
        }

        self.base.curve_type = CurveType::Line;

        let mut p2d = eval_pnt2d(&l.location(), &self.cylinder);
        if p2d.x() < 0.0 {
            p2d.set_x(p2d.x() + 2.0 * std::f64::consts::PI);
        }
        let mut signe = l.direction().dot(&axis);
        signe = if signe > 0.0 { 1.0 } else { -1.0 };
        let d2d = GpDir2d::new(0.0, signe).expect("Standard_ConstructionError: gp_Dir2d");

        self.base.lin = GpLin2d::from_pnt_dir(p2d, d2d);
        self.base.done = true;
    }

    /// `Project(const gp_Circ&)` (`cxx:122-157`).
    fn project_circ(&mut self, c: &GpCirc) {
        // "Check the circle's normal is parallel to the axis of cylinder. In
        // other cases, the projection is wrong." (`cxx:124-131`).
        let cyl_pos = self.cylinder.position();
        let circ_pos = c.position();
        if circ_pos
            .direction()
            .xyz()
            .crossed(cyl_pos.direction().xyz())
            .square_modulus()
            > precision::ANGULAR * precision::ANGULAR
        {
            return;
        }

        self.base.curve_type = CurveType::Line;

        let z_cyl = cyl_pos
            .x_direction()
            .crossed(cyl_pos.y_direction())
            .expect("Standard_ConstructionError: gp_Dir::Crossed");

        let u = cyl_pos
            .x_direction()
            .angle_with_ref(circ_pos.x_direction(), &z_cyl);

        let op = occt_core::gp::GpVec::from_pnts(&self.cylinder.location(), &c.location());
        let v = op.xyz().dot(cyl_pos.direction().xyz());

        let p2d1 = GpPnt2d::new(u, v);
        let d2d = if z_cyl.dot(&circ_pos.direction()) > 0.0 {
            GpDir2d::new(1.0, 0.0).expect("Standard_ConstructionError: gp_Dir2d")
        } else {
            GpDir2d::new(-1.0, 0.0).expect("Standard_ConstructionError: gp_Dir2d")
        };

        self.base.lin = GpLin2d::from_pnt_dir(p2d1, d2d);
        self.base.done = true;
    }

    /// `Project(const gp_Elips&)` (`cxx:159-164`): "the projection of an
    /// ellipse onto a cylinder is delegated to approximation", i.e. the body
    /// is empty and the projector stays not done.
    fn project_elips(&mut self, _e: &GpElips) {}

    /// `Project(const gp_Parab&)` (`cxx:167-171`): the base body.
    fn project_parab(&mut self, p: &GpParab) {
        self.base.curve_type = CurveType::OtherCurve;
        let _ = p;
    }

    /// `Project(const gp_Hypr&)` (`cxx:173-176`): the base body.
    fn project_hypr(&mut self, h: &GpHypr) {
        self.base.curve_type = CurveType::OtherCurve;
        let _ = h;
    }
}
