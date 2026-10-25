//! `ProjLib_Torus` — project circles onto a `gp_Torus`.
//!
//! Source: `ModelingData/TKGeomBase/ProjLib/ProjLib_Torus.hxx:27-83`,
//! `.cxx:28-193` (OCCT 8.0.0).
//!
//! A circle on a torus is either a V-isoline (its plane is parallel to the
//! torus axis and passes through the centre, `cxx:97-144`) or a U-isoline
//! (`cxx:146-171`); both give a `gp_Lin2d` in the `(U, V)` parameter space.
//! `EvalPnt2d` (`cxx:59-77`) only produces the `U` of a direction; the `V` is
//! filled in by the caller.

use occt_core::gp::{GpCirc, GpDir2d, GpElips, GpHypr, GpLin, GpLin2d, GpParab, GpPnt2d, GpTorus, GpVec, GpVec2d};
use occt_core::kernel::geomabs::CurveType;
use occt_core::precision;

use super::projector::{Projector, ProjectorState};

/// `gp::DY2d()` (`gp.hxx`).
const DY2D: GpDir2d = GpDir2d { x: 0.0, y: 1.0 };

/// `EvalPnt2d(Ve, To)` (`cxx:59-77`): the `U` of the direction `ve` in the
/// torus frame, with `V = 0` (the caller overwrites `V`).
fn eval_pnt2d(ve: &GpVec, to: &GpTorus) -> GpPnt2d {
    let pos = to.position();
    let x = ve.xyz().dot(pos.x_direction().xyz());
    let y = ve.xyz().dot(pos.y_direction().xyz());
    let u = if x.abs() > precision::PCONFUSION || y.abs() > precision::PCONFUSION {
        y.atan2(x)
    } else {
        0.0
    };
    GpPnt2d::new(u, 0.0)
}

/// `class ProjLib_Torus` (`hxx:27-83`).
pub struct Torus {
    base: ProjectorState,
    /// `myTorus` (`hxx:80`).
    torus: GpTorus,
}

impl Torus {
    /// `ProjLib_Torus()` (`cxx:30`).
    pub fn new() -> Self {
        Self {
            base: ProjectorState::new(),
            torus: GpTorus::default(),
        }
    }

    /// `ProjLib_Torus(const gp_Torus&)` (`cxx:36-39`).
    pub fn from_torus(to: &GpTorus) -> Self {
        let mut p = Self::new();
        p.init(to);
        p
    }

    /// `ProjLib_Torus(To, C)` (`cxx:43-47`).
    pub fn with_circ(to: &GpTorus, c: &GpCirc) -> Self {
        let mut p = Self::new();
        p.init(to);
        p.project_circ(c);
        p
    }

    /// `Init(To)` (`cxx:51-57`).
    pub fn init(&mut self, to: &GpTorus) {
        self.base.curve_type = CurveType::OtherCurve;
        self.torus = *to;
        self.base.periodic = false;
        self.base.done = false;
    }

    /// `myTorus` (`hxx:80`).
    pub fn torus(&self) -> &GpTorus {
        &self.torus
    }
}

impl Default for Torus {
    fn default() -> Self {
        Self::new()
    }
}

impl Projector for Torus {
    fn state(&self) -> &ProjectorState {
        &self.base
    }

    fn state_mut(&mut self) -> &mut ProjectorState {
        &mut self.base
    }

    /// `Project(const gp_Circ&)` (`cxx:82-174`).
    fn project_circ(&mut self, c: &GpCirc) {
        self.base.curve_type = CurveType::Line;

        let circ_pos = c.position();
        let torus_pos = self.torus.position();
        let xc = GpVec::from_xyz(circ_pos.x_direction().xyz());
        let yc = GpVec::from_xyz(circ_pos.y_direction().xyz());
        let xt = GpVec::from_xyz(torus_pos.x_direction().xyz());
        let yt = GpVec::from_xyz(torus_pos.y_direction().xyz());
        let zt = GpVec::from_xyz(torus_pos.direction().xyz());
        let oc = GpVec::from_pnts(&self.torus.location(), &c.location());

        let p1;
        let d2;

        if oc.magnitude() < precision::CONFUSION
            || circ_pos
                .direction()
                .is_parallel_tol(&torus_pos.direction(), precision::ANGULAR)
        {
            // Iso V (`cxx:99-144`).
            let mut p1 = eval_pnt2d(&xc, &self.torus); // evaluate U1
            let mut p2 = eval_pnt2d(&yc, &self.torus); // evaluate U2

            let mut z = oc.xyz().dot(torus_pos.direction().xyz());
            z /= self.torus.minor_radius();

            let mut v = if z > 1.0 {
                std::f64::consts::FRAC_PI_2 // simple protection against
                                            // calculation errors
            } else if z < -1.0 {
                -std::f64::consts::FRAC_PI_2
            } else {
                z.asin()
            };

            if c.radius() < self.torus.major_radius() {
                v = std::f64::consts::PI - v;
            } else if v < 0.0 {
                v += 2.0 * std::f64::consts::PI;
            }
            p1.set_y(v);
            p2.set_y(v);

            let mut v2d = GpVec2d::new(p2.x() - p1.x(), p2.y() - p1.y());
            // "Normally |P1.X() - P2.X()| = PI/2. If != PI/2, we crossed the
            // period => reverse the direction" (`cxx:131-135`).
            if (p1.x() - p2.x()).abs() > std::f64::consts::PI {
                v2d.reverse();
            }

            d2 = GpDir2d::new(v2d.x(), v2d.y()).expect("Standard_ConstructionError: gp_Dir2d");
            if p1.x() < 0.0 {
                p1.set_x(2.0 * std::f64::consts::PI + p1.x());
            }
            self.base.lin = GpLin2d::from_pnt_dir(p1, d2);
        } else {
            // Iso U -> U = angle(Xt, OC) (`cxx:146-171`).
            let xt_cross_yt = xt.xyz().crossed(yt.xyz());
            let mut u = xt.angle_with_ref(&oc, &GpVec::from_xyz(&xt_cross_yt));
            if u < 0.0 {
                u += 2.0 * std::f64::consts::PI;
            }

            let oc_cross_zt = oc.xyz().crossed(zt.xyz());
            let mut v1 = oc.angle_with_ref(&xc, &GpVec::from_xyz(&oc_cross_zt));
            if v1 < 0.0 {
                v1 += 2.0 * std::f64::consts::PI;
            }

            // Origin of the line.
            p1 = GpPnt2d::new(u, v1);

            // Direction of the line.
            let mut dir = DY2D;
            let xc_cross_yc = xc.xyz().crossed(yc.xyz());
            if oc_cross_zt.dot(&xc_cross_yc) < 0.0 {
                dir.reverse();
            }

            self.base.lin = GpLin2d::from_pnt_dir(p1, dir);
        }

        self.base.done = true;
    }

    /// `Project(const gp_Lin&)` (`cxx:176-179`): the base body.
    fn project_lin(&mut self, l: &GpLin) {
        self.base.curve_type = CurveType::OtherCurve;
        let _ = l;
    }

    /// `Project(const gp_Elips&)` (`cxx:181-184`): the base body.
    fn project_elips(&mut self, e: &GpElips) {
        self.base.curve_type = CurveType::OtherCurve;
        let _ = e;
    }

    /// `Project(const gp_Parab&)` (`cxx:186-189`): the base body.
    fn project_parab(&mut self, p: &GpParab) {
        self.base.curve_type = CurveType::OtherCurve;
        let _ = p;
    }

    /// `Project(const gp_Hypr&)` (`cxx:191-194`): the base body.
    fn project_hypr(&mut self, h: &GpHypr) {
        self.base.curve_type = CurveType::OtherCurve;
        let _ = h;
    }
}
