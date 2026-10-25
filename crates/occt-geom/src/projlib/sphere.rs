//! `ProjLib_Sphere` — project elementary curves onto a `gp_Sphere`.
//!
//! Source: `ModelingData/TKGeomBase/ProjLib/ProjLib_Sphere.hxx:31-88`,
//! `.cxx:33-249` (OCCT 8.0.0).
//!
//! Only circles are projected (`cxx:98-179`): a circle whose plane is normal
//! to the sphere axis and concentric with the sphere is a V-isoline
//! (`isIsoU`), one whose plane contains the axis is a U-isoline (`isIsoV`).
//! Both become a `gp_Lin2d` in the `(U, V)` parameter space, so the whole
//! adjustment of the seam (`SetInBounds`, `cxx:204-249`) is about moving that
//! line by whole periods.
//!
//! `EvalPnt2d` (`cxx:66-92`) is the `(U, V)` of a *direction* vector in the
//! sphere frame, not of a point: the callers pass the circle's `XDirection` /
//! `YDirection` so that the pole (`|V| == PI/2`) can be detected in `U`.

use occt_core::elib::clib;
use occt_core::gp::{
    GpAx2d, GpCirc, GpDir, GpDir2d, GpElips, GpHypr, GpLin, GpLin2d, GpParab, GpPnt2d, GpSphere,
    GpTrsf2d, GpVec, GpVec2d,
};
use occt_core::kernel::geomabs::CurveType;
use occt_core::precision;

use super::projector::dir2d_from_pnts;
use super::projector::{Projector, ProjectorState};

/// `gp::DX2d()` (`gp.hxx`) — the 2D X axis direction.
const DX2D: GpDir2d = GpDir2d { x: 1.0, y: 0.0 };
/// `gp::DY2d()` (`gp.hxx`).
const DY2D: GpDir2d = GpDir2d { x: 0.0, y: 1.0 };

/// `EvalPnt2d(P, Sp)` (`cxx:66-92`): `U = atan2(Y, X)` in `[0, 2*PI)` with the
/// pole folded onto `U = 0`, `V = asin(Z)` clamped to `[-1, 1]`.
fn eval_pnt2d(p: &GpDir, sp: &GpSphere) -> GpPnt2d {
    let pos = sp.position();
    let x = p.xyz().dot(pos.x_direction().xyz());
    let y = p.xyz().dot(pos.y_direction().xyz());
    let mut z = p.xyz().dot(pos.direction().xyz());

    let u = if x.abs() > precision::PCONFUSION || y.abs() > precision::PCONFUSION {
        clib::in_period(y.atan2(x), 0.0, 2.0 * std::f64::consts::PI)
    } else {
        0.0
    };

    if z > 1.0 {
        z = 1.0;
    } else if z < -1.0 {
        z = -1.0;
    }
    let v = z.asin();

    GpPnt2d::new(u, v)
}

/// `class ProjLib_Sphere` (`hxx:31-88`).
pub struct Sphere {
    base: ProjectorState,
    /// `mySphere` (`hxx:85`).
    sphere: GpSphere,
}

impl Sphere {
    /// `ProjLib_Sphere()` (`cxx:38`).
    pub fn new() -> Self {
        Self {
            base: ProjectorState::new(),
            sphere: GpSphere::default(),
        }
    }

    /// `ProjLib_Sphere(const gp_Sphere&)` (`cxx:42-45`).
    pub fn from_sphere(sp: &GpSphere) -> Self {
        let mut p = Self::new();
        p.init(sp);
        p
    }

    /// `ProjLib_Sphere(Sp, C)` (`cxx:49-53`).
    pub fn with_circ(sp: &GpSphere, c: &GpCirc) -> Self {
        let mut p = Self::new();
        p.init(sp);
        p.project_circ(c);
        p
    }

    /// `Init(Sp)` (`cxx:57-63`).
    pub fn init(&mut self, sp: &GpSphere) {
        self.base.curve_type = CurveType::OtherCurve;
        self.sphere = *sp;
        self.base.periodic = false;
        self.base.done = false;
    }

    /// `mySphere` (`hxx:85`).
    pub fn sphere(&self) -> &GpSphere {
        &self.sphere
    }

    /// `SetInBounds(U)` (`cxx:206-249`): bring the line's value at `U` into the
    /// `V` range `[-PI/2, PI/2]`, mirroring at the pole when needed and then
    /// wrapping `U` into `[0, 2*PI)`.
    ///
    /// Panics where OCCT raises `StdFail_NotDone` (`cxx:207`).
    pub fn set_in_bounds(&mut self, u: f64) {
        assert!(
            self.base.done,
            "StdFail_NotDone: ProjLib_Sphere::SetInBounds"
        );

        // "first set the y of the first point in -pi/2 pi/2" (`cxx:208-212`).
        let y = clib::line2d_value(&self.base.lin, u).y();
        let new_y = clib::in_period(y, -std::f64::consts::PI, std::f64::consts::PI);
        self.base.lin.translate_vec(&GpVec2d::new(0.0, new_y - y));

        let p = clib::line2d_value(&self.base.lin, u);
        const TOL: f64 = 1.0e-7;
        let d2d = *self.base.lin.direction();

        let axis_sign = if p.y() - std::f64::consts::FRAC_PI_2 > TOL
            || ((p.y() - std::f64::consts::FRAC_PI_2).abs() < TOL && d2d.is_equal(&DY2D, TOL))
        {
            1.0
        } else if p.y() + std::f64::consts::FRAC_PI_2 < -TOL
            || ((p.y() + std::f64::consts::FRAC_PI_2).abs() < TOL && d2d.is_opposite(&DY2D, TOL))
        {
            -1.0
        } else {
            return;
        };

        let mut trsf = GpTrsf2d::identity();
        let axis = GpAx2d::new(GpPnt2d::new(0.0, axis_sign * std::f64::consts::FRAC_PI_2), DX2D);
        trsf.set_mirror_ax2d(&axis);
        self.base.lin.transform(&trsf);

        self.base
            .lin
            .translate_vec(&GpVec2d::new(std::f64::consts::PI, 0.0));

        // "Now we need to adjust the U parameter" (`cxx:244-249`).
        let x = clib::line2d_value(&self.base.lin, u).x();
        let new_x = clib::in_period(x, 0.0, 2.0 * std::f64::consts::PI);
        self.base.lin.translate_vec(&GpVec2d::new(new_x - x, 0.0));
    }
}

impl Default for Sphere {
    fn default() -> Self {
        Self::new()
    }
}

impl Projector for Sphere {
    fn state(&self) -> &ProjectorState {
        &self.base
    }

    fn state_mut(&mut self) -> &mut ProjectorState {
        &mut self.base
    }

    /// `Project(const gp_Circ&)` (`cxx:98-179`).
    // `P2d1` / `P2d2` / `D2d` are declared up front in C++ (`cxx:121-122`) and
    // `myLin` is written from them unconditionally at `cxx:179`; the
    // "neither branch" case therefore keeps the default-constructed values,
    // which is why the initializers below are never read.
    #[allow(unused_assignments)]
    fn project_circ(&mut self, c: &GpCirc) {
        // "Check the validity: Xc & Yc must be perpendicular to Zs -> IsoV;
        // O, Zs is in the Plane O, Xc, Yc -> IsoU" (`cxx:103-106`).
        let o = self.sphere.location();
        let circ_pos = c.position();
        let xc = circ_pos.x_direction();
        let yc = circ_pos.y_direction();
        let zc = xc.crossed(yc).expect("Standard_ConstructionError: gp_Dir::Crossed");
        let sphere_pos = self.sphere.position();
        let xs = sphere_pos.x_direction();
        let ys = sphere_pos.y_direction();
        let zs = sphere_pos.direction();

        const TOL: f64 = precision::CONFUSION;

        let is_iso_u =
            zc.is_normal_tol(&zs, TOL) && o.is_equal_with_tolerance(&c.location(), TOL);
        let is_iso_v = xc.is_normal_tol(&zs, TOL) && yc.is_normal_tol(&zs, TOL);

        let mut p2d1 = GpPnt2d::new(0.0, 0.0);
        let mut p2d2 = GpPnt2d::new(0.0, 0.0);
        let mut d2d = GpDir2d::default();

        if is_iso_u {
            self.base.curve_type = CurveType::Line;

            p2d1 = eval_pnt2d(&xc, &self.sphere);
            p2d2 = eval_pnt2d(&yc, &self.sphere);

            if (p2d1.y() - std::f64::consts::FRAC_PI_2).abs() < precision::PCONFUSION
                || (p2d1.y() + std::f64::consts::FRAC_PI_2).abs() < precision::PCONFUSION
            {
                // "then P1 is on the apex of the sphere and U is undefined.
                // The value of U is given by P2d2.Y()." (`cxx:134-138`).
                p2d1.set_x(p2d2.x());
            } else if ((p2d1.x() - p2d2.x()).abs() - std::f64::consts::PI).abs()
                < precision::PCONFUSION
            {
                // "then we have U2 = U1 + PI; V2; we have to assume that
                // U1 = U2 so V2 = PI - V2" (`cxx:139-153`).
                p2d2.set_x(p2d1.x());
                if p2d2.y() < 0.0 {
                    p2d2.set_y(-std::f64::consts::PI - p2d2.y());
                } else {
                    p2d2.set_y(std::f64::consts::PI - p2d2.y());
                }
            } else {
                p2d2.set_x(p2d1.x());
            }

            let Some(dir) = dir2d_from_pnts(&p2d1, &p2d2) else {
                panic!("Standard_ConstructionError: gp_Dir2d");
            };
            d2d = dir;
            self.base.done = true;
        } else if is_iso_v {
            self.base.curve_type = CurveType::Line;

            // "P2d(U,V): first point of the PCurve." (`cxx:166-174`).
            let xs_cross_ys = xs.crossed(&ys).expect("Standard_ConstructionError: gp_Dir::Crossed");
            let mut u = xs.angle_with_ref(&xc, &xs_cross_ys);
            if u < 0.0 {
                u += 2.0 * std::f64::consts::PI;
            }
            let z = GpVec::from_pnts(&o, &c.location())
                .xyz()
                .dot(zs.xyz());
            let v = (z / self.sphere.radius()).asin();
            p2d1 = GpPnt2d::new(u, v);
            let dot = zc.dot(&xs_cross_ys);
            d2d = GpDir2d::new(dot, 0.0).expect("Standard_ConstructionError: gp_Dir2d");
            self.base.done = true;
        }
        self.base.lin = GpLin2d::from_pnt_dir(p2d1, d2d);
    }

    /// `Project(const gp_Lin&)` (`cxx:182-185`): the base body.
    fn project_lin(&mut self, l: &GpLin) {
        self.base.curve_type = CurveType::OtherCurve;
        let _ = l;
    }

    /// `Project(const gp_Elips&)` (`cxx:187-190`): the base body.
    fn project_elips(&mut self, e: &GpElips) {
        self.base.curve_type = CurveType::OtherCurve;
        let _ = e;
    }

    /// `Project(const gp_Parab&)` (`cxx:192-195`): the base body.
    fn project_parab(&mut self, p: &GpParab) {
        self.base.curve_type = CurveType::OtherCurve;
        let _ = p;
    }

    /// `Project(const gp_Hypr&)` (`cxx:197-200`): the base body.
    fn project_hypr(&mut self, h: &GpHypr) {
        self.base.curve_type = CurveType::OtherCurve;
        let _ = h;
    }
}
