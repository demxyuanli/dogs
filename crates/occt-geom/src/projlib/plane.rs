//! `ProjLib_Plane` — project elementary curves onto a `gp_Pln`'s parameter
//! plane.
//!
//! Source: `ModelingData/TKGeomBase/ProjLib/ProjLib_Plane.hxx:27-90`,
//! `.cxx:24-169` (OCCT 8.0.0).
//!
//! Every overload keeps the analytic type: a circle stays a `gp_Circ2d` with
//! the same radius and the projected axes (`cxx:111-122`), and so on. The two
//! static helpers `EvalPnt2d` / `EvalDir2d` (`cxx:87-97`) are the module-level
//! [`super::eval_pln_pnt2d`] / [`super::eval_pln_dir2d`] already used by the
//! pcurve code.

use occt_core::gp::{
    GpAx22d, GpCirc, GpDir, GpDir2d, GpElips, GpHypr, GpLin, GpLin2d, GpParab, GpPln, GpPnt2d,
};
use occt_core::kernel::geomabs::CurveType;

use super::projector::{Projector, ProjectorState};
use super::{eval_pln_dir2d, eval_pln_pnt2d};

/// `class ProjLib_Plane` (`hxx:27-90`).
pub struct Plane {
    base: ProjectorState,
    /// `myPlane` (`hxx:87`).
    plane: GpPln,
}

impl Plane {
    /// `ProjLib_Plane()` (`cxx:24`): every member default-initialized, so the
    /// state starts as `GeomAbs_BSplineCurve` / not done.
    pub fn new() -> Self {
        Self {
            base: ProjectorState::new(),
            plane: GpPln::default(),
        }
    }

    /// `ProjLib_Plane(const gp_Pln&)` (`cxx:31-34`).
    pub fn from_plane(pl: &GpPln) -> Self {
        let mut p = Self::new();
        p.init(pl);
        p
    }

    /// `ProjLib_Plane(Pl, L)` (`cxx:38-42`).
    pub fn with_line(pl: &GpPln, l: &GpLin) -> Self {
        let mut p = Self::new();
        p.init(pl);
        p.project_lin(l);
        p
    }

    /// `ProjLib_Plane(Pl, C)` (`cxx:46-50`).
    pub fn with_circ(pl: &GpPln, c: &GpCirc) -> Self {
        let mut p = Self::new();
        p.init(pl);
        p.project_circ(c);
        p
    }

    /// `ProjLib_Plane(Pl, E)` (`cxx:54-58`).
    pub fn with_elips(pl: &GpPln, e: &GpElips) -> Self {
        let mut p = Self::new();
        p.init(pl);
        p.project_elips(e);
        p
    }

    /// `ProjLib_Plane(Pl, P)` (`cxx:62-66`).
    pub fn with_parab(pl: &GpPln, pa: &GpParab) -> Self {
        let mut p = Self::new();
        p.init(pl);
        p.project_parab(pa);
        p
    }

    /// `ProjLib_Plane(Pl, H)` (`cxx:70-74`).
    pub fn with_hypr(pl: &GpPln, h: &GpHypr) -> Self {
        let mut p = Self::new();
        p.init(pl);
        p.project_hypr(h);
        p
    }

    /// `Init(Pl)` (`cxx:76-83`).
    pub fn init(&mut self, pl: &GpPln) {
        self.base.curve_type = CurveType::OtherCurve;
        self.base.done = false;
        self.base.periodic = false;
        self.plane = pl.clone();
    }

    /// `myPlane` (`hxx:87`).
    pub fn plane(&self) -> &GpPln {
        &self.plane
    }
}

impl Default for Plane {
    fn default() -> Self {
        Self::new()
    }
}

/// `EvalDir2d(D, Pl)` (`cxx:94-97`). Panics where `gp_Dir2d(X, Y)` raises
/// `Standard_ConstructionError`, i.e. when `d` is along the plane normal.
fn eval_dir2d(d: &GpDir, pl: &GpPln) -> GpDir2d {
    let (dx, dy) = eval_pln_dir2d(pl, d);
    GpDir2d::new(dx, dy).expect("Standard_ConstructionError: gp_Dir2d")
}

/// `gp_Ax22d(P2d, X2d, Y2d)` over the projected axes of `pos`. Panics where
/// the `gp_Ax22d` constructor raises `Standard_ConstructionError` (`X2d` and
/// `Y2d` not perpendicular in the plane frame).
fn projected_axes(pos: &occt_core::gp::GpAx2, pl: &GpPln) -> (GpPnt2d, GpDir2d, GpDir2d) {
    let p2d = eval_pln_pnt2d(pl, &pos.location());
    let x2d = eval_dir2d(pos.x_direction(), pl);
    let y2d = eval_dir2d(pos.y_direction(), pl);
    (p2d, x2d, y2d)
}

/// `gp_Ax22d(P2d, X2d, Y2d)`, panicking like the constructor.
fn axes22d(p2d: GpPnt2d, x2d: GpDir2d, y2d: GpDir2d) -> GpAx22d {
    GpAx22d::new(p2d, x2d, y2d).expect("Standard_ConstructionError: gp_Ax22d")
}

impl Projector for Plane {
    fn state(&self) -> &ProjectorState {
        &self.base
    }

    fn state_mut(&mut self) -> &mut ProjectorState {
        &mut self.base
    }

    /// `Project(const gp_Lin&)` (`cxx:102-107`).
    fn project_lin(&mut self, l: &GpLin) {
        self.base.curve_type = CurveType::Line;
        let p2d = eval_pln_pnt2d(&self.plane, &l.location());
        let d2d = eval_dir2d(&l.direction(), &self.plane);
        self.base.lin = GpLin2d::from_pnt_dir(p2d, d2d);
        self.base.done = true;
    }

    /// `Project(const gp_Circ&)` (`cxx:111-122`).
    fn project_circ(&mut self, c: &GpCirc) {
        self.base.curve_type = CurveType::Circle;

        let pos = c.position();
        let (p2d, x2d, y2d) = projected_axes(&pos, &self.plane);
        let ax = axes22d(p2d, x2d, y2d);

        self.base.circ = occt_core::gp::GpCirc2d::new(ax, c.radius());
        self.base.periodic = true;
        self.base.done = true;
    }

    /// `Project(const gp_Elips&)` (`cxx:126-139`).
    fn project_elips(&mut self, e: &GpElips) {
        self.base.curve_type = CurveType::Ellipse;

        let pos = *e.position();
        let (p2d, x2d, y2d) = projected_axes(&pos, &self.plane);
        let ax = axes22d(p2d, x2d, y2d);

        self.base.elips =
            occt_core::gp::GpElips2d::new(ax, e.major_radius(), e.minor_radius());
        self.base.periodic = true;
        self.base.done = true;
    }

    /// `Project(const gp_Parab&)` (`cxx:143-155`).
    fn project_parab(&mut self, p: &GpParab) {
        self.base.curve_type = CurveType::Parabola;

        let pos = *p.position();
        let (p2d, x2d, y2d) = projected_axes(&pos, &self.plane);
        let ax = axes22d(p2d, x2d, y2d);

        self.base.parab = occt_core::gp::GpParab2d::new(ax, p.focal);
        self.base.done = true;
    }

    /// `Project(const gp_Hypr&)` (`cxx:159-169`).
    fn project_hypr(&mut self, h: &GpHypr) {
        self.base.curve_type = CurveType::Hyperbola;

        let pos = *h.position();
        let (p2d, x2d, y2d) = projected_axes(&pos, &self.plane);
        let ax = axes22d(p2d, x2d, y2d);

        self.base.hypr =
            occt_core::gp::GpHypr2d::new(ax, h.major_radius, h.minor_radius);
        self.base.done = true;
    }
}
