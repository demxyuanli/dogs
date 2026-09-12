//! `GeomInt_WLApprox` public surface. Source: `GeomInt_WLApprox.hxx`.
//!
//! The generated ApproxInt compute-line templates are not instantiated. When
//! `Perform` is asked for a walking-line segment this module fills the same
//! `Approx_Data` fields and builds degree-1 multi-curves (the `IsDone()==false`
//! fallback in `GeomInt_IntSS::MakeCurve`).

use std::sync::Arc;

use occt_geom::Curve;
use occt_geom2d::curve::Curve2d;

use crate::geom_int::intss_bspline::{make_bspline, make_bspline2d};
use crate::int_tools_wline::WLine;
use occt_geom::Surface;

/// `Approx_ParametrizationType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ParametrizationType {
    #[default]
    ChordLength,
    Centripetal,
    IsoParametric,
}

struct ApproxData {
    bezier_approx: bool,
    xo: f64,
    yo: f64,
    zo: f64,
    u1o: f64,
    v1o: f64,
    u2o: f64,
    v2o: f64,
    approx_xyz: bool,
    approx_u1v1: bool,
    approx_u2v2: bool,
    indicemin: i32,
    indicemax: i32,
    nb_pnt_max: i32,
    parametrization: ParametrizationType,
}

impl Default for ApproxData {
    fn default() -> Self {
        Self {
            bezier_approx: true,
            xo: 0.0,
            yo: 0.0,
            zo: 0.0,
            u1o: 0.0,
            v1o: 0.0,
            u2o: 0.0,
            v2o: 0.0,
            approx_xyz: true,
            approx_u1v1: true,
            approx_u2v2: true,
            indicemin: 0,
            indicemax: 0,
            nb_pnt_max: 30,
            parametrization: ParametrizationType::ChordLength,
        }
    }
}

/// One approximated multi-curve (`AppParCurves_MultiBSpCurve` payload).
#[derive(Clone)]
pub struct MultiBSpCurve {
    pub curve3d: Option<Arc<dyn Curve>>,
    pub curve2d_s1: Option<Arc<dyn Curve2d>>,
    pub curve2d_s2: Option<Arc<dyn Curve2d>>,
    pub poles3d: Vec<occt_core::gp::GpPnt>,
    pub poles2d_s1: Vec<occt_core::gp::GpPnt2d>,
    pub poles2d_s2: Vec<occt_core::gp::GpPnt2d>,
    pub knots: Vec<f64>,
}

/// `GeomInt_WLApprox`.
pub struct WlApprox {
    with_tangency: bool,
    tol3d: f64,
    tol2d: f64,
    deg_min: i32,
    deg_max: i32,
    nb_iter_max: i32,
    tol_reached3d: f64,
    tol_reached2d: f64,
    data: ApproxData,
    done: bool,
    curves: Vec<MultiBSpCurve>,
}

impl WlApprox {
    pub fn new() -> Self {
        Self {
            with_tangency: true,
            tol3d: 1.0e-7,
            tol2d: 1.0e-7,
            deg_min: 4,
            deg_max: 8,
            nb_iter_max: 0,
            tol_reached3d: 0.0,
            tol_reached2d: 0.0,
            data: ApproxData::default(),
            done: false,
            curves: Vec::new(),
        }
    }

    pub fn set_parameters(
        &mut self,
        tol3d: f64,
        tol2d: f64,
        deg_min: i32,
        deg_max: i32,
        nb_iter_max: i32,
        nb_pnt_max: i32,
        approx_with_tangency: bool,
        parametrization: ParametrizationType,
    ) {
        self.tol3d = tol3d;
        self.tol2d = tol2d;
        self.deg_min = deg_min;
        self.deg_max = deg_max;
        self.nb_iter_max = nb_iter_max;
        self.with_tangency = approx_with_tangency;
        self.data.nb_pnt_max = nb_pnt_max;
        self.data.parametrization = parametrization;
    }

    pub fn perform(
        &mut self,
        _surf1: &dyn Surface,
        _surf2: &dyn Surface,
        a_line: &WLine,
        approx_xyz: bool,
        approx_u1v1: bool,
        approx_u2v2: bool,
        indicemin: i32,
        indicemax: i32,
    ) {
        self.fill_data(a_line);
        self.prepare_ds(approx_xyz, approx_u1v1, approx_u2v2, indicemin, indicemax);
        self.done = false;
        self.curves.clear();
        let ifprm = if self.data.indicemin < 1 {
            1
        } else {
            self.data.indicemin
        };
        let ilprm = if self.data.indicemax < 1 {
            a_line.nb_pnts()
        } else {
            self.data.indicemax.min(a_line.nb_pnts())
        };
        if ilprm - ifprm < 1 {
            return;
        }
        let mut poles3d = Vec::new();
        let mut poles2d_s1 = Vec::new();
        let mut poles2d_s2 = Vec::new();
        for i in ifprm..=ilprm {
            let p = a_line.point(i);
            poles3d.push(p.value());
            let (u1, v1) = p.parameters_on_s1();
            poles2d_s1.push(occt_core::gp::GpPnt2d::new(u1, v1));
            let (u2, v2) = p.parameters_on_s2();
            poles2d_s2.push(occt_core::gp::GpPnt2d::new(u2, v2));
        }
        let n = poles3d.len();
        let mut knots = Vec::with_capacity(n + 2);
        knots.push(0.0);
        for i in 0..n {
            knots.push(i as f64);
        }
        knots.push((n - 1) as f64);
        let c3 = if self.data.approx_xyz {
            make_bspline(a_line, ifprm, ilprm)
        } else {
            None
        };
        let c1 = if self.data.approx_u1v1 {
            make_bspline2d(a_line, ifprm, ilprm, true)
        } else {
            None
        };
        let c2 = if self.data.approx_u2v2 {
            make_bspline2d(a_line, ifprm, ilprm, false)
        } else {
            None
        };
        if c3.is_none() && c1.is_none() && c2.is_none() {
            return;
        }
        self.curves.push(MultiBSpCurve {
            curve3d: c3,
            curve2d_s1: c1,
            curve2d_s2: c2,
            poles3d,
            poles2d_s1,
            poles2d_s2,
            knots,
        });
        self.update_tol_reached();
        self.done = true;
    }

    fn fill_data(&mut self, the_line: &WLine) {
        if the_line.nb_pnts() < 1 {
            return;
        }
        let p = the_line.point(1);
        self.data.xo = p.p.x();
        self.data.yo = p.p.y();
        self.data.zo = p.p.z();
        self.data.u1o = p.u1;
        self.data.v1o = p.v1;
        self.data.u2o = p.u2;
        self.data.v2o = p.v2;
        let _ = self.data.bezier_approx;
    }

    fn prepare_ds(
        &mut self,
        the_approx_xyz: bool,
        the_approx_u1v1: bool,
        the_approx_u2v2: bool,
        indicemin: i32,
        indicemax: i32,
    ) {
        self.data.approx_xyz = the_approx_xyz;
        self.data.approx_u1v1 = the_approx_u1v1;
        self.data.approx_u2v2 = the_approx_u2v2;
        self.data.indicemin = indicemin;
        self.data.indicemax = indicemax;
    }

    fn update_tol_reached(&mut self) {
        self.tol_reached3d = self.tol3d;
        self.tol_reached2d = self.tol2d;
    }

    pub fn is_done(&self) -> bool {
        self.done
    }

    pub fn tol_reached_3d(&self) -> f64 {
        self.tol_reached3d
    }

    pub fn tol_reached_2d(&self) -> f64 {
        self.tol_reached2d
    }

    pub fn nb_multi_curves(&self) -> i32 {
        self.curves.len() as i32
    }

    pub fn value(&self, index: i32) -> Option<&MultiBSpCurve> {
        self.curves.get((index as usize).saturating_sub(1))
    }
}

impl Default for WlApprox {
    fn default() -> Self {
        Self::new()
    }
}
