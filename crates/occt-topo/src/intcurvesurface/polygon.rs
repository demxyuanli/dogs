//! `IntCurveSurface_ThePolygonOfHInter` and
//! `IntCurveSurface_ThePolygonToolOfHInter`.
//!
//! Source: `IntCurveSurface_ThePolygonOfHInter.hxx/.cxx` and
//! `IntCurveSurface_PolygonUtils.pxx` (TKGeomAlgo). The polygon is the polyline
//! approximation of the curve that `IntCurveSurface_TheInterferenceOfHInter`
//! weaves against the surface polyhedron, and the `CurveTool` side of
//! `Intf_InterferencePolygonPolyhedron`.
//!
//! Point storage is 1-based, exactly as in OCCT (`ThePnts(1..NbPntIn)`), so the
//! index arithmetic of `ApproxParamOnCurve` and the interference engine carries
//! over unchanged.

use occt_core::bnd::BndBox;
use occt_core::gp::{GpDir, GpLin, GpPnt, GpVec};
use occt_core::intf::IntfPolygon3dTool;
use occt_geom::Curve;

/// `IntCurveSurface_ThePolygonOfHInter` (`...hxx:29-101`).
#[derive(Clone)]
pub struct ThePolygonOfHInter {
    /// `NbPntIn` (`...hxx:96`).
    nb_pnt_in: usize,
    /// `ThePnts` (`...hxx:95`), 1-based: element 0 is unused.
    the_pnts: Vec<GpPnt>,
    /// `TheBnd` (`...hxx:93`).
    the_bnd: BndBox,
    /// `TheDeflection` (`...hxx:94`).
    the_deflection: f64,
    /// `ClosedPolygon` (`...hxx:97`).
    closed_polygon: bool,
    /// `Binf` (`...hxx:98`).
    binf: f64,
    /// `Bsup` (`...hxx:99`).
    bsup: f64,
    /// `myParams` (`...hxx:100`); `None` on the uniform ctor (`...cxx:77-84`).
    /// 1-based when present, same as `ThePnts`.
    my_params: Option<Vec<f64>>,
}

impl ThePolygonOfHInter {
    /// `ThePolygonOfHInter(Curve, NbPnt)` (`...cxx:27-37`): the whole natural
    /// range, with `NbPntIn = max(5, NbPnt)`.
    pub fn new(curve: &dyn Curve, nb_pnt: usize) -> Self {
        let nb_pnt_in = if nb_pnt < 5 { 5 } else { nb_pnt };
        let binf = curve.first_parameter();
        let bsup = curve.last_parameter();
        let mut me = Self {
            nb_pnt_in,
            the_pnts: vec![GpPnt::new(0.0, 0.0, 0.0); nb_pnt_in + 1],
            the_bnd: BndBox::new(),
            the_deflection: 0.0,
            closed_polygon: false,
            binf,
            bsup,
            my_params: None,
        };
        me.init_uniform(curve);
        me
    }

    /// `ThePolygonOfHInter(Curve, U1, U2, NbPnt)` (`...cxx:41-53`): an explicit
    /// parameter window.
    pub fn with_range(curve: &dyn Curve, u1: f64, u2: f64, nb_pnt: usize) -> Self {
        let nb_pnt_in = if nb_pnt < 5 { 5 } else { nb_pnt };
        let mut me = Self {
            nb_pnt_in,
            the_pnts: vec![GpPnt::new(0.0, 0.0, 0.0); nb_pnt_in + 1],
            the_bnd: BndBox::new(),
            the_deflection: 0.0,
            closed_polygon: false,
            binf: u1,
            bsup: u2,
            my_params: None,
        };
        me.init_uniform(curve);
        me
    }

    /// `ThePolygonOfHInter(Curve, Upars)` (`...cxx:57-66`): explicit sample
    /// parameters, remembered in `myParams`.
    pub fn with_params(curve: &dyn Curve, upars: &[f64]) -> Self {
        let nb_pnt_in = upars.len();
        let mut me = Self {
            nb_pnt_in,
            the_pnts: vec![GpPnt::new(0.0, 0.0, 0.0); nb_pnt_in + 1],
            the_bnd: BndBox::new(),
            the_deflection: 0.0,
            closed_polygon: false,
            binf: *upars.first().unwrap_or(&0.0),
            bsup: *upars.last().unwrap_or(&0.0),
            my_params: None,
        };
        me.init_with_params(curve, upars);
        me
    }

    /// `Init(Curve)` (`...cxx:72-83`) =
    /// `PolygonUtils::InitUniform` (`PolygonUtils.pxx:43-86`).
    fn init_uniform(&mut self, curve: &dyn Curve) {
        let du = (self.bsup - self.binf) / (self.nb_pnt_in - 1) as f64;
        let mut u = self.binf;
        for i in 1..=self.nb_pnt_in {
            let p = curve.d0(u);
            self.the_bnd.add_point(&p);
            self.the_pnts[i] = p;
            u += du;
        }
        self.the_deflection = 0.0;
        if self.nb_pnt_in > 3 {
            u = self.binf + du * 0.5;
            for i in 1..self.nb_pnt_in {
                let pm = curve.d0(u);
                let p1 = self.the_pnts[i];
                let p2 = self.the_pnts[i + 1];
                if let Ok(dir) = GpDir::from_vec(&GpVec::from_pnts(&p1, &p2)) {
                    let t = GpLin::from_pnt_dir(p1, dir).distance(&pm);
                    if t > self.the_deflection {
                        self.the_deflection = t;
                    }
                }
                u += du;
            }
            self.the_bnd.enlarge(1.5 * self.the_deflection);
        } else {
            self.the_bnd.enlarge(1e-10);
        }
        self.closed_polygon = false;
    }

    /// `Init(Curve, Upars)` (`...cxx:86-98`) =
    /// `PolygonUtils::InitWithParams` (`PolygonUtils.pxx:101-145`).
    fn init_with_params(&mut self, curve: &dyn Curve, upars: &[f64]) {
        let mut params = vec![0.0; self.nb_pnt_in + 1];
        for i in 1..=self.nb_pnt_in {
            params[i] = upars[i - 1];
            let p = curve.d0(upars[i - 1]);
            self.the_bnd.add_point(&p);
            self.the_pnts[i] = p;
        }
        self.the_deflection = 0.0;
        if self.nb_pnt_in > 3 {
            for i in 1..self.nb_pnt_in {
                let u = 0.5 * (upars[i - 1] + upars[i]);
                let pm = curve.d0(u);
                let p1 = self.the_pnts[i];
                let p2 = self.the_pnts[i + 1];
                if let Ok(dir) = GpDir::from_vec(&GpVec::from_pnts(&p1, &p2)) {
                    let t = GpLin::from_pnt_dir(p1, dir).distance(&pm);
                    if t > self.the_deflection {
                        self.the_deflection = t;
                    }
                }
            }
            self.the_bnd.enlarge(1.5 * self.the_deflection);
        } else {
            self.the_bnd.enlarge(1e-10);
        }
        self.my_params = Some(params);
        self.closed_polygon = false;
    }

    /// Line segment as a two-point polygon. Port-only shortcut for the
    /// `IntCurvesFace_Intersector` quick reject (`Intf_Tool::PolyhedronBox`);
    /// OCCT feeds `ThePolygonOfHInter` a `GeomAdaptor_Curve` on the line
    /// instead, which would sample `max(5, NbPnt)` points over the window.
    pub fn of_line(lin: &GpLin, t0: f64, t1: f64) -> Self {
        let p0 = occt_core::elib::clib::line_value(lin, t0);
        let p1 = occt_core::elib::clib::line_value(lin, t1);
        let mut the_bnd = BndBox::new();
        the_bnd.add_point(&p0);
        the_bnd.add_point(&p1);
        Self {
            nb_pnt_in: 2,
            the_pnts: vec![GpPnt::new(0.0, 0.0, 0.0), p0, p1],
            the_bnd,
            the_deflection: 0.0,
            closed_polygon: false,
            binf: t0,
            bsup: t1,
            my_params: Some(vec![0.0, t0, t1]),
        }
    }

    /// `ApproxParamOnCurve(Index, ParamOnLine)` (`...cxx:104-113` =
    /// `PolygonUtils.pxx:152-201`).
    pub fn approx_param_on_curve(&self, index: i32, param_on_line: f64) -> f64 {
        if !(0.0..=1.0).contains(&param_on_line) {
            return self.binf
                + (param_on_line * (self.bsup - self.binf)) / (self.nb_pnt_in - 1) as f64;
        }
        let mut index = index;
        let mut param_on_line = param_on_line;
        if index as usize == self.nb_pnt_in && param_on_line == 0.0 {
            index -= 1;
            param_on_line = 1.0;
        }
        let (du, u) = match &self.my_params {
            None => {
                let du = (self.bsup - self.binf) / (self.nb_pnt_in - 1) as f64;
                (du, self.binf + du * (index - 1) as f64)
            }
            Some(params) => {
                let du = params[index as usize + 1] - params[index as usize];
                (du, params[index as usize])
            }
        };
        u + du * param_on_line
    }

    /// `Bounding()` (`...hxx:46-48`).
    pub fn bounding(&self) -> &BndBox {
        &self.the_bnd
    }

    /// `DeflectionOverEstimation()` (`...hxx:50`).
    pub fn deflection_over_estimation(&self) -> f64 {
        self.the_deflection
    }

    /// `SetDeflectionOverEstimation(x)` (`...hxx:52-56`).
    pub fn set_deflection_over_estimation(&mut self, x: f64) {
        self.the_deflection = x;
        self.the_bnd.enlarge(x);
    }

    /// `Closed()` (`...hxx:60`).
    pub fn is_closed(&self) -> bool {
        self.closed_polygon
    }

    /// `Closed(flag)` (`...hxx:58`).
    pub fn set_closed(&mut self, flag: bool) {
        self.closed_polygon = flag;
    }

    /// `NbSegments()` (`...hxx:66`).
    pub fn nb_segments(&self) -> usize {
        self.nb_pnt_in - 1
    }

    /// `BeginOfSeg(Index)` (`...hxx:69`).
    pub fn begin_of_seg(&self, index: i32) -> GpPnt {
        self.the_pnts[index as usize]
    }

    /// `EndOfSeg(Index)` (`...hxx:72`).
    pub fn end_of_seg(&self, index: i32) -> GpPnt {
        self.the_pnts[index as usize + 1]
    }

    /// `InfParameter()` (`...hxx:75-76`).
    pub fn inf_parameter(&self) -> f64 {
        self.binf
    }

    /// `SupParameter()` (`...hxx:78-79`).
    pub fn sup_parameter(&self) -> f64 {
        self.bsup
    }

    /// Port-only accessor for the sampled point count (`NbPntIn`).
    pub fn nb_points(&self) -> usize {
        self.nb_pnt_in
    }
}

/// `IntCurveSurface_ThePolygonToolOfHInter` (`...hxx:31-67`).
pub struct ThePolygonToolOfHInter;

impl IntfPolygon3dTool for ThePolygonToolOfHInter {
    type Polygon3d = ThePolygonOfHInter;

    fn bounding(the_polyg: &Self::Polygon3d) -> &BndBox {
        the_polyg.bounding()
    }

    fn deflection_over_estimation(the_polyg: &Self::Polygon3d) -> f64 {
        the_polyg.deflection_over_estimation()
    }

    fn closed(the_polyg: &Self::Polygon3d) -> bool {
        the_polyg.is_closed()
    }

    fn nb_segments(the_polyg: &Self::Polygon3d) -> usize {
        the_polyg.nb_segments()
    }

    fn begin_of_seg(the_polyg: &Self::Polygon3d, index: usize) -> GpPnt {
        the_polyg.begin_of_seg(index as i32)
    }

    fn end_of_seg(the_polyg: &Self::Polygon3d, index: usize) -> GpPnt {
        the_polyg.end_of_seg(index as i32)
    }
}
