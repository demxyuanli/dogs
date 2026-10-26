//! `IntPatch_PolyLine` - the 2D polygon of an intersection line viewed on one
//! of the two surfaces.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_PolyLine.hxx` and
//! `.cxx`.
//!
//! Ported members (OCCT 8.0.0):
//! * default ctor (`cxx:29-34`) and `(InitDefle)` ctor (`cxx:37-42`); the
//!   default deflection is `Precision::PConfusion() * 100.` (`cxx:26`).
//! * `SetWLine` (`cxx:45-52`), `SetRLine` (`cxx:55-62`).
//! * `Prepare` (`cxx:65-151`): the bounding box and the deflection estimate,
//!   including the parabolic refinement of a too-large 3-point deflection.
//! * `ResetError` (`cxx:154-157`), `NbPoints` (`cxx:161-164`),
//!   `Point` (`cxx:168-220`).
//!
//! Translation notes:
//! * OCCT stores `occ::handle<IntPatch_WLine>` / `<IntPatch_RLine>`
//!   (`hxx:56-58`); the port borrows the line for the lifetime of the polygon,
//!   which is how OCCT's callers use it (`IntPatch_RstInt.cxx:462-489` builds a
//!   local `IntPatch_PolyLine` over the line it is processing).
//! * `typ` (`hxx:54`) is `IntPatch_IType`; the port spells the two live cases
//!   as the [`PolyLineSource`] variant.
//! * the `pnt` member (`hxx:53`) is written nowhere in `IntPatch_PolyLine.cxx`
//!   (it is a leftover of the pre-1993 implementation); it is not carried.

use occt_core::bnd::box2d::BndBox2d;
use occt_core::gp::GpPnt2d;
use occt_core::intf::IntfPolygon2d;
use occt_core::precision::Precision;

use crate::geom_int::RLine;
use crate::int_tools_wline::WLine;
use crate::intpatch_polygo::IntPatchPolygo;

/// `INITDEFLE` (`cxx:26`): `Precision::PConfusion() * 100.`
pub fn init_defle() -> f64 {
    Precision::PCONFUSION * 100.0
}

/// The `typ` + `wpoly` / `rpoly` triple of `IntPatch_PolyLine` (`hxx:54-58`).
#[derive(Clone, Copy)]
pub enum PolyLineSource<'a> {
    /// `typ == IntPatch_Walking`, `wpoly = Line`.
    WLine(&'a WLine),
    /// `typ == IntPatch_Restriction`, `rpoly = Line`.
    RLine(&'a RLine),
}

/// `IntPatch_PolyLine` (`hxx:32-59`).
pub struct IntPatchPolyLine<'a> {
    /// `myBox` (`Intf_Polygon2d::myBox`).
    my_box: BndBox2d,
    /// `myError` (`hxx:55`).
    my_error: f64,
    /// `onfirst` (`hxx:55`).
    on_first: bool,
    /// `typ` + `wpoly` / `rpoly`.
    source: Option<PolyLineSource<'a>>,
}

impl<'a> IntPatchPolyLine<'a> {
    /// `IntPatch_PolyLine()` (`cxx:29-34`): `IntPatch_Polygo(INITDEFLE)`,
    /// `onfirst(false)`.
    pub fn new() -> Self {
        Self::with_deflection(init_defle())
    }

    /// `IntPatch_PolyLine(const double InitDefle)` (`cxx:37-42`).
    pub fn with_deflection(init_defle: f64) -> Self {
        Self {
            my_box: BndBox2d::new(),
            my_error: init_defle,
            on_first: false,
            source: None,
        }
    }

    /// `SetWLine(OnFirst, Line)` (`cxx:45-52`).
    pub fn set_wline(&mut self, on_first: bool, line: &'a WLine) {
        self.source = Some(PolyLineSource::WLine(line));
        self.on_first = on_first;
        self.prepare();
    }

    /// `SetRLine(OnFirst, Line)` (`cxx:55-62`).
    pub fn set_rline(&mut self, on_first: bool, line: &'a RLine) {
        self.source = Some(PolyLineSource::RLine(line));
        self.on_first = on_first;
        self.prepare();
    }

    /// `Prepare()` (`cxx:65-151`).
    fn prepare(&mut self) {
        self.my_box.set_void();
        let n = self.nb_points();
        let eps_2 = self.my_error * self.my_error;

        let mut p1 = GpPnt2d::zero();
        let mut p2 = GpPnt2d::zero();
        if n >= 3 {
            p1 = self.point(1);
            p2 = self.point(2);
        }
        for i in 1..=n {
            let p3 = self.point(i);
            if i >= 3 {
                let v13 = p3.xy().subtracted(p1.xy());
                let v12 = p2.xy().subtracted(p1.xy());
                let d13_2 = v13.square_modulus();
                let mut d_2;
                if d13_2 > eps_2 {
                    // `d_2 = V13.CrossSquareMagnitude(V12) / d13_2;`
                    let cross = v13.crossed(&v12);
                    d_2 = cross * cross / d13_2;
                } else {
                    d_2 = eps_2;
                }
                if d_2 > self.my_error * self.my_error {
                    // try to compute deflection more precisely using parabola
                    // interpolation (`cxx:95-137`)
                    let v23 = p3.xy().subtracted(p2.xy());
                    let d12 = v12.modulus();
                    let d23 = v23.modulus();
                    let mut tm = d12 / (d12 + d23);
                    if tm > 0.1 && tm < 0.9 {
                        tm -= (tm - 0.5) * 0.6;
                        let tm1mtm = tm * (1.0 - tm);
                        // coefficients of the parabola
                        let ax = (tm * v13.x() - v12.x()) / tm1mtm;
                        let bx = (v12.x() - tm * tm * v13.x()) / tm1mtm;
                        let cx = p1.x();
                        let ay = (tm * v13.y() - v12.y()) / tm1mtm;
                        let by = (v12.y() - tm * tm * v13.y()) / tm1mtm;
                        let cy = p1.y();
                        // equations of lines P1-P2 and P2-P3
                        let a1 = v12.y() / d12;
                        let b1 = -v12.x() / d12;
                        let c1 = (p2.x() * p1.y() - p1.x() * p2.y()) / d12;
                        let a2 = v23.y() / d23;
                        let b2 = -v23.x() / d23;
                        let c2 = (p3.x() * p2.y() - p2.x() * p3.y()) / d23;
                        // points on the parabola with max deflection
                        let t1 = -0.5 * (a1 * bx + b1 * by) / (a1 * ax + b1 * ay);
                        let t2 = -0.5 * (a2 * bx + b2 * by) / (a2 * ax + b2 * ay);
                        let xt1 = ax * t1 * t1 + bx * t1 + cx;
                        let yt1 = ay * t1 * t1 + by * t1 + cy;
                        let xt2 = ax * t2 * t2 + bx * t2 + cx;
                        let yt2 = ay * t2 * t2 + by * t2 + cy;
                        // max deflection on segments P1-P2 and P2-P3
                        let mut d1 = (a1 * xt1 + b1 * yt1 + c1).abs();
                        let d2 = (a2 * xt2 + b2 * yt2 + c2).abs();
                        if d2 > d1 {
                            d1 = d2;
                        }
                        // select min deflection from linear and parabolic ones
                        if d1 * d1 < d_2 {
                            d_2 = d1 * d1;
                        }
                    }
                    if d_2 > self.my_error * self.my_error {
                        self.my_error = d_2.sqrt();
                    }
                }
                p1 = p2;
                p2 = p3;
            }
            self.my_box.add_point(&p3);
        }
        self.my_box.enlarge(self.my_error);
    }

    /// `ResetError()` (`cxx:154-157`).
    pub fn reset_error(&mut self) {
        self.my_error = init_defle();
    }

    /// `NbPoints()` (`cxx:161-164`):
    /// `(typ == IntPatch_Walking ? wpoly->NbPnts() : rpoly->NbPnts())`.
    pub fn nb_points(&self) -> usize {
        match self.source {
            Some(PolyLineSource::WLine(w)) => w.nb_pnts().max(0) as usize,
            Some(PolyLineSource::RLine(r)) => r.nb_pnts().max(0) as usize,
            None => 0,
        }
    }

    /// `Point(Index)` (`cxx:168-220`).
    ///
    /// The `1e-7 * (X - X1)` shift on the first and last point of a walking
    /// line (`cxx:180-187`, `:201-210`) is kept verbatim: it is what keeps the
    /// two end points off the restriction they sit on.
    pub fn point(&self, index: usize) -> GpPnt2d {
        let (x, y) = self.parameters(index);
        let mut dx = 0.0;
        let mut dy = 0.0;
        if let Some(PolyLineSource::WLine(w)) = self.source {
            let n = w.nb_pnts();
            if index == 1 {
                let (x1, y1) = self.parameters(2);
                dx = 0.0000001 * (x - x1);
                dy = 0.0000001 * (y - y1);
            } else if index as i32 == n {
                let (x1, y1) = self.parameters(index - 1);
                dx = 0.0000001 * (x - x1);
                dy = 0.0000001 * (y - y1);
            }
        }
        GpPnt2d::new(x + dx, y + dy)
    }

    /// The `ParametersOnS1` / `ParametersOnS2` read of `Point`
    /// (`cxx:171-215`) with the `onfirst` branch applied.
    fn parameters(&self, index: usize) -> (f64, f64) {
        match self.source {
            Some(PolyLineSource::WLine(w)) => {
                let p = w.point(index as i32);
                if self.on_first {
                    p.parameters_on_s1()
                } else {
                    p.parameters_on_s2()
                }
            }
            Some(PolyLineSource::RLine(r)) => {
                let p = r.point(index as i32);
                if self.on_first {
                    p.parameters_on_s1()
                } else {
                    p.parameters_on_s2()
                }
            }
            None => (0.0, 0.0),
        }
    }
}

impl Default for IntPatchPolyLine<'_> {
    fn default() -> Self {
        Self::new()
    }
}

impl IntfPolygon2d for IntPatchPolyLine<'_> {
    /// `Intf_Polygon2d::Bounding()`.
    fn bounding(&self) -> &BndBox2d {
        &self.my_box
    }

    /// `DeflectionOverEstimation()` (`lxx:24-27`): `myError`.
    fn deflection_over_estimation(&self) -> f64 {
        self.my_error
    }

    /// `NbSegments()` (`lxx:31-34`): `NbPoints() - 1`.
    fn nb_segments(&self) -> usize {
        self.nb_points().saturating_sub(1)
    }

    /// `Segment(Index, theBegin, theEnd)` (`lxx:38-43`).
    fn segment(&self, index: usize, begin: &mut GpPnt2d, end: &mut GpPnt2d) {
        *begin = self.point(index);
        *end = self.point(index + 1);
    }
}

impl IntPatchPolygo for IntPatchPolyLine<'_> {
    fn error(&self) -> f64 {
        self.my_error
    }

    fn nb_points(&self) -> usize {
        IntPatchPolyLine::nb_points(self)
    }

    fn point(&self, index: usize) -> GpPnt2d {
        IntPatchPolyLine::point(self, index)
    }
}
