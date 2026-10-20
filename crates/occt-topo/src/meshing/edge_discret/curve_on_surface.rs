//! `BRepMesh_CurveOnSurface` — `Curve` adapter evaluating a pcurve on its surface.

use super::prelude::*;
use super::*;

/// `Adaptor3d_CurveOnSurface` used by `BRepAdaptor_Curve(edge, face)`.
pub struct CurveOnSurface {
    pcurve: Arc<dyn Curve2d>,
    surface: Arc<dyn Surface>,
    first: f64,
    last: f64,
    /// `myCirc` (`cxx:1561-1564` plane arm; `cxx:1580-1733` iso arm).
    kpart_circ: Option<GpCirc>,
    /// `myLin` (`cxx:1565-1577`, the plane arm `D1(0)` projection).
    kpart_lin: Option<GpLin>,
    /// `Ellipse()` (`cxx:1399-1402`), reachable only on a plane surface.
    kpart_elips: Option<GpElips>,
    /// `Parabola()` (`cxx:1413-1416`), reachable only on a plane surface.
    kpart_parab: Option<GpParab>,
    /// `Hyperbola()` (`cxx:1406-1409`), reachable only on a plane surface.
    kpart_hypr: Option<GpHypr>,
    /// `Bezier()` (`cxx:1458-1485`) / `BSpline()` (`cxx:1489-1520`), reachable
    /// only on a plane surface; `EvalKPart`'s plane arm sets `myType = CTy`
    /// (`cxx:1560`), so a pole-classified pcurve answers these.
    kpart_poles: Option<Arc<dyn Curve>>,
    /// `myFirstSurf` from `EvalFirstLastSurf` (`cxx:1778-1783`).
    first_surf: Option<Arc<dyn Surface>>,
    /// `myLastSurf` from `EvalFirstLastSurf` (`cxx:1822-1827`).
    last_surf: Option<Arc<dyn Surface>>,
    /// `GeomLib_CurveOnSurfaceEvaluator` trim cache (`GeomLib.cxx:1009-1016`):
    /// `(FirstParam, LastParam, TrimCurve)`. `Mutex` keeps `Curve: Sync`.
    trim_cache: std::sync::Mutex<Option<(f64, f64, Arc<CurveOnSurface>)>>,
}

impl CurveOnSurface {
    /// Evaluate `surface(pcurve(t))` on `[first, last]`, with `EvalKPart` and
    /// `EvalFirstLastSurf` (`Adaptor3d_CurveOnSurface::Load` cxx:951-963).
    ///
    /// PORTED: Offset unwrap + `LocatePart_Offset` (BSpline + RevExt) +
    /// `LocatePart_RevExt` + `LocatePart` / `Locate1Coord` (surface+curve) /
    /// `Locate2Coord` (Arr+param) + UTrim/VTrim end patches for D1/D2
    /// (`cxx:1212-1266`, `1833-1866`) + `GeomAdaptor` LocalD1/IfUVBound/Span via
    /// `GeomRectangularTrimmedSurface` (`GeomAdaptor_Surface.cxx:1193-1195`).
    pub fn new(
        pcurve: Arc<dyn Curve2d>,
        surface: Arc<dyn Surface>,
        first: f64,
        last: f64,
    ) -> Self {
        let kpart = eval_k_part(pcurve.as_ref(), surface.as_ref());
        let (first_surf, last_surf) =
            super::super::cos_locate::eval_first_last_surf(pcurve.as_ref(), &surface, first, last);
        Self {
            pcurve,
            surface,
            first,
            last,
            kpart_circ: kpart.circ,
            kpart_lin: kpart.lin,
            kpart_elips: kpart.elips,
            kpart_parab: kpart.parab,
            kpart_hypr: kpart.hypr,
            kpart_poles: kpart.poles,
            first_surf,
            last_surf,
            // `cxx:1136`: seed outside the true range so the first Evaluate re-trims.
            trim_cache: std::sync::Mutex::new(None),
        }
    }

    /// `Adaptor3d_CurveOnSurface::Trim` (`cxx:1133-1140`) for the evaluator.
    fn trim_for_eval(&self, first: f64, last: f64) -> Arc<CurveOnSurface> {
        let mut cache = self.trim_cache.lock().unwrap_or_else(|e| e.into_inner());
        if let Some((f, l, trimmed)) = cache.as_ref() {
            if *f == first && *l == last {
                return trimmed.clone();
            }
        }
        let trimmed_pc: Arc<dyn Curve2d> = Arc::new(occt_geom2d::trimmed::Geom2dTrimmedCurve::new(
            Arc::from(self.pcurve.clone_dyn()),
            first,
            last,
        ));
        let trimmed = Arc::new(Self::new(
            trimmed_pc,
            self.surface.clone(),
            first,
            last,
        ));
        *cache = Some((first, last, trimmed.clone()));
        trimmed
    }

    fn end_surf(&self, u: f64) -> Option<&Arc<dyn Surface>> {
        // `Adaptor3d_CurveOnSurface::EvalD1` cxx:1208-1222.
        let tol = PCONFUSION / 10.0;
        if (u - self.first).abs() < tol {
            self.first_surf.as_ref()
        } else if (u - self.last).abs() < tol {
            self.last_surf.as_ref()
        } else {
            None
        }
    }
}

impl Curve for CurveOnSurface {
    fn d0(&self, u: f64) -> GpPnt {
        if let Some(ref c) = self.kpart_circ {
            return clib::circle_value(c, u);
        }
        if let Some(ref l) = self.kpart_lin {
            return clib::line_value(l, u);
        }
        let uv = self.pcurve.d0(u);
        self.surface.d0(uv.x(), uv.y())
    }

    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        if let Some(ref c) = self.kpart_circ {
            return clib::circle_d1(c, u);
        }
        if let Some(ref l) = self.kpart_lin {
            return clib::line_d1(l, u);
        }
        let (uv, duv) = self.pcurve.d1(u);
        let surf = self.end_surf(u).unwrap_or(&self.surface);
        let (p, su, sv) = surf.d1(uv.x(), uv.y());
        let tan = su
            .multiplied_scalar(duv.x())
            .added(&sv.multiplied_scalar(duv.y()));
        (p, tan)
    }

    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        if let Some(ref c) = self.kpart_circ {
            return clib::circle_d2(c, u);
        }
        if let Some(ref l) = self.kpart_lin {
            return clib::line_d2(l, u);
        }
        let (uv, duv, d2uv) = self.pcurve.d2(u);
        let surf = self.end_surf(u).unwrap_or(&self.surface);
        let (p, su, sv, suu, suv, svv) = surf.d2(uv.x(), uv.y());
        let d1 = su
            .multiplied_scalar(duv.x())
            .added(&sv.multiplied_scalar(duv.y()));
        let d2 = suu
            .multiplied_scalar(duv.x() * duv.x())
            .added(&suv.multiplied_scalar(2.0 * duv.x() * duv.y()))
            .added(&svv.multiplied_scalar(duv.y() * duv.y()))
            .added(&su.multiplied_scalar(d2uv.x()))
            .added(&sv.multiplied_scalar(d2uv.y()));
        (p, d1, d2)
    }

    fn first_parameter(&self) -> f64 {
        self.first
    }
    fn last_parameter(&self) -> f64 {
        self.last
    }
    fn is_periodic(&self) -> bool {
        if self.kpart_circ.is_some() {
            return true;
        }
        self.pcurve.is_periodic()
    }
    fn period(&self) -> f64 {
        if self.kpart_circ.is_some() {
            return 2.0 * std::f64::consts::PI;
        }
        self.pcurve.period()
    }
    fn circle_radius(&self) -> Option<f64> {
        self.kpart_circ.as_ref().map(|c| c.radius())
    }
    fn gp_circ(&self) -> Option<GpCirc> {
        self.kpart_circ.clone()
    }
    /// `Adaptor3d_CurveOnSurface::Line` (`cxx:1381-1386`): `myLin` as a real
    /// `gp_Lin`, so `EvalKPart`'s `GeomAbs_Line` cos keeps that type in
    /// `GeomAdaptor_Curve::GetType` (`extrema_pc/curve_tool.rs:36-56`).
    fn gp_line(&self) -> Option<GpLin> {
        self.kpart_lin
    }
    /// `Adaptor3d_CurveOnSurface::Ellipse` (`cxx:1399-1402`):
    /// `to3d(mySurface->Plane(), myCurve->Ellipse())`.
    fn gp_ellipse(&self) -> Option<GpElips> {
        self.kpart_elips.clone()
    }
    /// `Adaptor3d_CurveOnSurface::Parabola` (`cxx:1413-1416`).
    fn gp_parabola(&self) -> Option<GpParab> {
        self.kpart_parab.clone()
    }
    /// `Adaptor3d_CurveOnSurface::Hyperbola` (`cxx:1406-1409`).
    fn gp_hyperbola(&self) -> Option<GpHypr> {
        self.kpart_hypr.clone()
    }
    /// `Adaptor3d_CurveOnSurface::Bezier` (`cxx:1458-1485`): the plane-lifted
    /// pcurve poles, so `EvalKPart`'s `GeomAbs_BezierCurve` cos keeps that type
    /// in `GeomAdaptor_Curve::GetType` (`extrema_pc/curve_tool.rs:36-56`).
    fn bezier_poles(&self) -> Option<&[GpPnt]> {
        self.kpart_poles.as_ref()?.bezier_poles()
    }
    /// `Adaptor3d_CurveOnSurface::Bezier`'s `IsRational()` arm (`cxx:1475-1479`).
    fn bezier_weights(&self) -> Option<&[f64]> {
        self.kpart_poles.as_ref()?.bezier_weights()
    }
    /// `Adaptor3d_CurveOnSurface::BSpline` (`cxx:1489-1520`).
    fn bspline_poles(&self) -> Option<&[GpPnt]> {
        self.kpart_poles.as_ref()?.bspline_poles()
    }
    /// `Adaptor3d_CurveOnSurface::BSpline`'s `IsRational()` arm (`cxx:1509-1513`).
    fn bspline_weights(&self) -> Option<&[f64]> {
        self.kpart_poles.as_ref()?.bspline_weights()
    }
    /// `myCurve->Knots()` (`cxx:1506`).
    fn bspline_knots(&self) -> Option<&[f64]> {
        self.kpart_poles.as_ref()?.bspline_knots()
    }
    /// `myCurve->Degree()` (`cxx:1512`/`:1516`; `Adaptor3d_CurveOnSurface::Degree`
    /// `cxx:1418-1434` is `myCurve->Degree()`).
    fn nurbs_degree(&self) -> Option<usize> {
        self.kpart_poles.as_ref()?.nurbs_degree()
    }
    fn is_line(&self) -> bool {
        self.kpart_lin.is_some()
    }
    /// `Adaptor3d_CurveOnSurface::Resolution` (`Adaptor3d_CurveOnSurface.cxx:1364-1370`):
    /// `myCurve->Resolution(min(mySurface->UResolution(R3d), mySurface->VResolution(R3d)))`.
    /// The trait default is `GeomAdaptor_Curve::Resolution` (`curve.rs:133-137`), which the
    /// cxx overrides here; both surface side values come from the single faithful port of
    /// `GeomAdaptor_Surface::U/VResolution` (`approx_same_parameter.rs:386`, `:453`).
    fn resolution(&self, r3d: f64) -> f64 {
        let ru = occt_geom::approx_same_parameter::u_resolution(self.surface.as_ref(), r3d);
        let rv = occt_geom::approx_same_parameter::v_resolution(self.surface.as_ref(), r3d);
        self.pcurve.resolution(ru.min(rv))
    }
    fn continuity(&self) -> u8 {
        self.pcurve.continuity().min(self.surface.continuity())
    }

    /// `Adaptor3d_CurveOnSurface::Intervals` (`cxx:1119-1128` via `NbIntervals`
    /// `:1045-1114`): pcurve breaks plus roots on interior surface U/V knots.
    fn parameter_intervals(&self, continuity: u8) -> Vec<f64> {
        super::super::cos_intervals::curve_on_surface_intervals(
            self.pcurve.as_ref(),
            self.surface.as_ref(),
            continuity,
        )
    }

    fn nb_intervals(&self, continuity: u8) -> i32 {
        self.parameter_intervals(continuity)
            .len()
            .saturating_sub(1)
            .max(1) as i32
    }

    /// `GeomLib_CurveOnSurfaceEvaluator::Evaluate` (`GeomLib.cxx:1001-1047`).
    fn approx_eval(&self, first: f64, last: f64, u: f64, order: i32, out: &mut [f64]) -> i32 {
        let trimmed = self.trim_for_eval(first, last);
        match order {
            0 => {
                let p = trimmed.d0(u);
                out[0] = p.x();
                out[1] = p.y();
                out[2] = p.z();
            }
            1 => {
                let (_, v) = trimmed.d1(u);
                out[0] = v.x();
                out[1] = v.y();
                out[2] = v.z();
            }
            2 => {
                let (_, _, v) = trimmed.d2(u);
                out[0] = v.x();
                out[1] = v.y();
                out[2] = v.z();
            }
            _ => {
                out[0] = 0.0;
                out[1] = 0.0;
                out[2] = 0.0;
                return 3;
            }
        }
        0
    }

    fn transform(&mut self, _t: &occt_core::gp::GpTrsf) {}
    fn reverse(&mut self) {
        let tmp = self.first;
        self.first = self.last;
        self.last = tmp;
        if let Ok(mut cache) = self.trim_cache.lock() {
            *cache = None;
        }
    }
    fn clone_dyn(&self) -> Box<dyn Curve> {
        Box::new(Self {
            pcurve: Arc::from(self.pcurve.clone_dyn()),
            surface: self.surface.clone(),
            first: self.first,
            last: self.last,
            kpart_circ: self.kpart_circ.clone(),
            kpart_lin: self.kpart_lin,
            kpart_elips: self.kpart_elips.clone(),
            kpart_parab: self.kpart_parab.clone(),
            kpart_hypr: self.kpart_hypr.clone(),
            kpart_poles: self.kpart_poles.clone(),
            first_surf: self.first_surf.clone(),
            last_surf: self.last_surf.clone(),
            trim_cache: std::sync::Mutex::new(None),
        })
    }
}
