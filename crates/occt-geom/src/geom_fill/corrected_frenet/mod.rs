//! `GeomFill_CorrectedFrenet` — the Frenet law whose normal is rotated by a
//! torsion-minimising angle law.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_CorrectedFrenet.hxx` and
//! `GeomFill_CorrectedFrenet.cxx:304-1870` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `EvolAroundT`, `TLaw`, `AT`, `AN` and the OCC78 arrays are mutated through
//!   the shared handle (`SetCurve` / `SetInterval`), so they live in
//!   [`CorrectedState`] behind a `Mutex`; see
//!   [`crate::geom_fill::trihedron_law`] for the handle convention.
//! - `SetInterval` keeps `TLaw = EvolAroundT->Trim(...)` when the law is not a
//!   plain Frenet law (`cxx:964-974`).
//! - `EvaluateBestMode` (`cxx:977-1035`) returns [`GeomFillTrihedron`].

mod init;

use std::sync::{Arc, Mutex};

use occt_core::gp::{GpVec, GpXyz};
use occt_core::precision::PCONFUSION;

use crate::curve::Curve;
use crate::extrema_pc::{extrema_curve_type, GeomAbsCurveType};
use crate::geom_fill::frenet::Frenet;
use crate::geom_fill::trihedron::GeomFillTrihedron;
use crate::geom_fill::trihedron_law::{TrihedronD1, TrihedronD2, TrihedronLaw, TrihedronLawCurve};
use crate::law::function::LawFunction;

pub(crate) use init::CorrectedState;

/// `GeomFill_CorrectedFrenet` (`GeomFill_CorrectedFrenet.hxx:38-157`).
pub struct CorrectedFrenet {
    /// `myCurve` / `myTrimmed` (`GeomFill_TrihedronLaw.hxx:111-112`).
    pub(crate) data: TrihedronLawCurve,
    /// `frenet` (`hxx:145`).
    pub(crate) frenet: Frenet,
    /// `isFrenet`, `EvolAroundT`, `TLaw`, `AT`, `AN`, the OCC78 arrays
    /// (`hxx:146-155`).
    pub(crate) state: Mutex<CorrectedState>,
    /// `myForEvaluation` (`hxx:151`).
    for_evaluation: bool,
}

impl Default for CorrectedFrenet {
    fn default() -> Self {
        Self::new()
    }
}

impl CorrectedFrenet {
    /// `GeomFill_CorrectedFrenet` (`cxx:306-312`).
    pub fn new() -> Self {
        Self::new_for_evaluation(false)
    }

    /// `GeomFill_CorrectedFrenet(const bool ForEvaluation)`
    /// (`cxx:314-320`).
    pub fn new_for_evaluation(for_evaluation: bool) -> Self {
        Self {
            data: TrihedronLawCurve::new(),
            frenet: Frenet::new(),
            state: Mutex::new(CorrectedState::default()),
            for_evaluation,
        }
    }

    /// `myForEvaluation` (`hxx:151`).
    pub(crate) fn for_evaluation(&self) -> bool {
        self.for_evaluation
    }

    /// `myCurve` (`GeomFill_TrihedronLaw.hxx:111`).
    pub(crate) fn curve_arc(&self) -> Option<Arc<dyn Curve>> {
        self.data.curve()
    }

    /// `isFrenet` (`hxx:146`).
    fn is_frenet(&self) -> bool {
        self.state
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .is_frenet
    }

    fn set_is_frenet(&self, value: bool) {
        self.state
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .is_frenet = value;
    }

    /// `GeomFill_CorrectedFrenet::EvaluateBestMode` (`cxx:977-1035`).
    ///
    /// The `EvolAroundT` intervals and samples decide between Frenet,
    /// CorrectedFrenet and DiscreteTrihedron.
    pub fn evaluate_best_mode(&self) -> GeomFillTrihedron {
        let (evol, is_frenet) = {
            let state = self.state.lock().unwrap_or_else(|e| e.into_inner());
            (state.evol_around_t.clone(), state.is_frenet)
        };
        let Some(evol) = evol else {
            return GeomFillTrihedron::Frenet;
        };
        if is_frenet {
            return GeomFillTrihedron::Frenet;
        }

        const MAX_ANGLE: f64 = 3.0 * std::f64::consts::PI / 4.0;
        const MAX_TORSION: f64 = 100.0;
        const NB_SAMPLES: i32 = 10;

        let curve = match self.curve_arc() {
            Some(c) => c,
            None => return GeomFillTrihedron::Frenet,
        };
        let nb_int = evol.nb_intervals(occt_core::kernel::geomabs::Shape::CN).unwrap_or(1);
        let ints = match evol.intervals(occt_core::kernel::geomabs::Shape::CN) {
            Ok(v) => v,
            Err(_) => return GeomFillTrihedron::Frenet,
        };
        let mut prev = GpXyz::zero();
        let mut prev_vec = GpXyz::zero();
        let mut k = 1;
        for i in 1..=nb_int as usize {
            if i + 1 > ints.len() {
                break;
            }
            let tmin = ints[i - 1];
            let tmax = ints[i];
            let torsion = init::compute_torsion(tmin, curve.as_ref());
            if torsion.abs() > MAX_TORSION {
                return GeomFillTrihedron::DiscreteTrihedron;
            }
            let trimmed = match evol.trim(tmin, tmax, PCONFUSION / 2.0) {
                Ok(l) => l,
                Err(_) => return GeomFillTrihedron::Frenet,
            };
            let step = (tmax - tmin) / NB_SAMPLES as f64;
            for j in 0..=NB_SAMPLES {
                let u = tmin + j as f64 * step;
                let v = trimmed.value(u);
                let point = GpXyz::new(u, v, 0.0);
                if j != 0 {
                    let avec = GpXyz::new(
                        point.x - prev.x,
                        point.y - prev.y,
                        point.z - prev.z,
                    );
                    if k > 2 {
                        let ang = angle_2d(&prev_vec, &avec);
                        if ang.abs() > MAX_ANGLE {
                            return GeomFillTrihedron::DiscreteTrihedron;
                        }
                    }
                    prev_vec = avec;
                }
                prev = point;
                k += 1;
            }
        }
        GeomFillTrihedron::CorrectedFrenet
    }
}

/// `gp_Vec2d::Angle` for the two sample deltas of `EvaluateBestMode`
/// (`cxx:1019`).
fn angle_2d(a: &GpXyz, b: &GpXyz) -> f64 {
    let n1 = (a.x * a.x + a.y * a.y).sqrt();
    let n2 = (b.x * b.x + b.y * b.y).sqrt();
    if n1 <= occt_core::precision::RESOLUTION || n2 <= occt_core::precision::RESOLUTION {
        return 0.0;
    }
    let cos = ((a.x * b.x + a.y * b.y) / (n1 * n2)).clamp(-1.0, 1.0);
    cos.acos()
}

impl TrihedronLaw for CorrectedFrenet {
    /// `GeomFill_CorrectedFrenet::SetCurve` (`cxx:332-359`).
    fn set_curve(&self, curve: Arc<dyn Curve>) -> bool {
        self.data.set_curve(curve.clone());
        self.frenet.set_curve(curve.clone());
        let is_frenet = match extrema_curve_type(curve.as_ref()) {
            GeomAbsCurveType::Line
            | GeomAbsCurveType::Circle
            | GeomAbsCurveType::Ellipse
            | GeomAbsCurveType::Hyperbola
            | GeomAbsCurveType::Parabola => true,
            _ => {
                self.set_is_frenet(true);
                self.init();
                self.is_frenet()
            }
        };
        self.set_is_frenet(is_frenet);
        is_frenet
    }

    /// `GeomFill_CorrectedFrenet::Copy` (`cxx:322-330`).
    fn copy(&self) -> Arc<dyn TrihedronLaw> {
        let copy = CorrectedFrenet::new();
        if let Some(curve) = self.data.curve() {
            TrihedronLaw::set_curve(&copy, curve);
        }
        Arc::new(copy)
    }

    /// `GeomFill_CorrectedFrenet::D0` (`cxx:702-726`).
    fn d0(&self, param: f64) -> Result<(GpVec, GpVec, GpVec), String> {
        let (tangent, normal, _bi_normal) = self.frenet.d0(param)?;
        if self.is_frenet() {
            return Ok((tangent, normal, _bi_normal));
        }
        let angle_at = self.get_angle_at(param);
        let cross = tangent.crossed(&normal);
        let tcross = tangent.crossed(&cross);
        let normal = cross
            .multiplied_scalar(angle_at.sin())
            .added(&tcross.multiplied_scalar(1.0 - angle_at.cos()))
            .added(&normal);
        let bi_normal = tangent.crossed(&normal);
        Ok((tangent, normal, bi_normal))
    }

    /// `GeomFill_CorrectedFrenet::D1` (`cxx:732-785`).
    fn d1(&self, param: f64) -> Result<TrihedronD1, String> {
        let frame = self.frenet.d1(param)?;
        if self.is_frenet() {
            return Ok(frame);
        }
        let (_, d_angle_at) = match self.t_law_d1(param) {
            Some(v) => v,
            None => (0.0, 0.0),
        };
        let angle_at = self.get_angle_at(param);
        let (sina, cosa) = (angle_at.sin(), angle_at.cos());

        let mut normal = frame.normal;
        let mut d_normal = frame.dnormal;
        let cross = frame.tangent.crossed(&normal);
        let dcross = frame
            .dtangent
            .crossed(&normal)
            .added(&frame.tangent.crossed(&d_normal));
        let tcross = frame.tangent.crossed(&cross);
        let dtcross = frame
            .dtangent
            .crossed(&cross)
            .added(&frame.tangent.crossed(&dcross));

        let mut aux = dcross
            .multiplied_scalar(sina)
            .added(&cross.multiplied_scalar(cosa * d_angle_at));
        aux = dtcross
            .multiplied_scalar(1.0 - cosa)
            .added(&tcross.multiplied_scalar(sina * d_angle_at))
            .added(&aux);
        d_normal = d_normal.added(&aux);

        normal = cross
            .multiplied_scalar(sina)
            .added(&tcross.multiplied_scalar(1.0 - cosa))
            .added(&normal);
        let bi_normal = frame.tangent.crossed(&normal);
        let d_bi_normal = frame
            .dtangent
            .crossed(&normal)
            .added(&frame.tangent.crossed(&d_normal));
        Ok(TrihedronD1 {
            tangent: frame.tangent,
            dtangent: frame.dtangent,
            normal,
            dnormal: d_normal,
            binormal: bi_normal,
            dbinormal: d_bi_normal,
        })
    }

    /// `GeomFill_CorrectedFrenet::D2` (`cxx:787-900`).
    fn d2(&self, param: f64) -> Result<TrihedronD2, String> {
        let frame = self.frenet.d2(param)?;
        if self.is_frenet() {
            return Ok(frame);
        }
        let (_, d_angle_at, d2_angle_at) = match self.t_law_d2(param) {
            Some(v) => v,
            None => (0.0, 0.0, 0.0),
        };
        let angle_at = self.get_angle_at(param);
        let (sina, cosa) = (angle_at.sin(), angle_at.cos());

        let mut normal = frame.normal;
        let mut d_normal = frame.dnormal;
        let mut d2_normal = frame.d2normal;

        let cross = frame.tangent.crossed(&normal);
        let dcross = frame
            .dtangent
            .crossed(&normal)
            .added(&frame.tangent.crossed(&d_normal));
        let d2cross = frame
            .d2tangent
            .crossed(&normal)
            .added(&frame.dtangent.crossed(&d_normal).multiplied_scalar(2.0))
            .added(&frame.tangent.crossed(&d2_normal));

        let tcross = frame.tangent.crossed(&cross);
        let dtcross = frame
            .dtangent
            .crossed(&cross)
            .added(&frame.tangent.crossed(&dcross));
        let d2tcross = frame
            .d2tangent
            .crossed(&cross)
            .added(&frame.dtangent.crossed(&dcross).multiplied_scalar(2.0))
            .added(&frame.tangent.crossed(&d2cross));

        let mut aux = d2cross
            .multiplied_scalar(sina)
            .added(&dcross.multiplied_scalar(2.0 * cosa * d_angle_at))
            .added(&cross.multiplied_scalar(cosa * d2_angle_at - sina * d_angle_at * d_angle_at));
        aux = d2tcross
            .multiplied_scalar(1.0 - cosa)
            .added(&dtcross.multiplied_scalar(2.0 * sina * d_angle_at))
            .added(
                &tcross
                    .multiplied_scalar(cosa * d_angle_at * d_angle_at + sina * d2_angle_at),
            )
            .added(&aux);
        d2_normal = d2_normal.added(&aux);

        aux = dcross
            .multiplied_scalar(sina)
            .added(&cross.multiplied_scalar(cosa * d_angle_at));
        aux = dtcross
            .multiplied_scalar(1.0 - cosa)
            .added(&tcross.multiplied_scalar(sina * d_angle_at))
            .added(&aux);
        d_normal = d_normal.added(&aux);

        normal = cross
            .multiplied_scalar(sina)
            .added(&tcross.multiplied_scalar(1.0 - cosa))
            .added(&normal);
        let bi_normal = frame.tangent.crossed(&normal);
        let d_bi_normal = frame
            .dtangent
            .crossed(&normal)
            .added(&frame.tangent.crossed(&d_normal));
        let d2_bi_normal = frame
            .d2tangent
            .crossed(&normal)
            .added(&frame.dtangent.crossed(&d_normal).multiplied_scalar(2.0))
            .added(&frame.tangent.crossed(&d2_normal));
        Ok(TrihedronD2 {
            tangent: frame.tangent,
            dtangent: frame.dtangent,
            d2tangent: frame.d2tangent,
            normal,
            dnormal: d_normal,
            d2normal: d2_normal,
            binormal: bi_normal,
            dbinormal: d_bi_normal,
            d2binormal: d2_bi_normal,
        })
    }

    /// `GeomFill_CorrectedFrenet::NbIntervals` (`cxx:904-928`).
    fn nb_intervals(&self, s: u8) -> i32 {
        let nb_frenet = self.frenet.nb_intervals(s);
        if self.is_frenet() {
            return nb_frenet;
        }
        let t_law = self.t_law();
        let Some(t_law) = t_law else {
            return nb_frenet;
        };
        let shape = crate::geom_fill::snglr_func::shape_from_u8(s);
        let nb_law = t_law.nb_intervals(shape).unwrap_or(1);
        if nb_frenet == 1 {
            return nb_law;
        }
        let frenet_int = self.frenet.intervals(s);
        let law_int = match t_law.intervals(shape) {
            Ok(v) => v,
            Err(_) => return nb_frenet,
        };
        let fusion = crate::geom_lib::fuse_intervals(&frenet_int, &law_int, PCONFUSION, true);
        fusion.len() as i32 - 1
    }

    /// `GeomFill_CorrectedFrenet::Intervals` (`cxx:932-961`).
    fn intervals(&self, s: u8) -> Vec<f64> {
        if self.is_frenet() {
            return self.frenet.intervals(s);
        }
        let nb_frenet = self.frenet.nb_intervals(s);
        let t_law = self.t_law();
        let Some(t_law) = t_law else {
            return self.frenet.intervals(s);
        };
        let shape = crate::geom_fill::snglr_func::shape_from_u8(s);
        let frenet_int = self.frenet.intervals(s);
        let law_int = match t_law.intervals(shape) {
            Ok(v) => v,
            Err(_) => return frenet_int,
        };
        let _ = nb_frenet;
        crate::geom_lib::fuse_intervals(&frenet_int, &law_int, PCONFUSION, true)
    }

    /// `GeomFill_CorrectedFrenet::SetInterval` (`cxx:965-973`).
    fn set_interval(&self, first: f64, last: f64) -> Result<(), String> {
        self.data.set_interval(first, last)?;
        self.frenet.set_interval(first, last)?;
        if !self.is_frenet() {
            let evol = {
                let state = self.state.lock().unwrap_or_else(|e| e.into_inner());
                state.evol_around_t.clone()
            };
            if let Some(evol) = evol {
                if let Ok(trimmed) = (evol as Arc<dyn LawFunction>).trim(first, last, PCONFUSION / 2.0)
                {
                    self.state
                        .lock()
                        .unwrap_or_else(|e| e.into_inner())
                        .t_law = Some(trimmed);
                }
            }
        }
        Ok(())
    }

    /// `GeomFill_TrihedronLaw::GetInterval` (`cxx:62-67`).
    fn get_interval(&self) -> (f64, f64) {
        self.data.interval()
    }

    /// `GeomFill_CorrectedFrenet::GetAverageLaw` (`cxx:1037-1050`).
    fn get_average_law(&self) -> (GpVec, GpVec, GpVec) {
        if self.is_frenet() {
            return self.frenet.get_average_law();
        }
        let state = self.state.lock().unwrap_or_else(|e| e.into_inner());
        let at = state.at;
        let an = state.an;
        let bi_normal = at.crossed(&an);
        (at, an, bi_normal)
    }

    /// `GeomFill_CorrectedFrenet::IsConstant` (`cxx:1052-1057`).
    fn is_constant(&self) -> bool {
        match self.data.curve() {
            Some(c) => extrema_curve_type(c.as_ref()) == GeomAbsCurveType::Line,
            None => false,
        }
    }

    /// `GeomFill_CorrectedFrenet::IsOnlyBy3dCurve` (`cxx:1061-1064`).
    fn is_only_by_3d_curve(&self) -> bool {
        true
    }
}

impl CorrectedFrenet {
    /// `TLaw` (`hxx:148`).
    fn t_law(&self) -> Option<Arc<dyn LawFunction>> {
        self.state
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .t_law
            .clone()
    }

    /// `TLaw->D1` (`cxx:747`).
    fn t_law_d1(&self, param: f64) -> Option<(f64, f64)> {
        let law = self.t_law()?;
        Some(law.d1(param))
    }

    /// `TLaw->D2` (`cxx:815`).
    fn t_law_d2(&self, param: f64) -> Option<(f64, f64, f64)> {
        let law = self.t_law()?;
        Some(law.d2(param))
    }
}
