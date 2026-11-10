//! `GeomFill_CurveAndTrihedron` — the `GeomFill_LocationLaw` made of a spine
//! curve plus a `GeomFill_TrihedronLaw`.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_CurveAndTrihedron.hxx` and
//! `GeomFill_CurveAndTrihedron.cxx:34-342` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `myCurve` / `myTrimmed` (`hxx:70-71`) follow the same handle semantics as
//!   [`crate::geom_fill::trihedron_law::TrihedronLawCurve`]: `SetCurve` stores
//!   both handles, `SetInterval` calls the law's `SetInterval` and replaces
//!   `myTrimmed` with the trimmed adaptor, and every evaluation goes through
//!   `myTrimmed`. The port stores `myCurve` and the trimmed range behind a
//!   `Mutex`, so the non-const OCCT methods can take `&self`.
//! - `M.SetCols(V2, V3, V1)` is [`GpMat::set_cols`] with the (normal, binormal,
//!   tangent) column order of `cxx:100`; `M *= Trans` is
//!   [`GpMat::multiply`].
//! - `GeomLib::FuseIntervals` (`cxx:216`, `cxx:245`) is
//!   [`crate::geom_lib::fuse_intervals`] with
//!   `Precision::PConfusion() * 0.99`.
//! - `myCurve->GetType() == GeomAbs_Line` / `GeomAbs_Circle` becomes
//!   [`Curve::is_line`] / [`Curve::gp_circ`].

use std::sync::{Arc, Mutex};

use occt_core::gp::{GpMat, GpPnt2d, GpVec};
use occt_core::precision::PCONFUSION;

use crate::curve::Curve;
use crate::geom_fill::location_law::{LocationD1, LocationD2, LocationLaw};
use crate::geom_fill::trihedron_law::TrihedronLaw;
use crate::geom_lib::fuse_intervals;

/// The `myCurve` / `myTrimmed` / `Trans` / `WithTrans` state of
/// `GeomFill_CurveAndTrihedron` (`hxx:70-73`).
struct State {
    /// `myCurve` (`hxx:70`).
    curve: Option<Arc<dyn Curve>>,
    /// `myTrimmed` (`hxx:71`): the trimmed parameter range.
    trimmed: Option<(f64, f64)>,
    /// `Trans` (`hxx:72`).
    trans: GpMat,
    /// `WithTrans` (`hxx:73`).
    with_trans: bool,
}

/// `GeomFill_CurveAndTrihedron` (`GeomFill_CurveAndTrihedron.hxx:38-90`).
pub struct CurveAndTrihedron {
    /// `myLaw` (`hxx:68`).
    law: Arc<dyn TrihedronLaw>,
    state: Mutex<State>,
}

impl CurveAndTrihedron {
    /// `GeomFill_CurveAndTrihedron::GeomFill_CurveAndTrihedron` (`cxx:35-42`).
    pub fn new(trihedron: Arc<dyn TrihedronLaw>) -> Self {
        Self {
            law: trihedron,
            state: Mutex::new(State {
                curve: None,
                trimmed: None,
                trans: GpMat::identity(),
                with_trans: false,
            }),
        }
    }

    /// `myLaw` (`hxx:68`).
    pub fn law(&self) -> &Arc<dyn TrihedronLaw> {
        &self.law
    }

    /// The current trimmed range (`myTrimmed`), used by the laws that combine
    /// location intervals.
    fn trimmed_range(&self) -> (f64, f64) {
        self.state
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .trimmed
            .unwrap_or((0.0, 0.0))
    }

    fn curve_arc(&self) -> Option<Arc<dyn Curve>> {
        self.state
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .curve
            .clone()
    }
}

impl LocationLaw for CurveAndTrihedron {
    /// `GeomFill_CurveAndTrihedron::SetCurve` (`cxx:58-62`).
    fn set_curve(&self, curve: Arc<dyn Curve>) -> bool {
        let ok = self.law.set_curve(Arc::clone(&curve));
        let mut state = self.state.lock().unwrap_or_else(|e| e.into_inner());
        state.trimmed = Some((curve.first_parameter(), curve.last_parameter()));
        state.curve = Some(curve);
        ok
    }

    /// `GeomFill_CurveAndTrihedron::GetCurve` (`cxx:64-67`).
    fn get_curve(&self) -> Option<Arc<dyn Curve>> {
        self.curve_arc()
    }

    /// `GeomFill_CurveAndTrihedron::SetTrsf` (`cxx:71-89`).
    fn set_trsf(&self, transfo: &GpMat) {
        let mut state = self.state.lock().unwrap_or_else(|e| e.into_inner());
        state.trans = *transfo;
        let aux = GpMat::identity().subtract(transfo);
        let mut with_trans = false;
        for ii in 1..=3 {
            for jj in 1..=3 {
                if aux.value(ii, jj).abs() > 1e-14 {
                    with_trans = true;
                }
            }
            if with_trans {
                break;
            }
        }
        state.with_trans = with_trans;
    }

    /// `GeomFill_CurveAndTrihedron::Copy` (`cxx:46-55`).
    fn copy(&self) -> Arc<dyn LocationLaw> {
        let law = self.law.copy();
        let copy = CurveAndTrihedron::new(law);
        if let Some(curve) = self.curve_arc() {
            copy.set_curve(curve);
        }
        let state = self.state.lock().unwrap_or_else(|e| e.into_inner());
        copy.set_trsf(&state.trans);
        drop(state);
        Arc::new(copy)
    }

    /// `GeomFill_CurveAndTrihedron::D0` (`cxx:93-107`).
    fn d0(&self, param: f64) -> Result<(GpMat, GpVec), String> {
        let curve = self
            .curve_arc()
            .ok_or_else(|| "Standard_NullObject: GeomFill_CurveAndTrihedron::D0".to_string())?;
        let point = curve.d0(param);
        let v = GpVec::new(point.x(), point.y(), point.z());
        let (v1, v2, v3) = self.law.d0(param)?;
        let mut m = GpMat::from_cols(v2.xyz(), v3.xyz(), v1.xyz());
        let state = self.state.lock().unwrap_or_else(|e| e.into_inner());
        if state.with_trans {
            m = m.multiply(&state.trans);
        }
        drop(state);
        Ok((m, v))
    }

    /// `GeomFill_CurveAndTrihedron::D0(Param, M, V, Poles2d)` (`cxx:111-126`).
    fn d0_with_poles(&self, param: f64) -> Result<(GpMat, GpVec, Vec<GpPnt2d>), String> {
        let (m, v) = self.d0(param)?;
        Ok((m, v, Vec::new()))
    }

    /// `GeomFill_CurveAndTrihedron::D1` (`cxx:132-155`).
    fn d1(&self, param: f64) -> Result<LocationD1, String> {
        let curve = self
            .curve_arc()
            .ok_or_else(|| "Standard_NullObject: GeomFill_CurveAndTrihedron::D1".to_string())?;
        let (point, dv) = curve.d1(param);
        let v = GpVec::new(point.x(), point.y(), point.z());
        let frame = self.law.d1(param)?;
        let mut m = GpMat::from_cols(frame.normal.xyz(), frame.binormal.xyz(), frame.tangent.xyz());
        let mut dm = GpMat::from_cols(
            frame.dnormal.xyz(),
            frame.dbinormal.xyz(),
            frame.dtangent.xyz(),
        );
        let state = self.state.lock().unwrap_or_else(|e| e.into_inner());
        if state.with_trans {
            m = m.multiply(&state.trans);
            dm = dm.multiply(&state.trans);
        }
        drop(state);
        Ok(LocationD1 {
            m,
            v,
            dm,
            dv,
            poles2d: Vec::new(),
            dpoles2d: Vec::new(),
        })
    }

    /// `GeomFill_CurveAndTrihedron::D2` (`cxx:161-191`).
    fn d2(&self, param: f64) -> Result<LocationD2, String> {
        let curve = self
            .curve_arc()
            .ok_or_else(|| "Standard_NullObject: GeomFill_CurveAndTrihedron::D2".to_string())?;
        let (_point, dv, d2v) = curve.d2(param);
        let frame = self.law.d2(param)?;
        let mut m = GpMat::from_cols(frame.normal.xyz(), frame.binormal.xyz(), frame.tangent.xyz());
        let mut dm = GpMat::from_cols(
            frame.dnormal.xyz(),
            frame.dbinormal.xyz(),
            frame.dtangent.xyz(),
        );
        let mut d2m = GpMat::from_cols(
            frame.d2normal.xyz(),
            frame.d2binormal.xyz(),
            frame.d2tangent.xyz(),
        );
        let state = self.state.lock().unwrap_or_else(|e| e.into_inner());
        if state.with_trans {
            m = m.multiply(&state.trans);
            dm = dm.multiply(&state.trans);
            d2m = d2m.multiply(&state.trans);
        }
        drop(state);
        Ok(LocationD2 {
            m,
            v: GpVec::from_xyz(&occt_core::gp::GpXyz::new(_point.x(), _point.y(), _point.z())),
            dm,
            dv,
            d2m,
            d2v,
            poles2d: Vec::new(),
            dpoles2d: Vec::new(),
            d2poles2d: Vec::new(),
        })
    }

    /// `GeomFill_CurveAndTrihedron::NbIntervals` (`cxx:195-218`).
    fn nb_intervals(&self, s: u8) -> i32 {
        let curve = match self.curve_arc() {
            Some(c) => c,
            None => return self.law.nb_intervals(s),
        };
        let nb_sec = curve.nb_intervals(s);
        let nb_law = self.law.nb_intervals(s);
        if nb_sec == 1 {
            return nb_law;
        }
        if nb_law == 1 {
            return nb_sec;
        }
        let int_c = curve.parameter_intervals(s);
        let int_l = self.law.intervals(s);
        let inter = fuse_intervals(&int_c, &int_l, PCONFUSION * 0.99, false);
        inter.len() as i32 - 1
    }

    /// `GeomFill_CurveAndTrihedron::Intervals` (`cxx:222-251`).
    fn intervals(&self, s: u8) -> Vec<f64> {
        let curve = match self.curve_arc() {
            Some(c) => c,
            None => return self.law.intervals(s),
        };
        let nb_sec = curve.nb_intervals(s);
        let nb_law = self.law.nb_intervals(s);
        if nb_sec == 1 {
            return self.law.intervals(s);
        }
        if nb_law == 1 {
            return curve.parameter_intervals(s);
        }
        let int_c = curve.parameter_intervals(s);
        let int_l = self.law.intervals(s);
        fuse_intervals(&int_c, &int_l, PCONFUSION * 0.99, false)
    }

    /// `GeomFill_CurveAndTrihedron::SetInterval` (`cxx:255-259`).
    fn set_interval(&self, first: f64, last: f64) {
        let _ = self.law.set_interval(first, last);
        let mut state = self.state.lock().unwrap_or_else(|e| e.into_inner());
        state.trimmed = Some((first, last));
    }

    /// `GeomFill_CurveAndTrihedron::GetInterval` (`cxx:263-267`).
    fn get_interval(&self) -> (f64, f64) {
        self.trimmed_range()
    }

    /// `GeomFill_CurveAndTrihedron::GetDomain` (`cxx:271-275`).
    fn get_domain(&self) -> (f64, f64) {
        match self.curve_arc() {
            Some(c) => (c.first_parameter(), c.last_parameter()),
            None => (0.0, 0.0),
        }
    }

    /// `GeomFill_CurveAndTrihedron::GetMaximalNorm` (`cxx:283-286`): 1.
    fn get_maximal_norm(&self) -> f64 {
        1.0
    }

    /// `GeomFill_CurveAndTrihedron::GetAverageLaw` (`cxx:290-308`).
    fn get_average_law(&self) -> (GpMat, GpVec) {
        let (v1, v2, v3) = self.law.get_average_law();
        let am = GpMat::from_cols(v1.xyz(), v2.xyz(), v3.xyz());
        let curve = match self.curve_arc() {
            Some(c) => c,
            None => return (am, GpVec::zero()),
        };
        let (first, last) = self.get_interval();
        let delta = (last - first) / 10.0;
        let mut av = GpVec::zero();
        let mut u = first;
        for _ in 0..=10 {
            let p = curve.d0(u);
            av = av.added(&GpVec::new(p.x(), p.y(), p.z()));
            u += delta;
        }
        (am, av.divided(11.0))
    }

    /// `GeomFill_CurveAndTrihedron::IsTranslation` (`cxx:312-321`).
    fn is_translation(&self, error: &mut f64) -> bool {
        *error = 0.0;
        match self.curve_arc() {
            Some(curve) if curve.is_line() => {
                self.law.is_constant() || self.law.is_only_by_3d_curve()
            }
            _ => false,
        }
    }

    /// `GeomFill_CurveAndTrihedron::IsRotation` (`cxx:325-335`).
    fn is_rotation(&self, error: &mut f64) -> bool {
        *error = 0.0;
        match self.curve_arc() {
            Some(curve) if curve.gp_circ().is_some() => self.law.is_only_by_3d_curve(),
            _ => false,
        }
    }

    /// `GeomFill_CurveAndTrihedron::Rotation` (`cxx:339-342`):
    /// `myCurve->Circle().Location()`.
    fn rotation(&self) -> Result<occt_core::gp::GpPnt, String> {
        match self.curve_arc() {
            Some(curve) => match curve.gp_circ() {
                Some(circ) => Ok(circ.location()),
                None => Err("Standard_DomainError: GeomFill_CurveAndTrihedron::Rotation".into()),
            },
            None => Err("Standard_NullObject: GeomFill_CurveAndTrihedron::Rotation".into()),
        }
    }
}
