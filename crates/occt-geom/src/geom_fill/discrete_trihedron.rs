//! `GeomFill_DiscreteTrihedron` — the trihedron built by discrete transport:
//! the path is sampled, the first frame comes from `gp_Ax2(Origin, Tangent)`,
//! and every following frame is the previous one rotated by the angle between
//! the two tangents.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_DiscreteTrihedron.hxx` and
//! `GeomFill_DiscreteTrihedron.cxx:32-390` (OCCT 8.0.0): constructor
//! `cxx:36-41`, `Copy` `cxx:46-54`, `SetCurve` `cxx:58-85`, `Init`
//! `cxx:90-173`, `D0` `cxx:177-267`, `D1` `cxx:273-294`, `D2`
//! `cxx:298-333`, `NbIntervals` `cxx:337-340`, `Intervals` `cxx:344-348`,
//! `GetAverageLaw` `cxx:350-377`, `IsConstant` `cxx:381-384`,
//! `IsOnlyBy3dCurve` `cxx:387-390`.
//!
//! Translation notes:
//! - `myUseFrenet` / `myKnots` / `myTrihedrons` / `myPoint` (`hxx:110-114`)
//!   are written by the non-const `SetCurve` / `Init` through a shared handle;
//!   the port takes `&self` in those methods (see
//!   [`crate::geom_fill::trihedron_law`]), so they live behind one `Mutex`.
//! - `gp_Ax2` is [`GpAx2`]; `gp_Ax2(Origin, TangDir)` is
//!   [`GpAx2::from_axis`], `Rotated` / `SetDirection` keep their OCCT names.
//! - `gp_Vec::AngleWithRef` is
//!   [`crate::extrema_surf::elementary_curve_extrema::vec_angle_with_ref`].
//! - `gp_Dir(Tangent)` raising `Standard_ConstructionError` on a null vector
//!   becomes `GpDir::from_vec`'s `Err`.
//! - The `InsertBefore(i, NewParam)` + `i--` re-walk of `Init` (`cxx:159-160`)
//!   is the `while` loop below.

use std::sync::{Arc, Mutex};

use occt_core::gp::{GpAx1, GpAx2, GpDir, GpPnt, GpVec};
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::{CONFUSION, PCONFUSION, RESOLUTION};

use crate::curve::Curve;
use crate::extrema_pc::{extrema_curve_type, GeomAbsCurveType};
use crate::extrema_surf::vec_angle_with_ref;
use crate::geom_fill::frenet::Frenet;
use crate::geom_fill::trihedron_law::{TrihedronD1, TrihedronD2, TrihedronLaw, TrihedronLawCurve};

/// `TolConf` (`cxx:32`).
const TOL_CONF: f64 = CONFUSION;
/// `NbSamples` (`cxx:97`, `cxx:191`).
const NB_SAMPLES: usize = 10;

/// The mutable OCCT members of `GeomFill_DiscreteTrihedron`
/// (`GeomFill_DiscreteTrihedron.hxx:110-114`).
struct DiscreteState {
    /// `myUseFrenet` (`hxx:114`).
    use_frenet: bool,
    /// `myKnots` (`hxx:112`).
    knots: Vec<f64>,
    /// `myTrihedrons` (`hxx:111`).
    trihedrons: Vec<GpAx2>,
    /// `myPoint` (`hxx:110`): scratch of `D0`.
    point: GpPnt,
}

/// `GeomFill_DiscreteTrihedron` (`GeomFill_DiscreteTrihedron.hxx:36-116`).
pub struct DiscreteTrihedron {
    /// `myCurve` / `myTrimmed` (`GeomFill_TrihedronLaw.hxx:111-112`).
    data: TrihedronLawCurve,
    /// `myFrenet` (`hxx:113`).
    frenet: Frenet,
    /// `myUseFrenet` / `myKnots` / `myTrihedrons` / `myPoint`.
    state: Mutex<DiscreteState>,
}

impl Default for DiscreteTrihedron {
    fn default() -> Self {
        Self::new()
    }
}

impl DiscreteTrihedron {
    /// `GeomFill_DiscreteTrihedron()` (`cxx:36-41`): `myUseFrenet = false`,
    /// `myFrenet = new GeomFill_Frenet`, the two sequences empty.
    pub fn new() -> Self {
        Self {
            data: TrihedronLawCurve::new(),
            frenet: Frenet::new(),
            state: Mutex::new(DiscreteState {
                use_frenet: false,
                knots: Vec::new(),
                trihedrons: Vec::new(),
                point: GpPnt::new(0.0, 0.0, 0.0),
            }),
        }
    }

    /// `myUseFrenet` (`hxx:114`).
    pub fn use_frenet(&self) -> bool {
        self.state.lock().unwrap_or_else(|e| e.into_inner()).use_frenet
    }

    /// `myTrimmed` (`hxx`, inherited from the base).
    fn trimmed(&self) -> Result<Arc<dyn Curve>, String> {
        self.data.curve_or_err("GeomFill_DiscreteTrihedron")
    }

    /// `GeomFill_DiscreteTrihedron::Init` (`cxx:90-173`): the sample knots and
    /// the transported frames.
    fn init(&self) -> Result<(), String> {
        let curve = self.trimmed()?;
        let nb_intervals = curve.nb_intervals(Shape::CN as u8);
        let bounds = curve.parameter_intervals(Shape::CN as u8);
        if nb_intervals <= 0 || bounds.len() < 2 {
            return Err("GeomFill_DiscreteTrihedron::Init: no CN interval".into());
        }

        let mut knots: Vec<f64> = Vec::new();
        for i in 1..=nb_intervals as usize {
            let delta = (bounds[i] - bounds[i - 1]) / NB_SAMPLES as f64;
            for j in 0..NB_SAMPLES {
                knots.push(bounds[i - 1] + j as f64 * delta);
            }
        }
        knots.push(bounds[nb_intervals as usize]);

        let origin = GpPnt::new(0.0, 0.0, 0.0);
        let mut trihedrons: Vec<GpAx2> = Vec::new();
        let mut i = 1usize;
        while i <= knots.len() {
            let param = knots[i - 1];
            let (pnt, tangent) = curve.d1(param);
            let mut tangent = tangent;
            let norm = tangent.magnitude();
            if norm < TOL_CONF {
                let sub_delta = (knots[i] - knots[i - 1]) / NB_SAMPLES as f64;
                let sub_delta = if sub_delta < PCONFUSION {
                    knots[i] - knots[i - 1]
                } else {
                    sub_delta
                };
                let sub_pnt = curve.d0(param + sub_delta);
                tangent = GpVec::from_pnts(&pnt, &sub_pnt);
            }
            // `TangDir = Tangent; Tangent = TangDir;` (`cxx:131-132`).
            let tang_dir = GpDir::from_vec(&tangent)
                .map_err(|_| "Standard_ConstructionError: gp_Dir(Tangent)".to_string())?;
            let tangent = GpVec::from_xyz(tang_dir.xyz());

            if i == 1 {
                trihedrons.push(GpAx2::from_axis(origin, tang_dir));
            } else {
                let last_axis = *trihedrons
                    .last()
                    .ok_or_else(|| "GeomFill_DiscreteTrihedron::Init: empty".to_string())?;
                let last_tangent = GpVec::from_xyz(last_axis.direction().xyz());
                let axis_of_rotation = last_tangent.crossed(&tangent);
                if axis_of_rotation.magnitude() <= RESOLUTION {
                    let scalar_product = last_tangent.dot(&tangent);
                    if scalar_product > 0.0 {
                        trihedrons.push(last_axis);
                    } else {
                        let new_param = (knots[i - 2] + knots[i - 1]) / 2.0;
                        if new_param - knots[i - 2] < RESOLUTION {
                            return Err(
                                "Standard_ConstructionError: GeomFill_DiscreteTrihedron : \
                                 impassable singularities on path curve"
                                    .into(),
                            );
                        }
                        knots.insert(i - 1, new_param);
                        i -= 1;
                    }
                } else {
                    let the_angle = vec_angle_with_ref(&last_tangent, &tangent, &axis_of_rotation);
                    let dir = GpDir::from_vec(&axis_of_rotation)
                        .map_err(|_| "Standard_ConstructionError: gp_Dir(Axis)".to_string())?;
                    let the_axis_of_rotation = GpAx1::new(origin, dir);
                    let mut new_axis = last_axis.rotated(&the_axis_of_rotation, the_angle);
                    new_axis.set_direction(tang_dir);
                    trihedrons.push(new_axis);
                }
            }
            i += 1;
        }

        let mut state = self.state.lock().unwrap_or_else(|e| e.into_inner());
        state.knots = knots;
        state.trihedrons = trihedrons;
        Ok(())
    }
}

impl TrihedronLaw for DiscreteTrihedron {
    /// `GeomFill_DiscreteTrihedron::SetCurve` (`cxx:58-85`).
    fn set_curve(&self, curve: Arc<dyn Curve>) -> bool {
        self.data.set_curve(curve.clone());
        match extrema_curve_type(curve.as_ref()) {
            GeomAbsCurveType::Circle
            | GeomAbsCurveType::Ellipse
            | GeomAbsCurveType::Hyperbola
            | GeomAbsCurveType::Parabola
            | GeomAbsCurveType::Line => {
                self.frenet.set_curve(curve);
                self.state.lock().unwrap_or_else(|e| e.into_inner()).use_frenet = true;
            }
            _ => {
                self.state.lock().unwrap_or_else(|e| e.into_inner()).use_frenet = false;
                let _ = self.init();
            }
        }
        self.use_frenet()
    }

    /// `GeomFill_DiscreteTrihedron::Copy` (`cxx:46-54`).
    fn copy(&self) -> Arc<dyn TrihedronLaw> {
        let opy = DiscreteTrihedron::new();
        if let Some(curve) = self.data.curve() {
            TrihedronLaw::set_curve(&opy, curve);
        }
        Arc::new(opy)
    }

    /// `GeomFill_DiscreteTrihedron::D0` (`cxx:177-267`).
    fn d0(&self, param: f64) -> Result<(GpVec, GpVec, GpVec), String> {
        let use_frenet = self.use_frenet();
        if use_frenet {
            return self.frenet.d0(param);
        }
        let curve = self.trimmed()?;
        let origin = GpPnt::new(0.0, 0.0, 0.0);

        let mut state = self.state.lock().unwrap_or_else(|e| e.into_inner());
        let tol_par = PCONFUSION;
        let n = state.knots.len();
        if n == 0 {
            return Err("GeomFill_DiscreteTrihedron::D0: Init not run".into());
        }

        // Binary search for <Param> in <myKnots> (`cxx:200-221`).
        let mut i1 = 1usize;
        let mut i2 = n;
        let mut index;
        loop {
            let i = (i1 + i2) / 2;
            if param <= state.knots[i - 1] {
                i2 = i;
            } else {
                i1 = i;
            }
            if i2 - i1 <= 1 {
                break;
            }
        }
        index = i1;
        if (param - state.knots[i2 - 1]).abs() < tol_par {
            index = i2;
        }

        let prev_param = state.knots[index - 1];
        let prev_axis = state.trihedrons[index - 1];
        let the_axis;
        if (param - prev_param).abs() < tol_par {
            the_axis = prev_axis;
        } else {
            let (pnt, mut tangent) = curve.d1(param);
            state.point = pnt;
            let norm = tangent.magnitude();
            if norm < TOL_CONF {
                let sub_delta = (state.knots[index] - param) / NB_SAMPLES as f64;
                let sub_delta = if sub_delta < PCONFUSION {
                    state.knots[index] - param
                } else {
                    sub_delta
                };
                let sub_pnt = curve.d0(param + sub_delta);
                tangent = GpVec::from_pnts(&state.point, &sub_pnt);
            }
            let tang_dir = GpDir::from_vec(&tangent)
                .map_err(|_| "Standard_ConstructionError: gp_Dir(Tangent)".to_string())?;
            let tangent = GpVec::from_xyz(tang_dir.xyz());
            let prev_tangent = GpVec::from_xyz(prev_axis.direction().xyz());
            let axis_of_rotation = prev_tangent.crossed(&tangent);
            if axis_of_rotation.magnitude() <= RESOLUTION {
                // "we assume that tangents can not be opposite" (`cxx:250-251`).
                the_axis = prev_axis;
            } else {
                let the_angle = vec_angle_with_ref(&prev_tangent, &tangent, &axis_of_rotation);
                let dir = GpDir::from_vec(&axis_of_rotation)
                    .map_err(|_| "Standard_ConstructionError: gp_Dir(Axis)".to_string())?;
                let the_axis_of_rotation = GpAx1::new(origin, dir);
                the_axis = prev_axis.rotated(&the_axis_of_rotation, the_angle);
            }
            let mut the_axis = the_axis;
            the_axis.set_direction(tang_dir);
            return Ok(frame(&the_axis));
        }
        Ok(frame(&the_axis))
    }

    /// `GeomFill_DiscreteTrihedron::D1` (`cxx:273-294`).
    fn d1(&self, param: f64) -> Result<TrihedronD1, String> {
        if self.use_frenet() {
            return self.frenet.d1(param);
        }
        let (tangent, normal, binormal) = self.d0(param)?;
        Ok(TrihedronD1 {
            tangent,
            dtangent: GpVec::zero(),
            normal,
            dnormal: GpVec::zero(),
            binormal,
            dbinormal: GpVec::zero(),
        })
    }

    /// `GeomFill_DiscreteTrihedron::D2` (`cxx:298-333`).
    fn d2(&self, param: f64) -> Result<TrihedronD2, String> {
        if self.use_frenet() {
            return self.frenet.d2(param);
        }
        let (tangent, normal, binormal) = self.d0(param)?;
        Ok(TrihedronD2 {
            tangent,
            dtangent: GpVec::zero(),
            d2tangent: GpVec::zero(),
            normal,
            dnormal: GpVec::zero(),
            d2normal: GpVec::zero(),
            binormal,
            dbinormal: GpVec::zero(),
            d2binormal: GpVec::zero(),
        })
    }

    /// `GeomFill_DiscreteTrihedron::NbIntervals` (`cxx:337-340`): always the
    /// `GeomAbs_CN` count of `myTrimmed`, whatever `S` is.
    fn nb_intervals(&self, _s: u8) -> i32 {
        match self.data.curve() {
            Some(c) => c.nb_intervals(Shape::CN as u8),
            None => 1,
        }
    }

    /// `GeomFill_DiscreteTrihedron::Intervals` (`cxx:344-348`).
    fn intervals(&self, _s: u8) -> Vec<f64> {
        match self.data.curve() {
            Some(c) => c.parameter_intervals(Shape::CN as u8),
            None => vec![0.0, 0.0],
        }
    }

    /// `GeomFill_TrihedronLaw::SetInterval` (`GeomFill_TrihedronLaw.cxx:57-60`).
    fn set_interval(&self, first: f64, last: f64) -> Result<(), String> {
        self.data.set_interval(first, last)
    }

    /// `GeomFill_TrihedronLaw::GetInterval` (`GeomFill_TrihedronLaw.cxx:62-67`).
    fn get_interval(&self) -> (f64, f64) {
        self.data.interval()
    }

    /// `GeomFill_DiscreteTrihedron::GetAverageLaw` (`cxx:350-377`).
    fn get_average_law(&self) -> (GpVec, GpVec, GpVec) {
        const NUM: usize = 20;
        let mut a_tangent = GpVec::zero();
        let mut a_normal = GpVec::zero();
        let mut a_bi_normal = GpVec::zero();
        let (first, last) = self.data.interval();
        let step = (last - first) / NUM as f64;
        for i in 0..=NUM {
            let mut param = first + i as f64 * step;
            if param > last {
                param = last;
            }
            if let Ok((t, n, bn)) = self.d0(param) {
                a_tangent = a_tangent.added(&t);
                a_normal = a_normal.added(&n);
                a_bi_normal = a_bi_normal.added(&bn);
            }
        }
        a_tangent = a_tangent.divided((NUM + 1) as f64);
        a_normal = a_normal.divided((NUM + 1) as f64);
        a_tangent = a_tangent.normalized();
        a_bi_normal = a_tangent.crossed(&a_normal).normalized();
        a_normal = a_bi_normal.crossed(&a_tangent);
        (a_tangent, a_normal, a_bi_normal)
    }

    /// `GeomFill_DiscreteTrihedron::IsConstant` (`cxx:381-384`).
    fn is_constant(&self) -> bool {
        match self.data.curve() {
            Some(c) => matches!(extrema_curve_type(c.as_ref()), GeomAbsCurveType::Line),
            None => false,
        }
    }

    /// `GeomFill_DiscreteTrihedron::IsOnlyBy3dCurve` (`cxx:387-390`): true.
    fn is_only_by_3d_curve(&self) -> bool {
        true
    }
}

/// `Tangent = theAxis.Direction(); Normal = theAxis.XDirection();
/// BiNormal = theAxis.YDirection();` (`cxx:264-266`).
fn frame(axis: &GpAx2) -> (GpVec, GpVec, GpVec) {
    (
        GpVec::from_xyz(axis.direction().xyz()),
        GpVec::from_xyz(axis.x_direction().xyz()),
        GpVec::from_xyz(axis.y_direction().xyz()),
    )
}
