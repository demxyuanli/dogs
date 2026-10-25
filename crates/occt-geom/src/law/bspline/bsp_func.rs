//! `Law_BSpFunc` — a `Law_BSpline` restricted to a parametric window.
//!
//! Source: `TKGeomAlgo/Law/Law_BSpFunc.hxx`/`Law_BSpFunc.cxx:33-406`
//! (OCCT 8.0.0).
//!
//! `Law_BSpFunc::Value` / `D1` / `D2` / `Bounds` / `Trim` are the
//! `Law_Function` virtuals and live in the `impl LawFunction` block at the end
//! of this file (`Law_Function.hxx:35-77`).

use std::sync::Arc;

use occt_core::bspl::locate::locate_parameter_knots_range;
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::PCONFUSION;

use super::LawBSpline;
use crate::law::function::LawFunction;
use crate::law::knot_split::LawBSplineKnotSplitting;

/// `LocalValue`/`LocalD1`/`LocalD2` so that the boundary knots stay exact.

#[derive(Debug, Clone, Default)]
pub struct LawBSpFunc {
    curv: Option<LawBSpline>,
    first: f64,
    last: f64,
}

impl LawBSpFunc {
    /// `Law_BSpFunc::Law_BSpFunc` (`Law_BSpFunc.cxx:35-38`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `Law_BSpFunc::Law_BSpFunc(C, First, Last)` (`cxx:44-48`).
    pub fn new_with_curve(curve: LawBSpline, first: f64, last: f64) -> Self {
        Self {
            curv: Some(curve),
            first,
            last,
        }
    }

    fn curve_ref(&self) -> &LawBSpline {
        // OCCT dereferences an unset handle here; the port panics with a
        // diagnosable message instead of reading through a null pointer.
        self.curv
            .as_ref()
            .expect("Law_BSpFunc: no curve set (Law_BSpFunc.cxx:56-406)")
    }

    /// `Law_BSpFunc::SetCurve` (`Law_BSpFunc.cxx:396-403`).
    pub fn set_curve(&mut self, curve: LawBSpline) {
        self.first = curve.first_parameter();
        self.last = curve.last_parameter();
        self.curv = Some(curve);
    }

    /// `Law_BSpFunc::Curve` (`Law_BSpFunc.cxx:387-392`).
    pub fn curve(&self) -> Option<&LawBSpline> {
        self.curv.as_ref()
    }

    /// `Law_BSpFunc::Curve` as a mutable handle accessor.
    pub fn curve_mut(&mut self) -> Option<&mut LawBSpline> {
        self.curv.as_mut()
    }

    /// `Law_BSpFunc::Bounds` (`Law_BSpFunc.cxx:383-387`).
    pub fn bounds(&self) -> (f64, f64) {
        (self.first, self.last)
    }

    /// The `Ideb`/`Ifin` window `Value`/`D1`/`D2` build for an end parameter
    /// (`Law_BSpFunc.cxx:256-282`, `:295-322`, `:336-358`).
    fn end_window(&self, x: f64) -> (i32, i32) {
        let curv = self.curve_ref();
        let pos_tol = PCONFUSION / 2.0;
        let nb_knots = curv.nb_knots() as i32;
        let (mut ideb, mut ifin) = curv.locate_u(x, pos_tol, false);
        if x == self.first {
            if ideb < 1 {
                ideb = 1;
            }
            if ideb >= ifin {
                ifin = ideb + 1;
            }
        }
        if x == self.last {
            if ifin > nb_knots {
                ifin = nb_knots;
            }
            if ideb >= ifin {
                ideb = ifin - 1;
            }
        }
        (ideb, ifin)
    }

    fn continuity_range(&self, s: Shape, what: &str) -> Result<i32, String> {
        Ok(match s {
            Shape::G1 | Shape::G2 => {
                return Err(format!("Standard_DomainError: {what}"))
            }
            Shape::C1 => 1,
            Shape::C2 => 2,
            Shape::C3 => 3,
            Shape::CN => self.curve_ref().degree(),
            _ => return Ok(-1),
        })
    }

    /// `Law_BSpFunc::NbIntervals` (`Law_BSpFunc.cxx:57-146`).
    fn nb_intervals_impl(&self, s: Shape) -> Result<i32, String> {
        let curv = self.curve_ref();
        let mut my_nb = 1;
        if s > curv.continuity() {
            let cont = self.continuity_range(s, "Law_BSpFunc::NbIntervals")?;
            if cont < 0 {
                return Ok(1);
            }
            let convector = LawBSplineKnotSplitting::new(curv, cont)?;
            let nb_int = convector.nb_splits() as i32 - 1;
            let inter = convector.splitting();
            let nb = curv.nb_knots() as i32;
            let (mut index1, new_first) =
                locate_parameter_knots_range(curv.knots(), self.first, curv.is_periodic(), 1, nb);
            let (mut index2, new_last) =
                locate_parameter_knots_range(curv.knots(), self.last, curv.is_periodic(), 1, nb);
            if let Some(k) = curv.knots().get(index1 as usize) {
                if (new_first - k).abs() < PCONFUSION {
                    index1 += 1;
                }
            }
            if let Some(k) = curv.knots().get((index2 - 1) as usize) {
                if new_last - k > PCONFUSION {
                    index2 += 1;
                }
            }
            my_nb = 1;
            for i in 0..nb_int.max(0) as usize {
                if let Some(&value) = inter.get(i) {
                    if value > index1 && value < index2 {
                        my_nb += 1;
                    }
                }
            }
        }
        Ok(my_nb)
    }

    /// `Law_BSpFunc::Intervals` (`Law_BSpFunc.cxx:152-249`).
    fn intervals_impl(&self, s: Shape) -> Result<Vec<f64>, String> {
        let curv = self.curve_ref();
        let mut my_nb = 1;
        let mut t: Vec<f64> = Vec::new();
        if s > curv.continuity() {
            let cont = self.continuity_range(s, "Law_BSpFunc_Curve::Intervals")?;
            if cont < 0 {
                // `GeomAbs_C0` (`cxx:165-167`): no split, the bounds below win.
                return Ok(vec![self.first, self.last]);
            }
            let convector = LawBSplineKnotSplitting::new(curv, cont)?;
            let nb_int = convector.nb_splits() as i32 - 1;
            let inter = convector.splitting();
            let nb = curv.nb_knots() as i32;
            let (mut index1, new_first) =
                locate_parameter_knots_range(curv.knots(), self.first, curv.is_periodic(), 1, nb);
            let (mut index2, new_last) =
                locate_parameter_knots_range(curv.knots(), self.last, curv.is_periodic(), 1, nb);
            if let Some(k) = curv.knots().get(index1 as usize) {
                if (new_first - k).abs() < PCONFUSION {
                    index1 += 1;
                }
            }
            if let Some(k) = curv.knots().get((index2 - 1) as usize) {
                if new_last - k > PCONFUSION {
                    index2 += 1;
                }
            }
            let mut kept: Vec<i32> = vec![index1];
            my_nb = 1;
            for i in 0..nb_int.max(0) as usize {
                if let Some(&value) = inter.get(i) {
                    if value > index1 && value < index2 {
                        my_nb += 1;
                        kept.push(value);
                    }
                }
            }
            kept.push(index2);
            t = kept
                .iter()
                .map(|&i| *curv.knots().get((i - 1) as usize).unwrap_or(&0.0))
                .collect();
        }
        t.resize(my_nb + 1, 0.0);
        t[0] = self.first;
        t[my_nb as usize] = self.last;
        Ok(t)
    }
}

impl LawFunction for LawBSpFunc {
    /// `Law_BSpFunc::Continuity` (`cxx:54-57`).
    fn continuity(&self) -> Shape {
        self.curve_ref().continuity()
    }

    /// `Law_BSpFunc::NbIntervals` (`cxx:57-146`).
    fn nb_intervals(&self, s: Shape) -> Result<i32, String> {
        self.nb_intervals_impl(s)
    }

    /// `Law_BSpFunc::Intervals` (`cxx:152-249`).
    fn intervals(&self, s: Shape) -> Result<Vec<f64>, String> {
        self.intervals_impl(s)
    }

    /// `Law_BSpFunc::Value` (`cxx:253-286`).
    fn value(&self, x: f64) -> f64 {
        let curv = self.curve_ref();
        if x == self.first || x == self.last {
            let (ideb, ifin) = self.end_window(x);
            *curv
                .local_value(x, ideb, ifin)
                .as_ref()
                .unwrap_or(&f64::NAN)
        } else {
            curv.value(x)
        }
    }

    /// `Law_BSpFunc::D1` (`cxx:292-327`).
    fn d1(&self, x: f64) -> (f64, f64) {
        let curv = self.curve_ref();
        if x == self.first || x == self.last {
            let (ideb, ifin) = self.end_window(x);
            curv.local_d1(x, ideb, ifin).unwrap_or((f64::NAN, f64::NAN))
        } else {
            curv.d1(x)
        }
    }

    /// `Law_BSpFunc::D2` (`cxx:333-366`).
    fn d2(&self, x: f64) -> (f64, f64, f64) {
        let curv = self.curve_ref();
        if x == self.first || x == self.last {
            let (ideb, ifin) = self.end_window(x);
            curv.local_d2(x, ideb, ifin)
                .unwrap_or((f64::NAN, f64::NAN, f64::NAN))
        } else {
            curv.d2(x)
        }
    }

    /// `Law_BSpFunc::Trim` (`cxx:371-378`): the same curve over the new
    /// bounds.
    fn trim(&self, p_first: f64, p_last: f64, _tol: f64) -> Result<Arc<dyn LawFunction>, String> {
        Ok(Arc::new(LawBSpFunc {
            curv: self.curv.clone(),
            first: p_first,
            last: p_last,
        }))
    }

    /// `Law_BSpFunc::Bounds` (`cxx:383-387`).
    fn bounds(&self) -> (f64, f64) {
        (self.first, self.last)
    }
}