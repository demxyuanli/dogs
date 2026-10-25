//! `ApproxInt_MultiLine` — the multi-line the approximation drivers fit,
//! over an intersection line plus the surface pair that can insert points.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/ApproxInt/ApproxInt_MultiLine.gxx:30-778`
//! (OCCT 8.0.0), instantiated per line type by the `#define TheLine …` blocks
//! of `BRepApprox_TheMultiLineOfApprox_0.cxx:16-29` and its `GeomInt` /
//! `IntPatch` siblings.
//!
//! The OCCT template is parameterised by two names:
//! - `Handle_TheLine` / `TheLine`: an `IntSurf_LineOn2S`-shaped container with
//!   `NbPnts()` and `Point(i)`, built either from the curves of a
//!   `BRepApprox_ApproxLine` or from a freshly collected point list. That is
//!   [`LineSource`] here.
//! - `TheSvSurfaces`: the `ApproxInt_SvSurfaces` implementation used to insert
//!   points (`Compute` / `SeekPoint` / `Tangency*`) and to carry the
//!   `myUseSolver` flag. That is [`SvSurfaces`](super::SvSurfaces), held as
//!   `Rc<RefCell<..>>` because OCCT reaches it through an untyped `void*`
//!   (`gxx:57-58`, `:330`) and mutates the flag from `const` methods.
//!
//! `Value` / `Tangency` write into `NCollection_Array1` arguments whose length
//! the caller picks (`gxx:206`, `:273`). The port answers the same values as
//! `value_p3d`, `value_p2d`, `tangency`, `tangency_2d`; the "length >= 2"
//! guards of `gxx:206-208` and `:271-276` become the length of the returned
//! `Vec`, which is `nbp2d`.

use std::cell::RefCell;
use std::rc::Rc;

use occt_core::approx::ApproxStatus;
use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::precision::SQUARE_CONFUSION;

use crate::int_tools_wline::PntOn2S;

use super::SvSurfaces;

/// The all-zero `IntSurf_PntOn2S` OCCT default-constructs at
/// `ApproxInt_MultiLine.gxx:407` and `:667`.
fn fresh_point() -> PntOn2S {
    PntOn2S {
        p: GpPnt::zero(),
        u1: 0.0,
        v1: 0.0,
        u2: 0.0,
        v2: 0.0,
    }
}

/// The `TheLine` half of the OCCT template: the point container the multiline
/// reads (`myLine->NbPnts()`, `myLine->Point(i)`) and the two ways OCCT builds
/// a new one (`new TheLine(vide1, false)` at `gxx:310`, `new
/// TheLine(ResultPntOn2SLine, false)` at `gxx:510`).
pub trait LineSource: Clone {
    /// `NbPnts()`.
    fn nb_pnts(&self) -> i32;

    /// `Point(i)` (`IntSurf_LineOn2S::Value(i)` for the point container, the
    /// pole evaluation for the curve-backed line): 1-based.
    fn point(&self, index: i32) -> PntOn2S;

    /// `new TheLine(IntSurf_LineOn2S, false)`: a line over `points`.
    fn from_points(points: Vec<PntOn2S>) -> Self;

    /// `new TheLine(new IntSurf_LineOn2S(), false)`: the empty line.
    fn empty() -> Self;
}

/// `ApproxInt_MultiLine` (`gxx:31-778`).
pub struct MultiLine<L: LineSource> {
    /// `PtrOnmySvSurfaces` (`gxx:37`): `NULL` in OCCT is `None` here.
    ptr_on_my_sv_surfaces: Option<Rc<RefCell<dyn SvSurfaces>>>,
    /// `myLine`.
    my_line: L,
    /// `indicemin`.
    indicemin: i32,
    /// `indicemax`.
    indicemax: i32,
    /// `nbp3d`.
    nbp3d: i32,
    /// `nbp2d`.
    nbp2d: i32,
    /// `myApproxU1V1`.
    my_approx_u1v1: bool,
    /// `myApproxU2V2`.
    my_approx_u2v2: bool,
    /// `p2donfirst`.
    p2donfirst: bool,
    /// `Xo`, `Yo`, `Zo`.
    xo: f64,
    yo: f64,
    zo: f64,
    /// `U1o`, `V1o`.
    u1o: f64,
    v1o: f64,
    /// `U2o`, `V2o`.
    u2o: f64,
    v2o: f64,
}

impl<L: LineSource> Clone for MultiLine<L> {
    fn clone(&self) -> Self {
        Self {
            ptr_on_my_sv_surfaces: self.ptr_on_my_sv_surfaces.clone(),
            my_line: self.my_line.clone(),
            indicemin: self.indicemin,
            indicemax: self.indicemax,
            nbp3d: self.nbp3d,
            nbp2d: self.nbp2d,
            my_approx_u1v1: self.my_approx_u1v1,
            my_approx_u2v2: self.my_approx_u2v2,
            p2donfirst: self.p2donfirst,
            xo: self.xo,
            yo: self.yo,
            zo: self.zo,
            u1o: self.u1o,
            v1o: self.v1o,
            u2o: self.u2o,
            v2o: self.v2o,
        }
    }
}

impl<L: LineSource> MultiLine<L> {
    /// `ApproxInt_MultiLine()` (`gxx:33-51`): every field at its initial value,
    /// `myLine` null (the port keeps an empty line).
    pub fn new() -> Self {
        Self {
            ptr_on_my_sv_surfaces: None,
            my_line: L::empty(),
            indicemin: 0,
            indicemax: 0,
            nbp3d: 0,
            nbp2d: 0,
            my_approx_u1v1: false,
            my_approx_u2v2: false,
            p2donfirst: true,
            xo: 0.0,
            yo: 0.0,
            zo: 0.0,
            u1o: 0.0,
            v1o: 0.0,
            u2o: 0.0,
            v2o: 0.0,
        }
    }

    /// `ApproxInt_MultiLine(line, svsurf, NbP3d, NbP2d, ApproxU1V1, ApproxU2V2,
    /// xo, yo, zo, u1o, v1o, u2o, v2o, P2DOnFirst, IndMin, IndMax)`
    /// (`gxx:55-94`).
    #[allow(clippy::too_many_arguments)]
    pub fn with_surfaces(
        line: L,
        svsurf: Rc<RefCell<dyn SvSurfaces>>,
        nb_p3d: i32,
        nb_p2d: i32,
        approx_u1v1: bool,
        approx_u2v2: bool,
        offsets: [f64; 7],
        p2d_on_first: bool,
        ind_min: i32,
        ind_max: i32,
    ) -> Self {
        let [xo, yo, zo, u1o, v1o, u2o, v2o] = offsets;
        Self {
            ptr_on_my_sv_surfaces: Some(svsurf),
            my_line: line,
            indicemin: ind_min.min(ind_max),
            indicemax: ind_min.max(ind_max),
            nbp3d: nb_p3d,
            nbp2d: nb_p2d,
            my_approx_u1v1: approx_u1v1,
            my_approx_u2v2: approx_u2v2,
            p2donfirst: p2d_on_first,
            xo,
            yo,
            zo,
            u1o,
            v1o,
            u2o,
            v2o,
        }
    }

    /// `ApproxInt_MultiLine(line, NbP3d, NbP2d, …)` (`gxx:98-137`): the arm
    /// without a surface pair, i.e. `PtrOnmySvSurfaces = nullptr`.
    #[allow(clippy::too_many_arguments)]
    pub fn without_surfaces(
        line: L,
        nb_p3d: i32,
        nb_p2d: i32,
        approx_u1v1: bool,
        approx_u2v2: bool,
        offsets: [f64; 7],
        p2d_on_first: bool,
        ind_min: i32,
        ind_max: i32,
    ) -> Self {
        let [xo, yo, zo, u1o, v1o, u2o, v2o] = offsets;
        Self {
            ptr_on_my_sv_surfaces: None,
            my_line: line,
            indicemin: ind_min.min(ind_max),
            indicemax: ind_min.max(ind_max),
            nbp3d: nb_p3d,
            nbp2d: nb_p2d,
            my_approx_u1v1: approx_u1v1,
            my_approx_u2v2: approx_u2v2,
            p2donfirst: p2d_on_first,
            xo,
            yo,
            zo,
            u1o,
            v1o,
            u2o,
            v2o,
        }
    }

    /// `FirstPoint()` (`gxx:142-145`).
    pub fn first_point(&self) -> i32 {
        self.indicemin
    }

    /// `LastPoint()` (`gxx:149-152`).
    pub fn last_point(&self) -> i32 {
        self.indicemax
    }

    /// `WhatStatus()` (`gxx:156-162`).
    pub fn what_status(&self) -> ApproxStatus {
        if self.ptr_on_my_sv_surfaces.is_some() {
            ApproxStatus::PointsAdded
        } else {
            ApproxStatus::NoPointsAdded
        }
    }

    /// `NbP3d()` (`gxx:166-169`).
    pub fn nb_p3d(&self) -> i32 {
        self.nbp3d
    }

    /// `NbP2d()` (`gxx:173-176`).
    pub fn nb_p2d(&self) -> i32 {
        self.nbp2d
    }

    /// `Value(Index, TabPnt)` (`gxx:180-183`): the 3D point offset by `(Xo,
    /// Yo, Zo)`.
    pub fn value_p3d(&self, index: i32) -> GpPnt {
        let a_p = self.my_line.point(index).value();
        GpPnt::new(a_p.x() + self.xo, a_p.y() + self.yo, a_p.z() + self.zo)
    }

    /// `Value(Index, TabPnt2d)` (`gxx:187-210`): as many entries as `NbP2d`,
    /// with the second one only when `NbP2d >= 2`.
    pub fn value_p2d(&self, index: i32) -> Vec<GpPnt2d> {
        let (u1, v1, u2, v2) = self.my_line.point(index).parameters();
        if self.nbp2d == 1 {
            let (u, v) = if self.p2donfirst {
                (u1 + self.u1o, v1 + self.v1o)
            } else {
                (u2 + self.u2o, v2 + self.v2o)
            };
            return vec![GpPnt2d::new(u, v)];
        }
        let mut out = vec![GpPnt2d::new(u1 + self.u1o, v1 + self.v1o)];
        if self.nbp2d >= 2 {
            out.push(GpPnt2d::new(u2 + self.u2o, v2 + self.v2o));
        }
        out
    }

    /// `Value(Index, TabPnt, TabPnt2d)` (`gxx:214-221`).
    pub fn value(&self, index: i32) -> (GpPnt, Vec<GpPnt2d>) {
        (self.value_p3d(index), self.value_p2d(index))
    }

    /// `Tangency(Index, TabVec)` (`gxx:225-240`): `false` without a surface
    /// pair, and the tangent is zeroed when the evaluator fails.
    pub fn tangency(&self, index: i32) -> (bool, GpVec) {
        let Some(sv) = &self.ptr_on_my_sv_surfaces else {
            return (false, GpVec::new(0.0, 0.0, 0.0));
        };
        let (u1, v1, u2, v2) = self.my_line.point(index).parameters();
        let mut tg = GpVec::new(0.0, 0.0, 0.0);
        let ret = sv
            .borrow_mut()
            .tangency(u1, v1, u2, v2, &mut tg);
        if !ret {
            tg = GpVec::new(0.0, 0.0, 0.0);
        }
        (ret, tg)
    }

    /// `Tangency(Index, TabVec2d)` (`gxx:244-289`): `TangencyOnSurf1` /
    /// `TangencyOnSurf2` depending on `NbP2d` and `P2DOnFirst`.
    pub fn tangency_2d(&self, index: i32) -> (bool, Vec<GpVec2d>) {
        let Some(sv) = &self.ptr_on_my_sv_surfaces else {
            return (false, self.zero_tangents_2d());
        };
        let (u1, v1, u2, v2) = self.my_line.point(index).parameters();
        let mut t1 = GpVec2d::new(0.0, 0.0);
        let mut t2 = GpVec2d::new(0.0, 0.0);
        let mut ret;
        if self.nbp2d == 1 {
            let mut b = sv.borrow_mut();
            ret = if self.p2donfirst {
                b.tangency_on_surf1(u1, v1, u2, v2, &mut t1)
            } else {
                b.tangency_on_surf2(u1, v1, u2, v2, &mut t1)
            };
        } else {
            let mut b = sv.borrow_mut();
            ret = b.tangency_on_surf1(u1, v1, u2, v2, &mut t1);
            if ret && self.nbp2d >= 2 {
                ret = ret && b.tangency_on_surf2(u1, v1, u2, v2, &mut t2);
            }
        }
        if !ret {
            t1 = GpVec2d::new(0.0, 0.0);
            if self.nbp2d >= 2 {
                t2 = GpVec2d::new(0.0, 0.0);
            }
        }
        let mut out = vec![t1];
        if self.nbp2d >= 2 {
            out.push(t2);
        }
        (ret, out)
    }

    /// `Tangency(Index, TabVec, TabVec2d)` (`gxx:293-298`): both must succeed.
    pub fn tangency_both(&self, index: i32) -> (bool, GpVec, Vec<GpVec2d>) {
        let (ret1, tg) = self.tangency(index);
        let (ret2, tg2d) = self.tangency_2d(index);
        (ret1 && ret2, tg, tg2d)
    }

    /// The zeroed array `Tangency(Index, TabVec2d)` leaves behind on failure.
    fn zero_tangents_2d(&self) -> Vec<GpVec2d> {
        let mut out = vec![GpVec2d::new(0.0, 0.0)];
        if self.nbp2d >= 2 {
            out.push(GpVec2d::new(0.0, 0.0));
        }
        out
    }

    /// `MakeMLBetween(Low, High, aNbPntsToInsert)` (`gxx:302-626`): insert
    /// points at constant 3D arc length between two indices, re-evaluating the
    /// parameters on the surface pair.
    ///
    /// `None` is the empty line OCCT returns; every non-empty result carries
    /// the surface pair, so the caller can go on inserting.
    pub fn make_ml_between(&self, low: i32, high: i32, a_nb_pnts_to_insert: i32) -> Self {
        let Some(sv) = self.ptr_on_my_sv_surfaces.clone() else {
            return self.empty_with(1, 1);
        };

        // `gxx:330-334`
        let a_save_use_solver = sv.borrow().get_use_solver();
        if !a_save_use_solver {
            sv.borrow_mut().set_use_solver(true);
        }

        // `gxx:335-343`
        let mut nb_pnts_to_insert = a_nb_pnts_to_insert;
        if nb_pnts_to_insert < (high - low) {
            nb_pnts_to_insert = high - low;
        }
        let mut nb_pnts = nb_pnts_to_insert + high - low + 1;
        let mut nb_pnts_min = high - low;
        nb_pnts_min += nb_pnts_min;
        if nb_pnts < nb_pnts_min {
            nb_pnts = nb_pnts_min;
        }

        // `gxx:357-365`: `NCollection_Array1(Low, High)` becomes `Vec` indexed
        // by `i - Low`.
        let span = (high - low + 1) as usize;
        let mut u1s = vec![0.0; span];
        let mut v1s = vec![0.0; span];
        let mut u2s = vec![0.0; span];
        let mut v2s = vec![0.0; span];
        let mut ac = vec![0.0; span];
        let idx = |i: i32| (i - low) as usize;

        // `gxx:370-376`
        let (u1, v1, u2, v2) = self.my_line.point(low).parameters();
        u1s[idx(low)] = u1;
        v1s[idx(low)] = v1;
        u2s[idx(low)] = u2;
        v2s[idx(low)] = v2;
        ac[idx(low)] = 0.0;

        // `gxx:392-400`: arc length in 3D.
        for i in (low + 1)..=high {
            let (u1, v1, u2, v2) = self.my_line.point(i).parameters();
            u1s[idx(i)] = u1;
            v1s[idx(i)] = v1;
            u2s[idx(i)] = u2;
            v2s[idx(i)] = v2;
            ac[idx(i)] = ac[idx(i - 1)]
                + self
                    .my_line
                    .point(i - 1)
                    .value()
                    .distance(&self.my_line.point(i).value());
        }

        // `gxx:405-510`
        let mut result: Vec<PntOn2S> = Vec::new();
        let ds = ac[idx(high)] / (nb_pnts - 1) as f64;
        let mut indice = low;
        let mut has_been_inserted = false;
        let dsmin = ds * 0.3;
        let smax = ac[idx(high)];

        // `for (i = 2, s = ds; (s < smax && Indice <= High - 1); i++, s += ds)`
        // (`gxx:418`): `i` is never read inside the body, only `s` drives it.
        let mut s = ds;
        while s < smax && indice <= high - 1 {
            while ac[idx(indice + 1)] <= s {
                if !has_been_inserted {
                    result.push(self.my_line.point(indice));
                }
                has_been_inserted = false;
                indice += 1;
                if indice == high {
                    break;
                }
            }
            if indice == high {
                break;
            }
            if !has_been_inserted && ac[idx(indice)] <= s {
                result.push(self.my_line.point(indice));
                has_been_inserted = true;
            }

            let a = s - ac[idx(indice)];
            let b = ac[idx(indice + 1)] - s;
            let nab = 1.0 / (a + b);

            let mut inserted = false;
            if a > dsmin && b > dsmin {
                let mut u1 = (u1s[idx(indice)] * b + u1s[idx(indice + 1)] * a) * nab;
                let mut v1 = (v1s[idx(indice)] * b + v1s[idx(indice + 1)] * a) * nab;
                let mut u2 = (u2s[idx(indice)] * b + u2s[idx(indice + 1)] * a) * nab;
                let mut v2 = (v2s[idx(indice)] * b + v2s[idx(indice + 1)] * a) * nab;
                let mut p = GpPnt::zero();
                let mut t = GpVec::new(0.0, 0.0, 0.0);
                let mut ts1 = GpVec2d::new(0.0, 0.0);
                let mut ts2 = GpVec2d::new(0.0, 0.0);
                if sv
                    .borrow_mut()
                    .compute(&mut u1, &mut v1, &mut u2, &mut v2, &mut p, &mut t, &mut ts1, &mut ts2)
                {
                    let mut p_on2s = fresh_point();
                    p_on2s.set_value(p, u1, v1, u2, v2);
                    result.push(p_on2s);
                    inserted = true;
                }
            }
            if !inserted {
                if b < 0.0 {
                    while ac[idx(indice + 1)] <= s {
                        if !has_been_inserted {
                            result.push(self.my_line.point(indice));
                        }
                        has_been_inserted = false;
                        indice += 1;
                        if indice == high {
                            break;
                        }
                    }
                    if indice == high {
                        break;
                    }
                    if !has_been_inserted && ac[idx(indice)] <= s {
                        result.push(self.my_line.point(indice));
                        has_been_inserted = true;
                    }
                } else {
                    s += dsmin - ds;
                }
            }

            s += ds;
        }

        result.push(self.my_line.point(high));

        // `gxx:511-568`: the a-posteriori check on the 2D turns.
        let temp = L::from_points(result.clone());
        let nb_temp = temp.nb_pnts();
        let (mut u1, mut v1, mut u2, mut v2) = temp.point(1).parameters();
        let mut p1a = GpPnt2d::new(u1, v1);
        let mut p2a = GpPnt2d::new(u2, v2);
        let (u1b, v1b, u2b, v2b) = temp.point(2).parameters();
        let mut p1b = GpPnt2d::new(u1b, v1b);
        let mut p2b = GpPnt2d::new(u2b, v2b);
        let mut code_erreur = 0;
        let mut i = 3;
        loop {
            if !(code_erreur == 0 && i <= nb_temp) {
                break;
            }
            let (u1c, v1c, u2c, v2c) = temp.point(i).parameters();
            let p1c = GpPnt2d::new(u1c, v1c);
            let mut du = p1b.x() - p1a.x();
            let mut dv = p1b.y() - p1a.y();
            let duv2 = 0.25 * (du * du + dv * dv);
            u1 = p1b.x() + du;
            v1 = p1b.y() + dv;
            du = p1c.x() - u1;
            dv = p1c.y() - v1;
            if du * du + dv * dv > duv2 {
                code_erreur = 1;
                break;
            }

            let p2c = GpPnt2d::new(u2c, v2c);
            du = p2b.x() - p2a.x();
            dv = p2b.y() - p2a.y();
            let duv2 = 0.25 * (du * du + dv * dv);
            u2 = p2b.x() + du;
            v2 = p2b.y() + dv;
            du = p2c.x() - u2;
            dv = p2c.y() - v2;
            if du * du + dv * dv > duv2 {
                code_erreur = 2;
                break;
            }

            p1a = p1b;
            p2a = p2b;
            p1b = p1c;
            p2b = p2c;
            i += 1;
        }

        // `gxx:585-626`
        if nb_temp >= nb_pnts_to_insert + high - low + 1 && code_erreur == 0 {
            sv.borrow_mut().set_use_solver(a_save_use_solver);
            let sv2 = if high - low > 10 {
                Some(sv)
            } else {
                None
            };
            return self.span_with(temp, sv2, 1, nb_temp);
        }
        sv.borrow_mut().set_use_solver(a_save_use_solver);
        self.empty_with(1, 1)
    }

    /// `MakeMLOneMorePoint(theLow, theHigh, theIndbad, theNewMultiLine)`
    /// (`gxx:634-751`): insert one point between `theIndbad - 1` and
    /// `theIndbad`. `false` when the surface pair refuses, when the mid-point
    /// does not move in parameter space or lies too close to a neighbour.
    pub fn make_ml_one_more_point(&self, the_low: i32, the_high: i32, the_indbad: i32) -> (bool, Self) {
        let Some(sv) = self.ptr_on_my_sv_surfaces.clone() else {
            return (false, self.empty_with(0, 0));
        };

        // `gxx:640-645`
        let a_save_use_solver = sv.borrow().get_use_solver();
        if !a_save_use_solver {
            sv.borrow_mut().set_use_solver(true);
        }

        // `gxx:647-680`
        const SQ_TOL3D: f64 = SQUARE_CONFUSION;
        let tolerance = 1.0e-8;

        let mut result: Vec<PntOn2S> = Vec::new();
        for indice in the_low..=the_high {
            result.push(self.my_line.point(indice));
        }

        let prev_pnt = self.my_line.point(the_indbad - 1).value();
        let cur_pnt = self.my_line.point(the_indbad).value();
        let (uprev1, vprev1, uprev2, vprev2) = self.my_line.point(the_indbad - 1).parameters();
        let (ucur1, vcur1, ucur2, vcur2) = self.my_line.point(the_indbad).parameters();
        let umid1 = (uprev1 + ucur1) / 2.0;
        let vmid1 = (vprev1 + vcur1) / 2.0;
        let umid2 = (uprev2 + ucur2) / 2.0;
        let vmid2 = (vprev2 + vcur2) / 2.0;

        let mut is_new_point_invalid = self.my_approx_u1v1
            && (ucur1 - umid1).abs() <= tolerance
            && (vcur1 - vmid1).abs() <= tolerance;
        let mut other_line_made = false;
        if !is_new_point_invalid {
            is_new_point_invalid = self.my_approx_u2v2
                && (ucur2 - umid2).abs() <= tolerance
                && (vcur2 - vmid2).abs() <= tolerance;
            let mut mid_point = fresh_point();
            if !is_new_point_invalid
                && sv
                    .borrow_mut()
                    .seek_point(umid1, vmid1, umid2, vmid2, &mut mid_point)
            {
                let new_pnt = mid_point.value();
                let sq_dist_new_prev = new_pnt.square_distance(&prev_pnt);
                let sq_dist_new_cur = new_pnt.square_distance(&cur_pnt);
                is_new_point_invalid = sq_dist_new_prev <= SQ_TOL3D || sq_dist_new_cur <= SQ_TOL3D;
                if !is_new_point_invalid {
                    let (unew1, vnew1, unew2, vnew2) = mid_point.parameters();
                    if self.my_approx_u1v1 {
                        let sq_dist_cur_mid1 = (ucur1 - umid1) * (ucur1 - umid1)
                            + (vcur1 - vmid1) * (vcur1 - vmid1);
                        let sq_dist_mid_new1 = (umid1 - unew1) * (umid1 - unew1)
                            + (vmid1 - vnew1) * (vmid1 - vnew1);
                        is_new_point_invalid = sq_dist_mid_new1 > sq_dist_cur_mid1;
                    }
                    if !is_new_point_invalid {
                        if self.my_approx_u2v2 {
                            let sq_dist_cur_mid2 = (ucur2 - umid2) * (ucur2 - umid2)
                                + (vcur2 - vmid2) * (vcur2 - vmid2);
                            let sq_dist_mid_new2 = (umid2 - unew2) * (umid2 - unew2)
                                + (vmid2 - vnew2) * (vmid2 - vnew2);
                            is_new_point_invalid = sq_dist_mid_new2 > sq_dist_cur_mid2;
                        }
                        if !is_new_point_invalid {
                            // `ResultPntOn2SLine->InsertBefore(theIndbad - theLow + 1, MidPoint)`
                            result.insert((the_indbad - the_low) as usize, mid_point);
                            other_line_made = true;
                        }
                    }
                }
            }
        }

        if !other_line_made {
            sv.borrow_mut().set_use_solver(a_save_use_solver);
            return (false, self.empty_with(0, 0));
        }

        let temp = L::from_points(result);
        sv.borrow_mut().set_use_solver(a_save_use_solver);
        let nb = self.nb_pnts_of(&temp);
        (true, self.span_with(temp, Some(sv), 1, nb))
    }

    /// `myLine->NbPnts()` of an already-built line.
    fn nb_pnts_of(&self, line: &L) -> i32 {
        line.nb_pnts()
    }

    /// The empty `TheLine(new IntSurf_LineOn2S(), false)` result of
    /// `gxx:309-326` / `gxx:608-625`, keeping this multiline's settings.
    fn empty_with(&self, ind_min: i32, ind_max: i32) -> Self {
        Self {
            ptr_on_my_sv_surfaces: None,
            my_line: L::empty(),
            indicemin: ind_min,
            indicemax: ind_max,
            nbp3d: self.nbp3d,
            nbp2d: self.nbp2d,
            my_approx_u1v1: self.my_approx_u1v1,
            my_approx_u2v2: self.my_approx_u2v2,
            p2donfirst: self.p2donfirst,
            xo: self.xo,
            yo: self.yo,
            zo: self.zo,
            u1o: self.u1o,
            v1o: self.v1o,
            u2o: self.u2o,
            v2o: self.v2o,
        }
    }

    /// The `ApproxInt_MultiLine(temp, svsurf, …, 1, indMax)` result shared by
    /// `gxx:587-602` and `gxx:738-751`.
    fn span_with(
        &self,
        line: L,
        sv: Option<Rc<RefCell<dyn SvSurfaces>>>,
        ind_min: i32,
        ind_max: i32,
    ) -> Self {
        Self {
            ptr_on_my_sv_surfaces: sv,
            my_line: line,
            indicemin: ind_min,
            indicemax: ind_max,
            nbp3d: self.nbp3d,
            nbp2d: self.nbp2d,
            my_approx_u1v1: self.my_approx_u1v1,
            my_approx_u2v2: self.my_approx_u2v2,
            p2donfirst: self.p2donfirst,
            xo: self.xo,
            yo: self.yo,
            zo: self.zo,
            u1o: self.u1o,
            v1o: self.v1o,
            u2o: self.u2o,
            v2o: self.v2o,
        }
    }
}

impl<L: LineSource> Default for MultiLine<L> {
    fn default() -> Self {
        Self::new()
    }
}

