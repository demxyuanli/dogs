//! `Contap_TheSearch` - the restriction-arc search of `Contap_Contour`.
//!
//! `Contap_TheSearch` is `IntStart_SearchOnBoundaries` instantiated with
//! `Contap_ArcFunction` / `Contap_HContTool` (`Contap_TheSearch_0.cxx:32-54`),
//! so this module ports `IntStart_SearchOnBoundaries.gxx:128-1231`:
//! `FindVertex` (`:129-164`), `SolInfo` (`:191-228`), `BoundedArc`
//! (`:229-798`), `ComputeBoundsfromInfinite` (`:802-877`), `PointProcess`
//! (`:880-1001`), `IsRegularity` (`:1005-1017`), `TreatLC` (`:1021-1125`) and
//! `Perform` (`:1127-1231`).
//!
//! Two arms of `BoundedArc` cannot be reached by this instantiation and are
//! left UNPORTED:
//!
//! * the exact `IntCurveSurface_HInter` arm (`gxx:375-476`, `:475-488`,
//!   `:756-770`, `:610-624`, `:670-681`): it runs only when
//!   `Func.Quadric().TypeQuadric() != GeomAbs_OtherSurface`, and
//!   `Contap_ArcFunction::Quadric()` (`Contap_ArcFunction.cxx:204-207`) answers
//!   the `myQuad` member (`hxx:79`) that no `Contap` code ever assigns, so it
//!   is the default-constructed `IntSurf_Quadric` (`IntSurf_Quadric.cxx:23-30`,
//!   `typ = GeomAbs_OtherSurface`). `IsIntCSdone` is therefore always `false`
//!   and `TypeConS` is never read.
//! * the cylinder-tangency body of `TreatLC` (`gxx:1050-1121`): it returns
//!   early at `:1052-1055` because `aQuadric.TypeQuadric() != GeomAbs_Cylinder`
//!   for the same reason. The guards that run (`:1023-1047`) are ported, and
//!   the function always answers `1`, which is its own "no tangent line found"
//!   exit code.
//!
//! `Contap_HContTool::HasBeenSeen` always answers `false` (`cxx:320-323`), so
//! the `else` arm of `Perform` (`gxx:1190-1229`) is dead code in OCCT as well;
//! the port does not transcribe it.

use std::cmp::Ordering;
use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_core::math_function_all_roots::FunctionAllRoots;
use occt_core::math_function_sample::FunctionSample;
use occt_core::precision::{CONFUSION, INFINITE, PCONFUSION};
use occt_math::brent_minimum::BrentMinimum;

use occt_topo::brep_adaptor_curve2d::Curve2dAdaptor;
use occt_topo::brep_lib_continuity;
use occt_topo::brep_top_adaptor_topol_tool::TopolTool;

use super::arc_function::ContapArcFunction;
use super::path_point::ContapThePathPointOfTheSearch;
use super::sample;
use super::segment::ContapTheSegmentOfTheSearch;

/// `RealLast()` (`Standard_Real.hxx`).
const REAL_LAST: f64 = f64::MAX;

/// `Contap_TheSearch` (`hxx:41-94`).
#[derive(Clone)]
pub struct ContapTheSearch {
    /// `done` (`hxx:91`).
    done: bool,
    /// `all` (`hxx:92`).
    all: bool,
    /// `sseg` (`hxx:93`).
    sseg: Vec<ContapTheSegmentOfTheSearch>,
    /// `spnt` (`hxx:94`).
    spnt: Vec<ContapThePathPointOfTheSearch>,
}

impl Default for ContapTheSearch {
    fn default() -> Self {
        Self::new()
    }
}

impl ContapTheSearch {
    /// `Contap_TheSearch()` (`gxx:1127-1131`).
    pub fn new() -> Self {
        Self {
            done: false,
            all: false,
            sseg: Vec::new(),
            spnt: Vec::new(),
        }
    }

    /// `Perform(F, Domain, TolBoundary, TolTangency, RecheckOnRegularity = false)`
    /// (`gxx:1135-1231`).
    pub fn perform(
        &mut self,
        func: &mut ContapArcFunction,
        domain: &mut TopolTool,
        tol_boundary: f64,
        tol_tangency: f64,
        recheck_on_regularity: bool,
    ) -> Result<(), String> {
        self.done = false;
        self.spnt.clear();
        self.sseg.clear();

        domain.init();

        let mut all = domain.more();

        while domain.more() {
            let Some(arc) = domain.value().cloned().map(Arc::new) else {
                return Err("Contap_TheSearch::Perform: domain arc missing".to_string());
            };
            // `gxx:1167`: `HasBeenSeen(A)` is always false (`cxx:320-323`), so
            // the "already known arc" arm (`gxx:1190-1229`) is unreachable.
            func.set_arc(arc.clone());
            find_vertex(&arc, domain, func, &mut self.spnt, tol_boundary)?;
            let (mut pdeb, mut pfin) = (arc.first_parameter(), arc.last_parameter());
            // `gxx:1170-1174`: `Precision::IsNegativeInfinite` /
            // `IsPositiveInfinite`.
            if pdeb == -INFINITE || pfin == INFINITE {
                let (d, f, _nb_echant) = compute_bounds_from_infinite(func);
                pdeb = d;
                pfin = f;
            }
            let arcsol = bounded_arc(
                &arc,
                domain,
                func,
                pdeb,
                pfin,
                &mut self.spnt,
                &mut self.sseg,
                tol_boundary,
                tol_tangency,
                recheck_on_regularity,
            )?;
            all = all && arcsol;

            domain.next();
        }

        self.all = all;
        self.done = true;
        Ok(())
    }

    /// `IsDone()` (`lxx:19-22`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `AllArcSolution()` (`lxx:24-32`): `Err` where OCCT raises
    /// `StdFail_NotDone`.
    pub fn all_arc_solution(&self) -> Result<bool, &'static str> {
        if !self.done {
            return Err("Contap_TheSearch::AllArcSolution: StdFail_NotDone");
        }
        Ok(self.all)
    }

    /// `NbPoints()` (`lxx:34-41`).
    pub fn nb_points(&self) -> Result<usize, &'static str> {
        if !self.done {
            return Err("Contap_TheSearch::NbPoints: StdFail_NotDone");
        }
        Ok(self.spnt.len())
    }

    /// `Point(Index)` (`lxx:43-51`), 1-based.
    pub fn point(&self, index: usize) -> Option<&ContapThePathPointOfTheSearch> {
        if !self.done {
            return None;
        }
        index.checked_sub(1).and_then(|i| self.spnt.get(i))
    }

    /// `NbSegments()` (`lxx:53-60`).
    pub fn nb_segments(&self) -> Result<usize, &'static str> {
        if !self.done {
            return Err("Contap_TheSearch::NbSegments: StdFail_NotDone");
        }
        Ok(self.sseg.len())
    }

    /// `Segment(Index)` (`lxx:62-70`), 1-based.
    pub fn segment(&self, index: usize) -> Option<&ContapTheSegmentOfTheSearch> {
        if !self.done {
            return None;
        }
        index.checked_sub(1).and_then(|i| self.sseg.get(i))
    }
}

/// `FindVertex(A, Domain, Func, pnt, Toler)` (`gxx:129-164`).
fn find_vertex(
    arc: &Arc<Curve2dAdaptor>,
    domain: &mut TopolTool,
    func: &mut ContapArcFunction,
    points: &mut Vec<ContapThePathPointOfTheSearch>,
    toler: f64,
) -> Result<(), String> {
    domain.initialize_curve(arc);
    domain.init_vertex_iterator();
    while domain.more_vertex() {
        let vtx = domain
            .vertex()
            .ok_or_else(|| "Contap_TheSearch::FindVertex: missing vertex".to_string())?;
        let param = sample::parameter(&vtx, arc);

        // `gxx:157-158`: the value is evaluated and only its magnitude tested.
        let mut valf = 0.0;
        func.eval(param, &mut valf);
        if valf.abs() <= toler {
            let itemp = func.get_state_number();
            let p = *func
                .valpoint(itemp as usize)
                .ok_or_else(|| "Contap_TheSearch::FindVertex: missing state point".to_string())?;
            points.push(ContapThePathPointOfTheSearch::on_vertex(
                p,
                toler,
                vtx,
                arc.clone(),
                param,
            ));
        }
        domain.next_vertex();
    }
    Ok(())
}

/// `BoundedArc(A, Domain, Pdeb, Pfin, Func, pnt, seg, TolBoundary,
/// TolTangency, Arcsol, RecheckOnRegularity)` (`gxx:229-798`). Answers
/// `Arcsol`.
#[allow(clippy::too_many_arguments)]
fn bounded_arc(
    arc: &Arc<Curve2dAdaptor>,
    domain: &mut TopolTool,
    func: &mut ContapArcFunction,
    pdeb: f64,
    pfin: f64,
    points: &mut Vec<ContapThePathPointOfTheSearch>,
    segments: &mut Vec<ContapTheSegmentOfTheSearch>,
    tol_boundary: f64,
    tol_tangency: f64,
    recheck_on_regularity: bool,
) -> Result<bool, String> {
    // `gxx:261`: `double EpsX = 1.e-10;`.
    let mut eps_x = 1.0e-10;

    // `gxx:268-270`: `int NbEchant = Func.NbSamples(); if (NbEchant < 100)
    // NbEchant = 100;`.
    let nb_echant = func.nb_samples().max(100);

    // `gxx:272-276`.
    let mut n_tol_tangency = tol_tangency;
    if (pfin - pdeb) < (tol_tangency * 10.0) {
        n_tol_tangency = (pfin - pdeb) * 0.1;
    }
    // `gxx:278-281`.
    if eps_x > (n_tol_tangency + n_tol_tangency) {
        eps_x = n_tol_tangency * 0.1;
    }

    // `gxx:292-333`: the rejection of an arc whose function never crosses zero.
    let mut minr = REAL_LAST;
    let mut maxr = -minr;
    let mut max_dr = -minr;
    let dur = (pfin - pdeb) * 0.2;
    let mut ur = pdeb;
    let mut rejection = true;
    for _ in 1..=6 {
        let mut f = 0.0;
        let mut d = 0.0;
        if func.eval_f_and_d(ur, &mut f, &mut d) {
            let mut d_abs = d.abs();
            d_abs *= dur + dur;
            if d_abs > max_dr {
                max_dr = d_abs;
            }
            let lminr = f - d_abs;
            let lmaxr = f + d_abs;
            if lminr < minr {
                minr = lminr;
            }
            if lmaxr > maxr {
                maxr = lmaxr;
            }
            if minr < 0.0 && maxr > 0.0 {
                rejection = false;
                break;
            }
        }
        ur += dur;
    }
    if rejection {
        let dur = 0.001 + max_dr + (maxr - minr) * 0.1;
        minr -= dur;
        maxr += dur;
        if minr < 0.0 && maxr > 0.0 {
            rejection = false;
        }
    }

    let mut arcsol = false;
    if rejection {
        return Ok(arcsol);
    }

    // `gxx:352-371`: widen the point tolerance on a fine arc, otherwise clamp
    // it back to `TolBoundary`.
    let mut maxdist = tol_boundary + tol_tangency;
    let echant0 = FunctionSample::new(pdeb, pfin, nb_echant);
    let mut aelargir = true;
    let mut i = 1;
    while i <= nb_echant && aelargir {
        let u = echant0.parameter(i);
        let mut dist = 0.0;
        if func.eval(u, &mut dist) {
            if dist > maxdist || -dist > maxdist {
                aelargir = false;
            }
        }
        i += 1;
    }
    if !(aelargir && maxdist < 0.01) {
        maxdist = tol_boundary;
    }

    // `gxx:375-476` (exact `IntCurveSurface_HInter`) is UNPORTED: it needs
    // `Func.Quadric().TypeQuadric() != GeomAbs_OtherSurface`, which never holds
    // for `Contap_ArcFunction` (see the module note).

    // `gxx:478-489`.
    let echant = FunctionSample::new(pdeb, pfin, nb_echant);
    let p_sol = FunctionAllRoots::new(func, &echant, eps_x, maxdist, maxdist);
    if !p_sol.is_done() {
        return Err("Contap_TheSearch::BoundedArc: math_FunctionAllRoots not done".to_string());
    }
    let nbp = p_sol.nb_points();

    // `gxx:493-560`: the regularity recheck.
    if recheck_on_regularity && nbp > 0 && is_regularity(domain) {
        let the_tol = 5.0e-4;
        let sol_again = FunctionAllRoots::new(func, &echant, eps_x, the_tol, the_tol);
        if !sol_again.is_done() {
            return Err(
                "Contap_TheSearch::BoundedArc: regularity math_FunctionAllRoots not done"
                    .to_string(),
            );
        }
        let nbi_again = sol_again.nb_intervals();
        if nbi_again > 0 {
            let nb_samples = 10;
            let delta = (pfin - pdeb) / f64::from(nb_samples);
            let global_tol = the_tol * 10.0;
            let mut sol_on_boundary = true;
            for step in 0..=nb_samples {
                let a_param = pdeb + f64::from(step) * delta;
                let mut a_value = 0.0;
                func.eval(a_param, &mut a_value);
                if a_value.abs() > global_tol {
                    sol_on_boundary = false;
                    break;
                }
            }
            if sol_on_boundary {
                for seg_i in 1..=nbi_again {
                    let (mut pardeb, mut parfin) = sol_again.interval(seg_i);
                    if (pardeb - pdeb).abs() <= PCONFUSION {
                        pardeb = pdeb;
                    }
                    if (parfin - pfin).abs() <= PCONFUSION {
                        parfin = pfin;
                    }
                    let (ideb, ifin) = sol_again.interval_state(seg_i);
                    let ptdeb = *func.valpoint(ideb as usize).ok_or_else(|| {
                        "Contap_TheSearch::BoundedArc: missing interval start".to_string()
                    })?;
                    let ptfin = *func.valpoint(ifin as usize).ok_or_else(|| {
                        "Contap_TheSearch::BoundedArc: missing interval end".to_string()
                    })?;

                    let mut newseg = ContapTheSegmentOfTheSearch::new();
                    newseg.set_value(arc.clone());
                    let ranged = point_process(domain, points, ptdeb, pardeb, arc, the_tol)?;
                    newseg.set_limit_point(&points[ranged - 1], true);
                    let rangef = point_process(domain, points, ptfin, parfin, arc, the_tol)?;
                    newseg.set_limit_point(&points[rangef - 1], false);
                    segments.push(newseg);
                }
                return Ok(true);
            }
        }
    }

    // `gxx:566-739`: the near-tangency treatment and the point collection.
    let mut nbi;
    if nbp > 0 {
        // `gxx:570-580`: `SolInfo` copies of the sorted point parameters.
        let mut a_si: Vec<(f64, i32)> = (1..=nbp).map(|k| (p_sol.point(k), k)).collect();
        a_si.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal));

        // `gxx:590`: `int ip = TreatLC(...)`; `1` means "no tangent line
        // found", i.e. the usual path below still runs.
        let ip = treat_lc(domain, tol_boundary, points)?;
        if ip != 0 {
            for k in 1..nbp as usize {
                let mut parap1 = a_si[k].0;
                let mut para = a_si[k - 1].0;

                let mut param = (para + parap1) * 0.5;
                let mut ym = 0.0;
                let mut yf = 0.0;
                let mut yl = 0.0;
                if func.eval(param, &mut ym) && ym.abs() < maxdist {
                    let sm = ym.signum();
                    let mut a_tang = func.eval(para, &mut yf) && func.eval(parap1, &mut yl);
                    if a_tang {
                        // `gxx:613-616`.
                        a_tang = yf.abs() < maxdist && yl.abs() < maxdist;
                    }
                    // `gxx:617-625` (exact-intersection line tangency) is
                    // UNPORTED: it needs `IsIntCSdone` / `TypeConS`.
                    let _ = sm;
                    if a_tang {
                        // `gxx:628-668`.
                        let mut a_tol = tol_boundary * 1000.0;
                        if a_tol > 0.001 {
                            a_tol = 0.001;
                        }
                        // `gxx:636-637`: clamp the pair to +-1.e9.
                        parap1 = if parap1.abs() < 1.0e9 {
                            parap1
                        } else if parap1 >= 0.0 {
                            1.0e9
                        } else {
                            -1.0e9
                        };
                        para = if para.abs() < 1.0e9 {
                            para
                        } else if para >= 0.0 {
                            1.0e9
                        } else {
                            -1.0e9
                        };

                        let a_nb_nodes = ((parap1 - para) / a_tol).ceil() as i32;
                        let mut a_val = REAL_LAST;
                        let mut a_val_max = 0.0;
                        let a_delta = (parap1 - para) / (f64::from(a_nb_nodes) + 1.0);
                        for ii in 0..=(a_nb_nodes + 1) {
                            let a_cur_par = if ii < a_nb_nodes + 1 {
                                para + f64::from(ii) * a_delta
                            } else {
                                parap1
                            };
                            let mut a_cur_val = 0.0;
                            if func.eval(a_cur_par, &mut a_cur_val) {
                                let an_abs_val = a_cur_val.abs();
                                if an_abs_val < a_val {
                                    a_val = an_abs_val;
                                    param = a_cur_par;
                                }
                                if an_abs_val > a_val_max {
                                    a_val_max = an_abs_val;
                                }
                            }
                        }
                        // `gxx:669-679` (`IsIntCSdone && aNbNodes > 1`) is
                        // UNPORTED for the same reason.
                        let _ = a_val_max;
                        a_si[k - 1].0 = pdeb - 1.0;
                        a_si[k].0 = param;
                    }
                }
            }

            for k in 1..=nbp as usize {
                let mut para = a_si[k - 1].0;
                if (para - pdeb) < eps_x || (pfin - para) < eps_x {
                    continue;
                }

                let mut dist = 0.0;
                if !func.eval(para, &mut dist) {
                    continue;
                }
                let dist = dist.abs();

                let mut an_indx: i32 = -1;
                let a_param = a_si[k - 1].0;
                if dist < maxdist {
                    // `gxx:700-706`: the `!IsIntCSdone` prefix always holds.
                    if (a_param - pdeb).abs() <= PCONFUSION
                        || (a_param - pfin).abs() <= PCONFUSION
                    {
                        an_indx = p_sol.point_state(a_si[k - 1].1);
                    }
                }

                // `gxx:708`.
                let mut a_pnt = if an_indx < 0 {
                    *func.last_computed_point()
                } else {
                    *func.valpoint(an_indx as usize).ok_or_else(|| {
                        "Contap_TheSearch::BoundedArc: missing point state".to_string()
                    })?
                };

                // `gxx:710-732`: pull the point onto the intersection line.
                if dist > 0.1 * CONFUSION {
                    let a_f_par = if k == 1 {
                        pdeb
                    } else {
                        (para + a_si[k - 2].0) / 2.0
                    };
                    let a_l_par = if k == nbp as usize {
                        pfin
                    } else {
                        (para + a_si[k].0) / 2.0
                    };
                    // `MinFunction` (`gxx:105-127`): `F(x)^2`, failure when the
                    // function cannot be evaluated.
                    let mut minf = |x: f64| {
                        let mut v = 0.0;
                        if func.eval(x, &mut v) {
                            Some(v * v)
                        } else {
                            None
                        }
                    };
                    let mut a_min = BrentMinimum::with_tolerance(CONFUSION);
                    a_min.perform_fallible(&mut minf, a_f_par, para, a_l_par);
                    if a_min.is_done() {
                        if let Some(loc) = a_min.location() {
                            para = loc;
                            if let Some(p2d) = sample::arc_value(arc, para) {
                                if let Some(s) = func.surface() {
                                    if let Some(p) = s.d0(p2d.x(), p2d.y()) {
                                        a_pnt = p;
                                    }
                                }
                            }
                        }
                    }
                }

                // `gxx:736`.
                let _range = point_process(domain, points, a_pnt, para, arc, tol_boundary)?;
            }
        }
    }

    // `gxx:744-752`.
    nbi = p_sol.nb_intervals();
    if !recheck_on_regularity && nbp > 0 {
        nbi = 0;
    }

    // `gxx:754-787`: one segment per solution interval.
    let mut pardeb = 0.0;
    let mut parfin = 0.0;
    for seg_i in 1..=nbi {
        let mut newseg = ContapTheSegmentOfTheSearch::new();
        newseg.set_value(arc.clone());
        let (pd, pf) = p_sol.interval(seg_i);
        pardeb = pd;
        parfin = pf;
        let (ideb, ifin) = p_sol.interval_state(seg_i);
        let ptdeb = *func
            .valpoint(ideb as usize)
            .ok_or_else(|| "Contap_TheSearch::BoundedArc: missing interval start".to_string())?;
        let ptfin = *func
            .valpoint(ifin as usize)
            .ok_or_else(|| "Contap_TheSearch::BoundedArc: missing interval end".to_string())?;

        let ranged = point_process(domain, points, ptdeb, pardeb, arc, tol_boundary)?;
        newseg.set_limit_point(&points[ranged - 1], true);
        let rangef = point_process(domain, points, ptfin, parfin, arc, tol_boundary)?;
        newseg.set_limit_point(&points[rangef - 1], false);
        segments.push(newseg);
    }

    // `gxx:789-794`.
    if nbi == 1
        && (pardeb - pdeb).abs() < PCONFUSION
        && (parfin - pfin).abs() < PCONFUSION
    {
        arcsol = true;
    }

    Ok(arcsol)
}

/// `ComputeBoundsfromInfinite(Func, PDeb, PFin, NbEchant)` (`gxx:802-877`).
fn compute_bounds_from_infinite(func: &mut ContapArcFunction) -> (f64, f64, i32) {
    let nb_echant = 100;

    let mut u0 = 0.0;
    let d_u = 0.001;
    let mut dist0 = 0.0;
    let mut dist1 = 0.0;

    func.eval(u0, &mut dist0);
    func.eval(u0 + d_u, &mut dist1);
    let mut d_dist = dist1 - dist0;

    if d_dist != 0.0 {
        u0 -= d_u * dist0 / d_dist;
        let pdeb;
        let pfin;

        let mut umin = u0 - 1.0e5;
        func.eval(umin, &mut dist0);
        func.eval(umin + d_u, &mut dist1);
        d_dist = dist1 - dist0;
        if d_dist != 0.0 {
            umin -= d_u * dist0 / d_dist;
        } else {
            umin -= 10.0;
        }

        let mut umax = u0 + 1.0e8;
        func.eval(umax, &mut dist0);
        func.eval(umax + d_u, &mut dist1);
        d_dist = dist1 - dist0;
        if d_dist != 0.0 {
            umax -= d_u * dist0 / d_dist;
        } else {
            umax += 10.0;
        }

        if umin > u0 {
            umin = u0 - 10.0;
        }
        if umax < u0 {
            umax = u0 + 10.0;
        }

        pfin = umax + 10.0 * (umax - umin);
        pdeb = umin - 10.0 * (umax - umin);
        (pdeb, pfin, nb_echant)
    } else {
        // `gxx:868-873`: possibly an arc entirely inside the quadric.
        (1.0e10, -1.0e10, nb_echant)
    }
}

/// `PointProcess(Pt, Para, A, Domain, pnt, Tol, Range)` (`gxx:880-1001`).
/// Answers `Range`, 1-based.
fn point_process(
    domain: &mut TopolTool,
    points: &mut Vec<ContapThePathPointOfTheSearch>,
    pt: GpPnt,
    para: f64,
    arc: &Arc<Curve2dAdaptor>,
    tol: f64,
) -> Result<usize, String> {
    let nbsol = points.len() as i32;

    domain.initialize_curve(arc);
    domain.init_vertex_iterator();
    let mut found = false;
    let mut goon = domain.more_vertex();
    let mut range = 0usize;

    while goon {
        let vtx = domain
            .vertex()
            .ok_or_else(|| "Contap_TheSearch::PointProcess: missing vertex".to_string())?;
        let dist = (para - sample::parameter(&vtx, arc)).abs();
        let toler = sample::tolerance(&vtx, arc);
        if dist <= toler {
            // `gxx:925-950`: look for the vertex in the solution list.
            let mut k = 1i32;
            let mut done = k > nbsol;
            while !done {
                let ptsol = &points[(k - 1) as usize];
                if !ptsol.is_new() {
                    let same_vertex = ptsol
                        .vertex()
                        .map(|v| TopolTool::identical(v, &vtx))
                        .unwrap_or(false);
                    if same_vertex
                        && ptsol.is_on_arc(arc)
                        && (ptsol.parameter() - para).abs() <= toler
                    {
                        done = true;
                    } else {
                        k += 1;
                        done = k > nbsol;
                    }
                } else {
                    k += 1;
                    done = k > nbsol;
                }
            }
            if k <= nbsol {
                // `gxx:952-954`: the vertex is already in the list.
                range = k as usize;
            } else {
                // `gxx:956-958`.
                points.push(ContapThePathPointOfTheSearch::on_vertex(
                    pt,
                    tol,
                    vtx,
                    arc.clone(),
                    para,
                ));
                range = points.len();
            }
            found = true;
            goon = false;
        } else {
            domain.next_vertex();
            goon = domain.more_vertex();
        }
    }

    if !found {
        // `gxx:971-987`: no vertex, so reuse an existing interior point.
        let mut found_internal = false;
        for k in 0..points.len() {
            let ptsol = &points[k];
            if !ptsol.is_on_arc(arc) || !ptsol.is_new() {
                continue;
            }
            if (ptsol.parameter() - para).abs() <= PCONFUSION {
                found_internal = true;
                range = k + 1;
            }
        }
        if !found_internal {
            // `gxx:990-998`.
            let mut local_tol = tol;
            local_tol *= 1000.0;
            if local_tol > 0.005 {
                local_tol = 0.005;
            }
            points.push(ContapThePathPointOfTheSearch::at(
                pt,
                local_tol,
                arc.clone(),
                para,
            ));
            range = points.len();
        }
    }

    Ok(range)
}

/// `IsRegularity(A, aDomain)` (`gxx:1005-1017`): `BRep_Tool::HasContinuity` of
/// the domain's current edge.
fn is_regularity(domain: &TopolTool) -> bool {
    domain.edge().map(brep_lib_continuity::has_continuity).unwrap_or(false)
}

/// `TreatLC(A, aDomain, aQuadric, TolBoundary, pnt)` (`gxx:1021-1125`).
///
/// Only the guards run for this instantiation: `Contap_ArcFunction::Quadric()`
/// never answers a cylinder (see the module note), so the `Standard_EXPORT`
/// body below the `GeomAbs_Cylinder` test (`gxx:1050-1121`) is UNPORTED and the
/// answer is always `1`, the function's own "no tangent line" code.
fn treat_lc(
    domain: &TopolTool,
    _tol_boundary: f64,
    _points: &mut Vec<ContapThePathPointOfTheSearch>,
) -> Result<i32, String> {
    let an_exit_code = 1;

    // `gxx:1028-1032`: no edge behind the domain.
    let Some(edge) = domain.edge() else {
        return Ok(an_exit_code);
    };

    // `gxx:1034-1038`: `BRep_Tool::Degenerated`.
    if occt_topo::brep_tool::BRepTool::is_degenerated(edge) {
        return Ok(an_exit_code);
    }

    // `gxx:1041-1048`: `BRepAdaptor_Curve::GetType() != GeomAbs_Line`.
    let is_line = occt_topo::brep_tool::BRepTool::edge_curve(edge)
        .map(|c| c.gp_line().is_some())
        .unwrap_or(false);
    if !is_line {
        return Ok(an_exit_code);
    }

    // `gxx:1050-1055`: `aQuadric.TypeQuadric() != GeomAbs_Cylinder` - always
    // true here, so `gxx:1057-1121` (the `Extrema_ExtCC` tangency test) is
    // UNPORTED.
    Ok(an_exit_code)
}
