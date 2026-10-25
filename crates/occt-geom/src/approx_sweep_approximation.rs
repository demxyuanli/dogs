//! `Approx_SweepApproximation` — approximation of the surface swept by a
//! section law (`Approx_SweepFunction`), with optional 2d curves.
//! Source: `ModelingData/TKGeomBase/Approx/Approx_SweepApproximation.cxx`
//! (constructor `cxx:60-68`, `Perform` `cxx:71-275`, `Approximation`
//! `cxx:283-408`, `Eval` `cxx:410-435`, `D0` `cxx:437-501`, `D1`
//! `cxx:503-583`, `D2` `cxx:585-678`, `SurfShape`/`Surface` `cxx:680-717`,
//! `MaxErrorOnSurf` `cxx:719-753`, `AverageErrorOnSurf` `cxx:755-783`,
//! `Curves2dShape` `cxx:785-798`, `Curve2d` `cxx:800-820`) and
//! `Approx_SweepApproximation.lxx`.
//!
//! Translation notes:
//! - `NCollection_Array2` fields are `Vec<Vec<_>>` indexed `[u][v]`, keeping
//!   OCCT's row = section pole, column = approximation pole layout
//!   (`cxx:319-320`).
//! - The mutable pole caches of `D0`/`D1`/`D2` are `RefCell`s so that the
//!   evaluator can be handed to `AdvApprox_ApproxAFunction` as an `Fn`.
//! - `StdFail_NotDone` throws of the accessors (`cxx:688-690`, `:706-708`,
//!   `:723-725`, `:758-760`, `:787-789`) become `Option::None`; the
//!   `Standard_DomainError` of `Curves2dShape` (`cxx:791-794`) likewise.

use std::cell::{Cell, RefCell};
use std::sync::Arc;

use occt_core::adv_approx::{ApproxAFunction, Cutting, DichoCutting, PrefAndRec};
use occt_core::approx::{ApproxSweepFunction, SweepD1, SweepD2};
use occt_core::gp::{GpGTrsf2d, GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::kernel::geomabs::Shape;

/// `Approx_SweepApproximation` (`Approx_SweepApproximation.hxx:50-227`).
pub struct ApproxSweepApproximation {
    func: Arc<dyn ApproxSweepFunction>,
    done: bool,
    num_1d_ss: i32,
    num_2d_ss: i32,
    num_3d_ss: i32,
    udeg: i32,
    vdeg: i32,
    deg_2d: i32,
    /// `tabPoles(1..Num3DSS, 1..NbPoles)` (`cxx:319`).
    tab_poles: Vec<Vec<GpPnt>>,
    /// `tabWeights(1..Num3DSS, 1..NbPoles)` (`cxx:320`).
    tab_weights: Vec<Vec<f64>>,
    tab_u_knots: Vec<f64>,
    tab_u_mults: Vec<i32>,
    tab_v_knots: Vec<f64>,
    tab_v_mults: Vec<i32>,
    tab_2d_knots: Vec<f64>,
    tab_2d_mults: Vec<i32>,
    seq_poles_2d: Vec<Vec<GpPnt2d>>,
    m_error_1d: Vec<f64>,
    a_error_1d: Vec<f64>,
    m_error_3d: Vec<f64>,
    a_error_3d: Vec<f64>,
    tab_2d_error: Vec<f64>,
    ave_2d_error: Vec<f64>,
    con_surf_err: Vec<f64>,
    /// `AAffin(1..Num2DSS)` (`cxx:144`).
    affin: Vec<GpGTrsf2d>,
    translation: GpVec,
}

/// Mutable caches and configuration of `D0`/`D1`/`D2` (`cxx:437-678`).
struct SweepEval<'a> {
    func: &'a dyn ApproxSweepFunction,
    num_1d_ss: i32,
    num_2d_ss: i32,
    num_3d_ss: i32,
    translation: GpVec,
    affin: Vec<GpGTrsf2d>,
    poles: RefCell<Vec<GpPnt>>,
    dpoles: RefCell<Vec<GpVec>>,
    d2poles: RefCell<Vec<GpVec>>,
    poles2d: RefCell<Vec<GpPnt2d>>,
    dpoles2d: RefCell<Vec<GpVec2d>>,
    d2poles2d: RefCell<Vec<GpVec2d>>,
    weights: RefCell<Vec<f64>>,
    dweights: RefCell<Vec<f64>>,
    d2weights: RefCell<Vec<f64>>,
    order: Cell<i32>,
    param: Cell<f64>,
    first: Cell<f64>,
    last: Cell<f64>,
}

impl<'a> SweepEval<'a> {
    /// Cache initialisation of the constructor (`cxx:62-67`).
    fn new(
        func: &'a dyn ApproxSweepFunction,
        num_1d_ss: i32,
        num_2d_ss: i32,
        num_3d_ss: i32,
        translation: GpVec,
        affin: Vec<GpGTrsf2d>,
    ) -> Self {
        Self {
            func,
            num_1d_ss,
            num_2d_ss,
            num_3d_ss,
            translation,
            affin,
            poles: RefCell::new(Vec::new()),
            dpoles: RefCell::new(Vec::new()),
            d2poles: RefCell::new(Vec::new()),
            poles2d: RefCell::new(Vec::new()),
            dpoles2d: RefCell::new(Vec::new()),
            d2poles2d: RefCell::new(Vec::new()),
            weights: RefCell::new(Vec::new()),
            dweights: RefCell::new(Vec::new()),
            d2weights: RefCell::new(Vec::new()),
            order: Cell::new(-1),
            param: Cell::new(0.0),
            first: Cell::new(1.0e100),
            last: Cell::new(-1.0e100),
        }
    }

    /// `Approx_SweepApproximation::Eval` (`cxx:410-435`): returns 0 on success,
    /// non-zero on failure.
    fn eval(&self, first: f64, last: f64, param: f64, derive: i32, out: &mut [f64]) -> i32 {
        let ok = match derive {
            0 => self.d0(param, first, last, out),
            1 => self.d1(param, first, last, out),
            2 => self.d2(param, first, last, out),
            _ => false,
        };
        if ok {
            0
        } else {
            1
        }
    }

    /// `Approx_SweepApproximation::D0` (`cxx:437-501`).
    fn d0(&self, param: f64, first: f64, last: f64, out: &mut [f64]) -> bool {
        if self.first.get() != first || last != self.last.get() {
            self.func.set_interval(first, last);
        }
        if param != self.param.get()
            || self.order.get() < 0
            || self.first.get() != first
            || last != self.last.get()
        {
            let (poles, poles2d, weights) = match self.func.d0(param, first, last) {
                Ok(v) => v,
                Err(_) => return false,
            };
            let mut mp = self.poles.borrow_mut();
            *mp = poles;
            {
                let mut mp2 = self.poles2d.borrow_mut();
                *mp2 = poles2d;
            }
            let mut mw = self.weights.borrow_mut();
            *mw = weights;
            // `cxx:461-467`: poles3d are multiplied by weight after translation.
            for ii in 0..self.num_1d_ss.max(0) as usize {
                mp[ii].coord = mp[ii].coord.subtracted(self.translation.xyz());
                mp[ii].coord = mp[ii].coord.multiply_scalar(mw[ii]);
            }
            // `cxx:469-471`: the transformation is applied to poles 2d.
            let mut mp2 = self.poles2d.borrow_mut();
            for ii in 0..self.num_2d_ss.max(0) as usize {
                self.affin[ii].transforms_xy(&mut mp2[ii].coord);
            }
            drop(mp2);
            self.first.set(first);
            self.last.set(last);
            self.order.set(0);
            self.param.set(param);
        }
        // `cxx:481-499`: extraction of results.
        let mp = self.poles.borrow();
        let mp2 = self.poles2d.borrow();
        let mw = self.weights.borrow();
        let mut index = 0usize;
        for ii in 0..self.num_1d_ss.max(0) as usize {
            out[index] = mw[ii];
            index += 1;
        }
        for ii in 0..self.num_2d_ss.max(0) as usize {
            out[index] = mp2[ii].x();
            out[index + 1] = mp2[ii].y();
            index += 2;
        }
        for ii in 0..self.num_3d_ss.max(0) as usize {
            out[index] = mp[ii].x();
            out[index + 1] = mp[ii].y();
            out[index + 2] = mp[ii].z();
            index += 3;
        }
        true
    }

    /// `Approx_SweepApproximation::D1` (`cxx:503-583`).
    fn d1(&self, param: f64, first: f64, last: f64, out: &mut [f64]) -> bool {
        if self.first.get() != first || last != self.last.get() {
            self.func.set_interval(first, last);
        }
        if param != self.param.get()
            || self.order.get() < 1
            || self.first.get() != first
            || last != self.last.get()
        {
            let d: SweepD1 = match self.func.d1(param, first, last) {
                Ok(v) => v,
                Err(_) => return false,
            };
            let mut mp = self.poles.borrow_mut();
            *mp = d.poles;
            let mut md = self.dpoles.borrow_mut();
            *md = d.dpoles;
            let mut mw = self.weights.borrow_mut();
            *mw = d.weights;
            let mut mdw = self.dweights.borrow_mut();
            *mdw = d.dweights;
            // `cxx:535-546`.
            for ii in 0..self.num_1d_ss.max(0) as usize {
                mp[ii].coord = mp[ii].coord.subtracted(self.translation.xyz());
                let a_weight = mw[ii];
                md[ii] = md[ii].multiplied_scalar(a_weight);
                let vaux = GpVec::from_xyz(&mp[ii].coord);
                md[ii] = md[ii].added(&vaux.multiplied_scalar(mdw[ii]));
                mp[ii].coord = mp[ii].coord.multiply_scalar(a_weight);
            }
            // `cxx:549-554`: transformation 2d on the suitable vectors.
            let mut mp2 = self.poles2d.borrow_mut();
            *mp2 = d.poles2d;
            let mut md2 = self.dpoles2d.borrow_mut();
            *md2 = d.dpoles2d;
            for ii in 0..self.num_2d_ss.max(0) as usize {
                let mut vcoord = md2[ii].coord;
                self.affin[ii].transforms_xy(&mut vcoord);
                md2[ii].coord = vcoord;
                self.affin[ii].transforms_xy(&mut mp2[ii].coord);
            }
            drop(md2);
            drop(mp2);
            self.first.set(first);
            self.last.set(last);
            self.order.set(1);
            self.param.set(param);
        }
        // `cxx:562-581`: extraction of results.
        let md = self.dpoles.borrow();
        let md2 = self.dpoles2d.borrow();
        let mdw = self.dweights.borrow();
        let mut index = 0usize;
        for ii in 0..self.num_1d_ss.max(0) as usize {
            out[index] = mdw[ii];
            index += 1;
        }
        for ii in 0..self.num_2d_ss.max(0) as usize {
            out[index] = md2[ii].x();
            out[index + 1] = md2[ii].y();
            index += 2;
        }
        for ii in 0..self.num_3d_ss.max(0) as usize {
            out[index] = md[ii].x();
            out[index + 1] = md[ii].y();
            out[index + 2] = md[ii].z();
            index += 3;
        }
        true
    }

    /// `Approx_SweepApproximation::D2` (`cxx:585-678`).
    fn d2(&self, param: f64, first: f64, last: f64, out: &mut [f64]) -> bool {
        if self.first.get() != first || last != self.last.get() {
            self.func.set_interval(first, last);
        }
        if param != self.param.get()
            || self.order.get() < 2
            || self.first.get() != first
            || last != self.last.get()
        {
            let d: SweepD2 = match self.func.d2(param, first, last) {
                Ok(v) => v,
                Err(_) => return false,
            };
            let mut mp = self.poles.borrow_mut();
            *mp = d.poles;
            let mut md = self.dpoles.borrow_mut();
            *md = d.dpoles;
            let mut md2p = self.d2poles.borrow_mut();
            *md2p = d.d2poles;
            let mut mw = self.weights.borrow_mut();
            *mw = d.weights;
            let mut mdw = self.dweights.borrow_mut();
            *mdw = d.dweights;
            let mut md2w = self.d2weights.borrow_mut();
            *md2w = d.d2weights;
            // `cxx:618-636`: multiply poles3d by the weight after translation.
            for ii in 0..self.num_1d_ss.max(0) as usize {
                mp[ii].coord = mp[ii].coord.subtracted(self.translation.xyz());
                // Second derivative.
                md2p[ii] = md2p[ii].multiplied_scalar(mw[ii]);
                let vaux = GpVec::from_xyz(md[ii].xyz());
                md2p[ii] = md2p[ii].added(&vaux.multiplied_scalar(2.0 * mdw[ii]));
                let vaux = GpVec::from_xyz(&mp[ii].coord);
                md2p[ii] = md2p[ii].added(&vaux.multiplied_scalar(md2w[ii]));
                // Then the remainder for the cache.
                md[ii] = md[ii].multiplied_scalar(mw[ii]);
                let vaux = GpVec::from_xyz(&mp[ii].coord);
                md[ii] = md[ii].added(&vaux.multiplied_scalar(mdw[ii]));
                mp[ii].coord = mp[ii].coord.multiply_scalar(mw[ii]);
            }
            // `cxx:638-648`: transformation to poles 2d.
            let mut mp2 = self.poles2d.borrow_mut();
            *mp2 = d.poles2d;
            let mut md2 = self.dpoles2d.borrow_mut();
            *md2 = d.dpoles2d;
            let mut md2p2 = self.d2poles2d.borrow_mut();
            *md2p2 = d.d2poles2d;
            for ii in 0..self.num_2d_ss.max(0) as usize {
                let mut vcoord = md2p2[ii].coord;
                self.affin[ii].transforms_xy(&mut vcoord);
                md2p2[ii].coord = vcoord;
                let mut vcoord = md2[ii].coord;
                self.affin[ii].transforms_xy(&mut vcoord);
                md2[ii].coord = vcoord;
                self.affin[ii].transforms_xy(&mut mp2[ii].coord);
            }
            drop(md2p2);
            drop(md2);
            drop(mp2);
            self.first.set(first);
            self.last.set(last);
            self.order.set(2);
            self.param.set(param);
        }
        // `cxx:656-676`: extraction of results.
        let md2p = self.d2poles.borrow();
        let md2p2 = self.d2poles2d.borrow();
        let md2w = self.d2weights.borrow();
        let mut index = 0usize;
        for ii in 0..self.num_1d_ss.max(0) as usize {
            out[index] = md2w[ii];
            index += 1;
        }
        for ii in 0..self.num_2d_ss.max(0) as usize {
            out[index] = md2p2[ii].x();
            out[index + 1] = md2p2[ii].y();
            index += 2;
        }
        for ii in 0..self.num_3d_ss.max(0) as usize {
            out[index] = md2p[ii].x();
            out[index + 1] = md2p[ii].y();
            out[index + 2] = md2p[ii].z();
            index += 3;
        }
        true
    }
}

impl ApproxSweepApproximation {
    /// `Approx_SweepApproximation(Func)` (`cxx:60-68`).
    pub fn new(func: Arc<dyn ApproxSweepFunction>) -> Self {
        Self {
            func,
            done: false,
            num_1d_ss: 0,
            num_2d_ss: 0,
            num_3d_ss: 0,
            udeg: 0,
            vdeg: 0,
            deg_2d: 0,
            tab_poles: Vec::new(),
            tab_weights: Vec::new(),
            tab_u_knots: Vec::new(),
            tab_u_mults: Vec::new(),
            tab_v_knots: Vec::new(),
            tab_v_mults: Vec::new(),
            tab_2d_knots: Vec::new(),
            tab_2d_mults: Vec::new(),
            seq_poles_2d: Vec::new(),
            m_error_1d: Vec::new(),
            a_error_1d: Vec::new(),
            m_error_3d: Vec::new(),
            a_error_3d: Vec::new(),
            tab_2d_error: Vec::new(),
            ave_2d_error: Vec::new(),
            con_surf_err: Vec::new(),
            affin: Vec::new(),
            translation: GpVec::zero(),
        }
    }

    /// `Approx_SweepApproximation::Perform` (`cxx:71-275`).
    #[allow(clippy::too_many_arguments)]
    pub fn perform(
        &mut self,
        first: f64,
        last: f64,
        tol3d: f64,
        bound_tol: f64,
        tol2d: f64,
        tol_angular: f64,
        continuity: Shape,
        degmax: i32,
        segmax: i32,
    ) {
        // (1) Characteristics of a section (`cxx:82-91`).
        let (nb_pol_sect, _nb_knot_sect, udeg) = self.func.section_shape();
        self.udeg = udeg;
        self.num_2d_ss = self.func.nb_2d_curves();
        self.tab_u_knots = self.func.knots();
        self.tab_u_mults = self.func.mults();

        // (2) Decomposition into sub-spaces (`cxx:93-130`).
        self.num_3d_ss = nb_pol_sect;
        let mut one_d_tol = vec![0.0f64; self.num_3d_ss.max(0) as usize];
        let mut three_d_tol = self.func.get_tolerance(bound_tol, tol3d, tol_angular);
        three_d_tol.resize(self.num_3d_ss.max(0) as usize, 0.0);
        let mut tol3d_min = tol3d;
        for ii in 0..self.num_3d_ss.max(0) as usize {
            if three_d_tol[ii] < tol3d_min {
                tol3d_min = three_d_tol[ii];
            }
        }
        if self.func.is_rational() {
            self.num_1d_ss = nb_pol_sect;
            let wmin = self.func.get_minimal_weight().unwrap_or_default();
            let size = self.func.maximal_section().unwrap_or(0.0);
            self.translation = self.func.barycentre_of_surf().map(|p| GpVec::from_xyz(&p.coord)).unwrap_or_else(|_| GpVec::zero());
            for ii in 0..self.num_3d_ss.max(0) as usize {
                // To take account of the error on the final result.
                let mut tol = three_d_tol[ii] / 2.0;
                one_d_tol[ii] = tol * wmin.get(ii).copied().unwrap_or(1.0) / size;
                tol *= wmin.get(ii).copied().unwrap_or(1.0); // Factor of projection.
                three_d_tol[ii] = tol.max(1.0e-20);
            }
        } else {
            self.num_1d_ss = 0;
        }

        // (2.2) Tolerance and transformation 2d (`cxx:132-160`).
        let mut two_d_tol: Vec<f64> = Vec::new();
        let mut the_3d2d_tol = 0.0f64;
        if self.num_2d_ss > 0 {
            self.affin = vec![GpGTrsf2d::identity(); self.num_2d_ss as usize];
            the_3d2d_tol = 0.9 * bound_tol; // 10% of security.
            for ii in 0..self.num_2d_ss as usize {
                let (tolu, tolv) = self
                    .func
                    .resolution(ii as i32 + 1, the_3d2d_tol)
                    .unwrap_or((0.0, 0.0));
                let res = if tolu > tolv {
                    self.affin[ii].set_value(1, 1, tolv / tolu);
                    tolv
                } else {
                    self.affin[ii].set_value(2, 2, tolu / tolv);
                    tolu
                };
                two_d_tol.push(tol2d.min(res));
            }
        }

        // (3) Check whether D2 / D1 are implemented (`cxx:188-226`).
        let mut cont = continuity;
        if cont >= Shape::C2 && self.func.d2(first, first, last).is_err() {
            cont = Shape::C1;
        }
        if cont == Shape::C1 && self.func.d1(first, first, last).is_err() {
            cont = Shape::C0;
        }
        // So that F is at least 20 times more exact than its approximation.
        self.func.set_tolerance(tol3d_min / 20.0, tol2d / 20.0);

        let _nb_interval_c2 = self.func.nb_intervals(Shape::C2 as u8);
        let nb_interval_c3 = self.func.nb_intervals(Shape::C3 as u8);

        let eval = SweepEval::new(
            &*self.func,
            self.num_1d_ss,
            self.num_2d_ss,
            self.num_3d_ss,
            self.translation,
            self.affin.clone(),
        );
        let callback = |f: f64, l: f64, p: f64, d: i32, out: &mut [f64]| eval.eval(f, l, p, d, out);

        let approx = if nb_interval_c3 > 1 {
            // (3.1) Approximation with preferential cut (`cxx:236-256`).
            let param_c2 = self.func.intervals(Shape::C2 as u8);
            let param_c3 = self.func.intervals(Shape::C3 as u8);
            match PrefAndRec::new(param_c2, param_c3, 5.0) {
                Ok(cut) => ApproxAFunction::approx(
                    self.num_1d_ss,
                    self.num_2d_ss,
                    self.num_3d_ss,
                    &one_d_tol,
                    &two_d_tol,
                    &three_d_tol,
                    first,
                    last,
                    cont,
                    degmax,
                    segmax,
                    &cut,
                    &callback,
                ),
                Err(_) => Ok(ApproxAFunction::default()),
            }
        } else {
            // (3.2) Approximation without preferential cut (`cxx:258-272`).
            let cut = DichoCutting;
            ApproxAFunction::approx(
                self.num_1d_ss,
                self.num_2d_ss,
                self.num_3d_ss,
                &one_d_tol,
                &two_d_tol,
                &three_d_tol,
                first,
                last,
                cont,
                degmax,
                segmax,
                &cut as &dyn Cutting,
                &callback,
            )
        };

        let approx = match approx {
            Ok(a) => a,
            Err(_) => {
                self.done = false;
                return;
            }
        };
        self.done = approx.has_result;
        if !self.done {
            return;
        }

        // --> Fill the fields of the surface (`cxx:311-408`).
        self.vdeg = approx.degree;
        let nb_poles = approx.nb_poles();
        let mut tab_poles = vec![vec![GpPnt::zero(); nb_poles]; self.num_3d_ss.max(0) as usize];
        let mut tab_weights = vec![vec![0.0f64; nb_poles]; self.num_3d_ss.max(0) as usize];
        if self.num_1d_ss == self.num_3d_ss {
            for ii in 0..self.num_3d_ss.max(0) as usize {
                let poles = approx.poles(ii);
                let weights = approx.poles_1d(ii);
                for jj in 0..nb_poles {
                    let wpoid = weights[jj];
                    let mut p = poles[jj];
                    // It is necessary to divide poles by weight.
                    p.coord = p.coord.multiply_scalar(1.0 / wpoid);
                    p.coord = p.coord.added(self.translation.xyz());
                    tab_poles[ii][jj] = p;
                    tab_weights[ii][jj] = wpoid;
                }
            }
        } else {
            for ii in 0..self.num_3d_ss.max(0) as usize {
                let poles = approx.poles(ii);
                for jj in 0..nb_poles {
                    tab_poles[ii][jj] = poles[jj];
                    tab_weights[ii][jj] = 1.0;
                }
            }
        }
        self.tab_poles = tab_poles;
        self.tab_weights = tab_weights;
        self.tab_v_knots = approx.knots.clone();
        self.tab_v_mults = approx.mults.clone();

        // --> Filling of curves 2D (`cxx:355-373`).
        if self.num_2d_ss > 0 {
            self.deg_2d = self.vdeg;
            self.tab_2d_knots = approx.knots.clone();
            self.tab_2d_mults = approx.mults.clone();
            self.seq_poles_2d.clear();
            for ii in 0..self.num_2d_ss as usize {
                let trsf_inv = match self.affin[ii].inverted() {
                    Ok(t) => t,
                    Err(_) => GpGTrsf2d::identity(),
                };
                let mut p2d = approx.poles_2d(ii).to_vec();
                // Do not forget to apply inverted homothety.
                for p in p2d.iter_mut() {
                    trsf_inv.transforms_xy(&mut p.coord);
                }
                self.seq_poles_2d.push(p2d);
            }
        }

        // ---> Filling of errors (`cxx:376-407`).
        self.m_error_3d = (0..self.num_3d_ss.max(0) as usize)
            .map(|ii| approx.max_error(3, ii))
            .collect();
        self.a_error_3d = (0..self.num_3d_ss.max(0) as usize)
            .map(|ii| approx.average_error(3, ii))
            .collect();
        if self.func.is_rational() {
            self.m_error_1d = (0..self.num_1d_ss.max(0) as usize)
                .map(|ii| approx.max_error(1, ii))
                .collect();
            self.a_error_1d = (0..self.num_1d_ss.max(0) as usize)
                .map(|ii| approx.average_error(1, ii))
                .collect();
        }
        if self.num_2d_ss > 0 {
            self.tab_2d_error = (0..self.num_2d_ss as usize)
                .map(|ii| approx.max_error(2, ii))
                .collect();
            self.ave_2d_error = (0..self.num_2d_ss as usize)
                .map(|ii| approx.average_error(2, ii))
                .collect();
            self.con_surf_err = (0..self.num_2d_ss as usize)
                .map(|ii| (self.tab_2d_error[ii] / two_d_tol[ii]) * bound_tol)
                .collect();
        }
    }

    /// `Approx_SweepApproximation::IsDone` (`lxx`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `UDegree()` (`hxx:108`).
    pub fn u_degree(&self) -> i32 {
        self.udeg
    }

    /// `VDegree()` (`hxx:110`).
    pub fn v_degree(&self) -> i32 {
        self.vdeg
    }

    /// `SurfShape` (`cxx:680-696`): `(u_degree, v_degree, nb_u_poles,
    /// nb_v_poles, nb_u_knots, nb_v_knots)`. `nb_u_poles` is OCCT's
    /// `tabPoles->ColLength()` (poles inside one section) and `nb_v_poles` its
    /// `RowLength()` (number of sections) (`cxx:692-693`).
    pub fn surf_shape(&self) -> Option<(i32, i32, usize, usize, usize, usize)> {
        if !self.done {
            return None;
        }
        let nb_u = self.tab_poles.first().map(|v| v.len()).unwrap_or(0);
        Some((
            self.udeg,
            self.vdeg,
            nb_u,
            self.tab_poles.len(),
            self.tab_u_knots.len(),
            self.tab_v_knots.len(),
        ))
    }

    /// `SurfPoles()` (`hxx:112`), indexed like OCCT's `tabPoles`: row = section
    /// (`v`, 0-based section index), column = pole (`u`) (`cxx:319-339`).
    pub fn surf_poles(&self) -> &[Vec<GpPnt>] {
        &self.tab_poles
    }

    /// `SurfWeights()` (`hxx:114`), same layout as [`Self::surf_poles`].
    pub fn surf_weights(&self) -> &[Vec<f64>] {
        &self.tab_weights
    }

    /// `SurfUKnots()` (`hxx:116`).
    pub fn surf_u_knots(&self) -> &[f64] {
        &self.tab_u_knots
    }

    /// `SurfVKnots()` (`hxx:118`).
    pub fn surf_v_knots(&self) -> &[f64] {
        &self.tab_v_knots
    }

    /// `SurfUMults()` (`hxx:120`).
    pub fn surf_u_mults(&self) -> &[i32] {
        &self.tab_u_mults
    }

    /// `SurfVMults()` (`hxx:122`).
    pub fn surf_v_mults(&self) -> &[i32] {
        &self.tab_v_mults
    }

    /// `MaxErrorOnSurf()` (`cxx:719-753`).
    pub fn max_error_on_surf(&self) -> Option<f64> {
        if !self.done {
            return None;
        }
        if self.func.is_rational() {
            let wmin = self.func.get_minimal_weight().unwrap_or_default();
            let size = self.func.maximal_section().unwrap_or(0.0);
            let mut max_error = 0.0f64;
            for ii in 0..self.num_3d_ss.max(0) as usize {
                let w = wmin.get(ii).copied().unwrap_or(1.0);
                let err = (size * self.m_error_1d[ii] + self.m_error_3d[ii]) / w;
                if err > max_error {
                    max_error = err;
                }
            }
            Some(max_error)
        } else {
            let mut max_error = 0.0f64;
            for ii in 0..self.num_3d_ss.max(0) as usize {
                if self.m_error_3d[ii] > max_error {
                    max_error = self.m_error_3d[ii];
                }
            }
            Some(max_error)
        }
    }

    /// `AverageErrorOnSurf()` (`cxx:755-783`).
    pub fn average_error_on_surf(&self) -> Option<f64> {
        if !self.done {
            return None;
        }
        let mut moy_error = 0.0f64;
        if self.func.is_rational() {
            let wmin = self.func.get_minimal_weight().unwrap_or_default();
            let size = self.func.maximal_section().unwrap_or(0.0);
            for ii in 0..self.num_3d_ss.max(0) as usize {
                let w = wmin.get(ii).copied().unwrap_or(1.0);
                moy_error += (size * self.a_error_1d[ii] + self.a_error_3d[ii]) / w;
            }
        } else {
            for ii in 0..self.num_3d_ss.max(0) as usize {
                moy_error += self.a_error_3d[ii];
            }
        }
        Some(moy_error / self.num_3d_ss.max(1) as f64)
    }

    /// `NbCurves2d()` (`hxx:128`).
    pub fn nb_curves_2d(&self) -> usize {
        self.seq_poles_2d.len()
    }

    /// `Curves2dShape()` (`cxx:785-798`).
    pub fn curves_2d_shape(&self) -> Option<(i32, usize, usize)> {
        if !self.done || self.seq_poles_2d.is_empty() {
            return None;
        }
        Some((self.deg_2d, self.seq_poles_2d[0].len(), self.tab_2d_knots.len()))
    }

    /// `Curve2d(Index, ...)` (`cxx:800-820`): the poles of one 2d curve.
    pub fn curve_2d(&self, index: usize) -> Option<&[GpPnt2d]> {
        if !self.done {
            return None;
        }
        self.seq_poles_2d.get(index).map(|v| v.as_slice())
    }

    /// `Curves2dDegree()` (`hxx:136`).
    pub fn curves_2d_degree(&self) -> i32 {
        self.deg_2d
    }

    /// `Curve2dPoles(Index)` (`hxx:139`).
    pub fn curve_2d_poles(&self, index: usize) -> Option<&[GpPnt2d]> {
        self.seq_poles_2d.get(index).map(|v| v.as_slice())
    }

    /// `Curves2dKnots()` (`hxx:141`).
    pub fn curves_2d_knots(&self) -> &[f64] {
        &self.tab_2d_knots
    }

    /// `Curves2dMults()` (`hxx:143`).
    pub fn curves_2d_mults(&self) -> &[i32] {
        &self.tab_2d_mults
    }

    /// `Max2dError(Index)` (`hxx:147-148`).
    pub fn max_2d_error(&self, index: usize) -> f64 {
        self.tab_2d_error[index]
    }

    /// `Average2dError(Index)` (`hxx:150-151`).
    pub fn average_2d_error(&self, index: usize) -> f64 {
        self.ave_2d_error[index]
    }

    /// `TolCurveOnSurf(Index)` (`hxx:153-155`).
    pub fn tol_curve_on_surf(&self, index: usize) -> f64 {
        self.con_surf_err[index]
    }
}
