//! `Bisector_BisecPC` (`Bisector_BisecPC.cxx:36-884`, `.hxx:49-218`): the
//! algorithmic bisector between a curve and a point.
//!
//! The point at parameter `u` is obtained from the corresponding point `P(u)`
//! on the support curve by
//!
//! ```text
//!                1   ||PP(u)||**2
//! F(u) = P(u) - --- * ----------- * N(u)
//!                2    (N(u)|PP(u))
//! ```
//!
//! where `N(u)` is the normal of the support curve. The class also builds the
//! tangent extensions at the ends when the support curve is concave, and it
//! limits the definition range through `ComputeIntervals`.

use std::sync::Arc;

use occt_core::elib::clib2d::{line_parameter_ax2d, line_value_ax2d};
use occt_core::gp::{GpAx2d, GpDir2d, GpPnt2d, GpTrsf2d, GpVec2d};
use occt_core::precision::{CONFUSION, INFINITE, PCONFUSION, REAL_SMALL};

use crate::curve::Curve2d;
use crate::geom2d_api_project_point_on_curve::ProjectPointOnCurve2d;

use super::{curve_is_closed, is_convex, null_curve, BisectorBisecAnaError, BisectorCurve};

/// `gp_Dir2d(X, Y)` (`gp_Dir2d.hxx:304-320`): OCCT raises
/// `Standard_ConstructionError` on a null vector, the port returns the zero
/// direction instead of throwing.
fn dir2d_or_zero(x: f64, y: f64) -> GpDir2d {
    GpDir2d::new(x, y).unwrap_or(GpDir2d { x: 0.0, y: 0.0 })
}

/// `class Bisector_BisecPC` (`Bisector_BisecPC.hxx:49-218`).
pub struct BisectorBisecPC {
    /// `curve` (`.hxx:198`): the support curve, owned as a copy.
    curve: Arc<dyn Curve2d>,
    /// `point` (`.hxx:199`).
    point: GpPnt2d,
    /// `sign` (`.hxx:200`).
    sign: f64,
    /// `startIntervals` (`.hxx:201`), 1-based like OCCT: entry `k` is
    /// `startIntervals(k)`.
    start_intervals: Vec<f64>,
    /// `endIntervals` (`.hxx:202`).
    end_intervals: Vec<f64>,
    /// `bisInterval` (`.hxx:203`).
    bis_interval: usize,
    /// `currentInterval` (`.hxx:204`).
    current_interval: i32,
    /// `shiftParameter` (`.hxx:205`).
    shift_parameter: f64,
    /// `distMax` (`.hxx:206`).
    dist_max: f64,
    /// `isEmpty` (`.hxx:207`).
    is_empty: bool,
    /// `isConvex` (`.hxx:208`).
    is_convex: bool,
    /// `pointStartBis` (`.hxx:209`).
    point_start_bis: GpPnt2d,
    /// `pointEndBis` (`.hxx:210`).
    point_end_bis: GpPnt2d,
    /// `extensionStart` (`.hxx:211`).
    extension_start: bool,
    /// `extensionEnd` (`.hxx:212`).
    extension_end: bool,
}

impl Default for BisectorBisecPC {
    fn default() -> Self {
        Self::new()
    }
}

impl BisectorBisecPC {
    /// `Bisector_BisecPC()` (`cxx:36-50`).
    pub fn new() -> Self {
        Self {
            curve: null_curve(),
            point: GpPnt2d::zero(),
            sign: 0.0,
            start_intervals: Vec::new(),
            end_intervals: Vec::new(),
            bis_interval: 0,
            current_interval: 0,
            shift_parameter: 0.0,
            dist_max: 0.0,
            is_empty: true,
            is_convex: false,
            point_start_bis: GpPnt2d::zero(),
            point_end_bis: GpPnt2d::zero(),
            extension_start: false,
            extension_end: false,
        }
    }

    /// `Bisector_BisecPC(Cu, P, Side, DistMax)` (`cxx:54-58`).
    pub fn with_dist_max(cu: &dyn Curve2d, p: &GpPnt2d, side: f64, dist_max: f64) -> Self {
        let mut b = Self::new();
        b.perform(cu, p, side, dist_max);
        b
    }

    /// `Bisector_BisecPC(Cu, P, Side, UMin, UMax)` (`cxx:64-82`).
    ///
    /// UNPORTED detail: OCCT leaves `distMax`, `isEmpty` and
    /// `currentInterval` uninitialized in this constructor (they are plain
    /// members and the default constructor is not chained), so reading them is
    /// undefined behaviour there. The port stores `distMax = 0`,
    /// `is_empty = false` and `current_interval = 0`; no caller of this
    /// constructor reads `IsEmpty` or `Distance`.
    pub fn from_interval(
        cu: &dyn Curve2d,
        p: &GpPnt2d,
        side: f64,
        u_min: f64,
        u_max: f64,
    ) -> Self {
        let mut b = Self::new();
        b.curve = Arc::from(cu.clone_dyn());
        b.point = *p;
        b.sign = side;
        b.start_intervals.push(u_min);
        b.end_intervals.push(u_max);
        b.bis_interval = 1;
        b.extension_start = false;
        b.extension_end = false;
        b.point_start_bis = b.eval_d0(u_min);
        b.point_end_bis = b.eval_d0(u_max);
        b.is_convex = is_convex(&*b.curve, b.sign);
        b.is_empty = false;
        b
    }

    /// `Perform(Cu, P, Side, DistMax)` (`cxx:88-132`).
    pub fn perform(&mut self, cu: &dyn Curve2d, p: &GpPnt2d, side: f64, dist_max: f64) {
        self.curve = Arc::from(cu.clone_dyn());
        self.point = *p;
        self.dist_max = dist_max;
        self.sign = side;
        self.is_convex = is_convex(&*self.curve, self.sign);

        self.compute_intervals();
        if self.is_empty {
            return;
        }

        self.bis_interval = 1;
        self.extension_start = false;
        self.extension_end = false;
        self.point_start_bis = self.eval_d0(self.start_first());
        self.point_end_bis = self.eval_d0(self.end_last());

        if !self.is_convex {
            if self
                .point
                .is_equal(&self.curve.value(self.curve.first_parameter()), CONFUSION)
            {
                self.extension_start = true;
                let u_first = self.start_first() - p.distance(&self.point_start_bis);
                self.start_intervals.insert(0, u_first);
                let second = self.start_value(2);
                self.end_intervals.insert(0, second);
                self.bis_interval = 2;
            } else if self
                .point
                .is_equal(&self.curve.value(self.curve.last_parameter()), CONFUSION)
            {
                self.extension_end = true;
                let u_last = self.end_last() + p.distance(&self.point_end_bis);
                self.start_intervals.push(self.end_last());
                self.end_intervals.push(u_last);
                self.bis_interval = 1;
            }
        }
    }

    // --- interval sequence helpers (1-based, like NCollection_Sequence) ------

    fn start_first(&self) -> f64 {
        *self.start_intervals.first().unwrap_or(&0.0)
    }

    fn start_value(&self, k: usize) -> f64 {
        self.start_intervals[k - 1]
    }

    fn end_last(&self) -> f64 {
        *self.end_intervals.last().unwrap_or(&0.0)
    }

    fn end_value(&self, k: usize) -> f64 {
        self.end_intervals[k - 1]
    }

    /// `IsExtendAtStart` (`cxx:136-139`).
    pub fn is_extend_at_start(&self) -> bool {
        self.extension_start
    }

    /// `IsExtendAtEnd` (`cxx:143-146`).
    pub fn is_extend_at_end(&self) -> bool {
        self.extension_end
    }

    /// `ReversedParameter(U)` (`cxx:158-161`).
    pub fn reversed_parameter(&self, u: f64) -> f64 {
        self.last_parameter() + self.first_parameter() - u
    }

    /// `Copy()` (`cxx:165-184`).
    pub fn copy_bisector(&self) -> Self {
        let mut c = Self::new();
        c.curve = Arc::from(self.curve.clone_dyn());
        c.point = self.point;
        c.sign = self.sign;
        c.start_intervals = self.start_intervals.clone();
        c.end_intervals = self.end_intervals.clone();
        c.bis_interval = self.bis_interval;
        c.current_interval = self.current_interval;
        c.shift_parameter = self.shift_parameter;
        c.dist_max = self.dist_max;
        c.is_empty = self.is_empty;
        c.is_convex = self.is_convex;
        c.extension_start = self.extension_start;
        c.extension_end = self.extension_end;
        c.point_start_bis = self.point_start_bis;
        c.point_end_bis = self.point_end_bis;
        c
    }

    /// `Transform(T)` (`cxx:188-195`).
    pub fn transform_curve(&mut self, t: &GpTrsf2d) {
        let mut c = self.curve.clone_dyn();
        c.transform(t);
        self.curve = Arc::from(c);
        self.point.transform(t);
        self.point_start_bis.transform(t);
        self.point_end_bis.transform(t);
    }

    /// `IsCN(N)` (`cxx:199-202`): `curve->IsCN(N + 1)`.
    pub fn is_cn(&self, n: i32) -> bool {
        self.curve.continuity() as i32 >= n + 1
    }

    /// `FirstParameter()` (`cxx:206-209`); `0` when the sequence is empty
    /// (OCCT raises `Standard_NullObject` there, but every caller checks
    /// `IsEmpty` first).
    pub fn first_parameter(&self) -> f64 {
        self.start_first()
    }

    /// `LastParameter()` (`cxx:213-216`); `0` when the sequence is empty, see
    /// [`Self::first_parameter`].
    pub fn last_parameter(&self) -> f64 {
        self.end_last()
    }

    /// `Continuity()` (`cxx:220-238`).
    pub fn continuity(&self) -> u8 {
        match self.curve.continuity() {
            2 => 0, // GeomAbs_C1 -> GeomAbs_C0
            4 => 2, // GeomAbs_C2 -> GeomAbs_C1
            5 => 4, // GeomAbs_C3 -> GeomAbs_C2
            6 => 6, // GeomAbs_CN -> GeomAbs_CN
            _ => 0,
        }
    }

    /// `NbIntervals()` (`cxx:242-245`).
    pub fn nb_intervals(&self) -> i32 {
        self.start_intervals.len() as i32
    }

    /// `IntervalFirst(I)` (`cxx:249-252`).
    pub fn interval_first(&self, i: i32) -> Result<f64, BisectorBisecAnaError> {
        if i < 1 || (i as usize) > self.start_intervals.len() {
            return Err(BisectorBisecAnaError::OutOfRange);
        }
        Ok(self.start_value(i as usize))
    }

    /// `IntervalLast(I)` (`cxx:256-259`).
    pub fn interval_last(&self, i: i32) -> Result<f64, BisectorBisecAnaError> {
        if i < 1 || (i as usize) > self.end_intervals.len() {
            return Err(BisectorBisecAnaError::OutOfRange);
        }
        Ok(self.end_value(i as usize))
    }

    /// `IntervalContinuity()` (`cxx:263-279`).
    pub fn interval_continuity(&self) -> u8 {
        self.continuity()
    }

    /// `IsClosed()` (`cxx:283-298`).
    pub fn is_closed(&self) -> bool {
        if curve_is_closed(&*self.curve) {
            if self.start_first() == self.curve.first_parameter()
                && self.end_value(1) == self.curve.last_parameter()
            {
                return true;
            }
        }
        false
    }

    /// `IsEmpty()` (`cxx:769-772`).
    pub fn is_empty(&self) -> bool {
        self.is_empty
    }

    /// `LinkBisCurve(U)` (`cxx:753-756`).
    pub fn link_bis_curve(&self, u: f64) -> f64 {
        u - self.shift_parameter
    }

    /// `LinkCurveBis(U)` (`cxx:760-763`).
    pub fn link_curve_bis(&self, u: f64) -> f64 {
        u + self.shift_parameter
    }

    /// `Extension(U, P, V1, V2, V3)` (`cxx:308-349`).
    fn extension(&self, u: f64, p: &mut GpPnt2d, v1: &mut GpVec2d, v2: &mut GpVec2d, v3: &mut GpVec2d) {
        v1.set_coord(0.0, 0.0);
        v2.set_coord(0.0, 0.0);
        v3.set_coord(0.0, 0.0);
        if u < self.start_value(self.bis_interval) {
            if self.point_start_bis.is_equal(&self.point, PCONFUSION) {
                *p = self.point_start_bis;
            } else {
                let d_u = u - self.start_value(self.bis_interval);
                let dir_ext = dir2d_or_zero(
                    self.point_start_bis.x() - self.point.x(),
                    self.point_start_bis.y() - self.point.y(),
                );
                p.set_coord(
                    self.point_start_bis.x() + d_u * dir_ext.x(),
                    self.point_start_bis.y() + d_u * dir_ext.y(),
                );
                v1.set_coord(dir_ext.x(), dir_ext.y());
            }
        } else if u > self.end_value(self.bis_interval) {
            if self.point_end_bis.is_equal(&self.point, PCONFUSION) {
                *p = self.point_end_bis;
            } else {
                let d_u = u - self.end_value(self.bis_interval);
                let dir_ext = dir2d_or_zero(
                    self.point.x() - self.point_end_bis.x(),
                    self.point.y() - self.point_end_bis.y(),
                );
                p.set_coord(
                    self.point_end_bis.x() + d_u * dir_ext.x(),
                    self.point_end_bis.y() + d_u * dir_ext.y(),
                );
                v1.set_coord(dir_ext.x(), dir_ext.y());
            }
        }
    }

    /// `Values(U, N, P, V1, V2, V3)` (`cxx:363-448`).
    ///
    /// UNPORTED detail: the OCCT `switch (N)` has arms for `N = 0, 1, 2` only
    /// (`cxx:386-396`), so `N = 3` (reached from `EvalD3`) leaves `PC`, `Tu`,
    /// `Tuu` and `T3u` uninitialized - undefined behaviour there. The port
    /// evaluates `D3` for every `N >= 2`, which is what the `N = 2` arm does and
    /// what `EvalD3` needs.
    pub fn values(
        &self,
        u: f64,
        n: i32,
        p: &mut GpPnt2d,
        v1: &mut GpVec2d,
        v2: &mut GpVec2d,
        v3: &mut GpVec2d,
    ) {
        if u < self.start_value(self.bis_interval) {
            self.extension(u, p, v1, v2, v3);
            return;
        } else if u > self.end_value(self.bis_interval) {
            self.extension(u, p, v1, v2, v3);
            return;
        }
        let u_on_curve = self.link_bis_curve(u);

        let mut tu = GpVec2d::zero();
        let mut tuu = GpVec2d::zero();
        let mut t3u = GpVec2d::zero();
        let pc;

        if n == 0 {
            let (pt, d1) = self.curve.d1(u_on_curve);
            pc = pt;
            tu = d1;
        } else if n == 1 {
            let (pt, d1, d2) = self.curve.d2(u_on_curve);
            pc = pt;
            tu = d1;
            tuu = d2;
        } else {
            let (pt, d1, d2, d3) = self.curve.d3(u_on_curve);
            pc = pt;
            tu = d1;
            tuu = d2;
            t3u = d3;
        }

        let a_ppc = GpVec2d::new(pc.x() - self.point.x(), pc.y() - self.point.y());
        let nor = GpVec2d::new(-tu.y(), tu.x());

        let square_ppc = a_ppc.square_magnitude();
        let nor_ppc = nor.dot(&a_ppc);
        let a1;

        if nor_ppc.abs() > REAL_SMALL && (nor_ppc * self.sign) < 0.0 {
            a1 = 0.5 * square_ppc / nor_ppc;
            p.set_coord(pc.x() - nor.x() * a1, pc.y() - nor.y() * a1);
        } else {
            return;
        }

        if n == 0 {
            return;
        }

        let nu = GpVec2d::new(-tuu.y(), tuu.x());
        let nu_ppc = nu.dot(&a_ppc);
        let tu_ppc = tu.dot(&a_ppc);
        let nor_ppce2 = nor_ppc * nor_ppc;
        let a2 = tu_ppc / nor_ppc - 0.5 * nu_ppc * square_ppc / nor_ppce2;

        v1.set_coord(
            tu.x() - a1 * nu.x() - a2 * nor.x(),
            tu.y() - a1 * nu.y() - a2 * nor.y(),
        );

        if n == 1 {
            return;
        }

        let nuu = GpVec2d::new(-t3u.y(), t3u.x());
        let nor_ppce4 = nor_ppce2 * nor_ppce2;
        let nuu_ppc = nuu.dot(&a_ppc);
        let tuu_ppc = tuu.dot(&a_ppc);

        let a21 = tuu_ppc / nor_ppc - tu_ppc * nu_ppc / nor_ppce2;
        let a22 = (0.5 * nuu_ppc * square_ppc + nu_ppc * tu_ppc) / nor_ppce2
            - nu_ppc * square_ppc * nor_ppc * nu_ppc / nor_ppce4;
        let a2u = a21 - a22;

        v2.set_coord(
            tuu.x() - 2.0 * a2 * nu.x() - a1 * nuu.x() - a2u * nor.x(),
            tuu.y() - 2.0 * a2 * nu.y() - a1 * nuu.y() - a2u * nor.y(),
        );
    }

    /// `Distance(U)` (`cxx:486-541`): the **squared** distance from the bisector
    /// point at `U` to the support curve, or `Precision::Infinite()` when the
    /// point is outside the bisector's domain.
    pub fn distance(&self, u: f64) -> f64 {
        let u_on_curve = self.link_bis_curve(u);
        let (pc, tan) = self.curve.d1(u_on_curve);
        let a_ppc = GpVec2d::new(pc.x() - self.point.x(), pc.y() - self.point.y());
        let nor = GpVec2d::new(-tan.y(), tan.x());

        let nor_nor = nor.square_magnitude();
        let square_mag_ppc = a_ppc.square_magnitude();
        let prosca = nor.dot(&a_ppc);

        if self.point.is_equal(&pc, CONFUSION) {
            if self.is_convex {
                return 0.0;
            } else {
                return INFINITE;
            }
        }

        if prosca.abs() < CONFUSION || (prosca * self.sign) > 0.0 {
            INFINITE
        } else {
            let a = 0.5 * square_mag_ppc / prosca;
            a * a * nor_nor
        }
    }

    /// `EvalD0(U)` (`cxx:545-551`).
    pub fn eval_d0(&self, u: f64) -> GpPnt2d {
        let mut p = self.point;
        let mut v1 = GpVec2d::zero();
        let mut v2 = GpVec2d::zero();
        let mut v3 = GpVec2d::zero();
        self.values(u, 0, &mut p, &mut v1, &mut v2, &mut v3);
        p
    }

    /// `EvalD1(U)` (`cxx:555-562`).
    pub fn eval_d1(&self, u: f64) -> (GpPnt2d, GpVec2d) {
        let mut p = self.point;
        let mut d1 = GpVec2d::zero();
        let mut v2 = GpVec2d::zero();
        let mut v3 = GpVec2d::zero();
        self.values(u, 1, &mut p, &mut d1, &mut v2, &mut v3);
        (p, d1)
    }

    /// `EvalD2(U)` (`cxx:566-577`).
    pub fn eval_d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        let mut p = self.point;
        let mut d1 = GpVec2d::zero();
        let mut d2 = GpVec2d::zero();
        let mut v3 = GpVec2d::zero();
        self.values(u, 2, &mut p, &mut d1, &mut d2, &mut v3);
        (p, d1, d2)
    }

    /// `EvalD3(U)` (`cxx:581-589`).
    pub fn eval_d3(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d, GpVec2d) {
        let mut p = self.point;
        let mut d1 = GpVec2d::zero();
        let mut d2 = GpVec2d::zero();
        let mut d3 = GpVec2d::zero();
        self.values(u, 3, &mut p, &mut d1, &mut d2, &mut d3);
        (p, d1, d2, d3)
    }

    /// `EvalDN(U, N)` (`cxx:594-618`), for `N = 1, 2, 3`.
    pub fn eval_dn(&self, u: f64, n: i32) -> GpVec2d {
        let mut p = self.point;
        let mut v1 = GpVec2d::zero();
        let mut v2 = GpVec2d::zero();
        let mut v3 = GpVec2d::zero();
        self.values(u, n, &mut p, &mut v1, &mut v2, &mut v3);
        match n {
            1 => v1,
            2 => v2,
            3 => v3,
            // OCCT throws Geom2d_UndefinedDerivative; the trait returns a zero
            // vector for that case (see `Curve2d::eval_dn`).
            _ => GpVec2d::zero(),
        }
    }

    /// `SearchBound(U1, U2)` (`cxx:622-647`).
    pub fn search_bound(&self, u1: f64, u2: f64) -> f64 {
        const TOL: f64 = PCONFUSION;
        let dist_max2 = self.dist_max * self.dist_max;
        let mut u11 = u1;
        let mut u22 = u2;
        let mut dist1 = self.distance(u11);
        let mut u_mid = 0.0;

        while (u22 - u11) > TOL {
            u_mid = 0.5 * (u22 + u11);
            let dist_mid = self.distance(u_mid);
            if (dist1 > dist_max2) == (dist_mid > dist_max2) {
                u11 = u_mid;
                dist1 = dist_mid;
            } else {
                u22 = u_mid;
            }
        }
        u_mid
    }

    /// `ComputeIntervals()` (`cxx:655-750`).
    pub fn compute_intervals(&mut self) {
        let mut u_proj = 0.0;
        self.is_empty = false;
        self.shift_parameter = 0.0;
        let mut ya_proj = false;
        let dist_max2 = self.dist_max * self.dist_max;

        let u1 = self.curve.first_parameter();
        let u2 = self.curve.last_parameter();
        let dist1 = self.distance(u1);
        let dist2 = self.distance(u2);
        let mut dist_proj = INFINITE;

        let proj = ProjectPointOnCurve2d::with_curve_range(&self.point, &self.curve, u1, u2);
        if proj.nb_points() > 0 {
            u_proj = proj.lower_distance_parameter().unwrap_or(0.0);
            dist_proj = self.distance(u_proj);
            ya_proj = true;
        }

        let u_start;
        let u_end;

        if dist1 < dist_max2 && dist2 < dist_max2 {
            if dist_proj > dist_max2 && ya_proj {
                self.is_empty = true;
            } else {
                self.start_intervals.push(u1);
                self.end_intervals.push(u2);
            }
            return;
        } else if dist1 > dist_max2 && dist2 > dist_max2 {
            if dist_proj < dist_max2 {
                u_start = self.search_bound(u1, u_proj);
                u_end = self.search_bound(u_proj, u2);
            } else {
                self.is_empty = true;
                return;
            }
        } else if dist1 < dist_max2 {
            u_start = u1;
            u_end = self.search_bound(u1, u2);
        } else {
            u_end = u2;
            u_start = self.search_bound(u1, u2);
        }
        self.start_intervals.push(u_start);
        self.end_intervals.push(u_end);

        if self.curve.is_periodic() {
            if self.start_intervals.len() > 1 {
                if self.end_last() == self.curve.last_parameter()
                    && self.start_first() == self.curve.first_parameter()
                {
                    self.start_intervals.remove(0);
                    self.end_intervals.pop();

                    self.shift_parameter = self.curve.period() - self.start_first();
                    let n = self.start_intervals.len();
                    for k in 0..n {
                        self.end_intervals[k] += self.shift_parameter;
                        self.start_intervals[k] += self.shift_parameter;
                    }
                    self.start_intervals[0] = 0.0;
                }
            }
        }
    }

    /// `Parameter(P)` (`cxx:776-815`).
    pub fn parameter(&self, p: &GpPnt2d) -> f64 {
        const TOL: f64 = CONFUSION;

        if p.is_equal(&self.point_start_bis, TOL) {
            return self.start_value(self.bis_interval);
        }
        if p.is_equal(&self.point_end_bis, TOL) {
            return self.end_value(self.bis_interval);
        }

        if self.extension_start {
            let axe = GpAx2d::new(
                self.point_start_bis,
                dir2d_or_zero(
                    self.point_start_bis.x() - p.x(),
                    self.point_start_bis.y() - p.y(),
                ),
            );
            let u = line_parameter_ax2d(&axe, p);
            let proj = line_value_ax2d(u, &axe);
            if proj.is_equal(p, TOL) && u < 0.0 {
                return u + self.start_value(self.bis_interval);
            }
        }
        if self.extension_end {
            let axe = GpAx2d::new(
                self.point_end_bis,
                dir2d_or_zero(
                    p.x() - self.point_end_bis.x(),
                    p.y() - self.point_end_bis.y(),
                ),
            );
            let u = line_parameter_ax2d(&axe, p);
            let proj = line_value_ax2d(u, &axe);
            if proj.is_equal(p, TOL) && u > 0.0 {
                return u + self.end_value(self.bis_interval);
            }
        }
        let mut u_on_curve = 0.0;
        let proj = ProjectPointOnCurve2d::with_curve_range(
            p,
            &self.curve,
            self.curve.first_parameter(),
            self.curve.last_parameter(),
        );
        if proj.nb_points() > 0 {
            u_on_curve = proj.lower_distance_parameter().unwrap_or(0.0);
        }
        self.link_curve_bis(u_on_curve)
    }
}

impl Curve2d for BisectorBisecPC {
    fn d0(&self, u: f64) -> GpPnt2d {
        self.eval_d0(u)
    }

    fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d) {
        self.eval_d1(u)
    }

    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        self.eval_d2(u)
    }

    fn d3(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d, GpVec2d) {
        self.eval_d3(u)
    }

    fn eval_dn(&self, u: f64, n: i32) -> GpVec2d {
        BisectorBisecPC::eval_dn(self, u, n)
    }

    fn first_parameter(&self) -> f64 {
        BisectorBisecPC::first_parameter(self)
    }

    fn last_parameter(&self) -> f64 {
        BisectorBisecPC::last_parameter(self)
    }

    /// `Bisector_BisecPC::IsPeriodic` (`cxx:304-307`): always false.
    fn is_periodic(&self) -> bool {
        false
    }

    fn continuity(&self) -> u8 {
        BisectorBisecPC::continuity(self)
    }

    fn transform(&mut self, t: &GpTrsf2d) {
        self.transform_curve(t);
    }

    /// `Bisector_BisecPC::Reverse` (`cxx:152-155`) raises
    /// `Standard_NotImplemented`; the port panics for the same reason.
    fn reverse(&mut self) {
        panic!("Standard_NotImplemented: Bisector_BisecPC::Reverse");
    }

    fn reversed_parameter(&self, u: f64) -> f64 {
        BisectorBisecPC::reversed_parameter(self, u)
    }

    fn clone_dyn(&self) -> Box<dyn Curve2d> {
        Box::new(self.copy_bisector())
    }

    /// `STANDARD_TYPE(Bisector_BisecPC)` downcast (`Bisector_Inter.cxx:279`).
    fn bisector_bisec_pc(&self) -> Option<&BisectorBisecPC> {
        Some(self)
    }
}

impl BisectorCurve for BisectorBisecPC {
    fn is_extend_at_start(&self) -> bool {
        BisectorBisecPC::is_extend_at_start(self)
    }

    fn is_extend_at_end(&self) -> bool {
        BisectorBisecPC::is_extend_at_end(self)
    }

    fn nb_intervals(&self) -> i32 {
        BisectorBisecPC::nb_intervals(self)
    }

    fn interval_first(&self, index: i32) -> Result<f64, BisectorBisecAnaError> {
        BisectorBisecPC::interval_first(self, index)
    }

    fn interval_last(&self, index: i32) -> Result<f64, BisectorBisecAnaError> {
        BisectorBisecPC::interval_last(self, index)
    }
}
