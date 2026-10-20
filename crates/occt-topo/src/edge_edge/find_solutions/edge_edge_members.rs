use super::super::prelude::*;
use super::super::*;

use occt_core::bnd::BndBox;
use occt_core::precision::{INFINITE, SQUARE_CONFUSION};

use super::boxes::{bnd_build_box, split_range_on_segments};
use super::curve_kind::GeomCurveKind;
use super::dist_pc::{dist_pc, find_dist_pc, DistState, ProjPc};
use super::parameters::find_parameters;
use super::resolution::{is_closed, resolution_of};

// ---------------------------------------------------------------------------
// Member helpers
// ---------------------------------------------------------------------------

impl EdgeEdge {
    /// `IntTools_EdgeEdge::FindSolutions(theRanges1, theRanges2, bSplit2)`
    /// (`IntTools_EdgeEdge.cxx:290-349`) followed by `MergeSolutions`
    /// (`Perform`, `:241-242`).
    ///
    /// Replaces the previous substitute ("the extrema of the distance within
    /// tolerance are the intersections"): OCCT recurses on the parameter boxes
    /// so that every crossing whose boxes stay in contact is reported, even
    /// when it is not a distance extremum of the whole range.
    pub(in crate::edge_edge) fn find_solutions(&mut self) {
        let (t11, t12) = (self.r1().first, self.r1().last);
        let (t21, t22) = (self.r2().first, self.r2().last);
        let c1 = self.curve1.clone().unwrap();
        let c2 = self.curve2.clone().unwrap();

        let mut b_is_closed2 = is_closed(&*c2, t21, t22, self.tol2, self.res2);
        if b_is_closed2 {
            let b1 = bnd_build_box(&*c1, t11, t12, self.tol1);
            let p = c2.d0(t21);
            b_is_closed2 = !b1.is_out(&p);
        }

        let mut ranges1: Vec<IntRange> = Vec::new();
        let mut ranges2: Vec<IntRange> = Vec::new();
        let b_split2;

        if !b_is_closed2 {
            let b1 = bnd_build_box(&*c1, t11, t12, self.tol1);
            let b2 = bnd_build_box(&*c2, t21, t22, self.tol2);
            self.find_solutions_range(
                IntRange::new_unchecked(t11, t12),
                &b1,
                IntRange::new_unchecked(t21, t22),
                &b2,
                &mut ranges1,
                &mut ranges2,
            );
            self.merge_solutions_boxes(&ranges1, &ranges2, false);
            return;
        }

        if self.check_coincidence(t11, t12, t21, t22, self.tol, self.res1) == 0 {
            // `cxx:320-325`: the whole range is one common part when
            // `CheckCoincidence` reports the patches as coincident.
            ranges1.push(self.r1());
            ranges2.push(self.r2());
            self.merge_solutions_boxes(&ranges1, &ranges2, false);
            return;
        }

        let nb1 = if is_closed(&*c1, t11, t12, self.tol1, self.res1) {
            2
        } else {
            1
        };
        let segments1 = split_range_on_segments(t11, t12, self.res1, nb1);
        let segments2 = split_range_on_segments(t21, t22, self.res2, 2);

        for r1_val in &segments1 {
            let b1 = bnd_build_box(&*c1, r1_val.first, r1_val.last, self.tol1);
            for r2_val in &segments2 {
                let b2 = bnd_build_box(&*c2, r2_val.first, r2_val.last, self.tol2);
                self.find_solutions_range(*r1_val, &b1, *r2_val, &b2, &mut ranges1, &mut ranges2);
            }
        }
        b_split2 = segments2.len() > 1;
        self.merge_solutions_boxes(&ranges1, &ranges2, b_split2);
    }

    /// `IntTools_EdgeEdge::FindSolutions(theR1, theBox1, theR2, theBox2, …)`
    /// (`IntTools_EdgeEdge.cxx:353-549`) — the recursion.
    #[allow(clippy::too_many_arguments)]
    fn find_solutions_range(
        &self,
        r1: IntRange,
        box1: &BndBox,
        r2: IntRange,
        box2: &BndBox,
        ranges1: &mut Vec<IntRange>,
        ranges2: &mut Vec<IntRange>,
    ) {
        let c1 = self.curve1.as_ref().unwrap();
        let c2 = self.curve2.as_ref().unwrap();
        let (mut t11_val, mut t12_val) = (r1.first, r1.last);
        let (mut t21_val, mut t22_val) = (r2.first, r2.last);
        let mut b1 = *box1;
        // `cxx:371`: `aB2` starts at the incoming box and is rebuilt at step 2;
        // the rebuilt value carries into the next iteration's step 1.
        let mut b2 = *box2;

        let mut b_out;
        let mut b_thin = false;
        let mut b_stop = false;
        let mut i_com = 1;
        loop {
            let (tb11, tb12) = (t11_val, t12_val);
            let (tb21, tb22) = (t21_val, t22_val);

            // 1. Find the parameters of the second edge inside the box of the
            //    first one.
            b_out = b1.is_out_box(&b2);
            if b_out {
                break;
            }
            b_thin = (t12_val - t11_val) < self.res1
                || (b1.is_x_thin(self.tol)
                    && b1.is_y_thin(self.tol)
                    && b1.is_z_thin(self.tol));
            match find_parameters(
                &**c2,
                self.kind2,
                tb21,
                tb22,
                self.tol2,
                self.res2,
                self.ptol2,
                self.res_coeff2,
                &b1,
                t21_val,
                t22_val,
            ) {
                None => {
                    b_out = true;
                }
                Some((t21, t22)) => {
                    t21_val = t21;
                    t22_val = t22;
                    b_out = false;
                }
            }
            if b_out || b_thin {
                break;
            }

            // 2. Build the box of the second edge and find the parameters of
            //    the first one inside it.
            b2 = bnd_build_box(&**c2, t21_val, t22_val, self.tol2);
            b_out = b1.is_out_box(&b2);
            if b_out {
                break;
            }
            b_thin = (t22_val - t21_val) < self.res2
                || (b2.is_x_thin(self.tol)
                    && b2.is_y_thin(self.tol)
                    && b2.is_z_thin(self.tol));
            match find_parameters(
                &**c1,
                self.kind1,
                tb11,
                tb12,
                self.tol1,
                self.res1,
                self.ptol1,
                self.res_coeff1,
                &b2,
                t11_val,
                t12_val,
            ) {
                None => {
                    b_out = true;
                }
                Some((t11, t12)) => {
                    t11_val = t11;
                    t12_val = t12;
                    b_out = false;
                }
            }
            if b_out || b_thin {
                break;
            }

            // 3. Check whether it makes sense to continue.
            let mut small_step1 = (tb12 - tb11) / 250.0;
            let mut small_step2 = (tb22 - tb21) / 250.0;
            if small_step1 < self.res1 {
                small_step1 = self.res1;
            }
            if small_step2 < self.res2 {
                small_step2 = self.res2;
            }
            if ((t11_val - tb11) < small_step1)
                && ((tb12 - t12_val) < small_step1)
                && ((t21_val - tb21) < small_step2)
                && ((tb22 - t22_val) < small_step2)
            {
                b_stop = true;
            } else {
                b1 = bnd_build_box(&**c1, t11_val, t12_val, self.tol1);
            }
            if b_stop {
                break;
            }
        }

        if b_out {
            // no intersection
            return;
        }
        if !b_thin {
            // check the curves for coincidence on the ranges
            i_com = self.check_coincidence(t11_val, t12_val, t21_val, t22_val, self.tol, self.res1);
            if i_com == 0 {
                b_thin = true;
            }
        }
        if b_thin {
            if i_com != 0 {
                // check the intermediate points
                let t1 = (t11_val + t12_val) * 0.5;
                let p1 = c1.d0(t1);
                let proj = ProjPc::new(&**c2, t21_val, t22_val);
                let b_sol = match proj.project(&p1) {
                    Some((d, _)) => d <= self.tol,
                    None => {
                        let t2 = (t21_val + t22_val) * 0.5;
                        let p2 = c2.d0(t2);
                        p1.distance(&p2) <= self.tol
                    }
                };
                if !b_sol {
                    return;
                }
            }
            // add the common part
            ranges1.push(IntRange::new_unchecked(t11_val, t12_val));
            ranges2.push(IntRange::new_unchecked(t21_val, t22_val));
            return;
        }

        if !self.is_intersection(t11_val, t12_val, t21_val, t22_val) {
            return;
        }

        // Split the ranges on segments and repeat.
        let b1_first = bnd_build_box(&**c1, t11_val, t12_val, self.tol1);
        let b1_sq_extent = b1_first.square_extent();
        let r2_val = IntRange::new_unchecked(t21_val, t22_val);
        let b2_seg = bnd_build_box(&**c2, t21_val, t22_val, self.tol2);
        let segments1 = split_range_on_segments(t11_val, t12_val, self.res1, 3);
        let nb1 = segments1.len();
        for r1_val in segments1 {
            let b1s = bnd_build_box(&**c1, r1_val.first, r1_val.last, self.tol1);
            if !b1s.is_out_box(&b2_seg)
                && (nb1 == 1 || b1s.square_extent() < b1_sq_extent)
            {
                self.find_solutions_range(r1_val, &b1s, r2_val, &b2_seg, ranges1, ranges2);
            }
        }
    }

    /// `IntTools_EdgeEdge::MergeSolutions` (`IntTools_EdgeEdge.cxx:675-779`):
    /// groups the parameter ranges found by the recursion and stores them.
    fn merge_solutions_boxes(
        &mut self,
        ranges1: &[IntRange],
        ranges2: &[IntRange],
        b_split2: bool,
    ) {
        let nb_cp = ranges1.len();
        if nb_cp == 0 {
            return;
        }
        let c1 = self.curve1.clone().unwrap();
        let c2 = self.curve2.clone().unwrap();
        let res1_val = resolution_of(&*c1, self.kind1, self.res_coeff1, self.tol);
        let res2_val = resolution_of(&*c2, self.kind2, self.res_coeff2, self.tol);
        let (t11, t12) = (self.r1().first, self.r1().last);
        let (t21, t22) = (self.r2().first, self.r2().last);
        let d_tr1 = 20.0 * res1_val;
        let d_tr2 = 20.0 * res2_val;
        let mut type_edge = false;
        let mut used = vec![false; nb_cp];

        let mut i = 0usize;
        while i < nb_cp {
            if used[i] {
                i += 1;
                continue;
            }
            let (mut ti11, mut ti12) = (ranges1[i].first, ranges1[i].last);
            let (mut ti21, mut ti22) = (ranges2[i].first, ranges2[i].last);
            used[i] = true;

            let mut j = i + 1;
            while j < nb_cp {
                if used[j] {
                    j += 1;
                    continue;
                }
                let (tj11, tj12) = (ranges1[j].first, ranges1[j].last);
                let (tj21, tj22) = (ranges2[j].first, ranges2[j].last);
                let mut b_cond = ((ti12 - tj11).abs() < d_tr1)
                    || (tj11 > ti11 && tj11 < ti12)
                    || (ti11 > tj11 && ti11 < tj12)
                    || (b_split2 && (tj12 - ti11).abs() < d_tr1);
                if b_cond && b_split2 {
                    b_cond = ((ti22.max(tj22) - ti21.min(tj21))
                        - ((ti22 - ti21) + (tj22 - tj21)))
                        .abs()
                        < d_tr2
                        || (tj21 > ti21 && tj21 < ti22)
                        || (ti21 > tj21 && ti21 < tj22);
                }
                if b_cond {
                    ti11 = ti11.min(tj11);
                    ti12 = ti12.max(tj12);
                    ti21 = ti21.min(tj21);
                    ti22 = ti22.max(tj22);
                    used[j] = true;
                    j += 1;
                } else if !b_split2 {
                    i = j;
                    break;
                } else {
                    j += 1;
                }
            }

            if ((t11 - ti11).abs() < self.res1 && (t12 - ti12).abs() < self.res1)
                || ((t21 - ti21).abs() < self.res2 && (t22 - ti22).abs() < self.res2)
            {
                type_edge = true;
                self.common_parts.clear();
            }

            self.add_solution_occt(ti11, ti12, ti21, ti22, type_edge);
            if type_edge {
                break;
            }
            if b_split2 {
                i += 1;
            }
        }
    }

    /// `IntTools_EdgeEdge::AddSolution` (`IntTools_EdgeEdge.cxx:780-825`).
    ///
    /// `mySwap` is honoured exactly as OCCT does: the stored range and the
    /// vertex parameters are expressed on the **caller's** first/second edge.
    fn add_solution_occt(&mut self, t11: f64, t12: f64, t21: f64, t22: f64, is_edge: bool) {
        if is_edge {
            let range = if !self.swapped {
                IntRange::new_unchecked(t11, t12)
            } else {
                IntRange::new_unchecked(t21, t22)
            };
            self.common_parts
                .push(CommonPrt::with(CommonPartType::Edge, range, None, Vec::new()));
            return;
        }
        let (t1, t2) = self.find_best_solution_ranges(t11, t12, t21, t22);
        // `cxx:810-819`: when swapped, `aT1` lives on the caller's *second*
        // edge and the vertex parameters are exchanged.
        let (e1, e2) = if !self.swapped { (t1, t2) } else { (t2, t1) };
        let Some(c1) = self.curve1.as_ref() else {
            return;
        };
        let p = c1.d0(t1);
        self.points.push(PntOn2Faces::new(0, 1, p, p, (e1, 0.0), (e2, 0.0)));
    }

    /// `IntTools_EdgeEdge::FindBestSolution` (`IntTools_EdgeEdge.cxx:826-901`):
    /// the parameter pair of the closest approach inside the solution range.
    fn find_best_solution_ranges(
        &self,
        t11: f64,
        t12: f64,
        t21: f64,
        t22: f64,
    ) -> (f64, f64) {
        const A_SOL_CRITERIA: f64 = 5.0e-16;
        const A_TOUCH_CRITERIA: f64 = 5.0e-13;
        let c1 = self.curve1.as_ref().unwrap();
        let c2 = self.curve2.as_ref().unwrap();
        let mut d_min = INFINITE;
        let mut t1 = 0.0;
        let mut t2_val = 0.0;
        let mut b_touch = false;
        let mut b_touch_confirm = false;
        let (mut t11_touch, mut t12_touch) = (t11, t12);
        let (mut t21_touch, mut t22_touch) = (t21, t22);
        let mut is_sol_found = false;

        let res1 = resolution_of(&**c1, self.kind1, self.res_coeff1, self.tol);
        let nb_s = 10;
        let ranges = split_range_on_segments(t11, t12, 3.0 * res1, nb_s);
        let proj = ProjPc::new(&**c2, t21, t22);

        for r1 in ranges {
            let mut st = DistState { d_max: self.tol, t1_max: 0.0, t2_max: 0.0 };
            let i_err = find_dist_pc(
                &proj,
                &**c1,
                r1.first,
                r1.last,
                A_SOL_CRITERIA,
                self.ptol1,
                &mut st,
                false,
            );
            let d = st.d_max;
            if i_err != 1 {
                if d < d_min {
                    t1 = st.t1_max;
                    t2_val = st.t2_max;
                    d_min = d;
                    is_sol_found = true;
                }
                if d < A_TOUCH_CRITERIA {
                    if b_touch {
                        t12_touch = st.t1_max;
                        t22_touch = st.t2_max;
                        b_touch_confirm = true;
                    } else {
                        t11_touch = st.t1_max;
                        t21_touch = st.t2_max;
                        b_touch = true;
                    }
                }
            }
        }

        if !is_sol_found || b_touch_confirm {
            t1 = (t11_touch + t12_touch) * 0.5;
            let (i_err, _d, t2) = dist_pc(&proj, &**c1, t1, A_SOL_CRITERIA, -1);
            t2_val = t2;
            if i_err == 1 {
                t2_val = (t21_touch + t22_touch) * 0.5;
            }
        }
        (t1, t2_val)
    }

    /// `IntTools_EdgeEdge::CheckCoincidence` (`IntTools_EdgeEdge.cxx:1150-1206`).
    ///
    /// `0` = the patches are coincident, `1` = a point of C1 cannot be
    /// projected, `2` = the distance is too big.
    fn check_coincidence(
        &self,
        t11: f64,
        t12: f64,
        t21: f64,
        t22: f64,
        criteria: f64,
        curve_res1: f64,
    ) -> i32 {
        let c1 = self.curve1.as_ref().unwrap();
        let c2 = self.curve2.as_ref().unwrap();
        let proj = ProjPc::new(&**c2, t21, t22);
        let nb = 10;
        let ranges = split_range_on_segments(t11, t12, curve_res1, nb);
        let nb1 = ranges.len() as i32;

        // 1. Express evaluation. The 7-argument `DistPC` writes the distance
        //    straight into `aDmax` (`cxx:1174`), so there is no extremum
        //    tracking here.
        let mut d_max_val = -1.0;
        let mut i_err = 0;
        for k in 0..nb1.saturating_sub(1) {
            let t1b = ranges[k as usize].last;
            let (e, d, _t2) = dist_pc(&proj, &**c1, t1b, criteria, 1);
            d_max_val = d;
            i_err = e;
            if i_err != 0 {
                return i_err;
            }
        }
        // If the ranges are shorter than `theCurveRes1`, there is no need for
        // the deep evaluation.
        if nb1 < nb {
            return i_err;
        }

        // 2. Deep evaluation.
        for k in 1..nb1.saturating_sub(1) {
            let (t1a, t1b) = (ranges[k as usize].first, ranges[k as usize].last);
            let mut st = DistState { d_max: d_max_val, t1_max: 0.0, t2_max: 0.0 };
            i_err = find_dist_pc(
                &proj,
                &**c1,
                t1a,
                t1b,
                criteria,
                curve_res1,
                &mut st,
                true,
            );
            d_max_val = st.d_max;
            if i_err != 0 {
                return i_err;
            }
        }
        i_err
    }

    /// `IntTools_EdgeEdge::IsIntersection` (`IntTools_EdgeEdge.cxx:1060-1146`).
    fn is_intersection(&self, t11: f64, t12: f64, t21: f64, t22: f64) -> bool {
        let c1 = self.curve1.as_ref().unwrap();
        let c2 = self.curve2.as_ref().unwrap();
        let coef;
        if ((t12 - t11) > 1.0e5 * self.res1) && ((t22 - t21) > 1.0e5 * self.res2) {
            coef = 5000.0;
        } else {
            let tr_min = ((t12 - t11) / self.res1).min((t22 - t21) / self.res2);
            let coef_tmp = tr_min / 100.0;
            coef = if coef_tmp < 1.0 { 1.0 } else { coef_tmp };
        }
        let mut criteria = coef * self.tol;
        criteria *= criteria;

        let (p11, v11) = c1.d1(t11);
        let (p12, v12) = c1.d1(t12);
        let (p21, v21) = c2.d1(t21);
        let (p22, v22) = c2.d1(t22);

        let b_small_11_21 = p11.square_distance(&p21) < criteria;
        let b_small_11_22 = p11.square_distance(&p22) < criteria;
        let b_small_12_21 = p12.square_distance(&p21) < criteria;
        let b_small_12_22 = p12.square_distance(&p22) < criteria;

        let mut b_ret = true;
        if (b_small_11_21 && b_small_12_22) || (b_small_11_22 && b_small_12_21) {
            if coef == 1.0 {
                return b_ret;
            }
            let an_angle_criteria = 5.0e-3;
            let mut an_angle1 = 0.0;
            let mut an_angle2 = 0.0;
            let v_ok = |v: &GpVec| v.dot(v) > SQUARE_CONFUSION;
            if v_ok(&v11) && v_ok(&v12) && v_ok(&v21) && v_ok(&v22) {
                if b_small_11_21 && b_small_12_22 {
                    an_angle1 = v11.angle(&v21);
                    an_angle2 = v12.angle(&v22);
                } else {
                    an_angle1 = v11.angle(&v22);
                    an_angle2 = v12.angle(&v21);
                }
            }
            if (an_angle1 < an_angle_criteria
                || (std::f64::consts::PI - an_angle1) < an_angle_criteria)
                || (an_angle2 < an_angle_criteria
                    || (std::f64::consts::PI - an_angle2) < an_angle_criteria)
            {
                let proj = ProjPc::new(&**c2, t21, t22);
                let mut st = DistState { d_max: INFINITE, t1_max: 0.0, t2_max: 0.0 };
                let i_err =
                    find_dist_pc(&proj, &**c1, t11, t12, self.tol, self.res1, &mut st, false);
                b_ret = i_err == 2;
            }
        }
        b_ret
    }
}
