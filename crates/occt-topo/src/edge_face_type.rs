//! MakeType / CheckTouch / touch refinement for EdgeFace.
use occt_core::precision::PCONFUSION;

use crate::edge_face::EdgeFace;
use crate::edge_face_kind::{curve_kind, curve_resolution, surface_sample_bounds, CurveKind};
use crate::intcurvesurface::perform_curve_surface;
use crate::inttools_data::{CommonPartType, CommonPrt};

impl EdgeFace {
    /// Classify a range as `TopAbs_EDGE` (whole-range coincidence) or
    /// `TopAbs_VERTEX` (touch / piercing). Port of `IntTools_EdgeFace::MakeType`.
    ///
    /// VERTEX keeps the original range and stores the touch parameter in
    /// [`CommonPrt::vertex_parameter1`].
    pub(crate) fn make_type(&mut self, cp: &mut CommonPrt) -> i32 {
        if cp.all_null_flag {
            cp.part_type = CommonPartType::Edge;
            return 0;
        }
        let af1 = cp.range.first;
        let al1 = cp.range.last;
        let curve = self.curve.clone().expect("curve set");
        let a_pf = curve.d0(af1);
        let a_pl = curve.d0(al1);
        let df1 = a_pf.distance(&a_pl);
        let a_cr = curve_resolution(curve.as_ref(), self.criteria);
        let is_whole_range =
            (af1 - self.range.first).abs() < a_cr && (al1 - self.range.last).abs() < a_cr;

        if df1 > self.criteria * 2.0 && is_whole_range {
            cp.part_type = CommonPartType::Edge;
            return 0;
        }
        if is_whole_range {
            let tm = 0.5 * (af1 + al1);
            if a_pf.distance(&curve.d0(tm)) > self.criteria * 2.0 {
                cp.part_type = CommonPartType::Edge;
                return 0;
            }
        }
        let mut tm = 0.5 * (af1 + al1);
        if !self.check_touch(cp, &mut tm) {
            tm = 0.5 * (af1 + al1);
        }
        cp.part_type = CommonPartType::Vertex;
        cp.vertex_parameter1 = Some(tm);
        0
    }

    /// Whether the range contains a touch point within `myCriteria`, and the
    /// touch parameter. Port of `IntTools_EdgeFace::CheckTouch`.
    fn check_touch(&self, cp: &CommonPrt, tx: &mut f64) -> bool {
        let a_tf = cp.range.first;
        let a_tl = cp.range.last;
        let curve = self.curve.clone().expect("curve set");
        let a_cr = curve_resolution(curve.as_ref(), self.criteria);
        if (a_tf - self.range.first).abs() < a_cr && (a_tl - self.range.last).abs() < a_cr {
            return false; // whole range: keep EDGE
        }

        let (min_d, max_d, min_t) = self.distance_profile(a_tf, a_tl, 32);
        // Extrema parallel case: the distance is nearly constant over the range.
        if max_d - min_d <= 0.05 * (min_d + self.criteria.max(1e-9)) {
            return false;
        }
        let mut a_dist2 = min_d * min_d;
        let mut a_tx = min_t;

        // Exact curve–surface intersection fallback (the `Extrema` aNbExt == 0
        // branch of `IntTools_EdgeFace::CheckTouch`): when the sampled profile
        // found no near-surface minimum, consult the exact intersector.
        if a_dist2 > self.criteria * self.criteria && a_tl > a_tf {
            let surface = self.surface.clone().expect("surface set");
            let (u0, u1, v0, v1) = surface_sample_bounds(surface.as_ref());
            if let Ok(hr) = perform_curve_surface(
                curve.as_ref(),
                surface.as_ref(),
                (a_tf, a_tl),
                (u0, u1, v0, v1),
            ) {
                for i in 0..hr.nb_points() {
                    let p = hr.point(i);
                    if p.param() >= a_tf && p.param() <= a_tl {
                        a_dist2 = 0.0;
                        a_tx = p.param();
                        break;
                    }
                }
            }
        }

        let b1 = self.distance_function(a_tf) + self.criteria;
        if b1 * b1 < a_dist2 {
            a_dist2 = b1 * b1;
            a_tx = a_tf;
        }
        let b2 = self.distance_function(a_tl) + self.criteria;
        if b2 * b2 < a_dist2 {
            a_dist2 = b2 * b2;
            a_tx = a_tl;
        }
        let bm = self.distance_function(0.5 * (a_tf + a_tl)) + self.criteria;
        if bm * bm < a_dist2 {
            a_dist2 = bm * bm;
            a_tx = 0.5 * (a_tf + a_tl);
        }

        if a_dist2 > self.criteria * self.criteria {
            return false;
        }
        *tx = a_tx;
        if (a_tx - a_tf).abs() < PCONFUSION {
            return true;
        }
        if (a_tx - a_tl).abs() < PCONFUSION {
            return true;
        }
        if a_tx > a_tf && a_tx < a_tl {
            return true;
        }
        false
    }

    /// Vertex-specific touch refinement. Port of `IntTools_EdgeFace::CheckTouchVertex`.
    fn check_touch_vertex(&self, cp: &CommonPrt, tx: &mut f64) -> bool {
        let a_tf = cp.range.first;
        let a_tl = cp.range.last;
        let curve = self.curve.clone().expect("curve set");
        let a_type = curve_kind(curve.as_ref());
        let a_eps_t = if a_type == CurveKind::Line { 9e-5 } else { 8e-5 };
        let a_tm = 0.5 * (a_tf + a_tl);
        let a_dist2 = {
            let d = self.distance_function(a_tm);
            d * d
        };
        if a_tl <= a_tf {
            return false;
        }
        let (min_d, max_d, min_t) = self.distance_profile(a_tf, a_tl, 32);
        if max_d - min_d <= 0.05 * (min_d + self.criteria.max(1e-9)) {
            return false;
        }
        let a_dist2_new = min_d * min_d;
        if a_dist2_new > a_dist2 {
            *tx = a_tm;
            return true;
        }
        if a_dist2_new > self.criteria * self.criteria {
            return false;
        }
        let a_tx = min_t;
        if (a_tx - a_tf).abs() < a_eps_t {
            return false;
        }
        if (a_tx - a_tl).abs() < a_eps_t {
            return false;
        }
        if a_tx > a_tf && a_tx < a_tl {
            *tx = a_tx;
            return true;
        }
        false
    }

    /// The line/cylinder and circle/plane special treatment: refine EDGE/VERTEX
    /// common parts into touch points when the range is only tangent.
    /// Port of the Line/Cylinder and Circle/Plane blocks in
    /// `IntTools_EdgeFace::Perform`.
    pub(crate) fn refine_touch_parts(&mut self) {
        for i in 0..self.common_parts.len() {
            let a_type = self.common_parts[i].part_type;
            let cp = self.common_parts[i].clone();
            let mut tx = 0.0;
            if a_type == CommonPartType::Edge {
                if self.check_touch(&cp, &mut tx) {
                    self.common_parts[i].part_type = CommonPartType::Vertex;
                    self.common_parts[i].vertex_parameter1 = Some(tx);
                }
            } else if a_type == CommonPartType::Vertex
                && self.check_touch_vertex(&cp, &mut tx)
            {
                self.common_parts[i].vertex_parameter1 = Some(tx);
            }
        }
    }
}
