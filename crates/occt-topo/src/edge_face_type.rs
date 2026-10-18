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
        let pf = curve.d0(af1);
        let pl = curve.d0(al1);
        let df1 = pf.distance(&pl);
        let cr = curve_resolution(curve.as_ref(), self.criteria);
        let is_whole_range =
            (af1 - self.range.first).abs() < cr && (al1 - self.range.last).abs() < cr;
        if df1 > self.criteria * 2.0 && is_whole_range {
            cp.part_type = CommonPartType::Edge;
            return 0;
        }
        if is_whole_range {
            let tm = 0.5 * (af1 + al1);
            if pf.distance(&curve.d0(tm)) > self.criteria * 2.0 {
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
        let tf = cp.range.first;
        let tl = cp.range.last;
        let curve = self.curve.clone().expect("curve set");
        let cr = curve_resolution(curve.as_ref(), self.criteria);
        if (tf - self.range.first).abs() < cr && (tl - self.range.last).abs() < cr {
            return false; // whole range: keep EDGE
        }

        let (min_d, max_d, min_t) = self.distance_profile(tf, tl, 32);
        // Extrema parallel case: the distance is nearly constant over the range.
        if max_d - min_d <= 0.05 * (min_d + self.criteria.max(1e-9)) {
            return false;
        }
        let mut dist2 = min_d * min_d;
        let mut a_tx = min_t;

        // Exact curve–surface intersection fallback (the `Extrema` aNbExt == 0
        // branch of `IntTools_EdgeFace::CheckTouch`): when the sampled profile
        // found no near-surface minimum, consult the exact intersector.
        if dist2 > self.criteria * self.criteria && tl > tf {
            let surface = self.surface.clone().expect("surface set");
            let (u0, u1, v0, v1) = surface_sample_bounds(surface.as_ref());
            if let Ok(hr) = perform_curve_surface(
                curve.as_ref(),
                surface.as_ref(),
                (tf, tl),
                (u0, u1, v0, v1),
            ) {
                for i in 0..hr.nb_points() {
                    let p = hr.point(i);
                    if p.param() >= tf && p.param() <= tl {
                        dist2 = 0.0;
                        a_tx = p.param();
                        break;
                    }
                }
            }
        }

        let b1 = self.distance_function(tf) + self.criteria;
        if b1 * b1 < dist2 {
            dist2 = b1 * b1;
            a_tx = tf;
        }
        let b2 = self.distance_function(tl) + self.criteria;
        if b2 * b2 < dist2 {
            dist2 = b2 * b2;
            a_tx = tl;
        }
        let bm = self.distance_function(0.5 * (tf + tl)) + self.criteria;
        if bm * bm < dist2 {
            dist2 = bm * bm;
            a_tx = 0.5 * (tf + tl);
        }

        if dist2 > self.criteria * self.criteria {
            return false;
        }
        *tx = a_tx;
        if (a_tx - tf).abs() < PCONFUSION {
            return true;
        }
        if (a_tx - tl).abs() < PCONFUSION {
            return true;
        }
        if a_tx > tf && a_tx < tl {
            return true;
        }
        false
    }

    /// Vertex-specific touch refinement. Port of `IntTools_EdgeFace::CheckTouchVertex`.
    fn check_touch_vertex(&self, cp: &CommonPrt, tx: &mut f64) -> bool {
        let tf = cp.range.first;
        let tl = cp.range.last;
        let curve = self.curve.clone().expect("curve set");
        let kind = curve_kind(curve.as_ref());
        let eps_t = if kind == CurveKind::Line { 9e-5 } else { 8e-5 };
        let tm = 0.5 * (tf + tl);
        let dist2 = {
            let d = self.distance_function(tm);
            d * d
        };
        if tl <= tf {
            return false;
        }
        let (min_d, max_d, min_t) = self.distance_profile(tf, tl, 32);
        if max_d - min_d <= 0.05 * (min_d + self.criteria.max(1e-9)) {
            return false;
        }
        let dist2_new = min_d * min_d;
        if dist2_new > dist2 {
            *tx = tm;
            return true;
        }
        if dist2_new > self.criteria * self.criteria {
            return false;
        }
        let a_tx = min_t;
        if (a_tx - tf).abs() < eps_t {
            return false;
        }
        if (a_tx - tl).abs() < eps_t {
            return false;
        }
        if a_tx > tf && a_tx < tl {
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
            let kind = self.common_parts[i].part_type;
            let cp = self.common_parts[i].clone();
            let mut tx = 0.0;
            if kind == CommonPartType::Edge {
                if self.check_touch(&cp, &mut tx) {
                    self.common_parts[i].part_type = CommonPartType::Vertex;
                    self.common_parts[i].vertex_parameter1 = Some(tx);
                }
            } else if kind == CommonPartType::Vertex
                && self.check_touch_vertex(&cp, &mut tx)
            {
                self.common_parts[i].vertex_parameter1 = Some(tx);
            }
        }
    }
}
