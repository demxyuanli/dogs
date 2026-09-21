use super::prelude::*;
use super::*;

// ---------------------------------------------------------------------------
// Phase 18c — full analytic edge/edge solvers (appended; the original `perform`
// dispatch above is untouched). These mirror the analytic branches of
// `IntTools_EdgeEdge::FindSolutions` at the level this port needs.
// ---------------------------------------------------------------------------

impl EdgeEdge {
    /// Select the best (non-degenerate) intersection candidate among
    /// `candidates`.
    ///
    /// Port of OCCT `IntTools_EdgeEdge::FindBestSolution`. Each candidate is a
    /// `(parameter on edge 1, parameter on edge 2, 3D point)`. A candidate is
    /// *degenerate* when the two curves' evaluated points at its parameters are
    /// farther apart than `tol` — such candidates are dropped. Among the rest,
    /// the candidate whose two curve points agree most closely (minimum
    /// residual distance) is returned. `None` when no non-degenerate candidate
    /// exists (or the edges are not prepared).
    pub fn find_best_solution(
        &self,
        candidates: &[(f64, f64, GpPnt)],
        tol: f64,
    ) -> Option<(f64, f64, GpPnt)> {
        let c1 = self.curve1.as_ref()?;
        let c2 = self.curve2.as_ref()?;
        let etol = tol.max(1e-9);
        let mut best: Option<(f64, usize)> = None;
        for (i, (t1, t2, _)) in candidates.iter().enumerate() {
            let residual = c1.d0(*t1).distance(&c2.d0(*t2));
            if residual > etol {
                continue; // degenerate — not an actual intersection
            }
            match best {
                Some((br, _)) if br <= residual => {}
                _ => best = Some((residual, i)),
            }
        }
        best.map(|(_, i)| candidates[i])
    }

    /// Full analytic line/line intersection.
    ///
    /// Mirrors the line/line branch of `IntTools_EdgeEdge::FindSolutions`:
    /// - parallel & coincident lines → the overlapping segment becomes a
    ///   [`CommonPrt`] (via [`is_coincident`](Self::is_coincident));
    /// - parallel & distinct lines → no intersection;
    /// - skew (non-parallel) lines → the closest-approach point, kept only when
    ///   the lines actually meet within tolerance and the point lies in both
    ///   edge ranges.
    ///
    /// The result vectors are cleared first (this is a self-contained parse).
    /// Returns `Ok(true)` when any intersection/coincidence was stored,
    /// `Ok(false)` otherwise. `Err` when an edge curve is not a line.
    pub fn compute_line_line_full(&mut self) -> Result<bool, String> {
        let c1 = self.curve1.clone().ok_or("EdgeEdge: edge1 has no curve")?;
        let c2 = self.curve2.clone().ok_or("EdgeEdge: edge2 has no curve")?;
        let l1 = line_of_curve(&*c1).ok_or("EdgeEdge: edge1 is not a line")?;
        let l2 = line_of_curve(&*c2).ok_or("EdgeEdge: edge2 is not a line")?;
        let tol = self.tol.max(1e-9);

        self.common_parts.clear();
        self.points.clear();

        let d1 = GpVec::from_xyz(l1.direction().xyz());
        let d2 = GpVec::from_xyz(l2.direction().xyz());
        if d1.cross_magnitude(&d2) <= tol {
            // Parallel: coincident → segment, distinct → empty.
            if self.is_coincident() {
                self.push_coincident_common_part();
                return Ok(true);
            }
            return Ok(false);
        }

        let Some((t, s, p)) = line_line_closest(&l1, &l2) else { return Ok(false) };
        let et1 = self.curve_to_edge_param(1, t);
        let et2 = self.curve_to_edge_param(2, s);
        if self.r1().contains(et1) && self.r2().contains(et2) {
            if c1.d0(t).distance(&c2.d0(s)) <= tol {
                self.points.push(PntOn2Faces::new(0, 1, p, p, (et1, 0.0), (et2, 0.0)));
                return Ok(true);
            }
        }
        Ok(false)
    }

    /// Full analytic circle/circle intersection.
    ///
    /// **OCCT has no such branch**: `IntTools_EdgeEdge::Perform`
    /// (`IntTools_EdgeEdge.cxx:185-243`) special-cases only line/line
    /// (`ComputeLineLine`, `:902-1058`) and sends every other pair — circles
    /// included — through `FindSolutions`/`MergeSolutions`. This helper and the
    /// `perform` dispatch that calls it are therefore a port-local fast path;
    /// the faithful path now exists ([`EdgeEdge::find_solutions`]), so removing
    /// the dispatch is a behaviour change that needs a gate run (board task
    /// R2-17). Kept for its direct tests until then.
    ///
    /// Mirrors the coplanar circle/circle branch of `IntTools_EdgeEdge`:
    /// - coincident circles → the overlapping arc becomes a [`CommonPrt`];
    /// - externally/internally separated circles → no intersection;
    /// - externally/internally tangent circles → one vertex hit;
    /// - intersecting circles → two vertex hits (radical-line construction).
    ///
    /// Non-coplanar circles fall back to the general solver
    /// ([`intersect_edges`](Self::intersect_edges)). The result vectors are
    /// cleared first (self-contained parse). Returns `Ok(true)` when any result
    /// was stored. `Err` when an edge curve is not a circle.
    pub fn compute_circle_circle_full(&mut self) -> Result<bool, String> {
        let c1 = self.curve1.clone().ok_or("EdgeEdge: edge1 has no curve")?;
        let c2 = self.curve2.clone().ok_or("EdgeEdge: edge2 has no curve")?;
        let (cc1, r1, n1) = circle_of_curve(&*c1).ok_or("EdgeEdge: edge1 is not a circle")?;
        let (cc2, r2, n2) = circle_of_curve(&*c2).ok_or("EdgeEdge: edge2 is not a circle")?;
        let tol = self.tol.max(1e-9);

        self.common_parts.clear();
        self.points.clear();

        // Non-coplanar circles → general (sampling/exact) solver.
        if n1.crossed(&n2).magnitude() > 1e-6 {
            return self.fallback_general();
        }

        // Coincident (same center & radius) → overlapping arc as a common part.
        if cc1.distance(&cc2) <= tol && (r1 - r2).abs() <= tol {
            if self.is_coincident() {
                self.push_coincident_common_part();
                return Ok(true);
            }
            return Ok(false);
        }

        let dist = GpVec::from_pnts(&cc1, &cc2).magnitude();
        if dist <= tol {
            return Ok(false); // concentric, distinct radii
        }
        let rsum = r1 + r2;
        let rdiff = (r1 - r2).abs();
        if dist > rsum + tol || dist < rdiff - tol {
            return Ok(false); // external or internal separation
        }

        // Radical line: one or two intersection points in the shared plane.
        let u = GpVec::from_pnts(&cc1, &cc2).divided(dist);
        let x = (r1 * r1 - r2 * r2 + dist * dist) / (2.0 * dist);
        let h2 = r1 * r1 - x * x;
        let h = if h2 > 0.0 { h2.sqrt() } else { 0.0 };
        let w = u.crossed(&n1);
        let base = cc1.translated_vec(&u.multiplied_scalar(x));
        let mut pts = vec![base.translated_vec(&w.multiplied_scalar(h))];
        if h > tol {
            pts.push(base.translated_vec(&w.multiplied_scalar(-h)));
        }

        let mut found = false;
        for p in pts {
            if let Some(hit) = self.point_hit_on_both(&p) {
                self.points.push(hit);
                found = true;
            }
        }
        Ok(found)
    }

    /// Merged output of every intersection result, for external consumers.
    ///
    /// Discrete vertex hits appear as `(t1, t2, point)`; each coincident common
    /// part contributes its two boundary points `(t1, t2, point)` (the edge-2
    /// parameter comes from projecting the boundary point onto curve 2).
    pub fn intersection_points(&self) -> Vec<(f64, f64, GpPnt)> {
        let mut out: Vec<(f64, f64, GpPnt)> = Vec::new();
        for p in &self.points {
            out.push((p.uv1.0, p.uv2.0, p.pnt1));
        }
        for cp in &self.common_parts {
            out.extend(self.common_part_span(cp));
        }
        out
    }

    /// The coincident parameter range on edge 1, or `None` when the two curves
    /// do not coincide (overlap) on their configured ranges.
    ///
    /// The projection-based [`coincident_ranges`](Self::coincident_ranges) only
    /// clips parameter windows, so the overlap is additionally verified by
    /// checking the 3D distance at its mid-point — parallel-but-distinct curves
    /// (whose projections overlap) are thereby excluded.
    pub fn coincident_range(&self) -> Option<IntRange> {
        let (a, b, _, _) = self.coincident_ranges()?;
        let c1 = self.curve1.as_ref()?;
        let c2 = self.curve2.as_ref()?;
        let tol = self.tol.max(1e-9);
        let mid = 0.5 * (a + b);
        let p1 = c1.d0(self.edge_to_curve_param(1, mid));
        let pr = geom_api::project_point_on_curve(&**c2, &p1, tol)?;
        if pr.distance > tol {
            return None;
        }
        Some(IntRange::new_unchecked(a, b))
    }

    // ---- helpers for the full solvers ----

    /// General fallback for non-coplanar circles (the exact/sampling solver).
    pub(super) fn fallback_general(&mut self) -> Result<bool, String> {
        let hits = self.intersect_edges();
        if hits.is_empty() {
            return Ok(false);
        }
        self.push_hits(hits);
        Ok(true)
    }

    /// Map a 3D point lying on both curves to a `PntOn2Faces` hit, projecting
    /// for the parameter on each curve and checking the edge ranges.
    pub(super) fn point_hit_on_both(&self, p: &GpPnt) -> Option<PntOn2Faces> {
        let c1 = self.curve1.as_ref()?;
        let c2 = self.curve2.as_ref()?;
        let tol = self.tol.max(1e-9);
        let pr1 = geom_api::project_point_on_curve(&**c1, p, tol)?;
        let pr2 = geom_api::project_point_on_curve(&**c2, p, tol)?;
        if pr1.distance > tol || pr2.distance > tol {
            return None;
        }
        let et1 = self.curve_to_edge_param(1, pr1.parameter);
        let et2 = self.curve_to_edge_param(2, pr2.parameter);
        if self.r1().contains(et1) && self.r2().contains(et2) {
            Some(PntOn2Faces::new(0, 1, *p, *p, (et1, 0.0), (et2, 0.0)))
        } else {
            None
        }
    }

    /// The two boundary points of a common part as `(t1, t2, point)` entries.
    pub(super) fn common_part_span(&self, cp: &CommonPrt) -> Vec<(f64, f64, GpPnt)> {
        let c1 = match &self.curve1 { Some(c) => c, None => return Vec::new() };
        let c2 = match &self.curve2 { Some(c) => c, None => return Vec::new() };
        let tol = self.tol.max(1e-9);
        let r = cp.range;
        let mut out = Vec::new();
        for t1 in [r.first, r.last] {
            let u1 = self.edge_to_curve_param(1, t1);
            let p = c1.d0(u1);
            if let Some(pr) = geom_api::project_point_on_curve(&**c2, &p, tol) {
                let t2 = self.curve_to_edge_param(2, pr.parameter);
                out.push((t1, t2, p));
            }
        }
        out
    }
}

/// Closest points between two infinite lines — `(t on line1, s on line2, point)`.
/// `None` when the lines are parallel. Standard closest-approach solve of the
/// two-parameter linear system (with `r = p1 - p2`).
pub(super) fn line_line_closest(l1: &GpLin, l2: &GpLin) -> Option<(f64, f64, GpPnt)> {
    let p1 = l1.location();
    let p2 = l2.location();
    let d1 = GpVec::from_xyz(l1.direction().xyz());
    let d2 = GpVec::from_xyz(l2.direction().xyz());
    let r = GpVec::from_pnts(&p2, &p1); // p1 - p2
    let a = d1.dot(&d1);
    let b = d1.dot(&d2);
    let c = d1.dot(&r);
    let e = d2.dot(&d2);
    let f = d2.dot(&r);
    let denom = a * e - b * b;
    if denom.abs() <= 1e-30 {
        return None; // parallel
    }
    let t = (b * f - c * e) / denom;
    let s = (a * f - b * c) / denom;
    Some((t, s, p1.translated_vec(&d1.multiplied_scalar(t))))
}
