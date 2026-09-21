use super::prelude::*;
use super::*;

/// The analytic family of an edge curve, used to choose the intersection path.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]

pub(super) enum CurveType {
    /// Geometrically a straight line (constant tangent direction).
    Line,
    /// Geometrically a planar circle (coplanar, equidistant samples).
    Circle,
    /// Anything else (BSpline, Bezier, ellipse, …).
    Other,
}

/// Line/line outcome of `IntTools_EdgeEdge::ComputeLineLine`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum LineLineKind {
    /// Overlapping coincident ranges → `TopAbs_EDGE` common part.
    Coincide,
    /// Parallel distinct, or coincident 3d line with no range overlap.
    Empty,
    /// Non-parallel intersection that may produce a vertex hit.
    Crossing,
}

/// Edge/edge intersection algorithm.
///
/// Mirrors the public surface of `IntTools_EdgeEdge`: set the two edges (with
/// optional parameter ranges), run [`perform`](Self::perform), then read the
/// common parts and/or vertex points.
pub struct EdgeEdge {
    pub(super) edge1: Option<Edge>,
    pub(super) edge2: Option<Edge>,
    pub(super) range1: Option<IntRange>,
    pub(super) range2: Option<IntRange>,
    pub(super) fuzzy: f64,
    /// `IntTools_EdgeEdge::myQuickCoincidenceCheck`.
    pub(super) quick_coincidence_check: bool,
    // Prepared geometry (filled by `prepare`).
    pub(super) curve1: Option<Arc<dyn Curve>>,
    pub(super) curve2: Option<Arc<dyn Curve>>,
    pub(super) ctype1: CurveType,
    pub(super) ctype2: CurveType,
    pub(super) tol1: f64,
    pub(super) tol2: f64,
    pub(super) tol: f64,
    // Results.
    pub(super) done: bool,
    pub(super) common_parts: Vec<CommonPrt>,
    pub(super) points: Vec<PntOn2Faces>,
}

impl EdgeEdge {
    /// Empty algorithm. Edges must be set with [`set_edge1`](Self::set_edge1) /
    /// [`set_edge2`](Self::set_edge2) before [`perform`](Self::perform).
    pub fn new() -> Self {
        Self {
            edge1: None,
            edge2: None,
            range1: None,
            range2: None,
            fuzzy: 0.0,
            quick_coincidence_check: false,
            curve1: None,
            curve2: None,
            ctype1: CurveType::Other,
            ctype2: CurveType::Other,
            tol1: 0.0,
            tol2: 0.0,
            tol: 0.0,
            done: false,
            common_parts: Vec::new(),
            points: Vec::new(),
        }
    }

    /// Construct with both edges already set (full ranges).
    pub fn with_edges(e1: Edge, e2: Edge) -> Self {
        let mut s = Self::new();
        s.set_edge1(e1);
        s.set_edge2(e2);
        s
    }

    /// Set the first edge (its full parameter range is used unless
    /// [`set_range1`](Self::set_range1) is called).
    pub fn set_edge1(&mut self, e: Edge) {
        self.edge1 = Some(e);
    }

    /// Set the second edge.
    pub fn set_edge2(&mut self, e: Edge) {
        self.edge2 = Some(e);
    }

    /// Restrict the first edge to `r` (in the edge's parameter space).
    pub fn set_range1(&mut self, r: IntRange) {
        self.range1 = Some(r);
    }

    /// Restrict the second edge to `r` (in the edge's parameter space).
    pub fn set_range2(&mut self, r: IntRange) {
        self.range2 = Some(r);
    }

    /// Fuzzy tolerance added (half to each edge) on top of the edge
    /// tolerances. Mirrors `SetFuzzyValue`.
    pub fn set_fuzzy_value(&mut self, v: f64) {
        self.fuzzy = v;
    }

    /// `IntTools_EdgeEdge::UseQuickCoincidenceCheck`.
    pub fn set_quick_coincidence_check(&mut self, b: bool) {
        self.quick_coincidence_check = b;
    }

    /// Run the intersection. `Err` when the input edges are missing or carry
    /// no curve geometry.
    pub fn perform(&mut self) -> Result<(), String> {
        self.common_parts.clear();
        self.points.clear();
        self.done = false;
        self.prepare()?;
        // OCCT `IntTools_EdgeEdge::Perform`: Line/Line is handled first and
        // returns; the quick-coincidence shortcut then applies to every other
        // pair (`IntTools_EdgeEdge.cxx:198-214`).
        if self.ctype1 == CurveType::Line && self.ctype2 == CurveType::Line {
            self.compute_line_line();
            self.done = true;
            return Ok(());
        }
        if self.quick_coincidence_check && self.is_coincident() {
            self.push_coincident_common_part();
            self.done = true;
            return Ok(());
        }
        match (self.ctype1, self.ctype2) {
            (CurveType::Circle, CurveType::Circle) => self.compute_circle_circle(),
            _ => self.find_solutions(),
        }
        self.done = true;
        Ok(())
    }

    /// True when [`perform`](Self::perform) succeeded.
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// The resolved range on the first edge (always set after `prepare`).
    pub(super) fn r1(&self) -> IntRange {
        self.range1.unwrap_or_default()
    }

    /// The resolved range on the second edge.
    pub(super) fn r2(&self) -> IntRange {
        self.range2.unwrap_or_default()
    }

    /// The common parts found (coincident overlapping ranges).
    pub fn common_parts(&self) -> &[CommonPrt] {
        &self.common_parts
    }

    /// The discrete intersection points found, as `(parameter on edge 1,
    /// parameter on edge 2)` pairs stored in the UV slots.
    pub fn points(&self) -> &[PntOn2Faces] {
        &self.points
    }

    // -----------------------------------------------------------------------
    // Preparation
    // -----------------------------------------------------------------------

    /// Resolve ranges, classify the curves and compute tolerances.
    ///
    /// Mirrors `IntTools_EdgeEdge::Prepare`: the (0,0) range is treated as the
    /// "unset" sentinel and replaced by the edge's full parameter range; the
    /// total tolerance is `tol1 + tol2` with half the fuzzy value added to each
    /// edge.
    pub(super) fn prepare(&mut self) -> Result<(), String> {
        let e1 = self.edge1.as_ref().cloned().ok_or("EdgeEdge: edge1 not set")?;
        let e2 = self.edge2.as_ref().cloned().ok_or("EdgeEdge: edge2 not set")?;
        self.curve1 = BRepTool::edge_curve_world(&e1);
        self.curve2 = BRepTool::edge_curve_world(&e2);
        let c1 = self.curve1.as_ref().ok_or("EdgeEdge: edge1 has no curve")?.clone();
        let c2 = self.curve2.as_ref().ok_or("EdgeEdge: edge2 has no curve")?.clone();

        let (f1, l1) = BRepTool::edge_parameters(&e1);
        let (f2, l2) = BRepTool::edge_parameters(&e2);
        let r1 = self.range1.unwrap_or_else(|| IntRange::new_unchecked(f1, l1));
        let r2 = self.range2.unwrap_or_else(|| IntRange::new_unchecked(f2, l2));
        let r1 = if r1.first == 0.0 && r1.last == 0.0 { IntRange::new_unchecked(f1, l1) } else { r1 };
        let r2 = if r2.first == 0.0 && r2.last == 0.0 { IntRange::new_unchecked(f2, l2) } else { r2 };
        self.range1 = Some(r1);
        self.range2 = Some(r2);

        self.ctype1 = classify(&*c1);
        self.ctype2 = classify(&*c2);

        let add = self.fuzzy / 2.0;
        self.tol1 = BRepTool::edge_tolerance(&e1) + add;
        self.tol2 = BRepTool::edge_tolerance(&e2) + add;
        self.tol = self.tol1 + self.tol2;
        Ok(())
    }

    /// Map a curve-space parameter `t` into the edge's parameter space.
    ///
    /// Only curves reparametrized to `[0, 1]` (GeomTrimmedCurve) differ from
    /// the edge range here; everything else uses a one-to-one parameter space.
    pub(super) fn curve_to_edge_param(&self, _idx: usize, t: f64) -> f64 {
        // The edge parameter range is the curve parameter range (the solver's
        // `r1`/`r2` are sub-ranges *of the edge*, not a re-parameterization of
        // the curve), so a curve parameter is already an edge parameter.
        t
    }

    /// Map an edge-space parameter `u` into the curve's parameter space.
    pub(super) fn edge_to_curve_param(&self, _idx: usize, u: f64) -> f64 {
        // See `curve_to_edge_param`: the edge and curve share the parameter
        // space, so an edge parameter is already a curve parameter.
        u
    }

    // -----------------------------------------------------------------------
    // Analytic branches
    // -----------------------------------------------------------------------

    /// Line/line case: a coincident overlap becomes a common part, otherwise
    /// the exact segment solver reports the single vertex hit.
    ///
    /// Coincidence is the analytic `IntTools_EdgeEdge::ComputeLineLine` test
    /// (parallel + distance, or both endpoints of edge 1 on line 2), not the
    /// sampling `IsCoincident` used for general curves. Collinear ranges that
    /// do not overlap return empty — they must not fall through to a vertex hit
    /// at a shared endpoint (`IntTools_EdgeEdge.cxx:902-993`).
    pub(super) fn compute_line_line(&mut self) {
        match self.line_line_kind() {
            LineLineKind::Coincide => self.push_coincident_common_part(),
            LineLineKind::Empty => {}
            LineLineKind::Crossing => self.push_hits(self.intersect_edges()),
        }
    }

    /// Analytic line/line classification from `ComputeLineLine`.
    pub(super) fn line_line_kind(&self) -> LineLineKind {
        let Some(c1) = self.curve1.as_ref() else {
            return LineLineKind::Empty;
        };
        let Some(c2) = self.curve2.as_ref() else {
            return LineLineKind::Empty;
        };
        let Some(l1) = line_of_curve(&**c1) else {
            return LineLineKind::Crossing;
        };
        let Some(l2) = line_of_curve(&**c2) else {
            return LineLineKind::Crossing;
        };
        let a_tol = self.tol * self.tol;
        let a_d1 = l1.direction();
        let a_d2 = l2.direction();
        let mut is_coincide = a_d1.angle(&a_d2) < ANGULAR;
        if is_coincide && l1.square_distance(&l2.location()) > a_tol {
            return LineLineKind::Empty;
        }
        let r1 = self.r1();
        let r2 = self.r2();
        let (a_t11, a_t12) = (r1.first, r1.last);
        let (a_t21, a_t22) = (r2.first, r2.last);
        let a_p11 = c1.d0(a_t11);
        let a_p12 = c1.d0(a_t12);
        if !is_coincide {
            let o2 = if a_t21.is_finite() && a_t22.is_finite() {
                c2.d0(0.5 * (a_t21 + a_t22))
            } else {
                l2.location()
            };
            let d2 = GpVec::from_xyz(a_d2.xyz());
            let a_vec1 = GpVec::from_pnts(&o2, &a_p11).crossed(&d2);
            let a_vec2 = GpVec::from_pnts(&o2, &a_p12).crossed(&d2);
            is_coincide =
                a_vec1.square_magnitude() <= a_tol && a_vec2.square_magnitude() <= a_tol;
            if !is_coincide && a_vec1.dot(&a_vec2) > 0.0 {
                return LineLineKind::Empty;
            }
        }
        if !is_coincide {
            return LineLineKind::Crossing;
        }
        let proj_tol = self.tol.max(1e-9);
        let Some(qa) = geom_api::project_point_on_curve(&**c2, &a_p11, proj_tol) else {
            return LineLineKind::Empty;
        };
        let Some(qb) = geom_api::project_point_on_curve(&**c2, &a_p12, proj_tol) else {
            return LineLineKind::Empty;
        };
        let (t21, t22) = (qa.parameter, qb.parameter);
        if (t21 > a_t22 && t22 > a_t22) || (t21 < a_t21 && t22 < a_t21) {
            return LineLineKind::Empty;
        }
        LineLineKind::Coincide
    }

    /// Circle/circle case: coincident circles become a common part, otherwise
    /// the radical-line solver reports the up-to-two vertex hits.
    pub(super) fn compute_circle_circle(&mut self) {
        if self.is_coincident() {
            self.push_coincident_common_part();
            return;
        }
        self.push_hits(self.intersect_edges());
    }

    /// Exact (or sampling) intersection of the two edges restricted to the
    /// configured ranges, reusing `crate::inttools::edge_edge_intersections`.
    pub(super) fn intersect_edges(&self) -> Vec<EdgeEdgeHit> {
        let e1 = match &self.edge1 { Some(e) => e, None => return Vec::new() };
        let e2 = match &self.edge2 { Some(e) => e, None => return Vec::new() };
        let r1 = self.r1();
        let r2 = self.r2();
        edge_edge_intersections(e1, e2, self.tol.max(1e-9))
            .into_iter()
            .filter(|h| r1.contains(h.u1) && r2.contains(h.u2))
            .collect()
    }

    pub(super) fn push_hits(&mut self, hits: Vec<EdgeEdgeHit>) {
        for h in hits {
            self.points.push(PntOn2Faces::new(0, 1, h.point, h.point, (h.u1, 0.0), (h.u2, 0.0)));
        }
    }

    // -----------------------------------------------------------------------
    // General curves
    // -----------------------------------------------------------------------

    /// General-curve case: the extrema of the distance; a pair whose distance is
    /// within tolerance is an intersection.
    ///
    /// **UNPORTED (audit A11 / task T-47 sub-item 3)** — OCCT's
    /// `IntTools_EdgeEdge::FindSolutions` does not enumerate `Extrema_ExtCC`
    /// solutions: it recurses on the parameter boxes
    /// (`IntTools_EdgeEdge.cxx:290-549`) with `BndBuildBox` (`:1410-1419`),
    /// `FindParameters` (`:553-671`), `IsIntersection` (`:1060-1146`),
    /// `CheckCoincidence` (`:1150-1206`) and `SplitRangeOnSegments`
    /// (`:1366-1406`), so that *every* crossing whose boxes stay in contact is
    /// reported even when it is not a distance extremum of the whole range. The
    /// port's substitute is the "extrema within tolerance" criterion above; the
    /// previous body additionally merged the sampled
    /// `inttools::edge_edge_intersections` solutions, which has no OCCT
    /// counterpart at all and was removed (the faithful extrema engine
    /// `Extrema_ExtCC`/`GGenExtCC` covers those cases, see T-66).
    pub(super) fn find_solutions(&mut self) {
        if self.is_coincident() {
            self.push_coincident_common_part();
            return;
        }
        let c1 = self.curve1.clone().unwrap();
        let c2 = self.curve2.clone().unwrap();
        let tol = self.tol;
        let r1 = self.r1();
        let r2 = self.r2();

        let mut solutions: Vec<(f64, f64, GpPnt)> = Vec::new();
        // `Extrema_ExtCC` runs over the *edge* ranges: OCCT's `IntTools_EdgeEdge`
        // holds `BRepAdaptor_Curve` objects, whose `FirstParameter`/`LastParameter`
        // are the edge's (`IntTools_EdgeEdge.cxx:94-95`, `:164-165`), and
        // `Extrema_ExtCC` forwards them to the engine (`Extrema_ExtCC.cxx:180`).
        // The curve's own range may be infinite (an edge on an infinite line).
        for p in curve_curve_extrema_all_range(
            &*c1,
            &*c2,
            r1.first,
            r1.last,
            r2.first,
            r2.last,
        ) {
            if p.distance <= tol.max(1e-7) {
                if let Some(sol) = self.find_parameters(p.u1, p.u2) {
                    solutions.push(sol);
                }
            }
        }
        self.merge_solutions(solutions);
    }

    /// Polish a candidate `(u1, u2)` parameter pair into an intersection.
    ///
    /// The seed already comes from a local extremum of the distance; Newton
    /// polish (`locate_extcc`) sharpens it, and the seed itself is accepted when
    /// the polish strays. Returns `(edge1 param, edge2 param, point)`.
    pub(super) fn find_parameters(&self, u1: f64, u2: f64) -> Option<(f64, f64, GpPnt)> {
        let c1 = self.curve1.as_ref()?;
        let c2 = self.curve2.as_ref()?;
        let tol = self.tol.max(1e-7);
        let (u1r, u2r, p1) = match locate_extcc(&**c1, &**c2, u1, u2) {
            Some(e) if e.distance <= tol => (e.u1, e.u2, e.p1),
            _ => {
                let p = c1.d0(u1);
                if p.distance(&c2.d0(u2)) <= tol {
                    (u1, u2, p)
                } else {
                    return None;
                }
            }
        };
        let eu1 = self.curve_to_edge_param(1, u1r);
        let eu2 = self.curve_to_edge_param(2, u2r);
        if self.r1().contains(eu1) && self.r2().contains(eu2) {
            Some((eu1, eu2, p1))
        } else {
            None
        }
    }

    /// Deduplicate and sort solutions by the first parameter, then store them
    /// as [`PntOn2Faces`]. Reuses `inttools_roots::remove_identical_roots` and
    /// `sort_roots` on the primary parameter.
    pub(super) fn merge_solutions(&mut self, solutions: Vec<(f64, f64, GpPnt)>) {
        if solutions.is_empty() {
            return;
        }
        let eps = self.tol.max(1e-7);
        let mut roots: Vec<IntRoot> = solutions
            .iter()
            .enumerate()
            .map(|(i, s)| IntRoot::new(i as i32, RootType::IsRoot, IntRange::new_unchecked(s.0, s.0)))
            .collect();
        remove_identical_roots(&mut roots, eps);
        sort_roots(&mut roots);
        let mut out: Vec<(f64, f64, GpPnt)> = Vec::new();
        for r in roots {
            let s = solutions[r.root_index() as usize];
            let dup = out
                .iter()
                .any(|o| (o.0 - s.0).abs() <= eps && o.2.distance(&s.2) <= eps);
            if !dup {
                out.push(s);
            }
        }
        for (t1, t2, p) in out {
            self.points.push(PntOn2Faces::new(0, 1, p, p, (t1, 0.0), (t2, 0.0)));
        }
    }

    // -----------------------------------------------------------------------
    // Coincidence
    // -----------------------------------------------------------------------

    /// Whether the two curves coincide on their configured ranges: more than
    /// half of 24 samples of edge 1 project onto edge 2 within tolerance.
    ///
    /// Port of `IntTools_EdgeEdge::IsCoincident`.
    pub(super) fn is_coincident(&self) -> bool {
        let c1 = match &self.curve1 { Some(c) => c, None => return false };
        let c2 = match &self.curve2 { Some(c) => c, None => return false };
        let tol = self.tol.max(1e-9);
        let r1 = self.r1();
        let r2 = self.r2();
        let (a11, a12) = (r1.first, r1.last);
        let (a21, a22) = (r2.first, r2.last);
        let t2lo = self.edge_to_curve_param(2, a21);
        let t2hi = self.edge_to_curve_param(2, a22);
        let n = 24usize;
        let mut cnt = 0usize;
        for i in 0..=n {
            let u1 = a11 + (a12 - a11) * i as f64 / n as f64;
            let p1 = c1.d0(self.edge_to_curve_param(1, u1));
            if let Some(pr) = geom_api::project_point_on_curve(&**c2, &p1, tol) {
                if pr.distance < tol && pr.parameter >= t2lo - tol && pr.parameter <= t2hi + tol {
                    cnt += 1;
                }
            }
        }
        cnt as f64 / (n + 1) as f64 > 0.5
    }

    /// The overlapping sub-range of two coincident edges.
    ///
    /// Returns `(t11, t12, t21, t22)` — the overlap on edge 1 and on edge 2 in
    /// edge-parameter space — or `None` when the ranges do not overlap.
    pub fn coincident_ranges(&self) -> Option<(f64, f64, f64, f64)> {
        let c1 = self.curve1.as_ref()?;
        let c2 = self.curve2.as_ref()?;
        let tol = self.tol.max(1e-9);
        let r1 = self.r1();
        let r2 = self.r2();
        let (a11, a12) = (r1.first, r1.last);
        let (a21, a22) = (r2.first, r2.last);
        // Endpoints of edge 1's range projected onto curve 2.
        let p1a = c1.d0(self.edge_to_curve_param(1, a11));
        let p1b = c1.d0(self.edge_to_curve_param(1, a12));
        let qa = geom_api::project_point_on_curve(&**c2, &p1a, tol)?;
        let qb = geom_api::project_point_on_curve(&**c2, &p1b, tol)?;
        let (lo2, hi2) = (qa.parameter.min(qb.parameter), qa.parameter.max(qb.parameter));
        // Clip to edge 2's range (curve space).
        let t2lo = self.edge_to_curve_param(2, a21);
        let t2hi = self.edge_to_curve_param(2, a22);
        let (clo, chi) = (lo2.max(t2lo), hi2.min(t2hi));
        if chi - clo <= tol {
            return None;
        }
        // Endpoints of the overlap, projected back onto curve 1.
        let qlo = c2.d0(clo);
        let qhi = c2.d0(chi);
        let rlo = geom_api::project_point_on_curve(&**c1, &qlo, tol)?;
        let rhi = geom_api::project_point_on_curve(&**c1, &qhi, tol)?;
        let (t1lo, t1hi) = (rlo.parameter.min(rhi.parameter), rlo.parameter.max(rhi.parameter));
        Some((
            self.curve_to_edge_param(1, t1lo),
            self.curve_to_edge_param(1, t1hi),
            self.curve_to_edge_param(2, clo),
            self.curve_to_edge_param(2, chi),
        ))
    }

    /// Store the overlapping range of coincident edges as a [`CommonPrt`].
    pub(super) fn push_coincident_common_part(&mut self) {
        if let Some((a, b, _, _)) = self.coincident_ranges() {
            if b - a > self.tol.max(1e-9) {
                self.common_parts.push(CommonPrt::with(
                    CommonPartType::Edge,
                    IntRange::new_unchecked(a, b),
                    None,
                    Vec::new(),
                ));
            }
        }
    }
}

impl Default for EdgeEdge {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Curve classification
// ---------------------------------------------------------------------------

/// Classify a curve as line / circle / other from geometric invariants.
pub(super) fn classify(c: &dyn Curve) -> CurveType {
    if line_of_curve(c).is_some() {
        CurveType::Line
    } else if circle_of_curve(c).is_some() {
        CurveType::Circle
    } else {
        CurveType::Other
    }
}

/// Extract a `GpLin` when the curve has a constant tangent direction (a real
/// line, possibly a straight BSpline). The unbounded parameter range of
/// `GeomLine` alone is not a sufficient test (hyperbolas are unbounded too),
/// so the direction is verified by sampling.
pub(super) fn line_of_curve(c: &dyn Curve) -> Option<GpLin> {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    let u0 = if a.is_finite() { a } else { 0.0 };
    let p = c.d0(u0);
    let d = c.d1(u0).1;
    let m = d.magnitude();
    if m < 1e-12 {
        return None;
    }
    if a.is_finite() && b.is_finite() && b > a {
        for i in 1..=6 {
            let t = c.d1(a + (b - a) * i as f64 / 6.0).1;
            if t.cross_magnitude(&d) > 1e-6 * m * t.magnitude().max(m) {
                return None;
            }
        }
    }
    occt_core::gp::GpDir::from_vec(&d).ok().map(|dir| GpLin::from_pnt_dir(p, dir))
}

/// Reconstruct `(center, radius, unit plane normal)` when every sample of the
/// bounded curve is coplanar and equidistant from a single center.
pub(super) fn circle_of_curve(c: &dyn Curve) -> Option<(GpPnt, f64, GpVec)> {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !(a.is_finite() && b.is_finite() && (b - a).abs() > 1e-15) {
        return None;
    }
    let n = 8;
    let pts: Vec<GpPnt> = (0..=n).map(|i| c.d0(a + (b - a) * i as f64 / n as f64)).collect();
    let (p0, p1, p2) = first_three_spanning(&pts)?;
    let center = circumcenter(&p0, &p1, &p2)?;
    let radius = p0.distance(&center);
    if radius <= 1e-30 {
        return None;
    }
    let nrm = GpVec::from_pnts(&p0, &p1).crossed(&GpVec::from_pnts(&p0, &p2));
    let nm = nrm.magnitude();
    if nm <= 1e-30 {
        return None;
    }
    let nrm = nrm.divided(nm);
    let scale = radius.max(1.0).max((b - a).abs());
    let tol = 1e-6 * scale;
    for p in pts {
        let v = GpVec::from_pnts(&center, &p);
        if v.dot(&nrm).abs() > tol {
            return None;
        }
        if (v.magnitude() - radius).abs() > tol {
            return None;
        }
    }
    Some((center, radius, nrm))
}

/// First three samples that span a plane (non-collinear), if any.
pub(super) fn first_three_spanning(pts: &[GpPnt]) -> Option<(GpPnt, GpPnt, GpPnt)> {
    let p0 = pts[0];
    let mut i1 = None;
    for (i, p) in pts.iter().enumerate().skip(1) {
        if GpVec::from_pnts(&p0, p).magnitude() > 1e-9 {
            i1 = Some(i);
            break;
        }
    }
    let i1 = i1?;
    let p1 = pts[i1];
    let d0 = GpVec::from_pnts(&p0, &p1);
    for p in pts.iter().skip(i1 + 1) {
        if GpVec::from_pnts(&p0, p).cross_magnitude(&d0) > 1e-9 * d0.magnitude().max(1e-9) {
            return Some((p0, p1, *p));
        }
    }
    None
}

pub(super) fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

pub(super) fn solve3(a: &[[f64; 3]; 3], rhs: &[f64; 3]) -> Option<[f64; 3]> {
    let d = det3(a);
    if d.abs() < 1e-20 {
        return None;
    }
    let mut x = [0.0; 3];
    for k in 0..3 {
        let mut m = *a;
        for i in 0..3 {
            m[i][k] = rhs[i];
        }
        x[k] = det3(&m) / d;
    }
    Some(x)
}

/// Circumcenter of three non-collinear points (perpendicular-bisector system).
pub(super) fn circumcenter(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> Option<GpPnt> {
    let d1 = GpVec::from_pnts(a, b);
    let d2 = GpVec::from_pnts(a, c);
    let n = d1.crossed(&d2);
    if n.magnitude() < 1e-30 {
        return None;
    }
    let n2 = |p: &GpPnt| p.coord.dot(&p.coord);
    let mat = [
        [d1.xyz().x, d1.xyz().y, d1.xyz().z],
        [d2.xyz().x, d2.xyz().y, d2.xyz().z],
        [n.xyz().x, n.xyz().y, n.xyz().z],
    ];
    let rhs = [0.5 * (n2(b) - n2(a)), 0.5 * (n2(c) - n2(a)), a.coord.dot(&n.xyz())];
    let o = solve3(&mat, &rhs)?;
    Some(GpPnt::new(o[0], o[1], o[2]))
}
