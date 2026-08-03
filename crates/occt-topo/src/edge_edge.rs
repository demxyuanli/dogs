//! Edge–edge intersection — port of `IntTools_EdgeEdge`.
//!
//! Orchestrates the intersection of two edges into two kinds of results:
//!
//! - **discrete vertex hits** — [`PntOn2Faces`] entries carrying the parameter
//!   on each edge and the 3D point (`points()`);
//! - **coincident common parts** — [`CommonPrt`] entries with [`CommonPartType::Edge`]
//!   describing an overlapping sub-range of two (nearly) coincident edges
//!   (`common_parts()`).
//!
//! Line–line and circle–circle are dispatched to the exact solver in
//! `crate::inttools::edge_edge_intersections` (segment–segment, radical line);
//! every other combination — BSpline/Bezier/general curves — is routed through
//! `occt_geom::extrema_cc::curve_curve_extrema_all`, whose zero-distance local
//! extrema are the intersection points, with the exact/sampling solver as a
//! complement. Coincidence is detected by sampling one curve and projecting the
//! samples onto the other (port of `IntTools_EdgeEdge::IsCoincident`).
//!
//! The task spec suggested `geom2d_api::project_point_on_curve` for parameter
//! refinement; that is a 2-D API, so the 3-D analog `geom_api::project_point_on_curve`
//! (and the extrema Newton polish `extrema_cc::locate_extcc`) is used instead.

use std::cmp::Ordering;
use std::sync::Arc;

use occt_core::gp::{GpLin, GpPnt, GpVec};
use occt_geom::extrema_cc::{curve_curve_extrema_all, locate_extcc};
use occt_geom::geom_api;
use occt_geom::Curve;

use crate::brep_tool::BRepTool;
use crate::inttools::{edge_edge_intersections, EdgeEdgeHit};
use crate::inttools_data::{CommonPartType, CommonPrt, IntRange, IntRoot, PntOn2Faces, RootType};
use crate::inttools_roots::{remove_identical_roots, sort_roots};
use crate::shape::Edge;

/// The analytic family of an edge curve, used to choose the intersection path.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CurveType {
    /// Geometrically a straight line (constant tangent direction).
    Line,
    /// Geometrically a planar circle (coplanar, equidistant samples).
    Circle,
    /// Anything else (BSpline, Bezier, ellipse, …).
    Other,
}

/// Edge/edge intersection algorithm.
///
/// Mirrors the public surface of `IntTools_EdgeEdge`: set the two edges (with
/// optional parameter ranges), run [`perform`](Self::perform), then read the
/// common parts and/or vertex points.
pub struct EdgeEdge {
    edge1: Option<Edge>,
    edge2: Option<Edge>,
    range1: Option<IntRange>,
    range2: Option<IntRange>,
    fuzzy: f64,
    // Prepared geometry (filled by `prepare`).
    curve1: Option<Arc<dyn Curve>>,
    curve2: Option<Arc<dyn Curve>>,
    ctype1: CurveType,
    ctype2: CurveType,
    tol1: f64,
    tol2: f64,
    tol: f64,
    // Results.
    done: bool,
    common_parts: Vec<CommonPrt>,
    points: Vec<PntOn2Faces>,
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

    /// Run the intersection. `Err` when the input edges are missing or carry
    /// no curve geometry.
    pub fn perform(&mut self) -> Result<(), String> {
        self.common_parts.clear();
        self.points.clear();
        self.done = false;
        self.prepare()?;
        match (self.ctype1, self.ctype2) {
            (CurveType::Line, CurveType::Line) => self.compute_line_line(),
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
    fn r1(&self) -> IntRange {
        self.range1.unwrap_or_default()
    }

    /// The resolved range on the second edge.
    fn r2(&self) -> IntRange {
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
    fn prepare(&mut self) -> Result<(), String> {
        let e1 = self.edge1.as_ref().cloned().ok_or("EdgeEdge: edge1 not set")?;
        let e2 = self.edge2.as_ref().cloned().ok_or("EdgeEdge: edge2 not set")?;
        self.curve1 = BRepTool::edge_curve(&e1);
        self.curve2 = BRepTool::edge_curve(&e2);
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
    fn curve_to_edge_param(&self, _idx: usize, t: f64) -> f64 {
        // The edge parameter range is the curve parameter range (the solver's
        // `r1`/`r2` are sub-ranges *of the edge*, not a re-parameterization of
        // the curve), so a curve parameter is already an edge parameter.
        t
    }

    /// Map an edge-space parameter `u` into the curve's parameter space.
    fn edge_to_curve_param(&self, _idx: usize, u: f64) -> f64 {
        // See `curve_to_edge_param`: the edge and curve share the parameter
        // space, so an edge parameter is already a curve parameter.
        u
    }

    // -----------------------------------------------------------------------
    // Analytic branches
    // -----------------------------------------------------------------------

    /// Line/line case: a coincident overlap becomes a common part, otherwise
    /// the exact segment solver reports the single vertex hit.
    fn compute_line_line(&mut self) {
        if self.is_coincident() {
            self.push_coincident_common_part();
            return;
        }
        self.push_hits(self.intersect_edges());
    }

    /// Circle/circle case: coincident circles become a common part, otherwise
    /// the radical-line solver reports the up-to-two vertex hits.
    fn compute_circle_circle(&mut self) {
        if self.is_coincident() {
            self.push_coincident_common_part();
            return;
        }
        self.push_hits(self.intersect_edges());
    }

    /// Exact (or sampling) intersection of the two edges restricted to the
    /// configured ranges, reusing `crate::inttools::edge_edge_intersections`.
    fn intersect_edges(&self) -> Vec<EdgeEdgeHit> {
        let e1 = match &self.edge1 { Some(e) => e, None => return Vec::new() };
        let e2 = match &self.edge2 { Some(e) => e, None => return Vec::new() };
        let r1 = self.r1();
        let r2 = self.r2();
        edge_edge_intersections(e1, e2, self.tol.max(1e-9))
            .into_iter()
            .filter(|h| r1.contains(h.u1) && r2.contains(h.u2))
            .collect()
    }

    fn push_hits(&mut self, hits: Vec<EdgeEdgeHit>) {
        for h in hits {
            self.points.push(PntOn2Faces::new(0, 1, h.point, h.point, (h.u1, 0.0), (h.u2, 0.0)));
        }
    }

    // -----------------------------------------------------------------------
    // General curves
    // -----------------------------------------------------------------------

    /// General-curve case: local extrema of the distance, plus the sampling
    /// solver as a complement; near-zero-distance pairs are the intersections.
    ///
    /// Mirrors `IntTools_EdgeEdge::FindSolutions` at the level this port needs:
    /// no bounding-box recursion, but the same "distance within tolerance is an
    /// intersection" criterion.
    fn find_solutions(&mut self) {
        if self.is_coincident() {
            self.push_coincident_common_part();
            return;
        }
        let c1 = self.curve1.clone().unwrap();
        let c2 = self.curve2.clone().unwrap();
        let e1 = self.edge1.clone().unwrap();
        let e2 = self.edge2.clone().unwrap();
        let tol = self.tol;
        let r1 = self.r1();
        let r2 = self.r2();

        let mut solutions: Vec<(f64, f64, GpPnt)> = Vec::new();
        for p in curve_curve_extrema_all(&*c1, &*c2) {
            if p.distance <= tol.max(1e-7) {
                if let Some(sol) = self.find_parameters(p.u1, p.u2) {
                    solutions.push(sol);
                }
            }
        }
        // Complement: the existing solver already handles line/circle/general
        // combinations; merge (and dedupe) whatever it reports.
        for h in edge_edge_intersections(&e1, &e2, tol.max(1e-9)) {
            if r1.contains(h.u1) && r2.contains(h.u2) {
                solutions.push((h.u1, h.u2, h.point));
            }
        }
        self.merge_solutions(solutions);
    }

    /// Polish a candidate `(u1, u2)` parameter pair into an intersection.
    ///
    /// The seed already comes from a local extremum of the distance; Newton
    /// polish (`locate_extcc`) sharpens it, and the seed itself is accepted when
    /// the polish strays. Returns `(edge1 param, edge2 param, point)`.
    fn find_parameters(&self, u1: f64, u2: f64) -> Option<(f64, f64, GpPnt)> {
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
    fn merge_solutions(&mut self, solutions: Vec<(f64, f64, GpPnt)>) {
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
    fn is_coincident(&self) -> bool {
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
    fn push_coincident_common_part(&mut self) {
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
fn classify(c: &dyn Curve) -> CurveType {
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
fn line_of_curve(c: &dyn Curve) -> Option<GpLin> {
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
fn circle_of_curve(c: &dyn Curve) -> Option<(GpPnt, f64, GpVec)> {
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
fn first_three_spanning(pts: &[GpPnt]) -> Option<(GpPnt, GpPnt, GpPnt)> {
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

fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

fn solve3(a: &[[f64; 3]; 3], rhs: &[f64; 3]) -> Option<[f64; 3]> {
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
fn circumcenter(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> Option<GpPnt> {
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;
    use std::sync::Arc;

    use occt_core::gp::{GpAx2, GpCirc, GpDir};
    use occt_geom::{GeomBSplineCurve, GeomCircle};

    use crate::builder::TopoBuilder;
    use crate::inttools::edge_edge_intersections;
    use crate::shape::TopoShape;
    use crate::tgeometry::GeometryRegistry;

    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&TopoShape::from_handle(c));
        }
    }

    fn dir(x: f64, y: f64, z: f64) -> GpDir {
        GpDir::new(x, y, z).expect("dir")
    }

    fn run(e1: &Edge, e2: &Edge, fuzzy: f64) -> EdgeEdge {
        let mut ee = EdgeEdge::new();
        ee.set_edge1(e1.clone());
        ee.set_edge2(e2.clone());
        ee.set_fuzzy_value(fuzzy);
        ee.perform().unwrap();
        ee
    }

    fn sorted_x(ee: &EdgeEdge) -> Vec<f64> {
        let mut xs: Vec<f64> = ee.points().iter().map(|p| p.pnt1.x()).collect();
        xs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
        xs
    }

    #[test]
    fn line_line_crossing_single_hit() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 2.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(0.0, 2.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let ee = run(&e1, &e2, 0.0);
        assert!(ee.is_done());
        assert!(ee.common_parts().is_empty(), "common: {:?}", ee.common_parts());
        let pts = ee.points();
        assert_eq!(pts.len(), 1, "points: {pts:?}");
        assert!(pts[0].pnt1.distance(&GpPnt::new(1.0, 1.0, 0.0)) < 1e-6);
        let d = 2.0f64.sqrt();
        assert!((pts[0].uv1.0 - d).abs() < 1e-6, "u1={}", pts[0].uv1.0);
        assert!((pts[0].uv2.0 - d).abs() < 1e-6, "u2={}", pts[0].uv2.0);
        clear_tree(&e1.0);
        clear_tree(&e2.0);
    }

    #[test]
    fn line_line_coincident_overlap_common_part() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(3.0, 0.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(4.0, 0.0, 0.0));
        let ee = run(&e1, &e2, 1e-7);
        assert_eq!(ee.common_parts().len(), 1, "common: {:?}", ee.common_parts());
        let cp = &ee.common_parts()[0];
        assert_eq!(cp.part_type, CommonPartType::Edge);
        assert!((cp.range.first - 1.0).abs() < 1e-6, "first={}", cp.range.first);
        assert!((cp.range.last - 3.0).abs() < 1e-6, "last={}", cp.range.last);
        clear_tree(&e1.0);
        clear_tree(&e2.0);
    }

    #[test]
    fn line_line_collinear_non_overlap_empty() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(2.0, 0.0, 0.0), &GpPnt::new(3.0, 0.0, 0.0));
        let ee = run(&e1, &e2, 1e-7);
        assert!(ee.points().is_empty(), "points: {:?}", ee.points());
        assert!(ee.common_parts().is_empty(), "common: {:?}", ee.common_parts());
        clear_tree(&e1.0);
        clear_tree(&e2.0);
    }

    #[test]
    fn line_line_parallel_distinct_empty() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(0.0, 1.0, 0.0), &GpPnt::new(2.0, 1.0, 0.0));
        let ee = run(&e1, &e2, 0.0);
        assert!(ee.points().is_empty());
        assert!(ee.common_parts().is_empty());
        clear_tree(&e1.0);
        clear_tree(&e2.0);
    }

    #[test]
    fn circle_circle_two_hits() {
        let b = TopoBuilder::new();
        let c1 = b.make_edge(Arc::new(GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0))), 0.0, 2.0 * PI);
        let ax2 = GpAx2::new(GpPnt::new(1.0, 0.0, 0.0), dir(0.0, 0.0, 1.0), dir(1.0, 0.0, 0.0)).unwrap();
        let c2 = b.make_edge(Arc::new(GeomCircle::new(GpCirc::new(ax2, 1.0))), 0.0, 2.0 * PI);
        let ee = run(&c1, &c2, 1e-7);
        assert!(ee.common_parts().is_empty(), "common: {:?}", ee.common_parts());
        let pts = ee.points();
        assert_eq!(pts.len(), 2, "points: {pts:?}");
        for p in pts {
            assert!((p.pnt1.x() - 0.5).abs() < 1e-6, "x={}", p.pnt1.x());
            assert!((p.pnt1.y().abs() - 0.75f64.sqrt()).abs() < 1e-4, "y={}", p.pnt1.y());
        }
        clear_tree(&c1.0);
        clear_tree(&c2.0);
    }

    #[test]
    fn separated_circles_empty() {
        let b = TopoBuilder::new();
        let c1 = b.make_edge(Arc::new(GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0))), 0.0, 2.0 * PI);
        let ax2 = GpAx2::new(GpPnt::new(5.0, 0.0, 0.0), dir(0.0, 0.0, 1.0), dir(1.0, 0.0, 0.0)).unwrap();
        let c2 = b.make_edge(Arc::new(GeomCircle::new(GpCirc::new(ax2, 1.0))), 0.0, 2.0 * PI);
        let ee = run(&c1, &c2, 1e-7);
        assert!(ee.points().is_empty(), "points: {:?}", ee.points());
        clear_tree(&c1.0);
        clear_tree(&c2.0);
    }

    #[test]
    fn line_bspline_crossing_two_hits() {
        let b = TopoBuilder::new();
        let bs = Arc::new(
            GeomBSplineCurve::new(
                vec![
                    GpPnt::new(0.0, 0.0, 0.0),
                    GpPnt::new(2.0, 2.0, 0.0),
                    GpPnt::new(4.0, 0.0, 0.0),
                ],
                vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
                2,
            )
            .unwrap(),
        );
        let bs_edge = b.make_edge(bs, 0.0, 1.0);
        let line = b.make_edge_segment(&GpPnt::new(0.0, 0.5, 0.0), &GpPnt::new(4.0, 0.5, 0.0));
        let ee = run(&line, &bs_edge, 1e-7);
        let pts = ee.points();
        assert_eq!(pts.len(), 2, "points: {pts:?}");
        let xs = sorted_x(&ee);
        // Quadratic Bezier (0,0)-(2,2)-(4,0): x=4t, y=4t(1-t); the line y=0.5
        // meets it at t=(1±1/√2)/2 → x=4t = 2∓√2.
        let xa = 2.0f64 - std::f64::consts::SQRT_2;
        let xb = 2.0f64 + std::f64::consts::SQRT_2;
        assert!((xs[0] - xa).abs() < 1e-3, "x0={} expected {xa}", xs[0]);
        assert!((xs[1] - xb).abs() < 1e-3, "x1={} expected {xb}", xs[1]);
        // The line parameter equals the x coordinate; the BSpline parameter
        // equals x/4.
        let mut u1s: Vec<f64> = pts.iter().map(|p| p.uv1.0).collect();
        u1s.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
        assert!((u1s[0] - xa).abs() < 1e-3, "u1={}", u1s[0]);
        let mut u2s: Vec<f64> = pts.iter().map(|p| p.uv2.0).collect();
        u2s.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
        assert!((u2s[0] - xa / 4.0).abs() < 1e-3, "u2={}", u2s[0]);
        assert!((u2s[1] - xb / 4.0).abs() < 1e-3, "u2={}", u2s[1]);
        clear_tree(&line.0);
        clear_tree(&bs_edge.0);
    }

    #[test]
    fn separated_edges_empty() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(0.0, 2.0, 0.0), &GpPnt::new(1.0, 2.0, 0.0));
        let ee = run(&e1, &e2, 0.0);
        assert!(ee.points().is_empty(), "points: {:?}", ee.points());
        assert!(ee.common_parts().is_empty());
        clear_tree(&e1.0);
        clear_tree(&e2.0);
    }

    #[test]
    fn matches_inttools_line_line() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 2.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(0.0, 2.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let ref_hits = edge_edge_intersections(&e1, &e2, 1e-9);
        let ee = run(&e1, &e2, 0.0);
        assert_eq!(ref_hits.len(), ee.points().len(), "ref={ref_hits:?} ee={:?}", ee.points());
        for (h, p) in ref_hits.iter().zip(ee.points()) {
            assert!((h.u1 - p.uv1.0).abs() < 1e-6, "u1 {} vs {}", h.u1, p.uv1.0);
            assert!((h.u2 - p.uv2.0).abs() < 1e-6, "u2 {} vs {}", h.u2, p.uv2.0);
            assert!(h.point.distance(&p.pnt1) < 1e-6);
        }
        clear_tree(&e1.0);
        clear_tree(&e2.0);
    }

    #[test]
    fn matches_inttools_circle_circle() {
        let b = TopoBuilder::new();
        let c1 = b.make_edge(Arc::new(GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0))), 0.0, 2.0 * PI);
        let ax2 = GpAx2::new(GpPnt::new(1.0, 0.0, 0.0), dir(0.0, 0.0, 1.0), dir(1.0, 0.0, 0.0)).unwrap();
        let c2 = b.make_edge(Arc::new(GeomCircle::new(GpCirc::new(ax2, 1.0))), 0.0, 2.0 * PI);
        let ref_hits = edge_edge_intersections(&c1, &c2, 1e-6);
        let ee = run(&c1, &c2, 1e-7);
        assert_eq!(ref_hits.len(), ee.points().len(), "ref={ref_hits:?} ee={:?}", ee.points());
        let mut ref_xs: Vec<f64> = ref_hits.iter().map(|h| h.point.x()).collect();
        let mut ee_xs: Vec<f64> = ee.points().iter().map(|p| p.pnt1.x()).collect();
        ref_xs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
        ee_xs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
        for (r, e) in ref_xs.iter().zip(&ee_xs) {
            assert!((r - e).abs() < 1e-6, "x {r} vs {e}");
        }
        clear_tree(&c1.0);
        clear_tree(&c2.0);
    }

    #[test]
    fn range_restriction_filters_hits() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 2.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(0.0, 2.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        // Restrict edge 1 to [1.9, 2.0]-ish parameters: the crossing at u=√2
        // is excluded.
        let mut ee = EdgeEdge::new();
        ee.set_edge1(e1.clone());
        ee.set_edge2(e2.clone());
        ee.set_range1(IntRange::new(2.0, 2.8).unwrap());
        ee.perform().unwrap();
        assert!(ee.points().is_empty(), "points: {:?}", ee.points());
        clear_tree(&e1.0);
        clear_tree(&e2.0);
    }
}

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
    fn fallback_general(&mut self) -> Result<bool, String> {
        let hits = self.intersect_edges();
        if hits.is_empty() {
            return Ok(false);
        }
        self.push_hits(hits);
        Ok(true)
    }

    /// Map a 3D point lying on both curves to a `PntOn2Faces` hit, projecting
    /// for the parameter on each curve and checking the edge ranges.
    fn point_hit_on_both(&self, p: &GpPnt) -> Option<PntOn2Faces> {
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
    fn common_part_span(&self, cp: &CommonPrt) -> Vec<(f64, f64, GpPnt)> {
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
fn line_line_closest(l1: &GpLin, l2: &GpLin) -> Option<(f64, f64, GpPnt)> {
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

// ---------------------------------------------------------------------------
// Phase 18c tests — the appended full solvers and best-solution selection.
// (Deliberately a separate test module so the original tests stay untouched.)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests_full {
    use super::*;
    use std::f64::consts::PI;
    use std::sync::Arc;

    use occt_core::gp::{GpAx2, GpCirc, GpDir};
    use occt_geom::GeomCircle;

    use crate::builder::TopoBuilder;
    use crate::shape::{Edge, TopoShape};
    use crate::tgeometry::GeometryRegistry;

    fn clear_tree(s: &TopoShape) {
        GeometryRegistry::global().clear_shape(s);
        let children = s.tshape.read().unwrap().children.clone();
        for c in children {
            clear_tree(&TopoShape::from_handle(c));
        }
    }

    fn dir(x: f64, y: f64, z: f64) -> GpDir {
        GpDir::new(x, y, z).expect("dir")
    }

    /// Set two edges and run `prepare`, so the full solvers can run directly.
    fn prepared_ee(e1: &Edge, e2: &Edge, fuzzy: f64) -> EdgeEdge {
        let mut ee = EdgeEdge::new();
        ee.set_edge1(e1.clone());
        ee.set_edge2(e2.clone());
        ee.set_fuzzy_value(fuzzy);
        ee.prepare().unwrap();
        ee
    }

    #[test]
    fn line_line_full_coincident_segment() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(3.0, 0.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(1.0, 0.0, 0.0), &GpPnt::new(4.0, 0.0, 0.0));
        let mut ee = prepared_ee(&e1, &e2, 1e-7);
        let found = ee.compute_line_line_full().unwrap();
        assert!(found, "coincident lines intersect");
        assert_eq!(ee.common_parts().len(), 1, "common: {:?}", ee.common_parts());
        // coincident_range covers the overlapping domain [1, 3].
        let cr = ee.coincident_range().unwrap();
        assert!((cr.first - 1.0).abs() < 1e-6, "first={}", cr.first);
        assert!((cr.last - 3.0).abs() < 1e-6, "last={}", cr.last);
        // intersection_points exposes the segment endpoints.
        let pts = ee.intersection_points();
        assert_eq!(pts.len(), 2, "span endpoints: {pts:?}");
        assert!((pts[0].0 - 1.0).abs() < 1e-6, "t1={}", pts[0].0);
        assert!((pts[1].0 - 3.0).abs() < 1e-6, "t1={}", pts[1].0);
        assert!(pts[0].2.distance(&GpPnt::new(1.0, 0.0, 0.0)) < 1e-6);
        assert!(pts[1].2.distance(&GpPnt::new(3.0, 0.0, 0.0)) < 1e-6);
        clear_tree(&e1.0);
        clear_tree(&e2.0);
    }

    #[test]
    fn line_line_full_crossing_single_point() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 2.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(0.0, 2.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let mut ee = prepared_ee(&e1, &e2, 0.0);
        let found = ee.compute_line_line_full().unwrap();
        assert!(found);
        assert_eq!(ee.points().len(), 1, "points: {:?}", ee.points());
        let p = &ee.points()[0];
        assert!(p.pnt1.distance(&GpPnt::new(1.0, 1.0, 0.0)) < 1e-6, "pnt={:?}", p.pnt1);
        let d = 2.0f64.sqrt();
        assert!((p.uv1.0 - d).abs() < 1e-6, "u1={}", p.uv1.0);
        assert!((p.uv2.0 - d).abs() < 1e-6, "u2={}", p.uv2.0);
        clear_tree(&e1.0);
        clear_tree(&e2.0);
    }

    #[test]
    fn line_line_full_parallel_distinct_empty() {
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(0.0, 1.0, 0.0), &GpPnt::new(2.0, 1.0, 0.0));
        let mut ee = prepared_ee(&e1, &e2, 0.0);
        let found = ee.compute_line_line_full().unwrap();
        assert!(!found);
        assert!(ee.points().is_empty());
        assert!(ee.common_parts().is_empty());
        assert!(ee.coincident_range().is_none());
        clear_tree(&e1.0);
        clear_tree(&e2.0);
    }

    #[test]
    fn circle_circle_full_tangent_single_point() {
        let b = TopoBuilder::new();
        // r = 1 at the origin and at (2,0,0): externally tangent at (1,0,0).
        let c1 = b.make_edge(Arc::new(GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0))), 0.0, 2.0 * PI);
        let ax2 = GpAx2::new(GpPnt::new(2.0, 0.0, 0.0), dir(0.0, 0.0, 1.0), dir(1.0, 0.0, 0.0)).unwrap();
        let c2 = b.make_edge(Arc::new(GeomCircle::new(GpCirc::new(ax2, 1.0))), 0.0, 2.0 * PI);
        let mut ee = prepared_ee(&c1, &c2, 1e-7);
        let found = ee.compute_circle_circle_full().unwrap();
        assert!(found, "tangent circles touch");
        let pts = ee.points();
        assert_eq!(pts.len(), 1, "points: {pts:?}");
        assert!(pts[0].pnt1.distance(&GpPnt::new(1.0, 0.0, 0.0)) < 1e-6, "pnt={:?}", pts[0].pnt1);
        clear_tree(&c1.0);
        clear_tree(&c2.0);
    }

    #[test]
    fn circle_circle_full_intersecting_two_points() {
        let b = TopoBuilder::new();
        let c1 = b.make_edge(Arc::new(GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0))), 0.0, 2.0 * PI);
        let ax2 = GpAx2::new(GpPnt::new(1.0, 0.0, 0.0), dir(0.0, 0.0, 1.0), dir(1.0, 0.0, 0.0)).unwrap();
        let c2 = b.make_edge(Arc::new(GeomCircle::new(GpCirc::new(ax2, 1.0))), 0.0, 2.0 * PI);
        let mut ee = prepared_ee(&c1, &c2, 1e-7);
        let found = ee.compute_circle_circle_full().unwrap();
        assert!(found);
        let pts = ee.points();
        assert_eq!(pts.len(), 2, "points: {pts:?}");
        for p in pts {
            assert!((p.pnt1.x() - 0.5).abs() < 1e-6, "x={}", p.pnt1.x());
            assert!((p.pnt1.y().abs() - 0.75f64.sqrt()).abs() < 1e-4, "y={}", p.pnt1.y());
            // Both parameter pairs must fall inside the full circle range.
            assert!(p.uv1.0 >= 0.0 && p.uv1.0 <= 2.0 * PI, "u1={}", p.uv1.0);
            assert!(p.uv2.0 >= 0.0 && p.uv2.0 <= 2.0 * PI, "u2={}", p.uv2.0);
        }
        // The two hits are distinct in parameter space.
        assert!((pts[0].uv1.0 - pts[1].uv1.0).abs() > 1e-3, "u1 pair");
        assert!((pts[0].uv2.0 - pts[1].uv2.0).abs() > 1e-3, "u2 pair");
        // intersection_points mirrors the discrete hits.
        assert_eq!(ee.intersection_points().len(), 2);
        clear_tree(&c1.0);
        clear_tree(&c2.0);
    }

    #[test]
    fn circle_circle_full_separated_empty() {
        let b = TopoBuilder::new();
        let c1 = b.make_edge(Arc::new(GeomCircle::new(GpCirc::new(GpAx2::standard(), 1.0))), 0.0, 2.0 * PI);
        let ax2 = GpAx2::new(GpPnt::new(5.0, 0.0, 0.0), dir(0.0, 0.0, 1.0), dir(1.0, 0.0, 0.0)).unwrap();
        let c2 = b.make_edge(Arc::new(GeomCircle::new(GpCirc::new(ax2, 1.0))), 0.0, 2.0 * PI);
        let mut ee = prepared_ee(&c1, &c2, 1e-7);
        let found = ee.compute_circle_circle_full().unwrap();
        assert!(!found, "separated circles do not intersect");
        assert!(ee.points().is_empty());
        assert!(ee.common_parts().is_empty());
        clear_tree(&c1.0);
        clear_tree(&c2.0);
    }

    #[test]
    fn find_best_solution_selects_minimum_residual() {
        let b = TopoBuilder::new();
        // X-line crossing a Y-line at (1, 0, 0): c2 runs through (1,0,0) at
        // its parameter 1 (segment (1,-1,0)→(1,1,0)).
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(1.0, -1.0, 0.0), &GpPnt::new(1.0, 1.0, 0.0));
        let ee = prepared_ee(&e1, &e2, 0.0);

        let exact = (1.0, 1.0, GpPnt::new(1.0, 0.0, 0.0));
        let far = (0.0, 0.0, GpPnt::new(0.0, 0.0, 0.0)); // residual ≈ √2 → dropped
        let best = ee.find_best_solution(&[far, exact], 1e-6).unwrap();
        assert_eq!(best.0, 1.0);
        assert_eq!(best.1, 1.0);
        assert!(best.2.distance(&GpPnt::new(1.0, 0.0, 0.0)) < 1e-12);

        // All-degenerate (or empty) candidate lists yield None.
        assert!(ee.find_best_solution(&[far, far], 1e-6).is_none());
        assert!(ee.find_best_solution(&[], 1e-6).is_none());
        clear_tree(&e1.0);
        clear_tree(&e2.0);
    }
}
