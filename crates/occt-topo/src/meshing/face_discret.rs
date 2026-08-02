//! Port of OCCT face_discret — Wave 1 BRepMesh.
//!
//! `FaceDiscret` turns a face's boundary (already discretized edges) plus an
//! interior deflection grid into a UV vertex set (`BRepMesh_FaceDiscret`),
//! `FaceChecker` decides whether a face is degenerate / self-intersecting and
//! therefore meshable (`BRepMesh_FaceChecker`), and `Classifier` answers
//! "is this UV point inside the face?" with a winding-number test
//! (`BRepMesh_Classifier` + `CSLib_Class2d`).
//!
//! Contract types (`MeshParameters`/`MeshFace`/`GeomTool`) are owned by the
//! sibling stubs (`parameters.rs`/`data_model.rs`/`geom_tool.rs`). Until those
//! stubs are filled this module uses the stand-ins from `edge_discret` and
//! carries a local `GeomTool` stand-in so the algorithms stay self-testable.

use occt_core::gp::GpPnt2d;

use super::edge_discret::MeshFace;
use super::parameters::MeshParameters;

// ---------------------------------------------------------------------------
// GeomTool (stand-in) — 2D segment-segment intersection
// ---------------------------------------------------------------------------
// ponytail: temporary stand-in for `super::geom_tool::GeomTool`; the
// `IntSegSeg` port below is the only piece `FaceChecker` needs.

/// State of a 2D segment intersection check (`BRepMesh_GeomTool::IntFlag`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntFlag {
    NoIntersection,
    Cross,
    EndPointTouch,
    PointOnSegment,
    Glued,
    Same,
}

/// 2D geometry helpers (`BRepMesh_GeomTool`).
pub struct GeomTool;

impl GeomTool {
    /// Intersect segments `(s1, e1)` and `(s2, e2)`. When `consider_end_touch`
    /// / `consider_point_on_segment` are false the touching cases degrade to
    /// `NoIntersection` (matching `BRepMesh_FaceChecker`'s call). The
    /// intersection point, when any, is written to `int_pnt`.
    pub fn int_seg_seg(
        s1: &GpPnt2d,
        e1: &GpPnt2d,
        s2: &GpPnt2d,
        e2: &GpPnt2d,
        consider_end_touch: bool,
        consider_point_on_segment: bool,
        int_pnt: &mut GpPnt2d,
    ) -> IntFlag {
        let v = e1.coord.subtracted(&s1.coord);
        let w = e2.coord.subtracted(&s2.coord);
        let denom = v.crossed(&w);

        if denom.abs() < 1e-14 {
            // Parallel: collinear overlap -> Same/Glued, else nothing.
            let r = s2.coord.subtracted(&s1.coord);
            let collinear = r.crossed(&v).abs() < 1e-12;
            if collinear && segments_overlap(s1, e1, s2, e2) {
                *int_pnt = *s2;
                return IntFlag::Same;
            }
            return IntFlag::NoIntersection;
        }

        // t along (s1->e1), s along (s2->e2).
        let r = s2.coord.subtracted(&s1.coord);
        let t = r.crossed(&w) / denom;
        let s = r.crossed(&v) / denom;

        if t < -1e-12 || t > 1.0 + 1e-12 || s < -1e-12 || s > 1.0 + 1e-12 {
            return IntFlag::NoIntersection;
        }

        let p = GpPnt2d::new(
            s1.x() + t * v.x,
            s1.y() + t * v.y,
        );
        *int_pnt = p;

        let at_end = t <= 1e-12 || (t - 1.0).abs() <= 1e-12;
        let bs_end = s <= 1e-12 || (s - 1.0).abs() <= 1e-12;
        if (at_end || bs_end) && consider_end_touch {
            return IntFlag::EndPointTouch;
        }
        if (at_end || bs_end) && !consider_end_touch {
            return IntFlag::NoIntersection;
        }
        if consider_point_on_segment {
            return IntFlag::PointOnSegment;
        }
        IntFlag::Cross
    }
}

/// Whether the projection intervals of two (parallel, collinear) segments overlap.
fn segments_overlap(a: &GpPnt2d, b: &GpPnt2d, c: &GpPnt2d, d: &GpPnt2d) -> bool {
    let (p0, p1) = if a.x() <= b.x() { (a.x(), b.x()) } else { (b.x(), a.x()) };
    let (q0, q1) = if c.x() <= d.x() { (c.x(), d.x()) } else { (d.x(), c.x()) };
    let horiz_overlap = p0 <= q1 + 1e-12 && q0 <= p1 + 1e-12;
    let (p0y, p1y) = if a.y() <= b.y() { (a.y(), b.y()) } else { (b.y(), a.y()) };
    let (q0y, q1y) = if c.y() <= d.y() { (c.y(), d.y()) } else { (d.y(), c.y()) };
    horiz_overlap && p0y <= q1y + 1e-12 && q0y <= p1y + 1e-12
}

// ---------------------------------------------------------------------------
// Classifier
// ---------------------------------------------------------------------------

/// Result of a 2D point-in-face test.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PointState {
    In,
    Out,
    On,
}

/// Classifies UV points against the wires of a discrete face.
///
/// Port of `BRepMesh_Classifier`: wires are registered as closed UV polygons and
/// [`Classifier::perform`] decides whether a point lies inside the face. The
/// first registered wire is the outer boundary; every subsequent wire is a
/// hole. A point on any boundary is reported as [`PointState::On`] (treated as
/// OUT by [`Classifier::is_inside`], matching OCCT).
#[derive(Debug, Default)]
pub struct Classifier {
    wires: Vec<Vec<GpPnt2d>>,
    tolerances: Vec<f64>,
}

impl Classifier {
    pub fn new() -> Self {
        Self::default()
    }

    /// Register a boundary wire. `tol_uv` is the parametric tolerance used for
    /// the "on boundary" decision; `_range_u` / `_range_v` are accepted for
    /// OCCT signature parity but not needed by the winding-number test.
    pub fn register_wire(
        &mut self,
        wire: &[GpPnt2d],
        tol_uv: (f64, f64),
        _range_u: (f64, f64),
        _range_v: (f64, f64),
    ) {
        if wire.len() < 2 {
            return;
        }
        self.wires.push(wire.to_vec());
        self.tolerances.push(tol_uv.0.max(tol_uv.1));
    }

    /// Number of registered wires.
    pub fn wires_nb(&self) -> usize {
        self.wires.len()
    }

    /// Classify a point against the registered wires. First wire is outer, the
    /// rest are holes. A point on any boundary is `On`.
    pub fn perform(&self, p: &GpPnt2d) -> PointState {
        let Some(outer) = self.wires.first() else {
            return PointState::Out;
        };
        let tol = self.tolerances.first().copied().unwrap_or(1e-7);
        match classify_polygon(p, outer, tol) {
            PointState::On => return PointState::On,
            PointState::Out => return PointState::Out,
            PointState::In => {}
        }
        for (i, hole) in self.wires.iter().enumerate().skip(1) {
            let tol = self.tolerances.get(i).copied().unwrap_or(1e-7);
            match classify_polygon(p, hole, tol) {
                PointState::On => return PointState::On,
                PointState::In => return PointState::Out,
                PointState::Out => {}
            }
        }
        PointState::In
    }

    /// True when the point lies strictly inside the face.
    pub fn is_inside(&self, p: &GpPnt2d) -> bool {
        self.perform(p) == PointState::In
    }
}

/// Classify a single polygon with the winding-number / ray-crossing rule.
/// Returns `On` when `p` lies within `tol` of any segment.
pub fn classify_polygon(p: &GpPnt2d, polygon: &[GpPnt2d], tol: f64) -> PointState {
    let n = polygon.len();
    if n < 3 {
        return PointState::Out;
    }
    for i in 0..n {
        let a = polygon[i];
        let b = polygon[(i + 1) % n];
        if point_segment_dist_2d(p, &a, &b) <= tol {
            return PointState::On;
        }
    }
    if winding_number(p, polygon) != 0 {
        PointState::In
    } else {
        PointState::Out
    }
}

/// Winding number of `p` around `polygon` (nonzero ⇒ inside).
fn winding_number(p: &GpPnt2d, polygon: &[GpPnt2d]) -> i32 {
    let n = polygon.len();
    let mut wn = 0i32;
    for i in 0..n {
        let a = polygon[i];
        let b = polygon[(i + 1) % n];
        if a.y() <= p.y() {
            if b.y() > p.y() && is_left(&a, &b, p) > 0.0 {
                wn += 1;
            }
        } else if b.y() <= p.y() && is_left(&a, &b, p) < 0.0 {
            wn -= 1;
        }
    }
    wn
}

/// Sign of the cross product `(b-a) × (p-a)`.
fn is_left(a: &GpPnt2d, b: &GpPnt2d, p: &GpPnt2d) -> f64 {
    (b.x() - a.x()) * (p.y() - a.y()) - (p.x() - a.x()) * (b.y() - a.y())
}

// ---------------------------------------------------------------------------
// FaceChecker
// ---------------------------------------------------------------------------

/// Checks whether a face is degenerate or has self-intersecting boundary wires.
///
/// Port of `BRepMesh_FaceChecker`: every boundary edge of every wire becomes a
/// 2D segment; pairs of segments (within a wire and across wires) are tested
/// for proper crossings via [`GeomTool::int_seg_seg`]. Degenerate faces (no
/// outer wire, fewer than three distinct points, zero UV area) are rejected by
/// [`FaceChecker::check`].
#[derive(Debug)]
pub struct FaceChecker<'a> {
    face: &'a MeshFace,
    /// Collected `(p1, p2, wire_index)` boundary segments.
    segments: Vec<(GpPnt2d, GpPnt2d, usize)>,
    /// Registered proper crossings: `(segment_a, segment_b, intersection)`.
    intersections: Vec<(usize, usize, GpPnt2d)>,
}

impl<'a> FaceChecker<'a> {
    /// `tolerance` is accepted for OCCT API parity
    /// (`BRepMesh_FaceChecker(face, parameters)`) but the O(n²) segment scan
    /// does not need it.
    pub fn new(face: &'a MeshFace, tolerance: f64) -> Self {
        let _ = tolerance;
        let mut segments = Vec::new();
        for (wi, wire) in std::iter::once(&face.outer_wire)
            .chain(face.inner_wires.iter())
            .enumerate()
        {
            // Closed polygon: every edge including the closing `last -> first`.
            let n = wire.len();
            for i in 0..n {
                segments.push((wire[i], wire[(i + 1) % n], wi));
            }
        }
        Self {
            face,
            segments,
            intersections: Vec::new(),
        }
    }

    /// Run the self-intersection check. Returns true when no boundary segment
    /// properly crosses another (the face is wire-consistent).
    pub fn perform(&mut self) -> bool {
        self.intersections.clear();
        let n = self.segments.len();
        for i in 0..n {
            for j in (i + 1)..n {
                let (a1, a2, wa) = self.segments[i];
                let (b1, b2, wb) = self.segments[j];
                // Consecutive segments inside the same wire share an endpoint;
                // that is a legitimate adjacency, not a crossing.
                if wa == wb && (j == i + 1 || (i == 0 && j == n - 1)) {
                    continue;
                }
                let mut ip = GpPnt2d::zero();
                let flag = GeomTool::int_seg_seg(
                    &a1,
                    &a2,
                    &b1,
                    &b2,
                    false,
                    false,
                    &mut ip,
                );
                if flag == IntFlag::Cross || flag == IntFlag::Same {
                    self.intersections.push((i, j, ip));
                }
            }
        }
        self.intersections.is_empty()
    }

    /// Intersecting segment pairs from the last [`Self::perform`].
    pub fn intersecting_segments(&self) -> &[(usize, usize, GpPnt2d)] {
        &self.intersections
    }

    /// Whether the face can be meshed: non-degenerate and (when checked) free
    /// of boundary self-intersections.
    pub fn is_meshable(&mut self) -> bool {
        !Self::is_degenerate(self.face) && self.perform()
    }

    /// Whether a face is degenerate: no outer wire, fewer than three distinct
    /// outer points, or zero UV area.
    pub fn is_degenerate(face: &MeshFace) -> bool {
        if face.outer_wire.len() < 3 {
            return true;
        }
        // Distinct points only.
        let mut distinct = 0usize;
        let mut prev: Option<GpPnt2d> = None;
        for p in &face.outer_wire {
            if prev.map_or(true, |q| q.distance(p) > 1e-9) {
                distinct += 1;
                prev = Some(*p);
            }
        }
        if distinct < 3 {
            return true;
        }
        // Zero (or degenerate) signed area.
        Self::outer_signed_area2(face).abs() < 1e-12
    }

    /// Twice the signed area of the outer wire polygon (shoelace formula).
    pub fn outer_signed_area2(face: &MeshFace) -> f64 {
        let n = face.outer_wire.len();
        let mut a2 = 0.0;
        for i in 0..n {
            let p = face.outer_wire[i];
            let q = face.outer_wire[(i + 1) % n];
            a2 += p.x() * q.y() - q.x() * p.y();
        }
        a2
    }

    /// Human-readable check: `Ok(())` when the face is meshable, otherwise an
    /// `Err` describing the failure.
    pub fn check(face: &MeshFace) -> Result<(), String> {
        if face.outer_wire.is_empty() {
            return Err("face has no outer wire".into());
        }
        if Self::is_degenerate(face) {
            return Err("face is degenerate (fewer than 3 distinct points or zero area)".into());
        }
        let mut checker = FaceChecker::new(face, 1e-7);
        if !checker.perform() {
            let n = checker.intersections.len();
            return Err(format!("face boundary has {n} self-intersection(s)"));
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// FaceDiscret
// ---------------------------------------------------------------------------

/// Discretizes a face into a set of UV vertices.
///
/// Port of `BRepMesh_FaceDiscret`: the boundary UV points (from the already
/// discretized edges) are kept verbatim; the interior is covered by a uniform
/// UV grid whose resolution follows the deflection, filtered by the
/// [`Classifier`] so only points strictly inside the face survive.
#[derive(Debug, Clone)]
pub struct FaceDiscret {
    params: MeshParameters,
}

impl FaceDiscret {
    pub fn new(params: MeshParameters) -> Self {
        Self { params }
    }

    /// The parameters in effect.
    pub fn parameters(&self) -> &MeshParameters {
        &self.params
    }

    /// Discretize a face into UV vertices: boundary plus interior grid.
    pub fn discretize_face(&self, face: &MeshFace) -> Vec<GpPnt2d> {
        let mut pts = Self::boundary_points(face);
        pts.extend(Self::interior_grid(face, self.params.deflection, 4096));
        dedup_points(pts)
    }

    /// UV points of every boundary wire, in order (outer then holes).
    pub fn boundary_points(face: &MeshFace) -> Vec<GpPnt2d> {
        let mut out = face.outer_wire.clone();
        for hole in &face.inner_wires {
            out.extend_from_slice(hole);
        }
        out
    }

    /// Interior UV grid across the outer wire's bounding box. Grid resolution
    /// is `deflection` per cell; the total point budget is capped by
    /// `max_points`. Points outside the face (or inside a hole) are dropped.
    pub fn interior_grid(
        face: &MeshFace,
        deflection: f64,
        max_points: usize,
    ) -> Vec<GpPnt2d> {
        let (umin, umax, vmin, vmax) = Self::uv_bounds(&face.outer_wire);
        let (du, dv) = (umax - umin, vmax - vmin);
        if du <= 0.0 || dv <= 0.0 {
            return Vec::new();
        }
        let def = deflection.max(1e-9);
        let mut nu = (du / def).ceil() as usize + 1;
        let mut nv = (dv / def).ceil() as usize + 1;
        // Enforce the point budget: scale the finer axis down if needed.
        let budget = max_points.max(4) as f64;
        let scale = (budget / (nu as f64 * nv as f64)).sqrt().min(1.0);
        nu = (nu as f64 * scale).ceil().max(2.0) as usize;
        nv = (nv as f64 * scale).ceil().max(2.0) as usize;

        let mut classifier = Classifier::new();
        classifier.register_wire(&face.outer_wire, (1e-9, 1e-9), (umin, umax), (vmin, vmax));
        for hole in &face.inner_wires {
            classifier.register_wire(hole, (1e-9, 1e-9), (umin, umax), (vmin, vmax));
        }

        let mut out = Vec::with_capacity(nu * nv);
        for j in 0..nv {
            for i in 0..nu {
                let u = umin + du * i as f64 / (nu - 1) as f64;
                let v = vmin + dv * j as f64 / (nv - 1) as f64;
                let p = GpPnt2d::new(u, v);
                if classifier.is_inside(&p) {
                    out.push(p);
                }
            }
        }
        out
    }

    /// Axis-aligned UV bounds of a closed boundary polygon.
    pub fn uv_bounds(wire: &[GpPnt2d]) -> (f64, f64, f64, f64) {
        let mut umin = f64::INFINITY;
        let mut umax = f64::NEG_INFINITY;
        let mut vmin = f64::INFINITY;
        let mut vmax = f64::NEG_INFINITY;
        for p in wire {
            umin = umin.min(p.x());
            umax = umax.max(p.x());
            vmin = vmin.min(p.y());
            vmax = vmax.max(p.y());
        }
        if !(umin.is_finite() && umax.is_finite() && vmin.is_finite() && vmax.is_finite()) {
            return (0.0, 1.0, 0.0, 1.0);
        }
        (umin, umax, vmin, vmax)
    }
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Distance from `p` to the segment `a..b`.
fn point_segment_dist_2d(p: &GpPnt2d, a: &GpPnt2d, b: &GpPnt2d) -> f64 {
    let ab = b.coord.subtracted(&a.coord);
    let len2 = ab.square_modulus();
    if len2 <= f64::EPSILON {
        return p.distance(a);
    }
    let ap = p.coord.subtracted(&a.coord);
    let t = (ap.dot(&ab) / len2).clamp(0.0, 1.0);
    let proj = a.coord.added(&ab.multiplied(t));
    p.coord.subtracted(&proj).modulus()
}

/// Remove consecutive duplicates (within `1e-9`) from a UV point list.
fn dedup_points(pts: Vec<GpPnt2d>) -> Vec<GpPnt2d> {
    let mut out: Vec<GpPnt2d> = Vec::with_capacity(pts.len());
    for p in pts {
        if out.last().map_or(true, |q| q.distance(&p) > 1e-9) {
            out.push(p);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use occt_core::gp::GpPnt2d;

    fn square() -> Vec<GpPnt2d> {
        vec![
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(1.0, 0.0),
            GpPnt2d::new(1.0, 1.0),
            GpPnt2d::new(0.0, 1.0),
        ]
    }

    #[test]
    fn planar_face_classification_inside_outside() {
        let poly = square();
        assert_eq!(classify_polygon(&GpPnt2d::new(0.5, 0.5), &poly, 1e-9), PointState::In);
        assert_eq!(classify_polygon(&GpPnt2d::new(1.5, 0.5), &poly, 1e-9), PointState::Out);
        assert_eq!(classify_polygon(&GpPnt2d::new(0.5, -0.1), &poly, 1e-9), PointState::Out);
        // On-boundary within tolerance -> On.
        assert_eq!(classify_polygon(&GpPnt2d::new(0.0, 0.5), &poly, 1e-9), PointState::On);
        assert_eq!(classify_polygon(&GpPnt2d::new(0.5, 1.0), &poly, 1e-9), PointState::On);
    }

    #[test]
    fn classifier_handles_hole() {
        let outer = square();
        let hole = vec![
            GpPnt2d::new(0.25, 0.25),
            GpPnt2d::new(0.75, 0.25),
            GpPnt2d::new(0.75, 0.75),
            GpPnt2d::new(0.25, 0.75),
        ];
        let mut c = Classifier::new();
        c.register_wire(&outer, (1e-9, 1e-9), (0.0, 1.0), (0.0, 1.0));
        c.register_wire(&hole, (1e-9, 1e-9), (0.0, 1.0), (0.0, 1.0));
        assert_eq!(c.perform(&GpPnt2d::new(0.5, 0.5)), PointState::Out, "point in hole");
        assert_eq!(c.perform(&GpPnt2d::new(0.1, 0.1)), PointState::In, "point in ring");
        assert_eq!(c.perform(&GpPnt2d::new(0.5, 1.1)), PointState::Out, "point outside outer");
    }

    #[test]
    fn degenerate_face_detected() {
        let mut face = MeshFace::new(square());
        assert!(!FaceChecker::is_degenerate(&face), "unit square is not degenerate");
        assert!(FaceChecker::check(&face).is_ok());

        // Collapsed outer wire -> degenerate.
        face.outer_wire = vec![
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(1.0, 0.0),
            GpPnt2d::new(0.0, 0.0),
        ];
        assert!(FaceChecker::is_degenerate(&face));
        assert!(FaceChecker::check(&face).is_err());

        // Empty outer wire.
        face.outer_wire = Vec::new();
        assert!(FaceChecker::is_degenerate(&face));
        assert!(FaceChecker::check(&face).is_err());
    }

    #[test]
    fn face_checker_finds_self_intersection() {
        // Bowtie: two triangles sharing a crossing diagonal. Its two lobes
        // cancel in the shoelace sum, so it is zero-area AND self-intersecting.
        let bowtie = vec![
            GpPnt2d::new(0.0, 0.0),
            GpPnt2d::new(1.0, 1.0),
            GpPnt2d::new(1.0, 0.0),
            GpPnt2d::new(0.0, 1.0),
        ];
        let face = MeshFace::new(bowtie);
        let mut checker = FaceChecker::new(&face, 1e-7);
        assert!(!checker.perform(), "bowtie must self-intersect");
        assert!(!checker.intersecting_segments().is_empty());
        assert!(FaceChecker::check(&face).is_err(), "self-intersecting face must not be meshable");

        // Clean square passes.
        let face = MeshFace::new(square());
        let mut checker = FaceChecker::new(&face, 1e-7);
        assert!(checker.perform());
        assert!(checker.intersecting_segments().is_empty());
        assert!(FaceChecker::check(&face).is_ok());
    }

    #[test]
    fn pentagram_self_intersects_but_is_not_degenerate() {
        // Pentagram traced in star order: nonzero area yet self-intersecting.
        let star = vec![
            GpPnt2d::new(0.0, -1.0),
            GpPnt2d::new(0.587_785_252_3, 0.809_016_994_4),
            GpPnt2d::new(-0.951_056_516_3, -0.309_016_994_4),
            GpPnt2d::new(0.951_056_516_3, -0.309_016_994_4),
            GpPnt2d::new(-0.587_785_252_3, 0.809_016_994_4),
        ];
        let face = MeshFace::new(star);
        assert!(!FaceChecker::is_degenerate(&face), "star has nonzero area");
        let mut checker = FaceChecker::new(&face, 1e-7);
        assert!(!checker.perform(), "pentagram must self-intersect");
        assert!(FaceChecker::check(&face).is_err());
    }

    #[test]
    fn int_seg_seg_cross_vs_touch() {
        let mut ip = GpPnt2d::zero();
        let s1 = GpPnt2d::new(0.0, 0.0);
        let e1 = GpPnt2d::new(2.0, 0.0);
        let s2 = GpPnt2d::new(1.0, -1.0);
        let e2 = GpPnt2d::new(1.0, 1.0);
        assert_eq!(
            GeomTool::int_seg_seg(&s1, &e1, &s2, &e2, false, false, &mut ip),
            IntFlag::Cross
        );
        assert!(ip.distance(&GpPnt2d::new(1.0, 0.0)) < 1e-9);

        // Endpoint touch without the flag degrades to NoIntersection.
        let s3 = GpPnt2d::new(2.0, 1.0);
        let e3 = GpPnt2d::new(2.0, -1.0);
        assert_eq!(
            GeomTool::int_seg_seg(&s1, &e1, &s3, &e3, false, false, &mut ip),
            IntFlag::NoIntersection
        );
        // Parallel disjoint.
        let s4 = GpPnt2d::new(0.0, 1.0);
        let e4 = GpPnt2d::new(2.0, 1.0);
        assert_eq!(
            GeomTool::int_seg_seg(&s1, &e1, &s4, &e4, false, false, &mut ip),
            IntFlag::NoIntersection
        );
    }

    #[test]
    fn face_discret_produces_boundary_and_interior() {
        let face = MeshFace::new(square());
        let params = MeshParameters {
            deflection: 0.25,
            ..MeshParameters::default()
        };
        let fd = FaceDiscret::new(params);
        let pts = fd.discretize_face(&face);

        // Boundary corners are present.
        for c in square() {
            assert!(
                pts.iter().any(|p| p.distance(&c) < 1e-9),
                "missing corner {c:?}"
            );
        }
        // Interior points exist and are inside the unit square.
        let interior = pts
            .iter()
            .filter(|p| p.x() > 0.05 && p.x() < 0.95 && p.y() > 0.05 && p.y() < 0.95)
            .count();
        assert!(interior > 0, "no interior points generated");
        // No point outside the square bounds.
        assert!(pts.iter().all(|p| p.x() >= -1e-9 && p.x() <= 1.0 + 1e-9 && p.y() >= -1e-9 && p.y() <= 1.0 + 1e-9));
    }
}
