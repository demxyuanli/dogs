//! 2D point-in-face classifier — a port of `IntTools_FClass2d` (TKBO).
//!
//! Classifies a 2D point `(u, v)` in the parameter domain of a face against
//! the face's UV boundary. The face boundary is sampled edge-by-edge into 2D
//! polygons (pcurves preferred, sampled projection as fallback), chained into
//! closed rings by UV continuity, and decomposed into one outer ring plus hole
//! rings (largest |area| ring = outer). A point is:
//!
//! * `On`  — within `tol` of a boundary ring segment;
//! * `In`  — inside the outer ring and outside every hole ring;
//! * `Out` — otherwise.
//!
//! Periodic surfaces (cylinder, cone, sphere, torus) fold out-of-range `u`/`v`
//! by whole periods (`AdjustPeriodic` semantics) before classifying, and try
//! the periodic images in order so `perform((u + 2π, v))` agrees with
//! `perform((u, v))`.
//!
//! Source: `IntTools_FClass2d.hxx/.cxx` (TKBO). The UV rings are built by
//! sampling each edge pcurve (via [`crate::pcurve_full::make_pcurve_full`],
//! falling back to [`crate::brep_surface::edge_pcurve_on_face`]); edge
//! orientation within a wire is recovered by UV-continuity chaining, since this
//! port's flat `TShape` child tree drops the per-edge `TopAbs_Orientation`.

use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::geom::polygon_ops::{point_in_polygon2d, polygon_area2d};
use occt_core::gp::GpPnt2d;
use occt_core::precision::SQUARE_CONFUSION;
use occt_geom::Surface;
use occt_geom2d::curve::Curve2d;

use crate::abs::Orientation;
use crate::brep_surface::{edge_pcurve_on_face, face_uv_bounds};
use crate::pcurve::{pc_curve_kind, CurveKind};
use crate::pcurve_full::make_pcurve_full;
use crate::shape::{Edge, Face};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{edges_of_wire, wires_of_face};

/// State of a 2D point relative to the face region. Source: `TopAbs_State`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FaceState {
    In,
    Out,
    On,
    Unknown,
}

/// A face boundary region decomposed into an outer ring and hole rings.
#[derive(Debug, Clone, Default)]
pub struct FaceRegion {
    pub outer: Vec<GpPnt2d>,
    pub holes: Vec<Vec<GpPnt2d>>,
}

impl FaceRegion {
    /// Classify a 2D point against the region.
    ///
    /// `On` when the point is within `on_tol` of a boundary ring segment; `In`
    /// when it is inside the outer ring and outside every hole; `Out` otherwise.
    /// An empty region (no boundary) classifies every point `In`.
    fn point_in_region(&self, p: &GpPnt2d, on_tol: f64) -> FaceState {
        if self.outer.is_empty() {
            return FaceState::In;
        }
        if point_near_ring(&self.outer, p, on_tol) {
            return FaceState::On;
        }
        for hole in &self.holes {
            if point_near_ring(hole, p, on_tol) {
                return FaceState::On;
            }
        }
        if !point_in_polygon2d(&self.outer, p) {
            return FaceState::Out;
        }
        for hole in &self.holes {
            if point_in_polygon2d(hole, p) {
                return FaceState::Out;
            }
        }
        FaceState::In
    }
}

/// 2D point-in-face classifier. Source: `IntTools_FClass2d`.
#[derive(Debug, Clone)]
pub struct FClass2d {
    /// UV tolerance used for the `On` boundary test.
    tol: f64,
    /// The face whose boundary was sampled.
    face: Face,
    /// The sampled boundary region (None when the face has no boundary wires).
    region: Option<FaceRegion>,
    /// Whether the face surface is periodic in `u` / `v`.
    is_u_periodic: bool,
    is_v_periodic: bool,
    /// Surface periods (0.0 for non-periodic directions).
    u_period: f64,
    v_period: f64,
    /// Bounding box of the sampled boundary rings (`IntTools_FClass2d::Umin`…).
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
    /// Whether the face is a "hole" (its largest boundary loop winds clockwise
    /// in the face's UV). Source: `IntTools_FClass2d::IsHole`.
    my_is_hole: bool,
}

impl FClass2d {
    /// Build a classifier for `face` using the UV tolerance `tol`
    /// (`IntTools_FClass2d(F, TolUV)`).
    pub fn new(face: &Face, tol: f64) -> Result<Self, String> {
        let mut c = FClass2d {
            tol,
            face: Face::new(),
            region: None,
            is_u_periodic: false,
            is_v_periodic: false,
            u_period: 0.0,
            v_period: 0.0,
            umin: f64::INFINITY,
            umax: f64::NEG_INFINITY,
            vmin: f64::INFINITY,
            vmax: f64::NEG_INFINITY,
            my_is_hole: true,
        };
        c.init(face, tol)?;
        Ok(c)
    }

    /// (Re)initialize the classifier from `face` and tolerance `tol`
    /// (`IntTools_FClass2d::Init`).
    ///
    /// Samples every boundary wire into UV rings. Bad wires (open, degenerate,
    /// or with a zero-area ring) are skipped rather than failing the call, so a
    /// face with several wires still initializes when one is unusable.
    pub fn init(&mut self, face: &Face, tol: f64) -> Result<(), String> {
        self.tol = tol;
        self.face = face.clone();
        let Some(surf) = GeometryRegistry::global().face_surface(&face.0) else {
            return Err("FClass2d::init: face has no registered surface".into());
        };
        self.is_u_periodic = surf.is_u_periodic();
        self.is_v_periodic = surf.is_v_periodic();

        // Surface period from the natural parametric range; 2π fallback for an
        // unbounded periodic direction (cylinder/cone/sphere/torus use 2π).
        let (su0, su1, sv0, sv1) = face_uv_bounds(face);
        self.u_period = if self.is_u_periodic {
            if su0.is_finite() && su1.is_finite() {
                su1 - su0
            } else {
                2.0 * PI
            }
        } else {
            0.0
        };
        self.v_period = if self.is_v_periodic {
            if sv0.is_finite() && sv1.is_finite() {
                sv1 - sv0
            } else {
                2.0 * PI
            }
        } else {
            0.0
        };

        let mut rings: Vec<(Vec<GpPnt2d>, f64)> = Vec::new();
        let mut umin = f64::INFINITY;
        let mut umax = f64::NEG_INFINITY;
        let mut vmin = f64::INFINITY;
        let mut vmax = f64::NEG_INFINITY;

        for wire in wires_of_face(face) {
            let edges = edges_of_wire(&wire);
            let mut polylines: Vec<Vec<GpPnt2d>> = Vec::new();
            for edge in &edges {
                let or = edge.orientation();
                if or != Orientation::Forward && or != Orientation::Reversed {
                    continue;
                }
                let pl = edge_points(edge, face);
                if pl.len() < 2 {
                    continue;
                }
                polylines.push(pl);
            }
            if polylines.is_empty() {
                continue;
            }
            let ring = match chain_ring(&polylines, self.u_period, self.v_period) {
                Some(r) => dedup_ring(&r),
                None => continue,
            };
            if ring.len() < 3 {
                continue;
            }
            for p in &ring {
                umin = umin.min(p.x());
                umax = umax.max(p.x());
                vmin = vmin.min(p.y());
                vmax = vmax.max(p.y());
            }
            let area = polygon_area2d(&ring);
            if area.abs() < SQUARE_CONFUSION {
                continue;
            }
            rings.push((ring, area));
        }

        self.umin = umin;
        self.umax = umax;
        self.vmin = vmin;
        self.vmax = vmax;

        // Outer ring: the one with the largest |area|; every other ring is a
        // hole. This is winding-independent, so it stays correct for faces
        // whose UV loops happen to wind clockwise (e.g. inward-normal planes).
        let mut region = FaceRegion::default();
        if let Some((outer_idx, _)) = rings
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.1.abs().total_cmp(&b.1.abs()))
        {
            region.outer = rings[outer_idx].0.clone();
            for (i, (ring, _)) in rings.iter().enumerate() {
                if i != outer_idx {
                    region.holes.push(ring.clone());
                }
            }
        }
        self.region = if region.outer.is_empty() { None } else { Some(region) };

        // A face is a "hole" when its largest boundary loop winds clockwise
        // (material on the outside of the loop). An empty face counts as a hole.
        self.my_is_hole = match rings
            .iter()
            .max_by(|a, b| a.1.abs().total_cmp(&b.1.abs()))
        {
            Some((_, area)) => *area <= 0.0,
            None => true,
        };
        Ok(())
    }

    /// Classify the 2D point `puv` (`IntTools_FClass2d::Perform`).
    ///
    /// Periodic surface coordinates are folded by whole periods before
    /// classifying, and the periodic images are tried in order.
    pub fn perform(&self, puv: GpPnt2d) -> FaceState {
        self.perform_recadre(puv, true)
    }

    /// `Perform` with explicit `RecadreOnPeriodic` control.
    pub fn perform_recadre(&self, puv: GpPnt2d, recadre_on_periodic: bool) -> FaceState {
        self.perform_internal(puv, self.tol, recadre_on_periodic)
    }

    /// State of the infinite (far bottom-left) point: `In` when the face has
    /// no boundary (a closed periodic face such as a full sphere/torus),
    /// otherwise the classification of the UV-domain-outside corner point
    /// (`IntTools_FClass2d::PerformInfinitePoint`).
    pub fn perform_infinite_point(&self) -> FaceState {
        if !self.umin.is_finite()
            || !self.umax.is_finite()
            || !self.vmin.is_finite()
            || !self.vmax.is_finite()
        {
            return FaceState::In;
        }
        let p = GpPnt2d::new(
            self.umin - (self.umax - self.umin),
            self.vmin - (self.vmax - self.vmin),
        );
        self.perform_recadre(p, false)
    }

    /// Test whether `puv` lies on the face boundary restriction within `tol`
    /// (`IntTools_FClass2d::TestOnRestriction`). `On` when within `tol` of a
    /// boundary ring, `In`/`Out` otherwise.
    pub fn test_on_restriction(&self, puv: GpPnt2d, tol: f64) -> FaceState {
        self.perform_internal(puv, tol, true)
    }

    /// Whether the face is a "hole" (`IntTools_FClass2d::IsHole`).
    pub fn is_hole(&self) -> bool {
        self.my_is_hole
    }

    /// Classify `puv` against the sampled boundary region using the UV
    /// tolerance `on_tol`.
    pub fn point_in_region(&self, puv: GpPnt2d) -> FaceState {
        match &self.region {
            Some(r) => r.point_in_region(&puv, self.tol),
            None => FaceState::In,
        }
    }

    /// The sampled outer ring, if any (exposed for diagnostics/tests).
    pub fn outer_ring(&self) -> Option<&[GpPnt2d]> {
        self.region.as_ref().map(|r| r.outer.as_slice())
    }

    /// The sampled hole rings, if any.
    pub fn hole_rings(&self) -> &[Vec<GpPnt2d>] {
        self.region.as_ref().map(|r| r.holes.as_slice()).unwrap_or(&[])
    }

    /// `Perform` core: the periodic-image search loop from the C++ source.
    fn perform_internal(&self, puv: GpPnt2d, on_tol: f64, recadre: bool) -> FaceState {
        if self.region.is_none() {
            return FaceState::In;
        }
        let mut u = puv.x();
        let mut v = puv.y();
        let (mut uu, mut vv) = (u, v);
        if recadre {
            if self.is_u_periodic && self.u_period > 0.0 {
                uu = adjust_periodic(u, self.umin, self.umax, self.u_period).0;
            }
            if self.is_v_periodic && self.v_period > 0.0 {
                vv = adjust_periodic(v, self.vmin, self.vmax, self.v_period).0;
            }
        }
        let mut urecadre = false;
        let mut vrecadre = false;
        let mut a_status = FaceState::Unknown;
        for _ in 0..64 {
            a_status = self
                .region
                .as_ref()
                .map(|r| r.point_in_region(&GpPnt2d::new(u, v), on_tol))
                .unwrap_or(FaceState::In);

            if !recadre || (!self.is_u_periodic && !self.is_v_periodic) {
                return a_status;
            }
            if a_status == FaceState::In || a_status == FaceState::On {
                return a_status;
            }
            if !urecadre {
                u = uu;
                urecadre = true;
            } else if self.is_u_periodic {
                u += self.u_period;
            }
            if u > self.umax || !self.is_u_periodic {
                if !vrecadre {
                    v = vv;
                    vrecadre = true;
                } else if self.is_v_periodic {
                    v += self.v_period;
                }
                u = uu;
                if v > self.vmax || !self.is_v_periodic {
                    return a_status;
                }
            }
        }
        a_status
    }
}

// ---------------------------------------------------------------------------
// Edge sampling and ring chaining
// ---------------------------------------------------------------------------

/// Sample the UV pcurve of `edge` on `face` into a polyline (in the edge's
/// natural curve direction).
///
/// Preference: a pcurve stored on the edge for this face (matching
/// `BRep_Tool::CurveOnSurface`), then [`make_pcurve_full`], then a sampled
/// projection via [`edge_pcurve_on_face`]. Empty when the edge has no usable
/// curve or a degenerate parameter range.
fn edge_points(edge: &Edge, face: &Face) -> Vec<GpPnt2d> {
    let (a, b) = GeometryRegistry::global().edge_parameters(&edge.0);
    if !a.is_finite() || !b.is_finite() || b - a < 1e-15 {
        return Vec::new();
    }
    let face_key = GeometryRegistry::shape_key(&face.0);
    if let Some(pc) = GeometryRegistry::global().edge_pcurve(&edge.0, face_key) {
        return sample_pcurve(pc.as_ref(), a, b);
    }
    if let Ok(pc) = make_pcurve_full(edge, face) {
        return sample_pcurve(pc.as_ref(), a, b);
    }
    edge_pcurve_on_face(edge, face, 32)
}

/// Sample a pcurve uniformly over the edge range `[a, b]`.
fn sample_pcurve(pc: &dyn Curve2d, a: f64, b: f64) -> Vec<GpPnt2d> {
    let n = sample_count(pc);
    (0..n)
        .map(|i| {
            let t = a + (b - a) * i as f64 / (n.max(1) - 1) as f64;
            pc.d0(t)
        })
        .collect()
}

/// Number of samples for a pcurve: a straight UV line needs only its endpoints,
/// a circle pcurve 48 samples, anything else (curved, non-isoparametric) 32.
fn sample_count(pc: &dyn Curve2d) -> usize {
    match pc_curve_kind(pc) {
        CurveKind::Line => 2,
        CurveKind::Circle => 48,
        _ => 32,
    }
}

/// Chain edge polylines into a single closed ring by UV continuity.
///
/// The first polyline is taken in its natural direction; every following edge
/// is appended (forward or reversed) so its start connects to the current ring
/// end, shifting by whole periods when the surface is periodic. Among the
/// candidate connections the one with the fewest period shifts wins (ties go to
/// the natural direction). Returns `None` when the polylines do not form a
/// closed chain.
fn chain_ring(polylines: &[Vec<GpPnt2d>], u_per: f64, v_per: f64) -> Option<Vec<GpPnt2d>> {
    if polylines.is_empty() {
        return None;
    }
    let mut ring: Vec<GpPnt2d> = polylines[0].clone();
    let mut used = vec![false; polylines.len()];
    used[0] = true;

    for _ in 0..(polylines.len() * 2 + 16) {
        if used.iter().all(|&u| u) {
            break;
        }
        let cur_end = *ring.last().expect("ring is non-empty");
        let mut best: Option<(usize, bool, i64, i64)> = None;
        for i in 0..polylines.len() {
            if used[i] {
                continue;
            }
            let pl = &polylines[i];
            let fwd_start = *pl.first().expect("polyline is non-empty");
            let fwd_end = *pl.last().expect("polyline is non-empty");
            // Natural direction.
            if let Some((ku, kv)) = wrap_offset(&fwd_start, &cur_end, u_per, v_per) {
                if is_better(best, (i, false, ku, kv), u_per, v_per) {
                    best = Some((i, false, ku, kv));
                }
            }
            // Reversed direction.
            if let Some((ku, kv)) = wrap_offset(&fwd_end, &cur_end, u_per, v_per) {
                if is_better(best, (i, true, ku, kv), u_per, v_per) {
                    best = Some((i, true, ku, kv));
                }
            }
        }
        let (idx, rev, ku, kv) = best?;
        let pl = &polylines[idx];
        if rev {
            ring.extend(pl.iter().rev().map(|p| shift_pnt(p, ku, kv, u_per, v_per)).skip(1));
        } else {
            ring.extend(pl.iter().map(|p| shift_pnt(p, ku, kv, u_per, v_per)).skip(1));
        }
        used[idx] = true;
    }

    if !used.iter().all(|&u| u) {
        return None;
    }
    // Closure: the last point must coincide with the first (periodically).
    let start = *ring.first()?;
    let end = *ring.last()?;
    if wrap_offset(&start, &end, u_per, v_per).is_none() {
        return None;
    }
    Some(ring)
}

/// Whether candidate `(idx, rev, ku, kv)` beats `best` for the next chain link:
/// fewer period shifts wins; ties prefer the natural (non-reversed) direction.
fn is_better(
    best: Option<(usize, bool, i64, i64)>,
    cand: (usize, bool, i64, i64),
    _u_per: f64,
    _v_per: f64,
) -> bool {
    match best {
        None => true,
        Some((_, brev, bku, bkv)) => {
            let bcost = bku.unsigned_abs() + bkv.unsigned_abs();
            let cost = cand.2.unsigned_abs() + cand.3.unsigned_abs();
            cost < bcost || (cost == bcost && !cand.1 && brev)
        }
    }
}

/// The integer period offsets `(ku, kv)` such that `a + (ku*u_per, kv*v_per)`
/// equals `b` within tolerance; `None` when no such offsets exist.
fn wrap_offset(a: &GpPnt2d, b: &GpPnt2d, u_per: f64, v_per: f64) -> Option<(i64, i64)> {
    let du = b.x() - a.x();
    let dv = b.y() - a.y();
    let tol = 1e-6;
    let ku = if u_per > 0.0 && u_per.is_finite() {
        let k = (du / u_per).round();
        if (du - k * u_per).abs() <= tol * (1.0 + u_per.abs()) {
            k as i64
        } else {
            return None;
        }
    } else if du.abs() <= tol {
        0
    } else {
        return None;
    };
    let kv = if v_per > 0.0 && v_per.is_finite() {
        let k = (dv / v_per).round();
        if (dv - k * v_per).abs() <= tol * (1.0 + v_per.abs()) {
            k as i64
        } else {
            return None;
        }
    } else if dv.abs() <= tol {
        0
    } else {
        return None;
    };
    Some((ku, kv))
}

/// Translate a point by the periodic offset `(ku, kv)`.
fn shift_pnt(p: &GpPnt2d, ku: i64, kv: i64, u_per: f64, v_per: f64) -> GpPnt2d {
    let du = if u_per > 0.0 && u_per.is_finite() {
        ku as f64 * u_per
    } else {
        0.0
    };
    let dv = if v_per > 0.0 && v_per.is_finite() {
        kv as f64 * v_per
    } else {
        0.0
    };
    GpPnt2d::new(p.x() + du, p.y() + dv)
}

/// Remove consecutive duplicate points (degenerate zero-length segments).
fn dedup_ring(pts: &[GpPnt2d]) -> Vec<GpPnt2d> {
    let mut out: Vec<GpPnt2d> = Vec::with_capacity(pts.len());
    for &p in pts {
        if let Some(&last) = out.last() {
            if p.distance(&last) < 1e-9 {
                continue;
            }
        }
        out.push(p);
    }
    out
}

/// Fold `u` into `[umin, umax]` by adding/subtracting whole periods
/// (`GeomInt::AdjustPeriodic` semantics). Returns `(folded, offset)`.
fn adjust_periodic(u: f64, umin: f64, umax: f64, period: f64) -> (f64, f64) {
    if period <= 0.0
        || !period.is_finite()
        || !umin.is_finite()
        || !umax.is_finite()
        || umax <= umin
    {
        return (u, 0.0);
    }
    if u >= umin && u <= umax {
        return (u, 0.0);
    }
    let du = if u < umin {
        ((umin - u) / period).ceil() * period
    } else {
        -((u - umax) / period).ceil() * period
    };
    (u + du, du)
}

/// Minimum distance from `p` to any segment of a closed ring.
fn point_near_ring(ring: &[GpPnt2d], p: &GpPnt2d, tol: f64) -> bool {
    let n = ring.len();
    if n < 2 {
        return false;
    }
    for i in 0..n {
        let a = &ring[i];
        let b = &ring[(i + 1) % n];
        if point_segment_distance(p, a, b) <= tol {
            return true;
        }
    }
    false
}

/// Distance from a point to the segment `[a, b]`.
fn point_segment_distance(p: &GpPnt2d, a: &GpPnt2d, b: &GpPnt2d) -> f64 {
    let dx = b.x() - a.x();
    let dy = b.y() - a.y();
    let len2 = dx * dx + dy * dy;
    if len2 < 1e-24 {
        return p.distance(a);
    }
    let t = (((p.x() - a.x()) * dx + (p.y() - a.y()) * dy) / len2).clamp(0.0, 1.0);
    let qx = a.x() + t * dx;
    let qy = a.y() + t * dy;
    ((p.x() - qx).powi(2) + (p.y() - qy).powi(2)).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_extrema::test_box::unit_box;
    use crate::brep_surface::{face_is_planar, surface_closest_params};
    use crate::builder::TopoBuilder;
    use crate::builder_face::build_face_with_holes;
    use crate::primitives::BRepPrimCylinder;
    use crate::topo_tools_full::faces_of;
    use occt_core::gp::{GpAx2, GpAx3, GpCylinder, GpDir, GpPln, GpPnt, GpPnt2d};
    use occt_geom::{GeomCylinder, GeomPlane, Surface};

    fn p2(x: f64, y: f64) -> GpPnt2d {
        GpPnt2d::new(x, y)
    }

    /// UV parameters of a 3D point on a face surface (for tests).
    fn project_uv(face: &Face, p: &GpPnt) -> GpPnt2d {
        let surf = GeometryRegistry::global().face_surface(&face.0).expect("face surface");
        let (u, v) = surface_closest_params(surf.as_ref(), p, 32, 32);
        p2(u, v)
    }

    fn square_edges(b: &TopoBuilder, pts: &[GpPnt; 4]) -> Vec<Edge> {
        (0..4).map(|i| b.make_edge_segment(&pts[i], &pts[(i + 1) % 4])).collect()
    }

    #[test]
    fn unit_box_bottom_face_classifies_square() {
        let ub = unit_box();
        let face = &ub.faces[0]; // bottom (z = 0), UV square [0,1]²
        let cl = FClass2d::new(face, 1e-6).expect("classifier");
        assert_eq!(cl.perform(p2(0.5, 0.5)), FaceState::In);
        assert_eq!(cl.perform(p2(-0.5, 0.5)), FaceState::Out);
        assert_eq!(cl.perform(p2(0.5, 0.0)), FaceState::On);
        // The open bottom face has no closed periodic boundary: the infinite
        // point is Out.
        assert_eq!(cl.perform_infinite_point(), FaceState::Out);
    }

    #[test]
    fn face_with_square_hole() {
        let b = TopoBuilder::new();
        let outer = [
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
        ];
        let hole = [
            GpPnt::new(0.75, 0.75, 0.0),
            GpPnt::new(0.75, 0.25, 0.0),
            GpPnt::new(0.25, 0.25, 0.0),
            GpPnt::new(0.25, 0.75, 0.0),
        ];
        let face = build_face_with_holes(&square_edges(&b, &outer), &[square_edges(&b, &hole)])
            .expect("face with hole");
        let cl = FClass2d::new(&face, 1e-6).expect("classifier");

        // Hole centre -> Out; hole left-edge midpoint -> On; annulus -> In;
        // in the annulus near the outer corner -> In; outside the outer square
        // (projection clamps the grid to [-1,1], so use a direct UV point) -> Out.
        assert_eq!(cl.perform(project_uv(&face, &GpPnt::new(0.5, 0.5, 0.0))), FaceState::Out);
        assert_eq!(cl.perform(project_uv(&face, &GpPnt::new(0.25, 0.5, 0.0))), FaceState::On);
        assert_eq!(cl.perform(project_uv(&face, &GpPnt::new(0.9, 0.9, 0.0))), FaceState::In);
        assert_eq!(cl.perform(project_uv(&face, &GpPnt::new(0.05, 0.05, 0.0))), FaceState::In);
        assert_eq!(cl.perform(p2(1.5, 0.5)), FaceState::Out);
    }

    /// A clean cylinder lateral face: two cap circles plus two *distinct* seam
    /// edges whose stored pcurves sit at u = 0 (up) and u = 2π (down). The
    /// port's flat wire model cannot store the reversed seam orientation, so
    /// the pcurves are attached explicitly to make the UV ring a clean
    /// `[0, 2π] × [0, h]` rectangle.
    fn clean_cylinder_lateral(radius: f64, height: f64) -> Face {
        let b = TopoBuilder::new();
        let ax = GpAx3::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), &GpDir::new(1.0, 0.0, 0.0).unwrap())
            .expect("cylinder axis");
        let surface: Arc<dyn Surface> = Arc::new(GeomCylinder::new(
            GpCylinder::new(ax, radius).expect("cylinder radius"),
        ));

        let bottom = GpPnt::new(radius, 0.0, 0.0);
        let top = GpPnt::new(radius, 0.0, height);
        let ax2_bot = GpAx2::new(GpPnt::zero(), GpDir::new(0.0, 0.0, 1.0).unwrap(), GpDir::new(1.0, 0.0, 0.0).unwrap())
            .expect("bottom circle axis");
        let ax2_top = GpAx2::new(GpPnt::new(0.0, 0.0, height), GpDir::new(0.0, 0.0, 1.0).unwrap(), GpDir::new(1.0, 0.0, 0.0).unwrap())
            .expect("top circle axis");
        let bottom_circle = b.make_edge_circle(&ax2_bot, radius, 0.0, 2.0 * PI);
        let top_circle = b.make_edge_circle(&ax2_top, radius, 0.0, 2.0 * PI);
        let seam_up = b.make_edge_segment(&bottom, &top);
        let seam_down = b.make_edge_segment(&bottom, &top);

        let lateral = b.make_face(surface, &[]);
        let face_key = GeometryRegistry::shape_key(&lateral.0);
        let pc_up: Arc<dyn Curve2d> = Arc::new(occt_geom2d::Geom2dLine::from_pnt_dir(
            p2(0.0, 0.0),
            occt_core::gp::GpDir2d::new(0.0, 1.0).unwrap(),
        ));
        let pc_down: Arc<dyn Curve2d> = Arc::new(occt_geom2d::Geom2dLine::from_pnt_dir(
            p2(2.0 * PI, height),
            occt_core::gp::GpDir2d::new(0.0, -1.0).unwrap(),
        ));
        GeometryRegistry::global().set_edge_pcurve(&seam_up.0, face_key, pc_up);
        GeometryRegistry::global().set_edge_pcurve(&seam_down.0, face_key, pc_down);

        let wire = b.make_wire(&[bottom_circle, seam_up, top_circle, seam_down]);
        let mut lateral = lateral;
        b.add_wire(&mut lateral, &wire);
        lateral
    }

    #[test]
    fn cylinder_lateral_periodic_folding() {
        let radius = 1.0;
        let height = 2.0;
        let face = clean_cylinder_lateral(radius, height);
        let cl = FClass2d::new(&face, 1e-6).expect("classifier");
        let h = height;

        // Interior point is In, and its (u + 2π) periodic image classifies the
        // same.
        let u = 0.5;
        let v = 0.25 * h;
        assert_eq!(cl.perform(p2(u, v)), FaceState::In);
        assert_eq!(cl.perform(p2(u + 2.0 * PI, v)), FaceState::In);
        // A domain-outside left point folds by +2π into the face -> In.
        assert_eq!(cl.perform(p2(-PI, v)), FaceState::In);
        assert_eq!(cl.perform(p2(-PI, v)), cl.perform(p2(PI, v)));
        // Outside above the top edge (v not periodic) stays Out in both images.
        assert_eq!(cl.perform(p2(u, 1.5 * h)), FaceState::Out);
        assert_eq!(cl.perform(p2(u + 2.0 * PI, 1.5 * h)), FaceState::Out);

        // The lateral face is open (bounded caps): the infinite point is Out.
        assert_eq!(cl.perform_infinite_point(), FaceState::Out);
    }

    #[test]
    fn cylinder_primitive_lateral_periodic_equality() {
        let c = BRepPrimCylinder::make_cylinder(1.0, 2.0);
        let lateral = faces_of(&c.solid.0).into_iter().find(|f| !face_is_planar(f)).expect("lateral");
        let cl = FClass2d::new(&lateral, 1e-6).expect("classifier");
        let u = 0.5;
        let v = 0.5;
        // The periodic-image equality holds even for the primitive's (seam-
        // duplicated) lateral wire.
        assert_eq!(cl.perform(p2(u, v)), cl.perform(p2(u + 2.0 * PI, v)));
        assert_eq!(cl.perform_infinite_point(), FaceState::Out);
    }

    #[test]
    fn is_hole_reflects_winding() {
        let b = TopoBuilder::new();
        // CCW ring in the standard-plane UV -> not a hole.
        let ccw = [
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
        ];
        let pln = GpPln::new(GpAx3::standard());
        let surf: Arc<dyn Surface> = Arc::new(GeomPlane::new(pln.clone()));
        let f_ccw = crate::builder_face::make_face_from_wire(&square_edges(&b, &ccw), Some(surf))
            .expect("ccw face");
        let cl = FClass2d::new(&f_ccw, 1e-6).expect("classifier");
        assert!(!cl.is_hole(), "CCW ring is a bounded face, not a hole");

        // CW ring in the same UV -> a hole.
        let cw = [
            GpPnt::new(0.0, 0.0, 0.0),
            GpPnt::new(0.0, 1.0, 0.0),
            GpPnt::new(1.0, 1.0, 0.0),
            GpPnt::new(1.0, 0.0, 0.0),
        ];
        let surf: Arc<dyn Surface> = Arc::new(GeomPlane::new(pln));
        let f_cw = crate::builder_face::make_face_from_wire(&square_edges(&b, &cw), Some(surf))
            .expect("cw face");
        let cl = FClass2d::new(&f_cw, 1e-6).expect("classifier");
        assert!(cl.is_hole(), "CW ring reports a hole face");
    }

    #[test]
    fn infinite_point_closed_periodic_face_is_in() {
        // A full sphere face has no boundary wires -> the classifier sees an
        // empty boundary and the infinite point is In (closed periodic face).
        let sph = crate::primitives::BRepPrimSphere::make_sphere(2.0);
        let face = faces_of(&sph.solid.0).into_iter().next().expect("sphere face");
        let cl = FClass2d::new(&face, 1e-6).expect("classifier");
        assert_eq!(cl.perform_infinite_point(), FaceState::In);
        // Any 2D point is In (no restriction).
        assert_eq!(cl.perform(p2(1.0, 0.5)), FaceState::In);
    }

    #[test]
    fn test_on_restriction_detects_boundary() {
        let ub = unit_box();
        let face = &ub.faces[0];
        let cl = FClass2d::new(face, 1e-6).expect("classifier");
        // A point exactly on the boundary is On; an interior point is In.
        assert_eq!(cl.test_on_restriction(p2(0.5, 0.0), 1e-6), FaceState::On);
        assert_eq!(cl.test_on_restriction(p2(0.5, 0.5), 1e-6), FaceState::In);
        assert_eq!(cl.test_on_restriction(p2(1.5, 0.5), 1e-6), FaceState::Out);
    }
}
