//! `BRepMesh_CurveTessellator` (`gcpnts`-backed adaptive chord tessellation).

// ---------------------------------------------------------------------------
// CurveTessellator
// ---------------------------------------------------------------------------

use super::prelude::*;
use super::*;

/// Flattens a parametric curve into a deflection-bounded polyline.
///
/// Port of `BRepMesh_CurveTessellator` / `GCPnts_TangentialDeflection`: the curve
/// is seeded with at least `min_points` uniform samples and every segment is
/// then refined while its chord (evaluated at the parameter midpoint) deviates
/// from the true curve by more than `deflection`, or while the angle between the
/// segment end tangents exceeds `angular_deflection` (disabled when infinite).
/// Circles use `PerformCircular`. Other types follow
/// `GCPnts_TangentialDeflection::initialize` (`cxx:415-453`).
#[derive(Clone)]
pub struct CurveTessellator {
    curve: Arc<dyn Curve>,
    deflection: f64,
    angular_deflection: f64,
    min_points: usize,
    min_size: f64,
    pub(super) points: Vec<GpPnt>,
    pub(super) params: Vec<f64>,
}

struct TessCurve<'a>(&'a dyn Curve);

impl CurveSecondDeriv for TessCurve<'_> {
    fn point(&self, u: f64) -> GpPnt {
        self.0.d0(u)
    }
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        self.0.d2(u)
    }
}

impl CurveTessellator {
    /// Tessellate `curve` over its natural parameter range.
    pub fn new(curve: Arc<dyn Curve>, deflection: f64, min_points: usize) -> Self {
        let (a, b) = (curve.first_parameter(), curve.last_parameter());
        Self::from_range(curve, a, b, deflection, min_points)
    }

    /// Tessellate `curve` over `[first, last]` so every chord deviates from the
    /// curve by at most `deflection`. At least `min_points` points are emitted.
    pub fn from_range(
        curve: Arc<dyn Curve>,
        first: f64,
        last: f64,
        deflection: f64,
        min_points: usize,
    ) -> Self {
        Self::from_range_angular(curve, first, last, deflection, f64::INFINITY, min_points)
    }

    /// Tessellate `curve` over `[first, last]` bounding both the chord deviation
    /// (`deflection`) and the angular deviation between consecutive end tangents
    /// (`angular_deflection`, `GCPnts_TangentialDeflection`'s angular term).
    pub fn from_range_angular(
        curve: Arc<dyn Curve>,
        first: f64,
        last: f64,
        deflection: f64,
        angular_deflection: f64,
        min_points: usize,
    ) -> Self {
        Self::from_range_angular_min(
            curve,
            first,
            last,
            deflection,
            angular_deflection,
            min_points,
            CONFUSION,
        )
    }

    /// Same as `from_range_angular` with OCCT `Initialize(..., theMinLen)`.
    pub fn from_range_angular_min(
        curve: Arc<dyn Curve>,
        first: f64,
        last: f64,
        deflection: f64,
        angular_deflection: f64,
        min_points: usize,
        min_size: f64,
    ) -> Self {
        let mut t = Self {
            curve,
            deflection: deflection.max(1e-12),
            angular_deflection,
            min_points: min_points.max(2),
            min_size: min_size.max(CONFUSION),
            points: Vec::new(),
            params: Vec::new(),
        };
        t.build(first, last);
        t
    }

    /// Number of tessellation points.
    pub fn points_nb(&self) -> usize {
        self.points.len()
    }

    /// The `index`-th point (0-based). `None` out of range.
    pub fn point(&self, index: usize) -> Option<GpPnt> {
        self.points.get(index).copied()
    }

    /// The parameter of the `index`-th point (0-based). `None` out of range.
    pub fn parameter(&self, index: usize) -> Option<f64> {
        self.params.get(index).copied()
    }

    /// `(point, parameter)` of the `index`-th solution.
    pub fn value(&self, index: usize) -> Option<(GpPnt, f64)> {
        match (self.points.get(index), self.params.get(index)) {
            (Some(&p), Some(&u)) => Some((p, u)),
            _ => None,
        }
    }

    /// All points.
    pub fn points(&self) -> &[GpPnt] {
        &self.points
    }

    /// All parameters.
    pub fn params(&self) -> &[f64] {
        &self.params
    }

    /// `GCPnts_TangentialDeflection::AddPoint` (`cxx:458-491`).
    pub fn add_point(&mut self, pnt: GpPnt, param: f64, is_replace: bool) -> usize {
        let tol = PCONFUSION;
        let nb = self.params.len();
        for i in 0..nb {
            let dist = self.params[i] - param;
            if dist.abs() <= tol {
                if is_replace {
                    if i < self.points.len() {
                        self.points[i] = pnt;
                    }
                    self.params[i] = param;
                }
                return i;
            } else if dist > tol {
                let pi = i.min(self.points.len());
                self.points.insert(pi, pnt);
                self.params.insert(i, param);
                return i;
            }
        }
        self.points.push(pnt);
        self.params.push(param);
        self.params.len() - 1
    }

    /// `BRepMesh_CurveTessellator::addInternalVertices` (`cxx:193-208`).
    ///
    /// `BRep_Tool::Parameter(V,E)` walks `BRep_TVertex` PointOnCurve, which this
    /// port does not store. Fallback: project `BRep_Tool::Pnt` onto the 3D curve.
    pub fn add_internal_vertices(&mut self, edge: &Edge) {
        let stored = edge
            .0
            .tshape
            .read()
            .expect("poisoned TShape lock")
            .children
            .clone();
        for child in stored {
            if child.shape_type() != ShapeType::Vertex {
                continue;
            }
            if child.orientation() != Orientation::Internal {
                continue;
            }
            let vertex = Vertex(child);
            let pnt = BRepTool::vertex_point_world(&vertex);
            let Some(param) = vertex_parameter_on_edge(&vertex, edge) else {
                continue;
            };
            self.add_point(pnt, param, true);
        }
    }

    /// Measured worst chord deviation. Each segment's chord is sampled at
    /// `samples` interior curve parameters; the furthest distance to the chord
    /// is returned. `samples == 0` samples only the segment midpoint.
    pub fn max_chord_error(&self, samples: usize) -> f64 {
        let mut max_dev = 0.0f64;
        let k = samples.max(1);
        for i in 0..self.points.len().saturating_sub(1) {
            let pa = &self.points[i];
            let pb = &self.points[i + 1];
            let (ua, ub) = if self.params.len() == self.points.len() {
                (self.params[i], self.params[i + 1])
            } else {
                (0.0, 1.0)
            };
            for j in 1..=k {
                let u = ua + (ub - ua) * j as f64 / k as f64;
                let d = point_segment_dist(&self.curve.d0(u), pa, pb);
                if d > max_dev {
                    max_dev = d;
                }
            }
        }
        max_dev
    }

    fn build(&mut self, first: f64, last: f64) {
        let (a, b) = if first.is_finite() && last.is_finite() && last >= first {
            (first, last)
        } else {
            // Unbounded range (e.g. a bare line): fall back to the natural [0,1].
            (0.0, 1.0)
        };
        if b - a < 1e-12 {
            self.points = vec![self.curve.d0(a)];
            self.params = vec![a];
            return;
        }

        // `GCPnts_TangentialDeflection::initialize` dispatches `GeomAbs_Circle`
        // to `PerformCircular` (uniform `ArcAngularStep`). Binary-splitting a
        // 4-point seed on a full circle yields 48 samples at the default 0.5 rad
        // angle; OCCT walks `ceil(span / Du)` instead.
        if let Some(radius) = self.curve.circle_radius() {
            self.build_circular(a, b, radius);
            return;
        }

        // `GCPnts_TangentialDeflection::initialize` (`cxx:415-453`).
        // BSpline intervals are the adaptor-trimmed `NbIntervals(CN)` set
        // (`GeomAdaptor_Curve.cxx:371-413`), not the untrimmed unique knots.
        let adaptor = TessCurve(self.curve.as_ref());
        let two_poles = self
            .curve
            .bspline_poles()
            .or_else(|| self.curve.bezier_poles())
            .is_some_and(|p| p.len() == 2);
        if self.curve.is_line() || two_poles {
            let (params, points) = perform_linear(&adaptor, a, b, self.min_points);
            self.params = params;
            self.points = points;
            return;
        }

        let ang = if self.angular_deflection.is_finite() && self.angular_deflection > 0.0 {
            self.angular_deflection.max(ANGULAR)
        } else {
            std::f64::consts::PI
        };
        let intervals = if let (Some(knots), Some(deg)) =
            (self.curve.bspline_knots(), self.curve.nurbs_degree())
        {
            occt_core::bspl::adaptor_intervals(
                knots,
                deg,
                self.curve.is_periodic(),
                6,
                a,
                b,
                self.curve.resolution(CONFUSION).min(PCONFUSION),
            )
        } else {
            self.curve.parameter_intervals(6)
        };
        let degree_min_nb = self
            .curve
            .nurbs_degree()
            .map(|d| (d + 1).max(self.min_points))
            .unwrap_or(self.min_points);
        let (params, points) = perform_tangential_curve(
            &adaptor,
            a,
            b,
            ang,
            self.deflection,
            self.min_points,
            PCONFUSION,
            self.min_size,
            &intervals,
            degree_min_nb,
        );
        self.params = params;
        self.points = points;
    }

    /// `GCPnts_TangentialDeflection::PerformCircular` + `ArcAngularStep`.
    /// `Initialize` passes `myMinLen` (`cxx:357-360`, `cxx:413`).
    ///
    /// `GeomAdaptor_Curve::load` (`cxx:252-254`) unwraps `Geom_TrimmedCurve` to
    /// the basis, so Circle `U` is radians. Our STEP trim remaps that interval
    /// onto `[0, 1]`; `aDiff` for `ceil(aDiff / Du)` is the basis angle span,
    /// while stored parameters stay in `[first, last]` for SameParameter
    /// Tessellate2d / KeepParam.
    fn build_circular(&mut self, first: f64, last: f64, radius: f64) {
        let ang = if self.angular_deflection.is_finite() && self.angular_deflection > 0.0 {
            self.angular_deflection
        } else {
            std::f64::consts::PI
        };
        let mut du =
            super::super::range_splitter::arc_angular_step(radius, self.deflection, ang, self.min_size);
        let param_span = last - first;
        let angle_span = if let Some((bf, bl)) = self.curve.trimmed_basis_range() {
            let full = self.curve.last_parameter() - self.curve.first_parameter();
            if full.abs() > 1e-16 {
                (bl - bf).abs() * param_span / full
            } else {
                param_span
            }
        } else {
            param_span
        };
        if du <= 1e-12 {
            du = angle_span.abs();
        }
        let diff = angle_span.abs();
        let mut nb = (diff / du).ceil().min(1.0e6) as i32;
        nb = nb.max(self.min_points as i32 - 1).max(1);
        let du_param = param_span / nb as f64;
        self.params.clear();
        self.points.clear();
        let mut u = first;
        for _ in 0..nb {
            self.params.push(u);
            self.points.push(self.curve.d0(u));
            u += du_param;
        }
        self.params.push(last);
        self.points.push(self.curve.d0(last));
    }
}
