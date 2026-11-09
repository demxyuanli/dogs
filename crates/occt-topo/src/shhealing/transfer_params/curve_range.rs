//! `CurveRange`, the `CopyRanges` transfer and its parameter helpers.

use super::prelude::*;

/// `GeomAdaptor_Curve(C, first, last)` (`GeomAdaptor_Curve.cxx:239-245`,
/// `cxx:679-691`): the same parameter space as `curve`, restricted to
/// `[first, last]`. `Value(u)` is `curve.d0(u)`, never a rescaled one, which is
/// what `NextProject` / `Project` expect for the `GeomAdaptor_Curve GAC(C3d,
/// first, last)` built at `Proj.cxx:400` (and for the COS adaptor built from a
/// representation range at `Proj.cxx:462-464`).
pub(super) struct CurveRange {
    pub(super) curve: Arc<dyn Curve>,
    pub(super) first: f64,
    pub(super) last: f64,
}

impl CurveRange {
    pub(super) fn new(curve: Arc<dyn Curve>, first: f64, last: f64) -> Self {
        Self { curve, first, last }
    }
}

impl Curve for CurveRange {
    fn d0(&self, u: f64) -> GpPnt {
        self.curve.d0(u)
    }
    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        self.curve.d1(u)
    }
    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        self.curve.d2(u)
    }
    fn first_parameter(&self) -> f64 {
        self.first
    }
    fn last_parameter(&self) -> f64 {
        self.last
    }
    fn is_periodic(&self) -> bool {
        self.curve.is_periodic()
    }
    fn period(&self) -> f64 {
        self.curve.period()
    }
    fn continuity(&self) -> u8 {
        self.curve.continuity()
    }
    fn circle_radius(&self) -> Option<f64> {
        self.curve.circle_radius()
    }
    fn gp_circ(&self) -> Option<occt_core::gp::GpCirc> {
        self.curve.gp_circ()
    }
    fn gp_ellipse(&self) -> Option<occt_core::gp::GpElips> {
        self.curve.gp_ellipse()
    }
    fn is_line(&self) -> bool {
        self.curve.is_line()
    }
    fn nurbs_degree(&self) -> Option<usize> {
        self.curve.nurbs_degree()
    }
    fn bspline_poles(&self) -> Option<&[GpPnt]> {
        self.curve.bspline_poles()
    }
    fn bspline_weights(&self) -> Option<&[f64]> {
        self.curve.bspline_weights()
    }
    fn bezier_poles(&self) -> Option<&[GpPnt]> {
        self.curve.bezier_poles()
    }
    fn resolution(&self, r3d: f64) -> f64 {
        self.curve.resolution(r3d)
    }
    fn transform(&mut self, _t: &occt_core::gp::GpTrsf) {}
    fn reverse(&mut self) {}
    fn clone_dyn(&self) -> Box<dyn Curve> {
        Box::new(CurveRange::new(self.curve.clone(), self.first, self.last))
    }
}

/// `ShapeBuild_Edge::IsPeriodic(handle(Geom_Curve))`
/// (`ShapeBuild_Edge.cxx:164-183`): unwrap `Geom_OffsetCurve` /
/// `Geom_TrimmedCurve` down to the basis, then `IsPeriodic`. Returns the basis
/// period with the basis `FirstParameter()` / `LastParameter()` the caller
/// reads at `cxx:301-303`.
///
/// The `Geom_OffsetCurve` unwrap (`cxx:170-177`) has no `Curve` counterpart in
/// the port (no offset query), so only the trim unwrap runs
/// (`Curve::untrimmed_basis`, `GeomAdaptor_Curve.cxx:252-254`).
pub(super) fn periodic_basis(curve: &Arc<dyn Curve>) -> Option<(f64, f64, f64)> {
    let (basis, first, last) = match curve.untrimmed_basis() {
        Some((basis, f, l)) => (basis, f, l),
        None => (curve.clone(), curve.first_parameter(), curve.last_parameter()),
    };
    if basis.is_periodic() {
        Some((basis.period(), first, last))
    } else {
        None
    }
}

/// `ShapeBuild_Edge::IsPeriodic(handle(Geom2d_Curve))`
/// (`ShapeBuild_Edge.cxx:185-204`) plus the range read at `cxx:311-317`.
pub(super) fn periodic_basis_2d(curve: &Arc<dyn Curve2d>) -> Option<(f64, f64, f64)> {
    let mut basis: &dyn Curve2d = curve.as_ref();
    // `Geom2d_OffsetCurve` / `Geom2d_TrimmedCurve` unwrap (`cxx:191-201`).
    while let Some(inner) = basis.trimmed_basis().or_else(|| basis.offset_basis()) {
        basis = inner;
    }
    if basis.is_periodic() {
        Some((basis.period(), basis.first_parameter(), basis.last_parameter()))
    } else {
        None
    }
}

/// `ShapeBuild_Edge::CopyRanges(toedge, fromedge, alpha, beta)`
/// (`ShapeBuild_Edge.cxx:206-334`).
///
/// Walks the `fromedge` representations (`cxx:226-232`) and writes the matching
/// representation of `toedge` (`cxx:264-283`): the `BRep_Curve3D` range
/// (`cxx:285-289`, `cxx:326-331`) and every `BRep_CurveOnSurface` range
/// (`cxx:310-318`, `cxx:527-532`). A periodic basis whose shifted range leaves
/// its own range is re-shifted by a whole period, and `SameRange` /
/// `SameParameter` are cleared on `toedge` (`cxx:291-328`).
///
/// The port matches representations by the face key of the pcurve map, which IS
/// the surface identity `surface != toGC->Surface() || L != toGC->Location()`
/// tested at `cxx:281`.
pub(crate) fn copy_ranges(to: &Edge, from: &Edge, alpha: f64, beta: f64) {
    let reg = GeometryRegistry::global();
    let Some(from_geom) = reg.edge_geom(&from.0) else {
        return;
    };
    let Some(mut to_geom) = reg.edge_geom(&to.0) else {
        return;
    };

    // `cxx:234-245`: the `BRep_Curve3D` representation (`isC3d`).
    {
        let first = from_geom.first; // `cxx:285-286`
        let last = from_geom.last;
        let len = last - first;
        let mut new_first = first + alpha * len; // `cxx:289-290`
        let mut new_last = first + beta * len;
        if let Some((period, basis_first, basis_last)) = periodic_basis(&to_geom.curve) {
            // `cxx:321-328`.
            if ((new_first - basis_first).abs() > PCONFUSION && new_first < basis_first)
                || new_first >= basis_last
            {
                let shift = adjust_by_period(new_first, 0.5 * (basis_first + basis_last), period);
                new_first += shift;
                new_last += shift;
                reg.set_same_range(&to.0, false);
                reg.set_same_parameter(&to.0, false);
            }
        }
        to_geom.first = new_first;
        to_geom.last = new_last;
    }

    // `cxx:248-282`: every `BRep_CurveOnSurface` representation.
    let from_keys: Vec<usize> = from_geom.pcurves.keys().copied().collect();
    for key in from_keys {
        if !to_geom.pcurves.contains_key(&key) {
            continue; // `cxx:281`
        }
        let (first, last) = from_geom
            .pcurve_ranges
            .get(&key)
            .copied()
            .unwrap_or((from_geom.first, from_geom.last));
        let len = last - first;
        let mut new_first = first + alpha * len;
        let mut new_last = first + beta * len;
        let pcurve = to_geom.pcurves.get(&key).and_then(|v| v.first().cloned());
        if let Some((period, basis_first, basis_last)) = pcurve.and_then(|pc| periodic_basis_2d(&pc))
        {
            if ((new_first - basis_first).abs() > PCONFUSION && new_first < basis_first)
                || new_first >= basis_last
            {
                let shift = adjust_by_period(new_first, 0.5 * (basis_first + basis_last), period);
                new_first += shift;
                new_last += shift;
                reg.set_same_range(&to.0, false);
                reg.set_same_parameter(&to.0, false);
            }
        }
        to_geom.pcurve_ranges.insert(key, (new_first, new_last));
    }

    reg.set_edge(&to.0, to_geom);
}

/// `CorrectParameter` (`ShapeAnalysis_TransferParametersProj.cxx:255-281`).
///
/// Fully ported (batch 89, task T-15): the `Geom2d_BSplineCurve` knot snap
/// (`cxx:268-279`) runs through the new `Curve2d::bspline_knots2d` query — the
/// blocker recorded here was that `Curve2d` exposed no knot sequence.
pub(super) fn correct_parameter(c2d: &Arc<dyn Curve2d>, param: f64) -> f64 {
    correct_parameter_of(c2d.as_ref(), param)
}

/// `CorrectParameter` (`Proj.cxx:256-281`) on a borrowed curve, so the two
/// unwrap arms can recurse exactly as OCCT does.
///
/// * `Geom2d_TrimmedCurve` / `Geom2d_OffsetCurve` (`cxx:258-266`): recurse on
///   `BasisCurve()`.
/// * `Geom2d_BSplineCurve` (`cxx:268-279`): return the **first** knot within
///   `Precision::PConfusion()` of `param` (OCCT walks `Knot(j)` for
///   `j = FirstUKnotIndex()..LastUKnotIndex()`, i.e. the distinct knots; the
///   port's 2D B-spline is never periodic, `Curve2d::is_periodic` keeps its
///   default `false`, so that range is every distinct knot).
/// * Everything else: `param` unchanged.
pub(super) fn correct_parameter_of(c2d: &dyn Curve2d, param: f64) -> f64 {
    if let Some(basis) = c2d.trimmed_basis() {
        return correct_parameter_of(basis, param);
    }
    if let Some(basis) = c2d.offset_basis() {
        return correct_parameter_of(basis, param);
    }
    if let Some(knots) = c2d.bspline_knots2d() {
        let (distinct, _mults) = occt_core::bspl::unperiodize::distinct_knots_and_mults(knots);
        for valknot in distinct {
            if (valknot - param).abs() < PCONFUSION {
                return valknot;
            }
        }
    }
    param
}

/// `Precision::IsInfinite(p.X()) || IsInfinite(p.Y()) || IsInfinite(p.Z())`
/// (`Proj.cxx:319-321`, `:325-327`, `:342-343`, `:348-351`).
pub(super) fn point_is_infinite(p: &GpPnt) -> bool {
    Precision::is_infinite(p.x()) || Precision::is_infinite(p.y()) || Precision::is_infinite(p.z())
}
