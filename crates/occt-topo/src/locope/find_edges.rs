//! `LocOpe_FindEdges` (`LocOpe_FindEdges.hxx:29-57`, `.cxx:37-409`,
//! `.lxx:23-70`).
//!
//! `Set(FFrom, FTo)` walks the edges of both shapes and pairs them when their
//! curves are the same analytic object: `cxx:56-68` unwraps a
//! `Geom_TrimmedCurve` and keeps only `Geom_Line` / `Geom_Circle` /
//! `Geom_Ellipse` / `Geom_BSplineCurve` / `Geom_BezierCurve` on the `FFrom`
//! side, `cxx:74-88` requires the very same `DynamicType` on the `FTo` side,
//! and `cxx:91-407` then compares the two curves per type. The pairing is
//! recorded in the two parallel lists that the iterator walks `.lxx:34-70`.

use std::sync::Arc;

use occt_core::elib::clib;
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::Curve;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::geom_abs_curve_type::{curve_kind, GeomAbsCurveType};
use crate::shape::{Edge, TopoShape};
use crate::topexp::Explorer;

pub struct FindEdges {
    my_f_from: Option<TopoShape>,
    my_f_to: Option<TopoShape>,
    my_l_from: Vec<Edge>,
    my_l_to: Vec<Edge>,
    index: usize,
}

impl Default for FindEdges {
    fn default() -> Self {
        Self::new()
    }
}

impl FindEdges {
    /// `LocOpe_FindEdges()` (`.lxx:23`): default-constructed, no lists.
    pub fn new() -> Self {
        Self {
            my_f_from: None,
            my_f_to: None,
            my_l_from: Vec::new(),
            my_l_to: Vec::new(),
            index: 0,
        }
    }

    /// `LocOpe_FindEdges(FFrom, FTo)` (`.lxx:27-30`).
    pub fn with_shapes(from: &TopoShape, to: &TopoShape) -> Self {
        let mut f = Self::new();
        f.set(from, to);
        f
    }

    /// `Set(FFrom, FTo)` (`.cxx:37-409`).
    pub fn set(&mut self, from: &TopoShape, to: &TopoShape) {
        self.my_f_from = Some(from.clone());
        self.my_f_to = Some(to.clone());
        self.my_l_from.clear();
        self.my_l_to.clear();

        let edges_from = explorer_edges(from);
        let edges_to = explorer_edges(to);

        for edgf in &edges_from {
            let Some(cf) = BRepTool::edge_curve_world(edgf) else {
                continue;
            };
            let (ff, lf) = BRepTool::edge_parameters(edgf);
            // `cxx:59-63`: `Geom_TrimmedCurve -> BasisCurve()`.
            let cf = basis_of(&cf);
            let tf = curve_kind(cf.as_ref());
            // `cxx:64-70`.
            if !is_compared_kind(tf) {
                continue;
            }
            for edgt in &edges_to {
                let Some(ct) = BRepTool::edge_curve_world(edgt) else {
                    continue;
                };
                let (ft, lt) = BRepTool::edge_parameters(edgt);
                // `cxx:80-85`.
                let ct = basis_of(&ct);
                let tt = curve_kind(ct.as_ref());
                // `cxx:86-89`: `if (Tt != Tf) continue;`
                if tt != tf {
                    continue;
                }
                let tol = CONFUSION;
                let same = match tt {
                    GeomAbsCurveType::Line => same_line(cf.as_ref(), ct.as_ref(), ff, lf, ft, lt, tol),
                    GeomAbsCurveType::Circle => {
                        same_circle(cf.as_ref(), ct.as_ref(), ff, lf, ft, lt, tol)
                    }
                    GeomAbsCurveType::Ellipse => {
                        same_ellipse(cf.as_ref(), ct.as_ref(), ff, lf, ft, lt, tol)
                    }
                    GeomAbsCurveType::BSplineCurve => {
                        same_bspline(cf.as_ref(), ct.as_ref(), edgt, tol)
                    }
                    GeomAbsCurveType::BezierCurve => {
                        same_bezier(cf.as_ref(), ct.as_ref(), tol)
                    }
                    _ => false,
                };
                if same {
                    self.my_l_from.push(edgf.clone());
                    self.my_l_to.push(edgt.clone());
                    break;
                }
            }
        }
    }

    /// `InitIterator()` (`.lxx:34-38`).
    pub fn init_iterator(&mut self) {
        self.index = 0;
    }

    /// `More()` (`.lxx:42-45`).
    pub fn more(&self) -> bool {
        self.index < self.my_l_from.len()
    }

    /// `EdgeFrom()` (`.lxx:49-52`).
    pub fn edge_from(&self) -> Option<&Edge> {
        self.my_l_from.get(self.index)
    }

    /// `EdgeTo()` (`.lxx:56-59`).
    pub fn edge_to(&self) -> Option<&Edge> {
        self.my_l_to.get(self.index)
    }

    /// `Next()` (`.lxx:63-67`): both iterators advance together.
    pub fn next(&mut self) {
        self.index += 1;
    }
}

fn explorer_edges(shape: &TopoShape) -> Vec<Edge> {
    let mut out = Vec::new();
    let mut ex = Explorer::new(shape, ShapeType::Edge);
    while ex.more() {
        if let Some(e) = Edge::wrap(ex.current().clone()) {
            out.push(e);
        }
        ex.next();
    }
    out
}

/// The `Geom_TrimmedCurve -> BasisCurve()` unwrap of `cxx:59-63` /
/// `cxx:80-84`: a trimmed curve is compared through its basis, every other
/// curve is compared as it is.
fn basis_of(c: &Arc<dyn Curve>) -> Arc<dyn Curve> {
    match c.trim_basis() {
        Some(b) => b,
        None => c.clone(),
    }
}

/// The five `DynamicType()` values `cxx:64-70` accepts.
fn is_compared_kind(t: GeomAbsCurveType) -> bool {
    matches!(
        t,
        GeomAbsCurveType::Line
            | GeomAbsCurveType::Circle
            | GeomAbsCurveType::Ellipse
            | GeomAbsCurveType::BSplineCurve
            | GeomAbsCurveType::BezierCurve
    )
}

/// `cxx:93-115`: both endpoints of the `FFrom` interval must lie on the `FTo`
/// line within its own parameter range, then both points are re-evaluated at
/// the `FTo` parameters and compared with `Confusion^2`.
fn same_line(
    cf: &dyn Curve,
    ct: &dyn Curve,
    ff: f64,
    lf: f64,
    ft: f64,
    lt: f64,
    tol: f64,
) -> bool {
    let (Some(lif), Some(lit)) = (cf.gp_line(), ct.gp_line()) else {
        return false;
    };
    let p1 = clib::line_value(&lif, ff);
    let p2 = clib::line_value(&lif, lf);
    let prm1 = clib::line_parameter(&lit, &p1);
    let prm2 = clib::line_parameter(&lit, &p2);
    if !(prm1 >= ft - tol && prm1 <= lt + tol && prm2 >= ft - tol && prm2 <= lt + tol) {
        return false;
    }
    let tol2 = tol * tol;
    let pt = clib::line_value(&lit, prm1);
    if pt.square_distance(&p1) > tol2 {
        return false;
    }
    let pt = clib::line_value(&lit, prm2);
    pt.square_distance(&p2) <= tol2
}

/// `cxx:116-176`: same radius / location, then the start point is folded into
/// the `FTo` period and the second parameter is walked in the direction that
/// matches the `FFrom` tangent (`cxx:138-158`).
fn same_circle(
    cf: &dyn Curve,
    ct: &dyn Curve,
    ff: f64,
    lf: f64,
    ft: f64,
    lt: f64,
    tol: f64,
) -> bool {
    let (Some(cif), Some(cit)) = (cf.gp_circ(), ct.gp_circ()) else {
        return false;
    };
    if (cif.radius() - cit.radius()).abs() > tol
        || cif.location().square_distance(&cit.location()) > tol * tol
    {
        return false;
    }
    let (p1, tgf) = clib::circle_d1(&cif, ff);
    let p2 = clib::circle_value(&cif, lf);

    let mut prm1 = clib::circle_parameter(&cit.position(), &p1);
    let tol2d = PCONFUSION;
    if (prm1 - ft).abs() <= tol2d {
        prm1 = ft;
    }
    prm1 = clib::in_period(prm1, ft, ft + 2.0 * std::f64::consts::PI);
    let (_, tgt) = clib::circle_d1(&cit, prm1);

    let mut prm2 = clib::circle_parameter(&cit.position(), &p2);
    if tgt.dot(&tgf) > 0.0 {
        while prm2 <= prm1 {
            prm2 += 2.0 * std::f64::consts::PI;
        }
    } else {
        if (prm1 - ft).abs() <= occt_core::precision::ANGULAR {
            prm1 += 2.0 * std::f64::consts::PI;
        }
        while prm2 >= prm1 {
            prm2 -= 2.0 * std::f64::consts::PI;
        }
    }
    prm1 >= ft - tol && prm1 <= lt + tol && prm2 >= ft - tol && prm2 <= lt + tol
}

/// `cxx:178-236`: the circle comparison with `MajorRadius` / `MinorRadius`
/// instead of `Radius`, and without the initial `PConfusion` snap (`cxx:193`
/// goes straight to `InPeriod`).
fn same_ellipse(
    cf: &dyn Curve,
    ct: &dyn Curve,
    ff: f64,
    lf: f64,
    ft: f64,
    lt: f64,
    tol: f64,
) -> bool {
    let (Some(cif), Some(cit)) = (cf.gp_ellipse(), ct.gp_ellipse()) else {
        return false;
    };
    if (cif.major_radius - cit.major_radius).abs() > tol
        || (cif.minor_radius - cit.minor_radius).abs() > tol
        || cif.location().square_distance(&cit.location()) > tol * tol
    {
        return false;
    }
    let (p1, tgf) = clib::ellipse_d1(&cif, ff);
    let p2 = clib::ellipse_value(&cif, lf);

    let mut prm1 = clib::ellipse_parameter(&cit.position(), cit.major_radius, cit.minor_radius, &p1);
    prm1 = clib::in_period(prm1, ft, ft + 2.0 * std::f64::consts::PI);
    let (_, tgt) = clib::ellipse_d1(&cit, prm1);

    let mut prm2 = clib::ellipse_parameter(&cit.position(), cit.major_radius, cit.minor_radius, &p2);
    if tgt.dot(&tgf) > 0.0 {
        while prm2 <= prm1 {
            prm2 += 2.0 * std::f64::consts::PI;
        }
    } else {
        if (prm1 - ft).abs() <= occt_core::precision::ANGULAR {
            prm1 += 2.0 * std::f64::consts::PI;
        }
        while prm2 >= prm1 {
            prm2 -= 2.0 * std::f64::consts::PI;
        }
    }
    prm1 >= ft - tol && prm1 <= lt + tol && prm2 >= ft - tol && prm2 <= lt + tol
}

/// `cxx:234-335`: pole count, knot count, poles (against `FTo`'s edge
/// tolerance), knots, multiplicities, rationality and weights.
///
/// The reference compares knots with `(Kf(k) - Kt(k)) > Tol` (signed, not
/// absolute, `cxx:282`) and multiplicities through a `double` conversion
/// (`cxx:286`); both quirks are kept.
fn same_bspline(cf: &dyn Curve, ct: &dyn Curve, edgt: &Edge, tol: f64) -> bool {
    let (Some(bf), Some(bt)) = (cf.as_bspline_curve(), ct.as_bspline_curve()) else {
        return false;
    };
    if bf.nb_poles() != bt.nb_poles() {
        return false;
    }
    if bf.nb_knots() != bt.nb_knots() {
        return false;
    }
    let tol3d = BRepTool::edge_tolerance(edgt);
    for (pf, pt) in bf.poles.iter().zip(bt.poles.iter()) {
        if pf.distance(pt) > tol3d {
            return false;
        }
    }
    for (kf, kt) in bf.distinct_knots().iter().zip(bt.distinct_knots().iter()) {
        if (kf - kt) > tol {
            return false;
        }
    }
    for (mf, mt) in
        bf.multiplicities().iter().zip(bt.multiplicities().iter())
    {
        if (*mf as f64 - *mt as f64).abs() > tol {
            return false;
        }
    }
    if bf.is_rational() != bt.is_rational() {
        return false;
    }
    if bf.is_rational() {
        let (Some(wf), Some(wt)) = (cf.bspline_weights(), ct.bspline_weights()) else {
            return false;
        };
        for (a, b) in wf.iter().zip(wt.iter()) {
            if (a - b).abs() > tol {
                return false;
            }
        }
    }
    true
}

/// `cxx:336-406`: pole count, poles (against `Tol`), rationality and weights.
fn same_bezier(cf: &dyn Curve, ct: &dyn Curve, tol: f64) -> bool {
    let (Some(pf), Some(pt)) = (cf.bezier_poles(), ct.bezier_poles()) else {
        return false;
    };
    if pf.len() != pt.len() {
        return false;
    }
    for (a, b) in pf.iter().zip(pt.iter()) {
        if a.distance(b) > tol {
            return false;
        }
    }
    // `cxx:364-377`: `IsRational()` must agree.
    match (cf.bezier_weights(), ct.bezier_weights()) {
        (None, None) => true,
        (Some(wf), Some(wt)) => {
            // `cxx:381-393`: both weight arrays, `|Wf - Wt| <= Tol`.
            wf.iter().zip(wt.iter()).all(|(a, b)| (a - b).abs() <= tol)
        }
        _ => false,
    }
}
