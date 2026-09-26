//! The `Extrema_ExtPElC` dispatch and the small shared helpers of the two
//! reduced point/surface engines ([`ExtremaExtPExtS`], [`ExtremaExtPRevS`]).
//!
//! `Extrema_ExtPExtS::PerformExtPElC` (`Extrema_ExtPExtS.cxx:538-563`) and
//! `Extrema_ExtPRevS::PerformExtPElC` (`Extrema_ExtPRevS.cxx:106-131`) both
//! send the five elementary curve types to `Extrema_ExtPElC::Perform` with the
//! primitive's full natural range - line / hyperbola / parabola on
//! `(-Infinite, +Infinite)` (`Extrema_ExtPElC.cxx:58-81`, `:294-389`,
//! `:402-480`), circle / ellipse on `[0, 2*M_PI]` (`:92-190`, `:203-281`) -
//! and leave the solver not-done for every other type
//! (`Extrema_ExtPExtS.cxx:560-562`).
//!
//! The analytic solvers of `Extrema_ExtPElC` are already ported in
//! `crate::extrema_pc` (`extrema_pc/poly_roots.rs:249-383`); the public
//! `point_curve_extrema_all` (`extrema_pc/point_curve.rs:132-150`) dispatches
//! to them over exactly those natural ranges. It does **not** carry `myIsMin`
//! (`Extrema_ExtPElC.cxx:184`, `:277`, `:381`, `:473`), which is
//! recomputed here with OCCT's own per-type rule.
//!
//! Two deviations, both confined to the engines' input:
//!
//! * `point_curve_extrema_all` deduplicates and **sorts by distance**
//!   (`extrema_pc/point_curve.rs:6-19`), where `Extrema_ExtPElC` keeps its root
//!   order. The solution *set* is the same; only the order in which the engines
//!   consume it differs, which matters solely where the `NbExtMax = 4` cap
//!   truncates (`Extrema_ExtPExtS.cxx:270`, `:329`).
//! * For a hyperbola / parabola whose polynomial solve fails, OCCT leaves
//!   `IsDone()` false (`Extrema_ExtPElC.cxx:349-352`); the ported solvers always
//!   reduce the roots they find, so `done` is reported true.

use std::sync::Arc;

use super::prelude::*;
use super::*;

/// One `Extrema_ExtPElC` solution: `Extrema_POnCurv`'s `Parameter` /
/// `Value` (`Extrema_ExtPElC.hxx`) plus `myIsMin` and `mySqDist`.
#[derive(Clone, Copy)]
pub(crate) struct ExtPelcSolution {
    pub u: f64,
    pub point: GpPnt,
    pub sq_dist: f64,
    pub is_min: bool,
}

/// `Extrema_ExtPElC::IsDone()` together with its solutions. `done == false`
/// mirrors the early `return` s of `Extrema_ExtPElC::Perform` that leave
/// `IsDone()` false.
pub(crate) struct ExtPelcResult {
    pub done: bool,
    pub sols: Vec<ExtPelcSolution>,
}

/// `Extrema_ExtPExtS::ProjectPnt` (`Extrema_ExtPExtS.cxx:55-63`): project
/// `point` onto the plane `(plane_loc, plane_dir)` along `dir`.
pub(crate) fn project_pnt(
    plane_loc: &GpPnt,
    plane_dir: &GpDir,
    dir: &GpVec,
    point: &GpPnt,
) -> GpPnt {
    let po = GpVec::from_pnts(point, plane_loc);
    let mut alpha = po.dot(dir);
    alpha /= dir.dot(&GpVec::from_xyz(plane_dir.xyz()));
    point.translated_vec(&dir.multiplied_scalar(alpha))
}

/// `IsOriginalPnt` (`Extrema_ExtPExtS.cxx:67-77`, `Extrema_ExtPRevS.cxx:168-178`):
/// true when `p` is farther than `Precision::Confusion()` from all `nb`
/// already stored points.
pub(crate) fn is_original_pnt(p: &GpPnt, points: &[ExtremaPOnSurf], nb: usize) -> bool {
    for pt in &points[..nb] {
        if pt.value().distance(p) <= CONFUSION {
            return false;
        }
    }
    true
}

/// `gp_Pnt::IsEqual(theOther, Precision::Confusion())`.
pub(crate) fn pnt_equal(a: &GpPnt, b: &GpPnt) -> bool {
    a.distance(b) <= CONFUSION
}

/// `gp_Pln::SquareDistance(gp_Pnt)` for the plane `(loc, normal)`.
pub(crate) fn plane_square_distance(loc: &GpPnt, normal: &GpDir, p: &GpPnt) -> f64 {
    let d = GpVec::from_pnts(loc, p).dot(&GpVec::from_xyz(normal.xyz()));
    d * d
}

/// `gp_Vec::AngleWithRef(V2, VRef)` (`gp_Vec.cxx:96-113`): normalise both
/// vectors and take the signed `gp_Dir::AngleWithRef` (`gp_Dir.cxx:55-84`).
/// A zero vector raises `Standard_ConstructionError` in OCCT and returns 0.0
/// here, the convention `clib::ellipse_parameter` already uses
/// (`crates/occt-core/src/elib/clib.rs:463-468`).
pub(crate) fn vec_angle_with_ref(v1: &GpVec, v2: &GpVec, vref: &GpVec) -> f64 {
    let (Ok(d1), Ok(d2), Ok(dr)) = (
        GpDir::from_vec(v1),
        GpDir::from_vec(v2),
        GpDir::from_vec(vref),
    ) else {
        return 0.0;
    };
    d1.angle_with_ref(&d2, &dr)
}

/// The component of `p - loc` perpendicular to the axis `axis`: the vector
/// whose magnitude is `gp_Lin(axis).Distance(p)`, the quantity
/// `Extrema_ExtPElC::Perform` tests for circle / ellipse
/// (`Extrema_ExtPElC.cxx:132-136`, `:242-250`).
pub(crate) fn point_to_axis(p: &GpPnt, loc: &GpPnt, axis: &GpAx1) -> GpVec {
    let v = GpVec::from_pnts(loc, p);
    let n = GpVec::from_xyz(axis.direction().xyz());
    v.subtracted(&n.multiplied_scalar(v.dot(&n)))
}

/// `Adaptor3d_Curve::GetType()` membership for the five types
/// `Extrema_ExtPElC` handles (`Extrema_ExtPExtS.cxx:567-585`,
/// `Extrema_ExtPRevS.cxx:135-150`). The elementary queries of a trimmed view
/// forward to its basis, exactly like `GeomAdaptor_Curve::load`
/// (`GeomAdaptor_Curve.cxx:251-255`).
pub(crate) fn is_elementary_curve(c: &dyn Curve) -> bool {
    c.gp_line().is_some()
        || c.is_line()
        || c.gp_circ().is_some()
        || c.gp_ellipse().is_some()
        || c.gp_hyperbola().is_some()
        || c.gp_parabola().is_some()
}

/// The part of `gp_Ax2` the reduced engines read from `GetPosition`:
/// `Direction()` for the planarity tests and `Location()` / `Direction()`
/// for `ProjectPnt` (`Extrema_ExtPExtS.cxx:509-534`,
/// `Extrema_ExtPRevS.cxx:38-79`).
#[derive(Clone, Copy)]
pub(crate) struct CurvePlane {
    pub loc: GpPnt,
    pub dir: GpDir,
}

impl Default for CurvePlane {
    /// `gp_Ax2()` is `gp::XOY()` (`Extrema_ExtPExtS.cxx:531-532`,
    /// `Extrema_ExtPRevS.cxx:57`, `:77`).
    fn default() -> Self {
        let a = GpAx2::standard();
        CurvePlane {
            loc: a.location(),
            dir: a.direction(),
        }
    }
}

/// `GeomAdaptor_Curve::load` (`GeomAdaptor_Curve.cxx:251-255`) unwraps a
/// `Geom_TrimmedCurve` and keeps the analytic basis, so `Extrema_ExtPElC` is
/// driven in **basis** parameters. This port's `GeomTrimmedCurve` view remaps
/// the trim onto `[0, 1]` (`crates/occt-geom/src/trimmed.rs:44-81`), and the
/// surface that owns it uses that view; the solver's parameters therefore have
/// to be carried back into the view the surface evaluates.
pub(crate) struct CurveMap {
    /// The analytic basis curve handed to `point_curve_extrema_all`.
    pub analytic: Arc<dyn Curve>,
    /// Analytic basis range.
    pub b1: f64,
    pub b2: f64,
    /// The stored curve's own range.
    pub c1: f64,
    pub c2: f64,
    /// Whether the stored curve is a trimmed view whose parameters have to be
    /// remapped; an untrimmed curve is its own analytic basis and both
    /// parameterisations are identical (its range may even be infinite, so an
    /// affine remap would be undefined).
    trimmed: bool,
}

impl CurveMap {
    /// Builds the map for the curve the surface stores.
    pub(crate) fn of(c: &Arc<dyn Curve>) -> Self {
        match c.untrimmed_basis() {
            Some((basis, b1, b2)) => CurveMap {
                analytic: basis,
                b1,
                b2,
                c1: c.first_parameter(),
                c2: c.last_parameter(),
                trimmed: true,
            },
            None => {
                let (c1, c2) = (c.first_parameter(), c.last_parameter());
                CurveMap {
                    analytic: c.clone(),
                    b1: c1,
                    b2: c2,
                    c1,
                    c2,
                    trimmed: false,
                }
            }
        }
    }

    /// Analytic parameter -> stored-curve parameter.
    pub(crate) fn to_curve(&self, u: f64) -> f64 {
        if !self.trimmed {
            return u;
        }
        let den = self.b2 - self.b1;
        if den.abs() <= 0.0 {
            return self.c1;
        }
        self.c1 + (u - self.b1) / den * (self.c2 - self.c1)
    }
}

/// `Extrema_ExtPElC::Perform` for the five elementary types, with the bounds
/// `Extrema_ExtPExtS` / `Extrema_ExtPRevS` use
/// (`Extrema_ExtPExtS.cxx:538-563`). `c` must be the analytic basis curve
/// ([CurveMap::analytic]).
pub(crate) fn ext_pelc_perform(c: &dyn Curve, p: &GpPnt, tol: f64) -> ExtPelcResult {
    let not_done = || ExtPelcResult {
        done: false,
        sols: Vec::new(),
    };

    // Line (Extrema_ExtPElC.cxx:58-81): one solution, always a minimum; IsDone
    // is false when the parameter falls outside [Uinf, Usup], which here is
    // (-Infinite, +Infinite).
    if c.gp_line().is_some() || c.is_line() {
        let pairs = crate::extrema_pc::point_curve_extrema_all(c, p);
        let done = !pairs.is_empty();
        return ExtPelcResult {
            done,
            sols: pairs
                .into_iter()
                .map(|e| sol_of(p, e.u1, e.p2, true))
                .collect(),
        };
    }

    // Circle (Extrema_ExtPElC.cxx:92-190): IsDone is false when P projects onto
    // the circle's axis (cxx:133-136); NoSol == 0, the near point, is the
    // minimum and NoSol == 1 the maximum (cxx:184).
    if let Some(gc) = c.gp_circ() {
        if point_to_axis(p, &gc.location(), &gc.axis()).magnitude() < tol {
            return not_done();
        }
        let pairs = crate::extrema_pc::point_curve_extrema_all(c, p);
        // point_curve_extrema_all sorts ascending by distance, so index 0 is
        // the near point.
        let sols = pairs
            .into_iter()
            .enumerate()
            .map(|(i, e)| sol_of(p, e.u1, e.p2, i == 0))
            .collect();
        return ExtPelcResult { done: true, sols };
    }

    // Ellipse (Extrema_ExtPElC.cxx:203-281): done unless P lies on the axis of
    // a circle-like ellipse (cxx:244-250); myIsMin compares against C(Us + 0.1)
    // (cxx:276-277).
    if let Some(e) = c.gp_ellipse() {
        let on_axis = point_to_axis(p, &e.location(), e.axis()).magnitude() < tol;
        if on_axis && (e.major_radius() - e.minor_radius()).abs() < tol {
            return not_done();
        }
        let pairs = crate::extrema_pc::point_curve_extrema_all(c, p);
        let sols = pairs
            .into_iter()
            .map(|pair| {
                let is_min = p.square_distance(&pair.p2)
                    < p.square_distance(&clib::ellipse_value(&e, pair.u1 + 0.1));
                sol_of(p, pair.u1, pair.p2, is_min)
            })
            .collect();
        return ExtPelcResult { done: true, sols };
    }

    // Hyperbola (Extrema_ExtPElC.cxx:294-389): myIsMin compares against
    // C(Us + 1) (cxx:381). The polynomial solver failing (cxx:349-352) is the
    // only not-done path; the port's quartic arm always reduces the roots it
    // finds (extrema_pc/poly_roots.rs:321-351), so IsDone is reported true.
    if let Some(h) = c.gp_hyperbola() {
        let pairs = crate::extrema_pc::point_curve_extrema_all(c, p);
        let sols = pairs
            .into_iter()
            .map(|pair| {
                let is_min = p.square_distance(&pair.p2)
                    < p.square_distance(&clib::hyperbola_value(&h, pair.u1 + 1.0));
                sol_of(p, pair.u1, pair.p2, is_min)
            })
            .collect();
        return ExtPelcResult { done: true, sols };
    }

    // Parabola (Extrema_ExtPElC.cxx:402-480): myIsMin compares against
    // C(Us + 1) (cxx:473); same not-done remark as the hyperbola.
    if let Some(pa) = c.gp_parabola() {
        let pairs = crate::extrema_pc::point_curve_extrema_all(c, p);
        let sols = pairs
            .into_iter()
            .map(|pair| {
                let is_min = p.square_distance(&pair.p2)
                    < p.square_distance(&clib::parabola_value(&pa, pair.u1 + 1.0));
                sol_of(p, pair.u1, pair.p2, is_min)
            })
            .collect();
        return ExtPelcResult { done: true, sols };
    }

    // default: return; (Extrema_ExtPExtS.cxx:560-562,
    // Extrema_ExtPRevS.cxx:128-130).
    not_done()
}

/// `Extrema_POnCurv(Us, Cu)` + `mySqDist` / `myIsMin`.
fn sol_of(p: &GpPnt, u: f64, point: GpPnt, is_min: bool) -> ExtPelcSolution {
    ExtPelcSolution {
        u,
        point,
        sq_dist: p.square_distance(&point),
        is_min,
    }
}
