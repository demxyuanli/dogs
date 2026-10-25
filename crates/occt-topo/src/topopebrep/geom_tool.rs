//! `TopOpeBRep_GeomTool` — builds the DS curve and p-curves of an intersection
//! line for the DS filler.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRep/TopOpeBRep_GeomTool.hxx:30-59`,
//! `.cxx:31-188`.
//!
//! Translation notes:
//! - `L.Curve()` / `L.Curve(min, max)` are the non-const `LineInter` accessors
//!   (`TopOpeBRep_LineInter.cxx:480-508`, `:512-527`), so `l` is `&mut`.
//! - `Standard_ProgramError("TopOpeBRep_GeomTool::MakePrivateCurves")` is the
//!   default arm of both switches; the port returns that message in `Err`.
//! - `MakeBSpline1fromWALKING3d` / `...2d` collect the walking points through
//!   `TopOpeBRep_WPointInterIterator` (`cxx:152-160`, `:172-185`).
//! - `S1` / `S2` are unused by OCCT (the parameters are commented out in the
//!   definition, `cxx:40-42`); the port omits them.

use std::sync::Arc;

use occt_geom2d::curve::Curve2d;

use crate::topopebrep::ds::Curve as DsCurve;
use crate::topopebrep::line_inter::LineInter;
use crate::topopebrep::tool::curve_tool::CurveTool;
use crate::topopebrep::type_line_curve::TypeLineCurve;
use crate::topopebrep::w_point_inter_iterator::WPointInterIterator;

/// `class TopOpeBRep_GeomTool` (`hxx:30-59`).
pub struct GeomTool;

impl GeomTool {
    /// `MakeCurves(min, max, L, S1, S2, C, PC1, PC2)` (`cxx:37-91`): the DS
    /// curve `c` and the p-curves `pc1` / `pc2` of line `l`, whose bounds are
    /// `[min, max]`.
    pub fn make_curves(
        min: f64,
        max: f64,
        l: &mut LineInter,
        c: &mut DsCurve,
        pc1: &mut Option<Arc<dyn Curve2d>>,
        pc2: &mut Option<Arc<dyn Curve2d>>,
    ) -> Result<(), String> {
        let mut is_walk = false;
        let c3d;

        match l.type_line_curve() {
            TypeLineCurve::Walking => {
                // make BSplines of degree 1
                c3d = Self::make_bspline1_from_walking3d(l);
                *pc1 = Self::make_bspline1_from_walking2d(l, 1);
                *pc2 = Self::make_bspline1_from_walking2d(l, 2);
                if pc1.is_some() {
                    c.set_curve1(pc1.clone());
                }
                if pc2.is_some() {
                    c.set_curve2(pc2.clone());
                }
                is_walk = true;
            }
            TypeLineCurve::Line | TypeLineCurve::Circle | TypeLineCurve::Ellipse => {
                c3d = l.curve();
            }
            TypeLineCurve::Parabola | TypeLineCurve::Hyperbola => {
                c3d = l.curve_trimmed(min, max); // Trimmed
            }
            _ => {
                return Err("TopOpeBRep_GeomTool::MakePrivateCurves".to_string());
            }
        }

        let tol = c.tolerance();
        c.define_curve(c3d, tol, is_walk);
        c.set_range(min, max);
        Ok(())
    }

    /// `MakeCurve(min, max, L, C3D)` (`cxx:97-144`): the 3D curve of `l`.
    ///
    /// The conic arm repeats the `typeline` tests inside the `CIRCLE`-to-
    /// `HYPERBOLA` case, exactly as OCCT does (`cxx:118-134`). `C3D` is a
    /// handle that can stay null in OCCT, so the answer is an `Option`; the
    /// throwing default arm is the `Err`.
    pub fn make_curve(
        min: f64,
        max: f64,
        l: &mut LineInter,
    ) -> Result<Option<Arc<dyn occt_geom::curve::Curve>>, String> {
        let typeline = l.type_line_curve();

        let c3d = match typeline {
            TypeLineCurve::Walking => Self::make_bspline1_from_walking3d(l),
            TypeLineCurve::Line => l.curve(),
            TypeLineCurve::Circle | TypeLineCurve::Ellipse | TypeLineCurve::Parabola
            | TypeLineCurve::Hyperbola => {
                if typeline == TypeLineCurve::Circle {
                    l.curve()
                } else if typeline == TypeLineCurve::Ellipse {
                    l.curve()
                } else if typeline == TypeLineCurve::Parabola {
                    l.curve_trimmed(min, max) // Trimmed
                } else {
                    // TypeLineCurve::Hyperbola
                    l.curve_trimmed(min, max) // Trimmed
                }
            }
            _ => return Err("TopOpeBRep_GeomTool::MakePrivateCurves".to_string()),
        };

        Ok(c3d)
    }

    /// `MakeBSpline1fromWALKING3d(L)` (`cxx:148-162`): a degree-1 B-spline
    /// through the walking 3D points of `l`.
    pub fn make_bspline1_from_walking3d(l: &mut LineInter) -> Option<Arc<dyn occt_geom::curve::Curve>> {
        let nbpoints = l.nb_w_point();
        let mut points = Vec::with_capacity(nbpoints.max(0) as usize);
        let mut it = WPointInterIterator::of_line(l);
        it.init();
        while it.more() {
            points.push(it.current_wp().ok()?.value());
            it.next();
        }
        CurveTool::make_bspline1_from_pnt(&points)
    }

    /// `MakeBSpline1fromWALKING2d(L, SI)` (`cxx:166-187`): a degree-1 2D
    /// B-spline through the walking points of `l` on shape `si` (1 or 2).
    ///
    /// `values2d` is filled for every walking point; an `SI` outside `{1, 2}`
    /// leaves the entry untouched, as OCCT does (`cxx:178-185`), so such an
    /// `SI` yields the default `GpPnt2d` for those slots.
    pub fn make_bspline1_from_walking2d(l: &mut LineInter, si: i32) -> Option<Arc<dyn Curve2d>> {
        let nbpoints = l.nb_w_point();
        let mut points = vec![occt_core::gp::GpPnt2d::new(0.0, 0.0); nbpoints.max(0) as usize];
        let mut it = WPointInterIterator::of_line(l);
        it.init();
        let mut ip = 0usize;
        while it.more() {
            if si == 1 {
                points[ip] = it.current_wp().ok()?.value_on_s1();
            } else if si == 2 {
                points[ip] = it.current_wp().ok()?.value_on_s2();
            }
            ip += 1;
            it.next();
        }
        CurveTool::make_bspline1_from_pnt2d(&points)
    }
}
