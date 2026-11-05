//! `BRepFill_TrimEdgeTool::IntersectWith` (`cxx:272-628`).
//!
//! The 3D points answered are `gp_Pnt(bisector parameter, edge1 parameter,
//! edge2 parameter)`: `X` is the parameter on `myBis`, `Y` the parameter on
//! `Edge1`, `Z` the parameter on `Edge2` (`cxx:620-627`).
//!
//! Step correspondence:
//! * `cxx:285-298` — the two edge p-curves and their adaptors,
//! * `cxx:300-303` — the first intersection pass at `Precision::Confusion()`,
//! * `cxx:305-371` — the tolerance / iteration counts, the "both curves are far
//!   from the bisector line" early exit,
//! * `cxx:372-390` — the `SeanceDeRattrapage` retry loop with a growing
//!   tolerance,
//! * `cxx:402-468` — the `0 / 1`-point fix-ups that project the dangling side
//!   onto the bisector point,
//! * `cxx:470-537` — the pairing of the two sequences on the bisector
//!   parameter,
//! * `cxx:539-617` — the `GeomAbs_Intersection` restriction to the initial
//!   edge's parameter range,
//! * `cxx:620-627` — the `Z` coordinates.

use std::f64::consts::PI;

use occt_core::elib::clib2d::parameter_lin2d;
use occt_core::gp::GpPnt;
use occt_core::kernel::geomabs::{CurveType, JoinType};
use occt_geom2d::adaptor_curve::Geom2dAdaptorCurve;

use crate::abs::ShapeType;
use crate::geom_abs_curve_type::GeomAbsCurveType;
use crate::shape::{Edge, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topopebrep::build::tools2d::{first_vertex, last_vertex};
use crate::topopebrep::tool::adaptor::CurveAdaptor;

use super::eval_parameters::{eval_parameters, eval_parameters_bis};
use super::{Error, TrimEdgeTool};

impl TrimEdgeTool {
    /// `IntersectWith(Edge1, Edge2, InitShape1, InitShape2, End1, End2,
    /// theJoinType, IsOpenResult, Params)` (`cxx:272-628`).
    ///
    /// `Params` is cleared first (`cxx:283`) and answered ordered by increasing
    /// bisector parameter.
    #[allow(clippy::too_many_arguments)]
    pub fn intersect_with(
        &self,
        edge1: &Edge,
        edge2: &Edge,
        init_shape1: &TopoShape,
        init_shape2: &TopoShape,
        end1: Option<&Vertex>,
        end2: Option<&Vertex>,
        join_type: JoinType,
        is_open_result: bool,
    ) -> Result<Vec<GpPnt>, Error> {
        // `cxx:283`: `Params.Clear();`.
        let mut params: Vec<GpPnt> = Vec::new();

        // `cxx:285-298`: the p-curves and their adaptors.
        let ac1 = edge_adaptor(edge1)?;
        let ac2 = edge_adaptor(edge2)?;

        // `cxx:300-303`: the first pass.
        let mut points2: Vec<GpPnt> = Vec::new();
        eval_parameters(&self.bis, &ac1, &mut params)?;
        eval_parameters(&self.bis, &ac2, &mut points2)?;

        // `cxx:305-313`.
        let mut seance_de_rattrapage: i32 = 0;
        let mut tol_init = 1.0e-9;
        let mut nn: i32 = 7;
        if !is_line_or_circle(&ac1) || !is_line_or_circle(&ac2) {
            tol_init = 1.0e-8;
            nn = 6;
        }

        // `cxx:317-371`: no intersection at all and a line bisector — give up
        // when both curves sit farther than `dmax` from the line.
        if params.is_empty() && points2.is_empty() && self.bis.get_type() == Some(CurveType::Line) {
            let mut dmax = tol_init;
            for _ in 0..nn {
                dmax *= 10.0;
            }
            dmax *= dmax;
            let an_l = self.bis.line().ok_or(Error::NullCurve)?;
            let is_far1 = far_from_line(&an_l, &ac1, &self.bis, dmax)?;
            let is_far2 = far_from_line(&an_l, &ac2, &self.bis, dmax)?;
            if is_far1 && is_far2 {
                return Ok(params);
            }
        }

        // `cxx:372-390`: retry with a growing tolerance while the two
        // sequences disagree in length (or are both empty).
        while seance_de_rattrapage < nn
            && (points2.len() != params.len() || (points2.is_empty() && params.is_empty()))
        {
            params.clear();
            points2.clear();
            tol_init *= 10.0;
            eval_parameters_bis(&self.bis, &ac1, &mut params, tol_init)?;
            eval_parameters_bis(&self.bis, &ac2, &mut points2, tol_init)?;
            seance_de_rattrapage += 1;
        }

        // `cxx:402-433`: `Params` empty, one point on `Edge2`.
        if params.is_empty() && points2.len() == 1 {
            let dmax = 0.25 * self.offset * self.offset;
            let t_bis = points2[0].x();
            let pbis = self.bis.value(t_bis).ok_or(Error::NullCurve)?;
            let mut t = ac1.first_parameter();
            let pc = ac1.value(t).ok_or(Error::NullCurve)?;
            let dmin = pc.square_distance(&pbis);
            let mut p = GpPnt::new(t_bis, t, 0.0);
            if dmin < dmax {
                params.push(p);
            }
            t = ac1.last_parameter();
            let pc = ac1.value(t).ok_or(Error::NullCurve)?;
            let dmin1 = pc.square_distance(&pbis);
            if dmin > dmin1 && dmin1 < dmax {
                p.set_y(t);
                if params.is_empty() {
                    params.push(p);
                } else {
                    params[0] = p;
                }
            }
        } else if params.len() == 1 && points2.is_empty() {
            // `cxx:434-467`: `Params` has one point, `Edge2` none.
            let dmax = 0.25 * self.offset * self.offset;
            let t_bis = params[0].x();
            let pbis = self.bis.value(t_bis).ok_or(Error::NullCurve)?;
            let mut t = ac2.first_parameter();
            let pc = ac2.value(t).ok_or(Error::NullCurve)?;
            let dmin = pc.square_distance(&pbis);
            let mut p = GpPnt::new(t_bis, t, 0.0);
            if dmin < dmax {
                points2.push(p);
            }
            t = ac2.last_parameter();
            let pc = ac2.value(t).ok_or(Error::NullCurve)?;
            let dmin1 = pc.square_distance(&pbis);
            if dmin > dmin1 && dmin1 < dmax {
                p.set_y(t);
                if points2.is_empty() {
                    points2.push(p);
                } else {
                    points2[0] = p;
                }
            }
        }

        // `cxx:479-488`: a single pair is taken as is, with the second edge's
        // parameter in `Z`.
        if params.len() == 1 && points2.len() == 1 {
            params[0].set_z(points2[0].y());
            return Ok(params);
        }

        // `cxx:490-527`: drop the pairs whose bisector parameters differ by
        // less than `TolInit`, and the ones a `Tol` apart in the wrong order.
        let tol = 4.0 * 100.0 * occt_core::precision::PCONFUSION;
        let mut i = 0usize;
        while i < params.len().min(points2.len()) {
            let p1 = params[i];
            let p2 = points2[i];
            let p1x_p2x = (p1.x() - p2.x()).abs();
            if p1x_p2x > tol {
                if p1x_p2x > tol_init {
                    i += 1;
                } else if p1.x() < p2.x() {
                    params.remove(i);
                } else {
                    points2.remove(i);
                }
            } else {
                i += 1;
            }
        }

        // `cxx:529-536`: the two sequences end up the same length.
        if params.len() > points2.len() {
            params.truncate(points2.len());
        } else if params.len() < points2.len() {
            points2.truncate(params.len());
        }

        let nb_points = params.len();

        // `cxx:539-582`: the initial parameter range of `InitShape1` when it is
        // a circle and the join is an intersection.
        let mut init_fpar = -f64::MAX;
        let mut init_lpar = f64::MAX;
        if nb_points > 1
            && join_type == JoinType::Intersection
            && init_shape1.shape_type() != ShapeType::Vertex
            && init_shape2.shape_type() != ShapeType::Vertex
        {
            let init_edge1 = Edge::wrap(init_shape1.clone()).ok_or(Error::TypeMismatch)?;
            let mut to_extend_first_par = true;
            let mut to_extend_last_par = true;
            if is_open_result {
                // `cxx:554-561`.
                let v1 = first_vertex(&init_edge1.0);
                let v2 = last_vertex(&init_edge1.0);
                if same_as_any(&v1, end1, end2) {
                    to_extend_first_par = false;
                }
                if same_as_any(&v2, end1, end2) {
                    to_extend_last_par = false;
                }
            }
            // `cxx:563-581`.
            let ic1 = CurveAdaptor::of_edge(&init_edge1);
            if ic1.get_type() == GeomAbsCurveType::Circle {
                let delta = 2.0 * PI - ic1.last_parameter() + ic1.first_parameter();
                if to_extend_first_par && to_extend_last_par {
                    init_fpar = ac1.first_parameter() + delta / 2.0;
                } else if to_extend_first_par {
                    init_fpar = ac1.first_parameter() + delta;
                } else if to_extend_last_par {
                    init_fpar = ac1.first_parameter();
                }
                init_lpar = init_fpar + ic1.last_parameter() - ic1.first_parameter();
            }
        }

        // `cxx:584-617`: keep the minimal-parameter point plus every point
        // inside the initial range.
        if nb_points > 1 && join_type == JoinType::Intersection {
            let mut imin = 0usize;
            for i in 1..nb_points {
                if params[i].x() < params[imin].x() {
                    imin = i;
                }
            }
            let params_copy = params.clone();
            let points2_copy = points2.clone();
            params.clear();
            points2.clear();
            for (i, p) in params_copy.iter().enumerate() {
                if imin == i || (p.y() >= init_fpar && p.y() <= init_lpar) {
                    params.push(*p);
                    points2.push(points2_copy[i]);
                }
            }
        }

        // `cxx:620-627`: `Z` is the parameter on the second edge.
        let nb_points = params.len();
        for i in 0..nb_points {
            params[i].set_z(points2[i].y());
        }
        Ok(params)
    }
}

/// `BRep_Tool::CurveOnSurface(E, C, S, L, First, Last)` (`BRep_Tool.cxx:301-315`)
/// for one edge: the first CurveOnSurface representation's p-curve and COS
/// range.
///
/// The port's representations are keyed by surface handle
/// (`GeometryRegistry::edge_pcurve_reps`), so "first" is the first entry of
/// that list rather than OCCT's representation order; both edges of
/// `BRepFill_OffsetWire` carry one p-curve on the common face.
fn edge_adaptor(edge: &Edge) -> Result<Geom2dAdaptorCurve, Error> {
    let mut reps = GeometryRegistry::global().edge_pcurve_reps(&edge.0);
    if reps.is_empty() {
        return Err(Error::NullCurve);
    }
    let (_surf, pc, first, last) = reps.remove(0);
    Ok(Geom2dAdaptorCurve::with_range(pc, first, last))
}

/// `cxx:307-310`: `AC.GetType() != GeomAbs_Circle && AC.GetType() != GeomAbs_Line`.
fn is_line_or_circle(ac: &Geom2dAdaptorCurve) -> bool {
    matches!(
        ac.get_type(),
        Some(CurveType::Circle) | Some(CurveType::Line)
    )
}

/// `cxx:334-349` / `cxx:351-366`: whether every endpoint of `ac` inside the
/// bisector's range sits farther than `dmax` (squared) from the bisector line.
///
/// Both endpoints are tested in sequence, the first one assigning `d` and the
/// second one taking the minimum (`cxx:338-347`).
fn far_from_line(
    an_l: &occt_core::gp::GpLin2d,
    ac: &Geom2dAdaptorCurve,
    bis: &Geom2dAdaptorCurve,
    dmax: f64,
) -> Result<bool, Error> {
    let mut d = f64::MAX;
    let mut a_p = ac.value(ac.first_parameter()).ok_or(Error::NullCurve)?;
    let mut par = parameter_lin2d(an_l, &a_p);
    if par >= bis.first_parameter() && par <= bis.last_parameter() {
        d = an_l.square_distance(&a_p);
    }
    a_p = ac.value(ac.last_parameter()).ok_or(Error::NullCurve)?;
    par = parameter_lin2d(an_l, &a_p);
    if par >= bis.first_parameter() && par <= bis.last_parameter() {
        d = an_l.square_distance(&a_p).min(d);
    }
    Ok(d > dmax)
}

/// `V.IsSame(End1) || V.IsSame(End2)` (`cxx:556-561`), where a null
/// `TopoDS_Vertex` handle answers `false` (`IsSame` on a null handle).
fn same_as_any(v: &Option<TopoShape>, end1: Option<&Vertex>, end2: Option<&Vertex>) -> bool {
    let Some(v) = v else {
        return false;
    };
    let mut same = false;
    if let Some(e) = end1 {
        same |= crate::topo_tools::is_same(v, &e.0);
    }
    if let Some(e) = end2 {
        same |= crate::topo_tools::is_same(v, &e.0);
    }
    same
}
