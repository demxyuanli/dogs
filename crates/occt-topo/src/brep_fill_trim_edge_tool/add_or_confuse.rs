//! `BRepFill_TrimEdgeTool::AddOrConfuse` (`cxx:636-734`).
//!
//! When the arc extremity was not among the intersections answered by
//! `IntersectWith`, the first (or last) point of the bisector is projected on
//! both parallels and prepended (or appended) to the sequence.
//!
//! The `OCCT_DEBUG` printouts of `cxx:684-719` are not ported; the failure
//! branches they announce are the plain `return`s of `cxx:697`, `:705`, `:712`
//! and `:720`.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::precision::CONFUSION;
use occt_geom2d::adaptor_curve::Geom2dAdaptorCurve;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::geom2d_api_project_point_on_curve::ProjectPointOnCurve2d;

use crate::shape::Edge;
use crate::tgeometry::GeometryRegistry;

use super::{Error, TrimEdgeTool};

impl TrimEdgeTool {
    /// `AddOrConfuse(Start, Edge1, Edge2, Params)` (`cxx:636-734`).
    ///
    /// `Params` keeps the order of `IntersectWith`, so `Start` prepends and the
    /// other end appends.
    pub fn add_or_confuse(
        &self,
        start: bool,
        edge1: &Edge,
        edge2: &Edge,
        params: &mut Vec<GpPnt>,
    ) -> Result<(), Error> {
        // `constexpr double Tol = 10 * Precision::Confusion();` (`cxx:642`).
        let tol = 10.0 * CONFUSION;

        // `cxx:645-652`: the first edge's p-curve and adaptor.
        let (c1, f, l) = edge_pcurve(edge1)?;
        let ac1 = Geom2dAdaptorCurve::with_range(Arc::clone(&c1), f, l);

        // `cxx:654-662`: the bisector extremity.
        let pbis = if start {
            self.bis.value(self.bis.first_parameter())
        } else {
            self.bis.value(self.bis.last_parameter())
        }
        .ok_or(Error::NullCurve)?;

        // `cxx:664-677`: is the extremity already an intersection point?
        let mut to_proj = true;
        if !params.is_empty() {
            let u = if start {
                params[0].y()
            } else {
                params[params.len() - 1].y()
            };
            let p = ac1.value(u).ok_or(Error::NullCurve)?;
            to_proj = !pbis.is_equal(&p, tol);
        }

        if to_proj {
            // `cxx:687-690`: the second edge's p-curve.
            let (c2, f2, l2) = edge_pcurve(edge2)?;

            // `cxx:690-691`.
            let projector1 = ProjectPointOnCurve2d::with_curve_range(&pbis, &c1, f, l);
            let projector2 = ProjectPointOnCurve2d::with_curve_range(&pbis, &c2, f2, l2);

            // `cxx:693-698`.
            if projector1.nb_points() == 0 {
                return Ok(());
            }
            if !nearest_is(&projector1, &pbis, tol) {
                return Ok(());
            }
            // `cxx:707-713`.
            if projector2.nb_points() == 0 {
                return Ok(());
            }
            if !nearest_is(&projector2, &pbis, tol) {
                return Ok(());
            }

            // `cxx:721-731`.
            let u1 = projector1
                .lower_distance_parameter()
                .ok_or(Error::NullCurve)?;
            let u2 = projector2
                .lower_distance_parameter()
                .ok_or(Error::NullCurve)?;
            let mut p_int = GpPnt::new(0.0, u1, u2);
            if start {
                p_int.set_x(self.bis.first_parameter());
                params.insert(0, p_int);
            } else {
                p_int.set_x(self.bis.last_parameter());
                params.push(p_int);
            }
        }
        Ok(())
    }
}

/// `BRep_Tool::CurveOnSurface(E, C, S, L, f, l)` (`BRep_Tool.cxx:301-315`): the
/// first CurveOnSurface representation of `edge`, as the p-curve and its COS
/// range.
fn edge_pcurve(edge: &Edge) -> Result<(std::sync::Arc<dyn Curve2d>, f64, f64), Error> {
    let mut reps = GeometryRegistry::global().edge_pcurve_reps(&edge.0);
    if reps.is_empty() {
        return Err(Error::NullCurve);
    }
    let (_surf, pc, first, last) = reps.remove(0);
    Ok((pc, first, last))
}

/// `Projector.NearestPoint().IsEqual(PBis, Tol)` (`cxx:700`, `cxx:714`): a
/// projection that is not done has no nearest point and answers `false`.
fn nearest_is(
    projector: &ProjectPointOnCurve2d,
    pbis: &GpPnt2d,
    tol: f64,
) -> bool {
    projector
        .nearest_point()
        .map(|p| p.is_equal(pbis, tol))
        .unwrap_or(false)
}
