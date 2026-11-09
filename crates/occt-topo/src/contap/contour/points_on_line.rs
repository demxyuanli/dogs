//! `FindLine` (`Contap_Contour.cxx:1973-2011`) and `PutPointsOnLine`
//! (`cxx:2013-2098`): attach the boundary start points found by
//! `Contap_TheSearch` to the analytical line or circle they lie on.

use std::sync::Arc;

use occt_core::elib::clib;
use occt_core::gp::{GpDir, GpPnt, GpVec};
use occt_core::precision::{CONFUSION, REAL_EPSILON};

use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::intsurf_transition::{make_transition, IntSurfTransition};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::super::line::ContapLine;
use super::super::point::ContapPoint;
use super::super::surf_props;
use super::super::the_search::ContapTheSearch;
use super::super::types::ContapIType;
use super::super::tools::hcurve2d_d1;
use super::TOLPETIT;

/// The outputs of `FindLine` (`cxx:1973-2011`) other than its result.
struct FindLineResult {
    /// `Ptref`, written by `Contap_SurfProps::Normale`.
    ptref: GpPnt,
    /// `Paramin`.
    paramin: f64,
    /// `Tgmin`.
    tgmin: GpVec,
    /// `Norm`.
    norm: GpVec,
    /// `Ptmin.SquareDistance(Ptref) <= Tolpetit`.
    goon: bool,
}

/// `FindLine(Line, Surf, Pt2d, Ptref, Paramin, Tgmin, Norm)` (`cxx:1973-2011`).
///
/// `Dismin` starts at `RealLast()` and every candidate overwrites it, so the
/// "nearest" test of `cxx:2004-2009` always selects the only candidate; the
/// minimum selection is kept as a direct assignment.
fn find_line(line: &ContapLine, surf: &SurfaceAdaptor, pt2d: (f64, f64)) -> Option<FindLineResult> {
    let (ptref, norm) = surf_props::normale(surf, pt2d.0, pt2d.1)?;

    let (paramin, ptmin, tgmin) = match line.type_contour() {
        ContapIType::Lin => {
            let lin = line.line().ok()?;
            let para = clib::line_parameter(&lin, &ptref);
            let (pt, tg) = clib::line_d1(&lin, para);
            let _dist = pt.distance(&ptref) + norm.dot(&GpVec::from_xyz(lin.direction().xyz())).abs();
            (para, pt, tg)
        }
        ContapIType::Circle => {
            let cir = line.circle().ok()?;
            let pos = cir.position();
            let para = clib::circle_parameter(&pos, &ptref);
            let (pt, tg) = clib::circle_d1(&cir, para);
            // `Norm.Dot(tg / cir.Radius())`, `|tg| = Radius`.
            let _dist = pt.distance(&ptref) + (norm.dot(&tg) / cir.radius()).abs();
            (para, pt, tg)
        }
        // `cxx:1998` falls into the `Contap_Circle` arm for every non-line and
        // `Contap_Line::Circle()` raises `Standard_DomainError` for a line that
        // is neither a straight line nor a circle. `PutPointsOnLine` is only
        // entered with the analytical arm's lines and circles, so the port
        // declines instead of raising.
        _ => return None,
    };

    Some(FindLineResult {
        ptref,
        paramin,
        tgmin,
        norm,
        goon: ptmin.square_distance(&ptref) <= TOLPETIT,
    })
}

/// `PutPointsOnLine(solrst, Surf, slin)` (`cxx:2013-2098`).
pub(crate) fn put_points_on_line(
    solrst: &ContapTheSearch,
    surf: &SurfaceAdaptor,
    slin: &mut [ContapLine],
) {
    let nb_points = solrst.nb_points().unwrap_or(0);
    let nb_lin = slin.len();

    for l in 0..nb_lin {
        for i in 1..=nb_points {
            let Some(pstart) = solrst.point(i as usize) else {
                continue;
            };
            let thearc: Arc<Curve2dAdaptor> = pstart.arc().clone();
            let theparam = pstart.parameter();

            let (pt2d, d2d) = hcurve2d_d1(&thearc, theparam);

            let Some(found) = find_line(&slin[l], surf, (pt2d.x(), pt2d.y())) else {
                continue;
            };
            if !found.goon {
                continue;
            }
            let ptonsurf = found.ptref;
            let paramlin = found.paramin;
            let vectg = found.tgmin;
            let normale = found.norm;

            // `Adaptor3d_HSurfaceTool::D1(Surf, u, v, bidpt, d1u, d1v)`.
            let Some((_bidpt, d1u, d1v, _nrm)) = surf_props::deriv_and_norm(surf, pt2d.x(), pt2d.y())
            else {
                continue;
            };

            let mut ppoint = ContapPoint::with_uv(ptonsurf, pt2d.x(), pt2d.y());
            let mut t_line = IntSurfTransition::new();
            let mut t_arc = IntSurfTransition::new();

            if normale.magnitude() < REAL_EPSILON {
                t_line.set_undecided();
                t_arc.set_undecided();
            } else {
                // Small test for cone tips and sphere poles (`cxx:2062-2078`).
                let mut tgtrst;
                if d2d.y().abs() <= CONFUSION {
                    tgtrst = d1v.crossed(&normale);
                    if d2d.x() < 0.0 {
                        // `tgtrst.Reverse()`
                        tgtrst = tgtrst.reversed();
                    }
                } else {
                    tgtrst = GpVec::new(
                        d2d.x() * d1u.x() + d2d.y() * d1v.x(),
                        d2d.x() * d1u.y() + d2d.y() * d1v.y(),
                        d2d.x() * d1u.z() + d2d.y() * d1v.z(),
                    );
                }
                let ndir = GpDir::from_vec(&normale)
                    .unwrap_or_else(|_| GpDir::new(1.0, 0.0, 0.0).expect("+X"));
                make_transition(&vectg, &tgtrst, &ndir, &mut t_line, &mut t_arc);
            }

            ppoint.set_arc(thearc.clone(), theparam, t_line, t_arc);
            ppoint.set_parameter(paramlin);
            if !pstart.is_new() {
                if let Some(v) = pstart.vertex() {
                    ppoint.set_vertex(v.clone());
                }
            }
            slin[l].add_vertex(ppoint);
        }
    }
}
