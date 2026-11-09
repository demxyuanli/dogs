//! `KeepInsidePoints` (`Contap_Contour.cxx:515-554`): drop the interior points
//! that already sit on a solution restriction.

use occt_core::gp::GpPnt2d;
use occt_core::precision::CONFUSION;

use super::super::iwalking::types::InteriorPoint;
use super::super::search_inside::ContapTheSearchInside;
use super::super::surf_function::ContapSurfFunction;
use super::super::the_search::ContapTheSearch;
use super::super::tools::hcont_project;

/// `KeepInsidePoints(solins, solrst, Func, seqpins)` (`cxx:515-554`).
pub(crate) fn keep_inside_points(
    solins: &ContapTheSearchInside,
    solrst: &ContapTheSearch,
    func: &mut ContapSurfFunction,
) -> Vec<InteriorPoint> {
    let mut seqpins = Vec::new();
    let Some(surf) = func.surface().cloned() else {
        return seqpins;
    };

    let nba = solrst.nb_segments().unwrap_or(0);
    let nbp = solins.nb_points();

    for indp in 1..=nbp {
        let Some(pti) = solins.value(indp) else {
            continue;
        };
        let (u, v) = pti.parameters();
        let toproj = GpPnt2d::new(u, v);
        let mut tokeep = true;
        for inda in 1..=nba {
            let Some(seg) = solrst.segment(inda) else {
                continue;
            };
            let thearc = seg.curve();
            if let Some((_paramproj, pproj)) = hcont_project(thearc, &toproj) {
                // `cxx:542`: the 3D point of the projected parameters.
                let pprojete = surf.value(pproj.x(), pproj.y()).unwrap_or_default();
                if pti.value().distance(&pprojete) <= CONFUSION {
                    tokeep = false;
                    break;
                }
            }
        }
        if tokeep {
            seqpins.push(pti.clone());
        }
    }

    seqpins
}
