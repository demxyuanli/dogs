//! ImpPrm `Perform` tail: restriction vertices on the implicit surface.
//! Source: `IntPatch_ImpPrmIntersection.cxx:1766-1806` — `PutVertexOnLine` on
//! every section line, dead-line removal, walking lines moved to the end.

use occt_core::precision::CONFUSION;
use occt_geom::Surface;

use crate::geom_int::{GeomIntLine, TopolTool};
use crate::intpatch_rst_int::put_vertex_on_line;

/// `PutVertexOnLine` over every line, then drop the lines whose point count
/// falls to two or less and re-append the walking lines (`cxx:1766-1806`).
pub(crate) fn place_restriction_vertices(
    slin: &mut Vec<GeomIntLine>,
    s1: &dyn Surface,
    d1: &TopolTool,
    s2: &dyn Surface,
    d2: &TopolTool,
    reversed: bool,
    tol_tang: f64,
) {
    let mut i = 0usize;
    while i < slin.len() {
        // `occ::down_cast<IntPatch_PointLine>(slin(i))` (`cxx:1768`): only the
        // walking and restriction lines are point lines.
        if !matches!(
            slin[i],
            GeomIntLine::Walking(_) | GeomIntLine::Restriction(_)
        ) {
            i += 1;
            continue;
        }
        let mut line = slin.remove(i);
        // `cxx:1771-1783` — the domain belongs to `Surf1` (`!reversed`) or to
        // `Surf2` (`reversed`).
        let result = if !reversed {
            let mut dom = d1.brep_tool(s1);
            put_vertex_on_line(&mut line, s1, &mut dom, s2, true, tol_tang)
        } else {
            let mut dom = d2.brep_tool(s2);
            put_vertex_on_line(&mut line, s2, &mut dom, s1, false, tol_tang)
        };
        // OCCT ignores the outcome: `PutVertexOnLine` is `void` and raises for
        // a non-point line, which cannot reach this loop.
        let _ = result;

        // `cxx:1785-1800`: `NbPnts() < 2`, or two coincident points.
        let nb_pnts = match &line {
            GeomIntLine::Walking(w) => w.nb_pnts(),
            GeomIntLine::Restriction(r) => r.points.len() as i32,
            _ => 0,
        };
        if nb_pnts < 2 {
            continue;
        }
        if nb_pnts == 2 && first_two_same(&line) {
            continue;
        }

        if matches!(line, GeomIntLine::Walking(_)) {
            // `slin.Append(aWL); slin.Remove(i);` (`cxx:1802-1805`).
            slin.push(line);
        } else {
            slin.insert(i, line);
            i += 1;
        }
    }
}

/// `aL->Point(1).IsSame(aL->Point(2), Precision::Confusion())`: with the
/// default `theTol2D = -1` only the 3D distance is compared
/// (`IntSurf_PntOn2S.cxx:85-97`).
fn first_two_same(line: &GeomIntLine) -> bool {
    match line {
        GeomIntLine::Walking(w) => w.points[0].is_same(&w.points[1], CONFUSION, -1.0),
        GeomIntLine::Restriction(r) => r.points[0].is_same(&r.points[1], CONFUSION, -1.0),
        _ => false,
    }
}
