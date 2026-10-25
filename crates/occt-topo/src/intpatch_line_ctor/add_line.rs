//! `static void AddLine(L, i, j, TypeS1, TypeS2, TabIndex, slin)`
//! (`IntPatch_LineConstructor.cxx:592-880`): copy the `[i, j]` vertex range of
//! `L` into a new line of the same kind, re-numbering its vertices and
//! appending the coincident vertices OCCT duplicates through
//! `AppendSameVertex*`.
//!
//! Translation notes:
//! - The two `GeomAbs_SurfaceType` parameters (`cxx:594-597`) are unnamed and
//!   unused in OCCT; they are dropped.
//! - The `IntPatch_Line` constructors differ only in the base state, which the
//!   port copies verbatim (`IntPatchLine`), so the three-way
//!   `Undecided` / `Touch` / `In-Out` switch collapses onto one assignment.
//! - `WLine->ComputeVertexParameters(Precision::Confusion())` (`cxx:688`) and
//!   `RLine->ComputeVertexParameters(Precision::Confusion())` (`cxx:737`) are
//!   ported as `intpatch_vparams::compute_vertex_parameters_wline` /
//!   `_rline`.
//! - `wlig->SetArcOnS1(WLine->GetArcOnS1())` (`cxx:667-672`) copies the arc
//!   handles stored on the walking line itself.
//! - The cloned line carries the default creation way `IntPatch_WLUnknown`
//!   (`IntPatch_WLine.cxx:40`), which `AddLine` never overrides.

use crate::geom_int::{ALine, GLine, GLineKind, GeomIntLine, RLine};
use crate::int_tools_wline::WLine;
use crate::intpatch_vparams::{
    compute_vertex_parameters_rline, compute_vertex_parameters_wline,
};
use occt_core::precision::CONFUSION;

use super::helpers::{
    append_same_vertex_a, append_same_vertex_g, append_same_vertex_r, append_same_vertex_w,
};

/// `AddLine` (`cxx:592-880`). `i` / `j` are 1-based vertex indices,
/// `tab_index` is OCCT's 1-based `TabIndex`.
pub(crate) fn add_line(
    l: &GeomIntLine,
    i: i32,
    j: i32,
    tab_index: &mut [i32],
    slin: &mut Vec<GeomIntLine>,
) {
    let index_first_vertex = 1;
    let mut index_last_vertex = 2;
    if i == j {
        index_last_vertex = 1;
    }

    match l {
        GeomIntLine::Analytic(a) => {
            let mut alig = ALine::from_curve(a.curve.clone());
            alig.base = a.base;
            alig.add_vertex(a.vertex(i).clone());
            index_last_vertex += append_same_vertex_a(&mut alig, a, i, tab_index);
            if i != j {
                alig.add_vertex(a.vertex(j).clone());
                index_last_vertex += append_same_vertex_a(&mut alig, a, j, tab_index);
            }
            alig.set_first_point(index_first_vertex);
            alig.set_last_point(index_last_vertex);
            slin.push(GeomIntLine::Analytic(alig));
        }
        GeomIntLine::Walking(w) => {
            let mut wlig = WLine::new();
            wlig.base = w.base;
            let param_min_on_line = w.vertex(i).param_on_line as i32;
            let param_max_on_line = w.vertex(j).param_on_line as i32;
            for k in param_min_on_line..=param_max_on_line {
                if let Some(p) = w.points.get((k - 1) as usize) {
                    wlig.points.push(*p);
                }
            }
            if let Some(a) = w.arc_on_s1.clone() {
                wlig.set_arc_on_s1(a);
            }
            if let Some(a) = w.arc_on_s2.clone() {
                wlig.set_arc_on_s2(a);
            }
            let mut vtx = w.vertex(i).clone();
            vtx.param_on_line = 1.0;
            wlig.vertices.push(vtx);
            index_last_vertex += append_same_vertex_w(&mut wlig, w, i, 1.0, tab_index);
            if i != j {
                let mut vtx = w.vertex(j).clone();
                vtx.param_on_line = wlig.nb_pnts() as f64;
                wlig.vertices.push(vtx);
                let nb_pnts = wlig.nb_pnts() as f64;
                index_last_vertex += append_same_vertex_w(&mut wlig, w, j, nb_pnts, tab_index);
            }
            wlig.set_first_point(index_first_vertex);
            wlig.set_last_point(index_last_vertex);
            wlig.set_period(w.u1_period, w.v1_period, w.u2_period, w.v2_period);
            compute_vertex_parameters_wline(&mut wlig, CONFUSION);
            slin.push(GeomIntLine::Walking(wlig));
        }
        GeomIntLine::Restriction(r) => {
            // `cxx:696-698`: the two indices are re-initialised for RLine.
            let index_first_vertex = 1;
            let mut index_last_vertex = 2;
            // The RLine constructor plus `SetArcOnS1/2` (`cxx:710-718`) keep
            // the source's arc state; clearing the vertices gives the same
            // starting point.
            let mut rlig: RLine = r.clone();
            rlig.vertices.clear();
            rlig.has_first_point = false;
            rlig.has_last_point = false;

            rlig.add_vertex(r.vertex(i).clone());
            index_last_vertex += append_same_vertex_r(&mut rlig, r, i, tab_index);
            for k in (i + 1)..j {
                rlig.add_vertex(r.vertex(k).clone());
                index_last_vertex += 1;
            }
            if i != j {
                rlig.add_vertex(r.vertex(j).clone());
                index_last_vertex += append_same_vertex_r(&mut rlig, r, j, tab_index);
            }
            rlig.set_first_point(index_first_vertex);
            rlig.set_last_point(index_last_vertex);
            compute_vertex_parameters_rline(&mut rlig);
            slin.push(GeomIntLine::Restriction(rlig));
        }
        GeomIntLine::Geometric(g) => {
            let mut glig = GLine::new(g.kind.clone());
            glig.base = g.base;
            glig.add_vertex(g.vertex(i).clone());
            index_last_vertex += append_same_vertex_g(&mut glig, g, i, 0.0, tab_index);
            if i != j {
                let periodic = matches!(g.kind, GLineKind::Circ(_) | GLineKind::Elips(_));
                if periodic && i > j {
                    let mut vtx = g.vertex(j).clone();
                    vtx.param_on_line += std::f64::consts::PI * 2.0;
                    glig.add_vertex(vtx);
                    index_last_vertex += append_same_vertex_g(
                        &mut glig,
                        g,
                        j,
                        std::f64::consts::PI * 2.0,
                        tab_index,
                    );
                } else {
                    glig.add_vertex(g.vertex(j).clone());
                    index_last_vertex += append_same_vertex_g(&mut glig, g, j, 0.0, tab_index);
                }
            }
            glig.set_first_point(index_first_vertex);
            glig.set_last_point(index_last_vertex);
            slin.push(GeomIntLine::Geometric(glig));
        }
    }
}
