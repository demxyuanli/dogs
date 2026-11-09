//! `BRepOffset::CollapseSingularities` (`BRepOffset.cxx:206-365`).
//!
//! Only `Geom_BSplineSurface` is processed (`cxx:211-218`); any other surface
//! is answered unchanged. For a B-spline surface the four sides of the pole
//! grid (`U=0`, `V=0`, `U=1`, `V=1`, `cxx:249-257`) are tested: a side whose
//! poles all coincide with their centre of gravity is already collapsed
//! (`cxx:308-311`). A side that is *singular* -- every pole within
//! `thePrecision` of that centre -- or that fits exactly one vertex of a
//! degenerated edge (`cxx:313-340`) has all its poles set to the centre
//! (`cxx:342-355`). A degenerated point that fits such a side is dropped from
//! the candidate list so a single point cannot collapse two sides (`cxx:322-335`).
//!
//! Port notes:
//!
//! * `aPoles.LowerRow()` / `LowerCol()` are always 1 in OCCT's `Array2`
//!   accessors here, and `RowLength()` / `ColLength()` are
//!   `nb_poles_u()` / `nb_poles_v()` of the port's pole grid.
//! * `NCollection_List::Remove(anIterator)` (`cxx:326-327`) is the rebuild of
//!   the candidate vectors in [`remove_fitted`], which drops the same entries in
//!   the same order.
//! * `theSurface->Copy()` (`cxx:347`) is [`Surface::osculating_bspline`], which
//!   answers a clone for a `Geom_BSplineSurface`.
//! * `gp::Resolution()` is `RealSmall()` (`DBL_MIN`), never `Precision::Angular`;
//!   the port uses [`REAL_SMALL`].

use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_core::precision::REAL_SMALL;
use occt_geom::surface::Surface;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::shape::Face;
use crate::topo_tools::{edge_vertices, explorer_shapes};

/// `TopExp::Vertices(anEdge, aV1, aV2)` followed by the
/// `!aV1.IsSame(aV2)` test of `cxx:230-235`: the singular vertex of a
/// degenerated edge, with its tolerance.
fn degenerated_vertex(edge: &crate::shape::Edge) -> Option<(GpPnt, f64)> {
    // `if (!BRep_Tool::Degenerated(anEdge)) continue;` (`cxx:226-229`).
    if !BRepTool::is_degenerated(edge) {
        return None;
    }
    let (v1, v2) = edge_vertices(edge);
    let (v1, v2) = (v1?, v2?);
    // `if (!aV1.IsSame(aV2)) continue;` (`cxx:232-235`).
    if !v1.0.same_tshape(&v2.0) {
        return None;
    }
    // `aDegenPnt.Append(BRep_Tool::Pnt(aV1));
    //  aDegenTol.Append(BRep_Tool::Tolerance(aV1));` (`cxx:237-238`).
    Some((BRepTool::vertex_point(&v1), BRepTool::vertex_tolerance(&v1)))
}

/// `aDegenPnt.Remove(aDegPntIt)` / `aDegenTol.Remove(aDegTolIt)` (`cxx:326-327`)
/// for the entries whose `isDegenerated` flag is set: the kept entries keep
/// their relative order, so the two vectors stay aligned.
fn remove_fitted(
    degen: &mut Vec<(GpPnt, f64)>,
    is_degenerated: &[bool],
) -> usize {
    let mut a_nb_fit = 0;
    let mut kept: Vec<(GpPnt, f64)> = Vec::with_capacity(degen.len());
    for (i_degen, &deg) in is_degenerated.iter().enumerate() {
        if deg {
            a_nb_fit += 1;
        } else {
            kept.push(degen[i_degen].clone());
        }
    }
    *degen = kept;
    a_nb_fit
}

/// `BRepOffset::CollapseSingularities(theSurface, theFace, thePrecision)`
/// (`BRepOffset.cxx:206-365`).
pub fn collapse_singularities(
    the_surface: &Arc<dyn Surface>,
    the_face: &Face,
    the_precision: f64,
) -> Arc<dyn Surface> {
    // `occ::handle<Standard_Type> aType = theSurface->DynamicType();
    //  if (aType != STANDARD_TYPE(Geom_BSplineSurface)) return theSurface;`
    // (`cxx:212-218`): only B-spline surfaces are treated, so a trimmed
    // B-spline (its own type) is answered unchanged.
    if !the_surface.is_bspline_surface() {
        return the_surface.clone();
    }
    let Some(bsp) = the_surface.osculating_bspline() else {
        return the_surface.clone();
    };

    // `cxx:221-239`: the singularities, i.e. the vertices of degenerated edges.
    let mut a_degen: Vec<(GpPnt, f64)> = Vec::new();
    for edge_shape in explorer_shapes(&the_face.0, ShapeType::Edge) {
        let Some(edge) = crate::shape::Edge::wrap(edge_shape) else {
            continue;
        };
        if let Some(d) = degenerated_vertex(&edge) {
            a_degen.push(d);
        }
    }

    // `cxx:245`: `const NCollection_Array2<gp_Pnt>& aPoles = aBSpline->Poles();`
    let nb_rows = bsp.nb_poles_u();
    let nb_cols = bsp.nb_poles_v();
    if nb_rows == 0 || nb_cols == 0 {
        return the_surface.clone();
    }

    // `cxx:250-257`: the four sides, {U=0; V=0; U=1; V=1}.
    let row_start = [0usize, 0, nb_rows - 1, 0];
    let col_start = [0usize, 0, 0, nb_cols - 1];
    let row_step: [usize; 4] = [0, 1, 0, 1];
    let col_step: [usize; 4] = [1, 0, 1, 0];
    let nb_steps = [nb_cols, nb_rows, nb_cols, nb_rows];

    // `cxx:247`: set as soon as a side is collapsed.
    let mut a_copy: Option<occt_geom::bspline_surface::GeomBSplineSurface> = None;

    for i_side in 0..4 {
        // `cxx:260-268`: the centre of gravity of the side poles.
        let mut a_sum = [0.0f64; 3];
        for i_pole in 0..nb_steps[i_side] {
            let p = bsp.pole(
                row_start[i_side] + i_pole * row_step[i_side],
                col_start[i_side] + i_pole * col_step[i_side],
            );
            a_sum[0] += p.x();
            a_sum[1] += p.y();
            a_sum[2] += p.z();
        }
        let n = nb_steps[i_side] as f64;
        let a_center = GpPnt::from_xyz(&occt_core::gp::GpXyz::new(a_sum[0] / n, a_sum[1] / n, a_sum[2] / n));

        // `cxx:270-279`: the two flags and the per-degenerated-point fit array.
        let mut is_collapsed = true;
        let mut is_singular = true;
        let mut is_degenerated = vec![true; a_degen.len()];

        // `cxx:280-307`.
        for i_pole in 0..nb_steps[i_side] {
            let a_pole = bsp.pole(
                row_start[i_side] + i_pole * row_step[i_side],
                col_start[i_side] + i_pole * col_step[i_side],
            );
            // `double aDistCG = aCenter.Distance(aPole);` (`cxx:286`).
            let a_dist_cg = a_center.distance(&a_pole);
            if a_dist_cg > REAL_SMALL {
                is_collapsed = false;
            }
            if a_dist_cg > the_precision {
                is_singular = false;
            }
            // `cxx:299-306`: the distances to the degenerated points.
            for (i_degen, d) in a_degen.iter().enumerate() {
                if is_degenerated[i_degen] && d.0.distance(&a_pole) >= d.1 {
                    is_degenerated[i_degen] = false;
                }
            }
        }

        // `cxx:308-311`.
        if is_collapsed {
            continue;
        }

        // `cxx:315-340`.
        if !is_singular {
            let a_nb_fit = remove_fitted(&mut a_degen, &is_degenerated);
            // If the side fits more than one degenerated vertex, do not
            // collapse it, to be on the safe side (`cxx:337-339`).
            is_singular = a_nb_fit == 1;
        }

        // `cxx:342-355`: do collapse.
        if is_singular {
            if a_copy.is_none() {
                a_copy = Some(bsp.clone());
            }
            let copy = a_copy.as_mut().expect("copy");
            for i_pole in 0..nb_steps[i_side] {
                copy.set_pole(
                    row_start[i_side] + i_pole * row_step[i_side],
                    col_start[i_side] + i_pole * col_step[i_side],
                    a_center.clone(),
                );
            }
        }
    }

    // `cxx:358-361`.
    match a_copy {
        Some(copy) => Arc::new(copy),
        None => the_surface.clone(),
    }
}
