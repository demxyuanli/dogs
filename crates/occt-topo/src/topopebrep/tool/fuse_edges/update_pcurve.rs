//! `UpdatePCurve` (`TopOpeBRepTool_FuseEdges.cxx:885-1038`).
//!
//! Rebuilds, on the fused edge, every p-curve the chain's edges carried on the
//! surfaces of the faces they were shared by.
//!
//! The loop is driven by the **index** overload of
//! `BRep_Tool::CurveOnSurface(E, C2d, S, loc, f, l, iedg)`
//! (`BRep_Tool.cxx:488-530`): it walks the edge's `CurveOnSurface`
//! representations, counting a representation on a closed surface as two
//! entries (`PCurve` then `PCurve2`, `cxx:512-520`), and stops when the index
//! runs past the last one. [`GeometryRegistry::edge_pcurve_reps`] answers the
//! same sequence (forward p-curve then reversed for a seam), so iterating it
//! is the same walk.
//!
//! A representation is only rebuilt when one of the old edge's faces sits on
//! the **same** `Geom_Surface` handle (`cxx:910-924`, handle identity, not
//! `IsSame`). For a seam the two p-curves are written back together
//! (`cxx:939-956`); otherwise the p-curve is
//!
//! 1. unwrapped one `Geom2d_TrimmedCurve` level (`cxx:958-961`),
//! 2. concatenated with the p-curves of the other chain edges when its span is
//!    shorter than the new 3D range (`cxx:963-1006`) — that is
//!    `Geom2dConvert_CompCurveToBSplineCurve` with the chain's p-curves added
//!    after it, followed by a `BSplCLib::Reparametrize` onto `[ef, el]` when
//!    the concatenation did not land there,
//! 3. written back with `B.UpdateEdge` (`cxx:1007`).
//!
//! Finally the COS range is adapted to the old edge's orientation
//! (`cxx:1015-1023`) and, when a p-curve was rebuilt, the new edge is made
//! same-parameter again (`cxx:1029-1033`).
//!
//! `false` is the C++ `false` of `cxx:987` ("cannot merge pcurves"), which
//! leaves that chain unfused in `BuildListResultEdges` (`cxx:373-377`).
//! The `Standard_ConstructionError` of `Concat.Add` (`cxx:987-990` is the
//! `return false`; the raise lives inside `Geom2dConvert`) also answers
//! `false` here.

use std::sync::Arc;

use occt_core::precision::PCONFUSION;
use occt_core::convert::ParameterisationType;
use occt_geom2d::comp_curve_to_bspline::CompCurveToBSplineCurve2d;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::brep_tool::BRepTool;
use crate::shape::{Edge, TopoShape};
use crate::tgeometry::GeometryRegistry;

use super::FuseEdges;

/// `Geom2d_BoundedCurve` down-cast (`cxx:967`, `cxx:979`): `Geom2d_Line`,
/// `Circle`, `Ellipse`, `Parabola` and `Hyperbola` are the unbounded types, so
/// everything else is treated as bounded.
fn is_bounded_2d(c: &dyn Curve2d) -> bool {
    c.is_bspline2d() || c.is_bezier2d() || c.trimmed_basis().is_some()
}

/// `C = down_cast<Geom2d_TrimmedCurve>(C)->BasisCurve()` (`cxx:947-951`,
/// `cxx:958-961`): one level, not a loop.
fn untrim_2d(c: Arc<dyn Curve2d>) -> Arc<dyn Curve2d> {
    match c.trimmed_basis() {
        Some(b) => Arc::from(b.clone_dyn()),
        None => c,
    }
}

/// `down_cast<Geom2d_BoundedCurve>(C)` or `new Geom2d_TrimmedCurve(C, f, l)`
/// (`cxx:966-971`, `cxx:979-984`).
fn as_bounded_2d(c: Arc<dyn Curve2d>, first: f64, last: f64) -> Arc<dyn Curve2d> {
    if is_bounded_2d(c.as_ref()) {
        c
    } else {
        Arc::new(Geom2dTrimmedCurve::new(c, first, last))
    }
}

impl FuseEdges {
    /// `UpdatePCurve(theOldEdge, theNewEdge, theLstEdg)` (`cxx:885-1038`).
    pub(super) fn update_pcurve(
        &self,
        the_old_edge: &Edge,
        the_new_edge: &mut Edge,
        the_lst_edg: &[TopoShape],
    ) -> bool {
        let reg = GeometryRegistry::global();
        let old_shape = the_old_edge.0.clone();
        // `cxx:900`: `myMapEdgLstFac.FindFromKey(theOldEdge)`.
        let lmap_fac: Vec<TopoShape> = self
            .my_map_edg_lst_fac
            .find_from_key(&old_shape)
            .map(|v| v.to_vec())
            .unwrap_or_default();

        let (ef, el) = BRepTool::edge_parameters(the_new_edge);
        let mut pcurve_rebuilt = false;

        // `cxx:903`: the `CurveOnSurface(..., iedg)` sequence.
        for (surf, _pc, cf, cl) in reg.edge_pcurve_reps(&old_shape) {
            let surf_ptr = Arc::as_ptr(&surf) as *const () as usize;
            // `cxx:910-924`: a face of the old edge on that same surface handle.
            let mut same_surf: Option<TopoShape> = None;
            for face in &lmap_fac {
                let Some(fs) = reg.face_surface(face) else {
                    continue;
                };
                if Arc::as_ptr(&fs) as *const () as usize == surf_ptr {
                    same_surf = Some(face.clone());
                    break;
                }
            }
            let Some(face) = same_surf else {
                continue;
            };
            let face_key = GeometryRegistry::shape_key(&face);

            // The p-curves of the old edge on that surface, forward then
            // reversed for a seam (`BRep_TEdge::myCurves` order).
            let pcs = reg.edge_pcurves(&old_shape, face_key);
            let Some(first_pc) = pcs.first().cloned() else {
                continue;
            };

            if pcs.len() >= 2 {
                // `cxx:939-956`: closed surface. `aFEdge.Reverse()` +
                // `CurveOnSurface(aFEdge, aFFace, cf, cl)` reads `PCurve2`.
                let c1 = untrim_2d(first_pc);
                let c2 = untrim_2d(pcs[1].clone());
                reg.set_edge_pcurves(&the_new_edge.0, face_key, vec![c1, c2]);
            } else {
                // `cxx:958-961`
                let mut curv2d = untrim_2d(first_pc);
                let (f_c2d, l_c2d) = (curv2d.first_parameter(), curv2d.last_parameter());
                // `cxx:966`: `l - f + 2. * Epsilon(l - f) < el - ef`.
                if l_c2d - f_c2d + 2.0 * f64::EPSILON * (l_c2d - f_c2d) < el - ef {
                    // `cxx:968-971`: `Concat(bcurve)`.
                    let bcurve = as_bounded_2d(curv2d.clone(), cf, cl);
                    let Some(mut concat) = CompCurveToBSplineCurve2d::from_curve(
                        bcurve.as_ref(),
                        ParameterisationType::TgtThetaOver2,
                    ) else {
                        // UNPORTED inside Geom2dConvert: the basis conversion
                        // OCCT performs in the constructor raises
                        // `Standard_ConstructionError`.
                        return false;
                    };
                    // `cxx:972-990`: add the other chain edges' p-curves.
                    for e in the_lst_edg.iter().skip(1) {
                        // `cxx:975-984`
                        let Some((c, first, last)) = pcurve_on_surface(e, face_key) else {
                            return false;
                        };
                        let bc = as_bounded_2d(c, first, last);
                        // `cxx:986-990`: `Concat.Add(BC, Precision::PConfusion())`.
                        match concat.add(bc.as_ref(), PCONFUSION, true) {
                            Ok(true) => {}
                            // OCCT raises inside `Geom2dConvert` on `Err`.
                            _ => return false,
                        }
                    }
                    // `cxx:991`: `Curv2d = Concat.BSplineCurve();`
                    let Some(mut bs) = concat.into_curve() else {
                        return false;
                    };
                    // `cxx:993-1005`: force the range onto `[ef, el]`.
                    let (first, last) = (bs.first_parameter(), bs.last_parameter());
                    if (first - ef).abs() > PCONFUSION || (last - el).abs() > PCONFUSION {
                        // `BSplCLib::Reparametrize(ef, el, Knots)` then
                        // `bc->SetKnots(Knots)`.
                        let mut knots = bs.bspline_distinct_knots_mults().map(|(k, _)| k);
                        if let Some(k) = knots.as_mut() {
                            occt_core::bspl::knots::reparametrize(ef, el, k);
                            let _ = bs.set_knots(k);
                        }
                    }
                    curv2d = Arc::new(bs);
                    pcurve_rebuilt = true;
                }
                // `cxx:1007`: `B.UpdateEdge(theNewEdge, Curv2d, Surf, loc, tol)`.
                reg.set_edge_pcurve(&the_new_edge.0, face_key, curv2d);
            }

            // `cxx:1015-1023`: adapt the edge range to the p-curve range.
            if the_old_edge.0.orientation().is_reversed() {
                reg.set_ranges_all(&the_new_edge.0, cl - el + ef, cl);
            } else {
                reg.set_ranges_all(&the_new_edge.0, cf, cf + el - ef);
            }
        }

        // `cxx:1029-1033`: force same parameter again.
        if pcurve_rebuilt {
            reg.set_same_parameter(&the_new_edge.0, false);
            crate::brep_lib_same_parameter::same_parameter_edge_inplace(
                the_new_edge,
                BRepTool::edge_tolerance(the_new_edge),
            );
        }
        true
    }
}

/// `BRep_Tool::CurveOnSurface(E, S, loc, f, l)` (`BRep_Tool.cxx:301-316`): the
/// first `CurveOnSurface` representation for the surface behind `face_key`,
/// with that representation's range.
fn pcurve_on_surface(e: &TopoShape, face_key: usize) -> Option<(Arc<dyn Curve2d>, f64, f64)> {
    let reg = GeometryRegistry::global();
    let some_pc = reg.edge_pcurves(e, face_key).into_iter().next()?;
    let (f, l) = reg
        .pcurve_range(e, face_key)
        .unwrap_or_else(|| reg.edge_parameters(e));
    Some((some_pc, f, l))
}
