//! `Draft_Modification::Perform`, the surface / surface arm of the new-edge
//! block (`Draft_Modification_1.cxx:961-1331`).
//!
//! `GeomInt_IntSS` cuts the two faces' new surfaces; the arm then picks the
//! intersection line that starts at the old edge's first vertex, following the
//! reference's `Extrema_ExtPC` scan (`cxx:984-1125`) and its `detrompeur`
//! disambiguation of the "first vertex image" (`cxx:1002-1084`).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `IntSS::Perform` + the "no line" exit | `_1.cxx:961-978` | [`intersect_edge_faces`] |
//! | line selection by `Extrema_ExtPC` | `_1.cxx:984-1125` | [`intersect_edge_faces`] |
//! | chosen line, orientation, pcurves | `_1.cxx:1127-1157` | [`intersect_edge_faces`] |
//! | `Geom_BSplineCurve` gluing arm | `_1.cxx:1159-1325` | [`super::edges_spline`] |
//! | `Einf.Tolerance(max(.., TolReached3d))` | `_1.cxx:1329` | [`intersect_edge_faces`] |
//!
//! # Translation notes
//!
//! * `i2s.Line(1)->DynamicType() != STANDARD_TYPE(Geom_BSplineCurve)`
//!   (`cxx:985`) is `Curve::as_bspline_curve().is_none()`.
//! * `Extrema_ExtPC(pfv, TheCurve)` over `GeomAdaptor_Curve(i2s.Line(i))`
//!   (`cxx:990-991`) is `extrema_pc::extrema_ext_pc_range_done` over that
//!   line's own range.
//! * `Extrema_ExtPC::TrimmedSquareDistances` (`cxx:1103`) is recomputed from
//!   the two endpoint evaluations, its value when no extremum lies inside the
//!   range.
//! * The reference tests `theSurf`'s type but down-casts `S2` in the
//!   `ElSLib::Parameters` calls (`cxx:1050-1067`); the port keeps both, so the
//!   type test and the cast stay on the same handles as the reference.

use std::sync::Arc;

use occt_core::gp::GpVec2d;
use occt_core::precision::{Precision, CONFUSION};
use occt_geom::curve::Curve;
use occt_geom::Surface;
use occt_geom2d::curve::Curve2d;

use crate::brep_tool::BRepTool;
use crate::geom_int::IntSS;
use crate::shape::Edge;
use crate::topo_tools::edge_vertices;

use super::parameters::{
    surface_parameters_cylinder, surface_parameters_plane,
};
use super::DraftModification;

/// `Draft_Modification_1.cxx:961-1331`.
pub(super) fn intersect_edge_faces(
    dm: &mut DraftModification,
    the_edge: &Edge,
    c: &Arc<dyn Curve>,
    detrompeur: i32,
    ii: usize,
) -> Result<Arc<dyn Curve>, ()> {
    // `cxx:961-962`: the surfaces are re-read from the maps, so a trimmed
    // surface unwrapped by the `KPart` probe comes back as it was stored.
    let (first_face, second_face) = {
        let einf = dm.my_e_map.find_from_index(ii);
        (einf.first_face().cloned(), einf.second_face().cloned())
    };
    let (Some(first_face), Some(second_face)) = (first_face, second_face) else {
        return dm.fail_edge(the_edge);
    };
    let s1 = dm
        .my_f_map
        .find_from_key(&first_face.0)
        .and_then(|fi| fi.geometry())
        .cloned();
    let s2 = dm
        .my_f_map
        .find_from_key(&second_face.0)
        .and_then(|fi| fi.geometry())
        .cloned();
    let (Some(s1), Some(s2)) = (s1, s2) else {
        return dm.fail_edge(the_edge);
    };

    // `cxx:970-977`
    let mut i2s = IntSS::new();
    i2s.perform(s1.clone(), s2.clone(), CONFUSION, true, false, false);
    if !i2s.is_done() || i2s.nb_lines() <= 0 {
        return dm.fail_edge(the_edge);
    }

    // `cxx:853-862` are still in scope in the reference; recomputed here.
    let (Some(fv), Some(lv)) = edge_vertices(the_edge) else {
        return dm.fail_edge(the_edge);
    };
    let prmfv = BRepTool::parameter_of_vertex(&the_edge.0, &fv.0);
    let (pfv, d1fv) = c.d1(prmfv);
    // `cxx:858-861`: the last-vertex point, only read by the `detrompeur`
    // disambiguation below.
    let prmlv = BRepTool::parameter_of_vertex(&the_edge.0, &lv.0);
    let (plv, _d1lv) = c.d1(prmlv);

    // `cxx:979-983`
    let mut glob2_min = f64::MAX; // `RealLast()`
    let mut imin = 0i32;
    let mut pmin = 0.0f64;

    let first_is_bspline = i2s
        .line(1)
        .map(|l| l.curve.as_bspline_curve().is_some())
        .unwrap_or(false);

    let new_c: Option<Arc<dyn Curve>>;

    if !first_is_bspline {
        // `cxx:986-1157`
        let mut dist2_min = f64::MAX; // declared before the line loop
        let mut s1 = s1.clone();
        let mut s2 = s2.clone();
        let mut detrompeur = detrompeur;

        for i in 1..=i2s.nb_lines() {
            let Some(line) = i2s.line(i) else {
                continue;
            };
            let the_curve = line.curve.clone();
            let (sols, is_done) = occt_geom::extrema_pc::extrema_ext_pc_range_done(
                the_curve.as_ref(),
                &pfv,
                the_curve.first_parameter(),
                the_curve.last_parameter(),
            );
            let mut locpmin = 0.0f64;

            if is_done {
                // `cxx:995-1001`
                if !sols.is_empty() {
                    dist2_min = sols[0].sq_dist;
                    locpmin = sols[0].u;
                }
                // `cxx:1002-1084`
                if sols.len() == 2 && dist2_min > Precision::SQUARE_CONFUSION {
                    let d1_2 = sols[0].sq_dist;
                    let d2_2 = sols[1].sq_dist;
                    if d1_2 > 1.21 * d2_2 {
                        dist2_min = sols[1].sq_dist;
                        locpmin = sols[1].u;
                    } else if d2_2 > 1.21 * d1_2 {
                        dist2_min = sols[0].sq_dist;
                        locpmin = sols[0].u;
                    } else {
                        let pfvpar = sols[0].u;
                        let plvpar = sols[1].u;
                        let pfvprim = the_curve.d0(pfvpar);
                        let plvprim = the_curve.d0(plvpar);

                        // `cxx:1027-1046`
                        let the_surf: Arc<dyn Surface> = if detrompeur == 1 {
                            if let Some(basis) = s1.rectangular_trimmed_basis() {
                                s1 = basis;
                            }
                            s1.clone()
                        } else {
                            if let Some(basis) = s2.rectangular_trimmed_basis() {
                                s2 = basis;
                            }
                            s2.clone()
                        };

                        // `cxx:1047-1070`
                        if detrompeur != 0 && detrompeur != 4 {
                            let (mut ul, mut vl, mut uf, mut vf) = (0.0, 0.0, 0.0, 0.0);
                            let (mut ufprim, mut ulprim, mut vfprim, mut vlprim) =
                                (0.0, 0.0, 0.0, 0.0);
                            if the_surf.gp_pln().is_some() {
                                let pl = s2.gp_pln().expect("Geom_Plane");
                                (ul, vl) = surface_parameters_plane(&pl, &plv);
                                (uf, vf) = surface_parameters_plane(&pl, &pfv);
                                (ulprim, vlprim) = surface_parameters_plane(&pl, &plvprim);
                                (ufprim, vfprim) = surface_parameters_plane(&pl, &pfvprim);
                            } else if the_surf.gp_cylinder().is_some() {
                                let cy = s2.gp_cylinder().expect("Geom_CylindricalSurface");
                                (ul, vl) = surface_parameters_cylinder(&cy, &plv);
                                (uf, vf) = surface_parameters_cylinder(&cy, &pfv);
                                (ulprim, vlprim) = surface_parameters_cylinder(&cy, &plvprim);
                                (ufprim, vfprim) = surface_parameters_cylinder(&cy, &pfvprim);
                            } else {
                                detrompeur = 4;
                            }

                            // `cxx:1072-1082`
                            if detrompeur == 1 || detrompeur == 2 {
                                let v1 = GpVec2d::new(ul - ufprim, vl - vfprim);
                                let norm = GpVec2d::new(vf - vfprim, ufprim - uf);
                                let v2 = GpVec2d::new(ulprim - ufprim, vlprim - vfprim);
                                if v1.dot(&norm) * v2.dot(&norm) < 0.0 {
                                    dist2_min = sols[1].sq_dist;
                                    locpmin = sols[1].u;
                                }
                            }
                        }
                    }
                }

                // `cxx:1086-1121`
                if sols.len() == 1 || sols.len() > 2 || detrompeur == 4 {
                    dist2_min = sols[0].sq_dist;
                    locpmin = sols[0].u;
                    for j in 2..=sols.len() {
                        let dist2 = sols[j - 1].sq_dist;
                        if dist2 < dist2_min {
                            dist2_min = dist2;
                            locpmin = sols[j - 1].u;
                        }
                    }
                } else if sols.is_empty() {
                    // `cxx:1099-1119`: `TrimmedSquareDistances`.
                    let first = the_curve.first_parameter();
                    let last = the_curve.last_parameter();
                    let dist1_2 = pfv.square_distance(&the_curve.d0(first));
                    let dist2_2 = pfv.square_distance(&the_curve.d0(last));
                    if dist1_2 < dist2_2 {
                        dist2_min = dist1_2;
                        locpmin = first;
                    } else {
                        dist2_min = dist2_2;
                        locpmin = last;
                    }
                }

                // `cxx:1123-1127`
                if dist2_min < glob2_min {
                    glob2_min = dist2_min;
                    imin = i;
                    pmin = locpmin;
                }
            }
        }

        // `cxx:1129-1132`
        if imin == 0 {
            return dm.fail_edge(the_edge);
        }

        // `cxx:1134-1140`
        let Some(mut chosen) = i2s.line(imin).map(|l| l.curve.clone()) else {
            return dm.fail_edge(the_edge);
        };
        let (_, newd1) = chosen.d1(pmin);
        let ya_rev = d1fv.dot(&newd1) < 0.0;
        if ya_rev {
            chosen = Arc::from(chosen.reversed());
        }

        // `cxx:1142-1156`
        let first_pc = i2s.line_on_s1(imin).map(|p| {
            if ya_rev {
                Arc::from(p.reversed()) as Arc<dyn Curve2d>
            } else {
                p.clone()
            }
        });
        let second_pc = i2s.line_on_s2(imin).map(|p| {
            if ya_rev {
                Arc::from(p.reversed()) as Arc<dyn Curve2d>
            } else {
                p.clone()
            }
        });
        {
            let ei = dm.my_e_map.change_from_index(ii);
            if let Some(pc) = first_pc {
                *ei.change_first_pc() = Some(pc);
            }
            if let Some(pc) = second_pc {
                *ei.change_second_pc() = Some(pc);
            }
        }
        new_c = Some(chosen);
    } else {
        // `cxx:1159-1325`
        new_c = super::edges_spline::glue_bspline_lines(&i2s, &pfv, &d1fv, &s1);
    }
    // `cxx:1329`
    let tol3d = i2s.tol_reached_3d();
    {
        let ei = dm.my_e_map.change_from_index(ii);
        let tol = ei.tolerance().max(tol3d);
        ei.set_tolerance(tol);
    }

    match new_c {
        Some(nc) => Ok(nc),
        None => dm.fail_edge(the_edge),
    }
}
