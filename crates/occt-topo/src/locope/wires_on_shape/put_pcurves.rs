//! File-local `PutPCurves(Efrom, Eto, myShape)` (`LocOpe_WiresOnShape.cxx:913-1288`).
//!
//! Two cases:
//!
//! * `Efrom` is shared by two distinct faces: `Efrom` and `Eto` are different
//!   edges glued along the same face boundary, so a fresh p-curve of `Efrom` is
//!   built on each face that does not have one yet (`cxx:950-1097`).
//! * `Efrom` and `Eto` are the two representations of a closed (seam) edge on a
//!   single face: the existing p-curves are re-read on the forward and reversed
//!   edge and both are stored, after the seam-parameter translation
//!   (`cxx:1099-1285`).

use std::sync::Arc;

use occt_core::gp::{GpTrsf2d, GpVec2d};
use occt_core::precision::{self, PCONFUSION};
use occt_geom2d::curve::Curve2d;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d;
use crate::brep_offset::update_edge_pcurves;
use crate::brep_tool::BRepTool;
use crate::brep_tools_uv::uv_bounds_face;
use crate::brep_tools_wire_explorer::edge_vertices_cum_ori;
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{edges_of, explorer_shapes, is_same};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::projection::project_vertex_2d_on_edge;

/// `gp_Geom2d_Curve::Translated(gp_Vec2d)` (`cxx:1237`, `cxx:1242`, ...).
fn translated(curve: &Arc<dyn Curve2d>, dx: f64, dy: f64) -> Arc<dyn Curve2d> {
    let mut t = GpTrsf2d::identity();
    t.set_translation_vec(&GpVec2d::new(dx, dy));
    Arc::from(curve.transformed(&t))
}

/// `BRep_Tool::CurveOnSurface(E.Oriented(theOrient), F.Oriented(FORWARD), f, l)`
/// (`cxx:1108-1141`): the seam representation selected by the edge orientation.
fn curve_on_surface_oriented(
    edge: &Edge,
    face: &Face,
    orient: Orientation,
) -> Option<(Arc<dyn Curve2d>, f64, f64)> {
    let e = Edge(edge.0.clone().oriented(orient));
    let f = Face(face.0.clone().oriented(Orientation::Forward));
    boptools_2d::curve_on_surface_range(&e, &f)
}

/// `PutPCurves(Efrom, Eto, myShape)` (`LocOpe_WiresOnShape.cxx:913-1288`).
///
/// Answers `Err` where the reference throws `Standard_ConstructionError`
/// (`cxx:932`, `cxx:1105`).
pub(super) fn put_pcurves(
    e_from: &Edge,
    e_to: &Edge,
    my_shape: &TopoShape,
) -> Result<(), String> {
    // `cxx:917-928`: every face of the shape that carries `Eto` as an edge.
    let mut l_faces: Vec<Face> = Vec::new();
    for face_shape in explorer_shapes(my_shape, ShapeType::Face) {
        let Some(face) = Face::wrap(face_shape) else {
            continue;
        };
        for e in edges_of(&face.0) {
            if is_same(&e.0, &e_to.0) {
                l_faces.push(face.clone());
                break;
            }
        }
    }

    // `cxx:930-934`.
    if l_faces.len() != 1 && l_faces.len() != 2 {
        return Err(
            "LocOpe_WiresOnShape::PutPCurves: Standard_ConstructionError".to_string(),
        );
    }
    // `cxx:937-940`: a free boundary, left to `PutPCurve`.
    if l_faces.len() == 1 {
        return Ok(());
    }

    // `cxx:950`: two different faces -> project `Efrom` on each of them.
    if !is_same(&l_faces[0].0, &l_faces[1].0) {
        for fac in &l_faces {
            // `cxx:958-961`.
            if boptools_2d::curve_on_surface_range(e_from, fac).is_some() {
                continue;
            }
            // `cxx:962-971`: the basis surface, and a plane is skipped.
            let Some(mut s) = BRepTool::face_surface_world(fac) else {
                continue;
            };
            if let Some(basis) = s.rectangular_trimmed_basis() {
                s = basis;
            }
            if s.gp_pln().is_some() {
                continue;
            }

            // `cxx:972`.
            let (umin, umax, vmin, vmax) = uv_bounds_face(fac);

            // `cxx:973-980`.
            let Some(c) = BRepTool::edge_curve_world(e_from) else {
                continue;
            };
            // `cxx:981-983`: `Geom_TrimmedCurve(C, f, l)` is skipped - the port's
            // trimmed curve rescales to [0, 1] - and the projection window below
            // carries `[f, l]` explicitly instead.
            let (f, l) = BRepTool::edge_parameters(e_from);

            // `cxx:985`.
            let Some(s_world) = BRepTool::face_surface_world(fac) else {
                continue;
            };

            // `cxx:987-994`: `GeomAdaptor_Surface Gas(S, Umin, Umax, Vmin, Vmax)`
            // and the two resolutions at `tol3d`.
            let tol3d = BRepTool::edge_tolerance(e_from).max(BRepTool::face_tolerance(fac));
            let gas = SurfaceAdaptor::of_surface_with_range(s_world.clone(), umin, umax, vmin, vmax);
            let tol_u = gas.u_resolution(tol3d);
            let tol_v = gas.v_resolution(tol3d);
            let mut tol2d = tol_u.max(tol_v);

            // `cxx:996-999`: `GeomProjLib::Curve2d(C, S, Umin, Umax, Vmin, Vmax,
            // tol2d)`; a null answer returns from the whole helper (the
            // reference's `return`, not `continue`).
            let Some(mut c2d) =
                occt_geom::geom_proj_lib::curve_2d_range(&c, f, l, &s_world, &mut tol2d)
            else {
                return Ok(());
            };

            // `cxx:1001-1002`.
            let pf = c2d.d0(f);
            let pl = c2d.d0(l);

            // `cxx:1004-1048`: the U period shift, against `[Umin, Umax]`.
            if s_world.is_u_periodic() {
                let up = s_world.u_period();
                let tolu = PCONFUSION;
                let mut nbtra = 0_i32;
                let mut the_umin = pf.x().min(pl.x());
                let mut the_umax = pf.x().max(pl.x());
                if the_umin < umin - tolu {
                    while the_umin < umin - tolu {
                        the_umin += up;
                        the_umax += up;
                        nbtra += 1;
                    }
                } else if the_umax > umax + tolu {
                    while the_umax > umax + tolu {
                        the_umax -= up;
                        the_umin -= up;
                        nbtra -= 1;
                    }
                }
                if nbtra != 0 {
                    c2d = translated(&c2d, nbtra as f64 * up, 0.0);
                }
            }

            // `cxx:1050-1094`: the V period shift.
            if s_world.is_v_periodic() {
                let vp = s_world.v_period();
                let tolv = PCONFUSION;
                let mut nbtra = 0_i32;
                let mut the_vmin = pf.y().min(pl.y());
                let mut the_vmax = pf.y().max(pl.y());
                if the_vmin < vmin - tolv {
                    while the_vmin < vmin - tolv {
                        the_vmin += vp;
                        the_vmax += vp;
                        nbtra += 1;
                    }
                } else if the_vmax > vmax + tolv {
                    while the_vmax > vmax + tolv {
                        the_vmax -= vp;
                        the_vmin -= vp;
                        nbtra -= 1;
                    }
                }
                if nbtra != 0 {
                    c2d = translated(&c2d, 0.0, nbtra as f64 * vp);
                }
            }

            // `cxx:1096`.
            crate::brep_offset::update_edge_pcurve(
                e_from,
                &c2d,
                fac,
                BRepTool::edge_tolerance(e_from),
            );
        }
        return Ok(());
    }

    // `cxx:1099-1105`: the closed-edge case, one face carrying both `Efrom` and
    // `Eto`.
    let fac = l_faces[0].clone();
    if !BRepTool::is_closed_edge_face(e_to, &fac) {
        return Err(
            "LocOpe_WiresOnShape::PutPCurves: Standard_ConstructionError".to_string(),
        );
    }

    // `cxx:1107-1141`.
    let c2dff = curve_on_surface_oriented(e_from, &fac, Orientation::Forward);
    let c2dfr = curve_on_surface_oriented(e_from, &fac, Orientation::Reversed);
    let c2dtf = curve_on_surface_oriented(e_to, &fac, Orientation::Forward);
    let c2dtr = curve_on_surface_oriented(e_to, &fac, Orientation::Reversed);

    // `cxx:1143-1147`: both points are read at the range first parameter of the
    // last `CurveOnSurface` call (`c2dtr`).
    let (Some((c2dtf, _, _)), Some((c2dtr, f_cos, _))) = (c2dtf.as_ref(), c2dtr.as_ref()) else {
        return Ok(());
    };
    let ptf = c2dtf.d0(*f_cos);
    let ptr = c2dtr.d0(*f_cos);
    let iso_u = (ptf.y() - ptr.y()).abs() < precision::epsilon(ptf.x());

    // `cxx:1150-1168`: the 3d tangents at the `Efrom` start / the projected
    // parameter, and the 2d parameter of `Efrom`'s first vertex.
    let (f_from, _) = BRepTool::edge_parameters(e_from);
    let Some(c_from) = BRepTool::edge_curve_world(e_from) else {
        return Ok(());
    };
    let (_, d1f) = c_from.d1(f_from);

    let (Some(first_vertex), _) = edge_vertices_cum_ori(e_from) else {
        return Ok(());
    };
    let vtx_param = BRepTool::parameter_of_vertex(&e_from.0, &first_vertex.0);
    let ba_curve2d = Curve2dAdaptor::of_edge_face(e_from, &fac);
    let Some(p2d) = ba_curve2d.d0(vtx_param) else {
        return Ok(());
    };
    // `cxx:1168`: `Project(FirstVertex, p2d, Eto, Fac)`.
    let prmproj = project_vertex_2d_on_edge(&p2d, e_to, &fac);

    let Some(c_to) = BRepTool::edge_curve_world(e_to) else {
        return Ok(());
    };
    let (_, d1t) = c_to.d1(prmproj);
    let same_ori = d1t.dot(&d1f) > 0.0;

    // `cxx:1181-1221`.
    let mut c2dff = c2dff.map(|(c, _, _)| c);
    let mut c2dfr = c2dfr.map(|(c, _, _)| c);
    if c2dff.is_none() && c2dfr.is_none() {
        let Some(mut s) = BRepTool::face_surface_world(&fac) else {
            return Ok(());
        };
        if let Some(basis) = s.rectangular_trimmed_basis() {
            s = basis;
        }
        let Some(c) = BRepTool::edge_curve_world(e_from) else {
            return Ok(());
        };
        let (umin, umax, vmin, vmax) = uv_bounds_face(&fac);
        let tol3d = BRepTool::edge_tolerance(e_from).max(BRepTool::face_tolerance(&fac));
        let gas = SurfaceAdaptor::of_surface_with_range(s.clone(), umin, umax, vmin, vmax);
        let mut tol2d = gas.u_resolution(tol3d).max(gas.v_resolution(tol3d));
        let (f, l) = BRepTool::edge_parameters(e_from);
        let Some(c2d) = occt_geom::geom_proj_lib::curve_2d_range(&c, f, l, &s, &mut tol2d)
        else {
            return Ok(());
        };
        c2dff = Some(c2d.clone());
        c2dfr = Some(c2d);
    } else if c2dfr.is_none() {
        c2dfr = c2dff.clone();
    } else if c2dff.is_none() {
        c2dff = c2dfr.clone();
    }
    let (Some(mut c2dff), Some(mut c2dfr)) = (c2dff, c2dfr) else {
        return Ok(());
    };

    // `cxx:1224-1227`: re-read the `Efrom` range here - the range variables were
    // overwritten by the `Eto` lookups above.
    let (f_from, _) = BRepTool::edge_parameters(e_from);
    let p2f = c2dff.d0(f_from);
    let p2r = c2dfr.d0(f_from);

    // `cxx:1230-1284`: the seam translation, in U or V.
    if iso_u {
        if same_ori {
            if (ptf.x() - p2f.x()).abs() > precision::epsilon(ptf.x()) {
                c2dff = translated(&c2dff, ptf.x() - p2f.x(), 0.0);
            }
            if (ptr.x() - p2r.x()).abs() > precision::epsilon(ptr.x()) {
                c2dfr = translated(&c2dfr, ptr.x() - p2r.x(), 0.0);
            }
        } else {
            if (ptr.x() - p2f.x()).abs() > precision::epsilon(ptr.x()) {
                c2dff = translated(&c2dff, ptr.x() - p2f.x(), 0.0);
            }
            if (ptf.x() - p2r.x()).abs() > precision::epsilon(ptf.x()) {
                c2dfr = translated(&c2dfr, ptf.x() - p2r.x(), 0.0);
            }
        }
    } else if same_ori {
        if (ptf.y() - p2f.y()).abs() > precision::epsilon(ptf.y()) {
            c2dff = translated(&c2dff, 0.0, ptf.y() - p2f.y());
        }
        if (ptr.y() - p2r.y()).abs() > precision::epsilon(ptr.y()) {
            c2dfr = translated(&c2dfr, 0.0, ptr.y() - p2r.y());
        }
    } else {
        if (ptr.y() - p2f.y()).abs() > precision::epsilon(ptr.y()) {
            c2dff = translated(&c2dff, 0.0, ptr.y() - p2f.y());
        }
        if (ptf.y() - p2r.y()).abs() > precision::epsilon(ptf.y()) {
            c2dfr = translated(&c2dfr, 0.0, ptf.y() - p2r.y());
        }
    }

    // `cxx:1285`.
    update_edge_pcurves(e_from, &c2dff, &c2dfr, &fac, BRepTool::edge_tolerance(e_from));
    Ok(())
}
