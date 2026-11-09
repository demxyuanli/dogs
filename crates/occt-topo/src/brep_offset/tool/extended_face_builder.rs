//! File-local `MakeFace` (`BRepOffset_Tool.cxx:2580-2919`) and
//! `UpdatePCurves` (`:3230-3263`), the face-construction half of
//! `EnLargeFace` (`:3437-3634`).
//!
//! `MakeFace` is a private builder rather than a `BRepLib_MakeFace` because it
//! has to keep the enlarged bounds as they are, including the infinite ones:
//! the infinite directions get no vertex and no edge, and the face is left with
//! an open (one-wire) or seam-only (two-wire) boundary.
//!
//! Port notes:
//! * `S->IsKind(STANDARD_TYPE(Geom_ElementarySurface))` (`cxx:2698`) is
//!   [`is_elementary_surface`]: the five `gp_*` hooks with the trimmed and
//!   offset wrappers excluded, because those delegate their hooks to the basis
//!   while `IsKind` is false for them.
//! * A null `UIso` / `VIso` handle is never dereferenced here; OCCT would fault
//!   on `BRepOffset_Tool::Gabarit` at `cxx:2694-2723`, and the port just skips
//!   the degeneracy test for that side.

use std::sync::Arc;

use occt_core::gp::{GpDir2d, GpPnt2d};
use occt_core::precision::CONFUSION;
use occt_geom::{Curve, Surface};
use occt_geom2d::{Curve2d, Geom2dLine};

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::edges_of;

use super::enlarge_geometry::TOL_APEX;
use super::BRepOffsetTool;

/// `S->IsKind(STANDARD_TYPE(Geom_ElementarySurface))` (`cxx:2698`): the five
/// elementary surfaces, without the wrappers that delegate their `gp_*` hooks.
fn is_elementary_surface(s: &Arc<dyn Surface>) -> bool {
    if s.rectangular_trimmed_basis().is_some() || s.is_offset_surface() {
        return false;
    }
    s.gp_pln().is_some()
        || s.gp_cylinder().is_some()
        || s.gp_cone().is_some()
        || s.gp_sphere().is_some()
        || s.gp_torus().is_some()
}

/// `B.MakeEdge(E, C, tol)` when `curve` is set, `B.MakeEdge(E)` otherwise
/// (`cxx:2737-2743`). `B.MakeEdge(E, C)` stores the curve's own range, which
/// `B.Range(E, First, Last)` then overwrites (`cxx:2762`).
fn new_edge(
    builder: &TopoBuilder,
    curve: Option<&Arc<dyn Curve>>,
    first: f64,
    last: f64,
) -> Edge {
    match curve {
        Some(c) => builder.make_edge(c.clone(), first, last),
        None => Edge(builder.make_shape(ShapeType::Edge)),
    }
}

/// `B.UpdateEdge(E, C1, C2, F, tol)` / `B.UpdateEdge(E, C, F, tol)` for pcurves
/// of type `Geom2d_Line`, in OCCT's order (first curve first). A null curve is
/// not stored (`BRep_Builder.cxx:149`); two nulls are the explicit removal of
/// `B.UpdateEdge(E, NullPCurve, NullPCurve, F, tol)`.
fn attach_lines(
    reg: &GeometryRegistry,
    edge: &Edge,
    face: &Face,
    face_key: usize,
    curves: (Option<&Geom2dLine>, Option<&Geom2dLine>),
) {
    let mut pcs: Vec<Arc<dyn Curve2d>> = Vec::new();
    if let Some(c) = curves.0 {
        pcs.push(Arc::new(c.clone()));
    }
    if let Some(c) = curves.1 {
        pcs.push(Arc::new(c.clone()));
    }
    if pcs.is_empty() {
        reg.remove_pcurves_on_surface(&edge.0, &face.0);
    } else {
        reg.set_edge_pcurves(&edge.0, face_key, pcs);
    }
}

/// `MakeFace(S, Um, UM, Vm, VM, uclosed, vclosed, isVminDegen, isVmaxDegen, F)`
/// (`cxx:2580-2919`).
///
/// `is_vmin_degen` / `is_vmax_degen` carry the cone-apex flags `EnLargeFace`
/// computed; they are extended here by the `Gabarit(VIso) <= TolApex` test of
/// `cxx:2712-2724`.
pub(super) fn make_extended_face(
    s: &Arc<dyn Surface>,
    um: f64,
    umax: f64,
    vm: f64,
    vmax: f64,
    u_closed: bool,
    v_closed: bool,
    is_vmin_degen: bool,
    is_vmax_degen: bool,
) -> Face {
    let umin = um;
    let vmin = vm;

    // `cxx:2596-2602`.
    let umin_inf = occt_core::precision::Precision::is_negative_infinite(umin);
    let umax_inf = occt_core::precision::Precision::is_positive_infinite(umax);
    let vmin_inf = occt_core::precision::Precision::is_negative_infinite(vmin);
    let vmax_inf = occt_core::precision::Precision::is_positive_infinite(vmax);

    // `cxx:2604-2625`: the apex of a conical surface makes that V bound
    // degenerate.
    let mut vmin_degen = is_vmin_degen;
    let mut vmax_degen = is_vmax_degen;
    let the_surf = s.rectangular_trimmed_basis().unwrap_or_else(|| s.clone());
    if let Some(cone) = the_surf.gp_cone() {
        let apex = cone.apex();
        let (_, vapex) = occt_core::elib::slib::cone_parameters(
            &cone.position(),
            cone.radius(),
            cone.semi_angle(),
            &apex,
        );
        if (vmin - vapex).abs() <= CONFUSION {
            vmin_degen = true;
        }
        if (vmax - vapex).abs() <= CONFUSION {
            vmax_degen = true;
        }
    }

    let builder = TopoBuilder::new();
    let tol = CONFUSION;

    // `cxx:2630-2673`. OCCT starts from four null vertices; a bound that is
    // infinite keeps its vertex null, which is why a null one is only ever
    // added under the matching `!...inf` guard below.
    let mut v00: Option<Vertex> = None;
    let mut v10: Option<Vertex> = None;
    let mut v11: Option<Vertex> = None;
    let mut v01: Option<Vertex> = None;
    if !umin_inf {
        if !vmin_inf {
            v00 = Some(builder.make_vertex(s.value(umin, vmin), tol));
        }
        if !vmax_inf {
            v01 = Some(builder.make_vertex(s.value(umin, vmax), tol));
        }
    }
    if !umax_inf {
        if !vmin_inf {
            v10 = Some(builder.make_vertex(s.value(umax, vmin), tol));
        }
        if !vmax_inf {
            v11 = Some(builder.make_vertex(s.value(umax, vmax), tol));
        }
    }

    if u_closed {
        v10 = v00.clone();
        v11 = v01.clone();
    }
    if v_closed {
        v01 = v00.clone();
        v11 = v10.clone();
    }
    if vmin_degen {
        v10 = v00.clone();
    }
    if vmax_degen {
        v11 = v01.clone();
    }

    // `cxx:2678-2695`: the pcurves of the enlarged edges.
    let dir_y = GpDir2d::new(0.0, 1.0).unwrap_or_default();
    let dir_x = GpDir2d::new(1.0, 0.0).unwrap_or_default();
    let l_umin = if !umin_inf {
        Some(Geom2dLine::from_pnt_dir(GpPnt2d::new(umin, 0.0), dir_y))
    } else {
        None
    };
    let l_umax = if !umax_inf {
        Some(Geom2dLine::from_pnt_dir(GpPnt2d::new(umax, 0.0), dir_y))
    } else {
        None
    };
    let l_vmin = if !vmin_inf {
        Some(Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, vmin), dir_x))
    } else {
        None
    };
    let l_vmax = if !vmax_inf {
        Some(Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, vmax), dir_x))
    } else {
        None
    };

    // `cxx:2697-2725`: the 3D isolines, only on an elementary surface.
    let has_iso = is_elementary_surface(s);
    let mut c_umin: Option<Arc<dyn Curve>> = None;
    let mut c_umax: Option<Arc<dyn Curve>> = None;
    let mut c_vmin: Option<Arc<dyn Curve>> = None;
    let mut c_vmax: Option<Arc<dyn Curve>> = None;
    if has_iso {
        if !umin_inf {
            c_umin = s.u_iso_curve(umin);
        }
        if !umax_inf {
            c_umax = s.u_iso_curve(umax);
        }
        if !vmin_inf {
            c_vmin = s.v_iso_curve(vmin);
            if let Some(c) = &c_vmin {
                if BRepOffsetTool::gabarit(c.as_ref()) <= TOL_APEX {
                    vmin_degen = true;
                }
            }
        }
        if !vmax_inf {
            c_vmax = s.v_iso_curve(vmax);
            if let Some(c) = &c_vmax {
                if BRepOffsetTool::gabarit(c.as_ref()) <= TOL_APEX {
                    vmax_degen = true;
                }
            }
        }
    }

    // `cxx:2728`.
    let mut face = builder.make_face_with_tolerance(s.clone(), &[], tol);
    let face_key = GeometryRegistry::shape_key(&face.0);
    let reg = GeometryRegistry::global();

    // `cxx:2732-2763`.
    let mut eumin: Option<Edge> = None;
    if !umin_inf {
        let mut edge = new_edge(&builder, c_umin.as_ref(), vmin, vmax);
        if u_closed {
            attach_lines(reg, &edge, &face, face_key, (l_umax.as_ref(), l_umin.as_ref()));
        } else {
            attach_lines(reg, &edge, &face, face_key, (None, l_umin.as_ref()));
        }
        if !vmin_inf {
            if let Some(v) = &v00 {
                let mut vv = v.0.clone();
                vv.set_orientation(Orientation::Forward);
                builder.add(&mut edge.0, &vv);
            }
        }
        if !vmax_inf {
            if let Some(v) = &v01 {
                let mut vv = v.0.clone();
                vv.set_orientation(Orientation::Reversed);
                builder.add(&mut edge.0, &vv);
            }
        }
        // `B.Range(eumin, VMin, VMax)` (`cxx:2762`).
        reg.set_ranges_all(&edge.0, vmin, vmax);
        eumin = Some(edge);
    }

    // `cxx:2765-2793`.
    let mut eumax: Option<Edge> = None;
    if !umax_inf {
        if u_closed {
            eumax = eumin.clone();
        } else {
            let mut edge = new_edge(&builder, c_umax.as_ref(), vmin, vmax);
            attach_lines(reg, &edge, &face, face_key, (None, l_umax.as_ref()));
            if !vmin_inf {
                if let Some(v) = &v10 {
                    let mut vv = v.0.clone();
                    vv.set_orientation(Orientation::Forward);
                    builder.add(&mut edge.0, &vv);
                }
            }
            if !vmax_inf {
                if let Some(v) = &v11 {
                    let mut vv = v.0.clone();
                    vv.set_orientation(Orientation::Reversed);
                    builder.add(&mut edge.0, &vv);
                }
            }
            reg.set_ranges_all(&edge.0, vmin, vmax);
            eumax = Some(edge);
        }
    }

    // `cxx:2795-2828`.
    let mut evmin: Option<Edge> = None;
    if !vmin_inf {
        let curve = if has_iso && !vmin_degen {
            c_vmin.as_ref()
        } else {
            None
        };
        let mut edge = new_edge(&builder, curve, umin, umax);
        if v_closed {
            attach_lines(reg, &edge, &face, face_key, (l_vmin.as_ref(), l_vmax.as_ref()));
        } else {
            attach_lines(reg, &edge, &face, face_key, (None, l_vmin.as_ref()));
        }
        if !umin_inf {
            if let Some(v) = &v00 {
                let mut vv = v.0.clone();
                vv.set_orientation(Orientation::Forward);
                builder.add(&mut edge.0, &vv);
            }
        }
        if !umax_inf {
            if let Some(v) = &v10 {
                let mut vv = v.0.clone();
                vv.set_orientation(Orientation::Reversed);
                builder.add(&mut edge.0, &vv);
            }
        }
        reg.set_ranges_all(&edge.0, umin, umax);
        if vmin_degen {
            reg.set_degenerated(&edge.0, true);
        }
        evmin = Some(edge);
    }

    // `cxx:2830-2863`.
    let mut evmax: Option<Edge> = None;
    if !vmax_inf {
        if v_closed {
            evmax = evmin.clone();
        } else {
            let curve = if has_iso && !vmax_degen {
                c_vmax.as_ref()
            } else {
                None
            };
            let mut edge = new_edge(&builder, curve, umin, umax);
            attach_lines(reg, &edge, &face, face_key, (None, l_vmax.as_ref()));
            if !umin_inf {
                if let Some(v) = &v01 {
                    let mut vv = v.0.clone();
                    vv.set_orientation(Orientation::Forward);
                    builder.add(&mut edge.0, &vv);
                }
            }
            if !umax_inf {
                if let Some(v) = &v11 {
                    let mut vv = v.0.clone();
                    vv.set_orientation(Orientation::Reversed);
                    builder.add(&mut edge.0, &vv);
                }
            }
            reg.set_ranges_all(&edge.0, umin, umax);
            if vmax_degen {
                reg.set_degenerated(&edge.0, true);
            }
            evmax = Some(edge);
        }
    }

    // `cxx:2866-2867`: `eumin` and `evmax` are reversed before the wires are
    // built; `eumax` and `evmin` stay FORWARD.
    let mut eumin = eumin;
    if let Some(e) = eumin.as_mut() {
        e.0.set_orientation(Orientation::Reversed);
    }
    let mut evmax = evmax;
    if let Some(e) = evmax.as_mut() {
        e.0.set_orientation(Orientation::Reversed);
    }

    // `cxx:2869-2917`.
    if !umin_inf && !umax_inf && vmin_inf && vmax_inf {
        // Two wires in U.
        let mut w = Wire::new();
        if let Some(e) = &eumin {
            builder.add_edge(&mut w, e);
        }
        builder.add_wire(&mut face, &w);
        let mut w = Wire::new();
        if let Some(e) = &eumax {
            builder.add_edge(&mut w, e);
        }
        builder.add_wire(&mut face, &w);
        face.0.set_closed(u_closed);
    } else if umin_inf && umax_inf && !vmin_inf && !vmax_inf {
        // Two wires in V.
        let mut w = Wire::new();
        if let Some(e) = &evmin {
            builder.add_edge(&mut w, e);
        }
        builder.add_wire(&mut face, &w);
        let mut w = Wire::new();
        if let Some(e) = &evmax {
            builder.add_edge(&mut w, e);
        }
        builder.add_wire(&mut face, &w);
        face.0.set_closed(v_closed);
    } else if !umin_inf || !umax_inf || !vmin_inf || !vmax_inf {
        // One wire.
        let mut w = Wire::new();
        if let Some(e) = &eumin {
            builder.add_edge(&mut w, e);
        }
        if let Some(e) = &evmin {
            builder.add_edge(&mut w, e);
        }
        if let Some(e) = &eumax {
            builder.add_edge(&mut w, e);
        }
        if let Some(e) = &evmax {
            builder.add_edge(&mut w, e);
        }
        builder.add_wire(&mut face, &w);
        w.0.set_closed(!umin_inf && !umax_inf && !vmin_inf && !vmax_inf);
        face.0.set_closed(u_closed && v_closed);
    }

    face
}

/// `UpdatePCurves(F, BF)` (`cxx:3230-3263`): move every pcurve of `F` onto the
/// enlarged face `BF`, dropping `F`'s copy.
///
/// `B.UpdateEdge(CE, NullPCurve, F, tol)` / `(CE, NullPCurve, NullPCurve, F,
/// tol)` remove the representation keyed by `F`'s surface; the following
/// `UpdateEdge(.., BF, ..)` stores it under `BF`'s (enlarged) surface. The
/// caller re-registers that surface on `F` afterwards, which is why the
/// pcurves end up visible on both faces, as in OCCT where the two faces share
/// the same `Geom_Surface` handle.
///
/// `CE.Orientation(FORWARD)` (`cxx:3244`) is a no-op for the port: pcurve slots
/// are keyed by the edge `TShape`, not by an oriented view.
pub(super) fn update_pcurves(face: &Face, bf: &Face) {
    let reg = GeometryRegistry::global();
    let f_key = GeometryRegistry::shape_key(&face.0);
    let bf_key = GeometryRegistry::shape_key(&bf.0);
    for edge in edges_of(&face.0) {
        let ce = edge;
        let pcurves = reg.edge_pcurves(&ce.0, f_key);
        // `if (!C2.IsNull())` (`cxx:3245`).
        let Some(c2) = pcurves.first().cloned() else {
            continue;
        };
        // `BRep_Tool::CurveOnSurface(CE, F, f, l)`: the COS range.
        let Some((f, l)) = reg.pcurve_range(&ce.0, f_key) else {
            // Port guard: `f`/`l` would stay uninitialised in OCCT.
            continue;
        };
        let closed = BRepTool::is_closed_edge_face(&ce, face);
        // `B.UpdateEdge(CE, NullPCurve, F, tol)` / `(.., NullPCurve,
        // NullPCurve, ..)`.
        reg.remove_pcurves_on_surface(&ce.0, &face.0);
        if closed {
            // `cxx:3249-3253`: the reversed-edge overload answers the second
            // pcurve.
            match pcurves.get(1).cloned() {
                Some(c2r) => reg.set_edge_pcurves(&ce.0, bf_key, vec![c2, c2r]),
                None => reg.set_edge_pcurve(&ce.0, bf_key, c2),
            }
        } else {
            reg.set_edge_pcurve(&ce.0, bf_key, c2);
        }
        // `B.Range(CE, f, l)` (`cxx:3260`).
        reg.set_ranges_all(&ce.0, f, l);
    }
}
