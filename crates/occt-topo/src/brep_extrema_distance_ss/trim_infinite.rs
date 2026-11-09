//! `TRIM_INFINIT_EDGE` and `TRIM_INFINIT_FACE`.

use super::prelude::*;


/// `TRIM_INFINIT_EDGE` (`cxx:102-245`).
///
/// Returns `(res_edge, is_trim1, is_trim2)`; `res_edge` defaults to `s2`
/// (`cxx:113`).
pub(super) fn trim_infinit_edge(
    s1: &Edge,
    s2: &Edge,
) -> (Option<Edge>, bool, bool) {
    let mut b_is_trim1 = false;
    let mut b_is_trim2 = false;
    // `cxx:108-111`.
    if BRepTool::is_degenerated(s1) || BRepTool::is_degenerated(s2) {
        return (None, b_is_trim1, b_is_trim2);
    }
    // `cxx:114-115`.
    let (a_first1, a_last1) = BRepTool::edge_parameters(s1);
    let (a_first2, a_last2) = BRepTool::edge_parameters(s2);
    let p_curv1 = BRepTool::edge_curve(s1);
    let p_curv2 = BRepTool::edge_curve(s2);

    // `cxx:117-121`.
    if Precision::is_infinite(a_first1)
        && Precision::is_infinite(a_last1)
        && Precision::is_infinite(a_first2)
        && Precision::is_infinite(a_last2)
    {
        return (None, b_is_trim1, b_is_trim2);
    }

    // `cxx:123-126`.
    let mut u_min = 0.0;
    let mut u_max = 0.0;
    let mut b_umin = false;
    let mut b_umax = false;

    // `cxx:128-158`: pick the infinite curve and its known-bound side.
    let mut p_curv: Option<Arc<dyn Curve>> = None;
    if let Some(c) = &p_curv1 {
        if Precision::is_infinite(a_first1) || Precision::is_infinite(a_last1) {
            p_curv = Some(c.clone());
            b_is_trim1 = true;
            if !Precision::is_infinite(a_first1) {
                b_umin = true;
                u_min = a_first1;
            } else if !Precision::is_infinite(a_last1) {
                b_umax = true;
                u_max = a_last1;
            }
        }
    }
    if p_curv.is_none() {
        if let Some(c) = &p_curv2 {
            if Precision::is_infinite(a_first2) || Precision::is_infinite(a_last2) {
                p_curv = Some(c.clone());
                b_is_trim2 = true;
                if !Precision::is_infinite(a_first2) {
                    b_umin = true;
                    u_min = a_first2;
                } else if !Precision::is_infinite(a_last2) {
                    b_umax = true;
                    u_max = a_last2;
                }
            }
        }
    }

    let p_curv = match p_curv {
        Some(c) => c,
        // `cxx:159`: neither edge was infinite.
        None => return (None, b_is_trim1, b_is_trim2),
    };

    // `cxx:161-181`: the other edge's bounding box, then its 8 corners.
    let other = if b_is_trim1 { s2 } else { s1 };
    let mut a_edge_box = BndBox::new();
    brep_bnd_lib::add(&other.0, &mut a_edge_box, false);
    let (xmin, ymin, zmin, xmax, ymax, zmax) = match a_edge_box.get() {
        Some(v) => v,
        None => return (None, b_is_trim1, b_is_trim2),
    };
    let corners = [
        GpPnt::new(xmin, ymin, zmin),
        GpPnt::new(xmin, ymax, zmin),
        GpPnt::new(xmin, ymax, zmax),
        GpPnt::new(xmin, ymin, zmax),
        GpPnt::new(xmax, ymax, zmin),
        GpPnt::new(xmax, ymax, zmax),
        GpPnt::new(xmax, ymin, zmax),
        GpPnt::new(xmax, ymin, zmin),
    ];

    // `cxx:183-198`: `GeomAPI_ProjectPointOnCurve`'s `LowerDistanceParameter`
    // for each corner.
    let mut arr_u = [0.0f64; 8];
    for (i, c) in corners.iter().enumerate() {
        arr_u[i] = occt_geom::geom_api::project_point_on_curve(p_curv.as_ref(), c, Precision::CONFUSION)
            .map(|p| p.parameter)
            .unwrap_or(0.0);
    }

    // `cxx:200-207`.
    if !b_umin {
        u_min = arr_u[0];
    }
    if !b_umax {
        u_max = arr_u[0];
    }
    // `cxx:209-223`.
    for a_u in arr_u {
        if a_u < u_min {
            u_min = a_u;
        } else if a_u > u_max {
            u_max = a_u;
        }
    }

    // `cxx:225-233`.
    let mut tol = Precision::CONFUSION;
    if b_is_trim1 {
        tol = BRepTool::edge_tolerance(s1);
    } else if b_is_trim2 {
        tol = BRepTool::edge_tolerance(s2);
    }

    // `cxx:235-239`.
    let eps_u = p_curv.resolution(3.0 * tol);
    if (u_min - u_max).abs() < eps_u {
        u_min -= eps_u;
        u_max += eps_u;
    }

    // `cxx:241-242`: `new Geom_TrimmedCurve(pCurv, Umin, Umax)` then
    // `BRepBuilderAPI_MakeEdge(result)`.
    let trimmed = occt_geom::trimmed::GeomTrimmedCurve::new(p_curv, u_min, u_max);
    let b = TopoBuilder::new();
    let edge = b.make_edge(Arc::new(trimmed), u_min, u_max);
    (Some(edge), b_is_trim1, b_is_trim2)
}

/// `TRIM_INFINIT_FACE` (`cxx:249-484`).
///
/// Returns `(res_face, is_infinite)`; `res_face` is the second face argument
/// when the pair has no infinite face (`cxx:284`).
pub(super) fn trim_infinit_face(s1: &TopoShape, s2: &TopoShape) -> (Option<Face>, bool) {
    // `cxx:252`.
    let mut b_is_infinit = false;

    // `cxx:254-283`: orient the (edge, face) pair.
    let (a_e, a_f) = match (s1.shape_type(), s2.shape_type()) {
        (ShapeType::Edge, ShapeType::Face) => {
            let e = match Edge::wrap(s1.clone()) {
                Some(e) => e,
                None => return (None, b_is_infinit),
            };
            // `cxx:264-267`.
            if BRepTool::is_degenerated(&e) {
                return (None, b_is_infinit);
            }
            match Face::wrap(s2.clone()) {
                Some(f) => (e, f),
                None => return (None, b_is_infinit),
            }
        }
        (ShapeType::Face, ShapeType::Edge) => {
            let e = match Edge::wrap(s2.clone()) {
                Some(e) => e,
                None => return (None, b_is_infinit),
            };
            // `cxx:272-275`.
            if BRepTool::is_degenerated(&e) {
                return (None, b_is_infinit);
            }
            match Face::wrap(s1.clone()) {
                Some(f) => (e, f),
                None => return (None, b_is_infinit),
            }
        }
        // `cxx:278-283`.
        _ => return (None, false),
    };

    // `cxx:285-287`.
    let p_surf = match BRepTool::face_surface(&a_f) {
        Some(s) => s,
        None => return (None, b_is_infinit),
    };
    // `cxx:288`.
    let b_restrict = BRepTool::natural_restriction(&a_f);

    // `cxx:290-293`.
    let mut u_min = f64::MAX;
    let mut u_max = f64::MIN;
    let mut v_min = f64::MAX;
    let mut v_max = f64::MIN;
    let (mut b_umin, mut b_umax, mut b_vmin, mut b_vmax) = (false, false, false, false);
    let mut b_is_trim = false;

    // `cxx:295-337`.
    if b_restrict {
        let (u1, u2) = p_surf.u_range();
        let (v1, v2) = p_surf.v_range();
        if Precision::is_infinite(u1) {
            b_is_trim = true;
        } else {
            u_min = u1;
            b_umin = true;
        }
        if Precision::is_infinite(u2) {
            b_is_trim = true;
        } else {
            u_max = u2;
            b_umax = true;
        }
        if Precision::is_infinite(v1) {
            b_is_trim = true;
        } else {
            v_min = v1;
            b_vmin = true;
        }
        if Precision::is_infinite(v2) {
            b_is_trim = true;
        } else {
            v_max = v2;
            b_vmax = true;
        }
    } else {
        // `cxx:340-347`.
        let (u1, u2, v1, v2) = BRepTool::uv_bounds(&a_f);
        if Precision::is_infinite(u1)
            && Precision::is_infinite(u2)
            && Precision::is_infinite(v1)
            && Precision::is_infinite(v2)
        {
            b_is_trim = true;
        }
    }

    if !b_is_trim {
        // `cxx:478-482`.
        return (Some(a_f), false);
    }

    // `cxx:351-355`.
    let mut a_edge_box = BndBox::new();
    brep_bnd_lib::add(&a_e.0, &mut a_edge_box, false);
    if a_edge_box.is_whole() {
        return (None, b_is_infinit);
    }
    let (xmin, ymin, zmin, xmax, ymax, zmax) = match a_edge_box.get() {
        Some(v) => v,
        None => return (None, b_is_infinit),
    };
    let corners = [
        GpPnt::new(xmin, ymin, zmin),
        GpPnt::new(xmin, ymax, zmin),
        GpPnt::new(xmin, ymax, zmax),
        GpPnt::new(xmin, ymin, zmax),
        GpPnt::new(xmax, ymax, zmin),
        GpPnt::new(xmax, ymax, zmax),
        GpPnt::new(xmax, ymin, zmax),
        GpPnt::new(xmax, ymin, zmin),
    ];

    // `cxx:372-410`: `GeomAPI_ProjectPointOnSurf`'s
    // `LowerDistanceParameters`; failed projections leave `(0, 0)`.
    let mut arr_u = [0.0f64; 8];
    let mut arr_v = [0.0f64; 8];
    for (i, c) in corners.iter().enumerate() {
        if let Some(p) =
            occt_geom::geom_api::project_point_on_surface(p_surf.as_ref(), c, Precision::CONFUSION)
        {
            arr_u[i] = p.u;
            arr_v[i] = p.v;
        }
    }

    // `cxx:412-426`.
    if !b_umin {
        u_min = arr_u[0];
    }
    if !b_umax {
        u_max = arr_u[0];
    }
    if !b_vmin {
        v_min = arr_v[0];
    }
    if !b_vmax {
        v_max = arr_v[0];
    }
    // `cxx:428-453`.
    for i in 0..8 {
        let a_u = arr_u[i];
        if a_u < u_min {
            u_min = a_u;
        } else if a_u > u_max {
            u_max = a_u;
        }
        let a_v = arr_v[i];
        if a_v < v_min {
            v_min = a_v;
        } else if a_v > v_max {
            v_max = a_v;
        }
    }

    // `cxx:455-470`.
    let tol = BRepTool::face_tolerance(&a_f);
    let eps_u = occt_geom::approx_same_parameter::u_resolution(p_surf.as_ref(), 3.0 * tol);
    if (u_min - u_max).abs() < eps_u {
        u_min -= eps_u;
        u_max += eps_u;
    }
    let eps_v = occt_geom::approx_same_parameter::v_resolution(p_surf.as_ref(), 3.0 * tol);
    if (v_min - v_max).abs() < eps_v {
        v_min -= eps_v;
        v_max += eps_v;
    }

    // `cxx:471-475`.
    let trimmed = match occt_geom::rectangular_trimmed::GeomRectangularTrimmedSurface::uv(
        p_surf, u_min, u_max, v_min, v_max,
    ) {
        Ok(t) => t,
        Err(_) => return (None, b_is_infinit),
    };
    let b = TopoBuilder::new();
    let face = b.make_face_with_tolerance(Arc::new(trimmed), &[], Precision::CONFUSION);
    b_is_infinit = true;
    (Some(face), b_is_infinit)
}
