//! Analytic ImpImp pair dispatch (`iTT = iT1*10+iT2`).
//! Source: `IntPatch_ImpImpIntersection::Perform` switch and IntPP..IntSpSp.

use occt_core::gp::{GpAx1, GpPln};
use occt_core::precision::ANGULAR as PREC_ANG;
use occt_geom::intana::{
    quadric_quadric_cone_cone, quadric_quadric_cone_torus, quadric_quadric_cylinder_cone,
    quadric_quadric_cylinder_cylinder, quadric_quadric_cylinder_sphere,
    quadric_quadric_cylinder_torus, quadric_quadric_plane_cone, quadric_quadric_plane_cylinder,
    quadric_quadric_plane_sphere, quadric_quadric_plane_torus, quadric_quadric_planes,
    quadric_quadric_sphere_cone, quadric_quadric_sphere_sphere, quadric_quadric_sphere_torus,
    quadric_quadric_torus_torus, TorusIntersection,
};

use super::bounds;
use super::cycy::{cy_cy_no_geometric, CylUv};
use super::glines::{from_circles, from_quadric_intersection, PairOutcome};
use super::quad::{axis_distance, dist_point_axis, ImplicitQuad};

const TOL_ANG: f64 = 1.0e-8;

pub(crate) fn intersect_pair(
    q1: &ImplicitQuad,
    q2: &ImplicitQuad,
    tol_tang: f64,
    uv1: CylUv,
    uv2: CylUv,
) -> PairOutcome {
    let i_tt = q1.code() * 10 + q2.code();
    let reversed = q1.code() > q2.code();
    match i_tt {
        11 => int_pp(q1, q2, tol_tang),
        12 | 21 => int_p_cy(q1, q2, reversed, tol_tang),
        13 | 31 => int_p_co(q1, q2, reversed, tol_tang),
        14 | 41 => int_p_sp(q1, q2, reversed),
        15 | 51 => int_p_to(q1, q2, reversed, tol_tang),
        22 => int_cy_cy(q1, q2, tol_tang, uv1, uv2),
        23 | 32 => int_cy_co(q1, q2, reversed, tol_tang),
        24 | 42 => int_cy_sp(q1, q2, reversed, tol_tang),
        33 => int_co_co(q1, q2, tol_tang),
        34 | 43 => int_co_sp(q1, q2, reversed, tol_tang),
        44 => int_sp_sp(q1, q2, tol_tang),
        25 | 52 => int_cy_to(q1, q2, reversed, tol_tang),
        35 | 53 => int_co_to(q1, q2, reversed, tol_tang),
        45 | 54 => int_sp_to(q1, q2, reversed, tol_tang),
        55 => int_to_to(q1, q2, tol_tang),
        _ => PairOutcome::Fail,
    }
}

fn from_torus(r: TorusIntersection, seam: Option<&ImplicitQuad>) -> PairOutcome {
    match r {
        TorusIntersection::Fail => PairOutcome::Fail,
        TorusIntersection::Empty => PairOutcome::Empty,
        TorusIntersection::Same => PairOutcome::Same,
        TorusIntersection::Circles(c) => from_circles(c, seam),
    }
}

fn axes_normal(a: &occt_core::gp::GpDir, b: &occt_core::gp::GpDir) -> bool {
    (std::f64::consts::FRAC_PI_2 - a.angle(b)).abs() <= PREC_ANG
}

fn plane_of<'a>(q: &'a ImplicitQuad, other: &'a ImplicitQuad, reversed: bool) -> &'a GpPln {
    if reversed {
        other.as_plane().expect("plane")
    } else {
        q.as_plane().expect("plane")
    }
}

fn other_of<'a>(q: &'a ImplicitQuad, other: &'a ImplicitQuad, reversed: bool) -> &'a ImplicitQuad {
    if reversed {
        q
    } else {
        other
    }
}

fn int_pp(q1: &ImplicitQuad, q2: &ImplicitQuad, tol_tang: f64) -> PairOutcome {
    let p1 = q1.as_plane().expect("plane");
    let p2 = q2.as_plane().expect("plane");
    from_quadric_intersection(quadric_quadric_planes(p1, p2, TOL_ANG, tol_tang), None)
}

fn int_p_cy(q1: &ImplicitQuad, q2: &ImplicitQuad, reversed: bool, tol_tang: f64) -> PairOutcome {
    let pl = plane_of(q1, q2, reversed);
    let cy_q = other_of(q1, q2, reversed);
    let cy = cy_q.as_cylinder().expect("cylinder");
    from_quadric_intersection(
        quadric_quadric_plane_cylinder(pl, cy, TOL_ANG, tol_tang),
        Some(cy_q),
    )
}

fn int_p_co(q1: &ImplicitQuad, q2: &ImplicitQuad, reversed: bool, tol_tang: f64) -> PairOutcome {
    let pl = plane_of(q1, q2, reversed);
    let co_q = other_of(q1, q2, reversed);
    let co = co_q.as_cone().expect("cone");
    from_quadric_intersection(
        quadric_quadric_plane_cone(pl, co, TOL_ANG, tol_tang),
        Some(co_q),
    )
}

fn int_p_sp(q1: &ImplicitQuad, q2: &ImplicitQuad, reversed: bool) -> PairOutcome {
    let pl = plane_of(q1, q2, reversed);
    let sp_q = other_of(q1, q2, reversed);
    let sp = sp_q.as_sphere().expect("sphere");
    from_quadric_intersection(quadric_quadric_plane_sphere(pl, sp), Some(sp_q))
}

fn int_p_to(q1: &ImplicitQuad, q2: &ImplicitQuad, reversed: bool, tol_tang: f64) -> PairOutcome {
    let pl = plane_of(q1, q2, reversed);
    let to_q = other_of(q1, q2, reversed);
    let to = to_q.as_torus().expect("torus");
    let r = quadric_quadric_plane_torus(pl, to, tol_tang);
    let pl_ax = pl.axis();
    let seam = if axes_normal(pl_ax.direction(), to.axis().direction()) {
        None
    } else {
        Some(to_q)
    };
    from_torus(r, seam)
}

fn int_cy_cy(
    q1: &ImplicitQuad,
    q2: &ImplicitQuad,
    tol_tang: f64,
    uv1: CylUv,
    uv2: CylUv,
) -> PairOutcome {
    let c1 = q1.as_cylinder().expect("cylinder");
    let c2 = q2.as_cylinder().expect("cylinder");
    let qi = quadric_quadric_cylinder_cylinder(c1, c2, tol_tang);
    let out = from_quadric_intersection(qi, Some(q1));
    let parallel = c1.axis().direction().is_parallel(c2.axis().direction());
    if matches!(out, PairOutcome::Empty) && !parallel {
        let a2d = (tol_tang / c1.radius().max(c2.radius()).max(1.0)).min(1.0e-4);
        match cy_cy_no_geometric(c1, c2, uv1, uv2, tol_tang, a2d) {
            Ok(r) => r,
            Err(true) => PairOutcome::Infinite,
            Err(false) => PairOutcome::Fail,
        }
    } else {
        out
    }
}

fn int_cy_to(q1: &ImplicitQuad, q2: &ImplicitQuad, reversed: bool, tol_tang: f64) -> PairOutcome {
    let (cy_q, to_q) = if reversed { (q2, q1) } else { (q1, q2) };
    let cy = cy_q.as_cylinder().expect("cylinder");
    let to = to_q.as_torus().expect("torus");
    from_torus(quadric_quadric_cylinder_torus(cy, to, tol_tang), None)
}

fn int_co_to(q1: &ImplicitQuad, q2: &ImplicitQuad, reversed: bool, tol_tang: f64) -> PairOutcome {
    let (co_q, to_q) = if reversed { (q2, q1) } else { (q1, q2) };
    let co = co_q.as_cone().expect("cone");
    let to = to_q.as_torus().expect("torus");
    from_torus(quadric_quadric_cone_torus(co, to, tol_tang), None)
}

fn int_sp_to(q1: &ImplicitQuad, q2: &ImplicitQuad, reversed: bool, tol_tang: f64) -> PairOutcome {
    let (sp_q, to_q) = if reversed { (q2, q1) } else { (q1, q2) };
    let sp = sp_q.as_sphere().expect("sphere");
    let to = to_q.as_torus().expect("torus");
    from_torus(quadric_quadric_sphere_torus(sp, to, tol_tang), None)
}

fn int_to_to(q1: &ImplicitQuad, q2: &ImplicitQuad, tol_tang: f64) -> PairOutcome {
    let t1 = q1.as_torus().expect("torus");
    let t2 = q2.as_torus().expect("torus");
    from_torus(quadric_quadric_torus_torus(t1, t2, tol_tang), Some(q1))
}

fn int_cy_co(q1: &ImplicitQuad, q2: &ImplicitQuad, reversed: bool, tol_tang: f64) -> PairOutcome {
    let (cy_q, co_q) = if reversed { (q2, q1) } else { (q1, q2) };
    let cy = cy_q.as_cylinder().expect("cylinder");
    let co = co_q.as_cone().expect("cone");
    if !axes_coaxial(&cy.axis(), &co.axis()) {
        return bounds::from_cyl_quad(q1, q2, cy, co_q, Some(co), false, tol_tang);
    }
    from_quadric_intersection(quadric_quadric_cylinder_cone(cy, co, tol_tang), Some(cy_q))
}

fn int_cy_sp(q1: &ImplicitQuad, q2: &ImplicitQuad, _reversed: bool, tol_tang: f64) -> PairOutcome {
    let (cy_q, sp_q) = match (q1.as_cylinder(), q2.as_sphere()) {
        (Some(_), Some(_)) => (q1, q2),
        _ => (q2, q1),
    };
    let cy = cy_q.as_cylinder().expect("cylinder");
    let sp = sp_q.as_sphere().expect("sphere");
    if dist_point_axis(&sp.location(), &cy.axis()) > 1e-9 {
        return bounds::from_cyl_quad(q1, q2, cy, sp_q, None, true, tol_tang);
    }
    from_quadric_intersection(quadric_quadric_cylinder_sphere(cy, sp, tol_tang), Some(cy_q))
}

fn int_co_co(q1: &ImplicitQuad, q2: &ImplicitQuad, tol_tang: f64) -> PairOutcome {
    let c1 = q1.as_cone().expect("cone");
    let c2 = q2.as_cone().expect("cone");
    let qi = quadric_quadric_cone_cone(c1, c2, TOL_ANG, tol_tang);
    let out = from_quadric_intersection(qi, Some(q1));
    if matches!(out, PairOutcome::Empty) {
        let coaxial = axes_coaxial(&c1.axis(), &c2.axis());
        let parallel = c1.axis().direction().is_parallel(c2.axis().direction());
        let same_apex = c1.apex().square_distance(&c2.apex()) < tol_tang * tol_tang;
        if coaxial || (parallel && (c1.semi_angle() - c2.semi_angle()).abs() < TOL_ANG) || same_apex
        {
            out
        } else {
            bounds::from_cone_quad(q1, q2, c1, q2, false, false, tol_tang)
        }
    } else {
        out
    }
}

fn int_co_sp(q1: &ImplicitQuad, q2: &ImplicitQuad, _reversed: bool, tol_tang: f64) -> PairOutcome {
    let (sp_q, co_q) = match (q1.as_sphere(), q2.as_cone()) {
        (Some(_), Some(_)) => (q1, q2),
        _ => (q2, q1),
    };
    let sp = sp_q.as_sphere().expect("sphere");
    let co = co_q.as_cone().expect("cone");
    if dist_point_axis(&sp.location(), &co.axis()) > 1e-9 {
        return bounds::from_cone_quad(q1, q2, co, sp_q, false, false, tol_tang);
    }
    from_quadric_intersection(quadric_quadric_sphere_cone(sp, co, tol_tang), Some(co_q))
}

fn int_sp_sp(q1: &ImplicitQuad, q2: &ImplicitQuad, tol_tang: f64) -> PairOutcome {
    let s1 = q1.as_sphere().expect("sphere");
    let s2 = q2.as_sphere().expect("sphere");
    from_quadric_intersection(quadric_quadric_sphere_sphere(s1, s2, tol_tang), Some(q1))
}

fn axes_coaxial(a: &GpAx1, b: &GpAx1) -> bool {
    a.direction().is_parallel(b.direction()) && axis_distance(a, b) < 1e-14
}
