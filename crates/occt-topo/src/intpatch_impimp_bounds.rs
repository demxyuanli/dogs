//! `ProcessBounds` plus wrapping of `IntAna_IntQuadQuad` into `IntPatch_ALine`.
//! Source: `IntPatch_ImpImpIntersection.cxx` ProcessBounds / IntCySp..IntCoSp / ExploreCurve.

use std::sync::Arc;

use occt_core::gp::{GpCone, GpPnt};
use occt_core::precision::PCONFUSION;
use occt_geom::intana::{IntAnaCurve, IntAnaQuadric, IntQuadQuad};

use crate::geom_int::{ALine, GeomIntLine};
use crate::int_tools_wline::PatchPoint;

use super::glines::PairOutcome;
use super::quad::ImplicitQuad;

/// `ProcessBounds`: share ALine endpoints with existing Analytic lines.
pub(crate) fn process_bounds(
    alig: &mut ALine,
    slin: &mut [GeomIntLine],
    q1: &ImplicitQuad,
    q2: &ImplicitQuad,
    procf: &mut bool,
    ptf: GpPnt,
    first: f64,
    procl: &mut bool,
    ptl: GpPnt,
    last: f64,
    multpoint: &mut bool,
    tol: f64,
) {
    let mut j = if *procf && *procl { slin.len() } else { 0 };
    while j < slin.len() {
        let Some(aligold) = slin[j].as_aline_mut() else {
            j += 1;
            continue;
        };
        let mut k = 1;
        while k <= aligold.nb_vertex() {
            let mut ptsol = *aligold.vertex(k);
            if !*procf && ptf.distance(&ptsol.p) <= tol {
                ptsol.tolerance = tol;
                if !ptsol.is_multiple {
                    *multpoint = true;
                    ptsol.is_multiple = true;
                    aligold.replace(k, ptsol);
                }
                ptsol.param_on_line = first;
                alig.add_vertex(ptsol);
                alig.set_first_point(alig.nb_vertex());
                *procf = true;
                ptsol = *aligold.vertex(k);
            }
            if !*procl && ptl.distance(&ptsol.p) <= tol {
                ptsol.tolerance = tol;
                if !ptsol.is_multiple {
                    *multpoint = true;
                    ptsol.is_multiple = true;
                    aligold.replace(k, ptsol);
                }
                ptsol.param_on_line = last;
                alig.add_vertex(ptsol);
                alig.set_last_point(alig.nb_vertex());
                *procl = true;
                ptsol = *aligold.vertex(k);
            }
            if *procf && *procl {
                k = aligold.nb_vertex() + 1;
            } else {
                k += 1;
            }
        }
        if *procf && *procl {
            j = slin.len();
        } else {
            j += 1;
        }
    }

    if !*procf && !*procl {
        let (u1, v1) = q1.parameters(&ptf);
        let (u2, v2) = q2.parameters(&ptf);
        let mut ptsol = PatchPoint::new(ptf, first, u1, v1, u2, v2);
        ptsol.tolerance = tol;
        if ptf.distance(&ptl) <= tol {
            ptsol.is_multiple = true;
            *multpoint = true;
            alig.add_vertex(ptsol);
            alig.set_first_point(alig.nb_vertex());
            ptsol.param_on_line = last;
            alig.add_vertex(ptsol);
            alig.set_last_point(alig.nb_vertex());
        } else {
            alig.add_vertex(ptsol);
            alig.set_first_point(alig.nb_vertex());
            let (u1, v1) = q1.parameters(&ptl);
            let (u2, v2) = q2.parameters(&ptl);
            let mut ptsol = PatchPoint::new(ptl, last, u1, v1, u2, v2);
            ptsol.tolerance = tol;
            alig.add_vertex(ptsol);
            alig.set_last_point(alig.nb_vertex());
        }
    } else if !*procf {
        let (u1, v1) = q1.parameters(&ptf);
        let (u2, v2) = q2.parameters(&ptf);
        let mut ptsol = PatchPoint::new(ptf, first, u1, v1, u2, v2);
        ptsol.tolerance = tol;
        alig.add_vertex(ptsol);
        alig.set_first_point(alig.nb_vertex());
    } else if !*procl {
        let (u1, v1) = q1.parameters(&ptl);
        let (u2, v2) = q2.parameters(&ptl);
        let mut ptsol = PatchPoint::new(ptl, last, u1, v1, u2, v2);
        ptsol.tolerance = tol;
        alig.add_vertex(ptsol);
        alig.set_last_point(alig.nb_vertex());
    }
}

/// `ExploreCurve` — split an IntAna curve at the cone apex.
pub(crate) fn explore_curve(cone: &GpCone, crv: IntAnaCurve, tol: f64) -> Vec<IntAnaCurve> {
    let sq_tol = tol * tol;
    let apex = cone.apex();
    let (mut t1, t2) = crv.domain();
    let params = crv.find_parameter(&apex);
    if params.is_empty() {
        return vec![crv];
    }
    let mut out = Vec::new();
    for mut prm in params {
        if prm - t1 < PCONFUSION {
            continue;
        }
        let mut is_last = false;
        if t2 - prm < PCONFUSION {
            prm = t2;
            is_last = true;
        }
        let p = crv.value(prm);
        if p.square_distance(&apex) < sq_tol {
            let mut c1 = crv.clone();
            c1.set_domain(t1, prm);
            t1 = prm;
            out.push(c1);
        }
        if is_last {
            break;
        }
    }
    if out.is_empty() {
        return vec![crv];
    }
    if t2 - t1 > PCONFUSION {
        let mut c1 = crv;
        c1.set_domain(t1, t2);
        out.push(c1);
    }
    out
}

fn wrap_int_quad_quad(
    q1: &ImplicitQuad,
    q2: &ImplicitQuad,
    iqq: IntQuadQuad,
    cone_split: Option<&GpCone>,
    always_closed: bool,
    tol: f64,
) -> PairOutcome {
    if !iqq.is_done() {
        return PairOutcome::Fail;
    }
    if iqq.identical_elements() {
        return PairOutcome::Same;
    }
    if iqq.nb_curve() == 0 && iqq.nb_pnt() == 0 {
        return PairOutcome::Empty;
    }
    let mut lines = Vec::new();
    let mut points = Vec::new();
    for i in 1..=iqq.nb_pnt() {
        let p = iqq.point(i);
        let (u1, v1) = q1.parameters(&p);
        let (u2, v2) = q2.parameters(&p);
        let mut pt = PatchPoint::new(p, 0.0, u1, v1, u2, v2);
        pt.tolerance = tol;
        points.push(pt);
    }
    let mut mult = false;
    for i in 1..=iqq.nb_curve() {
        let c0 = iqq.curve(i).clone();
        let pieces = if let Some(co) = cone_split {
            explore_curve(co, c0, 10.0 * tol)
        } else {
            vec![c0]
        };
        for c in pieces {
            let (first, last) = c.domain();
            let first_open = c.is_first_open();
            let last_open = c.is_last_open();
            let ptf = if !first_open {
                c.value(first)
            } else {
                GpPnt::zero()
            };
            let ptl = if !last_open {
                c.value(last)
            } else {
                GpPnt::zero()
            };
            let mut alig = ALine::from_curve(Arc::new(c));
            let mut procf = if always_closed { false } else { first_open };
            let mut procl = if always_closed { false } else { last_open };
            process_bounds(
                &mut alig,
                &mut lines,
                q1,
                q2,
                &mut procf,
                ptf,
                first,
                &mut procl,
                ptl,
                last,
                &mut mult,
                tol,
            );
            lines.push(GeomIntLine::Analytic(alig));
        }
    }
    PairOutcome::Result { lines, points }
}

pub(crate) fn from_cyl_quad(
    q1: &ImplicitQuad,
    q2: &ImplicitQuad,
    cyl: &occt_core::gp::GpCylinder,
    quad: &ImplicitQuad,
    cone_split: Option<&GpCone>,
    always_closed: bool,
    tol: f64,
) -> PairOutcome {
    let Some(ana) = to_ana_quadric(quad) else {
        return PairOutcome::Fail;
    };
    let iqq = IntQuadQuad::cylinder_quad(cyl, &ana, tol);
    wrap_int_quad_quad(q1, q2, iqq, cone_split, always_closed, tol)
}

pub(crate) fn from_cone_quad(
    q1: &ImplicitQuad,
    q2: &ImplicitQuad,
    cone: &GpCone,
    quad: &ImplicitQuad,
    cone_split: bool,
    always_closed: bool,
    tol: f64,
) -> PairOutcome {
    let Some(ana) = to_ana_quadric(quad) else {
        return PairOutcome::Fail;
    };
    let iqq = IntQuadQuad::cone_quad(cone, &ana, tol);
    let split = if cone_split { Some(cone) } else { None };
    wrap_int_quad_quad(q1, q2, iqq, split, always_closed, tol)
}

fn to_ana_quadric(q: &ImplicitQuad) -> Option<IntAnaQuadric> {
    match q {
        ImplicitQuad::Plane(p) => Some(IntAnaQuadric::from_plane(p)),
        ImplicitQuad::Cylinder(c) => Some(IntAnaQuadric::from_cylinder(c)),
        ImplicitQuad::Cone(c) => Some(IntAnaQuadric::from_cone(c)),
        ImplicitQuad::Sphere(s) => Some(IntAnaQuadric::from_sphere(s)),
        ImplicitQuad::Torus(_) => None,
    }
}
