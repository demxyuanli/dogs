//! `BRepFeat` statics: edge sampling, barycentres, parametric extents, face
//! containment, the `FaceUntil` face extension and `Tool`.
//!
//! Source: `ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat.hxx:87-119`,
//! `BRepFeat.cxx:59-808`.

use std::sync::Arc;

use occt_core::bnd::BndBox;
use occt_core::elib::slib;
use occt_core::gp::{GpPnt, GpTrsf2d, GpVec2d};
use occt_core::precision::{self};
use occt_geom::curve::Curve;
use occt_geom::extrema_pc::extrema_ext_pc_range_done;
use occt_geom::gcpnts::QuasiUniformDeflection;
use occt_geom::geom_proj_lib::curve_2d_range;
use occt_geom::rectangular_trimmed::GeomRectangularTrimmedSurface;
use occt_geom2d::adaptor_curve::Geom2dAdaptorCurve;
use occt_geom2d::curve::Curve2d;

use crate::abs::{Orientation, ShapeType};
use crate::brep_bnd_lib;
use crate::brep_lib_make_face;
use crate::brep_tool::BRepTool;
use crate::brep_top_adaptor_fclass2d::BrepTopAdaptorFClass2d;
use crate::builder::TopoBuilder;
use crate::fclass2d::FaceState;
use crate::locope::{self, BuildShape, CSIntersector, ShapeSet};
use crate::shape::{Face, Shell, Solid, TopoShape};
use crate::topo_tools;

/// `#define NECHANTBARYC 11` (`BRepFeat.cxx:57`).
const NECHANTBARYC: i32 = 11;

/// `RealLast()` (`Standard_Real.hxx`).
const REAL_LAST: f64 = f64::MAX;
/// `RealFirst()` (`Standard_Real.hxx`).
const REAL_FIRST: f64 = -f64::MAX;

/// `BRepFeat::SampleEdges(S, Pt)` (`BRepFeat.cxx:61-64`): forwarded to
/// `LocOpe::SampleEdges`.
pub fn sample_edges(shape: &TopoShape) -> Vec<GpPnt> {
    locope::sample_edges(shape)
}

/// `BRepFeat::Barycenter(S, B)` (`BRepFeat.cxx:71-112`): the average of the
/// `NECHANTBARYC - 1` interior samples of every non-degenerated edge plus
/// every vertex, each counted once.
///
/// A null `BRep_Tool::Curve` dereferences in the reference (`cxx:91-92`); the
/// port skips such an edge and documents the guard here.
pub fn barycenter(shape: &TopoShape) -> GpPnt {
    let mut the_map = ShapeSet::new();
    let mut bar = [0.0f64; 3];
    let mut nbp = 0usize;

    for edg in topo_tools::edges_of(shape) {
        if !the_map.add(&edg.0) {
            continue;
        }
        if BRepTool::is_degenerated(&edg) {
            continue;
        }
        let (f, l) = BRepTool::edge_parameters(&edg);
        let Some(c) = BRepTool::edge_curve_world(&edg) else {
            continue;
        };
        for i in 1..NECHANTBARYC {
            let prm = ((NECHANTBARYC - i) as f64 * f + i as f64 * l) / NECHANTBARYC as f64;
            let p = c.d0(prm);
            bar[0] += p.x();
            bar[1] += p.y();
            bar[2] += p.z();
            nbp += 1;
        }
    }

    // `// Adds every vertex` (`cxx:100-108`).
    for v in topo_tools::vertices_of(shape) {
        if the_map.add(&v.0) {
            let p = BRepTool::vertex_point(&v);
            bar[0] += p.x();
            bar[1] += p.y();
            bar[2] += p.z();
            nbp += 1;
        }
    }

    let d = nbp as f64;
    GpPnt::new(bar[0] / d, bar[1] / d, bar[2] / d)
}

/// The `Extrema_ExtPC` projection of `p` onto `cc` over the whole curve range,
/// picking the solution of minimal square distance (`BRepFeat.cxx:147-165`,
/// `:176-189`, `:272-286`, `:303-321`). `None` is the reference's
/// `!IsDone() || NbExt() < 1` branch.
fn project_min_param(cc: &Arc<dyn Curve>, p: &GpPnt) -> Option<f64> {
    let (sols, done) = extrema_ext_pc_range_done(
        cc.as_ref(),
        p,
        cc.first_parameter(),
        cc.last_parameter(),
    );
    if !done || sols.is_empty() {
        return None;
    }
    let mut kmin = 0usize;
    let mut dist2_min = sols[0].sq_dist;
    for (k, sol) in sols.iter().enumerate().skip(1) {
        if sol.sq_dist < dist2_min {
            dist2_min = sol.sq_dist;
            kmin = k;
        }
    }
    Some(sols[kmin].u)
}

/// `BRepFeat::ParametricBarycenter(S, CC)` (`BRepFeat.cxx:119-195`): the mean
/// parameter on `CC` of the projected edge samples and vertices.
pub fn parametric_barycenter(shape: &TopoShape, cc: &Arc<dyn Curve>) -> f64 {
    let mut the_map = ShapeSet::new();
    let mut parbar = 0.0f64;
    let mut nbp = 0usize;

    for edg in topo_tools::edges_of(shape) {
        if !the_map.add(&edg.0) {
            continue;
        }
        if BRepTool::is_degenerated(&edg) {
            continue;
        }
        let (f, l) = BRepTool::edge_parameters(&edg);
        let Some(c) = BRepTool::edge_curve_world(&edg) else {
            continue;
        };
        for i in 1..NECHANTBARYC {
            let prm = ((NECHANTBARYC - i) as f64 * f + i as f64 * l) / NECHANTBARYC as f64;
            let pone = c.d0(prm);
            if let Some(prmp) = project_min_param(cc, &pone) {
                nbp += 1;
                parbar += prmp;
            }
        }
    }

    // `// Adds every vertex` (`cxx:169-191`).
    for v in topo_tools::vertices_of(shape) {
        if the_map.add(&v.0) {
            let pone = BRepTool::vertex_point(&v);
            if project_min_param(cc, &pone).is_some() {
                nbp += 1;
            }
        }
    }

    parbar / nbp as f64
}

/// `BRepFeat::ParametricMinMax(...)` (`BRepFeat.cxx:203-332`).
///
/// Returns `(prmin, prmax, prbmin, prbmax, flag)`. `prmin` / `prmax` come from
/// the `LocOpe_CSIntersector` hits on `CC` (`cxx:211-245`); `prbmin` /
/// `prbmax` from the projected edge samples and vertices (`cxx:247-331`).
pub fn parametric_min_max(
    shape: &TopoShape,
    cc: &Arc<dyn Curve>,
    ori: bool,
) -> (f64, f64, f64, f64, bool) {
    let mut asi = CSIntersector::with_shape(shape);
    let _ = asi.perform_curves(&[Some(Arc::clone(cc))]);

    let prmin;
    let prmax;
    let flag;

    if let Some(n) = asi.nb_points(1).ok().filter(|n| *n >= 1) {
        let first = asi.point(1, 1).expect("point 1 of a non-empty hit list");
        let last = asi.point(1, n).expect("point n of a non-empty hit list");
        if !ori {
            prmin = first.parameter().min(last.parameter());
            prmax = first.parameter().max(last.parameter());
        } else if first.orientation() == Orientation::Forward {
            prmin = first.parameter();
            prmax = last.parameter();
        } else {
            prmax = first.parameter();
            prmin = last.parameter();
        }
        flag = true;
    } else {
        prmax = REAL_FIRST;
        prmin = REAL_LAST;
        flag = false;
    }

    let mut the_map = ShapeSet::new();
    let mut prbmin = REAL_LAST;
    let mut prbmax = REAL_FIRST;

    for edg in topo_tools::edges_of(shape) {
        if !the_map.add(&edg.0) {
            continue;
        }
        if BRepTool::is_degenerated(&edg) {
            continue;
        }
        let (f, l) = BRepTool::edge_parameters(&edg);
        let Some(c) = BRepTool::edge_curve_world(&edg) else {
            continue;
        };
        for i in 1..NECHANTBARYC {
            let prm = ((NECHANTBARYC - i) as f64 * f + i as f64 * l) / NECHANTBARYC as f64;
            let pone = c.d0(prm);
            if let Some(prmp) = project_min_param(cc, &pone) {
                if prmp <= prbmin {
                    prbmin = prmp;
                }
                if prmp >= prbmax {
                    prbmax = prmp;
                }
            }
        }
    }

    // `// Adds every vertex` (`cxx:301-331`).
    for v in topo_tools::vertices_of(shape) {
        if the_map.add(&v.0) {
            let pone = BRepTool::vertex_point(&v);
            if let Some(prmp) = project_min_param(cc, &pone) {
                if prmp <= prbmin {
                    prbmin = prmp;
                }
                if prmp >= prbmax {
                    prbmax = prmp;
                }
            }
        }
    }

    (prmin, prmax, prbmin, prbmax, flag)
}

/// `static bool IsIn(BRepTopAdaptor_FClass2d& FC, const Geom2dAdaptor_Curve& AC)`
/// (`BRepFeat.cxx:337-351`): every quasi-uniform sample of `AC` must be
/// `TopAbs_IN` / `TopAbs_ON`, i.e. not `TopAbs_OUT`.
fn is_in(fc: &BrepTopAdaptorFClass2d, ac: &Geom2dAdaptorCurve) -> bool {
    let def = 100.0 * precision::CONFUSION;
    let Some(curve) = ac.curve() else {
        return true;
    };
    let mut qu = QuasiUniformDeflection::new();
    qu.initialize_2d(
        curve.as_ref(),
        ac.first_parameter(),
        ac.last_parameter(),
        def,
        1,
    );
    for i in 1..=qu.nb_points() {
        let Some(param) = qu.parameter(i) else {
            continue;
        };
        let Some(p) = ac.value(param) else {
            continue;
        };
        if fc.perform(p, false) == FaceState::Out {
            return false;
        }
    }
    true
}

/// `static void PutInBoundsU(...)` (`BRepFeat.cxx:364-406`): shifts the 2D
/// curve by whole U periods so its sample span sits inside `[umin, umax]`.
fn put_in_bounds_u(
    umin: f64,
    umax: f64,
    eps: f64,
    period: f64,
    f: f64,
    l: f64,
    c2d: &mut Arc<dyn Curve2d>,
) {
    let pf = c2d.value(f);
    let pl = c2d.value(l);
    let pm = c2d.value(0.34 * f + 0.66 * l);
    let mut min_c = pf.x().min(pl.x()).min(pm.x());
    let mut max_c = pf.x().max(pl.x()).max(pm.x());
    let mut du = 0.0f64;
    if min_c < umin - eps {
        du = (((umin - min_c) / period) as i32 + 1) as f64 * period;
    }
    if min_c > umax + eps {
        du = -((((min_c - umax) / period) as i32) + 1) as f64 * period;
    }
    if du != 0.0 {
        *c2d = translate_2d(c2d, du, 0.0);
        min_c += du;
        max_c += du;
    }
    // `// Ajuste au mieux la courbe dans le domaine.` (`cxx:393-404`).
    if max_c > umax + 100.0 * eps {
        let d1 = max_c - umax;
        let d2 = umin - min_c + period;
        if d2 < d1 {
            du = -period;
        }
        if du != 0.0 {
            *c2d = translate_2d(c2d, du, 0.0);
        }
    }
}

/// `static void PutInBoundsV(...)` (`BRepFeat.cxx:416-462`), the V twin of
/// [`put_in_bounds_u`].
fn put_in_bounds_v(
    vmin: f64,
    vmax: f64,
    eps: f64,
    period: f64,
    f: f64,
    l: f64,
    c2d: &mut Arc<dyn Curve2d>,
) {
    let pf = c2d.value(f);
    let pl = c2d.value(l);
    let pm = c2d.value(0.34 * f + 0.66 * l);
    let mut min_c = pf.y().min(pl.y()).min(pm.y());
    let mut max_c = pf.y().max(pl.y()).max(pm.y());
    let mut dv = 0.0f64;
    if min_c < vmin - eps {
        dv = (((vmin - min_c) / period) as i32 + 1) as f64 * period;
    }
    if min_c > vmax + eps {
        dv = -((((min_c - vmax) / period) as i32) + 1) as f64 * period;
    }
    if dv != 0.0 {
        *c2d = translate_2d(c2d, 0.0, dv);
        min_c += dv;
        max_c += dv;
    }
    if max_c > vmax + 100.0 * eps {
        let d1 = max_c - vmax;
        let d2 = vmin - min_c + period;
        if d2 < d1 {
            dv = -period;
        }
        if dv != 0.0 {
            *c2d = translate_2d(c2d, 0.0, dv);
        }
    }
}

/// `C2d->Translate(gp_Vec2d)` (`Geom2d_Curve::Translate`): the port's curves
/// are immutable handles, so an exact translated copy replaces the handle.
fn translate_2d(c: &Arc<dyn Curve2d>, dx: f64, dy: f64) -> Arc<dyn Curve2d> {
    let mut t = GpTrsf2d::identity();
    t.set_translation_vec(&GpVec2d::new(dx, dy));
    Arc::from(c.transformed(&t))
}

/// `BRepFeat::IsInside(F1, F2)` (`BRepFeat.cxx:467-520`): every edge of `F1`,
/// projected to the surface of `F2`, stays inside the `F2` region.
pub fn is_inside(f1: &Face, f2: &Face) -> bool {
    let Some(s) = BRepTool::face_surface(f2) else {
        // `BRep_Tool::Surface(F2)` answers a null handle in the reference and
        // the first `S->IsUPeriodic()` dereferences it.
        return false;
    };
    let (umin, umax, vmin, vmax) = BRepTool::uv_bounds(f2);

    let mut uperiod = 0.0;
    let mut flagu = false;
    if s.is_u_periodic() {
        flagu = true;
        uperiod = s.u_period();
    }
    let mut vperiod = 0.0;
    let mut flagv = false;
    if s.is_v_periodic() {
        flagv = true;
        vperiod = s.v_period();
    }

    let a_local_shape = f2.0.oriented(Orientation::Forward);
    let fc = BrepTopAdaptorFClass2d::new(&Face(a_local_shape), precision::CONFUSION);

    for edge in topo_tools::edges_of(&f1.0) {
        let (f1p, l1p) = BRepTool::edge_parameters(&edge);
        let Some(c0) = BRepTool::edge_curve(&edge) else {
            continue;
        };
        let mut tol = precision::CONFUSION;
        let Some(c) = curve_2d_range(&c0, f1p, l1p, &s, &mut tol) else {
            continue;
        };
        let mut c = c;
        if flagu || flagv {
            let eps = BRepTool::edge_tolerance(&edge);
            if flagu {
                put_in_bounds_u(umin, umax, eps, uperiod, f1p, l1p, &mut c);
            }
            if flagv {
                put_in_bounds_v(vmin, vmax, eps, vperiod, f1p, l1p, &mut c);
            }
        }
        let ac = Geom2dAdaptorCurve::with_range(c, f1p, l1p);
        if !is_in(&fc, &ac) {
            return false;
        }
    }
    true
}

/// `BRepFeat::IsInOut(FC, AC)` (`BRepFeat.hxx:105-106`) — declared in the
/// header, defined as the file-local `IsIn` above; exposed under its public
/// name.
pub fn is_in_out(fc: &BrepTopAdaptorFClass2d, ac: &Geom2dAdaptorCurve) -> bool {
    is_in(fc, ac)
}

/// `BRepFeat::FaceUntil(Sbase, FUntil)` (`BRepFeat.cxx:524-637`): extends a
/// single-face `FUntil` on a plane / cylinder / cone so it covers `Sbase` plus
/// a `10 * sqrt(square extent)` margin, and answers `None` for the
/// reference's `FUntil.Nullify()` arms.
pub fn face_until(sbase: &TopoShape, f_until: &Face) -> Option<Face> {
    let mut b = BndBox::new();
    brep_bnd_lib::add(sbase, &mut b, false);
    let (x, y, z) = match b.get() {
        Some((xmin, ymin, zmin, xmax, ymax, zmax)) => {
            ([xmin, xmax], [ymin, ymax], [zmin, zmax])
        }
        // `B.Get(...)` on a void box raises in the reference; no margin can be
        // derived, so the face is nullified like the non-analytic arms.
        None => return None,
    };
    let diam = 10.0 * b.square_extent().sqrt();

    let Some(s) = BRepTool::face_surface(f_until) else {
        return None;
    };
    let s = s.rectangular_trimmed_basis().unwrap_or(s);

    let str: Option<GeomRectangularTrimmedSurface> = if let Some(pln) = s.gp_pln() {
        let (mut umin, mut umax) = (REAL_LAST, -REAL_LAST);
        let (mut vmin, mut vmax) = (REAL_LAST, -REAL_LAST);
        for xi in x {
            for yj in y {
                for zk in z {
                    let (u, v) = slib::plane_parameters(&pln.position(), &GpPnt::new(xi, yj, zk));
                    umin = umin.min(u);
                    umax = umax.max(u);
                    vmin = vmin.min(v);
                    vmax = vmax.max(v);
                }
            }
        }
        umin -= diam;
        umax += diam;
        vmin -= diam;
        vmax += diam;
        GeomRectangularTrimmedSurface::uv_sense(s.clone(), umin, umax, vmin, vmax, true, true)
            .ok()
    } else if let Some(cyl) = s.gp_cylinder() {
        let (mut vmin, mut vmax) = (REAL_LAST, -REAL_LAST);
        for xi in x {
            for yj in y {
                for zk in z {
                    let (_u, v) =
                        slib::cylinder_parameters(&cyl.position(), &GpPnt::new(xi, yj, zk));
                    vmin = vmin.min(v);
                    vmax = vmax.max(v);
                }
            }
        }
        vmin -= diam;
        vmax += diam;
        GeomRectangularTrimmedSurface::one_param(s.clone(), vmin, vmax, false).ok()
    } else if let Some(cone) = s.gp_cone() {
        let (mut vmin, mut vmax) = (REAL_LAST, -REAL_LAST);
        for xi in x {
            for yj in y {
                for zk in z {
                    let (_u, v) = slib::cone_parameters(
                        &cone.position(),
                        cone.radius(),
                        cone.semi_angle(),
                        &GpPnt::new(xi, yj, zk),
                    );
                    vmin = vmin.min(v);
                    vmax = vmax.max(v);
                }
            }
        }
        vmin -= diam;
        vmax += diam;
        GeomRectangularTrimmedSurface::one_param(s.clone(), vmin, vmax, false).ok()
    } else {
        None
    };

    let str = str?;
    Some(brep_lib_make_face::make_face_from_surface(
        Arc::new(str),
        precision::CONFUSION,
    ))
}

/// `BRepFeat::Tool(SRef, Fac, Orf)` (`BRepFeat.cxx:643-712`): the solid of the
/// `SRef` faces oriented so that the face `Fac` carries the orientation `Orf`.
///
/// Returns `None` for the reference's null-solid arms: a non-shell / non-solid
/// `LocOpe_BuildShape` result, or a solid without exactly one shell.
pub fn tool(sref: &TopoShape, fac: &Face, orf: Orientation) -> Option<Solid> {
    let faces: Vec<TopoShape> = topo_tools::faces_of(sref)
        .into_iter()
        .map(|f| f.0)
        .collect();

    let mut bs = BuildShape::new();
    if bs.with_faces(&faces).is_err() {
        return None;
    }
    let res = bs.shape()?.clone();

    let mut sh = match res.shape_type() {
        ShapeType::Shell => Some(Shell(res.clone())),
        ShapeType::Solid => {
            let shells = topo_tools::shapes_of(&res, ShapeType::Shell);
            // `exp.Next(); if (exp.More()) Sh.Nullify();` (`cxx:669-673`).
            if shells.len() == 1 {
                Some(Shell(shells[0].clone()))
            } else {
                None
            }
        }
        _ => None,
    }?;

    // `Sh.Orientation(TopAbs_FORWARD);` (`cxx:683`).
    sh.0.set_orientation(Orientation::Forward);

    let mut orient = Orientation::Forward;
    for f in topo_tools::faces_of(&sh.0) {
        if topo_tools::is_same(&f.0, &fac.0) {
            orient = f.0.orientation();
            break;
        }
    }

    let reverse = (orient == fac.0.orientation() && orf == Orientation::Reversed)
        || (orient != fac.0.orientation() && orf == Orientation::Forward);
    if reverse {
        sh.0.reverse();
    }

    Some(TopoBuilder::new().make_solid(&[sh]))
}
