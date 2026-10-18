//! `BRepLib::SameParameter` for a shape (forced, mutable).
//! Source: `BRepLib.cxx` `InternalSameParameter` 913-1008, edge
//! `SameParameter` 1251-1739, `ComputeTol` 1070-1188.

use std::collections::HashSet;

use occt_core::precision::{CONFUSION, INFINITE, PCONFUSION, Precision};
use occt_core::kernel::geomabs::Shape;
use occt_geom::approx_same_parameter::{u_resolution, v_resolution};
use occt_geom::geom_lib::knots_reparametrize;
use occt_geom::{ApproxSameParameter, Curve, CurvilinearParameter, Surface};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::c0_bspline_to_c1_bspline_curve;
use std::sync::Arc;

use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};
use crate::tgeometry::{GeometryRegistry, VertexGeom};
use crate::topo_tools_full::{edge_vertices, edges_of, faces_of};

const NCONTROL: i32 = 22;

/// `EvalTol` (`BRepLib.cxx:1034-1066`): probe five parameters on the pcurve,
/// project the surface points onto the 3d curve with `Extrema_LocateExtPC`,
/// and widen `tolbail` to the reached distance. Returns `true` when more than
/// two projections succeed.
fn eval_tol(
    pc: &dyn Curve2d,
    surf: &dyn Surface,
    c3d: &dyn Curve,
    first: f64,
    last: f64,
    tol: f64,
) -> Option<f64> {
    let mut tolbail = tol;
    let mut ok = 0i32;
    for i in 1..=5 {
        let t_frac = i as f64 / 6.0;
        let t = (1.0 - t_frac) * first + t_frac * last;
        let puv = pc.d0(t);
        let p = surf.d0(puv.x(), puv.y());
        if let Some((_tp, q)) =
            crate::int_tools_vertex_line::extrema_locate_ext_pc(c3d, &p, t, first, last)
        {
            let dist2 = p.square_distance(&q);
            if dist2 > tolbail * tolbail {
                tolbail = dist2.sqrt();
            }
            ok += 1;
        }
    }
    if ok > 2 {
        Some(tolbail)
    } else {
        None
    }
}

/// Periodic `SetOrigin` realignment after `C0BSplineToC1BSplineCurve`
/// (`BRepLib.cxx:1419-1441`, `:1471-1493`).
fn realign_periodic_origin(bs2d: &mut occt_geom2d::Geom2dBSplineCurve, origin_point: &occt_core::gp::GpPnt2d) {
    if !bs2d.is_periodic() {
        return;
    }
    let new_origin = bs2d.d0(bs2d.first_parameter());
    if (origin_point.x() - new_origin.x()).abs() <= PCONFUSION
        && (origin_point.y() - new_origin.y()).abs() <= PCONFUSION
    {
        return;
    }
    let (knots, _) = bs2d.distinct_knots_and_mults();
    for (idx, &k) in knots.iter().enumerate() {
        let p = bs2d.d0(k);
        if (origin_point.x() - p.x()).abs() > PCONFUSION
            || (origin_point.y() - p.y()).abs() > PCONFUSION
        {
            continue;
        }
        let _ = bs2d.set_origin((idx + 1) as i32);
        break;
    }
}

/// Knot-span ratio "bad BSpline" test (`BRepLib.cxx:1537-1581`).
fn is_bad_bspline_knots(bs2d: &occt_geom2d::Geom2dBSplineCurve, error: f64, tol: f64) -> bool {
    let cont = bs2d.continuity();
    if cont <= 0 || error <= (1.0e-3_f64).max(tol) {
        return false;
    }
    let (knots, _) = bs2d.distinct_knots_and_mults();
    let nb_knots = knots.len();
    if nb_knots < 3 {
        return false;
    }
    const CRIT_RATIO: f64 = 10.0;
    let mut dtprev = knots[1] - knots[0];
    let mut dtmin = dtprev;
    let mut is_bad = false;
    for j in 1..nb_knots - 1 {
        let dtcur = knots[j + 1] - knots[j];
        dtmin = dtmin.min(dtcur);
        if !is_bad {
            let dtratio = if dtcur > dtprev {
                dtcur / dtprev
            } else if dtcur > 0.0 {
                dtprev / dtcur
            } else {
                f64::INFINITY
            };
            if dtratio > CRIT_RATIO {
                is_bad = true;
            }
        }
        dtprev = dtcur;
    }
    if is_bad {
        // `cxx:1575-1579`: clear IsBad when the smallest span is still above
        // the curve resolution at `max(1e-3, theTolerance)`.
        let dtcur = bs2d.resolution((1.0e-3_f64).max(tol));
        if dtmin < dtcur {
            is_bad = false;
        }
    }
    is_bad
}

/// `BRepLib.cxx:1070-1188`.
fn compute_tol(
    c3d: &dyn Curve,
    c2d: &dyn Curve2d,
    surf: &dyn Surface,
    first: f64,
    last: f64,
    nbp: i32,
) -> f64 {
    let nbp = nbp.max(1) as usize;
    let mut dist = vec![-1.0; nbp + 10];
    let (uf, ul) = surf.u_range();
    let (vf, vl) = surf.v_range();
    let du = 0.01 * (ul - uf);
    let dv = 0.01 * (vl - vf);
    let is_u_per = surf.is_u_periodic();
    let is_v_per = surf.is_v_periodic();
    let dsdu = 1.0 / u_resolution(surf, 1.0).max(1.0e-30);
    let dsdv = 1.0 / v_resolution(surf, 1.0).max(1.0e-30);
    let mut d2: f64 = 0.0;
    let mut dapp: f64 = -1.0;
    for i in 0..=nbp {
        let t = i as f64 / nbp as f64;
        let u = first * (1.0 - t) + last * t;
        let pc3d = c3d.d0(u);
        let puv = c2d.d0(u);
        if !is_u_per {
            if puv.x() < uf - du {
                dapp = dapp.max(dsdu * (uf - puv.x()));
                continue;
            } else if puv.x() > ul + du {
                dapp = dapp.max(dsdu * (puv.x() - ul));
                continue;
            }
        }
        if !is_v_per {
            if puv.y() < vf - dv {
                dapp = dapp.max(dsdv * (vf - puv.y()));
                continue;
            } else if puv.y() > vl + dv {
                dapp = dapp.max(dsdv * (puv.y() - vl));
                continue;
            }
        }
        let pcons = surf.d0(puv.x(), puv.y());
        if Precision::is_infinite(pcons.x())
            || Precision::is_infinite(pcons.y())
            || Precision::is_infinite(pcons.z())
        {
            return INFINITE;
        }
        let temp = pc3d.square_distance(&pcons);
        if i + 1 < dist.len() {
            dist[i + 1] = temp;
        }
        d2 = d2.max(temp);
    }
    if Precision::is_infinite(d2) {
        return d2;
    }
    d2 = d2.sqrt();
    if dapp > d2 {
        return dapp;
    }
    let mut n1 = 0i32;
    let mut n2 = 0i32;
    for &di in dist.iter().skip(1) {
        if di > 0.0 {
            if di < 1.0 {
                n1 += 1;
            } else {
                n2 += 1;
            }
        }
    }
    let mut n3 = 0i32;
    if n1 > n2 && n2 != 0 {
        n3 = 100 * n2 / (n1 + n2);
    }
    let mut ana = false;
    let mut d2_ana: f64 = 0.0;
    if n3 < 10 && n3 != 0 {
        ana = true;
        for &di in dist.iter().skip(1) {
            if di > 0.0 && di < 1.0 {
                d2_ana = d2_ana.max(di);
            }
        }
    }
    d2 = if !ana { 1.5 * d2 } else { 1.5 * d2_ana.sqrt() };
    d2.max(1.0e-7)
}

/// `BRepLib::SameParameter` four-argument core (`BRepLib.cxx:1251-1739`):
/// walk every `CurveOnSurface` representation (`cxx:1263-1291`, `:1352-1718`).
fn same_parameter_edge(edge: &Edge, tol: f64) -> f64 {
    let reg = GeometryRegistry::global();
    let Some(mut g) = reg.edge_geom(&edge.0) else {
        return -1.0;
    };
    if g.same_parameter {
        return -1.0;
    }
    if g.degenerated {
        return -1.0;
    }
    let c3d = g.curve.clone();
    let mut f3d = g.first;
    let mut l3d = g.last;
    if !c3d.is_periodic() {
        let udeb = c3d.first_parameter();
        let ufin = c3d.last_parameter();
        if udeb > f3d {
            f3d = udeb;
        }
        if l3d > ufin {
            l3d = ufin;
        }
    }
    // `GetCurve3d` Cos list (`cxx:1263-1291`): every stored representation.
    let reps = reg.edge_pcurve_reps(&edge.0);
    if reps.is_empty() {
        g.first = f3d;
        g.last = l3d;
        g.same_range = true;
        reg.set_edge(&edge.0, g);
        return -1.0;
    }
    let mut is_same_p = true;
    let mut maxdist = 0.0;
    let an_edge_tol = g.tolerance;
    let same_range = g.same_range;
    let mut ya_pcu = false;
    const BIG_ERROR: f64 = 1.0e10;
    // `cxx:1378`: TolSameRange = max(GAC.Resolution(theTolerance), PConfusion).
    // Passed to `GeomLib::SameRange` (`cxx:1389`, `:1678`) and used by the
    // Curvilinear reparam range check (`cxx:1610-1611`).
    let tol_same_range = crate::int_tools_vertex_line::adaptor_resolution(c3d.as_ref(), tol)
        .max(occt_core::precision::PCONFUSION);
    // Group Cos by representation key so seam pairs stay together when rewritten.
    let mut by_key: std::collections::BTreeMap<
        usize,
        (std::sync::Arc<dyn Surface>, f64, f64, Vec<std::sync::Arc<dyn Curve2d>>),
    > = std::collections::BTreeMap::new();
    for (surf, pc, first_on, last_on) in reps {
        let fk = std::sync::Arc::as_ptr(&surf) as *const () as usize;
        by_key
            .entry(fk)
            .and_modify(|e| e.3.push(pc.clone()))
            .or_insert_with(|| (surf, first_on, last_on, vec![pc]));
    }
    for (fk, (surf, first_on, last_on, pcs)) in by_key {
        let mut out_pcs: Vec<std::sync::Arc<dyn Curve2d>> = Vec::with_capacity(pcs.len());
        let mut pcs_changed = false;
        for pc in &pcs {
            ya_pcu = true;
            // `cxx:1387-1391`: `GeomLib::SameRange` when `!SameRange`.
            let mut cur_pc = if !same_range {
                let remapped = crate::shhealing::geom_lib_same_range(
                    tol_same_range,
                    pc.clone(),
                    first_on,
                    last_on,
                    f3d,
                    l3d,
                );
                if !std::sync::Arc::ptr_eq(&remapped, pc) {
                    pcs_changed = true;
                }
                remapped
            } else {
                pc.clone()
            };
            let mut error =
                compute_tol(c3d.as_ref(), cur_pc.as_ref(), surf.as_ref(), f3d, l3d, NCONTROL);
            if error > BIG_ERROR {
                maxdist = error;
                out_pcs.push(cur_pc);
                break;
            }
            // `cxx:1404-1622`: C0 BSpline pcurve -> C1 via
            // `Geom2dConvert::C0BSplineToC1BSplineCurve`.
            let mut goodpc = true;
            let mut is_bsp = false;
            let mut is_ana = false;
            if cur_pc.is_bspline2d() && cur_pc.continuity() == 0 {
                if let Some(mut bs2d) = cur_pc.bspline_copy2d() {
                    let u_resol = u_resolution(surf.as_ref(), tol);
                    let v_resol = v_resolution(surf.as_ref(), tol);
                    let tol_conf2d = u_resol.min(v_resol).max(PCONFUSION);
                    let bs2dsov = bs2d.clone();
                    let f_c0 = bs2d.first_parameter();
                    let l_c0 = bs2d.last_parameter();
                    let mut repar = true;
                    let origin_point = bs2d.d0(f_c0);
                    if c0_bspline_to_c1_bspline_curve(&mut bs2d, tol_conf2d).is_ok() {
                        is_bsp = true;
                        realign_periodic_origin(&mut bs2d, &origin_point);
                        if bs2d.continuity() == 0 {
                            // `cxx:1444-1506`: EvalTol then Tol2dbail retry.
                            if let Some(tolbail) = eval_tol(
                                cur_pc.as_ref(),
                                surf.as_ref(),
                                c3d.as_ref(),
                                f3d,
                                l3d,
                                tol,
                            ) {
                                bs2d = bs2dsov.clone();
                                let u_resbail = u_resolution(surf.as_ref(), tolbail);
                                let v_resbail = v_resolution(surf.as_ref(), tolbail);
                                let mut tol2dbail = u_resbail.min(v_resbail);
                                let origin_point = bs2d.d0(bs2d.first_parameter());
                                let nbp = bs2d.nb_poles();
                                let mut d = INFINITE;
                                if nbp >= 2 {
                                    if let Some(poles) = bs2d.poles2d() {
                                        let mut p = poles[0];
                                        for ip in 1..nbp {
                                            let p1 = poles[ip];
                                            d = d.min(p.square_distance(&p1));
                                            p = p1;
                                        }
                                        d = d.sqrt() * 0.1;
                                    }
                                }
                                tol2dbail = tol2dbail.min(d).max(tol_conf2d);
                                let _ = c0_bspline_to_c1_bspline_curve(&mut bs2d, tol2dbail);
                                realign_periodic_origin(&mut bs2d, &origin_point);
                                if bs2d.continuity() == 0 {
                                    goodpc = true;
                                    bs2d = bs2dsov.clone();
                                    repar = false;
                                }
                            } else {
                                goodpc = false;
                            }
                        }
                        if goodpc {
                            if repar {
                                let (mut knots, _) = bs2d.distinct_knots_and_mults();
                                knots_reparametrize(f_c0, l_c0, &mut knots);
                                let _ = bs2d.set_knots(&knots);
                                let new_pc: Arc<dyn Curve2d> = Arc::new(bs2d.clone());
                                let error1 = compute_tol(
                                    c3d.as_ref(),
                                    new_pc.as_ref(),
                                    surf.as_ref(),
                                    f3d,
                                    l3d,
                                    NCONTROL,
                                );
                                if error1 > error {
                                    // Keep `cur_pc` as the pre-convert curve
                                    // (`cxx:1525-1529`).
                                    bs2d = bs2dsov.clone();
                                    is_ana = true;
                                } else {
                                    error = error1;
                                    cur_pc = new_pc;
                                    pcs_changed = true;
                                }
                            }
                            // `cxx:1537-1621`: knot-span ratio + Curvilinear.
                            if is_bad_bspline_knots(&bs2d, error, tol) {
                                let mut cont = match bs2d.continuity() {
                                    0 => Shape::C0,
                                    1 => Shape::G1,
                                    2 => Shape::C1,
                                    3 => Shape::G2,
                                    4 => Shape::C2,
                                    5 => Shape::C3,
                                    _ => Shape::CN,
                                };
                                if cont > Shape::C2 {
                                    cont = Shape::C2;
                                }
                                let mut maxdeg = bs2d.degree() as i32;
                                if maxdeg == 1 {
                                    maxdeg = 14;
                                }
                                let pc_arc: Arc<dyn Curve2d> = Arc::new(bs2d.clone());
                                let surf_arc = surf.clone();
                                let cos: Arc<dyn Curve> = Arc::new(
                                    crate::meshing::edge_discret::CurveOnSurface::new(
                                        pc_arc.clone(),
                                        surf_arc.clone(),
                                        f3d,
                                        l3d,
                                    ),
                                );
                                let app = CurvilinearParameter::from_curve_on_surface(
                                    pc_arc,
                                    surf_arc.clone(),
                                    cos,
                                    (1.0e-3_f64).max(tol),
                                    cont,
                                    maxdeg,
                                    10,
                                    |c2d, s| {
                                        Arc::new(
                                            crate::meshing::edge_discret::CurveOnSurface::new(
                                                c2d, s, f3d, l3d,
                                            ),
                                        )
                                    },
                                );
                                if (app.done || app.has_result) && app.curve2d1.is_some() {
                                    let mut new_bs = app.curve2d1.unwrap();
                                    if (new_bs.first_parameter() - f_c0).abs() > tol_same_range
                                        || (new_bs.last_parameter() - l_c0).abs() > tol_same_range
                                    {
                                        let (mut knots, _) = new_bs.distinct_knots_and_mults();
                                        knots_reparametrize(f_c0, l_c0, &mut knots);
                                        let _ = new_bs.set_knots(&knots);
                                    }
                                    bs2d = new_bs;
                                    cur_pc = Arc::new(bs2d.clone());
                                    pcs_changed = true;
                                }
                            }
                        }
                    }
                }
            }
            if !goodpc {
                out_pcs.push(cur_pc);
                continue;
            }
            let tol_val = if is_ana && is_bsp { 1.0e-7 } else { tol };
            let same_p = ApproxSameParameter::new(
                c3d.as_ref(),
                cur_pc.as_ref(),
                surf.as_ref(),
                f3d,
                l3d,
                tol_val,
            );
            if same_p.same_parameter {
                maxdist = maxdist.max(same_p.tol_reached);
            } else if same_p.done {
                maxdist = maxdist.max(same_p.tol_reached.min(error));
                // `cxx:1648-1671`: replace pcurve when Approx rebuilt it.
                if let Some(ref new_pc) = same_p.curve2d {
                    if same_p.tol_reached <= error {
                        cur_pc = new_pc.clone();
                        pcs_changed = true;
                    }
                }
            } else {
                // `cxx:1678-1684`: Approx failed - still SameRange onto 3d domain.
                let remapped = crate::shhealing::geom_lib_same_range(
                    tol_same_range,
                    pc.clone(),
                    first_on,
                    last_on,
                    f3d,
                    l3d,
                );
                if !std::sync::Arc::ptr_eq(&remapped, &cur_pc) {
                    pcs_changed = true;
                }
                cur_pc = remapped;
                is_same_p = false;
            }
            if !is_same_p {
                // `BRepCheck::PrecCurve` / `PrecSurface` (`BRepCheck.cxx:70-129`):
                // ellipse/cone epsilon specials are not ported; default is `RealEpsilon`.
                let prec = f64::EPSILON;
                let cur_tol = an_edge_tol + prec;
                if cur_tol >= error {
                    maxdist = maxdist.max(an_edge_tol);
                    is_same_p = true;
                }
            }
            out_pcs.push(cur_pc);
        }
        if pcs_changed && !out_pcs.is_empty() {
            // `fk` is the surface-pointer Cos key; `repr_key` leaves it unchanged
            // when it is not a face shape key (`tgeometry.rs:587-593`).
            reg.set_edge_pcurves(&edge.0, fk, out_pcs);
            reg.set_pcurve_range(&edge.0, fk, f3d, l3d);
        }
    }
    g.first = f3d;
    g.last = l3d;
    g.same_range = true;
    let mut new_tol = -1.0;
    if is_same_p {
        if ya_pcu {
            maxdist = maxdist.max(CONFUSION);
            new_tol = maxdist;
            g.tolerance = maxdist;
        }
        g.same_parameter = true;
    }
    reg.set_edge(&edge.0, g);
    new_tol
}

fn update_v_tol(edge: &Edge, new_tol: f64) {
    if new_tol <= 0.0 {
        return;
    }
    let reg = GeometryRegistry::global();
    let (v1, v2) = edge_vertices(edge);
    for v in [v1, v2].into_iter().flatten() {
        let mut vg = reg.vertex_geom(&v.0).unwrap_or(VertexGeom {
            point: Default::default(),
            tolerance: 0.0,
        });
        if new_tol > vg.tolerance {
            vg.tolerance = new_tol;
            reg.set_vertex(&v.0, vg);
        }
    }
}

/// `BRepLib::SameParameter(theEdge, theTolerance)` (`BRepLib.cxx:1237-1246`):
/// the four argument overload with `IsUseOldEdge = true` (the edge is modified
/// in place) followed by `UpdateVTol` (`cxx:1242-1245`). Walks every Cos
/// representation (`cxx:1263-1291`).
pub fn same_parameter_edge_inplace(edge: &Edge, tolerance: f64) {
    let new_tol = same_parameter_edge(edge, tolerance);
    update_v_tol(edge, new_tol);
}

/// Face-scoped entry used by callers that already hold a face; Cos walk is
/// still every representation on the edge (`BRepLib.cxx:1251-1739`).
pub fn same_parameter_edge_inplace_on_face(edge: &Edge, _face: &Face, tolerance: f64) {
    same_parameter_edge_inplace(edge, tolerance);
}

fn update_tolerances(shape: &TopoShape) {
    let reg = GeometryRegistry::global();
    for face in faces_of(shape) {
        let ft = BRepTool::face_surface(&face)
            .map(|_| reg.face_tolerance(&face.0))
            .unwrap_or(0.0);
        for e in edges_of(&face.0) {
            if let Some(mut g) = reg.edge_geom(&e.0) {
                if g.tolerance < ft {
                    g.tolerance = ft;
                    reg.set_edge(&e.0, g);
                }
            }
            let et = reg.edge_geom(&e.0).map(|g| g.tolerance).unwrap_or(0.0);
            let (v1, v2) = edge_vertices(&e);
            for v in [v1, v2].into_iter().flatten() {
                let mut vg = reg.vertex_geom(&v.0).unwrap_or(VertexGeom {
                    point: Default::default(),
                    tolerance: 0.0,
                });
                if vg.tolerance < et {
                    vg.tolerance = et;
                    reg.set_vertex(&v.0, vg);
                }
            }
        }
    }
}

/// `BRepLib::SameParameter(S, Tolerance, forced)` (`cxx:1016`).
pub fn same_parameter_shape(shape: &TopoShape, tolerance: f64, forced: bool) {
    let reg = GeometryRegistry::global();
    let mut done: HashSet<usize> = HashSet::new();
    for e in edges_of(shape) {
        let key = GeometryRegistry::shape_key(&e.0);
        if !done.insert(key) {
            continue;
        }
        if let Some(mut g) = reg.edge_geom(&e.0) {
            if forced && (g.same_range || g.same_parameter) {
                g.same_range = false;
                g.same_parameter = false;
                reg.set_edge(&e.0, g);
            }
        }
        // `cxx:1251-1739`: every Cos on the edge, not one call-site face.
        let new_tol = same_parameter_edge(&e, tolerance);
        update_v_tol(&e, new_tol);
    }
    // Unported: planar `GetEdgeTol` (`cxx:973-1002`). Offset faces skip it.
    update_tolerances(shape);
}

/// `BRepLib::SameParameter(F, Confusion, true)` from `BRepLib_MakeFace.cxx:865`.
pub fn same_parameter_face(face: &Face, tolerance: f64, forced: bool) {
    same_parameter_shape(&face.0, tolerance, forced);
}
