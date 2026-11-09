//! `BRepOffset_Inter2d::ExtentEdge` (`BRepOffset_Inter2d.cxx:1193-1676`) with
//! the two file-local helpers it drives:
//!
//! * `EvaluateMaxSegment` (`:1078-1097`) - [`evaluate_max_segment`]
//! * `ExtendPCurve` (`:1100-1187`) - [`extend_pcurve`]
//!
//! The reference works on `NE = E.EmptyCopied()` (`cxx:1197`, `:1207`), which
//! `TopoDS_Shape::EmptyCopy` (`TopoDS_Shape.hxx:292`) defines as "the same
//! geometry and no sub-shapes": every `BRep_CurveRepresentation` of `E` comes
//! along. `ExtentEdge` then, for every CurveOnSurface representation:
//!
//! 1. prolongs a bounded pcurve (`ExtendPCurve`) or widens a periodic / closed
//!    one (`:1231-1262`),
//! 2. clips the resulting `[First, Last]` back into the surface's parametric
//!    rectangle by intersecting the pcurve with the four boundary lines
//!    (`:1265-1379`) and writes it with `BRep_GCurve::SetRange` (`:1381`),
//! 3. keeps the representation with the largest `FirstPar` and the one with the
//!    smallest `LastPar` (`:1386-1401`).
//!
//! Afterwards the 3D curve is rebuilt (`:1404-1603`): a present `Curve3D`
//! keeps its type and only its range is repaired (through
//! `GeomAPI_ProjectPointOnCurve` for the bounded case), a missing one is built
//! from the surviving pcurve (`GeomLib::BuildCurve3d`, or a straight
//! `Geom_Line` when both pcurve bounds are infinite), and every other pcurve is
//! re-projected onto its surface with `GeomProjLib::Curve2d`. An edge without
//! any pcurve is prolonged in 3D instead (`:1605-1672`).
//!
//! Port mapping of the representation list:
//!
//! * `BRep_TEdge::ChangeCurves()` holds one `BRep_CurveOnSurface` per
//!   `Geom_Surface` (a `BRep_CurveOnClosedSurface` carries two pcurves); the
//!   port stores the pcurves of one surface in one slot
//!   (`EdgePcurves::curves`, `tshape.rs:635-681`). A "representation" here is
//!   one surface slot: `NbPCurves` (`cxx:1226`) counts slots and `PCurve2`
//!   (`cxx:1245`) is `curves[1]`.
//! * `CurveRep->Location()` is the identity location - the port keys a slot by
//!   the surface pointer and keeps no per-representation location
//!   (`GeometryRegistry::repr_key`, `tgeometry/registry_faces.rs:28-37`), so
//!   `MinLoc` (`cxx:1406`) is the edge's own location.
//! * `CurveRep->SetRange` / `BRep_Builder::Range(edge, first, last)` touch the
//!   slot's `ranges` entry exactly as `BRep_GCurve::SetRange` does; the
//!   `Only3d = false` overload is [`GeometryRegistry::set_ranges_all`].
//! * `BRep_TEdge::ChangeCurves()` is insertion ordered while the port's slots
//!   are a map, so the representation order is pinned by surface address. It
//!   only matters for the `MinPC` tie-break below, i.e. when two
//!   representations answer the same `[First, Last]`.
//!
//! Port guards (no OCCT branch changed):
//!
//! * `BRep_Tool::Curve` answers null in the "no pcurves" arm (`cxx:1607`
//!   dereferences it): the edge is returned unchanged.
//! * `gp_Dir` / `gp_Dir2d` construction from a zero derivative raises in OCCT;
//!   the arm answers `false` instead.

use std::sync::Arc;

use occt_core::gp::{GpDir, GpDir2d, GpPnt2d, GpVec2d};
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::{Precision, ANGULAR, CONFUSION, PCONFUSION};
use occt_geom::geom_lib::geom_lib_build_curve3d;
use occt_geom::geom_proj_lib::curve_2d_range;
use occt_geom::geom_api::project_point_on_curve;
use occt_geom::{Curve, GeomLine, Surface};
use occt_geom2d::comp_curve_to_bspline::CompCurveToBSplineCurve2d;
use occt_geom2d::geom2d_int::Geom2dIntGInter;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;
use occt_geom2d::{Curve2d, Geom2dLine};

use occt_core::convert::ParameterisationType;

use crate::brep_lib_same_parameter::same_parameter_edge_inplace;
use crate::brep_tool::BRepTool;
use crate::meshing::edge_discret::CurveOnSurface;
use crate::shape::{Edge, TopoShape};
use crate::tgeometry::GeometryRegistry;

/// `RealLast()` (`Standard_Real.hxx`), the initial `LastParOnPC`.
const REAL_LAST: f64 = f64::MAX;
/// `RealFirst()` (`Standard_Real.hxx`), the initial `FirstParOnPC`.
const REAL_FIRST: f64 = -f64::MAX;

/// `IsKind(STANDARD_TYPE(Geom2d_BoundedCurve))`: `Geom2d_BoundedCurve.hxx`
/// documents exactly three concrete subclasses - `Geom2d_BezierCurve`,
/// `Geom2d_BSplineCurve` and `Geom2d_TrimmedCurve` - and the port's
/// `is_bezier2d` / `is_bspline2d` / `trimmed_basis` are the matching RTTI
/// queries (`occt-geom2d/src/curve.rs:84-114`, `:170-173`).
fn is_bounded_curve2d(c: &dyn Curve2d) -> bool {
    c.trimmed_basis().is_some() || c.is_bezier2d() || c.is_bspline2d()
}

/// `IsKind(STANDARD_TYPE(Geom_BoundedCurve))`: `Geom_BoundedCurve.hxx`
/// documents `Geom_BezierCurve`, `Geom_BSplineCurve` and
/// `Geom_TrimmedCurve` as its concrete subclasses.
fn is_bounded_curve3d(c: &dyn Curve) -> bool {
    c.is_geom_trimmed() || (!c.is_geom_trimmed() && c.bspline_poles().is_some()) || c.bezier_poles().is_some()
}

/// `IsKind(STANDARD_TYPE(Geom_BoundedSurface))`: `Geom_BoundedSurface.hxx`
/// documents `Geom_BezierSurface`, `Geom_BSplineSurface` and
/// `Geom_RectangularTrimmedSurface` (`Geom_BoundedSurface.hxx:41-49`); an
/// elementary surface derives from `Geom_ElementarySurface`, not from this
/// class (`:50`).
fn is_bounded_surface(s: &dyn Surface) -> bool {
    s.rectangular_trimmed_basis().is_some()
        || s.is_bspline_surface()
        || s.osculating_bspline().is_some()
}

/// `evaluateMaxSegment` (`cxx:1078-1097`): `30 + max(surface knots, pcurve
/// knots)`.
fn evaluate_max_segment(surface: &Arc<dyn Surface>, pcurve: &Arc<dyn Curve2d>) -> i32 {
    let mut nb_s_knots = 0usize;
    if surface.is_bspline_surface() {
        let nu = surface.bspline_surface_uknots().map_or(0, |k| k.len());
        let nv = surface.bspline_surface_vknots().map_or(0, |k| k.len());
        nb_s_knots = nu.max(nv);
    }
    let nb_c2d_knots = pcurve.bspline_nb_knots().unwrap_or(0);
    (30 + nb_s_knots.max(nb_c2d_knots)) as i32
}

/// The straight line through the poles of a two-pole Bezier / B-spline
/// (`cxx:1119-1138`).
fn line_through_poles(poles: &[GpPnt2d]) -> Option<Arc<dyn Curve2d>> {
    let dir = GpDir2d::from_vec2d(&GpVec2d::new(poles[1].x() - poles[0].x(), poles[1].y() - poles[0].y())).ok()?;
    Some(Arc::new(Geom2dLine::from_pnt_dir(poles[0], dir)))
}

/// `ExtendPCurve(PCurve, anEf, anEl, a2Offset, NewPCurve)` (`cxx:1100-1187`).
///
/// `None` is the reference's `false` return; `Some(NewPCurve)` is `true`.
fn extend_pcurve(
    pcurve: &Arc<dyn Curve2d>,
    an_ef: f64,
    an_el: f64,
    a2_offset: f64,
) -> Option<Arc<dyn Curve2d>> {
    // `NewPCurve = aPCurve;` then, for a trimmed curve, `BasisCurve()`
    // (`cxx:1105-1110`).
    let new_pcurve: Arc<dyn Curve2d> = match pcurve.trimmed_basis() {
        Some(b) => Arc::from(b.clone_dyn()),
        None => pcurve.clone(),
    };

    let first_par = new_pcurve.first_parameter();
    let last_par = new_pcurve.last_parameter();

    // `cxx:1113-1140`: a bounded curve that stops short of the requested
    // range and is a straight segment degenerates into an infinite line.
    if is_bounded_curve2d(new_pcurve.as_ref())
        && (first_par > an_ef - a2_offset || last_par < an_el + a2_offset)
    {
        if new_pcurve.is_bezier2d() {
            if let Some(poles) = new_pcurve.poles2d() {
                if poles.len() == 2 {
                    return line_through_poles(&poles);
                }
            }
        } else if new_pcurve.is_bspline2d() {
            let nb_knots = new_pcurve.bspline_nb_knots().unwrap_or(0);
            if nb_knots == 2 {
                if let Some(poles) = new_pcurve.poles2d() {
                    if poles.len() == 2 {
                        return line_through_poles(&poles);
                    }
                }
            }
        }
    }

    // `cxx:1141-1149`: the curve is not prolonged on either bound, so a
    // segment is appended at each end that stops short.
    let first_par = pcurve.first_parameter();
    let last_par = pcurve.last_parameter();
    let tr_curve: Arc<dyn Curve2d> =
        Arc::new(Geom2dTrimmedCurve::new(pcurve.clone(), first_par, last_par));
    let mut comp = CompCurveToBSplineCurve2d::from_curve(
        tr_curve.as_ref(),
        ParameterisationType::RationalC1,
    )?;
    let a_tol = CONFUSION;
    let a_delta = a2_offset.max(1.0);

    if first_par > an_ef - a2_offset {
        // `aPCurve->D1(FirstPar, aPBnd, aVBnd); aDBnd.SetXY(aVBnd.XY());`
        let (a_p_bnd, a_v_bnd) = pcurve.d1(first_par);
        let a_d_bnd = GpDir2d::from_vec2d(&a_v_bnd).ok()?;
        let a_p_beg = a_p_bnd
            .translated_vec(&GpVec2d::from_dir2d(&a_d_bnd).multiplied_scalar(-a_delta));
        let a_lin = Geom2dLine::from_pnt_dir(a_p_beg, a_d_bnd);
        let a_segment: Arc<dyn Curve2d> =
            Arc::new(Geom2dTrimmedCurve::new(Arc::new(a_lin), 0.0, a_delta));
        if !comp.add(a_segment.as_ref(), a_tol, true).ok()? {
            return None;
        }
    }

    if last_par < an_el + a2_offset {
        let (a_p_beg, a_v_bnd) = pcurve.d1(last_par);
        let a_d_bnd = GpDir2d::from_vec2d(&a_v_bnd).ok()?;
        let a_lin = Geom2dLine::from_pnt_dir(a_p_beg, a_d_bnd);
        let a_segment: Arc<dyn Curve2d> =
            Arc::new(Geom2dTrimmedCurve::new(Arc::new(a_lin), 0.0, a_delta));
        if !comp.add(a_segment.as_ref(), a_tol, true).ok()? {
            return None;
        }
    }

    // `NewPCurve = aCompCurve.BSplineCurve();` (`cxx:1185`).
    let bspl = comp.into_curve()?;
    Some(Arc::new(bspl))
}

/// Stable surface-slot key of a `CurveOnSurface` representation.
fn slot_key(surface: &Arc<dyn Surface>) -> usize {
    Arc::as_ptr(surface) as *const () as usize
}

/// `BRep_GCurve::SetRange(First, Last)` on one representation.
fn set_rep_range(edge: &TopoShape, key: usize, first: f64, last: f64) {
    edge.tshape
        .write()
        .unwrap()
        .edge_pcurves_mut()
        .ranges
        .insert(key, (first, last));
}

/// `BRep_GCurve::First()` / `Last()` of one representation.
fn rep_range(edge: &TopoShape, key: usize) -> Option<(f64, f64)> {
    let ts = edge.tshape.read().unwrap();
    ts.edge_pcurves().and_then(|g| g.ranges.get(&key).copied())
}

/// `CurveRep->PCurve(NewPCurve)` / `PCurve2(NewPCurve)`: replace the pcurve
/// list of one representation and leave its `[First, Last]` untouched.
fn set_rep_curves(edge: &TopoShape, key: usize, curves: Vec<Arc<dyn Curve2d>>) {
    edge.tshape
        .write()
        .unwrap()
        .edge_pcurves_mut()
        .curves
        .set(key, curves);
}

/// `BRepOffset_Inter2d::ExtentEdge(E, NE, theOffset)` (`cxx:1193-1676`).
///
/// `None` is the reference's `false` return (the caller propagates it out of
/// `ConnexIntByInt`).
pub fn extent_edge(e: &Edge, offset: f64) -> Option<Edge> {
    let reg = GeometryRegistry::global();
    let (an_ef, an_el) = BRepTool::edge_parameters(e);
    let a2_offset = 2.0 * offset.abs();

    // `NE = TopoDS::Edge(E.EmptyCopied())` (`cxx:1197`, `:1207`).
    let ne = e.0.empty_copied();
    if let Some(g) = reg.edge_geom(&e.0) {
        reg.set_edge(&ne, g);
    }

    // One entry per CurveOnSurface representation: the surface, its pcurves
    // (forward-then-reversed for a closed surface) and the COS range.
    let mut reps: Vec<(Arc<dyn Surface>, Vec<Arc<dyn Curve2d>>, (f64, f64))> = Vec::new();
    for (surf, pc, f, l) in reg.edge_pcurve_reps(&e.0) {
        let key = slot_key(&surf);
        match reps.iter_mut().find(|r| slot_key(&r.0) == key) {
            Some(slot) => slot.1.push(pc),
            None => reps.push((surf, vec![pc], (f, l))),
        }
    }
    reps.sort_by_key(|r| slot_key(&r.0));

    let mut nb_pcurves = 0usize;
    let mut first_par_on_pc = REAL_FIRST;
    let mut last_par_on_pc = REAL_LAST;
    let mut min_pc: Option<Arc<dyn Curve2d>> = None;
    let mut min_surf: Option<Arc<dyn Surface>> = None;

    for (surf, pcs, _range) in &reps {
        nb_pcurves += 1;
        let key = slot_key(surf);
        let mut curves = pcs.clone();
        let mut the_curve = curves[0].clone();
        let mut first_par = the_curve.first_parameter();
        let mut last_par = the_curve.last_parameter();

        // `cxx:1231-1252`.
        if is_bounded_curve2d(the_curve.as_ref())
            && (first_par > an_ef - a2_offset || last_par < an_el + a2_offset)
        {
            if let Some(new_pc) = extend_pcurve(&the_curve, an_ef, an_el, a2_offset) {
                curves[0] = new_pc.clone();
                the_curve = new_pc.clone();
                first_par = new_pc.first_parameter();
                last_par = new_pc.last_parameter();
                // `cxx:1242-1249`: the second pcurve of a closed surface.
                if curves.len() > 1 {
                    if let Some(new_pc2) = extend_pcurve(&curves[1], an_ef, an_el, a2_offset) {
                        curves[1] = new_pc2;
                    }
                }
            }
        } else if the_curve.is_periodic() {
            let mut delta = (the_curve.period() - (an_el - an_ef)) * 0.5;
            delta *= 0.95;
            first_par = an_ef - delta;
            last_par = an_el + delta;
        } else if is_closed2d(&the_curve) {
            last_par -= 0.05 * (last_par - first_par);
        }
        let rep_dirty = !Arc::ptr_eq(&curves[0], &pcs[0])
            || (curves.len() > 1 && !Arc::ptr_eq(&curves[1], &pcs[1]));
        if rep_dirty {
            set_rep_curves(&ne, key, curves.clone());
        }
        // `theCurve = CurveRep->PCurve();` (`cxx:1264`): the possibly
        // prolonged pcurve.
        let the_curve = curves[0].clone();

        // `cxx:1265-1279`: the four boundary lines of the surface's
        // parametric rectangle.
        let (u_min, u_max) = surf.u_range();
        let (v_min, v_max) = surf.v_range();
        let mut bound_lines: Vec<Arc<dyn Curve2d>> = Vec::new();
        if !Precision::is_infinite(v_min) {
            bound_lines.push(Arc::new(Geom2dLine::from_pnt_dir(
                GpPnt2d::new(0.0, v_min),
                GpDir2d::new(1.0, 0.0).ok()?,
            )));
        }
        if !Precision::is_infinite(u_min) {
            bound_lines.push(Arc::new(Geom2dLine::from_pnt_dir(
                GpPnt2d::new(u_min, 0.0),
                GpDir2d::new(0.0, 1.0).ok()?,
            )));
        }
        if !Precision::is_infinite(v_max) {
            bound_lines.push(Arc::new(Geom2dLine::from_pnt_dir(
                GpPnt2d::new(0.0, v_max),
                GpDir2d::new(1.0, 0.0).ok()?,
            )));
        }
        if !Precision::is_infinite(u_max) {
            bound_lines.push(Arc::new(Geom2dLine::from_pnt_dir(
                GpPnt2d::new(u_max, 0.0),
                GpDir2d::new(0.0, 1.0).ok()?,
            )));
        }

        // `cxx:1295-1337`: every intersection of the pcurve with a boundary
        // line that lies inside the rectangle.
        let mut params: Vec<f64> = Vec::new();
        let d_curve = Geom2dIntGInter::domain_over_range(
            the_curve.as_ref(),
            the_curve.first_parameter(),
            the_curve.last_parameter(),
            PCONFUSION,
        );
        for bound in &bound_lines {
            let mut int_cc = Geom2dIntGInter::new();
            let d_line = Geom2dIntGInter::domain_over_range(
                bound.as_ref(),
                bound.first_parameter(),
                bound.last_parameter(),
                PCONFUSION,
            );
            int_cc.perform(
                the_curve.as_ref(),
                &d_curve,
                bound.as_ref(),
                &d_line,
                PCONFUSION,
                PCONFUSION,
            );
            if !int_cc.is_done() {
                continue;
            }
            let in_box = |p: &GpPnt2d| {
                p.x() >= u_min && p.x() <= u_max && p.y() >= v_min && p.y() <= v_max
            };
            for j in 1..=int_cc.nb_points() {
                let ip = int_cc.point(j);
                if in_box(ip.value()) {
                    params.push(ip.param_on_first());
                }
            }
            for j in 1..=int_cc.nb_segments() {
                let is = int_cc.segment(j);
                if is.has_first_point() {
                    let ip = is.first_point();
                    if in_box(ip.value()) {
                        params.push(ip.param_on_first());
                    }
                }
                if is.has_last_point() {
                    let ip = is.last_point();
                    if in_box(ip.value()) {
                        params.push(ip.param_on_first());
                    }
                }
            }
        }

        // `cxx:1338-1379`.
        if !params.is_empty() {
            if params.len() == 1 {
                let pnt_first = the_curve.d0(first_par);
                if pnt_first.x() >= u_min
                    && pnt_first.x() <= u_max
                    && pnt_first.y() >= v_min
                    && pnt_first.y() <= v_max
                {
                    if last_par > params[0] {
                        last_par = params[0];
                    }
                } else if first_par < params[0] {
                    first_par = params[0];
                }
            } else {
                let mut fpar = REAL_LAST;
                let mut lpar = REAL_FIRST;
                for p in &params {
                    if *p < fpar {
                        fpar = *p;
                    }
                    if *p > lpar {
                        lpar = *p;
                    }
                }
                if first_par < fpar {
                    first_par = fpar;
                }
                if last_par > lpar {
                    last_par = lpar;
                }
            }
        }

        // `cxx:1381`.
        set_rep_range(&ne, key, first_par, last_par);

        // `cxx:1386-1401`.
        if first_par > first_par_on_pc {
            first_par_on_pc = first_par;
            min_pc = Some(the_curve.clone());
            min_surf = Some(surf.clone());
        }
        if last_par < last_par_on_pc {
            last_par_on_pc = last_par;
            min_pc = Some(the_curve);
            min_surf = Some(surf.clone());
        }
    }

    // `double f, l; C3d = BRep_Tool::Curve(NE, f, l);` (`cxx:1404-1405`).
    let (c3d_f, c3d_l) = BRepTool::edge_parameters(&Edge(ne.clone()));
    let c3d = reg.edge_curve(&ne);

    if nb_pcurves > 0 {
        // `MinLoc = E.Location() * MinLoc;` (`cxx:1406`).
        let min_loc = e.0.location();
        let (Some(min_pc), Some(min_surf)) = (min_pc.clone(), min_surf.clone()) else {
            // `NbPCurves` non-zero always sets `MinPC` (`cxx:1386-1401`).
            return Some(Edge(ne));
        };

        if let Some(c3d) = c3d {
            let mut f = c3d_f;
            let mut l = c3d_l;
            // `cxx:1411-1466`.
            if is_closed2d(&min_pc) {
                f = first_par_on_pc;
                l = last_par_on_pc;
            } else if c3d.is_periodic() {
                let mut delta = (c3d.period() - (l - f)) * 0.5;
                delta *= 0.95;
                f -= delta;
                l += delta;
            } else if is_closed3d(c3d.as_ref()) {
                l -= 0.05 * (l - f);
            } else {
                f = first_par_on_pc;
                l = last_par_on_pc;
                let min_loc_trsf = min_loc.transformation();
                if !Precision::is_infinite(first_par_on_pc) {
                    let p2d = min_pc.d0(first_par_on_pc);
                    let p1 = min_surf.value(p2d.x(), p2d.y()).transformed(&min_loc_trsf);
                    if let Some(proj) = project_point_on_curve(c3d.as_ref(), &p1, 0.0) {
                        f = proj.parameter;
                    }
                }
                if !Precision::is_infinite(last_par_on_pc) {
                    let p2d = min_pc.d0(last_par_on_pc);
                    let p2 = min_surf.value(p2d.x(), p2d.y()).transformed(&min_loc_trsf);
                    if let Some(proj) = project_point_on_curve(c3d.as_ref(), &p2, 0.0) {
                        l = proj.parameter;
                    }
                }
            }
            // `BB.Range(NE, f, l)` (`cxx:1464`).
            reg.set_ranges_all(&ne, f, l);
            if !Precision::is_infinite(f) && !Precision::is_infinite(l) {
                same_parameter_edge_inplace(&Edge(ne.clone()), CONFUSION);
            }
        } else if !reg.is_degenerated_edge(&e.0) {
            // No 3D curve (`cxx:1471-1603`).
            let min_loc_trsf = min_loc.transformation();
            let min_surf: Arc<dyn Surface> = Arc::from(min_surf.transformed(&min_loc_trsf));
            let mut max_deviation = 0.0;
            let mut new_c3d: Option<Arc<dyn Curve>> = None;

            if Precision::is_infinite(first_par_on_pc) || Precision::is_infinite(last_par_on_pc)
            {
                // `cxx:1477-1501`: an infinite pcurve on a plane, or a
                // y-parallel one on a cylinder / cone, is a straight 3D line.
                if min_pc.is_geom2d_line() {
                    let mut is_line = min_surf.gp_pln().is_some();
                    if !is_line
                        && (min_surf.gp_cylinder().is_some() || min_surf.gp_cone().is_some())
                    {
                        if let Some(the_line) = min_pc.gp_lin2d() {
                            let dy2d = GpDir2d::new(0.0, 1.0).ok()?;
                            if the_line.direction().is_parallel(&dy2d, ANGULAR) {
                                is_line = true;
                            }
                        }
                    }
                    if is_line {
                        let p2d_1 = min_pc.d0(0.0);
                        let p2d_2 = min_pc.d0(1.0);
                        let p1 = min_surf.value(p2d_1.x(), p2d_1.y());
                        let p2 = min_surf.value(p2d_2.x(), p2d_2.y());
                        let vec = occt_core::gp::GpVec::from_pnts(&p1, &p2);
                        let dir = GpDir::new(vec.x(), vec.y(), vec.z()).ok()?;
                        new_c3d = Some(Arc::new(GeomLine::from_pnt_dir(p1, dir)));
                    }
                }
            } else {
                // `cxx:1504-1524`.
                let cos: Arc<dyn Curve> = Arc::new(CurveOnSurface::new(
                    min_pc.clone(),
                    min_surf.clone(),
                    first_par_on_pc,
                    last_par_on_pc,
                ));
                let max_segment = evaluate_max_segment(&min_surf, &min_pc);
                let built = geom_lib_build_curve3d(
                    CONFUSION,
                    cos.as_ref(),
                    &min_surf,
                    &min_pc,
                    first_par_on_pc,
                    last_par_on_pc,
                    Shape::C1,
                    14,
                    max_segment,
                );
                new_c3d = built.curve;
                max_deviation = built.max_deviation;
            }

            // `BB.UpdateEdge(NE, C3d, max_deviation);` (`cxx:1525`). A null
            // handle leaves the existing (null) `BRep_Curve3D` representation
            // in place, which is the state `NE` was copied in.
            if let Some(c) = &new_c3d {
                reg.set_edge_curve3d(&ne, c.clone(), max_deviation);
            }

            // `cxx:1527-1595`: re-project the other pcurves onto their
            // surfaces through the rebuilt 3D curve.
            let mut projection_success = true;
            if nb_pcurves > 1 {
                // The pcurves `NE` carries now (the main loop may have
                // prolonged them); `CurveRep->PCurve()` is re-read per
                // representation in the reference too (`cxx:1538`).
                let mut reps_now: Vec<(Arc<dyn Surface>, Vec<Arc<dyn Curve2d>>)> = Vec::new();
                for (surf, pc, _f, _l) in reg.edge_pcurve_reps(&ne) {
                    let k = slot_key(&surf);
                    match reps_now.iter_mut().find(|r| slot_key(&r.0) == k) {
                        Some(slot) => slot.1.push(pc),
                        None => reps_now.push((surf, vec![pc])),
                    }
                }
                reps_now.sort_by_key(|r| slot_key(&r.0));

                for (surf, pcs) in &reps_now {
                    let key = slot_key(surf);
                    let the_curve = pcs[0].clone();
                    // `theCurve == MinPC && theSurf == MinSurf && theLoc ==
                    // MinLoc` (`cxx:1541`): handle identity.
                    if Arc::ptr_eq(&the_curve, &min_pc) && Arc::ptr_eq(surf, &min_surf) {
                        continue;
                    }
                    let (first_par, last_par) = rep_range(&ne, key).unwrap_or((c3d_f, c3d_l));
                    if (first_par - first_par_on_pc).abs() > PCONFUSION
                        || (last_par - last_par_on_pc).abs() > PCONFUSION
                    {
                        let the_loc = e.0.location();
                        let the_surf: Arc<dyn Surface> =
                            Arc::from(surf.transformed(&the_loc.transformation()));

                        if the_curve.is_geom2d_line() && is_bounded_surface(the_surf.as_ref()) {
                            if let Some(the_line) = the_curve.gp_lin2d() {
                                let dx2d = GpDir2d::new(1.0, 0.0).ok()?;
                                let dy2d = GpDir2d::new(0.0, 1.0).ok()?;
                                if the_line.direction().is_parallel(&dx2d, ANGULAR)
                                    || the_line.direction().is_parallel(&dy2d, ANGULAR)
                                {
                                    let (u1, u2) = the_surf.u_range();
                                    let (v1, v2) = the_surf.v_range();
                                    let origin = the_line.location();
                                    if (origin.x() - u1).abs() <= CONFUSION
                                        || (origin.x() - u2).abs() <= CONFUSION
                                        || (origin.y() - v1).abs() <= CONFUSION
                                        || (origin.y() - v2).abs() <= CONFUSION
                                    {
                                        same_parameter_edge_inplace(&Edge(ne.clone()), CONFUSION);
                                        break;
                                    }
                                }
                            }
                        }
                        if let Some(c) = &new_c3d {
                            if first_par_on_pc < last_par_on_pc {
                                let mut tol = PCONFUSION;
                                match curve_2d_range(
                                    c,
                                    first_par_on_pc,
                                    last_par_on_pc,
                                    &the_surf,
                                    &mut tol,
                                ) {
                                    Some(proj_pcurve) => {
                                        let mut curves = pcs.clone();
                                        curves[0] = proj_pcurve;
                                        set_rep_curves(&ne, key, curves);
                                    }
                                    None => projection_success = false,
                                }
                            } else {
                                return None;
                            }
                        } else {
                            return None;
                        }
                    }
                }
            }

            // `cxx:1594-1602`.
            if projection_success {
                reg.set_ranges_all(&ne, first_par_on_pc, last_par_on_pc);
            } else {
                reg.set_edge_range(&ne, first_par_on_pc, last_par_on_pc);
                same_parameter_edge_inplace(&Edge(ne.clone()), CONFUSION);
            }
        }
    } else {
        // No pcurves (`cxx:1605-1672`).
        let Some(c3d) = c3d else {
            // Port guard: `cxx:1607` dereferences the null handle.
            return Some(Edge(ne));
        };
        let mut first_par = c3d.first_parameter();
        let mut last_par = c3d.last_parameter();

        if is_bounded_curve3d(c3d.as_ref())
            && (first_par > an_ef - a2_offset || last_par < an_el + a2_offset)
        {
            let tr_curve: Arc<dyn Curve> =
                Arc::new(occt_geom::trimmed::GeomTrimmedCurve::new(
                    c3d.clone(),
                    first_par,
                    last_par,
                ));
            let mut comp = occt_geom::convert_bspl::CompCurveToBSplineCurve::from_curve(
                tr_curve.as_ref(),
                ParameterisationType::RationalC1,
            )
            .ok()?;
            let a_tol = CONFUSION;
            let a_delta = a2_offset.max(1.0);

            if first_par > an_ef - a2_offset {
                let (a_p_bnd, a_v_bnd) = c3d.d1(first_par);
                let a_d_bnd = GpDir::new(a_v_bnd.x(), a_v_bnd.y(), a_v_bnd.z()).ok()?;
                let a_p_beg = a_p_bnd.translated_vec(&occt_core::gp::GpVec::new(
                    -a_delta * a_d_bnd.x(),
                    -a_delta * a_d_bnd.y(),
                    -a_delta * a_d_bnd.z(),
                ));
                let a_lin = GeomLine::from_pnt_dir(a_p_beg, a_d_bnd);
                let a_segment: Arc<dyn Curve> = Arc::new(occt_geom::trimmed::GeomTrimmedCurve::new(
                    Arc::new(a_lin),
                    0.0,
                    a_delta,
                ));
                if !comp
                    .add(a_segment.as_ref(), a_tol, true, true, 0)
                    .ok()?
                {
                    return Some(Edge(ne));
                }
            }
            if last_par < an_el + a2_offset {
                let (a_p_beg, a_v_bnd) = c3d.d1(last_par);
                let a_d_bnd = GpDir::new(a_v_bnd.x(), a_v_bnd.y(), a_v_bnd.z()).ok()?;
                let a_lin = GeomLine::from_pnt_dir(a_p_beg, a_d_bnd);
                let a_segment: Arc<dyn Curve> = Arc::new(occt_geom::trimmed::GeomTrimmedCurve::new(
                    Arc::new(a_lin),
                    0.0,
                    a_delta,
                ));
                if !comp
                    .add(a_segment.as_ref(), a_tol, true, true, 0)
                    .ok()?
                {
                    return Some(Edge(ne));
                }
            }
            let bspl = comp.into_curve()?;
            let new_c3d: Arc<dyn Curve> = Arc::new(bspl);
            first_par = new_c3d.first_parameter();
            last_par = new_c3d.last_parameter();
            reg.set_edge_curve3d(&ne, new_c3d, CONFUSION);
        } else if c3d.is_periodic() {
            let mut delta = (c3d.period() - (an_el - an_ef)) * 0.5;
            delta *= 0.95;
            first_par = an_ef - delta;
            last_par = an_el + delta;
        } else if is_closed3d(c3d.as_ref()) {
            last_par -= 0.05 * (last_par - first_par);
        }

        // `BB.Range(NE, FirstPar, LastPar);` (`cxx:1670`).
        reg.set_ranges_all(&ne, first_par, last_par);
    }

    Some(Edge(ne))
}

/// `C3d->IsClosed()` (`Geom_Curve::IsClosed`): endpoint coincidence within
/// `Precision::Confusion()`.
fn is_closed3d(c: &dyn Curve) -> bool {
    let (first, last) = (c.first_parameter(), c.last_parameter());
    if Precision::is_infinite(first) || Precision::is_infinite(last) {
        return false;
    }
    c.d0(first).distance(&c.d0(last)) <= CONFUSION
}

/// `Geom2d_Curve::IsClosed()`: `Value(First).Distance(Value(Last)) <=
/// Precision::Confusion()`.
fn is_closed2d(c: &Arc<dyn Curve2d>) -> bool {
    let (first, last) = (c.first_parameter(), c.last_parameter());
    if Precision::is_infinite(first) || Precision::is_infinite(last) {
        return false;
    }
    c.d0(first).distance(&c.d0(last)) <= CONFUSION
}
