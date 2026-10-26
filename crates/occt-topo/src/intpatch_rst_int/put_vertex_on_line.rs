//! `IntPatch_RstInt::PutVertexOnLine` (`IntPatch_RstInt.cxx:433-1355`).
//!
//! The driver: for every restriction arc of `Domain` (the 2d domain of `Surf`),
//! sample it with an `IntPatch_PolyArc`, intersect it with the 2d polygon of the
//! line (`Intf_InterferencePolygon2d`), refine each section point on the
//! arc/surface pair with `IntPatch_CurvIntSurf` (`IntPatch_CSFunction`), project
//! the refined point back onto the line, and write it there as a vertex.
//!
//! Translation notes:
//! * `PLin` holds a live `occ::handle<IntPatch_PolyLine>` on `L` and reads the
//!   line's *points* through it. The vertices are written while the polygon is
//!   alive, so the port snapshots the line into `snapshot` (`GeomIntLine` is
//!   `Clone`) and builds the polygon on the snapshot. Only the vertex sequence
//!   changes here, which the polygon does not read.
//! * `Nbvtx` is updated in place (`cxx:980-1226` appends to the sequence) and is
//!   returned by [`super::vertex_write`]'s arms.
//! * `p2dFirst` / `p2dLast` (`cxx:573-587`) are written by the source and never
//!   read; the `D0` calls are kept, the values dropped.
//! * `nbTreated` (`cxx:725`, `:979`, `:1229`) can reach 1 at most, so the
//!   `nbTreated == 2` arm is unreachable in OCCT as well; it is kept verbatim.
//! * `debug_polygon2d` / `NumeroEdge` (`cxx:548`, `:703-709`) only feed the
//!   debug dump and are not translated.

use occt_core::bnd::BndBox2d;
use occt_core::elib::clib::in_period;
use occt_core::gp::{GpDir, GpPnt, GpPnt2d, GpVec};
use occt_core::intf::{IntfInterferencePolygon2d, IntfPolygon2d};
use occt_core::precision::{Precision, CONFUSION};
use occt_geom::Surface;
use occt_geom2d::curve::Curve2d;

use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_surface::{classify_surface, SurfaceKind};
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::geom_int::{GeomIntLine, IntPatchIType};
use crate::int_tools_wline::{u_resolution, v_resolution, PatchPoint, TransType};
use crate::intpatch_curv_int_surf::IntPatchCurvIntSurf;
use crate::intpatch_cs_function::IntPatchCSFunction;
use crate::intpatch_h_inter_tool as hinter;
use crate::intpatch_poly_arc::IntPatchPolyArc;
use crate::intpatch_poly_line::IntPatchPolyLine;
use crate::intpatch_vparams;
use crate::intsurf_transition::{make_transition, IntSurfTransition};

use super::geom_helpers::{
    coincide_on_arc, find_parameter, get_line_point2d, recadre, tol_3d_arc, tol_3d_vertex,
    verify_tgline,
};
use super::vertex_write::{append_new_vertex, update_existing_vertex, VertexUpdate};

/// `myInfinite` (`IntPatch_RstInt.cxx:31`): `RealLast()`.
const MY_INFINITE: f64 = f64::MAX;

/// `PutVertexOnLine(L, Surf, Domain, OtherSurf, OnFirst, Tol)`
/// (`cxx:433-1355`).
///
/// `Surf` is the surface whose restriction `Domain` describes, `OtherSurf` the
/// surface the line is intersected with; the restriction arcs of `Domain` are
/// intersected with `OtherSurf`.
pub fn put_vertex_on_line(
    line: &mut GeomIntLine,
    surf: &dyn Surface,
    domain: &mut TopolTool,
    other_surf: &dyn Surface,
    on_first: bool,
    tol: f64,
) -> Result<(), String> {
    let two_pi = std::f64::consts::PI + std::f64::consts::PI;

    let mut locpt: Vec<GpPnt> = Vec::new();
    let mut locpt2: Vec<GpPnt2d> = Vec::new();

    let typ_l = line.arc_type();
    let mut nbvtx = line.nb_vertex();

    // `tolPLin` (`cxx:459-461`).
    let tol_plin = u_resolution(surf, CONFUSION)
        .max(v_resolution(surf, CONFUSION))
        .min(CONFUSION);

    // `IntPatch_PolyLine PLin(tolPLin); PLin.SetWLine/SetRLine(OnFirst, L)`
    // (`cxx:462`, `:479-492`).
    let snapshot = line.clone();
    let mut plin = IntPatchPolyLine::with_deflection(tol_plin);
    match typ_l {
        IntPatchIType::Walking => plin.set_wline(
            on_first,
            snapshot
                .as_wline()
                .ok_or_else(|| "Standard_DomainError".to_string())?,
        ),
        IntPatchIType::Restriction => plin.set_rline(
            on_first,
            snapshot
                .as_rline()
                .ok_or_else(|| "Standard_DomainError".to_string())?,
        ),
        _ => return Err("Standard_DomainError".to_string()),
    }
    if !domain.has3d() {
        // `cxx:493-497`: don't use the computed deflection in the mode of pure
        // geometric intersection.
        plin.reset_error();
    }

    // `cxx:499-524`.
    let surface_is_u_closed = surf.is_u_closed();
    let surface_is_v_closed = surf.is_v_closed();
    let o_surface_is_u_closed = other_surf.is_u_closed();
    let o_surface_is_v_closed = other_surf.is_v_closed();
    let possibly_closed = surface_is_u_closed
        || surface_is_v_closed
        || o_surface_is_u_closed
        || o_surface_is_v_closed;
    let mut tol_u_closed = 0.0;
    let mut tol_v_closed = 0.0;
    let mut tol_ou_closed = 0.0;
    let mut tol_ov_closed = 0.0;
    if possibly_closed {
        if surface_is_u_closed {
            let (u0, u1) = surf.u_range();
            tol_u_closed = (u1 - u0) * 0.01;
        }
        if surface_is_v_closed {
            let (v0, v1) = surf.v_range();
            tol_v_closed = (v1 - v0) * 0.01;
        }
        if o_surface_is_u_closed {
            let (u0, u1) = other_surf.u_range();
            tol_ou_closed = (u1 - u0) * 0.01;
        }
        if o_surface_is_v_closed {
            let (v0, v1) = other_surf.v_range();
            tol_ov_closed = (v1 - v0) * 0.01;
        }
    }

    // `cxx:531-546`.
    let surf1: &dyn Surface = if on_first { surf } else { other_surf };
    let surf2: &dyn Surface = if on_first { other_surf } else { surf };
    let type_s1 = classify_surface(surf1);
    let type_s2 = classify_surface(surf2);
    let surfacetype = if on_first { type_s1 } else { type_s2 };
    let surface_is_periodic = matches!(
        surfacetype,
        SurfaceKind::Cylinder | SurfaceKind::Cone | SurfaceKind::Torus | SurfaceKind::Sphere
    );
    let surface_is_bi_periodic = surfacetype == SurfaceKind::Torus;

    // `Domain->Init()` / `More()` / `Value()` / `Next()` (`cxx:549-553`): walk
    // the restriction arcs of the domain. Each arc is materialised once so its
    // address is a stable `IntPatch_Point::ArcOnS1/2` identity, the convention
    // `intpatch_line_ctor::wline_to_rline::domain_arcs` established.
    let arcs: Vec<std::sync::Arc<Curve2dAdaptor>> = {
        let mut out = Vec::new();
        domain.init();
        while domain.more() {
            if let Some(c) = domain.value() {
                out.push(std::sync::Arc::new(c.clone()));
            }
            domain.next();
        }
        out
    };
    let mut arc_iter = arcs.iter();

    domain.init();
    while domain.more() {
        let arc_handle: std::sync::Arc<Curve2dAdaptor> = arc_iter
            .next()
            .ok_or_else(|| "Standard_DomainError".to_string())?
            .clone();
        // `cxx:553`: `arc = Domain->Value();`. The domain arc is materialised
        // once per call (`domain_arcs` convention of
        // `intpatch_line_ctor::wline_to_rline`) so its address is a stable
        // `IntPatch_Point::ArcOnS1/2` identity.
        let arc = arc_handle.as_ref();
        let arc_id = crate::int_tools_wline::arc_id(&arc_handle);
        let arc_curve: &dyn Curve2d = arc
            .curve()
            .ok_or_else(|| "Standard_DomainError".to_string())?
            .as_ref();

        // `cxx:556`: use the tolerance of this edge if possible.
        let edge_tol = tol_3d_arc(arc, domain, tol);

        // `IntPatch_HInterTool::Bounds` (`cxx:558-566`).
        let mut pfirst = arc.first_parameter();
        let mut plast = arc.last_parameter();
        if Precision::is_negative_infinite(pfirst) {
            pfirst = -MY_INFINITE;
        }
        if Precision::is_positive_infinite(plast) {
            plast = MY_INFINITE;
        }

        // `cxx:573-587`: the source stores the arc end points and never reads
        // them; the calls are kept, the values dropped.
        domain.initialize_curve(arc);
        domain.init_vertex_iterator();
        while domain.more_vertex() {
            if let Some(vtx) = domain.vertex() {
                let prm = hinter::parameter(&vtx, arc);
                if (prm - pfirst).abs() < occt_core::precision::PCONFUSION {
                    let _ = arc.d0(pfirst);
                } else if (prm - plast).abs() < occt_core::precision::PCONFUSION {
                    let _ = arc.d0(plast);
                }
            }
            domain.next_vertex();
        }

        let bplin: BndBox2d = plin.bounding().clone();
        let mut offset_v = 0.0;
        let mut offset_u = 0.0;

        // `switch (arc->GetType())` (`cxx:593-652`).
        let mut nb_echant: usize;
        if arc_curve.is_line() {
            nb_echant = 10;
            let (a_xmin, a_ymin, ..) = bplin
                .get()
                .ok_or_else(|| "IntPatch_PolyLine: void box".to_string())?;
            let a_lin = arc_curve
                .gp_lin2d()
                .ok_or_else(|| "Standard_DomainError".to_string())?;
            let a_loc = a_lin.location();
            let a_dir = *a_lin.direction();
            // Here, we consider rectangular axis-aligned domain only.
            let is_along_u = a_dir.x().abs() > a_dir.y().abs();
            if surface_is_periodic && !is_along_u {
                // Shift along U-direction.
                let a_new_location = in_period(a_loc.x(), a_xmin, a_xmin + two_pi);
                offset_u = a_new_location - a_loc.x();
            } else if surface_is_bi_periodic && is_along_u {
                // Shift along V-direction.
                let a_new_location = in_period(a_loc.y(), a_ymin, a_ymin + two_pi);
                offset_v = a_new_location - a_loc.y();
            }
        } else if arc_curve.is_bezier2d() {
            nb_echant = 3 + arc_curve.bezier_nb_poles().unwrap_or(0);
            if nb_echant < 10 {
                nb_echant = 10;
            } else if nb_echant > 50 {
                nb_echant = 50;
            }
        } else if arc_curve.bspline_nb_knots().is_some() {
            // `szv:const double nbs = ...` (`cxx:634-647`).
            let nbs = (arc_curve.bspline_nb_knots().unwrap_or(0) as f64
                * arc_curve.bspline_degree().unwrap_or(0) as f64)
                * (plast - pfirst)
                / (arc_curve.last_parameter() - arc_curve.first_parameter());
            nb_echant = if nbs < 2.0 { 2 } else { nbs as usize };
            if nb_echant < 10 {
                nb_echant = 10;
            } else if nb_echant > 50 {
                nb_echant = 50;
            }
        } else {
            nb_echant = 25;
        }

        // `cxx:654-671`.
        let mut bplin_ext = bplin.clone();
        if surface_is_periodic {
            let (xmin, ymin, xmax, ymax) = bplin_ext
                .get()
                .ok_or_else(|| "IntPatch_PolyLine: void box".to_string())?;
            let g = bplin_ext.gap();
            bplin_ext.set_void();
            bplin_ext.update(xmin - two_pi, ymin, xmax + two_pi, ymax);
            bplin_ext.set_gap(g);
        }
        if surface_is_bi_periodic {
            let (xmin, ymin, xmax, ymax) = bplin_ext
                .get()
                .ok_or_else(|| "IntPatch_PolyLine: void box".to_string())?;
            let g = bplin_ext.gap();
            bplin_ext.set_void();
            bplin_ext.update(xmin, ymin - two_pi, xmax, ymax + two_pi);
            bplin_ext.set_gap(g);
        }

        let mut brise = IntPatchPolyArc::new(arc_curve, nb_echant, pfirst, plast, &bplin_ext)?;

        let mut indice_offset_bi_periodic = 0;
        let mut indice_offset_periodic = 0;
        let a_ref_ou = offset_u;
        let a_ref_ov = offset_v;

        loop {
            if indice_offset_bi_periodic == 1 {
                offset_v = a_ref_ov - two_pi;
            } else if indice_offset_bi_periodic == 2 {
                offset_v = a_ref_ov + two_pi;
            }
            loop {
                if indice_offset_periodic == 1 {
                    offset_u = a_ref_ou - two_pi;
                } else if indice_offset_periodic == 2 {
                    offset_u = a_ref_ou + two_pi;
                }

                brise.set_offset(offset_u, offset_v);

                let mut commun = IntfInterferencePolygon2d::new();
                commun.perform_two(&plin, &brise);
                locpt.clear();
                locpt2.clear();

                // Only the section points are wanted; tangent zones are skipped
                // (`cxx:715-719`).
                let a_nb_section_pts = commun.base.nb_section_points();
                for i in 1..=a_nb_section_pts {
                    let a_w1 = commun.base.pnt_value(i).param_on_first();
                    let a_w2 = commun.base.pnt_value(i).param_on_second();

                    let mut nb_treated = 0;

                    let (u, v) = get_line_point2d(line, a_w1 + 1.0, !on_first)?;

                    let mut par = a_w2.trunc();
                    let mut irang = par as i32 + 1;
                    if irang as usize == brise.nb_points() {
                        irang -= 1;
                        par = 1.0;
                    } else {
                        par = (a_w2 - par).abs();
                    }
                    let w = (1.0 - par) * brise.parameter(irang as usize)
                        + par * brise.parameter(irang as usize + 1);

                    let thefunc = IntPatchCSFunction::new(other_surf, arc_curve, surf);
                    // MSV: extend the UV bounds to not miss a solution near the
                    // boundary (`cxx:750`).
                    let marg_coef = 0.004;
                    let int_cs = IntPatchCurvIntSurf::new(u, v, w, &thefunc, edge_tol, marg_coef);

                    let mut refined = false;
                    let mut u2 = 0.0;
                    let mut v2 = 0.0;
                    let mut paramarc = 0.0;
                    let mut ptsommet = GpPnt::zero();

                    if int_cs.is_done() && !int_cs.is_empty()? {
                        ptsommet = *int_cs.point()?;
                        let (ru2, rv2) = int_cs.parameter_on_surface()?;
                        u2 = ru2;
                        v2 = rv2;
                        let an_old_pnt = other_surf.d0(u, v);
                        let a_new_pnt = other_surf.d0(u2, v2);
                        let a_tol_conf = CONFUSION.max(edge_tol);
                        if an_old_pnt.square_distance(&a_new_pnt) < a_tol_conf * a_tol_conf {
                            u2 = u;
                            v2 = v;
                        }
                        paramarc = int_cs.parameter_on_curve()?;
                        refined = true;
                    }

                    if !refined {
                        continue;
                    }

                    let mut duplicate = false;
                    for j in 0..locpt.len() {
                        if ptsommet.distance(&locpt[j]) <= edge_tol {
                            if possibly_closed {
                                locpt2[j] = GpPnt2d::new(u, v);
                                if (o_surface_is_u_closed && (u - u2).abs() > tol_ou_closed)
                                    || (o_surface_is_v_closed && (v - v2).abs() > tol_ov_closed)
                                {
                                    continue;
                                }
                            }
                            duplicate = true;
                            break;
                        }
                    }
                    if duplicate {
                        continue;
                    }

                    let param_approx_on_line = a_w1 as i32 + 1;

                    let (p2d, d2d) = arc
                        .d1(paramarc)
                        .ok_or_else(|| "Standard_DomainError".to_string())?;
                    let mut u1 = p2d.x();
                    let mut v1 = p2d.y();
                    if typ_l == IntPatchIType::Walking && surface_is_periodic {
                        if on_first {
                            let wlin = line
                                .as_wline()
                                .ok_or_else(|| "Standard_DomainError".to_string())?;
                            recadre(
                                type_s1,
                                type_s2,
                                wlin,
                                param_approx_on_line,
                                &mut u1,
                                &mut v1,
                                &mut u2,
                                &mut v2,
                            );
                        } else {
                            let wlin = line
                                .as_wline()
                                .ok_or_else(|| "Standard_DomainError".to_string())?;
                            recadre(
                                type_s1,
                                type_s2,
                                wlin,
                                param_approx_on_line,
                                &mut u2,
                                &mut v2,
                                &mut u1,
                                &mut v1,
                            );
                        }
                    }
                    locpt.push(ptsommet);
                    locpt2.push(GpPnt2d::new(u2, v2));

                    let fp = find_parameter(
                        line,
                        other_surf,
                        edge_tol,
                        &ptsommet,
                        &GpPnt2d::new(u2, v2),
                        param_approx_on_line,
                        on_first,
                    )?;
                    let mut found = fp.is_some();
                    let (paramline, mut tgline) = fp.unwrap_or((0.0, GpVec::zero()));

                    if typ_l == IntPatchIType::Walking && found && possibly_closed {
                        // check in 2d (`cxx:824-845`).
                        if surface_is_u_closed || surface_is_v_closed {
                            let (lu, lv) = get_line_point2d(line, paramline, on_first)?;
                            if (surface_is_u_closed && (lu - u1).abs() > tol_u_closed)
                                || (surface_is_v_closed && (lv - v1).abs() > tol_v_closed)
                            {
                                found = false;
                            }
                        }
                        if found && (o_surface_is_u_closed || o_surface_is_v_closed) {
                            let (lu, lv) = get_line_point2d(line, paramline, !on_first)?;
                            if (o_surface_is_u_closed && (lu - u2).abs() > tol_ou_closed)
                                || (o_surface_is_v_closed && (lv - v2).abs() > tol_ov_closed)
                            {
                                found = false;
                            }
                        }
                    }
                    if !found {
                        continue;
                    }

                    let (mut vtx_on_arc, mut vtxarc) =
                        coincide_on_arc(&ptsommet, arc, surf, edge_tol, domain);
                    let mut vtx_tol = if vtx_on_arc {
                        let t = vtxarc
                            .as_ref()
                            .map(|v| tol_3d_vertex(v, domain, 0.0))
                            .unwrap_or(0.0);
                        t.max(edge_tol)
                    } else {
                        edge_tol
                    };

                    // It is necessary to test that the point does not already
                    // exist: it can be already a point on the arc, but on a
                    // different arc (`cxx:866-938`).
                    let mut ivtx = 0;
                    let mut dmin = f64::MAX;
                    let mut ptline = PatchPoint::new(GpPnt::zero(), 0.0, 0.0, 0.0, 0.0, 0.0);
                    for j in 1..=nbvtx {
                        let rptline = line
                            .vertex(j)
                            .cloned()
                            .ok_or_else(|| "Standard_DomainError".to_string())?;
                        let a_point_on_rst_still_exist = if on_first {
                            rptline.is_on_dom_s1() && rptline.arc_id(true) == arc_id
                        } else {
                            rptline.is_on_dom_s2() && rptline.arc_id(false) == arc_id
                        };
                        if a_point_on_rst_still_exist {
                            continue;
                        }
                        if possibly_closed {
                            if surface_is_u_closed || surface_is_v_closed {
                                let (ru, rv) = if on_first {
                                    rptline.parameters_on_s1()
                                } else {
                                    rptline.parameters_on_s2()
                                };
                                if (surface_is_u_closed && (ru - u1).abs() > tol_u_closed)
                                    || (surface_is_v_closed && (rv - v1).abs() > tol_v_closed)
                                {
                                    continue;
                                }
                            }
                            if o_surface_is_u_closed || o_surface_is_v_closed {
                                let (ru, rv) = if on_first {
                                    rptline.parameters_on_s2()
                                } else {
                                    rptline.parameters_on_s1()
                                };
                                if (o_surface_is_u_closed && (ru - u2).abs() > tol_ou_closed)
                                    || (o_surface_is_v_closed && (rv - v2).abs() > tol_ov_closed)
                                {
                                    continue;
                                }
                            }
                        }
                        let dist = ptsommet.distance(&rptline.p);
                        let dt = vtx_tol.max(rptline.tolerance);
                        if dist < dmin {
                            if dist <= dt {
                                ptline = rptline;
                                ivtx = j;
                                if surfacetype == SurfaceKind::Cone {
                                    ivtx = 0;
                                }
                            } else {
                                // cancel the previous solution because this
                                // point is better but its tolerance is not
                                // large enough.
                                ivtx = 0;
                            }
                            dmin = dist;
                        }
                    }

                    if ivtx != 0 && ptline.tolerance > vtx_tol {
                        vtx_tol = ptline.tolerance;
                        if !vtx_on_arc {
                            // now we should repeat the attempt to coincide with
                            // a bound of the arc (`cxx:944-956`).
                            let (c, va) = coincide_on_arc(&ptsommet, arc, surf, vtx_tol, domain);
                            vtx_on_arc = c;
                            vtxarc = va;
                            if vtx_on_arc {
                                let t = vtxarc
                                    .as_ref()
                                    .map(|v| tol_3d_vertex(v, domain, 0.0))
                                    .unwrap_or(0.0);
                                if t > vtx_tol {
                                    vtx_tol = t;
                                }
                            }
                        }
                    }

                    if typ_l == IntPatchIType::Walking {
                        let wlin = line
                            .as_wline()
                            .ok_or_else(|| "Standard_DomainError".to_string())?;
                        verify_tgline(wlin, paramline as i32, edge_tol, &mut tgline);
                    }

                    let (_ptbid, d1u, d1v) = surf.d1(u1, v1);
                    let tgrst = d1u
                        .multiplied_scalar(d2d.x())
                        .added(&d1v.multiplied_scalar(d2d.y()));

                    let normsurf = d1u.crossed(&d1v);
                    let (transline, _transarc) =
                        if normsurf.magnitude() < occt_core::precision::RESOLUTION {
                            let mut t1 = IntSurfTransition::new();
                            let mut t2 = IntSurfTransition::new();
                            t1.set_value_with_type(
                                true,
                                crate::intsurf_transition::IntSurfTypeTrans::Undecided,
                            );
                            t2.set_value_with_type(
                                true,
                                crate::intsurf_transition::IntSurfTypeTrans::Undecided,
                            );
                            (t1, t2)
                        } else {
                            let mut t1 = IntSurfTransition::new();
                            let mut t2 = IntSurfTransition::new();
                            let normale = GpDir::new(normsurf.x(), normsurf.y(), normsurf.z())
                                .map_err(|e| e.to_string())?;
                            make_transition(&tgline, &tgrst, &normale, &mut t1, &mut t2);
                            (t1, t2)
                        };

                    nb_treated += 1;

                    let update = VertexUpdate {
                        ptsommet,
                        vtx_tol,
                        u1,
                        v1,
                        u2,
                        v2,
                        param_line: paramline,
                        param_arc: paramarc,
                        trans_line: trans_type_of(&transline),
                        vtx_arc: vtxarc.clone(),
                        vtx_on_arc,
                        arc_id,
                        on_first,
                    };
                    if ivtx == 0 {
                        append_new_vertex(line, &update)?;
                        nbvtx += 1;
                    } else {
                        nbvtx =
                            update_existing_vertex(line, &update, ptline, ivtx, nbvtx, domain)?;
                    }
                    if nb_treated == 2 && typ_l == IntPatchIType::Walking {
                        // We processed a tangent zone and both ends have been
                        // treated, so mark the WLine as having the arc
                        // (`cxx:1229-1241`). `nb_treated` increments once per
                        // section point, so this arm is unreachable in OCCT
                        // too; it is kept verbatim.
                        if let Some(w) = line.as_wline_mut() {
                            if on_first {
                                w.set_arc_on_s1(arc_handle.clone());
                            } else {
                                w.set_arc_on_s2(arc_handle.clone());
                            }
                        }
                    }
                }

                indice_offset_periodic += 1;
                // `while (SurfaceIsPeriodic && IndiceOffsetPeriodic <= 2)`
                // (`cxx:1245`).
                if !(surface_is_periodic && indice_offset_periodic <= 2) {
                    break;
                }
            }
            indice_offset_bi_periodic += 1;
            // `while (SurfaceIsBiPeriodic && IndiceOffsetBiPeriodic <= 2)`
            // (`cxx:1248`).
            if !(surface_is_bi_periodic && indice_offset_bi_periodic <= 2) {
                break;
            }
        }

        domain.next();
    }

    // `cxx:1255-1354`: re-scale the line and recompute the vertex parameters.
    if typ_l == IntPatchIType::Walking {
        let mut pu1 = 0.0;
        let mut pv1 = 0.0;
        let mut pu2 = 0.0;
        let mut pv2 = 0.0;
        match type_s1 {
            SurfaceKind::Cylinder | SurfaceKind::Cone | SurfaceKind::Sphere => {
                pu1 = two_pi;
            }
            SurfaceKind::Torus => {
                pu1 = two_pi;
                pv1 = two_pi;
            }
            _ => {
                if surf1.is_u_periodic() {
                    pu1 = surf1.u_period();
                } else if surf1.is_u_closed() {
                    let (u0, u1) = surf1.u_range();
                    pu1 = u1 - u0;
                }
                if surf1.is_v_periodic() {
                    pv1 = surf1.v_period();
                } else if surf1.is_v_closed() {
                    let (v0, v1) = surf1.v_range();
                    pv1 = v1 - v0;
                }
            }
        }
        match type_s2 {
            SurfaceKind::Cylinder | SurfaceKind::Cone | SurfaceKind::Sphere => {
                pu2 = two_pi;
            }
            SurfaceKind::Torus => {
                pu2 = two_pi;
                pv2 = two_pi;
            }
            _ => {
                if surf2.is_u_periodic() {
                    pu2 = surf2.u_period();
                } else if surf2.is_u_closed() {
                    let (u0, u1) = surf2.u_range();
                    pu2 = u1 - u0;
                }
                if surf2.is_v_periodic() {
                    pv2 = surf2.v_period();
                } else if surf2.is_v_closed() {
                    let (v0, v1) = surf2.v_range();
                    pv2 = v1 - v0;
                }
            }
        }
        if let Some(w) = line.as_wline_mut() {
            w.set_period(pu1, pv1, pu2, pv2);
            intpatch_vparams::wline::compute_vertex_parameters(w, tol);
        }
    } else if let Some(r) = line.as_rline_mut() {
        intpatch_vparams::rline::compute_vertex_parameters_rline(r);
    }

    Ok(())
}

/// `IntSurf_Transition` -> the port's collapsed [`TransType`] slot (see
/// `int_tools_wline::PatchPoint`).
fn trans_type_of(t: &IntSurfTransition) -> TransType {
    match t.transition_type() {
        crate::intsurf_transition::IntSurfTypeTrans::In => TransType::In,
        crate::intsurf_transition::IntSurfTypeTrans::Out => TransType::Out,
        crate::intsurf_transition::IntSurfTypeTrans::Touch => TransType::Touch,
        crate::intsurf_transition::IntSurfTypeTrans::Undecided => TransType::Undecided,
    }
}
