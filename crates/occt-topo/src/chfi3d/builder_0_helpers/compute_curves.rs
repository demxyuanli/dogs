//! `ChFi3d_ComputeCurves` — the intersection curve of two surfaces, with a
//! p-curve on each, used by the "cut the fillet by a plane" extremity.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:3694-4295`
//! (OCCT 8.0.0).
//!
//! Arms, in source order:
//!
//! - cylinder ∩ plane (`cxx:3746-3864`): the `IntAna_QuadQuadGeo` analytic
//!   result, reoriented along `Pdeb -> Pfin`, p-curves by
//!   `ChFi3d_ProjectPCurv`, and a UV translation when a surface is a cylinder,
//! - plane ∩ plane (`cxx:3865-3887`): the analytic line,
//! - otherwise `GeomInt_IntSS` (`cxx:3889-4070`), trimmed to the two known
//!   extremities unless `whole_curv`,
//! - and, when that fails, the `IntWalk_PWalking` + `GeomInt_WLApprox` path
//!   (`cxx:4072-4294`), in [`super::compute_curves_walking`].
//!
//! Translation notes:
//! - `NCollection_Array1<double> Pardeb / Parfin` (indices `1..4`) become
//!   `[f64; 4]` counted from `0`,
//! - `C3d`, `Pc1`, `Pc2` and `tolreached` are OCCT out-parameters; the port
//!   answers them in [`ComputeCurves`],
//! - OCCT closes the analytic arms with
//!   `C3d = new Geom_TrimmedCurve(C3d, Udeb, Ufin)`. The port keeps the
//!   *untrimmed* analytic curve plus the `(Udeb, Ufin)` range in
//!   [`ComputeCurves::pardeb`] / [`ComputeCurves::parfin`]: `Geom_TrimmedCurve`
//!   evaluates its basis at the basis parameters (`Geom_TrimmedCurve.cxx:214`,
//!   `:258-267`) and only reports the trim as its range, which is exactly the
//!   pair the port carries. It is also the convention `ChFi3d_ComputeArete`
//!   already uses (untrimmed curve + `Pardeb` / `Parfin`), and the p-curves
//!   below are built on basis parameters either way,
//! - `trsfsurf` (`cxx:3600-3656`) is [`trsfsurf`]: the adaptor hands out the
//!   `Geom_Surface` it already holds, so the `BRepAdaptor_Surface` transform
//!   arm of OCCT's `down_cast` is folded into [`SurfaceAdaptor::surface`],
//! - `GeomAPI_ProjectPointOnCurve` is [`point_curve_extrema_all`].

use std::sync::Arc;

use occt_core::elib::clib;
use occt_core::gp::{GpLin, GpPnt, GpTrsf2d, GpVec, GpVec2d};
use occt_core::precision::{ANGULAR, PCONFUSION, CONFUSION};
use occt_geom::curve::Curve;
use occt_geom::extrema::ExtremaPair;
use occt_geom::extrema_pc::point_curve_extrema_all;
use occt_geom::intana::{
    quadric_quadric_plane_cylinder, quadric_quadric_planes, QuadricIntersection,
};
use occt_geom::surface::Surface;
use occt_geom::{GeomCircle, GeomEllipse, GeomLine};
use occt_geom::{GeomRectangularTrimmedSurface};
use occt_geom2d::curve::Curve2d;

use crate::chfi3d::builder_0_helpers::eval_tol_reached::ch_fi3d_eval_tol_reached;
use crate::chfi3d::builder_0_helpers::project_pcurv::ch_fi3d_project_pcurv;
use crate::chfi3d::builder_0_helpers::reparam_pcurv::ch_fi3d_reparam_pcurv;
use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::geom_int::IntSS;
use crate::topopebrep::tool::adaptor::{GeomAdaptorCurve, SurfaceAdaptor};

/// `tolap` of the `GeomInt_IntSS` call (`cxx:3901`).
const TOLAP: f64 = 2.0e-7;

/// The out-parameters of [`ch_fi3d_compute_curves`].
#[derive(Default, Clone)]
pub struct ComputeCurves {
    /// `C3d`.
    pub c3d: Option<Arc<dyn Curve>>,
    /// `Pc1`.
    pub pc1: Option<Arc<dyn Curve2d>>,
    /// `Pc2`.
    pub pc2: Option<Arc<dyn Curve2d>>,
    /// The first parameter of the `Geom_TrimmedCurve` OCCT would build around
    /// [`Self::c3d`]: `Udeb` in the analytic arms, `Uf` in the `GeomInt_IntSS`
    /// arm.
    pub pardeb: f64,
    /// The last parameter of that trimmed curve.
    pub parfin: f64,
    /// `tolreached`.
    pub tolreached: f64,
}

/// `trsfsurf(HS, dom)` (`cxx:3600-3656`): the `Geom_Surface` under the
/// adaptor, without a `Geom_RectangularTrimmedSurface` layer and trimmed to the
/// adaptor's parameter box, clamped to the surface's own bounds when that
/// parameter is not periodic.
pub fn trsfsurf(hs: &SurfaceAdaptor) -> Option<Arc<dyn Surface>> {
    let mut res = hs.surface()?.clone();
    // `occ::down_cast<Geom_RectangularTrimmedSurface>(res)` then
    // `res = tr->BasisSurface();` (`cxx:3616-3620`).
    if let Some(basis) = res.rectangular_trimmed_basis() {
        res = basis;
    }
    // `double U1 = HS->FirstUParameter(), ...` (`cxx:3622-3625`).
    let mut u1 = hs.first_u_parameter();
    let mut u2 = hs.last_u_parameter();
    let mut v1 = hs.first_v_parameter();
    let mut v2 = hs.last_v_parameter();
    // `cxx:3628-3651`.
    let (su1, su2) = res.u_range();
    let (sv1, sv2) = res.v_range();
    if !res.is_u_periodic() {
        if u1 < su1 {
            u1 = su1;
        }
        if u2 > su2 {
            u2 = su2;
        }
    }
    if !res.is_v_periodic() {
        if v1 < sv1 {
            v1 = sv1;
        }
        if v2 > sv2 {
            v2 = sv2;
        }
    }
    // `res = new Geom_RectangularTrimmedSurface(res, U1, U2, V1, V2);`
    // (`cxx:3652`). A degenerate box raises `Standard_ConstructionError` in
    // OCCT; the port answers the untrimmed surface instead of throwing.
    Some(match GeomRectangularTrimmedSurface::uv(res.clone(), u1, u2, v1, v2) {
        Ok(trimmed) => Arc::new(trimmed),
        Err(_) => res,
    })
}

/// `ChFi3d_ComputeCurves(S1, S2, Pardeb, Parfin, C3d, Pc1, Pc2, tol3d, tol2d,
/// tolreached, wholeCurv)` (`cxx:3694-4295`).
///
/// Answers `false` where OCCT returns `false` (no usable intersection) and
/// fills `out` where it returns `true`.
#[allow(clippy::too_many_arguments)]
pub fn ch_fi3d_compute_curves(
    s1: &SurfaceAdaptor,
    s2: &SurfaceAdaptor,
    pardeb: [f64; 4],
    parfin: [f64; 4],
    out: &mut ComputeCurves,
    tol3d: f64,
    tol2d: f64,
    whole_curv: bool,
) -> bool {
    let value = |s: &SurfaceAdaptor, u: f64, v: f64| s.value(u, v).unwrap_or_default();
    // `cxx:3706-3709`.
    let pdeb1 = value(s1, pardeb[0], pardeb[1]);
    let pfin1 = value(s1, parfin[0], parfin[1]);
    let pdeb2 = value(s2, pardeb[2], pardeb[3]);
    let pfin2 = value(s2, parfin[2], parfin[3]);

    // `cxx:3711-3734`: `distrefdeb`, `distreffin`, `pdeb`, `pfin`, `distref`.
    let mut distrefdeb = pdeb1.distance(&pdeb2);
    let mut distreffin = pfin1.distance(&pfin2);
    if distrefdeb < tol3d {
        distrefdeb = tol3d;
    }
    if distreffin < tol3d {
        distreffin = tol3d;
    }
    let pdeb = midpoint(&pdeb1, &pdeb2);
    let pfin = midpoint(&pfin1, &pfin2);
    let mut distref = 0.005 * pdeb.distance(&pfin);
    if distref < distrefdeb {
        distref = distrefdeb;
    }
    if distref < distreffin {
        distref = distreffin;
    }

    // `gp_Vec Vint, Vref(pdeb, pfin); double tolr1, tolr2; tolr1 = tolr2 =
    // tolreached = tol3d;` (`cxx:3741-3746`).
    let vref = GpVec::from_pnts(&pdeb, &pfin);
    out.tolreached = tol3d;
    let mut tolr1 = tol3d;
    let mut tolr2 = tol3d;

    let t1 = s1.get_type();
    let t2 = s2.get_type();
    if (t1 == GeomAbsSurfaceType::Cylinder && t2 == GeomAbsSurfaceType::Plane)
        || (t1 == GeomAbsSurfaceType::Plane && t2 == GeomAbsSurfaceType::Cylinder)
    {
        if let Some(result) = cylinder_plane_arm(
            s1, s2, t1, t2, &pdeb, &pfin, &vref, pardeb, tol3d, tol2d, out, &mut tolr1,
            &mut tolr2,
        ) {
            return result;
        }
    } else if t1 == GeomAbsSurfaceType::Plane && t2 == GeomAbsSurfaceType::Plane {
        if let Some(result) = plane_plane_arm(s1, s2, &pdeb, &pfin, &vref, tol3d, out) {
            return result;
        }
    } else {
        // `here GeomInt is approached.` (`cxx:3890`).
        if int_ss_arm(
            s1, s2, &pdeb, &pfin, &pdeb1, &pdeb2, &pfin1, &pfin2, pardeb, parfin, distref,
            distrefdeb, distreffin, tol2d, whole_curv, out,
        ) {
            return true;
        }
    }

    // `cxx:4072-4294`: the walking path.
    super::compute_curves_walking::walking_arm(
        s1, s2, &pdeb, &pfin, pardeb, parfin, tol3d, tol2d, out,
    )
}

/// `cxx:3746-3864`: cylinder ∩ plane.
///
/// `Some(result)` when the analytic intersection was usable, `None` when OCCT
/// would fall through to the next arm (`isIntDone` false, i.e. the `if
/// (isIntDone)` block was skipped and the `else` chain starts at `cxx:3889`).
#[allow(clippy::too_many_arguments)]
fn cylinder_plane_arm(
    s1: &SurfaceAdaptor,
    s2: &SurfaceAdaptor,
    t1: GeomAbsSurfaceType,
    t2: GeomAbsSurfaceType,
    pdeb: &GpPnt,
    pfin: &GpPnt,
    vref: &GpVec,
    pardeb: [f64; 4],
    tol3d: f64,
    tol2d: f64,
    out: &mut ComputeCurves,
    tolr1: &mut f64,
    tolr2: &mut f64,
) -> Option<bool> {
    // `cxx:3747-3761`.
    let (pl, cyl) = if t1 == GeomAbsSurfaceType::Plane {
        (
            s1.surface().and_then(|s| s.gp_pln())?,
            s2.surface().and_then(|s| s.gp_cylinder())?,
        )
    } else {
        (
            s2.surface().and_then(|s| s.gp_pln())?,
            s1.surface().and_then(|s| s.gp_cylinder())?,
        )
    };
    // `IntAna_QuadQuadGeo ImpKK(pl, cyl, Precision::Angular(), tol3d);`
    // (`cxx:3762`); that constructor leaves `H` at `0` (`IntAna_IntConicQuad`).
    let inter = quadric_quadric_plane_cylinder(&pl, &cyl, ANGULAR, tol3d, 0.0);

    // `bool isIntDone = ImpKK.IsDone(); if (TypeInter() == IntAna_Ellipse)
    // { isIntDone = (aMajorR < 100000.0 * aMinorR); }` (`cxx:3763-3770`).
    let mut is_int_done = !matches!(
        inter,
        QuadricIntersection::None | QuadricIntersection::NotDone
    );
    if let QuadricIntersection::Ellipse(e) = &inter {
        is_int_done = e.major_radius() < 100000.0 * e.minor_radius();
    }
    if !is_int_done {
        return None;
    }

    // `cxx:3775-3812`: the `switch (ImpKK.TypeInter())` arms.
    let mut c1line = false;
    let mut udeb;
    let mut ufin;
    let mut vint = GpVec::default();
    let c3d: Arc<dyn Curve> = match &inter {
        QuadricIntersection::Line(_) | QuadricIntersection::TwoLines(..) => {
            c1line = true;
            let candidates: Vec<GpLin> = match &inter {
                QuadricIntersection::Line(l) => vec![*l],
                QuadricIntersection::TwoLines(a, b) => vec![*a, *b],
                _ => unreachable!(),
            };
            // `for (int ilin = 1; ilin <= nbsol; ilin++) { C1 =
            // ImpKK.Line(ilin); Udeb = ElCLib::Parameter(C1, pdeb); ... if
            // (ptest.Distance(pdeb) < tol3d) break; }` (`cxx:3780-3790`). With
            // no early break `C1` keeps the last candidate.
            let mut chosen = candidates[0];
            let mut udeb_hit = 0.0;
            for cand in &candidates {
                chosen = *cand;
                udeb_hit = clib::line_parameter(&chosen, pdeb);
                if clib::line_value(&chosen, udeb_hit).distance(pdeb) < tol3d {
                    break;
                }
            }
            udeb = udeb_hit;
            ufin = clib::line_parameter(&chosen, pfin);
            vint = clib::line_d1(&chosen, udeb).1;
            Arc::new(GeomLine::new(chosen))
        }
        QuadricIntersection::Circle(c) => {
            udeb = clib::circle_parameter(&c.position(), pdeb);
            ufin = clib::circle_parameter(&c.position(), pfin);
            vint = clib::circle_d1(c, udeb).1;
            Arc::new(GeomCircle::new(c.clone()))
        }
        QuadricIntersection::Ellipse(e) => {
            udeb = clib::parameter_elips(e, pdeb);
            ufin = clib::parameter_elips(e, pfin);
            vint = clib::ellipse_d1(e, udeb).1;
            Arc::new(GeomEllipse::new(*e))
        }
        // `default: break;` (`cxx:3810-3811`): `C3d` stays null and the
        // `HC->Load` two lines below raises in OCCT.
        _ => return Some(false),
    };

    // `if (Vint.Dot(Vref) < 0) { C3d->Reverse(); if (c1line) {Udeb = -Udeb;
    // Ufin = -Ufin;} else {Udeb = 2*PI - Udeb; Ufin = 2*PI - Ufin;} }`
    // (`cxx:3813-3828`).
    let mut c3d = c3d;
    if vint.dot(vref) < 0.0 {
        c3d = Arc::from(c3d.reversed());
        if c1line {
            udeb = -udeb;
            ufin = -ufin;
        } else {
            udeb = 2.0 * std::f64::consts::PI - udeb;
            ufin = 2.0 * std::f64::consts::PI - ufin;
        }
    }
    // `if (!c1line) ElCLib::AdjustPeriodic(0., 2*PI, Precision::Angular(),
    // Udeb, Ufin);` (`cxx:3829-3832`).
    if !c1line {
        clib::adjust_periodic(
            0.0,
            2.0 * std::f64::consts::PI,
            ANGULAR,
            &mut udeb,
            &mut ufin,
        );
    }

    // `GeomAdaptor_Curve HC; HC.Load(C3d, Udeb, Ufin);` (`cxx:3833-3835`).
    let adapted: Arc<dyn Curve> = Arc::new(GeomAdaptorCurve::with_range(c3d.clone(), udeb, ufin));
    let s1_geom = s1.surface().cloned()?;
    let s2_geom = s2.surface().cloned()?;

    // `ChFi3d_ProjectPCurv(HC, S1, Pc1, tol3d, tolr1);` (`cxx:3836`).
    match ch_fi3d_project_pcurv(&adapted, &s1_geom, &mut out.pc1, tol3d) {
        Ok(Some(t)) => *tolr1 = t,
        Ok(None) => {}
        Err(_) => return Some(false),
    }
    // `if (S1->GetType() == GeomAbs_Cylinder) { Pc1->Value(Udeb).Coord(x, y);
    // x = Pardeb(1) - x; y = Pardeb(2) - y; if (fabs(x) >= tol2d || fabs(y) >=
    // tol2d) Pc1->Translate(gp_Vec2d(x, y)); }` (`cxx:3838-3846`).
    if t1 == GeomAbsSurfaceType::Cylinder {
        if let Some(pc1) = out.pc1.as_ref() {
            let p2d = pc1.value(udeb);
            let x = pardeb[0] - p2d.x();
            let y = pardeb[1] - p2d.y();
            if x.abs() >= tol2d || y.abs() >= tol2d {
                out.pc1 = Some(translate_2d(pc1, x, y));
            }
        }
    }
    // `ChFi3d_ProjectPCurv(HC, S2, Pc2, tol3d, tolr2);` (`cxx:3847`).
    match ch_fi3d_project_pcurv(&adapted, &s2_geom, &mut out.pc2, tol3d) {
        Ok(Some(t)) => *tolr2 = t,
        Ok(None) => {}
        Err(_) => return Some(false),
    }
    // `cxx:3848-3857`.
    if t2 == GeomAbsSurfaceType::Cylinder {
        if let Some(pc2) = out.pc2.as_ref() {
            let p2d = pc2.value(udeb);
            let x = pardeb[2] - p2d.x();
            let y = pardeb[3] - p2d.y();
            if x.abs() >= tol2d || y.abs() >= tol2d {
                out.pc2 = Some(translate_2d(pc2, x, y));
            }
        }
    }

    // `C3d = new Geom_TrimmedCurve(C3d, Udeb, Ufin);` (`cxx:3858`), kept as
    // the pair (untrimmed curve, range).
    out.c3d = Some(c3d);
    out.pardeb = udeb;
    out.parfin = ufin;
    // `tolreached = 1.5 * max(tolr1, tolr2); tolreached =
    // min(tolreached, ChFi3d_EvalTolReached(S1, Pc1, S2, Pc2, C3d));`
    // (`cxx:3859-3861`).
    let mut tolreached = 1.5 * tolr1.max(*tolr2);
    if let (Some(pc1), Some(pc2), Some(c3d_now)) =
        (out.pc1.clone(), out.pc2.clone(), out.c3d.clone())
    {
        let sampled: Arc<dyn Curve> =
            Arc::new(GeomAdaptorCurve::with_range(c3d_now, udeb, ufin));
        tolreached = tolreached.min(ch_fi3d_eval_tol_reached(
            s1_geom.as_ref(),
            pc1.as_ref(),
            s2_geom.as_ref(),
            pc2.as_ref(),
            sampled.as_ref(),
        ));
    }
    out.tolreached = tolreached;
    Some(true)
}

/// `cxx:3865-3887`: plane ∩ plane.
///
/// `None` when `LInt.IsDone()` is false (OCCT then leaves the analytic chain
/// and reaches the `GeomInt_IntSS` arm), `Some(result)` otherwise.
fn plane_plane_arm(
    s1: &SurfaceAdaptor,
    s2: &SurfaceAdaptor,
    pdeb: &GpPnt,
    pfin: &GpPnt,
    vref: &GpVec,
    tol3d: f64,
    out: &mut ComputeCurves,
) -> Option<bool> {
    // `IntAna_QuadQuadGeo LInt(S1->Plane(), S2->Plane(),
    // Precision::Angular(), tol3d);` (`cxx:3867-3868`).
    let pl1 = s1.surface().and_then(|s| s.gp_pln())?;
    let pl2 = s2.surface().and_then(|s| s.gp_pln())?;
    let inter = quadric_quadric_planes(&pl1, &pl2, ANGULAR, tol3d);
    // `if (LInt.IsDone())` (`cxx:3869`).
    let QuadricIntersection::Line(l) = inter else {
        return None;
    };
    // `cxx:3870-3879`.
    let mut udeb = clib::line_parameter(&l, pdeb);
    let mut ufin = clib::line_parameter(&l, pfin);
    let vint = clib::line_d1(&l, udeb).1;
    let mut c3d: Arc<dyn Curve> = Arc::new(GeomLine::new(l));
    if vint.dot(vref) < 0.0 {
        c3d = Arc::from(c3d.reversed());
        udeb = -udeb;
        ufin = -ufin;
    }
    // `cxx:3880-3887`.
    let adapted: Arc<dyn Curve> = Arc::new(GeomAdaptorCurve::with_range(c3d.clone(), udeb, ufin));
    let s1_geom = s1.surface().cloned()?;
    let s2_geom = s2.surface().cloned()?;
    match ch_fi3d_project_pcurv(&adapted, &s1_geom, &mut out.pc1, tol3d) {
        Ok(_) => {}
        Err(_) => return Some(false),
    }
    match ch_fi3d_project_pcurv(&adapted, &s2_geom, &mut out.pc2, tol3d) {
        Ok(_) => {}
        Err(_) => return Some(false),
    }
    // `C3d = new Geom_TrimmedCurve(C3d, Udeb, Ufin); return true;`
    // (`cxx:3888-3889`). OCCT leaves `tolreached` untouched here.
    out.c3d = Some(c3d);
    out.pardeb = udeb;
    out.parfin = ufin;
    Some(true)
}

/// `cxx:3889-4070`: the `GeomInt_IntSS` arm.
#[allow(clippy::too_many_arguments)]
fn int_ss_arm(
    s1: &SurfaceAdaptor,
    s2: &SurfaceAdaptor,
    pdeb: &GpPnt,
    pfin: &GpPnt,
    pdeb1: &GpPnt,
    pdeb2: &GpPnt,
    pfin1: &GpPnt,
    pfin2: &GpPnt,
    pardeb: [f64; 4],
    parfin: [f64; 4],
    distref: f64,
    distrefdeb: f64,
    distreffin: f64,
    tol2d: f64,
    whole_curv: bool,
    out: &mut ComputeCurves,
) -> bool {
    // `trsfsurf(S1, dom1)` / `trsfsurf(S2, dom2)` (`cxx:3892-3893`), then
    // `if (!gs1.IsNull() && !gs2.IsNull())` (`cxx:3896`).
    let (Some(gs1), Some(gs2)) = (trsfsurf(s1), trsfsurf(s2)) else {
        return false;
    };

    // `double tolap = 2.e-7; inter.Perform(gs1, gs2, tolap, true, true, true);`
    // (`cxx:3901-3903`).
    let mut inter = IntSS::new();
    inter.perform(gs1, gs2, TOLAP, true, true, true);
    if !inter.is_done() {
        return false;
    }

    // `for (int ilin = 1; ilin <= nbl; ilin++)` (`cxx:3928-4068`).
    for ilin in 1..=inter.nb_lines() {
        // `if (inter.HasLineOnS1(ilin) && inter.HasLineOnS2(ilin))`
        // (`cxx:3930`).
        if !(inter.has_line_on_s1(ilin) && inter.has_line_on_s2(ilin)) {
            continue;
        }
        let (Some(line), Some(pc1), Some(pc2)) = (
            inter.line(ilin).cloned(),
            inter.line_on_s1(ilin).cloned(),
            inter.line_on_s2(ilin).cloned(),
        ) else {
            continue;
        };
        let mut c3d = line.curve;
        let mut pc1 = pc1;
        let mut pc2 = pc2;

        // `cxx:3931-4026`: the two end parameters.
        let (uf, ul) = if whole_curv {
            // `Uf = C3d->FirstParameter(); Ul = C3d->LastParameter();`
            // (`cxx:3939-3940`).
            (c3d.first_parameter(), c3d.last_parameter())
        } else {
            let mut failed_f = false;
            let mut failed_l = false;
            let mut uf = 0.0;
            let mut ul = 0.0;

            // `proj.Init(pdeb1, C3d); if (proj.NbPoints() == 0 &&
            // distrefdeb > Precision::Confusion()) proj.Perform(pdeb2);`
            // (`cxx:3949-3953`).
            let mut proj = point_curve_extrema_all(c3d.as_ref(), pdeb1);
            if proj.is_empty() && distrefdeb > CONFUSION {
                proj = point_curve_extrema_all(c3d.as_ref(), pdeb2);
            }
            if proj.is_empty() {
                failed_f = true;
            } else {
                uf = nearest(&proj).u2;
            }
            // `proj.Perform(pfin1); if (proj.NbPoints() == 0 && distreffin >
            // Precision::Confusion()) proj.Perform(pfin2);` (`cxx:3961-3965`).
            let mut proj = point_curve_extrema_all(c3d.as_ref(), pfin1);
            if proj.is_empty() && distreffin > CONFUSION {
                proj = point_curve_extrema_all(c3d.as_ref(), pfin2);
            }
            if proj.is_empty() {
                failed_l = true;
            } else {
                ul = nearest(&proj).u2;
            }

            if failed_f && failed_l {
                // `cxx:3975-3979`.
                uf = c3d.first_parameter();
                ul = c3d.last_parameter();
            } else if failed_f || failed_l {
                // `cxx:3980-3994`.
                let uok = if failed_f { ul } else { uf };
                let u1 = c3d.first_parameter();
                let u2 = c3d.last_parameter();
                let uok = if (uok - u1).abs() > (uok - u2).abs() {
                    u1
                } else {
                    u2
                };
                if failed_f {
                    uf = uok;
                } else {
                    ul = uok;
                }
            } else if (uf - ul).abs() < PCONFUSION {
                // `cxx:3996-4002`.
                continue;
            }

            // `cxx:4004-4022`: `C3d->IsPeriodic() && !(failedF && failedL)`.
            if c3d.is_periodic() && !(failed_f && failed_l) {
                let ptestdeb = c3d.d0(uf);
                let ptestfin = c3d.d0(ul);
                let dir_new = GpVec::from_pnts(&ptestdeb, &ptestfin);
                let dir_old = c3d.d1(uf).1;
                if dir_old.dot(&dir_new) < 0.0 {
                    std::mem::swap(&mut uf, &mut ul);
                }
            }
            (uf, ul)
        };

        // `C3d = new Geom_TrimmedCurve(C3d, Uf, Ul); Pc1 = new
        // Geom2d_TrimmedCurve(Pc1, Uf, Ul); Pc2 = new Geom2d_TrimmedCurve(Pc2,
        // Uf, Ul);` (`cxx:4023-4025`). The port keeps the untrimmed curve and
        // reports `(Uf, Ul)`; the p-curves keep their own (basis) parameters,
        // which is what `Geom2d_TrimmedCurve` does in OCCT as well.
        let mut ptestdeb = c3d.d0(uf);
        let mut ptestfin = c3d.d0(ul);

        // `cxx:4027-4042`: is it necessary to invert?
        let mut dist_deb_to_deb = ptestdeb.distance(pdeb);
        let dist_deb_to_fin = ptestdeb.distance(pfin);
        let mut dist_fin_to_fin = ptestfin.distance(pfin);
        let dist_fin_to_deb = ptestfin.distance(pdeb);
        if dist_deb_to_deb > dist_deb_to_fin && dist_fin_to_fin > dist_fin_to_deb {
            c3d = Arc::from(c3d.reversed());
            pc1 = Arc::from(pc1.reversed());
            pc2 = Arc::from(pc2.reversed());
            ptestdeb = c3d.d0(uf);
            ptestfin = c3d.d0(ul);
            dist_deb_to_deb = ptestdeb.distance(pdeb);
            dist_fin_to_fin = ptestfin.distance(pfin);
        }
        // `if (DistDebToDeb < distref && DistFinToFin < distref)`
        // (`cxx:4043`).
        if !(dist_deb_to_deb < distref && dist_fin_to_fin < distref) {
            continue;
        }

        // `Uf = C3d->FirstParameter(); Ul = C3d->LastParameter();
        // ChFi3d_ReparamPcurv(Uf, Ul, Pc1); ChFi3d_ReparamPcurv(Uf, Ul,
        // Pc2);` (`cxx:4044-4048`).
        let uf = c3d.first_parameter();
        let ul = c3d.last_parameter();
        let mut pc1_opt = Some(pc1);
        let mut pc2_opt = Some(pc2);
        if ch_fi3d_reparam_pcurv(uf, ul, &mut pc1_opt).is_err() {
            return false;
        }
        if ch_fi3d_reparam_pcurv(uf, ul, &mut pc2_opt).is_err() {
            return false;
        }
        // `cxx:4049-4062`: move each p-curve onto `Pardeb` at `Uf`.
        pc1_opt = shift_at(pc1_opt, uf, pardeb[0], pardeb[1], tol2d);
        pc2_opt = shift_at(pc2_opt, uf, pardeb[2], pardeb[3], tol2d);

        // `tolreached = ChFi3d_EvalTolReached(S1, Pc1, S2, Pc2, C3d);`
        // (`cxx:4063`).
        let (Some(s1_geom), Some(s2_geom)) = (s1.surface().cloned(), s2.surface().cloned()) else {
            return false;
        };
        let (Some(pc1), Some(pc2)) = (pc1_opt.clone(), pc2_opt.clone()) else {
            return false;
        };
        let sampled: Arc<dyn Curve> = Arc::new(GeomAdaptorCurve::with_range(c3d.clone(), uf, ul));
        out.tolreached = ch_fi3d_eval_tol_reached(
            s1_geom.as_ref(),
            pc1.as_ref(),
            s2_geom.as_ref(),
            pc2.as_ref(),
            sampled.as_ref(),
        );
        out.c3d = Some(c3d);
        out.pc1 = pc1_opt;
        out.pc2 = pc2_opt;
        out.pardeb = uf;
        out.parfin = ul;
        return true;
    }
    false
}

/// The `Pc->Value(Uf).Coord(x, y); x = Pardeb(i) - x; ... Pc->Translate(...)`
/// couple of `cxx:4049-4062`.
fn shift_at(
    pc: Option<Arc<dyn Curve2d>>,
    uf: f64,
    target_u: f64,
    target_v: f64,
    tol2d: f64,
) -> Option<Arc<dyn Curve2d>> {
    let pc = pc?;
    let p2d = pc.value(uf);
    let x = target_u - p2d.x();
    let y = target_v - p2d.y();
    if x.abs() > tol2d || y.abs() > tol2d {
        Some(translate_2d(&pc, x, y))
    } else {
        Some(pc)
    }
}

fn midpoint(a: &GpPnt, b: &GpPnt) -> GpPnt {
    GpPnt::new(
        0.5 * (a.x() + b.x()),
        0.5 * (a.y() + b.y()),
        0.5 * (a.z() + b.z()),
    )
}

/// `proj.LowerDistanceParameter()` (`cxx:3960`, `:3971`).
fn nearest(solutions: &[ExtremaPair]) -> &ExtremaPair {
    let mut best = 0usize;
    for i in 1..solutions.len() {
        if solutions[i].distance < solutions[best].distance {
            best = i;
        }
    }
    &solutions[best]
}

/// `Pcurv->Translate(gp_Vec2d(x, y))`: the port's `Curve2d` mutates through a
/// fresh object.
fn translate_2d(curve: &Arc<dyn Curve2d>, x: f64, y: f64) -> Arc<dyn Curve2d> {
    let mut moved = curve.as_ref().clone_dyn();
    let mut t = GpTrsf2d::identity();
    t.set_translation_vec(&GpVec2d::new(x, y));
    moved.transform(&t);
    Arc::from(moved)
}
