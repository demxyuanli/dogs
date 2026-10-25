//! `ChFi3d_Builder::StoreData` — copy of an approximation result into a
//! `ChFiDS_SurfData`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_6.cxx:629-1005`
//! (OCCT 8.0.0).
//!
//! Divergences that are not heuristics:
//! - the member function becomes a free function over `&mut SurfData` and
//!   `&mut DataStructure` (`myDS->ChangeDS()`), so the caller keeps owning the
//!   builder; `S1` / `S2` are the port's [`SurfaceAdaptor`] and `S2` is an
//!   `Option` for the null `Adaptor3d_Surface` handle OCCT tests,
//! - `occ::handle<ChFiDS_SurfData>& Data` is `&mut SurfData`; the
//!   `ChangeVertex*` / `ChangeInterference*` references OCCT keeps open are
//!   re-read through [`SurfData::vertex_mut`] / [`SurfData::interference_mut`],
//! - the `const_cast` OCCT performs on `approx.SurfUKnots()` (`cxx:667-668`)
//!   becomes an owned `Vec` reparametrised in place; the approximate's own
//!   knot array is not modified, which is the only observable difference
//!   (`BRepBlend_AppSurface` keeps its arrays private),
//! - `ChFi3d_CheckSameParameter` takes a `GeomAdaptor_Curve` loaded with
//!   `(Crv3d, par1, par2)` (`cxx:657-658`, `:762`), i.e. the iso curve trimmed
//!   to the p-curve range; [`GeomTrimmedCurveBasis`] carries exactly that,
//! - `TopAbs_Face` `Orientation(FORWARD)` on a copy of the adaptor's face
//!   (`cxx:778-779`) does not change the p-curve lookup, but is reproduced,
//! - what OCCT reports by throwing `Standard_DomainError` from
//!   `ChFiDS_CommonPoint::Arc` / `Geom2dAdaptor_Curve::Value` is `Err` here;
//!   `Ok(false)` is the `return false` of the degenerated orientation loop.

use std::sync::Arc;

use occt_core::bspl::banded_interp::knot_sequence;
use occt_core::bspl::knots::reparametrize;
use occt_core::gp::{GpDir2d, GpPnt2d};
use occt_core::precision::CONFUSION;
use occt_geom::bspline_surface::{is_really_rational, GeomBSplineSurface};
use occt_geom::geom_lib_extend_surf::extend_bspline_surface;
use occt_geom::surface::Surface as GeomSurface;
use occt_geom::trimmed::GeomTrimmedCurveBasis;
use occt_geom2d::bspline_curve::Geom2dBSplineCurve;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;

use crate::abs::Orientation;
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_blend::app_surface::AppSurface;
use crate::brep_blend::line::Line;
use crate::chfi_ds::surf_data::SurfData;
use crate::shape::{Edge, Face};
use crate::topopebrep::ds::{Curve as DsCurve, DataStructure, Surface as DsSurface};
use crate::topopebrep::tool::adaptor::{CurveAdaptor, SurfaceAdaptor};

use crate::chfi3d::builder_helpers::{check_same_parameter, comp_param, fil_common_point, trsf_trans};

/// `gp::DY2d()` (`gp.cxx`): the direction of the `Geom2d_Line` OCCT puts on
/// the fillet surface (`cxx:753`, `:827`).
fn dy2d() -> GpDir2d {
    GpDir2d::new(0.0, 1.0).expect("DY2d")
}

/// `new Geom2d_BSplineCurve(approx.Curve2dPoles(I),
/// approx.Curves2dKnots(), approx.Curves2dMults(), approx.Curves2dDegree())`
/// (`cxx:757-759`, `:829-831`): the p-curve of the fillet on one face.
fn approx_curve_2d(approx: &AppSurface, index: i32) -> Option<Geom2dBSplineCurve> {
    let poles = approx.curve_2d_poles(index);
    let degree = approx.curves_2d_degree();
    if degree < 1 {
        return None;
    }
    let flat = knot_sequence(approx.curves_2d_knots(), approx.curves_2d_mults(), degree);
    let xs = poles.iter().map(|p| p.x()).collect();
    let ys = poles.iter().map(|p| p.y()).collect();
    Geom2dBSplineCurve::new(xs, ys, flat, degree as usize).ok()
}

/// A `TopoDS_Face` with the `FORWARD` orientation (`cxx:778-779`, `:794-795`,
/// `:845-846`, `:861-862`).
fn forward_face(f: &Face) -> Face {
    let mut face = f.clone();
    face.0.set_orientation(Orientation::Forward);
    face
}

/// `ChFi3d_CommonPoint::Arc()` as the `TopoDS_Edge` the adaptors load.
fn arc_edge(data: &SurfData, first: bool, on_s: i32) -> Result<Edge, String> {
    let arc = data
        .vertex(first, on_s)
        .arc()
        .ok_or_else(|| "ChFiDS_CommonPoint::Arc: Standard_DomainError".to_string())?;
    Edge::wrap(arc.clone()).ok_or_else(|| "ChFi3d_Builder::StoreData: arc is not an edge".to_string())
}

/// The `Gd` / `Gf` blocks of `StoreData`: `CompParam` on the arc of one common
/// point, then the tolerance and the arc it stores back (`cxx:775-806`,
/// `:842-873`).
#[allow(clippy::too_many_arguments)]
fn comp_arc_param(
    data: &mut SurfData,
    first: bool,
    on_s: i32,
    face: &Face,
    pcurve_on_face: &dyn Curve2d,
    tol_c: f64,
) -> Result<f64, String> {
    let edge = arc_edge(data, first, on_s)?;
    let mut brc = Curve2dAdaptor::new();
    brc.initialize(&edge, face);
    let carc = CurveAdaptor::of_edge(&edge);

    let (prefarc, preftg) = {
        let v = data.vertex(first, on_s);
        let prefarc = v
            .parameter_on_arc()
            .ok_or_else(|| "ChFiDS_CommonPoint::ParameterOnArc: Standard_DomainError".to_string())?;
        (prefarc, v.parameter())
    };
    let mut uarc = prefarc;
    let mut utg = preftg;
    comp_param(&brc, pcurve_on_face, &mut uarc, &mut utg, prefarc, preftg)?;

    let tolcheck = carc
        .d0(uarc)
        .unwrap_or_default()
        .distance(data.vertex(first, on_s).point());
    let tra = data
        .vertex(first, on_s)
        .transition_on_arc()
        .unwrap_or(Orientation::Forward);
    data.vertex_mut(first, on_s)
        .set_arc(tol_c + tolcheck, &edge.0, uarc, tra);
    Ok(utg)
}

/// `ChFi3d_Builder::StoreData(Data, approx, lin, S1, S2, Or1, Gd1, Gd2, Gf1,
/// Gf2, Reversed)` (`ChFi3d_Builder_6.cxx:629-1005`).
#[allow(clippy::too_many_arguments)]
pub fn store_data(
    data: &mut SurfData,
    approx: &AppSurface,
    lin: &Line,
    s1: &SurfaceAdaptor,
    s2: Option<&SurfaceAdaptor>,
    or1: Orientation,
    gd1: bool,
    gd2: bool,
    gf1: bool,
    gf2: bool,
    reversed: bool,
    ds: &mut DataStructure,
) -> Result<bool, String> {
    // Small control tools (`cxx:646-649`): one `GeomAdaptor_Curve` reused for
    // both sides, loaded with `(Crv3d, par1, par2)` below.
    let (tolget3d, tolget2d) = approx.tol_reached();
    let mut tol_c1 = tolget3d + approx.tol_curve_on_surf(1);
    let mut tol_c2 = 0.0;
    if s2.is_some() {
        tol_c2 = tolget3d + approx.tol_curve_on_surf(2);
    }

    // By default parametric space is created using a square surface
    // (`cxx:663-673`).
    let kv = approx.surf_v_knots();
    let larg = kv.last().copied().unwrap_or(0.0) - kv.first().copied().unwrap_or(0.0);
    let mut kku = approx.surf_u_knots().to_vec();
    reparametrize(0.0, larg, &mut kku);

    let udeg = approx.u_degree();
    let vdeg = approx.v_degree();
    if udeg < 1 || vdeg < 1 {
        return Err("Geom_BSplineSurface: degree must be >= 1".to_string());
    }
    let flat_u = knot_sequence(&kku, approx.surf_u_mults(), udeg);
    let flat_v = knot_sequence(kv, approx.surf_v_mults(), vdeg);
    let poles: Vec<Vec<occt_core::gp::GpPnt>> =
        approx.surf_poles().iter().map(|row| row.to_vec()).collect();
    // `new Geom_BSplineSurface(approx.SurfPoles(), approx.SurfWeights(), ...)`
    // (`cxx:669-676`): the constructor keeps the weight grid only when
    // `Rational()` finds it really rational (`Geom_BSplineSurface.cxx:110-138`,
    // `:270-272`), so a constant weight array stays non-rational.
    let weights = approx.surf_weights();
    let weight_grid = if is_really_rational(&weights) {
        Some(weights.iter().map(|row| row.to_vec()).collect())
    } else {
        None
    };
    let surf = GeomBSplineSurface::from_parts(
        poles,
        weight_grid,
        flat_u,
        flat_v,
        udeg as usize,
        vdeg as usize,
        false,
        false,
    );

    // Extension of the surface (`cxx:675-710`). The guards read the poles of
    // the *original* surface: OCCT's `ExtendSurfByLength` replaces the handle
    // instead of mutating it, so `Surf` still points at the unextended one
    // until the `down_cast` below.
    let nu0 = surf.nb_poles_u();
    let nv0 = surf.nb_poles_v();
    let p11 = surf.pole(0, 0);
    let p21 = surf.pole(nu0 - 1, 0);
    let p12 = surf.pole(0, nv0 - 1);
    let p22 = surf.pole(nu0 - 1, nv0 - 1);

    let length1 = data.first_extension_value();
    let length2 = data.last_extension_value();
    let eps = tolget3d.max(2.0 * CONFUSION);
    let mut ext1 = false;
    let mut ext2 = false;
    let mut bnd = surf.clone();
    if length1 > eps {
        // To avoid extending the surface with a singular boundary.
        if p11.distance(&p21) > eps {
            bnd = extend_bspline_surface(bnd, length1, 1, false, false)
                .ok_or_else(|| "GeomLib::ExtendSurfByLength failed".to_string())?;
            ext1 = true;
        }
    }
    if length2 > eps {
        if p12.distance(&p22) > eps {
            bnd = extend_bspline_surface(bnd, length2, 1, false, true)
                .ok_or_else(|| "GeomLib::ExtendSurfByLength failed".to_string())?;
            ext2 = true;
        }
    }

    // `Surf = down_cast<Geom_BSplineSurface>(aBndSurf)` (`cxx:710`).
    let mut surf = bnd;

    // Correction of the surface on the extremities (`cxx:712-728`).
    let (nu, nv) = (surf.nb_poles_u(), surf.nb_poles_v());
    if !ext1 {
        let p11 = lin.start_point_on_first().value();
        let p21 = lin.start_point_on_second().value();
        surf.set_pole(0, 0, *p11);
        surf.set_pole(nu - 1, 0, *p21);
    }
    if !ext2 {
        let p12 = lin.end_point_on_first().value();
        let p22 = lin.end_point_on_second().value();
        surf.set_pole(0, nv - 1, *p12);
        surf.set_pole(nu - 1, nv - 1, *p22);
    }

    // `Data->ChangeSurf(DStr.AddSurface(TopOpeBRepDS_Surface(Surf,
    // tolget3d)))` (`cxx:730`).
    let surf_index = ds.add_surface(&DsSurface::with_tolerance(
        Arc::new(surf.clone()) as Arc<dyn GeomSurface>,
        tolget3d,
    ));
    data.set_surf(surf_index);

    let (u_first, u_last) = surf.u_range();
    let (v_first, v_last) = surf.v_range();
    let mut uon1 = u_first;
    let mut uon2 = u_last;
    let mut ion1 = 1i32;
    let mut ion2 = 2i32;
    if reversed {
        uon1 = u_last;
        uon2 = u_first;
        ion1 = 2;
        ion2 = 1;
    }

    // The SurfData is filled in what concerns S1 (`cxx:750-821`).
    let crv3d1 = surf
        .u_iso_curve(uon1)
        .ok_or_else(|| "Geom_BSplineSurface::UIso failed".to_string())?;
    let pori1 = GpPnt2d::new(uon1, 0.0);
    let p_curve_on_surf = Arc::new(Geom2dLine::from_pnt_dir(pori1, dy2d())) as Arc<dyn Curve2d>;
    let p_curve_on_face = Arc::new(
        approx_curve_2d(approx, ion1)
            .ok_or_else(|| "BRepBlend_AppSurface::Curve2d failed".to_string())?,
    ) as Arc<dyn Curve2d>;
    let par1 = p_curve_on_face.first_parameter();
    let par2 = p_curve_on_face.last_parameter();
    let s1_surface = s1
        .surface()
        .ok_or_else(|| "BRepAdaptor_Surface: no surface".to_string())?;
    let s1_face = s1
        .face()
        .ok_or_else(|| "BRepAdaptor_Surface: no face".to_string())?
        .clone();
    let chc1 = GeomTrimmedCurveBasis::new(crv3d1.clone(), par1, par2);
    let (ok, tolcheck) = check_same_parameter(
        &chc1,
        p_curve_on_face.as_ref(),
        s1_surface.as_ref(),
        tol_c1,
    );
    if !ok {
        tol_c1 = tolcheck;
    }
    let index1_of_curve =
        ds.add_curve(&mut DsCurve::with_curve(Some(crv3d1.clone()), tol_c1, false));

    let mut pppdeb: f64;
    let mut pppfin: f64;
    if gd1 {
        let forwfac = forward_face(&s1_face);
        pppdeb = comp_arc_param(data, true, ion1, &forwfac, p_curve_on_face.as_ref(), tol_c1)?;
    } else {
        pppdeb = v_first;
    }
    if gf1 {
        let forwfac = forward_face(&s1_face);
        pppfin = comp_arc_param(data, false, ion1, &forwfac, p_curve_on_face.as_ref(), tol_c1)?;
    } else {
        pppfin = v_last;
    }
    let tra_on1 = trsf_trans(if reversed {
        lin.transition_on_s2()?
    } else {
        lin.transition_on_s1()?
    });
    {
        let fint1 = data.interference_mut(1);
        fint1.set_first_parameter(pppdeb);
        fint1.set_last_parameter(pppfin);
        fint1.set_interference(
            index1_of_curve,
            tra_on1,
            Some(p_curve_on_face.clone()),
            Some(p_curve_on_surf.clone()),
        );
    }

    // SurfData is filled in what concerns S2 (`cxx:823-896`).
    let crv3d2 = surf
        .u_iso_curve(uon2)
        .ok_or_else(|| "Geom_BSplineSurface::UIso failed".to_string())?;
    let pori2 = GpPnt2d::new(uon2, 0.0);
    let p_curve_on_surf = Arc::new(Geom2dLine::from_pnt_dir(pori2, dy2d())) as Arc<dyn Curve2d>;
    let mut p_curve_on_face2: Option<Arc<dyn Curve2d>> = None;
    if let Some(s2) = s2 {
        let pc = Arc::new(
            approx_curve_2d(approx, ion2)
                .ok_or_else(|| "BRepBlend_AppSurface::Curve2d failed".to_string())?,
        ) as Arc<dyn Curve2d>;
        let s2_surface = s2
            .surface()
            .ok_or_else(|| "BRepAdaptor_Surface: no surface".to_string())?;
        let chc2 = GeomTrimmedCurveBasis::new(crv3d2.clone(), par1, par2);
        let (ok, tolcheck) =
            check_same_parameter(&chc2, pc.as_ref(), s2_surface.as_ref(), tol_c2);
        if !ok {
            tol_c2 = tolcheck;
        }
        p_curve_on_face2 = Some(pc);
    }
    let index2_of_curve =
        ds.add_curve(&mut DsCurve::with_curve(Some(crv3d2.clone()), tol_c2, false));

    if gd2 {
        let s2 = s2.ok_or_else(|| "ChFi3d_Builder::StoreData: S2 is null".to_string())?;
        let forwfac = forward_face(
            s2.face()
                .ok_or_else(|| "BRepAdaptor_Surface: no face".to_string())?,
        );
        let pc = p_curve_on_face2
            .as_ref()
            .ok_or_else(|| "ChFi3d_Builder::StoreData: no p-curve on S2".to_string())?;
        pppdeb = comp_arc_param(data, true, ion2, &forwfac, pc.as_ref(), tol_c2)?;
    } else {
        pppdeb = v_first;
    }
    if gf2 {
        let s2 = s2.ok_or_else(|| "ChFi3d_Builder::StoreData: S2 is null".to_string())?;
        let forwfac = forward_face(
            s2.face()
                .ok_or_else(|| "BRepAdaptor_Surface: no face".to_string())?,
        );
        let pc = p_curve_on_face2
            .as_ref()
            .ok_or_else(|| "ChFi3d_Builder::StoreData: no p-curve on S2".to_string())?;
        pppfin = comp_arc_param(data, false, ion2, &forwfac, pc.as_ref(), tol_c2)?;
    } else {
        pppfin = v_last;
    }
    let fint2_first = pppdeb;
    let fint2_last = pppfin;
    if s2.is_some() {
        let tra_on2 = trsf_trans(if reversed {
            lin.transition_on_s1()?
        } else {
            lin.transition_on_s2()?
        });
        let fint2 = data.interference_mut(2);
        fint2.set_first_parameter(fint2_first);
        fint2.set_last_parameter(fint2_last);
        fint2.set_interference(
            index2_of_curve,
            tra_on2,
            p_curve_on_face2.clone(),
            Some(p_curve_on_surf.clone()),
        );
    } else {
        let fint2 = data.interference_mut(2);
        fint2.set_first_parameter(fint2_first);
        fint2.set_last_parameter(fint2_last);
        fint2.set_interference(
            index2_of_curve,
            Orientation::Forward,
            None,
            Some(p_curve_on_surf.clone()),
        );
    }

    // The orientation of the fillet in relation to the faces is evaluated
    // (`cxx:898-967`).
    let (sref, p_curve_on_face_orient): (&SurfaceAdaptor, Arc<dyn Curve2d>) = if reversed {
        let s2 = s2.ok_or_else(|| "ChFi3d_Builder::StoreData: S2 is null".to_string())?;
        (
            s2,
            data.interference(2)
                .p_curve_on_face()
                .cloned()
                .ok_or_else(|| "ChFiDS_FaceInterference: no p-curve on face".to_string())?,
        )
    } else {
        (
            s1,
            data.interference(1)
                .p_curve_on_face()
                .cloned()
                .ok_or_else(|| "ChFiDS_FaceInterference: no p-curve on face".to_string())?,
        )
    };
    let a_delta = v_last - v_first;
    let mut a_denom = 2.0f64;
    loop {
        let a_deltav = a_delta / a_denom;
        let a_param = v_first + a_deltav;
        let puv = p_curve_on_face_orient.d0(a_param);
        let (_, du1, dv1) = sref
            .d1(puv.x(), puv.y())
            .ok_or_else(|| "Adaptor3d_Surface::D1 failed".to_string())?;
        let mut du1 = du1.crossed(&dv1);
        if or1 == Orientation::Reversed {
            du1.reverse();
        }
        let (_, du2, dv2) = surf.d1(u_first, a_param);
        let du2 = du2.crossed(&dv2);
        if du1.magnitude() <= tolget3d || du2.magnitude() <= tolget3d {
            a_denom += 1.0;
            if a_deltav.abs() <= tolget2d {
                return Ok(false);
            }
            continue;
        }
        data.set_orientation(if du1.dot(&du2) > 0.0 {
            Orientation::Forward
        } else {
            Orientation::Reversed
        });
        break;
    }

    // `ChFi3d_FilCommonPoint` for the four extremities (`cxx:968-998`).
    if !gd1 {
        let mut cp = data.vertex(true, ion1).clone();
        fil_common_point(
            lin.start_point_on_first(),
            lin.transition_on_s1()?,
            true,
            &mut cp,
            tol_c1,
        );
        *data.vertex_mut(true, ion1) = cp;
    }
    if !gf1 {
        let mut cp = data.vertex(false, ion1).clone();
        fil_common_point(
            lin.end_point_on_first(),
            lin.transition_on_s1()?,
            false,
            &mut cp,
            tol_c1,
        );
        *data.vertex_mut(false, ion1) = cp;
    }
    if !gd2 && s2.is_some() {
        let mut cp = data.vertex(true, ion2).clone();
        fil_common_point(
            lin.start_point_on_second(),
            lin.transition_on_s2()?,
            true,
            &mut cp,
            tol_c2,
        );
        *data.vertex_mut(true, ion2) = cp;
    }
    if !gf2 && s2.is_some() {
        let mut cp = data.vertex(false, ion2).clone();
        fil_common_point(
            lin.end_point_on_second(),
            lin.transition_on_s2()?,
            false,
            &mut cp,
            tol_c2,
        );
        *data.vertex_mut(false, ion2) = cp;
    }

    // Parameters on ElSpine (`cxx:999-1003`).
    let nbp = lin.nb_points();
    let first = lin
        .point(1)
        .ok_or_else(|| "BRepBlend_Line::Point: Standard_OutOfRange".to_string())?;
    let last = lin
        .point(nbp)
        .ok_or_else(|| "BRepBlend_Line::Point: Standard_OutOfRange".to_string())?;
    data.set_first_spine_param(first.parameter());
    data.set_last_spine_param(last.parameter());
    Ok(true)
}
