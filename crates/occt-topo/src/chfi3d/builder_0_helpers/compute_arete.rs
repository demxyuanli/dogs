//! `ChFi3d_ComputeArete` — the 3D curve and the p-curve of one edge of a
//! fillet, between two common points.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:1969-2152`
//! (OCCT 8.0.0).
//!
//! The three arms, in source order:
//!
//! - `|UV1.X - UV2.X| <= tol2d`: the p-curve follows a U-iso of `Surf`. `IFlag
//!   == 0` also builds the 3D curve as that iso (`cxx:2004-2038`),
//! - `|UV1.Y - UV2.Y| <= tol2d`: the same on a V-iso (`cxx:2043-2082`),
//! - `IFlag == 0` otherwise (`cxx:2084-2139`): a Bezier p-curve through the two
//!   UV points, tangent to the two arcs when both common points are on one,
//!   then `GeomLib::BuildCurve3d` lifts it to 3D,
//! - `IFlag > 0` otherwise (`cxx:2141-2151`): `ChFi3d_ProjectPCurv` + a
//!   translation so that the p-curve starts on `UV1`.
//!
//! Translation notes:
//! - `C3d`, `Pcurv`, `Pardeb`, `Parfin` and `tolreached` are OCCT out
//!   parameters (the last two are also *in* for `IFlag > 0`, where `cxx:2143`
//!   loads them into `hc`); the port takes them as `&mut`,
//! - `Geom_TrimmedCurve::BasisCurve()` + `ElCLib::AdjustPeriodic`
//!   (`cxx:2017-2030`, `:2055-2070`) use the port's
//!   [`Curve::untrimmed_basis`] / [`occt_core::elib::clib::adjust_periodic`];
//!   `C3d->Reverse()` is `Arc::from(curve.reversed())`, because the port keeps
//!   curves behind `Arc<dyn Curve>` (mutation goes through a fresh object),
//! - `Adaptor3d_CurveOnSurface Cs(AHC, AHS)` (`cxx:2131-2137`) is
//!   [`CurveOnSurface`], and `GeomLib::BuildCurve3d` is
//!   [`geom_lib_build_curve3d`] with the OCCT defaults `C1` / 14 / 16.

use std::sync::Arc;

use occt_core::bnd::BndBox2d;
use occt_core::elib::clib;
use occt_core::gp::{GpDir2d, GpPnt2d, GpVec2d};
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::PCONFUSION;
use occt_geom::curve::Curve;
use occt_geom::curve_on_surface::CurveOnSurface;
use occt_geom::geom_lib::geom_lib_build_curve3d;
use occt_geom::surface::Surface;
use occt_geom2d::bezier_curve::Geom2dBezierCurve;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;

use crate::chfi_ds::CommonPoint;
use crate::chfi3d::builder_0_helpers::compute_pcurv::ch_fi3d_compute_pcurv;
use crate::chfi3d::builder_0_helpers::mkbound::build_p_curve_vec3d;
use crate::chfi3d::builder_0_helpers::project_pcurv::ch_fi3d_project_pcurv;
use crate::shape::Edge;
use crate::topopebrep::tool::adaptor::{CurveAdaptor, GeomAdaptorCurve, SurfaceAdaptor};

/// The out-parameters of [`ch_fi3d_compute_arete`].
#[derive(Default, Clone)]
pub struct ComputeArete {
    /// `C3d`.
    pub c3d: Option<Arc<dyn Curve>>,
    /// `Pcurv`.
    pub p_curve: Option<Arc<dyn Curve2d>>,
    /// `Pardeb`.
    pub pardeb: f64,
    /// `Parfin`.
    pub parfin: f64,
    /// `tolreached`.
    pub tolreached: f64,
}

/// `ChFi3d_ComputeArete(P1, UV1, P2, UV2, Surf, C3d, Pcurv, Pardeb, Parfin,
/// tol3d, tol2d, tolreached, IFlag)` (`cxx:1984-2152`).
///
/// `inout` carries the incoming `C3d` / `Pcurv` / `Pardeb` / `Parfin` (OCCT
/// only reads them on the `IFlag > 0` arm) and receives the results.
#[allow(clippy::too_many_arguments)]
pub fn ch_fi3d_compute_arete(
    p1: &CommonPoint,
    uv1: &GpPnt2d,
    p2: &CommonPoint,
    uv2: &GpPnt2d,
    surf: &Arc<dyn Surface>,
    inout: &mut ComputeArete,
    tol3d: f64,
    tol2d: f64,
    iflag: i32,
) -> Result<(), String> {
    // `tolreached = tol3d;` (`cxx:2004`).
    inout.tolreached = tol3d;

    if (uv1.x() - uv2.x()).abs() <= tol2d {
        // `Pardeb = UV1.Y(); Parfin = UV2.Y(); C3d = Surf->UIso(UV1.X());`
        // (`cxx:2008-2011`).
        if iflag == 0 {
            inout.pardeb = uv1.y();
            inout.parfin = uv2.y();
            let mut c = surf
                .u_iso_curve(uv1.x())
                .ok_or_else(|| "ChFi3d_ComputeArete: Geom_Surface::UIso is null".to_string())?;
            if inout.pardeb > inout.parfin {
                inout.pardeb = c.reversed_parameter(inout.pardeb);
                inout.parfin = c.reversed_parameter(inout.parfin);
                c = Arc::from(c.reversed());
            }
            // `tc = down_cast<Geom_TrimmedCurve>(C3d); if (!tc.IsNull()) {
            // C3d = tc->BasisCurve(); if (C3d->IsPeriodic())
            // ElCLib::AdjustPeriodic(...); }` (`cxx:2018-2031`).
            if let Some((basis, _, _)) = c.untrimmed_basis() {
                c = basis;
                if c.is_periodic() {
                    clib::adjust_periodic(
                        c.first_parameter(),
                        c.last_parameter(),
                        tol2d,
                        &mut inout.pardeb,
                        &mut inout.parfin,
                    );
                }
            }
            inout.c3d = Some(c);
        }
        if iflag != 1 {
            // `hs->Load(Surf); hc->Load(C3d, Pardeb, Parfin);
            // ChFi3d_ComputePCurv(aHCurve, UV1, UV2, Pcurv, hs, Pardeb,
            // Parfin, tol3d, tolreached, false);` (`cxx:2032-2038`).
            compute_p_curve(inout, uv1, uv2, surf, tol3d)?;
        } else {
            // `Pcurv = new Geom2d_Line(UV1, gp_Vec2d(UV1, UV2));`
            // (`cxx:2041-2042`).
            inout.p_curve = Some(line_through(uv1, uv2)?);
        }
        return Ok(());
    }

    if (uv1.y() - uv2.y()).abs() <= tol2d {
        // `Pardeb = UV1.X(); Parfin = UV2.X(); C3d = Surf->VIso(UV1.Y());`
        // (`cxx:2048-2051`).
        if iflag == 0 {
            inout.pardeb = uv1.x();
            inout.parfin = uv2.x();
            let mut c = surf
                .v_iso_curve(uv1.y())
                .ok_or_else(|| "ChFi3d_ComputeArete: Geom_Surface::VIso is null".to_string())?;
            if inout.pardeb > inout.parfin {
                inout.pardeb = c.reversed_parameter(inout.pardeb);
                inout.parfin = c.reversed_parameter(inout.parfin);
                c = Arc::from(c.reversed());
            }
            // `cxx:2057-2070`.
            if let Some((basis, _, _)) = c.untrimmed_basis() {
                c = basis;
                if c.is_periodic() {
                    clib::adjust_periodic(
                        c.first_parameter(),
                        c.last_parameter(),
                        tol2d,
                        &mut inout.pardeb,
                        &mut inout.parfin,
                    );
                }
            }
            inout.c3d = Some(c);
        }
        if iflag != 1 {
            // `cxx:2071-2078`.
            compute_p_curve(inout, uv1, uv2, surf, tol3d)?;
        } else {
            // `cxx:2079-2082`.
            inout.p_curve = Some(line_through(uv1, uv2)?);
        }
        return Ok(());
    }

    if iflag == 0 {
        // `if (P1.IsVertex() || P2.IsVertex() || !P1.IsOnArc() ||
        // !P2.IsOnArc())` (`cxx:2088-2096`): a straight Bezier.
        let straight = || -> Result<Arc<dyn Curve2d>, String> {
            Ok(Arc::new(
                Geom2dBezierCurve::new(vec![*uv1, *uv2]).map_err(|e| e.to_string())?,
            ))
        };
        let p_curve: Arc<dyn Curve2d> = if p1.is_vertex()
            || p2.is_vertex()
            || !p1.is_on_arc()
            || !p2.is_on_arc()
        {
            straight()?
        } else {
            // `BRepAdaptor_Curve C1(P1.Arc()); C1.D1(P1.ParameterOnArc(), Pp,
            // Vv1); C1.Initialize(P2.Arc()); C1.D1(P2.ParameterOnArc(), Pp,
            // Vv2);` (`cxx:2098-2105`).
            let (_, vv1) = tangent_of_arc(p1)?;
            let (_, vv2) = tangent_of_arc(p2)?;
            // `hs->Load(Surf); Pcurv = ChFi3d_BuildPCurve(hs, UV1, Vv1, UV2,
            // Vv2, true);` (`cxx:2106-2107`).
            let surf_adaptor = SurfaceAdaptor::of_surface(Arc::clone(surf));
            let pc = build_p_curve_vec3d(&surf_adaptor, uv1, &vv1, uv2, &vv2, true)?;
            // `Bnd_Box2d bs; double umin, umax, vmin, vmax;
            // Surf->Bounds(...); bs.Update(...); bs.SetGap(PConfusion());`
            // (`cxx:2112-2117`).
            let (umin, umax) = surf.u_range();
            let (vmin, vmax) = surf.v_range();
            let mut bs = BndBox2d::new();
            bs.update(umin, vmin, umax, vmax);
            bs.set_gap(PCONFUSION);
            // `for (ii = 1; ii <= 4 && aIN; ii++) if
            // (bs.IsOut(down_cast<Geom2d_BezierCurve>(Pcurv)->Pole(ii)))
            // { aIN = false; Pcurv = new Geom2d_BezierCurve(UV1, UV2); }`
            // (`cxx:2118-2128`).
            let mut a_in = true;
            for pole in pc.poles.iter() {
                if !a_in {
                    break;
                }
                if bs.is_out(pole) {
                    a_in = false;
                }
            }
            if a_in {
                Arc::new(pc)
            } else {
                straight()?
            }
        };
        // `Geom2dAdaptor_Curve AC(Pcurv); Adaptor3d_CurveOnSurface Cs(AHC,
        // AHS); Pardeb = Cs.FirstParameter(); Parfin = Cs.LastParameter();`
        // (`cxx:2130-2137`).
        let cons = CurveOnSurface::new(Arc::clone(&p_curve), Arc::clone(surf));
        inout.pardeb = cons.first_parameter();
        inout.parfin = cons.last_parameter();
        // `GeomLib::BuildCurve3d(tol3d, Cs, Pardeb, Parfin, C3d, tolreached,
        // avtol);` (`cxx:2138-2139`).
        let built = geom_lib_build_curve3d(
            tol3d,
            &cons,
            surf,
            &p_curve,
            inout.pardeb,
            inout.parfin,
            Shape::C1,
            14,
            16,
        );
        inout.c3d = built.curve;
        inout.tolreached = built.max_deviation;
        inout.p_curve = Some(p_curve);
        return Ok(());
    }

    // `hs->Load(Surf); hc->Load(C3d, Pardeb, Parfin);
    // ChFi3d_ProjectPCurv(hc, hs, Pcurv, tol3d, tolreached);` (`cxx:2141-2144`).
    let c3d = inout
        .c3d
        .clone()
        .ok_or_else(|| "ChFi3d_ComputeArete: IFlag > 0 without a C3d".to_string())?;
    let adaptor: Arc<dyn Curve> =
        Arc::new(GeomAdaptorCurve::with_range(c3d, inout.pardeb, inout.parfin));
    if let Some(tolreached) = ch_fi3d_project_pcurv(&adaptor, surf, &mut inout.p_curve, tol3d)? {
        inout.tolreached = tolreached;
    }
    // `gp_Pnt2d p2d = Pcurv->Value(Pardeb); if (!UV1.IsEqual(p2d,
    // PConfusion())) { gp_Vec2d v2d(p2d, UV1); Pcurv->Translate(v2d); }`
    // (`cxx:2145-2151`).
    let pc = inout
        .p_curve
        .as_ref()
        .ok_or_else(|| "ChFi3d_ComputeArete: IFlag > 0 without a p-curve".to_string())?;
    let p2d = pc.value(inout.pardeb);
    if !uv1.is_equal(&p2d, PCONFUSION) {
        // `gp_Vec2d v2d(p2d, UV1); Pcurv->Translate(v2d);`
        let mut moved = pc.as_ref().clone_dyn();
        let mut t = occt_core::gp::GpTrsf2d::identity();
        t.set_translation_vec(&GpVec2d::new(uv1.x() - p2d.x(), uv1.y() - p2d.y()));
        moved.transform(&t);
        inout.p_curve = Some(Arc::from(moved));
    }
    Ok(())
}

/// The `hs->Load(Surf); hc->Load(C3d, Pardeb, Parfin); ChFi3d_ComputePCurv(...)`
/// pair of the two iso arms (`cxx:2032-2038`, `:2071-2078`).
fn compute_p_curve(
    inout: &mut ComputeArete,
    uv1: &GpPnt2d,
    uv2: &GpPnt2d,
    surf: &Arc<dyn Surface>,
    tol3d: f64,
) -> Result<(), String> {
    let c3d = inout
        .c3d
        .clone()
        .ok_or_else(|| "ChFi3d_ComputeArete: no C3d for ChFi3d_ComputePCurv".to_string())?;
    let adaptor = GeomAdaptorCurve::with_range(c3d, inout.pardeb, inout.parfin);
    let tolreached = ch_fi3d_compute_pcurv(
        &adaptor,
        uv1,
        uv2,
        &mut inout.p_curve,
        surf.as_ref(),
        inout.pardeb,
        inout.parfin,
        tol3d,
        false,
    )?;
    inout.tolreached = tolreached;
    Ok(())
}

/// `BRepAdaptor_Curve C(arc); C.D1(par, P, D1)` on the arc of a common point
/// (`cxx:2098-2105`).
fn tangent_of_arc(p: &CommonPoint) -> Result<(occt_core::gp::GpPnt, occt_core::gp::GpVec), String> {
    let arc = p
        .arc()
        .ok_or_else(|| "ChFi3d_ComputeArete: a common point without an arc".to_string())?;
    let par = p
        .parameter_on_arc()
        .ok_or_else(|| "ChFi3d_ComputeArete: a common point without an arc parameter".to_string())?;
    CurveAdaptor::of_edge(&Edge(arc.clone()))
        .d1(par)
        .ok_or_else(|| "ChFi3d_ComputeArete: the arc has no curve".to_string())
}

/// `new Geom2d_Line(UV1, gp_Vec2d(UV1, UV2))` (`cxx:2041`, `:2080`).
fn line_through(uv1: &GpPnt2d, uv2: &GpPnt2d) -> Result<Arc<dyn Curve2d>, String> {
    let dir = GpDir2d::new(uv2.x() - uv1.x(), uv2.y() - uv1.y())
        .map_err(|_| "ChFi3d_ComputeArete: gp_Dir2d of a null vector".to_string())?;
    Ok(Arc::new(Geom2dLine::new(occt_core::gp::GpAx2d::new(
        *uv1, dir,
    ))))
}
