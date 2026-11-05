//! `ChFi3d_ComputesIntPC` - refine the intersection of two p-curves of two
//! fillet surfaces.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:4382-4447`
//! (OCCT 8.0.0): the six-argument body (`cxx:4389-4398`) that discards the
//! point, and the seven-argument one it forwards to (`cxx:4402-4447`).
//!
//! The two incoming parameters are only the seed of the search: each p-curve is
//! re-adapted over `[U - 0.1, U + 0.1]` (clamped to `0.05` of the p-curve's own
//! span), and `Extrema_LocateExtCC` moves them onto the local extremum. The
//! result only replaces the seeds when it is closer than they were.
//!
//! Translation notes:
//! - `Geom2dAdaptor_Curve(PCurve, U - delt, U + delt)` is
//!   [`occt_geom2d::adaptor_trim::TrimmedAdaptor2d`], and
//!   `Adaptor3d_CurveOnSurface(hc2d, HS)` is
//!   [`occt_geom::curve_on_surface::CurveOnSurface`],
//! - `Extrema_LocateExtCC` is [`occt_geom::extrema_cc::locate_extcc`]; its
//!   `None` is `!ext.IsDone()`,
//! - `P.SetXYZ(0.5 * (p3d1.XYZ() + p3d2.XYZ()))` is
//!   [`GpPnt::from_xyz`] over the summed [`GpXyz`],
//! - every null p-curve or surface OCCT dereferences is the `Err` of this port.

use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_geom::curve_on_surface::CurveOnSurface;
use occt_geom::extrema_cc::locate_extcc;
use occt_geom2d::adaptor_trim::TrimmedAdaptor2d;

use crate::chfi_ds::FaceInterference;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `P.SetXYZ(0.5 * (p3d1.XYZ() + p3d2.XYZ()))`.
fn midpoint(p3d1: &GpPnt, p3d2: &GpPnt) -> GpPnt {
    GpPnt::from_xyz(&p3d1.xyz().added(p3d2.xyz()).multiplied(0.5))
}

/// `ChFi3d_ComputesIntPC(Fi1, Fi2, HS1, HS2, UInt1, UInt2)` (`cxx:4389-4398`).
///
/// The point the seven-argument body computes is discarded by this overload,
/// exactly as in OCCT.
pub fn ch_fi3d_computes_int_pc(
    fi1: &FaceInterference,
    fi2: &FaceInterference,
    hs1: &SurfaceAdaptor,
    hs2: &SurfaceAdaptor,
    uint1: &mut f64,
    uint2: &mut f64,
) -> Result<(), String> {
    ch_fi3d_computes_int_pc_point(fi1, fi2, hs1, hs2, uint1, uint2).map(|_| ())
}

/// `ChFi3d_ComputesIntPC(Fi1, Fi2, HS1, HS2, UInt1, UInt2, P)` (`cxx:4402-4447`).
///
/// Answers the point `P` of the seven-argument body; `uint1` / `uint2` are
/// in/out like the OCCT references.
pub fn ch_fi3d_computes_int_pc_point(
    fi1: &FaceInterference,
    fi2: &FaceInterference,
    hs1: &SurfaceAdaptor,
    hs2: &SurfaceAdaptor,
    uint1: &mut f64,
    uint2: &mut f64,
) -> Result<GpPnt, String> {
    // `double x, y, distref2; Fi1.PCurveOnSurf()->Value(UInt1).Coord(x, y);
    // gp_Pnt p3d1 = HS1->Value(x, y);` (`cxx:4414-4416`).
    let pc1 = fi1
        .p_curve_on_surf()
        .cloned()
        .ok_or_else(|| "ChFi3d_ComputesIntPC: Fi1 has no p-curve on surf".to_string())?;
    let pc2 = fi2
        .p_curve_on_surf()
        .cloned()
        .ok_or_else(|| "ChFi3d_ComputesIntPC: Fi2 has no p-curve on surf".to_string())?;
    let s1 = hs1
        .surface()
        .cloned()
        .ok_or_else(|| "ChFi3d_ComputesIntPC: HS1 has no surface".to_string())?;
    let s2 = hs2
        .surface()
        .cloned()
        .ok_or_else(|| "ChFi3d_ComputesIntPC: HS2 has no surface".to_string())?;

    let uv = pc1.value(*uint1);
    let p3d1 = s1.d0(uv.x(), uv.y());
    let uv = pc2.value(*uint2);
    let p3d2 = s2.d0(uv.x(), uv.y());
    // `distref2 = p3d1.SquareDistance(p3d2); P.SetXYZ(0.5 * (...));`
    // (`cxx:4419-4420`).
    let distref2 = p3d1.square_distance(&p3d2);
    let mut p = midpoint(&p3d1, &p3d2);

    // `double delt1 = std::min(0.1, 0.05 * (Fi1.LastParameter() -
    // Fi1.FirstParameter()));` (`cxx:4422`).
    let delt1 = 0.1f64.min(0.05 * (fi1.last_parameter() - fi1.first_parameter()));
    let delt2 = 0.1f64.min(0.05 * (fi2.last_parameter() - fi2.first_parameter()));

    // `Geom2dAdaptor_Curve hc2d1(Fi1.PCurveOnSurf(), UInt1 - delt1, UInt1 +
    // delt1); Adaptor3d_CurveOnSurface cons1(hc2d1, HS1);` (`cxx:4424-4427`).
    let cons1 = CurveOnSurface::new(
        Arc::new(TrimmedAdaptor2d::new(
            Arc::clone(&pc1),
            *uint1 - delt1,
            *uint1 + delt1,
        )) as Arc<dyn occt_geom2d::curve::Curve2d>,
        s1,
    );
    let cons2 = CurveOnSurface::new(
        Arc::new(TrimmedAdaptor2d::new(
            Arc::clone(&pc2),
            *uint2 - delt2,
            *uint2 + delt2,
        )) as Arc<dyn occt_geom2d::curve::Curve2d>,
        s2,
    );

    // `Extrema_LocateExtCC ext(cons1, cons2, UInt1, UInt2); if (ext.IsDone())`
    // (`cxx:4430-4432`).
    let Some(ext) = locate_extcc(&cons1, &cons2, *uint1, *uint2) else {
        return Ok(p);
    };
    // `double dist2 = ext.SquareDistance(); if (dist2 < distref2)`
    // (`cxx:4434-4436`).
    let dist2 = ext.distance * ext.distance;
    if dist2 < distref2 {
        // `ext.Point(ponc1, ponc2); UInt1 = ponc1.Parameter(); UInt2 =
        // ponc2.Parameter();` (`cxx:4437-4439`).
        *uint1 = ext.u1;
        *uint2 = ext.u2;
        // `gp_Pnt Pnt1 = ponc1.Value(); gp_Pnt Pnt2 = ponc2.Value();
        // P.SetXYZ(0.5 * (Pnt1.XYZ() + Pnt2.XYZ()));` (`cxx:4440-4443`).
        p = midpoint(&ext.p1, &ext.p2);
    }
    Ok(p)
}
