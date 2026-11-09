//! The two tool classes `Contap_Contour` calls on its 2D restrictions:
//! `Contap_HContTool` and `Contap_HCurve2dTool`.
//!
//! Ported from `Contap_HContTool.cxx:230-305` (`NbSamplesOnArc`, `Project`) and
//! `Contap_HCurve2dTool.cxx:19-56` (`Value`, `D1`, `Resolution`). The rest of
//! `Contap_HContTool` (`cxx:19-180` sample counts, `:307-371` the
//! `HasBeenSeen`/`NbPoints`/`Value`/`IsVertex`/`Vertex`/`NbSegments` stubs) has
//! no caller in the `Contap_Contour` chain: `:307-371` answer a constant or
//! throw `Standard_OutOfRange` in OCCT, and `:19-228` feed
//! `IntStart_SearchInside` only through `TheTopolTool` (`TopolTool::nb_samples_u`
//! and friends) in the ported code.
//!
//! The projection is `Extrema_EPCOfExtPC2d(P, *C, 20, 1e-8, 1e-5)` followed by
//! a scan for the smallest `SquareDistance` (`cxx:283-303`).

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::extrema2d::EPcOfExtPc2d;

use crate::brep_adaptor_curve2d::Curve2dAdaptor;

/// `Contap_HCurve2dTool::Value(C, U)` (`Contap_HCurve2dTool.cxx:19-22`).
pub fn hcurve2d_value(c: &Curve2dAdaptor, u: f64) -> GpPnt2d {
    c.d0(u).unwrap_or_else(GpPnt2d::zero)
}

/// `Contap_HCurve2dTool::D1(C, U, P, T)` (`Contap_HCurve2dTool.cxx:34-42`).
pub fn hcurve2d_d1(c: &Curve2dAdaptor, u: f64) -> (GpPnt2d, GpVec2d) {
    c.d1(u)
        .unwrap_or_else(|| (GpPnt2d::zero(), GpVec2d::new(0.0, 0.0)))
}

/// `Contap_HCurve2dTool::Resolution(C, Ruv)` (`Contap_HCurve2dTool.cxx:50-56`).
pub fn hcurve2d_resolution(c: &Curve2dAdaptor, ruv: f64) -> f64 {
    c.resolution(ruv).unwrap_or(ruv)
}

/// `Contap_HContTool::NbSamplesOnArc(A)` (`Contap_HContTool.cxx:230-259`).
pub fn nb_samples_on_arc(a: &Curve2dAdaptor) -> i32 {
    let Some(curve) = a.curve() else {
        // No p-curve loaded: `GetType()` answers `GeomAbs_OtherCurve` on an
        // empty `Geom2dAdaptor_Curve`, which is the `default` arm.
        return 10;
    };
    let c: &dyn Curve2d = curve.as_ref();
    let nbs_on_c = if c.is_line() {
        2.0
    } else if c.is_geom2d_circle()
        || c.gp_elips2d().is_some()
        || c.gp_hypr2d().is_some()
        || c.gp_parab2d().is_some()
    {
        10.0
    } else if c.is_bezier2d() {
        // `A->NbPoles()`.
        c.bezier_nb_poles().unwrap_or(0) as f64
    } else if c.is_bspline2d() {
        // `2 + A->NbKnots() * A->Degree()`.
        match (c.bspline_distinct_knots_mults(), c.bspline2d()) {
            (Some((_knots, mults)), Some(bs)) => 2.0 + (mults.len() * bs.degree()) as f64,
            _ => 10.0,
        }
    } else {
        10.0
    };
    nbs_on_c as i32
}

/// `Contap_HContTool::Project(C, P, Paramproj, Ptproj)`
/// (`Contap_HContTool.cxx:269-303`).
///
/// `false` where OCCT returns `false`: the extrema search failed, or it found
/// no extremum at all.
pub fn hcont_project(c: &Curve2dAdaptor, p: &GpPnt2d) -> Option<(f64, GpPnt2d)> {
    let curve = c.curve()?;
    let nb_u = 20;
    let eps_x = 1.0e-8;
    let tol = 1.0e-5;

    let mut extrema = EPcOfExtPc2d::new(p, curve.as_ref(), nb_u, eps_x, tol);
    extrema.perform(p);
    if !extrema.is_done() {
        return None;
    }
    let nb_ext = extrema.nb_ext();
    if nb_ext == 0 {
        return None;
    }
    let mut index_min = 1usize;
    let mut dist2 = extrema.square_distance(1);
    for i in 2..=nb_ext {
        if extrema.square_distance(i) < dist2 {
            index_min = i;
            dist2 = extrema.square_distance(i);
        }
    }
    let _ = dist2;
    Some(extrema.point(index_min))
}
