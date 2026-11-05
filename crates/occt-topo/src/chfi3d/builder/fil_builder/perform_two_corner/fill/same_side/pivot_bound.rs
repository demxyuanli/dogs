//! The pivot boundary and the pivot p-curve of the `!sameparam` filling
//! (`ChFi3d_FilBuilder_C2.cxx:605-665`).
//!
//! When the two common points are not at the same parameter on the pivot arc,
//! the filling takes a fourth boundary, `GeomFill_SimpleBound` over
//! `Hpivot->Trim(min(parCP1,parCP2), max(parCP1,parCP2), tolesp)`
//! (`cxx:604-607`), and `CompleteData` takes the p-curve of the pivot arc on
//! one of the two incident faces (`cxx:608-665`).

use std::sync::Arc;

use occt_core::bspl::knots::reparametrize;
use occt_geom::geom_fill::boundary::{Boundary, SimpleBound};
use occt_geom::{Curve, GeomTrimmedCurve};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::geom2d_convert::{curve_to_bspline_curve_bspl, DEFAULT_PARAMETERISATION};
use occt_geom2d::trimmed::Geom2dTrimmedCurve;
use occt_geom2d::Geom2dBSplineCurve;

use crate::abs::ShapeType;
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::chfi3d::builder::ChFi3dBuilder;
use crate::shape::Edge;
use crate::topo_tools::{explorer_shapes, is_same};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::super::CornerCtx;

/// What `cxx:604-665` leaves for `fil.Init(...)` and `CompleteData(...)`:
/// `Bpiv`, `HFaPiv` and `PCurveOnPiv`.
pub(super) struct PivotBound {
    /// `Bpiv`
    pub boundary: Arc<dyn Boundary>,
    /// `HFaPiv`
    pub surface: SurfaceAdaptor,
    /// `PCurveOnPiv`
    pub p_curve: Option<Arc<dyn Curve2d>>,
}

/// `cxx:604-665`, the `else` arm of `if (sameparam)`.
///
/// `b1` is the `B1` boundary the caller already built: the body reads
/// `B1->Bounds` / `B1->D1` to decide which face carries the pivot p-curve.
/// `ctx.fa_piv`, `ctx.pctrans`, `ctx.reset_cp1` and `ctx.reset_cp2` are the
/// out-parameters of the block.
pub(super) fn pivot_bound(
    base: &ChFi3dBuilder,
    ctx: &mut CornerCtx,
    b1: &Arc<dyn Boundary>,
) -> Result<PivotBound, String> {
    let hpivot = ctx
        .hpivot
        .as_ref()
        .ok_or_else(|| "PerformTwoCorner: null Hpivot".to_string())?;
    let min_par = ctx.par_cp1.min(ctx.par_cp2);
    let max_par = ctx.par_cp1.max(ctx.par_cp2);

    // `occ::handle<Adaptor3d_Curve> HPivTrim = Hpivot->Trim(min(parCP1,
    // parCP2), max(parCP1, parCP2), tolesp); Bpiv = new
    // GeomFill_SimpleBound(HPivTrim, tolapp3d, 2.e-4);` (`cxx:605-607`).
    //
    // `GeomAdaptor_Curve::Trim` answers the adaptor itself when the window is
    // its own range and a `Geom_TrimmedCurve` over `Curve()` otherwise, so the
    // untrimmed arm keeps the adaptor's curve as it stands.
    let basis = hpivot
        .curve()
        .cloned()
        .ok_or_else(|| "PerformTwoCorner: null Hpivot curve".to_string())?;
    let trimmed: Arc<dyn Curve> =
        if min_par == hpivot.first_parameter() && max_par == hpivot.last_parameter() {
            basis
        } else {
            Arc::new(GeomTrimmedCurve::new(basis, min_par, max_par))
        };
    let boundary: Arc<dyn Boundary> = Arc::new(SimpleBound::new(trimmed, base.tolapp3d, 2.0e-4));

    // `Hpivot->D1(parCP1, bidon, dArc); B1->Bounds(fb1, lb1); B1->D1(lb1,
    // bidon, dcf); bool pivotverslebas = dArc.Dot(dcf) <= 0.;`
    // (`cxx:610-617`).
    let (_, d_arc) = hpivot
        .d1(ctx.par_cp1)
        .ok_or_else(|| "PerformTwoCorner: null Hpivot value".to_string())?;
    let (_, lb1) = b1.bounds();
    let (_, dcf) = b1.d1(lb1);
    let pivot_vers_le_bas = d_arc.dot(&dcf) <= 0.0;
    // `bool pcfalenvers = (parCP1 > parCP2);` (`cxx:618`).
    let pcfalenvers = ctx.par_cp1 > ctx.par_cp2;

    // `if ((pivotverslebas && !pcfalenvers) || (!pivotverslebas &&
    // pcfalenvers)) { FaPiv = FF2; HFaPiv = HBRS2; resetcp2 = true; } else {
    // FaPiv = FF1; HFaPiv = HBRS1; resetcp1 = true; }` (`cxx:619-628`).
    let (fa_piv, surface) = if (pivot_vers_le_bas && !pcfalenvers)
        || (!pivot_vers_le_bas && pcfalenvers)
    {
        ctx.reset_cp2 = true;
        (ctx.ff2.clone(), SurfaceAdaptor::of_face(&ctx.ff2))
    } else {
        ctx.reset_cp1 = true;
        (ctx.ff1.clone(), SurfaceAdaptor::of_face(&ctx.ff1))
    };
    // `FaPiv.Orientation(TopAbs_FORWARD);` (`cxx:629`).
    ctx.fa_piv = Some(crate::shape::Face(fa_piv.0.oriented(crate::abs::Orientation::Forward)));
    let fa_piv = ctx
        .fa_piv
        .clone()
        .ok_or_else(|| "PerformTwoCorner: null FaPiv".to_string())?;
    let pivot = Edge(ctx.pivot.clone());

    // `pcpivot.Initialize(pivot, FaPiv);` (`cxx:630`).
    let pcpivot = Curve2dAdaptor::of_edge_face(&pivot, &fa_piv);

    // `TopExp_Explorer Expl; for (Expl.Init(FaPiv, TopAbs_EDGE); ...) if
    // (Expl.Current().IsSame(pivot)) { pctrans = Expl.Current().Orientation();
    // break; }` (`cxx:632-641`).
    for e in explorer_shapes(&fa_piv.0, ShapeType::Edge) {
        if is_same(&e, &ctx.pivot) {
            ctx.pctrans = e.orientation();
            break;
        }
    }

    let pcurve = pcpivot
        .curve()
        .cloned()
        .ok_or_else(|| "PerformTwoCorner: null BRepAdaptor_Curve2d curve".to_string())?;

    // `if (pcpivot.GetType() != GeomAbs_BSplineCurve) { PCurveOnPiv =
    // CurveToBSplineCurve(new Geom2d_TrimmedCurve(pcpivot.Curve(), min, max));
    // } else { PCurveOnPiv = SplitBSplineCurve(down_cast<Geom2d_BSplineCurve>
    // (pcpivot.Curve()), min, max, tol2d); }` (`cxx:642-656`).
    //
    // `bspline2d()` is the `down_cast` of `cxx:652`, so a p-curve of another
    // type takes the `CurveToBSplineCurve` arm exactly when `GetType()` is not
    // `GeomAbs_BSplineCurve`.
    let mut bs: Geom2dBSplineCurve = match pcurve.bspline2d() {
        None => {
            let trc = Geom2dTrimmedCurve::new(Arc::clone(&pcurve), min_par, max_par);
            curve_to_bspline_curve_bspl(&trc, DEFAULT_PARAMETERISATION)
                .ok_or_else(|| "PerformTwoCorner: Geom2dConvert::CurveToBSplineCurve".to_string())?
        }
        Some(exact) => {
            // `Geom2dConvert::SplitBSplineCurve(C, U1, U2, Tol)` copies `C`,
            // runs `Segment(U1, U2)` and keeps `SameOrientation = true` (its
            // default), which is `Geom2dBSplineCurve::segment`.
            let mut copy = exact.clone();
            copy.segment(min_par, max_par, base.tol2d)
                .map_err(|e| format!("PerformTwoCorner: SplitBSplineCurve: {e}"))?;
            copy
        }
    };

    // `NCollection_Array1<double> kk(PCurveOnPiv->Knots());
    // BSplCLib::Reparametrize(0., 1., kk); PCurveOnPiv->SetKnots(kk);`
    // (`cxx:657-659`).
    let (mut kk, _mults) = bs.distinct_knots_and_mults();
    reparametrize(0.0, 1.0, &mut kk);
    bs.set_knots(&kk)
        .map_err(|e| format!("PerformTwoCorner: Geom2d_BSplineCurve::SetKnots: {e}"))?;

    // `if (pcfalenvers) { PCurveOnPiv->Reverse(); pctrans =
    // TopAbs::Reverse(pctrans); }` (`cxx:660-664`).
    if pcfalenvers {
        bs.reverse();
        ctx.pctrans = ctx.pctrans.reversed();
    }

    Ok(PivotBound {
        boundary,
        surface,
        p_curve: Some(Arc::new(bs) as Arc<dyn Curve2d>),
    })
}
