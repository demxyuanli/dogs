//! The static helpers of `BRepFill_TrimEdgeTool`.
//!
//! Source: `ModelingAlgorithms/TKBool/BRepFill/BRepFill_TrimEdgeTool.cxx:104-269`.
//!
//! * `Bubble` (`cxx:104-126`) — order the sequence of points by increasing `X`
//!   (increasing bisector parameter).
//! * `EvalParameters` (`cxx:128-197`) — intersect an edge p-curve with the
//!   bisector at `Precision::Confusion()`.
//! * `EvalParametersBis` (`cxx:201-267`) — the same intersection with a caller
//!   supplied tolerance, used by the `SeanceDeRattrapage` retry loop.
//!
//! Unlike `BRepFill_TrimSurfaceTool::EvalParameters`
//! (`BRepFill_TrimSurfaceTool.cxx:236-286`) neither of these removes the
//! duplicates: `IntersectWith` pairs the two sequences itself (`cxx:490-537`).

use occt_core::gp::GpPnt;
use occt_core::precision::CONFUSION;
use occt_geom2d::adaptor_curve::Geom2dAdaptorCurve;
use occt_geom2d::geom2d_int::Geom2dIntGInter;

use super::Error;

/// `Bubble(Seq)` (`cxx:106-126`): order by increasing `X`, exchanging
/// neighbours until no swap happens. The `for` bound is the **original** length
/// (`cxx:111`), as in OCCT.
pub fn bubble(seq: &mut [GpPnt]) {
    let nb_points = seq.len();
    let mut invert = true;
    while invert {
        invert = false;
        for i in 0..nb_points.saturating_sub(1) {
            let p1 = seq[i];
            let p2 = seq[i + 1];
            if p2.x() < p1.x() {
                seq.swap(i, i + 1);
                invert = true;
            }
        }
    }
}

/// `EvalParameters(Bis, AC, Params)` (`cxx:128-197`): intersect `AC` with
/// `Bis` at `Precision::Confusion()` and append
/// `gp_Pnt(bisector parameter, edge parameter, 0)` for every intersection.
pub fn eval_parameters(
    bis: &Geom2dAdaptorCurve,
    ac: &Geom2dAdaptorCurve,
    params: &mut Vec<GpPnt>,
) -> Result<(), Error> {
    // `constexpr double Tol = Precision::Confusion();` (`cxx:132`).
    let tol = CONFUSION;
    append_intersections(bis, ac, params, tol, tol)
}

/// `EvalParametersBis(Bis, AC, Params, Tol)` (`cxx:201-267`): the same
/// intersection with `TolConf == Tol == Tol` (`cxx:207`, `cxx:212`).
pub fn eval_parameters_bis(
    bis: &Geom2dAdaptorCurve,
    ac: &Geom2dAdaptorCurve,
    params: &mut Vec<GpPnt>,
    tol: f64,
) -> Result<(), Error> {
    append_intersections(bis, ac, params, tol, tol)
}

/// `cxx:140-195` / `cxx:212-265`: `Intersector = Geom2dInt_GInter(CAC, CBis,
/// TolConf, Tol)` followed by the point loop, the segment loop and `Bubble`.
fn append_intersections(
    bis: &Geom2dAdaptorCurve,
    ac: &Geom2dAdaptorCurve,
    params: &mut Vec<GpPnt>,
    tol_conf: f64,
    tol: f64,
) -> Result<(), Error> {
    let mut inter = Geom2dIntGInter::new();
    perform(&mut inter, ac, bis, tol_conf, tol);

    // `cxx:147-150` / `cxx:218-221`.
    if !inter.is_done() {
        return Err(Error::IntersectionNotDone);
    }

    // `cxx:152-163` / `cxx:223-234`: one point per intersection; `X` is the
    // parameter on the **second** adaptor (the bisector).
    let nb_points = inter.nb_points();
    for i in 1..=nb_points {
        let u1 = inter.point(i).param_on_second();
        let u2 = inter.point(i).param_on_first();
        params.push(GpPnt::new(u1, u2, 0.0));
    }

    // `cxx:165-195` / `cxx:236-265`: one point per overlapping segment.
    let nb_segments = inter.nb_segments();
    for i in 1..=nb_segments {
        let seg = inter.segment(i);
        let u1 = seg.first_point().param_on_second();
        let u_last = seg.last_point().param_on_second();
        if (u1 - bis.first_parameter()).abs() <= tol
            && (u_last - bis.last_parameter()).abs() <= tol
        {
            // `cxx:173-181`: the whole bisector is covered — both endpoints.
            params.push(GpPnt::new(u1, seg.first_point().param_on_first(), 0.0));
            params.push(GpPnt::new(u_last, seg.last_point().param_on_first(), 0.0));
        } else {
            // `cxx:182-192`: the middle of the segment.
            let mut u1 = u1;
            u1 += seg.last_point().param_on_second();
            u1 /= 2.0;
            let mut u2 = seg.first_point().param_on_first();
            u2 += seg.last_point().param_on_first();
            u2 /= 2.0;
            params.push(GpPnt::new(u1, u2, 0.0));
        }
    }

    // `cxx:196-197` / `cxx:266-267`.
    bubble(params);
    Ok(())
}

/// `Geom2dInt_GInter(A1, A2, TolConf, Tol)` (`Geom2dInt_GInter_0.cxx:29-89`):
/// both domains are `ComputeDomain(adaptor, max(TolConf, Tol))` over the
/// adaptor's own `myFirst` / `myLast`, which `load` leaves untouched
/// (`Geom2dAdaptor_Curve.cxx:250-258`).
fn perform(
    inter: &mut Geom2dIntGInter,
    ac: &Geom2dAdaptorCurve,
    bis: &Geom2dAdaptorCurve,
    tol_conf: f64,
    tol: f64,
) {
    let tol_domain = tol.max(tol_conf);
    let c1 = ac.curve().expect("loaded edge p-curve");
    let c2 = bis.curve().expect("loaded bisector");
    let d1 = Geom2dIntGInter::domain_over_range(
        c1.as_ref(),
        ac.first_parameter(),
        ac.last_parameter(),
        tol_domain,
    );
    let d2 = Geom2dIntGInter::domain_over_range(
        c2.as_ref(),
        bis.first_parameter(),
        bis.last_parameter(),
        tol_domain,
    );
    inter.perform(c1.as_ref(), &d1, c2.as_ref(), &d2, tol_conf, tol);
}
