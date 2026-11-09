//! `LocOpe_Pipe::Curves` / `LocOpe_Pipe::BarycCurve`
//! (`LocOpe_Pipe.cxx:329-426`, `:449-549`).
//!
//! Both members walk the wire `BRepFill_Pipe::PipeLine(P1)` sweeps, convert
//! every edge to a B-spline, reparametrise its knot vector onto
//! `[i, i + 1]` and stitch the pieces into one curve whose knots are the
//! distinct knots of the pieces minus the shared ends. The stitching body is
//! identical in the two members (`cxx:344-366` / `cxx:471-491` and
//! `cxx:377-424` / `cxx:503-548`), so it lives here once; the two callers in
//! `mod.rs` keep their own empty-sequence branches, which differ.

use occt_core::bspl::knots::reparametrize;
use occt_core::convert::ParameterisationType;
use occt_core::gp::GpPnt;
use occt_core::precision::PCONFUSION;
use occt_geom::convert_bspl::curve_to_bspline_curve;
use occt_geom::{Curve, GeomBSplineCurve};

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::shape::Edge;
use crate::topo_tools::explorer_shapes;

use super::Pipe;

/// `cxx:341-366` / `cxx:468-491`: the B-spline pieces of the wire swept by
/// `p1`, each reparametrised onto `[i, i + 1]`, together with the maximum
/// degree met (`MaxDeg`).
///
/// The `BRep_Tool::Curve(edge, p1, p2)` pair maps to
/// [`BRepTool::edge_curve`] plus [`BRepTool::edge_parameters`], the reference's
/// local-parameterisation form (`BRep_Tool.cxx:619-650`).
pub(super) fn collect_bsplines(
    pipe: &mut Pipe,
    p1: &GpPnt,
) -> Result<(Vec<GeomBSplineCurve>, usize), String> {
    // `cxx:342-343` / `cxx:469-470`.
    let mut max_deg = 0usize;
    let mut seq: Vec<GeomBSplineCurve> = Vec::new();

    // `cxx:344` / `cxx:471`.
    let wire = pipe.pipe_mut().pipe_line(p1)?;

    // `cxx:347-348` / `cxx:473-474`.
    for sub in explorer_shapes(&wire.0, ShapeType::Edge) {
        let edge = Edge::wrap(sub)
            .ok_or_else(|| "Standard_TypeMismatch: LocOpe_Pipe::Curves".to_string())?;
        // `cxx:349` / `cxx:476`: `BRep_Tool::Curve(edge, p1, p2)`.
        let c1 = BRepTool::edge_curve(&edge)
            .ok_or_else(|| "Standard_NullObject: LocOpe_Pipe::Curves".to_string())?;
        let (p1_edge, p2_edge) = BRepTool::edge_parameters(&edge);
        // `cxx:350` / `cxx:477`: `GeomConvert::CurveToBSplineCurve(C1)` with the
        // default `Convert_TgtThetaOver2` parameterisation
        // (`GeomConvert.hxx:270-272`).
        let mut c = curve_to_bspline_curve(&*c1, ParameterisationType::TgtThetaOver2)
            .map_err(|e| format!("GeomConvert::CurveToBSplineCurve: {e:?}"))?;
        // `cxx:352-355` / `cxx:479-482`: `if (C.IsNull()) continue;`. The
        // reference's null answer is defensive: `CurveToBSplineCurve` raises
        // for every curve it cannot convert, and a missing edge curve is the
        // null `C1` case, which the port answers as an error above.
        // `cxx:356` / `cxx:483`.
        max_deg = max_deg.max(c.degree());
        // `cxx:357`: `P1 = C->Value(p2)` is dead in the reference - the outer
        // loop re-reads `Spt(i)` into `P1` and the edge loop never reads it
        // after this assignment. Kept as the same evaluation.
        let _next_start = c.d0(p2_edge);
        // `cxx:358-361` / `cxx:484-487`.
        if p1_edge != c.first_parameter() || p2_edge != c.last_parameter() {
            c.segment(p1_edge, p2_edge, PCONFUSION)
                .map_err(str::to_string)?;
        }
        // `cxx:362-365` / `cxx:488-491`: `BSplCLib::Reparametrize(seq.Length(),
        // seq.Length() + 1, Tkn)` on the **distinct** knots, then `SetKnots`.
        let mut tkn = c.distinct_knots().to_vec();
        let first = seq.len() as f64;
        reparametrize(first, first + 1.0, &mut tkn);
        c.set_knots(&tkn).map_err(str::to_string)?;
        // `cxx:366` / `cxx:492`.
        seq.push(c);
    }

    Ok((seq, max_deg))
}

/// `cxx:377-424` / `cxx:503-548`: raise every piece to `max_deg`, then join
/// them, dropping the shared end pole and knot of each following piece.
pub(super) fn merge_bsplines(
    seq: &mut [GeomBSplineCurve],
    max_deg: usize,
) -> Result<GeomBSplineCurve, String> {
    // `cxx:368-370` / `cxx:495-497`.
    let nb_curv = seq.len();
    debug_assert!(nb_curv >= 1, "the callers handle the empty sequence");

    // `cxx:376-383` / `cxx:501-508`.
    let mut nbp = 0usize;
    let mut nbkn = 0usize;
    for bsp in seq.iter_mut() {
        bsp.increase_degree(max_deg).map_err(str::to_string)?;
        nbp += bsp.nb_poles();
        // `Geom_BSplineCurve::NbKnots()` is `myKnots.Length()`, i.e. the
        // **distinct** knot count (`Geom_BSplineCurve_1.cxx:380-383`).
        nbkn += bsp.distinct_knots().len();
    }
    // `cxx:386-387` / `cxx:508-509`.
    nbp -= nb_curv - 1;
    nbkn -= nb_curv - 1;

    // `cxx:388-389` / `cxx:510-511`.
    let mut tkn: Vec<f64> = Vec::with_capacity(nbkn);
    let mut tmu: Vec<i32> = Vec::with_capacity(nbkn);
    let mut tpol: Vec<GpPnt> = Vec::with_capacity(nbp);

    // `cxx:391-403` / `cxx:514-526`: the first piece contributes all its poles
    // and knots, and its last multiplicity loses one occurrence.
    {
        let bsp = &seq[0];
        for k in 0..bsp.nb_poles() {
            tpol.push(bsp.pole(k).clone());
        }
        for k in 0..bsp.distinct_knots().len() {
            tkn.push(bsp.distinct_knots()[k]);
            tmu.push(bsp.multiplicity((k + 1) as i32));
        }
        if let Some(last) = tmu.last_mut() {
            *last -= 1;
        }
    }

    // `cxx:405-420` / `cxx:528-545`: every following piece drops its first
    // pole and its first knot.
    for j in 1..nb_curv {
        let bsp = &seq[j];
        for k in 1..bsp.nb_poles() {
            tpol.push(bsp.pole(k).clone());
        }
        let nk = bsp.distinct_knots().len();
        for k in 1..nk {
            tkn.push(bsp.distinct_knots()[k]);
            tmu.push(bsp.multiplicity((k + 1) as i32));
        }
        if let Some(last) = tmu.last_mut() {
            *last -= 1;
        }
    }
    // `cxx:421` / `cxx:546`.
    if let Some(last) = tmu.last_mut() {
        *last += 1;
    }

    // `cxx:422` / `cxx:547`: `new Geom_BSplineCurve(Tpol, Tkn, Tmu, MaxDeg)`.
    GeomBSplineCurve::from_poles_knots_mults(tpol, tkn, tmu, max_deg)
        .map_err(str::to_string)
}
