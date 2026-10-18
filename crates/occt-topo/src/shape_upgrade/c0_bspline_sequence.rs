//! `ShapeUpgrade::C0BSplineToSequenceOfC1BSplineCurve`
//! (`ShapeUpgrade.cxx:40-130` for `Geom_BSplineCurve`, `:174-194` for
//! `Geom2d_BSplineCurve`).
//!
//! A C0 B-spline is "unusable in topology algorithms"
//! (`IGESToBRep_TopoCurve.cxx:995-998`), so the IGES reader splits it into
//! C1 pieces: `IGESToBRep_TopoCurve::TransferTopoBasicCurve` calls it at
//! `IGESToBRep_TopoCurve.cxx:1022` and `Transfer2dTopoBasicCurve` at `:1238`,
//! both through the pluggable `ShapeAlgo` container hook
//! (`ShapeAlgo_AlgoContainer.cxx:636-650`).
//!
//! The split is exact: it cuts at the knots whose multiplicity reaches the
//! degree and copies a contiguous pole/weight run per piece, so the geometry
//! does not change (`ShapeUpgrade.hxx:49-53`). OCCT's other splitter,
//! `GeomConvert::C0BSplineToArrayOfC1BSplineCurve`, is the tolerance-driven
//! one and is not the one the IGES reader reaches.
//!
//! The 2D overload (`cxx:174-194`) lifts the curve with `BSplineCurve2dTo3d`
//! (`cxx:133-152`), runs the 3D split, and drops each piece back with
//! `BSplineCurve3dTo2d` (`cxx:154-172`); the port shares one body for both.

use std::sync::Arc;

use occt_core::bspl::banded_interp::knot_sequence;
use occt_core::bspl::knots as bspl_knots;
use occt_core::bspl::locate::{first_u_knot_index, flat_index, last_u_knot_index};
use occt_core::gp::GpPnt;
use occt_core::precision::REAL_SMALL;
use occt_geom::bspline_curve::GeomBSplineCurve;
use occt_geom::convert_bspl::{bspline_copy, is_bspline};
use occt_geom::curve::Curve;
use occt_geom2d::bspline_curve::Geom2dBSplineCurve;
use occt_geom2d::curve::Curve2d;

/// `GeomAbs_C1` (`GeomAbs_Shape.hxx`), the minimum continuity the split leaves
/// behind. `ShapeUpgrade.cxx:43` and `:181` test `IsCN(1)`, which for a
/// B-spline is `Continuity() >= GeomAbs_C1`
/// (`Geom_BSplineCurve_1.cxx:34-60`); `IGESToBRep_TopoCurve.cxx:1016` and
/// `:1232` test the same bound before calling in.
pub const GEOMABS_C1: u8 = 2;

/// One B-spline produced by the split, in the constructor terms of
/// `Geom_BSplineCurve(Poles, Weights, Knots, Mults, Degree)`
/// (`Geom_BSplineCurve.cxx:177-225`): distinct knots plus multiplicities and
/// weights already checked by `Rational` (`cxx:95-105`); `None` weights are
/// OCCT's non-rational result.
struct Piece {
    /// 3D poles; a 2D input is lifted with `z = 0` (`cxx:143-147`).
    poles: Vec<GpPnt>,
    weights: Option<Vec<f64>>,
    knots: Vec<f64>,
    mults: Vec<i32>,
    degree: i32,
}

impl Piece {
    /// The `Geom_BSplineCurve` the C++ code appends to `seqBS`
    /// (`ShapeUpgrade.cxx:122-126`).
    fn to_3d(&self) -> Result<GeomBSplineCurve, &'static str> {
        let flat = knot_sequence(&self.knots, &self.mults, self.degree);
        match &self.weights {
            Some(w) => GeomBSplineCurve::rational(self.poles.clone(), w.clone(), flat, self.degree as usize),
            None => GeomBSplineCurve::new(self.poles.clone(), flat, self.degree as usize),
        }
    }

    /// `BSplineCurve3dTo2d` (`ShapeUpgrade.cxx:154-172`): the 2D pole is the
    /// 3D pole's `(X, Y)`.
    fn to_2d(&self) -> Result<Geom2dBSplineCurve, &'static str> {
        let flat = knot_sequence(&self.knots, &self.mults, self.degree);
        let xs: Vec<f64> = self.poles.iter().map(|p| p.x()).collect();
        let ys: Vec<f64> = self.poles.iter().map(|p| p.y()).collect();
        Geom2dBSplineCurve::from_flat(xs, ys, self.weights.clone(), flat, self.degree as usize, false)
    }
}

/// `Geom_BSplineCurve::IsCN`-style gate plus `Rational`
/// (`ShapeUpgrade.cxx:42-47`, `:179-182`): the split runs only for a B-spline
/// that is not already C1. `continuity()` is the curve's `mySmooth`.
fn is_c0_bspline(continuity: u8) -> bool {
    continuity < GEOMABS_C1
}

/// `ShapeUpgrade::C0BSplineToSequenceOfC1BSplineCurve(BS, seqBS)`
/// (`ShapeUpgrade.cxx:40-130`). `None` is OCCT's `false`: the curve is C1
/// already or is not a `Geom_BSplineCurve`.
pub fn c0_bspline_to_sequence_of_c1_bspline_curve(bs: &dyn Curve) -> Option<Vec<Arc<dyn Curve>>> {
    // `cxx:42-47`: the `IsKind(STANDARD_TYPE(Geom_BSplineCurve))` down_cast
    // plus `BS->IsCN(1)`. A `Geom_TrimmedCurve` forwards `IsCN` to its basis
    // and is excluded by `is_bspline`, exactly as the down_cast excludes it.
    if !is_bspline(bs) || !is_c0_bspline(bs.continuity()) {
        return None;
    }
    // `cxx:49-50,56-58`: work on `BS->Copy()`.
    let mut work = bspline_copy(bs).ok()?;
    // `cxx:50`: to have NbPoles = sum of the multiplicities - degree - 1.
    work.set_not_periodic();
    let (uknots, umults) = work.distinct_knots_and_mults();
    let pieces = split_c0_bspline(
        work.degree as i32,
        &work.knots,
        &work.poles,
        work.weights.as_deref(),
        &uknots,
        &umults,
    )?;
    let mut out = Vec::with_capacity(pieces.len());
    for piece in &pieces {
        out.push(Arc::new(piece.to_3d().ok()?) as Arc<dyn Curve>);
    }
    Some(out)
}

/// `ShapeUpgrade::C0BSplineToSequenceOfC1BSplineCurve(BS, seqBS)` for
/// `Geom2d_BSplineCurve` (`ShapeUpgrade.cxx:174-194`), through the
/// `BSplineCurve2dTo3d` / `BSplineCurve3dTo2d` round trip, which keeps the
/// periodic flag (`cxx:149-152`, `:169-172`) and re-tests `IsCN(1)` on the
/// lifted curve (`cxx:178`).
pub fn c0_bspline2d_to_sequence_of_c1_bspline_curve(bs: &dyn Curve2d) -> Option<Vec<Arc<dyn Curve2d>>> {
    if !bs.is_bspline2d() || !is_c0_bspline(bs.continuity()) {
        return None;
    }
    let poles2d = bs.poles2d()?;
    let degree = bs.bspline_degree()? as i32;
    // `cxx:183`: `BSplineCurve2dTo3d(BS)` (`cxx:133-152`), a `Geom_BSplineCurve`
    // over the same poles lifted with `z = 0`, the same weights, distinct knots,
    // multiplicities and `IsPeriodic()` flag.
    let mut poles: Vec<GpPnt> = poles2d.iter().map(|p| GpPnt::new(p.x(), p.y(), 0.0)).collect();
    let mut flat = bs.bspline_knots2d()?.to_vec();
    let mut weights = bs.bspline_weights2d().map(|w| w.to_vec());
    if bs.is_periodic() {
        // The 6-ary `Geom2d_BSplineCurve` constructor keeps the periodic
        // representation (`cxx:149-152`), and the 3D split then unperiodises
        // its own work copy (`cxx:50`).
        let mut work = Geom2dBSplineCurve::from_flat(
            poles.iter().map(|p| p.x()).collect(),
            poles.iter().map(|p| p.y()).collect(),
            weights,
            flat,
            degree as usize,
            true,
        )
        .ok()?;
        work.set_not_periodic();
        poles = work
            .xs
            .iter()
            .zip(work.ys.iter())
            .map(|(x, y)| GpPnt::new(*x, *y, 0.0))
            .collect();
        flat = work.knots.clone();
        weights = work.weights.clone();
    }
    // `cxx:54-57`: `Knots()`/`Multiplicities()` of the unperiodised work copy.
    let (uknots, umults) = bspl_knots::unique_knots_mults(&flat);
    let pieces = split_c0_bspline(degree, &flat, &poles, weights.as_deref(), &uknots, &umults)?;
    let mut out = Vec::with_capacity(pieces.len());
    for piece in &pieces {
        out.push(Arc::new(piece.to_2d().ok()?) as Arc<dyn Curve2d>);
    }
    Some(out)
}

/// `ShapeUpgrade.cxx:59-127`: the split loop. `flat` is `BS->KnotSequence()`
/// (knots repeated by multiplicity, non-periodic after `SetNotPeriodic`),
/// `uknots`/`umults` are `BS->Knots()`/`BS->Multiplicities()`, and `weights`
/// is `BS->WeightsArray()` (`None` for the unit weights of a non-rational
/// curve, `Geom_BSplineCurve.hxx:801-805`).
///
/// `None` is a piece the `Geom_BSplineCurve` constructor would reject.
fn split_c0_bspline(
    degree: i32,
    flat: &[f64],
    poles: &[GpPnt],
    weights: Option<&[f64]>,
    uknots: &[f64],
    umults: &[i32],
) -> Option<Vec<Piece>> {
    let nb_knots = uknots.len();
    if nb_knots == 0 || flat.is_empty() {
        return None;
    }
    // `cxx:60`: `StartKnotIndex = BS->FirstUKnotIndex()`.
    let last_uk = last_u_knot_index(degree, umults).max(1);
    let mut start = first_u_knot_index(degree, umults).max(1);
    let mut pieces = Vec::new();
    // `cxx:61`: `for (EndKnotIndex = StartKnotIndex + 1; End <= LastUKnotIndex; End++)`.
    while start + 1 <= last_uk {
        let mut end = start + 1;
        // `cxx:63-66`: skip knots whose multiplicity stays below the degree,
        // unless this is the last one.
        while end < last_uk && umults[(end - 1) as usize] < degree {
            end += 1;
        }
        // `cxx:69-71`.
        let start_flat = flat_index(degree, start, umults, false);
        let mut end_flat = flat_index(degree, end, umults, false);
        end_flat -= umults[(end - 1) as usize] - 1;

        // `cxx:73-91`: rebuild the knots and multiplicities of the span
        // `[StartFlatIndex - degree, EndFlatIndex + degree]` from the flat
        // sequence; runs of knots within `gp::Resolution()` merge.
        let mut temp_knots = vec![0.0f64; nb_knots];
        let mut temp_mults = vec![1i32; nb_knots];
        let mut temp_knot_index = 1i32;
        temp_knots[0] = *flat.get((start_flat - degree - 1) as usize)?;
        let mut j = start_flat - degree + 1;
        while j <= end_flat + degree {
            let a = *flat.get((j - 1) as usize)?;
            let b = *flat.get((j - 2) as usize)?;
            if (a - b).abs() <= REAL_SMALL {
                temp_mults[(temp_knot_index - 1) as usize] += 1;
            } else {
                if temp_knot_index as usize >= nb_knots {
                    return None;
                }
                temp_knot_index += 1;
                temp_knots[(temp_knot_index - 1) as usize] = a;
            }
            j += 1;
        }

        // `cxx:92-99`: raise a unit end multiplicity so the piece is clamped.
        let mut temp_start = 1i32;
        let mut temp_end = temp_knot_index;
        if nb_knots > 1 {
            if temp_mults[(temp_start - 1) as usize] == 1 {
                temp_start += 1;
                temp_mults[(temp_start - 1) as usize] += 1;
            }
            if temp_mults[(temp_end - 1) as usize] == 1 {
                temp_end -= 1;
                temp_mults[(temp_end - 1) as usize] += 1;
            }
        }
        if temp_end < temp_start {
            return None;
        }

        // `cxx:101-109`: slice the knots and multiplicities of the piece.
        let new_nb_knots = (temp_end - temp_start + 1) as usize;
        let mut new_mults = Vec::with_capacity(new_nb_knots);
        let mut new_knots = Vec::with_capacity(new_nb_knots);
        for k in 1..=new_nb_knots as i32 {
            let idx = (k + temp_start - 2) as usize;
            new_mults.push(*temp_mults.get(idx)?);
            new_knots.push(*temp_knots.get(idx)?);
        }

        // `cxx:111-120`: `BsplCLib::NbPoles(deg, false, newMults)` poles and
        // weights starting at `PoleIndex = StartFlatIndex - degree`.
        let new_nb_poles = bspl_knots::nb_poles(degree, false, &new_mults);
        if new_nb_poles <= 0 {
            return None;
        }
        let new_nb_poles = new_nb_poles as usize;
        let first_pole = (start_flat - degree - 1) as usize;
        if first_pole + new_nb_poles > poles.len() {
            return None;
        }
        let piece_poles = poles[first_pole..first_pole + new_nb_poles].to_vec();
        // `WeightsArray()` is always sized to `NbPoles`, so the slice exists
        // for a rational curve and is all ones otherwise.
        let piece_weights = match weights {
            Some(w) => {
                if first_pole + new_nb_poles > w.len() {
                    return None;
                }
                let sliced = w[first_pole..first_pole + new_nb_poles].to_vec();
                // `cxx:206-209`: `CheckRational` re-tests the sliced weights.
                if is_rational(&sliced) {
                    Some(sliced)
                } else {
                    None
                }
            }
            None => None,
        };
        pieces.push(Piece {
            poles: piece_poles,
            weights: piece_weights,
            knots: new_knots,
            mults: new_mults,
            degree,
        });

        // `cxx:128`: `StartKnotIndex = EndKnotIndex;`.
        start = end;
    }
    Some(pieces)
}

/// `Rational` (`Geom_BSplineCurve.cxx:97-106`, `Geom2d_BSplineCurve.cxx:95-107`):
/// the weights are rational when two neighbors differ by more than
/// `gp::Resolution()`.
fn is_rational(weights: &[f64]) -> bool {
    weights
        .windows(2)
        .any(|w| (w[0] - w[1]).abs() > REAL_SMALL)
}
