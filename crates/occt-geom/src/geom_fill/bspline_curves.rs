//! `GeomFill_BSplineCurves` — the B-spline surface bounded by two, three or
//! four curves.
//!
//! Source: `GeomFill_BSplineCurves.cxx:37-597` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `occ::handle<Geom_BSplineCurve>` becomes an owned [`GeomBSplineCurve`];
//!   `C->Copy()` is `Clone` and `C->Reversed()` is a clone run through
//!   `Curve::reverse` ([`reversed_curve`]).
//! - `StartPoint()` / `EndPoint()` are `Value(FirstParameter())` /
//!   `Value(LastParameter())` (`Geom_BSplineCurve_1.cxx:141-144`).
//! - `SetSameDistribution` keeps the OCCT sequence: reparametrize the shorter
//!   span, then insert each curve's knots into the other. OCCT builds the
//!   result with `BSplCLib::InsertKnots` on the raw pole arrays; the port calls
//!   `GeomBSplineCurve::InsertKnots`, which runs the same routine and keeps the
//!   rational flag.
//! - `GeomFill_Filling` is a handle to the abstract filling; the port takes a
//!   `&dyn Filling` built from one of [`crate::geom_fill::stretch::Stretch`],
//!   [`crate::geom_fill::coons::Coons`] or [`crate::geom_fill::curved::Curved`].

use occt_core::gp::GpPnt;
use occt_core::precision::{CONFUSION, PCONFUSION};

use crate::bspline_curve::GeomBSplineCurve;
use crate::bspline_surface::GeomBSplineSurface;
use crate::curve::Curve;
use crate::geom_fill::coons::Coons;
use crate::geom_fill::curved::Curved;
use crate::geom_fill::filling::Filling;
use crate::geom_fill::stretch::Stretch;

/// `GeomFill_FillingStyle` (`GeomFill_BSplineCurves.hxx:26-31`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FillingStyle {
    Stretch,
    Coons,
    Curved,
}

/// `GeomFill_BSplineCurves` (`GeomFill_BSplineCurves.hxx:33-80`).
#[derive(Clone, Default)]
pub struct BSplineCurves {
    surface: Option<GeomBSplineSurface>,
}

/// `Geom_BSplineCurve::Reversed()` (`Geom_BSplineCurve_1.cxx:160-168`).
pub(crate) fn reversed_curve(curve: &GeomBSplineCurve) -> GeomBSplineCurve {
    let mut copy = curve.clone();
    Curve::reverse(&mut copy);
    copy
}

/// `StartPoint()` (`Geom_BSplineCurve_1.cxx:141-144`).
fn start_point(curve: &GeomBSplineCurve) -> GpPnt {
    curve.d0(curve.first_parameter())
}

/// `EndPoint()` (`Geom_BSplineCurve_1.cxx:146-149`).
fn end_point(curve: &GeomBSplineCurve) -> GpPnt {
    curve.d0(curve.last_parameter())
}

/// `Knot(FirstUKnotIndex)` (`Geom_BSplineCurve_1.cxx:334-344`).
fn first_knot(curve: &GeomBSplineCurve) -> f64 {
    let (knots, mults) = curve.distinct_knots_and_mults();
    let index = occt_core::bspl::locate::first_u_knot_index(curve.degree() as i32, &mults);
    let index = index.clamp(1, knots.len() as i32) as usize;
    knots[index - 1]
}

/// `Knot(LastUKnotIndex)` (`Geom_BSplineCurve_1.cxx:334-344`).
fn last_knot(curve: &GeomBSplineCurve) -> f64 {
    let (knots, mults) = curve.distinct_knots_and_mults();
    let index = occt_core::bspl::locate::last_u_knot_index(curve.degree() as i32, &mults);
    let index = index.clamp(1, knots.len() as i32) as usize;
    knots[index - 1]
}

/// `Arrange` (`GeomFill_BSplineCurves.cxx:56-127`): put `C1..C4` in the order
/// `CC1 = bottom, CC2 = right, CC3 = top (reversed), CC4 = left (reversed)`.
///
/// A curve degenerate to a point is placed before the curvature leaves the
/// shared corner.
pub(crate) fn arrange_4(
    c1: &GeomBSplineCurve,
    c2: &GeomBSplineCurve,
    c3: &GeomBSplineCurve,
    c4: &GeomBSplineCurve,
    tol: f64,
) -> Option<[GeomBSplineCurve; 4]> {
    let mut gc = [c1.clone(), c2.clone(), c3.clone(), c4.clone()];
    // `cxx:74-123`: `i` is 1-based over the last three curves.
    for i in 1..=3usize {
        let mut trouve = false;
        // `cxx:78-92`: a degenerate curve that matches the last endpoint first.
        let mut j = i;
        while j <= 3 && !trouve {
            if start_point(&gc[j]).distance(&end_point(&gc[j])) < tol
                && start_point(&gc[j]).distance(&end_point(&gc[i - 1])) < tol
            {
                gc.swap(i, j);
                trouve = true;
            }
            j += 1;
        }
        // `cxx:95-114`: otherwise an ordinary curve, reversed when its end
        // matches.
        if !trouve {
            let mut j = i;
            while j <= 3 && !trouve {
                if start_point(&gc[j]).distance(&end_point(&gc[i - 1])) < tol {
                    gc.swap(i, j);
                    trouve = true;
                } else if end_point(&gc[j]).distance(&end_point(&gc[i - 1])) < tol {
                    gc[j] = reversed_curve(&gc[j]);
                    gc.swap(i, j);
                    trouve = true;
                }
                j += 1;
            }
        }
        // `cxx:116-121`.
        if !trouve {
            return None;
        }
    }
    // `cxx:123-127`.
    let cc3 = reversed_curve(&gc[2]);
    let cc4 = reversed_curve(&gc[3]);
    Some([gc[0].clone(), gc[1].clone(), cc3, cc4])
}

/// `SetSameDistribution` (`GeomFill_BSplineCurves.cxx:131-235`): make the two
/// curves share one knot vector and answer the common pole count.
pub(crate) fn set_same_distribution(
    c1: &mut GeomBSplineCurve,
    c2: &mut GeomBSplineCurve,
) -> Result<usize, String> {
    let mut k1 = c1.distinct_knots().to_vec();
    let m1 = c1.multiplicities().to_vec();
    let mut k2 = c2.distinct_knots().to_vec();
    let m2 = c2.multiplicities().to_vec();
    // `cxx:147-166`.
    let k11 = k1[0];
    let k12 = k1[k1.len() - 1];
    let k21 = k2[0];
    let k22 = k2[k2.len() - 1];
    if (k12 - k11) > (k22 - k21) {
        occt_core::bspl::knots::reparametrize(k11, k12, &mut k2);
        c2.set_knots(&k2)
            .map_err(|_| "Geom_BSplineCurve::SetKnots".to_string())?;
    } else if (k12 - k11) < (k22 - k21) {
        occt_core::bspl::knots::reparametrize(k21, k22, &mut k1);
        c1.set_knots(&k1)
            .map_err(|_| "Geom_BSplineCurve::SetKnots".to_string())?;
    } else if (k12 - k11).abs() > PCONFUSION {
        occt_core::bspl::knots::reparametrize(k11, k12, &mut k2);
        c2.set_knots(&k2)
            .map_err(|_| "Geom_BSplineCurve::SetKnots".to_string())?;
    }
    // `cxx:168-235`: insert each curve's knots into the other; a
    // `PrepareInsertKnots` failure is `Standard_ConstructionError`.
    c1.insert_knots(&k2, Some(&m2), PCONFUSION, false)
        .map_err(|_| "Standard_ConstructionError: GeomFill_BSplineCurves".to_string())?;
    c2.insert_knots(&k1, Some(&m1), PCONFUSION, false)
        .map_err(|_| "Standard_ConstructionError: GeomFill_BSplineCurves".to_string())?;
    Ok(c1.nb_poles())
}

/// `GeomFill_Filling` result assembled from the selected style.
fn build_filling(
    style: FillingStyle,
    p: [&[GpPnt]; 4],
    w: Option<[Vec<f64>; 4]>,
) -> Result<Box<dyn Filling>, String> {
    match (style, w) {
        (FillingStyle::Stretch, Some(w)) => Ok(Box::new(Stretch::new_rational(
            p[0], p[1], p[2], p[3], &w[0], &w[1], &w[2], &w[3],
        )?)),
        (FillingStyle::Stretch, None) => Ok(Box::new(Stretch::new(p[0], p[1], p[2], p[3])?)),
        // `cxx:347`: `GeomFill_Coons(P1, P4, P3, P2)`.
        (FillingStyle::Coons, Some(w)) => Ok(Box::new(Coons::new_rational(
            p[0], p[3], p[2], p[1], &w[0], &w[3], &w[2], &w[1],
        )?)),
        (FillingStyle::Coons, None) => Ok(Box::new(Coons::new(p[0], p[3], p[2], p[1])?)),
        (FillingStyle::Curved, Some(w)) => Ok(Box::new(Curved::new_rational(
            p[0], p[1], p[2], p[3], &w[0], &w[1], &w[2], &w[3],
        )?)),
        (FillingStyle::Curved, None) => Ok(Box::new(Curved::new(p[0], p[1], p[2], p[3])?)),
    }
}

impl BSplineCurves {
    /// `GeomFill_BSplineCurves(C1, C2, C3, C4, Type)`
    /// (`GeomFill_BSplineCurves.cxx:276-401`).
    pub fn new_4(
        c1: &GeomBSplineCurve,
        c2: &GeomBSplineCurve,
        c3: &GeomBSplineCurve,
        c4: &GeomBSplineCurve,
        style: FillingStyle,
    ) -> Result<Self, String> {
        // `cxx:282-287`.
        let mut curves = arrange_4(c1, c2, c3, c4, CONFUSION).ok_or_else(|| {
            "Standard_ConstructionError: GeomFill_BSplineCurves: Courbes non jointives".to_string()
        })?;
        // `cxx:290-314`: match the degrees along each direction.
        let deg1 = curves[0].degree();
        let deg2 = curves[1].degree();
        let deg3 = curves[2].degree();
        let deg4 = curves[3].degree();
        let deg_u = deg1.max(deg3);
        let deg_v = deg2.max(deg4);
        if deg1 < deg_u {
            curves[0]
                .increase_degree(deg_u)
                .map_err(|e| format!("Standard_ConstructionError: {e}"))?;
        }
        if deg2 < deg_v {
            curves[1]
                .increase_degree(deg_v)
                .map_err(|e| format!("Standard_ConstructionError: {e}"))?;
        }
        if deg3 < deg_u {
            curves[2]
                .increase_degree(deg_u)
                .map_err(|e| format!("Standard_ConstructionError: {e}"))?;
        }
        if deg4 < deg_v {
            curves[3]
                .increase_degree(deg_v)
                .map_err(|e| format!("Standard_ConstructionError: {e}"))?;
        }
        // `cxx:316-317`.
        let (u_side, v_side) = curves.split_at_mut(2);
        let nb_u_poles = set_same_distribution(&mut u_side[0], &mut v_side[0])?;
        let nb_v_poles = set_same_distribution(&mut u_side[1], &mut v_side[1])?;
        // `cxx:319-324`.
        if style == FillingStyle::Coons && (nb_u_poles < 4 || nb_v_poles < 4) {
            return Err(
                "Standard_ConstructionError: GeomFill_BSplineCurves: invalid filling style".into(),
            );
        }
        // `cxx:326-338`.
        let p1 = curves[0].poles.clone();
        let p2 = curves[1].poles.clone();
        let p3 = curves[2].poles.clone();
        let p4 = curves[3].poles.clone();
        let is_rat = curves[0].is_rational()
            || curves[1].is_rational()
            || curves[2].is_rational()
            || curves[3].is_rational();
        let weights = if is_rat {
            Some([
                weights_of(&curves[0]),
                weights_of(&curves[1]),
                weights_of(&curves[2]),
                weights_of(&curves[3]),
            ])
        } else {
            None
        };
        // `cxx:339-368`.
        let filling = build_filling(style, [&p1, &p2, &p3, &p4], weights)?;
        // `cxx:370-401`.
        let u_knots = curves[0].distinct_knots().to_vec();
        let u_mults = curves[0].multiplicities().to_vec();
        let v_knots = curves[1].distinct_knots().to_vec();
        let v_mults = curves[1].multiplicities().to_vec();
        let deg_u = curves[0].degree();
        let deg_v = curves[1].degree();
        let mut surface = GeomBSplineSurface::from_poles_knots_mults(
            filling.poles().to_vec(),
            u_knots,
            v_knots,
            u_mults,
            v_mults,
            deg_u,
            deg_v,
        )?;
        if filling.is_rational() {
            surface.weights = filling.weights().map(|w| w.to_vec());
        }
        Ok(Self {
            surface: Some(surface),
        })
    }

    /// `GeomFill_BSplineCurves(C1, C2, C3, Type)`
    /// (`GeomFill_BSplineCurves.cxx:406-441`): the fourth curve is the
    /// degree-1 segment closing `C1` and `C3` through `C2`'s knot span.
    pub fn new_3(
        c1: &GeomBSplineCurve,
        c2: &GeomBSplineCurve,
        c3: &GeomBSplineCurve,
        style: FillingStyle,
    ) -> Result<Self, String> {
        let tol = CONFUSION * CONFUSION;
        let mut poles = [GpPnt::zero(), GpPnt::zero()];
        // `cxx:417-425`.
        if start_point(c1).square_distance(&start_point(c2)) > tol
            && start_point(c1).square_distance(&end_point(c2)) > tol
        {
            poles[0] = start_point(c1);
        } else {
            poles[0] = end_point(c1);
        }
        // `cxx:427-435`.
        if start_point(c3).square_distance(&start_point(c2)) > tol
            && start_point(c3).square_distance(&end_point(c2)) > tol
        {
            poles[1] = start_point(c3);
        } else {
            poles[1] = end_point(c3);
        }
        // `cxx:436-439`.
        let knots = [first_knot(c2), last_knot(c2)];
        let mults = [2i32, 2];
        let c4 = GeomBSplineCurve::from_poles_knots_mults(poles.to_vec(), knots.to_vec(), mults.to_vec(), 1)
            .map_err(|e| format!("Standard_ConstructionError: GeomFill_BSplineCurves ({e})"))?;
        Self::new_4(c1, c2, c3, &c4, style)
    }

    /// `GeomFill_BSplineCurves(C1, C2, Type)`
    /// (`GeomFill_BSplineCurves.cxx:444-597`).
    pub fn new_2(
        c1: &GeomBSplineCurve,
        c2: &GeomBSplineCurve,
        style: FillingStyle,
    ) -> Result<Self, String> {
        let mut cc1 = c1.clone();
        let mut cc2 = c2.clone();
        let deg1 = cc1.degree();
        let deg2 = cc2.degree();
        let is_rat = cc1.is_rational() || cc2.is_rational();
        if style != FillingStyle::Curved {
            // `cxx:458-516`.
            let deg_u = deg1.max(deg2);
            if cc1.degree() < deg_u {
                cc1.increase_degree(deg_u)
                    .map_err(|e| format!("Standard_ConstructionError: {e}"))?;
            }
            if cc2.degree() < deg_u {
                cc2.increase_degree(deg_u)
                    .map_err(|e| format!("Standard_ConstructionError: {e}"))?;
            }
            let nb_poles = set_same_distribution(&mut cc1, &mut cc2)?;
            let p1 = cc1.poles.clone();
            let p2 = cc2.poles.clone();
            let mut poles = vec![vec![GpPnt::zero(); 2]; nb_poles];
            for i in 0..nb_poles {
                poles[i][0] = p1[i];
                poles[i][1] = p2[i];
            }
            let u_knots = cc1.distinct_knots().to_vec();
            let u_mults = cc1.multiplicities().to_vec();
            // `cxx:486-489`.
            let v_knots = vec![0.0, 1.0];
            let v_mults = vec![2i32, 2];
            let mut surface = if is_rat {
                // `cxx:492-512`: the V direction is never periodic.
                GeomBSplineSurface::from_poles_knots_mults(
                    poles,
                    u_knots,
                    v_knots,
                    u_mults,
                    v_mults,
                    cc1.degree(),
                    1,
                )?
            } else {
                // `cxx:514-516`.
                GeomBSplineSurface::from_poles_knots_mults(
                    poles,
                    u_knots,
                    v_knots,
                    u_mults,
                    v_mults,
                    cc1.degree(),
                    1,
                )?
            };
            if is_rat {
                let w1 = weights_of(&cc1);
                let w2 = weights_of(&cc2);
                let mut weights = vec![vec![0.0; 2]; nb_poles];
                for i in 0..nb_poles {
                    weights[i][0] = w1[i];
                    weights[i][1] = w2[i];
                }
                surface.weights = Some(weights);
            }
            return Ok(Self {
                surface: Some(surface),
            });
        }
        // `cxx:518-596`: `GeomFill_CurvedStyle` joins the two curves end to
        // end and extrudes one along the other.
        let eps = CONFUSION;
        let mut ok = false;
        if is_equal_eps(&start_point(&cc1), &start_point(&cc2), eps) {
            ok = true;
        } else if is_equal_eps(&start_point(&cc1), &end_point(&cc2), eps) {
            Curve::reverse(&mut cc2);
            ok = true;
        } else if is_equal_eps(&end_point(&cc1), &start_point(&cc2), eps) {
            // `cxx:533`: OCCT reverses the caller's `C1` here; the port reverses
            // the local copy so the argument stays untouched.
            Curve::reverse(&mut cc1);
            ok = true;
        } else if is_equal_eps(&end_point(&cc1), &end_point(&cc2), eps) {
            Curve::reverse(&mut cc1);
            Curve::reverse(&mut cc2);
            ok = true;
        }
        if !ok {
            return Err(
                "Standard_OutOfRange: GeomFill_BSplineCurves: Courbes non jointives".into(),
            );
        }
        let p1 = cc1.poles.clone();
        let p2 = cc2.poles.clone();
        let filling: Box<dyn Filling> = if is_rat {
            let w1 = weights_of(&cc1);
            let w2 = weights_of(&cc2);
            Box::new(Curved::new_two_rational(&p1, &p2, &w1, &w2)?)
        } else {
            Box::new(Curved::new_two(&p1, &p2)?)
        };
        // `cxx:553-558`.
        let u_knots = cc1.distinct_knots().to_vec();
        let u_mults = cc1.multiplicities().to_vec();
        let v_knots = cc2.distinct_knots().to_vec();
        let v_mults = cc2.multiplicities().to_vec();
        let mut surface = GeomBSplineSurface::from_poles_knots_mults(
            filling.poles().to_vec(),
            u_knots,
            v_knots,
            u_mults,
            v_mults,
            deg1,
            deg2,
        )?;
        if filling.is_rational() {
            surface.weights = filling.weights().map(|w| w.to_vec());
        }
        Ok(Self {
            surface: Some(surface),
        })
    }

    /// `GeomFill_BSplineCurves::Surface` (`GeomFill_BSplineCurves.hxx:60`).
    pub fn surface(&self) -> Option<&GeomBSplineSurface> {
        self.surface.as_ref()
    }

    /// The built surface, as the class's handle would be.
    pub fn into_surface(self) -> Result<GeomBSplineSurface, String> {
        self.surface
            .ok_or_else(|| "GeomFill_BSplineCurves: no surface".to_string())
    }
}

/// `WeightsArray()` for a curve that may be polynomial: OCCT would dereference
/// a null handle there, and the call sites only reach it when some curve is
/// rational, so the polynomial ones answer unit weights.
fn weights_of(curve: &GeomBSplineCurve) -> Vec<f64> {
    match &curve.weights {
        Some(w) => w.clone(),
        None => vec![1.0; curve.poles.len()],
    }
}

fn is_equal_eps(a: &GpPnt, b: &GpPnt, eps: f64) -> bool {
    a.distance(b) <= eps
}
