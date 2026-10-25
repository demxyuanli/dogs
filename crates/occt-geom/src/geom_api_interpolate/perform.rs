//! `GeomAPI_Interpolate::PerformNonPeriodic` / `PerformPeriodic`.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/GeomAPI/GeomAPI_Interpolate.cxx:434-808`
//! (OCCT 8.0.0).
//!
//! Both drivers fill a flat knot vector, a contact-order array and a pole array
//! whose entries alternate between the interpolation points and the end
//! tangents (each tangent adds one contact row of order 1), then solve the
//! whole system with `BSplCLib::Interpolate`
//! ([`occt_core::bspl::banded_interp::interpolate_contact`]).
//!
//! The local arrays are 1-based here, as in OCCT, so the index arithmetic of
//! the original is preserved: `parameters`, `flatknots`, `mults`,
//! `contact_order_array` and `poles` all start at slot 1.

use super::tangents::{build_periodic_tangent, build_tangents};
use super::*;
use occt_core::bspl::banded_interp::interpolate_contact;

#[inline]
fn get1(arr: &[f64], i: i32) -> f64 {
    arr[(i - 1) as usize]
}

#[inline]
fn set1(arr: &mut [f64], i: i32, v: f64) {
    arr[(i - 1) as usize] = v;
}

#[inline]
fn get1_mults(arr: &[i32], i: i32) -> i32 {
    arr[(i - 1) as usize]
}

#[inline]
fn set1_mults(arr: &mut [i32], i: i32, v: i32) {
    arr[(i - 1) as usize] = v;
}

/// `poles(index)` of the 3-dimensional pole array.
#[inline]
fn get_pole(poles: &[f64], i: i32) -> GpPnt {
    let base = ((i - 1) * 3) as usize;
    GpPnt::new(poles[base], poles[base + 1], poles[base + 2])
}

#[inline]
fn set_pole(poles: &mut [f64], i: i32, p: &GpPnt) {
    let base = ((i - 1) * 3) as usize;
    poles[base] = p.x();
    poles[base + 1] = p.y();
    poles[base + 2] = p.z();
}

/// `poles(index) = gp_Pnt(tangent.Coord(1..3))` — the tangent rows of the
/// interpolation system are stored as points (`cxx:523-527`, `:548-552`,
/// `:588-592`).
#[inline]
fn set_pole_tangent(poles: &mut [f64], i: i32, v: &GpVec) {
    let base = ((i - 1) * 3) as usize;
    poles[base] = v.x();
    poles[base + 1] = v.y();
    poles[base + 2] = v.z();
}

fn pole_vec(poles: &[f64], n: usize) -> Vec<GpPnt> {
    (1..=n as i32).map(|i| get_pole(poles, i)).collect()
}

impl GeomAPIInterpolate {
    /// `GeomAPI_Interpolate::PerformNonPeriodic`
    /// (`GeomAPI_Interpolate.cxx:613-808`).
    pub(super) fn perform_non_periodic(&mut self) {
        let num_points = self.points.len() as i32;
        let mut num_distinct_knots = num_points;
        let mut num_poles = num_points;
        let degree: i32;

        if num_poles == 2 && !self.tangent_request {
            degree = 1;
        } else if num_poles == 3 && !self.tangent_request {
            degree = 2;
            num_distinct_knots = 2;
        } else {
            degree = 3;
            num_poles += 2;
            if self.tangent_request {
                // `cxx:637-646`: interior flags only.
                for ii in 1..num_points - 1 {
                    if self.tangent_flags[ii as usize] {
                        num_poles += 1;
                    }
                }
            }
        }

        let np = num_poles as usize;
        let mut parameters = vec![0.0f64; np];
        let mut flatknots = vec![0.0f64; np + degree as usize + 1];
        let mut mults = vec![0i32; num_distinct_knots as usize];
        let mut knots = vec![0.0f64; num_distinct_knots as usize];
        let mut contact = vec![0i32; np];
        let mut poles = vec![0.0f64; np * 3];

        // `cxx:656-660`: the clamped ends.
        for ii in 1..=degree + 1 {
            set1(&mut flatknots, ii, get1(&self.parameters, 1));
            set1(&mut flatknots, ii + num_poles, get1(&self.parameters, num_points));
        }
        for ii in 2..num_distinct_knots {
            set1_mults(&mut mults, ii, 1);
        }
        set1_mults(&mut mults, 1, degree + 1);
        set1_mults(&mut mults, num_distinct_knots, degree + 1);

        match degree {
            1 => {
                // `cxx:672-681`.
                for ii in 1..=num_poles {
                    set_pole(&mut poles, ii, &self.points[(ii - 1) as usize]);
                }
                if let Ok(c) = GeomBSplineCurve::from_poles_knots_mults(
                    pole_vec(&poles, np),
                    self.parameters.clone(),
                    mults.clone(),
                    1,
                ) {
                    self.curve = Some(c);
                    self.is_done = true;
                }
            }
            2 => {
                // `cxx:682-699`.
                set1(&mut knots, 1, get1(&self.parameters, 1));
                set1(&mut knots, 2, get1(&self.parameters, 3));
                for ii in 1..=num_poles {
                    set_pole(&mut poles, ii, &self.points[(ii - 1) as usize]);
                }
                let params: Vec<f64> = (1..=num_points).map(|i| get1(&self.parameters, i)).collect();
                if interpolate_contact(2, &flatknots, &params, &contact, &mut poles, 3).is_ok() {
                    if let Ok(c) = GeomBSplineCurve::from_poles_knots_mults(
                        pole_vec(&poles, np),
                        knots.clone(),
                        mults.clone(),
                        2,
                    ) {
                        self.curve = Some(c);
                        self.is_done = true;
                    }
                }
            }
            3 => {
                // `cxx:701-803`.
                if num_points >= 3
                    && build_tangents(
                        &self.points,
                        &mut self.tangents,
                        &mut self.tangent_flags,
                        &self.parameters,
                    )
                    .is_err()
                {
                    return;
                }
                set1_mults(&mut contact, 2, 1);
                set1(&mut parameters, 1, get1(&self.parameters, 1));
                set1(&mut parameters, 2, get1(&self.parameters, 1));
                set_pole(&mut poles, 1, &self.points[0]);
                set_pole_tangent(&mut poles, 2, &self.tangents[0]);

                let mut mult_index = 2i32;
                let mut index = 3i32;
                let mut index1 = 2i32;
                let mut index2 = 2i32;
                let mut index3 = degree + 2;
                if self.tangent_request {
                    // `cxx:730-761`: `ii` runs over the interior parameters.
                    for ii in 2..num_points {
                        set1(&mut parameters, index, get1(&self.parameters, ii));
                        set_pole(&mut poles, index, &self.points[(index2 - 1) as usize]);
                        set1(&mut flatknots, index3, get1(&self.parameters, ii));
                        index += 1;
                        index3 += 1;
                        if self.tangent_flags[(index1 - 1) as usize] {
                            let v = get1_mults(&mults, mult_index) + 1;
                            set1_mults(&mut mults, mult_index, v);
                            contact[(index - 1) as usize] = 1;
                            set1(&mut flatknots, index3, get1(&self.parameters, ii));
                            set1(&mut parameters, index, get1(&self.parameters, ii));
                            let tangent = self.tangents[(ii - 1) as usize];
                            set_pole_tangent(&mut poles, index, &tangent);
                            index += 1;
                            index3 += 1;
                        }
                        mult_index += 1;
                        index1 += 1;
                        index2 += 1;
                    }
                } else {
                    // `cxx:762-781`.
                    index1 = 2;
                    for ii in 1..=num_points {
                        set1(&mut parameters, index1, get1(&self.parameters, ii));
                        index1 += 1;
                    }
                    index = 3;
                    for ii in 2..num_points {
                        set_pole(&mut poles, index, &self.points[(ii - 1) as usize]);
                        index += 1;
                    }
                    index = degree + 1;
                    for ii in 1..=num_points {
                        set1(&mut flatknots, index, get1(&self.parameters, ii));
                        index += 1;
                    }
                }

                // `cxx:783-796`: the last tangent row and the closing point.
                let tangent = self.tangents[(num_points - 1) as usize];
                set_pole_tangent(&mut poles, num_poles - 1, &tangent);
                set1_mults(&mut contact, num_poles - 1, 1);
                set1(&mut parameters, num_poles, get1(&self.parameters, num_points));
                set1(&mut parameters, num_poles - 1, get1(&self.parameters, num_points));
                set_pole(&mut poles, num_poles, &self.points[(num_points - 1) as usize]);

                let params: Vec<f64> =
                    (1..=num_poles).map(|i| get1(&parameters, i)).collect();
                if interpolate_contact(3, &flatknots, &params, &contact, &mut poles, 3).is_ok() {
                    let params: Vec<f64> =
                        (1..=num_points).map(|i| get1(&self.parameters, i)).collect();
                    if let Ok(c) = GeomBSplineCurve::from_poles_knots_mults(
                        pole_vec(&poles, np),
                        params,
                        mults.clone(),
                        3,
                    ) {
                        self.curve = Some(c);
                        self.is_done = true;
                    }
                }
            }
            _ => {}
        }
    }

    /// `GeomAPI_Interpolate::PerformPeriodic`
    /// (`GeomAPI_Interpolate.cxx:434-610`).
    pub(super) fn perform_periodic(&mut self) {
        let num_points = self.points.len() as i32;
        let period = get1(&self.parameters, self.parameters.len() as i32)
            - get1(&self.parameters, 1);
        let mut num_poles = num_points + 1;

        if num_points == 2 && !self.tangent_request {
            // `cxx:448-467`: a degree-1 periodic curve through two points.
            let mut deg1_mults = vec![0i32; num_poles as usize];
            for m in deg1_mults.iter_mut() {
                *m = 1;
            }
            if let Ok(c) = GeomBSplineCurve::from_poles_knots_mults_periodic(
                self.points.clone(),
                self.parameters.clone(),
                deg1_mults,
                1,
                true,
            ) {
                self.curve = Some(c);
                self.is_done = true;
            }
            return;
        }

        let num_distinct_knots = num_points + 1;
        let half_order = 2i32;
        let degree = 3i32;
        num_poles += 2;
        if self.tangent_request {
            // `cxx:474-483`: every flag from the second point on, including the
            // last one.
            for ii in 1..num_points {
                if self.tangent_flags[ii as usize] {
                    num_poles += 1;
                }
            }
        }

        let np = num_poles as usize;
        let mut parameters = vec![0.0f64; np];
        let mut flatknots = vec![0.0f64; np + degree as usize + 1];
        let mut mults = vec![0i32; num_distinct_knots as usize];
        let mut contact = vec![0i32; np];
        let mut poles = vec![0.0f64; np * 3];

        // `cxx:490-499`: the period-extended end knots.
        for ii in 1..=half_order {
            set1(&mut flatknots, ii, get1(&self.parameters, num_points) - period);
            set1(&mut flatknots, ii + half_order, get1(&self.parameters, 1));
            set1(&mut flatknots, num_poles + ii, get1(&self.parameters, num_points + 1));
            set1(
                &mut flatknots,
                num_poles + half_order + ii,
                get1(&self.parameters, half_order) + period,
            );
        }
        for ii in 2..num_distinct_knots {
            set1_mults(&mut mults, ii, 1);
        }
        set1_mults(&mut mults, 1, half_order);
        set1_mults(&mut mults, num_distinct_knots, half_order);

        if num_points >= 3
            && build_periodic_tangent(
                &self.points,
                &mut self.tangents,
                &mut self.tangent_flags,
                &self.parameters,
            )
            .is_err()
        {
            return;
        }
        set1_mults(&mut contact, 2, 1);
        set1(&mut parameters, 1, get1(&self.parameters, 1));
        set1(&mut parameters, 2, get1(&self.parameters, 1));
        set_pole(&mut poles, 1, &self.points[0]);
        set_pole_tangent(&mut poles, 2, &self.tangents[0]);

        let mut mult_index = 2i32;
        let mut index = 3i32;
        let mut index1 = degree + 2;
        if self.tangent_request {
            // `cxx:533-558`: `ii` runs from the second point to the last one.
            for ii in 2..=num_points {
                set1(&mut parameters, index, get1(&self.parameters, ii));
                set1(&mut flatknots, index1, get1(&self.parameters, ii));
                set_pole(&mut poles, index, &self.points[(ii - 1) as usize]);
                index += 1;
                index1 += 1;
                if self.tangent_flags[(ii - 1) as usize] {
                    let v = get1_mults(&mults, mult_index) + 1;
                    set1_mults(&mut mults, mult_index, v);
                    contact[(index - 1) as usize] = 1;
                    set1(&mut parameters, index, get1(&self.parameters, ii));
                    set1(&mut flatknots, index1, get1(&self.parameters, ii));
                    let tangent = self.tangents[(ii - 1) as usize];
                    set_pole_tangent(&mut poles, index, &tangent);
                    index += 1;
                    index1 += 1;
                }
                mult_index += 1;
            }
        } else {
            // `cxx:559-580`.
            index = degree + 1;
            index1 = 2;
            for ii in 1..=num_points + 1 {
                set1(&mut parameters, index1, get1(&self.parameters, ii));
                set1(&mut flatknots, index, get1(&self.parameters, ii));
                index += 1;
                index1 += 1;
            }
            index = 3;
            for ii in 2..=num_points {
                set_pole(&mut poles, index, &self.points[(ii - 1) as usize]);
                index += 1;
            }
        }

        // `cxx:581-596`: the wrap-around tangent row and closing point.
        set1_mults(&mut contact, num_poles - 1, 1);
        set1(
            &mut parameters,
            num_poles - 1,
            get1(&self.parameters, num_points + 1),
        );
        set_pole_tangent(&mut poles, num_poles - 1, &self.tangents[0]);
        set1(
            &mut parameters,
            num_poles,
            get1(&self.parameters, num_points + 1),
        );
        set_pole(&mut poles, num_poles, &self.points[0]);

        let params: Vec<f64> = (1..=num_poles).map(|i| get1(&parameters, i)).collect();
        if interpolate_contact(degree as usize, &flatknots, &params, &contact, &mut poles, 3)
            .is_err()
        {
            return;
        }

        // `cxx:604-608`: `myCurve = new Geom_BSplineCurve(poles(1), 1,
        // num_poles - 2)` — the last two solved poles close the curve and are
        // not part of the periodic representation.
        let kept = num_poles - 2;
        let newpoles: Vec<GpPnt> = (1..=kept).map(|i| get_pole(&poles, i)).collect();
        let knot_values: Vec<f64> = (1..=num_distinct_knots)
            .map(|i| get1(&self.parameters, i))
            .collect();
        let mut_knots: Vec<i32> = (1..=num_distinct_knots)
            .map(|i| get1_mults(&mults, i))
            .collect();
        if let Ok(c) = GeomBSplineCurve::from_poles_knots_mults_periodic(
            newpoles,
            knot_values,
            mut_knots,
            degree as usize,
            true,
        ) {
            self.curve = Some(c);
            self.is_done = true;
        }
    }
}




