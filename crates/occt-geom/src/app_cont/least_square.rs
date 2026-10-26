//! `AppCont_LeastSquare` - continuous least-squares fit of an
//! [`AppCont_Function`](crate::app_cont::Function) by one Bezier polynomial.
//!
//! Source: `ModelingData/TKGeomBase/AppCont/AppCont_LeastSquare.cxx:28-611`
//! and `AppCont_LeastSquare.hxx:24-67` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `math_Matrix` stays [`MathMatrix`] (1-based `Value(row, col)`), so the
//!   `myPoints` / `myPoles` / `myVB` indexing of the `cxx` is preserved
//!   verbatim. `math_Vector myParam` becomes a `Vec<f64>` indexed from 0 with a
//!   `- 1` on the 1-based reads.
//! - The five `AppCont_ContMatrices` tables are evaluated at run time by
//!   [`cont_matrices`](crate::app_cont::cont_matrices); see that module's note
//!   for how the tables were rebuilt.
//! - `aTabP` / `aTabP2d` are sized `max(nb, 1)` in OCCT because
//!   `NCollection_Array1` cannot be empty; the port keeps the same allocation
//!   so the `Value` / `D1` output slice lengths match the OCCT call.
//! - `PeriodicityInfo` (`AppCont_LeastSquare.hxx:28-32`) is
//!   [`PeriodicityInfo`].

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_math::{MathMatrix, gauss_legendre};

use crate::app_cont::cont_matrices;
use crate::app_cont::function::Function;
use crate::app_par_curves::constraint::Constraint;
use crate::app_par_curves::multi_curve::MultiCurve;
use crate::app_par_curves::multi_point::MultiPoint;

/// `PeriodicityInfo` (`AppCont_LeastSquare.hxx:28-32`).
#[derive(Debug, Clone, Copy)]
pub struct PeriodicityInfo {
    /// `isPeriodic`.
    pub is_periodic: bool,
    /// `myPeriod`.
    pub my_period: f64,
}

/// `AppCont_LeastSquare` (`AppCont_LeastSquare.hxx:34-67`).
///
/// `myParam` and `myPerInfo` (`hxx:62-64`) are private and read only inside
/// the constructor, so the port keeps them as constructor locals instead of
/// struct fields; nothing outside `cxx:91-515` observes them.
#[derive(Debug, Clone)]
pub struct LeastSquare {
    /// `mySCU` (`hxx:59`).
    my_scu: MultiCurve,
    /// `myPoints` (`hxx:60`).
    my_points: MathMatrix,
    /// `myPoles` (`hxx:61`).
    my_poles: MathMatrix,
    /// `myVB` (`hxx:63`).
    my_vb: MathMatrix,
    /// `myDone` (`hxx:65`).
    my_done: bool,
    /// `myDegre` (`hxx:66`).
    my_degre: i32,
    /// `myNbdiscret` (`hxx:67`).
    my_nbdiscret: i32,
    /// `myNbP` (`hxx:67`).
    my_nb_p: i32,
    /// `myNbP2d` (`hxx:67`).
    my_nb_p2d: i32,
}

impl LeastSquare {
    /// `AppCont_LeastSquare(SSP, U0, U1, FirstCons, LastCons, Deg, NbPoints)`
    /// (`cxx:91-515`).
    pub fn new<F: Function + ?Sized>(
        ssp: &F,
        u0: f64,
        u1: f64,
        first_cons: Constraint,
        last_cons: Constraint,
        deg: i32,
        nb_points: i32,
    ) -> Self {
        let (nb_p, nb_p2d) = ssp.get_number_of_points();
        let classe = deg + 1;
        let cl1 = deg;
        let nbcol = (3 * nb_p + 2 * nb_p2d) as usize;
        let first_p = 1i32;
        let last_p = nb_points;

        let my_scu = MultiCurve::with_nb_poles(classe as usize);
        let mut my_points = MathMatrix::with_init(1, nb_points as usize, 1, nbcol, 0.0);
        let mut my_poles = MathMatrix::with_init(1, classe as usize, 1, nbcol, 0.0);
        let mut my_param = vec![0.0f64; nb_points as usize];
        let mut my_per_info = vec![
            PeriodicityInfo {
                is_periodic: false,
                my_period: 0.0,
            };
            nbcol
        ];
        let mut my_done = false;
        let mut b = MathMatrix::with_init(1, classe as usize, 1, nbcol, 0.0);
        let mut bdeb = 1i32;
        let mut bfin = classe;
        let mut my_first_c = first_cons;
        let mut my_last_c = last_cons;

        for a_dim_idx in 1..=(3 * nb_p + 2 * nb_p2d) {
            let (is_periodic, my_period) = ssp.period_information(a_dim_idx);
            my_per_info[(a_dim_idx - 1) as usize] = PeriodicityInfo {
                is_periodic,
                my_period,
            };
        }

        let n3 = nb_p.max(1) as usize;
        let n2 = nb_p2d.max(1) as usize;
        let mut a_tab_p = vec![GpPnt::zero(); n3];
        let mut a_tab_p2d = vec![GpPnt2d::zero(); n2];
        let mut a_tab_v = vec![GpVec::zero(); n3];
        let mut a_tab_v2d = vec![GpVec2d::default(); n2];

        if my_first_c == Constraint::TangencyPoint && !ssp.d1(u0, &mut a_tab_v2d, &mut a_tab_v) {
            my_first_c = Constraint::PassPoint;
        }
        if my_last_c == Constraint::TangencyPoint && !ssp.d1(u1, &mut a_tab_v2d, &mut a_tab_v) {
            my_last_c = Constraint::PassPoint;
        }

        // Compute control points params on which approximation will be built.
        let (gauss_p, gauss_w) = gauss_legendre(-1.0, 1.0, nb_points as usize);
        let mut the_weights = vec![0.0f64; nb_points as usize];
        let d_u = 0.5 * (u1 - u0);
        let half = (nb_points + 1) / 2;
        for i in first_p..=last_p {
            let u = 0.5 * (u1 + u0) + d_u * gauss_p[(i - 1) as usize];
            if i <= half {
                let pos = last_p - i + 1;
                my_param[(pos - 1) as usize] = u;
                the_weights[(pos - 1) as usize] = 0.5 * gauss_w[(i - 1) as usize];
            } else {
                let pos = i - half;
                my_param[(pos - 1) as usize] = u;
                the_weights[(pos - 1) as usize] = 0.5 * gauss_w[(i - 1) as usize];
            }
        }

        // Compute control points.
        for i in first_p..=last_p {
            let u = my_param[(i - 1) as usize];
            ssp.value(u, &mut a_tab_p2d, &mut a_tab_p);

            let mut i2 = 1i32;
            for j in 1..=nb_p {
                let p = a_tab_p[(j - 1) as usize];
                my_points.set_value(i as usize, i2 as usize, p.x());
                my_points.set_value(i as usize, (i2 + 1) as usize, p.y());
                my_points.set_value(i as usize, (i2 + 2) as usize, p.z());
                i2 += 3;
            }
            for j in 1..=nb_p2d {
                let p = a_tab_p2d[(j - 1) as usize];
                my_points.set_value(i as usize, i2 as usize, p.x());
                my_points.set_value(i as usize, (i2 + 1) as usize, p.y());
                i2 += 2;
            }
        }

        // Fix possible "period jump".
        let a_max_dim = 3 * nb_p + 2 * nb_p2d;
        for a_dim_idx in 1..=a_max_dim {
            let ud = a_dim_idx as usize;
            let info = my_per_info[ud - 1];
            if info.is_periodic
                && (my_points.value(1, ud) - my_points.value(2, ud)).abs()
                    > info.my_period / 2.01
                && (my_points.value(2, ud) - my_points.value(3, ud)).abs()
                    < info.my_period / 2.01
            {
                let period_mult = if my_points.value(1, ud) < my_points.value(2, ud) {
                    1.0
                } else {
                    -1.0
                };
                let new_param = my_points.value(1, ud) + period_mult * info.my_period;
                my_points.set_value(1, ud, new_param);
            }
        }
        for a_pnt_idx in 1..nb_points {
            for a_dim_idx in 1..=a_max_dim {
                let ud = a_dim_idx as usize;
                let info = my_per_info[ud - 1];
                let p0 = my_points.value(a_pnt_idx as usize, ud);
                let p1 = my_points.value((a_pnt_idx + 1) as usize, ud);
                if info.is_periodic && (p0 - p1).abs() > info.my_period / 2.01 {
                    let period_mult = if p0 > p1 { 1.0 } else { -1.0 };
                    let new_param = p1 + period_mult * info.my_period;
                    my_points
                        .set_value((a_pnt_idx + 1) as usize, ud, new_param);
                }
            }
        }

        let my_vb = cont_matrices::v_bernstein(classe, nb_points);

        // Processing of the right-hand side:
        let mut tmppoints = vec![0.0f64; nbcol];
        for c in 1..=classe {
            for t in tmppoints.iter_mut() {
                *t = 0.0;
            }
            for i in 1..=nb_points {
                let coeff = the_weights[(i - 1) as usize] * my_vb.value(c as usize, i as usize);
                for j in 1..=nbcol as i32 {
                    tmppoints[(j - 1) as usize] +=
                        my_points.value(i as usize, j as usize) * coeff;
                }
            }
            for k in 1..=nbcol as i32 {
                let v = b.value(c as usize, k as usize) + tmppoints[(k - 1) as usize];
                b.set_value(c as usize, k as usize, v);
            }
        }

        if my_first_c == Constraint::NoConstraint && my_last_c == Constraint::NoConstraint {
            let inv_m = cont_matrices::inv_m_matrix(classe);
            // Direct computation of poles:
            for i in 1..=classe {
                for j in 1..=classe {
                    let ib_ij = inv_m.value(i as usize, j as usize);
                    for k in 1..=nbcol as i32 {
                        let v = my_poles.value(i as usize, k as usize)
                            + ib_ij * b.value(j as usize, k as usize);
                        my_poles.set_value(i as usize, k as usize, v);
                    }
                }
            }
        } else {
            let m = cont_matrices::m_matrix(classe);
            let mut a_fix_p2d = vec![GpPnt2d::zero(); n2];
            let mut a_fix_p = vec![GpPnt::zero(); n3];

            if my_first_c == Constraint::PassPoint || my_first_c == Constraint::TangencyPoint {
                ssp.value(u0, &mut a_tab_p2d, &mut a_tab_p);
                Self::fix_single_border_point(
                    ssp,
                    u0,
                    u0,
                    u1,
                    &mut a_fix_p2d,
                    &mut a_fix_p,
                );

                let mut i2 = 1i32;
                for k in 1..=nb_p {
                    let picked = if a_fix_p[(k - 1) as usize].distance(&a_tab_p[(k - 1) as usize])
                        > 0.1
                    {
                        a_fix_p[(k - 1) as usize]
                    } else {
                        a_tab_p[(k - 1) as usize]
                    };
                    my_poles.set_value(1, i2 as usize, picked.x());
                    my_poles.set_value(1, (i2 + 1) as usize, picked.y());
                    my_poles.set_value(1, (i2 + 2) as usize, picked.z());
                    i2 += 3;
                }
                for k in 1..=nb_p2d {
                    let picked = if a_fix_p2d[(k - 1) as usize]
                        .distance(&a_tab_p2d[(k - 1) as usize])
                        > 0.1
                    {
                        a_fix_p2d[(k - 1) as usize]
                    } else {
                        a_tab_p2d[(k - 1) as usize]
                    };
                    my_poles.set_value(1, i2 as usize, picked.x());
                    my_poles.set_value(1, (i2 + 1) as usize, picked.y());
                    i2 += 2;
                }

                for a_dim_idx in 1..=a_max_dim {
                    let ud = a_dim_idx as usize;
                    let info = my_per_info[ud - 1];
                    if info.is_periodic
                        && (my_poles.value(1, ud) - my_points.value(1, ud)).abs()
                            > info.my_period / 2.01
                    {
                        let mult = if my_poles.value(1, ud) < my_points.value(1, ud) {
                            1.0
                        } else {
                            -1.0
                        };
                        let v = my_poles.value(1, ud) + mult * info.my_period;
                        my_poles.set_value(1, ud, v);
                    }
                }
            }

            if my_last_c == Constraint::PassPoint || my_last_c == Constraint::TangencyPoint {
                ssp.value(u1, &mut a_tab_p2d, &mut a_tab_p);
                Self::fix_single_border_point(
                    ssp,
                    u1,
                    u0,
                    u1,
                    &mut a_fix_p2d,
                    &mut a_fix_p,
                );

                let mut i2 = 1i32;
                for k in 1..=nb_p {
                    let picked = if a_fix_p[(k - 1) as usize].distance(&a_tab_p[(k - 1) as usize])
                        > 0.1
                    {
                        a_fix_p[(k - 1) as usize]
                    } else {
                        a_tab_p[(k - 1) as usize]
                    };
                    my_poles.set_value(classe as usize, i2 as usize, picked.x());
                    my_poles.set_value(classe as usize, (i2 + 1) as usize, picked.y());
                    my_poles.set_value(classe as usize, (i2 + 2) as usize, picked.z());
                    i2 += 3;
                }
                for k in 1..=nb_p2d {
                    let picked = if a_fix_p2d[(k - 1) as usize]
                        .distance(&a_tab_p2d[(k - 1) as usize])
                        > 0.1
                    {
                        a_fix_p2d[(k - 1) as usize]
                    } else {
                        a_tab_p2d[(k - 1) as usize]
                    };
                    my_poles.set_value(classe as usize, i2 as usize, picked.x());
                    my_poles.set_value(classe as usize, (i2 + 1) as usize, picked.y());
                    i2 += 2;
                }

                for a_dim_idx in 1..=2 {
                    let ud = a_dim_idx as usize;
                    let info = my_per_info[ud - 1];
                    let last_pole = my_poles.value(classe as usize, ud);
                    let last_point = my_points.value(nb_points as usize, ud);
                    if info.is_periodic && (last_pole - last_point).abs() > info.my_period / 2.01 {
                        let mult = if last_pole < last_point { 1.0 } else { -1.0 };
                        my_poles
                            .set_value(classe as usize, ud, last_pole + mult * info.my_period);
                    }
                }
            }

            if my_first_c == Constraint::PassPoint {
                bdeb = 2;
                // Update the right-hand side:
                for i in 1..=classe {
                    let coeff = m.value(i as usize, 1);
                    for k in 1..=nbcol as i32 {
                        let v = b.value(i as usize, k as usize)
                            - my_poles.value(1, k as usize) * coeff;
                        b.set_value(i as usize, k as usize, v);
                    }
                }
            }

            if my_last_c == Constraint::PassPoint {
                bfin = cl1;
                for i in 1..=classe {
                    let coeff = m.value(i as usize, classe as usize);
                    for k in 1..=nbcol as i32 {
                        let v = b.value(i as usize, k as usize)
                            - my_poles.value(classe as usize, k as usize) * coeff;
                        b.set_value(i as usize, k as usize, v);
                    }
                }
            }

            if my_first_c == Constraint::TangencyPoint {
                // Fix the second pole:
                bdeb = 3;
                ssp.d1(u0, &mut a_tab_v2d, &mut a_tab_v);

                let mut i2 = 1i32;
                let coeff = (u1 - u0) / deg as f64;
                for k in 1..=nb_p {
                    let i2plus1 = i2 + 1;
                    let i2plus2 = i2 + 2;
                    let v = a_tab_v[(k - 1) as usize];
                    let p0 = my_poles.value(1, i2 as usize);
                    let p1 = my_poles.value(1, i2plus1 as usize);
                    let p2 = my_poles.value(1, i2plus2 as usize);
                    my_poles.set_value(2, i2 as usize, p0 + v.x() * coeff);
                    my_poles.set_value(2, i2plus1 as usize, p1 + v.y() * coeff);
                    my_poles.set_value(2, i2plus2 as usize, p2 + v.z() * coeff);
                    i2 += 3;
                }
                for k in 1..=nb_p2d {
                    let i2plus1 = i2 + 1;
                    let v = a_tab_v2d[(k - 1) as usize];
                    let p0 = my_poles.value(1, i2 as usize);
                    let p1 = my_poles.value(1, i2plus1 as usize);
                    my_poles.set_value(2, i2 as usize, p0 + v.x() * coeff);
                    my_poles.set_value(2, i2plus1 as usize, p1 + v.y() * coeff);
                    i2 += 2;
                }

                for i in 1..=classe {
                    let coeff1 = m.value(i as usize, 1);
                    let coeff2 = m.value(i as usize, 2);
                    for k in 1..=nbcol as i32 {
                        let v = b.value(i as usize, k as usize)
                            - my_poles.value(1, k as usize) * coeff1
                            - my_poles.value(2, k as usize) * coeff2;
                        b.set_value(i as usize, k as usize, v);
                    }
                }
            }

            if my_last_c == Constraint::TangencyPoint {
                bfin = classe - 2;
                ssp.d1(u1, &mut a_tab_v2d, &mut a_tab_v);
                let mut i2 = 1i32;
                let coeff = (u1 - u0) / deg as f64;
                for k in 1..=nb_p {
                    let i2plus1 = i2 + 1;
                    let i2plus2 = i2 + 2;
                    let v = a_tab_v[(k - 1) as usize];
                    let p0 = my_poles.value(classe as usize, i2 as usize);
                    let p1 = my_poles.value(classe as usize, i2plus1 as usize);
                    let p2 = my_poles.value(classe as usize, i2plus2 as usize);
                    my_poles.set_value(cl1 as usize, i2 as usize, p0 - v.x() * coeff);
                    my_poles
                        .set_value(cl1 as usize, i2plus1 as usize, p1 - v.y() * coeff);
                    my_poles
                        .set_value(cl1 as usize, i2plus2 as usize, p2 - v.z() * coeff);
                    i2 += 3;
                }
                for k in 1..=nb_p2d {
                    let i2plus1 = i2 + 1;
                    let v = a_tab_v2d[(k - 1) as usize];
                    let p0 = my_poles.value(classe as usize, i2 as usize);
                    let p1 = my_poles.value(classe as usize, i2plus1 as usize);
                    my_poles.set_value(cl1 as usize, i2 as usize, p0 - v.x() * coeff);
                    my_poles
                        .set_value(cl1 as usize, i2plus1 as usize, p1 - v.y() * coeff);
                    i2 += 2;
                }

                for i in 1..=classe {
                    let coeff1 = m.value(i as usize, classe as usize);
                    let coeff2 = m.value(i as usize, cl1 as usize);
                    for k in 1..=nbcol as i32 {
                        let v = b.value(i as usize, k as usize)
                            - my_poles.value(classe as usize, k as usize) * coeff1
                            - my_poles.value(cl1 as usize, k as usize) * coeff2;
                        b.set_value(i as usize, k as usize, v);
                    }
                }
            }

            if bdeb <= bfin {
                let mut b2 = MathMatrix::with_init(
                    bdeb as usize,
                    bfin as usize,
                    1,
                    b.col_count(),
                    0.0,
                );

                for i in bdeb..=bfin {
                    for j in 1..=classe {
                        let coeff = m.value(i as usize, j as usize);
                        for k in 1..=nbcol as i32 {
                            let v = b2.value(i as usize, k as usize)
                                + b.value(j as usize, k as usize) * coeff;
                            b2.set_value(i as usize, k as usize, v);
                        }
                    }
                }

                // Resolution:
                // ===========
                let dim = (bfin - bdeb + 1) as usize;
                let ibp: MathMatrix = if bdeb == 2 && bfin == classe - 1 && classe <= 26 {
                    cont_matrices::ibp_matrix(classe)
                } else if bdeb == 3 && bfin == classe - 2 && classe <= 26 {
                    cont_matrices::ibt_matrix(classe)
                } else {
                    // `MP.Transposed() * MP` (`cxx:487-497`), accumulated
                    // directly since `math_Matrix::Transposed` is square-only
                    // here and `MP` is not.
                    let mut ibp1 =
                        MathMatrix::with_init(bdeb as usize, bfin as usize, bdeb as usize, bfin as usize, 0.0);
                    for i in bdeb..=bfin {
                        for j in bdeb..=bfin {
                            let mut sum = 0.0;
                            for k in 1..=classe {
                                sum += m.value(k as usize, i as usize)
                                    * m.value(k as usize, j as usize);
                            }
                            ibp1.set_value(i as usize, j as usize, sum);
                        }
                    }
                    ibp1
                        .invert()
                        .expect("Standard_DimensionError: AppCont_LeastSquare IBP");
                    ibp1
                };
                debug_assert_eq!(ibp.row_count(), dim);

                my_done = true;
                for i in bdeb..=bfin {
                    for j in bdeb..=bfin {
                        let ibp_ij = ibp.value(i as usize, j as usize);
                        for k in 1..=nbcol as i32 {
                            let v = my_poles.value(i as usize, k as usize)
                                + ibp_ij * b2.value(j as usize, k as usize);
                            my_poles.set_value(i as usize, k as usize, v);
                        }
                    }
                }
            }
        }

        Self {
            my_scu,
            my_points,
            my_poles,
            my_vb,
            my_done,
            my_degre: deg,
            my_nbdiscret: nb_points,
            my_nb_p: nb_p,
            my_nb_p2d: nb_p2d,
        }
    }

    /// `FixSingleBorderPoint` (`cxx:30-87`).
    fn fix_single_border_point<F: Function + ?Sized>(
        ssp: &F,
        the_u: f64,
        the_u0: f64,
        the_u1: f64,
        the_fix2d: &mut [GpPnt2d],
        the_fix: &mut [GpPnt],
    ) {
        let a_max_iter = 15;
        let mut a_tab_p = vec![GpPnt::zero(); the_fix.len().max(1)];
        let mut a_prev_p = vec![GpPnt::zero(); the_fix.len().max(1)];
        let mut a_tab_p2d = vec![GpPnt2d::zero(); the_fix2d.len().max(1)];
        let mut a_prev_p2d = vec![GpPnt2d::zero(); the_fix2d.len().max(1)];
        let a_mult = if (the_u - the_u0) > (the_u1 - the_u) {
            1.0
        } else {
            -1.0
        };
        let a_start_param = the_u;
        let mut a_prev_dist = 1.0;
        let mut a_curr_dist = 1.0;

        let mut du = -(the_u1 - the_u0) / 2.0 * a_mult;
        let eps = f64::EPSILON;
        let dec = 0.1;
        for an_iter in 1..a_max_iter {
            du *= dec;
            let a_curr_param = a_start_param + du;
            ssp.value(a_curr_param, &mut a_tab_p2d, &mut a_tab_p);

            // from second iteration
            if an_iter > 1 {
                a_curr_dist = 0.0;
                for j in 0..a_tab_p.len() {
                    a_curr_dist += a_tab_p[j].distance(&a_prev_p[j]);
                }
                for j in 0..a_tab_p2d.len() {
                    a_curr_dist += a_tab_p2d[j].distance(&a_prev_p2d[j]);
                }
                // from the third iteration
                if an_iter > 2 && a_curr_dist / a_prev_dist > 10.0 {
                    break;
                }
            }
            a_prev_p.copy_from_slice(&a_tab_p);
            a_prev_p2d.copy_from_slice(&a_tab_p2d);
            a_prev_dist = a_curr_dist;
            if a_prev_dist <= eps {
                break;
            }
        }
        the_fix2d.copy_from_slice(&a_prev_p2d[..the_fix2d.len()]);
        the_fix.copy_from_slice(&a_prev_p[..the_fix.len()]);
    }

    /// `Value()` (`cxx:518-548`): the poles as a `MultiCurve`.
    pub fn value(&mut self) -> &MultiCurve {
        let ideb = 1;
        let ifin = self.my_degre + 1;
        for i in ideb..=ifin {
            let mut mpole = MultiPoint::with_size(self.my_nb_p as usize, self.my_nb_p2d as usize);
            let mut j2 = 1i32;
            for j in 1..=self.my_nb_p {
                let pt = GpPnt::new(
                    self.my_poles.value(i as usize, j2 as usize),
                    self.my_poles.value(i as usize, (j2 + 1) as usize),
                    self.my_poles.value(i as usize, (j2 + 2) as usize),
                );
                mpole.set_point(j, &pt);
                j2 += 3;
            }
            let nb_p = self.my_nb_p;
            for j in (nb_p + 1)..=(nb_p + self.my_nb_p2d) {
                let pt2d = GpPnt2d::new(
                    self.my_poles.value(i as usize, j2 as usize),
                    self.my_poles.value(i as usize, (j2 + 1) as usize),
                );
                mpole.set_point2d(j, &pt2d);
                j2 += 2;
            }
            self.my_scu.set_value(i, &mpole);
        }
        &self.my_scu
    }

    /// `Error(F, MaxE3d, MaxE2d)` (`cxx:551-606`).
    pub fn error(&self) -> (f64, f64, f64) {
        let classe = self.my_degre + 1;
        let mut f = 0.0f64;
        let mut max_e3d = 0.0f64;
        let mut max_e2d = 0.0f64;
        let ncol = self.my_points.col_count();
        let mut my_points = self.my_points.clone();
        let mut tmppoles = vec![0.0f64; ncol];

        for c in 1..=classe {
            for k in 1..=ncol as i32 {
                tmppoles[(k - 1) as usize] = self.my_poles.value(c as usize, k as usize);
            }
            for i in 1..=self.my_nbdiscret {
                let coeff = self.my_vb.value(c as usize, i as usize);
                for j in 1..=ncol as i32 {
                    let v = my_points.value(i as usize, j as usize)
                        - tmppoles[(j - 1) as usize] * coeff;
                    my_points.set_value(i as usize, j as usize, v);
                }
            }
        }

        for i in 1..=self.my_nbdiscret {
            let mut i2 = 1i32;
            for _j in 1..=self.my_nb_p {
                let e1 = my_points.value(i as usize, i2 as usize);
                let e2 = my_points.value(i as usize, (i2 + 1) as usize);
                let e3 = my_points.value(i as usize, (i2 + 2) as usize);
                let err3d = e1 * e1 + e2 * e2 + e3 * e3;
                max_e3d = max_e3d.max(err3d);
                f += err3d;
                i2 += 3;
            }
            for _j in 1..=self.my_nb_p2d {
                let e1 = my_points.value(i as usize, i2 as usize);
                let e2 = my_points.value(i as usize, (i2 + 1) as usize);
                let err2d = e1 * e1 + e2 * e2;
                max_e2d = max_e2d.max(err2d);
                f += err2d;
                i2 += 2;
            }
        }

        (f, max_e3d.sqrt(), max_e2d.sqrt())
    }

    /// `IsDone()` (`cxx:611-614`).
    pub fn is_done(&self) -> bool {
        self.my_done
    }
}
