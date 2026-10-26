//! `AppDef_Variational::Optimization` and `AssemblingConstraints` — one
//! finite-element solve plus the linear constraint assembly.
//!
//! Source: `AppDef_Variational.cxx:1593-1689` (`Optimization`) and
//! `cxx:2993-3313` (`AssemblingConstraints`) (OCCT 8.0.0).
//!
//! Translation notes:
//! - `math_Vector G0(0, MxDeg)` / `math_Matrix H(0, MxDeg, 0, MxDeg)` keep the
//!   zero lower bound, so the linear forms handed to
//!   [`Assembly::add_constraint`] start at index 0 as in OCCT.

use std::cell::RefCell;
use std::rc::Rc;

use occt_math::fem_tool::{Assembly, Curve};
use occt_math::matrix::MathMatrix;
use occt_math::vector::MathVector;

use super::Variational;

impl Variational {
    /// `Optimization(J, A, ToAssemble, EpsDeg, Curve, Parameters)`
    /// (`cxx:1593-1689`).
    pub(super) fn optimization(
        &mut self,
        a: &mut Assembly,
        to_assemble: bool,
        eps_deg: f64,
        curve: &Rc<RefCell<Curve>>,
        parameters: &MathVector,
    ) {
        let (mx_deg, nb_elm, nb_dim) = {
            let c = curve.borrow();
            (
                c.base().work_degree(),
                c.nb_elements(),
                c.dimension() as i32,
            )
        };

        let ass_table = a.assembly_table().to_vec();
        let mut h = MathMatrix::new(0, mx_deg as usize, 0, mx_deg as usize);
        let mut g = MathVector::new(0, mx_deg as usize);
        let mut sol = MathVector::new(1, a.nb_glob_var());

        let nb_constr = self.nb_pass_points + self.nb_tang_points + self.nb_curv_points;
        let cblong = *self.smooth_criterion.est_length();

        // Updating the assembly (`cxx:1615-1638`).
        if to_assemble {
            a.nullify_matrix();
        }
        a.nullify_vector();

        for el in 1..=nb_elm {
            if to_assemble {
                self.smooth_criterion.hessian(el, 1, 1, &mut h);
                for dim in 1..=nb_dim {
                    a.add_matrix(el, dim as usize, dim as usize, &h);
                }
            }
            for dim in 1..=nb_dim {
                self.smooth_criterion.gradient(el, dim, &mut g);
                a.add_vector(el, dim as usize, &g);
            }
        }

        // Solution of the system (`cxx:1640-1650`).
        if to_assemble {
            if nb_constr != 0 {
                self.assembling_constraints(curve, parameters, cblong, a);
            }
            a.solve();
        }
        a.solution(&mut sol);

        // Updating J (`cxx:1652-1654`).
        self.smooth_criterion.set_curve(curve);
        self.smooth_criterion.input_vector(&sol, &ass_table);

        // Updating the curve and reducing the degree (`cxx:1658-1688`).
        let mut c = curve.borrow_mut();
        if nb_constr == 0 {
            for el in 1..=nb_elm {
                let _ = c.reduce_degree(el, eps_deg);
            }
        } else {
            let tab_int: Vec<f64> = c.knots().to_vec();
            let mut icnt = 1i32;
            let p0 = parameters.lower() as i32 - self.first_point;
            for el in 1..=nb_elm {
                while icnt < nb_constr
                    && parameters.value(
                        (p0 + self.typ_constraints[(2 * icnt - 1) as usize]) as usize,
                    ) <= tab_int[el]
                {
                    icnt += 1;
                }
                let point = (p0 + self.typ_constraints[(2 * icnt - 1) as usize]) as usize;
                if parameters.value(point) <= tab_int[el] || parameters.value(point) >= tab_int[el + 1]
                {
                    let _ = c.reduce_degree(el, eps_deg);
                } else if c.degree(el) < mx_deg {
                    c.set_degree(el, mx_deg);
                }
            }
        }
    }

    /// `AssemblingConstraints(Curve, Parameters, CBLONG, A)`
    /// (`cxx:2993-3313`).
    fn assembling_constraints(
        &self,
        curve: &Rc<RefCell<Curve>>,
        parameters: &MathVector,
        cblong: f64,
        a: &mut Assembly,
    ) {
        let (mx_deg, nb_elm, nb_dim) = {
            let c = curve.borrow();
            (
                c.base().work_degree(),
                c.nb_elements(),
                c.dimension() as i32,
            )
        };
        let n_theta = 6 * self.nb_p3d + 2 * self.nb_p2d;
        let nb_constr = self.nb_pass_points + self.nb_tang_points + self.nb_curv_points;

        let ng3d = 3 * nb_constr + 3 * self.nb_tang_points + 5 * self.nb_curv_points;
        let ng2d = 2 * nb_constr + 2 * self.nb_tang_points + 3 * self.nb_curv_points;
        let nbeg2d = ng3d * self.nb_p3d;
        let ngpc1 = nb_constr + self.nb_tang_points + self.nb_curv_points;
        let mut ntang3d = 3 * ngpc1;
        let mut ntang2d = 2 * ngpc1;

        let intervals: Vec<f64> = curve.borrow().knots().to_vec();
        let order = curve.borrow().base().niv_constr() + 1;
        let base = curve.borrow().base().clone();

        let mut g0 = vec![0.0f64; mx_deg as usize + 1];
        let mut g1 = vec![0.0f64; mx_deg as usize + 1];
        let mut g2 = vec![0.0f64; mx_deg as usize + 1];

        let p0 = parameters.lower() as i32 - self.first_point;

        a.nullify_constraint();

        let mut curel = 1usize;
        let mut npass = 0i32;
        let mut ipnt = -1i32;
        let mut ityp = 0i32;

        for i in 1..=nb_constr {
            ipnt += 2;
            ityp += 2;

            let point = self.typ_constraints[ipnt as usize];
            let typ_of_constr = self.typ_constraints[ityp as usize];

            let mut t = parameters.value((p0 + point) as usize);

            let mut el = curel;
            while el <= nb_elm {
                el += 1;
                if t <= intervals[el] {
                    curel = el - 1;
                    break;
                }
            }

            let u_first = intervals[curel];
            let u_last = intervals[curel + 1];
            let coeff = (u_last - u_first) / 2.0;
            let c0 = (u_last + u_first) / 2.0;

            t = (t - c0) / coeff;

            if typ_of_constr == 0 {
                base.d0(t, &mut g0);
                for k in 1..order as usize {
                    let mfact = coeff.powi(k as i32);
                    g0[k] *= mfact;
                    g0[k + order as usize] *= mfact;
                }
            } else if typ_of_constr == 1 {
                base.d1(t, &mut g0, &mut g1);
                for k in 1..order as usize {
                    let mfact = coeff.powi(k as i32);
                    g0[k] *= mfact;
                    g0[k + order as usize] *= mfact;
                    g1[k] *= mfact;
                    g1[k + order as usize] *= mfact;
                }
                let mfact = 1.0 / coeff;
                for k in 0..=mx_deg as usize {
                    g1[k] *= mfact;
                }
            } else {
                base.d2(t, &mut g0, &mut g1, &mut g2);
                for k in 1..order as usize {
                    let mfact = coeff.powi(k as i32);
                    g0[k] *= mfact;
                    g0[k + order as usize] *= mfact;
                    g1[k] *= mfact;
                    g1[k + order as usize] *= mfact;
                    g2[k] *= mfact;
                    g2[k + order as usize] *= mfact;
                }
                let mfact = 1.0 / coeff;
                let mfact1 = mfact / coeff;
                for k in 0..=mx_deg as usize {
                    g1[k] *= mfact;
                    g2[k] *= mfact1;
                }
            }

            npass += 1;

            let mut j = nb_dim * (point - self.first_point);
            let mut n0 = npass;
            let mut curdim = 0i32;
            for _pnt in 1..=self.nb_p3d {
                let mut index_of_constraint = n0;
                for k in 1..=3 {
                    curdim += 1;
                    let lf = linear_form(&g0, 1.0);
                    a.add_constraint(
                        index_of_constraint as usize,
                        curel,
                        curdim as usize,
                        &lf,
                        self.tab_points[(j + k) as usize],
                    );
                    index_of_constraint += ngpc1;
                }
                j += 3;
                n0 += ng3d;
            }

            n0 = npass + nbeg2d;
            for _pnt in 1..=self.nb_p2d {
                let mut index_of_constraint = n0;
                for k in 1..=2 {
                    curdim += 1;
                    let lf = linear_form(&g0, 1.0);
                    a.add_constraint(
                        index_of_constraint as usize,
                        curel,
                        curdim as usize,
                        &lf,
                        self.tab_points[(j + k) as usize],
                    );
                    index_of_constraint += ngpc1;
                }
                j += 2;
                n0 += ng2d;
            }

            if typ_of_constr == 1 {
                // Tangency: unit tangency is imposed in the arc-length scale.
                npass += 1;
                n0 = npass;
                j = 2 * nb_dim * (i - 1);
                curdim = 0;
                for _pnt in 1..=self.nb_p3d {
                    let mut index_of_constraint = n0;
                    for k in 1..=3 {
                        curdim += 1;
                        let lf = linear_form(&g1, 1.0);
                        a.add_constraint(
                            index_of_constraint as usize,
                            curel,
                            curdim as usize,
                            &lf,
                            cblong * self.tab_constraints[(j + k) as usize],
                        );
                        index_of_constraint += ngpc1;
                    }
                    n0 += ng3d;
                    j += 6;
                }

                n0 = npass + nbeg2d;
                for _pnt in 1..=self.nb_p2d {
                    let mut index_of_constraint = n0;
                    for k in 1..=2 {
                        curdim += 1;
                        let lf = linear_form(&g1, 1.0);
                        a.add_constraint(
                            index_of_constraint as usize,
                            curel,
                            curdim as usize,
                            &lf,
                            cblong * self.tab_constraints[(j + k) as usize],
                        );
                        index_of_constraint += ngpc1;
                    }
                    n0 += ng2d;
                    j += 4;
                }
            }

            if typ_of_constr == 2 {
                // Curvature: tangency plus the second-derivative rows.
                npass += 1;
                n0 = npass;
                j = 2 * nb_dim * (i - 1);
                curdim = 0;
                for _pnt in 1..=self.nb_p3d {
                    let mut index_of_constraint = n0;
                    for k in 1..=3 {
                        curdim += 1;
                        let lf = linear_form(&g1, 1.0);
                        a.add_constraint(
                            index_of_constraint as usize,
                            curel,
                            curdim as usize,
                            &lf,
                            cblong * self.tab_constraints[(j + k) as usize],
                        );
                        index_of_constraint += ngpc1;
                    }
                    n0 += ng3d;
                    j += 6;
                }

                n0 = npass + nbeg2d;
                for _pnt in 1..=self.nb_p2d {
                    let mut index_of_constraint = n0;
                    for k in 1..=2 {
                        curdim += 1;
                        let lf = linear_form(&g1, 1.0);
                        a.add_constraint(
                            index_of_constraint as usize,
                            curel,
                            curdim as usize,
                            &lf,
                            cblong * self.tab_constraints[(j + k) as usize],
                        );
                        index_of_constraint += ngpc1;
                    }
                    n0 += ng2d;
                    j += 4;
                }

                j = 2 * nb_dim * (i - 1) + 3;
                let mut jt = n_theta * (i - 1);
                let mut index_of_constraint = ntang3d + 1;
                curdim = 0;
                for _pnt in 1..=self.nb_p3d {
                    let mut r1 = 0.0;
                    let mut r2 = 0.0;
                    for k in 1..=3 {
                        r1 += self.tab_constraints[(j + k) as usize]
                            * self.ttheta[(jt + k) as usize];
                        r2 += self.tab_constraints[(j + k) as usize]
                            * self.ttheta[(jt + 3 + k) as usize];
                    }
                    r1 *= cblong * cblong;
                    r2 *= cblong * cblong;
                    for k in 1..=3 {
                        curdim += 1;
                        if k > 1 {
                            r1 = 0.0;
                            r2 = 0.0;
                        }
                        let lf1 =
                            linear_form(&g2, self.tfthet[(jt + k) as usize]);
                        a.add_constraint(
                            index_of_constraint as usize,
                            curel,
                            curdim as usize,
                            &lf1,
                            r1,
                        );
                        let lf2 =
                            linear_form(&g2, self.tfthet[(jt + 3 + k) as usize]);
                        a.add_constraint(
                            (index_of_constraint + 1) as usize,
                            curel,
                            curdim as usize,
                            &lf2,
                            r2,
                        );
                    }
                    index_of_constraint += ng3d;
                    j += 6;
                    jt += 6;
                }

                j -= 1;
                index_of_constraint = nbeg2d + ntang2d + 1;
                for _pnt in 1..=self.nb_p2d {
                    let mut r1 = 0.0;
                    for k in 1..=2 {
                        r1 += self.tab_constraints[(j + k) as usize]
                            * self.ttheta[(jt + k) as usize];
                    }
                    r1 *= cblong * cblong;
                    for k in 1..=2 {
                        curdim += 1;
                        if k > 1 {
                            r1 = 0.0;
                        }
                        let lf = linear_form(&g2, self.tfthet[(jt + k) as usize]);
                        a.add_constraint(
                            index_of_constraint as usize,
                            curel,
                            curdim as usize,
                            &lf,
                            r1,
                        );
                    }
                    index_of_constraint += ng2d;
                    j += 4;
                    jt += 2;
                }

                ntang3d += 2;
                ntang2d += 1;
            }
        }
    }
}

/// `math_Vector(0, MxDeg)` filled from a zero-based slice and scaled
/// (`myTfthet->Value(...) * V2` in `cxx:3273`).
fn linear_form(values: &[f64], scalar: f64) -> MathVector {
    let mut v = MathVector::new(0, values.len().saturating_sub(1));
    for (i, value) in values.iter().enumerate() {
        v.set_value(i, value * scalar);
    }
    v
}
