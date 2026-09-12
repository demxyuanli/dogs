//! `math_FunctionRoots` NEWCODE path.
//! Source: `math_FunctionRoots.cxx:38-252` (helpers) and `:261-743` (ctor).
//! `MATH_FUNCTIONROOTS_OLDCODE` is not compiled in OCCT.

use crate::math_bracketed_root::{
    bracketed_root, BRACKETED_ROOT_NB_ITERATIONS, BRACKETED_ROOT_ZEPS,
};
use crate::math_fn::{MathFunction, MathFunctionWithDerivative};

const ITMAX: i32 = 100;
const EPSEPS: f64 = 2e-14;

/// `math_FunctionRoots` (`math_FunctionRoots.hxx:31-83`).
#[derive(Debug, Clone)]
pub struct FunctionRoots {
    done: bool,
    all_null: bool,
    sol: Vec<f64>,
    nb_state_sol: Vec<i32>,
}

impl FunctionRoots {
    /// `math_FunctionRoots::math_FunctionRoots` NEWCODE (`cxx:261-743`).
    pub fn new<F: MathFunctionWithDerivative>(
        f: &mut F,
        a: f64,
        b: f64,
        nb_sample: i32,
        eps_x: f64,
        eps_f: f64,
        eps_null: f64,
        k: f64,
    ) -> Self {
        let mut sol = Vec::new();
        let mut nb_state_sol = Vec::new();
        let eps_x_ctor = eps_x;
        let mut x0 = a;
        let mut xn = b;
        let mut n = nb_sample;
        if b < a {
            x0 = b;
            xn = a;
        }
        n *= 2;
        if n < 20 {
            n = 20;
        }
        let mut eps_x = eps_x;
        let delta_u = x0.abs() + xn.abs();
        let n_eps_x = 0.0000000001 * delta_u;
        if eps_x < n_eps_x {
            eps_x = n_eps_x;
        }
        let mut x = x0;
        let dx = (xn - x0) / (n as f64);
        let mut ptrval = vec![0.0; (n as usize) + 1];
        let mut nvalid: i32 = -1;
        let mut aux = 0.0;
        for _ in 0..=n {
            if x > xn {
                x = xn;
            }
            let ok = f.value(x, &mut aux);
            if ok {
                nvalid += 1;
                ptrval[nvalid as usize] = aux - k;
            }
            x += dx;
        }
        if nvalid < n {
            return Self {
                done: false,
                all_null: false,
                sol,
                nb_state_sol,
            };
        }
        let mut all_null = true;
        let mut i = 0;
        while all_null && i <= n {
            if ptrval[i as usize] > eps_null || ptrval[i as usize] < -eps_null {
                all_null = false;
            }
            i += 1;
        }
        if !all_null {
            let tol = eps_x;
            x = x0;
            for i in 0..n {
                let ip1 = i + 1;
                let mut x2 = x + dx;
                if x2 > xn {
                    x2 = xn;
                }
                if ptrval[i as usize] < 0.0 {
                    if ptrval[ip1 as usize] > 0.0 {
                        solve(
                            f,
                            k,
                            x,
                            ptrval[i as usize],
                            x2,
                            ptrval[ip1 as usize],
                            tol,
                            n_eps_x,
                            &mut sol,
                            &mut nb_state_sol,
                        );
                    }
                } else if ptrval[ip1 as usize] < 0.0 {
                    solve(
                        f,
                        k,
                        x,
                        ptrval[i as usize],
                        x2,
                        ptrval[ip1 as usize],
                        tol,
                        n_eps_x,
                        &mut sol,
                        &mut nb_state_sol,
                    );
                }
                x += dx;
            }
            for i in 0..=n {
                if ptrval[i as usize] == 0.0 {
                    let mut xx = x0 + (i as f64) * dx;
                    if xx > xn {
                        xx = xn;
                    }
                    let mut u0 = dx * 0.5;
                    let mut u1 = xx + u0;
                    u0 += xx;
                    if u0 < x0 {
                        u0 = x0;
                    }
                    if u0 > xn {
                        u0 = xn;
                    }
                    if u1 < x0 {
                        u1 = x0;
                    }
                    if u1 > xn {
                        u1 = xn;
                    }
                    let mut y0 = 0.0;
                    let mut y1 = 0.0;
                    f.value(u0, &mut y0);
                    y0 -= k;
                    f.value(u1, &mut y1);
                    y1 -= k;
                    if y0 * y1 < 0.0 {
                        solve(f, k, u0, y0, u1, y1, tol, n_eps_x, &mut sol, &mut nb_state_sol);
                    } else if y0 != 0.0 || y1 != 0.0 {
                        append_root(&mut sol, &mut nb_state_sol, xx, f, k, n_eps_x);
                    }
                }
            }
            if ptrval[0] <= eps_f && ptrval[0] >= -eps_f {
                append_root(&mut sol, &mut nb_state_sol, x0, f, k, n_eps_x);
            }
            if ptrval[n as usize] <= eps_f && ptrval[n as usize] >= -eps_f {
                append_root(&mut sol, &mut nb_state_sol, xn, f, k, n_eps_x);
            }
            let majdx = 5.0 * dx;
            let mut im1 = 0;
            let mut ip1 = 2;
            let mut xm = x0 + dx;
            for i in 1..n {
                let mut rediscr = false;
                if xm > xn {
                    xm = xn;
                }
                if ptrval[i as usize] > 0.0 {
                    if ptrval[im1 as usize] > ptrval[i as usize]
                        && ptrval[ip1 as usize] > ptrval[i as usize]
                    {
                        let mut xm1 = xm - dx;
                        if xm1 < x0 {
                            xm1 = x0;
                        }
                        let mut ym = 0.0;
                        let mut dym = 0.0;
                        f.values(xm1, &mut ym, &mut dym);
                        ym -= k;
                        if dym < -1e-10 || dym > 1e-10 {
                            let t = ym / dym;
                            if t < majdx && t > -majdx {
                                rediscr = true;
                            }
                        }
                        if !rediscr {
                            let mut xp1 = xm + dx;
                            if xp1 > xn {
                                xp1 = xn;
                            }
                            f.values(xp1, &mut ym, &mut dym);
                            ym -= k;
                            if dym < -1e-10 || dym > 1e-10 {
                                let t = ym / dym;
                                if t < majdx && t > -majdx {
                                    rediscr = true;
                                }
                            }
                        }
                    }
                } else if ptrval[i as usize] < 0.0
                    && ptrval[im1 as usize] < ptrval[i as usize]
                    && ptrval[ip1 as usize] < ptrval[i as usize]
                {
                    let mut xm1 = xm - dx;
                    if xm1 < x0 {
                        xm1 = x0;
                    }
                    let mut ym = 0.0;
                    let mut dym = 0.0;
                    f.values(xm1, &mut ym, &mut dym);
                    ym -= k;
                    if dym > 1e-10 || dym < -1e-10 {
                        let t = ym / dym;
                        if t < majdx && t > -majdx {
                            rediscr = true;
                        }
                    }
                    if !rediscr {
                        let mut xm1 = xm - dx;
                        if xm1 < x0 {
                            xm1 = x0;
                        }
                        f.values(xm1, &mut ym, &mut dym);
                        ym -= k;
                        if dym > 1e-10 || dym < -1e-10 {
                            let t = ym / dym;
                            if t < majdx && t > -majdx {
                                rediscr = true;
                            }
                        }
                    }
                }
                if rediscr {
                    rediscr_extrema(
                        f,
                        k,
                        eps_x_ctor,
                        eps_f,
                        n_eps_x,
                        tol,
                        x0,
                        xn,
                        xm,
                        dx,
                        ptrval[im1 as usize],
                        ptrval[ip1 as usize],
                        &mut sol,
                        &mut nb_state_sol,
                    );
                }
                xm += dx;
                im1 += 1;
                ip1 += 1;
            }
        }
        Self {
            done: true,
            all_null,
            sol,
            nb_state_sol,
        }
    }

    pub fn is_done(&self) -> bool {
        self.done
    }

    pub fn is_all_null(&self) -> bool {
        self.all_null
    }

    pub fn nb_solutions(&self) -> i32 {
        self.sol.len() as i32
    }

    /// `Value(Nieme)` — 1-based like OCCT.
    pub fn value(&self, nieme: i32) -> f64 {
        self.sol[(nieme as usize) - 1]
    }

    /// `StateNumber(Nieme)` — 1-based like OCCT.
    pub fn state_number(&self, nieme: i32) -> i32 {
        self.nb_state_sol[(nieme as usize) - 1]
    }
}

/// `DerivFunction` (`cxx:38-49`) plus the extrema rediscretization (`cxx:558-741`).
fn rediscr_extrema<F: MathFunctionWithDerivative>(
    f: &mut F,
    k: f64,
    eps_x_ctor: f64,
    eps_f: f64,
    n_eps_x: f64,
    tol: f64,
    x0_dom: f64,
    xn_dom: f64,
    xm: f64,
    dx: f64,
    f0_sample: f64,
    f3_sample: f64,
    sol: &mut Vec<f64>,
    nb_state_sol: &mut Vec<i32>,
) {
    let mut x0 = xm - dx;
    let mut x3 = xm + dx;
    if x0 < x0_dom {
        x0 = x0_dom;
    }
    if x3 > xn_dom {
        x3 = xn_dom;
    }
    let mut a_sol_x1 = 0.0;
    let mut a_sol_x2 = 0.0;
    let mut a_val1 = 0.0;
    let mut a_val2 = 0.0;
    let mut a_der1 = 0.0;
    let mut a_der2 = 0.0;
    let mut is_sol1 = false;
    let mut is_sol2 = false;
    {
        let mut der_f = DerivFunction { inner: f };
        let a_br = bracketed_root(
            &mut der_f,
            x0,
            x3,
            eps_x_ctor,
            BRACKETED_ROOT_NB_ITERATIONS,
            BRACKETED_ROOT_ZEPS,
        );
        if a_br.done {
            a_sol_x1 = a_br.root;
            let mut val = 0.0;
            f.value(a_sol_x1, &mut val);
            a_val1 = val.abs();
            if a_val1 < eps_f {
                is_sol1 = true;
                a_der1 = a_br.value;
            }
        }
    }
    let r = 0.61803399;
    let c = 1.0 - r;
    let tol_cr = n_eps_x * 10.0;
    let mut f0 = f0_sample;
    let mut f3 = f3_sample;
    let recherche_minimum = f0 > 0.0;
    let (mut x1, mut x2) = if (x3 - xm).abs() > (x0 - xm).abs() {
        (xm, xm + c * (x3 - xm))
    } else {
        (xm - c * (xm - x0), xm)
    };
    let mut f1 = 0.0;
    let mut f2 = 0.0;
    f.value(x1, &mut f1);
    f1 -= k;
    f.value(x2, &mut f2);
    f2 -= k;
    let tol_x = 0.001 * n_eps_x;
    while (x3 - x0).abs() > tol_cr * (x1.abs() + x2.abs()) && (x1 - x2).abs() > tol_x {
        if recherche_minimum {
            if f2 < f1 {
                x0 = x1;
                x1 = x2;
                x2 = r * x1 + c * x3;
                f0 = f1;
                f1 = f2;
                f.value(x2, &mut f2);
                f2 -= k;
            } else {
                x3 = x2;
                x2 = x1;
                x1 = r * x2 + c * x0;
                f3 = f2;
                f2 = f1;
                f.value(x1, &mut f1);
                f1 -= k;
            }
        } else if f2 > f1 {
            x0 = x1;
            x1 = x2;
            x2 = r * x1 + c * x3;
            f0 = f1;
            f1 = f2;
            f.value(x2, &mut f2);
            f2 -= k;
        } else {
            x3 = x2;
            x2 = x1;
            x1 = r * x2 + c * x0;
            f3 = f2;
            f2 = f1;
            f.value(x1, &mut f1);
            f1 -= k;
        }
        if f1 * f0 < 0.0 {
            solve(f, k, x0, f0, x1, f1, tol, n_eps_x, sol, nb_state_sol);
        }
        if f2 * f3 < 0.0 {
            solve(f, k, x2, f2, x3, f3, tol, n_eps_x, sol, nb_state_sol);
        }
    }
    if (recherche_minimum && f1 < f2) || (!recherche_minimum && f1 > f2) {
        if f1.abs() < eps_f {
            is_sol2 = true;
            a_sol_x2 = x1;
            a_val2 = f1.abs();
        }
    } else if f2.abs() < eps_f {
        is_sol2 = true;
        a_sol_x2 = x2;
        a_val2 = f2.abs();
    }
    if is_sol1 && is_sol2 {
        if a_val2 - a_val1 > eps_f {
            append_root(sol, nb_state_sol, a_sol_x1, f, k, n_eps_x);
        } else if a_val1 - a_val2 > eps_f {
            append_root(sol, nb_state_sol, a_sol_x2, f, k, n_eps_x);
        } else {
            a_der1 = a_der1.abs();
            f.derivative(a_sol_x2, &mut a_der2);
            a_der2 = a_der2.abs();
            if a_der1 < a_der2 {
                append_root(sol, nb_state_sol, a_sol_x1, f, k, n_eps_x);
            } else {
                append_root(sol, nb_state_sol, a_sol_x2, f, k, n_eps_x);
            }
        }
    } else if is_sol1 {
        append_root(sol, nb_state_sol, a_sol_x1, f, k, n_eps_x);
    } else if is_sol2 {
        append_root(sol, nb_state_sol, a_sol_x2, f, k, n_eps_x);
    }
}

/// `DerivFunction` (`cxx:38-49`): `Value` forwards to `Derivative`.
struct DerivFunction<'a, F: MathFunctionWithDerivative> {
    inner: &'a mut F,
}

impl<F: MathFunctionWithDerivative> MathFunction for DerivFunction<'_, F> {
    fn value(&mut self, x: f64, fval: &mut f64) -> bool {
        self.inner.derivative(x, fval)
    }
}

/// `AppendRoot` (`cxx:51-108`). `k` is unused in OCCT.
fn append_root<F: MathFunctionWithDerivative>(
    sol: &mut Vec<f64>,
    nb_state_sol: &mut Vec<i32>,
    x: f64,
    f: &mut F,
    _k: f64,
    d_x: f64,
) {
    let n = sol.len();
    let mut t = 0.0;
    if n == 0 {
        sol.push(x);
        f.value(x, &mut t);
        nb_state_sol.push(f.get_state_number());
        return;
    }
    let mut i = 1;
    let mut pl = n + 1;
    while i <= n {
        t = sol[i - 1];
        if t >= x {
            pl = i;
            i = n;
        }
        if (x - t).abs() <= d_x {
            pl = 0;
            i = n;
        }
        i += 1;
    }
    if pl > n {
        sol.push(x);
        f.value(x, &mut t);
        nb_state_sol.push(f.get_state_number());
    } else if pl > 0 {
        sol.insert(pl - 1, x);
        f.value(x, &mut t);
        nb_state_sol.insert(pl - 1, f.get_state_number());
    }
}

/// `Solve` Brent + Newton polish (`cxx:110-252`).
fn solve<F: MathFunctionWithDerivative>(
    f: &mut F,
    k: f64,
    x1: f64,
    y1: f64,
    x2: f64,
    y2: f64,
    tol: f64,
    d_x: f64,
    sol: &mut Vec<f64>,
    nb_state_sol: &mut Vec<i32>,
) {
    let tols2 = 0.5 * tol;
    let mut a = x1;
    let mut b = x2;
    let mut c = x2;
    let mut d = 0.0;
    let mut e = 0.0;
    let mut fa = y1;
    let mut fb = y2;
    let mut fc = y2;
    for _iter in 1..=ITMAX {
        if (fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0) {
            c = a;
            fc = fa;
            e = b - a;
            d = e;
        }
        if fc.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        let tol1 = EPSEPS * b.abs() + tols2;
        let xm = 0.5 * (c - b);
        if xm.abs() < tol1 || fb == 0.0 {
            let mut xp = b;
            let mut itern = 5;
            let mut ok;
            loop {
                let mut yp = 0.0;
                let mut dp = 0.0;
                ok = f.values(xp, &mut yp, &mut dp);
                if ok {
                    ok = false;
                    if dp > 1e-10 || dp < -1e-10 {
                        xp -= (yp - k) / dp;
                    }
                    if xp <= x2 && xp >= x1 {
                        f.value(xp, &mut yp);
                        yp -= k;
                        if yp.abs() < fb.abs() {
                            b = xp;
                            fb = yp;
                            ok = true;
                        }
                    }
                }
                itern -= 1;
                if !(ok && itern >= 0) {
                    break;
                }
            }
            append_root(sol, nb_state_sol, b, f, k, d_x);
            return;
        }
        if e.abs() >= tol1 && fa.abs() > fb.abs() {
            let s = fb / fa;
            let (p, mut q) = if a == c {
                let p = xm * s;
                (p + p, 1.0 - s)
            } else {
                let q = fa / fc;
                let r = fb / fc;
                (
                    s * ((xm + xm) * q * (q - r) - (b - a) * (r - 1.0)),
                    (q - 1.0) * (r - 1.0) * (s - 1.0),
                )
            };
            if p > 0.0 {
                q = -q;
            }
            let p = p.abs();
            let min1 = 3.0 * xm * q - (tol1 * q).abs();
            let min2 = (e * q).abs();
            if (p + p) < if min1 < min2 { min1 } else { min2 } {
                e = d;
                d = p / q;
            } else {
                d = xm;
                e = d;
            }
        } else {
            d = xm;
            e = d;
        }
        a = b;
        fa = fb;
        if d.abs() > tol1 {
            b += d;
        } else if xm >= 0.0 {
            b += tol1.abs();
        } else {
            b += -tol1.abs();
        }
        f.value(b, &mut fb);
        fb -= k;
    }
}
