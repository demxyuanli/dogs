//! `math_BrentMinimum` — Brent's method for the minimum of a function of one
//! variable, no derivative required.
//!
//! Source: `math_BrentMinimum.cxx` (`Perform` `:80-205`), `math_BrentMinimum.lxx`
//! (`IsSolutionReached`, accessors) and `math_BrentMinimum.hxx` (defaults
//! `NbIterations = 100`, `ZEPS = 1.0e-12`).
//!
//! The port keeps OCCT's field names (`a`, `b`, `x`, `fx`, `fv`, `fw`, `XTol`,
//! `EPSZ`) and its `Done`/`iter`/`Itermax`/`myF` state, so the control flow below
//! can be read side by side with the `.cxx`.

/// `CGOLD` (`math_BrentMinimum.cxx:21`): `0.5 * (3 - sqrt(5))`.
const CGOLD: f64 = 0.3819660;

/// `math_BrentMinimum` (`math_BrentMinimum.hxx:32-127`).
#[derive(Debug, Clone)]
pub struct BrentMinimum {
    a: f64,
    b: f64,
    x: f64,
    fx: f64,
    fv: f64,
    fw: f64,
    xtol: f64,
    epsz: f64,
    done: bool,
    iter: i32,
    itermax: i32,
    /// `myF`: true when `F(Bx)` was supplied to the constructor, so `Perform`
    /// must not evaluate it again (`math_BrentMinimum.cxx:92-100`).
    my_f: bool,
}

impl BrentMinimum {
    /// `math_BrentMinimum(TolX, NbIterations = 100, ZEPS = 1.0e-12)`
    /// (`math_BrentMinimum.cxx:38-54`).
    pub fn new(tol_x: f64, nb_iterations: i32, zeps: f64) -> Self {
        Self {
            a: 0.0,
            b: 0.0,
            x: 0.0,
            fx: 0.0,
            fv: 0.0,
            fw: 0.0,
            xtol: tol_x,
            epsz: zeps,
            done: false,
            iter: 0,
            itermax: nb_iterations,
            my_f: false,
        }
    }

    /// Defaults of `math_BrentMinimum.hxx:38-40`.
    pub fn with_tolerance(tol_x: f64) -> Self {
        Self::new(tol_x, 100, 1.0e-12)
    }

    /// `math_BrentMinimum(TolX, Fbx, NbIterations = 100, ZEPS = 1.0e-12)`
    /// (`math_BrentMinimum.cxx:58-75`): `F(Bx)` is already known.
    pub fn with_fbx(tol_x: f64, fbx: f64, nb_iterations: i32, zeps: f64) -> Self {
        Self {
            a: 0.0,
            b: 0.0,
            x: 0.0,
            fx: fbx,
            fv: 0.0,
            fw: 0.0,
            xtol: tol_x,
            epsz: zeps,
            done: false,
            iter: 0,
            itermax: nb_iterations,
            my_f: true,
        }
    }

    /// `IsSolutionReached` (`math_BrentMinimum.lxx:19-23`).
    fn is_solution_reached(&self) -> bool {
        let two_tol = 2.0 * (self.xtol * self.x.abs() + self.epsz);
        self.x <= two_tol + self.a && self.x >= self.b - two_tol
    }

    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `Location()` (`math_BrentMinimum.lxx:40-44`); `None` mirrors
    /// `StdFail_NotDone` when the minimization did not converge.
    pub fn location(&self) -> Option<f64> {
        if self.done {
            Some(self.x)
        } else {
            None
        }
    }

    /// `Minimum()` (`math_BrentMinimum.lxx:46-50`).
    pub fn minimum(&self) -> Option<f64> {
        if self.done {
            Some(self.fx)
        } else {
            None
        }
    }

    /// `NbIterations()` (`math_BrentMinimum.lxx:52-56`).
    pub fn nb_iterations(&self) -> Option<i32> {
        if self.done {
            Some(self.iter)
        } else {
            None
        }
    }

    /// `math_BrentMinimum::Perform(F, Ax, Bx, Cx)` (`math_BrentMinimum.cxx:80-205`).
    ///
    /// Brent minimization on the bracketing triplet `(ax, bx, cx)` with `bx`
    /// between `ax` and `cx` and `F(bx)` below both. Returns `false` where OCCT
    /// leaves `Done = false` (a failed function evaluation, or `Itermax`
    /// iterations without reaching `IsSolutionReached`).
    pub fn perform<F: Fn(f64) -> f64>(&mut self, f: &F, ax: f64, bx: f64, cx: f64) -> bool {
        let mut e = 0.0f64;
        let mut d = f64::MAX; // RealLast()
        let mut etemp;
        let mut p;
        let mut q;
        let mut r;
        let mut u;

        self.a = if ax < cx { ax } else { cx };
        self.b = if ax > cx { ax } else { cx };
        // `x = w = v = bx;` — OCCT's `x`/`w`/`v` are members, so
        // `IsSolutionReached` always reads the current `x`.
        self.x = bx;
        let mut w = bx;
        let mut v = bx;
        if !self.my_f {
            self.fx = f(self.x);
        }
        self.fw = self.fx;
        self.fv = self.fx;

        for iter in 1..=self.itermax {
            let xm = 0.5 * (self.a + self.b);
            let tol1 = self.xtol * self.x.abs() + self.epsz;
            let tol2 = 2.0 * tol1;
            if self.is_solution_reached() {
                self.iter = iter;
                self.done = true;
                return true;
            }
            if e.abs() > tol1 {
                r = (self.x - w) * (self.fx - self.fv);
                q = (self.x - v) * (self.fx - self.fw);
                p = (self.x - v) * q - (self.x - w) * r;
                q = 2.0 * (q - r);
                if q > 0.0 {
                    p = -p;
                }
                q = q.abs();
                etemp = e;
                e = d;
                if p.abs() >= (0.5 * q * etemp).abs()
                    || p <= q * (self.a - self.x)
                    || p >= q * (self.b - self.x)
                {
                    e = if self.x >= xm {
                        self.a - self.x
                    } else {
                        self.b - self.x
                    };
                    d = CGOLD * e;
                } else {
                    d = p / q;
                    u = self.x + d;
                    if u - self.a < tol2 || self.b - u < tol2 {
                        d = copysign(tol1, xm - self.x);
                    }
                }
            } else {
                e = if self.x >= xm {
                    self.a - self.x
                } else {
                    self.b - self.x
                };
                d = CGOLD * e;
            }
            u = if d.abs() >= tol1 {
                self.x + d
            } else {
                self.x + copysign(tol1, d)
            };
            let fu = f(u);
            if fu <= self.fx {
                if u >= self.x {
                    self.a = self.x;
                } else {
                    self.b = self.x;
                }
                // SHFT(v, w, x, u); SHFT(fv, fw, fx, fu)
                v = w;
                w = self.x;
                self.x = u;
                self.fv = self.fw;
                self.fw = self.fx;
                self.fx = fu;
            } else {
                if u < self.x {
                    self.a = u;
                } else {
                    self.b = u;
                }
                if fu <= self.fw || w == self.x {
                    v = w;
                    w = u;
                    self.fv = self.fw;
                    self.fw = fu;
                } else if fu <= self.fv || v == self.x || v == w {
                    v = u;
                    self.fv = fu;
                }
            }
        }
        self.done = false;
        false
    }
}

/// `copysign` with OCCT's `std::copysign` semantics (magnitude of `mag`, sign of
/// `sign`).
fn copysign(mag: f64, sign: f64) -> f64 {
    if sign < 0.0 {
        -mag.abs()
    } else {
        mag.abs()
    }
}
