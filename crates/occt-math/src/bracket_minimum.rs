//! `math_BracketMinimum` — bracket the minimum of a 1-D function.
//!
//! Source: `math_BracketMinimum.cxx` (`Perform` `:62-211`, `LimitAndMayBeSwap`
//! `:34-60`), `math_BracketMinimum.hxx:37-136` and `math_BracketMinimum.lxx`
//! (`SetLimits` `:46-51`, `SetFA`/`SetFB` `:53-63`, `Limited` `:65-68`).
//!
//! The port keeps OCCT's field names (`Ax`/`Bx`/`Cx`, `FAx`/`FBx`/`FCx`,
//! `myLeft`/`myRight`, `myIsLimited`, `myFA`/`myFB`, `Done`) so the control flow
//! below can be read side by side with the `.cxx`.

/// `GOLD` (`math_BracketMinimum.cxx:20`).
const GOLD: f64 = 1.618_034;
/// `GLIMIT` (`math_BracketMinimum.cxx:22`).
const GLIMIT: f64 = 100.0;
/// `TINY` (`math_BracketMinimum.cxx:23`).
const TINY: f64 = 1.0e-20;
/// `Precision::PConfusion()` (`Precision.hxx:334`).
const P_CONFUSION: f64 = 1.0e-9;
/// `Precision::Infinite()` (`Precision.hxx:371`).
const INFINITE: f64 = 2.0e100;

/// `math_BracketMinimum` (`math_BracketMinimum.hxx:37-136`).
#[derive(Debug, Clone)]
pub struct BracketMinimum {
    /// `Done` (`math_BracketMinimum.hxx:124`).
    done: bool,
    /// `Ax` / `Bx` / `Cx` (`math_BracketMinimum.hxx:125-127`).
    ax: f64,
    bx: f64,
    cx: f64,
    /// `FAx` / `FBx` / `FCx` (`math_BracketMinimum.hxx:128-130`).
    fax: f64,
    fbx: f64,
    fcx: f64,
    /// `myLeft` / `myRight` (`math_BracketMinimum.hxx:131-132`).
    left: f64,
    right: f64,
    /// `myIsLimited` (`math_BracketMinimum.hxx:133`).
    is_limited: bool,
    /// `myFA` / `myFB` (`math_BracketMinimum.hxx:134-135`).
    my_fa: bool,
    my_fb: bool,
}

impl BracketMinimum {
    /// `math_BracketMinimum(const double A, const double B)`
    /// (`math_BracketMinimum.lxx:19-33`): stores A and B only, does **not**
    /// perform the job.
    pub fn new(a: f64, b: f64) -> Self {
        Self {
            done: false,
            ax: a,
            bx: b,
            cx: 0.0,
            fax: 0.0,
            fbx: 0.0,
            fcx: 0.0,
            left: -INFINITE,
            right: INFINITE,
            is_limited: false,
            my_fa: false,
            my_fb: false,
        }
    }

    /// `SetLimits(theLeft, theRight)` (`math_BracketMinimum.lxx:46-51`).
    pub fn set_limits(&mut self, the_left: f64, the_right: f64) {
        self.left = the_left;
        self.right = the_right;
        self.is_limited = true;
    }

    /// `SetFA(theValue)` (`math_BracketMinimum.lxx:53-57`).
    pub fn set_fa(&mut self, the_value: f64) {
        self.fax = the_value;
        self.my_fa = true;
    }

    /// `SetFB(theValue)` (`math_BracketMinimum.lxx:59-63`).
    pub fn set_fb(&mut self, the_value: f64) {
        self.fbx = the_value;
        self.my_fb = true;
    }

    /// `IsDone()` (`math_BracketMinimum.lxx:35-38`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `Values(A, B, C)` (`math_BracketMinimum.cxx:271-278`); `None` mirrors
    /// `StdFail_NotDone` when the bracket was not found.
    pub fn values(&self) -> Option<(f64, f64, f64)> {
        if self.done {
            Some((self.ax, self.bx, self.cx))
        } else {
            None
        }
    }

    /// `FunctionValues(FA, FB, FC)` (`math_BracketMinimum.cxx:280-287`); `None`
    /// mirrors `StdFail_NotDone`.
    pub fn function_values(&self) -> Option<(f64, f64, f64)> {
        if self.done {
            Some((self.fax, self.fbx, self.fcx))
        } else {
            None
        }
    }

    /// `Limited(theValue)` (`math_BracketMinimum.lxx:65-68`).
    fn limited(&self, the_value: f64) -> f64 {
        if the_value < self.left {
            self.left
        } else if the_value > self.right {
            self.right
        } else {
            the_value
        }
    }

    /// Write the working triplet back into the object's members.
    fn commit(&mut self, ax: f64, bx: f64, cx: f64, fax: f64, fbx: f64, fcx: f64) {
        self.ax = ax;
        self.bx = bx;
        self.cx = cx;
        self.fax = fax;
        self.fbx = fbx;
        self.fcx = fcx;
    }

    /// `math_BracketMinimum::Perform(F)` (`math_BracketMinimum.cxx:62-211`).
    ///
    /// The objective is infallible here (the port's closures return a value),
    /// so the `F.Value == false` early returns of the `.cxx` collapse to the
    /// single `LimitAndMayBeSwap` failure (`cxx:42-45`). The triplet is carried
    /// in locals and committed only when `Done` is set, so the early returns
    /// leave `done == false` exactly as the `.cxx` does.
    pub fn perform<F: Fn(f64) -> f64>(&mut self, f: &F) {
        self.done = false;
        let lambda = GOLD;

        let (mut ax, mut bx) = (self.ax, self.bx);
        let (mut fax, mut fbx) = (self.fax, self.fbx);
        let mut cx: f64;
        let mut fcx: f64;

        // cxx:70-85
        if !self.my_fa {
            fax = f(ax);
        }
        if !self.my_fb {
            fbx = f(bx);
        }
        if fbx > fax {
            std::mem::swap(&mut ax, &mut bx);
            std::mem::swap(&mut fax, &mut fbx);
        }

        // get next prob after (A, B) — cxx:92-109
        cx = bx + lambda * (bx - ax);
        if self.is_limited {
            match limit_and_may_be_swap(f, ax, bx, fbx, cx, self.left, self.right) {
                Some((b_new, c_new, fb_new, fc_new)) => {
                    bx = b_new;
                    cx = c_new;
                    fbx = fb_new;
                    fcx = fc_new;
                }
                None => return,
            }
        } else {
            fcx = f(cx);
        }

        // cxx:111-209
        while fbx > fcx {
            let r = (bx - ax) * (fbx - fcx);
            let q = (bx - cx) * (fbx - fax);
            let mut u = bx
                - ((bx - cx) * q - (bx - ax) * r)
                    / (2.0 * sign((q - r).abs().max(TINY), q - r));
            let mut ulim = bx + GLIMIT * (cx - bx);
            if self.is_limited {
                ulim = self.limited(ulim);
            }

            let mut fu;
            if (bx - u) * (u - cx) > 0.0 {
                // u is between B and C
                fu = f(u);
                if fu < fcx {
                    // solution is found (B, u, c)
                    ax = bx;
                    bx = u;
                    fax = fbx;
                    fbx = fu;
                    self.commit(ax, bx, cx, fax, fbx, fcx);
                    self.done = true;
                    return;
                } else if fu > fbx {
                    // solution is found (A, B, u)
                    cx = u;
                    fcx = fu;
                    self.commit(ax, bx, cx, fax, fbx, fcx);
                    self.done = true;
                    return;
                }
                // get next prob after (B, C)
                u = cx + lambda * (cx - bx);
                if self.is_limited {
                    match limit_and_may_be_swap(f, bx, cx, fcx, u, self.left, self.right) {
                        Some((b_new, c_new, fb_new, fc_new)) => {
                            // SHFT(Ax, Bx, Cx, u); SHFT(FAx, FBx, FCx, fu);
                            // (the call may have swapped the old Cx and u).
                            ax = bx;
                            bx = b_new;
                            cx = c_new;
                            fax = fbx;
                            fbx = fb_new;
                            fcx = fc_new;
                            continue;
                        }
                        None => return,
                    }
                } else {
                    fu = f(u);
                }
            } else if (cx - u) * (u - ulim) > 0.0 {
                // u is beyond C but between C and limit
                fu = f(u);
            } else if (u - ulim) * (ulim - cx) >= 0.0 {
                // u is beyond limit
                u = ulim;
                fu = f(u);
            } else {
                // u tends to approach to the side of A, reset it to the next
                // prob after (B, C)
                u = cx + GOLD * (cx - bx);
                if self.is_limited {
                    match limit_and_may_be_swap(f, bx, cx, fcx, u, self.left, self.right) {
                        Some((b_new, c_new, fb_new, fc_new)) => {
                            ax = bx;
                            bx = b_new;
                            cx = c_new;
                            fax = fbx;
                            fbx = fb_new;
                            fcx = fc_new;
                            continue;
                        }
                        None => return,
                    }
                } else {
                    fu = f(u);
                }
            }
            // SHFT(Ax, Bx, Cx, u); SHFT(FAx, FBx, FCx, fu); — cxx:207-208
            ax = bx;
            bx = cx;
            cx = u;
            fax = fbx;
            fbx = fcx;
            fcx = fu;
        }
        self.commit(ax, bx, cx, fax, fbx, fcx);
        self.done = true;
    }
}

/// `math_BracketMinimum::LimitAndMayBeSwap` (`math_BracketMinimum.cxx:34-60`).
///
/// Returns the (possibly swapped) `(B, C, FB, FC)`; `None` when `C` collapsed
/// onto `B` within `Precision::PConfusion()` (`cxx:42-45`) or when the function
/// evaluation failed (the port's objective is infallible, so only the former).
fn limit_and_may_be_swap<F: Fn(f64) -> f64>(
    f: &F,
    the_a: f64,
    the_b: f64,
    the_fb: f64,
    the_c: f64,
    the_left: f64,
    the_right: f64,
) -> Option<(f64, f64, f64, f64)> {
    let the_c = if the_c < the_left {
        the_left
    } else if the_c > the_right {
        the_right
    } else {
        the_c
    };
    if (the_b - the_c).abs() < P_CONFUSION {
        return None;
    }
    let the_fc = f(the_c);
    // check that B is between A and C
    if (the_a - the_b) * (the_b - the_c) < 0.0 {
        // swap B and C
        Some((the_c, the_b, the_fc, the_fb))
    } else {
        Some((the_b, the_c, the_fb, the_fc))
    }
}

/// `SIGN(a, b)` (`math_BracketMinimum.cxx:28`): magnitude of `a`, sign of `b`.
fn sign(a: f64, b: f64) -> f64 {
    if b > 0.0 {
        a.abs()
    } else {
        -a.abs()
    }
}
