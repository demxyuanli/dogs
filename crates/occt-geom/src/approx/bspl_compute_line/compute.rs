//! `Approx_BSplComputeLine::Compute` (`Approx_BSplComputeLine.gxx:896-1080`):
//! one degree of the degree-escalating loop.

use crate::approx::line_tool::LineTool;
use crate::app_par_curves::bsp_gradient::BspGradient;
use crate::app_par_curves::constraint::Constraint;
use crate::app_par_curves::gradient::Gradient;
use crate::app_par_curves::least_square::LeastSquare;
use crate::app_par_curves::multi_bspline_curve::MultiBSpCurve;
use crate::approx::bspl_compute_line::BsplComputeLine;
use occt_math::MathVector;

/// `Cons >= AppParCurves_TangencyPoint` (`AppParCurves_Constraint.hxx:22-30`).
fn at_least_tangency(c: Constraint) -> bool {
    matches!(c, Constraint::TangencyPoint | Constraint::CurvaturePoint)
}

/// What one `deg` iteration does next: `return true` or the `continue` of the
/// degree loop.
enum DegOutcome {
    Done,
    Next,
}

/// `Compute(Line, fpt, lpt, Para, Knots, Mults)` (`gxx:896-1080`).
///
/// `Mults` is an in-out array in OCCT (the `!myhasmults` arm rewrites it from
/// the degree and the continuity), so the port takes `&mut [i32]`.
pub(super) fn compute<L: LineTool>(
    this: &mut BsplComputeLine,
    line: &L,
    fpt: i32,
    lpt: i32,
    para: &mut MathVector,
    knots: &[f64],
    mults: &mut [i32],
) -> bool {
    let nbp = lpt - fpt + 1;
    this.my_lambda1 = 0.0;
    this.my_lambda2 = 0.0;

    let mut a_params = MathVector::with_init(para.lower(), para.upper(), 0.0);

    for deg in this.my_degremin..=this.my_degremax {
        // `aParams = Para` (`gxx:917`).
        for i in para.lower()..=para.upper() {
            a_params.set_value(i, para.value(i));
        }

        // `gxx:919-942`: the pole count follows the multiplicities, which the
        // `!myhasmults` arm derives from the degree and the continuity.
        let nbpoles = if !this.my_has_mults {
            let last = mults.len() - 1;
            mults[0] = deg + 1;
            mults[last] = deg + 1;
            let multinter = if this.my_cont == -1 {
                1
            } else {
                1.max(deg - this.my_cont)
            };
            let mut n = deg + 1;
            for m in mults.iter_mut().take(last).skip(1) {
                *m = multinter;
                n += multinter;
            }
            n
        } else {
            let mut n = -deg - 1;
            for m in mults.iter() {
                n += *m;
            }
            n
        };

        match compute_at_degree(this, line, fpt, lpt, &mut a_params, knots, mults, deg, nbpoles, nbp)
        {
            DegOutcome::Done => return true,
            DegOutcome::Next => {}
        }
    }

    false
}

/// The body of one `deg` iteration (`gxx:943-1079`); `nbpoles` was established
/// by the caller.
#[allow(clippy::too_many_arguments)]
fn compute_at_degree<L: LineTool>(
    this: &mut BsplComputeLine,
    line: &L,
    fpt: i32,
    lpt: i32,
    a_params: &mut MathVector,
    knots: &[f64],
    mults: &[i32],
    deg: i32,
    nbpoles: i32,
    nbp: i32,
) -> DegOutcome {
    // `gxx:943-952`: four independent tests. `CurvaturePoint` matches only the
    // third/fourth one, so a constrained end adds exactly one either way.
    let mut nbpolestocompare = nbpoles;
    if this.real_first_c == Constraint::TangencyPoint {
        nbpolestocompare += 1;
    }
    if this.real_last_c == Constraint::TangencyPoint {
        nbpolestocompare += 1;
    }
    if this.real_first_c == Constraint::CurvaturePoint {
        nbpolestocompare += 1;
    }
    if this.real_last_c == Constraint::CurvaturePoint {
        nbpolestocompare += 1;
    }
    if nbpolestocompare > nbp {
        // `gxx:954-958`: not enough points, fall back to interpolation.
        super::interpol::interpol(this, line);
        this.tolreached = true;
        return DegOutcome::Done;
    }

    let mut my_scu = MultiBSpCurve::with_nb_poles(nbpoles.max(0) as usize);
    let mut the_tol3d = 0.0f64;
    let mut the_tol2d = 0.0f64;
    let mydone;

    if this.my_squares {
        let mut sq = LeastSquare::with_knots_parameters(
            line,
            knots,
            mults,
            fpt,
            lpt,
            this.real_first_c,
            this.real_last_c,
            a_params,
            nbpoles,
        );
        mydone = sq.is_done();
        if !mydone {
            return DegOutcome::Next;
        }
        my_scu = sq.bspline_value().clone();
        let (_fv, t3, t2) = sq.error();
        the_tol3d = t3;
        the_tol2d = t2;
    } else if nbpoles != deg + 1 {
        // `gxx:977-1014`: the tangency multipliers are searched only at the
        // lowest degree, then divided back out for the shared constructor.
        if deg == this.my_degremin
            && (at_least_tangency(this.real_first_c)
                || at_least_tangency(this.real_last_c))
        {
            let mut thefitt = LeastSquare::with_knots_parameters(
                line,
                knots,
                mults,
                fpt,
                lpt,
                this.real_first_c,
                this.real_last_c,
                a_params,
                nbpoles,
            );
            this.my_lambda1 = thefitt.first_lambda() * deg as f64;
            this.my_lambda2 = thefitt.last_lambda() * deg as f64;
        }
        let l1 = this.my_lambda1 / deg as f64;
        let l2 = this.my_lambda2 / deg as f64;
        let constraints = this.my_constraints.clone();

        let mut grad = BspGradient::with_lambdas(
            line,
            fpt,
            lpt,
            &constraints,
            a_params,
            knots,
            mults,
            deg,
            this.my_tol3d,
            this.my_tol2d,
            this.my_iter_max,
            l1,
            l2,
        );
        mydone = grad.is_done();
        if !mydone {
            return DegOutcome::Next;
        }
        my_scu = grad.value();
        the_tol3d = grad.max_error3d();
        the_tol2d = grad.max_error2d();
    } else {
        // `gxx:1015-1030`: as many poles as the degree, so the Bezier driver.
        let constraints = this.my_constraints.clone();
        let mut grad2 = Gradient::new(
            line,
            fpt,
            lpt,
            &constraints,
            a_params,
            deg,
            this.my_tol3d,
            this.my_tol2d,
            this.my_iter_max,
        );
        mydone = grad2.is_done();
        if !mydone {
            return DegOutcome::Next;
        }
        if grad2.value().nb_curves() == 0 {
            return DegOutcome::Next;
        }
        my_scu = MultiBSpCurve::from_curve(&grad2.value(), knots, mults);
        the_tol3d = grad2.max_error3d();
        the_tol2d = grad2.max_error2d();
    }

    // `gxx:1031-1040`: parameters outside `[0, 1]` are not stored.
    let mut save = true;
    for i in a_params.lower()..=a_params.upper() {
        let v = a_params.value(i);
        if v <= -0.000001 || v >= 1.000001 {
            save = false;
            break;
        }
    }

    if mydone && the_tol3d <= this.my_tol3d && the_tol2d <= this.my_tol2d {
        // Stockage de la multicurve approximee (`gxx:1044-1060`).
        this.tolreached = true;
        this.the_multi_bsp_curve = my_scu;
        this.current_tol3d = the_tol3d;
        this.current_tol2d = the_tol2d;
        if save {
            store_parameters(this, a_params);
        }
        return DegOutcome::Done;
    }

    // `gxx:1063-1076`: keep the best result so far and go on.
    if the_tol3d <= this.current_tol3d && the_tol2d <= this.current_tol2d {
        this.the_multi_bsp_curve = my_scu;
        this.current_tol3d = the_tol3d;
        this.current_tol2d = the_tol2d;
        if save {
            store_parameters(this, a_params);
        }
    }

    DegOutcome::Next
}

/// `myParameters->SetValue(i, aParams(i))` (`gxx:1054-1057`, `:1068-1074`).
fn store_parameters(this: &mut BsplComputeLine, a_params: &MathVector) {
    let Some((lower, vals)) = this.my_parameters.as_mut() else {
        return;
    };
    for i in a_params.lower()..=a_params.upper() {
        let slot = i as i32 - *lower;
        if slot >= 0 && (slot as usize) < vals.len() {
            vals[slot as usize] = a_params.value(i);
        }
    }
}





