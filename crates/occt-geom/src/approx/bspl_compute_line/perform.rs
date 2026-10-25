//! `Approx_BSplComputeLine::Perform` and `FindRealConstraints`
//! (`Approx_BSplComputeLine.gxx:636-788`, `:1206-1281`).

use occt_math::MathVector;

use crate::approx::line_tool::LineTool;
use crate::app_par_curves::constraint::{Constraint, ConstraintCouple};
use crate::approx::bspl_compute_line::BsplComputeLine;

/// `Cons >= AppParCurves_TangencyPoint` (`AppParCurves_Constraint.hxx:22-30`).
pub(super) fn at_least_tangency(c: Constraint) -> bool {
    matches!(c, Constraint::TangencyPoint | Constraint::CurvaturePoint)
}

/// `Perform(Line)` (`gxx:636-788`).
pub(super) fn perform<L: LineTool>(this: &mut BsplComputeLine, line: &L) {
    find_real_constraints(this, line);

    let thefirstpt = line.first_point();
    let thelastpt = line.last_point();
    let myfirstpt = thefirstpt;
    let mylastpt = thelastpt;

    // `gxx:656-658`.
    this.my_constraints = vec![
        ConstraintCouple::from_parts(myfirstpt, this.real_first_c),
        ConstraintCouple::from_parts(mylastpt, this.real_last_c),
    ];

    // `gxx:660-672`.
    let mut the_param = MathVector::with_init(
        thefirstpt.max(0) as usize,
        thelastpt.max(0) as usize,
        0.0,
    );
    if this.my_first_param.is_none() {
        let slice = super::parameters::compute_parameters(this, line, thefirstpt, thelastpt);
        for i in thefirstpt..=thelastpt {
            the_param.set_value(i as usize, slice[(i - thefirstpt) as usize]);
        }
    } else {
        let (lower, vals) = this.my_first_param.clone().unwrap();
        for i in lower..=lower + vals.len() as i32 - 1 {
            the_param.set_value(
                (i + thefirstpt - 1) as usize,
                vals[(i - lower) as usize],
            );
        }
    }
    this.my_parameters = Some((
        the_param.lower() as i32,
        (the_param.lower()..=the_param.upper())
            .map(|i| the_param.value(i))
            .collect(),
    ));

    let mut nbknots = 2;
    this.alldone = false;

    if !this.my_cut {
        // Case where no additional knots are desired (`gxx:684-708`).
        if !this.my_has_knots {
            // `theknots(1) = 0.; theknots(2) = 1.;` (`gxx:690-694`); `themults`
            // is left unwritten by OCCT and fully overwritten in `Compute`
            // because `myhasmults` is false.
            let theknots = vec![0.0, 1.0];
            let mut themults = vec![0i32; 2];
            this.alldone = super::compute::compute(
                this,
                line,
                myfirstpt,
                mylastpt,
                &mut the_param,
                &theknots,
                &mut themults,
            );
        } else if !this.my_has_mults {
            let (_, knots) = this.my_knots.clone().unwrap();
            let mut themults = vec![0i32; knots.len()];
            this.alldone = super::compute::compute(
                this,
                line,
                myfirstpt,
                mylastpt,
                &mut the_param,
                &knots,
                &mut themults,
            );
        } else {
            let (_, knots) = this.my_knots.clone().unwrap();
            let (_, mut mults) = this.my_mults.clone().unwrap();
            this.alldone = super::compute::compute(
                this,
                line,
                myfirstpt,
                mylastpt,
                &mut the_param,
                &knots,
                &mut mults,
            );
        }
    } else {
        // Iterate from user knots or from a Bezier (`gxx:710-782`).
        let mut finish = false;
        let mut begin = true;
        while !finish {
            this.current_tol3d = super::REAL_LAST;
            this.current_tol2d = super::REAL_LAST;

            if this.my_has_knots && begin {
                let (_, knots) = this.my_knots.clone().unwrap();
                if !this.my_has_mults {
                    // The user gives starting knots, the multiplicities follow
                    // from the continuity (`gxx:724-737`).
                    let mut the_mults = vec![0i32; knots.len()];
                    this.alldone = super::compute::compute(
                        this,
                        line,
                        myfirstpt,
                        mylastpt,
                        &mut the_param,
                        &knots,
                        &mut the_mults,
                    );
                } else {
                    // The user gives knots and multiplicities (`gxx:739-748`).
                    let (_, mut mults) = this.my_mults.clone().unwrap();
                    this.alldone = super::compute::compute(
                        this,
                        line,
                        myfirstpt,
                        mylastpt,
                        &mut the_param,
                        &knots,
                        &mut mults,
                    );
                }
                begin = false;
            } else {
                // No starting knots: spread `nbknots` knots over the
                // parametrisation (`gxx:752-777`).
                let mut theknots = vec![0.0f64; nbknots as usize];
                let mut themults = vec![0i32; nbknots as usize];
                theknots[0] = 0.0;
                theknots[nbknots as usize - 1] = 1.0;
                for i in 2..=nbknots - 1 {
                    let l = (mylastpt - myfirstpt) as f64 * (i - 1) as f64
                        / (nbknots - 1) as f64;
                    let ll = l as i32;
                    let a = l - ll as f64;
                    let p1 = the_param.value((ll + myfirstpt) as usize);
                    let p2 = the_param.value((ll + 1 + myfirstpt) as usize);
                    theknots[i as usize - 1] = (1.0 - a) * p1 + a * p2;
                }

                this.alldone = super::compute::compute(
                    this,
                    line,
                    myfirstpt,
                    mylastpt,
                    &mut the_param,
                    &theknots,
                    &mut themults,
                );
            }

            if !this.alldone {
                nbknots += 1;
            } else {
                finish = true;
            }
        }
    }
}

/// `FindRealConstraints(Line)` (`gxx:1206-1281`).
fn find_real_constraints<L: LineTool>(this: &mut BsplComputeLine, line: &L) {
    this.real_first_c = this.my_first_c;
    this.real_last_c = this.my_last_c;
    let nb_p3d = line.nb_p3d();
    let nb_p2d = line.nb_p2d();
    let thefirstpt = line.first_point();
    let thelastpt = line.last_point();

    if at_least_tangency(this.my_first_c) {
        let ok = tangency_available(line, thefirstpt, nb_p3d, nb_p2d);

        this.real_first_c = Constraint::PassPoint;
        if ok {
            this.real_first_c = Constraint::TangencyPoint;
            if this.my_first_c == Constraint::CurvaturePoint {
                // `gxx:1236-1248` repeats the same `Tangency` call.
                let ok = tangency_available(line, thefirstpt, nb_p3d, nb_p2d);
                if ok {
                    this.real_first_c = Constraint::CurvaturePoint;
                }
            }
        }
    }

    if at_least_tangency(this.my_last_c) {
        let ok = tangency_available(line, thelastpt, nb_p3d, nb_p2d);

        this.real_last_c = Constraint::PassPoint;
        if ok {
            this.real_last_c = Constraint::TangencyPoint;
            if this.my_last_c == Constraint::CurvaturePoint {
                // `gxx:1265-1277` repeats the same `Tangency` call.
                let ok = tangency_available(line, thelastpt, nb_p3d, nb_p2d);
                if ok {
                    this.real_last_c = Constraint::CurvaturePoint;
                }
            }
        }
    }
}

/// The `LineTool::Tangency` overload selected by the point dimensions
/// (`gxx:1224-1229`, `:1254-1259`).
fn tangency_available<L: LineTool>(line: &L, index: i32, nb_p3d: i32, nb_p2d: i32) -> bool {
    if nb_p3d != 0 && nb_p2d != 0 {
        line.tangency(index).is_some()
    } else if nb_p2d != 0 {
        line.tangency_2d(index).is_some()
    } else if nb_p3d != 0 {
        line.tangency_3d(index).is_some()
    } else {
        false
    }
}




