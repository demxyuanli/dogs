//! `Approx_BSplComputeLine::Interpol` (`Approx_BSplComputeLine.gxx:1284-1422`):
//! the constrained interpolation fallback, also used when the degree loop runs
//! out of comparison poles (`gxx:954-958`).

use occt_math::MathVector;

use crate::approx::line_tool::LineTool;
use crate::app_par_curves::constraint::Constraint;
use crate::app_par_curves::least_square::LeastSquare;
use crate::app_par_curves::multi_bspline_curve::MultiBSpCurve;
use crate::approx::bspl_compute_line::BsplComputeLine;

/// `Interpol(Line)` (`gxx:1284-1422`).
pub(super) fn interpol<L: LineTool>(this: &mut BsplComputeLine, line: &L) {
    // `gxx:1285-1288`: `deg` is the local constant 3, continuity is fixed to 2.
    let deg = 3;
    this.my_cont = 2;
    let thefirstpt = line.first_point();
    let thelastpt = line.last_point();

    let mut the_param = MathVector::with_init(
        thefirstpt.max(0) as usize,
        thelastpt.max(0) as usize,
        0.0,
    );
    if this.my_first_param.is_none() {
        // `gxx:1292-1295`.
        let slice = super::parameters::compute_parameters(this, line, thefirstpt, thelastpt);
        for i in thefirstpt..=thelastpt {
            the_param.set_value(i as usize, slice[(i - thefirstpt) as usize]);
        }
    } else {
        // `gxx:1296-1302`: `TheParam(i + Thefirstpt - 1) = myfirstParam(i)`.
        let (lower, vals) = this.my_first_param.clone().unwrap();
        for i in lower..=lower + vals.len() as i32 - 1 {
            the_param.set_value((i + thefirstpt - 1) as usize, vals[(i - lower) as usize]);
        }
    }

    let mut cons = Constraint::TangencyPoint;
    let mut lambda1;
    let mut lambda2;

    let nbpoints = thelastpt - thefirstpt + 1;

    if nbpoints == 2 {
        // `gxx:1311-1326`: a straight segment is interpolated, not approximated.
        cons = Constraint::NoConstraint;
        let mydeg = 1;
        let mut lsq = LeastSquare::with_parameters(
            line,
            thefirstpt,
            thelastpt,
            cons,
            cons,
            &the_param,
            mydeg + 1,
        );
        this.alldone = lsq.is_done();
        let the_knots = vec![
            the_param.value(thefirstpt as usize),
            the_param.value(thelastpt as usize),
        ];
        let the_mults = vec![mydeg + 1, mydeg + 1];
        this.the_multi_bsp_curve =
            MultiBSpCurve::from_curve(&lsq.bezier_value(), &the_knots, &the_mults);
        let (_fv, t3, t2) = lsq.error();
        this.current_tol3d = t3;
        this.current_tol2d = t2;
    } else {
        // `gxx:1328-1331`: one interior knot per point, two extra poles.
        let nbpoles = nbpoints + 2;
        let nbknots = nbpoints;

        // `gxx:1332-1346`.
        let mut theknots = vec![0.0f64; nbknots as usize];
        theknots[0] = the_param.value(thefirstpt as usize);
        theknots[nbknots as usize - 1] = the_param.value(thelastpt as usize);
        let mut the_mults = vec![0i32; nbknots as usize];
        the_mults[0] = deg + 1;
        the_mults[nbknots as usize - 1] = deg + 1;
        let low = the_param.lower() as i32;
        for i in 2..=nbknots - 1 {
            theknots[i as usize - 1] = the_param.value((i + low - 1) as usize);
            the_mults[i as usize - 1] = 1;
        }

        let nb_p = (3 * line.nb_p3d() + 2 * line.nb_p2d()).max(1) as usize;
        let mut v1 = MathVector::with_init(1, nb_p, 0.0);
        let mut v2 = MathVector::with_init(1, nb_p, 0.0);

        if nbpoints == 3 || nbpoints == 4 {
            // `gxx:1349-1359`.
            v1 = super::tangency::first_tangency_vector(this, line, thefirstpt);
            lambda1 = super::tangency::search_first_lambda(
                line,
                &the_param,
                &theknots,
                &v1,
                thefirstpt,
            );
            v2 = super::tangency::last_tangency_vector(this, line, thelastpt);
            lambda2 = super::tangency::search_last_lambda(
                line,
                &the_param,
                &theknots,
                &v2,
                thelastpt,
            );
            lambda1 /= deg as f64;
            lambda2 /= deg as f64;
        } else {
            // `gxx:1360-1402`: fit the first/last `nnp` points with a Bezier
            // and take the tangent there.
            let nnp = nbpoints.min(9);
            let nnpol = nnp;
            let lastp = thelastpt.min(thefirstpt + nnp - 1);
            let mut sq1 = LeastSquare::new(line, thefirstpt, lastp, cons, cons, nnpol);
            let mut p1 = MathVector::with_init(thefirstpt.max(0) as usize, lastp.max(0) as usize, 0.0);
            for i in thefirstpt..=lastp {
                p1.set_value(i as usize, the_param.value(i as usize));
            }
            sq1.perform_parameters(&p1);
            let c1 = sq1.bezier_value();
            v1 = super::tangency::tangency_vector(line, &c1, 0.0);

            let firstp = thefirstpt.max(thelastpt - nnp + 1);
            if firstp == thefirstpt && lastp == thelastpt {
                v2 = super::tangency::tangency_vector(line, &c1, 1.0);
            } else {
                let mut sq2 = LeastSquare::new(line, firstp, thelastpt, cons, cons, nnpol);
                let mut p2 =
                    MathVector::with_init(firstp.max(0) as usize, thelastpt.max(0) as usize, 0.0);
                for i in firstp..=thelastpt {
                    p2.set_value(i as usize, the_param.value(i as usize));
                }
                sq2.perform_parameters(&p2);
                let c2 = sq2.bezier_value();
                v2 = super::tangency::tangency_vector(line, &c2, 1.0);
            }

            let span = theknots[nbknots as usize - 1] - theknots[0];
            lambda1 = 1.0 / deg as f64 * (theknots[1] - theknots[0]) / span;
            lambda2 = 1.0 / deg as f64
                * (theknots[nbknots as usize - 1] - theknots[nbknots as usize - 2])
                / span;
        }

        if this.my_periodic {
            // `gxx:1404-1408`.
            for i in v1.lower()..=v1.upper() {
                let val = 0.5 * (v1.value(i) + v2.value(i));
                v1.set_value(i, val);
            }
            v2 = v1.clone();
        }

        // `gxx:1409-1416`.
        let mut sq = LeastSquare::with_knots(
            line,
            &theknots,
            &the_mults,
            thefirstpt,
            thelastpt,
            cons,
            cons,
            nbpoles,
        );
        sq.perform_tangents(&the_param, &v1, &v2, lambda1, lambda2);
        this.alldone = sq.is_done();
        this.the_multi_bsp_curve = sq.bspline_value().clone();
        let (_fv, t3, t2) = sq.error();
        this.current_tol3d = t3;
        this.current_tol2d = t2;
        this.tolreached = true;
    }

    // `gxx:1417-1421`.
    this.my_parameters = Some((
        the_param.lower() as i32,
        (the_param.lower()..=the_param.upper())
            .map(|i| the_param.value(i))
            .collect(),
    ));
}




