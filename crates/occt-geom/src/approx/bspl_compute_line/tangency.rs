//! `Approx_BSplComputeLine` tangency helpers: `FirstTangencyVector`,
//! `LastTangencyVector`, `TangencyVector`, `SearchFirstLambda` and
//! `SearchLastLambda` (`Approx_BSplComputeLine.gxx:197-490`, `:1427-1458`).

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_math::MathVector;

use crate::approx::line_tool::LineTool;
use crate::app_par_curves::constraint::Constraint;
use crate::app_par_curves::least_square::LeastSquare;
use crate::app_par_curves::multi_curve::MultiCurve;
use crate::approx::bspl_compute_line::BsplComputeLine;

/// `3 * NbP3d + 2 * NbP2d` (`gxx:1344`).
fn nb_of_floats<L: LineTool>(line: &L) -> usize {
    (3 * line.nb_p3d() + 2 * line.nb_p2d()).max(1) as usize
}

/// `FirstTangencyVector(Line, index, V)` (`gxx:197-279`).
pub(super) fn first_tangency_vector<L: LineTool>(this: &BsplComputeLine, line: &L, index: i32) -> MathVector {
    let nb_p3d = line.nb_p3d();
    let nb_p2d = line.nb_p2d();
    let mut v = MathVector::with_init(1, nb_of_floats(line), 0.0);

    let tang = if nb_p3d != 0 && nb_p2d != 0 {
        line.tangency(index)
    } else if nb_p2d != 0 {
        line.tangency_2d(index).map(|v2| (Vec::new(), v2))
    } else if nb_p3d != 0 {
        line.tangency_3d(index).map(|v3| (v3, Vec::new()))
    } else {
        None
    };

    if let Some((tab_v, tab_v2d)) = tang {
        // `gxx:222-243`.
        let mut j = 1usize;
        for t in tab_v.iter() {
            v.set_value(j, t.x());
            v.set_value(j + 1, t.y());
            v.set_value(j + 2, t.z());
            j += 3;
        }
        let mut j = (nb_p3d * 3 + 1) as usize;
        for t in tab_v2d.iter() {
            v.set_value(j, t.x());
            v.set_value(j + 1, t.y());
            j += 2;
        }
    } else {
        // Search for a tangent vector by construction of a parabola
        // (`gxx:246-277`).
        let nbpoles = 3;
        let mut mypar = MathVector::with_init(index.max(0) as usize, (index + 2).max(0) as usize, 0.0);
        let params = super::parameters::compute_parameters(this, line, index, index + 2);
        for (k, val) in params.iter().enumerate() {
            mypar.set_value(index as usize + k, *val);
        }
        let mut lsq = LeastSquare::with_parameters(
            line,
            index,
            index + 2,
            Constraint::PassPoint,
            Constraint::PassPoint,
            &mypar,
            nbpoles,
        );
        let c = lsq.bezier_value();

        let mut j = 1usize;
        for i in 1..=nb_p3d {
            let (_p, my_v) = c.d1_3d(i, 0.0);
            v.set_value(j, my_v.x());
            v.set_value(j + 1, my_v.y());
            v.set_value(j + 2, my_v.z());
            j += 3;
        }
        let mut j = (nb_p3d * 3 + 1) as usize;
        for i in nb_p3d + 1..=nb_p3d + nb_p2d {
            let (_p2d, my_v2d) = c.d1_2d(i, 0.0);
            v.set_value(j, my_v2d.x());
            v.set_value(j + 1, my_v2d.y());
            j += 2;
        }
    }

    v
}

/// `LastTangencyVector(Line, index, V)` (`gxx:283-364`).
pub(super) fn last_tangency_vector<L: LineTool>(this: &BsplComputeLine, line: &L, index: i32) -> MathVector {
    let nb_p3d = line.nb_p3d();
    let nb_p2d = line.nb_p2d();
    let mut v = MathVector::with_init(1, nb_of_floats(line), 0.0);

    let tang = if nb_p3d != 0 && nb_p2d != 0 {
        line.tangency(index)
    } else if nb_p2d != 0 {
        line.tangency_2d(index).map(|v2| (Vec::new(), v2))
    } else if nb_p3d != 0 {
        line.tangency_3d(index).map(|v3| (v3, Vec::new()))
    } else {
        None
    };

    if let Some((tab_v, tab_v2d)) = tang {
        let mut j = 1usize;
        for t in tab_v.iter() {
            v.set_value(j, t.x());
            v.set_value(j + 1, t.y());
            v.set_value(j + 2, t.z());
            j += 3;
        }
        let mut j = (nb_p3d * 3 + 1) as usize;
        for t in tab_v2d.iter() {
            v.set_value(j, t.x());
            v.set_value(j + 1, t.y());
            j += 2;
        }
    } else {
        // `gxx:331-362`: parabola over the last three points, `D1` at `u = 1`.
        let nbpoles = 3;
        let mut mypar = MathVector::with_init(
            (index - 2).max(0) as usize,
            index.max(0) as usize,
            0.0,
        );
        let params = super::parameters::compute_parameters(this, line, index - 2, index);
        for (k, val) in params.iter().enumerate() {
            mypar.set_value((index - 2) as usize + k, *val);
        }
        let mut lsq = LeastSquare::with_parameters(
            line,
            index - 2,
            index,
            Constraint::PassPoint,
            Constraint::PassPoint,
            &mypar,
            nbpoles,
        );
        let c = lsq.bezier_value();

        let mut j = 1usize;
        for i in 1..=nb_p3d {
            let (_p, my_v) = c.d1_3d(i, 1.0);
            v.set_value(j, my_v.x());
            v.set_value(j + 1, my_v.y());
            v.set_value(j + 2, my_v.z());
            j += 3;
        }
        let mut j = (nb_p3d * 3 + 1) as usize;
        for i in nb_p3d + 1..=nb_p3d + nb_p2d {
            let (_p2d, my_v2d) = c.d1_2d(i, 1.0);
            v.set_value(j, my_v2d.x());
            v.set_value(j + 1, my_v2d.y());
            j += 2;
        }
    }

    v
}

/// `TangencyVector(Line, C, U, V)` (`gxx:1427-1458`).
pub(super) fn tangency_vector<L: LineTool>(line: &L, c: &MultiCurve, u: f64) -> MathVector {
    let nb_p3d = line.nb_p3d();
    let nb_p2d = line.nb_p2d();
    let mut v = MathVector::with_init(1, nb_of_floats(line), 0.0);

    let mut j = 1usize;
    for i in 1..=nb_p3d {
        let (_p, my_v) = c.d1_3d(i, u);
        v.set_value(j, my_v.x());
        v.set_value(j + 1, my_v.y());
        v.set_value(j + 2, my_v.z());
        j += 3;
    }
    let mut j = (nb_p3d * 3 + 1) as usize;
    for i in nb_p3d + 1..=nb_p3d + nb_p2d {
        let (_p2d, my_v2d) = c.d1_2d(i, u);
        v.set_value(j, my_v2d.x());
        v.set_value(j + 1, my_v2d.y());
        j += 2;
    }

    v
}

/// `SearchFirstLambda(Line, aPar, Theknots, V, index)` (`gxx:369-429`).
pub(super) fn search_first_lambda<L: LineTool>(
    line: &L,
    a_par: &MathVector,
    theknots: &[f64],
    v: &MathVector,
    index: i32,
) -> f64 {
    let nb_p3d = line.nb_p3d();

    let (tab_p1, tab_p12d) = line.value(index);
    let (tab_p2, tab_p22d) = line.value(index + 1);

    let u1 = a_par.value(index as usize);
    let u2 = a_par.value((index + 1) as usize);
    let low = v.lower();
    let nbknots = theknots.len();

    let (lambda, s) = if nb_p3d != 0 {
        let p1: GpPnt = tab_p1[0];
        let p2: GpPnt = tab_p2[0];
        let p1p2 = GpVec::new(p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z());
        let my_v = GpVec::new(v.value(low), v.value(low + 1), v.value(low + 2));
        (
            p1p2.magnitude() / (my_v.magnitude() * (u2 - u1)),
            if p1p2.dot(&my_v) > 0.0 { 1.0 } else { -1.0 },
        )
    } else {
        let p12d: GpPnt2d = tab_p12d[0];
        let p22d: GpPnt2d = tab_p22d[0];
        let p1p2 = GpVec2d::new(p22d.x() - p12d.x(), p22d.y() - p12d.y());
        let my_v = GpVec2d::new(v.value(low), v.value(low + 1));
        (
            p1p2.magnitude() / (my_v.magnitude() * (u2 - u1)),
            if p1p2.dot(&my_v) > 0.0 { 1.0 } else { -1.0 },
        )
    };

    (s * lambda) * (theknots[1] - theknots[0]) / (theknots[nbknots - 1] - theknots[0])
}

/// `SearchLastLambda(Line, aPar, Theknots, V, index)` (`gxx:433-493`).
pub(super) fn search_last_lambda<L: LineTool>(
    line: &L,
    a_par: &MathVector,
    theknots: &[f64],
    v: &MathVector,
    index: i32,
) -> f64 {
    let nb_p3d = line.nb_p3d();
    let nb_p2d = line.nb_p2d();

    let (tab_p, tab_p2d) = line.value(index - 1);
    let (tab_p2, tab_p22d) = line.value(index);

    let u1 = a_par.value((index - 1) as usize);
    let u2 = a_par.value(index as usize);
    let low = v.lower();
    let nbknots = theknots.len();

    let (lambda, s) = if nb_p3d != 0 {
        let p1: GpPnt = tab_p[0];
        let p2: GpPnt = tab_p2[0];
        let p1p2 = GpVec::new(p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z());
        let my_v = GpVec::new(v.value(low), v.value(low + 1), v.value(low + 2));
        (
            p1p2.magnitude() / (my_v.magnitude() * (u2 - u1)),
            if p1p2.dot(&my_v) > 0.0 { 1.0 } else { -1.0 },
        )
    } else {
        let _ = nb_p2d;
        let p12d: GpPnt2d = tab_p2d[0];
        let p22d: GpPnt2d = tab_p22d[0];
        let p1p2 = GpVec2d::new(p22d.x() - p12d.x(), p22d.y() - p12d.y());
        let my_v = GpVec2d::new(v.value(low), v.value(low + 1));
        (
            p1p2.magnitude() / (my_v.magnitude() * (u2 - u1)),
            if p1p2.dot(&my_v) > 0.0 { 1.0 } else { -1.0 },
        )
    };

    (s * lambda) * (theknots[nbknots - 1] - theknots[nbknots - 2])
        / (theknots[nbknots - 1] - theknots[0])
}




