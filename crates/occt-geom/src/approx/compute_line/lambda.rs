//! `Approx_ComputeLine` tangency helpers: `FirstTangencyVector`,
//! `LastTangencyVector`, `SearchFirstLambda` and `SearchLastLambda`.
//!
//! Source: `ModelingData/TKGeomBase/Approx/Approx_ComputeLine.gxx:430-715`.
//!
//! Translation notes:
//! - `V` is a `math_Vector` whose lower bound the caller fixes at 1; the
//!   `V(j)` writes below keep OCCT's `j` running from 1 (`gxx:458`,
//!   `:541`), while `SearchFirstLambda` / `SearchLastLambda` index from
//!   `V.Lower()` (`gxx:636`, `:694`).
//! - `AppParCurves_MultiCurve::D1` is
//!   [`MultiCurve::d1_3d`](crate::app_par_curves::multi_curve::MultiCurve::d1_3d)
//!   / `d1_2d`.
//! - `gp_Vec2d` has no `Dot` in `occt-core`, so `P1P2.Dot(myV)`
//!   (`gxx:651`, `:709`) is the `gp_XY::Dot` formula written out.

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};

use crate::approx::line_tool::LineTool;
use crate::app_par_curves::constraint::Constraint;
use crate::app_par_curves::least_square::LeastSquare;
use crate::approx::compute_line::ComputeLine;
use occt_math::MathVector;

impl ComputeLine {
    /// `FirstTangencyVector(Line, index, V)` (`gxx:430-512`).
    pub(crate) fn first_tangency_vector<L: LineTool>(
        &self,
        line: &L,
        index: i32,
        v: &mut MathVector,
    ) {
        let nb_p3d = line.nb_p3d();
        let nb_p2d = line.nb_p2d();

        let ok = if nb_p3d != 0 && nb_p2d != 0 {
            let tab = line.tangency(index);
            if let Some((tab_v, tab_v2d)) = tab {
                copy_vectors(v, nb_p3d, &tab_v, &tab_v2d);
                true
            } else {
                false
            }
        } else if nb_p2d != 0 {
            let tab = line.tangency_2d(index);
            if let Some(tab_v2d) = tab {
                copy_vectors(v, nb_p3d, &[], &tab_v2d);
                true
            } else {
                false
            }
        } else if nb_p3d != 0 {
            let tab = line.tangency_3d(index);
            if let Some(tab_v) = tab {
                copy_vectors(v, nb_p3d, &tab_v, &[]);
                true
            } else {
                false
            }
        } else {
            false
        };

        if !ok {
            // Search for a tangent vector by construction of a parabola.
            let nbpoles = 3;
            let mypar = self.compute_parameters(line, index, index + 2);
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
            write_parabola_tangent(v, &c, nb_p3d, nb_p2d, 0.0);
        }
    }

    /// `LastTangencyVector(Line, index, V)` (`gxx:514-596`).
    pub(crate) fn last_tangency_vector<L: LineTool>(
        &self,
        line: &L,
        index: i32,
        v: &mut MathVector,
    ) {
        let nb_p3d = line.nb_p3d();
        let nb_p2d = line.nb_p2d();

        let ok = if nb_p3d != 0 && nb_p2d != 0 {
            let tab = line.tangency(index);
            if let Some((tab_v, tab_v2d)) = tab {
                copy_vectors(v, nb_p3d, &tab_v, &tab_v2d);
                true
            } else {
                false
            }
        } else if nb_p2d != 0 {
            let tab = line.tangency_2d(index);
            if let Some(tab_v2d) = tab {
                copy_vectors(v, nb_p3d, &[], &tab_v2d);
                true
            } else {
                false
            }
        } else if nb_p3d != 0 {
            let tab = line.tangency_3d(index);
            if let Some(tab_v) = tab {
                copy_vectors(v, nb_p3d, &tab_v, &[]);
                true
            } else {
                false
            }
        } else {
            false
        };

        if !ok {
            // Search for a tangent vector by construction of a parabola.
            let nbpoles = 3;
            let mypar = self.compute_parameters(line, index - 2, index);
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
            write_parabola_tangent(v, &c, nb_p3d, nb_p2d, 1.0);
        }
    }

    /// `SearchFirstLambda(Line, TheParam, V, index)` (`gxx:597-656`):
    /// `dq/dw = lambda * V = (p2 - p1) / (u2 - u1)`.
    pub(crate) fn search_first_lambda<L: LineTool>(
        &self,
        line: &L,
        the_param: &MathVector,
        v: &MathVector,
        index: i32,
    ) -> f64 {
        let nb_p3d = line.nb_p3d();
        let nb_p2d = line.nb_p2d();

        let (tab_p1, tab_p12d) = line_values(line, index, nb_p3d, nb_p2d);
        let (tab_p2, tab_p22d) = line_values(line, index + 1, nb_p3d, nb_p2d);

        let u1 = the_param.value(index as usize);
        let u2 = the_param.value((index + 1) as usize);
        let low = v.lower();

        let (lambda, s) = if nb_p3d != 0 {
            let p1 = tab_p1[0];
            let p2 = tab_p2[0];
            let p1p2 = GpVec::from_pnts(&p1, &p2);
            let my_v = GpVec::new(v.value(low), v.value(low + 1), v.value(low + 2));
            let lambda = p1p2.magnitude() / (my_v.magnitude() * (u2 - u1));
            let s = if p1p2.dot(&my_v) > 0.0 { 1.0 } else { -1.0 };
            (lambda, s)
        } else {
            let p12d = tab_p12d[0];
            let p22d = tab_p22d[0];
            let p1p2 = vec2d_from_pnts(&p12d, &p22d);
            let my_v = GpVec2d::new(v.value(low), v.value(low + 1));
            let lambda = p1p2.magnitude() / (my_v.magnitude() * (u2 - u1));
            let dot = p1p2.x() * my_v.x() + p1p2.y() * my_v.y();
            let s = if dot > 0.0 { 1.0 } else { -1.0 };
            (lambda, s)
        };
        s * lambda
    }

    /// `SearchLastLambda(Line, TheParam, V, index)` (`gxx:657-715`).
    pub(crate) fn search_last_lambda<L: LineTool>(
        &self,
        line: &L,
        the_param: &MathVector,
        v: &MathVector,
        index: i32,
    ) -> f64 {
        let nb_p3d = line.nb_p3d();
        let nb_p2d = line.nb_p2d();

        let (tab_p, tab_p2d) = line_values(line, index - 1, nb_p3d, nb_p2d);
        let (tab_p2, tab_p22d) = line_values(line, index, nb_p3d, nb_p2d);

        let u1 = the_param.value((index - 1) as usize);
        let u2 = the_param.value(index as usize);
        let low = v.lower();

        let (lambda, s) = if nb_p3d != 0 {
            let p1 = tab_p[0];
            let p2 = tab_p2[0];
            let p1p2 = GpVec::from_pnts(&p1, &p2);
            let my_v = GpVec::new(v.value(low), v.value(low + 1), v.value(low + 2));
            let lambda = p1p2.magnitude() / (my_v.magnitude() * (u2 - u1));
            let s = if p1p2.dot(&my_v) > 0.0 { 1.0 } else { -1.0 };
            (lambda, s)
        } else {
            let p12d = tab_p2d[0];
            let p22d = tab_p22d[0];
            let p1p2 = vec2d_from_pnts(&p12d, &p22d);
            let my_v = GpVec2d::new(v.value(low), v.value(low + 1));
            let lambda = p1p2.magnitude() / (my_v.magnitude() * (u2 - u1));
            let dot = p1p2.x() * my_v.x() + p1p2.y() * my_v.y();
            let s = if dot > 0.0 { 1.0 } else { -1.0 };
            (lambda, s)
        };
        s * lambda
    }
}

/// `LineTool::Value(Line, i, tabP, tabP2d)` dispatched on the column counts
/// (`gxx:618-630`, `:676-688`).
fn line_values<L: LineTool>(
    line: &L,
    index: i32,
    nb_p3d: i32,
    nb_p2d: i32,
) -> (Vec<GpPnt>, Vec<GpPnt2d>) {
    if nb_p3d != 0 && nb_p2d != 0 {
        line.value(index)
    } else if nb_p2d != 0 {
        (Vec::new(), line.value_2d(index))
    } else if nb_p3d != 0 {
        (line.value_3d(index), Vec::new())
    } else {
        (Vec::new(), Vec::new())
    }
}

/// The `if (Ok)` copy loops of `FirstTangencyVector` / `LastTangencyVector`
/// (`gxx:456-474`, `:538-558`).
fn copy_vectors(v: &mut MathVector, nb_p3d: i32, tab_v: &[GpVec], tab_v2d: &[GpVec2d]) {
    if nb_p3d != 0 {
        let mut j = 1i32;
        for vec in tab_v {
            v.set_value(j as usize, vec.x());
            v.set_value((j + 1) as usize, vec.y());
            v.set_value((j + 2) as usize, vec.z());
            j += 3;
        }
    }
    if !tab_v2d.is_empty() {
        let mut j = nb_p3d * 3 + 1;
        for vec in tab_v2d {
            v.set_value(j as usize, vec.x());
            v.set_value((j + 1) as usize, vec.y());
            j += 2;
        }
    }
}

/// The `else` branch shared by both tangency-vector helpers: approximate the
/// parabola `index-2..index` with a 3-pole `PassPoint` least square and read
/// its derivative at `u` (`gxx:494-511`, `:577-594`).
fn write_parabola_tangent(
    v: &mut MathVector,
    c: &crate::app_par_curves::multi_curve::MultiCurve,
    nb_p3d: i32,
    nb_p2d: i32,
    u: f64,
) {
    let mut j = 1i32;
    for i in 1..=nb_p3d {
        let (_, my_v) = c.d1_3d(i, u);
        v.set_value(j as usize, my_v.x());
        v.set_value((j + 1) as usize, my_v.y());
        v.set_value((j + 2) as usize, my_v.z());
        j += 3;
    }
    j = nb_p3d * 3 + 1;
    for i in nb_p3d + 1..=nb_p3d + nb_p2d {
        let (_, my_v2d) = c.d1_2d(i, u);
        v.set_value(j as usize, my_v2d.x());
        v.set_value((j + 1) as usize, my_v2d.y());
        j += 2;
    }
}

/// `gp_Vec2d(P1, P2)` (`gp_Vec2d.hxx`).
fn vec2d_from_pnts(p1: &GpPnt2d, p2: &GpPnt2d) -> GpVec2d {
    GpVec2d::new(p2.x() - p1.x(), p2.y() - p1.y())
}





