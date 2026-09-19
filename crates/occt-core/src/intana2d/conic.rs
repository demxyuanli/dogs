//! Conic defined by its implicit quadratic equation
//! `A.X**2 + B.Y**2 + 2.C.X*Y + 2.D.X + 2.E.Y + F = 0`.
//! Source: `IntAna2d_Conic.hxx` + `IntAna2d_Conic.cxx`.
use crate::gp::{GpAx2d, GpCirc2d, GpElips2d, GpHypr2d, GpLin2d, GpParab2d, GpXY};
use crate::precision::REAL_SMALL;

/// The six values `(T11, T12, T13, T21, T22, T23)` of
/// `gp_Trsf2d::SetTransformation(const gp_Ax2d&)` (`gp_Trsf2d.cxx:72-84`): the
/// change of basis from the basic frame to the frame `theA`
/// (`gp_Trsf2d.hxx:78-82`). `gp_Trsf2d::Value(Row, Col)`
/// (`gp_Trsf2d.hxx:301-312`) is `scale * matrix[Row-1][Col-1]` for `Col < 3` and
/// `loc.Coord(Row)` for `Col == 3`.
fn trsf2d_values(the_a: &GpAx2d) -> (f64, f64, f64, f64, f64, f64) {
    let v1x = the_a.vdir.x;
    let v1y = the_a.vdir.y;
    let lx = the_a.loc.x();
    let ly = the_a.loc.y();
    // `SetCol(1, V1)`, `SetCol(2, (-V1.Y(), V1.X()))`, then `Transpose()`.
    let t11 = v1x;
    let t12 = v1y;
    let t21 = -v1y;
    let t22 = v1x;
    // `loc` starts as `Location`, is multiplied by the transposed matrix, then reversed.
    let t13 = -(t11 * lx + t12 * ly);
    let t23 = -(t21 * lx + t22 * ly);
    (t11, t12, t13, t21, t22, t23)
}

/// `gp_Ax22d::XAxis()` — location plus X direction (`gp_Elips2d.hxx:211`,
/// `gp_Hypr2d.hxx:302`, `gp_Parab2d.hxx:167`).
fn x_axis22(pos: &crate::gp::GpAx22d) -> GpAx2d {
    GpAx2d::new(pos.point, pos.vxdir)
}

/// `gp_Lin2d::Coefficients` (`gp_Lin2d.hxx:91-96`).
fn lin2d_coefficients(l: &GpLin2d) -> (f64, f64, f64) {
    let a = l.pos.vdir.y;
    let b = -l.pos.vdir.x;
    let c = -(a * l.pos.loc.x() + b * l.pos.loc.y());
    (a, b, c)
}

/// `gp_Circ2d::Coefficients` (`gp_Circ2d.hxx:255-271`).
fn circ2d_coefficients(c: &GpCirc2d) -> (f64, f64, f64, f64, f64, f64) {
    let xc = c.pos.point.x();
    let yc = c.pos.point.y();
    (
        1.0,
        1.0,
        0.0,
        -xc,
        -yc,
        xc * xc + yc * yc - c.radius * c.radius,
    )
}

/// `gp_Elips2d::Coefficients` (`gp_Elips2d.cxx:25-63`).
fn elips2d_coefficients(e: &GpElips2d) -> (f64, f64, f64, f64, f64, f64) {
    let d_min = e.minor_radius * e.minor_radius;
    let d_maj = e.major_radius * e.major_radius;
    if d_min <= REAL_SMALL && d_maj <= REAL_SMALL {
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    }
    let (t11, t12, t13, t21, t22, t23) = trsf2d_values(&x_axis22(&e.pos));
    if d_min <= REAL_SMALL {
        (
            t11 * t11,
            t12 * t12,
            t11 * t12,
            t11 * t13,
            t12 * t13,
            t13 * t13 - d_maj,
        )
    } else {
        (
            (t11 * t11 / d_maj) + (t21 * t21 / d_min),
            (t12 * t12 / d_maj) + (t22 * t22 / d_min),
            (t11 * t12 / d_maj) + (t21 * t22 / d_min),
            (t11 * t13 / d_maj) + (t21 * t23 / d_min),
            (t12 * t13 / d_maj) + (t22 * t23 / d_min),
            (t13 * t13 / d_maj) + (t23 * t23 / d_min) - 1.0,
        )
    }
}

/// `gp_Hypr2d::Coefficients` (`gp_Hypr2d.cxx:25-62`).
fn hypr2d_coefficients(h: &GpHypr2d) -> (f64, f64, f64, f64, f64, f64) {
    let d_min = h.minor_radius * h.minor_radius;
    let d_maj = h.major_radius * h.major_radius;
    if d_min <= REAL_SMALL && d_maj <= REAL_SMALL {
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    }
    let (t11, t12, t13, t21, t22, t23) = trsf2d_values(&x_axis22(&h.pos));
    if d_min <= REAL_SMALL {
        (
            t11 * t11,
            t12 * t12,
            t11 * t12,
            t11 * t13,
            t12 * t13,
            t13 * t13 - d_maj,
        )
    } else {
        (
            (t11 * t11 / d_maj) - (t21 * t21 / d_min),
            (t12 * t12 / d_maj) - (t22 * t22 / d_min),
            (t11 * t12 / d_maj) - (t21 * t22 / d_min),
            (t11 * t13 / d_maj) - (t21 * t23 / d_min),
            (t12 * t13 / d_maj) - (t22 * t23 / d_min),
            (t13 * t13 / d_maj) - (t23 * t23 / d_min) - 1.0,
        )
    }
}

/// `gp_Parab2d::Coefficients` (`gp_Parab2d.cxx:44-62`).
fn parab2d_coefficients(p: &GpParab2d) -> (f64, f64, f64, f64, f64, f64) {
    let focal_param = 2.0 * p.focal;
    let (t11, t12, t13, t21, t22, t23) = trsf2d_values(&x_axis22(&p.pos));
    (
        t21 * t21,
        t22 * t22,
        t21 * t22,
        (t21 * t23) - (focal_param * t11),
        (t22 * t23) - (focal_param * t12),
        (t23 * t23) - (2.0 * focal_param * t13),
    )
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IntAna2dConic {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub d: f64,
    pub e: f64,
    pub f: f64,
}

impl IntAna2dConic {
    /// Ctor from `gp_Lin2d` (`cxx:22-31`): `a = b = c = 0`, then
    /// `f = 2 * f_lin`.
    pub fn from_lin2d(l: &GpLin2d) -> Self {
        let (d, e, f) = lin2d_coefficients(l);
        Self {
            a: 0.0,
            b: 0.0,
            c: 0.0,
            d,
            e,
            f: 2.0 * f,
        }
    }

    /// Ctor from `gp_Circ2d` (`cxx:33-37`).
    pub fn from_circ2d(c: &GpCirc2d) -> Self {
        let (a, b, c_, d, e, f) = circ2d_coefficients(c);
        Self {
            a,
            b,
            c: c_,
            d,
            e,
            f,
        }
    }

    /// Ctor from `gp_Elips2d` (`cxx:39-43`).
    pub fn from_elips2d(el: &GpElips2d) -> Self {
        let (a, b, c, d, e, f) = elips2d_coefficients(el);
        Self { a, b, c, d, e, f }
    }

    /// Ctor from `gp_Parab2d` (`cxx:45-48`).
    pub fn from_parab2d(p: &GpParab2d) -> Self {
        let (a, b, c, d, e, f) = parab2d_coefficients(p);
        Self { a, b, c, d, e, f }
    }

    /// Ctor from `gp_Hypr2d` (`cxx:50-53`).
    pub fn from_hypr2d(h: &GpHypr2d) -> Self {
        let (a, b, c, d, e, f) = hypr2d_coefficients(h);
        Self { a, b, c, d, e, f }
    }

    /// `Coefficients` (`cxx:114-124`).
    pub fn coefficients(&self) -> (f64, f64, f64, f64, f64, f64) {
        (self.a, self.b, self.c, self.d, self.e, self.f)
    }

    /// `Value(X, Y)` (`cxx:92-98`).
    pub fn value(&self, x: f64, y: f64) -> f64 {
        let (a, b, c, d, e, f) = self.coefficients();
        a * x * x + b * y * y + 2.0 * c * x * y + 2.0 * d * x + 2.0 * e * y + f
    }

    /// `Grad(X, Y)` (`cxx:100-106`).
    pub fn grad(&self, x: f64, y: f64) -> GpXY {
        let (a, b, c, d, e, _f) = self.coefficients();
        GpXY::new(
            2.0 * a * x + 2.0 * c * y + 2.0 * d,
            2.0 * b * y + 2.0 * c * x + 2.0 * e,
        )
    }

    /// `ValAndGrad(X, Y, Val, Grd)` (`cxx:108-112`).
    pub fn val_and_grad(&self, x: f64, y: f64) -> (f64, GpXY) {
        let (la, lb, lc, ld, le, lf) = self.coefficients();
        let grd = GpXY::new(
            2.0 * la * x + 2.0 * lc * y + 2.0 * ld,
            2.0 * lb * y + 2.0 * lc * x + 2.0 * le,
        );
        let val = la * x * x + lb * y * y + 2.0 * lc * x * y + 2.0 * ld * x + 2.0 * le * y + lf;
        (val, grd)
    }

    /// `NewCoefficients(A, B, C, D, E, F, Axis)` (`cxx:55-91`): the same conic
    /// written in the local frame `Dir1`. The in/out parameters `A..F` become
    /// the return tuple.
    pub fn new_coefficients(&self, dir1: &GpAx2d) -> (f64, f64, f64, f64, f64, f64) {
        let (a, b, c, d, e, f) = self.coefficients();
        // x = t11 X + t12 Y + t13 ; y = t21 X + t22 Y + t23
        let t11 = dir1.vdir.x;
        let t21 = dir1.vdir.y;
        let t13 = dir1.loc.x();
        let t23 = dir1.loc.y();
        let t22 = t11;
        let t12 = -t21;

        let a1 = t11 * (a * t11 + 2.0 * c * t21) + b * t21 * t21;
        let b1 = t12 * (a * t12 + 2.0 * c * t22) + b * t22 * t22;
        let c1 = t12 * (a * t11 + c * t21) + t22 * (c * t11 + b * t21);
        let d1 = t11 * (d + a * t13) + t21 * (e + c * t13) + t23 * (c * t11 + b * t21);
        let e1 = t12 * (d + a * t13) + t22 * (e + c * t13) + t23 * (c * t12 + b * t22);
        let f1 = f + t13 * (2.0 * d + a * t13) + t23 * (2.0 * e + 2.0 * c * t13 + b * t23);
        (a1, b1, c1, d1, e1, f1)
    }
}
