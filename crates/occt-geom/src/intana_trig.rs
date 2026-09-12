//! `math_TrigonometricFunctionRoots` (5-coeff) used by `IntAna_IntQuadQuad`.
//! Equation: `A cos^2 + 2 B cos sin + C cos + D sin + E = 0`.

use std::f64::consts::PI;

use occt_core::precision::PCONFUSION;

use super::quartic_roots;

const EPS: f64 = 1.5e-12;

/// Result of the trigonometric polynomial solver.
#[derive(Debug, Clone)]
pub(crate) enum TrigRoots {
    Fail,
    Infinite,
    Values(Vec<f64>),
}

pub(crate) fn trig_function_roots(
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    e: f64,
    inf: f64,
    sup: f64,
) -> TrigRoots {
    let depi = 2.0 * PI;
    let my_inf = inf;
    let mut delta = sup - inf;
    if delta > depi {
        delta = depi;
    }
    let md = inf / depi;

    if a.abs() <= EPS && b.abs() <= EPS {
        if c.abs() <= EPS {
            if d.abs() <= EPS {
                return if e.abs() <= EPS {
                    TrigRoots::Infinite
                } else {
                    TrigRoots::Values(Vec::new())
                };
            }
            let aa = -e / d;
            if aa.abs() > 1.0 {
                return TrigRoots::Values(Vec::new());
            }
            let z1 = aa.asin();
            let z2 = PI - z1;
            return TrigRoots::Values(collect_in_range(&[z1, z2], my_inf, delta, md, depi, inf, sup));
        }
        if d.abs() <= EPS {
            let aa = -e / c;
            if aa.abs() > 1.0 {
                return TrigRoots::Values(Vec::new());
            }
            let z1 = aa.acos();
            return TrigRoots::Values(collect_in_range(
                &[z1, -z1],
                my_inf,
                delta,
                md,
                depi,
                inf,
                sup,
            ));
        }
        let aa = e - c;
        let bb = 2.0 * d;
        let cc = e + c;
        let mut zer = super::quadratic_roots(aa, bb, cc);
        zer.retain(|t| t.is_finite());
        return wrap_tan_half(&zer, a, b, c, d, e, my_inf, delta, md, depi, inf, sup);
    }

    if a.abs() <= EPS && e.abs() <= EPS {
        if c.abs() <= EPS {
            let mut zer = vec![0.0, PI];
            let aa = -d / (b * 2.0);
            if aa.abs() <= 1.0 + PCONFUSION {
                if aa >= 1.0 {
                    zer.push(0.0);
                    zer.push(0.0);
                } else if aa <= -1.0 {
                    zer.push(PI);
                    zer.push(PI);
                } else {
                    let z = aa.acos();
                    zer.push(z);
                    zer.push(depi - z);
                }
            }
            return TrigRoots::Values(collect_in_range(&zer, my_inf, delta, md, depi, inf, sup));
        }
        if d.abs() <= EPS {
            let mut zer = vec![PI / 2.0, 1.5 * PI];
            let aa = -c / (b * 2.0);
            if aa.abs() <= 1.0 + PCONFUSION {
                if aa >= 1.0 {
                    zer.push(PI / 2.0);
                    zer.push(PI / 2.0);
                } else if aa <= -1.0 {
                    zer.push(1.5 * PI);
                    zer.push(1.5 * PI);
                } else {
                    let z = aa.asin();
                    zer.push(z);
                    zer.push(PI - z);
                }
            }
            return TrigRoots::Values(collect_in_range(&zer, my_inf, delta, md, depi, inf, sup));
        }
    }

    let mut ko = [a - c + e, 2.0 * d - 4.0 * b, 2.0 * e - 2.0 * a, 4.0 * b + 2.0 * d, a + c + e];
    let mut zer = Vec::new();
    for _ in 0..4 {
        zer = quartic_roots(ko[0], ko[1], ko[2], ko[3], ko[4]);
        zer.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
        let mut rescale = false;
        for w in zer.windows(2) {
            if (w[1] - w[0]).abs() < EPS {
                let qw = w[1];
                let va = ko[3] + qw * (2.0 * ko[2] + qw * (3.0 * ko[1] + qw * (4.0 * ko[0])));
                if va.abs() > EPS {
                    rescale = true;
                    break;
                }
            }
        }
        if !rescale {
            break;
        }
        for k in &mut ko {
            *k *= 0.0001;
        }
    }
    wrap_tan_half(&zer, a, b, c, d, e, my_inf, delta, md, depi, inf, sup)
}

fn collect_in_range(
    raw: &[f64],
    my_inf: f64,
    delta: f64,
    md: f64,
    depi: f64,
    inf: f64,
    sup: f64,
) -> Vec<f64> {
    let mut out = Vec::new();
    for &z0 in raw {
        let mut z = z0;
        if z <= -EPS {
            z = depi - z.abs();
        }
        z += md.trunc() * depi;
        let x = z - my_inf;
        if x >= -delta.abs() * f64::EPSILON && x <= delta + delta.abs() * f64::EPSILON {
            let mut t = z;
            if t < inf {
                t = inf;
            }
            if t > sup {
                t = sup;
            }
            push_unique(&mut out, t);
        }
    }
    out
}

fn wrap_tan_half(
    zer: &[f64],
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    e: f64,
    my_inf: f64,
    delta: f64,
    md: f64,
    depi: f64,
    inf: f64,
    sup: f64,
) -> TrigRoots {
    let mut sols = Vec::new();
    let supm = (sup - inf) * 0.01;
    for &z in zer {
        let mut teta = z.atan();
        teta += teta;
        if z <= -EPS {
            teta = depi - teta.abs();
        }
        teta += md.trunc() * depi;
        if teta - my_inf < 0.0 {
            teta += depi;
        }
        let x = teta - my_inf;
        if x >= -delta.abs() * f64::EPSILON && x <= delta + delta.abs() * f64::EPSILON {
            let polished = newton_trig(a, b, c, d, e, teta);
            let dn = polished - teta;
            let teta = if dn.abs() > supm { teta } else { polished };
            push_unique(&mut sols, teta);
        }
    }
    if sols.len() < 4 && (a - c + e).abs() <= EPS {
        let teta = PI + md.trunc() * depi;
        let x = teta - my_inf;
        if x >= -delta.abs() * f64::EPSILON && x <= delta + delta.abs() * f64::EPSILON {
            push_unique(&mut sols, teta);
        }
    }
    sols.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    TrigRoots::Values(sols)
}

fn newton_trig(a: f64, b: f64, c: f64, d: f64, e: f64, mut u: f64) -> f64 {
    for _ in 0..10 {
        let (si, co) = u.sin_cos();
        let f = a * co * co + 2.0 * b * co * si + c * co + d * si + e;
        let df = -2.0 * a * co * si + 2.0 * b * (co * co - si * si) - c * si + d * co;
        if df.abs() < 1e-18 {
            break;
        }
        u -= f / df;
    }
    u
}

fn push_unique(v: &mut Vec<f64>, t: f64) {
    if v.iter().any(|x| (*x - t).abs() <= EPS) {
        return;
    }
    let mut i = 0;
    while i < v.len() && v[i] < t {
        i += 1;
    }
    v.insert(i, t);
}
