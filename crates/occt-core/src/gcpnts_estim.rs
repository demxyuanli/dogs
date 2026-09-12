//! `GCPnts_TangentialDeflection::EstimDefl` (`cxx:960-1010`).
//!
//! Uses `math_BrentMinimum` (`math_BrentMinimum.cxx:80-187`, `.lxx:17-21`)
//! and, when `!IsDone`, `math_PSO` (`math_PSO.cxx:58-268`) then Brent again.

use crate::gp::GpPnt;
use crate::precision::RESOLUTION;

/// Squared distance from `p` to the infinite line `p1..p2`
/// (`GCPnts_DistFunction` / `gp_Lin::SquareDistance`).
fn line_sq_dist(p1: &GpPnt, p2: &GpPnt, p: &GpPnt) -> f64 {
    let vx = p2.x() - p1.x();
    let vy = p2.y() - p1.y();
    let vz = p2.z() - p1.z();
    let len2 = vx * vx + vy * vy + vz * vz;
    let (ax, ay, az) = if len2 > RESOLUTION {
        (vx, vy, vz)
    } else {
        return p.square_distance(p1);
    };
    let apx = p.x() - p1.x();
    let apy = p.y() - p1.y();
    let apz = p.z() - p1.z();
    let t = (apx * ax + apy * ay + apz * az) / (ax * ax + ay * ay + az * az);
    let dx = apx - t * ax;
    let dy = apy - t * ay;
    let dz = apz - t * az;
    dx * dx + dy * dy + dz * dz
}

/// `EstimDefl` (`cxx:960-1010`).
pub(super) fn estim_defl<C: super::CurveSecondDeriv>(
    curve: &C,
    u1: f64,
    u2: f64,
    u_tol: f64,
    curve_span: f64,
) -> (f64, f64) {
    let p1 = curve.point(u1);
    let mut p2 = curve.point(u2);
    // `GCPnts_DistFunction.cxx:27-36`: coincident ends use U1+0.01*(U2-U1).
    if p1.square_distance(&p2) <= RESOLUTION {
        p2 = curve.point(u1 + 0.01 * (u2 - u1));
    }
    let f = |u: f64| -> Option<f64> {
        // `DistFunction::Value` (`cxx:43-45`): reject outside [U1,U2] so Brent
        // `!OK` leaves `Done=false` and EstimDefl falls through to PSO.
        if u < u1 || u > u2 {
            return None;
        }
        Some(-line_sq_dist(&p1, &p2, &curve.point(u)))
    };
    let mid = 0.5 * (u1 + u2);
    let denom = u1.abs() + u2.abs();
    let a_rel_tol = 1.0e-3_f64.max(if denom > 0.0 {
        2.0 * u_tol / denom
    } else {
        1.0e-3
    });
    // `math_BrentMinimum(aRelTol, 100, myUTol)`: ZEPS is `myUTol` (`cxx:972-974`).
    if let Some((loc, minv)) = brent_minimum(f, u1, mid, u2, a_rel_tol, u_tol, 100) {
        return ((-minv).max(0.0).sqrt(), loc);
    }
    let du = curve_span.abs().max(u_tol);
    let step = (0.1 * du).max(100.0 * u_tol);
    // `RealToInt(32*(U2-U1)/Du)` (`cxx:985`).
    let n_particles = 8i32.max((32.0 * (u2 - u1) / du) as i32).max(1) as usize;
    let (value, t) = pso_1d(f, u1, u2, step, n_particles, 100);
    if let Some((loc, minv)) = brent_minimum(
        f,
        (t - step).max(u1),
        t,
        (t + step).min(u2),
        a_rel_tol,
        u_tol,
        100,
    ) {
        return ((-minv).max(0.0).sqrt(), loc);
    }
    ((-value).max(0.0).sqrt(), t)
}

/// `math_BrentMinimum::Perform` (`math_BrentMinimum.cxx:80-187`).
/// `IsSolutionReached` is `.lxx:17-21`, not the Numerical Recipes mid test.
/// `None` is `Done=false` (max iter or `F.Value` failed).
fn brent_minimum<F: Fn(f64) -> Option<f64>>(
    f: F,
    ax: f64,
    bx: f64,
    cx: f64,
    x_tol: f64,
    zeps: f64,
    max_iter: usize,
) -> Option<(f64, f64)> {
    const CGOLD: f64 = 0.3819660;
    let mut a = if ax < cx { ax } else { cx };
    let mut b = if ax > cx { ax } else { cx };
    let mut x = bx;
    let mut w = bx;
    let mut v = bx;
    let fx = f(x)?;
    let mut fx = fx;
    let mut fw = fx;
    let mut fv = fx;
    let mut e: f64 = 0.0;
    let mut d = f64::MAX;
    for _iter in 1..=max_iter {
        let xm = 0.5 * (a + b);
        let tol1 = x_tol * x.abs() + zeps;
        let two_tol = 2.0 * (x_tol * x.abs() + zeps);
        if x <= two_tol + a && x >= b - two_tol {
            return Some((x, fx));
        }
        let mut u;
        if e.abs() > tol1 {
            let r = (x - w) * (fx - fv);
            let mut q = (x - v) * (fx - fw);
            let mut p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if q > 0.0 {
                p = -p;
            }
            q = q.abs();
            let etemp = e;
            e = d;
            if p.abs() >= (0.5 * q * etemp).abs() || p <= q * (a - x) || p >= q * (b - x) {
                e = if x >= xm { a - x } else { b - x };
                d = CGOLD * e;
            } else {
                d = p / q;
                u = x + d;
                let tol2 = 2.0 * tol1;
                if u - a < tol2 || b - u < tol2 {
                    d = copysign_occt(tol1, xm - x);
                }
            }
        } else {
            e = if x >= xm { a - x } else { b - x };
            d = CGOLD * e;
        }
        u = if d.abs() >= tol1 {
            x + d
        } else {
            x + copysign_occt(tol1, d)
        };
        let fu = f(u)?;
        if fu <= fx {
            if u >= x {
                a = x;
            } else {
                b = x;
            }
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
        } else {
            if u < x {
                a = u;
            } else {
                b = u;
            }
            if fu <= fw || w == x {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            } else if fu <= fv || v == x || v == w {
                v = u;
                fv = fu;
            }
        }
    }
    None
}

fn copysign_occt(mag: f64, sgn: f64) -> f64 {
    mag.copysign(sgn)
}

/// `math_BullardGenerator` (`math_BullardGenerator.hxx:27-51`). Seed 1.
struct Bullard {
    hi: u32,
    lo: u32,
}

impl Bullard {
    fn new() -> Self {
        let mut g = Self { hi: 1, lo: 0 };
        g.hi = 1;
        g.lo = 1 ^ 0x49616E42;
        g
    }

    fn next_int(&mut self) -> u32 {
        self.hi = self.hi.wrapping_shr(2).wrapping_add(self.hi.wrapping_shl(2));
        self.hi = self.hi.wrapping_add(self.lo);
        self.lo = self.lo.wrapping_add(self.hi);
        self.hi
    }

    fn next_real(&mut self) -> f64 {
        self.next_int() as f64 / u32::MAX as f64
    }
}

#[derive(Clone)]
struct Particle {
    position: f64,
    velocity: f64,
    best_position: f64,
    distance: f64,
    best_distance: f64,
}

impl Particle {
    fn new() -> Self {
        Self {
            position: 0.0,
            velocity: 0.0,
            best_position: 0.0,
            distance: f64::MAX,
            best_distance: f64::MAX,
        }
    }
}

/// 1-variable `math_PSO::Perform` + `performPSOWithGivenParticles`
/// (`math_PSO.cxx:58-268`).
fn pso_1d<F: Fn(f64) -> Option<f64>>(
    f: F,
    low: f64,
    upp: f64,
    step: f64,
    nb_particles: usize,
    nb_iter: usize,
) -> (f64, f64) {
    const BORDER_DIVISOR: f64 = 1.0e4;
    let span = upp - low;
    let min_uv = low + span / BORDER_DIVISOR;
    let max_uv = upp - span / BORDER_DIVISOR;
    let mut pool = vec![Particle::new(); nb_particles.max(1)];

    let mut curr = min_uv;
    let grid_step = step.max(1.0e-15);
    loop {
        let value = f(curr).unwrap_or(f64::MAX);
        // `GetWorstParticle` = first `max_element` (`math_PSOParticlesPool.cxx:65-68`).
        let worst_i = first_max_distance(&pool);
        if value < pool[worst_i].distance {
            pool[worst_i].position = curr;
            pool[worst_i].best_position = curr;
            pool[worst_i].distance = value;
            pool[worst_i].best_distance = value;
        }
        curr += grid_step;
        if curr > max_uv {
            break;
        }
    }

    let mut rng = Bullard::new();
    for p in &mut pool {
        p.velocity = step * (rng.next_real() - 0.5) * 2.0;
    }

    // `GetBestParticle` = first `min_element` (`math_PSOParticlesPool.cxx:58-61`).
    let best_i = first_min_distance(&pool);
    let mut best_pos = pool[best_i].position;
    let mut best_dist = pool[best_i].distance;
    let term_vel = step / 2048.0;

    // `for (aStep = 1; aStep < myNbIter; ++aStep)` (`cxx:169`).
    let mut a_step = 1i32;
    while a_step < nb_iter as i32 {
        let mut min_vel = f64::MAX;
        for p in &mut pool {
            let ksi1 = rng.next_real();
            let ksi2 = rng.next_real();
            const RETENT: f64 = 0.72900;
            const PERSON: f64 = 1.49445;
            const SOCIAL: f64 = 1.49445;
            p.velocity = p.velocity * RETENT
                + (p.best_position - p.position) * (PERSON * ksi1)
                + (best_pos - p.position) * (SOCIAL * ksi2);
            p.position += p.velocity;
            p.position = p.position.clamp(min_uv, max_uv);
            min_vel = min_vel.min(p.velocity.abs());
            p.distance = f(p.position).unwrap_or(f64::MAX);
            if p.distance < p.best_distance {
                p.best_distance = p.distance;
                p.best_position = p.position;
                if p.distance < best_dist {
                    best_dist = p.distance;
                    best_pos = p.position;
                }
            }
        }
        if min_vel <= term_vel {
            const MIN_STEPS: i32 = 16;
            if a_step > MIN_STEPS {
                break;
            }
            for p in &mut pool {
                let ksi = rng.next_real();
                if p.position == min_uv || p.position == max_uv {
                    p.velocity = if p.position == min_uv {
                        step * ksi
                    } else {
                        -step * ksi
                    };
                } else {
                    p.velocity = step * (ksi - 0.5) * 2.0;
                }
            }
        }
        a_step += 1;
    }
    (best_dist, best_pos)
}

fn first_max_distance(pool: &[Particle]) -> usize {
    let mut i = 0;
    for k in 1..pool.len() {
        if pool[k].distance > pool[i].distance {
            i = k;
        }
    }
    i
}

fn first_min_distance(pool: &[Particle]) -> usize {
    let mut i = 0;
    for k in 1..pool.len() {
        if pool[k].distance < pool[i].distance {
            i = k;
        }
    }
    i
}
