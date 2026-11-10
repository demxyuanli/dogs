//! `IntAna_IntLinTorus` (`line_torus_intersect`).

use super::super::prelude::*;

use occt_core::math_direct_poly_roots::DirectPolynomialRoots;

// ---------------------------------------------------------------------------
// IntAna_IntLinTorus
// ---------------------------------------------------------------------------

/// Result of [`line_torus_intersect`] — port of `IntAna_IntLinTorus`
/// (`IntAna_IntLinTorus.hxx:29-66`), which carries a `done` flag beside the
/// points and throws `StdFail_NotDone` from every accessor when it is false
/// (`IntAna_IntLinTorus.lxx:23-44`).
#[derive(Debug, Clone)]
pub struct LineTorusIntersection {
    done: bool,
    points: Vec<GpPnt>,
}

impl LineTorusIntersection {
    /// `IsDone` (`IntAna_IntLinTorus.lxx:18-21`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// The `NbPoints` / `Value(i)` payload (`IntAna_IntLinTorus.lxx:23-44`),
    /// in the order the quartic roots were accepted (`cxx:109-115`); empty when
    /// `!is_done()`.
    pub fn points(&self) -> &[GpPnt] {
        &self.points
    }
}

/// Intersection of a line with a torus. Port of `IntAna_IntLinTorus::Perform`
/// (`IntAna_IntLinTorus.cxx:43-133`): a quartic in the line parameter in the
/// torus reference frame, with each root verified by re-evaluation on the
/// torus. `done` is the OCCT flag: `false` when the quartic solver declined
/// (`cxx:126-131`) or when it returned roots that were all rejected on
/// back-substitution (`cxx:116-120`).
pub fn line_torus_intersect(l: &GpLin, t: &GpTorus) -> LineTorusIntersection {
    let pl = l.location();
    let dl = l.direction();
    let tor_loc = t.location();
    // Reparametrize so the line location is nearest the torus location.
    let param_of_new_pl = GpVec::from_pnts(&pl, &tor_loc).dot(&GpVec::from_xyz(dl.xyz()));
    let new_pl = pl.translated_vec(&GpVec::from_xyz(dl.xyz()).multiplied_scalar(param_of_new_pl));

    // Express the line in the torus reference frame.
    let pos = t.position();
    let (xd, yd, zd) = (
        GpVec::from_xyz(pos.x_direction().xyz()),
        GpVec::from_xyz(pos.y_direction().xyz()),
        GpVec::from_xyz(pos.direction().xyz()),
    );
    let v = GpVec::from_pnts(&pos.location(), &new_pl);
    let (x0, y0, z0) = (v.dot(&xd), v.dot(&yd), v.dot(&zd));
    let dv = GpVec::from_xyz(dl.xyz());
    let (x1, y1, z1) = (dv.dot(&xd), dv.dot(&yd), dv.dot(&zd));

    let r = t.major_radius();
    let r2 = r * r;
    let rr = t.minor_radius();
    let rr2 = rr * rr;

    let a = x1 * x1 + y1 * y1 + z1 * z1;
    let b = 2.0 * (x1 * x0 + y1 * y0 + z1 * z0);
    let c = x0 * x0 + y0 * y0 + z0 * z0 - (r2 + rr2);

    let a4 = a * a;
    let a3 = 2.0 * a * b;
    let a2 = 2.0 * a * c + 4.0 * r2 * z1 * z1 + b * b;
    let a1 = 2.0 * b * c + 8.0 * r2 * z1 * z0;
    let a0 = c * c + 4.0 * r2 * (z0 * z0 - rr2);

    let mut out = Vec::new();
    // `IntAna_IntLinTorus.cxx:88-133`. The `Value(1 .. NbSolutions)` list is
    // walked through the ported solver itself, not through the
    // `quartic_roots` adapter, because that adapter folds `!IsDone()` into an
    // empty root list while OCCT distinguishes the two: a declined solver
    // (`cxx:126-131`) still sets `done = false` where a completed solver with
    // zero roots (`cxx:121-125`) sets `done = true`.
    let mdpr = DirectPolynomialRoots::new5(a4, a3, a2, a1, a0);
    if std::env::var("SFS_PROBE").is_ok() {
        eprintln!(
            "LTPROBE a4={a4:.6} a3={a3:.6} a2={a2:.6} a1={a1:.6} a0={a0:.6} newpl={param_of_new_pl:.9}"
        );
    }
    if !mdpr.is_done() {
        return LineTorusIntersection { done: false, points: Vec::new() };
    }
    let n = mdpr.nb_solutions();
    if std::env::var("SFS_PROBE").is_ok() {
        eprint!("LTPROBE n={n}");
        for i in 1..=n {
            eprint!(" root{i}={:.9}", mdpr.value(i));
        }
        eprintln!();
    }
    let mut nb_bad_sol = 0;
    for i in 1..=n {
        let tt = mdpr.value(i) + param_of_new_pl;
        let p = clib::line_value(l, tt);
        // `IntAna_IntLinTorus.cxx:98-106`: re-parameterise the candidate on the
        // torus (`ElSLib::Parameters`) and re-evaluate it (`ElSLib::Value`); a
        // root is kept only when the **square distance** between the line point
        // and the torus point is `<= 1e-10` (OCCT counts the others as bad
        // solutions). The previous body tested an invented implicit-equation
        // residual `|(rho-R)^2 + z^2 - r^2| < 1e-7` and deduplicated points with
        // a further `1e-7` (audit A15) — OCCT stores every valid root as is.
        let (u, v) = slib::torus_parameters(&t.position(), r, rr, &p);
        let p_sol_t = slib::torus_value(t, u, v);
        if p_sol_t.square_distance(&p) > 1.0e-10 {
            nb_bad_sol += 1;
            continue;
        }
        out.push(p);
    }
    if n > 0 && out.is_empty() && nb_bad_sol == n {
        // `cxx:116-120`: every root failed the back-substitution check, so the
        // intersection is reported as not done rather than as an empty set.
        return LineTorusIntersection { done: false, points: Vec::new() };
    }
    // `cxx:104-108` stores each accepted root at `thePoint[nbsolvalid]`, i.e. in
    // root-index order; OCCT never reorders the point list here. An earlier
    // `out.sort_by(|a, b| a.x() ...)` had no counterpart in the `.cxx` (audit
    // `occt-geom.md` section 9 item (b)); the consumers that need an order sort
    // on their own (`intcurvesurface/solvers.rs:210` by curve parameter).
    LineTorusIntersection { done: true, points: out }
}
