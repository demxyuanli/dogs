//! The `Law` static helpers: `Law::MixBnd`, `Law::MixTgt`, `Law::Scale`,
//! `Law::ScaleCub` and `Law::Reparametrize`.
//!
//! Source: `Law.hxx`/`Law.cxx:29-394` (OCCT 8.0.0).
//!
//! Translation notes:
//! - The `Knots` / `Mults` arrays of `MixBnd` / `MixTgt` are indexed 1-based in
//!   OCCT; the port keeps that convention so `Knots(Index)`, `Mults(1)` and the
//!   `Knots.Lower()` / `Knots.Upper()` reads map one-to-one.
//! - `BSplCLib::Interpolate` answers an `inversion_problem` code that both
//!   `Law::MixBnd` (`cxx:75-83`) and `Law::MixTgt` (`cxx:152`) ignore: the
//!   partially inverted poles are still used. The port ignores the `Err` the
//!   same way.
//! - `Adaptor3d_Curve` becomes `&dyn crate::curve::Curve`; only `Value`,
//!   `FirstParameter` and `LastParameter` are used.

use occt_core::bspl::banded_interp::{interpolate_contact, schoenberg_points};
use occt_core::precision::CONFUSION;

use crate::curve::Curve;
use crate::law::bspline::{LawBSpline, LawBSpFunc};
use crate::law::function::{LawFunction, LawLinear};
use crate::law::interpolate::LawInterpolate;

/// `Law::MixBnd(Degree, Knots, Mults, Lin)` (`Law.cxx:47-84`): the poles of the
/// blending basis obtained by interpolating `Lin` at the Schoenberg points of
/// the repeated knot vector.
pub fn law_mix_bnd_poles(
    degree: i32,
    knots: &[f64],
    mults: &[i32],
    lin: &LawLinear,
) -> Vec<f64> {
    // `cxx:56-64`: build the repeated knot vector.
    let mut nbfk = 0usize;
    for &m in mults {
        nbfk += m.max(0) as usize;
    }
    let mut fk = Vec::with_capacity(nbfk);
    for (i, &m) in mults.iter().enumerate() {
        for _ in 1..=m {
            fk.push(knots[i]);
        }
    }
    // `cxx:65-66`.
    let nbpol = nbfk as i32 - degree - 1;
    let nbpol_usize = nbpol.max(0) as usize;
    let par = schoenberg_points(degree as usize, &fk, nbpol_usize);
    // `cxx:67-73`.
    let mut pol: Vec<f64> = par.iter().map(|&p| lin.value(p)).collect();
    let contact = vec![0i32; nbpol_usize];
    let _inversion = interpolate_contact(
        degree as usize,
        &fk,
        &par,
        &contact,
        &mut pol,
        1,
    );
    // `cxx:79-83`: the second and second-to-last poles are copied from their
    // neighbours when the basis is long enough.
    if nbpol >= 4 {
        pol[1] = pol[0];
        pol[nbpol_usize - 2] = pol[nbpol_usize - 1];
    }
    pol
}

/// `Law::MixBnd(Lin)` (`Law.cxx:29-45`): cubic blending on the four knots
/// `f, 0.75 f + 0.25 l, 0.25 f + 0.75 l, l` with multiplicities `4, 1, 1, 4`.
pub fn law_mix_bnd(lin: &LawLinear) -> Result<LawBSpFunc, String> {
    let (f, l) = lin.bounds();
    let knots = [f, 0.75 * f + 0.25 * l, 0.25 * f + 0.75 * l, l];
    let mults = [4i32, 1, 1, 4];
    let pol = law_mix_bnd_poles(3, &knots, &mults, lin);
    let bs = LawBSpline::new(pol, knots.to_vec(), mults.to_vec(), 3, false)?;
    let mut bsf = LawBSpFunc::new();
    bsf.set_curve(bs);
    Ok(bsf)
}

/// `eval1` (`Law.cxx:86-121`): the tangent blending basis, zero on one side of
/// `piv` and cubic on the other.
fn mix_tgt_eval(p: f64, first: f64, last: f64, piv: f64, nulr: bool) -> f64 {
    if (nulr && p >= piv) || (!nulr && p <= piv) {
        return 0.0;
    }
    if nulr {
        let mut a = piv - first;
        a *= a;
        a = 1.0 / a;
        let mut b = p - first;
        a *= b;
        b = piv - p;
        a *= b;
        a *= b;
        a
    } else {
        let mut a = last - piv;
        a *= a;
        a = 1.0 / a;
        let mut b = last - p;
        a *= b;
        b = p - piv;
        a *= b;
        a *= b;
        a
    }
}

/// `Law::MixTgt(Degree, Knots, Mults, NulOnTheRight, Index)`
/// (`Law.cxx:123-156`): the poles of the tangent blending basis, interpolated
/// from [`mix_tgt_eval`] at the Schoenberg points.
pub fn law_mix_tgt(
    degree: i32,
    knots: &[f64],
    mults: &[i32],
    nul_on_the_right: bool,
    index: usize,
) -> Vec<f64> {
    let first = knots[0];
    let last = knots[knots.len() - 1];
    // `cxx:130`: `Knots(Index)`, 1-based.
    let piv = knots[index - 1];
    let mut nbfk = 0usize;
    for &m in mults {
        nbfk += m.max(0) as usize;
    }
    let mut fk = Vec::with_capacity(nbfk);
    for (i, &m) in mults.iter().enumerate() {
        for _ in 1..=m {
            fk.push(knots[i]);
        }
    }
    let nbpol = nbfk as i32 - degree - 1;
    let nbpol_usize = nbpol.max(0) as usize;
    let par = schoenberg_points(degree as usize, &fk, nbpol_usize);
    let mut pol: Vec<f64> = par
        .iter()
        .map(|&p| mix_tgt_eval(p, first, last, piv, nul_on_the_right))
        .collect();
    let contact = vec![0i32; nbpol_usize];
    let _inversion = interpolate_contact(degree as usize, &fk, &par, &contact, &mut pol, 1);
    pol
}

/// `Law::Reparametrize` (`Law.cxx:158-268`): the law $\lambda(t)$ that maps the
/// chord-length abscissa of `curve` onto `[First, Last]`.
///
/// `rev` reverses the curve direction and negates the given end tangents;
/// `has_df` / `has_dl` constrain the two end derivatives.
#[allow(clippy::too_many_arguments)]
pub fn law_reparametrize(
    curve: &dyn Curve,
    first: f64,
    last: f64,
    has_df: bool,
    has_dl: bool,
    d_first: f64,
    d_last: f64,
    rev: bool,
    nb_points: usize,
) -> Result<LawBSpline, String> {
    // `cxx:171-183`.
    let mut dd_first = d_first;
    let mut dd_last = d_last;
    if has_df && rev {
        dd_first = -d_first;
    }
    if has_dl && rev {
        dd_last = -d_last;
    }
    // `cxx:184-186`: 1-based work arrays.
    let n2 = 2 * nb_points;
    let mut cumdist = vec![0.0; n2 + 1];
    let mut ucourbe = vec![0.0; n2 + 1];
    let u1 = curve.first_parameter();
    let u2 = curve.last_parameter();
    let mut u;
    let du;
    let mut p1;
    if !rev {
        // `cxx:188-193`.
        p1 = curve.value(u1);
        u = u1;
        du = (u2 - u1) / (n2 as f64 - 1.0);
    } else {
        // `cxx:194-199`.
        p1 = curve.value(u2);
        u = u2;
        du = (u1 - u2) / (n2 as f64 - 1.0);
    }
    let mut length = 0.0;
    for i in 1..=n2 {
        // `cxx:200-207`.
        let p2 = curve.value(u);
        length += p2.distance(&p1);
        cumdist[i] = length;
        ucourbe[i] = u;
        u += du;
        p1 = p2;
    }
    // `cxx:208-215`.
    if rev {
        ucourbe[n2] = u1;
    } else {
        ucourbe[n2] = u2;
    }

    let mut point = vec![0.0; nb_points + 1];
    let mut param = vec![0.0; nb_points + 1];
    // `cxx:220-229`.
    let d_corde = length / (nb_points as f64 - 1.0);
    let mut corde = d_corde;
    let mut index = 1usize;
    let fac = 1.0 / (nb_points as f64 - 1.0);
    point[1] = ucourbe[1];
    param[1] = first;
    point[nb_points] = ucourbe[n2];
    param[nb_points] = last;
    // `cxx:231-246`.
    for i in 2..nb_points {
        while index <= n2 && cumdist[index] < corde {
            index += 1;
        }
        if index > n2 {
            return Err("Standard_OutOfRange: Law::Reparametrize".into());
        }
        let alpha = (corde - cumdist[index - 1]) / (cumdist[index] - cumdist[index - 1]);
        let uu = ucourbe[index - 1] + alpha * (ucourbe[index] - ucourbe[index - 1]);
        point[i] = uu;
        param[i] = ((nb_points - i) as f64 * first + (i - 1) as f64 * last) * fac;
        corde = i as f64 * d_corde;
    }
    // `cxx:247`.
    let mut inter = LawInterpolate::new_with_parameters(
        point[1..=nb_points].to_vec(),
        param[1..=nb_points].to_vec(),
        false,
        1.0e-9,
    )?;
    // `cxx:248-261`.
    if has_df || has_dl {
        let mut tgt = vec![0.0; nb_points];
        let mut flag = vec![false; nb_points];
        if has_df {
            flag[0] = true;
            tgt[0] = dd_first;
        }
        if has_dl {
            flag[nb_points - 1] = true;
            tgt[nb_points - 1] = dd_last;
        }
        inter.load(&tgt, &flag)?;
    }
    // `cxx:262-267`.
    inter.perform()?;
    if !inter.is_done() {
        return Err("Standard_Failure: Law::Reparametrize echec interpolation".into());
    }
    inter.curve().cloned()
}

/// `eval2` (`Law.cxx:270-294`): the piecewise `Scale` blending value.
#[allow(clippy::too_many_arguments)]
fn scale_eval(
    p: f64,
    mil: f64,
    has_first: bool,
    has_last: bool,
    bs1: Option<&LawBSpline>,
    bs2: Option<&LawBSpline>,
) -> f64 {
    if has_first && p < mil {
        match bs1 {
            Some(bs) => bs.value(p),
            None => f64::NAN,
        }
    } else if has_last && p > mil {
        match bs2 {
            Some(bs) => bs.value(p),
            None => f64::NAN,
        }
    } else {
        1.0
    }
}

/// `Law::Scale` (`Law.cxx:296-352`): the cubic law that stays at `VFirst` up to
/// the middle of `[First, Last]`, is `1` there, and reaches `VLast` at `Last`.
pub fn law_scale(
    first: f64,
    last: f64,
    has_f: bool,
    has_l: bool,
    v_first: f64,
    v_last: f64,
) -> Result<LawBSpline, String> {
    let milieu = 0.5 * (first + last);
    // `cxx:305-315`: 1-based `knot` / `mult`.
    let knot = [first, milieu, last];
    let mult = [4i32, 2, 4];
    // `cxx:311-315`: the repeated knot vector fed to the Schoenberg build.
    let fknot = [
        first, first, first, first, milieu, milieu, last, last, last, last,
    ];
    let mut bs1 = None;
    let mut bs2 = None;
    // `cxx:317-339`: 1-based `pbs` / `kbs` / `mbs`.
    if has_f {
        let pbs = [v_first, v_first, 1.0, 1.0];
        let kbs = [first, milieu];
        let mbs = [4i32, 4];
        bs1 = Some(LawBSpline::new(pbs.to_vec(), kbs.to_vec(), mbs.to_vec(), 3, false)?);
    }
    if has_l {
        let pbs = [1.0, 1.0, v_last, v_last];
        let kbs = [milieu, last];
        let mbs = [4i32, 4];
        bs2 = Some(LawBSpline::new(pbs.to_vec(), kbs.to_vec(), mbs.to_vec(), 3, false)?);
    }
    // `cxx:340-349`.
    let par = schoenberg_points(3, &fknot, 6);
    let mut pol: Vec<f64> = par
        .iter()
        .map(|&p| scale_eval(p, milieu, has_f, has_l, bs1.as_ref(), bs2.as_ref()))
        .collect();
    let contact = vec![0i32; 6];
    let _inversion = interpolate_contact(3, &fknot, &par, &contact, &mut pol, 1);
    // `cxx:350-351`.
    LawBSpline::new(pol, knot.to_vec(), mult.to_vec(), 3, false)
}

/// `Law::ScaleCub` (`Law.cxx:354-394`): the three-span cubic law with poles
/// `(VFirst, VFirst, 1, VLast, VLast)`.
pub fn law_scale_cub(
    first: f64,
    last: f64,
    has_f: bool,
    has_l: bool,
    v_first: f64,
    v_last: f64,
) -> Result<LawBSpline, String> {
    let milieu = 0.5 * (first + last);
    // `cxx:364-372`: `pol(1..=5)`, `knot` / `mult` 1-based.
    let mut pol = [0.0f64; 5];
    let knot = [first, milieu, last];
    let mult = [4i32, 1, 4];
    // `cxx:375-390`.
    if has_f {
        pol[0] = v_first;
        pol[1] = v_first;
    } else {
        pol[0] = 1.0;
        pol[1] = 1.0;
    }
    if has_l {
        pol[3] = v_last;
        pol[4] = v_last;
    } else {
        pol[3] = 1.0;
        pol[4] = 1.0;
    }
    pol[2] = 1.0;
    // `cxx:392`.
    LawBSpline::new(pol.to_vec(), knot.to_vec(), mult.to_vec(), 3, false)
}

/// `Law_Interpolate` with `Precision::Confusion()` as used by `Law_Interpol`
/// (`Law_Interpol.cxx:63`); re-exported so callers do not restate the
/// tolerance.
pub const LAW_INTERPOL_TOLERANCE: f64 = CONFUSION;
