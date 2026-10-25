//! The file-static helpers of `GeomFill_ConstrainedFilling.cxx`.
//!
//! Source: `GeomFill_ConstrainedFilling.cxx:55-86` (`inqadd`), `:90-101`
//! (`mklin`), `:103-214` (`sortbounds`), `:216-288` (`coonscnd`) and
//! `:290-332` (`killcorners`) (OCCT 8.0.0).
//!
//! Translation notes:
//! - `occ::handle<GeomFill_Boundary>*` becomes `&mut [Arc<dyn Boundary>]`; the
//!   `bound[i] <-> bound[j]` handle swaps are `slice::swap` on the same array,
//!   so the corner bookkeeping that refers to positions keeps its meaning.
//! - `bound[i]->Points(pf, pl)` / `Bounds(fi, li)` write through out-parameters
//!   in OCCT; the port takes the `(first, last)` pairs the trait returns.
//! - `Norm` / `D1Norm` can raise `Standard_Failure` on a boundary without
//!   normals; the trait returns `Result` and the helpers propagate it. OCCT
//!   only reaches those calls behind `HasNormals()`.

use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::gp::GpXyz;

use crate::geom_fill::boundary::Boundary;
use crate::geom_fill::corner_state::CornerState;
use crate::geom_fill::tgt_field::TgtField;
use crate::law::function::{LawFunction, LawLinear};
use crate::law::services::law_scale_cub;

/// `inqadd` (`cxx:57-86`): the knots to insert so that the blending function's
/// domain matches the requested `d1` / `d2`, written into `k`/`m` from index 0.
/// Returns `NbAdd`.
pub fn inqadd(d1: f64, d2: f64, k: &mut [f64; 2], m: &mut [i32; 2], deg: i32, tolk: f64) -> usize {
    // `cxx:61`.
    m[0] = deg - 2;
    m[1] = deg - 2;
    if d1 != 1.0 && d2 != 1.0 {
        if (d1 + d2 - 1.0).abs() < tolk {
            // `cxx:67-70`.
            k[0] = 0.5 * (d1 + 1.0 - d2);
            1
        } else {
            // `cxx:72-77`.
            k[0] = d1.min(1.0 - d2);
            k[1] = d1.max(1.0 - d2);
            2
        }
    } else if d1 != 1.0 {
        // `cxx:79-82`.
        k[0] = d1;
        1
    } else if d2 != 1.0 {
        // `cxx:84-87`.
        k[0] = d2;
        1
    } else {
        // `cxx:88`.
        0
    }
}

/// `mklin` (`cxx:90-101`): the law seen as a straight line. A `Law_Linear` is
/// returned as is (`occ::down_cast`); anything else is replaced by the line
/// through its end values.
pub fn mklin(func: &dyn LawFunction) -> LawLinear {
    if let Some(lin) = func.as_linear() {
        return *lin;
    }
    // `cxx:94-99`.
    let mut fu = LawLinear::new();
    let (d, f) = func.bounds();
    fu.set(d, func.value(d), f, func.value(f));
    fu
}

/// `sortbounds` (`cxx:103-214`): order the boundaries head to tail, flag the
/// ones to reverse and record the corner gaps / tangent and normal angles.
pub fn sort_bounds(
    nb: usize,
    bound: &mut [Arc<dyn Boundary>],
    rev: &mut [bool],
    stat: &mut [CornerState],
) -> Result<(), String> {
    // `cxx:115`.
    rev[0] = false;
    for i in 0..nb - 1 {
        // `cxx:119-126`.
        let (mut pf, mut pl) = bound[i].points();
        if rev[i] {
            std::mem::swap(&mut pf, &mut pl);
        }
        for j in i + 1..nb {
            // `cxx:129-152`.
            let (qf, ql) = bound[j].points();
            let df = qf.distance(&pl);
            let dl = ql.distance(&pl);
            if df < dl {
                if df < stat[i + 1].gap() {
                    bound.swap(i + 1, j);
                    stat[i + 1].set_gap(df);
                    rev[i + 1] = false;
                }
            } else if dl < stat[i + 1].gap() {
                bound.swap(i + 1, j);
                stat[i + 1].set_gap(dl);
                rev[i + 1] = true;
            }
        }
    }
    // `cxx:158-167`.
    let (mut pf, mut pl) = bound[nb - 1].points();
    if rev[nb - 1] {
        std::mem::swap(&mut pf, &mut pl);
    }
    let (qf, _) = bound[0].points();
    stat[0].set_gap(pl.distance(&qf));

    // `cxx:169-213`: the corner angles.
    for i in 0..nb {
        let next = (i + 1) % nb;
        let (fi, li) = bound[i].bounds();
        let li = if rev[i] { fi } else { li };
        let (_, mut tgi) = bound[i].d1(li);
        if rev[i] {
            tgi.reverse();
        }
        let (fn_, ln) = bound[next].bounds();
        let fn_ = if rev[next] { ln } else { fn_ };
        let (_, mut tgn) = bound[next].d1(fn_);
        if rev[next] {
            tgn.reverse();
        }
        let ang = PI - tgi.angle(&tgn);
        stat[next].set_tgt_ang(ang);
        if bound[i].has_normals() && bound[next].has_normals() {
            stat[next].constraint();
            let nori = bound[i].norm(li)?;
            let norn = bound[next].norm(fn_)?;
            stat[next].set_nor_ang(nori.angle(&norn));
        }
    }
    Ok(())
}

/// `coonscnd` (`cxx:216-288`): the Coons compatibility test on every
/// constrained corner. (`tga` is an unnamed parameter in OCCT, kept out.)
pub fn coonscnd(
    nb: usize,
    bound: &[Arc<dyn Boundary>],
    rev: &[bool],
    stat: &mut [CornerState],
    mintg: &[f64],
) -> Result<(), String> {
    let fact_normalization = 100.0;
    for i in 0..nb {
        if stat[i].has_constraint() {
            let ip = (i + nb - 1) % nb;
            // `cxx:231-232`.
            let mut tolang = bound[ip].tolang().min(bound[i].tolang());
            let mut an = stat[i].nor_ang();
            let mut twist = false;
            if an >= 0.5 * PI {
                twist = true;
                an = PI - an;
            }
            if an > tolang {
                // `cxx:239-242`.
                stat[i].do_kill(0.0);
            } else {
                // `cxx:244-284`.
                let fact = 0.5 * 27.0 / 4.0;
                tolang *= mintg[ip].min(mintg[i]) * fact * fact_normalization;
                let (fp, lp) = bound[ip].bounds();
                let lp = if rev[ip] { fp } else { lp };
                let (_, tgp) = bound[ip].d1(lp);
                let (_, dnorp) = bound[ip].d1_norm(lp)?;
                let (fi, li) = bound[i].bounds();
                let fi = if rev[i] { li } else { fi };
                let (_, tgi) = bound[i].d1(fi);
                let (_, dnori) = bound[i].d1_norm(fi)?;
                let scal1 = tgp.dot(&dnori);
                let mut scal2 = tgi.dot(&dnorp);
                if !twist {
                    scal2 *= -1.0;
                }
                let scal1 = (scal1 + scal2).abs();
                if scal1 > tolang {
                    // `cxx:279-281`.
                    let killfactor = tolang / scal1;
                    stat[i].do_kill(killfactor);
                }
            }
        }
    }
    Ok(())
}

/// `killcorners` (`cxx:290-332`): reparametrize the boundaries whose corner was
/// killed and rescale the tangency fields.
pub fn kill_corners(
    nb: usize,
    bound: &[Arc<dyn Boundary>],
    rev: &[bool],
    nrev: &[bool],
    stat: &[CornerState],
    tga: &[Option<Arc<dyn TgtField>>],
) -> Result<(), String> {
    for i in 0..nb {
        let inext = (i + 1) % nb;
        let (fnul, fscal, lnul, lscal) = if !nrev[i] {
            // `cxx:306-310`.
            let (fnul, fscal) = stat[i].is_to_kill();
            let (lnul, lscal) = stat[inext].is_to_kill();
            (fnul, fscal, lnul, lscal)
        } else {
            // `cxx:311-315`.
            let (lnul, lscal) = stat[i].is_to_kill();
            let (fnul, fscal) = stat[inext].is_to_kill();
            (fnul, fscal, lnul, lscal)
        };
        if fnul || lnul {
            // `cxx:320`.
            bound[i].reparametrize(0.0, 1.0, fnul, lnul, fscal, lscal, rev[i])?;
            // `cxx:324-328`.
            let scalable = tga[i].as_ref().map_or(false, |t| t.is_scalable());
            if bound[i].has_normals() && scalable {
                let bs = law_scale_cub(0.0, 1.0, fnul, lnul, fscal, lscal)?;
                if let Some(t) = tga[i].as_ref() {
                    t.scale(&bs);
                }
            }
        }
    }
    Ok(())
}

/// `gp_XYZ::Multiplied(-1) + Added + Multiplied(coef)`: one end term of the
/// corner corrections `v[4]` of `PerformS1` (`cxx:1164-1169`).
pub fn end_term(a: &GpXyz, b: &GpXyz, coef: f64) -> GpXyz {
    let v = a.multiplied(-1.0);
    v.added(b).multiplied(coef)
}
