//! `Adaptor3d_TopolTool` sampling: `ComputeSamplePoints`, `Analyse`,
//! `SamplePnts`, `BSplSamplePnts` and `IsUniformSampling`.
//!
//! Source: `Adaptor3d_TopolTool.cxx:638-749` (`Analyse`), `:750-898`
//! (`ComputeSamplePoints`), `:1019-1130` (`SamplePnts`) and `:1133-1646`
//! (`BSplSamplePnts`).

use occt_geom::approx_same_parameter::u_resolution;
use occt_geom::approx_same_parameter::v_resolution;
use occt_geom::Surface;

use crate::geom_bnd_lib_surface3d::{surface_kind, GeomAbsSurfaceType};

use crate::geom_int::finite_uv;

/// `myMinPnts` (`Adaptor3d_TopolTool.cxx:635`).
const MY_MIN_PNTS: usize = 4;

/// `aMaxPnts` (`Adaptor3d_TopolTool.cxx:1141`).
const A_MAX_PNTS: usize = 1001;

/// The two parameter arrays `SamplePnts` leaves in `myUPars` / `myVPars`.
#[derive(Debug, Clone, Default)]
pub struct SampleParameters {
    pub u_pars: Vec<f64>,
    pub v_pars: Vec<f64>,
}

impl SampleParameters {
    pub fn nb_u(&self) -> i32 {
        self.u_pars.len() as i32
    }

    pub fn nb_v(&self) -> i32 {
        self.v_pars.len() as i32
    }
}

/// `Adaptor3d_TopolTool::Analyse` (`cxx:638-749`).
///
/// `poles` is the U-major pole grid (`array2.Value(i, j)` is `poles[i-1][j-1]`).
fn analyse(poles: &[Vec<occt_core::gp::GpPnt>], nbup: usize, nbvp: usize) -> (i32, i32) {
    let at = |i: usize, j: usize| -> occt_core::gp::GpPnt {
        poles.get(i - 1).and_then(|row| row.get(j - 1)).copied().unwrap_or_default()
    };
    let second_diff = |a: &occt_core::gp::GpPnt,
                       b: &occt_core::gp::GpPnt,
                       c: &occt_core::gp::GpPnt|
     -> occt_core::gp::GpVec {
        occt_core::gp::GpVec::new(
            c.x() - b.x() - b.x() + a.x(),
            c.y() - b.y() - b.y() + a.y(),
            c.z() - b.z() - b.z() + a.z(),
        )
    };

    let mut sh: i32 = 1;
    let mut nbch = 0usize;

    if nbvp > 2 {
        for i in 2..nbup {
            let (a, b, c) = (at(i, 1), at(i, 2), at(i, 3));
            let mut vi = second_diff(&a, &b, &c);
            let mut loc_nbch = 0usize;
            for j in 3..nbvp {
                let (a1, b1, c1) = (at(i, j - 1), at(i, j), at(i, j + 1));
                let vip1 = second_diff(&a1, &b1, &c1);
                let pd = vi.dot(&vip1);
                vi = vip1;
                if pd > 1.0e-7 || pd < -1.0e-7 {
                    if pd > 0.0 {
                        if sh == -1 {
                            sh = 1;
                            loc_nbch += 1;
                        }
                    } else if sh == 1 {
                        sh = -1;
                        loc_nbch += 1;
                    }
                }
            }
            if loc_nbch > nbch {
                nbch = loc_nbch;
            }
        }
    }
    let nbsv = nbch + 5;

    nbch = 0;
    if nbup > 2 {
        for j in 2..nbvp {
            let (a, b, c) = (at(1, j), at(2, j), at(3, j));
            let mut vi = second_diff(&a, &b, &c);
            let mut loc_nbch = 0usize;
            for i in 3..nbup {
                let (a1, b1, c1) = (at(i - 1, j), at(i, j), at(i + 1, j));
                let vip1 = second_diff(&a1, &b1, &c1);
                let pd = vi.dot(&vip1);
                vi = vip1;
                if pd > 1.0e-7 || pd < -1.0e-7 {
                    if pd > 0.0 {
                        if sh == -1 {
                            sh = 1;
                            loc_nbch += 1;
                        }
                    } else if sh == 1 {
                        sh = -1;
                        loc_nbch += 1;
                    }
                }
            }
            if loc_nbch > nbch {
                nbch = loc_nbch;
            }
        }
    }
    let nbsu = nbch + 5;

    (nbsu as i32, nbsv as i32)
}

/// `Adaptor3d_TopolTool::ComputeSamplePoints` (`cxx:750-898`).
pub(crate) fn compute_sample_points(s: &dyn Surface) -> (i32, i32) {
    let (uinf, usup) = finite_uv(s.u_range().0, s.u_range().1);
    let (vinf, vsup) = finite_uv(s.v_range().0, s.v_range().1);

    let kind = surface_kind(s);
    let (mut nbsu, mut nbsv) = match kind {
        GeomAbsSurfaceType::Plane => (2i32, 2i32),
        GeomAbsSurfaceType::BezierSurface => {
            let (nbu, nbv) = bezier_pole_counts(s);
            (3 + nbu, 3 + nbv)
        }
        GeomAbsSurfaceType::BSplineSurface => {
            let (nbu, nbv) = bspline_sample_counts(s);
            (nbu.max(4), nbv.max(4))
        }
        GeomAbsSurfaceType::Cylinder
        | GeomAbsSurfaceType::Cone
        | GeomAbsSurfaceType::Sphere
        | GeomAbsSurfaceType::Torus
        | GeomAbsSurfaceType::SurfaceOfRevolution
        | GeomAbsSurfaceType::SurfaceOfExtrusion => (15, 15),
        _ => (10, 10),
    };

    if nbsu < 6 {
        nbsu = 6;
    }
    if nbsv < 6 {
        nbsv = 6;
    }

    if kind == GeomAbsSurfaceType::BSplineSurface {
        if nbsu > 8 || nbsv > 8 {
            if let Some(bs) = s.osculating_bspline() {
                let (nu, nv) = analyse(&bs.poles, bs.poles.len(), bs.poles.first().map_or(0, Vec::len));
                nbsu = nu;
                nbsv = nv;
            }
        }
        // Check anisotropy (`cxx:867-881`).
        let u_len = length_or_zero(usup - uinf) / u_resolution(s, 1.0);
        let v_len = length_or_zero(vsup - vinf) / v_resolution(s, 1.0);
        if v_len != 0.0 {
            let ratio = u_len / v_len;
            if ratio >= 10.0 {
                nbsu *= 2;
                nbsu = nbsu.min(50);
            } else if ratio <= 0.1 {
                nbsv *= 2;
                nbsv = nbsv.min(50);
            }
        }
    } else if kind == GeomAbsSurfaceType::BezierSurface && (nbsu > 8 || nbsv > 8) {
        if let Some((poles, nbu, nbv)) = s.bezier_control_grid() {
            // `bezier_control_grid` is row-major `poles[i * nbv + j]`; `Analyse`
            // indexes `array2.Value(i, j)`.
            let mut grid: Vec<Vec<occt_core::gp::GpPnt>> = Vec::with_capacity(nbu);
            for i in 0..nbu {
                let mut row = Vec::with_capacity(nbv);
                for j in 0..nbv {
                    row.push(*poles.get(i * nbv + j).unwrap_or(&occt_core::gp::GpPnt::zero()));
                }
                grid.push(row);
            }
            let (nu, nv) = analyse(&grid, nbu, nbv);
            nbsu = nu;
            nbsv = nv;
        }
    }

    (nbsu, nbsv)
}

fn length_or_zero(d: f64) -> f64 {
    if d.abs() < 1.0e-300 {
        0.0
    } else {
        d
    }
}

fn bezier_pole_counts(s: &dyn Surface) -> (i32, i32) {
    if let Some((_, nu, nv)) = s.bezier_control_grid() {
        return (nu as i32, nv as i32);
    }
    if let Some(b) = s.rectangular_trimmed_basis() {
        if let Some((_, nu, nv)) = b.bezier_control_grid() {
            return (nu as i32, nv as i32);
        }
    }
    (s.nb_u_poles(), s.nb_v_poles())
}

fn bspline_sample_counts(s: &dyn Surface) -> (i32, i32) {
    match s.osculating_bspline() {
        Some(bs) => (
            bs.uknots().len() as i32 * bs.deg_u() as i32,
            bs.vknots().len() as i32 * bs.deg_v() as i32,
        ),
        None => (6, 6),
    }
}

/// `Adaptor3d_TopolTool::IsUniformSampling` (`cxx:1649-1654`).
pub(crate) fn is_uniform_sampling(s: &dyn Surface) -> bool {
    surface_kind(s) != GeomAbsSurfaceType::BSplineSurface
}

/// `Adaptor3d_TopolTool::SamplePnts(theDefl, theNUmin, theNVmin)`
/// (`cxx:1019-1130`).
pub(crate) fn sample_pnts(
    s: &dyn Surface,
    the_defl: f64,
    the_nu_min: i32,
    the_nv_min: i32,
) -> SampleParameters {
    if surface_kind(s) == GeomAbsSurfaceType::BSplineSurface {
        return bspl_sample_pnts(s, the_defl, the_nu_min, the_nv_min);
    }

    let (uinf, usup) = finite_uv(s.u_range().0, s.u_range().1);
    let (vinf, vsup) = finite_uv(s.v_range().0, s.v_range().1);

    let (nbsu, nbsv) = compute_sample_points(s);

    let mut u_pars = vec![0.0; nbsu.max(1) as usize];
    let mut v_pars = vec![0.0; nbsv.max(1) as usize];
    let (nu, nv) = (nbsu as usize, nbsv as usize);

    if nu >= 2 {
        let dt = (usup - uinf) / (nu as f64 - 1.0);
        u_pars[0] = uinf;
        u_pars[nu - 1] = usup;
        let mut t = uinf + dt;
        for p in u_pars.iter_mut().take(nu - 1).skip(1) {
            *p = t;
            t += dt;
        }
    } else if nu == 1 {
        u_pars[0] = uinf;
    }

    if nv >= 2 {
        let dt = (vsup - vinf) / (nv as f64 - 1.0);
        v_pars[0] = vinf;
        v_pars[nv - 1] = vsup;
        let mut t = vinf + dt;
        for p in v_pars.iter_mut().take(nv - 1).skip(1) {
            *p = t;
            t += dt;
        }
    } else if nv == 1 {
        v_pars[0] = vinf;
    }

    SampleParameters { u_pars, v_pars }
}

/// `BSplCLib::FirstUKnotIndex` (`BSplCLib.cxx:111-122`).
fn first_u_knot_index(degree: usize, mults: &[i32]) -> usize {
    let mut index = 1usize;
    let mut sigma = mults.first().copied().unwrap_or(1).max(0) as usize;
    while sigma <= degree {
        index += 1;
        sigma += mults.get(index - 1).copied().unwrap_or(1).max(0) as usize;
        if index > mults.len() {
            break;
        }
    }
    index
}

/// `BSplCLib::LastUKnotIndex` (`BSplCLib.cxx:126-137`).
fn last_u_knot_index(degree: usize, mults: &[i32]) -> usize {
    let mut index = mults.len();
    let mut sigma = mults.last().copied().unwrap_or(1).max(0) as usize;
    while sigma <= degree && index > 1 {
        index -= 1;
        sigma += mults.get(index - 1).copied().unwrap_or(1).max(0) as usize;
    }
    index
}

/// Flattened knot vector: `UKnot(i)` of `Geom_BSplineSurface`.
fn flat_knots(knots: &[f64], mults: &[i32]) -> Vec<f64> {
    let mut out = Vec::new();
    for (k, m) in knots.iter().zip(mults.iter()) {
        for _ in 0..(*m).max(0) {
            out.push(*k);
        }
    }
    out
}

/// `Adaptor3d_TopolTool::BSplSamplePnts` (`cxx:1133-1646`).
fn bspl_sample_pnts(
    s: &dyn Surface,
    the_defl: f64,
    the_nu_min: i32,
    the_nv_min: i32,
) -> SampleParameters {
    let Some(bs) = s.osculating_bspline() else {
        return sample_pnts_uniform_counts(s);
    };

    let (uinf, usup) = s.u_range();
    let (vinf, vsup) = s.v_range();

    let u_knots = flat_knots(bs.uknots(), bs.umults());
    let v_knots = flat_knots(bs.vknots(), bs.vmults());

    let deg_u = bs.deg_u();
    let deg_v = bs.deg_v();
    let periodic_u = bs.u_periodic();
    let periodic_v = bs.v_periodic();

    let mut ui1 = if periodic_u { 1 } else { first_u_knot_index(deg_u, bs.umults()) };
    let mut ui2 = if periodic_u { u_knots.len() } else { last_u_knot_index(deg_u, bs.umults()) };
    let mut vi1 = if periodic_v { 1 } else { first_u_knot_index(deg_v, bs.vmults()) };
    let mut vi2 = if periodic_v { v_knots.len() } else { last_u_knot_index(deg_v, bs.vmults()) };

    // `cxx:1152-1187`: shrink the knot window to the parameter window.
    for i in ui1..ui2 {
        if i >= 1 && i + 1 <= u_knots.len() && uinf >= u_knots[i - 1] && uinf < u_knots[i] {
            ui1 = i;
            break;
        }
    }
    let mut i = ui2;
    while i > ui1 {
        if i >= 2 && i <= u_knots.len() && usup <= u_knots[i - 1] && usup > u_knots[i - 2] {
            ui2 = i;
            break;
        }
        i -= 1;
    }
    for i in vi1..vi2 {
        if i >= 1 && i + 1 <= v_knots.len() && vinf >= v_knots[i - 1] && vinf < v_knots[i] {
            vi1 = i;
            break;
        }
    }
    let mut i = vi2;
    while i > vi1 {
        if i >= 2 && i <= v_knots.len() && vsup <= v_knots[i - 1] && vsup > v_knots[i - 2] {
            vi2 = i;
            break;
        }
        i -= 1;
    }

    let mut nbsu = (ui2 - ui1 + 1) as i32;
    nbsu += (nbsu - 1) * (deg_u as i32 - 1);
    let mut nbsv = (vi2 - vi1 + 1) as i32;
    nbsv += (nbsv - 1) * (deg_v as i32 - 1);

    let mut b_u_uniform = false;
    let mut b_v_uniform = false;

    if nbsu < the_nu_min || nbsv < the_nv_min {
        if nbsu < nbsv {
            let mut a_nb = (nbsv as f64 * (the_nu_min as f64) / (nbsu as f64)) as i32;
            a_nb = a_nb.min(30);
            if a_nb > nbsv {
                b_v_uniform = true;
            }
            if b_v_uniform {
                nbsv = a_nb;
            }
        } else {
            let mut a_nb = (nbsu as f64 * (the_nv_min as f64) / (nbsv as f64)) as i32;
            a_nb = a_nb.min(30);
            if a_nb > nbsu {
                b_u_uniform = true;
            }
            if b_u_uniform {
                nbsu = a_nb;
            }
        }
    }

    if nbsu < the_nu_min {
        nbsu = the_nu_min;
        b_u_uniform = true;
    } else if nbsu > A_MAX_PNTS as i32 {
        nbsu = A_MAX_PNTS as i32;
        b_u_uniform = true;
    }
    if nbsv < the_nv_min {
        nbsv = the_nv_min;
        b_v_uniform = true;
    } else if nbsv > A_MAX_PNTS as i32 {
        nbsv = A_MAX_PNTS as i32;
        b_v_uniform = true;
    }

    let nbsu = nbsu.max(1) as usize;
    let nbsv = nbsv.max(1) as usize;

    // `cxx:1243-1335`: fill the sample parameters.
    let mut an_u_pars = vec![0.0f64; nbsu];
    let mut an_v_pars = vec![0.0f64; nbsv];
    let mut an_u_flg = vec![false; nbsu];
    let mut a_v_flg = vec![false; nbsv];

    if b_u_uniform {
        uniform_fill(&mut an_u_pars, uinf, usup);
    } else {
        let nbi = deg_u;
        if nbi > 0 {
            let mut k = 0usize;
            let mut t1 = uinf;
            let mut i = ui1 + 1;
            while i <= ui2 {
                let t2 = if i == ui2 {
                    usup
                } else if i <= u_knots.len() {
                    u_knots[i - 1]
                } else {
                    usup
                };
                let dt = (t2 - t1) / nbi as f64;
                for _ in 0..nbi {
                    if k < nbsu {
                        an_u_pars[k] = t1;
                    }
                    k += 1;
                    t1 += dt;
                }
                t1 = t2;
                i += 1;
            }
            if k < nbsu {
                an_u_pars[k] = t1;
            }
        }
    }

    if b_v_uniform {
        uniform_fill(&mut an_v_pars, vinf, vsup);
    } else {
        let nbi = deg_v;
        if nbi > 0 {
            let mut k = 0usize;
            let mut t1 = vinf;
            let mut i = vi1 + 1;
            while i <= vi2 {
                let t2 = if i == vi2 {
                    vsup
                } else if i <= v_knots.len() {
                    v_knots[i - 1]
                } else {
                    vsup
                };
                let dt = (t2 - t1) / nbi as f64;
                for _ in 0..nbi {
                    if k < nbsv {
                        an_v_pars[k] = t1;
                    }
                    k += 1;
                    t1 += dt;
                }
                t1 = t2;
                i += 1;
            }
            if k < nbsv {
                an_v_pars[k] = t1;
            }
        }
    }

    // `cxx:1336-1415`: keep only the parameters the deflection analysis needs.
    let a_defl2 = (the_defl * the_defl).max(1.0e-9);
    let tol = (0.01 * a_defl2).max(1.0e-9);

    if nbsu >= 2 {
        an_u_flg[0] = true;
        an_u_flg[nbsu - 1] = true;
    }
    // `j` is declared outside the analysis loops in OCCT (`cxx:1145`) and the
    // `myNbSamplesU == 2` branch below reuses the value the loop left there.
    let mut u_j_last = 0usize;
    for i in 0..nbsv {
        let t1 = an_v_pars[i];
        let mut j = 0usize;
        let mut b_cont = true;
        while j + 1 < nbsu && b_cont {
            if an_u_flg[j + 1] {
                j += 1;
                continue;
            }
            let p1 = s.d0(an_u_pars[j], t1);
            let mut k = j + 2;
            while k < nbsu {
                let p2 = s.d0(an_u_pars[k], t1);
                if p1.square_distance(&p2) <= tol {
                    k += 1;
                    continue;
                }
                let lin = match Line3::new(&p1, &p2) {
                    Some(l) => l,
                    None => {
                        k += 1;
                        continue;
                    }
                };
                let mut ok = true;
                for l in (j + 1)..k {
                    if an_u_flg[l] {
                        ok = false;
                        break;
                    }
                    let pp = s.d0(an_u_pars[l], t1);
                    if lin.square_distance(&pp) <= a_defl2 {
                        continue;
                    }
                    ok = false;
                    break;
                }
                if !ok {
                    j = k - 1;
                    an_u_flg[j] = true;
                    break;
                }
                if an_u_flg[k] {
                    j = k;
                    break;
                }
                k += 1;
            }
            if k + 1 >= nbsu {
                b_cont = false;
            }
        }
        u_j_last = j;
    }

    let mut my_nb_u = an_u_flg.iter().filter(|f| **f).count();

    if my_nb_u < MY_MIN_PNTS {
        if my_nb_u == 2 {
            let nn = nbsu / MY_MIN_PNTS;
            if nn < nbsu {
                an_u_flg[nn] = true;
            }
            if nn > 0 && nbsu > nn {
                an_u_flg[nbsu - 1 - nn] = true;
            }
        } else {
            // `cxx:1437-1449`: insert into the bigger segment.
            let mut i = 2usize;
            while i <= nbsu && !an_u_flg[i - 1] {
                i += 1;
            }
            // `while (!anUFlg(i++)) ;` leaves `i` one past the first set flag.
            i += 1;
            let j = if i < nbsu / 2 {
                ((i + (nbsu - i) / 2).min(nbsu - 1)).max(1)
            } else {
                (i / 2).max(2)
            };
            if j <= nbsu {
                an_u_flg[j - 1] = true;
            }
        }
        // `cxx:1451-1452`: the ==2 branch leaves `j` from the analysis loop.
        if my_nb_u == 2 && u_j_last < nbsu {
            an_u_flg[u_j_last] = true;
        }
        my_nb_u = MY_MIN_PNTS;
    }

    if nbsv >= 2 {
        a_v_flg[0] = true;
        a_v_flg[nbsv - 1] = true;
    }
    let mut v_j_last = 0usize;
    for i in 0..nbsu {
        let t1 = an_u_pars[i];
        let mut j = 0usize;
        let mut b_cont = true;
        while j + 1 < nbsv && b_cont {
            if a_v_flg[j + 1] {
                j += 1;
                continue;
            }
            let p1 = s.d0(t1, an_v_pars[j]);
            let mut k = j + 2;
            while k < nbsv {
                let p2 = s.d0(t1, an_v_pars[k]);
                if p1.square_distance(&p2) <= tol {
                    k += 1;
                    continue;
                }
                let lin = match Line3::new(&p1, &p2) {
                    Some(l) => l,
                    None => {
                        k += 1;
                        continue;
                    }
                };
                let mut ok = true;
                for l in (j + 1)..k {
                    if a_v_flg[l] {
                        ok = false;
                        break;
                    }
                    let pp = s.d0(t1, an_v_pars[l]);
                    if lin.square_distance(&pp) <= a_defl2 {
                        continue;
                    }
                    ok = false;
                    break;
                }
                if !ok {
                    j = k - 1;
                    a_v_flg[j] = true;
                    break;
                }
                if a_v_flg[k] {
                    j = k;
                    break;
                }
                k += 1;
            }
            if k + 1 >= nbsv {
                b_cont = false;
            }
        }
        v_j_last = j;
    }

    let mut my_nb_v = a_v_flg.iter().filter(|f| **f).count();

    if my_nb_v < MY_MIN_PNTS {
        if my_nb_v == 2 {
            let nn = nbsv / MY_MIN_PNTS;
            if nn < nbsv {
                a_v_flg[nn] = true;
            }
            if nn > 0 && nbsv > nn {
                a_v_flg[nbsv - 1 - nn] = true;
            }
        } else {
            let mut i = 2usize;
            while i <= nbsv && !a_v_flg[i - 1] {
                i += 1;
            }
            i += 1;
            let j = if i < nbsv / 2 {
                ((i + (nbsv - i) / 2).min(nbsv - 1)).max(1)
            } else {
                (i / 2).max(2)
            };
            if j <= nbsv {
                a_v_flg[j - 1] = true;
            }
        }
        if my_nb_v == 2 && v_j_last < nbsv {
            a_v_flg[v_j_last] = true;
        }
        my_nb_v = MY_MIN_PNTS;
    }

    // `cxx:1573-1622`: the `theNUmin` / `theNVmin` override.
    let u_flag = (my_nb_u as i32) < the_nu_min;
    let u_pars = if u_flag {
        an_u_pars.clone()
    } else {
        an_u_flg
            .iter()
            .enumerate()
            .filter(|(_, f)| **f)
            .map(|(i, _)| an_u_pars[i])
            .collect()
    };

    let v_flag = (my_nb_v as i32) < the_nv_min;
    let v_pars = if v_flag {
        an_v_pars.clone()
    } else {
        a_v_flg
            .iter()
            .enumerate()
            .filter(|(_, f)| **f)
            .map(|(i, _)| an_v_pars[i])
            .collect()
    };

    SampleParameters { u_pars, v_pars }
}

fn uniform_fill(pars: &mut [f64], t1: f64, t2: f64) {
    let n = pars.len();
    if n == 0 {
        return;
    }
    if n == 1 {
        pars[0] = t1;
        return;
    }
    let dt = (t2 - t1) / (n as f64 - 1.0);
    pars[0] = t1;
    pars[n - 1] = t2;
    let mut t = t1 + dt;
    for p in pars.iter_mut().take(n - 1).skip(1) {
        *p = t;
        t += dt;
    }
}

/// Uniform sampling over `ComputeSamplePoints`, used when a BSpline surface
/// cannot be viewed as `Geom_BSplineSurface`.
fn sample_pnts_uniform_counts(s: &dyn Surface) -> SampleParameters {
    let (uinf, usup) = finite_uv(s.u_range().0, s.u_range().1);
    let (vinf, vsup) = finite_uv(s.v_range().0, s.v_range().1);
    let (nbsu, nbsv) = compute_sample_points(s);
    let mut u_pars = vec![0.0; nbsu.max(1) as usize];
    let mut v_pars = vec![0.0; nbsv.max(1) as usize];
    uniform_fill(&mut u_pars, uinf, usup);
    uniform_fill(&mut v_pars, vinf, vsup);
    SampleParameters { u_pars, v_pars }
}

/// `gp_Lin` built from two points, with `square_distance` of `gp_Lin`.
struct Line3 {
    p: occt_core::gp::GpPnt,
    d: [f64; 3],
}

impl Line3 {
    fn new(p1: &occt_core::gp::GpPnt, p2: &occt_core::gp::GpPnt) -> Option<Self> {
        let v = occt_core::gp::GpVec::new(p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z());
        let m = (v.x() * v.x() + v.y() * v.y() + v.z() * v.z()).sqrt();
        if m < 1.0e-300 {
            return None;
        }
        Some(Self {
            p: *p1,
            d: [v.x() / m, v.y() / m, v.z() / m],
        })
    }

    fn square_distance(&self, p: &occt_core::gp::GpPnt) -> f64 {
        let w = [p.x() - self.p.x(), p.y() - self.p.y(), p.z() - self.p.z()];
        let t = w[0] * self.d[0] + w[1] * self.d[1] + w[2] * self.d[2];
        let q = [w[0] - t * self.d[0], w[1] - t * self.d[1], w[2] - t * self.d[2]];
        q[0] * q[0] + q[1] * q[1] + q[2] * q[2]
    }
}
