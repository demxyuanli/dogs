//! `BSplCLib::BuildKnots`.
//! Source: `BSplCLib.cxx:1555-1754`. Knot slices are 0-based; `index` is 1-based.

/// Local knot window of length `2 * degree` (`BuildKnots`).
pub fn build_knots(degree: i32, index: i32, periodic: bool, knots: &[f64], mults: Option<&[i32]>) -> Vec<f64> {
    let deg = degree.max(0) as usize;
    let mut lk = vec![0.0; 2 * deg.max(1)];
    if let Some(mults) = mults {
        build_knots_mults(degree, index, periodic, knots, mults, &mut lk);
    } else {
        // Default arm (`cxx:1672-1682`); cases 1-6 are the same window.
        let mut j = index - degree;
        let deg2 = degree << 1;
        for i in 0..deg2 {
            j += 1;
            lk[i as usize] = knots[(j - 1) as usize];
        }
    }
    lk
}

fn knot_1(knots: &[f64], i: i32) -> f64 {
    let n = knots.len() as i32;
    if n == 0 {
        return 0.0;
    }
    let mut k = i;
    if k < 1 {
        k = 1;
    }
    if k > n {
        k = n;
    }
    knots[(k - 1) as usize]
}

fn mult_1(mults: &[i32], i: i32) -> i32 {
    let n = mults.len() as i32;
    if n == 0 {
        return 0;
    }
    let mut k = i;
    if k < 1 {
        k = 1;
    }
    if k > n {
        k = n;
    }
    mults[(k - 1) as usize]
}

/// Mults arm (`cxx:1685-1752`).
fn build_knots_mults(
    degree: i32,
    index: i32,
    periodic: bool,
    knots: &[f64],
    mults: &[i32],
    lk: &mut [f64],
) {
    let deg1 = degree - 1;
    let k_lower = 1i32;
    let k_upper = knots.len() as i32;
    let m_lower = 1i32;
    let m_upper = mults.len() as i32;
    let mut dknot = 0.0;
    let mut ilow = index;
    let mut mlow = 0;
    let mut iupp = index + 1;
    let mut mupp = 0;
    let mut loffset = 0.0;
    let mut uoffset = 0.0;
    let mut getlow = true;
    let mut getupp = true;
    if periodic {
        dknot = knot_1(knots, k_upper) - knot_1(knots, k_lower);
        if iupp > m_upper {
            iupp = m_lower + 1;
            uoffset = dknot;
        }
    }
    for i in 0..degree {
        if getlow {
            mlow += 1;
            if mlow > mult_1(mults, ilow) {
                mlow = 1;
                ilow -= 1;
                getlow = ilow >= m_lower;
                if periodic && !getlow {
                    ilow = m_upper - 1;
                    loffset = dknot;
                    getlow = true;
                }
            }
            if getlow {
                lk[(deg1 - i) as usize] = knot_1(knots, ilow) - loffset;
            }
        }
        if getupp {
            mupp += 1;
            if mupp > mult_1(mults, iupp) {
                mupp = 1;
                iupp += 1;
                getupp = iupp <= m_upper;
                if periodic && !getupp {
                    iupp = m_lower + 1;
                    uoffset = dknot;
                    getupp = true;
                }
            }
            if getupp {
                lk[(degree + i) as usize] = knot_1(knots, iupp) + uoffset;
            }
        }
    }
}
