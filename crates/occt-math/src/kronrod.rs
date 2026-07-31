//! Kronrod adaptive integration. Source: `math_KronrodSingleIntegration.cxx`
use crate::MathStatus;

/// Kronrod 15-point rule for adaptive integration on [a, b].
const GK15_NODES: [f64; 15] = [-0.9914553711208126,-0.9491079123427585,-0.8648644233597691,-0.7415311855993945,-0.5860872354676911,-0.4058451513773972,-0.20778495500789848,0.0,0.20778495500789848,0.4058451513773972,0.5860872354676911,0.7415311855993945,0.8648644233597691,0.9491079123427585,0.9914553711208126];
const GK15_WG: [f64; 4] = [0.4179591836734694,0.3818300505051189,0.27970539148927664,0.1294849661688697]; // Gauss 7-point weights
const GK15_WK: [f64; 15] = [0.022935322010529224,0.06309209262997856,0.10479001032225019,0.14065325971552592,0.1690047266392679,0.19035057806478542,0.2044329400752989,0.20948214108472782,0.2044329400752989,0.19035057806478542,0.1690047266392679,0.14065325971552592,0.10479001032225019,0.06309209262997856,0.022935322010529224];

/// Integrate f over [a, b] adaptively to tolerance. Returns (integral, error_estimate).
pub fn integrate_adaptive<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64, tol: f64, max_depth: usize) -> Result<(f64, f64), MathStatus> {
    integrate_recursive(f, a, b, tol, max_depth, 0)
}

fn integrate_recursive<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64, tol: f64, max_depth: usize, depth: usize) -> Result<(f64, f64), MathStatus> {
    if depth > max_depth { return Err(MathStatus::TooManyIterations); }
    let xm = 0.5*(b+a); let xl = 0.5*(b-a);
    let mut g7 = 0.0f64; let mut k15 = 0.0f64;
    for i in 0..15 {
        let x = xm + xl * GK15_NODES[i];
        let fx = f(x);
        k15 += fx * GK15_WK[i];
        if i % 2 == 0 { /* 7-point subset */ }
    }
    // 7-point subset from 15-point nodes
    g7 = (0..7).map(|i| f(xm + xl * GK15_NODES[2*i+1]) * GK15_WG[i/2]).sum();
    g7 *= xl; k15 *= xl;

    let error = (200.0 * (g7 - k15).abs()).cbrt();
    if error < tol { return Ok((k15, error)); }

    let mid = 0.5*(a+b);
    let (left, el) = integrate_recursive(f, a, mid, tol*0.5, max_depth, depth+1)?;
    let (right, er) = integrate_recursive(f, mid, b, tol*0.5, max_depth, depth+1)?;
    Ok((left+right, el+er))
}

##[cfg(test)]
#mod tests {
#    use super::*;
#    #[test]
#    fn integrate_sin() {
#        let (v, _) = integrate_adaptive(&|x| x.sin(), 0.0, std::f64::consts::PI, 1e-10, 20).unwrap();
#        assert!((v-2.0).abs() < 1e-10);
#    }
#}
