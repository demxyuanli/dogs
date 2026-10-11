//! `AdvApp2Var_ApproxF2var::mmapptt_` (`AdvApp2Var_ApproxF2var.cxx:8313-8641`).
//!
//! Jacobi orders -1, 0, 1, and 2 read commons `MMAPGSS`, `MMAPGS0`,
//! `MMAPGS1`, and `MMAPGS2`.

include!("gauss_gss_xj.rs");
include!("gauss_gss_0.rs");
include!("gauss_gs0_xj.rs");
include!("gauss_gs0_0.rs");
include!("gauss_gs1_xj.rs");
include!("gauss_gs1_0.rs");
include!("gauss_gs2_xj.rs");
include!("gauss_gs2_0.rs");

fn allowed_points(nbpnts: i32) -> bool {
    matches!(nbpnts, 8 | 10 | 15 | 20 | 25 | 30 | 40 | 50 | 61)
}

fn copy_block(src: &[f64], src_at: usize, dst: &mut [f64], dst_at: usize, count: usize) -> bool {
    if src_at + count > src.len() || dst_at + count > dst.len() {
        return false;
    }
    dst[dst_at..dst_at + count].copy_from_slice(&src[src_at..src_at + count]);
    true
}

/// Gauss weights for one direction. `columns` is `JDeg - 2 * order - 1`.
pub(super) fn gauss_weights(
    ndgjac: i32,
    nbpnts: i32,
    jordre: i32,
    columns: i32,
) -> Result<Vec<f64>, String> {
    let infdg = (jordre + 1) << 1;
    if !allowed_points(nbpnts) || jordre < -1 || jordre > 2 || ndgjac >= nbpnts || ndgjac < infdg {
        return Err("AdvApp2Var_Context : Error in FORTRAN".into());
    }
    let (table_x, table_0) = if jordre == -1 {
        // `L1000` (`AdvApp2Var_ApproxF2var.cxx:8476-8505`). `infdg` is 0, so
        // the loop bound `ndgjac` matches `mxjac`.
        (GSLXJS.as_slice(), GSL0JS.as_slice())
    } else if jordre == 0 {
        // `L2000` (`AdvApp2Var_ApproxF2var.cxx:8510-8539`).
        (GSLXJ0.as_slice(), GSL0J0.as_slice())
    } else if jordre == 1 {
        (GSLXJ1.as_slice(), GSL0J1.as_slice())
    } else if jordre == 2 {
        // `L4000` (`AdvApp2Var_ApproxF2var.cxx:8580-8609`).
        (GSLXJ2.as_slice(), GSL0J2.as_slice())
    } else {
        return Err("AdvApp2Var_Context : Error in FORTRAN".into());
    };
    let half = nbpnts / 2;
    let dim = half + 1;
    let mxjac = ndgjac - infdg;
    let mut iptdb = 0;
    if nbpnts > 8 {
        iptdb += (8 - infdg) << 2;
    }
    if nbpnts > 10 {
        iptdb += (10 - infdg) * 5;
    }
    if nbpnts > 15 {
        iptdb += (15 - infdg) * 7;
    }
    if nbpnts > 20 {
        iptdb += (20 - infdg) * 10;
    }
    if nbpnts > 25 {
        iptdb += (25 - infdg) * 12;
    }
    if nbpnts > 30 {
        iptdb += (30 - infdg) * 15;
    }
    if nbpnts > 40 {
        iptdb += (40 - infdg) * 20;
    }
    if nbpnts > 50 {
        iptdb += (50 - infdg) * 25;
    }
    let mut ipdb0 = 1;
    if nbpnts > 15 {
        ipdb0 += (14 - infdg) / 2 + 1;
    }
    if nbpnts > 25 {
        ipdb0 += (24 - infdg) / 2 + 1;
    }
    let len = (dim * columns).max(0) as usize;
    let mut cgauss = vec![0.0; len];
    let count = half as usize;
    for kjac in 0..=mxjac {
        let iptt = iptdb + kjac * half + 1;
        let dst = (kjac * dim + 1) as usize;
        if !copy_block(table_x, (iptt - 1) as usize, &mut cgauss, dst, count) {
            return Err("AdvApp2Var_Context : Error in FORTRAN".into());
        }
    }
    if nbpnts % 2 == 1 {
        let mut iptt = ipdb0;
        let mut kjac = 0;
        while kjac <= mxjac {
            let dst = (kjac * dim) as usize;
            let src = (iptt - 1) as usize;
            if dst >= cgauss.len() || src >= table_0.len() {
                return Err("AdvApp2Var_Context : Error in FORTRAN".into());
            }
            cgauss[dst] = table_0[src];
            iptt += 1;
            kjac += 2;
        }
    }
    Ok(cgauss)
}
