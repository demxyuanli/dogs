//! `ComputePatches` cut decision (`AdvApp2Var_ApproxAFunc2Var.cxx:383-429`).
//!
//! `InitGrid(1)` stores one U span and one V span
//! (`AdvApp2Var_ApproxAFunc2Var.cxx:249-253`), so `NbPatchInU` and
//! `NbPatchInV` are 1 (`AdvApp2Var_Network.cxx:167-176`).

use crate::adv_approx::{Cutting, PrefAndRec};
use crate::Surface;

/// `UChoice.Value` / `VChoice.Value` then the `NumDec` assignments
/// (`AdvApp2Var_ApproxAFunc2Var.cxx:388-429`).
pub(super) fn authorize_cuts(
    surf: &dyn Surface,
    u0: f64,
    u1: f64,
    v0: f64,
    v1: f64,
    nb_u: i32,
    nb_v: i32,
    max_patches: i32,
) -> Result<(i32, f64, f64), String> {
    // `GeomConvert_ApproxSurface.cxx:381-398`: recommended cuts are the C2
    // breaks, preferred cuts are the C3 breaks, `Weight` stays 5.
    // `GeomAbs_C2` is 4 and `GeomAbs_C3` is 5.
    let u_choice = PrefAndRec::new(surf.u_intervals(4), surf.u_intervals(5), 5.0)
        .map_err(|_| "PrefAndRec : Weight is too small".to_string())?;
    let v_choice = PrefAndRec::new(surf.v_intervals(4), surf.v_intervals(5), 5.0)
        .map_err(|_| "PrefAndRec : Weight is too small".to_string())?;
    let mut u_dec = 0.0;
    let mut v_dec = 0.0;
    let u_more = u_choice.value(u0, u1, &mut u_dec);
    let v_more = v_choice.value(v0, v1, &mut v_dec);
    Ok((num_dec(nb_u, nb_v, max_patches, u_more, v_more), u_dec, v_dec))
}

fn num_dec(nb_u: i32, nb_v: i32, max_patches: i32, u_more: bool, v_more: bool) -> i32 {
    let nb_patch = nb_u * nb_v;
    let mut num_dec = 0;
    if (nb_patch + nb_v) <= max_patches && (nb_patch + nb_u) > max_patches && u_more {
        num_dec = 1;
    }
    if (nb_patch + nb_v) > max_patches && (nb_patch + nb_u) <= max_patches && v_more {
        num_dec = 2;
    }
    if (nb_patch + nb_v) <= max_patches && (nb_patch + nb_u) <= max_patches {
        if u_more {
            num_dec = 3;
        }
        if nb_v > nb_u && v_more {
            num_dec = 4;
        }
    }
    if (nb_u + 1) * (nb_v + 1) <= max_patches {
        if !u_more && !v_more {
            num_dec = 0;
        }
        if u_more && !v_more {
            num_dec = 3;
        }
        if !u_more && v_more {
            num_dec = 4;
        }
        if u_more && v_more {
            num_dec = 5;
        }
    }
    num_dec
}
