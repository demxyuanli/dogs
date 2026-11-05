//! `ChFi3d_SearchFD` — find the pair of surf data that faces each other.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:4549-4663`
//! (OCCT 8.0.0).
//!
//! The walk starts at `(ind1, ind2)` and grows the two windows one step at a
//! time along `sens1` / `sens2` until `ChFi3d_IsInFront` accepts a pair or both
//! windows are exhausted (`fini1` / `fini2`). The first pair that only answered
//! `visavis` is remembered (`visavisok`) and answers when nothing better shows
//! up.
//!
//! Translation notes:
//! - `i1` / `i2` / `p1` / `p2` / `face` / `sameside` / `jf1` / `jf2` are
//!   out-parameters in OCCT and are the fields of [`SearchFdAnswer`]; the
//!   `bool` return is [`SearchFdAnswer::found`],
//! - `TopoDS_Vertex Vtx;` (`cxx:4573`) is never set, and both calls pass
//!   `Check2dDistance = false` (`cxx:4597-4598`, `:4638-4639`), so the vertex is
//!   never read: the port hands `None`,
//! - the `for` loops advance by `sens` and OCCT relies on `sens` being `+1` /
//!   `-1`; that is what `ChFi3d_IndexOfSurfData` answers, and the port keeps the
//!   same arithmetic instead of re-deriving a range.

use crate::chfi3d::builder_0_helpers::is_in_front::{ch_fi3d_is_in_front, IsInFrontState};
use crate::chfi_ds::Stripe;
use crate::shape::Face;
use crate::topopebrep::ds::HDataStructure;

/// The out-parameters of `ChFi3d_SearchFD` plus its `bool` return.
#[derive(Debug, Clone)]
pub struct SearchFdAnswer {
    /// `found`
    pub found: bool,
    /// `i1`
    pub i1: i32,
    /// `i2`
    pub i2: i32,
    /// `p1`
    pub p1: f64,
    /// `p2`
    pub p2: f64,
    /// `face`
    pub face: Option<Face>,
    /// `sameside`
    pub sameside: bool,
    /// `jf1`
    pub jf1: i32,
    /// `jf2`
    pub jf2: i32,
}

impl Default for SearchFdAnswer {
    fn default() -> Self {
        Self {
            found: false,
            i1: 0,
            i2: 0,
            p1: 0.0,
            p2: 0.0,
            face: None,
            sameside: false,
            jf1: 0,
            jf2: 0,
        }
    }
}

/// `ChFi3d_SearchFD(DStr, cd1, cd2, sens1, sens2, i1, i2, p1, p2, ind1, ind2,
/// face, sameside, jf1, jf2)` (`cxx:4549-4663`).
#[allow(clippy::too_many_arguments)]
pub fn ch_fi3d_search_fd(
    ds: &HDataStructure,
    cd1: &Stripe,
    cd2: &Stripe,
    sens1: i32,
    sens2: i32,
    ind1: i32,
    ind2: i32,
) -> Result<SearchFdAnswer, String> {
    // `cxx:4566-4573`.
    let mut ans = SearchFdAnswer::default();
    // `int id1 = ind1, id2 = ind2;` — the `for` initials, read every round and
    // never written (`cxx:4567`).
    let id1 = ind1;
    let id2 = ind2;
    let mut if1 = ind1;
    let mut if2 = ind2;
    let l1 = cd1.set_of_surf_data().len() as i32;
    let l2 = cd2.set_of_surf_data().len() as i32;
    let mut fini1 = false;
    let mut fini2 = false;
    let mut visavisok = false;

    // `while (!found) { ... }` (`cxx:4575-4661`). The state carries `p1` / `p2`
    // across the calls, exactly like the OCCT out-parameters.
    let mut st = IsInFrontState::new(ans.p1, ans.p2);
    loop {
        // `for (i = id1; (i*sens1) <= (if1*sens1) && !found && !fini2; i += sens1)`
        // (`cxx:4577`).
        let mut i = id1;
        while i * sens1 <= if1 * sens1 && !ans.found && !fini2 {
            st.ok = false;
            ch_fi3d_is_in_front(ds, cd1, cd2, i, if2, sens1, sens2, &mut st, None, false, false)?;
            absorb(&mut ans, &st);
            if st.ok {
                ans.i1 = i;
                ans.i2 = if2;
                ans.found = true;
            } else if st.visavis && !visavisok {
                visavisok = true;
                ans.i1 = i;
                ans.i2 = if2;
            }
            i += sens1;
        }

        // `if (!fini1) { if1 += sens1; if (if1 < 1 || if1 > l1) { if1 -= sens1;
        //  fini1 = true; } }` (`cxx:4605-4613`).
        if !fini1 {
            if1 += sens1;
            if if1 < 1 || if1 > l1 {
                if1 -= sens1;
                fini1 = true;
            }
        }

        // `for (i = id2; (i*sens2) <= (if2*sens2) && !found && !fini1; i += sens2)`
        // (`cxx:4616`).
        let mut i = id2;
        while i * sens2 <= if2 * sens2 && !ans.found && !fini1 {
            st.ok = false;
            ch_fi3d_is_in_front(ds, cd1, cd2, if1, i, sens1, sens2, &mut st, None, false, false)?;
            absorb(&mut ans, &st);
            if st.ok {
                ans.i1 = if1;
                ans.i2 = i;
                ans.found = true;
            } else if st.visavis && !visavisok {
                visavisok = true;
                ans.i1 = if1;
                ans.i2 = i;
            }
            i += sens2;
        }

        // `if (!fini2) { if2 += sens2; if (if2 < 1 || if2 > l2) { if2 -= sens2;
        //  fini2 = true; } }` (`cxx:4648-4656`).
        if !fini2 {
            if2 += sens2;
            if if2 < 1 || if2 > l2 {
                if2 -= sens2;
                fini2 = true;
            }
        }

        // `if (fini1 && fini2) break;` (`cxx:4657-4660`).
        if fini1 && fini2 {
            break;
        }
    }
    Ok(ans)
}

/// Copies the out-parameters that every `ChFi3d_IsInFront` call writes
/// (`cxx:4592-4598` / `:4635-4639`): `p1` / `p2` / `face` / `sameside` /
/// `jf1` / `jf2`. `i1` / `i2` are assigned by the caller only on the two
/// OCCT branches that write them.
fn absorb(ans: &mut SearchFdAnswer, st: &IsInFrontState) {
    ans.p1 = st.p1;
    ans.p2 = st.p2;
    ans.face = st.face.clone();
    ans.sameside = st.sameside;
    ans.jf1 = st.jf1;
    ans.jf2 = st.jf2;
}
