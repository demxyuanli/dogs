//! The four `ChFi3d_Builder_2.cxx` file statics that edit the stripe's
//! `SurfData` sequence.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_2.cxx:2047-2204`
//! (OCCT 8.0.0).
//!
//! `PerformSetOfSurfOnElSpine` walks a band along the guideline and, at every
//! step, splices the newly computed `SurfData` into
//! `Stripe->ChangeSetOfSurfData()->ChangeSequence()`:
//!
//! | OCCT | here |
//! |---|---|
//! | `ChFi3d_Purge` (`cxx:2047-2101`) | [`ch_fi3d_purge`] |
//! | `InsertAfter` (`cxx:2108-2132`) | [`insert_after`] |
//! | `RemoveSD` (`cxx:2136-2172`) | [`remove_sd`] |
//! | `InsertBefore` (`cxx:2179-2203`) | [`insert_before`] |
//!
//! Divergences that are not heuristics:
//! - `NCollection_Sequence<handle(ChFiDS_SurfData)>` is the `Vec<Arc<SurfData>>`
//!   of [`crate::chfi_ds::Stripe::change_set_of_surf_data`], so `Prepend` /
//!   `Append` / `InsertBefore` / `InsertAfter` are `Vec::insert` / `push` and
//!   the 1-based positions of `NCollection_Sequence` become `i - 1`,
//! - `handle == handle` is `Arc::ptr_eq`, and the null handles the statics test
//!   with `IsNull()` are `None`,
//! - `Seq.Remove(iprev, inext)` removes the **inclusive** 1-based range, which
//!   is `Vec::drain(iprev - 1 .. inext)`,
//! - `ChangeVertex(...)` / `ChangeInterference(...)` on the `SurfData` handle go
//!   through [`std::sync::Arc::make_mut`] on the caller's slot, the same
//!   convention [`super::compute_data_face_face`] uses for shared data.

use std::sync::Arc;

use crate::chfi_ds::{CommonPoint, StripeHandle, SurfData};

/// `ChFi3d_Purge(Stripe, SD, VRef, isfirst, ons, intf, intl)`
/// (`_2.cxx:2047-2101`).
///
/// Drops the part of the band that fell in the extension zone: either the whole
/// `SurfData` (`cxx:2065-2079`) or its contact on `ons`, which is collapsed onto
/// `v_ref` (`cxx:2080-2100`).
///
/// `sd` is the handle OCCT passes by reference; the port looks the same entry
/// up in the stripe's sequence (by [`Arc::ptr_eq`]) so that the `Change*` calls
/// land on the object the stripe holds instead of a detached clone.
///
/// Answers the handle as it now lives in the sequence, or `None` when the
/// section was removed from it. OCCT keeps the caller's handle alive in both
/// cases; the lookup answers `Err` when the entry is gone without being removed,
/// which OCCT cannot hit (it owns the handle).
pub fn ch_fi3d_purge(
    stripe: &StripeHandle,
    sd: &Arc<SurfData>,
    v_ref: &CommonPoint,
    is_first: bool,
    ons: i32,
    intf: &mut i32,
    intl: &mut i32,
) -> Result<Option<Arc<SurfData>>, String> {
    // `if (isfirst) intf = 1; else intl = 1;` (`cxx:2055-2063`).
    if is_first {
        *intf = 1;
    } else {
        *intl = 1;
    }
    // `int opp = 3 - ons;` (`cxx:2064`).
    let opp = 3 - ons;
    // `if (!SD->Vertex(isfirst, opp).IsOnArc() || SD->TwistOnS1() ||
    // SD->TwistOnS2())` (`cxx:2065-2066`).
    let drop_section = !sd.vertex(is_first, opp).is_on_arc() || sd.twist_on_s1() || sd.twist_on_s2();
    let mut guard = stripe
        .write()
        .map_err(|_| "ChFi3d_Purge: stripe lock poisoned".to_string())?;
    let seq = guard.change_set_of_surf_data();
    if drop_section {
        // `if (isfirst) Seq.Remove(1); else Seq.Remove(Seq.Length());`
        // (`cxx:2067-2077`).
        if is_first {
            if !seq.is_empty() {
                seq.remove(0);
            }
        } else {
            seq.pop();
        }
        return Ok(None);
    }

    let Some(pos) = seq.iter().position(|value| Arc::ptr_eq(value, sd)) else {
        return Err("ChFi3d_Purge: the SurfData is not in the stripe".into());
    };
    let data = Arc::make_mut(&mut seq[pos]);

    // `if (ons == 1) SD->ChangeIndexOfS1(0); else SD->ChangeIndexOfS2(0);`
    // (`cxx:2080-2087`).
    if ons == 1 {
        data.set_index_of_s1(0);
    } else {
        data.set_index_of_s2(0);
    }

    // `SD->ChangeVertex(!isfirst, ons) = VRef; SD->ChangeVertex(isfirst, ons) =
    // VRef;` (`cxx:2089-2090`).
    *data.vertex_mut(!is_first, ons) = v_ref.clone();
    *data.vertex_mut(is_first, ons) = v_ref.clone();

    // `ChFiDS_FaceInterference& fi = SD->ChangeInterference(ons);`
    // (`cxx:2092`).
    let fi = data.interference_mut(ons);
    // `if (isfirst) fi.SetFirstParameter(fi.LastParameter()); else
    // fi.SetLastParameter(fi.FirstParameter());` (`cxx:2093-2099`).
    if is_first {
        let last = fi.last_parameter();
        fi.set_first_parameter(last);
    } else {
        let first = fi.first_parameter();
        fi.set_last_parameter(first);
    }
    // `fi.SetLineIndex(0);` (`cxx:2100`).
    fi.set_line_index(0);
    // The `Arc::make_mut` above may have replaced the sequence entry with a
    // fresh handle (the caller still holds the old one); hand the live one back.
    let current = seq[pos].clone();
    Ok(Some(current))
}

/// `InsertAfter(Stripe, Ref, Item)` (`_2.cxx:2108-2132`).
///
/// `Standard_Failure("InsertAfter : twice the same surfdata.")`
/// (`cxx:2111-2114`) is the `Err` of the port.
pub fn insert_after(
    stripe: &StripeHandle,
    reference: Option<&Arc<SurfData>>,
    item: &Arc<SurfData>,
) -> Result<(), String> {
    // `if (Ref == Item) throw ...` (`cxx:2111-2114`).
    if reference.is_some_and(|reference| Arc::ptr_eq(reference, item)) {
        return Err("InsertAfter : twice the same surfdata.".into());
    }
    let Ok(mut guard) = stripe.write() else {
        return Err("InsertAfter : stripe lock poisoned".into());
    };
    let seq = guard.change_set_of_surf_data();
    // `if (Seq.IsEmpty() || Ref.IsNull()) Seq.Prepend(Item);`
    // (`cxx:2120-2123`).
    if seq.is_empty() || reference.is_none() {
        seq.insert(0, item.clone());
        return Ok(());
    }
    let Some(reference) = reference else {
        return Ok(());
    };
    // `for (i = 1; i <= Seq.Length(); i++) if (Seq.Value(i) == Ref)
    // { Seq.InsertAfter(i, Item); break; }` (`cxx:2124-2131`).
    for i in 0..seq.len() {
        if Arc::ptr_eq(&seq[i], reference) {
            seq.insert(i + 1, item.clone());
            break;
        }
    }
    Ok(())
}

/// `RemoveSD(Stripe, Prev, Next)` (`_2.cxx:2136-2172`): drops the `SurfData`
/// strictly between `Prev` and `Next`.
pub fn remove_sd(
    stripe: &StripeHandle,
    prev: Option<&Arc<SurfData>>,
    next: Option<&Arc<SurfData>>,
) {
    let Ok(mut guard) = stripe.write() else {
        return;
    };
    let seq = guard.change_set_of_surf_data();
    // `if (Seq.IsEmpty()) return;` (`cxx:2142-2145`).
    if seq.is_empty() {
        return;
    }
    // `int iprev = 0, inext = 0;` (`cxx:2146`).
    let mut iprev = 0usize;
    let mut inext = 0usize;
    // `for (int i = 1; i <= Seq.Length(); i++) { ... }` (`cxx:2147-2158`).
    for (k, value) in seq.iter().enumerate() {
        if prev.is_some_and(|prev| Arc::ptr_eq(value, prev)) {
            iprev = k + 2;
        }
        if next.is_some_and(|next| Arc::ptr_eq(value, next)) {
            inext = k;
            break;
        }
    }
    // `if (Prev.IsNull()) iprev = 1; if (Next.IsNull()) inext = Seq.Length();`
    // (`cxx:2159-2165`).
    if prev.is_none() {
        iprev = 1;
    }
    if next.is_none() {
        inext = seq.len();
    }
    // `if (iprev <= inext) Seq.Remove(iprev, inext);` (`cxx:2166-2171`).
    // `iprev` stays `0` when `Prev` is neither null nor in the sequence; OCCT
    // would then raise `Standard_OutOfRange` from `Remove(0, ...)`, which the
    // port skips (`iprev >= 1`).
    if iprev <= inext && iprev >= 1 {
        seq.drain(iprev - 1..inext);
    }
}

/// `InsertBefore(Stripe, Ref, Item)` (`_2.cxx:2179-2203`).
///
/// `Standard_Failure("InsertBefore : twice the same surfdata.")`
/// (`cxx:2182-2185`) is the `Err` of the port.
pub fn insert_before(
    stripe: &StripeHandle,
    reference: Option<&Arc<SurfData>>,
    item: &Arc<SurfData>,
) -> Result<(), String> {
    // `if (Ref == Item) throw ...` (`cxx:2182-2185`).
    if reference.is_some_and(|reference| Arc::ptr_eq(reference, item)) {
        return Err("InsertBefore : twice the same surfdata.".into());
    }
    let Ok(mut guard) = stripe.write() else {
        return Err("InsertBefore : stripe lock poisoned".into());
    };
    let seq = guard.change_set_of_surf_data();
    // `if (Seq.IsEmpty() || Ref.IsNull()) Seq.Append(Item);`
    // (`cxx:2191-2194`).
    if seq.is_empty() || reference.is_none() {
        seq.push(item.clone());
        return Ok(());
    }
    let Some(reference) = reference else {
        return Ok(());
    };
    // `for (i = 1; i <= Seq.Length(); i++) if (Seq.Value(i) == Ref)
    // { Seq.InsertBefore(i, Item); break; }` (`cxx:2195-2202`).
    for i in 0..seq.len() {
        if Arc::ptr_eq(&seq[i], reference) {
            seq.insert(i, item.clone());
            break;
        }
    }
    Ok(())
}
