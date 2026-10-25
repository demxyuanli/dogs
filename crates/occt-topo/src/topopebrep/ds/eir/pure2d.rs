//! The 2d reduction helpers of `TopOpeBRepDS_EIR`.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_EIR.cxx` (OCCT 8.0.0):
//! - `FUN_keepl3dF` (`:690-717`)
//! - `FUN_reducepure2dI0` (`:719-772`)
//! - `FUN_reducepure2dI` (`:774-798`)
//!
//! `FUN_keepl3dF` takes two unnamed parameters in the C++ definition
//! (`cxx:690-691`: an `int` and an `HDataStructure` handle); the port drops the
//! handle, which OCCT never reads, and keeps the index as `_eix`.

use crate::abs::Orientation;

use super::super::export::fds_idata;
use super::super::interference_iterator::InterferenceList;
use super::super::kind::Kind;
use super::super::tki::TKI;
use super::super::transition::Transition;

/// `FUN_keepl3dF(EIX, HDS, l3dF, LR3dFE, l3dFkeep)` (`cxx:690-717`).
///
/// For every `(Kind, G)` group of `l3dF` with no group at all in `LR3dFE`, the
/// face interferences are kept: no reduced 3d interference exists for that
/// geometry.
pub fn fun_keep_l3d_f(
    _eix: i32,
    l3d_f: &InterferenceList,
    lr3d_fe: &InterferenceList,
    l3d_f_keep: &mut InterferenceList,
) {
    let mut tki = TKI::new();
    tki.fill_on_geometry(l3d_f);

    let mut tki_r = TKI::new();
    tki_r.fill_on_geometry(lr3d_fe);

    let mut order: Vec<(Kind, i32)> = Vec::new();
    tki.init();
    while tki.more() {
        if let Some((k, g, _)) = tki.value() {
            order.push((k, g));
        }
        tki.next();
    }

    for (k, g) in order {
        let is_bound = tki_r.is_bound(k, g);
        if !is_bound {
            l3d_f_keep.extend(tki.interferences(k, g).iter().cloned());
        }
    }
}

/// `FUN_reducepure2dI0(LI, RLI)` (`cxx:719-772`).
///
/// `xpu210798 : CTS21216 (e7,G=P3,S=E9,FTRASHA=f8,f10)` — for a section edge,
/// `2dI1=(REVERSED(ftrasha1),G,ES)` and `2dI2=(FORWARD(ftrasha2),G,ES)` reduce
/// to `2dR=(IN(ftrasha1,ftrasha2))`, which `I1` is rewritten into; `LI` is
/// emptied and `RLI` receives `I1`.
pub fn fun_reduce_pure2d_i0(li: &mut InterferenceList, rli: &mut InterferenceList) {
    let Some(h1) = li.first().cloned() else {
        return;
    };
    let Some(h2) = li.last().cloned() else {
        return;
    };

    let (o1, tsb1, isb1, _tsa1, isa1, _gt1, g1, _st1, _s1) = {
        let i1 = h1.read().expect("interference lock");
        let o = i1
            .transition()
            .orientation(crate::abs::State::In, crate::abs::ShapeType::Face);
        let (tsb, isb, tsa, isa, gt, g, st, s) = fds_idata(&i1);
        (o, tsb, isb, tsa, isa, gt, g, st, s)
    };

    let (o2, tsb2, isb2, _tsa2, isa2, _gt2, g2, _st2, _s2) = {
        let i2 = h2.read().expect("interference lock");
        let o = i2
            .transition()
            .orientation(crate::abs::State::In, crate::abs::ShapeType::Face);
        let (tsb, isb, tsa, isa, gt, g, st, s) = fds_idata(&i2);
        (o, tsb, isb, tsa, isa, gt, g, st, s)
    };

    if isb1 == isb2 {
        return; // xpu170898 FOR+REV ->INT/EXT (cto009B1(es6,p1,f9))
    }

    // recall : G1==G2 && ST1==ST2==EDGE
    let mut ok = g1 == g2;
    ok = ok
        && tsb1 == crate::abs::ShapeType::Face
        && tsb1 == tsb2
        && isb1 == isa1
        && isb2 == isa2;
    if !ok {
        return;
    }

    let int12 = o1 == Orientation::Reversed && o2 == Orientation::Forward;
    let int21 = o2 == Orientation::Reversed && o1 == Orientation::Forward;
    if !(int12 || int21) {
        return;
    }

    let mut new_t = Transition::from_orientation(Orientation::Internal);
    let bef = if int12 { isb1 } else { isb2 };
    let aft = if int21 { isb1 } else { isb2 };
    new_t.set_index_before(bef);
    new_t.set_index_after(aft);
    h1.write().expect("interference lock").set_transition(new_t);
    rli.push(h1);
    li.clear();
}

/// `FUN_reducepure2dI(LI, RLI)` (`cxx:774-798`).
///
/// Applies [`fun_reduce_pure2d_i0`] to the `(Kind == POINT)` groups of exactly
/// two elements; `LI` keeps every group in `TKI` order afterwards
/// (`cxx:795-797`).
pub fn fun_reduce_pure2d_i(li: &mut InterferenceList, rli: &mut InterferenceList) {
    let mut tki = TKI::new();
    tki.fill_on_geometry(li);

    let mut order: Vec<(Kind, i32)> = Vec::new();
    tki.init();
    while tki.more() {
        if let Some((k, g, _)) = tki.value() {
            order.push((k, g));
        }
        tki.next();
    }

    let mut new_li: InterferenceList = Vec::new();
    for (k, g) in order {
        let mut loi: InterferenceList = tki.interferences(k, g).to_vec();
        let mut rloi: InterferenceList = Vec::new();
        let nloi = loi.len();
        let ok = nloi == 2 && k == Kind::Point;
        if ok {
            fun_reduce_pure2d_i0(&mut loi, &mut rloi);
        }
        rli.extend(rloi);
        new_li.extend(loi);
    }
    li.clear();
    li.extend(new_li);
}
