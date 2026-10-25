//! The OCCT process-wide `GLOBAL_*` state shared by
//! `TopOpeBRepBuild_ffsfs.cxx` and `TopOpeBRepBuild_GridSS.cxx`.
//!
//! OCCT keeps these in file-scope `Standard_EXPORT` globals; the port keeps
//! one thread-local record per process so the two source files can drive the
//! same state without threading a context parameter through the whole
//! `Builder` grid.
//!
//! | OCCT | declaration |
//! |---|---|
//! | `GLOBAL_lfr1` | `GridSS.cxx:827`, `ffsfs.cxx` |
//! | `GLOBAL_lfrtoprocess` | `GridSS.cxx:829` |
//! | `GLOBAL_revownsplfacori` | `GridSS.cxx:807` |
//! | `GLOBAL_SplitAnc` | `ffsfs.cxx:230-231` |
//! | `STATIC_motheropedef` / `STATIC_Gmotherope` | `ffsfs.cxx:55-56` |
//! | `stabuild_IMELF1/2`, `stabuild_IDMEALF1/2`, `stabuild_IMEF`, `static_CONF1/2` | `GridSS.cxx:718-729` |

use std::cell::RefCell;

use crate::abs::{Orientation, ShapeType, State};
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, shape_map_key, ShapeMapKey};
use crate::topopebrep::build::g_topo::GTopo;
use crate::topopebrep::ds::kind::Config;

use super::Builder;

#[derive(Default)]
struct SsGlobals {
    /// `GLOBAL_lfr1`.
    lfr1: Option<Vec<TopoShape>>,
    /// `GLOBAL_lfrtoprocess`.
    lmfrtoprocess: bool,
    /// `GLOBAL_revownsplfacori`.
    revownsplfacori: bool,
    /// `GLOBAL_SplitAnc`; `None` is OCCT's `nullptr`.
    split_anc: Option<Vec<(ShapeMapKey, TopoShape, i32)>>,
    /// `STATIC_motheropedef`.
    motheropedef: bool,
    /// `STATIC_Gmotherope`.
    motherope: GTopo,
    /// `stabuild_IMELF1`.
    imelf1: Vec<(ShapeMapKey, TopoShape)>,
    /// `stabuild_IMELF2`.
    imelf2: Vec<(ShapeMapKey, TopoShape)>,
    /// `stabuild_IDMEALF1`: split shape -> ancestor faces.
    idmealf1: Vec<(ShapeMapKey, Vec<TopoShape>)>,
    /// `stabuild_IDMEALF2`.
    idmealf2: Vec<(ShapeMapKey, Vec<TopoShape>)>,
    /// `static_CONF1`.
    conf1: Config,
    /// `static_CONF2`.
    conf2: Config,
}

thread_local! {
    static SS_GLOBALS: RefCell<SsGlobals> = RefCell::new(SsGlobals::default());
}

/// `GLOBAL_lfr1->Clear()` (`ffsfs.cxx:363`, `:450`, `GridSS.cxx`).
pub fn lfr1_clear() {
    SS_GLOBALS.with(|g| g.borrow_mut().lfr1 = Some(Vec::new()));
}

/// `GLOBAL_lfr1->Append(S)` (`GridSS.cxx:968`).
pub fn lfr1_append(s: &TopoShape) {
    SS_GLOBALS.with(|g| {
        let mut g = g.borrow_mut();
        g.lfr1.get_or_insert_with(Vec::new).push(s.clone());
    });
}

/// `GLOBAL_lfr1` as a list; `None` before the first allocation.
pub fn lfr1_snapshot() -> Vec<TopoShape> {
    SS_GLOBALS.with(|g| g.borrow().lfr1.clone().unwrap_or_default())
}

/// `GLOBAL_lfrtoprocess`.
pub fn lfrtoprocess() -> bool {
    SS_GLOBALS.with(|g| g.borrow().lmfrtoprocess)
}

/// Assign `GLOBAL_lfrtoprocess`.
pub fn set_lfrtoprocess(v: bool) {
    SS_GLOBALS.with(|g| g.borrow_mut().lmfrtoprocess = v);
}

/// `GLOBAL_revownsplfacori`.
pub fn revownsplfacori() -> bool {
    SS_GLOBALS.with(|g| g.borrow().revownsplfacori)
}

/// Assign `GLOBAL_revownsplfacori`.
pub fn set_revownsplfacori(v: bool) {
    SS_GLOBALS.with(|g| g.borrow_mut().revownsplfacori = v);
}

/// `GLOBAL_SplitAnc` allocation (`ffsfs.cxx:374-381`): allocate when null,
/// then `Clear()`.
pub fn split_anc_prepare() {
    SS_GLOBALS.with(|g| g.borrow_mut().split_anc = Some(Vec::new()));
}

/// `GLOBAL_SplitAnc->Clear()` (`ffsfs.cxx:381`, `:613`).
pub fn split_anc_clear() {
    SS_GLOBALS.with(|g| {
        if let Some(m) = g.borrow_mut().split_anc.as_mut() {
            m.clear();
        }
    });
}

/// `GLOBAL_SplitAnc != nullptr`.
pub fn split_anc_defined() -> bool {
    SS_GLOBALS.with(|g| g.borrow().split_anc.is_some())
}

/// `GLOBAL_SplitAnc->IsBound(F)`.
pub fn split_anc_is_bound(f: &TopoShape) -> bool {
    let key = shape_map_key(f);
    SS_GLOBALS.with(|g| {
        g.borrow()
            .split_anc
            .as_ref()
            .map(|m| m.iter().any(|(k, _, _)| *k == key))
            .unwrap_or(false)
    })
}

/// `GLOBAL_SplitAnc->Find(F)`; `None` when unbound.
pub fn split_anc_find(f: &TopoShape) -> Option<i32> {
    let key = shape_map_key(f);
    SS_GLOBALS.with(|g| {
        g.borrow().split_anc.as_ref().and_then(|m| {
            m.iter()
                .find(|(k, _, _)| *k == key)
                .map(|(_, _, i)| *i)
        })
    })
}

/// `GLOBAL_SplitAnc->Bind(F, I)`.
pub fn split_anc_bind(f: &TopoShape, i: i32) {
    let key = shape_map_key(f);
    SS_GLOBALS.with(|g| {
        let mut g = g.borrow_mut();
        if let Some(m) = g.split_anc.as_mut() {
            if let Some(e) = m.iter_mut().find(|(k, _, _)| *k == key) {
                e.2 = i;
            } else {
                m.push((key, f.clone(), i));
            }
        }
    });
}

/// `FUN_setmotherope(G)` (`ffsfs.cxx:57-61`).
pub fn set_motherope(g: &GTopo) {
    SS_GLOBALS.with(|s| {
        let mut s = s.borrow_mut();
        s.motherope = g.clone();
        s.motheropedef = true;
    });
}

/// `FUN_unsetmotherope()` (`ffsfs.cxx:63-66`).
pub fn unset_motherope() {
    SS_GLOBALS.with(|s| s.borrow_mut().motheropedef = false);
}

/// `FUN_ismotheropedef()` (`ffsfs.cxx:68-71`).
pub fn is_motherope_defined() -> bool {
    SS_GLOBALS.with(|s| s.borrow().motheropedef)
}

/// `FUN_motherope()` (`ffsfs.cxx:73-76`).
pub fn motherope() -> GTopo {
    SS_GLOBALS.with(|s| s.borrow().motherope.clone())
}

/// `FUNBUILD_MAPSUBSHAPES(LOFS, T, IM)` (`GridSS.cxx:655-665`) into the given
/// map, appending only shapes not already present (an
/// `NCollection_IndexedMap::Add`).
fn map_sub_shapes_into(
    lofs: &[TopoShape],
    t: ShapeType,
    im: &mut Vec<(ShapeMapKey, TopoShape)>,
) {
    for s in lofs {
        for sub in explorer_shapes(s, t) {
            let key = shape_map_key(&sub);
            if !im.iter().any(|(k, _)| *k == key) {
                im.push((key, sub));
            }
        }
    }
}

/// `FUNBUILD_MAPANCSPLSHAPES(B, S, STATE, IDM)` (`GridSS.cxx:668-691`).
fn map_anc_split_shapes_into(
    b: &Builder,
    s: &TopoShape,
    state: State,
    idm: &mut Vec<(ShapeMapKey, Vec<TopoShape>)>,
) {
    if !b.is_split(s, state) {
        return;
    }
    for sps in b.splits(s, state).to_vec() {
        let key = shape_map_key(&sps);
        if let Some(e) = idm.iter_mut().find(|(k, _)| *k == key) {
            e.1.push(s.clone());
        } else {
            idm.push((key, vec![s.clone()]));
        }
    }
}

/// `FUNBUILD_MAPANCSPLSHAPES(B, S, IDM)` (`GridSS.cxx:693-702`).
fn map_anc_split_shapes_into_in_out(
    b: &Builder,
    s: &TopoShape,
    idm: &mut Vec<(ShapeMapKey, Vec<TopoShape>)>,
) {
    map_anc_split_shapes_into(b, s, State::In, idm);
    map_anc_split_shapes_into(b, s, State::Out, idm);
}

/// `FUNBUILD_ANCESTORRANKPREPARE(B, LF1, LF2, CONF1, CONF2)`
/// (`GridSS.cxx:731-745`).
pub fn build_ancestor_rank_prepare(
    b: &Builder,
    lf1: &[TopoShape],
    lf2: &[TopoShape],
    conf1: Config,
    conf2: Config,
) {
    SS_GLOBALS.with(|g| {
        let mut g = g.borrow_mut();
        g.conf1 = conf1;
        g.conf2 = conf2;
        g.imelf1.clear();
        g.imelf2.clear();
        g.idmealf1.clear();
        g.idmealf2.clear();
        map_sub_shapes_into(lf1, ShapeType::Edge, &mut g.imelf1);
        map_sub_shapes_into(lf2, ShapeType::Edge, &mut g.imelf2);
        // `cxx:743-744`: the ancestor maps are built from the edge maps, which
        // are traversed in insertion order.
        let imelf1 = g.imelf1.clone();
        let imelf2 = g.imelf2.clone();
        for (_, e) in &imelf1 {
            map_anc_split_shapes_into_in_out(b, e, &mut g.idmealf1);
        }
        for (_, e) in &imelf2 {
            map_anc_split_shapes_into_in_out(b, e, &mut g.idmealf2);
        }
    });
}

/// `FUNBUILD_ANCESTORRANKGET(B, f, of1, of2)` (`GridSS.cxx:748-776`).
///
/// `stabuild_IMEF` is reset per call (`cxx:753` reuses the static map
/// without `Clear()`, which is equivalent to a fresh `MapShapes` on the same
/// face for the `Contains` tests below).
pub fn build_ancestor_rank_get(f: &TopoShape) -> (bool, bool) {
    let edges = explorer_shapes(f, ShapeType::Edge);
    let mut of1 = false;
    for e in &edges {
        let key = shape_map_key(e);
        of1 = SS_GLOBALS.with(|g| g.borrow().idmealf1.iter().any(|(k, _)| *k == key));
        if of1 {
            break;
        }
    }
    let mut of2 = false;
    for e in &edges {
        let key = shape_map_key(e);
        of2 = SS_GLOBALS.with(|g| g.borrow().idmealf2.iter().any(|(k, _)| *k == key));
        if of2 {
            break;
        }
    }
    (of1, of2)
}

/// `FUNBUILD_ORIENTLOFS(B, TB1, TB2, LOFS)` (`GridSS.cxx:779-803`).
pub fn build_orient_lofs(tb1: State, tb2: State, lofs: &mut [TopoShape]) {
    let r12 = Builder::reverse(tb1, tb2);
    let r21 = Builder::reverse(tb2, tb1);
    for f in lofs.iter_mut() {
        let (of1, of2) = build_ancestor_rank_get(f);
        let rf = if of1 && !of2 {
            r12
        } else if of2 && !of1 {
            r21
        } else {
            false
        };
        let newori = Builder::orient(f.orientation(), rf);
        f.set_orientation(newori);
    }
}

/// `static_CONF1` / `static_CONF2` as written by
/// `FUNBUILD_ANCESTORRANKPREPARE`; exposed for callers that mirror OCCT's
/// static reads.
pub fn ancestor_rank_confs() -> (Config, Config) {
    SS_GLOBALS.with(|g| {
        let g = g.borrow();
        (g.conf1, g.conf2)
    })
}

/// `TopAbs::Complement` on an orientation, used by `GSplitFaceSFS`
/// (`GridSS.cxx:926`, `:958`).
pub fn complement(o: Orientation) -> Orientation {
    match o {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        Orientation::Internal => Orientation::External,
        Orientation::External => Orientation::Internal,
    }
}
