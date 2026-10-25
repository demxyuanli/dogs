//! `TopOpeBRepTool::RegularizeWires` / `Regularize` / `RegularizeFace` and the
//! wire classification they share.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepTool/TopOpeBRepTool_RegularizeW.cxx`
//! (OCCT 8.0.0, 887 lines): `RegularizeWires` (`:86-133`), `Regularize`
//! (`:135-157`), `FUN_tool_ClassifW` (`:676-850`) and `RegularizeFace`
//! (`:853-887`).
//!
//! UNPORTED:
//! - `FUN_tool_mkboundedF`, `FUN_tool_FindAPointInTheFace`,
//!   `FUN_tool_GetFP2d`, `FUN_tool_classiwithp2d`, `FUN_tool_mkBnd2d` and the
//!   first `FUN_tool_ClassifW` (`cxx:298-673`) are inside a `/* ... */`
//!   comment in OCCT (the classification moved to `TopOpeBRepTool_CLASSI`).
//! - the `OCCT_DEBUG` printers (`FUN_adds`, `FUN_tool_coutsta`,
//!   `FUN_tool_MakeFaces`) only print.
//! - `Regularize` constructs an unused `TopOpeBRepTool_REGUW REGUW(theFace)`
//!   (`cxx:140`) whose constructor only stores `Fref`; it has no other effect.

use std::collections::HashMap;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::shape::Face;
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, is_equal, shape_map_key, ShapeMapKey};
use crate::topopebrep::tool::classi::{
    fun_add_owlw, Classi, DIFF, ONE_IN_TWO, SAME, TWO_IN_ONE, UNKNOWN,
};
use crate::topopebrep::tool::geometry::tol_uv;
use crate::topopebrep::tool::reguw::Reguw;
use crate::topopebrep::tool::wire_to_faces::wire_to_face;

/// `NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>`
/// (`cxx:87-89`, `:676-681`): keyed by `TopoDS_Shape::IsEqual`, insertion
/// ordered because `ClassifW` walks the map.
#[derive(Clone, Default)]
pub struct WireMap {
    pairs: Vec<(TopoShape, Vec<TopoShape>)>,
}

impl WireMap {
    /// `WireMap()`.
    pub fn new() -> Self {
        Self::default()
    }

    /// `IsBound(K)`.
    pub fn is_bound(&self, k: &TopoShape) -> bool {
        self.pairs.iter().any(|(x, _)| is_equal(x, k))
    }

    /// `Find(K)`.
    pub fn find(&self, k: &TopoShape) -> Option<&Vec<TopoShape>> {
        self.pairs
            .iter()
            .find(|(x, _)| is_equal(x, k))
            .map(|(_, v)| v)
    }

    /// `ChangeFind(K)`.
    pub fn change_find(&mut self, k: &TopoShape) -> Option<&mut Vec<TopoShape>> {
        self.pairs
            .iter_mut()
            .find(|(x, _)| is_equal(x, k))
            .map(|(_, v)| v)
    }

    /// `Bind(K, I)`.
    pub fn bind(&mut self, k: &TopoShape, item: Vec<TopoShape>) {
        if let Some((_, v)) = self.pairs.iter_mut().find(|(x, _)| is_equal(x, k)) {
            *v = item;
        } else {
            self.pairs.push((k.clone(), item));
        }
    }

    /// `UnBind(K)`.
    pub fn unbind(&mut self, k: &TopoShape) {
        self.pairs.retain(|(x, _)| !is_equal(x, k));
    }

    /// `Extent()`.
    pub fn extent(&self) -> usize {
        self.pairs.len()
    }

    /// The `k`-th key, 1-based (`FindKey(k)`).
    pub fn key_at(&self, k: usize) -> Option<&TopoShape> {
        self.pairs.get(k.checked_sub(1)?).map(|(x, _)| x)
    }

    /// The `k`-th item, 1-based (`FindFromIndex(k)`).
    pub fn item_at(&self, k: usize) -> Option<&Vec<TopoShape>> {
        self.pairs.get(k.checked_sub(1)?).map(|(_, v)| v)
    }

    /// `Iterator` over `(key, item)`.
    pub fn iter(&self) -> impl Iterator<Item = &(TopoShape, Vec<TopoShape>)> {
        self.pairs.iter()
    }
}

/// `TopOpeBRepTool::RegularizeWires(theFace, mapoldWnewW, ESplits)`
/// (`cxx:86-133`).
///
/// `map_ow_nw` binds each wire of `theFace` to its splits (an empty list when
/// the wire is unchanged) and is updated in place; `e_splits` binds each split
/// edge to its splits.
pub fn regularize_wires(
    the_face: &TopoShape,
    map_ow_nw: &mut Vec<(TopoShape, Vec<TopoShape>)>,
    e_splits: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
) -> bool {
    // `cxx:93-95`: a null face answers false; `TopoShape::is_null` is always
    // false in this port, so the guard cannot fire.
    if the_face.is_null() {
        return false;
    }
    let a_face = the_face.oriented(Orientation::Forward);

    let mut reguw = Reguw::new(&a_face);

    // `REGUW.SetOwNw` / `SetEsplits` (`cxx:99-100`)
    let ow_nw: HashMap<ShapeMapKey, Vec<TopoShape>> = map_ow_nw
        .iter()
        .map(|(w, l)| (shape_map_key(w), l.clone()))
        .collect();
    reguw.set_ow_nw(&ow_nw);
    reguw.set_e_splits(e_splits);

    // `cxx:103-122`
    for w in explorer_shapes(&a_face, ShapeType::Wire) {
        reguw.init(&w);
        if !reguw.map_s() {
            return false;
        }
        if !reguw.split_eds() {
            return false;
        }
        if !reguw.regu_all() {
            return false;
        }
    }

    // `REGUW.GetEsplits(ESplits)` / `GetOwNw(mapoldWnewW)` (`cxx:124-125`).
    *e_splits = reguw.e_splits().clone();
    let after = reguw.ow_nw().clone();
    for (w, l) in map_ow_nw.iter_mut() {
        if let Some(new_l) = after.get(&shape_map_key(w)) {
            *l = new_l.clone();
        }
    }
    true
}

/// `TopOpeBRepTool::Regularize(theFace, aListOfFaces, ESplits)` (`cxx:135-157`).
pub fn regularize(
    the_face: &TopoShape,
    a_list_of_faces: &mut Vec<TopoShape>,
    e_splits: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
) -> bool {
    a_list_of_faces.clear();
    let mut map_old_w_new_w: Vec<(TopoShape, Vec<TopoShape>)> = Vec::new();
    let mut regu = regularize_wires(the_face, &mut map_old_w_new_w, e_splits);
    if regu {
        regu = regularize_face(the_face, &map_old_w_new_w, a_list_of_faces);
    }
    regu
}

/// `FUN_tool_ClassifW(F, mapOwNw, mapWlow)` (`cxx:676-850`).
///
/// `map_wlow` binds each new face's boundary wire to the wires that describe
/// its holes (an empty list when the face has none).
pub fn fun_tool_classif_w(
    f: &TopoShape,
    map_ow_nw: &[(TopoShape, Vec<TopoShape>)],
    map_wlow: &mut WireMap,
) -> bool {
    let tol_f = BRepTool::face_tolerance(&Face(f.clone()));
    let toluv = tol_uv(&Face(f.clone()), tol_f);
    let f_for = f.oriented(Orientation::Forward);
    let mut classi = Classi::new();
    classi.init_2d(&f_for);

    let n_old_w = map_ow_nw.len();

    // noldW = 1 (`cxx:696-701`)
    if n_old_w == 1 {
        let low = &map_ow_nw[0].1;
        // `Classilist` binds keys taken from `low`, so the wire shapes used as
        // keys here are known to be the wires of `low`.
        let mut hm: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
        if !classi.classilist(low, &mut hm) {
            return false;
        }
        for w in low {
            if let Some(holes) = hm.get(&shape_map_key(w)) {
                map_wlow.bind(w, holes.clone());
            }
        }
        return true;
    }

    // noldW > 1 (`cxx:703-846`)
    let mut l_ows: Vec<TopoShape> = Vec::new();
    let null: Vec<TopoShape> = Vec::new();
    for (owi, low) in map_ow_nw {
        l_ows.push(owi.clone());
        let mut lwresu: Vec<TopoShape> = Vec::new();
        fun_add_owlw(owi, low, &mut lwresu);
        for w in lwresu {
            map_wlow.bind(&w, null.clone());
        }
    }

    let mut map_done: Vec<TopoShape> = Vec::new();
    let n_ow = n_old_w;
    // `cxx:715-716`: `nite` is never incremented in OCCT, so the loop ends
    // through the `lOws.Extent() == 0` break or an explicit return.
    let nite = 0;
    let nitemax = (n_ow * (n_ow - 1) / 2) as i32;
    while nite <= nitemax {
        let n_ow = l_ows.len();
        if n_ow == 0 {
            break;
        }

        let ow1 = l_ows[0].clone();
        let isb1 = map_wlow.is_bound(&ow1)
            || !map_done.iter().any(|d| is_equal(d, &ow1));
        if !isb1 {
            continue;
        }

        let lw1 = map_ow_nw
            .iter()
            .find(|(k, _)| is_equal(k, &ow1))
            .map(|(_, l)| l.clone())
            .unwrap_or_default();

        if n_ow == 1 {
            // all wires of <mapWs> have been treated, except the last one
            break;
        }

        let mut out_all = false;
        let mut ow2: Option<TopoShape> = None;
        let mut sta12 = UNKNOWN;
        for cand in l_ows.iter().skip(1) {
            ow2 = Some(cand.clone());
            let isb2 = map_wlow.is_bound(cand)
                || !map_done.iter().any(|d| is_equal(d, cand));
            if !isb2 {
                continue;
            }
            let stageneric = classi.classi_bnd_2d(&ow1, cand, toluv, true);
            sta12 = classi.classip_2d(&ow1, cand, stageneric);
            if sta12 == DIFF {
                out_all = true;
                continue;
            } else if sta12 == UNKNOWN || sta12 == SAME {
                return false;
            }
            break;
        }

        if out_all {
            let mut ldone: Vec<TopoShape> = Vec::new();
            fun_add_owlw(&ow1, &lw1, &mut ldone);
            for w in ldone {
                if !map_done.iter().any(|d| is_equal(d, &w)) {
                    map_done.push(w);
                }
            }
            l_ows.remove(0);
            continue;
        }

        let Some(ow2) = ow2 else {
            // no candidate was classified: `sta12` stayed UNKNOWN.
            return false;
        };
        let lw2 = map_ow_nw
            .iter()
            .find(|(k, _)| is_equal(k, &ow2))
            .map(|(_, l)| l.clone())
            .unwrap_or_default();

        let mut lw1r: Vec<TopoShape> = Vec::new();
        fun_add_owlw(&ow1, &lw1, &mut lw1r);
        let mut lw2r: Vec<TopoShape> = Vec::new();
        fun_add_owlw(&ow2, &lw2, &mut lw2r);
        let mut lgre: Vec<TopoShape> = Vec::new();
        let mut lsma: Vec<TopoShape> = Vec::new();
        if sta12 == ONE_IN_TWO {
            lgre.extend(lw2r.iter().cloned());
            lsma.extend(lw1r.iter().cloned());
        }
        if sta12 == TWO_IN_ONE {
            lgre.extend(lw1r.iter().cloned());
            lsma.extend(lw2r.iter().cloned());
        }

        for wsma in &lsma {
            let isbsma = map_wlow.is_bound(wsma)
                || !map_done.iter().any(|d| is_equal(d, wsma));
            if !isbsma {
                continue;
            }
            for wgre in &lgre {
                let isbgre = map_wlow.is_bound(wgre)
                    || !map_done.iter().any(|d| is_equal(d, wgre));
                if !isbgre {
                    continue;
                }

                let stageneric = classi.classi_bnd_2d(wsma, wgre, toluv, true);
                let sta = classi.classip_2d(wsma, wgre, stageneric);
                if sta == DIFF {
                    continue;
                } else if sta == ONE_IN_TWO {
                    // wsma IN wgre
                    let holes_sma = map_wlow.find(wsma).cloned().unwrap_or_default();
                    if let Some(gre) = map_wlow.change_find(wgre) {
                        gre.extend(holes_sma);
                    }
                    map_wlow.unbind(wsma);
                } else if sta == TWO_IN_ONE {
                    // wgre IN wsma
                    let holes_gre = map_wlow.find(wgre).cloned().unwrap_or_default();
                    if let Some(sma) = map_wlow.change_find(wsma) {
                        sma.extend(holes_gre);
                    }
                    map_wlow.unbind(wgre);
                } else {
                    return false;
                }
            }
        }
        l_ows.remove(0);
    }
    true
}

/// `TopOpeBRepTool::RegularizeFace(theFace, mapoldWnewW, newFaces)`
/// (`cxx:853-887`).
pub fn regularize_face(
    the_face: &TopoShape,
    map_old_w_new_w: &[(TopoShape, Vec<TopoShape>)],
    new_faces: &mut Vec<TopoShape>,
) -> bool {
    let mut map_wlow = WireMap::new();

    let a_face = the_face.oriented(Orientation::Forward);
    if !fun_tool_classif_w(&a_face, map_old_w_new_w, &mut map_wlow) {
        return false;
    }

    // `TopOpeBRepTool_TOOL::WireToFace(theFace, mapWlow, newFaces)`
    let pairs: Vec<(TopoShape, Vec<TopoShape>)> = map_wlow
        .iter()
        .map(|(w, l)| (w.clone(), l.clone()))
        .collect();
    wire_to_face(the_face, &pairs, new_faces);
    true
}
