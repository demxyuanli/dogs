//! `TopOpeBRep_vprdeg.cxx:105-379` — the degenerated-edge data collection that
//! `FillLineVPonR` runs before it processes any VPoint.
//!
//! Sources (OCCT 8.0.0, `ModelingAlgorithms/TKBool/TopOpeBRep`):
//! - `FUN_FillVof12(L, pDS)` (`vprdeg.cxx:105-127`)
//! - `FUN_addmapve(mapve, v, e)` (`cxx:129-177`)
//! - `FUN_GetdgData(pDS, L, F1, F2, datamap)` (`cxx:179-379`)
//!
//! `FUN_GetdgData`'s purpose (`cxx:188-189`) is the map
//! `{vertex -> (closing edge, degenerated edge)}` for shapes of the same rank,
//! plus the same-domain completion `FillShapesSameDomain`.

use std::collections::HashMap;

use crate::abs::ShapeType;
use crate::shape::TopoShape;
use crate::topo_tools::{is_same, map_shapes_and_ancestors, shape_map_key};
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::faces_filler::DegenEdgeData;
use crate::topopebrep::faces_intersector::FacesIntersector;
use crate::topopebrep::tool::shape_tool::closed;

/// `NCollection_DataMap<TopoDS_Shape, int>` with `Bind` / `Find` semantics
/// (`Find` on an unbound key answers `0`, as OCCT's `DataMap::Find` inserts the
/// default value).
struct ShapeRank(Vec<(TopoShape, i32)>);

impl ShapeRank {
    fn new() -> Self {
        Self(Vec::new())
    }

    fn bind(&mut self, s: &TopoShape, rank: i32) {
        match self.0.iter_mut().find(|(k, _)| is_same(k, s)) {
            Some((_, r)) => *r = rank,
            None => self.0.push((s.clone(), rank)),
        }
    }

    fn find(&self, s: &TopoShape) -> i32 {
        self.0
            .iter()
            .find(|(k, _)| is_same(k, s))
            .map(|(_, r)| *r)
            .unwrap_or(0)
    }
}

/// `NCollection_DataMap<TopoDS_Shape, TopoDS_Shape>` with `Bind` / `IsBound` /
/// `Find`.
struct ShapeToShape(Vec<(TopoShape, TopoShape)>);

impl ShapeToShape {
    fn new() -> Self {
        Self(Vec::new())
    }

    fn bind(&mut self, k: &TopoShape, v: &TopoShape) {
        match self.0.iter_mut().find(|(k0, _)| is_same(k0, k)) {
            Some((_, v0)) => *v0 = v.clone(),
            None => self.0.push((k.clone(), v.clone())),
        }
    }

    fn is_bound(&self, k: &TopoShape) -> bool {
        self.0.iter().any(|(k0, _)| is_same(k0, k))
    }

    fn find(&self, k: &TopoShape) -> Option<&TopoShape> {
        self.0.iter().find(|(k0, _)| is_same(k0, k)).map(|(_, v)| v)
    }

    /// `Iterator` over the bindings, in insertion order.
    fn iter(&self) -> impl Iterator<Item = (&TopoShape, &TopoShape)> {
        self.0.iter().map(|(k, v)| (k, v))
    }
}

/// `FUN_FillVof12(L, pDS)` (`vprdeg.cxx:105-127`): every VPoint on both
/// surfaces that is a vertex on each side makes its two vertices same-domain.
pub fn fun_fill_v_of12(
    hds: &mut HDataStructure,
    facint: &mut FacesIntersector,
    line_index: i32,
) -> Result<(), String> {
    let nb = facint.change_line(line_index).nb_v_point();
    for k in 1..=nb {
        let vp = match facint.change_line(line_index).v_point(k) {
            Some(v) => v.clone(),
            None => continue,
        };
        if vp.shape_index() != 3 {
            continue;
        }
        let isvon1 = vp.is_vertex_on_s1()?;
        let isvon2 = vp.is_vertex_on_s2()?;
        if !(isvon1 && isvon2) {
            continue;
        }
        let v1 = vp.vertex_on_s1()?.clone();
        let v2 = vp.vertex_on_s2()?.clone();
        hds.change_ds().fill_shapes_same_domain(&v1, &v2, true)?;
    }
    Ok(())
}

/// `FUN_addmapve(mapve, v, e)` (`vprdeg.cxx:129-177`): the symmetric binding
/// `{v -> (e...), e -> (v...)}`.
fn fun_add_map_ve(mapve: &mut DegenEdgeData, v: &TopoShape, e: &TopoShape) {
    let visb = mapve.is_bound(v);
    let eisb = mapve.is_bound(e);
    if !visb && !eisb {
        mapve.bind(v, std::slice::from_ref(e));
        mapve.bind(e, std::slice::from_ref(v));
    } else if visb && !eisb {
        mapve.change_append(v, e);
        mapve.bind(e, std::slice::from_ref(v));
    } else if !visb && eisb {
        mapve.change_append(e, v);
        mapve.bind(v, std::slice::from_ref(e));
    } else {
        let found = mapve
            .find(v)
            .map(|le| le.iter().any(|x| is_same(x, e)))
            .unwrap_or(false);
        if !found {
            mapve.change_append(v, e);
            mapve.change_append(e, v);
        }
    }
}

/// `FUN_GetdgData(pDS, L, F1, F2, datamap)` (`vprdeg.cxx:179-379`).
pub fn fun_get_dg_data(
    hds: &mut HDataStructure,
    facint: &mut FacesIntersector,
    line_index: i32,
    f1: &TopoShape,
    f2: &TopoShape,
    datamap: &mut DegenEdgeData,
) -> Result<(), String> {
    let mut shaperk = ShapeRank::new(); // {shape, rank=1,2}
    let mut mapvec = DegenEdgeData::default(); // {(v,lec),(ec,lv)}
    let mut mapved = DegenEdgeData::default(); // {(v,led),(ed,lv)}
    let mut mapvvsd = ShapeToShape::new(); // {(v,vsd),(vsd,v)}

    let nb = facint.change_line(line_index).nb_v_point();
    for k in 1..=nb {
        let vp = match facint.change_line(line_index).v_point(k) {
            Some(v) => v.clone(),
            None => continue,
        };
        let isv1 = vp.is_vertex(1)?;
        let isv2 = vp.is_vertex(2)?;
        if !(isv1 || isv2) {
            continue;
        }

        let sind = vp.shape_index();
        let mut v = if isv1 {
            vp.vertex(1).cloned()
        } else {
            vp.vertex(2).cloned()
        };
        for i in 1..=2 {
            let f = if i == 1 { f1 } else { f2 };

            if vp.is_vertex(i)? {
                v = vp.vertex(i).cloned();
                if let Some(vv) = &v {
                    shaperk.bind(vv, i);
                }
            }

            let ison = (sind == i) || (sind == 3);
            if !ison {
                continue;
            }
            let Some(e) = vp.edge(i).cloned() else {
                continue;
            };
            shaperk.bind(&e, i);

            let isdg = crate::brep_tool::BRepTool::is_degenerated(&crate::shape::Edge(e.clone()));
            let iscl = if isdg { false } else { closed(&e, f) };
            if let Some(vv) = &v {
                if isdg {
                    fun_add_map_ve(&mut mapved, vv, &e);
                }
                if iscl {
                    fun_add_map_ve(&mut mapvec, vv, &e);
                }
            }
        }
    }

    // filling up map mapvvsd (`cxx:248-275`)
    for (key, value) in mapved.iter() {
        if key.shape_type() != ShapeType::Vertex {
            continue;
        }
        let v = key;
        let rkv = shaperk.find(v);
        for e in value {
            let rke = shaperk.find(e);
            if rke != rkv {
                // `TopExp_Explorer ex(e, TopAbs_VERTEX); ex.Current()`
                let Some(vsd) = first_vertex(e) else {
                    continue;
                };
                mapvvsd.bind(v, &vsd);
                mapvvsd.bind(&vsd, v);
            }
        }
    }

    // `cxx:277-361`: the degenerated edge arm
    for (key, value) in mapved.iter() {
        let dge = key;
        let rk = shaperk.find(dge);
        let f = if rk == 1 { f1 } else { f2 };
        if dge.shape_type() != ShapeType::Edge {
            continue;
        }
        let Some(v) = first_vertex(dge) else {
            continue;
        };
        let rkv = shaperk.find(&v);
        let hassd = mapvvsd.is_bound(&v);
        let vsd: Option<TopoShape> = if hassd {
            mapvvsd.find(&v).cloned()
        } else {
            None
        };

        let mut hasecl = false;
        let mut cle: Option<TopoShape> = None;
        let mut isbv = mapvec.is_bound(&v);
        let isbvsd = match &vsd {
            Some(vs) => mapvec.is_bound(vs),
            None => false,
        };
        if !isbv && !isbvsd {
            // `cxx:301-330`: look for a closing edge of `f` through `v`.
            let mapve: HashMap<_, Vec<TopoShape>> =
                map_shapes_and_ancestors(f, ShapeType::Vertex, ShapeType::Edge);
            if let Some(iteds) = mapve.get(&shape_map_key(&v)) {
                for ee in iteds {
                    if is_same(ee, dge) {
                        continue;
                    }
                    if !closed(ee, f) {
                        continue;
                    }
                    isbv = true;
                    cle = Some(ee.clone());
                    hasecl = true;
                    break;
                }
            }
        }
        if !hasecl && (isbv || isbvsd) {
            // `cxx:330-349`: `ite.Initialize(mapvec.Find(v))` — OCCT reads the
            // list keyed by `v` on both arms even when only `vsd` is bound, and
            // the local `vv` it computes there is shadowed unused at `cxx:354`.
            if let Some(list) = mapvec.find(&v) {
                for e in list {
                    let rke = shaperk.find(e);
                    if rke == rk {
                        cle = Some(e.clone());
                        hasecl = true;
                        break;
                    }
                }
            }
        }
        if !hasecl {
            continue;
        }
        let vv = if rkv == rk {
            v.clone()
        } else {
            match vsd.clone() {
                Some(s) => s,
                None => continue,
            }
        };
        let cle = match cle {
            Some(c) => c,
            None => continue,
        };
        datamap.bind(&vv, &[cle, dge.clone()]);
    }

    // filling sdm shapes (`cxx:363-377`)
    for (_, value) in mapvvsd.iter() {
        let v = value;
        let Some(ov) = mapvvsd.find(v).cloned() else {
            continue;
        };
        let rkv = shaperk.find(v);
        let v1 = if rkv == 1 { v.clone() } else { ov.clone() };
        let v2 = if rkv == 2 { v.clone() } else { ov };
        hds.change_ds().fill_shapes_same_domain(&v1, &v2, true)?;
    }
    Ok(())
}

/// `TopExp_Explorer(e, TopAbs_VERTEX)`'s `Current()`.
fn first_vertex(e: &TopoShape) -> Option<TopoShape> {
    let mut ex = crate::topexp::Explorer::new(e, ShapeType::Vertex);
    if ex.more() {
        Some(ex.current().clone())
    } else {
        None
    }
}
