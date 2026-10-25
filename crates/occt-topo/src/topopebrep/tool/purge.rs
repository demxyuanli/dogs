//! `TopOpeBRepTool_PURGE.cxx` — wire closure diagnostics used before faces
//! are rebuilt from a `WireEdgeSet`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepTool/TopOpeBRepTool_PURGE.cxx`
//! (OCCT 8.0.0).
//!
//! Ported:
//! - `FUN_getv` (`cxx:75-101`)
//! - `FUN_addtomap` (`cxx:58-73`)
//! - `FUN_tool_ClosedW` (`cxx:103-296`), the wire closure test called by
//!   `Builder::GWESMakeFaces` (`TopOpeBRepBuild_makefaces.cxx:685`) through
//!   `Builder::MakeWire`'s `Closed` flag.
//! - `TopOpeBRepTool::PurgeClosingEdges(Fin, FF, MWisOld, MshNOK)`
//!   (`cxx:300-401`) as [`purge_closing_edges`] and its `LOF` overload
//!   (`cxx:405-427`) as [`purge_closing_edges_faces`]; `MWisOld` is unused in
//!   OCCT (`cxx:305-306`).
//! - `FUN_connexX` (`cxx:677-770`) as [`fun_connex_x`].
//! - `TopOpeBRepTool::CorrectONUVISO` (`cxx:774-970`) as
//!   [`correct_on_uv_iso`].
//! - `TopOpeBRepTool::MakeFaces` (`cxx:974-1038`) as [`make_faces`].
//!
//! UNPORTED:
//! - `FUN_correctDegeneratedE` (`cxx:454-559`) and `FUN_tool_reguUV`
//!   (`cxx:561-674`) are commented out in OCCT (`/* ... */`).
//! - `FUN_correctClosingE` (`cxx:431-452`) is commented out in OCCT.
//! - `FUN_translate` / `FUN_tool_ttranslate` / `FUN_tool_translate`
//!   (`cxx:1040-1073`) are commented out in OCCT.

use std::collections::{HashMap, HashSet};

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{explorer_shapes, is_same_occurrence, shape_map_key, ShapeMapKey};
use crate::topopebrep::tool::closed::{par_e, vertices};
use crate::topopebrep::tool::corr_iso::{bind_fy, CorrIso, DECREASE, INCREASE};
use crate::topopebrep::tool::iso::{iso_on_clo_c2df, uvf};

/// `FUN_getv(orivine, e, v)` (`cxx:75-101`): the vertex of `e` whose stored
/// orientation is `orivine`, reading `e` `FORWARD` when `e` itself is
/// `INTERNAL`/`EXTERNAL` (`cxx:79-89`).
pub fn fun_getv(orivine: Orientation, e: &TopoShape) -> Option<TopoShape> {
    let oe = e.orientation();
    let eo = if oe == Orientation::Internal || oe == Orientation::External {
        e.oriented(Orientation::Forward)
    } else {
        e.clone()
    };
    explorer_shapes(&eo, ShapeType::Vertex)
        .into_iter()
        .find(|v| v.orientation() == orivine)
}

/// `FUN_addtomap(map, key, item)` (`cxx:58-73`).
fn fun_addtomap(map: &mut HashMap<ShapeMapKey, Vec<TopoShape>>, key: &TopoShape, item: &TopoShape) {
    map.entry(shape_map_key(key))
        .or_default()
        .push(item.clone());
}

/// `FUN_tool_ClosedW(W)` (`cxx:103-296`).
///
/// A wire is closed when its `FORWARD`/`REVERSED` vertices pair up, and every
/// vertex that appears on a single oriented edge is cancelled either by an
/// `INTERNAL` vertex of the same point or by the opposite orientation on one
/// other edge. An unoriented (`INTERNAL`/`EXTERNAL`) edge makes every one of
/// its vertices fail the test (`cxx:249-253`).
pub fn fun_tool_closed_w(w: &TopoShape) -> bool {
    // `cxx:110-111`
    let mut map_v_fine: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
    let mut map_v_rine: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
    let mut map_v_iine: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();

    // `cxx:113-147`
    for e in explorer_shapes(w, ShapeType::Edge) {
        let oe = e.orientation();
        let is_non_o = oe == Orientation::Internal || oe == Orientation::External;
        let eo = if is_non_o {
            e.oriented(Orientation::Forward)
        } else {
            e.clone()
        };
        for v in explorer_shapes(&eo, ShapeType::Vertex) {
            let oriv = v.orientation();
            if oriv == Orientation::Forward {
                fun_addtomap(&mut map_v_fine, &v, &e);
            }
            if oriv == Orientation::Reversed {
                fun_addtomap(&mut map_v_rine, &v, &e);
            }
            if oriv == Orientation::Internal {
                fun_addtomap(&mut map_v_iine, &v, &e);
            }
        }
    }

    // `cxx:149-152`
    if map_v_fine.is_empty() {
        return false;
    }

    // `cxx:154-197`
    let mut map_v_ok: HashSet<ShapeMapKey> = HashSet::new();
    let mut map_v_only_fine: HashSet<ShapeMapKey> = HashSet::new();
    let fine_keys: Vec<ShapeMapKey> = map_v_fine.keys().cloned().collect();
    for key in fine_keys {
        if map_v_iine.contains_key(&key) {
            map_v_ok.insert(key);
            continue;
        }
        if map_v_rine.contains_key(&key) {
            map_v_ok.insert(key);
            continue;
        }
        map_v_only_fine.insert(key);
    }
    let mut map_v_only_rine: HashSet<ShapeMapKey> = HashSet::new();
    let rine_keys: Vec<ShapeMapKey> = map_v_rine.keys().cloned().collect();
    for key in rine_keys {
        if map_v_ok.contains(&key) {
            continue;
        }
        if map_v_iine.contains_key(&key) {
            map_v_ok.insert(key);
            continue;
        }
        map_v_only_rine.insert(key);
    }

    // `cxx:199-294`: `nmap` walks the FORWARD-only set first, then the
    // REVERSED-only one.
    for v_fine in [true, false] {
        let to_check: Vec<ShapeMapKey> = if v_fine {
            map_v_only_fine.iter().cloned().collect()
        } else {
            map_v_only_rine.iter().cloned().collect()
        };
        for key in to_check {
            // `cxx:218-232`
            let eds = if v_fine {
                map_v_fine.get(&key)
            } else {
                map_v_rine.get(&key)
            };
            let eds = match eds {
                Some(l) => l,
                None => continue,
            };
            if eds.len() > 1 {
                return false;
            }
            let e = &eds[0];

            // `cxx:234-240`
            let ov_ori = if v_fine {
                Orientation::Reversed
            } else {
                Orientation::Forward
            };
            let ov = match fun_getv(ov_ori, e) {
                Some(v) => v,
                None => return false,
            };

            // `cxx:249-253`
            let oe = e.orientation();
            if oe == Orientation::Forward || oe == Orientation::Reversed {
                return false;
            }
            // `cxx:254-257`
            let ov_key = shape_map_key(&ov);
            if !map_v_ok.contains(&ov_key) {
                return false;
            }
            // `cxx:259-263`
            if map_v_iine.contains_key(&ov_key) {
                continue;
            }
            // `cxx:264-268`
            if !map_v_rine.contains_key(&ov_key) {
                return false;
            }
            // `cxx:269-273`
            if !map_v_rine.contains_key(&ov_key) {
                return false;
            }
            // `cxx:275-284`
            let eds_ov_fine = map_v_fine.get(&ov_key).map(|l| l.len()).unwrap_or(0);
            let eds_ov_rine = map_v_rine.get(&ov_key).map(|l| l.len()).unwrap_or(0);
            if eds_ov_fine > 1 {
                continue;
            }
            if eds_ov_rine > 1 {
                continue;
            }
            // `cxx:285-292`
            if let Some(first) = map_v_fine.get(&ov_key).and_then(|l| l.first()) {
                if crate::topo_tools::is_equal(first, e) {
                    return false;
                }
            }
            if let Some(first) = map_v_rine.get(&ov_key).and_then(|l| l.first()) {
                if crate::topo_tools::is_equal(first, e) {
                    return false;
                }
            }
        }
    }
    true
}

/// `MshNOK.Add(S)` (`NCollection_IndexedMap<Shape>::Add`): a shape already in
/// the list keeps its index, so a repeat is dropped.
fn indexed_map_add(m: &mut Vec<TopoShape>, s: &TopoShape) {
    let key = shape_map_key(s);
    if m.iter().any(|x| shape_map_key(x) == key) {
        return;
    }
    m.push(s.clone());
}

/// `TopOpeBRepTool::PurgeClosingEdges(Fin, FF, MWisOld, MshNOK)`
/// (`cxx:300-401`).
///
/// `Fin` is the ancestor face, `FF` a split face of it. A wire of `FF` whose
/// pcurves do not close up along the periodic direction that `Fin` closes in is
/// flagged: the faulty closing edges, the wire and `FF` itself land in
/// `msh_nok`, which [`make_faces`] then reads to rebuild `FF`.
///
/// `MWisOld` is unnamed and unused in OCCT (`cxx:305-306`).
pub fn purge_closing_edges(
    fin: &TopoShape,
    ff: &TopoShape,
    msh_nok: &mut Vec<TopoShape>,
) -> bool {
    let corr = CorrIso::new(fin);
    let tol_f = BRepTool::face_tolerance(&Face(fin.clone()));

    // `cxx:325-328`
    let uperiod = corr.ref_closed(1);
    let vperiod = corr.ref_closed(2);
    let uclosed = uperiod.is_some();
    let vclosed = vperiod.is_some();
    if !uclosed && !vclosed {
        return false;
    }

    // `cxx:330-334`
    let in_u = uclosed;
    let gas = corr.gas_ref();
    let xmin = if in_u {
        gas.first_u_parameter()
    } else {
        gas.first_v_parameter()
    };
    let xper = if in_u {
        uperiod.unwrap_or(0.0)
    } else {
        vperiod.unwrap_or(0.0)
    };
    let tolx = if in_u {
        corr.tol(1, tol_f)
    } else {
        corr.tol(2, tol_f)
    };

    for w in explorer_shapes(ff, ShapeType::Wire) {
        // `cxx:342-347`
        let mut corr = CorrIso::new(fin);
        corr.init(&w);
        if corr.uv_closed() {
            continue;
        }

        // `cxx:349-373`
        let mut c_eds: Vec<TopoShape> = Vec::new();
        let eds = corr.eds().to_vec();
        for e in eds.iter() {
            let mut closing =
                BRepTool::is_closed_edge_face(&Edge(e.clone()), &Face(fin.clone()));
            if !closing {
                // `cxx:357-358` (`xpu231198`): the pcurve was modified, the
                // `IsClosed` flag on the edge no longer knows.
                let Some(c2df) = corr.uv_rep(e) else {
                    return false;
                };
                if iso_on_clo_c2df(&c2df, in_u, xmin, xper, tolx) {
                    closing = true;
                }
            }
            if closing {
                c_eds.push(e.clone());
            }
        }

        // `cxx:375-379`
        let nc_e = c_eds.len();
        let nopurge = nc_e <= 1;
        if nopurge {
            return true;
        }

        // `cxx:381-395`
        if let Some(fy_eds) = corr.purge_fy_closing_e(&c_eds) {
            for f in fy_eds.iter() {
                indexed_map_add(msh_nok, f);
            }
            indexed_map_add(msh_nok, &w);
            indexed_map_add(msh_nok, ff);
        }
    }
    true
}

/// `TopOpeBRepTool::PurgeClosingEdges(Fin, LOF, MWisOld, MshNOK)`
/// (`cxx:405-427`): the `FF` loop over a list of split faces.
pub fn purge_closing_edges_faces(
    fin: &TopoShape,
    lof: &[TopoShape],
    msh_nok: &mut Vec<TopoShape>,
) -> bool {
    let closed = crate::topopebrep::tool::geometry::fun_tool_is_closed_s(&Face(fin.clone()));
    if !closed {
        return true;
    }
    for ff in lof.iter() {
        if !purge_closing_edges(fin, ff, msh_nok) {
            return false;
        }
    }
    true
}

/// `FUN_connexX(onU, CORRISO, EdstoCheck, fyEds)` (`cxx:677-770`).
///
/// `Fref` is `x`-periodic. Finds the edge with **both** bounds faulty and, for
/// each of its vertices, the connex edge whose vertex sits one period away in
/// `x`; that edge is what must be translated, and the sign says which way.
pub fn fun_connex_x(
    on_u: bool,
    corr: &CorrIso,
    eds_to_check: &[TopoShape],
) -> Vec<(TopoShape, i32)> {
    let mut fy_eds: Vec<(TopoShape, i32)> = Vec::new();

    // `cxx:687-693`
    let tol_f = BRepTool::face_tolerance(&Face(corr.f_ref().clone()));
    let index = if on_u { 1 } else { 2 };
    let Some(xperiod) = corr.ref_closed(index) else {
        return fy_eds;
    };
    let xtol = corr.tol(index, tol_f);

    // `cxx:699-702`: `fy` has its two UV bounds non-connexed.
    let Some((fy, _ify)) = corr.edge_with_faulty_uv_first(eds_to_check, 2) else {
        return fy_eds;
    };
    let Some(c2df) = corr.uv_rep(&fy) else {
        return fy_eds;
    };

    let (vfy1, vfy2) = vertices(&fy);
    for ii in 1..=2 {
        // `cxx:714-719`
        let vff = if ii == 1 { vfy1.clone() } else { vfy2.clone() };
        let Some(vff) = vff else {
            continue;
        };
        let Some(parvff) = par_e(ii, &fy) else {
            continue;
        };
        let Some(uvff) = uvf(parvff, &c2df) else {
            continue;
        };

        // `cxx:721-726`
        let Some(loe) = corr.connexity(&vff) else {
            return fy_eds;
        };

        for ee in loe.iter() {
            let (vee1, vee2) = vertices(ee);
            for ive in 1..=2 {
                let ve = if ive == 1 { vee1.clone() } else { vee2.clone() };
                let Some(ve) = ve else {
                    continue;
                };
                if !is_same_occurrence(&ve, &vff) {
                    continue;
                }
                if ive == ii {
                    continue;
                }
                let Some(c2dfe) = corr.uv_rep(ee) else {
                    return fy_eds;
                };
                let Some(paree) = par_e(ive, ee) else {
                    continue;
                };
                let Some(uve) = uvf(paree, &c2dfe) else {
                    continue;
                };

                // `cxx:756-767` (`xxtrsl`)
                let dxx = if on_u {
                    uve.x() - uvff.x()
                } else {
                    uve.y() - uvff.y()
                };
                let isper = (xperiod - dxx.abs()).abs() < xtol;
                if !isper {
                    continue;
                }

                let recadre = if dxx > 0.0 { INCREASE } else { DECREASE };
                bind_fy(&mut fy_eds, &fy, recadre);
            }
        }
    }
    fy_eds
}

/// `TopOpeBRepTool::CorrectONUVISO(Fin, Fsp)` (`cxx:774-970`).
///
/// `Fin` is the face on the `x`-periodic surface, `Fsp` a split face of it.
/// First tries to UV-connect the edges that sit one period away from their
/// neighbour ([`fun_connex_x`]); then translates the pcurves that fall outside
/// the `Fin` window. Either way, on success `Fsp` gets the corrected
/// representations and `true` comes back.
pub fn correct_on_uv_iso(fin: &TopoShape, fsp: &mut TopoShape) -> bool {
    let tol_f = BRepTool::face_tolerance(&Face(fin.clone()));
    let mut corr = CorrIso::new(fin);

    // `cxx:787-792`
    let uperiod = corr.ref_closed(1);
    let vperiod = corr.ref_closed(2);
    if uperiod.is_none() && vperiod.is_none() {
        return false;
    }

    // `cxx:794-799`
    corr.init(fsp);
    if corr.uv_closed() {
        return true;
    }

    // `cxx:804-828`: 1. connexity among all edges of `Fsp`.
    for i in 1..=2 {
        let on_u = i == 1;
        let to_check = corr.eds().to_vec();
        let fy_eds = fun_connex_x(on_u, &corr, &to_check);
        if fy_eds.is_empty() {
            continue;
        }
        if !corr.trsl_uv(on_u, &fy_eds) {
            continue;
        }
        if !corr.uv_closed() {
            continue;
        }
        return match corr.get_new_s() {
            Some(new_s) => {
                *fsp = new_s;
                true
            }
            None => false,
        };
    }

    // `cxx:831-860`: 2. the `x`-2D representations are in
    // `[xfirst, xfirst + xperiod]`.
    for i in 1..=2 {
        let on_u = i == 1;
        let Some(_xper) = corr.ref_closed(i) else {
            continue;
        };
        let mut tolx = corr.tol(i, tol_f);
        tolx *= 1.0e2; // `cxx:841` (BUC60380)
        let eds = corr.eds().to_vec();
        let fy_eds = corr.edges_out_of_bounds_uv(&eds, on_u, tolx);
        if fy_eds.is_empty() {
            continue;
        }
        if !corr.trsl_uv(on_u, &fy_eds) {
            return false;
        }
        if !corr.uv_closed() {
            continue;
        }
        return match corr.get_new_s() {
            Some(new_s) => {
                *fsp = new_s;
                true
            }
            None => false,
        };
    }
    false
}

/// `TopOpeBRepTool::MakeFaces(Fin, LOF, MshNOK, LOFF)` (`cxx:974-1038`).
///
/// Every face of `lof` that [`purge_closing_edges`] did not flag is kept as is;
/// a flagged face is rebuilt from `Fin`'s surface ([`TopoDS_Shape::EmptyCopied`],
/// `cxx:992-993`), dropping the flagged wires and, inside a flagged wire, the
/// flagged edges. A wire emptied that way is dropped (`cxx:1026-1029`).
pub fn make_faces(fin: &TopoShape, lof: &[TopoShape], msh_nok: &[TopoShape]) -> Vec<TopoShape> {
    let mut loff: Vec<TopoShape> = Vec::new();
    let bb = TopoBuilder::new();
    let is_nok = |s: &TopoShape| msh_nok.iter().any(|x| shape_map_key(x) == shape_map_key(s));

    for ff in lof.iter() {
        let valid = !is_nok(ff);
        if valid {
            loff.push(ff.clone());
            continue;
        }

        // `cxx:992-993`: `Fin.EmptyCopied()` keeps the surface, the location
        // and the tolerance (`BRep_TFace.cxx:36-43`), not the wires.
        let Some(surface) = BRepTool::face_surface_world(&Face(fin.clone())) else {
            continue;
        };
        let tol = BRepTool::face_tolerance(&Face(fin.clone()));
        let mut new_face = bb.make_face_with_tolerance(surface, &[], tol).0;
        new_face.set_location(fin.location());

        for w in explorer_shapes(ff, ShapeType::Wire) {
            let valid = !is_nok(&w);
            if valid {
                bb.add(&mut new_face, &w);
                continue;
            }

            let mut new_wire = bb.make_shape(ShapeType::Wire);
            let mut ne = 0;
            for e in explorer_shapes(&w, ShapeType::Edge) {
                if is_nok(&e) {
                    continue;
                }
                bb.add(&mut new_wire, &e);
                ne += 1;
            }
            if ne == 0 {
                continue;
            }
            let closed = fun_tool_closed_w(&new_wire);
            new_wire.set_closed(closed);
            bb.add(&mut new_face, &new_wire);
        }
        loff.push(new_face);
    }
    loff
}
