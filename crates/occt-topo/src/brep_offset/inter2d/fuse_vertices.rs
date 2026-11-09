//! `BRepOffset_Inter2d::FuseVertices` (`BRepOffset_Inter2d.cxx:2335-2393`)
//! with its file-local helper `MakeChain` (`:2311-2332`).
//!
//! `Store` records every pair of coincident intersection vertices in
//! `theDMVV` (`cxx:366-379`); `FuseVertices` walks those lists into connected
//! chains (`MakeChain`) and replaces every chain of two or more vertices by a
//! single vertex built with `BOPTools_AlgoTools::MakeVertex(list)`
//! (`cxx:2360`):
//!
//! 1. for every old member, every edge of `theAsDes->Ascendant(old)`
//!    (`cxx:2370`) gets that member's parameter written with
//!    `BRep_Builder::UpdateVertex(new, T, E, TolE)` (`cxx:2381`),
//! 2. `theAsDes->Replace(old, new)` (`cxx:2384`),
//! 3. when `old` has an image in `theImageVV`, that image gains `new`
//!    (`cxx:2386-2390`).
//!
//! `BRep_Tool::Parameter(V, E, T)` answering `false` (`cxx:2379-2380`) is the
//! `None` of [`BRepTool::parameter_of_vertex_checked`]; the reference returns
//! `false` there, which the port does as well.
//!
//! Port note: `MakeChain`'s `Seek` is [`Dmvv::seek`] and the `FindKey(i)` walk
//! of `cxx:2345` is [`Dmvv::keys`], which keeps OCCT's key order.

use std::collections::HashSet;

use crate::abs::Orientation;
use crate::algo_tools_range::make_vertex_from_list;
use crate::brep_algo::AsDes;
use crate::brep_algo_image::BRepAlgoImage;
use crate::brep_builder_update_vertex::update_vertex_on_edge;
use crate::brep_tool::BRepTool;
use crate::shape::TopoShape;
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{shape_map_key, ShapeMapKey};

use super::Dmvv;

/// `MakeChain(theV, theDMVV, theMDone, theChain)` (`cxx:2311-2332`): append
/// `theV` and, recursively, every vertex reachable from it through `theDMVV`,
/// visiting each vertex once (`theMDone`).
fn make_chain(
    v: &TopoShape,
    dmvv: &Dmvv,
    done: &mut HashSet<ShapeMapKey>,
    chain: &mut Vec<TopoShape>,
) {
    // `if (theMDone.Add(theV))` (`cxx:2318`).
    if !done.insert(shape_map_key(v)) {
        return;
    }
    // `theChain.Append(theV);` (`cxx:2320`).
    chain.push(v.clone());
    // `theDMVV.Seek(theV)` (`cxx:2321-2329`); an unbound key answers the null
    // pointer, i.e. nothing to recurse into.
    if let Some(list) = dmvv.seek(v) {
        for value in list {
            make_chain(value, dmvv, done, chain);
        }
    }
}

/// `BRepOffset_Inter2d::FuseVertices(theDMVV, theAsDes, theImageVV)`
/// (`cxx:2335-2393`).
///
/// `false` is the reference's early return at `cxx:2380` (no parameter for an
/// ascendant edge of an old vertex).
pub fn fuse_vertices(
    dmvv: &Dmvv,
    as_des: &mut AsDes,
    image_vv: &mut BRepAlgoImage,
) -> bool {
    let reg = GeometryRegistry::global();
    // `NCollection_Map<TopoDS_Shape> aMVDone;` (`cxx:2344`).
    let mut mv_done: HashSet<ShapeMapKey> = HashSet::new();

    // `for (i = 1; i <= aNb; ++i)` (`cxx:2346`).
    for i in 0..dmvv.extent() {
        let v_key = dmvv.keys()[i].clone();

        // `cxx:2349-2351`: the chain of vertices.
        let mut chain: Vec<TopoShape> = Vec::new();
        make_chain(&v_key, dmvv, &mut mv_done, &mut chain);
        // `if (aLVChain.Extent() < 2) continue;` (`cxx:2353-2356`).
        if chain.len() < 2 {
            continue;
        }

        // `BOPTools_AlgoTools::MakeVertex(aLVChain, aVNew)` (`cxx:2360`).
        let Ok(v_new) = make_vertex_from_list(&chain) else {
            // Port guard: `MakeVertex` cannot fail for a non-empty chain.
            return false;
        };
        // `TopoDS::Vertex(aVNew.Oriented(TopAbs_INTERNAL))` (`cxx:2362`).
        let v_new_int = v_new.oriented(Orientation::Internal);

        // `cxx:2365-2390`: every old member.
        for v_old in &chain {
            // `TopoDS::Vertex(aVOld.Oriented(TopAbs_INTERNAL))` (`cxx:2369`).
            let v_old_int = v_old.oriented(Orientation::Internal);
            // `theAsDes->Ascendant(aVOld)` (`cxx:2370`), copied because
            // `Replace` below rebinds the lists.
            let edges: Vec<TopoShape> = as_des.ascendant(v_old).to_vec();

            // `cxx:2372-2383`.
            for edge in &edges {
                let tol_e = reg.edge_tolerance(edge);
                // `if (!BRep_Tool::Parameter(aVOldInt, aE, aT)) return false;`
                // (`cxx:2379-2380`).
                let Some(t) = BRepTool::parameter_of_vertex_checked(edge, &v_old_int) else {
                    return false;
                };
                // `aBB.UpdateVertex(aVNewInt, aT, aE, aTolE);` (`cxx:2381`).
                let _ = update_vertex_on_edge(&v_new_int, t, edge, tol_e);
            }

            // `theAsDes->Replace(aVOld, aVNew);` (`cxx:2384`).
            as_des.replace(v_old, &v_new);

            // `cxx:2386-2390`: the image chain of `aVOld` gains `aVNew`.
            if image_vv.is_image(v_old) {
                if let Ok(pro_vertex) = image_vv.image_from(v_old) {
                    let _ = image_vv.add(&pro_vertex, &v_new.oriented(Orientation::Forward));
                }
            }
        }
    }
    true
}
