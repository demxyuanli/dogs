//! File-local statics of `LocOpe_Spliter.cxx`.
//!
//! | OCCT | lines |
//! |---|---|
//! | `static void RebuildWires(Ledges, PW)` | `620-631` |
//! | `static void Put(S, theMap)` | `634-652` |
//! | `static void Select(Ebase, lsubs)` | `657-719` |

use crate::brep_tool::BRepTool;
use crate::locope::build_wires::BuildWires;
use crate::locope::shape_maps::ShapeValueMap;
use crate::locope::wires_on_shape::WiresOnShape;
use crate::shape::{Edge, TopoShape};
use crate::topo_tools::children_keep_orientation;

/// `static void Put(S, theMap)` (`cxx:634-652`): bind `S` to `[S]` and recurse
/// into its children with `TopoDS_Iterator` (orientation of the child kept).
pub(super) fn put(s: &TopoShape, the_map: &mut ShapeValueMap<Vec<TopoShape>>) {
    // `cxx:639-642`.
    if the_map.is_bound(s) {
        return;
    }
    // `cxx:643-645`: `Bind(S, list); theMap(S).Append(S)`.
    the_map.bind(s, vec![s.clone()]);
    // `cxx:646-650`.
    for child in children_keep_orientation(s) {
        put(&child, the_map);
    }
}

/// `static void RebuildWires(ledge, PW)` (`cxx:620-631`): replace `ledge` by the
/// wires `LocOpe_BuildWires` chains out of it. The `Standard_ConstructionError`
/// of `cxx:628-631` becomes the `Err`.
pub(super) fn rebuild_wires(
    ledges: &[TopoShape],
    pw: &WiresOnShape,
) -> Result<Vec<TopoShape>, String> {
    let mut the_build = BuildWires::new();
    the_build.perform(ledges, pw)?;
    if !the_build.is_done() {
        return Err("LocOpe_Spliter: LocOpe_BuildWires not done".to_string());
    }
    Ok(the_build.result().unwrap_or(&[]).to_vec())
}

/// `static void Select(Ebase, lsubs)` (`cxx:657-719`): keep the single edge of
/// `lsubs` whose point at parameter `(f+l)/2` is the closest to the midpoint of
/// `Ebase`; an empty answer means "no candidate projects" (`cxx:699-702`).
pub(super) fn select(ebase: &Edge, lsubs: &mut Vec<TopoShape>) {
    // `cxx:665-674`: the reference curve, transformed by its location.
    let Some(curve) = BRepTool::edge_curve_world(ebase) else {
        lsubs.clear();
        return;
    };
    let (f, l) = BRepTool::edge_parameters(ebase);
    // `cxx:675`: the point at the middle parameter.
    let pt = curve.d0((f + l) / 2.0);

    // `cxx:677-698`: nearest projection over each candidate's own range.
    let mut dmin = f64::MAX;
    let mut imin = 0usize;
    for (i, s) in lsubs.iter().enumerate() {
        let Some(edg) = Edge::wrap(s.clone()) else {
            continue;
        };
        let Some(c) = BRepTool::edge_curve_world(&edg) else {
            continue;
        };
        let (ef, el) = BRepTool::edge_parameters(&edg);
        // `proj.Init(Pt, C, f, l)` + `NbPoints() > 0` + `LowerDistance()`.
        if let Some(e) = occt_geom::extrema_pc::point_curve_extrema_range(c.as_ref(), &pt, ef, el) {
            if e.distance < dmin {
                imin = i + 1;
                dmin = e.distance;
            }
        }
    }

    // `cxx:699-718`.
    if imin == 0 {
        lsubs.clear();
    } else {
        let keep = lsubs[imin - 1].clone();
        lsubs.clear();
        lsubs.push(keep);
    }
}
