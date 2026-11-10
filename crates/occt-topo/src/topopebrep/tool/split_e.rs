//! `TopOpeBRepTool_TOOL::SplitE` and the two helpers it needs.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepTool/TopOpeBRepTool_TOOL.cxx:436-469`
//! (`FUN_tool_sortVonE`), `:473-528` (`SplitE`), `:1319-1335` (`TrslUVModifE`)
//! and `TopOpeBRepTool_TOPOLOGY.cxx:1573-1623` (`FUN_ds_CopyEdge`,
//! `FUN_ds_Parameter`) (OCCT 8.0.0).
//!
//! `SplitE` cuts an edge that carries more than two vertices into one edge per
//! consecutive vertex pair: the vertices are sorted by parameter, then each
//! interval becomes a copy of the edge (`FUN_ds_CopyEdge`) with the two
//! vertices re-added at their parameters (`FUN_ds_Parameter`), oriented like
//! the input edge.
//!
//! The two `FUN_ds_*` helpers live in `TopOpeBRepTool_TOPOLOGY.cxx` only
//! because `TopOpeBRepDS_BuildTool` would close a dependency cycle
//! (`cxx:1569-1570`); their bodies are `BuildTool::CopyEdge`
//! (`TopOpeBRepDS_BuildTool.cxx:293-317`) and `BuildTool::Parameter`
//! (`:950-993`), so the port calls those instead of duplicating them.

use std::sync::Arc;

use occt_core::gp::{GpTrsf2d, GpVec2d};
use occt_geom2d::curve::Curve2d;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Face, TopoShape};
use crate::topo_tools::explorer_shapes;
use crate::topopebrep::ds::build_tool::BuildTool;

/// `FUN_tool_sortVonE(lov, E)` (`TOOL.cxx:436-469`): sort the vertices of
/// `lov` by their parameter on `E`.
///
/// OCCT keys an `IndexedMap<double>` by parameter (`mappar.Add(par)`) and a
/// `DataMap<int, shape>` by that index (`mapiv.Bind(iv, v)`), then sorts a
/// copy of the parameter array with `std::sort` and re-emits the shapes.
/// `IndexedMap::Add` returns the existing index for a parameter already in the
/// map, so a duplicate parameter keeps its first slot and the later shape
/// overwrites it — the port reproduces that. `std::sort` is not stable, so
/// OCCT leaves the order of equal parameters unspecified; the port sorts
/// stably, which is one of the orders `std::sort` may produce.
pub fn fun_tool_sort_von_e(lov: &mut Vec<TopoShape>, e: &TopoShape) {
    // cxx:439-446
    let mut order: Vec<(f64, usize)> = Vec::new();
    let mut by_index: Vec<(usize, TopoShape)> = Vec::new();
    for v in lov.iter() {
        let par = BRepTool::parameter_of_vertex(e, v);
        let iv = match order.iter().find(|(p, _)| *p == par) {
            Some((_, i)) => *i,
            None => {
                let i = by_index.len() + 1;
                order.push((par, i));
                by_index.push((i, v.clone()));
                i
            }
        };
        if let Some(slot) = by_index.iter_mut().find(|(i, _)| *i == iv) {
            slot.1 = v.clone();
        }
    }

    // cxx:448-458
    let mut pars: Vec<f64> = order.iter().map(|(p, _)| *p).collect();
    pars.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    // cxx:460-466
    let mut newlov: Vec<TopoShape> = Vec::with_capacity(pars.len());
    for par in pars {
        let Some((_, iv)) = order.iter().find(|(p, _)| *p == par) else {
            continue;
        };
        if let Some((_, v)) = by_index.iter().find(|(i, _)| i == iv) {
            newlov.push(v.clone());
        }
    }

    // cxx:467-468
    lov.clear();
    lov.extend(newlov);
}

/// `TopOpeBRepTool_TOOL::SplitE(Eanc, Splits)` (`TOOL.cxx:473-528`).
///
/// Returns the OCCT `bool`: `false` when `Eanc` has two vertices or fewer, or
/// when its vertex list is empty (`Splits` is left untouched then).
pub fn split_e(eanc: &TopoShape, splits: &mut Vec<TopoShape>) -> bool {
    // cxx:476-479
    let o_eanc = eanc.orientation();
    let efor = eanc.oriented(Orientation::Forward);

    // cxx:480-488
    let mut lov: Vec<TopoShape> = explorer_shapes(&efor, ShapeType::Vertex);

    // cxx:489-492
    if lov.len() <= 2 {
        return false;
    }

    // cxx:493
    fun_tool_sort_von_e(&mut lov, &efor);

    // cxx:495-505
    let mut it = lov.into_iter();
    let Some(mut v0) = it.next() else {
        return false;
    };

    // cxx:507-527
    for mut v in it {
        let par0 = BRepTool::parameter_of_vertex(&efor, &v0);
        let par = BRepTool::parameter_of_vertex(&efor, &v);

        // cxx:516-517: `TopoDS_Edge ed; FUN_ds_CopyEdge(EFOR, ed)`.
        let _ed = BuildTool::new().copy_edge(&efor);

        // cxx:518-521: `v0.Orientation(FORWARD); BB.Add(ed, v0);
        // FUN_ds_Parameter(ed, v0, par0)`.
        let mut v0f = v0.clone();
        v0f.set_orientation(Orientation::Forward);
        let mut ed = BuildTool::new().copy_edge(&efor);
        TopoBuilder::new().add(&mut ed, &v0f);
        BuildTool::new().parameter(&ed, &v0f, par0);

        // cxx:522-524: `v.Orientation(REVERSED); BB.Add(ed, v);
        // FUN_ds_Parameter(ed, v, par)`.
        v.set_orientation(Orientation::Reversed);
        TopoBuilder::new().add(&mut ed, &v);
        BuildTool::new().parameter(&ed, &v, par);

        // cxx:525: `Splits.Append(ed.Oriented(oEanc))`.
        splits.push(ed.oriented(o_eanc));

        // cxx:526: `v0 = v` — the next interval starts at the current vertex,
        // carrying the orientation given to it above.
        v0 = v;
    }
    true
}

/// `TopOpeBRepTool_TOOL::TrslUVModifE(t2d, F, E)` (`TOOL.cxx:1319-1335`).
///
/// Returns `false` when `E` has no p-curve on `F` (the OCCT
/// `if (PC.IsNull()) return false`); otherwise the translated p-curve is
/// written back through `BRep_Builder::UpdateEdge(E, PC, F, tol)`
/// (`cxx:1331-1333`), which `BuildTool::PCurve` reproduces
/// (`TopOpeBRepDS_BuildTool.cxx:1188-1239`).
///
/// Deviation: `FC2D_CurveOnSurface` reports `tol = BRep_Tool::Tolerance(E)`
/// when the p-curve is stored (`TopOpeBRepTool_2d.cxx:162-168`), and
/// `BRep_Builder::UpdateEdge` only raises the edge tolerance to `tol`
/// (`BRep_Builder.cxx:1205-1240`); `BuildTool::PCurve` writes the p-curve
/// without touching the tolerance, so the tolerance side of the OCCT call has
/// no effect here — the value can only ever equal the edge's own tolerance,
/// which the edge already has.
pub fn trsl_uv_modif_e(t2d: &GpVec2d, f: &TopoShape, e: &TopoShape) -> bool {
    let Some(ff) = Face::wrap(f.clone()) else {
        return false;
    };
    let Some(ee) = crate::shape::Edge::wrap(e.clone()) else {
        return false;
    };
    // `cxx:1321-1327`: `FC2D_CurveOnSurface(E, F, f, l, tol)`.
    let Some((pc, _first, _last)) = crate::boptools_2d::curve_on_surface_range(&ee, &ff) else {
        return false;
    };
    // `cxx:1329`: `PC->Translate(t2d)`.
    let mut tr = GpTrsf2d::identity();
    tr.set_translation_vec(t2d);
    let moved: Arc<dyn Curve2d> = Arc::from(pc.transformed(&tr));
    // `cxx:1331-1333`
    BuildTool::new().p_curve(f, e, &moved);
    true
}
