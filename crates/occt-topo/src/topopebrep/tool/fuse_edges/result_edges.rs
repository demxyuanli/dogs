//! `BuildListResultEdges` (`TopOpeBRepTool_FuseEdges.cxx:272-380`).
//!
//! For every chain built by `BuildListEdges`, the first edge is replaced by a
//! single edge built over the chain's support between the two **outer**
//! vertices:
//!
//! - the outer vertices are taken with the chain's traversal in mind
//!   (`cxx:305-315`): for a `REVERSED` first edge `VL` is the first vertex of
//!   the first edge and `VF` the last vertex of the last one, otherwise the
//!   other way round;
//! - the support curve is `BRep_Tool::Curve(OldEdge)` moved into world space by
//!   the edge location (`cxx:316-321`) and unwrapped one `Geom_TrimmedCurve`
//!   level (`cxx:324-327`) — "otherwise the makeedge will fail";
//! - `BRepLib_MakeEdge ME(C, VF, VL)` (`cxx:334`) is the `Init` overload with
//!   the two vertices, i.e. [`crate::brep_lib_make_edge::init_3d_projecting`];
//! - the p-curves are then rebuilt on the new edge by
//!   [`super::update_pcurve::update_pcurve`] (`cxx:373-377`); a chain whose
//!   p-curves cannot be merged is simply not bound, which leaves that chain
//!   unfused.
//!
//! UNPORTED (`cxx:334-363`): when `ME.IsDone()` is false OCCT copies the curve,
//! calls `GeomLib::ExtendCurveToPoint` twice (once per end) and retries
//! `ME.Init(ExtC, VF, VL)`, throwing `Standard_ConstructionError` when even
//! that fails. `GeomLib::ExtendCurveToPoint` (`GeomLib.cxx:1268-1411`) is not
//! ported (it needs `PLib::HermiteCoefficients` and `GeomLib::ComputeLambda`),
//! so the failure arm answers [`FuseEdgesError::FusionFailed`] directly.

use std::sync::Arc;

use occt_geom::Curve;

use crate::abs::Orientation;
use crate::brep_lib_make_edge::init_3d_projecting;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::Edge;

use super::{cumulated_vertex, FuseEdges, FuseEdgesError};

impl FuseEdges {
    /// `BuildListResultEdges()` (`cxx:272-380`).
    pub(super) fn build_list_result_edges(&mut self) -> Result<(), FuseEdgesError> {
        // `cxx:284-289`: nothing to do without chains.
        if self.my_map_lst_edg.is_empty() {
            return Ok(());
        }

        self.my_map_edg.clear();
        let keys: Vec<i32> = self.my_map_lst_edg.keys().copied().collect();

        for i_lst in keys {
            let Some(lmap_edg) = self.my_map_lst_edg.get(&i_lst).cloned() else {
                continue;
            };
            if lmap_edg.is_empty() {
                continue;
            }
            // `cxx:303`: `const TopoDS_Edge& OldEdge = TopoDS::Edge(LmapEdg.First());`
            let Some(old_edge) = Edge::wrap(lmap_edg[0].clone()) else {
                continue;
            };

            // `cxx:305-315`: the outer vertices, oriented along the chain.
            let last_shape = lmap_edg.last().cloned().unwrap();
            let (vf, vl) = if old_edge.orientation() == Orientation::Reversed {
                // `VL = FirstVertex(LmapEdg.First(), true)`,
                // `VF = LastVertex(LmapEdg.Last(), true)`.
                (
                    cumulated_vertex(&last_shape, true),
                    cumulated_vertex(&lmap_edg[0], false),
                )
            } else {
                // `VF = FirstVertex(LmapEdg.First(), true)`,
                // `VL = LastVertex(LmapEdg.Last(), true)`.
                (
                    cumulated_vertex(&lmap_edg[0], false),
                    cumulated_vertex(&last_shape, true),
                )
            };
            // OCCT dereferences the null vertex in `BRep_Tool::Pnt` and raises
            // `Standard_NullObject`; the port keeps the chain unfused.
            let (Some(vf), Some(vl)) = (vf, vl) else {
                continue;
            };

            // `cxx:316-321`: `C = BRep_Tool::Curve(OldEdge, loc, f, l)`,
            // moved into world space when the edge carries a location.
            let Some(curve) = BRepTool::edge_curve(&old_edge) else {
                continue;
            };
            let loc = old_edge.location().clone();
            // `cxx:324-327`: "if the curve is trimmed we get the basis curve to
            // fit the new vertices otherwise the makeedge will fail".
            let curve: Arc<dyn Curve> = {
                let world: Arc<dyn Curve> = if loc.is_identity() {
                    curve
                } else {
                    Arc::from(curve.transformed(&loc.transformation()))
                };
                match world.untrimmed_basis() {
                    Some((basis, _, _)) => basis,
                    None => world,
                }
            };

            // `cxx:334`: `BRepLib_MakeEdge ME(C, VF, VL)` — `Init(C, V1, V2)`
            // (`BRepLib_MakeEdge.cxx:542-572`) followed by
            // `Init(C, V1, V2, p1, p2)` (`cxx:606-796`).
            let Some(build) =
                init_3d_projecting(&curve, Some(&crate::shape::Vertex(vf.clone())), Some(&crate::shape::Vertex(vl.clone())))
            else {
                // UNPORTED `cxx:334-363` (`GeomLib::ExtendCurveToPoint`).
                return Err(FuseEdgesError::FusionFailed);
            };

            // `cxx:366`: `NewEdge = ME.Edge()` — `BRepLib_MakeEdge.cxx:780-795`:
            // `MakeEdge(E, C, preci)`, `Add(E, V1)`, `Add(E, V2)`, `Range`,
            // `Degenerated`, `Closed`.
            let b = TopoBuilder::new();
            let mut new_edge = b.make_edge(curve, build.first, build.last);
            crate::brep_lib_make_edge::add_endpoint_vertices(&b, &mut new_edge, &build);
            if build.closed {
                new_edge.0.set_closed(true);
            }
            if build.degenerated {
                crate::tgeometry::GeometryRegistry::global().set_degenerated(&new_edge.0, true);
            }

            // `cxx:373-377`
            if self.update_pcurve(&old_edge, &mut new_edge, &lmap_edg) {
                self.my_map_edg.insert(i_lst, new_edge.0);
            }
        }

        self.my_result_edges_done = true;
        Ok(())
    }
}
