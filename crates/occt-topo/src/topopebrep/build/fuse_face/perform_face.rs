//! `TopOpeBRepBuild_FuseFace::PerformFace` (`cxx:126-582`) — fuse the
//! co-surfacic faces of `myLRF` that are connected across a shared edge.
//!
//! Level 1 groups the faces (`cxx:214-245`), level 2 groups the wires of one
//! face family (`cxx:256-284`), level 3 groups the edges of one wire family
//! (`cxx:286-490`). A family of more than one shape is then rebuilt: the edges
//! are classified by their occurrences' orientations (`cxx:328-426`), the
//! kept edges are chained back into wires (`cxx:428-487`), and the wires —
//! plus one single-edge wire per internal / external / modified edge
//! (`cxx:505-535`) — are added to a face built on the family's surface
//! (`cxx:492-503`).

use occt_core::precision::CONFUSION;
use occt_geom::Surface;

use crate::abs::{Orientation, ShapeType};
use crate::brep_check::Analyzer;
use crate::brep_lib_make_wire::MakeWire;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::topo_tools::{explorer_shapes, shape_map_key};

use super::group_shape::group_shape;
use super::shape_list_map::ShapeListMap;
use super::FuseFace;

/// `BRepLib_MakeWire(E)` (`BRepLib_MakeWire.cxx`, the edge constructor):
/// a one-edge wire, null — i.e. not added at all (`BRep_Builder::Add` skips a
/// null shape, `BRep_Builder.cxx`) — when the wire could not be built.
fn single_edge_wire(shape: &TopoShape) -> Option<Wire> {
    let edge = Edge::wrap(shape.clone())?;
    let mut mw = MakeWire::new();
    mw.add(&edge).ok()?;
    if !mw.is_done() {
        return None;
    }
    Some(mw.wire())
}

/// `BRepLib_MakeWire MW; MW.Add(L)` followed by the "astuce pour contourner
/// Wire Not Closed", `BRepLib_MakeWire MW1(W); W = MW1.Wire()`
/// (`cxx:434-452`, `:739-757`). `None` is `!MW.IsDone()`.
pub(super) fn make_wire_from_list(edges: &[TopoShape]) -> Option<Wire> {
    let mut mw = MakeWire::new();
    for e in edges {
        let edge = Edge::wrap(e.clone())?;
        // `BRepLib_MakeWire::Add(const NCollection_List<TopoDS_Shape>&)` adds
        // every edge in list order; the error state is sticky.
        let _ = mw.add(&edge);
    }
    if !mw.is_done() {
        return None;
    }
    // Rebuild once more, edge by edge, like `MW1(W)`.
    let mut mw1 = MakeWire::new();
    for e in explorer_shapes(&mw.wire().0, ShapeType::Edge) {
        if let Some(edge) = Edge::wrap(e) {
            let _ = mw1.add(&edge);
        }
    }
    Some(mw1.wire())
}

impl FuseFace {
    /// `PerformFace()` (`cxx:126-582`).
    pub(crate) fn perform_face(&mut self) {
        self.my_modified = false;
        self.my_lff.clear();

        // `cxx:136-147`.
        if self.my_lrf.is_empty() {
            self.no_fuse();
            return;
        }
        // `cxx:149-161`.
        if self.my_lrf.len() == 1 {
            self.no_fuse();
            return;
        }

        let lrf = self.my_lrf.clone();
        let (mut ori3d_forward, mut ori3d_reversed) = (false, false);
        let mut mylist: Vec<TopoShape> = Vec::with_capacity(lrf.len());
        // `cxx:172-201`.
        for fac in &lrf {
            match fac.orientation() {
                Orientation::Forward => ori3d_forward = true,
                Orientation::Reversed => ori3d_reversed = true,
                _ => {}
            }
            // `cxx:184-197`: an invalid reconstructed face aborts the fusion.
            if !Analyzer::new(fac.clone()).is_valid() {
                self.no_fuse();
                return;
            }
            // `cxx:198`: the family is grouped with every face FORWARD.
            mylist.push(fac.oriented(Orientation::Forward));
        }

        // `cxx:203-212`: mixed 3d orientations cannot be fused.
        if ori3d_forward && ori3d_reversed {
            self.no_fuse();
            return;
        }

        // `cxx:214-245`: the face families, `Keep_Edge = false`.
        let mut map_fac_l_fac = ShapeListMap::new();
        group_shape(&mylist, false, &mut map_fac_l_fac);
        if map_fac_l_fac.is_empty() {
            self.no_fuse();
            return;
        }
        if lrf.len() == map_fac_l_fac.extent() {
            self.no_fuse();
            return;
        }

        // `cxx:247-560`: rebuild one face per family.
        for (fac, lfac) in map_fac_l_fac.pairs() {
            if lfac.len() == 1 {
                // `cxx:554-558`: a family of one keeps its representative.
                self.my_lff.push(fac.clone());
                continue;
            }
            if !self.fuse_face_family(fac, lfac, ori3d_reversed) {
                self.no_fuse();
                return;
            }
        }

        // `cxx:562-576`.
        if self.my_lff.is_empty() {
            self.no_fuse();
            return;
        }

        // `cxx:578-579`.
        self.my_modified = true;
        self.my_done = true;
    }

    /// The `n11 != 1` arm of `PerformFace` (`cxx:256-553`): rebuild one face
    /// from the wires of a face family. `false` is an OCCT early `return`.
    fn fuse_face_family(
        &mut self,
        fac: &TopoShape,
        lfac: &[TopoShape],
        ori3d_reversed: bool,
    ) -> bool {
        // `cxx:258-268`: every wire of the family.
        let mut lwir: Vec<TopoShape> = Vec::new();
        for fac1 in lfac {
            lwir.extend(explorer_shapes(fac1, ShapeType::Wire));
        }

        // `cxx:270-284`: the wire families, `Keep_Edge = false`.
        let mut map_wir_l_wir = ShapeListMap::new();
        group_shape(&lwir, false, &mut map_wir_l_wir);
        if map_wir_l_wir.is_empty() {
            return false;
        }

        // `cxx:286-288`: per-family output lists.
        let mut my_face_lie: Vec<TopoShape> = Vec::new();
        let mut my_face_lee: Vec<TopoShape> = Vec::new();
        let mut my_face_lme: Vec<TopoShape> = Vec::new();
        let mut my_face_lw: Vec<Wire> = Vec::new();

        // `cxx:289-489`: one wire family at a time.
        for (wir, lwir1) in map_wir_l_wir.pairs() {
            if lwir1.len() == 1 {
                // `cxx:488`: a wire that shares no edge is kept as it is.
                if let Some(w) = Wire::wrap(wir.clone()) {
                    my_face_lw.push(w);
                }
                continue;
            }

            // `cxx:295-313`: the edges of the family, then their families with
            // `Keep_Edge = true` (the shared-edge aggregation).
            let mut ledg: Vec<TopoShape> = Vec::new();
            for wir1 in lwir1 {
                ledg.extend(explorer_shapes(wir1, ShapeType::Edge));
            }
            let mut map_edg_l_edg = ShapeListMap::new();
            group_shape(&ledg, true, &mut map_edg_l_edg);
            if map_edg_l_edg.is_empty() {
                return false;
            }

            // `cxx:329-426`: classify every shared edge by its occurrences.
            let mut my_wire_le: Vec<TopoShape> = Vec::new();
            for (edg, ledg1) in map_edg_l_edg.pairs() {
                let (mut ori_reversed, mut ori_forward) = (false, false);
                let (mut ori_internal, mut ori_external) = (false, false);
                for edg1 in ledg1 {
                    match edg1.orientation() {
                        Orientation::Reversed => {
                            // `cxx:341-356`: two REVERSED occurrences mean the
                            // faces are not consistently oriented.
                            if ori_reversed {
                                return false;
                            }
                            ori_reversed = true;
                        }
                        Orientation::Forward => {
                            if ori_forward {
                                return false;
                            }
                            ori_forward = true;
                        }
                        Orientation::Internal => ori_internal = true,
                        Orientation::External => ori_external = true,
                    }
                }

                // `cxx:388-425`: forward / reversed first, then internal,
                // then external.
                if ori_reversed || ori_forward {
                    if ori_reversed && ori_forward {
                        let e = edg.oriented(Orientation::Internal);
                        self.my_lme.push(e.clone());
                        my_face_lme.push(e);
                    } else if ori_reversed {
                        my_wire_le.push(edg.oriented(Orientation::Reversed));
                    } else {
                        my_wire_le.push(edg.oriented(Orientation::Forward));
                    }
                } else if ori_internal {
                    let e = edg.oriented(Orientation::Internal);
                    self.my_lie.push(e.clone());
                    my_face_lie.push(e);
                } else if ori_external {
                    let e = edg.oriented(Orientation::External);
                    self.my_lee.push(e.clone());
                    my_face_lee.push(e);
                }
            }

            // `cxx:428-487`: rebuild as many wires as the edge list holds, one
            // `BRepLib_MakeWire` pass per wire.
            let mut number1 = my_wire_le.len();
            while number1 > 0 {
                let Some(w) = make_wire_from_list(&my_wire_le) else {
                    return false;
                };
                my_face_lw.push(w.clone());

                // `cxx:456-469`: mark the edges of the new wire.
                let mut m: Vec<_> = Vec::new();
                let mut nb = 0usize;
                for edg3 in explorer_shapes(&w.0, ShapeType::Edge) {
                    let k = shape_map_key(&edg3);
                    if !m.contains(&k) {
                        m.push(k);
                    }
                    nb += 1;
                }

                if nb == number1 {
                    number1 = 0;
                } else {
                    // `cxx:470-484`: keep the edges the wire did not consume.
                    let mut list_edge: Vec<TopoShape> = Vec::new();
                    for edg2 in &my_wire_le {
                        let k = shape_map_key(edg2);
                        if !m.contains(&k) {
                            m.push(k);
                            list_edge.push(edg2.clone());
                        }
                    }
                    my_wire_le = list_edge;
                    number1 = my_wire_le.len();
                }
            }
        }

        // `cxx:492-498`: the family's surface, untrimmed when it is a
        // `Geom_RectangularTrimmedSurface`.
        let Some(facref) = Face::wrap(fac.clone()) else {
            return false;
        };
        let Some(surface) = BRepTool::face_surface(&facref) else {
            return false;
        };
        let surface: std::sync::Arc<dyn Surface> =
            surface.rectangular_trimmed_basis().unwrap_or(surface);

        // `cxx:499-506`: the outer and inner wires.
        let mut wires: Vec<Wire> = my_face_lw;
        // `cxx:508-524`: internal and external edges as one-edge wires.
        for edg1 in &my_face_lie {
            if let Some(w) = single_edge_wire(edg1) {
                wires.push(w);
            }
        }
        for edg1 in &my_face_lee {
            if let Some(w) = single_edge_wire(edg1) {
                wires.push(w);
            }
        }
        // `cxx:525-535`: the modified edges only when internal connections are
        // kept (`Init(..., CXM == 2)`).
        if self.my_internal {
            for edg1 in &my_face_lme {
                if let Some(w) = single_edge_wire(edg1) {
                    wires.push(w);
                }
            }
        }

        // `cxx:499`: `BRepLib_MakeFace MF(S, Precision::Confusion())`, then one
        // `MF.Add(wir)` per wire (`cxx:500-535`). `MF.IsDone()` (`cxx:537-548`)
        // can only be false for a null surface, which the `face_surface` guard
        // above already answered.
        let b = TopoBuilder::new();
        let mut f = b.make_face_with_tolerance(surface, &wires, CONFUSION);

        // `cxx:550-553`.
        if ori3d_reversed {
            f.0.reverse();
        }
        self.my_lff.push(f.0);
        true
    }
}
