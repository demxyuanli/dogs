//! `TopOpeBRepBuild_Builder` — the edge builder driver of
//! `TopOpeBRepBuild_makeedges.cxx`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_makeedges.cxx:79-228`
//! (OCCT 8.0.0).
//!
//! `GEDBUMakeEdges` walks the edge/vertex iteration of a
//! `TopOpeBRepBuild_EdgeBuilder` and, for every edge with more than one vertex
//! loop, copies `EF` and adds the loop vertices, then splits the result at its
//! internal vertices (`TopOpeBRepTool_TOOL::SplitE`) and appends it to `LOE`.
//!
//! `myDataStructure->HasSameDomain(V)` is asked with `FindKeep = false` in
//! OCCT (`TopOpeBRepDS_HDataStructure::HasSameDomain(const TopoDS_Shape&)`,
//! default `FindKeep = Standard_True`); the port passes `true`.

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, is_equal};
use crate::topopebrep::tool::split_e::split_e;

use super::super::edge_builder::EdgeBuilder;
use super::Builder;

impl Builder {
    /// `GEDBUMakeEdges(EF, EDBU, LOE)` (`makeedges.cxx:79-228`).
    pub fn gedbu_make_edges(
        &self,
        ef: &TopoShape,
        edbu: &mut EdgeBuilder,
        loe: &mut Vec<TopoShape>,
    ) {
        edbu.init_edge();
        while edbu.more_edge() {
            // `cxx:108-115`: count the vertex loops of the current edge.
            let mut nloop = 0;
            edbu.init_vertex();
            while edbu.more_vertex() {
                nloop += 1;
                edbu.next_vertex();
            }
            // `cxx:116-118`: 0 or 1 vertex on the current edge => drop it.
            if nloop <= 1 {
                edbu.next_edge();
                continue;
            }

            let mut new_edge = self.my_build_tool.copy_edge(ef);

            let mut n_vf = 0;
            let mut n_vr = 0;
            let mut tosplit = false;

            edbu.init_vertex();
            while edbu.more_vertex() {
                let Some(v_raw) = edbu.vertex() else {
                    edbu.next_vertex();
                    continue;
                };
                let vori = v_raw.orientation();

                // `cxx:131-137`: with a same-domain vertex, take its reference.
                let mut v = v_raw.clone();
                let hassd = self
                    .data_structure()
                    .map(|h| h.has_same_domain(&v, true))
                    .unwrap_or(false);
                if hassd {
                    if let Some(hds) = self.data_structure() {
                        let iref = hds.same_domain_reference(&v);
                        if let Some(vref) = hds.shape(iref, true) {
                            v = vref;
                            v.set_orientation(vori);
                        }
                    }
                }

                let ori_v = v.orientation();
                if ori_v == Orientation::External {
                    edbu.next_vertex();
                    continue;
                }

                // `cxx:144-174`
                let par_v = edbu.parameter();
                let mut equafound = false;
                for ve in explorer_shapes(&new_edge, ShapeType::Vertex) {
                    let ori_ve = ve.orientation();
                    if is_equal(&v, &ve) {
                        equafound = true;
                        break;
                    } else if ori_ve == Orientation::Forward || ori_ve == Orientation::Reversed {
                        if ori_v == ori_ve {
                            equafound = true;
                            break;
                        }
                    } else if ori_ve == Orientation::Internal
                        || ori_ve == Orientation::External
                    {
                        let par_ve = BRepTool::parameter_of_vertex(&new_edge, &ve);
                        if par_v == par_ve {
                            equafound = true;
                            break;
                        }
                    }
                }

                if !equafound {
                    // `cxx:179-192`: `VF` / `VR` only feed the "Closed" bit of
                    // OCCT, which no later statement reads; the port keeps the
                    // counters only.
                    if vori == Orientation::Forward {
                        n_vf += 1;
                    }
                    if vori == Orientation::Reversed {
                        n_vr += 1;
                    }
                    if ori_v == Orientation::Internal {
                        tosplit = true;
                    }
                    self.my_build_tool.add_edge_vertex(&mut new_edge, &v);
                    self.my_build_tool.parameter(&new_edge, &v, par_v);
                }

                edbu.next_vertex();
            }

            let addedge = n_vf == 1 && n_vr == 1;
            if addedge {
                if tosplit {
                    let mut loe_sub: Vec<TopoShape> = Vec::new();
                    let ok = split_e(&new_edge, &mut loe_sub);
                    if !ok {
                        tosplit = false;
                    } else {
                        loe.extend(loe_sub);
                    }
                }
                if !tosplit {
                    loe.push(new_edge);
                }
            }

            edbu.next_edge();
        }
    }
}
