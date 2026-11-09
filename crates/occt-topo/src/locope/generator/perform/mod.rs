//! `LocOpe_Generator::Perform` (`LocOpe_Generator.cxx:71-439`): the seed faces,
//! the face-fusion search, the connecting edges and the per-vertex fusion
//! decisions.
//!
//! Stages 4 to 6 live in [`super::edge_fusion`], [`super::face_fusion`] and
//! [`super::assemble`].

mod assemble;
mod edge_fusion;
mod face_fusion;

use crate::abs::ShapeType;
use crate::locope::generated_shape::GeneratedShape;
use crate::locope::shape_maps::{ShapeDataMap, ShapeSet, ShapeValueMap};
use crate::shape::TopoShape;
use crate::topo_tools::{is_same, AncestorMap};

/// The locals of `Perform` (`cxx:90-147`) plus the base shape and the
/// `GeneratedShape` query interface. Every stage mutates them in place, exactly
/// as the reference's locals live across the whole body.
pub(super) struct Stage<'a> {
    /// `myShape` (`cxx:71`).
    pub(super) shape: TopoShape,
    /// The `G` argument (`cxx:71`).
    pub(super) g: &'a dyn GeneratedShape,
    /// `theLeft` (`cxx:91`): the faces on the left of the generating wire.
    pub(super) left: ShapeSet,
    /// `GEdg` (`cxx:94`): the generating edges.
    pub(super) g_edg: ShapeSet,
    /// `GVtx` (`cxx:94`): the endpoints of the generating edges.
    pub(super) g_vtx: ShapeSet,
    /// `theEFMap` (`cxx:133-135`).
    pub(super) ef_map: AncestorMap,
    /// `theEEMap` (`cxx:137-138`): edge -> the twin edges it fuses with.
    pub(super) ee_map: ShapeValueMap<Vec<TopoShape>>,
    /// `theFFMap` (`cxx:139-140`): base face -> the faces fused into it.
    pub(super) ff_map: ShapeValueMap<Vec<TopoShape>>,
    /// `toRemove` (`cxx:142`).
    pub(super) to_remove: ShapeSet,
    /// `DontFuse` (`cxx:269`): edge -> the face it must *not* fuse with.
    pub(super) dont_fuse: ShapeDataMap,
    /// `RebuildFace` (`cxx:268`).
    pub(super) rebuild_face: Vec<TopoShape>,
    /// `mapTreated` (`cxx:269`).
    pub(super) map_treated: ShapeSet,
}

impl<'a> Stage<'a> {
    /// Builds the stage and runs `theEFMap` construction, which the reference
    /// does between the seed loop and the fusion search (`cxx:133-135`).
    pub(super) fn new(shape: TopoShape, g: &'a dyn GeneratedShape) -> Self {
        let ef_map = AncestorMap::map_shapes_and_ancestors(&shape, ShapeType::Edge, ShapeType::Face);
        Self {
            shape,
            g,
            left: ShapeSet::new(),
            g_edg: ShapeSet::new(),
            g_vtx: ShapeSet::new(),
            ef_map,
            ee_map: ShapeValueMap::new(),
            ff_map: ShapeValueMap::new(),
            to_remove: ShapeSet::new(),
            dont_fuse: ShapeDataMap::new(),
            rebuild_face: Vec::new(),
            map_treated: ShapeSet::new(),
        }
    }

    /// `cxx:90-1090`: the whole `Perform` body up to the assembly.
    pub(super) fn run(&mut self, mod_shapes: &mut ShapeValueMap<Vec<TopoShape>>) {
        // `cxx:90-131`.
        self.seed_left(mod_shapes);
        // `cxx:144-265`.
        self.search_face_fusions();
        // `cxx:231-265`.
        self.add_connectivity_edges();
        // `cxx:267-439`.
        self.classify_vertices();
        // `cxx:454-522`.
        edge_fusion::fuse_edges(self, mod_shapes);
        // `cxx:526-1090`.
        face_fusion::fuse_faces(self, mod_shapes);
    }

    /// `cxx:1133-1218`.
    pub(super) fn finish(
        &self,
        mod_shapes: &mut ShapeValueMap<Vec<TopoShape>>,
    ) -> Option<TopoShape> {
        assemble::build(self, mod_shapes)
    }

    /// `cxx:90-131`: the faces on the left of the generating edges become seeds,
    /// and every one of them gets an empty descendant list.
    fn seed_left(&mut self, mod_shapes: &mut ShapeValueMap<Vec<TopoShape>>) {
        let ledges = self.g.generating_edges().to_vec();
        for edg in &ledges {
            // `cxx:100-108`.
            self.g_edg.add(edg);
            for vtx in crate::topo_tools::explorer_shapes(edg, ShapeType::Vertex) {
                self.g_vtx.add(&vtx);
            }
            // `cxx:110-130`: the first face of `myShape` that carries `edg` in
            // the same orientation, then stop - the reference `break`s out of
            // the face loop once one is added.
            for fac in crate::topo_tools::explorer_shapes(&self.shape, ShapeType::Face) {
                let matched = crate::topo_tools::explorer_shapes(&fac, ShapeType::Edge)
                    .iter()
                    .any(|e| is_same(e, edg) && e.orientation() == edg.orientation());
                if !matched {
                    continue;
                }
                self.left.add(&fac);
                if !mod_shapes.is_bound(&fac) {
                    mod_shapes.bind(&fac, Vec::new());
                }
                break;
            }
        }
    }

    /// `cxx:144-224`: for every generating edge whose first ancestor is not on
    /// the left, fuse the generated face into it when both are the same plane.
    fn search_face_fusions(&mut self) {
        for edg in self.g_edg.iter().to_vec() {
            let Some(ancestors) = self.ef_map.find_from_key(&edg) else {
                // `cxx:148-151`.
                continue;
            };
            // `cxx:152-163`: the first ancestor outside `theLeft`. None means
            // the edge is internal to the shell or on a free border.
            let mut target = None;
            for anc in ancestors {
                if !self.left.contains(anc) {
                    target = Some(anc.clone());
                    break;
                }
            }
            let Some(fac) = target else {
                continue;
            };
            // `cxx:165-168`.
            let Some(facbis) = self.g.generated_face(&edg) else {
                // `cxx:218-221`: a face generated by `edg` is only marked, and
                // the reference's `myFFromE` binding is commented out.
                continue;
            };
            if !super::statics::to_fuse_faces(&fac, &facbis) {
                continue;
            }
            // `cxx:169-193`: a face already fused with `facbis` takes `fac`
            // into its list.
            let mut found = None;
            for (key, value) in self.ff_map.iter() {
                if is_same(key, &fac) {
                    continue;
                }
                if value.iter().any(|v| is_same(v, &facbis)) {
                    found = Some(key.clone());
                    break;
                }
            }
            if let Some(key) = found {
                if let Some(list) = self.ff_map.change(&key) {
                    list.push(fac.clone());
                }
                self.to_remove.add(&fac);
                self.to_remove.add(&edg);
                continue;
            }
            // `cxx:196-216`.
            if !self.ff_map.is_bound(&fac) {
                self.ff_map.bind(&fac, Vec::new());
            }
            let already = self
                .ff_map
                .find(&fac)
                .is_some_and(|l| l.iter().any(|v| is_same(v, &facbis)));
            if !already {
                if let Some(list) = self.ff_map.change(&fac) {
                    list.push(facbis.clone());
                }
            }
            self.to_remove.add(&edg);
            self.to_remove.add(&facbis);
        }
    }

    /// `cxx:231-265`: an edge shared by two faces of the same fusion group is a
    /// connection edge and goes.
    fn add_connectivity_edges(&mut self) {
        let entries: Vec<(TopoShape, Vec<TopoShape>)> =
            self.ff_map.iter().map(|(k, v)| (k.clone(), v.clone())).collect();
        for (_, value) in entries {
            for itl in &value {
                for ed in crate::topo_tools::explorer_shapes(itl, ShapeType::Edge) {
                    // `cxx:241-244`.
                    if self.to_remove.contains(&ed) {
                        continue;
                    }
                    // `cxx:245-261`.
                    let mut shared = false;
                    for itl2 in &value {
                        if is_same(itl2, itl) {
                            continue;
                        }
                        if crate::topo_tools::explorer_shapes(itl2, ShapeType::Edge)
                            .iter()
                            .any(|e| is_same(e, &ed))
                        {
                            shared = true;
                            break;
                        }
                    }
                    if shared {
                        self.to_remove.add(&ed);
                    }
                }
            }
        }
    }

    /// `cxx:267-439`: per edge of every fusion group, decide for each generating
    /// vertex whether the edge fuses with the twin edge the vertex generated.
    fn classify_vertices(&mut self) {
        let entries: Vec<(TopoShape, Vec<TopoShape>)> =
            self.ff_map.iter().map(|(k, v)| (k.clone(), v.clone())).collect();
        for (fac, _) in entries {
            for edg in crate::topo_tools::explorer_shapes(&fac, ShapeType::Edge) {
                // `cxx:281-285`.
                if self.map_treated.contains(&edg) {
                    continue;
                }
                self.map_treated.add(&edg);
                // `cxx:287-292`.
                let vertices = crate::topo_tools::explorer_shapes(&edg, ShapeType::Vertex);
                for vtx in vertices {
                    if !self.g_vtx.contains(&vtx) {
                        continue;
                    }
                    self.classify_one(&fac, &edg, &vtx);
                }
            }
        }
        // `cxx:435-439`: a face to rebuild gets an empty fusion list.
        for fac in self.rebuild_face.clone() {
            self.ff_map.bind(&fac, Vec::new());
        }
    }

    /// The body of the `cxx:287-428` vertex loop for one `(fac, edg, vtx)`.
    fn classify_one(&mut self, fac: &TopoShape, edg: &TopoShape, vtx: &TopoShape) {
        // `cxx:294-299`.
        let Some(edgbis) = self.g.generated_edge(vtx) else {
            return;
        };
        let Some(edgbis_edge) = crate::shape::Edge::wrap(edgbis.clone()) else {
            return;
        };
        if crate::brep_tool::BRepTool::is_degenerated(&edgbis_edge)
            || !super::statics::to_fuse_edges(edg, &edgbis)
        {
            return;
        }
        // `cxx:302`: `BRepTools::IsReallyClosed(edg, fac)`.
        let really_closed = match crate::shape::Face::wrap(fac.clone()) {
            Some(f) => crate::brep_tools::is_really_closed(&edge_of(edg), &f),
            None => return,
        };
        if really_closed {
            self.classify_really_closed(fac, edg, vtx, &edgbis);
        } else {
            self.classify_open(fac, edg, vtx, &edgbis);
        }
    }

    /// `cxx:302-357`: the edge closes on the face, so the twin chain is built
    /// once for `edg` and each vertex either fuses or is bound in `DontFuse`.
    fn classify_really_closed(
        &mut self,
        fac: &TopoShape,
        edg: &TopoShape,
        vtx: &TopoShape,
        edgbis: &TopoShape,
    ) {
        if self.ee_map.is_bound(edg) {
            return;
        }
        self.ee_map.bind(edg, vec![edgbis.clone()]);
        // `cxx:308`: always true for a double edge.
        self.to_remove.add(edgbis);
        let mut fuse_edge = true;
        let (_, vl) = super::edge_vertex_shapes(edg);
        // `cxx:312`: `ConnectLast = Vl.IsSame(vtx)`.
        let connect_last = vl.as_ref().is_some_and(|l| is_same(l, vtx));
        // `cxx:313-351`.
        for eee in crate::topo_tools::explorer_shapes(
            &fac.oriented(crate::abs::Orientation::Forward),
            ShapeType::Edge,
        ) {
            let orient = eee.orientation();
            if is_same(&eee, edg) {
                continue;
            }
            let (vf, vl) = super::edge_vertex_shapes(&eee);
            let touches = |v: &Option<TopoShape>| v.as_ref().is_some_and(|x| is_same(x, vtx));
            if !(touches(&vf) || touches(&vl)) || self.to_remove.contains(&eee) {
                continue;
            }
            fuse_edge = false;
            // `cxx:325-348`.
            let forward_arm = (touches(&vf) && orient == crate::abs::Orientation::Forward)
                || (touches(&vl) && orient == crate::abs::Orientation::Reversed);
            let bound = if forward_arm == connect_last {
                crate::abs::Orientation::Forward
            } else {
                crate::abs::Orientation::Reversed
            };
            self.dont_fuse
                .bind(edg, &fac.oriented(bound));
            break;
        }
        // `cxx:352-355`.
        if fuse_edge {
            self.to_remove.add(vtx);
        }
    }

    /// `cxx:359-428`: the open edge case, where the fusion depends on the
    /// *ancestor* faces.
    fn classify_open(
        &mut self,
        fac: &TopoShape,
        edg: &TopoShape,
        vtx: &TopoShape,
        edgbis: &TopoShape,
    ) {
        // `cxx:362-368`.
        if !self.ee_map.is_bound(edg) {
            self.ee_map.bind(edg, Vec::new());
        }
        if let Some(list) = self.ee_map.change(edg) {
            list.push(edgbis.clone());
        }
        // `cxx:369-385`: the reference re-scans the list it just appended to
        // through a `const NCollection_List&` alias of the same map entry, so
        // `Lit.Value().IsSame(edgbis)` is always true and the second `Append`
        // never runs. The port keeps the single append.
        // `cxx:387-388`.
        let ancestors = self
            .ef_map
            .find_from_key(edg)
            .map(|a| a.to_vec())
            .unwrap_or_default();
        let fuse_edge = super::statics::to_fuse_edge_face_vertex(edg, fac, vtx, &self.to_remove);
        if !fuse_edge {
            // `cxx:389-391`.
            self.dont_fuse.bind(edg, fac);
            return;
        }
        // `cxx:393-426`.
        for anc in ancestors {
            if is_same(&anc, fac) {
                continue;
            }
            if self.ff_map.is_bound(&anc) {
                let ok = super::statics::to_fuse_edge_face_vertex(edg, &anc, vtx, &self.to_remove);
                if ok {
                    self.to_remove.add(vtx);
                } else {
                    if self.to_remove.contains(vtx) {
                        self.to_remove.remove(vtx);
                    }
                    self.dont_fuse.bind(edg, &anc);
                }
            } else {
                // `cxx:415-423`: a face to rebuild.
                self.rebuild_face.push(anc.clone());
                if self.to_remove.contains(vtx) {
                    self.to_remove.remove(vtx);
                }
                self.dont_fuse.bind(edg, &anc);
            }
            break;
        }
    }
}

/// `TopoDS::Edge` on a raw shape, panicking like the reference's unchecked
/// down cast only when the shape is not an edge.
fn edge_of(e: &TopoShape) -> crate::shape::Edge {
    crate::shape::Edge::wrap(e.clone()).expect("LocOpe_Generator: expected an edge")
}
