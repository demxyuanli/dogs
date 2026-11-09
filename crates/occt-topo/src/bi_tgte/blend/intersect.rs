//! `BiTgte_Blend::Intersect` (`BiTgte_Blend.cxx:2511-2664`): intersect one
//! newly built offset face with every offset already in `myMapSF`, and
//! propagate the tubes over the restrictions the intersection reached.
//!
//! Control flow, all from the `.cxx`:
//!
//! * `:2520` `Box1 = MapSBox(Face)`.
//! * `:2525-2526` `F1surBordLibre`: the offset comes from an edge that is a
//!   stopping face.
//! * `:2528-2540` walk `myMapSF`; `Box1.IsOut(MapSBox(F2))` rejects the pair,
//!   `Inter.IsDone(Face, F2)` skips an already computed pair.
//! * `:2546-2561` two tubes on free borders are not intersected.
//! * `:2566-2573` a tube is not intersected with a neighbouring face whose
//!   wire already holds it (`IsInFace`).
//! * `:2575` `Inter.FaceInter(Face, F2, myInitOffsetFace)`.
//! * `:2586-2661` for every common descendant edge and each of its two
//!   vertices, `IsOnRestriction` on both faces decides whether `Add` records
//!   the restriction in `myEdges` and whether the outer `JenRajoute` loop of
//!   `ComputeCenters` runs another pass.
//!
//! Deviations:
//! - `myMapSF` is [`OffsetMap`], which keeps insertion order so that both the
//!   key shape and the walk order are available; the reference's
//!   `NCollection_DataMap` walks in hash-bucket order. The key set looked at is
//!   the same; only the order of the pair-wise intersections differs.
//! - `MapSBox(Face)` raises `Standard_NoSuchObject` in the reference; the port
//!   answers `false` for the impossible unbound case.

use std::collections::{HashMap, HashSet};

use occt_core::bnd::box3d::BndBox;

use crate::abs::ShapeType;
use crate::brep_offset::{BRepOffsetInter3d, BRepOffsetOffset};
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{edge_vertices, shape_map_key, ShapeMapKey};

use super::statics::{add, is_in_face, is_on_restriction};
use super::Blend;

impl Blend {
    /// `bool BiTgte_Blend::Intersect(Init, Face, MapSBox, OF1, Inter)`
    /// (`cxx:2511-2664`).
    pub(super) fn intersect(
        &mut self,
        init: &TopoShape,
        face: &Face,
        map_sbox: &HashMap<ShapeMapKey, BndBox>,
        of1: &BRepOffsetOffset,
        inter: &mut BRepOffsetInter3d,
    ) -> bool {
        let mut jen_rajoute = false;

        // `cxx:2520`.
        let Some(box1) = map_sbox.get(&shape_map_key(&face.0)).cloned() else {
            return false;
        };

        // `cxx:2525-2526`.
        let f1_sur_bord_libre = of1
            .initial_shape()
            .map(|s| s.shape_type() == ShapeType::Edge && self.my_stop_faces.contains(&shape_map_key(s)))
            .unwrap_or(false);

        // `cxx:2528-2534`: the snapshot keeps `myEdges` free for `Add`.
        let entries: Vec<(TopoShape, BRepOffsetOffset)> = self
            .my_map_sf
            .iter()
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();

        for (it_key, of2) in entries.into_iter() {
            let Some(f2) = of2.face().cloned() else {
                continue;
            };

            // `cxx:2536-2539`.
            match map_sbox.get(&shape_map_key(&f2.0)) {
                Some(b) if box1.is_out_box(b) => continue,
                Some(_) => {}
                None => continue,
            }

            // `cxx:2541-2544`.
            if inter.is_done(face, &f2) {
                continue;
            }

            // `cxx:2546-2561`.
            let f2_sur_bord_libre = of2
                .initial_shape()
                .map(|s| s.shape_type() == ShapeType::Edge && self.my_stop_faces.contains(&shape_map_key(s)))
                .unwrap_or(false);
            if f1_sur_bord_libre && f2_sur_bord_libre {
                continue;
            }

            // `cxx:2566-2573`.
            if init.shape_type() == ShapeType::Edge
                && it_key.shape_type() == ShapeType::Face
                && is_in_face(&Edge(init.clone()), &Face(it_key.clone()))
            {
                continue;
            }

            // `cxx:2575`.
            inter.face_inter(face, &f2, &self.my_init_offset_face);

            // `cxx:2585-2661`.
            let mut l_int: Vec<TopoShape> = Vec::new();
            if !inter
                .as_des()
                .has_common_descendant(&face.0, &f2.0, &mut l_int)
            {
                continue;
            }
            let mut done: HashSet<ShapeMapKey> = HashSet::new();
            for cur_e in l_int {
                let cur_e = Edge(cur_e);
                let (v1, v2) = edge_vertices(&cur_e);
                // `cxx:2593-2658`: the same body for `V1` and `V2`.
                for v in [v1, v2].into_iter().flatten() {
                    // `cxx:2596` / `:2624`: `if (Done.Add(V1))`.
                    if !done.insert(shape_map_key(&v.0)) {
                        continue;
                    }
                    let mut e1: Option<Edge> = None;
                    let mut e2: Option<Edge> = None;
                    let is_on_r1 = is_on_restriction(&v, &cur_e, face, &mut e1);
                    let is_on_r2 = is_on_restriction(&v, &cur_e, &f2, &mut e2);
                    let both = is_on_r1 && is_on_r2;

                    if is_on_r1 && !self.my_stop_faces.contains(&shape_map_key(init)) {
                        if let Some(e1) = e1.as_ref() {
                            add(e1, &mut self.my_edges, init, of1, &self.my_analyse, both);
                            jen_rajoute = true;
                        }
                    }
                    if is_on_r2 && !self.my_stop_faces.contains(&shape_map_key(&it_key)) {
                        if let Some(e2) = e2.as_ref() {
                            add(e2, &mut self.my_edges, &it_key, &of2, &self.my_analyse, both);
                            jen_rajoute = true;
                        }
                    }
                }
            }
        }

        jen_rajoute
    }
}
