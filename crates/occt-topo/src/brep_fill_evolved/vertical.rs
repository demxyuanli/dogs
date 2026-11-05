//! `BRepFill_Evolved::VerticalPerform` (`BRepFill_Evolved.cxx:1418-1520`):
//! the sweep of `myProfile` along the planar spine `Sp` when the profile is a
//! vertical translation of its first parallel.
//!
//! The reference offsets the spine by `DistanceToOZ` of the first profile
//! vertex (`cxx:1456-1459`) and records that parallel as the base wire. Each
//! profile edge then translates the current base by `Altitud(V2) - Altitud(V1)`
//! with a `BRepSweep_Prism` (`cxx:1481`), appends the prism faces to
//! `myShape`, and records in `myMap` what the prism generated from every shape
//! of the base for that edge (`cxx:1489-1518`). `Base` becomes the top of the
//! prism, so the next edge starts where the previous one ended.
//!
//! `myIsDone` is left to the caller (`PrivatePerform`), which sets it.

use occt_core::gp::GpVec;
use occt_core::kernel::geomabs::JoinType;

use crate::abs::{Orientation, ShapeType};
use crate::brep_fill_offset_ancestors::BRepFillOffsetAncestors;
use crate::brep_fill_offset_wire::BRepFillOffsetWire;
use crate::brep_mat2d_bisecting_locus::BRepMat2dBisectingLocus;
use crate::brep_mat2d_link_topo_bilo::BRepMat2dLinkTopoBilo;
use crate::brep_sweep::Prism;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::builder::TopoBuilder;
use crate::shape::{Face, TopoShape, Wire};
use crate::topo_tools::{edges_of, explorer_shapes};

use super::evolved::BRepFillEvolved;
use super::shape_map::ShapeMap;
use super::side::{altitud, distance_to_oz, edge_vertices_oriented, FILL_CONFUSION};

impl BRepFillEvolved {
    /// `VerticalPerform(Sp, Pr, Locus, Link, Join)` (`cxx:1418-1520`).
    ///
    /// `Join` is passed through to `BRepFill_OffsetWire::PerformWithBiLo`
    /// (`cxx:1458`). The prism constructor raises in the reference when
    /// `Alt2 - Alt1` vanishes or the base is null; the port answers `Err` with
    /// the same condition.
    pub fn vertical_perform(
        &mut self,
        sp: &Face,
        pr: &Wire,
        locus: &BRepMat2dBisectingLocus,
        link: &mut BRepMat2dLinkTopoBilo,
        join: JoinType,
    ) -> Result<(), String> {
        // `cxx:1427-1434`.
        self.my_spine = Some(Face(sp.0.oriented(Orientation::Forward)));
        self.my_profile = Some(pr.clone());
        self.my_map.clear();

        let builder = TopoBuilder::new();
        let mut shape = builder.make_compound_of(&[]).0;
        self.my_shape = Some(shape.clone());

        let my_spine = self.my_spine.clone().expect("spine was just set");
        let my_profile = self.my_profile.clone().expect("profile was just set");

        // `cxx:1437-1443`.
        let mut paral = BRepFillOffsetWire::new();
        let mut off_anc = BRepFillOffsetAncestors::new();
        let mut first = true;
        let mut base: Option<TopoShape> = None;

        let mut prof_exp = WireExplorer::new();
        prof_exp.init(&my_profile, None);
        while prof_exp.more() {
            let Some(e) = prof_exp.current().cloned() else {
                prof_exp.next();
                continue;
            };

            // `cxx:1448-1449`.
            let (v1, v2) = edge_vertices_oriented(&e);
            let (Some(v1), Some(v2)) = (v1, v2) else {
                prof_exp.next();
                continue;
            };
            let alt1 = altitud(&v1);
            let alt2 = altitud(&v2);

            // `cxx:1451-1479`.
            if first {
                let mut offset = distance_to_oz(&v1);
                if offset.abs() < FILL_CONFUSION {
                    offset = 0.0;
                }
                paral.perform_with_bilo(&my_spine, offset, locus, link, join, alt1);
                off_anc.perform(&mut paral);

                let paral_shape = paral.shape().cloned();
                for an_edge in paral_shape.as_ref().map(edges_of).unwrap_or_default() {
                    let ae = off_anc.ancestor(&an_edge.0).clone();
                    // `cxx:1467-1477`: bind `V1` for every source of the
                    // parallel, then append the parallel edge.
                    let inner = self
                        .my_map
                        .entry_or_bind(&ae, ShapeMap::<Vec<TopoShape>>::new);
                    if !inner.is_bound(&v1.0) {
                        inner.bind(&v1.0, Vec::new());
                    }
                    inner
                        .find_mut(&v1.0)
                        .expect("V1 was just bound")
                        .push(an_edge.0.clone());
                }
                base = paral_shape;
                first = false;
            }

            // `cxx:1481-1488`: translate the current base to the altitude of
            // `V2`; the prism is not a copy (`Copy = false`).
            let Some(current_base) = base.clone() else {
                return Err("BRepFill_Evolved::VerticalPerform: null base wire".into());
            };
            let prism = Prism::from_vec(
                current_base,
                &GpVec::new(0.0, 0.0, alt2 - alt1),
                false,
                true,
            )?;
            base = prism.last_shape();

            if let Some(prism_shape) = prism.shape() {
                for f in explorer_shapes(&prism_shape, ShapeType::Face) {
                    builder.add(&mut shape, &f);
                }
            }

            // `cxx:1490-1518`: for every source shape of `myMap`, record what
            // the prism generated from its `V1` shapes at `V2` and on `E`.
            let keys: Vec<TopoShape> = self.my_map.keys().cloned().collect();
            for key in keys {
                // OCCT reads `it.Value()(V1)`, which raises when `V1` is
                // unbound; the walk above binds it for every key, so the port
                // reads the bound list and keeps an empty one otherwise.
                let lof = self
                    .my_map
                    .list(&key, &v1.0)
                    .cloned()
                    .unwrap_or_default();

                let inner = self
                    .my_map
                    .entry_or_bind(&key, ShapeMap::<Vec<TopoShape>>::new);
                if !inner.is_bound(&v2.0) {
                    inner.bind(&v2.0, Vec::new());
                }
                if !inner.is_bound(&e.0) {
                    inner.bind(&e.0, Vec::new());
                }

                for os in lof {
                    let last = prism.last_shape_of(&os);
                    let swept = prism.shape_of(&os);
                    let inner = self
                        .my_map
                        .entry_or_bind(&key, ShapeMap::<Vec<TopoShape>>::new);
                    if let Some(last) = last {
                        inner
                            .find_mut(&v2.0)
                            .expect("V2 was just bound")
                            .push(last);
                    }
                    if let Some(swept) = swept {
                        inner.find_mut(&e.0).expect("E was just bound").push(swept);
                    }
                }
            }

            prof_exp.next();
        }

        self.my_shape = Some(shape);
        Ok(())
    }
}
