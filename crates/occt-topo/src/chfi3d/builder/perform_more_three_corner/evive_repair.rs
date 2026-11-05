//! The `Evive` repair of `PerformMoreThreeCorner` (`cxx:1679-1782`).
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_CnCrn.cxx:1679-1782`
//! (OCCT 8.0.0).
//!
//! When the corner vertex sits on an edge that has no vertex at its extremity,
//! the sharp edge of the ring is not the edge the stripe sees: the common point
//! of the neighbouring surf data carries an arc that ends at the corner vertex
//! and is tangent (`|angle - pi| < 0.01`) to the sharp edge. That arc replaces
//! `Evive(ic)`, its two faces replace the `Fvive` pair of the slot, and the
//! parameter of the common point on the arc is copied into `p` for the two
//! neighbouring slots with `i = 1` (cts16288).
//!
//! Translation notes:
//! - `cp1` / `cp2` are `ChFiDS_CommonPoint` **copies** in OCCT (`cxx:1693`),
//!   taken from the mutable `ChangeVertex(...)` reference; the port reads the
//!   same cell through [`crate::chfi_ds::SurfData::vertex`] and clones it,
//! - the two `DStr.AddShape(F)` calls of one OCCT branch always add the same
//!   face, so the port adds it once and writes the answered index to both
//!   `numfa` cells.

use std::f64::consts::PI;

use crate::chfi3d::builder_helpers::{angle_edge, cherche_vertex, edge_common_faces};
use crate::chfi_ds::{ChFiDSMap, CommonPoint};
use crate::shape::{Edge, TopoShape};
use crate::topo_tools::is_same;
use crate::topopebrep::ds::HDataStructure;

use super::surf_index::indices;
use super::MoreThreeCorner;

impl MoreThreeCorner {
    /// `cxx:1679-1782`, entered only when `!deuxconges && !isOnSameDiff`
    /// (`cxx:1683`).
    pub(super) fn repair_evive(
        &mut self,
        ds: &mut HDataStructure,
        ef_map: &ChFiDSMap,
    ) -> Result<(), String> {
        let nedge = self.nedge as usize;
        for ic in 0..nedge {
            if !self.sharp[ic] {
                continue;
            }
            let (icplus, icmoins) = indices(nedge, ic);
            let arc = self.edge_at(ic)?;

            // `cxx:1712-1728` / `:1757-1775`: the common point of the
            // neighbouring surf data, and the `Evive` rewrite it may cause.
            self.try_arc(ic, icplus, icmoins, &arc, true, ds, ef_map)?;
            self.try_arc(ic, icmoins, icplus, &arc, false, ds, ef_map)?;
        }
        Ok(())
    }

    /// One arm of the `cxx:1694-1740` (`neighbour = icplus`, `jfp = 3 - jf`) or
    /// `cxx:1739-1780` (`neighbour = icmoins`, `jfp = jf`) block.
    ///
    /// `other` is the second neighbour of `ic`, used by the `Fvive` test that
    /// decides which of the two faces of the arc is written.
    #[allow(clippy::too_many_arguments)]
    fn try_arc(
        &mut self,
        ic: usize,
        neighbour: usize,
        other: usize,
        arc: &Edge,
        plus_side: bool,
        ds: &mut HDataStructure,
        ef_map: &ChFiDSMap,
    ) -> Result<(), String> {
        if self.sharp[neighbour] {
            return Ok(());
        }
        let isfirst = self.sens[neighbour] == 1;
        let jfp = if plus_side {
            // `jfp = 3 - jf.Value(icplus);` (`cxx:1698`) — a function-level
            // local of `PerformMoreThreeCorner`, written even when the arc test
            // below fails, and read again by the `SearchFD` stage (`cxx:1950`,
            // `:1964`) and the extremity fill (`cxx:2057`).
            self.jfp = 3 - self.jf[neighbour];
            self.jfp
        } else {
            self.jf[neighbour]
        };
        let cp = self.common_point(neighbour, isfirst, jfp)?;
        let Some(cp_arc_shape) = cp.arc().cloned() else {
            return Ok(());
        };
        // `ChFi3d_cherche_vertex(Arc, cp1.Arc(), Vcom, trouve);` (`cxx:1705`).
        let cp_arc = Edge::wrap(cp_arc_shape.clone())
            .ok_or_else(|| "PerformMoreThreeCorner: the common point arc is not an edge".to_string())?;
        let mut angedg = PI;
        if let Some(vcom) = cherche_vertex(arc, &cp_arc) {
            angedg = angle_edge(&vcom, arc, &cp_arc)?.abs();
        }
        if is_same(&cp_arc_shape, &arc.0) || (angedg - PI).abs() >= 0.01 {
            return Ok(());
        }

        // `Evive.SetValue(ic, cp1.Arc());` (`cxx:1712`).
        self.evive[ic] = Some(cp_arc_shape.clone());
        // `ChFi3d_edge_common_faces(myEFMap(cp1.Arc()), F1, F2);` (`cxx:1713`).
        let (f1, f2) = edge_common_faces(ef_map.get_from_key(&cp_arc_shape))
            .ok_or_else(|| "PerformMoreThreeCorner: the arc has no face".to_string())?;

        // `if (!Fvive.Value(ic, neighbour).IsSame(F1) && !...IsSame(F2))
        //  { if (Fvive.Value(ic, other).IsSame(F2)) write F1 else write F2; }`
        // (`cxx:1714-1728` for the `icplus` arm, `:1758-1775` for `icmoins`).
        let neighbour_cell = self.fvive[ic][neighbour]
            .clone()
            .ok_or_else(|| "PerformMoreThreeCorner: Fvive is null".to_string())?;
        if !is_same(&neighbour_cell, &f1.0) && !is_same(&neighbour_cell, &f2.0) {
            let partner_cell = self.fvive[ic][other]
                .clone()
                .ok_or_else(|| "PerformMoreThreeCorner: Fvive is null".to_string())?;
            let own = if is_same(&partner_cell, &f2.0) {
                &f1
            } else {
                &f2
            };
            let shape = own.0.clone();
            let num = ds.change_ds().add_shape(&shape);
            self.set_fvive_pair(ic, neighbour, &shape, num);
        }

        // `samedge.SetValue(ic, true); p.SetValue(...); i.SetValue(ic,
        // neighbour, 1);` (`cxx:1731-1735` / `:1776-1779`).
        let param = cp
            .parameter_on_arc()
            .ok_or_else(|| "PerformMoreThreeCorner: the common point has no arc parameter".to_string())?;
        self.samedge[ic] = true;
        self.p[ic][neighbour] = param;
        self.p[ic][other] = param;
        self.i[ic][neighbour] = 1;
        Ok(())
    }

    /// `Fvive.SetValue(ic, n, F); numfa.SetValue(ic, n, DStr.AddShape(F));
    /// Fvive.SetValue(n, ic, F); numfa.SetValue(n, ic, DStr.AddShape(F));`
    fn set_fvive_pair(&mut self, ic: usize, n: usize, face: &TopoShape, num: i32) {
        self.fvive[ic][n] = Some(face.clone());
        self.fvive[n][ic] = Some(face.clone());
        self.numfa[ic][n] = num;
        self.numfa[n][ic] = num;
    }

    /// `CD.Value(ic)->SetOfSurfData()->Value(Index.Value(ic))->ChangeVertex(
    /// isfirst, jfp)` (`cxx:1696-1701`, `:1741-1745`), read as a copy.
    fn common_point(&self, ic: usize, isfirst: bool, jfp: i32) -> Result<CommonPoint, String> {
        self.common_point_at(ic, isfirst, jfp)
    }
}
