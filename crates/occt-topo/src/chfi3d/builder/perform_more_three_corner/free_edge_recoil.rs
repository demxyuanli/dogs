//! The first free edge and the minimum recoil distance of
//! `PerformMoreThreeCorner`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_CnCrn.cxx:1786-1830`
//! (OCCT 8.0.0).
//!
//! `libre(ic)` marks the slot whose `Evive` is one of the two free borders (the
//! `Evive` repair of `cxx:1679-1782` can have replaced it, hence the re-scan),
//! and `distmin` is the shortest chord of the ring, divided by 1.5 — the recoil
//! distance that the corner curves of the later stages may not exceed.
//!
//! Translation notes:
//! - `gp_Pnt som = BRep_Tool::Pnt(V1)` (`cxx:1800`), `pic` and `p2`
//!   (`cxx:1800-1801`) are only read by the `distmini` test of `cxx:1986`;
//!   they live in [`super::search_fd_stage`], which re-reads the vertex point
//!   there, so this stage does not keep a second copy.

use crate::brep_tool::BRepTool;
use crate::topo_tools::edge_vertices;

use super::{spine_edge_at, MoreThreeCorner};

impl MoreThreeCorner {
    /// `cxx:1786-1830`.
    pub(super) fn free_edge_and_recoil(&mut self) -> Result<(), String> {
        let nedge = self.nedge as usize;

        // `trouve = false; for (ic = 0; ic < nedge && !trouve; ic++) { ... }`
        // (`cxx:1786-1797`).
        for ic in 0..nedge {
            let ecom = self.edge_at(ic)?;
            if self.is_free_border(&ecom) {
                self.libre[ic] = true;
                break;
            }
        }

        // `distmin = 1.e30; for (ic = 0; ic < nedge; ic++) { ... }`
        // (`cxx:1799-1830`), with `bool distmini = false;` of `cxx:1799`.
        self.distmini = false;
        let mut distmin = 1.0e30_f64;
        for ic in 0..nedge {
            let edgemin = if self.sharp[ic] {
                self.edge_at(ic)?
            } else {
                let stripe = self.cd[ic].clone().ok_or_else(|| {
                    "PerformMoreThreeCorner: CD is null for a non-sharp edge".to_string()
                })?;
                spine_edge_at(&stripe, self.sens[ic] == 1)?
            };
            let (v, v2) = edge_vertices(&edgemin);
            let v = v.ok_or_else(|| {
                "PerformMoreThreeCorner: the recoil edge has no first vertex".to_string()
            })?;
            let v2 = v2.ok_or_else(|| {
                "PerformMoreThreeCorner: the recoil edge has no last vertex".to_string()
            })?;
            let dst = BRepTool::vertex_point(&v).distance(&BRepTool::vertex_point(&v2)) / 1.5;
            if dst < distmin {
                distmin = dst;
            }
        }
        self.distmin = distmin;
        Ok(())
    }
}
