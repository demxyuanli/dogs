//! The vertex-distance scan and the parameter offset of
//! `PerformMoreThreeCorner`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_CnCrn.cxx:2064-2196`
//! (OCCT 8.0.0).
//!
//! Two blocks:
//! - `cxx:2064-2107` walks the non-sharp slots, evaluates the fillet surface at
//!   the two pairs `(i(ic, icmoins), p(ic, icmoins))` and
//!   `(i(ic, icplus), p(ic, icplus))` and records the distance of each point to
//!   the corner vertex in `dist1` / `dist2`, tracking the maximum in
//!   `distance`;
//! - `cxx:2109-2196` moves the parameters away from the vertex by that maximum
//!   distance: a sharp slot converts the 3d distance to a parametric one with
//!   `BRepAdaptor_Curve::Resolution(0.01)` (zero for a regular tangent edge) and
//!   shifts `p(ic, icmoins)` along the edge, while a non-sharp slot shifts by
//!   `distance - dist` and clears the `oksea` flag of the side it invalidates.
//!
//! Translation notes:
//! - `dist1` / `dist2` are `NCollection_Array1<double>(0, size)`; the port keeps
//!   them in [`MoreThreeCorner::dist1`] / [`MoreThreeCorner::dist2`],
//! - `inters` is written at `cxx:2151` and `:2170` but its only reader is the
//!   `if (!inters)` block of `cxx:2004-2062`, which has already run: the port
//!   keeps the writes (on [`MoreThreeCorner::inters`]) so the state matches
//!   OCCT's, and notes that they are dead in OCCT too,
//! - `TopExp::FirstVertex(E)` (`cxx:2129`) is the first entry of
//!   `TopoDS_Iterator(E)`; `IsSame(V1)` is [`crate::topo_tools::is_same`],
//! - `BRepAdaptor_Curve::Resolution(0.01)` (`cxx:2122`) is
//!   [`crate::topopebrep::tool::adaptor::CurveAdaptor::resolution`],
//! - the `distance * 100 * Resolution` arithmetic and the `1.3 * dist`
//!   threshold are OCCT's own; neither is a port-side heuristic.

use crate::brep_tool::BRepTool;
use crate::shape::{TopoShape, Vertex};
use crate::topo_tools::{edge_vertices, is_same};
use crate::topopebrep::ds::HDataStructure;
use crate::topopebrep::tool::adaptor::CurveAdaptor;

use super::calcul::calcul_p2d_on_surf;
use super::surf_index::{indices, surf_index, SurfType};
use super::MoreThreeCorner;

impl MoreThreeCorner {
    /// `cxx:2064-2196`: the `dist1` / `dist2` scan followed by the parameter
    /// offset. `nconges` is `PerformMoreThreeCorner`'s own argument
    /// (`cxx:1237`), read by the `nconges != 1` branch of `cxx:2172`.
    pub(super) fn vertex_distance_offsets(
        &mut self,
        ds: &HDataStructure,
        vtx: &TopoShape,
        nconges: usize,
    ) -> Result<(), String> {
        let nedge = self.nedge as usize;
        let vertex = Vertex::wrap(vtx.clone())
            .ok_or_else(|| "PerformMoreThreeCorner: V1 is not a vertex".to_string())?;
        let sommet = BRepTool::vertex_point(&vertex);

        // `double distance = 0.;` (`cxx:2067`) — shared by the two blocks.
        let mut distance = 0.0f64;

        // `cxx:2069-2107`.
        if !self.deuxconges {
            for ic in 0..nedge {
                let (icplus, icmoins) = indices(nedge, ic);
                if self.sharp[ic] {
                    // `dist1.SetValue(ic, 0); dist2.SetValue(ic, 0);`
                    self.dist1[ic] = 0.0;
                    self.dist2[ic] = 0.0;
                    continue;
                }
                // `jfp = 3 - jf.Value(ic);` (`cxx:2081`) — the function-level
                // local again, not a fresh one.
                self.jfp = 3 - self.jf[ic];

                // `cxx:2082-2091`.
                let d1 = self.surface_distance_at(
                    ds,
                    ic,
                    self.i[ic][icmoins],
                    self.jfp,
                    self.p[ic][icmoins],
                    &sommet,
                )?;
                self.dist1[ic] = d1;
                if d1 > distance {
                    distance = d1;
                }

                // `cxx:2093-2103`.
                let d2 = self.surface_distance_at(
                    ds,
                    ic,
                    self.i[ic][icplus],
                    self.jf[ic],
                    self.p[ic][icplus],
                    &sommet,
                )?;
                self.dist2[ic] = d2;
                if d2 > distance {
                    distance = d2;
                }
            }
        }

        // `cxx:2109-2196`.
        if !self.deuxconges && !self.deuxcgnontg {
            for ic in 0..nedge {
                let (icplus, icmoins) = indices(nedge, ic);
                if self.sharp[ic] {
                    // `BRepAdaptor_Curve C(TopoDS::Edge(Evive.Value(ic)));`
                    // (`cxx:2119`).
                    let edge = self.edge_at(ic)?;
                    let curve = CurveAdaptor::of_edge(&edge);
                    // `ec = distance * 100 * C.Resolution(0.01);`
                    // (`cxx:2121-2128`).
                    let ec = if !self.tangentregul[ic] {
                        distance * 100.0 * curve.resolution(0.01)
                    } else {
                        0.0
                    };
                    // `if (TopExp::FirstVertex(...).IsSame(V1))` (`cxx:2129`).
                    let first = edge_vertices(&edge).0;
                    let at_first = first
                        .as_ref()
                        .is_some_and(|v| is_same(v, vtx));
                    let para = if at_first {
                        self.p[ic][icmoins] + ec
                    } else {
                        self.p[ic][icmoins] - ec
                    };
                    self.p[ic][icmoins] = para;
                    // it is necessary to be on to remain on the edge.
                    self.p[ic][icplus] = self.p[ic][icmoins];
                    continue;
                }
                if self.distmini {
                    continue;
                }

                // `dist = dist1.Value(ic);` + `cxx:2144-2165`.
                let dist = self.dist1[ic];
                if !self.oksea[icmoins] || (self.oksea[icmoins] && distance > 1.3 * dist) {
                    let ec = distance - dist;
                    if self.oksea[icmoins] {
                        self.oksea[icmoins] = false;
                        self.inters = false;
                    }
                    let para = if self.sens[ic] == 1 {
                        self.p[ic][icmoins] + ec
                    } else {
                        self.p[ic][icmoins] - ec
                    };
                    self.p[ic][icmoins] = para;
                }

                // `dist = dist2.Value(ic);` + `cxx:2167-2193`.
                let dist = self.dist2[ic];
                if !self.oksea[ic] || (self.oksea[ic] && distance > 1.3 * dist) {
                    if self.oksea[ic] {
                        self.oksea[ic] = false;
                        self.inters = false;
                    }
                    if nconges != 1 {
                        let parold = self.p[ic][icplus];
                        let parnew = self.p[ic][icmoins];
                        if self.sens[ic] == 1 {
                            if parnew > parold {
                                self.p[ic][icplus] = self.p[ic][icmoins];
                            }
                        } else if parnew < parold {
                            self.p[ic][icplus] = self.p[ic][icmoins];
                        }
                    }
                }
            }
        }
        Ok(())
    }

    /// `Calcul_P2dOnSurf(CD.Value(ic), jfposit, ind, param, p2)` followed by
    /// `DStr.Surface(SurfIndex(CD, ic, ind, ChFiSURFACE)).Surface()->D0(p2.X(),
    /// p2.Y(), pic)` and `sommet.Distance(pic)` (`cxx:2082-2087`,
    /// `:2093-2098`).
    fn surface_distance_at(
        &self,
        ds: &HDataStructure,
        slot: usize,
        ind: i32,
        jf: i32,
        param: f64,
        sommet: &occt_core::gp::GpPnt,
    ) -> Result<f64, String> {
        let handle = self.cd[slot]
            .clone()
            .ok_or_else(|| "PerformMoreThreeCorner: CD is null".to_string())?;
        let p2 = {
            let guard = handle
                .read()
                .map_err(|_| "PerformMoreThreeCorner: poisoned stripe".to_string())?;
            calcul_p2d_on_surf(&guard, jf, ind, param)?
        };
        let indice = surf_index(&self.cd, slot, ind as usize, SurfType::ChFiSurface)?;
        let surface = ds
            .surface(indice)
            .surface()
            .cloned()
            .ok_or_else(|| format!("PerformMoreThreeCorner: surface {indice} is null"))?;
        Ok(sommet.distance(&surface.d0(p2.x(), p2.y())))
    }
}
