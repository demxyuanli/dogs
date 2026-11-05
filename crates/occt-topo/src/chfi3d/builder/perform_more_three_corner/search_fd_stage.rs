//! The `ChFi3d_SearchFD` stage of `PerformMoreThreeCorner` and the extremity
//! fill that follows it.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_CnCrn.cxx:1832-2063`
//! (OCCT 8.0.0).
//!
//! For every non-sharp slot of the ring the stage looks for the pair of surf
//! data that face each other (`ChFi3d_SearchFD`), prunes the surf-data
//! sequences that the sliding window walked over (`RemoveSD`) and rewrites
//! `Index` to the surf data that carries the corner vertex, then stores the
//! intersection parameters in `p` / `i` / `jf`. When a slot has no
//! intersection at all, the second half of the stage fills the missing
//! parameters from the extremities of the p-curves instead.
//!
//! Translation notes:
//! - `som` (`cxx:1800`), `pic` and `p2` are declared once in the OCCT function
//!   and reused here; the port keeps them as locals,
//! - `Fproj` (`cxx:1318`) is an `NCollection_Sequence<TopoDS_Shape>` that the
//!   face walks append to; those faces are read from `cxx:2595` on, so the port
//!   stores them in [`MoreThreeCorner::fproj`] for that stage
//!   ([`super::multi_face_hermite`] and the `plate_*` stages),
//! - `strip` is assigned in all three arms of the `i1 != Index(ic)` block in
//!   OCCT (`cxx:1904`, `:1919`, `:1932`); the port clones the handle once,
//! - the `3 - jfp` face selector of the `i2 != Index(icplus)` block
//!   (`cxx:1950`, `:1964`) reads the **function-level** `jfp`
//!   ([`MoreThreeCorner::jfp`]), not a fresh value,
//! - OCCT's `NCollection_Array2::Value(ic, icplus)` is `[ic][icplus]` in the
//!   port (both dimensions `0..size`).

use std::f64::consts::PI;

use crate::brep_tool::BRepTool;
use crate::chfi3d::builder_0_helpers::{ch_fi3d_index_of_surf_data, ch_fi3d_search_fd};
use crate::chfi3d::builder_helpers::angle_edge;
use crate::shape::{Vertex, TopoShape};
use crate::topopebrep::ds::HDataStructure;

use super::calcul::{calcul_param, calcul_p2d_on_surf};
use super::remove_surf_data::remove_sd;
use super::surf_index::{indices, surf_index, SurfType};
use super::MoreThreeCorner;

impl MoreThreeCorner {
    /// `cxx:1832-2063`: the `ChFi3d_SearchFD` loop followed by the `!inters`
    /// extremity fill.
    pub(super) fn intersect_stripes(
        &mut self,
        ds: &HDataStructure,
        vtx: &TopoShape,
    ) -> Result<(), String> {
        let nedge = self.nedge as usize;
        let vertex = Vertex::wrap(vtx.clone())
            .ok_or_else(|| "PerformMoreThreeCorner: V1 is not a vertex".to_string())?;
        let som = BRepTool::vertex_point(&vertex);

        // `bool inters = true;` (`cxx:1833`).
        let mut inters = true;

        // `cxx:1834-2001`.
        for ic in 0..nedge {
            let (icplus, _icmoins) = indices(nedge, ic);
            if self.sharp[ic] || self.sharp[icplus] {
                // `oksea.SetValue(ic, false);` (`cxx:1838-1840`).
                self.oksea[ic] = false;
            } else {
                // `int jf1 = 0; int i1 = 0, i2 = 0; double pa1 = 0., pa2;`
                // (`cxx:1844-1846`).
                let mut jf1 = 0i32;
                let mut i1 = 0i32;
                let mut i2 = 0i32;
                let mut pa1 = 0.0f64;
                let mut pa2 = 0.0f64;

                // `angedg = abs(AngleEdge(V1, Evive(ic), Evive(icplus)));`
                // (`cxx:1851`).
                let e1 = self.edge_at(ic)?;
                let e2 = self.edge_at(icplus)?;
                let angedg = angle_edge(&vertex, &e1, &e2)?.abs();

                // `if (abs(angedg - pi) > 0.01) ok = SearchFD(...); else
                //  ok = false;` (`cxx:1853-1876`, cts60046: two tangent edges
                //  are not intersected).
                let ok = if (angedg - PI).abs() > 0.01 {
                    let cd1 = self.cd[ic]
                        .clone()
                        .ok_or_else(|| "PerformMoreThreeCorner: CD is null".to_string())?;
                    let cd2 = self.cd[icplus]
                        .clone()
                        .ok_or_else(|| "PerformMoreThreeCorner: CD is null".to_string())?;
                    let ans = {
                        let g1 = cd1
                            .read()
                            .map_err(|_| "PerformMoreThreeCorner: poisoned stripe".to_string())?;
                        let g2 = cd2
                            .read()
                            .map_err(|_| "PerformMoreThreeCorner: poisoned stripe".to_string())?;
                        ch_fi3d_search_fd(
                            ds,
                            &g1,
                            &g2,
                            self.sens[ic],
                            self.sens[icplus],
                            self.index[ic],
                            self.index[icplus],
                        )?
                    };
                    jf1 = ans.jf1;
                    i1 = ans.i1;
                    i2 = ans.i2;
                    pa1 = ans.p1;
                    pa2 = ans.p2;
                    ans.found
                } else {
                    false
                };

                if ok {
                    // `cxx:1879-1939`: the first stripe walked over surf data
                    // that cannot reach the corner vertex.
                    if i1 != self.index[ic] {
                        self.prune_stripe_index(ds, ic, &mut i1, jf1)?;
                    }
                    // `cxx:1941-1980`: the `icplus` mirror, with `3 - jfp`.
                    if i2 != self.index[icplus] {
                        let jfp = 3 - self.jfp;
                        self.prune_stripe_index(ds, icplus, &mut i2, jfp)?;
                    }

                    // `Calcul_P2dOnSurf(CD.Value(ic), jf1, i1, pa1, p2);
                    //  indice = SurfIndex(CD, ic, i1, ChFiSURFACE);
                    //  DStr.Surface(indice).Surface()->D0(p2.X(), p2.Y(), pic);`
                    // (`cxx:1981-1984`).
                    let p2 = {
                        let handle = self.cd[ic].clone().ok_or_else(|| {
                            "PerformMoreThreeCorner: CD is null".to_string()
                        })?;
                        let guard = handle
                            .read()
                            .map_err(|_| "PerformMoreThreeCorner: poisoned stripe".to_string())?;
                        calcul_p2d_on_surf(&guard, jf1, i1, pa1)?
                    };
                    let indice = surf_index(&self.cd, ic, i1 as usize, SurfType::ChFiSurface)?;
                    let surface = ds
                        .surface(indice)
                        .surface()
                        .cloned()
                        .ok_or_else(|| {
                            format!("PerformMoreThreeCorner: surface {indice} is null")
                        })?;
                    let pic = surface.d0(p2.x(), p2.y());
                    // `if (pic.Distance(som) > distmin) distmini = true;`
                    // (`cxx:1986-1989`).
                    if pic.distance(&som) > self.distmin {
                        self.distmini = true;
                    }

                    // `cxx:1990-1994`.
                    self.jf[ic] = jf1;
                    self.i[ic][icplus] = i1;
                    self.i[icplus][ic] = i2;
                    self.p[ic][icplus] = pa1;
                    self.p[icplus][ic] = pa2;
                }
                self.oksea[ic] = ok;
            }
            // `if (!oksea.Value(ic)) inters = false;` (`cxx:1996-1999`).
            if !self.oksea[ic] {
                inters = false;
            }
        }

        // `if (!inters) { ... }` (`cxx:2004-2062`): fill the parameters that
        // have no intersection from the p-curve extremities.
        if !inters {
            for ic in 0..nedge {
                let (icplus, icmoins) = indices(nedge, ic);
                let (icplus2, _) = indices(nedge, icplus);
                if self.oksea[ic] {
                    continue;
                }
                // `cxx:2010-2041`: the `ic` side.
                if self.sharp[ic] {
                    if !self.samedge[ic] {
                        let edge = self.edge_at(ic)?;
                        let para = BRepTool::parameter_of_vertex(&edge.0, vtx);
                        self.p[ic][icplus] = para;
                        self.i[ic][icplus] = 1;
                    }
                } else {
                    let isfirst = self.sens[ic] == 1;
                    let handle = self.cd[ic]
                        .clone()
                        .ok_or_else(|| "PerformMoreThreeCorner: CD is null".to_string())?;
                    let mut sense = 0;
                    let num = ch_fi3d_index_of_surf_data(vtx, &handle, &mut sense)?;
                    self.i[ic][icplus] = num;
                    let para = if self.oksea[icmoins] {
                        self.p[ic][icmoins]
                    } else {
                        let guard = handle.read().map_err(|_| {
                            "PerformMoreThreeCorner: poisoned stripe".to_string()
                        })?;
                        calcul_param(&guard, self.jf[ic], num, isfirst)?
                    };
                    self.p[ic][icplus] = para;
                }
                // `cxx:2043-2060`: the `icplus` side; this arm writes `jfp`.
                if self.sharp[icplus] {
                    if !self.samedge[icplus] {
                        let edge = self.edge_at(icplus)?;
                        let para = BRepTool::parameter_of_vertex(&edge.0, vtx);
                        self.p[icplus][ic] = para;
                        self.i[icplus][ic] = 1;
                    }
                } else {
                    let isfirst = self.sens[icplus] == 1;
                    let handle = self.cd[icplus]
                        .clone()
                        .ok_or_else(|| "PerformMoreThreeCorner: CD is null".to_string())?;
                    let mut sense = 0;
                    let num = ch_fi3d_index_of_surf_data(vtx, &handle, &mut sense)?;
                    self.i[icplus][ic] = num;
                    let para = if self.oksea[icplus] {
                        self.p[icplus][icplus2]
                    } else {
                        self.jfp = 3 - self.jf[icplus];
                        let guard = handle.read().map_err(|_| {
                            "PerformMoreThreeCorner: poisoned stripe".to_string()
                        })?;
                        calcul_param(&guard, self.jfp, num, isfirst)?
                    };
                    self.p[icplus][ic] = para;
                }
            }
        }

        self.inters = inters;
        Ok(())
    }

    /// The `i1 != Index(ic)` / `i2 != Index(icplus)` block
    /// (`cxx:1879-1939` and `:1941-1980`): append the faces the walk crossed to
    /// `Fproj`, drop the surf data between the old and the new index, re-read
    /// the index of the vertex on the stripe, and answer the new index.
    ///
    /// `face_selector` is `3 - jf` for the first stripe and `3 - jfp` for the
    /// second one, exactly as OCCT writes the two blocks.
    fn prune_stripe_index(
        &mut self,
        ds: &HDataStructure,
        slot: usize,
        ind: &mut i32,
        face_selector: i32,
    ) -> Result<(), String> {
        let handle = self.cd[slot]
            .clone()
            .ok_or_else(|| "PerformMoreThreeCorner: CD is null".to_string())?;
        let old = self.index[slot];
        let sens = self.sens[slot];
        // `if (sens(ic) == 1) { ideb = Index(ic); ifin = i1 - 1; } else
        //  { ideb = i1 + 1; ifin = Index(ic); }` (`cxx:1882-1889`).
        let (ideb, ifin) = if sens == 1 {
            (old, *ind - 1)
        } else {
            (*ind + 1, old)
        };
        let face_type = if (3 - face_selector) == 1 {
            SurfType::Face1
        } else {
            SurfType::Face2
        };
        // `if (i1 < Index(ic)) { for (nb = i1; nb <= Index(ic); nb++) ... }`
        // (`cxx:1891-1908`) and its `i1 > Index(ic)` mirror (`:1910-1927`).
        if *ind < old {
            for nb in *ind..=old {
                let iface = surf_index(&self.cd, slot, nb as usize, face_type)?;
                self.push_projected_face(ds, iface)?;
            }
        } else if *ind > old {
            for nb in old..=*ind {
                let iface = surf_index(&self.cd, slot, nb as usize, face_type)?;
                self.push_projected_face(ds, iface)?;
            }
        }
        // `RemoveSD(strip, ideb, ifin);` (`cxx:1931` / `:1975`).
        {
            let mut guard = handle
                .write()
                .map_err(|_| "PerformMoreThreeCorner: poisoned stripe".to_string())?;
            remove_sd(&mut guard, ideb.max(1) as usize, ifin.max(0) as usize);
        }
        // `num = IndexOfSurfData(V1, CD.Value(ic), sense);
        //  Index.SetValue(ic, num);` (`cxx:1933-1936`).
        let mut sense = 0;
        let num = ch_fi3d_index_of_surf_data(&self.vtx, &handle, &mut sense)?;
        self.index[slot] = num;
        *ind = num;
        Ok(())
    }

    /// `Fproj.Append(TopoDS::Face(myDS->Shape(iface)));` (`cxx:1907`, `:1922`,
    /// `:1958`, `:1973`).
    fn push_projected_face(&mut self, ds: &HDataStructure, iface: i32) -> Result<(), String> {
        let shape = ds.shape(iface, false).ok_or_else(|| {
            format!("PerformMoreThreeCorner: DS shape {iface} is null")
        })?;
        self.fproj.push(shape);
        Ok(())
    }
}
