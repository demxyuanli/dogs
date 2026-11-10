//! The multi-face `CurveHermite` stage of `PerformMoreThreeCorner`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_CnCrn.cxx:2542-2640`
//! (OCCT 8.0.0).
//!
//! Where the regular-edge stage (`cxx:2368-2541`) bridges two stripes across
//! the faces it collects itself, this stage walks the slots that the corner
//! intersection gave up on (`!oksea`) and bridges them across the faces that
//! `Fproj` accumulated during the `i1 != Index(ic)` / `i2 != Index(icplus)`
//! walks (`cxx:1832-2063`).
//!
//! Translation notes:
//! - `Ecom`, `Eproj`, `parcom` and `Fproj` are **function-scope** sequences in
//!   OCCT (`cxx:2544-2546`, `:1318`), cleared before the loop and then grown by
//!   every `ic`, so `CurveHermite` reads `Ecom.Value(nb)` / `parcom.Value(nb)`
//!   from the *cumulative* sequence, not from the entries this `ic` appended.
//!   The port keeps them on [`MoreThreeCorner`] and indexes them the same way,
//! - `Fproj.Remove(nbface)` / `Remove(1)` / `Prepend(face1)` keep the length,
//!   so `nbface` stays `Fproj.Length()` after the two adjustments,
//! - `TopoDS::Face(myDS->Shape(iface))` uses the `FindKeep = true` default
//!   (`TopOpeBRepDS_DataStructure.hxx:254`), the same reading as [`super`]'s
//!   `face_at`.

use super::curve_hermite::{curve_hermite, HermiteSide};
use super::surf_index::{indices, surf_index, SurfType};
use super::MoreThreeCorner;
use crate::brep_tool::BRepTool;
use crate::chfi3d::builder_c1::cherche_edge1;
use crate::shape::Face;
use crate::topo_tools::{edge_vertices, is_same};
use crate::topopebrep::ds::DataStructure;

impl MoreThreeCorner {
    /// `cxx:2542-2640`.
    pub(super) fn multi_face_hermite(&mut self, ds: &DataStructure) -> Result<(), String> {
        // `cxx:2544-2546`: the three function-scope sequences.
        self.ecom.clear();
        self.eproj.clear();
        self.parcom.clear();

        if self.deuxconges {
            return Ok(());
        }

        let nedge = self.nedge as usize;
        // `cxx:2548-2639`.
        for ic in 0..nedge {
            let (icplus, _icmoins) = indices(nedge, ic);
            // `cxx:2557`.
            if self.oksea[ic] {
                continue;
            }

            // `cxx:2559-2560`.
            let mut iface1 = self.numfa[ic][icplus];
            let mut iface2 = self.numfa[icplus][ic];

            // `cxx:2561-2571`.
            if !self.sharp[ic] {
                let surf_type = if self.jf[ic] == 1 {
                    SurfType::Face1
                } else {
                    SurfType::Face2
                };
                iface1 = surf_index(&self.cd, ic, self.i[ic][icplus] as usize, surf_type)?;
            }
            let face1 = self.face_at(ds, iface1)?.ok_or_else(|| {
                format!("PerformMoreThreeCorner: DStr.Shape({iface1}) is not a face")
            })?;

            // `cxx:2574-2586`: the right side swaps the two face slots instead
            // of flipping `jf`, which reaches the same surface as
            // `jfp = 3 - jf(icplus)`.
            if !self.sharp[icplus] {
                let surf_type = if self.jf[icplus] == 1 {
                    SurfType::Face2
                } else {
                    SurfType::Face1
                };
                iface2 = surf_index(&self.cd, icplus, self.i[icplus][ic] as usize, surf_type)?;
            }
            let face2 = self.face_at(ds, iface2)?.ok_or_else(|| {
                format!("PerformMoreThreeCorner: DStr.Shape({iface2}) is not a face")
            })?;

            // `cxx:2587`.
            if is_same(&face1.0, &face2.0) {
                continue;
            }

            // `cxx:2589-2593`.
            if self.fproj.is_empty() {
                self.fproj.push(face1.0.clone());
                self.fproj.push(face2.0.clone());
            }
            self.moresurf[ic] = true;
            let nbface = self.fproj.len();
            self.nbface = nbface as i32;

            // `cxx:2596-2603`.
            if !is_same(&self.fproj[nbface - 1], &face2.0) {
                self.fproj.remove(nbface - 1);
                self.fproj.push(face2.0.clone());
            }
            if !is_same(&self.fproj[0], &face1.0) {
                self.fproj.remove(0);
                self.fproj.insert(0, face1.0.clone());
            }

            // `cxx:2605-2611`: the shared edges between consecutive faces and
            // the parameter of their first vertex.
            for nb in 1..nbface {
                let f1 = Face::wrap(self.fproj[nb - 1].clone()).ok_or_else(|| {
                    format!("PerformMoreThreeCorner: Fproj({nb}) is not a face")
                })?;
                let f2 = Face::wrap(self.fproj[nb].clone()).ok_or_else(|| {
                    format!("PerformMoreThreeCorner: Fproj({}) is not a face", nb + 1)
                })?;
                let edge = cherche_edge1(&f1, &f2)?;
                let first = edge_vertices(&edge)
                    .0
                    .ok_or_else(|| "PerformMoreThreeCorner: Ecom edge has no first vertex".to_string())?;
                let para = BRepTool::parameter_of_vertex(&edge.0, &first.0);
                self.ecom.push(edge);
                self.parcom.push(para);
            }

            // `cxx:2612-2634`: the bridge itself. `Fproj` is handed over as the
            // `Face` argument, so the projection runs over all of it, while
            // `Ecom` / `parcom` are indexed from their first entry.
            let faces: Vec<Option<Face>> = self.fproj.iter().cloned().map(Face::wrap).collect();
            let left_edge = self.edge_at(ic)?;
            let right_edge = self.edge_at(icplus)?;
            // The stripe handles are cloned out first: `curve_hermite` also
            // takes `&mut self.eproj` / `&mut self.parcom`, so the sides must
            // not borrow `self`.
            let left_stripe = self.stripe_at(ic)?.clone();
            let right_stripe = self.stripe_at(icplus)?.clone();
            let left = HermiteSide {
                stripe: &left_stripe,
                jfp: self.jf[ic],
                index: self.i[ic][icplus],
                par: self.p[ic][icplus],
                sens: self.sens[ic],
                sharp: self.sharp[ic],
                evive: &left_edge,
            };
            let right = HermiteSide {
                stripe: &right_stripe,
                // `CurveHermite` applies the `3 - jf` flip itself
                // (`cxx:336`), so the right side gets the raw index here.
                jfp: self.jf[icplus],
                index: self.i[icplus][ic],
                par: self.p[icplus][ic],
                sens: self.sens[icplus],
                sharp: self.sharp[icplus],
                evive: &right_edge,
            };

            let mut pr = Vec::new();
            let mut cr = Vec::new();
            let mut error = self.hermite_error;
            curve_hermite(
                ds,
                &left,
                &right,
                nbface,
                &self.ecom,
                &faces,
                &mut pr,
                &mut cr,
                &mut self.eproj,
                &mut self.parcom,
                &mut error,
            )?;
            self.hermite_error = error;

            // `cxx:2635-2636`: repeat the last edge and parameter, so the next
            // `ic` sees a pair whose second half closes the loop.
            let last_edge = self
                .ecom
                .get(nbface - 2)
                .cloned()
                .ok_or_else(|| "PerformMoreThreeCorner: Ecom is shorter than nbface".to_string())?;
            self.ecom.push(last_edge);
            let last_par = *self
                .parcom
                .get(nbface - 2)
                .ok_or_else(|| "PerformMoreThreeCorner: parcom is shorter than nbface".to_string())?;
            self.parcom.push(last_par);
        }

        Ok(())
    }
}
