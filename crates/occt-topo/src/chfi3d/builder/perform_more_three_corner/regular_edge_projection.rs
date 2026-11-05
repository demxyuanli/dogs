//! The regular-edge projection stage of `PerformMoreThreeCorner`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_CnCrn.cxx:2368-2541`
//! (OCCT 8.0.0).
//!
//! For every regular edge of the corner ring the stage walks to the first
//! non-regular edge, collects the faces that separate the two, builds the
//! Hermite bridge with [`super::curve_hermite::curve_hermite`] and writes the
//! returned parameters and projections back into the `p` / `i` / `isG1` grids.
//!
//! Translation notes:
//! - `proj2d1` / `proj2d2` / `cproj1` / `cproj2` are the `NCollection_Array1(0,
//!   size)` locals of `cxx:2376-2379`; only the slots the stage writes are
//!   filled, the rest stay `None` (`cxx:2506-2512` reads `Value(ind)` /
//!   `Value(ind + 1)`), and their consumers are the `GeomPlate` stages from
//!   `cxx:3063` on ([`super::plate_intermediary`]), so the grids live on
//!   [`MoreThreeCorner`],
//! - `Eproj` (`Epj` at the call site) is filled by `CurveHermite` and never
//!   read again (`cxx:2439`); the port keeps it on the state as well so the
//!   OCCT out-parameter is visible,
//! - `nbface`, `error`, `trouve` and `ii` are function-scope locals of
//!   `PerformMoreThreeCorner` in OCCT (`cxx:2373-2374`, `:1309`), so they are
//!   carried on the state instead of being reset per `ic`,
//! - the walk to the first non-regular edge is the `while (!trouve)` loop of
//!   `cxx:2399-2412`: on a ring whose entries are all regular it never
//!   terminates, which is the source's own behaviour.

use std::sync::Arc;

use occt_geom::Curve;
use occt_geom2d::curve::Curve2d;

use super::curve_hermite::{curve_hermite, HermiteSide};
use super::surf_index::{indices, surf_index, SurfType};
use super::MoreThreeCorner;
use crate::shape::{Edge, Face};
use crate::topopebrep::ds::DataStructure;

impl MoreThreeCorner {
    /// `cxx:2368-2541`.
    pub(super) fn project_regular_edges(&mut self, ds: &DataStructure) -> Result<(), String> {
        let nedge = self.nedge as usize;
        let size = self.size;

        // `cxx:2373-2379`.
        self.nbface = 0;
        self.hermite_error = 0.0;
        self.proj2d1 = vec![None; size + 1];
        self.proj2d2 = vec![None; size + 1];
        self.cproj1 = vec![None; size + 1];
        self.cproj2 = vec![None; size + 1];

        if self.deuxconges {
            return Ok(());
        }

        // `cxx:2383-2540`. The OCCT `for (ic = 0; ...)` is driven by `ic`, and
        // the body advances it by `nbface - 1` before the loop's own `ic++`.
        let mut ic = 0usize;
        while ic < nedge {
            if !self.regul[ic] {
                ic += 1;
                continue;
            }

            // `cxx:2393-2396`.
            let (icplus, icmoins) = indices(nedge, ic);
            let (icplus2, _) = indices(nedge, icplus);
            let _ = icplus2;

            // `cxx:2397-2412`: the first non-regular edge after `icplus`.
            let mut indfin = icplus;
            let mut trouve = false;
            let mut ii = icplus;
            while !trouve {
                if !self.regul[ii] {
                    indfin = ii;
                    trouve = true;
                }
                ii = if ii == nedge - 1 { 0 } else { ii + 1 };
            }

            // `cxx:2413`.
            let (indfinplus, indfinmoins) = indices(nedge, indfin);
            let _ = indfinplus;

            // `cxx:2414-2436`: the first face of the strip.
            let mut lface: Vec<Option<Face>> = Vec::new();
            if !self.sharp[icmoins] {
                let ilin = if self.jf[icmoins] == 1 {
                    surf_index(&self.cd, icmoins, self.i[icmoins][ic] as usize, SurfType::Face1)?
                } else {
                    surf_index(&self.cd, icmoins, self.i[icmoins][ic] as usize, SurfType::Face2)?
                };
                lface.push(self.face_at(ds, ilin)?);
            } else {
                lface.push(self.fvive_face(ic, icmoins)?);
            }

            // `cxx:2430-2437`.
            let nbface = if indfin > icmoins {
                indfin - icmoins
            } else {
                nedge - (icmoins - indfin)
            };
            self.nbface = nbface as i32;

            // `cxx:2439-2461`: the edges between the two stripes, with the
            // parameters of their live slots.
            let mut ledge: Vec<Edge> = Vec::new();
            let mut seqpr: Vec<f64> = Vec::new();
            let mut ii = ic;
            for nf in 1..nbface {
                let (iiplus, _iimoins) = indices(nedge, ii);
                ledge.push(self.edge_at(ii)?);
                seqpr.push(self.p[ii][iiplus]);
                // `if (nf != nbface - 1)` (`cxx:2449-2452`).
                if nf != nbface - 1 {
                    lface.push(self.fvive_face(ii, iiplus)?);
                }
                ii = if ii == nedge - 1 { 0 } else { ii + 1 };
            }

            // `cxx:2462-2477`: the last face of the strip.
            if !self.sharp[indfin] {
                let jfp = 3 - self.jf[indfin];
                let ilin = if jfp == 1 {
                    surf_index(
                        &self.cd,
                        indfin,
                        self.i[indfin][indfinmoins] as usize,
                        SurfType::Face1,
                    )?
                } else {
                    surf_index(
                        &self.cd,
                        indfin,
                        self.i[indfin][indfinmoins] as usize,
                        SurfType::Face2,
                    )?
                };
                lface.push(self.face_at(ds, ilin)?);
            } else {
                lface.push(self.fvive_face(indfin, indfinmoins)?);
            }

            // `cxx:2479-2500`.
            let left_edge = self.edge_at(icmoins)?;
            let right_edge = self.edge_at(indfin)?;
            let left = HermiteSide {
                stripe: self.stripe_at(icmoins)?,
                jfp: self.jf[icmoins],
                index: self.i[icmoins][ic],
                par: self.p[icmoins][ic],
                sens: self.sens[icmoins],
                sharp: self.sharp[icmoins],
                evive: &left_edge,
            };
            let right = HermiteSide {
                stripe: self.stripe_at(indfin)?,
                // `jfp = 3 - jf.Value(indfin)` (`cxx:2463`).
                jfp: 3 - self.jf[indfin],
                index: self.i[indfin][indfinmoins],
                par: self.p[indfin][indfinmoins],
                sens: self.sens[indfin],
                sharp: self.sharp[indfin],
                evive: &right_edge,
            };

            let mut pr: Vec<Option<Arc<dyn Curve2d>>> = Vec::new();
            let mut cr: Vec<Option<Arc<dyn Curve>>> = Vec::new();
            let mut epj: Vec<Option<Edge>> = Vec::new();
            let mut param = seqpr.clone();
            let mut error = 0.0;
            curve_hermite(
                ds,
                &left,
                &right,
                nbface,
                &ledge,
                &lface,
                &mut pr,
                &mut cr,
                &mut epj,
                &mut param,
                &mut error,
            )?;
            self.hermite_error = error;

            // `cxx:2501-2521`: write the parameters and the projections back
            // along the same ring walk.
            let mut ii = ic;
            for ind in 1..nbface {
                let (iiplus, iimoins) = indices(nedge, ii);
                let par = *param
                    .get(ind)
                    .ok_or_else(|| format!("PerformMoreThreeCorner: param({ind}) is missing"))?;
                self.p[ii][iiplus] = par;
                self.p[ii][iimoins] = par;
                // `proj2d1.SetValue(ii, pr.Value(ind))` - OCCT raises
                // `Standard_OutOfRange` when the projection appended fewer
                // entries than `nbface`.
                self.proj2d1[ii] = pr.get(ind).cloned().ok_or_else(|| {
                    format!("PerformMoreThreeCorner: pr({ind}) is missing")
                })?;
                self.proj2d2[ii] = pr.get(ind + 1).cloned().ok_or_else(|| {
                    format!("PerformMoreThreeCorner: pr({}) is missing", ind + 1)
                })?;
                self.cproj1[ii] = cr.get(ind).cloned().ok_or_else(|| {
                    format!("PerformMoreThreeCorner: cr({ind}) is missing")
                })?;
                self.cproj2[ii] = cr.get(ind + 1).cloned().ok_or_else(|| {
                    format!("PerformMoreThreeCorner: cr({}) is missing", ind + 1)
                })?;
                ii = if ii == nedge - 1 { 0 } else { ii + 1 };
            }

            // `cxx:2522-2534`: both ends non-sharp - the whole arc becomes G1.
            if !self.sharp[icmoins] && !self.sharp[indfin] {
                let mut ii = icmoins;
                while ii != indfin {
                    self.is_g1[ii] = true;
                    ii = if ii == nedge - 1 { 0 } else { ii + 1 };
                }
            }

            // `cxx:2538`: `ic = ic + nbface - 1;`, then the loop's `ic++`.
            ic += nbface - 1;
            ic += 1;
        }

        Ok(())
    }
}
