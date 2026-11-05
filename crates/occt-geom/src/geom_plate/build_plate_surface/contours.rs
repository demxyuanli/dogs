//! `Disc2dContour` (`cxx:871-1019`) and `Disc3dContour` (`cxx:1024-1162`):
//! sample the constraint points (and three quarters of every interval between
//! two consecutive constraint parameters) either in the initial surface's
//! parameter space or in 3D, as points or as normals.
//!
//! Translation notes:
//! - Both methods are unused inside OCCT itself: `GeomPlate_BuildPlateSurface`
//!   never calls them and no other class does either (they are kept for API
//!   parity with the header).
//! - `myProj.Initialize(...)` (`cxx:892`, `:1040`) only stores the projection
//!   window; the port's [`BuildPlateSurface::project_point`] rebuilds the
//!   engine with the same stored `(tol_u, tol_v)` on every call (see
//!   `projection.rs`).
//! - The `ProjectedCurve` / `Curve2dOnSurf` / `ProjectPoint` cascade repeated
//!   four times per interval is [`BuildPlateSurface::constraint_parameter`],
//!   the same one `EcartContraintesMil` uses.
//! - An empty `myParCont(i)` sequence makes OCCT raise on `Value(1)`; the port
//!   reports an error instead.

use occt_core::gp::{GpVec, GpXY, GpXyz};

use super::BuildPlateSurface;

impl BuildPlateSurface {
    /// `Disc2dContour(nbp, Seq2d)` (`cxx:871-1019`): the constraint points in
    /// the initial surface's parameter space.
    pub fn disc2d_contour(&self, seq2d: &mut Vec<GpXY>) -> Result<(), String> {
        seq2d.clear();
        // `cxx:893-904`: the point constraints.
        for i in 0..self.pnt_cont.len() {
            if self.pnt(i).order() != -1 {
                let p2d = self.pnt(i).pnt2d_on_surf();
                seq2d.push(GpXY::new(p2d.x(), p2d.y()));
            }
        }
        // `cxx:905-1015`: the curve constraints.
        for i in 0..self.lin_cont.len() {
            if self.lin(i).order() == -1 {
                continue;
            }
            let nb_pt = self.par_len(i);
            if nb_pt == 0 {
                return Err(
                    "GeomPlate_BuildPlateSurface::Disc2dContour: empty parameter table".to_string(),
                );
            }
            // `cxx:910-930`: the first constraint point (`j = 0`).
            let (pu, pv) = self.first_parameter_2d(i)?;
            seq2d.push(GpXY::new(pu, pv));
            // `cxx:931-1014`.
            for j in 1..nb_pt - 1 {
                let uj = self.par(i)[j];
                let ujp1 = self.par(i)[j + 1];
                for u in [
                    (ujp1 + 3.0 * uj) / 4.0,
                    (ujp1 + uj) / 2.0,
                    (3.0 * ujp1 + uj) / 4.0,
                    ujp1,
                ] {
                    let (pu, pv) = self.parameter_2d_at(i, u)?;
                    seq2d.push(GpXY::new(pu, pv));
                }
            }
        }
        Ok(())
    }

    /// `Disc3dContour(nbp, iordre, Seq3d)` (`cxx:1024-1162`): the constraint
    /// points in 3D, or their normals when `iordre == 1`.
    pub fn disc3d_contour(&self, iordre: i32, seq3d: &mut Vec<GpXyz>) -> Result<(), String> {
        seq3d.clear();
        // `cxx:1050-1070`: the point constraints.
        for i in 0..self.pnt_cont.len() {
            if self.pnt(i).order() == -1 {
                continue;
            }
            if iordre == 0 {
                let p3d = self.pnt(i).d0();
                seq3d.push(*p3d.xyz());
            } else {
                let (_, v1h, v2h) = self.pnt(i).d1();
                seq3d.push(*v1h.crossed(&v2h).xyz());
            }
        }
        // `cxx:1073-1161`: the curve constraints.
        for i in 0..self.lin_cont.len() {
            if self.lin(i).order() == -1 {
                continue;
            }
            let nb_pt = self.par_len(i);
            if nb_pt == 0 {
                return Err(
                    "GeomPlate_BuildPlateSurface::Disc3dContour: empty parameter table".to_string(),
                );
            }
            // `cxx:1080-1097`: the first constraint point.
            seq3d.push(self.curve_value_3d(i, self.par(i)[0], iordre)?);
            // `cxx:1099-1159`.
            for j in 1..nb_pt - 1 {
                let uj = self.par(i)[j];
                let ujp1 = self.par(i)[j + 1];
                for u in [
                    (ujp1 + 3.0 * uj) / 4.0,
                    (ujp1 + uj) / 2.0,
                    (3.0 * ujp1 + uj) / 4.0,
                    ujp1,
                ] {
                    seq3d.push(self.curve_value_3d(i, u, iordre)?);
                }
            }
        }
        Ok(())
    }

    /// `cxx:1031-1034`: `iordre` is only ever 0 or 1; the `OCCT_DEBUG` message
    /// has no return value in OCCT either.
    fn curve_value_3d(&self, i: usize, u: f64, iordre: i32) -> Result<GpXyz, String> {
        if iordre == 0 {
            Ok(*self.lin(i).d0(u).xyz())
        } else {
            let (_, v1h, v2h): (occt_core::gp::GpPnt, GpVec, GpVec) = self.lin(i).d1(u)?;
            Ok(*v1h.crossed(&v2h).xyz())
        }
    }

    /// `cxx:910-927`: the 2D point at the first constraint parameter.
    fn first_parameter_2d(&self, i: usize) -> Result<(f64, f64), String> {
        self.parameter_2d_at(i, self.par(i)[0])
    }

    /// The `ProjectedCurve` / `Curve2dOnSurf` / `ProjectPoint` cascade of
    /// `cxx:911-927`, `:936-949`, `:954-967`, `:974-989` and `:996-1010`.
    fn parameter_2d_at(&self, i: usize, u: f64) -> Result<(f64, f64), String> {
        let cc = self.lin(i);
        if let Some(projected) = cc.projected_curve() {
            let p2d = projected.d0(u);
            return Ok((p2d.x(), p2d.y()));
        }
        if let Some(curve2d) = cc.curve2d_on_surf()? {
            let p2d = curve2d.d0(u);
            return Ok((p2d.x(), p2d.y()));
        }
        let p3d = cc.d0(u);
        let p2d = self.project_point(&p3d)?;
        Ok((p2d.x(), p2d.y()))
    }
}
