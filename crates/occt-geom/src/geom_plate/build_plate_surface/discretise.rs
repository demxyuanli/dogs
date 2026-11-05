//! `Discretise` (`GeomPlate_BuildPlateSurface.cxx:2153-2357`): samples every
//! constraint curve with the cosine distribution of `myNbPtsOnCur`, books the
//! parameters in `myParCont`, and keeps in `myPlateCont` only the points that
//! are neither duplicated across an intersection interval (`PntInter`) nor
//! inside an incompatible zone (`PntG1G1`). The points inside an incompatible
//! zone are replaced by a G0 `Plate_PinpointConstraint` on the spot
//! (`cxx:2284-2296`, `:2324-2334`).
//!
//! Translation notes:
//! - `AC2d` is a `Geom2dAdaptor_Curve` loaded only when `Curve2dOnSurf()` is
//!   non-null (`cxx:2190-2195`), so on a null p-curve OCCT keeps whatever the
//!   previous iteration loaded and the G0 arm (`:2285`) evaluates that stale
//!   adaptor. The port keeps the same last-loaded p-curve in `ac2d` and reports
//!   an error only when the very first evaluation has none.
//! - `GCPnts_AbscissaPoint::Length(AC2d, Uinit, U)` is the three-argument
//!   overload (`GCPnts_AbscissaPoint.cxx:342-347`) ->
//!   [`occt_geom2d::curve_ops::curve2d_length_between`].
//! - The `myParCont` / `myPlateCont` sequences are built as local vectors and
//!   assigned at the end, and the G0 pinpoint constraints are collected in
//!   `pinpoints` and loaded at the end: the relative order of the `Load` calls
//!   is unchanged (they all precede `LoadCurve`), while the port avoids
//!   borrowing `self` mutably inside the `&self` evaluation loop.

use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::gp::{GpXY, GpXyz};
use occt_geom2d::Curve2d;
use occt_geom2d::curve_ops::curve2d_length_between;
use occt_geom2d::geom2d_int::{GeomAbsCurveType2d, geom2d_curve_type};

use crate::law::{LawFunction, LawInterpol};
use crate::plate::PinpointConstraint;

use super::{BuildPlateSurface, constraint_pnt};

impl BuildPlateSurface {
    /// `Discretise(PntInter, PntG1G1)` (`cxx:2156-2357`).
    pub(super) fn discretise(
        &mut self,
        pnt_inter: &[Vec<f64>],
        pnt_g1g1: &[Vec<f64>],
    ) -> Result<(), String> {
        let nt_lin_cont = self.lin_cont.len();
        let mut par_cont: Vec<Vec<f64>> = vec![Vec::new(); nt_lin_cont];
        let mut plate_cont: Vec<Vec<f64>> = vec![Vec::new(); nt_lin_cont];
        let mut pinpoints: Vec<PinpointConstraint> = Vec::new();
        // `AC2d` (`cxx:2165`): kept across iterations as OCCT keeps the loaded
        // adaptor.
        let mut ac2d: Option<Arc<dyn Curve2d>> = None;
        for i in 0..nt_lin_cont {
            // `cxx:2183-2188`.
            let uinit = self.lin(i).first_parameter();
            let ufinal = self.lin(i).last_parameter();
            let c2d = self.lin(i).curve2d_on_surf()?;
            let acr = c2d.is_some();
            if let Some(c) = &c2d {
                ac2d = Some(c.clone());
            }
            // `cxx:2189-2222`: a law close to the curvilinear abscissa.
            let mut acrlaw = LawInterpol::new();
            let mut length2d = 0.0;
            if acr {
                let c = c2d.as_ref().expect("ACR implies a p-curve");
                let nbint = 20usize;
                let mut table: Vec<(f64, f64)> = Vec::with_capacity(nbint + 1);
                table.push((uinit, 0.0));
                length2d = curve2d_length_between(c.as_ref(), uinit, ufinal);
                for ii in 2..=nbint {
                    let u = uinit
                        + (ufinal - uinit)
                            * ((1.0 - (((ii - 1) as f64) * PI / (nbint as f64)).cos()) / 2.0);
                    table.push((u, curve2d_length_between(c.as_ref(), uinit, u)));
                }
                table.push((ufinal, length2d));
                acrlaw.set(&table, false)?;
            }
            // `cxx:2224-2226`.
            let nb_pnt_i = self.lin(i).nb_points().max(0) as usize;
            let nb_pt_inter = pnt_inter[i].len();
            let nb_pt_g1g1 = pnt_g1g1[i].len();
            for j in 1..=nb_pnt_i {
                // `cxx:2241-2254`: cosine distribution on the 2D abscissa.
                let inter = if j == nb_pnt_i {
                    // `cxx:2243`: to avoid a bug on Sun.
                    ufinal
                } else if acr {
                    let cur_length = length2d
                        * (1.0 - ((((j - 1) as f64) * PI / ((nb_pnt_i - 1) as f64)).cos())) / 2.0;
                    acrlaw.value(cur_length)
                } else {
                    uinit
                        + (ufinal - uinit)
                            * ((1.0 - ((((j - 1) as f64) * PI / ((nb_pnt_i - 1) as f64)).cos()))
                                / 2.0)
                };
                // `cxx:2255`.
                par_cont[i].push(inter);
                // `cxx:2256-2347`.
                if nb_pt_inter != 0 {
                    // `cxx:2258-2310`.
                    let mut l = 0usize;
                    while l + 1 < nb_pt_inter {
                        let interval = &pnt_inter[i];
                        if inter > interval[l] && inter < interval[l + 1] {
                            // `cxx:2264`: leave the loop without storing.
                            break;
                        } else if l + 1 >= nb_pt_inter - 1 {
                            // `cxx:2270`: the whole table was parsed and the
                            // point belongs to no common interval.
                            self.g1g1_step(
                                i,
                                inter,
                                &pnt_g1g1[i],
                                &mut plate_cont,
                                &mut pinpoints,
                                ac2d.as_deref(),
                            )?;
                            break;
                        }
                        l += 2;
                    }
                } else if nb_pt_g1g1 != 0 {
                    // `cxx:2315-2347`.
                    self.g1g1_step(
                        i,
                        inter,
                        &pnt_g1g1[i],
                        &mut plate_cont,
                        &mut pinpoints,
                        ac2d.as_deref(),
                    )?;
                } else {
                    // `cxx:2347-2354`: `!mySurfInitIsGive` and the p-curve is
                    // not a circle, or an interior sample.
                    let not_circle = match &c2d {
                        Some(c) => geom2d_curve_type(c.as_ref()) != GeomAbsCurveType2d::Circle,
                        None => true,
                    };
                    if (!self.surf_init_is_give && not_circle) || (j > 1 && j < nb_pnt_i) {
                        plate_cont[i].push(inter);
                    }
                }
            }
        }
        // `cxx:2167-2168`: the tables are rebuilt from scratch.
        self.par_cont = par_cont;
        self.plate_cont = plate_cont;
        for pc in pinpoints {
            self.plate.load_pinpoint(pc);
        }
        Ok(())
    }

    /// The incompatible-zone loop shared by `cxx:2274-2306` and
    /// `cxx:2319-2346`: `PntG1G1` is a flat list of `(lo, hi)` pairs; a point
    /// inside one is replaced by a G0 pinpoint constraint, and a point past the
    /// last pair is stored.
    fn g1g1_step(
        &self,
        i: usize,
        inter: f64,
        g1g1: &[f64],
        plate_cont: &mut [Vec<f64>],
        pinpoints: &mut Vec<PinpointConstraint>,
        ac2d: Option<&dyn Curve2d>,
    ) -> Result<(), String> {
        let mut k = 0usize;
        while k + 1 < g1g1.len() {
            if inter > g1g1[k] && inter < g1g1[k + 1] {
                // `cxx:2284-2296`: add a G0 constraint point.
                let c2d = ac2d.ok_or_else(|| {
                    "GeomPlate_BuildPlateSurface::Discretise: no 2D curve".to_string()
                })?;
                let p2d = c2d.d0(inter);
                let p3d = constraint_pnt(&self.lin(i).d0(inter));
                let pp = self.surf_d0(&p2d)?;
                let pdif = GpXyz::new(-pp.x() + p3d.x(), -pp.y() + p3d.y(), -pp.z() + p3d.z());
                pinpoints.push(PinpointConstraint::new(*p2d.xy(), pdif, 0, 0));
                return Ok(());
            } else if k + 1 >= g1g1.len() - 1 {
                // `cxx:2300-2301`: the point belongs to no G1 interval.
                plate_cont[i].push(inter);
                return Ok(());
            }
            k += 2;
        }
        Ok(())
    }
}
