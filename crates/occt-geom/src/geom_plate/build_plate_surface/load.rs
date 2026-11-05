//! `LoadCurve` / `LoadPoint` (`GeomPlate_BuildPlateSurface.cxx:2404-2584`):
//! turns the sampled constraint points stored in `myPlateCont` into `Plate`
//! point and G1 / G2 constraints.
//!
//! Translation notes:
//! - The `CC->D0/D1/D2` evaluations and the `mySurfInit->D0/D1/D2` evaluations
//!   are gathered before any `myPlate.Load` call, because the port owns the
//!   constraints by value (a `&self` borrow of `lin_cont` would conflict with the
//!   `&mut self.plate` the loads need). The per-point order
//!   (G0, then G1, then G2) is unchanged.
//! - `Plate_FreeGtoCConstraint(P2d, D1init, D1final)` (`cxx:2462`) takes the
//!   defaults `IncrementalLoad = 1.0` / `orientation = 0`
//!   (`Plate_FreeGtoCConstraint.hxx:37-41`).
//! - The G2 arm always loads `Plate_GtoCConstraint`; the `FreeGtoCConstraint`
//!   alternative is commented out in OCCT (`cxx:2500-2510`, `:2563-2581`).

use occt_core::gp::{GpPnt2d, GpVec, GpXY, GpXyz};

use crate::plate::{FreeGtoCConstraint, GtoCConstraint, PinpointConstraint, PlateD1, PlateD2};

use super::BuildPlateSurface;

/// The values of one constraint point that `LoadCurve` feeds to the plate.
struct CurveSample {
    /// `P2d.XY()`.
    uv: GpXY,
    /// `Pdif.XYZ()`.
    pdif: GpXyz,
    /// `(D1init, D1final)` when the G1 branch runs.
    g1: Option<(PlateD1, PlateD1)>,
    /// `(D1init, D1final, D2init, D2final)` when the G2 branch runs.
    g2: Option<(PlateD1, PlateD1, PlateD2, PlateD2)>,
}

/// The values of one point constraint that `LoadPoint` feeds to the plate.
struct PointSample {
    uv: GpXY,
    pdif: GpXyz,
    g1: Option<(PlateD1, PlateD1)>,
    g2: Option<(PlateD1, PlateD1, PlateD2, PlateD2)>,
}

impl BuildPlateSurface {
    /// `myPlateCont->Value(i).Value(j)` mapped to the plate parameter
    /// (`cxx:2427-2441`, `:2285-2289`, `:2327-2332`): `ProjectedCurve` first,
    /// then `Curve2dOnSurf`, then a projection of the 3D point.
    pub(super) fn plate_parameter(
        &self,
        i: usize,
        par: f64,
        p3d: &occt_core::gp::GpPnt,
    ) -> Result<GpPnt2d, String> {
        let cc = self.lin(i);
        if let Some(projected) = cc.projected_curve() {
            return Ok(projected.d0(par));
        }
        if let Some(curve2d) = cc.curve2d_on_surf()? {
            return Ok(curve2d.d0(par));
        }
        self.project_point(p3d)
    }

    /// `LoadCurve(NbBoucle, OrderMax)` (`cxx:2406-2516`).
    pub(super) fn load_curve(&mut self, nb_boucle: i32, order_max: i32) -> Result<(), String> {
        let n = self.lin_cont.len();
        let mut samples: Vec<(i32, Vec<CurveSample>)> = Vec::with_capacity(n);
        for i in 0..n {
            let cc = self.lin(i);
            if cc.order() == -1 {
                continue;
            }
            // `cxx:2418-2419`.
            let tang = cc.order().min(order_max);
            if tang == -1 {
                continue;
            }
            let nt = self.plate_len(i);
            let mut per_curve = Vec::with_capacity(nt);
            for j in 0..nt {
                let par = self.plate_at(i)[j];
                // `cxx:2425-2441`: the 3D point and its parameter.
                let p3d = cc.d0(par);
                let p2d = self.plate_parameter(i, par, &p3d)?;
                let pp = self.surf_d0(&p2d)?;
                let pdif = GpXyz::new(
                    -pp.x() + p3d.x(),
                    -pp.y() + p3d.y(),
                    -pp.z() + p3d.z(),
                );
                let mut sample = CurveSample {
                    uv: GpXY::new(p2d.x(), p2d.y()),
                    pdif,
                    g1: None,
                    g2: None,
                };
                // `cxx:2449-2488`: G1.
                if tang == 1 {
                    let (_, v1, v2) = cc.d1(par)?;
                    let (_, v3, v4) = self.surf_d1(&p2d)?;
                    sample.g1 = Some((
                        PlateD1::new(*v3.xyz(), *v4.xyz()),
                        PlateD1::new(*v1.xyz(), *v2.xyz()),
                    ));
                }
                // `cxx:2491-2513`: G2.
                if tang == 2 {
                    let (_, v1, v2, v5, v6, v7) = cc.d2(par)?;
                    let (_, v3, v4, v8, v9, v10) = self.surf_d2(&p2d)?;
                    sample.g2 = Some((
                        PlateD1::new(*v3.xyz(), *v4.xyz()),
                        PlateD1::new(*v1.xyz(), *v2.xyz()),
                        PlateD2::new(*v8.xyz(), *v9.xyz(), *v10.xyz()),
                        PlateD2::new(*v5.xyz(), *v6.xyz(), *v7.xyz()),
                    ));
                }
                per_curve.push(sample);
            }
            samples.push((tang, per_curve));
        }

        // `cxx:2445-2512`: the loads, in the order OCCT issues them.
        for (tang, per_curve) in samples {
            for sample in per_curve {
                self.plate
                    .load_pinpoint(PinpointConstraint::new(sample.uv, sample.pdif, 0, 0));
                if tang == 1 {
                    let (d1_init, d1_final) = sample.g1.expect("G1 sample");
                    if !self.free {
                        self.plate
                            .load_gto_c(&GtoCConstraint::new_g1(sample.uv, &d1_init, &d1_final));
                    } else if nb_boucle == 1 {
                        let free =
                            FreeGtoCConstraint::new_g1(sample.uv, &d1_init, &d1_final, 1.0, 0)?;
                        self.plate.load_free_gto_c(&free);
                    } else {
                        // `cxx:2470-2486`.
                        let mut normal = GpVec::from_xyz(d1_final.du())
                            .crossed(&GpVec::from_xyz(d1_final.dv()));
                        let norm = normal.magnitude();
                        if norm > 1.0e-12 {
                            normal = normal.divided(norm);
                        }
                        let uv = sample.uv;
                        let der_plate_u = GpVec::from_xyz(&self.prev_plate.evaluate_derivative(uv, 1, 0));
                        let der_plate_v = GpVec::from_xyz(&self.prev_plate.evaluate_derivative(uv, 0, 1));
                        // `DU.SetLinearForm(-(V3 + DerPlateU).Dot(Normal), Normal, DerPlateU)`.
                        let a_u = GpVec::from_xyz(d1_init.du()).added(&der_plate_u).dot(&normal);
                        let du = GpVec::zero().set_linear_form_2(-a_u, &normal, 1.0, &der_plate_u);
                        let a_v = GpVec::from_xyz(d1_init.dv()).added(&der_plate_v).dot(&normal);
                        let dv = GpVec::zero().set_linear_form_2(-a_v, &normal, 1.0, &der_plate_v);
                        self.plate.load_pinpoint(PinpointConstraint::new(
                            uv,
                            *du.xyz(),
                            1,
                            0,
                        ));
                        self.plate.load_pinpoint(PinpointConstraint::new(
                            uv,
                            *dv.xyz(),
                            0,
                            1,
                        ));
                    }
                }
                if tang == 2 {
                    let (d1_init, d1_final, d2_init, d2_final) = sample.g2.expect("G2 sample");
                    self.plate.load_gto_c(&GtoCConstraint::new_g2(
                        sample.uv,
                        &d1_init,
                        &d1_final,
                        &d2_init,
                        &d2_final,
                    ));
                }
            }
        }
        Ok(())
    }

    /// `LoadPoint(NbBoucle, OrderMax)` (`cxx:2522-2584`); the `NbBoucle`
    /// argument is unnamed in OCCT (`cxx:2522`).
    pub(super) fn load_point(&mut self, order_max: i32) -> Result<(), String> {
        let n = self.pnt_cont.len();
        let mut samples: Vec<PointSample> = Vec::with_capacity(n);
        for i in 0..n {
            let pc = self.pnt(i);
            // `cxx:2534-2540`.
            let p3d = pc.d0();
            let p2d = pc.pnt2d_on_surf();
            let pp = self.surf_d0(&p2d)?;
            let pdif = GpXyz::new(-pp.x() + p3d.x(), -pp.y() + p3d.y(), -pp.z() + p3d.z());
            let mut sample = PointSample {
                uv: GpXY::new(p2d.x(), p2d.y()),
                pdif,
                g1: None,
                g2: None,
            };
            // `cxx:2541-2542`.
            let tang = pc.order().min(order_max);
            // `cxx:2543-2560`.
            if tang == 1 {
                let (_, v1, v2) = pc.d1();
                let (_, v3, v4) = self.surf_d1(&p2d)?;
                sample.g1 = Some((
                    PlateD1::new(*v3.xyz(), *v4.xyz()),
                    PlateD1::new(*v1.xyz(), *v2.xyz()),
                ));
            }
            // `cxx:2562-2582`.
            if tang == 2 {
                let (_, v1, v2, v5, v6, v7) = pc.d2();
                let (_, v3, v4, v8, v9, v10) = self.surf_d2(&p2d)?;
                sample.g2 = Some((
                    PlateD1::new(*v3.xyz(), *v4.xyz()),
                    PlateD1::new(*v1.xyz(), *v2.xyz()),
                    PlateD2::new(*v8.xyz(), *v9.xyz(), *v10.xyz()),
                    PlateD2::new(*v5.xyz(), *v6.xyz(), *v7.xyz()),
                ));
            }
            samples.push(sample);
        }

        // `cxx:2539-2541`: every point constraint is loaded at G0, then at G1
        // or G2.
        for sample in samples {
            self.plate
                .load_pinpoint(PinpointConstraint::new(sample.uv, sample.pdif, 0, 0));
            if let Some((d1_init, d1_final)) = sample.g1 {
                if !self.free {
                    self.plate
                        .load_gto_c(&GtoCConstraint::new_g1(sample.uv, &d1_init, &d1_final));
                } else {
                    let free =
                        FreeGtoCConstraint::new_g1(sample.uv, &d1_init, &d1_final, 1.0, 0)?;
                    self.plate.load_free_gto_c(&free);
                }
            }
            if let Some((d1_init, d1_final, d2_init, d2_final)) = sample.g2 {
                self.plate.load_gto_c(&GtoCConstraint::new_g2(
                    sample.uv,
                    &d1_init,
                    &d1_final,
                    &d2_init,
                    &d2_final,
                ));
            }
        }
        Ok(())
    }

    /// `mySurfInit->D0(U, V, P)`.
    pub(super) fn surf_d0(&self, p2d: &GpPnt2d) -> Result<occt_core::gp::GpPnt, String> {
        let surf = self
            .surf_init
            .as_ref()
            .ok_or_else(|| "GeomPlate_BuildPlateSurface: no initial surface".to_string())?;
        Ok(surf.d0(p2d.x(), p2d.y()))
    }

    /// `mySurfInit->D1(U, V, P, V1, V2)`.
    pub(super) fn surf_d1(
        &self,
        p2d: &GpPnt2d,
    ) -> Result<(occt_core::gp::GpPnt, GpVec, GpVec), String> {
        let surf = self
            .surf_init
            .as_ref()
            .ok_or_else(|| "GeomPlate_BuildPlateSurface: no initial surface".to_string())?;
        Ok(surf.d1(p2d.x(), p2d.y()))
    }

    /// `mySurfInit->D2(U, V, P, V1, V2, V3, V4, V5)`.
    #[allow(clippy::type_complexity)]
    pub(super) fn surf_d2(&self, p2d: &GpPnt2d) -> Result<(occt_core::gp::GpPnt, GpVec, GpVec, GpVec, GpVec, GpVec), String> {
        let surf = self
            .surf_init
            .as_ref()
            .ok_or_else(|| "GeomPlate_BuildPlateSurface: no initial surface".to_string())?;
        Ok(surf.d2(p2d.x(), p2d.y()))
    }
}
