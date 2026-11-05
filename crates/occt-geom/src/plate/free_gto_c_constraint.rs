//! `Plate_FreeGtoCConstraint` — a weaker G1 / G2 / G3 constraint on the plate
//! that mixes pinpoint rows with `LinearScalarConstraint` rows along the
//! normal.
//!
//! Source: `Plate_FreeGtoCConstraint.hxx:24-79`,
//! `Plate_FreeGtoCConstraint.cxx:32-463`, `Plate_FreeGtoCConstraint.lxx`.
//!
//! Translation notes:
//! - The three OCCT constructors are kept as three functions; their control
//!   flow differs (`cxx:136-141` / `:284-294` / `:346-350` set two scalar
//!   rows and return when `normaleS` is null, while the G1 constructor has no
//!   such arm at all), so folding them into one function would hide a branch.
//! - The `IncrementalLoad != 1.` rotation block is identical in all three
//!   (`cxx:54-93`, `:144-182`, `:297-338`); it lives in
//!   [`incremental_rotation`]. `Ok(false)` is OCCT's silent `return`
//!   (counts left at zero), `Err` maps the `gp_Dir` / `gp_Trsf::SetRotation`
//!   construction errors that OCCT raises as `Standard_ConstructionError`.
//! - `gp_XYZ::Normalize()` raises on a null vector; the port returns `Err`
//!   from the constructors.

use occt_core::gp::{GpAx1, GpDir, GpPnt, GpTrsf, GpXY, GpXyz};
use occt_math::{MathGauss, MathMatrix, MathVector};

use crate::plate::{LinearScalarConstraint, PinpointConstraint, PlateD1, PlateD2, PlateD3};

/// `NORMIN` (`cxx:33`).
const NORMIN: f64 = 1.0e-10;
/// `COSMIN` (`cxx:34`).
const COSMIN: f64 = 1.0e-2;

/// `Plate_FreeGtoCConstraint` (`hxx:34-79`).
#[derive(Debug, Clone)]
pub struct FreeGtoCConstraint {
    pnt2d: GpXY,
    nb_pp_constraints: usize,
    nb_ls_constraints: usize,
    ppc: [PinpointConstraint; 5],
    lsc: [LinearScalarConstraint; 4],
}

impl FreeGtoCConstraint {
    /// `Plate_FreeGtoCConstraint(point2d, D1S, D1T, IncrementalLoad,
    /// orientation)` (`cxx:37-101`): the G1-only constructor.
    pub fn new_g1(
        point2d: GpXY,
        d1s: &PlateD1,
        d1t: &PlateD1,
        incremental_load: f64,
        orientation: i32,
    ) -> Result<Self, String> {
        let mut out = Self::empty(point2d);

        let mut normale = d1t.du().crossed(d1t.dv());
        if normale.modulus() < NORMIN {
            return Ok(out);
        }
        normale.normalize();

        if incremental_load != 1.0 {
            // `gp_XYZ N0 = D1S.Du ^ D1S.Dv; if (N0.Modulus() < NORMIN) return;`
            // (`cxx:56-60`).
            let mut n0 = d1s.du().crossed(d1s.dv());
            if n0.modulus() < NORMIN {
                return Ok(out);
            }
            n0.normalize();
            if !incremental_rotation(
                &n0,
                &mut normale,
                incremental_load,
                orientation,
                &mut [],
            )? {
                return Ok(out);
            }
        }

        let du = d1s.du().multiplied(-1.0);
        let dv = d1s.dv().multiplied(-1.0);

        out.lsc[0] = LinearScalarConstraint::from_ppc_coeff(
            PinpointConstraint::new(out.pnt2d, du, 1, 0),
            normale,
        );
        out.lsc[1] = LinearScalarConstraint::from_ppc_coeff(
            PinpointConstraint::new(out.pnt2d, dv, 0, 1),
            normale,
        );
        out.nb_ls_constraints = 2;
        Ok(out)
    }

    /// `Plate_FreeGtoCConstraint(point2d, D1S, D1T, D2S, D2T, IncrementalLoad,
    /// orientation)` (`cxx:105-251`): G1 + G2.
    #[allow(clippy::too_many_arguments)]
    pub fn new_g2(
        point2d: GpXY,
        d1s: &PlateD1,
        d1t: &PlateD1,
        d2s: &PlateD2,
        d2t: &PlateD2,
        incremental_load: f64,
        orientation: i32,
    ) -> Result<Self, String> {
        let mut out = Self::empty(point2d);

        // `Plate_D1 D1T = D1T0; Plate_D2 D2T = D2T0;` (`cxx:117-118`).
        let mut d1t_du = *d1t.du();
        let mut d1t_dv = *d1t.dv();
        let mut d2t_duu = *d2t.duu();
        let mut d2t_duv = *d2t.duv();
        let mut d2t_dvv = *d2t.dvv();

        let mut normale = d1t_du.crossed(&d1t_dv);
        if normale.modulus() < NORMIN {
            return Ok(out);
        }
        normale.normalize();

        // G1 Constraints: `normaleS = D1S.Du ^ D1S.Dv` (`cxx:127`).
        let mut normale_s = d1s.du().crossed(d1s.dv());
        if normale_s.modulus() < NORMIN {
            if incremental_load != 1.0 {
                return Ok(out);
            }
            let du = d1s.du().multiplied(-1.0);
            let dv = d1s.dv().multiplied(-1.0);

            out.lsc[0] = LinearScalarConstraint::from_ppc_coeff(
                PinpointConstraint::new(out.pnt2d, du, 1, 0),
                normale,
            );
            out.lsc[1] = LinearScalarConstraint::from_ppc_coeff(
                PinpointConstraint::new(out.pnt2d, dv, 0, 1),
                normale,
            );
            out.nb_ls_constraints = 2;
            return Ok(out);
        }
        normale_s.normalize();

        if incremental_load != 1.0 {
            // `gp_XYZ N0 = normaleS;` (`cxx:146`).
            let n0 = normale_s;
            if !incremental_rotation(
                &n0,
                &mut normale,
                incremental_load,
                orientation,
                &mut [
                    &mut d1t_du,
                    &mut d1t_dv,
                    &mut d2t_duu,
                    &mut d2t_duv,
                    &mut d2t_dvv,
                ],
            )? {
                return Ok(out);
            }
        }

        let cos_normales = normale.dot(&normale_s);
        if cos_normales.abs() < COSMIN {
            let du = d1s.du().multiplied(-1.0);
            let dv = d1s.dv().multiplied(-1.0);

            out.lsc[0] = LinearScalarConstraint::from_ppc_coeff(
                PinpointConstraint::new(out.pnt2d, du, 1, 0),
                normale,
            );
            out.lsc[1] = LinearScalarConstraint::from_ppc_coeff(
                PinpointConstraint::new(out.pnt2d, dv, 0, 1),
                normale,
            );
            out.nb_ls_constraints = 2;
            return Ok(out);
        }

        let invcos = 1.0 / cos_normales;

        let du = normale_s.multiplied(-normale.dot(d1s.du()) * invcos);
        let dv = normale_s.multiplied(-normale.dot(d1s.dv()) * invcos);

        out.ppc[0] = PinpointConstraint::new(out.pnt2d, du, 1, 0);
        out.ppc[1] = PinpointConstraint::new(out.pnt2d, dv, 0, 1);
        out.nb_pp_constraints = 2;

        // G2 Constraints (`cxx:207-250`).
        let su = d1s.du().added(&du);
        let sv = d1s.dv().added(&dv);

        let mut mat = MathMatrix::new(0, 1, 0, 1);
        mat.set_value(0, 0, su.dot(&d1t_du));
        mat.set_value(0, 1, su.dot(&d1t_dv));
        mat.set_value(1, 0, sv.dot(&d1t_du));
        mat.set_value(1, 1, sv.dot(&d1t_dv));
        let gauss = MathGauss::new(&mat);
        if !gauss.is_done() {
            return Ok(out);
        }

        let mut vec = MathVector::new(0, 1);
        vec.set_value(0, su.dot(&su));
        vec.set_value(1, su.dot(&sv));
        let mut sol = MathVector::new(0, 1);

        if gauss.solve(&vec, &mut sol).is_err() {
            return Ok(out);
        }
        let a = sol.value(0);
        let b = sol.value(1);

        vec.set_value(0, sv.dot(&su));
        vec.set_value(1, sv.dot(&sv));

        if gauss.solve(&vec, &mut sol).is_err() {
            return Ok(out);
        }
        let c = sol.value(0);
        let d = sol.value(1);

        let suu = d2t_duu
            .multiplied(a * a)
            .added(&d2t_duv.multiplied(2.0 * a * b))
            .added(&d2t_dvv.multiplied(b * b));
        let suv = d2t_duu
            .multiplied(a * c)
            .added(&d2t_duv.multiplied(a * d + b * c))
            .added(&d2t_dvv.multiplied(b * d));
        let svv = d2t_duu
            .multiplied(c * c)
            .added(&d2t_duv.multiplied(2.0 * c * d))
            .added(&d2t_dvv.multiplied(d * d));

        let duu = suu.subtracted(d2s.duu()).multiplied(incremental_load);
        let duv = suv.subtracted(d2s.duv()).multiplied(incremental_load);
        let dvv = svv.subtracted(d2s.dvv()).multiplied(incremental_load);

        out.lsc[0] = LinearScalarConstraint::from_ppc_coeff(
            PinpointConstraint::new(out.pnt2d, duu, 2, 0),
            normale,
        );
        out.lsc[1] = LinearScalarConstraint::from_ppc_coeff(
            PinpointConstraint::new(out.pnt2d, duv, 1, 1),
            normale,
        );
        out.lsc[2] = LinearScalarConstraint::from_ppc_coeff(
            PinpointConstraint::new(out.pnt2d, dvv, 0, 2),
            normale,
        );
        out.nb_ls_constraints = 3;
        Ok(out)
    }

    /// `Plate_FreeGtoCConstraint(point2d, D1S, D1T, D2S, D2T, D3S, D3T,
    /// IncrementalLoad, orientation)` (`cxx:255-463`): G1 + G2 + G3.
    #[allow(clippy::too_many_arguments)]
    pub fn new_g3(
        point2d: GpXY,
        d1s: &PlateD1,
        d1t: &PlateD1,
        d2s: &PlateD2,
        d2t: &PlateD2,
        d3s: &PlateD3,
        d3t: &PlateD3,
        incremental_load: f64,
        orientation: i32,
    ) -> Result<Self, String> {
        let mut out = Self::empty(point2d);

        // `Plate_D1 D1T = D1T0; Plate_D2 D2T = D2T0; Plate_D3 D3T = D3T0;`
        // (`cxx:269-270`).
        let mut d1t_du = *d1t.du();
        let mut d1t_dv = *d1t.dv();
        let mut d2t_duu = *d2t.duu();
        let mut d2t_duv = *d2t.duv();
        let mut d2t_dvv = *d2t.dvv();
        let mut d3t_duuu = *d3t.duuu();
        let mut d3t_duuv = *d3t.duuv();
        let mut d3t_duvv = *d3t.duvv();
        let mut d3t_dvvv = *d3t.dvvv();

        let mut normale = d1t_du.crossed(&d1t_dv);
        if normale.modulus() < NORMIN {
            return Ok(out);
        }
        normale.normalize();

        // G1 Constraints (`cxx:280-295`).
        let mut normale_s = d1s.du().crossed(d1s.dv());
        if normale_s.modulus() < NORMIN {
            if incremental_load != 1.0 {
                return Ok(out);
            }
            let du = d1s.du().multiplied(-1.0);
            let dv = d1s.dv().multiplied(-1.0);

            out.lsc[0] = LinearScalarConstraint::from_ppc_coeff(
                PinpointConstraint::new(out.pnt2d, du, 1, 0),
                normale,
            );
            out.lsc[1] = LinearScalarConstraint::from_ppc_coeff(
                PinpointConstraint::new(out.pnt2d, dv, 0, 1),
                normale,
            );
            out.nb_ls_constraints = 2;
            return Ok(out);
        }
        normale_s.normalize();

        if incremental_load != 1.0 {
            let n0 = normale_s;
            if !incremental_rotation(
                &n0,
                &mut normale,
                incremental_load,
                orientation,
                &mut [
                    &mut d1t_du,
                    &mut d1t_dv,
                    &mut d2t_duu,
                    &mut d2t_duv,
                    &mut d2t_dvv,
                    &mut d3t_duuu,
                    &mut d3t_duuv,
                    &mut d3t_duvv,
                    &mut d3t_dvvv,
                ],
            )? {
                return Ok(out);
            }
        }

        let cos_normales = normale.dot(&normale_s);
        if cos_normales.abs() < COSMIN {
            let du = d1s.du().multiplied(-1.0);
            let dv = d1s.dv().multiplied(-1.0);

            out.lsc[0] = LinearScalarConstraint::from_ppc_coeff(
                PinpointConstraint::new(out.pnt2d, du, 1, 0),
                normale,
            );
            out.lsc[1] = LinearScalarConstraint::from_ppc_coeff(
                PinpointConstraint::new(out.pnt2d, dv, 0, 1),
                normale,
            );
            out.nb_ls_constraints = 2;
            return Ok(out);
        }

        let invcos = 1.0 / cos_normales;

        let du = normale_s.multiplied(-normale.dot(d1s.du()) * invcos);
        let dv = normale_s.multiplied(-normale.dot(d1s.dv()) * invcos);

        out.ppc[0] = PinpointConstraint::new(out.pnt2d, du, 1, 0);
        out.ppc[1] = PinpointConstraint::new(out.pnt2d, dv, 0, 1);
        out.nb_pp_constraints = 2;

        // G2 Constraints (`cxx:361-404`).
        let su = d1s.du().added(&du);
        let sv = d1s.dv().added(&dv);

        let mut mat = MathMatrix::new(0, 1, 0, 1);
        mat.set_value(0, 0, su.dot(&d1t_du));
        mat.set_value(0, 1, su.dot(&d1t_dv));
        mat.set_value(1, 0, sv.dot(&d1t_du));
        mat.set_value(1, 1, sv.dot(&d1t_dv));
        let gauss = MathGauss::new(&mat);
        if !gauss.is_done() {
            return Ok(out);
        }

        let mut vec = MathVector::new(0, 1);
        vec.set_value(0, su.dot(&su));
        vec.set_value(1, su.dot(&sv));
        let mut sol = MathVector::new(0, 1);

        if gauss.solve(&vec, &mut sol).is_err() {
            return Ok(out);
        }
        let a = sol.value(0);
        let b = sol.value(1);

        vec.set_value(0, sv.dot(&su));
        vec.set_value(1, sv.dot(&sv));

        if gauss.solve(&vec, &mut sol).is_err() {
            return Ok(out);
        }
        let c = sol.value(0);
        let d = sol.value(1);

        let suu = d2t_duu
            .multiplied(a * a)
            .added(&d2t_duv.multiplied(2.0 * a * b))
            .added(&d2t_dvv.multiplied(b * b));
        let suv = d2t_duu
            .multiplied(a * c)
            .added(&d2t_duv.multiplied(a * d + b * c))
            .added(&d2t_dvv.multiplied(b * d));
        let svv = d2t_duu
            .multiplied(c * c)
            .added(&d2t_duv.multiplied(2.0 * c * d))
            .added(&d2t_dvv.multiplied(d * d));

        let duu = normale_s.multiplied(normale.dot(&suu.subtracted(d2s.duu())) * invcos);
        let duv = normale_s.multiplied(normale.dot(&suv.subtracted(d2s.duv())) * invcos);
        let dvv = normale_s.multiplied(normale.dot(&svv.subtracted(d2s.dvv())) * invcos);

        out.ppc[2] = PinpointConstraint::new(out.pnt2d, duu, 2, 0);
        out.ppc[3] = PinpointConstraint::new(out.pnt2d, duv, 1, 1);
        out.ppc[4] = PinpointConstraint::new(out.pnt2d, dvv, 0, 2);
        out.nb_pp_constraints = 5;

        // G3 Constraints (`cxx:406-462`).
        vec.set_value(0, d2s.duu().added(&duu).subtracted(&suu).dot(&su));
        vec.set_value(1, d2s.duu().added(&duu).subtracted(&suu).dot(&sv));
        if gauss.solve(&vec, &mut sol).is_err() {
            return Ok(out);
        }
        let b0uu = sol.value(0);
        let b1uu = sol.value(1);

        vec.set_value(0, d2s.duv().added(&duv).subtracted(&suv).dot(&su));
        vec.set_value(1, d2s.duv().added(&duv).subtracted(&suv).dot(&sv));
        if gauss.solve(&vec, &mut sol).is_err() {
            return Ok(out);
        }
        let b0uv = sol.value(0);
        let b1uv = sol.value(1);

        vec.set_value(0, d2s.dvv().added(&dvv).subtracted(&svv).dot(&su));
        vec.set_value(1, d2s.dvv().added(&dvv).subtracted(&svv).dot(&sv));
        if gauss.solve(&vec, &mut sol).is_err() {
            return Ok(out);
        }
        let b0vv = sol.value(0);
        let b1vv = sol.value(1);

        let mut suuu = d3t_duuu
            .multiplied(a * a * a)
            .added(&d3t_duuv.multiplied(3.0 * a * a * b))
            .added(&d3t_duvv.multiplied(3.0 * a * b * b))
            .added(&d3t_dvvv.multiplied(b * b * b));
        let mut suuv = d3t_duuu
            .multiplied(a * a * c)
            .added(&d3t_duuv.multiplied(a * a * d + 2.0 * a * b * c))
            .added(&d3t_duvv.multiplied(b * b * c + 2.0 * a * b * d))
            .added(&d3t_dvvv.multiplied(b * b * d));
        let mut suvv = d3t_duuu
            .multiplied(a * c * c)
            .added(&d3t_duuv.multiplied(b * c * c + 2.0 * a * c * d))
            .added(&d3t_duvv.multiplied(a * d * d + 2.0 * b * c * d))
            .added(&d3t_dvvv.multiplied(b * d * d));
        let mut svvv = d3t_duuu
            .multiplied(c * c * c)
            .added(&d3t_duuv.multiplied(3.0 * c * c * d))
            .added(&d3t_duvv.multiplied(3.0 * c * d * d))
            .added(&d3t_dvvv.multiplied(d * d * d));

        let a0u = a;
        let a1u = b;
        let a0v = c;
        let a1v = d;
        suuu = suuu
            .added(&d2t_duu.multiplied(3.0 * a0u * b0uu))
            .added(&d2t_duv.multiplied(3.0 * (a0u * b1uu + a1u * b0uu)))
            .added(&d2t_dvv.multiplied(3.0 * a1u * b1uu));
        suuv = suuv
            .added(&d2t_duu.multiplied(2.0 * a0u * b0uv + a0v * b0uu))
            .added(
                &d2t_duv
                    .multiplied(2.0 * (a0u * b1uv + a1u * b0uv) + a0v * b1uu + a1v * b0uu),
            )
            .added(&d2t_dvv.multiplied(2.0 * a1u * b1uv + a1v * b1uu));
        suvv = suvv
            .added(&d2t_duu.multiplied(a0u * b0vv + 2.0 * a0v * b0uv))
            .added(
                &d2t_duv
                    .multiplied(2.0 * (a0v * b1uv + a1v * b0uv) + a0u * b1vv + a1u * b0vv),
            )
            .added(&d2t_dvv.multiplied(2.0 * a1v * b1uv + a1u * b1vv));
        svvv = svvv
            .added(&d2t_duu.multiplied(3.0 * a0v * b0vv))
            .added(&d2t_duv.multiplied(3.0 * (a0v * b1vv + a1v * b0vv)))
            .added(&d2t_dvv.multiplied(3.0 * a1v * b1vv));

        let duuu = suuu.subtracted(d3s.duuu()).multiplied(incremental_load);
        let duuv = suuv.subtracted(d3s.duuv()).multiplied(incremental_load);
        let duvv = suvv.subtracted(d3s.duvv()).multiplied(incremental_load);
        let dvvv = svvv.subtracted(d3s.dvvv()).multiplied(incremental_load);

        out.lsc[0] = LinearScalarConstraint::from_ppc_coeff(
            PinpointConstraint::new(out.pnt2d, duuu, 3, 0),
            normale,
        );
        out.lsc[1] = LinearScalarConstraint::from_ppc_coeff(
            PinpointConstraint::new(out.pnt2d, duuv, 2, 1),
            normale,
        );
        out.lsc[2] = LinearScalarConstraint::from_ppc_coeff(
            PinpointConstraint::new(out.pnt2d, duvv, 1, 2),
            normale,
        );
        out.lsc[3] = LinearScalarConstraint::from_ppc_coeff(
            PinpointConstraint::new(out.pnt2d, dvvv, 0, 3),
            normale,
        );
        out.nb_ls_constraints = 4;
        Ok(out)
    }

    /// `nb_PPC()` (`lxx:20`).
    pub fn nb_ppc(&self) -> usize {
        self.nb_pp_constraints
    }

    /// `GetPPC(Index)` (`lxx:25`), 0-based.
    pub fn get_ppc(&self, index: usize) -> &PinpointConstraint {
        &self.ppc[index]
    }

    /// `nb_LSC()` (`lxx:35`).
    pub fn nb_lsc(&self) -> usize {
        self.nb_ls_constraints
    }

    /// `LSC(Index)` (`lxx:40`), 0-based.
    pub fn lsc(&self, index: usize) -> &LinearScalarConstraint {
        &self.lsc[index]
    }

    /// `pnt2d = point2d; nb_PPConstraints = 0; nb_LSConstraints = 0;`
    /// (`cxx:43-45`, `:113-115`, `:265-267`).
    fn empty(point2d: GpXY) -> Self {
        Self {
            pnt2d: point2d,
            nb_pp_constraints: 0,
            nb_ls_constraints: 0,
            ppc: [PinpointConstraint::default(); 5],
            lsc: [
                LinearScalarConstraint::new(),
                LinearScalarConstraint::new(),
                LinearScalarConstraint::new(),
                LinearScalarConstraint::new(),
            ],
        }
    }
}

/// The `IncrementalLoad != 1.` block shared by the three constructors
/// (`cxx:54-93`, `:144-182`, `:297-338`): rotate `normale` (and every
/// derivative vector in `vecs`) about `N0 x N1` by
/// `angle * (IncrementalLoad - 1.)`.
///
/// `n0` must already be normalized by the caller. `Ok(false)` is OCCT's bare
/// `return` for the `(s < 1.0e-2) && (c < 0.)` arm (`cxx:78-82`).
fn incremental_rotation(
    n0: &GpXyz,
    normale: &mut GpXyz,
    incremental_load: f64,
    orientation: i32,
    vecs: &mut [&mut GpXyz],
) -> Result<bool, String> {
    let mut n1 = *normale;
    if orientation != 0 {
        n1 = n1.multiplied(orientation as f64);
    }
    let mut c = n0.dot(&n1);
    if orientation == 0 && c < 0.0 {
        c *= -1.0;
        n1 = n1.multiplied(-1.0);
    }

    let s = n0.cross_magnitude(&n1);
    if s < 1.0e-2 && c < 0.0 {
        return Ok(false);
    }
    let angle = c.atan2(s);

    let mut d = n0.crossed(&n1);
    // `gp_XYZ::Normalize()` raises on a null vector (`gp_XYZ.hxx:598-604`).
    if !d.normalize() {
        return Err("Plate_FreeGtoCConstraint: gp_XYZ::Normalize on a null axis".to_string());
    }
    // `gp_Dir dir = gp_Dir(d);` (`cxx:87`).
    let dir = GpDir::from_xyz(&d).map_err(|e| e.to_string())?;
    let mut rota = GpTrsf::identity();
    let axe = GpAx1::new(GpPnt::new(0.0, 0.0, 0.0), dir);
    // `rota.SetRotation(Axe, angle * (IncrementalLoad - 1.));` (`cxx:90`).
    rota.set_rotation_ax1(&axe, angle * (incremental_load - 1.0))
        .map_err(|e| e.to_string())?;
    rota.transforms_xyz(normale);
    for v in vecs.iter_mut() {
        rota.transforms_xyz(v);
    }
    Ok(true)
}
