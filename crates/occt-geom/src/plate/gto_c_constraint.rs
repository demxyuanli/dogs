//! `Plate_GtoCConstraint` — a G1, G2 or G3 constraint on the plate at one
//! `(u, v)`.
//!
//! Source: `Plate_GtoCConstraint.hxx:24-90`, `Plate_GtoCConstraint.cxx:20-595`,
//! `Plate_GtoCConstraint.lxx`.
//!
//! Translation notes:
//! - The seven OCCT constructors (`cxx:29-39` copy, `:41-79` G1, `:81-127` G1
//!   with `nP`, `:129-212` G2, `:214-309` G2 with `nP`, `:311-450` G3,
//!   `:452-595` G3 with `nP`) share one body in OCCT: the only difference is
//!   how many derivative levels are filled and whether the constrained
//!   direction is `normaleS` or its projection `nSP` onto the plane
//!   orthogonal to `nP`. The port keeps that shape in
//!   [`GtoCConstraint::build`]; every early `return` of the original (`cxx:51`,
//!   `:61`, `:68`, `:98`, `:103`, `:110`, `:117`, ...) is a `return` from
//!   `build` before `nb_pp_constraints` is raised.
//! - The last G3 overload takes `nP` but never reads it (`cxx:460-461` passes
//!   an unnamed `const gp_XYZ&`), so it is `new_g3` with `n_p = None`.
//! - `math_Gauss(mat)` uses the default `MinPivot` (`math_Gauss.hxx:49`) and
//!   the two `gauss.Solve` calls of `cxx:190-200` / `:285-297` / `:377-382`
//!   reuse one solution vector, as in OCCT.

use occt_core::gp::{GpXY, GpXyz};
use occt_math::{MathGauss, MathMatrix, MathVector};

use crate::plate::{PlateD1, PlateD2, PlateD3, PinpointConstraint};

/// `NORMIN` (`cxx:26`).
const NORMIN: f64 = 1.0e-10;
/// `COSMIN` (`cxx:27`).
const COSMIN: f64 = 1.0e-2;

/// `Plate_GtoCConstraint` (`hxx:33-89`).
#[derive(Debug, Clone)]
pub struct GtoCConstraint {
    ppc: [PinpointConstraint; 9],
    d1_surf_init: PlateD1,
    pnt2d: GpXY,
    nb_pp_constraints: usize,
}

impl GtoCConstraint {
    /// `Plate_GtoCConstraint(point2d, D1S, D1T)` (`cxx:41-79`).
    pub fn new_g1(point2d: GpXY, d1s: &PlateD1, d1t: &PlateD1) -> Self {
        let mut out = Self::empty(point2d, d1s);
        out.build(d1s, d1t, None, None, None);
        out
    }

    /// `Plate_GtoCConstraint(point2d, D1S, D1T, nP)` (`cxx:81-127`).
    pub fn new_g1_with_normal(point2d: GpXY, d1s: &PlateD1, d1t: &PlateD1, n_p: &GpXyz) -> Self {
        let mut out = Self::empty(point2d, d1s);
        out.build(d1s, d1t, None, None, Some(n_p));
        out
    }

    /// `Plate_GtoCConstraint(point2d, D1S, D1T, D2S, D2T)` (`cxx:129-212`).
    pub fn new_g2(
        point2d: GpXY,
        d1s: &PlateD1,
        d1t: &PlateD1,
        d2s: &PlateD2,
        d2t: &PlateD2,
    ) -> Self {
        let mut out = Self::empty(point2d, d1s);
        out.build(d1s, d1t, Some((d2s, d2t)), None, None);
        out
    }

    /// `Plate_GtoCConstraint(point2d, D1S, D1T, D2S, D2T, nP)` (`cxx:214-309`).
    pub fn new_g2_with_normal(
        point2d: GpXY,
        d1s: &PlateD1,
        d1t: &PlateD1,
        d2s: &PlateD2,
        d2t: &PlateD2,
        n_p: &GpXyz,
    ) -> Self {
        let mut out = Self::empty(point2d, d1s);
        out.build(d1s, d1t, Some((d2s, d2t)), None, Some(n_p));
        out
    }

    /// `Plate_GtoCConstraint(point2d, D1S, D1T, D2S, D2T, D3S, D3T)`
    /// (`cxx:311-450`).
    pub fn new_g3(
        point2d: GpXY,
        d1s: &PlateD1,
        d1t: &PlateD1,
        d2s: &PlateD2,
        d2t: &PlateD2,
        d3s: &PlateD3,
        d3t: &PlateD3,
    ) -> Self {
        let mut out = Self::empty(point2d, d1s);
        out.build(d1s, d1t, Some((d2s, d2t)), Some((d3s, d3t)), None);
        out
    }

    /// `Plate_GtoCConstraint(..., D3S, D3T, nP)` (`cxx:452-595`): `nP` is
    /// accepted but unused, exactly as in OCCT.
    #[allow(clippy::too_many_arguments)]
    pub fn new_g3_with_normal(
        point2d: GpXY,
        d1s: &PlateD1,
        d1t: &PlateD1,
        d2s: &PlateD2,
        d2t: &PlateD2,
        d3s: &PlateD3,
        d3t: &PlateD3,
        _n_p: &GpXyz,
    ) -> Self {
        let mut out = Self::empty(point2d, d1s);
        out.build(d1s, d1t, Some((d2s, d2t)), Some((d3s, d3t)), None);
        out
    }

    /// `Plate_GtoCConstraint(const Plate_GtoCConstraint& ref)` (`cxx:29-39`):
    /// `nb_PPConstraints` slots of `myPPC` are copied, the rest stay default.
    pub fn copy_of(&self) -> Self {
        let mut ppc = [PinpointConstraint::default(); 9];
        for i in 0..self.nb_pp_constraints {
            ppc[i] = self.ppc[i];
        }
        Self {
            ppc,
            d1_surf_init: self.d1_surf_init,
            pnt2d: self.pnt2d,
            nb_pp_constraints: self.nb_pp_constraints,
        }
    }

    /// `nb_PPC()` (`lxx:20`).
    pub fn nb_ppc(&self) -> usize {
        self.nb_pp_constraints
    }

    /// `GetPPC(Index)` (`lxx:25`), 0-based.
    pub fn get_ppc(&self, index: usize) -> &PinpointConstraint {
        &self.ppc[index]
    }

    /// `D1SurfInit()` (`lxx:30`).
    pub fn d1_surf_init(&self) -> &PlateD1 {
        &self.d1_surf_init
    }

    /// The `myD1SurfInit(D1S); pnt2d = point2d; nb_PPConstraints = 0;` prefix of
    /// every constructor.
    fn empty(point2d: GpXY, d1s: &PlateD1) -> Self {
        Self {
            ppc: [PinpointConstraint::default(); 9],
            d1_surf_init: *d1s,
            pnt2d: point2d,
            nb_pp_constraints: 0,
        }
    }

    /// The shared body of the six non-copy constructors.
    ///
    /// `d2` / `d3` select how many derivative levels are filled; `n_p`
    /// selects the `nSP` direction instead of `normaleS` (`cxx:109-115`,
    /// `:245-250`).
    #[allow(clippy::too_many_arguments)]
    fn build(
        &mut self,
        d1s: &PlateD1,
        d1t: &PlateD1,
        d2: Option<(&PlateD2, &PlateD2)>,
        d3: Option<(&PlateD3, &PlateD3)>,
        n_p: Option<&GpXyz>,
    ) {
        // `gp_XYZ normale = D1T.Du ^ D1T.Dv;` (`cxx:48`).
        let mut normale = d1t.du().crossed(d1t.dv());
        if normale.modulus() < NORMIN {
            return;
        }
        normale.normalize();

        // G1 Constraints: `normaleS = D1S.Du ^ D1S.Dv` (`cxx:58`).
        let mut normale_s = d1s.du().crossed(d1s.dv());
        if normale_s.modulus() < NORMIN {
            return;
        }
        normale_s.normalize();

        // `nSP = normaleS - nP * (nP * normaleS)` (`cxx:109-115`).
        let mut direction = normale_s;
        if let Some(n_p) = n_p {
            let mut n_sp = normale_s.subtracted(&n_p.multiplied(n_p.dot(&normale_s)));
            if n_sp.modulus() < NORMIN {
                return;
            }
            n_sp.normalize();
            direction = n_sp;
        }

        let cos_normales = normale.dot(&direction);
        if cos_normales.abs() < COSMIN {
            return;
        }
        let invcos = 1.0 / cos_normales;

        let du = direction.multiplied(-normale.dot(d1s.du()) * invcos);
        let dv = direction.multiplied(-normale.dot(d1s.dv()) * invcos);

        self.ppc[0] = PinpointConstraint::new(self.pnt2d, du, 1, 0);
        self.ppc[1] = PinpointConstraint::new(self.pnt2d, dv, 0, 1);
        self.nb_pp_constraints = 2;

        let Some((d2s, d2t)) = d2 else {
            return;
        };

        // G2 Constraints (`cxx:171-211`).
        let su = d1s.du().added(&du);
        let sv = d1s.dv().added(&dv);

        let mut mat = MathMatrix::new(0, 1, 0, 1);
        mat.set_value(0, 0, su.dot(d1t.du()));
        mat.set_value(0, 1, su.dot(d1t.dv()));
        mat.set_value(1, 0, sv.dot(d1t.du()));
        mat.set_value(1, 1, sv.dot(d1t.dv()));
        let gauss = MathGauss::new(&mat);
        if !gauss.is_done() {
            return;
        }

        let mut vec = MathVector::new(0, 1);
        vec.set_value(0, su.dot(&su));
        vec.set_value(1, su.dot(&sv));
        let mut sol = MathVector::new(0, 1);

        if gauss.solve(&vec, &mut sol).is_err() {
            return;
        }
        let a = sol.value(0);
        let b = sol.value(1);

        vec.set_value(0, sv.dot(&su));
        vec.set_value(1, sv.dot(&sv));

        if gauss.solve(&vec, &mut sol).is_err() {
            return;
        }
        let c = sol.value(0);
        let d = sol.value(1);

        let suu = d2t
            .duu()
            .multiplied(a * a)
            .added(&d2t.duv().multiplied(2.0 * a * b))
            .added(&d2t.dvv().multiplied(b * b));
        let suv = d2t
            .duu()
            .multiplied(a * c)
            .added(&d2t.duv().multiplied(a * d + b * c))
            .added(&d2t.dvv().multiplied(b * d));
        let svv = d2t
            .duu()
            .multiplied(c * c)
            .added(&d2t.duv().multiplied(2.0 * c * d))
            .added(&d2t.dvv().multiplied(d * d));

        let duu = direction.multiplied(normale.dot(&suu.subtracted(d2s.duu())) * invcos);
        let duv = direction.multiplied(normale.dot(&suv.subtracted(d2s.duv())) * invcos);
        let dvv = direction.multiplied(normale.dot(&svv.subtracted(d2s.dvv())) * invcos);

        self.ppc[2] = PinpointConstraint::new(self.pnt2d, duu, 2, 0);
        self.ppc[3] = PinpointConstraint::new(self.pnt2d, duv, 1, 1);
        self.ppc[4] = PinpointConstraint::new(self.pnt2d, dvv, 0, 2);
        self.nb_pp_constraints = 5;

        let Some((d3s, d3t)) = d3 else {
            return;
        };

        // G3 Constraints (`cxx:397-449`).
        vec.set_value(0, d2s.duu().added(&duu).subtracted(&suu).dot(&su));
        vec.set_value(1, d2s.duu().added(&duu).subtracted(&suu).dot(&sv));
        if gauss.solve(&vec, &mut sol).is_err() {
            return;
        }
        let b0uu = sol.value(0);
        let b1uu = sol.value(1);

        vec.set_value(0, d2s.duv().added(&duv).subtracted(&suv).dot(&su));
        vec.set_value(1, d2s.duv().added(&duv).subtracted(&suv).dot(&sv));
        if gauss.solve(&vec, &mut sol).is_err() {
            return;
        }
        let b0uv = sol.value(0);
        let b1uv = sol.value(1);

        vec.set_value(0, d2s.dvv().added(&dvv).subtracted(&svv).dot(&su));
        vec.set_value(1, d2s.dvv().added(&dvv).subtracted(&svv).dot(&sv));
        if gauss.solve(&vec, &mut sol).is_err() {
            return;
        }
        let b0vv = sol.value(0);
        let b1vv = sol.value(1);

        let mut suuu = d3t
            .duuu()
            .multiplied(a * a * a)
            .added(&d3t.duuv().multiplied(3.0 * a * a * b))
            .added(&d3t.duvv().multiplied(3.0 * a * b * b))
            .added(&d3t.dvvv().multiplied(b * b * b));
        let mut suuv = d3t
            .duuu()
            .multiplied(a * a * c)
            .added(&d3t.duuv().multiplied(a * a * d + 2.0 * a * b * c))
            .added(&d3t.duvv().multiplied(b * b * c + 2.0 * a * b * d))
            .added(&d3t.dvvv().multiplied(b * b * d));
        let mut suvv = d3t
            .duuu()
            .multiplied(a * c * c)
            .added(&d3t.duuv().multiplied(b * c * c + 2.0 * a * c * d))
            .added(&d3t.duvv().multiplied(a * d * d + 2.0 * b * c * d))
            .added(&d3t.dvvv().multiplied(b * d * d));
        let mut svvv = d3t
            .duuu()
            .multiplied(c * c * c)
            .added(&d3t.duuv().multiplied(3.0 * c * c * d))
            .added(&d3t.duvv().multiplied(3.0 * c * d * d))
            .added(&d3t.dvvv().multiplied(d * d * d));

        let a0u = a;
        let a1u = b;
        let a0v = c;
        let a1v = d;
        suuu = suuu
            .added(&d2t.duu().multiplied(3.0 * a0u * b0uu))
            .added(&d2t.duv().multiplied(3.0 * (a0u * b1uu + a1u * b0uu)))
            .added(&d2t.dvv().multiplied(3.0 * a1u * b1uu));
        suuv = suuv
            .added(&d2t.duu().multiplied(2.0 * a0u * b0uv + a0v * b0uu))
            .added(
                &d2t.duv()
                    .multiplied(2.0 * (a0u * b1uv + a1u * b0uv) + a0v * b1uu + a1v * b0uu),
            )
            .added(&d2t.dvv().multiplied(2.0 * a1u * b1uv + a1v * b1uu));
        suvv = suvv
            .added(&d2t.duu().multiplied(a0u * b0vv + 2.0 * a0v * b0uv))
            .added(
                &d2t.duv()
                    .multiplied(2.0 * (a0v * b1uv + a1v * b0uv) + a0u * b1vv + a1u * b0vv),
            )
            .added(&d2t.dvv().multiplied(2.0 * a1v * b1uv + a1u * b1vv));
        svvv = svvv
            .added(&d2t.duu().multiplied(3.0 * a0v * b0vv))
            .added(&d2t.duv().multiplied(3.0 * (a0v * b1vv + a1v * b0vv)))
            .added(&d2t.dvv().multiplied(3.0 * a1v * b1vv));

        let duuu = direction.multiplied(normale.dot(&suuu.subtracted(d3s.duuu())) * invcos);
        let duuv = direction.multiplied(normale.dot(&suuv.subtracted(d3s.duuv())) * invcos);
        let duvv = direction.multiplied(normale.dot(&suvv.subtracted(d3s.duvv())) * invcos);
        let dvvv = direction.multiplied(normale.dot(&svvv.subtracted(d3s.dvvv())) * invcos);

        self.ppc[5] = PinpointConstraint::new(self.pnt2d, duuu, 3, 0);
        self.ppc[6] = PinpointConstraint::new(self.pnt2d, duuv, 2, 1);
        self.ppc[7] = PinpointConstraint::new(self.pnt2d, duvv, 1, 2);
        self.ppc[8] = PinpointConstraint::new(self.pnt2d, dvvv, 0, 3);
        self.nb_pp_constraints = 9;
    }
}
