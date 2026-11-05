//! `BlendFunc_ConstRadInv` — the inverse function of `BlendFunc_ConstRad`: a
//! solution on a restriction (curve on surface) of one of the two supports.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_ConstRadInv.hxx:27-79`,
//! `BlendFunc_ConstRadInv.cxx:25-646` (OCCT 8.0.0).
//!
//! The vector `X` is `(t, w, U, V)`: the parameter on the curve on surface, the
//! parameter on the guide line, and the parameters of the partner point on the
//! other surface. The four equations are the same as `BlendFunc_ConstRad`
//! (`Value` / `Derivatives` / `Values`, `cxx:153-646`):
//! 1. the mid-point of the two extremities lies in the section plane,
//! 2. the two extremities are the ends of a circular arc of radius `|ray1|`
//!    centred on the guide-line point.
//!
//! Translation notes:
//! - `csurf` is an `occ::handle<Adaptor2d_Curve2d>`; the port stores it as
//!   `Option<Arc<dyn Curve2d>>` because `Set(OnFirst, C)` can be handed the
//!   empty handle. OCCT dereferences it unguarded in `Value` / `Derivatives` /
//!   `Values`; the port answers `false` instead of faulting. `GetTolerance`
//!   keeps the null case at `0.`.
//! - `Derivatives` computes `dnplan = (theD * nplan + d2cur) / |d1cur|`
//!   (`cxx:251-252`) while `Values` computes
//!   `dnplan = (-nplan . d2cur * nplan + d2cur) / |d1cur|` (`cxx:449-450`).
//!   The two differ; both are kept verbatim.
//! - `resul.SetLinearForm(ray1, ns1, -1., pts2.XYZ(), -ray2, ns2, pts1.XYZ())`
//!   (`cxx:227`) is `ray1*ns1 - ray2*ns2 + (pts1 - pts2)`; the port writes the
//!   same sum.
//! - `BlendFunc::ComputeNormal` (`cxx:187-211`, `:278-305`, `:484-510`) leaves
//!   the vector untouched when `CSLib` cannot define the normal; that is what
//!   `if let Some(..)` reproduces. A throwing `Surf->DN` is the port's panics
//!   with the message `try_eval_dn` reports, exactly where OCCT lets the
//!   exception escape.

use std::sync::Arc;

use occt_core::gp::{GpPnt2d, GpVec};
use occt_core::precision::Precision;
use occt_geom::approx_same_parameter;
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_geom2d::curve::Curve2d;
use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use crate::blend::{blend_func, BlendFuncInv};

use super::constants::EPS;
use super::pnt_vec;

/// `class BlendFunc_ConstRadInv : public Blend_FuncInv` (`hxx:27-79`).
pub struct ConstRadInv {
    /// `surf1` (`hxx:68`).
    surf1: Arc<dyn Surface>,
    /// `surf2` (`hxx:69`).
    surf2: Arc<dyn Surface>,
    /// `curv` (`hxx:70`): the guide line.
    curv: Arc<dyn Curve>,
    /// `csurf` (`hxx:71`): the curve on surface; `Set(OnFirst, C)` writes it.
    csurf: Option<Arc<dyn Curve2d>>,
    /// `ray1` (`hxx:72`).
    ray1: f64,
    /// `ray2` (`hxx:73`).
    ray2: f64,
    /// `choix` (`hxx:74`).
    choix: i32,
    /// `first` (`hxx:75`): `true` when `csurf` lies on `surf1`.
    first: bool,
}

impl ConstRadInv {
    /// `BlendFunc_ConstRadInv(S1, S2, C)` (`cxx:25-36`): `csurf` stays null,
    /// `ray1 = ray2 = 0.`, `choix = 0`, `first = false`.
    pub fn new(surf1: Arc<dyn Surface>, surf2: Arc<dyn Surface>, curv: Arc<dyn Curve>) -> Self {
        Self {
            surf1,
            surf2,
            curv,
            csurf: None,
            ray1: 0.0,
            ray2: 0.0,
            choix: 0,
            first: false,
        }
    }

    /// `Set(R, Choix)` (`cxx:40-75`): the sign of each radius follows the
    /// fillet quadrant.
    pub fn set_radius(&mut self, r: f64, choix: i32) {
        self.choix = choix;
        match choix {
            1 | 2 => {
                self.ray1 = -r;
                self.ray2 = -r;
            }
            3 | 4 => {
                self.ray1 = r;
                self.ray2 = -r;
            }
            5 | 6 => {
                self.ray1 = r;
                self.ray2 = r;
            }
            7 | 8 => {
                self.ray1 = -r;
                self.ray2 = r;
            }
            _ => {
                self.ray1 = -r;
                self.ray2 = -r;
            }
        }
    }

    /// `BlendFunc::ComputeNormal(Surf, P, N)`: the point `P` is the one of
    /// `csurf` when the normal is asked on the same support, `(U, V)` on the
    /// partner otherwise (`cxx:186-197`).
    ///
    /// Answers `None` where OCCT's `ComputeNormal` returns `false` — the caller
    /// then leaves the vector untouched.
    fn compute_normal_on(
        surf: &dyn Surface,
        on_curve_surface: bool,
        pt2d: &GpPnt2d,
        x: &MathVector,
        which: &str,
    ) -> Option<GpVec> {
        let p = if on_curve_surface {
            GpPnt2d::new(pt2d.x(), pt2d.y())
        } else {
            GpPnt2d::new(x.value(3), x.value(4))
        };
        blend_func::compute_normal(surf, &p)
            .unwrap_or_else(|e| panic!("BlendFunc_ConstRadInv: {which}->DN failed: {e}"))
    }
}

impl MathFunctionSetWithDerivatives for ConstRadInv {
    /// `Blend_FuncInv::NbVariables()` = `4` (`Blend_FuncInv.cxx:18-21`).
    fn nb_variables(&self) -> usize {
        4
    }

    /// `NbEquations()` (`cxx:78-81`): `4`.
    fn nb_equations(&self) -> usize {
        4
    }

    /// `Value(X, F)` (`cxx:153-232`). A missing `csurf` answers `false`.
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        let Some(csurf) = self.csurf.clone() else {
            return false;
        };
        let (ptcur, d1cur) = self.curv.d1(x.value(2));
        let nplan = d1cur.normalized();
        let the_d = -nplan.dot(&pnt_vec(&ptcur));

        let pt2d = csurf.value(x.value(1));

        let (pts1, d1u1, d1v1, pts2, d1u2, d1v2);
        if self.first {
            let r1 = self.surf1.d1(pt2d.x(), pt2d.y());
            pts1 = r1.0;
            d1u1 = r1.1;
            d1v1 = r1.2;
            let r2 = self.surf2.d1(x.value(3), x.value(4));
            pts2 = r2.0;
            d1u2 = r2.1;
            d1v2 = r2.2;
        } else {
            let r1 = self.surf1.d1(x.value(3), x.value(4));
            pts1 = r1.0;
            d1u1 = r1.1;
            d1v1 = r1.2;
            let r2 = self.surf2.d1(pt2d.x(), pt2d.y());
            pts2 = r2.0;
            d1u2 = r2.1;
            d1v2 = r2.2;
        }

        let sum = pnt_vec(&pts1).added(&pnt_vec(&pts2));
        f.set_value(1, nplan.dot(&sum) / 2.0 + the_d);

        let mut ns1 = d1u1.crossed(&d1v1);
        if ns1.magnitude() < EPS {
            if let Some(n) =
                Self::compute_normal_on(self.surf1.as_ref(), self.first, &pt2d, x, "surf1")
            {
                ns1 = n;
            }
        }

        let mut ns2 = d1u2.crossed(&d1v2);
        if ns2.magnitude() < EPS {
            if let Some(n) =
                Self::compute_normal_on(self.surf2.as_ref(), !self.first, &pt2d, x, "surf2")
            {
                ns2 = n;
            }
        }

        let mut norm1 = nplan.crossed(&ns1).magnitude();
        let mut norm2 = nplan.crossed(&ns2).magnitude();
        if norm1 < EPS {
            norm1 = 1.0;
        }
        if norm2 < EPS {
            // Unsatisfactory, but it is not necessary to stop (`cxx:220`).
            norm2 = 1.0;
        }

        ns1 = GpVec::zero().set_linear_form_2(nplan.dot(&ns1) / norm1, &nplan, -1.0 / norm1, &ns1);
        ns2 = GpVec::zero().set_linear_form_2(nplan.dot(&ns2) / norm2, &nplan, -1.0 / norm2, &ns2);
        // `resul = ray1*ns1 - 1*pts2 + (-ray2)*ns2 + pts1` (`cxx:227`).
        let resul = GpVec::zero()
            .set_linear_form_3(
                self.ray1,
                &ns1,
                -1.0,
                &pnt_vec(&pts2),
                -self.ray2,
                &ns2,
            )
            .added(&pnt_vec(&pts1));
        f.set_value(2, resul.x());
        f.set_value(3, resul.y());
        f.set_value(4, resul.z());

        true
    }

    /// `Derivatives(X, D)` (`cxx:236-433`).
    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        let Some(csurf) = self.csurf.clone() else {
            return false;
        };
        let (ptcur, d1cur, d2cur) = self.curv.d2(x.value(2));
        let normtgcur = d1cur.magnitude();
        let nplan = d1cur.normalized();
        let the_d = -nplan.dot(&pnt_vec(&ptcur));

        // `dnplan.SetLinearForm(theD, nplan, d2cur); dnplan /= normtgcur;`
        // (`cxx:251-252`). Note that `theD` is not `-nplan . d2cur` here; see
        // the module note.
        let dnplan = GpVec::zero()
            .set_linear_form_2(the_d, &nplan, 1.0, &d2cur)
            .divided(normtgcur);

        let (pt2d, v2d) = csurf.d1(x.value(1));

        let (pts1, d1u1, d1v1, d2u1, d2v1, d2uv1, pts2, d1u2, d1v2, d2u2, d2v2, d2uv2);
        if self.first {
            let (p, u1, v1, uu1, vv1, uv1) = self.surf1.d2(pt2d.x(), pt2d.y());
            pts1 = p;
            d1u1 = u1;
            d1v1 = v1;
            d2u1 = uu1;
            d2v1 = vv1;
            d2uv1 = uv1;
            let (p, u2, v2, uu2, vv2, uv2) = self.surf2.d2(x.value(3), x.value(4));
            pts2 = p;
            d1u2 = u2;
            d1v2 = v2;
            d2u2 = uu2;
            d2v2 = vv2;
            d2uv2 = uv2;
            let temp = GpVec::zero().set_linear_form_2(v2d.x(), &d1u1, v2d.y(), &d1v1);
            d.set_value(1, 1, nplan.dot(&temp) / 2.0);
            let temp = pnt_vec(&pts1).added(&pnt_vec(&pts2)).multiplied_scalar(0.5)
                .subtracted(&pnt_vec(&ptcur));
            d.set_value(1, 2, dnplan.dot(&temp) - normtgcur);
            d.set_value(1, 3, nplan.dot(&d1u2) / 2.0);
            d.set_value(1, 4, nplan.dot(&d1v2) / 2.0);
        } else {
            let (p, u1, v1, uu1, vv1, uv1) = self.surf1.d2(x.value(3), x.value(4));
            pts1 = p;
            d1u1 = u1;
            d1v1 = v1;
            d2u1 = uu1;
            d2v1 = vv1;
            d2uv1 = uv1;
            let (p, u2, v2, uu2, vv2, uv2) = self.surf2.d2(pt2d.x(), pt2d.y());
            pts2 = p;
            d1u2 = u2;
            d1v2 = v2;
            d2u2 = uu2;
            d2v2 = vv2;
            d2uv2 = uv2;
            let temp = GpVec::zero().set_linear_form_2(v2d.x(), &d1u2, v2d.y(), &d1v2);
            d.set_value(1, 1, nplan.dot(&temp) / 2.0);
            let temp = pnt_vec(&pts1).added(&pnt_vec(&pts2)).multiplied_scalar(0.5)
                .subtracted(&pnt_vec(&ptcur));
            d.set_value(1, 2, dnplan.dot(&temp) - normtgcur);
            d.set_value(1, 3, nplan.dot(&d1u1) / 2.0);
            d.set_value(1, 4, nplan.dot(&d1v1) / 2.0);
        }

        let mut ns1 = d1u1.crossed(&d1v1);
        if ns1.magnitude() < EPS {
            if let Some(n) =
                Self::compute_normal_on(self.surf1.as_ref(), self.first, &pt2d, x, "surf1")
            {
                ns1 = n;
            }
        }

        let mut ns2 = d1u2.crossed(&d1v2);
        if ns2.magnitude() < EPS {
            if let Some(n) =
                Self::compute_normal_on(self.surf2.as_ref(), !self.first, &pt2d, x, "surf2")
            {
                ns2 = n;
            }
        }

        let ncrossns1 = nplan.crossed(&ns1);
        let ncrossns2 = nplan.crossed(&ns2);
        let mut norm1 = ncrossns1.magnitude();
        let mut norm2 = ncrossns2.magnitude();
        if norm1 < EPS {
            norm1 = 1.0;
        }
        if norm2 < EPS {
            norm2 = 1.0;
        }

        let ndotns1 = nplan.dot(&ns1);
        let ndotns2 = nplan.dot(&ns2);

        // Derived compared to u1 (`cxx:330-340`).
        let temp = d2u1.crossed(&d1v1).added(&d1u1.crossed(&d2uv1));
        let grosterme = ncrossns1.dot(&nplan.crossed(&temp)) / norm1 / norm1;
        let resul1 = GpVec::zero()
            .set_linear_form_3(
                -self.ray1 / norm1 * (grosterme * ndotns1 - nplan.dot(&temp)),
                &nplan,
                self.ray1 * grosterme / norm1,
                &ns1,
                -self.ray1 / norm1,
                &temp,
            )
            .added(&d1u1);

        // Derived compared to v1 (`cxx:342-351`).
        let temp = d2uv1.crossed(&d1v1).added(&d1u1.crossed(&d2v1));
        let grosterme = ncrossns1.dot(&nplan.crossed(&temp)) / norm1 / norm1;
        let resul2 = GpVec::zero()
            .set_linear_form_3(
                -self.ray1 / norm1 * (grosterme * ndotns1 - nplan.dot(&temp)),
                &nplan,
                self.ray1 * grosterme / norm1,
                &ns1,
                -self.ray1 / norm1,
                &temp,
            )
            .added(&d1v1);

        if self.first {
            d.set_value(2, 1, v2d.x() * resul1.x() + v2d.y() * resul2.x());
            d.set_value(3, 1, v2d.x() * resul1.y() + v2d.y() * resul2.y());
            d.set_value(4, 1, v2d.x() * resul1.z() + v2d.y() * resul2.z());
        } else {
            d.set_value(2, 3, resul1.x());
            d.set_value(3, 3, resul1.y());
            d.set_value(4, 3, resul1.z());

            d.set_value(2, 4, resul2.x());
            d.set_value(3, 4, resul2.y());
            d.set_value(4, 4, resul2.z());
        }

        // Derived compared to w (parameter on guideline); it is assumed that
        // the radius is constant (`cxx:369-390`).
        let grosterme = ncrossns1.dot(&dnplan.crossed(&ns1)) / norm1 / norm1;
        let resul1 = GpVec::zero().set_linear_form_3(
            -self.ray1 / norm1 * (grosterme * ndotns1 - dnplan.dot(&ns1)),
            &nplan,
            self.ray1 * ndotns1 / norm1,
            &dnplan,
            self.ray1 * grosterme / norm1,
            &ns1,
        );

        let grosterme = ncrossns2.dot(&dnplan.crossed(&ns2)) / norm2 / norm2;
        let resul2 = GpVec::zero().set_linear_form_3(
            self.ray2 / norm2 * (grosterme * ndotns2 - dnplan.dot(&ns2)),
            &nplan,
            -self.ray2 * ndotns2 / norm2,
            &dnplan,
            -self.ray2 * grosterme / norm2,
            &ns2,
        );

        d.set_value(2, 2, resul1.x() + resul2.x());
        d.set_value(3, 2, resul1.y() + resul2.y());
        d.set_value(4, 2, resul1.z() + resul2.z());

        // Derived compared to u2 (`cxx:394-402`).
        let temp = d2u2.crossed(&d1v2).added(&d1u2.crossed(&d2uv2));
        let grosterme = ncrossns2.dot(&nplan.crossed(&temp)) / norm2 / norm2;
        let resul1 = GpVec::zero()
            .set_linear_form_3(
                self.ray2 / norm2 * (grosterme * ndotns2 - nplan.dot(&temp)),
                &nplan,
                -self.ray2 * grosterme / norm2,
                &ns2,
                self.ray2 / norm2,
                &temp,
            )
            .subtracted(&d1u2);

        // Derived compared to v2 (`cxx:404-412`).
        let temp = d2uv2.crossed(&d1v2).added(&d1u2.crossed(&d2v2));
        let grosterme = ncrossns2.dot(&nplan.crossed(&temp)) / norm2 / norm2;
        let resul2 = GpVec::zero()
            .set_linear_form_3(
                self.ray2 / norm2 * (grosterme * ndotns2 - nplan.dot(&temp)),
                &nplan,
                -self.ray2 * grosterme / norm2,
                &ns2,
                self.ray2 / norm2,
                &temp,
            )
            .subtracted(&d1v2);

        if !self.first {
            d.set_value(2, 1, v2d.x() * resul1.x() + v2d.y() * resul2.x());
            d.set_value(3, 1, v2d.x() * resul1.y() + v2d.y() * resul2.y());
            d.set_value(4, 1, v2d.x() * resul1.z() + v2d.y() * resul2.z());
        } else {
            d.set_value(2, 3, resul1.x());
            d.set_value(3, 3, resul1.y());
            d.set_value(4, 3, resul1.z());

            d.set_value(2, 4, resul2.x());
            d.set_value(3, 4, resul2.y());
            d.set_value(4, 4, resul2.z());
        }

        true
    }

    /// `Values(X, F, D)` (`cxx:436-646`): `Value` and `Derivatives` inlined the
    /// way OCCT inlines them, including the different `dnplan` (see the module
    /// note).
    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        let Some(csurf) = self.csurf.clone() else {
            return false;
        };
        let (ptcur, d1cur, d2cur) = self.curv.d2(x.value(2));
        let normtgcur = d1cur.magnitude();
        let nplan = d1cur.normalized();
        let the_d = -nplan.dot(&pnt_vec(&ptcur));

        // `dnplan.SetLinearForm(-nplan.Dot(d2cur), nplan, d2cur); dnplan /=
        // normtgcur;` (`cxx:449-450`).
        let dnplan = GpVec::zero()
            .set_linear_form_2(-nplan.dot(&d2cur), &nplan, 1.0, &d2cur)
            .divided(normtgcur);

        let (pt2d, v2d) = csurf.d1(x.value(1));

        let (pts1, d1u1, d1v1, d2u1, d2v1, d2uv1, pts2, d1u2, d1v2, d2u2, d2v2, d2uv2);
        if self.first {
            let (p, u1, v1, uu1, vv1, uv1) = self.surf1.d2(pt2d.x(), pt2d.y());
            pts1 = p;
            d1u1 = u1;
            d1v1 = v1;
            d2u1 = uu1;
            d2v1 = vv1;
            d2uv1 = uv1;
            let (p, u2, v2, uu2, vv2, uv2) = self.surf2.d2(x.value(3), x.value(4));
            pts2 = p;
            d1u2 = u2;
            d1v2 = v2;
            d2u2 = uu2;
            d2v2 = vv2;
            d2uv2 = uv2;

            let temp = GpVec::zero().set_linear_form_2(v2d.x(), &d1u1, v2d.y(), &d1v1);
            d.set_value(1, 1, nplan.dot(&temp) / 2.0);
            let temp = pnt_vec(&pts1).added(&pnt_vec(&pts2)).multiplied_scalar(0.5)
                .subtracted(&pnt_vec(&ptcur));
            d.set_value(1, 2, dnplan.dot(&temp) - normtgcur);
            d.set_value(1, 3, nplan.dot(&d1u2) / 2.0);
            d.set_value(1, 4, nplan.dot(&d1v2) / 2.0);
        } else {
            let (p, u1, v1, uu1, vv1, uv1) = self.surf1.d2(x.value(3), x.value(4));
            pts1 = p;
            d1u1 = u1;
            d1v1 = v1;
            d2u1 = uu1;
            d2v1 = vv1;
            d2uv1 = uv1;
            let (p, u2, v2, uu2, vv2, uv2) = self.surf2.d2(pt2d.x(), pt2d.y());
            pts2 = p;
            d1u2 = u2;
            d1v2 = v2;
            d2u2 = uu2;
            d2v2 = vv2;
            d2uv2 = uv2;

            let temp = GpVec::zero().set_linear_form_2(v2d.x(), &d1u2, v2d.y(), &d1v2);
            d.set_value(1, 1, nplan.dot(&temp) / 2.0);
            let temp = pnt_vec(&pts1).added(&pnt_vec(&pts2)).multiplied_scalar(0.5)
                .subtracted(&pnt_vec(&ptcur));
            d.set_value(1, 2, dnplan.dot(&temp) - normtgcur);
            d.set_value(1, 3, nplan.dot(&d1u1) / 2.0);
            d.set_value(1, 4, nplan.dot(&d1v1) / 2.0);
        }

        let sum = pnt_vec(&pts1).added(&pnt_vec(&pts2));
        f.set_value(1, nplan.dot(&sum) / 2.0 + the_d);

        let mut ns1 = d1u1.crossed(&d1v1);
        if ns1.magnitude() < EPS {
            if let Some(n) =
                Self::compute_normal_on(self.surf1.as_ref(), self.first, &pt2d, x, "surf1")
            {
                ns1 = n;
            }
        }

        let mut ns2 = d1u2.crossed(&d1v2);
        if ns2.magnitude() < EPS {
            if let Some(n) =
                Self::compute_normal_on(self.surf2.as_ref(), !self.first, &pt2d, x, "surf2")
            {
                ns2 = n;
            }
        }

        let ncrossns1 = nplan.crossed(&ns1);
        let ncrossns2 = nplan.crossed(&ns2);
        let mut norm1 = ncrossns1.magnitude();
        let mut norm2 = ncrossns2.magnitude();
        if norm1 < EPS {
            norm1 = 1.0;
        }
        if norm2 < EPS {
            norm2 = 1.0;
        }

        let ndotns1 = nplan.dot(&ns1);
        let ndotns2 = nplan.dot(&ns2);

        // `temp.SetLinearForm(ndotns1/norm1, nplan, -1/norm1, ns1);
        //  resul1.SetLinearForm(ray1, temp, gp_Vec(pts2, pts1));
        //  temp.SetLinearForm(ndotns2/norm2, nplan, -1/norm2, ns2);
        //  resul1.Subtract(ray2 * temp);` (`cxx:534-537`).
        let temp = GpVec::zero().set_linear_form_2(ndotns1 / norm1, &nplan, -1.0 / norm1, &ns1);
        let mut resul1 =
            GpVec::zero().set_linear_form_2(self.ray1, &temp, 1.0, &GpVec::from_pnts(&pts2, &pts1));
        let temp = GpVec::zero().set_linear_form_2(ndotns2 / norm2, &nplan, -1.0 / norm2, &ns2);
        resul1 = resul1.subtracted(&temp.multiplied_scalar(self.ray2));

        f.set_value(2, resul1.x());
        f.set_value(3, resul1.y());
        f.set_value(4, resul1.z());

        // Derived compared to u1 (`cxx:545-554`).
        let temp = d2u1.crossed(&d1v1).added(&d1u1.crossed(&d2uv1));
        let grosterme = ncrossns1.dot(&nplan.crossed(&temp)) / norm1 / norm1;
        let resul1 = GpVec::zero()
            .set_linear_form_3(
                -self.ray1 / norm1 * (grosterme * ndotns1 - nplan.dot(&temp)),
                &nplan,
                self.ray1 * grosterme / norm1,
                &ns1,
                -self.ray1 / norm1,
                &temp,
            )
            .added(&d1u1);

        // Derived compared to v1 (`cxx:556-565`).
        let temp = d2uv1.crossed(&d1v1).added(&d1u1.crossed(&d2v1));
        let grosterme = ncrossns1.dot(&nplan.crossed(&temp)) / norm1 / norm1;
        let resul2 = GpVec::zero()
            .set_linear_form_3(
                -self.ray1 / norm1 * (grosterme * ndotns1 - nplan.dot(&temp)),
                &nplan,
                self.ray1 * grosterme / norm1,
                &ns1,
                -self.ray1 / norm1,
                &temp,
            )
            .added(&d1v1);

        if self.first {
            d.set_value(2, 1, v2d.x() * resul1.x() + v2d.y() * resul2.x());
            d.set_value(3, 1, v2d.x() * resul1.y() + v2d.y() * resul2.y());
            d.set_value(4, 1, v2d.x() * resul1.z() + v2d.y() * resul2.z());
        } else {
            d.set_value(2, 3, resul1.x());
            d.set_value(3, 3, resul1.y());
            d.set_value(4, 3, resul1.z());

            d.set_value(2, 4, resul2.x());
            d.set_value(3, 4, resul2.y());
            d.set_value(4, 4, resul2.z());
        }

        // Derived compared to w (parameter on guideline); it is assumed that
        // the radius is constant (`cxx:585-604`).
        let grosterme = ncrossns1.dot(&dnplan.crossed(&ns1)) / norm1 / norm1;
        let resul1 = GpVec::zero().set_linear_form_3(
            -self.ray1 / norm1 * (grosterme * ndotns1 - dnplan.dot(&ns1)),
            &nplan,
            self.ray1 * ndotns1 / norm1,
            &dnplan,
            self.ray1 * grosterme / norm1,
            &ns1,
        );

        let grosterme = ncrossns2.dot(&dnplan.crossed(&ns2)) / norm2 / norm2;
        let resul2 = GpVec::zero().set_linear_form_3(
            self.ray2 / norm2 * (grosterme * ndotns2 - dnplan.dot(&ns2)),
            &nplan,
            -self.ray2 * ndotns2 / norm2,
            &dnplan,
            -self.ray2 * grosterme / norm2,
            &ns2,
        );

        d.set_value(2, 2, resul1.x() + resul2.x());
        d.set_value(3, 2, resul1.y() + resul2.y());
        d.set_value(4, 2, resul1.z() + resul2.z());

        // Derived compared to u2 (`cxx:607-616`).
        let temp = d2u2.crossed(&d1v2).added(&d1u2.crossed(&d2uv2));
        let grosterme = ncrossns2.dot(&nplan.crossed(&temp)) / norm2 / norm2;
        let resul1 = GpVec::zero()
            .set_linear_form_3(
                self.ray2 / norm2 * (grosterme * ndotns2 - nplan.dot(&temp)),
                &nplan,
                -self.ray2 * grosterme / norm2,
                &ns2,
                self.ray2 / norm2,
                &temp,
            )
            .subtracted(&d1u2);

        // Derived compared to v2 (`cxx:618-627`).
        let temp = d2uv2.crossed(&d1v2).added(&d1u2.crossed(&d2v2));
        let grosterme = ncrossns2.dot(&nplan.crossed(&temp)) / norm2 / norm2;
        let resul2 = GpVec::zero()
            .set_linear_form_3(
                self.ray2 / norm2 * (grosterme * ndotns2 - nplan.dot(&temp)),
                &nplan,
                -self.ray2 * grosterme / norm2,
                &ns2,
                self.ray2 / norm2,
                &temp,
            )
            .subtracted(&d1v2);

        if !self.first {
            d.set_value(2, 1, v2d.x() * resul1.x() + v2d.y() * resul2.x());
            d.set_value(3, 1, v2d.x() * resul1.y() + v2d.y() * resul2.y());
            d.set_value(4, 1, v2d.x() * resul1.z() + v2d.y() * resul2.z());
        } else {
            d.set_value(2, 3, resul1.x());
            d.set_value(3, 3, resul1.y());
            d.set_value(4, 3, resul1.z());

            d.set_value(2, 4, resul2.x());
            d.set_value(3, 4, resul2.y());
            d.set_value(4, 4, resul2.z());
        }

        true
    }
}

impl BlendFuncInv for ConstRadInv {
    /// `Set(OnFirst, COnSurf)` (`cxx:77-81`).
    fn set(&mut self, on_first: bool, c_on_surf: Option<Arc<dyn Curve2d>>) {
        self.first = on_first;
        self.csurf = c_on_surf;
    }

    /// `GetTolerance(Tolerance, Tol)` (`cxx:85-97`).
    fn get_tolerance(&self, tolerance: &mut MathVector, tol: f64) {
        // An empty `csurf` has no resolution; OCCT never reaches this with a
        // null handle.
        tolerance.set_value(1, match &self.csurf {
            Some(c) => c.resolution(tol),
            None => 0.0,
        });
        tolerance.set_value(2, self.curv.resolution(tol));
        if self.first {
            tolerance.set_value(3, approx_same_parameter::u_resolution(self.surf2.as_ref(), tol));
            tolerance.set_value(4, approx_same_parameter::v_resolution(self.surf2.as_ref(), tol));
        } else {
            tolerance.set_value(3, approx_same_parameter::u_resolution(self.surf1.as_ref(), tol));
            tolerance.set_value(4, approx_same_parameter::v_resolution(self.surf1.as_ref(), tol));
        }
    }

    /// `GetBounds(InfBound, SupBound)` (`cxx:101-147`): the box of `csurf` and
    /// the guide line, plus the box of the partner surface widened by its own
    /// range when both ends are finite.
    fn get_bounds(&self, inf_bound: &mut MathVector, sup_bound: &mut MathVector) {
        let (cfirst, clast) = match &self.csurf {
            Some(c) => (c.first_parameter(), c.last_parameter()),
            None => (0.0, 0.0),
        };
        inf_bound.set_value(1, cfirst);
        inf_bound.set_value(2, self.curv.first_parameter());
        sup_bound.set_value(1, clast);
        sup_bound.set_value(2, self.curv.last_parameter());

        let partner = if self.first { &self.surf2 } else { &self.surf1 };
        let (umin, umax) = partner.u_range();
        let (vmin, vmax) = partner.v_range();
        inf_bound.set_value(3, umin);
        inf_bound.set_value(4, vmin);
        sup_bound.set_value(3, umax);
        sup_bound.set_value(4, vmax);

        if !Precision::is_infinite(inf_bound.value(3)) && !Precision::is_infinite(sup_bound.value(3))
        {
            let range = sup_bound.value(3) - inf_bound.value(3);
            inf_bound.set_value(3, inf_bound.value(3) - range);
            sup_bound.set_value(3, sup_bound.value(3) + range);
        }
        if !Precision::is_infinite(inf_bound.value(4)) && !Precision::is_infinite(sup_bound.value(4))
        {
            let range = sup_bound.value(4) - inf_bound.value(4);
            inf_bound.set_value(4, inf_bound.value(4) - range);
            sup_bound.set_value(4, sup_bound.value(4) + range);
        }
    }

    /// `IsSolution(Sol, Tol)` (`cxx:149-153`): the first value and the norm of
    /// the last three are within `tol`.
    fn is_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        let mut valsol = MathVector::new(1, 4);
        self.value(sol, &mut valsol);
        valsol.value(1).abs() <= tol
            && valsol.value(2) * valsol.value(2)
                + valsol.value(3) * valsol.value(3)
                + valsol.value(4) * valsol.value(4)
                <= tol * tol
    }
}
