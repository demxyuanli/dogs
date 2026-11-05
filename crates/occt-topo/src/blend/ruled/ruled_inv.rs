//! `BlendFunc_RuledInv` — the inverse function of `BlendFunc_Ruled`: a
//! solution on a restriction (curve on surface) of one of the two supports.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_RuledInv.hxx:25-73`,
//! `BlendFunc_RuledInv.cxx:21-508` (OCCT 8.0.0).
//!
//! The vector `X` is `(t, w, U, V)`: the parameter on the curve on surface, the
//! parameter on the guide line, and the parameters of the partner point on the
//! other surface. The four equations are those of `BlendFunc_Ruled`:
//! 1. `nplan . pts1 + theD`,
//! 2. `nplan . pts2 + theD`,
//! 3. `(pts2 - pts1) . nor1`,
//! 4. `(pts2 - pts1) . nor2`,
//!
//! with `nplan` the unit tangent of the guide line at `X(2)` and `theD` its
//! plane offset; `nor1` / `nor2` are the surface normals projected into the
//! section plane.
//!
//! Translation notes:
//! - `csurf` is an `occ::handle<Adaptor2d_Curve2d>`; the port stores it as
//!   `Option<Arc<dyn Curve2d>>` because `Set(OnFirst, C)` can be handed the
//!   empty handle. OCCT dereferences it unguarded in `Value` / `Derivatives` /
//!   `Values`; the port answers `false` instead of faulting. `GetTolerance`
//!   and `GetBounds` keep the null case at `0.`.
//! - `Derivatives` (`cxx:156-327`) and `Values` (`cxx:329-508`) repeat the same
//!   body in OCCT, `Values` adding the four residuals; the port writes it once
//!   with `f: Option<&mut MathVector>`.
//! - The `D(3, ..)` / `D(4, ..)` block depends on `first` in opposite ways:
//!   the curve-on-surface variable is `X(1)` when the restriction lies on
//!   `surf1`, `X(3)` / `X(4)` otherwise (`cxx:224-236`, `:262-272`).
//! - `IsSolution` tests all four residuals one by one (`cxx:107-112`).

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d, GpVec};
use occt_core::precision::Precision;
use occt_geom::approx_same_parameter;
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_geom2d::curve::Curve2d;
use occt_math::{MathFunctionSetWithDerivatives, MathMatrix, MathVector};

use crate::blend::BlendFuncInv;

/// `class BlendFunc_RuledInv : public Blend_FuncInv` (`hxx:25-73`).
pub struct RuledInv {
    /// `surf1` (`hxx:66`).
    surf1: Arc<dyn Surface>,
    /// `surf2` (`hxx:67`).
    surf2: Arc<dyn Surface>,
    /// `curv` (`hxx:68`): the guide line.
    curv: Arc<dyn Curve>,
    /// `csurf` (`hxx:69`): the curve on surface; `Set(OnFirst, C)` writes it.
    csurf: Option<Arc<dyn Curve2d>>,
    /// `first` (`hxx:70`): `true` when `csurf` lies on `surf1`.
    first: bool,
}

impl RuledInv {
    /// `BlendFunc_RuledInv(S1, S2, C)` (`cxx:21-29`): `csurf` stays null and
    /// `first = false`.
    pub fn new(surf1: Arc<dyn Surface>, surf2: Arc<dyn Surface>, curv: Arc<dyn Curve>) -> Self {
        Self {
            surf1,
            surf2,
            curv,
            csurf: None,
            first: false,
        }
    }

    /// The surface that is not restricted by `csurf` (`first` selects it).
    fn partner(&self) -> &Arc<dyn Surface> {
        if self.first {
            &self.surf2
        } else {
            &self.surf1
        }
    }

    /// `Value(X, F)` (`cxx:114-154`). A missing `csurf` answers `false`.
    fn compute_value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        let csurf = match &self.csurf {
            Some(c) => Arc::clone(c),
            None => return false,
        };

        let (ptcur, d1cur) = self.curv.d1(x.value(2));
        let nplan = d1cur.normalized();
        let the_d = -nplan.dot(&GpVec::from_xyz(ptcur.xyz()));

        let pt2d = csurf.value(x.value(1));

        let (pts1, d1u1, d1v1, pts2, d1u2, d1v2) = self.surface_data(x, &pt2d);

        let temp = GpVec::from_pnts(&pts1, &pts2);

        let ns1 = d1u1.crossed(&d1v1);
        let ns2 = d1u2.crossed(&d1v2);
        let norm1 = nplan.crossed(&ns1).magnitude();
        let norm2 = nplan.crossed(&ns2).magnitude();
        let ns1 = GpVec::zero().set_linear_form_2(
            nplan.dot(&ns1) / norm1,
            &nplan,
            -1.0 / norm1,
            &ns1,
        );
        let ns2 = GpVec::zero().set_linear_form_2(
            nplan.dot(&ns2) / norm2,
            &nplan,
            -1.0 / norm2,
            &ns2,
        );

        f.set_value(1, nplan.dot(&GpVec::from_xyz(pts1.xyz())) + the_d);
        f.set_value(2, nplan.dot(&GpVec::from_xyz(pts2.xyz())) + the_d);
        f.set_value(3, temp.dot(&ns1));
        f.set_value(4, temp.dot(&ns2));

        true
    }

    /// `Derivatives(X, D)` (`cxx:156-327`) and `Values(X, F, D)`
    /// (`cxx:329-508`). `f` is `None` for `Derivatives`.
    fn compute_values(
        &mut self,
        x: &MathVector,
        f: Option<&mut MathVector>,
        d: &mut MathMatrix,
    ) -> bool {
        let csurf = match &self.csurf {
            Some(c) => Arc::clone(c),
            None => return false,
        };

        let (ptcur, d1cur, d2cur) = self.curv.d2(x.value(2));

        let normtgcur = d1cur.magnitude();
        let nplan = d1cur.normalized();
        let the_d = -nplan.dot(&GpVec::from_xyz(ptcur.xyz()));

        let dnplan = GpVec::zero()
            .set_linear_form_2(-nplan.dot(&d2cur), &nplan, 1.0, &d2cur)
            .divide(normtgcur);

        let (pt2d, v2d) = csurf.d1(x.value(1));

        let (pts1, d1u1, d1v1, d2u1, d2v1, d2uv1, pts2, d1u2, d1v2, d2u2, d2v2, d2uv2) =
            self.surface_data2(x, &pt2d);

        let dpdt = if self.first {
            GpVec::zero().set_linear_form_2(v2d.x(), &d1u1, v2d.y(), &d1v1)
        } else {
            GpVec::zero().set_linear_form_2(v2d.x(), &d1u2, v2d.y(), &d1v2)
        };
        let p1p2 = GpVec::from_pnts(&pts1, &pts2);

        let from_guide = |p: &GpPnt| dnplan.dot(&GpVec::from_pnts(&ptcur, p)) - normtgcur;

        if self.first {
            d.set_value(1, 1, dpdt.dot(&nplan));
            d.set_value(1, 2, from_guide(&pts1));
            d.set_value(1, 3, 0.0);
            d.set_value(1, 4, 0.0);

            d.set_value(2, 1, 0.0);
            d.set_value(2, 2, from_guide(&pts2));
            d.set_value(2, 3, d1u2.dot(&nplan));
            d.set_value(2, 4, d1v2.dot(&nplan));
        } else {
            d.set_value(1, 1, 0.0);
            d.set_value(1, 2, from_guide(&pts1));
            d.set_value(1, 3, d1u1.dot(&nplan));
            d.set_value(1, 4, d1v1.dot(&nplan));

            d.set_value(2, 1, dpdt.dot(&nplan));
            d.set_value(2, 2, from_guide(&pts2));
            d.set_value(2, 3, 0.0);
            d.set_value(2, 4, 0.0);
        }

        let ns1 = d1u1.crossed(&d1v1);
        let ns2 = d1u2.crossed(&d1v2);
        let ncrossns1 = nplan.crossed(&ns1);
        let ncrossns2 = nplan.crossed(&ns2);
        let norm1 = ncrossns1.magnitude();
        let norm2 = ncrossns2.magnitude();

        let ndotns1 = nplan.dot(&ns1);
        let ndotns2 = nplan.dot(&ns2);

        let nor1 = GpVec::zero().set_linear_form_2(
            ndotns1 / norm1,
            &nplan,
            -1.0 / norm1,
            &ns1,
        );
        let nor2 = GpVec::zero().set_linear_form_2(
            ndotns2 / norm2,
            &nplan,
            -1.0 / norm2,
            &ns2,
        );

        if self.first {
            d.set_value(3, 3, d1u2.dot(&nor1));
            d.set_value(3, 4, d1v2.dot(&nor1));

            d.set_value(4, 1, -dpdt.dot(&nor2));
        } else {
            d.set_value(3, 1, dpdt.dot(&nor1));

            d.set_value(4, 3, -d1u1.dot(&nor2));
            d.set_value(4, 4, -d1v1.dot(&nor2));
        }

        let mut resul1;
        let mut resul2;
        let mut temp;

        // Derivative of nor1 with respect to u1 (`cxx:242-250`).
        temp = d2u1.crossed(&d1v1).added(&d1u1.crossed(&d2uv1));
        let mut grosterme = ncrossns1.dot(&nplan.crossed(&temp)) / norm1 / norm1;
        resul1 = GpVec::zero().set_linear_form_3(
            -(grosterme * ndotns1 - nplan.dot(&temp)) / norm1,
            &nplan,
            grosterme / norm1,
            &ns1,
            -1.0 / norm1,
            &temp,
        );

        // Derivative with respect to v1 (`cxx:254-262`).
        temp = d2uv1.crossed(&d1v1).added(&d1u1.crossed(&d2v1));
        grosterme = ncrossns1.dot(&nplan.crossed(&temp)) / norm1 / norm1;
        resul2 = GpVec::zero().set_linear_form_3(
            -(grosterme * ndotns1 - nplan.dot(&temp)) / norm1,
            &nplan,
            grosterme / norm1,
            &ns1,
            -1.0 / norm1,
            &temp,
        );

        if self.first {
            resul1 = GpVec::zero().set_linear_form_2(v2d.x(), &resul1, v2d.y(), &resul2);
            d.set_value(3, 1, p1p2.dot(&resul1) - dpdt.dot(&nor1));
        } else {
            d.set_value(3, 3, -d1u1.dot(&nor1) + p1p2.dot(&resul1));
            d.set_value(3, 4, -d1v1.dot(&nor1) + p1p2.dot(&resul2));
        }

        // Derivative of nor2 with respect to u2 (`cxx:276-282`).
        temp = d2u2.crossed(&d1v2).added(&d1u2.crossed(&d2uv2));
        grosterme = ncrossns2.dot(&nplan.crossed(&temp)) / norm2 / norm2;
        resul1 = GpVec::zero().set_linear_form_3(
            -(grosterme * ndotns2 - nplan.dot(&temp)) / norm2,
            &nplan,
            grosterme / norm2,
            &ns2,
            -1.0 / norm2,
            &temp,
        );

        // Derivative with respect to v2 (`cxx:286-292`).
        temp = d2uv2.crossed(&d1v2).added(&d1u2.crossed(&d2v2));
        grosterme = ncrossns2.dot(&nplan.crossed(&temp)) / norm2 / norm2;
        resul2 = GpVec::zero().set_linear_form_3(
            -(grosterme * ndotns2 - nplan.dot(&temp)) / norm2,
            &nplan,
            grosterme / norm2,
            &ns2,
            -1.0 / norm2,
            &temp,
        );

        if self.first {
            d.set_value(4, 3, d1u2.dot(&nor2) + p1p2.dot(&resul1));
            d.set_value(4, 4, d1v2.dot(&nor2) + p1p2.dot(&resul2));
        } else {
            resul1 = GpVec::zero().set_linear_form_2(v2d.x(), &resul1, v2d.y(), &resul2);
            d.set_value(4, 1, p1p2.dot(&resul1) + dpdt.dot(&nor2));
        }

        // Derivatives with respect to w, the parameter on the guide line
        // (`cxx:306-325`).
        grosterme = ncrossns1.dot(&dnplan.crossed(&ns1)) / norm1 / norm1;
        resul1 = GpVec::zero().set_linear_form_3(
            -(grosterme * ndotns1 - dnplan.dot(&ns1)) / norm1,
            &nplan,
            ndotns1 / norm1,
            &dnplan,
            grosterme / norm1,
            &ns1,
        );

        grosterme = ncrossns2.dot(&dnplan.crossed(&ns2)) / norm2 / norm2;
        resul2 = GpVec::zero().set_linear_form_3(
            -(grosterme * ndotns2 - dnplan.dot(&ns2)) / norm2,
            &nplan,
            ndotns2 / norm2,
            &dnplan,
            grosterme / norm2,
            &ns2,
        );

        d.set_value(3, 2, p1p2.dot(&resul1));
        d.set_value(4, 2, p1p2.dot(&resul2));

        if let Some(f) = f {
            f.set_value(1, nplan.dot(&GpVec::from_xyz(pts1.xyz())) + the_d);
            f.set_value(2, nplan.dot(&GpVec::from_xyz(pts2.xyz())) + the_d);
            f.set_value(3, p1p2.dot(&nor1));
            f.set_value(4, p1p2.dot(&nor2));
        }

        true
    }

    /// `D1` of the two surfaces for `Value`: `(pts1, d1u1, d1v1, pts2, d1u2,
    /// d1v2)`, where the restricted surface is evaluated at the point of
    /// `csurf` (`cxx:126-135`).
    #[allow(clippy::type_complexity)]
    fn surface_data(
        &self,
        x: &MathVector,
        pt2d: &GpPnt2d,
    ) -> (GpPnt, GpVec, GpVec, GpPnt, GpVec, GpVec) {
        let (u1, v1, u2, v2) = self.surface_parameters(x, pt2d);
        let (pts1, d1u1, d1v1) = self.surf1.d1(u1, v1);
        let (pts2, d1u2, d1v2) = self.surf2.d1(u2, v2);
        (pts1, d1u1, d1v1, pts2, d1u2, d1v2)
    }

    /// `D2` of the two surfaces for `Derivatives` / `Values` (`cxx:177-197`).
    #[allow(clippy::type_complexity)]
    fn surface_data2(
        &self,
        x: &MathVector,
        pt2d: &GpPnt2d,
    ) -> (
        GpPnt,
        GpVec,
        GpVec,
        GpVec,
        GpVec,
        GpVec,
        GpPnt,
        GpVec,
        GpVec,
        GpVec,
        GpVec,
        GpVec,
    ) {
        let (u1, v1, u2, v2) = self.surface_parameters(x, pt2d);
        let (pts1, d1u1, d1v1, d2u1, d2v1, d2uv1) = self.surf1.d2(u1, v1);
        let (pts2, d1u2, d1v2, d2u2, d2v2, d2uv2) = self.surf2.d2(u2, v2);
        (
            pts1, d1u1, d1v1, d2u1, d2v1, d2uv1, pts2, d1u2, d1v2, d2u2, d2v2, d2uv2,
        )
    }

    /// `(U1, V1, U2, V2)` with the restricted surface taken at `pt2d`
    /// (`cxx:126-135`, `:177-197`).
    fn surface_parameters(
        &self,
        x: &MathVector,
        pt2d: &GpPnt2d,
    ) -> (f64, f64, f64, f64) {
        if self.first {
            (pt2d.x(), pt2d.y(), x.value(3), x.value(4))
        } else {
            (x.value(3), x.value(4), pt2d.x(), pt2d.y())
        }
    }
}

impl MathFunctionSetWithDerivatives for RuledInv {
    /// `Blend_FuncInv::NbVariables()` = `4` (`Blend_FuncInv.cxx:18-21`).
    fn nb_variables(&self) -> usize {
        4
    }

    /// `NbEquations()` (`cxx:37-40`): `4`.
    fn nb_equations(&self) -> usize {
        4
    }

    /// `Value(X, F)` (`cxx:114-154`).
    fn value(&mut self, x: &MathVector, f: &mut MathVector) -> bool {
        self.compute_value(x, f)
    }

    /// `Derivatives(X, D)` (`cxx:156-327`).
    fn derivatives(&mut self, x: &MathVector, d: &mut MathMatrix) -> bool {
        self.compute_values(x, None, d)
    }

    /// `Values(X, F, D)` (`cxx:329-508`).
    fn values(&mut self, x: &MathVector, f: &mut MathVector, d: &mut MathMatrix) -> bool {
        self.compute_values(x, Some(f), d)
    }
}

impl BlendFuncInv for RuledInv {
    /// `Set(OnFirst, COnSurf)` (`cxx:31-35`).
    fn set(&mut self, on_first: bool, c_on_surf: Option<Arc<dyn Curve2d>>) {
        self.first = on_first;
        self.csurf = c_on_surf;
    }

    /// `GetTolerance(Tolerance, Tol)` (`cxx:42-57`).
    fn get_tolerance(&self, tolerance: &mut MathVector, tol: f64) {
        // An empty `csurf` has no resolution; OCCT never reaches this with a
        // null handle.
        tolerance.set_value(1, match &self.csurf {
            Some(c) => c.resolution(tol),
            None => 0.0,
        });
        tolerance.set_value(2, self.curv.resolution(tol));
        let partner = self.partner();
        tolerance.set_value(3, approx_same_parameter::u_resolution(partner.as_ref(), tol));
        tolerance.set_value(4, approx_same_parameter::v_resolution(partner.as_ref(), tol));
    }

    /// `GetBounds(InfBound, SupBound)` (`cxx:59-104`): the box of `csurf` and
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

        let partner = self.partner();
        let (umin, umax) = partner.u_range();
        let (vmin, vmax) = partner.v_range();
        inf_bound.set_value(3, umin);
        inf_bound.set_value(4, vmin);
        sup_bound.set_value(3, umax);
        sup_bound.set_value(4, vmax);

        for i in 3..=4 {
            if !Precision::is_infinite(inf_bound.value(i))
                && !Precision::is_infinite(sup_bound.value(i))
            {
                let range = sup_bound.value(i) - inf_bound.value(i);
                inf_bound.set_value(i, inf_bound.value(i) - range);
                sup_bound.set_value(i, sup_bound.value(i) + range);
            }
        }
    }

    /// `IsSolution(Sol, Tol)` (`cxx:106-112`): all four residuals within
    /// `tol`.
    fn is_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        let mut valsol = MathVector::new(1, 4);
        self.value(sol, &mut valsol);
        valsol.value(1).abs() <= tol
            && valsol.value(2).abs() <= tol
            && valsol.value(3).abs() <= tol
            && valsol.value(4).abs() <= tol
    }
}
