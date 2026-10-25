//! `GeomFill_LocationDraft` evaluation: the frame plus, while `Intersec` is
//! set, the 2d poles of the generatrix / stop-surface intersection.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_LocationDraft.cxx:233-645` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `D0(Param, M, V, Poles2d)` (`cxx:233-341`) evaluates the point on
//!   `myCurve` (`cxx:245`), the frame from `myLaw->D0` (`cxx:247`) and, when
//!   `Intersec` is set, seeds `(W0, U0, V0)` by linear interpolation over
//!   `myPoles2d` (`cxx:272-301`) and solves `GeomFill_FunctionDraft` with
//!   `math_NewtonFunctionSetRoot` (`cxx:318-321`). Tolerances are `1.e-5` on
//!   `X`, `1.e-7` on `F` and 100 iterations (`cxx:310-317`); `D1` uses
//!   `1.e-4` / `1.e-6` / 100 (`cxx:428-434`) and `D2` `1.e-4` / `1.e-6` / 150
//!   (`cxx:568-575`).
//! - `D1` adds `DerivT` + `Derivatives` and the `math_Gauss` solve
//!   `DEDX * DSDT = -DEDT` (`cxx:450-469`); `D2` adds `Deriv2X`, `Deriv2T`,
//!   `DerivTX` and the second `math_Gauss` solve
//!   `DEDX * D2SDT2 = -aT*DSDT - 2*D2EDTDX*DSDT - D2EDT2` (`cxx:597-641`).
//!   `math_Gauss` is [`MathMatrix::solve`]; when it reports singularity OCCT
//!   still returns `true` with the derivative poles untouched (`cxx:462-469`,
//!   `cxx:629-640`), which is what the `if let Ok(..)` arms below do.
//! - `Poles2d(1)` is the surface point `(U, V)` and `Poles2d(2)` the curve
//!   point `(W, Param)` (`cxx:328-331`, `cxx:445-450`, `cxx:586-591`).
//! - The generatrix / surface intersection is
//!   [`super::intersection::line_surface_intersections`]; a non-planar stop
//!   surface reports [`super::intersection::UNPORTED_NON_PLANE`].
//! - `IsEqual(a, b)` (`cxx:278`, `cxx:296`, ...) is `|a - b| < RealSmall()`
//!   (`Standard_Real.hxx:148-151`).

use std::sync::Arc;

use occt_core::gp::{GpDir, GpMat, GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::precision::REAL_SMALL;
use occt_math::{MathMatrix, MathNewtonFunctionSetRoot, MathVector};

use crate::curve::Curve;
use crate::geom_fill::function_draft::FunctionDraft;
use crate::geom_fill::location_draft::{LdState, LocationDraft};
use crate::geom_fill::location_law::{LocationD1, LocationD2};
use crate::geom_fill::tensor::Tensor;
use crate::geom_fill::trihedron_law::TrihedronLaw;
use crate::line::GeomLine;

/// `IsEqual(a, b)` (`Standard_Real.hxx:148-151`).
fn is_equal(a: f64, b: f64) -> bool {
    (a - b).abs() < REAL_SMALL
}

/// `math_Vector X(1, 3)` filled with `(W, U, V)`.
fn vector3(a: f64, b: f64, c: f64) -> MathVector {
    let mut v = MathVector::with_init(1, 3, 0.0);
    v.set_value(1, a);
    v.set_value(2, b);
    v.set_value(3, c);
    v
}

/// `math_Matrix * math_Vector` (`math_Matrix::Multiply(Vector)`).
fn mat_vec(m: &MathMatrix, v: &MathVector) -> MathVector {
    let mut r = MathVector::with_init(m.lower_row(), m.upper_row(), 0.0);
    for i in m.lower_row()..=m.upper_row() {
        let mut s = 0.0;
        for j in m.lower_col()..=m.upper_col() {
            s += m.value(i, j) * v.value(j);
        }
        r.set_value(i, s);
    }
    r
}

/// The `(W0, U0, V0)` seed of `cxx:272-301` (same code in `cxx:390-419` and
/// `cxx:531-559`).
fn seed(st: &LdState, param: f64) -> (f64, f64, f64) {
    let mut w0 = 0.0;
    let mut u0 = 0.0;
    let mut v0 = 0.0;
    let nb = 2 * st.nb_pts;
    let mut ii = 1usize;
    while ii < nb && st.poles2d[ii - 1].y() < param {
        ii += 2;
    }
    if ii < nb && !is_equal(st.poles2d[ii - 1].y(), param) {
        let t1 = st.poles2d[ii - 1].y();
        let t2 = st.poles2d[ii - 3].y();
        let param_t1 = (param - t1) / (t2 - t1);
        let t2_param = (t2 - param) / (t2 - t1);
        w0 = st.poles2d[ii - 3].x() * param_t1 + st.poles2d[ii - 1].x() * t2_param;
        u0 = st.poles2d[ii - 2].x() * param_t1 + st.poles2d[ii].x() * t2_param;
        v0 = st.poles2d[ii - 2].y() * param_t1 + st.poles2d[ii].y() * t2_param;
    } else if ii < nb {
        w0 = st.poles2d[ii - 1].x();
        u0 = st.poles2d[ii].x();
        v0 = st.poles2d[ii].y();
    }
    (w0, u0, v0)
}

impl LocationDraft {
    /// `GeomFill_FunctionDraft E(mySurf, G)` (`cxx:318`): the generatrix line
    /// through `p` along `d` against the stop surface.
    fn function_draft(
        &self,
        st: &LdState,
        p: &GpPnt,
        d: &GpVec,
    ) -> Result<FunctionDraft, String> {
        let surf = st
            .surf
            .clone()
            .ok_or_else(|| "Standard_NullObject: GeomFill_LocationDraft::FunctionDraft".to_string())?;
        let dir = GpDir::from_vec(d)
            .map_err(|_| "GeomFill_LocationDraft: null generatrix direction".to_string())?;
        let line: Arc<dyn Curve> = Arc::new(GeomLine::from_pnt_dir(*p, dir));
        Ok(FunctionDraft::new(surf, line))
    }

    /// `math_NewtonFunctionSetRoot Result(E, XTol, FTol, Iter);
    /// Result.Perform(E, X); Result.Root(R)` (`cxx:319-326`).
    fn newton_root(
        e: &mut FunctionDraft,
        w0: f64,
        u0: f64,
        v0: f64,
        tol_x: f64,
        tol_f: f64,
        iter: i32,
    ) -> Result<(f64, f64, f64), String> {
        let x = vector3(w0, u0, v0);
        let xtol = MathVector::with_init(1, 3, tol_x);
        let mut result = MathNewtonFunctionSetRoot::new(e, &xtol, tol_f, iter);
        result.perform(e, &x);
        if !result.is_done() {
            return Err("GeomFill_LocationDraft: intersection not found".to_string());
        }
        let r = result.root().map_err(|err| err.to_string())?;
        Ok((r.value(1), r.value(2), r.value(3)))
    }

    /// `D0(Param, M, V, Poles2d)` (`cxx:233-341`).
    pub(crate) fn d0_poles(&self, param: f64) -> Result<(GpMat, GpVec, Vec<GpPnt2d>), String> {
        let st = self.lock();
        let curve = st
            .curve
            .clone()
            .ok_or_else(|| "Standard_NullObject: GeomFill_LocationDraft::D0".to_string())?;
        let p = curve.d0(param);
        let v = GpVec::new(p.x(), p.y(), p.z());
        let (t, n, b) = st.law.d0(param)?;
        let mut m = GpMat::from_cols(n.xyz(), b.xyz(), t.xyz());
        if st.with_trans {
            m = m.multiply(&st.trans);
        }
        let mut poles = Vec::new();
        if st.intersec {
            // `D = cos(myAngle)*B + sin(myAngle)*N` (`cxx:262`).
            let d = b
                .multiplied_scalar(st.angle.cos())
                .added(&n.multiplied_scalar(st.angle.sin()));
            let (w0, u0, v0) = seed(&st, param);
            let mut e = self.function_draft(&st, &p, &d)?;
            let (rw, ru, rv) = Self::newton_root(&mut e, w0, u0, v0, 1.0e-5, 1.0e-7, 100)?;
            poles.push(GpPnt2d::new(ru, rv));
            poles.push(GpPnt2d::new(rw, param));
        }
        Ok((m, v, poles))
    }

    /// `D1` (`cxx:346-478`).
    pub(crate) fn d1_frame(&self, param: f64) -> Result<LocationD1, String> {
        let st = self.lock();
        let curve = st
            .curve
            .clone()
            .ok_or_else(|| "Standard_NullObject: GeomFill_LocationDraft::D1".to_string())?;
        let (p, dv) = curve.d1(param);
        let v = GpVec::new(p.x(), p.y(), p.z());
        let tr = st.law.d1(param)?;
        let mut m = GpMat::from_cols(tr.normal.xyz(), tr.binormal.xyz(), tr.tangent.xyz());
        let mut dm = GpMat::from_cols(tr.dnormal.xyz(), tr.dbinormal.xyz(), tr.dtangent.xyz());
        if st.with_trans {
            m = m.multiply(&st.trans);
            dm = dm.multiply(&st.trans);
        }
        let mut poles = Vec::new();
        let mut dpoles = Vec::new();
        if st.intersec {
            let d = tr
                .binormal
                .multiplied_scalar(st.angle.cos())
                .added(&tr.normal.multiplied_scalar(st.angle.sin()));
            let (w0, u0, v0) = seed(&st, param);
            let mut e = self.function_draft(&st, &p, &d)?;
            let (rw, ru, rv) = Self::newton_root(&mut e, w0, u0, v0, 1.0e-4, 1.0e-6, 100)?;
            poles.push(GpPnt2d::new(ru, rv));
            poles.push(GpPnt2d::new(rw, param));

            let trimmed = st
                .trimmed
                .clone()
                .ok_or_else(|| "Standard_NullObject: GeomFill_LocationDraft::D1".to_string())?;
            let mut dedt = MathVector::with_init(1, 3, 0.0);
            e.deriv_t(trimmed.as_ref(), param, rw, &tr.dnormal, st.angle, &mut dedt);
            let mut dedx = MathMatrix::with_init(1, 3, 1, 3, 0.0);
            e.derivatives(&vector3(rw, ru, rv), &mut dedx);
            if let Ok(dsdt) = dedx.solve(&dedt.opposite()) {
                dpoles.push(GpVec2d::new(dsdt.value(2), dsdt.value(3)));
                dpoles.push(GpVec2d::new(dsdt.value(1), 1.0));
            }
        }
        Ok(LocationD1 {
            m,
            v,
            dm,
            dv,
            poles2d: poles,
            dpoles2d: dpoles,
        })
    }

    /// `D2` (`cxx:483-645`).
    pub(crate) fn d2_frame(&self, param: f64) -> Result<LocationD2, String> {
        let st = self.lock();
        let curve = st
            .curve
            .clone()
            .ok_or_else(|| "Standard_NullObject: GeomFill_LocationDraft::D2".to_string())?;
        let (p, dv, d2v) = curve.d2(param);
        let v = GpVec::new(p.x(), p.y(), p.z());
        let tr = st.law.d2(param)?;
        let mut m = GpMat::from_cols(tr.normal.xyz(), tr.binormal.xyz(), tr.tangent.xyz());
        let mut dm = GpMat::from_cols(tr.dnormal.xyz(), tr.dbinormal.xyz(), tr.dtangent.xyz());
        let mut d2m = GpMat::from_cols(tr.d2normal.xyz(), tr.d2binormal.xyz(), tr.d2tangent.xyz());
        if st.with_trans {
            m = m.multiply(&st.trans);
            dm = dm.multiply(&st.trans);
            d2m = d2m.multiply(&st.trans);
        }
        let mut poles = Vec::new();
        let mut dpoles = Vec::new();
        let mut d2poles = Vec::new();
        if st.intersec {
            let d = tr
                .binormal
                .multiplied_scalar(st.angle.cos())
                .added(&tr.normal.multiplied_scalar(st.angle.sin()));
            let (w0, u0, v0) = seed(&st, param);
            let mut e = self.function_draft(&st, &p, &d)?;
            let (rw, ru, rv) = Self::newton_root(&mut e, w0, u0, v0, 1.0e-4, 1.0e-6, 150)?;
            poles.push(GpPnt2d::new(ru, rv));
            poles.push(GpPnt2d::new(rw, param));

            let trimmed = st
                .trimmed
                .clone()
                .ok_or_else(|| "Standard_NullObject: GeomFill_LocationDraft::D2".to_string())?;
            let x = vector3(rw, ru, rv);

            // `DEDT` / `DEDX` / `DSDT` (`cxx:594-611`).
            let mut dedt = MathVector::with_init(1, 3, 0.0);
            e.deriv_t(trimmed.as_ref(), param, rw, &tr.dnormal, st.angle, &mut dedt);
            let mut dedx = MathMatrix::with_init(1, 3, 1, 3, 0.0);
            e.derivatives(&x, &mut dedx);
            let mut dsdt = MathVector::with_init(1, 3, 0.0);
            if let Ok(sol) = dedx.solve(&dedt.opposite()) {
                dsdt = sol;
                dpoles.push(GpVec2d::new(dsdt.value(2), dsdt.value(3)));
                dpoles.push(GpVec2d::new(dsdt.value(1), 1.0));
            }

            // `D2EDX2`, `D2EDT2`, `D2EDTDX`, `aT = D2EDX2 * DSDT` (`cxx:613-627`).
            let mut d2edx2 = Tensor::new(3, 3, 3);
            e.deriv2_x(&x, &mut d2edx2);
            let mut d2edt2 = MathVector::with_init(1, 3, 0.0);
            e.deriv2_t(trimmed.as_ref(), param, rw, &tr.d2normal, st.angle, &mut d2edt2);
            let mut d2edtdx = MathMatrix::with_init(1, 3, 1, 3, 0.0);
            e.deriv_tx(&tr.dnormal, st.angle, &mut d2edtdx);
            let mut at = MathMatrix::with_init(1, 3, 1, 3, 0.0);
            d2edx2.multiply(&dsdt, &mut at);

            // `-aT*DSDT - 2*D2EDTDX*DSDT - D2EDT2` (`cxx:631`).
            let mut rhs = mat_vec(&at, &dsdt).multiplied_scalar(-1.0);
            rhs.add(&mat_vec(&d2edtdx, &dsdt).multiplied_scalar(-2.0));
            rhs.subtract(&d2edt2);
            if let Ok(d2sdt2) = dedx.solve(&rhs) {
                d2poles.push(GpVec2d::new(d2sdt2.value(2), d2sdt2.value(3)));
                d2poles.push(GpVec2d::new(d2sdt2.value(1), 0.0));
            }
        }
        Ok(LocationD2 {
            m,
            v,
            dm,
            dv,
            d2m,
            d2v,
            poles2d: poles,
            dpoles2d: dpoles,
            d2poles2d: d2poles,
        })
    }
}
