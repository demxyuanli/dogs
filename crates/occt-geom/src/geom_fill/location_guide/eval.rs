//! `GeomFill_LocationGuide::D0` / `D1` / `D2`.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_LocationGuide.cxx:584-1145`
//! (OCCT 8.0.0).
//!
//! Translation notes:
//! - `D0` (`cxx:584-650`) and `D0` with the `Poles2d` array (`cxx:657-727`)
//!   are two copies of the same body in OCCT; the array is unused, so the port
//!   has one body and `d0_with_poles` returns an empty vector.
//! - When `rotation` is set, `D0` solves the guide / revolution-surface system
//!   with `GeomFill_FunctionGuide` + `math_FunctionSetRoot` seeded by `InitX`
//!   and rotates the `N` / `B` columns about the tangent by the solution angle.
//! - `D1` (`cxx:733-891`) returns `false` immediately when `rotation` is set
//!   (`cxx:768-770`); the derivative branch below that line is commented out in
//!   OCCT and stays UNPORTED.
//! - `D2` (`cxx:896-1145`) has the same early `false` for `rotation`
//!   (`cxx:933-935`); its `WithTrans` block (`cxx:925-930`) runs before `M`,
//!   `DM` and `D2M` are written, so the untransformed law frame is what OCCT
//!   actually returns. Reproduced.

use occt_core::gp::{GpMat, GpVec};

use super::LocationGuide;
use crate::geom_fill::location_law::{LocationD1, LocationD2};
use crate::geom_fill::pipe_error::GeomFillPipeError;

impl LocationGuide {
    /// `GeomFill_LocationGuide::D0` (`cxx:584-650`).
    pub(crate) fn d0_frame(&self, param: f64) -> Result<(GpMat, GpVec), String> {
        let (curve, law, rotation, u_on_s, inf, sup, tol_res) = {
            let st = self.lock();
            let curve = st.curve.clone().ok_or_else(|| {
                "Standard_NullObject: GeomFill_LocationGuide::D0".to_string()
            })?;
            if st.sec.is_none() || st.guide.is_none() {
                return Err("Standard_NullObject: GeomFill_LocationGuide::D0".to_string());
            }
            (
                curve,
                st.law.clone(),
                st.rotation,
                st.first_s + (param - st.curve.as_ref().unwrap().first_parameter()) * st.ratio,
                st.inf,
                st.sup,
                st.tol_res,
            )
        };
        let p = curve.d0(param);
        let v = GpVec::new(p.x(), p.y(), p.z());
        let Ok((tangent, normal, binormal)) = law.d0(param) else {
            self.lock().status = law.error_status();
            return Err("GeomFill_LocationGuide::D0: trihedron law failed".into());
        };
        let mut m = GpMat::from_cols(normal.xyz(), binormal.xyz(), tangent.xyz());
        let with_trans = self.lock().with_trans;
        let trans = self.lock().trans;
        if with_trans {
            m = m.multiply(&trans);
        }
        if !rotation {
            return Ok((m, v));
        }

        // Seed (`cxx:609`) and solve (`cxx:620-624`).
        self.init_x(param);
        let st = self.lock();
        let Some(sec) = st.sec.clone() else {
            return Err("Standard_NullObject: GeomFill_LocationGuide::D0".into());
        };
        let Some(guide) = st.guide.clone() else {
            return Err("Standard_NullObject: GeomFill_LocationGuide::D0".into());
        };
        let seed = st.x;
        drop(st);

        let t = m.column(3);
        let b = m.column(2);
        let n = m.column(1);
        let solved = super::rotation::solve_guide(
            &sec,
            &guide,
            u_on_s,
            &p,
            &t,
            &n,
            seed,
            inf,
            sup,
            tol_res,
            100,
        );
        let Some((root, _)) = solved else {
            self.lock().status = GeomFillPipeError::ImpossibleContact;
            return Err("GeomFill_LocationGuide::D0: no result".into());
        };

        // `Rot.SetRotation(t, R(2)); b *= Rot; n *= Rot;`
        let mut rot = GpMat::identity();
        rot.set_rotation(&t, root[1]).map_err(|e| e.to_string())?;
        let b = b.multiplied_mat(&rot);
        let n = n.multiplied_mat(&rot);
        m = GpMat::from_cols(&n, &b, &t);
        Ok((m, v))
    }

    /// `GeomFill_LocationGuide::D1` (`cxx:733-891`).
    pub(crate) fn d1_frame(&self, param: f64) -> Result<LocationD1, String> {
        let st = self.lock();
        let curve = st
            .curve
            .clone()
            .ok_or_else(|| "Standard_NullObject: GeomFill_LocationGuide::D1".to_string())?;
        let law = st.law.clone();
        let rotation = st.rotation;
        let with_trans = st.with_trans;
        let trans = st.trans;
        drop(st);

        let (p, dv) = curve.d1(param);
        let v = GpVec::new(p.x(), p.y(), p.z());
        let Ok(frame) = law.d1(param) else {
            self.lock().status = law.error_status();
            return Err("GeomFill_LocationGuide::D1: trihedron law failed".into());
        };
        let mut m = GpMat::from_cols(frame.normal.xyz(), frame.binormal.xyz(), frame.tangent.xyz());
        let mut dm = GpMat::from_cols(
            frame.dnormal.xyz(),
            frame.dbinormal.xyz(),
            frame.dtangent.xyz(),
        );
        if with_trans {
            m = m.multiply(&trans);
            dm = dm.multiply(&trans);
        }
        if rotation {
            // `cxx:768-770`: OCCT gives up here; the rotation branch below is
            // commented out and stays UNPORTED.
            return Err("GeomFill_LocationGuide::D1: rotating location law".into());
        }
        Ok(LocationD1 {
            m,
            v,
            dm,
            dv,
            poles2d: Vec::new(),
            dpoles2d: Vec::new(),
        })
    }

    /// `GeomFill_LocationGuide::D2` (`cxx:896-1145`).
    pub(crate) fn d2_frame(&self, param: f64) -> Result<LocationD2, String> {
        let st = self.lock();
        let curve = st
            .curve
            .clone()
            .ok_or_else(|| "Standard_NullObject: GeomFill_LocationGuide::D2".to_string())?;
        let law = st.law.clone();
        let rotation = st.rotation;
        let with_trans = st.with_trans;
        let trans = st.trans;
        drop(st);

        let (p, dv, d2v) = curve.d2(param);
        let v = GpVec::new(p.x(), p.y(), p.z());
        let Ok(frame) = law.d2(param) else {
            self.lock().status = law.error_status();
            return Err("GeomFill_LocationGuide::D2: trihedron law failed".into());
        };
        let _ = with_trans;
        let _ = trans;
        if rotation {
            return Err("GeomFill_LocationGuide::D2: rotating location law".into());
        }
        // `cxx:1136-1143`: the non-rotation output, taken from the law without
        // the `WithTrans` product (that product at `cxx:925-930` runs on the
        // still-unwritten output parameters).
        let m = GpMat::from_cols(frame.normal.xyz(), frame.binormal.xyz(), frame.tangent.xyz());
        let dm = GpMat::from_cols(
            frame.dnormal.xyz(),
            frame.dbinormal.xyz(),
            frame.dtangent.xyz(),
        );
        let d2m = GpMat::from_cols(
            frame.d2normal.xyz(),
            frame.d2binormal.xyz(),
            frame.d2tangent.xyz(),
        );
        Ok(LocationD2 {
            m,
            v,
            dm,
            dv,
            d2m,
            d2v,
            poles2d: Vec::new(),
            dpoles2d: Vec::new(),
            d2poles2d: Vec::new(),
        })
    }
}
