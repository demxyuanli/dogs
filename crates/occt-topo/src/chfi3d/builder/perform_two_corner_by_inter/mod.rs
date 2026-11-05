//! `ChFi3d_Builder::PerformTwoCornerbyInter` - the two-stripe corner taken by
//! intersecting the two fillet surfaces.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_C2.cxx:108-806`
//! (OCCT 8.0.0). The body is split by OCCT line ranges:
//!
//! - this file - `cxx:108-305`: the two stripes and their `SurfData`, the two
//!   `ChFi3d_IsInFront` calls that locate the common face and the two tangent
//!   lines, the `pivot` arc and its `sameparam` test, the two
//!   `ChFi3d_BoundSurf` and the `ok1` / `ok2` test on `myEFMap(pivot)`,
//! - [`same_param`] - `cxx:306-446`: the two p-curves meet on the common face,
//!   so one `ChFi3d_ComputeCurves` between the two fillet surfaces answers both
//!   stripes,
//! - [`diff_param`] - `cxx:447-804`: the two p-curves do not meet, so the body
//!   finds which stripe is the "big" one, intersects it with the small one and
//!   then with the face at the end of the small one, and writes the two pieces.
//!
//! Translation notes:
//! - OCCT keeps `Fd1` / `Fd2` as `ChFiDS_SurfData&` into the stripes'
//!   sequences; the port copies them into [`TwoCorner::fd1`] / `fd2` and writes
//!   them back with [`TwoCorner::flush`], which is also what keeps the two
//!   writes coherent when both stripes are the same handle (`Corner2 ==
//!   Corner1`),
//! - the `bool` answer of OCCT is the `Ok(())` / `self.done` pair: every
//!   `done = false; return done;` arm answers `Ok(false)` from the phase that
//!   took it and leaves `self.done` false, and only the `done = true` of
//!   `cxx:802-804` sets it,
//! - the `if (!pivot.IsSame(CP2.Arc()))` block of `cxx:256-269` is not
//!   translated: `pivot = CP1.Arc()` (`cxx:249-250`) and the `!CP1.Arc()
//!   .IsSame(CP2.Arc())` test of `cxx:236-246` has already returned, so the
//!   condition is never true. The block needs `GeomLib::ExtendCurveToPoint`
//!   and `Extrema_ExtPC`, neither of which is reachable from here,
//! - a null handle OCCT dereferences (a missing stripe, a `SurfData` out of
//!   range, a common point without an arc, a `myEFMap` entry missing) is the
//!   `Err` of this port.

mod diff_param;
mod same_param;

use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_geom2d::curve::Curve2d;

use crate::chfi3d::builder_0_helpers::{
    ch_fi3d_bound_surf, ch_fi3d_index_of_surf_data, ch_fi3d_is_in_front, IsInFrontState,
};
use crate::chfi_ds::{Stripe, StripeHandle, SurfData};
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::is_same;
use crate::topopebrep::tool::adaptor::{CurveAdaptor, SurfaceAdaptor};

use super::ChFi3dBuilder;

/// The locals of `PerformTwoCornerbyInter` that outlive the prologue.
pub(super) struct TwoCorner {
    /// `Vtx`
    pub vtx: TopoShape,
    /// `pivot`
    pub pivot: TopoShape,
    /// `parCP1`
    pub par_cp1: f64,
    /// `parCP2`
    pub par_cp2: f64,
    /// `Corner1`
    pub corner1: StripeHandle,
    /// `Corner2`
    pub corner2: StripeHandle,
    /// `IFd1`
    pub ifd1: i32,
    /// `IFd2`
    pub ifd2: i32,
    /// `Fd1`
    pub fd1: SurfData,
    /// `Fd2`
    pub fd2: SurfData,
    /// `isfirst1`
    pub is_first1: bool,
    /// `isfirst2`
    pub is_first2: bool,
    /// `UIntPC1`
    pub u_int_pc1: f64,
    /// `UIntPC2`
    pub u_int_pc2: f64,
    /// `IFaCo1`
    pub i_fa_co1: i32,
    /// `IFaCo2`
    pub i_fa_co2: i32,
    /// `IFaArc1`
    pub i_fa_arc1: i32,
    /// `IFaArc2`
    pub i_fa_arc2: i32,
    /// `sameparam`
    pub same_param: bool,
    /// `HS1` - mutated in place by `Reduce` in both branches
    pub hs1: SurfaceAdaptor,
    /// `HS2`
    pub hs2: SurfaceAdaptor,
    /// `Pardeb(1 .. 4)`
    pub pardeb: [f64; 4],
    /// `Parfin(1 .. 4)`
    pub parfin: [f64; 4],
}

impl TwoCorner {
    /// Writes the two copied `SurfData` back into their stripes.
    ///
    /// `flush` writes both slots of one sequence in a single pass when the two
    /// stripes are the same handle, so the second write cannot clobber the
    /// first.
    pub(super) fn flush(&mut self) {
        if Arc::ptr_eq(&self.corner1, &self.corner2) {
            if let Ok(mut stripe) = self.corner1.write() {
                let seq = stripe.change_set_of_surf_data();
                if let Some(slot) = index_slot(seq, self.ifd1) {
                    *Arc::make_mut(slot) = self.fd1.clone();
                }
                if let Some(slot) = index_slot(seq, self.ifd2) {
                    *Arc::make_mut(slot) = self.fd2.clone();
                }
            }
            return;
        }
        for (stripe, index, fd) in [
            (&self.corner1, self.ifd1, &self.fd1),
            (&self.corner2, self.ifd2, &self.fd2),
        ] {
            if let Ok(mut guard) = stripe.write() {
                if let Some(slot) = index_slot(guard.change_set_of_surf_data(), index) {
                    *Arc::make_mut(slot) = fd.clone();
                }
            }
        }
    }
}

/// `Seq.Value(i)` with the 1-based `i` of `NCollection_Sequence`.
fn index_slot(seq: &mut [Arc<SurfData>], index: i32) -> Option<&mut Arc<SurfData>> {
    index
        .checked_sub(1)
        .and_then(|i| usize::try_from(i).ok())
        .and_then(move |i| seq.get_mut(i))
}

/// `Fi.PCurveOnSurf()` (`cxx:315`, `:320`, `:346`, `:517`), a null handle as
/// `Err`.
pub(super) fn p_curve_on_surf(fd: &SurfData, index: i32) -> Result<Arc<dyn Curve2d>, String> {
    fd.interference(index)
        .p_curve_on_surf()
        .cloned()
        .ok_or_else(|| {
            format!("PerformTwoCornerbyInter: Interference({index}) has no p-curve on surf")
        })
}

/// `HS->Value(U, V)`.
pub(super) fn surf_value(hs: &SurfaceAdaptor, u: f64, v: f64) -> Result<GpPnt, String> {
    hs.value(u, v)
        .ok_or_else(|| "PerformTwoCornerbyInter: the bounded surface is null".to_string())
}

/// `IndexPoint(isfirst, on_s)` of a stripe.
pub(super) fn index_point(
    stripe: &StripeHandle,
    is_first: bool,
    on_s: i32,
) -> Result<i32, String> {
    let guard = stripe
        .read()
        .map_err(|_| "PerformTwoCornerbyInter: poisoned stripe".to_string())?;
    Ok(guard.index_point(is_first, on_s))
}

/// `ChangePCurve(First)` (`ChFiDS_Stripe.cxx`): the p-curve slot of one end.
pub(super) fn p_curve_mut(stripe: &mut Stripe, is_first: bool) -> &mut Option<Arc<dyn Curve2d>> {
    if is_first {
        stripe.change_first_p_curve()
    } else {
        stripe.change_last_p_curve()
    }
}

/// `Corner->SetOfSurfData()->Value(index)` (`cxx:127-129`), `index` 1-based.
pub(super) fn surf_data(stripe: &StripeHandle, index: i32) -> Result<SurfData, String> {    let guard = stripe
        .read()
        .map_err(|_| "PerformTwoCornerbyInter: poisoned stripe".to_string())?;
    index
        .checked_sub(1)
        .and_then(|i| usize::try_from(i).ok())
        .and_then(|i| guard.set_of_surf_data().get(i))
        .map(|sd| sd.as_ref().clone())
        .ok_or_else(|| format!("PerformTwoCornerbyInter: SurfData {index} is out of range"))
}

/// `Corner->SetOfSurfData()->Length()` (`cxx:138`).
fn surf_data_len(stripe: &StripeHandle) -> Result<i32, String> {
    let guard = stripe
        .read()
        .map_err(|_| "PerformTwoCornerbyInter: poisoned stripe".to_string())?;
    i32::try_from(guard.set_of_surf_data().len())
        .map_err(|_| "PerformTwoCornerbyInter: the SurfData sequence is too long".to_string())
}

/// One `ChFi3d_IsInFront` call (`cxx:159-201`).
///
/// `p1` / `p2` are in/out in OCCT: the second call of the body is seeded with
/// what the first one answered.
#[allow(clippy::too_many_arguments)]
fn is_in_front(
    b: &ChFi3dBuilder,
    corner1: &StripeHandle,
    corner2: &StripeHandle,
    ifd1: i32,
    ifd2: i32,
    sens1: i32,
    sens2: i32,
    vtx: &TopoShape,
    p1: f64,
    p2: f64,
    check_2d_distance: bool,
    enlarge: bool,
) -> Result<IsInFrontState, String> {
    let cd1 = corner1
        .read()
        .map_err(|_| "PerformTwoCornerbyInter: poisoned stripe".to_string())?;
    let cd2 = corner2
        .read()
        .map_err(|_| "PerformTwoCornerbyInter: poisoned stripe".to_string())?;
    let mut st = IsInFrontState::new(p1, p2);
    ch_fi3d_is_in_front(
        &b.my_ds,
        &cd1,
        &cd2,
        ifd1,
        ifd2,
        sens1,
        sens2,
        &mut st,
        Some(vtx),
        check_2d_distance,
        enlarge,
    )?;
    Ok(st)
}

/// `TopoDS::Face(DStr.Shape(index))` (`cxx:273-274`).
pub(super) fn face_at(ds: &crate::topopebrep::ds::DataStructure, index: i32) -> Result<Face, String> {
    let shape = ds
        .shape(index, true)
        .ok_or_else(|| format!("PerformTwoCornerbyInter: DStr.Shape({index}) is null"))?;
    Face::wrap(shape).ok_or_else(|| {
        format!("PerformTwoCornerbyInter: DStr.Shape({index}) is not a TopoDS_Face")
    })
}

impl ChFi3dBuilder {
    /// `PerformTwoCornerbyInter(Index)` (`ChFi3d_Builder_C2.cxx:108-806`).
    ///
    /// Answers `Ok(())` with `self.done` set like OCCT's `done` field: `true`
    /// on both successful completions (including the
    /// `PerformMoreThreeCorner(Index, 2)` hand-over of `cxx:242-247`) and
    /// `false` on the five early returns.
    pub fn perform_two_corner_by_inter(&mut self, index: usize) -> Result<(), String> {
        // `done = false;` (`cxx:111`).
        self.done = false;

        // `const TopoDS_Vertex& Vtx = myVDataMap.FindKey(Index);`
        // (`cxx:112`).
        let vtx = self
            .my_v_data_map
            .find_key(index)
            .cloned()
            .ok_or_else(|| format!("PerformTwoCornerbyInter: myVDataMap has no key {index}"))?;
        // `It.Initialize(myVDataMap(Index)); ...& Corner1 = It.ChangeValue();`
        // (`cxx:120-122`).
        let stripes: Vec<StripeHandle> = self
            .my_v_data_map
            .get_from_index(index)
            .map(|s| s.to_vec())
            .ok_or_else(|| {
                format!("PerformTwoCornerbyInter: myVDataMap({index}) is out of range")
            })?;
        let corner1 = stripes
            .first()
            .cloned()
            .ok_or_else(|| "PerformTwoCornerbyInter: the vertex has no stripe".to_string())?;
        let mut sens1 = 0;
        // `int IFd1 = ChFi3d_IndexOfSurfData(Vtx, Corner1, Sens1);`
        // (`cxx:124`).
        let ifd1 = ch_fi3d_index_of_surf_data(&vtx, &corner1, &mut sens1)?;
        let fd1 = surf_data(&corner1, ifd1)?;

        // `It.Next(); ...& Corner2 = It.ChangeValue();` (`cxx:130-131`). OCCT
        // dereferences the end iterator of a one-element list; the port answers
        // the `Err` instead.
        let corner2 = stripes
            .get(1)
            .cloned()
            .ok_or_else(|| "PerformTwoCornerbyInter: the vertex has a single stripe".to_string())?;
        // `if (Corner2 == Corner1) { Sens2 = -1; IFd2 = Corner2->SetOfSurfData()
        // ->Length(); } else { IFd2 = ChFi3d_IndexOfSurfData(Vtx, Corner2,
        // Sens2); }` (`cxx:134-143`).
        let (sens2, ifd2) = if Arc::ptr_eq(&corner1, &corner2) {
            (-1, surf_data_len(&corner2)?)
        } else {
            let mut s = 0;
            let i = ch_fi3d_index_of_surf_data(&vtx, &corner2, &mut s)?;
            (s, i)
        };
        let fd2 = surf_data(&corner2, ifd2)?;

        // `bool isfirst1 = (Sens1 == 1); bool isfirst2 = (Sens2 == 1);`
        // (`cxx:152-153`).
        let is_first1 = sens1 == 1;
        let is_first2 = sens2 == 1;

        // `OkinterCC = ChFi3d_IsInFront(DStr, Corner1, Corner2, IFd1, IFd2,
        // Sens1, Sens2, UIntPC1, UIntPC2, FaCo, SameSide, IFaCo1, IFaCo2,
        // Okvisavis, Vtx, true);` (`cxx:159-175`). `UIntPC1` / `UIntPC2` are
        // uninitialized out-parameters in OCCT; the port seeds them with `0.0`.
        let mut st = is_in_front(
            self,
            &corner1,
            &corner2,
            ifd1,
            ifd2,
            sens1,
            sens2,
            &vtx,
            0.0,
            0.0,
            true,
            false,
        )?;
        // `if (!Okvisavis) { done = false; return done; }` (`cxx:176-182`).
        if !st.visavis {
            self.done = false;
            return Ok(());
        }
        // `if (!OkinterCC) { OkinterCC = ChFi3d_IsInFront(..., true, true); }`
        // (`cxx:183-203`).
        if !st.ok {
            st = is_in_front(
                self,
                &corner1,
                &corner2,
                ifd1,
                ifd2,
                sens1,
                sens2,
                &vtx,
                st.p1,
                st.p2,
                true,
                true,
            )?;
        }
        // `if (!Okvisavis) { ... }` (`cxx:205-213`).
        if !st.visavis {
            self.done = false;
            return Ok(());
        }
        // `if (!OkinterCC) { ... }` (`cxx:214-222`).
        if !st.ok {
            self.done = false;
            return Ok(());
        }

        // `int IFaArc1 = 3 - IFaCo1, IFaArc2 = 3 - IFaCo2;` (`cxx:223`).
        let i_fa_co1 = st.jf1;
        let i_fa_co2 = st.jf2;
        let i_fa_arc1 = 3 - i_fa_co1;
        let i_fa_arc2 = 3 - i_fa_co2;

        // `const ChFiDS_CommonPoint& CP1 = Fd1->ChangeVertex(isfirst1,
        // IFaArc1);` and `CP2` (`cxx:227-228`).
        let cp1 = fd1.vertex(is_first1, i_fa_arc1).clone();
        let cp2 = fd2.vertex(is_first2, i_fa_arc2).clone();
        // `if (!CP1.IsOnArc() || !CP2.IsOnArc()) { done = false; return done; }`
        // (`cxx:230-237`).
        if !cp1.is_on_arc() || !cp2.is_on_arc() {
            self.done = false;
            return Ok(());
        }
        // `if (!CP1.Arc().IsSame(CP2.Arc())) { done = true;
        // PerformMoreThreeCorner(Index, 2); return done; }` (`cxx:238-247`).
        let arc1 = cp1
            .arc()
            .cloned()
            .ok_or_else(|| "PerformTwoCornerbyInter: CP1 has no arc".to_string())?;
        let arc2 = cp2
            .arc()
            .cloned()
            .ok_or_else(|| "PerformTwoCornerbyInter: CP2 has no arc".to_string())?;
        if !is_same(&arc1, &arc2) {
            self.done = true;
            self.perform_more_three_corner(index, 2)?;
            return Ok(());
        }

        // `TopoDS_Edge pivot = CP1.Arc(); double parCP1 =
        // CP1.ParameterOnArc(); double parCP2 = CP2.ParameterOnArc();`
        // (`cxx:249-252`).
        let pivot = arc1;
        let par_cp1 = cp1
            .parameter_on_arc()
            .ok_or_else(|| "PerformTwoCornerbyInter: CP1 has no parameter on arc".to_string())?;
        let par_cp2 = cp2
            .parameter_on_arc()
            .ok_or_else(|| "PerformTwoCornerbyInter: CP2 has no parameter on arc".to_string())?;

        // `occ::handle<BRepAdaptor_Curve> Hpivot = new BRepAdaptor_Curve(pivot);
        // gp_Pnt psp1 = Hpivot->Value(parCP1); gp_Pnt psp2 =
        // Hpivot->Value(parCP2); double sameparam = (psp1.Distance(psp2) < 10.0
        // * tolapp3d);` (`cxx:253-271`).
        let hpivot = CurveAdaptor::of_edge(&Edge(pivot.clone()));
        let psp1 = hpivot
            .d0(par_cp1)
            .ok_or_else(|| "PerformTwoCornerbyInter: the pivot has no value".to_string())?;
        let psp2 = hpivot
            .d0(par_cp2)
            .ok_or_else(|| "PerformTwoCornerbyInter: the pivot has no value".to_string())?;
        let same_param = psp1.distance(&psp2) < 10.0 * self.tolapp3d;

        // `TopoDS_Face FF1 = TopoDS::Face(DStr.Shape(Fd1->Index(IFaArc1)));
        // TopoDS_Face FF2 = ...;` (`cxx:273-274`).
        let ff1 = face_at(self.my_ds.ds(), fd1.index(i_fa_arc1))?;
        let ff2 = face_at(self.my_ds.ds(), fd2.index(i_fa_arc2))?;
        // `for (Kt.Initialize(myEFMap(pivot)); Kt.More(); Kt.Next()) { TopoDS_Face
        // F = TopoDS::Face(Kt.Value()); if (!ok1 && FF1.IsSame(F)) ok1 = true; if
        // (!ok2 && FF2.IsSame(F)) ok2 = true; }` (`cxx:275-288`).
        let mut ok1 = false;
        let mut ok2 = false;
        let faces: Vec<TopoShape> = self
            .my_ef_map
            .get_from_key(&pivot)
            .to_vec();
        for f in &faces {
            if !ok1 && is_same(f, &ff1.0) {
                ok1 = true;
            }
            if !ok2 && is_same(f, &ff2.0) {
                ok2 = true;
            }
        }
        // `if (!ok1 || !ok2) { done = false; return done; }` (`cxx:289-296`).
        if !ok1 || !ok2 {
            self.done = false;
            return Ok(());
        }

        // `HS1 = ChFi3d_BoundSurf(DStr, Fd1, IFaCo1, IFaArc1); HS2 =
        // ChFi3d_BoundSurf(DStr, Fd2, IFaCo2, IFaArc2);` (`cxx:298-300`).
        let hs1 = {
            let ds = self.my_ds.change_ds();
            ch_fi3d_bound_surf(ds, &fd1, i_fa_co1, i_fa_arc1)?
        };
        let hs2 = {
            let ds = self.my_ds.change_ds();
            ch_fi3d_bound_surf(ds, &fd2, i_fa_co2, i_fa_arc2)?
        };

        let mut tc = TwoCorner {
            vtx,
            pivot,
            par_cp1,
            par_cp2,
            corner1,
            corner2,
            ifd1,
            ifd2,
            fd1,
            fd2,
            is_first1,
            is_first2,
            u_int_pc1: st.p1,
            u_int_pc2: st.p2,
            i_fa_co1,
            i_fa_co2,
            i_fa_arc1,
            i_fa_arc2,
            same_param,
            hs1,
            hs2,
            // `NCollection_Array1<double> Pardeb(1, 4), Parfin(1, 4);`
            // (`cxx:303`): every element the two branches read is written by
            // the branch first.
            pardeb: [0.0; 4],
            parfin: [0.0; 4],
        };

        // `if (sameparam) { ... } else { ... }` (`cxx:306-801`).
        let answer = if tc.same_param {
            same_param::same_param(self, &mut tc)
        } else {
            diff_param::diff_param(self, &mut tc)
        };
        tc.flush();
        match answer {
            // `done = true; return done;` (`cxx:802-806`).
            Ok(true) => {
                self.done = true;
                Ok(())
            }
            // One of the `done = false; return done;` arms of the two branches.
            Ok(false) => {
                self.done = false;
                Ok(())
            }
            Err(e) => {
                self.done = false;
                Err(e)
            }
        }
    }
}
