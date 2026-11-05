//! `ChFi3d_Builder::PerformMoreSurfdata`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_C1.cxx:3771-4531`
//! (OCCT 8.0.0), declared `hxx:784`.
//!
//! The body runs when [`super::ChFi3dBuilder::more_surfdata`] reported that the
//! corner involves more than one `ChFiDS_SurfData` on the stripe. It takes the
//! `SurfData` at the corner (`anInd`) and its neighbour (`anIndPrev`), finds the
//! two arcs where the end face meets the fillet, classifies the four common
//! points into "on an arc" / "not on an arc", and then:
//!
//! - [`sections`] (`cxx:3918-4156`): `ChFi3d_ComputeArete` plus up to two
//!   `GeomInt_IntSS` sections, producing `aCint1`, `aTrCracc` and `aCint2`,
//! - [`store`] (`cxx:4158-4521`): the DS storage — a point, a solid/surface
//!   interference, the removal of `SeqFil(anInd)`, the stripe's point indices,
//!   and the curve / surface / shape interferences of the three curves.
//!
//! Translation notes (each also at the site):
//!
//! - `NCollection_Sequence<Handle(SurfData)>& aSeqSurfData = aStripe->
//!   ChangeSetOfSurfData()->ChangeSequence()` is copied into [`MoreSurfdata::seq`]
//!   and written back by [`MoreSurfdata::flush`]; the body removes one entry
//!   (`cxx:4216`) and has no `return` after that,
//! - `myListStripe.First()` (`cxx:3903`, `:3907`) is read from the builder, not
//!   from the stripe taken out of `myVDataMap(Index)`: the source uses those two
//!   different handles and the port keeps that,
//! - `DStr.Shape(I)` uses the OCCT default `FindKeep = true`,
//! - `ChFi3d_IndexPointInDS` / `ChFi3d_FilPointInDS` / `ChFi3d_FilCurveInDS`
//!   are the port's `builder_0_helpers` functions.

mod sections;
mod store;

use std::sync::Arc;

use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_geom2d::curve::Curve2d;

use crate::chfi_ds::{CommonPoint, StripeHandle, SurfData};
use crate::chfi3d::builder_0_helpers::ch_fi3d_index_of_surf_data;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::is_same;

use super::ChFi3dBuilder;

/// `const double aTol3d = 1.e-4;` — declared in the body (`cxx:3777`).
pub(super) const TOL3D: f64 = 1.0e-4;

/// The locals of `PerformMoreSurfdata` that more than one stage reads.
pub(super) struct MoreSurfdata {
    /// `aStripe` — `myVDataMap(Index).First()`.
    pub stripe: StripeHandle,
    /// `aSeqSurfData` — copied; see the module note.
    pub seq: Vec<Arc<SurfData>>,
    /// `isFirst` (`cxx:3792`).
    pub is_first: bool,
    /// `anInd` (`cxx:3791`), 1-based like the OCCT sequence.
    pub index: usize,
    /// `anIndPrev` (`cxx:3811-3818`), 1-based.
    pub ind_prev: usize,
    /// `aFace` — the `FindFace` answer.
    pub face: Face,
    /// `aSurfPrev` — `BRep_Tool::Surface(aFace)`.
    pub surf_prev: Arc<dyn Surface>,
    /// `aNeighborFace` — `DStr.Shape(indSurface)`.
    pub neighbor_face: Face,
    /// `aCP1onArc`
    pub cp1_on_arc: CommonPoint,
    /// `aCP2onArc`
    pub cp2_on_arc: CommonPoint,
    /// `aCP2NotonArc`
    pub cp2_not_on_arc: CommonPoint,
    /// `is1stCP1OnArc` (`cxx:3883-3896`) — read only by the prologue; kept as
    /// the source's local for the reader of the storage block.
    #[allow(dead_code)]
    pub is1st_cp1_on_arc: bool,
    /// `is2ndCP1OnArc` (`cxx:3863-3878`)
    pub is2nd_cp1_on_arc: bool,
    /// `aCline` (`cxx:3978`)
    pub cline: Arc<dyn Curve>,
    /// `aPClineOnSurf` (`cxx:3981`)
    pub pc_line_on_surf: Option<Arc<dyn Curve2d>>,
    /// `aPClineOnFace` (`cxx:3982`)
    pub pc_line_on_face: Option<Arc<dyn Curve2d>>,
    /// `indLine` (`cxx:3977`)
    pub ind_line: i32,
}

impl MoreSurfdata {
    /// Writes the copied `aSeqSurfData` back into the stripe.
    pub(super) fn flush(&self) {
        if let Ok(mut stripe) = self.stripe.write() {
            *stripe.change_set_of_surf_data() = self.seq.clone();
        }
    }
}

/// `SeqFil(I)`, 1-based (`SeqFil.Value(I)`).
pub(super) fn fd_at(st: &MoreSurfdata, i: usize) -> Result<Arc<SurfData>, String> {
    i.checked_sub(1)
        .and_then(|k| st.seq.get(k))
        .cloned()
        .ok_or_else(|| "PerformMoreSurfdata: SeqFil(I) is out of range".to_string())
}

/// `SeqFil(anInd)`.
pub(super) fn fd_of(st: &MoreSurfdata) -> Result<Arc<SurfData>, String> {
    fd_at(st, st.index)
}

impl ChFi3dBuilder {
    /// `PerformMoreSurfdata(Index)` (`cxx:3771-4531`).
    ///
    /// The source's bare `return`s are `Ok(())` here; the throws (a null
    /// `FindFace` face, a missing p-curve, a null `DStr.Surface`) are `Err`.
    pub fn perform_more_surfdata(&mut self, index: usize) -> Result<(), String> {
        let Some(mut st) = prologue(self, index)? else {
            return Ok(());
        };
        let Some(sec) = sections::compute(self, &st)? else {
            st.flush();
            return Ok(());
        };
        store::write(self, &mut st, &sec)?;
        st.flush();
        Ok(())
    }
}

/// `cxx:3773-3915`: the stripe / sequence / arcs / common-point classification.
///
/// `Ok(None)` is the two bare `return`s of that range: no stripe at the vertex,
/// or no `SurfData` at `anIndPrev` with a common point on either arc.
fn prologue(b: &mut ChFi3dBuilder, index: usize) -> Result<Option<MoreSurfdata>, String> {
    let tol3d = TOL3D;
    // `const NCollection_List<Handle(ChFiDS_Stripe)>& aLOfStripe =
    //  myVDataMap(Index); if (aLOfStripe.IsEmpty()) return;` (`cxx:3775-3781`).
    let key = b.my_v_data_map.find_key(index).cloned();
    let first_stripe = b
        .my_v_data_map
        .find_from_index(index)
        .and_then(|l| l.first().cloned());
    let (Some(vtx), Some(stripe)) = (key, first_stripe) else {
        return Ok(None);
    };

    // `aStripe = aLOfStripe.First(); aSpine = aStripe->Spine();` (`cxx:3783-3784`).
    // `aSeqSurfData = aStripe->ChangeSetOfSurfData()->ChangeSequence();`
    // (`cxx:3786-3787`).
    let seq: Vec<Arc<SurfData>> = {
        let Ok(guard) = stripe.read() else {
            return Ok(None);
        };
        guard.set_of_surf_data().to_vec()
    };

    // `aSens = 0; anInd = ChFi3d_IndexOfSurfData(aVtx, aStripe, aSens);
    //  isFirst = (aSens == 1);` (`cxx:3790-3792`).
    let mut sens = 0;
    let Ok(an_ind) = ch_fi3d_index_of_surf_data(&vtx, &stripe, &mut sens) else {
        return Ok(None);
    };
    let Ok(an_ind) = usize::try_from(an_ind) else {
        return Ok(None);
    };
    let is_first = sens == 1;

    // `aSurfData = aSeqSurfData.Value(anInd); aCP1 = aSurfData->Vertex(isFirst, 1);
    //  aCP2 = aSurfData->Vertex(isFirst, 2);` (`cxx:3798-3801`).
    let fd = an_ind
        .checked_sub(1)
        .and_then(|k| seq.get(k))
        .cloned()
        .ok_or_else(|| "PerformMoreSurfdata: SeqFil(anInd) is out of range".to_string())?;
    let cp1 = fd.vertex(is_first, 1).clone();
    let cp2 = fd.vertex(is_first, 2).clone();

    // `FindFace(aVtx, aCP1, aCP2, aFace);` (`cxx:3808`), then
    // `aSurfPrev = BRep_Tool::Surface(aFace);` (`cxx:3809`) — both raise in OCCT
    // when `aFace` is null.
    let face = b
        .find_face(&vtx, &cp1, &cp2, None)
        .ok_or_else(|| "PerformMoreSurfdata: FindFace returned no face".to_string())?;
    let surf_prev = crate::brep_tool::BRepTool::face_surface(&face)
        .ok_or_else(|| "PerformMoreSurfdata: BRep_Tool::Surface(aFace) is null".to_string())?;

    // `if (aSens == 1) anIndPrev = anInd + 1; else anIndPrev = anInd - 1;`
    // (`cxx:3811-3818`).
    let ind_prev = if sens == 1 {
        an_ind + 1
    } else {
        an_ind
            .checked_sub(1)
            .ok_or_else(|| "PerformMoreSurfdata: anIndPrev is out of range".to_string())?
    };

    // `anArc1` / `anArc2`: `cxx:3823-3845`. The source assigns the candidate
    // before testing it, so on a miss the handle holds the last map entry.
    let mut arc1: Option<Edge> = None;
    for item in b.my_ve_map.get_from_key(&vtx) {
        let e = Edge::wrap(item.clone())
            .ok_or_else(|| "PerformMoreSurfdata: TopoDS::Edge on myVEMap(aVtx)".to_string())?;
        let hit = crate::chfi3d::builder_c1::contain_e(&face, &e);
        arc1 = Some(e);
        if hit {
            break;
        }
    }
    let mut arc2: Option<Edge> = None;
    for item in b.my_ve_map.get_from_key(&vtx) {
        let e = Edge::wrap(item.clone())
            .ok_or_else(|| "PerformMoreSurfdata: TopoDS::Edge on myVEMap(aVtx)".to_string())?;
        let hit = crate::chfi3d::builder_c1::contain_e(&face, &e)
            && !arc1.as_ref().is_some_and(|a| is_same(&a.0, &e.0));
        arc2 = Some(e);
        if hit {
            break;
        }
    }

    // `aSurfData = aSeqSurfData.Value(anIndPrev); aCP1 = ...Vertex(isFirst, 1);
    //  aCP2 = ...Vertex(isFirst, 2);` (`cxx:3859-3861`).
    let fd_prev = ind_prev
        .checked_sub(1)
        .and_then(|k| seq.get(k))
        .cloned()
        .ok_or_else(|| "PerformMoreSurfdata: SeqFil(anIndPrev) is out of range".to_string())?;
    let cp1_prev = fd_prev.vertex(is_first, 1).clone();
    let cp2_prev = fd_prev.vertex(is_first, 2).clone();

    // `cxx:3863-3878`: which of the two `anIndPrev` points is on an arc.
    let arc_is = |candidate: Option<&TopoShape>| -> bool {
        let Some(c) = candidate else {
            return false;
        };
        arc1.as_ref().is_some_and(|a| is_same(c, &a.0))
            || arc2.as_ref().is_some_and(|a| is_same(c, &a.0))
    };
    let cp1_prev_on_arc = cp1_prev.is_on_arc() && arc_is(cp1_prev.arc());
    let cp2_prev_on_arc = cp2_prev.is_on_arc() && arc_is(cp2_prev.arc());
    let (cp2_on_arc, cp2_not_on_arc, is2nd_cp1_on_arc) = if cp1_prev_on_arc {
        (cp1_prev, cp2_prev, true)
    } else if cp2_prev_on_arc {
        (cp2_prev, cp1_prev, false)
    } else {
        return Ok(None);
    };

    // `aSurfData = aSeqSurfData.Value(anInd); aCP1 / aCP2 = Vertex(isFirst, 1/2);`
    // (`cxx:3880-3882`).
    let cp1 = fd.vertex(is_first, 1).clone();
    let cp2 = fd.vertex(is_first, 2).clone();

    // `cxx:3884-3898`: the point of `anInd` that is close to `aCP2onArc`.
    let (cp1_on_arc, is1st_cp1_on_arc) = if cp1.point().distance(cp2_on_arc.point()) <= tol3d {
        (cp2, false)
    } else {
        (cp1, true)
    };
    if !cp1_on_arc.is_on_arc() {
        return Ok(None);
    }

    // `cxx:3900-3910`: `indSurface` from `myListStripe.First()`.
    let current = b
        .my_list_stripe
        .first()
        .cloned()
        .ok_or_else(|| "PerformMoreSurfdata: myListStripe is empty".to_string())?;
    let ind_surface = {
        let Ok(guard) = current.read() else {
            return Err("PerformMoreSurfdata: myListStripe.First() is locked".to_string());
        };
        let seq_cur = guard.set_of_surf_data();
        let fd_cur = an_ind
            .checked_sub(1)
            .and_then(|k| seq_cur.get(k))
            .ok_or_else(|| "PerformMoreSurfdata: SeqFil(anInd) is out of range".to_string())?;
        if is1st_cp1_on_arc {
            fd_cur.index_of_s1()
        } else {
            fd_cur.index_of_s2()
        }
    };
    // `aNeighborFace = TopoDS::Face(myDS->Shape(indSurface));` (`cxx:3910`).
    let neighbor_face = b
        .my_ds
        .ds()
        .shape(ind_surface, true)
        .and_then(Face::wrap)
        .ok_or_else(|| "PerformMoreSurfdata: DStr.Shape(indSurface) is not a face".to_string())?;

    // `cxx:3959-3982`: the line on `anInd` and its two p-curves.
    let fi = if is1st_cp1_on_arc {
        fd.interference_on_s1().clone()
    } else {
        fd.interference_on_s2().clone()
    };
    let ind_line = fi.line_index();
    let cline = b
        .my_ds
        .ds()
        .curve(ind_line)
        .curve()
        .cloned()
        .ok_or_else(|| "PerformMoreSurfdata: DStr.Curve(indLine) is null".to_string())?;
    let pc_line_on_surf = fi.p_curve_on_surf().cloned();
    let pc_line_on_face = fi.p_curve_on_face().cloned();

    Ok(Some(MoreSurfdata {
        stripe,
        seq,
        is_first,
        index: an_ind,
        ind_prev,
        face,
        surf_prev,
        neighbor_face,
        cp1_on_arc,
        cp2_on_arc,
        cp2_not_on_arc,
        is1st_cp1_on_arc,
        is2nd_cp1_on_arc,
        cline,
        pc_line_on_surf,
        pc_line_on_face,
        ind_line,
    }))
}
