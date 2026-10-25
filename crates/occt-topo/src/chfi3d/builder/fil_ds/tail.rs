//! The tail of one `j` iteration of `ChFi3d_FilDS`: the last curve of the
//! stripe, the fillet/fillet curve in the middle of it, the "queue de billard"
//! degenerated end and the trailing point/edge interferences.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:3160-3339`
//! (OCCT 8.0.0):
//!
//! - `if (j == SeqFil.Length()) { if (!isInDS2) { ... Icurv =
//!   CorDat->LastCurve(); ... } }` (`cxx:3160-3206`),
//! - `else { if (!isInDS2) { ChFi3d_ComputeArete(...); Icurv =
//!   DStr.AddCurve(Crv); ... } }` (`cxx:3208-3240`),
//! - the degenerated "queue de billard" end (`cxx:3244-3300`) using the
//!   `NumEdge` / `BoutdeVtx` state of the stripe scan,
//! - the point/edge interferences of whatever is left (`cxx:3302-3339`).
//!
//! Translation notes:
//! - `reglist.Append(regcout)` / `reglist.Append(regfilfil)` copy the
//!   `ChFiDS_Regul` (`NCollection_List::Append(const TheItemType&)`), so the
//!   `regfilfil.SetCurve` / `SetS1` of `cxx:3232-3233` only updates the
//!   builder-local `ChFiDS_Regul` and not the copy appended by `cxx:2936`;
//!   [`FilDsState::regfilfil`] keeps the same builder-local role,
//! - `ChFi3d_ComputeArete(..., IFlag = 0)` (`cxx:3214-3228`) writes `C3d`,
//!   `PCurv`, `Pardeb`, `Parfin` and `tolreached`; [`ComputeArete`] carries
//!   them, with the incoming `PCurv` / `Pardeb` / `Parfin` of the previous
//!   iteration like the OCCT locals,
//! - `TopOpeBRepDS_Curve(C3d, tolreached)` is
//!   [`Curve::with_curve`] with `IsWalk = false` (the two-argument OCCT
//!   constructor of `TopOpeBRepDS_Curve.cxx` does not set the walk flag),
//! - the `while (NumEdge <= spine->NbEdges() && !Trouve)` scan of
//!   `cxx:3260-3285` is translated with its exact loop exits: the inner vertex
//!   scan stops as soon as `Trouve` is set, `Arcspine` keeps the last edge
//!   tried, and `NumEdge` is only incremented when the edge did not answer,
//! - OCCT would dereference a null `TopoDS_Edge` on a spine with no edge; the
//!   port reports that as an `Err`.

use std::sync::Arc;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::chfi_ds::{FaceInterference, StripeHandle, SurfData};
use crate::chfi3d::builder_0_helpers::compute_arete::{
    ch_fi3d_compute_arete, ComputeArete,
};
use crate::chfi3d::builder_0_helpers::fil_interf::{
    ch_fi3d_contains, ch_fi3d_fil_curve_in_ds, ch_fi3d_fil_point_in_ds, ch_fi3d_fil_vertex_in_ds,
};
use crate::iterator::cumulated_children;
use crate::topo_tools::is_same;
use crate::topopebrep::ds::Curve;

use super::super::spine_of;
use super::{p_curve_value, with_stripe, FilDsState};
use crate::chfi3d::builder::ChFi3dBuilder;

impl ChFi3dBuilder {
    /// `if (j == SeqFil.Length()) { if (!isInDS2) { ... } }` (`cxx:3160-3206`).
    pub(super) fn last_curve(
        &mut self,
        st: &mut FilDsState,
        cor_dat: &StripeHandle,
        isurf: i32,
        et1: &mut Orientation,
    ) -> Result<(), String> {
        // `Icurv = CorDat->LastCurve();` (`cxx:3165`).
        st.icurv = with_stripe(cor_dat, |s| s.last_curve())?;

        // `if (Closed && !Singulier_en_Bout && (Ipoin1 != Ipoin2)) {
        // regcout.SetS2(Isurf, false); reglist.Append(regcout); }`
        // (`cxx:3166-3170`).
        if st.closed && !st.singulier_en_bout && st.ipo_in1 != st.ipo_in2 {
            st.regcout.set_s2(isurf, false);
            self.append_regul(st.regcout);
        }

        // `PCurv = CorDat->LastPCurve(); ET1 = TopAbs::Compose(ET1,
        // CorDat->LastPCurveOrientation()); CorDat->LastParameters(Pardeb,
        // Parfin);` (`cxx:3171-3174`).
        let last_ori = with_stripe(cor_dat, |s| s.last_p_curve_orientation())?;
        st.p_curv = with_stripe(cor_dat, |s| s.last_p_curve().cloned())?;
        *et1 = Orientation::compose(*et1, last_ori);
        let (pardeb, parfin) = with_stripe(cor_dat, |s| s.last_parameters())?;
        st.pardeb = pardeb;
        st.parfin = parfin;

        // `if (Li.IsEmpty()) { if (LastPCurveOrientation() == TopAbs_REVERSED)
        // {...} else {...} Li.Append(Interfp5); Li.Append(Interfp6); }`
        // (`cxx:3176-3190`).
        if self.ds().curve_interferences(st.icurv).is_empty() {
            let (o1, p1, o2, p2) = if last_ori == Orientation::Reversed {
                (
                    Orientation::Reversed,
                    st.parfin,
                    Orientation::Forward,
                    st.pardeb,
                )
            } else {
                (
                    Orientation::Forward,
                    st.pardeb,
                    Orientation::Reversed,
                    st.parfin,
                )
            };
            let interfp5 = ch_fi3d_fil_point_in_ds(o1, st.icurv, st.ipo_in1, p1, st.is_vertex1);
            let interfp6 = ch_fi3d_fil_point_in_ds(o2, st.icurv, st.ipo_in2, p2, st.is_vertex2);
            let li = self.my_ds.change_ds().change_curve_interferences(st.icurv);
            li.push(interfp5);
            li.push(interfp6);
        }

        // `Interfc1 = ChFi3d_FilCurveInDS(Icurv, Isurf, PCurv, ET1);
        // DStr.ChangeSurfaceInterferences(Isurf).Append(Interfc1);`
        // (`cxx:3191-3192`).
        let interfc1 = ch_fi3d_fil_curve_in_ds(st.icurv, isurf, st.p_curv.clone(), *et1);
        self.my_ds
            .change_ds()
            .change_surface_interferences(isurf)
            .push(interfc1.clone());

        // `if (Ipoin1 == Ipoin2) { ... TCurv.ChangeCurve().Nullify();
        // TCurv.SetSCI(Interfc1, bidinterf); }` (`cxx:3193-3199`).
        if st.ipo_in1 == st.ipo_in2 {
            let tcurv = self.my_ds.change_ds().change_curve(st.icurv);
            *tcurv.change_curve() = None;
            tcurv.set_sci(Some(interfc1), None);
        }
        Ok(())
    }

    /// `else { if (!isInDS2) { ... } }` (`cxx:3208-3240`): the fillet/fillet
    /// curve between this `SurfData` and the next one.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn middle_curve(
        &mut self,
        st: &mut FilDsState,
        fd: &SurfData,
        fi1: &FaceInterference,
        fi2: &FaceInterference,
        isurf: i32,
        et1: Orientation,
        tol3d: f64,
        tol2d: f64,
    ) -> Result<(), String> {
        // `Fd->InterferenceOnS1().PCurveOnSurf()->Value(Fd->
        // InterferenceOnS1().LastParameter())` and the S2 companion
        // (`cxx:3216-3219`).
        let uv1 = p_curve_value(
            fi1.p_curve_on_surf(),
            fi1.last_parameter(),
            "the S1 interference",
        )?;
        let uv2 = p_curve_value(
            fi2.p_curve_on_surf(),
            fi2.last_parameter(),
            "the S2 interference",
        )?;

        // `DStr.Surface(Fd->Surf()).Surface()` (`cxx:3220`).
        let surf = self
            .ds()
            .surface(fd.surf())
            .surface()
            .cloned()
            .ok_or_else(|| "ChFi3d_FilDS: the fillet surface is not in the DS".to_string())?;

        // `ChFi3d_ComputeArete(V3, UV1, V4, UV2, Surf, C3d, PCurv, Pardeb,
        // Parfin, tol3d, tol2d, tolreached, 0);` (`cxx:3214-3228`).
        let mut arete = ComputeArete {
            c3d: None,
            p_curve: st.p_curv.clone(),
            pardeb: st.pardeb,
            parfin: st.parfin,
            tolreached: 0.0,
        };
        ch_fi3d_compute_arete(
            &st.v3,
            &uv1,
            &st.v4,
            &uv2,
            &surf,
            &mut arete,
            tol3d,
            tol2d,
            0,
        )?;

        // `Crv = TopOpeBRepDS_Curve(C3d, tolreached); Icurv =
        // DStr.AddCurve(Crv);` (`cxx:3229-3230`).
        let c3d = arete
            .c3d
            .clone()
            .ok_or_else(|| "ChFi3d_FilDS: ChFi3d_ComputeArete built no 3D curve".to_string())?;
        let mut crv = Curve::with_curve(Some(c3d), arete.tolreached, false);
        st.icurv = self.my_ds.change_ds().add_curve(&mut crv);
        st.p_curv = arete.p_curve.clone();
        st.pardeb = arete.pardeb;
        st.parfin = arete.parfin;

        // `regfilfil.SetCurve(Icurv); regfilfil.SetS1(Isurf, false);`
        // (`cxx:3231-3233`).
        st.regfilfil.set_curve(st.icurv);
        st.regfilfil.set_s1(isurf, false);

        // `Interfp5 = ChFi3d_FilPointInDS(TopAbs_FORWARD, Icurv, Ipoin1,
        // Pardeb, isVertex1); ...Append(Interfp5); Interfp6 =
        // ChFi3d_FilPointInDS(TopAbs_REVERSED, Icurv, Ipoin2, Parfin,
        // isVertex2); ...Append(Interfp6);` (`cxx:3234-3237`).
        let interfp5 =
            ch_fi3d_fil_point_in_ds(Orientation::Forward, st.icurv, st.ipo_in1, st.pardeb, st.is_vertex1);
        let interfp6 =
            ch_fi3d_fil_point_in_ds(Orientation::Reversed, st.icurv, st.ipo_in2, st.parfin, st.is_vertex2);
        {
            let li = self.my_ds.change_ds().change_curve_interferences(st.icurv);
            li.push(interfp5);
            li.push(interfp6);
        }

        // `Interfc1 = ChFi3d_FilCurveInDS(Icurv, Isurf, PCurv, ET1);
        // DStr.ChangeSurfaceInterferences(Isurf).Append(Interfc1);`
        // (`cxx:3238-3239`).
        let interfc1 = ch_fi3d_fil_curve_in_ds(st.icurv, isurf, st.p_curv.clone(), et1);
        self.my_ds
            .change_ds()
            .change_surface_interferences(isurf)
            .push(interfc1);
        Ok(())
    }

    /// `if (Degene) { ... }` (`cxx:3244-3300`): the "queue de billard" end.
    pub(super) fn degenerated_end(
        &mut self,
        st: &mut FilDsState,
        cor_dat: &StripeHandle,
        seq_fil: &[Arc<SurfData>],
        j: usize,
    ) -> Result<(), String> {
        // `bool Vertex = (V3.IsVertex()) && (V4.IsVertex()); if (!Vertex) {}
        // else { ... }` (`cxx:3247-3251`).
        if !(st.v3.is_vertex() && st.v4.is_vertex()) {
            return Ok(());
        }

        let spine = spine_of(cor_dat)
            .ok_or_else(|| "ChFi3d_FilDS: the stripe has no spine".to_string())?;
        let nb_edges = {
            let guard = spine
                .read()
                .map_err(|_| "ChFi3d_FilDS: poisoned spine".to_string())?;
            guard.nb_edges()
        };
        // `BoutdeVtx = V3.Vertex();` (`cxx:3258`).
        let bout = st
            .v3
            .vertex()
            .cloned()
            .ok_or_else(|| "ChFi3d_FilDS: V3 is a vertex without a shape".to_string())?;
        st.bout_de_vtx = Some(bout.clone());

        // `bool Trouve = false; TopoDS_Edge Arcspine; TopAbs_Orientation OVtx =
        // TopAbs_FORWARD; while (NumEdge <= spine->NbEdges() && !Trouve)
        // { Arcspine = spine->Edges(NumEdge); for (ex.Init(Arcspine.Oriented(
        // TopAbs_FORWARD), TopAbs_VERTEX); ex.More() && (!Trouve); ex.Next())
        // { if (BoutdeVtx.IsSame(ex.Current())) { OVtx = ex.Current().
        // Orientation(); if (Closed && (NumEdge == 1)) { Trouve =
        // (spine->NbEdges() == 1); } else { Trouve = true; } } } if (!Trouve)
        // { NumEdge++; } }` (`cxx:3255-3285`).
        let mut trouve = false;
        let mut arcspine: Option<crate::shape::TopoShape> = None;
        let mut ovtx = Orientation::Forward;
        while st.num_edge <= nb_edges && !trouve {
            let e = {
                let guard = spine
                    .read()
                    .map_err(|_| "ChFi3d_FilDS: poisoned spine".to_string())?;
                guard.edges(st.num_edge).cloned().ok_or_else(|| {
                    format!("ChFi3d_FilDS: spine->Edges({}) is out of range", st.num_edge)
                })?
            };
            for vv in cumulated_children(&e.oriented(Orientation::Forward)) {
                if vv.shape_type() != ShapeType::Vertex {
                    continue;
                }
                if is_same(&vv, &bout) {
                    ovtx = vv.orientation();
                    if st.closed && st.num_edge == 1 {
                        trouve = nb_edges == 1;
                    } else {
                        trouve = true;
                    }
                    if trouve {
                        break;
                    }
                }
            }
            arcspine = Some(e);
            if !trouve {
                st.num_edge += 1; // Go to the next edge
            }
        }
        let arcspine = arcspine.ok_or_else(|| {
            "ChFi3d_FilDS: the spine has no edge to remove at the degenerated end".to_string()
        })?;

        // `int IArcspine = DStr.AddShape(Arcspine); int IVtx; if (j ==
        // SeqFil.Length()) { IVtx = CorDat->IndexLastPointOnS1(); } else {
        // IVtx = DStr.AddShape(BoutdeVtx); }` (`cxx:3286-3293`).
        let i_arcspine = self.my_ds.change_ds().add_shape(&arcspine);
        let ivtx = if j == seq_fil.len() {
            with_stripe(cor_dat, |s| s.index_last_point_on_s1())?
        } else {
            self.my_ds.change_ds().add_shape(&bout)
        };

        // `OVtx = TopAbs::Reverse(OVtx); double parVtx =
        // BRep_Tool::Parameter(BoutdeVtx, Arcspine); interfv =
        // ChFi3d_FilVertexInDS(OVtx, IArcspine, IVtx, parVtx);
        // DStr.ChangeShapeInterferences(IArcspine).Append(interfv);`
        // (`cxx:3294-3298`).
        let ovtx = ovtx.reversed();
        let par_vtx = BRepTool::parameter_of_vertex(&arcspine, &bout);
        let interfv = ch_fi3d_fil_vertex_in_ds(ovtx, i_arcspine, ivtx, par_vtx);
        self.my_ds
            .change_ds()
            .change_shape_interferences_at(i_arcspine)
            .push(interfv);
        Ok(())
    }

    /// One "Processing of interference Point / Edges" block
    /// (`cxx:3302-3339`), for `side == 1` (`V3`) or `side == 2` (`V4`).
    pub(super) fn point_edge_interference(
        &mut self,
        st: &mut FilDsState,
        fd: &SurfData,
        side: i32,
    ) -> Result<(), String> {
        let v = if side == 1 {
            st.v3.clone()
        } else {
            st.v4.clone()
        };
        let ipoin = if side == 1 { st.ipo_in1 } else { st.ipo_in2 };

        // `if (V3.IsOnArc()) { if (!(V3.IsVertex() && Fd->IsOnCurve1())) {...}
        // }` (`cxx:3303-3320`).
        if !v.is_on_arc() {
            return Ok(());
        }
        let on_curve = if side == 1 {
            fd.is_on_curve1()
        } else {
            fd.is_on_curve2()
        };
        if v.is_vertex() && on_curve {
            return Ok(());
        }

        // `Iarc1 = DStr.AddShape(V3.Arc());` (`cxx:3307`).
        let arc = v
            .arc()
            .cloned()
            .ok_or_else(|| "ChFi3d_FilDS: the common point is on an arc without a shape".to_string())?;
        let iarc = self.my_ds.change_ds().add_shape(&arc);
        if side == 1 {
            st.iarc1 = iarc;
        } else {
            st.iarc2 = iarc;
        }

        // `if (!ChFi3d_Contains(DStr.ShapeInterferences(Iarc1), Iarc1, Ipoin1,
        // V3.IsVertex(), true)) { Interfpp = ChFi3d_FilPointInDS(V3.
        // TransitionOnArc(), Iarc1, Ipoin1, V3.ParameterOnArc(),
        // V3.IsVertex()); DStr.ChangeShapeInterferences(V3.Arc()).Append(
        // Interfpp); }` (`cxx:3308-3318`).
        if !ch_fi3d_contains(
            self.ds().shape_interferences_at(iarc, false),
            iarc,
            ipoin,
            v.is_vertex(),
            true,
        ) {
            let interfpp = ch_fi3d_fil_point_in_ds(
                v.transition_on_arc().ok_or_else(|| {
                    "ChFi3d_FilDS: the common point is on an arc without a transition".to_string()
                })?,
                iarc,
                ipoin,
                v.parameter_on_arc().ok_or_else(|| {
                    "ChFi3d_FilDS: the common point is on an arc without a parameter".to_string()
                })?,
                v.is_vertex(),
            );
            self.my_ds
                .change_ds()
                .change_shape_interferences(&arc)
                .push(interfpp);
        }
        Ok(())
    }
}
