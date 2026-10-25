//! The `for (j = 1; j <= SeqFil.Length(); j++)` body of `ChFi3d_FilDS`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:2655-3340`
//! (OCCT 8.0.0). Section anchors, in the order the body runs them:
//!
//! - `SeqFil(j)` / `Isurf` / `Ishape1` / `Ishape2` (`cxx:2658-2661`),
//! - the `isInDS1` / `isInDS2` end flags (`cxx:2663-2676`),
//! - the `SolidSurfaceInterference` (`cxx:2678-2687`),
//! - the "double interferences" of two consecutive `SurfData` sharing an arc
//!   (`cxx:2694-2720`),
//! - `V3` / `V4` and the `trafil1` / `trafil2` / `ET1` orientations
//!   (`cxx:2722-2754`),
//! - the four `CutEdge` calls (`cxx:2756-2777`),
//! - the initial processing of `j == 1` (`cxx:2779-2887`),
//! - the inter-fillet interference of `j > 1` (`cxx:2888-2952`),
//! - the fillet/faces interferences on S1 (`cxx:2954-3054`) and S2
//!   (`cxx:3056-3158`),
//! - the last curve of the stripe or, in the middle of it, the
//!   `ChFi3d_ComputeArete` fillet/fillet curve (`cxx:3160-3240`),
//! - the degenerated end ("queue de billard") and the point/edge
//!   interferences (`cxx:3242-3339`).
//!
//! Translation notes:
//! - `DStr.ChangeShapeInterferences(S)` / `ChangeCurveInterferences(I)` /
//!   `ChangeSurfaceInterferences(I)` are the `*_at` / shape-keyed mutators of
//!   [`crate::topopebrep::ds::DataStructure`]; every list that OCCT keeps live
//!   across a later append is re-read by index here, so the `ChFi3d_Contains`
//!   tests still see the interferences appended earlier in the same block,
//! - `Fd` is read-only in the whole body (it is `const` in `cxx:2658`), so the
//!   `SurfData` is only cloned out of the sequence,
//! - the `Interfp*` / `Interfc*` handles of `cxx:2593-2597` are single-use
//!   temporaries here; only `Interfc1`/`Interfc3`/`Interfc4` are also handed to
//!   `SetSCI`, and that clone is explicit,
//! - `TopOpeBRepDS_Curve& TCurv = DStr.ChangeCurve(Icurv); ... ChangeCurve()
//!   .Nullify(); TCurv.SetSCI(...)` is the sequence
//!   [`crate::topopebrep::ds::Curve::change_curve`] then
//!   [`crate::topopebrep::ds::Curve::set_sci`]; the `bidinterf` argument of
//!   `SetSCI` is a null handle, which the port spells `None`.

use std::sync::{Arc, RwLock};

use crate::abs::Orientation;
use crate::chfi_ds::{StripeHandle, SurfData};
use crate::chfi3d::builder_0_helpers::fil_interf::{
    ch_fi3d_contains, ch_fi3d_fil_point_in_ds, ch_fi3d_orientation, cut_edge,
};
use crate::topo_tools::is_same;
use crate::topopebrep::ds::{Interference, InterferenceHandle, Kind, Transition};

use super::{same_point, seq_at, with_stripe, FilDsState};
use crate::chfi3d::builder::ChFi3dBuilder;

impl ChFi3dBuilder {
    /// One `j` cell of `ChFi3d_FilDS` (`cxx:2655-3340`).
    #[allow(clippy::too_many_arguments)]
    pub(super) fn fil_ds_step(
        &mut self,
        st: &mut FilDsState,
        seq_fil: &[Arc<SurfData>],
        j: usize,
        cor_dat: &StripeHandle,
        solid_index: i32,
        tol3d: f64,
        tol2d: f64,
    ) -> Result<(), String> {
        let len = seq_fil.len();
        let fd = seq_at(seq_fil, j)?;

        // `Isurf = Fd->Surf(); Ishape1 = Fd->IndexOfS1(); Ishape2 =
        // Fd->IndexOfS2();` (`cxx:2659-2661`).
        let isurf = fd.surf();
        let ishape1 = fd.index_of_s1();
        let ishape2 = fd.index_of_s2();

        // `bool isInDS1 = false, isInDS2 = false; if (j <=
        // CorDat->IsInDS(true)) {...} if (SeqFil.Length() - j <
        // CorDat->IsInDS(false)) {...}` (`cxx:2666-2676`).
        let (in_ds_first, in_ds_last) =
            with_stripe(cor_dat, |s| (s.is_in_ds(true), s.is_in_ds(false)))?;
        let mut is_in_ds1 = false;
        let mut is_in_ds2 = false;
        let ji = j as i32;
        let ni = len as i32;
        if ji <= in_ds_first {
            is_in_ds1 = true;
            is_in_ds2 = ji + 1 <= in_ds_first;
        }
        if ni - ji < in_ds_last {
            is_in_ds2 = true;
            is_in_ds1 = is_in_ds1 || ni - ji + 1 < in_ds_last;
        }

        // `new TopOpeBRepDS_SolidSurfaceInterference(Transition(Fd->Orientation()),
        // SOLID, SolidIndex, SURFACE, Isurf); SolidInterfs.Append(SSI);`
        // (`cxx:2680-2687`).
        let ssi: InterferenceHandle = Arc::new(RwLock::new(Interference::solid_surface(
            Transition::from_orientation(fd.orientation()),
            Kind::Solid,
            solid_index,
            Kind::Surface,
            isurf,
        )));
        self.push_solid_interference(solid_index, ssi);

        // `const ChFiDS_FaceInterference& Fi1 = Fd->InterferenceOnS1();` and
        // the three companions (`cxx:2689-2692`). Copied out because the
        // `SurfData` is shared.
        let fi1 = fd.interference_on_s1().clone();
        let fi2 = fd.interference_on_s2().clone();
        let v1 = fd.vertex_first_on_s1().clone();
        let v2 = fd.vertex_first_on_s2().clone();

        // ---- Processing to manage double interferences (`cxx:2694-2720`).
        if j > 1 {
            if v1.is_on_arc() && st.v3.is_on_arc() {
                let a = v1.arc().cloned();
                let b = st.v3.arc().cloned();
                if let (Some(a), Some(b)) = (a, b) {
                    if is_same(&a, &b) {
                        // `ChFi3d_Contains(DStr.ShapeInterferences(Iarc1), Iarc1,
                        // Ipoin1) && (V1.TransitionOnArc() !=
                        // V3.TransitionOnArc())` (`cxx:2701-2702`).
                        let contains = ch_fi3d_contains(
                            self.ds().shape_interferences_at(st.iarc1, false),
                            st.iarc1,
                            st.ipo_in1,
                            false,
                            false,
                        );
                        let transition_differs =
                            v1.transition_on_arc() != st.v3.transition_on_arc();
                        if contains && transition_differs {
                            // `ChFi3d_FilPointInDS(V1.TransitionOnArc(), Iarc1,
                            // Ipoin1, V1.ParameterOnArc())` (`cxx:2704`).
                            let interfp1 = ch_fi3d_fil_point_in_ds(
                                v1.transition_on_arc().ok_or_else(|| {
                                    "ChFi3d_FilDS: V1 is on an arc without a transition"
                                        .to_string()
                                })?,
                                st.iarc1,
                                st.ipo_in1,
                                v1.parameter_on_arc().ok_or_else(|| {
                                    "ChFi3d_FilDS: V1 is on an arc without a parameter".to_string()
                                })?,
                                false,
                            );
                            // `DStr.ChangeShapeInterferences(V1.Arc()).Append(...)`
                            // (`cxx:2705`).
                            self.my_ds
                                .change_ds()
                                .change_shape_interferences(&a)
                                .push(interfp1);
                        }
                    }
                }
            }

            if v2.is_on_arc() && st.v4.is_on_arc() {
                let a = v2.arc().cloned();
                let b = st.v4.arc().cloned();
                if let (Some(a), Some(b)) = (a, b) {
                    if is_same(&a, &b) {
                        // `cxx:2713-2714`.
                        let contains = ch_fi3d_contains(
                            self.ds().shape_interferences_at(st.iarc2, false),
                            st.iarc2,
                            st.ipo_in2,
                            false,
                            false,
                        );
                        let transition_differs =
                            v2.transition_on_arc() != st.v4.transition_on_arc();
                        if contains && transition_differs {
                            let interfp2 = ch_fi3d_fil_point_in_ds(
                                v2.transition_on_arc().ok_or_else(|| {
                                    "ChFi3d_FilDS: V2 is on an arc without a transition"
                                        .to_string()
                                })?,
                                st.iarc2,
                                st.ipo_in2,
                                v2.parameter_on_arc().ok_or_else(|| {
                                    "ChFi3d_FilDS: V2 is on an arc without a parameter".to_string()
                                })?,
                                false,
                            );
                            // `cxx:2717`.
                            self.my_ds
                                .change_ds()
                                .change_shape_interferences(&a)
                                .push(interfp2);
                        }
                    }
                }
            }
        }

        // `V3 = Fd->VertexLastOnS1(); V4 = Fd->VertexLastOnS2();`
        // (`cxx:2722-2723`).
        st.v3 = fd.vertex_last_on_s1().clone();
        st.v4 = fd.vertex_last_on_s2().clone();

        // `if (Ishape1 != 0) {...} else {...}` (`cxx:2725-2752`).
        if ishape1 != 0 {
            if ishape1 > 0 {
                // `trafil1 = DStr.Shape(Ishape1).Orientation();` (`cxx:2729`).
                st.trafil1 = self.shape_orientation(ishape1)?;
            } else if let Some(or) =
                ch_fi3d_orientation(self.solid_interferences(solid_index), solid_index, -ishape1, false, false)
            {
                // `ChFi3d_Orientation(SolidInterfs, SolidIndex, -Ishape1,
                // trafil1);` (`cxx:2733`).
                st.trafil1 = or;
            }
            // `trafil1 = TopAbs::Compose(trafil1, Fd->Orientation()); trafil1 =
            // TopAbs::Compose(TopAbs::Reverse(Fi1.Transition()), trafil1);
            // trafil2 = TopAbs::Reverse(trafil1);` (`cxx:2735-2737`).
            st.trafil1 = Orientation::compose(st.trafil1, fd.orientation());
            st.trafil1 = Orientation::compose(fi1.transition().reversed(), st.trafil1);
            st.trafil2 = st.trafil1.reversed();
        } else {
            if ishape2 > 0 {
                st.trafil2 = self.shape_orientation(ishape2)?;
            } else if let Some(or) =
                ch_fi3d_orientation(self.solid_interferences(solid_index), solid_index, -ishape2, false, false)
            {
                st.trafil2 = or;
            }
            st.trafil2 = Orientation::compose(st.trafil2, fd.orientation());
            st.trafil2 = Orientation::compose(fi2.transition().reversed(), st.trafil2);
            st.trafil1 = st.trafil2.reversed();
        }

        // `ET1 = TopAbs::Reverse(trafil1);` (`cxx:2754`).
        let mut et1 = st.trafil1.reversed();

        // `if (V1.IsVertex() && Fd->IsOnCurve1()) { CutEdge(V1.Vertex(), Fd,
        // DStr, true, 1); }` and the three companions (`cxx:2758-2777`).
        if v1.is_vertex() && fd.is_on_curve1() {
            let vv1 = v1
                .vertex()
                .cloned()
                .ok_or_else(|| "ChFi3d_FilDS: V1 is a vertex without a shape".to_string())?;
            cut_edge(&vv1, &fd, self.my_ds.change_ds(), 1)?;
        }
        if v2.is_vertex() && fd.is_on_curve2() {
            let vv2 = v2
                .vertex()
                .cloned()
                .ok_or_else(|| "ChFi3d_FilDS: V2 is a vertex without a shape".to_string())?;
            cut_edge(&vv2, &fd, self.my_ds.change_ds(), 2)?;
        }
        if st.v3.is_vertex() && fd.is_on_curve1() {
            let vv3 = st
                .v3
                .vertex()
                .cloned()
                .ok_or_else(|| "ChFi3d_FilDS: V3 is a vertex without a shape".to_string())?;
            cut_edge(&vv3, &fd, self.my_ds.change_ds(), 1)?;
        }
        if st.v4.is_vertex() && fd.is_on_curve2() {
            let vv4 = st
                .v4
                .vertex()
                .cloned()
                .ok_or_else(|| "ChFi3d_FilDS: V4 is a vertex without a shape".to_string())?;
            cut_edge(&vv4, &fd, self.my_ds.change_ds(), 2)?;
        }

        if j == 1 {
            self.initial_processing(st, &fd, cor_dat, isurf, is_in_ds1, &mut et1)?;
        } else if !is_in_ds1 {
            self.inter_fillet(st, cor_dat, isurf, &fd, et1, tol3d)?;
        }

        self.fillet_face_interference(
            st, seq_fil, j, cor_dat, isurf, &fd, &fi1, 1, ishape1, is_in_ds2,
        )?;
        self.fillet_face_interference(
            st, seq_fil, j, cor_dat, isurf, &fd, &fi2, 2, ishape2, is_in_ds2,
        )?;

        // `ET1 = trafil1;` (`cxx:3160`).
        et1 = st.trafil1;

        if j == len {
            if !is_in_ds2 {
                self.last_curve(st, cor_dat, isurf, &mut et1)?;
            }
        } else if !is_in_ds2 {
            self.middle_curve(st, &fd, &fi1, &fi2, isurf, et1, tol3d, tol2d)?;
        }

        // `Degene = V3.Point().IsEqual(V4.Point(), 0);` (`cxx:3242`).
        st.degene = same_point(st.v3.point(), st.v4.point());

        if st.degene {
            self.degenerated_end(st, cor_dat, seq_fil, j)?;
        } else if !st.closed || j != len {
            self.point_edge_interference(st, &fd, 1)?;
            self.point_edge_interference(st, &fd, 2)?;
        }
        Ok(())
    }
}
