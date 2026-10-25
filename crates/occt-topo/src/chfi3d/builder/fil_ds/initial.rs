//! The initial processing (`j == 1`) and the inter-fillet interference
//! (`j > 1`) of `ChFi3d_FilDS`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx`:
//!
//! - `if (j == 1)` (`cxx:2779-2887`), the "End of the Initial Processing":
//!   the two first common points either form the "queue de billard" singular
//!   end or file their arc interferences, then the first curve of the stripe
//!   enters the DS with its two end points and its surface-curve interference,
//! - `else { // ---- Interference between Fillets ------ }` (`cxx:2888-2952`):
//!   the curve joining the previous `SurfData` to this one, either as a
//!   p-curve only (degenerated) or as a `ChFi3d_ComputePCurv`-projected
//!   p-curve on the fillet surface.
//!
//! Translation notes:
//! - `BoutdeVtx = V1.Vertex()` writes the outer `BoutdeVtx` of `cxx:2586`
//!   (`FilDsState::bout_de_vtx`) and is also used by the `j > 1` block, so the
//!   two functions share the state,
//! - `ChFi3d_ComputePCurv(UV1, UV2, PCurv, Pardeb, Parfin)` (`cxx:2926`) is the
//!   2D overload, [`ch_fi3d_compute_pcurv_uv`] with the default
//!   `reverse = false` of `ChFi3d_Builder_0.hxx:281-286`,
//! - `ChFi3d_ComputePCurv(TCurv.ChangeCurve(), ...)` (`cxx:2938`) reads the DS
//!   curve's 3D curve; a null handle would crash the OCCT adaptor, so the port
//!   reports it as an `Err`,
//! - the `reglist.Append(regfilfil)` of `cxx:2936` happens before the
//!   computation, exactly like OCCT, so a failing compute leaves the
//!   regulation list as OCCT would before throwing.

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::chfi_ds::{StripeHandle, SurfData};
use crate::chfi3d::builder_0_helpers::compute_pcurv::{
    ch_fi3d_compute_pcurv, ch_fi3d_compute_pcurv_uv,
};
use crate::chfi3d::builder_0_helpers::fil_interf::{
    ch_fi3d_contains, ch_fi3d_fil_curve_in_ds, ch_fi3d_fil_point_in_ds, ch_fi3d_fil_vertex_in_ds,
};

use super::super::spine_of;
use super::{p_curve_value, same_point, with_stripe, FilDsState};
use crate::chfi3d::builder::ChFi3dBuilder;

impl ChFi3dBuilder {
    /// `if (j == 1) { ... }` (`cxx:2779-2887`).
    #[allow(clippy::too_many_arguments)]
    pub(super) fn initial_processing(
        &mut self,
        st: &mut FilDsState,
        fd: &SurfData,
        cor_dat: &StripeHandle,
        isurf: i32,
        is_in_ds1: bool,
        et1: &mut Orientation,
    ) -> Result<(), String> {
        let v1 = fd.vertex_first_on_s1().clone();
        let v2 = fd.vertex_first_on_s2().clone();

        // `isVertex1 = V1.IsVertex(); isVertex2 = V2.IsVertex();
        // Singulier_en_Bout = (V1.Point().IsEqual(V2.Point(), 0));`
        // (`cxx:2781-2783`).
        st.is_vertex1 = v1.is_vertex();
        st.is_vertex2 = v2.is_vertex();
        st.singulier_en_bout = same_point(v1.point(), v2.point());

        if st.singulier_en_bout {
            // Queue de Billard (`cxx:2785-2816`).
            if v1.is_vertex() && v2.is_vertex() {
                st.is_vertex1 = true; // caution... (`cxx:2793`)
                st.is_vertex2 = true;
                // The edge is removed from spine starting on this vertex
                // (`cxx:2794-2798`).
                let spine = spine_of(cor_dat)
                    .ok_or_else(|| "ChFi3d_FilDS: the stripe has no spine".to_string())?;
                let arcspine = {
                    let guard = spine
                        .read()
                        .map_err(|_| "ChFi3d_FilDS: poisoned spine".to_string())?;
                    guard.edges(1).cloned().ok_or_else(|| {
                        "ChFi3d_FilDS: spine->Edges(1) is out of range".to_string()
                    })?
                };
                let bout = v1
                    .vertex()
                    .cloned()
                    .ok_or_else(|| "ChFi3d_FilDS: V1 is a vertex without a shape".to_string())?;
                st.bout_de_vtx = Some(bout.clone());
                let i_arcspine = self.my_ds.change_ds().add_shape(&arcspine);
                let ivtx = with_stripe(cor_dat, |s| s.index_first_point_on_s1())?;

                // `OVtx = ex.Current().Orientation(); ... OVtx =
                // TopAbs::Reverse(OVtx);` (`cxx:2802-2810`).
                let ovtx = super::vertex_orientation_in_edge(&arcspine, &bout).reversed();
                // `double parVtx = BRep_Tool::Parameter(BoutdeVtx,
                // Arcspine);` (`cxx:2811`).
                let par_vtx = BRepTool::parameter_of_vertex(&arcspine, &bout);
                // `ChFi3d_FilVertexInDS(OVtx, IArcspine, IVtx, parVtx)` and its
                // append (`cxx:2812-2814`).
                let interfv = ch_fi3d_fil_vertex_in_ds(ovtx, i_arcspine, ivtx, par_vtx);
                self.my_ds
                    .change_ds()
                    .change_shape_interferences_at(i_arcspine)
                    .push(interfv);
            }
        } else {
            // `if (V1.IsOnArc()) { Iarc1 = DStr.AddShape(V1.Arc()); if
            // (!ChFi3d_Contains(...)) {...} }` (`cxx:2819-2831`).
            if v1.is_on_arc() {
                let arc = v1
                    .arc()
                    .cloned()
                    .ok_or_else(|| "ChFi3d_FilDS: V1 is on an arc without a shape".to_string())?;
                let iarc1 = self.my_ds.change_ds().add_shape(&arc);
                st.iarc1 = iarc1;
                if !ch_fi3d_contains(
                    self.ds().shape_interferences_at(iarc1, false),
                    iarc1,
                    st.ipo_in1,
                    false,
                    false,
                ) {
                    let interfp1 = ch_fi3d_fil_point_in_ds(
                        v1.transition_on_arc().ok_or_else(|| {
                            "ChFi3d_FilDS: V1 is on an arc without a transition".to_string()
                        })?,
                        iarc1,
                        st.ipo_in1,
                        v1.parameter_on_arc().ok_or_else(|| {
                            "ChFi3d_FilDS: V1 is on an arc without a parameter".to_string()
                        })?,
                        st.is_vertex1,
                    );
                    self.my_ds
                        .change_ds()
                        .change_shape_interferences(&arc)
                        .push(interfp1);
                }
            }

            // `if (V2.IsOnArc()) { ... }` (`cxx:2833-2845`).
            if v2.is_on_arc() {
                let arc = v2
                    .arc()
                    .cloned()
                    .ok_or_else(|| "ChFi3d_FilDS: V2 is on an arc without a shape".to_string())?;
                let iarc2 = self.my_ds.change_ds().add_shape(&arc);
                st.iarc2 = iarc2;
                if !ch_fi3d_contains(
                    self.ds().shape_interferences_at(iarc2, false),
                    iarc2,
                    st.ipo_in2,
                    false,
                    false,
                ) {
                    let interfp2 = ch_fi3d_fil_point_in_ds(
                        v2.transition_on_arc().ok_or_else(|| {
                            "ChFi3d_FilDS: V2 is on an arc without a transition".to_string()
                        })?,
                        iarc2,
                        st.ipo_in2,
                        v2.parameter_on_arc().ok_or_else(|| {
                            "ChFi3d_FilDS: V2 is on an arc without a parameter".to_string()
                        })?,
                        st.is_vertex2,
                    );
                    self.my_ds
                        .change_ds()
                        .change_shape_interferences(&arc)
                        .push(interfp2);
                }
            }
        }

        // `if (!isInDS1) { ... }` (`cxx:2848-2886`).
        if !is_in_ds1 {
            // `ET1 = TopAbs::Compose(ET1,
            // CorDat->FirstPCurveOrientation());` (`cxx:2850`).
            let ori = with_stripe(cor_dat, |s| s.first_p_curve_orientation())?;
            *et1 = Orientation::compose(*et1, ori);
            // `Icurv = CorDat->FirstCurve();` (`cxx:2851`).
            st.icurv = with_stripe(cor_dat, |s| s.first_curve())?;
            // `if (Closed && !Singulier_en_Bout) { regcout.SetCurve(Icurv);
            // regcout.SetS1(Isurf, false); }` (`cxx:2852-2856`).
            if st.closed && !st.singulier_en_bout {
                st.regcout.set_curve(st.icurv);
                st.regcout.set_s1(isurf, false);
            }
            // `PCurv = CorDat->FirstPCurve(); CorDat->FirstParameters(Pardeb,
            // Parfin);` (`cxx:2857-2858`).
            st.p_curv = with_stripe(cor_dat, |s| s.first_p_curve().cloned())?;
            let (pardeb, parfin) = with_stripe(cor_dat, |s| s.first_parameters())?;
            st.pardeb = pardeb;
            st.parfin = parfin;

            // `NCollection_List<...>& Li = DStr.ChangeCurveInterferences(Icurv);
            // if (Li.IsEmpty()) { ... }` (`cxx:2860-2876`).
            if self.ds().curve_interferences(st.icurv).is_empty() {
                let (o1, p1, o2, p2) = if ori == Orientation::Reversed {
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
                let interfp1 =
                    ch_fi3d_fil_point_in_ds(o1, st.icurv, st.ipo_in1, p1, st.is_vertex1);
                let interfp2 =
                    ch_fi3d_fil_point_in_ds(o2, st.icurv, st.ipo_in2, p2, st.is_vertex2);
                let li = self.my_ds.change_ds().change_curve_interferences(st.icurv);
                li.push(interfp1);
                li.push(interfp2);
            }

            // `Interfc1 = ChFi3d_FilCurveInDS(Icurv, Isurf, PCurv, ET1);
            // DStr.ChangeSurfaceInterferences(Isurf).Append(Interfc1);`
            // (`cxx:2877-2878`).
            let interfc1 =
                ch_fi3d_fil_curve_in_ds(st.icurv, isurf, st.p_curv.clone(), *et1);
            self.my_ds
                .change_ds()
                .change_surface_interferences(isurf)
                .push(interfc1.clone());

            // `if (Ipoin1 == Ipoin2) { ... TCurv.ChangeCurve().Nullify();
            // TCurv.SetSCI(Interfc1, bidinterf); }` (`cxx:2879-2885`).
            if st.ipo_in1 == st.ipo_in2 {
                let tcurv = self.my_ds.change_ds().change_curve(st.icurv);
                *tcurv.change_curve() = None;
                tcurv.set_sci(Some(interfc1), None);
            }
        }
        Ok(())
    }

    /// `else { // ---- Interference between Fillets ------ }`
    /// (`cxx:2888-2952`).
    pub(super) fn inter_fillet(
        &mut self,
        st: &mut FilDsState,
        cor_dat: &StripeHandle,
        isurf: i32,
        fd: &SurfData,
        et1: Orientation,
        tol3d: f64,
    ) -> Result<(), String> {
        let fi1 = fd.interference_on_s1().clone();
        let fi2 = fd.interference_on_s2().clone();

        // `if (Degene && isVertex1) { NumEdge++; Arcspine =
        // spine->Edges(NumEdge); ... }` (`cxx:2895-2916`).
        if st.degene && st.is_vertex1 {
            st.num_edge += 1;
            let spine = spine_of(cor_dat)
                .ok_or_else(|| "ChFi3d_FilDS: the stripe has no spine".to_string())?;
            let arcspine = {
                let guard = spine
                    .read()
                    .map_err(|_| "ChFi3d_FilDS: poisoned spine".to_string())?;
                guard.edges(st.num_edge).cloned().ok_or_else(|| {
                    format!("ChFi3d_FilDS: spine->Edges({}) is out of range", st.num_edge)
                })?
            };
            let i_arcspine = self.my_ds.change_ds().add_shape(&arcspine);
            let bout = st.bout_de_vtx.clone().ok_or_else(|| {
                "ChFi3d_FilDS: BoutdeVtx is unset for the degenerated start".to_string()
            })?;
            let ivtx = self.my_ds.change_ds().add_shape(&bout);
            let ovtx = super::vertex_orientation_in_edge(&arcspine, &bout).reversed();
            let par_vtx = BRepTool::parameter_of_vertex(&arcspine, &bout);
            let interfv = ch_fi3d_fil_vertex_in_ds(ovtx, i_arcspine, ivtx, par_vtx);
            self.my_ds
                .change_ds()
                .change_shape_interferences_at(i_arcspine)
                .push(interfv);
        }

        // `gp_Pnt2d UV1 =
        // Fd->InterferenceOnS1().PCurveOnSurf()->Value(Fd->InterferenceOnS1().FirstParameter());`
        // and `UV2` (`cxx:2918-2921`).
        let uv1 = p_curve_value(
            fi1.p_curve_on_surf(),
            fi1.first_parameter(),
            "the S1 interference",
        )?;
        let uv2 = p_curve_value(
            fi2.p_curve_on_surf(),
            fi2.first_parameter(),
            "the S2 interference",
        )?;

        if st.degene {
            // `// pcurve is associated via SCI to TopOpeBRepDSCurve.
            // ChFi3d_ComputePCurv(UV1, UV2, PCurv, Pardeb, Parfin);`
            // (`cxx:2925-2931`).
            st.p_curv = Some(ch_fi3d_compute_pcurv_uv(
                &uv1,
                &uv2,
                st.pardeb,
                st.parfin,
                false,
            )?);
            let interfc1 =
                ch_fi3d_fil_curve_in_ds(st.icurv, isurf, st.p_curv.clone(), et1);
            self.my_ds
                .change_ds()
                .change_surface_interferences(isurf)
                .push(interfc1.clone());
            let tcurv = self.my_ds.change_ds().change_curve(st.icurv);
            *tcurv.change_curve() = None;
            tcurv.set_sci(Some(interfc1), None);
        } else {
            // `regfilfil.SetS2(Isurf, false); reglist.Append(regfilfil);`
            // (`cxx:2935-2936`).
            st.regfilfil.set_s2(isurf, false);
            self.append_regul(st.regfilfil);

            // `ChFi3d_ComputePCurv(TCurv.ChangeCurve(), UV1, UV2, PCurv,
            // DStr.Surface(Fd->Surf()).Surface(), Pardeb, Parfin, tol3d,
            // tolreached);` (`cxx:2938-2946`).
            let c3d = self
                .ds()
                .curve(st.icurv)
                .curve()
                .cloned()
                .ok_or_else(|| {
                    "ChFi3d_FilDS: the fillet/fillet curve has no 3D curve".to_string()
                })?;
            let surface = self
                .ds()
                .surface(fd.surf())
                .surface()
                .cloned()
                .ok_or_else(|| "ChFi3d_FilDS: the fillet surface is not in the DS".to_string())?;
            let tolreached = ch_fi3d_compute_pcurv(
                c3d.as_ref(),
                &uv1,
                &uv2,
                &mut st.p_curv,
                surface.as_ref(),
                st.pardeb,
                st.parfin,
                tol3d,
                false,
            )?;

            // `TCurv.Tolerance(std::max(TCurv.Tolerance(), tolreached));`
            // (`cxx:2947`).
            let tcurv = self.my_ds.change_ds().change_curve(st.icurv);
            let tol = tcurv.tolerance().max(tolreached);
            tcurv.set_tolerance(tol);

            // `Interfc1 = ChFi3d_FilCurveInDS(Icurv, Isurf, PCurv, ET1);
            // DStr.ChangeSurfaceInterferences(Isurf).Append(Interfc1);`
            // (`cxx:2948-2949`).
            let interfc1 =
                ch_fi3d_fil_curve_in_ds(st.icurv, isurf, st.p_curv.clone(), et1);
            self.my_ds
                .change_ds()
                .change_surface_interferences(isurf)
                .push(interfc1);
        }
        Ok(())
    }
}
