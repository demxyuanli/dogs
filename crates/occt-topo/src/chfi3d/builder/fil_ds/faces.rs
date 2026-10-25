//! The "Interference Fillets / Faces" blocks of `ChFi3d_FilDS`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:2954-3158`
//! (OCCT 8.0.0). OCCT writes the block twice, once for S1
//! (`cxx:2954-3054`) and once for S2 (`cxx:3056-3158`); the only differences
//! are the `...OnS1` / `...OnS2` accessors plus the extra "Pinch" branch of the
//! S2 block (`cxx:3116-3120`), so the port drives both from one function with
//! a `side` argument (`1` or `2`).
//!
//! Translation notes:
//! - `Fd->IndexOfS1()` is re-read at the top of each block (`cxx:2961`,
//!   `:3061`) although `cxx:2661` already read it, and the local is negated
//!   below (`cxx:2975`, `:3076`); the port keeps that local so the negation
//!   does not leak into the next `j`,
//! - `Ishape1 < 0` stores the interference on the *surface* list of `-Ishape1`
//!   and `Ishape1 > 0` on the *shape* list (`cxx:2977-2987`), which the port
//!   spells `change_surface_interferences` / `change_shape_interferences_at`,
//! - `NCollection_List<...>& Li = DStr.ChangeCurveInterferences(IcFil)` is a
//!   live list (`cxx:3036-3038`, `:3141-3143`): the second `ChFi3d_Contains`
//!   of `cxx:3047` / `:3150` sees the interference appended by the first test,
//!   so the port re-reads `curve_interferences` between the two appends,
//! - `gp_Pnt::IsEqual(Other, 1.0e-7)` and `(Other, 0)` are
//!   [`GpPnt::is_equal_with_tolerance`] with those tolerances,
//! - `DStr.Curve(I).Curve().IsNull()` on an unbound `I` answers the empty
//!   curve in the port ([`DataStructure::curve`]), the same answer
//!   `myEmptyCurve` gives.

use std::sync::Arc;

use crate::abs::Orientation;
use crate::chfi_ds::{CommonPoint, FaceInterference, Regul, StripeHandle, SurfData};
use crate::chfi3d::builder_0_helpers::fil_interf::{
    ch_fi3d_contains, ch_fi3d_fil_curve_in_ds, ch_fi3d_fil_point_in_ds, find_index_point,
};
use crate::chfi3d::builder_0_helpers::point_in_ds::ch_fi3d_index_point_in_ds;

use super::{with_stripe, FilDsState};
use crate::chfi3d::builder::ChFi3dBuilder;

/// `Fd->InterferenceOnS1()` / `...OnS2()` selected by `side`.
fn interference_on(sd: &SurfData, side: i32) -> &FaceInterference {
    if side == 1 {
        sd.interference_on_s1()
    } else {
        sd.interference_on_s2()
    }
}

/// `Fd->VertexFirstOnS1()` / `...OnS2()` selected by `side`.
fn vertex_first_on(sd: &SurfData, side: i32) -> &CommonPoint {
    if side == 1 {
        sd.vertex_first_on_s1()
    } else {
        sd.vertex_first_on_s2()
    }
}

/// `Fd->VertexLastOnS1()` / `...OnS2()` selected by `side`.
fn vertex_last_on(sd: &SurfData, side: i32) -> &CommonPoint {
    if side == 1 {
        sd.vertex_last_on_s1()
    } else {
        sd.vertex_last_on_s2()
    }
}

impl ChFi3dBuilder {
    /// One "Interference Fillets / Faces" block (`cxx:2954-3054` for
    /// `side == 1`, `cxx:3056-3158` for `side == 2`).
    #[allow(clippy::too_many_arguments)]
    pub(super) fn fillet_face_interference(
        &mut self,
        st: &mut FilDsState,
        seq_fil: &[Arc<SurfData>],
        j: usize,
        cor_dat: &StripeHandle,
        isurf: i32,
        fd: &SurfData,
        fi: &FaceInterference,
        side: i32,
        ishape: i32,
        is_in_ds2: bool,
    ) -> Result<(), String> {
        let len = seq_fil.len();
        // `IcFil1 = Fi1.LineIndex();` (`cxx:2955`) / `IcFil2 = ...`
        // (`cxx:3057`).
        let ic_fil = fi.line_index();
        if ic_fil == 0 {
            return Ok(());
        }

        // `Interfc3 = ChFi3d_FilCurveInDS(IcFil1, Isurf, Fi1.PCurveOnSurf(),
        // trafil1); DStr.ChangeSurfaceInterferences(Isurf).Append(Interfc3);`
        // (`cxx:2959-2960`).
        let trafil = if side == 1 { st.trafil1 } else { st.trafil2 };
        let interfc = ch_fi3d_fil_curve_in_ds(ic_fil, isurf, fi.p_curve_on_surf().cloned(), trafil);
        self.my_ds
            .change_ds()
            .change_surface_interferences(isurf)
            .push(interfc.clone());

        // `Ishape1 = Fd->IndexOfS1();` (`cxx:2961`).
        let mut ishape = ishape;

        // `TopOpeBRepDS_Curve& cc = DStr.ChangeCurve(IcFil1); if
        // (cc.Curve().IsNull()) { cc.SetSCI(Interfc3, bidinterf); }`
        // (`cxx:2963-2969`).
        if self.ds().curve(ic_fil).curve().is_none() {
            let cc = self.my_ds.change_ds().change_curve(ic_fil);
            cc.set_sci(Some(interfc), None);
        } else {
            // `ChFiDS_Regul regon1; regon1.SetCurve(IcFil1);
            // regon1.SetS1(Isurf, false);` (`cxx:2971-2974`).
            let mut regon = Regul::new();
            regon.set_curve(ic_fil);
            regon.set_s1(isurf, false);
            if ishape < 0 {
                // `Ishape1 = -Ishape1; regon1.SetS2(Ishape1, false); Interfc1
                // = ChFi3d_FilCurveInDS(IcFil1, Ishape1, Fi1.PCurveOnFace(),
                // Fi1.Transition()); DStr.ChangeSurfaceInterferences(Ishape1)
                // .Append(Interfc1);` (`cxx:2975-2981`).
                ishape = -ishape;
                regon.set_s2(ishape, false);
                let on_face =
                    ch_fi3d_fil_curve_in_ds(ic_fil, ishape, fi.p_curve_on_face().cloned(), fi.transition());
                self.my_ds
                    .change_ds()
                    .change_surface_interferences(ishape)
                    .push(on_face);
            } else if ishape > 0 {
                // `regon1.SetS2(Ishape1, true); Interfc1 =
                // ChFi3d_FilCurveInDS(...); DStr.ChangeShapeInterferences(
                // Ishape1).Append(Interfc1);` (`cxx:2983-2988`).
                regon.set_s2(ishape, true);
                let on_face =
                    ch_fi3d_fil_curve_in_ds(ic_fil, ishape, fi.p_curve_on_face().cloned(), fi.transition());
                self.my_ds
                    .change_ds()
                    .change_shape_interferences_at(ishape)
                    .push(on_face);
            }
            // `reglist.Append(regon1);` (`cxx:2989`).
            self.append_regul(regon);
        }

        // `int ipoin; bool isVertex = Fd->VertexLastOnS1().IsVertex();`
        // (`cxx:2991-2993`).
        let mut is_vertex = vertex_last_on(fd, side).is_vertex();
        let ipoin;

        if j == len {
            // `ipoin = CorDat->IndexLastPointOnS1();` (`cxx:2995-2998`).
            ipoin = with_stripe(cor_dat, |s| {
                if side == 1 {
                    s.index_last_point_on_s1()
                } else {
                    s.index_last_point_on_s2()
                }
            })?;
        } else if j == len - 1
            && self
                .ds()
                .curve(interference_on(&seq_fil[len - 1], side).line_index())
                .curve()
                .is_none()
        {
            // `j == (SeqFil.Length() - 1) && DStr.Curve(SeqFil.Last()->
            // InterferenceOnS1().LineIndex()).Curve().IsNull()`
            // (`cxx:2998-3002`).
            if st.closed {
                ipoin = with_stripe(cor_dat, |s| {
                    if side == 1 {
                        s.index_first_point_on_s1()
                    } else {
                        s.index_first_point_on_s2()
                    }
                })?;
                is_vertex = vertex_first_on(&seq_fil[0], side).is_vertex();
            } else {
                ipoin = with_stripe(cor_dat, |s| {
                    if side == 1 {
                        s.index_last_point_on_s1()
                    } else {
                        s.index_last_point_on_s2()
                    }
                })?;
                is_vertex = vertex_last_on(&seq_fil[len - 1], side).is_vertex();
            }
        } else if self.ds().curve(ic_fil).curve().is_none() {
            // Rotation (`cxx:3011-3015`): the fillet curve is a p-curve only.
            ipoin = if side == 1 { st.ipo_in1 } else { st.ipo_in2 };
            is_vertex = if side == 1 { st.is_vertex1 } else { st.is_vertex2 };
        } else if side == 2
            && fd
                .vertex_last_on_s2()
                .point()
                .is_equal_with_tolerance(fd.vertex_last_on_s1().point(), 0.0)
        {
            // Pinch (`cxx:3116-3120`), S2 only.
            ipoin = st.ipo_in1;
            is_vertex = st.is_vertex1;
        } else if (j == 1 || j == len - 1)
            && (vertex_last_on(fd, side)
                .point()
                .is_equal_with_tolerance(vertex_first_on(&seq_fil[0], side).point(), 1.0e-7)
                || vertex_last_on(fd, side)
                    .point()
                    .is_equal_with_tolerance(vertex_last_on(&seq_fil[len - 1], side).point(), 1.0e-7))
        {
            // Case of SurfData cut in "Triangular" way (`cxx:3016-3028`).
            ipoin = with_stripe(cor_dat, |s| {
                if side == 1 {
                    s.index_last_point_on_s1()
                } else {
                    s.index_last_point_on_s2()
                }
            })?;
        } else {
            // `else if (isInDS2 && findIndexPoint(DStr, Fd, 1, ipoin)) {} else
            // { ipoin = ChFi3d_IndexPointInDS(Fd->VertexLastOnS1(), DStr); }`
            // (`cxx:3028-3034`).
            let found = if is_in_ds2 {
                find_index_point(self.ds(), fd, side)?
            } else {
                None
            };
            ipoin = match found {
                Some(i) => i,
                None => {
                    let cp = vertex_last_on(fd, side).clone();
                    ch_fi3d_index_point_in_ds(&cp, self.my_ds.change_ds())
                }
            };
        }

        let ipo_in = if side == 1 { st.ipo_in1 } else { st.ipo_in2 };
        let is_vertex_in = if side == 1 { st.is_vertex1 } else { st.is_vertex2 };

        // `if (!ChFi3d_Contains(Li, IcFil1, Ipoin1)) { Interfp1 =
        // ChFi3d_FilPointInDS(TopAbs_FORWARD, IcFil1, Ipoin1,
        // Fi1.FirstParameter(), isVertex1); ...Append(Interfp1); }`
        // (`cxx:3039-3045`).
        if !ch_fi3d_contains(
            self.ds().curve_interferences(ic_fil),
            ic_fil,
            ipo_in,
            false,
            false,
        ) {
            let interfp = ch_fi3d_fil_point_in_ds(
                Orientation::Forward,
                ic_fil,
                ipo_in,
                fi.first_parameter(),
                is_vertex_in,
            );
            self.my_ds
                .change_ds()
                .change_curve_interferences(ic_fil)
                .push(interfp);
        }

        // `if (ipoin == Ipoin1 || !ChFi3d_Contains(Li, IcFil1, ipoin)) {
        // Interfp3 = ChFi3d_FilPointInDS(TopAbs_REVERSED, IcFil1, ipoin,
        // Fi1.LastParameter(), isVertex); ...Append(Interfp3); }`
        // (`cxx:3047-3052`).
        if ipoin == ipo_in
            || !ch_fi3d_contains(
                self.ds().curve_interferences(ic_fil),
                ic_fil,
                ipoin,
                false,
                false,
            )
        {
            let interfp = ch_fi3d_fil_point_in_ds(
                Orientation::Reversed,
                ic_fil,
                ipoin,
                fi.last_parameter(),
                is_vertex,
            );
            self.my_ds
                .change_ds()
                .change_curve_interferences(ic_fil)
                .push(interfp);
        }

        // `Ipoin1 = ipoin; isVertex1 = isVertex;` (`cxx:3053-3054`).
        if side == 1 {
            st.ipo_in1 = ipoin;
            st.is_vertex1 = is_vertex;
        } else {
            st.ipo_in2 = ipoin;
            st.is_vertex2 = is_vertex;
        }
        Ok(())
    }
}
