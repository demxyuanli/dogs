//! `PerformTwoCornerbyInter`'s `sameparam` branch -
//! `ChFi3d_Builder_C2.cxx:306-446` (OCCT 8.0.0).
//!
//! The two fillets meet on their common face, so the two tangent lines intersect
//! there. The body intersects the two bounded fillet surfaces into one curve
//! `Gc` with its two p-curves, and hands that single curve to both stripes: the
//! `SurfData` of each gets the same `CommonPoint`s (`Fd2` copies `Fd1`'s
//! outright) and the same point indices, so `Corner2` shares `Corner1`'s
//! points in the DS.
//!
//! OCCT reads `IFaCo1 == 1` / `== 2` to decide which of `Pardeb` / `Parfin` is
//! the p-curve pair of `HS1`; any other value leaves `Gc` null and the body
//! dereferences it, which the port answers as `Err`.

use std::sync::Arc;

use occt_core::bnd::BndBox;
use occt_geom::Curve;

use crate::abs::Orientation;
use crate::chfi3d::builder_0_helpers::{
    ch_fi3d_compute_curves, ch_fi3d_computes_int_pc, ch_fi3d_enlarge_box_edge,
    ch_fi3d_enlarge_box_surf_data, ch_fi3d_index_point_in_ds, ch_fi3d_set_point_tolerance,
    ComputeCurves,
};
use crate::chfi3d::builder_c2::reduce_surface;
use crate::topopebrep::ds::Curve as DsCurve;

use super::super::ChFi3dBuilder;
use super::{index_point, p_curve_mut, p_curve_on_surf, surf_value, TwoCorner};

/// `cxx:306-446`.
pub(super) fn same_param(b: &mut ChFi3dBuilder, tc: &mut TwoCorner) -> Result<bool, String> {
    let is_first1 = tc.is_first1;
    let is_first2 = tc.is_first2;
    let i_fa_co1 = tc.i_fa_co1;
    let i_fa_co2 = tc.i_fa_co2;
    let i_fa_arc1 = tc.i_fa_arc1;
    let i_fa_arc2 = tc.i_fa_arc2;
    let mut u_int_pc1 = tc.u_int_pc1;
    let mut u_int_pc2 = tc.u_int_pc2;

    // `ChFi3d_ComputesIntPC(Fd1->Interference(IFaCo1), Fd2->Interference(IFaCo2),
    // HS1, HS2, UIntPC1, UIntPC2);` (`cxx:310-317`).
    ch_fi3d_computes_int_pc(
        tc.fd1.interference(i_fa_co1),
        tc.fd2.interference(i_fa_co2),
        &tc.hs1,
        &tc.hs2,
        &mut u_int_pc1,
        &mut u_int_pc2,
    )?;
    // `UV = Fd1->Interference(IFaCo1).PCurveOnSurf()->Value(UIntPC1); Pardeb(1)
    // = UV.X(); Pardeb(2) = UV.Y();` (`cxx:314-319`).
    let uv = p_curve_on_surf(&tc.fd1, i_fa_co1)?.value(u_int_pc1);
    tc.pardeb[0] = uv.x();
    tc.pardeb[1] = uv.y();
    // `UV = Fd2->Interference(IFaCo2).PCurveOnSurf()->Value(UIntPC2); Pardeb(3)
    // = UV.X(); Pardeb(4) = UV.Y();` (`cxx:320-323`).
    let uv = p_curve_on_surf(&tc.fd2, i_fa_co2)?.value(u_int_pc2);
    tc.pardeb[2] = uv.x();
    tc.pardeb[3] = uv.y();
    // `gp_Pnt PFaCo = HS1->Value(Pardeb(1), Pardeb(2));` (`cxx:324`).
    let p_fa_co = surf_value(&tc.hs1, tc.pardeb[0], tc.pardeb[1])?;

    // `double UIntArc1 = Fd1->Interference(IFaArc1).Parameter(isfirst1);` and
    // `UIntArc2` (`cxx:327-328`).
    let mut u_int_arc1 = tc.fd1.interference(i_fa_arc1).parameter(is_first1);
    let mut u_int_arc2 = tc.fd2.interference(i_fa_arc2).parameter(is_first2);
    // `ChFi3d_ComputesIntPC(Fd1->Interference(IFaArc1), Fd2->Interference(IFaArc2),
    // HS1, HS2, UIntArc1, UIntArc2);` (`cxx:329-336`).
    ch_fi3d_computes_int_pc(
        tc.fd1.interference(i_fa_arc1),
        tc.fd2.interference(i_fa_arc2),
        &tc.hs1,
        &tc.hs2,
        &mut u_int_arc1,
        &mut u_int_arc2,
    )?;
    // `UV = Fd1->Interference(IFaArc1).PCurveOnSurf()->Value(UIntArc1);
    // Parfin(1) = UV.X(); Parfin(2) = UV.Y();` (`cxx:337-339`).
    let uv = p_curve_on_surf(&tc.fd1, i_fa_arc1)?.value(u_int_arc1);
    tc.parfin[0] = uv.x();
    tc.parfin[1] = uv.y();
    // `UV = Fd2->Interference(IFaArc2).PCurveOnSurf()->Value(UIntArc2);
    // Parfin(3) = UV.X(); Parfin(4) = UV.Y();` (`cxx:339-341`).
    let uv = p_curve_on_surf(&tc.fd2, i_fa_arc2)?.value(u_int_arc2);
    tc.parfin[2] = uv.x();
    tc.parfin[3] = uv.y();

    // `if (Fd1->Surf() == Fd2->Surf()) { Reduce(UIntPC1, UIntPC2, HS1, HS2); }`
    // (`cxx:343-346`).
    if tc.fd1.surf() == tc.fd2.surf() {
        reduce_surface(u_int_pc1, u_int_pc2, &mut tc.hs1, &mut tc.hs2);
    }

    // `double tolreached = tolapp3d; if (IFaCo1 == 1 && !ChFi3d_ComputeCurves(HS1,
    // HS2, Pardeb, Parfin, Gc, PGc1, PGc2, tolapp3d, tol2d, tolreached)) { ... }
    // else if (IFaCo1 == 2 && !ChFi3d_ComputeCurves(..., Parfin, Pardeb, ...))`
    // (`cxx:348-383`). The ten-argument calls keep `wholeCurv = true`.
    let mut computed = ComputeCurves {
        tolreached: b.tolapp3d,
        ..ComputeCurves::default()
    };
    if i_fa_co1 == 1 {
        if !ch_fi3d_compute_curves(
            &tc.hs1,
            &tc.hs2,
            tc.pardeb,
            tc.parfin,
            &mut computed,
            b.tolapp3d,
            b.tol2d,
            true,
        ) {
            return Ok(false);
        }
    } else if i_fa_co1 == 2 {
        if !ch_fi3d_compute_curves(
            &tc.hs1,
            &tc.hs2,
            tc.parfin,
            tc.pardeb,
            &mut computed,
            b.tolapp3d,
            b.tol2d,
            true,
        ) {
            return Ok(false);
        }
    }
    // `double WFirst = Gc->FirstParameter();` (`cxx:385`): a null `Gc` is the
    // OCCT null dereference.
    let gc: Arc<dyn Curve> = computed.c3d.clone().ok_or_else(|| {
        "PerformTwoCornerbyInter: ChFi3d_ComputeCurves left Gc null".to_string()
    })?;
    let w_first = gc.first_parameter();
    let w_last = gc.last_parameter();

    // `double tolpco = std::max(cpco1.Tolerance(), cpco2.Tolerance());
    // double tolparc = std::max(cparc1.Tolerance(), cparc2.Tolerance());`
    // (`cxx:392-395`).
    let tolpco = tc
        .fd1
        .vertex(is_first1, i_fa_co1)
        .tolerance()
        .max(tc.fd2.vertex(is_first2, i_fa_co2).tolerance());
    let mut tolparc = tc
        .fd1
        .vertex(is_first1, i_fa_arc1)
        .tolerance()
        .max(tc.fd2.vertex(is_first2, i_fa_arc2).tolerance());

    // `int ICurv = DStr.AddCurve(TopOpeBRepDS_Curve(Gc, tolreached));`
    // (`cxx:396`).
    let icurv = b.my_ds.change_ds().add_curve(&mut DsCurve::with_curve(
        Some(Arc::clone(&gc)),
        computed.tolreached,
        false,
    ));

    // `Corner1->SetParameters(isfirst1, WFirst, WLast); Corner1->SetCurve(ICurv,
    // isfirst1); Corner1->ChangePCurve(isfirst1) = PGc1;` (`cxx:398-399`).
    {
        let mut stripe = tc
            .corner1
            .write()
            .map_err(|_| "PerformTwoCornerbyInter: poisoned stripe".to_string())?;
        stripe.set_parameters(is_first1, w_first, w_last);
        stripe.set_curve(icurv, is_first1);
        *p_curve_mut(&mut stripe, is_first1) = computed.pc1.clone();
    }
    // `cpco1.Reset(); cpco1.SetPoint(PFaCo); cpco1.SetTolerance(std::max(
    // tolreached, tolpco)); Fd1->ChangeInterference(IFaCo1).SetParameter(UIntPC1,
    // isfirst1);` (`cxx:400-403`).
    {
        let cp = tc.fd1.vertex_mut(is_first1, i_fa_co1);
        cp.reset();
        cp.set_point(&p_fa_co);
        cp.set_tolerance(computed.tolreached.max(tolpco));
    }
    tc.fd1
        .interference_mut(i_fa_co1)
        .set_parameter(u_int_pc1, is_first1);
    // `tolparc = std::max(tolparc, tolreached); cparc1.SetTolerance(std::max(
    // tolparc, tolreached));` (`cxx:404-405`).
    tolparc = tolparc.max(computed.tolreached);
    tc.fd1
        .vertex_mut(is_first1, i_fa_arc1)
        .set_tolerance(tolparc.max(computed.tolreached));
    // `Ipoin1 = ChFi3d_IndexPointInDS(Fd1->Vertex(isfirst1, 1), DStr);
    // Corner1->SetIndexPoint(Ipoin1, isfirst1, 1);` (`cxx:406-407`).
    let ipoin1 = ch_fi3d_index_point_in_ds(
        tc.fd1.vertex(is_first1, 1),
        b.my_ds.change_ds(),
    );
    {
        let mut stripe = tc
            .corner1
            .write()
            .map_err(|_| "PerformTwoCornerbyInter: poisoned stripe".to_string())?;
        stripe.set_index_point(ipoin1, is_first1, 1);
    }
    // `Ipoin2 = ChFi3d_IndexPointInDS(Fd1->Vertex(isfirst1, 2), DStr);
    // Corner1->SetIndexPoint(Ipoin2, isfirst1, 2);` (`cxx:408-409`).
    let ipoin2 = ch_fi3d_index_point_in_ds(
        tc.fd1.vertex(is_first1, 2),
        b.my_ds.change_ds(),
    );
    {
        let mut stripe = tc
            .corner1
            .write()
            .map_err(|_| "PerformTwoCornerbyInter: poisoned stripe".to_string())?;
        stripe.set_index_point(ipoin2, is_first1, 2);
    }

    // `Corner2->SetParameters(isfirst2, WFirst, WLast); Corner2->SetCurve(ICurv,
    // isfirst2); Corner2->ChangePCurve(isfirst2) = PGc2;` (`cxx:411-413`).
    {
        let mut stripe = tc
            .corner2
            .write()
            .map_err(|_| "PerformTwoCornerbyInter: poisoned stripe".to_string())?;
        stripe.set_parameters(is_first2, w_first, w_last);
        stripe.set_curve(icurv, is_first2);
        *p_curve_mut(&mut stripe, is_first2) = computed.pc2.clone();
    }
    // `Fd2->ChangeInterference(IFaCo2).SetParameter(UIntPC2, isfirst2);`
    // (`cxx:414`).
    tc.fd2
        .interference_mut(i_fa_co2)
        .set_parameter(u_int_pc2, is_first2);
    // `Fd2->ChangeVertex(isfirst2, IFaCo2) = Fd1->Vertex(isfirst1, IFaCo1);
    // Fd2->ChangeVertex(isfirst2, IFaArc2) = Fd1->Vertex(isfirst1, IFaArc1);`
    // (`cxx:415-416`).
    let src_co = tc.fd1.vertex(is_first1, i_fa_co1).clone();
    *tc.fd2.vertex_mut(is_first2, i_fa_co2) = src_co;
    let src_arc = tc.fd1.vertex(is_first1, i_fa_arc1).clone();
    *tc.fd2.vertex_mut(is_first2, i_fa_arc2) = src_arc;
    // `if (IFaCo1 != IFaCo2) { Corner2->SetOrientation(TopAbs_REVERSED,
    // isfirst2); }` (`cxx:417-420`).
    if i_fa_co1 != i_fa_co2 {
        let mut stripe = tc
            .corner2
            .write()
            .map_err(|_| "PerformTwoCornerbyInter: poisoned stripe".to_string())?;
        stripe.set_orientation(Orientation::Reversed, is_first2);
    }
    // `Corner2->SetIndexPoint(Corner1->IndexPoint(isfirst1, IFaCo1), isfirst2,
    // IFaCo2); Corner2->SetIndexPoint(Corner1->IndexPoint(isfirst1, IFaArc1),
    // isfirst2, IFaArc2);` (`cxx:421-422`).
    let icc = index_point(&tc.corner1, is_first1, i_fa_co1)?;
    {
        let mut stripe = tc
            .corner2
            .write()
            .map_err(|_| "PerformTwoCornerbyInter: poisoned stripe".to_string())?;
        stripe.set_index_point(icc, is_first2, i_fa_co2);
    }
    let iarc = index_point(&tc.corner1, is_first1, i_fa_arc1)?;
    {
        let mut stripe = tc
            .corner2
            .write()
            .map_err(|_| "PerformTwoCornerbyInter: poisoned stripe".to_string())?;
        stripe.set_index_point(iarc, is_first2, i_fa_arc2);
    }

    // `Bnd_Box bco, barc; if (IFaCo1 == 1) ChFi3d_EnlargeBox(DStr, Corner1, Fd1,
    // bco, barc, isfirst1); else ChFi3d_EnlargeBox(DStr, Corner1, Fd1, barc, bco,
    // isfirst1);` and the `Corner2` / `Fd2` twin (`cxx:424-439`).
    let mut bco = BndBox::new();
    let mut barc = BndBox::new();
    {
        let stripe = tc
            .corner1
            .read()
            .map_err(|_| "PerformTwoCornerbyInter: poisoned stripe".to_string())?;
        let ds = b.my_ds.change_ds();
        if i_fa_co1 == 1 {
            ch_fi3d_enlarge_box_surf_data(
                ds,
                Some(&stripe),
                &tc.fd1,
                &mut bco,
                &mut barc,
                is_first1,
            )?;
        } else {
            ch_fi3d_enlarge_box_surf_data(
                ds,
                Some(&stripe),
                &tc.fd1,
                &mut barc,
                &mut bco,
                is_first1,
            )?;
        }
    }
    {
        let stripe = tc
            .corner2
            .read()
            .map_err(|_| "PerformTwoCornerbyInter: poisoned stripe".to_string())?;
        let ds = b.my_ds.change_ds();
        if i_fa_co2 == 1 {
            ch_fi3d_enlarge_box_surf_data(
                ds,
                Some(&stripe),
                &tc.fd2,
                &mut bco,
                &mut barc,
                is_first2,
            )?;
        } else {
            ch_fi3d_enlarge_box_surf_data(
                ds,
                Some(&stripe),
                &tc.fd2,
                &mut barc,
                &mut bco,
                is_first2,
            )?;
        }
    }

    // `const ChFiDS_CommonPoint& cparc = Fd1->Vertex(isfirst1, IFaArc1);
    // ChFi3d_EnlargeBox(cparc.Arc(), myEFMap(cparc.Arc()),
    // cparc.ParameterOnArc(), barc);` (`cxx:441-443`).
    let cparc = tc.fd1.vertex(is_first1, i_fa_arc1).clone();
    let arc = cparc
        .arc()
        .cloned()
        .ok_or_else(|| "PerformTwoCornerbyInter: cparc has no arc".to_string())?;
    let par = cparc.parameter_on_arc().ok_or_else(|| {
        "PerformTwoCornerbyInter: cparc has no parameter on arc".to_string()
    })?;
    let faces: Vec<crate::shape::TopoShape> = b
        .my_ef_map
        .get_from_key(&arc)
        .to_vec();
    ch_fi3d_enlarge_box_edge(&arc, &faces, par, &mut barc);

    // `ChFi3d_SetPointTolerance(DStr, barc, Corner1->IndexPoint(isfirst1,
    // IFaArc1)); ChFi3d_SetPointTolerance(DStr, bco, Corner1->IndexPoint(isfirst1,
    // IFaCo1));` (`cxx:444-445`).
    let ip_arc = index_point(&tc.corner1, is_first1, i_fa_arc1)?;
    ch_fi3d_set_point_tolerance(b.my_ds.change_ds(), &barc, ip_arc);
    let ip_co = index_point(&tc.corner1, is_first1, i_fa_co1)?;
    ch_fi3d_set_point_tolerance(b.my_ds.change_ds(), &bco, ip_co);

    // `done = true;` (`cxx:802`): the two branches share the single `done = true`
    // of the tail, reported through `Ok(true)`.
    Ok(true)
}
