//! `PerformTwoCornerbyInter`'s `else` branch -
//! `ChFi3d_Builder_C2.cxx:447-804` (OCCT 8.0.0).
//!
//! The two tangent lines do not meet on the common face, so the body splits the
//! two stripes into a "big" one and a "small" one and walks the intersection in
//! three steps:
//!
//! 1. `ChFi3d_ComputesIntPC` + `ChFi3d_IntCS` between the small one's tangency
//!    line (`Hctg`) and the big one's bounded surface give the end (`PMil`) of
//!    the intersection curve, and `ChFi3d_ComputeCurves` between the two
//!    surfaces gives `Gc` with its two p-curves; that updates the small stripe,
//! 2. the same curve is retried against the face `F` carried by
//!    `SmaFD->Index(IFaArcSma)`, this time with the p-curve *on that face*
//!    `(cxx:670-745)`, which is what the big stripe needs at its end,
//! 3. both `SurfData` ends, their interferences and the bounding boxes are
//!    written (`cxx:746-801`).
//!
//! Translation notes:
//! - `Big*` / `Sma*` are plain locals in OCCT that alias `HS1` / `HS2`,
//!   `Corner1` / `Corner2` and `Fd1` / `Fd2`; the port clones the two `Arc`
//!   handles (`hs`, `corner`) and borrows the two `SurfData` copies the caller
//!   lends it. `Reduce` mutating `SmaHS` / `BigHS` therefore only touches the
//!   local copies, which is all OCCT reads afterwards,
//! - `Gc`, `PGc1`, `PGc2` and `tolreached` are reused by both
//!   `ChFi3d_ComputeCurves` calls, so the second answer replaces the first,
//! - `done = false; return done;` is `Ok(false)`,
//! - a null handle OCCT dereferences (a tangency line without a curve, a
//!   p-curve missing on face or surf, a stripe whose `Arc` is lost) is the
//!   `Err` of this port.

use std::sync::Arc;

use occt_core::bnd::BndBox;
use occt_core::gp::GpPnt2d;
use occt_geom::curve::Curve as GeomCurve;
use occt_geom::curve_on_surface::CurveOnSurface;
use occt_geom2d::adaptor_trim::TrimmedAdaptor2d;
use occt_geom2d::curve::Curve2d;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_tool::BRepTool;
use crate::chfi3d::builder_0_helpers::{
    ch_fi3d_bound_fac, ch_fi3d_boite, ch_fi3d_compute_curves, ch_fi3d_computes_int_pc,
    ch_fi3d_enlarge_box_curve, ch_fi3d_enlarge_box_edge, ch_fi3d_enlarge_box_pcurve,
    ch_fi3d_enlarge_box_surf_data, ch_fi3d_fil_curve_in_ds, ch_fi3d_fil_point_in_ds,
    ch_fi3d_index_point_in_ds, ch_fi3d_int_cs, ch_fi3d_set_point_tolerance, ComputeCurves,
};
use crate::chfi3d::builder_c2::{reduce_curve, reduce_surface};
use crate::chfi_ds::{StripeHandle, SurfData};
use crate::shape::{Edge, TopoShape};
use crate::topo_tools::{explorer_shapes, is_same};
use crate::topopebrep::ds::Curve as DsCurve;
use crate::topopebrep::tool::adaptor::{GeomAdaptorCurve, SurfaceAdaptor};

use super::super::ChFi3dBuilder;
use super::{face_at, index_point, p_curve_on_surf, surf_value, TwoCorner};

/// `HS->Value(U, V)` of a 2d point.
fn value_at(hs: &SurfaceAdaptor, uv: &GpPnt2d) -> Result<occt_core::gp::GpPnt, String> {
    surf_value(hs, uv.x(), uv.y())
}

/// `cxx:447-804`.
pub(super) fn diff_param(b: &mut ChFi3dBuilder, tc: &mut TwoCorner) -> Result<bool, String> {
    // OCCT's `ChFiDS_SurfData& Fd1` / `Fd2` are references into the stripes'
    // sequences; the port keeps them as locals of this branch and writes them
    // back through `TwoCorner::flush`, so both ends of every `Ok(false)` path
    // still see what the branch wrote before it bailed out.
    let mut fd1 = std::mem::take(&mut tc.fd1);
    let mut fd2 = std::mem::take(&mut tc.fd2);
    let answer = diff_param_body(b, tc, &mut fd1, &mut fd2);
    tc.fd1 = fd1;
    tc.fd2 = fd2;
    answer
}

/// The body proper, with `SmaFD` / `BigFD` standing in for the two `SurfData`.
fn diff_param_body(
    b: &mut ChFi3dBuilder,
    tc: &mut TwoCorner,
    fd1: &mut SurfData,
    fd2: &mut SurfData,
) -> Result<bool, String> {
    let is_first1 = tc.is_first1;
    let is_first2 = tc.is_first2;
    let i_fa_co1 = tc.i_fa_co1;
    let i_fa_co2 = tc.i_fa_co2;
    let i_fa_arc1 = tc.i_fa_arc1;
    let i_fa_arc2 = tc.i_fa_arc2;
    let u_int_pc1 = tc.u_int_pc1;
    let u_int_pc2 = tc.u_int_pc2;

    // `bool parcrois = false; TopExp_Explorer Expl; for (Expl.Init(pivot
    // .Oriented(TopAbs_FORWARD), TopAbs_VERTEX); Expl.More(); Expl.Next()) { if
    // (Expl.Current().IsSame(Vtx)) { parcrois = (Expl.Current().Orientation() ==
    // TopAbs_FORWARD); break; } }` (`cxx:452-463`).
    let mut par_crois = false;
    for v in explorer_shapes(&tc.pivot.oriented(Orientation::Forward), ShapeType::Vertex) {
        if is_same(&v, &tc.vtx) {
            par_crois = v.orientation() == Orientation::Forward;
            break;
        }
    }

    // `if ((parcrois && parCP2 > parCP1) || (!parcrois && parCP2 < parCP1))`
    // (`cxx:472`).
    let big_is_2 = (par_crois && tc.par_cp2 > tc.par_cp1)
        || (!par_crois && tc.par_cp2 < tc.par_cp1);

    // The two assignments of `cxx:472-503`, `Big*` first.
    let big_cd: StripeHandle = if big_is_2 { tc.corner2.clone() } else { tc.corner1.clone() };
    let sma_cd: StripeHandle = if big_is_2 { tc.corner1.clone() } else { tc.corner2.clone() };
    let mut big_hs = if big_is_2 { tc.hs2.clone() } else { tc.hs1.clone() };
    let mut sma_hs = if big_is_2 { tc.hs1.clone() } else { tc.hs2.clone() };
    let (mut u_int_pc_big, mut u_int_pc_sma) = if big_is_2 {
        (u_int_pc2, u_int_pc1)
    } else {
        (u_int_pc1, u_int_pc2)
    };
    let (i_fa_co_big, i_fa_co_sma) = if big_is_2 {
        (i_fa_co2, i_fa_co1)
    } else {
        (i_fa_co1, i_fa_co2)
    };
    let (i_fa_arc_big, i_fa_arc_sma) = if big_is_2 {
        (i_fa_arc2, i_fa_arc1)
    } else {
        (i_fa_arc1, i_fa_arc2)
    };
    let (is_first_big, is_first_sma) = if big_is_2 {
        (is_first2, is_first1)
    } else {
        (is_first1, is_first2)
    };
    let (big_fd, sma_fd): (&mut SurfData, &mut SurfData) =
        if big_is_2 { (fd2, fd1) } else { (fd1, fd2) };

    // Intersection of the big with the small.

    // `ChFi3d_ComputesIntPC(SmaFD->Interference(IFaCoSma), BigFD->Interference(
    // IFaCoBig), SmaHS, BigHS, UIntPCSma, UIntPCBig);` (`cxx:510-516`).
    ch_fi3d_computes_int_pc(
        sma_fd.interference(i_fa_co_sma),
        big_fd.interference(i_fa_co_big),
        &sma_hs,
        &big_hs,
        &mut u_int_pc_sma,
        &mut u_int_pc_big,
    )?;
    // `UVi = BigFD->Interference(IFaCoBig).PCurveOnSurf()->Value(UIntPCBig);
    // Pardeb(3) = UVi.X(); Pardeb(4) = UVi.Y();` (`cxx:517-519`).
    let uv = p_curve_on_surf(big_fd, i_fa_co_big)?.value(u_int_pc_big);
    tc.pardeb[2] = uv.x();
    tc.pardeb[3] = uv.y();
    // `UVi = SmaFD->Interference(IFaCoSma).PCurveOnSurf()->Value(UIntPCSma);
    // Pardeb(1) = UVi.X(); Pardeb(2) = UVi.Y();` (`cxx:520-522`).
    let uv = p_curve_on_surf(sma_fd, i_fa_co_sma)?.value(u_int_pc_sma);
    tc.pardeb[0] = uv.x();
    tc.pardeb[1] = uv.y();
    // `gp_Pnt PFaCo = SmaHS->Value(UVi.X(), UVi.Y());` (`cxx:523`).
    let p_fa_co = value_at(&sma_hs, &uv)?;

    // Parfin (parameters of point PMil).
    // `const ChFiDS_FaceInterference& FiArcSma = SmaFD->Interference(IFaArcSma);`
    // (`cxx:526`).
    let fi_arc_sma = sma_fd.interference(i_fa_arc_sma).clone();
    // `occ::handle<Geom_Curve> ctg = DStr.Curve(FiArcSma.LineIndex()).Curve();`
    // (`cxx:527`).
    let ctg = b
        .my_ds
        .ds()
        .curve(fi_arc_sma.line_index())
        .curve()
        .cloned()
        .ok_or_else(|| {
            "PerformTwoCornerbyInter: the tangency line has no curve".to_string()
        })?;
    // `double temp, wi; if (isfirstSma) { wi = temp = FiArcSma.FirstParameter();
    // if (UIntPCSma < temp) temp = UIntPCSma; bid.Load(ctg, temp,
    // FiArcSma.LastParameter()); } else { wi = temp =
    // FiArcSma.LastParameter(); if (UIntPCSma > temp) temp = UIntPCSma;
    // bid.Load(ctg, FiArcSma.FirstParameter(), temp); }` (`cxx:530-549`).
    let (mut wi, mut hctg) = if is_first_sma {
        let mut temp = fi_arc_sma.first_parameter();
        let wi = temp;
        if u_int_pc_sma < temp {
            temp = u_int_pc_sma;
        }
        (
            wi,
            GeomAdaptorCurve::with_range(Arc::clone(&ctg), temp, fi_arc_sma.last_parameter()),
        )
    } else {
        let mut temp = fi_arc_sma.last_parameter();
        let wi = temp;
        if u_int_pc_sma > temp {
            temp = u_int_pc_sma;
        }
        (
            wi,
            GeomAdaptorCurve::with_range(Arc::clone(&ctg), fi_arc_sma.first_parameter(), temp),
        )
    };
    // `if (SmaFD->Surf() == BigFD->Surf()) { Reduce(UIntPCSma, UIntPCBig, SmaHS,
    // BigHS); Reduce(UIntPCSma, UIntPCBig, Hctg); }` (`cxx:550-554`).
    if sma_fd.surf() == big_fd.surf() {
        reduce_surface(u_int_pc_sma, u_int_pc_big, &mut sma_hs, &mut big_hs);
        reduce_curve(u_int_pc_sma, u_int_pc_big, &mut hctg);
    }
    // `if (!ChFi3d_IntCS(BigHS, Hctg, UVi, wi)) { done = false; return done; }`
    // (`cxx:555-562`).
    let mut uvi = GpPnt2d::new(0.0, 0.0);
    if !ch_fi3d_int_cs(&big_hs, &hctg, &mut uvi, &mut wi) {
        return Ok(false);
    }
    // `Parfin(3) = UVi.X(); Parfin(4) = UVi.Y();` (`cxx:563-564`).
    tc.parfin[2] = uvi.x();
    tc.parfin[3] = uvi.y();
    // `UVi = FiArcSma.PCurveOnSurf()->Value(wi); Parfin(1) = UVi.X();
    // Parfin(2) = UVi.Y();` (`cxx:565-567`).
    let uv = p_curve_on_surf(sma_fd, i_fa_arc_sma)?.value(wi);
    tc.parfin[0] = uv.x();
    tc.parfin[1] = uv.y();
    // `gp_Pnt PMil = SmaHS->Value(Parfin(1), Parfin(2));` (`cxx:568`).
    let p_mil = value_at(&sma_hs, &uv)?;

    // `if (!ChFi3d_ComputeCurves(SmaHS, BigHS, Pardeb, Parfin, Gc, PGc1, PGc2,
    // tolapp3d, tol2d, tolreached)) { done = false; return done; }`
    // (`cxx:570-588`).
    let mut computed = ComputeCurves::default();
    if !ch_fi3d_compute_curves(
        &sma_hs,
        &big_hs,
        tc.pardeb,
        tc.parfin,
        &mut computed,
        b.tolapp3d,
        b.tol2d,
        true,
    ) {
        return Ok(false);
    }
    // `double WFirst = Gc->FirstParameter(); double WLast =
    // Gc->LastParameter();` (`cxx:589-590`).
    let gc: Arc<dyn GeomCurve> = computed.c3d.clone().ok_or_else(|| {
        "PerformTwoCornerbyInter: ChFi3d_ComputeCurves left Gc null".to_string()
    })?;
    let mut w_first = gc.first_parameter();
    let mut w_last = gc.last_parameter();

    // `ChFiDS_CommonPoint& psmaco = SmaFD->ChangeVertex(isfirstSma, IFaCoSma);
    // ChFiDS_CommonPoint& pbigco = BigFD->ChangeVertex(isfirstBig, IFaCoBig);
    // double tolpco = std::max(psmaco.Tolerance(), pbigco.Tolerance());
    // ChFiDS_CommonPoint& psmamil = SmaFD->ChangeVertex(isfirstSma, IFaArcSma);
    // double tolpmil = psmamil.Tolerance();` (`cxx:592-597`).
    let tolpco = sma_fd
        .vertex(is_first_sma, i_fa_co_sma)
        .tolerance()
        .max(big_fd.vertex(is_first_big, i_fa_co_big).tolerance());
    let tolpmil = sma_fd.vertex(is_first_sma, i_fa_arc_sma).tolerance();
    // `int ICurv = DStr.AddCurve(TopOpeBRepDS_Curve(Gc, tolreached));`
    // (`cxx:598`).
    let icurv = b.my_ds.change_ds().add_curve(&mut DsCurve::with_curve(
        Some(Arc::clone(&gc)),
        computed.tolreached,
        false,
    ));

    // `SmaCD->SetParameters(isfirstSma, WFirst, WLast); SmaCD->SetCurve(ICurv,
    // isfirstSma); SmaCD->ChangePCurve(isfirstSma) = PGc1;` (`cxx:599-601`).
    {
        let mut stripe = sma_cd
            .write()
            .map_err(|_| "PerformTwoCornerbyInter: poisoned stripe".to_string())?;
        stripe.set_parameters(is_first_sma, w_first, w_last);
        stripe.set_curve(icurv, is_first_sma);
        *super::p_curve_mut(&mut stripe, is_first_sma) = computed.pc1.clone();
    }
    // `psmaco.Reset(); psmaco.SetPoint(PFaCo);
    // psmaco.SetTolerance(std::max(tolpco, tolreached));
    // SmaFD->ChangeInterference(IFaCoSma).SetParameter(UIntPCSma, isfirstSma);`
    // (`cxx:602-605`).
    {
        let cp = sma_fd.vertex_mut(is_first_sma, i_fa_co_sma);
        cp.reset();
        cp.set_point(&p_fa_co);
        cp.set_tolerance(tolpco.max(computed.tolreached));
    }
    sma_fd
        .interference_mut(i_fa_co_sma)
        .set_parameter(u_int_pc_sma, is_first_sma);
    // `psmamil.Reset(); psmamil.SetPoint(PMil);
    // psmamil.SetTolerance(std::max(tolpmil, tolreached));
    // SmaFD->ChangeInterference(IFaArcSma).SetParameter(wi, isfirstSma);`
    // (`cxx:606-609`).
    {
        let cp = sma_fd.vertex_mut(is_first_sma, i_fa_arc_sma);
        cp.reset();
        cp.set_point(&p_mil);
        cp.set_tolerance(tolpmil.max(computed.tolreached));
    }
    sma_fd
        .interference_mut(i_fa_arc_sma)
        .set_parameter(wi, is_first_sma);

    // `IpointCo = ChFi3d_IndexPointInDS(psmaco, DStr); SmaCD->SetIndexPoint(
    // IpointCo, isfirstSma, IFaCoSma); IpointMil =
    // ChFi3d_IndexPointInDS(psmamil, DStr); SmaCD->SetIndexPoint(IpointMil,
    // isfirstSma, IFaArcSma);` (`cxx:610-614`).
    let ipoint_co = ch_fi3d_index_point_in_ds(
        sma_fd.vertex(is_first_sma, i_fa_co_sma),
        b.my_ds.change_ds(),
    );
    {
        let mut stripe = sma_cd
            .write()
            .map_err(|_| "PerformTwoCornerbyInter: poisoned stripe".to_string())?;
        stripe.set_index_point(ipoint_co, is_first_sma, i_fa_co_sma);
    }
    let ipoint_mil = ch_fi3d_index_point_in_ds(
        sma_fd.vertex(is_first_sma, i_fa_arc_sma),
        b.my_ds.change_ds(),
    );
    {
        let mut stripe = sma_cd
            .write()
            .map_err(|_| "PerformTwoCornerbyInter: poisoned stripe".to_string())?;
        stripe.set_index_point(ipoint_mil, is_first_sma, i_fa_arc_sma);
    }
    // `if (IFaCoSma == 2) { SmaCD->SetOrientation(TopAbs_REVERSED, isfirstSma); }`
    // (`cxx:615-618`).
    if i_fa_co_sma == 2 {
        let mut stripe = sma_cd
            .write()
            .map_err(|_| "PerformTwoCornerbyInter: poisoned stripe".to_string())?;
        stripe.set_orientation(Orientation::Reversed, is_first_sma);
    }

    // `BigCD->SetIndexPoint(IpointCo, isfirstBig, IFaCoBig);
    // BigFD->ChangeVertex(isfirstBig, IFaCoBig) = psmaco;
    // BigFD->ChangeInterference(IFaCoBig).SetParameter(UIntPCBig, isfirstBig);`
    // (`cxx:620-623`).
    {
        let mut stripe = big_cd
            .write()
            .map_err(|_| "PerformTwoCornerbyInter: poisoned stripe".to_string())?;
        stripe.set_index_point(ipoint_co, is_first_big, i_fa_co_big);
    }
    let psma_co = sma_fd.vertex(is_first_sma, i_fa_co_sma).clone();
    *big_fd.vertex_mut(is_first_big, i_fa_co_big) = psma_co;
    big_fd
        .interference_mut(i_fa_co_big)
        .set_parameter(u_int_pc_big, is_first_big);

    // `Interfp = ChFi3d_FilPointInDS(TopAbs_FORWARD, ICurv, IpointCo, WFirst);
    // Li.Append(Interfp); Interfp = ChFi3d_FilPointInDS(TopAbs_REVERSED, ICurv,
    // IpointMil, WLast); Li.Append(Interfp);` (`cxx:625-630`).
    {
        let li = b.my_ds.change_ds().change_curve_interferences(icurv);
        li.push(ch_fi3d_fil_point_in_ds(
            Orientation::Forward,
            icurv,
            ipoint_co,
            w_first,
            false,
        ));
    }
    {
        let li = b.my_ds.change_ds().change_curve_interferences(icurv);
        li.push(ch_fi3d_fil_point_in_ds(
            Orientation::Reversed,
            icurv,
            ipoint_mil,
            w_last,
            false,
        ));
    }

    // `TopAbs_Orientation tra = BigFD->InterferenceOnS1().Transition();
    // TopAbs_Orientation ofac = DStr.Shape(BigFD->IndexOfS1()).Orientation();
    // TopAbs_Orientation ofil = BigFD->Orientation(); TopAbs_Orientation tracurv
    // = TopAbs::Compose(ofac, ofil); tracurv = TopAbs::Compose(tracurv, tra); if
    // (!isfirstBig) tracurv = TopAbs::Reverse(tracurv); if (IFaCoBig != 1)
    // tracurv = TopAbs::Reverse(tracurv);` (`cxx:633-644`).
    let tra = big_fd.interference_on_s1().transition();
    let ofac = b
        .my_ds
        .ds()
        .shape(big_fd.index_of_s1(), true)
        .ok_or_else(|| {
            "PerformTwoCornerbyInter: DStr.Shape(IndexOfS1) is null".to_string()
        })?
        .orientation();
    let ofil = big_fd.orientation();
    let mut tracurv = Orientation::compose(ofac, ofil);
    tracurv = Orientation::compose(tracurv, tra);
    if !is_first_big {
        tracurv = tracurv.reversed();
    }
    if i_fa_co_big != 1 {
        tracurv = tracurv.reversed();
    }

    // `int ISurf = BigFD->Surf(); Interfc = ChFi3d_FilCurveInDS(ICurv, ISurf,
    // PGc2, tracurv); DStr.ChangeSurfaceInterferences(ISurf).Append(Interfc);`
    // (`cxx:646-650`).
    let isurf = big_fd.surf();
    {
        let interfc = ch_fi3d_fil_curve_in_ds(icurv, isurf, computed.pc2.clone(), tracurv);
        b.my_ds
            .change_ds()
            .change_surface_interferences(isurf)
            .push(interfc);
    }

    // The tolerances of points are updated (beginning).
    // `Bnd_Box bco, bmil, barc; if (IFaCoSma == 1) ChFi3d_EnlargeBox(DStr, SmaCD,
    // SmaFD, bco, bmil, isfirstSma); else ChFi3d_EnlargeBox(DStr, SmaCD, SmaFD,
    // bmil, bco, isfirstSma); ChFi3d_EnlargeBox(BigHS, PGc2, WFirst, WLast, bco,
    // bmil);` (`cxx:653-663`).
    let mut bco = BndBox::new();
    let mut bmil = BndBox::new();
    let mut barc = BndBox::new();
    {
        let stripe = sma_cd
            .read()
            .map_err(|_| "PerformTwoCornerbyInter: poisoned stripe".to_string())?;
        let ds = b.my_ds.change_ds();
        if i_fa_co_sma == 1 {
            ch_fi3d_enlarge_box_surf_data(
                ds,
                Some(&stripe),
                sma_fd,
                &mut bco,
                &mut bmil,
                is_first_sma,
            )?;
        } else {
            ch_fi3d_enlarge_box_surf_data(
                ds,
                Some(&stripe),
                sma_fd,
                &mut bmil,
                &mut bco,
                is_first_sma,
            )?;
        }
    }
    ch_fi3d_enlarge_box_pcurve(
        big_hs
            .surface()
            .ok_or_else(|| "PerformTwoCornerbyInter: BigHS is null".to_string())?,
        computed
            .pc2
            .as_ref()
            .ok_or_else(|| "PerformTwoCornerbyInter: PGc2 is null".to_string())?,
        w_first,
        w_last,
        &mut bco,
        &mut bmil,
    );

    // Intersection of the big with the face at end.

    // `TopoDS_Face F = TopoDS::Face(DStr.Shape(SmaFD->Index(IFaArcSma)));
    // occ::handle<BRepAdaptor_Surface> HF = new BRepAdaptor_Surface(F);`
    // (`cxx:670-671`).
    let f = face_at(b.my_ds.ds(), sma_fd.index(i_fa_arc_sma))?;
    // `double fsma = FiArcSma.FirstParameter(); double lsma =
    // FiArcSma.LastParameter(); double deltSma = 0.05 * (lsma - fsma);`
    // (`cxx:672-674`).
    let mut fsma = fi_arc_sma.first_parameter();
    let mut lsma = fi_arc_sma.last_parameter();
    let delt_sma = 0.05 * (lsma - fsma);
    // `occ::handle<Geom2d_Curve> pcpc = SmaFD->Interference(IFaArcSma)
    // .PCurveOnFace(); fsma = std::max(pcpc->FirstParameter(), wi - deltSma);
    // lsma = std::min(pcpc->LastParameter(), wi + deltSma);` (`cxx:675-678`).
    let pcpc = fi_arc_sma.p_curve_on_face().cloned().ok_or_else(|| {
        "PerformTwoCornerbyInter: the arc interference has no p-curve on face".to_string()
    })?;
    fsma = pcpc.first_parameter().max(wi - delt_sma);
    lsma = pcpc.last_parameter().min(wi + delt_sma);
    // `if (lsma < fsma) { done = false; return done; }` (`cxx:679-683`).
    if lsma < fsma {
        return Ok(false);
    }
    // `occ::handle<Geom2dAdaptor_Curve> c2df = new Geom2dAdaptor_Curve(
    // SmaFD->Interference(IFaArcSma).PCurveOnFace(), fsma, lsma);
    // Adaptor3d_CurveOnSurface consf(c2df, HF);` (`cxx:684-685`).
    let hf = SurfaceAdaptor::of_face(&f);
    let c2df: Arc<dyn Curve2d> = Arc::new(TrimmedAdaptor2d::new(Arc::clone(&pcpc), fsma, lsma));
    let hf_surface = hf
        .surface()
        .cloned()
        .ok_or_else(|| "PerformTwoCornerbyInter: HF is null".to_string())?;
    let consf = CurveOnSurface::new(c2df, hf_surface);
    // `if (!ChFi3d_IntCS(BigHS, Hconsf, UVi, wi)) { done = false; return done; }`
    // (`cxx:687-695`).
    if !ch_fi3d_int_cs(&big_hs, &consf, &mut uvi, &mut wi) {
        return Ok(false);
    }
    // `Pardeb(3) = UVi.X(); Pardeb(4) = UVi.Y();` (`cxx:696-697`).
    tc.pardeb[2] = uvi.x();
    tc.pardeb[3] = uvi.y();
    // `UVi = SmaFD->Interference(IFaArcSma).PCurveOnFace()->Value(wi); Pardeb(1)
    // = UVi.X(); Pardeb(2) = UVi.Y(); gp_Pnt2d ppff1 = UVi;` (`cxx:698-700`).
    let uv = pcpc.value(wi);
    tc.pardeb[0] = uv.x();
    tc.pardeb[1] = uv.y();
    let ppff1 = uv;

    // Parfin (parameters of the point cpend).
    // `double ptg = BigFD->Interference(IFaArcBig).Parameter(isfirstBig); UVi =
    // BigFD->Interference(IFaArcBig).PCurveOnSurf()->Value(ptg); Parfin(3) =
    // UVi.X(); Parfin(4) = UVi.Y();` (`cxx:703-707`).
    let ptg = big_fd.interference(i_fa_arc_big).parameter(is_first_big);
    let uv = p_curve_on_surf(big_fd, i_fa_arc_big)?.value(ptg);
    tc.parfin[2] = uv.x();
    tc.parfin[3] = uv.y();
    // `ChFiDS_CommonPoint& cpend = BigFD->ChangeVertex(isfirstBig, IFaArcBig);
    // TopoDS_Edge etest = cpend.Arc(); if (BRep_Tool::IsClosed(etest, F))
    // etest.Reverse();` (`cxx:708-712`).
    let cpend = big_fd.vertex(is_first_big, i_fa_arc_big).clone();
    let mut etest = cpend.arc().cloned().ok_or_else(|| {
        "PerformTwoCornerbyInter: cpend has no arc".to_string()
    })?;
    if BRepTool::is_closed_edge_face(&Edge(etest.clone()), &f) {
        etest.set_orientation(etest.orientation().reversed());
    }
    // `BRepAdaptor_Curve2d arc(etest, F); UVi = arc.Value(cpend
    // .ParameterOnArc()); Parfin(1) = UVi.X(); Parfin(2) = UVi.Y(); gp_Pnt2d
    // ppff2 = UVi;` (`cxx:713-718`).
    let par_arc = cpend.parameter_on_arc().ok_or_else(|| {
        "PerformTwoCornerbyInter: cpend has no parameter on arc".to_string()
    })?;
    let (edge_pcurve, _, _) = curve_on_surface_range(&Edge(etest.clone()), &f).ok_or_else(|| {
        "PerformTwoCornerbyInter: etest has no p-curve on F".to_string()
    })?;
    let uv = edge_pcurve.value(par_arc);
    tc.parfin[0] = uv.x();
    tc.parfin[1] = uv.y();
    let ppff2 = uv;

    // Intersection.
    // `ChFi3d_Boite(ppff1, ppff2, uu1, uu2, vv1, vv2);
    // ChFi3d_BoundFac(*HF, uu1, uu2, vv1, vv2, true);` (`cxx:720-727`).
    let boite = ch_fi3d_boite(&ppff1, &ppff2);
    let hf = ch_fi3d_bound_fac(&hf, boite.u_min, boite.u_max, boite.v_min, boite.v_max);

    // `if (!ChFi3d_ComputeCurves(HF, BigHS, Pardeb, Parfin, Gc, PGc1, PGc2,
    // tolapp3d, tol2d, tolreached)) { done = false; return done; }`
    // (`cxx:729-744`). OCCT reuses the four out-parameters, so the second answer
    // replaces the first.
    let mut computed = ComputeCurves::default();
    if !ch_fi3d_compute_curves(
        &hf,
        &big_hs,
        tc.pardeb,
        tc.parfin,
        &mut computed,
        b.tolapp3d,
        b.tol2d,
        true,
    ) {
        return Ok(false);
    }

    // End of update of the BigCD and the DS.
    // `WFirst = Gc->FirstParameter(); WLast = Gc->LastParameter(); ICurv =
    // DStr.AddCurve(TopOpeBRepDS_Curve(Gc, tolreached));` (`cxx:747-748`).
    let gc: Arc<dyn GeomCurve> = computed.c3d.clone().ok_or_else(|| {
        "PerformTwoCornerbyInter: ChFi3d_ComputeCurves left Gc null".to_string()
    })?;
    w_first = gc.first_parameter();
    w_last = gc.last_parameter();
    let icurv = {
        let ds = b.my_ds.change_ds();
        ds.add_curve(&mut DsCurve::with_curve(
            Some(Arc::clone(&gc)),
            computed.tolreached,
            false,
        ))
    };
    // `cpend.SetTolerance(std::max(cpend.Tolerance(), tolreached));` (`cxx:749`).
    {
        let tol = big_fd.vertex(is_first_big, i_fa_arc_big).tolerance();
        big_fd
            .vertex_mut(is_first_big, i_fa_arc_big)
            .set_tolerance(tol.max(computed.tolreached));
    }
    // `IpointArc = ChFi3d_IndexPointInDS(cpend, DStr); BigCD->SetIndexPoint(
    // IpointArc, isfirstBig, IFaArcBig);` (`cxx:750-751`).
    let ipoint_arc = ch_fi3d_index_point_in_ds(
        big_fd.vertex(is_first_big, i_fa_arc_big),
        b.my_ds.change_ds(),
    );
    {
        let mut stripe = big_cd
            .write()
            .map_err(|_| "PerformTwoCornerbyInter: poisoned stripe".to_string())?;
        stripe.set_index_point(ipoint_arc, is_first_big, i_fa_arc_big);
    }
    // `Li7.Append(ChFi3d_FilPointInDS(TopAbs_FORWARD, ICurv, IpointMil, WFirst));
    // Li7.Append(ChFi3d_FilPointInDS(TopAbs_REVERSED, ICurv, IpointArc, WLast));`
    // (`cxx:755-758`).
    {
        let li7 = b.my_ds.change_ds().change_curve_interferences(icurv);
        li7.push(ch_fi3d_fil_point_in_ds(
            Orientation::Forward,
            icurv,
            ipoint_mil,
            w_first,
            false,
        ));
    }
    {
        let li7 = b.my_ds.change_ds().change_curve_interferences(icurv);
        li7.push(ch_fi3d_fil_point_in_ds(
            Orientation::Reversed,
            icurv,
            ipoint_arc,
            w_last,
            false,
        ));
    }
    // `Interfc = ChFi3d_FilCurveInDS(ICurv, ISurf, PGc2, tracurv);
    // DStr.ChangeSurfaceInterferences(ISurf).Append(Interfc);
    // BigCD->InDS(isfirstBig);` (`cxx:759-761`).
    {
        let interfc = ch_fi3d_fil_curve_in_ds(icurv, isurf, computed.pc2.clone(), tracurv);
        b.my_ds
            .change_ds()
            .change_surface_interferences(isurf)
            .push(interfc);
    }
    {
        let mut stripe = big_cd
            .write()
            .map_err(|_| "PerformTwoCornerbyInter: poisoned stripe".to_string())?;
        stripe.in_ds(is_first_big, 1);
    }

    // Finally the information on faces is placed in the DS.
    // `int IShape = DStr.AddShape(F);` (`cxx:767`).
    let ishape = b.my_ds.change_ds().add_shape(&f.0);
    // `if (SmaFD->Surf() == BigFD->Surf()) tracurv =
    // TopAbs::Compose(etest.Orientation(), cpend.TransitionOnArc()); else { for
    // (Exp.Init(F.Oriented(TopAbs_FORWARD), TopAbs_EDGE); Exp.More(); Exp.Next())
    // if (Exp.Current().IsSame(etest)) { tracurv =
    // TopAbs::Compose(Exp.Current().Orientation(),
    // cpend.TransitionOnArc()); break; } }` (`cxx:765-780`).
    let cpend_tra = big_fd
        .vertex(is_first_big, i_fa_arc_big)
        .transition_on_arc()
        .ok_or_else(|| {
            "PerformTwoCornerbyInter: cpend has no transition on arc".to_string()
        })?;
    if sma_fd.surf() == big_fd.surf() {
        tracurv = Orientation::compose(etest.orientation(), cpend_tra);
    } else {
        for e in explorer_shapes(&f.0.oriented(Orientation::Forward), ShapeType::Edge) {
            if is_same(&e, &etest) {
                tracurv = Orientation::compose(e.orientation(), cpend_tra);
                break;
            }
        }
    }
    // `Interfc = ChFi3d_FilCurveInDS(ICurv, IShape, PGc1, tracurv);
    // DStr.ChangeShapeInterferences(IShape).Append(Interfc);` (`cxx:781-782`).
    {
        let interfc = ch_fi3d_fil_curve_in_ds(icurv, ishape, computed.pc1.clone(), tracurv);
        b.my_ds
            .change_ds()
            .change_shape_interferences(&f.0)
            .push(interfc);
    }

    // The tolerances of points are updated (end).
    // `occ::handle<ChFiDS_Stripe> bidst; if (IFaCoBig == 1) ChFi3d_EnlargeBox(
    // DStr, bidst, BigFD, bco, barc, isfirstBig); else
    // ChFi3d_EnlargeBox(DStr, bidst, BigFD, barc, bco, isfirstBig);`
    // (`cxx:785-790`). `bidst` is a null handle.
    {
        let ds = b.my_ds.change_ds();
        if i_fa_co_big == 1 {
            ch_fi3d_enlarge_box_surf_data(ds, None, big_fd, &mut bco, &mut barc, is_first_big)?;
        } else {
            ch_fi3d_enlarge_box_surf_data(ds, None, big_fd, &mut barc, &mut bco, is_first_big)?;
        }
    }
    // `ChFi3d_EnlargeBox(BigHS, PGc2, WFirst, WLast, bmil, barc);
    // ChFi3d_EnlargeBox(HF, PGc1, WFirst, WLast, bmil, barc);
    // ChFi3d_EnlargeBox(Gc, WFirst, WLast, bmil, barc);` (`cxx:791-793`).
    let pc2 = computed
        .pc2
        .as_ref()
        .ok_or_else(|| "PerformTwoCornerbyInter: PGc2 is null".to_string())?;
    let pc1 = computed
        .pc1
        .as_ref()
        .ok_or_else(|| "PerformTwoCornerbyInter: PGc1 is null".to_string())?;
    ch_fi3d_enlarge_box_pcurve(
        big_hs
            .surface()
            .ok_or_else(|| "PerformTwoCornerbyInter: BigHS is null".to_string())?,
        pc2,
        w_first,
        w_last,
        &mut bmil,
        &mut barc,
    );
    ch_fi3d_enlarge_box_pcurve(
        hf.surface()
            .ok_or_else(|| "PerformTwoCornerbyInter: HF is null".to_string())?,
        pc1,
        w_first,
        w_last,
        &mut bmil,
        &mut barc,
    );
    ch_fi3d_enlarge_box_curve(&gc, w_first, w_last, &mut bmil, &mut barc);
    // `const ChFiDS_CommonPoint& cparc = BigFD->Vertex(isfirstBig, IFaArcBig);
    // ChFi3d_EnlargeBox(cparc.Arc(), myEFMap(cparc.Arc()),
    // cparc.ParameterOnArc(), barc);` (`cxx:794-796`).
    let cparc = big_fd.vertex(is_first_big, i_fa_arc_big).clone();
    let arc = cparc.arc().cloned().ok_or_else(|| {
        "PerformTwoCornerbyInter: cparc has no arc".to_string()
    })?;
    let par = cparc.parameter_on_arc().ok_or_else(|| {
        "PerformTwoCornerbyInter: cparc has no parameter on arc".to_string()
    })?;
    let faces: Vec<TopoShape> = b
        .my_ef_map
        .get_from_key(&arc)
        .to_vec();
    ch_fi3d_enlarge_box_edge(&arc, &faces, par, &mut barc);

    // `ChFi3d_SetPointTolerance(DStr, bco, SmaCD->IndexPoint(isfirstSma,
    // IFaCoSma)); ChFi3d_SetPointTolerance(DStr, bmil, SmaCD->IndexPoint(
    // isfirstSma, IFaArcSma)); ChFi3d_SetPointTolerance(DStr, barc,
    // BigCD->IndexPoint(isfirstBig, IFaArcBig));` (`cxx:800-802`).
    let ip_co_sma = index_point(&sma_cd, is_first_sma, i_fa_co_sma)?;
    ch_fi3d_set_point_tolerance(b.my_ds.change_ds(), &bco, ip_co_sma);
    let ip_arc_sma = index_point(&sma_cd, is_first_sma, i_fa_arc_sma)?;
    ch_fi3d_set_point_tolerance(b.my_ds.change_ds(), &bmil, ip_arc_sma);
    let ip_arc_big = index_point(&big_cd, is_first_big, i_fa_arc_big)?;
    ch_fi3d_set_point_tolerance(b.my_ds.change_ds(), &barc, ip_arc_big);

    // `done = true; return done;` (`cxx:803-805`).
    Ok(true)
}
