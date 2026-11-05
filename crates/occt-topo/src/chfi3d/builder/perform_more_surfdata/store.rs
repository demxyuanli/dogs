//! `PerformMoreSurfdata` — the DS storage block.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_C1.cxx:4158-4521`
//! (OCCT 8.0.0).
//!
//! Six interferences are appended here, in source order:
//!
//! 1. the curve point `aCP1onArc` on its own arc (`cxx:4198-4201`),
//! 2. a `TopOpeBRepDS_SolidSurfaceInterference` of the stripe's solid on the
//!    surface of `anInd` (`cxx:4203-4214`), followed by the removal of
//!    `SeqFil(anInd)` (`cxx:4216`) — after which `Value(anInd)` is the old
//!    `anIndPrev` entry, which is why every later `aSurfData->Surf()` is the
//!    neighbour's surface,
//! 3. `aCint2` on the neighbour's surface and on `aFace` (`cxx:4291-4343`),
//! 4. `aTrCracc` on the neighbour's surface and on the surface of `anInd`
//!    (`cxx:4386-4408`),
//! 5. `aCint1` on the surface of `anInd` and on `aFace` (`cxx:4412-4455`),
//! 6. `aCline` on the surface of `anInd` and on `aNeighborFace`
//!    (`cxx:4464-4521`).
//!
//! Translation notes:
//!
//! - `anOrSurf` starts as `aCurOrient` (`cxx:4282`) and is reversed by the
//!   orientation probes of `cxx:4317`, `:4380`, `:4402` and `:4433`; the port
//!   keeps it as one `mut` local so the accumulation is visible,
//! - `aPoint1` / `aPoint2` come from the `aFI` of `anInd` (`cxx:4288-4291`)
//!   through `aFI.PCurveOnFace()` and the face's surface,
//! - the `aTrCline` point choice at `cxx:4473-4483` uses a strict `<` while the
//!   two earlier ones use `<=`; the port keeps the difference.

use std::sync::Arc;

use occt_geom::curve::Curve;

use crate::abs::Orientation;

use crate::chfi3d::builder_0_helpers::{
    ch_fi3d_fil_curve_in_ds, ch_fi3d_fil_point_in_ds, ch_fi3d_index_point_in_ds,
};
use crate::shape::Face;
use crate::topopebrep::ds::{Curve as DsCurve, Interference, Kind, Point, Transition};

use super::super::ChFi3dBuilder;
use super::{sections::Section, MoreSurfdata, TOL3D};

/// `cxx:4158-4521`. `Ok(())` is the `if (indShape <= 0) return;` too.
pub(super) fn write(
    b: &mut ChFi3dBuilder,
    st: &mut MoreSurfdata,
    sec: &Section,
) -> Result<(), String> {
    let tol3d = TOL3D;
    // `aSurf` still holds `DStr.Surface(SeqFil(anIndPrev)->Surf()).Surface()`
    // from `cxx:4086`.
    let sd_prev = super::fd_at(st, st.ind_prev)?;
    let a_surf = b
        .my_ds
        .ds()
        .surface(sd_prev.surf())
        .surface()
        .cloned()
        .ok_or_else(|| "PerformMoreSurfdata: DStr.Surface(anIndPrev) is null".to_string())?;

    // `occ::handle<Geom2d_Curve> aPCraccS = GeomProjLib::Curve2d(aTrCracc, aSurf);`
    // (`cxx:4163`).
    let pcracc_s = occt_geom::geom_proj_lib::curve_2d(&sec.trcracc, &a_surf);

    // `if (is2ndCP1OnArc) { aFI = aSurfData->InterferenceOnS1(); indShape =
    //  aSurfData->IndexOfS1(); } else { ... S2 ... }` (`cxx:4165-4176`).
    let fi_prev = if st.is2nd_cp1_on_arc {
        sd_prev.interference_on_s1().clone()
    } else {
        sd_prev.interference_on_s2().clone()
    };
    let ind_shape_prev = if st.is2nd_cp1_on_arc {
        sd_prev.index_of_s1()
    } else {
        sd_prev.index_of_s2()
    };
    // `if (indShape <= 0) return;` (`cxx:4178-4181`).
    if ind_shape_prev <= 0 {
        return Ok(());
    }

    // `aCurOrient = DStr.Shape(indShape).Orientation(); aCurOrient =
    //  TopAbs::Compose(aCurOrient, aSurfData->Orientation()); aCurOrient =
    //  TopAbs::Compose(TopAbs::Reverse(aFI.Transition()), aCurOrient);`
    // (`cxx:4182-4186`).
    let shape_prev = b
        .my_ds
        .ds()
        .shape(ind_shape_prev, true)
        .ok_or_else(|| "PerformMoreSurfdata: DStr.Shape(indShape) is null".to_string())?;
    let mut cur_orient = shape_prev.orientation();
    cur_orient = Orientation::compose(cur_orient, sd_prev.orientation());
    cur_orient = Orientation::compose(fi_prev.transition().reversed(), cur_orient);

    // `aSurfData = aSeqSurfData.Value(anInd); ... indSurf1 = aSurfData->Surf();
    //  indArc1 = DStr.AddShape(aCP1onArc.Arc()); indSol = aStripe->SolidIndex();`
    // (`cxx:4189-4194`).
    let fd_ind = super::fd_of(st)?;
    let ind_surf1 = fd_ind.surf();
    let pt_cp1 = Point::with_tolerance(st.cp1_on_arc.point().clone(), st.cp1_on_arc.tolerance());
    let ind_cp1_on_arc = b.my_ds.change_ds().add_point(&pt_cp1);
    let arc1 = st
        .cp1_on_arc
        .arc()
        .cloned()
        .ok_or_else(|| "PerformMoreSurfdata: aCP1onArc has no arc".to_string())?;
    let ind_arc1 = b.my_ds.change_ds().add_shape(&arc1);
    let ind_sol = {
        let Ok(guard) = st.stripe.read() else {
            return Err("PerformMoreSurfdata: aStripe is locked".to_string());
        };
        guard.solid_index()
    };

    // `anInterfp1 = ChFi3d_FilPointInDS(aCP1onArc.TransitionOnArc(), indArc1,
    //  indCP1onArc, aCP1onArc.ParameterOnArc()); DStr.ChangeShapeInterferences(
    //  aCP1onArc.Arc()).Append(anInterfp1);` (`cxx:4198-4201`).
    let t_arc = st
        .cp1_on_arc
        .transition_on_arc()
        .ok_or_else(|| "PerformMoreSurfdata: aCP1onArc has no transition".to_string())?;
    let p_arc = st
        .cp1_on_arc
        .parameter_on_arc()
        .ok_or_else(|| "PerformMoreSurfdata: aCP1onArc has no parameter".to_string())?;
    let interf_p1 = ch_fi3d_fil_point_in_ds(t_arc, ind_arc1, ind_cp1_on_arc, p_arc, false);
    {
        let ds = b.my_ds.change_ds();
        ds.change_shape_interferences_at(ind_arc1).push(interf_p1);
    }

    // `SolidInterfs = DStr.ChangeShapeInterferences(indSol); SSI = new
    //  TopOpeBRepDS_SolidSurfaceInterference(TopOpeBRepDS_Transition(anOrSD1),
    //  SOLID, indSol, SURFACE, indSurf1); SolidInterfs.Append(SSI);`
    // (`cxx:4203-4214`).
    {
        let ssi = std::sync::Arc::new(std::sync::RwLock::new(Interference::solid_surface(
            Transition::from_orientation(sec.or_sd1),
            Kind::Solid,
            ind_sol,
            Kind::Surface,
            ind_surf1,
        )));
        let ds = b.my_ds.change_ds();
        ds.change_shape_interferences_at(ind_sol).push(ssi);
    }

    // `aSeqSurfData.Remove(anInd); if (!isFirst) anInd--; aSurfData =
    //  aSeqSurfData.Value(anInd);` (`cxx:4216-4224`).
    let pos = st
        .index
        .checked_sub(1)
        .ok_or_else(|| "PerformMoreSurfdata: anInd is not a sequence index".to_string())?;
    if pos >= st.seq.len() {
        return Err("PerformMoreSurfdata: aSeqSurfData.Remove(anInd) is out of range".to_string());
    }
    st.seq.remove(pos);
    if !st.is_first {
        st.index -= 1;
    }
    let fd_prev = super::fd_of(st)?;

    // `cxx:4230-4261`: the two point indices of the neighbour, on the side
    // `is2ndCP1OnArc` puts them.
    let (ind_cp2_on_arc, ind_cp2_not_on_arc) = {
        let i1 = ch_fi3d_index_point_in_ds(&st.cp2_on_arc, b.my_ds.change_ds());
        let i2 = ch_fi3d_index_point_in_ds(&st.cp2_not_on_arc, b.my_ds.change_ds());
        if st.is2nd_cp1_on_arc {
            if let Ok(mut stripe) = st.stripe.write() {
                stripe.set_index_point(i1, st.is_first, 1);
                stripe.set_index_point(i2, st.is_first, 2);
            }
            let Ok(guard) = st.stripe.read() else {
                return Err("PerformMoreSurfdata: aStripe is locked".to_string());
            };
            if st.is_first {
                (
                    guard.index_first_point_on_s1(),
                    guard.index_first_point_on_s2(),
                )
            } else {
                (
                    guard.index_last_point_on_s1(),
                    guard.index_last_point_on_s2(),
                )
            }
        } else {
            if let Ok(mut stripe) = st.stripe.write() {
                stripe.set_index_point(i1, st.is_first, 2);
                stripe.set_index_point(i2, st.is_first, 1);
            }
            let Ok(guard) = st.stripe.read() else {
                return Err("PerformMoreSurfdata: aStripe is locked".to_string());
            };
            if st.is_first {
                (
                    guard.index_first_point_on_s2(),
                    guard.index_first_point_on_s1(),
                )
            } else {
                (
                    guard.index_last_point_on_s2(),
                    guard.index_last_point_on_s1(),
                )
            }
        }
    };

    // `if (is2ndCP1OnArc) { aFI = aSurfData->InterferenceOnS1(); indShape =
    //  aSurfData->IndexOfS1(); } else { ... S2 ... }` (`cxx:4269-4278`).
    let fi_ind = if st.is2nd_cp1_on_arc {
        fd_prev.interference_on_s1().clone()
    } else {
        fd_prev.interference_on_s2().clone()
    };
    let ind_shape_ind = if st.is2nd_cp1_on_arc {
        fd_prev.index_of_s1()
    } else {
        fd_prev.index_of_s2()
    };

    // `TopAbs_Orientation anOrSurf = aCurOrient; anOrFace = aFace.Orientation();
    //  indaFace = DStr.AddShape(aFace); indPoint = indCP2onArc;`
    // (`cxx:4282-4285`).
    let mut or_surf = cur_orient;
    let mut or_face = st.face.orientation();
    let inda_face = b.my_ds.change_ds().add_shape(&st.face.0);

    // `aFI.PCurveOnFace()->D0(aFI.LastParameter(), aP2d); Stemp2 =
    //  BRep_Tool::Surface(TopoDS::Face(DStr.Shape(indShape))); Stemp2->D0(...,
    //  aPoint2); aFI.PCurveOnFace()->D0(aFI.FirstParameter(), aP2d); ...
    //  aPoint1;` (`cxx:4287-4291`).
    let stemp2 = {
        let ds = b.my_ds.ds();
        let s = ds
            .shape(ind_shape_ind, true)
            .ok_or_else(|| "PerformMoreSurfdata: DStr.Shape(indShape) is null".to_string())?;
        let f = Face::wrap(s)
            .ok_or_else(|| "PerformMoreSurfdata: DStr.Shape(indShape) is not a face".to_string())?;
        crate::brep_tool::BRepTool::face_surface(&f)
            .ok_or_else(|| "PerformMoreSurfdata: BRep_Tool::Surface(Stemp2) is null".to_string())?
    };
    let pc_face = fi_ind
        .p_curve_on_face()
        .cloned()
        .ok_or_else(|| "PerformMoreSurfdata: aFI has no p-curve on face".to_string())?;
    let p2d_last = pc_face.value(fi_ind.last_parameter());
    let point2 = stemp2.d0(p2d_last.x(), p2d_last.y());
    let p2d_first = pc_face.value(fi_ind.first_parameter());
    let point1 = stemp2.d0(p2d_first.x(), p2d_first.y());

    // `int indPoint = indCP2onArc;` (`cxx:4285`).
    let mut ind_point = ind_cp2_on_arc;
    // `int indPoint1; int indPoint2;` (`cxx:4264-4265`): assigned before every
    // read, on both paths.
    let mut ind_point1: i32;
    let mut ind_point2: i32;

    // `if (isDoSecondSection) { ... }` (`cxx:4293-4343`).
    if sec.is_do_second_section {
        let cint2 = sec
            .cint2
            .clone()
            .ok_or_else(|| "PerformMoreSurfdata: isDoSecondSection without aCint2".to_string())?;
        let pcint21 = sec
            .pcint21
            .clone()
            .ok_or_else(|| "PerformMoreSurfdata: isDoSecondSection without aPCint21".to_string())?;
        let pcint22 = sec
            .pcint22
            .clone()
            .ok_or_else(|| "PerformMoreSurfdata: isDoSecondSection without aPCint22".to_string())?;

        // `TopOpeBRepDS_Point tpoint(aPext, aTolex2); TopOpeBRepDS_Curve
        //  tcint2(aCint2, aTolex2); indPoint = DStr.AddPoint(tpoint); indCurve =
        //  DStr.AddCurve(tcint2);` (`cxx:4295-4299`).
        let tpoint = Point::with_tolerance(sec.pext.clone(), sec.tolex2);
        ind_point = b.my_ds.change_ds().add_point(&tpoint);
        let mut tcint2 = DsCurve::with_curve(Some(Arc::clone(&cint2)), sec.tolex2, false);
        let ind_curve = b.my_ds.change_ds().add_curve(&mut tcint2);

        // `aCint2->D0(First/LastParameter(), aPext1 / aPext2);` (`cxx:4301-4302`).
        let pext1 = cint2.d0(cint2.first_parameter());
        let pext2 = cint2.d0(cint2.last_parameter());
        // `cxx:4304-4313`.
        if pext1.distance(&sec.pext) <= tol3d {
            ind_point1 = ind_point;
            ind_point2 = ind_cp2_on_arc;
        } else {
            ind_point1 = ind_cp2_on_arc;
            ind_point2 = ind_point;
        }
        // `cxx:4316-4319`.
        if pext1.distance(&point2) > tol3d && pext2.distance(&point1) > tol3d {
            or_surf = or_surf.reversed();
        }
        // `cxx:4323-4327`.
        let interfp1 = ch_fi3d_fil_point_in_ds(
            Orientation::Forward,
            ind_curve,
            ind_point1,
            cint2.first_parameter(),
            false,
        );
        let interfp2 = ch_fi3d_fil_point_in_ds(
            Orientation::Reversed,
            ind_curve,
            ind_point2,
            cint2.last_parameter(),
            false,
        );
        {
            let ds = b.my_ds.change_ds();
            let list = ds.change_curve_interferences(ind_curve);
            list.push(interfp1);
            list.push(interfp2);
        }
        // `anInterfc = ChFi3d_FilCurveInDS(indCurve, aSurfData->Surf(), aPCint22,
        //  anOrSurf); DStr.ChangeSurfaceInterferences(aSurfData->Surf()).Append(
        //  anInterfc);` (`cxx:4329-4331`).
        let interfc = ch_fi3d_fil_curve_in_ds(
            ind_curve,
            fd_prev.surf(),
            Some(pcint22),
            or_surf,
        );
        b.my_ds
            .change_ds()
            .change_surface_interferences(fd_prev.surf())
            .push(interfc);
        // `if (anOrFace == anOrSD2) anOrFace = TopAbs::Reverse(anOrSurf); else
        //  anOrFace = anOrSurf;` (`cxx:4334-4342`).
        or_face = if or_face == sec.or_sd2 {
            or_surf.reversed()
        } else {
            or_surf
        };
        let interfc = ch_fi3d_fil_curve_in_ds(ind_curve, inda_face, Some(pcint21), or_face);
        b.my_ds
            .change_ds()
            .change_shape_interferences_at(inda_face)
            .push(interfc);
    }

    // `aTrCracc->D0(First/LastParameter(), aPext1 / aPext2);` then
    // `cxx:4347-4358`.
    let pext1 = sec.trcracc.d0(sec.trcracc.first_parameter());
    let pext2 = sec.trcracc.d0(sec.trcracc.last_parameter());
    if pext1.distance(st.cp2_not_on_arc.point()) <= tol3d {
        ind_point1 = ind_cp2_not_on_arc;
        ind_point2 = ind_point;
    } else {
        ind_point1 = ind_point;
        ind_point2 = ind_cp2_not_on_arc;
    }

    // `cxx:4361-4383`: the orientation of `aTrCracc`.
    let is_to_reverse = if sec.is_do_second_section {
        let cint2 = sec.cint2.as_ref().expect("isDoSecondSection");
        let p1 = sec.trcracc.d0(sec.trcracc.first_parameter());
        let p2 = sec.trcracc.d0(sec.trcracc.last_parameter());
        let p3 = cint2.d0(cint2.first_parameter());
        let p4 = cint2.d0(cint2.last_parameter());
        p1.distance(&p4) > tol3d && p2.distance(&p3) > tol3d
    } else {
        pext1.distance(&point2) > tol3d && pext2.distance(&point1) > tol3d
    };
    if is_to_reverse {
        or_surf = or_surf.reversed();
    }

    // `TopOpeBRepDS_Curve tct2(aTrCracc, aTolReached); indCurve =
    //  DStr.AddCurve(tct2);` (`cxx:4386-4389`).
    let ind_curve = {
        let mut tct2 = DsCurve::with_curve(Some(Arc::clone(&sec.trcracc)), sec.tol_reached, false);
        b.my_ds.change_ds().add_curve(&mut tct2)
    };
    // `cxx:4390-4395`.
    let interfp1 = ch_fi3d_fil_point_in_ds(
        Orientation::Forward,
        ind_curve,
        ind_point1,
        sec.trcracc.first_parameter(),
        false,
    );
    let interfp2 = ch_fi3d_fil_point_in_ds(
        Orientation::Reversed,
        ind_curve,
        ind_point2,
        sec.trcracc.last_parameter(),
        false,
    );
    {
        let ds = b.my_ds.change_ds();
        let list = ds.change_curve_interferences(ind_curve);
        list.push(interfp1);
        list.push(interfp2);
    }
    // `anInterfc = ChFi3d_FilCurveInDS(indCurve, aSurfData->Surf(), aPCraccS,
    //  anOrSurf);` (`cxx:4398-4399`).
    let interfc =
        ch_fi3d_fil_curve_in_ds(ind_curve, fd_prev.surf(), pcracc_s, or_surf);
    b.my_ds
        .change_ds()
        .change_surface_interferences(fd_prev.surf())
        .push(interfc);
    // `aStripe->InDS(isFirst);` (`cxx:4399`).
    if let Ok(mut stripe) = st.stripe.write() {
        stripe.in_ds(st.is_first, 1);
    }
    // `if (anOrSD1 == anOrSD2) anOrSurf = TopAbs::Reverse(anOrSurf);`
    // (`cxx:4401-4404`).
    if sec.or_sd1 == sec.or_sd2 {
        or_surf = or_surf.reversed();
    }
    // `anInterfc = ChFi3d_FilCurveInDS(indCurve, indSurf1, aPCurv1, anOrSurf);`
    // (`cxx:4406-4408`).
    let interfc = ch_fi3d_fil_curve_in_ds(
        ind_curve,
        ind_surf1,
        sec.pcurv1.clone(),
        or_surf,
    );
    b.my_ds
        .change_ds()
        .change_surface_interferences(ind_surf1)
        .push(interfc);

    // `aCint1->D0(aCint1->FirstParameter(), aPext1);` then `cxx:4415-4422`.
    let pext1 = sec.cint1.d0(sec.cint1.first_parameter());
    if pext1.distance(st.cp1_on_arc.point()) <= tol3d {
        ind_point1 = ind_cp1_on_arc;
        ind_point2 = ind_point;
    } else {
        ind_point1 = ind_point;
        ind_point2 = ind_cp1_on_arc;
    }

    // `cxx:4426-4436`: the orientation of `aCint1` against `aTrCracc`.
    let p1 = sec.cint1.d0(sec.cint1.first_parameter());
    let p2 = sec.cint1.d0(sec.cint1.last_parameter());
    let p3 = sec.trcracc.d0(sec.trcracc.first_parameter());
    let p4 = sec.trcracc.d0(sec.trcracc.last_parameter());
    if p1.distance(&p4) > tol3d && p2.distance(&p3) > tol3d {
        or_surf = or_surf.reversed();
    }

    // `TopOpeBRepDS_Curve aTCint1(aCint1, aTolex1); indCurve =
    //  DStr.AddCurve(aTCint1);` (`cxx:4438-4440`).
    let ind_curve = {
        let mut tcint1 = DsCurve::with_curve(Some(Arc::clone(&sec.cint1)), sec.tolex1, false);
        b.my_ds.change_ds().add_curve(&mut tcint1)
    };
    let interfp1 = ch_fi3d_fil_point_in_ds(
        Orientation::Forward,
        ind_curve,
        ind_point1,
        sec.cint1.first_parameter(),
        false,
    );
    let interfp2 = ch_fi3d_fil_point_in_ds(
        Orientation::Reversed,
        ind_curve,
        ind_point2,
        sec.cint1.last_parameter(),
        false,
    );
    {
        let ds = b.my_ds.change_ds();
        let list = ds.change_curve_interferences(ind_curve);
        list.push(interfp1);
        list.push(interfp2);
    }
    // `anInterfc = ChFi3d_FilCurveInDS(indCurve, indSurf1, aPCint12, anOrSurf);`
    // (`cxx:4446-4448`).
    let interfc = ch_fi3d_fil_curve_in_ds(
        ind_curve,
        ind_surf1,
        Some(Arc::clone(&sec.pcint12)),
        or_surf,
    );
    b.my_ds
        .change_ds()
        .change_surface_interferences(ind_surf1)
        .push(interfc);

    // `anOrFace = aFace.Orientation(); if (anOrFace == anOrSD1) anOrFace =
    //  TopAbs::Reverse(anOrSurf); else anOrFace = anOrSurf;` (`cxx:4452-4461`).
    or_face = st.face.orientation();
    or_face = if or_face == sec.or_sd1 {
        or_surf.reversed()
    } else {
        or_surf
    };
    let interfc = ch_fi3d_fil_curve_in_ds(
        ind_curve,
        inda_face,
        Some(Arc::clone(&sec.pcint11)),
        or_face,
    );
    b.my_ds
        .change_ds()
        .change_shape_interferences_at(inda_face)
        .push(interfc);

    // `aTrCline = new Geom_TrimmedCurve(aCline, First, Last); aTolerance =
    //  DStr.Curve(indLine).Tolerance(); aTct3(aTrCline, aTolerance); indCurve =
    //  DStr.AddCurve(aTct3);` (`cxx:4466-4472`).
    let tolerance = b.my_ds.ds().curve(st.ind_line).tolerance();
    let tr_cline: Arc<dyn Curve> = Arc::new(
        crate::topopebrep::tool::adaptor::GeomAdaptorCurve::with_range(
            Arc::clone(&st.cline),
            st.cline.first_parameter(),
            st.cline.last_parameter(),
        ),
    );
    let ind_curve = {
        let mut tct3 = DsCurve::with_curve(Some(Arc::clone(&tr_cline)), tolerance, false);
        b.my_ds.change_ds().add_curve(&mut tct3)
    };

    // `aTrCline->D0(aTrCline->FirstParameter(), aPext1); if (aPext1.Distance(
    //  aCP1onArc.Point()) < aTol3d) { ... } else { ... }` (`cxx:4474-4483`).
    let pext1 = tr_cline.d0(tr_cline.first_parameter());
    if pext1.distance(st.cp1_on_arc.point()) < tol3d {
        ind_point1 = ind_cp1_on_arc;
        ind_point2 = ind_cp2_not_on_arc;
    } else {
        ind_point1 = ind_cp2_not_on_arc;
        ind_point2 = ind_cp1_on_arc;
    }
    // `cxx:4486-4497`: the orientation of `aTrCline` against `aCint1`.
    let p1 = tr_cline.d0(tr_cline.first_parameter());
    let p2 = tr_cline.d0(tr_cline.last_parameter());
    let p3 = sec.cint1.d0(sec.cint1.first_parameter());
    let p4 = sec.cint1.d0(sec.cint1.last_parameter());
    if p1.distance(&p4) > tol3d && p2.distance(&p3) > tol3d {
        or_surf = or_surf.reversed();
    }
    let interfp1 = ch_fi3d_fil_point_in_ds(
        Orientation::Forward,
        ind_curve,
        ind_point1,
        tr_cline.first_parameter(),
        false,
    );
    let interfp2 = ch_fi3d_fil_point_in_ds(
        Orientation::Reversed,
        ind_curve,
        ind_point2,
        tr_cline.last_parameter(),
        false,
    );
    {
        let ds = b.my_ds.change_ds();
        let list = ds.change_curve_interferences(ind_curve);
        list.push(interfp1);
        list.push(interfp2);
    }
    // `anInterfc = ChFi3d_FilCurveInDS(indCurve, indSurf1, aPClineOnSurf,
    //  anOrSurf);` (`cxx:4504-4506`).
    let interfc = ch_fi3d_fil_curve_in_ds(
        ind_curve,
        ind_surf1,
        st.pc_line_on_surf.clone(),
        or_surf,
    );
    b.my_ds
        .change_ds()
        .change_surface_interferences(ind_surf1)
        .push(interfc);

    // `indShape = DStr.AddShape(aNeighborFace); anOrFace =
    //  aNeighborFace.Orientation(); if (anOrFace == anOrSD1) anOrFace =
    //  TopAbs::Reverse(anOrSurf); else anOrFace = anOrSurf;` (`cxx:4509-4519`).
    let ind_shape = b.my_ds.change_ds().add_shape(&st.neighbor_face.0);
    or_face = st.neighbor_face.orientation();
    or_face = if or_face == sec.or_sd1 {
        or_surf.reversed()
    } else {
        or_surf
    };
    let interfc = ch_fi3d_fil_curve_in_ds(ind_curve, ind_shape, st.pc_line_on_face.clone(), or_face);
    b.my_ds
        .change_ds()
        .change_shape_interferences_at(ind_shape)
        .push(interfc);

    Ok(())
}
