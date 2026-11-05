//! `PerformIntersectionAtEnd` prologue — `ChFi3d_Builder_C1.cxx:1828-2136`.
//!
//! Face / edge counts (`ChFi3d_nbface`, `ChFi3d_NbNotDegeneratedEdges`,
//! `ChFi3d_ChercheBordsLibres`), the sewing edge of the vertex
//! (`ChFi3d_CoutureOnVertex`), the regular-adjacent-edge analysis
//! (`ChFi3d_edge_common_faces` + `ChFi3d::IsTangentFaces`), the
//! `PerformMoreSurfdata` / `PerformMoreThreeCorner` dispatches
//! (`cxx:1998-2010`), the bounded fillet surface, the orientation of the end
//! curve and the `checkShrink` detection (`cxx:2088-2128`).

use std::sync::Arc;

use occt_core::gp::GpPnt2d;
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::{CONFUSION, PCONFUSION};

use crate::abs::Orientation;
use crate::chfi_ds::{FaceInterference, State, SurfData};
use crate::chfi3d::builder_0_helpers::{
    ch_fi3d_bound_surf, ch_fi3d_index_of_surf_data, ch_fi3d_index_point_in_ds,
};
use crate::chfi3d::builder_helpers::{
    cherche_bords_libres, cherche_element, couture_on_vertex, edge_common_faces,
    nb_not_degenerated_edges, nbface,
};
use crate::chfi3d::tangent_faces::is_tangent_faces;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools::is_same;

use super::super::ChFi3dBuilder;
use super::{EndIntersection, EndIntersectionExit};

/// What the prologue answers: the locals for the rest of the body, or the
/// corner hand-over that already returned in OCCT.
pub(super) enum PrologueOut {
    Continue(Box<EndIntersection>),
    Exit(EndIntersectionExit),
}

/// `SeqFil(i)`, 1-based like the OCCT sequence; out of range is the
/// out-of-bounds read OCCT would take.
fn seq_at(seq: &[Arc<SurfData>], i: i64) -> Result<&Arc<SurfData>, String> {
    usize::try_from(i)
        .ok()
        .and_then(|k| k.checked_sub(1))
        .and_then(|k| seq.get(k))
        .ok_or_else(|| "PerformIntersectionAtEnd: SeqFil index out of range".to_string())
}

pub(super) fn prologue(
    b: &mut ChFi3dBuilder,
    index: usize,
) -> Result<PrologueOut, String> {
    // `const TopoDS_Vertex& Vtx = myVDataMap.FindKey(Index);` and the iterator's
    // first value (`cxx:1845-1849`).
    let Some(vtx) = b.my_v_data_map.find_key(index).cloned() else {
        return Err("PerformIntersectionAtEnd: no vertex at Index".to_string());
    };
    let Some(stripe) = b
        .my_v_data_map
        .find_from_index(index)
        .and_then(|l| l.first())
        .cloned()
    else {
        return Err("PerformIntersectionAtEnd: empty stripe list at Index".to_string());
    };
    let vertex = Vertex::wrap(vtx.clone())
        .ok_or_else(|| "PerformIntersectionAtEnd: myVDataMap key is not a vertex".to_string())?;
    let spine = {
        let guard = stripe
            .read()
            .map_err(|_| "PerformIntersectionAtEnd: poisoned stripe".to_string())?;
        guard
            .spine()
            .cloned()
            .ok_or_else(|| "PerformIntersectionAtEnd: stripe has no spine".to_string())?
    };
    let mut seq: Vec<Arc<SurfData>> = {
        let guard = stripe
            .read()
            .map_err(|_| "PerformIntersectionAtEnd: poisoned stripe".to_string())?;
        guard.set_of_surf_data().to_vec()
    };

    // `int nface = ChFi3d_nbface(myVFMap(Vtx));` (`cxx:1858`).
    let nface = nbface(b.my_vf_map.get_from_key(&vtx));

    // `nbarete = ChFi3d_NbNotDegeneratedEdges(Vtx, myVEMap);
    //  ChFi3d_ChercheBordsLibres(myVEMap, Vtx, bordlibre, edgelibre1,
    //  edgelibre2);` (`cxx:1861-1862`).
    let mut nbarete = nb_not_degenerated_edges(&vertex, &b.my_ve_map);
    let (bord_libre, edge_libre1, edge_libre2) = cherche_bords_libres(&b.my_ve_map, &vertex);
    // `cxx:1863-1870`.
    if bord_libre {
        nbarete = (nbarete - 2) / 2 + 2;
    } else {
        nbarete /= 2;
    }

    // `for (ItF.Initialize(myVFMap(Vtx)); ItF.More() && !couture; ItF.Next())`
    // (`cxx:1877-1885`).
    let mut couture = false;
    let mut face_couture: Option<Face> = None;
    let mut edge_couture: Option<Edge> = None;
    for value in b.my_vf_map.get_from_key(&vtx) {
        if couture {
            break;
        }
        let Some(fcur) = Face::wrap(value.clone()) else {
            continue;
        };
        // `ChFi3d_CoutureOnVertex(fcur, Vtx, couture, edgecouture);`
        if let Some(e) = couture_on_vertex(&fcur, &vertex) {
            couture = true;
            edge_couture = Some(e);
            face_couture = Some(fcur);
        }
    }

    // `nbsurf = SeqFil.Length(); nbedge = spine->NbEdges();` (`cxx:1894-1895`).
    let nbsurf = seq.len();
    let nbedge = {
        let guard = spine
            .read()
            .map_err(|_| "PerformIntersectionAtEnd: poisoned spine".to_string())?;
        guard.nb_edges()
    };
    // `num = ChFi3d_IndexOfSurfData(Vtx, stripe, sens); isfirst = (sens == 1);`
    // (`cxx:1896-1897`).
    let mut sens = 0;
    let num = ch_fi3d_index_of_surf_data(&vtx, &stripe, &mut sens)?;
    let num = usize::try_from(num)
        .map_err(|_| "PerformIntersectionAtEnd: IndexOfSurfData out of range".to_string())?;
    let is_first = sens == 1;

    // `cxx:1898-1910`.
    let (edge_spine, mut num1, state) = {
        let guard = spine
            .read()
            .map_err(|_| "PerformIntersectionAtEnd: poisoned spine".to_string())?;
        if is_first {
            let e = guard
                .edges(1)
                .cloned()
                .ok_or_else(|| "PerformIntersectionAtEnd: spine has no first edge".to_string())?;
            (e, num as i64 + 1, guard.first_status())
        } else {
            let e = guard
                .edges(nbedge)
                .cloned()
                .ok_or_else(|| "PerformIntersectionAtEnd: spine has no last edge".to_string())?;
            (e, num as i64 - 1, guard.last_status())
        }
    };

    // `F1` / `F2` are only assigned inside the `nbsurf != nbedge && nbsurf != 1`
    // block (`cxx:1913`); the later stages re-read them from the DS
    // (`cxx:2146-2147`).
    let mut f1: Option<Face> = None;
    let mut f2: Option<Face> = None;
    let mut eadj1: Option<Edge> = None;
    let mut eadj2: Option<Edge> = None;
    let mut reg1 = false;
    let mut reg2 = false;

    // `if (nbsurf != nbedge && nbsurf != 1)` (`cxx:1911-1979`).
    if nbsurf != nbedge && nbsurf != 1 {
        // `ChFi3d_edge_common_faces(myEFMap(EdgeSpine), F1, F2);`
        let (common1, common2) = edge_common_faces(b.my_ef_map.get_from_key(&edge_spine))
            .ok_or_else(|| {
                "PerformIntersectionAtEnd: EdgeSpine has no face in myEFMap".to_string()
            })?;
        f1 = Some(common1.clone());
        f2 = Some(common2.clone());

        // `if (F1.IsSame(facecouture)) Eadj1 = edgecouture; else
        //  ChFi3d_cherche_element(Vtx, EdgeSpine, F1, Eadj1, Vbid1);`
        // (`cxx:1914-1921`).
        let e1 = match &face_couture {
            Some(fc) if is_same(&common1.0, &fc.0) => edge_couture.clone(),
            _ => None,
        };
        let e1 = match e1 {
            Some(e) => e,
            None => cherche_element(&vertex, &Edge(edge_spine.clone()), &common1)?.0,
        };
        let (fga, fdr) =
            edge_common_faces(b.my_ef_map.get_from_key(&e1.0)).ok_or_else(|| {
                "PerformIntersectionAtEnd: Eadj1 has no face in myEFMap".to_string()
            })?;
        // `reg1 = ChFi3d::IsTangentFaces(Eadj1, Fga, Fdr);` (`cxx:1925`).
        reg1 = is_tangent_faces(&e1, &fga, &fdr, Shape::G1)?;
        eadj1 = Some(e1);

        // The same for `F2` (`cxx:1927-1939`).
        let e2 = match &face_couture {
            Some(fc) if is_same(&common2.0, &fc.0) => edge_couture.clone(),
            _ => None,
        };
        let e2 = match e2 {
            Some(e) => e,
            None => cherche_element(&vertex, &Edge(edge_spine.clone()), &common2)?.0,
        };
        let (fga, fdr) =
            edge_common_faces(b.my_ef_map.get_from_key(&e2.0)).ok_or_else(|| {
                "PerformIntersectionAtEnd: Eadj2 has no face in myEFMap".to_string()
            })?;
        reg2 = is_tangent_faces(&e2, &fga, &fdr, Shape::G1)?;
        eadj2 = Some(e2);

        // `if (reg1 || reg2)` (`cxx:1942-1978`).
        if reg1 || reg2 {
            // `cp1 = SeqFil(num1)->ChangeVertex(isfirst, 1);` (copied, the
            // source only reads it) and `cp2` likewise (`cxx:1947-1948`).
            let cp1 = seq_at(&seq, num1)?.vertex(is_first, 1).clone();
            let cp2 = seq_at(&seq, num1)?.vertex(is_first, 2).clone();
            // `cxx:1949-1962`: `cp.Arc().IsSame(Eadj1) || cp.Arc().IsSame(Eadj2)`.
            let is_adjacent = |a: &TopoShape| {
                eadj1.as_ref().is_some_and(|e| is_same(a, &e.0))
                    || eadj2.as_ref().is_some_and(|e| is_same(a, &e.0))
            };
            let compoint1 = cp1.is_on_arc() && cp1.arc().is_some_and(is_adjacent);
            let compoint2 = cp2.is_on_arc() && cp2.arc().is_some_and(is_adjacent);
            // `if (compoint1 && compoint2) { SeqFil.Remove(num);
            //  num = ChFi3d_IndexOfSurfData(Vtx, stripe, sens); ... }`
            // (`cxx:1963-1977`). The removal is applied to the real sequence
            // because `ChFi3d_IndexOfSurfData` re-reads it.
            if compoint1 && compoint2 {
                if let Ok(mut guard) = stripe.write() {
                    guard.change_set_of_surf_data().remove(num - 1);
                }
                seq = {
                    let guard = stripe
                        .read()
                        .map_err(|_| "PerformIntersectionAtEnd: poisoned stripe".to_string())?;
                    guard.set_of_surf_data().to_vec()
                };
                let mut sens2 = 0;
                let num_new = ch_fi3d_index_of_surf_data(&vtx, &stripe, &mut sens2)?;
                let num_new = usize::try_from(num_new).map_err(|_| {
                    "PerformIntersectionAtEnd: IndexOfSurfData out of range".to_string()
                })?;
                if is_first {
                    num1 = num_new as i64 + 1;
                } else {
                    num1 = num_new as i64 - 1;
                }
                reg1 = false;
                reg2 = false;
            }
        }
    }

    // `Fd = SeqFil.ChangeValue(num); CV1 = Fd->ChangeVertex(isfirst, 1);`
    // (`cxx:1983-1985`). The source mutates `Fd` later; the port keeps the
    // stripe's sequence authoritative and copies the entry per stage.
    let num_v = num as i64;
    let fd = seq_at(&seq, num_v)?.clone();
    let cv1 = fd.vertex(is_first, 1).clone();
    let cv2 = fd.vertex(is_first, 2).clone();

    // `if (FindFace(Vtx, CV1, CV2, face))` (`cxx:1987-1997`).
    let face = b.find_face(&vtx, &cv1, &cv2, None);
    let mut one_corner = false;
    if let Some(f) = &face {
        if !couture {
            one_corner = true;
        } else if !face_couture.as_ref().is_some_and(|fc| is_same(&f.0, &fc.0)) {
            one_corner = true;
        }
    }

    // `if (onecorner) { if (MoreSurfdata(Index)) { PerformMoreSurfdata(Index);
    //  return; } }` (`cxx:1998-2005`).
    if one_corner && b.more_surfdata(index) {
        return Ok(PrologueOut::Exit(EndIntersectionExit::MoreSurfdata));
    }
    // `if (!onecorner && (reg1 || reg2) && !couture && state != ChFiDS_OnSame)
    //  { PerformMoreThreeCorner(Index, 1); return; }` (`cxx:2006-2010`).
    if !one_corner && (reg1 || reg2) && !couture && state != State::OnSame {
        return Ok(PrologueOut::Exit(EndIntersectionExit::MoreThreeCorner(1)));
    }

    // `HGs = ChFi3d_BoundSurf(DStr, Fd, 1, 2); Fi1 = Fd->InterferenceOnS1();
    //  Fi2 = Fd->InterferenceOnS2();` (`cxx:2011-2013`).
    let hgs = ch_fi3d_bound_surf(b.my_ds.change_ds(), &fd, 1, 2)?;
    let fi1: FaceInterference = fd.interference(1).clone();
    let fi2: FaceInterference = fd.interference(2).clone();

    // `orsurfdata = Fd->Orientation(); Isurf = Fd->Surf();` and `SDprev`
    // (`cxx:2021-2028`).
    let or_surf_data = fd.orientation();
    let isurf = fd.surf();
    let mut isurf_prev = 0;
    let mut sdprev: Option<Arc<SurfData>> = None;
    if num1 > 0 && num1 <= seq.len() as i64 {
        let sd = seq_at(&seq, num1)?.clone();
        isurf_prev = sd.surf();
        sdprev = Some(sd);
    }

    // `double tolpt = 1.e-4; double tolreached;` (`cxx:2031-2032`): OCCT leaves
    // `tolreached` uninitialised until `ChFi3d_ComputeCurves` writes it in the
    // per-face loop (`cxx:2597-3127`, ported in `face_loop`, which stores it
    // back through `EndIntersection::tol_reached`); the prologue starts it at 0.
    let tolpt = 1.0e-4;
    let tol_reached = 0.0;

    // `stripe->SetIndexPoint(ChFi3d_IndexPointInDS(CV1, DStr), isfirst, 1);`
    // (`cxx:2035-2036`).
    let ip1 = ch_fi3d_index_point_in_ds(&cv1, b.my_ds.change_ds());
    let ip2 = ch_fi3d_index_point_in_ds(&cv2, b.my_ds.change_ds());
    {
        let mut guard = stripe
            .write()
            .map_err(|_| "PerformIntersectionAtEnd: poisoned stripe".to_string())?;
        guard.set_index_point(ip1, is_first, 1);
        guard.set_index_point(ip2, is_first, 2);
    }

    // `Ishape1 = Fd->IndexOfS1(); if (Ishape1 != 0) { if (Ishape1 > 0) trafil1 =
    //  DStr.Shape(Ishape1).Orientation(); trafil1 = TopAbs::Compose(trafil1,
    //  Fd->Orientation()); trafil1 =
    //  TopAbs::Compose(TopAbs::Reverse(Fi1.Transition()), trafil1); }`
    // (`cxx:2041-2058`).
    let ishape1 = fd.index_of_s1();
    let mut trafil1 = Orientation::Forward;
    if ishape1 != 0 {
        if ishape1 > 0 {
            let shape = b.my_ds.change_ds().shape(ishape1, false).ok_or_else(|| {
                "PerformIntersectionAtEnd: DStr.Shape(IndexOfS1) is null".to_string()
            })?;
            trafil1 = shape.orientation();
        }
        trafil1 = Orientation::compose(trafil1, fd.orientation());
        trafil1 = Orientation::compose(fi1.transition().reversed(), trafil1);
    }
    // `if (!isfirst) orcourbe = trafil1; else orcourbe =
    //  TopAbs::Reverse(trafil1);` (`cxx:2071-2078`).
    let or_courbe = if !is_first {
        trafil1
    } else {
        trafil1.reversed()
    };

    // `P2d1 = Fi1.PCurveOnSurf()->Value(Fi1.Parameter(isfirst)); ...`
    // (`cxx:2096-2110`).
    let pc1 = fi1
        .p_curve_on_surf()
        .ok_or_else(|| "PerformIntersectionAtEnd: Fi1 has no p-curve on surf".to_string())?;
    let p2d1: GpPnt2d = pc1.value(fi1.parameter(is_first));
    let p2d2: GpPnt2d = pc1.value(fi1.parameter(!is_first));
    let a1 = hgs
        .d0(p2d1.x(), p2d1.y())
        .ok_or_else(|| "PerformIntersectionAtEnd: HGs evaluation failed".to_string())?;
    let a2 = hgs
        .d0(p2d2.x(), p2d2.y())
        .ok_or_else(|| "PerformIntersectionAtEnd: HGs evaluation failed".to_string())?;
    let fi1_length = a1.distance(&a2);
    let mut check_shrink = fi1_length <= CONFUSION;

    let pc2 = fi2
        .p_curve_on_surf()
        .ok_or_else(|| "PerformIntersectionAtEnd: Fi2 has no p-curve on surf".to_string())?;
    let p2d3: GpPnt2d = pc2.value(fi2.parameter(is_first));
    let p2d4: GpPnt2d = pc2.value(fi2.parameter(!is_first));
    let a1 = hgs
        .d0(p2d3.x(), p2d3.y())
        .ok_or_else(|| "PerformIntersectionAtEnd: HGs evaluation failed".to_string())?;
    let a2 = hgs
        .d0(p2d4.x(), p2d4.y())
        .ok_or_else(|| "PerformIntersectionAtEnd: HGs evaluation failed".to_string())?;
    let fi2_length = a1.distance(&a2);
    check_shrink = check_shrink || fi2_length <= CONFUSION;

    // `cxx:2112-2128`.
    let mut is_u_shrink = false;
    let mut check_shr_param = 0.0;
    if check_shrink {
        if (p2d2.y() - p2d4.y()).abs() <= PCONFUSION {
            is_u_shrink = false;
            check_shr_param = p2d2.y();
        } else if (p2d2.x() - p2d4.x()).abs() <= PCONFUSION {
            is_u_shrink = true;
            check_shr_param = p2d2.x();
        } else {
            check_shrink = false;
        }
    }

    Ok(PrologueOut::Continue(Box::new(EndIntersection {
        vtx,
        stripe,
        spine,
        seq,
        nface,
        nbarete,
        bord_libre,
        edge_libre1,
        edge_libre2,
        couture,
        face_couture,
        edge_couture,
        nbsurf,
        nbedge,
        num,
        num1,
        is_first,
        state,
        f1,
        f2,
        eadj1,
        eadj2,
        reg1,
        reg2,
        face,
        one_corner,
        hgs,
        fi1,
        fi2,
        cc: None,
        pc: None,
        ps: None,
        or_surf_data,
        isurf,
        isurf_prev,
        sdprev,
        tolpt,
        tol_reached,
        or_courbe,
        ishape1,
        trafil1,
        is_on_same1: false,
        is_on_same2: false,
        check_shrink,
        is_shrink: false,
        is_u_shrink,
        check_shr_param,
        prev_sd_param: 0.0,
        mid_p2d: GpPnt2d::zero(),
        mid_ipoint: 0,
        nbface: 0,
        nb: 1,
        edges: vec![None; super::NN],
        faces: vec![None; super::NN],
        tabedg: Vec::new(),
        find_on_f1: false,
        find_on_f2: false,
        one_intersection1: false,
        one_intersection2: false,
        face_sau: None,
        edge_sau: None,
        edge_spine,
        fi1_length,
        fi2_length,
    })))
}
