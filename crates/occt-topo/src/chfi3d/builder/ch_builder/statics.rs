//! The `ChFi3d_ChBuilder.cxx` file static `SearchCommonFaces` / the static
//! `ExtentSpineOnCommonFace` and the private member `ConexFaces`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_ChBuilder.cxx`
//! (OCCT 8.0.0): `SearchCommonFaces` (`cxx:67-93`), `ExtentSpineOnCommonFace`
//! (`cxx:105-186`), `ChFi3d_ChBuilder::ConexFaces` (`cxx:2283-2318`).
//!
//! Neither of the three has a fillet counterpart with the same body:
//! `SearchCommonFaces` is used by four `ChFi3d_ChBuilder` members only
//! (`SetDist`, `SetDists`, `SetDistAngle`, `ConexFaces`), `ExtentSpineOnCommonFace`
//! is called from `ExtentTwoCorner` / `ExtentThreeCorner` and `ConexFaces` from
//! the two `Extent*Corner` bodies — none of them exists in
//! `ChFi3d_FilBuilder.cxx`.
//!
//! Translation notes:
//! - `Spine->D1(Spine->Absc(V), tmp, tg1)` returns an undefined `gp_Pnt` /
//!   `gp_Vec` when `Absc` answers its `-1.` sentinel for a vertex that is not on
//!   the guideline (`ChFiDS_Spine.cxx:304-330`); the port answers `Err` there
//!   instead of extending the spine with a garbage tangent,
//! - `Spine->SetFirstParameter` / `SetFirstTgt` raise on a periodic guideline
//!   (`ChFiDS_Spine.cxx:353-392`); the port answers `Err`,
//! - `TopoDS::Face(It.Value())` of a shape of another type raises in OCCT; the
//!   port reports `Err`,
//! - a null `TopoDS_Face` (`F1.Nullify()`) is `None`; every caller of
//!   `SearchCommonFaces` hands both faces to something that dereferences them
//!   (`BRepAdaptor_Surface::Initialize`), so the callers report `Err` on `None`.

use crate::abs::Orientation;
use crate::brep_tools::is_really_closed;
use crate::chfi3d::side::concave_side;
use crate::chfi_ds::{ChFiDSMap, SpineHandle};
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::is_same;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `SearchCommonFaces(EFMap, E, F1, F2)` (`cxx:67-93`): the first two distinct
/// faces of the `EFMap` list of `e`.
///
/// The last block (`cxx:90-93`) answers `F2 = F1` when `E` is a closed edge on
/// `F1` and no second face was found, so a seam edge yields the same face
/// twice.
pub(super) fn search_common_faces(
    ef_map: &ChFiDSMap,
    e: &TopoShape,
) -> Result<(Option<Face>, Option<Face>), String> {
    // `TopoDS_Face Fc; F1.Nullify(); F2.Nullify();`
    // `for (It.Initialize(EFMap(E)); It.More(); It.Next())` (`cxx:72-89`).
    let mut f1: Option<Face> = None;
    let mut f2: Option<Face> = None;
    for shape in ef_map.get_from_key(e) {
        // `Fc = TopoDS::Face(It.Value());` (`cxx:79`).
        let fc = Face::wrap(shape.clone())
            .ok_or_else(|| "SearchCommonFaces: an ancestor of the edge is not a face".to_string())?;
        // `if (F1.IsNull()) F1 = Fc; else if (!Fc.IsSame(F1)) { F2 = Fc;
        // break; }` (`cxx:80-88`).
        match &f1 {
            None => f1 = Some(fc),
            Some(first) => {
                if !is_same(&fc.0, &first.0) {
                    f2 = Some(fc);
                    break;
                }
            }
        }
    }

    // `if (!F1.IsNull() && F2.IsNull() && BRepTools::IsReallyClosed(E, F1))
    // F2 = F1;` (`cxx:90-93`).
    if f2.is_none() {
        if let Some(first) = &f1 {
            let edge = Edge::wrap(e.clone()).ok_or_else(|| {
                "SearchCommonFaces: ChFiDS_Map key of an edge is not an edge".to_string()
            })?;
            if is_really_closed(&edge, first) {
                f2 = Some(first.clone());
            }
        }
    }
    Ok((f1, f2))
}

/// `ExtentSpineOnCommonFace(Spine1, Spine2, V, dis1, dis2, isfirst1, isfirst2)`
/// (`cxx:105-186`): extend the two guidelines of two chamfers by their common
/// face distance.
///
/// `isfirst(i)` is `false` when `Spine(i)` is oriented towards `V`. The three
/// guards of the body are kept literally: `tolesp = 1.e-7` (`cxx:112`), the
/// amplification `*3.` for `Spine1` (`cxx:158`) against `*1.5` for `Spine2`
/// (`cxx:174`), and the "extend only when the computed distance is positive"
/// test (`cxx:154`, `:171`).
pub(super) fn extent_spine_on_common_face(
    spine1: &SpineHandle,
    spine2: &SpineHandle,
    v: &TopoShape,
    dis1: f64,
    dis2: f64,
    isfirst1: bool,
    isfirst2: bool,
) -> Result<(), String> {
    // `double tolesp = 1.e-7;` (`cxx:112`).
    const TOLESP: f64 = 1.0e-7;

    // `Spine1->D1(Spine1->Absc(V), tmp, tg1); Spine2->D1(Spine2->Absc(V), tmp,
    // tg2); tg1.Normalize(); tg2.Normalize();` (`cxx:122-126`).
    let (mut tg1, mut tg2) = {
        let mut g1 = spine1
            .write()
            .map_err(|_| "ExtentSpineOnCommonFace: poisoned first spine".to_string())?;
        let a1 = g1.absc_of_vertex(v);
        let tg1 = g1
            .d1(a1)
            .ok_or_else(|| {
                "ExtentSpineOnCommonFace: Spine1->D1(Spine1->Absc(V))".to_string()
            })?
            .1;
        drop(g1);
        let mut g2 = spine2
            .write()
            .map_err(|_| "ExtentSpineOnCommonFace: poisoned second spine".to_string())?;
        let a2 = g2.absc_of_vertex(v);
        let tg2 = g2
            .d1(a2)
            .ok_or_else(|| {
                "ExtentSpineOnCommonFace: Spine2->D1(Spine2->Absc(V))".to_string()
            })?
            .1;
        (tg1, tg2)
    };
    tg1.normalize();
    tg2.normalize();
    // `if (isfirst1) tg1.Reverse(); if (isfirst2) tg2.Reverse();`
    // (`cxx:127-133`).
    if isfirst1 {
        tg1.reverse();
    }
    if isfirst2 {
        tg2.reverse();
    }

    // `cosalpha = tg1.Dot(tg2); sinalpha = sqrt(1 - cosalpha * cosalpha);`
    // (`cxx:137-138`).
    let cosalpha = tg1.dot(&tg2);
    let sinalpha = (1.0 - cosalpha * cosalpha).sqrt();

    // `temp = cosalpha + dis2 / dis1; if (|temp| > tolesp) { tga1 = sinalpha /
    // temp; d1plus = dis1 / tga1; }` (`cxx:141-146`).
    let mut d1plus = 0.0f64;
    let temp = cosalpha + dis2 / dis1;
    if temp.abs() > TOLESP {
        let tga1 = sinalpha / temp;
        d1plus = dis1 / tga1;
    }
    // `temp = cosalpha + dis1 / dis2; ...` (`cxx:147-152`).
    let mut d2plus = 0.0f64;
    let temp = cosalpha + dis1 / dis2;
    if temp.abs() > TOLESP {
        let tga2 = sinalpha / temp;
        d2plus = dis2 / tga2;
    }

    // `if (d1plus > 0.) { d1plus *= 3.; ... }` (`cxx:154-169`).
    if d1plus > 0.0 {
        d1plus *= 3.0;
        let mut guard = spine1
            .write()
            .map_err(|_| "ExtentSpineOnCommonFace: poisoned first spine".to_string())?;
        if isfirst1 {
            // `Spine1->SetFirstParameter(-d1plus); Spine1->SetFirstTgt(0.);`
            guard.set_first_parameter(-d1plus);
            guard
                .set_first_tgt(0.0)
                .ok_or_else(|| "ExtentSpineOnCommonFace: SetFirstTgt on a periodic spine".to_string())?;
        } else {
            // `param = Spine1->LastParameter(Spine1->NbEdges());
            // Spine1->SetLastParameter(d1plus + param);
            // Spine1->SetLastTgt(param);`
            let param = guard
                .last_parameter_of(guard.nb_edges() as i32)
                .ok_or_else(|| {
                    "ExtentSpineOnCommonFace: ChFiDS_Spine::LastParameter(NbEdges)".to_string()
                })?;
            guard.set_last_parameter(d1plus + param);
            guard.set_last_tgt(param).ok_or_else(|| {
                "ExtentSpineOnCommonFace: SetLastTgt on a periodic spine".to_string()
            })?;
        }
    }
    // `if (d2plus > 0.) { d2plus *= 1.5; ... }` (`cxx:170-185`).
    if d2plus > 0.0 {
        d2plus *= 1.5;
        let mut guard = spine2
            .write()
            .map_err(|_| "ExtentSpineOnCommonFace: poisoned second spine".to_string())?;
        if isfirst2 {
            guard.set_first_parameter(-d2plus);
            guard.set_first_tgt(0.0).ok_or_else(|| {
                "ExtentSpineOnCommonFace: SetFirstTgt on a periodic spine".to_string()
            })?;
        } else {
            let param = guard
                .last_parameter_of(guard.nb_edges() as i32)
                .ok_or_else(|| {
                    "ExtentSpineOnCommonFace: ChFiDS_Spine::LastParameter(NbEdges)".to_string()
                })?;
            guard.set_last_parameter(d2plus + param);
            guard.set_last_tgt(param).ok_or_else(|| {
                "ExtentSpineOnCommonFace: SetLastTgt on a periodic spine".to_string()
            })?;
        }
    }
    Ok(())
}

/// `ChFi3d_ChBuilder::ConexFaces(Spine, IEdge, F1, F2)` (`cxx:2283-2318`): the
/// two faces of `Spine->Edges(IEdge)` ordered so that `F1` carries the first
/// chamfer distance.
///
/// `F1` / `F2` are compared with the reference orientation of `Edges(1)`
/// (`cxx:2294-2301`) and swapped when the two `ConcaveSide` answers differ in
/// parity (`cxx:2309-2317`).
pub(super) fn cone_x_faces(
    ef_map: &ChFiDSMap,
    spine: &SpineHandle,
    i_edge: usize,
) -> Result<(Face, Face), String> {
    let (first_edge, iedge_edge) = {
        let guard = spine
            .read()
            .map_err(|_| "ConexFaces: poisoned spine".to_string())?;
        let first = guard.edges(1).cloned().ok_or_else(|| {
            "ConexFaces: ChFiDS_Spine::Edges(1) out of range".to_string()
        })?;
        let iedge = guard
            .edges(i_edge)
            .cloned()
            .ok_or_else(|| format!("ConexFaces: ChFiDS_Spine::Edges({i_edge}) out of range"))?;
        (first, iedge)
    };

    // `SearchCommonFaces(myEFMap, Spine->Edges(1), ff1, ff2);
    // ff1.Orientation(TopAbs_FORWARD); Sb1.Initialize(ff1);
    // ff2.Orientation(TopAbs_FORWARD); Sb2.Initialize(ff2);`
    // (`cxx:2294-2299`).
    let (ff1, ff2) = search_common_faces(ef_map, &first_edge)?;
    let (Some(ff1), Some(ff2)) = (ff1, ff2) else {
        return Err("ConexFaces: the first spine edge has no two faces".to_string());
    };
    let reff1 = Face(ff1.0.oriented(Orientation::Forward));
    let reff2 = Face(ff2.0.oriented(Orientation::Forward));
    let sb1 = SurfaceAdaptor::of_face(&reff1);
    let sb2 = SurfaceAdaptor::of_face(&reff2);
    // `RC = ChFi3d::ConcaveSide(Sb1, Sb2, Spine->Edges(1), tmp1, tmp2);`
    // (`cxx:2300`).
    let mut tmp1 = Orientation::Forward;
    let mut tmp2 = Orientation::Forward;
    let rc = concave_side(&sb1, &sb2, &Edge(first_edge), &mut tmp1, &mut tmp2)?;

    // `SearchCommonFaces(myEFMap, Spine->Edges(IEdge), f1, f2);
    // Sb1.Initialize(f1); Sb2.Initialize(f2); Choix =
    // ChFi3d::ConcaveSide(Sb1, Sb2, Spine->Edges(IEdge), tmp1, tmp2);`
    // (`cxx:2303-2307`) — no `Orientation(TopAbs_FORWARD)` here, unlike the
    // reference pair above.
    let (f1, f2) = search_common_faces(ef_map, &iedge_edge)?;
    let (Some(f1), Some(f2)) = (f1, f2) else {
        return Err("ConexFaces: the spine edge has no two faces".to_string());
    };
    let sb1 = SurfaceAdaptor::of_face(&f1);
    let sb2 = SurfaceAdaptor::of_face(&f2);
    let choix = concave_side(&sb1, &sb2, &Edge(iedge_edge), &mut tmp1, &mut tmp2)?;

    // `if (RC % 2 != Choix % 2) { F1 = f2; F2 = f1; } else { F1 = f1;
    // F2 = f2; }` (`cxx:2309-2317`).
    if rc % 2 != choix % 2 {
        Ok((f2, f1))
    } else {
        Ok((f1, f2))
    }
}
