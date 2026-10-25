//! Face-pair queries used to classify a corner: `ChFi3d_conexfaces`,
//! `ChFi3d_EdgeState` and `ChFi3d_evalconti`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:473-599`.

use occt_core::kernel::geomabs::Shape;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::chfi3d::side::concave_side;
use crate::chfi_ds::enums::State;
use crate::chfi_ds::map::ChFiDSMap;
use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::shape::{Edge, Face};
use crate::topo_tools::is_same;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `ChFi3d_conexfaces(E, F1, F2, EFMap)` (`cxx:473-499`): the two faces of
/// `E` from the edge→faces map.
///
/// The map holds one entry per edge *occurrence*, so a seam edge lists its face
/// twice; the second slot is then discarded and refilled while
/// `F2.IsSame(F1) && !BRep_Tool::IsClosed(E, F1)` (`cxx:492-498`).
///
/// `F1` is `None` when the map has no entry for `E` (OCCT leaves the output
/// null), `F2` is `None` for an open boundary or for a seam whose map entry
/// was consumed by the loop above.
pub fn conex_faces(e: &Edge, ef_map: &ChFiDSMap) -> (Option<Face>, Option<Face>) {
    let mut f1: Option<Face> = None;
    let mut f2: Option<Face> = None;
    for value in ef_map.get_from_key(&e.0) {
        let Some(face) = Face::wrap(value.clone()) else {
            // `TopoDS::Face` would raise `Standard_TypeMismatch`; the EF map is
            // filled from `TopAbs_FACE` ancestors only, so this cannot happen
            // for a map built by `ChFiDS_Map::Fill(E, FACE)`.
            continue;
        };
        match &f1 {
            None => f1 = Some(face),
            Some(first) => {
                f2 = Some(face);
                let second = f2.as_ref().expect("just assigned");
                if !is_same(&second.0, &first.0) || BRepTool::is_closed_edge_face(e, first) {
                    break;
                } else {
                    f2 = None;
                }
            }
        }
    }
    (f1, f2)
}

/// `ChFi3d_EdgeState(E[3], EFMap)` (`cxx:506-580`): the corner state of three
/// edges meeting at one vertex.
///
/// `FreeBoundary` as soon as one of the three edges has no opposite face
/// (`cxx:555-557`, the commented-out `F[0]` variant is dead in OCCT too).
/// Otherwise the three `ConcaveSide` calls fix the face orientations and the
/// comparison of `cxx:560-572` picks the state; `OnDiff` also covers a
/// `ConcaveSide` answer of `10` (one face ends in the extension of the other).
pub fn edge_state(edges: [&Edge; 3], ef_map: &ChFiDSMap) -> Result<State, String> {
    let (f1, f2) = conex_faces(edges[0], ef_map);
    let (f3, f4) = conex_faces(edges[1], ef_map);
    let (f5, f6) = conex_faces(edges[2], ef_map);

    if f2.is_none() || f4.is_none() || f6.is_none() {
        return Ok(State::FreeBoundary);
    }
    let f1 = f1.ok_or_else(|| "ChFi3d_EdgeState: vertex without faces".to_string())?;
    let f2 = f2.expect("checked above");
    let f3 = f3.ok_or_else(|| "ChFi3d_EdgeState: vertex without faces".to_string())?;
    let f4 = f4.expect("checked above");
    let f5 = f5.ok_or_else(|| "ChFi3d_EdgeState: vertex without faces".to_string())?;
    let f6 = f6.expect("checked above");

    let s1 = SurfaceAdaptor::of_face(&f1);
    let s2 = SurfaceAdaptor::of_face(&f2);
    let s3 = SurfaceAdaptor::of_face(&f3);
    let s4 = SurfaceAdaptor::of_face(&f4);
    let s5 = SurfaceAdaptor::of_face(&f5);
    let s6 = SurfaceAdaptor::of_face(&f6);

    let mut o01 = Orientation::Forward;
    let mut o02 = Orientation::Forward;
    let mut o11 = Orientation::Forward;
    let mut o12 = Orientation::Forward;
    let mut o21 = Orientation::Forward;
    let mut o22 = Orientation::Forward;
    let _ = concave_side(&s1, &s2, edges[0], &mut o01, &mut o02)?;
    // `cxx:558-559`: the two calls both write `i`, so what survives is the
    // second one — the one for `E[1]` — and that is what `cxx:564` tests.
    let i = concave_side(&s3, &s4, edges[1], &mut o11, &mut o12)?;
    let j = concave_side(&s5, &s6, edges[2], &mut o21, &mut o22)?;

    if o01 == o11 && o02 == o21 && o12 == o22 {
        Ok(State::AllSame)
    } else if o12 == o22 || i == 10 || j == 10 {
        Ok(State::OnDiff)
    } else {
        Ok(State::OnSame)
    }
}

/// `ChFi3d_evalconti(E, F1, F2)` (`cxx:581-599`): a cheap regularity probe.
///
/// `GeomAbs_G1` unless the two faces are the same *and* its surface is a cone,
/// a sphere or a torus, where the closing curve is `GeomAbs_CN`.
///
/// The `E` argument is unused in OCCT (`cxx:581` names it `/*E*/`), so it is
/// not part of the signature here.
pub fn eval_conti(f1: &Face, f2: &Face) -> Shape {
    let cont = Shape::G1;
    if !is_same(&f1.0, &f2.0) {
        return cont;
    }
    let forward = Face::wrap(f1.0.oriented(Orientation::Forward))
        .expect("orienting a face keeps it a face");
    let adaptor = SurfaceAdaptor::of_face(&forward);
    match adaptor.get_type() {
        GeomAbsSurfaceType::Cone
        | GeomAbsSurfaceType::Sphere
        | GeomAbsSurfaceType::Torus => Shape::CN,
        _ => cont,
    }
}
