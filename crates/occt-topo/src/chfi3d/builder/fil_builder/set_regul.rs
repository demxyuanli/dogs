//! `ChFi3d_FilBuilder::SetRegul` — the regularity encoding of the fillet
//! result, run at the end of `Compute`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_FilBuilder.cxx:2595-2634`
//! (OCCT 8.0.0). The chamfer override (`ChFi3d_ChBuilder.cxx:2218-2274`) is a
//! different body — it measures the normals of the two faces with
//! `BRepAdaptor_Surface` / `BRepAdaptor_Curve2d` instead of calling
//! `ChFi3d_evalconti` — so this module is the fillet translation only.
//!
//! Translation notes:
//! - `myCoup->NewEdges(reg.Curve())` takes the DS *curve* index
//!   (`TopOpeBRepBuild_HBuilder.hxx:105`, `int`); it is the same
//!   `NewEdges(int)` the port exposes as
//!   [`crate::topopebrep::build::HBuilder::new_edges`],
//! - `myCoup->Merged(myDS->Shape(reg.S1()), TopAbs_IN)` uses the default
//!   `FindKeep = true` of `TopOpeBRepDS_DataStructure::Shape(int, bool)`,
//! - `TopoDS::Edge` / `TopoDS::Face` of a shape of another type raise in OCCT;
//!   the port reports `Err`,
//! - `BRep_Builder::Continuity(E, F1, F2, C)` is
//!   [`crate::brep_lib_continuity::set_stored_continuity`], and
//!   `ChFi3d_evalconti` is [`crate::chfi3d::builder_helpers::eval_conti`].

use crate::abs::State as TopAbsState;
use crate::brep_lib_continuity::set_stored_continuity;
use crate::chfi3d::builder_helpers::eval_conti;
use crate::chfi3d::builder::ChFi3dBuilder;
use crate::shape::{Edge, Face, TopoShape};

/// `SetRegul()` (`cxx:2595-2634`).
///
/// For every regular curve of `myRegul`, the continuity `ChFi3d_evalconti`
/// answers for its two supporting faces is stored on the new edge, exactly as
/// the OCCT loop does; a curve that the cut did not produce (`NewEdges` empty)
/// or a face pair the cut did not build (`NewFaces` / `Merged` empty) is
/// skipped, as in OCCT.
pub(crate) fn set_regul(base: &ChFi3dBuilder) -> Result<(), String> {
    for reg in &base.my_regul {
        // `itc.Initialize(myCoup->NewEdges(reg.Curve()));`
        // `if (itc.More())` (`cxx:2606-2609`).
        let Some(new_edge) = base.my_coup.new_edges(reg.curve()).first().cloned() else {
            continue;
        };
        // `TopoDS_Edge E = TopoDS::Edge(itc.Value());` (`cxx:2609`).
        let e = Edge::wrap(new_edge)
            .ok_or_else(|| "SetRegul: the new edge of the regular curve is not an edge".to_string())?;

        // `if (reg.IsSurface1()) its1.Initialize(myCoup->NewFaces(reg.S1()));`
        // `else its1.Initialize(myCoup->Merged(myDS->Shape(reg.S1()),
        // TopAbs_IN));` (`cxx:2610-2618`).
        let faces1: &[TopoShape] = if reg.is_surface1() {
            match base.my_coup.new_faces(reg.s1()) {
                Some(faces) => faces,
                None => continue,
            }
        } else {
            let s1 = match base.my_ds.shape(reg.s1(), true) {
                Some(s) => s,
                None => continue,
            };
            base.my_coup.merged(&s1, TopAbsState::In)
        };

        // `if (reg.IsSurface2()) its2.Initialize(myCoup->NewFaces(reg.S2()));`
        // `else its2.Initialize(myCoup->Merged(myDS->Shape(reg.S2()),
        // TopAbs_IN));` (`cxx:2619-2627`).
        let faces2: &[TopoShape] = if reg.is_surface2() {
            match base.my_coup.new_faces(reg.s2()) {
                Some(faces) => faces,
                None => continue,
            }
        } else {
            let s2 = match base.my_ds.shape(reg.s2(), true) {
                Some(s) => s,
                None => continue,
            };
            base.my_coup.merged(&s2, TopAbsState::In)
        };

        // `if (its1.More() && its2.More())` (`cxx:2628`).
        let (Some(f1), Some(f2)) = (faces1.first().cloned(), faces2.first().cloned()) else {
            continue;
        };
        // `TopoDS_Face F1 = TopoDS::Face(its1.Value());`
        // `TopoDS_Face F2 = TopoDS::Face(its2.Value());` (`cxx:2629-2631`).
        let f1 = Face::wrap(f1)
            .ok_or_else(|| "SetRegul: the first support of the regular curve is not a face".to_string())?;
        let f2 = Face::wrap(f2)
            .ok_or_else(|| "SetRegul: the second support of the regular curve is not a face".to_string())?;
        // `GeomAbs_Shape cont = ChFi3d_evalconti(E, F1, F2);`
        // `B.Continuity(E, F1, F2, cont);` (`cxx:2632-2633`).
        let cont = eval_conti(&f1, &f2);
        set_stored_continuity(&e, &f1, &f2, cont as u8);
    }
    Ok(())
}
