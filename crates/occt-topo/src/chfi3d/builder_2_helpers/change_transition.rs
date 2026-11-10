//! `ChangeTransition` — flip the transition of a common point when its arc
//! does not cross the surface.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_2.cxx:592-627`
//! (OCCT 8.0.0), a file static of the anonymous namespace:
//!
//! ```text
//! static void ChangeTransition(const ChFiDS_CommonPoint& Precedant,
//!                              ChFiDS_CommonPoint&       Courant,
//!                              int                       FaceIndex,
//!                              const handle<TopOpeBRepDS_HDataStructure>& DS)
//! ```
//!
//! Both common points sit on the same support face, so OCCT only has to tell a
//! "cutting" arc (an arc and its reversed copy carry *different* p-curves on
//! `F`) from an ordinary one; on a cutting arc the transition is flipped only
//! when the arc's tangent is parallel to the preceding point's vector
//! (`cxx:611-618`), which is the `!` of what the surrounding code would call a
//! crossing.
//!
//! Translation notes:
//! - `getCurveOnSurface(E, F)` (`cxx:79-83`) is
//!   [`crate::boptools_2d::curve_on_surface`] (the p-curve handle without the
//!   parameter range), so `PCurve1 != PCurve2` — a *handle* comparison
//!   (`cxx:610`) — becomes [`Arc::ptr_eq`] over the two `Option`s,
//! - `Arc.Reversed()` is [`TopoShape::oriented`] with the flipped orientation,
//!   the same TShape, which is what `TopoDS_Shape::Reversed` keeps,
//! - `BRepAdaptor_Curve AC(Arc); AC.D1(u, P, tgarc)` is
//!   [`CurveAdaptor::of_edge`] then `d1`,
//! - `TopAbs::Reverse(T)` (`TopAbs.hxx:88-96`) is [`Orientation::reversed`],
//! - the OCCT null shape / null p-curve dereferences are the `Err` cases.

use std::sync::Arc;

use occt_core::gp::GpVec;
use occt_core::precision::CONFUSION;
use occt_geom2d::curve::Curve2d;

use crate::boptools_2d::curve_on_surface;
use crate::chfi_ds::CommonPoint;
use crate::shape::{Edge, Face, TopoShape};
use crate::topopebrep::ds::DataStructure;
use crate::topopebrep::tool::adaptor::CurveAdaptor;

/// True when the two handles stand for one p-curve, the `PCurve1 != PCurve2`
/// test of `cxx:610` (two null handles are equal, as in OCCT).
fn same_pcurve(a: Option<&Arc<dyn Curve2d>>, b: Option<&Arc<dyn Curve2d>>) -> bool {
    match (a, b) {
        (None, None) => true,
        (Some(a), Some(b)) => Arc::ptr_eq(a, b),
        _ => false,
    }
}

/// `ChangeTransition(Precedant, Courant, FaceIndex, DS)` (`cxx:592-627`).
pub fn change_transition(
    precedant: &CommonPoint,
    courant: &mut CommonPoint,
    face_index: i32,
    ds: &DataStructure,
) -> Result<(), String> {
    // `const TopoDS_Face& F = TopoDS::Face(DS->Shape(FaceIndex));` (`cxx:604`).
    let shape = ds
        .shape(face_index, false)
        .ok_or_else(|| "ChangeTransition: DS->Shape(FaceIndex) is null".to_string())?;
    let face = Face::wrap(shape)
        .ok_or_else(|| "ChangeTransition: DS->Shape(FaceIndex) is not a face".to_string())?;
    // `const TopoDS_Edge& Arc = Precedant.Arc();` (`cxx:605`).
    let arc = precedant
        .arc()
        .ok_or_else(|| "ChangeTransition: Precedant is not on an arc".to_string())?
        .clone();
    let arc_edge = Edge::wrap(arc.clone())
        .ok_or_else(|| "ChangeTransition: Precedant.Arc() is not an edge".to_string())?;

    // `cxx:607-609`.
    let pcurve1 = curve_on_surface(&arc_edge, &face);
    let rev_arc: TopoShape = arc.oriented(arc.orientation().reversed());
    let pcurve2 = Edge::wrap(rev_arc.clone())
        .and_then(|rev_edge| curve_on_surface(&rev_edge, &face));

    // `cxx:610-618`.
    let mut tochange = true;
    if !same_pcurve(pcurve1.as_ref(), pcurve2.as_ref()) {
        let parameter = precedant
            .parameter_on_arc()
            .ok_or_else(|| "ChangeTransition: Precedant has no arc parameter".to_string())?;
        let (_p, tgarc) = CurveAdaptor::of_edge(&arc_edge)
            .d1(parameter)
            .ok_or_else(|| "ChangeTransition: null curve on the arc".to_string())?;
        let vector: &GpVec = precedant
            .vector()
            .ok_or_else(|| "ChFiDS_CommonPoint::Vector (ChFi3d_Builder_2.cxx:617)".to_string())?;
        // `tgarc.IsParallel(Precedant.Vector(), Precision::Confusion())`.
        tochange = tgarc.is_parallel_ang(vector, CONFUSION);
    }

    // `cxx:620-626`.
    if tochange {
        let parameter = precedant
            .parameter_on_arc()
            .ok_or_else(|| "ChangeTransition: Precedant has no arc parameter".to_string())?;
        let transition = precedant
            .transition_on_arc()
            .ok_or_else(|| "ChangeTransition: Precedant has no arc transition".to_string())?;
        courant.set_arc(CONFUSION, &arc, parameter, transition.reversed());
    }
    Ok(())
}
