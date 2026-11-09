//! `PerformPlanes` (`BRepOffset_Tool.cxx:4537-4602`): the untrimmed
//! plane/plane shortcut of `BRepOffset_Tool::Inter3D`.
//!
//! Rather than running the pave filler, the two infinite planes are
//! intersected with `IntTools_FaceFace` directly and the single resulting
//! section curve becomes one edge:
//!
//! 1. `IntTools_FaceFace::SetParameters(true, true, true, Confusion)` then
//!    `Perform(theFace1, theFace2)` (`cxx:4547-4548`);
//! 2. the first (and, for a plane/plane pair, only) `IntTools_Curve` supplies
//!    the 3D curve, its tolerance and its bounds (`cxx:4563-4572`);
//! 3. the two bound vertices are made at the curve ends, the last one
//!    `REVERSED`, and added to the edge (`cxx:4573-4580`);
//! 4. the curve's two pcurves are attached to `theFace1` / `theFace2`
//!    (`cxx:4582-4583`) and the edge range is written (`cxx:4585`);
//! 5. `OrientSection` gives the two orientations, `OUT` reverses both
//!    (`cxx:4589-4595`), `BRepLib::SameParameter` runs on the new edge
//!    (`cxx:4597`) and the oriented copies are appended (`cxx:4600-4601`).
//!
//! Port guards (no OCCT counterpart): a failed/empty `IntTools_FaceFace`
//! answers two empty lists, and a pcurve that the intersection did not
//! produce is simply not attached - `OrientSection` then refuses the edge, as
//! it does for the same null handle in OCCT.

use occt_core::precision::CONFUSION;

use crate::abs::{Orientation, State};
use crate::brep_lib_same_parameter::same_parameter_edge_inplace;
use crate::builder::TopoBuilder;
use crate::int_face_face::FaceFace;
use crate::shape::{Face, TopoShape};
use crate::tgeometry::GeometryRegistry;

use super::BRepOffsetTool;

impl BRepOffsetTool {
    /// `PerformPlanes(theFace1, theFace2, theSide, theL1, theL2)`
    /// (`cxx:4537-4602`).
    pub fn perform_planes(
        face1: &Face,
        face2: &Face,
        side: State,
        l1: &mut Vec<TopoShape>,
        l2: &mut Vec<TopoShape>,
    ) {
        l1.clear();
        l2.clear();

        // `cxx:4547-4548`.
        let mut face_face = FaceFace::new();
        face_face.set_parameters(true, true, true, CONFUSION);
        face_face.set_face1(face1.clone());
        face_face.set_face2(face2.clone());
        if face_face.perform().is_err() || !face_face.is_done() {
            return;
        }
        // `cxx:4553-4558`.
        let result = face_face.result();
        if result.is_empty() {
            return;
        }

        // `cxx:4563-4572`: `BRep_Builder::MakeEdge(aE, aIC.Curve(),
        // aIC.Tolerance())` and `IntTools_Curve::Bounds`.
        let curve = result.curve(0);
        let builder = TopoBuilder::new();
        let mut edge = builder.make_edge(curve.curve.clone(), curve.range.first, curve.range.last);
        let p_first = curve.curve.d0(curve.range.first);
        let p_last = curve.curve.d0(curve.range.last);

        // `cxx:4574-4580`: two bound vertices, the last one reversed.
        let v_first = builder.make_vertex(p_first, curve.tolerance);
        let mut v_last = builder.make_vertex(p_last, curve.tolerance);
        v_last.0.set_orientation(Orientation::Reversed);
        builder.add(&mut edge.0, &v_first.0);
        builder.add(&mut edge.0, &v_last.0);

        // `cxx:4582-4583`: `BRep_Builder::UpdateEdge(aE, theC2D, theF,
        // aIC.Tolerance())` - the pcurve plus the edge tolerance.
        let reg = GeometryRegistry::global();
        reg.set_edge_tolerance(&edge.0, curve.tolerance);
        if let Some(pc) = &curve.pcurve1 {
            reg.set_edge_pcurve(&edge.0, GeometryRegistry::shape_key(&face1.0), pc.clone());
        }
        if let Some(pc) = &curve.pcurve2 {
            reg.set_edge_pcurve(&edge.0, GeometryRegistry::shape_key(&face2.0), pc.clone());
        }

        // `cxx:4589-4595`.
        let Some((o1, o2)) = Self::orient_section(&edge, face1, face2) else {
            return;
        };
        let (o1, o2) = if side == State::Out {
            (o1.reversed(), o2.reversed())
        } else {
            (o1, o2)
        };

        // `cxx:4597`.
        same_parameter_edge_inplace(&edge, CONFUSION);

        // `cxx:4600-4601`.
        l1.push(edge.0.oriented(o1));
        l2.push(edge.0.oriented(o2));
    }

    /// `IsInf(theVal)` (`cxx:4606-4611`): `theVal > TheInfini * 0.9` with
    /// `TheInfini = 1.e+7` (`cxx:121`).
    pub fn is_inf(value: f64) -> bool {
        value > 1.0e7 * 0.9
    }
}
