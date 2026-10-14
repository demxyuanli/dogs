//! `ShapeFix_Face::Perform` second-part loop and the `TopoDS_Iterator` helpers it uses, with the `ShapeBuild_ReShape` orientation restore. Source: `ShapeFix_Face.cxx:365-716`, `ShapeBuild_ReShape.cxx:317`.
use super::*;

/// Direct child edges of a wire **without** composing the wire orientation
/// (`TopoDS_Iterator` with `cumOri=false`). Used for `TranslateEdgeLoop`
/// SelectForwardSeam `EdgeO` (`cxx:700`).
pub(super) fn edges_stored_on_wire(wire: &Wire) -> Vec<Edge> {
    wire.0
        .tshape
        .read()
        .expect("poisoned TShape lock")
        .children
        .iter()
        .filter(|s| s.shape_type() == ShapeType::Edge)
        .cloned()
        .map(Edge)
        .collect()
}

/// Direct child wires of a face **without** composing the face orientation
/// (`TopoDS_Iterator(face, false)`) — the iterator both wire rounds of
/// `ShapeFix_Face::Perform` use: `ShapeFix_Face.cxx:400` and `:532`
/// (`for (TopoDS_Iterator iter(S, false); ...)`, `S = Context()->Apply(myFace)`).
/// The stored wire keeps the orientation `ShapeBuild_ReShape::applyImpl`
/// restores at `ShapeBuild_ReShape.cxx:317`; composing the face's orientation
/// into it would hand `ShapeFix_Wire` a wire whose `ShapeExtend_WireData::Init`
/// (`ShapeExtend_WireData.cxx:114-121`) prepends instead of appends, i.e. the
/// opposite `ShapeAnalysis_WireOrder` input.
pub(super) fn wires_stored_on_face(face: &Face) -> Vec<Wire> {
    face.0
        .tshape
        .read()
        .expect("poisoned TShape lock")
        .children
        .iter()
        .filter(|s| s.shape_type() == ShapeType::Wire)
        .cloned()
        .map(Wire)
        .collect()
}

/// `ShapeFix_Face::Perform`'s second-part loop, `FixOrientation` step
/// (`ShapeFix_Face.cxx:687-697`): the loop points the single `ShapeFix_Face` at
/// each face of `myResult` (`cxx:503-505`) and, for every one of them, runs
/// `FixOrientation(MapWires)` when `NeedFix(myFixOrientationMode)` holds. The
/// rebuilt face is recorded with `Context()->Replace(myFace, S)`
/// (`cxx:1639-1641`) and materialised later by `ShapeFix_Shape`'s
/// `Context()->Apply` (`ShapeFix_Shape.cxx:257`), which is why the caller has to
/// take the face back off the tool.
pub(super) fn fix_face_orientation(sff: &mut crate::shhealing::ShapeFixFace, face: &Face) -> Face {
    if !sff.fix_orientation_mode {
        return face.clone(); // cxx:691-697: `NeedFix(myFixOrientationMode)`.
    }
    sff.face = Some(face.clone()); // cxx:505: `myFace = TopoDS::Face(exp.Current())`.
    sff.fix_orientation();
    sff.face.clone().unwrap_or_else(|| face.clone())
}

/// `ShapeFix_Face::Perform`'s second-part loop, the wire steps that precede
/// `FixOrientation` (`ShapeFix_Face.cxx:676-685`): `FixWiresTwoCoincEdges`
/// (`cxx:676-679`) then, under `NeedFix(myFixIntersectingWiresMode)`,
/// `FixIntersectingWires` (`cxx:680-685`).
///
/// Both substitutes rebuild the face and record `Context()->Replace`
/// (`cxx:2894`, `cxx:2526`), so the current face has to be taken back off the
/// tool before `FixOrientation` runs on it.
pub(super) fn fix_face_wire_pairs(sff: &mut crate::shhealing::ShapeFixFace, face: &Face) -> Face {
    sff.face = Some(face.clone()); // cxx:505
    sff.fix_wires_two_coinc_edges(); // cxx:676-679, `ShapeExtend_DONE7`
    if sff.fix_intersecting_wires_mode {
        sff.fix_intersecting_wires(); // cxx:680-685, `ShapeExtend_DONE6`
    }
    sff.face.clone().unwrap_or_else(|| face.clone())
}

/// `ShapeFix_Face::Perform`'s second-part loop for one face of `myResult`
/// (`ShapeFix_Face.cxx:676-716`): the intersected-wire pair steps
/// (`cxx:676-685`), `FixOrientation(MapWires)` (`cxx:687-697`),
/// `FixAddNaturalBound` (`cxx:702-708`, whose `true` turns `NeedSplit` off) and
/// `FixSplitFace` (`cxx:710-716`, gated on `NeedFix(myFixSplitFaceMode)`,
/// `NeedSplit` and `MapWires.Extent() > 1`).
///
/// Returns `myResult` for this face, which `FixSplitFace` replaces by a
/// compound (`cxx:2993`); every other step keeps it a face.
///
/// UNPORTED: `NeedCheckSplitWire` / `SplitWire` (`cxx:629-674`) and the first
/// wire round (`cxx:365-480`, `ShapeFix_Wire::Perform`) run in the caller, on
/// the substituted wire data.
pub(super) fn fix_face_second_part(
    sff: &mut crate::shhealing::ShapeFixFace,
    face: &Face,
) -> crate::shape::TopoShape {
    let face = fix_face_wire_pairs(sff, face); // cxx:676-685
    let face = fix_face_orientation(sff, &face); // cxx:687-697
    // cxx:700: `BRepTools::Update(myFace)` is a no-op in the port.
    // cxx:702-708.
    let mut need_split = true;
    sff.face = Some(face.clone());
    if sff.fix_add_natural_bound() {
        need_split = false;
    }
    let face = sff.face.clone().unwrap_or(face);
    // cxx:710-716.
    if sff.fix_split_face_mode && need_split && sff.map_wires.len() > 1 {
        sff.face = Some(face.clone());
        if sff.fix_split_face() {
            if let Some(res) = sff.result.clone() {
                return res;
            }
        }
    }
    face.0.clone()
}
