//! `BRepOffset_Tool::CorrectOrientation` (`BRepOffset_Tool.cxx:4442-4505`) and
//! its file-local `IsInOut` (`:4422-4438`).
//!
//! `CorrectOrientation` walks every face of the initial shape `SI`, then every
//! image of that face in `InitOffset`, and only acts on offset faces whose
//! `AsDes` descendants contain a **new** (intersection) edge (`YaInt`,
//! `:4471-4478`). For those it builds a `BRepTopAdaptor_FClass2d` over the
//! initial face oriented `TopAbs_FORWARD` (`:4479-4480`), then walks the
//! descendants again: every new edge gets its p-curve on the offset face and,
//! when `Offset > 0` and the whole sampled p-curve lies `TopAbs_OUT` of the
//! initial face, the edge is reversed in place (`:4482-4502`).
//!
//! `IsInOut` (`:4422-4438`) samples the p-curve with
//! `GCPnts_QuasiUniformDeflection` at `Def = 100 * Precision::Confusion()`
//! and requires `FC.Perform(P) == S` at every sampled parameter. The `else`
//! arm of `Offset > 0` (`TopAbs_IN`, `:4496-4498`) is commented out in the
//! reference and therefore not ported.
//!
//! Port guards:
//!
//! * `AsDes->ChangeDescendant(OF)` (`:4463`) raises `Standard_NoSuchObject`
//!   for an offset face with no recorded descendants; the port skips it.
//! * `TopoDS::Face` / `TopoDS::Edge` (`:4456`, `:4469`) throw on a type
//!   mismatch; the port wraps through [`Face::wrap`] / [`Edge::wrap`] and
//!   skips the shape.
//! * `BRep_Tool::CurveOnSurface(E, OF, f, l)` (`:4488`) answers a null curve
//!   when the edge has no p-curve on that face; the port skips the edge
//!   instead of sampling a null adaptor.

use occt_core::precision::CONFUSION;
use occt_geom::gcpnts::QuasiUniformDeflection;
use occt_geom2d::Curve2d;

use crate::abs::{Orientation, ShapeType, State};
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_algo::AsDes;
use crate::brep_algo_image::BRepAlgoImage;
use crate::brep_offset::inter3d::ShapeIndexedMap;
use crate::brep_top_adaptor_fclass2d::BrepTopAdaptorFClass2d;
use crate::fclass2d::FaceState;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::explorer_shapes;

/// `TopAbs_State` of a [`FaceState`] (`TopAbs_State.hxx`).
fn state_of(s: FaceState) -> State {
    match s {
        FaceState::In => State::In,
        FaceState::Out => State::Out,
        FaceState::On => State::On,
        FaceState::Unknown => State::Unknown,
    }
}

/// File-local `IsInOut(FC, AC, S)` (`BRepOffset_Tool.cxx:4422-4438`).
///
/// `ac` is the `Geom2dAdaptor_Curve(CO2d, f, l)` of `:4488-4489`, i.e. the
/// p-curve plus the representation range.
fn is_in_out(
    fc: &BrepTopAdaptorFClass2d,
    ac: &dyn Curve2d,
    first: f64,
    last: f64,
    s: State,
) -> bool {
    // `constexpr double Def = 100 * Precision::Confusion()` (`:4426`).
    const DEF: f64 = 100.0 * CONFUSION;

    // `GCPnts_QuasiUniformDeflection QU(AC, Def)` (`:4427`): the `GeomAbs_C1`
    // default selects the derivative arm, `2` in the crate's encoding.
    let mut qu = QuasiUniformDeflection::new();
    qu.initialize_2d(ac, first, last, DEF, 2);

    for i in 1..=qu.nb_points() {
        let Some(param) = qu.parameter(i) else {
            continue;
        };
        // `gp_Pnt2d P = AC.Value(QU.Parameter(i))` (`:4431`).
        let p = ac.d0(param);
        // `if (FC.Perform(P) != S) return false` (`:4432`). `Perform` defaults
        // to `RecadreOnPeriodic = true` (`BRepTopAdaptor_FClass2d.hxx:40-41`).
        if state_of(fc.perform(p, true)) != s {
            return false;
        }
    }
    true
}

/// `BRepOffset_Tool::CorrectOrientation(SI, NewEdges, AsDes, InitOffset,
/// Offset)` (`BRepOffset_Tool.cxx:4442-4505`).
impl crate::brep_offset::BRepOffsetTool {
    pub fn correct_orientation(
        si: &TopoShape,
        new_edges: &ShapeIndexedMap,
        as_des: &mut AsDes,
        init_offset: &BRepAlgoImage,
        offset: f64,
    ) {
        // `cxx:4451`: `double f = 0., l = 0.` - the reference never reads the
        // previous iteration's value, so the pair lives with the p-curve lookup.
        // (Confirmed at `:4488`: `CurveOnSurface` always rewrites both.)

        for fi_shape in explorer_shapes(si, ShapeType::Face) {
            let Some(fi) = Face::wrap(fi_shape) else {
                continue;
            };
            // `cxx:4457`: `InitOffset.Image(FI)`; an unbound face answers an empty
            // list in `BRepAlgo_Image`.
            let lof = init_offset.image(&fi.0);

            for of_shape in lof {
                let Some(of) = Face::wrap(of_shape) else {
                    continue;
                };

                // `cxx:4463-4478`: does any descendant come from an intersection?
                let Some(loe) = as_des.change_descendant(&of.0) else {
                    continue;
                };
                let ya_int = loe.iter().any(|oe| new_edges.contains(oe));
                if !ya_int {
                    continue;
                }

                // `cxx:4479-4480`: the classifier over the initial face, forward.
                let a_local_face = Face(fi.0.oriented(Orientation::Forward));
                let fc = BrepTopAdaptorFClass2d::new(&a_local_face, CONFUSION);

                // `cxx:4482-4502`.
                for oe in loe.iter_mut() {
                    if !new_edges.contains(oe) {
                        continue;
                    }
                    let Some((co2d, first, last)) = curve_on_surface_range(&Edge(oe.clone()), &of)
                    else {
                        continue;
                    };
                    // `cxx:4488`'s out-parameters `f` / `l`.
                    let (f, l) = (first, last);
                    if offset > 0.0 && is_in_out(&fc, co2d.as_ref(), f, l, State::Out) {
                        // `OE.Reverse()` (`:4493`): `TopoDS_Shape::Reverse` flips
                        // the orientation of the element held by the AsDes list.
                        oe.reverse();
                    }
                }
            }
        }
    }
}
