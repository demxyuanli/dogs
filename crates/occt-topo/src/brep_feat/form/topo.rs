//! The topological route of `BRepFeat_Form::GlobalPerform`
//! (`BRepFeat_Form.cxx:653-811`).
//!
//! The `From` / `Until` shapes that are not part of the basis become tool
//! solids (`BRepFeat::Tool`), the generated shape is cut by their compound, the
//! resulting solids become `theGShape`, and `myMap` is rewritten with the
//! `Modified` images of the limit faces.

use crate::abs::{Orientation, ShapeType};
use crate::brep_builder_api::make_compound_shapes;
use crate::brep_check;
use crate::brep_feat::algo_api::cut;
use crate::brep_feat::perf_selection::PerfSelection;
use crate::brep_feat::statics::tool;
use crate::brep_feat::status_error::StatusError;
use crate::shape::{Face, TopoShape};
use crate::topo_tools::{explorer_shapes, is_same, shapes_of};

use super::descendants::update_descendants_bop;
use super::{BRepFeatForm, FormCurves};

/// `cxx:653-811`: answers the `theGShape` the caller hands to
/// `BRepFeat_Builder::Init`.
///
/// `Ok(None)` is a `NotDone()` early return of the reference; the caller then
/// leaves `GlobalPerform` without touching the result.
pub(super) fn topological_route(
    form: &mut BRepFeatForm,
    _curves: &mut dyn FormCurves,
    from_in_shape: bool,
    until_in_shape: bool,
    ffrom: Option<&Face>,
    funtil: Option<&Face>,
    oriffrom: Orientation,
    orifuntil: Orientation,
) -> Result<Option<TopoShape>, String> {
    // `cxx:663`.
    let mut the_g_shape: Option<TopoShape> = form.my_gshape.clone();

    // `cxx:670-706`: the tool compound.
    let mut tools: Vec<TopoShape> = Vec::new();
    if form.my_sfrom.is_some() || form.my_suntil.is_some() {
        if !from_in_shape {
            if let Some(sfrom) = form.my_sfrom.clone() {
                let fac = ffrom.cloned().unwrap_or_else(Face::new);
                if let Some(s) = tool(&sfrom, &fac, oriffrom) {
                    tools.push(s.0);
                }
            }
        }
        if !until_in_shape {
            if let Some(suntil) = form.my_suntil.clone() {
                // `cxx:686-703`: with a `From` shape that is the same shape the
                // Until tool is skipped, otherwise it is added.
                let same_as_from = match &form.my_sfrom {
                    Some(sfrom) => is_same(sfrom, &suntil),
                    None => false,
                };
                let add = match &form.my_sfrom {
                    Some(_) => !same_as_from,
                    None => true,
                };
                if add {
                    let fac = funtil.cloned().unwrap_or_else(Face::new);
                    if let Some(s) = tool(&suntil, &fac, orifuntil) {
                        tools.push(s.0);
                    }
                }
            }
        }
    }

    // `cxx:708-723`: the selection kind is dropped when its limit is part of
    // the basis.
    if form.my_perf_selection == PerfSelection::SelectionU && !until_in_shape {
        form.my_perf_selection = PerfSelection::NoSelection;
    } else if form.my_perf_selection == PerfSelection::SelectionFU
        && !from_in_shape
        && !until_in_shape
    {
        form.my_perf_selection = PerfSelection::NoSelection;
    } else if form.my_perf_selection == PerfSelection::SelectionShU && !until_in_shape {
        form.my_perf_selection = PerfSelection::NoSelection;
    }

    // `cxx:725-799`.
    let comp = make_compound_shapes(&tools).0;
    let comp_has_solids = !shapes_of(&comp, ShapeType::Solid).is_empty();
    if comp_has_solids && form.my_gshape.is_some() {
        let gshape = form
            .my_gshape
            .clone()
            .ok_or_else(|| "StdFail_NotDone".to_string())?;
        let mut trp = cut(&gshape, &comp);
        // `cxx:726-727`: `BRepAlgoAPI_Cut` then `trP.Shape()`. A failed cut
        // raises `StdFail_NotDone` in the reference; the port treats it as the
        // empty solid list, which is the `EmptyCutResult` arm right below.
        let cut_shape = if trp.build().is_ok() {
            trp.shape().cloned()
        } else {
            None
        };
        let solids = cut_shape
            .as_ref()
            .map(|s| explorer_shapes(s, ShapeType::Solid))
            .unwrap_or_default();

        // `cxx:730-735`.
        if solids.is_empty() {
            form.my_status_error = StatusError::EmptyCutResult;
            form.not_done();
            return Ok(None);
        }

        // `cxx:737-749`: only the solids are preserved.
        let new_g = make_compound_shapes(&solids).0;
        if !brep_check::is_valid(&new_g) {
            form.my_status_error = StatusError::InvShape;
            form.not_done();
            return Ok(None);
        }

        // `cxx:750-773`: the From faces.
        if !from_in_shape {
            if let Some(sfrom) = form.my_sfrom.clone() {
                for fac in explorer_shapes(&sfrom, ShapeType::Face) {
                    form.my_map.bind(&fac, Vec::new());
                    if trp.is_deleted(&fac) {
                        continue;
                    }
                    let mut list = trp.modified(&fac);
                    if list.is_empty() {
                        list.push(fac.clone());
                    }
                    form.my_map.bind(&fac, list);
                }
            }
        }

        // `cxx:774-796`: the Until faces. The reference tests
        // `myMap.IsEmpty()` instead of the list it just bound (`cxx:790`), so
        // the fallback of this arm is the whole-map emptiness; that is
        // reproduced as written.
        if !until_in_shape {
            if let Some(suntil) = form.my_suntil.clone() {
                for fac in explorer_shapes(&suntil, ShapeType::Face) {
                    form.my_map.bind(&fac, Vec::new());
                    if trp.is_deleted(&fac) {
                        continue;
                    }
                    let mut list = trp.modified(&fac);
                    if form.my_map.is_empty() {
                        list.push(fac.clone());
                    }
                    form.my_map.bind(&fac, list);
                }
            }
        }

        // `cxx:799`.
        update_descendants_bop(form, &trp, &new_g, true)?;
        the_g_shape = Some(new_g);
    }

    Ok(the_g_shape)
}
