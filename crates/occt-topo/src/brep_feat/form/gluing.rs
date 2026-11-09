//! The gluing attempt of `BRepFeat_Form::GlobalPerform` (`cxx:332-600`).
//!
//! Every face bound in `myGluedF` must still be reproduced by the cut of the
//! generated shape, and the pair must be inside one another; only then is the
//! `LocOpe_Gluer` allowed to rebuild the shape locally. Any mismatch falls back
//! to the topological route.

use crate::abs::{Orientation, ShapeType};
use crate::brep_builder_api::make_compound_shapes;
use crate::brep_check;
use crate::brep_feat::algo_api::cut;
use crate::brep_feat::statics::{is_inside, tool};
use crate::locope::operation::LocOpeOperation;
use crate::locope::{FindEdges, Gluer, ShapeValueMap};
use crate::shape::{Face, TopoShape};
use crate::topo_tools::{faces_of, is_same, shapes_of};

use super::descendants::update_descendants_bop;
use super::BRepFeatForm;

/// `cxx:332-600`. Mutates `the_ope` / `change_ope` and the gluer; the local
/// `Collage` of the reference dies with the block.
#[allow(clippy::too_many_arguments)]
pub(super) fn try_gluing(
    form: &mut BRepFeatForm,
    the_glue: &mut Gluer,
    the_ope: &mut i32,
    change_ope: &mut bool,
    from_in_shape: bool,
    until_in_shape: bool,
    ffrom: Option<&Face>,
    funtil: Option<&Face>,
    oriffrom: Orientation,
    orifuntil: Orientation,
) -> Result<(), String> {
    // `cxx:339`.
    let mut collage = true;

    // A null base has no faces in the reference; the port has no null
    // `TopoDS_Shape`, so an unset member is the empty group.
    let base = form
        .my_sbase
        .clone()
        .unwrap_or_else(|| TopoShape::new(ShapeType::Compound));

    // `cxx:340-359`: the tool compound.
    let mut comp_solids: Vec<TopoShape> = Vec::new();
    if let Some(sfrom) = form.my_sfrom.clone() {
        let fac = ffrom.cloned().unwrap_or_else(Face::new);
        if let Some(s) = tool(&sfrom, &fac, oriffrom) {
            comp_solids.push(s.0);
        }
    }
    if let Some(suntil) = form.my_suntil.clone() {
        let fac = funtil.cloned().unwrap_or_else(Face::new);
        if let Some(s) = tool(&suntil, &fac, orifuntil) {
            comp_solids.push(s.0);
        }
    }

    // `cxx:361-363`.
    let mut the_fe = FindEdges::new();
    let mut locmap: ShapeValueMap<Vec<TopoShape>> = ShapeValueMap::new();

    // `cxx:365`: the cut of the generated shape by the tools.
    if !comp_solids.is_empty() && form.my_gshape.is_some() {
        let gshape = form
            .my_gshape
            .clone()
            .unwrap_or_else(|| TopoShape::new(ShapeType::Compound));
        let comp = make_compound_shapes(&comp_solids).0;
        let mut trp = cut(&gshape, &comp);
        let built = trp.build().is_ok();
        let solids: Vec<TopoShape> = if built {
            trp.shape().map(|s| shapes_of(s, ShapeType::Solid)).unwrap_or_default()
        } else {
            Vec::new()
        };

        if solids.is_empty() {
            // `cxx:369-374`.
            *the_ope = 2;
            *change_ope = true;
            collage = false;
        } else {
            // `cxx:376-378`: only the solids are preserved.
            let the_g_shape = make_compound_shapes(&solids).0;
            if !brep_check::is_valid(&the_g_shape) {
                // `cxx:384-389`.
                *the_ope = 2;
                *change_ope = true;
                collage = false;
            } else {
                // ---- `cxx:392-510`: else X1 ----
                bind_cut_faces(
                    form,
                    &mut locmap,
                    &trp,
                    from_in_shape,
                    until_in_shape,
                )?;

                // `cxx:470`.
                update_descendants_bop(form, &trp, &the_g_shape, true)?;

                // `cxx:472-509`.
                the_glue.init(&base, &the_g_shape);
                let glued: Vec<(TopoShape, TopoShape)> =
                    form.my_glued_f.iter().map(|(k, v)| (k.clone(), v.clone())).collect();
                for (gl, glface_s) in glued {
                    let mut ldsc: Vec<TopoShape> = Vec::new();
                    if !trp.is_deleted(&gl) {
                        ldsc = trp.modified(&gl);
                        if ldsc.is_empty() {
                            ldsc.push(gl.clone());
                        }
                    }
                    let glface = Face::wrap(glface_s)
                        .ok_or_else(|| "Standard_TypeMismatch: TopoDS::Face".to_string())?;
                    for fac_s in ldsc {
                        let fac = Face::wrap(fac_s)
                            .ok_or_else(|| "Standard_TypeMismatch: TopoDS::Face".to_string())?;
                        // `cxx:493`: the arguments are `(fac, glface)`.
                        collage = is_inside(&fac, &glface);
                        if !collage {
                            *the_ope = 2;
                            *change_ope = true;
                            break;
                        } else {
                            the_glue.bind_face(&fac, &glface)?;
                            the_fe.set(&fac.0, &glface.0);
                            bind_edges(&mut the_fe, the_glue)?;
                        }
                    }
                }
            }
        }
    } else {
        // ---- `cxx:512-547`: no cut ----
        let gshape = form
            .my_gshape
            .clone()
            .unwrap_or_else(|| TopoShape::new(ShapeType::Compound));
        the_glue.init(&base, &gshape);
        let gshape_faces = faces_of(&gshape);
        let glued: Vec<(TopoShape, TopoShape)> =
            form.my_glued_f.iter().map(|(k, v)| (k.clone(), v.clone())).collect();
        for (glface_s, fac_s) in glued {
            let glface = Face::wrap(glface_s)
                .ok_or_else(|| "Standard_TypeMismatch: TopoDS::Face".to_string())?;
            let fac = Face::wrap(fac_s)
                .ok_or_else(|| "Standard_TypeMismatch: TopoDS::Face".to_string())?;
            if gshape_faces.iter().any(|f| is_same(&f.0, &glface.0)) {
                // `cxx:529`: here the arguments are `(glface, fac)`.
                collage = is_inside(&glface, &fac);
                if !collage {
                    *the_ope = 2;
                    *change_ope = true;
                    break;
                } else {
                    the_glue.bind_face(&glface, &fac)?;
                    the_fe.set(&glface.0, &fac.0);
                    bind_edges(&mut the_fe, the_glue)?;
                }
            }
        }
    }

    // `cxx:549-568`: gluing on the From faces.
    if from_in_shape && collage {
        if let Some(sfrom) = form.my_sfrom.clone() {
            for fac2 in faces_of(&sfrom) {
                let list = locmap.find(&fac2.0).cloned().unwrap_or_default();
                for fac1_s in list {
                    let fac1 = Face::wrap(fac1_s)
                        .ok_or_else(|| "Standard_TypeMismatch: TopoDS::Face".to_string())?;
                    the_fe.set(&fac1.0, &fac2.0);
                    the_glue.bind_face(&fac1, &fac2)?;
                    bind_edges(&mut the_fe, the_glue)?;
                }
            }
        }
    }

    // `cxx:570-592`: gluing on the Until faces.
    if until_in_shape && collage {
        if let Some(suntil) = form.my_suntil.clone() {
            for fac2 in faces_of(&suntil) {
                let list = locmap.find(&fac2.0).cloned().unwrap_or_default();
                for fac1_s in list {
                    let fac1 = Face::wrap(fac1_s)
                        .ok_or_else(|| "Standard_TypeMismatch: TopoDS::Face".to_string())?;
                    the_glue.bind_face(&fac1, &fac2)?;
                    the_fe.set(&fac1.0, &fac2.0);
                    bind_edges(&mut the_fe, the_glue)?;
                }
            }
        }
    }

    // `cxx:594-599`.
    let ope = the_glue.ope_type();
    if ope == LocOpeOperation::Invalid
        || (form.my_fuse && ope != LocOpeOperation::Fuse)
        || (!form.my_fuse && ope != LocOpeOperation::Cut)
        || !collage
    {
        *the_ope = 2;
        *change_ope = true;
    }

    Ok(())
}

/// `cxx:393-467`: bind the `myMap` / `locmap` images of the From and Until
/// faces after the cut.
fn bind_cut_faces(
    form: &mut BRepFeatForm,
    locmap: &mut ShapeValueMap<Vec<TopoShape>>,
    trp: &crate::brep_feat::algo_api::BooleanOperation,
    from_in_shape: bool,
    until_in_shape: bool,
) -> Result<(), String> {
    // `cxx:393-428`.
    if let Some(sfrom) = form.my_sfrom.clone() {
        for fac in faces_of(&sfrom) {
            if !from_in_shape {
                form.my_map.bind(&fac.0, Vec::new());
            } else {
                locmap.bind(&fac.0, Vec::new());
            }
            if trp.is_deleted(&fac.0) {
                continue;
            }
            let mut list = trp.modified(&fac.0);
            if list.is_empty() {
                list.push(fac.0.clone());
            }
            if !from_in_shape {
                form.my_map.bind(&fac.0, list);
            } else {
                locmap.bind(&fac.0, list);
            }
        }
    }
    // `cxx:430-467`.
    if let Some(suntil) = form.my_suntil.clone() {
        for fac in faces_of(&suntil) {
            if !until_in_shape {
                form.my_map.bind(&fac.0, Vec::new());
            } else {
                locmap.bind(&fac.0, Vec::new());
            }
            if trp.is_deleted(&fac.0) {
                continue;
            }
            let mut list = trp.modified(&fac.0);
            if list.is_empty() {
                list.push(fac.0.clone());
            }
            if !until_in_shape {
                form.my_map.bind(&fac.0, list);
            } else {
                locmap.bind(&fac.0, list);
            }
        }
    }
    Ok(())
}

/// `for (theFE.InitIterator(); theFE.More(); theFE.Next()) theGlue.Bind(...)`
/// (`cxx:502-506`, `:536-540`, `:561-565`, `:583-587`).
fn bind_edges(the_fe: &mut FindEdges, the_glue: &mut Gluer) -> Result<(), String> {
    the_fe.init_iterator();
    while the_fe.more() {
        let from = the_fe.edge_from().cloned();
        let to = the_fe.edge_to().cloned();
        if let (Some(from), Some(to)) = (from, to) {
            the_glue.bind_edge(&from, &to)?;
        }
        the_fe.next();
    }
    Ok(())
}
