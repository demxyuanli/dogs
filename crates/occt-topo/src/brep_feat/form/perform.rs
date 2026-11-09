//! `BRepFeat_Form::GlobalPerform` (`BRepFeat_Form.cxx:59-1240`).
//!
//! The dispatch of the two routes lives here:
//!
//! 1. the field validity check (`cxx:68-77`) and the choice of the route
//!    (`cxx:81-104`: gluing when faces were glued, topological otherwise);
//! 2. the `From` / `Until` containment tests against the basis (`cxx:109-158`)
//!    and the intersection scan that fixes `sens`, `FFrom`, `FUntil` and the two
//!    orientations (`cxx:160-325`);
//! 3. the gluing attempt ([`try_gluing`], `cxx:332-599`) and the `LocOpe_Gluer`
//!    result (`cxx:603-638`);
//! 4. the topological route ([`topological_route`], `cxx:653-811`);
//! 5. the tool-part selection ([`select_parts`], `cxx:812-1211`) and the final
//!    `PerformResult` (`cxx:1214-1237`).

use crate::abs::{Orientation, ShapeType};
use crate::brep_check;
use crate::brep_feat::builder::BRepFeatBuilder;
use crate::brep_feat::perf_selection::PerfSelection;
use crate::brep_feat::status_error::StatusError;
use crate::locope::{CSIntersector, Gluer};
use crate::shape::Face;
use crate::topo_tools::{explorer_shapes, is_same};

use super::descendants::update_descendants_gluer;
use super::gluing::try_gluing;
use super::selection::select_parts;
use super::topo::topological_route;
use super::{BRepFeatForm, FormCurves};

/// `GlobalPerform()` (`cxx:59-1240`).
pub fn global_perform(
    form: &mut BRepFeatForm,
    curves: &mut dyn FormCurves,
) -> Result<(), String> {
    // `cxx:68-77`.
    if !(form.my_sb_ok
        && form.my_gs_ok
        && form.my_sf_ok
        && form.my_su_ok
        && form.my_gf_ok
        && form.my_sk_ok
        && form.my_ps_ok)
    {
        form.my_status_error = StatusError::NotInitialized;
        form.not_done();
        return Ok(());
    }

    // `cxx:81-104`.
    let mut the_ope = 2;
    if form.my_just_feat && !form.my_fuse {
        form.my_status_error = StatusError::InvOption;
        form.not_done();
        return Ok(());
    } else if form.my_just_feat {
        the_ope = 2;
    } else if !form.my_glued_f.is_empty() {
        the_ope = 1;
    }
    let mut change_ope = false;

    // `cxx:109-134`.
    let mut from_in_shape = false;
    if let Some(sfrom) = form.my_sfrom.clone() {
        from_in_shape = true;
        let base_faces = match &form.my_sbase {
            Some(sbase) => explorer_shapes(sbase, ShapeType::Face),
            None => Vec::new(),
        };
        for ffrom in explorer_shapes(&sfrom, ShapeType::Face) {
            if !base_faces.iter().any(|x| is_same(x, &ffrom)) {
                from_in_shape = false;
                break;
            }
        }
    }

    // `cxx:136-158`.
    let mut until_in_shape = false;
    if let Some(suntil) = form.my_suntil.clone() {
        until_in_shape = true;
        let base_faces = match &form.my_sbase {
            Some(sbase) => explorer_shapes(sbase, ShapeType::Face),
            None => Vec::new(),
        };
        for funtil in explorer_shapes(&suntil, ShapeType::Face) {
            if !base_faces.iter().any(|x| is_same(x, &funtil)) {
                until_in_shape = false;
                break;
            }
        }
    }

    // `cxx:160-178`.
    let mut sens = 0;
    let scur = curves.curves();
    let mut orifuntil = Orientation::Internal;
    let mut oriffrom = Orientation::Internal;
    let mut f_from: Option<Face> = None;
    let mut f_until: Option<Face> = None;

    // `cxx:171-193`: the two intersection scans.
    let mut asi1 = CSIntersector::new();
    let mut asi2 = CSIntersector::new();
    if let Some(sfrom) = form.my_sfrom.clone() {
        asi1.init(&sfrom);
        let _ = asi1.perform_curves(&scur);
    }
    if let Some(suntil) = form.my_suntil.clone() {
        asi2.init(&suntil);
        let _ = asi2.perform_curves(&scur);
    }

    // `cxx:196-325`: find `sens`, `FFrom` and `FUntil`.
    for jj in 1..=scur.len() {
        if asi1.is_done() && asi2.is_done() {
            let nb_from = asi1.nb_points(jj).unwrap_or(0);
            if nb_from == 0 {
                continue;
            }
            let mf = asi1.point(jj, 1).map(|p| p.parameter()).unwrap_or(0.0);
            let m_f = asi1
                .point(jj, nb_from)
                .map(|p| p.parameter())
                .unwrap_or(0.0);
            let nb_until = asi2.nb_points(jj).unwrap_or(0);
            if nb_until == 0 {
                continue;
            }
            let mu = asi2.point(jj, 1).map(|p| p.parameter()).unwrap_or(0.0);
            let m_u = asi2
                .point(jj, nb_until)
                .map(|p| p.parameter())
                .unwrap_or(0.0);

            // `cxx:212`: a null curve answers "not periodic", where the
            // reference dereferences the handle.
            let periodic = scur
                .get(jj - 1)
                .and_then(|c| c.as_ref())
                .map(|c| c.is_periodic())
                .unwrap_or(false);
            if !periodic {
                let kf;
                let ku;
                if mu <= m_f && mf <= m_u {
                    // overlapping intervals
                    sens = 1;
                    kf = 1;
                    ku = nb_until;
                } else if mu > m_f {
                    if sens == -1 {
                        form.my_status_error = StatusError::IntervalOverlap;
                        form.not_done();
                        return Ok(());
                    }
                    sens = 1;
                    kf = 1;
                    ku = nb_until;
                } else {
                    if sens == 1 {
                        form.my_status_error = StatusError::IntervalOverlap;
                        form.not_done();
                        return Ok(());
                    }
                    sens = -1;
                    kf = nb_from;
                    ku = 1;
                }
                // `cxx:246-258`.
                if oriffrom == Orientation::Internal {
                    if let Ok(p) = asi1.point(jj, kf) {
                        let mut oript = p.orientation();
                        if oript.is_forward() || oript.is_reversed() {
                            if sens == -1 {
                                oript = oript.reversed();
                            }
                            oriffrom = oript.reversed();
                            f_from = Some(p.face().clone());
                        }
                    }
                }
                // `cxx:259-271`.
                if orifuntil == Orientation::Internal {
                    if let Ok(p) = asi2.point(jj, ku) {
                        let mut oript = p.orientation();
                        if oript.is_forward() || oript.is_reversed() {
                            if sens == -1 {
                                oript = oript.reversed();
                            }
                            orifuntil = oript;
                            f_until = Some(p.face().clone());
                        }
                    }
                }
                // The bindings shadow the reference's `kf`/`ku`, which are only
                // read inside this arm.
                let _ = (kf, ku);
            }
        } else if asi2.is_done() {
            // `cxx:273-321`.
            let nb_until = asi2.nb_points(jj).unwrap_or(0);
            if nb_until == 0 {
                continue;
            }
            if sens != 1 {
                let p1 = asi2.point(jj, 1).map(|p| p.parameter()).unwrap_or(0.0);
                let p2 = asi2
                    .point(jj, nb_until)
                    .map(|p| p.parameter())
                    .unwrap_or(0.0);
                if p1 * p2 <= 0.0 {
                    sens = 1;
                } else if p1 < 0.0 {
                    sens = -1;
                } else {
                    sens = 1;
                }
            }
            let ku = if sens == -1 { 1 } else { nb_until };
            if orifuntil == Orientation::Internal && sens != 0 {
                if let Ok(p) = asi2.point(jj, ku) {
                    let mut oript = p.orientation();
                    if oript.is_forward() || oript.is_reversed() {
                        if sens == -1 {
                            oript = oript.reversed();
                        }
                        orifuntil = oript;
                        f_until = Some(p.face().clone());
                    }
                }
            }
        } else {
            // `cxx:322-325`.
            sens = 1;
            break;
        }
    }

    // `cxx:328-599`: the gluing attempt.
    let mut the_glue = Gluer::new();
    if the_ope == 1 {
        try_gluing(
            form,
            &mut the_glue,
            &mut the_ope,
            &mut change_ope,
            from_in_shape,
            until_in_shape,
            f_from.as_ref(),
            f_until.as_ref(),
            oriffrom,
            orifuntil,
        )?;
    }

    // `cxx:603-638`.
    if the_ope == 1 {
        // A failed `Perform` leaves `IsDone` false, which is the fallback below.
        let _ = the_glue.perform();
        if the_glue.is_done() {
            if let Ok(shshs) = the_glue.resulting_shape().map(|s| s.clone()) {
                if brep_check::is_valid(&shshs) {
                    update_descendants_gluer(form, &the_glue)?;
                    form.my_new_edges = the_glue.edges().to_vec();
                    form.my_tgt_edges = the_glue.tgt_edges().to_vec();
                    form.done();
                    form.my_shape = Some(shshs);
                } else {
                    the_ope = 2;
                    change_ope = true;
                }
            } else {
                the_ope = 2;
                change_ope = true;
            }
        } else {
            the_ope = 2;
            change_ope = true;
        }
    }

    // `cxx:641-650`.
    if the_ope == 2 && change_ope && form.my_just_gluer {
        form.my_just_gluer = false;
        the_ope = 0;
    }

    // `cxx:653-811`.
    if the_ope == 2 {
        let the_g_shape = match topological_route(
            form,
            curves,
            from_in_shape,
            until_in_shape,
            f_from.as_ref(),
            f_until.as_ref(),
            oriffrom,
            orifuntil,
        )? {
            Some(s) => s,
            None => return Ok(()),
        };

        // `cxx:801-811`.
        let b_flag = form.my_perf_selection != PerfSelection::NoSelection;
        let sbase = form
            .my_sbase
            .clone()
            .ok_or_else(|| "Standard_ConstructionError".to_string())?;
        let mut the_builder = BRepFeatBuilder::new();
        the_builder.init_with_tool(&sbase, &the_g_shape);
        the_builder.set_operation_with_flag(if form.my_fuse { 1 } else { 0 }, b_flag);
        // The reference ignores a failed `Perform`; so does the port.
        let _ = the_builder.perform();
        let lshape = the_builder.parts_of_tool();

        // `cxx:812-1211`.
        if !select_parts(form, curves, &mut the_builder, &lshape, sens)? {
            return Ok(());
        }

        // `cxx:1214-1231`.
        if !form.my_just_feat {
            if b_flag {
                the_builder.perform_result()?;
                form.my_shape = Some(the_builder.shape().clone());
            } else {
                form.my_shape = Some(the_builder.shape().clone());
            }
            form.done();
        } else {
            form.done();
        }
    }

    // `cxx:1239`.
    form.my_status_error = StatusError::Ok;
    Ok(())
}
