//! `UpdateBlocksWithSharedVertices` using `IsVertexOnLine`.
//!
//! Source: `BOPAlgo_PaveFiller_6.cxx:3946-4052`.
//!
//! Non-destructive mode only. For each F/F pair with curves, collect old
//! (argument) vertices that sit On or In both faces, skip those that already
//! have an SD representative, and try `EstimatePaveOnCurve` (`IsVertexOnLine`)
//! on every section curve. A hit grows the vertex through `UpdateVertex` and
//! initialises pave blocks for it. Ends with
//! `UpdateCommonBlocksWithSDVertices`.

use std::collections::HashSet;

use crate::pave_common::update_common_blocks_with_sd_vertices;
use crate::pave_ff_pave_stick::estimate_pave_on_curve;
use crate::pave_filler::PaveFiller;
use crate::pave_force_ee::init_pave_blocks_for_vertex;
use crate::pave_update_sd::update_vertex;
use crate::shape::Vertex;
use crate::brep_tool::BRepTool;

/// `BOPAlgo_PaveFiller::UpdateBlocksWithSharedVertices` (`_6.cxx:3946`).
pub fn update_blocks_with_shared_vertices(f: &mut PaveFiller) -> Result<(), String> {
    if !f.non_destructive() {
        return Ok(());
    }
    let nb_ff = f.ds().interf_ff().len();
    if nb_ff == 0 {
        return Ok(());
    }
    let mut mf: HashSet<usize> = HashSet::new();
    let pairs: Vec<(usize, usize, usize)> = (0..nb_ff)
        .filter_map(|i| {
            let ff = &f.ds().interf_ff()[i];
            if ff.curves().is_empty() {
                return None;
            }
            let (n_f1, n_f2) = ff.indices();
            Some((i, n_f1, n_f2))
        })
        .collect();
    for (i, n_f1, n_f2) in pairs {
        if mf.insert(n_f1) {
            f.ds_mut().update_face_info_on(n_f1);
        }
        if mf.insert(n_f2) {
            f.ds_mut().update_face_info_on(n_f2);
        }
        let mut mi: HashSet<usize> = HashSet::new();
        let (on1, in1, on2, in2) = {
            let fi1 = f.ds().face_info(n_f1);
            let fi2 = f.ds().face_info(n_f2);
            (
                fi1.map(|x| x.verts_on().to_vec()).unwrap_or_default(),
                fi1.map(|x| x.verts_in().to_vec()).unwrap_or_default(),
                fi2.map(|x| x.verts_on().to_vec()).unwrap_or_default(),
                fi2.map(|x| x.verts_in().to_vec()).unwrap_or_default(),
            )
        };
        for &n_v in on1.iter().chain(in1.iter()) {
            if f.ds().is_new_shape(n_v) {
                continue;
            }
            if on2.contains(&n_v) || in2.contains(&n_v) {
                mi.insert(n_v);
            }
        }
        let nb_c = f.ds().interf_ff()[i].curves().len();
        let verts: Vec<usize> = mi.into_iter().collect();
        for j in 0..nb_c {
            let (tol_r3d, has_sd_skip) = {
                let nc = &f.ds().interf_ff()[i].curves()[j];
                (nc.tolerance().max(nc.tangential_tolerance()), ())
            };
            let _ = has_sd_skip;
            for &n_v in &verts {
                if f.ds().has_shape_sd(n_v).is_some() {
                    continue;
                }
                let on_curve = {
                    let nc = &f.ds().interf_ff()[i].curves()[j];
                    estimate_pave_on_curve(f.ds(), n_v, nc, tol_r3d)
                };
                if !on_curve {
                    continue;
                }
                let a_tol_v = f
                    .ds()
                    .shape(n_v)
                    .map(|s| BRepTool::vertex_tolerance(&Vertex(s.clone())))
                    .unwrap_or(0.0);
                let _ = update_vertex(f, n_v, a_tol_v)?;
                init_pave_blocks_for_vertex(f, n_v);
            }
        }
    }
    update_common_blocks_with_sd_vertices(f);
    Ok(())
}
