//! `BOPAlgo_PaveFiller::MakePCurves`, `UpdateVertices`, `Prepare`
//! (`BOPAlgo_PaveFiller_7.cxx:589/808/850`).
//!
//! Common-block IN/ON edges get a p-curve on each FaceInfo face. ON edges
//! that already carry a p-curve are skipped; otherwise a sibling common-block
//! member that already has a p-curve supplies the 2D geometry. Section edges
//! of F/F interferences are then updated (vertices compared 3D vs 2D).
//! `Prepare` builds missing p-curves of edges on planar faces unless
//! non-destructive mode is on.

use std::collections::HashSet;
use std::sync::Arc;

use occt_core::precision::PCONFUSION;
use occt_geom2d::curve::Curve2d;

use crate::abs::{Orientation, ShapeType};
use crate::algo_tools::D_TOLERANCE;
use crate::boptools_2d::{attach_existing_pcurve, build_pcurve_for_edge_on_face, has_curve_on_surface};
use crate::brep_surface::face_is_planar;
use crate::brep_tool::BRepTool;
use crate::iterator::ShapeIterator;
use crate::pave_blocks::PaveFillerLike;
use crate::pave_filler::PaveFiller;
use crate::pave_intersect::collect_pairs;
use crate::pcurve_full;
use crate::shape::{Edge, Face, TopoShape, Vertex};

use crate::topo_tools_full::edge_vertices;

/// `UpdateVertices` (`BOPAlgo_PaveFiller_7.cxx:808`).
pub fn update_vertices_pcurve(a_e: &Edge, a_f: &Face) {
    let mut a_ef = a_e.clone();
    a_ef.0.set_orientation(Orientation::Forward);
    let (v1, v2) = edge_vertices(&a_ef);
    let verts = [v1, v2];
    let Some(a_s) = BRepTool::face_surface(a_f) else {
        return;
    };
    let Some(a_c3d) = BRepTool::edge_curve(&a_ef) else {
        return;
    };
    let Some(a_c2d) = crate::boptools_2d::curve_on_surface(&a_ef, a_f) else {
        return;
    };
    let (t0, t1) = BRepTool::edge_parameters(&a_ef);
    let params = [t0, t1];
    for j in 0..2 {
        let Some(ref a_v) = verts[j] else {
            continue;
        };
        let a_tol_v2 = {
            let t = BRepTool::vertex_tolerance(a_v);
            t * t
        };
        let a_p3d = a_c3d.d0(params[j]);
        let a_p2d = a_c2d.d0(params[j]);
        let a_p3dx = a_s.d0(a_p2d.x(), a_p2d.y());
        let a_d2 = a_p3d.square_distance(&a_p3dx);
        if a_d2 > a_tol_v2 {
            let a_d = a_d2.sqrt();
            a_v.set_tolerance(a_d + D_TOLERANCE);
        }
    }
}

fn build_or_warn(f: &mut PaveFiller, edge: &Edge, face: &Face, n_e: usize, n_f: usize) {
    match build_pcurve_for_edge_on_face(edge, face) {
        Ok(_) => update_vertices_pcurve(edge, face),
        Err(e) => f.add_warning(format!(
            "make_pcurves: building p-curve failed for edge {n_e} on face {n_f}: {e}"
        )),
    }
}

/// `BOPAlgo_PaveFiller::MakePCurves`.
pub fn make_p_curves(f: &mut PaveFiller) -> Result<(), String> {
    if f.avoid_build_pcurve() {
        return Ok(());
    }
    let (on_s1, on_s2) = f.section_pcurve_on();
    if !on_s1 && !on_s2 {
        return Ok(());
    }
    let fi_pool = f.ds().face_info_pool().to_vec();
    for fi in &fi_pool {
        let n_f1 = fi.index();
        let Some(fs) = f.ds().shape(n_f1).cloned() else {
            continue;
        };
        let mut a_f1f = Face(fs);
        a_f1f.0.set_orientation(Orientation::Forward);
        for &(n_e, _, _) in fi.paves_in() {
            let Some(es) = f.ds().shape(n_e).cloned() else {
                continue;
            };
            build_or_warn(f, &Edge(es), &a_f1f, n_e, n_f1);
        }
        for &(n_e, t1, t2) in fi.paves_on() {
            let Some(es) = f.ds().shape(n_e).cloned() else {
                continue;
            };
            let a_e = Edge(es);
            if has_curve_on_surface(&a_e, &a_f1f) {
                continue;
            }
            let mut copied = false;
            if let Some(pb) = f.ds().pave_blocks(n_e).iter().find(|p| {
                let (a, b) = p.range();
                (a - t1).abs() <= PCONFUSION && (b - t2).abs() <= PCONFUSION
            }) {
                if let Some(cb) = f.ds().common_block(pb) {
                    if cb.pave_blocks().len() >= 2 {
                        for pbx in cb.pave_blocks() {
                            if pb_same(pbx, pb) {
                                continue;
                            }
                            let n_ex = pbx.original_edge();
                            let Some(exs) = f.ds().shape(n_ex).cloned() else {
                                continue;
                            };
                            let a_ex = Edge(exs);
                            if !has_curve_on_surface(&a_ex, &a_f1f) {
                                continue;
                            }
                            if let Ok(pc) = crate::boptools_2d::make_2d(&a_ex, &a_f1f) {
                                let _ = attach_existing_pcurve(&a_e, &a_f1f, pc);
                                copied = true;
                                break;
                            }
                        }
                    }
                }
            }
            if !copied {
                build_or_warn(f, &a_e, &a_f1f, n_e, n_f1);
            } else {
                update_vertices_pcurve(&a_e, &a_f1f);
            }
        }
    }

    let b_pc = [on_s1, on_s2];
    if b_pc[0] || b_pc[1] {
        let mut an_ef_pairs: HashSet<(usize, usize)> = HashSet::new();
        let ffs = f.ds().interf_ff().to_vec();
        for ff in ffs {
            if ff.curves().is_empty() {
                continue;
            }
            let (n_f0, n_f1) = ff.indices();
            let n_f = [n_f0, n_f1];
            let mut a_ff = [None, None];
            for m in 0..2 {
                if let Some(s) = f.ds().shape(n_f[m]).cloned() {
                    let mut face = Face(s);
                    face.0.set_orientation(Orientation::Forward);
                    a_ff[m] = Some(face);
                }
            }
            for nc in ff.curves() {
                for pb in nc.pave_blocks() {
                    let n_e = pb.edge();
                    let Some(es) = f.ds().shape(n_e).cloned() else {
                        continue;
                    };
                    let a_e = Edge(es);
                    for m in 0..2 {
                        if b_pc[m] && an_ef_pairs.insert((n_e, n_f[m])) {
                            if let Some(ref face) = a_ff[m] {
                                if has_curve_on_surface(&a_e, face) {
                                    update_vertices_pcurve(&a_e, face);
                                } else {
                                    build_or_warn(f, &a_e, face, n_e, n_f[m]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

fn pb_same(a: &crate::bopds::BopdsPaveBlock, b: &crate::bopds::BopdsPaveBlock) -> bool {
    crate::pave_ff_exist::pb_key(a) == crate::pave_ff_exist::pb_key(b)
}

/// `IsBasedOnPlane` (`BOPAlgo_PaveFiller_7.cxx:936`).
pub fn is_based_on_plane(a_f: &Face) -> bool {
    face_is_planar(a_f)
}

/// `BOPAlgo_PaveFiller::Prepare`.
pub fn prepare_pcurves_on_planes(f: &mut PaveFiller) -> Result<(), String> {
    if f.non_destructive() {
        return Ok(());
    }
    let mut a_mf: Vec<TopoShape> = Vec::new();
    let mut seen: HashSet<usize> = HashSet::new();
    let types = [ShapeType::Vertex, ShapeType::Edge, ShapeType::Face];
    for t in types {
        for (_, n_f) in collect_pairs(f.ds(), t, ShapeType::Face) {
            if !seen.insert(n_f) {
                continue;
            }
            let Some(fs) = f.ds().shape(n_f).cloned() else {
                continue;
            };
            let a_f = Face(fs.clone());
            if is_based_on_plane(&a_f) {
                a_mf.push(fs);
            }
        }
    }
    if a_mf.is_empty() {
        return Ok(());
    }
    for a_f in a_mf {
        let face = Face(a_f.clone());
        for sub in ShapeIterator::of_shape(&a_f) {
            if sub.shape_type() != ShapeType::Edge {
                continue;
            }
            let a_e = Edge(sub);
            if has_curve_on_surface(&a_e, &face) {
                continue;
            }
            match pcurve_full::make_pcurve_full(&a_e, &face) {
                Ok(pc) => {
                    let _ = attach_existing_pcurve(&a_e, &face, pc);
                }
                Err(e) => f.add_warning(format!("prepare: pcurve on plane failed: {e}")),
            }
        }
    }
    Ok(())
}

/// `BOPTools_AlgoTools2D::BuildPCurveForEdgeOnFace` wrapper used by MakePCurves.
pub fn ensure_pcurve(edge: &Edge, face: &Face) -> Result<Arc<dyn Curve2d>, String> {
    build_pcurve_for_edge_on_face(edge, face)
}

/// Keep `Vertex` reachable for UpdateVertices callers that pass DS shapes.
pub fn vertex_pair_of(edge: &Edge) -> (Option<Vertex>, Option<Vertex>) {
    edge_vertices(edge)
}
