//! `BOPAlgo_PaveFiller::MakePCurves`, `UpdateVertices`, `Prepare`
//! (`BOPAlgo_PaveFiller_7.cxx:589/808/850`).
//!
//! Common-block IN/ON edges get a p-curve on each FaceInfo face. ON edges
//! that already carry a p-curve are skipped; otherwise a sibling common-block
//! member that already has a p-curve supplies the 2D geometry. Section edges
//! of F/F interferences are then updated (vertices compared 3D vs 2D).
//! `Prepare` builds missing p-curves of edges on planar faces unless
//! non-destructive mode is on.
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use std::collections::HashSet;
use std::sync::Arc;

use occt_core::precision::PCONFUSION;
use occt_geom2d::curve::Curve2d;

use crate::abs::{Orientation, ShapeType};
use crate::algo_tools::D_TOLERANCE;
use crate::boptools_2d::{attach_existing_pcurve, build_pcurve_for_edge_on_face, has_curve_on_surface};
use crate::brep_surface::face_is_plane_type;
use crate::brep_tool::BRepTool;
use crate::iterator::ShapeIterator;
use crate::pave_blocks::PaveFillerLike;
use crate::pave_filler::PaveFiller;
use crate::pave_intersect::collect_pairs;
use crate::pcurve_full;
use crate::shape::{Edge, Face, TopoShape, Vertex};

use crate::topo_tools_full::edge_vertices;

/// `UpdateVertices` (`BOPAlgo_PaveFiller_7.cxx:808`).
pub fn update_vertices_pcurve(e: &Edge, f: &Face) {
    let mut ef = e.clone();
    ef.0.set_orientation(Orientation::Forward);
    let (v1, v2) = edge_vertices(&ef);
    let verts = [v1, v2];
    let Some(s) = BRepTool::face_surface(f) else {
        return;
    };
    let Some(c3d) = BRepTool::edge_curve(&ef) else {
        return;
    };
    let Some(c2d) = crate::boptools_2d::curve_on_surface(&ef, f) else {
        return;
    };
    let (t0, t1) = BRepTool::edge_parameters(&ef);
    let params = [t0, t1];
    for j in 0..2 {
        let Some(ref v) = verts[j] else {
            continue;
        };
        let tol_v2 = {
            let t = BRepTool::vertex_tolerance(v);
            t * t
        };
        let p3d = c3d.d0(params[j]);
        let p2d = c2d.d0(params[j]);
        let p3dx = s.d0(p2d.x(), p2d.y());
        let d2 = p3d.square_distance(&p3dx);
        if d2 > tol_v2 {
            let d = d2.sqrt();
            v.set_tolerance(d + D_TOLERANCE);
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
        let mut f1f = Face(fs);
        f1f.0.set_orientation(Orientation::Forward);
        for &(n_e, _, _) in fi.paves_in() {
            let Some(es) = f.ds().shape(n_e).cloned() else {
                continue;
            };
            build_or_warn(f, &Edge(es), &f1f, n_e, n_f1);
        }
        for &(n_e, t1, t2) in fi.paves_on() {
            let Some(es) = f.ds().shape(n_e).cloned() else {
                continue;
            };
            let e = Edge(es);
            if has_curve_on_surface(&e, &f1f) {
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
                            let ex = Edge(exs);
                            if !has_curve_on_surface(&ex, &f1f) {
                                continue;
                            }
                            if let Ok(pc) = crate::boptools_2d::make_2d(&ex, &f1f) {
                                let _ = attach_existing_pcurve(&e, &f1f, pc);
                                copied = true;
                                break;
                            }
                        }
                    }
                }
            }
            if !copied {
                build_or_warn(f, &e, &f1f, n_e, n_f1);
            } else {
                update_vertices_pcurve(&e, &f1f);
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
                    let e = Edge(es);
                    for m in 0..2 {
                        if b_pc[m] && an_ef_pairs.insert((n_e, n_f[m])) {
                            if let Some(ref face) = a_ff[m] {
                                if has_curve_on_surface(&e, face) {
                                    update_vertices_pcurve(&e, face);
                                } else {
                                    build_or_warn(f, &e, face, n_e, n_f[m]);
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

/// `IsBasedOnPlane` (`BOPAlgo_PaveFiller_7.cxx:936`) — the face surface, or its
/// basis when trimmed, is a `Geom_Plane`. Equivalent to
/// [`face_is_plane_type`]: the C++ down-cast to `Geom_Plane` is exactly the
/// `GeomAbs_Plane` adaptor type.
pub fn is_based_on_plane(f: &Face) -> bool {
    face_is_plane_type(f)
}

/// `BOPAlgo_PaveFiller::Prepare`.
pub fn prepare_pcurves_on_planes(f: &mut PaveFiller) -> Result<(), String> {
    if f.non_destructive() {
        return Ok(());
    }
    let mut mf: Vec<TopoShape> = Vec::new();
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
                mf.push(fs);
            }
        }
    }
    if mf.is_empty() {
        return Ok(());
    }
    for a_f in mf {
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
