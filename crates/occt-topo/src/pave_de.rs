//! Remainder of `BOPAlgo_PaveFiller_8.cxx`: `ProcessDE`, `FindPaveBlocks`,
//! `FillPaves`, `MakeSplitEdge`, `MakeSplitEdge1`, `AddSplitPoint`.
//!
//! Degenerated edges whose `ShapeInfo` carries `HasFlag(nF)` are split by
//! intersecting their 2D curve against every In/On/Section pave block that
//! passes through the degenerated vertex. Intersection tolerance is the
//! surface U/V resolution of the vertex tolerance (not a constant).

use occt_core::precision::PCONFUSION;
use occt_geom::Surface;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::geom2d_api::{intersect_curves, project_point_on_curve_2d as project_point_on_curve};

use crate::abs::{Orientation, ShapeType};
use crate::algo_tools::AlgoTools;
use crate::bopds::{BopdsPave, BopdsPaveBlock, BopdsShapeInfo};
use crate::boptools_2d;
use crate::brep_tool::BRepTool;
use crate::pave_blocks::PaveFillerLike;
use crate::shape::{Edge, Face, Vertex};
use crate::tgeometry::GeometryRegistry;

fn surface_uv_resolution(s: &dyn Surface, u: f64, v: f64, res: f64) -> (f64, f64) {
    let (_, du, dv) = s.d1(u, v);
    let ur = if du.magnitude() > 1.0e-12 {
        res / du.magnitude()
    } else {
        res
    };
    let vr = if dv.magnitude() > 1.0e-12 {
        res / dv.magnitude()
    } else {
        res
    };
    (ur, vr)
}

fn curve2d_is_line(c: &dyn Curve2d) -> bool {
    let (a, b) = (c.first_parameter(), c.last_parameter());
    if !a.is_finite() || !b.is_finite() || (b - a).abs() <= 1e-15 {
        return false;
    }
    let p0 = c.d0(a);
    let p1 = c.d0(0.5 * (a + b));
    let p2 = c.d0(b);
    let vx1 = p1.x() - p0.x();
    let vy1 = p1.y() - p0.y();
    let vx2 = p2.x() - p0.x();
    let vy2 = p2.y() - p0.y();
    let cross = (vx1 * vy2 - vy1 * vx2).abs();
    let scale = (vx1 * vx1 + vy1 * vy1).sqrt() * (vx2 * vx2 + vy2 * vy2).sqrt();
    cross < 1e-6 * scale.max(1e-12)
}

/// `BOPAlgo_PaveFiller::FindPaveBlocks`.
pub fn de_find_pave_blocks<F: PaveFillerLike>(
    f: &F,
    pave_index: usize,
    face_info_index: usize,
    found: &mut Vec<BopdsPaveBlock>,
) {
    let Some(fi) = f.ds().face_info(face_info_index) else {
        return;
    };
    let tuples: Vec<(usize, f64, f64)> = fi
        .paves_in()
        .iter()
        .chain(fi.paves_on().iter())
        .chain(fi.paves().iter())
        .copied()
        .collect();
    for (edge, _t1, _t2) in tuples {
        for pb in f.ds().pave_blocks(edge) {
            let (n1, n2) = pb.indices();
            if n1 == pave_index || n2 == pave_index {
                found.push(pb.clone());
            }
        }
    }
}

/// `AddSplitPoint`.
pub fn add_split_point(pbd: &mut BopdsPaveBlock, pave: BopdsPave, tol: f64) -> bool {
    let (td1, td2) = pbd.range();
    let t = pave.parameter();
    if t - td1 < tol || td2 - t < tol {
        return false;
    }
    if pbd.contains_parameter(t, tol).is_some() {
        return false;
    }
    pbd.append_ext_pave1(pave);
    true
}

/// `BOPAlgo_PaveFiller::FillPaves`.
pub fn de_fill_paves<F: PaveFillerLike>(
    f: &mut F,
    n_vd: usize,
    n_ed: usize,
    n_fd: usize,
    lpb_out: &[BopdsPaveBlock],
    pbd: &mut BopdsPaveBlock,
) {
    let Some(dv_s) = f.ds().shape(n_vd).cloned() else {
        return;
    };
    let Some(de_s) = f.ds().shape(n_ed).cloned() else {
        return;
    };
    let Some(df_s) = f.ds().shape(n_fd).cloned() else {
        return;
    };
    let dv = Vertex(dv_s);
    let de = Edge(de_s);
    let df = Face(df_s);
    let tol_v = BRepTool::vertex_tolerance(&dv);
    let Some(surf) = BRepTool::face_surface(&df) else {
        return;
    };
    let Ok(c2d_de) = boptools_2d::make_2d(&de, &df) else {
        return;
    };
    let (td1, td2) = BRepTool::edge_parameters(&de);
    let p_mid = c2d_de.d0(0.5 * (td1 + td2));
    let (u_res, v_res) = surface_uv_resolution(surf.as_ref(), p_mid.x(), p_mid.y(), tol_v);
    let tol_int = PCONFUSION.max(u_res.max(v_res));
    let p1 = c2d_de.d0(td1);
    let p2 = c2d_de.d0(td2);
    let b_u_dir = (p1.y() - p2.y()).abs() < PCONFUSION;
    let tol_cmp = PCONFUSION.max(if b_u_dir { u_res } else { v_res });
    let mut pave = BopdsPave::new(n_vd, 0.0);
    for pb in lpb_out {
        let n_e = if pb.edge() != 0 {
            pb.edge()
        } else {
            pb.original_edge()
        };
        if n_e == 0 || n_e >= f.ds().nb_shapes() {
            continue;
        }
        let Some(es) = f.ds().shape(n_e).cloned() else {
            continue;
        };
        let e = Edge(es);
        let Ok(c2d) = boptools_2d::make_2d(&e, &df) else {
            continue;
        };
        let hits = if curve2d_is_line(c2d.as_ref()) {
            intersect_curves(c2d_de.as_ref(), c2d.as_ref(), tol_int)
        } else {
            let (t1, t2) = BRepTool::edge_parameters(&e);
            let _ = (t1, t2);
            intersect_curves(c2d_de.as_ref(), c2d.as_ref(), tol_int)
        };
        if hits.is_empty() {
            let t = if n_vd == pb.pave1().index() {
                pb.pave1().parameter()
            } else {
                pb.pave2().parameter()
            };
            let p2d = c2d.d0(t);
            if let Some(proj) = project_point_on_curve(c2d_de.as_ref(), &p2d, tol_int) {
                pave.set_parameter(proj.parameter);
                add_split_point(pbd, pave, tol_cmp);
            }
        } else {
            for hit in hits {
                pave.set_parameter(hit.u1);
                add_split_point(pbd, pave, tol_cmp);
            }
        }
    }
}

/// `MakeSplitEdge1`.
pub fn make_split_edge1(
    a_e: &Edge,
    f: &Face,
    v1: &Vertex,
    p1: f64,
    v2: &Vertex,
    p2: f64,
) -> Result<Edge, String> {
    let tol = 1.0e-7;
    let e = AlgoTools::make_split_edge(a_e, Some(&v1.0), p1, Some(&v2.0), p2)?;
    if let Some(mut g) = GeometryRegistry::global().edge_geom(&e.0) {
        g.degenerated = true;
        g.tolerance = tol;
        GeometryRegistry::global().set_edge(&e.0, g);
    }
    let face_key = GeometryRegistry::shape_key(&f.0);
    if let Some(pc) = GeometryRegistry::global().edge_pcurve(&a_e.0, face_key) {
        GeometryRegistry::global().set_edge_pcurve(&e.0, face_key, pc);
    }
    Ok(e)
}

/// `BOPAlgo_PaveFiller::MakeSplitEdge` (ProcessDE branch).
pub fn de_make_split_edge<F: PaveFillerLike>(f: &mut F, n_de: usize, n_df: usize) -> Result<(), String> {
    let Some(de_s) = f.ds().shape(n_de).cloned() else {
        return Ok(());
    };
    let Some(df_s) = f.ds().shape(n_df).cloned() else {
        return Ok(());
    };
    let mut de = Edge(de_s);
    de.0.set_orientation(Orientation::Forward);
    let df = Face(df_s);
    let lpb = f.ds().pave_blocks(n_de).to_vec();
    let nb_pb = lpb.len();
    let mut clear_all = false;
    let mut splits: Vec<(usize, usize)> = Vec::new();
    for (i, a_pb) in lpb.iter().enumerate() {
        let (n_v1, n_v2) = a_pb.indices();
        let (t1, t2) = a_pb.range();
        if f.ds().is_new_shape(n_v1) || nb_pb > 1 {
            let Some(v1s) = f.ds().shape(n_v1).cloned() else {
                continue;
            };
            let Some(v2s) = f.ds().shape(n_v2).cloned() else {
                continue;
            };
            let mut v1 = Vertex(v1s);
            v1.0.set_orientation(Orientation::Forward);
            let mut v2 = Vertex(v2s);
            v2.0.set_orientation(Orientation::Reversed);
            let sp = make_split_edge1(&de, &df, &v1, t1, &v2, t2)?;
            let mut a_si = BopdsShapeInfo::new(sp.0);
            a_si.change_sub_shapes().extend_from_slice(&[n_v1, n_v2]);
            let n_sp = f.ds_mut().append_info(a_si);
            splits.push((i, n_sp));
        } else {
            clear_all = true;
            break;
        }
    }
    if clear_all {
        if let Some(si) = f.ds_mut().change_shape_info(n_de) {
            si.pb_reference = -1;
        }
        f.ds_mut().change_pave_blocks_mut(n_de).clear();
        return Ok(());
    }
    let blocks = f.ds_mut().change_pave_blocks_mut(n_de);
    for (i, n_sp) in splits {
        if let Some(pb) = blocks.get_mut(i) {
            pb.set_edge(n_sp);
        }
    }
    Ok(())
}

/// `BOPAlgo_PaveFiller::ProcessDE`.
pub fn de_process<F: PaveFillerLike>(f: &mut F) -> Result<(), String> {
    let n = f.ds().nb_source_shapes();
    for an_edge_index in 0..n {
        let Some(an_edge_info) = f.ds().shape_info(an_edge_index).cloned() else {
            continue;
        };
        if an_edge_info.shape_type() != ShapeType::Edge {
            continue;
        }
        let Some(n_f) = an_edge_info.flag() else {
            continue;
        };
        let Some(sif) = f.ds().shape_info(n_f).cloned() else {
            continue;
        };
        let n_v0 = an_edge_info.sub_shapes().first().copied().unwrap_or(0);
        let n_v = f.ds().get_same_domain_index(n_v0);
        match sif.shape_type() {
            ShapeType::Face => {
                let mut lpb_out: Vec<BopdsPaveBlock> = Vec::new();
                de_find_pave_blocks(f, n_v, n_f, &mut lpb_out);
                if !lpb_out.is_empty() {
                    let pbd = {
                        let blocks = f.ds_mut().change_pave_blocks_mut(an_edge_index);
                        let Some(first) = blocks.first().cloned() else {
                            continue;
                        };
                        first
                    };
                    let mut pbd = pbd;
                    de_fill_paves(f, n_v, an_edge_index, n_f, &lpb_out, &mut pbd);
                    let mut out = Vec::new();
                    pbd.update(&mut out, true);
                    let blocks = f.ds_mut().change_pave_blocks_mut(an_edge_index);
                    blocks.clear();
                    blocks.extend(out);
                }
                de_make_split_edge(f, an_edge_index, n_f)?;
            }
            ShapeType::Edge => {
                let Some(de) = f.ds().shape(an_edge_index).cloned() else {
                    continue;
                };
                let Some(vn) = f.ds().shape(n_v).cloned() else {
                    continue;
                };
                let e = AlgoTools::make_split_edge(
                    &Edge(de),
                    Some(&vn),
                    0.0,
                    Some(&vn),
                    0.0,
                )?;
                if let Some(mut g) = GeometryRegistry::global().edge_geom(&e.0) {
                    g.degenerated = true;
                    g.tolerance = occt_core::precision::CONFUSION;
                    GeometryRegistry::global().set_edge(&e.0, g);
                }
                let si = BopdsShapeInfo::new(e.0);
                let n_en = f.ds_mut().append_info(si);
                let blocks = f.ds_mut().change_pave_blocks_mut(an_edge_index);
                if let Some(a_pbd) = blocks.first_mut() {
                    a_pbd.set_edge(n_en);
                }
            }
            _ => {}
        }
    }
    Ok(())
}
