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
use occt_geom2d::geom2d_api::{intersect_curves, project_point_on_curve};

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
pub fn find_pave_blocks<F: PaveFillerLike>(
    f: &F,
    the_pave_index: usize,
    the_face_info_index: usize,
    the_found: &mut Vec<BopdsPaveBlock>,
) {
    let Some(fi) = f.ds().face_info(the_face_info_index) else {
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
            if n1 == the_pave_index || n2 == the_pave_index {
                the_found.push(pb.clone());
            }
        }
    }
}

/// `AddSplitPoint`.
pub fn add_split_point(the_pbd: &mut BopdsPaveBlock, the_pave: BopdsPave, the_tol: f64) -> bool {
    let (a_td1, a_td2) = the_pbd.range();
    let a_t = the_pave.parameter();
    if a_t - a_td1 < the_tol || a_td2 - a_t < the_tol {
        return false;
    }
    if the_pbd.contains_parameter(a_t, the_tol).is_some() {
        return false;
    }
    the_pbd.append_ext_pave1(the_pave);
    true
}

/// `BOPAlgo_PaveFiller::FillPaves`.
pub fn fill_paves<F: PaveFillerLike>(
    f: &mut F,
    n_vd: usize,
    n_ed: usize,
    n_fd: usize,
    a_lpb_out: &[BopdsPaveBlock],
    a_pbd: &mut BopdsPaveBlock,
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
    let a_dv = Vertex(dv_s);
    let a_de = Edge(de_s);
    let a_df = Face(df_s);
    let a_tol_v = BRepTool::vertex_tolerance(&a_dv);
    let Some(surf) = BRepTool::face_surface(&a_df) else {
        return;
    };
    let Ok(c2d_de) = boptools_2d::make_2d(&a_de, &a_df) else {
        return;
    };
    let (td1, td2) = BRepTool::edge_parameters(&a_de);
    let p_mid = c2d_de.d0(0.5 * (td1 + td2));
    let (a_u_res, a_v_res) = surface_uv_resolution(surf.as_ref(), p_mid.x(), p_mid.y(), a_tol_v);
    let a_tol_int = PCONFUSION.max(a_u_res.max(a_v_res));
    let p1 = c2d_de.d0(td1);
    let p2 = c2d_de.d0(td2);
    let b_u_dir = (p1.y() - p2.y()).abs() < PCONFUSION;
    let a_tol_cmp = PCONFUSION.max(if b_u_dir { a_u_res } else { a_v_res });
    let mut a_pave = BopdsPave::new(n_vd, 0.0);
    for a_pb in a_lpb_out {
        let n_e = if a_pb.edge() != 0 {
            a_pb.edge()
        } else {
            a_pb.original_edge()
        };
        if n_e == 0 || n_e >= f.ds().nb_shapes() {
            continue;
        }
        let Some(es) = f.ds().shape(n_e).cloned() else {
            continue;
        };
        let a_e = Edge(es);
        let Ok(a_c2d) = boptools_2d::make_2d(&a_e, &a_df) else {
            continue;
        };
        let hits = if curve2d_is_line(a_c2d.as_ref()) {
            intersect_curves(c2d_de.as_ref(), a_c2d.as_ref(), a_tol_int)
        } else {
            let (t1, t2) = BRepTool::edge_parameters(&a_e);
            let _ = (t1, t2);
            intersect_curves(c2d_de.as_ref(), a_c2d.as_ref(), a_tol_int)
        };
        if hits.is_empty() {
            let a_t = if n_vd == a_pb.pave1().index() {
                a_pb.pave1().parameter()
            } else {
                a_pb.pave2().parameter()
            };
            let a_p2d = a_c2d.d0(a_t);
            if let Some(proj) = project_point_on_curve(c2d_de.as_ref(), &a_p2d, a_tol_int) {
                a_pave.set_parameter(proj.parameter);
                add_split_point(a_pbd, a_pave, a_tol_cmp);
            }
        } else {
            for hit in hits {
                a_pave.set_parameter(hit.u1);
                add_split_point(a_pbd, a_pave, a_tol_cmp);
            }
        }
    }
}

/// `MakeSplitEdge1`.
pub fn make_split_edge1(
    a_e: &Edge,
    a_f: &Face,
    a_v1: &Vertex,
    a_p1: f64,
    a_v2: &Vertex,
    a_p2: f64,
) -> Result<Edge, String> {
    let a_tol = 1.0e-7;
    let e = AlgoTools::make_split_edge(a_e, Some(&a_v1.0), a_p1, Some(&a_v2.0), a_p2)?;
    if let Some(mut g) = GeometryRegistry::global().edge_geom(&e.0) {
        g.degenerated = true;
        g.tolerance = a_tol;
        GeometryRegistry::global().set_edge(&e.0, g);
    }
    let face_key = GeometryRegistry::shape_key(&a_f.0);
    if let Some(pc) = GeometryRegistry::global()
        .edge_geom(&a_e.0)
        .and_then(|g| g.get_pcurve(face_key))
    {
        GeometryRegistry::global().set_edge_pcurve(&e.0, face_key, pc);
    }
    Ok(e)
}

/// `BOPAlgo_PaveFiller::MakeSplitEdge` (ProcessDE branch).
pub fn make_split_edge_de<F: PaveFillerLike>(f: &mut F, n_de: usize, n_df: usize) -> Result<(), String> {
    let Some(de_s) = f.ds().shape(n_de).cloned() else {
        return Ok(());
    };
    let Some(df_s) = f.ds().shape(n_df).cloned() else {
        return Ok(());
    };
    let mut a_de = Edge(de_s);
    a_de.0.set_orientation(Orientation::Forward);
    let a_df = Face(df_s);
    let a_lpb = f.ds().pave_blocks(n_de).to_vec();
    let a_nb_pb = a_lpb.len();
    let mut clear_all = false;
    let mut splits: Vec<(usize, usize)> = Vec::new();
    for (i, a_pb) in a_lpb.iter().enumerate() {
        let (n_v1, n_v2) = a_pb.indices();
        let (a_t1, a_t2) = a_pb.range();
        if f.ds().is_new_shape(n_v1) || a_nb_pb > 1 {
            let Some(v1s) = f.ds().shape(n_v1).cloned() else {
                continue;
            };
            let Some(v2s) = f.ds().shape(n_v2).cloned() else {
                continue;
            };
            let mut a_v1 = Vertex(v1s);
            a_v1.0.set_orientation(Orientation::Forward);
            let mut a_v2 = Vertex(v2s);
            a_v2.0.set_orientation(Orientation::Reversed);
            let a_sp = make_split_edge1(&a_de, &a_df, &a_v1, a_t1, &a_v2, a_t2)?;
            let mut a_si = BopdsShapeInfo::new(a_sp.0);
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
pub fn process_de<F: PaveFillerLike>(f: &mut F) -> Result<(), String> {
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
        let Some(a_sif) = f.ds().shape_info(n_f).cloned() else {
            continue;
        };
        let n_v0 = an_edge_info.sub_shapes().first().copied().unwrap_or(0);
        let n_v = f.ds().get_same_domain_index(n_v0);
        match a_sif.shape_type() {
            ShapeType::Face => {
                let mut a_lpb_out: Vec<BopdsPaveBlock> = Vec::new();
                find_pave_blocks(f, n_v, n_f, &mut a_lpb_out);
                if !a_lpb_out.is_empty() {
                    let pbd = {
                        let blocks = f.ds_mut().change_pave_blocks_mut(an_edge_index);
                        let Some(first) = blocks.first().cloned() else {
                            continue;
                        };
                        first
                    };
                    let mut pbd = pbd;
                    fill_paves(f, n_v, an_edge_index, n_f, &a_lpb_out, &mut pbd);
                    let mut out = Vec::new();
                    pbd.update(&mut out, true);
                    let blocks = f.ds_mut().change_pave_blocks_mut(an_edge_index);
                    blocks.clear();
                    blocks.extend(out);
                }
                make_split_edge_de(f, an_edge_index, n_f)?;
            }
            ShapeType::Edge => {
                let Some(a_de) = f.ds().shape(an_edge_index).cloned() else {
                    continue;
                };
                let Some(a_vn) = f.ds().shape(n_v).cloned() else {
                    continue;
                };
                let a_e = AlgoTools::make_split_edge(
                    &Edge(a_de),
                    Some(&a_vn),
                    0.0,
                    Some(&a_vn),
                    0.0,
                )?;
                if let Some(mut g) = GeometryRegistry::global().edge_geom(&a_e.0) {
                    g.degenerated = true;
                    g.tolerance = occt_core::precision::CONFUSION;
                    GeometryRegistry::global().set_edge(&a_e.0, g);
                }
                let a_si = BopdsShapeInfo::new(a_e.0);
                let n_en = f.ds_mut().append_info(a_si);
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
