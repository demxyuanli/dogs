#[cfg(test)]
use super::prelude::*;
#[cfg(test)]
use super::*;

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Vertex/Vertex intersection: fuse coincident vertices of the arguments into
/// same-domain vertices.
///
/// Source: `BOPAlgo_PaveFiller::PerformVV` (`BOPAlgo_PaveFiller_1.cxx`). The
/// implementations live in the `pave_*` phase modules; these re-exports keep
/// the `pave_intersect::vertex_face` path working.
pub use crate::pave_vv::perform_vv;

/// Vertex/Vertex intersection restricted to the given pairs (used by the
/// repeat-intersection stage — `RepeatIntersection`).
pub(crate) use crate::pave_vv::perform_vv_pairs;

/// Vertex/Edge intersection: project every interfering vertex onto its edge and
/// insert paves at the projections, splitting the affected pave blocks.
///
/// Source: `BOPAlgo_PaveFiller::PerformVE` / `IntersectVE`
/// (`BOPAlgo_PaveFiller_2.cxx`).
pub use crate::pave_ve::perform_ve;

/// Vertex/Edge intersection restricted to the given pairs (used by the
/// repeat-intersection stage — `RepeatIntersection`).
pub(crate) use crate::pave_ve::perform_ve_pairs;

/// Edge/Edge intersection: intersect every interfering pair of pave blocks,
/// creating new vertices at the crossing points and recording coincident edge
/// common blocks.
///
/// Source: `BOPAlgo_PaveFiller::PerformEE` + `PerformNewVertices` /
/// `TreatNewVertices` (`BOPAlgo_PaveFiller_3.cxx`).
pub use crate::pave_ee_perform::perform_ee;

/// Vertex/Face intersection: classify every interfering vertex against its face
/// and record the vertices that project onto the face (inside or on the
/// boundary).
///
/// Source: `BOPAlgo_PaveFiller::PerformVF` + `TreatVerticesEE`
/// (`BOPAlgo_PaveFiller_4.cxx`).
pub use crate::pave_vf::perform_vf;

/// Vertex/Face intersection restricted to the given pairs (used by the
/// repeat-intersection stage — `RepeatIntersection`).
pub(crate) use crate::pave_vf::perform_vf_pairs;

// ---------------------------------------------------------------------------
// DS-only cores (kept for the `pave_intersect` unit tests)
//
// The live pipeline drives these through `crate::pave_ee_perform` /
// `crate::pave_split_blocks`, which own the filler handle; the DS-only cores
// below stay as the isolated-test entry points.
// ---------------------------------------------------------------------------

/// A new vertex to create at an edge/edge (or edge/face) intersection point.
#[cfg(test)]
#[derive(Debug, Clone)]
pub(crate) struct NewVertexSeed {
    /// 3-D point of the intersection.
    pub point: GpPnt,
    /// Tolerance of the new vertex.
    pub tol: f64,
    /// First hit edge index and the parameter on it.
    pub edge_a: usize,
    pub t_a: f64,
    /// Second hit edge index and the parameter on it.
    pub edge_b: usize,
    pub t_b: f64,
}

/// Cluster coincident seeds by point proximity.
#[cfg(test)]
pub(super) fn cluster_seeds(seeds: &[NewVertexSeed], fuzzy: f64) -> Vec<Vec<usize>> {
    let mut clusters: Vec<Vec<usize>> = Vec::new();
    let mut cluster_pts: Vec<GpPnt> = Vec::new();
    for (i, s) in seeds.iter().enumerate() {
        let tol = s.tol + fuzzy + CONFUSION;
        let tol2 = tol * tol;
        let mut placed = false;
        for (k, cpt) in cluster_pts.iter().enumerate() {
            if cpt.square_distance(&s.point) <= tol2 {
                clusters[k].push(i);
                placed = true;
                break;
            }
        }
        if !placed {
            clusters.push(vec![i]);
            cluster_pts.push(s.point);
        }
    }
    clusters
}

/// Fuse coincident new-vertex seeds, create one DS vertex per cluster and add
/// the vertex as an extra pave to every hit edge pave block.
///
/// Port of `BOPAlgo_PaveFiller::TreatNewVertices` + the append part of
/// `PerformNewVertices`. Returns the DS indices of the created vertices.
#[cfg(test)]
pub(crate) fn treat_new_vertices(
    ds: &mut BopdsDS,
    fuzzy: f64,
    seeds: &[NewVertexSeed],
) -> Result<Vec<usize>, String> {
    if seeds.is_empty() {
        return Ok(Vec::new());
    }
    let clusters = cluster_seeds(seeds, fuzzy);
    let mut out = Vec::new();
    for cluster in clusters {
        let mut cx: f64 = 0.0;
        let mut cy: f64 = 0.0;
        let mut cz: f64 = 0.0;
        let mut tol: f64 = 0.0;
        for &i in &cluster {
            let s = &seeds[i];
            cx += s.point.x();
            cy += s.point.y();
            cz += s.point.z();
            tol = tol.max(s.tol);
        }
        let n = cluster.len() as f64;
        let p = GpPnt::new(cx / n, cy / n, cz / n);
        // Reuse an existing source vertex at the same point (e.g. a box corner
        // of the other operand) instead of creating a coincident new one —
        // otherwise the two vertices split the same edge at the same parameter
        // into a zero-length block.
        let n_v = match find_source_vertex_at(ds, &p, tol + CONFUSION) {
            Some(existing) => existing,
            None => {
                let v = AlgoTools::make_new_vertex(&p, tol)?;
                ds.append(v)?
            }
        };
        for &i in &cluster {
            let s = &seeds[i];
            add_ext_pave(ds, s.edge_a, s.t_a, n_v);
            if s.edge_b != s.edge_a {
                add_ext_pave(ds, s.edge_b, s.t_b, n_v);
            }
        }
        out.push(n_v);
    }
    Ok(out)
}

/// The DS index of an existing source vertex within `tol` of `p`, if any.
#[cfg(test)]
pub(super) fn find_source_vertex_at(ds: &BopdsDS, p: &GpPnt, tol: f64) -> Option<usize> {
    let n = ds.nb_source_shapes();
    let tol2 = tol * tol;
    (0..n).find(|&i| {
        ds.shape_info(i).map(|si| si.shape_type() == ShapeType::Vertex).unwrap_or(false)
            && ds.shape(i).and_then(|s| {
                if s.shape_type() == ShapeType::Vertex {
                    Some(crate::brep_tool::BRepTool::vertex_point(&Vertex(s.clone())))
                } else {
                    None
                }
            })
            .map(|q| q.square_distance(p) <= tol2)
            .unwrap_or(false)
    })
}

/// Core of [`split_pave_blocks`]: split the pave blocks of the given edges that
/// carry extra paves, unifying block vertices when a block has no valid range.
///
/// When `edges` is empty every source edge with update-able pave blocks is
/// processed.
#[cfg(test)]
pub(crate) fn split_pave_blocks_impl(
    ds: &mut BopdsDS,
    ctx: &mut FillCtx,
    edges: &[usize],
) -> Result<(), String> {
    let to_split: Vec<usize> = if edges.is_empty() {
        (0..ds.nb_source_shapes())
            .filter(|&i| {
                ds.shape_info(i).map(|s| s.shape_type()) == Some(ShapeType::Edge)
                    && ds.pave_blocks(i).iter().any(|pb| pb.is_to_update())
            })
            .collect()
    } else {
        edges.to_vec()
    };
    let mut unify_fence: HashSet<(usize, usize)> = HashSet::new();

    for n_e in to_split {
        if !ds.has_pave_blocks(n_e) {
            continue;
        }
        let old = ds.pave_blocks(n_e).to_vec();
        let mut new_list: Vec<BopdsPaveBlock> = Vec::new();
        for mut pb in old {
            if !pb.is_to_update() {
                new_list.push(pb);
                continue;
            }
            let mut out: Vec<BopdsPaveBlock> = Vec::new();
            pb.update(&mut out, true);
            for mut new_pb in out {
                // Resolve bound indices to their SD representatives and recompute
                // the shrunk range for the elementary block.
                update_pb_with_sd_vertices(ds, &mut new_pb);
                fill_shrunk_data(ds, &mut new_pb);
                let (ts1, ts2, _) = new_pb.shrunk_data();
                let has_valid = (ts2 - ts1) > PCONFUSION;
                let splittable = new_pb.is_splittable();
                let b_check_dist = has_valid && !splittable;
                if !has_valid || b_check_dist {
                    let (n_v1, n_v2) = new_pb.indices();
                    if n_v1 == n_v2 {
                        // Same vertex on both bounds — nothing to unify.
                        continue;
                    }
                    // Decide whether the vertices interfere; when they do, the
                    // block has no valid range and they must be unified.
                    let mut unify = !has_valid;
                    if b_check_dist {
                        if let (Some(a), Some(b)) = (ds.shape(n_v1).cloned(), ds.shape(n_v2).cloned()) {
                            let pb = BRepTool::vertex_point(&Vertex(b));
                            if AlgoTools::compute_vv(&a, &pb, ctx.fuzzy) == 1 {
                                unify = true;
                            }
                        }
                    }
                    if unify {
                        let key = (n_v1.min(n_v2), n_v1.max(n_v2));
                        if unify_fence.insert(key) {
                            intersect_make_sd_vertices(ds, &[n_v1, n_v2], true)?;
                        }
                        continue;
                    }
                }
                new_list.push(new_pb);
            }
        }
        *ds.change_pave_blocks_mut(n_e) = new_list;
    }
    Ok(())
}
