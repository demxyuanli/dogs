use super::super::prelude::*;
use super::analytical::analytical_filter;
use super::intervals::{init_params_from_intervals, to_split_intervals, undefined_interval};
use super::set::ParamSet;
use super::splitter::RangeSplitter;
use super::surface_type::{mesh_surface_type, param_resolution};

/// Filters a parameter set so no two samples are closer than `min_dist`, then
/// keeps samples separated by at least `filter_dist` (except a trailing buffer).
/// Source: `BRepMesh_NURBSRangeSplitter::filterParameters`.
pub(in crate::meshing::range_splitter) fn filter_parameters(
    params: &ParamSet,
    min_dist: f64,
    filter_dist: f64,
) -> Vec<f64> {
    let mut arr: Vec<f64> = params.iter().copied().collect();
    if arr.is_empty() {
        return Vec::new();
    }
    // Mandatory pre-filtering using the minimal distance.
    let mut len = 1;
    for j in 1..arr.len() {
        if arr[j] - arr[len - 1] > min_dist {
            if len < j {
                arr[len] = arr[j];
            }
            len += 1;
        }
    }
    // Series filtering.
    let mut result = Vec::new();
    let mut last_added = arr[0];
    let mut last_candidate = last_added;
    let mut candidate_defined = false;
    result.push(last_added);
    // `cxx:614`: `for (j = 2; j < aParamLength; j++)` -- the last
    // pre-filtered sample is not series-filtered; it is appended after
    // the loop (`cxx:637`). Including it here promotes the trailing
    // candidate and duplicates the end parameter.
    let mut j = 1;
    while j + 1 < len {
        let val = arr[j];
        if val - last_added > filter_dist {
            if candidate_defined {
                last_added = last_candidate;
                candidate_defined = false;
                // Re-examine the same sample against the promoted candidate.
            } else {
                last_added = val;
                j += 1;
            }
            result.push(last_added);
        } else {
            last_candidate = val;
            candidate_defined = true;
            j += 1;
        }
    }
    result.push(arr[len - 1]);
    result
}

/// `IMeshTools_Parameters` MinSize after `BRepMesh_IncrementalMesh::initParameters`
/// (`discret_root.rs` / `IMeshTools` RelMinSize from deflection).
fn resolved_min_size(params: &MeshParameters) -> f64 {
    if params.min_size >= CONFUSION {
        params.min_size
    } else {
        let defl_int = if params.deflection_interior >= CONFUSION {
            params.deflection_interior
        } else {
            params.deflection
        };
        (MeshParameters::rel_min_size() * params.deflection.min(defl_int)).max(CONFUSION)
    }
}

/// Computes the grain (filter distance) and applies it to the source parameters.
/// Source: `BRepMesh_NURBSRangeSplitter::computeGrainAndFilterParameters`.
pub(in crate::meshing::range_splitter) fn compute_grain_and_filter(
    s: &dyn RangeSplitter,
    source: &ParamSet,
    tol2d: f64,
    range_diff: f64,
    delta: f64,
    params: &MeshParameters,
) -> Vec<f64> {
    let surface = match s.surface() {
        Some(sf) => sf.as_ref(),
        None => return Vec::new(),
    };
    let mut min_diff = PCONFUSION;
    if delta < 1.0 {
        min_diff /= delta;
    }
    // `cxx:556-557`: `U/VResolution(theParameters.MinSize)` after initParameters.
    let min_size = resolved_min_size(params);
    let min_size_2d = param_resolution(surface, min_size, true)
        .max(param_resolution(surface, min_size, false));
    min_diff = min_diff.max(min_size_2d);
    let diff_max_lim = 0.1 * range_diff;
    let diff_min_lim = (0.005 * range_diff).max(2.0 * tol2d);
    let diff = min_size_2d.max(diff_max_lim.min(diff_min_lim));
    filter_parameters(source, min_diff, diff)
}

/// Shared NURBS interior-node generation. Runs through `&dyn RangeSplitter` so
/// the `get_undefined_interval_nb` / `init_parameters` hooks dispatch to the
/// concrete splitter (NURBS, Undefined, BoundaryParams, Extrusion).
///
/// `BRepMesh_NURBSRangeSplitter::GenerateSurfaceNodes` (`cxx:317-398`).
pub(in crate::meshing::range_splitter) fn generate_nurbs_grid(
    s: &dyn RangeSplitter,
    params: &MeshParameters,
) -> Option<Vec<GpPnt2d>> {
    if !s.init_parameters() {
        return None;
    }
    let surface = s.surface()?;
    let range_u = s.range_u();
    let range_v = s.range_v();
    let delta = s.delta();
    let deflection = s.deflection();
    // `initParameters` always requests `GeomAbs_CN` (`cxx:441`).
    const GEOM_ABS_CN: u8 = 6;

    let mut u_params: ParamSet = s.parameters_u().cloned().unwrap_or_default();
    let mut v_params: ParamSet = s.parameters_v().cloned().unwrap_or_default();
    // `BRepMesh_NURBSRangeSplitter::initParameters` (`cxx:438-477`) fills CN
    // intervals. `BoundaryParamsRangeSplitter::initParameters` returns true
    // without that fill (`cxx` override) — only AddPoint samples.
    if s.seeds_cn_intervals() {
        // `initParameters` (`cxx:465-487`): fail the whole grid when either
        // direction adds no interval sample. INTERNAL edge params are already
        // in `parameters_u/v` via `grabParamsOfEdges` (`cxx:487`); do not
        // invent range-end samples (not in cxx).
        let u_intervals = undefined_interval(s, true, GEOM_ABS_CN, range_u);
        let v_intervals = undefined_interval(s, false, GEOM_ABS_CN, range_v);
        let is_split = to_split_intervals(surface.as_ref(), &u_intervals, &v_intervals);
        if std::env::var("NURBSGRID").is_ok() {
            eprintln!(
                "NURBSGRID intervals split={is_split} u_n={} v_n={} surf_u=[{:.17e},{:.17e}] surf_v=[{:.17e},{:.17e}]",
                u_intervals.len(),
                v_intervals.len(),
                surface.as_ref().u_range().0,
                surface.as_ref().u_range().1,
                surface.as_ref().v_range().0,
                surface.as_ref().v_range().1
            );
            eprint!("NURBSGRID u_int=");
            for (i, &u) in u_intervals.iter().enumerate() {
                if i > 0 {
                    eprint!(",");
                }
                eprint!("{u:.17e}");
            }
            eprintln!();
            eprint!("NURBSGRID v_int=");
            for (i, &v) in v_intervals.iter().enumerate() {
                if i > 0 {
                    eprint!(",");
                }
                eprint!("{v:.17e}");
            }
            eprintln!();
        }
        let mut u_iv = ParamSet::new();
        let mut v_iv = ParamSet::new();
        let dbg_none = std::env::var("NURBSGRID").is_ok();
        let fkey_dbg = s
            .base()
            .dface
            .as_ref()
            .map(|df| crate::tgeometry::GeometryRegistry::shape_key(&df.face().0))
            .unwrap_or(0);
        if !init_params_from_intervals(&u_intervals, range_u, is_split, &mut u_iv) {
            if dbg_none {
                eprintln!(
                    "NURBSGRID none key={fkey_dbg} why=u iv_n={} range=[{:.17e},{:.17e}]",
                    u_intervals.len(),
                    range_u.0,
                    range_u.1
                );
            }
            return None;
        }
        if !init_params_from_intervals(&v_intervals, range_v, is_split, &mut v_iv) {
            if dbg_none {
                eprintln!(
                    "NURBSGRID none key={fkey_dbg} why=v iv_n={} range=[{:.17e},{:.17e}]",
                    v_intervals.len(),
                    range_v.0,
                    range_v.1
                );
                for (i, &v) in v_intervals.iter().enumerate() {
                    eprint!("{}[{:.17e}]", if i == 0 { "" } else { "," }, v);
                }
                eprintln!();
            }
            return None;
        }
        for &p in u_iv.iter() {
            u_params.insert(p);
        }
        for &p in v_iv.iter() {
            v_params.insert(p);
        }
    }

    let tol_u = param_resolution(surface.as_ref(), deflection, true);
    let tol_v = param_resolution(surface.as_ref(), deflection, false);
    if std::env::var("NURBSGRID").is_ok() {
        // TEMP: tag each dump with the discrete face's shape key so a dump can be
        // paired with `FSTAT face=<i> key=<k>` (and with the probe's per-face rows).
        let fkey = s
            .base()
            .dface
            .as_ref()
            .map(|df| crate::tgeometry::GeometryRegistry::shape_key(&df.face().0))
            .unwrap_or(0);
        eprintln!(
            "NURBSGRID pre key={fkey} type={:?} u_params={} v_params={} ru=[{:.17e},{:.17e}] rv=[{:.17e},{:.17e}] tol_u={tol_u:.6e} tol_v={tol_v:.6e} defl={deflection:.6} delta=({:.17e},{:.17e})",
            mesh_surface_type(surface.as_ref()),
            u_params.len(),
            v_params.len(),
            range_u.0,
            range_u.1,
            range_v.0,
            range_v.1,
            delta.0,
            delta.1
        );
        let mut us: Vec<f64> = u_params.iter().copied().collect();
        us.sort_by(|a, b| a.partial_cmp(b).unwrap());
        eprint!("NURBSGRID u_params_sorted=");
        for (i, &u) in us.iter().enumerate() {
            if i > 0 {
                eprint!(",");
            }
            eprint!("{u:?}");
        }
        eprintln!();
        let mut vs: Vec<f64> = v_params.iter().copied().collect();
        vs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        eprint!("NURBSGRID v_params_sorted=");
        for (i, &v) in vs.iter().enumerate() {
            if i > 0 {
                eprint!(",");
            }
            eprint!("{v:.17e}");
        }
        eprintln!();
    }
    let mut u_seq =
        compute_grain_and_filter(s, &u_params, tol_u, range_u.1 - range_u.0, delta.0, params);
    let mut v_seq =
        compute_grain_and_filter(s, &v_params, tol_v, range_v.1 - range_v.0, delta.1, params);
    if std::env::var("NURBSGRID").is_ok() {
        eprintln!(
            "NURBSGRID grain u_seq={} v_seq={}",
            u_seq.len(),
            v_seq.len()
        );
        eprint!("NURBSGRID u_grain=");
        for (i, &u) in u_seq.iter().enumerate() {
            if i > 0 {
                eprint!(",");
            }
            eprint!("{u:.6}");
        }
        eprintln!();
        eprint!("NURBSGRID v_grain=");
        for (i, &v) in v_seq.iter().enumerate() {
            if i > 0 {
                eprint!(",");
            }
            eprint!("{v:.17e}");
        }
        eprintln!();
        if let Some(axis) = surface.as_ref().revolution_axis() {
            let loc = axis.location();
            let dir = axis.direction();
            eprint!("NURBSGRID v_rad=");
            for (i, &v) in v_seq.iter().enumerate() {
                let p = surface.as_ref().d0(0.0, v);
                let rel = occt_core::gp::GpVec::from_pnts(&loc, &p);
                let along = dir.xyz().multiplied(rel.xyz().dot(dir.xyz()));
                let rad = rel.xyz().subtracted(&along).modulus();
                if i > 0 {
                    eprint!(",");
                }
                eprint!("{rad:.6}");
            }
            eprintln!();
        }
    }

    // `AngleInterior` after `initParameters` (`2*Angle` when below Angular).
    let angle_interior = if params.angle_interior >= ANGULAR {
        params.angle_interior
    } else {
        2.0 * params.angle
    };
    let min_size = resolved_min_size(params);
    let mut u_remove: std::collections::HashSet<u64> = std::collections::HashSet::new();
    let mut v_remove: std::collections::HashSet<u64> = std::collections::HashSet::new();
    let mut u_fixed: std::collections::HashSet<u64> = std::collections::HashSet::new();
    let mut v_fixed: std::collections::HashSet<u64> = std::collections::HashSet::new();
    analytical_filter(
        surface.as_ref(),
        false,
        &v_seq,
        &mut u_seq,
        deflection,
        angle_interior,
        min_size,
        &mut u_remove,
        &mut v_fixed,
        &mut u_fixed,
    );
    analytical_filter(
        surface.as_ref(),
        true,
        &u_seq,
        &mut v_seq,
        deflection,
        angle_interior,
        min_size,
        &mut v_remove,
        &mut u_fixed,
        &mut v_fixed,
    );
    for k in &u_fixed {
        u_remove.remove(k);
    }
    for k in &v_fixed {
        v_remove.remove(k);
    }
    if std::env::var("NOAFRM").is_ok() {
        u_remove.clear();
        v_remove.clear();
    }

    let mut nodes = Vec::with_capacity(u_seq.len() * v_seq.len());
    for &u in &u_seq {
        if u_remove.contains(&u.to_bits()) {
            continue;
        }
        for &v in &v_seq {
            if v_remove.contains(&v.to_bits()) {
                continue;
            }
            nodes.push(GpPnt2d::new(u, v));
        }
    }
    if std::env::var("NURBSGRID").is_ok() {
        let nu = u_seq.iter().filter(|u| !u_remove.contains(&u.to_bits())).count();
        let nv = v_seq.iter().filter(|v| !v_remove.contains(&v.to_bits())).count();
        eprintln!(
            "NURBSGRID post u_kept={nu}/{} v_kept={nv}/{} nodes={} u_rm={} v_rm={} ang_int={angle_interior:.6}",
            u_seq.len(),
            v_seq.len(),
            nodes.len(),
            u_remove.len(),
            v_remove.len()
        );
        eprint!("NURBSGRID u_after_insert=");
        for (i, &u) in u_seq.iter().enumerate() {
            if i > 0 {
                eprint!(",");
            }
            eprint!("{u:.6}");
        }
        eprintln!();
        eprint!("NURBSGRID u_removed=");
        let mut rms: Vec<f64> = u_remove.iter().map(|b| f64::from_bits(*b)).collect();
        rms.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for (i, &u) in rms.iter().enumerate() {
            if i > 0 {
                eprint!(",");
            }
            eprint!("{u:.6}");
        }
        eprintln!();
        eprint!("NURBSGRID u_kept_vals=");
        let mut first = true;
        for &u in &u_seq {
            if u_remove.contains(&u.to_bits()) {
                continue;
            }
            if !first {
                eprint!(",");
            }
            first = false;
            eprint!("{u:.6}");
        }
        eprintln!();
    }
    Some(nodes)
}
