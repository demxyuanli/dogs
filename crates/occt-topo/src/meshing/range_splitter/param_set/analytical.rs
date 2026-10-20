use super::super::prelude::*;

/// Squared distance from `mid` to the infinite line through `p1..p2`.
/// `BRepMesh_GeomTool::SquareDeflectionOfSegment` (`hxx:181-193`) uses
/// `gp_Lin::SquareDistance`, not a clamped chord.
pub(in crate::meshing::range_splitter) fn sq_deflection_of_segment(
    p1: &GpPnt,
    p2: &GpPnt,
    mid: &GpPnt,
) -> f64 {
    crate::meshing::geom_tool::GeomTool::square_deflection_of_segment(p1, p2, mid)
}

/// `BRepMesh_NURBSRangeSplitter::AnalyticalFilter` (`cxx:30-215`): insert
/// midpoints that miss deflection / angle, and mark too-dense control params
/// for removal after both iso passes.
#[allow(clippy::too_many_arguments)]
pub(in crate::meshing::range_splitter) fn analytical_filter(
    surface: &dyn Surface,
    is_iso_u: bool,
    iso_params: &[f64],
    control_params: &mut Vec<f64>,
    deflection: f64,
    angle_interior: f64,
    min_size: f64,
    control_remove: &mut std::collections::HashSet<u64>,
    iso_forbidden: &mut std::collections::HashSet<u64>,
    control_forbidden: &mut std::collections::HashSet<u64>,
) {
    if control_params.len() < 2 || iso_params.is_empty() {
        return;
    }
    let sq_max_deflection = deflection * deflection;
    let sq_min_size = min_size * min_size;

    // OCCT: IsoU scans every iso-line; IsoV skips the two outer ones.
    let (start, end) = if is_iso_u {
        (0, iso_params.len())
    } else {
        (1, iso_params.len().saturating_sub(1))
    };

    // `cxx:74-75` builds `GeomAdaptor_Curve(UIso/VIso)` then `D1`.
    // Offset UIso/VIso is AdvApprox (`Geom_OffsetSurface.cxx:601-687`).
    for &iso_param in &iso_params[start..end] {
        let iso_curve = if is_iso_u {
            surface.u_iso_curve(iso_param)
        } else {
            surface.v_iso_curve(iso_param)
        };
        let iso_point = |t: f64| -> GpPnt {
            if let Some(ref c) = iso_curve {
                c.d0(t)
            } else if is_iso_u {
                surface.d0(iso_param, t)
            } else {
                surface.d0(t, iso_param)
            }
        };
        let iso_tangent = |t: f64| -> GpVec {
            if let Some(ref c) = iso_curve {
                c.d1(t).1
            } else {
                let (_, du, dv) = if is_iso_u {
                    surface.d1(iso_param, t)
                } else {
                    surface.d1(t, iso_param)
                };
                if is_iso_u {
                    dv
                } else {
                    du
                }
            }
        };

        let mut prev_param = control_params[0];
        let mut prev_pnt = iso_point(prev_param);
        let mut prev_vec = iso_tangent(prev_param);

        let mut j = 1usize;
        while j < control_params.len() {
            let curr_param = control_params[j];
            let curr_pnt = iso_point(curr_param);
            let curr_vec = iso_tangent(curr_param);

            let mid_param = 0.5 * (prev_param + curr_param);
            let mid_pnt = iso_point(mid_param);
            let sq_dist = sq_deflection_of_segment(&prev_pnt, &curr_pnt, &mid_pnt);
            let angle = if prev_vec.square_magnitude() > SQUARE_CONFUSION
                && curr_vec.square_magnitude() > SQUARE_CONFUSION
            {
                prev_vec.angle(&curr_vec).abs()
            } else {
                0.0
            };

            if (sq_dist > sq_max_deflection || angle > angle_interior) && sq_dist > sq_min_size {
                control_params.insert(j, mid_param);
                continue;
            }

            if (sq_dist < sq_max_deflection || angle < angle_interior)
                && control_params.len() > 3
                && j + 1 < control_params.len()
            {
                let next_param = control_params[j + 1];
                let next_pnt = iso_point(next_param);
                let next_vec = iso_tangent(next_param);
                let look_mid = iso_point(0.5 * (prev_param + next_param));
                let look_dist = sq_deflection_of_segment(&prev_pnt, &next_pnt, &look_mid);
                if look_dist < sq_max_deflection {
                    // `BRepMesh_NURBSRangeSplitter.cxx:165-167`: both
                    // `SquareMagnitude()` guards use `gp::Resolution()` =
                    // `RealSmall()` = `DBL_MIN` (`gp.hxx:60`).
                    let look_ok = prev_vec.square_magnitude() < REAL_SMALL
                        || next_vec.square_magnitude() < REAL_SMALL
                        || prev_vec.angle(&next_vec).abs() < angle_interior;
                    if look_ok {
                        control_remove.insert(curr_param.to_bits());
                        prev_param = next_param;
                        prev_pnt = next_pnt;
                        prev_vec = next_vec;
                        j += 2;
                        continue;
                    } else {
                        iso_forbidden.insert(iso_param.to_bits());
                        control_forbidden.insert(curr_param.to_bits());
                    }
                }
            }

            prev_param = curr_param;
            prev_pnt = curr_pnt;
            prev_vec = curr_vec;
            j += 1;
        }
    }
}
