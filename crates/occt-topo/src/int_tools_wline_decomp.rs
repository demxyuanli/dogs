//! `IntTools_WLineTool::DecompositionOfWLine`.
//!
//! Source: `IntTools_WLineTool.cxx:492`. Splits a walking line where samples
//! cross a periodic surface boundary so each resulting WLine stays in one
//! period. When only one piece is found the function returns `false` and
//! leaves `new_lines` empty (OCCT then uses the original WLine).

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_core::precision::{CONFUSION, PCONFUSION, RESOLUTION};
use occt_geom::Surface;
use occt_geom::geom_api::project_point_on_surface;

use crate::int_tools_wline::{
    adjust_by_neighbour, wline_adjust_periodic, find_point, is_point_on_boundary, u_period,
    u_resolution, v_period, v_resolution, PntOn2S, WLine,
};
use crate::shape::Face;

/// `IntTools_WLineTool::DecompositionOfWLine`.
///
/// `line_parts` are 1-based `(ifprm, ilprm)` ranges from `GeomInt_LineConstructor`.
/// When `avoid_line_constructor` is true the whole WLine is one part.
pub fn decomposition_of_wline(
    wline: &WLine,
    surface1: &dyn Surface,
    surface2: &dyn Surface,
    _the_face1: &Face,
    _the_face2: &Face,
    line_parts: &[(i32, i32)],
    avoid_l_constructor: bool,
    tol: f64,
    new_lines: &mut Vec<WLine>,
) -> bool {
    let b_ret = false;
    let nb_pnts_val = wline.nb_pnts();
    let b_avoid_line_constructor = avoid_l_constructor;
    if nb_pnts_val == 0 {
        return b_ret;
    }
    if !b_avoid_line_constructor && line_parts.is_empty() {
        return b_ret;
    }

    let mut an_array_of_lines: Vec<Vec<i32>> = vec![Vec::new(); (nb_pnts_val as usize) + 1];
    let mut an_array_of_line_type: Vec<i32> = vec![0; (nb_pnts_val as usize) + 1];
    let mut list_of_point_index: Vec<i32> = Vec::new();

    let mut nblines = 0i32;
    let tol_val = 0.5 * CONFUSION;
    let mut b_is_prev_point_on_boundary = false;
    let mut b_is_point_on_boundary = false;

    for pit in 1..=nb_pnts_val {
        let mut b_is_current_point_on_boundary = false;
        let point_val = wline.point(pit);
        for i in 0..2 {
            let ga = if i == 0 { surface1 } else { surface2 };
            let (umin, umax) = ga.u_range();
            let (vmin, vmax) = ga.v_range();
            let (u, v) = if i == 0 {
                point_val.parameters_on_s1()
            } else {
                point_val.parameters_on_s2()
            };
            for j in 0..2 {
                let isperiodic = if j == 0 {
                    ga.is_u_periodic()
                } else {
                    ga.is_v_periodic()
                };
                if !isperiodic {
                    continue;
                }
                let (resolution, period, alowerboundary, aupperboundary, parameter) =
                    if j == 0 {
                        (
                            u_resolution(ga, tol_val),
                            u_period(ga).unwrap_or(0.0),
                            umin,
                            umax,
                            u,
                        )
                    } else {
                        (
                            v_resolution(ga, tol_val),
                            v_period(ga).unwrap_or(0.0),
                            vmin,
                            vmax,
                            v,
                        )
                    };
                let (an_adjust_par, _anoffset) =
                    wline_adjust_periodic(parameter, alowerboundary, aupperboundary, period, 0.0);
                let mut b_is_on_first_boundary = true;
                b_is_point_on_boundary = is_point_on_boundary(
                    an_adjust_par,
                    alowerboundary,
                    aupperboundary,
                    resolution,
                    &mut b_is_on_first_boundary,
                );
                if b_is_point_on_boundary {
                    b_is_current_point_on_boundary = true;
                    break;
                }
            }
            if b_is_current_point_on_boundary {
                break;
            }
        }
        if b_is_current_point_on_boundary != b_is_prev_point_on_boundary {
            if !list_of_point_index.is_empty() {
                nblines += 1;
                an_array_of_lines[nblines as usize] = list_of_point_index.clone();
                an_array_of_line_type[nblines as usize] = i32::from(b_is_prev_point_on_boundary);
                list_of_point_index.clear();
            }
            b_is_prev_point_on_boundary = b_is_current_point_on_boundary;
        }
        list_of_point_index.push(pit);
    }
    if !list_of_point_index.is_empty() {
        nblines += 1;
        an_array_of_lines[nblines as usize] = list_of_point_index.clone();
        an_array_of_line_type[nblines as usize] = i32::from(b_is_prev_point_on_boundary);
        list_of_point_index.clear();
    }
    if nblines <= 1 {
        return b_ret;
    }

    // 2. Correct wlines.begin
    let mut an_array_of_line_ends: Vec<Vec<i32>> = vec![Vec::new(); (nblines as usize) + 1];
    let mut seq_of_pnt_on_2s: Vec<PntOn2S> = Vec::new();

    for i in 1..=nblines {
        if an_array_of_line_type[i as usize] != 0 {
            continue;
        }
        let list_of_index = &an_array_of_lines[i as usize];
        let mut list_of_fl_index: Vec<i32> = Vec::new();

        for j in 0..2 {
            let aneighbourindex = if j == 0 { i - 1 } else { i + 1 };
            if aneighbourindex < 1 || aneighbourindex > nblines {
                continue;
            }
            if an_array_of_line_type[aneighbourindex as usize] == 0 {
                continue;
            }
            let neighbour = &an_array_of_lines[aneighbourindex as usize];
            if neighbour.is_empty() {
                continue;
            }
            let an_index = if j == 0 {
                *neighbour.last().unwrap()
            } else {
                neighbour[0]
            };
            let point_val = *wline.point(an_index);
            let mut new_p = point_val;
            if list_of_index.len() < 2 {
                seq_of_pnt_on_2s.push(new_p);
                list_of_fl_index.push(seq_of_pnt_on_2s.len() as i32);
                continue;
            }
            let i_first = list_of_index[0];
            let i_last = *list_of_index.last().unwrap();

            for surfit in 0..2 {
                let ga = if surfit == 0 {
                    surface1
                } else {
                    surface2
                };
                let (umin, umax) = ga.u_range();
                let (vmin, vmax) = ga.v_range();
                let (u, v) = if surfit == 0 {
                    new_p.parameters_on_s1()
                } else {
                    new_p.parameters_on_s2()
                };
                let mut nbboundaries = 0i32;
                let mut b_is_near_boundary = false;
                let mut b_is_u_boundary = false;
                let mut b_is_first_boundary = false;

                for parit in 0..2 {
                    let isperiodic = if parit == 0 {
                        ga.is_u_periodic()
                    } else {
                        ga.is_v_periodic()
                    };
                    let resolution = if parit == 0 {
                        u_resolution(ga, tol_val)
                    } else {
                        v_resolution(ga, tol_val)
                    };
                    let alowerboundary = if parit == 0 { umin } else { vmin };
                    let aupperboundary = if parit == 0 { umax } else { vmax };
                    let parameter = if parit == 0 { u } else { v };
                    let mut b_is_on_first_boundary = true;

                    if !isperiodic {
                        b_is_point_on_boundary = is_point_on_boundary(
                            parameter,
                            alowerboundary,
                            aupperboundary,
                            resolution,
                            &mut b_is_on_first_boundary,
                        );
                        if b_is_point_on_boundary {
                            b_is_u_boundary = parit == 0;
                            b_is_first_boundary = b_is_on_first_boundary;
                            nbboundaries += 1;
                        }
                    } else {
                        let period = if parit == 0 {
                            u_period(ga).unwrap_or(0.0)
                        } else {
                            v_period(ga).unwrap_or(0.0)
                        };
                        let (an_adjust_par, _anoffset) =
                            wline_adjust_periodic(parameter, alowerboundary, aupperboundary, period, 0.0);
                        b_is_point_on_boundary = is_point_on_boundary(
                            an_adjust_par,
                            alowerboundary,
                            aupperboundary,
                            resolution,
                            &mut b_is_on_first_boundary,
                        );
                        if b_is_point_on_boundary {
                            b_is_u_boundary = parit == 0;
                            b_is_first_boundary = b_is_on_first_boundary;
                            nbboundaries += 1;
                        } else {
                            let mut an_epsilon = resolution * 100.0;
                            let part = (aupperboundary - alowerboundary) * 0.1;
                            an_epsilon = if an_epsilon > part { part } else { an_epsilon };
                            b_is_near_boundary = is_point_on_boundary(
                                an_adjust_par,
                                alowerboundary,
                                aupperboundary,
                                an_epsilon,
                                &mut b_is_on_first_boundary,
                            );
                        }
                    }
                }

                let mut b_compute_line_end = false;
                if nbboundaries == 2 {
                    b_compute_line_end = true;
                } else if nbboundaries == 1 {
                    let isperiodic = if b_is_u_boundary {
                        ga.is_u_periodic()
                    } else {
                        ga.is_v_periodic()
                    };
                    if isperiodic {
                        let alowerboundary = if b_is_u_boundary { umin } else { vmin };
                        let aupperboundary = if b_is_u_boundary { umax } else { vmax };
                        let period = if b_is_u_boundary {
                            u_period(ga).unwrap_or(0.0)
                        } else {
                            v_period(ga).unwrap_or(0.0)
                        };
                        let parameter = if b_is_u_boundary { u } else { v };
                        let (an_adjust_par, anoffset) =
                            wline_adjust_periodic(parameter, alowerboundary, aupperboundary, period, 0.0);
                        let adist = if b_is_first_boundary {
                            (an_adjust_par - alowerboundary).abs()
                        } else {
                            (an_adjust_par - aupperboundary).abs()
                        };
                        let mut another_par = if b_is_first_boundary {
                            aupperboundary - adist
                        } else {
                            alowerboundary + adist
                        };
                        another_par += anoffset;
                        let aneighbourpointindex = if j == 0 { i_first } else { i_last };
                        let neighbour_point = wline.point(aneighbourpointindex);
                        let (n_u1, n_v1) = if surfit == 0 {
                            neighbour_point.parameters_on_s1()
                        } else {
                            neighbour_point.parameters_on_s2()
                        };
                        let adist1 = if b_is_u_boundary {
                            (n_u1 - u).abs()
                        } else {
                            (n_v1 - v).abs()
                        };
                        let adist2 = if b_is_u_boundary {
                            (n_u1 - another_par).abs()
                        } else {
                            (n_v1 - another_par).abs()
                        };
                        b_compute_line_end = true;
                        let mut b_check_angle1 = false;
                        let mut b_check_angle2 = false;
                        let mut new_vec = GpVec2d::zero();
                        let anew_u = if b_is_u_boundary { another_par } else { u };
                        let anew_v = if b_is_u_boundary { v } else { another_par };

                        if (adist1 - adist2) > PCONFUSION && adist2 < (period / 4.0) {
                            b_check_angle1 = true;
                            new_vec = GpVec2d::new(anew_u - n_u1, anew_v - n_v1);
                            if new_vec.square_magnitude() < RESOLUTION {
                                new_p.set_uv_on(surfit == 0, anew_u, anew_v);
                                b_check_angle1 = false;
                            }
                        } else if adist1 < (period / 4.0) {
                            b_check_angle2 = true;
                            new_vec = GpVec2d::new(u - n_u1, v - n_v1);
                            if new_vec.square_magnitude() < RESOLUTION {
                                b_check_angle2 = false;
                            }
                        }

                        if b_check_angle1 || b_check_angle2 {
                            let mut anindexother = aneighbourpointindex;
                            while anindexother <= i_last && anindexother >= i_first {
                                anindexother = if j == 0 {
                                    anindexother + 1
                                } else {
                                    anindexother - 1
                                };
                                if anindexother > i_last || anindexother < i_first {
                                    break;
                                }
                                let prev = wline.point(anindexother);
                                let (n_u2, n_v2) = if surfit == 0 {
                                    prev.parameters_on_s1()
                                } else {
                                    prev.parameters_on_s2()
                                };
                                let vec_old = GpVec2d::new(n_u1 - n_u2, n_v1 - n_v2);
                                if vec_old.square_magnitude() <= RESOLUTION {
                                    continue;
                                }
                                let an_angle = new_vec.angle(&vec_old);
                                if an_angle.abs() < (std::f64::consts::PI * 0.25)
                                    && new_vec.dot(&vec_old) > 0.0
                                {
                                    if b_check_angle1 {
                                        let mut atmppoint = new_p;
                                        atmppoint.set_uv_on(surfit == 0, anew_u, anew_v);
                                        let (u1, v1, u2, v2) = atmppoint.parameters();
                                        let p1 = surface1.d0(u1, v1);
                                        let p2 = surface2.d0(u2, v2);
                                        let p0 = point_val.value();
                                        if p0.distance(&p1) <= tol_val
                                            && p0.distance(&p2) <= tol_val
                                            && p1.distance(&p2) <= tol_val
                                        {
                                            b_compute_line_end = false;
                                            new_p.set_uv_on(surfit == 0, anew_u, anew_v);
                                        }
                                    }
                                    if b_check_angle2 {
                                        b_compute_line_end = false;
                                    }
                                }
                                break;
                            }
                        }
                    }
                } else if b_is_near_boundary {
                    b_compute_line_end = true;
                }

                if b_compute_line_end {
                    let mut anewpoint = GpPnt2d::zero();
                    let mut found = false;
                    if b_is_near_boundary {
                        let (u1, v1, u2, v2) = new_p.parameters();
                        anewpoint = if surfit == 0 {
                            GpPnt2d::new(u1, v1)
                        } else {
                            GpPnt2d::new(u2, v2)
                        };
                        let aneighbourpointindex1 = if j == 0 { i_first } else { i_last };
                        let neighbour_point = wline.point(aneighbourpointindex1);
                        let (n_u1, n_v1) = if surfit == 0 {
                            neighbour_point.parameters_on_s1()
                        } else {
                            neighbour_point.parameters_on_s2()
                        };
                        let ap1 = GpPnt2d::new(n_u1, n_v1);
                        if ga.is_u_periodic() || ga.is_v_periodic() {
                            let ap2 = adjust_by_neighbour(&ap1, &anewpoint, ga);
                            if ap2.x() < umin || ap2.x() > umax || ap2.y() < vmin || ap2.y() > vmax
                            {
                                if let Some(p) = find_point(&ap1, &ap2, umin, umax, vmin, vmax) {
                                    anewpoint = p;
                                    found = true;
                                }
                            } else {
                                anewpoint = ap2;
                                new_p.set_uv_on(surfit == 0, anewpoint.x(), anewpoint.y());
                            }
                        }
                    } else {
                        let aneighbourpointindex1 = if j == 0 { i_first } else { i_last };
                        let neighbour_point = wline.point(aneighbourpointindex1);
                        let (n_u1, n_v1) = if surfit == 0 {
                            neighbour_point.parameters_on_s1()
                        } else {
                            neighbour_point.parameters_on_s2()
                        };
                        let ap1 = GpPnt2d::new(n_u1, n_v1);
                        let mut ap2 = GpPnt2d::new(n_u1, n_v1);
                        let mut aneighbourpointindex2 = aneighbourpointindex1;
                        while aneighbourpointindex2 <= i_last && aneighbourpointindex2 >= i_first {
                            aneighbourpointindex2 = if j == 0 {
                                aneighbourpointindex2 + 1
                            } else {
                                aneighbourpointindex2 - 1
                            };
                            if aneighbourpointindex2 > i_last || aneighbourpointindex2 < i_first {
                                break;
                            }
                            let prev = wline.point(aneighbourpointindex2);
                            let (n_u2, n_v2) = if surfit == 0 {
                                prev.parameters_on_s1()
                            } else {
                                prev.parameters_on_s2()
                            };
                            ap2 = GpPnt2d::new(n_u2, n_v2);
                            if ap1.square_distance(&ap2) > RESOLUTION {
                                break;
                            }
                        }
                        if let Some(p) = find_point(&ap2, &ap1, umin, umax, vmin, vmax) {
                            anewpoint = p;
                            found = true;
                        }
                    }

                    if found {
                        let criteria = tol;
                        let surface = if surfit == 0 {
                            surface1
                        } else {
                            surface2
                        };
                        let surface_other = if surfit == 0 {
                            surface2
                        } else {
                            surface1
                        };
                        let p3d_val = surface.d0(anewpoint.x(), anewpoint.y());
                        if let Some(proj) = project_point_on_surface(surface_other, &p3d_val, 0.0) {
                            if proj.distance < criteria {
                                let mut found_u = u;
                                let mut found_v = v;
                                found_u = proj.u;
                                found_v = proj.v;
                                let aneindex1 = if j == 0 { i_first } else { i_last };
                                let neighbour_point = wline.point(aneindex1);
                                let (n_un, n_vn) = if surfit == 0 {
                                    neighbour_point.parameters_on_s2()
                                } else {
                                    neighbour_point.parameters_on_s1()
                                };
                                let neighbour_2d = GpPnt2d::new(n_un, n_vn);
                                let an_adjusted_point = adjust_by_neighbour(
                                    &neighbour_2d,
                                    &GpPnt2d::new(found_u, found_v),
                                    surface_other,
                                );
                                found_u = an_adjusted_point.x();
                                found_v = an_adjusted_point.y();
                                // OCCT writes `X < umin && X > umax && Y < vmin && Y > vmax`
                                // (always false). Copied as written.
                                if (an_adjusted_point.x() < umin)
                                    && (an_adjusted_point.x() > umax)
                                    && (an_adjusted_point.y() < vmin)
                                    && (an_adjusted_point.y() > vmax)
                                {
                                    found_u = found_u.max(umin).min(umax);
                                    found_v = found_v.max(vmin).min(vmax);
                                    let p3d = surface_other.d0(found_u, found_v);
                                    if let Some(proj2) =
                                        project_point_on_surface(surface, &p3d, 0.0)
                                    {
                                        if proj2.distance < criteria {
                                            anewpoint = GpPnt2d::new(proj2.u, proj2.v);
                                        }
                                    }
                                }
                                if surfit == 0 {
                                    new_p.set_value(
                                        p3d_val,
                                        anewpoint.x(),
                                        anewpoint.y(),
                                        found_u,
                                        found_v,
                                    );
                                } else {
                                    new_p.set_value(
                                        p3d_val,
                                        found_u,
                                        found_v,
                                        anewpoint.x(),
                                        anewpoint.y(),
                                    );
                                }
                            }
                        }
                    }
                }
            }
            seq_of_pnt_on_2s.push(new_p);
            list_of_fl_index.push(seq_of_pnt_on_2s.len() as i32);
        }
        an_array_of_line_ends[i as usize] = list_of_fl_index;
    }

    // Split wlines.begin
    let nbiter = if b_avoid_line_constructor {
        1
    } else {
        line_parts.len().max(1)
    };
    for j in 1..=nbiter {
        let (ifprm, ilprm) = if b_avoid_line_constructor {
            (1, wline.nb_pnts())
        } else {
            line_parts[(j - 1) as usize]
        };
        let mut line_on_2s = WLine::new();
        for i in 1..=nblines {
            if an_array_of_line_type[i as usize] != 0 {
                continue;
            }
            let list_of_index = &an_array_of_lines[i as usize];
            let list_of_fl_index = &an_array_of_line_ends[i as usize];
            let mut bhasfirstpoint = list_of_fl_index.len() == 2;
            let mut bhaslastpoint = list_of_fl_index.len() == 2;
            if !bhasfirstpoint && !list_of_fl_index.is_empty() {
                bhasfirstpoint = i != 1;
            }
            if !bhaslastpoint && !list_of_fl_index.is_empty() {
                bhaslastpoint = i != nblines;
            }
            if list_of_index.is_empty() {
                continue;
            }
            let i_first = list_of_index[0];
            let i_last = *list_of_index.last().unwrap();
            let b_is_first_inside = ifprm >= i_first && ifprm <= i_last;
            let b_is_last_inside = ilprm >= i_first && ilprm <= i_last;

            if !b_is_first_inside && !b_is_last_inside {
                if ifprm < i_first && ilprm > i_last {
                    if bhasfirstpoint {
                        let pit = list_of_fl_index[0];
                        line_on_2s.add(seq_of_pnt_on_2s[(pit as usize).saturating_sub(1)]);
                    }
                    for &pit in list_of_index {
                        line_on_2s.add(*wline.point(pit));
                    }
                    if bhaslastpoint {
                        let pit = *list_of_fl_index.last().unwrap();
                        line_on_2s.add(seq_of_pnt_on_2s[(pit as usize).saturating_sub(1)]);
                    }
                    let aneighbour = i + 1;
                    let mut b_is_end_of_line = true;
                    if aneighbour <= nblines {
                        let list_of_neighbour = &an_array_of_lines[aneighbour as usize];
                        if an_array_of_line_type[aneighbour as usize] != 0
                            && list_of_neighbour.is_empty()
                        {
                            b_is_end_of_line = false;
                        }
                    }
                    if b_is_end_of_line {
                        if line_on_2s.nb_pnts() > 1 {
                            new_lines.push(std::mem::take(&mut line_on_2s));
                        }
                        line_on_2s = WLine::new();
                    }
                }
                continue;
            }
            if b_is_first_inside && b_is_last_inside {
                for &pit in list_of_index {
                    if pit < ifprm || pit > ilprm {
                        continue;
                    }
                    line_on_2s.add(*wline.point(pit));
                }
            } else {
                if b_is_first_inside {
                    for &pit in list_of_index {
                        if pit < ifprm {
                            continue;
                        }
                        line_on_2s.add(*wline.point(pit));
                    }
                    if bhaslastpoint {
                        let pit = *list_of_fl_index.last().unwrap();
                        line_on_2s.add(seq_of_pnt_on_2s[(pit as usize).saturating_sub(1)]);
                    }
                    let aneighbour = i + 1;
                    let mut b_is_end_of_line = true;
                    if aneighbour <= nblines {
                        let list_of_neighbour = &an_array_of_lines[aneighbour as usize];
                        if an_array_of_line_type[aneighbour as usize] != 0
                            && list_of_neighbour.is_empty()
                        {
                            b_is_end_of_line = false;
                        }
                    }
                    if b_is_end_of_line {
                        if line_on_2s.nb_pnts() > 1 {
                            new_lines.push(std::mem::take(&mut line_on_2s));
                        }
                        line_on_2s = WLine::new();
                    }
                }
                if b_is_last_inside {
                    if bhasfirstpoint {
                        let pit = list_of_fl_index[0];
                        line_on_2s.add(seq_of_pnt_on_2s[(pit as usize).saturating_sub(1)]);
                    }
                    for &pit in list_of_index {
                        if pit > ilprm {
                            continue;
                        }
                        line_on_2s.add(*wline.point(pit));
                    }
                }
            }
        }
        if line_on_2s.nb_pnts() > 1 {
            new_lines.push(line_on_2s);
        }
    }
    true
}
