//! `IntTools_WLineTool::DecompositionOfWLine`.
//!
//! Source: `IntTools_WLineTool.cxx:492`. Splits a walking line where samples
//! cross a periodic surface boundary so each resulting WLine stays in one
//! period. When only one piece is found the function returns `false` and
//! leaves `the_new_lines` empty (OCCT then uses the original WLine).

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_core::precision::{CONFUSION, PCONFUSION, RESOLUTION};
use occt_geom::Surface;
use occt_geom::geom_api::project_point_on_surface;

use crate::int_tools_wline::{
    adjust_by_neighbour, adjust_periodic, find_point, is_point_on_boundary, u_period,
    u_resolution, v_period, v_resolution, PntOn2S, WLine,
};
use crate::shape::Face;

/// `IntTools_WLineTool::DecompositionOfWLine`.
///
/// `line_parts` are 1-based `(ifprm, ilprm)` ranges from `GeomInt_LineConstructor`.
/// When `avoid_line_constructor` is true the whole WLine is one part.
pub fn decomposition_of_wline(
    the_wline: &WLine,
    the_surface1: &dyn Surface,
    the_surface2: &dyn Surface,
    _the_face1: &Face,
    _the_face2: &Face,
    line_parts: &[(i32, i32)],
    the_avoid_l_constructor: bool,
    the_tol: f64,
    the_new_lines: &mut Vec<WLine>,
) -> bool {
    let b_ret = false;
    let a_nb_pnts = the_wline.nb_pnts();
    let b_avoid_line_constructor = the_avoid_l_constructor;
    if a_nb_pnts == 0 {
        return b_ret;
    }
    if !b_avoid_line_constructor && line_parts.is_empty() {
        return b_ret;
    }

    let mut an_array_of_lines: Vec<Vec<i32>> = vec![Vec::new(); (a_nb_pnts as usize) + 1];
    let mut an_array_of_line_type: Vec<i32> = vec![0; (a_nb_pnts as usize) + 1];
    let mut a_list_of_point_index: Vec<i32> = Vec::new();

    let mut nblines = 0i32;
    let a_tol = 0.5 * CONFUSION;
    let mut b_is_prev_point_on_boundary = false;
    let mut b_is_point_on_boundary = false;

    for pit in 1..=a_nb_pnts {
        let mut b_is_current_point_on_boundary = false;
        let a_point = the_wline.point(pit);
        for i in 0..2 {
            let a_ga = if i == 0 { the_surface1 } else { the_surface2 };
            let (umin, umax) = a_ga.u_range();
            let (vmin, vmax) = a_ga.v_range();
            let (u, v) = if i == 0 {
                a_point.parameters_on_s1()
            } else {
                a_point.parameters_on_s2()
            };
            for j in 0..2 {
                let isperiodic = if j == 0 {
                    a_ga.is_u_periodic()
                } else {
                    a_ga.is_v_periodic()
                };
                if !isperiodic {
                    continue;
                }
                let (a_resolution, a_period, alowerboundary, aupperboundary, a_parameter) =
                    if j == 0 {
                        (
                            u_resolution(a_ga, a_tol),
                            u_period(a_ga).unwrap_or(0.0),
                            umin,
                            umax,
                            u,
                        )
                    } else {
                        (
                            v_resolution(a_ga, a_tol),
                            v_period(a_ga).unwrap_or(0.0),
                            vmin,
                            vmax,
                            v,
                        )
                    };
                let (an_adjust_par, _anoffset) =
                    adjust_periodic(a_parameter, alowerboundary, aupperboundary, a_period, 0.0);
                let mut b_is_on_first_boundary = true;
                b_is_point_on_boundary = is_point_on_boundary(
                    an_adjust_par,
                    alowerboundary,
                    aupperboundary,
                    a_resolution,
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
            if !a_list_of_point_index.is_empty() {
                nblines += 1;
                an_array_of_lines[nblines as usize] = a_list_of_point_index.clone();
                an_array_of_line_type[nblines as usize] = i32::from(b_is_prev_point_on_boundary);
                a_list_of_point_index.clear();
            }
            b_is_prev_point_on_boundary = b_is_current_point_on_boundary;
        }
        a_list_of_point_index.push(pit);
    }
    if !a_list_of_point_index.is_empty() {
        nblines += 1;
        an_array_of_lines[nblines as usize] = a_list_of_point_index.clone();
        an_array_of_line_type[nblines as usize] = i32::from(b_is_prev_point_on_boundary);
        a_list_of_point_index.clear();
    }
    if nblines <= 1 {
        return b_ret;
    }

    // 2. Correct wlines.begin
    let mut an_array_of_line_ends: Vec<Vec<i32>> = vec![Vec::new(); (nblines as usize) + 1];
    let mut a_seq_of_pnt_on_2s: Vec<PntOn2S> = Vec::new();

    for i in 1..=nblines {
        if an_array_of_line_type[i as usize] != 0 {
            continue;
        }
        let a_list_of_index = &an_array_of_lines[i as usize];
        let mut a_list_of_fl_index: Vec<i32> = Vec::new();

        for j in 0..2 {
            let aneighbourindex = if j == 0 { i - 1 } else { i + 1 };
            if aneighbourindex < 1 || aneighbourindex > nblines {
                continue;
            }
            if an_array_of_line_type[aneighbourindex as usize] == 0 {
                continue;
            }
            let a_neighbour = &an_array_of_lines[aneighbourindex as usize];
            if a_neighbour.is_empty() {
                continue;
            }
            let an_index = if j == 0 {
                *a_neighbour.last().unwrap()
            } else {
                a_neighbour[0]
            };
            let a_point = *the_wline.point(an_index);
            let mut a_new_p = a_point;
            if a_list_of_index.len() < 2 {
                a_seq_of_pnt_on_2s.push(a_new_p);
                a_list_of_fl_index.push(a_seq_of_pnt_on_2s.len() as i32);
                continue;
            }
            let i_first = a_list_of_index[0];
            let i_last = *a_list_of_index.last().unwrap();

            for surfit in 0..2 {
                let a_ga = if surfit == 0 {
                    the_surface1
                } else {
                    the_surface2
                };
                let (umin, umax) = a_ga.u_range();
                let (vmin, vmax) = a_ga.v_range();
                let (u, v) = if surfit == 0 {
                    a_new_p.parameters_on_s1()
                } else {
                    a_new_p.parameters_on_s2()
                };
                let mut nbboundaries = 0i32;
                let mut b_is_near_boundary = false;
                let mut b_is_u_boundary = false;
                let mut b_is_first_boundary = false;

                for parit in 0..2 {
                    let isperiodic = if parit == 0 {
                        a_ga.is_u_periodic()
                    } else {
                        a_ga.is_v_periodic()
                    };
                    let a_resolution = if parit == 0 {
                        u_resolution(a_ga, a_tol)
                    } else {
                        v_resolution(a_ga, a_tol)
                    };
                    let alowerboundary = if parit == 0 { umin } else { vmin };
                    let aupperboundary = if parit == 0 { umax } else { vmax };
                    let a_parameter = if parit == 0 { u } else { v };
                    let mut b_is_on_first_boundary = true;

                    if !isperiodic {
                        b_is_point_on_boundary = is_point_on_boundary(
                            a_parameter,
                            alowerboundary,
                            aupperboundary,
                            a_resolution,
                            &mut b_is_on_first_boundary,
                        );
                        if b_is_point_on_boundary {
                            b_is_u_boundary = parit == 0;
                            b_is_first_boundary = b_is_on_first_boundary;
                            nbboundaries += 1;
                        }
                    } else {
                        let a_period = if parit == 0 {
                            u_period(a_ga).unwrap_or(0.0)
                        } else {
                            v_period(a_ga).unwrap_or(0.0)
                        };
                        let (an_adjust_par, _anoffset) =
                            adjust_periodic(a_parameter, alowerboundary, aupperboundary, a_period, 0.0);
                        b_is_point_on_boundary = is_point_on_boundary(
                            an_adjust_par,
                            alowerboundary,
                            aupperboundary,
                            a_resolution,
                            &mut b_is_on_first_boundary,
                        );
                        if b_is_point_on_boundary {
                            b_is_u_boundary = parit == 0;
                            b_is_first_boundary = b_is_on_first_boundary;
                            nbboundaries += 1;
                        } else {
                            let mut an_epsilon = a_resolution * 100.0;
                            let a_part = (aupperboundary - alowerboundary) * 0.1;
                            an_epsilon = if an_epsilon > a_part { a_part } else { an_epsilon };
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
                        a_ga.is_u_periodic()
                    } else {
                        a_ga.is_v_periodic()
                    };
                    if isperiodic {
                        let alowerboundary = if b_is_u_boundary { umin } else { vmin };
                        let aupperboundary = if b_is_u_boundary { umax } else { vmax };
                        let a_period = if b_is_u_boundary {
                            u_period(a_ga).unwrap_or(0.0)
                        } else {
                            v_period(a_ga).unwrap_or(0.0)
                        };
                        let a_parameter = if b_is_u_boundary { u } else { v };
                        let (an_adjust_par, anoffset) =
                            adjust_periodic(a_parameter, alowerboundary, aupperboundary, a_period, 0.0);
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
                        let a_neighbour_point = the_wline.point(aneighbourpointindex);
                        let (n_u1, n_v1) = if surfit == 0 {
                            a_neighbour_point.parameters_on_s1()
                        } else {
                            a_neighbour_point.parameters_on_s2()
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
                        let mut a_new_vec = GpVec2d::zero();
                        let anew_u = if b_is_u_boundary { another_par } else { u };
                        let anew_v = if b_is_u_boundary { v } else { another_par };

                        if (adist1 - adist2) > PCONFUSION && adist2 < (a_period / 4.0) {
                            b_check_angle1 = true;
                            a_new_vec = GpVec2d::new(anew_u - n_u1, anew_v - n_v1);
                            if a_new_vec.square_magnitude() < RESOLUTION {
                                a_new_p.set_uv_on(surfit == 0, anew_u, anew_v);
                                b_check_angle1 = false;
                            }
                        } else if adist1 < (a_period / 4.0) {
                            b_check_angle2 = true;
                            a_new_vec = GpVec2d::new(u - n_u1, v - n_v1);
                            if a_new_vec.square_magnitude() < RESOLUTION {
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
                                let a_prev = the_wline.point(anindexother);
                                let (n_u2, n_v2) = if surfit == 0 {
                                    a_prev.parameters_on_s1()
                                } else {
                                    a_prev.parameters_on_s2()
                                };
                                let a_vec_old = GpVec2d::new(n_u1 - n_u2, n_v1 - n_v2);
                                if a_vec_old.square_magnitude() <= RESOLUTION {
                                    continue;
                                }
                                let an_angle = a_new_vec.angle(&a_vec_old);
                                if an_angle.abs() < (std::f64::consts::PI * 0.25)
                                    && a_new_vec.dot(&a_vec_old) > 0.0
                                {
                                    if b_check_angle1 {
                                        let mut atmppoint = a_new_p;
                                        atmppoint.set_uv_on(surfit == 0, anew_u, anew_v);
                                        let (u1, v1, u2, v2) = atmppoint.parameters();
                                        let p1 = the_surface1.d0(u1, v1);
                                        let p2 = the_surface2.d0(u2, v2);
                                        let p0 = a_point.value();
                                        if p0.distance(&p1) <= a_tol
                                            && p0.distance(&p2) <= a_tol
                                            && p1.distance(&p2) <= a_tol
                                        {
                                            b_compute_line_end = false;
                                            a_new_p.set_uv_on(surfit == 0, anew_u, anew_v);
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
                        let (u1, v1, u2, v2) = a_new_p.parameters();
                        anewpoint = if surfit == 0 {
                            GpPnt2d::new(u1, v1)
                        } else {
                            GpPnt2d::new(u2, v2)
                        };
                        let aneighbourpointindex1 = if j == 0 { i_first } else { i_last };
                        let a_neighbour_point = the_wline.point(aneighbourpointindex1);
                        let (n_u1, n_v1) = if surfit == 0 {
                            a_neighbour_point.parameters_on_s1()
                        } else {
                            a_neighbour_point.parameters_on_s2()
                        };
                        let ap1 = GpPnt2d::new(n_u1, n_v1);
                        if a_ga.is_u_periodic() || a_ga.is_v_periodic() {
                            let ap2 = adjust_by_neighbour(&ap1, &anewpoint, a_ga);
                            if ap2.x() < umin || ap2.x() > umax || ap2.y() < vmin || ap2.y() > vmax
                            {
                                if let Some(p) = find_point(&ap1, &ap2, umin, umax, vmin, vmax) {
                                    anewpoint = p;
                                    found = true;
                                }
                            } else {
                                anewpoint = ap2;
                                a_new_p.set_uv_on(surfit == 0, anewpoint.x(), anewpoint.y());
                            }
                        }
                    } else {
                        let aneighbourpointindex1 = if j == 0 { i_first } else { i_last };
                        let a_neighbour_point = the_wline.point(aneighbourpointindex1);
                        let (n_u1, n_v1) = if surfit == 0 {
                            a_neighbour_point.parameters_on_s1()
                        } else {
                            a_neighbour_point.parameters_on_s2()
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
                            let a_prev = the_wline.point(aneighbourpointindex2);
                            let (n_u2, n_v2) = if surfit == 0 {
                                a_prev.parameters_on_s1()
                            } else {
                                a_prev.parameters_on_s2()
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
                        let a_criteria = the_tol;
                        let a_surface = if surfit == 0 {
                            the_surface1
                        } else {
                            the_surface2
                        };
                        let a_surface_other = if surfit == 0 {
                            the_surface2
                        } else {
                            the_surface1
                        };
                        let a_p3d = a_surface.d0(anewpoint.x(), anewpoint.y());
                        if let Some(proj) = project_point_on_surface(a_surface_other, &a_p3d, 0.0) {
                            if proj.distance < a_criteria {
                                let mut found_u = u;
                                let mut found_v = v;
                                found_u = proj.u;
                                found_v = proj.v;
                                let aneindex1 = if j == 0 { i_first } else { i_last };
                                let a_neighbour_point = the_wline.point(aneindex1);
                                let (n_un, n_vn) = if surfit == 0 {
                                    a_neighbour_point.parameters_on_s2()
                                } else {
                                    a_neighbour_point.parameters_on_s1()
                                };
                                let a_neighbour_2d = GpPnt2d::new(n_un, n_vn);
                                let an_adjusted_point = adjust_by_neighbour(
                                    &a_neighbour_2d,
                                    &GpPnt2d::new(found_u, found_v),
                                    a_surface_other,
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
                                    let p3d = a_surface_other.d0(found_u, found_v);
                                    if let Some(proj2) =
                                        project_point_on_surface(a_surface, &p3d, 0.0)
                                    {
                                        if proj2.distance < a_criteria {
                                            anewpoint = GpPnt2d::new(proj2.u, proj2.v);
                                        }
                                    }
                                }
                                if surfit == 0 {
                                    a_new_p.set_value(
                                        a_p3d,
                                        anewpoint.x(),
                                        anewpoint.y(),
                                        found_u,
                                        found_v,
                                    );
                                } else {
                                    a_new_p.set_value(
                                        a_p3d,
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
            a_seq_of_pnt_on_2s.push(a_new_p);
            a_list_of_fl_index.push(a_seq_of_pnt_on_2s.len() as i32);
        }
        an_array_of_line_ends[i as usize] = a_list_of_fl_index;
    }

    // Split wlines.begin
    let nbiter = if b_avoid_line_constructor {
        1
    } else {
        line_parts.len().max(1)
    };
    for j in 1..=nbiter {
        let (ifprm, ilprm) = if b_avoid_line_constructor {
            (1, the_wline.nb_pnts())
        } else {
            line_parts[(j - 1) as usize]
        };
        let mut a_line_on_2s = WLine::new();
        for i in 1..=nblines {
            if an_array_of_line_type[i as usize] != 0 {
                continue;
            }
            let a_list_of_index = &an_array_of_lines[i as usize];
            let a_list_of_fl_index = &an_array_of_line_ends[i as usize];
            let mut bhasfirstpoint = a_list_of_fl_index.len() == 2;
            let mut bhaslastpoint = a_list_of_fl_index.len() == 2;
            if !bhasfirstpoint && !a_list_of_fl_index.is_empty() {
                bhasfirstpoint = i != 1;
            }
            if !bhaslastpoint && !a_list_of_fl_index.is_empty() {
                bhaslastpoint = i != nblines;
            }
            if a_list_of_index.is_empty() {
                continue;
            }
            let i_first = a_list_of_index[0];
            let i_last = *a_list_of_index.last().unwrap();
            let b_is_first_inside = ifprm >= i_first && ifprm <= i_last;
            let b_is_last_inside = ilprm >= i_first && ilprm <= i_last;

            if !b_is_first_inside && !b_is_last_inside {
                if ifprm < i_first && ilprm > i_last {
                    if bhasfirstpoint {
                        let pit = a_list_of_fl_index[0];
                        a_line_on_2s.add(a_seq_of_pnt_on_2s[(pit as usize).saturating_sub(1)]);
                    }
                    for &pit in a_list_of_index {
                        a_line_on_2s.add(*the_wline.point(pit));
                    }
                    if bhaslastpoint {
                        let pit = *a_list_of_fl_index.last().unwrap();
                        a_line_on_2s.add(a_seq_of_pnt_on_2s[(pit as usize).saturating_sub(1)]);
                    }
                    let aneighbour = i + 1;
                    let mut b_is_end_of_line = true;
                    if aneighbour <= nblines {
                        let a_list_of_neighbour = &an_array_of_lines[aneighbour as usize];
                        if an_array_of_line_type[aneighbour as usize] != 0
                            && a_list_of_neighbour.is_empty()
                        {
                            b_is_end_of_line = false;
                        }
                    }
                    if b_is_end_of_line {
                        if a_line_on_2s.nb_pnts() > 1 {
                            the_new_lines.push(std::mem::take(&mut a_line_on_2s));
                        }
                        a_line_on_2s = WLine::new();
                    }
                }
                continue;
            }
            if b_is_first_inside && b_is_last_inside {
                for &pit in a_list_of_index {
                    if pit < ifprm || pit > ilprm {
                        continue;
                    }
                    a_line_on_2s.add(*the_wline.point(pit));
                }
            } else {
                if b_is_first_inside {
                    for &pit in a_list_of_index {
                        if pit < ifprm {
                            continue;
                        }
                        a_line_on_2s.add(*the_wline.point(pit));
                    }
                    if bhaslastpoint {
                        let pit = *a_list_of_fl_index.last().unwrap();
                        a_line_on_2s.add(a_seq_of_pnt_on_2s[(pit as usize).saturating_sub(1)]);
                    }
                    let aneighbour = i + 1;
                    let mut b_is_end_of_line = true;
                    if aneighbour <= nblines {
                        let a_list_of_neighbour = &an_array_of_lines[aneighbour as usize];
                        if an_array_of_line_type[aneighbour as usize] != 0
                            && a_list_of_neighbour.is_empty()
                        {
                            b_is_end_of_line = false;
                        }
                    }
                    if b_is_end_of_line {
                        if a_line_on_2s.nb_pnts() > 1 {
                            the_new_lines.push(std::mem::take(&mut a_line_on_2s));
                        }
                        a_line_on_2s = WLine::new();
                    }
                }
                if b_is_last_inside {
                    if bhasfirstpoint {
                        let pit = a_list_of_fl_index[0];
                        a_line_on_2s.add(a_seq_of_pnt_on_2s[(pit as usize).saturating_sub(1)]);
                    }
                    for &pit in a_list_of_index {
                        if pit > ilprm {
                            continue;
                        }
                        a_line_on_2s.add(*the_wline.point(pit));
                    }
                }
            }
        }
        if a_line_on_2s.nb_pnts() > 1 {
            the_new_lines.push(a_line_on_2s);
        }
    }
    true
}
