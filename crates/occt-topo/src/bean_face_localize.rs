//! Localized bean/face intersection: `LocalizeSolutions` + `ComputeLocalized`.
//!
//! Source: `IntTools_BeanFaceIntersector.cxx:1369-2125` and `MergeSolutions`
//! at 2549. The previous 3-cell box cull is replaced by OCCT's recursive
//! range-sample localization.

use occt_core::bnd::BndBox;
use occt_core::precision::PCONFUSION;
use occt_geom::extrema_surf::curve_surface_extrema_all;

use crate::bean_face::BeanFaceIntersector;
use crate::bean_face_grid::{add_surface_to_box, get_surface_box};
use crate::bean_face_sample::{
    check_sampling, CurveRangeLocalizeData, CurveRangeSample, SurfaceRangeLocalizeData,
    SurfaceRangeSample,
};
use crate::int_tools_curve_box::add_curve_to_box;
use crate::inttools_data::IntRange;
use crate::meshing::range_splitter::{classify_surface as classify_surface_mesh, SurfaceType};

impl BeanFaceIntersector {
    /// Localized intersection for high-degree NURBS-like surfaces.
    /// Port of `ComputeLocalized`.
    pub(crate) fn compute_localized(&mut self) -> bool {
        let d_min_u = 10.0 * PCONFUSION;
        let d_min_v = d_min_u;
        let mut a_surface_data = SurfaceRangeLocalizeData::with_samples(3, 3, d_min_u, d_min_v);
        a_surface_data.remove_range_out_all();
        a_surface_data.clear_grid();

        let a_surface_range = SurfaceRangeSample::with_indices(0, 0, 0, 0);
        let mut fbox = a_surface_data
            .find_box(&a_surface_range)
            .unwrap_or_else(BndBox::new);
        let b_fbox_found = a_surface_data.find_box(&a_surface_range).is_some();

        let is_bspline = matches!(
            classify_surface_mesh(self.surface()),
            SurfaceType::BSplineSurface
        );
        if is_bspline && a_surface_data.has_grid() {
            if !b_fbox_found {
                fbox = get_surface_box(
                    self.surface(),
                    self.umin,
                    self.umax,
                    self.vmin,
                    self.vmax,
                    self.criteria,
                    &mut a_surface_data,
                );
                a_surface_data.add_box(a_surface_range, fbox);
            }
        } else if !b_fbox_found {
            add_surface_to_box(
                self.surface(),
                self.umin,
                self.umax,
                self.vmin,
                self.vmax,
                self.face_tolerance,
                &mut fbox,
            );
            a_surface_data.add_box(a_surface_range, fbox);
        }

        let mut ebox = BndBox::new();
        add_curve_to_box(
            self.curve(),
            self.first_parameter,
            self.last_parameter,
            self.bean_tolerance,
            &mut ebox,
        );

        if ebox.is_out_box(&fbox) {
            for i in 0..self.range_manager.len() {
                self.range_manager.set_flag(i, 1);
            }
            a_surface_data.clear_grid();
            return true;
        }

        let a_curve_range = {
            let mut c = CurveRangeSample::with_index(0);
            c.set_depth(0);
            c
        };
        let nb_sample_c = 3i32;
        let nb_sample_u = a_surface_data.nb_sample_u();
        let nb_sample_v = a_surface_data.nb_sample_v();
        let d_min_c = 10.0 * self.curve_resolution;

        let a_curve_data_tmp = CurveRangeLocalizeData::new(nb_sample_c, d_min_c);
        let a_surface_data_tmp =
            SurfaceRangeLocalizeData::with_samples(nb_sample_u, nb_sample_v, d_min_u, d_min_v);
        let (b_allow_c, b_allow_u, b_allow_v) = check_sampling(
            &a_curve_range,
            &a_surface_range,
            &a_curve_data_tmp,
            &a_surface_data_tmp,
            self.last_parameter - self.first_parameter,
            self.umax - self.umin,
            self.vmax - self.vmin,
        );

        let mut a_list_curve: Vec<CurveRangeSample> = Vec::new();
        let mut a_list_surface: Vec<SurfaceRangeSample> = Vec::new();
        let mut a_curve_data = CurveRangeLocalizeData::new(nb_sample_c, d_min_c);
        a_curve_data.add_box(a_curve_range, ebox);

        if !self.localize_solutions(
            a_curve_range,
            ebox,
            a_surface_range,
            fbox,
            &mut a_curve_data,
            &mut a_surface_data,
            &mut a_list_curve,
            &mut a_list_surface,
        ) {
            a_surface_data.clear_grid();
            return false;
        }

        let (a_list_curve_sort, a_list_surface_sort) =
            merge_solutions(&a_list_curve, &a_list_surface);

        let mut a_range_s_prev: Option<SurfaceRangeSample> = None;
        for (it_c, it_s) in a_list_curve_sort.iter().zip(a_list_surface_sort.iter()) {
            let mut a_range_c =
                IntRange::new_unchecked(self.first_parameter, self.last_parameter);
            if b_allow_c {
                a_range_c = it_c.get_range(self.first_parameter, self.last_parameter, nb_sample_c);
            }
            let mut a_range_u = IntRange::new_unchecked(self.umin, self.umax);
            if b_allow_u {
                a_range_u = it_s.get_range_u(self.umin, self.umax, nb_sample_u);
            }
            let mut a_range_v = IntRange::new_unchecked(self.vmin, self.vmax);
            if b_allow_v {
                a_range_v = it_s.get_range_v(self.vmin, self.vmax, nb_sample_v);
            }
            let anarg1 = a_range_c.first;
            let anarg2 = a_range_c.last;

            let mut n_min_index = self.range_manager.len() as i32;
            let mut n_max_index = -1i32;
            let inds1 = self.range_manager.get_indices(anarg1);
            for &n_index in &inds1 {
                let n = n_index as i32;
                n_min_index = n_min_index.min(n);
                n_max_index = n_max_index.max(n);
            }
            let mut b_found = false;
            if n_max_index >= 0 {
                for ind in n_min_index..=n_max_index {
                    if self.range_manager.flag(ind as usize) == 2 {
                        b_found = true;
                        break;
                    }
                }
            }
            if b_found {
                continue;
            }
            n_min_index = if n_max_index >= 0 {
                n_max_index
            } else {
                n_min_index
            };
            let inds2 = self.range_manager.get_indices(anarg2);
            for &n_index in &inds2 {
                let n = n_index as i32;
                n_min_index = n_min_index.min(n);
                n_max_index = n_max_index.max(n);
            }
            if n_max_index >= 0 {
                for ind in n_min_index..=n_max_index {
                    if self.range_manager.flag(ind as usize) == 2 {
                        b_found = true;
                        break;
                    }
                }
            }
            if b_found {
                continue;
            }

            let par_uf = a_range_u.first;
            let par_ul = a_range_u.last;
            let par_vf = a_range_v.first;
            let par_vl = a_range_v.last;
            let reuse = a_range_s_prev.map(|p| p.is_equal(it_s)).unwrap_or(false);
            let extrema_ok = self.gen_ext_cs_cell(
                anarg1,
                anarg2,
                par_uf,
                par_ul,
                par_vf,
                par_vl,
                reuse,
            );
            if !extrema_ok {
                let _ = self.range_manager.insert_range(anarg1, anarg2, 0);
            }
            a_range_s_prev = Some(*it_s);
        }

        if b_allow_c {
            for out_c in a_curve_data.list_range_out() {
                let a_range_c =
                    out_c.get_range(self.first_parameter, self.last_parameter, nb_sample_c);
                let _ = self
                    .range_manager
                    .insert_range(a_range_c.first, a_range_c.last, 1);
            }
        }
        self.compute_near_range_boundaries();
        a_surface_data.clear_grid();
        true
    }

    /// `LocalizeSolutions`.
    fn localize_solutions(
        &self,
        the_curve_range: CurveRangeSample,
        the_box_curve: BndBox,
        the_surface_range: SurfaceRangeSample,
        the_box_surface: BndBox,
        the_curve_data: &mut CurveRangeLocalizeData,
        the_surface_data: &mut SurfaceRangeLocalizeData,
        the_list_curve_range: &mut Vec<CurveRangeSample>,
        the_list_surface_range: &mut Vec<SurfaceRangeSample>,
    ) -> bool {
        let a_root_range_c = {
            let mut c = CurveRangeSample::with_index(0);
            c.set_depth(0);
            c
        };
        let a_root_range_s = SurfaceRangeSample::with_indices(0, 0, 0, 0);
        let mut a_main_box_c = the_box_curve;
        let mut a_main_box_s = the_box_surface;
        let mut b_main_box_found_s = false;
        let mut b_main_box_found_c = false;

        let mut a_list_curve_found: Vec<CurveRangeSample> = Vec::new();
        let mut a_list_surface_found: Vec<SurfaceRangeSample> = Vec::new();

        let a_range_c = the_curve_range.get_range(
            self.first_parameter,
            self.last_parameter,
            the_curve_data.nb_sample(),
        );
        let mut local_diff_c =
            (a_range_c.last - a_range_c.first) / the_curve_data.nb_sample() as f64;
        let mut a_cur_index_init =
            the_curve_range.range_index_deeper(the_curve_data.nb_sample());
        let mut a_list_c_to_avoid: Vec<i32> = Vec::new();
        let mut b_global_check_done = false;

        let mut a_cur_index_u =
            the_surface_range.range_index_u_deeper(the_surface_data.nb_sample_u());
        let mut a_cur_index_v_init =
            the_surface_range.range_index_v_deeper(the_surface_data.nb_sample_v());
        let a_range_v = the_surface_range.get_range_v(
            self.vmin,
            self.vmax,
            the_surface_data.nb_sample_v(),
        );
        let a_range_u = the_surface_range.get_range_u(
            self.umin,
            self.umax,
            the_surface_data.nb_sample_u(),
        );
        let mut a_cur_par_u = a_range_u.first;
        let mut a_local_diff_u =
            (a_range_u.last - a_range_u.first) / the_surface_data.nb_sample_u() as f64;
        let mut a_prev_par_u = a_cur_par_u;
        let mut a_local_diff_v =
            (a_range_v.last - a_range_v.first) / the_surface_data.nb_sample_v() as f64;

        let (b_allow_c, b_allow_u, b_allow_v) = check_sampling(
            &the_curve_range,
            &the_surface_range,
            the_curve_data,
            the_surface_data,
            local_diff_c,
            a_local_diff_u,
            a_local_diff_v,
        );
        if !b_allow_c && !b_allow_u && !b_allow_v {
            the_list_curve_range.push(the_curve_range);
            the_list_surface_range.push(the_surface_range);
            return true;
        }

        let mut a_new_range_c_template = CurveRangeSample::new();
        if !b_allow_c {
            a_new_range_c_template = the_curve_range;
            a_cur_index_init = the_curve_range.range_index();
            local_diff_c = a_range_c.last - a_range_c.first;
        } else {
            a_new_range_c_template.set_depth(the_curve_range.depth() + 1);
            a_new_range_c_template.set_range_index(a_cur_index_init);
        }

        let mut a_new_range_s_template = the_surface_range;
        if b_allow_u {
            a_new_range_s_template.set_depth_u(the_surface_range.depth_u() + 1);
        } else {
            a_cur_index_u = a_new_range_s_template.index_u();
            a_local_diff_u = a_range_u.last - a_range_u.first;
        }
        if b_allow_v {
            a_new_range_s_template.set_depth_v(the_surface_range.depth_v() + 1);
        } else {
            a_cur_index_v_init = the_surface_range.index_v();
            a_local_diff_v = a_range_v.last - a_range_v.first;
        }

        let mut b_has_out = false;
        let nb_u = if b_allow_u {
            the_surface_data.nb_sample_u()
        } else {
            1
        };
        let nb_v = if b_allow_v {
            the_surface_data.nb_sample_v()
        } else {
            1
        };
        let nb_c = if b_allow_c {
            the_curve_data.nb_sample()
        } else {
            1
        };

        for _u_it in 1..=nb_u {
            a_prev_par_u = a_cur_par_u;
            a_cur_par_u += a_local_diff_u;
            let mut a_cur_par_v = a_range_v.first;
            let mut a_prev_par_v = a_cur_par_v;
            let mut a_cur_index_v = a_cur_index_v_init;
            let mut b_has_out_v = false;

            for _v_it in 1..=nb_v {
                a_prev_par_v = a_cur_par_v;
                a_cur_par_v += a_local_diff_v;
                let mut a_new_range_s = a_new_range_s_template;
                if b_allow_u {
                    a_new_range_s.set_index_u(a_cur_index_u);
                }
                if b_allow_v {
                    a_new_range_s.set_index_v(a_cur_index_v);
                }
                a_cur_index_v += 1;

                if the_surface_data.is_range_out(&a_new_range_s) {
                    b_has_out_v = true;
                    continue;
                }

                let a_box_s = match the_surface_data.find_box(&a_new_range_s) {
                    Some(b) => b,
                    None => {
                        let mut box_s = BndBox::new();
                        let is_bspline = matches!(
                            classify_surface_mesh(self.surface()),
                            SurfaceType::BSplineSurface
                        );
                        if is_bspline && the_surface_data.has_grid() {
                            box_s = get_surface_box(
                                self.surface(),
                                a_prev_par_u,
                                a_cur_par_u,
                                a_prev_par_v,
                                a_cur_par_v,
                                self.criteria,
                                the_surface_data,
                            );
                        } else {
                            add_surface_to_box(
                                self.surface(),
                                a_prev_par_u,
                                a_cur_par_u,
                                a_prev_par_v,
                                a_cur_par_v,
                                self.criteria,
                                &mut box_s,
                            );
                        }
                        if !b_main_box_found_c {
                            if let Some(b) = the_curve_data.find_box(&a_root_range_c) {
                                a_main_box_c = b;
                                b_main_box_found_c = true;
                            }
                        }
                        if box_s.is_out_box(&a_main_box_c) {
                            the_surface_data.add_out_range(a_new_range_s);
                            b_has_out_v = true;
                            continue;
                        }
                        the_surface_data.add_box(a_new_range_s, box_s);
                        box_s
                    }
                };

                if a_box_s.is_out_box(&the_box_curve) {
                    b_has_out_v = true;
                    continue;
                }

                let mut a_list_of_box: Vec<BndBox> = Vec::new();
                let mut a_list_of_index: Vec<i32> = Vec::new();
                let mut b_has_out_c = false;
                let mut a_cur_par = a_range_c.first;
                let mut a_prev_par = a_range_c.first;
                let mut a_cur_range_c = a_new_range_c_template;
                let mut a_cur_index = a_cur_index_init;

                for t_it in 1..=nb_c {
                    a_prev_par = a_cur_par;
                    a_cur_par += local_diff_c;
                    let mut b_found = a_list_c_to_avoid.contains(&t_it);
                    if !b_found {
                        if b_allow_c {
                            a_cur_range_c.set_range_index(a_cur_index);
                        }
                        b_found = the_curve_data.is_range_out(&a_cur_range_c);
                    }
                    a_cur_index += 1;
                    if b_found {
                        b_has_out_c = true;
                        continue;
                    }

                    let a_box_c = match the_curve_data.find_box(&a_cur_range_c) {
                        Some(b) => b,
                        None => {
                            let mut box_c = BndBox::new();
                            add_curve_to_box(
                                self.curve(),
                                a_prev_par,
                                a_cur_par,
                                self.criteria,
                                &mut box_c,
                            );
                            if !b_main_box_found_s {
                                if let Some(b) = the_surface_data.find_box(&a_root_range_s) {
                                    a_main_box_s = b;
                                    b_main_box_found_s = true;
                                }
                            }
                            if box_c.is_out_box(&a_main_box_s) {
                                the_curve_data.add_out_range(a_cur_range_c);
                                b_has_out_c = true;
                                continue;
                            }
                            the_curve_data.add_box(a_cur_range_c, box_c);
                            box_c
                        }
                    };

                    if !b_global_check_done && a_box_c.is_out_box(&the_box_surface) {
                        a_list_c_to_avoid.push(t_it);
                        b_has_out_c = true;
                        continue;
                    }
                    if a_box_c.is_out_box(&a_box_s) {
                        b_has_out_v = true;
                        b_has_out_c = true;
                        continue;
                    }
                    a_list_of_index.push(t_it);
                    a_list_of_box.push(a_box_c);
                }
                b_global_check_done = true;
                if b_has_out_c {
                    b_has_out_v = true;
                }

                let mut a_new_range_c = a_new_range_c_template;
                let mut b_use_old_c = false;
                let mut b_use_old_s = false;
                let b_check_size = !b_has_out_c;

                for (t_it, a_box_c) in a_list_of_index.iter().zip(a_list_of_box.iter()) {
                    a_cur_index = a_cur_index_init + t_it - 1;
                    b_use_old_s = false;
                    if b_allow_c {
                        a_new_range_c.set_range_index(a_cur_index);
                    }
                    if b_check_size {
                        if the_curve_range.depth() == 0
                            || the_surface_range.depth_u() == 0
                            || the_surface_range.depth_v() == 0
                        {
                            b_has_out_c = true;
                            b_has_out_v = true;
                        } else if the_curve_range.depth() < 4
                            && the_surface_range.depth_u() < 4
                            && the_surface_range.depth_v() < 4
                        {
                            if !a_box_c.is_whole() && !a_box_s.is_whole() {
                                let a_diag_c = a_box_c.square_extent();
                                let a_diag_s = a_box_s.square_extent();
                                if a_diag_c < a_diag_s {
                                    if a_diag_c * 10.0 < a_diag_s {
                                        b_use_old_c = true;
                                        b_has_out_c = true;
                                        b_has_out_v = true;
                                        break;
                                    }
                                } else if a_diag_s * 10.0 < a_diag_c {
                                    b_use_old_s = true;
                                    b_has_out_c = true;
                                    b_has_out_v = true;
                                }
                            }
                        }
                    }
                    if !b_has_out_c {
                        a_list_curve_found.push(a_new_range_c);
                        a_list_surface_found.push(a_new_range_s);
                    } else {
                        if b_use_old_s && a_new_range_c == the_curve_range {
                            return false;
                        }
                        let surf_range = if b_use_old_s {
                            the_surface_range
                        } else {
                            a_new_range_s
                        };
                        let surf_box = if b_use_old_s {
                            the_box_surface
                        } else {
                            a_box_s
                        };
                        if !self.localize_solutions(
                            a_new_range_c,
                            *a_box_c,
                            surf_range,
                            surf_box,
                            the_curve_data,
                            the_surface_data,
                            the_list_curve_range,
                            the_list_surface_range,
                        ) {
                            return false;
                        }
                    }
                }

                if b_has_out_v
                    && b_use_old_c
                    && b_allow_c
                    && (b_allow_u || b_allow_v)
                    && !self.localize_solutions(
                        the_curve_range,
                        the_box_curve,
                        a_new_range_s,
                        a_box_s,
                        the_curve_data,
                        the_surface_data,
                        the_list_curve_range,
                        the_list_surface_range,
                    )
                {
                    return false;
                }
                let _ = a_box_s;
            }
            a_cur_index_u += 1;
            if b_has_out_v {
                b_has_out = true;
            }
        }

        if !b_has_out {
            the_list_curve_range.push(the_curve_range);
            the_list_surface_range.push(the_surface_range);
        } else {
            for (c, s) in a_list_curve_found.iter().zip(a_list_surface_found.iter()) {
                the_list_curve_range.push(*c);
                the_list_surface_range.push(*s);
            }
        }
        true
    }

    /// `Extrema_GenExtCS` analogue for one localized cell: 10 curve samples
    /// plus `curve_surface_extrema_all`, restricted to the cell.
    fn gen_ext_cs_cell(
        &mut self,
        anarg1: f64,
        anarg2: f64,
        par_uf: f64,
        par_ul: f64,
        par_vf: f64,
        par_vl: f64,
        _reuse: bool,
    ) -> bool {
        let mut any = false;
        let mut candidates: Vec<(f64, f64, f64, f64)> = Vec::new();
        let exts = curve_surface_extrema_all(self.curve(), self.surface());
        for e in &exts {
            if e.u1 < anarg1 - PCONFUSION || e.u1 > anarg2 + PCONFUSION {
                continue;
            }
            let u = e.u2;
            let v = e.v2.unwrap_or(0.0);
            if u < par_uf || u > par_ul || v < par_vf || v > par_vl {
                continue;
            }
            candidates.push((e.u1, u, v, e.distance * e.distance));
            any = true;
        }
        const NB: usize = 10;
        let span = anarg2 - anarg1;
        if span.abs() > PCONFUSION {
            let f = |t: f64| self.closest_params_dist(&self.curve_d0(t)).2;
            let (t_min, d_min) = crate::bean_face::golden_1d(&f, anarg1, anarg2, 1e-12);
            let (u, v, _) = self.closest_params_dist(&self.curve_d0(t_min));
            if u >= par_uf && u <= par_ul && v >= par_vf && v <= par_vl {
                candidates.push((t_min, u, v, d_min * d_min));
                any = true;
            }
            for i in 0..=NB {
                let t = anarg1 + span * i as f64 / NB as f64;
                let (u, v, d) = self.closest_params_dist(&self.curve_d0(t));
                if u < par_uf || u > par_ul || v < par_vf || v > par_vl {
                    continue;
                }
                candidates.push((t, u, v, d * d));
                any = true;
            }
        }
        if !any {
            return false;
        }
        let crit2 = self.criteria * self.criteria;
        let mut found_ext = false;
        for (t, mut u, mut v, sq) in candidates {
            if sq >= crit2 {
                continue;
            }
            found_ext = true;
            u = in_period_if(self.surface().is_u_periodic(), u, par_uf, par_ul);
            v = in_period_if(self.surface().is_v_periodic(), v, par_vf, par_vl);
            let mut t_adj = t;
            if self.curve().is_periodic() {
                t_adj = in_period(t, anarg1, anarg1 + self.curve().period());
            }
            if u < self.umin {
                u = self.umin;
            }
            if u > self.umax {
                u = self.umax;
            }
            if v < self.vmin {
                v = self.vmin;
            }
            if v > self.vmax {
                v = self.vmax;
            }
            let n = self.range_manager.len();
            self.compute_range_from_start_point(false, t_adj, u, v);
            self.compute_range_from_start_point(true, t_adj, u, v);
            if n == self.range_manager.len() {
                self.set_empty_result_range(t_adj);
            }
        }
        found_ext
    }
}

fn merge_solutions(
    list_curve: &[CurveRangeSample],
    list_surface: &[SurfaceRangeSample],
) -> (Vec<CurveRangeSample>, Vec<SurfaceRangeSample>) {
    let mut map_to_avoid: Vec<SurfaceRangeSample> = Vec::new();
    let mut curve_id_map: Vec<Vec<usize>> = Vec::new();
    let mut curve_range_vector: Vec<CurveRangeSample> = Vec::new();
    for (c, s) in list_curve.iter().zip(list_surface.iter()) {
        let id = curve_range_vector.len();
        curve_range_vector.push(*c);
        if let Some(idx) = map_to_avoid.iter().position(|x| x.is_equal(s)) {
            curve_id_map[idx].push(id);
        } else {
            map_to_avoid.push(*s);
            curve_id_map.push(vec![id]);
        }
    }
    let mut out_c = Vec::new();
    let mut out_s = Vec::new();
    for (i, surf) in map_to_avoid.iter().enumerate() {
        for &id in &curve_id_map[i] {
            out_s.push(*surf);
            out_c.push(curve_range_vector[id]);
        }
    }
    (out_c, out_s)
}

fn in_period(value: f64, first: f64, last: f64) -> f64 {
    let period = last - first;
    if period.abs() <= PCONFUSION {
        return value;
    }
    let mut v = value;
    while v < first {
        v += period;
    }
    while v > last {
        v -= period;
    }
    v
}

fn in_period_if(periodic: bool, value: f64, first: f64, last: f64) -> f64 {
    if periodic {
        in_period(value, first, last)
    } else {
        value
    }
}
