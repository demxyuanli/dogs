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
use crate::meshing::range_splitter::{mesh_surface_type as mesh_surface_type_mesh, MeshSurfaceType};

impl BeanFaceIntersector {
    /// Localized intersection for high-degree NURBS-like surfaces.
    /// Port of `ComputeLocalized`.
    pub(crate) fn compute_localized(&mut self) -> bool {
        let d_min_u = 10.0 * PCONFUSION;
        let d_min_v = d_min_u;
        let mut surface_data = SurfaceRangeLocalizeData::with_samples(3, 3, d_min_u, d_min_v);
        surface_data.remove_range_out_all();
        surface_data.clear_grid();

        let surface_range = SurfaceRangeSample::with_indices(0, 0, 0, 0);
        let mut fbox = surface_data
            .find_box(&surface_range)
            .unwrap_or_else(BndBox::new);
        let b_fbox_found = surface_data.find_box(&surface_range).is_some();

        let is_bspline = matches!(
            mesh_surface_type_mesh(self.surface()),
            MeshSurfaceType::BSplineSurface
        );
        if is_bspline && surface_data.has_grid() {
            if !b_fbox_found {
                fbox = get_surface_box(
                    self.surface(),
                    self.umin,
                    self.umax,
                    self.vmin,
                    self.vmax,
                    self.criteria,
                    &mut surface_data,
                );
                surface_data.add_box(surface_range, fbox);
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
            surface_data.add_box(surface_range, fbox);
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
            surface_data.clear_grid();
            return true;
        }

        let curve_range = {
            let mut c = CurveRangeSample::with_index(0);
            c.set_depth(0);
            c
        };
        let nb_sample_c = 3i32;
        let nb_sample_u = surface_data.nb_sample_u();
        let nb_sample_v = surface_data.nb_sample_v();
        let d_min_c = 10.0 * self.curve_resolution;

        let curve_data_tmp = CurveRangeLocalizeData::new(nb_sample_c, d_min_c);
        let surface_data_tmp =
            SurfaceRangeLocalizeData::with_samples(nb_sample_u, nb_sample_v, d_min_u, d_min_v);
        let (b_allow_c, b_allow_u, b_allow_v) = check_sampling(
            &curve_range,
            &surface_range,
            &curve_data_tmp,
            &surface_data_tmp,
            self.last_parameter - self.first_parameter,
            self.umax - self.umin,
            self.vmax - self.vmin,
        );

        let mut list_curve: Vec<CurveRangeSample> = Vec::new();
        let mut list_surface: Vec<SurfaceRangeSample> = Vec::new();
        let mut curve_data = CurveRangeLocalizeData::new(nb_sample_c, d_min_c);
        curve_data.add_box(curve_range, ebox);

        if !self.localize_solutions(
            curve_range,
            ebox,
            surface_range,
            fbox,
            &mut curve_data,
            &mut surface_data,
            &mut list_curve,
            &mut list_surface,
        ) {
            surface_data.clear_grid();
            return false;
        }

        let (list_curve_sort, list_surface_sort) =
            merge_solutions(&list_curve, &list_surface);

        let mut range_s_prev: Option<SurfaceRangeSample> = None;
        for (it_c, it_s) in list_curve_sort.iter().zip(list_surface_sort.iter()) {
            let mut range_c =
                IntRange::new_unchecked(self.first_parameter, self.last_parameter);
            if b_allow_c {
                range_c = it_c.get_range(self.first_parameter, self.last_parameter, nb_sample_c);
            }
            let mut range_u = IntRange::new_unchecked(self.umin, self.umax);
            if b_allow_u {
                range_u = it_s.get_range_u(self.umin, self.umax, nb_sample_u);
            }
            let mut range_v = IntRange::new_unchecked(self.vmin, self.vmax);
            if b_allow_v {
                range_v = it_s.get_range_v(self.vmin, self.vmax, nb_sample_v);
            }
            let anarg1 = range_c.first;
            let anarg2 = range_c.last;

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

            let par_uf = range_u.first;
            let par_ul = range_u.last;
            let par_vf = range_v.first;
            let par_vl = range_v.last;
            let reuse = range_s_prev.map(|p| p.is_equal(it_s)).unwrap_or(false);
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
            range_s_prev = Some(*it_s);
        }

        if b_allow_c {
            for out_c in curve_data.list_range_out() {
                let range_c =
                    out_c.get_range(self.first_parameter, self.last_parameter, nb_sample_c);
                let _ = self
                    .range_manager
                    .insert_range(range_c.first, range_c.last, 1);
            }
        }
        self.compute_near_range_boundaries();
        surface_data.clear_grid();
        true
    }

    /// `LocalizeSolutions`.
    fn localize_solutions(
        &self,
        curve_range: CurveRangeSample,
        box_curve: BndBox,
        surface_range: SurfaceRangeSample,
        box_surface: BndBox,
        curve_data: &mut CurveRangeLocalizeData,
        surface_data: &mut SurfaceRangeLocalizeData,
        list_curve_range: &mut Vec<CurveRangeSample>,
        list_surface_range: &mut Vec<SurfaceRangeSample>,
    ) -> bool {
        let root_range_c = {
            let mut c = CurveRangeSample::with_index(0);
            c.set_depth(0);
            c
        };
        let root_range_s = SurfaceRangeSample::with_indices(0, 0, 0, 0);
        let mut main_box_c = box_curve;
        let mut main_box_s = box_surface;
        let mut b_main_box_found_s = false;
        let mut b_main_box_found_c = false;

        let mut list_curve_found: Vec<CurveRangeSample> = Vec::new();
        let mut list_surface_found: Vec<SurfaceRangeSample> = Vec::new();

        let range_c = curve_range.get_range(
            self.first_parameter,
            self.last_parameter,
            curve_data.nb_sample(),
        );
        let mut local_diff_c =
            (range_c.last - range_c.first) / curve_data.nb_sample() as f64;
        let mut cur_index_init =
            curve_range.range_index_deeper(curve_data.nb_sample());
        let mut list_c_to_avoid: Vec<i32> = Vec::new();
        let mut b_global_check_done = false;

        let mut cur_index_u =
            surface_range.range_index_u_deeper(surface_data.nb_sample_u());
        let mut cur_index_v_init =
            surface_range.range_index_v_deeper(surface_data.nb_sample_v());
        let range_v = surface_range.get_range_v(
            self.vmin,
            self.vmax,
            surface_data.nb_sample_v(),
        );
        let range_u = surface_range.get_range_u(
            self.umin,
            self.umax,
            surface_data.nb_sample_u(),
        );
        let mut cur_par_u = range_u.first;
        let mut local_diff_u =
            (range_u.last - range_u.first) / surface_data.nb_sample_u() as f64;
        let mut prev_par_u = cur_par_u;
        let mut local_diff_v =
            (range_v.last - range_v.first) / surface_data.nb_sample_v() as f64;

        let (b_allow_c, b_allow_u, b_allow_v) = check_sampling(
            &curve_range,
            &surface_range,
            curve_data,
            surface_data,
            local_diff_c,
            local_diff_u,
            local_diff_v,
        );
        if !b_allow_c && !b_allow_u && !b_allow_v {
            list_curve_range.push(curve_range);
            list_surface_range.push(surface_range);
            return true;
        }

        let mut new_range_c_template = CurveRangeSample::new();
        if !b_allow_c {
            new_range_c_template = curve_range;
            cur_index_init = curve_range.range_index();
            local_diff_c = range_c.last - range_c.first;
        } else {
            new_range_c_template.set_depth(curve_range.depth() + 1);
            new_range_c_template.set_range_index(cur_index_init);
        }

        let mut new_range_s_template = surface_range;
        if b_allow_u {
            new_range_s_template.set_depth_u(surface_range.depth_u() + 1);
        } else {
            cur_index_u = new_range_s_template.index_u();
            local_diff_u = range_u.last - range_u.first;
        }
        if b_allow_v {
            new_range_s_template.set_depth_v(surface_range.depth_v() + 1);
        } else {
            cur_index_v_init = surface_range.index_v();
            local_diff_v = range_v.last - range_v.first;
        }

        let mut b_has_out = false;
        let nb_u = if b_allow_u {
            surface_data.nb_sample_u()
        } else {
            1
        };
        let nb_v = if b_allow_v {
            surface_data.nb_sample_v()
        } else {
            1
        };
        let nb_c = if b_allow_c {
            curve_data.nb_sample()
        } else {
            1
        };

        for _u_it in 1..=nb_u {
            prev_par_u = cur_par_u;
            cur_par_u += local_diff_u;
            let mut cur_par_v = range_v.first;
            let mut prev_par_v = cur_par_v;
            let mut cur_index_v = cur_index_v_init;
            let mut b_has_out_v = false;

            for _v_it in 1..=nb_v {
                prev_par_v = cur_par_v;
                cur_par_v += local_diff_v;
                let mut new_range_s = new_range_s_template;
                if b_allow_u {
                    new_range_s.set_index_u(cur_index_u);
                }
                if b_allow_v {
                    new_range_s.set_index_v(cur_index_v);
                }
                cur_index_v += 1;

                if surface_data.is_range_out(&new_range_s) {
                    b_has_out_v = true;
                    continue;
                }

                let a_box_s = match surface_data.find_box(&new_range_s) {
                    Some(b) => b,
                    None => {
                        let mut box_s = BndBox::new();
                        let is_bspline = matches!(
                            mesh_surface_type_mesh(self.surface()),
                            MeshSurfaceType::BSplineSurface
                        );
                        if is_bspline && surface_data.has_grid() {
                            box_s = get_surface_box(
                                self.surface(),
                                prev_par_u,
                                cur_par_u,
                                prev_par_v,
                                cur_par_v,
                                self.criteria,
                                surface_data,
                            );
                        } else {
                            add_surface_to_box(
                                self.surface(),
                                prev_par_u,
                                cur_par_u,
                                prev_par_v,
                                cur_par_v,
                                self.criteria,
                                &mut box_s,
                            );
                        }
                        if !b_main_box_found_c {
                            if let Some(b) = curve_data.find_box(&root_range_c) {
                                main_box_c = b;
                                b_main_box_found_c = true;
                            }
                        }
                        if box_s.is_out_box(&main_box_c) {
                            surface_data.add_out_range(new_range_s);
                            b_has_out_v = true;
                            continue;
                        }
                        surface_data.add_box(new_range_s, box_s);
                        box_s
                    }
                };

                if a_box_s.is_out_box(&box_curve) {
                    b_has_out_v = true;
                    continue;
                }

                let mut list_of_box: Vec<BndBox> = Vec::new();
                let mut list_of_index: Vec<i32> = Vec::new();
                let mut b_has_out_c = false;
                let mut cur_par = range_c.first;
                let mut prev_par = range_c.first;
                let mut cur_range_c = new_range_c_template;
                let mut cur_index = cur_index_init;

                for t_it in 1..=nb_c {
                    prev_par = cur_par;
                    cur_par += local_diff_c;
                    let mut b_found = list_c_to_avoid.contains(&t_it);
                    if !b_found {
                        if b_allow_c {
                            cur_range_c.set_range_index(cur_index);
                        }
                        b_found = curve_data.is_range_out(&cur_range_c);
                    }
                    cur_index += 1;
                    if b_found {
                        b_has_out_c = true;
                        continue;
                    }

                    let a_box_c = match curve_data.find_box(&cur_range_c) {
                        Some(b) => b,
                        None => {
                            let mut box_c = BndBox::new();
                            add_curve_to_box(
                                self.curve(),
                                prev_par,
                                cur_par,
                                self.criteria,
                                &mut box_c,
                            );
                            if !b_main_box_found_s {
                                if let Some(b) = surface_data.find_box(&root_range_s) {
                                    main_box_s = b;
                                    b_main_box_found_s = true;
                                }
                            }
                            if box_c.is_out_box(&main_box_s) {
                                curve_data.add_out_range(cur_range_c);
                                b_has_out_c = true;
                                continue;
                            }
                            curve_data.add_box(cur_range_c, box_c);
                            box_c
                        }
                    };

                    if !b_global_check_done && a_box_c.is_out_box(&box_surface) {
                        list_c_to_avoid.push(t_it);
                        b_has_out_c = true;
                        continue;
                    }
                    if a_box_c.is_out_box(&a_box_s) {
                        b_has_out_v = true;
                        b_has_out_c = true;
                        continue;
                    }
                    list_of_index.push(t_it);
                    list_of_box.push(a_box_c);
                }
                b_global_check_done = true;
                if b_has_out_c {
                    b_has_out_v = true;
                }

                let mut new_range_c = new_range_c_template;
                let mut b_use_old_c = false;
                let mut b_use_old_s = false;
                let b_check_size = !b_has_out_c;

                for (t_it, a_box_c) in list_of_index.iter().zip(list_of_box.iter()) {
                    cur_index = cur_index_init + t_it - 1;
                    b_use_old_s = false;
                    if b_allow_c {
                        new_range_c.set_range_index(cur_index);
                    }
                    if b_check_size {
                        if curve_range.depth() == 0
                            || surface_range.depth_u() == 0
                            || surface_range.depth_v() == 0
                        {
                            b_has_out_c = true;
                            b_has_out_v = true;
                        } else if curve_range.depth() < 4
                            && surface_range.depth_u() < 4
                            && surface_range.depth_v() < 4
                        {
                            if !a_box_c.is_whole() && !a_box_s.is_whole() {
                                let diag_c = a_box_c.square_extent();
                                let diag_s = a_box_s.square_extent();
                                if diag_c < diag_s {
                                    if diag_c * 10.0 < diag_s {
                                        b_use_old_c = true;
                                        b_has_out_c = true;
                                        b_has_out_v = true;
                                        break;
                                    }
                                } else if diag_s * 10.0 < diag_c {
                                    b_use_old_s = true;
                                    b_has_out_c = true;
                                    b_has_out_v = true;
                                }
                            }
                        }
                    }
                    if !b_has_out_c {
                        list_curve_found.push(new_range_c);
                        list_surface_found.push(new_range_s);
                    } else {
                        if b_use_old_s && new_range_c == curve_range {
                            return false;
                        }
                        let surf_range = if b_use_old_s {
                            surface_range
                        } else {
                            new_range_s
                        };
                        let surf_box = if b_use_old_s {
                            box_surface
                        } else {
                            a_box_s
                        };
                        if !self.localize_solutions(
                            new_range_c,
                            *a_box_c,
                            surf_range,
                            surf_box,
                            curve_data,
                            surface_data,
                            list_curve_range,
                            list_surface_range,
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
                        curve_range,
                        box_curve,
                        new_range_s,
                        a_box_s,
                        curve_data,
                        surface_data,
                        list_curve_range,
                        list_surface_range,
                    )
                {
                    return false;
                }
                let _ = a_box_s;
            }
            cur_index_u += 1;
            if b_has_out_v {
                b_has_out = true;
            }
        }

        if !b_has_out {
            list_curve_range.push(curve_range);
            list_surface_range.push(surface_range);
        } else {
            for (c, s) in list_curve_found.iter().zip(list_surface_found.iter()) {
                list_curve_range.push(*c);
                list_surface_range.push(*s);
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
