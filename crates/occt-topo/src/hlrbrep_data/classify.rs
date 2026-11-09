//! Point-in-face classification of one edge against the current hiding face -
//! `Classify` (`HLRBRep_Data.cxx:1993-2270`), `SimplClassify`
//! (`:2272-2328`) and `IsBadFace` (`:2655-2682`), plus the local
//! `REJECT1` helper (`cxx:1946-1991`).
//!
//! `Classify` projects one point of the tested edge, prunes it against the
//! hiding face with the sixteen-direction box test, shoots a ray from the eye
//! through the projected point (`myProj.Shoot`) and counts the crossings of
//! that ray with the face in front of the edge (`w < wLim`). Every crossing is
//! confirmed by `myClassifier->Classify(pnt2d, PConfusion())`, i.e. by the 2D
//! face classifier of `BRepTopAdaptor_TopolTool`, and `Level` counts the
//! crossings unless `LevelFlag` is false (then the first crossing ends the
//! search and the answer is just `IN`).
//!
//! The four `REJECT1` blocks of `Classify` (`cxx:2019-2044`, `:2052-2077`,
//! `:2080-2105`, `:2110-2176`) and the one of `SimplClassify`
//! (`cxx:2294-2320`) run the same chain - `InitMinMax(Precision::Infinite())`,
//! `UpdateMinMax`, `EnlargeMinMax(tol)`, `REJECT1`, `EncodeMinMax`, then the
//! sixteen-direction comparison against `iFaceMinMax` - so they share
//! [`HlrBRepData::rejected_point_z`].
//!
//! Deviation: `myClassifier` is an `occ::handle` in OCCT and the port stores
//! `Arc<TopolTool>`. `Classify` needs `&mut` on the tool (its `FClass2d` is
//! built lazily), so the handle is moved out for the point loop and
//! `Arc::make_mut` clones the tool when the `MST` map still holds a second
//! reference. The clone carries the same face and the same cache contents, so
//! the answer is unchanged; only the lazily-built `FClass2d` may be rebuilt
//! instead of reused. A missing classifier answers `UNKNOWN`, which is
//! `!= TopAbs_OUT` exactly as a null OCCT handle would crash.

use occt_core::gp::GpPnt2d;
use occt_core::precision::{CONFUSION, INFINITE};

use crate::abs::State;
use crate::fclass2d::FaceState;
use crate::hlr_algo::MinMaxIndices;
use crate::hlr_algo::{encode_min_max, enlarge_min_max, init_min_max, update_min_max};

use super::HlrBRepData;
use super::edge_iter::separated;

/// `REJECT1(theDeca, theTotMin, theTotMax, theSurD, theVertMin, theVertMax)`
/// (`cxx:1946-1991`): the sixteen scaled bounds packed into two
/// [`MinMaxIndices`], the first eight entries in `Min`, the last eight in
/// `Max`.
pub(crate) fn reject1(
    the_deca: &[f64; 16],
    the_tot_min: &[f64; 16],
    the_tot_max: &[f64; 16],
    the_sur_d: &[f64; 16],
) -> (MinMaxIndices, MinMaxIndices) {
    let mut vert_min = MinMaxIndices::default();
    let mut vert_max = MinMaxIndices::default();
    for i in 0..8 {
        vert_min.min[i] = ((the_deca[i] + the_tot_min[i]) * the_sur_d[i]) as i32;
        vert_max.min[i] = ((the_deca[i] + the_tot_max[i]) * the_sur_d[i]) as i32;
        vert_min.max[i] = ((the_deca[i + 8] + the_tot_min[i + 8]) * the_sur_d[i + 8]) as i32;
        vert_max.max[i] = ((the_deca[i + 8] + the_tot_max[i + 8]) * the_sur_d[i + 8]) as i32;
    }
    (vert_min, vert_max)
}

impl HlrBRepData {
    /// One `REJECT1` block: project `(x, y, z)` into the sixteen directions,
    /// grow by `tol`, rescale by `myDeca` / `mySurD` and test the result
    /// against `iFaceMinMax`. `true` means "rejected" (`cxx:2019-2044` and the
    /// three sibling blocks of `Classify`, `:2294-2320` of `SimplClassify`).
    fn rejected_point_z(&self, tol: f64, x: f64, y: f64, z: f64) -> bool {
        let mut tot_min = [0.0f64; 16];
        let mut tot_max = [0.0f64; 16];
        init_min_max(INFINITE, &mut tot_min, &mut tot_max);
        update_min_max(x, y, z, &mut tot_min, &mut tot_max);
        enlarge_min_max(tol, &mut tot_min, &mut tot_max);

        let (vert_min, vert_max) = reject1(&self.my_deca, &tot_min, &tot_max, &self.my_sur_d);
        let mut min_max_vert = MinMaxIndices::default();
        encode_min_max(&vert_min, &vert_max, &mut min_max_vert);

        match self.i_face_min_max {
            Some(face) => separated(&face, &min_max_vert),
            // `iFaceMinMax` is a raw pointer in OCCT; a null one would crash.
            None => false,
        }
    }

    /// `Classify(E, ED, LevelFlag, Level, param)` (`cxx:1993-2270`).
    ///
    /// `E` is only `(void)E` in OCCT (`cxx:1999`): the edge under test is `ED`,
    /// which every caller takes from `myEData(E)`, so the port derives it from
    /// the index and never aliases `&mut self`.
    pub fn classify(&mut self, e: i32, level_flag: bool, level: &mut i32, param: f64) -> State {
        let ed = e as usize;
        let tol = self.my_e_data[ed].tolerance() as f64;

        *level = 0;
        let mut state = State::Out;

        let sta;
        if level_flag {
            // `sta = param` (`cxx:2014-2017`).
            sta = param;
            let p = self.my_e_data[ed].geometry().value3d(sta);
            let (x, y, z) = self.my_proj.project_coords(&p);
            if self.rejected_point_z(tol, x, y, z) {
                return state;
            }
        } else {
            // `cxx:2046-2176`: start point, end point, then the 40/60 point of
            // the (unshrunk) parameter range.
            sta = self.my_e_data[ed].geometry().parameter3d(
                self.my_e_data[ed].geometry().first_parameter(),
            );
            let p = self.my_e_data[ed].geometry().value3d(sta);
            let (x, y, z) = self.my_proj.project_coords(&p);
            if self.rejected_point_z(tol, x, y, z) {
                return state;
            }

            let end = self.my_e_data[ed]
                .geometry()
                .parameter3d(self.my_e_data[ed].geometry().last_parameter());
            let p = self.my_e_data[ed].geometry().value3d(end);
            let (x, y, z) = self.my_proj.project_coords(&p);
            if self.rejected_point_z(tol, x, y, z) {
                return state;
            }

            // `sta = 0.4 * sta + 0.6 * end;` (`cxx:2107`).
            let sta_mid = 0.4 * sta + 0.6 * end;
            let p = self.my_e_data[ed].geometry().value3d(sta_mid);
            let (x, y, z) = self.my_proj.project_coords(&p);
            if self.rejected_point_z(tol, x, y, z) {
                return state;
            }
            // The intersection below uses `sta`; OCCT overwrites it in place.
            self.classify_ray(ed, sta_mid, level_flag, level, &mut state);
            return state;
        }

        self.classify_ray(ed, sta, level_flag, level, &mut state);
        state
    }

    /// The `cxx:2180-2267` tail of both `Classify` arms: shoot the ray, load
    /// the face into `myIntersector`, and count the crossings in front of the
    /// edge.
    fn classify_ray(&mut self, ed: usize, sta: f64, level_flag: bool, level: &mut i32, state: &mut State) {
        let psta = self.my_e_data[ed].geometry().value(sta);
        let plim = self.my_e_data[ed].geometry().value3d(sta);

        let l = match self.my_proj.shoot(psta.x(), psta.y()) {
            Ok(l) => l,
            // `HLRAlgo_Projector::Shoot` throws from `gp_Dir` on a null
            // direction (`Projector.cxx:347-358`).
            Err(_) => return,
        };
        let mut w_lim = occt_core::elib::clib::line_parameter(&l, &plim);
        self.my_intersector.perform_line(&l, w_lim);
        if !self.my_intersector.is_done() {
            return;
        }
        let nb_points = self.my_intersector.nb_points();
        if nb_points == 0 {
            return;
        }

        // `TolZ = myBigSize * 0.000001` and its two `iFaceTest` arms
        // (`cxx:2202-2211`).
        let mut tol_z = self.my_big_size * 0.000_001;
        if self.i_face_test {
            if !self.my_le_out_line && !self.my_le_internal {
                tol_z = self.my_big_size * 0.001;
            } else {
                tol_z = self.my_big_size * 0.01;
            }
        }
        w_lim -= tol_z;

        // `cxx:2215-2232`: the two periods, 0.0 when not periodic.
        let (period_u, umin, umax, period_v, vmin, vmax) = match self.i_face_geom {
            Some(idx) => {
                let s = self.my_f_data[idx as usize].geometry();
                let (pu, umin, umax) = if s.is_u_periodic() {
                    (s.u_period(), s.first_u_parameter(), s.last_u_parameter())
                } else {
                    (0.0, 0.0, 0.0)
                };
                let (pv, vmin, vmax) = if s.is_v_periodic() {
                    (s.v_period(), s.first_v_parameter(), s.last_v_parameter())
                } else {
                    (0.0, 0.0, 0.0)
                };
                (pu, umin, umax, pv, vmin, vmax)
            }
            None => return,
        };

        // The handle is moved out so `Arc::make_mut` can hand out `&mut` while
        // `myIntersector` / `myFData` are read (see the module notes).
        let mut classifier = self.my_classifier.take();
        for i in 0..nb_points {
            let cp = self.my_intersector.cs_point(i);
            let (mut u, mut v, w) = (cp.u(), cp.v(), cp.param());
            if w >= w_lim {
                continue;
            }
            if period_u > 0.0 {
                let (nu, _) = crate::geom_int::geom_int_adjust_periodic(u, umin, umax, period_u, 0.0);
                u = nu;
            }
            if period_v > 0.0 {
                let (nv, _) = crate::geom_int::geom_int_adjust_periodic(v, vmin, vmax, period_v, 0.0);
                v = nv;
            }
            let pnt2d = GpPnt2d::new(u, v);
            let st = match classifier.as_mut() {
                Some(arc) => std::sync::Arc::make_mut(arc).classify(pnt2d, CONFUSION, true),
                None => FaceState::Unknown,
            };
            if st != FaceState::Out {
                *state = State::In;
                *level += 1;
                if !level_flag {
                    break;
                }
            }
        }
        self.my_classifier = classifier;
    }

    /// `SimplClassify(E, ED, Nbp, p1, p2)` (`cxx:2272-2328`): `OUT` as soon as
    /// one of the `Nbp` evenly spaced interior points is rejected, `IN`
    /// otherwise. `E` is unused in OCCT too.
    pub fn simpl_classify(&mut self, e: i32, nbp: i32, p1: f64, p2: f64) -> State {
        let ed = e as usize;
        let tol = self.my_e_data[ed].tolerance() as f64;
        let dp = (p2 - p1) / (nbp as f64 + 1.0);

        let mut sta = p1 + dp;
        let mut i = 1;
        while i <= nbp {
            let p = self.my_e_data[ed].geometry().value3d(sta);
            let (x, y, z) = self.my_proj.project_coords(&p);
            if self.rejected_point_z(tol, x, y, z) {
                return State::Out;
            }
            i += 1;
            sta += dp;
        }
        State::In
    }

    /// `IsBadFace()` (`cxx:2655-2682`): a periodic direction whose face bounds
    /// span more than two periods is garbage data.
    pub fn is_bad_face(&self) -> bool {
        let Some(idx) = self.i_face_geom else {
            return false;
        };
        let s = self.my_f_data[idx as usize].geometry();
        if s.is_u_periodic() && s.u_period() * 2.0 < s.last_u_parameter() - s.first_u_parameter() {
            return true;
        }
        if s.is_v_periodic() && s.v_period() * 2.0 < s.last_v_parameter() - s.first_v_parameter() {
            return true;
        }
        false
    }
}
