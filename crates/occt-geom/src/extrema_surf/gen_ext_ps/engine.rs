use super::super::prelude::*;
use super::super::{ExtremaExtAlgo, ExtremaExtFlag};
use super::func_ps_norm::ExtremaFuncPsNorm;
use super::helpers::{correct_nb_samples, distinct_knots, fill_params};
use super::point_params::{ExtremaElementType, ExtremaPOnSurfParams};

use occt_core::bnd::BndSphere;
use occt_math::{MathFunctionSetRoot, MathVector};

// ---------------------------------------------------------------------------
// Extrema_GenExtPS
// ---------------------------------------------------------------------------

/// Extrema_GenExtPS (Extrema_GenExtPS.hxx:36-179).
pub struct ExtremaGenExtPs<'a> {
    /// myDone.
    done: bool,
    /// myInit.
    init: bool,
    /// myumin / myusup / myvmin / myvsup.
    umin: f64,
    usup: f64,
    vmin: f64,
    vsup: f64,
    /// myusample / myvsample.
    usample: i32,
    vsample: i32,
    /// mytolu / mytolv.
    tolu: f64,
    tolv: f64,
    /// myPoints (cxx:166), bounds 0..=usample+1 x 0..=vsample+1.
    points: Vec<Vec<ExtremaPOnSurfParams>>,
    /// mySphereArray (cxx:168); index 0..usample*vsample.
    sphere_array: Vec<BndSphere>,
    /// The (NoU, NoV) each mySphereArray entry was built for - the U/V carried
    /// by Bnd_Sphere (cxx:924, Bnd_Sphere.hxx), which the port's BndSphere does
    /// not store.
    sphere_uv: Vec<(i32, i32)>,
    /// !mySphereUBTree.IsNull() (cxx:859). The tree itself is an acceleration
    /// structure; see the module header.
    tree_built: bool,
    /// myF (cxx:169).
    f: ExtremaFuncPsNorm<'a>,
    /// myS (cxx:170).
    s: Option<&'a dyn Surface>,
    /// myFlag / myAlgo.
    flag: ExtremaExtFlag,
    algo: ExtremaExtAlgo,
    /// myUParams / myVParams (cxx:173-174). Vec index i is OCCT index i + 1.
    u_params: Option<Vec<f64>>,
    v_params: Option<Vec<f64>>,
    /// myFacePntParams (cxx:175), bounds 0..=usample x 0..=vsample.
    face_pnt_params: Vec<Vec<ExtremaPOnSurfParams>>,
    /// myUEdgePntParams (cxx:176), OCCT bounds 1..=usample-1 x 1..=vsample.
    u_edge_pnt_params: Vec<Vec<ExtremaPOnSurfParams>>,
    /// myVEdgePntParams (cxx:177), OCCT bounds 1..=usample x 1..=vsample-1.
    v_edge_pnt_params: Vec<Vec<ExtremaPOnSurfParams>>,
    /// myGridParam (cxx:178).
    grid_param: ExtremaPOnSurfParams,
}

impl<'a> Default for ExtremaGenExtPs<'a> {
    fn default() -> Self {
        Self::new()
    }
}

impl<'a> ExtremaGenExtPs<'a> {
    /// Extrema_GenExtPS() (cxx:216-231).
    pub fn new() -> Self {
        ExtremaGenExtPs {
            done: false,
            init: false,
            umin: 0.0,
            usup: 0.0,
            vmin: 0.0,
            vsup: 0.0,
            usample: 0,
            vsample: 0,
            tolu: 0.0,
            tolv: 0.0,
            points: Vec::new(),
            sphere_array: Vec::new(),
            sphere_uv: Vec::new(),
            tree_built: false,
            f: ExtremaFuncPsNorm::default(),
            s: None,
            flag: ExtremaExtFlag::MinMax,
            algo: ExtremaExtAlgo::Grad,
            u_params: None,
            v_params: None,
            face_pnt_params: Vec::new(),
            u_edge_pnt_params: Vec::new(),
            v_edge_pnt_params: Vec::new(),
            grid_param: ExtremaPOnSurfParams::new_default(),
        }
    }

    /// Extrema_GenExtPS(P, S, NbU, NbV, TolU, TolV, F, A) (cxx:239-253):
    /// Initialize then Perform.
    #[allow(clippy::too_many_arguments)]
    pub fn with_surface(
        p: &GpPnt,
        s: &'a dyn Surface,
        nb_u: i32,
        nb_v: i32,
        tol_u: f64,
        tol_v: f64,
        flag: ExtremaExtFlag,
        algo: ExtremaExtAlgo,
    ) -> Self {
        let mut e = ExtremaGenExtPs::new();
        e.flag = flag;
        e.algo = algo;
        e.initialize(s, nb_u, nb_v, tol_u, tol_v);
        e.perform(p);
        e
    }

    /// Extrema_GenExtPS(P, S, NbU, NbV, Umin, Usup, Vmin, Vsup, TolU, TolV, F, A)
    /// (cxx:255-273).
    #[allow(clippy::too_many_arguments)]
    pub fn with_window(
        p: &GpPnt,
        s: &'a dyn Surface,
        nb_u: i32,
        nb_v: i32,
        umin: f64,
        usup: f64,
        vmin: f64,
        vsup: f64,
        tol_u: f64,
        tol_v: f64,
        flag: ExtremaExtFlag,
        algo: ExtremaExtAlgo,
    ) -> Self {
        let mut e = ExtremaGenExtPs::new();
        e.flag = flag;
        e.algo = algo;
        e.initialize_window(s, nb_u, nb_v, umin, usup, vmin, vsup, tol_u, tol_v);
        e.perform(p);
        e
    }

    /// Extrema_GenExtPS::Initialize(S, NbU, NbV, TolU, TolV) (cxx:275-286).
    pub fn initialize(&mut self, s: &'a dyn Surface, nb_u: i32, nb_v: i32, tol_u: f64, tol_v: f64) {
        let (umin, usup) = s.u_range();
        let (vmin, vsup) = s.v_range();
        self.initialize_window(s, nb_u, nb_v, umin, usup, vmin, vsup, tol_u, tol_v);
    }

    /// Extrema_GenExtPS::Initialize(S, NbU, NbV, Umin, Usup, Vmin, Vsup, TolU,
    /// TolV) (cxx:288-319).
    #[allow(clippy::too_many_arguments)]
    pub fn initialize_window(
        &mut self,
        s: &'a dyn Surface,
        nb_u: i32,
        nb_v: i32,
        umin: f64,
        usup: f64,
        vmin: f64,
        vsup: f64,
        tol_u: f64,
        tol_v: f64,
    ) {
        self.s = Some(s);
        self.usample = nb_u;
        self.vsample = nb_v;
        self.tolu = tol_u;
        self.tolv = tol_v;
        self.umin = umin;
        self.usup = usup;
        self.vmin = vmin;
        self.vsup = vsup;

        if self.usample < 2 || self.vsample < 2 {
            // cxx:308-311 throws Standard_OutOfRange.
            panic!("Extrema_GenExtPS::Initialize: nb samples < 2");
        }

        self.f.initialize(s);

        // cxx:315-317: mySphereUBTree.Nullify(); myUParams.Nullify();
        // myVParams.Nullify();
        self.tree_built = false;
        self.sphere_array.clear();
        self.sphere_uv.clear();
        self.u_params = None;
        self.v_params = None;
        self.init = false;
    }

    /// Extrema_GenExtPS::GetGridPoints (cxx:371-454).
    fn grid_points(&mut self, s: &dyn Surface) {
        if s.is_offset_surface() {
            // cxx:375-378.
            if let Some(basis) = s.offset_basis_surface() {
                self.grid_points(&*basis);
            }
        } else if s.is_bspline_surface() {
            // cxx:380-390.
            if let (Some(uk), Some(vk)) = (s.bspline_surface_uknots(), s.bspline_surface_vknots()) {
                let ud = distinct_knots(uk);
                let vd = distinct_knots(vk);
                self.u_params = fill_params(&ud, s.u_degree(), self.umin, self.usup, self.usample);
                self.v_params = fill_params(&vd, s.v_degree(), self.vmin, self.vsup, self.vsample);
            }
        } else if s.is_bezier_surface() {
            // cxx:392-405.
            let (u1, u2) = s.u_range();
            let (v1, v2) = s.v_range();
            self.u_params = fill_params(&[u1, u2], s.u_degree(), self.umin, self.usup, self.usample);
            self.v_params = fill_params(&[v1, v2], s.v_degree(), self.vmin, self.vsup, self.vsample);
        } else if s.is_surface_of_revolution() || s.is_surface_of_linear_extrusion() {
            // cxx:407-444: parametric points follow the generating curve when it
            // is a B-spline or a Bezier.
            let basis = if s.is_surface_of_revolution() {
                s.revolution_basis_curve()
            } else {
                s.extrusion_basis_curve()
            };
            let Some(curve) = basis else {
                return;
            };
            let mut arr: Option<Vec<f64>> = None;
            let mut degree: i32 = 0;
            if let Some(knots) = curve.bspline_knots() {
                // BasisCurve()->GetType() == GeomAbs_BSplineCurve (cxx:412-420).
                arr = Some(distinct_knots(knots));
                degree = curve.nurbs_degree().unwrap_or(1) as i32;
            }
            if arr.is_none() && curve.bezier_poles().is_some() {
                // BasisCurve()->GetType() == GeomAbs_BezierCurve (cxx:421-431).
                arr = Some(vec![curve.first_parameter(), curve.last_parameter()]);
                degree = curve.nurbs_degree().unwrap_or(1) as i32;
            }
            let Some(arr) = arr else {
                return;
            };
            if s.is_surface_of_revolution() {
                self.v_params = fill_params(&arr, degree, self.vmin, self.vsup, self.vsample);
            } else {
                self.u_params = fill_params(&arr, degree, self.umin, self.usup, self.usample);
            }
        }
        // cxx:445-453: update the number of points in sample.
        if let Some(u) = &self.u_params {
            self.usample = u.len() as i32;
        }
        if let Some(v) = &self.v_params {
            self.vsample = v.len() as i32;
        }
    }

    /// Extrema_GenExtPS::ComputeEdgeParameters (cxx:460-526). OCCT returns a
    /// reference to myGridParam (or to one of the inputs); the port returns the
    /// value, which is what every caller stores.
    fn compute_edge_parameters(
        &mut self,
        is_u_edge: bool,
        param0: &ExtremaPOnSurfParams,
        param1: &ExtremaPOnSurfParams,
        point: &GpPnt,
        diff_tol: f64,
    ) -> ExtremaPOnSurfParams {
        let s = self.s.expect("Extrema_GenExtPS: surface not set");
        let sqr_dist01 = param0.value().square_distance(&param1.value());

        if sqr_dist01 <= diff_tol {
            // cxx:469-473: the points are confused.
            return *param0;
        }
        let diff_dist = (param0.sqr_distance() - param1.sqr_distance()).abs();
        if diff_dist >= sqr_dist01 - diff_tol {
            // cxx:476-491: the shortest distance is one of the nodes.
            return if param0.sqr_distance() > param1.sqr_distance() {
                *param1
            } else {
                *param0
            };
        }
        // cxx:493-524: the shortest distance is inside the edge.
        let pop = GpVec::from_pnts(&param0.value(), point);
        let pop1 = GpVec::from_pnts(&param0.value(), &param1.value());
        let ratio = pop.dot(&pop1) / sqr_dist01;
        let (au0, av0) = param0.parameter();
        let (au1, av1) = param1.parameter();
        let (mut u_par, mut v_par) = (au0, av0);
        if is_u_edge {
            u_par += ratio * (au1 - au0);
        } else {
            v_par += ratio * (av1 - av0);
        }
        self.grid_param.set_parameters(u_par, v_par, s.d0(u_par, v_par));
        let (i0, j0) = param0.indices();
        self.grid_param.set_element_type(if is_u_edge {
            ExtremaElementType::UIsoEdge
        } else {
            ExtremaElementType::VIsoEdge
        });
        self.grid_param
            .set_sqr_distance(point.square_distance(&self.grid_param.value()));
        self.grid_param.set_indices(i0, j0);
        self.grid_param
    }

    /// Extrema_GenExtPS::BuildGrid (cxx:528-753).
    fn build_grid(&mut self, point: &GpPnt) {
        let s = self.s.expect("Extrema_GenExtPS: surface not set");

        if !self.init {
            // cxx:534: build parametric grid for complex geometry.
            self.grid_points(s);
            let usample = self.usample as usize;
            let vsample = self.vsample as usize;

            // cxx:537-549: build grid in other cases (U).
            if self.u_params.is_none() {
                let mut pas_u = self.usup - self.umin;
                let mut u0 = pas_u / self.usample as f64 / 100.0;
                pas_u = (pas_u - u0) / (self.usample as f64 - 1.0);
                u0 = u0 / 2.0 + self.umin;
                let mut params = Vec::with_capacity(usample);
                let mut u = u0;
                for _ in 1..=usample {
                    params.push(u);
                    u += pas_u;
                }
                self.u_params = Some(params);
            }
            // cxx:551-564: build grid in other cases (V).
            if self.v_params.is_none() {
                let mut pas_v = self.vsup - self.vmin;
                let mut v0 = pas_v / self.vsample as f64 / 100.0;
                pas_v = (pas_v - v0) / (self.vsample as f64 - 1.0);
                v0 = v0 / 2.0 + self.vmin;
                let mut params = Vec::with_capacity(vsample);
                let mut v = v0;
                for _ in 1..=vsample {
                    params.push(v);
                    v += pas_v;
                }
                self.v_params = Some(params);
            }

            // cxx:566-567: myPoints.Resize(0, usample+1, 0, vsample+1, false).
            let nu = self.usample as usize + 2;
            let nv = self.vsample as usize + 2;
            self.points = vec![vec![ExtremaPOnSurfParams::new_default(); nv]; nu];

            // cxx:569-572: GeomGridEval_Surface::EvaluateGrid over the two
            // parameter arrays - the same points as Value(U_i, V_j).
            let u_params = self.u_params.clone().expect("u params");
            let v_params = self.v_params.clone().expect("v params");
            for no_u in 1..=self.usample as usize {
                for no_v in 1..=self.vsample as usize {
                    let p1 = s.d0(u_params[no_u - 1], v_params[no_v - 1]);
                    let mut param =
                        ExtremaPOnSurfParams::new(u_params[no_u - 1], v_params[no_v - 1], p1);
                    param.set_element_type(ExtremaElementType::Node);
                    param.set_indices(no_u as i32, no_v as i32);
                    self.points[no_u][no_v] = param;
                }
            }

            // cxx:591-593.
            let fu = self.usample as usize + 1;
            let fv = self.vsample as usize + 1;
            self.face_pnt_params = vec![vec![ExtremaPOnSurfParams::new_default(); fv]; fu];
            self.u_edge_pnt_params = vec![vec![ExtremaPOnSurfParams::new_default(); fv]; fu];
            self.v_edge_pnt_params = vec![vec![ExtremaPOnSurfParams::new_default(); fv]; fu];

            // cxx:595-607: fill the boundary with negative square distance (used
            // for the maximum search).
            for no_v in 0..=self.vsample as usize {
                self.points[0][no_v].set_sqr_distance(-1.0);
                self.points[self.usample as usize + 1][no_v].set_sqr_distance(-1.0);
            }
            for no_u in 1..=self.usample as usize {
                self.points[no_u][0].set_sqr_distance(-1.0);
                self.points[no_u][self.vsample as usize + 1].set_sqr_distance(-1.0);
            }

            self.init = true;
        }

        let usample = self.usample as usize;
        let vsample = self.vsample as usize;

        // cxx:612-622: step 1, distances to the nodes.
        for no_u in 1..=usample {
            for no_v in 1..=vsample {
                let d = point.square_distance(&self.points[no_u][no_v].value());
                self.points[no_u][no_v].set_sqr_distance(d);
            }
        }

        // cxx:624-752: for the minimum, compute distances to the mesh (edges
        // then faces).
        if self.flag == ExtremaExtFlag::Min || self.flag == ExtremaExtFlag::MinMax {
            let diff_tol = self.tolu + self.tolv;

            // cxx:631-660: step 2, distances to the edges.
            for no_u in 1..=usample {
                for no_v in 1..=vsample {
                    let param0 = self.points[no_u][no_v];
                    if no_u < usample {
                        let param1 = self.points[no_u + 1][no_v];
                        let edge =
                            self.compute_edge_parameters(true, &param0, &param1, point, diff_tol);
                        self.u_edge_pnt_params[no_u][no_v] = edge;
                    }
                    if no_v < vsample {
                        let param1 = self.points[no_u][no_v + 1];
                        let edge = self.compute_edge_parameters(
                            false, &param0, &param1, point, diff_tol,
                        );
                        self.v_edge_pnt_params[no_u][no_v] = edge;
                    }
                }
            }

            // cxx:662-738: step 3, distances to the faces.
            for no_u in 1..usample {
                for no_v in 1..vsample {
                    let ue0 = self.u_edge_pnt_params[no_u][no_v];
                    let ue1 = self.u_edge_pnt_params[no_u][no_v + 1];
                    let ve0 = self.v_edge_pnt_params[no_u][no_v];
                    let ve1 = self.v_edge_pnt_params[no_u + 1][no_v];

                    let mut sqr_dist01 = ue0.value().square_distance(&ue1.value());
                    let mut diff_dist = (ue0.sqr_distance() - ue1.sqr_distance()).abs();
                    let mut is_out = false;
                    if diff_dist >= sqr_dist01 - diff_tol {
                        // cxx:684-688: projection outside the face.
                        is_out = true;
                    } else {
                        sqr_dist01 = ve0.value().square_distance(&ve1.value());
                        diff_dist = (ve0.sqr_distance() - ve1.sqr_distance()).abs();
                        if diff_dist >= sqr_dist01 - diff_tol {
                            is_out = true;
                        }
                    }

                    if is_out {
                        // cxx:701-712: closest point on an edge.
                        let ue_min = if ue0.sqr_distance() < ue1.sqr_distance() {
                            ue0
                        } else {
                            ue1
                        };
                        let ve_min = if ve0.sqr_distance() < ve1.sqr_distance() {
                            ve0
                        } else {
                            ve1
                        };
                        let e_min = if ue_min.sqr_distance() < ve_min.sqr_distance() {
                            ue_min
                        } else {
                            ve_min
                        };
                        self.face_pnt_params[no_u][no_v] = e_min;
                    } else {
                        // cxx:713-736: closest point inside the face.
                        let (au0, _av0) = ue0.parameter();
                        let (au1, _av1) = ue1.parameter();
                        let u_par = 0.5 * (au0 + au1);
                        let (_bu0, bv0) = ve0.parameter();
                        let (_bu1, bv1) = ve1.parameter();
                        let v_par = 0.5 * (bv0 + bv1);

                        let mut param =
                            ExtremaPOnSurfParams::new(u_par, v_par, s.d0(u_par, v_par));
                        param.set_element_type(ExtremaElementType::Face);
                        param.set_sqr_distance(point.square_distance(&param.value()));
                        param.set_indices(no_u as i32, no_v as i32);
                        self.face_pnt_params[no_u][no_v] = param;
                    }
                }
            }

            // cxx:740-751: fill the boundary with RealLast.
            for no_v in 0..=self.vsample as usize {
                self.face_pnt_params[0][no_v].set_sqr_distance(f64::MAX);
                self.face_pnt_params[self.usample as usize][no_v].set_sqr_distance(f64::MAX);
            }
            for no_u in 1..self.usample as usize {
                self.face_pnt_params[no_u][0].set_sqr_distance(f64::MAX);
                self.face_pnt_params[no_u][self.vsample as usize].set_sqr_distance(f64::MAX);
            }
        }
    }

    /// Extrema_GenExtPS::BuildTree (cxx:856-931). The grid spheres are the
    /// selection input (mySphereArray, cxx:911-929); the NCollection_UBTree is
    /// not reproduced - see the module header.
    fn build_tree(&mut self) {
        if self.tree_built {
            // cxx:858-862.
            return;
        }
        let s = self.s.expect("Extrema_GenExtPS: surface not set");

        // cxx:864-878: raise the sampling for B-spline surfaces.
        if s.is_bspline_surface() {
            let uk = s
                .bspline_surface_uknots()
                .map(distinct_knots)
                .unwrap_or_default();
            let vk = s
                .bspline_surface_vknots()
                .map(distinct_knots)
                .unwrap_or_default();
            let u_value = s.u_degree() * uk.len() as i32;
            let v_value = s.v_degree() * vk.len() as i32;
            if u_value > self.usample {
                self.usample = u_value.min(300);
            }
            if v_value > self.vsample {
                self.vsample = v_value.min(300);
            }
        }
        // cxx:880.
        correct_nb_samples(
            s,
            self.umin,
            self.usup,
            &mut self.usample,
            self.vmin,
            self.vsup,
            &mut self.vsample,
        );

        // cxx:882-903: uniform parametric grid.
        let mut pas_u = self.usup - self.umin;
        let mut pas_v = self.vsup - self.vmin;
        let mut u0 = pas_u / self.usample as f64 / 100.0;
        let mut v0 = pas_v / self.vsample as f64 / 100.0;
        pas_u = (pas_u - u0) / (self.usample as f64 - 1.0);
        pas_v = (pas_v - v0) / (self.vsample as f64 - 1.0);
        u0 = u0 / 2.0 + self.umin;
        v0 = v0 / 2.0 + self.vmin;

        let mut u_params = Vec::with_capacity(self.usample as usize);
        let mut v_params = Vec::with_capacity(self.vsample as usize);
        let mut u = u0;
        for _ in 1..=self.usample {
            u_params.push(u);
            u += pas_u;
        }
        let mut v = v0;
        for _ in 1..=self.vsample {
            v_params.push(v);
            v += pas_v;
        }
        self.u_params = Some(u_params);
        self.v_params = Some(v_params);

        // cxx:905-929: build the sphere array.
        self.sphere_array.clear();
        self.sphere_uv.clear();
        let u_params = self.u_params.clone().unwrap();
        let v_params = self.v_params.clone().unwrap();
        for no_u in 1..=self.usample as usize {
            for no_v in 1..=self.vsample as usize {
                let p1 = s.d0(u_params[no_u - 1], v_params[no_v - 1]);
                // Bnd_Sphere(aP1.XYZ(), 0, NoU, NoV) (cxx:924).
                self.sphere_array.push(BndSphere::from_center_radius(p1, 0.0));
                self.sphere_uv.push((no_u as i32, no_v as i32));
            }
        }
        // cxx:930: aFiller.Fill().
        self.tree_built = true;
    }

    /// Extrema_GenExtPS::FindSolution (cxx:933-952).
    fn find_solution(&mut self, _p: &GpPnt, params: &ExtremaPOnSurfParams) {
        let mut tol = MathVector::new(1, 2);
        tol.set_value(1, self.tolu);
        tol.set_value(2, self.tolv);

        let mut uv = MathVector::new(1, 2);
        let (u, v) = params.parameter();
        uv.set_value(1, u);
        uv.set_value(2, v);

        let mut uv_inf = MathVector::new(1, 2);
        let mut uv_sup = MathVector::new(1, 2);
        uv_inf.set_value(1, self.umin);
        uv_inf.set_value(2, self.vmin);
        uv_sup.set_value(1, self.usup);
        uv_sup.set_value(2, self.vsup);

        // math_FunctionSetRoot S(myF, Tol); S.Perform(myF, UV, UVinf, UVsup);
        let mut solver = MathFunctionSetRoot::new(&self.f, &tol, 100);
        solver.perform_with_bounds(&mut self.f, &uv, &uv_inf, &uv_sup, false);

        self.done = true;
    }

    /// Extrema_GenExtPS::SetFlag (cxx:954-957).
    pub fn set_flag(&mut self, f: ExtremaExtFlag) {
        self.flag = f;
    }

    /// Extrema_GenExtPS::SetAlgo (cxx:959-966).
    pub fn set_algo(&mut self, a: ExtremaExtAlgo) {
        if self.algo != a {
            self.init = false;
        }
        self.algo = a;
    }

    /// The flag set by [ExtremaGenExtPs::set_flag].
    pub fn flag(&self) -> ExtremaExtFlag {
        self.flag
    }

    /// The algorithm set by [ExtremaGenExtPs::set_algo].
    pub fn algo(&self) -> ExtremaExtAlgo {
        self.algo
    }

    /// Extrema_GenExtPS::Perform (cxx:968-1151).
    pub fn perform(&mut self, p: &GpPnt) {
        self.done = false;
        self.f.set_point(*p);

        if self.algo == ExtremaExtAlgo::Grad {
            self.build_grid(p);
            let usample = self.usample;
            let vsample = self.vsample;

            if self.flag == ExtremaExtFlag::Min || self.flag == ExtremaExtFlag::MinMax {
                // cxx:978-1079 - minimums over the faces.
                for no_u in 1..usample as usize {
                    for no_v in 1..vsample as usize {
                        let param = self.face_pnt_params[no_u][no_v];
                        let mut is_min = false;
                        let an_elem_type = param.element_type();

                        if an_elem_type == ExtremaElementType::Face {
                            is_min = true;
                        } else {
                            let (i_u, i_v) = param.indices();
                            if an_elem_type == ExtremaElementType::UIsoEdge {
                                is_min = i_v == 1 || i_v == vsample;
                            } else if an_elem_type == ExtremaElementType::VIsoEdge {
                                is_min = i_u == 1 || i_u == usample;
                            } else if an_elem_type == ExtremaElementType::Node {
                                is_min =
                                    (i_u == 1 || i_u == usample) && (i_v == 1 || i_v == vsample);
                            }

                            if !is_min {
                                if an_elem_type == ExtremaElementType::UIsoEdge
                                    || (an_elem_type == ExtremaElementType::Node
                                        && (i_u == 1 || i_u == usample))
                                {
                                    // cxx:1025-1033: check the down face.
                                    let down = self.face_pnt_params[no_u][no_v - 1];
                                    if down.element_type() == an_elem_type {
                                        let (i_u2, i_v2) = down.indices();
                                        is_min = i_u == i_u2 && i_v == i_v2;
                                    }
                                } else if an_elem_type == ExtremaElementType::VIsoEdge
                                    || (an_elem_type == ExtremaElementType::Node
                                        && (i_v == 1 || i_v == vsample))
                                {
                                    // cxx:1034-1045: check the right face.
                                    let right = self.face_pnt_params[no_u - 1][no_v];
                                    if right.element_type() == an_elem_type {
                                        let (i_u2, i_v2) = right.indices();
                                        is_min = i_u == i_u2 && i_v == i_v2;
                                    }
                                } else if i_u == no_u as i32 && i_v == no_v as i32 {
                                    // cxx:1046-1069: the lower-left node.
                                    is_min = true;
                                    let others = [
                                        self.face_pnt_params[no_u][no_v - 1],
                                        self.face_pnt_params[no_u - 1][no_v - 1],
                                        self.face_pnt_params[no_u - 1][no_v],
                                    ];
                                    for other in others.iter() {
                                        if !is_min {
                                            break;
                                        }
                                        if other.element_type() == ExtremaElementType::Node {
                                            let (i_u2, i_v2) = other.indices();
                                            is_min = i_u == i_u2 && i_v == i_v2;
                                        } else {
                                            is_min = false;
                                        }
                                    }
                                }
                            }
                        }

                        if is_min {
                            self.find_solution(p, &param);
                        }
                    }
                }
            }

            if self.flag == ExtremaExtFlag::Max || self.flag == ExtremaExtFlag::MinMax {
                // cxx:1081-1111 - maximums over the 3x3 node neighborhood.
                for no_u in 1..=usample as usize {
                    for no_v in 1..=vsample as usize {
                        let main = self.points[no_u][no_v];
                        let p1 = self.points[no_u - 1][no_v - 1];
                        let p2 = self.points[no_u - 1][no_v];
                        let p3 = self.points[no_u - 1][no_v + 1];
                        let p4 = self.points[no_u][no_v - 1];
                        let p5 = self.points[no_u][no_v + 1];
                        let p6 = self.points[no_u + 1][no_v - 1];
                        let p7 = self.points[no_u + 1][no_v];
                        let p8 = self.points[no_u + 1][no_v + 1];

                        let dist = main.sqr_distance();
                        if p1.sqr_distance() <= dist
                            && p2.sqr_distance() <= dist
                            && p3.sqr_distance() <= dist
                            && p4.sqr_distance() <= dist
                            && p5.sqr_distance() <= dist
                            && p6.sqr_distance() <= dist
                            && p7.sqr_distance() <= dist
                            && p8.sqr_distance() <= dist
                        {
                            self.find_solution(p, &main);
                        }
                    }
                }
            }
        } else {
            // cxx:1113-1150 - the tree algorithm.
            self.build_tree();

            if self.flag == ExtremaExtFlag::Min || self.flag == ExtremaExtFlag::MinMax {
                let idx = self.select_sphere(p, true);
                let (nu, nv) = self.sphere_uv[idx];
                let u = self.u_params.as_ref().unwrap()[(nu - 1) as usize];
                let v = self.v_params.as_ref().unwrap()[(nv - 1) as usize];
                let s = self.s.expect("Extrema_GenExtPS: surface not set");
                let mut params = ExtremaPOnSurfParams::new(u, v, s.d0(u, v));
                params.set_sqr_distance(p.square_distance(&params.value()));
                params.set_indices(nu, nv);
                self.find_solution(p, &params);
            }
            if self.flag == ExtremaExtFlag::Max || self.flag == ExtremaExtFlag::MinMax {
                let idx = self.select_sphere(p, false);
                let (nu, nv) = self.sphere_uv[idx];
                let u = self.u_params.as_ref().unwrap()[(nu - 1) as usize];
                let v = self.v_params.as_ref().unwrap()[(nv - 1) as usize];
                let s = self.s.expect("Extrema_GenExtPS: surface not set");
                let mut params = ExtremaPOnSurfParams::new(u, v, s.d0(u, v));
                params.set_sqr_distance(p.square_distance(&params.value()));
                params.set_indices(nu, nv);
                self.find_solution(p, &params);
            }
        }
    }

    /// Bnd_SphereUBTreeSelectorMin/Max::Accept (cxx:104-169) over the whole
    /// mySphereArray: select the sphere minimizing (min == true) or maximizing
    /// Bnd_Sphere::Distance(theCheckPoint). The UBTree only prunes candidates, so
    /// the selected sphere is the same.
    fn select_sphere(&self, p: &GpPnt, min: bool) -> usize {
        let mut best = 0usize;
        let mut best_dist = self.sphere_array[0].distance(p);
        for i in 1..self.sphere_array.len() {
            let d = self.sphere_array[i].distance(p);
            if (min && d < best_dist) || (!min && d > best_dist) {
                best = i;
                best_dist = d;
            }
        }
        best
    }

    /// Extrema_GenExtPS::IsDone (cxx:1155-1158).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// Extrema_GenExtPS::NbExt (cxx:1162-1169). OCCT throws StdFail_NotDone when
    /// the search failed; the port reports the raw count.
    pub fn nb_ext(&self) -> usize {
        self.f.nb_ext()
    }

    /// Extrema_GenExtPS::SquareDistance(N) (cxx:1173-1181).
    pub fn square_distance(&self, n: usize) -> f64 {
        self.f.square_distance(n)
    }

    /// Extrema_GenExtPS::Point(N) (cxx:1185-1193).
    pub fn point(&self, n: usize) -> (f64, f64, GpPnt) {
        let ps = self.f.point(n);
        let (u, v) = ps.parameter();
        (u, v, ps.value())
    }
}
