//! `ApproxInt_ImpPrmSvSurfaces` — the analytic-parametric surface-pair
//! evaluator: one surface is a quadric (`IntSurf_Quadric`), the other a
//! parametric surface (`BRepAdaptor_Surface`), and every point of the
//! intersection line is found by a 2D root search on the parametric surface.
//!
//! Sources:
//! - `ModelingAlgorithms/TKGeomAlgo/ApproxInt/ApproxInt_ImpPrmSvSurfaces.gxx:1-1020`
//!   (`IsSingular` `:32-56`, `SingularProcessing` `:80-267`,
//!   `NonSingularProcessing` `:287-322`, the constructors `:325-352`, the four
//!   tangent accessors `:356-429`, `Compute` `:437-768`, `SeekPoint` `:776-861`,
//!   `FillInitialVectorOfSolution` `:864-1020`)
//! - `ModelingAlgorithms/TKTopAlgo/BRepApprox/BRepApprox_TheImpPrmSvSurfacesOfApprox_0.cxx:16-44`
//!   (the `#define ThePSurface BRepAdaptor_Surface` / `ThePSurfaceTool
//!   BRepApprox_SurfaceTool` / `TheISurface IntSurf_Quadric` / `TheISurfaceTool
//!   IntSurf_QuadricTool` / `TheZerImpFunc …` binding)
//! - `BRepApprox_TheImpPrmSvSurfacesOfApprox.hxx:39-165`
//!
//! The `.gxx` is one class; the port keeps it one class spread over three
//! files of the same module namespace:
//! - this file: the members, the constructors, `Compute` and the
//!   `ApproxInt_SvSurfaces` implementation;
//! - [`imp_prm_tangent_processing`](super::imp_prm_tangent_processing): the
//!   file-static `IsSingular` / `SingularProcessing` / `NonSingularProcessing`;
//! - [`imp_prm_initial_guess`](super::imp_prm_initial_guess): `SeekPoint` and
//!   `FillInitialVectorOfSolution`.
//!
//! The scalar function is [`ZerImpFunc`](crate::brep_approx::ZerImpFunc), the
//! port of `IntImp_ZerImpFunc.gxx` over the same two surfaces; `MyZerImpFunc`
//! keeps the member the OCCT header declares.
//!
//! Translation notes:
//! - `TheISurfaceTool::Value` / `Gradient` / `ValueAndGradient`
//!   (`IntSurf_QuadricTool.lxx:18-45`) are `IntSurf_Quadric::Distance` /
//!   `Gradient` / `ValAndGrad`; the port calls the methods directly, the way
//!   [`IntSurfQuadric`](crate::int_surf_quadric::IntSurfQuadric) exposes them.
//! - `math_FunctionSetRoot Rsnld(MyZerImpFunc)` (`gxx:552`) uses the default
//!   `NbIterations == 100` (`math_FunctionSetRoot.hxx:52`) and the
//!   four-argument `Perform` (`hxx:84`), i.e. `StopOnDivergent == false`.
//! - `MyParOnS2.SetCoord(tu2, tu2)` at `gxx:611` writes `tu2` for **both**
//!   coordinates. That is what OCCT does; the port keeps it and marks it.
//! - `Tg` / `Tguv1` / `Tguv2` are in/out and aliased by `aQuadTg` / `aPrmTg`
//!   (`gxx:448-449`) according to `MyImplicitFirst`; the port keeps the alias
//!   in a local and writes it back where OCCT writes the member.
//! - The two surfaces are read through `MyZerImpFunc` (`gxx:447-448`); the port
//!   clones them at the top of `Compute` because that method also writes
//!   `self`. The body never mutates either surface.

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::kernel::geomabs::SurfaceType;
use occt_core::precision::{ANGULAR, APPROXIMATION};
use occt_math::{MathFunctionSetRoot, MathVector};

use crate::approx_int::{SolverFlag, SvSurfaces};
use crate::brep_approx::imp_prm_tangent_processing::{
    is_singular, non_singular_processing, singular_processing,
};
use crate::brep_approx::{SurfaceTool, ZerImpFunc};
use crate::int_surf_quadric::IntSurfQuadric;
use crate::int_tools_wline::PntOn2S;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// The two 2D tangents ordered back into `(Tguv1, Tguv2)` according to
/// `MyImplicitFirst` (`gxx:448-449` aliases them the other way round).
fn split_tangents(quad_is_s1: bool, quad_tg: &GpVec2d, prm_tg: &GpVec2d) -> (GpVec2d, GpVec2d) {
    if quad_is_s1 {
        (*quad_tg, *prm_tg)
    } else {
        (*prm_tg, *quad_tg)
    }
}

/// `BRepApprox_TheImpPrmSvSurfacesOfApprox` (`hxx:39`).
#[derive(Clone)]
pub struct ImpPrmSvSurfaces {
    /// `MyParOnS1`.
    my_par_on_s1: GpPnt2d,
    /// `MyParOnS2`.
    my_par_on_s2: GpPnt2d,
    /// `MyPnt`.
    pub(crate) my_pnt: GpPnt,
    /// `MyTguv1`.
    my_tguv1: GpVec2d,
    /// `MyTguv2`.
    my_tguv2: GpVec2d,
    /// `MyTg`.
    my_tg: GpVec,
    /// `MyIsTangent`.
    my_is_tangent: bool,
    /// `MyHasBeenComputed`.
    pub(crate) my_has_been_computed: bool,
    /// `MyParOnS1bis`.
    my_par_on_s1bis: GpPnt2d,
    /// `MyParOnS2bis`.
    my_par_on_s2bis: GpPnt2d,
    /// `MyPntbis`.
    my_pntbis: GpPnt,
    /// `MyTguv1bis`.
    my_tguv1bis: GpVec2d,
    /// `MyTguv2bis`.
    my_tguv2bis: GpVec2d,
    /// `MyTgbis`.
    my_tgbis: GpVec,
    /// `MyIsTangentbis`.
    my_is_tangentbis: bool,
    /// `MyHasBeenComputedbis`.
    my_has_been_computedbis: bool,
    /// `MyImplicitFirst`.
    pub(crate) my_implicit_first: bool,
    /// `MyZerImpFunc`.
    pub(crate) my_zer_imp_func: ZerImpFunc,
    /// `ApproxInt_SvSurfaces::myUseSolver`.
    flag: SolverFlag,
}

impl ImpPrmSvSurfaces {
    /// `ApproxInt_ImpPrmSvSurfaces(const TheISurface& ISurf, const ThePSurface&
    /// PSurf)` (`gxx:325-338`): the quadric is surface 1.
    pub fn new_implicit_first(isurf: IntSurfQuadric, psurf: SurfaceAdaptor) -> Self {
        let mut this = Self::with_surfaces(psurf, isurf, true);
        this.set_use_solver(true);
        this
    }

    /// `ApproxInt_ImpPrmSvSurfaces(const ThePSurface& PSurf, const TheISurface&
    /// ISurf)` (`gxx:343-352`): the quadric is surface 2.
    pub fn new_param_first(psurf: SurfaceAdaptor, isurf: IntSurfQuadric) -> Self {
        let mut this = Self::with_surfaces(psurf, isurf, false);
        this.set_use_solver(true);
        this
    }

    /// The member initialisation list both constructors share
    /// (`gxx:327-331`, `:345-349`): everything false / zero except
    /// `MyImplicitFirst`, plus `MyZerImpFunc(PSurf, ISurf)`.
    fn with_surfaces(psurf: SurfaceAdaptor, isurf: IntSurfQuadric, implicit_first: bool) -> Self {
        Self {
            my_par_on_s1: GpPnt2d::new(0.0, 0.0),
            my_par_on_s2: GpPnt2d::new(0.0, 0.0),
            my_pnt: GpPnt::zero(),
            my_tguv1: GpVec2d::new(0.0, 0.0),
            my_tguv2: GpVec2d::new(0.0, 0.0),
            my_tg: GpVec::zero(),
            my_is_tangent: false,
            my_has_been_computed: false,
            my_par_on_s1bis: GpPnt2d::new(0.0, 0.0),
            my_par_on_s2bis: GpPnt2d::new(0.0, 0.0),
            my_pntbis: GpPnt::zero(),
            my_tguv1bis: GpVec2d::new(0.0, 0.0),
            my_tguv2bis: GpVec2d::new(0.0, 0.0),
            my_tgbis: GpVec::zero(),
            my_is_tangentbis: false,
            my_has_been_computedbis: false,
            my_implicit_first: implicit_first,
            my_zer_imp_func: ZerImpFunc::new(psurf, isurf),
            flag: SolverFlag::default(),
        }
    }

    /// `Compute(u1, v1, u2, v2, P, Tg, Tguv1, Tguv2)` (`gxx:437-768`).
    #[allow(clippy::too_many_arguments)]
    fn compute(
        &mut self,
        u1: &mut f64,
        v1: &mut f64,
        u2: &mut f64,
        v2: &mut f64,
        p: &mut GpPnt,
        tg: &mut GpVec,
        tguv1: &mut GpVec2d,
        tguv2: &mut GpVec2d,
    ) -> bool {
        let a_q_surf = self.my_zer_imp_func.isurface().clone();
        let a_p_surf = self.my_zer_imp_func.psurface().clone();
        let quad_is_s1 = self.my_implicit_first;

        // for square
        let a_null_value = APPROXIMATION * APPROXIMATION;
        let an_ang_tol = ANGULAR;

        let tu1 = *u1;
        let tu2 = *u2;
        let tv1 = *v1;
        let tv2 = *v2;

        if self.my_has_been_computed {
            if self.my_par_on_s1.x() == *u1
                && self.my_par_on_s1.y() == *v1
                && self.my_par_on_s2.x() == *u2
                && self.my_par_on_s2.y() == *v2
            {
                return self.my_is_tangent;
            } else if !self.my_has_been_computedbis {
                self.my_tgbis = self.my_tg;
                self.my_tguv1bis = self.my_tguv1;
                self.my_tguv2bis = self.my_tguv2;
                self.my_pntbis = self.my_pnt;
                self.my_par_on_s1bis = self.my_par_on_s1;
                self.my_par_on_s2bis = self.my_par_on_s2;
                self.my_is_tangentbis = self.my_is_tangent;
                self.my_has_been_computedbis = self.my_has_been_computed;
            }
        }

        if self.my_has_been_computedbis
            && self.my_par_on_s1bis.x() == *u1
            && self.my_par_on_s1bis.y() == *v1
            && self.my_par_on_s2bis.x() == *u2
            && self.my_par_on_s2bis.y() == *v2
        {
            let tv_vec = self.my_tg;
            let tv1 = self.my_tguv1;
            let tv2 = self.my_tguv2;
            let tp = self.my_pnt;
            let tp1 = self.my_par_on_s1;
            let tp2 = self.my_par_on_s2;
            let tb = self.my_is_tangent;

            self.my_tg = self.my_tgbis;
            self.my_tguv1 = self.my_tguv1bis;
            self.my_tguv2 = self.my_tguv2bis;
            self.my_pnt = self.my_pntbis;
            self.my_par_on_s1 = self.my_par_on_s1bis;
            self.my_par_on_s2 = self.my_par_on_s2bis;
            self.my_is_tangent = self.my_is_tangentbis;

            self.my_tgbis = tv_vec;
            self.my_tguv1bis = tv1;
            self.my_tguv2bis = tv2;
            self.my_pntbis = tp;
            self.my_par_on_s1bis = tp1;
            self.my_par_on_s2bis = tp2;
            self.my_is_tangentbis = tb;

            return self.my_is_tangent;
        }

        let mut x = MathVector::new(1, 2);
        let mut born_inf = MathVector::new(1, 2);
        let mut born_sup = MathVector::new(1, 2);
        let mut tolerance = MathVector::new(1, 2);
        // --- ThePSurfaceTool::GetResolution(aPSurf,Tolerance(1),Tolerance(2));
        tolerance.set_value(1, 1.0e-8);
        tolerance.set_value(2, 1.0e-8);
        let binfu = SurfaceTool::first_u_parameter(&a_p_surf);
        let binfv = SurfaceTool::first_v_parameter(&a_p_surf);
        let bsupu = SurfaceTool::last_u_parameter(&a_p_surf);
        let bsupv = SurfaceTool::last_v_parameter(&a_p_surf);
        born_inf.set_value(1, binfu);
        born_sup.set_value(1, bsupu);
        born_inf.set_value(2, binfv);
        born_sup.set_value(2, bsupv);
        let mut translation_u = 0.0_f64;
        let mut translation_v = 0.0_f64;

        if !self.fill_initial_vector_of_solution(
            *u1,
            *v1,
            *u2,
            *v2,
            binfu,
            bsupu,
            binfv,
            bsupv,
            &mut x,
            &mut translation_u,
            &mut translation_v,
        ) {
            self.my_is_tangent = false;
            self.my_is_tangentbis = false;
            self.my_has_been_computed = false;
            self.my_has_been_computedbis = false;
            return false;
        }

        let mut a_rsnld_is_done = false;
        let pour_tester_u = x.value(1);
        let pour_tester_v = x.value(2);
        if self.get_use_solver() {
            let mut rsnld = MathFunctionSetRoot::with_iterations(&self.my_zer_imp_func, 100);
            rsnld.set_tolerance(&tolerance);
            rsnld.perform_with_bounds(
                &mut self.my_zer_imp_func,
                &x,
                &born_inf,
                &born_sup,
                false,
            );
            a_rsnld_is_done = rsnld.is_done();
            if a_rsnld_is_done {
                x = rsnld.root();
            }
        }
        if !(a_rsnld_is_done || !self.get_use_solver()) {
            self.my_is_tangent = false;
            self.my_has_been_computed = false;
            self.my_has_been_computedbis = false;
            return false;
        }

        self.my_has_been_computed = true;

        let dist_avant_apres_u = (pour_tester_u - x.value(1)).abs();
        let dist_avant_apres_v = (pour_tester_v - x.value(2)).abs();

        self.my_pnt = SurfaceTool::value(&a_p_surf, x.value(1), x.value(2));
        *p = self.my_pnt;

        if dist_avant_apres_v > 0.001 || dist_avant_apres_u > 0.001 {
            // "ApproxInt_ImpImpSvSurfaces.gxx : Distance apres recadrage
            // Trop Grande"
            self.my_is_tangent = false;
            self.my_has_been_computed = false;
            self.my_has_been_computedbis = false;
            return false;
        }

        let mut d1u_prm = GpVec::zero();
        let mut d1v_prm = GpVec::zero();
        let (d1u_quad, d1v_quad) = if quad_is_s1 {
            *u2 = x.value(1) - translation_u;
            *v2 = x.value(2) - translation_v;

            if a_q_surf.type_quadric() != SurfaceType::Plane {
                while *u1 - tu1 > std::f64::consts::PI {
                    *u1 -= std::f64::consts::PI + std::f64::consts::PI;
                }
                while tu1 - *u1 > std::f64::consts::PI {
                    *u1 += std::f64::consts::PI + std::f64::consts::PI;
                }
            }

            self.my_par_on_s1.set_coord(tu1, tv1);
            self.my_par_on_s2.set_coord(tu2, tv2);

            SurfaceTool::d1(
                &a_p_surf,
                x.value(1),
                x.value(2),
                p,
                &mut d1u_prm,
                &mut d1v_prm,
            );
            let (a_p2, du, dv) = a_q_surf.d1(*u1, *v1);

            // Middle-point of P-P2 segment
            *p = GpPnt::from_xyz(&p.coord.added(&a_p2.coord).multiplied(0.5));
            (du, dv)
        } else {
            *u1 = x.value(1) - translation_u;
            *v1 = x.value(2) - translation_v;
            // aQSurf.Parameters(P, u2, v2);
            if a_q_surf.type_quadric() != SurfaceType::Plane {
                while *u2 - tu2 > std::f64::consts::PI {
                    *u2 -= std::f64::consts::PI + std::f64::consts::PI;
                }
                while tu2 - *u2 > std::f64::consts::PI {
                    *u2 += std::f64::consts::PI + std::f64::consts::PI;
                }
            }

            self.my_par_on_s1.set_coord(tu1, tv1);
            // OCCT writes `tu2` for both coordinates here (`gxx:611`).
            self.my_par_on_s2.set_coord(tu2, tu2);

            SurfaceTool::d1(
                &a_p_surf,
                x.value(1),
                x.value(2),
                p,
                &mut d1u_prm,
                &mut d1v_prm,
            );
            let (a_p2, du, dv) = a_q_surf.d1(*u2, *v2);

            // Middle-point of P-P2 segment
            *p = GpPnt::from_xyz(&p.coord.added(&a_p2.coord).multiplied(0.5));
            (du, dv)
        };

        self.my_pnt = *p;

        // Normals to the surfaces
        let mut a_normal_prm = d1u_prm.crossed(&d1v_prm);
        let mut a_normal_imp = a_q_surf.normale(&self.my_pnt);

        let a_sq_magn_prm = a_normal_prm.square_magnitude();
        let a_sq_magn_imp = a_normal_imp.square_magnitude();

        let mut is_prm_singular = false;
        let mut is_imp_singular = false;

        // `Tg` and the two aliases `aQuadTg` / `aPrmTg` (`gxx:448-449`).
        let mut a_tg = *tg;
        let mut a_quad_tg = if quad_is_s1 { *tguv1 } else { *tguv2 };
        let mut a_prm_tg = if quad_is_s1 { *tguv2 } else { *tguv1 };

        if is_singular(&d1u_prm, &d1v_prm, a_null_value, an_ang_tol) {
            is_prm_singular = true;
            if !singular_processing(
                &d1u_prm,
                &d1v_prm,
                true,
                a_null_value,
                an_ang_tol,
                &mut a_tg,
                &mut a_prm_tg,
            ) {
                self.my_is_tangent = false;
                self.my_has_been_computed = false;
                self.my_has_been_computedbis = false;
                return false;
            }
            self.my_tg = a_tg;
        } else {
            a_normal_prm = a_normal_prm.divided(a_sq_magn_prm.sqrt());
        }

        // Analogically for implicit surface
        if a_sq_magn_imp < a_null_value {
            is_imp_singular = true;

            if !singular_processing(
                &d1u_quad,
                &d1v_quad,
                !is_prm_singular,
                a_null_value,
                an_ang_tol,
                &mut a_tg,
                &mut a_quad_tg,
            ) {
                self.my_is_tangent = false;
                self.my_has_been_computed = false;
                self.my_has_been_computedbis = false;
                return false;
            }
            self.my_tg = a_tg;
        } else {
            a_normal_imp = a_normal_imp.divided(a_sq_magn_imp.sqrt());
        }

        if is_imp_singular && is_prm_singular {
            // All abnormal cases were processed above.
            let (uv1, uv2) = split_tangents(quad_is_s1, &a_quad_tg, &a_prm_tg);
            *tguv1 = uv1;
            *tguv2 = uv2;
            self.my_tguv1 = uv1;
            self.my_tguv2 = uv2;

            self.my_is_tangent = true;
            return self.my_is_tangent;
        } else if !(is_imp_singular || is_prm_singular) {
            // Processing pure non-singular case (the 3D and 2D tangents are
            // still not defined). aNormalImp and aNormalPrm are normalized
            // here, so |Tg| == 0 iff (aNormalImp || aNormalPrm).
            a_tg = a_normal_imp.crossed(&a_normal_prm);
        }

        let a_sq_magn_tg = a_tg.square_magnitude();

        if a_sq_magn_tg < a_null_value {
            self.my_is_tangent = false;
            self.my_has_been_computed = false;
            self.my_has_been_computedbis = false;
            return false;
        }

        // Normalize Tg vector
        a_tg = a_tg.divided(a_sq_magn_tg.sqrt());
        *tg = a_tg;
        self.my_tg = a_tg;

        if !is_prm_singular {
            // When isPrmSingular, aPrmTg has already been computed.
            if !non_singular_processing(
                &d1u_prm,
                &d1v_prm,
                &a_tg,
                a_null_value,
                an_ang_tol,
                &mut a_prm_tg,
            ) {
                self.my_is_tangent = false;
                self.my_has_been_computed = false;
                self.my_has_been_computedbis = false;
                return false;
            }
        }

        if !is_imp_singular {
            // When isImpSingular, aQuadTg has already been computed.
            if !non_singular_processing(
                &d1u_quad,
                &d1v_quad,
                &a_tg,
                a_null_value,
                an_ang_tol,
                &mut a_quad_tg,
            ) {
                self.my_is_tangent = false;
                self.my_has_been_computed = false;
                self.my_has_been_computedbis = false;
                return false;
            }
        }

        let (uv1, uv2) = split_tangents(quad_is_s1, &a_quad_tg, &a_prm_tg);
        *tguv1 = uv1;
        *tguv2 = uv2;
        self.my_tguv1 = uv1;
        self.my_tguv2 = uv2;

        self.my_is_tangent = true;

        true
    }
}

impl SvSurfaces for ImpPrmSvSurfaces {
    /// `Compute` (`gxx:437-768`).
    fn compute(
        &mut self,
        u1: &mut f64,
        v1: &mut f64,
        u2: &mut f64,
        v2: &mut f64,
        pt: &mut GpPnt,
        tg: &mut GpVec,
        tguv1: &mut GpVec2d,
        tguv2: &mut GpVec2d,
    ) -> bool {
        ImpPrmSvSurfaces::compute(self, u1, v1, u2, v2, pt, tg, tguv1, tguv2)
    }

    /// `Pnt` (`gxx:356-371`): `P = MyPnt` after the `Compute` call, whatever
    /// that call returned.
    fn pnt(&mut self, u1: f64, v1: f64, u2: f64, v2: f64, p: &mut GpPnt) {
        let mut a_p = GpPnt::zero();
        let mut a_t = GpVec::zero();
        let mut a_ts1 = GpVec2d::new(0.0, 0.0);
        let mut a_ts2 = GpVec2d::new(0.0, 0.0);
        let mut tu1 = u1;
        let mut tu2 = u2;
        let mut tv1 = v1;
        let mut tv2 = v2;
        let _ = ImpPrmSvSurfaces::compute(
            self,
            &mut tu1,
            &mut tv1,
            &mut tu2,
            &mut tv2,
            &mut a_p,
            &mut a_t,
            &mut a_ts1,
            &mut a_ts2,
        );
        *p = self.my_pnt;
    }

    /// `SeekPoint` (`gxx:776-861`).
    fn seek_point(&mut self, u1: f64, v1: f64, u2: f64, v2: f64, point: &mut PntOn2S) -> bool {
        ImpPrmSvSurfaces::seek_point(self, u1, v1, u2, v2, point)
    }

    /// `Tangency` (`gxx:374-390`).
    fn tangency(&mut self, u1: f64, v1: f64, u2: f64, v2: f64, tg: &mut GpVec) -> bool {
        let mut a_p = GpPnt::zero();
        let mut a_t = GpVec::zero();
        let mut a_ts1 = GpVec2d::new(0.0, 0.0);
        let mut a_ts2 = GpVec2d::new(0.0, 0.0);
        let mut tu1 = u1;
        let mut tu2 = u2;
        let mut tv1 = v1;
        let mut tv2 = v2;
        let t = ImpPrmSvSurfaces::compute(
            self,
            &mut tu1,
            &mut tv1,
            &mut tu2,
            &mut tv2,
            &mut a_p,
            &mut a_t,
            &mut a_ts1,
            &mut a_ts2,
        );
        *tg = self.my_tg;
        t
    }

    /// `TangencyOnSurf1` (`gxx:394-410`).
    fn tangency_on_surf1(&mut self, u1: f64, v1: f64, u2: f64, v2: f64, tg: &mut GpVec2d) -> bool {
        let mut a_p = GpPnt::zero();
        let mut a_t = GpVec::zero();
        let mut a_ts1 = GpVec2d::new(0.0, 0.0);
        let mut a_ts2 = GpVec2d::new(0.0, 0.0);
        let mut tu1 = u1;
        let mut tu2 = u2;
        let mut tv1 = v1;
        let mut tv2 = v2;
        let t = ImpPrmSvSurfaces::compute(
            self,
            &mut tu1,
            &mut tv1,
            &mut tu2,
            &mut tv2,
            &mut a_p,
            &mut a_t,
            &mut a_ts1,
            &mut a_ts2,
        );
        *tg = self.my_tguv1;
        t
    }

    /// `TangencyOnSurf2` (`gxx:414-430`).
    fn tangency_on_surf2(&mut self, u1: f64, v1: f64, u2: f64, v2: f64, tg: &mut GpVec2d) -> bool {
        let mut a_p = GpPnt::zero();
        let mut a_t = GpVec::zero();
        let mut a_ts1 = GpVec2d::new(0.0, 0.0);
        let mut a_ts2 = GpVec2d::new(0.0, 0.0);
        let mut tu1 = u1;
        let mut tu2 = u2;
        let mut tv1 = v1;
        let mut tv2 = v2;
        let t = ImpPrmSvSurfaces::compute(
            self,
            &mut tu1,
            &mut tv1,
            &mut tu2,
            &mut tv2,
            &mut a_p,
            &mut a_t,
            &mut a_ts1,
            &mut a_ts2,
        );
        *tg = self.my_tguv2;
        t
    }

    fn set_use_solver(&mut self, use_solver: bool) {
        self.flag.set_use_solver(use_solver);
    }

    fn get_use_solver(&self) -> bool {
        self.flag.get_use_solver()
    }
}
