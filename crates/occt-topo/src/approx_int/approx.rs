//! `ApproxInt_Approx` — the driver that fits a B-spline to the point set of a
//! `ApproxInt_MultiLine`, in Bezier pieces or in one B-spline.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/ApproxInt/ApproxInt_Approx.gxx:1-751`
//! (OCCT 8.0.0). The class members are declared in the instantiating headers
//! (`BRepApprox_Approx.hxx:45-91` carries the private `Approx_Data` struct and
//! the `myComputeLine` / `myComputeLineBezier` / `myBezToBSpl` / `myKnots`
//! members), so the private data layout below follows that header while the
//! bodies follow the `.gxx`.
//!
//! Translation notes:
//! - `myComputeLine` is `Approx_BSplComputeLine` over the multiline
//!   (`BRepApprox_TheComputeLineOfApprox_0.cxx:99` includes
//!   `Approx_BSplComputeLine.gxx`) and `myComputeLineBezier` is
//!   `Approx_ComputeLine` over it (`…TheComputeLineBezierOfApprox_0.cxx:69`
//!   includes `Approx_ComputeLine.gxx`); both are the generic drivers of
//!   [`occt_geom::approx`] parameterised by [`LineTool`].
//! - `TheWLine` / `ThePSurface` / `TheISurface` are template parameters of the
//!   `.gxx`; the port keeps the line generic over [`LineSource`] and leaves the
//!   surface-typed dispatch of `gxx:225-345` to the instantiation
//!   ([`crate::brep_approx::approx`]), which owns the `ThePSurfaceTool` calls.
//! - `void* const ptrsvsurf` (`gxx:333`, `:381`, `:539`, `:623`) is
//!   `NULL` on the plain `Perform` arms and the surface-pair evaluator on the
//!   others; the port models it as `Option<Rc<RefCell<dyn SvSurfaces>>>`, the
//!   same type [`MultiLine`] already carries.
//! - `NCollection_DynamicArray<int> myKnots` is 0-based (`Lower() == 0`,
//!   `Upper() == Length() - 1`), so `myKnots(kind)` is `my_knots[kind]`.
//! - `Perform()` (no arguments) is declared at `hxx` but never defined in the
//!   `.gxx`; it has no body to port and no caller in OCCT 8.0.0.
//! - `ComputeTrsf3d` / `ComputeTrsf2d` (`gxx:36-82`) are file-static helpers in
//!   OCCT; they stay file-private functions here.

use std::cell::RefCell;
use std::rc::Rc;

use occt_core::approx::ApproxParametrizationType;
use occt_core::gp::{GpPnt, GpPnt2d};
use occt_geom::app_par_curves::constraint::Constraint;
use occt_geom::app_par_curves::multi_bspline_curve::MultiBSpCurve;
use occt_geom::approx::bspl_compute_line::BsplComputeLine;
use occt_geom::approx::compute_line::ComputeLine;
use occt_geom::approx::line_tool::LineTool;
use occt_geom::approx::m_curves_to_bsp_curve::MCurvesToBSpCurve;
use occt_math::MathVector;

use crate::approx_int::approx_data::{compute_trsf2d, compute_trsf3d, line_values, ApproxData};
use crate::approx_int::knot_tools;
use crate::approx_int::multi_line::{LineSource, MultiLine};
use crate::approx_int::SvSurfaces;

/// `aMinNbPointsForApprox` (`gxx:28`): below this point count the driver
/// interpolates instead of cutting into Bezier pieces.
const A_MIN_NB_POINTS_FOR_APPROX: i32 = 5;

/// `RatioTol` (`gxx:31`).
const RATIO_TOL: f64 = 1.5;

/// `ApproxInt_Approx` (`gxx:164`).
#[derive(Debug, Clone)]
pub struct ApproxIntApprox {
    /// `myComputeLine` (`hxx`): `BRepApprox_TheComputeLineOfApprox`.
    my_compute_line: BsplComputeLine,
    /// `myComputeLineBezier`: `BRepApprox_TheComputeLineBezierOfApprox`.
    my_compute_line_bezier: ComputeLine,
    /// `myBezToBSpl`.
    my_bez_to_bspl: MCurvesToBSpCurve,
    /// `myWithTangency`.
    my_with_tangency: bool,
    /// `myTol3d`.
    my_tol3d: f64,
    /// `myTol2d`.
    my_tol2d: f64,
    /// `myDegMin`.
    my_deg_min: i32,
    /// `myDegMax`.
    my_deg_max: i32,
    /// `myNbIterMax`.
    my_nb_iter_max: i32,
    /// `myTolReached3d`.
    my_tol_reached3d: f64,
    /// `myTolReached2d`.
    my_tol_reached2d: f64,
    /// `myData`.
    my_data: ApproxData,
    /// `myKnots`.
    my_knots: Vec<i32>,
}

impl Default for ApproxIntApprox {
    fn default() -> Self {
        Self::new()
    }
}

impl ApproxIntApprox {
    /// `ApproxInt_Approx()` (`gxx:164-181`).
    pub fn new() -> Self {
        let mut my_compute_line = BsplComputeLine::with_settings(
            4,
            8,
            0.001,
            0.001,
            5,
            true,
            ApproxParametrizationType::ChordLength,
            false,
        );
        let my_compute_line_bezier = ComputeLine::with_settings(
            4,
            8,
            0.001,
            0.001,
            5,
            true,
            ApproxParametrizationType::ChordLength,
            false,
        );
        my_compute_line.set_continuity(2);
        // `// myComputeLineBezier.SetContinuity(2);` (`gxx:179`).
        Self {
            my_compute_line,
            my_compute_line_bezier,
            my_bez_to_bspl: MCurvesToBSpCurve::new(),
            my_with_tangency: true,
            my_tol3d: 0.001,
            my_tol2d: 0.001,
            my_deg_min: 4,
            my_deg_max: 8,
            my_nb_iter_max: 5,
            my_tol_reached3d: 0.0,
            my_tol_reached2d: 0.0,
            my_data: ApproxData::default(),
            my_knots: Vec::new(),
        }
    }

    /// `Parameters(Line, firstP, lastP, Par, TheParameters)` (`gxx:86-162`).
    pub fn parameters<L: LineTool>(
        line: &L,
        first_p: i32,
        last_p: i32,
        par: ApproxParametrizationType,
        the_parameters: &mut MathVector,
    ) {
        let nb_p3d = line.nb_p3d();
        let nb_p2d = line.nb_p2d();

        if par == ApproxParametrizationType::ChordLength
            || par == ApproxParametrizationType::Centripetal
        {
            the_parameters.set_value(first_p as usize, 0.0);
            for i in (first_p + 1)..=last_p {
                let (tab_p, tab_p2d) = line_values(line, i - 1, nb_p3d, nb_p2d);
                let (tab_pp, tab_pp2d) = line_values(line, i, nb_p3d, nb_p2d);

                let mut dist = 0.0;
                for j in 0..nb_p3d as usize {
                    dist += tab_pp[j].square_distance(&tab_p[j]);
                }
                for j in 0..nb_p2d as usize {
                    dist += tab_pp2d[j].square_distance(&tab_p2d[j]);
                }

                let dist = dist.sqrt();
                let prev = the_parameters.value((i - 1) as usize);
                let value = if par == ApproxParametrizationType::ChordLength {
                    prev + dist
                } else {
                    // `Par == Approx_Centripetal` (`gxx:147-150`).
                    prev + dist.sqrt()
                };
                the_parameters.set_value(i as usize, value);
            }
            let last = the_parameters.value(last_p as usize);
            for i in first_p..=last_p {
                let value = the_parameters.value(i as usize) / last;
                the_parameters.set_value(i as usize, value);
            }
        } else {
            for i in first_p..=last_p {
                let value =
                    ((i - first_p) as f64) / ((last_p - first_p) as f64);
                the_parameters.set_value(i as usize, value);
            }
        }
    }

    /// `Perform(theline, ApproxXYZ, ApproxU1V1, ApproxU2V2, indicemin,
    /// indicemax)` (`gxx:183-218`): the arm without a surface pair.
    #[allow(clippy::too_many_arguments)]
    pub fn perform_line<L: LineSource>(
        &mut self,
        theline: &L,
        approx_xyz: bool,
        approx_u1v1: bool,
        approx_u2v2: bool,
        indicemin: i32,
        indicemax: i32,
    ) {
        // Prepare DS (`gxx:191-192`).
        self.prepare_ds(approx_xyz, approx_u1v1, approx_u2v2, indicemin, indicemax);

        let nbpntbez = self.my_data.indicemax - self.my_data.indicemin;
        self.my_data.my_bezier_approx = nbpntbez >= A_MIN_NB_POINTS_FOR_APPROX;

        // Fill data structure (`gxx:200-201`).
        self.fill_data(theline);

        // Build knots (`gxx:203-210`).
        self.build_knots(theline, None);
        if self.my_knots.len() == 2
            && indicemax - indicemin > 2 * self.my_data.my_nb_pnt_max
        {
            // At least 3 knots for BrepApprox.
            let last = self.my_knots.len() - 1;
            self.my_knots[last] = (indicemax - indicemin) / 2;
            self.my_knots.push(indicemax);
        }

        self.my_compute_line.init(
            self.my_deg_min,
            self.my_deg_max,
            self.my_tol3d,
            self.my_tol2d,
            self.my_nb_iter_max,
            true,
            self.my_data.parametrization,
            false,
        );
        self.my_compute_line_bezier.init(
            self.my_deg_min,
            self.my_deg_max,
            self.my_tol3d,
            self.my_tol2d,
            self.my_nb_iter_max,
            true,
            self.my_data.parametrization,
            false,
        );

        self.build_curve(theline, None);
    }

    /// The body shared by `Perform(Surf1, Surf2, …)` (`gxx:311-344`) and
    /// `Perform(ISurf, PSurf, …)` (`gxx:360-395`): the surface-pair arms of the
    /// `.gxx`, with the pair evaluator already built by the caller.
    ///
    /// Both OCCT arms differ only in the evaluator they construct (and in
    /// `SetUseSolver(false)` for the analytic-parametric one), so they share
    /// this body verbatim.
    #[allow(clippy::too_many_arguments)]
    fn prepare_with_evaluator<L: LineSource>(
        &mut self,
        sv_surfaces: Rc<RefCell<dyn SvSurfaces>>,
        theline: &L,
        approx_xyz: bool,
        approx_u1v1: bool,
        approx_u2v2: bool,
        indicemin: i32,
        indicemax: i32,
    ) {
        // Prepare DS (`gxx:312-313`, `:361-362`).
        self.prepare_ds(approx_xyz, approx_u1v1, approx_u2v2, indicemin, indicemax);

        let nbpntbez = indicemax - indicemin;
        self.my_data.my_bezier_approx = nbpntbez >= A_MIN_NB_POINTS_FOR_APPROX;

        let cut = self.my_data.my_bezier_approx;

        // Fill data structure (`gxx:328-329`, `:385-386`).
        self.fill_data(theline);

        // Build knots (`gxx:335`, `:388`).
        self.build_knots(theline, Some(sv_surfaces.clone()));

        self.my_compute_line.init(
            self.my_deg_min,
            self.my_deg_max,
            self.my_tol3d,
            self.my_tol2d,
            self.my_nb_iter_max,
            cut,
            self.my_data.parametrization,
            false,
        );
        self.my_compute_line_bezier.init(
            self.my_deg_min,
            self.my_deg_max,
            self.my_tol3d,
            self.my_tol2d,
            self.my_nb_iter_max,
            cut,
            self.my_data.parametrization,
            false,
        );

        self.build_curve(theline, Some(sv_surfaces));
    }

    /// `Perform(Surf1, Surf2, …)`, non-analytical arm (`gxx:311-344`): the
    /// `ApproxInt_ThePrmPrmSvSurfaces` evaluator (`gxx:316`).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_with_surfaces<L: LineSource>(
        &mut self,
        sv_surfaces: Rc<RefCell<dyn SvSurfaces>>,
        theline: &L,
        approx_xyz: bool,
        approx_u1v1: bool,
        approx_u2v2: bool,
        indicemin: i32,
        indicemax: i32,
    ) {
        self.prepare_with_evaluator(
            sv_surfaces,
            theline,
            approx_xyz,
            approx_u1v1,
            approx_u2v2,
            indicemin,
            indicemax,
        );
    }

    /// `Perform(ISurf, PSurf, …)`, analytic-parametric arm (`gxx:355-395`):
    /// `myImpPrmSvSurfaces.SetUseSolver(false)` (`gxx:369`) before the shared
    /// body, so `FillInitialVectorOfSolution` takes the parameters as exact.
    #[allow(clippy::too_many_arguments)]
    pub fn perform_imp_prm<L: LineSource>(
        &mut self,
        sv_surfaces: Rc<RefCell<dyn SvSurfaces>>,
        theline: &L,
        approx_xyz: bool,
        approx_u1v1: bool,
        approx_u2v2: bool,
        indicemin: i32,
        indicemax: i32,
    ) {
        sv_surfaces.borrow_mut().set_use_solver(false);
        self.prepare_with_evaluator(
            sv_surfaces,
            theline,
            approx_xyz,
            approx_u1v1,
            approx_u2v2,
            indicemin,
            indicemax,
        );
    }

    /// `SetParameters(Tol3d, Tol2d, DegMin, DegMax, NbIterMax, NbPntMax,
    /// ApproxWithTangency, Parametrization)` (`gxx:398-428`).
    #[allow(clippy::too_many_arguments)]
    pub fn set_parameters(
        &mut self,
        tol3d: f64,
        tol2d: f64,
        deg_min: i32,
        deg_max: i32,
        nb_iter_max: i32,
        nb_pnt_max: i32,
        approx_with_tangency: bool,
        parametrization: ApproxParametrizationType,
    ) {
        self.my_data.my_nb_pnt_max = nb_pnt_max;
        self.my_with_tangency = approx_with_tangency;
        self.my_tol3d = tol3d / RATIO_TOL;
        self.my_tol2d = tol2d / RATIO_TOL;
        self.my_deg_min = deg_min;
        self.my_deg_max = deg_max;
        self.my_nb_iter_max = nb_iter_max;

        self.my_compute_line.init(
            self.my_deg_min,
            self.my_deg_max,
            self.my_tol3d,
            self.my_tol2d,
            self.my_nb_iter_max,
            true,
            parametrization,
            false,
        );
        self.my_compute_line_bezier.init(
            self.my_deg_min,
            self.my_deg_max,
            self.my_tol3d,
            self.my_tol2d,
            self.my_nb_iter_max,
            true,
            parametrization,
            false,
        );

        if !approx_with_tangency {
            self.my_compute_line
                .set_constraints(Constraint::PassPoint, Constraint::PassPoint);
            self.my_compute_line_bezier
                .set_constraints(Constraint::PassPoint, Constraint::PassPoint);
        }

        self.my_data.my_bezier_approx = true;
    }

    /// `NbMultiCurves()` (`gxx:432-435`).
    pub fn nb_multi_curves(&self) -> i32 {
        1
    }

    /// `UpdateTolReached()` (`gxx:438-454`).
    pub fn update_tol_reached(&mut self) {
        if self.my_data.my_bezier_approx {
            let nb_curves = self.my_compute_line_bezier.nb_multi_curves();
            for i_cur in 1..=nb_curves {
                let (tol3d, tol2d) = self.my_compute_line_bezier.error(i_cur as i32);
                self.my_tol_reached3d = self.my_tol_reached3d.max(tol3d);
                self.my_tol_reached2d = self.my_tol_reached2d.max(tol2d);
            }
        } else {
            let (tol3d, tol2d) = self.my_compute_line.error();
            self.my_tol_reached3d = tol3d;
            self.my_tol_reached2d = tol2d;
        }
    }

    /// `TolReached3d()` (`gxx:458-462`).
    pub fn tol_reached_3d(&self) -> f64 {
        self.my_tol_reached3d * RATIO_TOL
    }

    /// `TolReached2d()` (`gxx:464-468`).
    pub fn tol_reached_2d(&self) -> f64 {
        self.my_tol_reached2d * RATIO_TOL
    }

    /// `IsDone()` (`gxx:472-483`).
    pub fn is_done(&self) -> bool {
        if self.my_data.my_bezier_approx {
            self.my_compute_line_bezier.nb_multi_curves() > 0
        } else {
            self.my_compute_line.is_tolerance_reached()
        }
    }

    /// `Value(Index)` (`gxx:486-498`): `Index` is unused.
    pub fn value(&self) -> &MultiBSpCurve {
        if self.my_data.my_bezier_approx {
            self.my_bez_to_bspl.value()
        } else {
            self.my_compute_line.value()
        }
    }

    /// Zero both reached tolerances (`gxx:238`).
    pub fn reset_tol_reached(&mut self) {
        self.my_tol_reached3d = 0.0;
        self.my_tol_reached2d = 0.0;
    }

    /// `fillData(theline)` (`gxx:502-517`).
    fn fill_data<L: LineSource>(&mut self, theline: &L) {
        if self.my_data.approx_xyz {
            let (xo, yo, zo) = compute_trsf3d(theline);
            self.my_data.xo = xo;
            self.my_data.yo = yo;
            self.my_data.zo = zo;
        } else {
            self.my_data.xo = 0.0;
            self.my_data.yo = 0.0;
            self.my_data.zo = 0.0;
        }

        if self.my_data.approx_u1v1 {
            let (u1o, v1o) = compute_trsf2d(theline, true);
            self.my_data.u1o = u1o;
            self.my_data.v1o = v1o;
        } else {
            self.my_data.u1o = 0.0;
            self.my_data.v1o = 0.0;
        }

        if self.my_data.approx_u2v2 {
            let (u2o, v2o) = compute_trsf2d(theline, false);
            self.my_data.u2o = u2o;
            self.my_data.v2o = v2o;
        } else {
            self.my_data.u2o = 0.0;
            self.my_data.v2o = 0.0;
        }
    }

    /// `prepareDS(ApproxXYZ, ApproxU1V1, ApproxU2V2, indicemin, indicemax)`
    /// (`gxx:521-534`).
    fn prepare_ds(
        &mut self,
        the_approx_xyz: bool,
        the_approx_u1v1: bool,
        the_approx_u2v2: bool,
        the_indicemin: i32,
        the_indicemax: i32,
    ) {
        self.my_tol_reached3d = 0.0;
        self.my_tol_reached2d = 0.0;
        self.my_data.approx_u1v1 = the_approx_u1v1;
        self.my_data.approx_u2v2 = the_approx_u2v2;
        self.my_data.approx_xyz = the_approx_xyz;
        self.my_data.indicemin = the_indicemin;
        self.my_data.indicemax = the_indicemax;
        self.my_data.parametrization = self.my_compute_line_bezier.parametrization();
    }

    /// The `ApproxInt_TheMultiLine` construction shared by `buildKnots`
    /// (`gxx:548-564`), `buildCurve` (`gxx:638-653`) and the insertion helpers
    /// of `ApproxInt_MultiLine`: `thePtrSVSurf == NULL` is the constructor
    /// without a surface pair (`ApproxInt_MultiLine.gxx:98-137`).
    fn make_multi_line<L: LineSource>(
        &self,
        theline: &L,
        sv_surfaces: Option<Rc<RefCell<dyn SvSurfaces>>>,
        ind_min: i32,
        ind_max: i32,
    ) -> MultiLine<L> {
        let d = &self.my_data;
        let nb_p3d = if d.approx_xyz { 1 } else { 0 };
        let nb_p2d = (if d.approx_u1v1 { 1 } else { 0 }) + (if d.approx_u2v2 { 1 } else { 0 });
        let offsets = [d.xo, d.yo, d.zo, d.u1o, d.v1o, d.u2o, d.v2o];
        match sv_surfaces {
            Some(sv) => MultiLine::with_surfaces(
                theline.clone(),
                sv,
                nb_p3d,
                nb_p2d,
                d.approx_u1v1,
                d.approx_u2v2,
                offsets,
                d.approx_u1v1,
                ind_min,
                ind_max,
            ),
            None => MultiLine::without_surfaces(
                theline.clone(),
                nb_p3d,
                nb_p2d,
                d.approx_u1v1,
                d.approx_u2v2,
                offsets,
                d.approx_u1v1,
                ind_min,
                ind_max,
            ),
        }
    }

    /// `buildKnots(theline, thePtrSVSurf)` (`gxx:537-619`).
    fn build_knots<L: LineSource>(
        &mut self,
        theline: &L,
        the_ptr_sv_surf: Option<Rc<RefCell<dyn SvSurfaces>>>,
    ) {
        self.my_knots.clear();
        if !self.my_data.my_bezier_approx {
            self.my_knots.push(self.my_data.indicemin);
            self.my_knots.push(self.my_data.indicemax);
            return;
        }

        let lo = self.my_data.indicemin;
        let hi = self.my_data.indicemax;
        let a_test_line = self.make_multi_line(theline, the_ptr_sv_surf, lo, hi);

        let nbp3d = a_test_line.nb_p3d();
        let nbp2d = a_test_line.nb_p2d();
        let nb = (hi - lo + 1).max(0) as usize;
        let mut a_pnt_xyz = vec![GpPnt::zero(); nb];
        let mut a_pnt_u1v1 = vec![GpPnt2d::new(0.0, 0.0); nb];
        let mut a_pnt_u2v2 = vec![GpPnt2d::new(0.0, 0.0); nb];

        for i in lo..=hi {
            let idx = (i - lo) as usize;
            let (tab_pnt3d, tab_pnt2d) = line_values(&a_test_line, i, nbp3d, nbp2d);
            if nbp3d > 0 {
                a_pnt_xyz[idx] = tab_pnt3d[0];
            }
            if nbp2d > 1 {
                a_pnt_u1v1[idx] = tab_pnt2d[0];
                a_pnt_u2v2[idx] = tab_pnt2d[1];
            } else if nbp2d > 0 {
                if self.my_data.approx_u1v1 {
                    a_pnt_u1v1[idx] = tab_pnt2d[0];
                } else {
                    a_pnt_u2v2[idx] = tab_pnt2d[0];
                }
            }
        }

        let a_min_nb_pnts = self.my_data.my_nb_pnt_max;

        // Expected parametrization (`gxx:607-608`).
        let mut a_pars = MathVector::new(lo as usize, hi as usize);
        Self::parameters(
            &a_test_line,
            lo,
            hi,
            self.my_data.parametrization,
            &mut a_pars,
        );

        knot_tools::build_knots(
            &a_pnt_xyz,
            &a_pnt_u1v1,
            &a_pnt_u2v2,
            &a_pars,
            self.my_data.approx_xyz,
            self.my_data.approx_u1v1,
            self.my_data.approx_u2v2,
            a_min_nb_pnts,
            &mut self.my_knots,
        );
    }

    /// `buildCurve(theline, thePtrSVSurf)` (`gxx:622-751`).
    fn build_curve<L: LineSource>(
        &mut self,
        theline: &L,
        the_ptr_sv_surf: Option<Rc<RefCell<dyn SvSurfaces>>>,
    ) {
        if self.my_data.my_bezier_approx {
            self.my_bez_to_bspl.reset();
        }

        let mut kind = 0usize;
        let mut other_inter;
        loop {
            // Base cycle: iterate over knots (`gxx:633-635`).
            let imin = self.my_knots[kind];
            let imax = self.my_knots[kind + 1];
            let my_multi_line =
                self.make_multi_line(theline, the_ptr_sv_surf.clone(), imin, imax);

            if self.my_data.my_bezier_approx {
                self.my_compute_line_bezier.perform(&my_multi_line);
                if self.my_compute_line_bezier.nb_multi_curves() == 0 {
                    return;
                }
            } else {
                self.my_compute_line.perform(&my_multi_line);
            }

            self.update_tol_reached();

            let mut indice3d = 1;
            let mut indice2d1 = 2;
            let mut indice2d2 = 3;
            if !self.my_data.approx_xyz {
                indice2d1 -= 1;
                indice2d2 -= 1;
            }
            if !self.my_data.approx_u1v1 {
                indice2d2 -= 1;
            }

            let (xo, yo, zo, u1o, v1o, u2o, v2o) = (
                self.my_data.xo,
                self.my_data.yo,
                self.my_data.zo,
                self.my_data.u1o,
                self.my_data.v1o,
                self.my_data.u2o,
                self.my_data.v2o,
            );

            if self.my_data.approx_xyz {
                if self.my_data.my_bezier_approx {
                    let nb = self.my_compute_line_bezier.nb_multi_curves();
                    for nbmc in (1..=nb).rev() {
                        self.my_compute_line_bezier
                            .change_value(nbmc as i32)
                            .transform(indice3d, -xo, 1.0, -yo, 1.0, -zo, 1.0);
                    }
                } else {
                    self.my_compute_line
                        .change_value()
                        .base_mut()
                        .transform(indice3d, -xo, 1.0, -yo, 1.0, -zo, 1.0);
                }
            }
            if self.my_data.approx_u1v1 {
                if self.my_data.my_bezier_approx {
                    let nb = self.my_compute_line_bezier.nb_multi_curves();
                    for nbmc in (1..=nb).rev() {
                        self.my_compute_line_bezier
                            .change_value(nbmc as i32)
                            .transform2d(indice2d1, -u1o, 1.0, -v1o, 1.0);
                    }
                } else {
                    self.my_compute_line
                        .change_value()
                        .base_mut()
                        .transform2d(indice2d1, -u1o, 1.0, -v1o, 1.0);
                }
            }
            if self.my_data.approx_u2v2 {
                if self.my_data.my_bezier_approx {
                    let nb = self.my_compute_line_bezier.nb_multi_curves();
                    for nbmc in (1..=nb).rev() {
                        self.my_compute_line_bezier
                            .change_value(nbmc as i32)
                            .transform2d(indice2d2, -u2o, 1.0, -v2o, 1.0);
                    }
                } else {
                    self.my_compute_line
                        .change_value()
                        .base_mut()
                        .transform2d(indice2d2, -u2o, 1.0, -v2o, 1.0);
                }
            }

            other_inter = false;
            if self.my_data.my_bezier_approx {
                let nb = self.my_compute_line_bezier.nb_multi_curves();
                for nbmc in 1..=nb {
                    let mc = self.my_compute_line_bezier.value(nbmc as i32).clone();
                    self.my_bez_to_bspl.append(&mc);
                }
                kind += 1;
                if kind < self.my_knots.len() - 1 {
                    other_inter = true;
                }
            }

            if !other_inter {
                break;
            }
        }

        if self.my_data.my_bezier_approx {
            self.my_bez_to_bspl.perform();
        }
    }
}
