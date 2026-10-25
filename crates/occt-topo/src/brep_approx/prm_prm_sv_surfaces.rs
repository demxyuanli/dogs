//! `BRepApprox_ThePrmPrmSvSurfacesOfApprox` — the two-parametric surface
//! evaluator (`ApproxInt_PrmPrmSvSurfaces.gxx` instantiated for
//! `BRepAdaptor_Surface` / `BRepApprox_SurfaceTool`).
//!
//! Sources:
//! - `ModelingAlgorithms/TKGeomAlgo/ApproxInt/ApproxInt_PrmPrmSvSurfaces.gxx:28-373`
//! - `ModelingAlgorithms/TKTopAlgo/BRepApprox/BRepApprox_ThePrmPrmSvSurfacesOfApprox_0.cxx:19-47`
//!   (the `#define ThePSurface BRepAdaptor_Surface` / `ThePSurfaceTool
//!   BRepApprox_SurfaceTool` / `ApproxInt_TheInt2S
//!   BRepApprox_TheInt2SOfThePrmPrmSvSurfacesOfApprox` binding)
//! - `BRepApprox_ThePrmPrmSvSurfacesOfApprox.hxx:34-108`
//!
//! `MyIntersectionOn2S` is `TheInt2S` over the two surfaces,
//! [`crate::intwalk_int2s::TheInt2S`], which `IntImp_Int2S.gxx` provides with
//! the same `Perform(Param, Rsnld)` / `IsDone` / `IsEmpty` / `IsTangent` /
//! `Point` / `Direction` / `DirectionOnS1` / `DirectionOnS2` surface. OCCT holds
//! it as a member; the port rebuilds it per `Compute` because
//! `TheInt2S<'a>` borrows the surfaces while this struct must own them. The
//! object carries no state across `Perform` calls that `Compute` does not
//! re-read afterwards, so the rebuild is observationally the same.
//!
//! `TOLTANGENCY` (`gxx:18`) is `1e-10`, passed as the squared tolerance by
//! `TheInt2S::new` exactly as `AppParCurves_Constraint`-independent code does.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::precision::ANGULAR;
use occt_geom::surface::Surface;

use crate::approx_int::{SvSurfaces, SolverFlag};
use crate::brep_approx::SurfaceTool;
use crate::int_tools_wline::PntOn2S;
use crate::intwalk::TheInt2S;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `#define TOLTANGENCY 0.0000000001` (`gxx:18`).
const TOL_TANGENCY: f64 = 1.0e-10;

/// `BRepApprox_ThePrmPrmSvSurfacesOfApprox` (`hxx:40`).
pub struct PrmPrmSvSurfaces {
    /// `MyIntersectionOn2S.Function().AuxillarSurface1()`, the `Surf1` member.
    surf1: SurfaceAdaptor,
    /// `AuxillarSurface2()`.
    surf2: SurfaceAdaptor,
    /// The evaluators' world surfaces, kept so [`TheInt2S`] can borrow them.
    world1: Option<Arc<dyn Surface>>,
    world2: Option<Arc<dyn Surface>>,
    /// `MyIsTangent` (`gxx:32`).
    my_is_tangent: bool,
    /// `MyHasBeenComputed`.
    my_has_been_computed: bool,
    /// `MyIsTangentbis`.
    my_is_tangentbis: bool,
    /// `MyHasBeenComputedbis`.
    my_has_been_computedbis: bool,
    /// `MyTg`.
    my_tg: GpVec,
    /// `MyTguv1`.
    my_tguv1: GpVec2d,
    /// `MyTguv2`.
    my_tguv2: GpVec2d,
    /// `MyPnt`.
    my_pnt: GpPnt,
    /// `MyParOnS1` — the parameters the previous `Compute` was asked for.
    my_par_on_s1: GpPnt2d,
    /// `MyParOnS2`.
    my_par_on_s2: GpPnt2d,
    /// `MyTgbis` / `MyTguv1bis` / `MyTguv2bis` / `MyPntbis` /
    /// `MyParOnS1bis` / `MyParOnS2bis` (`gxx:66-74`).
    my_tgbis: GpVec,
    my_tguv1bis: GpVec2d,
    my_tguv2bis: GpVec2d,
    my_pntbis: GpPnt,
    my_par_on_s1bis: GpPnt2d,
    my_par_on_s2bis: GpPnt2d,
    /// `ApproxInt_SvSurfaces::myUseSolver`.
    flag: SolverFlag,
}

impl PrmPrmSvSurfaces {
    /// `ApproxInt_PrmPrmSvSurfaces(Surf1, Surf2)`
    /// (`BRepApprox_ThePrmPrmSvSurfacesOfApprox.hxx:38-39`, `gxx:28-37`).
    pub fn new(surf1: SurfaceAdaptor, surf2: SurfaceAdaptor) -> Self {
        let world1 = surf1.surface().cloned();
        let world2 = surf2.surface().cloned();
        Self {
            surf1,
            surf2,
            world1,
            world2,
            my_is_tangent: false,
            my_has_been_computed: false,
            my_is_tangentbis: false,
            my_has_been_computedbis: false,
            my_tg: GpVec::new(0.0, 0.0, 0.0),
            my_tguv1: GpVec2d::new(0.0, 0.0),
            my_tguv2: GpVec2d::new(0.0, 0.0),
            my_pnt: GpPnt::zero(),
            my_par_on_s1: GpPnt2d::new(0.0, 0.0),
            my_par_on_s2: GpPnt2d::new(0.0, 0.0),
            my_tgbis: GpVec::new(0.0, 0.0, 0.0),
            my_tguv1bis: GpVec2d::new(0.0, 0.0),
            my_tguv2bis: GpVec2d::new(0.0, 0.0),
            my_pntbis: GpPnt::zero(),
            my_par_on_s1bis: GpPnt2d::new(0.0, 0.0),
            my_par_on_s2bis: GpPnt2d::new(0.0, 0.0),
            flag: SolverFlag::default(),
        }
    }
}

impl SvSurfaces for PrmPrmSvSurfaces {
    /// `Compute` (`gxx:44-226`).
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
        let tu1 = *u1;
        let tu2 = *u2;
        let tv1 = *v1;
        let tv2 = *v2;

        // `gxx:59-78`: the previous call's answer, verbatim.
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

        // `gxx:79-111`: the one-before-last answer, swapped into the current
        // slots so the next call can hit the cache again.
        if self.my_has_been_computedbis
            && self.my_par_on_s1bis.x() == *u1
            && self.my_par_on_s1bis.y() == *v1
            && self.my_par_on_s2bis.x() == *u2
            && self.my_par_on_s2bis.y() == *v2
        {
            let tv = self.my_tg;
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

            self.my_tgbis = tv;
            self.my_tguv1bis = tv1;
            self.my_tguv2bis = tv2;
            self.my_pntbis = tp;
            self.my_par_on_s1bis = tp1;
            self.my_par_on_s2bis = tp2;
            self.my_is_tangentbis = tb;

            return self.my_is_tangent;
        }

        self.my_is_tangent = true;

        // `gxx:114-121`: `MyIntersectionOn2S.Perform(Param, Rsnld)`.
        let (Some(w1), Some(w2)) = (self.world1.as_ref(), self.world2.as_ref()) else {
            self.my_has_been_computed = false;
            self.my_has_been_computedbis = false;
            return false;
        };
        let mut inter = TheInt2S::new(w1.as_ref(), w2.as_ref(), TOL_TANGENCY);
        let _ = inter.perform([*u1, *v1, *u2, *v2]);
        if !inter.is_done() {
            self.my_has_been_computed = false;
            self.my_has_been_computedbis = false;
            return false;
        }
        if inter.is_empty() {
            self.my_is_tangent = false;
            self.my_has_been_computed = false;
            self.my_has_been_computedbis = false;
            return false;
        }
        self.my_has_been_computed = true;
        self.my_pnt = inter.point().value();
        *p = self.my_pnt;

        let (iu1, iv1, iu2, iv2) = inter.point().parameters();
        *u1 = iu1;
        *v1 = iv1;
        *u2 = iu2;
        *v2 = iv2;
        self.my_par_on_s1.set_coord(tu1, tv1);
        self.my_par_on_s2.set_coord(tu2, tv2);

        if inter.is_tangent() {
            self.my_is_tangent = false;
            self.my_has_been_computed = false;
            self.my_has_been_computedbis = false;
            return false;
        }
        let d3d = inter.direction();
        let d2d1 = inter.direction_on_s1();
        let d2d2 = inter.direction_on_s2();
        self.my_tg = GpVec::new(d3d.x(), d3d.y(), d3d.z());
        *tg = self.my_tg;
        self.my_tguv1 = GpVec2d::new(d2d1.x, d2d1.y);
        *tguv1 = self.my_tguv1;
        self.my_tguv2 = GpVec2d::new(d2d2.x, d2d2.y);
        *tguv2 = self.my_tguv2;

        // `gxx:172-197`
        tg.normalize();
        self.my_tg = *tg;

        let mut pbid = GpPnt::zero();
        let mut tu = GpVec::new(0.0, 0.0, 0.0);
        let mut tv = GpVec::new(0.0, 0.0, 0.0);
        SurfaceTool::d1(&self.surf1, *u1, *v1, &mut pbid, &mut tu, &mut tv);
        let tutu = tu.dot(&tu);
        let tvtv = tv.dot(&tv);
        let tutv = tu.dot(&tv);
        let tgtu = tg.dot(&tu);
        let tgtv = tg.dot(&tv);
        let dis = tutu * tvtv - tutv * tutv;
        if dis.abs() < ANGULAR {
            self.my_is_tangent = false;
            self.my_has_been_computed = false;
            self.my_has_been_computedbis = false;
            return false;
        }
        let delta_u = (tgtu * tvtv - tgtv * tutv) / dis;
        let delta_v = (tgtv * tutu - tgtu * tutv) / dis;
        tguv1.set_coord(delta_u, delta_v);
        self.my_tguv1 = *tguv1;

        // `gxx:204-226`
        SurfaceTool::d1(&self.surf2, *u2, *v2, &mut pbid, &mut tu, &mut tv);
        let tutu = tu.dot(&tu);
        let tvtv = tv.dot(&tv);
        let tutv = tu.dot(&tv);
        let tgtu = tg.dot(&tu);
        let tgtv = tg.dot(&tv);
        let dis = tutu * tvtv - tutv * tutv;
        if dis.abs() < ANGULAR {
            self.my_is_tangent = false;
            self.my_has_been_computed = false;
            self.my_has_been_computedbis = false;
            return false;
        }
        let delta_u = (tgtu * tvtv - tgtv * tutv) / dis;
        let delta_v = (tgtv * tutu - tgtu * tutv) / dis;
        tguv2.set_coord(delta_u, delta_v);
        self.my_tguv2 = *tguv2;

        true
    }

    /// `Pnt` (`gxx:232-247`).
    fn pnt(&mut self, u1: f64, v1: f64, u2: f64, v2: f64, p: &mut GpPnt) {
        let mut a_p = GpPnt::zero();
        let mut a_t = GpVec::new(0.0, 0.0, 0.0);
        let mut a_ts1 = GpVec2d::new(0.0, 0.0);
        let mut a_ts2 = GpVec2d::new(0.0, 0.0);
        let mut tu1 = u1;
        let mut tu2 = u2;
        let mut tv1 = v1;
        let mut tv2 = v2;
        let _ = self.compute(&mut tu1, &mut tv1, &mut tu2, &mut tv2, &mut a_p, &mut a_t, &mut a_ts1, &mut a_ts2);
        *p = self.my_pnt;
    }

    /// `SeekPoint` (`gxx:251-273`).
    fn seek_point(&mut self, u1: f64, v1: f64, u2: f64, v2: f64, point: &mut PntOn2S) -> bool {
        let mut a_p = GpPnt::zero();
        let mut a_t = GpVec::new(0.0, 0.0, 0.0);
        let mut a_ts1 = GpVec2d::new(0.0, 0.0);
        let mut a_ts2 = GpVec2d::new(0.0, 0.0);
        let mut tu1 = u1;
        let mut tu2 = u2;
        let mut tv1 = v1;
        let mut tv2 = v2;
        if !self.compute(&mut tu1, &mut tv1, &mut tu2, &mut tv2, &mut a_p, &mut a_t, &mut a_ts1, &mut a_ts2)
        {
            return false;
        }
        point.set_value(a_p, tu1, tv1, tu2, tv2);
        true
    }

    /// `Tangency` (`gxx:277-292`).
    fn tangency(&mut self, u1: f64, v1: f64, u2: f64, v2: f64, tg: &mut GpVec) -> bool {
        let mut a_p = GpPnt::zero();
        let mut a_t = GpVec::new(0.0, 0.0, 0.0);
        let mut a_ts1 = GpVec2d::new(0.0, 0.0);
        let mut a_ts2 = GpVec2d::new(0.0, 0.0);
        let mut tu1 = u1;
        let mut tu2 = u2;
        let mut tv1 = v1;
        let mut tv2 = v2;
        let t = self.compute(&mut tu1, &mut tv1, &mut tu2, &mut tv2, &mut a_p, &mut a_t, &mut a_ts1, &mut a_ts2);
        *tg = self.my_tg;
        t
    }

    /// `TangencyOnSurf1` (`gxx:296-312`).
    fn tangency_on_surf1(&mut self, u1: f64, v1: f64, u2: f64, v2: f64, tg: &mut GpVec2d) -> bool {
        let mut a_p = GpPnt::zero();
        let mut a_t = GpVec::new(0.0, 0.0, 0.0);
        let mut a_ts1 = GpVec2d::new(0.0, 0.0);
        let mut a_ts2 = GpVec2d::new(0.0, 0.0);
        let mut tu1 = u1;
        let mut tu2 = u2;
        let mut tv1 = v1;
        let mut tv2 = v2;
        let t = self.compute(&mut tu1, &mut tv1, &mut tu2, &mut tv2, &mut a_p, &mut a_t, &mut a_ts1, &mut a_ts2);
        *tg = self.my_tguv1;
        t
    }

    /// `TangencyOnSurf2` (`gxx:316-332`).
    fn tangency_on_surf2(&mut self, u1: f64, v1: f64, u2: f64, v2: f64, tg: &mut GpVec2d) -> bool {
        let mut a_p = GpPnt::zero();
        let mut a_t = GpVec::new(0.0, 0.0, 0.0);
        let mut a_ts1 = GpVec2d::new(0.0, 0.0);
        let mut a_ts2 = GpVec2d::new(0.0, 0.0);
        let mut tu1 = u1;
        let mut tu2 = u2;
        let mut tv1 = v1;
        let mut tv2 = v2;
        let t = self.compute(&mut tu1, &mut tv1, &mut tu2, &mut tv2, &mut a_p, &mut a_t, &mut a_ts1, &mut a_ts2);
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


