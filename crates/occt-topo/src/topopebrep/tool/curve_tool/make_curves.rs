//! `TopOpeBRepTool_CurveTool::MakeCurves` and the file-local helpers of
//! `TopOpeBRepTool_CurveTool.cxx`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepTool/TopOpeBRepTool_CurveTool.cxx`
//! (OCCT 8.0.0): `CheckApproxResults` `:141-163`, `CheckPCurve` `:166-227`,
//! `MakeCurve3DfromWLineApprox` `:230-251`,
//! `MakeCurve2DfromWLineApproxAndPlane` `:253-274`,
//! `MakeCurve2DfromWLineApprox` `:276-296`, `MakeCurves` `:300-895`.
//!
//! `MakeCurves` re-fits the degree-1 walking approximation of an intersection
//! edge (`C3D`, `PC1`, `PC2`) as B-splines on the two faces' surfaces, with
//! `BRepApprox_Approx` driving the fit over `BRepApprox_ApproxLine`.
//!
//! `#define CurveImprovement` is active (`cxx:68`), so the `nbpol > 100`
//! decimation (`cxx:358-545`) runs. OCCT writes the reduced degree-1 curves
//! back through a `const_cast` of its `const occ::handle&` inputs
//! (`cxx:516-528`, `:535-544`); the port cannot alias its `Arc` inputs, so it
//! answers them in [`MakeCurvesOutcome::improved_c3d`] and friends, and the
//! caller (`TopOpeBRepDS_BuildTool::ApproxCurves`) re-reads them on the failure
//! path exactly where OCCT re-reads `C3D` / `PC1` / `PC2`.
//!
//! UNPORTED: nothing.

use std::sync::Arc;

use occt_core::approx::ApproxParametrizationType;
use occt_core::elib::slib;
use occt_core::gp::{GpDir, GpDir2d, GpLin, GpLin2d, GpPln, GpPnt, GpPnt2d, GpVec};
use occt_core::precision::{APPROXIMATION, CONFUSION};
use occt_geom::app_par_curves::multi_bspline_curve::MultiBSpCurve;
use occt_geom::bspline_curve::GeomBSplineCurve;
use occt_geom::curve::Curve;
use occt_geom::geom_lib_check_bspline::{Check2dBSplineCurve, CheckBSplineCurve};
use occt_geom2d::bspline_curve::Geom2dBSplineCurve;
use occt_geom2d::curve::Curve2d;

use crate::brep_approx::{ApproxLine, BRepApproxApprox, SurfaceTool};
use crate::brep_tool::BRepTool;
use crate::shape::Face;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::CurveTool;

/// `MakeCurves` (`cxx:300-895`) result. OCCT returns the same information
/// through its reference parameters plus the `bool` status.
#[derive(Clone, Default)]
pub struct MakeCurvesOutcome {
    /// `MakeCurves`'s return value.
    pub done: bool,
    /// `C3Dnew` (`cxx:317`).
    pub c3d_new: Option<Arc<dyn Curve>>,
    /// `PC1new` (`cxx:318`).
    pub pc1_new: Option<Arc<dyn Curve2d>>,
    /// `PC2new` (`cxx:319`).
    pub pc2_new: Option<Arc<dyn Curve2d>>,
    /// `TolReached3d` (`cxx:320`); only written on success (`cxx:783`).
    pub tol_reached_3d: f64,
    /// `TolReached2d` (`cxx:321`); only written on success (`cxx:784`).
    pub tol_reached_2d: f64,
    /// The `CurveImprovement` rewrite of `C3D` (`cxx:516-518`), `None` when the
    /// block did not run. OCCT writes it through `const_cast`; see the module
    /// docs.
    pub improved_c3d: Option<Arc<dyn Curve>>,
    /// The `CurveImprovement` rewrite of `PC1` (`cxx:520-524`).
    pub improved_pc1: Option<Arc<dyn Curve2d>>,
    /// The `CurveImprovement` rewrite of `PC2` (`cxx:526-530`).
    pub improved_pc2: Option<Arc<dyn Curve2d>>,
}

impl CurveTool {
    /// `MakeCurves(parmin, parmax, C3D, PC1, PC2, S1, S2, C3Dnew, PC1new,
    /// PC2new, TolReached3d, TolReached2d)` (`cxx:300-895`).
    ///
    /// `s1` / `s2` are `TopoDS::Face(S1)` / `TopoDS::Face(S2)` (`cxx:326-327`
    /// of `CurveTool.cxx` in the `ApproxCurves` caller); the two adaptors are
    /// `BRepAdaptor_Surface(F, compminmaxUV = true)` (`cxx:687-689`).
    pub fn make_curves(
        &self,
        parmin: f64,
        parmax: f64,
        c3d: Option<Arc<dyn Curve>>,
        pc1: Option<Arc<dyn Curve2d>>,
        pc2: Option<Arc<dyn Curve2d>>,
        s1: &Face,
        s2: &Face,
    ) -> MakeCurvesOutcome {
        const TOLCHECK: f64 = 1.0e-7;
        const TOLANGCHECK: f64 = 1.0e-6;

        let mut out = MakeCurvesOutcome::default();

        // cxx:333-337
        let comp_c3d = self.geom_tool().comp_c3d();
        if !comp_c3d {
            return out;
        }

        // cxx:339-343
        let comp_pc1 = self.geom_tool().comp_pc1();
        let comp_pc2 = self.geom_tool().comp_pc2();
        let (mut tol3d, mut tol2d) = self.geom_tool().get_tolerances();
        let mut nb_pnt_max = self.geom_tool().nb_pnt_max();

        // cxx:345-346: `(int)` truncates toward zero, which `as i32` reproduces.
        let iparmin = parmin as i32;
        let mut iparmax = parmax as i32;

        // cxx:329-331: `HC3D` / `HPC1` / `HPC2` down-casts. A null handle is
        // dereferenced at `cxx:551` (`HC3D->NbPoles()`), so the port reports the
        // failure instead of reproducing the crash.
        let Some(c3d) = c3d else {
            return out;
        };
        let Some(h_c3d_poles) = c3d.bspline_poles() else {
            return out;
        };
        let nbpol = h_c3d_poles.len();

        //--------------------- cxx:358-545: `#ifdef CurveImprovement`
        if nbpol > 100 {
            let Some(improved) = decimate_poles(c3d.as_ref(), pc1.as_deref(), pc2.as_deref(), comp_pc1, comp_pc2, tol3d, tol2d) else {
                // `HPC1->Poles()` / `HPC2->Poles()` on a null down-cast throws
                // in OCCT (`cxx:361-363`); the port reports the failure.
                return out;
            };
            out.improved_c3d = Some(improved.c3d);
            out.improved_pc1 = improved.pc1;
            out.improved_pc2 = improved.pc2;
            iparmax = improved.nb_pol as i32;
        }

        //--------------------- cxx:547-549
        // `HC3D` is a handle *copy* taken at `cxx:329`, so OCCT keeps reading
        // the original poles here even when the block above rewrote `C3D`.
        let pol3d: Vec<GpPnt> = h_c3d_poles.to_vec();

        // cxx:552-566: `par` = chord-length parameters of the poles, aborted as
        // `IsBad` when two consecutive poles coincide.
        let mut par = vec![0.0f64; nbpol + 1];
        let mut is_bad = false;
        let mut prev = pol3d[0];
        for ip in 2..=nbpol {
            let dist = prev.distance(&pol3d[ip - 1]);
            if dist < CONFUSION {
                is_bad = true;
                break;
            }
            par[ip] = par[ip - 1] + dist;
            prev = pol3d[ip - 1];
        }

        // cxx:568-650
        if !is_bad {
            let mut curvature = vec![0.0f64; nbpol + 1];

            // cxx:573-611: curvature at the ends from the osculating circle.
            if nbpol > 3 {
                let mut ic = 1usize;
                let mut ip = 2usize;
                while ip <= nbpol - 1 {
                    curvature[ic] = osculating_curvature(&pol3d[ip - 2], &pol3d[ip - 1], &pol3d[ip]);
                    ic = nbpol;
                    ip += nbpol - 3;
                }
            } else if nbpol == 3 {
                let kc = osculating_curvature(&pol3d[0], &pol3d[1], &pol3d[2]);
                curvature[1] = kc;
                curvature[nbpol] = kc;
            } else {
                curvature[1] = 0.0;
                curvature[nbpol] = 0.0;
            }

            // cxx:613-639: interior curvature from the second difference.
            let mut dt = par[2] - par[1];
            let mut dx = (pol3d[1].x() - pol3d[0].x()) / dt;
            let mut dy = (pol3d[1].y() - pol3d[0].y()) / dt;
            let mut dz = (pol3d[1].z() - pol3d[0].z()) / dt;
            for ip in 2..=nbpol - 1 {
                dt = par[ip + 1] - par[ip];
                let dx1 = (pol3d[ip].x() - pol3d[ip - 1].x()) / dt;
                let dy1 = (pol3d[ip].y() - pol3d[ip - 1].y()) / dt;
                let dz1 = (pol3d[ip].z() - pol3d[ip - 1].z()) / dt;
                let d2t = 2.0 / (par[ip + 1] - par[ip - 1]);
                let d2x = d2t * (dx1 - dx);
                let d2y = d2t * (dy1 - dy);
                let d2z = d2t * (dz1 - dz);
                curvature[ip] = (d2x * d2x + d2y * d2y + d2z * d2z).sqrt();
                dx = dx1;
                dy = dy1;
                dz = dz1;
            }

            // cxx:641-667. `crit` is OCCT's own empirical bound (`cxx:641`).
            let crit = 1000.0;
            dt = par[2] - par[1];
            let kcprev = (curvature[2] - curvature[1]) / dt;
            let mut ksi = 0.0;
            for ip in 2..=nbpol - 1 {
                dt = par[ip + 1] - par[ip];
                let kc = (curvature[ip + 1] - curvature[ip]) / dt;
                ksi += (kc - kcprev).abs();
                if ksi > crit {
                    is_bad = true;
                    break;
                }
                // `cxx:661` reads `kc = kcprev;` — OCCT never advances
                // `kcprev`, so it keeps its `cxx:643` value for every
                // iteration. Reproduced as written.
            }
        }

        // cxx:673-680
        let mut degmax = 8;
        let mut parametrization = ApproxParametrizationType::ChordLength;
        if is_bad {
            let tt = (10.0 * tol3d).min(APPROXIMATION);
            tol2d = tt * tol2d / tol3d;
            tol3d = tt;
            nb_pnt_max = 40;
            degmax = 4;
            parametrization = ApproxParametrizationType::Centripetal;
        }

        // cxx:682-689
        let degmin = 4;
        let nitmax = 0;
        let with_tangency = true;
        let bas1 = SurfaceAdaptor::of_face_restricted(s1);
        let bas2 = SurfaceAdaptor::of_face_restricted(s2);

        // cxx:691-692
        let al = ApproxLine::new(Some(c3d.clone()), pc1.clone(), pc2.clone());

        // cxx:694-696
        let mut approx = BRepApproxApprox::new();
        approx.set_parameters(
            tol3d,
            tol2d,
            degmin,
            degmax,
            nitmax,
            nb_pnt_max,
            with_tangency,
            parametrization,
        );

        // cxx:698-714
        let s1_is_plane = bas1.get_type() == crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType::Plane;
        let s2_is_plane = bas2.get_type() == crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType::Plane;
        let plane_first = comp_c3d && comp_pc1 && s1_is_plane;
        let plane_second = comp_c3d && comp_pc2 && s2_is_plane;
        if plane_first {
            // The curve X,Y,Z and U2,V2 is approximated.
            approx.perform(&bas1, &bas2, &al, comp_c3d, false, comp_pc2, iparmin, iparmax);
        } else if plane_second {
            // The curve X,Y,Z and U1,V1 is approximated.
            approx.perform(&bas1, &bas2, &al, comp_c3d, comp_pc1, false, iparmin, iparmax);
        } else {
            approx.perform(&bas1, &bas2, &al, comp_c3d, comp_pc1, comp_pc2, iparmin, iparmax);
        }

        // cxx:716-724
        let mut done = approx.is_done() && check_approx_results(approx.value());

        if done {
            // cxx:726-758
            let amc = approx.value();
            let mut bs3d = None;
            let mut bs2d1 = None;
            let mut bs2d2 = None;
            if plane_first {
                bs3d = make_curve3d_from_wline_approx(amc);
                bs2d1 = make_curve2d_from_wline_approx_and_plane(amc, &SurfaceTool::plane(&bas1));
                if comp_pc2 {
                    bs2d2 = make_curve2d_from_wline_approx(amc, 2);
                }
            } else if plane_second {
                bs3d = make_curve3d_from_wline_approx(amc);
                if comp_pc1 {
                    bs2d1 = make_curve2d_from_wline_approx(amc, 2);
                }
                bs2d2 = make_curve2d_from_wline_approx_and_plane(amc, &SurfaceTool::plane(&bas2));
            } else {
                if comp_c3d {
                    bs3d = make_curve3d_from_wline_approx(amc);
                }
                if comp_pc1 {
                    bs2d1 = make_curve2d_from_wline_approx(amc, 2);
                }
                if comp_pc2 {
                    let i32_index = if comp_pc1 { 3 } else { 2 };
                    bs2d2 = make_curve2d_from_wline_approx(amc, i32_index);
                }
            }

            // cxx:761-768: a curve the fit could not build leaves the outcome
            // null, exactly like a null OCCT handle.
            if (comp_c3d && bs3d.is_none()) || (comp_pc1 && bs2d1.is_none()) || (comp_pc2 && bs2d2.is_none()) {
                done = false;
            }

            // cxx:761-768: the pcurves are checked against the faces bounds.
            if done && comp_pc1 {
                done = bs2d1
                    .as_ref()
                    .map(|c| check_pcurve(c, s1))
                    .unwrap_or(false);
            }
            if done && comp_pc2 {
                done = bs2d2
                    .as_ref()
                    .map(|c| check_pcurve(c, s2))
                    .unwrap_or(false);
            }

            if done {
                // cxx:783-784
                out.tol_reached_3d = approx.tol_reached_3d();
                out.tol_reached_2d = approx.tol_reached_2d();

                // cxx:788-845: `GeomLib_CheckBSplineCurve` /
                // `GeomLib_Check2dBSplineCurve` tangent fixing, with the flags
                // recomputed per curve just as OCCT does.
                if let Some(c) = bs3d.as_mut() {
                    let mut chk = CheckBSplineCurve::new(c, TOLCHECK, TOLANGCHECK);
                    let (bf, bl) = chk.need_tangent_fix();
                    chk.fix_tangent(c, bf, bl);
                }
                if let Some(c) = bs2d1.as_mut() {
                    let mut chk = Check2dBSplineCurve::new(c, TOLCHECK, TOLANGCHECK);
                    let (bf, bl) = chk.need_tangent_fix();
                    chk.fix_tangent(c, bf, bl);
                }
                if let Some(c) = bs2d2.as_mut() {
                    let mut chk = Check2dBSplineCurve::new(c, TOLCHECK, TOLANGCHECK);
                    let (bf, bl) = chk.need_tangent_fix();
                    chk.fix_tangent(c, bf, bl);
                }

                out.c3d_new = bs3d.map(|c| Arc::new(c) as Arc<dyn Curve>);
                out.pc1_new = bs2d1.map(|c| Arc::new(c) as Arc<dyn Curve2d>);
                out.pc2_new = bs2d2.map(|c| Arc::new(c) as Arc<dyn Curve2d>);
            }
        }

        // cxx:770-782: on failure every computed curve is nullified and the
        // tolerances keep the caller's values.
        if !done {
            out.c3d_new = None;
            out.pc1_new = None;
            out.pc2_new = None;
            return out;
        }

        out.done = true;
        out
    }
}

/// The `CurveImprovement` block (`cxx:358-545`) outputs.
struct Decimated {
    c3d: Arc<dyn Curve>,
    pc1: Option<Arc<dyn Curve2d>>,
    pc2: Option<Arc<dyn Curve2d>>,
    nb_pol: usize,
}

/// `#ifdef CurveImprovement` decimation (`cxx:358-545`): drop poles that stay
/// within `tol` of the chord between the previous kept pole and the next one,
/// then rebuild the three curves as degree-1 B-splines over the kept poles.
fn decimate_poles(
    c3d: &dyn Curve,
    pc1: Option<&dyn Curve2d>,
    pc2: Option<&dyn Curve2d>,
    comp_pc1: bool,
    comp_pc2: bool,
    tol3d: f64,
    tol2d: f64,
) -> Option<Decimated> {
    let pol_c3d: Vec<GpPnt> = c3d.bspline_poles()?.to_vec();
    let nbpol = pol_c3d.len();
    let pol_pc1: Option<Vec<GpPnt2d>> = if comp_pc1 {
        let (xs, ys) = pc1?.bspline_poles2d()?;
        Some(xs.iter().zip(ys.iter()).map(|(x, y)| GpPnt2d::new(*x, *y)).collect())
    } else {
        None
    };
    let pol_pc2: Option<Vec<GpPnt2d>> = if comp_pc2 {
        let (xs, ys) = pc2?.bspline_poles2d()?;
        Some(xs.iter().zip(ys.iter()).map(|(x, y)| GpPnt2d::new(*x, *y)).collect())
    } else {
        None
    };

    // cxx:365-369: distances are compared squared.
    let tol = (1.0e-10f64).max(100.0 * tol3d * tol3d);
    let tl2d = tol * (tol2d * tol2d) / (tol3d * tol3d);
    let def = tol;
    let def2d = tol2d;

    // 1-based like OCCT's `NCollection_Array1<bool> IsValid(1, nbpol)`.
    let mut is_valid = vec![true; nbpol + 1];
    let mut p = pol_c3d[0];
    let mut p1 = pol_pc1.as_ref().map(|v| v[0]);
    let mut p2 = pol_pc2.as_ref().map(|v| v[0]);
    let mut nb_pol = 1usize;

    for ip in 2..=nbpol {
        if !is_valid[ip] {
            continue;
        }
        let d = p.square_distance(&pol_c3d[ip - 1]);
        let d1 = match (&p1, &pol_pc1) {
            (Some(a), Some(v)) => a.square_distance(&v[ip - 1]),
            _ => 0.0,
        };
        let d2 = match (&p2, &pol_pc2) {
            (Some(a), Some(v)) => a.square_distance(&v[ip - 1]),
            _ => 0.0,
        };
        if !(d > tol || d1 > tl2d || d2 > tl2d) {
            is_valid[ip] = false;
            continue;
        }

        let dd = if ip < nbpol { p.distance(&pol_c3d[ip]) } else { 0.0 };
        if ip < nbpol && dd < 10.0 * tol {
            // cxx:406-411: distance to the chord `P -> PolC3D(ip + 1)`.
            let d = square_distance_to_line(&p, &pol_c3d[ip], &pol_c3d[ip - 1]);
            let d1 = match (&p1, &pol_pc1) {
                (Some(a), Some(v)) => square_distance_to_line2d(a, &v[ip], &v[ip - 1]),
                _ => 0.0,
            };
            let d2 = match (&p2, &pol_pc2) {
                (Some(a), Some(v)) => square_distance_to_line2d(a, &v[ip], &v[ip - 1]),
                _ => 0.0,
            };
            if d > def || d1 > def2d || d2 > def2d {
                nb_pol += 1;
                p = pol_c3d[ip - 1];
                p1 = pol_pc1.as_ref().map(|v| v[ip - 1]);
                p2 = pol_pc2.as_ref().map(|v| v[ip - 1]);
            } else {
                is_valid[ip] = false;
            }
        } else {
            // `cxx:431` writes `IsValid(ip + 1)`, which is past the array end
            // when `ip == nbpol` (the `dd == 0` / `ip == nbpol` case always
            // lands here). The write is clamped to the last pole: OCCT's
            // out-of-range store is not reproducible.
            if ip + 1 <= nbpol {
                is_valid[ip + 1] = false;
            }
            nb_pol += 1;
            p = pol_c3d[ip - 1];
            p1 = pol_pc1.as_ref().map(|v| v[ip - 1]);
            p2 = pol_pc2.as_ref().map(|v| v[ip - 1]);
        }
    }

    // cxx:451-455
    if nb_pol < 2 {
        is_valid[nbpol] = true;
        nb_pol += 1;
    }

    if nb_pol >= nbpol {
        return None;
    }

    // cxx:457-530
    let mut polc3d: Vec<GpPnt> = vec![GpPnt::zero(); nb_pol];
    let mut polpc1: Vec<GpPnt2d> = vec![GpPnt2d::new(0.0, 0.0); nb_pol];
    let mut polpc2: Vec<GpPnt2d> = vec![GpPnt2d::new(0.0, 0.0); nb_pol];
    let mut knots: Vec<f64> = vec![0.0; nb_pol];
    let mut mults: Vec<i32> = vec![1; nb_pol];
    mults[0] = 2;
    mults[nb_pol - 1] = 2;
    let mut count = 0usize;
    for ip in 1..=nbpol {
        if is_valid[ip] {
            count += 1;
            if count > nb_pol {
                break;
            }
            polc3d[count - 1] = pol_c3d[ip - 1];
            if let Some(v) = &pol_pc1 {
                polpc1[count - 1] = v[ip - 1];
            }
            if let Some(v) = &pol_pc2 {
                polpc2[count - 1] = v[ip - 1];
            }
            knots[count - 1] = count as f64;
        }
    }

    // cxx:503-530: the last pole always comes from `nbpol`.
    polc3d[nb_pol - 1] = pol_c3d[nbpol - 1];
    if let Some(v) = &pol_pc1 {
        polpc1[nb_pol - 1] = v[nbpol - 1];
    }
    if let Some(v) = &pol_pc2 {
        polpc2[nb_pol - 1] = v[nbpol - 1];
    }

    let c3d_new = GeomBSplineCurve::from_poles_knots_mults(polc3d, knots.clone(), mults.clone(), 1).ok()?;
    let pc1_new = if comp_pc1 {
        let xs: Vec<f64> = polpc1.iter().map(|p| p.x()).collect();
        let ys: Vec<f64> = polpc1.iter().map(|p| p.y()).collect();
        let flat = occt_core::bspl::banded_interp::knot_sequence(&knots, &mults, 1);
        Some(Arc::new(Geom2dBSplineCurve::new(xs, ys, flat, 1).ok()?) as Arc<dyn Curve2d>)
    } else {
        None
    };
    let pc2_new = if comp_pc2 {
        let xs: Vec<f64> = polpc2.iter().map(|p| p.x()).collect();
        let ys: Vec<f64> = polpc2.iter().map(|p| p.y()).collect();
        let flat = occt_core::bspl::banded_interp::knot_sequence(&knots, &mults, 1);
        Some(Arc::new(Geom2dBSplineCurve::new(xs, ys, flat, 1).ok()?) as Arc<dyn Curve2d>)
    } else {
        None
    };

    Some(Decimated {
        c3d: Arc::new(c3d_new),
        pc1: pc1_new,
        pc2: pc2_new,
        nb_pol,
    })
}

/// `1. / cir.Radius()` of `gce_MakeCirc(P1, P2, P3)` (`cxx:590-592`): `0.` when
/// the two chords are parallel (`cxx:583-589`), which is the same guard that
/// keeps `gce_MakeCirc` from failing.
///
/// `gp_Circ::Radius()` of the circle through three points is the triangle's
/// circumradius `|P1P2| |P2P3| |P3P1| / (2 |P1P2 x P2P3|)`; the port computes it
/// directly instead of building the `gp_Circ`.
fn osculating_curvature(p1: &GpPnt, p2: &GpPnt, p3: &GpPnt) -> f64 {
    let v1 = GpVec::from_pnts(p1, p2);
    let v2 = GpVec::from_pnts(p2, p3);
    if v1.is_parallel_ang(&v2, 1.0e-4) {
        return 0.0;
    }
    let twice_area = v1.crossed(&v2).magnitude();
    if twice_area <= 0.0 {
        return 0.0;
    }
    let r = v1.magnitude() * v2.magnitude() * p3.distance(p1) / (2.0 * twice_area);
    if r > 0.0 {
        1.0 / r
    } else {
        0.0
    }
}

/// `gp_Lin(P, Q).SquareDistance(R)` (`cxx:409`).
fn square_distance_to_line(p: &GpPnt, q: &GpPnt, r: &GpPnt) -> f64 {
    let v = GpVec::from_pnts(p, q);
    let Ok(dir) = GpDir::new(v.x(), v.y(), v.z()) else {
        return 0.0;
    };
    GpLin::from_pnt_dir(*p, dir).square_distance(r)
}

/// `gp_Lin2d` through `p` and `q`, then its squared distance to `r`
/// (`cxx:412-414`, `:419-421`).
fn square_distance_to_line2d(p: &GpPnt2d, q: &GpPnt2d, r: &GpPnt2d) -> f64 {
    let Ok(dir) = GpDir2d::new(q.x() - p.x(), q.y() - p.y()) else {
        return 0.0;
    };
    let d = GpLin2d::from_pnt_dir(*p, dir).distance(r);
    d * d
}

/// `CheckApproxResults(Approx)` (`cxx:141-163`).
fn check_approx_results(amc: &MultiBSpCurve) -> bool {
    let np = amc.base().nb_poles();
    let nc = amc.base().nb_curves();
    if np < 2 || nc < 1 {
        return false;
    }

    // check the knots for coincidence (`cxx:151-160`)
    let knots = amc.knots();
    for i in 0..knots.len().saturating_sub(1) {
        if knots[i + 1] - knots[i] <= epsilon_of(knots[i]) {
            return false;
        }
    }
    true
}

/// `Epsilon(double)` (`Standard_Real.hxx:242-246`): the distance from
/// `value` to the next representable double towards the infinity of its sign.
fn epsilon_of(value: f64) -> f64 {
    if value >= 0.0 {
        f64::from_bits(value.to_bits() + 1) - value
    } else {
        value - f64::from_bits(value.to_bits() - 1)
    }
}

/// `CheckPCurve(aPC, aFace)` (`cxx:166-227`): sample the pcurve inside
/// `[FirstParameter, LastParameter]` and refuse it when a sample leaves the
/// face's UV bounds grown by 1% of their span.
fn check_pcurve(pc: &dyn Curve2d, face: &Face) -> bool {
    const NPOINTS: i32 = 23;

    let (umin0, umax0, vmin0, vmax0) = BRepTool::uv_bounds(face);
    let tol_u = ((umax0 - umin0) * 0.01).max(CONFUSION);
    let tol_v = ((vmax0 - vmin0) * 0.01).max(CONFUSION);
    let fp = pc.first_parameter();
    let lp = pc.last_parameter();
    let step = (lp - fp) / f64::from(NPOINTS + 1);

    // cxx:190-204: adjust the domain for periodic surfaces around the pcurve's
    // mid parameter.
    let mut umin = umin0;
    let mut umax = umax0;
    let mut vmin = vmin0;
    let mut vmax = vmax0;
    let surface = BRepTool::face_surface_world(face);
    let mid = pc.value((fp + lp) * 0.5);
    if let Some(surf) = surface.as_ref() {
        if surf.is_u_periodic() {
            let per = surf.u_period();
            let mut nshift = ((mid.x() - umin) / per) as i32;
            if mid.x() < umin + per * f64::from(nshift) {
                nshift -= 1;
            }
            umin += per * f64::from(nshift);
            umax += per * f64::from(nshift);
        }
        if surf.is_v_periodic() {
            let per = surf.v_period();
            let mut nshift = ((mid.y() - vmin) / per) as i32;
            if mid.y() < vmin + per * f64::from(nshift) {
                nshift -= 1;
            }
            vmin += per * f64::from(nshift);
            vmax += per * f64::from(nshift);
        }
    }

    for i in 1..=NPOINTS {
        let p = fp + f64::from(i) * step;
        let pnt = pc.value(p);
        let (u, v) = (pnt.x(), pnt.y());
        if umin - u > tol_u || u - umax > tol_u || vmin - v > tol_v || v - vmax > tol_v {
            return false;
        }
    }
    true
}

/// `MakeCurve3DfromWLineApprox(Approx, 1)` (`cxx:230-251`): curve 1 of the fit
/// as a `Geom_BSplineCurve`.
fn make_curve3d_from_wline_approx(amc: &MultiBSpCurve) -> Option<GeomBSplineCurve> {
    let poles = amc.base().curve3d(1);
    GeomBSplineCurve::from_poles_knots_mults(
        poles,
        amc.knots().to_vec(),
        amc.multiplicities().to_vec(),
        amc.degree() as usize,
    )
    .ok()
}

/// `MakeCurve2DfromWLineApproxAndPlane(Approx, Pl)` (`cxx:253-274`): curve 1 of
/// the fit, mapped onto the plane's `(U, V)`.
fn make_curve2d_from_wline_approx_and_plane(
    amc: &MultiBSpCurve,
    pln: &GpPln,
) -> Option<Geom2dBSplineCurve> {
    let poles3d = amc.base().curve3d(1);
    let mut xs = Vec::with_capacity(poles3d.len());
    let mut ys = Vec::with_capacity(poles3d.len());
    for p in &poles3d {
        let (u, v) = slib::plane_parameters(&pln.pos, p);
        xs.push(u);
        ys.push(v);
    }
    let flat = occt_core::bspl::banded_interp::knot_sequence(
        amc.knots(),
        amc.multiplicities(),
        amc.degree(),
    );
    Geom2dBSplineCurve::new(xs, ys, flat, amc.degree() as usize).ok()
}

/// `MakeCurve2DfromWLineApprox(Approx, CI)` (`cxx:276-296`): curve `CI` of the
/// fit as a `Geom2d_BSplineCurve`.
fn make_curve2d_from_wline_approx(amc: &MultiBSpCurve, ci: i32) -> Option<Geom2dBSplineCurve> {
    let poles = amc.base().curve2d(ci);
    let xs: Vec<f64> = poles.iter().map(|p| p.x()).collect();
    let ys: Vec<f64> = poles.iter().map(|p| p.y()).collect();
    let flat = occt_core::bspl::banded_interp::knot_sequence(
        amc.knots(),
        amc.multiplicities(),
        amc.degree(),
    );
    Geom2dBSplineCurve::new(xs, ys, flat, amc.degree() as usize).ok()
}
