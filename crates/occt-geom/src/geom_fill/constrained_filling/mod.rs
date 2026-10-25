//! `GeomFill_ConstrainedFilling` — the transfinite filling with tangency
//! constraints.
//!
//! Source: `GeomFill_ConstrainedFilling.hxx:37-193`,
//! `GeomFill_ConstrainedFilling.cxx:364-1560` (OCCT 8.0.0).
//!
//! Module layout:
//! ```text
//! bound_sorts  the file-static helpers (inqadd, mklin, sortbounds,
//!              coonscnd, killcorners)
//! perform      PerformApprox / MatchKnots / PerformS0 / PerformS1 /
//!              PerformSurface
//! ```
//!
//! Translation notes:
//! - `ptch`, `tgalg[4]`, `curvpol[4]`, `tgtepol[4]`, `S0`, `S1` and `surf` are
//!   `occ::handle`s that start null; the port uses `Option` (`None` until
//!   `Init`/`Build` sets them).
//! - `S0` / `S1` are `NCollection_Array2<gp_Pnt>` indexed `(i, j)` with `i` over
//!   the bound-0 direction; the port keeps the same `s0[i][j]` layout (note
//!   OCCT's `RowLength()` is the column count — `NCollection_Array2.hxx:198`).
//! - `Init` does not clear `tgalg` in OCCT; the port keeps that behaviour and
//!   only assigns the entries the corresponding `Init` overload assigns
//!   (`cxx:430-441` for three bounds, `cxx:531-538` for four).
//! - `SetDomain` compares the boundary handle against `ptch->Bound(I)`; the port
//!   compares the shared pointer with [`std::sync::Arc::ptr_eq`].
//! - `CheckCoonsAlgPatch` / `CheckTgteField` / `CheckApprox` / `CheckResult`
//!   (`cxx:1465-1710`) are debug printers: they write to `std::cout` and draw
//!   through `Law_draw1dcurve`. They are **not ported**; the port draws
//!   nothing.

mod bound_sorts;
mod perform;

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpVec, GpXyz};

use bound_sorts::{coonscnd, kill_corners, sort_bounds};
use crate::bspline_surface::GeomBSplineSurface;
use crate::geom_fill::boundary::{Boundary, DegeneratedBound};
use crate::geom_fill::coons_alg_patch::CoonsAlgPatch;
use crate::geom_fill::corner_state::CornerState;
use crate::geom_fill::tgt_field::{TgtField, TgtOnCoons};
use crate::law::services::law_mix_bnd;

/// `GeomFill_ConstrainedFilling` (`GeomFill_ConstrainedFilling.hxx:53-193`).
pub struct ConstrainedFilling {
    /// `degmax` (`hxx:167`).
    degmax: i32,
    /// `segmax` (`hxx:168`).
    segmax: i32,
    /// `ptch` (`hxx:169`).
    ptch: Option<Arc<CoonsAlgPatch>>,
    /// `tgalg[4]` (`hxx:170`).
    tgalg: [Option<Arc<dyn TgtField>>; 4],
    /// `mig[4]` (`hxx:171`).
    mig: [f64; 4],
    /// `stcor[4]` (`hxx:172`).
    stcor: [CornerState; 4],
    /// `v[4]` (`hxx:173`).
    v: [GpVec; 4],
    /// `appdone` (`hxx:174`).
    appdone: bool,
    /// `degree[2]` (`hxx:175`).
    degree: [i32; 2],
    /// `curvpol[4]` (`hxx:176`).
    curvpol: [Option<Vec<GpPnt>>; 4],
    /// `tgtepol[4]` (`hxx:177`).
    tgtepol: [Option<Vec<GpPnt>>; 4],
    /// `mults[2]` (`hxx:178`).
    mults: [Vec<i32>; 2],
    /// `knots[2]` (`hxx:179`).
    knots: [Vec<f64>; 2],
    /// `ab[4]` (`hxx:180`).
    ab: [Vec<f64>; 4],
    /// `pq[4]` (`hxx:181`).
    pq: [Vec<f64>; 4],
    /// `dom[4]` (`hxx:182`).
    dom: [f64; 4],
    /// `ncpol[4]` (`hxx:183`).
    ncpol: [Option<Vec<GpPnt>>; 4],
    /// `ntpol[4]` (`hxx:184`).
    ntpol: [Option<Vec<GpPnt>>; 4],
    /// `nm[2]` (`hxx:185`).
    nm: [Vec<i32>; 2],
    /// `nk[2]` (`hxx:186`).
    nk: [Vec<f64>; 2],
    /// `ibound[2]` (`hxx:187`).
    ibound: [usize; 2],
    /// `ctr[2]` (`hxx:188`).
    ctr: [i32; 2],
    /// `nbd3` (`hxx:189`).
    nbd3: i32,
    /// `S0` (`hxx:190`).
    s0: Vec<Vec<GpPnt>>,
    /// `S1` (`hxx:191`).
    s1: Vec<Vec<GpPnt>>,
    /// `surf` (`hxx:192`).
    surf: Option<GeomBSplineSurface>,
}

impl ConstrainedFilling {
    /// `GeomFill_ConstrainedFilling(MaxDeg, MaxSeg)` (`cxx:367-380`).
    pub fn new(max_deg: i32, max_seg: i32) -> Self {
        Self {
            degmax: max_deg,
            segmax: max_seg,
            ptch: None,
            tgalg: [None, None, None, None],
            mig: [1.0; 4],
            stcor: [CornerState::new(); 4],
            v: [GpVec::zero(); 4],
            appdone: false,
            degree: [0; 2],
            curvpol: [None, None, None, None],
            tgtepol: [None, None, None, None],
            mults: [Vec::new(), Vec::new()],
            knots: [Vec::new(), Vec::new()],
            ab: [Vec::new(), Vec::new(), Vec::new(), Vec::new()],
            pq: [Vec::new(), Vec::new(), Vec::new(), Vec::new()],
            // `cxx:375`.
            dom: [1.0; 4],
            ncpol: [None, None, None, None],
            ntpol: [None, None, None, None],
            nm: [Vec::new(), Vec::new()],
            nk: [Vec::new(), Vec::new()],
            ibound: [0, 0],
            ctr: [0, 0],
            nbd3: 0,
            s0: Vec::new(),
            s1: Vec::new(),
            surf: None,
        }
    }

    /// The current algorithmic patch (`ptch`), which `Init` always creates.
    fn patch(&self) -> Result<&Arc<CoonsAlgPatch>, String> {
        self.ptch
            .as_ref()
            .ok_or_else(|| "Standard_Failure: GeomFill_ConstrainedFilling not initialized".into())
    }

    /// `GeomFill_ConstrainedFilling::Init(B1, B2, B3, NoCheck)`
    /// (`cxx:382-488`).
    pub fn init_3(
        &mut self,
        b1: Arc<dyn Boundary>,
        b2: Arc<dyn Boundary>,
        b3: Arc<dyn Boundary>,
        no_check: bool,
    ) -> Result<(), String> {
        // `cxx:395-402`.
        let mut bound = vec![b1, b2, b3];
        let mut rev = [false; 3];
        sort_bounds(3, &mut bound, &mut rev, &mut self.stcor)?;

        // `cxx:405`: on reoriente.
        rev[2] = !rev[2];

        // `cxx:413-417`.
        for i in 0..=2 {
            bound[i].reparametrize(0.0, 1.0, false, false, 1.0, 1.0, rev[i])?;
        }

        // `cxx:424-429`: the missing degenerate bound.
        let p1 = bound[1].value(1.0);
        let p2 = bound[2].value(1.0);
        let ppp = GpPnt::from_xyz(&p1.xyz().added(p2.xyz()).multiplied(0.5));
        let t3 = bound[1].tol3d().max(bound[2].tol3d());
        let db: Arc<dyn Boundary> = Arc::new(DegeneratedBound::new(ppp, 0.0, 1.0, t3, 10.0));
        let ptch = Arc::new(CoonsAlgPatch::new(
            bound[0].clone(),
            bound[1].clone(),
            db,
            bound[2].clone(),
        ));
        self.ptch = Some(ptch.clone());

        // `cxx:431-443`: `ttgalg` and `tgalg` share the created fields.
        let mut ttgalg: [Option<Arc<dyn TgtField>>; 4] = [None, None, None, None];
        if bound[0].has_normals() {
            let t: Arc<dyn TgtField> = Arc::new(TgtOnCoons::new(ptch.clone(), 0));
            ttgalg[0] = Some(t.clone());
            self.tgalg[0] = Some(t);
        }
        if bound[1].has_normals() {
            let t: Arc<dyn TgtField> = Arc::new(TgtOnCoons::new(ptch.clone(), 1));
            ttgalg[1] = Some(t.clone());
            self.tgalg[1] = Some(t);
        }
        if bound[2].has_normals() {
            let t: Arc<dyn TgtField> = Arc::new(TgtOnCoons::new(ptch.clone(), 3));
            ttgalg[2] = Some(t.clone());
            self.tgalg[3] = Some(t);
        }

        // `cxx:445-451`.
        for i in 0..=3 {
            self.mig[i] = 1.0;
            if self.tgalg[i].is_some() {
                self.min_tgte(i)?;
            }
        }

        if !no_check {
            // `cxx:456-464`.
            let nrev = [false, false, true];
            self.mig[2] = self.mig[3];
            coonscnd(3, &bound, &nrev, &mut self.stcor, &self.mig)?;
            kill_corners(3, &bound, &rev, &nrev, &self.stcor, &ttgalg)?;
        }
        // `cxx:466`.
        self.stcor[3] = self.stcor[2];

        // `cxx:468-483`: the tangent fields must not twist.
        self.check_and_mix()?;

        self.build()
    }

    /// `GeomFill_ConstrainedFilling::Init(B1, B2, B3, B4, NoCheck)`
    /// (`cxx:491-579`).
    pub fn init_4(
        &mut self,
        b1: Arc<dyn Boundary>,
        b2: Arc<dyn Boundary>,
        b3: Arc<dyn Boundary>,
        b4: Arc<dyn Boundary>,
        no_check: bool,
    ) -> Result<(), String> {
        // `cxx:504-511`.
        let mut bound = vec![b1, b2, b3, b4];
        let mut rev = [false; 4];
        sort_bounds(4, &mut bound, &mut rev, &mut self.stcor)?;

        // `cxx:514-515`.
        rev[2] = !rev[2];
        rev[3] = !rev[3];

        // `cxx:521-525`.
        for i in 0..=3 {
            bound[i].reparametrize(0.0, 1.0, false, false, 1.0, 1.0, rev[i])?;
        }

        // `cxx:531-538`.
        let ptch = Arc::new(CoonsAlgPatch::new(
            bound[0].clone(),
            bound[1].clone(),
            bound[2].clone(),
            bound[3].clone(),
        ));
        self.ptch = Some(ptch.clone());
        for i in 0..=3 {
            if bound[i].has_normals() {
                self.tgalg[i] = Some(Arc::new(TgtOnCoons::new(ptch.clone(), i)));
            }
        }

        // `cxx:542-549`.
        for i in 0..=3 {
            self.mig[i] = 1.0;
            if self.tgalg[i].is_some() {
                self.min_tgte(i)?;
            }
        }

        if !no_check {
            // `cxx:554-559`.
            let nrev = [false, false, true, true];
            coonscnd(4, &bound, &nrev, &mut self.stcor, &self.mig)?;
            kill_corners(4, &bound, &rev, &nrev, &self.stcor, &self.tgalg)?;
        }

        // `cxx:561-576`.
        self.check_and_mix()?;

        self.build()
    }

    /// The `CheckTgte` loop that replaces the blending functions by `Law::MixBnd`
    /// when a tangency field twists (`cxx:468-483` and `:561-576`).
    ///
    /// OCCT `occ::down_cast<Law_Linear>`s the two functions with no fallback, so
    /// a non-linear blending law is a `Standard_Failure` here as well.
    fn check_and_mix(&mut self) -> Result<(), String> {
        for i in 0..=3 {
            self.mig[i] = 1.0;
            if self.tgalg[i].is_some() && !self.check_tgte(i)? {
                let ptch = self.patch()?.clone();
                let (fu1, fu2) = ptch.funcs();
                let lin1 = fu1.as_linear().ok_or_else(|| {
                    "Standard_Failure: GeomFill_ConstrainedFilling Init Law_Linear cast".to_string()
                })?;
                let ffu1 = Arc::new(law_mix_bnd(lin1)?);
                let lin2 = fu2.as_linear().ok_or_else(|| {
                    "Standard_Failure: GeomFill_ConstrainedFilling Init Law_Linear cast".to_string()
                })?;
                let ffu2 = Arc::new(law_mix_bnd(lin2)?);
                ptch.set_func(ffu1, ffu2);
                break;
            }
        }
        Ok(())
    }

    /// `GeomFill_ConstrainedFilling::SetDomain(l, B)` (`cxx:583-604`).
    pub fn set_domain(&mut self, l: f64, b: &Arc<dyn Boundary>) -> Result<(), String> {
        let ptch = self.patch()?.clone();
        let value = 1.0f64.min(l.abs());
        if Arc::ptr_eq(ptch.bound(0), b) {
            self.dom[0] = value;
        } else if Arc::ptr_eq(ptch.bound(1), b) {
            self.dom[1] = value;
        } else if Arc::ptr_eq(ptch.bound(2), b) {
            self.dom[2] = value;
        } else if Arc::ptr_eq(ptch.bound(3), b) {
            self.dom[3] = value;
        }
        Ok(())
    }

    /// `GeomFill_ConstrainedFilling::ReBuild` (`cxx:608-616`).
    pub fn rebuild(&mut self) -> Result<(), String> {
        if !self.appdone {
            return Err(
                "Standard_Failure: GeomFill_ConstrainedFilling::ReBuild Approx non faite".into(),
            );
        }
        self.match_knots()?;
        self.perform_s0()?;
        self.perform_s1()?;
        self.perform_surface()?;
        Ok(())
    }

    /// `GeomFill_ConstrainedFilling::Boundary(I)` (`cxx:622-625`).
    pub fn boundary(&self, i: usize) -> Result<Arc<dyn Boundary>, String> {
        Ok(self.patch()?.bound(i).clone())
    }

    /// `GeomFill_ConstrainedFilling::Surface()` (`cxx:629-632`).
    pub fn surface(&self) -> Option<&GeomBSplineSurface> {
        self.surf.as_ref()
    }

    /// `GeomFill_ConstrainedFilling::Build()` (`cxx:636-690`).
    fn build(&mut self) -> Result<(), String> {
        for count in 0..2usize {
            // `cxx:638-640`.
            self.ibound[0] = count;
            self.ibound[1] = count + 2;
            self.ctr = [0, 0];
            self.nbd3 = 0;
            for ii in 0..2 {
                let b = self.patch()?.bound(self.ibound[ii]).clone();
                if b.has_normals() {
                    self.ctr[ii] = 2;
                } else if !b.is_degenerated() {
                    self.ctr[ii] = 1;
                }
                self.nbd3 += self.ctr[ii];
            }
            if self.nbd3 != 0 {
                self.perform_approx()?;
            }
        }
        // `cxx:666`.
        self.appdone = true;
        self.match_knots()?;
        self.perform_s0()?;
        self.perform_s1()?;
        self.perform_surface()?;
        Ok(())
    }

    /// `GeomFill_ConstrainedFilling::Eval(W, Ord, Result)` (`cxx:1414-1465`).
    ///
    /// `Result` is the `3 * nbd3` flat buffer `AdvApprox_EvaluatorFunction`
    /// hands over; `ctr[ii]` writes the bound value and, at `2`, the tangency
    /// field of each of the two approximated bounds.
    ///
    /// `tgalg[ii]->Value` can raise when the field's boundary has no normals;
    /// `Eval` reports that through a non-zero error code (`cxx:1464` returns
    /// `0` on every path, so the code is the port's only signal).
    pub(crate) fn eval(&self, w: f64, ord: i32, res: &mut [f64]) -> i32 {
        let Some(ptch) = self.ptch.as_ref() else {
            return 1;
        };
        let jmp = 3 * self.ctr[0] as usize;
        match ord {
            0 => {
                if self.ctr[0] != 0 {
                    put_xyz(ptch.bound(self.ibound[0]).value(w).xyz(), res, 0);
                }
                if self.ctr[0] == 2 {
                    let Some(t) = self.tgalg[self.ibound[0]].as_ref() else {
                        return 1;
                    };
                    let Ok(v) = t.value(w) else {
                        return 1;
                    };
                    put_xyz(v.xyz(), res, 3);
                }
                if self.ctr[1] != 0 {
                    put_xyz(ptch.bound(self.ibound[1]).value(w).xyz(), res, jmp);
                }
                if self.ctr[1] == 2 {
                    let Some(t) = self.tgalg[self.ibound[1]].as_ref() else {
                        return 1;
                    };
                    let Ok(v) = t.value(w) else {
                        return 1;
                    };
                    put_xyz(v.xyz(), res, jmp + 3);
                }
            }
            1 => {
                if self.ctr[0] != 0 {
                    let (_, vt) = ptch.bound(self.ibound[0]).d1(w);
                    put_xyz(vt.xyz(), res, 0);
                }
                if self.ctr[0] == 2 {
                    let Some(t) = self.tgalg[self.ibound[0]].as_ref() else {
                        return 1;
                    };
                    let Ok(v) = t.d1(w) else {
                        return 1;
                    };
                    put_xyz(v.xyz(), res, 3);
                }
                if self.ctr[1] != 0 {
                    let (_, vt) = ptch.bound(self.ibound[1]).d1(w);
                    put_xyz(vt.xyz(), res, jmp);
                }
                if self.ctr[1] == 2 {
                    let Some(t) = self.tgalg[self.ibound[1]].as_ref() else {
                        return 1;
                    };
                    let Ok(v) = t.d1(w) else {
                        return 1;
                    };
                    put_xyz(v.xyz(), res, jmp + 3);
                }
            }
            _ => {}
        }
        0
    }

    /// `GeomFill_ConstrainedFilling::CheckTgte(I)` (`cxx:1356-1390`): false
    /// when the frame `(tangent, normal, field)` flips orientation.
    fn check_tgte(&self, i: usize) -> Result<bool, String> {
        let bou = self.patch()?.bound(i).clone();
        if !bou.has_normals() {
            return Ok(true);
        }
        let tgalg = self.tgalg[i]
            .as_ref()
            .ok_or("Standard_Failure: GeomFill_ConstrainedFilling::CheckTgte")?;
        let ll = 1.0 / 12.0;
        let mut pmix = 0.0;
        for iu in 0..13 {
            let uu = iu as f64 * ll;
            let (_, tgte) = bou.d1(uu);
            let norm = bou.norm(uu)?;
            let vfield = tgalg.value(uu)?;
            if iu == 0 {
                pmix = vfield.dot(&tgte.crossed(&norm));
            } else {
                let pmixcur = vfield.dot(&tgte.crossed(&norm));
                if pmix * pmixcur < 0.0 {
                    return Ok(false);
                }
            }
        }
        Ok(true)
    }

    /// `GeomFill_ConstrainedFilling::MinTgte(I)` (`cxx:1394-1412`): the
    /// smallest magnitude of the tangency field along bound `I`.
    fn min_tgte(&mut self, i: usize) -> Result<(), String> {
        let bou = self.patch()?.bound(i).clone();
        if !bou.has_normals() {
            return Ok(());
        }
        let tgalg = self.tgalg[i]
            .as_ref()
            .ok_or("Standard_Failure: GeomFill_ConstrainedFilling::MinTgte")?;
        // `cxx:1400`: `minmag = RealLast()`.
        let mut minmag = f64::MAX;
        let ll = 0.02;
        for iu in 0..=30 {
            let uu = 0.2 + iu as f64 * ll;
            let vv = tgalg.value(uu)?;
            let temp = vv.square_magnitude();
            if temp < minmag {
                minmag = temp;
            }
        }
        // `cxx:1411`.
        self.mig[i] = minmag.sqrt();
        Ok(())
    }
}

/// `res[at..at+3] = xyz` — the `Coord(res[i], res[i+1], res[i+2])` of `Eval`.
fn put_xyz(xyz: &GpXyz, res: &mut [f64], at: usize) {
    if at + 2 < res.len() {
        res[at] = xyz.x();
        res[at + 1] = xyz.y();
        res[at + 2] = xyz.z();
    }
}
