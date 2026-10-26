//! `IntPatch_Intersection` — Geom-Geom / Geom-Param / Param-Param dispatch.
//! Source: `IntPatch_Intersection.cxx` Perform / GeomGeomPerfom /
//! GeomParamPerfom / ParamParamPerfom.

use occt_geom::Surface;

use crate::brep_surface::{classify_surface, SurfaceKind};
use crate::geom_int::{GeomIntLine, TopolTool};
use crate::int_tools_wline::PatchPoint;
use crate::intpatch::aline_to_wline::ALineToWLine;
use crate::intpatch::impimp::{try_set_quad, ImpImpIntersection};
use crate::intpatch::wline_tool;

use super::impprm::ImpPrmIntersection;
use super::prmprm::PrmPrmIntersection;
use super::trim_inf_surf;

/// `IntPatch_Intersection`.
#[derive(Clone)]
pub struct PatchIntersection {
    done: bool,
    empty: bool,
    tangent_faces: bool,
    opposite: bool,
    slin: Vec<GeomIntLine>,
    spnt: Vec<PatchPoint>,
    /// `myFleche` (`IntPatch_Intersection.cxx:1080`, `:185`).
    fleche: f64,
    /// `myUVMaxStep` (`:1084`, `:190`).
    uv_max_step: f64,
}

impl PatchIntersection {
    pub fn new() -> Self {
        Self {
            done: false,
            empty: true,
            tangent_faces: false,
            opposite: false,
            slin: Vec::new(),
            spnt: Vec::new(),
            fleche: 0.01,
            uv_max_step: 0.01,
        }
    }

    /// `SetTolerances(TolArc, TolTang, UVMaxStep, Fleche)`
    /// (`IntPatch_Intersection.cxx:136-173`), including its clamps.
    pub fn set_tolerances(&mut self, uv_max_step: f64, fleche: f64) {
        self.fleche = fleche.clamp(1.0e-3, 10.0);
        self.uv_max_step = uv_max_step.min(0.5);
    }

    /// `Perform(S1, D1, S2, D2, TolArc, TolTang)` — the 8-argument overload
    /// with its defaults `isGeomInt = true`, `theIsReqToKeepRLine = false`,
    /// `theIsReqToPostWLProc = true` (`IntPatch_Intersection.hxx:92-94`).
    pub fn perform(
        &mut self,
        s1: &dyn Surface,
        d1: &TopolTool,
        s2: &dyn Surface,
        d2: &TopolTool,
        tol_arc: f64,
        tol_tang: f64,
    ) {
        self.perform_with_options(s1, d1, s2, d2, tol_arc, tol_tang, true)
    }

    /// `Perform(..., theIsReqToPostWLProc)` (`IntPatch_Intersection.cxx:
    /// 1201-1372`): the dispatch, then the post-processing pass of
    /// `:1343-1372`.
    #[allow(clippy::too_many_arguments)]
    pub fn perform_with_options(
        &mut self,
        s1: &dyn Surface,
        d1: &TopolTool,
        s2: &dyn Surface,
        d2: &TopolTool,
        tol_arc: f64,
        tol_tang: f64,
        post_wl_proc: bool,
    ) {
        self.done = false;
        self.empty = true;
        self.tangent_faces = false;
        self.opposite = false;
        self.slin.clear();
        self.spnt.clear();

        let q1 = try_set_quad(s1).is_some();
        let q2 = try_set_quad(s2).is_some();
        let ip1 = is_impprm_quadric(s1);
        let ip2 = is_impprm_quadric(s2);

        if q1 && q2 {
            let mut imp = ImpImpIntersection::new();
            imp.perform(s1, d1, s2, d2, tol_arc, tol_tang);
            if imp.is_done() {
                self.take_impimp(&imp, s1, s2);
            } else {
                self.run_prmprm(s1, d1, s2, d2, tol_arc, tol_tang);
            }
        } else if ip1 != ip2 {
            let mut ip = ImpPrmIntersection::new();
            ip.perform(s1, d1, s2, d2, tol_arc, tol_tang, 0.01, 0.01);
            if ip.is_done() {
                self.done = true;
                self.empty = ip.is_empty();
                self.slin.extend(ip.lines().iter().cloned());
            } else {
                self.run_prmprm(s1, d1, s2, d2, tol_arc, tol_tang);
            }
        } else {
            self.run_prmprm(s1, d1, s2, d2, tol_arc, tol_tang);
        }

        // `if (!theIsReqToPostWLProc) return;` (`cxx:1341-1344`).
        if post_wl_proc {
            self.post_wline_processing(s1, d1, s2, d2);
        }
    }

    /// `for (i = Lower(); i <= slin.Upper(); i++)` of `cxx:1346-1371`
    /// (`:1626-1652` in the `LOfPnts` overload and `:2072-2098` in the
    /// start-point one are the same loop).
    ///
    /// Every `IntPatch_WLine` whose `IsPurgingAllowed()` is `true` is replaced
    /// by `IntPatch_WLineTool::ComputePurgedWLine`; when that returns null the
    /// original line stays. `InsertAfter(i, aRW); Remove(i)` leaves `aRW` at
    /// index `i`, so the loop moves on to the next line.
    fn post_wline_processing(
        &mut self,
        s1: &dyn Surface,
        d1: &TopolTool,
        s2: &dyn Surface,
        d2: &TopolTool,
    ) {
        let mut i = 0usize;
        while i < self.slin.len() {
            let Some(wl) = self.slin[i].as_wline() else {
                i += 1;
                continue;
            };
            if !wl.is_purging_allowed() {
                i += 1;
                continue;
            }
            if let Some(purged) = wline_tool::compute_purged_wline(wl, s1, s2, d1, d2) {
                self.slin[i] = GeomIntLine::Walking(purged);
            }
            i += 1;
        }
    }

    fn take_impimp(&mut self, imp: &ImpImpIntersection, s1: &dyn Surface, s2: &dyn Surface) {
        self.done = true;
        self.empty = imp.is_empty();
        self.tangent_faces = imp.tangent_faces();
        self.opposite = imp.opposite_faces();
        self.slin
            .extend(convert_alines(imp.lines(), s1, s2));
        for i in 1..=imp.nb_pnts() {
            if let Some(p) = imp.point(i) {
                self.spnt.push(p.clone());
            }
        }
    }

    /// `ParamParamPerfom(S1, D1, S2, D2, TolArc, TolTang, ListOfPnts, typs1,
    /// typs2)` (`IntPatch_Intersection.cxx:1658-1768`).
    ///
    /// Both domains finite: the `ClearFlag` overload with `ClearFlag = true`
    /// (`:1680`); `ListOfPnts` is empty on every call of this port, so the
    /// `LOfPnts` branch at `:1676-1679` does not run.
    ///
    /// Exactly one domain infinite: `FUN_GetMinMaxXYZPnt` on the finite surface,
    /// `FUN_TrimInfSurf` on the infinite one (`:1685-1708`), then the
    /// 8-argument `Perform` on the trimmed pair (see
    /// `intpatch_trim_inf_surf.rs`).
    ///
    /// Both domains infinite (`:1710-1750`). **UNPORTED**:
    /// `FUN_PL_Intersection` (`:450-963`) and `FUN_TrimBothSurf` (`:964-...`)
    /// have no Rust body, so this arm keeps passing the untrimmed pair.
    fn run_prmprm(
        &mut self,
        s1: &dyn Surface,
        d1: &TopolTool,
        s2: &dyn Surface,
        d2: &TopolTool,
        tol_arc: f64,
        tol_tang: f64,
    ) {
        let mut pp = PrmPrmIntersection::new();
        let d1_inf = d1.domain_is_infinite();
        let d2_inf = d2.domain_is_infinite();
        if !d1_inf && !d2_inf {
            pp.perform_clear_flag(
                s1,
                d1,
                s2,
                d2,
                tol_tang,
                tol_arc,
                self.fleche,
                self.uv_max_step,
                true,
            );
        } else if d1_inf ^ d2_inf {
            // `:1685-1708`: `AP = max(MU, MV)` of the *finite* surface, and the
            // infinite one is trimmed against the finite one's point cloud.
            let (finite, infinite, infinite_is_s1) = if d1_inf {
                (s2, s1, true)
            } else {
                (s1, s2, false)
            };
            let (u0, u1) = finite.u_range();
            let (v0, v1) = finite.v_range();
            let mu = u0.abs().max(u1.abs());
            let mv = v0.abs().max(v1.abs());
            let ap = mu.max(mv);
            let (pmin, pmax) = trim_inf_surf::get_min_max_xyz_pnt(finite);
            let trimmed = trim_inf_surf::trim_inf_surf(&pmin, &pmax, infinite, ap);
            let (ts1, ts2): (&dyn Surface, &dyn Surface) = if infinite_is_s1 {
                (&*trimmed, s2)
            } else {
                (s1, &*trimmed)
            };
            pp.perform(
                ts1,
                d1,
                ts2,
                d2,
                tol_tang,
                tol_arc,
                self.fleche,
                self.uv_max_step,
            );
        } else {
            pp.perform(s1, d1, s2, d2, tol_tang, tol_arc, self.fleche, self.uv_max_step);
        }
        self.done = pp.is_done();
        self.empty = pp.is_empty();
        self.slin.extend(pp.lines().iter().cloned());
    }

    pub fn is_done(&self) -> bool {
        self.done
    }

    pub fn is_empty(&self) -> bool {
        self.empty
    }

    pub fn tangent_faces(&self) -> bool {
        self.tangent_faces
    }

    pub fn opposite_faces(&self) -> bool {
        self.opposite
    }

    pub fn lines(&self) -> &[GeomIntLine] {
        &self.slin
    }

    pub fn nb_lines(&self) -> i32 {
        self.slin.len() as i32
    }

    pub fn nb_pnts(&self) -> i32 {
        self.spnt.len() as i32
    }
}

impl Default for PatchIntersection {
    fn default() -> Self {
        Self::new()
    }
}

fn convert_alines(lines: &[GeomIntLine], s1: &dyn Surface, s2: &dyn Surface) -> Vec<GeomIntLine> {
    let conv = ALineToWLine::new(s1, s2, 200);
    let mut out = Vec::with_capacity(lines.len());
    for line in lines {
        match line {
            GeomIntLine::Analytic(a) => {
                for wl in conv.make_wline(a) {
                    // `if (aWLine->NbPnts() > 1) theLines.Append(aWLine)`
                    // (`IntPatch_ALineToWLine.cxx:955-960`). The line carries
                    // `EnablePurging(false)` set at `:956`, so the
                    // post-processing pass skips it.
                    if wl.nb_pnts() > 1 {
                        out.push(GeomIntLine::Walking(wl));
                    }
                }
            }
            other => out.push(other.clone()),
        }
    }
    out
}

fn is_impprm_quadric(s: &dyn Surface) -> bool {
    matches!(
        classify_surface(s),
        SurfaceKind::Plane | SurfaceKind::Cylinder | SurfaceKind::Sphere | SurfaceKind::Cone
    )
}
