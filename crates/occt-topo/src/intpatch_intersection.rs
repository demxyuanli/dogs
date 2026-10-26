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

    /// `Perform(S1, D1, S2, D2, TolArc, TolTang)`.
    pub fn perform(
        &mut self,
        s1: &dyn Surface,
        d1: &TopolTool,
        s2: &dyn Surface,
        d2: &TopolTool,
        tol_arc: f64,
        tol_tang: f64,
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
                self.take_impimp(&imp, s1, d1, s2, d2);
                return;
            }
            self.run_prmprm(s1, d1, s2, d2, tol_arc, tol_tang);
            return;
        }
        if ip1 != ip2 {
            let mut ip = ImpPrmIntersection::new();
            ip.perform(s1, d1, s2, d2, tol_arc, tol_tang, 0.01, 0.01);
            if ip.is_done() {
                self.done = true;
                self.empty = ip.is_empty();
                self.slin.extend(ip.lines().iter().cloned());
                return;
            }
        }
        self.run_prmprm(s1, d1, s2, d2, tol_arc, tol_tang);
    }

    fn take_impimp(
        &mut self,
        imp: &ImpImpIntersection,
        s1: &dyn Surface,
        d1: &TopolTool,
        s2: &dyn Surface,
        d2: &TopolTool,
    ) {
        self.done = true;
        self.empty = imp.is_empty();
        self.tangent_faces = imp.tangent_faces();
        self.opposite = imp.opposite_faces();
        self.slin
            .extend(convert_alines(imp.lines(), s1, d1, s2, d2));
        for i in 1..=imp.nb_pnts() {
            if let Some(p) = imp.point(i) {
                self.spnt.push(p.clone());
            }
        }
    }

    /// `ParamParamPerfom(S1, D1, S2, D2, TolArc, TolTang, ListOfPnts, typs1,
    /// typs2)` (`IntPatch_Intersection.cxx:1658-1716`).
    ///
    /// Both domains finite: the `ClearFlag` overload with `ClearFlag = true`
    /// (`:1680`); `ListOfPnts` is empty on every call of this port, so the
    /// `LOfPnts` branch at `:1676-1679` does not run.
    ///
    /// One domain infinite: the 8-argument `Perform`
    /// (`IntPatch_PrmPrmIntersection.cxx:369-928`, the polyhedron body) on the
    /// *trimmed* surface (`:1692`, `:1708`). **UNPORTED**: `FUN_TrimInfSurf`
    /// (`:1712-...`) has no Rust body, so that branch keeps passing the
    /// untrimmed surface.
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
        if !d1.domain_is_infinite() && !d2.domain_is_infinite() {
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

fn convert_alines(
    lines: &[GeomIntLine],
    s1: &dyn Surface,
    d1: &TopolTool,
    s2: &dyn Surface,
    d2: &TopolTool,
) -> Vec<GeomIntLine> {
    let conv = ALineToWLine::new(s1, s2, 200);
    let mut out = Vec::with_capacity(lines.len());
    for line in lines {
        match line {
            GeomIntLine::Analytic(a) => {
                for wl in conv.make_wline(a) {
                    if let Some(purged) = wline_tool::compute_purged_wline(&wl, s1, s2, d1, d2) {
                        if purged.nb_pnts() >= 2 {
                            out.push(GeomIntLine::Walking(purged));
                        }
                    } else if wl.nb_pnts() >= 2 {
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
