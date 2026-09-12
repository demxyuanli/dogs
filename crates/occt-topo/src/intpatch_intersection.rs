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
        }
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
                self.spnt.push(*p);
            }
        }
    }

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
        pp.perform(s1, d1, s2, d2, tol_tang, tol_arc, 0.1, 0.01);
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
