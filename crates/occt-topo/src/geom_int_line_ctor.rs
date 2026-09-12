//! `GeomInt_LineConstructor`. Source: `GeomInt_LineConstructor.cxx` Perform
//! walking / dispatch, plus `Load` / `Part`.

use std::sync::Arc;

use occt_core::gp::GpPnt2d;
use occt_core::precision::{PCONFUSION, Precision};
use occt_geom::Surface;

use crate::brep_surface::{classify_surface, SurfaceKind};
use crate::fclass2d::FaceState;
use crate::geom_int::line_ctor_gline;
use crate::geom_int::line_ctor_rline;
use crate::geom_int::line_tool;
use crate::geom_int::quadric::{adjust_periodic_uv, surface_parameters};
use crate::geom_int::topol::TopolTool;
use crate::geom_int::types::{GeomIntLine, IntPatchIType};
use crate::int_tools_wline::WLineWay;

/// Splits an intersection line into in-domain parts.
pub struct LineConstructor {
    done: bool,
    seqp: Vec<f64>,
    dom1: Option<TopolTool>,
    dom2: Option<TopolTool>,
    hs1: Option<Arc<dyn Surface>>,
    hs2: Option<Arc<dyn Surface>>,
}

impl LineConstructor {
    pub fn new() -> Self {
        Self {
            done: false,
            seqp: Vec::new(),
            dom1: None,
            dom2: None,
            hs1: None,
            hs2: None,
        }
    }

    pub fn load(
        &mut self,
        d1: TopolTool,
        d2: TopolTool,
        s1: Arc<dyn Surface>,
        s2: Arc<dyn Surface>,
    ) {
        self.dom1 = Some(d1);
        self.dom2 = Some(d2);
        self.hs1 = Some(s1);
        self.hs2 = Some(s2);
    }

    pub fn is_done(&self) -> bool {
        self.done
    }

    pub fn nb_parts(&self) -> i32 {
        (self.seqp.len() / 2) as i32
    }

    /// 1-based part `[WFirst, WLast]`.
    pub fn part(&self, i: i32) -> Option<(f64, f64)> {
        if !self.done || i < 1 {
            return None;
        }
        let idx = (2 * i - 1) as usize;
        if idx + 1 > self.seqp.len() {
            return None;
        }
        Some((self.seqp[idx - 1], self.seqp[idx]))
    }

    pub(crate) fn seqp_mut(&mut self) -> &mut Vec<f64> {
        &mut self.seqp
    }

    pub(crate) fn set_done(&mut self, d: bool) {
        self.done = d;
    }

    pub(crate) fn surfaces(&self) -> Option<(&dyn Surface, &dyn Surface)> {
        Some((self.hs1.as_ref()?.as_ref(), self.hs2.as_ref()?.as_ref()))
    }

    pub(crate) fn domains(&self) -> Option<(&TopolTool, &TopolTool)> {
        Some((self.dom1.as_ref()?, self.dom2.as_ref()?))
    }

    pub(crate) fn surface_arcs(&self) -> Option<(Arc<dyn Surface>, Arc<dyn Surface>)> {
        Some((self.hs1.clone()?, self.hs2.clone()?))
    }

    pub(crate) fn domain_clones(&self) -> Option<(TopolTool, TopolTool)> {
        Some((self.dom1.clone()?, self.dom2.clone()?))
    }

    /// `GeomInt_LineConstructor::Perform`.
    pub fn perform(&mut self, line: &GeomIntLine) {
        let tol = PCONFUSION * 35.0;
        let typl = line.arc_type();
        if typl == IntPatchIType::Analytic {
            self.perform_analytic(line, tol);
            return;
        }
        if typl == IntPatchIType::Walking {
            self.perform_walking(line, tol);
            return;
        }
        if typl != IntPatchIType::Restriction {
            line_ctor_gline::perform_gline(self, line, tol);
            return;
        }
        line_ctor_rline::perform_restriction(self, line, tol);
    }

    fn perform_analytic(&mut self, line: &GeomIntLine, tol: f64) {
        self.seqp.clear();
        let Some(aline) = line.as_aline() else {
            self.done = true;
            return;
        };
        let Some((s1, s2)) = self.surface_arcs() else {
            self.done = false;
            return;
        };
        let Some((d1, d2)) = self.domain_clones() else {
            self.done = false;
            return;
        };
        let s1 = s1.as_ref();
        let s2 = s2.as_ref();
        let nbvtx = line_tool::nb_vertex(line);
        for i in 1..nbvtx {
            let firstp = line_tool::vertex(line, i).parameter_on_line();
            let lastp = line_tool::vertex(line, i + 1).parameter_on_line();
            if firstp != lastp {
                let pmid = (firstp + lastp) * 0.5;
                let p = aline.value(pmid);
                if let Some((mut u1, mut v1)) = surface_parameters(s1, &p) {
                    if let Some((mut u2, mut v2)) = surface_parameters(s2, &p) {
                        adjust_periodic_uv(s1, s2, &mut u1, &mut v1, &mut u2, &mut v2);
                        if d1.classify(GpPnt2d::new(u1, v1), tol) != FaceState::Out
                            && d2.classify(GpPnt2d::new(u2, v2), tol) != FaceState::Out
                        {
                            self.seqp.push(firstp);
                            self.seqp.push(lastp);
                        }
                    }
                }
            }
        }
        self.done = true;
    }

    fn perform_walking(&mut self, line: &GeomIntLine, tol: f64) {
        self.seqp.clear();
        let Some(wline) = line.as_wline() else {
            self.done = true;
            return;
        };
        let Some((s1, s2)) = self.surface_arcs() else {
            self.done = false;
            return;
        };
        let Some((d1, d2)) = self.domain_clones() else {
            self.done = false;
            return;
        };
        let s1 = s1.as_ref();
        let s2 = s2.as_ref();
        let nbvtx = line_tool::nb_vertex(line);
        for i in 1..nbvtx {
            let firstp = line_tool::vertex(line, i).parameter_on_line();
            let lastp = line_tool::vertex(line, i + 1).parameter_on_line();
            if firstp == lastp {
                continue;
            }
            if lastp != firstp + 1.0 {
                let pmid = ((firstp + lastp) / 2.0) as i32;
                let p = wline.point(pmid);
                let (mut u1, mut v1, mut u2, mut v2) = p.parameters();
                adjust_periodic_uv(s1, s2, &mut u1, &mut v1, &mut u2, &mut v2);
                if d1.classify(GpPnt2d::new(u1, v1), tol) != FaceState::Out
                    && d2.classify(GpPnt2d::new(u2, v2), tol) != FaceState::Out
                {
                    self.seqp.push(firstp);
                    self.seqp.push(lastp);
                }
            } else if wline.creating_way() == WLineWay::ImpPrm {
                let pf = wline.point(firstp as i32);
                let pl = wline.point(lastp as i32);
                let (mut u1, mut v1, mut u2, mut v2) = pf.parameters();
                adjust_periodic_uv(s1, s2, &mut u1, &mut v1, &mut u2, &mut v2);
                let (mut a_u21, mut a_v21, mut a_u22, mut a_v22) = pl.parameters();
                adjust_periodic_uv(s1, s2, &mut a_u21, &mut a_v21, &mut a_u22, &mut a_v22);
                u1 = 0.5 * (u1 + a_u21);
                v1 = 0.5 * (v1 + a_v21);
                u2 = 0.5 * (u2 + a_u22);
                v2 = 0.5 * (v2 + a_v22);
                if d1.classify(GpPnt2d::new(u1, v1), tol) != FaceState::Out
                    && d2.classify(GpPnt2d::new(u2, v2), tol) != FaceState::Out
                {
                    self.seqp.push(firstp);
                    self.seqp.push(lastp);
                }
            } else {
                let pf = wline.point(firstp as i32);
                let (mut u1, mut v1, mut u2, mut v2) = pf.parameters();
                adjust_periodic_uv(s1, s2, &mut u1, &mut v1, &mut u2, &mut v2);
                if d1.classify(GpPnt2d::new(u1, v1), tol) != FaceState::Out
                    && d2.classify(GpPnt2d::new(u2, v2), tol) != FaceState::Out
                {
                    let pl = wline.point(lastp as i32);
                    let (mut u1, mut v1, mut u2, mut v2) = pl.parameters();
                    adjust_periodic_uv(s1, s2, &mut u1, &mut v1, &mut u2, &mut v2);
                    if d1.classify(GpPnt2d::new(u1, v1), tol) != FaceState::Out
                        && d2.classify(GpPnt2d::new(u2, v2), tol) != FaceState::Out
                    {
                        self.seqp.push(firstp);
                        self.seqp.push(lastp);
                    }
                }
            }
        }
        let a_nb_parts = self.seqp.len() / 2;
        if a_nb_parts > 1 {
            let a_st1 = classify_surface(s1);
            let a_st2 = classify_surface(s2);
            let mut b_cond = false;
            if a_st1 == SurfaceKind::Plane
                && matches!(a_st2, SurfaceKind::Other)
            {
                b_cond = true;
            } else if a_st2 == SurfaceKind::Plane && matches!(a_st1, SurfaceKind::Other) {
                b_cond = true;
            }
            if b_cond {
                merge_connected_walking_parts(&mut self.seqp);
            }
        }
        let _ = Precision::CONFUSION;
        self.done = true;
    }
}

impl Default for LineConstructor {
    fn default() -> Self {
        Self::new()
    }
}

/// PKV 22.Apr.2002 block: collapse duplicate indices in `seqp`.
fn merge_connected_walking_parts(seqp: &mut Vec<f64>) {
    let mut a_map: Vec<i32> = Vec::new();
    let mut a_seq_tmp: Vec<f64> = Vec::new();
    for &lastp in seqp.iter() {
        let an_index = lastp as i32;
        if !a_map.contains(&an_index) {
            a_map.push(an_index);
            a_seq_tmp.push(lastp);
        } else if !a_seq_tmp.is_empty() {
            a_seq_tmp.pop();
        }
    }
    seqp.clear();
    let a_nb = a_seq_tmp.len() / 2;
    for i in 1..=a_nb {
        let jx = 2 * i;
        seqp.push(a_seq_tmp[jx - 2]);
        seqp.push(a_seq_tmp[jx - 1]);
    }
}
