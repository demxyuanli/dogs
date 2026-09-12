//! `GeomInt_IntSS`. Source: `GeomInt_IntSS.cxx` Perform / InternalPerform.

use std::sync::Arc;

use occt_core::gp::GpPnt2d;
use occt_core::precision::CONFUSION;
use occt_geom::{Curve, Surface};
use occt_geom2d::curve::Curve2d;

use crate::brep_surface::{classify_surface, SurfaceKind};
use crate::fclass2d::FaceState;
use crate::geom_int::intss_make;
use crate::geom_int::line_ctor::LineConstructor;
use crate::geom_int::topol::TopolTool;
use crate::geom_int::types::GeomIntLine;
use crate::int_tools_wline::WLine;
use crate::intpatch::{self, IntersectionCurve, PatchIntersection, SurfaceIntersection};

/// One `GeomInt_IntSS` result line (3D + pcurves on S1 / S2).
#[derive(Clone)]
pub struct IntSSLine {
    pub curve: Arc<dyn Curve>,
    pub pcurve1: Option<Arc<dyn Curve2d>>,
    pub pcurve2: Option<Arc<dyn Curve2d>>,
}

/// Surface/surface intersection (`GeomInt_IntSS`).
pub struct IntSS {
    pub(crate) hs1: Option<Arc<dyn Surface>>,
    pub(crate) hs2: Option<Arc<dyn Surface>>,
    pub(crate) l_construct: LineConstructor,
    pub(crate) sline: Vec<IntSSLine>,
    pub(crate) nbrestr: i32,
    pub(crate) tol_reached_2d: f64,
    pub(crate) tol_reached_3d: f64,
    pub(crate) same_surfaces: bool,
    done: bool,
    tangent_faces: bool,
}

impl IntSS {
    pub fn new() -> Self {
        Self {
            hs1: None,
            hs2: None,
            l_construct: LineConstructor::new(),
            sline: Vec::new(),
            nbrestr: 0,
            tol_reached_2d: 0.0,
            tol_reached_3d: 0.0,
            same_surfaces: false,
            done: false,
            tangent_faces: false,
        }
    }

    /// `Load` domains and surfaces onto the line constructor.
    pub fn load(
        &mut self,
        d1: TopolTool,
        d2: TopolTool,
        s1: Arc<dyn Surface>,
        s2: Arc<dyn Surface>,
    ) {
        self.same_surfaces = Arc::ptr_eq(&s1, &s2);
        self.hs1 = Some(s1.clone());
        self.hs2 = Some(s2.clone());
        self.l_construct.load(d1, d2, s1, s2);
    }

    /// `GeomInt_IntSS::Perform(S1, S2, Tol, Approx, ApproxS1, ApproxS2)`.
    pub fn perform(
        &mut self,
        s1: Arc<dyn Surface>,
        s2: Arc<dyn Surface>,
        tol: f64,
        approx: bool,
        approx_s1: bool,
        approx_s2: bool,
    ) {
        let d1 = TopolTool::from_surface(s1.as_ref());
        let d2 = TopolTool::from_surface(s2.as_ref());
        self.load(d1, d2, s1, s2);
        self.internal_perform(tol, approx, approx_s1, approx_s2, false, 0.0, 0.0, 0.0, 0.0);
    }

    /// `Perform` with a starting UV on each surface.
    pub fn perform_with_start(
        &mut self,
        s1: Arc<dyn Surface>,
        s2: Arc<dyn Surface>,
        tol: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        approx: bool,
        approx_s1: bool,
        approx_s2: bool,
    ) {
        let d1 = TopolTool::from_surface(s1.as_ref());
        let d2 = TopolTool::from_surface(s2.as_ref());
        self.load(d1, d2, s1, s2);
        self.internal_perform(tol, approx, approx_s1, approx_s2, true, u1, v1, u2, v2);
    }

    /// `Perform` using domains already installed by [`Self::load`].
    pub fn perform_loaded(
        &mut self,
        tol: f64,
        approx: bool,
        approx_s1: bool,
        approx_s2: bool,
    ) {
        self.internal_perform(tol, approx, approx_s1, approx_s2, false, 0.0, 0.0, 0.0, 0.0);
    }

    fn internal_perform(
        &mut self,
        tol: f64,
        approx: bool,
        approx_s1: bool,
        approx_s2: bool,
        use_start: bool,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) {
        self.tol_reached_2d = 0.0;
        self.tol_reached_3d = 0.0;
        self.nbrestr = 0;
        self.sline.clear();
        self.done = false;
        self.tangent_faces = false;

        let Some(hs1) = self.hs1.clone() else {
            return;
        };
        let Some(hs2) = self.hs2.clone() else {
            return;
        };
        let _tol_arc = tol;
        let _tol_tang = tol;
        let mut _deflection = 0.1;
        if classify_surface(hs1.as_ref()) == SurfaceKind::Other
            && classify_surface(hs2.as_ref()) == SurfaceKind::Other
        {
            _deflection /= 10.0;
        }
        let _ = _deflection;

        let (d1, d2) = if let Some(pair) = self.l_construct.domain_clones() {
            pair
        } else {
            let d1 = TopolTool::from_surface(hs1.as_ref());
            let d2 = TopolTool::from_surface(hs2.as_ref());
            self.l_construct
                .load(d1.clone(), d2.clone(), hs1.clone(), hs2.clone());
            (d1, d2)
        };

        if use_start {
            let st1 = d1.classify(GpPnt2d::new(u1, v1), tol);
            let st2 = d2.classify(GpPnt2d::new(u2, v2), tol);
            if st1 == FaceState::Out || st2 == FaceState::Out {
                self.done = true;
                return;
            }
        }

        if self.same_surfaces {
            match intpatch::surface_surface_intersection(hs1.as_ref(), hs1.as_ref(), tol) {
                SurfaceIntersection::Curves(ics) => {
                    self.emit_curves(&ics, approx, approx_s1, approx_s2, tol);
                }
                SurfaceIntersection::Coincident | SurfaceIntersection::None => {}
            }
        } else {
            let mut inter = PatchIntersection::new();
            inter.perform(hs1.as_ref(), &d1, hs2.as_ref(), &d2, tol, tol);
            if inter.is_done() {
                self.tangent_faces = inter.tangent_faces();
                if !inter.tangent_faces() {
                    for line in inter.lines() {
                        self.make_curve(line, approx, approx_s1, approx_s2, tol);
                    }
                }
                self.done = true;
                return;
            }
        }
        self.done = true;
        let _ = CONFUSION;
    }

    fn emit_curves(
        &mut self,
        ics: &[IntersectionCurve],
        approx: bool,
        approx_s1: bool,
        approx_s2: bool,
        tol: f64,
    ) {
        for ic in ics {
            let line = line_from_intersection_curve(ic);
            self.make_curve(&line, approx, approx_s1, approx_s2, tol);
        }
    }

    /// `GeomInt_IntSS::MakeCurve`.
    pub fn make_curve(
        &mut self,
        line: &GeomIntLine,
        approx: bool,
        approx_s1: bool,
        approx_s2: bool,
        tol: f64,
    ) {
        intss_make::make_curve(self, line, approx, approx_s1, approx_s2, tol);
    }

    pub(crate) fn bump_tol_2d(&mut self, tolpc: f64) {
        if tolpc > self.tol_reached_2d || self.tol_reached_2d == 0.0 {
            self.tol_reached_2d = tolpc;
        }
    }

    pub(crate) fn bump_tol_3d(&mut self, t: f64) {
        if t > self.tol_reached_3d || self.tol_reached_3d == 0.0 {
            self.tol_reached_3d = t;
        }
    }

    pub(crate) fn append_line(
        &mut self,
        curve: Arc<dyn Curve>,
        pcurve1: Option<Arc<dyn Curve2d>>,
        pcurve2: Option<Arc<dyn Curve2d>>,
    ) {
        self.sline.push(IntSSLine {
            curve,
            pcurve1,
            pcurve2,
        });
    }

    pub fn is_done(&self) -> bool {
        self.done
    }

    pub fn tangent_faces(&self) -> bool {
        self.tangent_faces
    }

    pub fn nb_lines(&self) -> i32 {
        self.sline.len() as i32
    }

    pub fn line(&self, index: i32) -> Option<&IntSSLine> {
        self.sline.get((index as usize).saturating_sub(1))
    }

    pub fn lines(&self) -> &[IntSSLine] {
        &self.sline
    }

    pub fn has_line_on_s1(&self, index: i32) -> bool {
        self.line(index)
            .map(|l| l.pcurve1.is_some())
            .unwrap_or(false)
    }

    pub fn has_line_on_s2(&self, index: i32) -> bool {
        self.line(index)
            .map(|l| l.pcurve2.is_some())
            .unwrap_or(false)
    }

    pub fn line_on_s1(&self, index: i32) -> Option<&Arc<dyn Curve2d>> {
        self.line(index).and_then(|l| l.pcurve1.as_ref())
    }

    pub fn line_on_s2(&self, index: i32) -> Option<&Arc<dyn Curve2d>> {
        self.line(index).and_then(|l| l.pcurve2.as_ref())
    }

    pub fn nb_boundaries(&self) -> i32 {
        self.nbrestr
    }

    pub fn boundary(&self, index: i32) -> Option<&IntSSLine> {
        if index < 1 || index > self.nbrestr {
            return None;
        }
        self.sline.get((index as usize).saturating_sub(1))
    }

    pub fn tol_reached_2d(&self) -> f64 {
        self.tol_reached_2d
    }

    pub fn tol_reached_3d(&self) -> f64 {
        self.tol_reached_3d
    }
}

impl Default for IntSS {
    fn default() -> Self {
        Self::new()
    }
}

fn line_from_intersection_curve(ic: &IntersectionCurve) -> GeomIntLine {
    GeomIntLine::Walking(WLine::from_intersection_curve(ic))
}
