//! `TopOpeBRep_FacesIntersector` — intersection of two faces as the
//! `TopOpeBRep` boolean engine sees it.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRep/TopOpeBRep_FacesIntersector.hxx:35-140`,
//! `.cxx:155-650`.
//!
//! Translation notes:
//! - `myIntersector` is [`PatchIntersection`]; `mySurface1` / `mySurface2` are
//!   [`SurfaceAdaptor`] (the port's `BRepAdaptor_Surface` slice) and stay
//!   `Option` because `BRepAdaptor_Surface` has no default state here. They are
//!   `Some` from the first `Perform` on.
//! - `myDomain1->Initialize(aSurf1)` (`cxx:194-197`) resolves the face from the
//!   adaptor (`BRepTopAdaptor_TopolTool.cxx:71-96`), forces it `FORWARD` and
//!   builds one p-curve per face edge; the port passes the already-`FORWARD`
//!   face to [`TopolTool::initialize`], which does the same.
//! - `Lines()`'s `myHAL` is a `NCollection_HArray1<TopOpeBRep_LineInter>` whose
//!   bounds move (allocated `0..n` at `cxx:379-380`, rebuilt `1..myLineNb` at
//!   `cxx:458`); [`LinesArray`] keeps the bounds next to the values.
//! - `GLOBAL_tolFF` (`cxx:213`) is a file-scope exported global that no ported
//!   caller of this class reads yet; it is not reproduced.
//! - `myLine` (`hxx:126`) is never touched by the translated methods and is
//!   not carried.
//! - The `newV == false` arm of `PrepareLines` is dead source (`cxx:391-401`
//!   assigns `newV = true` unconditionally and the other arm is commented
//!   out); only the `newV == true` arm is translated.

pub mod merged_wline;
pub mod rline_build;
pub mod wlines_to_arc;

use std::sync::Arc;

use occt_core::bnd::BndBox;
use occt_core::precision;
use occt_geom::Surface;

use crate::abs::{Orientation, ShapeType};
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::geom_int::GeomIntLine;
use crate::intpatch::PatchIntersection;
use crate::intpatch_line_ctor::LineConstructor;
use crate::shape::{Face, TopoShape};
use crate::topexp::Explorer;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;
use crate::topopebrep::tool::shape_tool::{surfaces_same_oriented, tolerance_of};
use crate::topopebrep::tool::tol::ftol_face_tolerances3d_boxes;
use crate::topopebrep::{LineInter, TypeLineCurve};

use wlines_to_arc::{merge_wlines_if_all_segments_along_restriction, test_wlines_to_an_arc};

/// `NCollection_HArray1<TopOpeBRep_LineInter>` (`hxx:125`).
pub struct LinesArray {
    /// `Lower()` of the OCCT array.
    pub lower: i32,
    /// `Upper()` of the OCCT array.
    pub upper: i32,
    /// The elements, `values[i - lower]`.
    pub values: Vec<LineInter>,
}

impl LinesArray {
    /// `new NCollection_HArray1<TopOpeBRep_LineInter>(aLower, anUpper)`: every
    /// element default-constructed.
    pub fn alloc(lower: i32, upper: i32) -> Self {
        let n = if upper >= lower { (upper - lower + 1) as usize } else { 0 };
        Self {
            lower,
            upper,
            values: (0..n).map(|_| LineInter::new()).collect(),
        }
    }

    /// 1-based `Value(i)` / `ChangeValue(i)`.
    pub fn value(&self, i: i32) -> Option<&LineInter> {
        let k = i - self.lower;
        if k < 0 {
            return None;
        }
        self.values.get(k as usize)
    }

    /// `ChangeValue(i)`.
    pub fn value_mut(&mut self, i: i32) -> Option<&mut LineInter> {
        let k = i - self.lower;
        if k < 0 {
            return None;
        }
        self.values.get_mut(k as usize)
    }
}

/// `class TopOpeBRep_FacesIntersector` (`hxx:35-140`).
pub struct FacesIntersector {
    /// `myIntersector` (`hxx:119`).
    intersector: PatchIntersection,
    /// `myIntersectionDone` (`hxx:120`).
    intersection_done: bool,
    /// `myTol1` (`hxx:121`).
    tol1: f64,
    /// `myTol2` (`hxx:122`).
    tol2: f64,
    /// `myForceTolerances` (`hxx:123`).
    force_tolerances: bool,
    /// `myHAL` (`hxx:125`).
    hal: Option<LinesArray>,
    /// `myLineIndex` (`hxx:127`).
    line_index: i32,
    /// `myLineFound` (`hxx:128`).
    line_found: bool,
    /// `myLineNb` (`hxx:129`).
    line_nb: i32,
    /// `myFace1` (`hxx:130`).
    face1: Option<Face>,
    /// `myFace2` (`hxx:131`).
    face2: Option<Face>,
    /// `mySurface1` (`hxx:132`).
    surface1: Option<SurfaceAdaptor>,
    /// `mySurface2` (`hxx:133`).
    surface2: Option<SurfaceAdaptor>,
    /// `mySurfaceType1` (`hxx:134`).
    surface_type1: GeomAbsSurfaceType,
    /// `mySurfaceType2` (`hxx:135`).
    surface_type2: GeomAbsSurfaceType,
    /// `mySurfacesSameOriented` (`hxx:136`).
    surfaces_same_oriented: bool,
    /// `myDomain1` (`hxx:137`).
    domain1: TopolTool,
    /// `myDomain2` (`hxx:138`).
    domain2: TopolTool,
    /// `myEdgeRestrictionMap` (`hxx:139`).
    edge_restriction_map: Vec<TopoShape>,
    /// `GLOBAL_tolFF` (`TopOpeBRep_FacesIntersector.cxx:37:1.e-7`,
    /// set to `max(tol1, tol2)` at `:212`). A file-scope export in OCCT that
    /// `TopOpeBRep_ProcessGR.cxx:301` reads through
    /// `TopOpeBRep_FacesFiller::LSameDomainERL`; kept on the intersector that
    /// owns the tolerances.
    global_tol_ff: f64,
}

impl Default for FacesIntersector {
    fn default() -> Self {
        Self::new()
    }
}

impl FacesIntersector {
    /// `TopOpeBRep_FacesIntersector()` (`cxx:155-165`).
    pub fn new() -> Self {
        let mut me = Self {
            intersector: PatchIntersection::new(),
            intersection_done: false,
            tol1: 0.0,
            tol2: 0.0,
            force_tolerances: false,
            hal: None,
            line_index: 1,
            line_found: false,
            line_nb: 0,
            face1: None,
            face2: None,
            surface1: None,
            surface2: None,
            surface_type1: GeomAbsSurfaceType::OtherSurface,
            surface_type2: GeomAbsSurfaceType::OtherSurface,
            surfaces_same_oriented: false,
            domain1: TopolTool::new(),
            domain2: TopolTool::new(),
            edge_restriction_map: Vec::new(),
            global_tol_ff: 1.0e-7,
        };
        me.reset_intersection();
        me.tol1 = precision::CONFUSION;
        me.tol2 = precision::CONFUSION;
        me.force_tolerances = false;
        me
    }

    /// `Perform(F1, F2, B1, B2)` (`cxx:166-266`).
    pub fn perform(&mut self, f1: &Face, f2: &Face, b1: &BndBox, b2: &BndBox) {
        self.reset_intersection();
        if !self.force_tolerances {
            self.shape_tolerances(f1, f2);
        }

        let mut face1 = f1.clone();
        face1.0.set_orientation(Orientation::Forward);
        let mut face2 = f2.clone();
        face2.0.set_orientation(Orientation::Forward);
        self.face1 = Some(face1.clone());
        self.face2 = Some(face2.clone());

        let s1 = SurfaceAdaptor::of_face(&face1);
        let s2 = SurfaceAdaptor::of_face(&face2);
        self.surface_type1 = s1.get_type();
        self.surface_type2 = s2.get_type();
        let a_surf1 = s1.surface().cloned();
        let a_surf2 = s2.surface().cloned();
        if let Some(a) = &a_surf1 {
            if self.domain1.initialize(a.clone(), &face1).is_err() {
                return;
            }
        }
        if let Some(a) = &a_surf2 {
            if self.domain2.initialize(a.clone(), &face2).is_err() {
                return;
            }
        }

        let mut deflection = 0.01;
        let mut max_uv = 0.01;
        if !self.force_tolerances {
            ftol_face_tolerances3d_boxes(
                b1,
                b2,
                &face1,
                &face2,
                &s1,
                &s2,
                &mut self.tol1,
                &mut self.tol2,
                &mut deflection,
                &mut max_uv,
            );
            self.tol1 = if self.tol1 > 1.0e-4 { 1.0e-4 } else { self.tol1 };
            self.tol2 = if self.tol2 > 1.0e-4 { 1.0e-4 } else { self.tol2 };
        }

        let tol1 = self.tol1;
        let tol2 = self.tol2;
        self.global_tol_ff = if tol1 > tol2 { tol1 } else { tol2 };
        self.surface1 = Some(s1.clone());
        self.surface2 = Some(s2.clone());

        let (Some(a1), Some(a2)) = (a_surf1, a_surf2) else {
            return;
        };
        // OCCT hands the very same `BRepTopAdaptor_TopolTool` to
        // `IntPatch_Intersection` through its `Adaptor3d_TopolTool` base.
        // The port's `PatchIntersection` takes the `geom_int` adaptor view
        // instead, so the classifier is rebuilt from the same face that
        // `myDomain1` / `myDomain2` were initialized with.
        let d1 = crate::geom_int::TopolTool::from_face(&face1, a1.as_ref());
        let d2 = crate::geom_int::TopolTool::from_face(&face2, a2.as_ref());
        self.intersector
            .perform(a1.as_ref(), &d1, a2.as_ref(), &d2, tol1, tol2);

        // xpu180998 : cto900Q1 — the intersection may be incomplete.
        if !self.intersector.is_done() {
            return;
        }

        self.prepare_lines();
        self.intersection_done = true;

        // mySurfacesSameOriented : a mettre dans IntPatch NYI
        if self.intersector.tangent_faces() {
            self.surfaces_same_oriented = surfaces_same_oriented(&s1, &s2);
        }

        // build the map of edges found as RESTRICTION
        let mut arcs: Vec<TopoShape> = Vec::new();
        self.init_line();
        while self.more_line() {
            if self.current_line().type_line_curve() == TypeLineCurve::Restriction {
                if let Some(a) = self.current_line().arc() {
                    arcs.push(a.clone());
                }
            }
            self.next_line();
        }
        for a in arcs {
            self.edge_restriction_map.push(a);
        }
    }

    /// `Perform(F1, F2)` (`cxx:268-272`).
    pub fn perform_faces(&mut self, f1: &Face, f2: &Face) {
        let b1 = BndBox::new();
        let b2 = BndBox::new();
        self.perform(f1, f2, &b1, &b2);
    }

    /// `IsEmpty()` (`cxx:276-302`).
    pub fn is_empty(&mut self) -> bool {
        if !self.intersection_done {
            return false;
        }
        let done = self.intersector.is_done();
        let empty = self.intersector.is_empty();
        if !done || empty {
            return true;
        }
        // ElemIntersector is done and is not empty: returns True if at least one
        // VPoint is found.
        let mut empty = true;
        self.init_line();
        while self.more_line() {
            empty = self.current_line().nb_v_point() == 0;
            if !empty {
                break;
            }
            self.next_line();
        }
        empty
    }

    /// `IsDone()` (`cxx:307-310`).
    pub fn is_done(&self) -> bool {
        self.intersection_done
    }

    /// `SameDomain()` (`cxx:314-325`): raises `Standard_ProgramError` before a
    /// successful `Perform`.
    pub fn same_domain(&self) -> Result<bool, String> {
        if !self.intersection_done {
            return Err("FacesIntersector : bad SameDomain".to_string());
        }
        Ok(self.intersector.tangent_faces())
    }

    /// `Face(Index)` (`cxx:331-346`): raises for an index other than 1 or 2.
    pub fn face(&self, index: i32) -> Result<&Face, String> {
        match index {
            1 => self
                .face1
                .as_ref()
                .ok_or_else(|| "TopOpeBRep_FacesIntersector::Face".to_string()),
            2 => self
                .face2
                .as_ref()
                .ok_or_else(|| "TopOpeBRep_FacesIntersector::Face".to_string()),
            _ => Err("TopOpeBRep_FacesIntersector::Face".to_string()),
        }
    }

    /// `GLOBAL_tolFF` (`FacesIntersector.cxx:37`, `:212`): the larger of the
    /// two working tolerances, `1.e-7` before the first `Perform`.
    pub fn global_tol_ff(&self) -> f64 {
        self.global_tol_ff
    }

    /// `SurfacesSameOriented()` (`cxx:349-357`).
    pub fn surfaces_same_oriented(&self) -> Result<bool, String> {
        if self.same_domain()? {
            return Ok(self.surfaces_same_oriented);
        }
        Err("FacesIntersector : bad SurfacesSameOriented".to_string())
    }

    /// `IsRestriction(E)` (`cxx:360-364`).
    pub fn is_restriction(&self, e: &TopoShape) -> bool {
        self.edge_restriction_map
            .iter()
            .any(|r| crate::topo_tools::is_same(r, e))
    }

    /// `Restrictions()` (`cxx:368-372`).
    pub fn restrictions(&self) -> &[TopoShape] {
        &self.edge_restriction_map
    }

    /// `PrepareLines()` (`cxx:376-485`).
    pub fn prepare_lines(&mut self) {
        self.line_nb = 0;
        let n = self.intersector.nb_lines();
        self.hal = Some(LinesArray::alloc(0, n));
        if n == 0 {
            return;
        }

        let tol1 = self.tol1;
        let a_surf1 = self.surface1.as_ref().and_then(|a| a.surface().cloned());
        let a_surf2 = self.surface2.as_ref().and_then(|a| a.surface().cloned());
        let (Some(a1), Some(a2)) = (a_surf1, a_surf2) else {
            return;
        };
        let a1: Arc<dyn Surface> = a1;
        let a2: Arc<dyn Surface> = a2;
        let sl: Vec<GeomIntLine> = self.intersector.lines().to_vec();

        let Self {
            domain1,
            domain2,
            hal,
            line_nb,
            ..
        } = self;

        let mut seq_of_result_lines: Vec<GeomIntLine> = Vec::new();
        for i in 1..=n {
            let mut line_constructor = LineConstructor::new(2);
            line_constructor.perform(
                &sl,
                &sl[(i - 1) as usize],
                a1.as_ref(),
                domain1,
                a2.as_ref(),
                domain2,
                tol1,
            );
            let mut seq_of_lines: Vec<GeomIntLine> = line_constructor.lines().to_vec();
            test_wlines_to_an_arc(
                &mut seq_of_lines,
                a1.as_ref(),
                domain1,
                a2.as_ref(),
                domain2,
                tol1,
            );
            seq_of_result_lines.extend(seq_of_lines);
        }

        *line_nb = seq_of_result_lines.len() as i32;

        merge_wlines_if_all_segments_along_restriction(
            &mut seq_of_result_lines,
            a1.as_ref(),
            domain1,
            a2.as_ref(),
            domain2,
            tol1,
        );
        *line_nb = seq_of_result_lines.len() as i32;

        if *line_nb > 0 {
            let mut values = Vec::with_capacity(seq_of_result_lines.len());
            for (index, l) in seq_of_result_lines.iter().enumerate() {
                let mut li = LineInter::new();
                li.set_line(l, a1.as_ref(), a2.as_ref());
                li.set_index((index + 1) as i32);
                values.push(li);
            }
            *hal = Some(LinesArray {
                lower: 1,
                upper: *line_nb,
                values,
            });
        }
    }

    /// `Lines()` (`cxx:487-490`).
    pub fn lines(&self) -> Option<&LinesArray> {
        self.hal.as_ref()
    }

    /// `NbLines()` (`cxx:494-497`).
    pub fn nb_lines(&self) -> i32 {
        self.line_nb
    }

    /// `InitLine()` (`cxx:501-505`).
    pub fn init_line(&mut self) {
        self.line_index = 1;
        self.find_line();
    }

    /// `FindLine()` (`cxx:509-531`).
    fn find_line(&mut self) {
        self.line_found = false;
        if !self.intersection_done {
            return;
        }
        while self.line_index <= self.line_nb {
            let ok = match self.hal.as_ref().and_then(|h| h.value(self.line_index)) {
                Some(l) => l.ok(),
                None => false,
            };
            self.line_found = ok;
            if ok {
                break;
            }
            self.line_index += 1;
        }
    }

    /// `MoreLine()` (`cxx:534-537`).
    pub fn more_line(&self) -> bool {
        self.line_found
    }

    /// `NextLine()` (`cxx:541-545`).
    pub fn next_line(&mut self) {
        self.line_index += 1;
        self.find_line();
    }

    /// `CurrentLine()` (`cxx:549-553`).
    pub fn current_line(&self) -> &LineInter {
        self.hal
            .as_ref()
            .and_then(|h| h.value(self.line_index))
            .expect("FacesIntersector :: CurrentLine")
    }

    /// `CurrentLineIndex()` (`cxx:557-560`).
    pub fn current_line_index(&self) -> i32 {
        self.line_index
    }

    /// `ChangeLine(IL)` (`cxx:564-568`).
    pub fn change_line(&mut self, il: i32) -> &mut LineInter {
        self.hal
            .as_mut()
            .and_then(|h| h.value_mut(il))
            .expect("FacesIntersector :: ChangeLine")
    }

    /// Two aliasing `ChangeLine` borrows, as `TopOpeBRep_FacesFiller`'s
    /// `VP_PositionOnL` / `PequalVPonR` walk the other lines of `FI` while
    /// holding the current one (`TopOpeBRep_FacesFiller.cxx:436-455`,
    /// `:626-663`). OCCT reaches them through the raw pointer it stores in
    /// `myFacesIntersector`; the port hands out a disjoint pair instead.
    ///
    /// `(None, None)` when either index is out of range or when both name the
    /// same line, which is exactly the `iOL == Lindex` skip of `cxx:441-444`.
    pub fn two_lines_mut(
        &mut self,
        i: i32,
        j: i32,
    ) -> (Option<&mut LineInter>, Option<&mut LineInter>) {
        if i == j {
            return (None, None);
        }
        let Some(hal) = self.hal.as_mut() else {
            return (None, None);
        };
        let ki = i - hal.lower;
        let kj = j - hal.lower;
        let n = hal.values.len() as i32;
        if ki < 0 || kj < 0 || ki >= n || kj >= n {
            return (None, None);
        }
        if ki < kj {
            let (a, b) = hal.values.split_at_mut(kj as usize);
            (Some(&mut a[ki as usize]), Some(&mut b[0]))
        } else {
            let (a, b) = hal.values.split_at_mut(ki as usize);
            (Some(&mut b[0]), Some(&mut a[kj as usize]))
        }
    }

    /// `ResetIntersection()` (`cxx:572-580`).
    fn reset_intersection(&mut self) {
        self.intersection_done = false;
        self.line_index = 1;
        self.line_nb = 0;
        self.edge_restriction_map.clear();
        self.line_found = false;
    }

    /// `ForceTolerances(Tol1, Tol2)` (`cxx:584-593`).
    pub fn force_tolerances(&mut self, tol1: f64, tol2: f64) {
        self.tol1 = tol1;
        self.tol2 = tol2;
        self.force_tolerances = true;
    }

    /// `GetTolerances(Tol1, Tol2)` (`cxx:597-601`).
    pub fn get_tolerances(&self) -> (f64, f64) {
        (self.tol1, self.tol2)
    }

    /// `ShapeTolerances(S1, S2)` (`cxx:605-626`): the `ToleranceMax` line is
    /// commented out in the source, so both tolerances are `Precision::Confusion`.
    fn shape_tolerances(&mut self, _s1: &TopoShape, _s2: &TopoShape) {
        self.tol1 = precision::CONFUSION;
        self.tol2 = self.tol1;
        self.force_tolerances = false;
    }

    /// `ToleranceMax(S, T)` (`cxx:628-646`).
    pub fn tolerance_max(&self, s: &TopoShape, t: ShapeType) -> f64 {
        let mut e = Explorer::new(s, t);
        if !e.more() {
            return precision::INTERSECTION;
        }
        let mut tol = f64::MIN;
        while e.more() {
            let v = tolerance_of(e.current()).unwrap_or(precision::CONFUSION);
            tol = if v > tol { v } else { tol };
            e.next();
        }
        tol
    }
}
