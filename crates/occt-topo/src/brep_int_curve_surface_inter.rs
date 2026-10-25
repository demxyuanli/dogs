//! `BRepIntCurveSurface_Inter` — a curve walked against every face of a shape.
//!
//! Source: `ModelingAlgorithms/TKTopAlgo/BRepIntCurveSurface/BRepIntCurveSurface_Inter.hxx:44-147`,
//! `.cxx:30-262`.
//!
//! Translation notes:
//! - `myCurve` is a `GeomAdaptor_Curve` handle; the port stores the
//!   `Arc<dyn Curve>` the adaptor wraps and reads `FirstParameter()` /
//!   `LastParameter()` from it.
//! - `myIntcs.Perform(myCurve, aSurfForFastClass)` (`cxx:179`) is
//!   [`crate::intcurvesurface::perform_curve_surface`]; `BRepAdaptor_Surface(F)`
//!   is built with the default `Restriction = true` (`BRepAdaptor_Surface.hxx`),
//!   so the UV window is the face box, [`finite_uv`].
//! - `BndLib_Add3dCurve::Add(*myCurve, 0., myCurveBox)` (`cxx:88`) is
//!   [`crate::int_tools_curve_box::add_curve_to_box`] with a zero tolerance.
//! - `aFaceBox.IsOut(myCurve->Line())` (`cxx:172`) is
//!   [`occt_core::bnd::box3d::BndBox::is_out_line`].
//! - `myFastClass->Initialize(aSurf)` + `Classify(Puv, Tol)` (`cxx:185-186`,
//!   `:216`) is [`crate::geom_int_topol::TopolTool`] built from the face.
//! - `StdFail_NotDone` (`cxx:216`, `:229`, ...) is returned as `Err`, since the
//!   accessors are total in the port.
//! - `myIndFace` is 1-based in OCCT; the port keeps that so `More()` reads
//!   `ind_face <= faces.len()` exactly.

use std::sync::Arc;

use occt_core::bnd::box3d::BndBox;
use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::precision::Precision;
use occt_geom::{Curve, Surface};

use crate::abs::State;
use crate::bbox_from_geometry::shape_bbox;
use crate::brep_tool::BRepTool;
use crate::fclass2d::FaceState;
use crate::geom_int::TopolTool;
use crate::int_curves_face::{finite_uv, Transition};
use crate::int_tools_curve_box::add_curve_to_box;
use crate::intcurvesurface::{self, HInterResult, IntersectionPoint};
use crate::shape::Face;
use crate::topo_tools::faces_of;

/// `class BRepIntCurveSurface_Inter` (`hxx:44-147`).
pub struct BRepIntCurveSurfaceInter {
    /// `myTolerance` (`hxx:124`).
    tolerance: f64,
    /// `myCurve` (`hxx:125`).
    curve: Option<Arc<dyn Curve>>,
    /// `myIntcs` (`hxx:126`): the result of the last `Perform`.
    intcs: HInterResult,
    /// `myCurrentindex` (`hxx:127`), 0 when there is no current point.
    current_index: usize,
    /// `myCurrentnbpoints` (`hxx:128`).
    current_nb_points: usize,
    /// `myFastClass` (`hxx:129`).
    fast_class: Option<TopolTool>,
    /// `myCurrentstate` (`hxx:130`).
    current_state: FaceState,
    /// `myCurrentU` (`hxx:131`).
    current_u: f64,
    /// `myCurrentV` (`hxx:132`).
    current_v: f64,
    /// `myCurveBox` (`hxx:133`).
    curve_box: BndBox,
    /// `myCurve->FirstParameter()` (`BRepIntCurveSurface_Inter.cxx:88`).
    curve_first: f64,
    /// `myCurve->LastParameter()` (`BRepIntCurveSurface_Inter.cxx:89`).
    curve_last: f64,
    /// `myIndFace` (`hxx:134`), 1-based; 0 before the first face.
    ind_face: usize,
    /// `myFaces` (`hxx:135`).
    faces: Vec<Face>,
    /// `myFaceBoxes` (`hxx:136`): `None` is the null `HArray1`.
    face_boxes: Option<Vec<BndBox>>,
}

impl Default for BRepIntCurveSurfaceInter {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepIntCurveSurfaceInter {
    /// `BRepIntCurveSurface_Inter()` (`cxx:30-35`): `myFastClass` is built and
    /// `Clear()` runs.
    pub fn new() -> Self {
        let mut i = Self {
            tolerance: 0.0,
            curve: None,
            intcs: HInterResult::default(),
            current_index: 0,
            current_nb_points: 0,
            fast_class: None,
            current_state: FaceState::Unknown,
            current_u: 0.0,
            current_v: 0.0,
            curve_box: BndBox::new(),
            curve_first: 0.0,
            curve_last: 0.0,
            ind_face: 0,
            faces: Vec::new(),
            face_boxes: None,
        };
        i.clear();
        i
    }

    /// `Clear()` (`cxx:60-68`).
    fn clear(&mut self) {
        self.current_index = 0;
        self.current_nb_points = 0;
        self.ind_face = 0;
        self.current_state = FaceState::Unknown;
        self.current_u = 0.0;
        self.current_v = 0.0;
    }

    /// `Load(S, Tol)` (`cxx:72-83`): reset and collect the `FACE` sub-shapes.
    pub fn load(&mut self, shape: &crate::shape::TopoShape, tol: f64) {
        self.clear();
        self.faces.clear();
        self.face_boxes = None;
        self.tolerance = tol;
        self.faces = faces_of(shape);
    }

    /// `Init(S, C, Tol)` (`cxx:38-44`).
    pub fn init_shape_curve(&mut self, shape: &crate::shape::TopoShape, curve: Arc<dyn Curve>, tol: f64) {
        self.load(shape, tol);
        self.init_curve(curve);
    }

    /// `Init(S, L, Tol)` (`cxx:48-56`): the line is wrapped into a
    /// `Geom_Line` / `GeomAdaptor_Curve` before `Init(C)`.
    pub fn init_shape_line(
        &mut self,
        shape: &crate::shape::TopoShape,
        line: &occt_core::gp::GpLin,
        tol: f64,
    ) {
        let curve = Arc::new(occt_geom::GeomLine::new(*line));
        self.load(shape, tol);
        self.init_curve(curve as Arc<dyn Curve>);
    }

    /// `Init(C)` (`cxx:87-98`): (re)start the walk of `myCurve` over the
    /// loaded faces, over the curve's own range.
    pub fn init_curve(&mut self, curve: Arc<dyn Curve>) {
        let first = curve.first_parameter();
        let last = curve.last_parameter();
        self.init_curve_range(curve, first, last);
    }

    /// `myCurve->Load(C, First, Last)` followed by `Find()` (`cxx:87-98`). The
    /// range is the adaptor's, which for `GeomAdaptor_Curve` loaded from
    /// `BRep_Tool::Curve(E, loc, f, l)` is the edge range `[f, l]`, not the
    /// underlying curve's.
    pub fn init_curve_range(&mut self, curve: Arc<dyn Curve>, first: f64, last: f64) {
        self.clear();
        self.curve_box.set_void();
        self.curve = Some(curve);
        self.curve_first = first;
        self.curve_last = last;
        if !Precision::is_infinite(first) && !Precision::is_infinite(last) {
            let c = self.curve.as_ref().expect("curve set above");
            add_curve_to_box(c.as_ref(), first, last, 0.0, &mut self.curve_box);
        }
        self.find();
    }

    /// `More()` (`cxx:102-105`).
    pub fn more(&self) -> bool {
        self.ind_face <= self.faces.len()
    }

    /// `Next()` (`cxx:109-116`).
    pub fn next(&mut self) {
        if self.current_nb_points > 0 {
            self.current_index += 1;
        }
        self.find();
    }

    /// `Find()` (`cxx:120-190`).
    fn find(&mut self) {
        if self.current_nb_points > 0 && self.current_index <= self.current_nb_points && self.find_point()
        {
            return;
        }

        self.current_nb_points = 0;
        self.current_index = 0;

        let mut i = self.ind_face + 1;
        while i <= self.faces.len() {
            let face = self.faces[i - 1].clone();
            if self.face_boxes.is_none() {
                self.face_boxes = Some(vec![BndBox::new(); self.faces.len()]);
            }
            let boxes = self.face_boxes.as_mut().expect("boxes allocated above");
            let face_box = &mut boxes[i - 1];
            if face_box.is_void() {
                *face_box = shape_bbox(&face.0);
                face_box.set_gap(self.tolerance);
            }
            let is_out = {
                let Some(curve) = self.curve.as_ref() else {
                    return;
                };
                match curve.gp_line() {
                    Some(lin) => face_box.is_out_line(&lin),
                    None => {
                        if !self.curve_box.is_void() {
                            face_box.is_out_box(&self.curve_box)
                        } else {
                            false
                        }
                    }
                }
            };
            if is_out {
                i += 1;
                continue;
            }
            let Some(surf) = BRepTool::face_surface(&face) else {
                i += 1;
                continue;
            };
            let (first, last) = (self.curve_first, self.curve_last);
            let (u0, u1, v0, v1) = finite_uv(&face);
            let Some(curve) = self.curve.clone() else {
                return;
            };
            let res = intcurvesurface::perform_curve_surface(
                curve.as_ref(),
                surf.as_ref(),
                (first, last),
                (u0, u1, v0, v1),
            );
            self.intcs = match res {
                Ok(r) => r,
                Err(_) => HInterResult::default(),
            };
            self.current_nb_points = self.intcs.nb_points();
            if self.current_nb_points == 0 {
                i += 1;
                continue;
            }
            self.fast_class = Some(TopolTool::from_face(&face, surf.as_ref()));
            self.ind_face = i;
            if self.find_point() {
                return;
            }
            self.current_nb_points = 0;
            i += 1;
        }

        if self.current_nb_points == 0 && i > self.faces.len() {
            self.ind_face = i;
        }
    }

    /// `FindPoint()` (`cxx:194-213`): the first point of the current face that
    /// the fast classifier places `IN` or `ON`.
    fn find_point(&mut self) -> bool {
        let mut j = if self.current_index == 0 {
            1
        } else {
            self.current_index
        };

        while j <= self.current_nb_points {
            let p = *self.intcs.point(j - 1);
            let u = p.u();
            let v = p.v();
            let uv = GpPnt2d::new(u, v);

            let st = match self.fast_class.as_ref() {
                Some(tool) => tool.classify(uv, self.tolerance),
                None => FaceState::Unknown,
            };
            self.current_state = st;
            if st == FaceState::On || st == FaceState::In {
                self.current_index = j;
                self.current_u = u;
                self.current_v = v;
                return true;
            }
            j += 1;
        }
        false
    }

    /// `Point()` (`cxx:217-232`): the current point with the classified
    /// `(U, V)`.
    pub fn point(&self) -> Result<IntersectionPoint, String> {
        if self.current_index == 0 {
            return Err("StdFail_NotDone".to_string());
        }
        let p = *self.intcs.point(self.current_index - 1);
        Ok(IntersectionPoint::new(
            p.param(),
            self.current_u,
            self.current_v,
            p.pnt(),
            p.state(),
        ))
    }

    /// `Pnt()` (`cxx:253-262`).
    pub fn pnt(&self) -> Result<GpPnt, String> {
        if self.current_index == 0 {
            return Err("StdFail_NotDone".to_string());
        }
        Ok(self.intcs.point(self.current_index - 1).pnt())
    }

    /// `U()` (`cxx:235-243`).
    pub fn u(&self) -> Result<f64, String> {
        if self.current_index == 0 {
            return Err("StdFail_NotDone".to_string());
        }
        Ok(self.current_u)
    }

    /// `V()` (`cxx:246-253`).
    pub fn v(&self) -> Result<f64, String> {
        if self.current_index == 0 {
            return Err("StdFail_NotDone".to_string());
        }
        Ok(self.current_v)
    }

    /// `W()` (`cxx:256-263`).
    pub fn w(&self) -> Result<f64, String> {
        if self.current_index == 0 {
            return Err("StdFail_NotDone".to_string());
        }
        Ok(self.intcs.point(self.current_index - 1).param())
    }

    /// `State()` (`cxx:266-274`): `IN` or `ON`.
    pub fn state(&self) -> Result<State, String> {
        if self.current_index == 0 {
            return Err("StdFail_NotDone".to_string());
        }
        Ok(match self.current_state {
            FaceState::In => State::In,
            FaceState::On => State::On,
            FaceState::Out => State::Out,
            FaceState::Unknown => State::Unknown,
        })
    }

    /// `Transition()` (`cxx:277-283`).
    pub fn transition(&self) -> Result<Transition, String> {
        if self.current_index == 0 {
            return Err("StdFail_NotDone".to_string());
        }
        Ok(match self.intcs.point(self.current_index - 1).state() {
            intcurvesurface::State::In => Transition::In,
            intcurvesurface::State::Out => Transition::Out,
            intcurvesurface::State::On | intcurvesurface::State::Unknown => Transition::Tangent,
        })
    }

    /// `Face()` (`cxx:286-289`).
    pub fn face(&self) -> Option<&Face> {
        if self.ind_face == 0 {
            return None;
        }
        self.faces.get(self.ind_face - 1)
    }
}
