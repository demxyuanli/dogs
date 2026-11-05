//! `BRepFill_Filling` — the N-side filling algorithm.
//!
//! Source: `TKBool/BRepFill/BRepFill_Filling.{hxx,cxx}` (OCCT 8.0.0).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | ctor `(Degree, NbPtsOnCur, NbIter, Anisotropie, Tol2d, Tol3d, TolAng, TolCurv, MaxDeg, MaxSegments)` | `cxx:142-171` | [`Filling::new`] |
//! | `SetConstrParam` / `SetResolParam` / `SetApproxParam` | `cxx:177-207` | the same names |
//! | `LoadInitSurface` | `cxx:209-213` | [`Filling::load_init_surface`] |
//! | the five `Add` overloads | `cxx:217-299` | `add` / `add_with_support` / `add_free` / `add_point` / `add_point_on_face` |
//! | `AddConstraints(SeqOfConstraints)` | `cxx:304-389` | [`Filling::add_constraints`] |
//! | `BuildWires` / `FindExtremitiesOfHoles` | `cxx:392-567` | [`build_wires`], [`find_extremities_of_holes`] |
//! | `Build` | `cxx:571-777` | `build.rs` |
//! | `IsDone` / `Face` / `Generated` | `cxx:782-806` | the same names |
//! | `G0Error` / `G1Error` / `G2Error` (+ the `Index` overloads) | `cxx:812-865` | the same names |
//!
//! # Translation notes
//!
//! * `myBuilder` is a `std::shared_ptr<GeomPlate_BuildPlateSurface>` reset at the
//!   start of `Build` (`cxx:574-582`), so it is an `Option` here; the error
//!   accessors answer `0.0` before `Build` where OCCT dereferences a null
//!   pointer.
//! * `myBoundary` / `myConstraints` / `myFreeConstraints` / `myPoints` are
//!   `NCollection_Sequence`s appended to by the `Add` overloads; they become
//!   `Vec`s with the same 1-based return values (the returned `int` is
//!   `Length()`, i.e. the index of the constraint just added).
//! * `myOldNewMap` keys on `TopTools_ShapeMapHasher` (TShape + location,
//!   orientation ignored) — see [`OldNewMap`].
//! * `BRepAdaptor_Curve(E)` is [`BRepTool::edge_curve_adaptor_world`] and
//!   `BRep_Tool::Surface(Face)` is [`BRepTool::face_surface_world`]; the port
//!   stores located geometry, so no `Transform(loc.Transformation())` copy step
//!   is needed on either.
//! * `Geom_BuildPlateSurface::Curves2d()` answers `Option<Arc<dyn Curve2d>>` per
//!   bound; a bound whose pcurve is missing keeps OCCT's
//!   `UpdateEdge(..., null, ...)` arm, i.e. the edge gets no pcurve on the
//!   result surface.
//!
//! # Gaps
//!
//! `BRepFill_ComputeCLine` (`BRepFill_ComputeCLine_0.cxx:19-22`, the
//! `Approx_ComputeCLine.gxx` instantiation over `BRepFill_MultiLine`) is
//! **UNPORTED**: OCCT 8.0 only `#include`s its header
//! (`BRepFill_TrimSurfaceTool.cxx:20`) and never instantiates the class, and the
//! port has no `approx::LineTool` implementation for `BRepFill_MultiLine`
//! (`occt-geom/src/approx/line_tool.rs` implements it for `AppDef_MultiLine`
//! only). Nothing in `BRepFill_Filling` reaches it.

mod build;
mod build_wires;
mod constraints;
mod curve_constraint;

pub use build_wires::{build_wires, find_extremities_of_holes, wire_from_list, OldNewMap};
pub use constraints::{shape_to_plate_order, EdgeFaceAndOrder, FaceAndOrder};
pub use curve_constraint::{from_curve, from_curve_on_surface};

use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_core::kernel::geomabs::Shape;
use occt_geom::curve_on_surface::CurveOnSurface;
use occt_geom::geom_plate::{BuildPlateSurface, CurveConstraint, PointConstraint};

use crate::brep_tool::BRepTool;
use crate::shape::{Face, TopoShape};

/// `class BRepFill_Filling` (`hxx:66-231`).
pub struct Filling {
    /// `myBuilder` (`hxx:210`).
    builder: Option<BuildPlateSurface>,
    /// `myBoundary` (`hxx:211`).
    boundary: Vec<EdgeFaceAndOrder>,
    /// `myConstraints` (`hxx:212`).
    constraints: Vec<EdgeFaceAndOrder>,
    /// `myFreeConstraints` (`hxx:213`).
    free_constraints: Vec<FaceAndOrder>,
    /// `myPoints` (`hxx:214`).
    points: Vec<PointConstraint>,
    /// `myOldNewMap` (`hxx:215`).
    old_new_map: OldNewMap,
    /// `myGenerated` (`hxx:217`).
    generated: Vec<TopoShape>,
    /// `myFace` (`hxx:218`).
    face: Option<Face>,
    /// `myInitFace` (`hxx:219`).
    init_face: Option<Face>,
    /// `myTol2d` (`hxx:220`).
    tol2d: f64,
    /// `myTol3d` (`hxx:221`).
    tol3d: f64,
    /// `myTolAng` (`hxx:222`).
    tol_ang: f64,
    /// `myTolCurv` (`hxx:223`).
    tol_curv: f64,
    /// `myMaxDeg` (`hxx:224`).
    max_deg: i32,
    /// `myMaxSegments` (`hxx:225`).
    max_segments: i32,
    /// `myDegree` (`hxx:226`).
    degree: i32,
    /// `myNbPtsOnCur` (`hxx:227`).
    nb_pts_on_cur: i32,
    /// `myNbIter` (`hxx:228`).
    nb_iter: i32,
    /// `myAnisotropie` (`hxx:229`).
    anisotropie: bool,
    /// `myIsInitFaceGiven` (`hxx:230`).
    is_init_face_given: bool,
    /// `myIsDone` (`hxx:231`).
    is_done: bool,
}

impl Default for Filling {
    /// The header defaults of `hxx:73-82`.
    fn default() -> Self {
        Self::new(3, 15, 2, false, 0.00001, 0.0001, 0.01, 0.1, 8, 9)
    }
}

impl Filling {
    /// `BRepFill_Filling(Degree, NbPtsOnCur, NbIter, Anisotropie, Tol2d, Tol3d,
    /// TolAng, TolCurv, MaxDeg, MaxSegments)` (`cxx:142-171`).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        degree: i32,
        nb_pts_on_cur: i32,
        nb_iter: i32,
        anisotropie: bool,
        tol2d: f64,
        tol3d: f64,
        tol_ang: f64,
        tol_curv: f64,
        max_deg: i32,
        max_segments: i32,
    ) -> Self {
        Self {
            builder: None,
            boundary: Vec::new(),
            constraints: Vec::new(),
            free_constraints: Vec::new(),
            points: Vec::new(),
            old_new_map: OldNewMap::new(),
            generated: Vec::new(),
            face: None,
            init_face: None,
            tol2d,
            tol3d,
            tol_ang,
            tol_curv,
            max_deg,
            max_segments,
            degree,
            nb_pts_on_cur,
            nb_iter,
            anisotropie,
            is_init_face_given: false,
            is_done: false,
        }
    }

    /// `SetConstrParam(Tol2d, Tol3d, TolAng, TolCurv)` (`cxx:177-186`).
    pub fn set_constr_param(&mut self, tol2d: f64, tol3d: f64, tol_ang: f64, tol_curv: f64) {
        self.tol2d = tol2d;
        self.tol3d = tol3d;
        self.tol_ang = tol_ang;
        self.tol_curv = tol_curv;
    }

    /// `SetResolParam(Degree, NbPtsOnCur, NbIter, Anisotropie)` (`cxx:188-197`).
    pub fn set_resol_param(
        &mut self,
        degree: i32,
        nb_pts_on_cur: i32,
        nb_iter: i32,
        anisotropie: bool,
    ) {
        self.degree = degree;
        self.nb_pts_on_cur = nb_pts_on_cur;
        self.nb_iter = nb_iter;
        self.anisotropie = anisotropie;
    }

    /// `SetApproxParam(MaxDeg, MaxSegments)` (`cxx:201-205`).
    pub fn set_approx_param(&mut self, max_deg: i32, max_segments: i32) {
        self.max_deg = max_deg;
        self.max_segments = max_segments;
    }

    /// `LoadInitSurface(aFace)` (`cxx:209-213`).
    pub fn load_init_surface(&mut self, a_face: &Face) {
        self.init_face = Some(a_face.clone());
        self.is_init_face_given = true;
    }

    /// `Add(anEdge, Order, IsBound)` (`cxx:217-233`): a constraint that also
    /// defines an edge of the face's wire, with no support face.
    pub fn add(&mut self, an_edge: &crate::shape::Edge, order: Shape, is_bound: bool) -> i32 {
        // `TopoDS_Face NullFace;` (`cxx:219`).
        let edge_face_and_order = EdgeFaceAndOrder::with_face(an_edge, None, order);
        if is_bound {
            self.boundary.push(edge_face_and_order);
            self.boundary.len() as i32
        } else {
            self.constraints.push(edge_face_and_order);
            (self.boundary.len() + self.free_constraints.len() + self.constraints.len()) as i32
        }
    }

    /// `Add(anEdge, Support, Order, IsBound)` (`cxx:237-253`): an edge
    /// constraint whose tangency / curvature is measured against `support`.
    pub fn add_with_support(
        &mut self,
        an_edge: &crate::shape::Edge,
        support: &Face,
        order: Shape,
        is_bound: bool,
    ) -> i32 {
        let edge_face_and_order = EdgeFaceAndOrder::with_face(an_edge, Some(support), order);
        if is_bound {
            self.boundary.push(edge_face_and_order);
            self.boundary.len() as i32
        } else {
            self.constraints.push(edge_face_and_order);
            (self.boundary.len() + self.free_constraints.len() + self.constraints.len()) as i32
        }
    }

    /// `Add(Support, Order)` (`cxx:257-263`): a "free constraint", a face with
    /// no edge of its own.
    pub fn add_free(&mut self, support: &Face, order: Shape) -> i32 {
        self.free_constraints.push(FaceAndOrder::with_face(support, order));
        (self.boundary.len() + self.free_constraints.len()) as i32
    }

    /// `Add(Point)` (`cxx:267-275`).
    pub fn add_point(&mut self, point: &GpPnt) -> Result<i32, String> {
        // `new GeomPlate_PointConstraint(Point, GeomAbs_C0, myTol3d)`.
        let a_pc = PointConstraint::from_point(*point, 0, self.tol3d)?;
        self.points.push(a_pc);
        Ok((self.boundary.len()
            + self.free_constraints.len()
            + self.constraints.len()
            + self.points.len()) as i32)
    }

    /// `Add(U, V, Support, Order)` (`cxx:281-298`).
    pub fn add_point_on_face(
        &mut self,
        u: f64,
        v: f64,
        support: &Face,
        order: Shape,
    ) -> Result<i32, String> {
        let surf = BRepTool::face_surface_world(support).ok_or_else(|| {
            "BRepFill_Filling::Add: the support face has no surface".to_string()
        })?;
        let a_pc = PointConstraint::from_surface(
            u,
            v,
            surf,
            shape_to_plate_order(order),
            self.tol3d,
            self.tol_ang,
            self.tol_curv,
        )?;
        self.points.push(a_pc);
        Ok((self.boundary.len()
            + self.free_constraints.len()
            + self.constraints.len()
            + self.points.len()) as i32)
    }

    /// `void AddConstraints(SeqOfConstraints)` (`cxx:304-389`): turn every
    /// record into a plate curve constraint and load it into the builder.
    ///
    /// `Err` is one of the two reference throws — `"Add"` for an edge with no
    /// curve-on-surface representation (`cxx:336-340`), and the
    /// `"Attempt to access to null object"` that `BRepAdaptor_Curve2d` raises
    /// on an edge without a pcurve on `CurFace` (`cxx:362-365`).
    fn add_constraints(
        builder: &mut BuildPlateSurface,
        seq: &[EdgeFaceAndOrder],
        nb_pts_on_cur: i32,
        tol3d: f64,
        tol_ang: f64,
        tol_curv: f64,
        init_face: Option<&Face>,
    ) -> Result<(), String> {
        for record in seq {
            let cur_edge = &record.edge;
            let cur_order = record.order;

            let constr: CurveConstraint = if let Some(cur_face) = record.face.as_ref() {
                // `cxx:356-378`.
                let surface = BRepTool::face_surface_world(cur_face).ok_or_else(|| {
                    "BRepFill_Filling::AddConstraints: the support face has no surface".to_string()
                })?;
                let (c2d, _f, _l) = crate::boptools_2d::curve_on_surface_range(cur_edge, cur_face)
                    .ok_or_else(|| "Attempt to access to null object".to_string())?;
                let cos = CurveOnSurface::new(c2d, surface);
                from_curve_on_surface(
                    Arc::new(cos),
                    shape_to_plate_order(cur_order),
                    nb_pts_on_cur,
                    tol3d,
                    tol_ang,
                    tol_curv,
                )?
            } else if cur_order == Shape::C0 {
                // `cxx:321-326`: `BRepAdaptor_Curve(E)` straight into the
                // free 3D-curve constructor.
                let h_curve = BRepTool::edge_curve_adaptor_world(cur_edge).ok_or_else(|| {
                    "BRepFill_Filling::AddConstraints: the edge has no 3D curve".to_string()
                })?;
                from_curve(h_curve, 0, nb_pts_on_cur, tol3d)?
            } else {
                // `cxx:327-353`: "Pas de representation Topologique / On
                // prend une representation Geometrique".
                let reps = crate::tgeometry::GeometryRegistry::global()
                    .edge_pcurve_reps(&cur_edge.0);
                let Some((surface, c2d, first, last)) = reps.into_iter().next() else {
                    return Err("Add".to_string());
                };
                // `cxx:341-343`: `Surface = Surface->Copy();
                // Surface->Transform(loc.Transformation());`. The port's
                // representation stores the already-located surface, so the
                // copy-and-transform is the identity here.
                let _ = (first, last);
                let cos = CurveOnSurface::new(c2d, surface);
                from_curve_on_surface(
                    Arc::new(cos),
                    shape_to_plate_order(cur_order),
                    nb_pts_on_cur,
                    tol3d,
                    tol_ang,
                    tol_curv,
                )?
            };

            let mut constr = constr;
            if let Some(init_face) = init_face {
                // `cxx:375-387`.
                if let Some((curve2d, first_par, last_par)) =
                    crate::boptools_2d::curve_on_surface_range(cur_edge, init_face)
                {
                    let trimmed = occt_geom2d::trimmed::Geom2dTrimmedCurve::new(
                        curve2d, first_par, last_par,
                    );
                    constr.set_curve2d_on_surf(Arc::new(trimmed));
                }
            }
            builder.add_curve_constraint(constr);
        }
        Ok(())
    }

    /// `IsDone()` (`cxx:782-785`).
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// `Face()` (`cxx:788-791`).
    pub fn face(&self) -> Option<&Face> {
        self.face.as_ref()
    }

    /// `Generated(S)` (`cxx:796-806`): the new edge made from old edge `S`, or
    /// an empty list.
    pub fn generated(&mut self, s: &TopoShape) -> &[TopoShape] {
        self.generated.clear();
        if let Some(found) = self.old_new_map.get(&crate::topo_tools::shape_map_key(s)) {
            self.generated.push(found.clone());
        }
        &self.generated
    }

    /// `G0Error()` (`cxx:812-815`).
    ///
    /// Before `Build` the reference dereferences its null `myBuilder`; the port
    /// answers `0.0`.
    pub fn g0_error(&self) -> f64 {
        self.builder.as_ref().map_or(0.0, |b| b.g0_error())
    }

    /// `G1Error()` (`cxx:820-823`), same null-`myBuilder` note.
    pub fn g1_error(&self) -> f64 {
        self.builder.as_ref().map_or(0.0, |b| b.g1_error())
    }

    /// `G2Error()` (`cxx:828-831`), same null-`myBuilder` note.
    pub fn g2_error(&self) -> f64 {
        self.builder.as_ref().map_or(0.0, |b| b.g2_error())
    }

    /// `G0Error(Index)` (`cxx:840-846`).
    pub fn g0_error_of(&mut self, index: usize) -> Result<f64, String> {
        let builder = self
            .builder
            .as_mut()
            .ok_or_else(|| "BRepFill_Filling::G0Error: Build has not run".to_string())?;
        builder.g0_error_of(index)
    }

    /// `G1Error(Index)` (`cxx:851-857`).
    pub fn g1_error_of(&mut self, index: usize) -> Result<f64, String> {
        let builder = self
            .builder
            .as_mut()
            .ok_or_else(|| "BRepFill_Filling::G1Error: Build has not run".to_string())?;
        builder.g1_error_of(index)
    }

    /// `G2Error(Index)` (`cxx:862-868`).
    pub fn g2_error_of(&mut self, index: usize) -> Result<f64, String> {
        let builder = self
            .builder
            .as_mut()
            .ok_or_else(|| "BRepFill_Filling::G2Error: Build has not run".to_string())?;
        builder.g2_error_of(index)
    }

    /// `myBoundary` length (`myBoundary.Length()`).
    pub fn nb_boundary(&self) -> usize {
        self.boundary.len()
    }
}
