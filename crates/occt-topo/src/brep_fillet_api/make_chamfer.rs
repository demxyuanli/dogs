//! `BRepFilletAPI_MakeChamfer` (`BRepFilletAPI_MakeChamfer.cxx:28-371`,
//! `BRepFilletAPI_MakeChamfer.hxx:48-297`).
//!
//! The front end over `ChFi3d_ChBuilder` (`myBuilder`). As in
//! [`super::MakeFillet`] every member forwards, and `Build` is the only one
//! that reaches `ChFi3d_Builder::Compute`.
//!
//! Translation notes:
//! - `myMap` (`cxx:283-288`) is an insertion-ordered `Vec` tested with
//!   [`crate::topo_tools::is_same`],
//! - `GetDist` / `Dists` / `GetDistAngle` (`cxx:57-60`, `:84-90`, `:114-117`)
//!   are `double&` out parameters in OCCT and became returned values in
//!   [`ChFi3dChBuilder`](crate::chfi3d::builder::ch_builder::ChFi3dChBuilder),
//! - `Dists` copies through two locals (`cxx:86-89`) so that a failing
//!   `myBuilder.Dists` leaves the caller's values untouched; the port's
//!   `Result` makes that arm explicit,
//! - `IsSymetric` / `IsTwoDistances` / `IsDistanceAngle` (`cxx:128-169`) all
//!   compare the contour's `ChFiDS_ChamfMethod`, which the port answers as a
//!   `Result` because `ChFi3d_ChBuilder::IsChamfer` raises on an out-of-range
//!   contour.

use std::sync::Arc;

use crate::abs::State;
use crate::brep_fillet_api::LocalOperation;
use crate::chfi3d::builder::ChFi3dChBuilder;
use crate::chfi_ds::{ChamfMethod, ChamfMode, CircSectionArray};
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools::{explorer_shapes, is_same};

/// `BRepFilletAPI_MakeChamfer` (`hxx:48-297`).
pub struct MakeChamfer {
    /// `myBuilder`.
    builder: ChFi3dChBuilder,
    /// `myShape`.
    shape: Option<TopoShape>,
    /// `myDone`.
    done: bool,
    /// `myGenerated`.
    generated: Vec<TopoShape>,
    /// `myMap` (`cxx:283-288`).
    map: Vec<TopoShape>,
}

impl MakeChamfer {
    /// `BRepFilletAPI_MakeChamfer(S)` (`cxx:29-31`): `ChFi3d_ChBuilder(S)`, i.e.
    /// `Ta = 1.0e-2` (`ChFi3d_ChBuilder.hxx:53`).
    pub fn new(s: &TopoShape) -> Self {
        Self {
            builder: ChFi3dChBuilder::new(s, 1.0e-2),
            shape: None,
            done: false,
            generated: Vec::new(),
            map: Vec::new(),
        }
    }

    /// `Add(Dis, E)` (`cxx:42-45`).
    pub fn add_with_dist(&mut self, dis: f64, e: &Edge) -> Result<(), String> {
        self.builder.add_with_dist(dis, e)
    }

    /// `SetDist(Dis, IC, F)` (`cxx:52-55`).
    pub fn set_dist(&mut self, dis: f64, ic: usize, f: &Face) -> Result<(), String> {
        self.builder.set_dist(dis, ic, f)
    }

    /// `GetDist(IC, Dis)` (`cxx:57-60`).
    pub fn get_dist(&self, ic: usize) -> Result<f64, String> {
        self.builder.get_dist(ic)
    }

    /// `Add(Dis1, Dis2, E, F)` (`cxx:64-70`).
    pub fn add_with_two_dists(
        &mut self,
        dis1: f64,
        dis2: f64,
        e: &Edge,
        f: &Face,
    ) -> Result<(), String> {
        self.builder.add_with_two_dists(dis1, dis2, e, f)
    }

    /// `SetDists(Dis1, Dis2, IC, F)` (`cxx:73-80`).
    pub fn set_dists(
        &mut self,
        dis1: f64,
        dis2: f64,
        ic: usize,
        f: &Face,
    ) -> Result<(), String> {
        self.builder.set_dists(dis1, dis2, ic, f)
    }

    /// `Dists(IC, Dis1, Dis2)` (`cxx:84-90`).
    pub fn dists(&self, ic: usize) -> Result<(f64, f64), String> {
        self.builder.dists(ic)
    }

    /// `AddDA(Dis, Angle, E, F)` (`cxx:94-100`).
    pub fn add_with_dist_angle(
        &mut self,
        dis: f64,
        angle: f64,
        e: &Edge,
        f: &Face,
    ) -> Result<(), String> {
        self.builder.add_with_dist_angle(dis, angle, e, f)
    }

    /// `SetDistAngle(Dis, Angle, IC, F)` (`cxx:104-110`).
    pub fn set_dist_angle(
        &mut self,
        dis: f64,
        angle: f64,
        ic: usize,
        f: &Face,
    ) -> Result<(), String> {
        self.builder.set_dist_angle(dis, angle, ic, f)
    }

    /// `GetDistAngle(IC, Dis, Angle)` (`cxx:114-117`).
    pub fn get_dist_angle(&self, ic: usize) -> Result<(f64, f64), String> {
        self.builder.get_dist_angle(ic)
    }

    /// `SetMode(theMode)` (`cxx:121-124`).
    pub fn set_mode(&mut self, mode: ChamfMode) {
        self.builder.set_mode(mode);
    }

    /// `IsSymetric(IC)` (`cxx:128-138`).
    pub fn is_symetric(&self, ic: usize) -> Result<bool, String> {
        Ok(self.builder.is_chamfer(ic)? == ChamfMethod::Sym)
    }

    /// `IsTwoDistances(IC)` (`cxx:143-154`).
    pub fn is_two_distances(&self, ic: usize) -> Result<bool, String> {
        Ok(self.builder.is_chamfer(ic)? == ChamfMethod::TwoDist)
    }

    /// `IsDistanceAngle(IC)` (`cxx:159-169`).
    pub fn is_distance_angle(&self, ic: usize) -> Result<bool, String> {
        Ok(self.builder.is_chamfer(ic)? == ChamfMethod::DistAngle)
    }

    /// `Builder()` (`cxx:267-270`).
    pub fn builder(&self) -> &crate::topopebrep::build::HBuilder {
        self.builder.base().builder()
    }

    /// `Build(theRange)` (`cxx:274-290`): `myBuilder.Compute()`, then `Done()`
    /// + `myShape = myBuilder.Shape()` when the computation succeeded, and the
    /// result's faces into `myMap`.
    ///
    /// `Err` carries the failure `Compute` reports; OCCT leaves `myDone` false
    /// in that case without rethrowing.
    pub fn build(&mut self) -> Result<(), String> {
        self.builder.compute()?;
        if self.builder.base().is_done() {
            self.done = true;
            self.shape = self.builder.base().shape()?;
            self.map.clear();
            if let Some(shape) = self.shape.clone() {
                self.map = explorer_shapes(&shape, crate::abs::ShapeType::Face);
            }
        }
        Ok(())
    }

    /// `BRepBuilderAPI_MakeShape::IsDone()` (`BRepBuilderAPI_MakeShape.hxx:47`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `BRepBuilderAPI_MakeShape::Shape()` (`BRepBuilderAPI_MakeShape.hxx:50`).
    pub fn shape(&self) -> Option<&TopoShape> {
        self.shape.as_ref()
    }

    /// `Generated(EorV)` (`cxx:305-308`).
    pub fn generated(&mut self, e_or_v: Option<&TopoShape>) -> &[TopoShape] {
        self.builder.base_mut().generated(e_or_v)
    }

    /// `Modified(F)` (`cxx:311-339`): the builder's `Splits(F, OUT)` then
    /// `Splits(F, IN)` then `Splits(F, ON)`, concatenated after clearing
    /// `myGenerated`.
    pub fn modified(&mut self, f: &TopoShape) -> &[TopoShape] {
        self.generated.clear();
        for state in [State::Out, State::In, State::On] {
            self.generated
                .extend(self.builder.base().builder().splits(f, state).iter().cloned());
        }
        &self.generated
    }

    /// `IsDeleted(F)` (`cxx:343-347`): `F` is neither in `myMap` nor split in
    /// any of the three states.
    pub fn is_deleted(&self, f: &TopoShape) -> bool {
        let builder = self.builder.base().builder();
        let in_map = self.map.iter().any(|s| is_same(s, f));
        !(in_map
            || builder.is_split(f, State::Out)
            || builder.is_split(f, State::In)
            || builder.is_split(f, State::On))
    }

    /// The wrapped builder, for the members `BRepFilletAPI_MakeChamfer` does not
    /// forward (the fault accessors belong to `BRepFilletAPI_MakeFillet` only).
    pub fn ch_builder(&self) -> &ChFi3dChBuilder {
        &self.builder
    }

    /// The mutable wrapped builder.
    pub fn ch_builder_mut(&mut self) -> &mut ChFi3dChBuilder {
        &mut self.builder
    }
}

impl LocalOperation for MakeChamfer {
    /// `Add(E)` (`cxx:36-39`).
    fn add(&mut self, e: &Edge) -> Result<(), String> {
        self.builder.add(e)
    }

    /// `ResetContour(IC)` (`cxx:174-177`).
    fn reset_contour(&mut self, ic: usize) -> Result<(), String> {
        self.builder.reset_contour(ic)
    }

    /// `NbContours()` (`cxx:180-183`).
    fn nb_contours(&self) -> usize {
        self.builder.base().nb_elements()
    }

    /// `Contour(E)` (`cxx:187-190`).
    fn contour(&self, e: &Edge) -> usize {
        self.builder.base().contains(e)
    }

    /// `NbEdges(I)` (`cxx:194-198`).
    fn nb_edges(&self, i: usize) -> usize {
        self.builder
            .base()
            .value(i)
            .and_then(|s| s.read().ok().map(|g| g.nb_edges()))
            .unwrap_or(0)
    }

    /// `Edge(I, J)` (`cxx:202-207`).
    fn edge(&self, i: usize, j: usize) -> Option<TopoShape> {
        let spine = self.builder.base().value(i)?;
        let guard = spine.read().ok()?;
        guard.edges(j).cloned()
    }

    /// `Remove(E)` (`cxx:211-214`).
    fn remove(&mut self, e: &Edge) -> Result<(), String> {
        self.builder.base_mut().remove(e);
        Ok(())
    }

    /// `Length(IC)` (`cxx:218-222`).
    fn length(&self, ic: usize) -> f64 {
        self.builder.base().length(ic).unwrap_or(-1.0)
    }

    /// `FirstVertex(IC)` (`cxx:225-228`).
    fn first_vertex(&self, ic: usize) -> Option<TopoShape> {
        self.builder.base().first_vertex(ic)
    }

    /// `LastVertex(IC)` (`cxx:232-235`).
    fn last_vertex(&self, ic: usize) -> Option<TopoShape> {
        self.builder.base().last_vertex(ic)
    }

    /// `Abscissa(IC, V)` (`cxx:240-243`).
    fn abscissa(&self, ic: usize, v: &Vertex) -> f64 {
        self.builder.base().abscissa(ic, v)
    }

    /// `RelativeAbscissa(IC, V)` (`cxx:247-250`).
    fn relative_abscissa(&self, ic: usize, v: &Vertex) -> f64 {
        self.builder.base().relative_abscissa(ic, v)
    }

    /// `ClosedAndTangent(IC)` (`cxx:254-257`).
    fn closed_and_tangent(&self, ic: usize) -> bool {
        self.builder.base().closed_and_tangent(ic)
    }

    /// `Closed(IC)` (`cxx:261-264`).
    fn closed(&self, ic: usize) -> bool {
        self.builder.base().closed(ic)
    }

    /// `Reset()` (`cxx:294-299`): `NotDone()` + `myBuilder.Reset()` +
    /// `myMap.Clear()`.
    fn reset(&mut self) {
        self.done = false;
        self.builder.base_mut().reset();
        self.map.clear();
    }

    /// `Simulate(IC)` (`cxx:352-355`).
    fn simulate(&mut self, ic: usize) -> Result<(), String> {
        self.builder.simulate(ic)
    }

    /// `NbSurf(IC)` (`cxx:359-362`).
    fn nb_surf(&self, ic: usize) -> usize {
        self.builder.nb_surf(ic)
    }

    /// `Sect(IC, IS)` (`cxx:366-371`).
    fn sect(&self, ic: usize, is: usize) -> Option<Arc<CircSectionArray>> {
        self.builder.sect(ic, is)
    }
}
