//! `BRepFilletAPI_MakeFillet` (`BRepFilletAPI_MakeFillet.cxx:30-545`,
//! `BRepFilletAPI_MakeFillet.hxx:52-385`).
//!
//! The class is a front end over `ChFi3d_FilBuilder` (`myBuilder`):
//! `Build` is the only member that does more than forward, and it is the only
//! one that reaches `ChFi3d_Builder::Compute`.
//!
//! Translation notes:
//! - the five `Add` overloads deliberately do **not** call
//!   `ChFi3d_FilBuilder::Add(Radius, E)` (`cxx:66-77` comments the call out):
//!   they `Add(E)`, then `Contains(E, IinC)`, then `SetRadius`. The port keeps
//!   that order, so [`MakeFillet::add_with_radius`] is not
//!   `ChFi3dFilBuilder::add_with_radius` (which forwards to `SetRadius(R, IC, E)`
//!   instead of `SetRadius(R, IC, IinC)` and would miss the `IinC`).
//! - `myMap.Add(ex.Current())` (`cxx:382`) is a `TopTools_ShapeMapHasher` set;
//!   the port keeps an insertion-ordered `Vec` tested with
//!   [`crate::topo_tools::is_same`],
//! - `NbSurfaces()` (`cxx:398-401`) dereferences `Builder()->DataStructure()`,
//!   which is null before `Build`; the port answers `None`,
//! - `NewFaces(I)` (`cxx:405-408`) casts away const on the builder's list; the
//!   port borrows it [read-only](HBuilder::new_faces),
//! - `NotDone()` / `Done()` are the `done` flag; `Reset()` clears it (`cxx:391-396`).

use std::sync::Arc;

use occt_core::gp::GpPnt2d;
use occt_core::gp::GpXY;
use occt_core::precision::CONFUSION;
use occt_core::kernel::geomabs::Shape;
use occt_geom::law::LawFunction;
use occt_geom::Surface as GeomSurface;

use crate::abs::State;
use crate::brep_fillet_api::LocalOperation;
use crate::chfi3d::builder::ChFi3dFilBuilder;
use crate::chfi3d::fillet_shape::FilletShape;
use crate::chfi_ds::{CircSectionArray, ErrorStatus};
use crate::shape::{Edge, TopoShape, Vertex};
use crate::topo_tools::{explorer_shapes, is_same};

/// `BRepFilletAPI_MakeFillet` (`hxx:52-385`).
pub struct MakeFillet {
    /// `myBuilder`.
    builder: ChFi3dFilBuilder,
    /// `myShape`.
    shape: Option<TopoShape>,
    /// `myDone`.
    done: bool,
    /// `myGenerated`.
    generated: Vec<TopoShape>,
    /// `myMap` (`TopTools_DataMapOfShapeInteger` in OCCT 8.0 is a
    /// `TopTools_ShapeMapHasher` set of the result's faces, `cxx:377-385`).
    map: Vec<TopoShape>,
}

impl MakeFillet {
    /// `BRepFilletAPI_MakeFillet(S, FShape)` (`cxx:32-35`):
    /// `ChFi3d_FilBuilder(S, FShape)`, i.e. `Ta = 1.0e-2`
    /// (`ChFi3d_FilBuilder.hxx:53-56`).
    pub fn new(s: &TopoShape, fshape: FilletShape) -> Self {
        Self {
            builder: ChFi3dFilBuilder::new(s, fshape, 1.0e-2),
            shape: None,
            done: false,
            generated: Vec::new(),
            map: Vec::new(),
        }
    }

    /// `SetParams(Tang, Tesp, T2d, TApp3d, TolApp2d, Fleche)` (`cxx:41-47`).
    pub fn set_params(
        &mut self,
        tang: f64,
        tesp: f64,
        t2d: f64,
        tapp3d: f64,
        tolapp2d: f64,
        fleche: f64,
    ) {
        self.builder
            .base_mut()
            .set_params(tang, tesp, t2d, tapp3d, tolapp2d, fleche);
    }

    /// `SetContinuity(InternalContinuity, AngleTol)` (`cxx:52-56`).
    pub fn set_continuity(&mut self, internal_continuity: Shape, angle_tol: f64) {
        self.builder
            .base_mut()
            .set_continuity(internal_continuity, angle_tol);
    }

    /// `Add(Radius, E)` (`cxx:66-77`): `Add(E)`, then a constant radius on the
    /// contour that now holds `E`.
    pub fn add_with_radius(&mut self, radius: f64, e: &Edge) -> Result<(), String> {
        self.builder.add(e)?;
        let (ic, i_in_c) = self.builder.base().contains_with_index(e);
        if ic != 0 {
            self.set_radius(radius, ic, i_in_c)?;
        }
        Ok(())
    }

    /// `Add(R1, R2, E)` (`cxx:80-90`).
    pub fn add_with_two_radii(&mut self, r1: f64, r2: f64, e: &Edge) -> Result<(), String> {
        self.builder.add(e)?;
        let (ic, i_in_c) = self.builder.base().contains_with_index(e);
        if ic != 0 {
            self.set_radius_pair(r1, r2, ic, i_in_c)?;
        }
        Ok(())
    }

    /// `Add(const handle(Law_Function)& L, E)` (`cxx:92-105`).
    pub fn add_with_law(
        &mut self,
        l: Arc<dyn LawFunction>,
        e: &Edge,
    ) -> Result<(), String> {
        self.builder.add(e)?;
        let (ic, i_in_c) = self.builder.base().contains_with_index(e);
        if ic != 0 {
            self.set_radius_law(l, ic, i_in_c)?;
        }
        Ok(())
    }

    /// `Add(const NCollection_Array1<gp_Pnt2d>& UandR, E)` (`cxx:107-117`).
    pub fn add_with_radius_array(&mut self, u_and_r: &[GpPnt2d], e: &Edge) -> Result<(), String> {
        self.builder.add(e)?;
        let (ic, i_in_c) = self.builder.base().contains_with_index(e);
        if ic != 0 {
            self.set_radius_array(u_and_r, ic, i_in_c)?;
        }
        Ok(())
    }

    /// `SetRadius(Radius, IC, IinC)` (`cxx:121-126`): the two end parameters
    /// `(0., Radius)` and `(1., Radius)`.
    pub fn set_radius(&mut self, radius: f64, ic: usize, i_in_c: usize) -> Result<(), String> {
        let i_in_c = i_in_c as i32;
        self.builder
            .set_radius_at(GpXY::new(0.0, radius), ic, i_in_c)?;
        self.builder
            .set_radius_at(GpXY::new(1.0, radius), ic, i_in_c)
    }

    /// `SetRadius(R1, R2, IC, IinC)` (`cxx:130-149`): the two radii are averaged
    /// when `|R1 - R2| < Precision::Confusion()`.
    pub fn set_radius_pair(
        &mut self,
        r1: f64,
        r2: f64,
        ic: usize,
        i_in_c: usize,
    ) -> Result<(), String> {
        let (r1, r2) = if (r1 - r2).abs() < CONFUSION {
            let r = (r1 + r2) * 0.5;
            (r, r)
        } else {
            (r1, r2)
        };
        let i_in_c = i_in_c as i32;
        self.builder.set_radius_at(GpXY::new(0.0, r1), ic, i_in_c)?;
        self.builder.set_radius_at(GpXY::new(1.0, r2), ic, i_in_c)
    }

    /// `SetRadius(const handle(Law_Function)& L, IC, IinC)` (`cxx:152-157`).
    pub fn set_radius_law(
        &mut self,
        l: Arc<dyn LawFunction>,
        ic: usize,
        i_in_c: usize,
    ) -> Result<(), String> {
        self.builder.set_radius_law(l, ic, i_in_c as i32)
    }

    /// `SetRadius(const NCollection_Array1<gp_Pnt2d>& UandR, IC, IinC)`
    /// (`cxx:161-188`): one point is a constant radius, two points are the two
    /// end radii, more points are renormalised by
    /// `(Ucur - Uf) / (Ul - Uf)` and forwarded one by one.
    pub fn set_radius_array(
        &mut self,
        u_and_r: &[GpPnt2d],
        ic: usize,
        i_in_c: usize,
    ) -> Result<(), String> {
        let Some(first) = u_and_r.first() else {
            return Ok(());
        };
        let last = u_and_r[u_and_r.len() - 1];
        match u_and_r.len() {
            1 => self.set_radius(first.y(), ic, i_in_c),
            2 => self.set_radius_pair(first.y(), last.y(), ic, i_in_c),
            _ => {
                let uf = first.x();
                let ul = last.x();
                for point in u_and_r {
                    let ucur = (point.x() - uf) / (ul - uf);
                    self.builder
                        .set_radius_at(GpXY::new(ucur, point.y()), ic, i_in_c as i32)?;
                }
                Ok(())
            }
        }
    }

    /// `IsConstant(IC)` (`cxx:190-193`).
    pub fn is_constant(&self, ic: usize) -> Result<bool, String> {
        self.builder.is_constant(ic)
    }

    /// `Radius(IC)` (`cxx:196-200`).
    pub fn radius(&self, ic: usize) -> Result<f64, String> {
        self.builder.radius(ic)
    }

    /// `IsConstant(IC, E)` (`cxx:211-214`).
    pub fn is_constant_of_edge(&self, ic: usize, e: &Edge) -> Result<bool, String> {
        self.builder.is_constant_of_edge(ic, e)
    }

    /// `Radius(IC, E)` (`cxx:217-221`).
    pub fn radius_of_edge(&self, ic: usize, e: &Edge) -> Result<f64, String> {
        self.builder.radius_of_edge(ic, e)
    }

    /// `SetRadius(Radius, IC, E)` (`cxx:223-228`).
    pub fn set_radius_on_edge(&mut self, radius: f64, ic: usize, e: &Edge) -> Result<(), String> {
        self.builder.set_radius_on_edge(radius, ic, e)
    }

    /// `GetBounds(IC, E, F, L)` (`cxx:230-235`): `None` is OCCT's `false` (no
    /// elementary law for `E`), `Some` is `true` plus the two out parameters.
    pub fn get_bounds(
        &self,
        ic: usize,
        e: &Edge,
    ) -> Result<Option<(f64, f64)>, String> {
        self.builder.get_bounds(ic, e)
    }

    /// `GetLaw(IC, E)` (`cxx:238-241`).
    pub fn get_law(
        &self,
        ic: usize,
        e: &Edge,
    ) -> Result<Option<Arc<dyn LawFunction>>, String> {
        self.builder.get_law(ic, e)
    }

    /// `SetLaw(IC, E, L)` (`cxx:244-251`).
    pub fn set_law(
        &self,
        ic: usize,
        e: &Edge,
        l: Arc<dyn LawFunction>,
    ) -> Result<(), String> {
        self.builder.set_law(ic, e, l)
    }

    /// `SetRadius(Radius, IC, V)` (`cxx:254-258`).
    pub fn set_radius_on_vertex(
        &mut self,
        radius: f64,
        ic: usize,
        v: &Vertex,
    ) -> Result<(), String> {
        self.builder.set_radius_on_vertex(radius, ic, v)
    }

    /// `SetFilletShape(FShape)` (`cxx:261-265`).
    pub fn set_fillet_shape(&mut self, fshape: FilletShape) {
        self.builder.set_fillet_shape(fshape);
    }

    /// `GetFilletShape()` (`cxx:267-272`).
    pub fn get_fillet_shape(&self) -> FilletShape {
        self.builder.get_fillet_shape()
    }

    /// `Build(theRange)` (`cxx:372-386`): `myBuilder.Compute()`, then `Done()`
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

    /// `Builder()` (`cxx:363-368`).
    pub fn builder(&self) -> &crate::topopebrep::build::HBuilder {
        self.builder.base().builder()
    }

    /// `NbSurfaces()` (`cxx:398-401`): `Builder()->DataStructure()->NbSurfaces()`.
    ///
    /// `None` is the null `DataStructure()` handle OCCT dereferences before the
    /// first `Build`.
    pub fn nb_surfaces(&self) -> Option<i32> {
        self.builder
            .base()
            .builder()
            .data_structure()
            .map(|ds| ds.nb_surfaces())
    }

    /// `NewFaces(I)` (`cxx:404-408`).
    pub fn new_faces(&self, i: i32) -> Option<&[TopoShape]> {
        self.builder.base().builder().new_faces(i)
    }

    /// `Generated(EorV)` (`cxx:437-440`): `myBuilder.Generated(EorV)`.
    pub fn generated(&mut self, e_or_v: Option<&TopoShape>) -> &[TopoShape] {
        self.builder.base_mut().generated(e_or_v)
    }

    /// `Modified(F)` (`cxx:443-473`): the builder's `Splits(F, OUT)` then
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

    /// `IsDeleted(F)` (`cxx:477-482`): `F` is neither in `myMap` nor split in
    /// any of the three states.
    pub fn is_deleted(&self, f: &TopoShape) -> bool {
        let builder = self.builder.base().builder();
        let in_map = self.map.iter().any(|s| is_same(s, f));
        !(in_map
            || builder.is_split(f, State::Out)
            || builder.is_split(f, State::In)
            || builder.is_split(f, State::On))
    }

    /// `NbFaultyContours()` (`cxx:485-488`).
    pub fn nb_faulty_contours(&self) -> usize {
        self.builder.base().nb_faulty_contours()
    }

    /// `FaultyContour(I)` (`cxx:491-494`).
    pub fn faulty_contour(&self, i: usize) -> usize {
        self.builder.base().faulty_contour(i)
    }

    /// `NbComputedSurfaces(IC)` (`cxx:497-500`).
    pub fn nb_computed_surfaces(&self, ic: usize) -> usize {
        self.builder.base().nb_computed_surfaces(ic)
    }

    /// `ComputedSurface(IC, IS)` (`cxx:503-510`).
    pub fn computed_surface(&self, ic: usize, is: usize) -> Option<Arc<dyn GeomSurface>> {
        self.builder.base().computed_surface(ic, is)
    }

    /// `NbFaultyVertices()` (`cxx:513-516`).
    pub fn nb_faulty_vertices(&self) -> usize {
        self.builder.base().nb_faulty_vertices()
    }

    /// `FaultyVertex(IV)` (`cxx:519-522`).
    pub fn faulty_vertex(&self, iv: usize) -> Option<Vertex> {
        self.builder.base().faulty_vertex(iv)
    }

    /// `HasResult()` (`cxx:526-529`).
    pub fn has_result(&self) -> bool {
        self.builder.base().has_result()
    }

    /// `BadShape()` (`cxx:533-536`).
    pub fn bad_shape(&self) -> Result<Option<TopoShape>, String> {
        self.builder.base().bad_shape()
    }

    /// `StripeStatus(IC)` (`cxx:539-542`).
    pub fn stripe_status(&self, ic: usize) -> Option<ErrorStatus> {
        self.builder.base().stripe_status(ic)
    }

    /// The wrapped builder, for the members `BRepFilletAPI_MakeFillet` does not
    /// forward.
    pub fn fil_builder(&self) -> &ChFi3dFilBuilder {
        &self.builder
    }

    /// The mutable wrapped builder.
    pub fn fil_builder_mut(&mut self) -> &mut ChFi3dFilBuilder {
        &mut self.builder
    }
}

impl LocalOperation for MakeFillet {
    /// `Add(E)` (`cxx:60-63`).
    fn add(&mut self, e: &Edge) -> Result<(), String> {
        self.builder.add(e)
    }

    /// `ResetContour(IC)` (`cxx:204-207`).
    fn reset_contour(&mut self, ic: usize) -> Result<(), String> {
        self.builder.reset_contour(ic)
    }

    /// `NbContours()` (`cxx:276-279`).
    fn nb_contours(&self) -> usize {
        self.builder.base().nb_elements()
    }

    /// `Contour(E)` (`cxx:283-286`).
    fn contour(&self, e: &Edge) -> usize {
        self.builder.base().contains(e)
    }

    /// `NbEdges(I)` (`cxx:290-295`).
    fn nb_edges(&self, i: usize) -> usize {
        self.builder
            .base()
            .value(i)
            .and_then(|s| s.read().ok().map(|g| g.nb_edges()))
            .unwrap_or(0)
    }

    /// `Edge(I, J)` (`cxx:299-303`).
    fn edge(&self, i: usize, j: usize) -> Option<TopoShape> {
        let spine = self.builder.base().value(i)?;
        let guard = spine.read().ok()?;
        guard.edges(j).cloned()
    }

    /// `Remove(E)` (`cxx:307-310`).
    fn remove(&mut self, e: &Edge) -> Result<(), String> {
        self.builder.base_mut().remove(e);
        Ok(())
    }

    /// `Length(IC)` (`cxx:314-317`).
    fn length(&self, ic: usize) -> f64 {
        self.builder.base().length(ic).unwrap_or(-1.0)
    }

    /// `FirstVertex(IC)` (`cxx:321-324`).
    fn first_vertex(&self, ic: usize) -> Option<TopoShape> {
        self.builder.base().first_vertex(ic)
    }

    /// `LastVertex(IC)` (`cxx:328-331`).
    fn last_vertex(&self, ic: usize) -> Option<TopoShape> {
        self.builder.base().last_vertex(ic)
    }

    /// `Abscissa(IC, V)` (`cxx:335-338`).
    fn abscissa(&self, ic: usize, v: &Vertex) -> f64 {
        self.builder.base().abscissa(ic, v)
    }

    /// `RelativeAbscissa(IC, V)` (`cxx:342-345`).
    fn relative_abscissa(&self, ic: usize, v: &Vertex) -> f64 {
        self.builder.base().relative_abscissa(ic, v)
    }

    /// `ClosedAndTangent(IC)` (`cxx:350-353`).
    fn closed_and_tangent(&self, ic: usize) -> bool {
        self.builder.base().closed_and_tangent(ic)
    }

    /// `Closed(IC)` (`cxx:357-360`).
    fn closed(&self, ic: usize) -> bool {
        self.builder.base().closed(ic)
    }

    /// `Reset()` (`cxx:391-396`): `NotDone()` + `myBuilder.Reset()` +
    /// `myMap.Clear()`.
    fn reset(&mut self) {
        self.done = false;
        self.builder.base_mut().reset();
        self.map.clear();
    }

    /// `Simulate(IC)` (`cxx:412-415`).
    fn simulate(&mut self, ic: usize) -> Result<(), String> {
        self.builder.simulate(ic)
    }

    /// `NbSurf(IC)` (`cxx:419-422`).
    fn nb_surf(&self, ic: usize) -> usize {
        self.builder.nb_surf(ic)
    }

    /// `Sect(IC, IS)` (`cxx:426-431`).
    fn sect(&self, ic: usize, is: usize) -> Option<Arc<CircSectionArray>> {
        self.builder.sect(ic, is)
    }
}
