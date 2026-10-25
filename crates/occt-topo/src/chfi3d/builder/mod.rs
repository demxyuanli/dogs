//! `ChFi3d_Builder` — the root class that computes the surfaces of the fillets
//! and chamfers and rebuilds the shape.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder.hxx:77-848`
//! (OCCT 8.0.0). This module ports the data layer of the class:
//!
//! - the constructor and the five acquisition maps (`ChFi3d_Builder_1.cxx:340-356`),
//! - `SetParams` / `SetContinuity` (`_1.cxx:358-381`),
//! - the result and fault accessors (`_1.cxx:384-556`),
//! - the contour accessors (`_1.cxx:1181-1356`),
//! - `FaceTangency` (`_1.cxx:564-623`),
//! - `UpdateTolesp` (`ChFi3d_Builder_C2.cxx:814-824`),
//! - `Reset` (`ChFi3d_Builder.cxx:924-948`) and `Generated`
//!   (`ChFi3d_Builder.cxx:950-977`),
//! - [`StoreData`](store_data::store_data) (`ChFi3d_Builder_6.cxx:629-1005`),
//!   the four `CompleteData` overloads ([`complete_data`], `_6.cxx:393-624`),
//!   the four `ComputeData` overloads ([`compute_data`], `_6.cxx:1009-1342`, and
//!   the face / face one, `_6.cxx:1676-2507`), the four `SimulData` overloads
//!   ([`simul_data`], `_6.cxx:1346-1669`, and the face / face one,
//!   `_6.cxx:2511-2703`) and the `UpdateLine` helper
//!   ([`update_line::update_line`], `_6.cxx:344-387`),
//! - [`ChFi3dBuilder::search_face`] (`ChFi3d_Builder_2.cxx:1523-1627`), the
//!   `SearchFace` member the face / face `ComputeData` walks its extremities
//!   with,
//! - [`ChFi3dBuilder::perform_extremity`]
//!   ([`perform_extremity`], `ChFi3d_Builder_1.cxx:730-885`) and
//!   [`ChFi3dBuilder::cone_x_faces`] (`ChFi3d_Builder_2.cxx:859-880`), with the
//!   `ReorderFaces` / `MakeOffsetEdge` / `TangentExtremity` / `TangentOnVertex`
//!   file statics of `_1.cxx` in [`perform_helpers`],
//! - [`ChFi3dBuilder::perform_element`] ([`perform_element`],
//!   `ChFi3d_Builder_1.cxx:887-1179`), the downstream / upstream tangency walk
//!   that fills a spine,
//! - [`ChFi3dBuilder::call_perform_surf`] ([`call_perform_surf`],
//!   `ChFi3d_Builder_2.cxx:635-818`) and
//!   [`ChFi3dBuilder::stripe_orientations`] (`_2.cxx:826-857`), the two drivers
//!   that choose the face pair / `Choix` handed to the subclass, with the two
//!   protected pure virtuals they dispatch to stated as
//!   [`ChFi3dSurfOps`] (`ChFi3d_Builder.hxx:203-221`, `:445-466`),
//! - [`ChFi3dBuilder::start_sol_el_spine`] ([`start_sol`],
//!   `ChFi3d_Builder_2.cxx:900-1108`), the first-section search along the
//!   guideline, with its file statics `BonVoisin` (`_2.cxx:213-318`) and
//!   `Projection` (`_2.cxx:320-370`) in [`builder_2_helpers::start_sol`] and
//!   `ChFi3d_InterPlaneEdge` (`ChFi3d_Builder_0.cxx:718-767`) in
//!   [`crate::chfi3d::builder_0_helpers`],
//! - [`ChFi3dBuilder::start_sol_with_ref`] ([`start_sol_obstacle`],
//!   `_2.cxx:1162-1520`), the overload that walks the last section of an
//!   already computed band onto the next face, with `ChFi3d_BuildPlane`
//!   (`_2.cxx:1113-1144`) and `isTangentToArc` (`_2.cxx:192-209`) in
//!   [`builder_2_helpers::build_plane`] / [`builder_2_helpers::is_tangent_to_arc`]
//!   and `BRep_Tool::Parameters` (`BRep_Tool.cxx:1661-1703`) in
//!   [`crate::brep_tool::BRepTool::parameters_of_vertex_on_face`],
//! - the four `SurfData`-sequence surgery statics of `_2.cxx:2047-2204`
//!   (`ChFi3d_Purge`, `InsertAfter`, `RemoveSD`, `InsertBefore`) in
//!   [`stripe_surgery`],
//! - [`ChFi3dBuilder::perform_set_of_surf_on_el_spine`]
//!   ([`perform_set_of_surf`], `_2.cxx:2207-3003`), the section-by-section walk
//!   of one guideline: the driver and its `Walk` state in
//!   `perform_set_of_surf/mod.rs`, the starting solution in `start.rs`
//!   (`:2283-2412`), one turn of the `while (!fini)` loop in `step.rs`
//!   (`:2416-2988`), the three obstacle arms in
//!   `perform_set_of_surf/obstacles.rs` and the splice in
//!   `perform_set_of_surf/splice.rs`.
//!
//! UNPORTED — the orchestration above the data layer, whose `BRepBlend` /
//! `BlendFunc` dependencies (the four walking drivers, the four `*Rad` surface
//! functions, `AppFuncRoot` / `AppSurface`) are now all ported:
//!
//! - `Compute` (`ChFi3d_Builder.cxx:178-680`), `PerformSingularCorner`
//!   (`:682-757`), `PerformFilletOnVertex` (`:759-922`), `ExtentAnalyse`
//!   (`:144-176`),
//! - `PerformFirstSection` (`ChFi3d_Builder.hxx:577-588`) is declared in
//!   [`ChFi3dSurfOps`] but has no implementor (pure virtual in OCCT, implemented
//!   by `ChFi3d_FilBuilder` / `ChFi3d_ChBuilder`), together with `SimulKPart`
//!   (`hxx:201`), `PerformTwoCornerbyInter`, `SplitKPart`,
//! - `PerformSetOfSurfOnElSpine` (`ChFi3d_Builder_2.cxx:2207-3003`) — ported,
//!   but still without a caller; its obstacle arms need the edge-contact
//!   `PerformSurf` / `SimulSurf` of the subclasses, so they stay on the
//!   `ChFi3d_Builder_NotImp.cxx` `throw` defaults for now —, `PerformOneCorner`,
//!   `IntersectMoreCorner`, `PerformMoreSurfdata`, `PerformIntersectionAtEnd`,
//!   `FindFace` (`ChFi3d_Builder_C1.cxx:4533-4600`), `MoreSurfdata`
//!   (`ChFi3d_Builder_C1.cxx:4601-...`).
//!
//! Divergences that are not heuristics (each is also stated at the call site):
//! - `occ::handle<T>` becomes an owned value (`HDataStructure`, `HBuilder`) or
//!   the shared-mutable [`SpineHandle`] / [`StripeHandle`] the crate uses for
//!   the transient data,
//! - OCCT raises (`Standard_NoSuchObject`, `Standard_ConstructionError`,
//!   `Standard_OutOfRange`) from `Shape` / `BadShape` / `FaceTangency` and from
//!   the accessors that index an unloaded guideline; the port answers `Result`
//!   where OCCT has a message, `Option` where it does not,
//! - a null `TopoDS_Shape` member (`myShapeResult`, `badShape`) has no
//!   counterpart in the port, so those members are `Option<TopoShape>`,
//! - `NCollection_List<T>` / `NCollection_Sequence<T>` become `Vec<T>`, and
//!   `NCollection_DataMap<TopoDS_Shape, T>` becomes a `Vec` of pairs addressed
//!   with `TopoDS_Shape::IsSame`, like [`crate::chfi_ds::ChFiDSMap`].

use std::sync::Arc;

use occt_core::kernel::geomabs::Shape;
use occt_geom::Surface as GeomSurface;

mod call_perform_surf;
mod complete_data;
mod compute_data;
mod compute_data_face_face;
mod perform_element;
mod perform_extremity;
mod perform_helpers;
mod perform_set_of_surf;
mod search_face;
mod simul_data;
mod simul_data_face_face;
mod start_sol;
mod start_sol_obstacle;
mod store_data;
mod stripe_surgery;
mod surf_ops;
mod update_line;

pub use store_data::store_data;
pub use surf_ops::ChFi3dSurfOps;
pub use update_line::update_line;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::chfi_ds::{ChFiDSMap, ErrorStatus, Regul, SpineHandle, StripeHandle, StripeMap};
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools::is_same;
use crate::topopebrep::build::HBuilder;
use crate::topopebrep::ds::{BuildTool, HDataStructure};
use crate::topopebrep::tool::out_curve_type::OutCurveType;

use super::tangent_faces::is_tangent_faces;

/// The `mkbuildtool()` helper of `ChFi3d_Builder_1.cxx:330-338`:
/// `TopOpeBRepTool_GeomTool(TopOpeBRepTool_BSPLINE1, true, false, false)` into
/// a `TopOpeBRepDS_BuildTool`, then `OverWrite(false)` / `Translate(false)`.
fn mk_build_tool() -> BuildTool {
    let mut bt = BuildTool::new();
    bt.change_geom_tool()
        .define(OutCurveType::BSpline1, true, false, false);
    bt.set_over_write(false);
    bt.set_translate(false);
    bt
}

/// `Spine()` of a stripe handle, with the `IsNull()` of the OCCT handle
/// answered as `None`.
fn spine_of(stripe: &StripeHandle) -> Option<SpineHandle> {
    let guard = stripe.read().ok()?;
    guard.spine().cloned()
}

/// The `j` of `sp->Edges(j)` that is `E`, `1`-based
/// (`ChFi3d_Builder_1.cxx:1242-1248`, `:1263-1272`).
fn edge_index_in_spine(spine: &SpineHandle, e: &TopoShape) -> Option<usize> {
    let guard = spine.read().ok()?;
    for j in 1..=guard.nb_edges() {
        if let Some(candidate) = guard.edges(j) {
            if is_same(candidate, e) {
                return Some(j);
            }
        }
    }
    None
}

/// Source: `ChFi3d_Builder.hxx:77-848`.
pub struct ChFi3dBuilder {
    // Data members, in the OCCT declaration order (`hxx:742-766`, `:842-846`).
    /// `tolappangle`
    tolappangle: f64,
    /// `tolesp`
    tolesp: f64,
    /// `tol2d`
    tol2d: f64,
    /// `tolapp3d`
    tolapp3d: f64,
    /// `tolapp2d`
    tolapp2d: f64,
    /// `fleche`
    fleche: f64,
    /// `myConti`
    my_conti: Shape,
    /// `myEFMap`
    my_ef_map: ChFiDSMap,
    /// `myESoMap`
    my_eso_map: ChFiDSMap,
    /// `myEShMap`
    my_esh_map: ChFiDSMap,
    /// `myVFMap`
    my_vf_map: ChFiDSMap,
    /// `myVEMap`
    my_ve_map: ChFiDSMap,
    /// `myDS`
    my_ds: HDataStructure,
    /// `myCoup`
    my_coup: HBuilder,
    /// `myListStripe`
    my_list_stripe: Vec<StripeHandle>,
    /// `myVDataMap`
    my_v_data_map: StripeMap,
    /// `myRegul`
    my_regul: Vec<Regul>,
    /// `badstripes`
    bad_stripes: Vec<StripeHandle>,
    /// `badvertices`
    bad_vertices: Vec<TopoShape>,
    /// `myEVIMap` — a `NCollection_DataMap<TopoDS_Shape, List<int>>`; the port
    /// keeps the insertion order in a `Vec` and looks keys up with `IsSame`.
    my_evi_map: Vec<(TopoShape, Vec<i32>)>,
    /// `myEdgeFirstFace`
    my_edge_first_face: Vec<(TopoShape, TopoShape)>,
    /// `done`
    done: bool,
    /// `hasresult`
    has_result: bool,
    /// `myShape`
    my_shape: TopoShape,
    /// `angular`
    angular: f64,
    /// `myGenerated`
    my_generated: Vec<TopoShape>,
    /// `myShapeResult`
    my_shape_result: Option<TopoShape>,
    /// `badShape`
    bad_shape: Option<TopoShape>,
}

impl ChFi3dBuilder {
    /// `ChFi3d_Builder(S, Ta)` (`ChFi3d_Builder_1.cxx:340-356`).
    ///
    /// The five maps are filled with `ChFiDS_Map::Fill`; the parameters are the
    /// OCCT defaults `(Ta, 1.0e-4, 1.e-5, 1.e-4, 1.e-5, 1.e-3)` and the internal
    /// continuity is `GeomAbs_C1` with `Ta` as angular tolerance.
    pub fn new(shape: &TopoShape, ta: f64) -> Self {
        let mut my_ef_map = ChFiDSMap::new();
        my_ef_map.fill(shape, ShapeType::Edge, ShapeType::Face);
        let mut my_eso_map = ChFiDSMap::new();
        my_eso_map.fill(shape, ShapeType::Edge, ShapeType::Solid);
        let mut my_esh_map = ChFiDSMap::new();
        my_esh_map.fill(shape, ShapeType::Edge, ShapeType::Shell);
        let mut my_vf_map = ChFiDSMap::new();
        my_vf_map.fill(shape, ShapeType::Vertex, ShapeType::Face);
        let mut my_ve_map = ChFiDSMap::new();
        my_ve_map.fill(shape, ShapeType::Vertex, ShapeType::Edge);
        let mut builder = Self {
            tolappangle: 0.0,
            tolesp: 0.0,
            tol2d: 0.0,
            tolapp3d: 0.0,
            tolapp2d: 0.0,
            fleche: 0.0,
            my_conti: Shape::C0,
            my_ef_map,
            my_eso_map,
            my_esh_map,
            my_vf_map,
            my_ve_map,
            my_ds: HDataStructure::new(),
            my_coup: HBuilder::new(mk_build_tool()),
            my_list_stripe: Vec::new(),
            my_v_data_map: StripeMap::new(),
            my_regul: Vec::new(),
            bad_stripes: Vec::new(),
            bad_vertices: Vec::new(),
            my_evi_map: Vec::new(),
            my_edge_first_face: Vec::new(),
            done: false,
            has_result: false,
            my_shape: shape.clone(),
            angular: 0.0,
            my_generated: Vec::new(),
            my_shape_result: None,
            bad_shape: None,
        };
        builder.set_params(ta, 1.0e-4, 1.0e-5, 1.0e-4, 1.0e-5, 1.0e-3);
        builder.set_continuity(Shape::C1, ta);
        builder
    }

    /// `SetParams(Tang, Tesp, T2d, TApp3d, TolApp2d, Fleche)`
    /// (`ChFi3d_Builder_1.cxx:358-372`).
    pub fn set_params(
        &mut self,
        tang: f64,
        tesp: f64,
        t2d: f64,
        tapp3d: f64,
        tolapp2d: f64,
        fleche: f64,
    ) {
        self.angular = tang;
        self.tolesp = tesp;
        self.tol2d = t2d;
        self.tolapp3d = tapp3d;
        self.tolapp2d = tolapp2d;
        self.fleche = fleche;
    }

    /// `SetContinuity(InternalContinuity, AngularTolerance)`
    /// (`ChFi3d_Builder_1.cxx:376-381`).
    pub fn set_continuity(&mut self, internal_continuity: Shape, angular_tolerance: f64) {
        self.my_conti = internal_continuity;
        self.tolappangle = angular_tolerance;
    }

    /// `IsDone()` (`ChFi3d_Builder_1.cxx:384-387`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `Shape()` (`ChFi3d_Builder_1.cxx:391-395`): `Standard_NoSuchObject` when
    /// the computation did not run, which the port answers as `Err`; the
    /// returned `Option` is the OCCT null-shape case.
    pub fn shape(&self) -> Result<Option<TopoShape>, String> {
        if !self.done {
            return Err("ChFi3d_Builder::Shape() - no result".to_string());
        }
        Ok(self.my_shape_result.clone())
    }

    /// `NbFaultyContours()` (`ChFi3d_Builder_1.cxx:399-402`).
    pub fn nb_faulty_contours(&self) -> usize {
        self.bad_stripes.len()
    }

    /// `FaultyContour(I)` (`ChFi3d_Builder_1.cxx:405-434`): the `I`-th faulty
    /// contour is the position of its stripe in `myListStripe`, `0` when the
    /// stripe is not there (the OCCT sentinel).
    pub fn faulty_contour(&self, i: usize) -> usize {
        let Some(stripe) = i.checked_sub(1).and_then(|k| self.bad_stripes.get(k)) else {
            return 0;
        };
        for (k, candidate) in self.my_list_stripe.iter().enumerate() {
            if Arc::ptr_eq(candidate, stripe) {
                return k + 1;
            }
        }
        0
    }

    /// `NbComputedSurfaces(IC)` (`ChFi3d_Builder_1.cxx:438-467`).
    pub fn nb_computed_surfaces(&self, ic: usize) -> usize {
        let Some(stripe) = self.stripe_at(ic) else {
            return 0;
        };
        let Ok(guard) = stripe.read() else {
            return 0;
        };
        if guard.spine().is_none() {
            return 0;
        }
        guard.set_of_surf_data().len()
    }

    /// `ComputedSurface(IC, IS)` (`ChFi3d_Builder_1.cxx:470-488`): the surface
    /// of the `IS`-th `SurfData` of the `IC`-th contour, read through the DS.
    pub fn computed_surface(&self, ic: usize, is: usize) -> Option<Arc<dyn GeomSurface>> {
        let stripe = self.stripe_at(ic)?;
        let guard = stripe.read().ok()?;
        let data = guard.set_of_surf_data().get(is.checked_sub(1)?)?;
        self.my_ds.surface(data.surf()).surface().cloned()
    }

    /// `NbFaultyVertices()` (`ChFi3d_Builder_1.cxx:491-494`).
    pub fn nb_faulty_vertices(&self) -> usize {
        self.bad_vertices.len()
    }

    /// `FaultyVertex(IV)` (`ChFi3d_Builder_1.cxx:497-513`), `IV` `1`-based.
    pub fn faulty_vertex(&self, iv: usize) -> Option<Vertex> {
        let shape = iv.checked_sub(1).and_then(|k| self.bad_vertices.get(k))?;
        Vertex::wrap(shape.clone())
    }

    /// `HasResult()` (`ChFi3d_Builder_1.cxx:517-520`).
    pub fn has_result(&self) -> bool {
        self.has_result
    }

    /// `BadShape()` (`ChFi3d_Builder_1.cxx:524-528`): `Standard_NoSuchObject`
    /// when no partial result was stored, which the port answers as `Err`.
    pub fn bad_shape(&self) -> Result<Option<TopoShape>, String> {
        if !self.has_result {
            return Err("ChFi3d_Builder::BadShape() - no result".to_string());
        }
        Ok(self.bad_shape.clone())
    }

    /// `StripeStatus(IC)` (`ChFi3d_Builder_1.cxx:532-548`).
    pub fn stripe_status(&self, ic: usize) -> Option<ErrorStatus> {
        let stripe = self.stripe_at(ic)?;
        let spine = spine_of(&stripe)?;
        let guard = spine.read().ok()?;
        Some(guard.error_status())
    }

    /// `Remove(E)` (`ChFi3d_Builder_1.cxx:1181-1197`): drops the first stripe
    /// one of whose edges is `E`.
    pub fn remove(&mut self, e: &Edge) {
        let mut found = None;
        for (k, stripe) in self.my_list_stripe.iter().enumerate() {
            let Some(spine) = spine_of(stripe) else {
                continue;
            };
            if edge_index_in_spine(&spine, &e.0).is_some() {
                found = Some(k);
                break;
            }
        }
        if let Some(k) = found {
            self.my_list_stripe.remove(k);
        }
    }

    /// `Value(I)` (`ChFi3d_Builder_1.cxx:1201-1210`).
    pub fn value(&self, i: usize) -> Option<SpineHandle> {
        let stripe = self.stripe_at(i)?;
        spine_of(&stripe)
    }

    /// `NbElements()` (`ChFi3d_Builder_1.cxx:1213-1226`): the stripes up to the
    /// first one whose `Spine()` is null.
    pub fn nb_elements(&self) -> usize {
        let mut n = 0;
        for stripe in &self.my_list_stripe {
            if spine_of(stripe).is_none() {
                break;
            }
            n += 1;
        }
        n
    }

    /// `Contains(E)` (`ChFi3d_Builder_1.cxx:1231-1252`): the `1`-based index of
    /// the contour holding `E`, `0` when no contour does.
    pub fn contains(&self, e: &Edge) -> usize {
        self.contains_with_index(e).0
    }

    /// `Contains(E, IndexInSpine)` (`ChFi3d_Builder_1.cxx:1255-1277`): the
    /// contour index and the `1`-based index of `E` inside that contour.
    pub fn contains_with_index(&self, e: &Edge) -> (usize, usize) {
        for (k, stripe) in self.my_list_stripe.iter().enumerate() {
            let Some(spine) = spine_of(stripe) else {
                break;
            };
            if let Some(j) = edge_index_in_spine(&spine, &e.0) {
                return (k + 1, j);
            }
        }
        (0, 0)
    }

    /// `Length(IC)` (`ChFi3d_Builder_1.cxx:1281-1290`):
    /// `Value(IC)->LastParameter(NbEdges())`. The `-1` of an out-of-range `IC`
    /// and the `Standard_OutOfRange` of a not-yet-`Load`ed guideline both
    /// answer `None`.
    pub fn length(&self, ic: usize) -> Option<f64> {
        if ic == 0 || ic > self.nb_elements() {
            return None;
        }
        let spine = self.value(ic)?;
        let guard = spine.read().ok()?;
        guard.last_parameter_of(guard.nb_edges() as i32)
    }

    /// `FirstVertex(IC)` (`ChFi3d_Builder_1.cxx:1293-1300`).
    pub fn first_vertex(&self, ic: usize) -> Option<TopoShape> {
        if ic == 0 || ic > self.nb_elements() {
            return None;
        }
        let spine = self.value(ic)?;
        let guard = spine.read().ok()?;
        guard.first_vertex()
    }

    /// `LastVertex(IC)` (`ChFi3d_Builder_1.cxx:1303-1311`).
    pub fn last_vertex(&self, ic: usize) -> Option<TopoShape> {
        if ic == 0 || ic > self.nb_elements() {
            return None;
        }
        let spine = self.value(ic)?;
        let guard = spine.read().ok()?;
        guard.last_vertex()
    }

    /// `Abscissa(IC, V)` (`ChFi3d_Builder_1.cxx:1315-1322`): `-1.` both when
    /// `IC` is out of range and when `V` is not on the contour, as in OCCT.
    pub fn abscissa(&self, ic: usize, v: &Vertex) -> f64 {
        if ic == 0 || ic > self.nb_elements() {
            return -1.0;
        }
        let Some(spine) = self.value(ic) else {
            return -1.0;
        };
        let Ok(guard) = spine.read() else {
            return -1.0;
        };
        guard.absc_of_vertex(&v.0)
    }

    /// `RelativeAbscissa(IC, V)` (`ChFi3d_Builder_1.cxx:1325-1334`).
    pub fn relative_abscissa(&self, ic: usize, v: &Vertex) -> f64 {
        if ic == 0 || ic > self.nb_elements() {
            return -1.0;
        }
        self.abscissa(ic, v) / self.length(ic).unwrap_or(-1.0)
    }

    /// `Closed(IC)` (`ChFi3d_Builder_1.cxx:1337-1345`).
    pub fn closed(&self, ic: usize) -> bool {
        if ic == 0 || ic > self.nb_elements() {
            return false;
        }
        let Some(spine) = self.value(ic) else {
            return false;
        };
        let Ok(guard) = spine.read() else {
            return false;
        };
        // `IsClosed()` raises on a guideline without edges; the port answers
        // `None`, which is "not closed" here.
        guard.is_closed().unwrap_or(false)
    }

    /// `ClosedAndTangent(IC)` (`ChFi3d_Builder_1.cxx:1348-1355`): OCCT asks
    /// `IsPeriodic()`, i.e. `firstState == ChFiDS_Closed`.
    pub fn closed_and_tangent(&self, ic: usize) -> bool {
        if ic == 0 || ic > self.nb_elements() {
            return false;
        }
        let Some(spine) = self.value(ic) else {
            return false;
        };
        let Ok(guard) = spine.read() else {
            return false;
        };
        guard.is_periodic()
    }

    /// `Builder()` (`ChFi3d_Builder_1.cxx:552-555`).
    pub fn builder(&self) -> &HBuilder {
        &self.my_coup
    }

    /// The protected `myEFMap` (`hxx:750`), read by the subclasses
    /// (`ChFi3d_ChBuilder.cxx:205`, `ChFi3d_FilBuilder_C2.cxx:434`).
    pub fn ef_map(&self) -> &ChFiDSMap {
        &self.my_ef_map
    }

    /// The protected `myESoMap` (`hxx:751`).
    pub fn eso_map(&self) -> &ChFiDSMap {
        &self.my_eso_map
    }

    /// The protected `myEShMap` (`hxx:752`).
    pub fn esh_map(&self) -> &ChFiDSMap {
        &self.my_esh_map
    }

    /// The protected `myVFMap` (`hxx:753`).
    pub fn vf_map(&self) -> &ChFiDSMap {
        &self.my_vf_map
    }

    /// The protected `myShape` (`hxx:842`), the input shape of the builder.
    pub fn input_shape(&self) -> &TopoShape {
        &self.my_shape
    }

    /// `FaceTangency(E0, E1, V)` (`ChFi3d_Builder_1.cxx:564-623`): the faces
    /// opposite to `E0` and `E1` are tangent through the faces starting at `V`.
    ///
    /// `Err` carries the `Standard_ConstructionError("ChFi3d_Builder:only 2
    /// faces")` OCCT raises when an edge has more than two faces, and the
    /// failure `ChFi3d::IsTangentFaces` may raise.
    pub fn face_tangency(&self, e0: &Edge, e1: &Edge, v: &Vertex) -> Result<bool, String> {
        // It is checked if the connection is not on a regular edge.
        let mut faces: Vec<Face> = Vec::with_capacity(2);
        for value in self.my_ef_map.get_from_key(&e1.0) {
            if faces.len() > 1 {
                return Err("ChFi3d_Builder:only 2 faces".to_string());
            }
            let face = Face::wrap(value.clone())
                .ok_or_else(|| "ChFi3d_Builder:only 2 faces".to_string())?;
            faces.push(face);
        }
        if faces.len() < 2 {
            return Ok(false);
        }
        // `ChFi3d::IsTangentFaces(E1, F[0], F[1])` with its `GeomAbs_G1`
        // default (`ChFi3d.hxx:44-48`).
        if is_tangent_faces(e1, &faces[0], &faces[1], Shape::G1)? {
            return Ok(false);
        }

        for value in self.my_ve_map.get_from_key(&v.0) {
            let Some(ec) = Edge::wrap(value.clone()) else {
                continue;
            };
            if is_same(&ec.0, &e0.0) || is_same(&ec.0, &e1.0) {
                continue;
            }
            let orientation = ec.0.orientation();
            if orientation == Orientation::Internal || orientation == Orientation::External {
                continue;
            }
            if BRepTool::is_degenerated(&ec) {
                continue;
            }
            let mut faces: Vec<Face> = Vec::with_capacity(2);
            for candidate in self.my_ef_map.get_from_key(&ec.0) {
                if faces.len() > 1 {
                    return Err("ChFi3d_Builder:only 2 faces".to_string());
                }
                let face = Face::wrap(candidate.clone())
                    .ok_or_else(|| "ChFi3d_Builder:only 2 faces".to_string())?;
                faces.push(face);
            }
            if faces.len() < 2 {
                return Ok(false);
            }
            if !is_tangent_faces(&ec, &faces[0], &faces[1], Shape::G1)? {
                return Ok(false);
            }
        }
        Ok(true)
    }

    /// `Reset()` (`ChFi3d_Builder.cxx:924-948`): back to the state of the last
    /// acquisition, dropping the stripes that lost their spine.
    pub fn reset(&mut self) {
        self.done = false;
        self.my_v_data_map.clear();
        self.my_regul.clear();
        self.my_evi_map.clear();
        self.bad_stripes.clear();
        self.bad_vertices.clear();

        let mut i = 0;
        while i < self.my_list_stripe.len() {
            if spine_of(&self.my_list_stripe[i]).is_some() {
                if let Ok(mut guard) = self.my_list_stripe[i].write() {
                    guard.reset();
                }
                i += 1;
            } else {
                self.my_list_stripe.remove(i);
            }
        }
    }

    /// `Generated(EouV)` (`ChFi3d_Builder.cxx:950-977`): the faces built by the
    /// topological operation for an edge or a vertex, through `myEVIMap`.
    ///
    /// `e_or_v` is an `Option` because a null `TopoDS_Shape` handle has no
    /// counterpart in the port.
    pub fn generated(&mut self, e_or_v: Option<&TopoShape>) -> &[TopoShape] {
        self.my_generated.clear();
        let Some(eou_v) = e_or_v else {
            return &self.my_generated;
        };
        if eou_v.shape_type() != ShapeType::Edge && eou_v.shape_type() != ShapeType::Vertex {
            return &self.my_generated;
        }
        let Some(list) = self.evi_indices(eou_v).map(|l| l.to_vec()) else {
            return &self.my_generated;
        };
        for i in list {
            if let Some(faces) = self.my_coup.new_faces(i) {
                self.my_generated.extend_from_slice(faces);
            }
        }
        &self.my_generated
    }

    /// `UpdateTolesp()` (`ChFi3d_Builder_C2.cxx:814-824`): lowers `tolesp` to
    /// the `tolesp` of every stripe's guideline.
    pub fn update_tolesp(&mut self) {
        for stripe in &self.my_list_stripe {
            let Some(spine) = spine_of(stripe) else {
                continue;
            };
            let Ok(guard) = spine.read() else {
                continue;
            };
            let current_stripe_tolesp = guard.tolesp();
            if self.tolesp > current_stripe_tolesp {
                self.tolesp = current_stripe_tolesp;
            }
        }
    }

    /// `IsBound(E)` / `Bind(E, L)` / `ChangeFind(E)` of `myEVIMap`
    /// (`ChFi3d_Builder_2.cxx:3801-3806`, `ChFi3d_Builder_C1.cxx:2817-2834`):
    /// the index list of `s`, created empty on first sight of `s`.
    pub fn change_evi(&mut self, s: &TopoShape) -> &mut Vec<i32> {
        if let Some(k) = self
            .my_evi_map
            .iter()
            .position(|(key, _)| is_same(key, s))
        {
            return &mut self.my_evi_map[k].1;
        }
        self.my_evi_map.push((s.clone(), Vec::new()));
        let last = self.my_evi_map.len() - 1;
        &mut self.my_evi_map[last].1
    }

    /// `myEdgeFirstFace.Bind(Ec, Face)` (`ChFi3d_Builder_1.cxx:927`, `:1040`,
    /// `:1146`).
    pub fn bind_edge_first_face(&mut self, e: &TopoShape, f: &TopoShape) {
        if let Some(k) = self
            .my_edge_first_face
            .iter()
            .position(|(key, _)| is_same(key, e))
        {
            self.my_edge_first_face[k].1 = f.clone();
            return;
        }
        self.my_edge_first_face.push((e.clone(), f.clone()));
    }

    /// `myEdgeFirstFace(Ec)` (`ChFi3d_Builder_1.cxx:1089`,
    /// `ChFi3d_Builder_2.cxx:836`, `:879`).
    pub fn edge_first_face(&self, e: &TopoShape) -> Option<&TopoShape> {
        self.my_edge_first_face
            .iter()
            .find(|(key, _)| is_same(key, e))
            .map(|(_, f)| f)
    }

    /// The `I`-th `myListStripe` entry, `I` `1`-based.
    fn stripe_at(&self, ic: usize) -> Option<StripeHandle> {
        self.my_list_stripe.get(ic.checked_sub(1)?).cloned()
    }

    /// `myEVIMap.Find(EouV)`, without creating the entry.
    fn evi_indices(&self, s: &TopoShape) -> Option<&[i32]> {
        self.my_evi_map
            .iter()
            .find(|(key, _)| is_same(key, s))
            .map(|(_, list)| list.as_slice())
    }
}
