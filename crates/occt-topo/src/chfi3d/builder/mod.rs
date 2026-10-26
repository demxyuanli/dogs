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
//!   `perform_set_of_surf/splice.rs`,
//! - [`ChFi3dBuilder::perform_set_of_k_part`] ([`perform_set_of_k_part`],
//!   `_2.cxx:3004-3281`), the first of the two walks `PerformSetOfSurf` runs:
//!   the particular-corner scan of the spine, its `SurfData` sequence and the
//!   `ElSpine` list `PerformSetOfKGen` then walks, with
//!   [`crate::chfi3d::builder_0_helpers::ch_fi3d_k_particular`]
//!   (`ChFi3d_Builder_0.cxx:604-683`), [`tgt_kp`] (`_2.cxx:373-398`) and
//!   [`crate::chfi3d::builder_0_helpers::ch_fi3d_solid_index`]
//!   (`ChFi3d_Builder_0.cxx:2300-2321`) under it,
//! - [`ChFi3dBuilder::perform_set_of_k_gen`] ([`perform_set_of_k_gen`],
//!   `_2.cxx:3298-3881`), the walk of every `ElSpine` of the spine through
//!   `PerformSetOfSurfOnElSpine`, then — outside the simultaneous run — the
//!   twisted-section repair (`:3318-3648`) with `ChFi3d_IntTraces`
//!   ([`crate::chfi3d::builder_0_helpers::ch_fi3d_int_traces`],
//!   `ChFi3d_Builder_0.cxx:1234-1447`), the six `ChFi3d_mkbound` overloads
//!   ([`crate::chfi3d::builder_0_helpers::mkbound_curve`] and its siblings,
//!   `:1750-1967`) and `ChFi3d_ReparamPcurv`
//!   ([`crate::chfi3d::builder_0_helpers::ch_fi3d_reparam_pcurv`],
//!   `:1472-...`) on the `GeomFill_ConstrainedFilling` patch, the point
//!   tolerance update (`:3650-3754`) and the `myEVIMap` edge / new-face links
//!   (`:3755-3878`),
//! - [`ChFi3dBuilder::perform_set_of_surf`] (`_2.cxx:3882-3900`): the solid
//!   index, `PerformSetOfKPart` when the spine is not split, then the two above
//!   and `ChFi3d_MakeExtremities` (`:1766-2046`, with `ChFi3d_SingularExtremity`
//!   `:1635-1743` and `IsFree` `:1744-1765`) in [`make_extremities`],
//! - `ChFi3d_PerformElSpine` (`ChFi3d_Builder_0.cxx:4797-5327`) as
//!   [`crate::chfi3d::builder_0_helpers::ch_fi3d_perform_el_spine`], the tail
//!   `PerformSetOfKPart` runs on every `ElSpine` once the KPart scan is done,
//!   with `ChFi3d_ApproxByC2` ([`crate::chfi3d::builder_helpers::approx_by_c2`])
//!   and `CurveCleaner`
//!   ([`crate::chfi3d::builder_0_helpers::curve_cleaner`]) under it,
//! - `ChFi3d_Builder::Compute` (`ChFi3d_Builder.cxx:178-680`, [`compute`]) with
//!   `CompleteDS` (`:80-136`, [`complete_ds`]) and `ExtentAnalyse`
//!   (`:144-174`, [`extent_analyse`]) — the pipeline that walks
//!   `PerformSetOfSurf` over every contour, dispatches the corners, indexes the
//!   solids / shells in the DS, cross-checks the stripes with
//!   [`crate::chfi3d::builder_0_helpers::ch_fi3d_stripe_edge_inter`]
//!   (`ChFi3d_Builder_0.cxx:3346-3443`), runs `ChFi3d_FilDS`, `CompleteDS`, the
//!   `myCoup` reconstruction and the `SameParameter` pass —,
//! - `PerformFilletOnVertex` (`:759-920`, [`perform_fillet_on_vertex`]) and
//!   `PerformSingularCorner` (`:682-757`, [`perform_singular_corner`]), the
//!   corner dispatch and the degenerated-corner arc builder; both reach the
//!   corner bodies through [`ChFi3dCornerOps`],
//!
//! UNPORTED — the corner bodies the two dispatchers call:
//!
//! - `PerformFirstSection` (`ChFi3d_Builder.hxx:577-588`) is declared in
//!   [`ChFi3dSurfOps`] but has no implementor (pure virtual in OCCT, implemented
//!   by `ChFi3d_FilBuilder` / `ChFi3d_ChBuilder`), together with `SimulKPart`
//!   (`hxx:201`, declared with the `ChFi3d_Builder_NotImp.cxx` answer as its
//!   default), `PerformTwoCornerbyInter`,
//! - `ChFi3d_FilDS` (`ChFi3d_Builder_0.cxx:2565-3351`, [`fil_ds`]): the whole
//!   fillet-side DS filling, blocked on `ChFi3d_FilPointInDS` /
//!   `FilVertexInDS` / `FilCurveInDS` / `Contains` / `CutEdge` /
//!   `findIndexPoint` / `ChFi3d_ComputePCurv`, none of which is ported,
//! - `SplitKPart` (`ChFi3d_Builder_SpKP.cxx:749-1160`) in [`split_k_part`]: it
//!   needs `Geom2dHatch_Hatcher` / `Geom2dHatch_Intersector` and the `Tri` /
//!   `FillSD` / `CpSD` statics of the same file,
//! - `ChFiKPart_ComputeData::ComputeCorner` — the three overloads are ported
//!   ([`crate::chfi3d::chfi_k_part`]), but their OCCT callers
//!   `PerformTwoCorner` / `PerformThreeCorner` (`ChFi3d_Builder_C1.cxx`) are
//!   not,
//! - `PerformSetOfSurfOnElSpine` (`ChFi3d_Builder_2.cxx:2207-3003`) — ported and
//!   now called from `PerformSetOfKGen`, but its obstacle arms need the
//!   edge-contact `PerformSurf` / `SimulSurf` of the subclasses, so they stay on
//!   the `ChFi3d_Builder_NotImp.cxx` `throw` defaults for now —,
//! - the corner bodies reached through [`ChFi3dCornerOps`]:
//!   `PerformOneCorner` (`ChFi3d_Builder_C1.cxx:611-1827`),
//!   `PerformIntersectionAtEnd` (`:1828-3770`), `PerformMoreSurfdata`
//!   (`:3771-4531`), `IntersectMoreCorner` (`:4708-...`), and
//!   `PerformMoreThreeCorner` (`:1237-1652`), all of them `ChFiKPart` /
//!   `BRepBlend` consumers; the six pure virtuals (`PerformTwoCorner`,
//!   `PerformThreeCorner`, the three `Extent*Corner`, `SetRegul`) stay in the
//!   subclasses (`ChFi3d_FilBuilder` / `ChFi3d_ChBuilder`).
//! - `ComputeCurve2d` (`ChFi3d_Builder_C1.cxx:456-483`) - UNPORTED with its
//!   `BRepAlgo_NormalProjection` dependency, see
//!   [`crate::chfi3d::builder_c1`].
//!
//! Ported from `ChFi3d_Builder_C1.cxx` (see [`crate::chfi3d::builder_c1`] for
//! the file statics):
//!
//! - [`ChFi3dBuilder::find_face`] (`:4533-4598`, both overloads),
//! - [`ChFi3dBuilder::more_surfdata`] (`:4601-4705`), the `const` member
//!   `PerformFilletOnVertex` branches on.
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

use occt_core::kernel::geomabs::Shape;

mod accessors;
mod call_perform_surf;
mod complete_data;
mod complete_ds;
mod compute;
mod compute_data;
mod compute_data_face_face;
mod corner_ops;
mod extent_analyse;
mod fil_ds;
mod find_face;
mod make_extremities;
mod more_surfdata;
mod perform_element;
mod perform_extremity;
mod perform_fillet_on_vertex;
mod perform_helpers;
mod perform_one_corner;
mod perform_set_of_k_part;
mod perform_set_of_k_gen;
mod perform_set_of_surf;
mod perform_singular_corner;
mod search_face;
mod simul_data;
mod simul_data_face_face;
mod start_sol;
mod start_sol_obstacle;
mod state;
mod store_data;
mod stripe_surgery;
mod split_k_part;
mod surf_ops;
mod update_line;

pub use corner_ops::ChFi3dCornerOps;
pub use store_data::store_data;
pub use surf_ops::ChFi3dSurfOps;
pub use update_line::update_line;

use crate::abs::ShapeType;
use crate::chfi_ds::{ChFiDSMap, Regul, SpineHandle, StripeHandle, StripeMap};
use crate::shape::TopoShape;
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

    // The result / fault / contour accessors live in accessors.rs and the
    // state transitions in state.rs; both are impl ChFi3dBuilder blocks in
    // child modules, so the struct fields stay private to builder.
}
