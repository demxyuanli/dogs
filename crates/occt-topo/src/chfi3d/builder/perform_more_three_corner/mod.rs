//! `ChFi3d_Builder::PerformMoreThreeCorner` - the corner of a vertex with more
//! than three edges.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_CnCrn.cxx:1237-3927`
//! (OCCT 8.0.0), plus the file statics it calls (`cxx:128-1230`).
//!
//! The OCCT body is one function of 2691 lines that fills a block of parallel
//! arrays (`CD` / `jf` / `Index` / `sens` / `Evive` / `Fvive` / `numfa` / `p` /
//! ...) and then walks them. The port keeps those arrays in
//! [`MoreThreeCorner`] so the walk can be split by OCCT line range, and takes
//! them a stage at a time:
//!
//! | stage | OCCT lines | module | state |
//! |---|---|---|---|
//! | vertex counts, sewing edge, edge/face parse, `regul` | `cxx:1237-1581` | [`topology_parse`] | ported |
//! | `deuxconges` / `deuxcgnontg` / `isOnSameDiff` flags | `cxx:1583-1676` | [`two_conges`] | ported |
//! | `Evive` repair from the neighbouring common points | `cxx:1679-1782` | [`evive_repair`] | ported |
//! | first free edge, minimum recoil distance | `cxx:1786-1830` | [`free_edge_recoil`] | ported |
//! | `ChFi3d_SearchFD` between stripes and the `p` / `i` fill | `cxx:1832-2063` | [`search_fd_stage`] | ported |
//! | `dist1` / `dist2` scan and the parameter offset | `cxx:2064-2196` | [`vertex_distance`] | ported |
//! | common-point edge limiting, free-border parameter | `cxx:2198-2366` | [`common_point_limit`] | ported |
//! | regular-edge projection, `CurveHermite` | `cxx:2368-2541` | [`regular_edge_projection`], [`curve_hermite`] | ported |
//! | multi-face `CurveHermite` over `Fproj` / `Ecom` | `cxx:2542-2640` | [`multi_face_hermite`] | ported |
//! | `PerformTwoCornerSameExt` at `nconges == 2` | `cxx:2642-2705` | [`two_corner_same_ext`] | ported |
//! | stripe `GeomPlate` constraints, living-edge points | `cxx:2707-2905` | [`plate_stripe_constraints`] | ported |
//! | intermediary connecting curves | `cxx:2907-3260` | [`plate_intermediary`] | ported |
//! | multi-face curves, tangent free borders | `cxx:3262-3488` | [`plate_multi_face`] | ported |
//! | `PSurf.Perform` / `MakeApprox` / DS surface | `cxx:3490-3626` | [`plate_solve_approx`] | ported |
//! | plate tails (p-curves, `ComputeArete`, wires) | `cxx:3628-3927` | [`plate_tails`] | ported |
//!
//! Translation notes:
//! - `TopoDS_Shape` has no null value in the port, so the `Evive` / `Fvive`
//!   cells are `Option<TopoShape>` and `CD` is `Option<StripeHandle>`,
//! - OCCT's 1-based `NCollection_Sequence` indices (`CD`, `jf`, `Index`,
//!   `numfa`, `p`, `i`, `indpoint`) are kept 1-based by the struct so that
//!   `SurfIndex` / `ChFi3d_IndexOfSurfData` answer the same values,
//! - `ChFi3d_cherche_edge` / `ChFi3d_cherche_face1` write their answer through
//!   out-parameters and keep the previous value when they find nothing; only a
//!   still-null answer raises `Standard_ConstructionError`
//!   (`ChFi3d_Builder_0.cxx:5458-5461`, `:5350-5353`). The port's helpers
//!   answer `Result` / `Option`, so the two call sites here restore that
//!   out-parameter behaviour explicitly ([`topology_parse::cherche_edge_into`]),
//! - `PerformMoreThreeCorner` is a plain member in OCCT (`hxx:564`, not
//!   virtual), so it is an inherent method of [`super::ChFi3dBuilder`].

mod calcul;
mod common_point_limit;
mod curve_hermite;
mod evive_repair;
mod free_edge_recoil;
mod multi_face_hermite;
mod plate_helpers;
mod plate_intermediary;
mod plate_multi_face;
mod plate_solve_approx;
mod plate_stripe_constraints;
mod plate_tails;
mod regular_edge_projection;
mod remove_surf_data;
mod search_fd_stage;
mod surface_normal;
mod surf_index;
mod topology_parse;
mod two_conges;
mod two_corner_same_ext;
mod vertex_distance;

use std::sync::Arc;

use occt_core::gp::GpVec;
use occt_geom::geom_plate::BuildPlateSurface;
use occt_geom::Curve;
use occt_geom2d::curve::Curve2d;

use crate::chfi_ds::{CommonPoint, StripeHandle};
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::is_same;

use crate::abs::Orientation;
use occt_geom::Surface;

use super::ChFi3dBuilder;

/// The locals of `PerformMoreThreeCorner` that the ported stages produce.
///
/// Fields are declared here in OCCT order so each stage reads the value the
/// body it translates wrote; the `dead_code` allow covers the cells OCCT itself
/// writes without a reader on every path.
#[allow(dead_code)]
pub struct MoreThreeCorner {
    /// `V1`
    pub vtx: TopoShape,
    /// `nedge` — after the free-border adjustment of `cxx:1256-1263`.
    pub nedge: i32,
    /// `size` — `nedge * 2`.
    pub size: usize,
    /// `bordlibre`
    pub bord_libre: bool,
    /// `edgelibre1`
    pub edge_libre1: Option<Edge>,
    /// `edgelibre2`
    pub edge_libre2: Option<Edge>,
    /// `droit` — the two free borders are opposite.
    pub droit: bool,
    /// `CD`
    pub cd: Vec<Option<StripeHandle>>,
    /// `jf`
    pub jf: Vec<i32>,
    /// `Index`
    pub index: Vec<i32>,
    /// `Indice`
    pub indice: Vec<i32>,
    /// `sens`
    pub sens: Vec<i32>,
    /// `indcurve3d`
    pub indcurve3d: Vec<i32>,
    /// `numfa`
    pub numfa: Vec<Vec<i32>>,
    /// `Order`
    pub order: Vec<i32>,
    /// `i`
    pub i: Vec<Vec<i32>>,
    /// `indpoint`
    pub indpoint: Vec<Vec<i32>>,
    /// `oksea`
    pub oksea: Vec<bool>,
    /// `sharp`
    pub sharp: Vec<bool>,
    /// `regul`
    pub regul: Vec<bool>,
    /// `p`
    pub p: Vec<Vec<f64>>,
    /// `errapp`
    pub errapp: Vec<f64>,
    /// `Evive`
    pub evive: Vec<Option<TopoShape>>,
    /// `Fvive`
    pub fvive: Vec<Vec<Option<TopoShape>>>,
    /// `ponctuel`
    pub ponctuel: Vec<bool>,
    /// `samedge`
    pub samedge: Vec<bool>,
    /// `moresurf`
    pub moresurf: Vec<bool>,
    /// `libre`
    pub libre: Vec<bool>,
    /// `tangentregul`
    pub tangentregul: Vec<bool>,
    /// `isG1`
    pub is_g1: Vec<bool>,
    /// `couture` — the vertex is on a sewing edge.
    pub couture: bool,
    /// `facecouture`
    pub face_couture: Option<Face>,
    /// `edgecouture`
    pub edge_couture: Option<Edge>,
    /// `SumFaceNormalAtV1`
    pub sum_face_normal_at_v1: GpVec,
    /// `deuxconges` — two fillets on two tangent living edges (`cxx:1586`).
    pub deuxconges: bool,
    /// `deuxcgnontg` — three edges, two fillets, not `deuxconges` (`cxx:1628`).
    pub deuxcgnontg: bool,
    /// `isOnSameDiff` — one stripe `OnSame` and one `OnDiff` (`cxx:1640`).
    pub is_on_same_diff: bool,
    /// `distmin` — the minimum recoil distance (`cxx:1799-1830`).
    pub distmin: f64,
    /// `distmini` — set when a `p` / `i` pair strays past `distmin`
    /// (`cxx:1799`, `:1986-1989`).
    pub distmini: bool,
    /// `jfp` — a function-level local of `PerformMoreThreeCorner`
    /// (`cxx:1308`), written by the `Evive` repair (`cxx:1698`) and by the
    /// extremity fill (`cxx:2057`), read by the `i1 != Index(ic)` /
    /// `i2 != Index(icplus)` face walks (`cxx:1950`, `:1964`).
    pub jfp: i32,
    /// `Fproj` (`cxx:1318`, a `NCollection_Sequence<TopoDS_Shape>`) — the faces
    /// the `i1 != Index` / `i2 != Index` walks append; read from `cxx:2595` on
    /// by [`multi_face_hermite`] and [`plate_multi_face`].
    pub fproj: Vec<TopoShape>,
    /// `inters` — every slot found an intersection (`cxx:1833`, `:1996-1999`).
    pub inters: bool,
    /// `dist1` — `NCollection_Array1<double>(0, size)` (`cxx:2065`): the
    /// distance from the corner vertex to the surface point of the
    /// `i(ic, icmoins)` / `p(ic, icmoins)` pair.
    pub dist1: Vec<f64>,
    /// `dist2` — the `i(ic, icplus)` / `p(ic, icplus)` mirror (`cxx:2066`).
    pub dist2: Vec<f64>,
    /// `tolcp` — the largest `ChFiDS_CommonPoint::Tolerance()` (< 1) the
    /// common-point limiting of `cxx:2201-2315` accepted; read by the
    /// `GeomPlate` section (`cxx:2875`, `:2889`), [`plate_stripe_constraints`].
    pub tolcp: f64,
    /// `nbface` — the function-scope counter of `cxx:2373`, written by the
    /// regular-edge projection (`cxx:2430-2437`, `:2538`) and read by the
    /// `GeomPlate` section (`cxx:2628`, `:2635`), [`plate_multi_face`].
    pub nbface: i32,
    /// `error` — the smallest projected-edge tolerance (`cxx:2374`,
    /// `:394-458`), read by the `GeomPlate` section, [`plate_intermediary`].
    pub hermite_error: f64,
    /// `proj2d1` (`cxx:2376`) — `NCollection_Array1(0, size)`; read at
    /// `cxx:3082` ([`plate_intermediary`]).
    pub proj2d1: Vec<Option<Arc<dyn Curve2d>>>,
    /// `proj2d2` (`cxx:2377`); read at `cxx:3069` ([`plate_intermediary`]).
    pub proj2d2: Vec<Option<Arc<dyn Curve2d>>>,
    /// `cproj1` (`cxx:2378`); read at `cxx:3076` / `:3083` ([`plate_intermediary`]).
    pub cproj1: Vec<Option<Arc<dyn Curve>>>,
    /// `cproj2` (`cxx:2379`); read at `cxx:3063` / `:3070` ([`plate_multi_face`]).
    pub cproj2: Vec<Option<Arc<dyn Curve>>>,
    /// `Eproj` (`cxx:2545`) — appended by `CurveHermite` (`cxx:441-471`) and
    /// read from `cxx:3294` on ([`plate_multi_face`]).
    pub eproj: Vec<Option<Edge>>,
    /// `Ecom` (`cxx:2544`) — the shared edges between consecutive `Fproj`
    /// faces, appended by the multi-face stage (`cxx:2606-2635`) and read by
    /// `CurveHermite` (`cxx:475-526`) from its first entry on.
    pub ecom: Vec<Edge>,
    /// `parcom` (`cxx:2546`) — the first-vertex parameter of every `Ecom`
    /// edge, indexed like `Ecom`.
    pub parcom: Vec<f64>,
    /// `IVtx` (`cxx:1310`) — the DS shape index of `V1`, written by the
    /// `deuxconges` arm of the living-edge points (`cxx:2865`) and read by the
    /// `deuxconges` arcs of the plate storage (`cxx:3589`).
    pub ivtx: i32,
    /// `Surf` (`cxx:3522`) — the BSpline the plate approximation produced.
    pub plate_surface: Option<Arc<dyn Surface>>,
    /// `Isurf` (`cxx:3542`) — the DS surface index of [`Self::plate_surface`].
    pub isurf: i32,
    /// `apperror` (`cxx:3525-3526`) — `Mapp.CriterionError() * 1.1`, read by
    /// the tails (`cxx:3698`, `:3870`).
    pub apperror: f64,
    /// `orsurfdata` (`cxx:3617`) — the orientation of the first stripe's surf
    /// data the plate is built for.
    pub or_surfdata: Orientation,
    /// `orplate` (`cxx:3618`) — the `PlateOrientation` answer.
    pub or_plate: Orientation,
    /// `isfirst` (`cxx:3598`) — `sens(0) == 1`, re-read by the p-curve
    /// orientation of the tails (`cxx:3658-3665`).
    pub is_first_plate: bool,
    /// `n3d` (`cxx:1310`) — the running curve counter of the `GeomPlate`
    /// section. The plate stages take it by reference; the tails re-read it
    /// from here because `cxx:3733` resets it to `0`.
    pub n3d_curves: usize,
}

impl Default for MoreThreeCorner {
    fn default() -> Self {
        Self {
            vtx: TopoShape::new(crate::abs::ShapeType::Vertex),
            nedge: 0,
            size: 0,
            bord_libre: false,
            edge_libre1: None,
            edge_libre2: None,
            droit: false,
            cd: Vec::new(),
            jf: Vec::new(),
            index: Vec::new(),
            indice: Vec::new(),
            sens: Vec::new(),
            indcurve3d: Vec::new(),
            numfa: Vec::new(),
            order: Vec::new(),
            i: Vec::new(),
            indpoint: Vec::new(),
            oksea: Vec::new(),
            sharp: Vec::new(),
            regul: Vec::new(),
            p: Vec::new(),
            errapp: Vec::new(),
            evive: Vec::new(),
            fvive: Vec::new(),
            ponctuel: Vec::new(),
            samedge: Vec::new(),
            moresurf: Vec::new(),
            libre: Vec::new(),
            tangentregul: Vec::new(),
            is_g1: Vec::new(),
            couture: false,
            face_couture: None,
            edge_couture: None,
            sum_face_normal_at_v1: GpVec::zero(),
            deuxconges: false,
            deuxcgnontg: false,
            is_on_same_diff: false,
            distmin: 0.0,
            distmini: false,
            jfp: 0,
            fproj: Vec::new(),
            inters: false,
            dist1: Vec::new(),
            dist2: Vec::new(),
            tolcp: 0.0,
            nbface: 0,
            hermite_error: 0.0,
            proj2d1: Vec::new(),
            proj2d2: Vec::new(),
            cproj1: Vec::new(),
            cproj2: Vec::new(),
            eproj: Vec::new(),
            ecom: Vec::new(),
            parcom: Vec::new(),
            ivtx: 0,
            plate_surface: None,
            isurf: 0,
            apperror: 0.0,
            or_surfdata: Orientation::Forward,
            or_plate: Orientation::Forward,
            is_first_plate: false,
            n3d_curves: 0,
        }
    }
}

impl MoreThreeCorner {
    /// The empty state; the arrays are sized once `nedge` is known.
    pub fn new() -> Self {
        Self::default()
    }

    /// The `NCollection_Array1(0, size)` / `NCollection_Array2(0, size, ...)`
    /// allocations of `cxx:1265-1288` together with the initialization loop of
    /// `cxx:1290-1304` (every cell zero / `false`, `indpoint` both columns `0`).
    pub fn allocate(&mut self, size: usize) {
        let n = size + 1;
        self.size = size;
        self.cd = vec![None; n];
        self.jf = vec![0; n];
        self.index = vec![0; n];
        self.indice = vec![0; n];
        self.sens = vec![0; n];
        self.indcurve3d = vec![0; n];
        self.numfa = vec![vec![0; n]; n];
        self.order = vec![0; n];
        self.i = vec![vec![0; n]; n];
        self.indpoint = vec![vec![0; 2]; n];
        self.oksea = vec![false; n];
        self.sharp = vec![false; n];
        self.regul = vec![false; n];
        self.p = vec![vec![0.0; n]; n];
        self.errapp = vec![0.0; n];
        self.evive = vec![None; n];
        self.fvive = vec![vec![None; n]; n];
        self.ponctuel = vec![false; n];
        self.samedge = vec![false; n];
        self.moresurf = vec![false; n];
        self.libre = vec![false; n];
        self.tangentregul = vec![false; n];
        self.is_g1 = vec![false; n];
        self.dist1 = vec![0.0; n];
        self.dist2 = vec![0.0; n];
    }

    /// The non-null cells of `Evive`, in array order — the exclusion list
    /// `ChFi3d_cherche_edge(V1, Evive, Fcur, Enext, VV)` compares against
    /// (`Evive` is handed over whole, and `IsSame` against a null shape is
    /// always false).
    pub(super) fn excluded_edges(&self) -> Vec<TopoShape> {
        self.evive.iter().flatten().cloned().collect()
    }

    /// `TopoDS::Edge(Evive.Value(ic))`.
    pub(super) fn edge_at(&self, ic: usize) -> Result<Edge, String> {
        let shape = self.evive[ic]
            .clone()
            .ok_or_else(|| "PerformMoreThreeCorner: Evive is null".to_string())?;
        Edge::wrap(shape)
            .ok_or_else(|| "PerformMoreThreeCorner: Evive is not an edge".to_string())
    }

    /// `CD.Value(slot)->SetOfSurfData()->Value(Index.Value(slot))->Vertex(
    /// isfirst, on_s)` read as a copy — the `ChangeVertex(isfirst, jfp)` of
    /// `cxx:1696-1701`, `:1741-1745`, `:2217-2223`, `:2224-2231`,
    /// `:2341`, `:2355`.
    pub(super) fn common_point_at(
        &self,
        slot: usize,
        isfirst: bool,
        on_s: i32,
    ) -> Result<CommonPoint, String> {
        let stripe = self.cd[slot]
            .clone()
            .ok_or_else(|| "PerformMoreThreeCorner: CD is null".to_string())?;
        let index = self.index[slot];
        if index < 1 {
            return Err("PerformMoreThreeCorner: Index is null".to_string());
        }
        let guard = stripe
            .read()
            .map_err(|_| "PerformMoreThreeCorner: poisoned stripe".to_string())?;
        let surf = guard
            .set_of_surf_data()
            .get(index as usize - 1)
            .ok_or_else(|| {
                format!("PerformMoreThreeCorner: SurfData {index} is out of range")
            })?;
        Ok(surf.vertex(isfirst, on_s).clone())
    }

    /// `CD.Value(slot)` as a borrowed stripe (`cxx:2479`, `:2486`).
    pub(super) fn stripe_at(&self, slot: usize) -> Result<&StripeHandle, String> {
        self.cd
            .get(slot)
            .and_then(|s| s.as_ref())
            .ok_or_else(|| format!("PerformMoreThreeCorner: CD({slot}) is null"))
    }

    /// `TopoDS::Face(DStr.Shape(index))` (`cxx:2426`, `:2476`): the default
    /// `FindKeep` of `TopOpeBRepDS_DataStructure::Shape` is `true`
    /// (`hxx:254`), so a shape that is not kept answers nothing here.
    pub(super) fn face_at(
        &self,
        ds: &crate::topopebrep::ds::DataStructure,
        index: i32,
    ) -> Result<Option<Face>, String> {
        let shape = ds
            .shape(index, true)
            .ok_or_else(|| format!("PerformMoreThreeCorner: DStr.Shape({index}) is null"))?;
        Ok(Face::wrap(shape))
    }

    /// `Fvive(row, column)` (`cxx:2429`, `:2450`, `:2477`) as a face.
    pub(super) fn fvive_face(&self, row: usize, column: usize) -> Result<Option<Face>, String> {
        let cell = self
            .fvive
            .get(row)
            .and_then(|r| r.get(column))
            .ok_or_else(|| format!("PerformMoreThreeCorner: Fvive({row}, {column}) is missing"))?;
        Ok(cell.clone().and_then(Face::wrap))
    }

    /// `E.IsSame(edgelibre1) || E.IsSame(edgelibre2)`; a null free border is
    /// never the same shape.
    pub(super) fn is_free_border(&self, e: &Edge) -> bool {        self.edge_libre1
            .as_ref()
            .is_some_and(|l| is_same(&e.0, &l.0))
            || self.edge_libre2
                .as_ref()
                .is_some_and(|l| is_same(&e.0, &l.0))
    }
}

/// `spine->Edges(1)` / `spine->Edges(spine->NbEdges())` — the edge of the
/// spine that carries the corner vertex.
pub(super) fn spine_edge_at(stripe: &StripeHandle, first: bool) -> Result<Edge, String> {
    let spine = {
        let guard = stripe
            .read()
            .map_err(|_| "PerformMoreThreeCorner: poisoned stripe".to_string())?;
        guard
            .spine()
            .cloned()
            .ok_or_else(|| "PerformMoreThreeCorner: the stripe has no spine".to_string())?
    };
    let guard = spine
        .read()
        .map_err(|_| "PerformMoreThreeCorner: poisoned spine".to_string())?;
    let i = if first { 1 } else { guard.nb_edges() };
    let shape = guard
        .edges(i)
        .cloned()
        .ok_or_else(|| format!("PerformMoreThreeCorner: spine edge {i} is missing"))?;
    Edge::wrap(shape)
        .ok_or_else(|| format!("PerformMoreThreeCorner: spine entry {i} is not an edge"))
}

impl ChFi3dBuilder {
    /// The three tables the `GeomPlate` section writes besides the DS — the
    /// data structure, `myEVIMap` and `myRegul` — borrowed together so that a
    /// stage can append to all of them at once.
    pub(super) fn plate_tables(
        &mut self,
    ) -> (
        &mut crate::topopebrep::ds::DataStructure,
        &mut Vec<(TopoShape, Vec<i32>)>,
        &mut Vec<crate::chfi_ds::Regul>,
    ) {
        (self.my_ds.change_ds(), &mut self.my_evi_map, &mut self.my_regul)
    }

    /// `PerformMoreThreeCorner(Jndex, nconges)` (`cxx:1237-3927`): every stage
    /// of the body is ported; the driver below walks them in the OCCT order.
    pub fn perform_more_three_corner(
        &mut self,
        index: usize,
        nconges: usize,
    ) -> Result<(), String> {
        let mut state = MoreThreeCorner::new();
        self.more_three_corner_topology(index, &mut state)?;

        // `cxx:1583-1676` and the `PerformOneCorner(Jndex, true)` it triggers
        // (`cxx:1671-1675`).
        state.two_conges(nconges)?;
        if state.is_on_same_diff {
            // The OCCT call does not return here: the `Evive` repair below is
            // then skipped only through its own `!isOnSameDiff` guard.
            self.perform_one_corner(index, true)?;
        }

        // `cxx:1679-1782`.
        if !state.deuxconges && !state.is_on_same_diff {
            state.repair_evive(&mut self.my_ds, &self.my_ef_map)?;
        }

        // `cxx:1786-1830`.
        state.free_edge_and_recoil()?;

        // `cxx:1832-2063`.
        let vtx = state.vtx.clone();
        state.intersect_stripes(&self.my_ds, &vtx)?;

        // `cxx:2064-2196`.
        state.vertex_distance_offsets(&self.my_ds, &vtx, nconges)?;

        // `cxx:2198-2366`.
        state.limit_by_common_points(&vtx, nconges)?;

        // `cxx:2368-2541`.
        state.project_regular_edges(self.my_ds.ds())?;

        // `cxx:2542-2640`.
        state.multi_face_hermite(self.my_ds.ds())?;

        // `cxx:2642-2705`: two fillets on the vertex that are not tangent - the
        // corner is settled by the intersection of the two stripes whenever the
        // two ends of that intersection land on the common points they share,
        // and `PerformTwoCornerSameExt` returns (the `return` of `cxx:2698`).
        if state.two_corner_same_ext(self.my_ds.change_ds(), nconges)? {
            return Ok(());
        }

        // `cxx:2707-2718`: the plate solver. `degree = 3`, `nbcurvpnt = 10`,
        // `nbiter = 1`, and the six-argument constructor leaves `TolCurv` at
        // its `0.1` default with `Anisotropie` false
        // (`GeomPlate_BuildPlateSurface.hxx:66-73`).
        let tolapp3d = self.tolapp3d;
        let angular = self.angular;
        let tol2d = self.tol2d;
        let mut psurf = BuildPlateSurface::new(3, 10, 1, tol2d, tolapp3d, angular, false)?;

        let (ds, evi, regul) = self.plate_tables();
        let mut n3d = 0usize;

        // `cxx:2716-2852`, `:2854-2905`, `:2907-3260`, `:3262-3407` and
        // `:3409-3488`: the constraints of the plate solver.
        state.plate_stripe_constraints(ds, &mut psurf, &mut n3d, tolapp3d, angular)?;
        // `tolcp` is the largest common-point tolerance of `cxx:2201-2315`.
        let tolcp = state.tolcp;
        state.plate_living_point_indices(ds, tolcp)?;
        state.plate_intermediary_curves(ds, &mut psurf, &mut n3d, nconges, tolapp3d, angular)?;
        state.plate_multi_face_curves(ds, &mut psurf, &mut n3d, tolapp3d, angular)?;
        state.plate_droit(ds, &mut psurf, &mut n3d, tolapp3d, angular)?;
        state.n3d_curves = n3d;

        // `cxx:3494-3626` and, when the plate failed, `cxx:3891-3926`.
        let plate_done = state.plate_solve_approx(ds, evi, &mut psurf)?;
        if plate_done {
            // `cxx:3628-3890`.
            state.plate_tails(ds, regul, &psurf)?;
        } else {
            // `else { done = false; hasresult = true; ... }` (`cxx:3891-3926`).
            // The two flags are only written on this arm; the success path
            // leaves them to `Compute`.
            state.plate_partial_tail(ds)?;
            self.done = false;
            self.has_result = true;
        }

        Ok(())
    }
}
