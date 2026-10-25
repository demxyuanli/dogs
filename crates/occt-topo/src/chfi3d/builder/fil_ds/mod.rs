//! `ChFi3d_FilDS` - fill the DS with one stripe's fillet surface data.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:2565-3341`
//! (OCCT 8.0.0). `ChFi3d_Builder::Compute` calls it once per stripe after the
//! `ChFi3d_StripeEdgeInter` sweep (`ChFi3d_Builder.cxx:391-395`) and it throws
//! straight out of `Compute` when it fails (there is no `try` around that
//! call).
//!
//! Layout of the translation:
//! - this file: the entry point, the `ChFiDS_Spine` / `SeqFil` prelude, the
//!   `FilDsState` locals of `cxx:2571-2607` and the small shared helpers,
//! - [`nullify`]: the "Nullify degenerated ChFi/Faces interferences" prelude
//!   (`cxx:2611-2653`),
//! - [`step`]: the driver of the `for (j = 1; ...)` loop (`cxx:2655-3340`),
//! - [`initial`]: the `j == 1` initial processing (`cxx:2779-2887`) and the
//!   `j > 1` inter-fillet interference (`cxx:2888-2952`),
//! - [`faces`]: the fillet/faces interferences on S1 and S2
//!   (`cxx:2954-3158`), including the end-point index choice,
//! - [`tail`]: the last curve / middle curve (`cxx:3160-3240`) and the
//!   degenerated end plus the point/edge interferences (`cxx:3242-3339`),
//! - [`crate::chfi3d::builder_0_helpers::fil_interf`]: `ChFi3d_FilPointInDS` /
//!   `FilVertexInDS` / `FilCurveInDS` / `ChFi3d_Contains` / `ChFi3d_Orientation`
//!   / `QueryAddVertexInEdge` / `CutEdge` / `findIndexPoint`,
//! - [`crate::chfi3d::builder_0_helpers::point_in_ds`]: `ChFi3d_IndexPointInDS`,
//! - [`crate::chfi3d::builder_0_helpers::compute_pcurv`] /
//!   [`crate::chfi3d::builder_0_helpers::compute_arete`]: `ChFi3d_ComputePCurv`
//!   and `ChFi3d_ComputeArete`.
//!
//! Translation notes:
//! - `DStr` is `self.my_ds.change_ds()` and `reglist` is `self.my_regul`, so
//!   the two OCCT parameters that live in the builder are not passed around,
//! - `Fd->ChangeVertex(...)` in the nullify prelude is `Arc::make_mut` on the
//!   stripe's sequence entry; every other use of `Fd` in the loop is read-only,
//!   exactly like the `const occ::handle<ChFiDS_SurfData>& Fd` of `cxx:2658`,
//! - the loop state of `cxx:2571-2607` (`Ipoin1`, `Degene`, `isVertex1`,
//!   `Icurv`, `trafil1`, `Pardeb`, `V3`, ...) is [`FilDsState`] instead of a
//!   long list of locals, because the body is split across functions,
//! - `spine->Edges(i)` that OCCT indexes past the end answers `None`; that is
//!   reported as an `Err` with the OCCT call site (OCCT would dereference a
//!   null `TopoDS_Edge`),
//! - `BRep_Tool::Parameter(V, E)` is [`BRepTool::parameter_of_vertex`] and
//!   `TopExp_Explorer(E.Oriented(FORWARD), TopAbs_VERTEX)` is
//!   [`cumulated_children`] over the edge, which yields the same stored
//!   orientations.

mod faces;
mod initial;
mod nullify;
mod step;
mod tail;

use std::sync::Arc;

use occt_geom2d::curve::Curve2d;

use crate::abs::{Orientation, ShapeType};
use crate::chfi_ds::{CommonPoint, Regul, Stripe, StripeHandle, SurfData};
use crate::iterator::cumulated_children;
use crate::shape::TopoShape;
use crate::topo_tools::is_same;

use super::{spine_of, ChFi3dBuilder};

/// The locals of `ChFi3d_FilDS` that survive one `j` iteration
/// (`cxx:2572-2609`).
pub(super) struct FilDsState {
    /// `Closed` (`cxx:2575`, assigned by `cxx:2579`): `spine->IsPeriodic()`.
    pub closed: bool,
    /// `Ipoin1` (`cxx:2583`).
    pub ipo_in1: i32,
    /// `Ipoin2` (`cxx:2584`).
    pub ipo_in2: i32,
    /// `NumEdge` (`cxx:2585`), 1-based.
    pub num_edge: usize,
    /// `BoutdeVtx` (`cxx:2586`).
    pub bout_de_vtx: Option<TopoShape>,
    /// `Icurv` (`cxx:2587`).
    pub icurv: i32,
    /// `Iarc1` (`cxx:2588`).
    pub iarc1: i32,
    /// `Iarc2` (`cxx:2588`).
    pub iarc2: i32,
    /// `trafil1` (`cxx:2589`).
    pub trafil1: Orientation,
    /// `trafil2` (`cxx:2589`).
    pub trafil2: Orientation,
    /// `Pardeb` (`cxx:2591`).
    pub pardeb: f64,
    /// `Parfin` (`cxx:2591`).
    pub parfin: f64,
    /// `PCurv` (`cxx:2599`).
    pub p_curv: Option<Arc<dyn Curve2d>>,
    /// `Degene` (`cxx:2576`).
    pub degene: bool,
    /// `isVertex1` (`cxx:2576`).
    pub is_vertex1: bool,
    /// `isVertex2` (`cxx:2576`).
    pub is_vertex2: bool,
    /// `Singulier_en_Bout` (`cxx:2576`).
    pub singulier_en_bout: bool,
    /// `V3` (`cxx:2608`).
    pub v3: CommonPoint,
    /// `V4` (`cxx:2609`).
    pub v4: CommonPoint,
    /// `regcout` (`cxx:2605`), "for closed and tangent CD".
    pub regcout: Regul,
    /// `regfilfil` (`cxx:2606`), "for connections Surf/Surf".
    pub regfilfil: Regul,
}

impl FilDsState {
    /// The in-class initial values of `cxx:2572-2609`.
    fn new(closed: bool, ipo_in1: i32, ipo_in2: i32) -> Self {
        Self {
            closed,
            ipo_in1,
            ipo_in2,
            num_edge: 1,
            bout_de_vtx: None,
            icurv: 0,
            iarc1: 0,
            iarc2: 0,
            trafil1: Orientation::Forward,
            trafil2: Orientation::Forward,
            pardeb: 0.0,
            parfin: 0.0,
            p_curv: None,
            degene: false,
            is_vertex1: false,
            is_vertex2: false,
            singulier_en_bout: false,
            v3: CommonPoint::new(),
            v4: CommonPoint::new(),
            regcout: Regul::new(),
            regfilfil: Regul::new(),
        }
    }
}

/// Read `stripe` under its lock. OCCT holds the handle throughout the call;
/// the port keeps the lock only for the access, so no lock is held while the
/// DS is mutated.
pub(super) fn with_stripe<T>(
    stripe: &StripeHandle,
    f: impl FnOnce(&Stripe) -> T,
) -> Result<T, String> {
    let guard = stripe
        .read()
        .map_err(|_| "ChFi3d_FilDS: poisoned stripe".to_string())?;
    Ok(f(&guard))
}

/// `SeqFil(j)` (`cxx:2658`), `1`-based.
pub(super) fn seq_at(seq: &[Arc<SurfData>], j: usize) -> Result<Arc<SurfData>, String> {
    let index = j
        .checked_sub(1)
        .ok_or_else(|| "ChFi3d_FilDS: SeqFil is 1-based, index 0 requested".to_string())?;
    seq.get(index)
        .cloned()
        .ok_or_else(|| format!("ChFi3d_FilDS: SeqFil has no element {j}"))
}

/// `SeqFil(Index)->ChangeVertex(First, OnS)` on the stripe's shared sequence
/// entries (`cxx:2641`, `:2647`). `index` is the 1-based `SeqFil` index.
///
/// OCCT mutates the `SurfData` behind the handle; the port writes the edited
/// common point back into the sequence entry with `Arc::make_mut`.
pub(super) fn write_vertex(
    seq: &mut [Arc<SurfData>],
    index: usize,
    first: bool,
    on_s: i32,
    point: CommonPoint,
) {
    let Some(slot) = index.checked_sub(1).and_then(|i| seq.get_mut(i)) else {
        return;
    };
    *Arc::make_mut(slot).vertex_mut(first, on_s) = point;
}

/// `gp_Pnt::IsEqual(Other, 0)` of the degenerate-end tests (`cxx:2783`,
/// `cxx:3242`): exact coincidence, expressed as the `Distance <= 0`
/// comparison OCCT performs.
pub(super) fn same_point(a: &occt_core::gp::GpPnt, b: &occt_core::gp::GpPnt) -> bool {
    a.is_equal_with_tolerance(b, 0.0)
}

/// `CorDat->FirstPCurve()` / `LastPCurve()` and their `...PyCurve` siblings are
/// `occ::handle<Geom2d_Curve>`; a null handle is `None` in the port and is
/// stored as such by `ChFi3d_FilCurveInDS`.
pub(super) fn p_curve_value(
    p_curve: Option<&Arc<dyn Curve2d>>,
    par: f64,
    what: &str,
) -> Result<occt_core::gp::GpPnt2d, String> {
    let curve = p_curve
        .ok_or_else(|| format!("ChFi3d_FilDS: {what} has no p-curve on the surface"))?;
    Ok(curve.value(par))
}

impl ChFi3dBuilder {
    /// `ChFi3d_FilDS(SolidIndex, CorDat, DStr, reglist, tol3d, tol2d)`
    /// (`ChFi3d_Builder_0.cxx:2565-3341`).
    pub fn fil_ds(
        &mut self,
        solid_index: i32,
        cor_dat: &StripeHandle,
        tol3d: f64,
        tol2d: f64,
    ) -> Result<(), String> {
        // `occ::handle<ChFiDS_Spine> spine = CorDat->Spine(); bool Closed =
        // false; if (!spine.IsNull()) Closed = spine->IsPeriodic();`
        // (`cxx:2567-2574`).
        let spine = spine_of(cor_dat);
        let closed = match &spine {
            Some(sp) => sp
                .read()
                .map_err(|_| "ChFi3d_FilDS: poisoned spine".to_string())?
                .is_periodic(),
            None => false,
        };

        // `int Ipoin1 = CorDat->IndexFirstPointOnS1(); int Ipoin2 =
        // CorDat->IndexFirstPointOnS2();` (`cxx:2584-2585`).
        let (ipo_in1, ipo_in2) = with_stripe(cor_dat, |s| {
            (s.index_first_point_on_s1(), s.index_first_point_on_s2())
        })?;

        let mut st = FilDsState::new(closed, ipo_in1, ipo_in2);

        // `const NCollection_Sequence<...>& SeqFil =
        // CorDat->SetOfSurfData()->Sequence();` (`cxx:2579-2580`) over the
        // nullify prelude (`cxx:2613-2653`).
        let seq_len = with_stripe(cor_dat, |s| s.set_of_surf_data().len())?;
        if seq_len > 1 {
            let mut seq = with_stripe(cor_dat, |s| s.set_of_surf_data().to_vec())?;
            for j in 1..=seq_len {
                for on_s in 1..=2 {
                    self.nullify_degenerated_interference(&mut seq, j, on_s)?;
                }
            }
            cor_dat
                .write()
                .map_err(|_| "ChFi3d_FilDS: poisoned stripe".to_string())?
                .change_set_of_surf_data()
                .clone_from(&seq);
        }

        // `for (j = 1; j <= SeqFil.Length(); j++)` (`cxx:2655`). The sequence is
        // re-read after the nullify prelude because `SeqFil(j)` must show the
        // common-point edits it made.
        let seq_fil = with_stripe(cor_dat, |s| s.set_of_surf_data().to_vec())?;
        for j in 1..=seq_fil.len() {
            self.fil_ds_step(&mut st, &seq_fil, j, cor_dat, solid_index, tol3d, tol2d)?;
        }
        Ok(())
    }

    /// `SolidInterfs` of `cxx:2599`: the interference list of the solid that
    /// owns the fillet. `ChFi3d_Orientation(SolidInterfs, ...)` reads it and
    /// `SolidInterfs.Append(SSI)` (`cxx:2687`) writes it.
    pub(super) fn solid_interferences(&self, solid_index: i32) -> &[crate::topopebrep::ds::InterferenceHandle] {
        self.my_ds.ds().shape_interferences_at(solid_index, false)
    }

    /// `SolidInterfs.Append(SSI)` (`cxx:2687`).
    pub(super) fn push_solid_interference(
        &mut self,
        solid_index: i32,
        interference: crate::topopebrep::ds::InterferenceHandle,
    ) {
        self.my_ds
            .change_ds()
            .change_shape_interferences_at(solid_index)
            .push(interference);
    }

    /// `reglist.Append(R)` (`cxx:2936`, `:2988`, `:3088`, `:3169`, `:3231`).
    pub(super) fn append_regul(&mut self, regul: Regul) {
        self.my_regul.push(regul);
    }

    /// The `DStr` of the builder, read-only.
    pub(super) fn ds(&self) -> &crate::topopebrep::ds::DataStructure {
        self.my_ds.ds()
    }

    /// `DStr.Shape(Ishape1).Orientation()` (`cxx:2729`, `:2743`).
    pub(super) fn shape_orientation(&self, i: i32) -> Result<Orientation, String> {
        self.my_ds
            .ds()
            .shape(i, false)
            .map(|s| s.orientation())
            .ok_or_else(|| format!("ChFi3d_FilDS: the DS shape {i} is not bound"))
    }
}

/// `TopExp_Explorer(E.Oriented(FORWARD), TopAbs_VERTEX)` looking for `vertex`
/// (`cxx:2802-2809`, `:2903-2910`, `:3263-3278`): the orientation the edge
/// stores for that vertex, `FORWARD` when the edge does not carry it (the
/// initial value of `OVtx`).
pub(super) fn vertex_orientation_in_edge(
    edge: &TopoShape,
    vertex: &TopoShape,
) -> Orientation {
    for vv in cumulated_children(&edge.oriented(Orientation::Forward)) {
        if vv.shape_type() != ShapeType::Vertex {
            continue;
        }
        if is_same(&vv, vertex) {
            return vv.orientation();
        }
    }
    Orientation::Forward
}
