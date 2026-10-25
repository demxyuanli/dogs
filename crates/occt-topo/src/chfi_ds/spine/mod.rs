//! `ChFiDS_Spine` — the guideline of a fillet/chamfer: the run of edges the
//! bisector walks on, plus the elementary spines of those edges.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiDS/ChFiDS_Spine.hxx:63-289`,
//! `ChFiDS_Spine.lxx:20-170`, `ChFiDS_Spine.cxx:35-932` (OCCT 8.0.0).
//!
//! The class is not an `Adaptor3d_Curve`; it only *exposes* the adaptor
//! queries (`Resolution`, `FirstParameter`, `LastParameter`, `Value`, `D1`,
//! `D2`) over its single reused `BRepAdaptor_Curve myCurve`. OCCT rewrites
//! `indexofcurve` / `myCurve` from `const` methods through `const_cast`
//! (`cxx:548-552`, `:584-590`, `:727-734`, `:776-781`, `:826-831`); the port
//! keeps the same single adaptor but asks for `&mut self` on those methods.
//!
//! The `abscissa` array is `NCollection_HArray1<double>` in OCCT, 1-based and
//! indexed by the 1-based edge number everywhere (`cxx:238-263`, `:450-480`,
//! `:567-578`); the port keeps a `Vec<f64>` and translates the index at the
//! accessors, so the public 1-based numbering is unchanged.
//!
//! Divergences that are not heuristics:
//! - `occ::handle<ChFiDS_ElSpine>` becomes `Arc<ElSpine>`,
//! - `NCollection_List` / `NCollection_Sequence` become `Vec`,
//! - the accessors that would dereference a null `abscissa` handle, an empty
//!   `spine` sequence or a null `myCurve` answer `Option` instead of raising,
//! - `Period` / `SetFirstTgt` / `SetLastTgt` throw `Standard_Failure` in OCCT
//!   (`cxx:336-341`, `:353-358`, `:372-377`); the port answers `Option`.
//!
//! Module layout:
//! ```text
//! mod          the class, `SpineClass`, the edge / elementary-spine lists, the
//!              parameters, `Index`, `Index(edge)`, and the private `abscissa`
//!              helpers
//! state        the extremity targets and reference, plus the status flags
//! curve_query  `Absc`, `Parameter`, `Value`, `D0`-`D2`, `Line`, `Circle`
//! ```

use std::sync::{Arc, RwLock};

use occt_core::elib::clib;
use occt_core::gp::{GpAx2, GpCirc, GpDir, GpLin, GpPnt, GpVec};
use occt_core::precision::CONFUSION;
use occt_geom::curve::Curve;
use occt_geom::gcpnts::{abscissa_point_with_guess, curve_length, curve_length_between};

use crate::geom_abs_curve_type::GeomAbsCurveType;
use crate::shape::{Edge, TopoShape};
use crate::topo_tools;
use crate::topopebrep::tool::adaptor::CurveAdaptor;

use super::chamf_spine::ChamfSpineData;
use super::el_spine::ElSpine;
use super::enums::{ChamfMode, ErrorStatus, State, TypeOfConcavity};
use super::fil_spine::FilSpineData;

mod curve_query;
mod state;

/// Which of the three `ChFiDS` guideline classes an object is: OCCT derives
/// `ChFiDS_FilSpine` and `ChFiDS_ChamfSpine` from `ChFiDS_Spine` and tells them
/// apart by RTTI (`IMPLEMENT_STANDARD_RTTIEXT`), and `ChFi3d` down-casts
/// (`ChFi3d_FilBuilder::PerformTwoCorner` etc.).
///
/// A `SpineHandle` is `Arc<RwLock<Spine>>`, so the handle cannot hold a
/// subclass: the subclass members live in [`Spine`] behind this discriminant
/// instead of behind a vtable.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpineClass {
    /// `ChFiDS_Spine`
    Base,
    /// `ChFiDS_FilSpine`
    Fil,
    /// `ChFiDS_ChamfSpine`
    Chamf,
}


/// `occ::handle<ChFiDS_Spine>`: OCCT hands the spine out by shared handle and
/// mutates it through it (`ChFiDS_Stripe::ChangeSpine`), so the port uses the
/// same shared-mutable handle as the other transient data of this crate
/// (`topopebrep/ds/interference.rs:33`).
pub type SpineHandle = Arc<RwLock<Spine>>;

/// Source: `ChFiDS_Spine.hxx:68-289`.
#[derive(Debug, Clone)]
pub struct Spine {
    /// `splitdone`, read and written by `super::fil_spine` too.
    pub(crate) splitdone: bool,
    /// `elspines`
    elspines: Vec<Arc<ElSpine>>,
    /// `offset_elspines`
    offset_elspines: Vec<Arc<ElSpine>>,
    /// `myMode`
    my_mode: ChamfMode,
    /// `myCurve` — `BRepAdaptor_Curve`, default-constructed (null edge) until
    /// `Load` / `SetCurrent`.
    my_curve: Option<CurveAdaptor>,
    /// `myOffsetCurve` — never initialized by this class in OCCT; kept for the
    /// member list.
    my_offset_curve: Option<CurveAdaptor>,
    /// `indexofcurve`
    index_of_curve: i32,
    /// `myTypeOfConcavity`
    my_type_of_concavity: TypeOfConcavity,
    /// `firstState`
    first_state: State,
    /// `lastState`
    last_state: State,
    /// `spine` — `NCollection_Sequence<TopoDS_Shape>` holding edges.
    spine: Vec<TopoShape>,
    /// `offsetspine`
    offset_spine: Vec<TopoShape>,
    /// `abscissa` — cumulative arc length per edge, 1-based in OCCT.
    abscissa: Option<Vec<f64>>,
    /// `offset_abscissa` — declared in OCCT (`hxx:270`) and never written by
    /// this class; kept for the member list.
    offset_abscissa: Option<Vec<f64>>,
    /// `tolesp`
    tolesp: f64,
    /// `firstparam`
    firstparam: f64,
    /// `lastparam`
    lastparam: f64,
    /// `firstprolon`
    firstprolon: bool,
    /// `lastprolon`
    lastprolon: bool,
    /// `firstistgt`
    firstistgt: bool,
    /// `lastistgt`
    lastistgt: bool,
    /// `firsttgtpar`
    firsttgtpar: f64,
    /// `lasttgtpar`
    lasttgtpar: f64,
    /// `hasfirsttgt`
    hasfirsttgt: bool,
    /// `haslasttgt`
    haslasttgt: bool,
    /// `firstori`
    firstori: GpPnt,
    /// `lastori`
    lastori: GpPnt,
    /// `firsttgt`
    firsttgt: GpVec,
    /// `lasttgt`
    lasttgt: GpVec,
    /// `valref`
    valref: f64,
    /// `hasref`
    hasref: bool,
    /// `errorstate`
    errorstate: ErrorStatus,
    /// Which subclass this object is (`SpineClass`).
    class: SpineClass,
    /// `ChFiDS_FilSpine::parandrad` / `::laws` (`ChFiDS_FilSpine.hxx:98-99`).
    /// OCCT keeps them in the derived class; a `RwLock<Spine>` cannot hold a
    /// derived class, so they live here (see [`SpineClass`]).
    pub(crate) fil: FilSpineData,
    /// `ChFiDS_ChamfSpine::d1` / `::d2` / `::angle` / `::mChamf`
    /// (`ChFiDS_ChamfSpine.hxx:53-56`), for the same reason.
    pub(crate) chamf: ChamfSpineData,
}

impl Default for Spine {
    fn default() -> Self {
        Self::new()
    }
}

impl Spine {
    /// `ChFiDS_Spine()` (`cxx:35-59`): `tolesp = Precision::Confusion()`, both
    /// ends `ChFiDS_OnSame`, concavity `ChFiDS_Other`, mode
    /// `ChFiDS_ClassicChamfer`.
    pub fn new() -> Self {
        Self::with_tolerance(CONFUSION)
    }

    /// `ChFiDS_Spine(const double Tol)` (`cxx:63-86`).
    pub fn with_tolerance(tol: f64) -> Self {
        Self {
            splitdone: false,
            elspines: Vec::new(),
            offset_elspines: Vec::new(),
            my_mode: ChamfMode::ClassicChamfer,
            my_curve: None,
            my_offset_curve: None,
            index_of_curve: 0,
            my_type_of_concavity: TypeOfConcavity::Other,
            first_state: State::OnSame,
            last_state: State::OnSame,
            spine: Vec::new(),
            offset_spine: Vec::new(),
            abscissa: None,
            offset_abscissa: None,
            tolesp: tol,
            firstparam: 0.0,
            lastparam: 0.0,
            firstprolon: false,
            lastprolon: false,
            firstistgt: false,
            lastistgt: false,
            firsttgtpar: 0.0,
            lasttgtpar: 0.0,
            hasfirsttgt: false,
            haslasttgt: false,
            firstori: GpPnt::zero(),
            lastori: GpPnt::zero(),
            firsttgt: GpVec::zero(),
            lasttgt: GpVec::zero(),
            valref: 0.0,
            hasref: false,
            errorstate: ErrorStatus::Ok,
            class: SpineClass::Base,
            fil: FilSpineData::new(),
            chamf: ChamfSpineData::new(),
        }
    }

    /// `ChFiDS_FilSpine()` (`ChFiDS_FilSpine.cxx:41`): `ChFiDS_Spine()` with
    /// the fillet members default-constructed.
    pub fn new_fil_spine() -> Self {
        Self::with_tolerance_fil(CONFUSION)
    }

    /// `ChFiDS_FilSpine(const double Tol)` (`ChFiDS_FilSpine.cxx:43-46`).
    pub fn with_tolerance_fil(tol: f64) -> Self {
        let mut spine = Self::with_tolerance(tol);
        spine.class = SpineClass::Fil;
        spine
    }

    /// `ChFiDS_ChamfSpine()` (`ChFiDS_ChamfSpine.cxx:27-33`): `ChFiDS_Spine()`,
    /// then `d1 = d2 = angle = 0`, `mChamf = ChFiDS_Sym` and
    /// `myMode = ChFiDS_ClassicChamfer` (the base constructor already writes
    /// that mode, `ChFiDS_Spine.cxx:57`; OCCT writes it twice).
    pub fn new_chamf_spine() -> Self {
        Self::with_tolerance_chamf(CONFUSION)
    }

    /// `ChFiDS_ChamfSpine(const double Tol)` (`ChFiDS_ChamfSpine.cxx:35-42`).
    pub fn with_tolerance_chamf(tol: f64) -> Self {
        let mut spine = Self::with_tolerance(tol);
        spine.class = SpineClass::Chamf;
        spine.my_mode = ChamfMode::ClassicChamfer;
        spine.chamf = ChamfSpineData::new();
        spine
    }

    /// The RTTI discriminant ([`SpineClass`]).
    pub fn class(&self) -> SpineClass {
        self.class
    }

    /// `ChFiDS_FilSpine`'s own members (`ChFiDS_FilSpine.hxx:98-99`).
    pub fn fil(&self) -> &FilSpineData {
        &self.fil
    }

    /// `ChFiDS_FilSpine`'s own members, for the writers in
    /// [`super::fil_spine`].
    pub fn fil_mut(&mut self) -> &mut FilSpineData {
        &mut self.fil
    }

    /// `ChFiDS_ChamfSpine`'s own members (`ChFiDS_ChamfSpine.hxx:53-56`).
    pub fn chamf(&self) -> &ChamfSpineData {
        &self.chamf
    }

    /// `ChFiDS_ChamfSpine`'s own members, for the writers in
    /// [`super::chamf_spine`].
    pub fn chamf_mut(&mut self) -> &mut ChamfSpineData {
        &mut self.chamf
    }

    /// `SetEdges(E)` (`lxx:124-127`): appends to the guideline.
    pub fn set_edges(&mut self, e: &TopoShape) {
        self.spine.push(e.clone());
    }

    /// `SetOffsetEdges(E)` (`lxx:131-134`).
    pub fn set_offset_edges(&mut self, e: &TopoShape) {
        self.offset_spine.push(e.clone());
    }

    /// `PutInFirst(E)` (`lxx:138-141`) — `InsertBefore(1, E)`.
    pub fn put_in_first(&mut self, e: &TopoShape) {
        self.spine.insert(0, e.clone());
    }

    /// `PutInFirstOffset(E)` (`lxx:144-147`).
    pub fn put_in_first_offset(&mut self, e: &TopoShape) {
        self.offset_spine.insert(0, e.clone());
    }

    /// `NbEdges()` (`lxx:104-107`).
    pub fn nb_edges(&self) -> usize {
        self.spine.len()
    }

    /// `Edges(I)` (`lxx:110-113`), `I` 1-based.
    pub fn edges(&self, i: usize) -> Option<&TopoShape> {
        if i == 0 {
            return None;
        }
        self.spine.get(i - 1)
    }

    /// `OffsetEdges(I)` (`lxx:116-119`), `I` 1-based.
    pub fn offset_edges(&self, i: usize) -> Option<&TopoShape> {
        if i == 0 {
            return None;
        }
        self.offset_spine.get(i - 1)
    }

    /// `SetFirstStatus(S)` (`lxx:27-30`).
    pub fn set_first_status(&mut self, s: State) {
        self.first_state = s;
    }

    /// `SetLastStatus(S)` (`lxx:33-36`).
    pub fn set_last_status(&mut self, s: State) {
        self.last_state = s;
    }

    /// `FirstStatus()` (`lxx:49-52`).
    pub fn first_status(&self) -> State {
        self.first_state
    }

    /// `LastStatus()` (`lxx:56-59`).
    pub fn last_status(&self) -> State {
        self.last_state
    }

    /// `SetStatus(S, IsFirst)` (`lxx:63-71`).
    pub fn set_status(&mut self, s: State, is_first: bool) {
        if is_first {
            self.first_state = s;
        } else {
            self.last_state = s;
        }
    }

    /// `Status(IsFirst)` (`lxx:76-82`).
    pub fn status(&self, is_first: bool) -> State {
        if is_first {
            self.first_state
        } else {
            self.last_state
        }
    }

    /// `GetTypeOfConcavity()` (`lxx:35-38`).
    pub fn type_of_concavity(&self) -> TypeOfConcavity {
        self.my_type_of_concavity
    }

    /// `SetTypeOfConcavity(theType)` (`lxx:21-24`).
    pub fn set_type_of_concavity(&mut self, the_type: TypeOfConcavity) {
        self.my_type_of_concavity = the_type;
    }

    /// `IsTangencyExtremity(IsFirst)` (`lxx:88-96`).
    pub fn is_tangency_extremity(&self, is_first: bool) -> bool {
        if is_first {
            self.firstistgt
        } else {
            self.lastistgt
        }
    }

    /// `SetTangencyExtremity(IsTangency, IsFirst)` (`lxx:83-89`).
    pub fn set_tangency_extremity(&mut self, is_tangency: bool, is_first: bool) {
        if is_first {
            self.firstistgt = is_tangency;
        } else {
            self.lastistgt = is_tangency;
        }
    }

    /// `AppendElSpine(Els)` (`cxx:89-92`); virtual in OCCT.
    pub fn append_el_spine(&mut self, els: Arc<ElSpine>) {
        self.elspines.push(els);
    }

    /// `AppendOffsetElSpine(Els)` (`cxx:96-99`); virtual in OCCT.
    pub fn append_offset_el_spine(&mut self, els: Arc<ElSpine>) {
        self.offset_elspines.push(els);
    }

    /// `ElSpine(IE)` (`cxx:108-116`): the mid-abscissa elementary spine of
    /// edge `IE`.
    pub fn el_spine_of_edge_index(&self, ie: i32) -> Option<Arc<ElSpine>> {
        let mut wmil = 0.5 * (self.first_parameter_of(ie)? + self.last_parameter_of(ie)?);
        if self.is_periodic() {
            if let Some(last) = self.last_parameter() {
                wmil = clib::in_period(wmil, self.first_parameter(), last);
            }
        }
        self.el_spine(wmil)
    }

    /// `ElSpine(const TopoDS_Edge&)` (`cxx:104-107`).
    pub fn el_spine_of_edge(&self, e: &TopoShape) -> Option<Arc<ElSpine>> {
        self.el_spine_of_edge_index(self.index_of_edge(e))
    }

    /// `ElSpine(W)` (`cxx:118-137`): the single elementary spine when there is
    /// one, otherwise the first whose `[FirstParameter, LastParameter]` covers
    /// `W`. OCCT returns a null handle when none does.
    pub fn el_spine(&self, w: f64) -> Option<Arc<ElSpine>> {
        if self.elspines.len() == 1 {
            return self.elspines.first().cloned();
        }
        for cur in &self.elspines {
            let uf = cur.first_parameter();
            let ul = cur.last_parameter();
            if uf <= w && w <= ul {
                return Some(cur.clone());
            }
        }
        None
    }

    /// `elspines`, for `ChFiDS_FilSpine::Law` which walks the list in step with
    /// `laws` (`ChFiDS_FilSpine.cxx:860-869`).
    pub fn el_spines(&self) -> &[Arc<ElSpine>] {
        &self.elspines
    }

    /// `ChangeElSpines()` (`cxx:144-147`).
    pub fn el_spines_mut(&mut self) -> &mut Vec<Arc<ElSpine>> {
        &mut self.elspines
    }

    /// `ChangeOffsetElSpines()` (`cxx:151-154`).
    pub fn offset_el_spines_mut(&mut self) -> &mut Vec<Arc<ElSpine>> {
        &mut self.offset_elspines
    }

    /// `Reset(AllData = false)` (`cxx:173-183`); virtual in OCCT.
    pub fn reset(&mut self, all_data: bool) -> Option<()> {
        self.splitdone = false;
        self.elspines.clear();
        if all_data {
            self.firstparam = 0.0;
            self.lastparam = self.abscissa_upper()?;
            self.firstprolon = false;
            self.lastprolon = false;
        }
        Some(())
    }

    /// `SplitDone()` (`cxx:165-168`).
    pub fn split_done(&self) -> bool {
        self.splitdone
    }

    /// `SplitDone(B)` (`cxx:157-160`).
    pub fn set_split_done(&mut self, b: bool) {
        self.splitdone = b;
    }

    /// `Load()` (`cxx:504-534`): one `BRepAdaptor_Curve::Initialize` per edge,
    /// cumulative `GCPnts_AbscissaPoint::Length`, then the `tolesp` clamp to
    /// `5e-5 * (umax - umin)`.
    pub fn load(&mut self) -> Option<()> {
        let len = self.spine.len();
        let mut abscissa = vec![0.0f64; len];
        let mut a1 = 0.0f64;
        for i in 0..len {
            // `cxx:516-518`: `myCurve.Initialize(spine.Value(i))` then
            // `a1 += GCPnts_AbscissaPoint::Length(myCurve)`.
            let adaptor = CurveAdaptor::of_edge(&Edge(self.spine[i].clone()));
            a1 += curve_length(&adaptor);
            abscissa[i] = a1;
        }
        self.abscissa = Some(abscissa);
        self.index_of_curve = 1;
        let first = self.spine.first()?.clone();
        self.my_curve = Some(CurveAdaptor::of_edge(&Edge(first)));

        // `cxx:523-532`.
        let umin = self.first_parameter();
        let umax = self.last_parameter()?;
        let new_tolesp = 5.0e-5 * (umax - umin);
        if self.tolesp > new_tolesp {
            self.tolesp = new_tolesp;
        }
        Some(())
    }

    /// `Resolution(R3d)` (`cxx:346-349`): the guideline answers `R3d`.
    pub fn resolution(&self, r3d: f64) -> f64 {
        r3d
    }

    /// `GetTolesp()` (`lxx:167-170`).
    pub fn tolesp(&self) -> f64 {
        self.tolesp
    }

    /// `IsClosed()` (`cxx:275-278`): first vertex `IsSame` last vertex. OCCT
    /// raises on an empty `spine`; the port answers `None`.
    pub fn is_closed(&self) -> Option<bool> {
        let first = self.first_vertex()?;
        let last = self.last_vertex()?;
        Some(topo_tools::is_same(&first, &last))
    }

    /// `FirstParameter()` (`cxx:186-193`): `firstparam` only when the guideline
    /// was extended at the start, `0.` otherwise.
    pub fn first_parameter(&self) -> f64 {
        if self.firstprolon {
            return self.firstparam;
        }
        0.0
    }

    /// `LastParameter()` (`cxx:197-204`).
    pub fn last_parameter(&self) -> Option<f64> {
        if self.lastprolon {
            return Some(self.lastparam);
        }
        self.abscissa_upper()
    }

    /// `SetFirstParameter(Par)` (`cxx:208-218`), `OCCT_DEBUG` block dropped.
    pub fn set_first_parameter(&mut self, par: f64) {
        self.firstprolon = true;
        self.firstparam = par;
    }

    /// `SetLastParameter(Par)` (`cxx:222-233`), `OCCT_DEBUG` block dropped.
    pub fn set_last_parameter(&mut self, par: f64) {
        self.lastprolon = true;
        self.lastparam = par;
    }

    /// `FirstParameter(IndexSpine)` (`cxx:237-244`): the total length of the
    /// arcs before `IndexSpine`, `0.` on the first one, `IndexSpine` 1-based.
    pub fn first_parameter_of(&self, index_spine: i32) -> Option<f64> {
        if index_spine == 1 {
            return Some(0.0);
        }
        self.abscissa_value(index_spine - 1)
    }

    /// `LastParameter(IndexSpine)` (`cxx:248-251`): the total length up to and
    /// including `IndexSpine`.
    pub fn last_parameter_of(&self, index_spine: i32) -> Option<f64> {
        self.abscissa_value(index_spine)
    }

    /// `Length(IndexSpine)` (`cxx:255-263`).
    pub fn length_of(&self, index_spine: i32) -> Option<f64> {
        if index_spine == 1 {
            return self.abscissa_value(index_spine);
        }
        Some(self.abscissa_value(index_spine)? - self.abscissa_value(index_spine - 1)?)
    }

    /// `IsPeriodic()` (`cxx:267-270`): `firstState == ChFiDS_Closed`.
    pub fn is_periodic(&self) -> bool {
        self.first_state == State::Closed
    }

    /// `Period()` (`cxx:335-342`): throws `Standard_Failure("Non-periodic
    /// Spine")` when the guideline is not closed; the port answers `None`.
    pub fn period(&self) -> Option<f64> {
        if !self.is_periodic() {
            return None;
        }
        self.abscissa_upper()
    }

    /// `FirstVertex()` (`cxx:283-290`): the stored `FORWARD` vertex of the
    /// first edge when that edge is `FORWARD`, else its stored `REVERSED` one
    /// (`TopExp::FirstVertex` / `LastVertex` with `CumOri = false`).
    pub fn first_vertex(&self) -> Option<TopoShape> {
        let e = self.spine.first()?;
        let (first, last) = topo_tools::edge_vertices(&Edge(e.clone()));
        let v = if e.orientation().is_forward() {
            first
        } else {
            last
        };
        v.map(|v| v.0)
    }

    /// `LastVertex()` (`cxx:294-301`).
    pub fn last_vertex(&self) -> Option<TopoShape> {
        let e = self.spine.last()?;
        let (first, last) = topo_tools::edge_vertices(&Edge(e.clone()));
        let v = if e.orientation().is_forward() {
            last
        } else {
            first
        };
        v.map(|v| v.0)
    }

    /// `Absc(const TopoDS_Vertex&)` (`cxx:304-330`): the abscissa of a vertex of
    /// the guideline, `-1.` when the vertex is not on it (the same sentinel as
    /// OCCT, see `ChFi3d_Builder_1.cxx:1315-1322`).
    pub fn absc_of_vertex(&self, v: &TopoShape) -> f64 {
        for i in 1..=self.spine.len() as i32 {
            let e = &self.spine[(i - 1) as usize];
            let (d, f) = topo_tools::edge_vertices(&Edge(e.clone()));
            let d = d.map(|v| v.0);
            let f = f.map(|v| v.0);
            let is_d = d.as_ref().is_some_and(|d| topo_tools::is_same(d, v));
            let is_f = f.as_ref().is_some_and(|f| topo_tools::is_same(f, v));
            let forward = e.orientation().is_forward();
            let par = if is_d {
                if forward {
                    self.first_parameter_of(i)
                } else {
                    self.last_parameter_of(i)
                }
            } else if is_f {
                if forward {
                    self.last_parameter_of(i)
                } else {
                    self.first_parameter_of(i)
                }
            } else {
                continue;
            };
            return par.unwrap_or(-1.0);
        }
        -1.0
    }

    /// `Index(W, Forward = true)` (`cxx:443-478`), 1-based.
    pub fn index(&self, w: f64, forward: bool) -> Option<i32> {
        let len = self.abscissa_len()?;
        let mut par = w;
        let last = self.abscissa_upper()?;
        let t = self.tolesp.max(CONFUSION);

        if self.is_periodic() && par.abs() >= t && (par - last).abs() >= t {
            par = clib::in_period(par, 0.0, last);
        }

        let mut f = 0.0f64;
        let mut l = 0.0f64;
        let mut ind = 1i32;
        while ind <= len {
            f = l;
            l = self.abscissa_value(ind)?;
            if par < l || ind == len {
                break;
            }
            ind += 1;
        }
        if forward && ind < len && (par - l).abs() < t {
            ind += 1;
        } else if !forward && ind > 1 && (par - f).abs() < t {
            ind -= 1;
        } else if forward && self.is_periodic() && ind == len && (par - l).abs() < t {
            ind = 1;
        } else if !forward && self.is_periodic() && ind == 1 && (par - f).abs() < t {
            ind = len;
        }
        Some(ind)
    }

    /// `Index(const TopoDS_Edge&)` (`cxx:483-493`): `0` when the edge is not
    /// part of the guideline (OCCT's sentinel).
    pub fn index_of_edge(&self, e: &TopoShape) -> i32 {
        for (i, s) in self.spine.iter().enumerate() {
            if topo_tools::is_same(s, e) {
                return i as i32 + 1;
            }
        }
        0
    }


    /// `Prepare(L, Ind)` (`cxx:618-699`), rewritten to return `Ind` and to take
    /// `L` by reference, exactly like the OCCT in/out parameters.
    fn prepare(&self, l: &mut f64) -> Option<i32> {
        let tol = self.tolesp.max(CONFUSION);
        let last = self.abscissa_upper()?;
        let len = self.abscissa_len()?;

        // `cxx:625-628`.
        if self.is_periodic() && l.abs() >= tol && (*l - last).abs() >= tol {
            *l = clib::in_period(*l, 0.0, last);
        }

        let mut ind;
        if self.hasfirsttgt && (*l <= self.firsttgtpar) {
            // `cxx:629-641`.
            if self.hasref && self.valref >= *l && (*l - self.firsttgtpar).abs() <= tol {
                ind = self.index(*l, true)?;
            } else {
                ind = -1;
                *l -= self.firsttgtpar;
            }
        } else if *l <= 0.0 {
            // `cxx:642-645`.
            ind = 1;
        } else if self.haslasttgt && (*l >= self.lasttgtpar) {
            // `cxx:646-656`.
            if self.hasref && self.valref <= *l && (*l - self.lasttgtpar).abs() <= tol {
                ind = self.index(*l, true)?;
            } else {
                ind = len + 1;
                *l -= self.lasttgtpar;
            }
        } else if *l >= last {
            // `cxx:657-660`.
            ind = len;
        } else {
            // `cxx:661-687`.
            ind = len;
            for i in 1..len {
                if *l < self.abscissa_value(i)? {
                    ind = i;
                    break;
                }
            }
            if self.hasref {
                if *l >= self.valref && ind != 1 {
                    if (*l - self.abscissa_value(ind - 1)?).abs() <= CONFUSION {
                        ind -= 1;
                    }
                } else if *l <= self.valref && ind != len {
                    if (*l - self.abscissa_value(ind)?).abs() <= CONFUSION {
                        ind += 1;
                    }
                }
            }
        }

        // `cxx:688-698`: local abscissa inside the elementary spine, measured
        // from the start the orientation walks from.
        if ind >= 1 && ind <= len {
            if self
                .spine
                .get(usize::try_from(ind - 1).ok()?)?
                .orientation()
                .is_reversed()
            {
                *l = self.abscissa_value(ind)? - *l;
            } else if ind != 1 {
                *l -= self.abscissa_value(ind - 1)?;
            }
        }
        Some(ind)
    }

    /// `indexofcurve = Index; myCurve.Initialize(spine.Value(Index))`
    /// (`cxx:856-862` and the five `const_cast` copies of it).
    fn ensure_current(&mut self, index: i32) -> Option<()> {
        if index != self.index_of_curve {
            self.index_of_curve = index;
            let edge = self.spine.get(usize::try_from(index - 1).ok()?)?.clone();
            self.my_curve = Some(CurveAdaptor::of_edge(&Edge(edge)));
        }
        self.my_curve.as_ref()?;
        Some(())
    }

    /// `abscissa->Value(I)`, `I` 1-based.
    fn abscissa_value(&self, i: i32) -> Option<f64> {
        let idx = usize::try_from(i).ok()?;
        if idx == 0 {
            return None;
        }
        self.abscissa.as_ref()?.get(idx - 1).copied()
    }

    /// `abscissa->Value(abscissa->Upper())`.
    fn abscissa_upper(&self) -> Option<f64> {
        self.abscissa.as_ref()?.last().copied()
    }

    /// `abscissa->Length()`.
    fn abscissa_len(&self) -> Option<i32> {
        i32::try_from(self.abscissa.as_ref()?.len()).ok()
    }
}
