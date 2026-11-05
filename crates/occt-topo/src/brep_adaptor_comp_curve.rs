//! `BRepAdaptor_CompCurve` — an `Adaptor3d_Curve` over a whole wire.
//!
//! Source: `ModelingData/TKBRep/BRepAdaptor/BRepAdaptor_CompCurve.cxx:39-600`
//! (OCCT 8.0.0), `.hxx:30-145`.
//!
//! The adaptor walks the wire with `BRepTools_WireExplorer`, keeps one
//! `BRepAdaptor_Curve` per non-degenerated edge plus a knot table, and maps the
//! global parameter `W` onto the current edge's own parameter with `Prepare`
//! (`cxx:497-568`). Two knot conventions exist:
//!
//! - `theIsAC == false` (the `BRepFill_SectionPlacement` constructor call): the
//!   knots are the edge ranks `0, 1, 2, ...`, so the global parameter is the
//!   "edge index" space of `BRepFill_LocationLaw`.
//! - `theIsAC == true`: `myKnots(i+1) = myKnots(i) + GCPnts_AbscissaPoint::Length`
//!   of that edge's adaptor (`cxx:139-140`).
//!
//! Translation notes:
//! - `myCurves` (`NCollection_HArray1<BRepAdaptor_Curve>`) and `myKnots`
//!   (`NCollection_HArray1<double>`, `1..NbEdge+1`) are the paired vectors
//!   `curves` / `edges` (index `i-1` is OCCT's `i`) and `knots` (same
//!   `i-1` shift; `knots[i]` is OCCT's `myKnots(i+1)`).
//! - `BRepAdaptor_Curve` is [`CurveAdaptor`]; its `Edge()` is the `edges`
//!   entry, because OCCT reads the orientation and the `BRep_Tool::Range` of
//!   the edge the adaptor was built from.
//! - `GCPnts_AbscissaPoint::Length(adaptor)` is
//!   [`occt_geom::gcpnts::curve_length`] over that adaptor.
//! - `TopExp::CommonVertex` / `LastVertex` are the existing
//!   [`crate::topo_tools::common_vertex`] and
//!   [`crate::brep_tools_wire_explorer::last_vertex_cum_ori`] (`CumOri = true`).
//! - `Adaptor3d_Curve::EvalD0..EvalDN` are the `Curve` trait methods; the
//!   chain rule `d^n` scaling of `cxx:384-407` is applied around them.
//! - `myWire.IsNull()` (`cxx:107` region) has no counterpart: `Wire` is never
//!   null here, so an empty wire simply leaves `curves` empty and the
//!   geometric queries raise the `Standard_NullObject` OCCT would hit on a
//!   null `myCurves` handle.

use occt_core::gp::{GpPnt, GpVec};
use occt_core::kernel::geomabs::CurveType;
use occt_geom::curve::Curve;
use occt_geom::gcpnts::curve_length;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::{last_vertex_cum_ori, WireExplorer};
use crate::shape::{Edge, Wire};
use crate::topo_tools::{common_vertex, is_same};
use crate::topopebrep::tool::adaptor::CurveAdaptor;

/// `BRepAdaptor_CompCurve` (`BRepAdaptor_CompCurve.hxx:30-145`).
#[derive(Clone)]
pub struct CompCurveAdaptor {
    /// `myWire` (`hxx:126`).
    wire: Option<Wire>,
    /// `TFirst` (`hxx:118`).
    t_first: f64,
    /// `TLast` (`hxx:119`).
    t_last: f64,
    /// `PTol` (`hxx:120`).
    p_tol: f64,
    /// `CurIndex` (`hxx:121`), 1-based; `-1` before `Initialize`.
    cur_index: i32,
    /// `Forward` (`hxx:122`).
    forward: bool,
    /// `IsbyAC` (`hxx:123`).
    is_by_ac: bool,
    /// `myCurves` (`hxx:124`).
    curves: Vec<CurveAdaptor>,
    /// `myCurves->Value(i).Edge()`.
    edges: Vec<Edge>,
    /// `myKnots` (`hxx:125`), index `i-1` is OCCT's `myKnots(i)`.
    knots: Vec<f64>,
}

impl Default for CompCurveAdaptor {
    /// `BRepAdaptor_CompCurve()` (`cxx:39-47`).
    fn default() -> Self {
        Self {
            wire: None,
            t_first: 0.0,
            t_last: 0.0,
            p_tol: 0.0,
            cur_index: -1,
            forward: false,
            is_by_ac: false,
            curves: Vec::new(),
            edges: Vec::new(),
            knots: Vec::new(),
        }
    }
}

impl CompCurveAdaptor {
    /// `BRepAdaptor_CompCurve()` (`cxx:39-47`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `BRepAdaptor_CompCurve(W, theIsAC = false)` (`cxx:49-58`).
    pub fn of_wire(w: &Wire, is_ac: bool) -> Self {
        let mut adaptor = Self {
            wire: Some(w.clone()),
            is_by_ac: is_ac,
            ..Self::default()
        };
        adaptor.initialize(w, is_ac);
        adaptor
    }

    /// `BRepAdaptor_CompCurve(W, theIsAC, First, Last, Tol)` (`cxx:61-75`).
    pub fn with_range(w: &Wire, is_ac: bool, first: f64, last: f64, tol: f64) -> Self {
        let mut adaptor = Self {
            wire: Some(w.clone()),
            t_first: first,
            t_last: last,
            p_tol: tol,
            is_by_ac: is_ac,
            ..Self::default()
        };
        adaptor.initialize_trimmed(w, is_ac, first, last, tol);
        adaptor
    }

    /// `ShallowCopy` (`cxx:79-100`): every member is copied, the per-edge
    /// adaptors are shallow copies of themselves (`CurveAdaptor` is a value
    /// type here, so a clone *is* the shallow copy).
    pub fn shallow_copy(&self) -> Self {
        self.clone()
    }

    /// `Initialize(W, AC)` (`cxx:102-172`).
    pub fn initialize(&mut self, w: &Wire, ac: bool) {
        self.wire = Some(w.clone());
        self.p_tol = 0.0;
        self.is_by_ac = ac;

        // `cxx:111-118`: count then materialize the non-degenerated edges.
        let mut edges: Vec<Edge> = Vec::new();
        let mut explorer = WireExplorer::of_wire(w);
        while explorer.more() {
            if let Some(e) = explorer.current_edge() {
                if !BRepTool::is_degenerated(&e) {
                    edges.push(e);
                }
            }
            explorer.next();
        }

        // `cxx:120-123`: `NbEdge == 0` returns with `myCurves` still null.
        if edges.is_empty() {
            return;
        }

        let nb_edge = edges.len();
        // `cxx:127`: `CurIndex = (NbEdge + 1) / 2`.
        self.cur_index = ((nb_edge + 1) / 2) as i32;
        self.knots = vec![0.0; nb_edge + 1];
        self.curves.clear();
        self.edges.clear();

        // `cxx:130-148`.
        for (ii, e) in edges.iter().enumerate() {
            let adaptor = CurveAdaptor::of_edge(e);
            if ac {
                // `myKnots->ChangeValue(ii + 1) += GCPnts_AbscissaPoint::Length(...)`.
                let length = curve_length(&adaptor);
                self.knots[ii + 1] = self.knots[ii] + length;
            } else {
                self.knots[ii + 1] = (ii + 1) as f64;
            }
            self.curves.push(adaptor);
            self.edges.push(e.clone());
        }

        // `cxx:149-170`: the direction of parsing.
        self.forward = true;
        if nb_edge > 2 || (nb_edge == 2 && !w.0.is_closed()) {
            let or = self.edges[0].0.orientation();
            let vi = common_vertex(&self.edges[0], &self.edges[1]);
            let vl = last_vertex_cum_ori(&self.edges[0]).map(|v| v.0);
            let common_is_last = match (&vi, &vl) {
                (Some(a), Some(b)) => is_same(a, b),
                _ => false,
            };
            if common_is_last {
                // `cxx:155-163`: the direction of parsing is preserved.
                if or == Orientation::Reversed {
                    self.forward = false;
                }
            } else {
                // `cxx:165-169`: the direction of parsing is reversed.
                if or != Orientation::Reversed {
                    self.forward = false;
                }
            }
        }

        // `cxx:171-172`.
        self.t_first = 0.0;
        self.t_last = self.knots[nb_edge];
    }

    /// `Initialize(W, AC, First, Last, Tol)` (`cxx:174-235`).
    pub fn initialize_trimmed(&mut self, w: &Wire, ac: bool, first: f64, last: f64, tol: f64) {
        self.initialize(w, ac);
        self.t_first = first;
        self.t_last = last;
        self.p_tol = tol;
        if self.curves.is_empty() {
            return;
        }

        // `cxx:187-193`.
        let mut i1 = self.cur_index;
        let mut i2 = self.cur_index;
        let mut f = self.t_first;
        let mut d = 0.0;
        self.prepare(&mut f, &mut d, &mut i1);
        let mut l = self.t_last;
        self.prepare(&mut l, &mut d, &mut i2);
        // `cxx:194`: small optimization.
        self.cur_index = (i1 + i2) / 2;

        if i1 == i2 {
            // `cxx:196-205`.
            let (a, b) = if l > f { (f, l) } else { (l, f) };
            let index = (i1 - 1) as usize;
            self.curves[index] = self.curves[index].trim(a, b);
        } else {
            // `cxx:207-235`.
            let index1 = (i1 - 1) as usize;
            let k = self.curves[index1].last_parameter();
            let (a, b) = if k > f { (f, k) } else { (k, f) };
            self.curves[index1] = self.curves[index1].trim(a, b);

            let index2 = (i2 - 1) as usize;
            let k = self.curves[index2].first_parameter();
            let (a, b) = if k <= l { (k, l) } else { (l, k) };
            self.curves[index2] = self.curves[index2].trim(a, b);
        }
    }

    /// `Wire()` (`cxx:237-240`).
    pub fn wire(&self) -> Option<&Wire> {
        self.wire.as_ref()
    }

    /// `Edge(U, E, UonE)` (`cxx:242-249`): the edge holding `U` and the local
    /// parameter on it.
    pub fn edge_at(&self, u: f64) -> Option<(Edge, f64)> {
        let (index, local, _) = self.locate(u)?;
        Some((self.edges[index].clone(), local))
    }

    /// `FirstParameter()` (`cxx:251-254`).
    pub fn first_parameter(&self) -> f64 {
        self.t_first
    }

    /// `LastParameter()` (`cxx:256-259`).
    pub fn last_parameter(&self) -> f64 {
        self.t_last
    }

    /// `Continuity()` (`cxx:261-268`).
    pub fn continuity(&self) -> u8 {
        if self.curves.len() > 1 {
            // `GeomAbs_C0`.
            return 0;
        }
        self.curves.first().map(|c| c.continuity()).unwrap_or(0)
    }

    /// `NbIntervals(S)` (`cxx:270-279`).
    pub fn nb_intervals(&self, continuity: u8) -> i32 {
        self.curves
            .iter()
            .map(|c| c.nb_intervals(continuity))
            .sum()
    }

    /// `Intervals(T, S)` (`cxx:281-335`), `NbIntervals(S) + 1` knots.
    pub fn intervals(&self, continuity: u8) -> Vec<f64> {
        let nb = self.nb_intervals(continuity);
        if nb <= 0 {
            return Vec::new();
        }
        let mut t = vec![0.0; (nb + 1) as usize];
        let mut kk = 0usize;

        // First curve: direction of parsing of the edge (`cxx:287-306`).
        let ti = self.curves[0].parameter_intervals(continuity);
        let (f, delta) = self.inv_prepare(0);
        let base = self.knots[0];
        if delta < 0.0 {
            for jj in (0..ti.len()).rev() {
                t[kk] = base + (ti[jj] - f) * delta;
                kk += 1;
            }
        } else {
            for value in ti.iter() {
                t[kk] = base + (value - f) * delta;
                kk += 1;
            }
        }

        // And the next ones (`cxx:308-334`).
        for ii in 1..self.curves.len() {
            let ti = self.curves[ii].parameter_intervals(continuity);
            let (f, delta) = self.inv_prepare(ii);
            let base = self.knots[ii];
            if delta < 0.0 {
                // `for (jj = Ti->Length() - 1; jj > 0; kk++, jj--)`.
                for jj in (0..ti.len().saturating_sub(1)).rev() {
                    t[kk] = base + (ti[jj] - f) * delta;
                    kk += 1;
                }
            } else {
                // `for (jj = 2; jj <= Ti->Length(); kk++, jj++)`.
                for jj in 1..ti.len() {
                    t[kk] = base + (ti[jj] - f) * delta;
                    kk += 1;
                }
            }
        }
        t
    }

    /// `Trim(First, Last, Tol)` (`cxx:337-345`): a fresh adaptor built from the
    /// same wire and the new limits.
    pub fn trim(&self, first: f64, last: f64, tol: f64) -> Option<Self> {
        let w = self.wire.as_ref()?;
        Some(Self::with_range(w, self.is_by_ac, first, last, tol))
    }

    /// `IsClosed()` (`cxx:347-350`).
    pub fn is_closed(&self) -> bool {
        self.wire.as_ref().map(|w| w.0.is_closed()).unwrap_or(false)
    }

    /// `IsPeriodic()` (`cxx:352-355`).
    pub fn is_periodic(&self) -> bool {
        false
    }

    /// `Period()` (`cxx:357-360`).
    pub fn period(&self) -> f64 {
        self.t_last - self.t_first
    }

    /// `Resolution(R3d)` (`cxx:410-422`).
    pub fn resolution(&self, r3d: f64) -> f64 {
        let mut res = 1.0e200_f64;
        for c in &self.curves {
            let r = c.resolution(r3d);
            if r < res {
                res = r;
            }
        }
        res
    }

    /// `GetType()` (`cxx:424-430`): always `GeomAbs_OtherCurve`.
    pub fn get_type(&self) -> CurveType {
        CurveType::OtherCurve
    }

    /// `Prepare(W, Delta, theCurIndex)` (`cxx:497-568`).
    fn prepare(&self, w: &mut f64, delta: &mut f64, index: &mut i32) {
        // `cxx:500-508`: the offset sign selects the side of a knot.
        let eps = if *w - self.t_first < self.t_last - *w {
            self.p_tol
        } else {
            -self.p_tol
        };
        let w_test = *w + eps;

        // `cxx:510-543`: find the index.
        let mut cur = *index as usize;
        let mut found = false;
        if self.knots[cur - 1] > w_test {
            let mut ii = cur - 1;
            while ii > 0 && !found {
                if self.knots[ii - 1] <= w_test {
                    cur = ii;
                    found = true;
                }
                ii -= 1;
            }
            if !found {
                // Out of limits.
                cur = 1;
            }
        } else if self.knots[cur] <= w_test {
            let mut ii = cur + 1;
            while ii <= self.curves.len() && !found {
                if self.knots[ii] > w_test {
                    cur = ii;
                    found = true;
                }
                ii += 1;
            }
            if !found {
                cur = self.curves.len();
            }
        }
        *index = cur as i32;

        // `cxx:545-549`: invert?
        let edge = &self.edges[cur - 1];
        let or = edge.0.orientation();
        let reverse = (self.forward && or == Orientation::Reversed)
            || (!self.forward && or != Orientation::Reversed);

        // `cxx:551-566`: the local parameter.
        let (f, l) = BRepTool::edge_parameters(edge);
        *delta = self.knots[cur] - self.knots[cur - 1];
        if *delta > self.p_tol * 1.0e-9 {
            *delta = (l - f) / *delta;
        }
        if reverse {
            *delta = -*delta;
            *w = l + (*w - self.knots[cur - 1]) * *delta;
        } else {
            *w = f + (*w - self.knots[cur - 1]) * *delta;
        }
    }

    /// `InvPrepare(index, First, Delta)` (`cxx:570-597`).
    fn inv_prepare(&self, index: usize) -> (f64, f64) {
        let edge = &self.edges[index];
        let or = edge.0.orientation();
        let reverse = (self.forward && or == Orientation::Reversed)
            || (!self.forward && or != Orientation::Reversed);

        let (f, l) = BRepTool::edge_parameters(edge);
        let mut delta = self.knots[index + 1] - self.knots[index];
        if l - f > self.p_tol * 1.0e-9 {
            delta /= l - f;
        }
        if reverse {
            (l, -delta)
        } else {
            (f, delta)
        }
    }

    /// The `Prepare` lookup shared by `Edge(U)` and every `EvalD*`:
    /// `(0-based curve index, local parameter, reparametrisation factor)`.
    fn locate(&self, u: f64) -> Option<(usize, f64, f64)> {
        if self.curves.is_empty() || self.cur_index < 1 {
            return None;
        }
        let mut w = u;
        let mut d = 0.0;
        let mut index = self.cur_index;
        self.prepare(&mut w, &mut d, &mut index);
        if index < 1 || index as usize > self.curves.len() {
            return None;
        }
        Some(((index - 1) as usize, w, d))
    }

    /// `EvalD0(theU)` (`cxx:362-367`).
    fn eval_d0(&self, u: f64) -> GpPnt {
        let (index, local, _) = self
            .locate(u)
            .expect("BRepAdaptor_CompCurve: no curve loaded (Standard_NullObject)");
        <CurveAdaptor as Curve>::d0(&self.curves[index], local)
    }

    /// `EvalD1(theU)` (`cxx:369-376`).
    fn eval_d1(&self, u: f64) -> (GpPnt, GpVec) {
        let (index, local, d) = self
            .locate(u)
            .expect("BRepAdaptor_CompCurve: no curve loaded (Standard_NullObject)");
        let (p, d1) = <CurveAdaptor as Curve>::d1(&self.curves[index], local);
        (p, d1.multiplied_scalar(d))
    }

    /// `EvalD2(theU)` (`cxx:378-387`).
    fn eval_d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        let (index, local, d) = self
            .locate(u)
            .expect("BRepAdaptor_CompCurve: no curve loaded (Standard_NullObject)");
        let (p, d1, d2) = <CurveAdaptor as Curve>::d2(&self.curves[index], local);
        (
            p,
            d1.multiplied_scalar(d),
            d2.multiplied_scalar(d * d),
        )
    }

    /// `EvalD3(theU)` (`cxx:389-400`).
    fn eval_d3(&self, u: f64) -> (GpPnt, GpVec, GpVec, GpVec) {
        let (index, local, d) = self
            .locate(u)
            .expect("BRepAdaptor_CompCurve: no curve loaded (Standard_NullObject)");
        let (p, d1, d2, d3) = <CurveAdaptor as Curve>::d3(&self.curves[index], local);
        (
            p,
            d1.multiplied_scalar(d),
            d2.multiplied_scalar(d * d),
            d3.multiplied_scalar(d * d * d),
        )
    }
}

/// `Adaptor3d_Curve` view: `BRepFill_SectionPlacement` hands the adaptor to
/// `GeomFill_SectionPlacement::Perform(Adaptor3d_Curve, Tol)`
/// (`BRepFill_SectionPlacement.cxx:170-171`), which the port expresses as
/// `&dyn Curve`.
impl Curve for CompCurveAdaptor {
    fn d0(&self, u: f64) -> GpPnt {
        self.eval_d0(u)
    }

    fn d1(&self, u: f64) -> (GpPnt, GpVec) {
        self.eval_d1(u)
    }

    fn d2(&self, u: f64) -> (GpPnt, GpVec, GpVec) {
        self.eval_d2(u)
    }

    fn d3(&self, u: f64) -> (GpPnt, GpVec, GpVec, GpVec) {
        self.eval_d3(u)
    }

    /// `EvalDN(theU, theN)` (`cxx:402-408`): `D{N}` scaled by `Delta^N`.
    fn eval_dn(&self, u: f64, n: i32) -> GpVec {
        let (index, local, d) = self
            .locate(u)
            .expect("BRepAdaptor_CompCurve: no curve loaded (Standard_NullObject)");
        self.curves[index]
            .eval_dn(local, n)
            .multiplied_scalar(d.powi(n))
    }

    fn first_parameter(&self) -> f64 {
        self.t_first
    }

    fn last_parameter(&self) -> f64 {
        self.t_last
    }

    fn continuity(&self) -> u8 {
        CompCurveAdaptor::continuity(self)
    }

    fn nb_intervals(&self, continuity: u8) -> i32 {
        CompCurveAdaptor::nb_intervals(self, continuity)
    }

    fn parameter_intervals(&self, continuity: u8) -> Vec<f64> {
        self.intervals(continuity)
    }

    fn resolution(&self, r3d: f64) -> f64 {
        CompCurveAdaptor::resolution(self, r3d)
    }

    fn is_periodic(&self) -> bool {
        false
    }

    fn period(&self) -> f64 {
        CompCurveAdaptor::period(self)
    }

    /// `Line()` (`cxx:432-435`) / `Circle()` (`cxx:437-440`) and the other
    /// analytic accessors forward to `myCurves->Value(1)`; `GetType()` stays
    /// `GeomAbs_OtherCurve`, so a caller that checks the type never reaches
    /// them, exactly as in OCCT.
    fn gp_line(&self) -> Option<occt_core::gp::GpLin> {
        self.curves.first().and_then(|c| c.gp_line())
    }

    fn gp_circ(&self) -> Option<occt_core::gp::GpCirc> {
        self.curves.first().and_then(|c| c.gp_circ())
    }

    fn gp_ellipse(&self) -> Option<occt_core::gp::GpElips> {
        self.curves.first().and_then(|c| c.gp_ellipse())
    }

    fn gp_hyperbola(&self) -> Option<occt_core::gp::GpHypr> {
        self.curves.first().and_then(|c| c.gp_hyperbola())
    }

    fn gp_parabola(&self) -> Option<occt_core::gp::GpParab> {
        self.curves.first().and_then(|c| c.gp_parabola())
    }

    /// `Degree()` (`cxx:457-460`).
    fn nurbs_degree(&self) -> Option<usize> {
        self.curves.first().and_then(|c| c.nurbs_degree())
    }

    /// `NbPoles()` (`cxx:466-469`).
    fn bspline_poles(&self) -> Option<&[GpPnt]> {
        self.curves.first().and_then(|c| c.bspline_poles())
    }

    /// `IsRational()` (`cxx:462-465`) through the weight slice.
    fn bspline_weights(&self) -> Option<&[f64]> {
        self.curves.first().and_then(|c| c.bspline_weights())
    }

    /// `Bezier()` (`cxx:477-480`).
    fn bezier_poles(&self) -> Option<&[GpPnt]> {
        self.curves.first().and_then(|c| c.bezier_poles())
    }

    /// `Bezier()` weights.
    fn bezier_weights(&self) -> Option<&[f64]> {
        self.curves.first().and_then(|c| c.bezier_weights())
    }

    /// `Adaptor3d_Curve` has no `Transform`.
    fn transform(&mut self, _t: &occt_core::gp::GpTrsf) {}

    /// `Adaptor3d_Curve` has no `Reverse`.
    fn reverse(&mut self) {}

    fn clone_dyn(&self) -> Box<dyn Curve> {
        Box::new(self.clone())
    }
}
