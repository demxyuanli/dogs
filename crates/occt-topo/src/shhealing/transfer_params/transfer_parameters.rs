//! `ShapeAnalysis_TransferParameters`.

use super::prelude::*;
use super::curve_range::*;

/// `ShapeAnalysis_TransferParameters` (`ShapeAnalysis_TransferParameters.cxx`).
#[derive(Clone)]
pub struct TransferParameters {
    /// `myFirst` / `myLast`: the 3D range (`hxx:85-86`).
    pub(super) first: f64,
    pub(super) last: f64,
    /// `myEdge` (`hxx:87`).
    pub(super) edge: Edge,
    /// `myMaxTolerance` (`hxx:88`).
    pub(super) max_tolerance: f64,
    /// `myFirst2d` / `myLast2d` (`hxx:92-93`).
    pub(super) first2d: f64,
    pub(super) last2d: f64,
    /// `myFace` (`hxx:94`).
    pub(super) face: Option<Face>,
    /// `myShift` - private in OCCT (`hxx:91`).
    pub(super) shift: f64,
    /// `myScale` - private in OCCT (`hxx:90`).
    pub(super) scale: f64,
}

impl Default for TransferParameters {
    fn default() -> Self {
        Self::new()
    }
}

impl TransferParameters {
    /// `ShapeAnalysis_TransferParameters()` (`cxx:29-33`): `myShift = 0`,
    /// `myScale = 1`. OCCT leaves `myFirst` / `myLast` / `myFirst2d` /
    /// `myLast2d` / `myMaxTolerance` uninitialised here (they are only set by
    /// `Init`, `cxx:45-74`, or `SetMaxTolerance`, `cxx:79-82`); the port uses
    /// the `Proj` class's own fallbacks (`Proj.cxx:52` sets `myMaxTolerance = 1`)
    /// plus `0` / `1` for the ranges, which every call site overwrites.
    pub fn new() -> Self {
        Self {
            first: 0.0,
            last: 1.0,
            edge: Edge::new(),
            max_tolerance: 1.0,
            first2d: 0.0,
            last2d: 0.0,
            face: None,
            shift: 0.0,
            scale: 1.0,
        }
    }

    /// `ShapeAnalysis_TransferParameters(E, F)` (`cxx:37-41`).
    pub fn with_edge_face(e: &Edge, f: Option<&Face>) -> Self {
        let mut t = Self::new();
        t.init(e, f);
        t
    }

    /// `Init` (`cxx:45-74`).
    ///
    /// OCCT leaves `myFirst` / `myLast` uninitialised when the edge carries no
    /// 3D curve (`cxx:48` declares `double l, f;` and `cxx:55` only writes them
    /// on success) and `myFirst2d` / `myLast2d` at `0` when `F.IsNull()`
    /// (`cxx:49`); the port keeps the constructor's `0` / `1` / `0` / `0` in
    /// those cases.
    pub fn init(&mut self, e: &Edge, f: Option<&Face>) {
        self.scale = 1.0; // `cxx:46-47`
        self.shift = 0.0;
        let mut first2d = 0.0;
        let mut last2d = 0.0;
        self.edge = e.clone(); // `cxx:52`

        // `sae.Curve3d(E, curve3d, f, l, false)` (`cxx:54-56`): `orient = false`,
        // so the range is the stored one, not swapped for a REVERSED edge
        // (`ShapeAnalysis_Edge.cxx:100-125`).
        let mut curve3d: Option<Arc<dyn Curve>> = None;
        if let Some(curve) = BRepTool::edge_curve(e) {
            let (a, b) = BRepTool::edge_parameters(e);
            self.first = a;
            self.last = b;
            curve3d = Some(curve);
        }

        // `sae.PCurve(E, F, curve2d, f2d, l2d, false)` (`cxx:58-62`), skipped for
        // a null face (`cxx:59`: "process free edges").
        let mut curve2d: Option<Arc<dyn Curve2d>> = None;
        if let Some(f) = f {
            if let Some((pc, a, b)) = curve_on_surface_oriented(e, f, false) {
                curve2d = Some(pc);
                first2d = a;
                last2d = b;
            }
        }
        self.first2d = first2d;
        self.last2d = last2d;
        self.face = f.cloned();

        if curve3d.is_none() || curve2d.is_none() {
            return; // `cxx:66-69`
        }
        let ln2d = last2d - first2d; // `cxx:71`
        let ln3d = self.last - self.first; // `cxx:72`
        self.scale = if ln3d <= REAL_SMALL {
            1.0 // `cxx:73`: `ln3d <= gp::Resolution()`
        } else {
            ln2d / ln3d
        };
        self.shift = first2d - self.first * self.scale; // `cxx:74`
    }

    /// `SetMaxTolerance` (`cxx:79-82`).
    pub fn set_max_tolerance(&mut self, maxtol: f64) {
        self.max_tolerance = maxtol;
    }

    /// `Perform(Params, To2d)` (`cxx:86-96`).
    pub fn perform(&self, params: &[f64], to2d: bool) -> Vec<f64> {
        params.iter().map(|&p| self.perform_param(p, to2d)).collect()
    }

    /// `Perform(Param, To2d)` (`cxx:100-112`).
    pub fn perform_param(&self, param: f64, to2d: bool) -> f64 {
        if to2d {
            self.shift + param * self.scale // `cxx:105`
        } else {
            -self.shift / self.scale + param / self.scale // `cxx:109`
        }
    }

    /// `TransferRange(newEdge, prevPar, currPar, Is2d)` (`cxx:116-146`).
    pub fn transfer_range(&self, new_edge: &mut Edge, prev_par: f64, curr_par: f64, is2d: bool) {
        if is2d {
            let span2d = self.last2d - self.first2d; // `cxx:124`
            let (tmp1, tmp2) = if prev_par > curr_par {
                (curr_par, prev_par)
            } else {
                (prev_par, curr_par)
            }; // `cxx:126-132`
            let alpha = (tmp1 - self.first2d) / span2d; // `cxx:133`
            let beta = (tmp2 - self.first2d) / span2d; // `cxx:134`
            copy_ranges(new_edge, &self.edge, alpha, beta); // `cxx:135`
        } else {
            let alpha = (prev_par - self.first) / (self.last - self.first); // `cxx:142`
            let beta = (curr_par - self.first) / (self.last - self.first); // `cxx:143`
            copy_ranges(new_edge, &self.edge, alpha, beta); // `cxx:144`
        }
    }

    /// `IsSameRange` (`cxx:150-153`).
    pub fn is_same_range(&self) -> bool {
        self.shift == 0.0 && self.scale == 1.0
    }
}
