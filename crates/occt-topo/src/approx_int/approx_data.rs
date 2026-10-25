//! The data block of `ApproxInt_Approx` plus the three file-static helpers
//! that fill it — split out of [`super::approx`] so that neither file crosses
//! the size budget. All four items stay in the `approx_int` namespace.
//!
//! Sources: `ModelingAlgorithms/TKGeomAlgo/ApproxInt/ApproxInt_Approx.gxx`
//! - `ComputeTrsf3d` (`:36-53`)
//! - `ComputeTrsf2d` (`:57-82`)
//! - `TheMultiLineTool::Value(Line, Index, …)` selection by `NbP3d` / `NbP2d`
//!   (`:115-128`, `:577-581`)
//! - `Approx_Data` (`BRepApprox_Approx.hxx:46-78`)

use occt_core::approx::ApproxParametrizationType;
use occt_core::gp::{GpPnt, GpPnt2d};
use occt_geom::approx::line_tool::LineTool;

use crate::approx_int::multi_line::LineSource;

/// `Approx_Data` (`BRepApprox_Approx.hxx:46-78`): the flags and offsets of the
/// current `Perform`.
#[derive(Debug, Clone)]
pub struct ApproxData {
    /// `myBezierApprox`.
    pub my_bezier_approx: bool,
    /// `Xo`, `Yo`, `Zo`.
    pub xo: f64,
    /// `Yo`.
    pub yo: f64,
    /// `Zo`.
    pub zo: f64,
    /// `U1o`.
    pub u1o: f64,
    /// `V1o`.
    pub v1o: f64,
    /// `U2o`.
    pub u2o: f64,
    /// `V2o`.
    pub v2o: f64,
    /// `ApproxXYZ`.
    pub approx_xyz: bool,
    /// `ApproxU1V1`.
    pub approx_u1v1: bool,
    /// `ApproxU2V2`.
    pub approx_u2v2: bool,
    /// `indicemin`.
    pub indicemin: i32,
    /// `indicemax`.
    pub indicemax: i32,
    /// `myNbPntMax`.
    pub my_nb_pnt_max: i32,
    /// `parametrization`.
    pub parametrization: ApproxParametrizationType,
}

impl Default for ApproxData {
    /// `Approx_Data()` (`hxx:48-70`).
    fn default() -> Self {
        Self {
            my_bezier_approx: true,
            xo: 0.0,
            yo: 0.0,
            zo: 0.0,
            u1o: 0.0,
            v1o: 0.0,
            u2o: 0.0,
            v2o: 0.0,
            approx_xyz: true,
            approx_u1v1: true,
            approx_u2v2: true,
            indicemin: 0,
            indicemax: 0,
            my_nb_pnt_max: 30,
            parametrization: ApproxParametrizationType::ChordLength,
        }
    }
}

/// `ComputeTrsf3d` (`gxx:36-53`).
pub(crate) fn compute_trsf3d<L: LineSource>(theline: &L) -> (f64, f64, f64) {
    let a_nb_pnts = theline.nb_pnts();
    let mut a_xmin = f64::MAX;
    let mut a_ymin = f64::MAX;
    let mut a_zmin = f64::MAX;
    for i in 1..=a_nb_pnts {
        let p = theline.point(i).value();
        a_xmin = p.x().min(a_xmin);
        a_ymin = p.y().min(a_ymin);
        a_zmin = p.z().min(a_zmin);
    }
    (-a_xmin, -a_ymin, -a_zmin)
}

/// `ComputeTrsf2d` (`gxx:57-82`).
pub(crate) fn compute_trsf2d<L: LineSource>(theline: &L, on_first: bool) -> (f64, f64) {
    let a_nb_pnts = theline.nb_pnts();
    let mut a_umin = f64::MAX;
    let mut a_vmin = f64::MAX;
    for i in 1..=a_nb_pnts {
        let (u1, v1, u2, v2) = theline.point(i).parameters();
        let (u, v) = if on_first { (u1, v1) } else { (u2, v2) };
        a_umin = u.min(a_umin);
        a_vmin = v.min(a_vmin);
    }
    (-a_umin, -a_vmin)
}

/// `TheMultiLineTool::Value(Line, Index, …)` for the three arities OCCT picks
/// by the `NbP3d` / `NbP2d` flags (`gxx:115-128`, `:577-581`).
pub(crate) fn line_values<L: LineTool>(
    line: &L,
    index: i32,
    nb_p3d: i32,
    nb_p2d: i32,
) -> (Vec<GpPnt>, Vec<GpPnt2d>) {
    if nb_p3d != 0 && nb_p2d != 0 {
        line.value(index)
    } else if nb_p2d != 0 {
        (Vec::new(), line.value_2d(index))
    } else if nb_p3d != 0 {
        (line.value_3d(index), Vec::new())
    } else {
        (Vec::new(), Vec::new())
    }
}
