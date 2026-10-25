//! `TopOpeBRepTool_TOOL` isoparametric-identification family.
//!
//! Source (OCCT 8.0.0, `ModelingAlgorithms/TKBool/TopOpeBRepTool`):
//! - `UVISO(PC, isoU, isoV, d2d, o2d)` (`TopOpeBRepTool_TOOL.cxx:1191-1222`)
//! - `UVISO(E, F, isoU, isoV, d2d, o2d)` (`cxx:1224-1243`)
//! - `ParISO(uv, E, F, par)` (`cxx:544-576`)
//!
//! `UVISO(C2DF, ...)` (`cxx:1245-1259`) forwards the stored `PC` to the curve
//! arm ([`uv_iso_c2df`]).
//! - `IsonCLO(PC, onU, xfirst, xperiod, xtol)` (`cxx:1263-1296`)
//! - `IsonCLO(C2DF, onU, xfirst, xperiod, xtol)` (`cxx:1297-1307`)
//! - `UVF(par, C2DF)` (`cxx:533-541`)
//!
//! `BASISCURVE2D(PC)` (`TopOpeBRepTool_GEOMETRY.cxx:35-50`, used at
//! `TOOL.cxx:1199`) unwraps `Geom2d_OffsetCurve` and `Geom2d_TrimmedCurve`
//! recursively before the `Geom2d_Line` type test, so an offset or trimmed
//! p-curve of a line still answers an iso direction. [`basis_curve2d`] is that
//! recursion; `Curve2d::gp_lin2d` alone resolves through
//! `Geom2d_TrimmedCurve` only (`occt-geom2d/src/trimmed.rs:170-172`).

use occt_core::gp::{GpDir2d, GpPnt2d};
use occt_core::precision::PCONFUSION;
use occt_geom2d::curve::Curve2d;

use crate::abs::Orientation;
use crate::shape::{Edge, Face};
use crate::shape::TopoShape;
use crate::topopebrep::tool::c2d::C2dF;

/// `BASISCURVE2D(C)` (`TopOpeBRepTool_GEOMETRY.cxx:35-50`): descend through
/// `Geom2d_OffsetCurve` / `Geom2d_TrimmedCurve` until the curve is neither.
pub fn basis_curve2d<'a>(c: &'a dyn Curve2d) -> &'a dyn Curve2d {
    if let Some(basis) = c.offset_basis() {
        return basis_curve2d(basis);
    }
    if let Some(basis) = c.trimmed_basis() {
        return basis_curve2d(basis);
    }
    c
}

/// Out-parameters of [`uv_iso_curve`]: `(isoU, isoV, d2d, o2d)`.
pub struct UvIso {
    pub iso_u: bool,
    pub iso_v: bool,
    pub dir: GpDir2d,
    pub origin: GpPnt2d,
}

/// `UVISO(const occ::handle<Geom2d_Curve>& PC, ...)` (`cxx:1191-1222`): the
/// curve is a `Geom2d_Line` whose direction is parallel to a parametric axis.
///
/// `isoU = |d2d.X()| < Precision::Parametric(Precision::Confusion())` and
/// `isoV` likewise on `.Y()`; `Precision::Parametric(P)` is `P * 0.01`
/// (`Precision.hxx:328`), i.e. [`PCONFUSION`].
pub fn uv_iso_curve(pc: &dyn Curve2d) -> Option<UvIso> {
    let lin = basis_curve2d(pc).gp_lin2d()?;
    let dir = *lin.direction();
    let iso_u = dir.x.abs() < PCONFUSION;
    let iso_v = dir.y.abs() < PCONFUSION;
    if !(iso_u || iso_v) {
        return None;
    }
    Some(UvIso {
        iso_u,
        iso_v,
        dir,
        origin: lin.location(),
    })
}

/// `UVISO(const TopoDS_Edge& E, const TopoDS_Face& F, ...)` (`cxx:1224-1243`).
///
/// OCCT reads `FC2D_HasOldCurveOnSurface` then `FC2D_EditableCurveOnSurface`
/// and attaches the freshly built curve when none was stored (`cxx:1232-1239`);
/// [`crate::boptools_2d::make_2d`] is the port's stored-else-built equivalent
/// of that cascade.
pub fn uv_iso_edge_face(e: &TopoShape, f: &TopoShape) -> Option<UvIso> {
    let pc = crate::boptools_2d::make_2d(&Edge(e.clone()), &Face(f.clone())).ok()?;
    uv_iso_curve(pc.as_ref())
}

/// `ParISO(const gp_Pnt2d& uv, const TopoDS_Edge& E, const TopoDS_Face& F,
/// double& par)` (`cxx:544-576`): the parameter of `uv` on the isoparametric
/// line of `E` on `F`.
///
/// `par` starts at `1.e7` and is only overwritten by the branches that match
/// the curve's own iso axis; both `if`s (not `else if`) are written out as in
/// OCCT.
pub fn par_iso(uv: &GpPnt2d, e: &TopoShape, f: &TopoShape) -> Option<f64> {
    let iso = uv_iso_edge_face(e, f)?;
    let mut par = 1.0e7;
    if iso.iso_u {
        par = uv.y() - iso.origin.y();
        if iso.dir.y < 0.0 {
            par = -par;
        }
    }
    if iso.iso_v {
        par = uv.x() - iso.origin.x();
        if iso.dir.x < 0.0 {
            par = -par;
        }
    }
    Some(par)
}

/// `UVISO(const TopOpeBRepTool_C2DF& C2DF, ...)` (`cxx:1245-1259`): the same
/// test on the p-curve the `C2DF` stores. A `C2DF` without a curve answers
/// `false` at the `PC.IsNull()` gate of the curve arm (`cxx:1197-1201`).
pub fn uv_iso_c2df(c2df: &C2dF) -> Option<UvIso> {
    let (curve, _, _, _) = c2df.pc();
    uv_iso_curve(curve?.as_ref())
}

/// `IsonCLO(PC, onU, xfirst, xperiod, xtol)` (`cxx:1263-1296`): is the iso
/// line of `PC` the one at parameter `xfirst`, modulo `xperiod`?
///
/// `onU` selects the axis the caller asks about; the direction of the line is
/// not read, only its origin (`cxx:1283-1295`). Both the closed and the
/// periodic seam position count as "on the closed line".
pub fn iso_on_clo(pc: &dyn Curve2d, on_u: bool, xfirst: f64, xperiod: f64, xtol: f64) -> bool {
    let iso = match uv_iso_curve(pc) {
        Some(iso) => iso,
        None => return false,
    };
    let on_x = (on_u && iso.iso_u) || ((!on_u) && iso.iso_v);
    if !on_x {
        return false;
    }
    let dxx = if on_u {
        (iso.origin.x() - xfirst).abs()
    } else {
        (iso.origin.y() - xfirst).abs()
    };
    let mut onclo = dxx < xtol;
    onclo = onclo || ((xperiod - dxx).abs() < xtol);
    onclo
}

/// `IsonCLO(const TopOpeBRepTool_C2DF& C2DF, ...)` (`cxx:1297-1307`).
pub fn iso_on_clo_c2df(
    c2df: &C2dF,
    on_u: bool,
    xfirst: f64,
    xperiod: f64,
    xtol: f64,
) -> bool {
    let (curve, _, _, _) = c2df.pc();
    match curve {
        Some(pc) => iso_on_clo(pc.as_ref(), on_u, xfirst, xperiod, xtol),
        None => false,
    }
}

/// `UVF(par, C2DF)` (`cxx:533-541`): the `(u, v)` point of the stored p-curve
/// at `par`. A `C2DF` without a curve has no point to answer.
pub fn uvf(par: f64, c2df: &C2dF) -> Option<GpPnt2d> {
    let (curve, _, _, _) = c2df.pc();
    Some(curve?.d0(par))
}

/// `FUN_tool_orientEinFFORWARD`-free helper kept private to this module: the
/// `M_FORWARD` test of `myFi.Orientation()` used by the `xxi` reversals in
/// `TopOpeBRepTool_mkTondgE::MkTonE` (`TopOpeBRepTool_mkTondgE.cxx:317`,
/// `:451`).
pub(crate) fn orientation_is_forward(f: &TopoShape) -> bool {
    f.orientation() == Orientation::Forward
}
