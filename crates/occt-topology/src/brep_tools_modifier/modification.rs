//! `BRepTools_Modification` (`BRepTools_Modification.hxx:41-144`): the abstract
//! description of the geometric changes [`super::Modifier`] applies to a shape.
//!
//! OCCT answers "has this sub-shape been modified?" with a `bool` return value
//! and hands the new geometry back through out parameters. The port answers
//! with `Option`: `None` is OCCT's `false`, `Some(info)` carries exactly what
//! the out parameters carried. The three virtuals whose OCCT body is the
//! default "not modified" (`BRepTools_Modification.cxx:20-40`) are trait
//! defaults here.
//!
//! UNPORTED (no data model in this port):
//! - `NewPolygon` (`BRepTools_Modification.hxx:87`) — `EdgeGeom` keeps no
//!   `Poly_Polygon3D` (`tgeometry::EdgeGeom`), the same gap recorded in
//!   `bop_bop.rs:135-137`.
//! - `NewPolygonOnTriangulation` (`BRepTools_Modification.hxx:91-94`) —
//!   `EdgePcurves` keeps no `Poly_PolygonOnTriangulation`.
//! Their consumers, `BRepTools_Modifier.cxx:318-332` (copy polygon on the
//! edge) and `cxx:548-585` (copy polygon on triangulation), are marked at the
//! corresponding place in `rebuild.rs`.

use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_core::kernel::geomabs::Shape as GeomAbsShape;
use occt_core::poly::Triangulation;
use occt_core::toploc::TopLocLocation;
use occt_geom::{Curve, Surface};
use occt_geom2d::curve::Curve2d;

use crate::shape::TopoShape;

/// `BRepTools_Modifier::NewSurfaceInfo` (`BRepTools_Modifier.hxx:87-94`): the
/// five out parameters of `NewSurface` (`BRepTools_Modification.hxx:58-64`).
#[derive(Clone)]
pub struct NewSurfaceInfo {
    /// `mySurface`.
    pub surface: Arc<dyn Surface>,
    /// `myLoc` — the location the modification reports for the new surface.
    /// The rebuilt face keeps the geometry in `myLoc.Predivided(S.Location())`
    /// (`BRepTools_Modifier.cxx:261`); see `rebuild.rs`.
    pub location: TopLocLocation,
    /// `myToler`.
    pub tolerance: f64,
    /// `myRevWires` — set when the modification reverses the normal of the
    /// surface, so the wires of the rebuilt face have to be reversed
    /// (`BRepTools_Modifier.cxx:258`, `cxx:390-393`).
    pub rev_wires: bool,
    /// `myRevFace` — set when the orientation of the modified face changes in
    /// the shells that contain it (`BRepTools_Modifier.cxx:264-267`).
    pub rev_face: bool,
}

/// `BRepTools_Modifier::NewCurveInfo` (`BRepTools_Modifier.hxx:80-85`): the
/// out parameters of `NewCurve` (`BRepTools_Modification.hxx:78-81`).
#[derive(Clone)]
pub struct NewCurveInfo {
    /// `myCurve`. OCCT allows a **null** curve here and then takes the
    /// no-3D-curve branch of `Rebuild` (`BRepTools_Modifier.cxx:295-301`), so
    /// the port keeps the handle-null case as `None`.
    pub curve: Option<Arc<dyn Curve>>,
    /// `myLoc`.
    pub location: TopLocLocation,
    /// `myToler`.
    pub tolerance: f64,
}

/// Out parameters of `NewPoint` (`BRepTools_Modification.hxx:102-104`).
#[derive(Clone, Copy)]
pub struct NewPointInfo {
    /// `P`.
    pub point: GpPnt,
    /// `Tol`.
    pub tolerance: f64,
}

/// Out parameters of `NewParameter` (`BRepTools_Modification.hxx:127-130`).
#[derive(Clone, Copy)]
pub struct NewParameterInfo {
    /// `P`.
    pub parameter: f64,
    /// `Tol`.
    pub tolerance: f64,
}

/// `BRepTools_Modification` (`BRepTools_Modification.hxx:41-145`).
///
/// Every method takes the **source** sub-shape; the rebuilt shapes are only
/// handed to `new_curve2d`, exactly as OCCT passes `NewE` / `NewF` to
/// `NewCurve2d` (`BRepTools_Modification.hxx:114-119`).
pub trait Modification {
    /// `NewSurface` (`BRepTools_Modification.hxx:58-64`), pure virtual.
    fn new_surface(&self, f: &TopoShape) -> Option<NewSurfaceInfo>;

    /// `NewTriangulation` (`BRepTools_Modification.hxx:68-69`). The OCCT body is
    /// the default "not modified" (`BRepTools_Modification.cxx:25-28`).
    fn new_triangulation(&self, f: &TopoShape) -> Option<Arc<Triangulation>> {
        let _ = f;
        None
    }

    /// `NewCurve` (`BRepTools_Modification.hxx:78-81`), pure virtual.
    fn new_curve(&self, e: &TopoShape) -> Option<NewCurveInfo>;

    /// `NewPoint` (`BRepTools_Modification.hxx:102`), pure virtual.
    fn new_point(&self, v: &TopoShape) -> Option<NewPointInfo>;

    /// `NewCurve2d` (`BRepTools_Modification.hxx:114-119`), pure virtual.
    ///
    /// The `Geom2d_Curve& C` out parameter is the return value. OCCT also
    /// writes `Tol`, but neither caller reads it: `BRepTools_Modifier.cxx:491`
    /// / `:501` / `:508` pass a literal `0.`, and `cxx:522` / `cxx:605`
    /// overwrite it from `NewParameter` before use, so the port drops it.
    fn new_curve2d(
        &self,
        e: &TopoShape,
        f: &TopoShape,
        new_e: &TopoShape,
        new_f: &TopoShape,
    ) -> Option<Arc<dyn Curve2d>>;

    /// `NewParameter` (`BRepTools_Modification.hxx:127-130`), pure virtual.
    fn new_parameter(&self, v: &TopoShape, e: &TopoShape) -> Option<NewParameterInfo>;

    /// `Continuity` (`BRepTools_Modification.hxx:137-142`), pure virtual: the continuity of `new_e`
    /// between `new_f1` and `new_f2`.
    fn continuity(
        &self,
        e: &TopoShape,
        f1: &TopoShape,
        f2: &TopoShape,
        new_e: &TopoShape,
        new_f1: &TopoShape,
        new_f2: &TopoShape,
    ) -> GeomAbsShape;
}

/// `u8` (the numeric value of `GeomAbs_Shape`) to `GeomAbsShape`, the inverse
/// of the `as u8` conversion [`crate::tshape_continuity::set_stored_continuity`]
/// takes. `GeomAbs_C0` = 0 … `GeomAbs_CN` = 6 (`geomabs.rs:43`).
pub fn geom_abs_shape_from_u8(v: u8) -> GeomAbsShape {
    match v {
        1 => GeomAbsShape::G1,
        2 => GeomAbsShape::C1,
        3 => GeomAbsShape::G2,
        4 => GeomAbsShape::C2,
        5 => GeomAbsShape::C3,
        6 => GeomAbsShape::CN,
        _ => GeomAbsShape::C0,
    }
}
