//! `ShapeFix_Face` subset: `CheckWire` (`ShapeFix_Face.cxx:1652-1718`), the
//! setup / wire-pair selection plus seam construction of `FixMissingSeam`
//! (`ShapeFix_Face.cxx:1722-2330`), `IsPeriodicConicalLoop`
//! (`cxx:3018-3098`) and `FixPeriodicDegenerated` (`cxx:3101-3259`).
//!
//! UNPORTED: `ShapeFix_Face::Perform` (`cxx:345-...`) beyond its
//! `FixPeriodicDegenerated` / `FixMissingSeam` tail (`cxx:482-498`) — the
//! wire-fixing first part (`cxx:365-480`, needs `ShapeFix_Wire::Perform`) and
//! the post-`FixMissingSeam` face loop (`cxx:500-...`). See
//! `specs/_a3n00_gap_analysis.md` §9.35.

mod prelude {
    pub(crate) use std::sync::Arc;
    pub(crate) use std::f64::consts::PI;
    pub(crate) use occt_core::bnd::BndBox2d;
    pub(crate) use occt_core::gp::{GpDir2d, GpPnt2d, GpVec2d};
    pub(crate) use occt_core::precision::{Precision, CONFUSION, INFINITE, PCONFUSION};
    pub(crate) use occt_geom::{GeomRectangularTrimmedSurface, Surface};
    pub(crate) use occt_geom2d::curve::Curve2d;
    pub(crate) use occt_geom2d::line::Geom2dLine;
    pub(crate) use crate::abs::{Orientation, ShapeType};
    pub(crate) use crate::boptools_2d;
    pub(crate) use crate::brep_tool::BRepTool;
    pub(crate) use crate::brep_tools::{self, add_uv_bounds_on_wire};
    pub(crate) use crate::builder::TopoBuilder;
    pub(crate) use crate::shape::{Edge, Face, TopoShape, Wire};
    pub(crate) use crate::shape_analysis::is_outer_bound;
    pub(crate) use crate::shape_fix_compose_shell::{
        compose_adjust_to_period, reverse_wire_data_on_face, CompositeSurface, ComposeShell, Parametrisation,
        ReShape,
    };
    pub(crate) use crate::shhealing::{adjust_by_period, first_vertex, fix_reorder_wire, last_vertex};
    pub(crate) use crate::tgeometry::GeometryRegistry;
    pub(crate) use crate::topo_tools::{edges_of_wire, wires_of_face};
}

use prelude::*;

/// `CheckWire(wire, face, dU, dV, isuopen, isvopen, isDeg)`
/// (`ShapeFix_Face.cxx:1652-1718`).
pub fn healing_check_wire(
    wire: &Wire,
    face: &Face,
    d_u: f64,
    d_v: f64,
) -> Option<(i32, i32, bool)> {
    let mut vec_x = 0.0f64;
    let mut vec_y = 0.0f64;
    let mut is_deg = true;
    for edge in edges_of_wire(wire) {
        if !BRepTool::is_degenerated(&edge) {
            is_deg = false; // cxx:1670-1673
        }
        let (c2d, f, l) = match boptools_2d::curve_on_surface_oriented(&edge, face, true) {
            Some(v) => v,
            None => {
                return None; // cxx:1676-1679
            }
        };
        let pl = c2d.d0(l);
        let pf = c2d.d0(f);
        vec_x += pl.x() - pf.x(); // cxx:1680
        vec_y += pl.y() - pf.y();
    }

    // cxx:1683-1715.
    let isuopen = {
        let delta = vec_x.abs() - d_u;
        if delta.abs() < 0.1 * d_u {
            if vec_x > 0.0 { 1 } else { -1 }
        } else {
            0
        }
    };
    let isvopen = {
        let delta = vec_y.abs() - d_v;
        if delta.abs() < 0.1 * d_v {
            if vec_y > 0.0 { 1 } else { -1 }
        } else {
            0
        }
    };
    if isuopen != 0 || isvopen != 0 {
        Some((isuopen, isvopen, is_deg))
    } else {
        None
    }
}

/// `IsSurfaceUVPeriodic(theSurf)` (`ShapeFix_Face.cxx:103-108`).
pub(in crate::shhealing) fn is_surface_uv_periodic(surf: &dyn Surface) -> bool {
    (surf.is_u_periodic() && surf.is_v_periodic()) || surf.gp_sphere().is_some()
}

/// `ShapeExtend_WireData::Wire()` (`ShapeExtend_WireData.cxx:651-684`): a fresh
/// FORWARD wire over `data`, then the non-manifold edges appended unchanged
/// (`cxx:678-683`), with the closed flag set when the first and last vertex
/// coincide (`cxx:668-677`).
pub(in crate::shhealing) fn wire_from_wire_data(data: &[Edge], non_manifold: &[Edge]) -> Wire {
    let builder = TopoBuilder::new();
    let mut wire = builder.make_wire(data);
    let is_manifold = data.iter().all(|e| {
        let o = e.0.orientation();
        o == Orientation::Forward || o == Orientation::Reversed
    });
    if is_manifold {
        let vf = data.first().and_then(first_vertex);
        let vl = data.last().and_then(last_vertex);
        if let (Some(vf), Some(vl)) = (vf, vl) {
            if crate::topo_tools::is_same(&vf.0, &vl.0) {
                wire.0.set_closed(true); // cxx:670-676 `W.Closed(true)`
            }
        }
    }
    for e in non_manifold {
        builder.add_edge(&mut wire, e); // cxx:678-683
    }
    wire
}

/// `IsPeriodicConicalLoop` out-parameters (`ShapeFix_Face.cxx:3018-3098`).
#[derive(Debug, Clone, Copy)]
struct PeriodicConicalLoop {
    min_u: f64,
    max_u: f64,
    min_v: f64,
    max_v: f64,
    is_u_decrease: bool,
}

/// `IsPeriodicConicalLoop(theSurf, theWire, theTolerance, theMinU, theMaxU,
/// theMinV, theMaxV, isUDecrease)` (`ShapeFix_Face.cxx:3018-3098`): true when
/// the wire belts the conical surface by exactly one period (the absolute sum
/// of the pcurve `dU` is `2*pi` and the U span covers the apex). `None` is
/// `false`.
///
/// `cxx:3027-3032` (null surface) is the caller's `BRep_Tool::Surface(myFace)`
/// test; `ShapeAnalysis_Edge::PCurve` here is
/// `boptools_2d::curve_on_surface_oriented` keyed by the face's surface
/// (`tgeometry/registry_faces.rs::repr_key`), which is the same handle OCCT looks the pcurve
/// up by.
fn is_periodic_conical_loop(
    wire: &Wire,
    face: &Face,
    tolerance: f64,
) -> Option<PeriodicConicalLoop> {
    let mut cumul_delta_u = 0.0f64;
    let mut cumul_delta_u_abs = 0.0f64;
    // `cxx:3037-3040`: `RealLast()` / `-RealLast()` (`Standard_Real.hxx:179`).
    let mut min_u = f64::MAX;
    let mut min_v = min_u;
    let mut max_u = -min_u;
    let mut max_v = max_u;

    // `cxx:3043`: `TopoDS_Iterator(theWire, false)`.
    for edge in edges_of_wire(wire) {
        let (c2d, p_first, p_last) = match boptools_2d::curve_on_surface_oriented(&edge, face, true)
        {
            Some(v) => v,
            None => return None, // cxx:3051-3054
        };

        let uv_first = c2d.d0(p_first);
        let uv_last = c2d.d0(p_last);

        let (u_first, u_last) = (uv_first.x(), uv_last.x());
        let (v_first, v_last) = (uv_first.y(), uv_last.y());

        // cxx:3062-3079.
        let cur_max_u = u_first.max(u_last);
        let cur_min_u = u_first.min(u_last);
        let cur_max_v = v_first.max(v_last);
        let cur_min_v = v_first.min(v_last);
        if cur_min_u < min_u {
            min_u = cur_min_u;
        }
        if cur_max_u > max_u {
            max_u = cur_max_u;
        }
        if cur_min_v < min_v {
            min_v = cur_min_v;
        }
        if cur_max_v > max_v {
            max_v = cur_max_v;
        }

        // cxx:3081-3084.
        let delta_u = u_last - u_first;
        cumul_delta_u += delta_u;
        cumul_delta_u_abs += delta_u.abs();
    }

    // cxx:3094-3098.
    let is_2pi_delta = (cumul_delta_u_abs - 2.0 * PI).abs() <= tolerance;
    let is_around_apex = (max_u - min_u).abs() > 2.0 * PI - tolerance;
    if is_2pi_delta && is_around_apex {
        Some(PeriodicConicalLoop {
            min_u,
            max_u,
            min_v,
            max_v,
            is_u_decrease: cumul_delta_u < 0.0,
        })
    } else {
        None
    }
}
/// `ShapeFix_Face` state needed by `FixMissingSeam` (`ShapeFix_Face.hxx`).
pub struct ShapeFixFace {
    pub face: Option<Face>,
    pub surf: Option<Arc<dyn Surface>>,
    pub status: i32,
    /// `ShapeFix_Root::myContext` (`ShapeBuild_ReShape`).
    ///
    /// `ShapeFix_ComposeShell` is handed the same handle
    /// (`SetContext(Context())`, `ShapeFix_ComposeShell.cxx:2259`).
    pub context: crate::shape_fix_compose_shell::SharedReShape,
    /// `ShapeFix_Root::MaxTolerance` (`FromSTEP.FixShape.MaxTolerance3d`).
    pub max_tol: f64,
    /// `ShapeFix_Face::myResult`.
    pub result: Option<TopoShape>,
    /// `myFixMissingSeamMode` (`ShapeFix_Face.cxx:137`, default `-1` = auto).
    pub fix_missing_seam_mode: bool,
    /// `myFixPeriodicDegenerated` (`ShapeFix_Face.cxx:144`, default `-1`).
    /// `NeedFix(-1)` is true (`ShapeFix_Root.lxx:101`) and the STEP
    /// `ShapeProcess` operator never sets this flag — it only sets
    /// `FixMissingSeamMode` (`ShapeProcess_OperLibrary.cxx:830`) — so the
    /// conic degenerate-apex fix always runs on the reader path.
    pub fix_periodic_degenerated_mode: bool,
    /// `myFixOrientationMode` (`ShapeFix_Face.cxx:141`, default `-1` = auto).
    /// `FromSTEP.FixShape` sets it from `FixOrientationMode`
    /// (`ShapeProcess_OperLibrary.cxx:828`), whose default is also `-1`, so
    /// `NeedFix` is true on the reader path.
    pub fix_orientation_mode: bool,
    /// `myFixIntersectingWiresMode` (`ShapeFix_Face.cxx:140`, default `-1` =
    /// auto). `ShapeProcess_OperLibrary.cxx:833` sets it from
    /// `FixIntersectingWiresMode`, whose default is also `-1`, so `NeedFix` is
    /// true on the reader path.
    pub fix_intersecting_wires_mode: bool,
    /// `myFixAddNaturalBoundMode` (`ShapeFix_Face.cxx:143`, default `-1`);
    /// `NeedFix` is true, so `isNeedAddNaturalBound` runs.
    pub fix_add_natural_bound_mode: bool,
    /// `myFixSplitFaceMode` (`ShapeFix_Face.cxx:142`, default `-1` = auto).
    /// `Perform` gates `FixSplitFace` on `NeedFix(myFixSplitFaceMode)`
    /// (`cxx:711`), so the face-splitting tail of the second-part loop runs.
    pub fix_split_face_mode: bool,
    /// `Perform`'s local `MapWires` (`ShapeFix_Face.cxx:689-691`): the outer
    /// wire -> inner wires nesting that `FixOrientation` computes and
    /// `FixSplitFace` consumes (`cxx:711-716`).
    ///
    /// Outer wire to the wires it contains. Written by `orient_several_wires`
    /// (`ShapeFix_Face.cxx:1537`, `:1542`, `:1579`, `:1583`).
    pub map_wires: Vec<(Wire, Vec<Wire>)>,
    /// `ShapeFix_Root::myFwd` (`ShapeFix_Face.cxx:218`):
    /// `theFace.Orientation() != TopAbs_REVERSED`. `FixOrientation` re-applies
    /// it to the rebuilt face (`cxx:1634-1636`).
    pub fwd: bool,
}

impl Default for ShapeFixFace {
    fn default() -> Self {
        Self {
            face: None,
            surf: None,
            status: 0,
            context: crate::shape_fix_compose_shell::SharedReShape::new(),
            max_tol: 1.0,
            result: None,
            fix_missing_seam_mode: true,
            fix_periodic_degenerated_mode: true,
            fix_orientation_mode: true,
            fix_intersecting_wires_mode: true,
            fix_add_natural_bound_mode: true,
            fix_split_face_mode: true,
            map_wires: Vec::new(),
            fwd: true,
        }
    }
}
impl ShapeFixFace {
    /// `ShapeFix_Face(face)` with an externally owned context: OCCT's
    /// `ShapeFix_Shape` / `ShapeFix_Shell` call `myFixFace->SetContext(Context())`
    /// (`ShapeFix_Shell.cxx:108`, `ShapeFix_Shape.cxx:202`) before `Init`, so
    /// every face of one `FixShape` pass shares one `ShapeBuild_ReShape`.
    pub fn with_face_and_context(
        face: &Face,
        context: crate::shape_fix_compose_shell::SharedReShape,
    ) -> Self {
        let mut s = Self::with_face(face);
        s.context = context;
        s
    }

    /// `ShapeFix_Face(face)` + `mySurf = ShapeAnalysis_Surface(surface)`
    /// (`ShapeFix_Face.cxx:181-224`).
    pub fn with_face(face: &Face) -> Self {
        Self {
            face: Some(face.clone()),
            surf: BRepTool::face_surface(face),
            status: 0,
            context: crate::shape_fix_compose_shell::SharedReShape::new(),
            max_tol: 1.0,
            result: None,
            fix_missing_seam_mode: true,
            fix_periodic_degenerated_mode: true,
            fix_orientation_mode: true,
            fix_intersecting_wires_mode: true,
            fix_add_natural_bound_mode: true,
            fix_split_face_mode: true,
            map_wires: Vec::new(),
            // `myFwd = (theFace.Orientation() != TopAbs_REVERSED)`
            // (`ShapeFix_Face.cxx:218`): INTERNAL/EXTERNAL count as forward.
            fwd: face.0.orientation() != Orientation::Reversed,
        }
    }

    /// `ShapeFix_Face::Perform` (`ShapeFix_Face.cxx:345-...`) reduced to the
    /// `FixPeriodicDegenerated` / `FixMissingSeam` tail (`cxx:482-498`):
    /// `myResult = myFace;` then, when `myFixPeriodicDegeneratedMode` is on,
    /// `FixPeriodicDegenerated()` (`cxx:486-489`), then `FixMissingSeam()`
    /// (`cxx:492-498`).
    ///
    /// UNPORTED: the wire-fixing first part (`cxx:365-480`, needs
    /// `ShapeFix_Wire::Perform`) and the post-`FixMissingSeam` face loop
    /// (`cxx:500-...`).
    pub fn perform_fix_missing_seam(&mut self) -> Option<TopoShape> {
        self.result = self.face.as_ref().map(|f| f.0.clone()); // cxx:482
        if self.fix_periodic_degenerated_mode {
            self.fix_periodic_degenerated(); // cxx:486-489
        }
        if self.fix_missing_seam_mode {
            // cxx:492-498.
            if self.fix_missing_seam() {
                self.status |= 0x0004; // ShapeExtend_DONE3
            }
        }
        self.result.clone()
    }

    /// `ShapeFix_Face::FixIntersectingWires` (`ShapeFix_Face.cxx:2821-2825`):
    /// builds `ShapeFix_IntersectionTool(Context(), Precision(), MaxTolerance())`
    /// and runs `ITool.FixIntersectingWires(myFace)`.
    ///
    /// The tool rebuilds the face and records `Context()->Replace(myFace, newface)`
    /// (`intersection_tool::fix_intersecting_wires`, `cxx:2526`), so `myFace` here
    /// follows the out-parameter update. `ShapeFix_Face::Perform` (`cxx:680`)
    /// consumes this; wiring that driver still needs the unported wire-fixing
    /// first part (`cxx:365-480`).
    pub fn fix_intersecting_wires(&mut self) -> bool {
        let mut face = match self.face.clone() {
            Some(f) => f,
            None => return false,
        };
        let done = crate::shhealing::fix_intersecting_wires(&mut face, &mut self.context);
        if done {
            self.face = Some(face);
        }
        done
    }

    /// `ShapeFix_Face::isNeedAddNaturalBound` (`ShapeFix_Face.cxx:1120-1160`).
    pub(in crate::shhealing) fn is_need_add_natural_bound(&self, oriented_wires: &[Wire]) -> bool {
        // cxx:1124-1127.
        if !self.fix_add_natural_bound_mode {
            return false;
        }
        // cxx:1130-1133: only a surface closed in both U and V (or a sphere)
        // can be missing its natural bounds.
        if let Some(surf) = &self.surf {
            if !is_surface_uv_periodic(surf.as_ref()) {
                return false;
            }
        }
        let face = match &self.face {
            Some(f) => f.clone(),
            None => return false,
        };
        // cxx:1135-1138.
        if is_outer_bound(&face) {
            return false;
        }
        // cxx:1140-1157: a seam or degenerated edge means this wire is the
        // outer one already, and only its orientation gets corrected.
        for w in oriented_wires {
            for e in edges_of_wire(w) {
                if BRepTool::is_degenerated(&e) {
                    return false;
                }
                if BRepTool::is_closed_edge_face(&e, &face) {
                    return false;
                }
            }
        }
        true
    }

}

mod missing_seam;
mod orientation;
mod periodic_degenerated;
