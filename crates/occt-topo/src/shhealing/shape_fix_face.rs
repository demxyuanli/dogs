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

use std::sync::Arc;

use std::f64::consts::PI;

use occt_core::bnd::BndBox2d;
use occt_core::gp::{GpDir2d, GpPnt2d, GpVec2d};
use occt_core::precision::{Precision, CONFUSION, INFINITE, PCONFUSION};
use occt_geom::{GeomRectangularTrimmedSurface, Surface};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d;
use crate::brep_tool::BRepTool;
use crate::brep_tools::{self, add_uv_bounds_on_wire};
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::shape_analysis::is_outer_bound;
use crate::shape_fix_compose_shell::{
    adjust_to_period, CompositeSurface, ComposeShell, Parametrisation, ReShape,
};
use crate::shhealing::{adjust_by_period, fix_reorder_wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{edges_of_wire, wires_of_face};

/// `CheckWire(wire, face, dU, dV, isuopen, isvopen, isDeg)`
/// (`ShapeFix_Face.cxx:1652-1718`).
pub fn check_wire(
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
        let a_delta = vec_x.abs() - d_u;
        if a_delta.abs() < 0.1 * d_u {
            if vec_x > 0.0 { 1 } else { -1 }
        } else {
            0
        }
    };
    let isvopen = {
        let a_delta = vec_y.abs() - d_v;
        if a_delta.abs() < 0.1 * d_v {
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
/// (`tgeometry.rs::repr_key`), which is the same handle OCCT looks the pcurve
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

    /// `ShapeFix_Face::FixPeriodicDegenerated` (`ShapeFix_Face.cxx:3101-3259`):
    /// a conical face whose single wire belts the surface by one period gets a
    /// degenerated apex edge and a second (apex) wire, so that the following
    /// `FixMissingSeam` sees the 2-wire configuration and stitches the seam.
    /// Without it a cone read from STEP keeps 1 edge and 0 V-span, and the
    /// face's discrete range is empty (`IMeshData_Failure`).
    ///
    /// `Precision()` in this routine is `ShapeFix_Root::Precision()`
    /// (`ShapeFix_Root.lxx:34`) = `Precision::Confusion()` by default
    /// (`ShapeFix_Root.cxx:26`; the STEP operator does not `SetPrecision`).
    pub fn fix_periodic_degenerated(&mut self) -> bool {
        // cxx:3105-3112.
        let mut face = match &self.face {
            Some(f) => f.clone(),
            None => return false,
        };
        {
            let applied = self.context.apply(&face.0);
            if applied.is_face() {
                face = Face(applied);
                self.face = Some(face.clone());
            }
        }

        // cxx:3116-3129: `TopoDS_Iterator(myFace, false)` — oriented wires only.
        let children: Vec<TopoShape> = {
            let ts = face.0.tshape.read().unwrap();
            ts.children.clone()
        };
        let mut wire_seq: Vec<Wire> = Vec::new();
        for child in children {
            let o = child.orientation();
            if child.shape_type() != ShapeType::Wire
                || (o != Orientation::Forward && o != Orientation::Reversed)
            {
                continue;
            }
            wire_seq.push(Wire(child));
        }

        // cxx:3131-3139: only a single wire on a conical surface is checked.
        let surf = match BRepTool::face_surface(&face) {
            Some(s) => s,
            None => return false,
        };
        let cone = match surf.gp_cone() {
            Some(c) => c,
            None => return false,
        };
        if wire_seq.len() != 1 {
            return false;
        }
        let mut sole_wire = wire_seq.remove(0);

        // cxx:3146-3161: does the wire belt the cone by one period?
        let lp = match is_periodic_conical_loop(&sole_wire, &face, CONFUSION) {
            Some(v) => v,
            None => return false,
        };

        // cxx:3167-3180: the base circle the cone was built from
        // (`VIso(0.0)` -> `Geom_Circle`, whose radius is `Radius`).
        let base_r = match surf.v_iso_curve(0.0).and_then(|c| c.circle_radius()) {
            Some(r) => r,
            None => return false,
        };
        let semi_angle = cone.semi_angle();
        if semi_angle.abs() <= CONFUSION {
            return false; // cxx:3177-3180: bad surface
        }

        // cxx:3183-3184: the V parameter of the apex.
        let apex_v = -(base_r / semi_angle.sin());

        // cxx:3187: `BRepBuilderAPI_MakeVertex(aConeSurf->Apex())`.
        let builder = TopoBuilder::new();
        let apex_pnt = cone.apex();

        // cxx:3198-3201: reject when the apex V sits on or between the wire's
        // V bounds — the 2D support line would not be consistent with the wire.
        if (apex_v - lp.min_v).abs() <= CONFUSION
            || (apex_v - lp.max_v).abs() <= CONFUSION
            || (apex_v < lp.max_v && apex_v > lp.min_v)
        {
            return false;
        }

        // cxx:3205-3223: the 2D apex support line plus the wire flip that puts
        // the apex before the wire along the seam direction. Exactly one of the
        // two branches holds, the third V relationship was rejected above.
        let (line_pnt, line_dir, flip_wire) = if apex_v < lp.min_v {
            (
                GpPnt2d::new(lp.min_u, apex_v),
                GpDir2d::new(1.0, 0.0).expect("dir"),
                !lp.is_u_decrease,
            )
        } else if apex_v > lp.max_v {
            (
                GpPnt2d::new(lp.max_u, apex_v),
                GpDir2d::new(-1.0, 0.0).expect("dir"),
                lp.is_u_decrease,
            )
        } else {
            return false;
        };

        // cxx:3226-3232: `UpdateEdge` (pcurve on `myFace` with the
        // `Precision()` tolerance), the degenerate vertex twice, the
        // `Degenerated` flag and `Range(E, First, Last)` (the edge's own range,
        // which `BRep_TEdge` shares between the 3D curve and the pcurves).
        let reg = GeometryRegistry::global();
        let face_key = GeometryRegistry::shape_key(&face.0);
        let line: Arc<dyn Curve2d> = Arc::new(Geom2dLine::from_pnt_dir(line_pnt, line_dir));
        let mut apex_edge = builder.make_shape(ShapeType::Edge);
        reg.set_edge_pcurve(&apex_edge, face_key, line);
        reg.set_edge_tolerance(&apex_edge, CONFUSION);
        reg.set_degenerated(&apex_edge, true);
        // cxx:3232 `Range(E, First, Last)`: `BRep_TEdge` carries a single range
        // shared by the 3D curve and every CurveOnSurface representation, so
        // the apex edge's pcurve range becomes `[0, dU]` as well.
        reg.set_edge_range(&apex_edge, 0.0, (lp.max_u - lp.min_u).abs());
        reg.set_pcurve_range(&apex_edge, face_key, 0.0, (lp.max_u - lp.min_u).abs());
        let mut apex = builder.make_vertex(apex_pnt, CONFUSION);
        apex.0.set_orientation(Orientation::Forward);
        builder.add(&mut apex_edge, &apex.0);
        apex.0.set_orientation(Orientation::Reversed);
        builder.add(&mut apex_edge, &apex.0);
        let apex_wire = builder.make_wire(&[Edge(apex_edge)]);

        // cxx:3239-3253: the new face over the old surface, carrying the old
        // face's location and orientation (`EmptyCopied`) with the two wires.
        if flip_wire {
            sole_wire.0.reverse();
        }
        let mut new_face = builder.make_face(surf.clone(), &[sole_wire, apex_wire]);
        new_face.0.set_location(face.0.location());
        new_face.0.set_orientation(face.0.orientation());

        // cxx:3255-3257.
        self.result = Some(new_face.0.clone());
        self.context.replace(&face.0, &new_face.0);
        true
    }

    /// `ShapeFix_Face::FixMissingSeam` (`ShapeFix_Face.cxx:1722-2330`), ported
    /// through the setup and wire-pair selection (`cxx:1722-1898`).
    ///
    /// UNPORTED: the seam construction and the `w2 != null` merge
    /// (`cxx:1899-2330`) — port them before wiring `ShapeFix_Face::Perform`.
    pub fn fix_missing_seam(&mut self) -> bool {
        let surf = match &self.surf {
            Some(s) => s.clone(),
            None => return false, // cxx:1724-1727
        };
        let mut face = match &self.face {
            Some(f) => f.clone(),
            None => return false,
        };
        // cxx:1729-1735.
        let uclosed = crate::pcurve_full::sa_is_u_closed(surf.as_ref(), CONFUSION);
        let vclosed = crate::pcurve_full::sa_is_v_closed(surf.as_ref(), CONFUSION);
        if !uclosed && !vclosed {
            return false;
        }

        // cxx:1737-1741: `myFace = TopoDS::Face(Context()->Apply(myFace))`.
        // This is what hands the face built by `FixPeriodicDegenerated`
        // (`cxx:3108-3111` applies the context as well) to the wire-pair
        // selection below.
        let applied = self.context.apply(&face.0);
        if applied.is_face() {
            face = Face(applied);
            self.face = Some(face.clone());
        }

        // cxx:1744-1751: a BSpline surface must be U- or V-periodic.
        if surf.is_bspline_surface() && !surf.is_u_periodic() && !surf.is_v_periodic() {
            return false;
        }

        // cxx:1753-1756.
        let (mut suf, mut sul) = surf.u_range();
        let (mut svf, mut svl) = surf.v_range();
        let (f_u1, f_u2, f_v1, f_v2) = brep_tools::uv_bounds(&face);

        // cxx:1758-1802: replace infinite surface bounds with the face's.
        if Precision::is_infinite(suf) || Precision::is_infinite(sul) {
            if Precision::is_infinite(suf) {
                suf = f_u1;
            }
            if Precision::is_infinite(sul) {
                sul = f_u2;
            }
            if (sul - suf).abs() < PCONFUSION {
                if Precision::is_infinite(suf) {
                    suf -= 1000.0;
                } else {
                    sul += 1000.0;
                }
            }
        }
        if Precision::is_infinite(svf) || Precision::is_infinite(svl) {
            if Precision::is_infinite(svf) {
                svf = f_v1;
            }
            if Precision::is_infinite(svl) {
                svl = f_v2;
            }
            if (svl - svf).abs() < PCONFUSION {
                if Precision::is_infinite(svf) {
                    svf -= 1000.0;
                } else {
                    svl += 1000.0;
                }
            }
        }

        // cxx:1804-1805.
        let u_range = (sul - suf).abs().min(INFINITE);
        let v_range = (svl - svf).abs().min(INFINITE);

        // cxx:1810-1822: `TopoDS_Iterator(myFace, false)` — oriented wires go
        // to `ws`, every other direct child (non-wire or non-oriented) to
        // `aSeqNonManif`.
        let children: Vec<TopoShape> = {
            let ts = face.0.tshape.read().unwrap();
            ts.children.clone()
        };
        let mut ws: Vec<Wire> = Vec::new();
        let mut non_manifold: Vec<TopoShape> = Vec::new();
        for child in children {
            let o = child.orientation();
            if child.shape_type() == ShapeType::Wire
                && (o == Orientation::Forward || o == Orientation::Reversed)
            {
                ws.push(Wire(child));
            } else {
                non_manifold.push(child);
            }
        }

        // cxx:1824-1881: select the wire pair.
        let mut w1: Option<Wire> = None;
        let mut w2: Option<Wire> = None;
        let mut ismodeu = 0i32;
        let mut ismodev = 0i32;
        let mut isdeg1 = 0usize;
        let mut isdeg2 = 0usize;
        let mut i = 0usize;
        while i < ws.len() {
            let wire = ws[i].clone();
            let (isuopen, isvopen, isdeg) = match check_wire(&wire, &face, u_range, v_range) {
                Some(v) => v,
                None => {
                    i += 1;
                    continue; // cxx:1831-1834
                }
            };
            if w1.is_none() {
                ismodeu = isuopen;
                ismodev = isvopen;
                isdeg1 = if isdeg { i + 1 } else { 0 };
                w1 = Some(wire);
            } else if w2.is_none() {
                if ismodeu == -isuopen && ismodev == -isvopen {
                    isdeg2 = if isdeg { i + 1 } else { 0 };
                    w2 = Some(wire);
                } else if ismodeu == isuopen && ismodev == isvopen {
                    // cxx:1849-1865.
                    isdeg2 = usize::from(isdeg);
                    w2 = Some(wire);
                    if isdeg1 != 0 {
                        if let Some(w) = w1.as_mut() {
                            w.0.reverse();
                        }
                        ismodeu = -ismodeu;
                        ismodev = -ismodev;
                    } else if let Some(w) = w2.as_mut() {
                        w.0.reverse();
                    }
                }
            } else if isdeg || isdeg1 != 0 || isdeg2 != 0 {
                // cxx:1868-1879.
                let rem = if isdeg {
                    i + 1
                } else if isdeg2 != 0 {
                    isdeg2
                } else {
                    isdeg1
                };
                ws.remove(rem - 1);
                w1 = None;
                w2 = None;
                i = 0;
                continue;
            }
            i += 1;
        }

        // cxx:1883-1893: a torus with MajorRadius < MinorRadius counts as
        // degenerated, unless a second wire exists.
        let torus = surf.gp_torus();
        let mut an_is_degenerated_tor = false;
        if let Some(t) = &torus {
            an_is_degenerated_tor = t.major_radius < t.minor_radius;
        }
        if an_is_degenerated_tor && w2.is_some() {
            an_is_degenerated_tor = false;
        }

        // cxx:1895-1898.
        if w1.is_none() {
            return false;
        }
        if w2.is_none() {
            // cxx:1899-1992: add a degenerated edge for spheres and BSpline
            // cone-like surfaces, then handle them as a regular pair.
            let mut builder_pt = GpPnt2d::new(0.0, 0.0);
            let dir;
            let a_range;
            if ismodeu != 0 && an_is_degenerated_tor {
                let t = torus.as_ref().expect("torus");
                let a_ra = t.major_radius;
                let a_ri = t.minor_radius;
                let a_phi = (-a_ra / a_ri).acos();
                builder_pt = GpPnt2d::new(
                    0.0,
                    if ismodeu > 0 { std::f64::consts::PI + a_phi } else { a_phi },
                );
                dir = GpDir2d::new(-(ismodeu as f64), 0.0).expect("dir");
                a_range = 2.0 * std::f64::consts::PI;
            } else if ismodeu != 0 && surf.gp_sphere().is_some() {
                builder_pt = GpPnt2d::new(
                    if ismodeu < 0 { 0.0 } else { 2.0 * std::f64::consts::PI },
                    ismodeu as f64 * 0.5 * std::f64::consts::PI,
                );
                dir = GpDir2d::new(-(ismodeu as f64), 0.0).expect("dir");
                a_range = 2.0 * std::f64::consts::PI;
            } else if ismodev != 0 && surf.is_bspline_surface() {
                let u_coord;
                if surf.d0(suf, svf).distance(&surf.d0(suf, (svf + svl) / 2.0)) < CONFUSION {
                    u_coord = suf;
                } else if surf.d0(sul, svf).distance(&surf.d0(sul, (svf + svl) / 2.0)) < CONFUSION {
                    u_coord = sul;
                } else {
                    return false; // cxx:1940-1943
                }
                builder_pt = GpPnt2d::new(u_coord, if ismodev < 0 { 0.0 } else { v_range });
                dir = GpDir2d::new(0.0, -(ismodev as f64)).expect("dir");
                a_range = v_range;
            } else if ismodeu != 0 && surf.is_bspline_surface() {
                let v_coord;
                if surf.d0(suf, svf).distance(&surf.d0((suf + sul) / 2.0, svf)) < CONFUSION {
                    v_coord = svf;
                } else if surf.d0(sul, svl).distance(&surf.d0((suf + sul) / 2.0, svl)) < CONFUSION {
                    v_coord = svl;
                } else {
                    return false; // cxx:1962-1965
                }
                builder_pt = GpPnt2d::new(if ismodeu < 0 { 0.0 } else { u_range }, v_coord);
                dir = GpDir2d::new(-(ismodeu as f64), 0.0).expect("dir");
                a_range = u_range;
            } else {
                return false; // cxx:1972-1975
            }

            // cxx:1977-1991.
            let builder = TopoBuilder::new();
            let line: Arc<dyn Curve2d> = Arc::new(Geom2dLine::from_pnt_dir(builder_pt, dir));
            let reg = GeometryRegistry::global();
            let face_key = GeometryRegistry::shape_key(&face.0);
            let mut edge = builder.make_shape(ShapeType::Edge);
            reg.set_degenerated(&edge, true);
            reg.set_edge_pcurve(&edge, face_key, line);
            reg.set_pcurve_range(&edge, face_key, 0.0, a_range);
            reg.set_edge_tolerance(&edge, CONFUSION);
            let mut v = builder.make_vertex(surf.d0(builder_pt.x(), builder_pt.y()), CONFUSION);
            v.0.set_orientation(Orientation::Forward);
            builder.add(&mut edge, &v.0);
            v.0.set_orientation(Orientation::Reversed);
            builder.add(&mut edge, &v.0);
            let new_w2 = builder.make_wire(&[Edge(edge)]);
            ws.push(new_w2.clone());
            w2 = Some(new_w2);
        }

        // cxx:1994-2006: UV bounds of the two wires on the face's surface.
        let wire_uv_bounds = |w: &Wire| -> (f64, f64, f64, f64) {
            let mut b = BndBox2d::new();
            add_uv_bounds_on_wire(&face, w, &mut b);
            match b.get() {
                Some((u0, v0, u1, v1)) => (u0, u1, v0, v1),
                None => (0.0, 0.0, 0.0, 0.0),
            }
        };
        let (w1, w2) = (w1.unwrap(), w2.unwrap());
        let m1 = wire_uv_bounds(&w1);
        let m2 = wire_uv_bounds(&w2);
        let coord = if ismodeu != 0 { 1usize } else { 0usize };
        let isneg = if ismodeu != 0 { ismodeu } else { -ismodev };
        let period = if ismodeu != 0 { u_range } else { v_range };
        // `m1[coord]` / `m2[coord]` in `ShapeFix_Face.cxx:1996-2127`: the arrays
        // are filled by `ShapeAnalysis::GetFaceUVBounds(F, UMin, UMax, VMin,
        // VMax)` (cxx:268-273), so `m[..][0]` holds the U bounds and `m[..][1]`
        // the V bounds, while `coord = (ismodeu ? 1 : 0)` — `coord == 1`
        // therefore selects **V**, not U. `pick(m, 1)` returns the V pair.
        let pick = |m: &(f64, f64, f64, f64), c: usize| -> (f64, f64) {
            if c == 1 { (m.2, m.3) } else { (m.0, m.1) }
        };
        let mut w1 = w1;
        let mut w2 = w2;
        if !vclosed || !uclosed || an_is_degenerated_tor {
            let (a0, a1) = pick(&m1, coord);
            let (b0, b1) = pick(&m2, coord);
            let delta_other = 0.5 * (b0 + b1) - 0.5 * (a0 + a1);
            if delta_other * (isneg as f64) < 0.0 {
                w1.0.reverse(); // cxx:2016-2020
                w2.0.reverse();
            }
        }

        // UNPORTED (cxx:2029-2032): `ShapeFix_Wire::FixReorder` — the port's
        // `fix_reorder_wire` only reports the 3D order, it does not rewrite the
        // stored edge list, so `w11`/`w21` keep the original order.
        let _ = fix_reorder_wire(&w1, &face);
        let w11 = w1.clone();
        let w21 = w2.clone();

        // cxx:2036-2066: rebuild the face, replacing w1/w2 with w11/w21 and
        // orienting the other wires as holes.
        let builder = TopoBuilder::new();
        let mut tmp_wires: Vec<Wire> = Vec::new();
        for w in &ws {
            let mut wire = w.clone();
            if crate::topo_tools_full::is_same(&wire.0, &w1.0) {
                wire = w11.clone();
            } else if crate::topo_tools_full::is_same(&wire.0, &w2.0) {
                wire = w21.clone();
            } else {
                // cxx:2052-2062.
                let cur_face = builder.make_face(surf.clone(), &[wire.clone()]);
                let mut cur_face = cur_face;
                cur_face.0.set_orientation(face.0.orientation());
                if is_outer_bound(&cur_face) {
                    wire.0.reverse();
                }
            }
            tmp_wires.push(wire);
        }
        let mut tmp_f = builder.make_face(surf.clone(), &tmp_wires);
        tmp_f.0.set_orientation(face.0.orientation());

        // cxx:2068-2136: torus-like FixShifted along the missing seam direction.
        let m1 = m1;
        let mut m1 = m1;
        let m2 = m2;
        let mut uf2 = suf;
        let mut vf2 = svf;
        if uclosed && vclosed && !an_is_degenerated_tor {
            let shiftw2 = adjust_by_period(
                0.5 * (pick(&m2, coord).0 + pick(&m2, coord).1),
                0.5 * (pick(&m1, coord).0 + pick(&m1, coord).1
                    + (isneg as f64) * (period + PCONFUSION)),
                period,
            );
            let (m2c0, m2c1) = pick(&m2, coord);
            let (m1c0, m1c1) = pick(&m1, coord);
            // cxx:2077-2078: `m1[coord][0] = min(...)`, `m1[coord][1] = max(...)`
            // — the update lands on whichever pair `coord` selects (V when
            // `ismodeu`), not on U.
            if coord == 1 {
                m1.2 = m1c0.min(m2c0 + shiftw2);
                m1.3 = m1c1.max(m2c1 + shiftw2);
            } else {
                m1.0 = m1c0.min(m2c0 + shiftw2);
                m1.1 = m1c1.max(m2c1 + shiftw2);
            }
            let tmp_face_wires = wires_of_face(&tmp_f);
            for w in &tmp_face_wires {
                if crate::topo_tools_full::is_same(&w.0, &w11.0) {
                    continue;
                }
                let shift = if crate::topo_tools_full::is_same(&w.0, &w21.0) {
                    shiftw2
                } else {
                    let mw = wire_uv_bounds(w);
                    adjust_by_period(
                        0.5 * (pick(&mw, coord).0 + pick(&mw, coord).1),
                        0.5 * (pick(&m1, coord).0 + pick(&m1, coord).1),
                        period,
                    )
                };
                if shift != 0.0 {
                    // cxx:2107 `V.SetCoord(coord + 1, shift)`: `coord == 1`
                    // writes gp_Vec2d's Y (V) component, `coord == 0` its X (U).
                    let vshift = if coord == 1 {
                        GpVec2d::new(0.0, shift)
                    } else {
                        GpVec2d::new(shift, 0.0)
                    };
                    let tkey = GeometryRegistry::shape_key(&tmp_f.0);
                    for e in edges_of_wire(w) {
                        if let Some((c2d, _a, _b)) = boptools_2d::curve_on_surface_range(&e, &tmp_f) {
                            let mut nc = c2d.clone_dyn();
                            let mut tr = occt_core::gp::GpTrsf2d::identity();
                            tr.set_translation_vec(&vshift);
                            nc.transform(&tr);
                            GeometryRegistry::global().set_edge_pcurve(&e.0, tkey, Arc::from(nc));
                        }
                    }
                }
            }
            // cxx:2122-2135.
            if pick(&m1, coord).1 - pick(&m1, coord).0 <= period {
                let other = 0.5 * (pick(&m1, coord).0 + pick(&m1, coord).1 - period);
                if ismodeu != 0 {
                    vf2 = other;
                } else {
                    uf2 = other;
                }
            }
        }

        // cxx:2138-2234: find the best place to insert the seam.
        let mut found_u = 0i32;
        let mut found_v = 0i32;
        let wd1_edges = edges_of_wire(&w11);
        let wd2_edges = edges_of_wire(&w21);
        let nb1 = wd1_edges.len();
        let nb2 = wd2_edges.len();
        let mut i1 = 1usize;
        while i1 <= nb1 + nb2 {
            let edge1 = if i1 <= nb1 {
                wd1_edges[i1 - 1].clone()
            } else {
                wd2_edges[i1 - nb1 - 1].clone()
            };
            let (c2d, f, l) = match boptools_2d::curve_on_surface_range(&edge1, &tmp_f) {
                Some(v) => v,
                None => return false, // cxx:2149-2152
            };
            let mut pos1 = c2d.d0(l);
            let mut skip_u = !uclosed;
            if uclosed && ismodeu != 0 {
                pos1 = GpPnt2d::new(pos1.x() + adjust_by_period(pos1.x(), suf, u_range), pos1.y());
                if found_u == 2 && pos1.x().abs() > uf2.abs() {
                    skip_u = true;
                } else if found_u == 0 || (found_u == 1 && pos1.x().abs() < uf2.abs()) {
                    found_u = 1;
                    uf2 = pos1.x();
                }
            }
            let mut skip_v = !vclosed;
            if vclosed && ismodeu == 0 {
                pos1 = GpPnt2d::new(pos1.x(), pos1.y() + adjust_by_period(pos1.y(), svf, v_range));
                if found_v == 2 && pos1.y().abs() > vf2.abs() {
                    skip_v = true;
                } else if found_v == 0 || (found_v == 1 && pos1.y().abs() < vf2.abs()) {
                    found_v = 1;
                    vf2 = pos1.y();
                }
            }
            if skip_u && skip_v {
                if i1 <= nb1 {
                    i1 += 1;
                    continue;
                } else {
                    break;
                }
            }
            for i2 in 1..=nb2 {
                if i1 > nb1 {
                    break;
                }
                let edge2 = wd2_edges[i2 - 1].clone();
                let (c2d2, f2, _l2) = match boptools_2d::curve_on_surface_range(&edge2, &tmp_f) {
                    Some(v) => v,
                    None => return false, // cxx:2198-2201
                };
                let mut pos2 = c2d2.d0(f2);
                if uclosed && ismodeu != 0 {
                    pos2 = GpPnt2d::new(pos2.x() + adjust_by_period(pos2.x(), pos1.x(), u_range), pos2.y());
                    if (pos2.x() - pos1.x()).abs() < PCONFUSION
                        && (found_u != 2 || pos1.x().abs() < uf2.abs())
                    {
                        found_u = 2;
                        uf2 = pos1.x();
                    }
                }
                if vclosed && ismodeu == 0 {
                    pos2 = GpPnt2d::new(pos2.x(), pos2.y() + adjust_by_period(pos2.y(), pos1.y(), v_range));
                    if (pos2.y() - pos1.y()).abs() < PCONFUSION
                        && (found_v != 2 || pos1.y().abs() < vf2.abs())
                    {
                        found_v = 2;
                        vf2 = pos1.y();
                    }
                }
            }
            i1 += 1;
        }

        // cxx:2226-2234.
        if uf2 < suf || uf2 > sul {
            uf2 += adjust_to_period(uf2, suf, suf + u_range);
        }
        if vf2 < svf || vf2 > svl {
            vf2 += adjust_to_period(vf2, svf, svf + v_range);
        }

        // cxx:2236-2261: fictive grid + ComposeShell inserts the seam.
        let rts: Arc<dyn Surface> = Arc::new(GeomRectangularTrimmedSurface::uv(
            surf.clone(),
            uf2,
            uf2 + u_range,
            vf2,
            vf2 + v_range,
        ));
        let grid = CompositeSurface::with_grid(vec![vec![rts.clone()]], Parametrisation::Natural);
        // cxx:2246-2250: re-add the non-manifold children.
        for s in &non_manifold {
            builder.add(&mut tmp_f.0, s);
        }
        let mut comp = ComposeShell::new();
        comp.init(grid, &tmp_f, CONFUSION); // cxx:2252-2253
        comp.set_closed_mode(true); // cxx:2258
        comp.set_context(self.context.clone()); // cxx:2259
        comp.set_max_tolerance(self.max_tol); // cxx:2260
        comp.perform(); // cxx:2261
        // cxx:2263-2264: reset mySurf to the trimmed surface.
        self.surf = Some(rts);
        // cxx:2266 myResult = CompShell.Result().
        let mut result = comp.result().cloned();
        // cxx:2270-2322: drop the small wires / faces ComposeShell can generate.
        //  * cxx:2273-2300: every wire of every result face goes through
        //    ShapeFix_Wire::FixSmall(true, Precision()); a wire that loses all
        //    its edges is discarded and a face that loses all its wires is
        //    removed (Context()->Remove).
        //  * cxx:2303-2321: when the result keeps more than one face
        //    (nbFaces > 1), FixSmallAreaWire(true) (cxx:2331-2394) discards
        //    every wire for which ShapeAnalysis_Wire::CheckSmallArea (cxx:2004)
        //    reports a null area; a face whose wires all disappear is removed.
        // UNPORTED: Context()->Apply (cxx:2302/2322) is expressed by rebuilding
        // the result here; BRepTools::Update (cxx:2320) is a no-op.
        if let Some(res) = &result {
            if res.shape_type() != ShapeType::Face {
                let builder = TopoBuilder::new();
                // First round: FixSmall.
                let mut round1: Vec<(Face, Vec<Wire>)> = Vec::new();
                for f in crate::topo_tools_full::faces_of(res) {
                    let mut kept: Vec<Wire> = Vec::new();
                    for w in wires_of_face(&f) {
                        let mut w2 = w.clone();
                        crate::shhealing::fix_small_all(&mut w2, &f, CONFUSION, CONFUSION, false);
                        if !edges_of_wire(&w2).is_empty() {
                            kept.push(w2);
                        }
                    }
                    if !kept.is_empty() {
                        round1.push((f, kept));
                    }
                }
                let nb_faces = round1.len();
                let mut kept_faces: Vec<Face> = Vec::new();
                for (f, wires) in round1 {
                    let mut kept_wires = wires;
                    if nb_faces > 1 {
                        kept_wires.retain(|w| !crate::shhealing::check_small_area(w, &f));
                    }
                    if kept_wires.is_empty() {
                        continue; // cxx:2377-2385 / 2294-2298: remove the face
                    }
                    let orig_n = wires_of_face(&f).len();
                    if kept_wires.len() == orig_n {
                        kept_faces.push(f);
                    } else if let Some(s) = BRepTool::face_surface(&f) {
                        let mut nf = builder.make_face(s, &kept_wires);
                        nf.0.set_orientation(f.0.orientation());
                        kept_faces.push(nf);
                    }
                }
                result = match kept_faces.len() {
                    0 => None,
                    1 => Some(kept_faces.remove(0).0.clone()),
                    _ => Some(builder.make_shell(&kept_faces).0),
                };
            }
        }
        self.result = result;
        // cxx:2268: Context()->Replace(myFace, myResult).
        if let Some(res) = &self.result {
            self.context.replace(&face.0, res);
        }
        true
    }
}