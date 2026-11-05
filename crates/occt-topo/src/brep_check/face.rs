//! `BRepCheck_Face` — face validity checks.
//!
//! Source: `ModelingAlgorithms/TKTopAlgo/BRepCheck/BRepCheck_Face.cxx`. Ported
//! here: `ClassifyWires` (`cxx:302-434`), `OrientationOfWires` (`cxx:436-557`),
//! `IntersectWires` (`cxx:169-296`), the static `Intersect`
//! (`cxx:626-797`), `IsInside` (`cxx:802-875`) and `CheckThin` (`cxx:877-952`).
//! `OrientationOfWires` is ported but not wired into [`face_check`] yet — see
//! [`ORIENTATION_OF_WIRES`].
//!
//! `BRepTopAdaptor_FClass2d` is [`crate::brep_top_adaptor_fclass2d`].

use std::sync::Arc;

use occt_core::bnd::BndBox2d;
use occt_core::gp::{GpPnt, GpPnt2d, GpVec2d};
use occt_core::intres2d::IntRes2dDomain;
use occt_core::precision::{Precision, CONFUSION, PCONFUSION};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::geom2d_int::Geom2dIntGInter;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::brep_top_adaptor_fclass2d::BrepTopAdaptorFClass2d;
use crate::fclass2d::FaceState;
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{is_same, wires_of_face};

use super::{wire_check, Status};

/// Whether [`face_check`] runs `BRepCheck_Face::OrientationOfWires`
/// (`BRepCheck_Face.cxx:436-557`). Off, because the port reports
/// `BadOrientationOfSubshape` where OCCT 8.0.0's `BRepCheck_Analyzer` reports
/// `NoError`: on the record-parity set two faces of `data/rev.step` come back
/// `BadOrientationOfSubshape`, which turns `VALID` to 0 for `rev`, `screw`,
/// `OffsetPlaneHoleEdge`, `linkrods` and `Shape-1` in the STEP-read path — all
/// five are `VALID 1` after an IGES round trip here and on the OCCT side in both
/// modes. A false `BadOrientationOfSubshape` rejects valid geometry, so the
/// branch stays unwired while it is unverified.
const ORIENTATION_OF_WIRES: bool = false;

pub(crate) fn face_check(f: &Face) -> Vec<Status> {
    if BRepTool::face_surface(f).is_none() {
        return vec![Status::NoSurface];
    }
    let wires = wires_of_face(f);
    if wires.is_empty() {
        return vec![Status::InvalidWire];
    }
    for w in &wires {
        if wire_check(w).iter().any(|s| *s != Status::NoError) {
            return vec![Status::InvalidWire];
        }
    }
    // `BRepCheck_Analyzer::Perform` runs `BRepCheck_Face::OrientationOfWires`
    // for a face whose edges and wires are all clean
    // (`BRepCheck_Analyzer.cxx:269-281`).
    if ORIENTATION_OF_WIRES {
        let st = orientation_of_wires(f);
        if st != Status::NoError {
            return vec![st];
        }
    }
    vec![Status::NoError]
}

/// `DetectWires` half of `BRepCheck_Face::IntersectWires`
/// (`BRepCheck_Face.cxx:169-296`).
///
/// Ported: the wire map (`cxx:197-216`, duplicate wire -> `RedundantWire`), the
/// 2D boxes (`cxx:219-247`) and the pairwise `Intersect` test
/// (`cxx:249-294` -> `BRepCheck_Face.cxx:626-797`).
fn intersect_wires(f: &Face) -> Status {
    let wires = wires_of_face(f);
    // `cxx:197-216`: `myMapImb` is keyed by `TopTools_ShapeMapHasher`, so a
    // wire met twice is `RedundantWire`.
    for i in 0..wires.len() {
        for j in (i + 1)..wires.len() {
            if crate::topo_tools::is_same(&wires[i].0, &wires[j].0) {
                return Status::RedundantWire;
            }
        }
    }
    // `cxx:224-247`: per-edge and per-wire 2D boxes on the face's surface.
    let mut edge_boxes: Vec<(TopoShape, BndBox2d)> = Vec::new();
    let mut wire_boxes: Vec<BndBox2d> = Vec::with_capacity(wires.len());
    for w in &wires {
        let mut box_w = BndBox2d::new();
        for e in crate::topo_tools::edges_of_wire(w) {
            let Some((pc, mut first, mut last)) = pcurve_on_face(&e, f) else {
                // `cxx:230`: `Geom2dAdaptor_Curve::Load` of a null handle
                // throws, which the analyzer turns into a failure
                // (`BRepCheck_Analyzer.cxx:284-302`).
                return Status::CheckFail;
            };
            // `cxx:232-239`: clamp the representation range into the adaptor.
            if pc.first_parameter() > first {
                first = pc.first_parameter();
            }
            if pc.last_parameter() < last {
                last = pc.last_parameter();
            }
            let mut box_e = BndBox2d::new();
            crate::bnd_lib_add2d::add_geom2d_range(pc.as_ref(), first, last, 0.0, &mut box_e);
            box_w.add_box(&box_e);
            edge_boxes.push((e.0.clone(), box_e));
        }
        wire_boxes.push(box_w);
    }

    // `cxx:249-294`: every wire pair, the outer box test first.
    for i in 0..wires.len().saturating_sub(1) {
        let box1 = wire_boxes[i];
        for j in (i + 1)..wires.len() {
            if !box1.is_void() && !wire_boxes[j].is_void() && box1.is_out_box(&wire_boxes[j]) {
                continue;
            }
            match wires_intersect(&wires[i], &wires[j], f, &edge_boxes) {
                Ok(true) => return Status::IntersectingWires,
                Ok(false) => {}
                Err(()) => return Status::CheckFail,
            }
        }
    }
    Status::NoError
}

/// `myMapImb` entry (`BRepCheck_Face.cxx:197-216` binds every wire to an empty
/// list; `cxx:365-368` re-binds a hole wire with the `REVERSED` orientation,
/// which `TopTools_ShapeMapHasher` still resolves to the same key because it
/// compares `IsSame`).
struct WireEntry {
    /// The map key.
    wire: Wire,
    /// The key's orientation: `Forward` normally, `Reversed` when
    /// `BRepTopAdaptor_FClass2d::PerformInfinitePoint` said the wire defines a
    /// hole (`cxx:362-368`).
    orientation: Orientation,
    /// The wires this one contains (`cxx:377-381`).
    contains: Vec<Wire>,
}

/// `BRepCheck_Face::ClassifyWires` (`BRepCheck_Face.cxx:302-434`).
///
/// Returns the status together with the filled `myMapImb`, which
/// [`orientation_of_wires`] reads. The throwaway `BRep_Builder` face of
/// `cxx:355-360` is [`BrepTopAdaptorFClass2d::new_with_wires`].
fn classify_wires(f: &Face) -> (Status, Vec<WireEntry>) {
    let st = intersect_wires(f);
    if st != Status::NoError {
        return (st, Vec::new());
    }

    let face = Face(f.0.oriented(Orientation::Forward));
    let wires = wires_of_face(&face);
    let nbwire = wires.len();
    if nbwire < 1 {
        return (Status::NoError, Vec::new());
    }

    let mut map: Vec<WireEntry> = wires
        .iter()
        .map(|w| WireEntry {
            wire: w.clone(),
            orientation: w.0.orientation(),
            contains: Vec::new(),
        })
        .collect();

    for wir1 in &wires {
        // `cxx:355-360`: a face made of this wire alone, on the same surface.
        let fclass2d = BrepTopAdaptorFClass2d::new_with_wires(
            &face,
            std::slice::from_ref(wir1),
            PCONFUSION,
        );
        let wire_bien_oriente = fclass2d.perform_infinite_point() != FaceState::Out;
        if wire_bien_oriente {
            // `cxx:365-368`: the given wire defines a hole.
            if let Some(e) = map.iter_mut().find(|e| is_same(&e.wire.0, &wir1.0)) {
                e.orientation = e.orientation.reversed();
            }
        }
        for wir2 in &wires {
            if is_same(&wir2.0, &wir1.0) {
                continue;
            }
            if is_inside(wir2, wire_bien_oriente, &fclass2d, &face) {
                // `cxx:379`.
                if let Some(e) = map.iter_mut().find(|e| is_same(&e.wire.0, &wir1.0)) {
                    e.contains.push(wir2.clone());
                }
            }
        }
    }

    // `cxx:386-411`: at most one wire may contain anything, and it must contain
    // every other wire.
    let mut nb_with_content = 0usize;
    for e in &map {
        if !e.contains.is_empty() {
            nb_with_content += 1;
            if nb_with_content > 1 {
                return (Status::InvalidImbricationOfWires, map);
            }
        }
    }
    if let Some(wext) = map.iter().find(|e| !e.contains.is_empty()) {
        // `cxx:414-417`.
        if wext.contains.len() != nbwire - 1 {
            return (Status::InvalidImbricationOfWires, map);
        }
    }

    (Status::NoError, map)
}

/// `BRepCheck_Face::OrientationOfWires` (`BRepCheck_Face.cxx:436-557`).
///
/// Not wired into [`face_check`] — see [`ORIENTATION_OF_WIRES`].
fn orientation_of_wires(f: &Face) -> Status {
    // WARNING: it is assumed that the edges of a wire are correctly oriented.
    let infinite = f.0.infinite();
    let (st, map) = classify_wires(f);
    if st != Status::NoError {
        return st;
    }

    let nbwire = map.len();
    // `cxx:472-491`.
    let wext: Option<&WireEntry> = if nbwire == 1 {
        if infinite {
            None
        } else {
            map.first()
        }
    } else {
        map.iter().find(|e| !e.contains.is_empty())
    };

    if wext.is_none() && !infinite {
        if nbwire > 0 {
            return Status::InvalidImbricationOfWires;
        }
        return Status::NoError;
    }

    let face_forward = Face(f.0.oriented(Orientation::Forward));
    for wir in wires_of_face(&face_forward) {
        let or = wir.0.orientation();
        let key_orientation = map
            .iter()
            .find(|e| is_same(&e.wire.0, &wir.0))
            .map(|e| e.orientation);

        let is_wext = wext
            .map(|w| is_same(&w.wire.0, &wir.0))
            .unwrap_or(false);
        if is_wext {
            let wext_or = wext.map(|w| w.orientation).unwrap_or(Orientation::Forward);
            if or != wext_or {
                // The exterior wire defines a hole.
                if check_thin(&wir, &face_forward) {
                    return Status::NoError;
                }
                return Status::BadOrientationOfSubshape;
            }
        } else {
            // `cxx:531-548`: a wire whose orientation matches its map key does
            // not define a hole.
            if key_orientation == Some(or) {
                return Status::BadOrientationOfSubshape;
            }
        }
    }

    Status::NoError
}

/// `IsInside(theWire, WireBienOriente, FClass2d, theFace)`
/// (`BRepCheck_Face.cxx:802-875`): one point of the wire's first usable edge,
/// classified in the UV of the throwaway face.
fn is_inside(
    wire: &Wire,
    wire_bien_oriente: bool,
    fclass2d: &BrepTopAdaptorFClass2d,
    face: &Face,
) -> bool {
    for e in crate::topo_tools::edges_of_wire(wire) {
        let Some((pc2d, first, last)) = pcurve_on_face(&e, face) else {
            // `cxx:818-820`: a null `CurveOnSurface` would dereference null in
            // OCCT; skip instead.
            continue;
        };

        let param;
        if !Precision::is_negative_infinite(first) && !Precision::is_positive_infinite(last) {
            param = (first + last) * 0.5;
            // `cxx:821-824`: skip an edge whose parametric range is too small.
            if (param - first).abs() < PCONFUSION {
                continue;
            }
            // `cxx:827-831`: skip an edge with no 3D curve.
            let (first3d, last3d) = BRepTool::edge_parameters(&e);
            let Some(curve3d) = BRepTool::edge_curve(&e) else {
                continue;
            };
            // `cxx:834-841`: the start point is evaluated at `aFirst` (the 2D
            // range start) while the middle point uses the 3D range — kept as
            // written in OCCT.
            let p0 = curve3d.d0(first);
            let pmid = curve3d.d0((first3d + last3d) / 2.0);
            if p0.distance(&pmid) < CONFUSION {
                continue;
            }
        } else if Precision::is_negative_infinite(first) && Precision::is_positive_infinite(last) {
            param = 0.0;
        } else if Precision::is_negative_infinite(first) {
            param = last - 1.0;
        } else {
            param = first + 1.0;
        }

        let p2d: GpPnt2d = pc2d.d0(param);
        let state = fclass2d.perform(p2d, false);
        return if wire_bien_oriente {
            state == FaceState::Out
        } else {
            state == FaceState::In
        };
    }
    false
}

/// `CheckThin(w, f)` (`BRepCheck_Face.cxx:877-952`): a 2-edge wire may still be
/// a hole when the two edges run in opposite directions or in opposite
/// orientations.
fn check_thin(w: &Wire, f: &Face) -> bool {
    let edges = crate::topo_tools::edges_of_wire(w);
    if edges.len() != 2 {
        return false;
    }
    let e1 = &edges[0];
    let e2 = &edges[1];

    let (v1, v2) = crate::topo_tools::edge_vertices(e1);
    let (v3, v4) = crate::topo_tools::edge_vertices(e2);
    let (Some(v1), Some(v2), Some(v3), Some(v4)) = (v1, v2, v3, v4) else {
        return false;
    };

    if crate::topo_tools::is_same(&v1.0, &v2.0) || crate::topo_tools::is_same(&v3.0, &v4.0) {
        return false;
    }

    let mut is_first = false;
    let mut is_last = false;
    if crate::topo_tools::is_same(&v1.0, &v3.0) || crate::topo_tools::is_same(&v1.0, &v4.0) {
        is_first = true;
    }
    if crate::topo_tools::is_same(&v2.0, &v3.0) || crate::topo_tools::is_same(&v2.0, &v4.0) {
        is_last = true;
    }
    if !is_first || !is_last {
        return false;
    }

    let e1or = e1.0.orientation();
    let e2or = e2.0.orientation();

    let (Some((pc1, f1, l1)), Some((pc2, f2, l2))) = (
        pcurve_on_face(e1, f),
        pcurve_on_face(e2, f),
    ) else {
        return false;
    };

    let d1 = (l1 - f1).abs() / 100.0;
    let d2 = (l2 - f2).abs() / 100.0;
    let m1 = (l1 + f1) * 0.5;
    let m2 = (l2 + f2) * 0.5;

    let p1f = pc1.d0(m1 - d1);
    let p1l = pc1.d0(m1 + d1);
    let p2f = pc2.d0(m2 - d2);
    let p2l = pc2.d0(m2 + d2);

    let vc1 = GpVec2d::new(p1l.x() - p1f.x(), p1l.y() - p1f.y());
    let vc2 = GpVec2d::new(p2l.x() - p2f.x(), p2l.y() - p2f.y());

    // `vc1 * vc2` is `gp_Vec2d::Dot`.
    vc1.dot(&vc2) < 0.0 || e1or != e2or
}

/// `BRep_Tool::CurveOnSurface(E, F, First, Last)` (`BRep_Tool.cxx:301-315`):
/// a `REVERSED` face reverses a local copy of the edge first, and
/// `cxx:354-361` then answers `PCurve2` for a reversed seam. The stored range
/// is *not* reordered here (that is `ShapeAnalysis_Edge::PCurve`'s extra
/// `orient` step, `ShapeAnalysis_Edge.cxx:201-206`).
fn pcurve_on_face(edge: &Edge, face: &Face) -> Option<(Arc<dyn Curve2d>, f64, f64)> {
    let face_reversed = face.0.orientation() == Orientation::Reversed;
    let local_reversed = edge.0.orientation().is_reversed() != face_reversed;
    let face_key = GeometryRegistry::shape_key(&face.0);
    let reg = GeometryRegistry::global();
    let pcs = reg.edge_pcurves(&edge.0, face_key);
    if pcs.is_empty() {
        // `BRep_Tool.cxx:367-372`: no stored representation -> `CurveOnPlane`
        // projection, planar surfaces only.
        return Some(crate::boptools_2d::curve_on_surface_oriented(edge, face, false)?);
    }
    let pc = if pcs.len() >= 2 && local_reversed {
        pcs[1].clone()
    } else {
        pcs[0].clone()
    };
    let (first, last) = reg
        .pcurve_range(&edge.0, face_key)
        .unwrap_or_else(|| reg.edge_parameters(&edge.0));
    if !(first.is_finite() && last.is_finite()) {
        return None;
    }
    Some((pc, first, last))
}

/// The static `Intersect(wir1, wir2, F, theMapEdgeBox)` of
/// `BRepCheck_Face.cxx:626-797`. `Err(())` is the OCCT exception path
/// (`BRep_Tool::Parameters` throwing `Standard_NoSuchObject`,
/// `BRep_Tool.cxx:1702`), which the analyzer records as a failure.
fn wires_intersect(
    wir1: &Wire,
    wir2: &Wire,
    f: &Face,
    edge_boxes: &[(TopoShape, BndBox2d)],
) -> Result<bool, ()> {
    // `cxx:630`: `Inter2dTol`.
    const INTER_2D_TOL: f64 = 1.0e-10;

    // `cxx:637-650`: common vertices of the two wires - the non-manifold case.
    let verts1 = crate::topo_tools::vertices_of(&wir1.0);
    let mut common: Vec<Vertex> = Vec::new();
    for v2 in crate::topo_tools::vertices_of(&wir2.0) {
        if verts1.iter().any(|v1| crate::topo_tools::is_same(&v1.0, &v2.0)) {
            common.push(v2);
        }
    }

    // `cxx:654-664`: the 3D points of those vertices on the surface.
    let Some(surf) = BRepTool::face_surface(f) else {
        return Err(());
    };
    let mut pnts: Vec<GpPnt> = Vec::with_capacity(common.len());
    for v in &common {
        let p2d = BRepTool::parameters_of_vertex_on_face(&v.0, f).map_err(|_| ())?;
        pnts.push(surf.d0(p2d.x(), p2d.y()));
    }

    let edges1 = crate::topo_tools::edges_of_wire(wir1);
    let edges2 = crate::topo_tools::edges_of_wire(wir2);
    for e1 in &edges1 {
        // `cxx:675-687`.
        let Some((c1, mut first1, mut last1)) = pcurve_on_face(e1, f) else {
            return Err(());
        };
        if c1.first_parameter() > first1 {
            first1 = c1.first_parameter();
        }
        if c1.last_parameter() < last1 {
            last1 = c1.last_parameter();
        }
        let mut box1 = edge_box(edge_boxes, e1);
        if box1.is_void() {
            crate::bnd_lib_add2d::add_geom2d_range(c1.as_ref(), first1, last1, 0.0, &mut box1);
        }
        for e2 in &edges2 {
            // `cxx:699`: a shared edge is not an intersection.
            if crate::topo_tools::is_same(&e1.0, &e2.0) {
                continue;
            }
            let Some((c2, mut first2, mut last2)) = pcurve_on_face(e2, f) else {
                return Err(());
            };
            if c2.first_parameter() > first2 {
                first2 = c2.first_parameter();
            }
            if c2.last_parameter() < last2 {
                last2 = c2.last_parameter();
            }
            let mut box2 = edge_box(edge_boxes, e2);
            if box2.is_void() {
                crate::bnd_lib_add2d::add_geom2d_range(c2.as_ref(), first2, last2, 0.0, &mut box2);
            }
            if box1.is_out_box(&box2) {
                continue;
            }
            // `cxx:723-731`.
            let (pfirst1, plast1) = BRepTool::uv_points(e1, f);
            let (pfirst2, plast2) = BRepTool::uv_points(e2, f);
            let d1 = IntRes2dDomain::bounded(
                &pfirst1,
                first1,
                INTER_2D_TOL,
                &plast1,
                last1,
                INTER_2D_TOL,
            );
            let d2 = IntRes2dDomain::bounded(
                &pfirst2,
                first2,
                INTER_2D_TOL,
                &plast2,
                last2,
                INTER_2D_TOL,
            );
            let mut inter = Geom2dIntGInter::new();
            inter.perform(
                c1.as_ref(),
                &d1,
                c2.as_ref(),
                &d2,
                INTER_2D_TOL,
                INTER_2D_TOL,
            );
            if !inter.is_done() {
                return Ok(true);
            }
            // `cxx:733-766`: every overlapping segment must end on a common
            // vertex, otherwise the wires cross.
            if inter.nb_segments() > 0 {
                if pnts.is_empty() {
                    return Ok(true);
                }
                let mut nb_coinc = 0usize;
                for i in 1..=inter.nb_segments() {
                    let seg = inter.segment(i);
                    if !seg.has_first_point() || !seg.has_last_point() {
                        return Ok(true);
                    }
                    let fp = *seg.first_point().value();
                    let lp = *seg.last_point().value();
                    let first3d = surf.d0(fp.x(), fp.y());
                    let last3d = surf.d0(lp.x(), lp.y());
                    for (j, p) in pnts.iter().enumerate() {
                        // `cxx:755`.
                        let tolv = BRepTool::vertex_tolerance(&common[j]);
                        if first3d.is_equal_with_tolerance(p, tolv)
                            || last3d.is_equal_with_tolerance(p, tolv)
                        {
                            nb_coinc += 1;
                            break;
                        }
                    }
                }
                return Ok(nb_coinc != inter.nb_segments());
            }
            // `cxx:767-792`: the same test for isolated intersection points,
            // with the small extra tolerance OCCT grants the intersection
            // point itself (`cxx:781`).
            if inter.nb_points() > 0 {
                if pnts.is_empty() {
                    return Ok(true);
                }
                let mut nb_coinc = 0usize;
                for i in 1..=inter.nb_points() {
                    let p2d = *inter.point(i).value();
                    let p = surf.d0(p2d.x(), p2d.y());
                    for (j, q) in pnts.iter().enumerate() {
                        let tolv = BRepTool::vertex_tolerance(&common[j]) + 1.0e-8;
                        if p.square_distance(q) <= tolv * tolv {
                            nb_coinc += 1;
                            break;
                        }
                    }
                }
                return Ok(nb_coinc != inter.nb_points());
            }
        }
    }
    Ok(false)
}

/// `theMapEdgeBox(edg)` lookup of `BRepCheck_Face.cxx:687-692`.
fn edge_box(edge_boxes: &[(TopoShape, BndBox2d)], e: &Edge) -> BndBox2d {
    for (s, b) in edge_boxes {
        if crate::topo_tools::is_same(s, &e.0) {
            return *b;
        }
    }
    BndBox2d::new()
}
