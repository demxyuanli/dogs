//! Wire and face helper queries: the wire edge order, finite UV bounds and the per-edge reversal data. Source: `ShapeExtend_WireData`, `ShapeAnalysis_Edge`, `ShapeFix_Wire` and `BRepToIGES_BRWire`.
use super::prelude::*;
use super::*;

/// The `ShapeExtend_WireData` edge list of a wire (`ShapeExtend_WireData.cxx:72-148`
/// with the default `chained = true`, `theManifold = true`): `TopoDS_Iterator`'s
/// composed orientation, and — when the wire itself is REVERSED — every edge
/// `Prepend`ed instead of appended (`cxx:114-121`), so the list comes out in the
/// wire's own traversal order. `INTERNAL`/`EXTERNAL` edges go to
/// `myNonmanifoldEdges` (`cxx:84-89`), which `TransferWire` never reads.
pub(super) fn wire_data_edges(wire: &crate::shape::Wire) -> Vec<Edge> {
    let mut stored: Vec<Edge> = edges_of_wire(wire)
        .into_iter()
        .filter(|e| {
            let o = e.0.orientation();
            o == crate::abs::Orientation::Forward || o == crate::abs::Orientation::Reversed
        })
        .collect();
    if wire.0.orientation() == crate::abs::Orientation::Reversed {
        stored.reverse();
    }
    stored
}

/// `ShapeAnalysis_Edge::FirstVertex(edge, CumOri = true)`
/// (`ShapeAnalysis_Edge.cxx:230-238`): the FORWARD vertex of the edge *as
/// traversed*, i.e. the stored REVERSED one on a REVERSED edge.
pub(super) fn wire_edge_first_vertex(e: &Edge) -> Option<Vertex> {
    let (f, l) = edge_vertices(e);
    if e.0.orientation() == crate::abs::Orientation::Reversed {
        l
    } else {
        f
    }
}

/// `ShapeAnalysis_Edge::LastVertex(edge, CumOri = true)`.
pub(super) fn wire_edge_last_vertex(e: &Edge) -> Option<Vertex> {
    let (f, l) = edge_vertices(e);
    if e.0.orientation() == crate::abs::Orientation::Reversed {
        f
    } else {
        l
    }
}

/// `BRepToIGES_BRWire::TransferWire` (`BRepToIGES_BRWire.cxx:704-726`) orders the
/// edges of a wire through `ShapeFix_Wire::FixReorder()` before transferring them
/// (`cxx:706`). That is the 3D `ShapeAnalysis_Wire::CheckOrder(sawo, myClosedMode,
/// true, false)` of `ShapeFix_Wire::FixReorder()` (`ShapeFix_Wire.cxx:487-534`)
/// followed by `FixReorder(sawo)` (`cxx:1351-1400`). `Perform` ignores its
/// `closed` argument (`ShapeAnalysis_WireOrder.cxx:205`), so only the edge
/// connectivity matters.
///
/// Returns the same list as `wire_data_edges` when the order is unchanged
/// (`status == 0`, `cxx:1360-1363`), when an edge has no vertex (FAIL2,
/// `ShapeAnalysis_Wire.cxx:617-621`) or when the order could not be applied
/// (`cxx:1372-1385`). `Ordered(i) < 0` is `ShapeExtend_WireData::Edge(signed)`
/// = the edge reversed (`ShapeExtend_WireData.cxx:583-591`).
pub(super) fn reordered_wire_edges(wire: &crate::shape::Wire) -> Vec<Edge> {
    let stored = wire_data_edges(wire);
    if stored.len() < 2 {
        return stored;
    }
    let mut order = crate::meshing::wire_order::WireOrder::new();
    for e in &stored {
        let (Some(v1), Some(v2)) = (wire_edge_first_vertex(e), wire_edge_last_vertex(e)) else {
            return stored;
        };
        order.add_edge_xyz(BRepTool::vertex_point(&v1), BRepTool::vertex_point(&v2));
    }
    order.perform();
    if order.status() == crate::meshing::wire_order::WireOrderStatus::Same
        || order.nb_edges() != stored.len()
    {
        return stored;
    }
    let mut out = Vec::with_capacity(stored.len());
    for i in 1..=stored.len() {
        let signed = order.ordered(i);
        if signed == 0 {
            return stored;
        }
        let mut e = stored[signed.unsigned_abs() as usize - 1].clone();
        if signed < 0 {
            e.0.reverse();
        }
        out.push(e);
    }
    out
}

/// Direct children of `s` whose type is `t`.
pub(super) fn children_of_type(s: &TopoShape, t: ShapeType) -> Vec<TopoShape> {
    s.tshape
        .read()
        .unwrap()
        .children
        .iter()
        .filter(|h| h.shape_type() == t)
        .cloned()
        .collect()
}

/// The UV box handed to `GeomToIGES_GeomSurface::TransferSurface`.
///
/// OCCT takes it straight from the face (`BRep_Tool::UVBounds`, i.e.
/// `BRepTools::AddUVBounds` -> the edge p-curve boxes, see
/// [`crate::brep_tools::add_uv_bounds`]). A void box is reported by
/// `BRepTools::UVBounds` as `0,0,0,0` (`BRepTools.cxx:64-80`) and a surface
/// whose natural range is infinite yields no box either; the revolution arms
/// here cannot use either, so both map to `[0,1]^2`.
pub(super) fn face_uv_bounds_finite(f: &Face) -> (f64, f64, f64, f64) {
    let mut b = occt_core::bnd::BndBox2d::new();
    crate::brep_tools::add_uv_bounds(f, &mut b);
    match b.get() {
        Some((u0, v0, u1, v1))
            if u0.is_finite() && u1.is_finite() && v0.is_finite() && v1.is_finite() =>
        {
            (u0, u1, v0, v1)
        }
        _ => (0.0, 1.0, 0.0, 1.0),
    }
}

/// `Geom_ElementarySurface::Bounds` (`Geom_ElementarySurface.cxx:82-89`): every
/// elementary surface reports the same "infinite" range, `U [0, 2*PI],
/// V [0, 2*PI]`. `TransferFace` mirrors the p-curves of a REVERSED face about
/// `U = 0.5*(U1+U2)` of the *reversed* surface's own bounds
/// (`BRepToIGES_BRShell.cxx:176-183`), so this is `2*aCenter` for every surface
/// class that has an `UReversed` branch.
pub(super) const ELEMENTARY_REVERSED_TWO_CENTER: f64 = 2.0 * std::f64::consts::PI;

/// `BRepToIGES_BRShell::TransferFace` (`BRepToIGES_BRShell.cxx:126-237`): a face
/// whose orientation is `TopAbs_REVERSED` is not transferred as it is. OCCT
/// copies the face, replaces its surface with `aSurf->UReversed()` after peeling
/// every `Geom_RectangularTrimmedSurface` (`:135-141`) and mirrors each p-curve
/// about `gp_Ax2d(gp_Pnt2d(0.5*(U1+U2), V1), gp_Dir2d(0,1))` (`:176-201`), then
/// transfers the copy.
///
/// The port transfers the original face (the copied wires only get re-oriented
/// to the originals, `:150-174`) and mirrors the p-curves it reads on the fly
/// (`transfer_edge_uv`), so this carries what both call sites need: the
/// U-reversed surface and the mirror's `2*aCenter`.
///
/// `None` means the REVERSED branch is not taken: either the face is not
/// REVERSED, or its surface class has no `Geom_Surface::UReversed` branch yet
/// (elementary surfaces, a non-periodic-U B-spline, a surface of revolution or
/// linear extrusion, and an offset whose basis can be U-reversed — see
/// `Surface::u_reversed`), in which case the face is transferred as if it were
/// FORWARD, exactly as an unimplemented OCCT branch would leave it.
pub(super) struct FaceReversal {
    /// `aSurf->UReversed()` (`cxx:135-142`), already trim-stripped.
    pub(super) surf: Arc<dyn Surface>,
    /// `2*aCenter` of `gp_Ax2d(gp_Pnt2d(aCenter, V1), gp_Dir2d(0,1))`: the mirror
    /// maps `U` to `2*aCenter - U` (`cxx:178-183`).
    pub(super) two_center: f64,
}

pub(super) fn face_reversal(f: &Face) -> Option<FaceReversal> {
    if f.0.orientation() != crate::abs::Orientation::Reversed {
        return None;
    }
    let surf = BRepTool::face_surface(f)?;
    // `cxx:135-141`: peel the rectangular trims first, "because pcurves will be
    // transformed, so trim will be shifted, accorded to new face bounds".
    let mut base: Arc<dyn Surface> = surf.clone();
    while let Some(b) = base.rectangular_trimmed_basis() {
        base = b;
    }
    let rev = base.u_reversed()?;
    // `cxx:176-182`: `aSurf->Bounds(U1, U2, V1, V2); aCenter = 0.5*(U1+U2)` on the
    // REVERSED surface, so the mirror is `U -> 2*aCenter - U`. Finite bounds
    // (B-spline, revolution, offset of those) use that sum; an unbounded
    // elementary plane keeps the `[0, 2*PI]` constant the other elementary
    // surfaces report (`Geom_ElementarySurface::Bounds`).
    let (ru0, ru1) = rev.u_range();
    let two_center = if ru0.is_finite() && ru1.is_finite() {
        ru0 + ru1
    } else {
        ELEMENTARY_REVERSED_TWO_CENTER
    };
    if crate::debug_flags::iges_trace_face() {
        eprintln!(
            "TRACEREV ptr={:?} surf={:?} two_center={two_center}",
            Arc::as_ptr(&f.0.tshape) as usize,
            crate::pcurve_full::classify_surface_kind(rev.as_ref()),
        );
    }
    Some(FaceReversal {
        surf: rev,
        two_center,
    })
}

/// 'GeomToIGES_GeomSurface::Length()' ('GeomToIGES_GeomSurface.cxx:1441-1443'):
/// the 'TheLength' the matching 'TransferSurface' overload stored for the
/// face's surface. Only two arms are not 1:
///
/// * 'TransferSurface(Geom_ConicalSurface)' ('cxx:777-823') sets it to the
///   generatrix span 'gen1.Distance(gen2)', and the generatrix is a unit
///   direction, so the value is '|V2 - V1|' over the (infinite-substituted) V
///   bounds;
/// * 'TransferSurface(Geom_SurfaceOfRevolution)' ('cxx:1111-1169') sets it the
///   same way only when the trim-unwrapped basis curve is a line.
///
/// Every other overload ('cxx:618', '714', '879', '947', '1047', '1132', '1253')
/// leaves it at 1.
pub(super) fn surface_transfer_length(surf: &dyn Surface, v0: f64, v1: f64) -> f64 {
    let bound = |v: f64, positive: bool| -> f64 {
        if positive {
            if occt_core::precision::Precision::is_positive_infinite(v) {
                occt_core::precision::INFINITE
            } else {
                v
            }
        } else if occt_core::precision::Precision::is_negative_infinite(v) {
            -occt_core::precision::INFINITE
        } else {
            v
        }
    };
    let (vv0, vv1) = (bound(v0, false), bound(v1, true));
    // 'cxx:737-741' (GeomToIGES_GeomSurface::TransferSurface(Geom_Cylinder)):
    // the generatrix is 'Geom_Line(gp_Pnt(R,0,0), gp::DZ())', so
    // 'TheLength = gen1.Distance(gen2) = |V2 - V1|'.
    if surf.gp_cylinder().is_some() {
        return (vv1 - vv0).abs();
    }
    // 'cxx:817-825' (GeomToIGES_GeomSurface::TransferSurface(Geom_ConicalSurface)):
    // the generatrix is 'Geom_Line(gp_Pnt(RefRadius,0,0), gp_Dir(sin a,0,cos a))',
    // a unit direction, so the same expression holds.
    if surf.gp_cone().is_some() {
        return (vv1 - vv0).abs();
    }
    if surf.is_surface_of_revolution() {
        if let Some(c) = surf.revolution_basis_curve() {
            let c: Arc<dyn Curve> = if c.is_geom_trimmed() {
                c.untrimmed_basis().map(|(b, _, _)| b).unwrap_or(c)
            } else {
                c
            };
            if c.is_line() {
                return (vv1 - vv0).abs();
            }
        }
    }
    1.0
}
