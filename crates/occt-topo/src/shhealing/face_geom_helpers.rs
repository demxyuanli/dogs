//! Face-level helpers shared by `ShapeFix_Face`'s second-part loop
//! (`ShapeFix_Face.cxx:500-717`): `IsSurfaceUVInfinite` (`cxx:94-103`),
//! `Shift2dWire` (`cxx:776-806`), `CutInterval` (`cxx:808-850`),
//! `FindBestInterval` (`cxx:853-866`) and the `EmptyCopied` face view.
//!
//! They live here rather than in `shape_fix_face.rs` because that file is
//! already past the repository's split threshold, and because
//! `FixAddNaturalBound`, `FixOrientation` and `FixSplitFace` each need a
//! subset.

use std::sync::Arc;

use occt_core::gp::{GpPnt2d, GpTrsf2d, GpVec2d};
use occt_core::precision::{self, PCONFUSION};
use occt_geom::Surface;

use crate::abs::Orientation;
use crate::boptools_2d::curve_on_surface_oriented;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Face, TopoShape, Wire};
use crate::tgeometry::{GeometryRegistry, VertexGeom};
use crate::topo_tools_full::edges_of_wire;

/// `IsSurfaceUVInfinite(theSurf)` (`ShapeFix_Face.cxx:94-103`):
/// `Precision::IsInfinite` on any of the four bounds.
pub(in crate::shhealing) fn is_surface_uv_infinite(surf: &dyn Surface) -> bool {
    let (umin, umax) = surf.u_range();
    let (vmin, vmax) = surf.v_range();
    precision::Precision::is_infinite(umin)
        || precision::Precision::is_infinite(umax)
        || precision::Precision::is_infinite(vmin)
        || precision::Precision::is_infinite(vmax)
}

/// `TopoDS_Shape::EmptyCopied` for a face: a fresh `TopoDS_Face` carrying the
/// same surface and face tolerance and none of the original children
/// (`TopoDS_Shape::EmptyCopied` -> `TShape::EmptyCopy`).
///
/// The location is copied explicitly because the port's `BRep_Builder::Add`
/// compensates a child's location through the target's
/// (`BRep_Builder.cxx:155-175`), exactly as `EmptyCopied` plus `B.Add` do.
pub(in crate::shhealing) fn empty_copied_face(face: &Face) -> Face {
    let mut new_face = match BRepTool::face_surface(face) {
        Some(surf) => TopoBuilder::new().make_face(surf, &[]),
        None => Face::new(),
    };
    new_face.0.set_orientation(Orientation::Forward);
    new_face.0.set_location(face.0.location());
    let reg = GeometryRegistry::global();
    if let Some(geom) = reg.face_geom(&face.0) {
        reg.set_face(&new_face.0, geom);
    }
    new_face
}

/// `Shift2dWire(w, f, vec, theSurface, recompute3d)`
/// (`ShapeFix_Face.cxx:776-806`): translate every pcurve of `w` on `f` by
/// `vec` and, when `recompute3d` holds, drop the 3D curve, rebuild it from the
/// shifted pcurve and move the edge's first vertex onto the new pcurve point
/// with tolerance `0.` (`B.UpdateVertex(V, P, 0.)`).
pub(in crate::shhealing) fn shift_2d_wire(
    w: &Wire,
    f: &Face,
    vec: &GpVec2d,
    surface: &dyn Surface,
    recompute3d: bool,
) {
    // cxx:781-783.
    let mut tr2d = GpTrsf2d::identity();
    tr2d.set_translation_vec(vec);
    let face_key = GeometryRegistry::shape_key(&f.0);
    let reg = GeometryRegistry::global();
    for edge in edges_of_wire(w) {
        // cxx:788-793: `sae.PCurve(edge, f, C2d, cf, cl, true)`; a failure is
        // `continue`, not an abort. The returned handle is the pcurve of the
        // edge's orientation slot, so writing the transformed curve back into
        // that same slot keeps a seam edge's other pcurve untouched - the C++
        // mutates the one `Geom2d_Curve` handle `PCurve` handed out.
        let reversed = edge.0.orientation() == Orientation::Reversed;
        let pcurves = reg.edge_pcurves(&edge.0, face_key);
        let slot = if pcurves.len() >= 2 && reversed { 1 } else { 0 };
        let Some((c2d, cf, _cl)) = curve_on_surface_oriented(&edge, f, true) else {
            continue;
        };
        if slot >= pcurves.len() {
            continue;
        }
        // cxx:795: `C2d->Transform(tr2d)`; the COS range is unchanged by a
        // rigid 2D translation.
        let mut updated = pcurves.clone();
        updated[slot] = Arc::from(c2d.transformed(&tr2d));
        let range = reg.pcurve_range(&edge.0, face_key);
        reg.set_edge_pcurves(&edge.0, face_key, updated);
        if let Some((first, last)) = range {
            reg.set_pcurve_range(&edge.0, face_key, first, last);
        }
        if !recompute3d {
            continue;
        }
        // cxx:797-802: `sbe.RemoveCurve3d(edge); sbe.BuildCurve3d(edge);`
        // then `B.UpdateVertex(sae.FirstVertex(edge),
        // theSurface->Value(C2d->Value(cf)), 0.)`.
        edge.0
            .tshape
            .write()
            .expect("poisoned TShape lock")
            .edge_core_mut()
            .curve = None;
        crate::shhealing::ShapeBuildEdge.build_curve3d(&edge);
        let Some(v1) = crate::shhealing::first_vertex(&edge) else {
            continue;
        };
        // `theSurface->Value(C2d->Value(cf))` on the transformed curve.
        let Some(p) = reg
            .edge_pcurve(&edge.0, face_key)
            .map(|pc| surface.d0(pc.d0(cf).x(), pc.d0(cf).y()))
        else {
            continue;
        };
        reg.set_vertex(&v1.0, VertexGeom { point: p, tolerance: 0.0 });
    }
}

/// `CutInterval(intervals, toAddI, period)` (`ShapeFix_Face.cxx:808-850`):
/// cuts `toAddI` out of `intervals`, splitting the interval it lands in.
/// Returns `false` only for an empty sequence (`cxx:810-813`).
pub(in crate::shhealing) fn cut_interval(
    intervals: &mut Vec<GpPnt2d>,
    to_add_i: &GpPnt2d,
    period: f64,
) -> bool {
    // cxx:810-813.
    if intervals.is_empty() {
        return false;
    }
    // cxx:815-845: try twice, aligning to the bottom then to the top.
    for j in 0..2 {
        let mut i = 0usize;
        while i < intervals.len() {
            let interval = intervals[i];
            // cxx:818-826.
            let shift = crate::shhealing::adjust_by_period(
                if j == 1 { to_add_i.x() } else { to_add_i.y() },
                0.5 * (interval.x() + interval.y()),
                period,
            );
            let to_add = GpPnt2d::new(to_add_i.x() + shift, to_add_i.y() + shift);
            // cxx:827-829.
            if to_add.y() <= interval.x() || to_add.x() >= interval.y() {
                i += 1;
                continue;
            }
            // cxx:830-844.
            if to_add.x() > interval.x() {
                if to_add.y() < interval.y() {
                    // cxx:833-836: `InsertBefore(i, interval)` then
                    // `ChangeValue(i + 1).SetX(toAdd.Y())`.
                    intervals.insert(i, interval);
                    intervals[i + 1] = GpPnt2d::new(to_add.y(), interval.y());
                }
                // cxx:837: `ChangeValue(i).SetY(toAdd.X())`.
                intervals[i] = GpPnt2d::new(interval.x(), to_add.x());
                // cxx:836's trailing `i++`: move past the split.
                i += 1;
            } else if to_add.y() < interval.y() {
                // cxx:839-842.
                intervals[i] = GpPnt2d::new(to_add.y(), interval.y());
                i += 1;
            } else {
                // cxx:843-846: `Remove(i--)`, so the same index is re-examined.
                intervals.remove(i);
            }
        }
    }
    true
}

/// `FindBestInterval(intervals)` (`ShapeFix_Face.cxx:853-866`): the middle of
/// the widest interval, `0.` when none is wider than `-1`.
pub(in crate::shhealing) fn find_best_interval(intervals: &[GpPnt2d]) -> f64 {
    let mut shift = 0.0;
    let mut max = -1.0;
    for interval in intervals.iter() {
        // cxx:859-862.
        if interval.y() - interval.x() <= max {
            continue;
        }
        max = interval.y() - interval.x();
        shift = interval.x() + 0.5 * max;
    }
    shift
}

/// `shift.XY().Modulus() < Precision::PConfusion()` (`ShapeFix_Face.cxx:1020`).
pub(in crate::shhealing) fn shift_is_negligible(shift: &GpPnt2d) -> bool {
    GpVec2d::new(shift.x(), shift.y()).magnitude() < PCONFUSION
}

/// `TopoDS_Iterator(face, false)` split by `ShapeFix_Face.cxx:892-904`:
/// `(oriented wires, everything else)`.
pub(in crate::shhealing) fn split_face_children(face: &Face) -> (Vec<Wire>, Vec<TopoShape>) {
    let children: Vec<TopoShape> = {
        let ts = face.0.tshape.read().expect("poisoned TShape lock");
        ts.children.clone()
    };
    let mut wires = Vec::new();
    let mut others = Vec::new();
    for child in children {
        let ori = child.orientation();
        if child.shape_type() == crate::abs::ShapeType::Wire
            && (ori == Orientation::Forward || ori == Orientation::Reversed)
        {
            wires.push(Wire(child));
        } else {
            others.push(child);
        }
    }
    (wires, others)
}
