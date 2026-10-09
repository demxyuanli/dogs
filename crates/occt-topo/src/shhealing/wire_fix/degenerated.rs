use super::*;

/// `ShapeFix_Wire::FixDegenerated` (`ShapeFix_Wire.cxx:2130-2205`), the pass
/// `ShapeFix_Wire::Perform` runs at `cxx:386-392` under `myFixDegeneratedMode`
/// (`FromSTEP.FixShape.FixDegeneratedMode`, `STEPControl_Controller.cxx:236`,
/// default `-1`). The read-time chain is `ShapeProcess_OperLibrary.cxx:846`
/// (`sfw->FixDegeneratedMode() = ctx->IntegerVal("FixDegeneratedMode", -1)`) ->
/// `ShapeFix_Face::Perform` -> its first wire loop (`ShapeFix_Face.cxx:423`
/// `theAdvFixWire->Perform`, `StatusDegenerated` read at `cxx:431`; the second
/// loop disables the mode again at `cxx:520`) -> `ShapeFix_Wire::Perform`.
/// `NeedFix` maps `-1` to true (`ShapeFix_Root.lxx:101-104` with the default
/// `def = true`, `ShapeFix_Root.hxx:112`), so OCCT does run this pass on a STEP
/// read.
///
/// `Shape-1`'s 16 `SPHERICAL_SURFACE` faces are the read-shape effect of this
/// pass: the file has 116 `EDGE_CURVE`s while the OCCT read shape has 132
/// edges, the 16 extra being one degenerated edge per spherical surface. Each
/// adds one vertex reference forward and the same one reversed (`+301 0 -301 0
/// *` in the OCCT read shape), so the read shape keeps the file's 61
/// `VERTEX_POINT`s with no new vertex -- `FixMissingSeam` would have created one
/// per edge (`ShapeFix_Face.cxx:1984` `B.MakeVertex(V, mySurf->Value(...))`) and
/// given 77. All 16 new edges carry the range `0 1.01893653209489`, which is
/// `vect2d.Magnitude()` for a sphere at `-22.5 -6.129 2.5` r=1.9 (the face's own
/// u-span; all 16 spheres are the same shape, so it is one value), not `2*pi`,
/// again ruling out `FixMissingSeam`'s sphere branch (`ShapeFix_Face.cxx:1920-1925`,
/// `aRange = 2*pi`; that pass is enabled at read time too,
/// `ShapeProcess_OperLibrary.cxx:830`, so it is the range, not the mode, that
/// excludes it). Wire and vertex counts are unchanged (63 wires, 61 vertices in
/// file and reference).
///
/// `ShapeAnalysis_Wire::CheckDegenerated` (`ShapeAnalysis_Wire.cxx:896-1116`).
///
/// `num` is the 1-based wire position of the edge to test (`0` = the last edge,
/// `cxx:903`). The singularity side reuses `SurfaceSingularities`
/// (`pcurve_full/pcurve_ranges.rs`), the port of `ShapeAnalysis_Surface::ComputeSingularities`
/// / `DegeneratedValues` / `IsDegenerated`.
pub(in crate::shhealing) enum DegeneratedCheck {
    /// Returned false without `ShapeExtend_FAIL2`: nothing to fix.
    None,
    /// `ShapeExtend_FAIL2` (`cxx:930-931`, `cxx:945-947`, `cxx:1047-1048`,
    /// `cxx:1093-1097`): `ShapeFix_Wire::FixDegenerated` drops the edge
    /// (`cxx:2148-2152`).
    Remove,
    /// `ShapeExtend_DONE` (`cxx:1110`): insert (`lack` = `ShapeExtend_DONE1`,
    /// `ShapeExtend_WireData::Add`) or replace (`ShapeExtend_DONE2`,
    /// `ShapeExtend_WireData::Set`) a degenerated edge spanning `p2d1`..`p2d2`.
    Found { p2d1: GpPnt2d, p2d2: GpPnt2d, lack: bool },
}

/// `ShapeAnalysis_Edge::HasPCurve(edge, face)` (`ShapeAnalysis_Edge.cxx:183-191`):
/// a `BRep_Tool::CurveOnSurface` representation exists.
pub(in crate::shhealing) fn has_pcurve(edge: &Edge, face: &Face) -> bool {
    curve_on_surface_oriented(edge, face, true).is_some()
}

pub(in crate::shhealing) fn check_degenerated(wire: &Wire, face: &Face, prec: f64, num: usize) -> DegeneratedCheck {
    let Some(surf) = BRepTool::face_surface(face) else {
        return DegeneratedCheck::None;
    };
    let edges = edges_of_wire(wire);
    let nb = edges.len();
    if nb < 1 || num > nb {
        return DegeneratedCheck::None;
    }
    // `cxx:903-905`.
    let n2 = if num > 0 { num } else { nb };
    let n1 = if n2 > 1 { n2 - 1 } else { nb };
    let n3 = if n2 < nb { n2 + 1 } else { 1 };
    let e1 = edges[n1 - 1].clone();
    let e2 = edges[n2 - 1].clone();
    let e3 = edges[n3 - 1].clone();

    // `cxx:913-934`: an edge already flagged degenerated *and* carrying a
    // pcurve is only re-checked for a pcurve that no longer spans the
    // singularity (OCC7630) -- that one is removable.
    if BRepTool::is_degenerated(&e2) && has_pcurve(&e2, face) {
        if has_pcurve(&e1, face) && has_pcurve(&e3, face) {
            let (Some(pc2), Some(pc1), Some(pc3)) = (
                pcurve_at(&e2, face),
                pcurve_at(&e1, face),
                pcurve_at(&e3, face),
            ) else {
                return DegeneratedCheck::None;
            };
            let p21 = pc2.0.d0(pc2.1);
            let p22 = pc2.0.d0(pc2.2);
            let p12 = pc1.4;
            let p31 = pc3.3;
            if (p12.distance(&p31) - p21.distance(&p22)).abs() > 2.0 * PCONFUSION {
                return DegeneratedCheck::Remove;
            }
        }
        return DegeneratedCheck::None;
    }

    // `cxx:938-948`: `n1 != n2` (OCC320, two sequences of degenerated edges on
    // separate surface bounds) and the previous edge is degenerated without a
    // pcurve.
    if n1 != n2 && BRepTool::is_degenerated(&e1) && !has_pcurve(&e1, face) {
        if BRepTool::is_degenerated(&e2) {
            return DegeneratedCheck::Remove;
        }
        return DegeneratedCheck::None;
    }

    // `cxx:950-970`.
    let (Some(vp), Some(v0), Some(v1v), Some(v2v)) = (
        first_vertex(&e1),
        last_vertex(&e1),
        first_vertex(&e2),
        last_vertex(&e2),
    ) else {
        return DegeneratedCheck::None;
    };
    let pp = BRepTool::vertex_point(&vp);
    let p0 = BRepTool::vertex_point(&v0);
    let p1 = BRepTool::vertex_point(&v1v);
    let p2 = BRepTool::vertex_point(&v2v);
    let tol1 = BRepTool::vertex_tolerance(&v1v);
    let prec_first = prec.min(tol1);
    let prec_fin = prec.max(tol1);
    let prec_vtx = if prec < tol1 { 2.0 * prec_fin } else { prec_fin };
    // OCCT's `forward` (`cxx:959`) is unused by
    // `ShapeAnalysis_Surface::DegeneratedValues` (`cxx:373-411`,
    // `const bool /*forward*/`).

    let sing = crate::pcurve_full::SurfaceSingularities::compute(surf.as_ref());
    let mut p2d1 = GpPnt2d::new(0.0, 0.0);
    let mut p2d2 = GpPnt2d::new(0.0, 0.0);
    let mut dgnr = false;
    // `cxx:974-991`: the edge's own wire is already closed on the singular
    // point.
    if p1.distance(&p2) <= prec_first {
        if let Some((a, b, _, _)) = sing.min_gap(&p1, prec_vtx) {
            p2d1 = a;
            p2d2 = b;
            dgnr = true;
            // `cxx:979-990`: do not turn a closed edge whose mid-point is away
            // from the singular point into a degenerated one.
            let (a3, b3) = BRepTool::edge_parameters(&e2);
            if a3.is_finite() && b3.is_finite() {
                if let Some(c3d) = BRepTool::edge_curve(&e2) {
                    let pm = c3d.d0(0.5 * (a3 + b3));
                    if pm.square_distance(&p1) > prec_vtx * prec_vtx {
                        dgnr = false;
                    }
                }
            }
        }
    }
    let mut lack = false;
    if !dgnr {
        // `cxx:995-999`.
        if n1 != n2
            && p1.distance(&pp) <= prec_first
            && sing.is_degenerated(&pp, prec_first)
            && !BRepTool::is_degenerated(&e1)
        {
            return DegeneratedCheck::None;
        }
        // `cxx:1007-1037`: a missing degenerated edge, found by the same
        // minimum-gap singularity search as `DegeneratedValues`.
        if p0.distance(&p1) <= prec_fin {
            if let Some((a, b, _, _)) = sing.min_gap(&p1, prec_vtx) {
                p2d1 = a;
                p2d2 = b;
                lack = true;
            }
        }
    }
    if !lack && !dgnr {
        if BRepTool::is_degenerated(&e2) && !has_pcurve(&e2, face) {
            return DegeneratedCheck::Remove;
        }
        return DegeneratedCheck::None;
    }

    // `cxx:1064-1086`: parametrize the new pcurve exactly from the end of the
    // previous edge's pcurve to the start of the next one.
    if lack || n1 != n2 {
        if let Some((_, _, _, _, pb)) = pcurve_at(&e1, face) {
            p2d1 = pb;
        }
        let target = if dgnr { &e3 } else { &e2 };
        if let Some((_, _, _, a, _)) = pcurve_at(target, face) {
            p2d2 = a;
        }
    }

    // `cxx:1089-1098`: the fix is postponed to `ShapeFix_Wire::FixLacking` when
    // the pcurve is no longer degenerate.
    if !is_degenerated_2d(surf.as_ref(), p2d1, p2d2, prec_vtx, 10.0) {
        if BRepTool::is_degenerated(&e2) {
            return DegeneratedCheck::Remove;
        }
        return DegeneratedCheck::None;
    }

    // `cxx:1101-1107`: the parametric space is already closed. `gp::Resolution()`
    // is `REAL_SMALL`.
    let ads_max = occt_geom::approx_same_parameter::u_resolution(surf.as_ref(), prec)
        .max(occt_geom::approx_same_parameter::v_resolution(surf.as_ref(), prec));
    if p2d1.distance(&p2d2) <= ads_max + occt_core::precision::REAL_SMALL {
        return DegeneratedCheck::None;
    }
    DegeneratedCheck::Found { p2d1, p2d2, lack }
}

/// `ShapeExtend_WireData::Add(edge, atnum)` (`ShapeExtend_WireData.cxx:259-281`):
/// `myEdges->InsertBefore(atnum, edge)` with a 1-based `atnum` (`0` appends).
/// `TShape::children` is the port's `myEdges`, so the new edge is inserted at
/// the same position instead of through `TopoBuilder::add_edge` (which is
/// guarded by `TShape::free()` and would skip a wire already bound to its
/// face, `builder.rs:32-52`).
pub(in crate::shhealing) fn wire_insert_edge_before(wire: &mut Wire, at: usize, edge: &Edge) {
    let mut t = wire.0.tshape.write().expect("poisoned TShape lock");
    let pos = if at == 0 {
        t.children.len()
    } else {
        (at - 1).min(t.children.len())
    };
    t.children.insert(pos, edge.0.clone());
    t.set_modified(true);
}

/// `ShapeExtend_WireData::Set(edge, num)` (`ShapeExtend_WireData.cxx:456-476`):
/// `myEdges->SetValue(num, edge)`, which ignores an out-of-range `num`.
pub(in crate::shhealing) fn wire_set_edge(wire: &mut Wire, at: usize, edge: &Edge) {
    let mut t = wire.0.tshape.write().expect("poisoned TShape lock");
    let len = t.children.len();
    let pos = if at == 0 { len.saturating_sub(1) } else { at - 1 };
    if pos < len {
        t.children[pos] = edge.0.clone();
        t.set_modified(true);
    }
}

/// `ShapeExtend_WireData::Set(edge, num)` for an edge read back through
/// [`edges_of_wire`]. `ShapeExtend_WireData::Init` stores the edges with the
/// wire's own orientation already composed (`ShapeExtend_WireData.cxx:114-121`)
/// and `ShapeExtend_WireData::Wire()` rebuilds a fresh FORWARD wire out of that
/// list (`ShapeExtend_WireData.cxx:651-685`). Our children stay under the
/// original wire, so undo the wire's own orientation before storing, exactly as
/// `fix_reorder_wire` does; `edges_of_wire` then re-composes to the value OCCT
/// would have.
pub(in crate::shhealing) fn wire_set_edge_composed(wire: &mut Wire, at: usize, edge: &Edge) {
    let mut e = edge.clone();
    if wire.0.orientation() == Orientation::Reversed {
        e.0.reverse();
    }
    wire_set_edge(wire, at, &e);
}

/// `ShapeExtend_WireData::Remove(num)` (`ShapeExtend_WireData.cxx:446-452`).
pub(in crate::shhealing) fn wire_remove_edge(wire: &mut Wire, at: usize) {
    let mut t = wire.0.tshape.write().expect("poisoned TShape lock");
    let len = t.children.len();
    let pos = if at == 0 { len.saturating_sub(1) } else { at - 1 };
    if pos < len {
        t.children.remove(pos);
        t.set_modified(true);
    }
}

/// Number of edges `ShapeExtend_WireData` holds (`ShapeExtend_WireData.cxx:576-579`).
pub(in crate::shhealing) fn wire_edges_nb(wire: &Wire) -> usize {
    edges_of_wire(wire).len()
}

/// `ShapeFix_Wire::FixDegenerated(const int num)`
/// (`ShapeFix_Wire.cxx:2130-2205`).
///
/// Returns `(ShapeExtend_DONE, ShapeExtend_DONE2)`; the `DONE2` bit is what the
/// no-argument driver uses to drop duplicated degenerated edges (`cxx:1050`).
pub(in crate::shhealing) fn fix_degenerated(wire: &mut Wire, face: &Face, prec: f64, num: usize) -> (bool, bool) {
    let nb = wire_edges_nb(wire);
    if nb < 1 {
        return (false, false);
    }
    match check_degenerated(wire, face, prec, num) {
        DegeneratedCheck::None => (false, false),
        // `cxx:2148-2152`: FAIL2 -> `WireData()->Remove(num)` + DONE3.
        DegeneratedCheck::Remove => {
            wire_remove_edge(wire, num);
            (true, false)
        }
        DegeneratedCheck::Found { p2d1, p2d2, lack } => {
            let Some(surf) = BRepTool::face_surface(face) else {
                return (false, false);
            };
            // `cxx:2160-2163`.
            let vect2d = GpVec2d::new(p2d2.x() - p2d1.x(), p2d2.y() - p2d1.y());
            let Ok(dir2d) = GpDir2d::from_vec2d(&vect2d) else {
                return (false, false);
            };
            let line2d: Arc<dyn Curve2d> = Arc::new(Geom2dLine::from_pnt_dir(p2d1, dir2d));
            let mag = vect2d.magnitude();

            // `cxx:2165-2168`: `B.MakeEdge(degEdge)`, `B.Degenerated(degEdge,
            // true)`, `B.UpdateEdge(degEdge, line2d, Face(), Precision::Confusion())`,
            // `B.Range(degEdge, Face(), 0., vect2d.Magnitude())`.
            //
            // `BRep_Builder::Degenerated` clears the 3D curve (`cxx:1082-1084`),
            // and this port's `EdgeGeom` cannot hold a null one, so the edge is
            // built over the curve `BRepAdaptor_Curve(edge, face)` produces once
            // the 3D curve is gone: the pcurve image on the surface
            // (`Adaptor3d_CurveOnSurface`, `meshing/edge_discret.rs`).
            let deg_curve: Arc<dyn Curve> = Arc::new(
                crate::meshing::edge_discret::CurveOnSurface::new(
                    line2d.clone(),
                    surf.clone(),
                    0.0,
                    mag,
                ),
            );
            let b = TopoBuilder::new();
            let mut deg = b.make_edge(deg_curve, 0.0, mag);
            let reg = GeometryRegistry::global();
            let face_key = GeometryRegistry::shape_key(&face.0);
            reg.set_degenerated(&deg.0, true);
            reg.set_edge_pcurve(&deg.0, face_key, line2d);
            reg.set_pcurve_range(&deg.0, face_key, 0.0, mag);

            // `cxx:2171-2186`.
            let n2 = if num > 0 { num } else { nb };
            let n1 = if n2 > 1 { n2 - 1 } else { nb };
            let n3 = if lack { n2 } else if n2 < nb { n2 + 1 } else { 1 };
            let edges = edges_of_wire(wire);
            if let Some(v1) = last_vertex(&edges[n1 - 1]) {
                let mut v = v1.0.clone();
                v.set_orientation(Orientation::Forward);
                b.add(&mut deg.0, &v);
            }
            if let Some(v2) = first_vertex(&edges[n3 - 1]) {
                let mut v = v2.0.clone();
                v.set_orientation(Orientation::Reversed);
                b.add(&mut deg.0, &v);
            }
            deg.0.set_orientation(Orientation::Forward);

            // `cxx:2188-2198`: `Add` (DONE1) inserts, `Set` (DONE2) replaces.
            if lack {
                wire_insert_edge_before(wire, n2, &deg);
            } else {
                wire_set_edge(wire, n2, &deg);
            }
            (true, !lack)
        }
    }
}

/// `ShapeFix_Wire::FixDegenerated()` (`ShapeFix_Wire.cxx:1034-1076`), the
/// driver `ShapeFix_Wire::Perform` calls at `cxx:386-392`.
pub fn fix_degenerated_all(wire: &mut Wire, face: &Face, prec: f64) -> bool {
    let mut done = false;
    let mut last_coded = -1i32;
    let mut prev_coded = 0i32;
    // `myClosedMode` is true (`ShapeFix_Wire.cxx:171`) so `stop = 0`
    // (`cxx:1045`) and the scan runs down to index 1.
    let mut i = wire_edges_nb(wire) as isize;
    while i > 0 {
        let idx = i as usize;
        let (d, coded2) = fix_degenerated(wire, face, prec, idx);
        done |= d;
        let coded = i32::from(coded2);
        if last_coded == -1 {
            last_coded = coded;
        }
        // `cxx:1050-1071`: PRO7226 -- drop a duplicated degenerated edge and
        // clear the flag of the edge that follows it.
        if coded == 1
            && (prev_coded == 1 || (i == 1 && last_coded == 1))
            && wire_edges_nb(wire) > 1
        {
            wire_remove_edge(wire, idx);
            if prev_coded == 0 {
                i = wire_edges_nb(wire) as isize;
            }
            if let Some(e) = edges_of_wire(wire).get(i as usize - 1) {
                GeometryRegistry::global().set_degenerated(&e.0, false);
            }
            // `B.Degenerated(sbwd->Edge(i++), false)` -- the `i++` is undone by
            // the loop's `i--`.
            i += 1;
            prev_coded = 0;
        } else {
            prev_coded = coded;
        }
        i -= 1;
    }
    done
}
