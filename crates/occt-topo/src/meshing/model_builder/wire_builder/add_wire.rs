use super::pcurve::pcurve_and_range;

use super::super::prelude::*;

/// Port of `BRepMesh_ShapeVisitor::addWire` (2D pcurve mode).
///
/// Orders the wire's edges into one connected chain
/// (`ShapeAnalysis_Wire::CheckOrder` → `ShapeAnalysis_WireOrder::Perform`) and,
/// for every non-EXTERNAL edge, registers its pcurve on the face and its
/// position/orientation in the wire chain. Returns `false` when the wire cannot
/// be ordered (a missing pcurve → `ShapeExtend_FAIL`); a wire whose edges needed
/// reversing sets the face `UNORIENTED_WIRE` status (`ShapeExtend_DONE3`).
pub(super) fn add_wire(
    model: &mut MeshModel,
    face_index: usize,
    wire: &Wire,
    edge_index: &mut HashMap<usize, usize>,
) -> bool {
    let face = match model.face(face_index) {
        Ok(f) => f.face().clone(),
        Err(_) => return false,
    };
    // `ShapeExtend_WireData(theWire, chained=true, manifold=false)` keeps only
    // FORWARD/REVERSED edges in the main list; INTERNAL/EXTERNAL are non-manifold.
    // `ShapeExtend_WireData::Init`: FORWARD/REVERSED only; a REVERSED wire
    // prepends each edge so the list is the reverse of iterator order
    // (`ShapeExtend_WireData.cxx:115-122`).
    //
    // OCCT's mesh model builder applies NO orientation fixup here:
    // `BRepMesh_ShapeVisitor::addWire` (`BRepMesh_ShapeVisitor.cxx:93-134`) only
    // runs `ShapeExtend_WireData(wire, true, false)` + `ShapeAnalysis_Wire::
    // CheckOrder` and stores the resulting chain as-is. `ShapeFix_Face::
    // FixOrientation` is a `ShapeFix` pass and is not part of the mesher, so the
    // wire's UV winding is left untouched (the clockwise/ccw split is the
    // classifier's job, not this loop's).
    //
    // Reversing a single clockwise wire here flipped face 433/435 of
    // `data/occ/acs10.stp` from `ShapeAnalysis_WireOrder`'s reference chain
    // (`e0,e15,e14,...,e1`) to its mirror (`e15,e0,...,e14`) and starved the
    // Delaunay cavity: 49/50 triangles instead of OCCT's 73/74.
    let work = wire.clone();
    let face_fwd = Face(face.0.oriented(Orientation::Forward));
    let mut stored: Vec<Edge> = edges_of_wire(&work)
        .into_iter()
        .filter(|e| {
            let o = e.0.orientation();
            o == Orientation::Forward || o == Orientation::Reversed
        })
        .collect();
    if work.0.orientation().is_reversed() {
        stored.reverse();
    }
    if stored.is_empty() {
        return false;
    }

    // `ShapeAnalysis_Wire::CheckOrder(..., isClosed=true, mode3d=false)`:
    // pcurve endpoints on `face.Oriented(FORWARD)` with `PCurve(..., orient=true)`
    // (reversed edges swap cf/cl).
    let mut order = WireOrder::new();
    for e in &stored {
        // Port of `ShapeAnalysis_Wire::CheckOrder` (`ShapeAnalysis_Wire.cxx:627-649`):
        //   ShapeAnalysis_Edge EA; EA.PCurve(E, face.Oriented(FORWARD), c2d, f, l);
        //   gp_Pnt2d (c2d->Value(f)), (c2d->Value(l));
        // `PCurve(..., orient=true)` (`ShapeAnalysis_Edge.cxx:192-208`) reads
        // `BRep_Tool::CurveOnSurface(E, F, f, l)` - the **COS representation's**
        // range - and swaps f/l on a REVERSED edge. Seam edges have two pcurves
        // (one per side) and each traversal must get its own side, so the pcurve
        // is selected by the edge's *actual* orientation.
        let (pc, a, b) = match pcurve_and_range(e, &face_fwd) {
            Some(v) => v,
            None => return false,
        };
        if !a.is_finite() || !b.is_finite() || (b - a).abs() < 1e-15 {
            return false;
        }
        // `ShapeAnalysis_Wire::CheckOrder` (`cxx:648-649`): `c2d->Value(f/l)`
        // after `PCurve(..., orient=true)` toggles f/l on a REVERSED edge.
        let (p_a, p_b) = (pc.d0(a), pc.d0(b));
        let (begin, end) = if e.0.orientation().is_reversed() {
            (p_b, p_a)
        } else {
            (p_a, p_b)
        };
        order.add_edge(begin, end);
    }
    order.perform();
    // `BRepMesh_ShapeVisitor::addWire` (`BRepMesh_ShapeVisitor.cxx:99-134`)
    // applies no chain-level reversal after `CheckOrder`: it only reads
    // `aOrderTool.Ordered(i)` and lets `ShapeExtend_WireData::Edge(signed)`
    // reverse each edge individually (`ShapeExtend_WireData.cxx:583-588`). The
    // wire's UV winding stays whatever the oriented edges give (hole wires stay
    // clockwise); the inside/outside split is done later by the classifier, not
    // by flipping the whole chain on the sign of a signed area.

    if order.status() == WireOrderStatus::Reversed {
        if let Ok(f) = model.face_mut(face_index) {
            f.set_status(MeshStatus::UNORIENTED_WIRE);
        }
    }
    if order.nb_edges() != stored.len() {
        return false;
    }

    // `ShapeExtend_WireData::Edge(signed)` (`cxx:583-588`): `Ordered < 0`
    // is `Edge.Reverse()`, then `AddPCurve` / `AddEdge` use that orientation
    // (`BRepMesh_ShapeVisitor.cxx:123-130`). On a seam `Reverse` selects
    // PCurve2 (`BRep_Tool.cxx:354-357`).
    let wire_index = model.add_wire(wire.clone());
    // TEMPORARY diagnostic (DELREGSKEL): the `ShapeExtend_WireData` list fed to
    // `ShapeAnalysis_WireOrder::Perform` and the `Ordered(i)` chain it returns, in
    // the `occt_probe --orderdbg` `ODBG` line format.
    if crate::debug_flags::del_reg_skel() {
        {
            let raw = work.0.tshape.read().expect("lock").children.clone();
            let rawlist = raw
                .iter()
                .map(|c| format!("{:?}:{:?}", c.shape_type(), c.orientation()))
                .collect::<Vec<_>>()
                .join(" ");
            eprintln!(
                "WIRE face={} wire={} wptr={:#x} work_ori={:?} raw_children=[{}]",
                face_index,
                wire_index,
                Arc::as_ptr(&work.0.tshape) as usize,
                work.0.orientation(),
                rawlist
            );
        }
        let list = stored
            .iter()
            .map(|e| format!("{:p}:{:?}", Arc::as_ptr(&e.0.tshape), e.0.orientation()))
            .collect::<Vec<_>>()
            .join(" ");
        let chain = (1..=stored.len())
            .map(|i| order.ordered(i).to_string())
            .collect::<Vec<_>>()
            .join(" ");
        eprintln!(
            "WIRE face={} wire={} ori={:?} stored=[{}] ordered=[{}]",
            face_index,
            wire_index,
            work.0.orientation(),
            list,
            chain
        );
        for (i, e) in stored.iter().enumerate() {
            let n_pc = GeometryRegistry::global()
                .edge_pcurves(&e.0, GeometryRegistry::shape_key(&face_fwd.0))
                .len();
            match pcurve_and_range(e, &face_fwd) {
                Some((pc, a, b)) => {
                    let (pa, pb) = (pc.d0(a), pc.d0(b));
                    let (bu, eu) = if e.0.orientation().is_reversed() {
                        (pb, pa)
                    } else {
                        (pa, pb)
                    };
                    eprintln!(
                        "WIRE   i={} ptr={:p} ori={:?} npc={} pc2d=({:.9},{:.9})->({:.9},{:.9}) trav=({:.9},{:.9})->({:.9},{:.9})",
                        i + 1,
                        Arc::as_ptr(&e.0.tshape),
                        e.0.orientation(),
                        n_pc,
                        pa.x(),
                        pa.y(),
                        pb.x(),
                        pb.y(),
                        bu.x(),
                        bu.y(),
                        eu.x(),
                        eu.y()
                    );
                }
                None => {
                    eprintln!(
                        "WIRE   i={} ptr={:p} ori={:?} npc={} no-pcurve",
                        i + 1,
                        Arc::as_ptr(&e.0.tshape),
                        e.0.orientation(),
                        n_pc
                    );
                }
            }
        }
    }
    for i in 1..=stored.len() {
        let signed = order.ordered(i);
        let e = &stored[signed.unsigned_abs() as usize - 1];
        let mut orientation = e.0.orientation();
        if orientation == Orientation::External {
            continue;
        }        // `ShapeExtend_WireData::Edge(signed)` is `Edge.Reverse()`. On a plane
        // that is the same pcurve walked backwards. On a seam `Reverse` selects
        // PCurve2 (`BRep_Tool.cxx:354-357`), leaving the CheckOrder iso; walk
        // that measured chord backwards instead.
        let mut reverse_walk = false;
        if signed < 0 {
            let face_key = GeometryRegistry::shape_key(&face_fwd.0);
            let n_pc = GeometryRegistry::global().edge_pcurves(&e.0, face_key).len();
            if n_pc >= 2 {
                reverse_walk = true;
            } else {
                orientation = orientation.reversed();
            }
        }
        let key = Arc::as_ptr(&e.0.tshape) as usize;
        let eidx = match edge_index.get(&key) {
            Some(&i) => i,
            None => {
                let i = model.add_edge(e.clone());
                edge_index.insert(key, i);
                i
            }
        };
        if let Ok(edge) = model.edge_mut(eidx) {
            // Keep the wire occurrence (location) on the pcurve: a prism side
            // reuses the profile `TShape` for bottom and top under different
            // locations (`BRep_Tool::CurveOnPlane` relative location).
            let occurrence = Edge(e.0.oriented(orientation));
            edge.add_pcurve(face_index, orientation, occurrence);
        }
        if let Ok(w) = model.wire_mut(wire_index) {
            w.add_edge_ordered(eidx, orientation, reverse_walk);
        }
    }
    if let Ok(f) = model.face_mut(face_index) {
        f.add_wire(wire_index);
    }
    true
}
