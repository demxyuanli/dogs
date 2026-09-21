use super::prelude::*;
use super::*;

/// Port of `BRepMesh_ShapeVisitor::addWire` (2D pcurve mode).
///
/// Orders the wire's edges into one connected chain
/// (`ShapeAnalysis_Wire::CheckOrder` → `ShapeAnalysis_WireOrder::Perform`) and,
/// for every non-EXTERNAL edge, registers its pcurve on the face and its
/// position/orientation in the wire chain. Returns `false` when the wire cannot
/// be ordered (a missing pcurve → `ShapeExtend_FAIL`); a wire whose edges needed
/// reversing sets the face `UNORIENTED_WIRE` status (`ShapeExtend_DONE3`).
fn add_wire(
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
    // `ShapeFix_Face::FixOrientation` single-wire branch (`cxx:1254-1271`):
    // `!IsOuterBound` → `ShapeExtend_WireData::Reverse(face)` (edge Reverse +
    // list Reverse + SwapSeam on FORWARD seams). STEP `FACE_BOUND(.F.)` leaves
    // the torus loop clockwise. Multi-wire faces keep hole winding.
    let mut work = wire.clone();
    let face_fwd = Face(face.0.oriented(Orientation::Forward));
    if wires_of_face(&face_fwd).len() == 1 && wire_area_2d(&work, &face_fwd) < 0.0 {
        reverse_wire_on_face(&mut work, &face_fwd);
    }
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
        // Pass the edge's *actual* orientation: a seam edge has two pcurves
        // (one per side) and must hand each traversal its own side, not the
        // forward pcurve twice.
        let pc = match make_pcurve_full(e, &face_fwd) {
            Ok(pc) => pc,
            Err(_) => return false,
        };
        let fwd = Edge(e.0.oriented(Orientation::Forward));
        let (a, b) = BRepTool::edge_parameters(&fwd);
        if !a.is_finite() || !b.is_finite() || b - a < 1e-15 {
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
    // CheckOrder follows the first stored chord. After `Reverse(face)` that
    // chord can still start a clockwise UV cycle. On a single-wire face,
    // flip the chain when its area is negative so the discrete walk is an
    // outer bound (`ShapeAnalysis::IsOuterBound`, `TotCross2D >= 0`).
    // Multi-wire holes stay clockwise.
    if wires_of_face(&face_fwd).len() == 1 && order.chain_area() < 0.0 {
        order.reverse_chain();
    }

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
    for i in 1..=stored.len() {
        let signed = order.ordered(i);
        let e = &stored[signed.unsigned_abs() as usize - 1];
        let mut orientation = e.0.orientation();
        if orientation == Orientation::External {
            continue;
        }
        // `ShapeExtend_WireData::Edge(signed)` is `Edge.Reverse()`. On a plane
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
            edge.add_pcurve(face_index, orientation);
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

/// `ShapeAnalysis::OuterWire` over an already-oriented wire list.
pub(crate) fn outer_of_wires(wires: &[Wire], face: &Face) -> Option<Wire> {
    if wires.is_empty() {
        return None;
    }
    for (i, w) in wires.iter().enumerate() {
        if i == wires.len() - 1 {
            return Some(w.clone());
        }
        if wire_area_2d(w, face) >= 0.0 {
            return Some(w.clone());
        }
    }
    None
}

/// `ShapeExtend_WireData(wire, chained=true)` edge list: FORWARD/REVERSED
/// only; a REVERSED wire prepends so iterator order is reversed
/// (`ShapeExtend_WireData.cxx:115-122`).
fn wire_edges_sewd(wire: &Wire) -> Vec<Edge> {
    let mut stored: Vec<Edge> = edges_of_wire(wire)
        .into_iter()
        .filter(|e| {
            let o = e.0.orientation();
            o == Orientation::Forward || o == Orientation::Reversed
        })
        .collect();
    if wire.0.orientation().is_reversed() {
        stored.reverse();
    }
    stored
}

fn edge_oriented_ends(e: &Edge) -> (Option<Vertex>, Option<Vertex>) {
    let (a, b) = edge_vertices(e);
    if e.0.orientation().is_reversed() {
        (b, a)
    } else {
        (a, b)
    }
}

/// `BRepTools_WireExplorer` order: follow `LastVertex` → next `FirstVertex`
/// (`CumOri=true`). When vertices are distinct TShapes at the same UV
/// point, `Next` with a face compares pcurve ends (`cxx:428-475`).
fn wire_edges_explorer(wire: &Wire, face: &Face) -> Vec<Edge> {
    let mut unused = wire_edges_sewd(wire);
    if unused.len() <= 1 {
        return unused;
    }
    let mut out = vec![unused.remove(0)];
    while !unused.is_empty() {
        let last = out.last().unwrap();
        let Some(lv) = edge_oriented_ends(last).1 else {
            out.extend(unused);
            break;
        };
        let by_vertex = unused.iter().position(|e| {
            edge_oriented_ends(e)
                .0
                .as_ref()
                .is_some_and(|v| Arc::ptr_eq(&v.0.tshape, &lv.0.tshape))
        });
        if let Some(i) = by_vertex {
            out.push(unused.remove(i));
            continue;
        }
        let Some((_, last_uv)) = edge_uv_ends(last, face) else {
            out.extend(unused);
            break;
        };
        match unused.iter().position(|e| {
            edge_uv_ends(e, face).is_some_and(|(u0, _)| uv_close(&u0, &last_uv))
        }) {
            Some(i) => out.push(unused.remove(i)),
            None => {
                out.extend(unused);
                break;
            }
        }
    }
    out
}

fn edge_uv_ends(e: &Edge, face: &Face) -> Option<(GpPnt2d, GpPnt2d)> {
    let pc = make_pcurve_full(e, face).ok()?;
    let fwd = Edge(e.0.oriented(Orientation::Forward));
    let (a, b) = BRepTool::edge_parameters(&fwd);
    if !a.is_finite() || !b.is_finite() {
        return None;
    }
    let (pa, pb) = (pc.d0(a), pc.d0(b));
    if e.0.orientation().is_reversed() {
        Some((pb, pa))
    } else {
        Some((pa, pb))
    }
}

fn uv_close(a: &GpPnt2d, b: &GpPnt2d) -> bool {
    let dx = a.x() - b.x();
    let dy = a.y() - b.y();
    dx * dx + dy * dy <= PCONFUSION * PCONFUSION
}

/// Pcurve samples in `WireExplorer` order (`BRepTopAdaptor_FClass2d`).
fn wire_uv_points(wire: &Wire, face: &Face) -> Vec<GpPnt2d> {
    let mut pts = Vec::new();
    for e in wire_edges_explorer(wire, face) {
        let Ok(pc) = make_pcurve_full(&e, face) else {
            continue;
        };
        let fwd = Edge(e.0.oriented(Orientation::Forward));
        let (a, b) = BRepTool::edge_parameters(&fwd);
        if !a.is_finite() || !b.is_finite() {
            continue;
        }
        let mut s = sample_pcurve(pc.as_ref(), a, b);
        if e.0.orientation().is_reversed() {
            s.reverse();
        }
        if !pts.is_empty() && !s.is_empty() {
            s.remove(0);
        }
        pts.extend(s);
    }
    pts
}

/// One-wire `BRepTopAdaptor_FClass2d` stand-in: `CSLib_Class2d` plus
/// `TabOrien` from the sampled `square` (`cxx:338`, `445`), not the
/// turning-angle `BRepMesh_Classifier` uses (coarse circles flip that sign).
fn wire_classifier(wire: &Wire, face: &Face) -> Option<(Class2d, bool, BndBox2d)> {
    let pts = wire_uv_points(wire, face);
    if pts.len() < 2 {
        return None;
    }
    let mut box2d = BndBox2d::new();
    for p in &pts {
        box2d.update_point(p.x(), p.y());
    }
    let (u0, v0, u1, v1) = box2d.get()?;
    let clas = Class2d::new(&pts, PCONFUSION, PCONFUSION, u0, v0, u1, v1);
    let mut area = 0.0;
    for i in 0..pts.len() {
        let p = pts[i];
        let q = pts[(i + 1) % pts.len()];
        area += (p.x() - q.x()) * (p.y() + q.y()) / 2.0;
    }
    let orient_ccw = area >= 0.0;
    Some((clas, orient_ccw, box2d))
}

fn classify_on_wire(clas: &Class2d, orient_ccw: bool, p: &GpPnt2d) -> PointState {
    match clas.si_dans(p) {
        Class2dResult::Uncertain => PointState::On,
        Class2dResult::Outside if orient_ccw => PointState::Out,
        Class2dResult::Inside if !orient_ccw => PointState::Out,
        _ => PointState::In,
    }
}

/// `ShapeExtend_WireData::Reverse(face)` (`cxx:483-572`): reverse the edge
/// list, reverse each edge, then `SwapSeam`. The rebuilt wire is FORWARD.
fn reverse_wire_on_face(wire: &mut Wire, face: &Face) {
    let mut edges = wire_edges_sewd(wire);
    for e in &mut edges {
        e.0.reverse();
    }
    edges.reverse();
    let b = TopoBuilder::new();
    let mut nw = b.make_wire(&edges);
    nw.0.set_closed(wire.0.closed());
    nw.0.set_orientation(Orientation::Forward);
    *wire = nw;
    swap_wire_seams(wire, face);
}

/// `ShapeFix_Face::FixOrientation` multi-wire branch (`cxx:1274-1602`):
/// classify each wire against the others; reverse when it is OUT and the
/// infinite point is IN, or IN (contained) and the infinite point is OUT.
fn fix_multiwire_orientation(face: &Face, wires: &mut [Wire]) {
    let n = wires.len();
    if n < 2 {
        return;
    }
    let mut boxes = Vec::with_capacity(n);
    let mut classifiers = Vec::with_capacity(n);
    for w in wires.iter() {
        match wire_classifier(w, face) {
            Some((c, orient, b)) => {
                classifiers.push(Some((c, orient)));
                boxes.push(b);
            }
            None => {
                classifiers.push(None);
                boxes.push(BndBox2d::new());
            }
        }
    }
    // SI: 0 unknown / 1 keep / 2 contained+inf OUT / 3 contained+inf IN
    let mut si = vec![0i32; n];
    let mut map_int: Vec<bool> = vec![false; n];
    for i in 0..n {
        let Some((clas, orient_ccw)) = classifiers[i].as_ref() else {
            continue;
        };
        // Infinite point: In on a clockwise wire (`TotCross2D < 0`).
        let staout_in = !*orient_ccw;
        let mut sta_in = false;
        let mut unknown = false;
        for j in 0..n {
            if i == j {
                continue;
            }
            if boxes[j].is_out_box(&boxes[i]) {
                continue;
            }
            let mut stb: Option<bool> = None;
            for e in wire_edges_sewd(&wires[j]) {
                let Ok(pc) = make_pcurve_full(&e, face) else {
                    continue;
                };
                let fwd = Edge(e.0.oriented(Orientation::Forward));
                let (a, b) = BRepTool::edge_parameters(&fwd);
                if !a.is_finite() || !b.is_finite() {
                    continue;
                }
                let unp = pc.d0((a + b) * 0.5);
                match classify_on_wire(clas, *orient_ccw, &unp) {
                    PointState::On => {}
                    st @ (PointState::In | PointState::Out) => {
                        let is_in = st == PointState::In;
                        match stb {
                            None => stb = Some(is_in),
                            Some(prev) if prev != is_in => {
                                unknown = true;
                                break;
                            }
                            _ => {}
                        }
                    }
                }
            }
            if unknown {
                break;
            }
            let Some(is_in) = stb else {
                continue;
            };
            if is_in == staout_in {
                sta_in = true;
            } else {
                map_int[j] = true;
            }
        }
        if unknown {
            continue;
        }
        if !sta_in {
            if staout_in {
                reverse_wire_on_face(&mut wires[i], face);
                si[i] = 1;
            } else {
                si[i] = 1;
            }
        } else if !staout_in {
            si[i] = 2;
        } else {
            si[i] = 3;
        }
    }
    for i in 0..n {
        if si[i] <= 1 {
            continue;
        }
        if !map_int[i] {
            if si[i] == 3 {
                reverse_wire_on_face(&mut wires[i], face);
            }
        } else if si[i] == 2 {
            reverse_wire_on_face(&mut wires[i], face);
        }
    }
}

/// `ShapeExtend_WireData::Reverse(face)` seam pass (`cxx:545-571`):
/// `ComputeSeams` then `SwapSeam` once per closed-surface edge. `SwapSeam`
/// returns immediately on a REVERSED occurrence (`cxx:518-520`).
fn swap_wire_seams(wire: &Wire, face: &Face) {
    let mut seen = HashMap::new();
    for e in edges_of_wire(wire) {
        if e.0.orientation().is_reversed() {
            continue;
        }
        let k = GeometryRegistry::shape_key(&e.0);
        if seen.insert(k, ()).is_none() {
            swap_seam_pcurves(&e, face);
        }
    }
}

/// `SwapSeam` (`ShapeExtend_WireData.cxx:511-543`): exchange PCurve/PCurve2
/// once, on the FORWARD occurrence only (`cxx:518-520`).
fn swap_seam_pcurves(edge: &Edge, face: &Face) {
    if edge.0.orientation().is_reversed() {
        return;
    }
    let face_key = GeometryRegistry::shape_key(&face.0);
    let mut pcs = GeometryRegistry::global().edge_pcurves(&edge.0, face_key);
    if pcs.len() < 2 {
        return;
    }
    pcs.swap(0, 1);
    GeometryRegistry::global().set_edge_pcurves(&edge.0, face_key, pcs);
}

/// `ShapeAnalysis::TotCross2D` — signed 2D area of a wire's pcurves (trapezoid
/// rule over sampled pcurve points, sequence reversed per REVERSED edge).
fn wire_area_2d(wire: &Wire, face: &Face) -> f64 {
    let mut totcross = 0.0;
    let mut uv0: Option<GpPnt2d> = None;
    let mut fuv = GpPnt2d::new(0.0, 0.0);
    let mut nbc = 0usize;
    for e in wire_edges_sewd(wire) {
        // `BRep_Tool::CurveOnSurface(edge, face)` picks PCurve2 on a reversed
        // seam (`BRep_Tool.cxx:354-357`); then TotCross2D reverses the samples.
        let Ok(pc) = make_pcurve_full(&e, face) else { continue };
        let fwd = Edge(e.0.oriented(Orientation::Forward));
        let (a, b) = BRepTool::edge_parameters(&fwd);
        if !a.is_finite() || !b.is_finite() {
            continue;
        }
        let mut pts = sample_pcurve(pc.as_ref(), a, b);
        if e.0.orientation().is_reversed() {
            pts.reverse();
        }
        nbc += 1;
        if nbc == 1 {
            fuv = pts[0];
            uv0 = Some(pts[0]);
        }
        for p in &pts {
            totcross += (fuv.x() - p.x()) * (fuv.y() + p.y()) / 2.0;
            fuv = *p;
        }
    }
    if let Some(u0) = uv0 {
        totcross += (fuv.x() - u0.x()) * (fuv.y() + u0.y()) / 2.0;
    }
    totcross
}

/// Uniform sample of a pcurve over `[a, b]` (endpoints + interior) with the
/// point count OCCT uses for both `ShapeAnalysis::TotCross2D`
/// (`ShapeAnalysis_Curve::GetSamplePoints`, `ShapeAnalysis_Curve.cxx:1317-1337`)
/// and the `BRepTopAdaptor_FClass2d` boundary polygon
/// (`BRepTopAdaptor_FClass2d.cxx:179-185`).
fn sample_pcurve(pc: &dyn Curve2d, a: f64, b: f64) -> Vec<GpPnt2d> {
    let n = sample_count_2d(pc, a, b);
    let mut pts = Vec::with_capacity(n);
    for i in 0..n {
        let t = a + (b - a) * i as f64 / (n - 1) as f64;
        pts.push(pc.d0(t));
    }
    pts
}

/// `Geom2dInt_Geom2dCurveTool::NbSamples(const Adaptor2d_Curve2d&)`
/// (`Geom2dInt_Geom2dCurveTool.cxx:73-91`) plus the `nbs *= 4` applied by the
/// two callers that share this count:
///
/// * `ShapeAnalysis_Curve.cxx:1325-1331` (`GetSamplePoints` for `TotCross2D`)
/// * `BRepTopAdaptor_FClass2d.cxx:182-184`
///
/// Both apply the `*4` at their own call site, so it stays here rather than in
/// [`crate::curve_sampling_2d::nb_samples`].
fn sample_count_2d(pc: &dyn Curve2d, a: f64, b: f64) -> usize {
    let mut nbs = crate::curve_sampling_2d::nb_samples(pc, a, b);
    if nbs > 2 {
        nbs *= 4;
    }
    // `nbPoints` never returns less than 2 (the `Geom2d_Line` arm).
    nbs.max(2)
}

/// Tool for building a discrete model from a topological shape.
/// Port of `BRepMesh_ModelBuilder` (+ the visiting logic of
/// `BRepMesh_ShapeVisitor`).
pub struct ModelBuilder;

impl ModelBuilder {
    /// Build the discrete model of `shape` under `params`.
    ///
    /// Mirrors `BRepMesh_ModelBuilder::performInternal` plus the
    /// `BRepMesh_ShapeVisitor` walk: every face gets its outer wire first
    /// (`ShapeAnalysis::OuterWire`), then its inner wires, each wire reordered
    /// into a connected chain (`ShapeAnalysis_Wire::CheckOrder`), with EXTERNAL
    /// edges skipped and shared edges deduplicated by TShape identity. Free edges
    /// (not bounding any face) are added as well.
    ///
    /// Returns `Err` when the shape is empty (void bounding box → `Message_Fail1`).
    pub fn build_model(shape: &TopoShape, params: &MeshParameters) -> Result<MeshModel, String> {
        let mut model = MeshModel::new(shape.clone());

        // Maximum size of the shape's bounding box (used by relative deflection).
        let bbox = shape_bbox(shape);
        if bbox.is_void() {
            return Err("BRepMesh_ModelBuilder::build_model: empty shape".to_string());
        }
        let max_size = if params.relative {
            box_max_dimension(&bbox).unwrap_or(0.0)
        } else {
            params.deflection.max(params.deflection_interior.max(0.0))
        };
        model.set_max_size(max_size);

        // Visit(Face): `IMeshTools_ShapeExplorer.cxx:101` stores only
        // `Oriented(FORWARD)` faces "to prevent inverse issue". A REVERSED
        // face would Compose hole wires to FORWARD and fill the hole.
        // Outer wire first; a failure on the outer wire fails the face, a
        // failure on an inner wire only marks it unoriented.
        let mut edge_index: HashMap<usize, usize> = HashMap::new();
        for f in faces_of(shape) {
            let f = Face(f.0.oriented(Orientation::Forward));
            let face_index = model.add_face(f.clone());

            let mut wires = wires_of_face(&f);
            if wires.len() > 1 {
                fix_multiwire_orientation(&f, &mut wires);
            }
            let outer = outer_of_wires(&wires, &f);
            if let Some(outer_wire) = &outer {
                if !add_wire(&mut model, face_index, outer_wire, &mut edge_index) {
                    if let Ok(fm) = model.face_mut(face_index) {
                        fm.set_status(MeshStatus::FAILURE);
                    }
                    continue;
                }
            }

            for w in &wires {
                if let Some(outer_wire) = &outer {
                    if Arc::ptr_eq(&w.0.tshape, &outer_wire.0.tshape) {
                        continue;
                    }
                }
                if !add_wire(&mut model, face_index, w, &mut edge_index) {
                    if let Ok(fm) = model.face_mut(face_index) {
                        fm.set_status(MeshStatus::UNORIENTED_WIRE);
                    }
                }
            }
        }

        // Free edges that do not bound any face.
        for e in edges_of(shape) {
            let key = Arc::as_ptr(&e.0.tshape) as usize;
            if !edge_index.contains_key(&key) {
                model.add_edge(e.clone());
            }
        }

        Ok(model)
    }
}
