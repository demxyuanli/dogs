//! `BRepClass3d_SolidClassifier` / `SolidExplorer` / `SClassifier`.
//!
//! Source: `BRepClass3d_SolidClassifier.cxx`, `BRepClass3d_SolidExplorer.cxx`
//! (`OtherSegment` at 493), `BRepClass3d_SClassifier.cxx` (`Perform` at 203,
//! `Trans` at 728). A probe line from the query point toward a point inside a
//! face is intersected with every face by [`crate::int_curves_face`]. The
//! nearest IN hit's CS transition decides IN/OUT (`Trans`); |W| <= Tol is ON.

use occt_core::gp::{GpDir, GpLin, GpPnt, GpVec};
use occt_core::math_bullard::BullardGenerator;
use occt_core::precision::CONFUSION;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::fclass2d::FaceState;
use crate::int_curves_face::{FaceIntersector, Transition};
use crate::iterator::cumulated_children;
use crate::shape::{Edge, Face, TopoShape};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::faces_of;
use std::collections::HashMap;

mod other_segment;
mod sclassifier_segments;

use sclassifier_segments::face_normal_cxx;

/// `myParamOnEdge` after `InitShape` (`BRepClass3d_SolidExplorer.cxx:905`).
const PARAM_ON_EDGE_INIT: f64 = 0.512345;

/// `NB_MAX_POINTS_PER_FACE` (`BRepClass3d_SClassifier.cxx:124`).
const NB_MAX_POINTS_PER_FACE: usize = 10;

/// Probe line produced by `SolidExplorer::OtherSegment`.
struct Segment {
    lin: GpLin,
    par: f64,
    flag: i32,
}

/// `BRepClass3d_SolidExplorer`.
pub struct SolidExplorer {
    shape: TopoShape,
    faces: Vec<Face>,
    /// `myMapEV` (`BRepClass3d_SolidExplorer::Init`, `cxx:930-982`): the edges
    /// and vertices that bound the shape's non-`INTERNAL`/`EXTERNAL` faces,
    /// with `INTERNAL`/`EXTERNAL` and degenerated edges skipped. Used for the
    /// ON test in `BRepClass3d_SClassifier::Perform` (`cxx:217-227`).
    map_ev: Vec<TopoShape>,
    first_face: i32,
    /// `myMapOfInter` (`InitShape`, `cxx:921-926`): one intersector per face,
    /// keyed by the face TShape (`GeometryRegistry::shape_key`). Rebound by
    /// `OtherSegment` (`cxx:537-542`).
    inters: HashMap<usize, FaceIntersector>,
    /// `myParamOnEdge` (`InitShape`, `cxx:905`); updated by the `OtherSegment`
    /// retry ladder (`cxx:700-784`).
    param_on_edge: f64,
}

impl SolidExplorer {
    pub fn load(shape: TopoShape) -> Self {
        let faces = faces_of(&shape);
        let map_ev = Self::edge_vertex_map(&shape);
        // `myMapOfInter.Bind(Face, ptr)` (`cxx:924-925`): a repeated key is
        // overwritten, as `NCollection_DataMap::Bind` does.
        let mut inters = HashMap::new();
        for f in &faces {
            inters.insert(
                GeometryRegistry::shape_key(&f.0),
                FaceIntersector::new(f.clone(), CONFUSION, true, false),
            );
        }
        Self {
            shape,
            faces,
            map_ev,
            first_face: 0,
            inters,
            param_on_edge: PARAM_ON_EDGE_INIT,
        }
    }

    /// The ON-test edge/vertex list (see [`SolidExplorer::map_ev`]).
    pub fn map_ev(&self) -> &[TopoShape] {
        &self.map_ev
    }

    /// `BRepClass3d_SolidExplorer::Init` (`cxx:930-982`): walk the shape's
    /// faces; for each face keep its edges unless the face or the edge is
    /// `INTERNAL`/`EXTERNAL` or the edge is degenerated; `TopExp::MapShapes(aE,
    /// myMapEV)` then adds the edge **and its vertices**, deduplicated. A vertex
    /// or edge that is an internal *child* of the solid is therefore absent,
    /// which is what keeps an internal vertex from classifying the query point
    /// as ON.
    fn edge_vertex_map(shape: &TopoShape) -> Vec<TopoShape> {
        fn push_unique(out: &mut Vec<TopoShape>, s: &TopoShape) {
            if !out.iter().any(|k| k.same_tshape(s)) {
                out.push(s.clone());
            }
        }
        let mut out: Vec<TopoShape> = Vec::new();
        for f in faces_of(shape) {
            let fo = f.0.orientation();
            if fo == Orientation::Internal || fo == Orientation::External {
                continue;
            }
            for w in crate::topo_tools_full::wires_of_face(&f) {
                for e in crate::topo_tools_full::edges_of_wire(&w) {
                    let eo = e.0.orientation();
                    if eo == Orientation::Internal || eo == Orientation::External {
                        continue;
                    }
                    if BRepTool::is_degenerated(&e) {
                        continue;
                    }
                    push_unique(&mut out, &e.0);
                    let (v1, v2) = crate::topo_tools_full::edge_vertices(&e);
                    if let Some(v) = v1 {
                        push_unique(&mut out, &v.0);
                    }
                    if let Some(v) = v2 {
                        push_unique(&mut out, &v.0);
                    }
                }
            }
        }
        out
    }

    pub fn shape(&self) -> &TopoShape {
        &self.shape
    }

    /// `Reject(P)` — solid without faces is the whole space.
    pub fn reject(&self, _p: &GpPnt) -> bool {
        self.faces.is_empty()
    }

    /// `Segment(P, L, Par)` (`BRepClass3d_SolidExplorer.cxx:1095-1101`):
    /// restart the face index, then `OtherSegment`.
    fn segment(&mut self, p: &GpPnt) -> Option<Segment> {
        self.first_face = 0;
        self.other_segment(p)
    }

    /// `GetFaceSegmentIndex()` (`cxx:795-797`): the face index of the last
    /// `Segment` / `OtherSegment` call.
    pub fn face_segment_index(&self) -> i32 {
        self.first_face
    }
}

/// `BRepClass3d_SClassifier::Perform` (`cxx:217-227`): a query point within
/// tolerance of a vertex or edge of `map_ev` — the explorer's `myMapEV` — is ON.
/// The list is built by [`SolidExplorer::edge_vertex_map`], which skips internal
/// children exactly as `BRepClass3d_SolidExplorer::Init` does.
fn map_ev_accepts_point(map_ev: &[TopoShape], p: &GpPnt) -> bool {
    map_ev.iter().any(|s| {
        if s.is_vertex() {
            let v = crate::shape::Vertex(s.clone());
            let t = BRepTool::vertex_tolerance(&v);
            BRepTool::vertex_point(&v).square_distance(p) < t * t
        } else if s.is_edge() {
            edge_accepts_point(&Edge(s.clone()), p)
        } else {
            false
        }
    })
}

/// Edge branch of `BRepClass3d_BndBoxTreeSelectorPoint::Accept`: any
/// `Extrema_ExtPC` solution over the edge range closer than the edge tolerance.
fn edge_accepts_point(e: &Edge, p: &GpPnt) -> bool {
    let Some(curve) = BRepTool::edge_curve(e) else {
        return false;
    };
    let t = BRepTool::edge_tolerance(e);
    let (f, l) = BRepTool::edge_parameters(e);
    occt_geom::extrema_pc::extrema_ext_pc_range(&*curve, p, f, l)
        .iter()
        .any(|sol| sol.sq_dist < t * t)
}

/// `BRepClass3d_SClassifier`.
pub struct SClassifier {
    pub(crate) state: FaceState,
    face: Option<Face>,
}

impl SClassifier {
    pub fn new() -> Self {
        Self {
            state: FaceState::Unknown,
            face: None,
        }
    }

    pub fn state(&self) -> FaceState {
        self.state
    }

    pub fn face(&self) -> Option<&Face> {
        self.face.as_ref()
    }

    /// `Perform(Explorer, P, Tol)`.
    pub fn perform(&mut self, expl: &mut SolidExplorer, p: &GpPnt, tol: f64) {
        self.face = None;
        if expl.reject(p) {
            self.state = FaceState::In;
            return;
        }
        if map_ev_accepts_point(expl.map_ev(), p) {
            self.state = FaceState::On;
            return;
        }
        self.state = FaceState::Out;
        // Segment loop (`cxx:254-516`), see `sclassifier_segments.rs`.
        self.classify_segments(expl, p, tol);
    }

    /// `BRepClass3d_SClassifier::PerformInfinitePoint` (`cxx:82-199`).
    /// Probe lines leave random points of each face along the inverted face
    /// normal. The nearest hit (minimum parameter over all faces) decides the
    /// state by its transition. `Tol` is unused in OCCT (`cxx:83`).
    pub fn perform_infinite_point(&mut self, expl: &mut SolidExplorer, _tol: f64) {
        if expl.reject(&GpPnt::zero()) {
            // `myState = 3` (`cxx:98`).
            self.state = FaceState::In;
            return;
        }
        // `myState = 2` (`cxx:111`): stays ON when no probe decides.
        self.state = FaceState::On;
        let mut rng = BullardGenerator::new();
        for _itry in 0..NB_MAX_POINTS_PER_FACE {
            for iface in 0..expl.faces.len() {
                // `aParam = 0.1 + 0.8 * NextReal()` (`cxx:134`).
                let a_param = 0.1 + 0.8 * rng.next_real();
                let Some(lin) = probe_line(&expl.faces[iface], a_param) else {
                    continue;
                };
                let mut state = FaceState::Out;
                let mut transition = Transition::Tangent;
                let mut parmin = f64::MAX;
                // `aSE.Intersector(CurFace)` (`cxx:153`) is the per-face
                // intersector kept by the explorer (`myMapOfInter`).
                for k in 0..expl.faces.len() {
                    let key = GeometryRegistry::shape_key(&expl.faces[k].0);
                    let Some(inter) = expl.inters.get_mut(&key) else {
                        continue;
                    };
                    inter.perform(&lin, -f64::MAX, parmin);
                    if !inter.is_done() || inter.nb_pnt() <= 0 {
                        continue;
                    }
                    let mut imin = 1;
                    for i in 2..=inter.nb_pnt() {
                        if inter.w_parameter(i) < inter.w_parameter(imin) {
                            imin = i;
                        }
                    }
                    parmin = inter.w_parameter(imin);
                    state = inter.state(imin);
                    transition = inter.transition(imin);
                }
                if state == FaceState::In {
                    // `_cxx:182-195`: Out => infinite point is IN, In => OUT.
                    match transition {
                        Transition::Out => {
                            self.state = FaceState::In;
                            return;
                        }
                        Transition::In => {
                            self.state = FaceState::Out;
                            return;
                        }
                        Transition::Tangent => {}
                    }
                }
            }
        }
    }
}

/// `BRepClass3d_SClassifier::PerformInfinitePoint` probe (`cxx:134-141`):
/// `FindAPointInTheFace` with parameter `a_param`, then `gp_Lin(aPoint, -aDN)`
/// where `aDN` is `FaceNormal` (`cxx:606-627`). Returns `None` where OCCT
/// `continue`s.
fn probe_line(face: &Face, a_param: f64) -> Option<GpLin> {
    let mut ap = GpPnt::new(0.0, 0.0, 0.0);
    let mut u = 0.0_f64;
    let mut v = 0.0_f64;
    let mut d1u = GpVec::new(0.0, 0.0, 0.0);
    let mut d1v = GpVec::new(0.0, 0.0, 0.0);
    if !other_segment::find_a_point_in_the_face(
        face,
        a_param,
        &mut ap,
        &mut u,
        &mut v,
        &mut d1u,
        &mut d1v,
    ) {
        return None;
    }
    // `-aDN` (`cxx:141`), with `aDN` from `FaceNormal` (`cxx:606-627`).
    let dn = face_normal_cxx(face, u, v)?;
    let x = dn.xyz();
    let dir = GpDir::from_vec(&GpVec::new(-x.x, -x.y, -x.z)).ok()?;
    Some(GpLin::from_pnt_dir(ap, dir))
}

impl Default for SClassifier {
    fn default() -> Self {
        Self::new()
    }
}

fn apply_trans(parmin: f64, tran: &mut Transition, state: &mut FaceState) {
    if parmin < 0.0 {
        *tran = match *tran {
            Transition::Out => Transition::In,
            Transition::In => Transition::Out,
            Transition::Tangent => Transition::Tangent,
        };
    }
    *state = if *tran == Transition::Out {
        FaceState::In
    } else {
        FaceState::Out
    };
}

/// `BRepClass3d_SolidClassifier`.
pub struct SolidClassifier {
    explorer: Option<SolidExplorer>,
    inner: SClassifier,
}

impl SolidClassifier {
    pub fn new() -> Self {
        Self {
            explorer: None,
            inner: SClassifier::new(),
        }
    }

    pub fn load(&mut self, shape: TopoShape) {
        self.explorer = Some(SolidExplorer::load(shape));
        self.inner = SClassifier::new();
    }

    /// `Perform(P, Tol)`.
    pub fn perform(&mut self, p: &GpPnt, tol: f64) {
        let Some(expl) = self.explorer.as_mut() else {
            self.inner.state = FaceState::Unknown;
            return;
        };
        self.inner.perform(expl, p, tol);
    }

    /// `PerformInfinitePoint(Tol)`.
    pub fn perform_infinite_point(&mut self, tol: f64) {
        let Some(expl) = self.explorer.as_mut() else {
            self.inner.state = FaceState::Unknown;
            return;
        };
        self.inner.perform_infinite_point(expl, tol);
    }

    pub fn state(&self) -> FaceState {
        self.inner.state()
    }

    /// Classify `p` in `shape` without reusing a loaded explorer.
    pub fn classify(shape: &TopoShape, p: &GpPnt, tol: f64) -> FaceState {
        let mut sc = Self::new();
        sc.load(shape.clone());
        sc.perform(p, tol);
        sc.state()
    }
}

impl Default for SolidClassifier {
    fn default() -> Self {
        Self::new()
    }
}

/// `BRepTools::OrientClosedSolid` — reverse the solid when infinity is IN.
pub fn orient_closed_solid(shape: &mut TopoShape) -> bool {
    let mut sc = SolidClassifier::new();
    sc.load(shape.clone());
    sc.perform_infinite_point(CONFUSION);
    match sc.state() {
        FaceState::In => {
            shape.reverse();
            true
        }
        FaceState::Out | FaceState::On => true,
        FaceState::Unknown => false,
    }
}

/// Orientation of `edge` inside `face` (`BRep_Tool` / `TopExp` composed ori).
pub fn orientation_of_edge_in_face(edge: &TopoShape, face: &TopoShape) -> Option<Orientation> {
    for w in cumulated_children(face) {
        if w.shape_type() != crate::abs::ShapeType::Wire {
            continue;
        }
        for e in cumulated_children(&w) {
            if e.shape_type() == crate::abs::ShapeType::Edge && e.same_tshape(edge) {
                return Some(e.orientation());
            }
        }
    }
    None
}

/// `BOPTools_AlgoTools::OrientFacesOnShell` reverse pass: shared edges of two
/// already-placed faces must have opposite orientations.
pub fn apply_orient_faces_on_shell(shell: &mut TopoShape) {
    use std::collections::HashMap;

    use crate::abs::ShapeType;
    use crate::bop_split_seam::is_closed_on_face;
    use crate::tgeometry::GeometryRegistry;
    use crate::topo_tools_full::{edges_of_wire, wires_of_face};

    let faces: Vec<TopoShape> = shell
        .tshape
        .read()
        .unwrap()
        .children
        .iter()
        .filter(|h| h.shape_type() == ShapeType::Face)
        .cloned()
        .collect();
    if faces.len() < 2 {
        return;
    }
    let mut ef: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut edge_of: HashMap<usize, TopoShape> = HashMap::new();
    for (fi, f) in faces.iter().enumerate() {
        let face = Face(f.clone());
        for w in wires_of_face(&face) {
            for e in edges_of_wire(&w) {
                if BRepTool::is_degenerated(&e) {
                    continue;
                }
                let k = GeometryRegistry::shape_key(&e.0);
                edge_of.entry(k).or_insert_with(|| e.0.clone());
                let v = ef.entry(k).or_default();
                if !v.contains(&fi) {
                    v.push(fi);
                }
            }
        }
    }
    let mut processed = vec![false; faces.len()];
    let mut out: Vec<TopoShape> = Vec::new();
    let mut working = faces.clone();
    for (ek, fis) in &ef {
        if fis.len() != 2 {
            continue;
        }
        let i1 = fis[0];
        let i2 = fis[1];
        let p1 = processed[i1];
        let p2 = processed[i2];
        if p1 && p2 {
            continue;
        }
        let Some(edge) = edge_of.get(ek) else {
            continue;
        };
        if !p1 && !p2 {
            processed[i1] = true;
            out.push(working[i1].clone());
        }
        let f1 = &working[i1];
        let f2 = &working[i2];
        let Some(or1) = orientation_of_edge_in_face(edge, f1) else {
            continue;
        };
        let Some(or2) = orientation_of_edge_in_face(edge, f2) else {
            continue;
        };
        if p1 && !p2 {
            if or1 == or2 {
                let ff1 = Face(f1.clone());
                let ff2 = Face(f2.clone());
                if !is_closed_on_face(&Edge(edge.clone()), &ff1)
                    && !is_closed_on_face(&Edge(edge.clone()), &ff2)
                {
                    working[i2].reverse();
                }
            }
            processed[i2] = true;
            out.push(working[i2].clone());
        } else if !p1 && p2 {
            if or1 == or2 {
                let ff1 = Face(f1.clone());
                let ff2 = Face(f2.clone());
                if !is_closed_on_face(&Edge(edge.clone()), &ff1)
                    && !is_closed_on_face(&Edge(edge.clone()), &ff2)
                {
                    working[i1].reverse();
                }
            }
            processed[i1] = true;
            out.push(working[i1].clone());
        }
    }
    for (i, f) in working.iter().enumerate() {
        if !processed[i] {
            out.push(f.clone());
        }
    }
    if let Ok(mut t) = shell.tshape.write() {
        t.children = out;
    }
}
