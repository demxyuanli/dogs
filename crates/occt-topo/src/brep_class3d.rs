//! `BRepClass3d_SolidClassifier` / `SolidExplorer` / `SClassifier`.
//!
//! Source: `BRepClass3d_SolidClassifier.cxx`, `BRepClass3d_SolidExplorer.cxx`
//! (`OtherSegment` at 493), `BRepClass3d_SClassifier.cxx` (`Perform` at 203,
//! `Trans` at 728). A probe line from the query point toward a point inside a
//! face is intersected with every face by [`crate::int_curves_face`]. The
//! nearest IN hit's CS transition decides IN/OUT (`Trans`); |W| <= Tol is ON.

use occt_core::bnd::BndBox;
use occt_core::gp::{GpDir, GpLin, GpPnt, GpPnt2d, GpVec};
use occt_core::precision::{CONFUSION, PCONFUSION};

use crate::abs::Orientation;
use crate::brep_extrema::closest_point_on_edge;
use crate::brep_surface::{face_uv_bounds, surface_closest_params, surface_normal};
use crate::brep_tool::BRepTool;
use crate::fclass2d::{FaceState, FClass2d};
use crate::int_curves_face::{FaceIntersector, Transition};
use crate::iterator::cumulated_children;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools_full::{edges_of, faces_of, vertices_of};

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
    first_face: i32,
}

impl SolidExplorer {
    pub fn load(shape: TopoShape) -> Self {
        let faces = faces_of(&shape);
        Self {
            shape,
            faces,
            first_face: 0,
        }
    }

    pub fn shape(&self) -> &TopoShape {
        &self.shape
    }

    /// `Reject(P)` — solid without faces is the whole space.
    pub fn reject(&self, _p: &GpPnt) -> bool {
        self.faces.is_empty()
    }

    pub fn reset_segment(&mut self) {
        self.first_face = 0;
    }

    /// `Segment` / `OtherSegment`. `flag`: 0 ok, 1 ON infinite face, 2 empty,
    /// 3 ON surface but OUT of face.
    fn other_segment(&mut self, p: &GpPnt) -> Option<Segment> {
        while (self.first_face as usize) < self.faces.len() {
            let face = self.faces[self.first_face as usize].clone();
            self.first_face += 1;
            let (u1, u2, v1, v2) = face_uv_bounds(&face);
            let uv_inf = !u1.is_finite() || !u2.is_finite() || !v1.is_finite() || !v2.is_finite();
            if uv_inf {
                if let Some(surf) = BRepTool::face_surface(&face) {
                    let (su, sv) = surface_closest_params(surf.as_ref(), p, 16, 16);
                    if surf.d0(su, sv).distance(p) <= CONFUSION {
                        if let Ok(cl) = FClass2d::new(&face, CONFUSION) {
                            let st = cl.perform(GpPnt2d::new(su, sv));
                            if st == FaceState::In || st == FaceState::On {
                                return Some(Segment {
                                    lin: GpLin::from_pnt_dir(*p, GpDir::default_dir()),
                                    par: 0.0,
                                    flag: 1,
                                });
                            }
                            return Some(Segment {
                                lin: GpLin::from_pnt_dir(*p, GpDir::default_dir()),
                                par: 0.0,
                                flag: 3,
                            });
                        }
                    }
                }
            } else if (u2 - u1).abs() < PCONFUSION || (v2 - v1).abs() < PCONFUSION {
                return Some(Segment {
                    lin: GpLin::from_pnt_dir(*p, GpDir::default_dir()),
                    par: 0.0,
                    flag: 2,
                });
            }
            let Some((ap, _u, _v)) = find_a_point_in_the_face(&face) else {
                continue;
            };
            let dist = p.distance(&ap);
            if dist <= CONFUSION {
                return Some(Segment {
                    lin: GpLin::from_pnt_dir(*p, GpDir::default_dir()),
                    par: 0.0,
                    flag: 1,
                });
            }
            if let Some(surf) = BRepTool::face_surface(&face) {
                let (su, sv) = surface_closest_params(surf.as_ref(), p, 16, 16);
                if surf.d0(su, sv).distance(p) <= CONFUSION {
                    if let Ok(cl) = FClass2d::new(&face, CONFUSION) {
                        if cl.perform(GpPnt2d::new(su, sv)) == FaceState::Out {
                            return Some(Segment {
                                lin: GpLin::from_pnt_dir(*p, GpDir::default_dir()),
                                par: 0.0,
                                flag: 3,
                            });
                        }
                    }
                }
            }
            let dir = GpVec::from_pnts(p, &ap);
            let Ok(d) = GpDir::from_vec(&dir) else {
                continue;
            };
            return Some(Segment {
                lin: GpLin::from_pnt_dir(*p, d),
                par: dist,
                flag: 0,
            });
        }
        None
    }
}

fn on_vertex_or_edge(shape: &TopoShape, p: &GpPnt, tol: f64) -> bool {
    for v in vertices_of(shape) {
        if BRepTool::vertex_point(&v).distance(p) <= tol {
            return true;
        }
    }
    for e in edges_of(shape) {
        let (_, q) = closest_point_on_edge(&e, p, 32);
        if q.distance(p) <= tol {
            return true;
        }
    }
    false
}

fn find_a_point_in_the_face(face: &Face) -> Option<(GpPnt, f64, f64)> {
    let surf = BRepTool::face_surface(face)?;
    let (u1, u2, v1, v2) = finite_uv_of_face(face)?;
    let cl = FClass2d::new(face, CONFUSION).ok()?;
    const PAR_T: f64 = 0.43213918;
    for k in 0..12 {
        let s = 0.45 * (0.75_f64).powi(k);
        let u = u1 + (u2 - u1) * (0.5 + (PAR_T - 0.5) * (1.0 - 2.0 * s));
        let v = v1 + (v2 - v1) * (0.5 + (PAR_T - 0.5) * s.max(0.1));
        if cl.perform(GpPnt2d::new(u, v)) == FaceState::In {
            return Some((surf.d0(u, v), u, v));
        }
    }
    let mut ctx = crate::int_tools_full::IntToolsContext::new();
    crate::algo_tools3d::point_in_face(face, &mut ctx)
        .ok()
        .map(|(p, uv)| (p, uv.x(), uv.y()))
}

fn finite_uv_of_face(face: &Face) -> Option<(f64, f64, f64, f64)> {
    let (u1, u2, v1, v2) = face_uv_bounds(face);
    if u1.is_finite()
        && u2.is_finite()
        && v1.is_finite()
        && v2.is_finite()
        && (u2 - u1).abs() > PCONFUSION
        && (v2 - v1).abs() > PCONFUSION
    {
        return Some((u1, u2, v1, v2));
    }
    let surf = BRepTool::face_surface(face)?;
    let mut ua = f64::INFINITY;
    let mut ub = f64::NEG_INFINITY;
    let mut va = f64::INFINITY;
    let mut vb = f64::NEG_INFINITY;
    for vtx in vertices_of(&face.0) {
        let p = BRepTool::vertex_point(&vtx);
        let (u, v) = surface_closest_params(surf.as_ref(), &p, 8, 8);
        ua = ua.min(u);
        ub = ub.max(u);
        va = va.min(v);
        vb = vb.max(v);
    }
    if !ua.is_finite() || !ub.is_finite() {
        return None;
    }
    if ub - ua < PCONFUSION {
        ub = ua + 1.0;
    }
    if vb - va < PCONFUSION {
        vb = va + 1.0;
    }
    Some((ua, ub, va, vb))
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
        if on_vertex_or_edge(expl.shape(), p, tol) {
            self.state = FaceState::On;
            return;
        }
        expl.reset_segment();
        self.state = FaceState::Out;
        let mut faulty = true;
        let mut tries = 0;
        while faulty && tries < expl.faces.len() + 2 {
            tries += 1;
            let Some(seg) = expl.other_segment(p) else {
                self.state = FaceState::Out;
                return;
            };
            if seg.flag == 1 {
                self.state = FaceState::On;
                return;
            }
            if seg.flag == 2 {
                self.state = FaceState::Out;
                return;
            }
            if seg.flag == 3 {
                continue;
            }
            faulty = false;
            let mut parmin = f64::MAX;
            let add_w = (10.0 * tol).max(0.01 * seg.par);
            for f in &expl.faces {
                let mut inter = FaceIntersector::new(f.clone(), tol, true, true);
                let box_add = add_to_param(&seg.lin, seg.par, &inter.bounding());
                let add = add_w.max(box_add);
                let min_w = -add_w;
                let max_w = (seg.par * 10.0).min(seg.par + add);
                inter.perform(&seg.lin, min_w, max_w);
                if !inter.is_done() {
                    continue;
                }
                for i in 1..=inter.nb_pnt() {
                    let w = inter.w_parameter(i);
                    if w.abs() >= parmin.abs() - PCONFUSION {
                        continue;
                    }
                    parmin = w;
                    let st = inter.state(i);
                    if parmin.abs() <= tol && inter.pnt(i).distance(p) <= tol {
                        self.state = FaceState::On;
                        self.face = Some(f.clone());
                        return;
                    }
                    if st == FaceState::In {
                        let mut tran = inter.transition(i);
                        if tran == Transition::Tangent {
                            continue;
                        }
                        apply_trans(parmin, &mut tran, &mut self.state);
                        self.face = Some(f.clone());
                    } else if st == FaceState::On {
                        faulty = true;
                        break;
                    }
                }
                if self.state == FaceState::On {
                    return;
                }
                if faulty {
                    break;
                }
            }
        }
    }

    /// `PerformInfinitePoint`.
    pub fn perform_infinite_point(&mut self, expl: &mut SolidExplorer, tol: f64) {
        self.face = None;
        if expl.reject(&GpPnt::zero()) {
            self.state = FaceState::In;
            return;
        }
        self.state = FaceState::Out;
        for f in &expl.faces {
            let Some((ap, u, v)) = find_a_point_in_the_face(f) else {
                continue;
            };
            let Some(surf) = BRepTool::face_surface(f) else {
                continue;
            };
            let n = surface_normal(surf.as_ref(), u, v);
            let Ok(dn) = GpDir::from_vec(&n) else {
                continue;
            };
            let mut dn = dn;
            if f.0.orientation() == Orientation::Reversed {
                dn = dn.reversed();
            }
            let lin = GpLin::from_pnt_dir(ap, dn.reversed());
            let mut parmin = f64::MAX;
            let mut found = false;
            let mut tran_keep = Transition::Tangent;
            for g in &expl.faces {
                let mut inter = FaceIntersector::new(g.clone(), tol, true, true);
                inter.perform(&lin, f64::NEG_INFINITY, parmin);
                if !inter.is_done() {
                    continue;
                }
                for i in 1..=inter.nb_pnt() {
                    let w = inter.w_parameter(i);
                    if w >= parmin {
                        continue;
                    }
                    if inter.state(i) != FaceState::In {
                        continue;
                    }
                    let t = inter.transition(i);
                    if t == Transition::Tangent {
                        continue;
                    }
                    parmin = w;
                    tran_keep = t;
                    found = true;
                    self.face = Some(g.clone());
                }
            }
            if found {
                // `_cxx:182-195`: Out => infinite point is IN, In => OUT.
                self.state = if tran_keep == Transition::Out {
                    FaceState::In
                } else {
                    FaceState::Out
                };
                return;
            }
        }
    }
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

fn add_to_param(lin: &GpLin, par: f64, box_: &BndBox) -> f64 {
    let Some((xmin, xmax, ymin, ymax, zmin, zmax)) = box_.get() else {
        return 0.0;
    };
    if !xmin.is_finite() || !xmax.is_finite() {
        return 0.0;
    }
    let loc = lin.location();
    let d = lin.direction();
    let corners = [
        GpPnt::new(xmin, ymin, zmin),
        GpPnt::new(xmax, ymin, zmin),
        GpPnt::new(xmin, ymax, zmin),
        GpPnt::new(xmax, ymax, zmin),
        GpPnt::new(xmin, ymin, zmax),
        GpPnt::new(xmax, ymin, zmax),
        GpPnt::new(xmin, ymax, zmax),
        GpPnt::new(xmax, ymax, zmax),
    ];
    let mut tmax = 0.0;
    for c in &corners {
        let v = GpVec::from_pnts(&loc, c);
        let t = v.dot(&GpVec::from_xyz(d.xyz())) - par;
        if t > tmax {
            tmax = t;
        }
    }
    tmax.max(0.0)
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
