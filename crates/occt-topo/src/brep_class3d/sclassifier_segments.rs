//! Segment loop of `BRepClass3d_SClassifier::Perform` (`BRepClass3d_SClassifier.cxx:
//! 254-516`) with its helpers: `FaceNormal` (`cxx:606-627`),
//! `GetNormalOnFaceBound` (`cxx:631-650`), `GetTransi` (`cxx:654-724`),
//! `GetAddToParam` (`cxx:569-602`) and the line selector
//! `BRepClass3d_BndBoxTreeSelectorLine::Accept` (`BRepClass3d_BndBoxTree.cxx:75-156`).
//!
//! The `BRepClass3d_BndBoxTree` is replaced by a linear scan over `myMapEV` that
//! applies the same `Accept` predicate. `RejectShell` and `RejectFace` always
//! return false in OCCT (`BRepClass3d_SolidExplorer.cxx:1037-1040`,
//! `:1084-1087`), so every face of every shell is visited.

use std::collections::HashMap;

use occt_core::gp::{GpDir, GpLin, GpPnt, GpPnt2d, GpVec};
use occt_core::precision::{ANGULAR, INFINITE, PCONFUSION};
use occt_geom::extrema_cc::curve_curve_extrema_all_range;
use occt_geom::extrema_pc::extrema_ext_pc_range;
use occt_geom::extrema_surf::{ExtPs, ExtremaExtFlag};
use occt_geom::GeomLine;

use crate::abs::Orientation;
use crate::boptools_2d::curve_on_surface_range;
use occt_core::bnd::BndBox;

use crate::brep_tool::BRepTool;
use crate::fclass2d::FaceState;
use crate::int_curves_face::Transition;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{edge_vertices, edges_of_wire, wires_of_face};

use super::{apply_trans, Segment, SClassifier, SolidExplorer};

/// `FaceNormal` (`cxx:606-627`): unit normal `D1U x D1V` at `(u, v)`, reversed
/// for REVERSED faces. `None` when the magnitude is not above `gp::Resolution()`.
pub(super) fn face_normal_cxx(face: &Face, u: f64, v: f64) -> Option<GpDir> {
    let surf = BRepTool::face_surface(face)?;
    let (_, du, dv) = surf.d1(u, v);
    let n = du.xyz().crossed(dv.xyz());
    let m = n.modulus();
    if m <= f64::MIN_POSITIVE {
        return None;
    }
    let mut dn = [n.x / m, n.y / m, n.z / m];
    if face.0.orientation() == Orientation::Reversed {
        dn = [-dn[0], -dn[1], -dn[2]];
    }
    GpDir::from_vec(&GpVec::new(dn[0], dn[1], dn[2])).ok()
}

/// `GetNormalOnFaceBound` (`cxx:631-650`).
fn get_normal_on_face_bound(edge: &Edge, face: &Face, param: f64) -> Option<GpDir> {
    let (c2d, first, last) = curve_on_surface_range(edge, face)?;
    if param < first || param > last {
        return None;
    }
    let uv = c2d.d0(param);
    face_normal_cxx(face, uv.x(), uv.y())
}

/// `GetTransi` (`cxx:654-724`). Returns the OCCT status (`1` OK, `0` skip,
/// `-1` probably faulty line) together with the transition.
fn get_transi(f1: &Face, f2: &Face, e: &Edge, param: f64, lin: &GpLin) -> (i32, Transition) {
    let mut tran = Transition::Tangent;
    let Some(nf1) = get_normal_on_face_bound(e, f1, param) else {
        return (-1, tran);
    };
    let Some(nf2) = get_normal_on_face_bound(e, f2, param) else {
        return (-1, tran);
    };
    let ldir = lin.direction();
    if ldir.dot(&nf1).abs() < ANGULAR || ldir.dot(&nf2).abs() < ANGULAR {
        // Line is orthogonal to a normal (`cxx:677-683`).
        return (-1, tran);
    }
    if nf1.is_parallel_tol(&nf2, ANGULAR) {
        let ang_d = nf1.dot(&ldir);
        if ang_d.abs() < ANGULAR {
            return (-1, tran);
        }
        tran = if ang_d > 0.0 {
            Transition::Out
        } else {
            Transition::In
        };
        return (1, tran);
    }
    // `gp_Dir ProjL = N ^ LDir ^ N`: projection of LDir on the plane of nf1/nf2.
    let n = nf1.xyz().crossed(nf2.xyz());
    let proj_xyz = n.crossed(ldir.xyz()).crossed(&n);
    // The constructor only fails for a zero vector, which the orthogonality
    // check above excludes; `(0, _)` is the skip status for that case.
    let Ok(proj) = GpDir::from_xyz(&proj_xyz) else {
        return (0, tran);
    };
    let fad = nf1.dot(&proj);
    let sad = nf2.dot(&proj);
    if fad < -ANGULAR && sad < -ANGULAR {
        tran = Transition::In;
    } else if fad > ANGULAR && sad > ANGULAR {
        tran = Transition::Out;
    } else {
        return (0, tran);
    }
    (1, tran)
}

/// `GetAddToParam` (`cxx:569-602`): how far the line must be prolonged so that
/// it reaches the bounding box `box_`. Returns `1e20` when a box corner lies
/// beyond `1e20` from the line origin.
fn get_add_to_param(lin: &GpLin, par: f64, box_: &BndBox) -> f64 {
    let Some((xmin, xmax, ymin, ymax, zmin, zmax)) = box_.get() else {
        return 0.0;
    };
    let loc = lin.location();
    let dir = GpVec::from_xyz(lin.direction().xyz());
    let mut best = par;
    for x in [xmin, xmax] {
        for y in [ymin, ymax] {
            for z in [zmin, zmax] {
                if (x - loc.x()).abs() >= 1.0e20
                    || (y - loc.xyz().y).abs() >= 1.0e20
                    || (z - loc.xyz().z).abs() >= 1.0e20
                {
                    return 1.0e20;
                }
                let t = GpVec::from_pnts(&loc, &GpPnt::new(x, y, z)).dot(&dir);
                if t > best {
                    best = t;
                }
            }
        }
    }
    best - par
}

/// Same shape occurrence: `TopTools_ShapeMapHasher` compares with `IsSame`,
/// i.e. TShape and location (`TopoDS_Shape.hxx:268-271`).
fn same_occurrence(a: &TopoShape, b: &TopoShape) -> bool {
    a.same_tshape(b) && a.location().is_equal(b.location())
}

/// `TopExp::MapShapesAndAncestors(shape, EDGE, FACE)` (`cxx:698-699`, `mapEF`):
/// each occurrence of an edge under a face is recorded, so a seam edge stored
/// twice in one face has that face twice (`TopExp.cxx:90-107`).
fn edge_face_ancestors(expl: &SolidExplorer) -> HashMap<usize, Vec<Face>> {
    let mut out: HashMap<usize, Vec<Face>> = HashMap::new();
    for face in &expl.faces {
        for wire in wires_of_face(face) {
            for edge in edges_of_wire(&wire) {
                out.entry(GeometryRegistry::shape_key(&edge.0))
                    .or_default()
                    .push(face.clone());
            }
        }
    }
    out
}

/// Hits of the line selector (`myEP` / `myVP`).
struct LineHits {
    /// Edge, parameter on the edge, parameter on the line.
    edges: Vec<(Edge, f64, f64)>,
    /// Vertex and its parameter on the line.
    vertices: Vec<(Vertex, f64)>,
}

/// `BRepClass3d_BndBoxTreeSelectorLine::Accept` applied to every `myMapEV`
/// object. The line is `Geom_Line(L)` loaded on `[-PConfusion, Par]`
/// (`BndBoxTree.hxx:93`).
fn select_line(map_ev: &[TopoShape], lin: &GpLin, par: f64) -> LineHits {
    let line = GeomLine::new(*lin);
    let mut hits = LineHits {
        edges: Vec::new(),
        vertices: Vec::new(),
    };
    for s in map_ev {
        if s.is_edge() {
            let e = Edge(s.clone());
            let Some(curve) = BRepTool::edge_curve(&e) else {
                continue;
            };
            let t = BRepTool::edge_tolerance(&e);
            let (f, l) = BRepTool::edge_parameters(&e);
            for pair in curve_curve_extrema_all_range(&*curve, &line, f, l, -PCONFUSION, par) {
                if pair.p1.square_distance(&pair.p2) < t * t {
                    hits.edges.push((e.clone(), pair.u1, pair.u2));
                }
            }
        } else if s.is_vertex() {
            let v = Vertex(s.clone());
            let t = BRepTool::vertex_tolerance(&v);
            let pv = BRepTool::vertex_point(&v);
            let sols = extrema_ext_pc_range(&line, &pv, -INFINITE, INFINITE);
            if let Some(first) = sols.first() {
                if first.sq_dist < t * t {
                    hits.vertices.push((v, first.u));
                }
            }
        }
    }
    hits
}

impl SClassifier {
    /// The segment loop of `Perform` (`cxx:254-516`). `myState` persists across
    /// iterations of the faulty-line loop, as in OCCT.
    pub(super) fn classify_segments(&mut self, expl: &mut SolidExplorer, p: &GpPnt, tol: f64) {
        // `isFaultyLine` and `anIndFace` (`cxx:254-255`).
        let mut faulty = true;
        let mut ind_face: i32 = 0;
        while faulty {
            // `Segment` on the first pass, `OtherSegment` afterwards (`cxx:259-266`).
            let Some(seg) = (if ind_face == 0 {
                expl.segment(p)
            } else {
                expl.other_segment(p)
            }) else {
                self.state = FaceState::Out;
                return;
            };
            // Faulty line when the face index does not advance (`cxx:268-278`).
            let cur_ind = expl.face_segment_index();
            if cur_ind > ind_face {
                ind_face = cur_ind;
            } else {
                self.state = FaceState::Out;
                return;
            }
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
            let mut parmin = f64::MAX; // `RealLast()` (`cxx:300`)
            let mut near_fault_par = f64::MAX; // (`cxx:302`)
            self.line_interference(expl, &seg, &mut parmin, &mut near_fault_par);
            self.scan_faces(expl, p, tol, &seg, &mut parmin, &mut faulty);

            // `cxx:511-515`.
            if near_fault_par != f64::MAX
                && parmin.abs() >= near_fault_par.abs() - PCONFUSION
            {
                faulty = true;
            }
        }
    }

    /// Line and edge / vertex interference (`cxx:305-361`).
    fn line_interference(
        &mut self,
        expl: &SolidExplorer,
        seg: &Segment,
        parmin: &mut f64,
        near_fault_par: &mut f64,
    ) {
        let hits = select_line(expl.map_ev(), &seg.lin, seg.par);
        if hits.edges.is_empty() && hits.vertices.is_empty() {
            return;
        }
        let mut lv_ints: Vec<TopoShape> = Vec::new();
        for (v, lp) in &hits.vertices {
            lv_ints.push(v.0.clone());
            if lp.abs() < near_fault_par.abs() {
                *near_fault_par = *lp;
            }
        }
        if hits.edges.is_empty() {
            return;
        }
        let ancestors = edge_face_ancestors(expl);
        for (e, param, lpar) in &hits.edges {
            let Some(ffs) = ancestors.get(&GeometryRegistry::shape_key(&e.0)) else {
                continue;
            };
            if ffs.len() != 2 {
                continue;
            }
            let f1 = &ffs[0];
            let f2 = &ffs[ffs.len() - 1];
            let (v1, v2) = edge_vertices(e);
            let hit_vertex = |v: &Option<Vertex>| {
                v.as_ref()
                    .is_some_and(|v| lv_ints.iter().any(|x| same_occurrence(x, &v.0)))
            };
            if hit_vertex(&v1) || hit_vertex(&v2) {
                continue;
            }
            let (status, mut tran) = get_transi(f1, f2, e, *param, &seg.lin);
            if status == 1 && lpar.abs() < parmin.abs() {
                *parmin = *lpar;
                apply_trans(*parmin, &mut tran, &mut self.state);
            } else if lpar.abs() < near_fault_par.abs() {
                *near_fault_par = *lpar;
            }
        }
    }

    /// Shell and face loop (`cxx:363-509`). Each face is intersected with the
    /// segment through its `myMapOfInter` intersector.
    fn scan_faces(
        &mut self,
        expl: &mut SolidExplorer,
        p: &GpPnt,
        tol: f64,
        seg: &Segment,
        parmin: &mut f64,
        faulty: &mut bool,
    ) {
        // `addW` / `AddW` (`cxx:380-381`).
        let add_w = (10.0 * tol).max(0.01 * seg.par);
        let nb_faces = expl.faces.len();
        'faces: for fi in 0..nb_faces {
            if *faulty {
                break;
            }
            let f = expl.faces[fi].clone();
            let key = GeometryRegistry::shape_key(&f.0);
            let Some(inter) = expl.inters.get_mut(&key) else {
                continue;
            };

            // Prolong the segment to the bounds of a finite face box (`cxx:383-397`).
            let mut add = add_w;
            let bbox = inter.bounding();
            if !bbox.is_void() && !bbox.is_whole() {
                add = add.max(get_add_to_param(&seg.lin, seg.par, &bbox));
            }
            let min_w = -add_w;
            let max_w = (seg.par * 10.0).min(seg.par + add);
            inter.perform(&seg.lin, min_w, max_w);
            if !inter.is_done() {
                continue;
            }

            if inter.nb_pnt() == 0 && inter.is_parallel() {
                // Distance from P to the surface (`cxx:404-441`).
                if let Some(surf) = BRepTool::face_surface(&f) {
                    let (u0, u1) = surf.u_range();
                    let (v0, v1) = surf.v_range();
                    let mut ext = ExtPs::new();
                    ext.initialize(surf.as_ref(), u0, u1, v0, v1, PCONFUSION, PCONFUSION);
                    ext.set_flag(ExtremaExtFlag::Min);
                    ext.perform(p);
                    if ext.is_done() && ext.nb_ext() > 0 {
                        let mut d = f64::MAX;
                        let mut indmin = 0usize;
                        for i in 1..=ext.nb_ext() {
                            let sq = ext.square_distance(i);
                            if sq < d {
                                d = sq;
                                indmin = i;
                            }
                        }
                        if indmin > 0 && d <= tol * tol {
                            let (u, v, _) = ext.point(indmin);
                            let st = inter.classify_uv_point(GpPnt2d::new(u, v));
                            if st == FaceState::In || st == FaceState::On {
                                self.state = FaceState::On;
                                self.face = Some(f.clone());
                                *parmin = 0.0;
                                break 'faces;
                            }
                        }
                    }
                }
            }

            // Intersection points (`cxx:444-488`).
            for i in 1..=inter.nb_pnt() {
                let w = inter.w_parameter(i);
                if w.abs() >= parmin.abs() - PCONFUSION {
                    continue;
                }
                *parmin = w;
                let a_state = inter.state(i);
                if parmin.abs() <= tol {
                    self.state = FaceState::On;
                    self.face = Some(f.clone());
                    break;
                } else if a_state == FaceState::In {
                    let mut tran = inter.transition(i);
                    if tran == Transition::Tangent {
                        // `cxx:463-470`: ignore this point.
                        continue;
                    }
                    apply_trans(*parmin, &mut tran, &mut self.state);
                    self.face = Some(f.clone());
                } else if a_state == FaceState::On {
                    *faulty = true;
                    break;
                }
            }
            if self.state == FaceState::On {
                break;
            }
        }
    }
}
