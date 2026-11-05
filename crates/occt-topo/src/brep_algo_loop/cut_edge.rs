//! `BRepAlgo_Loop::CutEdge` and its two statics (`BRepAlgo_Loop.cxx:80-227`,
//! `:790-947`).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `Bubble` | `cxx:80-122` | [`bubble`] |
//! | `UpdateClosedEdge` | `cxx:160-227` | [`update_closed_edge`] |
//! | `CutEdge` | `cxx:790-947` | [`cut_edge`] |

use occt_core::precision::CONFUSION;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{edge_vertices, is_equal, is_same};

/// `BRep_Tool::Parameter(V, E)`.
fn parameter(v: &TopoShape, e: &TopoShape) -> f64 {
    BRepTool::parameter_of_vertex(e, v)
}

/// `Bubble(E, Seq)` (`cxx:80-122`): remove the duplicated vertices of `Seq`
/// and order the rest by increasing parameter on `E`.
pub fn bubble(e: &TopoShape, seq: &mut Vec<TopoShape>) {
    // `cxx:82-96`: `Seq(i) == Seq(j)` is `IsEqual`, i.e. orientation included.
    let mut i = 0usize;
    while i + 1 < seq.len() {
        let mut j = i + 1;
        while j < seq.len() {
            if is_equal(&seq[i], &seq[j]) {
                seq.remove(j);
            } else {
                j += 1;
            }
        }
        i += 1;
    }

    // `cxx:98-121`: bubble sort on `BRep_Tool::Parameter(V, E)`.
    let nb_points = seq.len();
    let mut invert = true;
    while invert {
        invert = false;
        for k in 0..nb_points.saturating_sub(1) {
            if k + 1 >= seq.len() {
                break;
            }
            let v1 = seq[k].clone().oriented(Orientation::Internal);
            let v2 = seq[k + 1].clone().oriented(Orientation::Internal);
            let u1 = parameter(&v1, e);
            let u2 = parameter(&v2, e);
            if u2 < u1 {
                seq.swap(k, k + 1);
                invert = true;
            }
        }
    }
}

/// `UpdateClosedEdge(E, SV)` (`cxx:160-227`): drop from `SV` the intersection
/// vertices that sit on the closing vertex of a closed edge.
pub fn update_closed_edge(e: &TopoShape, sv: &mut Vec<TopoShape>) -> Option<TopoShape> {
    let (v1, _) = edge_vertices(&Edge(e.clone()));
    // `cxx:170` reads `BRep_Tool::Tolerance(V1)`; a null `V1` dereferences a
    // null handle there, so the port answers `0.0` instead of crashing.
    let tol = match &v1 {
        Some(v) => BRepTool::vertex_tolerance(v),
        None => 0.0,
    };

    if sv.is_empty() {
        return None;
    }

    let vb0 = sv.first().cloned().unwrap();
    let vb1 = sv.last().cloned().unwrap();
    let pc = match &v1 {
        Some(v) => BRepTool::vertex_point(v),
        None => return None,
    };

    let mut res: Option<TopoShape> = None;
    let mut on_start = false;
    let mut on_end = false;
    for (i, vb) in [&vb0, &vb1].into_iter().enumerate() {
        let p = BRepTool::vertex_point(&Vertex(vb.clone()));
        if p.is_equal_with_tolerance(&pc, tol) {
            res = Some(vb.clone());
            if i == 0 {
                on_start = true;
            } else {
                on_end = true;
            }
        }
    }

    if on_start && on_end {
        if is_same(&vb0, &vb1) {
            sv.remove(0);
            if !sv.is_empty() {
                sv.remove(sv.len() - 1);
            }
        }
    } else if on_start {
        sv.remove(0);
    } else if on_end {
        sv.remove(sv.len() - 1);
    }

    res
}

/// `CutEdge(E, VOnE, NE)` (`cxx:790-947`): split `E` at the vertices of `VOnE`.
pub fn cut_edge(
    e: &TopoShape,
    v_on_e: &[TopoShape],
    ne: &mut Vec<TopoShape>,
    face: &Face,
    _tol_conf: f64,
) {
    let we = e.clone().oriented(Orientation::Forward);

    let mut sv: Vec<TopoShape> = v_on_e.to_vec();

    // `cxx:810`: order the vertices by parameterisation.
    bubble(&we, &mut sv);

    let nb_ver = sv.len();
    if sv.is_empty() {
        ne.push(e.clone());
        return;
    }

    let (f, l) = BRepTool::edge_parameters(&Edge(we.clone()));
    let (mut vf, mut vl) = edge_vertices(&Edge(we.clone()));

    if nb_ver == 2 {
        let same_f = vf.as_ref().map(|v| is_equal(&sv[0], &v.0)).unwrap_or(false);
        let same_l = vl.as_ref().map(|v| is_equal(&sv[1], &v.0)).unwrap_or(false);
        if same_f && same_l {
            ne.push(e.clone());
            return;
        }
    }

    let closed_edge = match (&vf, &vl) {
        (Some(a), Some(b)) => is_same(&a.0, &b.0),
        _ => false,
    };
    if closed_edge {
        // `cxx:842-856`.
        if let Some(vcei) = update_closed_edge(&we, &mut sv) {
            vf = Some(Vertex(vcei.clone().oriented(Orientation::Forward)));
            vl = Some(Vertex(vcei.oriented(Orientation::Reversed)));
        }
        if let Some(v) = &vf {
            sv.insert(0, v.0.clone());
        }
        if let Some(v) = &vl {
            sv.push(v.0.clone());
        }
    } else {
        // `cxx:858-868`: make sure both ends of the edge are in `SV`.
        if let Some(v) = &vf {
            let same = sv.first().map(|x| is_same(x, &v.0)).unwrap_or(false);
            if !same {
                sv.insert(0, v.0.clone());
            }
        }
        if let Some(v) = &vl {
            let same = sv.last().map(|x| is_same(x, &v.0)).unwrap_or(false);
            if !same {
                sv.push(v.0.clone());
            }
        }
    }

    let reg = GeometryRegistry::global();
    let builder = TopoBuilder::new();

    // `cxx:870-919`.
    while !sv.is_empty() {
        while !sv.is_empty() && sv[0].orientation() != Orientation::Forward {
            sv.remove(0);
        }
        if sv.is_empty() {
            break;
        }
        let v1 = sv.remove(0);
        if sv.is_empty() {
            break;
        }
        if sv[0].orientation() == Orientation::Reversed {
            let v2 = sv.remove(0);

            // `cxx:891-916`: `WE.EmptyCopied()` keeps the curve
            // representations, which `B.Range` then re-trims to the piece.
            let mut new_edge = we.empty_copied();
            if let Some(g) = reg.edge_geom(&we) {
                reg.set_edge(&new_edge, g);
            }
            reg.copy_edge_pcurve_slots(&new_edge, &we);
            builder.add(&mut new_edge, &v1.oriented(Orientation::Forward));
            builder.add(&mut new_edge, &v2.oriented(Orientation::Reversed));

            let u1 = match &vf {
                Some(f_v) if is_same(&v1, &f_v.0) => f,
                _ => parameter(&v1.clone().oriented(Orientation::Internal), &we),
            };
            let u2 = match &vl {
                Some(l_v) if is_same(&v2, &l_v.0) => l,
                _ => parameter(&v2.clone().oriented(Orientation::Internal), &we),
            };
            reg.set_ranges_all(&new_edge, u1, u2);
            ne.push(new_edge.oriented(e.orientation()));
        }
    }

    // `cxx:921-946`: drop the pieces shorter than the 2D tolerance.
    let tol = 0.001;
    let mut k = 0usize;
    while k < ne.len() {
        let ee = Edge(ne[k].clone());
        let (fpar, lpar) = BRepTool::edge_parameters(&ee);
        if lpar - fpar <= CONFUSION {
            ne.remove(k);
            continue;
        }
        let (pf, pl) = BRepTool::uv_points(&ee, face);
        if pf.distance(&pl) <= tol && !BRepTool::is_closed(&ee.0) {
            ne.remove(k);
            continue;
        }
        k += 1;
    }
}
