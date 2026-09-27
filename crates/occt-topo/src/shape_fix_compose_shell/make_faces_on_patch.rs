//! ShapeFix_ComposeShell::MakeFacesOnPatch (ShapeFix_ComposeShell.cxx:2978-3271).

use std::sync::Arc;

use occt_core::gp::GpPnt2d;
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::Surface;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::fclass2d::{FaceState, FClass2d};
use crate::pcurve_full::surface_value_of_uv;
use crate::shape::{Edge, TopoShape, Vertex, Wire};
use crate::shhealing::transfer_params::copy_nm_vertex_face;
use crate::topo_tools_full::edges_of_wire;

use super::reshape::ReShape;
use super::shell::ComposeShell;

fn is_fwd_or_rev(o: Orientation) -> bool {
    matches!(o, Orientation::Forward | Orientation::Reversed)
}

/// The first edge of `wire` whose orientation is FORWARD or REVERSED.
fn first_oriented_edge(wire: &Wire) -> Option<Edge> {
    edges_of_wire(wire).into_iter().find(|e| is_fwd_or_rev(e.0.orientation()))
}

impl ComposeShell {
    /// `MakeFacesOnPatch(faces, surf, loops)` (`cxx:2978-3271`).
    pub fn make_faces_on_patch(
        &mut self,
        faces: &mut Vec<TopoShape>,
        surf: &Arc<dyn Surface>,
        loops: &mut Vec<TopoShape>,
    ) {
        let builder = TopoBuilder::new();
        let my_face = match self.face() {
            Some(f) => f.clone(),
            None => return,
        };

        // cxx:2984-3010: single loop.
        if loops.len() == 1 {
            if loops[0].shape_type() != ShapeType::Wire {
                return;
            }
            let wire = Wire(loops[0].clone());
            let new_face = builder.make_face(surf.clone(), &[wire]);
            if self.invert_edge_status() {
                // UNPORTED (cxx:2997-3006): `ShapeFix_Face::FixOrientation`.
                // `ShapeFix_ComposeShell::Perform` clears `myInvertEdgeStatus`
                // (`cxx:209`) before `MakeFacesOnPatch`, so the branch is not
                // taken on the `FixMissingSeam` path.
            }
            faces.push(new_face.0);
            return;
        }

        // cxx:3014-3016: pseudo-face.
        let pf = builder.make_face(surf.clone(), &[]);
        let surf_ref = BRepTool::face_surface(&pf).unwrap_or_else(|| surf.clone());
        let mut roots: Vec<TopoShape> = Vec::new();

        // cxx:3021-3143: find roots.
        for i in 0..loops.len() {
            let a_shape = loops[i].clone();
            if a_shape.shape_type() != ShapeType::Wire || !is_fwd_or_rev(a_shape.orientation()) {
                continue;
            }
            let wr = Wire(a_shape.clone());
            let ew = edges_of_wire(&wr);
            let mut k = 0usize;
            while k < ew.len() && !is_fwd_or_rev(ew[k].0.orientation()) {
                k += 1;
            }
            if k >= ew.len() {
                continue;
            }
            let (cw, cf, cl) = match boptools_2d::curve_on_surface_range(&ew[k], &pf) {
                Some(v) => v,
                None => continue,
            };
            let unp = cw.d0(0.5 * (cf + cl));

            let mut j = 0usize;
            while j < loops.len() {
                if i == j {
                    j += 1;
                    continue;
                }
                let a_shape2 = loops[j].clone();
                if a_shape2.shape_type() != ShapeType::Wire || !is_fwd_or_rev(a_shape2.orientation()) {
                    j += 1;
                    continue;
                }
                let w1 = Wire(a_shape2);
                let w1_edges: Vec<Edge> = edges_of_wire(&w1)
                    .into_iter()
                    .filter(|e| is_fwd_or_rev(e.0.orientation()))
                    .collect();
                if w1_edges.is_empty() {
                    j += 1;
                    continue;
                }
                let awtmp = builder.make_wire(&w1_edges);
                let fc = builder.make_face(surf.clone(), &[awtmp]);
                let clas = match FClass2d::new(&fc, PCONFUSION) {
                    Ok(c) => c,
                    Err(_) => {
                        j += 1;
                        continue;
                    }
                };
                let mut st_point = clas.perform(unp);
                if st_point == FaceState::On || st_point == FaceState::Unknown {
                    // cxx:3104-3130.
                    let mut eidx = k;
                    let mut a_cw = cw.clone();
                    loop {
                        st_point = clas.perform(a_cw.d0(cl));
                        eidx += 1;
                        if eidx >= ew.len() {
                            break;
                        }
                        if !is_fwd_or_rev(ew[eidx].0.orientation()) {
                            continue;
                        }
                        if let Some((c2d, _a, _b)) = boptools_2d::curve_on_surface_range(&ew[eidx], &pf) {
                            a_cw = c2d;
                        }
                    }
                }
                let st_infin = clas.perform_infinite_point();
                if st_point != st_infin {
                    break;
                }
                j += 1;
            }
            if j >= loops.len() {
                roots.push(a_shape);
            }
        }

        // cxx:3146-3156: remove roots from loops.
        let mut i = 0usize;
        while i < loops.len() {
            let mut removed = false;
            for r in &roots {
                if crate::topo_tools_full::is_same(&loops[i], r) {
                    loops.remove(i);
                    removed = true;
                    break;
                }
            }
            if !removed {
                i += 1;
            }
        }

        // cxx:3159-3170: lost wires become roots.
        if roots.is_empty() && !loops.is_empty() {
            roots.extend(loops.drain(..));
        }

        // cxx:3173-3270: iterate on loops.
        let n_roots_start = roots.len();
        let mut ri = 0usize;
        while ri < roots.len() {
            let reverse;
            let wire = Wire(roots[ri].clone());
            let fc = builder.make_face(surf.clone(), &[wire.clone()]);
            let clas = match FClass2d::new(&fc, PCONFUSION) {
                Ok(c) => c,
                Err(_) => {
                    ri += 1;
                    continue;
                }
            };
            reverse = clas.perform_infinite_point() == FaceState::In; // cxx:3181-3188

            // cxx:3190-3229: find holes.
            let mut holes: Vec<TopoShape> = Vec::new();
            let mut j = 0usize;
            while j < loops.len() {
                let unp: Option<GpPnt2d> = match loops[j].shape_type() {
                    ShapeType::Wire => {
                        let bw = Wire(loops[j].clone());
                        match first_oriented_edge(&bw) {
                            Some(ed) => match boptools_2d::curve_on_surface_range(&ed, &pf) {
                                Some((cw, cf, cl)) => Some(cw.d0(0.5 * (cf + cl))),
                                None => None,
                            },
                            None => None,
                        }
                    }
                    ShapeType::Vertex => {
                        let a_v = Vertex(loops[j].clone());
                        let a_p = BRepTool::vertex_point(&a_v);
                        Some(surface_value_of_uv(surf_ref.as_ref(), &a_p, CONFUSION))
                    }
                    _ => None,
                };
                let unp = match unp {
                    Some(p) => p,
                    None => {
                        j += 1;
                        continue;
                    }
                };
                let state = clas.perform(unp);
                if (state == FaceState::Out) == reverse {
                    holes.push(loops[j].clone());
                    loops.remove(j);
                } else {
                    j += 1;
                }
            }

            // cxx:3232-3249.
            let mut new_face = builder.make_face(surf.clone(), &[wire]);
            for h in &holes {
                if h.shape_type() == ShapeType::Vertex {
                    let a_v = Vertex(h.clone());
                    if let Some(a_new_v) = copy_nm_vertex_face(&a_v, &new_face, &my_face) {
                        self.context_mut().replace(&a_v.0, &a_new_v.0);
                        builder.add(&mut new_face.0, &a_new_v.0);
                    }
                } else {
                    builder.add(&mut new_face.0, h);
                }
            }
            faces.push(new_face.0);

            // cxx:3252-3269.
            if ri + 1 == roots.len() && !loops.is_empty() {
                let mut extra: Vec<TopoShape> = Vec::new();
                for s in loops.iter() {
                    if s.shape_type() == ShapeType::Wire && is_fwd_or_rev(s.orientation()) {
                        extra.push(s.clone());
                    }
                }
                roots.extend(extra);
                loops.clear();
                let _ = n_roots_start;
            }
            ri += 1;
        }
    }
}
