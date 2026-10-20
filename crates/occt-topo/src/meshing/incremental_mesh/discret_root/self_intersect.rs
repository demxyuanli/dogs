//! Self-intersecting wire healing helpers (`BRepMesh_ModelHealer`).

use super::*;

impl IncrementalMesh {
    /// `BRepMesh_ModelHealer::process` FaceChecker + `amplifyEdges`
    /// (`BRepMesh_ModelHealer.cxx:234-211`).
    pub(super) fn heal_self_intersecting_wires(&self, model: &mut MeshModel) -> Result<(), String> {
        let mut intersecting: HashMap<usize, std::collections::HashSet<usize>> = HashMap::new();
        for i in 0..model.faces_nb() {
            if let Some(edges) = self.face_intersecting_edges(model, i)? {
                if crate::debug_flags::mbdiag() {
                    let mut v: Vec<usize> = edges.iter().copied().collect();
                    v.sort_unstable();
                    eprintln!("MBDIAG selfint face={i} edges={v:?}");
                }
                intersecting.insert(i, edges);
            }
        }

        const AMP_ITERS: usize = 5;
        for _ in 0..AMP_ITERS {
            let mut edges: std::collections::HashSet<usize> = std::collections::HashSet::new();
            for set in intersecting.values() {
                edges.extend(set.iter().copied());
            }
            if edges.is_empty() {
                break;
            }
            // `EdgeAmplifier` (`cxx:50-83`): remesh only the intersecting
            // edges. Do not run `discretize_pcurves` on the whole model --
            // that rebuilds every pcurve with healer MinSize and wipes
            // `connectClosestPoints` ends (`Tessellate2d(false)` keeps them).
            self.amplify_remesh_edges(model, &edges)?;
            let faces: Vec<usize> = intersecting.keys().copied().collect();
            intersecting.clear();
            for i in faces {
                ModelHealer::fix_face_boundaries(model, i)?;
                if let Some(edges) = self.face_intersecting_edges(model, i)? {
                    intersecting.insert(i, edges);
                }
            }
        }

        for i in intersecting.keys().copied() {
            let f = model.face_mut(i)?;
            f.set_status(MeshStatus::SELF_INTERSECTING_WIRE);
            f.set_status(MeshStatus::FAILURE);
        }
        Ok(())
    }

    /// `BRepMesh_FaceChecker::Perform` plus the 2-edge / 2-point branch
    /// (`BRepMesh_ModelHealer.cxx:248-278`).
    pub(super) fn face_intersecting_edges(
        &self,
        model: &MeshModel,
        face_index: usize,
    ) -> Result<Option<std::collections::HashSet<usize>>, String> {
        // `BRepMesh_FaceChecker` / `SegmentsFiller` (`cxx:61-81`) walks the
        // discretized pcurve points, not a re-projected 3D polyline.
        let uv = match Self::face_pcurve_uv_face(model, face_index) {
            Some(uv) => uv,
            None => return Ok(None),
        };
        let edge_pts = Self::face_pcurve_edge_polylines(model, face_index);
        let mut checker = FaceChecker::from_pcurve_edges(&uv, &edge_pts, 1e-7);
        let mut edges = std::collections::HashSet::new();
        // `BRepMesh_ModelHealer.cxx:248-278`: FaceChecker failure marks only
        // GetIntersectingEdges; the 2-edge / 2-point branch marks those two.
        if !checker.perform() {
            edges.extend(checker.intersecting_edges().iter().copied());
        } else if let Some((e0, e1)) = Self::two_edge_two_point_wire(model, face_index) {
            edges.insert(e0);
            edges.insert(e1);
        }
        if edges.is_empty() {
            Ok(None)
        } else {
            Ok(Some(edges))
        }
    }

    /// One UV chain per wire from stored pcurve samples.
    /// Source: `BRepMesh_FaceChecker.cxx:61-81`.
    pub(super) fn face_pcurve_uv_face(model: &MeshModel, face_index: usize) -> Option<UvFace> {
        let face = model.face(face_index).ok()?;
        let mut wires: Vec<Vec<GpPnt2d>> = Vec::new();
        for &wi in face.wires() {
            let wire = model.wire(wi).ok()?;
            let mut chain = Vec::new();
            for j in 0..wire.edges_nb() {
                let ei = wire.edge(j).ok()?;
                let ori = wire.edge_orientation(j).ok()?;
                let reverse_walk = wire.edge_reverse_walk(j);
                let edge = model.edge(ei).ok()?;
                let Some(pc) = edge.pcurve_for(face_index, ori) else {
                    continue;
                };
                let pts = pc.points();
                let reverse = !pc.is_forward() ^ reverse_walk;
                let seq: Vec<GpPnt2d> = if reverse {
                    pts.iter().rev().copied().collect()
                } else {
                    pts.to_vec()
                };
                for p in seq {
                    if chain
                        .last()
                        .map(|q: &GpPnt2d| q.distance(&p) <= 1e-9)
                        .unwrap_or(false)
                    {
                        continue;
                    }
                    chain.push(p);
                }
            }
            if chain.len() >= 2 {
                wires.push(chain);
            }
        }
        let outer_wire = wires.first()?.clone();
        let inner_wires = wires.into_iter().skip(1).collect();
        Some(UvFace {
            outer_wire,
            inner_wires,
            deflection: face.deflection(),
            id: face_index,
        })
    }

    /// One polyline per edge pcurve in `GetPoint` order
    /// (`BRepMesh_FaceChecker.cxx:65-81`).
    pub(super) fn face_pcurve_edge_polylines(
        model: &MeshModel,
        face_index: usize,
    ) -> Vec<(usize, usize, Vec<GpPnt2d>)> {
        let mut out = Vec::new();
        let Ok(face) = model.face(face_index) else {
            return out;
        };
        for (wi_local, &wi) in face.wires().iter().enumerate() {
            let Ok(wire) = model.wire(wi) else {
                continue;
            };
            for j in 0..wire.edges_nb() {
                let Ok(ei) = wire.edge(j) else {
                    continue;
                };
                let Ok(ori) = wire.edge_orientation(j) else {
                    continue;
                };
                let Ok(edge) = model.edge(ei) else {
                    continue;
                };
                let Some(pc) = edge.pcurve_for(face_index, ori) else {
                    continue;
                };
                let pts = pc.points().to_vec();
                if pts.len() >= 2 {
                    out.push((wi_local, ei, pts));
                }
            }
        }
        out
    }

    pub(super) fn two_edge_two_point_wire(model: &MeshModel, face_index: usize) -> Option<(usize, usize)> {
        let face = model.face(face_index).ok()?;
        if face.wires_nb() != 1 {
            return None;
        }
        let wi = face.wire(0).ok()?;
        let w = model.wire(wi).ok()?;
        if w.edges_nb() != 2 {
            return None;
        }
        let (e0, e1) = (w.edge(0).ok()?, w.edge(1).ok()?);
        let (o0, o1) = (w.edge_orientation(0).ok()?, w.edge_orientation(1).ok()?);
        let pc0 = model
            .edge(e0)
            .ok()
            .and_then(|e| e.pcurve_for(face_index, o0).map(|p| p.parameters_nb()));
        let pc1 = model
            .edge(e1)
            .ok()
            .and_then(|e| e.pcurve_for(face_index, o1).map(|p| p.parameters_nb()));
        if pc0 == Some(2) && pc1 == Some(2) {
            Some((e0, e1))
        } else {
            None
        }
    }

}
