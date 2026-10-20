//! Edge amplification/remeshing (`IncrementalMesh::amplify_remesh_edges`).

use super::*;

impl IncrementalMesh {
    /// `BRepMesh_ModelHealer::EdgeAmplifier` (`cxx:40-84`).
    pub(super) fn amplify_remesh_edges(
        &self,
        model: &mut MeshModel,
        edges: &std::collections::HashSet<usize>,
    ) -> Result<(), String> {
        let mut list: Vec<usize> = edges.iter().copied().collect();
        list.sort_unstable();
        for ei in list {
            self.amplify_one_edge(model, ei)?;
        }
        Ok(())
    }

    /// One `EdgeAmplifier::operator()` (`cxx:50-83`): keep endpoints, insert
    /// interiors, `MinSize = Confusion`, min points = previous `ParametersNb`
    /// (plus one when a face outer wire has at most two edges).
    pub(super) fn amplify_one_edge(&self, model: &mut MeshModel, edge_index: usize) -> Result<(), String> {
        let mut points_nb = model.edge(edge_index)?.discretization().parameters_nb();
        let n_pc = model.edge(edge_index)?.pcurves_nb();
        for p in 0..n_pc {
            let face_index = model.edge(edge_index)?.pcurve(p)?.face();
            let face = model.face(face_index)?;
            if face.wires_nb() == 0 {
                continue;
            }
            let wi = face.wire(0)?;
            if model.wire(wi)?.edges_nb() <= 2 {
                points_nb += 1;
                break;
            }
        }
        let nd = (model.edge(edge_index)?.deflection() / 3.0).max(occt_core::precision::CONFUSION);
        {
            let e = model.edge_mut(edge_index)?;
            e.set_deflection(nd);
            e.set_status(MeshStatus::REMESH);
            e.clear(true);
        }
        let min_size = occt_core::precision::CONFUSION;
        let topo_edge = Edge(model.edge(edge_index)?.edge().0.oriented(Orientation::Forward));
        let (a, b) = BRepTool::edge_parameters(&topo_edge);
        if !(a.is_finite() && b.is_finite() && b - a >= 1e-15) {
            return Ok(());
        }
        let curve = model.edge(edge_index)?.curve().ok_or_else(|| {
            format!("EdgeAmplifier: edge {edge_index} has no 3D curve")
        })?;
        let (mut def, ang, degenerated, is_free, same_param, same_range, ori) = {
            let e = model.edge(edge_index)?;
            (
                (0.5 * e.deflection()).max(1e-9),
                0.5 * e.angular_deflection(),
                e.is_degenerated(),
                e.is_free(),
                e.same_param(),
                e.same_range(),
                e.edge().0.orientation(),
            )
        };
        if ori == Orientation::Internal {
            def = (0.5 * def).max(1e-9);
        }
        let mut cos_value: Option<(
            std::sync::Arc<dyn occt_geom2d::curve::Curve2d>,
            std::sync::Arc<dyn occt_geom::Surface>,
        )> = None;
        let (tess_curve, tess_first, tess_last) = if !is_free && !same_param {
            let (face_index, orientation) = {
                let pc = model.edge(edge_index)?.pcurve(0)?;
                (pc.face(), pc.orientation())
            };
            let topo_face = model.face(face_index)?.face().clone();
            let oriented = Edge(model.edge(edge_index)?.edge().0.oriented(orientation));
            match (
                crate::pcurve_full::make_pcurve_full(&oriented, &topo_face),
                BRepTool::face_surface_world(&topo_face),
            ) {
                (Ok(c2), Some(s)) => {
                    let (pf, pl) = crate::boptools_2d::curve_on_surface_range(&oriented, &topo_face)
                        .map(|(_, f, l)| (f, l))
                        .unwrap_or((a, b));
                    cos_value = Some((c2.clone(), s.clone()));
                    (
                        std::sync::Arc::new(CurveOnSurface::new(c2, s, pf, pl))
                            as std::sync::Arc<dyn occt_geom::Curve>,
                        pf,
                        pl,
                    )
                }
                _ => (curve.clone(), a, b),
            }
        } else {
            (curve.clone(), a, b)
        };
        let min_pts = tessellator_min_points(tess_curve.as_ref()).max(points_nb.max(2));
        let tess_edge = model.edge(edge_index)?.edge().clone();
        let mut tess = CurveTessellator::from_range_angular_min(
            tess_curve,
            tess_first,
            tess_last,
            def,
            ang,
            min_pts,
            min_size,
        );
        tess.add_internal_vertices(&tess_edge);
        let mut tess_params = tess.params().to_vec();
        if !is_free && same_param && same_range && tess_params.len() > 1 {
            let mut pcs = Vec::new();
            for p in 0..model.edge(edge_index)?.pcurves_nb() {
                let (face_index, orientation) = {
                    let pc = model.edge(edge_index)?.pcurve(p)?;
                    (pc.face(), pc.orientation())
                };
                let topo_face = model.face(face_index)?.face().clone();
                let oriented = Edge(model.edge(edge_index)?.edge().0.oriented(orientation));
                if let (Ok(c2), Some(s)) = (
                    crate::pcurve_full::make_pcurve_full(&oriented, &topo_face),
                    BRepTool::face_surface_world(&topo_face),
                ) {
                    pcs.push((c2, s));
                }
            }
            split_by_deflection2d(curve.as_ref(), &mut tess_params, &pcs, def, min_size);
        }
        if degenerated {
            return Ok(());
        }
        let edge_tol = BRepTool::edge_tolerance(&model.edge(edge_index)?.edge());
        let n = tess_params.len();
        for k in 1..n.saturating_sub(1) {
            let t = tess_params[k];
            if let Some((c2, s)) = cos_value.as_ref() {
                let p = curve.d0(t);
                if !curve_tessellator_value_ok(c2.as_ref(), s.as_ref(), t, &p, edge_tol) {
                    continue;
                }
            }
            let p = curve.d0(t);
            let pos = model.edge(edge_index)?.discretization().parameters_nb().saturating_sub(1);
            model
                .edge_mut(edge_index)?
                .discretization_mut()
                .insert_point(pos, p, t)?;
        }
        // `Tessellate2d(false)` (`cxx:307-325`): interiors only, keep snapped ends.
        let params3d = model.edge(edge_index)?.discretization().parameters().to_vec();
        let pts3d = model.edge(edge_index)?.discretization().points().to_vec();
        let n3 = params3d.len();
        if n3 < 3 {
            return Ok(());
        }
        let n_pc = model.edge(edge_index)?.pcurves_nb();
        for p in 0..n_pc {
            let face_index = model.edge(edge_index)?.pcurve(p)?.face();
            let orientation = model.edge(edge_index)?.pcurve(p)?.orientation();
            let topo_face = model.face(face_index)?.face().clone();
            let oriented = Edge(model.edge(edge_index)?.edge().0.oriented(orientation));
            let Ok(pc) = crate::pcurve_full::make_pcurve_full(&oriented, &topo_face) else {
                continue;
            };
            let Some(surface) = BRepTool::face_surface_world(&topo_face) else {
                continue;
            };
            let (pf, pl) = crate::boptools_2d::curve_on_surface_range(&oriented, &topo_face)
                .map(|(_, f, l)| (f, l))
                .filter(|(f, l)| f.is_finite() && l.is_finite())
                .unwrap_or((a, b));
            let mut provider = if same_param {
                EdgeParameterProvider::new(pf, pl)
            } else {
                let old_first = params3d.first().copied().unwrap_or(a);
                let old_last = params3d.last().copied().unwrap_or(b);
                EdgeParameterProvider::with_stored(pf, pl, old_first, old_last)
            };
            let cos = CurveOnSurface::new(pc.clone(), surface.clone(), pf, pl);
            for k in 1..n3 - 1 {
                let t = params3d[k];
                let u = provider.parameter_of(t, &pts3d[k], &cos);
                let pos = model.edge(edge_index)?.pcurve(p)?.parameters_nb().saturating_sub(1);
                model
                    .edge_mut(edge_index)?
                    .pcurve_mut(p)?
                    .insert_point(pos, pc.d0(u), t)?;
            }
        }
        Ok(())
    }

}
