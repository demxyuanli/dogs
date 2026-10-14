//! `StepToTopoDS_GeometricTool` (`IsSeamCurve` / `IsLikeSeam` / `PCurve`) with `TranslateEdge::MakePCurve` and the `TranslateEdgeLoop` seam-pcurve ordering (`cxx:699-734`).
use super::*;

impl<'a> Resolver<'a> {
    /// `StepToTopoDS_TranslateFace.cxx:479-491` — SolidWorks torus with a
    /// negative major radius is treated as a reversed `Face_Surface`.
    pub(in crate::step) fn step_surface_is_reversed(&self, surf_id: usize) -> bool {
        let Ok(rec) = self.record(surf_id) else {
            return false;
        };
        match rec.type_name.as_str() {
            "RECTANGULAR_TRIMMED_SURFACE" => rec
                .args
                .get(1)
                .and_then(|s| parse_ref(s))
                .map(|basis| self.step_surface_is_reversed(basis))
                .unwrap_or(false),
            "TOROIDAL_SURFACE" => rec
                .args
                .get(2)
                .and_then(|s| parse_f64(s).ok())
                .map(|maj| maj < 0.0)
                .unwrap_or(false),
            _ => false,
        }
    }

    /// True when `assoc_ref` is a `PCURVE` whose basis surface is `surf_ref`
    /// (`StepToTopoDS_GeometricTool.cxx:57-63` / `:94-97` / `:143-146`).
    pub(super) fn pcurve_lies_on(&self, assoc_ref: usize, surf_ref: usize) -> bool {
        matches!(self.resolve_pcurve(assoc_ref), Ok((basis, _)) if basis == surf_ref)
    }

    /// The raw STEP origin and direction ratios of a `PCURVE` whose underlying
    /// 2D curve is a `LINE` (`StepToTopoDS_GeometricTool.cxx:167-181`), as
    /// `(Pnt.X, Pnt.Y, Dir.X, Dir.Y)`. `None` for any other entity.
    pub(super) fn step_2d_line(&self, assoc_ref: usize) -> Option<(f64, f64, f64, f64)> {
        let rec = self.record(assoc_ref).ok()?;
        if rec.type_name != "PCURVE" {
            return None;
        }
        let dri = parse_ref(&rec.args[2])?;
        let dri_rec = self.record(dri).ok()?;
        if dri_rec.type_name != "DEFINITIONAL_REPRESENTATION" {
            return None;
        }
        let line_ref = *parse_ref_list(&dri_rec.args[1]).first()?;
        let line = self.record(line_ref).ok()?;
        if line.type_name != "LINE" {
            return None;
        }
        let pnt = self.record(parse_ref(&line.args[1])?).ok()?;
        let dir = self.record(parse_ref(&line.args[2])?).ok()?;
        if pnt.type_name != "CARTESIAN_POINT" || dir.type_name != "DIRECTION" {
            return None;
        }
        let p = parse_xy(&pnt.args[1]).ok()?;
        let d = parse_xy(&dir.args[1]).ok()?;
        Some((p.0, p.1, d.0, d.1))
    }

    /// `StepToTopoDS_GeometricTool::IsSeamCurve` (`cxx:78-121`). `nb_oe` is how
    /// many oriented edges of the loop reference the same step edge.
    pub(super) fn is_seam_curve(&self, curve_ref: usize, surf_ref: usize, nb_oe: usize) -> bool {
        if let Ok(rec) = self.record(curve_ref) {
            if rec.type_name == "SEAM_CURVE" {
                return true;
            }
        }
        let assoc = self
            .surface_curve_pcurves
            .borrow()
            .get(&curve_ref)
            .cloned()
            .unwrap_or_default();
        if assoc.len() != 2 {
            return false;
        }
        if !self.pcurve_lies_on(assoc[0], surf_ref) || !self.pcurve_lies_on(assoc[1], surf_ref) {
            return false;
        }
        // `cxx:104-116`: two oriented edges of the same wire share this edge.
        nb_oe == 2
    }

    /// `StepToTopoDS_GeometricTool::IsLikeSeam` (`cxx:133-218`): the two pcurves
    /// lie on the same surface but the edge is used once by the loop, and both
    /// are `LINE`s sharing an origin coordinate and a direction (CATIA BRep).
    pub(super) fn is_like_seam(&self, curve_ref: usize, surf_ref: usize, nb_oe: usize) -> bool {
        let assoc = self
            .surface_curve_pcurves
            .borrow()
            .get(&curve_ref)
            .cloned()
            .unwrap_or_default();
        if assoc.len() != 2 {
            return false;
        }
        if !self.pcurve_lies_on(assoc[0], surf_ref) || !self.pcurve_lies_on(assoc[1], surf_ref) {
            return false;
        }
        // `cxx:160-165`: the two oriented edges are not in the same wire.
        if nb_oe != 1 {
            return false;
        }
        let (Some(l1), Some(l2)) = (self.step_2d_line(assoc[0]), self.step_2d_line(assoc[1]))
        else {
            return false;
        };
        // `cxx:183-196`, `preci2d = Precision::PConfusion()`.
        let preci2d = occt_core::precision::PCONFUSION;
        let delta_x = (l1.0 - l2.0).abs();
        let delta_y = (l1.1 - l2.1).abs();
        let delta_dir_x = (l1.2 - l2.2).abs();
        let delta_dir_y = (l1.3 - l2.3).abs();
        if delta_x < preci2d || delta_y < preci2d {
            delta_dir_x < preci2d && delta_dir_y < preci2d
        } else {
            false
        }
    }

    /// `TranslateEdgeLoop.cxx:699-734`: order the two seam pcurves so the first
    /// entry is the FORWARD one. `ShapeAnalysis_Curve::SelectForwardSeam`, then
    /// flip when the cumulative edge/wire/face orientation is reversed.
    pub(super) fn order_seam_pcurves(
        &self,
        edge: &Edge,
        matched: &mut [Arc<dyn Curve2d>],
        face_ori: Orientation,
        wire_o: Orientation,
    ) {
        if matched.len() != 2 {
            return;
        }
        // WireO is FaceBound->Orientation (`cxx:272`), not the post-sameSense wire.
        let mut fwd =
            crate::pcurve_full::select_forward_seam(matched[0].as_ref(), matched[1].as_ref());
        if fwd != 0 {
            let edge_o = edge.0.orientation();
            let cumul = Orientation::compose(edge_o, wire_o);
            let cumul = Orientation::compose(cumul, face_ori);
            if cumul != Orientation::Forward {
                fwd = 3 - fwd;
            }
            if fwd == 2 {
                matched.swap(0, 1);
            }
        }
    }

    /// Match a wire edge's SURFACE_CURVE pcurve to `surf_ref` and attach its 2D
    /// curve to the edge for `face_key`. Source:
    /// `StepToTopoDS_GeometricTool::PCurve` + `StepToTopoDS_TranslateEdge::MakePCurve`.
    ///
    /// `nb_oe` is the number of oriented edges of the current loop that
    /// reference this edge; it selects the seam arm
    /// (`TranslateEdgeLoop.cxx:611-667`, `:689-767`, `:779-796`).
    pub(in crate::step) fn associate_edge_pcurve(
        &self,
        edge: &Edge,
        surf_ref: usize,
        face_key: usize,
        face_ori: Orientation,
        wire_o: Orientation,
        nb_oe: usize,
    ) -> Result<(), String> {
        let key = Arc::as_ptr(&edge.0.tshape) as usize;
        let Some(&curve_ref) = self.edge_curve_ref.borrow().get(&key) else {
            return Ok(());
        };
        let Some(pcurves) = self.surface_curve_pcurves.borrow().get(&curve_ref).cloned() else {
            return Ok(());
        };
        // Each associated_geometry entry is a pcurve_or_surface: a PCURVE (with
        // a 2D curve) or a SURFACE (intersection curve). Only PCURVEs carry the
        // pcurve; a SURFACE entry fails to resolve and is skipped. The raw 2D
        // curve is stored with its own parameterization; the meshing maps the
        // edge's 3D parameter onto it (`BRepMesh_EdgeParameterProvider`).
        let mut matched: Vec<Arc<dyn Curve2d>> = Vec::new();
        for pc in pcurves {
            if let Ok((basis_surf, c2d)) = self.resolve_pcurve(pc) {
                if basis_surf == surf_ref {
                    matched.push(c2d);
                }
            }
        }
        let is_seam = self.is_seam_curve(curve_ref, surf_ref, nb_oe);
        let is_like_seam = !is_seam && self.is_like_seam(curve_ref, surf_ref, nb_oe);
        let matched_n = matched.len();
        let stored: Vec<Arc<dyn Curve2d>> = if is_seam {
            // `cxx:736-767`: a seam edge carries two pcurves on the same face
            // (one per side of the seam, e.g. the cone's `u = 0` and `u = 2pi`),
            // stored forward-then-reversed so the wire assembly can hand each
            // traversal its own side.
            self.order_seam_pcurves(edge, &mut matched, face_ori, wire_o);
            matched
        } else if is_like_seam {
            // `cxx:748-767`: `else UpdateEdge(E, C2d2, Face, 0.)` - only the
            // forward pcurve is stored for CATIA-like seams.
            self.order_seam_pcurves(edge, &mut matched, face_ori, wire_o);
            matched.into_iter().take(1).collect()
        } else {
            // `cxx:645-667` walks every matching pcurve and keeps the last, then
            // `cxx:786-793` does `B.UpdateEdge(E, C2d, Face, 0.)` with it.
            matched.pop().into_iter().collect()
        };
        let stored_count = stored.len();
        if !stored.is_empty() {
            GeometryRegistry::global().set_edge_pcurves(&edge.0, face_key, stored);
        }
        Ok(())
    }

}
