use super::prelude::*;
use super::*;

impl AlgoTools {

    /// `CorrectTolerances` with an avoid-map of source V/E/F keys.
    pub fn correct_tolerances_avoid(shape: &TopoShape, avoid: &HashSet<usize>, tol: f64) {
        for edge in edges_of(shape) {
            if avoid.contains(&GeometryRegistry::shape_key(&edge.0)) {
                continue;
            }
            let e_tol = BRepTool::edge_tolerance(&edge);
            let curve = BRepTool::edge_curve(&edge);
            let (a, b) = BRepTool::edge_parameters(&edge);
            let kids: Vec<TopoShape> = edge
                .0
                .tshape
                .read()
                .unwrap()
                .children
                .iter()
                .filter(|h| h.shape_type() == ShapeType::Vertex)
                .cloned()
                .collect();
            if let (Some(curve), true, true) = (curve, a.is_finite(), b.is_finite()) {
                let pa = curve.d0(a);
                let pb = curve.d0(b);
                for (i, p_end) in [(0usize, &pa), (1, &pb)] {
                    if let Some(k) = kids.get(i) {
                        if avoid.contains(&GeometryRegistry::shape_key(k)) {
                            continue;
                        }
                        let v = Vertex(k.clone());
                        let want = BRepTool::vertex_point(&v).distance(p_end) + D_TOLERANCE;
                        if want > BRepTool::vertex_tolerance(&v) && want <= tol {
                            v.set_tolerance(want);
                        }
                    }
                }
            }
            for k in &kids {
                if avoid.contains(&GeometryRegistry::shape_key(k)) {
                    continue;
                }
                let v = Vertex(k.clone());
                let cur = BRepTool::vertex_tolerance(&v);
                if e_tol > cur && e_tol <= tol {
                    v.set_tolerance(e_tol);
                }
            }
        }
    }

    /// Find a parameter of `edge` at which the edge's 3-D curve is *off* `face`
    /// — a sampled point whose UV classification on `face` is `Out`
    /// (`BOPTools_AlgoTools::GetEdgeOff`, geometric form).
    ///
    /// Returns `None` when the whole edge lies on or inside the face.
    pub fn get_edge_off(edge: &Edge, face: &Face) -> Option<f64> {
        let surf = BRepTool::face_surface(face)?;
        let cl = FClass2d::new(face, 1e-7).ok()?;
        let (a, b) = BRepTool::edge_parameters(edge);
        if !a.is_finite() || !b.is_finite() {
            return None;
        }
        for i in 0..=32 {
            let t = a + (b - a) * i as f64 / 32.0;
            let p = AlgoTools::point_on_edge(edge, t).ok()?;
            // Faithful `GeomAPI_ProjectPointOnSurf` (`Extrema_ExtPS`, T-52); the
            // previous 32x32 grid was A1's substitute. A sample whose projection
            // is not done contributes nothing (OCCT's classifier has no (u, v)
            // to work with there).
            let Some(ps) = occt_geom::geom_api::project_point_on_surface(surf.as_ref(), &p, 1e-7)
            else {
                continue;
            };
            if cl.perform(GpPnt2d::new(ps.u, ps.v)) == FaceState::Out {
                return Some(t);
            }
        }
        None
    }

    // ---- private helpers ---------------------------------------------------

    /// Distinct `TShape` identity keys of the edges of `shape`. A face yields
    /// the edges of its boundary wires; an edge yields itself; any other shape
    /// yields its distinct edge sub-shapes.
    pub(super) fn shape_edge_keys(shape: &TopoShape) -> Vec<usize> {
        match shape.shape_type() {
            ShapeType::Edge => vec![GeometryRegistry::shape_key(shape)],
            ShapeType::Face => {
                let f = Face(shape.clone());
                let mut out = Vec::new();
                for w in wires_of_face(&f) {
                    for e in edges_of_wire(&w) {
                        out.push(GeometryRegistry::shape_key(&e.0));
                    }
                }
                out
            }
            _ => edges_of(shape).into_iter().map(|e| GeometryRegistry::shape_key(&e.0)).collect(),
        }
    }

    /// Wrap a connected group into a [`ConnexityBlock`] and compute its
    /// regularity (every shared edge used by exactly two shapes).
    pub(super) fn connexity_block(group: Vec<TopoShape>) -> ConnexityBlock {
        let mut block = ConnexityBlock::new();
        block.change_shapes_mut().extend(group);
        let edge_sets: Vec<Vec<usize>> =
            block.shapes().iter().map(AlgoTools::shape_edge_keys).collect();
        let mut counts: HashMap<usize, usize> = HashMap::new();
        for es in &edge_sets {
            for e in es {
                *counts.entry(*e).or_insert(0) += 1;
            }
        }
        block.set_regular(counts.values().all(|&c| c == 2));
        block
    }

    /// All boundary edges of a face (its wires' edges).
    pub(super) fn face_edges(face: &Face) -> Vec<Edge> {
        let mut out = Vec::new();
        for w in wires_of_face(face) {
            for e in edges_of_wire(&w) {
                out.push(e);
            }
        }
        out
    }

    /// Unit surface normal of `face` at the parameter point closest to `p`.
    pub(super) fn face_normal_at_point(face: &Face, p: &GpPnt) -> Option<GpVec> {
        let surf = BRepTool::face_surface(face)?;
        // Faithful `GeomAPI_ProjectPointOnSurf` (`Extrema_ExtPS`, T-52).
        let (u, v) = occt_geom::geom_api::project_point_on_surface(surf.as_ref(), p, 1e-7)
            .map(|ps| (ps.u, ps.v))?;
        let n = surface_normal(surf.as_ref(), u, v);
        if n.square_magnitude() < 1e-30 { None } else { Some(n) }
    }
}
/// 2D (p-curve) helpers mirroring `BOPTools_AlgoTools2D`.
pub struct AlgoTools2D;

impl AlgoTools2D {
    /// The p-curve of `edge` on `face` — a convenience alias of
    /// [`AlgoTools::make_pcurve`] matching the OCCT `EdgeToFace` role of
    /// producing the edge's UV representation on a face.
    pub fn edge_to_face(edge: &Edge, face: &Face) -> Result<Arc<dyn Curve2d>, String> {
        AlgoTools::make_pcurve(edge, face)
    }

    /// `BOPTools_AlgoTools2D::AdjustPCurveOnSurf`
    /// (`BOPTools_AlgoTools2D.cxx:247-400`): shift `curve` by whole periods of
    /// the surface so that the point at the **middle of the edge's range**
    /// `[first, last]` (OCCT's `aT = 0.5*(aFirst+aLast)`, *not* the pcurve's own
    /// range — a `Geom2d_Line` may declare an infinite one) lies inside the
    /// face's UV bounds.
    ///
    /// OCCT never trims the pcurve here: the result either is the same curve or
    /// the same curve translated by `(du, dv)` (`cxx:389-399`). The pieces are,
    /// in order:
    /// * a whole-period shift of `u2`/`v2` onto the bounds (`cxx:273-344`),
    ///   including the cylinder special case `dFi = MaxToleranceEdge(face)/R`
    ///   and the `(VMax - VMin) < aVPeriod` tie-break for `dv`;
    /// * a `BRepClass_FaceClassifier` cross-check of `(u2 + du, v2 + dv)` for
    ///   surfaces whose period is narrower than the face range
    ///   (`cxx:346-387`), which may flip `du`/`dv` by one more period;
    /// * the translation itself (`cxx:388-399`).
    ///
    /// The classifier here is the port's [`crate::fclass2d::FClass2d`]; when it
    /// cannot be built (a degenerate face) the cross-check is skipped, which
    /// leaves `du`/`dv` as the first stage computed them.
    pub fn adjust_pcurve_on_surf(
        curve: &Arc<dyn Curve2d>,
        face: &Face,
        first: f64,
        last: f64,
        _tol: f64,
    ) -> Result<Arc<dyn Curve2d>, String> {
        use crate::fclass2d::{FaceState, FClass2d};

        let surf = crate::brep_tool::BRepTool::face_surface(face)
            .ok_or("AdjustPCurveOnSurf: face has no surface")?;
        let (umin, umax, vmin, vmax) = crate::brep_tool::BRepTool::uv_bounds(face);
        let u_periodic = surf.is_u_periodic();
        let v_periodic = surf.is_v_periodic();
        let u_period = if u_periodic { surf.u_period() } else { 0.0 };
        let v_period = if v_periodic { surf.v_period() } else { 0.0 };
        let a_delta = occt_core::precision::PCONFUSION;

        let a_t = 0.5 * (first + last);
        let p = curve.d0(a_t);
        let (mut u2, mut v2) = (p.x(), p.y());

        // du (`cxx:273-315`)
        let mut du = 0.0;
        if u_periodic && u_period > 0.0 {
            if (u2 - umin).abs() < a_delta {
                u2 = umin;
            } else if (u2 - umin - u_period).abs() < a_delta {
                u2 = umin + u_period;
            }
            let (nu2, ndu) =
                crate::geom_int::adjust_periodic(u2, umin, umax, u_period, 0.0);
            u2 = nu2;
            du = ndu;
            if du == 0.0 {
                if let Some(radius) = cylinder_radius(surf.as_ref()) {
                    let a_tol = max_tolerance_edge(face);
                    let mut d_fi = a_tol / radius;
                    if d_fi < a_delta {
                        d_fi = a_delta;
                    }
                    let min_cond = umin - u2 > d_fi;
                    let max_cond = u2 - umax > d_fi;
                    if min_cond || max_cond {
                        du = if min_cond { u_period } else { -u_period };
                    }
                }
            }
        }

        // dv (`cxx:317-344`)
        let mut dv = 0.0;
        if v_periodic && v_period > 0.0 {
            let min_cond = vmin - v2 > a_delta;
            let max_cond = v2 - vmax > a_delta;
            if min_cond || max_cond {
                dv = if min_cond { v_period } else { -v_period };
            }
            if (vmax - vmin < v_period) && dv != 0.0 {
                let v_mid = 0.5 * (vmin + vmax);
                let d_vm = (v2 - v_mid).abs();
                let d_vr = (v2 + dv - v_mid).abs();
                if d_vm < d_vr {
                    dv = 0.0;
                }
            }
        }

        // Classifier cross-check (`cxx:346-387`)
        if let Ok(cl) = FClass2d::new(face, a_delta) {
            if u_periodic && u_period > 0.0 && (umax - umin - 2.0 * a_delta) > u_period {
                let u = u2 + du;
                if u > umin + a_delta + u_period || u < umax - a_delta - u_period {
                    if cl.perform(GpPnt2d::new(u, v2 + dv)) == FaceState::Out {
                        du += if u > umin + a_delta + u_period {
                            -u_period
                        } else {
                            u_period
                        };
                    }
                }
            }
            if v_periodic && v_period > 0.0 && (vmax - vmin - 2.0 * a_delta) > v_period {
                let u = u2 + du;
                let v = v2 + dv;
                if v > vmin + a_delta + v_period || v < vmax - a_delta - v_period {
                    if cl.perform(GpPnt2d::new(u, v)) == FaceState::Out {
                        dv += if v > vmin + a_delta + v_period {
                            -v_period
                        } else {
                            v_period
                        };
                    }
                }
            }
        }

        if du == 0.0 && dv == 0.0 {
            return Ok(curve.clone());
        }
        let mut tr = occt_core::gp::GpTrsf2d::identity();
        tr.set_translation_vec(&occt_core::gp::GpVec2d::new(du, dv));
        Ok(Arc::from(curve.transformed(&tr)))
    }
}

/// `BRep_Tool::MaxTolerance(face, TopAbs_EDGE)`.
fn max_tolerance_edge(face: &Face) -> f64 {
    crate::topo_tools_full::edges_of(&face.0)
        .iter()
        .map(|e| crate::brep_tool::BRepTool::edge_tolerance(e))
        .fold(0.0f64, f64::max)
}

/// The radius of a cylindrical surface, for the `dFi = aTol / aR` branch of
/// `AdjustPCurveOnSurf` (`cxx:294-313`).
fn cylinder_radius(surf: &dyn Surface) -> Option<f64> {
    surf.gp_cylinder().map(|c| c.radius())
}

/// Shape-set utilities mirroring `BOPTools_Set`.
pub struct BOPToolsSet;

impl BOPToolsSet {
    /// Deduplicate `shapes`, keeping the first occurrence of each distinct
    /// shape. Two shapes are the same when they share the same `TShape` data
    /// (`TopoShape::same_tshape` — pointer identity), matching OCCT's shape-map
    /// hashing by `TopoDS_Shape`.
    pub fn shape_list(shapes: &[TopoShape]) -> Vec<TopoShape> {
        let mut out: Vec<TopoShape> = Vec::new();
        for s in shapes {
            if !out.iter().any(|o| o.same_tshape(s)) {
                out.push(s.clone());
            }
        }
        out
    }

    /// Count how many shapes in `shapes` have type `kind` (`TopAbs_ShapeEnum`).
    pub fn type_count(shapes: &[TopoShape], kind: ShapeType) -> usize {
        shapes.iter().filter(|s| s.shape_type() == kind).count()
    }
}
