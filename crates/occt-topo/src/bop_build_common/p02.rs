use super::prelude::*;
use super::*;

/// 3-D state of `shape` relative to `solid`, from a single representative point
/// (`BOPTools_AlgoTools::ComputeStateByOnePoint`).
pub(super) fn compute_state_by_one_point(
    shape: &TopoShape,
    solid: &TopoShape,
    tol: f64,
) -> Result<FaceState, String> {
    match shape.shape_type() {
        ShapeType::Vertex => {
            let p = BRepTool::vertex_point(&Vertex(shape.clone()));
            Ok(point_solid_state(solid, &p, tol))
        }
        ShapeType::Edge => {
            let e = Edge(shape.clone());
            let (a, b) = BRepTool::edge_parameters(&e);
            let Some(curve) = BRepTool::edge_curve(&e) else {
                return Ok(FaceState::Unknown);
            };
            let p = if a.is_finite() && b.is_finite() {
                curve.d0(0.5 * (a + b))
            } else {
                return Ok(FaceState::Unknown);
            };
            Ok(point_solid_state(solid, &p, tol))
        }
        ShapeType::Face => {
            // Prefer a boundary edge not lying on the solid; fall back to the
            // face surface centre when every edge is on the solid.
            let f = Face(shape.clone());
            let solid_edges = edges_of(solid);
            let mut p: Option<GpPnt> = None;
            for e in edges_of(&f.0) {
                if BRepTool::is_degenerated(&e) {
                    continue;
                }
                if solid_edges.iter().any(|se| se.same_tshape(&e.0)) {
                    continue;
                }
                let (a, b) = BRepTool::edge_parameters(&e);
                if a.is_finite() && b.is_finite() {
                    if let Some(c) = BRepTool::edge_curve(&e) {
                        p = Some(c.d0(0.5 * (a + b)));
                        break;
                    }
                }
            }
            let p = match p {
                Some(p) => p,
                None => match face_sample_point(&f) {
                    Some(p) => p,
                    None => return Ok(FaceState::Unknown),
                },
            };
            Ok(point_solid_state(solid, &p, tol))
        }
        _ => {
            let kids = shape.tshape.read().unwrap().children.clone();
            for sub in kids {
                let st = compute_state_by_one_point(&sub, solid, tol)?;
                if st != FaceState::Unknown {
                    return Ok(st);
                }
            }
            Ok(FaceState::Unknown)
        }
    }
}

/// State of a point relative to a solid: `On` within `tol` of a boundary face,
/// otherwise `In`/`Out` by the parity test. Distances are measured against the
/// solid's faces only, so internal vertices/edges already children of the solid
/// do not make the point `On` (matching the OCCT `BRepClass3d` classifier used
/// by `ComputeState`).
pub(super) fn point_solid_state(solid: &TopoShape, p: &GpPnt, tol: f64) -> FaceState {
    let mut d = f64::INFINITY;
    for f in faces_of(solid) {
        let (_, q) = closest_point_on_face(&f, p, 16, 16);
        d = d.min(q.distance(p));
    }
    if d.is_finite() && d <= tol {
        return FaceState::On;
    }
    if is_inside(solid, p) {
        FaceState::In
    } else {
        FaceState::Out
    }
}

// ---------------------------------------------------------------------------
// Draft solid
// ---------------------------------------------------------------------------

/// Builds a draft solid from a (closed) shell, rebuilding the shell from the
/// face splits.
///
/// Mirrors `BOPAlgo_Builder::BuildDraftSolid` restricted to a single shell
/// argument: each face of the shell is replaced by its image splits (a split
/// whose orientation is inverted relative to the original face is reversed
/// first); `INTERNAL` faces are dropped (the caller collects them separately in
/// OCCT — this port does not expose the internal-face list). The rebuilt shell
/// is flagged closed and wrapped into a solid.
///
/// A solid input is treated as a collection of shells, mirroring the OCCT loop.
pub fn build_draft_solid<B: BopBuildOps>(
    f: &mut B,
    shell: &TopoShape,
) -> Result<TopoShape, String> {
    let bld = TopoBuilder::new();
    let solid_or = shell.orientation();
    let mut solid = Solid::new();
    solid.0.set_orientation(solid_or);

    let shells: Vec<TopoShape> = match shell.shape_type() {
        ShapeType::Solid => direct_children(shell)
            .into_iter()
            .filter(|c| c.shape_type() == ShapeType::Shell)
            .collect(),
        _ => vec![shell.clone()],
    };

    for sh in &shells {
        let mut new_shell = Shell::new();
        new_shell.0.set_orientation(sh.orientation());
        let mut i_flag = false;
        for child in direct_children(sh) {
            let or = child.orientation();
            if let Some(images) = f.history().image(&child) {
                for im in images {
                    let mut fx = im.clone();
                    if has_same_domain(f.ds(), &fx) {
                        if or == Orientation::Internal {
                            // Internal face: collected by the caller, not added.
                        } else {
                            if is_split_to_reverse(&fx, &child) {
                                reverse_orientation(&mut fx);
                            }
                            if add_draft_face(&bld, &mut new_shell.0, &fx) {
                                i_flag = true;
                            }
                        }
                    } else {
                        fx.set_orientation(or);
                        if or == Orientation::Internal {
                            // Internal face.
                        } else {
                            if add_draft_face(&bld, &mut new_shell.0, &fx) {
                                i_flag = true;
                            }
                        }
                    }
                }
            } else if or != Orientation::Internal {
                if add_draft_face(&bld, &mut new_shell.0, &child) {
                    i_flag = true;
                }
            }
        }
        if i_flag {
            new_shell.set_closed(!AlgoTools::is_open_shell(&new_shell.0));
            bld.add_shell(&mut solid, &new_shell);
        }
    }
    Ok(solid.0)
}

/// True when `shape` has a same-domain counterpart in the data structure.
pub(super) fn has_same_domain(ds: &BopdsDS, shape: &TopoShape) -> bool {
    ds.index(shape)
        .and_then(|i| ds.has_shape_sd(i))
        .is_some()
}

/// Adds `face` to the rebuilt shell, once. A face that carries no area
/// (its boundary collapses to fewer than three distinct vertices) cannot bound
/// a shell and is dropped — the coincident-source-face splitting can produce
/// such degenerate sliver pieces alongside the real split face. A face already
/// present in the shell is not added twice. Returns true when the face was
/// added.
pub(super) fn add_draft_face(bld: &TopoBuilder, shell: &mut TopoShape, face: &TopoShape) -> bool {
    if face_is_degenerate(&Face(face.clone())) {
        return false;
    }
    if set_contains(&direct_children(shell), face) {
        return false;
    }
    bld.add(shell, face);
    true
}

/// True when the face has no area: fewer than three distinct boundary-vertex
/// positions (a sliver / segment wire produced by an over-split).
pub(super) fn face_is_degenerate(face: &Face) -> bool {
    let mut keys: Vec<(i64, i64, i64)> = Vec::new();
    for v in vertices_of(&face.0) {
        let p = BRepTool::vertex_point(&v);
        let k = (
            (p.x() / 1e-6).round() as i64,
            (p.y() / 1e-6).round() as i64,
            (p.z() / 1e-6).round() as i64,
        );
        if !keys.contains(&k) {
            keys.push(k);
        }
        if keys.len() >= 3 {
            return false;
        }
    }
    true
}
