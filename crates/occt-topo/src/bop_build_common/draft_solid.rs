use super::prelude::*;
use super::*;

/// 3-D state of `shape` relative to `solid`, from a single representative point
/// (`BOPTools_AlgoTools::ComputeStateByOnePoint`, `BOPTools_AlgoTools.cxx:623-657`).
///
/// The per-shape-type control flow is shared with the public
/// [`crate::bopalgo_tools_class::compute_state_by_one_point`]; this wrapper only
/// adapts it to the `Result` shape the draft-solid builder uses.
pub(super) fn compute_state_by_one_point(
    shape: &TopoShape,
    solid: &TopoShape,
    tol: f64,
) -> Result<FaceState, String> {
    Ok(crate::bopalgo_tools_class::compute_state_by_one_point(
        shape, solid, tol,
    ))
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
