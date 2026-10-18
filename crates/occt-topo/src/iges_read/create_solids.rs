//! `CreateSolids` / `CollectSolids` (`ShapeFix_Solid.cxx:115-456`, `:577-640`).
//! `CreateOpenSolidMode` is 1 (`IGESData.cxx:119`). The one-shell arm stays in
//! `perform_iges_solid`.

use occt_core::precision::CONFUSION;

use crate::abs::ShapeType;
use crate::brep_class3d::SolidClassifier;
use crate::builder::TopoBuilder;
use crate::fclass2d::FaceState;
use crate::shape::{Shell, TopoShape};
use crate::shape_fix_compose_shell::{ReShape, SharedReShape};
use crate::topo_tools_full::{edges_of, faces_of, is_same, vertices_of, vertex_position};

struct ShellHoles {
    shell: TopoShape,
    holes: Vec<TopoShape>,
}

fn contains_shape(list: &[TopoShape], shape: &TopoShape) -> bool {
    list.iter().any(|item| is_same(item, shape))
}

fn find_holes(map: &[ShellHoles], shell: &TopoShape) -> Option<usize> {
    map.iter().position(|item| is_same(&item.shell, shell))
}

fn classify_infinite(shell: &TopoShape) -> FaceState {
    let mut trial = shell.clone();
    if !trial.free() {
        trial.set_free(true);
    }
    let solid = TopoBuilder::new().make_solid(&[Shell(trial)]);
    let mut sc = SolidClassifier::new();
    sc.load(solid.0);
    sc.perform_infinite_point(CONFUSION);
    sc.state()
}

fn status_of(state: FaceState) -> i32 {
    match state {
        FaceState::In => 1,
        FaceState::Out => 2,
        _ => 0,
    }
}

fn state_of(status: i32) -> FaceState {
    if status == 1 {
        FaceState::In
    } else if status == 2 {
        FaceState::Out
    } else {
        FaceState::Unknown
    }
}

/// `CollectSolids` (`cxx:115-302`). The status map is `theMapStatus`.
fn collect_solids(shells: &[TopoShape]) -> (Vec<ShellHoles>, Vec<(TopoShape, i32)>) {
    let mut map: Vec<ShellHoles> = shells
        .iter()
        .map(|shell| ShellHoles {
            shell: shell.clone(),
            holes: Vec::new(),
        })
        .collect();
    let mut holes = Vec::new();
    let mut status: Vec<(TopoShape, i32)> = Vec::new();
    for (i, shell1) in shells.iter().enumerate() {
        if !shell1.is_closed() || edges_of(shell1).is_empty() {
            continue;
        }
        let st = if let Some((_, found)) = status.iter().find(|(s, _)| is_same(s, shell1)) {
            *found
        } else {
            let st = status_of(classify_infinite(shell1));
            status.push((shell1.clone(), st));
            st
        };
        if st == 0 {
            continue;
        }
        let infin = state_of(st);
        let solid = {
            let mut trial = shell1.clone();
            if !trial.free() {
                trial.set_free(true);
            }
            TopoBuilder::new().make_solid(&[Shell(trial)]).0
        };
        let mut sc = SolidClassifier::new();
        sc.load(solid);
        for (j, shell2) in shells.iter().enumerate() {
            if i == j || !shell2.is_closed() || contains_shape(&holes, shell2) {
                continue;
            }
            if let Some(idx) = find_holes(&map, shell2) {
                if contains_shape(&map[idx].holes, shell1) {
                    continue;
                }
            }
            let mut point = FaceState::Unknown;
            let mut numon = 0;
            for vertex in vertices_of(shell2).into_iter().take(10) {
                if point != FaceState::Unknown && !(point == FaceState::On && numon < 3) {
                    break;
                }
                sc.perform(&vertex_position(&vertex), CONFUSION);
                point = sc.state();
                if point == FaceState::On {
                    numon += 1;
                }
            }
            if numon == 3 && point == FaceState::On {
                point = FaceState::Out;
            }
            if point != infin {
                if let Some(idx) = find_holes(&map, shell1) {
                    map[idx].holes.push(shell2.clone());
                }
                if let Some(pos) = holes.iter().position(|item| is_same(item, shell2)) {
                    holes.remove(pos);
                } else {
                    holes.push(shell2.clone());
                }
            }
        }
    }
    let mut nested = Vec::new();
    for item in &map {
        if contains_shape(&holes, &item.shell) || item.holes.is_empty() {
            continue;
        }
        for hole in &item.holes {
            if !contains_shape(&holes, hole) {
                continue;
            }
            if let Some(idx) = find_holes(&map, hole) {
                nested.extend(map[idx].holes.clone());
            }
        }
    }
    holes.retain(|hole| !contains_shape(&nested, hole));
    map.retain(|item| !contains_shape(&holes, &item.shell));
    let mut ordered = Vec::new();
    for shell in shells {
        if let Some(idx) = find_holes(&map, shell) {
            ordered.push(ShellHoles {
                shell: map[idx].shell.clone(),
                holes: map[idx].holes.clone(),
            });
        }
    }
    (ordered, status)
}

fn cached_state(status: &[(TopoShape, i32)], shell: &TopoShape) -> Option<FaceState> {
    status
        .iter()
        .find(|(item, _)| is_same(item, shell))
        .map(|(_, st)| state_of(*st))
}

/// Hole status (`cxx:383-386`): a bound `0` is `OUT`, unlike the outer shell.
fn hole_state(status: &[(TopoShape, i32)], shell: &TopoShape) -> Option<FaceState> {
    status.iter().find(|(item, _)| is_same(item, shell)).map(|(_, st)| {
        if *st == 1 {
            FaceState::In
        } else {
            FaceState::Out
        }
    })
}

fn explorer_solids(shape: &TopoShape) -> Vec<TopoShape> {
    match shape.shape_type() {
        ShapeType::Solid => vec![shape.clone()],
        ShapeType::CompSolid | ShapeType::Compound => {
            let children = shape.tshape.read().expect("lock").children.clone();
            children.iter().flat_map(explorer_solids).collect()
        }
        _ => Vec::new(),
    }
}

/// `CreateSolids` (`cxx:306-456`). `None` is `isDone == false`.
fn create_solids(shape: &TopoShape) -> Option<Vec<TopoShape>> {
    let shells = {
        let mut out = Vec::new();
        fn walk(shape: &TopoShape, out: &mut Vec<TopoShape>) {
            if shape.shape_type() == ShapeType::Shell {
                out.push(shape.clone());
                return;
            }
            let children = shape.tshape.read().expect("lock").children.clone();
            for child in children {
                walk(&child, out);
            }
        }
        walk(shape, &mut out);
        out
    };
    let (roots, status) = collect_solids(&shells);
    let mut shell_solid: Vec<(TopoShape, TopoShape)> = Vec::new();
    let mut done = false;
    for root in &roots {
        let shell = &root.shell;
        if !shell.is_closed() || edges_of(shell).is_empty() {
            shell_solid.push((shell.clone(), shell.clone()));
            done = true;
            continue;
        }
        let infin = cached_state(&status, shell).unwrap_or_else(|| classify_infinite(shell));
        let mut body = shell.clone();
        if infin == FaceState::In {
            done = true;
            if !body.free() {
                body.set_free(true);
            }
            body.reverse();
        }
        if !body.free() {
            body.set_free(true);
        }
        let mut solid = TopoBuilder::new().make_solid(&[Shell(body)]);
        for hole in &root.holes {
            let hole_state = hole_state(&status, hole).unwrap_or_else(|| classify_infinite(hole));
            let mut piece = hole.clone();
            if hole_state == FaceState::Out {
                if !piece.free() {
                    piece.set_free(true);
                }
                piece.reverse();
                done = true;
            }
            if !piece.free() {
                piece.set_free(true);
            }
            TopoBuilder::new().add_shell(&mut solid, &Shell(piece));
        }
        shell_solid.push((shell.clone(), solid.0));
    }
    let mut face_shells: Vec<(TopoShape, Vec<TopoShape>)> = Vec::new();
    for shell in &shells {
        for face in faces_of(shell) {
            if let Some(slot) = face_shells.iter_mut().find(|(key, _)| is_same(key, &face.0)) {
                if !contains_shape(&slot.1, shell) {
                    slot.1.push(shell.clone());
                }
            } else {
                face_shells.push((face.0, vec![shell.clone()]));
            }
        }
    }
    for (_, owners) in &face_shells {
        if owners.len() < 2 {
            continue;
        }
        done = shape.shape_type() != ShapeType::CompSolid || done;
        let mut solids = Vec::new();
        for owner in owners {
            if let Some((_, value)) = shell_solid.iter().find(|(key, _)| is_same(key, owner)) {
                solids.extend(explorer_solids(value));
            }
        }
        if solids.len() > 1 {
            let mut comp = TopoShape::new(ShapeType::CompSolid);
            comp.set_free(true);
            let builder = TopoBuilder::new();
            for solid in &solids {
                builder.add(&mut comp, solid);
            }
            for owner in owners {
                if let Some(slot) = shell_solid.iter_mut().find(|(key, _)| is_same(key, owner)) {
                    slot.1 = comp.clone();
                }
            }
        }
    }
    let mut out = Vec::new();
    for (_, value) in &shell_solid {
        if !contains_shape(&out, value) {
            out.push(value.clone());
        }
    }
    done = out.len() > 1 || done;
    if done { Some(out) } else { None }
}

/// `ShapeFix_Solid::Perform` multi-shell arm (`cxx:577-640`).
pub(in crate::iges_read) fn replace_multi_shell(shape: &TopoShape, ctx: &mut SharedReShape) {
    let Some(solids) = create_solids(shape) else {
        return;
    };
    if solids.len() == 1 {
        let result = if shape.shape_type() == ShapeType::Shell {
            let mut piece = solids[0].clone();
            if !piece.free() {
                piece.set_free(true);
            }
            TopoBuilder::new().make_solid(&[Shell(piece)]).0
        } else {
            solids[0].clone()
        };
        ctx.replace(shape, &result);
    } else if solids.len() > 1 {
        let mut compound = TopoShape::new(ShapeType::Compound);
        compound.set_free(true);
        let builder = TopoBuilder::new();
        for item in &solids {
            let piece = if shape.shape_type() == ShapeType::Shell {
                let mut open = item.clone();
                open.set_closed(false);
                if !open.free() {
                    open.set_free(true);
                }
                TopoBuilder::new().make_solid(&[Shell(open)]).0
            } else {
                item.clone()
            };
            builder.add(&mut compound, &piece);
        }
        ctx.replace(shape, &compound);
    }
}
