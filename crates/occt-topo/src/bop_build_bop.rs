//! `BOPAlgo_Builder::BuildBOP` — select faces by IN/OUT of the opposite group
//! and rebuild solids with `BuilderSolid`.
//!
//! Source: `BOPAlgo_Builder.cxx:479`. Used by `BOPAlgo_BOP::BuildShape` when
//! argument solids are open (or BuilderSolid left unused faces). Closed-solid
//! Fuse/Cut/Common still go through `BuildRC` / `BuildSolid`.
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use std::collections::{HashMap, HashSet};

use crate::abs::Orientation;
use crate::algo_tools_face::is_split_to_reverse;
use crate::bop_builder2::BopBuilder;
use crate::bop_occt_util::{explore, shape_key};
use crate::builder::TopoBuilder;
use crate::builder_solid::BuilderSolid;
use crate::fclass2d::FaceState;
use crate::int_tools_full::IntToolsContext;
use crate::shape::TopoShape;
use crate::topo_tools_full::faces_of;

fn ori_key(s: &TopoShape) -> (usize, u8) {
    let o = match s.orientation() {
        Orientation::Forward => 0,
        Orientation::Reversed => 1,
        Orientation::Internal => 2,
        Orientation::External => 3,
    };
    (shape_key(s), o)
}

fn reverse_copy(s: &TopoShape) -> TopoShape {
    let mut r = s.clone();
    r.reverse();
    r
}

fn collect_group_faces(
    b: &BopBuilder,
    group: &[TopoShape],
    ctx: &mut IntToolsContext,
) -> (Vec<TopoShape>, HashSet<usize>, HashSet<usize>) {
    let mut faces_ori: Vec<TopoShape> = Vec::new();
    let mut faces: HashSet<usize> = HashSet::new();
    let mut in_keys: HashSet<usize> = HashSet::new();
    for shape in group {
        for solid in explore(shape, crate::abs::ShapeType::Solid) {
            for f in faces_of(&solid) {
                if f.0.orientation() != Orientation::Forward
                    && f.0.orientation() != Orientation::Reversed
                {
                    continue;
                }
                let imgs = match b.history().image(&f.0) {
                    Some(list) if !list.is_empty() => list.to_vec(),
                    _ => vec![f.0.clone()],
                };
                for fim in imgs {
                    let mut oriented = fim.clone();
                    if is_split_to_reverse(&oriented, &f.0, ctx).unwrap_or(false) {
                        oriented.reverse();
                    }
                    faces_ori.push(oriented.clone());
                    faces.insert(shape_key(&oriented));
                }
            }
            if let Some(ins) = b.in_parts().get(&shape_key(&solid)) {
                for f in ins {
                    in_keys.insert(shape_key(f));
                }
            }
        }
    }
    (faces_ori, faces, in_keys)
}

/// `BOPAlgo_Builder::BuildBOP`.
pub fn build_bop(b: &mut BopBuilder) -> Result<(), String> {
    let obj_state = b.obj_state();
    let tools_state = b.tools_state();
    if !matches!(obj_state, FaceState::In | FaceState::Out)
        || !matches!(tools_state, FaceState::In | FaceState::Out)
    {
        return Err("BOPAlgo_AlertBOPNotSet".into());
    }
    let mut ctx = IntToolsContext::new();
    let objects = b.objects().to_vec();
    let tools = b.tools().to_vec();
    let (obj_ori, obj_faces, in_obj) = collect_group_faces(b, &objects, &mut ctx);
    let (tool_ori, tool_faces, in_tool) = collect_group_faces(b, &tools, &mut ctx);

    let is_objects_in = obj_state == FaceState::In;
    let is_tools_in = tools_state == FaceState::In;
    let avoid_in = !is_objects_in && !is_tools_in;
    let avoid_in_both = is_objects_in != is_tools_in;
    let same_ori_needed = obj_state == tools_state;

    let mut res_ori: Vec<TopoShape> = Vec::new();
    let mut res_fence: HashSet<usize> = HashSet::new();
    let mut avoid: HashSet<usize> = HashSet::new();
    let mut fence: HashSet<usize> = HashSet::new();
    let mut fence_ori: HashSet<(usize, u8)> = HashSet::new();

    for (i, map) in [obj_ori.as_slice(), tool_ori.as_slice()].into_iter().enumerate() {
        let opposite = if i == 0 { &tool_faces } else { &obj_faces };
        let in_map = if i == 0 { &in_obj } else { &in_tool };
        let in_opp = if i == 0 { &in_tool } else { &in_obj };
        let take_in = if i == 0 { is_objects_in } else { is_tools_in };
        for fim in map {
            let k = shape_key(fim);
            let is_in = in_map.contains(&k);
            let is_in_opp = in_opp.contains(&k);
            if avoid_in && (is_in || is_in_opp) {
                continue;
            }
            if avoid_in_both && is_in && is_in_opp {
                continue;
            }
            if !fence.insert(k) {
                if !opposite.contains(&k) {
                    if take_in != same_ori_needed {
                        avoid.insert(k);
                    }
                } else {
                    let same_ori = !fence_ori.insert(ori_key(fim));
                    if same_ori_needed == same_ori {
                        if res_fence.insert(k) {
                            res_ori.push(fim.clone());
                        }
                    } else {
                        avoid.insert(k);
                    }
                    continue;
                }
            }
            if !fence_ori.insert(ori_key(fim)) {
                continue;
            }
            if take_in == is_in_opp {
                if is_in {
                    res_ori.push(fim.clone());
                    res_ori.push(reverse_copy(fim));
                } else if take_in && !same_ori_needed {
                    res_ori.push(reverse_copy(fim));
                } else {
                    res_ori.push(fim.clone());
                }
                res_fence.insert(k);
            }
        }
    }

    let res_faces: Vec<TopoShape> = res_ori
        .into_iter()
        .filter(|f| !avoid.contains(&shape_key(f)))
        .collect();
    if res_faces.is_empty() {
        return Err("BOPAlgo_AlertBuilderFailed".into());
    }
    let mut bs = BuilderSolid::new();
    bs.set_shapes(res_faces.clone());
    bs.set_fuzzy(b.fuzzy_value());
    bs.perform()?;
    let obj_ori_keys: HashSet<(usize, u8)> = {
        let mut s = HashSet::new();
        for f in &obj_ori {
            s.insert(ori_key(f));
        }
        for f in &tool_ori {
            s.insert(ori_key(f));
        }
        s
    };
    let mut solids: Vec<TopoShape> = Vec::new();
    for area in bs.areas() {
        let ok = faces_of(area).iter().any(|f| obj_ori_keys.contains(&ori_key(&f.0)));
        if ok {
            solids.push(area.clone());
        }
    }
    if solids.is_empty() {
        return Err("BOPAlgo_AlertBuilderFailed".into());
    }
    let bb = TopoBuilder::new();
    b.set_result_shape(bb.make_compound_of(&solids).0);
    Ok(())
}

/// Map of IN faces used by BuildBOP (`myInParts`).
pub fn in_face_keys(in_parts: &HashMap<usize, Vec<TopoShape>>) -> HashSet<usize> {
    let mut s = HashSet::new();
    for faces in in_parts.values() {
        for f in faces {
            s.insert(shape_key(f));
        }
    }
    s
}
