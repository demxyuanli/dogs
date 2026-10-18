//! `BOPAlgo_Builder::BuildSplitSolids`.
//!
//! Source: `BOPAlgo_Builder_3.cxx:413-618`. For each source solid bound in
//! `theDraftSolids`:
//! * no IN faces → the draft itself is the image;
//! * IN faces present → draft faces plus each IN face FORWARD and REVERSED
//!   are fed to `BOPAlgo_SplitSolid` / [`crate::builder_solid::BuilderSolid`];
//! * the resulting areas are interned with [`crate::bop_tools_set::BopToolsSetMap`]
//!   so same-domain solids collapse to one representative (`myShapesSD`).
//!
//! This is the GF (General Fuse) split, with **no** operation-state filter.
//! BOP Cut/Common/Fuse filtering is `crate::bop_bop::build_shape`
//! (`BOPAlgo_BOP::BuildShape`), not `BuildResult`.

use std::collections::HashSet;

use crate::abs::ShapeType;
use crate::bop_fill_in3d::FillIn3dPartsResult;
use crate::bop_occt_util::{
    both_orientations, images_bound, nb_source, origins_append, shape_key, BopSolidHost,
};
use crate::bop_tools_set::{intern_solid_sd, BopToolsSet, BopToolsSetMap};
use crate::builder_solid::BuilderSolid;
use crate::shape::TopoShape;
use crate::topo_tools::{faces_of, ShapeMapKey};

/// `BOPAlgo_SplitSolid`: BuilderSolid plus the source solid it splits.
struct SplitSolid {
    solid: TopoShape,
    builder: BuilderSolid,
    warnings: Vec<String>,
}

impl SplitSolid {
    fn new(solid: TopoShape, shapes: Vec<TopoShape>, fuzzy: f64) -> Self {
        let mut builder = BuilderSolid::new();
        builder.set_shapes(shapes);
        builder.set_fuzzy(fuzzy);
        Self {
            solid,
            builder,
            warnings: Vec::new(),
        }
    }

    fn perform(&mut self) {
        match self.builder.perform() {
            Ok(()) => {
                self.warnings.extend(self.builder.warnings().iter().cloned());
            }
            Err(e) => {
                self.warnings.push(format!(
                    "BOPAlgo_SplitSolid: solid split failed: {e}"
                ));
            }
        }
    }

    fn areas(&self) -> &[TopoShape] {
        self.builder.areas()
    }
}

/// Shell-face set of a draft solid (`_3.cxx:494-499`).
fn draft_faces(draft: &TopoShape) -> Vec<TopoShape> {
    faces_of(draft).into_iter().map(|f| f.0).collect()
}

/// IN faces doubled FORWARD+REVERSED (`_3.cxx:502-511`).
fn doubled_in_faces(in_faces: &[TopoShape]) -> Vec<TopoShape> {
    let mut out = Vec::with_capacity(in_faces.len() * 2);
    for f in in_faces {
        out.extend(both_orientations(f));
    }
    out
}

/// `BOPAlgo_Builder::BuildSplitSolids` (`_3.cxx:413`).
pub fn build_split_solids_occt<B: BopSolidHost>(
    f: &mut B,
    fill: &FillIn3dPartsResult,
) -> Result<(), String> {
    let mut mst = BopToolsSetMap::new();
    let mut m_fence: HashSet<ShapeMapKey> = HashSet::new();

    // 0. Same-domain seeds for non-interfered solids.
    let n = nb_source(f.ds());
    for i in 0..n {
        let Some(si) = f.ds().shape_info(i) else { continue };
        if si.shape_type() != ShapeType::Solid {
            continue;
        }
        let s = si.shape().clone();
        if !m_fence.insert(shape_key(&s)) {
            continue;
        }
        if fill.draft_solids.contains_key(&shape_key(&s)) {
            continue;
        }
        let mut st = BopToolsSet::new();
        st.add(&s, ShapeType::Face);
        mst.add(st);
    }

    // 1. Build solids for interfered source solids.
    let mut solids_im: Vec<(TopoShape, Vec<TopoShape>)> = Vec::new();
    let mut vbs: Vec<SplitSolid> = Vec::new();
    let fuzzy = f.fuzzy_value();

    for i in 0..n {
        let Some(si) = f.ds().shape_info(i) else { continue };
        if si.shape_type() != ShapeType::Solid {
            continue;
        }
        let s = si.shape().clone();
        let Some(sd) = fill.draft_solids.get(&shape_key(&s)) else {
            continue;
        };
        let p_lfin = fill.in_parts.get(&shape_key(&s));
        if p_lfin.map(|v| v.is_empty()).unwrap_or(true) {
            solids_im.push((s, vec![sd.clone()]));
            continue;
        }
        let mut sfs = draft_faces(sd);
        sfs.extend(doubled_in_faces(p_lfin.unwrap()));
        vbs.push(SplitSolid::new(s, sfs, fuzzy));
    }

    // Sequential `BOPTools_Parallel::Perform`.
    for bs in &mut vbs {
        bs.perform();
        for w in &bs.warnings {
            f.add_warning(w.clone());
        }
        solids_im.push((bs.solid.clone(), bs.areas().to_vec()));
    }

    // Add new solids to images map (`_3.cxx:579-617`).
    for (s, lsr) in solids_im {
        if images_bound(f.history(), &s) {
            continue;
        }
        for sr in lsr {
            let (sx, b_flag_sd) = intern_solid_sd(&mut mst, &sr);
            f.history_mut().add_image(&s, sx.clone());
            origins_append(f.origins_mut(), &sx, s.clone());
            if b_flag_sd {
                f.bind_shapes_sd(sr, sx);
            }
        }
    }
    Ok(())
}
