//! `BRepFeat_Builder::RebuildFaces` (`BRepFeat_Builder.cxx:271-558`).
//!
//! The kept parts of the tool (`myShapes`) decide which images of a basis face
//! survive: an image edge that keeps a vertex inside `myShapes` stays, the rest
//! go to `aME`. A boundary edge whose images are partly kept and partly removed
//! is re-split by [`BRepFeatBuilder::rebuild_edge`], and the resulting edge set
//! is handed to `BOPAlgo_BuilderFace` to rebuild the face split. The areas are
//! then matched against the same-domain set of the kept faces, so identical
//! topologies collapse onto one representative.

use crate::abs::{Orientation, ShapeType};
use crate::algo_tools_face::is_split_to_reverse;
use crate::bop_split_seam::do_split_seam_on_face;
use crate::bop_tools_set::{BopToolsSet, BopToolsSetMap};
use crate::brep_tool::BRepTool;
use crate::builder_area::AreaBuilder;
use crate::builder_face::FaceBuilder;
use crate::int_tools_full::IntToolsContext;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{edges_of, explorer_shapes, shape_map_key};

use super::{BRepFeatBuilder, ShapeSet};

/// `RebuildFaces()` (`cxx:271-558`).
pub(super) fn rebuild_faces(b: &mut BRepFeatBuilder) -> Result<(), String> {
    // `cxx:288-298`: the edge sets of the kept faces.
    let mut a_mst = BopToolsSetMap::new();
    for a_s in b.my_shapes.iter() {
        if a_s.shape_type() == ShapeType::Face {
            a_mst.add(BopToolsSet::from_shape(a_s, ShapeType::Edge));
        }
    }

    // `cxx:282`: the edges already re-split.
    let mut a_me_split = ShapeSet::new();
    let mut ctx = IntToolsContext::new();

    // `cxx:300-556`.
    let a_nb_s = b.bop().ds().nb_source_shapes();
    for i in 0..a_nb_s {
        // `cxx:304-325`: the images of the tool shapes keep only the parts
        // inside `myShapes`.
        if b.bop().ds().rank(i) == 1 {
            let Some(a_s) = b.bop().ds().shape(i).cloned() else {
                continue;
            };
            if b.bop().history().has_image(&a_s) {
                let images = b
                    .bop()
                    .history()
                    .image(&a_s)
                    .map(|v| v.to_vec())
                    .unwrap_or_default();
                let mut kept: Vec<TopoShape> = Vec::new();
                for a_s_im in images {
                    if b.my_shapes.contains(&a_s_im) {
                        kept.push(a_s_im);
                    }
                }
                b.bop().history_mut().set_images(&a_s, kept);
            }
            continue;
        }

        // `cxx:326-332`.
        let Some(a_s) = b.bop().ds().shape(i).cloned() else {
            continue;
        };
        if a_s.shape_type() != ShapeType::Face {
            continue;
        }

        // `cxx:334-338`.
        if !b.bop().history().has_image(&a_s) {
            continue;
        }

        // `cxx:340-345`.
        let an_ori_f = a_s.orientation();
        let mut a_ff = a_s.clone();
        a_ff.set_orientation(Orientation::Forward);
        let a_ff_face = Face::wrap(a_ff.clone())
            .ok_or_else(|| "Standard_TypeMismatch: TopoDS::Face".to_string())?;

        // `cxx:347-348`: `PaveBlocksIn` / `PaveBlocksSc`.
        let (paves_in, paves_sc) = match b.bop().ds().face_info(i) {
            Some(fi) => (fi.paves_in().to_vec(), fi.paves().to_vec()),
            None => (Vec::new(), Vec::new()),
        };

        let mut a_le: Vec<TopoShape> = Vec::new();

        // `cxx:350-474`: the bounding edges.
        for be in edges_of(&a_ff) {
            let an_ori_e = be.0.orientation();
            let b_is_degenerated = BRepTool::is_degenerated(&be);
            let b_is_closed = BRepTool::is_closed_edge_face(&be, &a_ff_face);
            if !b.bop().history().has_image(&be.0) {
                // `cxx:471-473`.
                a_le.push(be.0.clone());
                continue;
            }

            // `cxx:359-406`.
            let a_le_im = b
                .bop()
                .history()
                .image(&be.0)
                .map(|v| v.to_vec())
                .unwrap_or_default();
            let mut b_rem = false;
            let mut b_im = false;
            let mut a_me = ShapeSet::new();
            let mut a_le_im_new: Vec<TopoShape> = Vec::new();
            for a_sx in &a_le_im {
                // `cxx:373-389`: the image keeps the edge when the edge or one
                // of its vertices is inside `myShapes`.
                let mut b_v_in_shapes = b.my_shapes.contains(a_sx);
                if !b_v_in_shapes {
                    for v in explorer_shapes(a_sx, ShapeType::Vertex) {
                        if b.my_shapes.contains(&v) {
                            b_v_in_shapes = true;
                            break;
                        }
                    }
                }
                if b_v_in_shapes {
                    b_im = true;
                    a_le_im_new.push(a_sx.clone());
                } else {
                    b_rem = true;
                    a_me.add(a_sx);
                }
            }

            // `cxx:407-411`.
            if !b_im {
                a_le.push(be.0.clone());
                continue;
            }

            // `cxx:412-428`.
            let mut a_le_im = a_le_im;
            if b_rem && b_im {
                if a_le_im.len() == 2 {
                    a_le.push(be.0.clone());
                    continue;
                }
                if !a_me_split.contains(&be.0) {
                    a_me_split.add(&be.0);
                    b.rebuild_edge(&be.0, &a_ff_face, &a_me, &mut a_le_im_new)?;
                    a_le_im = a_le_im_new;
                    if a_le_im.len() == 1 {
                        a_le.push(be.0.clone());
                        continue;
                    }
                }
            }

            // `cxx:429-470`.
            for a_sp in a_le_im {
                let mut a_sp = a_sp;

                // `cxx:433-438`.
                if b_is_degenerated {
                    a_sp.set_orientation(an_ori_e);
                    a_le.push(a_sp);
                    continue;
                }

                // `cxx:439-448`.
                if an_ori_e == Orientation::Internal {
                    a_sp.set_orientation(Orientation::Forward);
                    a_le.push(a_sp.clone());
                    a_sp.set_orientation(Orientation::Reversed);
                    a_le.push(a_sp);
                    continue;
                }

                // `cxx:449-461`.
                if b_is_closed {
                    let sp_edge = Edge(a_sp.clone());
                    if !BRepTool::is_closed_edge_face(&sp_edge, &a_ff_face) {
                        do_split_seam_on_face(&sp_edge, &a_ff_face);
                    }
                    a_sp.set_orientation(Orientation::Forward);
                    a_le.push(a_sp.clone());
                    a_sp.set_orientation(Orientation::Reversed);
                    a_le.push(a_sp);
                    continue;
                }

                // `cxx:463-470`. A failed analysis answers `false` here, which
                // is what the reference's `IsSplitToReverse` yields on a
                // geometry failure.
                a_sp.set_orientation(an_ori_e);
                let b_to_reverse =
                    is_split_to_reverse(&a_sp, &be.0, &mut ctx).unwrap_or(false);
                if b_to_reverse {
                    let ori = a_sp.orientation().reversed();
                    a_sp.set_orientation(ori);
                }
                a_le.push(a_sp);
            }
        }

        // `cxx:477-495`: the IN edges.
        for (n_sp, _, _) in paves_in {
            let Some(a_sp) = b.bop().ds().shape(n_sp).cloned() else {
                continue;
            };
            if b.my_removed.contains(&a_sp) {
                continue;
            }
            let mut a_sp = a_sp;
            a_sp.set_orientation(Orientation::Forward);
            a_le.push(a_sp.clone());
            a_sp.set_orientation(Orientation::Reversed);
            a_le.push(a_sp);
        }

        // `cxx:496-511`: the section edges.
        for (n_sp, _, _) in paves_sc {
            let Some(a_sp) = b.bop().ds().shape(n_sp).cloned() else {
                continue;
            };
            if b.my_removed.contains(&a_sp) {
                continue;
            }
            let mut a_sp = a_sp;
            a_sp.set_orientation(Orientation::Forward);
            a_le.push(a_sp.clone());
            a_sp.set_orientation(Orientation::Reversed);
            a_le.push(a_sp);
        }

        // `cxx:513-519`: `BOPAlgo_BuilderFace`.
        let mut a_bf = FaceBuilder::new();
        a_bf.set_face(&a_ff_face);
        a_bf.set_shapes(&a_le);
        crate::bop_builder_face_perform::perform_face_builder(&mut a_bf)?;

        // `cxx:522-551`.
        let a_lfr = a_bf.areas().to_vec();
        let mut a_lf_im: Vec<TopoShape> = Vec::new();
        for a_fr in a_lfr {
            let a_st = BopToolsSet::from_shape(&a_fr, ShapeType::Edge);
            let b_flag_sd = a_mst.contains(&a_st);
            let mut a_sx = a_mst.added(a_st).shape().clone();
            a_sx.set_orientation(an_ori_f);
            a_lf_im.push(a_sx.clone());

            // `cxx:540-546`: `myOrigins(aSx).Append(aF)`.
            let key = shape_map_key(&a_sx);
            b.bop()
                .origins_mut()
                .entry(key)
                .or_insert_with(Vec::new)
                .push(a_s.clone());

            // `cxx:547-550`.
            if b_flag_sd {
                b.bop().bind_shapes_sd(a_fr.clone(), a_sx.clone());
            }
        }

        // `cxx:523-524` and `cxx:552-555`.
        if a_lf_im.is_empty() {
            b.bop().history_mut().unbind_image(&a_s);
        } else {
            b.bop().history_mut().set_images(&a_s, a_lf_im);
        }
    }

    let _ = a_me_split;
    Ok(())
}
