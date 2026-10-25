//! `TopOpeBRepBuild_Builder::GWESMakeFaces`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_makefaces.cxx`
//! (OCCT 8.0.0, 85-246).
//!
//! Virtual (`Builder.hxx`), overridden by `Builder1::GWESMakeFaces`
//! (`Builder1.cxx:1373-1412`, ported in
//! `builder/builder1/gwes_make_faces.rs`); [`Builder::gwes_make_faces`] branches
//! on [`Builder::is_builder1`] where the vtable would land and runs
//! [`Builder::base_gwes_make_faces`] otherwise.
//!
//! The call is the whole face-building sequence of a `WireEdgeSet`:
//! `FaceBuilder` on the set, unclosed-wire correction, `GFABUMakeFaces`, the
//! closing-edge purge (`TopOpeBRepTool::PurgeClosingEdges` then
//! `TopOpeBRepTool::MakeFaces`, `TopOpeBRepTool_PURGE.cxx:405-427`,
//! `:974-1038`), the p-curve correction on a closed face
//! (`TopOpeBRepTool::CorrectONUVISO`) and finally
//! `RegularizeFaces` (`FREGU.cxx:39-131`).
//!
//! UNPORTED:
//! - the `OCCT_DEBUG` blocks (`DEBSHASET`, `GdumpSHA`, `WES.DumpSS`,
//!   `debwesmf`, `GetcontextNOPURGE` / `GetcontextNOCORRISO` / `debcorriso`).
//!   `topurge` and `corronISO` start at `true` (`cxx:114`, `:214`).

use std::collections::HashMap;

use crate::shape::Face;
use crate::shape::TopoShape;
use crate::topo_tools::ShapeMapKey;
use crate::topopebrep::build::face_builder::FaceBuilder;
use crate::topopebrep::build::wire_edge_set::WireEdgeSet;
use crate::topopebrep::ds::data_structure::IndexedShapeData;
use crate::topopebrep::ds::data_structure::IndexedShapeSet;
use crate::topopebrep::tool::geometry::fun_tool_is_closed_s;
use crate::topopebrep::tool::purge::{
    correct_on_uv_iso, make_faces as tool_make_faces, purge_closing_edges_faces,
};

use super::Builder;

impl Builder {
    /// `GWESMakeFaces(FF, WES, LOF)` (`cxx:85-246`).
    ///
    /// Virtual in OCCT (`Builder.hxx`), overridden by
    /// `Builder1::GWESMakeFaces` (`Builder1.cxx:1373-1412`); the branch on
    /// [`Builder::is_builder1`] is where the vtable would land.
    pub fn gwes_make_faces(
        &mut self,
        ff: &TopoShape,
        wes: &mut WireEdgeSet,
        lof: &mut Vec<TopoShape>,
    ) {
        if self.is_builder1 {
            self.b1_gwes_make_faces(ff, wes, lof);
            return;
        }
        self.base_gwes_make_faces(ff, wes, lof);
    }

    /// `TopOpeBRepBuild_Builder::GWESMakeFaces(FF, WES, LOF)` (`cxx:85-246`).
    ///
    /// `lof` is an out-list: OCCT appends to it without clearing (`cxx:168`),
    /// and the port does the same.
    pub(crate) fn base_gwes_make_faces(
        &mut self,
        ff: &TopoShape,
        wes: &mut WireEdgeSet,
        lof: &mut Vec<TopoShape>,
    ) {
        let force_class = true;
        let mut fabu = FaceBuilder::new();
        fabu.init_face_builder(wes, ff, force_class);

        // Wire checking (`cxx:107-160`): a wire that has no edge connected to
        // the others in the face's UV representation is repaired by removing
        // the lonesome closing edge.
        let mut topurge = true;

        if topurge {
            let mut map_pie = IndexedShapeSet::new();
            fabu.detect_pseudo_internal_edge(&mut map_pie);

            let mut map_vv_same_g: IndexedShapeData<TopoShape> = IndexedShapeData::new();
            let mut map_v_on1_edge: IndexedShapeData<TopoShape> = IndexedShapeData::new();
            fabu.detect_unclosed_wire(&mut map_vv_same_g, &mut map_v_on1_edge);

            let n_vv = map_vv_same_g.extent();
            let mut map_vv_ref: IndexedShapeData<TopoShape> = IndexedShapeData::new();
            if n_vv > 0 {
                // `cxx:128-158`: update the DS with the same-domain vertices
                // and fill `mapVVref`.
                for i in 1..=n_vv {
                    let Some(v) = map_vv_same_g.find_key(i).cloned() else {
                        continue;
                    };
                    let Some(hds) = self.my_data_structure.as_mut() else {
                        continue;
                    };
                    let mut hsdm = hds.has_same_domain(&v, false);
                    if !hsdm {
                        let rank_v = hds.change_ds().ancestor_rank_of(&v);
                        let Some(v_same_g) = map_vv_same_g.find_from_index(i).cloned() else {
                            continue;
                        };
                        // MSV Oct 4, 2001: prefer old vertex as
                        // SameDomainReference
                        let rank_v_same_g = hds.change_ds().ancestor_rank_of(&v_same_g);
                        let other_ref = rank_v_same_g != 0 && rank_v != 1;
                        let _ = if other_ref {
                            hds.change_ds()
                                .fill_shapes_same_domain(&v_same_g, &v, true)
                        } else {
                            hds.change_ds().fill_shapes_same_domain(&v, &v_same_g, true)
                        };
                        hsdm = hds.has_same_domain(&v, false);
                    }
                    if hsdm {
                        let i_ref = hds.same_domain_reference(&v);
                        let v_ref = hds.shape(i_ref, false);
                        if let Some(v_ref) = v_ref {
                            map_vv_ref.add(&v, v_ref);
                        }
                    }
                }
                fabu.correct_gclosed_wire(&map_vv_ref, &map_v_on1_edge);
                fabu.detect_unclosed_wire(&mut map_vv_same_g, &mut map_v_on1_edge);
            }
        }

        // `cxx:162-164`
        let mut m_wis_old: HashMap<ShapeMapKey, i32> = HashMap::new();
        self.g_fabumake_faces(ff, &mut fabu, lof, &mut m_wis_old);

        // 2. on the periodic face F: find the faulty shapes to avoid when
        // building up the shapes (`cxx:166-203`).
        let mut msh_nok: Vec<TopoShape> = Vec::new();
        if topurge {
            let fa = ff.clone();
            let puok = purge_closing_edges_faces(&fa, lof, &mut msh_nok);
            if !puok {
                // `cxx:185`: `throw Standard_Failure("TopOpeBRepBuild::GWESMakeFaces")`
                panic!("TopOpeBRepBuild::GWESMakeFaces");
            }
            topurge = !msh_nok.is_empty();
        }

        if topurge {
            let loff = tool_make_faces(ff, lof, &msh_nok);
            lof.clear();
            lof.extend(loff);
        }

        // 1. on the periodic face F: translate the p-curves of the edges into
        // F's UV bounds (`cxx:206-233`).
        let corron_iso = fun_tool_is_closed_s(&Face(ff.clone()));
        if corron_iso {
            let mut new_lof: Vec<TopoShape> = Vec::new();
            let ffa = ff.clone();
            for fa in lof.iter() {
                let mut fa = fa.clone();
                correct_on_uv_iso(&ffa, &mut fa);
                new_lof.push(fa);
            }
            lof.clear();
            lof.extend(new_lof);
        }

        // xpu280898: regularization after GFABUMakeFaces, purge processings
        // (`cxx:236-241`)
        let mut new_lof: Vec<TopoShape> = Vec::new();
        self.regularize_faces(ff, lof, &mut new_lof);
        lof.clear();
        lof.extend(new_lof);
    }
}
