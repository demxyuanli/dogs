//! The pairwise rebuild pass `IntersectFaces` (`cxx:7233-7311`) and the
//! invalid-split filter `GetInvalidEdges` (`cxx:7607-7676`).
//!
//! `IntersectFaces(theFInv, theFi, theFj, theLFInv, theLFImi, theLFImj,
//! theLFEi, theLFEj, theMECV, theMEToInt)` intersects two splits through
//! `BRepOffset_Tool::Inter3D` (`cxx:7250-7258`) and returns at once when the
//! faces do not cross (`cxx:7259-7262`). The common vertices of the two split
//! lists are added to `theMECV` only when there is more than one of them
//! (`cxx:7265-7275`), the origins of the parts shared with the invalid face are
//! looked up (`cxx:7277-7297`), and every intersection edge goes into
//! `theLFEi`, `theLFEj`, `myEdgesOrigins` (a `Bind`, i.e. a replace,
//! `cxx:7307-7309`) and `theMEToInt` (`cxx:7311`).
//!
//! `GetInvalidEdges(theVertsToAvoid, theMVBounds, theGF, theMEInv)` walks the
//! vertices of the builder result (`cxx:7620-7624`) and marks a vertex invalid
//! when it is neither a bounding vertex (`cxx:7639-7642`) nor an old shape
//! (`cxx:7645-7648`) and coincides with a vertex to avoid
//! (`cxx:7652-7660`); all the edges of such a vertex become invalid
//! (`cxx:7663-7670`).

use crate::abs::{ShapeType, State};
use crate::bop_builder2::BopBuilder;
use crate::brep_offset::inter3d::ShapeIndexedMap;
use crate::brep_offset::BRepOffsetTool;
use crate::int_tools_lines::compute_vv;
use crate::shape::{Face, TopoShape, Vertex};
use crate::topo_tools::{explorer_shapes, is_same, shape_map_key};

use super::common_edges::find_origins;
use super::statics::append_to_list;
use super::BRepOffsetBuildOffsetFaces;

impl BRepOffsetBuildOffsetFaces {
    /// `void BRepOffset_BuildOffsetFaces::IntersectFaces(const TopoDS_Shape&
    /// theFInv, const TopoDS_Shape& theFi, const TopoDS_Shape& theFj, const
    /// List& theLFInv, const List& theLFImi, const List& theLFImj, List&
    /// theLFEi, List& theLFEj, IndexedMap& theMECV, IndexedMap& theMEToInt)`
    /// (`cxx:7233-7311`).
    #[allow(clippy::too_many_arguments)]
    pub fn intersect_faces_of_pair(
        &mut self,
        the_f_inv: &TopoShape,
        the_fi: &TopoShape,
        the_fj: &TopoShape,
        the_lf_inv: &[TopoShape],
        the_lf_imi: &[TopoShape],
        the_lf_imj: &[TopoShape],
        the_lfei: &mut Vec<TopoShape>,
        the_lfej: &mut Vec<TopoShape>,
        the_mecv: &mut ShapeIndexedMap,
        the_me_to_int: &mut ShapeIndexedMap,
    ) {
        // `cxx:7244-7258`: `BRepOffset_Tool::Inter3D` with a null reference
        // edge and two null reference faces.
        let mut a_l_int1: Vec<TopoShape> = Vec::new();
        let mut a_l_int2: Vec<TopoShape> = Vec::new();
        BRepOffsetTool::inter3d(
            &Face(the_fi.clone()),
            &Face(the_fj.clone()),
            &mut a_l_int1,
            &mut a_l_int2,
            State::Out,
            None,
            None,
            None,
        );

        // `cxx:7259-7262`.
        if a_l_int1.is_empty() {
            return;
        }

        // `cxx:7265-7275`: the common vertices are used for trimming only when
        // there is more than one.
        let mut a_lcv: Vec<TopoShape> = Vec::new();
        super::statics::find_common_parts(
            the_lf_imi,
            the_lf_imj,
            &mut a_lcv,
            ShapeType::Vertex,
        );
        if a_lcv.len() > 1 {
            for a_cv in &a_lcv {
                the_mecv.add(a_cv);
            }
        }

        // `cxx:7277-7287`: the parts of the invalid splits each side shares
        // with the invalidity.
        let mut a_le_to_find_origins: Vec<TopoShape> = Vec::new();
        if !is_same(the_fi, the_f_inv) {
            super::statics::find_common_parts(
                the_lf_imi,
                the_lf_inv,
                &mut a_le_to_find_origins,
                ShapeType::Edge,
            );
        }
        if !is_same(the_fj, the_f_inv) {
            super::statics::find_common_parts(
                the_lf_imj,
                the_lf_inv,
                &mut a_le_to_find_origins,
                ShapeType::Edge,
            );
        }

        // `cxx:7289-7297`.
        let mut a_me_to_find_origins = ShapeIndexedMap::new();
        for a_ec in &a_le_to_find_origins {
            a_me_to_find_origins.add(a_ec);
        }
        let mut a_le_or_init: Vec<TopoShape> = Vec::new();
        find_origins(
            the_lf_imi,
            the_lf_imj,
            &a_me_to_find_origins,
            &self.my_edges_origins,
            &mut a_le_or_init,
        );

        // `cxx:7299-7311`.
        for a_e_int in &a_l_int1 {
            append_to_list(the_lfei, a_e_int);
            append_to_list(the_lfej, a_e_int);
            // `cxx:7307-7309`: `myEdgesOrigins->Bind` replaces the list.
            if !a_le_or_init.is_empty() {
                self.my_edges_origins
                    .insert(shape_map_key(a_e_int), a_le_or_init.clone());
            }
            the_me_to_int.add(a_e_int);
        }
    }

    /// `void BRepOffset_BuildOffsetFaces::GetInvalidEdges(const Map&
    /// theVertsToAvoid, const Map& theMVBounds, BOPAlgo_Builder& theGF,
    /// Map& theMEInv)` (`cxx:7607-7676`).
    pub fn get_invalid_edges(
        &self,
        the_verts_to_avoid: &ShapeIndexedMap,
        the_mv_bounds: &ShapeIndexedMap,
        the_gf: &BopBuilder,
        the_me_inv: &mut ShapeIndexedMap,
    ) {
        // `cxx:7612-7615`.
        if the_verts_to_avoid.extent() == 0 {
            return;
        }

        // `cxx:7620-7624`: `MapShapesAndAncestors(theGF.Shape(), VERTEX,
        // EDGE)`, kept in the order the reference's `FindKey(i)` walks.
        let a_res = the_gf.result();
        let mut a_dmve_keys: Vec<TopoShape> = Vec::new();
        let mut a_dmve: Vec<Vec<TopoShape>> = Vec::new();
        for a_e in explorer_shapes(a_res, ShapeType::Edge) {
            for a_v in explorer_shapes(&a_e, ShapeType::Vertex) {
                let key = shape_map_key(&a_v);
                match a_dmve_keys.iter().position(|k| shape_map_key(k) == key) {
                    Some(i) => a_dmve[i].push(a_e.clone()),
                    None => {
                        a_dmve_keys.push(a_v);
                        a_dmve.push(vec![a_e.clone()]);
                    }
                }
            }
        }

        let p_ds = the_gf.ds();

        // `cxx:7635-7671`.
        for (a_v, a_lve) in a_dmve_keys.iter().zip(a_dmve.iter()) {
            // `cxx:7639-7642`: the bounding vertices are never invalid.
            if the_mv_bounds.contains(a_v) {
                continue;
            }
            // `cxx:7645-7648`: an old shape is never invalid. `pDS->Index`
            // answers a negative index for an unregistered shape, which the
            // reference lets through.
            if let Some(n_v) = p_ds.index(a_v) {
                if !p_ds.is_new_shape(n_v) {
                    continue;
                }
            }

            // `cxx:7650-7660`: a vertex coinciding with a vertex to avoid is
            // invalid.
            let mut b_invalid = false;
            for a_v_avoid in the_verts_to_avoid.iter() {
                let i_flag = compute_vv(
                    &Vertex(a_v.clone()),
                    &Vertex(a_v_avoid.clone()),
                );
                if i_flag == 0 {
                    b_invalid = true;
                    break;
                }
            }

            // `cxx:7663-7670`.
            if b_invalid {
                for a_e in a_lve {
                    the_me_inv.add(a_e);
                }
            }
        }
    }
}
