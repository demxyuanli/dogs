//! `TopOpeBRepBuild_Builder::SplitSectionEdge` — one section edge through the
//! `G*` machinery.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_Section.cxx:928-1230`
//! (OCCT 8.0.0).
//!
//! The driver answers, for one section edge `EOR`, three split lists keyed by
//! the build state: `ON` (`GComUnsh(FACE, FACE)`), `IN` and `OUT`
//! (`GCutUnsh(FACE, FACE)`, the last two being the two permutations of one
//! table). Each is filled by the `G` pipeline: `GFillEdgePVS` on the forward
//! edge, a `PaveClassifier` over the resulting `PaveSet`, an `EdgeBuilder`
//! over the two, then `GEDBUMakeEdges`.
//!
//! Two shortcuts run before the `G` pipeline (`cxx:1025-1080`): an edge that
//! carries neither geometry nor a 3-D same-domain partner — or whose every
//! interference sits on a bound vertex — is answered by a copy of itself
//! (with its vertices rebound to their same-domain references when they have
//! one).
//!
//! UNPORTED, with its OCCT location: the `OCCT_DEBUG` dumps (`cxx:929-934`,
//! `:968-979`, `:981-987`, `:1085-1097`, `:1158-1162`, `:1195-1198`).

use crate::abs::{ShapeType, State};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, TopoShape, Vertex};
use crate::topo_tools::{edge_vertices, explorer_shapes, is_same};
use crate::topopebrep::build::edge_builder::EdgeBuilder;
use crate::topopebrep::build::g_tool::GTool;
use crate::topopebrep::build::pave_classifier::PaveClassifier;
use crate::topopebrep::build::pave_set::PaveSet;
use crate::topopebrep::ds::export::{fds_data, fds_has_same_domain2d, fds_has_same_domain3d};
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::tool::split_e::split_e;

use super::grid_ee::set_global_hassd;
use super::grid_ee::set_global_issp;
use super::Builder;

/// `#define TheIN (1)` (`Section.cxx:915`).
const THE_IN: i32 = 1;
/// `#define TheON (2)` (`Section.cxx:916`).
const THE_ON: i32 = 2;
/// `#define TheOUT (3)` (`Section.cxx:918`).
const THE_OUT: i32 = 3;

impl Builder {
    /// `SplitSectionEdge(EA)` (`Section.cxx:928-1230`).
    pub fn split_section_edge(&mut self, ea: &TopoShape) {
        let eor = ea.clone();
        let mut ef = eor.clone();
        ef.set_orientation(crate::abs::Orientation::Forward);

        let Some(hds) = self.data_structure() else {
            return;
        };
        let bds = hds.ds();
        // `cxx:941-943`: `AncestorRank(EF)`, set by `GMapShapes` (`Merge`).
        let rank_ef = bds.ancestor_rank_of(&ef);

        // `cxx:947-950`
        let hg = self
            .data_structure()
            .map(|h| h.has_geometry(&eor))
            .unwrap_or(false);
        let hsd3d = self
            .data_structure()
            .map(|h| fds_has_same_domain3d(h.ds(), &eor, false).0)
            .unwrap_or(false);
        let hsd2d = self
            .data_structure()
            .map(|h| fds_has_same_domain2d(h.ds(), &eor, false).0)
            .unwrap_or(false);

        // `cxx:953`
        let cond = hg || hsd3d;

        // `cxx:957-964`
        set_global_hassd(0);
        if hsd3d {
            set_global_hassd(3);
        }
        if hsd2d {
            set_global_hassd(2);
        }

        // `cxx:966-979`
        if self.my_split_section_edges_done {
            return;
        }

        // `cxx:989-1027`: every interference of EOR sits on a bound vertex.
        let all_gb1 = self.all_on_bound_vertices(&eor, rank_ef);

        // `cxx:1028-1030`
        let mut mke = !cond;
        mke = mke || !hg; // xpu271098 cto002J2(e9on)
        mke = mke || all_gb1; // xpu161198 BUC60382(e3on)

        if mke {
            self.mark_split(&eor, State::On, true);
            // `cxx:1037-1041`: `FUN_ds_CopyEdge(EOR, newEOR)`.
            let mut new_eor = self.my_build_tool.copy_edge(&eor);
            let mut hasnew_eor = false;
            for v in explorer_shapes(&eor, ShapeType::Vertex) {
                // `cxx:1047-1050`
                let par_v = BRepTool::parameter_of_vertex(&eor, &v);
                let (iv, ivref, hsd) = match self.data_structure() {
                    Some(hds) => {
                        let iv = hds.shape_index(&v, true);
                        let mut ivref = 0;
                        let mut hsd = false;
                        if iv != 0 {
                            hsd = hds.has_same_domain(&v, true);
                        }
                        if hsd {
                            ivref = hds.same_domain_reference(&v);
                        }
                        (iv, ivref, hsd)
                    }
                    None => (0, 0, false),
                };
                // `cxx:1057-1058`
                let setref = hsd && iv != ivref;
                let vref = self
                    .data_structure()
                    .and_then(|hds| hds.shape(ivref, true))
                    .and_then(|s| Vertex::wrap(s));
                let (Some(vref), true) = (vref, setref) else {
                    // `cxx:1060-1064`
                    self.my_build_tool.add_edge_vertex(&mut new_eor, &v);
                    self.my_build_tool.parameter(&new_eor, &v, par_v);
                    continue;
                };
                // `cxx:1067-1069`
                hasnew_eor = true;
                let mut vref = vref.0;
                vref.set_orientation(v.orientation());
                self.my_build_tool.add_edge_vertex(&mut new_eor, &vref);
                self.my_build_tool.parameter(&new_eor, &vref, par_v);
            }
            // `cxx:1071-1078`
            let new_list = if hasnew_eor { new_eor } else { eor.clone() };
            self.change_split(&eor, State::On).push(new_list);
            return;
        }

        // `cxx:1082-1083`: `GFindSamDom(EOR, LESD1, LESD2)` feeds only the
        // `OCCT_DEBUG` dumps at `cxx:1085-1097`, so nothing is computed here.

        // ----- (ON, SOLID) (`cxx:1098-1153`) -----
        {
            let g = GTool::g_com_unsh(ShapeType::Face, ShapeType::Face);
            self.my_edge_reference = Some(ef.clone());
            let mut pvs = PaveSet::new(&ef);

            set_global_issp(THE_ON);
            self.fill_edge_pvs(&ef, &g, &mut pvs);
            set_global_issp(0);

            let mut lem: Vec<TopoShape> = Vec::new();
            self.edge_builder_edges(&ef, &mut pvs, &mut lem);

            // xpu : 18-03-98 (cto 016 B2): split the edges of LEM that carry
            // internal vertices.
            let mut new_lem: Vec<TopoShape> = Vec::new();
            for esp in &lem {
                let mut lspe: Vec<TopoShape> = Vec::new();
                let ok = split_e(esp, &mut lspe);
                let nonwesp = !ok || lspe.len() < 2;
                if nonwesp {
                    new_lem.push(esp.clone());
                } else {
                    new_lem.extend(lspe);
                }
            }

            // `cxx:1147-1150`
            self.mark_split(&eor, State::On, true);
            let lon = self.change_split(&eor, State::On);
            lon.extend(new_lem);
        }

        // ----- (IN, OUT) 3d (`cxx:1155-1194`) -----
        {
            let g = GTool::g_cut_unsh(ShapeType::Face, ShapeType::Face).copy_permuted();
            self.my_edge_reference = Some(ef.clone());
            let mut pvs = PaveSet::new(&ef);

            set_global_issp(THE_IN);
            self.fill_edge_pvs(&ef, &g, &mut pvs);
            set_global_issp(0);

            let mut lemnc: Vec<TopoShape> = Vec::new();
            self.edge_builder_edges(&ef, &mut pvs, &mut lemnc);

            // `cxx:1189-1192`
            self.mark_split(&eor, State::In, true);
            let linn = self.change_split(&eor, State::In);
            linn.extend(lemnc);
        }

        // ----- (OUT, IN) 3d (`cxx:1196-1227`) -----
        {
            let g = GTool::g_cut_unsh(ShapeType::Face, ShapeType::Face);
            self.my_edge_reference = Some(ef.clone());
            let mut pvs = PaveSet::new(&ef);

            set_global_issp(THE_OUT);
            self.fill_edge_pvs(&ef, &g, &mut pvs);
            set_global_issp(0);

            let mut lem: Vec<TopoShape> = Vec::new();
            self.edge_builder_edges(&ef, &mut pvs, &mut lem);

            // `cxx:1224-1227`
            self.mark_split(&eor, State::Out, true);
            let lout = self.change_split(&eor, State::Out);
            lout.extend(lem);
        }

        // `cxx:1229`
        set_global_hassd(0);
    }

    /// `cxx:1098-1121`, `:1163-1186`, `:1199-1222`: the `PaveSet` -> classifier
    /// -> `EdgeBuilder` -> `GEDBUMakeEdges` pipeline shared by the three state
    /// blocks.
    fn edge_builder_edges(
        &self,
        ef: &TopoShape,
        pvs: &mut PaveSet,
        lem: &mut Vec<TopoShape>,
    ) {
        let mut vcl = PaveClassifier::new(ef);
        // `cxx:1114-1118`
        if pvs.has_equal_parameters() {
            vcl.set_first_parameter(pvs.equal_parameters());
        }
        let mut edbu = EdgeBuilder::new();
        edbu.init_edge_builder(pvs, &mut vcl, false);
        self.gedbu_make_edges(ef, &mut edbu, lem);
    }

    /// `cxx:989-1027`: every interference of `eor` is an edge/vertex
    /// interference whose vertex sits on one of the two end vertices of `eor`,
    /// either directly (same rank) or through its same-domain partner.
    fn all_on_bound_vertices(&self, eor: &TopoShape, rank_ef: i32) -> bool {
        let Some(hds) = self.data_structure() else {
            return false;
        };
        let bds = hds.ds();
        let edge = Edge(eor.clone());
        let (vf, vl) = edge_vertices(&edge);
        let on_end = |v: &TopoShape| {
            vf.as_ref().is_some_and(|x| is_same(&x.0, v))
                || vl.as_ref().is_some_and(|x| is_same(&x.0, v))
        };

        let mut all_gb1 = false;
        for handle in bds.shape_interferences(eor, true) {
            let it = handle.read().expect("poisoned Interference lock");
            let (gt, g, _, _) = fds_data(&it);
            if gt == Kind::Point {
                all_gb1 = false;
                break;
            }
            let rkg = bds.ancestor_rank(g);
            let Some(v) = hds.shape(g, true) else {
                all_gb1 = false;
                break;
            };
            if rkg == rank_ef {
                all_gb1 = on_end(&v);
            } else {
                // `cxx:1013-1021`: `FUN_ds_getoov(v, HDS, oov)`. OCCT compares
                // an unset `oov` when the call answers false; the port treats
                // that as "not on a bound".
                all_gb1 = match crate::topopebrep::ds::export::fun_ds_getoov_in_hds(hds, &v) {
                    Some(oov) => on_end(&oov),
                    None => false,
                };
            }
            if !all_gb1 {
                break;
            }
        }
        all_gb1
    }
}
