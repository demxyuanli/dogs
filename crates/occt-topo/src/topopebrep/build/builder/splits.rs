//! `TopOpeBRepBuild_Builder` — the split drivers.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_Builder.cxx:150-179`
//! (`AddIntersectionEdges`), `:925-1153` (`SplitEdge`, `SplitEdge1`,
//! `SplitEdge2`), `:1155-1533` (`SplitFace`, `SplitFace1`, `SplitFace2`),
//! `:1535-1711` (`SplitSolid`), `:1713-1874` (`SplitShapes`) (OCCT 8.0.0).
//!
//! A split driver collects, for one shape of one operand, the same-domain
//! partners on the other operand, fills a `PaveSet` / `WireEdgeSet` /
//! `ShellFaceSet` with the intersection geometry plus the untouched
//! sub-shapes, and hands the set to `EdgeBuilder` / `FaceBuilder` /
//! `SolidBuilder`. The built sub-shapes land in `ChangeMerged`, and every
//! same-domain shape of either operand is marked split and, when the
//! connectivity table of the operation asks for it, receives the same list in
//! `ChangeSplit`.
//!
//! `SplitEdge` and `SplitFace` dispatch on `TopOpeBRepBuild_GetcontextSF2()` in
//! OCCT (`cxx:927-935`, `:1158-1166`); that is a debug-only switch
//! (`#ifdef OCCT_DEBUG`), so a release build always runs the `1` variant and
//! the port does the same.

use std::sync::atomic::{AtomicI32, Ordering};

use occt_core::precision::CONFUSION;

use crate::abs::{Orientation, ShapeType, State};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Solid, TopoShape};
use crate::topo_tools::{edge_vertices, explorer_shapes};
use crate::topopebrep::ds::data_structure::DataStructure;
use crate::topopebrep::tool::shape_explorer::ShapeExplorer;
use crate::topopebrep::tool::solid_classifier::SolidClassifier;

use super::super::edge_builder::EdgeBuilder;
use super::super::face_builder::FaceBuilder;
use super::super::loop_set::LoopSetSource;
use super::super::pave_classifier::PaveClassifier;
use super::super::pave_set::PaveSet;
use super::super::shape_set::ShapeSetSource;
use super::super::shell_face_set::ShellFaceSet;
use super::super::solid_builder::SolidBuilder;
use super::super::wire_edge_set::WireEdgeSet;
use super::Builder;

/// `STATIC_SOLIDINDEX` (`Builder.cxx:1609`, `:1618`, `:1828-1835`): the
/// operand `SplitShapes` classifies a touched edge against.
static STATIC_SOLID_INDEX: AtomicI32 = AtomicI32::new(0);

impl Builder {
    /// `SplitEdge(E, ToBuild1, ToBuild2)` (`cxx:925-936`).
    pub fn split_edge(&mut self, e: &TopoShape, to_build1: State, to_build2: State) {
        self.split_edge1(e, to_build1, to_build2);
    }

    /// `SplitEdge1(Eoriented, ToBuild1, ToBuild2)` (`cxx:941-1079`).
    pub fn split_edge1(&mut self, e_oriented: &TopoShape, to_build1: State, to_build2: State) {
        // cxx:947-948: work on a FORWARD edge <Eforward>.
        let mut e_forward = e_oriented.clone();
        e_forward.set_orientation(Orientation::Forward);

        // cxx:950
        let tosplit = self.to_split(e_oriented, to_build1);
        if !tosplit {
            return;
        }

        // cxx:969-971
        let rev_ori1 = Self::reverse(to_build1, to_build2);
        let rev_ori2 = Self::reverse(to_build2, to_build1);
        let _ = (rev_ori1, rev_ori2);
        let connect_to1 = true;
        let connect_to2 = false;

        // cxx:975-977: the edges to split on both sides.
        let mut le1 = vec![e_forward.clone()];
        let mut le2: Vec<TopoShape> = Vec::new();
        self.find_same_domain(&mut le1, &mut le2);

        // cxx:1010-1017: a `PaveSet` on <Eforward> with every point or vertex
        // of the DS found on it.
        let mut pvs = PaveSet::new(&e_forward);
        let mut epit = match self.data_structure() {
            Some(hds) => hds.edge_points(&e_forward),
            None => return,
        };
        self.fill_vertex_set(&mut epit, to_build1, &mut pvs);

        // cxx:1018-1022
        let mut vcl = PaveClassifier::new(&e_forward);
        let equalpar = pvs.has_equal_parameters();
        if equalpar {
            let p = pvs.equal_parameters();
            vcl.set_first_parameter(p);
        }

        // cxx:1027-1041: mark <Eforward> as split <ToBuild1> before leaving,
        // even when there is no vertex to split at.
        self.mark_split(&e_forward, to_build1, true);
        pvs.init_loop();
        if !pvs.more_loop() {
            return;
        }

        // cxx:1044-1049
        let mut ebu = EdgeBuilder::with_loop_set(&mut pvs, &mut vcl, false);
        let mut edge_list = self.take_merged(&e_forward, to_build1);
        self.make_edges(&e_forward, &mut ebu, &mut edge_list);
        self.put_merged(&e_forward, to_build1, edge_list.clone());

        // cxx:1053-1076: connect the new edges to the LE1 / LE2 edges.
        for e_cur in &le1 {
            self.mark_split(e_cur, to_build1, true);
            let el = self.change_split(e_cur, to_build1);
            if connect_to1 {
                *el = edge_list.clone();
            }
        }
        for e_cur in &le2 {
            self.mark_split(e_cur, to_build2, true);
            let el = self.change_split(e_cur, to_build2);
            if connect_to2 {
                *el = edge_list.clone();
            }
        }
    }

    /// `SplitEdge2(Eoriented, ToBuild1, ToBuild2)` (`cxx:1083-1150`).
    pub fn split_edge2(&mut self, e_oriented: &TopoShape, to_build1: State, _to_build2: State) {
        // cxx:1088-1090
        let tosplit = self.to_split(e_oriented, to_build1);
        if !tosplit {
            return;
        }

        // cxx:1094-1096
        let mut e_forward = e_oriented.clone();
        self.my_build_tool
            .set_orientation(&mut e_forward, Orientation::Forward);

        // cxx:1110-1117
        let mut pvs = PaveSet::new(&e_forward);
        let mut epit = match self.data_structure() {
            Some(hds) => hds.edge_points(&e_forward),
            None => return,
        };
        self.fill_vertex_set(&mut epit, to_build1, &mut pvs);

        // cxx:1118-1122
        let mut vcl = PaveClassifier::new(&e_forward);
        let equalpar = pvs.has_equal_parameters();
        if equalpar {
            let p = pvs.equal_parameters();
            vcl.set_first_parameter(p);
        }

        // cxx:1128-1140
        self.mark_split(&e_forward, to_build1, true);
        pvs.init_loop();
        if !pvs.more_loop() {
            return;
        }

        // cxx:1143-1149
        let mut ebu = EdgeBuilder::with_loop_set(&mut pvs, &mut vcl, false);
        let mut edge_list = self.take_merged(&e_forward, to_build1);
        self.make_edges(&e_forward, &mut ebu, &mut edge_list);
        self.put_merged(&e_forward, to_build1, edge_list);
    }

    /// `SplitFace(Foriented, ToBuild1, ToBuild2)` (`cxx:1155-1167`).
    pub fn split_face(&mut self, f: &TopoShape, to_build1: State, to_build2: State) {
        self.split_face1(f, to_build1, to_build2);
    }

    /// `SplitFace1(Foriented, ToBuild1, ToBuild2)` (`cxx:1171-1300`).
    ///
    /// The connectivity table of `cxx:1176-1182` (common: connect to 1 and 2,
    /// fuse: both, cut: 1 only) is what `ConnectTo1 = true` / `ConnectTo2 =
    /// false` encode.
    pub fn split_face1(&mut self, f_oriented: &TopoShape, to_build1: State, to_build2: State) {
        // cxx:1184-1187
        let tosplit = self.to_split(f_oriented, to_build1);
        if !tosplit {
            return;
        }

        // cxx:1189-1192
        let mut rev_ori1 = Self::reverse(to_build1, to_build2);
        let mut rev_ori2 = Self::reverse(to_build2, to_build1);
        let connect_to1 = true;
        let connect_to2 = false;

        // cxx:1194-1196
        let mut f_forward = f_oriented.clone();
        self.my_build_tool
            .set_orientation(&mut f_forward, Orientation::Forward);

        // cxx:1198-1202
        let mut lf1 = vec![f_forward.clone()];
        let mut lf2: Vec<TopoShape> = Vec::new();
        self.find_same_domain(&mut lf1, &mut lf2);
        let n1 = lf1.len();
        let n2 = lf2.len();

        // cxx:1206-1214: with no same-domain partner on the other shape the
        // orientation is not reversed in `FillFace`.
        if n2 == 0 {
            rev_ori1 = false;
        }
        if n1 == 0 {
            rev_ori2 = false;
        }

        // cxx:1218
        let mut wes = WireEdgeSet::new(&f_forward);

        // cxx:1239-1252
        for f_cur in &lf1 {
            self.fill_face(f_cur, to_build1, &lf2, to_build2, &mut wes, rev_ori1);
        }
        for f_cur in &lf2 {
            self.fill_face(f_cur, to_build2, &lf1, to_build1, &mut wes, rev_ori2);
        }

        // cxx:1256
        self.add_intersection_edges(&f_forward, to_build1, rev_ori1, &mut wes);

        // cxx:1264-1268: `forceclass = false`.
        let mut fbu = FaceBuilder::with_wire_edge_set(&mut wes, &f_forward, false);

        // cxx:1272-1275
        let mut face_list = self.take_merged(&f_forward, to_build1);
        self.make_faces(&f_forward, &mut fbu, &mut face_list);
        self.put_merged(&f_forward, to_build1, face_list.clone());

        // cxx:1279-1297
        for f_cur in &lf1 {
            self.mark_split(f_cur, to_build1, true);
            let fl = self.change_split(f_cur, to_build1);
            if connect_to1 {
                *fl = face_list.clone();
            }
        }
        for f_cur in &lf2 {
            self.mark_split(f_cur, to_build2, true);
            let fl = self.change_split(f_cur, to_build2);
            if connect_to2 {
                *fl = face_list.clone();
            }
        }
    }

    /// `SplitFace2(Foriented, ToBuild1, ToBuild2)` (`cxx:1304-1531`).
    ///
    /// Splits the same-domain partners of both orientations into two face
    /// builders: `WES1` for the `<ToBuild1>` faces, `WES2` for `<ToBuild2>`.
    pub fn split_face2(&mut self, f_oriented: &TopoShape, to_build1: State, to_build2: State) {
        // cxx:1315-1318
        let tosplit = self.to_split(f_oriented, to_build1);
        if !tosplit {
            return;
        }

        // cxx:1322-1325
        let mut rev_ori1 = Self::reverse(to_build1, to_build2);
        let mut rev_ori2 = Self::reverse(to_build2, to_build1);
        let connect_to1 = true;
        let connect_to2 = false;

        // cxx:1328-1329
        let mut f_forward = f_oriented.clone();
        self.my_build_tool
            .set_orientation(&mut f_forward, Orientation::Forward);

        // cxx:1331-1337
        let mut lf1 = vec![f_forward.clone()];
        let mut lf2: Vec<TopoShape> = Vec::new();
        self.find_same_domain(&mut lf1, &mut lf2);
        let n1 = lf1.len();
        let n2 = lf2.len();

        // cxx:1352-1359
        if n2 == 0 {
            rev_ori1 = false;
        }
        if n1 == 0 {
            rev_ori2 = false;
        }

        // cxx:1361-1367
        let mut lfso = vec![f_forward.clone()];
        let mut lfoo: Vec<TopoShape> = Vec::new();
        self.find_same_domain_same_orientation(&mut lfso, &mut lfoo);

        // cxx:1372-1379: the partners of the operand <F> belongs to and the
        // partners of the other operand.
        let rank_f = self.shape_rank(f_oriented);
        let rank_x = if rank_f != 0 {
            if rank_f == 1 {
                2
            } else {
                1
            }
        } else {
            0
        };
        let mut lfso1: Vec<TopoShape> = Vec::new();
        let mut lfoo1: Vec<TopoShape> = Vec::new();
        let mut lfso2: Vec<TopoShape> = Vec::new();
        let mut lfoo2: Vec<TopoShape> = Vec::new();
        self.find_same_rank(&lfso, rank_f, &mut lfso1);
        self.find_same_rank(&lfoo, rank_f, &mut lfoo1);
        self.find_same_rank(&lfso, rank_x, &mut lfso2);
        self.find_same_rank(&lfoo, rank_x, &mut lfoo2);

        // cxx:1393-1396
        let tob1 = to_build1;
        let tob2 = to_build2;
        let tob1comp = if to_build1 == State::In {
            State::Out
        } else {
            State::In
        };
        let tob2comp = if to_build2 == State::In {
            State::Out
        } else {
            State::In
        };

        // ----- WES1 (`cxx:1400-1450`) -----
        let mut wes1 = WireEdgeSet::new(&f_forward);
        for f_cur in &lfso1 {
            self.fill_face(f_cur, tob1, &lf2, tob2, &mut wes1, rev_ori1);
        }
        for f_cur in &lfso2 {
            self.fill_face(f_cur, tob2, &lf1, tob1, &mut wes1, rev_ori2);
        }
        for f_cur in &lfoo1 {
            self.fill_face(f_cur, tob1comp, &lf2, to_build2, &mut wes1, !rev_ori1);
        }
        for f_cur in &lfoo2 {
            self.fill_face(f_cur, tob2comp, &lf1, to_build1, &mut wes1, !rev_ori2);
        }
        self.add_intersection_edges(&f_forward, to_build1, rev_ori1, &mut wes1);

        // cxx:1444-1449
        let mut fbu1 = FaceBuilder::with_wire_edge_set(&mut wes1, &f_forward, false);
        let mut face_list1 = self.take_merged(&f_forward, to_build1);
        self.make_faces(&f_forward, &mut fbu1, &mut face_list1);
        self.put_merged(&f_forward, to_build1, face_list1.clone());

        // cxx:1453-1462
        for f_cur in &lf1 {
            self.mark_split(f_cur, to_build1, true);
            let fl = self.change_split(f_cur, to_build1);
            if connect_to1 {
                *fl = face_list1.clone();
            }
        }

        // ----- WES2 (`cxx:1467-1515`) -----
        let mut wes2 = WireEdgeSet::new(&f_forward);
        for f_cur in &lfso1 {
            self.fill_face(f_cur, tob1comp, &lf2, tob2, &mut wes2, !rev_ori1);
        }
        for f_cur in &lfso2 {
            self.fill_face(f_cur, tob2comp, &lf1, tob1, &mut wes2, !rev_ori2);
        }
        for f_cur in &lfoo1 {
            self.fill_face(f_cur, tob1, &lf2, to_build2, &mut wes2, rev_ori1);
        }
        for f_cur in &lfoo2 {
            self.fill_face(f_cur, tob2, &lf1, to_build1, &mut wes2, rev_ori2);
        }
        self.add_intersection_edges(&f_forward, to_build2, rev_ori2, &mut wes2);

        // cxx:1510-1515
        let mut fbu2 = FaceBuilder::with_wire_edge_set(&mut wes2, &f_forward, false);
        let mut face_list2 = self.take_merged(&f_forward, to_build2);
        self.make_faces(&f_forward, &mut fbu2, &mut face_list2);
        self.put_merged(&f_forward, to_build2, face_list2.clone());

        // cxx:1520-1529
        for f_cur in &lf2 {
            self.mark_split(f_cur, to_build2, true);
            let fl = self.change_split(f_cur, to_build2);
            if connect_to2 {
                *fl = face_list2.clone();
            }
        }
    }

    /// `SplitSolid(S1oriented, ToBuild1, ToBuild2)` (`cxx:1535-1700`).
    pub fn split_solid(&mut self, s1_oriented: &TopoShape, to_build1: State, to_build2: State) {
        // cxx:1540-1562: a shell is split when any of its faces is.
        let is_shell = s1_oriented.shape_type() == ShapeType::Shell;
        let tosplit = if is_shell {
            explorer_shapes(s1_oriented, ShapeType::Face)
                .iter()
                .any(|sh| self.to_split(sh, to_build1))
        } else {
            self.to_split(s1_oriented, to_build1)
        };
        if !tosplit {
            return;
        }

        // cxx:1566-1569
        let mut rev_ori1 = Self::reverse(to_build1, to_build2);
        let mut rev_ori2 = Self::reverse(to_build2, to_build1);
        let connect_to1 = true;
        let connect_to2 = false;

        // cxx:1572-1574
        let mut s1_forward = s1_oriented.clone();
        self.my_build_tool
            .set_orientation(&mut s1_forward, Orientation::Forward);

        // cxx:1576-1580
        let mut ls1 = vec![s1_forward.clone()];
        let mut ls2: Vec<TopoShape> = Vec::new();
        self.find_same_domain(&mut ls1, &mut ls2);
        let n1 = ls1.len();
        let n2 = ls2.len();

        // cxx:1582-1591
        if n2 == 0 {
            rev_ori1 = false;
        }
        if n1 == 0 {
            rev_ori2 = false;
        }

        // cxx:1595
        let mut sfs = ShellFaceSet::new();

        // cxx:1609-1622
        STATIC_SOLID_INDEX.store(1, Ordering::SeqCst);
        for s_cur in &ls1 {
            self.fill_solid(s_cur, to_build1, &ls2, to_build2, &mut sfs, rev_ori1);
        }
        STATIC_SOLID_INDEX.store(2, Ordering::SeqCst);
        for s_cur in &ls2 {
            self.fill_solid(s_cur, to_build2, &ls1, to_build1, &mut sfs, rev_ori2);
        }

        // cxx:1626-1653: the new faces of the intersection surfaces.
        if let Some(hds) = self.data_structure() {
            if hds.nb_surfaces() > 0 {
                let mut ssurfaces = hds.solid_surfaces_of_shape(&s1_forward);
                while ssurfaces.more() {
                    let is = ssurfaces.geometry();
                    let ori = ssurfaces.orientation(to_build1);
                    let lnew_f: Vec<TopoShape> =
                        self.new_faces(is).map(|l| l.to_vec()).unwrap_or_default();
                    for a_face in lnew_f {
                        let mut a_face = a_face;
                        self.my_build_tool.set_orientation(&mut a_face, ori);
                        sfs.add_element(&a_face);
                    }
                    ssurfaces.next();
                }
            }
        }

        // cxx:1657-1672
        let mut sobu = SolidBuilder::with_shell_face_set(&mut sfs, false);
        let mut solid_list = self.take_merged(s1_oriented, to_build1);
        if is_shell {
            self.make_shells(&mut sobu, &mut solid_list);
        } else {
            self.make_solids(&mut sobu, &mut solid_list);
        }
        self.put_merged(s1_oriented, to_build1, solid_list.clone());

        // cxx:1676-1698
        for s_cur in &ls1 {
            self.mark_split(s_cur, to_build1, true);
            let sl = self.change_split(s_cur, to_build1);
            if connect_to1 {
                *sl = solid_list.clone();
            }
        }
        for s_cur in &ls2 {
            self.mark_split(s_cur, to_build2, true);
            let sl = self.change_split(s_cur, to_build2);
            if connect_to2 {
                *sl = solid_list.clone();
            }
        }
    }

    /// `SplitShapes(Ex, ToBuild1, ToBuild2, aSet, RevOri)` (`cxx:1713-1872`).
    pub fn split_shapes(
        &mut self,
        ex: &mut ShapeExplorer,
        to_build1: State,
        to_build2: State,
        a_set: &mut dyn ShapeSetSource,
        rev_ori: bool,
    ) {
        while ex.more() {
            let Some(a_shape) = ex.current().cloned() else {
                break;
            };

            // cxx:1727: the orientation the new shapes must carry.
            let newori = Self::orient(self.my_build_tool.orientation(&a_shape), rev_ori);

            let t = a_shape.shape_type();
            match t {
                ShapeType::Solid | ShapeType::Shell => {
                    self.split_solid(&a_shape, to_build1, to_build2)
                }
                ShapeType::Face => self.split_face(&a_shape, to_build1, to_build2),
                ShapeType::Edge => self.split_edge(&a_shape, to_build1, to_build2),
                // cxx:1743-1746
                _ => {
                    ex.next();
                    continue;
                }
            }

            if self.is_split(&a_shape, to_build1) {
                // cxx:1749-1794
                let ls: Vec<TopoShape> = self.splits(&a_shape, to_build1).to_vec();
                // A split edge of <ToBuild1> with no part of its own uses the
                // `<ON>` parts, filtered by the 2d-only `myListOfFace`.
                let is_lson = t == ShapeType::Edge && to_build1 == State::In && ls.is_empty();
                let parts: Vec<TopoShape> = if is_lson {
                    self.splits(&a_shape, State::On).to_vec()
                } else {
                    ls
                };
                for new_shape in parts {
                    let mut new_shape = new_shape;
                    self.my_build_tool.set_orientation(&mut new_shape, newori);
                    if is_lson {
                        // cxx:1777-1787 (2d pur)
                        let lof = self.my_list_of_face.clone();
                        let add = if lof.is_empty() {
                            true
                        } else {
                            self.keep_shape(&new_shape, &lof, to_build1)
                        };
                        if add {
                            a_set.add_start_element(&new_shape);
                        }
                    } else {
                        a_set.add_start_element(&new_shape);
                    }
                }
            } else {
                // cxx:1796-1868: no split part for <ToBuild1>, so the
                // <ToBuild1> parts of <aShape> are kept when it survives.
                let mut add = true;
                let isedge = t == ShapeType::Edge;
                let mut hs = false;
                let mut hg = false;
                if let Some(hds) = self.data_structure() {
                    hs = hds.has_shape(&a_shape, true);
                    hg = hds.has_geometry(&a_shape);
                }

                let mut testkeep = isedge && hs && !hg;

                // cxx:1806-1813: an edge absent from the DS but whose bounding
                // vertex is on the other shape is still classified
                // (`USA60299`).
                let mut istouched = isedge && !hs && !hg;
                if istouched {
                    if let (Some(edge), Some(hds)) =
                        (Edge::wrap(a_shape.clone()), self.data_structure())
                    {
                        istouched = fun_touched(hds.ds(), &edge);
                    }
                }
                testkeep = testkeep || istouched;

                if testkeep {
                    if !self.my_list_of_face.is_empty() {
                        // cxx:1818-1822 (2d pur)
                        let lof = self.my_list_of_face.clone();
                        add = self.keep_shape(&a_shape, &lof, to_build1);
                    } else {
                        // cxx:1826-1861 (3d): classify a fixed probe point of
                        // the edge against the other operand solid.
                        let sol = if STATIC_SOLID_INDEX.load(Ordering::SeqCst) == 1 {
                            self.my_shape2.clone()
                        } else {
                            self.my_shape1.clone()
                        };
                        match sol {
                            Some(sol) => {
                                let edge = Edge::wrap(a_shape.clone())
                                    .expect("SplitShapes: not an edge");
                                let c3d = BRepTool::edge_curve(&edge)
                                    .expect("SplitShapes no 3D curve on edge");
                                let (first, last) = BRepTool::edge_parameters(&edge);
                                // `tt` is OCCT's probe parameter (`cxx:1843`);
                                // it is not derived from the edge.
                                let tt = 0.127956477;
                                let par = (1.0 - tt) * first + tt * last;
                                let p3d = c3d.d0(par);
                                let solid =
                                    Solid::wrap(sol).expect("SplitShapes: not a solid");
                                let mut scl = SolidClassifier::new();
                                let state =
                                    scl.classify_solid(&solid, &p3d, CONFUSION);
                                add = state == to_build1;
                            }
                            // cxx:1858-1861
                            None => add = true,
                        }
                    }
                }

                if add {
                    let mut shape = a_shape.clone();
                    self.my_build_tool.set_orientation(&mut shape, newori);
                    a_set.add_element(&shape);
                }
            }

            ex.next();
        }
    }

    /// `AddIntersectionEdges(aFace, ToBuild1, RevOri1, WES)` (`cxx:150-179`).
    ///
    /// Every new edge of each curve of the face is oriented for `<ToBuild1>`
    /// and gets its pcurve on `aFace` before entering the set.
    fn add_intersection_edges(
        &self,
        a_face: &TopoShape,
        to_build1: State,
        rev_ori1: bool,
        wes: &mut WireEdgeSet,
    ) {
        let Some(hds) = self.data_structure() else {
            return;
        };
        let mut fcurves = hds.face_curves_of_shape(a_face);
        while fcurves.more() {
            let ic = fcurves.geometry();
            let lnew_e: Vec<TopoShape> = self.new_edges(ic).to_vec();
            for an_edge in lnew_e {
                let ori = fcurves.orientation(to_build1);
                let newori = Self::orient(ori, rev_ori1);
                if newori == Orientation::External {
                    continue;
                }
                let mut an_edge = an_edge;
                self.my_build_tool.set_orientation(&mut an_edge, newori);
                if let Some(pc) = fcurves.p_curve() {
                    self.my_build_tool.p_curve(a_face, &an_edge, &pc);
                }
                wes.add_start_element(&an_edge);
            }
            fcurves.next();
        }
    }

    /// `ChangeMerged(S, ToBuild)` with the list detached, so a `Make*` call can
    /// push into it while borrowing `self` mutably.
    ///
    /// OCCT hands the list by reference and `Make*` appends; detaching and
    /// putting it back keeps exactly that append behaviour.
    fn take_merged(&mut self, s: &TopoShape, to_build: State) -> Vec<TopoShape> {
        std::mem::take(self.change_merged(s, to_build))
    }

    /// Puts back a [`Builder::take_merged`] list.
    fn put_merged(&mut self, s: &TopoShape, to_build: State, l: Vec<TopoShape>) {
        *self.change_merged(s, to_build) = l;
    }
}

/// `FUN_touched(BDS, EOR)` (`cxx:1702-1709`): true when one of the two bounding
/// vertices of `EOR` is a shape of the DS.
fn fun_touched(bds: &DataStructure, eor: &Edge) -> bool {
    let (vf, vl) = edge_vertices(eor);
    let hvf = vf.map(|v| bds.has_shape(&v.0, true)).unwrap_or(false);
    let hvl = vl.map(|v| bds.has_shape(&v.0, true)).unwrap_or(false);
    hvf || hvl
}
