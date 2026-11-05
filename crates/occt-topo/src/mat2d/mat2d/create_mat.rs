//! `MAT2d_Mat2d::CreateMatOpen` (`MAT2d_Mat2d.cxx:122-859`) and
//! `MAT2d_Mat2d::CreateMat` (`cxx:852-1628`).
//!
//! Both functions are copies of one main loop (they are 738 and 777 lines of
//! near-identical C++). The port keeps one body and spells the differences out
//! where OCCT's two copies differ, with the line numbers of both:
//!
//! | difference | open | closed |
//! |---|---|---|
//! | initial `noofbisectors` | `noofedges - 1` (`:176`) | `noofedges` (`:913`) |
//! | initial bisector loop bound | `Number() - 1` (`:201`) | `Number()` (`:939`) |
//! | edge / bisector assignment | chain (`:221-233`) | ring (`:957-963`) |
//! | scanned edges | second to last (`:303`) | all (`:1074`) |
//! | anti-loop patch | absent | `:996-1033`, `:1396-1401` |
//! | semi-infinite root loop | `Number() - 1` (`:826`) | `Number()` (`:1599`) |
//!
//! Everything else - the removal loop, the step 2.5 parameter analysis, the
//! step 2.5 edge elimination and the step 3 root collection - is line for line
//! the same in both copies.

use occt_core::precision::INFINITE;

use super::{classify, Mat2dMat2d};
use crate::mat::{MatBisector, MatEdge, MatTList};
use crate::mat2d::Mat2dTool2d;

impl Mat2dMat2d {
    /// The shared body of `CreateMat` / `CreateMatOpen`.
    ///
    /// `noofbisectorstoremove` and `beginbisector` keep their OCCT
    /// initialisation (`cxx:134`, `:138`, `:177`) even though the first main
    /// loop iteration overwrites them; the algorithm relies on them being
    /// defined when the main loop is left through `interrupt` on its first pass.
    #[allow(unused_assignments)]
    pub(super) fn create(&mut self, tool: &mut Mat2dTool2d, closed: bool) {
        let mut interrupt = false;

        let mut noofbisectorstoremove: usize = 0;
        let mut beginbisector: i32 = 0;
        let mut noofbisectors: i32;

        // `cxx:140-146` / `:873-879`.
        let mut nb_iter_bis: i32 = 0;
        let even_nb_iter_bis: i32 = 10;
        let mut edge_numbers = vec![-1i32; (even_nb_iter_bis + 2) as usize];
        let mut to_nullify_noofbisectorstoremove = false;

        let noofedges = tool.number_of_items();
        let toleranceofconfusion = tool.tolerance_of_confusion();
        let number_max_of_ite = noofedges * noofedges;

        // `NCollection_Array1<int> firstarea(0, noofedges)` and friends.
        let mut firstarea = vec![0i32; (noofedges + 1) as usize];
        let mut lastarea = vec![0i32; (noofedges + 1) as usize];
        let mut noofarea = vec![0i32; (noofedges + 1) as usize];

        let mut parama = [0i32; 2];
        let mut paramb = [0i32; 2];

        // `cxx:900-904`: state of the closed-only anti-loop patch.
        let mut a_nb_of_narea1 = 0;
        let mut a_pref_narea = 0;
        let a_nb_max_narea1 = 10;
        let mut a_nb_elts = [0i32; 2];
        let mut a_count_elts = [0i32; 2];
        let mut is_break = false;

        // `cxx:170-177` / `:905-911`.
        self.bisector_one_to_remove.clear();
        self.bisector_two_to_remove.clear();
        self.type_of_bisector_to_remove.clear();
        self.bisector_map.clear();

        self.is_done = true;
        noofbisectors = if closed { noofedges } else { noofedges - 1 };
        beginbisector = 0;

        // `cxx:182-195` / `:920-931`: one `MAT_Edge` per circuit item.
        self.edge_list = MatTList::new();
        self.removed_edges.clear();
        for i in 0..noofedges {
            let mut edge = MatEdge::new();
            edge.set_edge_number(i + 1);
            edge.set_distance(-1.0);
            let id = self.graph.add_edge(edge);
            self.edge_list.back_add(id);
        }
        self.edge_list.loop_list();

        // `cxx:198-215` / `:935-952`: the bisectors issued from the contour.
        self.edge_list.first();
        let edge_count = self.edge_list.number() as i32;
        let bisector_count = if closed { edge_count } else { edge_count - 1 };
        for i in 0..bisector_count {
            let mut bisector = MatBisector::new();
            bisector.set_index_number(i);

            let first_edge = *self.edge_list.current();
            bisector.set_first_edge(Some(first_edge));
            let first_number = self.edg(first_edge).edge_number();
            let first_vector = tool.tangent_before(first_number, self.my_is_open_result);
            bisector.set_first_vector(first_vector);

            self.edge_list.next();

            let second_edge = *self.edge_list.current();
            bisector.set_second_edge(Some(second_edge));
            let second_number = self.edg(second_edge).edge_number();
            let mut dist = 0.0;
            let issue_point = tool.first_point(second_number, &mut dist);
            bisector.set_issue_point(issue_point);
            bisector.set_dist_issue_point(dist);
            let second_vector = tool.tangent_after(second_number, self.my_is_open_result);
            bisector.set_second_vector(second_vector);

            let id = self.graph.add_bisector(bisector);
            self.bisector_map.insert(i, id);
        }

        // `cxx:217-233` (open) / `:954-963` (closed): the two bisectors of
        // every edge. The open chain gives both ends the same bisector; the
        // closed ring shifts the first one by one.
        self.edge_list.first();
        if closed {
            for i in 0..edge_count {
                let first = self.map((i - 1 + noofbisectors).rem_euclid(noofbisectors));
                let second = self.map(i);
                let current = *self.edge_list.current();
                self.edg_mut(current).set_first_bisector(Some(first));
                self.edg_mut(current).set_second_bisector(Some(second));
                self.edge_list.next();
            }
        } else {
            let first = self.map(0);
            let current = *self.edge_list.current();
            self.edg_mut(current).set_first_bisector(Some(first));
            self.edg_mut(current).set_second_bisector(Some(first));
            self.edge_list.next();
            for i in 1..edge_count - 1 {
                let previous = self.map(i - 1);
                let current_bisector = self.map(i);
                let current = *self.edge_list.current();
                self.edg_mut(current).set_first_bisector(Some(previous));
                self.edg_mut(current).set_second_bisector(Some(current_bisector));
                self.edge_list.next();
            }
            let last = self.map(edge_count - 2);
            let current = *self.edge_list.current();
            self.edg_mut(current).set_first_bisector(Some(last));
            self.edg_mut(current).set_second_bisector(Some(last));
        }

        // step 2: the main loop (`cxx:238` / `:969`).
        let mut number_of_ite = 0;
        while (self.edge_list.number() as i32) > 1 {
            // step 2.1: `cxx:246-259` / `:977-990`.
            let a_nb_bis = noofbisectors - beginbisector;
            for i in beginbisector..noofbisectors {
                let id = self.map(i);
                let edge1 = self.first_edge_number(id);
                let edge2 = self.second_edge_number(id);
                let mut bisector = self.bis(id).clone();
                tool.create_bisector_of(&mut bisector, edge1, edge2);
                *self.bis_mut(id) = bisector;
                self.number_of_bisectors += 1;
            }

            // `cxx:996-1033`: the closed variant's patch against an infinite
            // loop caused by bad geometry.
            if closed {
                if a_nb_bis == 1 {
                    if a_pref_narea == 1 {
                        a_nb_of_narea1 += 1;
                        let start = self.map(beginbisector);
                        let edge1number = self.first_edge_number(start);
                        let edge2number = self.second_edge_number(start);
                        if a_nb_elts[0] == edge1number {
                            a_count_elts[0] += 1;
                        } else {
                            a_count_elts[0] = 0;
                            a_nb_elts[0] = edge1number;
                        }
                        if a_nb_elts[1] == edge2number {
                            a_count_elts[1] += 1;
                        } else {
                            a_count_elts[1] = 0;
                            a_nb_elts[1] = edge2number;
                        }
                        if a_nb_of_narea1 >= a_nb_max_narea1
                            && (a_count_elts[0] >= a_nb_max_narea1
                                || a_count_elts[1] >= a_nb_max_narea1)
                        {
                            is_break = true;
                        }
                    } else {
                        a_nb_of_narea1 = 0;
                        a_count_elts[0] = 0;
                        a_count_elts[1] = 0;
                    }
                }
                a_pref_narea = a_nb_bis;
            }

            // step 2 exit test: `cxx:264-270` / `:1037-1042`.
            if (self.edge_list.number() as i32) < 3 {
                break;
            }

            // loop 2: while there are bisectors to remove (`cxx:275` / `:1046`).
            loop {
                nb_iter_bis += 1;
                noofbisectorstoremove = 0;

                // `cxx:280-282` / `:1051-1052`.
                self.edge_list.first();
                if !closed {
                    self.edge_list.next();
                }

                // `cxx:287-301` / `:1059-1072`.
                if nb_iter_bis <= even_nb_iter_bis + 1 {
                    edge_numbers[nb_iter_bis as usize] = self.edge_list.number() as i32;
                } else {
                    for k in 1..=even_nb_iter_bis {
                        edge_numbers[k as usize] = edge_numbers[(k + 1) as usize];
                    }
                    edge_numbers[(even_nb_iter_bis + 1) as usize] =
                        self.edge_list.number() as i32;
                }
                if edge_numbers[(even_nb_iter_bis + 1) as usize] == edge_numbers[1] {
                    to_nullify_noofbisectorstoremove = true;
                }

                // step 2.2: `cxx:303-337` / `:1074-1106`.
                let scan_count = if closed {
                    self.edge_list.number() as i32
                } else {
                    self.edge_list.number() as i32 - 2
                };
                for _ in 0..scan_count {
                    self.intersect_current_edge(tool, beginbisector, &mut noofbisectorstoremove);
                    self.edge_list.next();
                }

                // loop 2 exit test: `cxx:340-347` / `:1110-1118`.
                if to_nullify_noofbisectorstoremove {
                    noofbisectorstoremove = 0;
                }
                if noofbisectorstoremove == 0 {
                    break;
                }

                // step 2.4: `cxx:352-458` / `:1122-1228`.
                for i in 0..noofbisectorstoremove {
                    self.annihilate(tool, i);
                }
            }

            // step 2.5: `cxx:469-607` / `:1238-1375`. The parameters of the
            // intersections decide which runs of edges leave the contour.
            self.edge_list.first();
            if !closed {
                self.edge_list.next();
            }

            let current = *self.edge_list.current();
            let firstbisector = self
                .edg(current)
                .first_bisector()
                .expect("MAT2d_Mat2d::CreateMat: edge without FirstBisector");
            let parameters = classify(self.bis(firstbisector), tool, toleranceofconfusion);
            parama[0] = parameters.0;
            paramb[0] = parameters.1;

            let mut narea: i32 = -1;
            let scan_count = if closed {
                self.edge_list.number() as i32
            } else {
                self.edge_list.number() as i32 - 2
            };
            for _ in 0..scan_count {
                let current_index = self.edge_list.index() as i32;
                let current = *self.edge_list.current();
                let secondbisector = self
                    .edg(current)
                    .second_bisector()
                    .expect("MAT2d_Mat2d::CreateMat: edge without SecondBisector");
                let parameters = classify(self.bis(secondbisector), tool, toleranceofconfusion);
                parama[1] = parameters.0;
                paramb[1] = parameters.1;

                // `cxx:575-601` / `:1345-1371`.
                if paramb[0] > 0 && parama[1] > 0 {
                    if narea < 0 {
                        narea += 1;
                        firstarea[narea as usize] = current_index;
                        lastarea[narea as usize] = firstarea[narea as usize];
                        noofarea[narea as usize] = 1;
                    } else if current_index == lastarea[narea as usize] + 1 {
                        lastarea[narea as usize] += 1;
                        noofarea[narea as usize] += 1;
                    } else {
                        narea += 1;
                        firstarea[narea as usize] = current_index;
                        lastarea[narea as usize] = firstarea[narea as usize];
                        noofarea[narea as usize] = 1;
                    }
                }
                parama[0] = parama[1];
                paramb[0] = paramb[1];
                self.edge_list.next();
            }

            // `cxx:608-617` / `:1377-1386`: a run that wraps around the ring is
            // merged into the first one.
            let mut compact = 0;
            if narea > 0
                && lastarea[narea as usize] == self.edge_list.number() as i32
                && firstarea[0] == 1
            {
                firstarea[0] = firstarea[narea as usize];
                noofarea[0] += noofarea[narea as usize];
                compact = noofarea[narea as usize];
                narea -= 1;
            }
            narea += 1;

            // `cxx:1395-1401`: the closed variant's anti-loop exit.
            if closed && narea == 1 && is_break {
                narea = 0;
            }

            // step 2.6: `cxx:623-629` / `:1403-1407`.
            if narea == 0 {
                interrupt = true;
                break;
            }

            // step 2.5 (continued): `cxx:636-728` / `:1409-1503`, elimination of
            // the runs and creation of the bisectors between the new
            // neighbours.
            beginbisector = noofbisectors;
            let mut shift = 0;
            let mut all = 0;
            if narea == 1 && noofarea[0] == self.edge_list.number() as i32 {
                all = 1;
            }

            for i in 0..narea {
                if i == 1 {
                    shift -= compact;
                }
                self.edge_list.first();
                if !closed {
                    self.edge_list.next();
                }
                let bracket = (firstarea[i as usize] - shift) as usize;
                let mut edgetoremove = *self.edge_list.brackets(bracket);

                // `cxx:657-669` / `:1431-1443`.
                let first_bisector = self
                    .edg(edgetoremove)
                    .first_bisector()
                    .expect("MAT2d_Mat2d::CreateMat: edge without FirstBisector");
                let intersection_point = self.edg(edgetoremove).intersection_point();
                self.bis_mut(first_bisector)
                    .set_end_point(intersection_point);
                let second_parameter = self.bis(first_bisector).second_parameter();
                self.bis_mut(first_bisector)
                    .set_first_parameter(second_parameter);
                let bisector = self.bis(first_bisector).clone();
                tool.trim_bisector(&bisector);

                // `cxx:671-677` / `:1445-1451`.
                let mut newbisector = MatBisector::new();
                newbisector.set_index_number(noofbisectors);
                newbisector.set_dist_issue_point(self.edg(edgetoremove).distance());
                newbisector.set_issue_point(intersection_point);
                newbisector.set_first_edge(Some(*self.edge_list.previous_item()));
                newbisector.add_bisector(first_bisector);
                let new_id = self.graph.add_bisector(newbisector);
                self.bisector_map.insert(noofbisectors, new_id);

                for j in 0..noofarea[i as usize] {
                    let current = *self.edge_list.current();
                    self.removed_edges.push(current);
                    self.edge_list.unlink();
                    self.edge_list.next();
                    shift += 1;

                    // `cxx:689-692` / `:1463-1466`.
                    if all == 0 || j + 1 != noofarea[i as usize] {
                        let second_bisector = self
                            .edg(edgetoremove)
                            .second_bisector()
                            .expect("MAT2d_Mat2d::CreateMat: edge without SecondBisector");
                        self.bis_mut(new_id).add_bisector(second_bisector);
                    }
                    let second_bisector = self
                        .edg(edgetoremove)
                        .second_bisector()
                        .expect("MAT2d_Mat2d::CreateMat: edge without SecondBisector");
                    let edge_point = self.edg(edgetoremove).intersection_point();
                    self.bis_mut(second_bisector).set_end_point(edge_point);
                    let first_parameter = self.bis(second_bisector).first_parameter();
                    self.bis_mut(second_bisector)
                        .set_second_parameter(first_parameter);
                    let bisector = self.bis(second_bisector).clone();
                    tool.trim_bisector(&bisector);

                    edgetoremove = *self.edge_list.current();
                }

                // `cxx:709-727` / `:1483-1501`.
                let current = *self.edge_list.current();
                self.bis_mut(new_id).set_second_edge(Some(current));
                let previous = *self.edge_list.previous_item();
                self.edg_mut(previous).set_second_bisector(Some(new_id));
                self.edg_mut(current).set_first_bisector(Some(new_id));

                let first_bisector = self
                    .bis(new_id)
                    .first_bisector()
                    .expect("MAT2d_Mat2d::CreateMat: new bisector has no FirstBisector");
                let number = self.bis(first_bisector).bisector_number();
                let first_vector = tool.tangent(number);
                self.bis_mut(new_id).set_first_vector(first_vector);

                let last_bisector = self
                    .bis(new_id)
                    .last_bisector()
                    .expect("MAT2d_Mat2d::CreateMat: new bisector has no LastBisector");
                let number = self.bis(last_bisector).bisector_number();
                let second_vector = tool.tangent(number);
                self.bis_mut(new_id).set_second_vector(second_vector);

                noofbisectors += 1;

                let previous = *self.edge_list.previous_item();
                self.edg_mut(previous).set_distance(-1.0);
                let current = *self.edge_list.current();
                self.edg_mut(current).set_distance(-1.0);

                let first_bisector = self
                    .edg(previous)
                    .first_bisector()
                    .expect("MAT2d_Mat2d::CreateMat: edge without FirstBisector");
                self.bis_mut(first_bisector).set_second_parameter(INFINITE);
                let second_bisector = self
                    .edg(current)
                    .second_bisector()
                    .expect("MAT2d_Mat2d::CreateMat: edge without SecondBisector");
                self.bis_mut(second_bisector)
                    .set_first_parameter(INFINITE);
            }

            // `cxx:735-741` / `:1507-1513`.
            if number_of_ite > number_max_of_ite {
                self.is_done = false;
                break;
            }
            number_of_ite += 1;
        }

        // step 3: `cxx:751-858` / `:1525-1627`.
        self.finish_create(tool, closed, interrupt, beginbisector, noofbisectors);
    }

    /// `MAT2d_Tool2d::IntersectBisector` on the current edge, and the three
    /// `Intersect` calls it triggers (`cxx:307-333` / `:1078-1104`).
    fn intersect_current_edge(
        &mut self,
        tool: &mut Mat2dTool2d,
        beginbisector: i32,
        noofbisectorstoremove: &mut usize,
    ) {
        let edge = *self.edge_list.current();
        if self.edg(edge).distance() != -1.0 {
            return;
        }
        let firstbisector = self
            .edg(edge)
            .first_bisector()
            .expect("MAT2d_Mat2d::CreateMat: edge without FirstBisector");
        let secondbisector = self
            .edg(edge)
            .second_bisector()
            .expect("MAT2d_Mat2d::CreateMat: edge without SecondBisector");

        let edges_one = self.edge_numbers(firstbisector);
        let edges_two = self.edge_numbers(secondbisector);
        let mut bisector_one = self.bis(firstbisector).clone();
        let mut bisector_two = self.bis(secondbisector).clone();
        let (distance, int_pnt) =
            tool.intersect_bisector(&mut bisector_one, edges_one, &mut bisector_two, edges_two);
        *self.bis_mut(firstbisector) = bisector_one;
        *self.bis_mut(secondbisector) = bisector_two;

        self.edg_mut(edge).set_distance(distance);
        // `IntPnt` is only assigned when the intersection is retained
        // (`MAT2d_Tool2d.cxx:1006`); OCCT then stores an uninitialised value,
        // the port leaves the edge's point untouched.
        if let Some(point) = int_pnt {
            self.edg_mut(edge).set_intersection_point(point);
        }

        if distance == INFINITE {
            if self.bis(firstbisector).index_number() >= beginbisector
                || self.bis(secondbisector).index_number() >= beginbisector
            {
                self.intersect(tool, 0, noofbisectorstoremove, firstbisector, secondbisector);
            }
        } else {
            if self.bis(firstbisector).index_number() >= beginbisector {
                self.intersect(tool, 1, noofbisectorstoremove, firstbisector, secondbisector);
            }
            if self.bis(secondbisector).index_number() >= beginbisector {
                self.intersect(tool, 2, noofbisectorstoremove, firstbisector, secondbisector);
            }
        }
    }

    /// step 2.4, one entry of the `...toremove` maps (`cxx:352-458` /
    /// `:1122-1228`): the sub-bisectors between the two recorded bisectors are
    /// nullified and their edges reinserted into the contour.
    fn annihilate(&mut self, tool: &mut Mat2dTool2d, i: usize) {
        let mut bisectortoremove = *self
            .bisector_one_to_remove
            .get(&i)
            .expect("MAT2d_Mat2d::CreateMat: bisectoronetoremove");
        let last_to_remove = *self
            .bisector_two_to_remove
            .get(&i)
            .expect("MAT2d_Mat2d::CreateMat: bisectortwotoremove");
        let tree_side = *self
            .type_of_bisector_to_remove
            .get(&i)
            .expect("MAT2d_Mat2d::CreateMat: typeofbisectortoremove");

        loop {
            self.number_of_bisectors -= 1;

            // `cxx:375-386` / `:1146-1155`.
            let mut currentbisectorlist = self.bis(bisectortoremove).list().clone();
            currentbisectorlist.first();
            let currentbisector = *currentbisectorlist.current();
            let mut previousedge = self
                .bis(currentbisector)
                .first_edge()
                .expect("MAT2d_Mat2d::CreateMat: sub-bisector without FirstEdge");
            self.edge_list.init(&previousedge);
            self.edg_mut(previousedge).set_distance(-1.0);
            let first_bisector = self
                .edg(previousedge)
                .first_bisector()
                .expect("MAT2d_Mat2d::CreateMat: edge without FirstBisector");
            self.bis_mut(first_bisector).set_second_parameter(INFINITE);
            let second_bisector = self
                .edg(previousedge)
                .second_bisector()
                .expect("MAT2d_Mat2d::CreateMat: edge without SecondBisector");
            self.bis_mut(second_bisector).set_first_parameter(INFINITE);

            // `cxx:388-424` / `:1159-1196`.
            while currentbisectorlist.more() {
                let currentbisector = *currentbisectorlist.current();
                let currentedge = self
                    .bis(currentbisector)
                    .second_edge()
                    .expect("MAT2d_Mat2d::CreateMat: sub-bisector without SecondEdge");
                self.edge_list.link_after(currentedge);
                self.edge_list.next();
                self.edg_mut(currentedge)
                    .set_first_bisector(Some(currentbisector));
                self.edg_mut(previousedge)
                    .set_second_bisector(Some(currentbisector));

                self.bis_mut(currentbisector).set_first_parameter(INFINITE);
                self.bis_mut(currentbisector).set_second_parameter(INFINITE);
                let bisector = self.bis(currentbisector).clone();
                tool.trim_bisector(&bisector);

                self.edg_mut(currentedge).set_distance(-1.0);
                let first_bisector = self
                    .edg(currentedge)
                    .first_bisector()
                    .expect("MAT2d_Mat2d::CreateMat: edge without FirstBisector");
                self.bis_mut(first_bisector).set_second_parameter(INFINITE);
                let second_bisector = self
                    .edg(currentedge)
                    .second_bisector()
                    .expect("MAT2d_Mat2d::CreateMat: edge without SecondBisector");
                self.bis_mut(second_bisector).set_first_parameter(INFINITE);

                previousedge = currentedge;
                currentbisectorlist.next();
            }

            // `cxx:426-428` / `:1197-1198`.
            let current = *self.edge_list.current();
            self.removed_edges.push(current);
            self.edge_list.unlink();

            // `cxx:430-436` / `:1200-1206`.
            if self.bis(bisectortoremove).bisector_number()
                == self.bis(last_to_remove).bisector_number()
            {
                break;
            }

            // `cxx:438-449` / `:1208-1220`.
            bisectortoremove = if tree_side == 1 {
                self.bis(bisectortoremove)
                    .first_bisector()
                    .expect("MAT2d_Mat2d::CreateMat: empty bisector list")
            } else {
                self.bis(bisectortoremove)
                    .last_bisector()
                    .expect("MAT2d_Mat2d::CreateMat: empty bisector list")
            };
        }
    }
}
