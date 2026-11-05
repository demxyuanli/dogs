//! step 3 of `MAT2d_Mat2d::CreateMatOpen` (`MAT2d_Mat2d.cxx:749-858`) and of
//! `MAT2d_Mat2d::CreateMat` (`cxx:1523-1627`): the two remaining bisectors are
//! merged when the contour collapsed to two edges, the semi-infinite roots are
//! collected, and the roots list is handed out.

use super::Mat2dMat2d;
use crate::mat::MatTList;
use crate::mat2d::Mat2dTool2d;

impl Mat2dMat2d {
    /// `cxx:1525-1627` (closed) / `:751-858` (open).
    ///
    /// `beginbisector` and `noofbisectors` are the values the main loop left
    /// behind; `interrupt` is the `interrupt` flag of step 2.6.
    pub(super) fn finish_create(
        &mut self,
        tool: &mut Mat2dTool2d,
        closed: bool,
        interrupt: bool,
        beginbisector: i32,
        noofbisectors: i32,
    ) {
        let mut beginbisector = beginbisector;
        let mut noofbisectors = noofbisectors;

        if interrupt {
            // `cxx:756-758` / `:1530-1532`.
            self.semi_infinite = true;
        } else {
            self.semi_infinite = false;

            // `cxx:768-817` / `:1548-1594`: the two last bisectors separate the
            // same edges, so they are either merged or both semi-infinite. The
            // branch is unreachable since the two-edge case leaves the main
            // loop (OCCT says so in the comment at `:773`); it is ported as is.
            if (self.edge_list.number() as i32) > 1 {
                self.edge_list.first();
                let edge = *self.edge_list.current();
                let first_bisector = self
                    .edg(edge)
                    .first_bisector()
                    .expect("MAT2d_Mat2d::CreateMat: edge without FirstBisector");
                let second_bisector = self
                    .edg(edge)
                    .second_bisector()
                    .expect("MAT2d_Mat2d::CreateMat: edge without SecondBisector");

                if self.bis(first_bisector).index_number() == noofbisectors - 1 {
                    // `cxx:780-792`.
                    let issue_point = self.bis(first_bisector).issue_point();
                    let bisector = self.bis(second_bisector).clone();
                    if tool.trim_bisector_by_point(&bisector, issue_point) {
                        if self.bis(second_bisector).end_point() == 0 {
                            self.bis_mut(second_bisector).set_end_point(issue_point);
                        }
                        let root = self.map(noofbisectors - 1);
                        self.bis_mut(root).add_bisector(second_bisector);
                    } else {
                        self.semi_infinite = true;
                    }
                } else {
                    // `cxx:794-808`.
                    let issue_point = self.bis(second_bisector).issue_point();
                    let bisector = self.bis(first_bisector).clone();
                    if tool.trim_bisector_by_point(&bisector, issue_point) {
                        if self.bis(first_bisector).end_point() == 0 {
                            self.bis_mut(first_bisector).set_end_point(issue_point);
                        }
                        let root = self.map(noofbisectors - 1);
                        self.bis_mut(root).add_bisector(first_bisector);
                    } else {
                        self.semi_infinite = true;
                    }
                }

                // `cxx:810-816`.
                if !self.semi_infinite {
                    self.number_of_bisectors -= 1;
                    let root = self.map(noofbisectors - 1);
                    let current = *self.edge_list.current();
                    self.bis_mut(root).set_second_edge(Some(current));
                    self.bis_mut(root).set_bisector_number(-1);
                }
            }
        }

        if self.semi_infinite {
            // `cxx:820-831` / `:1596-1607`: every remaining edge contributes its
            // second bisector as a root. The open variant stops one short.
            beginbisector = noofbisectors;
            self.edge_list.first();
            let count = if closed {
                self.edge_list.number() as i32
            } else {
                self.edge_list.number() as i32 - 1
            };
            for _ in 0..count {
                let edge = *self.edge_list.current();
                let second_bisector = self
                    .edg(edge)
                    .second_bisector()
                    .expect("MAT2d_Mat2d::CreateMat: edge without SecondBisector");
                self.bisector_map.insert(noofbisectors, second_bisector);
                noofbisectors += 1;
                self.edge_list.next();
            }
        }

        // `cxx:836-857` / `:1612-1626`.
        self.roots = MatTList::new();
        let last = self.map(noofbisectors - 1);
        if self.bis(last).bisector_number() == -1 {
            // The root is the fused bisector's own list, and the distance of its
            // first edge is the distance of the issue point of the fusion.
            self.roots = self.bis(last).list().clone();
            self.roots.first();
            let current = *self.roots.current();
            let first_edge = self
                .bis(current)
                .first_edge()
                .expect("MAT2d_Mat2d::CreateMat: root without FirstEdge");
            let distance = self.bis(last).dist_issue_point();
            self.edg_mut(first_edge).set_distance(distance);
        } else {
            for i in beginbisector..noofbisectors {
                let id = self.map(i);
                self.roots.back_add(id);
            }
        }
    }
}
