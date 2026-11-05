//! `MAT2d_Mat2d::LoadBisectorsToRemove` (`MAT2d_Mat2d.cxx:1632-1699`) and
//! `MAT2d_Mat2d::Intersect` (`cxx:1701-1826`).
//!
//! `Intersect` walks the descendants of the two bisectors that surround an edge
//! and records which branch of the tree has to be cut. `LoadBisectorsToRemove`
//! turns the two candidate branches into one entry of the three
//! `...toremove` maps, so that the next pass of the main loop can nullify the
//! subtree between them.

use occt_core::precision::INFINITE;

use super::Mat2dMat2d;
use crate::mat::BisectorId;
use crate::mat2d::Mat2dTool2d;

impl Mat2dMat2d {
    /// `Intersect(aTool, aSide, NoOfBisectorsToRemove, Bisector1, Bisector2)`
    /// (`cxx:1719-1826`).
    ///
    /// * `aside == 0`: `Bisector1` against the leftmost descendants of
    ///   `Bisector2`, and `Bisector2` against the rightmost descendants of
    ///   `Bisector1`.
    /// * `aside == 1`: `Bisector1` against its own leftmost and rightmost
    ///   descendants.
    /// * `aside == 2`: `Bisector2` against its own leftmost and rightmost
    ///   descendants.
    ///
    /// `distant` keeps its OCCT initialisation (`cxx:1776-1778`) so that a
    /// bisector whose `List()` is empty still contributes the infinite distance
    /// the branch was entered with.
    #[allow(unused_assignments)]
    pub(crate) fn intersect(
        &mut self,
        tool: &mut Mat2dTool2d,
        aside: i32,
        noofbisectortoremove: &mut usize,
        firstbisector: BisectorId,
        secondbisector: BisectorId,
    ) {
        let mut distance = [INFINITE, INFINITE];
        let mut firstbisectortoremove: [Option<BisectorId>; 2] = [None, None];
        let mut lastbisectortoremove: [Option<BisectorId>; 2] = [None, None];

        for bisectornumber in 0..2 {
            // `cxx:1734-1753`.
            firstbisectortoremove[bisectornumber] = if aside == 0 {
                if bisectornumber == 0 {
                    Some(secondbisector)
                } else {
                    Some(firstbisector)
                }
            } else if aside == 1 {
                Some(firstbisector)
            } else {
                Some(secondbisector)
            };
            let Some(mut lastbisector) = firstbisectortoremove[bisectornumber] else {
                continue;
            };

            let mut previousbisector;
            if aside == 0 {
                previousbisector = lastbisector;
            } else {
                if self.bis(lastbisector).list().is_empty() {
                    continue;
                }
                previousbisector = if bisectornumber == 0 {
                    self.bis(lastbisector)
                        .first_bisector()
                        .expect("MAT2d_Mat2d::Intersect: empty List")
                } else {
                    self.bis(lastbisector)
                        .last_bisector()
                        .expect("MAT2d_Mat2d::Intersect: empty List")
                };
            }

            let mut distant = distance[bisectornumber];
            // `cxx:1779-1811`.
            while !self.bis(previousbisector).list().is_empty() {
                previousbisector = if bisectornumber == 0 {
                    self.bis(previousbisector)
                        .first_bisector()
                        .expect("MAT2d_Mat2d::Intersect: empty List")
                } else {
                    self.bis(previousbisector)
                        .last_bisector()
                        .expect("MAT2d_Mat2d::Intersect: empty List")
                };

                if aside == 1 || (aside == 0 && bisectornumber == 0) {
                    // `cxx:1794-1798`: the intersection is written back into
                    // `firstbisector`; the previous bisector's `FirstParameter`
                    // is restored.
                    distant = self
                        .intersect_pair(tool, firstbisector, previousbisector, false, true)
                        .0;
                } else {
                    // `cxx:1799-1807`: the intersection is written back into
                    // `secondbisector`; the previous bisector's
                    // `SecondParameter` is restored.
                    distant = self
                        .intersect_pair(tool, previousbisector, secondbisector, true, false)
                        .0;
                }

                if distant < INFINITE {
                    distance[bisectornumber] = distant;
                    lastbisectortoremove[bisectornumber] = Some(lastbisector);
                }

                lastbisector = previousbisector;
            }
        }

        // `cxx:1813-1825`.
        self.load_bisectors_to_remove(
            noofbisectortoremove,
            distance[0],
            distance[1],
            firstbisectortoremove,
            lastbisectortoremove,
        );
    }

    /// `LoadBisectorsToRemove` (`cxx:1632-1699`).
    fn load_bisectors_to_remove(
        &mut self,
        noofbisectorstoremove: &mut usize,
        distance1: f64,
        distance2: f64,
        firstbisectortoremove: [Option<BisectorId>; 2],
        lastbisectortoremove: [Option<BisectorId>; 2],
    ) {
        // `cxx:1651-1661`.
        let index: i32 = if distance1 < INFINITE && distance2 == INFINITE {
            0
        } else if distance2 < INFINITE && distance1 == INFINITE {
            1
        } else {
            -1
        };
        if index == -1 {
            return;
        }
        let index = index as usize;
        let first = firstbisectortoremove[index]
            .expect("MAT2d_Mat2d::LoadBisectorsToRemove: null first bisector");
        let last = lastbisectortoremove[index]
            .expect("MAT2d_Mat2d::LoadBisectorsToRemove: null last bisector");

        // `cxx:1666-1680`.
        let mut found = *noofbisectorstoremove as i32;
        for j in 0..*noofbisectorstoremove {
            let one = *self
                .bisector_one_to_remove
                .get(&j)
                .expect("MAT2d_Mat2d::LoadBisectorsToRemove: bisectoronetoremove");
            if self.bis(one).bisector_number() == self.bis(first).bisector_number() {
                found = j as i32;
                let two = *self
                    .bisector_two_to_remove
                    .get(&j)
                    .expect("MAT2d_Mat2d::LoadBisectorsToRemove: bisectortwotoremove");
                if self.bis(two).bisector_number() < self.bis(last).bisector_number() {
                    found = -1;
                }
                break;
            }
        }

        // `cxx:1682-1698`.
        if found != -1 {
            let found = found as usize;
            self.bisector_one_to_remove.insert(found, first);
            self.bisector_two_to_remove.insert(found, last);
            self.type_of_bisector_to_remove
                .insert(found, index as i32 + 1);
            if found == *noofbisectorstoremove {
                *noofbisectorstoremove += 1;
            }
        }
    }
}
