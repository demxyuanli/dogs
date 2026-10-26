//! Chaining of start points into section lines:
//! `CheckCoupleAndGetAngle` / `CheckCoupleAndGetAngle2` / `CheckNextStartPoint`
//! (`IntPolyh_MaillageAffinage.cxx:3188-3328`), `StartPointsChain`
//! (`cxx:3333-3611`) and `GetNextChainStartPoint` (`cxx:3617-3780`).

use super::super::arrays::{ArrayOfSectionLines, ArrayOfTangentZones};
use super::super::edge_couple::IntPolyhCouple;
use super::super::section_line::IntPolyhSectionLine;
use super::super::start_point::IntPolyhStartPoint;
use super::{MaillageAffinage, MY_CONFUSION_PRECISION};

/// `CheckCoupleAndGetAngle(T1, T2, Angle, TTrianglesContacts)`
/// (`cxx:3188-3209`).
fn check_couple_and_get_angle(
    t1: i32,
    t2: i32,
    angle: &mut f64,
    contacts: &mut [IntPolyhCouple],
) -> bool {
    for couple in contacts.iter_mut() {
        if !couple.analyzed && couple.index1 == t1 && couple.index2 == t2 {
            couple.analyzed = true;
            *angle = couple.angle;
            return true;
        }
    }
    false
}

/// `CheckCoupleAndGetAngle2(T1, T2, T11, T22, theItCT11, theItCT22, Angle,
/// TTrianglesContacts)` (`cxx:3212-3268`).
///
/// Returns `(found, slot_ct11, slot_ct22, angle)`, the two slots being the
/// positions the `.cxx` hands back through its list iterators.
fn check_couple_and_get_angle2(
    t1: i32,
    t2: i32,
    t11: i32,
    t22: i32,
    contacts: &mut [IntPolyhCouple],
) -> (bool, Option<usize>, Option<usize>, f64) {
    let mut test1 = false;
    let mut test2 = false;
    let mut test3 = false;
    let mut slot_ct11 = None;
    let mut slot_ct22 = None;
    let mut angle = -2.0;

    for (i, couple) in contacts.iter_mut().enumerate() {
        if couple.analyzed {
            continue;
        }
        if couple.index1 == t1 {
            if couple.index2 == t2 {
                test1 = true;
                couple.analyzed = true;
                angle = couple.angle;
            } else if couple.index2 == t22 {
                test2 = true;
                slot_ct11 = Some(i);
                angle = couple.angle;
            }
        } else if couple.index1 == t11 && couple.index2 == t2 {
            test3 = true;
            slot_ct22 = Some(i);
            angle = couple.angle;
        }

        if test1 && test2 && test3 {
            break;
        }
    }

    (test1, slot_ct11, slot_ct22, angle)
}

/// `CheckNextStartPoint(SectionLine, TTangentZones, SP, Prepend)`
/// (`cxx:3277-3324`): 1 when the chain continues, i.e. the point is not a top.
fn check_next_start_point(
    section_line: &mut IntPolyhSectionLine,
    tangent_zones: &mut ArrayOfTangentZones,
    sp: &mut IntPolyhStartPoint,
    prepend: bool,
) -> i32 {
    let mut test = 1;
    if sp.e1 == -1 || sp.e2 == -1 {
        // The tops of the triangle are analysed: skip the ones already stored.
        let fin_ttz = tangent_zones.nb_items();
        let mut found = false;
        for i in 0..fin_ttz {
            let test_sp = tangent_zones.get(i);
            if (sp.u1 - test_sp.u1).abs() < MY_CONFUSION_PRECISION
                && (sp.v1 - test_sp.v1).abs() < MY_CONFUSION_PRECISION
                && (sp.u2 - test_sp.u2).abs() < MY_CONFUSION_PRECISION
                && (sp.v2 - test_sp.v2).abs() < MY_CONFUSION_PRECISION
            {
                found = true;
                break;
            }
        }
        if !found {
            // The top does not belong to the list of tangent zones.
            sp.chainlist = -1;
            let slot = tangent_zones.get_mut(fin_ttz);
            *slot = *sp;
            tangent_zones.increment_nb_items();
        }
        test = 0;
    } else if test != 0 {
        if prepend {
            section_line.prepend(sp);
        } else {
            let idx = section_line.nb_start_points();
            *section_line.change_value(idx) = *sp;
            section_line.increment_nb_start_points();
        }
    }
    test
}

impl MaillageAffinage<'_> {
    /// `StartPointsChain(TSectionLines, TTangentZones)` (`cxx:3333-3611`).
    pub fn start_points_chain(
        &mut self,
        section_lines: &mut ArrayOfSectionLines,
        tangent_zones: &mut ArrayOfTangentZones,
    ) -> i32 {
        // `cxx:3335-3337`: the couples are walked from the mesh; take the list
        // out so the `&self` queries below do not alias it.
        let mut contacts = std::mem::take(&mut self.contacts);

        for i_couple in 0..contacts.len() {
            if contacts[i_couple].analyzed {
                continue;
            }

            let mut section_line_index = section_lines.nb_items();
            // Fill the last section line if it is still empty (eap).
            if section_line_index > 0
                && section_lines.get(section_line_index - 1).nb_start_points() == 0
            {
                section_line_index -= 1;
            } else {
                section_lines.increment_nb_items();
            }

            let mut my_section_line = section_lines.get(section_line_index);
            if my_section_line.get_n() == 0 {
                my_section_line.init(10000); // initialisation of the array of start points
            }

            let t1i = contacts[i_couple].index1;
            let t2i = contacts[i_couple].index2;

            // Start points for the current couple.
            let mut sp1 = IntPolyhStartPoint::default();
            let mut sp2 = IntPolyhStartPoint::default();
            let nb_points = self.starting_points_research(t1i, t2i, &mut sp1, &mut sp2);
            contacts[i_couple].analyzed = true;

            if nb_points == 1 {
                // Particular case top/triangle or edge/edge.
                sp1.chainlist = section_line_index;
                sp1.angle = contacts[i_couple].angle;
                if check_next_start_point(&mut my_section_line, tangent_zones, &mut sp1, false) != 0
                {
                    // Chain of a side.
                    let mut sp11 = IntPolyhStartPoint::default();
                    if sp1.e1 >= 0 {
                        let next_triangle1 = {
                            let e = self.edges1.get(sp1.e1);
                            if e.triangle1 != t1i {
                                e.triangle1
                            } else {
                                e.triangle2
                            }
                        };

                        let mut angle = -2.0;
                        if check_couple_and_get_angle(
                            next_triangle1,
                            t2i,
                            &mut angle,
                            &mut contacts,
                        ) {
                            let nb_points11 =
                                self.next_starting_points_research(next_triangle1, t2i, &sp1, &mut sp11);
                            if nb_points11 == 1 {
                                sp11.chainlist = section_line_index;
                                sp11.angle = angle;

                                if check_next_start_point(
                                    &mut my_section_line,
                                    tangent_zones,
                                    &mut sp11,
                                    false,
                                ) != 0
                                {
                                    let mut end_chain_list = true;
                                    while end_chain_list {
                                        let mut sp_next1 = IntPolyhStartPoint::default();
                                        let test_sp1 = self.get_next_chain_start_point(
                                            &mut contacts,
                                            &sp11,
                                            &mut sp_next1,
                                            &mut my_section_line,
                                            tangent_zones,
                                            false,
                                        );
                                        if test_sp1 == 1 {
                                            sp_next1.chainlist = section_line_index;
                                            if check_next_start_point(
                                                &mut my_section_line,
                                                tangent_zones,
                                                &mut sp_next1,
                                                false,
                                            ) != 0
                                            {
                                                sp11 = sp_next1;
                                            } else {
                                                end_chain_list = false;
                                            }
                                        } else {
                                            end_chain_list = false; // no next point
                                        }
                                    }
                                }
                            } else if nb_points11 > 1 {
                                // The point goes to the tangent zones.
                                let slot = tangent_zones.get_mut(tangent_zones.nb_items());
                                *slot = sp11; // default list number = -1
                                tangent_zones.increment_nb_items();
                            }
                        }
                    }

                    // Chain of the other side.
                    let mut sp12 = IntPolyhStartPoint::default();
                    if sp1.e2 >= 0 {
                        let next_triangle2 = {
                            let e = self.edges2.get(sp1.e2);
                            if e.triangle1 != t2i {
                                e.triangle1
                            } else {
                                e.triangle2
                            }
                        };

                        let mut angle = -2.0;
                        if check_couple_and_get_angle(t1i, next_triangle2, &mut angle, &mut contacts)
                        {
                            let nb_points12 =
                                self.next_starting_points_research(t1i, next_triangle2, &sp1, &mut sp12);
                            if nb_points12 == 1 {
                                sp12.chainlist = section_line_index;
                                sp12.angle = angle;
                                let prepend = true; // eap

                                if check_next_start_point(
                                    &mut my_section_line,
                                    tangent_zones,
                                    &mut sp12,
                                    prepend,
                                ) != 0
                                {
                                    let mut end_chain_list = true;
                                    while end_chain_list {
                                        let mut sp_next1 = IntPolyhStartPoint::default();
                                        let test_sp1 = self.get_next_chain_start_point(
                                            &mut contacts,
                                            &sp12,
                                            &mut sp_next1,
                                            &mut my_section_line,
                                            tangent_zones,
                                            prepend,
                                        );
                                        if test_sp1 == 1 {
                                            sp_next1.chainlist = section_line_index;
                                            if check_next_start_point(
                                                &mut my_section_line,
                                                tangent_zones,
                                                &mut sp_next1,
                                                prepend,
                                            ) != 0
                                            {
                                                sp12 = sp_next1;
                                            } else {
                                                end_chain_list = false;
                                            }
                                        } else {
                                            end_chain_list = false; // no next point
                                        }
                                    }
                                }
                            } else if nb_points12 > 1 {
                                let slot = tangent_zones.get_mut(tangent_zones.nb_items());
                                *slot = sp12;
                                tangent_zones.increment_nb_items();
                            }
                        }
                    }
                }
            } else if nb_points == 2 {
                // The start points are input in the array.
                sp1.chainlist = section_line_index;
                sp1.angle = contacts[i_couple].angle;
                if check_next_start_point(&mut my_section_line, tangent_zones, &mut sp1, false) != 0
                {
                    // Chain of a side.
                    let mut end_chain_list = true;
                    while end_chain_list {
                        let mut sp_next2 = IntPolyhStartPoint::default();
                        let test_sp2 = self.get_next_chain_start_point(
                            &mut contacts,
                            &sp1,
                            &mut sp_next2,
                            &mut my_section_line,
                            tangent_zones,
                            false,
                        );
                        if test_sp2 == 1 {
                            sp_next2.chainlist = section_line_index;
                            if check_next_start_point(
                                &mut my_section_line,
                                tangent_zones,
                                &mut sp_next2,
                                false,
                            ) != 0
                            {
                                sp1 = sp_next2;
                            } else {
                                end_chain_list = false;
                            }
                        } else {
                            end_chain_list = false; // no next point
                        }
                    }
                }

                sp2.chainlist = section_line_index;
                sp2.angle = contacts[i_couple].angle;
                let prepend = true; // eap

                if check_next_start_point(&mut my_section_line, tangent_zones, &mut sp2, prepend) != 0
                {
                    // Chain of the other side.
                    let mut end_chain_list = true;
                    while end_chain_list {
                        let mut sp_next2 = IntPolyhStartPoint::default();
                        let test_sp2 = self.get_next_chain_start_point(
                            &mut contacts,
                            &sp2,
                            &mut sp_next2,
                            &mut my_section_line,
                            tangent_zones,
                            prepend,
                        );
                        if test_sp2 == 1 {
                            sp_next2.chainlist = section_line_index;
                            if check_next_start_point(
                                &mut my_section_line,
                                tangent_zones,
                                &mut sp_next2,
                                prepend,
                            ) != 0
                            {
                                sp2 = sp_next2;
                            } else {
                                end_chain_list = false;
                            }
                        } else {
                            end_chain_list = false; // no next point
                        }
                    }
                }
            } else if (3..7).contains(&nb_points) {
                // More than two start points: only the first one is stored.
                sp1.chainlist = section_line_index;
                check_next_start_point(&mut my_section_line, tangent_zones, &mut sp1, false);
            }

            section_lines.get_mut(section_line_index).copy_from(&my_section_line);
        }

        self.contacts = contacts;
        1
    }

    /// `GetNextChainStartPoint(SP, SPNext, MySectionLine, TTangentZones,
    /// Prepend)` (`cxx:3617-3780`).
    ///
    /// `contacts` replaces the `.cxx` member list so the caller keeps exclusive
    /// access while the `&self` queries run.
    #[allow(clippy::too_many_arguments)]
    fn get_next_chain_start_point(
        &self,
        contacts: &mut [IntPolyhCouple],
        sp: &IntPolyhStartPoint,
        sp_next: &mut IntPolyhStartPoint,
        section_line: &mut IntPolyhSectionLine,
        tangent_zones: &mut ArrayOfTangentZones,
        prepend: bool,
    ) -> i32 {
        let mut nb_points = 0;

        if sp.e1 >= 0 && sp.e2 == -2 {
            // The point is on an edge of T1.
            let next_triangle1 = {
                let e = self.edges1.get(sp.e1);
                if e.triangle1 != sp.t1 {
                    e.triangle1
                } else {
                    e.triangle2
                }
            };
            let mut angle = -2.0;
            if check_couple_and_get_angle(next_triangle1, sp.t2, &mut angle, contacts) {
                nb_points = self.next_starting_points_research(next_triangle1, sp.t2, sp, sp_next);
                if nb_points != 1 {
                    if nb_points > 1 {
                        // The new point is checked; `cxx:3644` keeps the
                        // returned count, so the chain still stops here.
                        check_next_start_point(section_line, tangent_zones, sp_next, prepend);
                    } else {
                        nb_points = 0;
                    }
                } else {
                    sp_next.angle = angle;
                }
            } else {
                nb_points = 0; // this couple does not intersect
            }
        } else if sp.e1 == -2 && sp.e2 >= 0 {
            // The point is on an edge of T2.
            let next_triangle2 = {
                let e = self.edges2.get(sp.e2);
                if e.triangle1 != sp.t2 {
                    e.triangle1
                } else {
                    e.triangle2
                }
            };
            let mut angle = -2.0;
            if check_couple_and_get_angle(sp.t1, next_triangle2, &mut angle, contacts) {
                nb_points = self.next_starting_points_research(sp.t1, next_triangle2, sp, sp_next);
                if nb_points != 1 {
                    if nb_points > 1 {
                        check_next_start_point(section_line, tangent_zones, sp_next, prepend);
                    } else {
                        nb_points = 0;
                    }
                } else {
                    sp_next.angle = angle;
                }
            } else {
                nb_points = 0;
            }
        } else if sp.e1 == -2 && sp.e2 == -2 {
            // No edge is touched or cut.
            nb_points = 0;
        } else if sp.e1 >= 0 && sp.e2 >= 0 {
            // The point is located on two edges.
            let next_triangle1 = {
                let e = self.edges1.get(sp.e1);
                if e.triangle1 != sp.t1 {
                    e.triangle1
                } else {
                    e.triangle2
                }
            };
            let next_triangle2 = {
                let e = self.edges2.get(sp.e2);
                if e.triangle1 != sp.t2 {
                    e.triangle1
                } else {
                    e.triangle2
                }
            };

            let (found, slot_ct11, slot_ct22, angle) =
                check_couple_and_get_angle2(next_triangle1, next_triangle2, sp.t1, sp.t2, contacts);
            if found {
                nb_points =
                    self.next_starting_points_research(next_triangle1, next_triangle2, sp, sp_next);
                if nb_points != 1 {
                    if nb_points > 1 {
                        // The new point is checked (the `.cxx` keeps the result
                        // of the call without branching on it).
                        check_next_start_point(section_line, tangent_zones, sp_next, prepend);
                    }
                    nb_points = 0;
                } else {
                    // NbPoints == 1.
                    sp_next.angle = angle;
                    if let Some(i) = slot_ct11 {
                        contacts[i].analyzed = true;
                    }
                    if let Some(i) = slot_ct22 {
                        contacts[i].analyzed = true;
                    }
                }
            } else {
                nb_points = 0;
            }
        } else if sp.e1 == -1 || sp.e2 == -1 {
            // The points are tops of a triangle; they go to the intermediary
            // array, which the `.cxx` does not do here either.
        }

        nb_points
    }
}
