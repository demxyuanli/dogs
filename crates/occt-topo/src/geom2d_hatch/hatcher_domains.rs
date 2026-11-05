//! `Geom2dHatch_Hatcher::ComputeDomains` — the domain state machine.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/Geom2dHatch/Geom2dHatch_Hatcher.cxx:1126-1806`
//! (OCCT 8.0.0). Split out of [`super::hatcher`] only to keep each file under
//! the 800-line limit; the `Hatcher` fields stay in that module and the two
//! `ComputeDomains` overloads are the whole content here.
//!
//! The machine walks the restriction points of one hatching in parameter order
//! and, for each point, decides from `StateBefore` / `StateAfter` and the
//! segment flags whether a domain starts, ends or is emitted on the spot. Every
//! bail-out writes `HatchGen_IncoherentParity` when no domain exists yet and
//! `IsDone(true)` otherwise (`cxx:1330-1341` and its twins); that shared tail
//! is [`Hatcher::abort_domain`].

use crate::abs::State;
use crate::hatch_gen::{Domain, ErrorStatus};

use super::fclass2d::Classifier;
use super::hatcher::Hatcher;

impl Hatcher {
    /// `ComputeDomains()` (`cxx:1126-1136`).
    pub fn compute_domains(&mut self) {
        for ind_h in 1..=self.my_nb_hatchings {
            if self.my_hatchings.contains_key(&ind_h) {
                self.compute_domains_ind(ind_h);
            }
        }
    }

    /// `ComputeDomains(IndH)` (`cxx:1144-1806`).
    pub fn compute_domains_ind(&mut self, ind_h: i32) {
        {
            let hatching = self
                .my_hatchings
                .get_mut(&ind_h)
                .expect("Geom2dHatch_Hatcher::ComputeDomains");
            hatching.clr_domains();
            hatching.set_is_done(false);
        }

        if !self.hatching(ind_h).trim_done() {
            self.trim_ind(ind_h);
        }

        if self.hatching(ind_h).status() != ErrorStatus::NoProblem {
            return;
        }

        let points = self.my_keep_points;
        let segments = self.my_keep_segments;
        let mut i_sav: i32 = 0;
        let mut sav_pnt = false;
        let mut nb_opened_segments: i32 = 0;
        let nb_pnt = self.hatching(ind_h).nb_points();

        if nb_pnt == 0 {
            // `cxx:1175-1186`: no restriction point, classify the point that
            // decides whether the whole hatching is inside the region.
            let classification_point = self.hatching(ind_h).classification_point();
            let classifier =
                Classifier::with_args(&mut self.my_elements, &classification_point, 0.0000001);
            if classifier.state() == State::In {
                let domain = Domain::new();
                self.my_hatchings
                    .get_mut(&ind_h)
                    .expect("Geom2dHatch_Hatcher::ComputeDomains")
                    .add_domain(&domain);
            }
            self.my_hatchings
                .get_mut(&ind_h)
                .expect("Geom2dHatch_Hatcher::ComputeDomains")
                .set_is_done(true);
            return;
        }

        for ipnt in 1..=nb_pnt {
            let no_domain = self.hatching(ind_h).nb_domains() == 0;
            let first_point = ipnt == 1;
            let last_point = ipnt == nb_pnt;

            let cur_pnt = self.hatching(ind_h).point(ipnt).clone();

            let mut state_before = cur_pnt.state_before();
            let mut state_after = cur_pnt.state_after();
            let segment_begin = cur_pnt.segment_beginning();
            let segment_end = cur_pnt.segment_end();

            let mut domain = Domain::new();

            // Initialisations dues au premier point.
            if first_point {
                sav_pnt = false;
                i_sav = 0;
                nb_opened_segments = 0;
                if segment_end && segment_begin {
                    if state_after == State::Unknown {
                        state_after = State::In;
                    }
                    if state_before == State::Unknown {
                        state_before = State::In;
                    }
                    if segments {
                        sav_pnt = true;
                        i_sav = 0;
                    }
                } else if segment_end {
                    if state_after == State::Unknown {
                        state_after = State::In;
                    }
                    if segments {
                        sav_pnt = true;
                        i_sav = 0;
                    }
                } else if segment_begin {
                    if state_before == State::Unknown {
                        state_before = State::In;
                    }
                    if state_before == State::In {
                        sav_pnt = true;
                        i_sav = 0;
                    }
                } else if state_before == State::In {
                    sav_pnt = true;
                    i_sav = 0;
                }
            }

            // Initialisations dues au dernier point.
            if last_point {
                if segment_end && segment_begin {
                    if state_after == State::Unknown {
                        state_after = State::In;
                    }
                    if state_before == State::Unknown {
                        state_before = State::In;
                    }
                } else if segment_end {
                    if state_after == State::Unknown {
                        state_after = State::In;
                    }
                } else if segment_begin && state_before == State::Unknown {
                    state_before = State::In;
                }
            }

            // Cas general.
            let mut to_append = false;

            if segment_end && segment_begin {
                if state_before != State::In && state_after != State::In {
                    self.my_hatchings
                        .get_mut(&ind_h)
                        .expect("Geom2dHatch_Hatcher::ComputeDomains")
                        .set_status(ErrorStatus::IncompatibleStates);
                    return;
                }

                if points {
                    if segments {
                        if !sav_pnt {
                            self.abort_domain(ind_h, no_domain);
                            return;
                        }

                        if i_sav != 0 {
                            domain.set_first_point(self.hatching(ind_h).point(i_sav as usize));
                        }

                        domain.set_second_point(&cur_pnt);
                        to_append = true;
                        sav_pnt = true;
                        i_sav = ipnt as i32;
                    } else {
                        let is_in_in =
                            state_before == State::In && state_after == State::In;
                        if sav_pnt && !is_in_in {
                            self.abort_domain(ind_h, no_domain);
                            return;
                        }

                        domain.set_points(&cur_pnt, &cur_pnt);
                        to_append = true;
                        sav_pnt = false;
                        i_sav = 0;
                    }
                }
            } else if segment_end {
                if segments {
                    if state_after == State::Out {
                        if !sav_pnt {
                            self.abort_domain(ind_h, no_domain);
                            return;
                        }

                        if i_sav != 0 {
                            domain.set_first_point(self.hatching(ind_h).point(i_sav as usize));
                        }

                        domain.set_second_point(&cur_pnt);
                        to_append = true;
                    } else if points {
                        if i_sav != 0 {
                            domain.set_first_point(self.hatching(ind_h).point(i_sav as usize));
                        }

                        domain.set_second_point(&cur_pnt);
                        to_append = true;
                        sav_pnt = true;
                        i_sav = ipnt as i32;
                    }
                } else if state_after == State::In {
                    sav_pnt = true;
                    i_sav = ipnt as i32;
                }

                nb_opened_segments -= 1;
            } else if segment_begin {
                if segments {
                    if state_before == State::Out {
                        sav_pnt = true;
                        i_sav = ipnt as i32;
                    } else if points {
                        if !sav_pnt {
                            self.abort_domain(ind_h, no_domain);
                            return;
                        }

                        if i_sav != 0 {
                            domain.set_first_point(self.hatching(ind_h).point(i_sav as usize));
                        }

                        domain.set_second_point(&cur_pnt);
                        to_append = true;
                        sav_pnt = true;
                        i_sav = ipnt as i32;
                    }
                } else if state_before == State::In {
                    if !sav_pnt {
                        self.abort_domain(ind_h, no_domain);
                        return;
                    }

                    if i_sav != 0 {
                        domain.set_first_point(self.hatching(ind_h).point(i_sav as usize));
                    }

                    domain.set_second_point(&cur_pnt);
                    to_append = true;

                    // `cxx:1494-1502`: kept as OCCT left it, the saved point
                    // moves on instead of being cleared.
                    sav_pnt = true;
                    i_sav = ipnt as i32;
                }

                nb_opened_segments += 1;
            } else {
                // `cxx:1510-1529`: two points delimiting an iso with OUT/OUT
                // transitions are turned into IN transitions.
                if state_before == State::Out && state_after == State::Out && nb_pnt == 2 {
                    if first_point {
                        state_after = State::In;
                    } else {
                        state_before = State::In;
                    }
                }

                if state_before == State::Out && state_after == State::Out {
                    if sav_pnt {
                        self.abort_domain(ind_h, no_domain);
                        return;
                    }

                    if points {
                        domain.set_points(&cur_pnt, &cur_pnt);
                        to_append = true;
                        sav_pnt = true;
                        i_sav = ipnt as i32;
                    }
                } else if state_before == State::Out && state_after == State::In {
                    sav_pnt = true;
                    i_sav = ipnt as i32;
                } else if state_before == State::In && state_after == State::Out {
                    if !sav_pnt {
                        self.abort_domain(ind_h, no_domain);
                        return;
                    }

                    if i_sav != 0 {
                        domain.set_first_point(self.hatching(ind_h).point(i_sav as usize));
                    }

                    domain.set_second_point(&cur_pnt);
                    to_append = true;
                    sav_pnt = false;
                    i_sav = 0;
                } else if state_before == State::In && state_after == State::In {
                    if points {
                        if nb_opened_segments == 0 {
                            if !sav_pnt {
                                self.abort_domain(ind_h, no_domain);
                                // `cxx:1602`: here OCCT continues the point
                                // loop instead of returning.
                                continue;
                            }

                            if i_sav != 0 {
                                domain.set_first_point(self.hatching(ind_h).point(i_sav as usize));
                            }

                            domain.set_second_point(&cur_pnt);
                            to_append = true;
                            sav_pnt = true;
                            i_sav = ipnt as i32;
                        } else if segments {
                            if !sav_pnt {
                                self.abort_domain(ind_h, no_domain);
                                return;
                            }

                            if i_sav != 0 {
                                domain.set_first_point(self.hatching(ind_h).point(i_sav as usize));
                            }

                            domain.set_second_point(&cur_pnt);
                            to_append = true;
                            sav_pnt = true;
                            i_sav = ipnt as i32;
                        } else {
                            if sav_pnt {
                                self.abort_domain(ind_h, no_domain);
                                return;
                            }

                            domain.set_points(&cur_pnt, &cur_pnt);
                            to_append = true;
                            sav_pnt = false;
                            i_sav = 0;
                        }
                    }
                } else {
                    self.my_hatchings
                        .get_mut(&ind_h)
                        .expect("Geom2dHatch_Hatcher::ComputeDomains")
                        .set_status(ErrorStatus::IncompatibleStates);
                    return;
                }
            }

            // Ajout du domaine.
            if to_append {
                self.my_hatchings
                    .get_mut(&ind_h)
                    .expect("Geom2dHatch_Hatcher::ComputeDomains")
                    .add_domain(&domain);
            }

            // Traitement lie au dernier point.
            if last_point {
                domain.set_points_infinite();
                to_append = false;

                if segment_end && segment_begin {
                    if segments {
                        if !sav_pnt {
                            self.abort_domain(ind_h, no_domain);
                            return;
                        }

                        if i_sav != 0 {
                            domain.set_first_point(self.hatching(ind_h).point(i_sav as usize));
                        }

                        to_append = true;
                    }
                } else if segment_end {
                    if state_after == State::In {
                        if !sav_pnt {
                            self.abort_domain(ind_h, no_domain);
                            return;
                        }

                        if i_sav != 0 {
                            domain.set_first_point(self.hatching(ind_h).point(i_sav as usize));
                        }

                        to_append = true;
                    }
                } else if segment_begin {
                    if segments {
                        if !sav_pnt {
                            self.abort_domain(ind_h, no_domain);
                            return;
                        }

                        if i_sav != 0 {
                            domain.set_first_point(self.hatching(ind_h).point(i_sav as usize));
                        }

                        to_append = true;
                    }
                } else if state_after == State::In {
                    if !sav_pnt {
                        self.abort_domain(ind_h, no_domain);
                        return;
                    }

                    if i_sav != 0 {
                        domain.set_first_point(self.hatching(ind_h).point(i_sav as usize));
                    }

                    to_append = true;
                }

                if to_append {
                    self.my_hatchings
                        .get_mut(&ind_h)
                        .expect("Geom2dHatch_Hatcher::ComputeDomains")
                        .add_domain(&domain);
                }
            }
        }

        self.my_hatchings
            .get_mut(&ind_h)
            .expect("Geom2dHatch_Hatcher::ComputeDomains")
            .set_is_done(true);
    }

    /// The `if (!SavPnt) { if (NoDomain) Status(IncoherentParity); else
    /// IsDone(true); return; }` tail repeated through `cxx:1330-1790`.
    fn abort_domain(&mut self, ind_h: i32, no_domain: bool) {
        let hatching = self
            .my_hatchings
            .get_mut(&ind_h)
            .expect("Geom2dHatch_Hatcher::ComputeDomains");
        if no_domain {
            hatching.set_status(ErrorStatus::IncoherentParity);
        } else {
            hatching.set_is_done(true);
        }
    }
}
