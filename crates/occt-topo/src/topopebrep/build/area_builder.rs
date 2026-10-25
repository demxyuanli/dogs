//! `TopOpeBRepBuild_AreaBuilder` — reconstructs the areas of a `LoopSet` out of
//! its loops and a loop classifier.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_AreaBuilder.hxx:34-121`,
//! `.cxx:30-430` (OCCT 8.0.0).
//!
//! An area is a list of loops. `InitAreaBuilder` (`cxx:124-330`) walks the
//! loops of the `LoopSet`: a boundary loop is inserted into the first area it
//! is `IN`, otherwise it goes to `boundaryloops`; a block loop either absorbs
//! the area it is `IN` (ejecting the loops outside it) or opens a new area and
//! pulls in the boundary loops that are `IN` it.
//!
//! `Area1dBuilder`'s `InitAreaBuilder` (`Area1dBuilder.cxx:67-298`) repeats the
//! base body verbatim apart from `OCCT_DEBUG` output, and its
//! `ADD_Loop_TO_LISTOFLoop` / `REM_Loop_FROM_LISTOFLoop` /
//! `ADD_LISTOFLoop_TO_LISTOFLoop` (`Area1dBuilder.cxx:302-413`) only append and
//! remove, exactly like the base ones (`cxx:401-430`); the port therefore keeps
//! a single implementation here and lets `Area1dBuilder` delegate.

use crate::abs::State;

use super::loop_classifier::LoopClassifier;
use super::loop_enum::LoopEnum;
use super::loop_set::LoopSetSource;
use super::loop_shape::LoopHandle;

/// `class TopOpeBRepBuild_AreaBuilder` (`hxx:51`).
pub struct AreaBuilder {
    /// `myArea` (`hxx:116`).
    my_area: Vec<Vec<LoopHandle>>,
    /// `myAreaIterator` (`hxx:117`) — a cursor into `my_area`.
    my_area_iter: usize,
    /// `myLoopIterator` (`hxx:118`) — a cursor into the current area.
    my_loop_iter: usize,
    /// `myUNKNOWNRaise` (`hxx:119`).
    my_unknown_raise: bool,
}

impl AreaBuilder {
    /// `TopOpeBRepBuild_AreaBuilder()` (`cxx:30-33`).
    pub fn new() -> Self {
        Self {
            my_area: Vec::new(),
            my_area_iter: 0,
            my_loop_iter: 0,
            my_unknown_raise: false, // no raise if UNKNOWN state found
        }
    }

    /// `TopOpeBRepBuild_AreaBuilder(LS, LC, ForceClass)` (`cxx:36-43`).
    pub fn with_loop_set(
        ls: &mut dyn LoopSetSource,
        lc: &mut dyn LoopClassifier,
        force_class: bool,
    ) -> Self {
        let mut b = Self::new();
        b.init_area_builder(ls, lc, force_class);
        b
    }

    /// `CompareLoopWithListOfLoop(LC, L, LOL, what)` (`cxx:59-103`).
    ///
    /// `OUT` as soon as one tested loop of `LOL` answers `OUT`; `OUT` for an
    /// empty list.
    pub fn compare_loop_with_list_of_loop(
        lc: &mut dyn LoopClassifier,
        l: &LoopHandle,
        lol: &[LoopHandle],
        what: LoopEnum,
    ) -> State {
        if lol.is_empty() {
            return State::Out;
        }

        let mut state = State::Unknown;
        for cur_l in lol {
            let totest = match what {
                LoopEnum::AnyLoop => true,
                LoopEnum::Boundary => cur_l.read().expect("poisoned Loop lock").is_shape(),
                LoopEnum::Block => !cur_l.read().expect("poisoned Loop lock").is_shape(),
            };
            if totest {
                state = lc.compare(l, cur_l);
                if state == State::Out {
                    // <L> is out of at least one Loop of <LOL>: stop
                    break;
                }
            }
        }
        state
    }

    /// `Atomize(state, newstate)` (`cxx:112-122`).
    ///
    /// `myUNKNOWNRaise` is `false` for every builder of this family (it is set
    /// only by the constructor `cxx:30-33` and never changed), so the
    /// `Standard_DomainError` arm is unreachable and `state` takes `newstate`.
    fn atomize(&self, state: &mut State, newstate: State) {
        if self.my_unknown_raise {
            // `Standard_DomainError_Raise_if(state == UNKNOWN, ...)`
            assert!(*state != State::Unknown, "AreaBuilder : Position Unknown");
        } else {
            *state = newstate;
        }
    }

    /// `InitAreaBuilder(LS, LC, ForceClass)` (`cxx:124-330`).
    pub fn init_area_builder(
        &mut self,
        ls: &mut dyn LoopSetSource,
        lc: &mut dyn LoopClassifier,
        force_class: bool,
    ) {
        // boundaryloops : list of boundary loops out of the areas.
        let mut boundaryloops: Vec<LoopHandle> = Vec::new();

        self.my_area.clear();

        ls.init_loop();
        while ls.more_loop() {
            // process a new loop : L is the new current Loop
            let Some(l) = ls.loop_() else {
                break;
            };
            let boundary_l = l.read().expect("poisoned Loop lock").is_shape();

            // L = Shape and ForceClass  : L is handled as a block
            // L = Shape and !ForceClass : L is handled as a pure Shape
            // L = !Shape                : L is handled as a block
            let traitercommeblock = !boundary_l || force_class;
            if !traitercommeblock {
                // the loop L is a boundary loop:
                // - try to insert it in an existing area, such as L is inside
                //   all the block loops. Only block loops of the area are
                //   compared.
                // - if L could not be inserted, store it in boundaryloops.
                let mut loop_inside = false;
                let mut i = 0;
                while i < self.my_area.len() {
                    if self.my_area[i].is_empty() {
                        i += 1;
                        continue;
                    }
                    let mut state =
                        Self::compare_loop_with_list_of_loop(lc, &l, &self.my_area[i], LoopEnum::Block);
                    if state == State::Unknown {
                        self.atomize(&mut state, State::In);
                    }
                    loop_inside = state == State::In;
                    if loop_inside {
                        break;
                    }
                    i += 1;
                }

                if loop_inside {
                    Self::add_loop_to_list_of_loop(&l, &mut self.my_area[i]);
                } else {
                    Self::add_loop_to_list_of_loop(&l, &mut boundaryloops);
                }
            } else {
                // the loop L is a block loop
                let mut loop_inside = false;
                let mut i = 0;
                while i < self.my_area.len() {
                    if self.my_area[i].is_empty() {
                        i += 1;
                        continue;
                    }
                    let mut state =
                        Self::compare_loop_with_list_of_loop(lc, &l, &self.my_area[i], LoopEnum::AnyLoop);
                    if state == State::Unknown {
                        self.atomize(&mut state, State::In);
                    }
                    loop_inside = state == State::In;
                    if loop_inside {
                        break;
                    }
                    i += 1;
                }

                if loop_inside {
                    // the loop L is IN theArea:
                    // - stop area scan, insert L in theArea.
                    // - remove from the area all the loops outside L
                    // - make a new area with them, unless they are all boundary
                    // - if they are all boundary put them back in boundaryLoops
                    let a_area = i;
                    let mut all_shape = true;
                    let mut removed_loops: Vec<LoopHandle> = Vec::new();
                    let mut j = 0;
                    while j < self.my_area[a_area].len() {
                        let cur = self.my_area[a_area][j].clone();
                        let mut state = lc.compare(&cur, &l);
                        if state == State::Unknown {
                            self.atomize(&mut state, State::In); // not OUT
                        }
                        let loop_outside = state == State::Out;
                        if loop_outside {
                            // remove the loop from the area
                            Self::add_loop_to_list_of_loop(&cur, &mut removed_loops);
                            all_shape =
                                all_shape && cur.read().expect("poisoned Loop lock").is_shape();
                            Self::rem_loop_from_list_of_loop(j, &mut self.my_area[a_area]);
                        } else {
                            j += 1;
                        }
                    }
                    // insert the loop in the area
                    Self::add_loop_to_list_of_loop(&l, &mut self.my_area[a_area]);
                    if !removed_loops.is_empty() {
                        if all_shape {
                            Self::add_list_of_loop_to_list_of_loop(
                                &removed_loops,
                                &mut boundaryloops,
                            );
                        } else {
                            // make a new area with the removed loops
                            self.my_area.push(Vec::new());
                            let last = self.my_area.len() - 1;
                            Self::add_list_of_loop_to_list_of_loop(
                                &removed_loops,
                                &mut self.my_area[last],
                            );
                        }
                    }
                } else {
                    // the loop L is not IN any area:
                    // - create a new area with L.
                    // - insert boundary loops that are IN the new area (and
                    //   remove them from boundaryloops)
                    self.my_area.push(Vec::new());
                    let new_area0 = self.my_area.len() - 1;
                    Self::add_loop_to_list_of_loop(&l, &mut self.my_area[new_area0]);

                    let mut j = 0;
                    while j < boundaryloops.len() {
                        let cur = boundaryloops[j].clone();
                        let mut state = lc.compare(&cur, &l);
                        if state == State::Unknown {
                            self.atomize(&mut state, State::In);
                        }
                        let ashapeinside = state == State::In;
                        let mut ablockinside = false;
                        if ashapeinside {
                            let mut state = lc.compare(&l, &cur);
                            if state == State::Unknown {
                                self.atomize(&mut state, State::In);
                            }
                            ablockinside = state == State::In;
                        }
                        if ashapeinside && ablockinside {
                            Self::add_loop_to_list_of_loop(&cur, &mut self.my_area[new_area0]);
                            Self::rem_loop_from_list_of_loop(j, &mut boundaryloops);
                        } else {
                            j += 1;
                        }
                    }
                }
            }

            ls.next_loop();
        }

        self.init_area();
    }

    /// `InitArea()` (`cxx:334-341`): returns the number of areas.
    pub fn init_area(&mut self) -> i32 {
        self.my_area_iter = 0;
        self.init_loop();
        self.my_area.len() as i32
    }

    /// `MoreArea()` (`cxx:345-349`).
    pub fn more_area(&self) -> bool {
        self.my_area_iter < self.my_area.len()
    }

    /// `NextArea()` (`cxx:353-357`).
    pub fn next_area(&mut self) {
        self.my_area_iter += 1;
        self.init_loop();
    }

    /// `InitLoop()` (`cxx:361-375`): returns the number of loops of the current
    /// area, and leaves a null iterator when there is no current area.
    pub fn init_loop(&mut self) -> i32 {
        if self.more_area() {
            let n = self.my_area[self.my_area_iter].len();
            self.my_loop_iter = 0;
            n as i32
        } else {
            // Create an empty ListIteratorOfListOfLoop
            self.my_loop_iter = 0;
            0
        }
    }

    /// `MoreLoop()` (`cxx:379-383`).
    pub fn more_loop(&self) -> bool {
        self.more_area() && self.my_loop_iter < self.my_area[self.my_area_iter].len()
    }

    /// `NextLoop()` (`cxx:387-390`).
    pub fn next_loop(&mut self) {
        self.my_loop_iter += 1;
    }

    /// `Loop()` (`cxx:394-398`).
    pub fn loop_(&self) -> Option<&LoopHandle> {
        if !self.more_area() {
            return None;
        }
        self.my_area[self.my_area_iter].get(self.my_loop_iter)
    }

    /// `ADD_Loop_TO_LISTOFLoop(L, LOL)` (`cxx:401-408`).
    fn add_loop_to_list_of_loop( l: &LoopHandle, lol: &mut Vec<LoopHandle>) {
        lol.push(l.clone());
    }

    /// `REM_Loop_FROM_LISTOFLoop(ITLOL, LOL)` (`cxx:412-418`).
    ///
    /// `NCollection_List::Remove(iterator)` removes the element the iterator
    /// points at and leaves the iterator on the next one, which the index-based
    /// cursor reproduces by not advancing `it`.
    fn rem_loop_from_list_of_loop( it: usize, lol: &mut Vec<LoopHandle>) {
        if it < lol.len() {
            lol.remove(it);
        }
    }

    /// `ADD_LISTOFLoop_TO_LISTOFLoop(LOL1, LOL2)` (`cxx:422-430`).
    fn add_list_of_loop_to_list_of_loop( lol1: &[LoopHandle], lol2: &mut Vec<LoopHandle>) {
        lol2.extend(lol1.iter().cloned());
    }
}

impl Default for AreaBuilder {
    fn default() -> Self {
        Self::new()
    }
}
