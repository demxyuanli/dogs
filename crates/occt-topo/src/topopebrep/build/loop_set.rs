//! `TopOpeBRepBuild_LoopSet` — a set of loops explored through
//! `InitLoop` / `MoreLoop` / `NextLoop` / `Loop`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_LoopSet.hxx:28-50`,
//! `.cxx:22-61`.

use super::loop_shape::LoopHandle;

/// The `InitLoop` / `MoreLoop` / `NextLoop` / `Loop` sequence of
/// `TopOpeBRepBuild_LoopSet` (`hxx:37-45`).
///
/// `TopOpeBRepBuild_PaveSet` (`PaveSet.hxx:33-52`) derives from `LoopSet` and
/// overrides all four methods, and the area builders take a `LoopSet&`
/// (`AreaBuilder.hxx:44-58`); the port keeps that polymorphism in a trait.
pub trait LoopSetSource {
    /// `ChangeListOfLoop()` (`hxx:35`).
    fn change_list_of_loop(&mut self) -> &mut Vec<LoopHandle>;

    /// `InitLoop()` (`hxx:37`).
    fn init_loop(&mut self);

    /// `MoreLoop()` (`hxx:39`).
    fn more_loop(&self) -> bool;

    /// `NextLoop()` (`hxx:41`).
    fn next_loop(&mut self);

    /// `Loop()` (`hxx:43`).
    fn loop_(&self) -> Option<LoopHandle>;
}

/// `class TopOpeBRepBuild_LoopSet` (`hxx:29`).
#[derive(Debug, Clone, Default)]
pub struct LoopSet {
    /// `myListOfLoop` (`hxx:46`).
    my_list_of_loop: Vec<LoopHandle>,
    /// `myLoopIterator` (`hxx:47`) — a cursor into `myListOfLoop`.
    my_loop_iterator: usize,
    /// `myLoopIndex` (`hxx:48`).
    my_loop_index: i32,
    /// `myNbLoop` (`hxx:49`).
    my_nb_loop: i32,
}

impl LoopSet {
    /// `TopOpeBRepBuild_LoopSet()` (`cxx:22-26`): `myLoopIndex = 1`,
    /// `myNbLoop = 0`.
    pub fn new() -> Self {
        Self {
            my_list_of_loop: Vec::new(),
            my_loop_iterator: 0,
            my_loop_index: 1,
            my_nb_loop: 0,
        }
    }

    /// `myLoopIndex`.
    pub fn loop_index(&self) -> i32 {
        self.my_loop_index
    }

    /// `myNbLoop`, set by the last `InitLoop()`.
    pub fn nb_loop(&self) -> i32 {
        self.my_nb_loop
    }
}

impl LoopSetSource for LoopSet {
    /// `ChangeListOfLoop()` (`cxx:57-60`).
    fn change_list_of_loop(&mut self) -> &mut Vec<LoopHandle> {
        &mut self.my_list_of_loop
    }

    /// `InitLoop()` (`cxx:30-35`).
    fn init_loop(&mut self) {
        self.my_loop_iterator = 0;
        self.my_loop_index = 1;
        self.my_nb_loop = self.my_list_of_loop.len() as i32;
    }

    /// `MoreLoop()` (`cxx:39-42`).
    fn more_loop(&self) -> bool {
        self.my_loop_iterator < self.my_list_of_loop.len()
    }

    /// `NextLoop()` (`cxx:46-50`).
    fn next_loop(&mut self) {
        self.my_loop_index += 1;
        self.my_loop_iterator += 1;
    }

    /// `Loop()` (`cxx:54-58`).
    fn loop_(&self) -> Option<LoopHandle> {
        self.my_list_of_loop.get(self.my_loop_iterator).cloned()
    }
}
