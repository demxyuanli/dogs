//! `TopOpeBRepBuild_GIter` — iterator on the true cases of a `GTopo`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_GIter.hxx:26-52`,
//! `.cxx:17-94`.
//!
//! OCCT keeps `void* mypG` (`hxx:50`) and dereferences it through the `MYGTOPO`
//! macro (`cxx:17`). The iterator only reads the `GTopo`, which is a value
//! type, so the port stores the `GTopo` it was initialized with.

use crate::abs::State;

use super::g_topo::GTopo;

/// `class TopOpeBRepBuild_GIter` (`hxx:29`).
#[derive(Debug, Clone, Default)]
pub struct GIter {
    /// `myII` (`hxx:49`).
    my_ii: i32,
    /// `mypG` (`hxx:50`) — `None` before the first `Init(G)`.
    my_g: Option<GTopo>,
}

impl GIter {
    /// `TopOpeBRepBuild_GIter()` (`cxx:21-25`).
    pub fn new() -> Self {
        Self {
            my_ii: 0,
            my_g: None,
        }
    }

    /// `TopOpeBRepBuild_GIter(const GTopo& G)` (`cxx:27-33`).
    pub fn with_topo(g: &GTopo) -> Self {
        let mut i = Self::new();
        i.init_with(g);
        i
    }

    /// `Find()` (`cxx:35-47`): skip to the next `II` whose case is true.
    fn find(&mut self) {
        while self.my_ii <= 8 {
            let b = self.g().value_index(self.my_ii);
            if b {
                break;
            }
            self.my_ii += 1;
        }
    }

    fn g(&self) -> &GTopo {
        self.my_g.as_ref().expect("GIter::Init(G) not called")
    }

    /// `Init()` (`cxx:49-53`).
    pub fn init(&mut self) {
        self.my_ii = 0;
        self.find();
    }

    /// `Init(const GTopo& G)` (`cxx:55-60`).
    pub fn init_with(&mut self, g: &GTopo) {
        self.my_g = Some(g.clone());
        self.init();
    }

    /// `More()` (`cxx:62-73`).
    pub fn more(&self) -> bool {
        if self.my_ii <= 8 {
            self.g().value_index(self.my_ii)
        } else {
            false
        }
    }

    /// `Next()` (`cxx:75-79`).
    pub fn next(&mut self) {
        self.my_ii += 1;
        self.find();
    }

    /// `Current(s1, s2)` (`cxx:81-94`): `(UNKNOWN, UNKNOWN)` when the iterator
    /// is exhausted.
    pub fn current(&self) -> (State, State) {
        if !self.more() {
            return (State::Unknown, State::Unknown);
        }
        let (i1, i2) = self.g().index(self.my_ii);
        (self.g().g_state(i1), self.g().g_state(i2))
    }

    /// `Dump(OS)` (`cxx:86-94`) — debug output only.
    pub fn dump(&self) {}
}
