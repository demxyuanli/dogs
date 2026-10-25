//! `TopOpeBRepBuild_GTopo` — the 3x3 `IN`/`ON`/`OUT` case table of a boolean
//! operation between two shapes of given types.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_GTopo.hxx:26-122`,
//! `.cxx:23-473`.

use crate::abs::{ShapeType, State};
use crate::topopebrep::ds::kind::Config;

/// `class TopOpeBRepBuild_GTopo` (`hxx:32`).
#[derive(Debug, Clone)]
pub struct GTopo {
    /// `mycases[3][3]` (`hxx:113`) — indexed `[i1][i2]` with
    /// `TopAbs_IN = 0`, `TopAbs_ON = 1`, `TopAbs_OUT = 2` (`cxx:224-243`).
    my_cases: [[bool; 3]; 3],
    /// `myt1` (`hxx:114`).
    my_t1: ShapeType,
    /// `myt2` (`hxx:115`).
    my_t2: ShapeType,
    /// `myConfig1` (`hxx:116`).
    my_config1: Config,
    /// `myConfig2` (`hxx:117`).
    my_config2: Config,
    /// `myReverseForce` (`hxx:118`).
    my_reverse_force: bool,
    /// `myReverseValue` (`hxx:119`).
    my_reverse_value: bool,
}

impl GTopo {
    /// `TopOpeBRepBuild_GTopo()` (`cxx:23-26`).
    pub fn new() -> Self {
        let mut g = Self {
            my_cases: [[false; 3]; 3],
            my_t1: ShapeType::Shape,
            my_t2: ShapeType::Shape,
            my_config1: Config::UnshGeometry,
            my_config2: Config::UnshGeometry,
            my_reverse_force: false,
            my_reverse_value: false,
        };
        g.reset();
        g
    }

    /// `TopOpeBRepBuild_GTopo(ii..oo, t1, t2, C1, C2)` (`cxx:30-55`).
    #[allow(clippy::too_many_arguments)]
    pub fn with_params(
        ii: bool,
        in_: bool,
        io: bool,
        ni: bool,
        nn: bool,
        no: bool,
        oi: bool,
        on: bool,
        oo: bool,
        t1: ShapeType,
        t2: ShapeType,
        c1: Config,
        c2: Config,
    ) -> Self {
        let mut g = Self::new();
        g.set(ii, in_, io, ni, nn, no, oi, on, oo);
        g.my_t1 = t1;
        g.my_t2 = t2;
        g.my_config1 = c1;
        g.my_config2 = c2;
        g
    }

    /// `Reset()` (`cxx:59-67`).
    pub fn reset(&mut self) {
        self.my_t1 = ShapeType::Shape;
        self.my_t2 = ShapeType::Shape;
        self.my_cases = [[false; 3]; 3];
        self.my_config1 = Config::UnshGeometry;
        self.my_config2 = Config::UnshGeometry;
        self.my_reverse_force = false;
        self.my_reverse_value = false;
    }

    /// `Set(ii..oo)` (`cxx:71-90`).
    #[allow(clippy::too_many_arguments)]
    pub fn set(
        &mut self,
        ii: bool,
        in_: bool,
        io: bool,
        ni: bool,
        nn: bool,
        no: bool,
        oi: bool,
        on: bool,
        oo: bool,
    ) {
        self.my_cases[0][0] = ii;
        self.my_cases[0][1] = in_;
        self.my_cases[0][2] = io;
        self.my_cases[1][0] = ni;
        self.my_cases[1][1] = nn;
        self.my_cases[1][2] = no;
        self.my_cases[2][0] = oi;
        self.my_cases[2][1] = on;
        self.my_cases[2][2] = oo;
    }

    /// `Type(t1, t2)` (`cxx:94-98`).
    pub fn type_(&self) -> (ShapeType, ShapeType) {
        (self.my_t1, self.my_t2)
    }

    /// `ChangeType(t1, t2)` (`cxx:102-106`).
    pub fn change_type(&mut self, t1: ShapeType, t2: ShapeType) {
        self.my_t1 = t1;
        self.my_t2 = t2;
    }

    /// `Config1()` (`cxx:110-113`).
    pub fn config1(&self) -> Config {
        self.my_config1
    }

    /// `Config2()` (`cxx:117-120`).
    pub fn config2(&self) -> Config {
        self.my_config2
    }

    /// `ChangeConfig(C1, C2)` (`cxx:124-128`).
    pub fn change_config(&mut self, c1: Config, c2: Config) {
        self.my_config1 = c1;
        self.my_config2 = c2;
    }

    /// `Value(i1, i2)` (`cxx:132-136`).
    pub fn value(&self, i1: i32, i2: i32) -> bool {
        self.my_cases[i1 as usize][i2 as usize]
    }

    /// `Value(s1, s2)` (`cxx:140-146`).
    pub fn value_state(&self, s1: State, s2: State) -> bool {
        let i1 = self.g_index(s1);
        let i2 = self.g_index(s2);
        self.my_cases[i1 as usize][i2 as usize]
    }

    /// `Value(II)` (`cxx:150-156`).
    pub fn value_index(&self, ii: i32) -> bool {
        let (i1, i2) = self.index(ii);
        self.value(i1, i2)
    }

    /// `ChangeValue(i1, i2, b)` (`cxx:160-163`).
    pub fn change_value(&mut self, i1: i32, i2: i32, b: bool) {
        self.my_cases[i1 as usize][i2 as usize] = b;
    }

    /// `ChangeValue(s1, s2, b)` (`cxx:167-172`).
    pub fn change_value_state(&mut self, s1: State, s2: State, b: bool) {
        let i1 = self.g_index(s1);
        let i2 = self.g_index(s2);
        self.my_cases[i1 as usize][i2 as usize] = b;
    }

    /// `GIndex(s)` (`cxx:176-197`) — `Standard_ProgramError("GIndex : bad
    /// input")` for `UNKNOWN`.
    pub fn g_index(&self, s: State) -> i32 {
        match s {
            State::In => 0,
            State::On => 1,
            State::Out => 2,
            State::Unknown => panic!("GIndex : bad input"),
        }
    }

    /// `GState(i)` (`cxx:201-221`) — `Standard_ProgramError("GState : bad
    /// input")` outside `0..2`.
    pub fn g_state(&self, i: i32) -> State {
        match i {
            0 => State::In,
            1 => State::On,
            2 => State::Out,
            _ => panic!("GState : bad input"),
        }
    }

    /// `Index(II, i1, i2)` (`cxx:225-271`). An out-of-range `II` leaves both
    /// answers untouched in OCCT; the port answers `(0, 0)`.
    pub fn index(&self, ii: i32) -> (i32, i32) {
        match ii {
            0 => (0, 0),
            1 => (0, 1),
            2 => (0, 2),
            3 => (1, 0),
            4 => (1, 1),
            5 => (1, 2),
            6 => (2, 0),
            7 => (2, 1),
            8 => (2, 2),
            _ => (0, 0),
        }
    }

    /// `DumpVal(OS, s1, s2)` (`cxx:275-281`).
    pub fn dump_val(&self, s1: State, s2: State) -> String {
        format!("{}", self.value_state(s1, s2))
    }

    /// `DumpType(OS)` (`cxx:285-290`).
    pub fn dump_type(&self) -> String {
        format!("{}/{}", self.my_t1.to_str(), self.my_t2.to_str())
    }

    /// `static DumpSSB(OS, s1, s2, b)` (`cxx:294-302`).
    pub fn dump_ssb(s1: State, s2: State, b: bool) -> String {
        format!("{} {} : {}", s1.to_str(), s2.to_str(), b)
    }

    /// `Dump(OS, s)` (`cxx:306-368`) — debug output only.
    pub fn dump(&self) {}

    /// `StatesON(s1, s2)` (`cxx:372-393`): the states of `s1`/`s2` under the
    /// `ON` cases, `Standard_ProgramError("Gtopo : StatesON incorrect")` when
    /// `(ON, IN)`/`(ON, OUT)` or `(IN, ON)`/`(OUT, ON)` are both false.
    pub fn states_on(&self) -> (State, State) {
        let mut s1 = State::Unknown;
        if self.value_state(State::On, State::In) {
            s1 = State::In;
        } else if self.value_state(State::On, State::Out) {
            s1 = State::Out;
        }

        let mut s2 = State::Unknown;
        if self.value_state(State::In, State::On) {
            s2 = State::In;
        } else if self.value_state(State::Out, State::On) {
            s2 = State::Out;
        }

        if s1 == State::Unknown || s2 == State::Unknown {
            panic!("Gtopo : StatesON incorrect");
        }

        (s1, s2)
    }

    /// `IsToReverse1()` (`cxx:397-421`).
    pub fn is_to_reverse1(&self) -> bool {
        if self.my_reverse_force {
            return self.my_reverse_value;
        }
        let (s1, s2) = self.states_on();
        if s1 == State::In && s2 == State::In {
            false
        } else {
            s1 == State::In
        }
    }

    /// `IsToReverse2()` (`cxx:425-448`).
    pub fn is_to_reverse2(&self) -> bool {
        if self.my_reverse_force {
            return self.my_reverse_value;
        }
        let (s1, s2) = self.states_on();
        if s1 == State::In && s2 == State::In {
            false
        } else {
            s2 == State::In
        }
    }

    /// `SetReverse(rev)` (`cxx:452-456`).
    pub fn set_reverse(&mut self, rev: bool) {
        self.my_reverse_force = true;
        self.my_reverse_value = rev;
    }

    /// `Reverse()` (`cxx:460-467`) — `Standard_ProgramError("GTopo::ReverseValue
    /// undefined")` when the value was never forced.
    pub fn reverse(&self) -> bool {
        if self.my_reverse_force {
            return self.my_reverse_value;
        }
        panic!("GTopo::ReverseValue undefined");
    }

    /// `CopyPermuted()` (`cxx:471-489`): swaps the two shapes and transposes
    /// the case table (`if (myReverseForce) g.SetReverse(myReverseValue)`).
    pub fn copy_permuted(&self) -> Self {
        let mut g = GTopo::new();

        g.change_type(self.my_t2, self.my_t1);
        g.change_config(self.my_config2, self.my_config1);
        for i in 0..3 {
            for j in 0..3 {
                g.change_value(j as i32, i as i32, self.value(i as i32, j as i32));
            }
        }
        if self.my_reverse_force {
            g.set_reverse(self.my_reverse_value);
        }

        g
    }
}

impl Default for GTopo {
    fn default() -> Self {
        Self::new()
    }
}
