//! Port of `Intf_SectionPoint`
//! (`src/ModelingAlgorithms/TKGeomAlgo/Intf/Intf_SectionPoint.hxx/.cxx/.lxx`)
//! and the `Intf_PIType` enum (`Intf_PIType.hxx`).

use crate::gp::{GpPnt, GpPnt2d};

/// `Intf_PIType` (`Intf_PIType.hxx:20-28`).
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub enum IntfPIType {
    /// `Intf_EXTERNAL`
    External,
    /// `Intf_FACE`
    Face,
    /// `Intf_EDGE`
    Edge,
    /// `Intf_VERTEX`
    Vertex,
}

/// `Intf_SectionPoint` (`Intf_SectionPoint.hxx:30-122`).
#[derive(Clone, Copy, Debug)]
pub struct IntfSectionPoint {
    /// `myPnt` (`hxx:111`)
    my_pnt: GpPnt,
    /// `DimenObje` (`hxx:112`)
    dimen_obje: IntfPIType,
    /// `IndexO1` (`hxx:113`)
    index_o1: i32,
    /// `IndexO2` (`hxx:114`)
    index_o2: i32,
    /// `ParamObje` (`hxx:115`)
    param_obje: f64,
    /// `DimenTool` (`hxx:116`)
    dimen_tool: IntfPIType,
    /// `IndexT1` (`hxx:117`)
    index_t1: i32,
    /// `IndexT2` (`hxx:118`)
    index_t2: i32,
    /// `ParamTool` (`hxx:119`)
    param_tool: f64,
    /// `Incide` (`hxx:120`)
    incide: f64,
}

impl Default for IntfSectionPoint {
    /// `Intf_SectionPoint()` (`Intf_SectionPoint.cxx:136-148`).
    fn default() -> Self {
        Self::new()
    }
}

impl IntfSectionPoint {
    /// `Intf_SectionPoint()` (`Intf_SectionPoint.cxx:136-148`).
    pub fn new() -> Self {
        Self {
            my_pnt: GpPnt::new(0.0, 0.0, 0.0),
            dimen_obje: IntfPIType::External,
            index_o1: 0,
            index_o2: 0,
            param_obje: 0.0,
            dimen_tool: IntfPIType::External,
            index_t1: 0,
            index_t2: 0,
            param_tool: 0.0,
            incide: 0.0,
        }
    }

    /// `Intf_SectionPoint(Where, DimeO, AddrO1, AddrO2, ParamO, DimeT, AddrT1,
    /// AddrT2, ParamT, Incid)` (`Intf_SectionPoint.cxx:152-174`).
    #[allow(clippy::too_many_arguments)]
    pub fn with_3d(
        where_: &GpPnt,
        dime_o: IntfPIType,
        addr_o1: i32,
        addr_o2: i32,
        param_o: f64,
        dime_t: IntfPIType,
        addr_t1: i32,
        addr_t2: i32,
        param_t: f64,
        incid: f64,
    ) -> Self {
        Self {
            my_pnt: *where_,
            dimen_obje: dime_o,
            index_o1: addr_o1,
            index_o2: addr_o2,
            param_obje: param_o,
            dimen_tool: dime_t,
            index_t1: addr_t1,
            index_t2: addr_t2,
            param_tool: param_t,
            incide: incid,
        }
    }

    /// `Intf_SectionPoint(Where, DimeO, AddrO1, ParamO, DimeT, AddrT1, ParamT,
    /// Incid)` (2D, `Intf_SectionPoint.cxx:176-198`). The 2D point is lifted to
    /// `(x, y, 0)` and the addresses go to `IndexO2` / `IndexT2`.
    pub fn with_2d(
        where_: &GpPnt2d,
        dime_o: IntfPIType,
        addr_o1: i32,
        param_o: f64,
        dime_t: IntfPIType,
        addr_t1: i32,
        param_t: f64,
        incid: f64,
    ) -> Self {
        Self {
            my_pnt: GpPnt::new(where_.x(), where_.y(), 0.0),
            dimen_obje: dime_o,
            index_o1: 0,
            index_o2: addr_o1,
            param_obje: param_o,
            dimen_tool: dime_t,
            index_t1: 0,
            index_t2: addr_t1,
            param_tool: param_t,
            incide: incid,
        }
    }

    /// `Pnt()` (`Intf_SectionPoint.cxx:23-27`).
    pub fn pnt(&self) -> &GpPnt {
        &self.my_pnt
    }

    /// `ParamOnFirst()` (`SectionPoint.lxx:17-20`).
    pub fn param_on_first(&self) -> f64 {
        f64::from(self.index_o2 - 1) + self.param_obje
    }

    /// `ParamOnSecond()` (`SectionPoint.lxx:22-25`).
    pub fn param_on_second(&self) -> f64 {
        f64::from(self.index_t2 - 1) + self.param_tool
    }

    /// `TypeOnFirst()` (`SectionPoint.lxx:27-30`).
    pub fn type_on_first(&self) -> IntfPIType {
        self.dimen_obje
    }

    /// `TypeOnSecond()` (`SectionPoint.lxx:32-35`).
    pub fn type_on_second(&self) -> IntfPIType {
        self.dimen_tool
    }

    /// `InfoFirst(Dim, Add1, Add2, Param)` (`Intf_SectionPoint.cxx:31-37`).
    pub fn info_first_2(&self) -> (IntfPIType, i32, i32, f64) {
        (self.dimen_obje, self.index_o1, self.index_o2, self.param_obje)
    }

    /// `InfoFirst(Dim, Add, Param)` (`Intf_SectionPoint.cxx:41-47`).
    pub fn info_first(&self) -> (IntfPIType, i32, f64) {
        (self.dimen_obje, self.index_o2, self.param_obje)
    }

    /// `InfoSecond(Dim, Add1, Add2, Param)` (`Intf_SectionPoint.cxx:51-57`).
    pub fn info_second_2(&self) -> (IntfPIType, i32, i32, f64) {
        (self.dimen_tool, self.index_t1, self.index_t2, self.param_tool)
    }

    /// `InfoSecond(Dim, Add, Param)` (`Intf_SectionPoint.cxx:61-67`).
    pub fn info_second(&self) -> (IntfPIType, i32, f64) {
        (self.dimen_tool, self.index_t2, self.param_tool)
    }

    /// `Incidence()` (`Intf_SectionPoint.cxx:69-72`).
    pub fn incidence(&self) -> f64 {
        self.incide
    }

    /// `IsEqual(Other)` (`Intf_SectionPoint.lxx:37-40`).
    pub fn is_equal(&self, other: &Self) -> bool {
        self.dimen_obje == other.dimen_obje
            && self.index_o1 == other.index_o1
            && self.index_o2 == other.index_o2
            && self.dimen_tool == other.dimen_tool
            && self.index_t1 == other.index_t1
            && self.index_t2 == other.index_t2
    }

    /// `IsOnSameEdge(Other)` (`Intf_SectionPoint.cxx:75-134`).
    pub fn is_on_same_edge(&self, other: &Self) -> bool {
        let mut is_on = false;
        if self.dimen_obje == IntfPIType::Edge {
            if other.dimen_obje == IntfPIType::Edge {
                is_on = self.index_o1 == other.index_o1 && self.index_o2 == other.index_o2;
            } else if other.dimen_obje == IntfPIType::Vertex {
                is_on = self.index_o1 == other.index_o1 || self.index_o2 == other.index_o1;
            }
        } else if self.dimen_obje == IntfPIType::Vertex {
            if other.dimen_obje == IntfPIType::Edge {
                is_on = self.index_o1 == other.index_o1 || self.index_o1 == other.index_o2;
            } else if other.dimen_obje == IntfPIType::Vertex {
                is_on = self.index_t1 == other.index_t1;
            }
        }
        if !is_on {
            if self.dimen_tool == IntfPIType::Edge {
                if other.dimen_tool == IntfPIType::Edge {
                    is_on = self.index_t1 == other.index_t1 && self.index_t2 == other.index_t2;
                } else if other.dimen_tool == IntfPIType::Vertex {
                    is_on = self.index_t1 == other.index_t1 || self.index_t2 == other.index_t1;
                }
            } else if self.dimen_tool == IntfPIType::Vertex {
                if other.dimen_tool == IntfPIType::Edge {
                    is_on = self.index_t1 == other.index_t1 || self.index_t1 == other.index_t2;
                } else if other.dimen_tool == IntfPIType::Vertex {
                    is_on = self.index_t1 == other.index_t1;
                }
            }
        }
        is_on
    }

    /// `Merge(Other)` (`Intf_SectionPoint.cxx:202-233`).
    pub fn merge(&mut self, other: &mut Self) {
        other.my_pnt = self.my_pnt;
        if self.dimen_obje >= other.dimen_obje {
            other.dimen_obje = self.dimen_obje;
            other.index_o1 = self.index_o1;
            other.index_o2 = self.index_o2;
            other.param_obje = self.param_obje;
        } else {
            self.dimen_obje = other.dimen_obje;
            self.index_o1 = other.index_o1;
            self.index_o2 = other.index_o2;
            self.param_obje = other.param_obje;
        }
        if self.dimen_tool >= other.dimen_tool {
            other.dimen_tool = self.dimen_tool;
            other.index_t1 = self.index_t1;
            other.index_t2 = self.index_t2;
            other.param_tool = self.param_tool;
        } else {
            self.dimen_tool = other.dimen_tool;
            self.index_t1 = other.index_t1;
            self.index_t2 = other.index_t2;
            self.param_tool = other.param_tool;
        }
    }
}

impl PartialEq for IntfSectionPoint {
    fn eq(&self, other: &Self) -> bool {
        self.is_equal(other)
    }
}
