//! Port of `Intf_TangentZone`
//! (`src/ModelingAlgorithms/TKGeomAlgo/Intf/Intf_TangentZone.hxx/.cxx/.lxx`).

use super::section_point::IntfSectionPoint;

/// `Intf_TangentZone` (`Intf_TangentZone.hxx:31-103`).
#[derive(Clone, Debug)]
pub struct IntfTangentZone {
    /// `Result` (`hxx:99`)
    result: Vec<IntfSectionPoint>,
    /// `ParamOnFirstMin` (`hxx:100`)
    param_on_first_min: f64,
    /// `ParamOnFirstMax` (`hxx:101`)
    param_on_first_max: f64,
    /// `ParamOnSecondMin` (`hxx:102`)
    param_on_second_min: f64,
    /// `ParamOnSecondMax` (`hxx:103`)
    param_on_second_max: f64,
}

impl Default for IntfTangentZone {
    /// `Intf_TangentZone()` (`Intf_TangentZone.cxx:25-30`).
    fn default() -> Self {
        Self::new()
    }
}

impl IntfTangentZone {
    /// `Intf_TangentZone()` (`Intf_TangentZone.cxx:25-30`).
    pub fn new() -> Self {
        Self {
            result: Vec::new(),
            param_on_first_min: f64::MAX,
            param_on_first_max: -f64::MAX,
            param_on_second_min: f64::MAX,
            param_on_second_max: -f64::MAX,
        }
    }

    /// `NumberOfPoints()` (`TangentZone.lxx:25-28`).
    pub fn number_of_points(&self) -> usize {
        self.result.len()
    }

    /// `GetPoint(Index)` (`Intf_TangentZone.cxx:226-229`). OCCT uses 1-based
    /// indexing.
    pub fn get_point(&self, index: usize) -> &IntfSectionPoint {
        &self.result[index - 1]
    }

    /// `Append(Pi)` (`Intf_TangentZone.cxx:37-57`).
    pub fn append_point(&mut self, pi: &IntfSectionPoint) {
        self.result.push(*pi);
        if self.param_on_first_min > pi.param_on_first() {
            self.param_on_first_min = pi.param_on_first();
        }
        if self.param_on_second_min > pi.param_on_second() {
            self.param_on_second_min = pi.param_on_second();
        }
        if self.param_on_first_max < pi.param_on_first() {
            self.param_on_first_max = pi.param_on_first();
        }
        if self.param_on_second_max < pi.param_on_second() {
            self.param_on_second_max = pi.param_on_second();
        }
    }

    /// `Append(Tzi)` (`Intf_TangentZone.cxx:64-71`).
    pub fn append_zone(&mut self, tzi: &Self) {
        for ipi in 1..=tzi.number_of_points() {
            let p = *tzi.get_point(ipi);
            self.polygon_insert(&p);
        }
    }

    /// `Insert(Pi)` (`Intf_TangentZone.cxx:79-118`). The OCCT body is commented
    /// out and always returns `false`.
    pub fn insert(&mut self, _pi: &IntfSectionPoint) -> bool {
        false
    }

    /// `PolygonInsert(Pi)` (`Intf_TangentZone.cxx:126-171`).
    pub fn polygon_insert(&mut self, pi: &IntfSectionPoint) {
        let nbp_tz = self.number_of_points();
        if nbp_tz == 0 {
            self.append_point(pi);
            return;
        }
        if pi.param_on_first() >= self.param_on_first_max {
            self.append_point(pi);
        } else if pi.param_on_first() >= self.param_on_first_min {
            self.insert_before(1, pi);
        } else {
            self.append_point(pi);
        }
    }

    /// `InsertAfter(Index, Pi)` (`Intf_TangentZone.cxx:175-195`).
    pub fn insert_after(&mut self, index: usize, pi: &IntfSectionPoint) {
        self.result.insert(index, *pi);
        if self.param_on_first_min > pi.param_on_first() {
            self.param_on_first_min = pi.param_on_first();
        }
        if self.param_on_second_min > pi.param_on_second() {
            self.param_on_second_min = pi.param_on_second();
        }
        if self.param_on_first_max < pi.param_on_first() {
            self.param_on_first_max = pi.param_on_first();
        }
        if self.param_on_second_max < pi.param_on_second() {
            self.param_on_second_max = pi.param_on_second();
        }
    }

    /// `InsertBefore(Index, Pi)` (`Intf_TangentZone.cxx:199-219`).
    pub fn insert_before(&mut self, index: usize, pi: &IntfSectionPoint) {
        self.result.insert(index - 1, *pi);
        if self.param_on_first_min > pi.param_on_first() {
            self.param_on_first_min = pi.param_on_first();
        }
        if self.param_on_second_min > pi.param_on_second() {
            self.param_on_second_min = pi.param_on_second();
        }
        if self.param_on_first_max < pi.param_on_first() {
            self.param_on_first_max = pi.param_on_first();
        }
        if self.param_on_second_max < pi.param_on_second() {
            self.param_on_second_max = pi.param_on_second();
        }
    }

    /// `IsEqual(Other)` (`Intf_TangentZone.cxx:232-247`).
    pub fn is_equal(&self, other: &Self) -> bool {
        if self.result.len() != other.result.len() {
            return false;
        }
        for i in 0..self.result.len() {
            if !self.result[i].is_equal(&other.result[i]) {
                return false;
            }
        }
        true
    }

    /// `Contains(ThePI)` (`Intf_TangentZone.cxx:252-262`).
    pub fn contains(&self, the_pi: &IntfSectionPoint) -> bool {
        self.result.iter().any(|p| the_pi.is_equal(p))
    }

    /// `ParamOnFirst(paraMin, paraMax)` (`TangentZone.lxx:32-36`).
    pub fn param_on_first(&self) -> (f64, f64) {
        (self.param_on_first_min, self.param_on_first_max)
    }

    /// `ParamOnSecond(paraMin, paraMax)` (`TangentZone.lxx:40-44`).
    pub fn param_on_second(&self) -> (f64, f64) {
        (self.param_on_second_min, self.param_on_second_max)
    }

    /// `InfoFirst(segMin, paraMin, segMax, paraMax)`
    /// (`Intf_TangentZone.cxx:266-274`).
    pub fn info_first(&self) -> (i32, f64, i32, f64) {
        let (mut para_min, mut para_max) = self.param_on_first();
        let seg_min = para_min.trunc() as i32;
        para_min -= f64::from(seg_min);
        let seg_max = para_max.trunc() as i32;
        para_max -= f64::from(seg_max);
        (seg_min, para_min, seg_max, para_max)
    }

    /// `InfoSecond(segMin, paraMin, segMax, paraMax)`
    /// (`Intf_TangentZone.cxx:278-286`).
    pub fn info_second(&self) -> (i32, f64, i32, f64) {
        let (mut para_min, mut para_max) = self.param_on_second();
        let seg_min = para_min.trunc() as i32;
        para_min -= f64::from(seg_min);
        let seg_max = para_max.trunc() as i32;
        para_max -= f64::from(seg_max);
        (seg_min, para_min, seg_max, para_max)
    }

    /// `RangeContains(ThePI)` (`Intf_TangentZone.cxx:288-296`).
    pub fn range_contains(&self, the_pi: &IntfSectionPoint) -> bool {
        let (a, b) = self.param_on_first();
        let (c, d) = self.param_on_second();
        a <= the_pi.param_on_first()
            && the_pi.param_on_first() <= b
            && c <= the_pi.param_on_second()
            && the_pi.param_on_second() <= d
    }

    /// `HasCommonRange(Other)` (`Intf_TangentZone.cxx:300-311`).
    pub fn has_common_range(&self, other: &Self) -> bool {
        let (a1, b1) = self.param_on_first();
        let (a2, b2) = self.param_on_second();
        let (c1, d1) = other.param_on_first();
        let (c2, d2) = other.param_on_second();

        ((c1 <= a1 && a1 <= d1) || (c1 <= b1 && b1 <= d1) || (a1 <= c1 && c1 <= b1))
            && ((c2 <= a2 && a2 <= d2) || (c2 <= b2 && b2 <= d2) || (a2 <= c2 && c2 <= b2))
    }
}

impl PartialEq for IntfTangentZone {
    fn eq(&self, other: &Self) -> bool {
        self.is_equal(other)
    }
}
