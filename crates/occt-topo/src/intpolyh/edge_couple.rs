//! `IntPolyh_Edge` (`IntPolyh_Edge.hxx:26-83`) and `IntPolyh_Couple`
//! (`IntPolyh_Couple.hxx:31-94`) - the two small index records of the
//! triangulation intersection.

/// `IntPolyh_Edge` (`hxx:26`): an edge between two nodes, shared by two
/// triangles. Unset indices are `-1`, as in the `.hxx:34-37`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IntPolyhEdge {
    pub point1: i32,
    pub point2: i32,
    pub triangle1: i32,
    pub triangle2: i32,
}

impl Default for IntPolyhEdge {
    fn default() -> Self {
        Self {
            point1: -1,
            point2: -1,
            triangle1: -1,
            triangle2: -1,
        }
    }
}

impl IntPolyhEdge {
    /// `IntPolyh_Edge(thePoint1, thePoint2, theTriangle1, theTriangle2)`
    /// (`hxx:41-50`).
    pub fn new(point1: i32, point2: i32, triangle1: i32, triangle2: i32) -> Self {
        Self {
            point1,
            point2,
            triangle1,
            triangle2,
        }
    }
}

/// `IntPolyh_Couple` (`hxx:31`): an unordered pair of interfering triangles
/// plus the angle between their normals.
#[derive(Debug, Clone, Copy)]
pub struct IntPolyhCouple {
    pub index1: i32,
    pub index2: i32,
    pub analyzed: bool,
    pub angle: f64,
}

impl Default for IntPolyhCouple {
    /// `IntPolyh_Couple()` (`hxx:38-45`): `myAngle = -2.0`.
    fn default() -> Self {
        Self {
            index1: -1,
            index2: -1,
            analyzed: false,
            angle: -2.0,
        }
    }
}

impl IntPolyhCouple {
    /// `IntPolyh_Couple(theTriangle1, theTriangle2, theAngle = -2.0)`
    /// (`hxx:48-54`).
    pub fn new(index1: i32, index2: i32, angle: f64) -> Self {
        Self {
            index1,
            index2,
            analyzed: false,
            angle,
        }
    }

    /// `IsEqual(theOther)` (`hxx:81-85`): the pair is unordered.
    pub fn is_equal(&self, other: &Self) -> bool {
        self.index1 == other.index1 && self.index2 == other.index2
            || self.index1 == other.index2 && self.index2 == other.index1
    }

    /// The `std::hash<IntPolyh_Couple>` specialisation (`hxx:113-116`): both
    /// indices are sorted before hashing, matching `IsEqual`.
    pub fn hash_key(&self) -> (i32, i32) {
        if self.index1 > self.index2 {
            (self.index2, self.index1)
        } else {
            (self.index1, self.index2)
        }
    }
}
